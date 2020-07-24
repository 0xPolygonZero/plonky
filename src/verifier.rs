use anyhow::{bail, ensure, Result};

use crate::partition::get_subgroup_shift;

use crate::gates::evaluate_all_constraints;
use crate::halo::verify_ipa;
use crate::plonk_proof::OldProof;
use crate::plonk_util::{eval_poly, halo_g, halo_n, halo_n_mul, halo_s, pedersen_hash, powers, reduce_with_powers};
use crate::util::{ceil_div_usize, log2_strict};
use crate::{blake_hash_usize_to_curve, fft_precompute, ifft_with_precomputation_power_of_2, msm_execute_parallel, msm_precompute, AffinePoint, Circuit, FftPrecomputation, Field, HaloCurve, MsmPrecomputation, Proof, GRID_WIDTH, NUM_ROUTED_WIRES, NUM_WIRES};

pub const SECURITY_BITS: usize = 128;

#[derive(Debug, Clone)]
pub struct VerificationKey<C: HaloCurve> {
    pub c_constants: Vec<AffinePoint<C>>,
    pub c_s_sigmas: Vec<AffinePoint<C>>,
    pub degree: usize,
    pub num_public_inputs: usize,
    pub security_bits: usize,
    pub pedersen_g_msm_precomputation: Option<MsmPrecomputation<C>>,
    pub fft_precomputation: Option<FftPrecomputation<C::ScalarField>>,
}

impl<C: HaloCurve> From<Circuit<C>> for VerificationKey<C> {
    fn from(circuit: Circuit<C>) -> Self {
        circuit.to_vk()
    }
}

/// Verifies a proof `proof` and some old proofs G points for a given verification key.
/// If `verify_g` is `true`, the function completely verifies the proof, including the
/// linear time check of the G point.
/// If `verify_g` is `false`, returns an `OldProof` to be checked in a later verification.
pub fn verify_proof<C: HaloCurve, InnerC: HaloCurve<BaseField = C::ScalarField>>(
    public_inputs: &[C::ScalarField],
    proof: &Proof<C>,
    old_proofs: &[OldProof<C>],
    vk: &VerificationKey<C>,
    verify_g: bool,
) -> Result<Option<OldProof<C>>> {
    // Verify that the proof parameters are valid.
    check_proof_parameters(proof)?;

    if proof.o_public_inputs.is_some() {
        // Check public inputs.
        verify_public_inputs(public_inputs, &proof, vk.num_public_inputs)?;
    } else {
        ensure!(
            proof.c_pis_quotient.is_some(),
            "Public input information missing from the proof."
        );
    }

    // Observe the transcript and generate the associated challenge points using Fiat-Shamir.
    let challs = proof.get_challenges(public_inputs)?;

    // Check the old proofs' openings.
    verify_old_proof_evaluation(old_proofs, &proof, challs.zeta)?;

    let degree = vk.degree;

    let constraint_terms = evaluate_all_constraints::<C, InnerC>(
        &proof.o_local.o_constants,
        &proof.o_local.o_wires,
        &proof.o_right.o_wires,
        &proof.o_below.o_wires,
    );

    // Evaluate zeta^degree.
    let zeta_power_d = challs.zeta.exp_usize(degree);
    // Evaluate Z_H(zeta).
    let one = <C::ScalarField as Field>::ONE;
    let zero_of_zeta = zeta_power_d - one;

    // Evaluate L_1(zeta) = (zeta^degree - 1) / (degree * (zeta - 1)).
    let lagrange_1_eval =
        zero_of_zeta / (C::ScalarField::from_canonical_usize(degree) * (challs.zeta - one));

    // Get z(zeta), z(g.zeta) from the proof openings.
    let (z_x, z_gx) = (proof.o_local.o_plonk_z, proof.o_right.o_plonk_z);
    // Evaluate the L_1(x) (Z(x) - 1) vanishing term.
    let vanishing_z_1_term = lagrange_1_eval * (z_x - one);

    // Compute Z(zeta) f'(zeta) - Z(g * zeta) g'(zeta), which should vanish on H.
    let mut f_prime = one;
    let mut g_prime = one;
    for i in 0..NUM_ROUTED_WIRES {
        let k_i = get_subgroup_shift::<C::ScalarField>(i);
        let s_id = k_i * challs.zeta;
        let beta_s_id = challs.beta * s_id;
        let beta_s_sigma = challs.beta * proof.o_local.o_plonk_sigmas[i];
        let f_prime_part = proof.o_local.o_wires[i] + beta_s_id + challs.gamma;
        let g_prime_part = proof.o_local.o_wires[i] + beta_s_sigma + challs.gamma;
        f_prime = f_prime * f_prime_part;
        g_prime = g_prime * g_prime_part;
    }
    let vanishing_v_shift_term = f_prime * z_x - g_prime * z_gx;

    let vanishing_terms = [
        vec![vanishing_z_1_term],
        vec![vanishing_v_shift_term],
        constraint_terms,
    ]
    .concat();

    // Compute t(zeta).
    let computed_t_opening = reduce_with_powers(&vanishing_terms, challs.alpha) / zero_of_zeta;

    // Compute the purported opening of t(zeta).
    let purported_t_opening = reduce_with_powers(&proof.o_local.o_plonk_t, zeta_power_d);

    // If the two values differ, the proof is invalid.
    if computed_t_opening != purported_t_opening {
        bail!("Incorrect opening of the t polynomial.");
    }

    let subgroup_generator_n = C::ScalarField::primitive_root_of_unity(log2_strict(vk.degree));

    // Verify that the purported opening of the public input quotient polynomial is valid.
    if let Some(purported_pis_quotient_opening) = proof.o_local.o_pi_quotient {
        let num_public_input_gates = ceil_div_usize(vk.num_public_inputs, NUM_WIRES);
        // Compute the denominator `prod_{pi \in PI} (X - pi)`.
        let pis_quotient_denominator = (0..num_public_input_gates)
            .fold(C::ScalarField::ONE, |acc, i| {
                acc * (challs.zeta - subgroup_generator_n.exp_usize(2 * i))
            });
        // Compute the numerator `WirePoly - PIPoly`.
        let pis_quotient_numerator = (0..NUM_WIRES).fold(C::ScalarField::ZERO, |acc, i| {
            acc + proof.o_local.o_wires[i] * challs.alpha.exp_usize(i)
        }) - eval_poly(
            &public_inputs_to_polynomial(
                public_inputs,
                challs.alpha,
                vk.degree,
                vk.fft_precomputation.as_ref(),
            ),
            challs.zeta,
        );
        let computed_pis_quotient_opening = pis_quotient_numerator / pis_quotient_denominator;

        if computed_pis_quotient_opening != purported_pis_quotient_opening {
            bail!("Incorrect opening of the public inputs quotient polynomial.");
        }
    }

    // Verify polynomial commitment openings.
    let pedersen_h = blake_hash_usize_to_curve(vk.degree);
    let u_curve = blake_hash_usize_to_curve(vk.degree + 1);
    ensure!(
        verify_all_ipas::<C>(
            &vk.c_constants,
            &vk.c_s_sigmas,
            vk.num_public_inputs,
            subgroup_generator_n,
            u_curve,
            pedersen_h,
            &proof,
            old_proofs,
            challs.u,
            challs.v,
            challs.u_scaling,
            challs.zeta,
            &challs.halo_us,
            challs.schnorr_challenge,
            vk.security_bits,
        ),
        "Invalid IPA proof."
    );

    if verify_g {
        let pedersen_g: Vec<_> = (0..vk.degree).map(blake_hash_usize_to_curve::<C>).collect();
        let w = 8; // TODO: Should really be set dynamically based on MSM size.
        let pedersen_g_msm_precomputation =
            msm_precompute(&AffinePoint::batch_to_projective(&pedersen_g), w);

        // Verify that `self.halo_g = <s, G>`.
        if proof.halo_g
            == pedersen_hash(&halo_s(&challs.halo_us), &pedersen_g_msm_precomputation).to_affine()
        {
            Ok(None)
        } else {
            bail!("Invalid G point.");
        }
    } else {
        Ok(Some(OldProof {
            halo_g: proof.halo_g,
            halo_us: challs.halo_us,
        }))
    }
}

/// Verify all IPAs in the given proof using a reduction to a single polynomial.
fn verify_all_ipas<C: HaloCurve>(
    c_constants: &[AffinePoint<C>],
    c_s_sigmas: &[AffinePoint<C>],
    num_public_inputs: usize,
    subgroup_generator_n: C::ScalarField,
    u_curve: AffinePoint<C>,
    pedersen_h: AffinePoint<C>,
    proof: &Proof<C>,
    old_proofs: &[OldProof<C>],
    u: C::ScalarField,
    v: C::ScalarField,
    u_scaling: C::ScalarField,
    zeta: C::ScalarField,
    halo_us: &[C::ScalarField],
    schnorr_challenge: C::ScalarField,
    security_bits: usize,
) -> bool {
    // Reduce all polynomial commitments to a single one, i.e. a random combination of them.
    let c_all: Vec<AffinePoint<C>> = [
        c_constants,
        c_s_sigmas,
        &proof.c_wires,
        &[proof.c_plonk_z],
        &proof.c_plonk_t,
        &old_proofs.iter().map(|p| p.halo_g).collect::<Vec<_>>(),
        &proof.c_pis_quotient.map(|c| vec![c]).unwrap_or_default(),
    ]
    .concat();
    let powers_of_u = powers(u, c_all.len());
    let actual_scalars = powers_of_u
        .iter()
        .map(|u_pow| halo_n::<C>(&u_pow.to_canonical_bool_vec()[..security_bits]))
        .collect::<Vec<_>>();
    let precomputation = msm_precompute(&AffinePoint::batch_to_projective(&c_all), 8);
    let c_reduction = msm_execute_parallel(&precomputation, &actual_scalars);

    // For each opening set, we do a similar reduction, using the actual scalars above.
    let opening_set_reductions: Vec<C::ScalarField> = proof
        .all_opening_sets()
        .iter()
        .map(|opening_set| C::ScalarField::inner_product(&opening_set.to_vec(), &actual_scalars))
        .collect();

    // Then, we reduce the above opening set reductions to a single value.
    let reduced_opening = reduce_with_powers(&opening_set_reductions, v);

    let u_prime =
        halo_n_mul(&u_scaling.to_canonical_bool_vec()[..security_bits], u_curve).to_projective();

    let num_public_input_gates = ceil_div_usize(num_public_inputs, NUM_WIRES);
    let points = [
        if proof.o_public_inputs.is_some() {
            (0..2 * num_public_input_gates)
                .step_by(2)
                .map(|i| subgroup_generator_n.exp_usize(i))
                .collect::<Vec<_>>()
        } else {
            vec![]
        },
        vec![
            zeta,
            zeta * subgroup_generator_n,
            zeta * subgroup_generator_n.exp_usize(GRID_WIDTH),
        ],
    ]
    .concat();
    let halo_bs = points
        .iter()
        .map(|&p| halo_g(p, &halo_us))
        .collect::<Vec<_>>();
    let halo_b = reduce_with_powers(&halo_bs, v);
    verify_ipa::<C>(
        &proof.halo_l,
        &proof.halo_r,
        proof.halo_g,
        c_reduction,
        reduced_opening,
        halo_b,
        halo_us,
        u_prime,
        pedersen_h,
        schnorr_challenge,
        proof.schnorr_proof,
    )
}

/// Verifies that the purported public inputs in a proof match a given set of scalars.
fn verify_public_inputs<C: HaloCurve>(
    public_inputs: &[C::ScalarField],
    proof: &Proof<C>,
    num_public_inputs: usize,
) -> Result<()> {
    if public_inputs.len() != num_public_inputs {
        bail!("Incorrect number of public inputs.")
    }
    if let Some(proof_pis) = &proof.o_public_inputs {
        for i in 0..num_public_inputs {
            // If the value `v` doesn't match the corresponding wire in the `PublicInputGate`, return false.
            if public_inputs[i] != proof_pis[i / NUM_WIRES].o_wires[i % NUM_WIRES] {
                bail!("{}-th public input is incorrect", i);
            }
        }
    }
    Ok(())
}

/// Verifies that the purported openings at `zeta` for the old proofs in `old_proofs` is valid.
fn verify_old_proof_evaluation<C: HaloCurve>(
    old_proofs: &[OldProof<C>],
    proof: &Proof<C>,
    zeta: C::ScalarField,
) -> Result<()> {
    if old_proofs.len() != proof.o_local.o_old_proofs.len() {
        bail!("Incorrect number of old proofs opening.")
    }
    for (i, p) in old_proofs.iter().enumerate() {
        // If the value `v` doesn't match the corresponding wire in the `PublicInputGate`, return false.
        if halo_g(zeta, &p.halo_us) != proof.o_local.o_old_proofs[i] {
            bail!("{}-th old proof opening is incorrect", i);
        }
    }
    Ok(())
}

/// Check that the parameters in a proof are well-formed, i.e,
/// that curve points are on the curve, and field elements are in range.
/// Panics otherwise.
fn check_proof_parameters<C: HaloCurve>(proof: &Proof<C>) -> Result<()> {
    let Proof {
        c_wires,
        c_plonk_z,
        c_plonk_t,
        halo_l,
        halo_r,
        halo_g,
        schnorr_proof,
        ..
    } = proof;
    // Verify that the curve points are valid.
    ensure!(
        c_wires.iter().all(|p| p.is_valid()),
        "A wire polynomial commitment is not on the curve."
    );
    ensure!(
        c_plonk_z.is_valid(),
        "The Z polynomial commitment is not on the curve."
    );
    ensure!(
        c_plonk_t.iter().all(|p| p.is_valid()),
        "A t polynomial commitment is not on the curve."
    );
    ensure!(
        halo_l.iter().all(|p| p.is_valid()),
        "A Halo left point is not on the curve."
    );
    ensure!(
        halo_r.iter().all(|p| p.is_valid()),
        "A Halo right point is not on the curve."
    );
    ensure!(halo_g.is_valid(), "The Halo G point is not on the curve.");
    // Verify that the field elements are valid.
    ensure!(
        proof.all_opening_sets().iter().all(|v| {
            v.to_vec().iter().all(|x| {
                <C::ScalarField as Field>::is_valid_canonical_u64(&x.to_canonical_u64_vec())
            })
        }),
        "An opening element is not in the field."
    );

    // Verify that the Halo vectors have same length.
    ensure!(
        halo_l.len() == halo_r.len(),
        "The Halo L and R vecor don't have the same length."
    );

    // Verify that the Schnorr protocol data are valid.
    ensure!(
        &schnorr_proof.r.is_valid(),
        "The Z polynomial commitment is not on the curve."
    );
    ensure!(
        <C::ScalarField as Field>::is_valid_canonical_u64(&schnorr_proof.z1.to_canonical_u64_vec()),
        "The first element in the Schnorr proof is not in the field."
    );
    ensure!(
        <C::ScalarField as Field>::is_valid_canonical_u64(&schnorr_proof.z2.to_canonical_u64_vec()),
        "The second element in the Schnorr proof is not in the field."
    );

    Ok(())
}

/// Reshapes a vector of field elements to a zero-padded matrix of height `NUM_WIRES`.
fn public_inputs_to_wires_vec<F: Field>(public_inputs: &[F]) -> Vec<Vec<F>> {
    let n = public_inputs.len();
    let mut wires = (0..NUM_WIRES)
        .map(|i| {
            (i..n)
                .step_by(NUM_WIRES)
                .flat_map(|j| vec![public_inputs[j], F::ZERO])
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    if n % NUM_WIRES > 0 {
        for i in n % NUM_WIRES..NUM_WIRES {
            wires[i].extend(&[F::ZERO; 2]);
        }
    }
    wires
}

fn public_inputs_to_polynomial<F: Field>(
    public_inputs: &[F],
    alpha: F,
    degree: usize,
    fft_precomputation: Option<&FftPrecomputation<F>>,
) -> Vec<F> {
    let pis_wires_vec = public_inputs_to_wires_vec(public_inputs);
    let mut scaled_wires_vec = (0..pis_wires_vec[0].len())
        .map(|i| {
            (0..pis_wires_vec.len())
                .map(|j| pis_wires_vec[j][i] * alpha.exp_usize(j))
                .fold(F::ZERO, |acc, x| acc + x)
        })
        .collect::<Vec<_>>();
    for _ in scaled_wires_vec.len()..degree {
        scaled_wires_vec.push(F::ZERO);
    }

    fft_precomputation.map_or_else(
        || ifft_with_precomputation_power_of_2(&scaled_wires_vec, &fft_precompute(degree)),
        |precomp| ifft_with_precomputation_power_of_2(&scaled_wires_vec, precomp),
    )
}
