use anyhow::{anyhow, bail, ensure, Result};

use crate::partition::get_subgroup_shift;
use crate::plonk_challenger::Challenger;
use crate::plonk_gates::evaluate_all_constraints;
use crate::plonk_util::{build_s, halo_g, halo_n, powers, reduce_with_powers};
use crate::util::ceil_div_usize;
use crate::{
    msm_execute_parallel, msm_precompute, AffinePoint, Circuit, Curve, Field, HaloCurve,
    ProjectivePoint, Proof, SchnorrProof, GRID_WIDTH, NUM_ROUTED_WIRES, NUM_WIRES,
};

const SECURITY_BITS: usize = 128;

pub struct VerificationKey<C: Curve> {
    selector_commitments: Vec<AffinePoint<C>>,
    sigma_commitments: Vec<AffinePoint<C>>,
    degree: usize,
    degree_log: usize,
}

pub fn verify_proof_circuit<C: HaloCurve, InnerC: HaloCurve<BaseField = C::ScalarField>>(
    public_inputs: &[C::ScalarField],
    proof: &Proof<C>,
    circuit: &Circuit<C>,
) -> Result<()> {
    // Verify that the proof parameters are valid.
    check_proof_parameters(proof)?;

    // Check public inputs.
    verify_public_inputs(public_inputs, proof)?;

    // Observe the transcript and generate the associated challenge points using Fiat-Shamir.
    let challs = get_challenges(proof, Challenger::new(SECURITY_BITS))?;

    let degree = circuit.degree();

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

    // Verify polynomial commitment openings.
    if !verify_all_ipas::<C>(
        &circuit,
        &proof,
        challs.u,
        challs.v,
        challs.u_scaling,
        challs.zeta,
        challs.ipa_challenges,
        challs.schnorr_challenge,
    ) {
        bail!("Invalid IPA proof.");
    }
    Ok(())
}

// pub fn verify_proof_vk<C: Curve>(
//     public_inputs: &[C::ScalarField],
//     proof: &Proof<C>,
//     vk: &VerificationKey<C>,
// ) -> Result<bool> {
//     let Proof {
//         c_wires,
//         c_plonk_z,
//         c_plonk_t,
//         o_public_inputs,
//         o_local,
//         o_right,
//         o_below,
//         halo_l,
//         halo_r,
//         halo_g,
//     } = proof;
//     // Verify that the proof parameters are valid.
//     check_proof_parameters(proof);

//     // Check public inputs.
//     if !verify_public_inputs(public_inputs, proof) {
//         return Ok(false);
//     }

//     // Observe the transcript and generate the associated challenge points using Fiat-Shamir.
//     let challs = get_challenges(proof, Challenger::new(SECURITY_BITS));

//     // Evaluate zeta^degree.
//     let mut zeta_power_d = challs.zeta.exp_usize(vk.degree_pow);
//     // Evaluate Z_H(zeta).
//     let one = <C::ScalarField as Field>::ONE;
//     let z_of_zeta = zeta_power_d - one;
//     // Evaluate L_1(zeta) = (zeta^degree - 1) / (degree * (zeta - 1)).
//     let lagrange_1_eval =
//         z_of_zeta / (C::ScalarField::from_canonical_usize(vk.degree_pow) * (challs.zeta - one));

//     // Get z(zeta), z(g.zeta) from the proof openings.
//     let (z_x, z_gx) = (proof.o_local.o_plonk_z, proof.o_right.o_plonk_z);
//     // Compute Z(zeta) f'(zeta) - Z(g * zeta) g'(zeta), which should vanish on H.
//     let mut f_prime = one;
//     let mut g_prime = one;
//     for i in 0..NUM_ROUTED_WIRES {
//         let k_i = get_subgroup_shift::<C::ScalarField>(i);
//         let s_id = k_i * challs.zeta;
//         let beta_s_id = challs.beta * s_id;
//         let beta_s_sigma = challs.beta * o_local.o_plonk_sigmas[i];
//         let f_prime_part = o_local.o_wires[i] + beta_s_id + challs.gamma;
//         let g_prime_part = o_local.o_wires[i] + beta_s_sigma + challs.gamma;
//         f_prime = f_prime * f_prime_part;
//         g_prime = g_prime * g_prime_part;
//     }
//     let vanishing_v_shift_term = f_prime * z_x - g_prime * z_gx;

//     // Evaluate the L_1(x) (Z(x) - 1) vanishing term.
//     let vanishing_z_1_term = lagrange_1_eval * (z_x - one);

//     // TODO: Evaluate constraint polynomial
//     let constraint_term = one;

//     // Compute t(zeta).
//     let computed_t_opening = reduce_with_powers(
//         &[vanishing_z_1_term, vanishing_v_shift_term, constraint_term],
//         challs.alpha,
//     );
//     // Compute the purported opening of t(zeta).
//     let purported_t_opening = reduce_with_powers(&proof.o_local.o_plonk_t, zeta_power_d);

//     // If the two values differ, the proof is invalid.
//     if computed_t_opening != purported_t_opening {
//         return Ok(false);
//     }

//     // Verify polynomial commitment openings.
//     // let (u_l, u_r) = verify_all_ipas::<C, InnerC>(&proof, u, v, x, ipa_challenges);
//     todo!()
// }

/// Verify all IPAs in the given proof using a reduction to a single polynomial.
fn verify_all_ipas<C: HaloCurve>(
    circuit: &Circuit<C>,
    proof: &Proof<C>,
    u: C::ScalarField,
    v: C::ScalarField,
    u_scaling: C::ScalarField,
    zeta: C::ScalarField,
    ipa_challenges: Vec<C::ScalarField>,
    schnorr_challenge: C::ScalarField,
) -> bool {
    // Reduce all polynomial commitments to a single one, i.e. a random combination of them.
    let c_all: Vec<AffinePoint<C>> = [
        circuit.c_constants.clone(),
        circuit.c_s_sigmas.clone(),
        proof.c_wires.clone(),
        vec![proof.c_plonk_z],
        proof.c_plonk_t.clone(),
    ]
    .concat();
    let powers_of_u = powers(u, c_all.len());
    let actual_scalars = powers_of_u
        .iter()
        .map(|u_pow| halo_n::<C>(&u_pow.to_canonical_bool_vec()[..circuit.security_bits]))
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

    let u_curve = C::convert(u_scaling) * circuit.u.to_projective();

    let num_public_input_gates = ceil_div_usize(circuit.num_public_inputs, NUM_WIRES);
    let b = circuit.build_halo_b(
        &[
            (0..2 * num_public_input_gates)
                .step_by(2)
                .map(|i| circuit.subgroup_generator_n.exp_usize(i))
                .collect::<Vec<_>>(),
            vec![
                zeta,
                zeta * circuit.subgroup_generator_n,
                zeta * circuit.subgroup_generator_n.exp_usize(GRID_WIDTH),
            ],
        ]
        .concat(),
        v,
    );
    verify_ipa::<C>(
        proof,
        c_reduction,
        reduced_opening,
        b,
        ipa_challenges,
        u_curve,
        circuit.pedersen_h,
        proof.halo_g,
        schnorr_challenge,
        proof.schnorr_proof,
    )
}

/// Verify the final IPA.
fn verify_ipa<C: HaloCurve>(
    proof: &Proof<C>,
    commitment: ProjectivePoint<C>,
    value: C::ScalarField,
    b: Vec<C::ScalarField>,
    ipa_challenges: Vec<C::ScalarField>,
    u_curve: ProjectivePoint<C>,
    pedersen_h: AffinePoint<C>,
    halo_g_curve: AffinePoint<C>,
    schnorr_challenge: C::ScalarField,
    schnorr_proof: SchnorrProof<C>,
) -> bool {
    // Now we begin IPA verification by computing P' and u' as in Protocol 1 of Bulletproofs.
    // In Protocol 1 we compute u' = [x] u, but we leverage to endomorphism, instead computing
    // u' = [n(x)] u.

    // Compute [c] [n(x)] u = [c] u'.
    let u_n_x_c = C::convert(value) * u_curve;
    let p_prime = commitment + u_n_x_c;

    // Compute Q as defined in the Halo paper.
    let mut points = proof.halo_l.clone();
    points.extend(proof.halo_r.iter());
    let mut scalars = ipa_challenges
        .iter()
        .map(|u| u.square())
        .collect::<Vec<_>>();
    scalars.extend(
        ipa_challenges
            .iter()
            .map(|chal| chal.multiplicative_inverse_assuming_nonzero().square()),
    );
    let precomputation = msm_precompute(&AffinePoint::batch_to_projective(&points), 8);
    let q = msm_execute_parallel(&precomputation, &scalars) + p_prime;
    dbg!(q.to_affine());

    // Performing ZK opening protocol.
    // let b = halo_g(point, &ipa_challenges);
    dbg!(build_s(&ipa_challenges));
    let b = C::ScalarField::inner_product(&build_s(&ipa_challenges), &b);
    dbg!(b);
    dbg!(
        (C::convert(schnorr_challenge) * q + schnorr_proof.r).to_affine(),
        (C::convert(schnorr_proof.z1) * (halo_g_curve.to_projective() + C::convert(b) * u_curve)
            + C::convert(schnorr_proof.z2) * pedersen_h.to_projective())
        .to_affine()
    );
    C::convert(schnorr_challenge) * q + schnorr_proof.r
        == C::convert(schnorr_proof.z1) * (halo_g_curve.to_projective() + C::convert(b) * u_curve)
            + C::convert(schnorr_proof.z2) * pedersen_h.to_projective()
}

/// Verifies that the purported public inputs in a proof match a given set of scalars.
fn verify_public_inputs<C: Curve>(
    public_inputs: &[C::ScalarField],
    proof: &Proof<C>,
) -> Result<()> {
    for (i, &v) in public_inputs.iter().enumerate() {
        // If the value `v` doesn't match the corresponding wire in the `PublicInputGate`, return false.
        if v != proof.o_public_inputs[i / NUM_WIRES].o_wires[i % NUM_WIRES] {
            bail!("{}-th public input is incorrect", i);
        }
    }
    Ok(())
}

/// Check that the parameters in a proof are well-formed, i.e,
/// that curve points are on the curve, and field elements are in range.
/// Panics otherwise.
fn check_proof_parameters<C: Curve>(proof: &Proof<C>) -> Result<()> {
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

#[derive(Debug, Clone)]
struct ProofChallenge<C: Curve> {
    beta: C::ScalarField,
    gamma: C::ScalarField,
    alpha: C::ScalarField,
    zeta: C::ScalarField,
    v: C::ScalarField,
    u: C::ScalarField,
    u_scaling: C::ScalarField,
    ipa_challenges: Vec<C::ScalarField>,
    schnorr_challenge: C::ScalarField,
}

// Computes all challenges used in the proof verification.
fn get_challenges<C: Curve>(
    proof: &Proof<C>,
    mut challenger: Challenger<C::BaseField>,
) -> Result<ProofChallenge<C>> {
    let error_msg = "Conversion from base to scalar field failed.";
    challenger.observe_affine_points(&proof.c_wires);
    let (beta_bf, gamma_bf) = challenger.get_2_challenges();
    let beta = C::try_convert_b2s(beta_bf).map_err(|_| anyhow!(error_msg))?;
    let gamma = C::try_convert_b2s(gamma_bf).map_err(|_| anyhow!(error_msg))?;
    challenger.observe_affine_point(proof.c_plonk_z);
    let alpha_bf = challenger.get_challenge();
    let alpha = C::try_convert_b2s(alpha_bf).map_err(|_| anyhow!(error_msg))?;
    challenger.observe_affine_points(&proof.c_plonk_t);
    let zeta_bf = challenger.get_challenge();
    let zeta = C::try_convert_b2s(zeta_bf).map_err(|_| anyhow!(error_msg))?;
    for os in proof.all_opening_sets().iter() {
        for &f in os.to_vec().iter() {
            challenger.observe_element(C::try_convert_s2b(f).map_err(|_| anyhow!(error_msg))?);
        }
    }
    let (v_bf, u_bf) = challenger.get_2_challenges();
    let v = C::try_convert_b2s(v_bf).map_err(|_| anyhow!(error_msg))?;
    let u = C::try_convert_b2s(u_bf).map_err(|_| anyhow!(error_msg))?;

    let u_scaling_bf = challenger.get_challenge();
    let u_scaling = C::try_convert_b2s(u_scaling_bf).map_err(|_| anyhow!(error_msg))?;

    // Compute IPA challenges.
    let mut ipa_challenges = Vec::new();
    for i in 0..proof.halo_l.len() {
        challenger.observe_affine_points(&[proof.halo_l[i], proof.halo_r[i]]);
        let l_challenge = challenger.get_challenge();
        ipa_challenges.push(C::try_convert_b2s(l_challenge).map_err(|_| anyhow!(error_msg))?);
    }

    // Compute challenge for Schnorr protocol.
    challenger.observe_affine_point(proof.schnorr_proof.r);
    let schnorr_challenge_bf = challenger.get_challenge();
    let schnorr_challenge =
        C::try_convert_b2s(schnorr_challenge_bf).map_err(|_| anyhow!(error_msg))?;

    Ok(ProofChallenge {
        beta,
        gamma,
        alpha,
        zeta,
        v,
        u,
        u_scaling,
        ipa_challenges,
        schnorr_challenge,
    })
}
