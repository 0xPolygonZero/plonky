use crate::plookup::{eval_l_i, padded, SECURITY_BITS};
use crate::proof::{PlookupProof, PlookupProofChallenge};
use anyhow::{ensure, Result};
use plonky::halo::verify_ipa;
use plonky::plonk_util::{halo_g, halo_n, halo_n_mul, powers, reduce_with_powers};
use plonky::util::log2_strict;
use plonky::{
    blake_hash_usize_to_curve, fft_precompute, ifft_with_precomputation_power_of_2,
    msm_execute_parallel, msm_precompute, AffinePoint, Field, HaloCurve, PolynomialCommitment,
};

/// Verifies that a proof is valid for a set `t`.
/// TODO: The verifier should have some auxiliary knowledge of `c_t`. For now, it is stored in the `proof`.
pub fn verify<C: HaloCurve>(t: &[C::ScalarField], proof: &PlookupProof<C>) -> Result<()> {
    let n = proof.n;
    let t = padded(t, n + 1);
    let fft_precomputation = fft_precompute(n + 1);
    let gs = (0..2 * n + 2)
        .map(|i| blake_hash_usize_to_curve::<C>(i))
        .collect::<Vec<_>>();
    let h = blake_hash_usize_to_curve(2 * n + 2);
    let u_curve = blake_hash_usize_to_curve(2 * n + 3);
    let t_coeffs = ifft_with_precomputation_power_of_2(&t, &fft_precomputation);
    let msm_precomputation = msm_precompute(&AffinePoint::batch_to_projective(&gs[..n + 1]), 8);
    let c_t = PolynomialCommitment::coeffs_to_commitment(&t_coeffs, &msm_precomputation, h, false);
    ensure!(c_t.to_affine() == proof.c_t, "Incorrect table commitment");

    let challs = proof.get_challenges()?;
    let PlookupProofChallenge {
        beta,
        gamma,
        zeta,
        alpha,
        ..
    } = challs;
    let generator = C::ScalarField::primitive_root_of_unity(log2_strict(n + 1));
    let beta1 = challs.beta + C::ScalarField::ONE;
    let gamma_beta1 = challs.gamma * beta1;

    let vanishing_z1_term =
        eval_l_i(n + 1, 0, generator, zeta) * (proof.openings.z.local - C::ScalarField::ONE);
    let vanishing_shift_term = (zeta - generator.exp_usize(n))
        * proof.openings.z.local
        * beta1
        * (gamma + proof.openings.f.local)
        * (gamma_beta1 + proof.openings.t.local + beta * proof.openings.t.right)
        - (zeta - generator.exp_usize(n))
            * proof.openings.z.right
            * (gamma_beta1 + proof.openings.h1.local + beta * proof.openings.h1.right)
            * (gamma_beta1 + proof.openings.h2.local + beta * proof.openings.h2.right);
    let eval_last = eval_l_i(n + 1, n, generator, zeta);
    let vanishing_hs_term = eval_last * (proof.openings.h1.local - proof.openings.h2.right);
    let vanishing_last_term = eval_last * (proof.openings.z.local - C::ScalarField::ONE);

    let numerator = reduce_with_powers(
        &[
            vanishing_z1_term,
            vanishing_shift_term,
            vanishing_hs_term,
            vanishing_last_term,
        ],
        alpha,
    );
    let denominator = zeta.exp_usize(n + 1) - C::ScalarField::ONE;

    let purpoted_quotient_opening = numerator / denominator;

    ensure!(
        purpoted_quotient_opening == proof.openings.quotient.local,
        "Incorrect quotient opening"
    );

    let c_all = vec![
        proof.c_f,
        proof.c_t,
        proof.c_h1,
        proof.c_h2,
        proof.c_z,
        proof.c_quotient,
    ];

    ensure!(
        verify_all_ipas(
            &c_all,
            generator,
            u_curve,
            h,
            proof,
            challs.u,
            challs.v,
            challs.u_scaling,
            challs.zeta,
            &challs.halo_us,
            challs.schnorr_challenge,
            SECURITY_BITS,
        ),
        "Invalid IPA proof."
    );

    // TODO: Add option to verify `halo_g` point.
    Ok(())
}

/// Verify all IPAs in the given proof using a reduction to a single polynomial.
fn verify_all_ipas<C: HaloCurve>(
    c_all: &[AffinePoint<C>],
    subgroup_generator_n: C::ScalarField,
    u_curve: AffinePoint<C>,
    pedersen_h: AffinePoint<C>,
    proof: &PlookupProof<C>,
    u: C::ScalarField,
    v: C::ScalarField,
    u_scaling: C::ScalarField,
    zeta: C::ScalarField,
    halo_us: &[C::ScalarField],
    schnorr_challenge: C::ScalarField,
    security_bits: usize,
) -> bool {
    // Reduce all polynomial commitments to a single one, i.e. a random combination of them.
    let powers_of_u = powers(u, c_all.len());
    let actual_scalars = powers_of_u
        .iter()
        .map(|u_pow| halo_n::<C>(&u_pow.to_canonical_bool_vec()[..security_bits]))
        .collect::<Vec<_>>();
    let precomputation = msm_precompute(&AffinePoint::batch_to_projective(&c_all), 8);
    let c_reduction = msm_execute_parallel(&precomputation, &actual_scalars);

    // For each opening set, we do a similar reduction, using the actual scalars above.
    let opening_set_reductions = vec![
        C::ScalarField::inner_product(&actual_scalars, &proof.openings.local()),
        C::ScalarField::inner_product(&actual_scalars, &proof.openings.right()),
    ];

    // Then, we reduce the above opening set reductions to a single value.
    let reduced_opening = reduce_with_powers(&opening_set_reductions, v);

    let u_prime =
        halo_n_mul(&u_scaling.to_canonical_bool_vec()[..security_bits], u_curve).to_projective();

    let points = [zeta, zeta * subgroup_generator_n];
    let halo_bs = points
        .iter()
        .map(|&p| halo_g(p, &halo_us))
        .collect::<Vec<_>>();
    let halo_b = reduce_with_powers(&halo_bs, v);
    verify_ipa::<C>(
        &proof.halo_proof.halo_l,
        &proof.halo_proof.halo_r,
        proof.halo_proof.halo_g,
        c_reduction,
        reduced_opening,
        halo_b,
        halo_us,
        u_prime,
        pedersen_h,
        schnorr_challenge,
        proof.halo_proof.schnorr_proof,
    )
}
