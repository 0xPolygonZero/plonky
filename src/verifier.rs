use crate::fft::{fft_precompute, ifft_with_precomputation_power_of_2};
use crate::plonk_challenger::Challenger;
use crate::plonk_util::{eval_poly, reduce_with_powers};
use crate::partition::get_subgroup_shift;
use crate::{AffinePoint, Curve, Field, Proof, NUM_ROUTED_WIRES};
use anyhow::Result;

const SECURITY_BITS: usize = 128;

pub struct VerificationKey<C: Curve> {
    selector_commitments: Vec<AffinePoint<C>>,
    sigma_commitments: Vec<AffinePoint<C>>,
    degree_log: usize,
    degree_pow: usize,
}

pub fn verify_proof<C: Curve>(
    public_inputs: &[C::ScalarField],
    proof: &Proof<C>,
    vk: &VerificationKey<C>,
) -> Result<bool> {
    let Proof {
        c_wires,
        c_plonk_z,
        c_plonk_t,
        o_public_inputs,
        o_local,
        o_right,
        o_below,
        halo_l,
        halo_r,
        halo_g,
    } = proof;
    // Verify that the curve points are valid.
    assert!(c_wires.iter().all(|p| p.is_valid()));
    assert!(c_plonk_z.is_valid());
    assert!(c_plonk_t.iter().all(|p| p.is_valid()));
    assert!(halo_l.iter().all(|p| p.is_valid()));
    assert!(halo_r.iter().all(|p| p.is_valid()));
    assert!(halo_g.is_valid());
    // Verify that the field elements are valid.
    assert!(proof.all_opening_sets().iter().all(|v| {
        v.to_vec()
            .iter()
            .all(|x| <C::ScalarField as Field>::is_valid_canonical_u64(&x.to_canonical_u64_vec()))
    }));
    // Verify that the public input elements are valid.
    assert!(public_inputs
        .iter()
        .all(|x| <C::ScalarField as Field>::is_valid_canonical_u64(&x.to_canonical_u64_vec())));

    // Verify that the Halo vectors have same length.
    assert!(halo_l.len() == halo_r.len());

    // Observe the transcript and generate the associated challenge points using Fiat-Shamir.
    let mut challenger = Challenger::new(SECURITY_BITS);
    challenger.observe_affine_points(&proof.c_wires);
    let (beta_bf, gamma_bf) = challenger.get_2_challenges();
    let beta = C::try_convert_b2s(beta_bf).expect("Improbable");
    let gamma = C::try_convert_b2s(gamma_bf).expect("Improbable");
    challenger.observe_affine_point(proof.c_plonk_z);
    let alpha_bf = challenger.get_challenge();
    let alpha = C::try_convert_b2s(alpha_bf).expect("Improbable");
    challenger.observe_affine_points(&proof.c_plonk_t);
    let zeta_bf = challenger.get_challenge();
    let zeta = C::try_convert_b2s(zeta_bf).expect("Improbable");
    proof.all_opening_sets().iter().for_each(|os| {
        os.to_vec().iter().for_each(|&f| {
            challenger.observe_element(C::try_convert_s2b(f).expect("Improbable"));
        })
    });
    let (v, u, x) = challenger.get_3_challenges();

    // Evaluate zeta^degree.
    let mut zeta_power_d = zeta.exp_usize(vk.degree_pow);
    // Evaluate Z_H(zeta).
    let one = <C::ScalarField as Field>::ONE;
    let z_of_zeta = zeta_power_d - one;
    // Evaluate L_1(zeta) = (zeta^degree - 1) / (degree * (zeta - 1)).
    let lagrange_1_eval = z_of_zeta / (C::ScalarField::from_canonical_usize(vk.degree_pow) * (zeta - one));

    let pi_poly = public_input_polynomial(&public_inputs, vk.degree_pow);
    let pi_poly_zeta = eval_poly(&pi_poly, zeta);

    // Get z(zeta), z(g.zeta) from the proof openings.
    let (z_x, z_gx) = (proof.o_local.o_plonk_z, proof.o_right.o_plonk_z);
    // Compute Z(zeta) f'(zeta) - Z(g * zeta) g'(zeta), which should vanish on H.
    let mut f_prime = one;
    let mut g_prime = one;
    for i in 0..NUM_ROUTED_WIRES {
        let k_i = get_subgroup_shift::<C::ScalarField>(i);
        let s_id = k_i * zeta;
        let beta_s_id = beta * s_id;
        let beta_s_sigma = beta * o_local.o_plonk_sigmas[i];
        let f_prime_part = o_local.o_wires[i] + beta_s_id + gamma;
        let g_prime_part = o_local.o_wires[i] + beta_s_sigma + gamma;
        f_prime = f_prime * f_prime_part;
        g_prime = g_prime * g_prime_part;
    }
    let vanishing_v_shift_term = f_prime * z_x - g_prime * z_gx;
    
    // Evaluate the L_1(x) (Z(x) - 1) vanishing term.
    let vanishing_z_1_term =  lagrange_1_eval * (z_x - one);

    // TODO: Evaluate constraint polynomial
    let constraint_term = one;

    // Compute t(zeta).
    let computed_t_opening = reduce_with_powers(&[vanishing_z_1_term, vanishing_v_shift_term, constraint_term], alpha);
    // Compute the purported opening of t(zeta).
    let purported_t_opening = reduce_with_powers(&proof.o_local.o_plonk_t, zeta_power_d);

    // If the two values differ, the proof is invalid.
    if computed_t_opening != purported_t_opening {
        return Ok(false);
    }


    // Compute IPA challenges.
    let mut ipa_challenges = Vec::new();
    for i in 0..proof.halo_l.len() {
        challenger.observe_affine_points(&[proof.halo_l[i], proof.halo_r[i]]);
        let l_challenge = challenger.get_challenge();
        ipa_challenges.push(l_challenge);
    }

    // Verify polynomial commitment openings.
    todo!()
}

fn public_input_polynomial<F: Field>(public_input: &[F], degree: usize) -> Vec<F> {
    let mut values = vec![F::ZERO; degree];
    (0..public_input.len()).for_each(|i| values[i] = public_input[i]);
    let fft_precomputation = fft_precompute(degree);
    ifft_with_precomputation_power_of_2(&values, &fft_precomputation)
}
