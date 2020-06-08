use crate::fft::{fft_precompute, ifft_with_precomputation_power_of_2};
use crate::partition::get_subgroup_shift;
use crate::plonk_challenger::Challenger;
use crate::plonk_util::{eval_poly, reduce_with_powers};
use crate::{AffinePoint, Curve, Field, HaloCurve, Proof, NUM_ROUTED_WIRES};
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
    let (v_bf, u_bf, x_bf) = challenger.get_3_challenges();
    let v = C::try_convert_b2s(v_bf).expect("Improbable");
    let u = C::try_convert_b2s(u_bf).expect("Improbable");
    let x = C::try_convert_b2s(x_bf).expect("Improbable");

    // Evaluate zeta^degree.
    let mut zeta_power_d = zeta.exp_usize(vk.degree_pow);
    // Evaluate Z_H(zeta).
    let one = <C::ScalarField as Field>::ONE;
    let z_of_zeta = zeta_power_d - one;
    // Evaluate L_1(zeta) = (zeta^degree - 1) / (degree * (zeta - 1)).
    let lagrange_1_eval =
        z_of_zeta / (C::ScalarField::from_canonical_usize(vk.degree_pow) * (zeta - one));

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
    let vanishing_z_1_term = lagrange_1_eval * (z_x - one);

    // TODO: Evaluate constraint polynomial
    let constraint_term = one;

    // Compute t(zeta).
    let computed_t_opening = reduce_with_powers(
        &[vanishing_z_1_term, vanishing_v_shift_term, constraint_term],
        alpha,
    );
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
        ipa_challenges.push(C::try_convert_b2s(l_challenge).expect("Improbable"));
    }

    // Verify polynomial commitment openings.
    // let (u_l, u_r) = verify_all_ipas::<C, InnerC>(&proof, u, v, x, ipa_challenges);
    todo!()
}

fn public_input_polynomial<F: Field>(public_input: &[F], degree: usize) -> Vec<F> {
    let mut values = vec![F::ZERO; degree];
    (0..public_input.len()).for_each(|i| values[i] = public_input[i]);
    let fft_precomputation = fft_precompute(degree);
    ifft_with_precomputation_power_of_2(&values, &fft_precomputation)
}

// /// Verify all IPAs in the given proof. Return `(u_l, u_r)`, which roughly correspond to `u` and
// /// `u^{-1}` in the Halo paper, respectively.
// fn verify_all_ipas<C: HaloCurve>(
//     proof: &Proof<C>,
//     u: C::ScalarField,
//     v: C::ScalarField,
//     x: C::ScalarField,
//     ipa_challenges: Vec<C::ScalarField>,
// ) -> (Vec<C::ScalarField>, Vec<C::ScalarField>) {
//     // Reduce all polynomial commitments to a single one, i.e. a random combination of them.
//     // TODO: Configure the actual constants and permutations of whatever circuit we wish to verify.
//     // For now, we use a dummy point for each of those polynomial commitments.
//     let dummy_point = builder.constant_affine_point(InnerC::GENERATOR_AFFINE);
//     let c_constants = vec![dummy_point; NUM_CONSTANTS];
//     let c_s_sigmas = vec![dummy_point; NUM_ROUTED_WIRES];
//     let c_all: Vec<AffinePointTarget> = [
//         c_constants,
//         c_s_sigmas,
//         proof.c_wires.clone(),
//         vec![proof.c_plonk_z],
//         proof.c_plonk_t.clone(),
//     ]
//     .concat();
//     let mut c_reduction_muls = Vec::new();
//     let powers_of_u = powers_recursive(builder, u, c_all.len());
//     for (&c, &power) in c_all.iter().zip(powers_of_u.iter()) {
//         let mul = CurveMulOp {
//             scalar: power,
//             point: c,
//         };
//         c_reduction_muls.push(mul);
//     }
//     let c_reduction_msm_result = builder.curve_msm_endo::<InnerC>(&c_reduction_muls);
//     let actual_scalars = c_reduction_msm_result.actual_scalars;
//     let c_reduction = c_reduction_msm_result.msm_result;

//     // For each opening set, we do a similar reduction, using the actual scalars above.
//     let opening_set_reductions: Vec<Target> = proof
//         .all_opening_sets()
//         .iter()
//         .map(|opening_set| {
//             reduce_with_coefficients(builder, &opening_set.to_vec(), &actual_scalars)
//         })
//         .collect();

//     // Then, we reduce the above opening set reductions to a single value.
//     let reduced_opening = reduce_with_powers_recursive(builder, &opening_set_reductions, v);

//     verify_ipa::<C, InnerC>(
//         builder,
//         proof,
//         c_reduction,
//         reduced_opening,
//         x,
//         ipa_challenges,
//     )
// }

// /// Verify the final IPA. Return `(u_l, u_r)`, which roughly correspond to `u` and `u^{-1}` in the
// /// Halo paper, respectively.
// fn verify_ipa<C: HaloCurve>(
//     proof: &Proof,
//     p: AffinePoint<C>,
//     c: C::ScalarField,
//     x: C::ScalarField,
//     ipa_challenges: Vec<C::ScalarField>,
// ) -> (Vec<C::ScalarField>, Vec<C::ScalarField>) {
//     // Now we begin IPA verification by computing P' and u' as in Protocol 1 of Bulletproofs.
//     // In Protocol 1 we compute u' = [x] u, but we leverage to endomorphism, instead computing
//     // u' = [n(x)] u.
//     let u = C::GENERATOR_AFFINE;
//     let u_prime = builder
//         .curve_mul_endo::<InnerC>(CurveMulOp {
//             scalar: x,
//             point: u,
//         })
//         .mul_result;

//     // Compute [c] [n(x)] u = [c] u'.
//     let u_n_x_c = builder.curve_mul::<InnerC>(CurveMulOp {
//         scalar: c,
//         point: u_prime,
//     });
//     let p_prime = builder.curve_add::<InnerC>(p, u_n_x_c);

//     // Compute Q as defined in the Halo paper.
//     let mut q_muls = Vec::new();
//     for i in 0..proof.halo_l_i.len() {
//         let l_i = proof.halo_l_i[i];
//         let r_i = proof.halo_r_i[i];
//         let l_challenge = ipa_challenges[i];
//         let r_challenge = builder.inv(l_challenge);
//         q_muls.push(CurveMulOp {
//             scalar: l_challenge,
//             point: l_i,
//         });
//         q_muls.push(CurveMulOp {
//             scalar: r_challenge,
//             point: r_i,
//         });
//     }
//     let q_msm_result = builder.curve_msm_endo::<InnerC>(&q_muls);
//     let _q = builder.curve_add::<InnerC>(p_prime, q_msm_result.msm_result);

//     // TODO: Perform ZK opening protocol.

//     // Take the square roots of the actual scalars as u_l and u_r, which will be used elsewhere
//     // in the protocol.
//     let mut u_l = Vec::new();
//     let mut u_r = Vec::new();
//     for (i, &scalar) in q_msm_result.actual_scalars.iter().enumerate() {
//         let sqrt = builder.deterministic_square_root(scalar);
//         if i % 2 == 0 {
//             u_l.push(sqrt);
//         } else {
//             u_r.push(sqrt);
//         }
//     }

//     (u_l, u_r)
// }
