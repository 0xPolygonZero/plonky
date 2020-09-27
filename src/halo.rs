use crate::plonk_challenger::Challenger;
use crate::plonk_util::{halo_n, halo_n_mul, powers, reduce_with_powers};
use crate::util::log2_strict;
use crate::{msm_execute_parallel, msm_parallel, msm_precompute, AffinePoint, Curve, Field, HaloCurve, PolynomialCommitment, ProjectivePoint, SchnorrProof};
use anyhow::Result;
use rayon::prelude::*;

pub struct OpeningProof<C: HaloCurve> {
    pub halo_l: Vec<AffinePoint<C>>,
    pub halo_r: Vec<AffinePoint<C>>,
    pub halo_g: AffinePoint<C>,
    pub schnorr_proof: SchnorrProof<C>,
}

#[allow(clippy::too_many_arguments)]
pub fn batch_opening_proof<C: HaloCurve>(
    polynomials_coeffs: &[&[C::ScalarField]],
    commitments: &[PolynomialCommitment<C>],
    opening_points: &[C::ScalarField],
    pedersen_g: &[AffinePoint<C>],
    pedersen_h: ProjectivePoint<C>,
    u_curve: AffinePoint<C>,
    u: C::ScalarField,
    v: C::ScalarField,
    u_scaling: C::ScalarField,
    degree: usize,
    security_bits: usize,
    challenger: &mut Challenger<C::BaseField>,
) -> Result<OpeningProof<C>> {
    // Normally we would reduce these lists using powers of u, but for the sake of efficiency
    // (particularly in the recursive verifier) we instead use n(u^i) for each u^i, where n is
    // the injective function related to the Halo endomorphism. Here we compute n(u^i).
    let actual_scalars: Vec<C::ScalarField> = powers(u, polynomials_coeffs.len())
        .iter()
        .map(|u_power| halo_n::<C>(&u_power.to_canonical_bool_vec()[..security_bits]))
        .collect();

    // Reduce the coefficient list to a single set of polynomial coefficients.
    let mut reduced_coeffs = vec![C::ScalarField::ZERO; degree];
    for (i, coeffs) in polynomials_coeffs.iter().enumerate() {
        for (j, &c) in coeffs.iter().enumerate() {
            reduced_coeffs[j] = reduced_coeffs[j] + actual_scalars[i] * c;
        }
    }

    let u_prime =
        halo_n_mul(&u_scaling.to_canonical_bool_vec()[..security_bits], u_curve).to_projective();

    // Make a list of all polynomials' commitment randomness and coefficients, to be reduced later.
    // This must match the order of OpeningSet::to_vec.
    let all_randomness = commitments.iter().map(|c| c.randomness).collect::<Vec<_>>();

    // Final IPA proof.
    let mut halo_a = reduced_coeffs;
    // The Halo b vector is a random combination of the powers of all opening points.
    let mut halo_b = build_halo_b::<C>(opening_points, v, degree);

    let mut halo_g = AffinePoint::batch_to_projective(pedersen_g);
    let mut halo_l = Vec::new();
    let mut halo_r = Vec::new();
    let mut randomness = C::ScalarField::inner_product(&actual_scalars, &all_randomness);
    let degree_pow = log2_strict(degree);
    for j in (1..=degree_pow).rev() {
        let n = 1 << j;
        let middle = n / 2;

        debug_assert_eq!(halo_a.len(), n);
        debug_assert_eq!(halo_b.len(), n);
        debug_assert_eq!(halo_g.len(), n);

        let a_lo = &halo_a[..middle];
        let a_hi = &halo_a[middle..];
        let b_lo = &halo_b[..middle];
        let b_hi = &halo_b[middle..];
        let g_lo = &halo_g[..middle];
        let g_hi = &halo_g[middle..];

        let window_size = 8;

        // We may need to re-generate L_i/R_i a few times with different blinding factors until
        // we get a challenge r such that n(r) is square.
        let u_j = loop {
            let l_j_blinding_factor = C::ScalarField::rand();
            let r_j_blinding_factor = C::ScalarField::rand();

            // L_i = <a_lo, G_hi> + [l_j] H + [<a_lo, b_hi>] U.
            let halo_l_j = msm_parallel(a_lo, g_hi, window_size)
                + C::convert(l_j_blinding_factor) * pedersen_h
                + C::convert(C::ScalarField::inner_product(a_lo, b_hi)) * u_prime;
            // R_i = <a_hi, G_lo> + [r_j] H + [<a_hi, b_lo>] U.
            let halo_r_j = msm_parallel(a_hi, g_lo, window_size)
                + C::convert(r_j_blinding_factor) * pedersen_h
                + C::convert(C::ScalarField::inner_product(a_hi, b_lo)) * u_prime;

            let mut challenger_fork = challenger.clone();
            challenger_fork.observe_proj_points(&[halo_l_j, halo_r_j]);
            let r_bf = challenger_fork.get_challenge();
            let r_sf = r_bf.try_convert::<C::ScalarField>()?;
            let r_bits = &r_sf.to_canonical_bool_vec()[..security_bits];
            let u_j_squared = halo_n::<C>(r_bits);

            if let Some(u_j) = u_j_squared.square_root() {
                let u_squared_inv = u_j_squared.multiplicative_inverse().expect("Improbable");

                halo_l.push(halo_l_j);
                halo_r.push(halo_r_j);
                randomness = randomness
                    + u_j_squared * l_j_blinding_factor
                    + u_squared_inv * r_j_blinding_factor;
                *challenger = challenger_fork;

                break u_j;
            }
        };
        let u_j_inv = u_j.multiplicative_inverse().expect("Improbable");

        halo_a = C::ScalarField::add_slices(&u_j_inv.scale_slice(a_hi), &u_j.scale_slice(a_lo));
        halo_b = C::ScalarField::add_slices(&u_j_inv.scale_slice(b_lo), &u_j.scale_slice(b_hi));
        halo_g = g_lo
            .into_par_iter()
            .zip(g_hi)
            .map(|(&g_lo_i, &g_hi_i)| msm_parallel(&[u_j_inv, u_j], &[g_lo_i, g_hi_i], 4))
            .collect();
    }

    debug_assert_eq!(halo_g.len(), 1);
    let halo_g = halo_g[0].to_affine();

    debug_assert_eq!(halo_a.len(), 1);
    debug_assert_eq!(halo_b.len(), 1);
    let schnorr_proof = schnorr_protocol(
        halo_a[0], halo_b[0], halo_g, randomness, u_prime, pedersen_h, challenger,
    );

    Ok(OpeningProof {
        halo_g,
        halo_l: ProjectivePoint::batch_to_affine(&halo_l),
        halo_r: ProjectivePoint::batch_to_affine(&halo_r),
        schnorr_proof,
    })
}

fn build_halo_b<C: Curve>(
    points: &[C::ScalarField],
    v: C::ScalarField,
    degree: usize,
) -> Vec<C::ScalarField> {
    let power_points = points
        .iter()
        .map(|&p| powers(p, degree))
        .collect::<Vec<_>>();
    (0..degree)
        .map(|i| reduce_with_powers(&power_points.iter().map(|v| v[i]).collect::<Vec<_>>(), v))
        .collect()
}

fn schnorr_protocol<C: HaloCurve>(
    halo_a: C::ScalarField,
    halo_b: C::ScalarField,
    halo_g: AffinePoint<C>,
    randomness: C::ScalarField,
    u_curve: ProjectivePoint<C>,
    pedersen_h: ProjectivePoint<C>,
    challenger: &mut Challenger<C::BaseField>,
) -> SchnorrProof<C> {
    let (d, s) = (C::ScalarField::rand(), C::ScalarField::rand());
    let r_curve = C::convert(d) * (halo_g.to_projective() + C::convert(halo_b) * u_curve)
        + C::convert(s) * pedersen_h;

    challenger.observe_proj_point(r_curve);
    let chall_bf = challenger.get_challenge();
    let chall = chall_bf
        .try_convert::<C::ScalarField>()
        .expect("Improbable");
    let z1 = halo_a * chall + d;
    let z2 = randomness * chall + s;
    SchnorrProof {
        r: r_curve.to_affine(),
        z1,
        z2,
    }
}

/// Verify the final IPA.
#[allow(clippy::too_many_arguments)]
pub fn verify_ipa<C: HaloCurve>(
    halo_l: &[AffinePoint<C>],
    halo_r: &[AffinePoint<C>],
    halo_g: AffinePoint<C>,
    commitment: ProjectivePoint<C>,
    value: C::ScalarField,
    halo_b: C::ScalarField,
    halo_us: &[C::ScalarField],
    u_prime: ProjectivePoint<C>,
    pedersen_h: AffinePoint<C>,
    schnorr_challenge: C::ScalarField,
    schnorr_proof: SchnorrProof<C>,
) -> bool {
    // Now we begin IPA verification by computing P' and u' as in Protocol 1 of Bulletproofs.
    // In Protocol 1 we compute u' = [x] u, but we leverage to endomorphism, instead computing
    // u' = [n(x)] u.

    // Compute [c] [n(x)] u = [c] u'.
    let u_n_x_c = C::convert(value) * u_prime;
    let p_prime = commitment + u_n_x_c;

    // Compute Q as defined in the Halo paper.
    let mut points = halo_l.to_vec();
    points.extend(halo_r.iter());
    let mut scalars = halo_us.iter().map(|u| u.square()).collect::<Vec<_>>();
    scalars.extend(
        halo_us
            .iter()
            .map(|chal| chal.multiplicative_inverse_assuming_nonzero().square()),
    );
    let precomputation = msm_precompute(&AffinePoint::batch_to_projective(&points), 8);
    let q = msm_execute_parallel(&precomputation, &scalars) + p_prime;

    // Performing ZK opening protocol.
    C::convert(schnorr_challenge) * q + schnorr_proof.r
        == C::convert(schnorr_proof.z1) * (halo_g.to_projective() + C::convert(halo_b) * u_prime)
            + C::convert(schnorr_proof.z2) * pedersen_h.to_projective()
}
