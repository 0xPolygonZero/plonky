#![allow(clippy::many_single_char_names)]
#![allow(clippy::too_many_arguments)]
use crate::openings::open_all_polynomials;
use crate::proof::PlookupProof;
use anyhow::Result;
use plonky::halo::batch_opening_proof;
use plonky::plonk_challenger::Challenger;
use plonky::plonk_util::reduce_with_powers;
use plonky::polynomial::Polynomial;
use plonky::util::log2_strict;
use plonky::{blake_hash_usize_to_curve, fft_precompute, msm_precompute, AffinePoint, Field, HaloCurve};

pub const SECURITY_BITS: usize = 128;

/// Computes a proof that `t` is a sub(multi)set of `f`, following the Plookup protocol (https://ia.cr/2020/315).
pub fn prove<C: HaloCurve>(f: &[C::ScalarField], t: &[C::ScalarField]) -> Result<PlookupProof<C>> {
    let (n, f, t) = pad_inputs(f, t);

    // Compute the `s` multiset of the Plookup protocol.
    let mut s = f.clone();
    s.extend_from_slice(&t);
    sort_by(&mut s, &t);

    // Initiate the Fiat-Shamir challenger.
    let mut challenger = Challenger::new(SECURITY_BITS);

    // FFT precomputation on the cyclic subgroup of order `n+1`.
    let fft_precomputation = fft_precompute(n + 1);

    // Compute the polynomials corresponding to `f`, `t`, `h1` and `h2`.
    let f_padded = padded(&f, n + 1);
    let f_poly = Polynomial::from_evaluations(&f_padded, &fft_precomputation);
    let t_poly = Polynomial::from_evaluations(&t, &fft_precomputation);
    let h1_poly = Polynomial::from_evaluations(&s[..n + 1], &fft_precomputation);
    let h2_poly = Polynomial::from_evaluations(&s[n..], &fft_precomputation);

    // Curve points used in the IPA.
    let gs = (0..2 * n + 2)
        .map(blake_hash_usize_to_curve::<C>)
        .collect::<Vec<_>>();
    let h = blake_hash_usize_to_curve(2 * n + 2);
    let u_curve = blake_hash_usize_to_curve(2 * n + 3);
    let msm_precomputation = msm_precompute(&AffinePoint::batch_to_projective(&gs[..n + 1]), 8);
    let msm_precomputation_2n2 = msm_precompute(&AffinePoint::batch_to_projective(&gs), 8);

    // Commit to all polynomials.
    let c_f = f_poly.commit(&msm_precomputation, h, true);
    let c_t = t_poly.commit(&msm_precomputation, h, false);
    let c_h1 = h1_poly.commit(&msm_precomputation, h, false);
    let c_h2 = h2_poly.commit(&msm_precomputation, h, false);

    // Observe the commitments to get verifier challenges.
    // `beta` and `gamma` are used to construct the Plookup grand product.
    challenger.observe_affine_points(&[
        c_f.to_affine(),
        c_t.to_affine(),
        c_h1.to_affine(),
        c_h2.to_affine(),
    ]);
    let (beta_bf, gamma_bf) = challenger.get_2_challenges();
    let beta = C::BaseField::try_convert(&beta_bf).unwrap();
    let gamma = C::BaseField::try_convert(&gamma_bf).unwrap();

    // Compute and commit to the `z` grand product polynomial of Plookup.
    // `f` is a subset of `t` iff this polynomial is well-formed and its last value is 1.
    let z_values = grand_polynomial(&f, &t, &s, beta, gamma);
    let z_poly = Polynomial::from_evaluations(&z_values, &fft_precomputation);
    let c_z = z_poly.commit(&msm_precomputation, h, true);

    // Observe the commitment to get a verifier challenge.
    // `alpha` is used to batch all vanishing polynomials.
    challenger.observe_affine_point(c_z.to_affine());
    let alpha_bf = challenger.get_challenge();
    let alpha = C::BaseField::try_convert(&alpha_bf).unwrap();

    // Compute the coefficients of the "vanishing polynomial".
    // `f` is a subset of `t` iff this polynomial vanishes on subgroup `H`.
    let vanishing_poly = vanishing_polynomial(
        &f_poly, &t_poly, &h1_poly, &h2_poly, &z_poly, beta, gamma, alpha, n,
    );

    // Commit to the quotient of the vanishing polynomial by `x^(n+1) - 1`. The quotient exists
    // because the vanishing polynomial vanishes on `H`.
    let mut quotient_poly = vanishing_poly.divide_by_z_h(n + 1);
    quotient_poly.trim();
    assert!(quotient_poly.len() <= 2 * n + 1);
    quotient_poly.pad(2 * n + 2);
    let c_quotient = quotient_poly.commit(&msm_precomputation_2n2, h, true);

    // Observe the commitment to get a verifier challenge.
    // `zeta` is the point at which we'll open all polynomials.
    challenger.observe_affine_point(c_quotient.to_affine());
    let zeta_bf = challenger.get_challenge();
    let zeta: C::ScalarField = C::BaseField::try_convert(&zeta_bf).unwrap();

    // Open all polynomials.
    let generator = C::ScalarField::primitive_root_of_unity(log2_strict(n + 1));
    let openings = open_all_polynomials(
        &f_poly,
        &t_poly,
        &h1_poly,
        &h2_poly,
        &z_poly,
        &quotient_poly,
        zeta,
        generator,
    );

    // Observe the openings to get verifier challenges.
    // `v`, `u` and `u_scaling` are used in the IPA proof.
    let openings_bf: Vec<_> = openings
        .to_vec()
        .into_iter()
        .map(|f| {
            C::try_convert_s2b(f)
                .expect("For now, we assume that all opened values fit in both fields")
        })
        .collect();
    challenger.observe_elements(&openings_bf);
    let (v_bf, u_bf, u_scaling_bf) = challenger.get_3_challenges();
    let v = v_bf.try_convert::<C::ScalarField>()?;
    let u = u_bf.try_convert::<C::ScalarField>()?;
    let u_scaling = u_scaling_bf.try_convert::<C::ScalarField>()?;

    let commitments = vec![c_f, c_t, c_h1, c_h2, c_z, c_quotient];

    let coeffs = vec![
        padded(&f_poly[..], 2 * n + 2),
        padded(&t_poly[..], 2 * n + 2),
        padded(&h1_poly[..], 2 * n + 2),
        padded(&h2_poly[..], 2 * n + 2),
        padded(&z_poly[..], 2 * n + 2),
        padded(&quotient_poly[..], 2 * n + 2),
    ];

    // Compute the Halo opening proof.
    let halo_proof = batch_opening_proof(
        &coeffs.iter().map(|c| c.as_slice()).collect::<Vec<_>>(),
        &commitments,
        &[zeta, zeta * generator],
        &gs,
        h.to_projective(),
        u_curve,
        u,
        v,
        u_scaling,
        2 * n + 2,
        SECURITY_BITS,
        &mut challenger,
    )?;

    Ok(PlookupProof::from((commitments, openings, halo_proof, n)))
}

fn pad_inputs<F: Field>(f: &[F], t: &[F]) -> (usize, Vec<F>, Vec<F>) {
    let d = t.len();
    let f = if f.len() + 1 < d {
        padded(f, d - 1)
    } else {
        f.to_vec()
    };
    let n = f.len().next_power_of_two() - 1;

    let f = padded(&f, n);
    let t = padded(t, n + 1);
    (n, f, t)
}

/// Sorts a vector `f` by another `t` in-place. `f` is sorted by `t` if `f` is a subset of `t`,
/// and elements of `f` appear in the same order as elements of `t`.
fn sort_by<F: Field>(f: &mut [F], t: &[F]) {
    f.sort_by(|a, b| {
        let i = t.iter().position(|x| x == a).unwrap();
        let j = t.iter().position(|x| x == b).unwrap();
        i.cmp(&j)
    });
}

/// Computes the Plookup grand product polynomial.
fn grand_polynomial<F: Field>(f: &[F], t: &[F], s: &[F], beta: F, gamma: F) -> Vec<F> {
    let n = f.len();
    let mut values = vec![F::ONE];

    let beta1 = beta + F::ONE;
    let gamma_beta1 = gamma * beta1;
    let mut beta1_pow = beta1;
    let mut prod_a = gamma + f[0];
    let mut prod_b = gamma_beta1 + t[0] + beta * t[1];
    let mut prod_c = (gamma_beta1 + s[0] + beta * s[1]) * (gamma_beta1 + s[n] + beta * s[n + 1]);
    for i in 1..n {
        values.push((beta1_pow * prod_a * prod_b) / prod_c);

        beta1_pow = beta1_pow * beta1;
        prod_a = prod_a * (gamma + f[i]);
        prod_b = prod_b * (gamma_beta1 + t[i] + beta * t[i + 1]);
        prod_c = prod_c
            * (gamma_beta1 + s[i] + beta * s[i + 1])
            * (gamma_beta1 + s[n + i] + beta * s[n + i + 1]);
    }
    values.push(F::ONE);
    values
}

/// Computes the Plookup vanishing polynomial.
fn vanishing_polynomial<F: Field>(
    f_poly: &Polynomial<F>,
    t_poly: &Polynomial<F>,
    h1_poly: &Polynomial<F>,
    h2_poly: &Polynomial<F>,
    z_poly: &Polynomial<F>,
    beta: F,
    gamma: F,
    alpha: F,
    n: usize,
) -> Polynomial<F> {
    let order = 4 * (n + 1);
    let generator_4 = F::primitive_root_of_unity(log2_strict(order));
    let fft_precomp4 = fft_precompute(order);
    let z_4_values = z_poly.eval_domain(&fft_precomp4);
    let f_4_values = f_poly.eval_domain(&fft_precomp4);
    let t_4_values = t_poly.eval_domain(&fft_precomp4);
    let h1_4_values = h1_poly.eval_domain(&fft_precomp4);
    let h2_4_values = h2_poly.eval_domain(&fft_precomp4);

    let beta1 = beta + F::ONE;
    let gamma_beta1 = gamma * beta1;

    let values = F::cyclic_subgroup_known_order(generator_4, order)
        .into_iter()
        .enumerate()
        .map(|(i, x)| {
            let vanishing_z1_term =
                eval_l_i(n + 1, 0, generator_4.exp_usize(4), x) * (z_4_values[i] - F::ONE);
            let next = (i + 4) % order;
            let vanishing_shift_term = (x - generator_4.exp_usize(4 * n))
                * z_4_values[i]
                * beta1
                * (gamma + f_4_values[i])
                * (gamma_beta1 + t_4_values[i] + beta * t_4_values[next])
                - (x - generator_4.exp_usize(4 * n))
                    * z_4_values[next]
                    * (gamma_beta1 + h1_4_values[i] + beta * h1_4_values[next])
                    * (gamma_beta1 + h2_4_values[i] + beta * h2_4_values[next]);
            let eval_last = eval_l_i(n + 1, n, generator_4.exp_usize(4), x);
            let vanishing_hs_term = eval_last * (h1_4_values[i] - h2_4_values[next]);
            let vanishing_last_term = eval_last * (z_4_values[i] - F::ONE);

            if cfg!(debug_assertions) && i % 4 == 0 {
                for &x in &[
                    vanishing_z1_term,
                    vanishing_shift_term,
                    vanishing_hs_term,
                    vanishing_last_term,
                ] {
                    assert_eq!(x, F::ZERO);
                }
            }

            reduce_with_powers(
                &[
                    vanishing_z1_term,
                    vanishing_shift_term,
                    vanishing_hs_term,
                    vanishing_last_term,
                ],
                alpha,
            )
        })
        .collect::<Vec<_>>();
    Polynomial::from_evaluations(&values, &fft_precomp4)
}

/// Computes the evaluation of the `i`-th Lagrange basis polynomial of order `n` on a point `x`.
/// Uses the formula `L_i(x) = w^i(x^n - 1) / n(x-w^i)`.
pub fn eval_l_i<F: Field>(n: usize, i: usize, generator: F, x: F) -> F {
    let g = generator.exp_usize(i);
    if x == g {
        F::ZERO
    } else {
        g * (x.exp_usize(n) - F::ONE) / (F::from_canonical_usize(n) * (x - g))
    }
}

/// Zero-pad a slice to have length `n`.
pub fn padded<F: Field>(s: &[F], n: usize) -> Vec<F> {
    let mut pad = s.to_vec();
    pad.extend((s.len()..n).map(|_| F::ZERO));
    pad
}

#[cfg(test)]
mod tests {
    use crate::plookup::sort_by;
    use plonky::{Field, TweedledeeBase};

    #[test]
    fn test_sort_by() {
        type F = TweedledeeBase;
        let mut f = vec![F::FIVE, F::TWO, F::ONE];
        let t = vec![F::ONE, F::TWO, F::THREE, F::FOUR, F::FIVE];
        sort_by(&mut f, &t);
        assert_eq!(f, vec![F::ONE, F::TWO, F::FIVE]);
    }
}
