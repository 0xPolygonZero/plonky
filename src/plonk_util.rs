use crate::partition::get_subgroup_shift;
use crate::witness::Witness;
use crate::{ifft_with_precomputation_power_of_2, msm_execute_parallel, AffinePoint, CircuitBuilder, Curve, FftPrecomputation, Field, HaloCurve, MsmPrecomputation, Polynomial, PolynomialCommitment, ProjectivePoint, Target, NUM_ROUTED_WIRES};
use rayon::prelude::*;

/// Evaluate the polynomial which vanishes on any multiplicative subgroup of a given order `n`.
pub(crate) fn eval_zero_poly<F: Field>(n: usize, x: F) -> F {
    // Z(x) = x^n - 1
    x.exp_usize(n) - F::ONE
}

/// Evaluate the Lagrange basis `L_1` with `L_1(1) = 1`, and `L_1(x) = 0` for other members of an
/// order `n` multiplicative subgroup.
pub(crate) fn eval_l_1<F: Field>(n: usize, x: F) -> F {
    if x.is_one() {
        // The code below would divide by zero, since we have (x - 1) in both the numerator and
        // denominator.
        return F::ONE;
    }

    // L_1(x) = (x^n - 1) / (n * (x - 1))
    //        = Z(x) / (n * (x - 1))
    eval_zero_poly(n, x) / (F::from_canonical_usize(n) * (x - F::ONE))
}

/// Computes a sum of terms weighted by powers of alpha.
pub fn reduce_with_powers<F: Field>(terms: &[F], alpha: F) -> F {
    let mut sum = F::ZERO;
    for &term in terms.iter().rev() {
        sum = sum * alpha + term;
    }
    sum
}

/// Computes a sum of terms weighted by powers of alpha.
pub(crate) fn reduce_with_powers_recursive<C: HaloCurve>(
    builder: &mut CircuitBuilder<C>,
    terms: &[Target<C::ScalarField>],
    alpha: Target<C::ScalarField>,
) -> Target<C::ScalarField> {
    let mut sum = builder.zero_wire();
    for &term in terms.iter().rev() {
        sum = builder.mul_add(sum, alpha, term);
    }
    sum
}

/// Compute `n(x)` for a given `x`, where `n` is the injective function related to the Halo
/// endomorphism.
pub fn halo_n<C: HaloCurve>(s_bits: &[bool]) -> C::ScalarField {
    // This is based on Algorithm 2 of the Halo paper, except that we start with (a, b) = (0, 0).

    debug_assert_eq!(s_bits.len() % 2, 0, "Number of scalar bits must be even");

    let zero = C::ScalarField::ZERO;
    let mut a = zero;
    let mut b = zero;

    for s_bits_chunk in s_bits.chunks(2) {
        let bit_lo = s_bits_chunk[0];
        let bit_hi = s_bits_chunk[1];

        let sign = if bit_lo {
            C::ScalarField::ONE
        } else {
            C::ScalarField::NEG_ONE
        };
        let (c, d) = if bit_hi { (sign, zero) } else { (zero, sign) };

        a = a.double() + c;
        b = b.double() + d;
    }

    a * C::ZETA_SCALAR + b
}

/// Compute `[n(s)].P` for a given `s`, where `n` is the injective function related to the Halo
/// endomorphism.
pub fn halo_n_mul<C: HaloCurve>(s_bits: &[bool], p: AffinePoint<C>) -> AffinePoint<C> {
    // This is based on Algorithm 1 of the Halo paper, except that we start with Acc = O.

    debug_assert_eq!(s_bits.len() % 2, 0, "Number of scalar bits must be even");

    let p_p = p.to_projective();
    let p_n = -p_p;
    let endo_p_p = p.endomorphism().to_projective();
    let endo_p_n = -endo_p_p;

    let mut acc = ProjectivePoint::<C>::ZERO;

    for s_bits_chunk in s_bits.chunks(2) {
        let bit_lo = s_bits_chunk[0];
        let bit_hi = s_bits_chunk[1];

        let s = if bit_hi {
            if bit_lo {
                endo_p_p
            } else {
                endo_p_n
            }
        } else if bit_lo {
            p_p
        } else {
            p_n
        };
        acc = acc.double() + s;
    }

    acc.to_affine()
}

pub fn eval_poly<F: Field>(coeffs: &[F], x: F) -> F {
    let mut ans = F::ZERO;
    let mut x_pow = F::ONE;
    for &c in coeffs {
        ans = ans + (c * x_pow);
        x_pow = x_pow * x;
    }
    ans
}

/// Compute `[x^0, x^1, ..., x^(n - 1)]`.
pub fn powers<F: Field>(x: F, n: usize) -> Vec<F> {
    let mut powers = Vec::new();
    let mut current = F::ONE;
    for i in 0..n {
        if i != 0 {
            current = current * x;
        }
        powers.push(current);
    }
    powers
}

/// Compute `[x^0, x^1, ..., x^(n - 1)]`.
pub(crate) fn powers_recursive<C: HaloCurve>(
    builder: &mut CircuitBuilder<C>,
    x: Target<C::ScalarField>,
    n: usize,
) -> Vec<Target<C::ScalarField>> {
    let mut powers = Vec::new();
    let mut current = builder.one_wire();
    for i in 0..n {
        if i != 0 {
            current = builder.mul(current, x);
        }
        powers.push(current);
    }
    powers
}

/// Returns the evaluation of a list of polynomials at a point.
pub(crate) fn eval_polys<F: Field>(polys: &[Polynomial<F>], powers: &[F]) -> Vec<F> {
    polys.iter().map(|p| p.eval_from_power(powers)).collect()
}

/// Zero-pad a list of `n` polynomial coefficients to a length of `8n`, which is the degree at
/// which we do most polynomial arithmetic.
pub(crate) fn pad_to_8n<F: Field>(coeffs: &[F]) -> Vec<F> {
    let n = coeffs.len();

    let mut result = coeffs.to_vec();
    while result.len() < 8 * n {
        result.push(F::ZERO);
    }
    result
}

pub(crate) fn values_to_polynomials<F: Field>(
    values_vec: &[Vec<F>],
    fft_precomputation: &FftPrecomputation<F>,
) -> Vec<Polynomial<F>> {
    values_vec
        .par_iter()
        .map(|values| Polynomial::from_evaluations(&values, fft_precomputation))
        .collect()
}

pub(crate) fn polynomials_to_values_padded<F: Field>(
    polys_vec: &[Polynomial<F>],
    fft_precomputation: &FftPrecomputation<F>,
) -> Vec<Vec<F>> {
    polys_vec
        .par_iter()
        .map(|poly| {
            let padded_poly = poly.padded(poly.len() * 8);
            padded_poly.eval_domain(fft_precomputation)
        })
        .collect()
}

/// Like `pedersen_commit`, but with no blinding factor.
pub fn pedersen_hash<C: Curve>(
    xs: &[C::ScalarField],
    pedersen_g_msm_precomputation: &MsmPrecomputation<C>,
) -> ProjectivePoint<C> {
    msm_execute_parallel(pedersen_g_msm_precomputation, xs)
}

#[allow(dead_code)]
fn pedersen_commit<C: Curve>(
    xs: &[C::ScalarField],
    opening: C::ScalarField,
    h: AffinePoint<C>,
    pedersen_g_msm_precomputation: &MsmPrecomputation<C>,
) -> ProjectivePoint<C> {
    // TODO: Couldn't get this working with *.
    let h = h.to_projective();
    let mul_precomputation = h.mul_precompute();
    let blinding_term = h.mul_with_precomputation(opening, mul_precomputation);

    msm_execute_parallel(pedersen_g_msm_precomputation, xs) + blinding_term
}

pub fn commit_polynomials<C: Curve>(
    polynomials: &[Polynomial<C::ScalarField>],
    msm_precomputation: &MsmPrecomputation<C>,
    blinding_point: AffinePoint<C>,
    blinding: bool,
) -> Vec<PolynomialCommitment<C>> {
    PolynomialCommitment::coeffs_vec_to_commitments(
        polynomials
            .iter()
            .map(|p| p.coeffs())
            .collect::<Vec<_>>()
            .as_slice(),
        msm_precomputation,
        blinding_point,
        blinding,
    )
}

// Generate Z, which is used in Plonk's permutation argument.
pub fn permutation_polynomial<F: Field>(
    degree: usize,
    subgroup: &[F],
    witness: &Witness<F>,
    sigma_values: &[Vec<F>],
    beta: F,
    gamma: F,
) -> Vec<F> {
    let mut plonk_z_points = vec![F::ONE];
    let k_is = (0..NUM_ROUTED_WIRES)
        .map(get_subgroup_shift::<F>)
        .collect::<Vec<_>>();
    for i in 1..degree {
        let x = subgroup[i - 1];
        let mut numerator = F::ONE;
        let mut denominator = F::ONE;
        for j in 0..NUM_ROUTED_WIRES {
            let wire_value = witness.get_indices(i - 1, j);
            let k_i = k_is[j];
            let s_id = k_i * x;
            let s_sigma = sigma_values[j][8 * (i - 1)];
            numerator = numerator * (wire_value + beta * s_id + gamma);
            denominator = denominator * (wire_value + beta * s_sigma + gamma);
        }
        let last = *plonk_z_points.last().unwrap();
        plonk_z_points.push(last * numerator / denominator);
    }
    plonk_z_points
}

pub fn sigma_polynomials<F: Field>(
    sigma: Vec<usize>,
    degree: usize,
    subgroup_generator: F,
) -> Vec<Vec<F>> {
    sigma
        .chunks(degree)
        .map(|chunk| {
            chunk
                .par_iter()
                .map(|&x| {
                    get_subgroup_shift::<F>(x / degree) * subgroup_generator.exp_usize(x % degree)
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

/// Given polynomials `[p_0,...,p_k]` of degree `degree` and `alpha \in F`, returns `\sum_{i=0}^k alpha^i p_i`.
pub(crate) fn scale_polynomials<F: Field>(
    polynomials: Vec<Polynomial<F>>,
    alpha: F,
    degree: usize,
) -> Polynomial<F> {
    let alpha_powers = powers(alpha, polynomials.len());
    Polynomial::from(
        (0..degree)
            .map(|i| {
                (0..polynomials.len())
                    .map(|j| polynomials[j][i] * alpha_powers[j])
                    .fold(F::ZERO, |acc, x| acc + x)
            })
            .collect::<Vec<_>>(),
    )
}

#[allow(dead_code)]
pub(crate) fn polynomial_degree_plus_1<F: Field>(
    points: &[F],
    fft_precomputation: &FftPrecomputation<F>,
) -> usize {
    let coeffs = ifft_with_precomputation_power_of_2(&points, fft_precomputation);
    coeffs.iter().rev().skip_while(|c| c.is_zero()).count()
}

// TODO: Maybe a streaming version using an `Iterator` would be faster and wouldn't require as much memory for large circuits.
// TODO: Optimize this.
pub fn halo_s<F: Field>(us: &[F]) -> Vec<F> {
    let n = 1 << us.len();
    let mut res = vec![F::ONE; n];
    let us_inv = F::batch_multiplicative_inverse(us);

    for (j, (&u, &u_inv)) in us.iter().rev().zip(us_inv.iter().rev()).enumerate() {
        for (i, x) in res.iter_mut().enumerate() {
            if i & (1 << j) == 0 {
                *x = *x * u_inv;
            } else {
                *x = *x * u;
            }
        }
    }
    res
}

/// Evaluate `g(X, {u_i})` as defined in the Halo paper.
pub fn halo_g<F: Field>(x: F, us: &[F]) -> F {
    let mut product = F::ONE;
    let mut x_power = x;
    for &u_i in us.iter().rev() {
        let u_i_inv = u_i.multiplicative_inverse_assuming_nonzero();
        let term = u_i * x_power + u_i_inv;
        product = product * term;
        x_power = x_power.square();
    }
    product
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{CircuitBuilder, Curve, Field, PartialWitness, Tweedledee};

    #[test]
    fn test_halo_n() {
        type C = Tweedledee;
        type SF = <Tweedledee as Curve>::ScalarField;
        let p = C::convert(SF::rand()) * C::GENERATOR_PROJECTIVE;
        let r = SF::rand();
        let res = C::convert(halo_n::<C>(&r.to_canonical_bool_vec()[..128])) * p;
        let p = p.to_affine();
        assert_eq!(
            res.to_affine(),
            halo_n_mul::<C>(&r.to_canonical_bool_vec()[..128], p)
        )
    }

    #[test]
    fn test_permutation_polynomial() {
        let mut builder = CircuitBuilder::<Tweedledee>::new(128);
        let one = builder.one_wire();
        let t = builder.add_virtual_target();
        let t_sq = builder.square(t);
        let quad = builder.add_many(&[one, t, t_sq]);
        let seven =
            builder.constant_wire(<Tweedledee as Curve>::ScalarField::from_canonical_usize(7));
        let res = builder.sub(quad, seven);
        builder.assert_zero(res);
        let mut partial_witness = PartialWitness::new();
        partial_witness.set_target(t, <Tweedledee as Curve>::ScalarField::TWO);
        let circuit = builder.build();
        let witness = circuit.generate_witness(partial_witness);
        let beta = <Tweedledee as Curve>::ScalarField::rand();
        let gamma = <Tweedledee as Curve>::ScalarField::rand();
        let plonk_z_points_n = permutation_polynomial(
            circuit.degree(),
            &circuit.subgroup_n,
            &witness,
            &circuit.s_sigma_values_8n,
            beta,
            gamma,
        );
        // Verify that the permutation polynomial is well-formed.
        let k_is = (0..NUM_ROUTED_WIRES)
            .map(get_subgroup_shift::<<Tweedledee as Curve>::ScalarField>)
            .collect::<Vec<_>>();
        let wire_values = &witness.transpose();
        for (i, &x) in circuit.subgroup_n.iter().enumerate() {
            let (z_x, z_gz) = (
                plonk_z_points_n[i],
                plonk_z_points_n[(i + 1) % circuit.degree()],
            );
            let mut f_prime = <Tweedledee as Curve>::ScalarField::ONE;
            let mut g_prime = <Tweedledee as Curve>::ScalarField::ONE;
            for j in 0..NUM_ROUTED_WIRES {
                let wire_value = wire_values[j][i];
                let k_i = k_is[j];
                let s_id = k_i * x;
                let s_sigma = circuit.s_sigma_values_8n[j][8 * i];
                f_prime = f_prime * (wire_value + beta * s_id + gamma);
                g_prime = g_prime * (wire_value + beta * s_sigma + gamma);
            }
            let vanishing_v_shift_term = f_prime * z_x - g_prime * z_gz;
            assert_eq!(
                vanishing_v_shift_term,
                <Tweedledee as Curve>::ScalarField::ZERO
            );
        }
    }

    #[test]
    fn test_s_vector_g_function() {
        type F = <Tweedledee as Curve>::ScalarField;
        let us = (0..10).map(|_| F::rand()).collect::<Vec<_>>();
        let x = F::rand();
        assert_eq!(
            F::inner_product(&halo_s(&us), &powers(x, 1 << 10)),
            halo_g(x, &us)
        );
    }
}
