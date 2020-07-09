use crate::partition::get_subgroup_shift;
use crate::witness::Witness;
use crate::{divide_by_z_h, fft_with_precomputation_power_of_2, ifft_with_precomputation_power_of_2, msm_execute_parallel, msm_parallel, polynomial_division, AffinePoint, CircuitBuilder, Curve, FftPrecomputation, Field, HaloCurve, MsmPrecomputation, ProjectivePoint, Target, NUM_ROUTED_WIRES};
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
pub(crate) fn reduce_with_powers<F: Field>(terms: &[F], alpha: F) -> F {
    let mut sum = F::ZERO;
    for &term in terms.iter().rev() {
        sum = sum * alpha + term;
    }
    sum
}

/// Computes a sum of terms weighted by powers of alpha.
pub(crate) fn reduce_with_powers_recursive<C: HaloCurve>(
    builder: &mut CircuitBuilder<C>,
    terms: &[Target],
    alpha: Target,
) -> Target {
    let mut sum = builder.zero_wire();
    for &term in terms.iter().rev() {
        sum = builder.mul_add(sum, alpha, term);
    }
    sum
}

/// Compute `n(x)` for a given `x`, where `n` is the injective function related to the Halo
/// endomorphism.
pub(crate) fn halo_n<C: HaloCurve>(s_bits: &[bool]) -> C::ScalarField {
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
pub(crate) fn halo_n_mul<C: HaloCurve>(s_bits: &[bool], p: AffinePoint<C>) -> AffinePoint<C> {
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
        } else {
            if bit_lo {
                p_p
            } else {
                p_n
            }
        };
        acc = acc.double() + s;
    }

    acc.to_affine()
}

pub(crate) fn eval_poly<F: Field>(coeffs: &[F], x: F) -> F {
    let mut ans = F::ZERO;
    let mut x_pow = F::ONE;
    for &c in coeffs {
        ans = ans + (c * x_pow);
        x_pow = x_pow * x;
    }
    ans
}

/// Compute `[x^0, x^1, ..., x^(n - 1)]`.
pub(crate) fn powers<F: Field>(x: F, n: usize) -> Vec<F> {
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
    x: Target,
    n: usize,
) -> Vec<Target> {
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
pub(crate) fn eval_coeffs<F: Field>(coeffs: &[Vec<F>], powers: &[F]) -> Vec<F> {
    coeffs
        .iter()
        .map(|c| F::inner_product(c, &powers))
        .collect()
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

pub(crate) fn values_to_coeffs<F: Field>(
    values_vec: &[Vec<F>],
    fft_precomputation: &FftPrecomputation<F>,
) -> Vec<Vec<F>> {
    values_vec
        .iter()
        .map(|values| ifft_with_precomputation_power_of_2(values, fft_precomputation))
        .collect()
}

pub(crate) fn coeffs_to_values<F: Field>(
    coefficients_vec: &[Vec<F>],
    fft_precomputation: &FftPrecomputation<F>,
) -> Vec<Vec<F>> {
    coefficients_vec
        .iter()
        .map(|coeffs| fft_with_precomputation_power_of_2(coeffs, fft_precomputation))
        .collect()
}

pub(crate) fn coeffs_to_values_padded<F: Field>(
    coefficients_vec: &[Vec<F>],
    fft_precomputation: &FftPrecomputation<F>,
) -> Vec<Vec<F>> {
    coefficients_vec
        .iter()
        .map(|coeffs| fft_with_precomputation_power_of_2(&pad_to_8n(coeffs), fft_precomputation))
        .collect()
}

/// Like `pedersen_commit`, but with no blinding factor.
pub fn pedersen_hash<C: Curve>(
    xs: &[C::ScalarField],
    pedersen_g_msm_precomputation: &MsmPrecomputation<C>,
) -> ProjectivePoint<C> {
    msm_execute_parallel(pedersen_g_msm_precomputation, xs)
}

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
        .map(|j| get_subgroup_shift::<F>(j))
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

/// Computes the commitments (without randomness) `Comm(L_i)` for `i=0..num`, where `L_i` is the `i-th`
/// Lagrange basis polynomial on the subgroup of order `order` of `C`'s scalar field,
/// i.e., `L_i(g^j) = delta_ij`, with `g` a `order`-th primitive root of unity.
pub fn precompute_lagrange_commitments<C: Curve>(
    order_log: usize,
    num: usize,
    msm_precomputation: &MsmPrecomputation<C>,
) -> Vec<ProjectivePoint<C>> {
    let order = 1 << order_log;
    let g = C::ScalarField::primitive_root_of_unity(order_log);
    let mut xn_minus_one = vec![C::ScalarField::ZERO; order + 1];
    xn_minus_one[order] = C::ScalarField::ONE;
    xn_minus_one[0] = C::ScalarField::NEG_ONE;

    (0..num)
        .into_par_iter()
        .map(|i| {
            let a = g.exp_u32(2 * i as u32);
            // We divide by X - a.
            let denom = vec![-a, C::ScalarField::ONE];
            let mut d = polynomial_division(&xn_minus_one, &denom);
            debug_assert_eq!(d.1, vec![C::ScalarField::ZERO]);
            let factor = C::ScalarField::from_canonical_usize(order)
                .multiplicative_inverse_assuming_nonzero()
                * a;
            d.0.iter_mut().for_each(|x| {
                *x = *x * factor;
            });
            msm_execute_parallel(msm_precomputation, &d.0)
        })
        .collect()
}

/// Compute the commitment of the wire polynomials on the `PublicInputGate`s.
/// Returns a vector of length `NUM_WIRES` where the `i-th` value is the commitment of
/// the polynomial interpolating the `i-th` wire of the `PublicInputGate`s and then having
/// zero values on the remaining gates.
/// If `step_by_two` is true, we assume that the public input gates are at indices
/// `0, 2, ..., 2*num_public_gates-2`. Otherwise, we assume they are at indices
/// `0,1,...,num_public_gates-1`.
pub fn pis_commitments<C: Curve>(
    wires: &[Vec<C::ScalarField>],
    precomputed_lagrange_commitments: &[ProjectivePoint<C>],
    num_public_gates: usize,
    step_by_two: bool,
) -> Vec<ProjectivePoint<C>> {
    wires
        .par_iter()
        .map(|wire_vec| {
            msm_parallel(
                &(if step_by_two {
                    (0..num_public_gates)
                        .map(|i| wire_vec[2 * i])
                        .collect::<Vec<_>>()
                } else {
                    wire_vec[..num_public_gates].to_vec()
                }),
                precomputed_lagrange_commitments,
                8,
            )
        })
        .collect()
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::util::log2_ceil;
    use crate::{blake_hash_usize_to_curve, fft_precompute, msm_precompute, Circuit, CircuitBuilder, Curve, Field, PartialWitness, Tweedledee, NUM_WIRES};
    use std::time::Instant;

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
            .map(|j| get_subgroup_shift::<<Tweedledee as Curve>::ScalarField>(j))
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

    fn manually_compute_pis_commitments<C: Curve>(
        wires: &[C::ScalarField],
        msm_precomputation: &MsmPrecomputation<C>,
        fft_precomputation: &FftPrecomputation<C::ScalarField>,
        order: usize,
    ) -> ProjectivePoint<C> {
        let mut padded_wires = vec![C::ScalarField::ZERO; order];
        (0..wires.len()).for_each(|i| {
            padded_wires[2 * i] = wires[i];
        });
        let coeffs = ifft_with_precomputation_power_of_2(&padded_wires, fft_precomputation);
        msm_execute_parallel(msm_precomputation, &coeffs)
    }

    #[test]
    fn test_public_input_commitments() {
        let noww = Instant::now();
        type C = Tweedledee;
        type F = <C as Curve>::ScalarField;
        let order_log = 14;
        let order = 1 << order_log;
        let num = 10;
        let now = Instant::now();
        let pedersen_g: Vec<_> = (0..order).map(blake_hash_usize_to_curve::<C>).collect();
        let w = 11;
        let msm_precomputation = msm_precompute(&AffinePoint::batch_to_projective(&pedersen_g), w);
        let fft_precomputation = fft_precompute(order);
        dbg!(now.elapsed());

        let now = Instant::now();
        let precomputed_lagrange_commitments =
            precompute_lagrange_commitments(order_log, num, &msm_precomputation);
        dbg!(now.elapsed());

        let wires = (0..NUM_WIRES)
            .map(|_| (0..num).map(|_| F::rand()).collect::<Vec<_>>())
            .collect::<Vec<_>>();

        let now = Instant::now();
        let pis_commitments =
            pis_commitments(&wires, &precomputed_lagrange_commitments, num, false);
        dbg!(now.elapsed());

        let now = Instant::now();
        let manual_pis_commitments = wires
            .into_iter()
            .map(|w| {
                manually_compute_pis_commitments(
                    &w,
                    &msm_precomputation,
                    &fft_precomputation,
                    order,
                )
            })
            .collect::<Vec<_>>();
        dbg!(now.elapsed());

        pis_commitments
            .into_iter()
            .zip(manual_pis_commitments)
            .for_each(|(a, b)| {
                assert_eq!(a, b);
            });
        dbg!(noww.elapsed());
    }
}
