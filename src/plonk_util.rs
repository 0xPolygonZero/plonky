use crate::{
    fft_with_precomputation_power_of_2, ifft_with_precomputation_power_of_2, pedersen_hash,
    AffinePoint, CircuitBuilder, Curve, FftPrecomputation, Field, HaloCurve, MsmPrecomputation,
    ProjectivePoint, Target,
};

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
    let one = C::ScalarField::ONE;
    let two = C::ScalarField::TWO;

    let mut a = zero;
    let mut b = zero;

    for s_bits_chunk in s_bits.chunks(2) {
        let bit_lo = s_bits_chunk[0];
        let bit_hi = s_bits_chunk[1];

        let sign = two * C::ScalarField::from_canonical_bool(bit_lo) - one;
        let (c, d) = if bit_hi { (zero, sign) } else { (sign, zero) };

        a = a.double() + c;
        b = b.double() + d;
    }

    a * C::ZETA_SCALAR + b
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

pub(crate) fn coeffs_to_commitments<C: Curve>(
    coefficients_vec: &[Vec<C::ScalarField>],
    msm_precomputation: &MsmPrecomputation<C>,
) -> Vec<AffinePoint<C>> {
    let projs: Vec<_> = coefficients_vec
        .iter()
        .map(|coeffs| pedersen_hash(coeffs, msm_precomputation))
        .collect();
    ProjectivePoint::batch_to_affine(&projs)
}
