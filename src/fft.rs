use rayon::prelude::*;

use crate::Bls12Scalar;

/// Permutes `arr` such that each index is mapped to its reverse in binary.
fn reverse_index_bits<T: Copy>(arr: Vec<T>) -> Vec<T> {
    let n = arr.len();
    let n_power = log2_strict(n);

    let mut result = Vec::with_capacity(n);
    for i in 0..n {
        result.push(arr[reverse_bits(i, n_power)]);
    }
    result
}

fn reverse_bits(n: usize, num_bits: usize) -> usize {
    let mut result = 0;
    for i in 0..num_bits {
        let i_rev = num_bits - i - 1;
        result |= (n >> i & 1) << i_rev;
    }
    result
}

/// Computes `ceil(log_2(n))`.
fn log2_ceil(n: usize) -> usize {
    let mut exp = 0;
    while 1 << exp < n {
        exp += 1;
    }
    exp
}

/// Computes `log_2(n)`, panicking if `n` is not a power of two.
fn log2_strict(n: usize) -> usize {
    let mut exp = log2_ceil(n);
    assert_eq!(1 << exp, n, "Input not a power of 2");
    exp
}

pub struct FftPrecomputation {
    /// For each layer index i, stores the cyclic subgroup corresponding to the evaluation domain of
    /// layer i. The indices within these subgroup vectors are bit-reversed.
    subgroups_rev: Vec<Vec<Bls12Scalar>>,
}

pub fn fft(coefficients: &[Bls12Scalar]) -> Vec<Bls12Scalar> {
    let precomputation = fft_precompute(coefficients.len());
    fft_with_precomputation(coefficients, &precomputation)
}

pub fn fft_precompute(degree: usize) -> FftPrecomputation {
    let degree_pow = log2_ceil(degree);

    let mut subgroups_rev = Vec::new();
    for i in 0..=degree_pow {
        let g_i = Bls12Scalar::primitive_root_of_unity(i);
        let subgroup = Bls12Scalar::cyclic_subgroup_known_order(g_i, 1 << i);
        let subgroup_rev = reverse_index_bits(subgroup);
        subgroups_rev.push(subgroup_rev);
    }

    FftPrecomputation { subgroups_rev }
}

pub fn fft_with_precomputation(
    coefficients: &[Bls12Scalar],
    precomputation: &FftPrecomputation,
) -> Vec<Bls12Scalar> {
    let degree = coefficients.len();
    let degree_padded = 1 << log2_ceil(degree);

    if degree == degree_padded {
        fft_with_precomputation_power_of_2(coefficients, precomputation)
    } else {
        let mut coefficients_padded = Vec::with_capacity(degree_padded);
        for &c in coefficients {
            coefficients_padded.push(c);
        }
        for _i in degree..degree_padded {
            coefficients_padded.push(Bls12Scalar::ZERO);
        }
        fft_with_precomputation_power_of_2(&coefficients_padded, precomputation)
    }
}

pub fn ifft_with_precomputation_power_of_2(
    points: &[Bls12Scalar],
    precomputation: &FftPrecomputation,
) -> Vec<Bls12Scalar> {
    let n = points.len();
    let n_inv = Bls12Scalar::from_canonical_usize(n).multiplicative_inverse().unwrap();
    let mut result = fft_with_precomputation(points, precomputation);
    // TODO: Could do this in-place with swaps.
    let mut coefficients = Vec::with_capacity(n);
    coefficients.push(result[0] * n_inv);
    for x in result.into_iter().skip(1).rev() {
        coefficients.push(x * n_inv);
    }
    coefficients
}

pub fn fft_with_precomputation_power_of_2(
    coefficients: &[Bls12Scalar],
    precomputation: &FftPrecomputation,
) -> Vec<Bls12Scalar> {
    let degree = coefficients.len();
    let half_degree = coefficients.len() >> 1;
    let degree_pow = log2_strict(degree);

    // In the base layer, we're just evaluating "degree 0 polynomials", i.e. the coefficients
    // themselves.
    let coefficients_vec = coefficients.to_vec();
    let mut evaluations = reverse_index_bits(coefficients_vec);

    for i in 1..=degree_pow {
        // In layer i, we're evaluating a series of polynomials, each at 2^i points. In practice
        // we evaluate a pair of points together, so we have 2^(i - 1) pairs.
        let points_per_poly = 1 << i;
        let pairs_per_poly = 1 << (i - 1);

        let pair_indices: Vec<usize> = (0..half_degree).collect();
        let new_evaluations = pair_indices.par_chunks(2000).flat_map(|chunk| {
            let mut new_evaluations_chunk = Vec::new();
            for pair_index in chunk {
                let poly_index = pair_index / pairs_per_poly;
                let pair_index_within_poly = pair_index % pairs_per_poly;

                let child_index_0 = poly_index * points_per_poly + pair_index_within_poly;
                let child_index_1 = child_index_0 + pairs_per_poly;

                let even = evaluations[child_index_0];
                let odd = evaluations[child_index_1];

                let point_0 = precomputation.subgroups_rev[i][pair_index_within_poly * 2];
                let product = point_0 * odd;
                new_evaluations_chunk.push(even + product);
                new_evaluations_chunk.push(even - product);
            }
            new_evaluations_chunk
        }).collect();
        evaluations = new_evaluations;
    }

    // Reorder so that evaluations' indices correspond to (g_0, g_1, g_2, ...)
    reverse_index_bits(evaluations)
}

#[cfg(test)]
mod tests {
    use crate::{Bls12Scalar, ifft_with_precomputation_power_of_2, fft_precompute, fft_with_precomputation};
    use crate::fft::{log2_ceil, log2_strict, reverse_bits, reverse_index_bits};

    #[test]
    fn fft_and_ifft() {
        let degree = 200;
        let degree_padded = log2_ceil(degree);
        let mut coefficients = Vec::new();
        for i in 0..degree {
            coefficients.push(Bls12Scalar::from_canonical_usize(i * 1337 % 100));
        }

        let precomputation = fft_precompute(degree);
        let points = fft_with_precomputation(&coefficients, &precomputation);
        assert_eq!(points, evaluate_naive(&coefficients));

        let interpolated_coefficients = ifft_with_precomputation_power_of_2(&points, &precomputation);
        for i in 0..degree {
            assert_eq!(interpolated_coefficients[i], coefficients[i]);
        }
        for i in degree..degree_padded {
            assert_eq!(interpolated_coefficients[i], Bls12Scalar::ZERO);
        }
    }

    #[test]
    fn test_reverse_bits() {
        assert_eq!(reverse_bits(0b00110101, 8), 0b10101100);
        assert_eq!(reverse_index_bits(vec!["a", "b"]), vec!["a", "b"]);
        assert_eq!(reverse_index_bits(vec!["a", "b", "c", "d"]), vec!["a", "c", "b", "d"]);
    }

    fn evaluate_naive(coefficients: &[Bls12Scalar]) -> Vec<Bls12Scalar> {
        let degree = coefficients.len();
        let degree_padded = 1 << log2_ceil(degree);

        let mut coefficients_padded = Vec::with_capacity(degree_padded);
        for c in coefficients {
            coefficients_padded.push(*c);
        }
        for _i in degree..degree_padded {
            coefficients_padded.push(Bls12Scalar::ZERO);
        }
        evaluate_naive_power_of_2(&coefficients_padded)
    }

    fn evaluate_naive_power_of_2(coefficients: &[Bls12Scalar]) -> Vec<Bls12Scalar> {
        let degree = coefficients.len();
        let degree_pow = log2_strict(degree);

        let g = Bls12Scalar::primitive_root_of_unity(degree_pow);
        let powers_of_g = Bls12Scalar::cyclic_subgroup_known_order(g, degree);

        powers_of_g.into_iter().map(|x| evaluate_at_naive(&coefficients, x)).collect()
    }

    fn evaluate_at_naive(coefficients: &[Bls12Scalar], point: Bls12Scalar) -> Bls12Scalar {
        let mut sum = Bls12Scalar::ZERO;
        let mut point_power = Bls12Scalar::ONE;
        for &c in coefficients {
            sum = sum + c * point_power;
            point_power = point_power * point;
        }
        sum
    }
}
