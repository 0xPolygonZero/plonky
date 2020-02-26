use chashmap::CHashMap;

use lazy_static::lazy_static;

use crate::{Bls12Base, Bls12Scalar};

lazy_static! {
    static ref SUBGROUPS_BY_ORDER_POWER: CHashMap<usize, Vec<Bls12Scalar>> = CHashMap::new();
}

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

fn get_subgroup(order_power: usize) -> Vec<Bls12Scalar> {
    match SUBGROUPS_BY_ORDER_POWER.get(&order_power) {
        Some(subgroup) => subgroup.clone(),
        None => {
            let subgroup = Bls12Scalar::cyclic_subgroup_unknown_order(
                Bls12Scalar::primitive_root_of_unity(order_power));
            SUBGROUPS_BY_ORDER_POWER.insert(order_power, subgroup.clone());
            subgroup
        }
    }
}

pub fn fft(coefficients: &[Bls12Scalar]) -> Vec<Bls12Scalar> {
    let degree = coefficients.len();
    let degree_padded = 1 << log2_ceil(degree);

    if degree == degree_padded {
        fft_power_of_2(coefficients)
    } else {
        let mut coefficients_padded = Vec::with_capacity(degree_padded);
        for &c in coefficients {
            coefficients_padded.push(c);
        }
        for _i in degree..degree_padded {
            coefficients_padded.push(Bls12Scalar::ZERO);
        }
        fft_power_of_2(&coefficients_padded)
    }
}

pub fn fft_power_of_2(coefficients: &[Bls12Scalar]) -> Vec<Bls12Scalar> {
    let degree = coefficients.len();
    let degree_pow = log2_strict(degree);

    // In the base layer, we're just evaluating "degree 0 polynomials", i.e. the coefficients
    // themselves.
    let coefficients_vec = coefficients.to_vec();
    let mut evaluations = reverse_index_bits(coefficients_vec);

    for i in 1..=degree_pow {
        let chunk_size = 1 << i;
        let half_chunk_size = 1 << (i - 1);
        let num_chunks = 1 << (degree_pow - i);

        // TODO: This stuff can be precomputed.
        let g_i = Bls12Scalar::primitive_root_of_unity(i);
        let powers_of_g_i = Bls12Scalar::cyclic_subgroup_known_order(g_i, 1 << i);
        let powers_of_g_i_rev = reverse_index_bits(powers_of_g_i);

        let mut new_evaluations = Vec::new();
        for j in 0..num_chunks {
            let first_chunk_index = j * chunk_size;
            for k in 0..half_chunk_size {
                let index_0 = first_chunk_index + k * 2;
                let index_1 = index_0 + 1;
                let child_index_0 = first_chunk_index + k;
                let child_index_1 = child_index_0 + half_chunk_size;

                let even = evaluations[child_index_0];
                let odd = evaluations[child_index_1];

                let point_0 = powers_of_g_i_rev[k * 2];
                let product = point_0 * odd;
                new_evaluations.push(even + product);
                new_evaluations.push(even - product);
            }
        }
        evaluations = new_evaluations;
    }

    // Reorder so that evaluations' indices correspond to (g_0, g_1, g_2, ...)
    reverse_index_bits(evaluations)
}

#[cfg(test)]
mod tests {
    use rand::RngCore;
    use rand::rngs::OsRng;

    use crate::{Bls12Scalar, fft};
    use crate::fft::{log2_strict, reverse_bits, reverse_index_bits, log2_ceil};

    #[test]
    fn test_fft() {
        let degree = 200;
        let mut coefficients = Vec::new();
        for i in 0..degree {
            coefficients.push(Bls12Scalar::from_canonical_u64(i * 1337 % 100));
        }
        assert_eq!(fft(&coefficients), evaluate_naive(coefficients));
    }

    #[test]
    fn test_reverse_bits() {
        assert_eq!(reverse_bits(0b00110101, 8), 0b10101100);
        assert_eq!(reverse_index_bits(vec!["a", "b"]), vec!["a", "b"]);
        assert_eq!(reverse_index_bits(vec!["a", "b", "c", "d"]), vec!["a", "c", "b", "d"]);
    }

    fn evaluate_naive(coefficients: Vec<Bls12Scalar>) -> Vec<Bls12Scalar> {
        let degree = coefficients.len();
        let degree_padded = 1 << log2_ceil(degree);

        let mut coefficients_padded = Vec::with_capacity(degree_padded);
        for c in coefficients {
            coefficients_padded.push(c);
        }
        for _i in degree..degree_padded {
            coefficients_padded.push(Bls12Scalar::ZERO);
        }
        evaluate_naive_power_of_2(coefficients_padded)
    }

    fn evaluate_naive_power_of_2(coefficients: Vec<Bls12Scalar>) -> Vec<Bls12Scalar> {
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
