use std::collections::HashMap;

use crate::{Bls12Scalar, G1AffinePoint, G1ProjectivePoint};

pub struct MsmPrecomputation {
    /// For each generator (in the order they were passed to `msm_precompute`), contains a vector
    /// of powers, i.e. [(2^w)^i] for i < DIGITS.
    // TODO: Use compressed coordinates here.
    powers_per_generator: Vec<Vec<G1AffinePoint>>,
}

pub fn msm_precompute(generators: &Vec<G1ProjectivePoint>, w: usize) -> MsmPrecomputation {
    MsmPrecomputation {
        powers_per_generator: generators.iter()
            .map(|g| precompute_single_generator(g, w))
            .collect()
    }
}

fn precompute_single_generator(g: &G1ProjectivePoint, w: usize) -> Vec<G1AffinePoint> {
    let digits = (Bls12Scalar::BITS + w - 1) / w;
    let mut powers = Vec::with_capacity(digits);
    powers.push(g.to_affine());
    for i in 1..digits {
        let mut power_i_proj = powers[i - 1].to_projective();
        for _j in 0..w {
            power_i_proj = power_i_proj.double();
        }
        powers.push(power_i_proj.to_affine());
    }
    powers
}

pub fn msm_execute(
    precomputation: &MsmPrecomputation,
    scalars: &Vec<Bls12Scalar>,
    w: usize,
) -> G1ProjectivePoint {
    assert_eq!(precomputation.powers_per_generator.len(), scalars.len());
    let digits = (Bls12Scalar::BITS + w - 1) / w;
    let base = 1 << w;

    // This is a variant of Yao's method, adapted to the multi-scalar setting. Because we use
    // extremely large windows, the repeated scans in Yao's method could be more expensive than the
    // actual group operations. To avoid this, we store a multimap from each possible digit to the
    // positions in which that digit occurs in the scalars. These positions have the form (i, j),
    // where i is the index of the generator and j is an index into the digits of the scalar
    // associated with that generator.
    let mut digit_occurrences: Vec<Vec<(usize, usize)>> = Vec::with_capacity(digits);
    for i in 0..digits {
        digit_occurrences.push(Vec::new());
    }
    for (i, scalar) in scalars.iter().enumerate() {
        let digits = to_digits(scalar, w);
        for (j, &digit) in digits.iter().enumerate() {
            digit_occurrences[digit].push((i, j));
        }
    }

    let mut y = G1ProjectivePoint::ZERO;
    let mut u = G1ProjectivePoint::ZERO;

    for digit in (1..base).rev() {
        for &(i, j) in &digit_occurrences[digit] {
            u = u + precomputation.powers_per_generator[i][j];
        }
        y = y + u;
    }

    y
}

pub(crate) fn to_digits(x: &Bls12Scalar, w: usize) -> Vec<usize> {
    let num_digits = (Bls12Scalar::BITS + w - 1) / w;

    // Convert x to a bool array.
    let x_canonical = x.to_canonical();
    let mut x_bits = [false; Bls12Scalar::BITS];
    for i in 0..Bls12Scalar::BITS {
        x_bits[i] = (x_canonical[i / 64] >> (i as u64 % 64) & 1) != 0;
    }

    let mut digits = Vec::with_capacity(num_digits);
    for i in 0..num_digits {
        let mut digit = 0;
        for j in ((i * w)..((i + 1) * w).min(Bls12Scalar::BITS)).rev() {
            digit <<= 1;
            digit |= x_bits[j] as usize;
        }
        digits.push(digit);
    }
    digits
}

#[cfg(test)]
mod tests {
    use crate::{Bls12Scalar, to_digits, G1_GENERATOR, msm_precompute, msm_execute};

    #[test]
    fn test_to_digits() {
        let x_canonical = [
            0b1010101010101010101010101010101010101010101010101010101010101010,
            0b1100110011001100110011001100110011001100110011001100110011001100,
            0b1111000011110000111100001111000011110000111100001111000011110000,
            0b0000111111111111111111111111111111111111111111111111111111111111];
        let x = Bls12Scalar::from_canonical(x_canonical);
        assert_eq!(x.to_canonical(), x_canonical);
        assert_eq!(to_digits(&x, 17), vec![
            0b01010101010101010,
            0b10101010101010101,
            0b01010101010101010,
            0b11001010101010101,
            0b01100110011001100,
            0b00110011001100110,
            0b10011001100110011,
            0b11110000110011001,
            0b01111000011110000,
            0b00111100001111000,
            0b00011110000111100,
            0b11111111111111110,
            0b11111111111111111,
            0b11111111111111111,
            0b00011111111111111,
        ]);
    }

    #[test]
    fn test_msm() {
        let w = 5;

        let generator_1 = G1_GENERATOR;
        let generator_2 = generator_1 + generator_1;
        let generator_3 = generator_1 + generator_2;

        let scalar_1 = Bls12Scalar::from_canonical([11111111, 22222222, 33333333, 44444444]);
        let scalar_2 = Bls12Scalar::from_canonical([22222222, 22222222, 33333333, 44444444]);
        let scalar_3 = Bls12Scalar::from_canonical([33333333, 22222222, 33333333, 44444444]);

        let generators = vec![generator_1, generator_2, generator_3];
        let scalars = vec![scalar_1, scalar_2, scalar_3];

        let precomputation = msm_precompute(&generators, w);
        let result_msm = msm_execute(&precomputation, &scalars, w);

        let result_naive = scalar_1 * generator_1 + scalar_2 * generator_2 + scalar_3 * generator_3;
        assert_eq!(result_msm, result_naive);
    }
}
