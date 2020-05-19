use std::time::Instant;

use rayon::prelude::*;

use crate::{affine_multisummation_best, AffinePoint, Curve, Field, ProjectivePoint};

/// In Yao's method, we compute an affine summation for each digit. In a parallel setting, it would
/// be easiest to assign individual summations to threads, but this would be sub-optimal because
/// multi-summations can be more efficient than repeating individual summations (see
/// `affine_multisummation_best`). Thus we divide digits into large chunks, and assign chunks of
/// digits to threads. Note that there is a delicate balance here, as large chunks can result in
/// uneven distributions of work among threads.
const DIGITS_PER_CHUNK: usize = 80;

pub struct MsmPrecomputation<C: Curve> {
    /// For each generator (in the order they were passed to `msm_precompute`), contains a vector
    /// of powers, i.e. [(2^w)^i] for i < DIGITS.
    // TODO: Use compressed coordinates here.
    powers_per_generator: Vec<Vec<AffinePoint<C>>>,

    /// The window size.
    w: usize,
}

pub fn msm_precompute<C: Curve>(generators: &[ProjectivePoint<C>], w: usize) -> MsmPrecomputation<C> {
    MsmPrecomputation {
        powers_per_generator: generators.into_par_iter()
            .map(|&g| precompute_single_generator(g, w))
            .collect(),
        w,
    }
}

fn precompute_single_generator<C: Curve>(g: ProjectivePoint<C>, w: usize) -> Vec<AffinePoint<C>> {
    let digits = (C::ScalarField::BITS + w - 1) / w;
    let mut powers: Vec<ProjectivePoint<C>> = Vec::with_capacity(digits);
    powers.push(g);
    for i in 1..digits {
        let mut power_i_proj = powers[i - 1];
        for _j in 0..w {
            power_i_proj = power_i_proj.double();
        }
        powers.push(power_i_proj);
    }
    ProjectivePoint::batch_to_affine(&powers)
}

pub fn msm_execute<C: Curve>(
    precomputation: &MsmPrecomputation<C>,
    scalars: &[C::ScalarField],
) -> ProjectivePoint<C> {
    assert_eq!(precomputation.powers_per_generator.len(), scalars.len());
    let w = precomputation.w;
    let digits = (C::ScalarField::BITS + w - 1) / w;
    let base = 1 << w;

    // This is a variant of Yao's method, adapted to the multi-scalar setting. Because we use
    // extremely large windows, the repeated scans in Yao's method could be more expensive than the
    // actual group operations. To avoid this, we store a multimap from each possible digit to the
    // positions in which that digit occurs in the scalars. These positions have the form (i, j),
    // where i is the index of the generator and j is an index into the digits of the scalar
    // associated with that generator.
    let mut digit_occurrences: Vec<Vec<(usize, usize)>> = Vec::with_capacity(digits);
    for _i in 0..base {
        digit_occurrences.push(Vec::new());
    }
    for (i, scalar) in scalars.iter().enumerate() {
        let digits = to_digits::<C>(scalar, w);
        for (j, &digit) in digits.iter().enumerate() {
            digit_occurrences[digit].push((i, j));
        }
    }

    let mut y = ProjectivePoint::ZERO;
    let mut u = ProjectivePoint::ZERO;

    for digit in (1..base).rev() {
        for &(i, j) in &digit_occurrences[digit] {
            u = u + precomputation.powers_per_generator[i][j];
        }
        y = y + u;
    }

    y
}

pub fn msm_execute_parallel<C: Curve>(
    precomputation: &MsmPrecomputation<C>,
    scalars: &[C::ScalarField],
    w: usize,
) -> ProjectivePoint<C> {
    assert_eq!(precomputation.powers_per_generator.len(), scalars.len());
    let digits = (C::ScalarField::BITS + w - 1) / w;
    let base = 1 << w;

    // This is a variant of Yao's method, adapted to the multi-scalar setting. Because we use
    // extremely large windows, the repeated scans in Yao's method could be more expensive than the
    // actual group operations. To avoid this, we store a multimap from each possible digit to the
    // positions in which that digit occurs in the scalars. These positions have the form (i, j),
    // where i is the index of the generator and j is an index into the digits of the scalar
    // associated with that generator.
    let mut digit_occurrences: Vec<Vec<(usize, usize)>> = Vec::with_capacity(digits);
    for _i in 0..base {
        digit_occurrences.push(Vec::new());
    }
    for (i, scalar) in scalars.iter().enumerate() {
        let digits = to_digits::<C>(scalar, w);
        for (j, &digit) in digits.iter().enumerate() {
            digit_occurrences[digit].push((i, j));
        }
    }

    // For each digit, we add up the powers associated with all occurrences that digit.
    let digits: Vec<usize> = (0..base).collect();
    let start = Instant::now();
    let digit_acc: Vec<ProjectivePoint<C>> = digits.par_chunks(DIGITS_PER_CHUNK)
        .flat_map(|chunk| {
            let summations: Vec<Vec<AffinePoint<C>>> = chunk.iter()
                .map(|&digit| digit_occurrences[digit].iter()
                    .map(|&(i, j)| precomputation.powers_per_generator[i][j])
                    .collect())
                .collect();
            affine_multisummation_best(summations)
        })
        .collect();
    println!("Computing the per-digit summations (in parallel) took {}s", start.elapsed().as_secs_f64());

    let start = Instant::now();
    let mut y = ProjectivePoint::ZERO;
    let mut u = ProjectivePoint::ZERO;
    for digit in (1..base).rev() {
        u = u + digit_acc[digit];
        y = y + u;
    }
    println!("Final summation (sequential) {}s", start.elapsed().as_secs_f64());
    y
}

pub(crate) fn to_digits<C: Curve>(x: &C::ScalarField, w: usize) -> Vec<usize> {
    let scalar_bits = C::ScalarField::BITS;
    let num_digits = (scalar_bits + w - 1) / w;

    // Convert x to a bool array.
    let x_canonical = x.to_canonical_u64_vec();
    let mut x_bits = Vec::with_capacity(scalar_bits);
    for i in 0..scalar_bits {
        x_bits.push((x_canonical[i / 64] >> (i as u64 % 64) & 1) != 0);
    }

    let mut digits = Vec::with_capacity(num_digits);
    for i in 0..num_digits {
        let mut digit = 0;
        for j in ((i * w)..((i + 1) * w).min(scalar_bits)).rev() {
            digit <<= 1;
            digit |= x_bits[j] as usize;
        }
        digits.push(digit);
    }
    digits
}

#[cfg(test)]
mod tests {
    use crate::{Bls12377, Bls12377Scalar, Curve, msm_execute, msm_precompute, to_digits};

    #[test]
    fn test_to_digits() {
        let x_canonical = [
            0b1010101010101010101010101010101010101010101010101010101010101010,
            0b1100110011001100110011001100110011001100110011001100110011001100,
            0b1111000011110000111100001111000011110000111100001111000011110000,
            0b0000111111111111111111111111111111111111111111111111111111111111];
        let x = Bls12377Scalar::from_canonical(x_canonical);
        assert_eq!(x.to_canonical(), x_canonical);
        assert_eq!(to_digits::<Bls12377>(&x, 17), vec![
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

        let generator_1 = Bls12377::GENERATOR_PROJECTIVE;
        let generator_2 = generator_1 + generator_1;
        let generator_3 = generator_1 + generator_2;

        let scalar_1 = Bls12377Scalar::from_canonical([11111111, 22222222, 33333333, 44444444]);
        let scalar_2 = Bls12377Scalar::from_canonical([22222222, 22222222, 33333333, 44444444]);
        let scalar_3 = Bls12377Scalar::from_canonical([33333333, 22222222, 33333333, 44444444]);

        let generators = vec![generator_1, generator_2, generator_3];
        let scalars = vec![scalar_1, scalar_2, scalar_3];

        let precomputation = msm_precompute(&generators, w);
        let result_msm = msm_execute(&precomputation, &scalars);

        let result_naive = scalar_1 * generator_1 + scalar_2 * generator_2 + scalar_3 * generator_3;
        assert_eq!(result_msm, result_naive);
    }
}
