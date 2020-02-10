use std::ops::Mul;

use chashmap::CHashMap;

use lazy_static::lazy_static;

use crate::{Bls12Scalar, G1ProjectivePoint};

const WINDOW_BITS: usize = 4;
const BASE: usize = 1 << WINDOW_BITS;
const DIGITS: usize = (Bls12Scalar::BITS + WINDOW_BITS - 1) / WINDOW_BITS;

/// Precomputed state used for Bls12Scalar x G1ProjectivePoint multiplications,
/// specific to a particular generator.
#[derive(Copy, Clone)]
struct G1GeneratorPrecomputations {
    /// [(2^w)^i] g for each i < DIGITS
    powers: [G1ProjectivePoint; DIGITS],
}

// TODO: Use compressed coordinates in the cache.
lazy_static! {
    static ref G1_MUL_PRECOMPUTATIONS: CHashMap<G1ProjectivePoint, G1GeneratorPrecomputations> = CHashMap::new();
}

fn get_precomputation(g: G1ProjectivePoint) -> G1GeneratorPrecomputations {
    match G1_MUL_PRECOMPUTATIONS.get(&g) {
        Some(x) => *x,
        None => {
            let precomputation = precompute(g);
            G1_MUL_PRECOMPUTATIONS.insert(g, precomputation);
            precomputation
        }
    }
}

fn precompute(g: G1ProjectivePoint) -> G1GeneratorPrecomputations {
    let mut powers = [G1ProjectivePoint::ZERO; DIGITS];
    powers[0] = g;
    for i in 1..DIGITS {
        powers[i] = powers[i - 1];
        for _j in 0..WINDOW_BITS {
            powers[i] = powers[i] + powers[i];
        }
    }
    G1GeneratorPrecomputations { powers }
}

impl Mul<G1ProjectivePoint> for Bls12Scalar {
    type Output = G1ProjectivePoint;

    fn mul(self, rhs: G1ProjectivePoint) -> Self::Output {
        // Yao's method; see https://koclab.cs.ucsb.edu/teaching/ecc/eccPapers/Doche-ch09.pdf
        let precomputed_powers = get_precomputation(rhs).powers;

        let digits = to_digits(&self);

        let mut y = G1ProjectivePoint::ZERO;
        let mut u = G1ProjectivePoint::ZERO;
        let mut ops = 0; // todo
        for j in (1..BASE).rev() {
            for (i, &digit) in digits.iter().enumerate() {
                if digit == j as u64 {
                    u = u + precomputed_powers[i];
                    ops += 1;
                }
            }
            y = y + u;
            ops += 1;
        }
        y
    }
}

fn to_digits(x: &Bls12Scalar) -> [u64; DIGITS] {
    debug_assert!(64 % WINDOW_BITS == 0,
                  "For simplicity, only power-of-two window sizes are handled for now");
    let digits_per_u64 = 64 / WINDOW_BITS;
    let mut digits = [0; DIGITS];
    for (i, limb) in x.limbs.iter().enumerate() {
        for j in 0..digits_per_u64 {
            digits[i * digits_per_u64 + j] = (limb >> (j * WINDOW_BITS)) % BASE as u64;
        }
    }

    digits
}

#[cfg(test)]
mod tests {
    use crate::{G1_GENERATOR, Bls12Scalar};

    #[test]
    fn mul() {
        let ten = Bls12Scalar::from(10);
        let product = ten * G1_GENERATOR;
        let sum = G1_GENERATOR + G1_GENERATOR + G1_GENERATOR + G1_GENERATOR + G1_GENERATOR
            + G1_GENERATOR + G1_GENERATOR + G1_GENERATOR + G1_GENERATOR + G1_GENERATOR;
        assert_eq!(product, sum);
    }
}
