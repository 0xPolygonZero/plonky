use std::ops::Mul;

use crate::{Bls12Scalar, G1ProjectivePoint, G1AffinePoint, affine_summation_batch_inversion};

const WINDOW_BITS: usize = 4;
const BASE: usize = 1 << WINDOW_BITS;
const DIGITS: usize = (Bls12Scalar::BITS + WINDOW_BITS - 1) / WINDOW_BITS;

/// Precomputed state used for Bls12Scalar x G1ProjectivePoint multiplications,
/// specific to a particular generator.
#[derive(Clone)]
pub struct G1GeneratorPrecomputations {
    /// [(2^w)^i] g for each i < DIGITS.
    powers: Vec<G1AffinePoint>,
}

impl G1ProjectivePoint {
    pub fn mul_precompute(&self) -> G1GeneratorPrecomputations {
        let mut powers_proj = Vec::with_capacity(DIGITS);
        powers_proj.push(*self);
        for i in 1..DIGITS {
            let mut power_i_proj = powers_proj[i - 1];
            for _j in 0..WINDOW_BITS {
                power_i_proj = power_i_proj.double();
            }
            powers_proj.push(power_i_proj);
        }

        let powers = G1ProjectivePoint::batch_to_affine(&powers_proj);
        G1GeneratorPrecomputations { powers }
    }

    pub fn mul_with_precomputation(
        &self,
        scalar: Bls12Scalar,
        precomputation: G1GeneratorPrecomputations,
    ) -> Self {
        // Yao's method; see https://koclab.cs.ucsb.edu/teaching/ecc/eccPapers/Doche-ch09.pdf
        let precomputed_powers = precomputation.powers;

        let digits = to_digits(&scalar);

        let mut y = G1ProjectivePoint::ZERO;
        let mut u = G1ProjectivePoint::ZERO;
        for j in (1..BASE).rev() {
            let mut u_summands = Vec::new();
            for (i, &digit) in digits.iter().enumerate() {
                if digit == j as u64 {
                    u_summands.push(precomputed_powers[i]);
                }
            }
            u = u + affine_summation_batch_inversion(u_summands);
            y = y + u;
        }
        y
    }
}

impl Mul<G1ProjectivePoint> for Bls12Scalar {
    type Output = G1ProjectivePoint;

    fn mul(self, rhs: G1ProjectivePoint) -> Self::Output {
        let precomputation = rhs.mul_precompute();
        rhs.mul_with_precomputation(self, precomputation)
    }
}

fn to_digits(x: &Bls12Scalar) -> [u64; DIGITS] {
    debug_assert!(64 % WINDOW_BITS == 0,
                  "For simplicity, only power-of-two window sizes are handled for now");
    let digits_per_u64 = 64 / WINDOW_BITS;
    let mut digits = [0; DIGITS];
    for (i, limb) in x.to_canonical().iter().enumerate() {
        for j in 0..digits_per_u64 {
            digits[i * digits_per_u64 + j] = (limb >> (j * WINDOW_BITS)) % BASE as u64;
        }
    }

    digits
}
