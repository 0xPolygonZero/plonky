use std::ops::Mul;

use crate::{affine_summation_batch_inversion, AffinePoint, Curve, CurveScalar, Field, ProjectivePoint};

const WINDOW_BITS: usize = 4;
const BASE: usize = 1 << WINDOW_BITS;

fn digits_per_scalar<C: Curve>() -> usize {
    (C::ScalarField::BITS + WINDOW_BITS - 1) / WINDOW_BITS
}

/// Precomputed state used for scalar x ProjectivePoint multiplications,
/// specific to a particular generator.
#[derive(Clone)]
pub struct MultiplicationPrecomputation<C: Curve> {
    /// [(2^w)^i] g for each i < digits_per_scalar.
    powers: Vec<AffinePoint<C>>,
}

impl<C: Curve> ProjectivePoint<C> {
    pub fn mul_precompute(&self) -> MultiplicationPrecomputation<C> {
        let num_digits = digits_per_scalar::<C>();
        let mut powers_proj = Vec::with_capacity(num_digits);
        powers_proj.push(*self);
        for i in 1..num_digits {
            let mut power_i_proj = powers_proj[i - 1];
            for _j in 0..WINDOW_BITS {
                power_i_proj = power_i_proj.double();
            }
            powers_proj.push(power_i_proj);
        }

        let powers = ProjectivePoint::batch_to_affine(&powers_proj);
        MultiplicationPrecomputation { powers }
    }

    pub fn mul_with_precomputation(
        &self,
        scalar: C::ScalarField,
        precomputation: MultiplicationPrecomputation<C>,
    ) -> Self {
        // Yao's method; see https://koclab.cs.ucsb.edu/teaching/ecc/eccPapers/Doche-ch09.pdf
        let precomputed_powers = precomputation.powers;

        let digits = to_digits::<C>(&scalar);

        let mut y = ProjectivePoint::ZERO;
        let mut u = ProjectivePoint::ZERO;
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

impl<C: Curve> Mul<ProjectivePoint<C>> for CurveScalar<C> {
    type Output = ProjectivePoint<C>;

    fn mul(self, rhs: ProjectivePoint<C>) -> Self::Output {
        let precomputation = rhs.mul_precompute();
        rhs.mul_with_precomputation(self.0, precomputation)
    }
}

#[allow(clippy::assertions_on_constants)]
fn to_digits<C: Curve>(x: &C::ScalarField) -> Vec<u64> {
    debug_assert!(
        64 % WINDOW_BITS == 0,
        "For simplicity, only power-of-two window sizes are handled for now"
    );
    let digits_per_u64 = 64 / WINDOW_BITS;
    let mut digits = Vec::with_capacity(digits_per_scalar::<C>());
    for limb in x.to_canonical_u64_vec() {
        for j in 0..digits_per_u64 {
            digits.push((limb >> (j * WINDOW_BITS) as u64) % BASE as u64);
        }
    }

    digits
}
