// Since we work with elliptic curve groups etc., nearly all the operations are suspicious to
// Clippy.
#![allow(clippy::suspicious_arithmetic_impl)]

// We have tons of bigint literals in Montgomery form, which won't be readable with or without underscores.
#![allow(clippy::unreadable_literal)]

pub use bigint_arithmetic::*;
pub use bls12_377_base::*;
pub use bls12_377_curve::*;
pub use bls12_377_scalar::*;
pub use curve::*;
pub use curve_adds::*;
pub use curve_msm::*;
pub use curve_multiplication::*;
pub use curve_summations::*;
pub use fft::*;
pub use field::*;
pub use plonk::*;
pub use poly_commit::*;
pub use tweedledee_base::*;
pub use tweedledee_curve::*;
pub use tweedledum_base::*;
pub use tweedledum_curve::*;

mod bigint_arithmetic;
mod bigint_inverse;
mod curve_adds;
mod bls12_377_base;
mod bls12_377_curve;
mod bls12_377_scalar;
mod curve_summations;
mod fft;
mod field;
mod curve;
mod curve_msm;
mod curve_multiplication;
mod plonk;
mod plonk_gates;
mod plonk_recursion;
mod poly_commit;
mod tweedledee_base;
mod tweedledee_curve;
mod tweedledum_base;
mod tweedledum_curve;

#[cfg(test)]
mod conversions;
