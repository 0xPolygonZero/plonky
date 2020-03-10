// Since we work with elliptic curve groups etc., nearly all the operations are suspicious to
// Clippy.
#![allow(clippy::suspicious_arithmetic_impl)]

// We have tons of bigint literals in Montgomery form, which won't be readable with or without underscores.
#![allow(clippy::unreadable_literal)]

pub use bigint_arithmetic::*;
pub use bls12_377_adds::*;
pub use bls12_377_base::*;
pub use bls12_377_curve::*;
pub use bls12_377_scalar::*;
pub use bls12_377_summations::*;
pub use fft::*;
pub use field::*;
pub use group::*;
pub use group_msm::*;
pub use group_multiplication::*;
pub use plonk::*;
pub use poly_commit::*;

mod bigint_arithmetic;
mod bls12_377_adds;
mod bls12_377_base;
mod bls12_377_curve;
mod bls12_377_scalar;
mod bls12_377_summations;
mod fft;
mod field;
mod group;
mod group_msm;
mod group_multiplication;
mod plonk;
mod poly_commit;

#[cfg(test)]
mod conversions;
