// Since we work with elliptic curve groups etc., nearly all the operations are suspicious to
// Clippy.
#![allow(clippy::suspicious_arithmetic_impl)]

// We have tons of bigint literals in Montgomery form, which won't be readable with or without underscores.
#![allow(clippy::unreadable_literal)]

// Constant equality checks are usually due to macro-based loop unrolling.
#![warn(clippy::eq_op)]

pub use fft::*;
pub use field::*;
pub use group::*;
pub use group_adds::*;
pub use group_msm::*;
pub use group_multiplication::*;
pub use group_summations::*;
pub use plonk::*;
pub use poly_commit::*;

mod conversions;
mod field;
mod fft;
mod group;
mod group_adds;
mod group_msm;
mod group_multiplication;
mod group_summations;
mod plonk;
mod poly_commit;
