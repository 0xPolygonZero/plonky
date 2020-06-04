// Since we work with elliptic curve groups etc., nearly all the operations are suspicious to
// Clippy.
#![allow(clippy::suspicious_arithmetic_impl)]
// We have tons of bigint literals in Montgomery form, which won't be readable with or without underscores.
#![allow(clippy::unreadable_literal)]

pub use bigint_arithmetic::*;
pub use bigint_inverse::*;
pub use bls12_377_base::*;
pub use bls12_377_curve::*;
pub use bls12_377_scalar::*;
pub use conversions::*;
pub use curve::*;
pub use curve_adds::*;
pub use curve_msm::*;
pub use curve_multiplication::*;
pub use curve_summations::*;
pub use fft::*;
pub use field::*;
pub use hash_to_curve::*;
pub use mds::*;
pub use plonk::*;
pub use plonk_gates::*;
pub use plonk_proof::*;
pub use plonk_recursion::*;
// pub use poly_commit::*;
pub use circuit_builder::*;
pub use partition::*;
pub use polynomial::*;
pub use pseudorandom::*;
pub use rescue::*;
pub use serialization::*;
pub use target::*;
pub use tweedledee_base::*;
pub use tweedledee_curve::*;
pub use tweedledum_base::*;
pub use tweedledum_curve::*;
pub use witness::*;

mod bigint_arithmetic;
mod bigint_inverse;
mod bls12_377_base;
mod bls12_377_curve;
mod bls12_377_scalar;
mod conversions;
mod curve;
mod curve_adds;
mod curve_msm;
mod curve_multiplication;
mod curve_summations;
mod fft;
mod field;
mod hash_to_curve;
mod mds;
mod plonk;
mod plonk_challenger;
mod plonk_gates;
mod plonk_proof;
mod plonk_recursion;
mod plonk_util;
// mod poly_commit;
mod circuit_builder;
mod partition;
mod polynomial;
mod pseudorandom;
mod rescue;
mod serialization;
mod target;
mod tweedledee_base;
mod tweedledee_curve;
mod tweedledum_base;
mod tweedledum_curve;
mod util;
mod witness;
