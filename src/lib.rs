// Since we work with elliptic curve groups etc., nearly all the operations are suspicious to
// Clippy.
#![allow(clippy::suspicious_arithmetic_impl)]
// We have tons of bigint literals in Montgomery form, which won't be readable with or without underscores.
#![allow(clippy::unreadable_literal)]
#![feature(associated_type_bounds)]
// This is annoying and often wrong.
#![allow(clippy::needless_range_loop)]

// Required for generic low-level functions on small arrays.
#![feature(const_generics)]
// Unfortunatly it makes rustc complain, so we include
#![allow(incomplete_features)]


pub use bigint::*;
pub use circuit_bigint::*;
pub use circuit_builder::*;
pub use circuit_curve::*;
pub use circuit_foreign_field::*;
pub use circuit_ordering::*;
pub use conversions::*;
pub use curve::*;
pub use fft::*;
pub use field::*;
pub use gates::*;
pub use hash_to_curve::*;
pub use mds::*;
pub use partition::*;
pub use plonk::*;
pub use plonk_proof::*;
pub use plonk_recursion::*;
pub use poly_commit::*;
pub use polynomial::*;
pub use pseudorandom::*;
pub use rescue::*;
pub use serialization::*;
pub use target::*;
pub use verifier::*;
pub use witness::*;

mod bigint;
mod circuit_bigint;
mod circuit_builder;
mod circuit_curve;
mod circuit_foreign_field;
mod circuit_ordering;
mod conversions;
mod curve;
mod fft;
mod field;
mod gates;
pub mod halo;
mod hash_to_curve;
mod mds;
mod partition;
mod plonk;
pub mod plonk_challenger;
mod plonk_proof;
mod plonk_recursion;
pub mod plonk_util;
pub mod poly_commit;
pub mod polynomial;
mod pseudorandom;
mod rescue;
mod serialization;
mod target;
pub mod util;
mod verifier;
mod witness;

#[macro_use]
extern crate log;
