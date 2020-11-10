use std::convert::Infallible;
use std::hash::Hash;
use std::marker::PhantomData;
use std::rc::Rc;

pub use constraint_polynomial::*;
pub use gate::*;
pub use generator::*;
pub use prover::*;
pub use verifier::*;
pub use witness::*;

use crate::{Field, GateInstance, Gate2, Wire};

mod constraint_polynomial;
mod gate;
mod generator;
mod partitions;
mod prover;
mod witness;
mod verifier;

pub struct CircuitBuilder2<F: Field> {
    gates: Vec<Rc<dyn Gate2<F>>>,
    gate_instances: Vec<GateInstance<F>>,
}

impl<F: Field> CircuitBuilder2<F> {
    pub fn new() -> Self {
        CircuitBuilder2 {
            gates: Vec::new(),
            gate_instances: Vec::new(),
        }
    }

    /// Adds a gate to the circuit, and returns its index.
    pub fn add_gate(&mut self, gate_instance: GateInstance<F>) -> usize {
        let index = self.gate_instances.len();
        self.gate_instances.push(gate_instance);
        index
    }

    /// Shorthand for `generate_copy` and `assert_equal`.
    /// Both elements must be routable, otherwise this method will panic.
    pub fn copy(&mut self, x: Target2<F>, y: Target2<F>) {
        self.generate_copy(x, y);
        self.assert_equal(x, y);
    }

    /// Adds a generator which will copy `x` to `y`.
    pub fn generate_copy(&mut self, x: Target2<F>, y: Target2<F>) {
        todo!();
    }

    /// Uses Plonk's permutation argument to require that two elements be equal.
    /// Both elements must be routable, otherwise this method will panic.
    pub fn assert_equal(&mut self, x: Target2<F>, y: Target2<F>) {
        assert!(x.is_routable());
        assert!(y.is_routable());
    }

    /// Returns a routable target with a value of 0.
    pub fn zero(&mut self) -> Target2<F> {
        self.constant(F::ZERO)
    }

    /// Returns a routable target with a value of 1.
    pub fn one(&mut self) -> Target2<F> {
        self.constant(F::ONE)
    }

    /// Returns a routable target with a value of 2.
    pub fn two(&mut self) -> Target2<F> {
        self.constant(F::TWO)
    }

    /// Returns a routable target with a value of `ORDER - 1`.
    pub fn neg_one(&mut self) -> Target2<F> {
        self.constant(F::NEG_ONE)
    }

    /// Returns a routable target with the given constant value.
    pub fn constant(&mut self, c: F) -> Target2<F> {
        todo!()
    }
}

/// A location in the witness.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum Target2<F: Field> {
    Wire(Wire),
    PublicInput { index: usize },
    VirtualAdviceTarget { index: usize },
    // Trick taken from https://github.com/rust-lang/rust/issues/32739#issuecomment-627765543.
    _Field(Infallible, PhantomData<F>),
}

impl<F: Field> Target2<F> {
    pub fn wire(gate: usize, input: usize) -> Self {
        Self::Wire(Wire { gate, input })
    }

    pub fn is_routable(&self) -> bool {
        match self {
            Target2::Wire(wire) => wire.is_routable(),
            Target2::PublicInput { .. } => true,
            Target2::VirtualAdviceTarget { .. } => false,
            Target2::_Field(_, _) => unreachable!(),
        }
    }
}
