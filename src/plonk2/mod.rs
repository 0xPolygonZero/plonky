use std::collections::HashSet;
use std::convert::Infallible;
use std::hash::Hash;
use std::marker::PhantomData;

pub use constraint_polynomial::*;
pub use gate::*;
pub use generator::*;
pub use prover::*;
pub use verifier::*;
pub use witness::*;

use crate::{Field, GateInstance, GateRef, Wire};

mod constraint_polynomial;
mod gate;
mod generator;
mod partitions;
mod prover;
mod witness;
mod verifier;

#[derive(Copy, Clone)]
pub struct CircuitConfig {
    pub num_wires: usize,
    pub num_routed_wires: usize,
}

impl CircuitConfig {
    pub fn advice_wires(&self) -> usize {
        self.num_wires - self.num_routed_wires
    }
}

pub struct CircuitBuilder2<F: Field> {
    config: CircuitConfig,
    gates: HashSet<GateRef<F>>,
    gate_instances: Vec<GateInstance<F>>,
    generators: Vec<Box<dyn WitnessGenerator2<F>>>,
}

impl<F: Field> CircuitBuilder2<F> {
    pub fn new(config: CircuitConfig) -> Self {
        CircuitBuilder2 {
            config,
            gates: HashSet::new(),
            gate_instances: Vec::new(),
            generators: Vec::new(),
        }
    }

    /// Adds a gate to the circuit, and returns its index.
    pub fn add_gate(&mut self, gate_instance: GateInstance<F>) -> usize {
        // If we haven't seen a gate of this type before, check that it's compatible with our
        // circuit configuration, then register it.
        if !self.gates.contains(&gate_instance.gate_type) {
            let gate = gate_instance.gate_type.clone();
            self.check_gate_compatibility(&gate);
            self.gates.insert(gate);
        }

        let index = self.gate_instances.len();
        self.gate_instances.push(gate_instance);
        index
    }

    fn check_gate_compatibility(&self, gate: &GateRef<F>) {
        assert!(gate.0.min_wires(self.config) <= self.config.num_wires);
    }

    /// Shorthand for `generate_copy` and `assert_equal`.
    /// Both elements must be routable, otherwise this method will panic.
    pub fn route(&mut self, src: Target2<F>, dst: Target2<F>) {
        self.generate_copy(src, dst);
        self.assert_equal(src, dst);
    }

    /// Adds a generator which will copy `src` to `dst`.
    pub fn generate_copy(&mut self, src: Target2<F>, dst: Target2<F>) {
        self.add_generator(CopyGenerator { src, dst });
    }

    /// Uses Plonk's permutation argument to require that two elements be equal.
    /// Both elements must be routable, otherwise this method will panic.
    pub fn assert_equal(&mut self, x: Target2<F>, y: Target2<F>) {
        assert!(x.is_routable());
        assert!(y.is_routable());
    }

    pub fn add_generator<G: WitnessGenerator2<F>>(&mut self, generator: G) {
        self.generators.push(Box::new(generator));
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
