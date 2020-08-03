use crate::{NUM_ROUTED_WIRES, NUM_WIRES};
use num::BigUint;

/// A sort of proxy wire, in the context of routing and witness generation. It is not an actual
/// witness element (i.e. wire) itself, but it can be copy-constrained to wires, listed as a
/// dependency in generators, etc.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct VirtualTarget {
    pub index: usize,
}

/// Represents a wire in the circuit.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct Wire {
    /// The index of the associated gate.
    pub gate: usize,
    /// The index of the gate input wherein this wire is inserted.
    pub input: usize,
}

impl Wire {
    pub fn is_routable(&self) -> bool {
        self.input < NUM_ROUTED_WIRES
    }
}

/// A routing target.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum Target {
    VirtualTarget(VirtualTarget),
    Wire(Wire),
}

#[derive(Clone)]
/// A `Target` with a (inclusive) known upper bound.
pub struct BoundedTarget {
    pub target: Target,
    /// An inclusive upper bound on this number.
    pub max: BigUint,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct PublicInput {
    pub index: usize,
}

/// See `PublicInputGate` for an explanation of how we make public inputs routable.
impl PublicInput {
    pub fn original_wire(&self) -> Wire {
        let gate = self.index / NUM_WIRES * 2;
        let input = self.index % NUM_WIRES;
        Wire { gate, input }
    }

    pub fn routable_target(&self) -> Target {
        let Wire {
            mut gate,
            mut input,
        } = self.original_wire();
        if input >= NUM_ROUTED_WIRES {
            gate += 1;
            input -= NUM_ROUTED_WIRES;
        }
        Target::Wire(Wire { gate, input })
    }
}
