use crate::{NUM_ROUTED_WIRES, NUM_WIRES, HaloCurve};
use num::BigUint;
use std::marker::PhantomData;

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

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct AffinePointTarget<C: HaloCurve> {
    pub x: Target,
    pub y: Target,
    phantom: PhantomData<C>,
}

impl<C: HaloCurve> AffinePointTarget<C> {
    pub fn new(x: Target, y: Target) -> Self {
        AffinePointTarget { x, y, phantom: PhantomData }
    }

    pub fn to_vec(&self) -> Vec<Target> {
        vec![self.x, self.y]
    }
}

/// Represents a scalar * point multiplication operation.
pub struct CurveMulOp<C: HaloCurve> {
    pub scalar: Target,
    pub point: AffinePointTarget<C>,
    phantom: PhantomData<C>,
}

impl<C: HaloCurve> CurveMulOp<C> {
    pub fn new(scalar: Target, point: AffinePointTarget<C>) -> Self {
        CurveMulOp { scalar, point, phantom: PhantomData }
    }
}

pub struct CurveMulEndoResult<C: HaloCurve> {
    pub mul_result: AffinePointTarget<C>,
    pub actual_scalar: Target,
    pub(crate) phantom: PhantomData<C>,
}

pub struct CurveMsmEndoResult<C: HaloCurve> {
    pub msm_result: AffinePointTarget<C>,
    /// While `msm` computes a sum of `[s] P` terms, `msm_endo` computes a sum of `[n(s)] P` terms
    /// for some injective `n`. Here we return each `n(s)`, i.e., the scalar by which the point was
    /// actually multiplied.
    pub actual_scalars: Vec<Target>,
    pub(crate) phantom: PhantomData<C>,
}
