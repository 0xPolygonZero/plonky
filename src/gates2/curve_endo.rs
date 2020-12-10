use std::marker::PhantomData;

use crate::{CircuitConfig, ConstraintPolynomial, Curve, Gate2, SimpleGenerator, WitnessGenerator2};

/// Performs a step of Halo's accumulate-with-endomorphism loop.
pub struct CurveEndoGate2<C: Curve> {
    _phantom: PhantomData<C>,
}

impl<C: Curve> CurveEndoGate2<C> {
    pub const WIRE_ADDEND_X: usize = 0;
    pub const WIRE_ADDEND_Y: usize = 1;
    pub const WIRE_SCALAR_BIT_0: usize = 2;
    pub const WIRE_SCALAR_BIT_1: usize = 3;
    pub const WIRE_GROUP_ACC_X: usize = 4;
    pub const WIRE_GROUP_ACC_Y: usize = 5;
    pub const WIRE_INVERSE: usize = 6;
}

impl<C: Curve> Gate2<C::BaseField> for CurveEndoGate2<C> {
    fn id(&self) -> String {
        "CurveEndoGate".into()
    }

    fn constraints(&self, _config: CircuitConfig) -> Vec<ConstraintPolynomial<C::BaseField>> {
        unimplemented!()
    }

    fn generators(
        &self,
        _config: CircuitConfig,
        gate_index: usize,
        _local_constants: Vec<C::BaseField>,
        _next_constants: Vec<<C as Curve>::BaseField>,
    ) -> Vec<Box<dyn WitnessGenerator2<C::BaseField>>> {
        unimplemented!()
    }
}

struct CurveEndoGateGenerator<C: Curve> {
    gate_index: usize,
    _phantom: PhantomData<C>,
}

impl<C: Curve> SimpleGenerator<C::BaseField> for CurveEndoGateGenerator<C> {}
