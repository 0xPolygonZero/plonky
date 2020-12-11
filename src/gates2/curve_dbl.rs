use std::marker::PhantomData;

use crate::{CircuitConfig, ConstraintPolynomial, Curve, Field, Gate2, PartialWitness2, SimpleGenerator, Target2, Wire, WitnessGenerator2};

pub struct CurveDblGate2<C: Curve> {
    _phantom: PhantomData<C>,
}

impl<C: Curve> CurveDblGate2<C> {
    pub const WIRE_X_OLD: usize = 0;
    pub const WIRE_Y_OLD: usize = 1;
    pub const WIRE_X_NEW: usize = 2;
    pub const WIRE_Y_NEW: usize = 3;
    pub const WIRE_INVERSE: usize = 4;
    pub const WIRE_LAMBDA: usize = 5;
}

impl<C: Curve> Gate2<C::BaseField> for CurveDblGate2<C> {
    fn id(&self) -> String {
        "CurveDblGate".into()
    }

    fn constraints(&self, _config: CircuitConfig) -> Vec<ConstraintPolynomial<C::BaseField>> {
        unimplemented!()
    }

    fn generators(
        &self,
        _config: CircuitConfig,
        gate_index: usize,
        _local_constants: Vec<C::BaseField>,
        _next_constants: Vec<C::BaseField>,
    ) -> Vec<Box<dyn WitnessGenerator2<C::BaseField>>> {
        let gen = CurveDblGateGenerator::<C> { gate_index, _phantom: PhantomData };
        vec![Box::new(gen)]
    }
}

struct CurveDblGateGenerator<C: Curve> {
    gate_index: usize,
    _phantom: PhantomData<C>,
}

impl<C: Curve> SimpleGenerator<C::BaseField> for CurveDblGateGenerator<C> {
    fn dependencies(&self) -> Vec<Target2<C::BaseField>> {
        vec![
            Target2::Wire(Wire { gate: self.gate_index, input: CurveDblGate2::<C>::WIRE_X_OLD }),
            Target2::Wire(Wire { gate: self.gate_index, input: CurveDblGate2::<C>::WIRE_Y_OLD }),
        ]
    }

    fn run_once(
        &mut self,
        witness: &PartialWitness2<C::BaseField>,
    ) -> PartialWitness2<C::BaseField> {
        let x_old_wire = Wire { gate: self.gate_index, input: CurveDblGate2::<C>::WIRE_X_OLD };
        let y_old_wire = Wire { gate: self.gate_index, input: CurveDblGate2::<C>::WIRE_Y_OLD };
        let x_new_wire = Wire { gate: self.gate_index, input: CurveDblGate2::<C>::WIRE_X_NEW };
        let y_new_wire = Wire { gate: self.gate_index, input: CurveDblGate2::<C>::WIRE_Y_NEW };
        let inverse_wire = Wire { gate: self.gate_index, input: CurveDblGate2::<C>::WIRE_INVERSE };
        let lambda_wire = Wire { gate: self.gate_index, input: CurveDblGate2::<C>::WIRE_LAMBDA };

        let x_old = witness.get_wire(x_old_wire);
        let y_old = witness.get_wire(y_old_wire);
        let inverse = y_old.double().multiplicative_inverse().expect("y = 0");
        let lambda = x_old.square().triple() * inverse;
        let x_new = lambda.square() - x_old.double();
        let y_new = lambda * (x_old - x_new) - y_old;

        let mut result = PartialWitness2::new();
        result.set_wire(inverse_wire, inverse);
        result.set_wire(lambda_wire, lambda);
        result.set_wire(x_new_wire, x_new);
        result.set_wire(y_new_wire, y_new);
        result
    }
}
