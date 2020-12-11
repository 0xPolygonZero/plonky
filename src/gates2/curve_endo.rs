use std::marker::PhantomData;

use crate::{CircuitConfig, ConstraintPolynomial, Curve, Field, Gate2, PartialWitness2, SimpleGenerator, Target2, Wire, WitnessGenerator2, HaloCurve};

/// Performs a step of Halo's accumulate-with-endomorphism loop.
pub struct CurveEndoGate2<C: HaloCurve> {
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> CurveEndoGate2<C> {
    pub const WIRE_ADDEND_X: usize = 0;
    pub const WIRE_ADDEND_Y: usize = 1;
    pub const WIRE_SCALAR_BIT_0: usize = 2;
    pub const WIRE_SCALAR_BIT_1: usize = 3;
    pub const WIRE_GROUP_ACC_X: usize = 4;
    pub const WIRE_GROUP_ACC_Y: usize = 5;
    pub const WIRE_INVERSE: usize = 6;
    pub const WIRE_LAMBDA: usize = 7;
}

impl<C: HaloCurve> Gate2<C::BaseField> for CurveEndoGate2<C> {
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
        let gen = CurveEndoGateGenerator::<C> { gate_index, _phantom: PhantomData };
        vec![Box::new(gen)]
    }
}

struct CurveEndoGateGenerator<C: HaloCurve> {
    gate_index: usize,
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> SimpleGenerator<C::BaseField> for CurveEndoGateGenerator<C> {
    fn dependencies(&self) -> Vec<Target2<C::BaseField>> {
        vec![
            Target2::Wire(Wire { gate: self.gate_index, input: CurveEndoGate2::<C>::WIRE_ADDEND_X }),
            Target2::Wire(Wire { gate: self.gate_index, input: CurveEndoGate2::<C>::WIRE_ADDEND_Y }),
            Target2::Wire(Wire { gate: self.gate_index, input: CurveEndoGate2::<C>::WIRE_SCALAR_BIT_0 }),
            Target2::Wire(Wire { gate: self.gate_index, input: CurveEndoGate2::<C>::WIRE_SCALAR_BIT_1 }),
            Target2::Wire(Wire { gate: self.gate_index, input: CurveEndoGate2::<C>::WIRE_GROUP_ACC_X }),
            Target2::Wire(Wire { gate: self.gate_index, input: CurveEndoGate2::<C>::WIRE_GROUP_ACC_Y }),
        ]
    }

    fn run_once(
        &mut self,
        witness: &PartialWitness2<C::BaseField>,
    ) -> PartialWitness2<C::BaseField> {
        let addend_x_wire = Wire { gate: self.gate_index, input: CurveEndoGate2::<C>::WIRE_ADDEND_X };
        let addend_y_wire = Wire { gate: self.gate_index, input: CurveEndoGate2::<C>::WIRE_ADDEND_Y };
        let scalar_bit_0_wire = Wire { gate: self.gate_index, input: CurveEndoGate2::<C>::WIRE_SCALAR_BIT_0 };
        let scalar_bit_1_wire = Wire { gate: self.gate_index, input: CurveEndoGate2::<C>::WIRE_SCALAR_BIT_1 };
        let group_acc_old_x_wire = Wire { gate: self.gate_index, input: CurveEndoGate2::<C>::WIRE_GROUP_ACC_X };
        let group_acc_old_y_wire = Wire { gate: self.gate_index, input: CurveEndoGate2::<C>::WIRE_GROUP_ACC_Y };
        let group_acc_new_x_wire = Wire { gate: self.gate_index + 1, input: CurveEndoGate2::<C>::WIRE_GROUP_ACC_X };
        let group_acc_new_y_wire = Wire { gate: self.gate_index + 1, input: CurveEndoGate2::<C>::WIRE_GROUP_ACC_Y };
        let inverse_wire = Wire { gate: self.gate_index, input: CurveEndoGate2::<C>::WIRE_INVERSE };
        let lambda_wire = Wire { gate: self.gate_index, input: CurveEndoGate2::<C>::WIRE_LAMBDA };

        // Load input values.
        let addend_x = witness.get_wire(addend_x_wire);
        let addend_y = witness.get_wire(addend_y_wire);
        let scalar_bit_0 = witness.get_wire(scalar_bit_0_wire);
        let scalar_bit_1 = witness.get_wire(scalar_bit_1_wire);
        let group_acc_old_x = witness.get_wire(group_acc_old_x_wire);
        let group_acc_old_y = witness.get_wire(group_acc_old_y_wire);

        // Compute S_i as defined in Halo.
        let mut s_i_x = addend_x;
        if scalar_bit_0 == C::BaseField::ONE {
            s_i_x = s_i_x * C::ZETA;
        }
        let mut s_i_y = addend_y;
        if scalar_bit_1 == C::BaseField::ZERO {
            s_i_y = -s_i_y;
        }

        // Compute group_acc_new = group_acc_old_x + s_i.
        let dx = group_acc_old_x - s_i_x;
        let dy = group_acc_old_y - s_i_y;
        let inverse = dx.multiplicative_inverse().expect("x_1 = x_2");
        let lambda = dy * inverse;
        let group_acc_new_x = lambda.square() - group_acc_old_x - s_i_x;
        let group_acc_new_y = lambda * (group_acc_old_x - group_acc_new_x) - group_acc_old_y;

        let mut result = PartialWitness2::new();
        result.set_wire(inverse_wire, inverse);
        result.set_wire(lambda_wire, lambda);
        result.set_wire(group_acc_new_x_wire, group_acc_new_x);
        result.set_wire(group_acc_new_y_wire, group_acc_new_y);
        result
    }
}
