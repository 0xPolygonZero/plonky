use serde::export::PhantomData;

use crate::{CircuitConfig, ConstraintPolynomial, Curve, Field, Gate2, PartialWitness2, SimpleGenerator, Target2, Wire, WitnessGenerator2};

pub struct CurveAddGate2<C: Curve> {
    _phantom: PhantomData<C>,
}

impl<C: Curve> CurveAddGate2<C> {
    pub const WIRE_GROUP_ACC_X: usize = 0;
    pub const WIRE_GROUP_ACC_Y: usize = 1;
    pub const WIRE_SCALAR_ACC_OLD: usize = 2;
    pub const WIRE_SCALAR_ACC_NEW: usize = 3;
    pub const WIRE_ADDEND_X: usize = 4;
    pub const WIRE_ADDEND_Y: usize = 5;
    pub const WIRE_SCALAR_BIT: usize = 6;
    pub const WIRE_INVERSE: usize = 7;
    pub const WIRE_LAMBDA: usize = 8;
}

impl<C: Curve> Gate2<C::BaseField> for CurveAddGate2<C> {
    fn id(&self) -> String {
        "CurveAddGate".into()
    }

    fn constraints(
        &self,
        config: CircuitConfig,
    ) -> Vec<ConstraintPolynomial<<C as Curve>::BaseField>> {
        // Notation:
        // - p1 is the accumulator;
        // - p2 is the addend;
        // - p3 = p1 + p2;
        // - p4 = if scalar_bit { p3 } else { p1 }

        let x1 = ConstraintPolynomial::<C::BaseField>::local_wire_value(Self::WIRE_GROUP_ACC_X);
        let y1 = ConstraintPolynomial::<C::BaseField>::local_wire_value(Self::WIRE_GROUP_ACC_Y);
        let x4 = ConstraintPolynomial::<C::BaseField>::next_wire_value(Self::WIRE_GROUP_ACC_X);
        let y4 = ConstraintPolynomial::<C::BaseField>::next_wire_value(Self::WIRE_GROUP_ACC_Y);
        let scalar_acc_old = ConstraintPolynomial::<C::BaseField>::local_wire_value(Self::WIRE_SCALAR_ACC_OLD);
        let scalar_acc_new = ConstraintPolynomial::<C::BaseField>::local_wire_value(Self::WIRE_SCALAR_ACC_NEW);
        let x2 = ConstraintPolynomial::<C::BaseField>::local_wire_value(Self::WIRE_ADDEND_X);
        let y2 = ConstraintPolynomial::<C::BaseField>::local_wire_value(Self::WIRE_ADDEND_Y);
        let scalar_bit = ConstraintPolynomial::<C::BaseField>::local_wire_value(Self::WIRE_SCALAR_BIT);
        let inverse = ConstraintPolynomial::<C::BaseField>::local_wire_value(Self::WIRE_INVERSE);
        let lambda = ConstraintPolynomial::<C::BaseField>::local_wire_value(Self::WIRE_LAMBDA);

        let computed_lambda = (&y1 - &y2) * &inverse;
        let x3 = lambda.square() - &x1 - &x2;
        // We subtract x4 instead of x3 in order to minimize degree. This will give an incorrect
        // result for y3 if x3 != x4, which happens when scalar_bit = 0, but in that case y3 will
        // be ignored (i.e. multiplied by zero), so we're okay.
        let y3 = &lambda * (&x1 - &x4) - &y1;

        let not_scalar_bit = ConstraintPolynomial::constant_usize(1) - &scalar_bit;
        let computed_x4 = &scalar_bit * &x3 + &not_scalar_bit * &x1;
        let computed_y4 = &scalar_bit * &y3 + &not_scalar_bit * &y1;

        vec![
            &computed_lambda - &lambda,
            &computed_x4 - &x4,
            &computed_y4 - &y4,
            &scalar_acc_new - (scalar_acc_old.double() + &scalar_bit),
            &scalar_bit * &not_scalar_bit,
            &inverse * (&x1 - &x2) - C::BaseField::ONE,
        ]
    }

    fn generators(
        &self,
        config: CircuitConfig,
        gate_index: usize,
        local_constants: Vec<<C as Curve>::BaseField>,
        next_constants: Vec<<C as Curve>::BaseField>,
    ) -> Vec<Box<dyn WitnessGenerator2<C::BaseField>>> {
        let gen = CurveAddGateGenerator::<C> { gate_index, _phantom: PhantomData };
        vec![Box::new(gen)]
    }
}

struct CurveAddGateGenerator<C: Curve> {
    gate_index: usize,
    _phantom: PhantomData<C>,
}

impl<C: Curve> SimpleGenerator<C::BaseField> for CurveAddGateGenerator<C> {
    fn dependencies(&self) -> Vec<Target2<C::BaseField>> {
        vec![
            Target2::Wire(Wire {
                gate: self.gate_index,
                input: CurveAddGate2::<C>::WIRE_GROUP_ACC_X,
            }),
            Target2::Wire(Wire {
                gate: self.gate_index,
                input: CurveAddGate2::<C>::WIRE_GROUP_ACC_Y,
            }),
            Target2::Wire(Wire {
                gate: self.gate_index,
                input: CurveAddGate2::<C>::WIRE_SCALAR_ACC_OLD,
            }),
            Target2::Wire(Wire {
                gate: self.gate_index,
                input: CurveAddGate2::<C>::WIRE_ADDEND_X,
            }),
            Target2::Wire(Wire {
                gate: self.gate_index,
                input: CurveAddGate2::<C>::WIRE_ADDEND_Y,
            }),
        ]
    }

    fn run_once(
        &mut self,
        witness: &PartialWitness2<C::BaseField>,
    ) -> PartialWitness2<C::BaseField> {
        // Notation:
        // - p1 is the accumulator;
        // - p2 is the addend;
        // - p3 = p1 + p2;
        // - p4 = if scalar_bit { p3 } else { p1 }

        let x1_target = Wire {
            gate: self.gate_index,
            input: CurveAddGate2::<C>::WIRE_GROUP_ACC_X,
        };
        let y1_target = Wire {
            gate: self.gate_index,
            input: CurveAddGate2::<C>::WIRE_GROUP_ACC_Y,
        };
        let x4_target = Wire {
            gate: self.gate_index + 1,
            input: CurveAddGate2::<C>::WIRE_GROUP_ACC_X,
        };
        let y4_target = Wire {
            gate: self.gate_index + 1,
            input: CurveAddGate2::<C>::WIRE_GROUP_ACC_Y,
        };
        let scalar_acc_old_target = Wire {
            gate: self.gate_index,
            input: CurveAddGate2::<C>::WIRE_SCALAR_ACC_OLD,
        };
        let scalar_acc_new_target = Wire {
            gate: self.gate_index,
            input: CurveAddGate2::<C>::WIRE_SCALAR_ACC_NEW,
        };
        let x2_target = Wire {
            gate: self.gate_index,
            input: CurveAddGate2::<C>::WIRE_ADDEND_X,
        };
        let y2_target = Wire {
            gate: self.gate_index,
            input: CurveAddGate2::<C>::WIRE_ADDEND_Y,
        };
        let scalar_bit_target = Wire {
            gate: self.gate_index,
            input: CurveAddGate2::<C>::WIRE_SCALAR_BIT,
        };
        let inverse_target = Wire {
            gate: self.gate_index,
            input: CurveAddGate2::<C>::WIRE_INVERSE,
        };
        let lambda_target = Wire {
            gate: self.gate_index,
            input: CurveAddGate2::<C>::WIRE_LAMBDA,
        };

        let x1 = witness.get_wire(x1_target);
        let y1 = witness.get_wire(y1_target);

        let scalar_acc_old = witness.get_wire(scalar_acc_old_target);

        let x2 = witness.get_wire(x2_target);
        let y2 = witness.get_wire(y2_target);

        let scalar_bit = witness.get_wire(scalar_bit_target);
        debug_assert!(scalar_bit.is_zero() || scalar_bit.is_one());

        let scalar_acc_new = scalar_acc_old.double() + scalar_bit;

        let dx = x1 - x2;
        let dy = y1 - y2;
        let inverse = dx.multiplicative_inverse().expect("x_1 = x_2");
        let lambda = dy * inverse;
        let x3 = lambda.square() - x1 - x2;
        let y3 = lambda * (x1 - x3) - y1;

        let (x4, y4) = if scalar_bit.is_one() {
            (x3, y3)
        } else {
            (x1, y1)
        };

        let mut result = PartialWitness2::new();
        result.set_wire(x4_target, x4);
        result.set_wire(y4_target, y4);
        result.set_wire(scalar_acc_new_target, scalar_acc_new);
        result.set_wire(inverse_target, inverse);
        result.set_wire(lambda_target, lambda);
        result
    }
}
