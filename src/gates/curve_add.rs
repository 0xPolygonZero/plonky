use std::marker::PhantomData;

use crate::gates::{assert_binary_recursively, assert_inverses_recursively, Gate};
use crate::{CircuitBuilder, Curve, Field, HaloCurve, PartialWitness, Target, Wire, WitnessGenerator};

/// A gate which performs incomplete point addition, conditioned on an input bit. In order to
/// facilitate MSMs which use this gate, it also adds the bit to an accumulator.
///
/// `C` is the curve whose points are being added.
pub struct CurveAddGate<C: HaloCurve, InnerC: Curve<BaseField = C::ScalarField>> {
    pub index: usize,
    _phantom_oc: PhantomData<C>,
    _phantom_ic: PhantomData<InnerC>,
}

impl<C: HaloCurve, InnerC: Curve<BaseField = C::ScalarField>> CurveAddGate<C, InnerC> {
    pub fn new(index: usize) -> Self {
        CurveAddGate {
            index,
            _phantom_oc: PhantomData,
            _phantom_ic: PhantomData,
        }
    }

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

impl<C: HaloCurve, InnerC: Curve<BaseField = C::ScalarField>> Gate<C> for CurveAddGate<C, InnerC> {
    const NAME: &'static str = "CurveAddGate";

    const PREFIX: &'static [bool] = &[true, false, true, false, true];

    fn evaluate_unfiltered(
        _local_constant_values: &[InnerC::BaseField],
        local_wire_values: &[InnerC::BaseField],
        right_wire_values: &[InnerC::BaseField],
        _below_wire_values: &[InnerC::BaseField],
    ) -> Vec<InnerC::BaseField> {
        // Notation:
        // - p1 is the accumulator;
        // - p2 is the addend;
        // - p3 = p1 + p2;
        // - p4 = if scalar_bit { p3 } else { p1 }

        let x1 = local_wire_values[Self::WIRE_GROUP_ACC_X];
        let y1 = local_wire_values[Self::WIRE_GROUP_ACC_Y];
        let x4 = right_wire_values[Self::WIRE_GROUP_ACC_X];
        let y4 = right_wire_values[Self::WIRE_GROUP_ACC_Y];

        let scalar_acc_old = local_wire_values[Self::WIRE_SCALAR_ACC_OLD];
        let scalar_acc_new = local_wire_values[Self::WIRE_SCALAR_ACC_NEW];
        let x2 = local_wire_values[Self::WIRE_ADDEND_X];
        let y2 = local_wire_values[Self::WIRE_ADDEND_Y];
        let scalar_bit = local_wire_values[Self::WIRE_SCALAR_BIT];
        let inverse = local_wire_values[Self::WIRE_INVERSE];
        let lambda = local_wire_values[Self::WIRE_LAMBDA];

        let computed_lambda = (y1 - y2) * inverse;
        let x3 = lambda.square() - x1 - x2;
        // We subtract x4 instead of x3 in order to minimize degree. This will give an incorrect
        // result for y3 if x3 != x4, which happens when scalar_bit = 0, but in that case y3 will
        // be ignored (i.e. multiplied by zero), so we're okay.
        let y3 = lambda * (x1 - x4) - y1;

        let not_scalar_bit = InnerC::BaseField::ONE - scalar_bit;
        let computed_x4 = scalar_bit * x3 + not_scalar_bit * x1;
        let computed_y4 = scalar_bit * y3 + not_scalar_bit * y1;

        vec![
            computed_lambda - lambda,
            computed_x4 - x4,
            computed_y4 - y4,
            scalar_acc_new - (scalar_acc_old.double() + scalar_bit),
            scalar_bit * not_scalar_bit,
            inverse * (x1 - x2) - InnerC::BaseField::ONE,
        ]
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        _local_constant_values: &[Target<C::ScalarField>],
        local_wire_values: &[Target<C::ScalarField>],
        right_wire_values: &[Target<C::ScalarField>],
        _below_wire_values: &[Target<C::ScalarField>],
    ) -> Vec<Target<C::ScalarField>> {
        // Notation:
        // - p1 is the accumulator;
        // - p2 is the addend;
        // - p3 = p1 + p2;
        // - p4 = if scalar_bit { p3 } else { p1 }

        let x1 = local_wire_values[Self::WIRE_GROUP_ACC_X];
        let y1 = local_wire_values[Self::WIRE_GROUP_ACC_Y];
        let x4 = right_wire_values[Self::WIRE_GROUP_ACC_X];
        let y4 = right_wire_values[Self::WIRE_GROUP_ACC_Y];

        let scalar_acc_old = local_wire_values[Self::WIRE_SCALAR_ACC_OLD];
        let scalar_acc_new = local_wire_values[Self::WIRE_SCALAR_ACC_NEW];
        let x2 = local_wire_values[Self::WIRE_ADDEND_X];
        let y2 = local_wire_values[Self::WIRE_ADDEND_Y];
        let scalar_bit = local_wire_values[Self::WIRE_SCALAR_BIT];
        let inverse = local_wire_values[Self::WIRE_INVERSE];
        let lambda = local_wire_values[Self::WIRE_LAMBDA];

        let x1_minus_x2 = builder.sub(x1, x2);
        let x1_plus_x2 = builder.add(x1, x2);
        let y1_minus_y2 = builder.sub(y1, y2);

        let computed_lambda = builder.mul(y1_minus_y2, inverse);
        let x3 = builder.mul_sub(lambda, lambda, x1_plus_x2);
        let x1_minus_x4 = builder.sub(x1, x4);
        // We subtract x4 instead of x3 in order to minimize degree. This will give an incorrect
        // result for y3 if x3 != x4, which happens when scalar_bit = 0, but in that case y3 will
        // be ignored (i.e. multiplied by zero), so we're okay.
        let y3 = builder.mul_sub(lambda, x1_minus_x4, y1);

        let not_scalar_bit = builder.not(scalar_bit);
        let x1_conditioned = builder.mul(x1, not_scalar_bit);
        let y1_conditioned = builder.mul(y1, not_scalar_bit);
        let computed_x4 = builder.mul_add(scalar_bit, x3, x1_conditioned);
        let computed_y4 = builder.mul_add(scalar_bit, y3, y1_conditioned);

        let double_scalar_acc_old = builder.double(scalar_acc_old);
        let computed_scalar_acc_new = builder.add(double_scalar_acc_old, scalar_bit);

        vec![
            builder.sub(computed_lambda, lambda),
            builder.sub(computed_x4, x4),
            builder.sub(computed_y4, y4),
            builder.sub(computed_scalar_acc_new, scalar_acc_new),
            assert_binary_recursively(builder, scalar_bit),
            assert_inverses_recursively(builder, inverse, x1_minus_x2),
        ]
    }
}

impl<C: HaloCurve, InnerC: Curve<BaseField = C::ScalarField>> WitnessGenerator<C::ScalarField>
    for CurveAddGate<C, InnerC>
{
    fn dependencies(&self) -> Vec<Target<C::ScalarField>> {
        vec![
            Target::Wire(Wire {
                gate: self.index,
                input: Self::WIRE_GROUP_ACC_X,
            }),
            Target::Wire(Wire {
                gate: self.index,
                input: Self::WIRE_GROUP_ACC_Y,
            }),
            Target::Wire(Wire {
                gate: self.index,
                input: Self::WIRE_SCALAR_ACC_OLD,
            }),
            Target::Wire(Wire {
                gate: self.index,
                input: Self::WIRE_ADDEND_X,
            }),
            Target::Wire(Wire {
                gate: self.index,
                input: Self::WIRE_ADDEND_Y,
            }),
        ]
    }

    fn generate(
        &self,
        _constants: &Vec<Vec<C::ScalarField>>,
        witness: &PartialWitness<InnerC::BaseField>,
    ) -> PartialWitness<InnerC::BaseField> {
        // Notation:
        // - p1 is the accumulator;
        // - p2 is the addend;
        // - p3 = p1 + p2;
        // - p4 = if scalar_bit { p3 } else { p1 }

        let x1_target = Wire {
            gate: self.index,
            input: Self::WIRE_GROUP_ACC_X,
        };
        let x4_target = Wire {
            gate: self.index + 1,
            input: Self::WIRE_GROUP_ACC_X,
        };
        let y1_target = Wire {
            gate: self.index,
            input: Self::WIRE_GROUP_ACC_Y,
        };
        let y4_target = Wire {
            gate: self.index + 1,
            input: Self::WIRE_GROUP_ACC_Y,
        };
        let scalar_acc_old_target = Wire {
            gate: self.index,
            input: Self::WIRE_SCALAR_ACC_OLD,
        };
        let scalar_acc_new_target = Wire {
            gate: self.index,
            input: Self::WIRE_SCALAR_ACC_NEW,
        };
        let x2_target = Wire {
            gate: self.index,
            input: Self::WIRE_ADDEND_X,
        };
        let y2_target = Wire {
            gate: self.index,
            input: Self::WIRE_ADDEND_Y,
        };
        let scalar_bit_target = Wire {
            gate: self.index,
            input: Self::WIRE_SCALAR_BIT,
        };
        let inverse_target = Wire {
            gate: self.index,
            input: Self::WIRE_INVERSE,
        };
        let lambda_target = Wire {
            gate: self.index,
            input: Self::WIRE_LAMBDA,
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

        let mut result = PartialWitness::new();
        result.set_wire(x4_target, x4);
        result.set_wire(y4_target, y4);
        result.set_wire(scalar_acc_new_target, scalar_acc_new);
        result.set_wire(inverse_target, inverse);
        result.set_wire(lambda_target, lambda);
        result
    }
}

#[cfg(test)]
mod tests {
    use crate::{test_gate_low_degree, CurveAddGate, Tweedledee, Tweedledum};

    test_gate_low_degree!(
        low_degree_CurveAddGate,
        Tweedledum,
        CurveAddGate<Tweedledum, Tweedledee>
    );
}
