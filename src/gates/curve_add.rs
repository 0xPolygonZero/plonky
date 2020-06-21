use std::marker::PhantomData;

use crate::{CircuitBuilder, Curve, Field, HaloCurve, PartialWitness, Target, Wire, WitnessGenerator};
use crate::gates::{assert_binary_recursively, assert_inverses_recursively, Gate};

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

    const PREFIX: &'static [bool] = &[true, false, true, false, false, false];

    fn evaluate_unfiltered(
        _local_constant_values: &[InnerC::BaseField],
        local_wire_values: &[InnerC::BaseField],
        right_wire_values: &[InnerC::BaseField],
        _below_wire_values: &[InnerC::BaseField],
    ) -> Vec<InnerC::BaseField> {
        let x1 = local_wire_values[Self::WIRE_GROUP_ACC_X];
        let y1 = local_wire_values[Self::WIRE_GROUP_ACC_Y];
        let x3 = right_wire_values[Self::WIRE_GROUP_ACC_X];
        let y3 = right_wire_values[Self::WIRE_GROUP_ACC_Y];

        let scalar_acc_old = local_wire_values[Self::WIRE_SCALAR_ACC_OLD];
        let scalar_acc_new = local_wire_values[Self::WIRE_SCALAR_ACC_NEW];
        let x2 = local_wire_values[Self::WIRE_ADDEND_X];
        let y2 = local_wire_values[Self::WIRE_ADDEND_Y];
        let scalar_bit = local_wire_values[Self::WIRE_SCALAR_BIT];
        let inverse = local_wire_values[Self::WIRE_INVERSE];
        let lambda = local_wire_values[Self::WIRE_LAMBDA];

        let computed_lambda = (y1 - y2) * inverse;
        let computed_x3 = lambda.square() - x1 - x2;
        let computed_y3 = lambda * (x1 - x3) - y1;

        vec![
            computed_lambda - lambda,
            computed_x3 - x3,
            computed_y3 - y3,
            scalar_acc_new - scalar_acc_old.double() + scalar_bit,
            scalar_bit * (scalar_bit - InnerC::BaseField::ONE),
            inverse * (x1 - x2) - InnerC::BaseField::ONE,
        ]
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        _local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        _below_wire_values: &[Target],
    ) -> Vec<Target> {
        let x1 = local_wire_values[Self::WIRE_GROUP_ACC_X];
        let y1 = local_wire_values[Self::WIRE_GROUP_ACC_Y];
        let x3 = right_wire_values[Self::WIRE_GROUP_ACC_X];
        let y3 = right_wire_values[Self::WIRE_GROUP_ACC_Y];

        let scalar_acc_old = local_wire_values[Self::WIRE_SCALAR_ACC_OLD];
        let scalar_acc_new = local_wire_values[Self::WIRE_SCALAR_ACC_NEW];
        let x2 = local_wire_values[Self::WIRE_ADDEND_X];
        let y2 = local_wire_values[Self::WIRE_ADDEND_Y];
        let scalar_bit = local_wire_values[Self::WIRE_SCALAR_BIT];
        let inverse = local_wire_values[Self::WIRE_INVERSE];
        let lambda = local_wire_values[Self::WIRE_LAMBDA];

        let x1_minus_x2 = builder.sub(x1, x2);
        let x1_plus_x2 = builder.add(x1, x2);
        let x1_minus_x3 = builder.sub(x1, x3);
        let y1_minus_y2 = builder.sub(y1, y2);

        let computed_lambda = builder.mul(y1_minus_y2, inverse);
        let computed_x3 = builder.mul_sub(lambda, lambda, x1_plus_x2);
        let computed_y3 = builder.mul_sub(lambda, x1_minus_x3, y1);

        let double_scalar_acc_old = builder.double(scalar_acc_old);
        let computed_scalar_acc_new = builder.add(double_scalar_acc_old, scalar_bit);

        vec![
            builder.sub(computed_lambda, lambda),
            builder.sub(computed_x3, x3),
            builder.sub(computed_y3, y3),
            builder.sub(computed_scalar_acc_new, scalar_acc_new),
            assert_binary_recursively(builder, scalar_bit),
            assert_inverses_recursively(builder, inverse, x1_minus_x2),
        ]
    }
}

impl<C: HaloCurve, InnerC: Curve<BaseField = C::ScalarField>> WitnessGenerator<C::ScalarField>
for CurveAddGate<C, InnerC>
{
    fn dependencies(&self) -> Vec<Target> {
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
            Target::Wire(Wire {
                gate: self.index,
                input: Self::WIRE_SCALAR_BIT,
            }),
        ]
    }

    fn generate(
        &self,
        _constants: &Vec<Vec<C::ScalarField>>,
        witness: &PartialWitness<InnerC::BaseField>,
    ) -> PartialWitness<InnerC::BaseField> {
        let group_acc_old_x_target = Wire {
            gate: self.index,
            input: Self::WIRE_GROUP_ACC_X,
        };
        let group_acc_new_x_target = Wire {
            gate: self.index + 1,
            input: Self::WIRE_GROUP_ACC_X,
        };
        let group_acc_old_y_target = Wire {
            gate: self.index,
            input: Self::WIRE_GROUP_ACC_Y,
        };
        let group_acc_new_y_target = Wire {
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
        let addend_x_target = Wire {
            gate: self.index,
            input: Self::WIRE_ADDEND_X,
        };
        let addend_y_target = Wire {
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

        let group_acc_old_x = witness.get_wire(group_acc_old_x_target);
        let group_acc_old_y = witness.get_wire(group_acc_old_y_target);

        let scalar_acc_old = witness.get_wire(scalar_acc_old_target);

        let addend_x = witness.get_wire(addend_x_target);
        let addend_y = witness.get_wire(addend_y_target);

        let scalar_bit = witness.get_wire(scalar_bit_target);
        debug_assert!(scalar_bit.is_zero() || scalar_bit.is_one());

        let scalar_acc_new = scalar_acc_old.double() + scalar_bit;

        // Here's where our abstraction leaks a bit. Although we already have the sum, we need to
        // redo part of the computation in order to populate the purported inverse wire.
        let dx = group_acc_old_x - addend_x;
        let dy = group_acc_old_y - addend_y;
        let inverse = dx.multiplicative_inverse().expect("x_1 = x_2");
        let lambda = dy * inverse;

        let x_3 = lambda.square() - group_acc_old_x - addend_x;
        let y_3 = lambda * (group_acc_old_x - x_3) - group_acc_old_y;
        let (group_acc_new_x, group_acc_new_y) = if scalar_bit.is_one() {
            (x_3, y_3)
        } else {
            (group_acc_old_x, group_acc_old_y)
        };

        let mut result = PartialWitness::new();
        result.set_wire(group_acc_new_x_target, group_acc_new_x);
        result.set_wire(group_acc_new_y_target, group_acc_new_y);
        result.set_wire(scalar_acc_new_target, scalar_acc_new);
        result.set_wire(inverse_target, inverse);
        result.set_wire(lambda_target, lambda);
        result
    }
}

#[cfg(test)]
mod tests {
    use crate::{CurveAddGate, test_gate_low_degree, Tweedledee, Tweedledum};

    test_gate_low_degree!(
        low_degree_CurveAddGate,
        Tweedledum,
        CurveAddGate<Tweedledum, Tweedledee>
    );
}