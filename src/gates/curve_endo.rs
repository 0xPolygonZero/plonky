use std::marker::PhantomData;

use crate::{AffinePoint, CircuitBuilder, Field, GRID_WIDTH, HaloCurve, PartialWitness, Target, Wire, WitnessGenerator};
use crate::gates::{assert_binary_recursively, assert_inverses_recursively, Gate};

/// A gate which performs an iteration of an simultaneous doubling MSM loop, employing the
/// endomorphism described in the Halo paper. `C` is the curve of the inner proof.
pub struct CurveEndoGate<C: HaloCurve, InnerC: HaloCurve<BaseField = C::ScalarField>> {
    pub index: usize,
    _phantom_oc: PhantomData<C>,
    _phantom_ic: PhantomData<InnerC>,
}

impl<C: HaloCurve, InnerC: HaloCurve<BaseField = C::ScalarField>> CurveEndoGate<C, InnerC> {
    pub fn new(index: usize) -> Self {
        CurveEndoGate {
            index,
            _phantom_oc: PhantomData,
            _phantom_ic: PhantomData,
        }
    }

    pub const WIRE_GROUP_ACC_X: usize = 0;
    pub const WIRE_GROUP_ACC_Y: usize = 1;
    pub const WIRE_SCALAR_ACC_UNSIGNED: usize = 2;
    pub const WIRE_SCALAR_ACC_SIGNED: usize = 3;
    pub const WIRE_ADDEND_X: usize = 4;
    pub const WIRE_ADDEND_Y: usize = 5;
    pub const WIRE_SCALAR_BIT_0: usize = 6;
    pub const WIRE_SCALAR_BIT_1: usize = 7;
    pub const WIRE_INVERSE: usize = 8;
}

impl<C: HaloCurve, InnerC: HaloCurve<BaseField = C::ScalarField>> Gate<C>
for CurveEndoGate<C, InnerC>
{
    const NAME: &'static str = "CurveEndoGate";

    const PREFIX: &'static [bool] = &[true, true];

    fn evaluate_unfiltered(
        _local_constant_values: &[InnerC::BaseField],
        local_wire_values: &[InnerC::BaseField],
        right_wire_values: &[InnerC::BaseField],
        below_wire_values: &[InnerC::BaseField],
    ) -> Vec<InnerC::BaseField> {
        let one = InnerC::BaseField::ONE;

        let x1 = local_wire_values[Self::WIRE_GROUP_ACC_X];
        let y1 = local_wire_values[Self::WIRE_GROUP_ACC_Y];
        let x_in = local_wire_values[Self::WIRE_ADDEND_X];
        let y_in = local_wire_values[Self::WIRE_ADDEND_Y];
        let x3 = right_wire_values[Self::WIRE_GROUP_ACC_X];
        let y3 = right_wire_values[Self::WIRE_GROUP_ACC_Y];
        let scalar_acc_unsigned_old = local_wire_values[Self::WIRE_SCALAR_ACC_UNSIGNED];
        let scalar_acc_unsigned_new = below_wire_values[Self::WIRE_SCALAR_ACC_UNSIGNED];
        let scalar_acc_signed_old = local_wire_values[Self::WIRE_SCALAR_ACC_SIGNED];
        let scalar_acc_signed_new = below_wire_values[Self::WIRE_SCALAR_ACC_SIGNED];
        let scalar_bit_0 = local_wire_values[Self::WIRE_SCALAR_BIT_0];
        let scalar_bit_1 = local_wire_values[Self::WIRE_SCALAR_BIT_1];
        let inverse = local_wire_values[Self::WIRE_INVERSE];

        // Conditionally apply the endo and conditionally negate in order to get (x2, y2), which is
        // the actual point we want to add to the accumulator.
        let x2 = ((InnerC::ZETA - one) * scalar_bit_1 + one) * x_in;
        let y2 = (scalar_bit_0.double() - one) * y_in;

        let lambda = (y1 - y2) * inverse;
        let computed_x3 = lambda.square() - x1 - x2;
        let computed_y3 = lambda * (x1 - x3) - y1;

        // This is based on Algorithm 2 in the Halo paper.
        let signed_limb_multiplier = (InnerC::ZETA - one) * scalar_bit_1 + one;
        let signed_limb = (scalar_bit_0.double() - one) * signed_limb_multiplier;

        vec![
            computed_x3 - x3,
            computed_y3 - y3,
            scalar_acc_unsigned_new - (scalar_acc_unsigned_old.quadruple()
                + scalar_bit_1.double()
                + scalar_bit_0),
            scalar_acc_signed_new - (scalar_acc_signed_old.double() + signed_limb),
            scalar_bit_0 * (scalar_bit_0 - one),
            scalar_bit_1 * (scalar_bit_1 - one),
            inverse * (x1 - x2) - one,
        ]
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        _local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        below_wire_values: &[Target],
    ) -> Vec<Target> {
        let one = builder.one_wire();
        let two = builder.two_wire();
        let four = builder.constant_wire_u32(4);
        let zeta = builder.constant_wire(InnerC::ZETA);

        let x1 = local_wire_values[Self::WIRE_GROUP_ACC_X];
        let y1 = local_wire_values[Self::WIRE_GROUP_ACC_Y];
        let x_in = local_wire_values[Self::WIRE_ADDEND_X];
        let y_in = local_wire_values[Self::WIRE_ADDEND_Y];
        let x3 = right_wire_values[Self::WIRE_GROUP_ACC_X];
        let y3 = right_wire_values[Self::WIRE_GROUP_ACC_Y];
        let scalar_acc_unsigned_old = local_wire_values[Self::WIRE_SCALAR_ACC_UNSIGNED];
        let scalar_acc_unsigned_new = below_wire_values[Self::WIRE_SCALAR_ACC_UNSIGNED];
        let scalar_acc_signed_old = local_wire_values[Self::WIRE_SCALAR_ACC_SIGNED];
        let scalar_acc_signed_new = below_wire_values[Self::WIRE_SCALAR_ACC_SIGNED];
        let scalar_bit_0 = local_wire_values[Self::WIRE_SCALAR_BIT_0];
        let scalar_bit_1 = local_wire_values[Self::WIRE_SCALAR_BIT_1];
        let inverse = local_wire_values[Self::WIRE_INVERSE];

        // Conditionally apply the endo and conditionally negate in order to get (x2, y2), which is
        // the actual point we want to add to the accumulator.
        let zeta_minus_one = builder.sub(zeta, one);
        let x_in_multiplier = builder.mul_add(zeta_minus_one, scalar_bit_1, one);
        let x2 = builder.mul(x_in_multiplier, x_in);
        let y_in_multiplier = builder.mul_sub(scalar_bit_0, two, one);
        let y2 = builder.mul(y_in_multiplier, y_in);

        let x1_minus_x2 = builder.sub(x1, x2);
        let x1_plus_x2 = builder.add(x1, x2);
        let x1_minus_x3 = builder.sub(x1, x3);
        let y1_minus_y2 = builder.sub(y1, y2);

        let lambda = builder.mul(y1_minus_y2, inverse);
        let computed_x3 = builder.mul_sub(lambda, lambda, x1_plus_x2);
        let computed_y3 = builder.mul_sub(lambda, x1_minus_x3, y1);

        let unsigned_limb = builder.mul_add(scalar_bit_1, two, scalar_bit_0);
        let computed_scalar_acc_unsigned_new =
            builder.mul_add(scalar_acc_unsigned_old, four, unsigned_limb);

        // This is based on Algorithm 2 in the Halo paper.
        // TODO: Wrong zeta used here?
        let signed_limb_multiplier = builder.mul_add(zeta_minus_one, scalar_bit_1, one);
        let signed_limb_sign = builder.mul_sub(scalar_bit_0, two, one);
        let signed_limb = builder.mul(signed_limb_sign, signed_limb_multiplier);
        let computed_scalar_acc_signed_new =
            builder.mul_add(scalar_acc_signed_old, two, signed_limb);

        vec![
            builder.sub(computed_x3, x3),
            builder.sub(computed_y3, y3),
            builder.sub(computed_scalar_acc_unsigned_new, scalar_acc_unsigned_new),
            builder.sub(computed_scalar_acc_signed_new, scalar_acc_signed_new),
            assert_binary_recursively(builder, scalar_bit_0),
            assert_binary_recursively(builder, scalar_bit_1),
            assert_inverses_recursively(builder, inverse, x1_minus_x2),
        ]
    }
}

impl<C: HaloCurve, InnerC: HaloCurve<BaseField = C::ScalarField>> WitnessGenerator<C::ScalarField>
for CurveEndoGate<C, InnerC>
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
                input: Self::WIRE_SCALAR_ACC_UNSIGNED,
            }),
            Target::Wire(Wire {
                gate: self.index,
                input: Self::WIRE_SCALAR_ACC_SIGNED,
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
                input: Self::WIRE_SCALAR_BIT_0,
            }),
            Target::Wire(Wire {
                gate: self.index,
                input: Self::WIRE_SCALAR_BIT_1,
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

        let scalar_acc_unsigned_old_target = Wire {
            gate: self.index,
            input: Self::WIRE_GROUP_ACC_Y,
        };
        let scalar_acc_unsigned_new_target = Wire {
            gate: self.index + GRID_WIDTH,
            input: Self::WIRE_GROUP_ACC_Y,
        };
        let scalar_acc_signed_old_target = Wire {
            gate: self.index,
            input: Self::WIRE_GROUP_ACC_Y,
        };
        let scalar_acc_signed_new_target = Wire {
            gate: self.index + GRID_WIDTH,
            input: Self::WIRE_GROUP_ACC_Y,
        };

        let addend_x_target = Wire {
            gate: self.index,
            input: Self::WIRE_ADDEND_X,
        };
        let addend_y_target = Wire {
            gate: self.index,
            input: Self::WIRE_ADDEND_Y,
        };
        let scalar_bit_0_target = Wire {
            gate: self.index,
            input: Self::WIRE_SCALAR_BIT_0,
        };
        let scalar_bit_1_target = Wire {
            gate: self.index,
            input: Self::WIRE_SCALAR_BIT_1,
        };
        let inverse_target = Wire {
            gate: self.index,
            input: Self::WIRE_INVERSE,
        };

        let group_acc_old_x = witness.get_wire(group_acc_old_x_target);
        let group_acc_old_y = witness.get_wire(group_acc_old_y_target);
        let group_acc_old = AffinePoint::<InnerC>::nonzero(group_acc_old_x, group_acc_old_y);

        let scalar_acc_unsigned_old = witness.get_wire(scalar_acc_unsigned_old_target);
        let scalar_acc_signed_old = witness.get_wire(scalar_acc_signed_old_target);

        let scalar_bit_0 = witness.get_wire(scalar_bit_0_target);
        let scalar_bit_1 = witness.get_wire(scalar_bit_1_target);
        debug_assert!(scalar_bit_0.is_zero() || scalar_bit_0.is_one());
        debug_assert!(scalar_bit_1.is_zero() || scalar_bit_1.is_one());

        let p_x = witness.get_wire(addend_x_target);
        let p_y = witness.get_wire(addend_y_target);

        let mut s_i_x = p_x;
        if scalar_bit_0 == InnerC::BaseField::ONE {
            s_i_x = s_i_x * InnerC::ZETA;
        }
        let mut s_i_y = p_y;
        if scalar_bit_1 == InnerC::BaseField::ZERO {
            s_i_y = -s_i_y;
        }
        let s_i = AffinePoint::nonzero(s_i_x, s_i_y);
        let group_acc_new = group_acc_old + s_i;

        let scalar_acc_unsigned_new =
            scalar_acc_unsigned_old.quadruple() + scalar_bit_0 + scalar_bit_1.double();

        // This is based on Algorithm 2 in the Halo paper.
        let mut scalar_acc_signed_limb = if scalar_bit_0 == InnerC::BaseField::ONE {
            InnerC::BaseField::ONE
        } else {
            InnerC::BaseField::NEG_ONE
        };
        if scalar_bit_1 == InnerC::BaseField::ONE {
            scalar_acc_signed_limb = scalar_acc_signed_limb * InnerC::ZETA;
        }
        let scalar_acc_signed_new = scalar_acc_signed_old.double() + scalar_acc_signed_limb;

        // Here's where our abstraction leaks a bit. Although we already have the sum, we need to
        // redo part of the computation in order to populate the purported inverse wire.
        let dx = group_acc_old_x - p_x;
        let _dy = group_acc_old_y - p_y;
        let inverse = dx.multiplicative_inverse().expect("x_1 = x_2");

        let mut result = PartialWitness::new();
        result.set_wire(group_acc_new_x_target, group_acc_new.x);
        result.set_wire(group_acc_new_y_target, group_acc_new.y);
        result.set_wire(scalar_acc_unsigned_new_target, scalar_acc_unsigned_new);
        result.set_wire(scalar_acc_signed_new_target, scalar_acc_signed_new);
        result.set_wire(inverse_target, inverse);

        result
    }
}

#[cfg(test)]
mod tests {
    use crate::{CurveEndoGate, test_gate_low_degree, Tweedledee, Tweedledum};

    test_gate_low_degree!(
        low_degree_CurveEndoGate,
        Tweedledum,
        CurveEndoGate<Tweedledum, Tweedledee>
    );
}