use std::marker::PhantomData;

use crate::gates::{assert_inverses_recursively, Gate};
use crate::{CircuitBuilder, Curve, Field, HaloCurve, PartialWitness, Target, Wire, WitnessGenerator};

/// A curve which performs point doubling.
pub struct CurveDblGate<C: HaloCurve, InnerC: Curve<BaseField = C::ScalarField>> {
    pub index: usize,
    _phantom_oc: PhantomData<C>,
    _phantom_ic: PhantomData<InnerC>,
}

impl<C: HaloCurve, InnerC: Curve<BaseField = C::ScalarField>> CurveDblGate<C, InnerC> {
    pub fn new(index: usize) -> Self {
        CurveDblGate {
            index,
            _phantom_oc: PhantomData,
            _phantom_ic: PhantomData,
        }
    }

    pub const WIRE_X_OLD: usize = 0;
    pub const WIRE_Y_OLD: usize = 1;
    pub const WIRE_X_NEW: usize = 2;
    pub const WIRE_Y_NEW: usize = 3;
    pub const WIRE_INVERSE: usize = 4;
    pub const WIRE_LAMBDA: usize = 5;
}

impl<C: HaloCurve, InnerC: Curve<BaseField = C::ScalarField>> Gate<C> for CurveDblGate<C, InnerC> {
    const NAME: &'static str = "CurveDblGate";
    const DEGREE: usize = 3;
    const NUM_CONSTANTS: usize = 0;

    type Constraints = ();
    const PREFIX: &'static [bool] = &[true, false, true, true, true];

    fn evaluate_unfiltered(
        _local_constant_values: &[InnerC::BaseField],
        local_wire_values: &[InnerC::BaseField],
        _right_wire_values: &[InnerC::BaseField],
        _below_wire_values: &[InnerC::BaseField],
    ) -> Vec<InnerC::BaseField> {
        let x_old = local_wire_values[Self::WIRE_X_OLD];
        let y_old = local_wire_values[Self::WIRE_Y_OLD];
        let x_new = local_wire_values[Self::WIRE_X_NEW];
        let y_new = local_wire_values[Self::WIRE_Y_NEW];
        let inverse = local_wire_values[Self::WIRE_INVERSE];
        let lambda = local_wire_values[Self::WIRE_LAMBDA];

        let computed_lambda_numerator = x_old.square().triple() + InnerC::A;
        let computed_lambda = computed_lambda_numerator * inverse;
        let computed_x_new = lambda.square() - x_old.double();
        let computed_y_new = lambda * (x_old - x_new) - y_old;

        vec![
            // Verify that computed_lambda matches lambda.
            computed_lambda - lambda,
            // Verify that computed_x_new matches x_new.
            computed_x_new - x_new,
            // Verify that computed_y_new matches y_new.
            computed_y_new - y_new,
            // Verify that 2 * y_old times its purported inverse is 1.
            y_old.double() * inverse - InnerC::BaseField::ONE,
        ]
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        _local_constant_values: &[Target<C::ScalarField>],
        local_wire_values: &[Target<C::ScalarField>],
        _right_wire_values: &[Target<C::ScalarField>],
        _below_wire_values: &[Target<C::ScalarField>],
    ) -> Vec<Target<C::ScalarField>> {
        let _one = builder.one_wire();
        let three = builder.constant_wire_u32(3);
        let a = builder.constant_wire(InnerC::A);

        let x_old = local_wire_values[Self::WIRE_X_OLD];
        let y_old = local_wire_values[Self::WIRE_Y_OLD];
        let x_new = local_wire_values[Self::WIRE_X_NEW];
        let y_new = local_wire_values[Self::WIRE_Y_NEW];
        let inverse = local_wire_values[Self::WIRE_INVERSE];
        let lambda = local_wire_values[Self::WIRE_LAMBDA];

        let two_x_old = builder.double(x_old);
        let two_y_old = builder.double(y_old);
        let x_old_squared = builder.square(x_old);
        let three_x_old_squared = builder.mul(three, x_old_squared);
        let computed_lambda_numerator = builder.add(three_x_old_squared, a);
        let computed_lambda = builder.mul(computed_lambda_numerator, inverse);
        let lambda_squared = builder.square(lambda);
        let computed_x_new = builder.sub(lambda_squared, two_x_old);
        let delta_x = builder.sub(x_old, x_new);
        let lambda_times_delta_x = builder.mul(lambda, delta_x);
        let computed_y_new = builder.sub(lambda_times_delta_x, y_old);

        vec![
            // Verify that computed_lambda matches lambda.
            builder.sub(computed_lambda, lambda),
            // Verify that computed_x_new matches x_new.
            builder.sub(computed_x_new, x_new),
            // Verify that computed_y_new matches y_new.
            builder.sub(computed_y_new, y_new),
            // Verify that 2 * y_old times its purported inverse is 1.
            assert_inverses_recursively(builder, two_y_old, inverse),
        ]
    }
}

impl<C: HaloCurve, InnerC: Curve<BaseField = C::ScalarField>> WitnessGenerator<C::ScalarField>
    for CurveDblGate<C, InnerC>
{
    fn dependencies(&self) -> Vec<Target<C::ScalarField>> {
        vec![
            Target::Wire(Wire {
                gate: self.index,
                input: Self::WIRE_X_OLD,
            }),
            Target::Wire(Wire {
                gate: self.index,
                input: Self::WIRE_Y_OLD,
            }),
        ]
    }

    fn generate(
        &self,
        _constants: &[Vec<C::ScalarField>],
        witness: &PartialWitness<InnerC::BaseField>,
    ) -> PartialWitness<InnerC::BaseField> {
        let x_old_target = Wire {
            gate: self.index,
            input: Self::WIRE_X_OLD,
        };
        let y_old_target = Wire {
            gate: self.index,
            input: Self::WIRE_Y_OLD,
        };
        let x_new_target = Wire {
            gate: self.index,
            input: Self::WIRE_X_NEW,
        };
        let y_new_target = Wire {
            gate: self.index,
            input: Self::WIRE_Y_NEW,
        };
        let inverse_target = Wire {
            gate: self.index,
            input: Self::WIRE_INVERSE,
        };
        let lambda_target = Wire {
            gate: self.index,
            input: Self::WIRE_LAMBDA,
        };

        let x_old = witness.get_wire(x_old_target);
        let y_old = witness.get_wire(y_old_target);

        let inverse = y_old.double().multiplicative_inverse().expect("y = 0");
        let lambda = x_old.square().triple() * inverse;
        let x_new = lambda.square() - x_old.double();
        let y_new = lambda * (x_old - x_new) - y_old;

        let mut result = PartialWitness::new();
        result.set_wire(inverse_target, inverse);
        result.set_wire(lambda_target, lambda);
        result.set_wire(x_new_target, x_new);
        result.set_wire(y_new_target, y_new);
        result
    }
}

#[cfg(test)]
mod tests {
    use crate::{test_gate_low_degree, CurveDblGate, Tweedledee, Tweedledum};

    test_gate_low_degree!(
        low_degree_CurveDblGate,
        Tweedledum,
        CurveDblGate<Tweedledum, Tweedledee>
    );
}
