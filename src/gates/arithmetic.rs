use std::marker::PhantomData;

use crate::gates::gate_collection::{GateCollection, GatePrefixes};
use crate::gates::Gate;
use crate::{CircuitBuilder, HaloCurve, PartialWitness, Target, Wire, WitnessGenerator};

/// A gate which can be configured to perform various arithmetic. In particular, it computes
///
/// ```text
/// output := const_0 * multiplicand_0 * multiplicand_1 + const_1 * addend
/// ```
pub struct ArithmeticGate<C: HaloCurve> {
    pub index: usize,
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> ArithmeticGate<C> {
    pub fn new(index: usize) -> Self {
        ArithmeticGate {
            index,
            _phantom: PhantomData,
        }
    }

    pub const WIRE_MULTIPLICAND_0: usize = 0;
    pub const WIRE_MULTIPLICAND_1: usize = 1;
    pub const WIRE_ADDEND: usize = 2;
    pub const WIRE_OUTPUT: usize = 3;
}

impl<C: HaloCurve> Gate<C> for ArithmeticGate<C> {
    fn name(&self) -> &'static str {
        "ArithmeticGate"
    }
    fn degree(&self) -> usize {
        3
    }
    fn num_constants(&self) -> usize {
        2
    }

    fn evaluate_unfiltered(
        &self,
        gates: &GateCollection<C>,
        local_constant_values: &[C::ScalarField],
        local_wire_values: &[C::ScalarField],
        _right_wire_values: &[C::ScalarField],
        _below_wire_values: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        let prefix_len = gates.prefix(self).len();
        let const_0 = local_constant_values[prefix_len];
        let const_1 = local_constant_values[prefix_len + 1];
        let multiplicand_0 = local_wire_values[Self::WIRE_MULTIPLICAND_0];
        let multiplicand_1 = local_wire_values[Self::WIRE_MULTIPLICAND_1];
        let addend = local_wire_values[Self::WIRE_ADDEND];
        let output = local_wire_values[Self::WIRE_OUTPUT];
        let computed_output = const_0 * multiplicand_0 * multiplicand_1 + const_1 * addend;
        vec![computed_output - output]
    }

    fn evaluate_unfiltered_recursively(
        &self,
        builder: &mut CircuitBuilder<C>,
        gates: &GateCollection<C>,
        local_constant_values: &[Target<C::ScalarField>],
        local_wire_values: &[Target<C::ScalarField>],
        _right_wire_values: &[Target<C::ScalarField>],
        _below_wire_values: &[Target<C::ScalarField>],
    ) -> Vec<Target<C::ScalarField>> {
        let prefix_len = gates.prefix(self).len();
        let const_0 = local_constant_values[prefix_len];
        let const_1 = local_constant_values[prefix_len + 1];
        let multiplicand_0 = local_wire_values[Self::WIRE_MULTIPLICAND_0];
        let multiplicand_1 = local_wire_values[Self::WIRE_MULTIPLICAND_1];
        let addend = local_wire_values[Self::WIRE_ADDEND];
        let output = local_wire_values[Self::WIRE_OUTPUT];

        let product_term = builder.mul_many(&[const_0, multiplicand_0, multiplicand_1]);
        let addend_term = builder.mul(const_1, addend);
        let computed_output = builder.add_many(&[product_term, addend_term]);
        vec![builder.sub(computed_output, output)]
    }
}

impl<C: HaloCurve> WitnessGenerator<C::ScalarField> for ArithmeticGate<C> {
    fn dependencies(&self) -> Vec<Target<C::ScalarField>> {
        vec![
            Target::Wire(Wire {
                gate: self.index,
                input: Self::WIRE_MULTIPLICAND_0,
            }),
            Target::Wire(Wire {
                gate: self.index,
                input: Self::WIRE_MULTIPLICAND_1,
            }),
            Target::Wire(Wire {
                gate: self.index,
                input: Self::WIRE_ADDEND,
            }),
        ]
    }

    fn generate(
        &self,
        prefixes: &GatePrefixes,
        constants: &[Vec<C::ScalarField>],
        witness: &PartialWitness<C::ScalarField>,
    ) -> PartialWitness<C::ScalarField> {
        let multiplicand_0_target = Wire {
            gate: self.index,
            input: Self::WIRE_MULTIPLICAND_0,
        };
        let multiplicand_1_target = Wire {
            gate: self.index,
            input: Self::WIRE_MULTIPLICAND_1,
        };
        let addend_target = Wire {
            gate: self.index,
            input: Self::WIRE_ADDEND,
        };
        let output_target = Wire {
            gate: self.index,
            input: Self::WIRE_OUTPUT,
        };

        let prefix_len = prefixes
            .get(self.name())
            .expect(&format!("Gate {} not found.", self.name()))
            .len();
        let const_0 = constants[self.index][prefix_len];
        let const_1 = constants[self.index][prefix_len + 1];

        let multiplicand_0 = witness.get_wire(multiplicand_0_target);
        let multiplicand_1 = witness.get_wire(multiplicand_1_target);
        let addend = witness.get_wire(addend_target);

        let output = const_0 * multiplicand_0 * multiplicand_1 + const_1 * addend;

        let mut result = PartialWitness::new();
        result.set_wire(output_target, output);
        result
    }
}

#[cfg(test)]
mod tests {
    use crate::{test_gate_low_degree, ArithmeticGate, Tweedledee, Tweedledum};

    test_gate_low_degree!(
        low_degree_ArithmeticGate,
        Tweedledum,
        Tweedledee,
        ArithmeticGate<Tweedledum>
    );
}
