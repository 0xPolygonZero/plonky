use std::marker::PhantomData;

use crate::gates::Gate;
use crate::{CircuitBuilder, HaloCurve, PartialWitness, Target, Wire, WitnessGenerator};

/// A gate which takes a single constant parameter and outputs that value.
pub struct ConstantGate<C: HaloCurve> {
    pub index: usize,
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> ConstantGate<C> {
    pub const WIRE_OUTPUT: usize = 0;

    pub fn new(index: usize) -> Self {
        ConstantGate {
            index,
            _phantom: PhantomData,
        }
    }
}

impl<C: HaloCurve> Gate<C> for ConstantGate<C> {
    const NAME: &'static str = "ConstantGate";

    const PREFIX: &'static [bool] = &[true, false, true, true, false];

    fn evaluate_unfiltered(
        local_constant_values: &[C::ScalarField],
        local_wire_values: &[C::ScalarField],
        _right_wire_values: &[C::ScalarField],
        _below_wire_values: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        let c = local_constant_values[Self::PREFIX.len()];
        let out = local_wire_values[Self::WIRE_OUTPUT];
        vec![c - out]
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        local_constant_values: &[Target<C::ScalarField>],
        local_wire_values: &[Target<C::ScalarField>],
        _right_wire_values: &[Target<C::ScalarField>],
        _below_wire_values: &[Target<C::ScalarField>],
    ) -> Vec<Target<C::ScalarField>> {
        let c = local_constant_values[Self::PREFIX.len()];
        let out = local_wire_values[Self::WIRE_OUTPUT];
        vec![builder.sub(c, out)]
    }
}

impl<C: HaloCurve> WitnessGenerator<C::ScalarField> for ConstantGate<C> {
    fn dependencies(&self) -> Vec<Target<C::ScalarField>> {
        Vec::new()
    }

    fn generate(
        &self,
        constants: &Vec<Vec<C::ScalarField>>,
        _witness: &PartialWitness<C::ScalarField>,
    ) -> PartialWitness<C::ScalarField> {
        let constants = &constants[self.index];
        let c = constants[Self::PREFIX.len()];
        let mut result = PartialWitness::new();
        result.set_wire(
            Wire {
                gate: self.index,
                input: Self::WIRE_OUTPUT,
            },
            c,
        );
        result
    }
}

#[cfg(test)]
mod tests {
    use crate::{test_gate_low_degree, ConstantGate, Tweedledum};

    test_gate_low_degree!(
        low_degree_ConstantGate,
        Tweedledum,
        ConstantGate<Tweedledum>
    );
}
