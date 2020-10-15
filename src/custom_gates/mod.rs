use crate::gates::gate_collection::{GateCollection, GatePrefixes};
use crate::{CircuitBuilder, Curve, Field, Gate, HaloCurve, PartialWitness, Target, WitnessGenerator, NUM_CONSTANTS};
use std::sync::Arc;

pub mod multivariate_polynomial;
pub mod tree_builder;

#[allow(clippy::type_complexity)]
#[derive(Clone)]
pub struct CustomGateCore<C: HaloCurve> {
    pub name: &'static str,
    pub degree: usize,
    pub num_constants: usize,
    pub evaluate: Arc<
        dyn Fn(
                &[C::ScalarField],
                &[C::ScalarField],
                &[C::ScalarField],
                &[C::ScalarField],
            ) -> Vec<C::ScalarField>
            + Sync
            + Send,
    >,
    pub evaluate_recursively: Arc<
        dyn Fn(
                &mut CircuitBuilder<C>,
                &[Target<C::ScalarField>],
                &[Target<C::ScalarField>],
                &[Target<C::ScalarField>],
                &[Target<C::ScalarField>],
            ) -> Vec<Target<C::ScalarField>>
            + Sync
            + Send,
    >,
    pub dependencies: Arc<dyn Fn(usize) -> Vec<Target<C::ScalarField>> + Send + Sync>,
    pub generate: Arc<
        dyn Fn(
                usize,
                usize,
                &[Vec<C::ScalarField>],
                &PartialWitness<C::ScalarField>,
            ) -> PartialWitness<C::ScalarField>
            + Sync
            + Send,
    >,
}

impl<C: HaloCurve> CustomGateCore<C> {
    pub fn at(self, index: usize) -> CustomGate<C> {
        CustomGate { index, core: self }
    }
}

pub struct CustomGate<C: HaloCurve> {
    pub index: usize,
    pub core: CustomGateCore<C>,
}

impl<C: HaloCurve> Gate<C> for CustomGate<C> {
    fn name(&self) -> &'static str {
        self.core.name
    }

    fn degree(&self) -> usize {
        self.core.degree
    }

    fn num_constants(&self) -> usize {
        self.core.num_constants
    }

    fn evaluate_unfiltered(
        &self,
        gates: &GateCollection<C>,
        local_constant_values: &[<C as Curve>::ScalarField],
        local_wire_values: &[<C as Curve>::ScalarField],
        right_wire_values: &[<C as Curve>::ScalarField],
        below_wire_values: &[<C as Curve>::ScalarField],
    ) -> Vec<<C as Curve>::ScalarField> {
        let prefix_len = gates.prefix(self).len();
        let mut constants = [C::ScalarField::ZERO; NUM_CONSTANTS];
        constants[..(NUM_CONSTANTS - prefix_len)]
            .copy_from_slice(&local_constant_values[prefix_len..]);
        (self.core.evaluate)(
            &constants,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        )
    }

    fn evaluate_unfiltered_recursively(
        &self,
        builder: &mut CircuitBuilder<C>,
        gates: &GateCollection<C>,
        local_constant_values: &[Target<<C as Curve>::ScalarField>],
        local_wire_values: &[Target<<C as Curve>::ScalarField>],
        right_wire_values: &[Target<<C as Curve>::ScalarField>],
        below_wire_values: &[Target<<C as Curve>::ScalarField>],
    ) -> Vec<Target<<C as Curve>::ScalarField>> {
        let prefix_len = gates.prefix(self).len();
        let zero = builder.zero_wire();
        let mut constants = [zero; NUM_CONSTANTS];
        constants[..(NUM_CONSTANTS - prefix_len)]
            .copy_from_slice(&local_constant_values[prefix_len..]);
        (self.core.evaluate_recursively)(
            builder,
            &constants,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        )
    }
}

impl<C: HaloCurve> WitnessGenerator<C::ScalarField> for CustomGate<C> {
    fn dependencies(&self) -> Vec<Target<<C as Curve>::ScalarField>> {
        (self.core.dependencies)(self.index)
    }

    fn generate(
        &self,
        prefixes: &GatePrefixes,
        constants: &[Vec<<C as Curve>::ScalarField>],
        witness: &PartialWitness<<C as Curve>::ScalarField>,
    ) -> PartialWitness<<C as Curve>::ScalarField> {
        let prefix_len = prefixes
            .get(self.name())
            .unwrap_or_else(|| panic!("Gate {} not found.", self.name()))
            .len();
        (self.core.generate)(self.index, prefix_len, constants, witness)
    }
}

#[cfg(test)]
mod tests {
    use crate::custom_gates::CustomGateCore;
    use crate::{PartialWitness, Target, Tweedledee, Wire, NUM_WIRES};
    use std::sync::Arc;

    #[test]
    fn test_custom_gate_construction() {
        let _arithm_gate = CustomGateCore::<Tweedledee> {
            name: "CustomArithmeticGate",
            degree: 3,
            num_constants: 2,
            evaluate: Arc::new(move |constants, local, _right, _below| {
                let c0 = constants[0];
                let c1 = constants[1];

                let x = local[0];
                let y = local[1];
                let z = local[2];

                let output = local[3];

                vec![c0 * x * y + c1 * z - output]
            }),
            evaluate_recursively: Arc::new(move |builder, constants, local, _right, _below| {
                let c0 = constants[0];
                let c1 = constants[1];

                let x = local[0];
                let y = local[1];
                let z = local[2];
                let output = local[3];

                let product_term = builder.mul_many(&[c0, x, y]);
                let add_term = builder.mul(c1, z);
                let computed_output = builder.add(product_term, add_term);
                vec![builder.sub(output, computed_output)]
            }),
            dependencies: Arc::new(|i| {
                vec![
                    Target::Wire(Wire { gate: i, input: 0 }),
                    Target::Wire(Wire { gate: i, input: 1 }),
                    Target::Wire(Wire { gate: i, input: 2 }),
                ]
            }),
            generate: Arc::new(move |index, _prefix_len, constants, witness| {
                let targets = (0..NUM_WIRES)
                    .map(|i| {
                        Target::Wire(Wire {
                            gate: index,
                            input: i,
                        })
                    })
                    .collect::<Vec<_>>();
                let x = witness.get_target(targets[0]);
                let y = witness.get_target(targets[1]);
                let z = witness.get_target(targets[2]);

                let c0 = constants[index][0];
                let c1 = constants[index][1];

                let output = c0 * x * y + c1 * z;

                let mut result = PartialWitness::new();
                result.set_target(targets[3], output);
                result
            }),
        };
    }
}
