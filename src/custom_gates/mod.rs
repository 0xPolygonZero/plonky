use crate::gates::gate_collection::{GateCollection, GatePrefixes};
use crate::{CircuitBuilder, Curve, Gate, HaloCurve, PartialWitness, Target, WitnessGenerator};
use std::sync::Arc;

pub mod multivariate_polynomial;
pub mod tree_builder;
// mod gate_macro;

pub struct CustomGate<C: HaloCurve> {
    index: usize,
    name: &'static str,
    degree: usize,
    num_constants: usize,
    evaluate: Arc<
        dyn Fn(
                &[C::ScalarField],
                &[C::ScalarField],
                &[C::ScalarField],
                &[C::ScalarField],
            ) -> Vec<C::ScalarField>
            + Sync
            + Send,
    >,
    evaluate_recursively: Arc<
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
    dependencies: Vec<Target<C::ScalarField>>,
    generate: Arc<
        dyn Fn(
                &[Target<C::ScalarField>],
                &[Vec<C::ScalarField>],
                &PartialWitness<C::ScalarField>,
            ) -> PartialWitness<C::ScalarField>
            + Sync
            + Send,
    >,
}

impl<C: HaloCurve> Gate<C> for CustomGate<C> {
    fn name(&self) -> &'static str {
        self.name
    }

    fn degree(&self) -> usize {
        self.degree
    }

    fn num_constants(&self) -> usize {
        self.num_constants
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
        let constants = &local_constant_values[prefix_len..];
        (self.evaluate)(
            constants,
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
        let constants = &local_constant_values[prefix_len..];
        (self.evaluate_recursively)(
            builder,
            constants,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        )
    }
}

impl<C: HaloCurve> WitnessGenerator<C::ScalarField> for CustomGate<C> {
    fn dependencies(&self) -> Vec<Target<<C as Curve>::ScalarField>> {
        self.dependencies.clone()
    }

    fn generate(
        &self,
        prefixes: &GatePrefixes,
        constants: &[Vec<<C as Curve>::ScalarField>],
        witness: &PartialWitness<<C as Curve>::ScalarField>,
    ) -> PartialWitness<<C as Curve>::ScalarField> {
        unimplemented!()
    }
}
