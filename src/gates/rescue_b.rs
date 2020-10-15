use std::marker::PhantomData;

use crate::gates::gate_collection::{GateCollection, GatePrefixes};
use crate::gates::Gate;
use crate::{mds_matrix, CircuitBuilder, Field, HaloCurve, PartialWitness, Target, Wire, WitnessGenerator, RESCUE_SPONGE_WIDTH};

/// The second step of Rescue, i.e. the one with the `x^5` layer.
pub struct RescueStepBGate<C: HaloCurve> {
    pub index: usize,
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> RescueStepBGate<C> {
    pub fn new(index: usize) -> Self {
        RescueStepBGate {
            index,
            _phantom: PhantomData,
        }
    }

    /// Returns the index of the `i`th accumulator wire.
    pub fn wire_acc(i: usize) -> usize {
        i
    }
}

impl<C: HaloCurve> Gate<C> for RescueStepBGate<C> {
    fn name(&self) -> &'static str {
        "RescueStepBGate"
    }
    fn degree(&self) -> usize {
        5
    }
    fn num_constants(&self) -> usize {
        4
    }

    fn evaluate_unfiltered(
        &self,
        gates: &GateCollection<C>,
        local_constant_values: &[C::ScalarField],
        local_wire_values: &[C::ScalarField],
        right_wire_values: &[C::ScalarField],
        _below_wire_values: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        let ins: Vec<C::ScalarField> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| local_wire_values[Self::wire_acc(i)])
            .collect();

        let exps: Vec<C::ScalarField> = ins.into_iter().map(|n| n.exp_usize(5)).collect();

        let outs: Vec<C::ScalarField> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| right_wire_values[Self::wire_acc(i)])
            .collect();

        let mds = mds_matrix::<C::ScalarField>(RESCUE_SPONGE_WIDTH);

        let prefix_len = gates.prefix(self).len();
        let mut constraints = Vec::new();
        for i in 0..RESCUE_SPONGE_WIDTH {
            let mut computed_out_i = local_constant_values[prefix_len + i];
            for j in 0..RESCUE_SPONGE_WIDTH {
                computed_out_i = computed_out_i + mds.get(i, j) * exps[j];
            }
            constraints.push(computed_out_i - outs[i]);
        }
        constraints
    }

    fn evaluate_unfiltered_recursively(
        &self,
        builder: &mut CircuitBuilder<C>,
        gates: &GateCollection<C>,
        local_constant_values: &[Target<C::ScalarField>],
        local_wire_values: &[Target<C::ScalarField>],
        right_wire_values: &[Target<C::ScalarField>],
        _below_wire_values: &[Target<C::ScalarField>],
    ) -> Vec<Target<C::ScalarField>> {
        let ins: Vec<Target<C::ScalarField>> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| local_wire_values[Self::wire_acc(i)])
            .collect();

        let exps: Vec<Target<C::ScalarField>> = ins
            .into_iter()
            .map(|n| builder.exp_constant_usize(n, 5))
            .collect();

        let outs: Vec<Target<C::ScalarField>> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| right_wire_values[Self::wire_acc(i)])
            .collect();

        let mds = mds_matrix::<C::ScalarField>(RESCUE_SPONGE_WIDTH);

        let prefix_len = gates.prefix(self).len();
        let mut constraints = Vec::new();
        for i in 0..RESCUE_SPONGE_WIDTH {
            let mut computed_out_i = local_constant_values[prefix_len + i];
            for j in 0..RESCUE_SPONGE_WIDTH {
                let mds_entry = builder.constant_wire(mds.get(i, j));
                computed_out_i = builder.mul_add(mds_entry, exps[i], computed_out_i);
            }
            constraints.push(builder.sub(computed_out_i, outs[i]));
        }
        constraints
    }
}

impl<C: HaloCurve> WitnessGenerator<C::ScalarField> for RescueStepBGate<C> {
    fn dependencies(&self) -> Vec<Target<C::ScalarField>> {
        (0..RESCUE_SPONGE_WIDTH)
            .map(|i| {
                Target::Wire(Wire {
                    gate: self.index,
                    input: Self::wire_acc(i),
                })
            })
            .collect()
    }

    fn generate(
        &self,
        prefixes: &GatePrefixes,
        constants: &[Vec<C::ScalarField>],
        witness: &PartialWitness<C::ScalarField>,
    ) -> PartialWitness<C::ScalarField> {
        let constants = &constants[self.index];

        let ins: Vec<C::ScalarField> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| {
                witness.get_wire(Wire {
                    gate: self.index,
                    input: Self::wire_acc(i),
                })
            })
            .collect();

        let exps: Vec<C::ScalarField> = ins.iter().map(|n| n.exp_usize(5)).collect();

        let mds = mds_matrix::<C::ScalarField>(RESCUE_SPONGE_WIDTH);

        let prefix_len = prefixes
            .get(self.name())
            .unwrap_or_else(|| panic!("Gate {} not found.", self.name()))
            .len();
        let mut result = PartialWitness::new();
        for i in 0..RESCUE_SPONGE_WIDTH {
            let mut out_i = constants[prefix_len + i];
            for j in 0..RESCUE_SPONGE_WIDTH {
                out_i = out_i + mds.get(i, j) * exps[j];
            }
            let wire_out_i = Wire {
                gate: self.index + 1,
                input: Self::wire_acc(i),
            };
            result.set_wire(wire_out_i, out_i);
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use crate::{test_gate_low_degree, RescueStepBGate, Tweedledee, Tweedledum};

    test_gate_low_degree!(
        low_degree_RescueStepBGate,
        Tweedledum,
        Tweedledee,
        RescueStepBGate<Tweedledum>
    );
}
