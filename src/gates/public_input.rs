use std::marker::PhantomData;

use crate::gates::gate_collection::{GateCollection, GatePrefixes};
use crate::gates::Gate;
use crate::{CircuitBuilder, HaloCurve, PartialWitness, Target, Wire, WitnessGenerator, NUM_ADVICE_WIRES, NUM_ROUTED_WIRES, NUM_WIRES};

/// A gate for receiving public inputs. These gates will be placed at static indices and the wire
/// polynomials will always be opened at those indices.
///
/// Because our gate arity is 11 but only 6 of the wires are routed, it may seem as though each gate
/// can only receive 6 public inputs. To work around this, we place a BufferGate immediately after
/// each PublicInputGate, and have the PublicInputGate copy its 5 non-routed wires to routed wires
/// of the BufferGate.
pub struct PublicInputGate<C: HaloCurve> {
    pub index: usize,
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> PublicInputGate<C> {
    pub fn new(index: usize) -> Self {
        PublicInputGate {
            index,
            _phantom: PhantomData,
        }
    }
}

impl<C: HaloCurve> Gate<C> for PublicInputGate<C> {
    fn name(&self) -> &'static str {
        "PublicInputGate"
    }
    fn degree(&self) -> usize {
        1
    }
    fn num_constants(&self) -> usize {
        0
    }

    fn evaluate_unfiltered(
        &self,
        gates: &GateCollection<C>,
        _local_constant_values: &[C::ScalarField],
        local_wire_values: &[C::ScalarField],
        right_wire_values: &[C::ScalarField],
        _below_wire_values: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        // This ensures that advice wires' values are copied to the following buffer gate.
        // TODO: Consider enforcing this via copy constraints, in which case there would be nothing to do here.
        (0..NUM_ADVICE_WIRES)
            .map(|i| local_wire_values[NUM_ROUTED_WIRES + i] - right_wire_values[i])
            .collect()
    }

    fn evaluate_unfiltered_recursively(
        &self,
        builder: &mut CircuitBuilder<C>,
        gates: &GateCollection<C>,
        _local_constant_values: &[Target<C::ScalarField>],
        local_wire_values: &[Target<C::ScalarField>],
        right_wire_values: &[Target<C::ScalarField>],
        _below_wire_values: &[Target<C::ScalarField>],
    ) -> Vec<Target<C::ScalarField>> {
        // This ensures that advice wires' values are copied to the following buffer gate.
        // TODO: Consider enforcing this via copy constraints, in which case there would be nothing to do here.
        (0..NUM_ADVICE_WIRES)
            .map(|i| {
                builder.sub(
                    local_wire_values[NUM_ROUTED_WIRES + i],
                    right_wire_values[i],
                )
            })
            .collect()
    }
}

impl<C: HaloCurve> WitnessGenerator<C::ScalarField> for PublicInputGate<C> {
    fn dependencies(&self) -> Vec<Target<C::ScalarField>> {
        Vec::new()
    }

    fn generate(
        &self,
        prefixes: &GatePrefixes,
        _constants: &[Vec<C::ScalarField>],
        witness: &PartialWitness<C::ScalarField>,
    ) -> PartialWitness<C::ScalarField> {
        let targets: Vec<Target<C::ScalarField>> = (0..NUM_WIRES)
            .map(|i| {
                Target::Wire(Wire {
                    gate: self.index,
                    input: i,
                })
            })
            .collect();

        let mut result = PartialWitness::new();
        for i_advice in 0..NUM_ADVICE_WIRES {
            let i_wire = NUM_ROUTED_WIRES + i_advice;
            if witness.contains_target(targets[i_wire]) {
                let value = witness.get_target(targets[i_wire]);
                result.set_wire(
                    Wire {
                        gate: self.index + 1,
                        input: i_advice,
                    },
                    value,
                );
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use crate::{test_gate_low_degree, PublicInputGate, Tweedledee, Tweedledum};

    test_gate_low_degree!(
        low_degree_PublicInputGate,
        Tweedledum,
        Tweedledee,
        PublicInputGate<Tweedledum>
    );
}
