use std::marker::PhantomData;

use crate::{CircuitBuilder, Field, HaloCurve, mds_matrix, PartialWitness, RESCUE_SPONGE_WIDTH, Target, Wire, WitnessGenerator};
use crate::gates::Gate;

/// The first step of Rescue, i.e. the one with the `x^(1/5)` layer.
pub struct RescueStepAGate<C: HaloCurve> {
    pub index: usize,
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> RescueStepAGate<C> {
    pub fn new(index: usize) -> Self {
        RescueStepAGate {
            index,
            _phantom: PhantomData,
        }
    }

    /// Returns the index of the `i`th accumulator wire.
    pub fn wire_acc(i: usize) -> usize {
        return i;
    }

    /// Returns the index of the `i`th root wire.
    pub fn wire_root(i: usize) -> usize {
        return RESCUE_SPONGE_WIDTH + i;
    }
}

impl<C: HaloCurve> Gate<C> for RescueStepAGate<C> {
    const NAME: &'static str = "RescueStepAGate";

    const PREFIX: &'static [bool] = &[false, false];

    fn evaluate_unfiltered(
        local_constant_values: &[C::ScalarField],
        local_wire_values: &[C::ScalarField],
        right_wire_values: &[C::ScalarField],
        _below_wire_values: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        let ins: Vec<C::ScalarField> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| local_wire_values[Self::wire_acc(i)])
            .collect();
        let outs: Vec<C::ScalarField> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| right_wire_values[Self::wire_acc(i)])
            .collect();
        let roots: Vec<C::ScalarField> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| local_wire_values[Self::wire_root(i)])
            .collect();

        let mds = mds_matrix::<C::ScalarField>(RESCUE_SPONGE_WIDTH);

        let mut constraints = Vec::new();
        for i in 0..RESCUE_SPONGE_WIDTH {
            constraints.push(roots[i].exp_usize(5) - ins[i]);

            let mut computed_out_i = local_constant_values[Self::PREFIX.len() + i];
            for j in 0..RESCUE_SPONGE_WIDTH {
                computed_out_i = computed_out_i + mds.get(i, j) * roots[j];
            }
            constraints.push(computed_out_i - outs[i]);
        }
        constraints
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        _below_wire_values: &[Target],
    ) -> Vec<Target> {
        let ins: Vec<Target> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| local_wire_values[Self::wire_acc(i)])
            .collect();

        let outs: Vec<Target> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| right_wire_values[Self::wire_acc(i)])
            .collect();

        let roots: Vec<Target> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| local_wire_values[Self::wire_root(i)])
            .collect();

        let mds = mds_matrix::<C::ScalarField>(RESCUE_SPONGE_WIDTH);

        let mut constraints = Vec::new();
        for i in 0..RESCUE_SPONGE_WIDTH {
            let computed_in_i = builder.exp_constant_usize(roots[i], 5);
            constraints.push(builder.sub(computed_in_i, ins[i]));

            let mut computed_out_i = local_constant_values[Self::PREFIX.len() + i];
            for j in 0..RESCUE_SPONGE_WIDTH {
                let mds_entry = builder.constant_wire(mds.get(i, j));
                computed_out_i = builder.mul_add(mds_entry, roots[j], computed_out_i);
            }
            constraints.push(builder.sub(computed_out_i, outs[i]));
        }
        constraints
    }
}

impl<C: HaloCurve> WitnessGenerator<C::ScalarField> for RescueStepAGate<C> {
    fn dependencies(&self) -> Vec<Target> {
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
        constants: &Vec<Vec<C::ScalarField>>,
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

        let roots: Vec<C::ScalarField> = ins.iter().map(|n| n.kth_root_u32(5)).collect();

        let mds = mds_matrix::<C::ScalarField>(RESCUE_SPONGE_WIDTH);

        let mut result = PartialWitness::new();
        for i in 0..RESCUE_SPONGE_WIDTH {
            let wire_root_i = Wire {
                gate: self.index,
                input: Self::wire_root(i),
            };
            result.set_wire(wire_root_i, roots[i]);

            let mut out_i = constants[Self::PREFIX.len() + i];
            for j in 0..RESCUE_SPONGE_WIDTH {
                out_i = out_i + mds.get(i, j) * roots[j];
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
    use crate::{ArithmeticGate, test_gate_low_degree, Tweedledum};

    test_gate_low_degree!(
        low_degree_ArithmeticGate,
        Tweedledum,
        ArithmeticGate<Tweedledum>
    );
}