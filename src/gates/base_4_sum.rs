use std::marker::PhantomData;

use crate::gates::Gate;
use crate::{CircuitBuilder, Curve, Field, HaloCurve, PartialWitness, Target, WitnessGenerator, NUM_ROUTED_WIRES, NUM_WIRES};

/// A gate for accumulating base-4 limbs.
pub struct Base4SumGate<C: Curve> {
    pub index: usize,
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> Base4SumGate<C> {
    pub fn new(index: usize) -> Self {
        Base4SumGate {
            index,
            _phantom: PhantomData,
        }
    }

    pub const WIRE_ACC_OLD: usize = 0;
    pub const WIRE_ACC_NEW: usize = 1;
    pub const NUM_LIMBS: usize = NUM_WIRES - 2;
    pub const NUM_ROUTED_LIMBS: usize = NUM_ROUTED_WIRES - 2;

    /// Returns the index of the `i`th limb wire.
    pub fn wire_limb(i: usize) -> usize {
        2 + i
    }
}

impl<C: HaloCurve> Gate<C> for Base4SumGate<C> {
    const NAME: &'static str = "Base4SumGate";

    const PREFIX: &'static [bool] = &[true, false, false, false];

    fn evaluate_unfiltered(
        _local_constant_values: &[C::ScalarField],
        local_wire_values: &[C::ScalarField],
        _right_wire_values: &[C::ScalarField],
        _below_wire_values: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        let acc_old = local_wire_values[Self::WIRE_ACC_OLD];
        let acc_new = local_wire_values[Self::WIRE_ACC_NEW];
        let limbs: Vec<C::ScalarField> = (0..Self::NUM_LIMBS)
            .map(|i| local_wire_values[Self::wire_limb(i)])
            .collect();

        let mut computed_acc_new = acc_old;
        for &limb in &limbs {
            computed_acc_new = computed_acc_new.quadruple() + limb;
        }

        let mut constraints = vec![computed_acc_new - acc_new];
        for limb in limbs {
            let mut product = C::ScalarField::ONE;
            for j in 0..4 {
                product = product * (limb - C::ScalarField::from_canonical_usize(j));
            }
            constraints.push(product);
        }
        constraints
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        _local_constant_values: &[Target<C::ScalarField>],
        local_wire_values: &[Target<C::ScalarField>],
        _right_wire_values: &[Target<C::ScalarField>],
        _below_wire_values: &[Target<C::ScalarField>],
    ) -> Vec<Target<C::ScalarField>> {
        let four = builder.constant_wire_u32(4);

        let acc_old = local_wire_values[Self::WIRE_ACC_OLD];
        let acc_new = local_wire_values[Self::WIRE_ACC_NEW];
        let limbs: Vec<Target<C::ScalarField>> = (0..Self::NUM_LIMBS)
            .map(|i| local_wire_values[Self::wire_limb(i)])
            .collect();

        let mut computed_acc_new = acc_old;
        for &limb in &limbs {
            let shifted_acc_new = builder.mul(computed_acc_new, four);
            computed_acc_new = builder.add(shifted_acc_new, limb);
        }

        let mut constraints = vec![builder.sub(computed_acc_new, acc_new)];
        for limb in limbs {
            let mut product = builder.one_wire();
            for j in 0..4 {
                let j_target = builder.constant_wire_u32(j);
                let limb_minus_j = builder.sub(limb, j_target);
                product = builder.mul(product, limb_minus_j);
            }
            constraints.push(product);
        }
        constraints
    }
}

impl<C: HaloCurve> WitnessGenerator<C::ScalarField> for Base4SumGate<C> {
    fn dependencies(&self) -> Vec<Target<C::ScalarField>> {
        Vec::new()
    }

    fn generate(
        &self,
        _constants: &[Vec<C::ScalarField>],
        _witness: &PartialWitness<C::ScalarField>,
    ) -> PartialWitness<C::ScalarField> {
        // For base 4 decompositions, we don't do any witness generation on a per-gate level.
        // Instead, we have a single generator which generates values for an entire decomposition.
        PartialWitness::new()
    }
}

#[cfg(test)]
mod tests {
    use crate::{test_gate_low_degree, Base4SumGate, Tweedledum};

    test_gate_low_degree!(
        low_degree_Base4SumGate,
        Tweedledum,
        Base4SumGate<Tweedledum>
    );
}
