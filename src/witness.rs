use crate::gates::gate_collection::{GateCollection, GatePrefixes};
use crate::util::transpose;
use crate::{biguint_to_field, biguint_to_limbs, field_to_biguint, AffinePoint, AffinePointTarget, BigIntTarget, Curve, Field, ForeignFieldTarget, OrderingTarget, Target, Wire, LIMB_BITS, NUM_ADVICE_WIRES, NUM_ROUTED_WIRES, NUM_WIRES};
use num::{BigUint, Zero};
use std::{cmp::Ordering, collections::HashMap};

#[derive(Debug)]
pub struct PartialWitness<F: Field> {
    wire_values: HashMap<Target<F>, F>,
}

impl<F: Field> Default for PartialWitness<F> {
    fn default() -> Self {
        PartialWitness::new()
    }
}

impl<F: Field> PartialWitness<F> {
    pub fn new() -> Self {
        PartialWitness {
            wire_values: HashMap::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.wire_values.is_empty()
    }

    pub fn contains_target(&self, target: Target<F>) -> bool {
        self.wire_values.contains_key(&target)
    }

    pub fn contains_wire(&self, wire: Wire) -> bool {
        self.contains_target(Target::Wire(wire))
    }

    pub fn contains_all_targets(&self, targets: &[Target<F>]) -> bool {
        targets.iter().all(|&t| self.contains_target(t))
    }

    pub fn all_populated_targets(&self) -> Vec<Target<F>> {
        self.wire_values.keys().cloned().collect()
    }

    pub fn get_target(&self, target: Target<F>) -> F {
        self.wire_values[&target]
    }

    pub fn get_targets(&self, targets: &[Target<F>]) -> Vec<F> {
        targets.iter().map(|&t| self.get_target(t)).collect()
    }

    pub fn get_point_target<InnerC: Curve<BaseField = F>>(
        &self,
        target: AffinePointTarget<InnerC>,
    ) -> AffinePoint<InnerC> {
        let x = self.get_target(target.x);
        let y = self.get_target(target.y);
        AffinePoint::nonzero(x, y)
    }

    pub fn get_ordering_target(&self, target: OrderingTarget<F>) -> Ordering {
        let OrderingTarget { lt, eq, gt } = target;
        let lt_eq_gt = &self.get_targets(&[lt, eq, gt]);

        if lt_eq_gt == &[F::ONE, F::ZERO, F::ZERO] {
            return Ordering::Less;
        }
        if lt_eq_gt == &[F::ZERO, F::ONE, F::ZERO] {
            return Ordering::Equal;
        }
        if lt_eq_gt == &[F::ZERO, F::ZERO, F::ONE] {
            return Ordering::Greater;
        }

        panic!("Invalid ordering values")
    }

    pub fn set_ordering_target(&mut self, target: OrderingTarget<F>, value: Ordering) {
        let OrderingTarget { lt, eq, gt } = target;
        let values = match value {
            Ordering::Less => [F::ONE, F::ZERO, F::ZERO],
            Ordering::Equal => [F::ZERO, F::ONE, F::ZERO],
            Ordering::Greater => [F::ZERO, F::ZERO, F::ONE],
        };
        self.set_targets(&[lt, eq, gt], &values);
    }

    pub fn get_bigint_target(&self, target: &BigIntTarget<F>) -> BigUint {
        let mut result = BigUint::zero();
        for (i, &limb) in target.limbs.iter().enumerate() {
            let limb_value = field_to_biguint(self.get_target(limb));
            result += limb_value << (i * LIMB_BITS);
        }
        result
    }

    pub fn set_bigint_target(&mut self, target: &BigIntTarget<F>, value: &BigUint) {
        let mut value_limbs = biguint_to_limbs(value);

        debug_assert!(
            value_limbs.len() <= target.limbs.len(),
            "Not enough limbs to fit the given value"
        );

        while value_limbs.len() < target.limbs.len() {
            value_limbs.push(F::ZERO);
        }

        self.set_targets(&target.limbs, &value_limbs);
    }

    pub fn get_foreign_field_target<FF: Field>(&self, target: &ForeignFieldTarget<F, FF>) -> FF {
        biguint_to_field(self.get_bigint_target(&target.value))
    }

    pub fn set_foreign_field_target<FF: Field>(
        &mut self,
        target: &ForeignFieldTarget<F, FF>,
        value: FF,
    ) {
        self.set_bigint_target(&target.value, &field_to_biguint(value))
    }

    pub fn get_wire(&self, wire: Wire) -> F {
        self.get_target(Target::Wire(wire))
    }

    pub fn set_target(&mut self, target: Target<F>, value: F) {
        let opt_old_value = self.wire_values.insert(target, value);
        if let Some(old_value) = opt_old_value {
            debug_assert_eq!(
                old_value, value,
                "Target {:?} was set twice with different values",
                target
            );
        }
    }

    pub fn set_targets(&mut self, targets: &[Target<F>], values: &[F]) {
        debug_assert_eq!(targets.len(), values.len());
        targets
            .iter()
            .zip(values.iter())
            .for_each(|(&target, &value)| self.set_target(target, value))
    }

    pub fn set_point_target<InnerC: Curve<BaseField = F>>(
        &mut self,
        point_target: AffinePointTarget<InnerC>,
        point: AffinePoint<InnerC>,
    ) {
        self.set_target(point_target.x, point.x);
        self.set_target(point_target.y, point.y);
    }

    pub fn set_point_targets<InnerC: Curve<BaseField = F>>(
        &mut self,
        point_targets: &[AffinePointTarget<InnerC>],
        points: &[AffinePoint<InnerC>],
    ) {
        debug_assert_eq!(point_targets.len(), points.len());
        point_targets
            .iter()
            .zip(points.iter())
            .for_each(|(&point_target, &point)| self.set_point_target(point_target, point))
    }

    pub fn set_wire(&mut self, wire: Wire, value: F) {
        self.set_target(Target::Wire(wire), value);
    }

    pub fn extend(&mut self, other: PartialWitness<F>) {
        for (target, value) in other.wire_values {
            self.set_target(target, value);
        }
    }

    /// Replace all `PublicInput`-type targets by their corresponding `Wire`-type targets
    /// in the partial witness.
    pub(crate) fn replace_public_inputs(&mut self, offset: usize) {
        let new_pis = self
            .wire_values
            .iter()
            .filter_map(|(t, v)| {
                if let Target::PublicInput(pi) = t {
                    Some((Target::Wire(pi.original_wire(offset)), *v))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        self.wire_values
            .retain(|t, _| !matches!(t, Target::PublicInput(_)));
        self.wire_values.extend(new_pis);
    }

    /// Looks through the keys and looks for targets in `BufferGate`s following a `PublicInputGate`.
    /// If some are found, add the corresponding non-routable wire in the `PublicInputGate` to the
    /// partial witness.
    pub(crate) fn copy_buffer_to_pi_gate(&mut self, offset: usize) {
        let pis_wires = self
            .wire_values
            .iter()
            .filter_map(|(t, &v)| match t {
                Target::Wire(Wire { gate: n, input: i })
                    if ((*n > offset) && ((n - offset) % 2 == 1) && (*i < NUM_ADVICE_WIRES)) =>
                {
                    Some((
                        Target::Wire(Wire {
                            gate: n - 1,
                            input: NUM_ROUTED_WIRES + i,
                        }),
                        v,
                    ))
                }
                _ => None,
            })
            .collect::<Vec<_>>();
        self.wire_values.extend(pis_wires);
    }
}

#[derive(Debug, Clone)]
pub struct Witness<F: Field> {
    wire_values: Vec<Vec<F>>,
}

impl<F: Field> Witness<F> {
    pub fn new(wire_values: Vec<Vec<F>>) -> Self {
        Self { wire_values }
    }

    pub fn get(&self, wire: Wire) -> F {
        self.wire_values[wire.gate][wire.input]
    }

    pub fn get_indices(&self, i: usize, j: usize) -> F {
        self.wire_values[i][j]
    }

    pub fn transpose(&self) -> Vec<Vec<F>> {
        transpose(&self.wire_values)
    }

    /// Converts a `PartialWitness` to a a `Witness`.
    /// The partial witness should be sufficiently preprocessed, e.g., it should contain copy constraints.
    pub fn from_partial(pw: &PartialWitness<F>, degree: usize) -> Self {
        let mut wire_values: Vec<Vec<F>> = Vec::new();
        for i in 0..degree {
            let mut gate_i_wires = Vec::new();
            for j in 0..NUM_WIRES {
                let wire = Wire { gate: i, input: j };
                let value = if pw.contains_wire(wire) {
                    pw.get_wire(wire)
                } else {
                    // In our circuit model, a lot of wires are unused. We just set them to zero.
                    F::ZERO
                };
                gate_i_wires.push(value);
            }
            wire_values.push(gate_i_wires);
        }
        Witness::new(wire_values)
    }
}

pub trait WitnessGenerator<F: Field>: 'static + Sync {
    fn dependencies(&self) -> Vec<Target<F>>;

    /// Given a partial witness, return any newly generated values. The caller will merge them in.
    fn generate(
        &self,
        prefixes: &GatePrefixes,
        constants: &[Vec<F>],
        witness: &PartialWitness<F>,
    ) -> PartialWitness<F>;
}
