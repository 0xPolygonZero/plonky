use crate::util::transpose;
use crate::{biguint_to_field, biguint_to_limbs, field_to_biguint, AffinePoint, AffinePointTarget, BigIntTarget, Curve, Field, ForeignFieldTarget, OrderingTarget, PublicInput, Target, Wire, LIMB_BITS, NUM_WIRES};
use num::{BigUint, Zero};
use std::{cmp::Ordering, collections::HashMap};

pub struct PartialWitness<F: Field> {
    wire_values: HashMap<Target, F>,
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

    pub fn contains_target(&self, target: Target) -> bool {
        self.wire_values.contains_key(&target)
    }

    pub fn contains_wire(&self, wire: Wire) -> bool {
        self.contains_target(Target::Wire(wire))
    }

    pub fn contains_all_targets(&self, targets: &[Target]) -> bool {
        targets.iter().all(|&t| self.contains_target(t))
    }

    pub fn all_populated_targets(&self) -> Vec<Target> {
        self.wire_values.keys().cloned().collect()
    }

    pub fn get_target(&self, target: Target) -> F {
        self.wire_values[&target]
    }

    pub fn get_targets(&self, targets: &[Target]) -> Vec<F> {
        targets.iter().map(|&t| self.get_target(t)).collect()
    }

    pub fn get_point_target<InnerC: Curve<BaseField = F>>(
        &self,
        target: AffinePointTarget,
    ) -> AffinePoint<InnerC> {
        let x = self.get_target(target.x);
        let y = self.get_target(target.y);
        AffinePoint::nonzero(x, y)
    }

    pub fn get_ordering_target(&self, target: OrderingTarget) -> Ordering {
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

    pub fn set_ordering_target(&mut self, target: OrderingTarget, value: Ordering) {
        let OrderingTarget { lt, eq, gt } = target;
        let values = match value {
            Ordering::Less => [F::ONE, F::ZERO, F::ZERO],
            Ordering::Equal => [F::ZERO, F::ONE, F::ZERO],
            Ordering::Greater => [F::ZERO, F::ZERO, F::ONE],
        };
        self.set_targets(&[lt, eq, gt], &values);
    }

    pub fn get_bigint_target(&self, target: &BigIntTarget) -> BigUint {
        let mut result = BigUint::zero();
        for (i, &limb) in target.limbs.iter().enumerate() {
            let limb_value = field_to_biguint(self.get_target(limb));
            result += limb_value << i * LIMB_BITS;
        }
        result
    }

    pub fn set_bigint_target(&mut self, target: &BigIntTarget, value: &BigUint) {
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

    pub fn get_foreign_field_target<FF: Field>(&self, target: &ForeignFieldTarget) -> FF {
        biguint_to_field(self.get_bigint_target(&target.value))
    }

    pub fn set_foreign_field_target<FF: Field>(&mut self, target: &ForeignFieldTarget, value: FF) {
        self.set_bigint_target(&target.value, &field_to_biguint(value))
    }

    pub fn get_wire(&self, wire: Wire) -> F {
        self.get_target(Target::Wire(wire))
    }

    pub fn set_target(&mut self, target: Target, value: F) {
        let opt_old_value = self.wire_values.insert(target, value);
        if let Some(old_value) = opt_old_value {
            debug_assert_eq!(
                old_value, value,
                "Target {:?} was set twice with different values",
                target
            );
        }
    }

    pub fn set_targets(&mut self, targets: &[Target], values: &[F]) {
        debug_assert_eq!(targets.len(), values.len());
        targets
            .iter()
            .zip(values.iter())
            .for_each(|(&target, &value)| self.set_target(target, value))
    }

    pub fn set_point_target<InnerC: Curve<BaseField = F>>(
        &mut self,
        point_target: AffinePointTarget,
        point: AffinePoint<InnerC>,
    ) {
        self.set_target(point_target.x, point.x);
        self.set_target(point_target.y, point.y);
    }

    pub fn set_point_targets<InnerC: Curve<BaseField = F>>(
        &mut self,
        point_targets: &[AffinePointTarget],
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

    pub fn set_public_input(&mut self, pi: PublicInput, value: F) {
        self.set_wire(pi.original_wire(), value);
        self.set_target(pi.routable_target(), value);
    }

    pub fn extend(&mut self, other: PartialWitness<F>) {
        for (target, value) in other.wire_values {
            self.set_target(target, value);
        }
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
    fn dependencies(&self) -> Vec<Target>;

    /// Given a partial witness, return any newly generated values. The caller will merge them in.
    fn generate(&self, constants: &Vec<Vec<F>>, witness: &PartialWitness<F>) -> PartialWitness<F>;
}
