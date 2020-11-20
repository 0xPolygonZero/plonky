use crate::{DeterministicGate, Field, ConstraintPolynomial};

/// A gate which takes as inputs limbs of sum small base, verifies that each limb is in `[0, base)`,
/// and outputs the weighted sum `limb[0] + base limb[1] + base^2 limb[2] + ...`.
pub struct LimbSumGate {
    base: usize,
    num_limbs: usize,
}

impl<F: Field> DeterministicGate<F> for LimbSumGate {
    fn id(&self) -> String {
        format!("LimbSumGate-{}x{}", self.base, self.num_limbs)
    }

    fn outputs(&self) -> Vec<(usize, ConstraintPolynomial<F>)> {
        // We compute `out = limb[0] + base * limb[1] + base^2 * limb[2] + ...`.
        let out = (0..self.num_limbs).map(|i| {
            let limb = ConstraintPolynomial::local_wire_value(i);
            let weight = F::from_canonical_usize(self.base).exp_usize(i);
            limb * weight
        }).sum();

        vec![(self.num_limbs, out)]
    }

    fn additional_constraints(&self) -> Vec<ConstraintPolynomial<F>> {
        // For each limb,
        (0..self.num_limbs).map(|i| {
            let limb = ConstraintPolynomial::local_wire_value(i);

            // Assert that this limb is in `[0, base)` by enforcing that
            // `limb (limb - 1) .. (limb - (base - 1)) = 0`.
            (0..self.base).map(|possible_value| {
                &limb - possible_value
            }).product()
        }).collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::{CircuitBuilder2, TweedledumBase};

    fn valid() {
        let mut builder = CircuitBuilder2::<TweedledumBase>::new();
        let zero = builder.zero();
        let one = builder.one();
        let two = builder.two();
        todo!()
    }
}