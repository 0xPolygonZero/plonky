use crate::{CircuitConfig, ConstraintPolynomial, DeterministicGate, DeterministicGateAdapter, Field, GateRef};

/// A gate which takes as inputs limbs of sum small base, verifies that each limb is in `[0, base)`,
/// and outputs the weighted sum `limb[0] + base limb[1] + base^2 limb[2] + ...`.
pub struct LimbSumGate {
    base: usize,
    num_limbs: usize,
}

impl LimbSumGate {
    pub fn get_ref<F: Field>(base: usize, num_limbs: usize) -> GateRef<F> {
        let gate = LimbSumGate { base, num_limbs };
        GateRef::new(DeterministicGateAdapter::new(gate))
    }
}

impl<F: Field> DeterministicGate<F> for LimbSumGate {
    fn id(&self) -> String {
        format!("LimbSumGate[base={}, num_limbs={}]", self.base, self.num_limbs)
    }

    fn outputs(&self, _config: CircuitConfig) -> Vec<(usize, ConstraintPolynomial<F>)> {
        // We compute `out = limb[0] + base * limb[1] + base^2 * limb[2] + ...`.
        let out = (0..self.num_limbs).map(|i| {
            let limb = ConstraintPolynomial::local_wire_value(i);
            let weight = F::from_canonical_usize(self.base).exp_usize(i);
            limb * weight
        }).sum();

        vec![(self.num_limbs, out)]
    }

    fn additional_constraints(&self, _config: CircuitConfig) -> Vec<ConstraintPolynomial<F>> {
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
    use crate::{CircuitBuilder2, CircuitConfig, TweedledumBase};

    fn valid() {
        let config = CircuitConfig { num_wires: 3, num_routed_wires: 3, security_bits: 128 };
        let mut builder = CircuitBuilder2::<TweedledumBase>::new(config);
        let zero = builder.zero();
        let one = builder.one();
        let two = builder.two();
        todo!()
    }
}
