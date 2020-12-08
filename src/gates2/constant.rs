use crate::{ConstraintPolynomial, DeterministicGate, Field, CircuitConfig, GateRef, DeterministicGateAdapter};

/// A gate which takes a single constant parameter and outputs that value.
pub struct ConstantGate2;

impl ConstantGate2 {
    pub fn get_ref<F: Field>() -> GateRef<F> {
        GateRef::new(DeterministicGateAdapter::new(ConstantGate2))
    }

    pub const CONST_INPUT: usize = 0;

    pub const WIRE_OUTPUT: usize = 0;
}

impl<F: Field> DeterministicGate<F> for ConstantGate2 {
    fn id(&self) -> String {
        "ConstantGate".into()
    }

    fn outputs(&self, _config: CircuitConfig) -> Vec<(usize, ConstraintPolynomial<F>)> {
        let out = ConstraintPolynomial::local_constant(Self::CONST_INPUT);
        vec![(Self::WIRE_OUTPUT, out)]
    }
}
