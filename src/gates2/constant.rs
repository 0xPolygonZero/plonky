use crate::{Field, DeterministicGate, ConstraintPolynomial};

/// A gate which takes a single constant parameter and outputs that value.
pub struct ConstantGate2;

impl<F: Field> DeterministicGate<F> for ConstantGate2 {
    fn id(&self) -> String {
        "ConstantGate".into()
    }

    fn outputs(&self) -> Vec<(usize, ConstraintPolynomial<F>)> {
        let out = ConstraintPolynomial::local_constant(0);
        vec![(0, out)]
    }
}
