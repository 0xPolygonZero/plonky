use crate::{CircuitConfig, ConstraintPolynomial, DeterministicGate, Field};

/// A gate which doesn't perform any arithmetic, but just acts as a buffer for receiving data.
/// Some gates, such as the Rescue round gate, "output" their results using one of the next gate's
/// "input" wires. The last such gate has no next gate of the same type, so we add a buffer gate
/// for receiving the last gate's output.
pub struct BufferGate2;

impl<F: Field> DeterministicGate<F> for BufferGate2 {
    fn id(&self) -> String {
        "Buffer".into()
    }

    fn outputs(&self, _config: CircuitConfig) -> Vec<(usize, ConstraintPolynomial<F>)> {
        Vec::new()
    }
}
