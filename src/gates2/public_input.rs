use crate::{Field, Gate2, WitnessGenerator2, ConstraintPolynomial};

struct PublicInputGate2 {
    pis_per_gate: usize,
    routed_pis_per_gate: usize,
}
