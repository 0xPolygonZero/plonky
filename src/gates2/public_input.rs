use crate::{CircuitConfig, ConstraintPolynomial, Field, Gate2, GateRef, WitnessGenerator2};

/// A gate for receiving public inputs. These gates will be placed at static indices and the wire
/// polynomials will always be opened at those indices.
///
/// Each `PublicInputGate` can receive `max(num_wires, 2 * num_routed_wires)` public inputs. If all
/// wires are routed, each `PublicInputGate` simply receives `num_wires` public inputs, but if not,
/// it gets a bit more complex. We place a `BufferGate` after each `PublicInputGate`, then "copy"
/// any non-routed public inputs from the `PublicInputGate` to the routed wires of the following
/// `BufferGate`. These routed `BufferGate` wires can then be used to route public inputs elsewhere.
pub struct PublicInputGate2;

impl<F: Field> Gate2<F> for PublicInputGate2 {
    fn id(&self) -> String {
        "PublicInputGate".into()
    }

    fn constraints(&self, config: CircuitConfig) -> Vec<ConstraintPolynomial<F>> {
        let routed_pis = config.num_routed_wires;
        let non_routed_pis = config.advice_wires().min(routed_pis);

        // For each non-routed PI, we "copy" that PI to the following gate, which should be a
        // `BufferGate` just for receiving these these PIs and making them routable.
        (0..non_routed_pis).map(|i| {
            let non_routed_pi_wire = ConstraintPolynomial::local_wire_value(routed_pis + i);
            let routed_receiving_wire = ConstraintPolynomial::next_wire_value(i);
            non_routed_pi_wire - routed_receiving_wire
        }).collect()
    }

    fn generators(
        &self,
        _config: CircuitConfig,
        _gate_index: usize,
        _local_constants: Vec<F>,
        _next_constants: Vec<F>,
    ) -> Vec<Box<dyn WitnessGenerator2<F>>> {
        // CircuitBuilder handles copying public input values around.
        Vec::new()
    }
}

impl PublicInputGate2 {
    pub fn get_ref<F: Field>() -> GateRef<F> {
        GateRef::new(PublicInputGate2)
    }
}
