use crate::{apply_mds, CircuitConfig, ConstraintPolynomial, Field, Gate2, PartialWitness2, SimpleGenerator, Target2, Wire, WitnessGenerator2, apply_mds_constraint_polys};

/// Implements a round of the Rescue permutation, modified with a different key schedule to reduce
/// the number of constants involved.
#[derive(Copy, Clone)]
pub struct ModifiedRescueGate {
    width: usize,
    alpha: usize,
}

impl ModifiedRescueGate {
    /// Returns the index of the `i`th accumulator wire. These act as both input and output wires.
    pub fn wire_acc(&self, i: usize) -> usize {
        debug_assert!(i < self.width);
        i
    }

    /// Returns the index of the `i`th (purported) root wire.
    pub fn wire_root(&self, i: usize) -> usize {
        debug_assert!(i < self.width);
        self.width + i
    }
}

impl<F: Field> Gate2<F> for ModifiedRescueGate {
    fn id(&self) -> String {
        format!("ModifiedRescueGate[width={}, alpha={}]", self.width, self.alpha)
    }

    fn constraints(&self, _config: CircuitConfig) -> Vec<ConstraintPolynomial<F>> {
        let w = self.width;

        // Load the input layer variables.
        let layer_0 = (0..w)
            .map(|i| ConstraintPolynomial::local_wire_value(self.wire_acc(i)))
            .collect::<Vec<_>>();

        // Load the (purported) alpha'th root layer variables.
        let layer_1 = (0..w)
            .map(|i| ConstraintPolynomial::local_wire_value(self.wire_root(i)))
            .collect::<Vec<_>>();

        // Compute the input layer from the alpha'th root layer, so that we can verify the
        // (purported) roots.
        let computed_inputs = layer_1.iter()
            .map(|x| x.exp(self.alpha))
            .collect::<Vec<_>>();

        // Apply the MDS matrix.
        let layer_2 = apply_mds_constraint_polys(layer_1);

        // Add a constant to the first element.
        let mut layer_3 = layer_2;
        layer_3[0] = &layer_3[0] + ConstraintPolynomial::local_constant(0);

        // Raise to the alpha'th power.
        let layer_4 = layer_3.iter()
            .map(|x| x.exp(self.alpha))
            .collect::<Vec<_>>();

        // Apply the MDS matrix.
        let layer_5 = apply_mds_constraint_polys(layer_4);

        // Add a constant to the first element.
        let mut layer_6 = layer_5;
        layer_6[0] = &layer_6[0] + ConstraintPolynomial::local_constant(1);

        let mut constraints = Vec::new();
        for i in 0..w {
            // Check that the computed input matches the actual input.
            constraints.push(&computed_inputs[i] - &layer_0[i]);

            // Check that the computed output matches the actual output.
            let actual_out = ConstraintPolynomial::next_wire_value(self.wire_acc(i));
            constraints.push(&layer_6[i] - actual_out);
        }
        constraints
    }

    fn generators(
        &self,
        _config: CircuitConfig,
        gate_index: usize,
        local_constants: Vec<F>,
        _next_constants: Vec<F>,
    ) -> Vec<Box<dyn WitnessGenerator2<F>>> {
        let gen = ModifiedRescueGenerator::<F> {
            gate: *self,
            gate_index,
            constants: local_constants.clone(),
        };
        vec![Box::new(gen)]
    }
}

struct ModifiedRescueGenerator<F: Field> {
    gate: ModifiedRescueGate,
    gate_index: usize,
    constants: Vec<F>,
}

impl<F: Field> SimpleGenerator<F> for ModifiedRescueGenerator<F> {
    fn dependencies(&self) -> Vec<Target2<F>> {
        (0..self.gate.width)
            .map(|i| Target2::Wire(Wire { gate: self.gate_index, input: self.gate.wire_acc(i) }))
            .collect()
    }

    fn run_once(&mut self, witness: &PartialWitness2<F>) -> PartialWitness2<F> {
        let w = self.gate.width;

        // Load inputs.
        let layer_0 = (0..w)
            .map(|i| witness.get_wire(
                Wire { gate: self.gate_index, input: self.gate.wire_acc(i) }))
            .collect::<Vec<_>>();

        // Take alpha'th roots.
        let layer_1 = layer_0.iter()
            .map(|x| x.kth_root_usize(self.gate.alpha))
            .collect::<Vec<_>>();
        let layer_roots = layer_1.clone();

        // Apply MDS matrix.
        let layer_2 = apply_mds(layer_1);

        // Add a constant to the first element.
        let mut layer_3 = layer_2;
        layer_3[0] = layer_3[0] + self.constants[0];

        // Raise to the alpha'th power.
        let layer_4 = layer_3.iter()
            .map(|x| x.exp_usize(self.gate.alpha))
            .collect::<Vec<_>>();

        // Apply MDS matrix.
        let layer_5 = apply_mds(layer_4);

        // Add a constant to the first element.
        let mut layer_6 = layer_5;
        layer_6[0] = layer_6[0] + self.constants[1];

        let mut result = PartialWitness2::new();
        for i in 0..w {
            // Set the i'th root wire.
            result.set_wire(
                Wire { gate: self.gate_index, input: self.gate.wire_root(i) },
                layer_roots[i]);
            // Set the i'th output wire.
            result.set_wire(
                Wire { gate: self.gate_index + 1, input: self.gate.wire_acc(i) },
                layer_6[i]);
        }
        result
    }
}
