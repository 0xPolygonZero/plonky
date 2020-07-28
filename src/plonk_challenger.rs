use std::marker::PhantomData;

use crate::{rescue_permutation, AffinePoint, AffinePointTarget, CircuitBuilder, Curve, Field, HaloCurve, ProjectivePoint, Target, RESCUE_SPONGE_RATE, RESCUE_SPONGE_WIDTH};

/// Observes prover messages, and generates challenges by hashing the transcript.
#[derive(Clone)]
pub struct Challenger<F: Field> {
    sponge_state: Vec<F>,
    input_buffer: Vec<F>,
    output_buffer: Vec<F>,
    security_bits: usize,
}

/// Observes prover messages, and generates verifier challenges based on the transcript.
///
/// The implementation is roughly based on a duplex sponge with a Rescue permutation. Note that in
/// each round, our sponge can absorb an arbitrary number of prover messages and generate an
/// arbitrary number of verifier challenges. This might appear to diverge from the duplex sponge
/// design, but it can be viewed as a duplex sponge whose inputs are sometimes zero (when we perform
/// multiple squeezes) and whose outputs are sometimes ignored (when we perform multiple
/// absorptions). Thus the security properties of a duplex sponge still apply to our design.
impl<F: Field> Challenger<F> {
    pub fn new(security_bits: usize) -> Challenger<F> {
        Challenger {
            sponge_state: vec![F::ZERO; RESCUE_SPONGE_WIDTH],
            input_buffer: Vec::new(),
            output_buffer: Vec::new(),
            security_bits,
        }
    }

    pub fn observe_element(&mut self, element: F) {
        // Any buffered outputs are now invalid, since they wouldn't reflect this input.
        self.output_buffer.clear();

        self.input_buffer.push(element);
    }

    pub fn observe_elements(&mut self, elements: &[F]) {
        for &element in elements {
            self.observe_element(element);
        }
    }

    pub fn observe_affine_point<C: Curve<BaseField = F>>(&mut self, point: AffinePoint<C>) {
        debug_assert!(!point.zero);
        self.observe_element(point.x);
        self.observe_element(point.y);
    }

    pub fn observe_affine_points<C: Curve<BaseField = F>>(&mut self, points: &[AffinePoint<C>]) {
        for &point in points {
            self.observe_affine_point(point);
        }
    }

    pub fn observe_proj_point<C: Curve<BaseField = F>>(&mut self, point: ProjectivePoint<C>) {
        self.observe_affine_point(point.to_affine());
    }

    pub fn observe_proj_points<C: Curve<BaseField = F>>(&mut self, points: &[ProjectivePoint<C>]) {
        self.observe_affine_points(&ProjectivePoint::batch_to_affine(points));
    }

    pub fn get_challenge(&mut self) -> F {
        self.absorb_buffered_inputs();

        if self.output_buffer.is_empty() {
            // Evaluate the permutation to produce `r` new outputs.
            self.sponge_state = rescue_permutation(&self.sponge_state, self.security_bits);
            self.output_buffer = self.sponge_state[0..RESCUE_SPONGE_RATE].to_vec();
        }

        self.output_buffer
            .pop()
            .expect("Output buffer should be non-empty")
    }

    pub fn get_2_challenges(&mut self) -> (F, F) {
        (self.get_challenge(), self.get_challenge())
    }

    pub fn get_3_challenges(&mut self) -> (F, F, F) {
        (
            self.get_challenge(),
            self.get_challenge(),
            self.get_challenge(),
        )
    }

    pub fn get_n_challenges(&mut self, n: usize) -> Vec<F> {
        (0..n).map(|_| self.get_challenge()).collect()
    }

    /// Absorb any buffered inputs. After calling this, the input buffer will be empty.
    fn absorb_buffered_inputs(&mut self) {
        for input_chunk in self.input_buffer.chunks(RESCUE_SPONGE_RATE) {
            // Add the inputs to our sponge state.
            for (i, &input) in input_chunk.iter().enumerate() {
                self.sponge_state[i] = self.sponge_state[i] + input;
            }

            // Apply the permutation.
            self.sponge_state = rescue_permutation(&self.sponge_state, self.security_bits);
        }

        self.output_buffer = self.sponge_state[0..RESCUE_SPONGE_RATE].to_vec();

        self.input_buffer.clear();
    }
}

/// A recursive version of `Challenger`.
pub(crate) struct RecursiveChallenger<C: HaloCurve> {
    sponge_state: Vec<Target>,
    input_buffer: Vec<Target>,
    output_buffer: Vec<Target>,
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> RecursiveChallenger<C> {
    pub(crate) fn new(builder: &mut CircuitBuilder<C>) -> RecursiveChallenger<C> {
        let zero = builder.zero_wire();
        RecursiveChallenger {
            sponge_state: vec![zero; RESCUE_SPONGE_WIDTH],
            input_buffer: Vec::new(),
            output_buffer: Vec::new(),
            _phantom: PhantomData,
        }
    }

    pub(crate) fn observe_element(&mut self, target: Target) {
        // Any buffered outputs are now invalid, since they wouldn't reflect this input.
        self.output_buffer.clear();

        self.input_buffer.push(target);
    }

    pub(crate) fn observe_elements(&mut self, targets: &[Target]) {
        for &target in targets {
            self.observe_element(target);
        }
    }

    pub(crate) fn observe_affine_point(&mut self, point: AffinePointTarget) {
        self.observe_element(point.x);
        self.observe_element(point.y);
    }

    pub(crate) fn observe_affine_points(&mut self, points: &[AffinePointTarget]) {
        for &point in points {
            self.observe_affine_point(point);
        }
    }

    pub(crate) fn get_challenge(&mut self, builder: &mut CircuitBuilder<C>) -> Target {
        self.absorb_buffered_inputs(builder);

        if self.output_buffer.is_empty() {
            // Evaluate the permutation to produce `r` new outputs.
            self.sponge_state = builder.rescue_permutation(&self.sponge_state);
            self.output_buffer = self.sponge_state[0..RESCUE_SPONGE_RATE].to_vec();
        }

        self.output_buffer
            .pop()
            .expect("Output buffer should be non-empty")
    }

    pub(crate) fn get_2_challenges(&mut self, builder: &mut CircuitBuilder<C>) -> (Target, Target) {
        (self.get_challenge(builder), self.get_challenge(builder))
    }

    pub(crate) fn get_3_challenges(
        &mut self,
        builder: &mut CircuitBuilder<C>,
    ) -> (Target, Target, Target) {
        (
            self.get_challenge(builder),
            self.get_challenge(builder),
            self.get_challenge(builder),
        )
    }

    #[allow(dead_code)]
    pub(crate) fn get_n_challenges(
        &mut self,
        builder: &mut CircuitBuilder<C>,
        n: usize,
    ) -> Vec<Target> {
        (0..n).map(|_| self.get_challenge(builder)).collect()
    }

    /// Absorb any buffered inputs. After calling this, the input buffer will be empty.
    fn absorb_buffered_inputs(&mut self, builder: &mut CircuitBuilder<C>) {
        for input_chunk in self.input_buffer.chunks(RESCUE_SPONGE_RATE) {
            // Add the inputs to our sponge state.
            for (i, &input) in input_chunk.iter().enumerate() {
                self.sponge_state[i] = builder.add(self.sponge_state[i], input);
            }

            // Apply the permutation.
            self.sponge_state = builder.rescue_permutation(&self.sponge_state);
        }

        self.output_buffer = self.sponge_state[0..RESCUE_SPONGE_RATE].to_vec();

        self.input_buffer.clear();
    }
}

#[cfg(test)]
mod tests {
    use crate::plonk_challenger::{Challenger, RecursiveChallenger};
    use crate::{CircuitBuilder, Curve, Field, PartialWitness, Target, Tweedledum};

    /// Tests for consistency between `Challenger` and `RecursiveChallenger`.
    #[test]
    fn test_consistency() {
        type C = Tweedledum;
        type SF = <C as Curve>::ScalarField;

        // These are mostly arbitrary, but we want to test some rounds with enough inputs/outputs to
        // trigger multiple absorptions/squeezes.
        let num_inputs_per_round = vec![2, 5, 3];
        let num_outputs_per_round = vec![1, 2, 4];

        // Generate random input messages.
        let inputs_per_round: Vec<Vec<SF>> = num_inputs_per_round
            .iter()
            .map(|&n| (0..n).map(|_| SF::rand()).collect::<Vec<_>>())
            .collect();

        let mut challenger = Challenger::new(128);
        let mut outputs_per_round: Vec<Vec<SF>> = Vec::new();
        for (r, inputs) in inputs_per_round.iter().enumerate() {
            challenger.observe_elements(inputs);
            outputs_per_round.push(challenger.get_n_challenges(num_outputs_per_round[r]));
        }

        let mut builder = CircuitBuilder::<C>::new(128);
        let mut recursive_challenger = RecursiveChallenger::new(&mut builder);
        let mut recursive_outputs_per_round: Vec<Vec<Target>> = Vec::new();
        for (r, inputs) in inputs_per_round.iter().enumerate() {
            recursive_challenger.observe_elements(&builder.constant_wires(inputs));
            recursive_outputs_per_round.push(
                recursive_challenger.get_n_challenges(&mut builder, num_outputs_per_round[r]),
            );
        }
        let circuit = builder.build();
        let witness = circuit.generate_partial_witness(PartialWitness::new());
        let recursive_output_values_per_round: Vec<Vec<SF>> = recursive_outputs_per_round
            .iter()
            .map(|outputs| witness.get_targets(outputs))
            .collect();

        assert_eq!(outputs_per_round, recursive_output_values_per_round);
    }
}
