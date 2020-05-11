use crate::{Field, apply_mds, PRF};
use crate::util::ceil_div_usize;

pub struct RescuePrf {
    security_bits: usize,
}

impl Default for RescuePrf {
    fn default() -> Self {
        RescuePrf {
            security_bits: 128,
        }
    }
}

impl<F: Field> PRF<F> for RescuePrf {
    fn rand(&self, x: F) -> F {
        rescue_hash_1_to_1(x, self.security_bits)
    }
}

pub fn rescue_hash_1_to_1<F: Field>(input: F, security_bits: usize) -> F {
    rescue_hash_n_to_1(vec![input], security_bits)
}

pub fn rescue_hash_n_to_1<F: Field>(inputs: Vec<F>, security_bits: usize) -> F {
    rescue_sponge(inputs, 1, security_bits)[0]
}

pub fn rescue_sponge<F: Field>(
    inputs: Vec<F>,
    num_outputs: usize,
    security_bits: usize,
) -> Vec<F> {
    // This is mostly arbitrary, but we wouldn't want a huge width as the MDS layer could get
    // expensive.
    let rate = 2;
    let capacity = 1;
    let width = rate + capacity;

    let mut state = vec![F::ZERO; width];

    // Absorb all input chunks.
    for input_chunk in inputs.chunks(rate) {
        for i in 0..input_chunk.len() {
            state[i] = state[i] + input_chunk[i];
        }
        state = rescue_permutation(state, security_bits);
    }

    // Squeeze until we have the desired number of outputs.
    let mut outputs = Vec::new();
    loop {
        for i in 0..rate {
            outputs.push(state[i]);
            if outputs.len() == num_outputs {
                return outputs;
            }
        }
        state = rescue_permutation(state, security_bits);
    }
}

pub fn rescue_permutation<F: Field>(mut state: Vec<F>, security_bits: usize) -> Vec<F> {
    let width = state.len();
    let constants = generate_rescue_constants(width, security_bits);

    for (step_a_constants, step_b_constants) in constants {
        // Step A.
        state = state.iter().map(|x| x.exp(F::ALPHA)).collect();
        state = apply_mds(state);
        state = add_vecs(state, step_a_constants);

        // Step B.
        state = state.iter().map(|x| x.kth_root(F::ALPHA)).collect();
        state = apply_mds(state);
        state = add_vecs(state, step_b_constants);
    }

    state
}

fn add_vecs<F: Field>(a: Vec<F>, b: Vec<F>) -> Vec<F> {
    a.iter().zip(b.iter())
        .map(|(a_i, b_i)| *a_i + *b_i)
        .collect()
}

pub(crate) fn generate_rescue_constants<F: Field>(
    width: usize,
    security_bits: usize,
) -> Vec<(Vec<F>, Vec<F>)> {
    // TODO: This should use deterministic randomness.
    let mut constants = Vec::new();
    for _i in 0..recommended_rounds::<F>(width, security_bits) {
        let mut step_a_constants = Vec::new();
        for _k in 0..width {
            step_a_constants.push(F::rand());
        }

        let mut step_b_constants = Vec::new();
        for _k in 0..width {
            step_b_constants.push(F::rand());
        }

        constants.push((step_a_constants, step_b_constants));
    }
    constants
}

pub(crate) fn recommended_rounds<F: Field>(
    width: usize,
    security_bits: usize,
) -> usize {
    ceil_div_usize(security_bits, 2 * width).max(10)
}
