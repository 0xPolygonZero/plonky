use crate::Field;
use crate::util::ceil_div_usize;

pub(crate) fn generate_rescue_constants<F: Field>(width: usize) -> Vec<(Vec<F>, Vec<F>)> {
    // TODO: This should use deterministic randomness.
    let mut constants = Vec::new();
    for _i in 0..recommended_rounds::<F>(width) {
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

pub(crate) fn recommended_rounds<F: Field>(width: usize) -> usize {
    recommended_rounds_for_security_bits::<F>(width, 128)
}

fn recommended_rounds_for_security_bits<F: Field>(
    width: usize,
    security_bits: usize,
) -> usize {
    ceil_div_usize(security_bits, 2 * width).max(10)
}
