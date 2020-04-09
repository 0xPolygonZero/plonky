use crate::Field;

pub(crate) fn generate_rescue_constants<F: Field>(width: usize) -> Vec<Vec<F>> {
    let mut constants = Vec::new();
    for _i in 0..recommended_rounds::<F>(width) {
        let mut round_constants = Vec::new();
        for _j in 0..2 {
            for _k in 0..width {
                // TODO: This should use deterministic randomness.
                round_constants.push(F::rand());
            }
        }
        constants.push(round_constants);
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
    ceil_div(security_bits, 2 * width).max(10)
}

fn ceil_div(a: usize, b: usize) -> usize {
    (a + b - 1) / b
}
