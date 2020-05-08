use crate::Field;

/// Apply an MDS matrix to the given state vector.
pub(crate) fn apply_mds<F: Field>(inputs: Vec<F>) -> Vec<F> {
    let n = inputs.len();
    let mut result = vec![F::ZERO; n];
    for r in 0..n {
        for c in 0..n {
            result[r] = result[r] + mds::<F>(n, r, c) * inputs[c];
        }
    }
    result
}

/// Returns entry `(r, c)` of an `n` by `n` MDS matrix.
pub(crate) fn mds<F: Field>(n: usize, r: usize, c: usize) -> F {
    debug_assert!(r < n);
    debug_assert!(c < n);

    // We use a Cauchy matrix with x_r = n + r, y_c = c.
    let x = F::from_canonical_usize(n + r);
    let y = F::from_canonical_usize(c);
    (x - y).multiplicative_inverse().unwrap()
}
