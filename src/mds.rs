use crate::Field;

/// Returns entry `(r, c)` of an `n` by `n` MDS matrix.
pub(crate) fn mds<F: Field>(n: usize, r: usize, c: usize) -> F {
    debug_assert!(r < n);
    debug_assert!(c < n);

    // We use a Cauchy matrix with x_r = n + r, y_c = c.
    let x = F::from_canonical_usize(n + r);
    let y = F::from_canonical_usize(c);
    (x - y).multiplicative_inverse().unwrap()
}
