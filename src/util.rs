// TODO: Can this impl usize?
pub(crate) fn ceil_div_usize(a: usize, b: usize) -> usize {
    (a + b - 1) / b
}

/// Computes `ceil(log_2(n))`.
pub(crate) fn log2_ceil(n: usize) -> usize {
    n.next_power_of_two().trailing_zeros() as usize
}

/// Computes `log_2(n)`, panicking if `n` is not a power of two.
pub(crate) fn log2_strict(n: usize) -> usize {
    assert!(n.is_power_of_two(), "Not a power of two");
    log2_ceil(n)
}
