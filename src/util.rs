// TODO: Can this impl usize?
pub(crate) fn ceil_div_usize(a: usize, b: usize) -> usize {
    (a + b - 1) / b
}
