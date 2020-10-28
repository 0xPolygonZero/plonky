#[inline]
#[unroll_for_loops]
fn mul_4_4(a: [u64; 4], b: [u64; 4]) -> [u64; 8] {
    // Grade school multiplication. To avoid carrying at each of
    // O(n^2) steps, we first add each intermediate product to a
    // 128-bit accumulator, then propagate carries at the end.
    let mut acc128 = [0u128; 4 + 4];

    for i in 0..4 {
        for j in 0..4 {
            let a_i_b_j = a[i] as u128 * b[j] as u128;
            // Add the less significant chunk to the less significant
            // accumulator.
            acc128[i + j] += a_i_b_j as u64 as u128;
            // Add the more significant chunk to the more significant
            // accumulator.
            acc128[i + j + 1] += a_i_b_j >> 64;
        }
    }

    let mut acc = [0u64; 8];
    acc[0] = acc128[0] as u64;
    let mut carry = false;
    for i in 1..8 {
        let last_chunk_big = (acc128[i - 1] >> 64) as u64;
        let curr_chunk_small = acc128[i] as u64;
        // Note that last_chunk_big won't get anywhere near 2^64,
        // since it's essentially a carry from some additions in the
        // previous phase, so we can add the carry bit to it without
        // fear of overflow.
        let result = curr_chunk_small.overflowing_add(
            last_chunk_big + carry as u64);
        acc[i] += result.0;
        carry = result.1;
    }
    debug_assert!(!carry);
    acc
}

/// This modular arithmetic representation is based on Daira's
/// adaptation of the tricks for moduli of the form 2^b + c.
/// Source: https://hackmd.io/drzN-z-_So28zDLhK2tegw
pub trait DairaRepr {

}
