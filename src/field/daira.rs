use unroll::unroll_for_loops;
use crate::{add_no_overflow, sub, cmp, mul2};

// FIXME: These functions are copypasta from monty.rs
#[inline]
fn mul_add_cy_in(a: u64, b: u64, c: u64, cy_in: u64) -> (u64, u64) {
    let t = (a as u128) * (b as u128) + (c as u128) + (cy_in as u128);
    ((t >> 64) as u64, t as u64)
}

#[inline]
fn mul_add(a: u64, b: u64, c: u64) -> (u64, u64) {
    let t = (a as u128) * (b as u128) + (c as u128);
    ((t >> 64) as u64, t as u64)
}

#[unroll_for_loops]
fn sqr_4(a: [u64; 4]) -> [u64; 8] {
    let mut res = [u64; 8];

    // Calculate the off-diagonal part of the square
    for i in 0 .. 4 {
        let mut hi_in = 0u64;
        for j in i+1 .. 4 {
            let (hi_out, lo) = mul_add(a[j], a[i], hi_in);
            res[i + j] = lo;
            hi_in = hi_out;
        }
        res[i + 4] = hi_in;
    }
    res = mul2(res); // NB: Top and bottom words are zero

    // Calculate and add in the diagonal part
    let mut hi_in = 0u64;
    for i in 0 .. 4 {
        let (hi_out, lo) = mul_add_cy_in(a[i], a[i], res[2*i], hi_in);
        res[2*i] = lo;
        let (t, cy) = res[2*i + 1].overflowing_add(hi_out);
        res[2*i + 1] = t;
        hi_in = cy as u64;
    }
    res
}

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
