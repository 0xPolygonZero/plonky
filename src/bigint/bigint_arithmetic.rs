use std::cmp::Ordering;
use std::cmp::Ordering::{Equal, Greater, Less};

use rand::rngs::OsRng;
use rand::Rng;
use unroll::unroll_for_loops;

/// This module provides functions for big integer arithmetic using little-endian encoded u64
/// arrays.
#[unroll_for_loops]
pub(crate) fn cmp_4_4(a: [u64; 4], b: [u64; 4]) -> Ordering {
    for i in (0..4).rev() {
        if a[i] < b[i] {
            return Less;
        }
        if a[i] > b[i] {
            return Greater;
        }
    }
    Equal
}

/// Public only for use with Criterion benchmarks.
#[unroll_for_loops]
pub fn cmp_6_6(a: [u64; 6], b: [u64; 6]) -> Ordering {
    for i in (0..6).rev() {
        if a[i] < b[i] {
            return Less;
        }
        if a[i] > b[i] {
            return Greater;
        }
    }
    Equal
}

/// Computes `a + b`. Assumes that there is no overflow; this is verified only in debug builds.
#[unroll_for_loops]
pub(crate) fn add_4_4_no_overflow(a: [u64; 4], b: [u64; 4]) -> [u64; 4] {
    let mut carry = false;
    let mut sum = [0; 4];
    for i in 0..4 {
        // First add a[i] + b[i], then add in the carry.
        let result1 = a[i].overflowing_add(b[i]);
        let result2 = result1.0.overflowing_add(carry as u64);
        sum[i] = result2.0;
        // If either sum overflowed, set the carry bit.
        carry = result1.1 | result2.1;
    }
    debug_assert_eq!(carry, false, "Addition overflowed");
    sum
}

/// Computes `a + b`. Assumes that there is no overflow; this is verified only in debug builds.
#[unroll_for_loops]
pub(crate) fn add_6_6_no_overflow(a: [u64; 6], b: [u64; 6]) -> [u64; 6] {
    let mut carry = false;
    let mut sum = [0; 6];
    for i in 0..6 {
        // First add a[i] + b[i], then add in the carry.
        let result1 = a[i].overflowing_add(b[i]);
        let result2 = result1.0.overflowing_add(carry as u64);
        sum[i] = result2.0;
        // If either sum overflowed, set the carry bit.
        carry = result1.1 | result2.1;
    }
    debug_assert_eq!(carry, false, "Addition overflowed");
    sum
}

/// Computes `a - b`. Assumes `a >= b`, otherwise the behavior is undefined.
#[unroll_for_loops]
pub(crate) fn sub_4_4(a: [u64; 4], b: [u64; 4]) -> [u64; 4] {
    let mut borrow = false;
    let mut difference = [0; 4];
    for i in 0..4 {
        let result1 = a[i].overflowing_sub(b[i]);
        let result2 = result1.0.overflowing_sub(borrow as u64);
        difference[i] = result2.0;
        // If either difference underflowed, set the borrow bit.
        borrow = result1.1 | result2.1;
    }

    debug_assert!(!borrow, "a < b: {:?} < {:?}", a, b);
    difference
}

/// Computes `a - b`. Assumes `a >= b`, otherwise the behavior is undefined.
#[unroll_for_loops]
pub(crate) fn sub_6_6(a: [u64; 6], b: [u64; 6]) -> [u64; 6] {
    let mut borrow = false;
    let mut difference = [0; 6];

    for i in 0..6 {
        let result1 = a[i].overflowing_sub(b[i]);
        let result2 = result1.0.overflowing_sub(borrow as u64);
        difference[i] = result2.0;
        // If either difference underflowed, set the borrow bit.
        borrow = result1.1 | result2.1;
    }

    debug_assert!(!borrow, "a < b: {:?} < {:?}", a, b);
    difference
}

#[allow(dead_code)]
#[unroll_for_loops]
pub(crate) fn mul_4_4(a: [u64; 4], b: [u64; 4]) -> [u64; 8] {
    // Grade school multiplication. To avoid carrying at each of O(n^2) steps, we first add each
    // intermediate product to a 128-bit accumulator, then propagate carries at the end.
    let mut acc128 = [0u128; 4 + 4];

    for i in 0..4 {
        for j in 0..4 {
            let a_i_b_j = a[i] as u128 * b[j] as u128;
            // Add the less significant chunk to the less significant accumulator.
            acc128[i + j] += a_i_b_j as u64 as u128;
            // Add the more significant chunk to the more significant accumulator.
            acc128[i + j + 1] += a_i_b_j >> 64;
        }
    }

    let mut acc = [0u64; 8];
    acc[0] = acc128[0] as u64;
    let mut carry = false;
    for i in 1..8 {
        let last_chunk_big = (acc128[i - 1] >> 64) as u64;
        let curr_chunk_small = acc128[i] as u64;
        // Note that last_chunk_big won't get anywhere near 2^64, since it's essentially a carry
        // from some additions in the previous phase, so we can add the carry bit to it without
        // fear of overflow.
        let result = curr_chunk_small.overflowing_add(last_chunk_big + carry as u64);
        acc[i] += result.0;
        carry = result.1;
    }
    debug_assert!(!carry);
    acc
}

/// Public only for use with Criterion benchmarks.
#[unroll_for_loops]
pub fn mul_6_6(a: [u64; 6], b: [u64; 6]) -> [u64; 6 + 6] {
    // Grade school multiplication. To avoid carrying at each of O(n^2) steps, we first add each
    // intermediate product to a 128-bit accumulator, then propagate carries at the end.
    let mut acc128 = [0u128; 6 + 6];

    for i in 0..6 {
        for j in 0..6 {
            let a_i_b_j = a[i] as u128 * b[j] as u128;
            // Add the less significant chunk to the less significant accumulator.
            acc128[i + j] += a_i_b_j as u64 as u128;
            // Add the more significant chunk to the more significant accumulator.
            acc128[i + j + 1] += a_i_b_j >> 64;
        }
    }

    let mut acc = [0u64; 6 + 6];
    acc[0] = acc128[0] as u64;
    let mut carry = false;
    for i in 1..(6 + 6) {
        let last_chunk_big = (acc128[i - 1] >> 64) as u64;
        let curr_chunk_small = acc128[i] as u64;
        // Note that last_chunk_big won't get anywhere near 2^64, since it's essentially a carry
        // from some additions in the previous phase, so we can add the carry bit to it without
        // fear of overflow.
        let result = curr_chunk_small.overflowing_add(last_chunk_big + carry as u64);
        acc[i] += result.0;
        carry = result.1;
    }
    debug_assert!(!carry);
    acc
}

#[inline(always)]
pub(crate) fn is_even_4(x: [u64; 4]) -> bool {
    x[0] & 1 == 0
}

#[inline(always)]
pub(crate) fn is_even_6(x: [u64; 6]) -> bool {
    x[0] & 1 == 0
}

#[inline(always)]
pub(crate) fn is_odd_4(x: [u64; 4]) -> bool {
    x[0] & 1 != 0
}

#[inline(always)]
pub(crate) fn is_odd_6(x: [u64; 6]) -> bool {
    x[0] & 1 != 0
}

/// Shift in the direction of increasing significance by one bit. Equivalent to integer
/// multiplication by two. Assumes that the input's most significant bit is clear.
#[unroll_for_loops]
pub(crate) fn mul2_6(x: [u64; 6]) -> [u64; 6] {
    debug_assert_eq!(x[5] >> 63, 0, "Most significant bit should be clear");

    let mut result = [0; 6];
    result[0] = x[0] << 1;
    for i in 1..6 {
        result[i] = x[i] << 1 | x[i - 1] >> 63;
    }
    result
}

/// Shift in the direction of increasing significance by one bit. Equivalent to integer
/// division by two.
#[unroll_for_loops]
pub(crate) fn div2_4(x: [u64; 4]) -> [u64; 4] {
    let mut result = [0; 4];
    for i in 0..3 {
        result[i] = x[i] >> 1 | x[i + 1] << 63;
    }
    result[3] = x[3] >> 1;
    result
}

/// Shift in the direction of increasing significance by one bit. Equivalent to integer
/// division by two.
#[unroll_for_loops]
pub(crate) fn div2_6(x: [u64; 6]) -> [u64; 6] {
    let mut result = [0; 6];
    for i in 0..5 {
        result[i] = x[i] >> 1 | x[i + 1] << 63;
    }
    result[5] = x[5] >> 1;
    result
}

pub(crate) fn rand_range_6(limit_exclusive: [u64; 6]) -> [u64; 6] {
    rand_range_6_from_rng(limit_exclusive, &mut OsRng)
}

// Same as `rand_range_6` but specifying a RNG (useful when dealing with seeded RNGs).
pub(crate) fn rand_range_6_from_rng<R: Rng>(limit_exclusive: [u64; 6], rng: &mut R) -> [u64; 6] {
    // Our approach is to repeatedly generate random u64 arrays until one of them happens to be
    // within the limit. This could take a lot of attempts if the limit has many leading zero bits,
    // though. It is more efficient to generate n-bit random numbers, where n is the number of bits
    // in the limit.
    let bits_to_strip = limit_exclusive[5].leading_zeros();

    let mut limbs = [0; 6];

    loop {
        for limb_i in &mut limbs {
            *limb_i = rng.next_u64();
        }
        limbs[5] >>= bits_to_strip;

        if cmp_6_6(limbs, limit_exclusive) == Less {
            return limbs;
        }
    }
}

pub(crate) fn rand_range_4(limit_exclusive: [u64; 4]) -> [u64; 4] {
    rand_range_4_from_rng(limit_exclusive, &mut OsRng)
}

// Same as `rand_range_4` but specifying a RNG (useful when dealing with seeded RNGs).
pub(crate) fn rand_range_4_from_rng<R: Rng>(limit_exclusive: [u64; 4], rng: &mut R) -> [u64; 4] {
    // Our approach is to repeatedly generate random u64 arrays until one of them happens to be
    // within the limit. This could take a lot of attempts if the limit has many leading zero bits,
    // though. It is more efficient to generate n-bit random numbers, where n is the number of bits
    // in the limit.
    let bits_to_strip = limit_exclusive[3].leading_zeros();

    let mut limbs = [0; 4];

    loop {
        for limb_i in &mut limbs {
            *limb_i = rng.next_u64();
        }
        limbs[3] >>= bits_to_strip;

        if cmp_4_4(limbs, limit_exclusive) == Less {
            return limbs;
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::conversions::u64_slice_to_biguint;
    use crate::{div2_6, mul_6_6};

    #[test]
    fn test_mul_6_6() {
        let a = [
            11111111u64,
            22222222,
            33333333,
            44444444,
            55555555,
            66666666,
        ];
        let b = [
            77777777u64,
            88888888,
            99999999,
            11111111,
            22222222,
            33333333,
        ];
        assert_eq!(
            u64_slice_to_biguint(&mul_6_6(a, b)),
            u64_slice_to_biguint(&a) * u64_slice_to_biguint(&b)
        );
    }

    #[test]
    fn test_div2() {
        assert_eq!(div2_6([40, 0, 0, 0, 0, 0]), [20, 0, 0, 0, 0, 0]);

        assert_eq!(
            div2_6([
                15668009436471190370,
                3102040391300197453,
                4166322749169705801,
                3518225024268476800,
                11231577158546850254,
                226224965816356276
            ]),
            [
                17057376755090370993,
                10774392232504874534,
                2083161374584852900,
                1759112512134238400,
                5615788579273425127,
                113112482908178138
            ]
        );
    }
}
