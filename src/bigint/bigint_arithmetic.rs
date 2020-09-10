use std::cmp::Ordering;
use std::cmp::Ordering::{Equal, Greater, Less};

use rand::rngs::OsRng;
use rand::Rng;
use unroll::unroll_for_loops;

/// This module provides functions for big integer arithmetic using little-endian encoded u64
/// arrays.
#[unroll_for_loops]
pub(crate) fn cmp<const N: usize>(a: [u64; N], b: [u64; N]) -> Ordering {
    for i in (0..N).rev() {
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
pub(crate) fn add_no_overflow<const N: usize>(a: [u64; N], b: [u64; N]) -> [u64; N] {
    let mut carry = false;
    let mut sum = [0; N];
    for i in 0..N {
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
pub(crate) fn sub<const N: usize>(a: [u64; N], b: [u64; N]) -> [u64; N] {
    let mut borrow = false;
    let mut difference = [0; N];
    for i in 0..N {
        let result1 = a[i].overflowing_sub(b[i]);
        let result2 = result1.0.overflowing_sub(borrow as u64);
        difference[i] = result2.0;
        // If either difference underflowed, set the borrow bit.
        borrow = result1.1 | result2.1;
    }

    debug_assert!(!borrow, "a < b: {:?} < {:?}", a, b);
    difference
}

#[inline(always)]
pub(crate) fn is_even<const N: usize>(x: [u64; N]) -> bool {
    x[0] & 1 == 0
}

#[inline(always)]
pub(crate) fn is_odd<const N: usize>(x: [u64; N]) -> bool {
    x[0] & 1 != 0
}

/// Shift in the direction of increasing significance by one bit. Equivalent to integer
/// multiplication by two. Assumes that the input's most significant bit is clear.
#[unroll_for_loops]
pub(crate) fn mul2<const N: usize>(x: [u64; N]) -> [u64; N] {
    debug_assert_eq!(x[N-1] >> 63, 0, "Most significant bit should be clear");

    let mut result = [0; N];
    result[0] = x[0] << 1;
    for i in 1..N {
        result[i] = x[i] << 1 | x[i - 1] >> 63;
    }
    result
}

/// Shift in the direction of increasing significance by one bit. Equivalent to integer
/// division by two.
#[unroll_for_loops]
pub(crate) fn div2<const N: usize>(x: [u64; N]) -> [u64; N] {
    let mut result = [0; N];
    for i in 0..N-1 {
        result[i] = x[i] >> 1 | x[i + 1] << 63;
    }
    result[N-1] = x[N-1] >> 1;
    result
}

pub(crate) fn rand_range<const N: usize>(limit_exclusive: [u64; N]) -> [u64; N] {
    rand_range_from_rng(limit_exclusive, &mut OsRng)
}

// Same as `rand_range` but specifying a RNG (useful when dealing with seeded RNGs).
pub(crate) fn rand_range_from_rng<R: Rng, const N: usize>(limit_exclusive: [u64; N], rng: &mut R) -> [u64; N] {
    // Our approach is to repeatedly generate random u64 arrays until one of them happens to be
    // within the limit. This could take a lot of attempts if the limit has many leading zero bits,
    // though. It is more efficient to generate n-bit random numbers, where n is the number of bits
    // in the limit.
    let bits_to_strip = limit_exclusive[N-1].leading_zeros();

    let mut limbs = [0; N];

    loop {
        for limb_i in &mut limbs {
            *limb_i = rng.next_u64();
        }
        limbs[N-1] >>= bits_to_strip;

        if cmp(limbs, limit_exclusive) == Less {
            return limbs;
        }
    }
}

#[macro_export]
macro_rules! one_array {
    [$type:ty; $N:expr] => {
        {
            let mut tmp: [$type; $N] = [0 as $type; $N];
            tmp[0] = 1;
            tmp
        }
    };
}

#[cfg(test)]
mod tests {
    use crate::div2;

    #[test]
    fn test_div2() {
        assert_eq!(div2([40, 0, 0, 0, 0, 0]), [20, 0, 0, 0, 0, 0]);

        assert_eq!(
            div2([
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
