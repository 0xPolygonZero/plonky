//! This module contains field arithmetic implementations for BLS12's base field and scalar field,
//! with an emphasis on performance.
//!
//! Field elements are represented as `u64` arrays, which works well with modern x86 systems which
//! natively support multiplying two `u64`s to obtain a `u128`. All encodings in the public API are
//! little-endian.
//!
//! We use fixed-length arrays so that there is no need for heap allocation or bounds checking.
//! Unfortunately, this means that some methods need to be duplicated to handle various array
//! lengths. They can be rewritten to use const generics when that feature becomes stable.

use std::cmp::Ordering;
use std::convert::TryInto;
use std::ops::Mul;

/// An element of the BLS12 group's base field.
#[derive(Copy, Clone)]
pub struct Bls12Base {
    /// The limbs in little-endian form.
    limbs: [u64; 6],
}

/// An element of the BLS12 group's scalar field.
#[derive(Copy, Clone)]
pub struct Bls12Scalar {
    /// The limbs in little-endian form.
    limbs: [u64; 4],
}

impl Bls12Base {
    pub const ZERO: Self = Self { limbs: [0; 6] };

    // The order of the field.
    pub const ORDER: [u64; 6] = [13402431016077863595, 2210141511517208575, 7435674573564081700,
        7239337960414712511, 5412103778470702295, 1873798617647539866];

    // Precomputed R for the Barrett reduction algorithm.
    const BARRET_FACTOR: [u64; 6] = [17027978386419893992, 5649138592172459777, 3421924034565217767,
        11848418460761227941, 4080332095855958760, 2837504485842123031];

    const BARRET_K: usize = 381;
}

impl Bls12Scalar {
    pub const ZERO: Self = Self { limbs: [0; 4] };

    // The order of the field.
    pub const ORDER: [u64; 4] = [18446744069414584321, 6034159408538082302, 3691218898639771653, 8353516859464449352];

    // Precomputed R for the Barrett reduction algorithm.
    const BARRET_CONSTANT: [u64; 4] = [5808762262936312036, 15654811016218471260, 1021603728894469044, 10183805594867568095];

    const BARRET_K: usize = 255;
}

impl Mul<Bls12Base> for Bls12Base {
    type Output = Bls12Base;

    fn mul(self, rhs: Bls12Base) -> Bls12Base {
        // First we do a widening multiplication.
        let product = mul_6_6(self.limbs, rhs.limbs);
//        println!("x: {:?}", product);

        // Then, to make it a modular multiplication, we apply the Barrett reduction algorithm.
        // See https://www.nayuki.io/page/barrett-reduction-algorithm
        let product_r = mul_12_6(product, Self::BARRET_FACTOR);
        println!("xr: {:?} * {:?} = {:?}", product, Self::BARRET_FACTOR, product_r);

        // Shift left to divide by 4^k.
        let mut product_r_shifted = [0u64; 6];
        for i in 0..6 {
            let shift_total_bits = Self::BARRET_K * 2;
            let shift_words = shift_total_bits / 64;
            let shift_bits = shift_total_bits as u64 % 64;
            product_r_shifted[i] = product_r[i + shift_words] >> shift_bits
                | product_r[i + shift_words + 1] << (64 - shift_bits);
        }

        let product_r_shifted_n = mul_6_6(product_r_shifted, Self::ORDER);
        let t_12 = sub_12_12(product, product_r_shifted_n);

        // t should fit in m + 1 bits, where m is the modulus length in bits.
        debug_assert!(t_12[6] == 0 || t_12[6] == 1,
                      "t should fit in m + 1 bits: {:?}", t_12);
        for i in 7..12 {
            debug_assert_eq!(t_12[i], 0,
                             "t should fit in m + 1 bits: {:?}", t_12);
        }
        // For efficiency, truncate t down to 7 `u64`s.
        let t_7: [u64; 7] = (&t_12[0..7]).try_into().unwrap();

        let result_7 = if cmp_7_6(t_7, Self::ORDER) == Ordering::Less {
            t_7
        } else {
            sub_7_6(t_7, Self::ORDER)
        };
        // The difference should fit into 6 bits; truncate the most significant bit.
        debug_assert_eq!(t_7[6], 0);
        let result_6 = (&result_7[0..6]).try_into().unwrap();
        Self { limbs: result_6 }
    }
}

fn cmp_6_6(a: [u64; 6], b: [u64; 6]) -> Ordering {
    for i in (0..6).rev() {
        if a[i] < b[i] {
            return Ordering::Less;
        }
        if a[i] > b[i] {
            return Ordering::Greater;
        }
    }
    Ordering::Equal
}

fn cmp_7_6(a: [u64; 7], b: [u64; 6]) -> Ordering {
    if a[6] != 0 {
        return Ordering::Greater;
    }
    for i in (0..6).rev() {
        if a[i] < b[i] {
            return Ordering::Less;
        }
        if a[i] > b[i] {
            return Ordering::Greater;
        }
    }
    Ordering::Equal
}

fn add_6_6(a: [u64; 6], b: [u64; 6]) -> [u64; 7] {
    let mut carry = false;
    let mut sum = [0; 7];
    for i in 0..6 {
        // First add a[i] + b[i], then add in the carry.
        let result1 = a[i].overflowing_add(b[i]);
        let result2 = result1.0.overflowing_add(carry as u64);
        sum[i] = result2.0;
        // If either sum overflowed, set the carry bit.
        carry = result1.1 | result2.1;
    }
    sum[6] = carry as u64;
    sum
}

/// Compute `a - b`. Assumes `a >= b`, otherwise the behavior is undefined.
fn sub_12_12(a: [u64; 12], b: [u64; 12]) -> [u64; 12] {
    let mut borrow = false;
    let mut difference = [0; 12];
    for i in 0..12 {
        let result1 = a[i].overflowing_sub(b[i]);
        let result2 = result1.0.overflowing_sub(borrow as u64);
        difference[i] = result2.0;
        // If either difference underflowed, set the borrow bit.
        borrow = result1.1 | result2.1;
    }
    debug_assert!(!borrow, "a < b: {:?} < {:?}", a, b);
    difference
}

/// Compute `a - b`. Assumes `a >= b`, otherwise the behavior is undefined.
fn sub_6_6(a: [u64; 6], b: [u64; 6]) -> [u64; 6] {
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

/// Compute `a - b`. Assumes `a >= b`, otherwise the behavior is undefined.
fn sub_7_6(a: [u64; 7], b: [u64; 6]) -> [u64; 7] {
    let mut borrow = false;
    let mut difference = [0; 7];
    for i in 0..6 {
        let result1 = a[i].overflowing_sub(b[i]);
        let result2 = result1.0.overflowing_sub(borrow as u64);
        difference[i] = result2.0;
        // If either difference underflowed, set the borrow bit.
        borrow = result1.1 | result2.1;
    }

    // For the last digit, we do `a - carry`, since `b` has been fully consumed already.
    let result = a[6].overflowing_sub(borrow as u64);
    difference[6] = result.0;
    debug_assert!(!result.1, "a < b: {:?} < {:?}", a, b);

    difference
}

fn mul_6_6(a: [u64; 6], b: [u64; 6]) -> [u64; 12] {
    // Grade school multiplication. To avoid carrying at each of O(n^2) steps, we first add each
    // intermediate product to a 128-bit accumulator, then propagate carries at the end.
    let mut acc128 = [0u128; 12];

    for i in 0..6 {
        for j in 0..6 {
            let a_i_b_j = a[i] as u128 * b[j] as u128;
            // Add the less significant chunk to the less significant accumulator.
            acc128[i + j] += a_i_b_j as u64 as u128;
            // Add the more significant chunk to the more significant accumulator.
            acc128[i + j + 1] += a_i_b_j >> 64;
        }
    }

    let mut acc = [0u64; 12];
    acc[0] = acc128[0] as u64;
    let mut carry = false;
    for i in 1..12 {
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

fn mul_12_6(a: [u64; 12], b: [u64; 6]) -> [u64; 18] {
    // Grade school multiplication. To avoid carrying at each of O(n^2) steps, we first add each
    // intermediate product to a 128-bit accumulator, then propagate carries at the end.
    let mut acc128 = [0u128; 18];

    for i in 0..12 {
        for j in 0..6 {
            let a_i_b_j = a[i] as u128 * b[j] as u128;
            // Add the least significant chunk to the less significant accumulator.
            acc128[i + j] += a_i_b_j as u64 as u128;
            // Add the more significant chunk to the more significant accumulator.
            acc128[i + j + 1] += a_i_b_j >> 64;
        }
    }

    let mut acc = [0u64; 18];
    acc[0] = acc128[0] as u64;
    let mut carry = false;
    for i in 1..18 {
        let last_chunk_big = (acc128[i - 1] >> 64) as u64;
        let curr_chunk_small = acc128[i] as u64;
        // Note that last_chunk_big won't get anywhere near 2^64, since it's essentially a carry
        // from some additions in the previous phase, so we can add the carry bit to it without
        // fear of overflow.
        let result = curr_chunk_small.overflowing_add(last_chunk_big + carry as u64);
        acc[i] += result.0;
        carry = result.1;
    }
    acc
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use num::BigUint;

    use crate::field::{Bls12Base, mul_12_6, mul_6_6};

    fn u64_slice_to_biguint(n: &[u64]) -> BigUint {
        let mut bytes_le = Vec::new();
        for n_i in n {
            for j in 0..8 {
                bytes_le.push((n_i >> j * 8) as u8);
            }
        }
        BigUint::from_bytes_le(&bytes_le)
    }

    #[test]
    fn convert_single_u64() {
        let biguint = u64_slice_to_biguint(&[12379813738877118345]);
        assert_eq!(biguint, BigUint::from_str("12379813738877118345").unwrap());
    }

    #[test]
    fn test_mul_6_6() {
        let a = [11111111u64, 22222222, 33333333, 44444444, 55555555, 66666666];
        let b = [77777777u64, 88888888, 99999999, 11111111, 22222222, 33333333];
        assert_eq!(
            u64_slice_to_biguint(&mul_6_6(a, b)),
            u64_slice_to_biguint(&a) * u64_slice_to_biguint(&b));
    }

    #[test]
    fn test_mul_12_6() {
        let a = [11111111u64, 22222222, 33333333, 44444444, 55555555, 66666666,
            77777777, 88888888, 99999999, 11111111, 22222222, 33333333];
        let b = [77777777u64, 88888888, 99999999, 11111111, 22222222, 33333333];
        assert_eq!(
            u64_slice_to_biguint(&mul_12_6(a, b)),
            u64_slice_to_biguint(&a) * u64_slice_to_biguint(&b));
    }

    #[test]
    fn test_mul_bls12_base() {
        let a = [1, 2, 3, 4, 0, 0];
        let b = [3, 4, 5, 6, 0, 0];

        let a_blsbase = Bls12Base { limbs: a };
        let b_blsbase = Bls12Base { limbs: b };
        let a_biguint = u64_slice_to_biguint(&a);
        let b_biguint = u64_slice_to_biguint(&b);
        let order_biguint = u64_slice_to_biguint(&Bls12Base::ORDER);

        println!("{} {}",
                 u64_slice_to_biguint(&(a_blsbase * b_blsbase).limbs),
                 &a_biguint * &b_biguint % &order_biguint);

        assert_eq!(
            u64_slice_to_biguint(&(a_blsbase * b_blsbase).limbs),
            a_biguint * b_biguint % order_biguint);
    }
}
