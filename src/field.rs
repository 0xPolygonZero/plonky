//! This module contains field arithmetic implementations for BLS12's base field and scalar field,
//! with an emphasis on performance.
//!
//! Field elements are represented as `u64` arrays, which works well with modern x86 systems which
//! natively support multiplying two `u64`s to obtain a `u128`. All encodings in the public API are
//! little-endian.
//!
//! We use fixed-length arrays so that there is no need for heap allocation or bounds checking.
//! Unfortunately, this means that we need several variants of each function to handle different
//! array sizes. For now, we use macros to generate these variants. This API clutter can be removed
//! in the future when const generics become stable.

use std::cmp::Ordering;
use std::convert::TryInto;
use std::ops::{Mul, Add, Neg, Sub, Div};

/// An element of the BLS12 group's base field.
#[derive(Copy, Clone, Eq, PartialEq)]
pub struct Bls12Base {
    /// The limbs in little-endian form.
    pub limbs: [u64; 6],
}

/// An element of the BLS12 group's scalar field.
#[derive(Copy, Clone, Eq, PartialEq)]
pub struct Bls12Scalar {
    /// The limbs in little-endian form.
    pub limbs: [u64; 4],
}

impl Bls12Base {
    pub const ZERO: Self = Self { limbs: [0; 6] };

    // The order of the field.
    pub const ORDER: [u64; 6] = [13402431016077863595, 2210141511517208575, 7435674573564081700,
        7239337960414712511, 5412103778470702295, 1873798617647539866];

    // Precomputed R for the Barrett reduction algorithm.
    const BARRETT_R: [u64; 6] = [17027978386419893992, 5649138592172459777, 3421924034565217767,
        11848418460761227941, 4080332095855958760, 2837504485842123031];

    const BARRETT_K: usize = 381;
}

impl Bls12Scalar {
    pub const ZERO: Self = Self { limbs: [0; 4] };

    // The order of the field.
    pub const ORDER: [u64; 4] = [18446744069414584321, 6034159408538082302, 3691218898639771653, 8353516859464449352];

    // Precomputed R for the Barrett reduction algorithm.
    const BARRETT_R: [u64; 4] = [5808762262936312036, 15654811016218471260, 1021603728894469044, 10183805594867568095];

    const BARRETT_K: usize = 255;
}

impl Add<Bls12Scalar> for Bls12Scalar {
    type Output = Bls12Scalar;

    fn add(self, rhs: Bls12Scalar) -> Bls12Scalar {
        // First we do a widening addition, then we reduce if necessary.
        let sum = add_4_4(self.limbs, rhs.limbs);
        let result_5 = match cmp_5_4(sum, Bls12Scalar::ORDER) {
            Ordering::Less => sum,
            _ => sub_5_4(sum, Bls12Scalar::ORDER)
        };
        debug_assert_eq!(result_5[4], 0, "reduced sum should fit in 4 u64s");
        let limbs = (&result_5[0..4]).try_into().unwrap();
        Bls12Scalar { limbs }
    }
}

impl Add<Bls12Base> for Bls12Base {
    type Output = Bls12Base;

    fn add(self, rhs: Bls12Base) -> Bls12Base {
        // First we do a widening addition, then we reduce if necessary.
        let sum = add_6_6(self.limbs, rhs.limbs);
        let result_7 = match cmp_7_6(sum, Bls12Base::ORDER) {
            Ordering::Less => sum,
            _ => sub_7_6(sum, Bls12Base::ORDER)
        };
        debug_assert_eq!(result_7[6], 0, "reduced sum should fit in 6 u64s");
        let limbs = (&result_7[0..6]).try_into().unwrap();
        Bls12Base { limbs }
    }
}

impl Sub<Bls12Base> for Bls12Base {
    type Output = Bls12Base;

    fn sub(self, rhs: Bls12Base) -> Self::Output {
        self + -rhs
    }
}

impl Mul<Bls12Scalar> for Bls12Scalar {
    type Output = Bls12Scalar;

    fn mul(self, rhs: Bls12Scalar) -> Bls12Scalar {
        // First we do a widening multiplication, then a modular reduction.
        let product = mul_4_4(self.limbs, rhs.limbs);
        let limbs = barrett_reduction_bls12_scalar(product);
        Bls12Scalar { limbs }
    }
}

impl Mul<Bls12Base> for Bls12Base {
    type Output = Bls12Base;

    fn mul(self, rhs: Bls12Base) -> Bls12Base {
        // First we do a widening multiplication, then a modular reduction.
        let product = mul_6_6(self.limbs, rhs.limbs);
        let limbs = barrett_reduction_bls12_base(product);
        Bls12Base { limbs }
    }
}

impl Mul<u64> for Bls12Base {
    type Output = Bls12Base;

    fn mul(self, rhs: u64) -> Self::Output {
        let rhs_field = Bls12Base{limbs: [rhs, 0x0, 0x0, 0x0, 0x0, 0x0]};
        self * rhs_field
    }
}

impl Neg for Bls12Base {
    type Output = Bls12Base;

    fn neg(self) -> Self::Output {
        if self == Bls12Base::ZERO {
            Bls12Base::ZERO
        }else {
            Bls12Base{limbs: sub_6_6(Bls12Base::ORDER, self.limbs)}
        }
    }
}

macro_rules! barrett_reduction {
    (
        $name:ident,
        $n:expr,
        $mul_n_n:ident,
        $mul_2n_n:ident,
        $cmp_n_plus_1_n:ident,
        $sub_2n_2n:ident,
        $sub_n_plus_1_n:ident,
        $modulus:expr,
        $barrett_r:expr,
        $barrett_k:expr
    ) => {
        fn $name(x: [u64; 2 * $n]) -> [u64; $n] {
            // Then, to make it a modular multiplication, we apply the Barrett reduction algorithm.
            // See https://www.nayuki.io/page/barrett-reduction-algorithm
            let x_r = $mul_2n_n(x, $barrett_r);

            // Shift left to divide by 4^k.
            let mut x_r_shifted = [0u64; $n];
            for i in 0..$n {
                let shift_total_bits = $barrett_k * 2;
                let shift_words = shift_total_bits / 64;
                let shift_bits = shift_total_bits as u64 % 64;
                x_r_shifted[i] = x_r[i + shift_words] >> shift_bits
                    | x_r[i + shift_words + 1] << (64 - shift_bits);
            }

            let x_r_shifted_n = $mul_n_n(x_r_shifted, $modulus);
            let t_2n = $sub_2n_2n(x, x_r_shifted_n);

            // t should fit in m + 1 bits, where m is the modulus length in bits.
            debug_assert!(t_2n[$n] == 0 || t_2n[$n] == 1,
                          "t should fit in m + 1 bits: {:?}", t_2n);
            for i in ($n + 1)..(2 * $n) {
                debug_assert_eq!(t_2n[i], 0,
                                 "t should fit in m + 1 bits: {:?}", t_2n);
            }
            // For efficiency, truncate t down to `$n + 1` `u64`s.
            let t_n_plus_1: [u64; $n + 1] = (&t_2n[0..($n + 1)]).try_into().unwrap();

            let result_n_plus_1 = if $cmp_n_plus_1_n(t_n_plus_1, $modulus) == Ordering::Less {
                t_n_plus_1
            } else {
                $sub_n_plus_1_n(t_n_plus_1, $modulus)
            };

            // The difference should fit into 6 bits; truncate the most significant bit.
            debug_assert_eq!(t_n_plus_1[$n], 0);
            (&result_n_plus_1[0..$n]).try_into().unwrap()
        }
    }
}

macro_rules! cmp_symmetric {
    ($name:ident, $len:expr) => {
        cmp_asymmetric!($name, $len, $len);
    }
}

/// Generates comparison functions for `u64` arrays.
macro_rules! cmp_asymmetric {
    ($name:ident, $a_len:expr, $b_len:expr) => {
        fn $name(a: [u64; $a_len], b: [u64; $b_len]) -> Ordering {
            // If any of the "a only" bits are set, then a is greater, as b's associated bit is
            // implicitly zero.
            for i in $b_len..$a_len {
                if a[i] != 0 {
                    return Ordering::Greater;
                }
            }

            for i in (0..$b_len).rev() {
                if a[i] < b[i] {
                    return Ordering::Less;
                }
                if a[i] > b[i] {
                    return Ordering::Greater;
                }
            }

            Ordering::Equal
        }
    }
}

/// Generates addition functions for `u64` arrays.
macro_rules! add_symmetric {
    ($name:ident, $len:expr) => {
        /// Computes `a + b`.
        fn $name(a: [u64; $len], b: [u64; $len]) -> [u64; $len + 1] {
            let mut carry = false;
            let mut sum = [0; $len + 1];
            for i in 0..$len {
                // First add a[i] + b[i], then add in the carry.
                let result1 = a[i].overflowing_add(b[i]);
                let result2 = result1.0.overflowing_add(carry as u64);
                sum[i] = result2.0;
                // If either sum overflowed, set the carry bit.
                carry = result1.1 | result2.1;
            }
            sum[$len] = carry as u64;
            sum
        }
    }
}

/// Generates subtraction functions for `u64` arrays where the operand lengths are equal.
macro_rules! sub_symmetric {
    ($name:ident, $len:expr) => {
        /// Computes `a - b`. Assumes `a >= b`, otherwise the behavior is undefined.
        fn $name(a: [u64; $len], b: [u64; $len]) -> [u64; $len] {
            debug_assert!($len == $len);

            let mut borrow = false;
            let mut difference = [0; $len];
            for i in 0..$len {
                let result1 = a[i].overflowing_sub(b[i]);
                let result2 = result1.0.overflowing_sub(borrow as u64);
                difference[i] = result2.0;
                // If either difference underflowed, set the borrow bit.
                borrow = result1.1 | result2.1;
            }

            debug_assert!(!borrow, "a < b: {:?} < {:?}", a, b);
            difference
        }
    }
}

/// Generates subtraction functions for `u64` arrays where the length of the first operand exceeds
/// that of the second.
macro_rules! sub_asymmetric {
    ($name:ident, $a_len:expr, $b_len:expr) => {
        /// Computes `a - b`. Assumes `a >= b`, otherwise the behavior is undefined.
        fn $name(a: [u64; $a_len], b: [u64; $b_len]) -> [u64; $a_len] {
            debug_assert!($a_len > $b_len);

            let mut borrow = false;
            let mut difference = [0; $a_len];

            for i in 0..$b_len {
                let result1 = a[i].overflowing_sub(b[i]);
                let result2 = result1.0.overflowing_sub(borrow as u64);
                difference[i] = result2.0;
                // If either difference underflowed, set the borrow bit.
                borrow = result1.1 | result2.1;
            }

            // At this point `b` has been fully consumed, so we just subtract carry bits from digits
            // of `a`.
            for i in $b_len..$a_len {
                let result = a[i].overflowing_sub(borrow as u64);
                difference[i] = result.0;
                borrow = result.1;
            }

            debug_assert!(!borrow, "a < b: {:?} < {:?}", a, b);
            difference
        }
    }
}

macro_rules! mul_symmetric {
    ($name:ident, $len:expr) => {
        mul_asymmetric!($name, $len, $len);
    }
}

macro_rules! mul_asymmetric {
    ($name:ident, $a_len:expr, $b_len:expr) => {
        pub fn $name(a: [u64; $a_len], b: [u64; $b_len]) -> [u64; $a_len + $b_len] {
            // Grade school multiplication. To avoid carrying at each of O(n^2) steps, we first add each
            // intermediate product to a 128-bit accumulator, then propagate carries at the end.

            // TODO: It may be faster to do this in one pass. We can minimize carry propagation by
            // enumerating the intermediate products in a least-significant-to-most-significant manner.

            let mut acc128 = [0u128; $a_len + $b_len];

            for i in 0..$a_len {
                for j in 0..$b_len {
                    let a_i_b_j = a[i] as u128 * b[j] as u128;
                    // Add the less significant chunk to the less significant accumulator.
                    acc128[i + j] += a_i_b_j as u64 as u128;
                    // Add the more significant chunk to the more significant accumulator.
                    acc128[i + j + 1] += a_i_b_j >> 64;
                }
            }

            let mut acc = [0u64; $a_len + $b_len];
            acc[0] = acc128[0] as u64;
            let mut carry = false;
            for i in 1..($a_len + $b_len) {
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
    }
}

add_symmetric!(add_4_4, 4);
add_symmetric!(add_6_6, 6);

sub_symmetric!(sub_4_4, 4);
sub_symmetric!(sub_6_6, 6);
sub_symmetric!(sub_8_8, 8);
sub_symmetric!(sub_12_12, 12);
sub_asymmetric!(sub_7_6, 7, 6);
sub_asymmetric!(sub_5_4, 5, 4);

mul_symmetric!(mul_4_4, 4);
mul_symmetric!(mul_6_6, 6);
mul_asymmetric!(mul_8_4, 8, 4);
mul_asymmetric!(mul_12_6, 12, 6);

cmp_symmetric!(cmp_4_4, 4);
cmp_symmetric!(cmp_6_6, 6);
cmp_asymmetric!(cmp_5_4, 5, 4);
cmp_asymmetric!(cmp_7_6, 7, 6);

barrett_reduction!(barrett_reduction_bls12_scalar, 4, mul_4_4, mul_8_4, cmp_5_4, sub_8_8, sub_5_4,
Bls12Scalar::ORDER, Bls12Scalar::BARRETT_R, Bls12Scalar::BARRETT_K);
barrett_reduction!(barrett_reduction_bls12_base, 6, mul_6_6, mul_12_6, cmp_7_6, sub_12_12, sub_7_6,
Bls12Base::ORDER, Bls12Base::BARRETT_R, Bls12Base::BARRETT_K);

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use num::BigUint;

    use crate::field::{Bls12Base, mul_12_6, mul_6_6};
    use crate::u64_slice_to_biguint;

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

        assert_eq!(
            u64_slice_to_biguint(&(a_blsbase * b_blsbase).limbs),
            a_biguint * b_biguint % order_biguint);
    }
}
