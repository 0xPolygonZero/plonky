//! This module implements field arithmetic for BLS12-377's base field.

use std::cmp::Ordering;
use std::cmp::Ordering::Less;
use std::ops::{Add, Div, Mul, Neg, Sub};

use rand::RngCore;
use rand::rngs::OsRng;
use unroll::unroll_for_loops;

/// An element of the BLS12 group's base field.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct Bls12Base {
    /// Montgomery representation, encoded with little-endian u64 limbs.
    pub limbs: [u64; 6],
}

impl Bls12Base {
    pub const ZERO: Self = Self { limbs: [0; 6] };
    pub const ONE: Self = Self {
        limbs: [202099033278250856, 5854854902718660529, 11492539364873682930, 8885205928937022213,
            5545221690922665192, 39800542322357402]
    };
    pub const TWO: Self = Self {
        limbs: [404198066556501712, 11709709805437321058, 4538334656037814244, 17770411857874044427,
            11090443381845330384, 79601084644714804]
    };
    pub const THREE: Self = Self {
        limbs: [606297099834752568, 17564564708155981587, 16030874020911497174, 8208873713101515024,
            16635665072767995577, 119401626967072206]
    };

    pub const BITS: usize = 377;

    /// The order of the field:
    /// 258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177
    pub const ORDER: [u64; 6] = [9586122913090633729, 1660523435060625408, 2230234197602682880,
        1883307231910630287, 14284016967150029115, 121098312706494698];

    pub const ORDER_X2: [u64; 6] = [725501752471715842, 3321046870121250817, 4460468395205365760,
        3766614463821260574, 10121289860590506614, 242196625412989397];

    /// R in the context of the Montgomery reduction, i.e. 2^384 % |F|.
    pub(crate) const R: [u64; 6] = [202099033278250856, 5854854902718660529, 11492539364873682930,
        8885205928937022213, 5545221690922665192, 39800542322357402];

    /// R^2 in the context of the Montgomery reduction, i.e. 2^384^2 % |F|.
    pub(crate) const R2: [u64; 6] = [13224372171368877346, 227991066186625457, 2496666625421784173,
        13825906835078366124, 9475172226622360569, 30958721782860680];

    /// R^3 in the context of the Montgomery reduction, i.e. 2^384^3 % |F|.
    pub(crate) const R3: [u64; 6] = [6349885463227391520, 16505482940020594053, 3163973454937060627,
        7650090842119774734, 4571808961100582073, 73846176275226021];

    /// In the context of Montgomery multiplication, µ = -|F|^-1 mod 2^64.
    const MU: u64 = 9586122913090633727;

    pub fn from_canonical(c: [u64; 6]) -> Self {
        // We compute M(c, R^2) = c * R^2 * R^-1 = c * R.
        Self { limbs: Self::montgomery_multiply(c, Self::R2) }
    }

    pub fn from_canonical_u64(c: u64) -> Self {
        Self::from_canonical([c, 0, 0, 0, 0, 0])
    }

    pub fn to_canonical(&self) -> [u64; 6] {
        // Let x * R = self. We compute M(x * R, 1) = x * R * R^-1 = x.
        Self::montgomery_multiply(self.limbs, [1, 0, 0, 0, 0, 0])
    }

    pub fn is_zero(&self) -> bool {
        *self == Self::ZERO
    }

    pub fn num_bits(&self) -> usize {
        let mut n = 0;
        for (i, limb) in self.to_canonical().iter().enumerate() {
            for j in 0..64 {
                if (limb >> j & 1) != 0 {
                    n = i * 64 + j + 1;
                }
            }
        }
        n
    }

    pub fn double(&self) -> Self {
        let result = Self::mul2(self.limbs);
        let limbs = if cmp_6_6(result, Self::ORDER) == Less {
            result
        } else {
            sub_6_6(result, Self::ORDER)
        };
        Self { limbs }
    }

    pub fn triple(&self) -> Self {
        // First compute (unreduced) self * 3 via a double and add, then reduce. We might need to do
        // two comparisons, but at least we'll always do a single subtraction.

        let mut sum = Self::mul2(self.limbs);
        sum = add_6_6_no_overflow(sum, self.limbs);
        let limbs = if cmp_6_6(sum, Self::ORDER) == Less {
            sum
        } else if cmp_6_6(sum, Self::ORDER_X2) == Less {
            sub_6_6(sum, Self::ORDER)
        } else {
            sub_6_6(sum, Self::ORDER_X2)
        };
        Self { limbs }
    }

    pub fn square(&self) -> Self {
        // TODO: Some intermediate products are the redundant, so this can be made faster.
        *self * *self
    }

    pub fn cube(&self) -> Self {
        self.square() * *self
    }

    pub fn multiplicative_inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            // Let x R = self. We compute M((x R)^-1, R^3) = x^-1 R^-1 R^3 R^-1 = x^-1 R.
            let self_r_inv = Self::nonzero_multiplicative_inverse_canonical(self.limbs);
            Some(Self { limbs: Self::montgomery_multiply(self_r_inv, Self::R3) })
        }
    }

    fn nonzero_multiplicative_inverse_canonical(a: [u64; 6]) -> [u64; 6] {
        // Based on Algorithm 16 of "Efficient Software-Implementation of Finite Fields with
        // Applications to Cryptography".

        let zero = [0, 0, 0, 0, 0, 0];
        let one = [1, 0, 0, 0, 0, 0];

        let mut u = a;
        let mut v = Self::ORDER;
        let mut b = one;
        let mut c = zero;

        while u != one && v != one {
            while Self::is_even(u) {
                u = Self::div2(u);
                if Self::is_odd(b) {
                    b = add_6_6_no_overflow(b, Self::ORDER);
                }
                b = Self::div2(b);
            }

            while Self::is_even(v) {
                v = Self::div2(v);
                if Self::is_odd(c) {
                    c = add_6_6_no_overflow(c, Self::ORDER);
                }
                c = Self::div2(c);
            }

            if cmp_6_6(u, v) == Less {
                v = sub_6_6(v, u);
                if cmp_6_6(c, b) == Less {
                    c = add_6_6_no_overflow(c, Self::ORDER);
                }
                c = sub_6_6(c, b);
            } else {
                u = sub_6_6(u, v);
                if cmp_6_6(b, c) == Less {
                    b = add_6_6_no_overflow(b, Self::ORDER);
                }
                b = sub_6_6(b, c);
            }
        }

        if u == one {
            b
        } else {
            c
        }
    }

    fn is_even(x: [u64; 6]) -> bool {
        x[0] & 1 == 0
    }

    fn is_odd(x: [u64; 6]) -> bool {
        x[0] & 1 == 1
    }

    /// Shift in the direction of increasing significance by 1. Equivalent to integer
    /// division by two.
    #[unroll_for_loops]
    fn div2(x: [u64; 6]) -> [u64; 6] {
        let mut result = [0; 6];
        for i in 0..5 {
            result[i] = x[i] >> 1 | x[i + 1] << 63;
        }
        result[5] = x[5] >> 1;
        result
    }

    /// Shift in the direction of increasing significance by 1. Equivalent to integer
    /// multiplication by two. Assumes that the input's most significant bit is clear.
    #[unroll_for_loops]
    fn mul2(x: [u64; 6]) -> [u64; 6] {
        debug_assert_eq!(x[5] >> 63, 0, "Most significant bit should be clear");

        let mut result = [0; 6];
        result[0] = x[0] << 1;
        for i in 1..6 {
            result[i] = x[i] << 1 | x[i - 1] >> 63;
        }
        result
    }

    pub fn batch_multiplicative_inverse_opt(x: &[Self]) -> Vec<Option<Self>> {
        let n = x.len();
        let mut x_nonzero = Vec::with_capacity(n);
        let mut index_map = Vec::with_capacity(n);

        for x_i in x {
            if !x_i.is_zero() {
                index_map.push(x_nonzero.len());
                x_nonzero.push(*x_i);
            } else {
                // Push an intentionally out-of-bounds index to make sure it isn't used.
                index_map.push(n);
            }
        }

        let x_inv = Self::batch_multiplicative_inverse(&x_nonzero);

        let mut result = Vec::with_capacity(n);
        for i in 0..n {
            if !x[i].is_zero() {
                result.push(Some(x_inv[index_map[i]]));
            } else {
                result.push(None);
            }
        }
        result
    }

    pub fn batch_multiplicative_inverse(x: &[Self]) -> Vec<Self> {
        // This is Montgomery's trick. At a high level, we invert the product of the given field
        // elements, then derive the individual inverses from that via multiplication.

        let n = x.len();
        if n == 0 {
            return Vec::new();
        }

        let mut a = Vec::with_capacity(n);
        a.push(x[0]);
        for i in 1..n {
            a.push(a[i - 1] * x[i]);
        }

        let mut a_inv = vec![Self::ZERO; n];
        a_inv[n - 1] = a[n - 1].multiplicative_inverse().expect("No inverse");
        for i in (0..n - 1).rev() {
            a_inv[i] = x[i + 1] * a_inv[i + 1];
        }

        let mut x_inv = Vec::with_capacity(n);
        x_inv.push(a_inv[0]);
        for i in 1..n {
            x_inv.push(a[i - 1] * a_inv[i]);
        }
        x_inv
    }

    // TODO: replace with a CSPRNG
    pub fn rand() -> Self {
        let mut limbs = [0; 6];

        for limb_i in &mut limbs {
            *limb_i = OsRng.next_u64();
        }

        // Remove a few of the most significant bits to ensure we're in range.
        limbs[5] >>= 4;

        Self { limbs }
    }

    #[unroll_for_loops]
    fn montgomery_multiply(a: [u64; 6], b: [u64; 6]) -> [u64; 6] {
        // Interleaved Montgomery multiplication, as described in Algorithm 2 of
        // https://eprint.iacr.org/2017/1057.pdf

        // Note that in the loop below, to avoid explicitly shifting c, we will treat i as the least
        // significant digit and wrap around.
        let mut c = [0u64; 7];

        for i in 0..6 {
            // Add a[i] b to c.
            let mut carry = 0;
            for j in 0..6 {
                let result = c[(i + j) % 7] as u128 + a[i] as u128 * b[j] as u128 + carry as u128;
                c[(i + j) % 7] = result as u64;
                carry = (result >> 64) as u64;
            }
            c[(i + 6) % 7] += carry;

            // q = u c mod r = u c[0] mod r.
            let q = Bls12Base::MU.wrapping_mul(c[i]);

            // C += N q
            carry = 0;
            for j in 0..6 {
                let result = c[(i + j) % 7] as u128 + q as u128 * Self::ORDER[j] as u128 + carry as u128;
                c[(i + j) % 7] = result as u64;
                carry = (result >> 64) as u64;
            }
            c[(i + 6) % 7] += carry;

            debug_assert_eq!(c[i], 0);
        }

        let mut result = [c[6], c[0], c[1], c[2], c[3], c[4]];
        // Final conditional subtraction.
        if cmp_6_6(result, Self::ORDER) != Ordering::Less {
            result = sub_6_6(result, Self::ORDER);
        }
        result
    }
}

impl Add<Bls12Base> for Bls12Base {
    type Output = Bls12Base;

    fn add(self, rhs: Bls12Base) -> Bls12Base {
        // First we do a widening addition, then we reduce if necessary.
        let sum = add_6_6_no_overflow(self.limbs, rhs.limbs);
        let limbs = if cmp_6_6(sum, Bls12Base::ORDER) == Less {
            sum
        } else {
            sub_6_6(sum, Bls12Base::ORDER)
        };
        Bls12Base { limbs }
    }
}

impl Sub<Bls12Base> for Bls12Base {
    type Output = Bls12Base;

    fn sub(self, rhs: Bls12Base) -> Self::Output {
        let limbs = if cmp_6_6(self.limbs, rhs.limbs) == Ordering::Less {
            // Underflow occurs, so we compute the difference as `self + (-rhs)`.
            add_6_6_no_overflow(self.limbs, (-rhs).limbs)
        } else {
            // No underflow, so it's faster to subtract directly.
            sub_6_6(self.limbs, rhs.limbs)
        };
        Self { limbs }
    }
}

impl Mul<Bls12Base> for Bls12Base {
    type Output = Bls12Base;

    fn mul(self, rhs: Bls12Base) -> Bls12Base {
        Self { limbs: Self::montgomery_multiply(self.limbs, rhs.limbs) }
    }
}

impl Mul<u64> for Bls12Base {
    type Output = Bls12Base;

    fn mul(self, rhs: u64) -> Self::Output {
        // TODO: Do a 6x1 multiplication instead of padding to 6x6.
        let rhs_field = Bls12Base::from_canonical_u64(rhs);
        self * rhs_field
    }
}

impl Div<Bls12Base> for Bls12Base {
    type Output = Bls12Base;

    fn div(self, rhs: Bls12Base) -> Bls12Base {
        self * rhs.multiplicative_inverse().expect("No inverse")
    }
}

impl Neg for Bls12Base {
    type Output = Bls12Base;

    fn neg(self) -> Bls12Base {
        if self == Bls12Base::ZERO {
            Bls12Base::ZERO
        } else {
            Bls12Base { limbs: sub_6_6(Bls12Base::ORDER, self.limbs) }
        }
    }
}

#[unroll_for_loops]
pub fn cmp_6_6(a: [u64; 6], b: [u64; 6]) -> Ordering {
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

/// Computes `a + b`. Assumes that there is no overflow; this is verified only in debug builds.
#[unroll_for_loops]
fn add_6_6_no_overflow(a: [u64; 6], b: [u64; 6]) -> [u64; 6] {
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

#[cfg(test)]
mod tests {
    use crate::conversions::u64_slice_to_biguint;
    use crate::{Bls12Base, mul_6_6};

    #[test]
    fn test_mul_6_6() {
        let a = [11111111u64, 22222222, 33333333, 44444444, 55555555, 66666666];
        let b = [77777777u64, 88888888, 99999999, 11111111, 22222222, 33333333];
        assert_eq!(
            u64_slice_to_biguint(&mul_6_6(a, b)),
            u64_slice_to_biguint(&a) * u64_slice_to_biguint(&b));
    }

    #[test]
    fn bls12base_to_and_from_canonical() {
        let a = [1, 2, 3, 4, 0, 0];
        let a_biguint = u64_slice_to_biguint(&a);
        let order_biguint = u64_slice_to_biguint(&Bls12Base::ORDER);
        let r_biguint = u64_slice_to_biguint(&Bls12Base::R);

        let a_bls12base = Bls12Base::from_canonical(a);
        assert_eq!(u64_slice_to_biguint(&a_bls12base.limbs),
                   &a_biguint * &r_biguint % &order_biguint);
        assert_eq!(u64_slice_to_biguint(&a_bls12base.to_canonical()), a_biguint);
    }

    #[test]
    fn mul_bls12_base() {
        let a = [1, 2, 3, 4, 0, 0];
        let b = [3, 4, 5, 6, 0, 0];

        let a_biguint = u64_slice_to_biguint(&a);
        let b_biguint = u64_slice_to_biguint(&b);
        let order_biguint = u64_slice_to_biguint(&Bls12Base::ORDER);

        let a_blsbase = Bls12Base::from_canonical(a);
        let b_blsbase = Bls12Base::from_canonical(b);

        assert_eq!(
            u64_slice_to_biguint(&(a_blsbase * b_blsbase).to_canonical()),
            a_biguint * b_biguint % order_biguint);
    }

    #[test]
    fn test_bls12_rand() {
        let random_element = Bls12Base::rand();

        for i in 0..4 {
            assert_ne!(random_element.limbs[i], 0x0);
        }
    }

    #[test]
    fn negation() {
        for i in 0..25 {
            let i_blsbase = Bls12Base::from_canonical_u64(i);
            assert_eq!(i_blsbase + -i_blsbase, Bls12Base::ZERO);
        }
    }

    #[test]
    fn multiplicative_inverse() {
        for i in 0..25 {
            let i_blsbase = Bls12Base::from_canonical_u64(i);
            let i_inv_blsbase = i_blsbase.multiplicative_inverse();
            if i == 0 {
                assert!(i_inv_blsbase.is_none());
            } else {
                assert_eq!(i_blsbase * i_inv_blsbase.unwrap(), Bls12Base::ONE);
            }
        }
    }

    #[test]
    fn batch_multiplicative_inverse() {
        let mut x = Vec::new();
        for i in 1..25 {
            x.push(Bls12Base::from_canonical_u64(i));
        }

        let x_inv = Bls12Base::batch_multiplicative_inverse(&x);
        assert_eq!(x.len(), x_inv.len());

        for (x_i, x_i_inv) in x.into_iter().zip(x_inv) {
            assert_eq!(x_i * x_i_inv, Bls12Base::ONE);
        }
    }

    #[test]
    fn test_div2() {
        assert_eq!(Bls12Base::div2([40, 0, 0, 0, 0, 0]), [20, 0, 0, 0, 0, 0]);

        assert_eq!(
            Bls12Base::div2(
                [15668009436471190370, 3102040391300197453, 4166322749169705801, 3518225024268476800, 11231577158546850254, 226224965816356276]),
            [17057376755090370993, 10774392232504874534, 2083161374584852900, 1759112512134238400, 5615788579273425127, 113112482908178138]);
    }

    #[test]
    fn num_bits() {
        assert_eq!(Bls12Base::from_canonical_u64(0b10101).num_bits(), 5);
        assert_eq!(Bls12Base::from_canonical_u64(u64::max_value()).num_bits(), 64);
        assert_eq!(Bls12Base::from_canonical([0, 1, 0, 0, 0, 0]).num_bits(), 64 + 1);
        assert_eq!(Bls12Base::from_canonical([0, 0, 0, 0, 0, 1]).num_bits(), 64 * 5 + 1);
        assert_eq!(Bls12Base::from_canonical([0, 0, 0, 0, 0, 0b10101]).num_bits(), 64 * 5 + 5)
    }
}
