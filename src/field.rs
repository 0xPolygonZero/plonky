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
use std::collections::HashSet;
use std::convert::TryInto;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::str::FromStr;

use num::{BigUint, FromPrimitive};
use num::integer::Integer;
use num::traits::cast::ToPrimitive;
use num::traits::identities::One;
use rand::RngCore;
use rand::rngs::OsRng;
use unroll::unroll_for_loops;

use crate::conversions::{biguint_to_u64_vec, u64_slice_to_biguint};
use crate::num_util::modinv;

/// An element of the BLS12 group's base field.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct Bls12Base {
    /// Montgomery representation, encoded with little-endian u64 limbs.
    pub limbs: [u64; 6],
}

/// An element of the BLS12 group's scalar field.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct Bls12Scalar {
    /// Montgomery representation, encoded with little-endian u64 limbs.
    pub limbs: [u64; 4],
}

impl Bls12Base {
    pub const ZERO: Self = Self { limbs: [0; 6] };
    pub const ONE: Self = Self { limbs: [1, 0, 0, 0, 0, 0] };

    pub const BITS: usize = 377;

    /// The order of the field:
    /// 258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177
    pub const ORDER: [u64; 6] = [9586122913090633729, 1660523435060625408, 2230234197602682880,
        1883307231910630287, 14284016967150029115, 121098312706494698];

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

    fn from_canonical(c: [u64; 6]) -> Self {
        // We compute M(c, R^2) = c * R^2 * R^-1 = c * R.
        Self { limbs: Self::montgomery_multiply(c, Self::R2) }
    }

    fn from_canonical_u64(c: u64) -> Self {
        Self::from_canonical([c, 0, 0, 0, 0, 0])
    }

    fn to_canonical(&self) -> [u64; 6] {
        // Let x * R = self. We compute M(x * R, 1) = x * R * R^-1 = x.
        Self::montgomery_multiply(self.limbs, [1, 0, 0, 0, 0, 0])
    }

    pub fn double(&self) -> Self {
        // TODO: Shift instead of adding.
        *self + *self
    }

    pub fn square(&self) -> Self {
        // TODO: Some intermediate products are the redundant, so this can be made faster.
        *self * *self
    }

    pub fn multiplicative_inverse(&self) -> Option<Self> {
        // TODO: This assumes canonical representation. Need to M-multiply by R^3.
        let self_biguint = u64_slice_to_biguint(&self.limbs);
        let order_biguint = u64_slice_to_biguint(&Self::ORDER);
        let opt_inverse_biguint = modinv(self_biguint, order_biguint);

        opt_inverse_biguint.map(|inverse_biguint| Self {
            limbs: biguint_to_u64_vec(inverse_biguint, 6).as_slice().try_into().unwrap()
        })
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

impl Bls12Scalar {
    pub const ZERO: Self = Self { limbs: [0; 4] };
    pub const ONE: Self = Self { limbs: [9015221291577245683, 8239323489949974514, 1646089257421115374, 958099254763297437] };
    pub const TWO: Self = Self { limbs: [17304940830682775525, 10017539527700119523, 14770643272311271387, 570918138838421475] };
    pub const THREE: Self = Self { limbs: [7147916296078753751, 11795755565450264533, 9448453213491875784, 183737022913545514] };

    pub const BITS: usize = 253;

    /// The order of the field:
    /// 8444461749428370424248824938781546531375899335154063827935233455917409239041
    pub const ORDER: [u64; 4] = [725501752471715841, 6461107452199829505, 6968279316240510977, 1345280370688173398];

    /// R in the context of the Montgomery reduction, i.e. 2^256 % |F|.
    pub(crate) const R: [u64; 4] =
        [9015221291577245683, 8239323489949974514, 1646089257421115374, 958099254763297437];

    /// R^2 in the context of the Montgomery reduction, i.e. 2^256^2 % |F|.
    pub(crate) const R2: [u64; 4] =
        [2726216793283724667, 14712177743343147295, 12091039717619697043, 81024008013859129];

    /// R^3 in the context of the Montgomery reduction, i.e. 2^256^3 % |F|.
    pub(crate) const R3: [u64; 4] =
        [7656847007262524748, 7083357369969088153, 12818756329091487507, 432872940405820890];

    /// In the context of Montgomery multiplication, µ = -|F|^-1 mod 2^64.
    const MU: u64 = 725501752471715839;

    pub const TWO_ADICITY: usize = 47;

    /// Generator of [1, order).
    const GENERATOR: Bls12Scalar = Bls12Scalar { limbs: [
        1855201571499933546, 8511318076631809892, 6222514765367795509, 1122129207579058019] };

    fn from_canonical(c: [u64; 4]) -> Self {
        // We compute M(c, R^2) = c * R^2 * R^-1 = c * R.
        Self { limbs: Self::montgomery_multiply(c, Self::R2) }
    }

    fn from_canonical_u64(c: u64) -> Self {
        Self::from_canonical([c, 0, 0, 0])
    }

    fn to_canonical(&self) -> [u64; 4] {
        // Let x * R = self. We compute M(x * R, 1) = x * R * R^-1 = x.
        Self::montgomery_multiply(self.limbs, [1, 0, 0, 0])
    }

    /// Computes a `2^n_power`th primitive root of unity.
    pub fn primitive_root_of_unity(n_power: usize) -> Bls12Scalar {
        assert!(n_power <= 47);
        let t = BigUint::from_str("60001509534603559531609739528203892656505753216962260608619555").unwrap();
        let base_root = Self::GENERATOR.exp(t);
        base_root.exp(BigUint::from_u64(1u64 << 47u64 - n_power as u64).unwrap())
    }

    pub fn cyclic_subgroup_unknown_order(generator: Bls12Scalar) -> Vec<Bls12Scalar> {
        let mut subgroup_vec = Vec::new();
        let mut subgroup_set = HashSet::new();
        let mut current = generator;
        loop {
            if !subgroup_set.insert(current) {
                break;
            }
            subgroup_vec.push(current);
            current = current * generator;
        }
        subgroup_vec
    }

    pub fn cyclic_subgroup_known_order(generator: Bls12Scalar, order: usize) -> Vec<Bls12Scalar> {
        let mut subgroup = Vec::new();
        let mut current = generator;
        for i in 0..order {
            subgroup.push(current);
            current = current * generator;
        }
        subgroup
    }

    pub fn generator_order(generator: Bls12Scalar) -> usize {
        Self::cyclic_subgroup_unknown_order(generator).len()
    }

    pub fn exp(&self, power: BigUint) -> Bls12Scalar {
        let mut current = *self;
        let mut product = Bls12Scalar::ONE;
        for byte in power.to_bytes_le() {
            for i in 0..8 {
                if (byte >> i & 1) != 0 {
                    product = product * current;
                }
                current = current.square();
            }
        }
        product
    }

    pub fn square(&self) -> Self {
        // TODO: Some intermediate products are the redundant, so this can be made faster.
        *self * *self
    }

    pub fn multiplicative_inverse(&self) -> Option<Self> {
        // Let x R = self. We compute M((x R)^-1, R^3) = x^-1 R^-1 R^3 R^-1 = x^-1 R.
        // We use BigUints for now, since we don't care much about inverse performance.
        let self_biguint = u64_slice_to_biguint(&self.limbs);
        let order_biguint = u64_slice_to_biguint(&Self::ORDER);
        let opt_inverse_biguint = modinv(self_biguint, order_biguint);

        opt_inverse_biguint.map(|inverse_biguint| Self {
            limbs: Self::montgomery_multiply(
                biguint_to_u64_vec(inverse_biguint, 4).as_slice().try_into().unwrap(),
                Bls12Scalar::R3
            )
        })
    }

    // TODO: replace with a CSPRNG
    pub fn rand() -> Bls12Scalar {
        let mut limbs = [0; 4];

        for i in 0..4 {
            limbs[i] = OsRng.next_u64();
        }

        Bls12Scalar { limbs }
    }

    #[unroll_for_loops]
    fn montgomery_multiply(a: [u64; 4], b: [u64; 4]) -> [u64; 4] {
        // Interleaved Montgomery multiplication, as described in Algorithm 2 of
        // https://eprint.iacr.org/2017/1057.pdf

        // Note that in the loop below, to avoid explicitly shifting c, we will treat i as the least
        // significant digit and wrap around.
        let mut c = [0u64; 5];

        for i in 0..4 {
            // Add a[i] b to c.
            let mut carry = 0;
            for j in 0..4 {
                let result = c[(i + j) % 5] as u128 + a[i] as u128 * b[j] as u128 + carry as u128;
                c[(i + j) % 5] = result as u64;
                carry = (result >> 64) as u64;
            }
            c[(i + 4) % 5] += carry;

            // q = u c mod r = u c[0] mod r.
            let q = Bls12Scalar::MU.wrapping_mul(c[i]);

            // C += N q
            carry = 0;
            for j in 0..4 {
                let result = c[(i + j) % 5] as u128 + q as u128 * Self::ORDER[j] as u128 + carry as u128;
                c[(i + j) % 5] = result as u64;
                carry = (result >> 64) as u64;
            }
            c[(i + 4) % 5] += carry;

            debug_assert_eq!(c[i], 0);
        }

        let mut result = [c[4], c[0], c[1], c[2]];
        // Final conditional subtraction.
        if cmp_4_4(result, Self::ORDER) != Ordering::Less {
            result = sub_4_4(result, Self::ORDER);
        }
        result
    }
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
        Self { limbs: Self::montgomery_multiply(self.limbs, rhs.limbs) }
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
        let rhs_field = Bls12Base { limbs: [rhs, 0x0, 0x0, 0x0, 0x0, 0x0] };
        self * rhs_field
    }
}

impl Div<Bls12Base> for Bls12Base {
    type Output = Bls12Base;

    fn div(self, rhs: Bls12Base) -> Bls12Base {
        self * rhs.multiplicative_inverse().expect("No inverse")
    }
}

impl Div<Bls12Scalar> for Bls12Scalar {
    type Output = Bls12Scalar;

    fn div(self, rhs: Bls12Scalar) -> Bls12Scalar {
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

impl Neg for Bls12Scalar {
    type Output = Bls12Scalar;

    fn neg(self) -> Bls12Scalar {
        if self == Bls12Scalar::ZERO {
            Bls12Scalar::ZERO
        } else {
            Bls12Scalar { limbs: sub_4_4(Bls12Scalar::ORDER, self.limbs) }
        }
    }
}

macro_rules! shl {
    ($name:ident, $in_len:expr, $out_len:expr) => {
        /// Shift each bit `n` bits to the left, i.e., in the direction of decreasing significance.
        #[unroll_for_loops]
        fn $name(a: [u64; $in_len], n: usize) -> [u64; $out_len] {
            let mut result = [0u64; $out_len];
            for i in 0..$out_len {
                let shift_words = n / 64;
                let shift_bits = n as u64 % 64;
                result[i] = a[i + shift_words] >> shift_bits
                    | a[i + shift_words + 1] << (64 - shift_bits);
            }
            result
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
        #[unroll_for_loops]
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
        #[unroll_for_loops]
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
        #[unroll_for_loops]
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
        #[unroll_for_loops]
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
        #[unroll_for_loops]
        pub fn $name(a: [u64; $a_len], b: [u64; $b_len]) -> [u64; $a_len + $b_len] {
            // Grade school multiplication. To avoid carrying at each of O(n^2) steps, we first add each
            // intermediate product to a 128-bit accumulator, then propagate carries at the end.
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

macro_rules! div_asymmetric {
    ($name:ident, $a_len:expr, $b_len:expr) => {
        /// Integer division. Returns (quotient, remainder).
        pub fn $name(a: [u64; $a_len], b: [u64; $b_len]) -> ([u64; $a_len], [u64; $a_len]) {
            // For now, we're not too interested in optimizing division speed, so we just use num's
            // implementation.
            let a_biguint = u64_slice_to_biguint(&a);
            let b_biguint = u64_slice_to_biguint(&b);
            let (q_biguint, r_biguint) = a_biguint.div_rem(&b_biguint);
            let q = biguint_to_u64_vec(q_biguint, $a_len).as_slice().try_into().unwrap();
            let r = biguint_to_u64_vec(r_biguint, $a_len).as_slice().try_into().unwrap();
            (q, r)
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

div_asymmetric!(div_4_4, 4, 4);
div_asymmetric!(div_6_6, 6, 6);
div_asymmetric!(div_8_4, 8, 4);
div_asymmetric!(div_12_6, 12, 6);

cmp_symmetric!(cmp_4_4, 4);
cmp_symmetric!(cmp_6_6, 6);
cmp_asymmetric!(cmp_5_4, 5, 4);
cmp_asymmetric!(cmp_7_6, 7, 6);

shl!(shl_12_4, 12, 4);
shl!(shl_18_6, 18, 6);

#[cfg(test)]
mod tests {
    use num::{BigUint, FromPrimitive, One, Zero};

    use crate::conversions::u64_slice_to_biguint;
    use crate::field::{Bls12Base, Bls12Scalar, mul_12_6, mul_6_6};

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
    fn test_to_and_from_canonical() {
        let a = [1, 2, 3, 4, 0, 0];
        let a_biguint = u64_slice_to_biguint(&a);
        let order_biguint = u64_slice_to_biguint(&Bls12Base::ORDER);
        let r_biguint = u64_slice_to_biguint(&Bls12Base::R);

        let a_blsbase = Bls12Base::from_canonical(a);
        assert_eq!(u64_slice_to_biguint(&a_blsbase.limbs),
                   &a_biguint * &r_biguint % &order_biguint);
        assert_eq!(u64_slice_to_biguint(&a_blsbase.to_canonical()), a_biguint);
    }

    #[test]
    fn test_mul_bls12_base() {
        let a = [1, 2, 3, 4, 0, 0];
        let b = [3, 4, 5, 6, 0, 0];

        let a_biguint = u64_slice_to_biguint(&a);
        let b_biguint = u64_slice_to_biguint(&b);
        let order_biguint = u64_slice_to_biguint(&Bls12Base::ORDER);
        let r_biguint = u64_slice_to_biguint(&Bls12Base::R);

        let a_blsbase = Bls12Base::from_canonical(a);
        let b_blsbase = Bls12Base::from_canonical(b);

        assert_eq!(
            u64_slice_to_biguint(&(a_blsbase * b_blsbase).to_canonical()),
            a_biguint * b_biguint % order_biguint);
    }

    #[test]
    fn test_bls12_rand() {
        let random_element = Bls12Scalar::rand();

        for i in 0..4 {
            assert!(random_element.limbs[i] != 0x0);
        }
    }

    #[test]
    fn exp() {
        assert_eq!(Bls12Scalar::THREE.exp(BigUint::zero()), Bls12Scalar::ONE);
        assert_eq!(Bls12Scalar::THREE.exp(BigUint::one()), Bls12Scalar::THREE);
        assert_eq!(Bls12Scalar::THREE.exp(BigUint::from_u8(2).unwrap()), Bls12Scalar::from_canonical_u64(9));
        assert_eq!(Bls12Scalar::THREE.exp(BigUint::from_u8(3).unwrap()), Bls12Scalar::from_canonical_u64(27));
    }

    #[test]
    fn roots_of_unity() {
        for n_power in 0..10 {
            let n = 1 << n_power as u64;
            let root = Bls12Scalar::primitive_root_of_unity(n_power);

            assert_eq!(root.exp(BigUint::from_u64(n).unwrap()), Bls12Scalar::ONE);

            if n > 1 {
                assert_ne!(root.exp(BigUint::from_u64(n - 1).unwrap()), Bls12Scalar::ONE)
            }
        }
    }

    #[test]
    fn primitive_root_order() {
        for n_power in 0..10 {
            let root = Bls12Scalar::primitive_root_of_unity(n_power);
            let order = Bls12Scalar::generator_order(root);
            assert_eq!(order, 1 << n_power, "2^{}'th primitive root", n_power);
        }
    }
}
