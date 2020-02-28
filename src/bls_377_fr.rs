//! This module implements field arithmetic for BLS12-377's scalar field.

use std::cmp::Ordering;
use std::cmp::Ordering::Less;
use std::collections::HashSet;
use std::ops::{Add, Div, Mul, Neg, Sub};

use rand::RngCore;
use rand::rngs::OsRng;
use unroll::unroll_for_loops;

/// An element of the BLS12 group's scalar field.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct Bls12Scalar {
    /// Montgomery representation, encoded with little-endian u64 limbs.
    pub limbs: [u64; 4],
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

    /// t = (r - 1) / 2^s =
    /// 60001509534603559531609739528203892656505753216962260608619555
    const T: Self = Self { limbs: [725501752471715841, 6461107452199829505, 6968279316240510977, 1345280370688042326] };

    pub const TWO_ADICITY: usize = 47;

    /// Generator of [1, order).
    const GENERATOR: Bls12Scalar = Bls12Scalar {
        limbs: [
            1855201571499933546, 8511318076631809892, 6222514765367795509, 1122129207579058019]
    };

    pub fn from_canonical(c: [u64; 4]) -> Self {
        // We compute M(c, R^2) = c * R^2 * R^-1 = c * R.
        Self { limbs: Self::montgomery_multiply(c, Self::R2) }
    }

    pub fn from_canonical_u64(c: u64) -> Self {
        Self::from_canonical([c, 0, 0, 0])
    }

    pub fn from_canonical_usize(c: usize) -> Self {
        Self::from_canonical_u64(c as u64)
    }

    pub fn to_canonical(&self) -> [u64; 4] {
        // Let x * R = self. We compute M(x * R, 1) = x * R * R^-1 = x.
        Self::montgomery_multiply(self.limbs, [1, 0, 0, 0])
    }

    pub fn is_zero(&self) -> bool {
        *self == Self::ZERO
    }

    /// Computes a `2^n_power`th primitive root of unity.
    pub fn primitive_root_of_unity(n_power: usize) -> Bls12Scalar {
        assert!(n_power <= Self::TWO_ADICITY);
        let base_root = Self::GENERATOR.exp(Self::T);
        base_root.exp(Self::from_canonical_u64(1u64 << Self::TWO_ADICITY as u64 - n_power as u64))
    }

    pub fn cyclic_subgroup_unknown_order(generator: Bls12Scalar) -> Vec<Bls12Scalar> {
        let mut subgroup_vec = Vec::new();
        let mut subgroup_set = HashSet::new();
        let mut current = Bls12Scalar::ONE;
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
        let mut current = Bls12Scalar::ONE;
        for _i in 0..order {
            subgroup.push(current);
            current = current * generator;
        }
        subgroup
    }

    pub fn generator_order(generator: Bls12Scalar) -> usize {
        Self::cyclic_subgroup_unknown_order(generator).len()
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

    pub fn exp(&self, power: Bls12Scalar) -> Bls12Scalar {
        let power_bits = power.num_bits();
        let mut current = *self;
        let mut product = Bls12Scalar::ONE;

        for (i, limb) in power.to_canonical().iter().enumerate() {
            for j in 0..64 {
                // If we've gone through all the 1 bits already, no need to keep squaring.
                let bit_index = i * 64 + j;
                if bit_index == power_bits {
                    return product;
                }

                if (limb >> j & 1) != 0 {
                    product = product * current;
                }
                current = current.square();
            }
        }

        product
    }

    pub fn exp_usize(&self, power: usize) -> Bls12Scalar {
        self.exp(Self::from_canonical_usize(power))
    }

    pub fn square(&self) -> Self {
        // TODO: Some intermediate products are the redundant, so this can be made faster.
        *self * *self
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

    fn nonzero_multiplicative_inverse_canonical(a: [u64; 4]) -> [u64; 4] {
        // Based on Algorithm 16 of "Efficient Software-Implementation of Finite Fields with
        // Applications to Cryptography".

        let zero = [0, 0, 0, 0];
        let one = [1, 0, 0, 0];

        let mut u = a;
        let mut v = Self::ORDER;
        let mut b = one;
        let mut c = zero;

        while u != one && v != one {
            while Self::is_even(u) {
                u = Self::div2(u);
                if Self::is_odd(b) {
                    b = add_4_4_no_overflow(b, Self::ORDER);
                }
                b = Self::div2(b);
            }

            while Self::is_even(v) {
                v = Self::div2(v);
                if Self::is_odd(c) {
                    c = add_4_4_no_overflow(c, Self::ORDER);
                }
                c = Self::div2(c);
            }

            if cmp_4_4(u, v) == Less {
                v = sub_4_4(v, u);
                if cmp_4_4(c, b) == Less {
                    c = add_4_4_no_overflow(c, Self::ORDER);
                }
                c = sub_4_4(c, b);
            } else {
                u = sub_4_4(u, v);
                if cmp_4_4(b, c) == Less {
                    b = add_4_4_no_overflow(b, Self::ORDER);
                }
                b = sub_4_4(b, c);
            }
        }

        if u == one {
            b
        } else {
            c
        }
    }

    fn is_even(x: [u64; 4]) -> bool {
        x[0] & 1 == 0
    }

    fn is_odd(x: [u64; 4]) -> bool {
        x[0] & 1 == 1
    }

    /// Shift left (in the direction of increasing significance) by 1. Equivalent to integer
    /// division by two.
    #[unroll_for_loops]
    fn div2(x: [u64; 4]) -> [u64; 4] {
        let mut result = [0; 4];
        for i in 0..3 {
            result[i] = x[i] >> 1 | x[i + 1] << 63;
        }
        result[3] = x[3] >> 1;
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
    pub fn rand() -> Bls12Scalar {
        let mut limbs = [0; 4];

        for limb_i in &mut limbs {
            *limb_i = OsRng.next_u64();
        }

        // Remove a few of the most significant bits to ensure we're in range.
        limbs[3] >>= 4;

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
        let sum = add_4_4_no_overflow(self.limbs, rhs.limbs);
        let limbs = if cmp_4_4(sum, Bls12Scalar::ORDER) == Less {
            sum
        } else {
            sub_4_4(sum, Bls12Scalar::ORDER)
        };
        Bls12Scalar { limbs }
    }
}

impl Sub<Bls12Scalar> for Bls12Scalar {
    type Output = Bls12Scalar;

    fn sub(self, rhs: Bls12Scalar) -> Self::Output {
        let limbs = if cmp_4_4(self.limbs, rhs.limbs) == Ordering::Less {
            // Underflow occurs, so we compute the difference as `self + (-rhs)`.
            add_4_4_no_overflow(self.limbs, (-rhs).limbs)
        } else {
            // No underflow, so it's faster to subtract directly.
            sub_4_4(self.limbs, rhs.limbs)
        };
        Self { limbs }
    }
}

impl Mul<Bls12Scalar> for Bls12Scalar {
    type Output = Bls12Scalar;

    fn mul(self, rhs: Bls12Scalar) -> Bls12Scalar {
        Self { limbs: Self::montgomery_multiply(self.limbs, rhs.limbs) }
    }
}

impl Div<Bls12Scalar> for Bls12Scalar {
    type Output = Bls12Scalar;

    fn div(self, rhs: Bls12Scalar) -> Bls12Scalar {
        self * rhs.multiplicative_inverse().expect("No inverse")
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

#[unroll_for_loops]
pub fn cmp_4_4(a: [u64; 4], b: [u64; 4]) -> Ordering {
    for i in (0..4).rev() {
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
fn add_4_4_no_overflow(a: [u64; 4], b: [u64; 4]) -> [u64; 4] {
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

/// Computes `a - b`. Assumes `a >= b`, otherwise the behavior is undefined.
#[unroll_for_loops]
fn sub_4_4(a: [u64; 4], b: [u64; 4]) -> [u64; 4] {
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

#[unroll_for_loops]
pub fn mul_4_4(a: [u64; 4], b: [u64; 4]) -> [u64; 8] {
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

#[cfg(test)]
mod tests {
    use crate::conversions::u64_slice_to_biguint;
    use crate::{Bls12Scalar, mul_6_6};

    #[test]
    fn test_mul_6_6() {
        let a = [11111111u64, 22222222, 33333333, 44444444, 55555555, 66666666];
        let b = [77777777u64, 88888888, 99999999, 11111111, 22222222, 33333333];
        assert_eq!(
            u64_slice_to_biguint(&mul_6_6(a, b)),
            u64_slice_to_biguint(&a) * u64_slice_to_biguint(&b));
    }

    #[test]
    fn bls12scalar_to_and_from_canonical() {
        let a = [1, 2, 3, 4];
        let a_biguint = u64_slice_to_biguint(&a);
        let order_biguint = u64_slice_to_biguint(&Bls12Scalar::ORDER);
        let r_biguint = u64_slice_to_biguint(&Bls12Scalar::R);

        let a_bls12scalar = Bls12Scalar::from_canonical(a);
        assert_eq!(u64_slice_to_biguint(&a_bls12scalar.limbs),
                   &a_biguint * &r_biguint % &order_biguint);
        assert_eq!(u64_slice_to_biguint(&a_bls12scalar.to_canonical()), a_biguint);
    }

    #[test]
    fn mul_bls12_scalar() {
        let a = [1, 2, 3, 4];
        let b = [3, 4, 5, 6];

        let a_biguint = u64_slice_to_biguint(&a);
        let b_biguint = u64_slice_to_biguint(&b);
        let order_biguint = u64_slice_to_biguint(&Bls12Scalar::ORDER);

        let a_blsbase = Bls12Scalar::from_canonical(a);
        let b_blsbase = Bls12Scalar::from_canonical(b);

        assert_eq!(
            u64_slice_to_biguint(&(a_blsbase * b_blsbase).to_canonical()),
            a_biguint * b_biguint % order_biguint);
    }

    #[test]
    fn test_bls12_rand() {
        let random_element = Bls12Scalar::rand();

        for i in 0..4 {
            assert_ne!(random_element.limbs[i], 0x0);
        }
    }

    #[test]
    fn exp() {
        assert_eq!(Bls12Scalar::THREE.exp(Bls12Scalar::ZERO), Bls12Scalar::ONE);
        assert_eq!(Bls12Scalar::THREE.exp(Bls12Scalar::ONE), Bls12Scalar::THREE);
        assert_eq!(Bls12Scalar::THREE.exp(Bls12Scalar::from_canonical_u64(2)), Bls12Scalar::from_canonical_u64(9));
        assert_eq!(Bls12Scalar::THREE.exp(Bls12Scalar::from_canonical_u64(3)), Bls12Scalar::from_canonical_u64(27));
    }

    #[test]
    fn negation() {
        for i in 0..25 {
            let i_blsscalar = Bls12Scalar::from_canonical_u64(i);
            assert_eq!(i_blsscalar + -i_blsscalar, Bls12Scalar::ZERO);
        }
    }

    #[test]
    fn multiplicative_inverse() {
        for i in 0..25 {
            let i_blsscalar = Bls12Scalar::from_canonical_u64(i);
            let i_inv_blsscalar = i_blsscalar.multiplicative_inverse();
            if i == 0 {
                assert!(i_inv_blsscalar.is_none());
            } else {
                assert_eq!(i_blsscalar * i_inv_blsscalar.unwrap(), Bls12Scalar::ONE);
            }
        }
    }

    #[test]
    fn batch_multiplicative_inverse() {
        let mut x = Vec::new();
        for i in 1..25 {
            x.push(Bls12Scalar::from_canonical_u64(i));
        }

        let x_inv = Bls12Scalar::batch_multiplicative_inverse(&x);
        assert_eq!(x.len(), x_inv.len());

        for (x_i, x_i_inv) in x.into_iter().zip(x_inv) {
            assert_eq!(x_i * x_i_inv, Bls12Scalar::ONE);
        }
    }

    #[test]
    fn test_div2() {
        assert_eq!(Bls12Scalar::div2([40, 0, 0, 0]), [20, 0, 0, 0]);

        assert_eq!(
            Bls12Scalar::div2(
                [15668009436471190370, 3102040391300197453, 4166322749169705801, 3518225024268476800]),
            [17057376755090370993, 10774392232504874534, 2083161374584852900, 1759112512134238400]);
    }

    #[test]
    fn num_bits() {
        assert_eq!(Bls12Scalar::from_canonical_u64(0b10101).num_bits(), 5);
        assert_eq!(Bls12Scalar::from_canonical_u64(u64::max_value()).num_bits(), 64);
        assert_eq!(Bls12Scalar::from_canonical([0, 1, 0, 0]).num_bits(), 64 + 1);
        assert_eq!(Bls12Scalar::from_canonical([0, 0, 0, 1]).num_bits(), 64 * 3 + 1);
        assert_eq!(Bls12Scalar::from_canonical([0, 0, 0, 0b10101]).num_bits(), 64 * 3 + 5)
    }

    #[test]
    fn roots_of_unity() {
        for n_power in 0..10 {
            let n = 1 << n_power as u64;
            let root = Bls12Scalar::primitive_root_of_unity(n_power);

            assert_eq!(root.exp(Bls12Scalar::from_canonical_u64(n)), Bls12Scalar::ONE);

            if n > 1 {
                assert_ne!(root.exp(Bls12Scalar::from_canonical_u64(n - 1)), Bls12Scalar::ONE)
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
