//! This module implements field arithmetic for BLS12-377's scalar field.

use std::cmp::Ordering::Less;
use std::convert::TryInto;
use std::ops::{Add, Div, Mul, Neg, Sub};
use rand::Rng;

use unroll::unroll_for_loops;

use crate::{add_4_4_no_overflow, cmp_4_4, Field, sub_4_4, field_to_biguint, rand_range_4, rand_range_4_from_rng};
use crate::nonzero_multiplicative_inverse_4;
use std::cmp::Ordering;
use std::fmt;
use std::fmt::{Display, Formatter};

/// An element of the BLS12 group's scalar field.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub struct Bls12377Scalar {
    /// Montgomery representation, encoded with little-endian u64 limbs.
    pub limbs: [u64; 4],
}

impl Bls12377Scalar {
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

    pub fn from_canonical(c: [u64; 4]) -> Self {
        // We compute M(c, R^2) = c * R^2 * R^-1 = c * R.
        Self { limbs: Self::montgomery_multiply(c, Self::R2) }
    }

    pub fn to_canonical(&self) -> [u64; 4] {
        // Let x * R = self. We compute M(x * R, 1) = x * R * R^-1 = x.
        Self::montgomery_multiply(self.limbs, [1, 0, 0, 0])
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
            let q = Self::MU.wrapping_mul(c[i]);

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
        if cmp_4_4(result, Self::ORDER) != Less {
            result = sub_4_4(result, Self::ORDER);
        }
        result
    }
}

impl Add<Bls12377Scalar> for Bls12377Scalar {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        // First we do a widening addition, then we reduce if necessary.
        let sum = add_4_4_no_overflow(self.limbs, rhs.limbs);
        let limbs = if cmp_4_4(sum, Self::ORDER) == Less {
            sum
        } else {
            sub_4_4(sum, Self::ORDER)
        };
        Self { limbs }
    }
}

impl Sub<Bls12377Scalar> for Bls12377Scalar {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let limbs = if cmp_4_4(self.limbs, rhs.limbs) == Less {
            // Underflow occurs, so we compute the difference as `self + (-rhs)`.
            add_4_4_no_overflow(self.limbs, (-rhs).limbs)
        } else {
            // No underflow, so it's faster to subtract directly.
            sub_4_4(self.limbs, rhs.limbs)
        };
        Self { limbs }
    }
}

impl Mul<Self> for Bls12377Scalar {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self { limbs: Self::montgomery_multiply(self.limbs, rhs.limbs) }
    }
}

impl Div<Bls12377Scalar> for Bls12377Scalar {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        self * rhs.multiplicative_inverse().expect("No inverse")
    }
}

impl Neg for Bls12377Scalar {
    type Output = Self;

    fn neg(self) -> Self {
        if self == Self::ZERO {
            Self::ZERO
        } else {
            Self { limbs: sub_4_4(Self::ORDER, self.limbs) }
        }
    }
}

impl Field for Bls12377Scalar {
    const BITS: usize = 253;
    const BYTES: usize = 32;

    const ZERO: Self = Self { limbs: [0; 4] };
    const ONE: Self = Self { limbs: [9015221291577245683, 8239323489949974514, 1646089257421115374, 958099254763297437] };
    const TWO: Self = Self { limbs: [17304940830682775525, 10017539527700119523, 14770643272311271387, 570918138838421475] };
    const THREE: Self = Self { limbs: [7147916296078753751, 11795755565450264533, 9448453213491875784, 183737022913545514] };
    const FOUR: Self = Self { limbs: [16163137587655999434, 1588334981690687431, 11094542470912991159, 1141836277676842951] };
    const FIVE: Self = Self { limbs: [6006113053051977660, 3366551019440832441, 5772352412093595556, 754655161751966990] };
    const NEG_ONE: Self = Self { limbs: [10157024534604021774, 16668528035959406606, 5322190058819395602, 387181115924875961] };

    const MULTIPLICATIVE_SUBGROUP_GENERATOR: Self = Self { limbs: [1855201571499933546, 8511318076631809892, 6222514765367795509, 1122129207579058019] };

    /// x^11 is a permutation in this field.
    const ALPHA: Self = Self { limbs: [1855201571499933546, 8511318076631809892, 6222514765367795509, 1122129207579058019] };

    const TWO_ADICITY: usize = 47;

    /// 60001509534603559531609739528203892656505753216962260608619555
    const T: Self = Self { limbs: [725501752471715841, 6461107452199829505, 6968279316240510977, 1345280370688042326] };

    fn to_canonical_u64_vec(&self) -> Vec<u64> {
        self.to_canonical().to_vec()
    }

    fn from_canonical_u64_vec(v: Vec<u64>) -> Self {
        Self::from_canonical(v[..].try_into().unwrap())
    }

    fn from_canonical_u64(n: u64) -> Self {
        Self::from_canonical([n, 0, 0, 0])
    }

    fn is_valid_canonical_u64(v: &Vec<u64>) -> bool {
        v.len() == 4 && cmp_4_4(v[..].try_into().unwrap(), Self::ORDER) == Less
    }

    fn multiplicative_inverse_assuming_nonzero(&self) -> Self {
        // Let x R = self. We compute M((x R)^-1, R^3) = x^-1 R^-1 R^3 R^-1 = x^-1 R.
        let self_r_inv = nonzero_multiplicative_inverse_4(self.limbs, Self::ORDER);
        Self { limbs: Self::montgomery_multiply(self_r_inv, Self::R3) }
    }

    fn rand() -> Self {
        Self {
            limbs: rand_range_4(Self::ORDER),
        }
    }

    fn rand_from_rng<R: Rng>(rng: &mut R) -> Self {
        Self {
            limbs: rand_range_4_from_rng(Self::ORDER, rng),
        }
    }


}

impl Ord for Bls12377Scalar {
    fn cmp(&self, other: &Self) -> Ordering {
        self.cmp_helper(other)
    }
}

impl PartialOrd for Bls12377Scalar {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Display for Bls12377Scalar {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        field_to_biguint(*self).fmt(f)
    }
}


#[cfg(test)]
mod tests {
    use crate::{Bls12377Scalar, Field};
    use crate::conversions::u64_slice_to_biguint;
    use crate::test_square_root;

    #[test]
    fn bls12scalar_to_and_from_canonical() {
        let a = [1, 2, 3, 4];
        let a_biguint = u64_slice_to_biguint(&a);
        let order_biguint = u64_slice_to_biguint(&Bls12377Scalar::ORDER);
        let r_biguint = u64_slice_to_biguint(&Bls12377Scalar::R);

        let a_bls12scalar = Bls12377Scalar::from_canonical(a);
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
        let order_biguint = u64_slice_to_biguint(&Bls12377Scalar::ORDER);

        let a_blsbase = Bls12377Scalar::from_canonical(a);
        let b_blsbase = Bls12377Scalar::from_canonical(b);

        assert_eq!(
            u64_slice_to_biguint(&(a_blsbase * b_blsbase).to_canonical()),
            a_biguint * b_biguint % order_biguint);
    }

    #[test]
    fn test_bls12_rand() {
        let random_element = Bls12377Scalar::rand();

        for i in 0..4 {
            assert_ne!(random_element.limbs[i], 0x0);
        }
    }

    #[test]
    fn exp() {
        assert_eq!(Bls12377Scalar::THREE.exp(Bls12377Scalar::ZERO), Bls12377Scalar::ONE);
        assert_eq!(Bls12377Scalar::THREE.exp(Bls12377Scalar::ONE), Bls12377Scalar::THREE);
        assert_eq!(Bls12377Scalar::THREE.exp(Bls12377Scalar::from_canonical_u64(2)), Bls12377Scalar::from_canonical_u64(9));
        assert_eq!(Bls12377Scalar::THREE.exp(Bls12377Scalar::from_canonical_u64(3)), Bls12377Scalar::from_canonical_u64(27));
    }

    #[test]
    fn negation() {
        for i in 0..25 {
            let i_blsscalar = Bls12377Scalar::from_canonical_u64(i);
            assert_eq!(i_blsscalar + -i_blsscalar, Bls12377Scalar::ZERO);
        }
    }

    #[test]
    fn multiplicative_inverse() {
        for i in 0..25 {
            let i_blsscalar = Bls12377Scalar::from_canonical_u64(i);
            let i_inv_blsscalar = i_blsscalar.multiplicative_inverse();
            if i == 0 {
                assert!(i_inv_blsscalar.is_none());
            } else {
                assert_eq!(i_blsscalar * i_inv_blsscalar.unwrap(), Bls12377Scalar::ONE);
            }
        }
    }

    #[test]
    fn batch_multiplicative_inverse() {
        let mut x = Vec::new();
        for i in 1..25 {
            x.push(Bls12377Scalar::from_canonical_u64(i));
        }

        let x_inv = Bls12377Scalar::batch_multiplicative_inverse(&x);
        assert_eq!(x.len(), x_inv.len());

        for (x_i, x_i_inv) in x.into_iter().zip(x_inv) {
            assert_eq!(x_i * x_i_inv, Bls12377Scalar::ONE);
        }
    }

    #[test]
    fn num_bits() {
        assert_eq!(Bls12377Scalar::from_canonical_u64(0b10101).num_bits(), 5);
        assert_eq!(Bls12377Scalar::from_canonical_u64(u64::max_value()).num_bits(), 64);
        assert_eq!(Bls12377Scalar::from_canonical([0, 1, 0, 0]).num_bits(), 64 + 1);
        assert_eq!(Bls12377Scalar::from_canonical([0, 0, 0, 1]).num_bits(), 64 * 3 + 1);
        assert_eq!(Bls12377Scalar::from_canonical([0, 0, 0, 0b10101]).num_bits(), 64 * 3 + 5)
    }

    #[test]
    fn roots_of_unity() {
        for n_power in 0..10 {
            let n = 1 << n_power as u64;
            let root = Bls12377Scalar::primitive_root_of_unity(n_power);

            assert_eq!(root.exp(Bls12377Scalar::from_canonical_u64(n)), Bls12377Scalar::ONE);

            if n > 1 {
                assert_ne!(root.exp(Bls12377Scalar::from_canonical_u64(n - 1)), Bls12377Scalar::ONE)
            }
        }
    }

    #[test]
    fn primitive_root_order() {
        for n_power in 0..10 {
            let root = Bls12377Scalar::primitive_root_of_unity(n_power);
            let order = Bls12377Scalar::generator_order(root);
            assert_eq!(order, 1 << n_power, "2^{}'th primitive root", n_power);
        }
    }

    test_square_root!(Bls12377Scalar);
}
