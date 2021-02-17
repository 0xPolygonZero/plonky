use rand::Rng;
use std::cmp::Ordering;
use std::cmp::Ordering::Less;
use std::convert::TryInto;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::fmt;
use std::fmt::{Debug, Display, Formatter};

use crate::{cmp, field_to_biguint,
            rand_range, rand_range_from_rng,
            MontyRepr, Field};

/// An element of the Pallas group's base field.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Default)]
pub struct PallasBase {
    /// Montgomery representation, encoded with little-endian u64 limbs.
    pub limbs: [u64; 4],
}

impl MontyRepr for PallasBase {
    /// The order of the field: 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
    const ORDER: [u64; 4] = [
        0x992d30ed00000001,
        0x224698fc094cf91b,
        0x0,
        0x4000000000000000,
    ];

    /// Twice the order of the field
    #[allow(dead_code)]
    const ORDER_X2: [u64; 4] = [
        0x325a61da00000002,
        0x448d31f81299f237,
        0x0,
        0x8000000000000000,
    ];

    /// R in the context of the Montgomery reduction, i.e. 2^256 % |F|.
    const R: [u64; 4] = [
        0x34786d38fffffffd,
        0x992c350be41914ad,
        0xffffffffffffffff,
        0x3fffffffffffffff
    ];

    /// R^2 in the context of the Montgomery reduction, i.e. 2^(256*2) % |F|.
    const R2: [u64; 4] = [
        0x8c78ecb30000000f,
        0xd7d30dbd8b0de0e7,
        0x7797a99bc3c95d18,
        0x96d41af7b9cb714
    ];

    /// R^3 in the context of the Montgomery reduction, i.e. 2^(256*3) % |F|.
    const R3: [u64; 4] = [
        0xf185a5993a9e10f9,
        0xf6a68f3b6ac5b1d1,
        0xdf8d1014353fd42c,
        0x2ae309222d2d9910
    ];

    /// In the context of Montgomery multiplication, µ = -|F|^-1 mod 2^64.
    const MU: u64 = 0x992d30ecffffffff;
}

impl PallasBase {
    pub fn from_canonical(c: [u64; 4]) -> Self {
        Self { limbs: Self::from_monty(c) }
    }

    pub fn to_canonical(&self) -> [u64; 4] {
        Self::to_monty(self.limbs)
    }
}

impl Add<PallasBase> for PallasBase {
    type Output = Self;

    fn add(self, rhs: PallasBase) -> Self::Output {
        Self { limbs: Self::monty_add(self.limbs, rhs.limbs) }
    }
}

impl Sub<PallasBase> for PallasBase {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self { limbs: Self::monty_sub(self.limbs, rhs.limbs) }
    }
}

impl Mul<PallasBase> for PallasBase {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self { limbs: Self::monty_multiply(self.limbs, rhs.limbs) }
    }
}

impl Div<PallasBase> for PallasBase {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        self * rhs.multiplicative_inverse().expect("No inverse")
    }
}

impl Neg for PallasBase {
    type Output = Self;

    fn neg(self) -> Self {
        Self { limbs: Self::monty_neg(self.limbs) }
    }
}

impl Field for PallasBase {
    const BITS: usize = 255;
    const BYTES: usize = 32;
    const ZERO: Self = Self { limbs: <Self as MontyRepr>::ZERO };
    const ONE: Self = Self { limbs: <Self as MontyRepr>::ONE };
    const TWO: Self = Self {
        limbs: [
            0xcfc3a984fffffff9,
            0x1011d11bbee5303e,
            0xffffffffffffffff,
            0x3fffffffffffffff
        ],
    };
    const THREE: Self = Self {
        limbs: [
            0x6b0ee5d0fffffff5,
            0x86f76d2b99b14bd0,
            0xfffffffffffffffe,
            0x3fffffffffffffff
        ],
    };
    const FOUR: Self = Self {
        limbs: [
            0x65a221cfffffff1,
            0xfddd093b747d6762,
            0xfffffffffffffffd,
            0x3fffffffffffffff
        ],
    };
    const FIVE: Self = Self {
        limbs: [
            0xa1a55e68ffffffed,
            0x74c2a54b4f4982f3,
            0xfffffffffffffffd,
            0x3fffffffffffffff
        ],
    };
    const NEG_ONE: Self = Self {
        limbs: [0x64b4c3b400000004, 0x891a63f02533e46e, 0, 0]
    };

    const MULTIPLICATIVE_SUBGROUP_GENERATOR: Self = Self::FIVE;

    const ALPHA: Self = Self::FIVE;

    const TWO_ADICITY: usize = 32;

    /// T = (ORDER - 1) / 2^TWO_ADICITY  in Monty form
    const T: Self = Self {
        limbs: [
            0x992d30ed00000001,
            0x224698fc094cf91b,
            0x0,
            0x3fffffff00000000
        ],
    };

    fn to_canonical_u64_vec(&self) -> Vec<u64> {
        self.to_canonical().to_vec()
    }

    fn from_canonical_u64_vec(v: Vec<u64>) -> Self {
        Self::from_canonical(v[..].try_into().unwrap())
    }

    fn from_canonical_u64(n: u64) -> Self {
        Self::from_canonical([n, 0, 0, 0])
    }

    fn is_valid_canonical_u64(v: &[u64]) -> bool {
        v.len() == 4 && cmp(v[..].try_into().unwrap(), Self::ORDER) == Less
    }

    fn multiplicative_inverse_assuming_nonzero(&self) -> Self {
        Self {
            limbs: Self::monty_inverse(self.limbs)
        }
    }

    fn rand() -> Self {
        Self {
            limbs: rand_range(Self::ORDER),
        }
    }

    fn rand_from_rng<R: Rng>(rng: &mut R) -> Self {
        Self {
            limbs: rand_range_from_rng(Self::ORDER, rng),
        }
    }

    #[inline(always)]
    fn square(&self) -> Self {
        Self {
            limbs: Self::monty_square(self.limbs),
        }
    }
}

impl Ord for PallasBase {
    fn cmp(&self, other: &Self) -> Ordering {
        self.cmp_helper(other)
    }
}

impl PartialOrd for PallasBase {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Display for PallasBase {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}", field_to_biguint(*self))
    }
}

impl Debug for PallasBase {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "PallasBase {}", field_to_biguint(*self))
    }
}

#[cfg(test)]
mod tests {
    use crate::test_arithmetic;
    use crate::Field;
    use crate::PallasBase;
    use crate::MontyRepr; // This is just to access ORDER_X2.

    #[test]
    fn primitive_root_order() {
        for n_power in 0..10 {
            let root = PallasBase::primitive_root_of_unity(n_power);
            let order = PallasBase::generator_order(root);
            assert_eq!(order, 1 << n_power, "2^{}'th primitive root", n_power);
        }
    }

    #[test]
    fn valid_canonical_vec() {
        let small = <PallasBase as Field>::ONE.to_canonical_u64_vec();
        assert!(PallasBase::is_valid_canonical_u64(&small));

        let big = PallasBase::ORDER_X2.to_vec();
        assert_eq!(PallasBase::is_valid_canonical_u64(&big), false);

        let limbs = vec![1, 2, 3, 4, 5];
        assert_eq!(PallasBase::is_valid_canonical_u64(&limbs), false);
    }

    test_arithmetic!(crate::PallasBase);
}
