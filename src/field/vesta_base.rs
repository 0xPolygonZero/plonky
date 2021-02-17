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

/// An element of the Vesta group's base field.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Default)]
pub struct VestaBase {
    /// Montgomery representation, encoded with little-endian u64 limbs.
    pub limbs: [u64; 4],
}

impl MontyRepr for VestaBase {
    /// The order of the field: 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001
    const ORDER: [u64; 4] = [
        0x8c46eb2100000001,
        0x224698fc0994a8dd,
        0x0,
        0x4000000000000000
    ];

    /// Twice the order of the field:
    #[allow(dead_code)]
    const ORDER_X2: [u64; 4] = [
        0x188dd64200000002,
        0x448d31f8132951bb,
        0x0,
        0x8000000000000000
    ];

    /// R in the context of the Montgomery reduction, i.e. 2^256 % |F|.
    const R: [u64; 4] = [
        0x5b2b3e9cfffffffd,
        0x992c350be3420567,
        0xffffffffffffffff,
        0x3fffffffffffffff
    ];

    /// R^2 in the context of the Montgomery reduction, i.e. 2^(256*2) % |F|.
    const R2: [u64; 4] = [
        0xfc9678ff0000000f,
        0x67bb433d891a16e3,
        0x7fae231004ccf590,
        0x96d41af7ccfdaa9
    ];

    /// R^3 in the context of the Montgomery reduction, i.e. 2^(256*3) % |F|.
    const R3: [u64; 4] = [
        0x8b421c249dae4c,
        0xe13bda50dba41326,
        0x88fececb8e15cb63,
        0x7dd97a06e6792c8
    ];

    /// In the context of Montgomery multiplication, µ = -|F|^-1 mod 2^64.
    const MU: u64 = 0x8c46eb20ffffffff;
}

impl VestaBase {
    pub fn from_canonical(c: [u64; 4]) -> Self {
        Self { limbs: Self::from_monty(c) }
    }

    pub fn to_canonical(&self) -> [u64; 4] {
        Self::to_monty(self.limbs)
    }
}

impl Add<VestaBase> for VestaBase {
    type Output = Self;

    fn add(self, rhs: VestaBase) -> Self::Output {
        Self { limbs: Self::monty_add(self.limbs, rhs.limbs) }
    }
}

impl Sub<VestaBase> for VestaBase {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self { limbs: Self::monty_sub(self.limbs, rhs.limbs) }
    }
}

impl Mul<VestaBase> for VestaBase {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self { limbs: Self::monty_multiply(self.limbs, rhs.limbs) }
    }
}

impl Div<VestaBase> for VestaBase {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        self * rhs.multiplicative_inverse().expect("No inverse")
    }
}

impl Neg for VestaBase {
    type Output = Self;

    fn neg(self) -> Self {
        Self { limbs: Self::monty_neg(self.limbs) }
    }
}

impl Field for VestaBase {
    const BITS: usize = 255;
    const BYTES: usize = 32;
    const ZERO: Self = Self { limbs: <Self as MontyRepr>::ZERO };
    const ONE: Self = Self { limbs: <Self as MontyRepr>::ONE };
    const TWO: Self = Self {
        limbs: [
            0x2a0f9218fffffff9,
            0x1011d11bbcef61f1,
            0xffffffffffffffff,
            0x3fffffffffffffff
        ],
    };
    const THREE: Self = Self {
        limbs: [
            0xf8f3e594fffffff5,
            0x86f76d2b969cbe7a,
            0xfffffffffffffffe,
            0x3fffffffffffffff
        ],
    };
    const FOUR: Self = Self {
        limbs: [
            0xc7d83910fffffff1,
            0xfddd093b704a1b04,
            0xfffffffffffffffd,
            0x3fffffffffffffff
        ],
    };
    const FIVE: Self = Self {
        limbs: [
            0x96bc8c8cffffffed,
            0x74c2a54b49f7778e,
            0xfffffffffffffffd,
            0x3fffffffffffffff
        ],
    };
    const NEG_ONE: Self = Self {
        limbs: [0x311bac8400000004, 0x891a63f02652a376, 0, 0]
    };

    const MULTIPLICATIVE_SUBGROUP_GENERATOR: Self = Self::FIVE;

    const ALPHA: Self = Self::FIVE;

    const TWO_ADICITY: usize = 32;

    /// T = (ORDER - 1) / 2^TWO_ADICITY  in Monty form
    const T: Self = Self {
        limbs: [
            0x8c46eb2100000001,
            0x224698fc0994a8dd,
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

impl Ord for VestaBase {
    fn cmp(&self, other: &Self) -> Ordering {
        self.cmp_helper(other)
    }
}

impl PartialOrd for VestaBase {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Display for VestaBase {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}", field_to_biguint(*self))
    }
}

impl Debug for VestaBase {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "VestaBase({})", field_to_biguint(*self))
    }
}

#[cfg(test)]
mod tests {
    use crate::test_arithmetic;
    use crate::Field;
    use crate::VestaBase;

    #[test]
    fn primitive_root_order() {
        for n_power in 0..10 {
            let root = VestaBase::primitive_root_of_unity(n_power);
            let order = VestaBase::generator_order(root);
            assert_eq!(order, 1 << n_power, "2^{}'th primitive root", n_power);
        }
    }

    test_arithmetic!(crate::VestaBase);
}
