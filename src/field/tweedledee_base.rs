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

/// An element of the Tweedledee group's base field.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Default)]
pub struct TweedledeeBase {
    /// Montgomery representation, encoded with little-endian u64 limbs.
    pub limbs: [u64; 4],
}

impl MontyRepr for TweedledeeBase {
    /// The order of the field: 28948022309329048855892746252171976963322203655954433126947083963168578338817
    const ORDER: [u64; 4] = [
        9524180637049683969,
        255193519543715529,
        0,
        4611686018427387904,
    ];

    /// Twice the order of the field: 57896044618658097711785492504343953926644407311908866253894167926337156677634
    #[allow(dead_code)]
    const ORDER_X2: [u64; 4] = [
        601617200389816322,
        510387039087431059,
        0,
        9223372036854775808,
    ];

    /// R in the context of the Montgomery reduction, i.e. 2^256 % |F|.
    const R: [u64; 4] = [
        8320946236270051325,
        17681163515078405027,
        18446744073709551615,
        4611686018427387903,
    ];

    /// R^2 in the context of the Montgomery reduction, i.e. 2^(256*2) % |F|.
    const R2: [u64; 4] = [
        9625875206237061136,
        9085631154807722544,
        17636350113745641634,
        56485833733595155,
    ];

    /// R^3 in the context of the Montgomery reduction, i.e. 2^(256*3) % |F|.
    const R3: [u64; 4] = [
        11971961131424865118,
        6311318431551332850,
        14638507591886519234,
        739379759776372087,
    ];

    /// In the context of Montgomery multiplication, µ = -|F|^-1 mod 2^64.
    const MU: u64 = 9524180637049683967;
}

impl TweedledeeBase {
    pub fn from_canonical(c: [u64; 4]) -> Self {
        Self { limbs: Self::from_monty(c) }
    }

    pub fn to_canonical(&self) -> [u64; 4] {
        Self::to_monty(self.limbs)
    }
}

impl Add<TweedledeeBase> for TweedledeeBase {
    type Output = Self;

    fn add(self, rhs: TweedledeeBase) -> Self::Output {
        Self { limbs: Self::monty_add(self.limbs, rhs.limbs) }
    }
}

impl Sub<TweedledeeBase> for TweedledeeBase {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self { limbs: Self::monty_sub(self.limbs, rhs.limbs) }
    }
}

impl Mul<TweedledeeBase> for TweedledeeBase {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self { limbs: Self::monty_multiply(self.limbs, rhs.limbs) }
    }
}

impl Div<TweedledeeBase> for TweedledeeBase {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        self * rhs.multiplicative_inverse().expect("No inverse")
    }
}

impl Neg for TweedledeeBase {
    type Output = Self;

    fn neg(self) -> Self {
        Self { limbs: Self::monty_neg(self.limbs) }
    }
}

impl Field for TweedledeeBase {
    const BITS: usize = 255;
    const BYTES: usize = 32;
    const ZERO: Self = Self { limbs: <Self as MontyRepr>::ZERO };
    const ONE: Self = Self { limbs: <Self as MontyRepr>::ONE };
    const TWO: Self = Self {
        limbs: [
            7117711835490418681,
            16660389436903542909,
            18446744073709551615,
            4611686018427387903,
        ],
    };
    const THREE: Self = Self {
        limbs: [
            5914477434710786037,
            15639615358728680791,
            18446744073709551615,
            4611686018427387903,
        ],
    };
    const FOUR: Self = Self {
        limbs: [
            4711243033931153393,
            14618841280553818673,
            18446744073709551615,
            4611686018427387903,
        ],
    };
    const FIVE: Self = Self {
        limbs: [
            3508008633151520749,
            13598067202378956555,
            18446744073709551615,
            4611686018427387903,
        ],
    };
    const NEG_ONE: Self = Self {
        limbs: [1203234400779632644, 1020774078174862118, 0, 0],
    };

    const MULTIPLICATIVE_SUBGROUP_GENERATOR: Self = Self::FIVE;

    const ALPHA: Self = Self::FIVE;

    const TWO_ADICITY: usize = 34;

    /// 1684996666696914987166688442938726917102595538363933628829375605749
    const T: Self = Self {
        limbs: [
            9524180637049683969,
            255193519543715529,
            0,
            4611686017353646080,
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

impl Ord for TweedledeeBase {
    fn cmp(&self, other: &Self) -> Ordering {
        self.cmp_helper(other)
    }
}

impl PartialOrd for TweedledeeBase {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Display for TweedledeeBase {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}", field_to_biguint(*self))
    }
}

impl Debug for TweedledeeBase {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "TweedledeeBase {}", field_to_biguint(*self))
    }
}

#[cfg(test)]
mod tests {
    use crate::test_arithmetic;
    use crate::Field;
    use crate::TweedledeeBase;
    use crate::MontyRepr; // This is just to access ORDER_X2.

    #[test]
    fn primitive_root_order() {
        for n_power in 0..10 {
            let root = TweedledeeBase::primitive_root_of_unity(n_power);
            let order = TweedledeeBase::generator_order(root);
            assert_eq!(order, 1 << n_power, "2^{}'th primitive root", n_power);
        }
    }

    #[test]
    fn valid_canonical_vec() {
        let small = TweedledeeBase::ONE.to_canonical_u64_vec();
        assert!(TweedledeeBase::is_valid_canonical_u64(&small));

        let big = TweedledeeBase::ORDER_X2.to_vec();
        assert_eq!(TweedledeeBase::is_valid_canonical_u64(&big), false);

        let limbs = vec![1, 2, 3, 4, 5];
        assert_eq!(TweedledeeBase::is_valid_canonical_u64(&limbs), false);
    }

    test_arithmetic!(crate::TweedledeeBase);
}
