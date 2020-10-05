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

/// An element of the Tweedledum group's base field.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Default)]
pub struct TweedledumBase {
    /// Montgomery representation, encoded with little-endian u64 limbs.
    pub limbs: [u64; 4],
}

impl MontyRepr for TweedledumBase {
    /// The order of the field: 28948022309329048855892746252171976963322203655955319056773317069363642105857
    const ORDER: [u64; 4] = [
        11619397960441266177,
        255193519591741881,
        0,
        4611686018427387904,
    ];

    /// Twice the order of the field: 57896044618658097711785492504343953926644407311910638113546634138727284211714
    #[allow(dead_code)]
    const ORDER_X2: [u64; 4] = [
        4792051847172980738,
        510387039183483763,
        0,
        9223372036854775808,
    ];

    /// R in the context of the Montgomery reduction, i.e. 2^256 % |F|.
    const R: [u64; 4] = [
        2035294266095304701,
        17681163514934325971,
        18446744073709551615,
        4611686018427387903,
    ];

    /// R^2 in the context of the Montgomery reduction, i.e. 2^(256*2) % |F|.
    const R2: [u64; 4] = [
        2885853259929485328,
        10494584067553537908,
        15959394653775906393,
        56485833754855950,
    ];

    /// R^3 in the context of the Montgomery reduction, i.e. 2^(256*3) % |F|.
    const R3: [u64; 4] = [
        11023471670160566071,
        18013763770685241468,
        7203328081223416457,
        2412999303287602290,
    ];

    /// In the context of Montgomery multiplication, µ = -|F|^-1 mod 2^64.
    const MU: u64 = 11619397960441266175;
}

impl TweedledumBase {
    pub fn from_canonical(c: [u64; 4]) -> Self {
        Self { limbs: Self::from_monty(c) }
    }

    pub fn to_canonical(&self) -> [u64; 4] {
        Self::to_monty(self.limbs)
    }
}

impl Add<TweedledumBase> for TweedledumBase {
    type Output = Self;

    fn add(self, rhs: TweedledumBase) -> Self::Output {
        Self { limbs: Self::monty_add(self.limbs, rhs.limbs) }
    }
}

impl Sub<TweedledumBase> for TweedledumBase {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self { limbs: Self::monty_sub(self.limbs, rhs.limbs) }
    }
}

impl Mul<TweedledumBase> for TweedledumBase {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self { limbs: Self::monty_multiply(self.limbs, rhs.limbs) }
    }
}

impl Div<TweedledumBase> for TweedledumBase {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        self * rhs.multiplicative_inverse().expect("No inverse")
    }
}

impl Neg for TweedledumBase {
    type Output = Self;

    fn neg(self) -> Self {
        Self { limbs: Self::monty_neg(self.limbs) }
    }
}

impl Field for TweedledumBase {
    const BITS: usize = 255;
    const BYTES: usize = 32;
    const ZERO: Self = Self { limbs: <Self as MontyRepr>::ZERO };
    const ONE: Self = Self { limbs: <Self as MontyRepr>::ONE };
    const TWO: Self = Self {
        limbs: [
            10897934645458894841,
            16660389436567358444,
            18446744073709551615,
            4611686018427387903,
        ],
    };
    const THREE: Self = Self {
        limbs: [
            1313830951112933365,
            15639615358200390918,
            18446744073709551615,
            4611686018427387903,
        ],
    };
    const FOUR: Self = Self {
        limbs: [
            10176471330476523505,
            14618841279833423391,
            18446744073709551615,
            4611686018427387903,
        ],
    };
    const FIVE: Self = Self {
        limbs: [
            592367636130562029,
            13598067201466455865,
            18446744073709551615,
            4611686018427387903,
        ],
    };
    const NEG_ONE: Self = Self {
        limbs: [9584103694345961476, 1020774078366967526, 0, 0],
    };

    const MULTIPLICATIVE_SUBGROUP_GENERATOR: Self = Self::FIVE;

    const ALPHA: Self = Self::FIVE;

    const TWO_ADICITY: usize = 33;

    /// 3369993333393829974333376885877453834205191076727970393464588218993
    const T: Self = Self {
        limbs: [
            11619397960441266177,
            255193519591741881,
            0,
            4611686016279904256,
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

impl Ord for TweedledumBase {
    fn cmp(&self, other: &Self) -> Ordering {
        self.cmp_helper(other)
    }
}

impl PartialOrd for TweedledumBase {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Display for TweedledumBase {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}", field_to_biguint(*self))
    }
}

impl Debug for TweedledumBase {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "TweedledumBase({})", field_to_biguint(*self))
    }
}

#[cfg(test)]
mod tests {
    use crate::test_arithmetic;
    use crate::Field;
    use crate::TweedledumBase;

    #[test]
    fn primitive_root_order() {
        for n_power in 0..10 {
            let root = TweedledumBase::primitive_root_of_unity(n_power);
            let order = TweedledumBase::generator_order(root);
            assert_eq!(order, 1 << n_power, "2^{}'th primitive root", n_power);
        }
    }

    test_arithmetic!(crate::TweedledumBase);
}
