use rand::Rng;
use std::cmp::Ordering;
use std::cmp::Ordering::{Less, Equal};
use std::convert::TryInto;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::fmt;
use std::fmt::{Debug, Display, Formatter};

use crate::{cmp, field_to_biguint,
            rand_range, rand_range_from_rng,
            DairaRepr, Field};

/// An element of the Tweedledee group's base field.
#[derive(Copy, Clone, Hash, Default)]
pub struct TweedledeeBase {
    /// Daira representation, encoded with little-endian u64 limbs.
    pub limbs: [u64; 4],
}

impl DairaRepr for TweedledeeBase {
    /// ORDER = 2^254 + C
    const ORDER: [u64; 4] = [
        0x842cafd400000001,
        0x38aa127696286c9,
        0x0,
        0x4000000000000000
    ];

    const C: [u64; 2] = [0x842cafd400000001, 0x38aa127696286c9];

    /// Multiplication lookup table: ith element is (i+1)*C^2
    const C_SQR_TBL: [[u64; 4]; 15] = [
    [0x8595fa800000001, 0x27e16a565c689523, 0x3f4c0e6fcb03aa0a, 0xc8ad910688601],
    [0x10b2bf5000000002, 0x4fc2d4acb8d12a46, 0x7e981cdf96075414, 0x1915b220d10c02],
    [0x190c1ef800000003, 0x77a43f031539bf69, 0xbde42b4f610afe1e, 0x25a08b31399203],
    [0x21657ea000000004, 0x9f85a95971a2548c, 0xfd3039bf2c0ea828, 0x322b6441a21804],
    [0x29bede4800000005, 0xc76713afce0ae9af, 0x3c7c482ef7125232, 0x3eb63d520a9e06],
    [0x32183df000000006, 0xef487e062a737ed2, 0x7bc8569ec215fc3c, 0x4b411662732407],
    [0x3a719d9800000007, 0x1729e85c86dc13f5, 0xbb14650e8d19a647, 0x57cbef72dbaa08],
    [0x42cafd4000000008, 0x3f0b52b2e344a918, 0xfa60737e581d5051, 0x6456c883443009],
    [0x4b245ce800000009, 0x66ecbd093fad3e3b, 0x39ac81ee2320fa5b, 0x70e1a193acb60b],
    [0x537dbc900000000a, 0x8ece275f9c15d35e, 0x78f8905dee24a465, 0x7d6c7aa4153c0c],
    [0x5bd71c380000000b, 0xb6af91b5f87e6881, 0xb8449ecdb9284e6f, 0x89f753b47dc20d],
    [0x64307be00000000c, 0xde90fc0c54e6fda4, 0xf790ad3d842bf879, 0x96822cc4e6480e],
    [0x6c89db880000000d, 0x6726662b14f92c7, 0x36dcbbad4f2fa284, 0xa30d05d54ece10],
    [0x74e33b300000000e, 0x2e53d0b90db827ea, 0x7628ca1d1a334c8e, 0xaf97dee5b75411],
    [0x7d3c9ad80000000f, 0x56353b0f6a20bd0d, 0xb574d88ce536f698, 0xbc22b7f61fda12]
    ];

    /// K = 2*C, M = ORDER, K_M = K*M
    const K_M: [u64; 6] = [
        0x10b2bf5000000002, 0x4fc2d4acb8d12a46, 0x7e981cdf96075414,
        0x801915b220d10c02, 0xc21657ea00000000, 0x1c55093b4b14364
    ];
}

impl TweedledeeBase {
    pub fn from_canonical(c: [u64; 4]) -> Self {
        debug_assert!(cmp(c, Self::ORDER) == Less, "Given non-canonical input: {:?}", c);
        Self { limbs: c }
    }

    pub fn to_canonical(&self) -> [u64; 4] {
        Self::daira_to_canonical(self.limbs)
    }
}

impl Add<TweedledeeBase> for TweedledeeBase {
    type Output = Self;

    fn add(self, rhs: TweedledeeBase) -> Self {
        Self { limbs: Self::daira_add(self.limbs, rhs.limbs) }
    }
}

impl Sub<TweedledeeBase> for TweedledeeBase {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self { limbs: Self::daira_sub(self.limbs, rhs.limbs) }
    }
}

impl Mul<TweedledeeBase> for TweedledeeBase {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self { limbs: Self::daira_multiply(self.limbs, rhs.limbs) }
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
        Self { limbs: Self::daira_neg(self.limbs) }
    }
}

impl Field for TweedledeeBase {
    const BITS: usize = 255;
    const BYTES: usize = 32;
    const ZERO: Self = Self { limbs: <Self as DairaRepr>::ZERO };
    const ONE: Self = Self { limbs: <Self as DairaRepr>::ONE };
    const TWO: Self = Self { limbs: [2u64, 0, 0, 0] };
    const THREE: Self = Self { limbs: [3u64, 0, 0, 0] };
    const FOUR: Self = Self { limbs: [4u64, 0, 0, 0] };
    const FIVE: Self = Self { limbs: [5u64, 0, 0, 0] };

    const NEG_ONE: Self = Self {
        limbs: [0x842cafd400000000, 0x38aa127696286c9, 0x0, 0x4000000000000000]
    };

    const MULTIPLICATIVE_SUBGROUP_GENERATOR: Self = Self::FIVE;

    const ALPHA: Self = Self::FIVE;

    const TWO_ADICITY: usize = 34;

    /// 1684996666696914987166688442938726917102595538363933628829375605749
    const T: Self = Self {
        limbs: [0xda58a1b2610b2bf5, 0xe2a849, 0x0, 0x10000000]
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
            limbs: Self::daira_inverse(self.limbs)
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
            limbs: Self::daira_square(self.limbs),
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

/// NB: We need to implement this so that is_zero and == call
/// to_canonical (via cmp_helper) before comparing numbers.
impl PartialEq for TweedledeeBase {
    fn eq(&self, other: &Self) -> bool {
        self.cmp_helper(other) == Equal
    }
}

impl Eq for TweedledeeBase { }


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
        let small = <TweedledeeBase as Field>::ONE.to_canonical_u64_vec();
        assert!(TweedledeeBase::is_valid_canonical_u64(&small));

        let limbs = vec![1, 2, 3, 4, 5];
        assert_eq!(TweedledeeBase::is_valid_canonical_u64(&limbs), false);
    }

    test_arithmetic!(crate::TweedledeeBase);
}
