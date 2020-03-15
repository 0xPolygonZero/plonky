use std::cmp::Ordering::Less;
use std::convert::TryInto;
use std::ops::{Add, Div, Mul, Neg, Sub};

use unroll::unroll_for_loops;

use crate::{add_4_4_no_overflow, cmp_4_4, Field, rand_range_4, sub_4_4, TwoAdicField};
use crate::bigint_inverse::nonzero_multiplicative_inverse_4;

/// An element of the Tweedledee group's base field.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct TweedledeeBase {
    /// Montgomery representation, encoded with little-endian u64 limbs.
    pub limbs: [u64; 4],
}

impl TweedledeeBase {
    /// The order of the field: 28948022309329048855892746252171976963322203655954433126947083963168578338817
    const ORDER: [u64; 4] = [9524180637049683969, 255193519543715529, 0, 4611686018427387904];

    /// Twice the order of the field:
    const ORDER_X2: [u64; 4] = [601617200389816322, 510387039087431059, 0, 9223372036854775808];

    /// R in the context of the Montgomery reduction, i.e. 2^256 % |F|.
    const R: [u64; 4] = [8320946236270051325, 17681163515078405027, 18446744073709551615, 4611686018427387903];

    /// R^2 in the context of the Montgomery reduction, i.e. 2^(256*2) % |F|.
    const R2: [u64; 4] = [9625875206237061136, 9085631154807722544, 17636350113745641634, 56485833733595155];

    /// R^3 in the context of the Montgomery reduction, i.e. 2^(256*3) % |F|.
    const R3: [u64; 4] = [11971961131424865118, 6311318431551332850, 14638507591886519234, 739379759776372087];

    /// In the context of Montgomery multiplication, µ = -|F|^-1 mod 2^64.
    const MU: u64 = 9524180637049683967;

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

impl Add<TweedledeeBase> for TweedledeeBase {
    type Output = Self;

    fn add(self, rhs: TweedledeeBase) -> Self::Output {
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

impl Sub<TweedledeeBase> for TweedledeeBase {
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

impl Mul<TweedledeeBase> for TweedledeeBase {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self { limbs: Self::montgomery_multiply(self.limbs, rhs.limbs) }
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
        if self == Self::ZERO {
            Self::ZERO
        } else {
            Self { limbs: sub_4_4(Self::ORDER, self.limbs) }
        }
    }
}

impl Field for TweedledeeBase {
    const BITS: usize = 255;
    const ZERO: Self = Self { limbs: [0; 4] };
    const ONE: Self = Self { limbs: [8320946236270051325, 17681163515078405027, 18446744073709551615, 4611686018427387903] };
    const TWO: Self = Self { limbs: [7117711835490418681, 16660389436903542909, 18446744073709551615, 4611686018427387903] };
    const THREE: Self = Self { limbs: [5914477434710786037, 15639615358728680791, 18446744073709551615, 4611686018427387903] };

    fn to_canonical_vec(&self) -> Vec<u64> {
        self.to_canonical().to_vec()
    }

    fn from_canonical_vec(v: Vec<u64>) -> Self {
        Self::from_canonical(v[..].try_into().unwrap())
    }

    fn from_canonical_u64(n: u64) -> Self {
        Self::from_canonical([n, 0, 0, 0])
    }

    fn multiplicative_inverse_assuming_nonzero(&self) -> Self {
        // Let x R = self. We compute M((x R)^-1, R^3) = x^-1 R^-1 R^3 R^-1 = x^-1 R.
        let self_r_inv = nonzero_multiplicative_inverse_4(self.limbs, Self::ORDER);
        Self { limbs: Self::montgomery_multiply(self_r_inv, Self::R3) }
    }

    fn rand() -> Self {
        Self { limbs: rand_range_4(Self::ORDER) }
    }
}

impl TwoAdicField for TweedledeeBase {
    const TWO_ADICITY: usize = 34;

    fn primitive_root_of_unity(n_power: usize) -> Self {
        assert!(n_power <= Self::TWO_ADICITY);
        todo!()
    }
}
