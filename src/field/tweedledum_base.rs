use rand::Rng;
use std::cmp::Ordering::Less;
use std::convert::TryInto;
use std::ops::{Add, Div, Mul, Neg, Sub};

use unroll::unroll_for_loops;

use crate::nonzero_multiplicative_inverse_4;
use crate::{add_4_4_no_overflow, cmp_4_4, field_to_biguint, rand_range_4, rand_range_4_from_rng, sub_4_4, Field};
use std::cmp::Ordering;
use std::fmt;
use std::fmt::{Debug, Display, Formatter};

/// An element of the Tweedledum group's base field.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Default)]
pub struct TweedledumBase {
    /// Montgomery representation, encoded with little-endian u64 limbs.
    pub limbs: [u64; 4],
}

impl TweedledumBase {
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

    pub fn from_canonical(c: [u64; 4]) -> Self {
        // We compute M(c, R^2) = c * R^2 * R^-1 = c * R.
        Self {
            limbs: Self::montgomery_multiply(c, Self::R2),
        }
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
                let result =
                    c[(i + j) % 5] as u128 + q as u128 * Self::ORDER[j] as u128 + carry as u128;
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

impl Add<TweedledumBase> for TweedledumBase {
    type Output = Self;

    fn add(self, rhs: TweedledumBase) -> Self::Output {
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

impl Sub<TweedledumBase> for TweedledumBase {
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

impl Mul<TweedledumBase> for TweedledumBase {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self {
            limbs: Self::montgomery_multiply(self.limbs, rhs.limbs),
        }
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
        if self == Self::ZERO {
            Self::ZERO
        } else {
            Self {
                limbs: sub_4_4(Self::ORDER, self.limbs),
            }
        }
    }
}

impl Field for TweedledumBase {
    const BITS: usize = 255;
    const BYTES: usize = 32;
    const ZERO: Self = Self { limbs: [0; 4] };
    const ONE: Self = Self { limbs: Self::R };
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

    fn is_valid_canonical_u64(v: &Vec<u64>) -> bool {
        v.len() == 4 && cmp_4_4(v[..].try_into().unwrap(), Self::ORDER) == Less
    }

    fn multiplicative_inverse_assuming_nonzero(&self) -> Self {
        // Let x R = self. We compute M((x R)^-1, R^3) = x^-1 R^-1 R^3 R^-1 = x^-1 R.
        let self_r_inv = nonzero_multiplicative_inverse_4(self.limbs, Self::ORDER);
        Self {
            limbs: Self::montgomery_multiply(self_r_inv, Self::R3),
        }
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
    use crate::test_square_root;
    use crate::{Field, TweedledumBase};
    use std::io::{Result, Read, BufReader};
    use std::fs::File;
    use crate::serialization::FromBytes;
    use std::ops::{Add,Sub,Mul,Div,Neg};

    #[test]
    fn primitive_root_order() {
        for n_power in 0..10 {
            let root = TweedledumBase::primitive_root_of_unity(n_power);
            let order = TweedledumBase::generator_order(root);
            assert_eq!(order, 1 << n_power, "2^{}'th primitive root", n_power);
        }
    }

    test_square_root!(TweedledumBase);

    /* Read 4 bytes from f and interpret little endian as i32. */
    fn read_i32<R: Read>(f: &mut R) -> Result<i32> {
        let mut i32_buf = [0; 4];
        f.read_exact(&mut i32_buf)?;
        Ok(i32::from_le_bytes(i32_buf))
    }

    // TODO: There is probably a cleaner way to express this.
    // TODO: I don't understand the signature for parameter f.
    fn read_vec<R: Read>(mut f: &mut R, nelts: usize) -> Result<Vec<TweedledumBase>> {
        let mut res = Vec::new();
        for _i in 0 .. nelts {
            res.push(TweedledumBase::read(&mut f)?);
        }
        Ok(res)
    }

    fn run_test_cases<F>(op_str: &str, op: BinFn) -> Result<()>
        where BinFn: Fn(TweedledumBase, TweedledumBase) -> TweedledumBase
    {
        let file_prefix = "src/field/test-data/tweedledum_";
        let input_file = File::open([file_prefix, op_str].concat())?;
        let mut reader = BufReader::new(input_file);

        // TODO: Check whether I can just pass a mut here rather than &mut.
        let bytes_per_elt = read_i32(&mut reader)?;
        assert_eq!(bytes_per_elt as usize, TweedledumBase::BYTES,
                   "mismatch in expected size");
        let n_input_elts = read_i32(&mut reader)? as usize;
        let n_outputs_per_op = read_i32(&mut reader)?;
        assert_eq!(n_outputs_per_op, 1, "unexpected value for #outputs/op");

        let inputs = read_vec(&mut reader, n_input_elts)?;

        for i in 0 .. n_input_elts {
            // Iterator over inputs rotated right by i places. Since
            // cycle().skip(i) rotates left by i, we need to rotate by
            // n_input_elts - i.
            let shifted_inputs = inputs.iter().cycle().skip(n_input_elts - i);
            // Calculate pointwise operation
            let output = inputs.iter().zip(shifted_inputs).map(
                |(&x, &y)| op(x, y));
            // Read expected outputs
            let expected_output = read_vec(&mut reader, n_input_elts)?;
            // Compare expected outputs with actual outputs
            assert!(output.zip(expected_output.iter()).all(|(x, &y)| x == y),
                    "output differs from expected at rotation {}", i);
        }
        Ok(())
    }

    #[test]
    fn addition() -> Result<()> {
        run_test_cases("add", TweedledumBase::add)
    }

    #[test]
    fn subtraction() -> Result<()> {
        run_test_cases("sub", TweedledumBase::sub)
    }

        //run_test_cases("neg", TweedledumBase::neg);

    #[test]
    fn multiplication() -> Result<()> {
        run_test_cases("mul", TweedledumBase::mul)
    }

    #[test]
    fn square() -> Result<()> {
        run_test_cases("sqr", TweedledumBase::mul)
    }

    #[test]
    fn division() -> Result<()> {
        run_test_cases("div", TweedledumBase::div)
    }
}
