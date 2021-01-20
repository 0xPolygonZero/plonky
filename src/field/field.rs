use std::cmp::{Ordering, min};
use std::cmp::Ordering::Equal;
use std::collections::HashSet;
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::ops::{Add, Div, Mul, Neg, Sub};

use anyhow::{Error, Result};
use num::{BigUint, Integer, One, Zero};
use rand::Rng;
use serde::{de::DeserializeOwned, Serialize};

use crate::{biguint_to_field, Curve, field_to_biguint, ProjectivePoint};

pub trait Field:
    'static
    + Sized
    + Copy
    + Ord
    + Hash
    + Send
    + Sync
    + Debug
    + Display
    + Default
    + Neg<Output = Self>
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<Self, Output = Self>
    + Div<Self, Output = Self>
    + Serialize
    + DeserializeOwned
{
    const BITS: usize;
    const BYTES: usize;

    const ZERO: Self;
    const ONE: Self;
    const TWO: Self;
    const THREE: Self;
    const FOUR: Self;
    const FIVE: Self;
    const NEG_ONE: Self;

    const MULTIPLICATIVE_SUBGROUP_GENERATOR: Self;

    /// An element `a` such that `x^a` is a permutation in this field. Although not strictly
    /// required, the smallest such `a` should be configured to minimize the cost of evaluating the
    /// monomial.
    const ALPHA: Self;

    const TWO_ADICITY: usize;

    /// `T = (ORDER - 1) / 2^TWO_ADICITY`
    const T: Self;

    fn to_canonical_u64_vec(&self) -> Vec<u64>;

    fn to_canonical_u32_vec(&self) -> Vec<u32> {
        let mut limbs = Vec::new();
        for u64_limb in self.to_canonical_u64_vec() {
            limbs.push(u64_limb as u32);
            limbs.push(u64_limb.overflowing_shr(32).0 as u32);
        }
        limbs
    }

    fn to_canonical_u8_vec(&self) -> Vec<u8> {
        let mut limbs = Vec::new();
        for u64_limb in self.to_canonical_u64_vec() {
            limbs.push(u64_limb as u8);
            limbs.push(u64_limb.overflowing_shr(8).0 as u8);
            limbs.push(u64_limb.overflowing_shr(16).0 as u8);
            limbs.push(u64_limb.overflowing_shr(24).0 as u8);
            limbs.push(u64_limb.overflowing_shr(32).0 as u8);
            limbs.push(u64_limb.overflowing_shr(40).0 as u8);
            limbs.push(u64_limb.overflowing_shr(48).0 as u8);
            limbs.push(u64_limb.overflowing_shr(56).0 as u8);
        }
        limbs
    }

    fn from_canonical_u8_vec(u8_limbs: Vec<u8>) -> Result<Self> {
        let mut u64_chunks = Vec::new();
        for u8_chunk in u8_limbs.chunks(8) {
            u64_chunks.push(
                (u8_chunk[7] as u64) << 56
                    | (u8_chunk[6] as u64) << 48
                    | (u8_chunk[5] as u64) << 40
                    | (u8_chunk[4] as u64) << 32
                    | (u8_chunk[3] as u64) << 24
                    | (u8_chunk[2] as u64) << 16
                    | (u8_chunk[1] as u64) << 8
                    | (u8_chunk[0] as u64),
            );
        }

        if Self::is_valid_canonical_u64(&u64_chunks) {
            Ok(Self::from_canonical_u64_vec(u64_chunks))
        } else {
            Err(Error::msg("Out of range"))
        }
    }

    fn to_canonical_bool_vec(&self) -> Vec<bool> {
        let mut limbs = Vec::new();
        for u64_limb in self.to_canonical_u64_vec() {
            for i in 0..64 {
                limbs.push((u64_limb.overflowing_shr(i).0 & 1) != 0);
            }
        }
        limbs
    }

    fn from_canonical_u64_vec(v: Vec<u64>) -> Self;

    fn from_canonical_u32_vec(u32_limbs: Vec<u32>) -> Self {
        let mut u64_chunks = Vec::new();
        for u32_chunk in u32_limbs.chunks(2) {
            u64_chunks.push((u32_chunk[1] as u64) << 32 | u32_chunk[0] as u64);
        }
        Self::from_canonical_u64_vec(u64_chunks)
    }

    fn from_canonical_u64(n: u64) -> Self;

    fn from_canonical_u32(n: u32) -> Self {
        Self::from_canonical_u64(n as u64)
    }

    fn from_canonical_usize(n: usize) -> Self {
        Self::from_canonical_u64(n as u64)
    }

    fn from_canonical_bool(b: bool) -> Self {
        if b {
            Self::ONE
        } else {
            Self::ZERO
        }
    }

    fn is_valid_canonical_u64(v: &[u64]) -> bool;

    #[inline(always)]
    fn is_zero(&self) -> bool {
        *self == Self::ZERO
    }

    #[inline(always)]
    fn is_one(&self) -> bool {
        *self == Self::ONE
    }

    #[inline(always)]
    fn is_nonzero(&self) -> bool {
        *self != Self::ZERO
    }

    #[inline]
    fn multiplicative_inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(self.multiplicative_inverse_assuming_nonzero())
        }
    }

    fn multiplicative_inverse_assuming_nonzero(&self) -> Self;

    #[inline(always)]
    fn square(&self) -> Self {
        *self * *self
    }

    #[inline(always)]
    fn cube(&self) -> Self {
        self.square() * *self
    }

    #[inline(always)]
    fn double(&self) -> Self {
        *self * Self::TWO
    }

    #[inline(always)]
    fn triple(&self) -> Self {
        *self * Self::THREE
    }

    #[inline(always)]
    fn quadruple(&self) -> Self {
        *self * Self::FOUR
    }

    fn scale_slice(&self, slice: &[Self]) -> Vec<Self> {
        slice.iter().map(|&x| *self * x).collect()
    }

    fn scale_proj_point_slice<C>(&self, slice: &[ProjectivePoint<C>]) -> Vec<ProjectivePoint<C>>
    where
        C: Curve<ScalarField = Self>,
    {
        slice.iter().map(|&p| C::convert(*self) * p).collect()
    }

    fn add_slices(a: &[Self], b: &[Self]) -> Vec<Self> {
        assert_eq!(a.len(), b.len());
        a.iter()
            .zip(b.iter())
            .map(|(&a_i, &b_i)| a_i + b_i)
            .collect()
    }

    fn inner_product(a: &[Self], b: &[Self]) -> Self {
        assert_eq!(a.len(), b.len());
        let mut sum = Self::ZERO;
        for (&a_i, &b_i) in a.iter().zip(b.iter()) {
            sum = sum + a_i * b_i;
        }
        sum
    }

    fn batch_multiplicative_inverse_opt<F: Field>(x: &[F]) -> Vec<Option<F>> {
        let n = x.len();
        let mut x_nonzero = Vec::with_capacity(n);
        let mut index_map = Vec::with_capacity(n);

        for x_i in x {
            if !x_i.is_zero() {
                index_map.push(x_nonzero.len());
                x_nonzero.push(*x_i);
            } else {
                // Push an intentionally out-of-bounds index to make sure it isn't used.
                index_map.push(n);
            }
        }

        let x_inv = F::batch_multiplicative_inverse(&x_nonzero);

        let mut result = Vec::with_capacity(n);
        for i in 0..n {
            if !x[i].is_zero() {
                result.push(Some(x_inv[index_map[i]]));
            } else {
                result.push(None);
            }
        }
        result
    }

    fn batch_multiplicative_inverse<F: Field>(x: &[F]) -> Vec<F> {
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

        let mut a_inv = vec![F::ZERO; n];
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

    fn cyclic_subgroup_unknown_order(generator: Self) -> Vec<Self> {
        let mut subgroup = Vec::new();
        let mut current = generator;
        // Start with the identity, thn add the generator until we reach the identity again.
        subgroup.push(Self::ONE);
        while current != Self::ONE {
            subgroup.push(current);
            current = current * generator;
        }
        subgroup
    }

    fn cyclic_subgroup_known_order(generator: Self, order: usize) -> Vec<Self> {
        let mut subgroup = Vec::new();
        let mut current = Self::ONE;
        for _i in 0..order {
            subgroup.push(current);
            current = current * generator;
        }
        subgroup
    }

    fn generator_order(generator: Self) -> usize
    where
        Self: Hash,
    {
        Self::cyclic_subgroup_unknown_order(generator).len()
    }

    fn exp(&self, power: Self) -> Self {
        let mut power_bits = power.num_bits();
        let mut current = *self;
        let mut product = Self::ONE;

        for limb in power.to_canonical_u64_vec().iter() {
            // To minimize branching, it's better to iterate up to the min than to conditionally
            // break inside the loop.
            for j in 0..min(64, power_bits) {
                if (limb >> j & 1) != 0 {
                    product = product * current;
                }
                current = current.square();
            }
            if power_bits >= 64 {
                power_bits -= 64;
            } else {
                break;
            }
        }
        product
    }

    fn exp_u32(&self, power: u32) -> Self {
        self.exp(Self::from_canonical_u32(power))
    }

    fn exp_usize(&self, power: usize) -> Self {
        self.exp(Self::from_canonical_usize(power))
    }

    fn kth_root_u32(&self, k: u32) -> Self {
        self.kth_root(Self::from_canonical_u32(k))
    }

    /// Computes `x^(1/k)`. Assumes that `x^k` is a permutation in this field; undefined behavior
    /// otherwise.
    fn kth_root(&self, k: Self) -> Self {
        // By Fermat's little theorem, x^p = x and x^(p - 1) = 1, so x^(p + n(p - 1)) = x for any n.
        // Our assumption that the k'th root operation is a permutation implies gcd(p - 1, k) = 1,
        // so there exists some n such that p + n(p - 1) is a multiple of k. Once we find such an n,
        // we can rewrite the above as
        //    x^((p + n(p - 1))/k)^k = x,
        // implying that x^((p + n(p - 1))/k) is a k'th root of x.

        let p_minus_1_bu = field_to_biguint(Self::NEG_ONE);
        let k_bu = field_to_biguint(k);
        let mut n = BigUint::zero();
        let mut numerator_bu = &p_minus_1_bu + BigUint::one();

        while n < k_bu {
            // We can safely increment first, thus skipping the check for n=0, since n=0 will never
            // satisfy the relation above.
            n += BigUint::one();
            numerator_bu += &p_minus_1_bu;

            if numerator_bu.is_multiple_of(&k_bu) {
                let power_bu = numerator_bu.div_floor(&k_bu).mod_floor(&p_minus_1_bu);
                return self.exp(biguint_to_field(power_bu));
            }
        }

        panic!(
            "x^{} and x^(1/{}) are not permutations in this field, or we have a bug!",
            k, k
        );
    }

    fn is_quadratic_residue(&self) -> bool {
        if self.is_zero() {
            return true;
        }
        // This is based on Euler's criterion.
        let power = biguint_to_field(field_to_biguint(Self::NEG_ONE) / 2u8);
        let exp = self.exp(power);
        if exp == Self::ONE {
            return true;
        }
        if exp == Self::NEG_ONE {
            return false;
        }
        panic!("Number theory is a lie!")
    }

    /// The number of bits in the binary encoding of this field element.
    fn num_bits(&self) -> usize {
        // Search for the most significant nonzero limb.
        let limbs = self.to_canonical_u64_vec();
        for (i, &limb) in limbs.iter().rev().enumerate() {
            if limb != 0 {
                // This can be understood as the size of all limbs (limbs.len() * 64), minus the
                // leading zeros from all-zero limbs (i * 64), minus the leading zeros from this
                // limb (limb.leading_zeros()).
                return (limbs.len() - i) * 64 - limb.leading_zeros() as usize;
            }
        }
        0
    }

    /// Like `Ord::cmp`. We can't implement `Ord` directly due to Rust's trait coherence rules, so
    /// instead we provide this helper which implementations can use to trivially implement `Ord`.
    fn cmp_helper(&self, other: &Self) -> Ordering {
        let self_limbs = self.to_canonical_u64_vec().into_iter();
        let other_limbs = other.to_canonical_u64_vec().into_iter();

        let mut result = Equal;
        for (self_limb, other_limb) in self_limbs.zip(other_limbs) {
            let limb_ordering = self_limb.cmp(&other_limb);
            if limb_ordering != Equal {
                result = limb_ordering;
            }
        }
        result
    }

    fn rand() -> Self;

    fn rand_from_rng<R: Rng>(rng: &mut R) -> Self;

    /// Computes a `2^n_power`th primitive root of unity.
    fn primitive_root_of_unity(n_power: usize) -> Self {
        assert!(n_power <= Self::TWO_ADICITY);
        let base_root = Self::MULTIPLICATIVE_SUBGROUP_GENERATOR.exp(Self::T);
        base_root.exp(Self::from_canonical_u64(
            1u64 << (Self::TWO_ADICITY as u64 - n_power as u64),
        ))
    }

    /// If this is a quadratic residue, return an arbitrary (but deterministic) one of its square
    /// roots, otherwise return `None`.
    /// Inspired by implementation in https://github.com/scipr-lab/zexe/blob/85bae796a411077733ddeefda042d02f4b4772e5/algebra-core/src/fields/arithmetic.rs
    fn square_root(&self) -> Option<Self> {
        if self.is_zero() {
            Some(*self)
        } else if self.is_quadratic_residue() {
            let mut z = Self::MULTIPLICATIVE_SUBGROUP_GENERATOR.exp(Self::T);
            let mut w = self.exp((Self::T - Self::ONE) / Self::TWO);
            let mut x = w * *self;
            let mut b = x * w;

            let mut v = Self::TWO_ADICITY as usize;

            while !b.is_one() {
                let mut k = 0usize;
                let mut b2k = b;
                while !b2k.is_one() {
                    b2k = b2k.square();
                    k += 1;
                }
                let j = v - k - 1;
                w = z;
                for _ in 0..j {
                    w = w.square();
                }

                z = w.square();
                b = b * z;
                x = x * w;
                v = k;
            }
            Some(x)
        } else {
            None
        }
    }

    /// Return this field element re-encoded as an element of `F` if it fits, or `Err` if not.
    fn try_convert<F: Field>(&self) -> Result<F> {
        F::from_canonical_u8_vec(self.to_canonical_u8_vec())
    }

    fn try_convert_all<F: Field>(values: &[Self]) -> Result<Vec<F>> {
        values.iter().map(|f| f.try_convert::<F>()).collect()
    }
}

#[cfg(test)]
pub mod field_tests {
    use std::io::Result;

    use num::{BigUint, One, Zero};

    use crate::{biguint_to_field, Field, field_to_biguint};
    use crate::util::ceil_div_usize;

    /// Generates a series of non-negative integers less than
    /// `modulus` which cover a range of values and which will
    /// generate lots of carries, especially at `word_bits` word
    /// boundaries.
    pub fn test_inputs(modulus: &num::BigUint, word_bits: usize) -> Vec<num::BigUint> {
        assert!(word_bits == 32 || word_bits == 64);
        let modwords = ceil_div_usize(modulus.bits() as usize, word_bits);
        // Start with basic set close to zero: 0 .. 10
        const BIGGEST_SMALL: u32 = 10;
        let smalls: Vec<_> = (0..BIGGEST_SMALL).map(BigUint::from).collect();
        // ... and close to MAX: MAX - x
        let word_max = (BigUint::one() << word_bits) - 1u32;
        let bigs = smalls.iter().map(|x| &word_max - x).collect();
        let one_words = [smalls, bigs].concat();
        // For each of the one word inputs above, create a new one at word i.
        // TODO: Create all possible `modwords` combinations of those
        let multiple_words = (1..modwords)
            .flat_map(|i| {
                one_words
                    .iter()
                    .map(|x| x << (word_bits * i))
                    .collect::<Vec<BigUint>>()
            })
            .collect();
        let basic_inputs = [one_words, multiple_words].concat();

        // Biggest value that will fit in `modwords` words
        let maxval = (BigUint::one() << (modwords * word_bits)) - 1u32;
        // Inputs 'difference from' maximum value
        let diff_max = basic_inputs
            .iter()
            .map(|x| &maxval - x)
            .filter(|x| x < modulus)
            .collect();
        // Inputs 'difference from' modulus value
        let diff_mod = basic_inputs
            .iter()
            .filter(|x| *x < modulus && !x.is_zero())
            .map(|x| modulus - x)
            .collect();
        let basics = basic_inputs
            .into_iter()
            .filter(|x| x < modulus)
            .collect::<Vec<BigUint>>();
        [basics, diff_max, diff_mod].concat()

        // // There should be a nicer way to express the code above; something
        // // like this (and removing collect() calls from diff_max and diff_mod):
        // basic_inputs.into_iter()
        //     .chain(diff_max)
        //     .chain(diff_mod)
        //     .filter(|x| x < &modulus)
        //     .collect()
    }

    /// Apply the unary functions `op` and `expected_op`
    /// coordinate-wise to the inputs from `test_inputs(modulus,
    /// word_bits)` and panic if the two resulting vectors differ.
    pub fn run_unaryop_test_cases<F, UnaryOp, ExpectedOp>(
        modulus: &BigUint,
        word_bits: usize,
        op: UnaryOp,
        expected_op: ExpectedOp,
    ) -> Result<()>
    where
        F: Field,
        UnaryOp: Fn(F) -> F,
        ExpectedOp: Fn(BigUint) -> BigUint,
    {
        let inputs = test_inputs(modulus, word_bits);
        let expected = inputs.iter().map(|x| expected_op(x.clone()));
        let output = inputs
            .iter()
            .map(|x| field_to_biguint(op(biguint_to_field(x.clone()))));
        // Compare expected outputs with actual outputs
        assert!(
            output.zip(expected).all(|(x, y)| x == y),
            "output differs from expected"
        );
        Ok(())
    }

    /// Apply the binary functions `op` and `expected_op` to each pair
    /// in `zip(inputs, rotate_right(inputs, i))` where `inputs` is
    /// `test_inputs(modulus, word_bits)` and `i` ranges from 0 to
    /// `inputs.len()`.  Panic if the two functions ever give
    /// different answers.
    pub fn run_binaryop_test_cases<F, BinaryOp, ExpectedOp>(
        modulus: &BigUint,
        word_bits: usize,
        op: BinaryOp,
        expected_op: ExpectedOp,
    ) -> Result<()>
    where
        F: Field,
        BinaryOp: Fn(F, F) -> F,
        ExpectedOp: Fn(BigUint, BigUint) -> BigUint,
    {
        let inputs = test_inputs(modulus, word_bits);

        for i in 0..inputs.len() {
            // Iterator over inputs rotated right by i places. Since
            // cycle().skip(i) rotates left by i, we need to rotate by
            // n_input_elts - i.
            let shifted_inputs = inputs.iter().cycle().skip(inputs.len() - i);
            // Calculate pointwise operations
            let expected = inputs
                .iter()
                .zip(shifted_inputs.clone())
                .map(|(x, y)| expected_op(x.clone(), y.clone()));
            let output = inputs.iter().zip(shifted_inputs).map(|(x, y)| {
                field_to_biguint(op(biguint_to_field(x.clone()), biguint_to_field(y.clone())))
            });
            // Compare expected outputs with actual outputs
            assert!(
                output.zip(expected).all(|(x, y)| x == y),
                "output differs from expected at rotation {}",
                i
            );
        }
        Ok(())
    }
}

#[macro_export]
macro_rules! test_arithmetic {
    ($field:ty) => {
        mod arithmetic {
            use crate::{biguint_to_field, field_tests, field_to_biguint, Field};

            use num::{BigUint, Zero};
            use std::io::Result;
            use std::ops::{Add, Div, Mul, Neg, Sub};

            /// Return the modulus of the type `Fld`.
            fn field_modulus<F: Field>() -> BigUint {
                field_to_biguint(F::NEG_ONE) + 1u32
            }

            // Can be 32 or 64; doesn't have to be computer's actual word
            // bits. Choosing 32 gives more tests...
            const WORD_BITS: usize = 32;

            #[test]
            fn arithmetic_addition() -> Result<()> {
                let modulus = field_modulus::<$field>();
                field_tests::run_binaryop_test_cases(&modulus, WORD_BITS, <$field>::add, |x, y| {
                    let z = x + y;
                    if z < modulus {
                        z
                    } else {
                        z - &modulus
                    }
                })
            }

            #[test]
            fn arithmetic_subtraction() -> Result<()> {
                let modulus = field_modulus::<$field>();
                field_tests::run_binaryop_test_cases(&modulus, WORD_BITS, <$field>::sub, |x, y| {
                    if x >= y {
                        x - y
                    } else {
                        &modulus - y + x
                    }
                })
            }

            #[test]
            fn arithmetic_negation() -> Result<()> {
                let modulus = field_modulus::<$field>();
                field_tests::run_unaryop_test_cases(&modulus, WORD_BITS, <$field>::neg, |x| {
                    if x.is_zero() {
                        BigUint::zero()
                    } else {
                        &modulus - x
                    }
                })
            }

            #[test]
            fn arithmetic_multiplication() -> Result<()> {
                let modulus = field_modulus::<$field>();
                field_tests::run_binaryop_test_cases(&modulus, WORD_BITS, <$field>::mul, |x, y| {
                    x * y % &modulus
                })
            }

            #[test]
            fn arithmetic_square() -> Result<()> {
                let modulus = field_modulus::<$field>();
                field_tests::run_unaryop_test_cases(
                    &modulus, WORD_BITS,
                    |x| <$field>::square(&x),
                    |x| &x * &x % &modulus)
            }

            #[test]
            #[ignore]
            fn arithmetic_division() -> Result<()> {
                // This test takes ages to finish so is #[ignore]d by default.
                // TODO: Re-enable and reimplement when
                // https://github.com/rust-num/num-bigint/issues/60 is finally resolved.
                let modulus = field_modulus::<$field>();
                field_tests::run_binaryop_test_cases(
                    &modulus,
                    WORD_BITS,
                    // Need to help the compiler infer the type of y here
                    |x: $field, y: $field| {
                        // TODO: Work out how to check that div() panics
                        // appropriately when given a zero divisor.
                        if !y.is_zero() {
                            <$field>::div(x, y)
                        } else {
                            <$field>::ZERO
                        }
                    },
                    |x, y| {
                        // yinv = y^-1 (mod modulus)
                        let exp = &modulus - BigUint::from(2u32);
                        let yinv = y.modpow(&exp, &modulus);
                        // returns 0 if y was 0
                        x * yinv % &modulus
                    },
                )
            }

            #[test]
            fn square_root() -> Result<()> {
                // We don't use run_{unary,binary}op_test_cases here because
                // we're just testing 'internal consistency'.
                // TODO: Test that it fails appropriately for non-quadratic
                // residues.
                // NB: Could calculate modular sqrt with BigUint with
                // x^{(p+1)/2} (mod p) but this will be slow as with modular
                // inverse above.
                let modulus = field_modulus::<$field>();
                const WORD_BITS: usize = 32;
                let inputs = field_tests::test_inputs(&modulus, WORD_BITS)
                    .iter()
                    .map(|x| biguint_to_field::<$field>(x.clone()))
                    .collect::<Vec<_>>();
                let squares = inputs.iter().map(|&x| x.square()).collect::<Vec<_>>();
                let roots = squares
                    .iter()
                    .map(|xx| xx.square_root().unwrap())
                    .collect::<Vec<_>>();
                assert!(roots.iter().zip(inputs).all(|(&x, y)| x == y || x == -y));
                Ok(())
            }

            /// Basic GCD algo; replace with Binary GCD for better performance.
            fn gcd(x: BigUint, y: u32) -> u32 {
                // TODO: This function probably belongs somewhere else, but at
                // present is only used in kth_root_consistent_with_exp below.

                use num::ToPrimitive;

                let r = x % y;
                if r.is_zero() {
                    y
                } else {
                    let mut a = y;
                    let mut b = r.to_u32().unwrap();
                    while b != 0u32 {
                        let t = a % b;
                        a = b;
                        b = t;
                    }
                    a
                }
            }

            #[test]
            fn kth_root_consistent_with_exp() {
                // We only test degrees that are coprime q-1 as these are
                // the ones that give rise to permutations.
                let degs = (3u32..101).step_by(2)
                    .filter(|&x| gcd(field_to_biguint(<$field>::NEG_ONE), x) == 1u32);
                for deg in degs {
                    let num = <$field>::rand();
                    assert_eq!(num, num.exp_u32(deg).kth_root_u32(deg));
                }
            }
        }
    };
}
