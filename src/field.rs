use std::collections::HashSet;
use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{Add, Div, Mul, Neg, Shl, Sub};

use num::{BigUint, FromPrimitive, Integer, One};

use crate::{biguint_to_field, field_to_biguint};

pub trait Field: 'static + Sized + Copy + Eq + Hash + Send + Sync + Debug
+ Neg<Output=Self>
+ Add<Self, Output=Self>
+ Sub<Self, Output=Self>
+ Mul<Self, Output=Self>
+ Div<Self, Output=Self> {
    const BITS: usize;

    const ZERO: Self;
    const ONE: Self;
    const TWO: Self;
    const THREE: Self;
    const FOUR: Self;
    const FIVE: Self;
    const NEG_ONE: Self;

    const MULTIPLICATIVE_SUBGROUP_GENERATOR: Self;

    fn to_canonical_u64_vec(&self) -> Vec<u64>;

    fn to_canonical_u32_vec(&self) -> Vec<u32> {
        let mut limbs = Vec::new();
        for u64_limb in self.to_canonical_u64_vec() {
            limbs.push(u64_limb as u32);
            limbs.push(u64_limb.overflowing_shr(32).0 as u32);
        }
        limbs
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
        Self::from_canonical_u64(if b { 1 } else { 0 })
    }

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
        let mut subgroup_vec = Vec::new();
        let mut subgroup_set = HashSet::new();
        let mut current = Self::ONE;
        loop {
            if !subgroup_set.insert(current) {
                break;
            }
            subgroup_vec.push(current);
            current = current * generator;
        }
        subgroup_vec
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

    fn generator_order(generator: Self) -> usize where Self: Hash {
        Self::cyclic_subgroup_unknown_order(generator).len()
    }

    fn exp(&self, power: Self) -> Self {
        let power_bits = power.num_bits();
        let mut current = *self;
        let mut product = Self::ONE;

        for (i, limb) in power.to_canonical_u64_vec().iter().enumerate() {
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

    fn exp_u32(&self, power: u32) -> Self {
        self.exp(Self::from_canonical_u32(power))
    }

    /// Assumes x^k is a permutation in this field; undefined behavior otherwise.
    fn kth_root_u32(&self, k: u32) -> Self {
        // By Fermat's little theorem, x^p = x and x^(p - 1) = 1, so x^(p + n(p - 1)) = x for any n.
        // Our assumption that the k'th root operation is a permutation implies gcd(p - 1, k) = 1,
        // so there exists some n such that p + n(p - 1) is a multiple of k. Once we find such an n,
        // we can rewrite the above as
        //    x^((p + n(p - 1))/k)^k = x,
        // implying that x^((p + n(p - 1))/k) is a k'th root of x.
        let p_minus_1_bu = field_to_biguint(Self::NEG_ONE);
        let p_bu = &p_minus_1_bu + BigUint::one();
        let k_bu = BigUint::from_u32(k).unwrap();
        for n in 0..k {
            let numerator_bu = &p_bu + BigUint::from_u32(n).unwrap() * &p_minus_1_bu;
            if numerator_bu.is_multiple_of(&k_bu) {
                let power_bu = numerator_bu.div_floor(&k_bu).mod_floor(&p_minus_1_bu);
                return self.exp(biguint_to_field(power_bu));
            }
        }
        panic!("x^{} and x^(1/{}) are not permutations in this field, or we have a bug!", k, k);
    }

    fn is_quadratic_residue(&self) -> bool {
        // This is based on Euler's criterion.
        let power = biguint_to_field(field_to_biguint(Self::NEG_ONE).shl(1));
        let exp = self.exp(power);
        if exp == Self::ONE {
            return true;
        }
        if exp == Self::NEG_ONE {
            return false;
        }
        panic!("Number theory is a lie!")
    }

    /// Generate a list of unique quadratic nonresidues in this field. This is guaranteed to be
    /// deterministic, but the behavior is undefined beyond that.
    fn generate_quadratic_nonresidues(n: usize) -> Vec<Self> {
        let mut a = Self::TWO;
        let mut residues = Vec::new();
        while residues.len() < n {
            if !a.is_quadratic_residue() {
                residues.push(a);
            }
            a = a + Self::ONE;
        }
        residues
    }

    fn num_bits(&self) -> usize {
        let mut n = 0;
        for (i, limb) in self.to_canonical_u64_vec().iter().enumerate() {
            for j in 0..64 {
                if (limb >> j & 1) != 0 {
                    n = i * 64 + j + 1;
                }
            }
        }
        n
    }

    fn rand() -> Self;
}

pub trait TwoAdicField: Field {
    const TWO_ADICITY: usize;

    /// `T = (ORDER - 1) / 2^TWO_ADICITY`
    const T: Self;

    /// Computes a `2^n_power`th primitive root of unity.
    fn primitive_root_of_unity(n_power: usize) -> Self {
        assert!(n_power <= Self::TWO_ADICITY);
        let base_root = Self::MULTIPLICATIVE_SUBGROUP_GENERATOR.exp(Self::T);
        base_root.exp(Self::from_canonical_u64(1u64 << (Self::TWO_ADICITY as u64 - n_power as u64)))
    }
}
