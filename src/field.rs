use std::collections::HashSet;
use std::hash::Hash;
use std::ops::{Add, Div, Mul, Neg, Sub};

pub trait Field: Sized + Copy + Eq + Send + Sync
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

    fn to_canonical_vec(&self) -> Vec<u64>;

    fn from_canonical_vec(v: Vec<u64>) -> Self;

    fn from_canonical_u64(n: u64) -> Self;

    fn from_canonical_usize(n: usize) -> Self {
        Self::from_canonical_u64(n as u64)
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        *self == Self::ZERO
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

    fn cyclic_subgroup_unknown_order(generator: Self) -> Vec<Self> where Self: Hash {
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

    fn rand() -> Self;
}

pub trait TwoAdicField: Field {
    const TWO_ADICITY: usize;

    /// Computes a `2^n_power`th primitive root of unity.
    fn primitive_root_of_unity(n_power: usize) -> Self;
}
