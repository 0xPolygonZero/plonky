use crate::{CircuitBuilder, Field, HaloCurve, Target};
use std::collections::HashMap;
use std::iter::FromIterator;

/// Struct holding a multivariate polynomial \sum c_e.x^e as the map `e => c`.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct MultivariatePolynomial<F: Field, const N: usize>(HashMap<[usize; N], F>);

impl<F: Field, const N: usize> MultivariatePolynomial<F, N> {
    /// Takes a slice of coefficients and returns the corresponding polynomial.
    pub fn from_coeffs(coeffs: &[([usize; N], F)]) -> Self {
        Self(HashMap::from_iter(coeffs.iter().copied()))
    }

    /// Zero polynomial with length `len`.
    /// `len = 1` is the standard representation, but sometimes it's useful to set `len > 1`
    /// to have polynomials with uniform length.
    pub fn zero() -> Self {
        Self(HashMap::from_iter(vec![([0; N], F::ZERO)]))
    }

    pub fn iter(&self) -> std::collections::hash_map::Iter<[usize; N], F> {
        self.0.iter()
    }

    pub fn iter_mut(&mut self) -> std::collections::hash_map::IterMut<[usize; N], F> {
        self.0.iter_mut()
    }

    pub fn is_zero(&self) -> bool {
        self.0.values().all(|c| c.is_zero())
    }

    /// Number of coefficients held by the polynomial.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Degree of the polynomial. Panics if polynomial is empty.
    pub fn degree(&self) -> usize {
        self.0.keys().map(|e| e.iter().sum()).max().unwrap()
    }

    /// Degree of the polynomial + 1.
    #[allow(dead_code)]
    fn degree_plus_one(&self) -> usize {
        self.0
            .keys()
            .map(|e| e.iter().sum::<usize>() + 1)
            .max()
            .unwrap_or(0)
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Evaluates the polynomial at a point `x`.
    pub fn eval(&self, xs: [F; N]) -> F {
        self.iter().fold(F::ZERO, |acc, (es, &c)| {
            acc + c * xs
                .iter()
                .zip(es)
                // Evaluate monomial.
                .fold(F::ONE, |acc, (x, &e)| acc * x.exp_usize(e))
        })
    }

    /// Evaluates the polynomial at a target `t`.
    /// TODO: This isn't optimal at all. A multivariate Horner's method would be better.
    pub fn eval_recursively<C: HaloCurve<ScalarField = F>>(
        &self,
        builder: &mut CircuitBuilder<C>,
        ts: [Target<F>; N],
    ) -> Target<F> {
        let zero = builder.constant_wire(F::ZERO);
        self.iter().fold(zero, |acc, (es, &c)| {
            let terms = ts
                .iter()
                .zip(es)
                .map(|(&t, &e)| builder.exp_constant(t, F::from_canonical_usize(e)))
                .collect::<Vec<_>>();
            let monomial = builder.mul_many(&terms);
            let c_target = builder.constant_wire(c);
            let scaled_monomial = builder.mul(c_target, monomial);
            builder.add(acc, scaled_monomial)
        })
    }

    /// Leading coefficient.
    pub fn lead(&self) -> F {
        self.iter()
            .max_by(|x, y| x.0.iter().sum::<usize>().cmp(&y.0.iter().sum::<usize>()))
            .map(|(_e, &c)| c)
            .unwrap_or(F::ZERO)
    }

    /// Negates the polynomial's coefficients.
    pub fn neg(&self) -> Self {
        Self::from_coeffs(&self.iter().map(|(&e, &x)| (e, -x)).collect::<Vec<_>>())
    }

    /// Multiply the polynomial's coefficients by a scalar.
    #[allow(dead_code)]
    fn scalar_mul(&self, c: F) -> Self {
        Self::from_coeffs(&self.iter().map(|(&e, &x)| (e, c * x)).collect::<Vec<_>>())
    }

    /// Removes all zero coefficients.
    pub fn trim(&mut self) {
        self.0.retain(|_e, c| c.is_nonzero());
    }

    /// Polynomial addition.
    pub fn add(&self, other: &Self) -> Self {
        let mut coeffs = self.0.clone();
        other.iter().for_each(|(&e, &c)| {
            let v = coeffs.entry(e).or_insert(F::ZERO);
            *v = *v + c;
        });
        Self(coeffs)
    }

    /// Polynomial substraction.
    pub fn sub(&self, other: &Self) -> Self {
        let mut coeffs = self.0.clone();
        other.iter().for_each(|(&e, &c)| {
            let v = coeffs.entry(e).or_insert(F::ZERO);
            *v = *v - c;
        });
        Self(coeffs)
    }

    /// Polynomial multiplication.
    pub fn mul(&self, other: &Self) -> Self {
        fn add_coeffs<const N: usize>(e1: &[usize; N], e2: &[usize; N]) -> [usize; N] {
            let mut res = [0; N];
            for i in 0..N {
                res[i] = e1[i] + e2[i];
            }
            res
        }

        let mut res = HashMap::new();
        for (e1, &c1) in self.iter() {
            for (e2, &c2) in other.iter() {
                let v = res.entry(add_coeffs(e1, e2)).or_insert(F::ZERO);
                *v = *v + c1 * c2;
            }
        }

        Self(res)
    }

    pub fn exp(&self, power: F) -> Self {
        let power_bits = power.num_bits();
        let mut current = self.clone();
        let mut product = Self::constant(F::ONE);

        for (i, limb) in power.to_canonical_u64_vec().iter().enumerate() {
            for j in 0..64 {
                // If we've gone through all the 1 bits already, no need to keep squaring.
                let bit_index = i * 64 + j;
                if bit_index == power_bits {
                    return product;
                }

                if (limb >> j & 1) != 0 {
                    product = product.mul(&current);
                }
                current = current.mul(&current);
            }
        }

        product
    }

    /// Returns `true` iff the polynomial uses the `i`-th variable.
    pub fn has_var(&self, i: usize) -> bool {
        self.0.iter().any(|(e, c)| (e[i] > 0) && c.is_nonzero())
    }

    /// Returns the `i`-th variable, i.e., the polynomial `X_i`.
    pub fn var(i: usize) -> Self {
        let mut e = [0; N];
        e[i] = 1;
        Self::from_coeffs(&[(e, F::ONE)])
    }

    /// Returns the constant polynomial with coefficient `c`.
    pub fn constant(c: F) -> Self {
        Self::from_coeffs(&[([0; N], c)])
    }
}

// TODO: Add tests
