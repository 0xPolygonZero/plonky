use crate::Field;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::iter::FromIterator;
use std::ops::{Index, IndexMut, RangeBounds};
use std::slice::{Iter, IterMut, SliceIndex};

/// Struct holding a mulivariate polynomial \sum c_e.x^e as the map `e => c`.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct MultivariatePolynomial<F: Field, const N: usize>(HashMap<[usize; N], F>);

// impl<F: Field, const N: usize> Parti
//     fn eq(&self, other: &Self) -> bool {
//         self.hashmap() == other.hashmap()
//     }
// }

// impl<F: Field, const N: usize> Eq for MultivariatePolynomial<F, N> {}
//
// impl<F, I, const N: usize> Index<I> for MultivariatePolynomial<F, N>
// where
//     F: Field,
//     I: SliceIndex<[(F, [usize; N])]>,
// {
//     type Output = I::Output;
//
//     /// Indexing on the coefficients.
//     fn index(&self, index: I) -> &Self::Output {
//         &self.0[index]
//     }
// }

// impl<F, I, const N: usize> IndexMut<I> for MultivariatePolynomial<F, N>
// where
//     F: Field,
//     I: SliceIndex<[(F, [usize; N])]>,
// {
//     fn index_mut(&mut self, index: I) -> &mut <Self as Index<I>>::Output {
//         &mut self.0[index]
//     }
// }

// impl<F: Field, const N: usize> From<Vec<(F, [usize; N])>> for MultivariatePolynomial<F, N> {
//     /// Takes a vector of coefficients and returns the corresponding polynomial.
//     fn from(coeffs: Vec<(F, [usize; N])>) -> Self {
//         Self(coeffs)
//     }
// }

impl<F: Field, const N: usize> MultivariatePolynomial<F, N> {
    /// Takes a slice of coefficients and returns the corresponding polynomial.
    pub fn from_coeffs(coeffs: &[([usize; N], F)]) -> Self {
        Self(HashMap::from_iter(coeffs.iter().copied()))
    }
    //
    // /// Returns the coefficient vector.
    // pub fn coeffs(&self) -> &[([usize; N], F)] {
    //     &self.0.iter().copied().collect::<Vec<_>>()
    // }
    //
    // /// Empty polynomial;
    // pub fn empty() -> Self {
    //     Self(HashMap::new())
    // }

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

    // /// Evaluates the polynomial at a point `x`, given the list of powers of `x`.
    // /// Assumes that `self.len() == x_pow.len()`.
    // pub fn eval_from_power(&self, x_pow: &[F]) -> F {
    //     F::inner_product(&self[..], x_pow)
    // }

    // /// Evaluates the polynomial on subgroup of `F^*` with a given FFT precomputation.
    // pub fn eval_domain(&self, fft_precomputation: &FftPrecomputation<F>) -> Vec<F> {
    //     let domain_size = fft_precomputation.size();
    //     if self.len() < domain_size {
    //         // Need to pad the polynomial to have the same length as the domain.
    //         fft_with_precomputation(&self.padded(domain_size).coeffs(), fft_precomputation)
    //     } else {
    //         fft_with_precomputation(&self.coeffs(), fft_precomputation)
    //     }
    // }

    // /// Computes the interpolating polynomial of a list of `values` on a subgroup of `F^*`.
    // pub fn from_evaluations(values: &[F], fft_precomputation: &FftPrecomputation<F>) -> Self {
    //     Self(ifft_with_precomputation_power_of_2(
    //         values,
    //         fft_precomputation,
    //     ))
    // }

    /// Leading coefficient.
    pub fn lead(&self) -> F {
        self.iter()
            .max_by(|x, y| x.0.iter().sum::<usize>().cmp(&y.0.iter().sum::<usize>()))
            .map(|(e, &c)| c)
            .unwrap_or(F::ZERO)
    }

    // /// Negates the polynomial's coefficients.
    // fn neg(&self) -> Self {
    //     Self(self.iter().map(|&x| -x).collect())
    // }

    /// Multiply the polynomial's coefficients by a scalar.
    fn scalar_mul(&self, c: F) -> Self {
        Self::from_coeffs(&self.iter().map(|(&e, &x)| (e, c * x)).collect::<Vec<_>>())
    }

    /// Removes all zero coefficients.
    pub fn trim(&mut self) {
        self.0.retain(|e,c| c.is_nonzero());
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
                *v = *v + c1*c2;
            }
        }

        Self(res)
    }
}

// fn compose<F: Field, const M: usize, const N: usize>(p: MultivariatePolynomial<F, M>, polys: [MultivariatePolynomial<F, N>; M]) -> MultivariatePolynomial<F, N> {
//     todo!()
// }
