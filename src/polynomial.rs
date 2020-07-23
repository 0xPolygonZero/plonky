use crate::{fft_precompute, fft_with_precomputation, ifft_with_precomputation_power_of_2, util::log2_ceil, FftPrecomputation, Field};
use std::ops::{Index, IndexMut, RangeBounds};
use std::slice::{Iter, IterMut, SliceIndex};

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Polynomial<F: Field>(Vec<F>);

impl<F: Field> From<Vec<F>> for Polynomial<F> {
    fn from(coeffs: Vec<F>) -> Self {
        Self(coeffs)
    }
}

impl<F, I> Index<I> for Polynomial<F>
where
    F: Field,
    I: SliceIndex<[F]>,
{
    type Output = I::Output;

    fn index(&self, index: I) -> &Self::Output {
        &self.0[index]
    }
}

impl<F, I> IndexMut<I> for Polynomial<F>
where
    F: Field,
    I: SliceIndex<[F]>,
{
    fn index_mut(&mut self, index: I) -> &mut <Self as Index<I>>::Output {
        &mut self.0[index]
    }
}

impl<F: Field> Polynomial<F> {
    fn from_coeffs(coeffs: &[F]) -> Self {
        Self(coeffs.to_vec())
    }

    fn coeffs(&self) -> Vec<F> {
        self.0.clone()
    }

    fn empty() -> Self {
        Self(Vec::new())
    }

    fn zero(len: usize) -> Self {
        Self(vec![F::ZERO; len])
    }

    fn iter(&self) -> Iter<F> {
        self.0.iter()
    }

    fn iter_mut(&mut self) -> IterMut<F> {
        self.0.iter_mut()
    }

    fn is_zero(&self) -> bool {
        self.0.iter().all(|x| x.is_zero())
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn drain<R: RangeBounds<usize>>(&mut self, range: R) {
        self.0.drain(range);
    }

    fn eval(&self, x: F) -> F {
        self.iter().rev().fold(F::ZERO, |acc, &c| acc * x + c)
    }

    fn eval_domain(&self, fft_precomputation: &FftPrecomputation<F>) -> Vec<F> {
        fft_with_precomputation(&self.coeffs(), fft_precomputation)
    }

    fn from_evaluations(values: &[F], fft_precomputation: &FftPrecomputation<F>) -> Self {
        Self(ifft_with_precomputation_power_of_2(
            values,
            fft_precomputation,
        ))
    }

    fn degree(&self) -> usize {
        (0usize..self.0.len())
            .rev()
            .find(|&i| self.0[i].is_nonzero())
            .expect("Zero polynomial")
    }

    fn lead(&self) -> F {
        self.iter()
            .rev()
            .find(|x| x.is_nonzero())
            .map_or(F::ZERO, |x| *x)
    }

    fn rev(&self) -> Self {
        let d = self.degree();
        Self(self.0[..=d].iter().rev().copied().collect())
    }

    fn neg(&self) -> Self {
        Self(self.iter().map(|&x| -x).collect())
    }

    fn trim(&mut self) {
        if !self.is_zero() {
            self.0.drain(self.degree() + 1..);
        }
    }

    fn add(&self, other: &Self) -> Self {
        let (mut a, mut b) = (self.clone(), other.clone());
        if a.len() < b.len() {
            a.pad(b.len())
        } else if b.len() < a.len() {
            b.pad(a.len());
        }
        Self(a.iter().zip(b.iter()).map(|(&x, &y)| x + y).collect())
    }

    fn pad(&mut self, len: usize) {
        assert!(self.len() <= len);
        self.0.extend((self.len()..len).map(|_| F::ZERO));
    }

    fn mul(&self, b: &Self) -> Self {
        if self.is_zero() || b.is_zero() {
            return Self::zero(1);
        }
        let a_deg = self.degree();
        let b_deg = b.degree();
        let mut a_pad = self.clone();
        a_pad.pad(a_deg + b_deg + 1);
        let mut b_pad = b.clone();
        b_pad.pad(a_deg + b_deg + 1);

        let precomputation = fft_precompute(a_deg + b_deg + 1);
        let a_evals = fft_with_precomputation(&a_pad.0, &precomputation);
        let b_evals = fft_with_precomputation(&b_pad.0, &precomputation);

        let mul_evals: Vec<F> = a_evals
            .iter()
            .zip(b_evals.iter())
            .map(|(&pa, &pb)| pa * pb)
            .collect();
        ifft_with_precomputation_power_of_2(&mul_evals, &precomputation).into()
    }

    pub fn polynomial_long_division(&self, b: &Self) -> (Self, Self) {
        let (a_degree, b_degree) = (self.degree(), b.degree());
        if self.is_zero() {
            (Self::zero(1), Self::empty())
        } else if b.is_zero() {
            panic!("Division by zero polynomial");
        } else if a_degree < b_degree {
            (Self::zero(1), self.clone())
        } else {
            // Now we know that self.degree() >= divisor.degree();
            let mut quotient = Self::zero(a_degree - b_degree + 1);
            let mut remainder = self.clone();
            // Can unwrap here because we know self is not zero.
            let divisor_leading_inv = b.lead().multiplicative_inverse_assuming_nonzero();
            while !remainder.is_zero() && remainder.degree() >= b_degree {
                let cur_q_coeff = remainder.lead() * divisor_leading_inv;
                let cur_q_degree = remainder.degree() - b_degree;
                quotient[cur_q_degree] = cur_q_coeff;

                for (i, &div_coeff) in b.iter().enumerate() {
                    remainder[cur_q_degree + i] =
                        remainder[cur_q_degree + i] - (cur_q_coeff * div_coeff);
                }
                remainder.trim();
            }
            (quotient, remainder)
        }
    }

    /// Computes the inverse of `self` modulo `x^n`.
    fn inv_mod_xn(&self, n: usize) -> Self {
        assert!(self[0].is_nonzero(), "Inverse doesn't exist.");
        let mut h = self.clone();
        if h.len() < n {
            h.pad(n);
        }
        let mut a = Self::empty();
        a.0.push(h[0].multiplicative_inverse_assuming_nonzero());
        for i in 0..log2_ceil(n) {
            let l = 1 << i;
            let h0 = h[..l].to_vec().into();
            let mut h1: Polynomial<F> = h[l..].to_vec().into();
            let mut c = a.mul(&h0);
            if l == c.len() {
                c = Self::zero(1);
            } else {
                c.drain(0..l);
            }
            h1.trim();
            let mut tmp = a.mul(&h1);
            tmp = tmp.add(&c);
            tmp.iter_mut().for_each(|x| *x = -(*x));
            tmp.trim();
            let b = &a.mul(&tmp)[..l];
            a.0.extend_from_slice(b);
        }
        a.drain(n..);
        a
    }

    /// Returns `(q,r)` the quotient and remainder of the polynomial division of `a` by `b`.
    /// Algorithm from http://people.csail.mit.edu/madhu/ST12/scribe/lect06.pdf
    pub fn polynomial_division(&self, b: &Self) -> (Self, Self) {
        let (a_degree, b_degree) = (self.degree(), b.degree());
        if self.is_zero() {
            (Self::zero(1), Self::empty())
        } else if b.is_zero() {
            panic!("Division by zero polynomial");
        } else if a_degree < b_degree {
            (Self::zero(1), self.clone())
        } else {
            let rev_b = b.rev();
            let rev_b_inv = rev_b.inv_mod_xn(a_degree - b_degree + 1);
            let rev_q: Polynomial<F> = rev_b_inv
                .mul(&self.rev()[..=a_degree - b_degree].to_vec().into())[..=a_degree - b_degree]
                .to_vec()
                .into();
            let mut q = rev_q.rev();
            let mut qb = q.mul(b);
            qb.trim();
            qb.pad(self.len());
            let mut r = self.add(&qb.neg());
            q.trim();
            r.trim();
            (q, r)
        }
    }

    // Divides a polynomial `a` by `Z_H = X^n - 1`. Assumes `Z_H | a`, otherwise result is meaningless.
    pub fn divide_by_z_h(&self, n: usize) -> Self {
        let mut a_trim = self.clone();
        a_trim.trim();
        let g = F::MULTIPLICATIVE_SUBGROUP_GENERATOR;
        let mut g_pow = F::ONE;
        // Multiply the i-th coefficient of `a` by `g^i`. Then `new_a(w^j) = old_a(g.w^j)`.
        a_trim.iter_mut().for_each(|x| {
            *x = (*x) * g_pow;
            g_pow = g * g_pow;
        });
        let d = a_trim.degree();
        let root = F::primitive_root_of_unity(log2_ceil(d + 1));
        let precomputation = fft_precompute(d + 1);
        // Equals to the evaluation of `a` on `{g.w^i}`.
        let mut a_eval = a_trim.eval_domain(&precomputation);
        // Compute the denominators `1/(g^n.w^(n*i) - 1)` using batch inversion.
        let denominator_g = g.exp_usize(n);
        let root_n = root.exp_usize(n);
        let mut root_pow = F::ONE;
        let denominators = (0..a_eval.len())
            .map(|i| {
                if i != 0 {
                    root_pow = root_pow * root_n;
                }
                denominator_g * root_pow - F::ONE
            })
            .collect::<Vec<_>>();
        let denominators_inv = F::batch_multiplicative_inverse(&denominators);
        // Divide every element of `a_eval` by the corresponding denominator.
        // Then, `a_eval` is the evaluation of `a/Z_H` on `{g.w^i}`.
        a_eval
            .iter_mut()
            .zip(denominators_inv.iter())
            .for_each(|(x, &d)| {
                *x = (*x) * d;
            });
        // `p` is the interpolating polynomial of `a_eval` on `{w^i}`.
        let mut p = Self::from_evaluations(&a_eval, &precomputation);
        // We need to scale it by `g^(-i)` to get the interpolating polynomial of `a_eval` on `{g.w^i}`,
        // a.k.a `a/Z_H`.
        let g_inv = g.multiplicative_inverse_assuming_nonzero();
        let mut g_inv_pow = F::ONE;
        p.iter_mut().for_each(|x| {
            *x = (*x) * g_inv_pow;
            g_inv_pow = g_inv_pow * g_inv;
        });
        p
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{Field, TweedledeeBase};
    use rand::{thread_rng, Rng};
    use std::time::Instant;

    fn evaluate_at_naive<F: Field>(coefficients: &[F], point: F) -> F {
        let mut sum = F::ZERO;
        let mut point_power = F::ONE;
        for &c in coefficients {
            sum = sum + c * point_power;
            point_power = point_power * point;
        }
        sum
    }

    #[test]
    fn test_polynomial_multiplication() {
        type F = TweedledeeBase;
        let mut rng = thread_rng();
        let (a_deg, b_deg) = (rng.gen_range(1, 10_000), rng.gen_range(1, 10_000));
        let a = Polynomial((0..a_deg).map(|_| F::rand()).collect());
        let b = Polynomial((0..b_deg).map(|_| F::rand()).collect());
        let m1 = a.mul(&b);
        let m2 = a.mul(&b);
        for _ in 0..1000 {
            let x = F::rand();
            assert_eq!(m1.eval(x), a.eval(x) * b.eval(x));
            assert_eq!(m2.eval(x), a.eval(x) * b.eval(x));
        }
    }

    #[test]
    fn test_inv_mod_xn() {
        type F = TweedledeeBase;
        let mut rng = thread_rng();
        let a_deg = rng.gen_range(1, 1_000);
        let n = rng.gen_range(1, 1_000);
        let a = Polynomial((0..a_deg).map(|_| F::rand()).collect());
        let b = a.inv_mod_xn(n);
        let mut m = a.mul(&b);
        m.drain(n..);
        m.trim();
        assert_eq!(m, Polynomial(vec![F::ONE]));
    }

    #[test]
    fn test_polynomial_long_division() {
        type F = TweedledeeBase;
        let mut rng = thread_rng();
        let (a_deg, b_deg) = (rng.gen_range(1, 10_000), rng.gen_range(1, 10_000));
        let a = Polynomial((0..a_deg).map(|_| F::rand()).collect());
        let b = Polynomial((0..b_deg).map(|_| F::rand()).collect());
        let (q, r) = a.polynomial_long_division(&b);
        for _ in 0..1000 {
            let x = F::rand();
            assert_eq!(a.eval(x), b.eval(x) * q.eval(x) + r.eval(x));
        }
    }

    #[test]
    fn test_polynomial_division() {
        type F = TweedledeeBase;
        let mut rng = thread_rng();
        let (a_deg, b_deg) = (rng.gen_range(1, 10_000), rng.gen_range(1, 10_000));
        let a = Polynomial((0..a_deg).map(|_| F::rand()).collect());
        let b = Polynomial((0..b_deg).map(|_| F::rand()).collect());
        let (q, r) = a.polynomial_division(&b);
        for _ in 0..1000 {
            let x = F::rand();
            assert_eq!(a.eval(x), b.eval(x) * q.eval(x) + r.eval(x));
        }
    }

    #[test]
    fn test_division_by_z_h() {
        type F = TweedledeeBase;
        let mut rng = thread_rng();
        let a_deg = rng.gen_range(1, 10_000);
        let n = rng.gen_range(1, a_deg);
        let mut a = Polynomial((0..a_deg).map(|_| F::rand()).collect());
        a.trim();
        let z_h = {
            let mut z_h_vec = vec![F::ZERO; n + 1];
            z_h_vec[n] = F::ONE;
            z_h_vec[0] = F::NEG_ONE;
            Polynomial(z_h_vec)
        };
        let m = a.mul(&z_h);
        let now = Instant::now();
        let mut a_test = m.divide_by_z_h(n);
        a_test.trim();
        println!("Division time: {:?}", now.elapsed());
        assert_eq!(a, a_test);
    }

    // Test to see which polynomial division method is faster for divisions of the type
    // `(X^n - 1)/(X - a)
    #[test]
    fn test_division_linear() {
        type F = TweedledeeBase;
        let mut rng = thread_rng();
        let l = 14;
        let n = 1 << l;
        let g = F::primitive_root_of_unity(l);
        let xn_minus_one = {
            let mut xn_min_one_vec = vec![F::ZERO; n + 1];
            xn_min_one_vec[n] = F::ONE;
            xn_min_one_vec[0] = F::NEG_ONE;
            Polynomial(xn_min_one_vec)
        };

        let a = g.exp_usize(rng.gen_range(0, n));
        let denom = Polynomial(vec![-a, F::ONE]);
        let now = Instant::now();
        xn_minus_one.polynomial_division(&denom);
        println!("Division time: {:?}", now.elapsed());
        let now = Instant::now();
        xn_minus_one.polynomial_long_division(&denom);
        println!("Division time: {:?}", now.elapsed());
    }
}
