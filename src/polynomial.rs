use crate::{
    fft_precompute, fft_with_precomputation, ifft_with_precomputation_power_of_2, util::log2_ceil,
    Field,
};
use std::time::Instant;

// Store polynomial as a list of coefficient starting with the constant coefficient.
type Polynomial<F: Field> = Vec<F>;

fn degree<F: Field>(a: &Polynomial<F>) -> usize {
    let mut ans = 0;
    a.iter().enumerate().for_each(|(i, x)| {
        if x.is_nonzero() {
            ans = i;
        }
    });
    ans
}

fn lead<F: Field>(a: &Polynomial<F>) -> F {
    let mut ans = F::ZERO;
    a.iter().for_each(|&x| {
        if x.is_nonzero() {
            ans = x;
        }
    });
    ans
}

fn is_zero<F: Field>(a: &Polynomial<F>) -> bool {
    a.iter().all(|&x| x.is_zero())
}

fn rev<F: Field>(a: &Polynomial<F>) -> Polynomial<F> {
    assert!(a.last().unwrap().is_nonzero());
    let mut r = a.to_vec();
    r.reverse();
    r
}

fn neg<F: Field>(a: &Polynomial<F>) -> Polynomial<F> {
    a.iter().map(|&x| -x).collect()
}

fn trim<F: Field>(a: &mut Polynomial<F>) {
    let d = degree(a);
    a.drain(d + 1..);
}

pub fn polynomial_addition<F: Field>(a: &Polynomial<F>, b: &Polynomial<F>) -> Polynomial<F> {
    a.iter().zip(b.iter()).map(|(&ca, &cb)| ca + cb).collect()
}

pub fn polynomial_multiplication<F: Field>(a: &Polynomial<F>, b: &Polynomial<F>) -> Polynomial<F> {
    let a_deg = a.len();
    let b_deg = b.len();
    let mut a_pad = a.clone();
    a_pad.append(&mut vec![F::ZERO; b_deg]);
    let mut b_pad = b.clone();
    b_pad.append(&mut vec![F::ZERO; a_deg]);
    let precomputation = fft_precompute(a_deg + b_deg);
    let a_evals = fft_with_precomputation(&a_pad, &precomputation);
    let b_evals = fft_with_precomputation(&b_pad, &precomputation);

    let mul_evals: Vec<F> = a_evals
        .iter()
        .zip(b_evals.iter())
        .map(|(&pa, &pb)| pa * pb)
        .collect();
    ifft_with_precomputation_power_of_2(&mul_evals, &precomputation)
}

pub fn polynomial_long_division<F: Field>(
    a: &Polynomial<F>,
    b: &Polynomial<F>,
) -> (Polynomial<F>, Polynomial<F>) {
    let (a_degree, b_degree) = (degree(a), degree(b));
    if is_zero(a) {
        (Vec::new(), Vec::new())
    } else if is_zero(b) {
        panic!("Division by zero polynomial");
    } else if a_degree < b_degree {
        (Vec::new(), a.to_vec())
    } else {
        // Now we know that self.degree() >= divisor.degree();
        let mut quotient = vec![F::ZERO; a_degree - b_degree + 1];
        let mut remainder = a.to_vec();
        // Can unwrap here because we know self is not zero.
        let divisor_leading_inv = lead(b).multiplicative_inverse_assuming_nonzero();
        while !is_zero(&remainder) && degree(&remainder) >= b_degree {
            let cur_q_coeff = *remainder.last().unwrap() * divisor_leading_inv;
            let cur_q_degree = remainder.len() - 1 - b_degree;
            quotient[cur_q_degree] = cur_q_coeff;

            for (i, &div_coeff) in b.iter().enumerate() {
                remainder[cur_q_degree + i] =
                    remainder[cur_q_degree + i] - (cur_q_coeff * div_coeff);
            }
            while let Some(true) = remainder.last().map(|c| c.is_zero()) {
                remainder.pop();
            }
        }
        (quotient, remainder)
    }
}

fn inv_mod_xn<F: Field>(h: &Polynomial<F>, n: usize) -> Polynomial<F> {
    assert!(h[0].is_nonzero(), "Inverse doesn't exist.");
    let mut hh = h.clone();
    if h.len() < n {
        for _ in 0..n - h.len() {
            hh.push(F::ZERO);
        }
    }
    let mut a: Polynomial<F> = Vec::new();
    a.push(hh[0].multiplicative_inverse_assuming_nonzero());
    for i in 0..log2_ceil(n) {
        let l = 1 << i;
        let mut h0 = hh[..l].to_vec();
        let mut h1 = hh[l..].to_vec();
        let mut c = polynomial_multiplication(&a, &h0);
        c.drain(0..l);
        trim(&mut h1);
        let mut tmp = polynomial_multiplication(&h1, &a);
        tmp = polynomial_addition(&tmp, &c);
        tmp.iter_mut().for_each(|x| *x = -(*x));
        trim(&mut tmp);
        let mut b = polynomial_multiplication(&a, &tmp)[..l].to_vec();
        a.append(&mut b);
    }
    a.drain(n..);
    a
}

// Algorithm from http://people.csail.mit.edu/madhu/ST12/scribe/lect06.pdf
pub fn polynomial_division<F: Field>(
    a: &Polynomial<F>,
    b: &Polynomial<F>,
) -> (Polynomial<F>, Polynomial<F>) {
    let deg_a = degree(a);
    let deg_b = degree(b);
    let rev_b = rev(b);
    let rev_b_inv = inv_mod_xn(&rev_b, deg_a - deg_b + 1);
    let rev_q = polynomial_multiplication(&rev_b_inv, &rev(a)[..=deg_a - deg_b].to_vec())
        [..=deg_a - deg_b]
        .to_vec();
    let q = rev(&rev_q);
    let qb = polynomial_multiplication(&q, &b);
    let r = polynomial_addition(&a, &neg(&qb));
    (q, r)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{Field, TweedledeeBase};

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
        let a: Vec<F> = (0..128000).map(|_| F::rand()).collect();
        let b: Vec<F> = (0..16000).map(|_| F::rand()).collect();
        let m = polynomial_multiplication(&a, &b);
        for _ in 0..1000 {
            let x = F::rand();
            assert_eq!(
                evaluate_at_naive(&m, x),
                evaluate_at_naive(&a, x) * evaluate_at_naive(&b, x)
            );
        }
    }

    #[test]
    fn test_inv_mod_xn() {
        type F = TweedledeeBase;
        let a: Vec<F> = (0..128000).map(|_| F::rand()).collect();
        let b = inv_mod_xn(&a, 200_000);
        let m = polynomial_multiplication(&a, &b);
        assert_eq!(m[0], F::ONE);
        m[1..200_000].iter().for_each(|&x| {
            assert_eq!(x, F::ZERO);
        });
    }

    #[test]
    fn test_polynomial_long_division() {
        type F = TweedledeeBase;
        let a: Vec<F> = (0..12800).map(|_| F::rand()).collect();
        let b: Vec<F> = (0..1600).map(|_| F::rand()).collect();
        let (q, r) = polynomial_long_division(&a, &b);
        for _ in 0..1000 {
            let x = F::rand();
            assert_eq!(
                evaluate_at_naive(&a, x),
                evaluate_at_naive(&b, x) * evaluate_at_naive(&q, x) + evaluate_at_naive(&r, x)
            );
        }
    }

    #[test]
    fn test_polynomial_division_small() {
        type F = TweedledeeBase;
        let a: Vec<F> = (0..12_000).map(|_| F::rand()).collect();
        let b: Vec<F> = (0..1_000).map(|_| F::rand()).collect();
        let (mut q, mut r) = polynomial_division(&a, &b);
        assert!(degree(&r) < degree(&b));
        let (mut ql, mut rl) = polynomial_long_division(&a, &b);
        trim(&mut q);
        trim(&mut r);
        trim(&mut ql);
        trim(&mut rl);
        assert_eq!(&q, &ql);
        assert_eq!(&r, &rl);
        for _ in 0..1000 {
            let x = F::rand();
            assert_eq!(
                evaluate_at_naive(&a, x),
                evaluate_at_naive(&b, x) * evaluate_at_naive(&q, x) + evaluate_at_naive(&r, x)
            );
        }
    }

    #[test]
    fn test_polynomial_division_large() {
        type F = TweedledeeBase;
        let a: Vec<F> = (0..128_000).map(|_| F::rand()).collect();
        let b: Vec<F> = (0..16_000).map(|_| F::rand()).collect();
        let now = Instant::now();
        let (q, r) = polynomial_division(&a, &b);
        println!("Division time: {:?}", now.elapsed());
        assert!(degree(&r) < degree(&b));
        for _ in 0..1000 {
            let x = F::rand();
            assert_eq!(
                evaluate_at_naive(&a, x),
                evaluate_at_naive(&b, x) * evaluate_at_naive(&q, x) + evaluate_at_naive(&r, x)
            );
        }
    }
}
