use crate::{Field, Curve, CircuitBuilder, Target};

/// Evaluate the polynomial which vanishes on any multiplicative subgroup of a given order `n`.
pub(crate) fn eval_zero_poly<F: Field>(n: usize, x: F) -> F {
    // Z(x) = x^n - 1
    x.exp_usize(n) - F::ONE
}

/// Evaluate the Lagrange basis `L_1` with `L_1(1) = 1`, and `L_1(x) = 0` for other members of an
/// order `n` multiplicative subgroup.
pub(crate) fn eval_l_1<F: Field>(n: usize, x: F) -> F {
    // L_1(x) = (x^n - 1) / (n * (x - 1))
    //        = Z(x) / (n * (x - 1))
    eval_zero_poly(n, x) / (F::from_canonical_usize(n) * (x - F::ONE))
}

/// Computes a sum of terms weighted by powers of alpha.
pub(crate) fn reduce_with_powers<F: Field>(
    terms: &[F],
    alpha: F,
) -> F {
    let mut sum = F::ZERO;
    for &term in terms.iter().rev() {
        sum = sum * alpha + term;
    }
    sum
}

/// Computes a sum of terms weighted by powers of alpha.
pub(crate) fn reduce_with_powers_recursive<C: Curve>(
    builder: &mut CircuitBuilder<C>,
    terms: &[Target],
    alpha: Target,
) -> Target {
    let mut sum = builder.zero_wire();
    for &term in terms.iter().rev() {
        sum = builder.mul_add(sum, alpha, term);
    }
    sum
}

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

pub fn polynomial_division<F: Field>(a: &Polynomial<F>, b: &Polynomial<F>) -> (Vec<F>, Vec<F>) {
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
            let cur_q_degree = degree(&remainder) - b_degree;
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

#[cfg(test)]
mod test {
    use super::*;
    use crate::{Bls12377Scalar, Field};

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
    fn test_polynomial_division() {
        type F = Bls12377Scalar;
        let a: Vec<F> = (0..200).map(|_| F::rand()).collect();
        let b: Vec<F> = (0..100).map(|_| F::rand()).collect();
        let (q, r) = polynomial_division(&a, &b);
        for _i in 0..1000 {
            let x = F::rand();
            assert_eq!(evaluate_at_naive(&a, x), evaluate_at_naive(&b, x)*evaluate_at_naive(&q, x) + evaluate_at_naive(&r, x));
        }
    }
}