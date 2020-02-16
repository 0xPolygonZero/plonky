use std::cmp::Ordering;
use std::mem;
use std::ops::{Neg, Shr};

use num::Integer;
use num::traits::NumRef;
use num::traits::RefNum;

/// See https://github.com/rust-num/num-integer/pull/10/files

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct GcdResult<T> {
    /// Greatest common divisor.
    pub gcd: T,
    /// Coefficients such that: gcd(a, b) = c1*a + c2*b
    pub c1: T, pub c2: T,
}

/// Calculate greatest common divisor and the corresponding coefficients.
pub fn extended_gcd<T: Integer + NumRef>(a: T, b: T) -> GcdResult<T>
    where for<'a> &'a T: RefNum<T>
{
    // Euclid's extended algorithm
    let (mut s, mut old_s) = (T::zero(), T::one());
    let (mut t, mut old_t) = (T::one(), T::zero());
    let (mut r, mut old_r) = (b, a);

    while r != T::zero() {
        let quotient = &old_r / &r;
        old_r = old_r - &quotient * &r; mem::swap(&mut old_r, &mut r);
        old_s = old_s - &quotient * &s; mem::swap(&mut old_s, &mut s);
        old_t = old_t - quotient * &t; mem::swap(&mut old_t, &mut t);
    }

    let _quotients = (t, s); // == (a, b) / gcd

    GcdResult { gcd: old_r, c1: old_s, c2: old_t }
}

/// Find the standard representation of a (mod n).
pub fn normalize<T: Integer + NumRef>(a: T, n: &T) -> T {
    let a = a % n;
    match a.cmp(&T::zero()) {
        Ordering::Less => a + n,
        _ => a,
    }
}

/// Calculate the inverse of a (mod n).
pub fn inverse<T: Integer + NumRef + Clone>(a: T, n: &T) -> Option<T>
    where for<'a> &'a T: RefNum<T>
{
    let GcdResult { gcd, c1: c, .. } = extended_gcd(a, n.clone());
    if gcd == T::one() {
        Some(normalize(c, n))
    } else {
        None
    }
}

/// Calculate base^exp (mod modulus).
pub fn powm<T>(base: &T, exp: &T, modulus: &T) -> T
    where T: Integer + NumRef + Clone + Neg<Output = T> + Shr<i32, Output = T>,
          for<'a> &'a T: RefNum<T>
{
    let zero = T::zero();
    let one = T::one();
    let two = &one + &one;
    let mut exp = exp.clone();
    let mut result = one.clone();
    let mut base = base % modulus;
    if exp < zero {
        exp = -exp;
        base = inverse(base, modulus).unwrap();
    }
    while exp > zero {
        if &exp % &two == one {
            result = (result * &base) % modulus;
        }
        exp = exp >> 1;
        base = (&base * &base) % modulus;
    }
    result
}
