#![allow(clippy::many_single_char_names)]
use std::cmp::Ordering::Less;

use crate::{add_4_4_no_overflow, add_6_6_no_overflow, cmp_4_4, cmp_6_6, div2_4, div2_6, is_even_4, is_even_6, is_odd_4, is_odd_6, sub_4_4, sub_6_6};

pub(crate) fn nonzero_multiplicative_inverse_4(a: [u64; 4], order: [u64; 4]) -> [u64; 4] {
    // Based on Algorithm 16 of "Efficient Software-Implementation of Finite Fields with
    // Applications to Cryptography".

    let zero = [0, 0, 0, 0];
    let one = [1, 0, 0, 0];

    let mut u = a;
    let mut v = order;
    let mut b = one;
    let mut c = zero;

    while u != one && v != one {
        while is_even_4(u) {
            u = div2_4(u);
            if is_odd_4(b) {
                b = add_4_4_no_overflow(b, order);
            }
            b = div2_4(b);
        }

        while is_even_4(v) {
            v = div2_4(v);
            if is_odd_4(c) {
                c = add_4_4_no_overflow(c, order);
            }
            c = div2_4(c);
        }

        if cmp_4_4(u, v) == Less {
            v = sub_4_4(v, u);
            if cmp_4_4(c, b) == Less {
                c = add_4_4_no_overflow(c, order);
            }
            c = sub_4_4(c, b);
        } else {
            u = sub_4_4(u, v);
            if cmp_4_4(b, c) == Less {
                b = add_4_4_no_overflow(b, order);
            }
            b = sub_4_4(b, c);
        }
    }

    if u == one {
        b
    } else {
        c
    }
}

pub(crate) fn nonzero_multiplicative_inverse_6(a: [u64; 6], order: [u64; 6]) -> [u64; 6] {
    // Based on Algorithm 16 of "Efficient Software-Implementation of Finite Fields with
    // Applications to Cryptography".

    let zero = [0, 0, 0, 0, 0, 0];
    let one = [1, 0, 0, 0, 0, 0];

    let mut u = a;
    let mut v = order;
    let mut b = one;
    let mut c = zero;

    while u != one && v != one {
        while is_even_6(u) {
            u = div2_6(u);
            if is_odd_6(b) {
                b = add_6_6_no_overflow(b, order);
            }
            b = div2_6(b);
        }

        while is_even_6(v) {
            v = div2_6(v);
            if is_odd_6(c) {
                c = add_6_6_no_overflow(c, order);
            }
            c = div2_6(c);
        }

        if cmp_6_6(u, v) == Less {
            v = sub_6_6(v, u);
            if cmp_6_6(c, b) == Less {
                c = add_6_6_no_overflow(c, order);
            }
            c = sub_6_6(c, b);
        } else {
            u = sub_6_6(u, v);
            if cmp_6_6(b, c) == Less {
                b = add_6_6_no_overflow(b, order);
            }
            b = sub_6_6(b, c);
        }
    }

    if u == one {
        b
    } else {
        c
    }
}
