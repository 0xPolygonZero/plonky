use std::cmp::Ordering::Less;

use crate::{add_no_overflow, cmp, div2, is_even, is_odd, sub, one_array};

pub(crate) fn nonzero_multiplicative_inverse<const N: usize>(a: [u64; N], order: [u64; N]) -> [u64; N] {
    // Based on Algorithm 16 of "Efficient Software-Implementation of Finite Fields with
    // Applications to Cryptography".

    let zero = [0; N];
    let one = one_array![u64; N];

    let mut u = a;
    let mut v = order;
    let mut b = one;
    let mut c = zero;

    while u != one && v != one {
        while is_even(u) {
            u = div2(u);
            if is_odd(b) {
                b = add_no_overflow(b, order);
            }
            b = div2(b);
        }

        while is_even(v) {
            v = div2(v);
            if is_odd(c) {
                c = add_no_overflow(c, order);
            }
            c = div2(c);
        }

        if cmp(u, v) == Less {
            v = sub(v, u);
            if cmp(c, b) == Less {
                c = add_no_overflow(c, order);
            }
            c = sub(c, b);
        } else {
            u = sub(u, v);
            if cmp(b, c) == Less {
                b = add_no_overflow(b, order);
            }
            b = sub(b, c);
        }
    }

    if u == one {
        b
    } else {
        c
    }
}
