use std::cmp::Ordering::Less;
use unroll::unroll_for_loops;

use crate::{add_no_overflow, sub, cmp, mul2, nonzero_multiplicative_inverse};

// TODO: These low-level functions should be collected together in their own
// module.
#[inline]
fn mul_add_cy_in(a: u64, b: u64, c: u64, cy_in: u64) -> (u64, u64) {
    let t = (a as u128) * (b as u128) + (c as u128) + (cy_in as u128);
    ((t >> 64) as u64, t as u64)
}

#[inline]
fn mul_add(a: u64, b: u64, c: u64) -> (u64, u64) {
    let t = (a as u128) * (b as u128) + (c as u128);
    ((t >> 64) as u64, t as u64)
}


pub trait MontyRepr {
    /// The order of the field
    const ORDER: [u64; 4];
    /// Twice the order of the field
    const ORDER_X2: [u64; 4];
    /// R in the context of the Montgomery reduction, i.e. 2^256 % |F|.
    const R: [u64; 4];
    /// R^2 in the context of the Montgomery reduction, i.e. 2^(256*2) % |F|.
    const R2: [u64; 4];
    /// R^3 in the context of the Montgomery reduction, i.e. 2^(256*3) % |F|.
    const R3: [u64; 4];
    /// In the context of Montgomery multiplication, µ = -|F|^-1 mod 2^64.
    const MU: u64;

    const ZERO: [u64; 4] = [0u64; 4];
    const ONE: [u64; 4] = Self::R;

    fn monty_add(lhs: [u64; 4], rhs: [u64; 4]) -> [u64; 4] {
        // First we do a widening addition, then we reduce if necessary.
        let sum = add_no_overflow(lhs, rhs);
        if cmp(sum, Self::ORDER) == Less {
            sum
        } else {
            sub(sum, Self::ORDER)
        }
    }

    fn monty_sub(lhs: [u64; 4], rhs: [u64; 4]) -> [u64; 4] {
        if cmp(lhs, rhs) == Less {
            // Underflow occurs, so we compute the difference as `self + (-rhs)`.
            add_no_overflow(lhs, Self::monty_neg(rhs))
        } else {
            // No underflow, so it's faster to subtract directly.
            sub(lhs, rhs)
        }
    }

    fn monty_neg(limbs: [u64; 4]) -> [u64; 4] {
        // FIXME: Should have a MontyRepr::ZERO or 'is_zero()'
        if limbs == Self::ZERO {
            Self::ZERO
        } else {
            sub(Self::ORDER, limbs)
        }
    }

    #[unroll_for_loops]
    fn monty_multiply(a: [u64; 4], b: [u64; 4]) -> [u64; 4] {
        // Interleaved Montgomery multiplication, as described in Algorithm 2 of
        // https://eprint.iacr.org/2017/1057.pdf

        // Note that in the loop below, to avoid explicitly shifting c, we will treat i as the least
        // significant digit and wrap around.
        let mut c = [0u64; 5];

        for i in 0..4 {
            // Add a[i] b to c.
            let mut carry = 0;
            for j in 0..4 {
                let result = c[(i + j) % 5] as u128 + a[i] as u128 * b[j] as u128 + carry as u128;
                c[(i + j) % 5] = result as u64;
                carry = (result >> 64) as u64;
            }
            c[(i + 4) % 5] += carry;

            // q = u c mod r = u c[0] mod r.
            let q = Self::MU.wrapping_mul(c[i]);

            // C += N q
            carry = 0;
            for j in 0..4 {
                let result =
                    c[(i + j) % 5] as u128 + q as u128 * Self::ORDER[j] as u128 + carry as u128;
                c[(i + j) % 5] = result as u64;
                carry = (result >> 64) as u64;
            }
            c[(i + 4) % 5] += carry;

            debug_assert_eq!(c[i], 0);
        }

        let mut result = [c[4], c[0], c[1], c[2]];
        // Final conditional subtraction.
        if cmp(result, Self::ORDER) != Less {
            result = sub(result, Self::ORDER);
        }
        result
    }

    #[unroll_for_loops]
    fn monty_square(a: [u64; 4]) -> [u64; 4] {
        let mut c = [0u64; 4];
        let mut hi = 0u64;

        for i in 0 .. 4 {
            // u holds the off-diagonal part of the square calculation. Only the
            // first N - i elements are used.
            let mut u = [0u64; 4];
            let mut hi_in = 0u64;
            for j in i+1 .. 4 {
                let (hi_out, lo) = mul_add(a[j], a[i], hi_in);
                u[j - (i+1)] = lo;
                hi_in = hi_out;
            }
            u[4 - (i+1)] = hi_in;
            u = mul2(u);
            let (mut hi_in, lo) = mul_add(a[i], a[i], c[i]);
            c[i] = lo;
            for j in i+1 .. 4 {
                // c[j] = c[j] + u[j] + hi_in
                let (t, cy1) = c[j].overflowing_add(hi_in);
                let (t, cy2) = u[j - (i+1)].overflowing_add(t);
                c[j] = t;
                hi_in = (cy1 as u64) + (cy2 as u64);
            }

            let (t, cy1) = hi.overflowing_add(hi_in);
            let (t, cy2) = u[4 - (i+1)].overflowing_add(t);
            hi = t;
            debug_assert_eq!(cy1 | cy2, false);

            let m = c[0].wrapping_mul(Self::MU);
            let (mut hi_in, lo) = mul_add(Self::ORDER[0], m, c[0]);
            debug_assert_eq!(lo, 0u64);
            for j in 1 .. 4 {
                let (hi_out, lo) = mul_add_cy_in(Self::ORDER[j], m, c[j], hi_in);
                c[j - 1] = lo;
                hi_in = hi_out;
            }
            let (t, cy) = hi.overflowing_add(hi_in);
            c[4 - 1] = t;
            hi = cy as u64;
        }
        debug_assert_eq!(hi, 0u64);
        // Final conditional subtraction.
        if cmp(c, Self::ORDER) != Less {
            c = sub(c, Self::ORDER);
        }
        c
    }

    fn monty_inverse(limbs: [u64; 4]) -> [u64; 4] {
        // Let x R = self. We compute M((x R)^-1, R^3) = x^-1 R^-1 R^3 R^-1 =
        // x^-1 R.
        // TODO: There are faster ways to do this (for example, see McIvor,
        // McLoone and McCanny (2004)).
        let self_r_inv = nonzero_multiplicative_inverse(limbs, Self::ORDER);
        Self::monty_multiply(self_r_inv, Self::R3)
    }

    fn from_monty(c: [u64; 4]) -> [u64; 4] {
        // We compute M(c, R^2) = c * R^2 * R^-1 = c * R.
        Self::monty_multiply(c, Self::R2)
    }

    fn to_monty(limbs: [u64; 4]) -> [u64; 4] {
        // Let x * R = limbs. We compute M(x * R, 1) = x * R * R^-1 = x.
        Self::monty_multiply(limbs, [1, 0, 0, 0])
    }
}
