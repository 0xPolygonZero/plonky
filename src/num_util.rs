use num::{BigInt, BigUint, One, Zero};

use num::bigint::ToBigInt;

fn egcd(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
    if a.is_zero() {
        (b, BigInt::zero(), BigInt::one())
    } else {
        let (g, y, x) = egcd(&b % &a, a.clone());
        return (g, x - (b / a) * &y, y)
    }
}

// TODO: Delete this poor implementation once it's fully replaced.
pub fn modinv(a: BigUint, m: BigUint) -> Option<BigUint> {
    let a = a.to_bigint().unwrap();
    let m = m.to_bigint().unwrap();
    let (g, x, _y) = egcd(a, m.clone());
    if g.is_one() {
        let result_bigint = (x % &m + &m) % m;
        Some((result_bigint).to_biguint().unwrap())
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use num::{BigUint, FromPrimitive, ToPrimitive};
    use crate::num_util::modinv;

    #[test]
    fn test_modinv() {
        let m = 197;
        let m_biguint = BigUint::from_u32(m).unwrap();

        for x in 0..m {
            let x_biguint = BigUint::from_u32(x).unwrap();
            let x_inv_biguint = modinv(x_biguint.clone(), m_biguint.clone());
            if x == 0 {
                assert!(x_inv_biguint.is_none());
            } else {
                let x_inv = x_inv_biguint.expect("No inverse").to_u32().unwrap();
                assert_eq!(x * x_inv % m, 1);
            }
        }
    }
}
