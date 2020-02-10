use num::BigUint;

pub fn u64_slice_to_biguint(n: &[u64]) -> BigUint {
    let mut bytes_le = Vec::new();
    for n_i in n {
        for j in 0..8 {
            bytes_le.push((n_i >> j * 8) as u8);
        }
    }
    BigUint::from_bytes_le(&bytes_le)
}

pub fn biguint_to_u64_vec(biguint: BigUint, len: usize) -> Vec<u64> {
    let mut vec = Vec::new();
    for u32_pair in biguint.to_u32_digits().chunks(2) {
        let mut limb = u32_pair[0] as u64;
        if u32_pair.len() == 2 {
            limb += (u32_pair[1] as u64) << 32;
        }
        vec.push(limb)
    }

    debug_assert!(vec.len() <= len);

    while vec.len() < len {
        vec.push(0);
    }

    vec
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use num::BigUint;

    use crate::conversions::u64_slice_to_biguint;

    #[test]
    fn convert_single_u64() {
        let biguint = u64_slice_to_biguint(&[12379813738877118345]);
        assert_eq!(biguint, BigUint::from_str("12379813738877118345").unwrap());
    }
}
