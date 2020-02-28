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
