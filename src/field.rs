use std::convert::TryInto;
use std::ops::Mul;

/// An element of the BLS12 group's base field.
#[derive(Copy, Clone)]
pub struct Bls12Base {
    /// The limbs in little-endian form.
    limbs: [u64; 6],
}

/// An element of the BLS12 group's scalar field.
#[derive(Copy, Clone)]
pub struct Bls12Scalar {
    /// The limbs in little-endian form.
    limbs: [u64; 4],
}

impl Bls12Base {
    pub const ZERO: Bls12Base = Bls12Base { limbs: [0; 6] };

    // The order of the field.
    pub const ORDER: [u64; 6] = [13402431016077863595, 2210141511517208575, 7435674573564081700, 7239337960414712511, 5412103778470702295, 1873798617647539866];

    // Precomputed R for the Barrett reduction algorithm.
    const BARRET_FACTOR: [u64; 6] = [17027978386419893992, 5649138592172459777, 3421924034565217767, 11848418460761227941, 4080332095855958760, 2837504485842123031];
    const BARRET_K: usize = 381;

    pub fn mul(a: Bls12Base, b: Bls12Base) -> Bls12Base {
        let product = mul_6_6(a.limbs, b.limbs);
        // This is the Barrett reduction algorithm.
        let product_r = mul_12_6(product, Bls12Base::BARRET_FACTOR);
        let product_r_shifted = todo!();
        let product_r_shifted_n = mul_6_6(product_r_shifted, Bls12Base::ORDER);
        let result = sub_12x64(product, product_r_shifted_n);

        // The 6 higher-order limbs should all be 0 after the subtraction. Truncate them off.
        for i in 6..12 {
            assert_eq!(result[i], 0);
        }
        let result_slice = &result[0..6];
        let limbs: [u64; 6] = result_slice.try_into().unwrap();
        Bls12Base { limbs }
    }
}


impl Bls12Scalar {
    pub const ZERO: Bls12Scalar = Bls12Scalar { limbs: [0; 4] };

    // The order of the field.
    pub const ORDER: [u64; 4] = [18446744069414584321, 6034159408538082302, 3691218898639771653, 8353516859464449352];

    // Precomputed R for the Barrett reduction algorithm.
    const BARRET_CONSTANT: [u64; 4] = [5808762262936312036, 15654811016218471260, 1021603728894469044, 10183805594867568095];
    const BARRET_K: usize = 255;

    pub fn mul(a: Bls12Scalar, b: Bls12Scalar) -> Bls12Scalar {
        let mut acc = Bls12Scalar::ZERO;
        todo!();
//        mac_4x64(acc.limbs, a.limbs, b.limbs);
        acc
    }
}

impl Mul<Bls12Base> for Bls12Base {
    type Output = Bls12Base;

    fn mul(self, rhs: Bls12Base) -> Bls12Base {
        todo!()
    }
}

fn sub_12x64(a: [u64; 12], b: [u64; 12]) -> [u64; 12] {
    todo!()
}

fn mul_6_6(a: [u64; 6], b: [u64; 6]) -> [u64; 12] {
    // Grade school multiplication.
    let mut acc128 = [0u128; 12];

    for i in 0..6 {
        for j in 0..6 {
            let a_i_b_j = a[i] as u128 * b[j] as u128;
            // Add the least significant chunk to the less significant accumulator.
            acc128[i + j] += a_i_b_j as u64 as u128;
            // Add the more significant chunk to the more significant accumulator.
            acc128[i + j + 1] = a_i_b_j >> 64;
        }
    }

    let mut acc = [0u64; 12];
    acc[0] = acc128[0] as u64;
    let mut carry = false;
    for i in 1..12 {
        let last_chunk_big = (acc[i - 1] >> 64) as u64;
        let curr_chunk_small = acc[i] as u64;
        // Note that last_chunk_big won't get anywhere near 2^64, since it's essentially a carry
        // from some additions in the previous phase, so we can add the carry bit to it without
        // fear of overflow.
        let result = curr_chunk_small.overflowing_add(last_chunk_big + carry as u64);
        acc[i] += result.0;
        carry = result.1;
    }
    acc
}

fn mul_12_6(a: [u64; 12], b: [u64; 6]) -> [u64; 18] {
    // Grade school multiplication.
    let mut acc128 = [0u128; 18];

    for i in 0..12 {
        for j in 0..6 {
            let a_i_b_j = a[i] as u128 * b[j] as u128;
            // Add the least significant chunk to the less significant accumulator.
            acc128[i + j] += a_i_b_j as u64 as u128;
            // Add the more significant chunk to the more significant accumulator.
            acc128[i + j + 1] = a_i_b_j >> 64;
        }
    }

    let mut acc = [0u64; 18];
    acc[0] = acc128[0] as u64;
    let mut carry = false;
    for i in 1..18 {
        let last_chunk_big = (acc[i - 1] >> 64) as u64;
        let curr_chunk_small = acc[i] as u64;
        // Note that last_chunk_big won't get anywhere near 2^64, since it's essentially a carry
        // from some additions in the previous phase, so we can add the carry bit to it without
        // fear of overflow.
        let result = curr_chunk_small.overflowing_add(last_chunk_big + carry as u64);
        acc[i] += result.0;
        carry = result.1;
    }
    acc
}

#[cfg(test)]
mod tests {
    #[test]
    fn todo() {
        todo!()
    }
}
