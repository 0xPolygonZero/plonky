use crate::{rescue_sponge, AffinePoint, Curve, Field};
use blake3;

pub fn hash_u32_to_curve<C: Curve>(seed: u32, security_bits: usize) -> AffinePoint<C> {
    let seed_f = C::BaseField::from_canonical_u32(seed);
    hash_base_field_to_curve(seed_f, security_bits)
}

pub fn hash_usize_to_curve<C: Curve>(seed: usize, security_bits: usize) -> AffinePoint<C> {
    let seed_f = C::BaseField::from_canonical_usize(seed);
    hash_base_field_to_curve(seed_f, security_bits)
}

pub fn blake_field<F: Field>(iter: u8, seed: F) -> (F, bool) {
    // Number of bytes required to make a field element.
    let byte_length = (F::BITS / 8) + if F::BITS % 8 == 0 { 0 } else { 1 };
    // Number of Blake hashes necessary to get `byte_length` bytes.
    let n = (byte_length / blake3::OUT_LEN)
        + if byte_length % blake3::OUT_LEN == 0 {
            0
        } else {
            1
        };
    // Bytes version of the field element.
    let mut bytes = seed.to_canonical_u8_vec();
    (0..3).for_each(|_| bytes.push(0));
    // Add the `iter` value to get a different result at each iteration.
    bytes[byte_length] = iter;
    let mut hash_container = vec![0; n * blake3::OUT_LEN];
    // Loop index.
    let mut j = 0;
    loop {
        // Add the loop index to get a different result if the hash fails.
        bytes[byte_length+1] = j;
        if j > 0 {
        // Reset the container to zeros.
            hash_container.iter_mut().for_each(|x| *x = 0);
        }
        // Compute `n` different Blake hashes and copy them to the container.
        for i in 0..n {
            bytes[byte_length+2] = i as u8;
            &hash_container[(i * blake3::OUT_LEN)..((i + 1) * blake3::OUT_LEN)]
                .copy_from_slice(blake3::hash(&bytes).as_bytes());
        }
        // Try to convert the hash to a field element.
        let x = F::from_canonical_u8_vec(hash_container[..byte_length].to_vec());
        if let Ok(good) = x {
            // Compute a final hash to get the `y_neg` boolean value
            // bytes.push(n as u8);
            bytes[byte_length+2] = n as u8;
            let y_neg = blake3::hash(&bytes).as_bytes()[0] & 1 == 1;
            return (good, y_neg);
        } else {
            j += 1;
        }
    }
}

pub fn blake_hash_base_field_to_curve<C: Curve>(seed: C::BaseField) -> AffinePoint<C> {
    // Based on the MapToGroup method of BLS.
    let mut i = 0;
    loop {
        // Let (x, y_neg) = H(seed, i).
        let (x, y_neg) = blake_field(i, seed);

        // We compute x^3 + a x + b, then check if it's a square in the field. If it is (which
        // occurs with a probability of ~0.5), we have found a point on the curve.
        let square_candidate = x.cube() + C::A * x + C::B;
        if let Some(mut y) = square_candidate.square_root() {
            if y_neg {
                y = -y;
            }
            return AffinePoint::nonzero(x, y);
        }

        i += 1;
    }
}

// TODO: This is rather slow! Should use ChaCha20 or something instead of Rescue.
pub fn hash_base_field_to_curve<C: Curve>(
    mut seed: C::BaseField,
    security_bits: usize,
) -> AffinePoint<C> {
    // Based on the MapToGroup method of BLS.
    let mut i = 0;
    loop {
        // Let (x, y_neg) = H(seed, i).
        let inputs = vec![seed, C::BaseField::from_canonical_u32(i)];
        let outputs = rescue_sponge(inputs, 2, security_bits);
        let x = outputs[0];
        let y_neg = outputs[1].to_canonical_bool_vec()[0];

        // We compute x^3 + a x + b, then check if it's a square in the field. If it is (which
        // occurs with a probability of ~0.5), we have found a point on the curve.
        let square_candidate = x.cube() + C::A * x + C::B;
        if let Some(mut y) = square_candidate.square_root() {
            if y_neg {
                y = -y;
            }
            return AffinePoint::nonzero(x, y);
        }

        i += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        hash_u32_to_curve, rescue_sponge, Field, Tweedledee, TweedledeeBase, Tweedledum,
        TweedledumBase,
    };
    use std::time::Instant;

    #[test]
    fn test_rescue_deterministic() {
        let inputs = [TweedledumBase::ZERO, TweedledumBase::ZERO];
        let outputs1 = rescue_sponge(inputs.to_vec(), 2, 128);
        let outputs2 = rescue_sponge(inputs.to_vec(), 2, 128);
        assert_eq!(outputs1, outputs2);
    }

    #[test]
    fn test_hash_u32_to_point() {
        // Just make sure it runs with no errors and is deterministic.
        for i in 0..5 {
            let x = hash_u32_to_curve::<Tweedledum>(i, 128);
            assert_eq!(x, hash_u32_to_curve::<Tweedledum>(i, 128));
        }
    }

    #[test]
    fn test_hash_rescue_10000() {
        let N = 10000;
        let points: Vec<_> = (0..N).map(|_| TweedledeeBase::rand()).collect();
        let now = Instant::now();
        for &p in points.iter() {
            hash_base_field_to_curve::<Tweedledee>(p, 128);
        }
        println!("Elapsed: {:.2?}", now.elapsed());
    }
    #[test]
    fn test_hash_blake_10000() {
        let N = 1000000;
        let points: Vec<_> = (0..N).map(|_| TweedledeeBase::rand()).collect();
        let now = Instant::now();
        for &p in points.iter() {
            blake_hash_base_field_to_curve::<Tweedledee>(p);
        }
        println!("Elapsed: {:.2?}", now.elapsed());
    }
    #[test]
    fn test_hash_blake_deterministic() {
        let N = 10000;
        let points: Vec<_> = (0..N).map(|_| TweedledeeBase::rand()).collect();
        let now = Instant::now();
        for &p in points.iter() {
            let x = blake_hash_base_field_to_curve::<Tweedledee>(p);
            let y = blake_hash_base_field_to_curve::<Tweedledee>(p);
            assert_eq!(x, y);
        }
        println!("Elapsed: {:.2?}", now.elapsed());
    }
}
