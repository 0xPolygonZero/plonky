use crate::{rescue_sponge, util::ceil_div_usize, AffinePoint, Curve, Field};
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
    let mut hasher = blake3::Hasher::new();
    // Number of bytes required to make a field element.
    let byte_length = F::BYTES;
    // Bytes version of the field element.
    let mut bytes = seed.to_canonical_u8_vec();
    (0..2).for_each(|_| bytes.push(0));
    // Add the `iter` value to get a different result at each iteration.
    bytes[byte_length] = iter;
    // One extra-byte for `y_neg`.
    let mut hash_container = vec![0; byte_length + 1];
    // Loop index.
    let mut j = 0;
    loop {
        // Add the loop index to get a different result if the hash fails.
        bytes[byte_length + 1] = j;
        if j > 0 {
            // Reset the container to zeros.
            hash_container.iter_mut().for_each(|x| *x = 0);
            // Reset the hasher.
            hasher.reset();
        }
        // Fill the container with an extended hash.
        hasher
            .update(&bytes)
            .finalize_xof()
            .fill(&mut hash_container);
        hash_container[byte_length-1] >>= 8*F::BYTES - F::BITS;
        // Try to convert the hash to a field element.
        let x = F::from_canonical_u8_vec(hash_container[..byte_length].to_vec());
        if let Ok(good) = x {
            // Use the extra-byte in `hash_container` to deduce `y_neg`.
            let y_neg = hash_container.last().unwrap() & 1 == 1;
            return (good, y_neg);
        } else {
            j += 1;
        }
    }
}

pub fn blake_hash_usize_to_curve<C: Curve>(seed: usize) -> AffinePoint<C> {
    blake_hash_base_field_to_curve(C::BaseField::from_canonical_usize(seed))
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
        TweedledumBase, Bls12377, Bls12377Base
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
    fn test_hash_blake_deterministic() {
        let n = 10000;
        let points_dee: Vec<_> = (0..n).map(|_| TweedledeeBase::rand()).collect();
        let points_dum: Vec<_> = (0..n).map(|_| TweedledumBase::rand()).collect();
        let points_377: Vec<_> = (0..n).map(|_| Bls12377Base::rand()).collect();
        let now = Instant::now();
        for i in 0..n {
            let x = blake_hash_base_field_to_curve::<Tweedledee>(points_dee[i]);
            let y = blake_hash_base_field_to_curve::<Tweedledee>(points_dee[i]);
            assert_eq!(x, y);
            let x = blake_hash_base_field_to_curve::<Tweedledum>(points_dum[i]);
            let y = blake_hash_base_field_to_curve::<Tweedledum>(points_dum[i]);
            assert_eq!(x, y);
            let x = blake_hash_base_field_to_curve::<Bls12377>(points_377[i]);
            let y = blake_hash_base_field_to_curve::<Bls12377>(points_377[i]);
            assert_eq!(x, y);
        }
        println!("Elapsed: {:.2?}", now.elapsed());
    }
}
