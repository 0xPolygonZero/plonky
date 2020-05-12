use crate::{rescue_sponge, AffinePoint, Curve, Field};

pub fn hash_u32_to_curve<C: Curve>(seed: u32, security_bits: usize) -> AffinePoint<C> {
    let seed_f = C::BaseField::from_canonical_u32(seed);
    hash_base_field_to_curve(seed_f, security_bits)
}

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
    use crate::{hash_u32_to_curve, rescue_sponge, Field, Tweedledum, TweedledumBase};

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
}
