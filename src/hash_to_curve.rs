use crate::{AffinePoint, Curve, Field, PRFBasedPRG, RescuePrf, PRG};

pub fn hash_u32_to_curve<C: Curve>(seed: u32) -> AffinePoint<C> {
    hash_base_field_to_curve(C::BaseField::from_canonical_u32(seed))
}

pub fn hash_base_field_to_curve<C: Curve>(mut seed: C::BaseField) -> AffinePoint<C> {
    // This is the try-and-increment method from https://eprint.iacr.org/2009/226.pdf
    let prf = RescuePrf::default();
    let mut prg = PRFBasedPRG::seeded(prf, seed);
    loop {
        let x = prg.next_field();
        let y_neg = prg.next_bool();

        // We compute x^3 + a x + b, then check if it's a square in the field. If it is (which
        // occurs with a probability of ~0.5), we have found a point on the curve.
        let square_candidate = x.cube() + C::A * x + C::B;
        if let Some(mut y) = square_candidate.square_root() {
            if y_neg {
                y = -y;
            }
            return AffinePoint::nonzero(x, y)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{hash_u32_to_curve, Tweedledum};

    #[test]
    fn test_hash_u32_to_point() {
        // Just make sure it runs with no errors.
        for i in 0..5 {
            hash_u32_to_curve::<Tweedledum>(i);
        }
    }
}
