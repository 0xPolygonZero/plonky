use crate::{AffinePoint, Curve, Field, HaloCurve, PallasBase, VestaBase};
use serde::{Deserialize, Serialize};

#[derive(Debug, Copy, Clone, Eq, PartialEq, Default, Serialize, Deserialize)]
pub struct Pallas;

impl Curve for Pallas {
    type BaseField = PallasBase;
    type ScalarField = VestaBase;

    const A: PallasBase = PallasBase::ZERO;
    const B: PallasBase = PallasBase::FIVE;

    const GENERATOR_AFFINE: AffinePoint<Self> = AffinePoint {
        x: PallasBase::NEG_ONE,
        y: PallasBase::TWO,
        zero: false,
    };
}

impl HaloCurve for Pallas {
    // See https://github.com/zcash/pasta/blob/fb448f35380263143160b6f190cca110858b0473/amicable.sage#L149
    // for how to find zeta and "zeta scalar"
    const ZETA: Self::BaseField = PallasBase {
        limbs: [
            0xfbdfd7aa9e65eac8, 0x0cd4d654e50025fb,
            0xd59892a33785b99a, 0x2a27fb62585e8789
        ],
    };
    const ZETA_SCALAR: Self::ScalarField = VestaBase {
        limbs: [
            0x410e7d207feeeee3, 0x6afdf14fd8fa2279,
            0xfd3d8a04eca4d4d7, 0x2de2d60777dba4ef
        ],
    };
}

#[cfg(test)]
mod tests {
    use crate::curve::{Curve, HaloCurve, ProjectivePoint};
    use crate::{Field, Pallas};

    /// A simple, somewhat inefficient implementation of multiplication which is used as a reference
    /// for correctness.
    fn mul_naive(
        lhs: <Pallas as Curve>::ScalarField,
        rhs: ProjectivePoint<Pallas>,
    ) -> ProjectivePoint<Pallas> {
        let mut g = rhs;
        let mut sum = ProjectivePoint::ZERO;
        for limb in lhs.to_canonical().iter() {
            for j in 0..64 {
                if (limb >> j & 1u64) != 0u64 {
                    sum = sum + g;
                }
                g = g.double();
            }
        }
        sum
    }

    #[test]
    fn test_endomorphism_pallas() {
        type C = Pallas;
        let g = C::convert(<C as Curve>::ScalarField::rand()) * C::GENERATOR_PROJECTIVE;
        let g = g.to_affine();
        let h = g.endomorphism();
        assert_eq!(
            h,
            mul_naive(Pallas::ZETA_SCALAR, g.to_projective()).to_affine()
        );
    }

    #[test]
    fn is_safe_curve() {
        type C = Pallas;
        assert!(
           C::is_safe_curve()
        );
    }
}
