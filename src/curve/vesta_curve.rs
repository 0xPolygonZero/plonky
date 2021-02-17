use crate::{AffinePoint, Curve, Field, HaloCurve, PallasBase, VestaBase};
use serde::{Deserialize, Serialize};

#[derive(Debug, Copy, Clone, Eq, PartialEq, Default, Serialize, Deserialize)]
pub struct Vesta;

impl Curve for Vesta {
    type BaseField = VestaBase;
    type ScalarField = PallasBase;

    const A: VestaBase = VestaBase::ZERO;
    const B: VestaBase = VestaBase::FIVE;

    const GENERATOR_AFFINE: AffinePoint<Self> = AffinePoint {
        x: VestaBase::NEG_ONE,
        y: VestaBase::TWO,
        zero: false,
    };
}

impl HaloCurve for Vesta {
    const ZETA: Self::BaseField = VestaBase {
        limbs: [
            0x410e7d207feeeee3, 0x6afdf14fd8fa2279,
            0xfd3d8a04eca4d4d7, 0x2de2d60777dba4ef
        ]
    };
    const ZETA_SCALAR: Self::ScalarField = PallasBase {
        limbs: [
            0xfbdfd7aa9e65eac8, 0x0cd4d654e50025fb,
            0xd59892a33785b99a, 0x2a27fb62585e8789
        ]
    };
}

#[cfg(test)]
mod tests {
    use crate::curve::{Curve, HaloCurve, ProjectivePoint};
    use crate::{Vesta, Field};

    /// A simple, somewhat inefficient implementation of multiplication which is used as a reference
    /// for correctness.
    fn mul_naive(
        lhs: <Vesta as Curve>::ScalarField,
        rhs: ProjectivePoint<Vesta>,
    ) -> ProjectivePoint<Vesta> {
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
    fn test_endomorphism_vesta() {
        type C = Vesta;
        let g = C::convert(<C as Curve>::ScalarField::rand()) * C::GENERATOR_PROJECTIVE;
        let g = g.to_affine();
        let h = g.endomorphism();
        assert_eq!(
            h,
            mul_naive(C::ZETA_SCALAR, g.to_projective()).to_affine()
        );
    }

    #[test]
    fn is_safe_curve() {
        type C = Vesta;
        assert!(
            C::is_safe_curve()
        );
    }
}
