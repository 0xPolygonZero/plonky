use crate::{AffinePoint, Curve, Field, HaloCurve, TweedledeeBase, TweedledumBase};
use serde::{Deserialize, Serialize};

#[derive(Debug, Copy, Clone, Eq, PartialEq, Default, Serialize, Deserialize)]
pub struct Tweedledee;

impl Curve for Tweedledee {
    type BaseField = TweedledeeBase;
    type ScalarField = TweedledumBase;

    const A: TweedledeeBase = TweedledeeBase::ZERO;
    const B: TweedledeeBase = TweedledeeBase::FIVE;

    const GENERATOR_AFFINE: AffinePoint<Self> = AffinePoint {
        x: TweedledeeBase::NEG_ONE,
        y: TweedledeeBase::TWO,
        zero: false,
    };
}

impl HaloCurve for Tweedledee {
    const ZETA: Self::BaseField = TweedledeeBase {
        limbs: [
            1444470991491022206,
            3301226169728360777,
            72516509137424193,
            708688398506307241,
        ],
    };
    const ZETA_SCALAR: Self::ScalarField = TweedledumBase {
        limbs: [
            13597504620482004229,
            16590497220115833568,
            15137822970486674306,
            1901757351910266741,
        ],
    };
}

#[cfg(test)]
mod tests {
    use crate::curve::{Curve, HaloCurve, ProjectivePoint};
    use crate::{Field, Tweedledee};

    /// A simple, somewhat inefficient implementation of multiplication which is used as a reference
    /// for correctness.
    fn mul_naive(
        lhs: <Tweedledee as Curve>::ScalarField,
        rhs: ProjectivePoint<Tweedledee>,
    ) -> ProjectivePoint<Tweedledee> {
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
    fn test_endomorphism_tweedledee() {
        type C = Tweedledee;
        let g = C::convert(<C as Curve>::ScalarField::rand()) * C::GENERATOR_PROJECTIVE;
        let g = g.to_affine();
        let h = g.endomorphism();
        assert_eq!(
            h,
            mul_naive(Tweedledee::ZETA_SCALAR, g.to_projective()).to_affine()
        );
    }

    #[test]
    fn is_safe_curve() {
        type C = Tweedledee;
        assert!(
           C::is_safe_curve()
        );
    }
}
