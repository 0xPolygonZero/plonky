use crate::{AffinePoint, Curve, Field, HaloEndomorphismCurve, TweedledeeBase, TweedledumBase};

#[derive(Debug, Copy, Clone)]
pub struct Tweedledum;

impl Curve for Tweedledum {
    type BaseField = TweedledumBase;
    type ScalarField = TweedledeeBase;

    const A: TweedledumBase = TweedledumBase::ZERO;
    // B = 7
    const B: TweedledumBase = TweedledumBase {
        limbs: [
            18317648394857742309,
            11556519044732520811,
            18446744073709551615,
            4611686018427387903,
        ]
    };

    const GENERATOR_AFFINE: AffinePoint<Self> = AffinePoint {
        x: TweedledumBase::ONE,
        y: TweedledumBase {
            limbs: [
                12815994359195135157,
                12442237869110527732,
                9256472484777506843,
                1114242145010923164,
            ]
        },
        zero: false,
    };
}

impl HaloEndomorphismCurve for Tweedledum {
    const ZETA: Self::BaseField = TweedledumBase {
        limbs: [
            7605997034305223424,
            3132214451552427455,
            3308921103222877309,
            2709928666517121162,
        ]
    };
    const ZETA_SCALAR: Self::ScalarField = TweedledeeBase {
        limbs: [
            9282944046338294407,
            16421485501699768486,
            18374227564572127422,
            3902997619921080662,
        ]
    };
}

#[cfg(test)]
mod tests {
    use crate::curve::{AffinePoint, Curve, HaloEndomorphismCurve, ProjectivePoint};
    use crate::Tweedledum;

    /// A simple, somewhat inefficient implementation of multiplication which is used as a reference
    /// for correctness.
    fn mul_naive(
        lhs: <Tweedledum as Curve>::ScalarField,
        rhs: ProjectivePoint<Tweedledum>,
    ) -> ProjectivePoint<Tweedledum> {
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
    fn test_endomorphism_tweedledum() {
        let g = Tweedledum::GENERATOR_AFFINE;
        assert!(g.is_valid());
        let h = AffinePoint::<Tweedledum> {
            x: g.x * Tweedledum::ZETA,
            y: g.y,
            zero: false,
        };
        assert_eq!(
            h,
            mul_naive(Tweedledum::ZETA_SCALAR, g.to_projective()).to_affine()
        );
    }
}
