use crate::{AffinePoint, Curve, Field, HaloEndomorphismCurve, ProjectivePoint, TweedledeeBase, TweedledumBase};

#[derive(Debug, Copy, Clone)]
pub struct Tweedledum;

impl Curve for Tweedledum {
    type BaseField = TweedledumBase;
    type ScalarField = TweedledeeBase;

    const A: TweedledumBase = TweedledumBase::ZERO;
    const B: TweedledumBase = TweedledumBase::ONE;

    const GENERATOR_AFFINE: AffinePoint<Self> = AffinePoint {
        x: TweedledumBase::NEG_ONE,
        y: TweedledumBase::TWO,
        zero: false,
    };
}

impl HaloEndomorphismCurve for Tweedledum {
    // TODO: Configure zeta
    const ZETA: Self::BaseField = TweedledumBase::ZERO;
}
