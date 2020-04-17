use crate::{AffinePoint, Curve, Field, HaloEndomorphismCurve, TweedledeeBase, TweedledumBase};

#[derive(Debug, Copy, Clone)]
pub struct Tweedledee;

impl Curve for Tweedledee {
    type BaseField = TweedledeeBase;
    type ScalarField = TweedledumBase;

    const A: TweedledeeBase = TweedledeeBase::ZERO;
    const B: TweedledeeBase = TweedledeeBase::ONE;

    const GENERATOR_AFFINE: AffinePoint<Self> = AffinePoint {
        x: TweedledeeBase::NEG_ONE,
        y: TweedledeeBase::TWO,
        zero: false,
    };
}

impl HaloEndomorphismCurve for Tweedledee {
    // TODO: Configure zeta
    const ZETA: Self::BaseField = TweedledeeBase::ZERO;
}
