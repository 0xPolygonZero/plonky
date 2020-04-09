use crate::{AffinePoint, Curve, Field, HaloEndomorphismCurve, ProjectivePoint, TweedledeeBase, TweedledumBase};

#[derive(Debug, Copy, Clone)]
pub struct Tweedledee;

pub const TWEEDLEDEE_GENERATOR_AFFINE: AffinePoint<Tweedledee> = AffinePoint {
    x: TweedledeeBase::NEG_ONE,
    y: TweedledeeBase::TWO,
    zero: false,
};

pub const TWEEDLEDEE_GENERATOR_PROJECTIVE: ProjectivePoint<Tweedledee> = ProjectivePoint {
    x: TweedledeeBase::NEG_ONE,
    y: TweedledeeBase::TWO,
    z: TweedledeeBase::ONE,
    zero: false,
};

impl Curve for Tweedledee {
    type BaseField = TweedledeeBase;
    type ScalarField = TweedledumBase;

    const A: TweedledeeBase = TweedledeeBase::ZERO;
    const B: TweedledeeBase = TweedledeeBase::ONE;
}

impl HaloEndomorphismCurve for Tweedledee {
    // TODO: Configure zeta
    const ZETA: Self::BaseField = TweedledeeBase::ZERO;
}
