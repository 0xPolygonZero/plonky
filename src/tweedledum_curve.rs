use crate::{Curve, Field, TweedledeeBase, TweedledumBase, AffinePoint, ProjectivePoint, HaloEndomorphismCurve};

#[derive(Debug, Copy, Clone)]
pub struct Tweedledum;

pub const TWEEDLEDUM_GENERATOR_AFFINE: AffinePoint<Tweedledum> = AffinePoint {
    x: TweedledumBase::NEG_ONE,
    y: TweedledumBase::TWO,
    zero: false,
};

pub const TWEEDLEDUM_GENERATOR_PROJECTIVE: ProjectivePoint<Tweedledum> = ProjectivePoint {
    x: TweedledumBase::NEG_ONE,
    y: TweedledumBase::TWO,
    z: TweedledumBase::ONE,
    zero: false,
};

impl Curve for Tweedledum {
    type BaseField = TweedledumBase;
    type ScalarField = TweedledeeBase;

    const A: TweedledumBase = TweedledumBase::ZERO;
    const B: TweedledumBase = TweedledumBase::ONE;
}

impl HaloEndomorphismCurve for Tweedledum {
    // TODO: Configure zeta
    const ZETA: Self::BaseField = TweedledumBase::ZERO;
}
