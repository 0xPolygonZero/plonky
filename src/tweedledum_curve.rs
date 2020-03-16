use crate::{Curve, Field, TweedledeeBase, TweedledumBase};

#[derive(Debug, Copy, Clone)]
pub struct Tweedledum;

impl Curve for Tweedledum {
    type BaseField = TweedledumBase;
    type ScalarField = TweedledeeBase;

    const A: TweedledumBase = TweedledumBase::ZERO;
    const B: TweedledumBase = TweedledumBase::ONE;
}
