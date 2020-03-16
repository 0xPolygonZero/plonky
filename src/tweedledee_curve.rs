use crate::{Curve, Field, TweedledeeBase, TweedledumBase};

#[derive(Debug, Copy, Clone)]
pub struct Tweedledee;

impl Curve for Tweedledee {
    type BaseField = TweedledeeBase;
    type ScalarField = TweedledumBase;

    const A: TweedledeeBase = TweedledeeBase::ZERO;
    const B: TweedledeeBase = TweedledeeBase::ONE;
}
