use crate::{Curve, Field, TweedledeeBase, TweedledumBase};

#[derive(Debug, Copy, Clone)]
pub struct Bls12377;

impl Curve for Bls12377 {
    type BaseField = TweedledumBase;
    type ScalarField = TweedledeeBase;

    const A: TweedledumBase = TweedledumBase::ZERO;
    const B: TweedledumBase = TweedledumBase::ONE;
}
