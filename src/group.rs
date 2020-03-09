use crate::Field;

pub trait Group {
    type BaseField: Field;
    type ScalarField: Field;
}
