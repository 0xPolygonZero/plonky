use std::convert::Infallible;
use std::marker::PhantomData;

use crate::{Field, Wire};

/// A location in the witness.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum Target2<F: Field> {
    Wire(Wire),
    PublicInput { index: usize },
    VirtualAdviceTarget { index: usize },
    // Trick taken from https://github.com/rust-lang/rust/issues/32739#issuecomment-627765543.
    _Field(Infallible, PhantomData<F>),
}

impl<F: Field> Target2<F> {
    pub fn wire(gate: usize, input: usize) -> Self {
        Self::Wire(Wire { gate, input })
    }

    pub fn is_routable(&self) -> bool {
        match self {
            Target2::Wire(wire) => wire.is_routable(),
            Target2::PublicInput { .. } => true,
            Target2::VirtualAdviceTarget { .. } => false,
            Target2::_Field(_, _) => unreachable!(),
        }
    }
}
