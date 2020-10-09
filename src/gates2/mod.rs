use std::rc::Rc;

use crate::{Field, ConstraintPolynomial, WitnessGenerator2, Target2, PartialWitness2};

pub mod arithmetic_gate;

pub trait Gate2<F: Field> {
    fn id(&self) -> String;

    fn constraints(&self) -> Vec<ConstraintPolynomial<F>>;

    fn generators(&self, index: usize) -> Vec<Box<dyn WitnessGenerator2<F>>>;
}

/// A deterministic gate. Each entry in `outputs()` describes how that output is evaluated; this is
/// used to create both the constraint set and the generator set.
pub trait DeterministicGate<F: Field> {
    fn id(&self) -> String;

    fn outputs(&self) -> Vec<(usize, ConstraintPolynomial<F>)>;
}

impl<F: Field> Gate2<F> for dyn DeterministicGate<F> {
    fn id(&self) -> String {
        self.id()
    }

    fn constraints(&self) -> Vec<ConstraintPolynomial<F>> {
        self.outputs().into_iter()
            .map(|(i, out)| out - ConstraintPolynomial::local_wire_value(i))
            .collect()
    }

    fn generators(&self, index: usize) -> Vec<Box<dyn WitnessGenerator2<F>>> {
        struct OutputGenerator<F: Field> {
            i: usize,
            out: ConstraintPolynomial<F>,
        }

        impl<F: Field> WitnessGenerator2<F> for OutputGenerator<F> {
            fn watch_list(&self) -> Vec<Target2<F>> {
                todo!()
            }

            fn run(&mut self, witness: &PartialWitness2<F>) -> (PartialWitness2<F>, bool) {
                // ConstraintPolynomial::evaluate_all(&self.outputs)
                todo!()
            }
        }

        self.outputs()
            .into_iter()
            .map(|(i, out)| {
                let b: Box::<dyn WitnessGenerator2<F> + 'static> =
                    Box::new(OutputGenerator { i, out });
                b
            })
            .collect()
    }
}

pub struct GateInstance<F: Field> {
    gate_type: Rc<dyn Gate2<F>>,
    constants: Vec<F>,
}
