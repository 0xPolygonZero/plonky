use std::hash::{Hash, Hasher};
use std::iter;
use std::rc::Rc;

mod arithmetic;
mod limb_sum;

pub use arithmetic::*;
pub use limb_sum::*;

use crate::{ConstraintPolynomial, EvaluationVars, Field, PartialWitness2, SimpleGenerator, Target2, Wire, WitnessGenerator2};

pub(crate) struct GateWrapper<F: Field>(Rc<dyn Gate2<F>>);

impl<F: Field> GateWrapper<F> {
    pub fn new<G: Gate2<F>>(gate: G) -> GateWrapper<F> {
        GateWrapper(Rc::new(gate))
    }
}

impl<F: Field> PartialEq for GateWrapper<F> {
    fn eq(&self, other: &Self) -> bool {
        self.0.id() == other.0.id()
    }
}

impl<F: Field> Hash for GateWrapper<F> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.id().hash(state)
    }
}

impl<F: Field> Eq for GateWrapper<F> {}

pub trait Gate2<F: Field>: 'static {
    fn id(&self) -> String;

    /// A set of expressions which must evaluate to zero.
    fn constraints(&self) -> Vec<ConstraintPolynomial<F>>;

    fn generators(
        &self,
        gate_index: usize,
        local_constants: Vec<F>,
        next_constants: Vec<F>,
    ) -> Vec<Box<dyn WitnessGenerator2<F>>>;

    fn max_constant_index(&self) -> Option<usize>;

    fn max_wire_input_index(&self) -> Option<usize>;
}

/// A deterministic gate. Each entry in `outputs()` describes how that output is evaluated; this is
/// used to create both the constraint set and the generator set.
pub trait DeterministicGate<F: Field>: 'static {
    /// A unique identifier for this gate.
    fn id(&self) -> String;

    fn outputs(&self) -> Vec<(usize, ConstraintPolynomial<F>)>;

    /// Additional constraints to be enforced, besides the (automatically generated) ones that
    /// constraint output values.
    fn additional_constraints(&self) -> Vec<ConstraintPolynomial<F>> {
        Vec::new()
    }
}

impl<F: Field, DG: DeterministicGate<F>> Gate2<F> for DG {
    fn id(&self) -> String {
        self.id()
    }

    fn constraints(&self) -> Vec<ConstraintPolynomial<F>> {
        // For each output, we add a constraint of the form `out - expression = 0`,
        // then we append any additional constraints that the gate defines.
        self.outputs().into_iter()
            .map(|(i, out)| out - ConstraintPolynomial::local_wire_value(i))
            .chain(self.additional_constraints().into_iter())
            .collect()
    }

    fn generators(
        &self,
        gate_index: usize,
        local_constants: Vec<F>,
        next_constants: Vec<F>,
    ) -> Vec<Box<dyn WitnessGenerator2<F>>> {
        struct OutputGenerator<F: Field> {
            gate_index: usize,
            input_index: usize,
            out: ConstraintPolynomial<F>,
            local_constants: Vec<F>,
            next_constants: Vec<F>,
        }

        impl<F: Field> SimpleGenerator<F> for OutputGenerator<F> {
            fn dependencies(&self) -> Vec<Target2<F>> {
                self.out.dependencies(self.gate_index)
                    .into_iter()
                    .map(Target2::Wire)
                    .collect()
            }

            fn run_once(&mut self, witness: &PartialWitness2<F>) -> PartialWitness2<F> {
                let mut local_wire_values = Vec::new();
                let mut next_wire_values = Vec::new();

                // Get an exclusive upper bound on the largest input index in this constraint.
                let input_limit_exclusive = self.out.max_wire_input_index()
                    .map_or(0, |i| i + 1);

                for input in 0..input_limit_exclusive {
                    let local_wire = Wire { gate: self.gate_index, input };
                    let next_wire = Wire { gate: self.gate_index + 1, input };

                    // Lookup the values if they exist. If not, we can just insert a zero, knowing
                    // that it will not be used. (If it was used, it would have been included in our
                    // dependencies, and this generator would not have run yet.)
                    let local_value = witness.try_get(Target2::Wire(local_wire)).unwrap_or(F::ZERO);
                    let next_value = witness.try_get(Target2::Wire(next_wire)).unwrap_or(F::ZERO);

                    local_wire_values.push(local_value);
                    next_wire_values.push(next_value);
                }

                let vars = EvaluationVars {
                    local_constants: &self.local_constants,
                    next_constants: &self.next_constants,
                    local_wire_values: &local_wire_values,
                    next_wire_values: &next_wire_values,
                };

                let result_wire = Wire { gate: self.gate_index, input: self.input_index };
                let result_value = self.out.evaluate(vars);
                PartialWitness2::singleton(Target2::Wire(result_wire), result_value)
            }
        }

        self.outputs()
            .into_iter()
            .map(|(input_index, out)| {
                let og = OutputGenerator {
                    gate_index,
                    input_index,
                    out,
                    local_constants: local_constants.clone(),
                    next_constants: next_constants.clone(),
                };
                // We need the type system to treat this as a boxed `WitnessGenerator2<F>`, rather
                // than a boxed `OutputGenerator<F>`.
                let b: Box::<dyn WitnessGenerator2<F>> = Box::new(og);
                b
            })
            .collect()
    }

    fn max_constant_index(&self) -> Option<usize> {
        self.outputs().into_iter()
            .map(|(_i, out)| out.max_constant_index())
            .filter_map(|out_max| out_max)
            .max()
    }

    fn max_wire_input_index(&self) -> Option<usize> {
        self.outputs().into_iter()
            // For each output, we consider both the output wire and the wires it depends on.
            .flat_map(|(i, out)| out.max_wire_input_index()
                .into_iter()
                .chain(iter::once(i)))
            .max()
    }
}

/// A gate along with any constants used to configure it.
pub struct GateInstance<F: Field> {
    gate_type: GateWrapper<F>,
    constants: Vec<F>,
}
