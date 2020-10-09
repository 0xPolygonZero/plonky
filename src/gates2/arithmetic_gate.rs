use std::marker::PhantomData;
use std::rc::Rc;

use crate::{CircuitBuilder2, ConstraintPolynomial, DeterministicGate, Field, Gate2, GateInstance, Target2};

pub const ID: &'static str = "ARITHMETIC";
pub const CONST_PRODUCT_WEIGHT: usize = 0;
pub const CONST_ADDEND_WEIGHT: usize = 1;
pub const WIRE_MULTIPLICAND_0: usize = 0;
pub const WIRE_MULTIPLICAND_1: usize = 1;
pub const WIRE_ADDEND: usize = 2;
pub const WIRE_OUTPUT: usize = 3;

/// A gate which can be configured to perform various arithmetic. In particular, it computes
///
/// ```text
/// output := const_product_weight * multiplicand_0 * multiplicand_1
///         + const_addend_weight * addend
/// ```
struct ArithmeticGate<F: Field> {
    _phantom: PhantomData<F>
}

impl<F: Field> DeterministicGate<F> for ArithmeticGate<F> {
    fn id(&self) -> String {
        "ArithmeticGate".into()
    }

    fn outputs(&self) -> Vec<(usize, ConstraintPolynomial<F>)> {
        let const_0 = ConstraintPolynomial::local_constant(CONST_PRODUCT_WEIGHT);
        let const_1 = ConstraintPolynomial::local_constant(CONST_ADDEND_WEIGHT);
        let multiplicand_0 = ConstraintPolynomial::local_wire_value(WIRE_MULTIPLICAND_0);
        let multiplicand_1 = ConstraintPolynomial::local_wire_value(WIRE_MULTIPLICAND_1);
        let addend = ConstraintPolynomial::local_wire_value(WIRE_ADDEND);

        let out = const_0 * multiplicand_0 * &multiplicand_1 + const_1 * &addend;
        vec![(WIRE_OUTPUT, out)]
    }
}

fn gate_type<F: Field>() -> Rc<Gate2<F>> {
    todo!()
}

fn add<F: Field>(builder: &mut CircuitBuilder2<F>, x: Target2<F>, y: Target2<F>) -> Target2<F> {
    let gate_type = gate_type();
    let constants = vec![F::ONE, F::ONE];
    let gate = builder.add_gate(GateInstance { gate_type, constants });
    let one = builder.one();

    builder.copy(x, Target2::wire(gate, WIRE_MULTIPLICAND_0));
    builder.copy(one, Target2::wire(gate, WIRE_MULTIPLICAND_1));
    builder.copy(y, Target2::wire(gate, WIRE_ADDEND));

    Target2::wire(gate, WIRE_OUTPUT)
}
