use crate::{CircuitBuilder2, ConstraintPolynomial, DeterministicGate, DeterministicGateAdapter, Field, GateInstance, GateRef, Target2, CircuitConfig};

/// A gate which can be configured to perform various arithmetic. In particular, it computes
///
/// ```text
/// output := product_weight * multiplicand_0 * multiplicand_1
///         + addend_weight * addend
/// ```
///
/// where `product_weight` and `addend_weight` are constants, and the other variables are wires.
#[derive(Eq, PartialEq, Hash)]
struct ArithmeticGate2;

impl ArithmeticGate2 {
    pub fn get_ref<F: Field>() -> GateRef<F> {
        GateRef::new(DeterministicGateAdapter::new(ArithmeticGate2))
    }
}

impl ArithmeticGate2 {
    pub const CONST_PRODUCT_WEIGHT: usize = 0;
    pub const CONST_ADDEND_WEIGHT: usize = 1;
    pub const WIRE_MULTIPLICAND_0: usize = 0;
    pub const WIRE_MULTIPLICAND_1: usize = 1;
    pub const WIRE_ADDEND: usize = 2;
    pub const WIRE_OUTPUT: usize = 3;

    /// Computes `x y + z`.
    pub fn mul_add<F: Field>(
        builder: &mut CircuitBuilder2<F>,
        x: Target2<F>,
        y: Target2<F>,
        z: Target2<F>,
    ) -> Target2<F> {
        let gate_type = ArithmeticGate2::get_ref();
        let constants = vec![F::ONE, F::ONE];
        let gate = builder.add_gate(GateInstance { gate_type, constants });

        builder.route(x, Target2::wire(gate, Self::WIRE_MULTIPLICAND_0));
        builder.route(y, Target2::wire(gate, Self::WIRE_MULTIPLICAND_1));
        builder.route(z, Target2::wire(gate, Self::WIRE_ADDEND));

        Target2::wire(gate, Self::WIRE_OUTPUT)
    }

    /// Computes `x y`.
    pub fn mul<F: Field>(
        builder: &mut CircuitBuilder2<F>,
        x: Target2<F>,
        y: Target2<F>,
    ) -> Target2<F> {
        let zero = builder.zero();
        Self::mul_add(builder, x, y, zero)
    }

    /// Computes `x + y`.
    pub fn add<F: Field>(
        builder: &mut CircuitBuilder2<F>,
        x: Target2<F>,
        y: Target2<F>,
    ) -> Target2<F> {
        let one = builder.one();
        Self::mul_add(builder, x, one, y)
    }
}

impl<F: Field> DeterministicGate<F> for ArithmeticGate2 {
    fn id(&self) -> String {
        "ArithmeticGate".into()
    }

    fn outputs(&self, _config: CircuitConfig) -> Vec<(usize, ConstraintPolynomial<F>)> {
        let const_0 = ConstraintPolynomial::local_constant(Self::CONST_PRODUCT_WEIGHT);
        let const_1 = ConstraintPolynomial::local_constant(Self::CONST_ADDEND_WEIGHT);
        let multiplicand_0 = ConstraintPolynomial::local_wire_value(Self::WIRE_MULTIPLICAND_0);
        let multiplicand_1 = ConstraintPolynomial::local_wire_value(Self::WIRE_MULTIPLICAND_1);
        let addend = ConstraintPolynomial::local_wire_value(Self::WIRE_ADDEND);

        let out = const_0 * multiplicand_0 * &multiplicand_1 + const_1 * &addend;
        vec![(Self::WIRE_OUTPUT, out)]
    }
}

#[cfg(test)]
mod tests {
    use crate::{CircuitBuilder2, CircuitConfig, TweedledumBase};
    use crate::gates2::arithmetic::ArithmeticGate2;

    fn add() {
        let config = CircuitConfig { num_wires: 3, num_routed_wires: 3 };
        let mut builder = CircuitBuilder2::<TweedledumBase>::new(config);
        let one = builder.one();
        let two = builder.two();
        let sum = ArithmeticGate2::add(&mut builder, one, one);
        todo!()
    }
}
