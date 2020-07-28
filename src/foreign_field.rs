use crate::{HaloCurve, CircuitBuilder, BigIntTarget, Field, field_to_biguint};
use num::{BigUint, One};

/// Represents an element of a field other than the native field.
#[derive(Clone)]
pub struct ForeignFieldTarget {
    pub value: BigIntTarget,
}

impl ForeignFieldTarget {
    pub fn zero() -> Self {
        ForeignFieldTarget { value: BigIntTarget::zero() }
    }
}

impl<C: HaloCurve> CircuitBuilder<C> {
    pub fn constant_foreign_field<FF: Field>(&mut self, constant: FF) -> ForeignFieldTarget {
        let value = self.constant_bigint(&field_to_biguint(constant));
        ForeignFieldTarget { value }
    }

    pub fn foreign_field_add_many<FF: Field>(
        &mut self,
        terms: &[ForeignFieldTarget],
    ) -> ForeignFieldTarget {
        let term_values = terms.iter()
            .map(|ff| ff.value.clone())
            .collect::<Vec<_>>();
        let sum = self.bigint_add_many(&term_values);
        self.reduce::<FF>(&sum)
    }

    pub fn foreign_field_add<FF: Field>(
        &mut self,
        x: &ForeignFieldTarget,
        y: &ForeignFieldTarget,
    ) -> ForeignFieldTarget {
        self.foreign_field_add_many::<FF>(&[x.clone(), y.clone()])
    }

    pub fn foreign_field_mul<FF: Field>(
        &mut self,
        lhs: &ForeignFieldTarget,
        rhs: &ForeignFieldTarget,
    ) -> ForeignFieldTarget {
        let product = self.bigint_mul(&lhs.value, &rhs.value);
        self.reduce::<FF>(&product)
    }

    /// Returns `x % |FF|` as a `ForeignFieldTarget`.
    fn reduce<FF: Field>(&mut self, x: &BigIntTarget) -> ForeignFieldTarget {
        let order = field_to_biguint(FF::NEG_ONE) + BigUint::one();
        let order_target = self.constant_bigint(&order);
        let value = self.bigint_rem(&x, &order_target);
        ForeignFieldTarget { value }
    }
}

#[cfg(test)]
mod tests {
    use crate::{CircuitBuilder, PartialWitness, Tweedledum, Field, Curve};
    use num::{BigUint, FromPrimitive};

    #[test]
    fn test_foreign_field_add() {
        type C = Tweedledum;
        type FF = <C as Curve>::ScalarField;

        let x_value = FF::rand();
        let y_value = FF::rand();
        let expected_z_value = x_value + y_value;

        let mut builder = CircuitBuilder::<C>::new(128);
        let x = builder.constant_foreign_field(x_value);
        let y = builder.constant_foreign_field(y_value);
        let z = builder.foreign_field_add::<FF>(&x, &y);
        let circuit = builder.build();

        let witness = circuit.generate_partial_witness(PartialWitness::new());
        let actual_z_value: FF = witness.get_foreign_field_target(&z);
        assert_eq!(actual_z_value, expected_z_value);
    }

    #[test]
    fn test_foreign_field_mul() {
        type C = Tweedledum;
        type FF = <C as Curve>::ScalarField;

        let x_value = FF::rand();
        let y_value = FF::rand();
        let expected_z_value = x_value * y_value;

        let mut builder = CircuitBuilder::<C>::new(128);
        let x = builder.constant_foreign_field(x_value);
        let y = builder.constant_foreign_field(y_value);
        let z = builder.foreign_field_mul::<FF>(&x, &y);
        let circuit = builder.build();

        let witness = circuit.generate_partial_witness(PartialWitness::new());
        let actual_z_value: FF = witness.get_foreign_field_target(&z);
        assert_eq!(actual_z_value, expected_z_value);
    }
}
