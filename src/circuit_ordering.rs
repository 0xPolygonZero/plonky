use crate::gates::gate_collection::GatePrefixes;
use crate::{CircuitBuilder, Field, HaloCurve, PartialWitness, Target, WitnessGenerator};
use std::cmp::Ordering;

#[derive(Copy, Clone)]
pub struct OrderingTarget<F: Field> {
    pub lt: Target<F>,
    pub eq: Target<F>,
    pub gt: Target<F>,
}

impl<C: HaloCurve> CircuitBuilder<C> {
    pub fn constant_ordering(&mut self, ordering: Ordering) -> OrderingTarget<C::ScalarField> {
        match ordering {
            Ordering::Less => self.ordering_lt(),
            Ordering::Equal => self.ordering_eq(),
            Ordering::Greater => self.ordering_gt(),
        }
    }

    pub fn ordering_lt(&mut self) -> OrderingTarget<C::ScalarField> {
        let _false = self.zero_wire();
        let _true = self.one_wire();
        OrderingTarget {
            lt: _true,
            eq: _false,
            gt: _false,
        }
    }

    pub fn ordering_eq(&mut self) -> OrderingTarget<C::ScalarField> {
        let _false = self.zero_wire();
        let _true = self.one_wire();
        OrderingTarget {
            lt: _false,
            eq: _true,
            gt: _false,
        }
    }

    pub fn ordering_gt(&mut self) -> OrderingTarget<C::ScalarField> {
        let _false = self.zero_wire();
        let _true = self.one_wire();
        OrderingTarget {
            lt: _false,
            eq: _false,
            gt: _true,
        }
    }

    pub fn add_virtual_ordering_target(
        &mut self,
        validate: bool,
    ) -> OrderingTarget<C::ScalarField> {
        let lt = self.add_virtual_target();
        let eq = self.add_virtual_target();
        let gt = self.add_virtual_target();

        let ordering = OrderingTarget { lt, eq, gt };
        if validate {
            self.ordering_assert_valid(ordering);
        }
        ordering
    }

    /// Adds a generator to generate `ordering` by comparing `lhs` and `rhs`.
    pub(crate) fn add_ordering_generator(
        &mut self,
        ordering: OrderingTarget<C::ScalarField>,
        lhs: Target<C::ScalarField>,
        rhs: Target<C::ScalarField>,
    ) {
        struct OrderingGenerator<F: Field> {
            ordering: OrderingTarget<F>,
            lhs: Target<F>,
            rhs: Target<F>,
        }

        impl<F: Field> WitnessGenerator<F> for OrderingGenerator<F> {
            fn dependencies(&self) -> Vec<Target<F>> {
                vec![self.lhs, self.rhs]
            }

            fn generate(
                &self,
                _prefixes: &GatePrefixes,
                _constants: &[Vec<F>],
                witness: &PartialWitness<F>,
            ) -> PartialWitness<F> {
                let lhs = witness.get_target(self.lhs);
                let rhs = witness.get_target(self.rhs);

                let mut result = PartialWitness::new();
                result.set_ordering_target(self.ordering, lhs.cmp(&rhs));
                result
            }
        }

        self.add_generator(OrderingGenerator { ordering, lhs, rhs });
    }

    pub fn ordering_assert_valid(&mut self, ordering: OrderingTarget<C::ScalarField>) {
        let OrderingTarget { lt, eq, gt } = ordering;

        self.assert_binary(lt);
        self.assert_binary(eq);
        self.assert_binary(gt);

        // Two of the three targets must be zero, so each product must be zero.
        let lt_eq = self.mul(lt, eq);
        let lt_gt = self.mul(lt, gt);
        let eq_gt = self.mul(eq, gt);
        self.assert_zero(lt_eq);
        self.assert_zero(lt_gt);
        self.assert_zero(eq_gt);

        // The remaining target must be one, so the sum must be one.
        let sum = self.add_many(&[lt, eq, gt]);
        self.assert_one(sum);
    }
}
