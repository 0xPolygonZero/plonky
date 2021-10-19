use crate::util::ceil_div_usize;
use crate::{biguint_to_field, field_to_biguint, util::pad_to_multiple_usize, Base4SumGate, BoundedTarget, CircuitBuilder, Field, HaloCurve, OrderingTarget, PartialWitness, Target, WitnessGenerator};
use num::{BigUint, Integer, One, Zero};

/// We use 86-bit limbs so that
/// - Any ~256 bit field element can be encoded as three limbs
/// - Each limb product will take at most 172 bits
/// - If our native field size is at least ~256 bits, we can accumulate many limb products witout
///   overflowing the native field
pub(crate) const LIMB_DIBITS: usize = 43;
pub(crate) const LIMB_BITS: usize = LIMB_DIBITS * 2;

/// Targets representing an unsigned big integer with limbs defined over a field `F`.
#[derive(Clone)]
pub struct BigIntTarget<F: Field> {
    pub limbs: Vec<Target<F>>,
    /// An inclusive upper bound on this number.
    pub max: BigUint,
}

impl<F: Field> BigIntTarget<F> {
    pub fn new_bounded(limbs: Vec<Target<F>>, max: BigUint) -> Self {
        Self { limbs, max }
    }

    pub fn new_unbounded(limbs: Vec<Target<F>>) -> Self {
        let bits = LIMB_BITS * limbs.len();
        let one = BigUint::one();
        let max = (&one << bits) - one;
        Self { limbs, max }
    }

    pub fn zero() -> Self {
        Self {
            limbs: Vec::new(),
            max: BigUint::zero(),
        }
    }

    pub fn num_limbs(&self) -> usize {
        self.limbs.len()
    }

    pub fn get_limb(&self, index: usize) -> Target<F> {
        self.limbs[index]
    }

    pub fn get_bounded_limb(&self, index: usize) -> BoundedTarget<F> {
        // We shift self.max to get the max value of this limb AND any more significant limbs.
        let max_high_limbs = &self.max >> (LIMB_BITS * index);
        let max_any_limb = (BigUint::one() << LIMB_BITS) - BigUint::one();
        let max_this_limb = max_high_limbs.min(max_any_limb);
        BoundedTarget {
            target: self.get_limb(index),
            max: max_this_limb,
        }
    }

    fn get_bounded_limb_or_default(
        &self,
        index: usize,
        default: BoundedTarget<F>,
    ) -> BoundedTarget<F> {
        if index < self.num_limbs() {
            self.get_bounded_limb(index)
        } else {
            default
        }
    }

    /// Return `(first, rest)`, where `first` is this bigint's least significant limb, and `rest` is
    /// a bigint consisting of the remaining limbs.
    fn split_smallest_limb(&self) -> (Target<F>, Self) {
        let first = self.get_limb(0);
        let rest_limbs = self.limbs[1..].to_vec();
        let rest_max = &self.max >> LIMB_BITS;
        let rest = BigIntTarget {
            limbs: rest_limbs,
            max: rest_max,
        };
        (first, rest)
    }
}

impl<F: Field> From<BoundedTarget<F>> for BigIntTarget<F> {
    fn from(bounded_target: BoundedTarget<F>) -> Self {
        Self {
            limbs: vec![bounded_target.target],
            max: bounded_target.max,
        }
    }
}

pub(crate) fn biguint_to_limbs<F: Field>(biguint: &BigUint) -> Vec<F> {
    let num_limbs = ceil_div_usize(biguint.bits() as usize, LIMB_BITS);
    let base = BigUint::one() << LIMB_BITS;
    (0..num_limbs)
        .map(|i| biguint_to_field((biguint >> (i * LIMB_BITS)) % &base))
        .collect()
}

impl<C: HaloCurve> CircuitBuilder<C> {
    pub fn add_virtual_bigint_target(
        &mut self,
        max: &BigUint,
        validate: bool,
    ) -> BigIntTarget<C::ScalarField> {
        let num_limbs = ceil_div_usize(max.bits() as usize, LIMB_BITS);
        let limbs = self.add_virtual_targets(num_limbs);

        if validate {
            // Check that we have a valid bigint encoding, with each limb being in the proper range.
            for &limb in &limbs {
                self.assert_dibit_length(limb, LIMB_DIBITS);
            }
        }

        BigIntTarget {
            limbs,
            max: max.clone(),
        }
    }

    pub fn constant_bigint(&mut self, value: &BigUint) -> BigIntTarget<C::ScalarField> {
        let limbs = biguint_to_limbs(value)
            .into_iter()
            .map(|limb| self.constant_wire(limb))
            .collect();
        BigIntTarget {
            limbs,
            max: value.clone(),
        }
    }

    pub fn bigint_cmp(
        &mut self,
        x: &BigIntTarget<C::ScalarField>,
        y: &BigIntTarget<C::ScalarField>,
    ) -> OrderingTarget<C::ScalarField> {
        // We test each pair of limbs for equality. Then we find the most significant pair that
        // does not match, and compare those limbs. A roughly similar approach was described in
        // https://github.com/mir-protocol/r1cs-workshop/blob/master/workshop.pdf

        // Zero-pad the inputs if needed.
        let num_limbs = x.limbs.len().max(y.limbs.len());
        let x = self.bigint_pad_limbs(x, num_limbs);
        let y = self.bigint_pad_limbs(y, num_limbs);

        // These track the most significant pair of limbs that differed, if any.
        let mut x_diff = self.zero_wire();
        let mut y_diff = self.zero_wire();

        for i in 0..num_limbs {
            let x_i = x.limbs[i];
            let y_i = y.limbs[i];
            let equal = self.is_equal(x_i, y_i);
            x_diff = self.select(equal, x_diff, x_i);
            y_diff = self.select(equal, y_diff, y_i);
        }

        self.limb_cmp(x_diff, y_diff)
    }

    fn limb_cmp(
        &mut self,
        x: Target<C::ScalarField>,
        y: Target<C::ScalarField>,
    ) -> OrderingTarget<C::ScalarField> {
        let ordering = self.add_virtual_ordering_target(true);
        let OrderingTarget { gt, eq, lt } = ordering;
        self.add_ordering_generator(ordering, x, y);

        // Check that eq == (x == y).
        let is_equal = self.is_equal(x, y);
        self.copy(eq, is_equal);

        // Now we want to check that if lt == 1, x <= y, and if gt == 1, x >= y.
        // (Strict equality is not required here since we already checked the eq case.)
        // We will do this by computing
        //     r = lt * (y - x) + gt * (x - y)
        //       = lt * (y - x) - gt * (y - x)
        // and range checking r to ensure that no underflow occurs.
        let delta = self.sub(y, x);
        let gt_delta = self.mul(gt, delta);
        let r = self.mul_sub(lt, delta, gt_delta);

        // We have some flexibility in what upper bound to use for the range check. The max
        // number of dibits must be at least LIMB_DIBITS, since r can legitimately be that large
        // without underflow. We will pad to a multiple of Base4Gate::NUM_LIMBS, since
        // assert_dibit_length is more efficient in that case. We will still detect any underflow,
        // since that would result in r consuming at least roughly F_dibits - LIMB_DIBITS, which
        // for any reasonably large field, will exceed pad(LIMB_DIBITS).
        let max_dibits = pad_to_multiple_usize(LIMB_DIBITS, Base4SumGate::<C>::NUM_LIMBS);
        self.assert_dibit_length(r, max_dibits);

        ordering
    }

    pub fn bigint_add(
        &mut self,
        x: &BigIntTarget<C::ScalarField>,
        y: &BigIntTarget<C::ScalarField>,
    ) -> BigIntTarget<C::ScalarField> {
        self.bigint_add_many(&[x.clone(), y.clone()])
    }

    pub fn bigint_add_many(
        &mut self,
        terms: &[BigIntTarget<C::ScalarField>],
    ) -> BigIntTarget<C::ScalarField> {
        let num_limbs = terms
            .iter()
            .map(BigIntTarget::num_limbs)
            .max()
            .expect("No operands");

        let mut carry: BoundedTarget<C::ScalarField> = self.zero_bounded_target();
        let mut result_limbs = Vec::new();
        for i in 0..num_limbs {
            let mut bounded_limbs = vec![carry];
            for term in terms {
                if term.num_limbs() > i {
                    bounded_limbs.push(term.get_bounded_limb(i));
                }
            }

            let sum_of_limbs = self.sum_limbs(&bounded_limbs);

            // Unless the number of terms is astronomical, we should have two limbs at most.
            assert!(sum_of_limbs.num_limbs() <= 2);

            // The first limb (or zero if there are no limbs) becomes a limb of the result.
            result_limbs.push(
                sum_of_limbs
                    .limbs
                    .get(0)
                    .cloned()
                    .unwrap_or_else(|| self.zero_wire()),
            );

            // The second limb (or zero if there isn't one) becomes our carry.
            carry = sum_of_limbs.get_bounded_limb_or_default(1, self.zero_bounded_target());
        }

        if !carry.max.is_zero() {
            result_limbs.push(carry.target);
        }

        let max = terms.iter().map(|t| &t.max).sum();
        BigIntTarget {
            limbs: result_limbs,
            max,
        }
    }

    fn sum_limbs(
        &mut self,
        limbs: &[BoundedTarget<C::ScalarField>],
    ) -> BigIntTarget<C::ScalarField> {
        let nonzero_limbs: Vec<BoundedTarget<C::ScalarField>> =
            limbs.iter().cloned().filter(|l| !l.max.is_zero()).collect();

        if nonzero_limbs.is_empty() {
            return BigIntTarget::zero();
        }

        if nonzero_limbs.len() == 1 {
            return nonzero_limbs[0].clone().into();
        }

        // Convert to a Vec of unbounded limbs.
        let nonzero_limbs: Vec<Target<C::ScalarField>> = nonzero_limbs
            .iter()
            .map(|bounded_limb| bounded_limb.target)
            .collect();

        // Compute an overall bound based on each limb's bound.
        let mut max = BigUint::zero();
        for limb in limbs {
            max += &limb.max;
        }

        let sum = self.add_many(&nonzero_limbs);
        self.target_to_bigint(&BoundedTarget { target: sum, max })
    }

    /// Split the given bounded target into a `BigIntTarget`.
    fn target_to_bigint(
        &mut self,
        input: &BoundedTarget<C::ScalarField>,
    ) -> BigIntTarget<C::ScalarField> {
        struct SplitGenerator<F: Field> {
            input: BoundedTarget<F>,
            output: BigIntTarget<F>,
        }

        impl<F: Field> WitnessGenerator<F> for SplitGenerator<F> {
            fn dependencies(&self) -> Vec<Target<F>> {
                vec![self.input.target]
            }

            fn generate(
                &self,
                _constants: &[Vec<F>],
                witness: &PartialWitness<F>,
            ) -> PartialWitness<F> {
                let value = witness.get_target(self.input.target);

                let mut result = PartialWitness::new();
                result.set_bigint_target(&self.output, &field_to_biguint(value));
                result
            }
        }

        let output = self.add_virtual_bigint_target(&input.max, true);
        self.add_generator(SplitGenerator {
            input: input.clone(),
            output: output.clone(),
        });

        // Check that a weighted sum of the limbs matches the original input.
        let joined = self.bigint_to_target(&output);
        self.copy(joined.target, input.target);

        output
    }

    /// Join a `BigIntTarget` into a `BoundedTarget`.
    fn bigint_to_target(
        &mut self,
        bigint: &BigIntTarget<C::ScalarField>,
    ) -> BoundedTarget<C::ScalarField> {
        let mut sum = self.zero_wire();
        let limb_multiplier = self.constant_wire(C::ScalarField::TWO.exp_usize(LIMB_BITS));
        for &limb in bigint.limbs.iter().rev() {
            sum = self.mul_add(sum, limb_multiplier, limb);
        }
        BoundedTarget {
            target: sum,
            max: bigint.max.clone(),
        }
    }

    pub fn bigint_mul(
        &mut self,
        x: &BigIntTarget<C::ScalarField>,
        y: &BigIntTarget<C::ScalarField>,
    ) -> BigIntTarget<C::ScalarField> {
        let x_n = x.num_limbs();
        let y_n = y.num_limbs();

        let mut result_digits = Vec::new();
        let mut carry: BigIntTarget<C::ScalarField> = BigIntTarget::zero();

        // We want to enumerate the cartesian product of x's and y's limbs, ordered from least
        // significant to most significant. We will do so by enumerating the possible "shifts",
        // measured in limbs. We will add up the (shifted) intermediate products as we go.
        for shift in 0..=(x_n + y_n - 2) {
            let mut sum_of_limb_products = self.bigint_to_target(&carry);

            // We want to enumerate valid pairs of indices such that x_index + y_index = shift.
            for x_index in 0..x_n {
                for y_index in 0..y_n {
                    if x_index + y_index == shift {
                        let x_limb = x.get_bounded_limb(x_index);
                        let y_limb = y.get_bounded_limb(y_index);
                        sum_of_limb_products =
                            self.bounded_mul_add(&x_limb, &y_limb, &sum_of_limb_products);
                    }
                }
            }

            let sum_of_limb_products_bigint = self.target_to_bigint(&sum_of_limb_products);
            let (first, rest) = sum_of_limb_products_bigint.split_smallest_limb();
            result_digits.push(first);
            carry = rest;
        }

        // Add any remaining carry digits.
        for i in 0..carry.num_limbs() {
            result_digits.push(carry.get_limb(i));
        }

        let max = &x.max * &y.max;
        BigIntTarget {
            limbs: result_digits,
            max,
        }
    }

    pub fn bigint_div(
        &mut self,
        x: &BigIntTarget<C::ScalarField>,
        y: &BigIntTarget<C::ScalarField>,
    ) -> BigIntTarget<C::ScalarField> {
        let (div, _rem) = self.bigint_div_rem(x, y);
        div
    }

    pub fn bigint_rem(
        &mut self,
        x: &BigIntTarget<C::ScalarField>,
        y: &BigIntTarget<C::ScalarField>,
    ) -> BigIntTarget<C::ScalarField> {
        let (_div, rem) = self.bigint_div_rem(x, y);
        rem
    }

    /// Returns `(x / y, x % y)`.
    pub fn bigint_div_rem(
        &mut self,
        x: &BigIntTarget<C::ScalarField>,
        y: &BigIntTarget<C::ScalarField>,
    ) -> (BigIntTarget<C::ScalarField>, BigIntTarget<C::ScalarField>) {
        struct DivRemGenerator<F: Field> {
            x: BigIntTarget<F>,
            y: BigIntTarget<F>,
            div: BigIntTarget<F>,
            rem: BigIntTarget<F>,
        }

        impl<F: Field> WitnessGenerator<F> for DivRemGenerator<F> {
            fn dependencies(&self) -> Vec<Target<F>> {
                [self.x.limbs.as_slice(), self.y.limbs.as_slice()].concat()
            }

            fn generate(
                &self,
                _constants: &[Vec<F>],
                witness: &PartialWitness<F>,
            ) -> PartialWitness<F> {
                let x = witness.get_bigint_target(&self.x);
                let y = witness.get_bigint_target(&self.y);
                let (div, rem) = x.div_rem(&y);

                let mut result = PartialWitness::new();
                result.set_bigint_target(&self.div, &div);
                result.set_bigint_target(&self.rem, &rem);
                result
            }
        }

        let max_rem = &y.max - BigUint::one();
        let div = self.add_virtual_bigint_target(&x.max, true);
        let rem = self.add_virtual_bigint_target(&max_rem, true);

        self.add_generator(DivRemGenerator {
            x: x.clone(),
            y: y.clone(),
            div: div.clone(),
            rem: rem.clone(),
        });

        // Check that x = div * y + rem.
        let div_y = self.bigint_mul(&div, y);
        let div_y_plus_rem = self.bigint_add(&div_y, &rem);
        self.copy_bigint(x, &div_y_plus_rem);

        // Check that rem < y.
        let cmp_rem_y = self.bigint_cmp(&rem, y).lt;
        self.assert_one(cmp_rem_y);

        (div, rem)
    }

    /// Assert that the two given bigints encode the same integer.
    pub fn copy_bigint(
        &mut self,
        lhs: &BigIntTarget<C::ScalarField>,
        rhs: &BigIntTarget<C::ScalarField>,
    ) {
        // The number of limbs may differ, in which case we assert equality for any limb indices
        // which are valid for both bigints, then assert that any "extra" limbs (present in one
        // bigint but not the other) are zero.

        let min_limbs = lhs.num_limbs().min(rhs.num_limbs());
        for i in 0..min_limbs {
            self.copy(lhs.get_limb(i), rhs.get_limb(i));
        }

        for i in min_limbs..lhs.num_limbs() {
            self.assert_zero(lhs.get_limb(i));
        }
        for i in min_limbs..rhs.num_limbs() {
            self.assert_zero(rhs.get_limb(i));
        }
    }

    fn bigint_pad_limbs(
        &mut self,
        x: &BigIntTarget<C::ScalarField>,
        num_limbs: usize,
    ) -> BigIntTarget<C::ScalarField> {
        assert!(x.limbs.len() <= num_limbs);
        let mut result = x.clone();
        result.limbs.resize(num_limbs, self.zero_wire());
        result
    }
}

#[cfg(test)]
mod tests {
    use crate::{CircuitBuilder, PartialWitness, Tweedledum};
    use num::{BigUint, FromPrimitive, Integer};

    #[test]
    fn test_bigint_add() {
        let x_value = BigUint::from_u128(22222222222222222222222222222222222222).unwrap();
        let y_value = BigUint::from_u128(33333333333333333333333333333333333333).unwrap();
        let expected_z_value = &x_value + &y_value;

        let mut builder = CircuitBuilder::<Tweedledum>::new(128);
        let x = builder.constant_bigint(&x_value);
        let y = builder.constant_bigint(&y_value);
        let z = builder.bigint_add(&x, &y);
        let circuit = builder.build();

        let witness = circuit.generate_partial_witness(PartialWitness::new());
        let actual_z_value = witness.get_bigint_target(&z);
        assert_eq!(actual_z_value, expected_z_value);
    }

    #[test]
    fn test_bigint_mul() {
        let x_value = BigUint::from_u128(123123123123123123123123123123123123).unwrap();
        let y_value = BigUint::from_u128(456456456456456456456456456456456456).unwrap();
        let expected_z_value = &x_value * &y_value;

        let mut builder = CircuitBuilder::<Tweedledum>::new(128);
        let x = builder.constant_bigint(&x_value);
        let y = builder.constant_bigint(&y_value);
        let z = builder.bigint_mul(&x, &y);
        let circuit = builder.build();

        let witness = circuit.generate_partial_witness(PartialWitness::new());
        let actual_z_value = witness.get_bigint_target(&z);
        assert_eq!(actual_z_value, expected_z_value);
    }

    #[test]
    fn test_bigint_div_rem() {
        let x_value = BigUint::from_u128(456456456456456456456456456456456456).unwrap();
        let y_value = BigUint::from_u128(123123123123123123123123123123123123).unwrap();
        let (expected_div_value, expected_rem_value) = x_value.div_rem(&y_value);

        let mut builder = CircuitBuilder::<Tweedledum>::new(128);
        let x = builder.constant_bigint(&x_value);
        let y = builder.constant_bigint(&y_value);
        let (div, rem) = builder.bigint_div_rem(&x, &y);
        let circuit = builder.build();

        let witness = circuit.generate_partial_witness(PartialWitness::new());
        let actual_div_value = witness.get_bigint_target(&div);
        let actual_rem_value = witness.get_bigint_target(&rem);
        assert_eq!(actual_div_value, expected_div_value);
        assert_eq!(actual_rem_value, expected_rem_value);
    }
}
