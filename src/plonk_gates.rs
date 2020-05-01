//! For reference, here is our gate prefix tree:
//!
//! ```text
//! 00000 <unused> (TODO: Use for "loading" several constants in one gate?)
//! 00001 CurveAddGate
//! 00010 CurveDblGate
//! 00011 CurveEndoGate
//! 00100 Base4SumGate
//! 00101 PublicInputGate
//! 0011* BufferGate
//! 01*** ArithmeticGate
//! 10*** RescueStepAGate
//! 11*** RescueStepBGate
//! ```

use std::marker::PhantomData;

use crate::{AffinePoint, Circuit, CircuitBuilder, Curve, Field, GRID_WIDTH, HaloEndomorphismCurve, NUM_ADVICE_WIRES, NUM_ROUTED_WIRES, NUM_WIRES, PartialWitness, Target, Wire, WitnessGenerator};
use crate::mds::mds;

pub(crate) fn evaluate_all_constraints<C: HaloEndomorphismCurve>(
    local_constant_values: &[C::BaseField],
    local_wire_values: &[C::BaseField],
    right_wire_values: &[C::BaseField],
    below_wire_values: &[C::BaseField],
) -> Vec<C::BaseField> {
    let constraint_sets_per_gate = vec![
        CurveAddGate::<C>::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        CurveDblGate::<C>::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        CurveEndoGate::<C>::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        Base4SumGate::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        PublicInputGate::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        BufferGate::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        ArithmeticGate::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        RescueStepAGate::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        RescueStepBGate::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
    ];

    let mut unified_constraint_set = vec![];
    for constraint_sets in constraint_sets_per_gate {
        while unified_constraint_set.len() < constraint_sets.len() {
            unified_constraint_set.push(C::BaseField::ZERO);
        }
        for i in 0..constraint_sets.len() {
            unified_constraint_set[i] = unified_constraint_set[i] + constraint_sets[i];
        }
    }
    unified_constraint_set
}

pub(crate) fn evaluate_all_constraints_recursively<C: HaloEndomorphismCurve>(
    builder: &mut CircuitBuilder<C::BaseField>,
    local_constant_values: &[Target],
    local_wire_values: &[Target],
    right_wire_values: &[Target],
    below_wire_values: &[Target],
) -> Vec<Target> {
    let constraint_sets_per_gate = vec![
        CurveAddGate::<C>::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        CurveDblGate::<C>::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        CurveEndoGate::<C>::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        Base4SumGate::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        PublicInputGate::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        BufferGate::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        ArithmeticGate::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        RescueStepAGate::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        RescueStepBGate::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
    ];

    let mut unified_constraint_set = vec![];
    for constraint_set in constraint_sets_per_gate {
        while unified_constraint_set.len() < constraint_set.len() {
            unified_constraint_set.push(builder.zero_wire());
        }
        for i in 0..constraint_set.len() {
            unified_constraint_set[i] = builder.add(unified_constraint_set[i], constraint_set[i]);
        }
    }
    unified_constraint_set
}

/// Computes `x * (x - 1)`, which should vanish iff `x` is binary.
fn assert_binary_recursively<F: Field>(builder: &mut CircuitBuilder<F>, x: Target) -> Target {
    let one = builder.one_wire();
    let x_minus_one = builder.sub(x, one);
    builder.mul(x, x_minus_one)
}

/// Computes `x * y - 1`, which should vanish iff `x` and `y` are inverses.
fn assert_inverses_recursively<F: Field>(
    builder: &mut CircuitBuilder<F>,
    x: Target,
    y: Target,
) -> Target {
    let one = builder.one_wire();
    let x_y = builder.mul(x, y);
    builder.sub(x_y, one)
}

pub trait Gate<F: Field>: WitnessGenerator<F> {
    const NAME: &'static str;

    /// In order to combine the constraints of various gate types into a unified constraint set, we
    /// assign each gate type a binary prefix such that no two prefixes overlap.
    const PREFIX: &'static [bool];

    fn evaluate_filtered(local_constant_values: &[F],
                         local_wire_values: &[F],
                         right_wire_values: &[F],
                         below_wire_values: &[F],
    ) -> Vec<F> {
        let filter = Self::evaluate_prefix_filter(local_constant_values);
        let unfiltered = Self::evaluate_unfiltered(
            local_constant_values, local_wire_values, right_wire_values, below_wire_values);
        unfiltered.into_iter().map(|u| filter * u).collect()
    }

    fn evaluate_filtered_recursively(builder: &mut CircuitBuilder<F>,
                                     local_constant_values: &[Target],
                                     local_wire_values: &[Target],
                                     right_wire_values: &[Target],
                                     below_wire_values: &[Target],
    ) -> Vec<Target> {
        let filter = Self::evaluate_prefix_filter_recursively(builder, local_constant_values);
        let unfiltered = Self::evaluate_unfiltered_recursively(
            builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values);
        unfiltered.into_iter().map(|u| builder.mul(filter, u)).collect()
    }

    fn evaluate_prefix_filter(local_constant_values: &[F]) -> F {
        let mut product = F::ONE;
        for (i, &bit) in Self::PREFIX.iter().enumerate() {
            let c = local_constant_values[i];
            if bit {
                product = product * c;
            } else {
                product = product * (F::ONE - c);
            }
        }
        product
    }

    fn evaluate_prefix_filter_recursively(
        builder: &mut CircuitBuilder<F>,
        local_constant_values: &[Target],
    ) -> Target {
        let one = builder.one_wire();
        let mut product = one;
        for (i, &bit) in Self::PREFIX.iter().enumerate() {
            let c = local_constant_values[i];
            let term = if bit {
                c
            } else {
                builder.sub(one, c)
            };
            product = builder.mul(product, term);
        }
        product
    }

    /// Evaluate the constraints implied by this gate at the given challenge point.
    ///
    /// For example, if the gate computes `c = a * b`, this should return `[c(x) - a(x) * b(x)]`,
    /// where `x` is the challenge point.
    fn evaluate_unfiltered(local_constant_values: &[F],
                           local_wire_values: &[F],
                           right_wire_values: &[F],
                           below_wire_values: &[F]) -> Vec<F>;

    /// Like the other `evaluate` method, but in the context of a recursive circuit.
    fn evaluate_unfiltered_recursively(builder: &mut CircuitBuilder<F>,
                                       local_constant_values: &[Target],
                                       local_wire_values: &[Target],
                                       right_wire_values: &[Target],
                                       below_wire_values: &[Target]) -> Vec<Target>;
}

/// A gate for receiving public inputs. These gates will be placed at static indices and the wire
/// polynomials will always be opened at those indices.
///
/// Because our gate arity is 11 but only 6 of the wires are routed, it may seem as though each gate
/// can only receive 6 public inputs. To work around this, we place a BufferGate immediately after
/// each PublicInputGate, and have the PublicInputGate copy its 5 non-routed wires to routed wires
/// of the BufferGate.
pub(crate) struct PublicInputGate {
    pub index: usize,
    /// Make the constructor private.
    _private: (),
}

impl PublicInputGate {
    pub fn new(index: usize) -> Self {
        PublicInputGate { index, _private: () }
    }
}

impl<F: Field> Gate<F> for PublicInputGate {
    const NAME: &'static str = "PublicInputGate";

    const PREFIX: &'static [bool] = &[false, false, true, false, true];

    fn evaluate_unfiltered(
        local_constant_values: &[F],
        local_wire_values: &[F],
        right_wire_values: &[F],
        below_wire_values: &[F],
    ) -> Vec<F> {
        // This ensures that advice wires' values are copied to the following buffer gate.
        // TODO: Consider enforcing this via copy constraints, in which case there would be nothing to do here.
        (0..NUM_ADVICE_WIRES).map(|i| {
            local_wire_values[NUM_ROUTED_WIRES + i] - right_wire_values[i]
        }).collect()
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<F>,
        local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        below_wire_values: &[Target],
    ) -> Vec<Target> {
        // This ensures that advice wires' values are copied to the following buffer gate.
        // TODO: Consider enforcing this via copy constraints, in which case there would be nothing to do here.
        (0..NUM_ADVICE_WIRES).map(|i| {
            builder.sub(local_wire_values[NUM_ROUTED_WIRES + i], right_wire_values[i])
        }).collect()
    }
}

impl<F: Field> WitnessGenerator<F> for PublicInputGate {
    fn dependencies(&self) -> Vec<Target> {
        (0..NUM_WIRES)
            .map(|i| Target::Wire(Wire { gate: self.index, input: i }))
            .collect()
    }

    fn generate(&self, circuit: Circuit<F>, witness: &PartialWitness<F>) -> PartialWitness<F> {
        let self_as_generator: &dyn WitnessGenerator<F> = self;
        let targets: Vec<Target> = self_as_generator.dependencies();

        let mut result = PartialWitness::new();
        for i_advice in 0..NUM_ADVICE_WIRES {
            let i_wire = NUM_ROUTED_WIRES + i_advice;
            let value = witness.get_target(targets[i_wire]);
            result.set_wire(Wire { gate: self.index + 1, input: i_advice }, value);
        }
        result
    }
}

/// A gate which doesn't perform any arithmetic, but just acts as a buffer for receiving data. This
/// is used in a couple ways:
/// * Some gates, such as the Rescue round gate, "output" their results using one of the next gate's
///   "input" wires. The last such gate has no next gate of the same type, so we add a buffer gate
///   for receiving the last gate's output.
/// * The first constant value configured for this gate will be proxied to its `WIRE_BUFFER_CONST`
///   wire; this allows us to create routable constant wires.
pub(crate) struct BufferGate {
    pub index: usize,
    /// Make the constructor private.
    _private: (),
}

impl BufferGate {
    pub fn new(index: usize) -> Self {
        BufferGate { index, _private: () }
    }

    pub const WIRE_BUFFER_0: usize = 0;
    pub const WIRE_BUFFER_1: usize = 1;
    pub const WIRE_BUFFER_2: usize = 2;
    pub const WIRE_BUFFER_3: usize = 3;
    pub const WIRE_BUFFER_4: usize = 4;
    pub const WIRE_BUFFER_CONST: usize = 5;
}

impl<F: Field> Gate<F> for BufferGate {
    const NAME: &'static str = "BufferGate";

    const PREFIX: &'static [bool] = &[false, false, true, true];

    fn evaluate_unfiltered(
        local_constant_values: &[F],
        local_wire_values: &[F],
        right_wire_values: &[F],
        below_wire_values: &[F],
    ) -> Vec<F> {
        let wire_value = local_wire_values[Self::WIRE_BUFFER_CONST];
        let const_value = local_constant_values[<Self as Gate<F>>::PREFIX.len()];
        vec![wire_value - const_value]
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<F>,
        local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        below_wire_values: &[Target],
    ) -> Vec<Target> {
        let wire_value = local_wire_values[Self::WIRE_BUFFER_CONST];
        let const_value = local_constant_values[<Self as Gate<F>>::PREFIX.len()];
        vec![builder.sub(wire_value, const_value)]
    }
}

impl<F: Field> WitnessGenerator<F> for BufferGate {
    fn dependencies(&self) -> Vec<Target> {
        Vec::new()
    }

    fn generate(&self, circuit: Circuit<F>, _witness: &PartialWitness<F>) -> PartialWitness<F> {
        let buffer_const_target = Wire { gate: self.index, input: Self::WIRE_BUFFER_CONST };

        let mut witness = PartialWitness::new();
        let const_value = circuit.gate_constants[self.index][<Self as Gate<F>>::PREFIX.len()];
        witness.set_wire(buffer_const_target, const_value);
        witness
    }
}

/// A gate which performs incomplete point addition, conditioned on an input bit. In order to
/// facilitate MSMs which use this gate, it also adds the bit to an accumulator.
///
/// `C` is the curve whose points are being added.
pub(crate) struct CurveAddGate<C: Curve> {
    pub index: usize,
    _phantom: PhantomData<C>,
}

impl<C: Curve> CurveAddGate<C> {
    pub fn new(index: usize) -> Self {
        CurveAddGate { index, _phantom: PhantomData }
    }

    pub const WIRE_GROUP_ACC_X: usize = 0;
    pub const WIRE_GROUP_ACC_Y: usize = 1;
    pub const WIRE_SCALAR_ACC_OLD: usize = 2;
    pub const WIRE_SCALAR_ACC_NEW: usize = 3;
    pub const WIRE_ADDEND_X: usize = 4;
    pub const WIRE_ADDEND_Y: usize = 5;
    pub const WIRE_SCALAR_BIT: usize = 6;
    pub const WIRE_INVERSE: usize = 7;
}

impl<C: Curve> Gate<C::BaseField> for CurveAddGate<C> {
    const NAME: &'static str = "CurveAddGate";

    const PREFIX: &'static [bool] = &[false, false, false, false, true];

    fn evaluate_unfiltered(
        local_constant_values: &[C::BaseField],
        local_wire_values: &[C::BaseField],
        right_wire_values: &[C::BaseField],
        below_wire_values: &[C::BaseField],
    ) -> Vec<C::BaseField> {
        let x1 = local_wire_values[Self::WIRE_GROUP_ACC_X];
        let y1 = local_wire_values[Self::WIRE_GROUP_ACC_Y];
        let x3 = right_wire_values[Self::WIRE_GROUP_ACC_X];
        let y3 = right_wire_values[Self::WIRE_GROUP_ACC_Y];

        let scalar_acc_old = local_wire_values[Self::WIRE_SCALAR_ACC_OLD];
        let scalar_acc_new = local_wire_values[Self::WIRE_SCALAR_ACC_NEW];
        let x2 = local_wire_values[Self::WIRE_ADDEND_X];
        let y2 = local_wire_values[Self::WIRE_ADDEND_Y];
        let scalar_bit = local_wire_values[Self::WIRE_SCALAR_BIT];
        let inverse = local_wire_values[Self::WIRE_INVERSE];

        let lambda = (y1 - y2) * inverse;
        let computed_x3 = lambda.square() - x1 - x2;
        let computed_y3 = lambda * (x1 - x3) - y1;

        vec![
            computed_x3 - x3,
            computed_y3 - y3,
            scalar_acc_new - scalar_acc_old.double() + scalar_bit,
            scalar_bit * (scalar_bit - C::BaseField::ONE),
            inverse * (x1 - x2) - C::BaseField::ONE,
        ]
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C::BaseField>,
        local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        below_wire_values: &[Target],
    ) -> Vec<Target> {
        let x1 = local_wire_values[Self::WIRE_GROUP_ACC_X];
        let y1 = local_wire_values[Self::WIRE_GROUP_ACC_Y];
        let x3 = right_wire_values[Self::WIRE_GROUP_ACC_X];
        let y3 = right_wire_values[Self::WIRE_GROUP_ACC_Y];

        let scalar_acc_old = local_wire_values[Self::WIRE_SCALAR_ACC_OLD];
        let scalar_acc_new = local_wire_values[Self::WIRE_SCALAR_ACC_NEW];
        let x2 = local_wire_values[Self::WIRE_ADDEND_X];
        let y2 = local_wire_values[Self::WIRE_ADDEND_Y];
        let scalar_bit = local_wire_values[Self::WIRE_SCALAR_BIT];
        let inverse = local_wire_values[Self::WIRE_INVERSE];

        let x1_minus_x2 = builder.sub(x1, x2);
        let x1_plus_x2 = builder.sub(x1, x2);
        let x1_minus_x3 = builder.sub(x1, x3);
        let y1_minus_y2 = builder.sub(y1, y2);

        let lambda = builder.mul(y1_minus_y2, inverse);
        let lambda_squared = builder.square(lambda);
        let computed_x3 = builder.sub(lambda_squared, x1_plus_x2);
        let computed_y3 = builder.mul_sub(lambda, x1_minus_x3, y1);

        let double_scalar_acc_old = builder.double(scalar_acc_old);
        let computed_scalar_acc_new = builder.add(double_scalar_acc_old, scalar_bit);

        vec![
            builder.sub(computed_x3, x3),
            builder.sub(computed_y3, y3),
            builder.sub(computed_scalar_acc_new, scalar_acc_new),
            assert_binary_recursively(builder, scalar_bit),
            assert_inverses_recursively(builder, inverse, x1_minus_x2),
        ]
    }
}

impl<C: Curve> WitnessGenerator<C::BaseField> for CurveAddGate<C> {
    fn dependencies(&self) -> Vec<Target> {
        vec![
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_GROUP_ACC_X }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_GROUP_ACC_Y }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_SCALAR_ACC_OLD }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_ADDEND_X }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_ADDEND_Y }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_SCALAR_BIT }),
        ]
    }

    fn generate(&self, circuit: Circuit<C::BaseField>, witness: &PartialWitness<C::BaseField>) -> PartialWitness<C::BaseField> {
        let group_acc_old_x_target = Wire { gate: self.index, input: Self::WIRE_GROUP_ACC_X };
        let group_acc_new_x_target = Wire { gate: self.index + 1, input: Self::WIRE_GROUP_ACC_X };
        let group_acc_old_y_target = Wire { gate: self.index, input: Self::WIRE_GROUP_ACC_Y };
        let group_acc_new_y_target = Wire { gate: self.index + 1, input: Self::WIRE_GROUP_ACC_Y };
        let scalar_acc_old_target = Wire { gate: self.index, input: Self::WIRE_SCALAR_ACC_OLD };
        let scalar_acc_new_target = Wire { gate: self.index, input: Self::WIRE_SCALAR_ACC_NEW };
        let addend_x_target = Wire { gate: self.index, input: Self::WIRE_ADDEND_X };
        let addend_y_target = Wire { gate: self.index, input: Self::WIRE_ADDEND_Y };
        let scalar_bit_target = Wire { gate: self.index, input: Self::WIRE_SCALAR_BIT };
        let inverse_target = Wire { gate: self.index, input: Self::WIRE_INVERSE };

        let group_acc_old_x = witness.get_wire(group_acc_old_x_target);
        let group_acc_old_y = witness.get_wire(group_acc_old_y_target);
        let group_acc_old = AffinePoint::<C>::nonzero(group_acc_old_x, group_acc_old_y);

        let scalar_acc_old = witness.get_wire(scalar_acc_old_target);

        let addend_x = witness.get_wire(addend_x_target);
        let addend_y = witness.get_wire(addend_y_target);
        let addend = AffinePoint::<C>::nonzero(addend_x, addend_y);

        let scalar_bit = witness.get_wire(scalar_bit_target);
        debug_assert!(scalar_bit.is_zero() || scalar_bit.is_one());

        let mut group_acc_new = group_acc_old;
        if scalar_bit.is_one() {
            group_acc_new = (group_acc_new + addend).to_affine();
        }

        let scalar_acc_new = scalar_acc_old.double() + scalar_bit;

        // Here's where our abstraction leaks a bit. Although we already have the sum, we need to
        // redo part of the computation in order to populate the purported inverse wire.
        let dx = group_acc_old_x - addend_x;
        let inverse = dx.multiplicative_inverse().expect("x_1 = x_2");

        let mut result = PartialWitness::new();
        result.set_wire(group_acc_new_x_target, group_acc_new.x);
        result.set_wire(group_acc_new_y_target, group_acc_new.y);
        result.set_wire(scalar_acc_new_target, scalar_acc_new);
        result.set_wire(inverse_target, inverse);
        result
    }
}

/// A curve which performs point doubling.
pub(crate) struct CurveDblGate<C: Curve> {
    pub index: usize,
    _phantom: PhantomData<C>,
}

impl<C: Curve> CurveDblGate<C> {
    pub fn new(index: usize) -> Self {
        CurveDblGate { index, _phantom: PhantomData }
    }

    pub const WIRE_X_OLD: usize = 0;
    pub const WIRE_Y_OLD: usize = 1;
    pub const WIRE_X_NEW: usize = 2;
    pub const WIRE_Y_NEW: usize = 3;
    pub const WIRE_INVERSE: usize = 4;
}

impl<C: Curve> Gate<C::BaseField> for CurveDblGate<C> {
    const NAME: &'static str = "CurveDblGate";

    const PREFIX: &'static [bool] = &[false, false, false, true, false];

    fn evaluate_unfiltered(
        local_constant_values: &[C::BaseField],
        local_wire_values: &[C::BaseField],
        right_wire_values: &[C::BaseField],
        below_wire_values: &[C::BaseField],
    ) -> Vec<C::BaseField> {
        let x_old = local_wire_values[Self::WIRE_X_OLD];
        let y_old = local_wire_values[Self::WIRE_Y_OLD];
        let x_new = local_wire_values[Self::WIRE_X_NEW];
        let y_new = local_wire_values[Self::WIRE_Y_NEW];
        let inverse = local_wire_values[Self::WIRE_INVERSE];

        let lambda_numerator = x_old.square().triple() + C::A;
        let lambda = lambda_numerator * inverse;
        let computed_x_new = lambda.square() - x_old.double();
        let computed_y_new = lambda * (x_old - x_new) - y_old;

        vec![
            // Verify that computed_x_new matches x_new.
            computed_x_new - x_new,
            // Verify that computed_y_new matches y_new.
            computed_y_new - y_new,
            // Verify that 2 * y_old times its purported inverse is 1.
            y_old.double() * inverse - C::BaseField::ONE,
        ]
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C::BaseField>,
        local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        below_wire_values: &[Target],
    ) -> Vec<Target> {
        let one = builder.one_wire();
        let three = builder.constant_wire_u32(3);
        let a = builder.constant_wire(C::A);

        let x_old = local_wire_values[Self::WIRE_X_OLD];
        let y_old = local_wire_values[Self::WIRE_Y_OLD];
        let x_new = local_wire_values[Self::WIRE_X_NEW];
        let y_new = local_wire_values[Self::WIRE_Y_NEW];
        let inverse = local_wire_values[Self::WIRE_INVERSE];

        let two_x_old = builder.double(x_old);
        let two_y_old = builder.double(y_old);
        let x_old_squared = builder.square(x_old);
        let three_x_old_squared = builder.mul(three, x_old_squared);
        let lambda_numerator = builder.add(three_x_old_squared, a);
        let lambda = builder.mul(lambda_numerator, inverse);
        let lambda_squared = builder.square(lambda);
        let computed_x_new = builder.sub(lambda_squared, two_x_old);
        let delta_x = builder.sub(x_old, x_new);
        let lambda_times_delta_x = builder.mul(lambda, delta_x);
        let computed_y_new = builder.sub(lambda_times_delta_x, y_old);

        vec![
            // Verify that computed_x_new matches x_new.
            builder.sub(computed_x_new, x_new),
            // Verify that computed_y_new matches y_new.
            builder.sub(computed_y_new, y_new),
            // Verify that 2 * y_old times its purported inverse is 1.
            assert_inverses_recursively(builder, two_y_old, inverse),
        ]
    }
}

impl<C: Curve> WitnessGenerator<C::BaseField> for CurveDblGate<C> {
    fn dependencies(&self) -> Vec<Target> {
        vec![
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_X_OLD }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_Y_OLD }),
        ]
    }

    fn generate(&self, circuit: Circuit<C::BaseField>, witness: &PartialWitness<C::BaseField>) -> PartialWitness<C::BaseField> {
        let x_old_target = Wire { gate: self.index, input: Self::WIRE_X_OLD };
        let y_old_target = Wire { gate: self.index, input: Self::WIRE_Y_OLD };
        let x_new_target = Wire { gate: self.index, input: Self::WIRE_X_NEW };
        let y_new_target = Wire { gate: self.index, input: Self::WIRE_Y_NEW };
        let inverse_target = Wire { gate: self.index, input: Self::WIRE_INVERSE };

        let x_old = witness.get_wire(x_old_target);
        let y_old = witness.get_wire(y_old_target);
        let old = AffinePoint::<C>::nonzero(x_old, y_old);
        let new = old.double();

        // Here's where our abstraction leaks a bit. Although we already have the result, we need to
        // redo part of the computation in order to populate the purported inverse wire.
        let inverse = y_old.double().multiplicative_inverse().expect("y = 0");

        let mut result = PartialWitness::new();
        result.set_wire(inverse_target, inverse);
        result.set_wire(x_new_target, new.x);
        result.set_wire(y_new_target, new.y);
        result
    }
}

/// A gate which performs an iteration of an simultaneous doubling MSM loop, employing the
/// endomorphism described in the Halo paper. `C` is the curve of the inner proof.
pub(crate) struct CurveEndoGate<C: HaloEndomorphismCurve> {
    pub index: usize,
    _phantom: PhantomData<C>,
}

impl<C: HaloEndomorphismCurve> CurveEndoGate<C> {
    pub fn new(index: usize) -> Self {
        CurveEndoGate { index, _phantom: PhantomData }
    }

    pub const WIRE_GROUP_ACC_X: usize = 0;
    pub const WIRE_GROUP_ACC_Y: usize = 1;
    pub const WIRE_SCALAR_ACC_UNSIGNED: usize = 2;
    pub const WIRE_SCALAR_ACC_SIGNED: usize = 3;
    pub const WIRE_ADDEND_X: usize = 4;
    pub const WIRE_ADDEND_Y: usize = 5;
    pub const WIRE_SCALAR_BIT_0: usize = 6;
    pub const WIRE_SCALAR_BIT_1: usize = 7;
    pub const WIRE_INVERSE: usize = 8;
}

impl<C: HaloEndomorphismCurve> Gate<C::BaseField> for CurveEndoGate<C> {
    const NAME: &'static str = "CurveEndoGate";

    const PREFIX: &'static [bool] = &[false, false, false, true, true];

    fn evaluate_unfiltered(
        local_constant_values: &[C::BaseField],
        local_wire_values: &[C::BaseField],
        right_wire_values: &[C::BaseField],
        below_wire_values: &[C::BaseField],
    ) -> Vec<C::BaseField> {
        unimplemented!()
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C::BaseField>,
        local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        below_wire_values: &[Target],
    ) -> Vec<Target> {
        unimplemented!()
    }
}

impl<C: HaloEndomorphismCurve> WitnessGenerator<C::BaseField> for CurveEndoGate<C> {
    fn dependencies(&self) -> Vec<Target> {
        vec![
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_GROUP_ACC_X }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_GROUP_ACC_Y }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_SCALAR_ACC_UNSIGNED }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_SCALAR_ACC_SIGNED }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_ADDEND_X }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_ADDEND_Y }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_SCALAR_BIT_0 }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_SCALAR_BIT_1 }),
        ]
    }

    fn generate(&self, circuit: Circuit<C::BaseField>, witness: &PartialWitness<C::BaseField>) -> PartialWitness<C::BaseField> {
        let group_acc_old_x_target = Wire { gate: self.index, input: Self::WIRE_GROUP_ACC_X };
        let group_acc_new_x_target = Wire { gate: self.index + 1, input: Self::WIRE_GROUP_ACC_X };
        let group_acc_old_y_target = Wire { gate: self.index, input: Self::WIRE_GROUP_ACC_Y };
        let group_acc_new_y_target = Wire { gate: self.index + 1, input: Self::WIRE_GROUP_ACC_Y };

        let scalar_acc_unsigned_old_target = Wire { gate: self.index, input: Self::WIRE_GROUP_ACC_Y };
        let scalar_acc_unsigned_new_target = Wire { gate: self.index + GRID_WIDTH, input: Self::WIRE_GROUP_ACC_Y };
        let scalar_acc_signed_old_target = Wire { gate: self.index, input: Self::WIRE_GROUP_ACC_Y };
        let scalar_acc_signed_new_target = Wire { gate: self.index + GRID_WIDTH, input: Self::WIRE_GROUP_ACC_Y };

        let addend_x_target = Wire { gate: self.index, input: Self::WIRE_ADDEND_X };
        let addend_y_target = Wire { gate: self.index, input: Self::WIRE_ADDEND_Y };
        let scalar_bit_0_target = Wire { gate: self.index, input: Self::WIRE_SCALAR_BIT_0 };
        let scalar_bit_1_target = Wire { gate: self.index, input: Self::WIRE_SCALAR_BIT_1 };
        let inverse_target = Wire { gate: self.index, input: Self::WIRE_INVERSE };

        let group_acc_old_x = witness.get_wire(group_acc_old_x_target);
        let group_acc_old_y = witness.get_wire(group_acc_old_y_target);
        let group_acc_old = AffinePoint::<C>::nonzero(group_acc_old_x, group_acc_old_y);

        let scalar_acc_unsigned_old = witness.get_wire(scalar_acc_unsigned_old_target);
        let scalar_acc_signed_old = witness.get_wire(scalar_acc_signed_old_target);

        let scalar_bit_0 = witness.get_wire(scalar_bit_0_target);
        let scalar_bit_1 = witness.get_wire(scalar_bit_1_target);
        debug_assert!(scalar_bit_0.is_zero() || scalar_bit_0.is_one());
        debug_assert!(scalar_bit_1.is_zero() || scalar_bit_1.is_one());

        let p_x = witness.get_wire(addend_x_target);
        let p_y = witness.get_wire(addend_y_target);

        let mut s_i_x = p_x;
        if scalar_bit_0 == C::BaseField::ONE {
            s_i_x = s_i_x * C::ZETA;
        }
        let mut s_i_y = p_y;
        if scalar_bit_1 == C::BaseField::ZERO {
            s_i_y = -s_i_y;
        }
        let s_i = AffinePoint::nonzero(s_i_x, s_i_y);
        let group_acc_new = group_acc_old + s_i;

        let scalar_acc_unsigned_new = scalar_acc_unsigned_old.quadruple()
            + scalar_bit_0 + scalar_bit_1.double();

        // This is based on Algorithm 2 in the Halo paper.
        let mut scalar_acc_signed_limb = if scalar_bit_0 == C::BaseField::ONE {
            C::BaseField::ONE
        } else {
            C::BaseField::NEG_ONE
        };
        if scalar_bit_1 == C::BaseField::ONE {
            scalar_acc_signed_limb = scalar_acc_signed_limb * C::ZETA;
        }
        let scalar_acc_signed_new = scalar_acc_signed_old.double() + scalar_acc_signed_limb;

        // Here's where our abstraction leaks a bit. Although we already have the sum, we need to
        // redo part of the computation in order to populate the purported inverse wire.
        let dx = group_acc_old_x - p_x;
        let dy = group_acc_old_y - p_y;
        let inverse = dx.multiplicative_inverse().expect("x_1 = x_2");

        let mut result = PartialWitness::new();
        result.set_wire(group_acc_new_x_target, group_acc_new.x);
        result.set_wire(group_acc_new_y_target, group_acc_new.y);
        result.set_wire(scalar_acc_unsigned_new_target, scalar_acc_unsigned_new);
        result.set_wire(scalar_acc_signed_new_target, scalar_acc_signed_new);
        result.set_wire(inverse_target, inverse);

        result
    }
}

/// The first step of Rescue, i.e. the one with the `x^(1/5)` layer.
pub(crate) struct RescueStepAGate<F: Field> {
    pub index: usize,
    _phantom: PhantomData<F>,
}

impl<F: Field> RescueStepAGate<F> {
    pub fn new(index: usize) -> Self {
        RescueStepAGate { index, _phantom: PhantomData }
    }

    pub const WIRE_INPUT_0: usize = 0;
    pub const WIRE_INPUT_1: usize = 1;
    pub const WIRE_INPUT_2: usize = 2;
    pub const WIRE_OUTPUT_0: usize = 3;
    pub const WIRE_OUTPUT_1: usize = 4;
    pub const WIRE_OUTPUT_2: usize = 5;
    pub const WIRE_ROOT_0: usize = 6;
    pub const WIRE_ROOT_1: usize = 7;
    pub const WIRE_ROOT_2: usize = 8;
}

impl<F: Field> Gate<F> for RescueStepAGate<F> {
    const NAME: &'static str = "RescueStepAGate";

    const PREFIX: &'static [bool] = &[true, false];

    fn evaluate_unfiltered(
        local_constant_values: &[F],
        local_wire_values: &[F],
        right_wire_values: &[F],
        below_wire_values: &[F],
    ) -> Vec<F> {
        let in_0 = local_wire_values[Self::WIRE_INPUT_0];
        let in_1 = local_wire_values[Self::WIRE_INPUT_1];
        let in_2 = local_wire_values[Self::WIRE_INPUT_2];
        let out_0 = local_wire_values[Self::WIRE_OUTPUT_0];
        let out_1 = local_wire_values[Self::WIRE_OUTPUT_1];
        let out_2 = local_wire_values[Self::WIRE_OUTPUT_2];
        let root_0 = local_wire_values[Self::WIRE_ROOT_0];
        let root_1 = local_wire_values[Self::WIRE_ROOT_1];
        let root_2 = local_wire_values[Self::WIRE_ROOT_2];

        let computed_out_0 = mds::<F>(3, 0, 0) * root_0 + mds::<F>(3, 0, 1) * root_1 + mds::<F>(3, 0, 2) * root_2 + local_constant_values[Self::PREFIX.len()];
        let computed_out_1 = mds::<F>(3, 1, 0) * root_0 + mds::<F>(3, 1, 1) * root_1 + mds::<F>(3, 1, 2) * root_2 + local_constant_values[Self::PREFIX.len() + 1];
        let computed_out_2 = mds::<F>(3, 2, 0) * root_0 + mds::<F>(3, 2, 1) * root_1 + mds::<F>(3, 2, 2) * root_2 + local_constant_values[Self::PREFIX.len() + 2];

        vec![
            root_0.exp_u32(5) - in_0,
            root_1.exp_u32(5) - in_1,
            root_2.exp_u32(5) - in_2,
            computed_out_0 - out_0,
            computed_out_1 - out_1,
            computed_out_2 - out_2,
        ]
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<F>,
        local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        below_wire_values: &[Target],
    ) -> Vec<Target> {
        let in_0 = local_wire_values[Self::WIRE_INPUT_0];
        let in_1 = local_wire_values[Self::WIRE_INPUT_1];
        let in_2 = local_wire_values[Self::WIRE_INPUT_2];
        let out_0 = local_wire_values[Self::WIRE_OUTPUT_0];
        let out_1 = local_wire_values[Self::WIRE_OUTPUT_1];
        let out_2 = local_wire_values[Self::WIRE_OUTPUT_2];
        let root_0 = local_wire_values[Self::WIRE_ROOT_0];
        let root_1 = local_wire_values[Self::WIRE_ROOT_1];
        let root_2 = local_wire_values[Self::WIRE_ROOT_2];

        let computed_in_0 = builder.exp_constant_usize(root_0, 5);
        let computed_in_1 = builder.exp_constant_usize(root_1, 5);
        let computed_in_2 = builder.exp_constant_usize(root_2, 5);

        let mds_00 = builder.constant_wire(mds::<F>(3, 0, 0));
        let mds_01 = builder.constant_wire(mds::<F>(3, 0, 1));
        let mds_02 = builder.constant_wire(mds::<F>(3, 0, 2));
        let mds_10 = builder.constant_wire(mds::<F>(3, 1, 0));
        let mds_11 = builder.constant_wire(mds::<F>(3, 1, 1));
        let mds_12 = builder.constant_wire(mds::<F>(3, 1, 2));
        let mds_20 = builder.constant_wire(mds::<F>(3, 2, 0));
        let mds_21 = builder.constant_wire(mds::<F>(3, 2, 1));
        let mds_22 = builder.constant_wire(mds::<F>(3, 2, 2));

        let mds_00_root_0 = builder.mul(mds_00, root_0);
        let mds_01_root_1 = builder.mul(mds_01, root_1);
        let mds_02_root_2 = builder.mul(mds_02, root_2);
        let mds_10_root_0 = builder.mul(mds_10, root_0);
        let mds_11_root_1 = builder.mul(mds_11, root_1);
        let mds_12_root_2 = builder.mul(mds_12, root_2);
        let mds_20_root_0 = builder.mul(mds_20, root_0);
        let mds_21_root_1 = builder.mul(mds_21, root_1);
        let mds_22_root_2 = builder.mul(mds_22, root_2);

        let computed_out_0 = builder.add_many(&[mds_00_root_0, mds_01_root_1, mds_02_root_2, local_constant_values[Self::PREFIX.len()]]);
        let computed_out_1 = builder.add_many(&[mds_10_root_0, mds_11_root_1, mds_12_root_2, local_constant_values[Self::PREFIX.len() + 1]]);
        let computed_out_2 = builder.add_many(&[mds_20_root_0, mds_21_root_1, mds_22_root_2, local_constant_values[Self::PREFIX.len() + 2]]);

        vec![
            builder.sub(computed_in_0, in_0),
            builder.sub(computed_in_1, in_1),
            builder.sub(computed_in_2, in_2),
            builder.sub(computed_out_0, out_0),
            builder.sub(computed_out_1, out_1),
            builder.sub(computed_out_2, out_2),
        ]
    }
}

impl<F: Field> WitnessGenerator<F> for RescueStepAGate<F> {
    fn dependencies(&self) -> Vec<Target> {
        vec![
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_INPUT_0 }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_INPUT_1 }),
        ]
    }

    fn generate(&self, circuit: Circuit<F>, witness: &PartialWitness<F>) -> PartialWitness<F> {
        let constants = &circuit.gate_constants[self.index];

        let in_0_target = Wire { gate: self.index, input: Self::WIRE_INPUT_0 };
        let in_1_target = Wire { gate: self.index, input: Self::WIRE_INPUT_1 };
        let in_2_target = Wire { gate: self.index, input: Self::WIRE_INPUT_2 };

        let root_0_target = Wire { gate: self.index, input: Self::WIRE_ROOT_0 };
        let root_1_target = Wire { gate: self.index, input: Self::WIRE_ROOT_1 };
        let root_2_target = Wire { gate: self.index, input: Self::WIRE_ROOT_2 };

        let out_0_target = Wire { gate: self.index + 1, input: Self::WIRE_INPUT_0 };
        let out_1_target = Wire { gate: self.index + 1, input: Self::WIRE_INPUT_1 };
        let out_2_target = Wire { gate: self.index + 1, input: Self::WIRE_INPUT_2 };

        let in_0 = witness.get_wire(in_0_target);
        let in_1 = witness.get_wire(in_1_target);
        let in_2 = witness.get_wire(in_2_target);

        let root_0 = in_0.kth_root_u32(5);
        let root_1 = in_1.kth_root_u32(5);
        let root_2 = in_2.kth_root_u32(5);

        let out_0 = mds::<F>(3, 0, 0) * root_0 + mds::<F>(3, 0, 1) * root_1 + mds::<F>(3, 0, 2) * root_2 + constants[Self::PREFIX.len()];
        let out_1 = mds::<F>(3, 1, 0) * root_0 + mds::<F>(3, 1, 1) * root_1 + mds::<F>(3, 1, 2) * root_2 + constants[Self::PREFIX.len() + 1];
        let out_2 = mds::<F>(3, 2, 0) * root_0 + mds::<F>(3, 2, 1) * root_1 + mds::<F>(3, 2, 2) * root_2 + constants[Self::PREFIX.len() + 2];

        let mut result = PartialWitness::new();
        result.set_wire(root_0_target, root_0);
        result.set_wire(root_1_target, root_1);
        result.set_wire(root_2_target, root_2);
        result.set_wire(out_0_target, out_0);
        result.set_wire(out_1_target, out_1);
        result.set_wire(out_2_target, out_2);
        result
    }
}

/// The second step of Rescue, i.e. the one with the `x^5` layer.
pub(crate) struct RescueStepBGate<F: Field> {
    pub index: usize,
    _phantom: PhantomData<F>,
}

impl<F: Field> RescueStepBGate<F> {
    pub fn new(index: usize) -> Self {
        RescueStepBGate { index, _phantom: PhantomData }
    }

    pub const WIRE_INPUT_0: usize = 0;
    pub const WIRE_INPUT_1: usize = 1;
    pub const WIRE_INPUT_2: usize = 2;
    pub const WIRE_OUTPUT_0: usize = 3;
    pub const WIRE_OUTPUT_1: usize = 4;
    pub const WIRE_OUTPUT_2: usize = 5;
}

impl<F: Field> Gate<F> for RescueStepBGate<F> {
    const NAME: &'static str = "RescueStepBGate";

    const PREFIX: &'static [bool] = &[true, true];

    fn evaluate_unfiltered(
        local_constant_values: &[F],
        local_wire_values: &[F],
        right_wire_values: &[F],
        below_wire_values: &[F],
    ) -> Vec<F> {
        let in_0 = local_wire_values[Self::WIRE_INPUT_0];
        let in_1 = local_wire_values[Self::WIRE_INPUT_1];
        let in_2 = local_wire_values[Self::WIRE_INPUT_2];
        let out_0 = local_wire_values[Self::WIRE_OUTPUT_0];
        let out_1 = local_wire_values[Self::WIRE_OUTPUT_1];
        let out_2 = local_wire_values[Self::WIRE_OUTPUT_2];

        let in_0_exp_5 = in_0.exp_u32(5);
        let in_1_exp_5 = in_1.exp_u32(5);
        let in_2_exp_5 = in_2.exp_u32(5);

        let computed_out_0 = mds::<F>(3, 0, 0) * in_0_exp_5 + mds::<F>(3, 0, 1) * in_1_exp_5 + mds::<F>(3, 0, 2) * in_2_exp_5 + local_constant_values[Self::PREFIX.len()];
        let computed_out_1 = mds::<F>(3, 1, 0) * in_0_exp_5 + mds::<F>(3, 1, 1) * in_1_exp_5 + mds::<F>(3, 1, 2) * in_2_exp_5 + local_constant_values[Self::PREFIX.len() + 1];
        let computed_out_2 = mds::<F>(3, 2, 0) * in_0_exp_5 + mds::<F>(3, 2, 1) * in_1_exp_5 + mds::<F>(3, 2, 2) * in_2_exp_5 + local_constant_values[Self::PREFIX.len() + 2];

        vec![
            computed_out_0 - out_0,
            computed_out_1 - out_1,
            computed_out_2 - out_2,
        ]
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<F>,
        local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        below_wire_values: &[Target],
    ) -> Vec<Target> {
        let in_0 = local_wire_values[Self::WIRE_INPUT_0];
        let in_1 = local_wire_values[Self::WIRE_INPUT_1];
        let in_2 = local_wire_values[Self::WIRE_INPUT_2];
        let out_0 = local_wire_values[Self::WIRE_OUTPUT_0];
        let out_1 = local_wire_values[Self::WIRE_OUTPUT_1];
        let out_2 = local_wire_values[Self::WIRE_OUTPUT_2];

        let in_0_exp_5 = builder.exp_constant_usize(in_0, 5);
        let in_1_exp_5 = builder.exp_constant_usize(in_1, 5);
        let in_2_exp_5 = builder.exp_constant_usize(in_2, 5);

        let mds_00 = builder.constant_wire(mds::<F>(3, 0, 0));
        let mds_01 = builder.constant_wire(mds::<F>(3, 0, 1));
        let mds_02 = builder.constant_wire(mds::<F>(3, 0, 2));
        let mds_10 = builder.constant_wire(mds::<F>(3, 1, 0));
        let mds_11 = builder.constant_wire(mds::<F>(3, 1, 1));
        let mds_12 = builder.constant_wire(mds::<F>(3, 1, 2));
        let mds_20 = builder.constant_wire(mds::<F>(3, 2, 0));
        let mds_21 = builder.constant_wire(mds::<F>(3, 2, 1));
        let mds_22 = builder.constant_wire(mds::<F>(3, 2, 2));

        let mds_00_in_0_exp_5 = builder.mul(mds_00, in_0_exp_5);
        let mds_01_in_1_exp_5 = builder.mul(mds_01, in_1_exp_5);
        let mds_02_in_2_exp_5 = builder.mul(mds_02, in_2_exp_5);
        let mds_10_in_0_exp_5 = builder.mul(mds_10, in_0_exp_5);
        let mds_11_in_1_exp_5 = builder.mul(mds_11, in_1_exp_5);
        let mds_12_in_2_exp_5 = builder.mul(mds_12, in_2_exp_5);
        let mds_20_in_0_exp_5 = builder.mul(mds_20, in_0_exp_5);
        let mds_21_in_1_exp_5 = builder.mul(mds_21, in_1_exp_5);
        let mds_22_in_2_exp_5 = builder.mul(mds_22, in_2_exp_5);

        let computed_out_0 = builder.add_many(&[mds_00_in_0_exp_5, mds_01_in_1_exp_5, mds_02_in_2_exp_5, local_constant_values[Self::PREFIX.len()]]);
        let computed_out_1 = builder.add_many(&[mds_10_in_0_exp_5, mds_11_in_1_exp_5, mds_12_in_2_exp_5, local_constant_values[Self::PREFIX.len() + 1]]);
        let computed_out_2 = builder.add_many(&[mds_20_in_0_exp_5, mds_21_in_1_exp_5, mds_22_in_2_exp_5, local_constant_values[Self::PREFIX.len() + 2]]);

        vec![
            builder.sub(computed_out_0, out_0),
            builder.sub(computed_out_1, out_1),
            builder.sub(computed_out_2, out_2),
        ]
    }
}

impl<F: Field> WitnessGenerator<F> for RescueStepBGate<F> {
    fn dependencies(&self) -> Vec<Target> {
        vec![
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_INPUT_0 }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_INPUT_1 }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_INPUT_2 }),
        ]
    }

    fn generate(&self, circuit: Circuit<F>, witness: &PartialWitness<F>) -> PartialWitness<F> {
        let constants = &circuit.gate_constants[self.index];

        let in_0_target = Wire { gate: self.index, input: Self::WIRE_INPUT_0 };
        let in_1_target = Wire { gate: self.index, input: Self::WIRE_INPUT_1 };
        let in_2_target = Wire { gate: self.index, input: Self::WIRE_INPUT_2 };

        let out_0_target = Wire { gate: self.index + 1, input: Self::WIRE_INPUT_0 };
        let out_1_target = Wire { gate: self.index + 1, input: Self::WIRE_INPUT_1 };
        let out_2_target = Wire { gate: self.index + 1, input: Self::WIRE_INPUT_2 };

        let in_0 = witness.get_wire(in_0_target);
        let in_1 = witness.get_wire(in_1_target);
        let in_2 = witness.get_wire(in_2_target);

        let exp_0 = in_0.exp_u32(5);
        let exp_1 = in_1.exp_u32(5);
        let exp_2 = in_2.exp_u32(5);

        let out_0 = mds::<F>(3, 0, 0) * exp_0 + mds::<F>(3, 0, 1) * exp_1 + mds::<F>(3, 0, 2) * exp_2 + constants[Self::PREFIX.len()];
        let out_1 = mds::<F>(3, 1, 0) * exp_0 + mds::<F>(3, 1, 1) * exp_1 + mds::<F>(3, 1, 2) * exp_2 + constants[Self::PREFIX.len() + 1];
        let out_2 = mds::<F>(3, 2, 0) * exp_0 + mds::<F>(3, 2, 1) * exp_1 + mds::<F>(3, 2, 2) * exp_2 + constants[Self::PREFIX.len() + 2];

        let mut result = PartialWitness::new();
        result.set_wire(out_0_target, out_0);
        result.set_wire(out_1_target, out_1);
        result.set_wire(out_2_target, out_2);
        result
    }
}

/// A gate for accumulating base-4 limbs.
pub(crate) struct Base4SumGate {
    pub index: usize,
    /// Make the constructor private.
    _private: (),
}

impl Base4SumGate {
    pub fn new(index: usize) -> Self {
        Base4SumGate { index, _private: () }
    }

    pub const WIRE_ACC_OLD: usize = 0;
    pub const WIRE_ACC_NEW: usize = 1;
    pub const WIRE_LIMB_0: usize = 2;
    pub const NUM_LIMBS: usize = 9;
}

impl<F: Field> Gate<F> for Base4SumGate {
    const NAME: &'static str = "Base4SumGate";

    // TODO: Need to reduce this prefix length by 1 bit in order to keep everything within degree 8n.
    const PREFIX: &'static [bool] = &[false, false, true, false, false];

    fn evaluate_unfiltered(
        local_constant_values: &[F],
        local_wire_values: &[F],
        right_wire_values: &[F],
        below_wire_values: &[F],
    ) -> Vec<F> {
        let acc_old = local_wire_values[Self::WIRE_ACC_OLD];
        let acc_new = local_wire_values[Self::WIRE_ACC_NEW];
        let limbs: Vec<F> = (0..Self::NUM_LIMBS).map(
            |i| local_wire_values[Self::WIRE_LIMB_0 + i]
        ).collect();

        let mut computed_acc_new = acc_old;
        for &limb in &limbs {
            computed_acc_new = computed_acc_new.quadruple() + limb;
        }

        let mut constraints = vec![computed_acc_new - acc_new];
        for limb in limbs {
            let mut product = F::ONE;
            for j in 0..4 {
                product = product * (limb - F::from_canonical_usize(j));
            }
            constraints.push(product);
        }
        constraints
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<F>,
        local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        below_wire_values: &[Target],
    ) -> Vec<Target> {
        let four = builder.constant_wire_u32(4);

        let acc_old = local_wire_values[Self::WIRE_ACC_OLD];
        let acc_new = local_wire_values[Self::WIRE_ACC_NEW];
        let limbs: Vec<Target> = (0..Self::NUM_LIMBS).map(
            |i| local_wire_values[Self::WIRE_LIMB_0 + i]
        ).collect();

        let mut computed_acc_new = acc_old;
        for &limb in &limbs {
            let shifted_acc_new = builder.mul(computed_acc_new, four);
            computed_acc_new = builder.add(shifted_acc_new, limb);
        }

        let mut constraints = vec![builder.sub(computed_acc_new, acc_new)];
        for limb in limbs {
            let mut product = builder.one_wire();
            for j in 0..4 {
                let j_target = builder.constant_wire_u32(j);
                let limb_minus_j = builder.sub(limb, j_target);
                product = builder.mul(product, limb_minus_j);
            }
            constraints.push(product);
        }
        constraints
    }
}

impl<F: Field> WitnessGenerator<F> for Base4SumGate {
    fn dependencies(&self) -> Vec<Target> {
        Vec::new()
    }

    fn generate(&self, circuit: Circuit<F>, _witness: &PartialWitness<F>) -> PartialWitness<F> {
        // For base 4 decompositions, we don't do any witness generation on a per-gate level.
        // Instead, we have a single generator which generates values for an entire decomposition.
        PartialWitness::new()
    }
}

/// A gate which can be configured to perform various arithmetic. In particular, it computes
///
/// ```text
/// output := const_0 * multiplicand_0 * multiplicand_1 + const_1 * addend + const_2
/// ```
pub(crate) struct ArithmeticGate<F: Field> {
    pub index: usize,
    _phantom: PhantomData<F>,
}

impl<F: Field> ArithmeticGate<F> {
    pub fn new(index: usize) -> Self {
        ArithmeticGate { index, _phantom: PhantomData }
    }

    pub const WIRE_MULTIPLICAND_0: usize = 0;
    pub const WIRE_MULTIPLICAND_1: usize = 1;
    pub const WIRE_ADDEND: usize = 2;
    pub const WIRE_OUTPUT: usize = 3;
}

impl<F: Field> Gate<F> for ArithmeticGate<F> {
    const NAME: &'static str = "ArithmeticGate";

    const PREFIX: &'static [bool] = &[false, true];

    fn evaluate_unfiltered(
        local_constant_values: &[F],
        local_wire_values: &[F],
        right_wire_values: &[F],
        below_wire_values: &[F],
    ) -> Vec<F> {
        unimplemented!()
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<F>,
        local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        below_wire_values: &[Target],
    ) -> Vec<Target> {
        unimplemented!()
    }
}

impl<F: Field> WitnessGenerator<F> for ArithmeticGate<F> {
    fn dependencies(&self) -> Vec<Target> {
        vec![
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_MULTIPLICAND_0 }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_MULTIPLICAND_1 }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_ADDEND }),
        ]
    }

    fn generate(&self, circuit: Circuit<F>, witness: &PartialWitness<F>) -> PartialWitness<F> {
        let multiplicand_0_target = Wire { gate: self.index, input: Self::WIRE_MULTIPLICAND_0 };
        let multiplicand_1_target = Wire { gate: self.index, input: Self::WIRE_MULTIPLICAND_1 };
        let addend_target = Wire { gate: self.index, input: Self::WIRE_ADDEND };
        let output_target = Wire { gate: self.index, input: Self::WIRE_OUTPUT };

        let const_0 = circuit.gate_constants[self.index][Self::PREFIX.len()];
        let const_1 = circuit.gate_constants[self.index][Self::PREFIX.len() + 1];
        let const_2 = circuit.gate_constants[self.index][Self::PREFIX.len() + 2];

        let multiplicand_0 = witness.get_wire(multiplicand_0_target);
        let multiplicand_1 = witness.get_wire(multiplicand_1_target);
        let addend = witness.get_wire(addend_target);

        let output = const_0 * multiplicand_0 * multiplicand_1 + const_1 * addend + const_2;

        let mut result = PartialWitness::new();
        result.set_wire(output_target, output);
        result
    }
}
