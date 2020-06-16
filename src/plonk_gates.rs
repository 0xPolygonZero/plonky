//! For reference, here is our gate prefix tree:
//!
//! ```text
//! 101001 PublicInputGate
//! 101000 CurveAddGate
//! 10111* CurveDblGate
//! 11**** CurveEndoGate
//! 1000** Base4SumGate
//! 101010 BufferGate
//! 10110* ConstantGate
//! 1001** ArithmeticGate
//! 00**** RescueStepAGate
//! 01**** RescueStepBGate
//! ```
//!
//! The `*`s above represent constants which are not used in the gate prefix, and are thus available
//! for gate configuration.

use std::marker::PhantomData;

use crate::{AffinePoint, CircuitBuilder, Curve, Field, GRID_WIDTH, HaloCurve, NUM_ADVICE_WIRES, NUM_ROUTED_WIRES, NUM_WIRES, PartialWitness, Target, Wire, WitnessGenerator, mds_matrix};

pub const RESCUE_SPONGE_WIDTH: usize = 4;

pub fn evaluate_all_constraints<C: HaloCurve, InnerC: HaloCurve<BaseField=C::ScalarField>>(
    local_constant_values: &[C::ScalarField],
    local_wire_values: &[C::ScalarField],
    right_wire_values: &[C::ScalarField],
    below_wire_values: &[C::ScalarField],
) -> Vec<C::ScalarField> {
    let constraint_sets_per_gate = vec![
        CurveAddGate::<C, InnerC>::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        CurveDblGate::<C, InnerC>::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        CurveEndoGate::<C, InnerC>::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        Base4SumGate::<C>::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        PublicInputGate::<C>::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        BufferGate::<C>::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        ConstantGate::<C>::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        ArithmeticGate::<C>::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        RescueStepAGate::<C>::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        RescueStepBGate::<C>::evaluate_filtered(local_constant_values, local_wire_values, right_wire_values, below_wire_values),
    ];

    let mut unified_constraint_set = vec![];
    for constraint_sets in constraint_sets_per_gate {
        while unified_constraint_set.len() < constraint_sets.len() {
            unified_constraint_set.push(C::ScalarField::ZERO);
        }
        for i in 0..constraint_sets.len() {
            unified_constraint_set[i] = unified_constraint_set[i] + constraint_sets[i];
        }
    }
    unified_constraint_set
}

pub fn evaluate_all_constraints_recursively<C: HaloCurve, InnerC: HaloCurve<BaseField=C::ScalarField>>(
    builder: &mut CircuitBuilder<C>,
    local_constant_values: &[Target],
    local_wire_values: &[Target],
    right_wire_values: &[Target],
    below_wire_values: &[Target],
) -> Vec<Target> {
    let constraint_sets_per_gate = vec![
        CurveAddGate::<C, InnerC>::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        CurveDblGate::<C, InnerC>::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        CurveEndoGate::<C, InnerC>::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        Base4SumGate::<C>::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        PublicInputGate::<C>::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        BufferGate::<C>::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        ConstantGate::<C>::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        ArithmeticGate::<C>::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        RescueStepAGate::<C>::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
        RescueStepBGate::<C>::evaluate_filtered_recursively(builder, local_constant_values, local_wire_values, right_wire_values, below_wire_values),
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
fn assert_binary_recursively<C: HaloCurve>(builder: &mut CircuitBuilder<C>, x: Target) -> Target {
    let one = builder.one_wire();
    let x_minus_one = builder.sub(x, one);
    builder.mul(x, x_minus_one)
}

/// Computes `x * y - 1`, which should vanish iff `x` and `y` are inverses.
fn assert_inverses_recursively<C: HaloCurve>(
    builder: &mut CircuitBuilder<C>,
    x: Target,
    y: Target,
) -> Target {
    let one = builder.one_wire();
    let x_y = builder.mul(x, y);
    builder.sub(x_y, one)
}

pub trait Gate<C: HaloCurve>: WitnessGenerator<C::ScalarField> {
    const NAME: &'static str;

    /// In order to combine the constraints of various gate types into a unified constraint set, we
    /// assign each gate type a binary prefix such that no two prefixes overlap.
    const PREFIX: &'static [bool];

    fn evaluate_filtered(local_constant_values: &[C::ScalarField],
                         local_wire_values: &[C::ScalarField],
                         right_wire_values: &[C::ScalarField],
                         below_wire_values: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        let filter = Self::evaluate_prefix_filter(local_constant_values);
        let unfiltered = Self::evaluate_unfiltered(
            local_constant_values, local_wire_values, right_wire_values, below_wire_values);
        unfiltered.into_iter().map(|u| filter * u).collect()
    }

    fn evaluate_filtered_recursively(builder: &mut CircuitBuilder<C>,
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

    fn evaluate_prefix_filter(local_constant_values: &[C::ScalarField]) -> C::ScalarField {
        let mut product = C::ScalarField::ONE;
        for (i, &bit) in Self::PREFIX.iter().enumerate() {
            let c = local_constant_values[i];
            if bit {
                product = product * c;
            } else {
                product = product * (C::ScalarField::ONE - c);
            }
        }
        product
    }

    fn evaluate_prefix_filter_recursively(
        builder: &mut CircuitBuilder<C>,
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
    fn evaluate_unfiltered(local_constant_values: &[C::ScalarField],
                           local_wire_values: &[C::ScalarField],
                           right_wire_values: &[C::ScalarField],
                           below_wire_values: &[C::ScalarField]) -> Vec<C::ScalarField>;

    /// Like the other `evaluate` method, but in the context of a recursive circuit.
    fn evaluate_unfiltered_recursively(builder: &mut CircuitBuilder<C>,
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
pub struct PublicInputGate<C: HaloCurve> {
    pub index: usize,
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> PublicInputGate<C> {
    pub fn new(index: usize) -> Self {
        PublicInputGate { index, _phantom: PhantomData }
    }
}

impl<C: HaloCurve> Gate<C> for PublicInputGate<C> {
    const NAME: &'static str = "PublicInputGate";

    const PREFIX: &'static [bool] = &[true, false, true, false, false, true];

    fn evaluate_unfiltered(
        _local_constant_values: &[C::ScalarField],
        local_wire_values: &[C::ScalarField],
        right_wire_values: &[C::ScalarField],
        _below_wire_values: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        // This ensures that advice wires' values are copied to the following buffer gate.
        // TODO: Consider enforcing this via copy constraints, in which case there would be nothing to do here.
        (0..NUM_ADVICE_WIRES).map(|i| {
            local_wire_values[NUM_ROUTED_WIRES + i] - right_wire_values[i]
        }).collect()
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        _local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        _below_wire_values: &[Target],
    ) -> Vec<Target> {
        // This ensures that advice wires' values are copied to the following buffer gate.
        // TODO: Consider enforcing this via copy constraints, in which case there would be nothing to do here.
        (0..NUM_ADVICE_WIRES).map(|i| {
            builder.sub(local_wire_values[NUM_ROUTED_WIRES + i], right_wire_values[i])
        }).collect()
    }
}

impl<C: HaloCurve> WitnessGenerator<C::ScalarField> for PublicInputGate<C> {
    fn dependencies(&self) -> Vec<Target> {
        (0..NUM_WIRES)
            .map(|i| Target::Wire(Wire { gate: self.index, input: i }))
            .collect()
    }

    fn generate(&self, _constants: &Vec<Vec<C::ScalarField>>, witness: &PartialWitness<C::ScalarField>) -> PartialWitness<C::ScalarField> {
        let self_as_generator: &dyn WitnessGenerator<C::ScalarField> = self;
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

/// A gate which doesn't perform any arithmetic, but just acts as a buffer for receiving data.
/// Some gates, such as the Rescue round gate, "output" their results using one of the next gate's
/// "input" wires. The last such gate has no next gate of the same type, so we add a buffer gate
/// for receiving the last gate's output.
pub struct BufferGate<C: HaloCurve> {
    pub index: usize,
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> BufferGate<C> {
    pub fn new(index: usize) -> Self {
        BufferGate { index, _phantom: PhantomData }
    }
}

impl<C: HaloCurve> Gate<C> for BufferGate<C> {
    const NAME: &'static str = "BufferGate";

    const PREFIX: &'static [bool] = &[true, false, true, false, true, false];

    fn evaluate_unfiltered(
        _local_constant_values: &[C::ScalarField],
        _local_wire_values: &[C::ScalarField],
        _right_wire_values: &[C::ScalarField],
        _below_wire_values: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        Vec::new()
    }

    fn evaluate_unfiltered_recursively(
        _builder: &mut CircuitBuilder<C>,
        _local_constant_values: &[Target],
        _local_wire_values: &[Target],
        _right_wire_values: &[Target],
        _below_wire_values: &[Target],
    ) -> Vec<Target> {
        Vec::new()
    }
}

impl<C: HaloCurve> WitnessGenerator<C::ScalarField> for BufferGate<C> {
    fn dependencies(&self) -> Vec<Target> {
        Vec::new()
    }

    fn generate(&self, _constants: &Vec<Vec<C::ScalarField>>, _witness: &PartialWitness<C::ScalarField>) -> PartialWitness<C::ScalarField> {
        PartialWitness::new()
    }
}

/// A gate which takes a single constant parameter and outputs that value.
pub struct ConstantGate<C: HaloCurve> {
    pub index: usize,
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> ConstantGate<C> {
    pub const WIRE_OUTPUT: usize = 0;

    pub fn new(index: usize) -> Self {
        ConstantGate { index, _phantom: PhantomData }
    }
}

impl<C: HaloCurve> Gate<C> for ConstantGate<C> {
    const NAME: &'static str = "ConstantGate";

    const PREFIX: &'static [bool] = &[true, false, true, true, false];

    fn evaluate_unfiltered(
        local_constant_values: &[C::ScalarField],
        local_wire_values: &[C::ScalarField],
        _right_wire_values: &[C::ScalarField],
        _below_wire_values: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        let c = local_constant_values[Self::PREFIX.len()];
        let out = local_wire_values[Self::WIRE_OUTPUT];
        vec![c - out]
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        local_constant_values: &[Target],
        local_wire_values: &[Target],
        _right_wire_values: &[Target],
        _below_wire_values: &[Target],
    ) -> Vec<Target> {
        let c = local_constant_values[Self::PREFIX.len()];
        let out = local_wire_values[Self::WIRE_OUTPUT];
        vec![builder.sub(c, out)]
    }
}

impl<C: HaloCurve> WitnessGenerator<C::ScalarField> for ConstantGate<C> {
    fn dependencies(&self) -> Vec<Target> {
        Vec::new()
    }

    fn generate(
        &self,
        constants: &Vec<Vec<C::ScalarField>>,
        _witness: &PartialWitness<C::ScalarField>,
    ) -> PartialWitness<C::ScalarField> {
        let constants = &constants[self.index];
        let c = constants[Self::PREFIX.len()];
        let mut result = PartialWitness::new();
        result.set_wire(Wire { gate: self.index, input: Self::WIRE_OUTPUT }, c);
        result
    }
}

/// A gate which performs incomplete point addition, conditioned on an input bit. In order to
/// facilitate MSMs which use this gate, it also adds the bit to an accumulator.
///
/// `C` is the curve whose points are being added.
pub struct CurveAddGate<C: HaloCurve, InnerC: Curve<BaseField=C::ScalarField>> {
    pub index: usize,
    _phantom_oc: PhantomData<C>,
    _phantom_ic: PhantomData<InnerC>,
}

impl<C: HaloCurve, InnerC: Curve<BaseField=C::ScalarField>>
CurveAddGate<C, InnerC> {
    pub fn new(index: usize) -> Self {
        CurveAddGate { index, _phantom_oc: PhantomData, _phantom_ic: PhantomData }
    }

    pub const WIRE_GROUP_ACC_X: usize = 0;
    pub const WIRE_GROUP_ACC_Y: usize = 1;
    pub const WIRE_SCALAR_ACC_OLD: usize = 2;
    pub const WIRE_SCALAR_ACC_NEW: usize = 3;
    pub const WIRE_ADDEND_X: usize = 4;
    pub const WIRE_ADDEND_Y: usize = 5;
    pub const WIRE_SCALAR_BIT: usize = 6;
    pub const WIRE_INVERSE: usize = 7;
    pub const WIRE_LAMBDA: usize = 8;
}

impl<C: HaloCurve, InnerC: Curve<BaseField=C::ScalarField>>
Gate<C> for CurveAddGate<C, InnerC> {
    const NAME: &'static str = "CurveAddGate";

    const PREFIX: &'static [bool] = &[true, false, true, false, false, false];

    fn evaluate_unfiltered(
        _local_constant_values: &[InnerC::BaseField],
        local_wire_values: &[InnerC::BaseField],
        right_wire_values: &[InnerC::BaseField],
        _below_wire_values: &[InnerC::BaseField],
    ) -> Vec<InnerC::BaseField> {
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
        let lambda = local_wire_values[Self::WIRE_LAMBDA];

        let computed_lambda = (y1 - y2) * inverse;
        let computed_x3 = lambda.square() - x1 - x2;
        let computed_y3 = lambda * (x1 - x3) - y1;

        vec![
            computed_lambda - lambda,
            computed_x3 - x3,
            computed_y3 - y3,
            scalar_acc_new - scalar_acc_old.double() + scalar_bit,
            scalar_bit * (scalar_bit - InnerC::BaseField::ONE),
            inverse * (x1 - x2) - InnerC::BaseField::ONE,
        ]
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        _local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        _below_wire_values: &[Target],
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
        let lambda = local_wire_values[Self::WIRE_LAMBDA];

        let x1_minus_x2 = builder.sub(x1, x2);
        let x1_plus_x2 = builder.add(x1, x2);
        let x1_minus_x3 = builder.sub(x1, x3);
        let y1_minus_y2 = builder.sub(y1, y2);

        let computed_lambda = builder.mul(y1_minus_y2, inverse);
        let computed_x3 = builder.mul_sub(lambda, lambda, x1_plus_x2);
        let computed_y3 = builder.mul_sub(lambda, x1_minus_x3, y1);

        let double_scalar_acc_old = builder.double(scalar_acc_old);
        let computed_scalar_acc_new = builder.add(double_scalar_acc_old, scalar_bit);

        vec![
            builder.sub(computed_lambda, lambda),
            builder.sub(computed_x3, x3),
            builder.sub(computed_y3, y3),
            builder.sub(computed_scalar_acc_new, scalar_acc_new),
            assert_binary_recursively(builder, scalar_bit),
            assert_inverses_recursively(builder, inverse, x1_minus_x2),
        ]
    }
}

impl<C: HaloCurve, InnerC: Curve<BaseField=C::ScalarField>>
WitnessGenerator<C::ScalarField> for CurveAddGate<C, InnerC> {
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

    fn generate(&self, _constants: &Vec<Vec<C::ScalarField>>, witness: &PartialWitness<InnerC::BaseField>) -> PartialWitness<InnerC::BaseField> {
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
        let lambda_target = Wire { gate: self.index, input: Self::WIRE_LAMBDA };

        let group_acc_old_x = witness.get_wire(group_acc_old_x_target);
        let group_acc_old_y = witness.get_wire(group_acc_old_y_target);

        let scalar_acc_old = witness.get_wire(scalar_acc_old_target);

        let addend_x = witness.get_wire(addend_x_target);
        let addend_y = witness.get_wire(addend_y_target);

        let scalar_bit = witness.get_wire(scalar_bit_target);
        debug_assert!(scalar_bit.is_zero() || scalar_bit.is_one());

        let scalar_acc_new = scalar_acc_old.double() + scalar_bit;

        // Here's where our abstraction leaks a bit. Although we already have the sum, we need to
        // redo part of the computation in order to populate the purported inverse wire.
        let dx = group_acc_old_x - addend_x;
        let dy = group_acc_old_y - addend_y;
        let inverse = dx.multiplicative_inverse().expect("x_1 = x_2");
        let lambda = dy * inverse;

        let x_3 = lambda.square() - group_acc_old_x - addend_x;
        let y_3 = lambda * (group_acc_old_x - x_3) - group_acc_old_y;
        let (group_acc_new_x, group_acc_new_y) = if scalar_bit.is_one() {
            (x_3, y_3)
        } else {
            (group_acc_old_x, group_acc_old_y)
        };

        let mut result = PartialWitness::new();
        result.set_wire(group_acc_new_x_target, group_acc_new_x);
        result.set_wire(group_acc_new_y_target, group_acc_new_y);
        result.set_wire(scalar_acc_new_target, scalar_acc_new);
        result.set_wire(inverse_target, inverse);
        result.set_wire(lambda_target, lambda);
        result
    }
}

/// A curve which performs point doubling.
pub struct CurveDblGate<C: HaloCurve, InnerC: Curve<BaseField=C::ScalarField>> {
    pub index: usize,
    _phantom_oc: PhantomData<C>,
    _phantom_ic: PhantomData<InnerC>,
}

impl<C: HaloCurve, InnerC: Curve<BaseField=C::ScalarField>>
CurveDblGate<C, InnerC> {
    pub fn new(index: usize) -> Self {
        CurveDblGate { index, _phantom_oc: PhantomData, _phantom_ic: PhantomData }
    }

    pub const WIRE_X_OLD: usize = 0;
    pub const WIRE_Y_OLD: usize = 1;
    pub const WIRE_X_NEW: usize = 2;
    pub const WIRE_Y_NEW: usize = 3;
    pub const WIRE_INVERSE: usize = 4;
    pub const WIRE_LAMBDA: usize = 5;
}

impl<C: HaloCurve, InnerC: Curve<BaseField=C::ScalarField>>
Gate<C> for CurveDblGate<C, InnerC> {
    const NAME: &'static str = "CurveDblGate";

    const PREFIX: &'static [bool] = &[true, false, true, true, true];

    fn evaluate_unfiltered(
        _local_constant_values: &[InnerC::BaseField],
        local_wire_values: &[InnerC::BaseField],
        _right_wire_values: &[InnerC::BaseField],
        _below_wire_values: &[InnerC::BaseField],
    ) -> Vec<InnerC::BaseField> {
        let x_old = local_wire_values[Self::WIRE_X_OLD];
        let y_old = local_wire_values[Self::WIRE_Y_OLD];
        let x_new = local_wire_values[Self::WIRE_X_NEW];
        let y_new = local_wire_values[Self::WIRE_Y_NEW];
        let inverse = local_wire_values[Self::WIRE_INVERSE];
        let lambda = local_wire_values[Self::WIRE_LAMBDA];

        let computed_lambda_numerator = x_old.square().triple() + InnerC::A;
        let computed_lambda = computed_lambda_numerator * inverse;
        let computed_x_new = lambda.square() - x_old.double();
        let computed_y_new = lambda * (x_old - x_new) - y_old;

        vec![
            // Verify that computed_lambda matches lambda.
            computed_lambda - lambda,
            // Verify that computed_x_new matches x_new.
            computed_x_new - x_new,
            // Verify that computed_y_new matches y_new.
            computed_y_new - y_new,
            // Verify that 2 * y_old times its purported inverse is 1.
            y_old.double() * inverse - InnerC::BaseField::ONE,
        ]
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        _local_constant_values: &[Target],
        local_wire_values: &[Target],
        _right_wire_values: &[Target],
        _below_wire_values: &[Target],
    ) -> Vec<Target> {
        let _one = builder.one_wire();
        let three = builder.constant_wire_u32(3);
        let a = builder.constant_wire(InnerC::A);

        let x_old = local_wire_values[Self::WIRE_X_OLD];
        let y_old = local_wire_values[Self::WIRE_Y_OLD];
        let x_new = local_wire_values[Self::WIRE_X_NEW];
        let y_new = local_wire_values[Self::WIRE_Y_NEW];
        let inverse = local_wire_values[Self::WIRE_INVERSE];
        let lambda = local_wire_values[Self::WIRE_LAMBDA];

        let two_x_old = builder.double(x_old);
        let two_y_old = builder.double(y_old);
        let x_old_squared = builder.square(x_old);
        let three_x_old_squared = builder.mul(three, x_old_squared);
        let computed_lambda_numerator = builder.add(three_x_old_squared, a);
        let computed_lambda = builder.mul(computed_lambda_numerator, inverse);
        let lambda_squared = builder.square(lambda);
        let computed_x_new = builder.sub(lambda_squared, two_x_old);
        let delta_x = builder.sub(x_old, x_new);
        let lambda_times_delta_x = builder.mul(lambda, delta_x);
        let computed_y_new = builder.sub(lambda_times_delta_x, y_old);

        vec![
            // Verify that computed_lambda matches lambda.
            builder.sub(computed_lambda, lambda),
            // Verify that computed_x_new matches x_new.
            builder.sub(computed_x_new, x_new),
            // Verify that computed_y_new matches y_new.
            builder.sub(computed_y_new, y_new),
            // Verify that 2 * y_old times its purported inverse is 1.
            assert_inverses_recursively(builder, two_y_old, inverse),
        ]
    }
}

impl<C: HaloCurve, InnerC: Curve<BaseField=C::ScalarField>>
WitnessGenerator<C::ScalarField> for CurveDblGate<C, InnerC> {
    fn dependencies(&self) -> Vec<Target> {
        vec![
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_X_OLD }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_Y_OLD }),
        ]
    }

    fn generate(&self, _constants: &Vec<Vec<C::ScalarField>>, witness: &PartialWitness<InnerC::BaseField>) -> PartialWitness<InnerC::BaseField> {
        let x_old_target = Wire { gate: self.index, input: Self::WIRE_X_OLD };
        let y_old_target = Wire { gate: self.index, input: Self::WIRE_Y_OLD };
        let x_new_target = Wire { gate: self.index, input: Self::WIRE_X_NEW };
        let y_new_target = Wire { gate: self.index, input: Self::WIRE_Y_NEW };
        let inverse_target = Wire { gate: self.index, input: Self::WIRE_INVERSE };
        let lambda_target = Wire { gate: self.index, input: Self::WIRE_LAMBDA };

        let x_old = witness.get_wire(x_old_target);
        let y_old = witness.get_wire(y_old_target);

        let inverse = y_old.double().multiplicative_inverse().expect("y = 0");
        let lambda = x_old.square().triple() * inverse;
        let x_new = lambda.square() - x_old.double();
        let y_new = lambda * (x_old - x_new) - y_old;

        let mut result = PartialWitness::new();
        result.set_wire(inverse_target, inverse);
        result.set_wire(lambda_target, lambda);
        result.set_wire(x_new_target, x_new);
        result.set_wire(y_new_target, y_new);
        result
    }
}

/// A gate which performs an iteration of an simultaneous doubling MSM loop, employing the
/// endomorphism described in the Halo paper. `C` is the curve of the inner proof.
pub struct CurveEndoGate<C: HaloCurve, InnerC: HaloCurve<BaseField=C::ScalarField>> {
    pub index: usize,
    _phantom_oc: PhantomData<C>,
    _phantom_ic: PhantomData<InnerC>,
}

impl<C: HaloCurve, InnerC: HaloCurve<BaseField=C::ScalarField>>
CurveEndoGate<C, InnerC> {
    pub fn new(index: usize) -> Self {
        CurveEndoGate { index, _phantom_oc: PhantomData, _phantom_ic: PhantomData }
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

impl<C: HaloCurve, InnerC: HaloCurve<BaseField=C::ScalarField>>
Gate<C> for CurveEndoGate<C, InnerC> {
    const NAME: &'static str = "CurveEndoGate";

    const PREFIX: &'static [bool] = &[true, true];

    fn evaluate_unfiltered(
        _local_constant_values: &[InnerC::BaseField],
        local_wire_values: &[InnerC::BaseField],
        right_wire_values: &[InnerC::BaseField],
        below_wire_values: &[InnerC::BaseField],
    ) -> Vec<InnerC::BaseField> {
        let one = InnerC::BaseField::ONE;

        let x1 = local_wire_values[Self::WIRE_GROUP_ACC_X];
        let y1 = local_wire_values[Self::WIRE_GROUP_ACC_Y];
        let x_in = local_wire_values[Self::WIRE_ADDEND_X];
        let y_in = local_wire_values[Self::WIRE_ADDEND_Y];
        let x3 = right_wire_values[Self::WIRE_GROUP_ACC_X];
        let y3 = right_wire_values[Self::WIRE_GROUP_ACC_Y];
        let scalar_acc_unsigned_old = local_wire_values[Self::WIRE_SCALAR_ACC_UNSIGNED];
        let scalar_acc_unsigned_new = below_wire_values[Self::WIRE_SCALAR_ACC_UNSIGNED];
        let scalar_acc_signed_old = local_wire_values[Self::WIRE_SCALAR_ACC_SIGNED];
        let scalar_acc_signed_new = below_wire_values[Self::WIRE_SCALAR_ACC_SIGNED];
        let scalar_bit_0 = local_wire_values[Self::WIRE_SCALAR_BIT_0];
        let scalar_bit_1 = local_wire_values[Self::WIRE_SCALAR_BIT_1];
        let inverse = local_wire_values[Self::WIRE_INVERSE];

        // Conditionally apply the endo and conditionally negate in order to get (x2, y2), which is
        // the actual point we want to add to the accumulator.
        let x2 = ((InnerC::ZETA - one) * scalar_bit_1 + one) * x_in;
        let y2 = (scalar_bit_0.double() - one) * y_in;

        let lambda = (y1 - y2) * inverse;
        let computed_x3 = lambda.square() - x1 - x2;
        let computed_y3 = lambda * (x1 - x3) - y1;

        // This is based on Algorithm 2 in the Halo paper.
        let signed_limb_multiplier = (InnerC::ZETA - one) * scalar_bit_1 + one;
        let signed_limb = (scalar_bit_0.double() - one) * signed_limb_multiplier;

        vec![
            computed_x3 - x3,
            computed_y3 - y3,
            scalar_acc_unsigned_new - scalar_acc_unsigned_old.quadruple() + scalar_bit_1.double() + scalar_bit_0,
            scalar_acc_signed_new - scalar_acc_signed_old.double() + signed_limb,
            scalar_bit_0 * (scalar_bit_0 - one),
            scalar_bit_1 * (scalar_bit_1 - one),
            inverse * (x1 - x2) - one,
        ]
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        _local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        below_wire_values: &[Target],
    ) -> Vec<Target> {
        let one = builder.one_wire();
        let two = builder.two_wire();
        let four = builder.constant_wire_u32(4);
        let zeta = builder.constant_wire(InnerC::ZETA);

        let x1 = local_wire_values[Self::WIRE_GROUP_ACC_X];
        let y1 = local_wire_values[Self::WIRE_GROUP_ACC_Y];
        let x_in = local_wire_values[Self::WIRE_ADDEND_X];
        let y_in = local_wire_values[Self::WIRE_ADDEND_Y];
        let x3 = right_wire_values[Self::WIRE_GROUP_ACC_X];
        let y3 = right_wire_values[Self::WIRE_GROUP_ACC_Y];
        let scalar_acc_unsigned_old = local_wire_values[Self::WIRE_SCALAR_ACC_UNSIGNED];
        let scalar_acc_unsigned_new = below_wire_values[Self::WIRE_SCALAR_ACC_UNSIGNED];
        let scalar_acc_signed_old = local_wire_values[Self::WIRE_SCALAR_ACC_SIGNED];
        let scalar_acc_signed_new = below_wire_values[Self::WIRE_SCALAR_ACC_SIGNED];
        let scalar_bit_0 = local_wire_values[Self::WIRE_SCALAR_BIT_0];
        let scalar_bit_1 = local_wire_values[Self::WIRE_SCALAR_BIT_1];
        let inverse = local_wire_values[Self::WIRE_INVERSE];

        // Conditionally apply the endo and conditionally negate in order to get (x2, y2), which is
        // the actual point we want to add to the accumulator.
        let zeta_minus_one = builder.sub(zeta, one);
        let x_in_multiplier = builder.mul_add(zeta_minus_one, scalar_bit_1, one);
        let x2 = builder.mul(x_in_multiplier, x_in);
        let y_in_multiplier = builder.mul_sub(scalar_bit_0, two, one);
        let y2 = builder.mul(y_in_multiplier, y_in);

        let x1_minus_x2 = builder.sub(x1, x2);
        let x1_plus_x2 = builder.add(x1, x2);
        let x1_minus_x3 = builder.sub(x1, x3);
        let y1_minus_y2 = builder.sub(y1, y2);

        let lambda = builder.mul(y1_minus_y2, inverse);
        let computed_x3 = builder.mul_sub(lambda, lambda, x1_plus_x2);
        let computed_y3 = builder.mul_sub(lambda, x1_minus_x3, y1);

        let unsigned_limb = builder.mul_add(scalar_bit_1, two, scalar_bit_0);
        let computed_scalar_acc_unsigned_new = builder.mul_add(scalar_acc_unsigned_old, four, unsigned_limb);

        // This is based on Algorithm 2 in the Halo paper.
        // TODO: Wrong zeta used here?
        let signed_limb_multiplier = builder.mul_add(zeta_minus_one, scalar_bit_1, one);
        let signed_limb_sign = builder.mul_sub(scalar_bit_0, two, one);
        let signed_limb = builder.mul(signed_limb_sign, signed_limb_multiplier);
        let computed_scalar_acc_signed_new = builder.mul_add(scalar_acc_signed_old, two, signed_limb);

        vec![
            builder.sub(computed_x3, x3),
            builder.sub(computed_y3, y3),
            builder.sub(computed_scalar_acc_unsigned_new, scalar_acc_unsigned_new),
            builder.sub(computed_scalar_acc_signed_new, scalar_acc_signed_new),
            assert_binary_recursively(builder, scalar_bit_0),
            assert_binary_recursively(builder, scalar_bit_1),
            assert_inverses_recursively(builder, inverse, x1_minus_x2),
        ]
    }
}

impl<C: HaloCurve, InnerC: HaloCurve<BaseField=C::ScalarField>>
WitnessGenerator<C::ScalarField> for CurveEndoGate<C, InnerC> {
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

    fn generate(&self, _constants: &Vec<Vec<C::ScalarField>>, witness: &PartialWitness<InnerC::BaseField>) -> PartialWitness<InnerC::BaseField> {
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
        let group_acc_old = AffinePoint::<InnerC>::nonzero(group_acc_old_x, group_acc_old_y);

        let scalar_acc_unsigned_old = witness.get_wire(scalar_acc_unsigned_old_target);
        let scalar_acc_signed_old = witness.get_wire(scalar_acc_signed_old_target);

        let scalar_bit_0 = witness.get_wire(scalar_bit_0_target);
        let scalar_bit_1 = witness.get_wire(scalar_bit_1_target);
        debug_assert!(scalar_bit_0.is_zero() || scalar_bit_0.is_one());
        debug_assert!(scalar_bit_1.is_zero() || scalar_bit_1.is_one());

        let p_x = witness.get_wire(addend_x_target);
        let p_y = witness.get_wire(addend_y_target);

        let mut s_i_x = p_x;
        if scalar_bit_0 == InnerC::BaseField::ONE {
            s_i_x = s_i_x * InnerC::ZETA;
        }
        let mut s_i_y = p_y;
        if scalar_bit_1 == InnerC::BaseField::ZERO {
            s_i_y = -s_i_y;
        }
        let s_i = AffinePoint::nonzero(s_i_x, s_i_y);
        let group_acc_new = group_acc_old + s_i;

        let scalar_acc_unsigned_new = scalar_acc_unsigned_old.quadruple()
            + scalar_bit_0 + scalar_bit_1.double();

        // This is based on Algorithm 2 in the Halo paper.
        let mut scalar_acc_signed_limb = if scalar_bit_0 == InnerC::BaseField::ONE {
            InnerC::BaseField::ONE
        } else {
            InnerC::BaseField::NEG_ONE
        };
        if scalar_bit_1 == InnerC::BaseField::ONE {
            scalar_acc_signed_limb = scalar_acc_signed_limb * InnerC::ZETA;
        }
        let scalar_acc_signed_new = scalar_acc_signed_old.double() + scalar_acc_signed_limb;

        // Here's where our abstraction leaks a bit. Although we already have the sum, we need to
        // redo part of the computation in order to populate the purported inverse wire.
        let dx = group_acc_old_x - p_x;
        let _dy = group_acc_old_y - p_y;
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
pub struct RescueStepAGate<C: HaloCurve> {
    pub index: usize,
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> RescueStepAGate<C> {
    pub fn new(index: usize) -> Self {
        RescueStepAGate { index, _phantom: PhantomData }
    }

    /// Returns the index of the `i`th accumulator wire.
    pub fn wire_acc(i: usize) -> usize {
        return i;
    }

    /// Returns the index of the `i`th root wire.
    pub fn wire_root(i: usize) -> usize {
        return RESCUE_SPONGE_WIDTH + i;
    }
}

impl<C: HaloCurve> Gate<C> for RescueStepAGate<C> {
    const NAME: &'static str = "RescueStepAGate";

    const PREFIX: &'static [bool] = &[false, false];

    fn evaluate_unfiltered(
        local_constant_values: &[C::ScalarField],
        local_wire_values: &[C::ScalarField],
        right_wire_values: &[C::ScalarField],
        _below_wire_values: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        let ins: Vec<C::ScalarField> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| local_wire_values[Self::wire_acc(i)])
            .collect();
        let outs: Vec<C::ScalarField> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| right_wire_values[Self::wire_acc(i)])
            .collect();
        let roots: Vec<C::ScalarField> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| local_wire_values[Self::wire_root(i)])
            .collect();

        let mds = mds_matrix::<C::ScalarField>(RESCUE_SPONGE_WIDTH);

        let mut constraints = Vec::new();
        for i in 0..RESCUE_SPONGE_WIDTH {
            constraints.push(roots[i].exp_usize(5) - ins[i]);

            let mut computed_out_i = local_constant_values[Self::PREFIX.len() + i];
            for j in 0..RESCUE_SPONGE_WIDTH {
                computed_out_i = computed_out_i + mds.get(i, j) * roots[j];
            }
            constraints.push(computed_out_i - outs[i]);
        }
        constraints
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        _below_wire_values: &[Target],
    ) -> Vec<Target> {
        let ins: Vec<Target> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| local_wire_values[Self::wire_acc(i)])
            .collect();

        let outs: Vec<Target> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| right_wire_values[Self::wire_acc(i)])
            .collect();

        let roots: Vec<Target> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| local_wire_values[Self::wire_root(i)])
            .collect();

        let mds = mds_matrix::<C::ScalarField>(RESCUE_SPONGE_WIDTH);

        let mut constraints = Vec::new();
        for i in 0..RESCUE_SPONGE_WIDTH {
            let computed_in_i = builder.exp_constant_usize(roots[i], 5);
            constraints.push(builder.sub(computed_in_i, ins[i]));

            let mut computed_out_i = local_constant_values[Self::PREFIX.len() + i];
            for j in 0..RESCUE_SPONGE_WIDTH {
                let mds_entry = builder.constant_wire(mds.get(i, j));
                computed_out_i = builder.mul_add(mds_entry, roots[j], computed_out_i);
            }
            constraints.push(builder.sub(computed_out_i, outs[i]));
        }
        constraints
    }
}

impl<C: HaloCurve> WitnessGenerator<C::ScalarField> for RescueStepAGate<C> {
    fn dependencies(&self) -> Vec<Target> {
        (0..RESCUE_SPONGE_WIDTH)
            .map(|i| Target::Wire(Wire { gate: self.index, input: Self::wire_acc(i) }))
            .collect()
    }

    fn generate(&self, constants: &Vec<Vec<C::ScalarField>>, witness: &PartialWitness<C::ScalarField>) -> PartialWitness<C::ScalarField> {
        let constants = &constants[self.index];

        let ins: Vec<C::ScalarField> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| witness.get_wire(Wire { gate: self.index, input: Self::wire_acc(i) }))
            .collect();

        let roots: Vec<C::ScalarField> = ins.iter()
            .map(|n| n.kth_root_u32(5))
            .collect();

        let mds = mds_matrix::<C::ScalarField>(RESCUE_SPONGE_WIDTH);

        let mut result = PartialWitness::new();
        for i in 0..RESCUE_SPONGE_WIDTH {
            let wire_root_i = Wire { gate: self.index, input: Self::wire_root(i) };
            result.set_wire(wire_root_i, roots[i]);

            let mut out_i = constants[Self::PREFIX.len() + i];
            for j in 0..RESCUE_SPONGE_WIDTH {
                out_i = out_i + mds.get(i, j) * roots[j];
            }
            let wire_out_i = Wire { gate: self.index + 1, input: Self::wire_acc(i) };
            result.set_wire(wire_out_i, out_i);
        }
        result
    }
}

/// The second step of Rescue, i.e. the one with the `x^5` layer.
pub struct RescueStepBGate<C: HaloCurve> {
    pub index: usize,
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> RescueStepBGate<C> {
    pub fn new(index: usize) -> Self {
        RescueStepBGate { index, _phantom: PhantomData }
    }

    /// Returns the index of the `i`th accumulator wire.
    pub fn wire_acc(i: usize) -> usize {
        return i;
    }
}

impl<C: HaloCurve> Gate<C> for RescueStepBGate<C> {
    const NAME: &'static str = "RescueStepBGate";

    const PREFIX: &'static [bool] = &[false, true];

    fn evaluate_unfiltered(
        local_constant_values: &[C::ScalarField],
        local_wire_values: &[C::ScalarField],
        right_wire_values: &[C::ScalarField],
        _below_wire_values: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        let ins: Vec<C::ScalarField> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| local_wire_values[Self::wire_acc(i)])
            .collect();

        let exps: Vec<C::ScalarField> = ins.into_iter()
            .map(|n| n.exp_usize(5))
            .collect();

        let outs: Vec<C::ScalarField> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| right_wire_values[Self::wire_acc(i)])
            .collect();

        let mds = mds_matrix::<C::ScalarField>(RESCUE_SPONGE_WIDTH);

        let mut constraints = Vec::new();
        for i in 0..RESCUE_SPONGE_WIDTH {
            let mut computed_out_i = local_constant_values[Self::PREFIX.len() + i];
            for j in 0..RESCUE_SPONGE_WIDTH {
                computed_out_i = computed_out_i + mds.get(i, j) * exps[i];
            }
            constraints.push(computed_out_i - outs[i]);
        }
        constraints
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        _below_wire_values: &[Target],
    ) -> Vec<Target> {
        let ins: Vec<Target> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| local_wire_values[Self::wire_acc(i)])
            .collect();

        let exps: Vec<Target> = ins.into_iter()
            .map(|n| builder.exp_constant_usize(n, 5))
            .collect();

        let outs: Vec<Target> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| right_wire_values[Self::wire_acc(i)])
            .collect();

        let mds = mds_matrix::<C::ScalarField>(RESCUE_SPONGE_WIDTH);

        let mut constraints = Vec::new();
        for i in 0..RESCUE_SPONGE_WIDTH {
            let mut computed_out_i = local_constant_values[Self::PREFIX.len() + i];
            for j in 0..RESCUE_SPONGE_WIDTH {
                let mds_entry = builder.constant_wire(mds.get(i, j));
                computed_out_i = builder.mul_add(mds_entry, exps[i], computed_out_i);
            }
            constraints.push(builder.sub(computed_out_i, outs[i]));
        }
        constraints
    }
}

impl<C: HaloCurve> WitnessGenerator<C::ScalarField> for RescueStepBGate<C> {
    fn dependencies(&self) -> Vec<Target> {
        (0..RESCUE_SPONGE_WIDTH)
            .map(|i| Target::Wire(Wire { gate: self.index, input: Self::wire_acc(i) }))
            .collect()
    }

    fn generate(&self, constants: &Vec<Vec<C::ScalarField>>, witness: &PartialWitness<C::ScalarField>) -> PartialWitness<C::ScalarField> {
        let constants = &constants[self.index];

        let ins: Vec<C::ScalarField> = (0..RESCUE_SPONGE_WIDTH)
            .map(|i| witness.get_wire(Wire { gate: self.index, input: Self::wire_acc(i) }))
            .collect();

        let exps: Vec<C::ScalarField> = ins.iter()
            .map(|n| n.exp_usize(5))
            .collect();

        let mds = mds_matrix::<C::ScalarField>(RESCUE_SPONGE_WIDTH);

        let mut result = PartialWitness::new();
        for i in 0..RESCUE_SPONGE_WIDTH {
            let mut out_i = constants[Self::PREFIX.len() + i];
            for j in 0..RESCUE_SPONGE_WIDTH {
                out_i = out_i + mds.get(i, j) * exps[j];
            }
            let wire_out_i = Wire { gate: self.index + 1, input: Self::wire_acc(i) };
            result.set_wire(wire_out_i, out_i);
        }
        result
    }
}

/// A gate for accumulating base-4 limbs.
pub struct Base4SumGate<C: Curve> {
    pub index: usize,
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> Base4SumGate<C> {
    pub fn new(index: usize) -> Self {
        Base4SumGate { index, _phantom: PhantomData }
    }

    pub const WIRE_ACC_OLD: usize = 0;
    pub const WIRE_ACC_NEW: usize = 1;
    pub const WIRE_LIMB_0: usize = 2;
    pub const NUM_LIMBS: usize = 7;
}

impl<C: HaloCurve> Gate<C> for Base4SumGate<C> {
    const NAME: &'static str = "Base4SumGate";

    const PREFIX: &'static [bool] = &[true, false, false, false];

    fn evaluate_unfiltered(
        _local_constant_values: &[C::ScalarField],
        local_wire_values: &[C::ScalarField],
        _right_wire_values: &[C::ScalarField],
        _below_wire_values: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        let acc_old = local_wire_values[Self::WIRE_ACC_OLD];
        let acc_new = local_wire_values[Self::WIRE_ACC_NEW];
        let limbs: Vec<C::ScalarField> = (0..Self::NUM_LIMBS).map(
            |i| local_wire_values[Self::WIRE_LIMB_0 + i]
        ).collect();

        let mut computed_acc_new = acc_old;
        for &limb in &limbs {
            computed_acc_new = computed_acc_new.quadruple() + limb;
        }

        let mut constraints = vec![computed_acc_new - acc_new];
        for limb in limbs {
            let mut product = C::ScalarField::ONE;
            for j in 0..4 {
                product = product * (limb - C::ScalarField::from_canonical_usize(j));
            }
            constraints.push(product);
        }
        constraints
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        _local_constant_values: &[Target],
        local_wire_values: &[Target],
        _right_wire_values: &[Target],
        _below_wire_values: &[Target],
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

impl<C: HaloCurve> WitnessGenerator<C::ScalarField> for Base4SumGate<C> {
    fn dependencies(&self) -> Vec<Target> {
        Vec::new()
    }

    fn generate(&self, _constants: &Vec<Vec<C::ScalarField>>, _witness: &PartialWitness<C::ScalarField>) -> PartialWitness<C::ScalarField> {
        // For base 4 decompositions, we don't do any witness generation on a per-gate level.
        // Instead, we have a single generator which generates values for an entire decomposition.
        PartialWitness::new()
    }
}

/// A gate which can be configured to perform various arithmetic. In particular, it computes
///
/// ```text
/// output := const_0 * multiplicand_0 * multiplicand_1 + const_1 * addend
/// ```
pub struct ArithmeticGate<C: HaloCurve> {
    pub index: usize,
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> ArithmeticGate<C> {
    pub fn new(index: usize) -> Self {
        ArithmeticGate { index, _phantom: PhantomData }
    }

    pub const WIRE_MULTIPLICAND_0: usize = 0;
    pub const WIRE_MULTIPLICAND_1: usize = 1;
    pub const WIRE_ADDEND: usize = 2;
    pub const WIRE_OUTPUT: usize = 3;
}

impl<C: HaloCurve> Gate<C> for ArithmeticGate<C> {
    const NAME: &'static str = "ArithmeticGate";

    const PREFIX: &'static [bool] = &[true, false, false, true];

    fn evaluate_unfiltered(
        local_constant_values: &[C::ScalarField],
        local_wire_values: &[C::ScalarField],
        _right_wire_values: &[C::ScalarField],
        _below_wire_values: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        let const_0 = local_constant_values[Self::PREFIX.len()];
        let const_1 = local_constant_values[Self::PREFIX.len() + 1];
        let multiplicand_0 = local_wire_values[Self::WIRE_MULTIPLICAND_0];
        let multiplicand_1 = local_wire_values[Self::WIRE_MULTIPLICAND_1];
        let addend = local_wire_values[Self::WIRE_ADDEND];
        let output = local_wire_values[Self::WIRE_OUTPUT];
        let computed_output = const_0 * multiplicand_0 * multiplicand_1 + const_1 * addend;
        vec![computed_output - output]
    }

    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        local_constant_values: &[Target],
        local_wire_values: &[Target],
        _right_wire_values: &[Target],
        _below_wire_values: &[Target],
    ) -> Vec<Target> {
        let const_0 = local_constant_values[Self::PREFIX.len()];
        let const_1 = local_constant_values[Self::PREFIX.len() + 1];
        let multiplicand_0 = local_wire_values[Self::WIRE_MULTIPLICAND_0];
        let multiplicand_1 = local_wire_values[Self::WIRE_MULTIPLICAND_1];
        let addend = local_wire_values[Self::WIRE_ADDEND];
        let output = local_wire_values[Self::WIRE_OUTPUT];

        let product_term = builder.mul_many(&[const_0, multiplicand_0, multiplicand_1]);
        let addend_term = builder.mul(const_1, addend);
        let computed_output = builder.add_many(&[product_term, addend_term]);
        vec![builder.sub(computed_output, output)]
    }
}

impl<C: HaloCurve> WitnessGenerator<C::ScalarField> for ArithmeticGate<C> {
    fn dependencies(&self) -> Vec<Target> {
        vec![
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_MULTIPLICAND_0 }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_MULTIPLICAND_1 }),
            Target::Wire(Wire { gate: self.index, input: Self::WIRE_ADDEND }),
        ]
    }

    fn generate(&self, constants: &Vec<Vec<C::ScalarField>>, witness: &PartialWitness<C::ScalarField>) -> PartialWitness<C::ScalarField> {
        let multiplicand_0_target = Wire { gate: self.index, input: Self::WIRE_MULTIPLICAND_0 };
        let multiplicand_1_target = Wire { gate: self.index, input: Self::WIRE_MULTIPLICAND_1 };
        let addend_target = Wire { gate: self.index, input: Self::WIRE_ADDEND };
        let output_target = Wire { gate: self.index, input: Self::WIRE_OUTPUT };

        let const_0 = constants[self.index][Self::PREFIX.len()];
        let const_1 = constants[self.index][Self::PREFIX.len() + 1];

        let multiplicand_0 = witness.get_wire(multiplicand_0_target);
        let multiplicand_1 = witness.get_wire(multiplicand_1_target);
        let addend = witness.get_wire(addend_target);

        let output = const_0 * multiplicand_0 * multiplicand_1 + const_1 * addend;

        let mut result = PartialWitness::new();
        result.set_wire(output_target, output);
        result
    }
}

/// Test that a gate's constraints are within degree 8n, including the gate prefix filter.
#[macro_export]
macro_rules! test_gate_low_degree {
    ($method:ident, $curve:ty, $gate:ty) => {
        #[test]
        #[ignore] // Too slow to run regularly.
        #[allow(non_snake_case)]
        fn $method() {
            type C = $curve;
            type SF = <C as $crate::curve::Curve>::ScalarField;

            let n = 256;
            let fft_precomputation_n = $crate::fft::fft_precompute::<SF>(n);
            let fft_precomputation_8n = $crate::fft::fft_precompute::<SF>(8 * n);
            let fft_precomputation_16n = $crate::fft::fft_precompute::<SF>(16 * n);

            // Generate random constant and wire polynomials.
            let mut constant_values_n: Vec<Vec<SF>> = vec![Vec::new(); $crate::plonk::NUM_CONSTANTS];
            let mut wire_values_n: Vec<Vec<SF>> = vec![Vec::new(); $crate::plonk::NUM_WIRES];
            for i in 0..n {
                for points in constant_values_n.iter_mut() {
                    points.push(<SF as $crate::field::Field>::rand())
                }
                for points in wire_values_n.iter_mut() {
                    points.push(<SF as $crate::field::Field>::rand())
                }
            }

            // Low-degree extend them to 16n values.
            let mut constant_coeffs_16n = $crate::plonk_util::values_to_coeffs(&constant_values_n, &fft_precomputation_n);
            let mut wire_coeffs_16n = $crate::plonk_util::values_to_coeffs(&wire_values_n, &fft_precomputation_n);
            for coeffs in constant_coeffs_16n.iter_mut().chain(wire_coeffs_16n.iter_mut()) {
                while coeffs.len() < 16 * n {
                    coeffs.push(<SF as $crate::field::Field>::ZERO);
                }
            }
            let constant_values_16n: Vec<Vec<SF>> = constant_coeffs_16n.iter()
                .map(|coeffs| $crate::fft::fft_with_precomputation_power_of_2(coeffs, &fft_precomputation_16n))
                .collect();
            let wire_values_16n: Vec<Vec<SF>> = wire_coeffs_16n.iter()
                .map(|coeffs| $crate::fft::fft_with_precomputation_power_of_2(coeffs, &fft_precomputation_16n))
                .collect();

            // Make sure each extended polynomial is still degree n.
            for values_16n in constant_values_16n.iter().chain(wire_values_16n.iter()) {
                assert!($crate::plonk_util::polynomial_degree(values_16n, &fft_precomputation_16n) < n);
            }

            let constant_values_16n_t = $crate::util::transpose(&constant_values_16n);
            let wire_values_16n_t = $crate::util::transpose(&wire_values_16n);

            // Evaluate constraints at each of our 16n points.
            let mut constraint_values_16n: Vec<Vec<SF>> = Vec::new();
            for i in 0..16 * n {
                let constraints: Vec<SF> = <$gate as $crate::plonk_gates::Gate<C>>::evaluate_filtered(
                    &constant_values_16n_t[i],
                    &wire_values_16n_t[i],
                    &wire_values_16n_t[(i + 16) % (16 * n)],
                    &wire_values_16n_t[(i + 16 * $crate::plonk::GRID_WIDTH) % (16 * n)]);
                for (j, &c) in constraints.iter().enumerate() {
                    if constraint_values_16n.len() <= j {
                        constraint_values_16n.push(Vec::new());
                    }
                    constraint_values_16n[j].push(c);
                }
            }

            // Check that the degree of each constraint is within the limit.
            let constraint_degrees = constraint_values_16n.iter()
                .map(|c| $crate::plonk_util::polynomial_degree(c, &fft_precomputation_16n))
                .collect::<Vec<_>>();
            let max_degree_excl = (crate::plonk::QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER + 1) * n;
            for (i, &deg) in constraint_degrees.iter().enumerate() {
                assert!(deg < max_degree_excl,
                "Constraint at index {} has degree {}; should be less than {}n = {}",
                i,
                deg,
                crate::plonk::QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER + 1,
                max_degree_excl);
            }
        }
    };
}

#[cfg(test)]
mod tests {
    use crate::{Tweedledum, Tweedledee, PublicInputGate, CurveAddGate, CurveDblGate, CurveEndoGate, Base4SumGate, BufferGate, ConstantGate, ArithmeticGate, RescueStepAGate, RescueStepBGate};

    test_gate_low_degree!(low_degree_PublicInputGate, Tweedledum, PublicInputGate<Tweedledum>);
    test_gate_low_degree!(low_degree_CurveAddGate, Tweedledum, CurveAddGate<Tweedledum, Tweedledee>);
    test_gate_low_degree!(low_degree_CurveDblGate, Tweedledum, CurveDblGate<Tweedledum, Tweedledee>);
    test_gate_low_degree!(low_degree_CurveEndoGate, Tweedledum, CurveEndoGate<Tweedledum, Tweedledee>);
    test_gate_low_degree!(low_degree_Base4SumGate, Tweedledum, Base4SumGate<Tweedledum>);
    test_gate_low_degree!(low_degree_BufferGate, Tweedledum, BufferGate<Tweedledum>);
    test_gate_low_degree!(low_degree_ConstantGate, Tweedledum, ConstantGate<Tweedledum>);
    test_gate_low_degree!(low_degree_ArithmeticGate, Tweedledum, ArithmeticGate<Tweedledum>);
    test_gate_low_degree!(low_degree_RescueStepAGate, Tweedledum, RescueStepAGate<Tweedledum>);
    test_gate_low_degree!(low_degree_RescueStepBGate, Tweedledum, RescueStepBGate<Tweedledum>);
}
