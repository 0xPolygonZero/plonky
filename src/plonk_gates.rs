use std::marker::PhantomData;

use crate::{AffinePoint, Circuit, CircuitBuilder, Curve, Field, Wire, GRID_WIDTH, HaloEndomorphismCurve, PartialWitness, Target, WitnessGenerator};
use crate::mds::mds;

pub trait Gate<F: Field>: WitnessGenerator<F> {
    /// In order to combine the constraints of various gate types into a unified constraint set, we
    /// assign each gate type a binary prefix such that no two prefixes overlap.
    const PREFIX: &'static [bool];

    /// Evaluate the constraints implied by this gate at the given challenge point.
    ///
    /// For example, if the gate computes `c = a * b`, this should return `[c(x) - a(x) * b(x)]`,
    /// where `x` is the challenge point.
    fn evaluate(&self,
                local_constant_values: Vec<F>,
                local_wire_values: Vec<F>,
                right_wire_values: Vec<F>,
                below_wire_values: Vec<F>) -> Vec<F>;

    /// Like the other `evaluate` method, but in the context of a recursive circuit.
    fn evaluate_recursively(&self,
                            builder: &mut CircuitBuilder<F>,
                            local_constant_values: Vec<Target>,
                            local_wire_values: Vec<Target>,
                            right_wire_values: Vec<Target>,
                            below_wire_values: Vec<Target>) -> Vec<Target>;
}

/// A gate which doesn't perform any arithmetic, but just acts as a buffer for receiving data. This
/// is used in a couple ways:
/// * Public inputs can be "received" via `WIRE_BUFFER_PI`.
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
    pub const WIRE_BUFFER_PI: usize = 4;
    pub const WIRE_BUFFER_CONST: usize = 5;
}

impl<F: Field> Gate<F> for BufferGate {
    const PREFIX: &'static [bool] = &[false, false, true, true];

    fn evaluate(
        &self,
        local_constant_values: Vec<F>,
        local_wire_values: Vec<F>,
        right_wire_values: Vec<F>,
        below_wire_values: Vec<F>,
    ) -> Vec<F> {
        unimplemented!()
    }

    fn evaluate_recursively(
        &self,
        builder: &mut CircuitBuilder<F>,
        local_constant_values: Vec<Target>,
        local_wire_values: Vec<Target>,
        right_wire_values: Vec<Target>,
        below_wire_values: Vec<Target>,
    ) -> Vec<Target> {
        unimplemented!()
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
    const PREFIX: &'static [bool] = &[false, false, false, false, true];

    fn evaluate(
        &self,
        local_constant_values: Vec<C::BaseField>,
        local_wire_values: Vec<C::BaseField>,
        right_wire_values: Vec<C::BaseField>,
        below_wire_values: Vec<C::BaseField>,
    ) -> Vec<C::BaseField> {
        let local_group_acc_x = local_wire_values[Self::WIRE_GROUP_ACC_X];
        let local_group_acc_y = local_wire_values[Self::WIRE_GROUP_ACC_Y];
        let right_group_acc_x = right_wire_values[Self::WIRE_GROUP_ACC_X];
        let right_group_acc_y = right_wire_values[Self::WIRE_GROUP_ACC_Y];

        let scalar_acc_old = local_wire_values[Self::WIRE_SCALAR_ACC_OLD];
        let scalar_acc_new = local_wire_values[Self::WIRE_SCALAR_ACC_NEW];
        let addend_x = local_wire_values[Self::WIRE_ADDEND_X];
        let addend_y = local_wire_values[Self::WIRE_ADDEND_Y];
        let scalar_bit = local_wire_values[Self::WIRE_SCALAR_BIT];
        let inverse = local_wire_values[Self::WIRE_INVERSE];

        let dx = local_group_acc_x - addend_x;
        let dy = local_group_acc_y - addend_y;
        let lambda = dy * inverse;
        let sum_x = lambda.square() - local_group_acc_x - addend_x;
        let sum_y = lambda * dx - local_group_acc_y;

        vec![
            sum_x - right_group_acc_x,
            sum_y - right_group_acc_y,
            scalar_acc_new - scalar_acc_old.double() + scalar_bit,
            scalar_bit * (scalar_bit - C::BaseField::ONE),
            inverse * dx - C::BaseField::ONE,
        ]
    }

    fn evaluate_recursively(
        &self,
        builder: &mut CircuitBuilder<C::BaseField>,
        local_constant_values: Vec<Target>,
        local_wire_values: Vec<Target>,
        right_wire_values: Vec<Target>,
        below_wire_values: Vec<Target>,
    ) -> Vec<Target> {
        unimplemented!()
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
    const PREFIX: &'static [bool] = &[false, false, false, true, false];

    fn evaluate(
        &self,
        local_constant_values: Vec<C::BaseField>,
        local_wire_values: Vec<C::BaseField>,
        right_wire_values: Vec<C::BaseField>,
        below_wire_values: Vec<C::BaseField>,
    ) -> Vec<C::BaseField> {
        unimplemented!()
    }

    fn evaluate_recursively(
        &self,
        builder: &mut CircuitBuilder<C::BaseField>,
        local_constant_values: Vec<Target>,
        local_wire_values: Vec<Target>,
        right_wire_values: Vec<Target>,
        below_wire_values: Vec<Target>,
    ) -> Vec<Target> {
        unimplemented!()
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
    const PREFIX: &'static [bool] = &[false, false, false, true, true];

    fn evaluate(
        &self,
        local_constant_values: Vec<C::BaseField>,
        local_wire_values: Vec<C::BaseField>,
        right_wire_values: Vec<C::BaseField>,
        below_wire_values: Vec<C::BaseField>,
    ) -> Vec<C::BaseField> {
        unimplemented!()
    }

    fn evaluate_recursively(
        &self,
        builder: &mut CircuitBuilder<C::BaseField>,
        local_constant_values: Vec<Target>,
        local_wire_values: Vec<Target>,
        right_wire_values: Vec<Target>,
        below_wire_values: Vec<Target>,
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
    const PREFIX: &'static [bool] = &[true, false];

    fn evaluate(
        &self,
        local_constant_values: Vec<F>,
        local_wire_values: Vec<F>,
        right_wire_values: Vec<F>,
        below_wire_values: Vec<F>,
    ) -> Vec<F> {
        unimplemented!()
    }

    fn evaluate_recursively(
        &self,
        builder: &mut CircuitBuilder<F>,
        local_constant_values: Vec<Target>,
        local_wire_values: Vec<Target>,
        right_wire_values: Vec<Target>,
        below_wire_values: Vec<Target>,
    ) -> Vec<Target> {
        unimplemented!()
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

        let out_0 = mds::<F>(3, 0, 0) * root_0 + mds::<F>(3, 0, 1) * root_1 + mds::<F>(3, 0, 2) * root_2 + constants[Self::PREFIX.len() + 0];
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
    const PREFIX: &'static [bool] = &[true, true];

    fn evaluate(
        &self,
        local_constant_values: Vec<F>,
        local_wire_values: Vec<F>,
        right_wire_values: Vec<F>,
        below_wire_values: Vec<F>,
    ) -> Vec<F> {
        unimplemented!()
    }

    fn evaluate_recursively(
        &self,
        builder: &mut CircuitBuilder<F>,
        local_constant_values: Vec<Target>,
        local_wire_values: Vec<Target>,
        right_wire_values: Vec<Target>,
        below_wire_values: Vec<Target>,
    ) -> Vec<Target> {
        unimplemented!()
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

        let out_0 = mds::<F>(3, 0, 0) * exp_0 + mds::<F>(3, 0, 1) * exp_1 + mds::<F>(3, 0, 2) * exp_2 + constants[Self::PREFIX.len() + 0];
        let out_1 = mds::<F>(3, 1, 0) * exp_0 + mds::<F>(3, 1, 1) * exp_1 + mds::<F>(3, 1, 2) * exp_2 + constants[Self::PREFIX.len() + 1];
        let out_2 = mds::<F>(3, 2, 0) * exp_0 + mds::<F>(3, 2, 1) * exp_1 + mds::<F>(3, 2, 2) * exp_2 + constants[Self::PREFIX.len() + 2];

        let mut result = PartialWitness::new();
        result.set_wire(out_0_target, out_0);
        result.set_wire(out_1_target, out_1);
        result.set_wire(out_2_target, out_2);
        result
    }
}

pub(crate) struct Base4SumGate {
    pub index: usize,
    /// Make the constructor private.
    _private: (),
}

impl Base4SumGate {
    pub fn new(index: usize) -> Self {
        Base4SumGate { index, _private: () }
    }

    pub const WIRE_ACC: usize = 0;
    pub const WIRE_LIMB_0: usize = 1;
    pub const WIRE_LIMB_1: usize = 2;
    pub const WIRE_LIMB_2: usize = 3;
    pub const WIRE_LIMB_3: usize = 4;
    pub const WIRE_LIMB_4: usize = 5;
    pub const WIRE_LIMB_5: usize = 6;
    pub const WIRE_LIMB_6: usize = 7;
    pub const WIRE_LIMB_7: usize = 8;
}

impl<F: Field> Gate<F> for Base4SumGate {
    const PREFIX: &'static [bool] = &[false, false, true, false, false];

    fn evaluate(
        &self,
        local_constant_values: Vec<F>,
        local_wire_values: Vec<F>,
        right_wire_values: Vec<F>,
        below_wire_values: Vec<F>,
    ) -> Vec<F> {
        unimplemented!()
    }

    fn evaluate_recursively(
        &self,
        builder: &mut CircuitBuilder<F>,
        local_constant_values: Vec<Target>,
        local_wire_values: Vec<Target>,
        right_wire_values: Vec<Target>,
        below_wire_values: Vec<Target>,
    ) -> Vec<Target> {
        unimplemented!()
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

/// A "multiply and add" gate. In particular, it computes
///
/// ```text
/// output := const_0 * multiplicand_0 * multiplicand_1 + const_1 * addend + const_2
/// ```
pub(crate) struct MaddGate<F: Field> {
    pub index: usize,
    _phantom: PhantomData<F>,
}

impl<F: Field> MaddGate<F> {
    pub fn new(index: usize) -> Self {
        MaddGate { index, _phantom: PhantomData }
    }

    pub const WIRE_MULTIPLICAND_0: usize = 0;
    pub const WIRE_MULTIPLICAND_1: usize = 1;
    pub const WIRE_ADDEND: usize = 2;
    pub const WIRE_OUTPUT: usize = 3;
}

impl<F: Field> Gate<F> for MaddGate<F> {
    const PREFIX: &'static [bool] = &[false, true];

    fn evaluate(
        &self,
        local_constant_values: Vec<F>,
        local_wire_values: Vec<F>,
        right_wire_values: Vec<F>,
        below_wire_values: Vec<F>,
    ) -> Vec<F> {
        unimplemented!()
    }

    fn evaluate_recursively(
        &self,
        builder: &mut CircuitBuilder<F>,
        local_constant_values: Vec<Target>,
        local_wire_values: Vec<Target>,
        right_wire_values: Vec<Target>,
        below_wire_values: Vec<Target>,
    ) -> Vec<Target> {
        unimplemented!()
    }
}

impl<F: Field> WitnessGenerator<F> for MaddGate<F> {
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

        let const_0 = circuit.gate_constants[self.index][Self::PREFIX.len() + 0];
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
