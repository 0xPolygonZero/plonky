use crate::{GateInput, WitnessGenerator, Field, PartialWitness, AffinePoint, HaloEndomorphismCurve, Curve};
use std::marker::PhantomData;

pub(crate) trait Gate<F: Field>: 'static + WitnessGenerator<F> {
    const ID: usize;
}

pub(crate) struct NoopGate { index: usize }

impl NoopGate {
    const WIRE_BUFFER_0: usize = 0;
    const WIRE_BUFFER_1: usize = 1;
    const WIRE_BUFFER_2: usize = 2;
    const WIRE_BUFFER_3: usize = 3;
}

impl<F: Field> Gate<F> for NoopGate {
    const ID: usize = 0;
}

impl<F: Field> WitnessGenerator<F> for NoopGate {
    fn dependencies(&self) -> Vec<GateInput> {
        Vec::new()
    }

    fn generate(&self, _witness: &PartialWitness<F>) -> PartialWitness<F> {
        // This gate does not generate any witness values.
        PartialWitness::new()
    }
}

/// A gate which performs incomplete point addition, conditioned on an input bit. In order to
/// facilitate MSMs which use this gate, it also adds the bit to an accumulator. `C` is the curve
/// of the inner proof.
pub(crate) struct CurveAddGate<C: Curve> {
    index: usize,
    _phantom: PhantomData<C>,
}

impl<C: Curve> CurveAddGate<C> {
    const WIRE_GROUP_ACC_X: usize = 0;
    const WIRE_GROUP_ACC_Y: usize = 1;
    const WIRE_SCALAR_ACC_OLD: usize = 2;
    const WIRE_SCALAR_ACC_NEW: usize = 3;
    const WIRE_ADDEND_X: usize = 4;
    const WIRE_ADDEND_Y: usize = 5;
    const WIRE_SCALAR_BIT: usize = 6;
    const WIRE_INVERSE: usize = 7;
}

impl<C: Curve> Gate<C::BaseField> for CurveAddGate<C> {
    const ID: usize = 1;
}

impl<C: Curve> WitnessGenerator<C::BaseField> for CurveAddGate<C> {
    fn dependencies(&self) -> Vec<GateInput> {
        vec![
            GateInput { gate: self.index, input: Self::WIRE_GROUP_ACC_X },
            GateInput { gate: self.index, input: Self::WIRE_GROUP_ACC_Y },
            GateInput { gate: self.index, input: Self::WIRE_SCALAR_ACC_OLD },
            GateInput { gate: self.index, input: Self::WIRE_ADDEND_X },
            GateInput { gate: self.index, input: Self::WIRE_ADDEND_Y },
            GateInput { gate: self.index, input: Self::WIRE_SCALAR_BIT },
        ]
    }

    fn generate(&self, witness: &PartialWitness<<C as Curve>::BaseField>) -> PartialWitness<<C as Curve>::BaseField> {
        let group_acc_old_x_target = GateInput { gate: self.index, input: Self::WIRE_GROUP_ACC_X };
        let group_acc_new_x_target = GateInput { gate: self.index + 1, input: Self::WIRE_GROUP_ACC_X };
        let group_acc_old_y_target = GateInput { gate: self.index, input: Self::WIRE_GROUP_ACC_Y };
        let group_acc_new_y_target = GateInput { gate: self.index + 1, input: Self::WIRE_GROUP_ACC_Y };
        let scalar_acc_old_target = GateInput { gate: self.index, input: Self::WIRE_SCALAR_ACC_OLD };
        let scalar_acc_new_target = GateInput { gate: self.index, input: Self::WIRE_SCALAR_ACC_NEW };
        let addend_x_target = GateInput { gate: self.index, input: Self::WIRE_ADDEND_X };
        let addend_y_target = GateInput { gate: self.index, input: Self::WIRE_ADDEND_Y };
        let scalar_bit_target = GateInput { gate: self.index, input: Self::WIRE_SCALAR_BIT };
        let inverse_target = GateInput { gate: self.index, input: Self::WIRE_INVERSE };

        let group_acc_old_x = witness.wire_values[&group_acc_old_x_target];
        let group_acc_old_y = witness.wire_values[&group_acc_old_y_target];
        let group_acc_old = AffinePoint::<C>::nonzero(group_acc_old_x, group_acc_old_y);

        let scalar_acc_old = witness.wire_values[&scalar_acc_old_target];

        let addend_x = witness.wire_values[&addend_x_target];
        let addend_y = witness.wire_values[&addend_y_target];
        let addend = AffinePoint::<C>::nonzero(addend_x, addend_y);

        let scalar_bit = witness.wire_values[&scalar_bit_target];
        debug_assert!(scalar_bit.is_zero() || scalar_bit.is_one());

        let mut group_acc_new = group_acc_old;
        if scalar_bit.is_one() {
            group_acc_new = (group_acc_new + addend).to_affine();
        }

        let scalar_acc_new = scalar_acc_old.double() + scalar_bit;

        // Here's where our abstraction leaks a bit. Although we already have the sum, we need to
        // redo part of the computation in order to populate the purported inverse wire.
        let dx = group_acc_old_x - addend_x;
        let dy = group_acc_old_y - addend_y;
        let inverse = dx.multiplicative_inverse().expect("x_1 = x_2");

        let mut result = PartialWitness::new();
        result.wire_values.insert(group_acc_new_x_target, group_acc_new.x);
        result.wire_values.insert(group_acc_new_y_target, group_acc_new.y);
        result.wire_values.insert(scalar_acc_new_target, scalar_acc_new);
        result.wire_values.insert(inverse_target, inverse);
        result
    }
}

/// A curve which performs point doubling.
pub(crate) struct CurveDblGate<C: Curve> {
    index: usize,
    _phantom: PhantomData<C>,
}

impl<C: Curve> CurveDblGate<C> {
    const WIRE_X_OLD: usize = 0;
    const WIRE_Y_OLD: usize = 1;
    const WIRE_X_NEW: usize = 2;
    const WIRE_Y_NEW: usize = 3;
    const WIRE_INVERSE: usize = 4;
}

impl<C: Curve> Gate<C::BaseField> for CurveDblGate<C> {
    const ID: usize = 2;
}

impl<C: Curve> WitnessGenerator<C::BaseField> for CurveDblGate<C> {
    fn dependencies(&self) -> Vec<GateInput> {
        vec![
            GateInput { gate: self.index, input: Self::WIRE_X_OLD },
            GateInput { gate: self.index, input: Self::WIRE_Y_OLD },
        ]
    }

    fn generate(&self, witness: &PartialWitness<<C as Curve>::BaseField>) -> PartialWitness<<C as Curve>::BaseField> {
        let x_old_target = GateInput { gate: self.index, input: Self::WIRE_X_OLD };
        let y_old_target = GateInput { gate: self.index, input: Self::WIRE_Y_OLD };
        let x_new_target = GateInput { gate: self.index, input: Self::WIRE_X_NEW };
        let y_new_target = GateInput { gate: self.index, input: Self::WIRE_Y_NEW };
        let inverse_target = GateInput { gate: self.index, input: Self::WIRE_INVERSE };

        let x_old = witness.wire_values[&x_old_target];
        let y_old = witness.wire_values[&y_old_target];
        let old = AffinePoint::<C>::nonzero(x_old, y_old);
        let new = old.double();

        // Here's where our abstraction leaks a bit. Although we already have the result, we need to
        // redo part of the computation in order to populate the purported inverse wire.
        let inverse = y_old.double().multiplicative_inverse().expect("y = 0");

        let mut result = PartialWitness::new();
        result.wire_values.insert(inverse_target, inverse);
        result.wire_values.insert(x_new_target, new.x);
        result.wire_values.insert(y_new_target, new.y);
        result
    }
}

/// A gate which performs an iteration of an simultaneous doubling MSM loop, employing the
/// endomorphism described in the Halo paper. `C` is the curve of the inner proof.
pub(crate) struct CurveEndoGate<C: HaloEndomorphismCurve> {
    index: usize,
    _phantom: PhantomData<C>,
}

impl<C: HaloEndomorphismCurve> CurveEndoGate<C> {
    const WIRE_GROUP_ACC_X: usize = 0;
    const WIRE_GROUP_ACC_Y: usize = 1;
    const WIRE_SCALAR_ACC_UNSIGNED: usize = 2;
    const WIRE_SCALAR_ACC_SIGNED: usize = 3;
    const WIRE_ADDEND_X: usize = 4;
    const WIRE_ADDEND_Y: usize = 5;
    const WIRE_SCALAR_BIT_0: usize = 6;
    const WIRE_SCALAR_BIT_1: usize = 7;
    const WIRE_INVERSE: usize = 8;
}

impl<C: HaloEndomorphismCurve> Gate<C::BaseField> for CurveEndoGate<C> {
    const ID: usize = 3;
}

impl<C: HaloEndomorphismCurve> WitnessGenerator<C::BaseField> for CurveEndoGate<C> {
    fn dependencies(&self) -> Vec<GateInput> {
        vec![
            GateInput { gate: self.index, input: Self::WIRE_GROUP_ACC_X },
            GateInput { gate: self.index, input: Self::WIRE_GROUP_ACC_Y },
            GateInput { gate: self.index, input: Self::WIRE_SCALAR_ACC_UNSIGNED },
            GateInput { gate: self.index, input: Self::WIRE_SCALAR_ACC_SIGNED },
            GateInput { gate: self.index, input: Self::WIRE_ADDEND_X },
            GateInput { gate: self.index, input: Self::WIRE_ADDEND_Y },
            GateInput { gate: self.index, input: Self::WIRE_SCALAR_BIT_0 },
            GateInput { gate: self.index, input: Self::WIRE_SCALAR_BIT_1 },
        ]
    }

    fn generate(&self, witness: &PartialWitness<C::BaseField>) -> PartialWitness<C::BaseField> {
        let group_acc_old_x_target = GateInput { gate: self.index, input: Self::WIRE_GROUP_ACC_X };
        let group_acc_new_x_target = GateInput { gate: self.index + 1, input: Self::WIRE_GROUP_ACC_X };
        let group_acc_old_y_target = GateInput { gate: self.index, input: Self::WIRE_GROUP_ACC_Y };
        let group_acc_new_y_target = GateInput { gate: self.index + 1, input: Self::WIRE_GROUP_ACC_Y };

        let scalar_acc_unsigned_old_target = GateInput { gate: self.index, input: Self::WIRE_GROUP_ACC_Y };
        let scalar_acc_unsigned_new_target = GateInput { gate: self.index + 1, input: Self::WIRE_GROUP_ACC_Y };
        let scalar_acc_signed_old_target = GateInput { gate: self.index, input: Self::WIRE_GROUP_ACC_Y };
        let scalar_acc_signed_new_target = GateInput { gate: self.index + 1, input: Self::WIRE_GROUP_ACC_Y };

        let addend_x_target = GateInput { gate: self.index, input: Self::WIRE_ADDEND_X };
        let addend_y_target = GateInput { gate: self.index, input: Self::WIRE_ADDEND_Y };
        let scalar_bit_0_target = GateInput { gate: self.index, input: Self::WIRE_SCALAR_BIT_0 };
        let scalar_bit_1_target = GateInput { gate: self.index, input: Self::WIRE_SCALAR_BIT_1 };
        let inverse_target = GateInput { gate: self.index, input: Self::WIRE_INVERSE };

        let group_acc_old_x = witness.wire_values[&group_acc_old_x_target];
        let group_acc_old_y = witness.wire_values[&group_acc_old_y_target];
        let group_acc_old = AffinePoint::<C>::nonzero(group_acc_old_x, group_acc_old_y);

        let scalar_acc_unsigned_old = witness.wire_values[&scalar_acc_unsigned_old_target];
        let scalar_acc_signed_old = witness.wire_values[&scalar_acc_signed_old_target];

        let scalar_bit_0 = witness.wire_values[&scalar_bit_0_target];
        let scalar_bit_1 = witness.wire_values[&scalar_bit_1_target];
        debug_assert!(scalar_bit_0.is_zero() || scalar_bit_0.is_one());
        debug_assert!(scalar_bit_1.is_zero() || scalar_bit_1.is_one());

        let p_x = witness.wire_values[&addend_x_target];
        let p_y = witness.wire_values[&addend_y_target];

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
        result.wire_values.insert(group_acc_new_x_target, group_acc_new.x);
        result.wire_values.insert(group_acc_new_y_target, group_acc_new.y);
        result.wire_values.insert(scalar_acc_unsigned_new_target, scalar_acc_unsigned_new);
        result.wire_values.insert(scalar_acc_signed_new_target, scalar_acc_signed_new);
        result.wire_values.insert(inverse_target, inverse);

        result
    }
}

pub(crate) struct RescueGate { index: usize }

impl RescueGate {
    const WIRE_INPUT_0: usize = 0;
    const WIRE_INPUT_1: usize = 1;
    const WIRE_INPUT_2: usize = 2;
    const WIRE_INPUT_3: usize = 3;
    const WIRE_STEP_0: usize = 4;
    const WIRE_STEP_1: usize = 5;
    const WIRE_STEP_2: usize = 6;
    const WIRE_STEP_3: usize = 7;

    const MDS: [[u64; 4]; 4] = [[2, 3, 1, 1], [1, 2, 3, 1], [1, 1, 2, 3], [3, 1, 1, 2]];

    fn mds<F: Field>() -> [[F; 4]; 4] {
        let mut result = [[F::ZERO; 4]; 4];
        for r in 0..4 {
            for c in 0..4 {
                result[r][c] = F::from_canonical_u64(Self::MDS[r][c]);
            }
        }
        result
    }
}

impl<F: Field> Gate<F> for RescueGate {
    const ID: usize = 4;
}

impl<F: Field> WitnessGenerator<F> for RescueGate {
    fn dependencies(&self) -> Vec<GateInput> {
        vec![
            GateInput { gate: self.index, input: Self::WIRE_INPUT_0 },
            GateInput { gate: self.index, input: Self::WIRE_INPUT_1 },
            GateInput { gate: self.index, input: Self::WIRE_INPUT_2 },
            GateInput { gate: self.index, input: Self::WIRE_INPUT_3 },
        ]
    }

    fn generate(&self, witness: &PartialWitness<F>) -> PartialWitness<F> {
        let in_0_target = GateInput { gate: self.index, input: Self::WIRE_INPUT_0 };
        let in_1_target = GateInput { gate: self.index, input: Self::WIRE_INPUT_1 };
        let in_2_target = GateInput { gate: self.index, input: Self::WIRE_INPUT_2 };
        let in_3_target = GateInput { gate: self.index, input: Self::WIRE_INPUT_3 };

        let step_0_target = GateInput { gate: self.index, input: Self::WIRE_STEP_0 };
        let step_1_target = GateInput { gate: self.index, input: Self::WIRE_STEP_1 };
        let step_2_target = GateInput { gate: self.index, input: Self::WIRE_STEP_2 };
        let step_3_target = GateInput { gate: self.index, input: Self::WIRE_STEP_3 };

        let out_0_target = GateInput { gate: self.index + 1, input: Self::WIRE_INPUT_0 };
        let out_1_target = GateInput { gate: self.index + 1, input: Self::WIRE_INPUT_1 };
        let out_2_target = GateInput { gate: self.index + 1, input: Self::WIRE_INPUT_2 };
        let out_3_target = GateInput { gate: self.index + 1, input: Self::WIRE_INPUT_3 };

        let in_0 = witness.wire_values[&in_0_target];
        let in_1 = witness.wire_values[&in_1_target];
        let in_2 = witness.wire_values[&in_2_target];
        let in_3 = witness.wire_values[&in_3_target];

        let in_0_cubed = in_0.cube();
        let in_1_cubed = in_1.cube();
        let in_2_cubed = in_2.cube();
        let in_3_cubed = in_3.cube();

        let mds = Self::mds::<F>();
        let step_0 = mds[0][0] * in_0_cubed + mds[0][1] * in_1_cubed + mds[0][2] * in_2_cubed + mds[0][3] * in_3_cubed;
        let step_1 = mds[1][0] * in_0_cubed + mds[1][1] * in_1_cubed + mds[1][2] * in_2_cubed + mds[1][3] * in_3_cubed;
        let step_2 = mds[2][0] * in_0_cubed + mds[2][1] * in_1_cubed + mds[2][2] * in_2_cubed + mds[2][3] * in_3_cubed;
        let step_3 = mds[3][0] * in_0_cubed + mds[3][1] * in_1_cubed + mds[3][2] * in_2_cubed + mds[3][3] * in_3_cubed;

        let step_0_cubed = step_0.cube();
        let step_1_cubed = step_1.cube();
        let step_2_cubed = step_2.cube();
        let step_3_cubed = step_3.cube();

        let out_0 = mds[0][0] * step_0_cubed + mds[0][1] * step_1_cubed + mds[0][2] * step_2_cubed + mds[0][3] * step_3_cubed;
        let out_1 = mds[1][0] * step_0_cubed + mds[1][1] * step_1_cubed + mds[1][2] * step_2_cubed + mds[1][3] * step_3_cubed;
        let out_2 = mds[2][0] * step_0_cubed + mds[2][1] * step_1_cubed + mds[2][2] * step_2_cubed + mds[2][3] * step_3_cubed;
        let out_3 = mds[3][0] * step_0_cubed + mds[3][1] * step_1_cubed + mds[3][2] * step_2_cubed + mds[3][3] * step_3_cubed;

        let mut result = PartialWitness::new();

        // Write step wire values.
        result.wire_values.insert(step_0_target, step_0);
        result.wire_values.insert(step_1_target, step_1);
        result.wire_values.insert(step_2_target, step_2);
        result.wire_values.insert(step_3_target, step_3);

        // Write output wire values.
        result.wire_values.insert(out_0_target, out_0);
        result.wire_values.insert(out_1_target, out_1);
        result.wire_values.insert(out_2_target, out_2);
        result.wire_values.insert(out_3_target, out_3);

        result
    }
}

pub(crate) struct Base4SumGate { index: usize }

impl Base4SumGate {
    const WIRE_ACC: usize = 0;
    const WIRE_LIMB_0: usize = 1;
    const WIRE_LIMB_1: usize = 2;
    const WIRE_LIMB_2: usize = 3;
    const WIRE_LIMB_3: usize = 4;
    const WIRE_LIMB_4: usize = 5;
    const WIRE_LIMB_5: usize = 6;
    const WIRE_LIMB_6: usize = 7;
    const WIRE_LIMB_7: usize = 8;
}

impl<F: Field> Gate<F> for Base4SumGate {
    const ID: usize = 5;
}

impl<F: Field> WitnessGenerator<F> for Base4SumGate {
    fn dependencies(&self) -> Vec<GateInput> {
        Vec::new()
    }

    fn generate(&self, _witness: &PartialWitness<F>) -> PartialWitness<F> {
        // For base 4 decompositions, we don't do any witness generation on a per-gate level.
        // Instead, we have a single generator which generates values for an entire decomposition.
        PartialWitness::new()
    }
}

/// A "multiply, add, and rescale" gate.
pub(crate) struct MaddGate<F: Field> {
    index: usize,
    scalar: F
}

impl<F: Field> MaddGate<F> {
    const WIRE_MULTIPLICAND_0: usize = 0;
    const WIRE_MULTIPLICAND_1: usize = 1;
    const WIRE_ADDEND: usize = 2;
    const WIRE_OUTPUT: usize = 3;
}

impl<F: Field> Gate<F> for MaddGate<F> {
    const ID: usize = 6;
}

impl<F: Field> WitnessGenerator<F> for MaddGate<F> {
    fn dependencies(&self) -> Vec<GateInput> {
        vec![
            GateInput { gate: self.index, input: Self::WIRE_MULTIPLICAND_0 },
            GateInput { gate: self.index, input: Self::WIRE_MULTIPLICAND_1 },
            GateInput { gate: self.index, input: Self::WIRE_ADDEND },
        ]
    }

    fn generate(&self, witness: &PartialWitness<F>) -> PartialWitness<F> {
        let multiplicand_0_target = GateInput { gate: self.index, input: Self::WIRE_MULTIPLICAND_0 };
        let multiplicand_1_target = GateInput { gate: self.index, input: Self::WIRE_MULTIPLICAND_1 };
        let addend_target = GateInput { gate: self.index, input: Self::WIRE_ADDEND };
        let output_target = GateInput { gate: self.index, input: Self::WIRE_OUTPUT };

        let multiplicand_0 = witness.wire_values[&multiplicand_0_target];
        let multiplicand_1 = witness.wire_values[&multiplicand_1_target];
        let addend = witness.wire_values[&addend_target];

        let output = multiplicand_0 * multiplicand_1 + addend;

        let mut result = PartialWitness::new();
        result.wire_values.insert(output_target, output);
        result
    }
}
