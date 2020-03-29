use crate::{GateInput, WitnessGenerator, Field, PartialWitness};

trait Gate<F: Field>: WitnessGenerator<F> {
    const ID: usize;
}

struct NoopGate { index: usize }

impl NoopGate {
    const BUFFER_0: usize = 0;
    const BUFFER_1: usize = 1;
    const BUFFER_2: usize = 2;
    const BUFFER_3: usize = 3;
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

struct MsmGate { index: usize }

impl MsmGate {
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

impl<F: Field> Gate<F> for MsmGate {
    const ID: usize = 1;
}

impl<F: Field> WitnessGenerator<F> for MsmGate {
    fn dependencies(&self) -> Vec<GateInput> {
        todo!()
    }

    fn generate(&self, witness: &PartialWitness<F>) -> PartialWitness<F> {
        todo!()
    }
}

struct RescueGate { index: usize }

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
    const ID: usize = 2;
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

struct Base4SumGate { index: usize }

impl Base4SumGate {
    const ACC: usize = 0;
    const LIMB_0: usize = 1;
    const LIMB_1: usize = 2;
    const LIMB_2: usize = 3;
    const LIMB_3: usize = 4;
    const LIMB_4: usize = 5;
    const LIMB_5: usize = 6;
    const LIMB_6: usize = 7;
    const LIMB_7: usize = 8;
}

impl<F: Field> Gate<F> for Base4SumGate {
    const ID: usize = 3;
}

impl<F: Field> WitnessGenerator<F> for Base4SumGate {
    fn dependencies(&self) -> Vec<GateInput> {
        todo!()
    }

    fn generate(&self, witness: &PartialWitness<F>) -> PartialWitness<F> {
        todo!()
    }
}

struct MiscGate { index: usize }

impl MiscGate {
    const MULTIPLICAND_0: usize = 0;
    const MULTIPLICAND_1: usize = 1;
    const ADDEND: usize = 2;
}

impl<F: Field> Gate<F> for MiscGate {
    const ID: usize = 4;
}

impl<F: Field> WitnessGenerator<F> for MiscGate {
    fn dependencies(&self) -> Vec<GateInput> {
        unimplemented!()
    }

    fn generate(&self, witness: &PartialWitness<F>) -> PartialWitness<F> {
        unimplemented!()
    }
}
