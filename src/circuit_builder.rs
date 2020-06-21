use std::collections::{BTreeMap, HashMap};

use crate::gates::*;
use crate::plonk_util::{coeffs_to_values_padded, sigma_polynomials, values_to_coeffs};
use crate::poly_commit::PolynomialCommitment;
use crate::util::{ceil_div_usize, log2_strict, transpose};
use crate::{blake_hash_base_field_to_curve, blake_hash_usize_to_curve, fft_precompute, generate_rescue_constants, msm_precompute, AffinePoint, AffinePointTarget, Circuit, Curve, CurveMsmEndoResult, CurveMulEndoResult, CurveMulOp, Field, HaloCurve, PartialWitness, PublicInput, Target, TargetPartitions, VirtualTarget, Wire, WitnessGenerator, NUM_CONSTANTS, NUM_WIRES};

pub struct CircuitBuilder<C: HaloCurve> {
    pub(crate) security_bits: usize,
    public_input_index: usize,
    virtual_target_index: usize,
    gate_counts: BTreeMap<&'static str, usize>,
    gate_constants: Vec<Vec<C::ScalarField>>,
    copy_constraints: Vec<(Target, Target)>,
    generators: Vec<Box<dyn WitnessGenerator<C::ScalarField>>>,
    constant_wires: HashMap<C::ScalarField, Target>,
}

impl<C: HaloCurve> CircuitBuilder<C> {
    pub fn new(security_bits: usize) -> Self {
        CircuitBuilder {
            security_bits,
            public_input_index: 0,
            virtual_target_index: 0,
            gate_counts: BTreeMap::new(),
            gate_constants: Vec::new(),
            copy_constraints: Vec::new(),
            generators: Vec::new(),
            constant_wires: HashMap::new(),
        }
    }

    pub fn stage_public_input(&mut self) -> PublicInput {
        let index = self.public_input_index;
        self.public_input_index += 1;
        PublicInput { index }
    }

    pub fn stage_public_inputs(&mut self, n: usize) -> Vec<PublicInput> {
        (0..n).map(|_i| self.stage_public_input()).collect()
    }

    /// Add `PublicInputGate`s which enable public inputs to be routed. Should be called after all
    /// `stage_public_input[s]` calls, but before any gates are added.
    pub fn route_public_inputs(&mut self) {
        debug_assert_eq!(
            self.num_gates(),
            0,
            "Must be called before any gates are added"
        );
        let num_pi_gates = ceil_div_usize(self.public_input_index, NUM_WIRES);
        for i in 0..num_pi_gates {
            self.add_gate_no_constants(PublicInputGate::new(i * 2));
            self.add_gate_no_constants(BufferGate::new(i * 2 + 1));
        }
    }

    pub fn add_virtual_target(&mut self) -> Target {
        let index = self.virtual_target_index;
        self.virtual_target_index += 1;
        Target::VirtualTarget(VirtualTarget { index })
    }

    pub fn add_virtual_targets(&mut self, n: usize) -> Vec<Target> {
        (0..n).map(|_i| self.add_virtual_target()).collect()
    }

    pub fn add_virtual_point_target(&mut self) -> AffinePointTarget {
        let x = self.add_virtual_target();
        let y = self.add_virtual_target();
        AffinePointTarget { x, y }
    }

    pub fn add_virtual_point_targets(&mut self, n: usize) -> Vec<AffinePointTarget> {
        (0..n).map(|_i| self.add_virtual_point_target()).collect()
    }

    pub fn zero_wire(&mut self) -> Target {
        self.constant_wire(C::ScalarField::ZERO)
    }

    pub fn one_wire(&mut self) -> Target {
        self.constant_wire(C::ScalarField::ONE)
    }

    pub fn two_wire(&mut self) -> Target {
        self.constant_wire(C::ScalarField::TWO)
    }

    pub fn neg_one_wire(&mut self) -> Target {
        self.constant_wire(C::ScalarField::NEG_ONE)
    }

    pub fn constant_wire(&mut self, c: C::ScalarField) -> Target {
        if self.constant_wires.contains_key(&c) {
            self.constant_wires[&c]
        } else {
            let result = self.create_constant_wire(c);
            self.constant_wires.insert(c, result);
            result
        }
    }

    pub fn constant_wire_u32(&mut self, c: u32) -> Target {
        self.constant_wire(C::ScalarField::from_canonical_u32(c))
    }

    fn create_constant_wire(&mut self, c: C::ScalarField) -> Target {
        // We will create a ConstantGate and pass c as its first (and only) constant, which will
        // cause it to populate its output wire with the same value c.
        let gate = self.num_gates();
        self.add_gate(ConstantGate::new(gate), vec![c]);
        Target::Wire(Wire {
            gate,
            input: ConstantGate::<C>::WIRE_OUTPUT,
        })
    }

    pub fn constant_affine_point<InnerC: Curve<BaseField = C::ScalarField>>(
        &mut self,
        point: AffinePoint<InnerC>,
    ) -> AffinePointTarget {
        assert!(!point.zero);
        AffinePointTarget {
            x: self.constant_wire(point.x),
            y: self.constant_wire(point.y),
        }
    }

    pub fn assert_zero(&mut self, x: Target) {
        let zero = self.zero_wire();
        self.copy(x, zero);
    }

    pub fn assert_binary(&mut self, x: Target) {
        let _zero = self.zero_wire();
        let one = self.one_wire();

        let x_minus_1 = self.sub(x, one);
        let product = self.mul(x, x_minus_1);
        self.assert_zero(product);
    }

    pub fn add(&mut self, x: Target, y: Target) -> Target {
        let zero = self.zero_wire();
        if x == zero {
            return y;
        }
        if y == zero {
            return x;
        }

        let one = self.one_wire();
        let index = self.num_gates();
        self.add_gate(
            ArithmeticGate::new(index),
            vec![C::ScalarField::ONE, C::ScalarField::ONE],
        );
        self.copy(
            x,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_0,
            }),
        );
        self.copy(
            one,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_1,
            }),
        );
        self.copy(
            y,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_ADDEND,
            }),
        );
        Target::Wire(Wire {
            gate: index,
            input: ArithmeticGate::<C>::WIRE_OUTPUT,
        })
    }

    pub fn add_many(&mut self, terms: &[Target]) -> Target {
        let mut sum = self.zero_wire();
        for term in terms {
            sum = self.add(sum, *term);
        }
        sum
    }

    pub fn double(&mut self, x: Target) -> Target {
        self.add(x, x)
    }

    pub fn sub(&mut self, x: Target, y: Target) -> Target {
        let zero = self.zero_wire();
        if y == zero {
            return x;
        }

        let one = self.one_wire();
        let index = self.num_gates();
        self.add_gate(
            ArithmeticGate::new(index),
            vec![C::ScalarField::ONE, C::ScalarField::NEG_ONE],
        );
        self.copy(
            x,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_0,
            }),
        );
        self.copy(
            one,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_1,
            }),
        );
        self.copy(
            y,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_ADDEND,
            }),
        );
        Target::Wire(Wire {
            gate: index,
            input: ArithmeticGate::<C>::WIRE_OUTPUT,
        })
    }

    pub fn mul(&mut self, x: Target, y: Target) -> Target {
        let one = self.one_wire();
        if x == one {
            return y;
        }
        if y == one {
            return x;
        }

        let zero = self.zero_wire();
        let index = self.num_gates();
        self.add_gate(
            ArithmeticGate::new(index),
            vec![C::ScalarField::ONE, C::ScalarField::ZERO],
        );
        self.copy(
            x,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_0,
            }),
        );
        self.copy(
            y,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_1,
            }),
        );
        self.copy(
            zero,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_ADDEND,
            }),
        );
        Target::Wire(Wire {
            gate: index,
            input: ArithmeticGate::<C>::WIRE_OUTPUT,
        })
    }

    pub fn mul_many(&mut self, terms: &[Target]) -> Target {
        let mut product = self.one_wire();
        for term in terms {
            product = self.mul(product, *term);
        }
        product
    }

    pub fn square(&mut self, x: Target) -> Target {
        self.mul(x, x)
    }

    /// Note: This assumes the most significant bit of each scalar is unset. This occurs with high
    /// probability if the field size is slightly larger than a power of two and the inputs are
    /// uniformly random.
    pub fn deterministic_square_root(&mut self, x: Target) -> Target {
        // Assume x != 0. Let y, z be the square roots of x. Since y + z = |F|, and |F| is odd, the
        // parity of y and z must differ, so we can enforce determinism by checking for a certain
        // parity bit state. We chose a parity bit of 0 since this also works for the x = 0 case.

        struct SqrtGenerator {
            x: Target,
            x_sqrt: Target,
        }

        impl<F: Field> WitnessGenerator<F> for SqrtGenerator {
            fn dependencies(&self) -> Vec<Target> {
                vec![self.x]
            }

            fn generate(
                &self,
                _constants: &Vec<Vec<F>>,
                witness: &PartialWitness<F>,
            ) -> PartialWitness<F> {
                let x_value = witness.get_target(self.x);
                let mut x_sqrt_value = x_value.square_root().expect("Not square");

                if x_sqrt_value.to_canonical_bool_vec()[0] {
                    // The parity bit is 1; we want the other square root.
                    x_sqrt_value = -x_sqrt_value;
                    debug_assert!(!x_sqrt_value.to_canonical_bool_vec()[0]);
                }

                let mut result = PartialWitness::new();
                result.set_target(self.x_sqrt, x_sqrt_value);
                result
            }
        }

        let x_sqrt = self.add_virtual_target();
        self.add_generator(SqrtGenerator { x, x_sqrt });

        // We assume each most significant bit is unset; see the note in the method doc.
        let f_bits = C::ScalarField::BITS - 1;
        assert_eq!(
            f_bits, 254,
            "We currently only handle fields of size 2^254 + epsilon"
        );
        let (bits, dibits) = self.split_binary_and_base_4(x_sqrt, 2, 126);

        // Verify that x_sqrt * x_sqrt = x.
        let x_sqrt_squared = self.square(x_sqrt);
        self.copy(x_sqrt_squared, x);

        // Verify that the parity bit is 0, and the other bit is binary.
        self.assert_zero(bits[0]);
        self.assert_binary(bits[1]);

        // Verify the decomposition by taking a weighted sum of the limbs. Since the first bit is
        // always zero, we start with the second bit (scaled by its weight of 2) and then add the
        // base 4 limbs.
        let mut sum = self.double(bits[1]);
        for chunk in dibits.chunks(Base4SumGate::<C>::NUM_LIMBS) {
            assert_eq!(
                chunk.len(),
                Base4SumGate::<C>::NUM_LIMBS,
                "Should not have a partial chunk"
            );

            let index = self.num_gates();
            self.add_gate_no_constants(Base4SumGate::new(index));
            for i in 0..chunk.len() {
                self.copy(
                    sum,
                    Target::Wire(Wire {
                        gate: index,
                        input: Base4SumGate::<C>::WIRE_ACC_OLD,
                    }),
                );
                self.copy(
                    chunk[i],
                    Target::Wire(Wire {
                        gate: index,
                        input: Base4SumGate::<C>::WIRE_LIMB_0 + i,
                    }),
                );
                sum = Target::Wire(Wire {
                    gate: index,
                    input: Base4SumGate::<C>::WIRE_ACC_NEW,
                })
            }
        }
        self.copy(sum, x);

        x_sqrt
    }

    /// Compute `x^power`, where `power` is a constant.
    pub fn exp_constant(&mut self, x: Target, power: C::ScalarField) -> Target {
        let power_bits = power.num_bits();
        let mut current = x;
        let mut product = self.one_wire();

        for (i, limb) in power.to_canonical_u64_vec().iter().enumerate() {
            for j in 0..64 {
                // If we've gone through all the 1 bits already, no need to keep squaring.
                let bit_index = i * 64 + j;
                if bit_index == power_bits {
                    return product;
                }

                if (limb >> j & 1) != 0 {
                    product = self.mul(product, current);
                }
                current = self.square(current);
            }
        }

        product
    }

    /// Compute `x^power`, where `power` is a constant `usize`.
    pub fn exp_constant_usize(&mut self, x: Target, power: usize) -> Target {
        self.exp_constant(x, C::ScalarField::from_canonical_usize(power))
    }

    pub fn inv(&mut self, x: Target) -> Target {
        struct InverseGenerator {
            x: Target,
            x_inv: Target,
        }

        impl<F: Field> WitnessGenerator<F> for InverseGenerator {
            fn dependencies(&self) -> Vec<Target> {
                vec![self.x]
            }

            fn generate(
                &self,
                _constants: &Vec<Vec<F>>,
                witness: &PartialWitness<F>,
            ) -> PartialWitness<F> {
                let x_value = witness.get_target(self.x);
                let x_inv_value = x_value.multiplicative_inverse().expect("x = 0");

                let mut result = PartialWitness::new();
                result.set_target(self.x_inv, x_inv_value);
                result
            }
        }

        let x_inv = self.add_virtual_target();
        self.add_generator(InverseGenerator { x, x_inv });

        // Enforce that x * x_inv = 1.
        let product = self.mul(x, x_inv);
        let one = self.one_wire();
        self.copy(product, one);

        x_inv
    }

    pub fn div(&mut self, x: Target, y: Target) -> Target {
        let y_inv = self.inv(y);
        self.mul(x, y_inv)
    }

    /// Multiply and add; i.e. computes `x * y + z`.
    pub fn mul_add(&mut self, x: Target, y: Target, z: Target) -> Target {
        let index = self.num_gates();
        self.add_gate(
            ArithmeticGate::new(index),
            vec![C::ScalarField::ONE, C::ScalarField::ONE],
        );
        self.copy(
            x,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_0,
            }),
        );
        self.copy(
            y,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_1,
            }),
        );
        self.copy(
            z,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_ADDEND,
            }),
        );
        Target::Wire(Wire {
            gate: index,
            input: ArithmeticGate::<C>::WIRE_OUTPUT,
        })
    }

    /// Multiply and subtract; i.e. computes `x * y - z`.
    pub fn mul_sub(&mut self, x: Target, y: Target, z: Target) -> Target {
        let index = self.num_gates();
        self.add_gate(
            ArithmeticGate::new(index),
            vec![C::ScalarField::ONE, C::ScalarField::NEG_ONE],
        );
        self.copy(
            x,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_0,
            }),
        );
        self.copy(
            y,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_1,
            }),
        );
        self.copy(
            z,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_ADDEND,
            }),
        );
        Target::Wire(Wire {
            gate: index,
            input: ArithmeticGate::<C>::WIRE_OUTPUT,
        })
    }

    /// Computes `-x`.
    pub fn neg(&mut self, x: Target) -> Target {
        let neg_one = self.neg_one_wire();
        self.mul(x, neg_one)
    }

    /// Splits `x` into its binary representation. Note that this method merely adds a generator to
    /// populate the bit wires; it does not enforce constraints to verify the decomposition.
    fn split_binary(&mut self, x: Target, num_bits: usize) -> Vec<Target> {
        let (bits, _dibits) = self.split_binary_and_base_4(x, num_bits, 0);
        bits
    }

    /// Splits `x` into a combination of binary and base 4 limbs. Note that this method merely adds
    /// a generator to populate the limb wires; it does not enforce constraints to verify the
    /// decomposition.
    fn split_binary_and_base_4(
        &mut self,
        x: Target,
        num_bits: usize,
        num_dibits: usize,
    ) -> (Vec<Target>, Vec<Target>) {
        struct SplitGenerator {
            x: Target,
            bits: Vec<Target>,
            dibits: Vec<Target>,
        }

        impl<F: Field> WitnessGenerator<F> for SplitGenerator {
            fn dependencies(&self) -> Vec<Target> {
                vec![self.x]
            }

            fn generate(
                &self,
                _constants: &Vec<Vec<F>>,
                witness: &PartialWitness<F>,
            ) -> PartialWitness<F> {
                let x = witness.get_target(self.x);
                let x_bits = x.to_canonical_bool_vec();

                let mut result = PartialWitness::new();
                for i in 0..self.bits.len() {
                    result.set_target(self.bits[i], F::from_canonical_bool(x_bits[i]));
                }
                for i in 0..self.dibits.len() {
                    let bit_1 = x_bits[self.bits.len() + i * 2];
                    let bit_2 = x_bits[self.bits.len() + i * 2 + 1];
                    let dibit = if bit_1 { 1 } else { 0 } + if bit_2 { 2 } else { 0 };
                    result.set_target(self.dibits[i], F::from_canonical_u32(dibit));
                }
                result
            }
        }

        let bits = self.add_virtual_targets(num_bits);
        let dibits = self.add_virtual_targets(num_dibits);
        let generator = SplitGenerator {
            x,
            bits: bits.clone(),
            dibits: dibits.clone(),
        };
        self.add_generator(generator);
        (bits, dibits)
    }

    pub fn rescue_hash_n_to_1(&mut self, inputs: &[Target]) -> Target {
        self.rescue_sponge(inputs, 1)[0]
    }

    pub fn rescue_hash_n_to_2(&mut self, inputs: &[Target]) -> (Target, Target) {
        let outputs = self.rescue_sponge(inputs, 2);
        (outputs[0], outputs[1])
    }

    pub fn rescue_hash_n_to_3(&mut self, inputs: &[Target]) -> (Target, Target, Target) {
        let outputs = self.rescue_sponge(inputs, 3);
        (outputs[0], outputs[1], outputs[2])
    }

    pub fn rescue_sponge(&mut self, inputs: &[Target], num_outputs: usize) -> Vec<Target> {
        let zero = self.zero_wire();
        let mut state = vec![zero; RESCUE_SPONGE_WIDTH];

        // Absorb all input chunks.
        for input_chunk in inputs.chunks(RESCUE_SPONGE_WIDTH - 1) {
            for i in 0..input_chunk.len() {
                state[i] = self.add(state[i], input_chunk[i]);
            }
            state = self.rescue_permutation(&state);
        }

        // Squeeze until we have the desired number of outputs.
        let mut outputs = Vec::new();
        loop {
            for i in 0..(RESCUE_SPONGE_WIDTH - 1) {
                outputs.push(state[i]);
                if outputs.len() == num_outputs {
                    return outputs;
                }
            }
            state = self.rescue_permutation(&state);
        }
        outputs
    }

    pub fn rescue_permutation(&mut self, inputs: &[Target]) -> Vec<Target> {
        assert_eq!(inputs.len(), RESCUE_SPONGE_WIDTH);

        // Route the input wires.
        for i in 0..RESCUE_SPONGE_WIDTH {
            self.copy(
                inputs[i],
                Target::Wire(Wire {
                    gate: self.num_gates(),
                    input: RescueStepAGate::<C>::wire_acc(i),
                }),
            );
        }

        let all_constants = generate_rescue_constants(RESCUE_SPONGE_WIDTH, self.security_bits);
        for (a_constants, b_constants) in all_constants.into_iter() {
            let a_index = self.num_gates();
            let a_gate = RescueStepAGate::new(a_index);
            self.add_gate(a_gate, a_constants);

            let b_index = self.num_gates();
            let b_gate = RescueStepBGate::new(b_index);
            self.add_gate(b_gate, b_constants);
        }

        // Use a BufferGate to receive the final accumulator states.
        let gate = self.num_gates();
        self.add_gate_no_constants(BufferGate::new(gate));
        (0..RESCUE_SPONGE_WIDTH)
            .map(|i| {
                Target::Wire(Wire {
                    gate,
                    input: RescueStepBGate::<C>::wire_acc(i),
                })
            })
            .collect()
    }

    /// Assert that a given coordinate pair is on the curve `C`.
    pub fn curve_assert_valid<InnerC: Curve<BaseField = C::ScalarField>>(
        &mut self,
        p: AffinePointTarget,
    ) {
        // Recall the short Weierstrass equation: y^2 = x^3 + a*x + b.
        let a = self.constant_wire(InnerC::A);
        let b = self.constant_wire(InnerC::B);
        let y_squared = self.square(p.y);
        let x_cubed = self.exp_constant_usize(p.x, 3);
        let a_x_plus_b = self.mul_add(a, p.x, b);
        let rhs = self.add(x_cubed, a_x_plus_b);
        self.copy(y_squared, rhs);
    }

    pub fn curve_neg<InnerC: Curve<BaseField = C::ScalarField>>(
        &mut self,
        p: AffinePointTarget,
    ) -> AffinePointTarget {
        let neg_y = self.neg(p.y);
        AffinePointTarget { x: p.x, y: neg_y }
    }

    pub fn curve_add<InnerC: Curve<BaseField = C::ScalarField>>(
        &mut self,
        _p_1: AffinePointTarget,
        _p_2: AffinePointTarget,
    ) -> AffinePointTarget {
        let add_index = self.num_gates();
        self.add_gate_no_constants(CurveAddGate::<C, InnerC>::new(add_index));
        let buffer_index = self.num_gates();
        self.add_gate_no_constants(BufferGate::new(buffer_index));

        // TODO: Wiring.

        let result_x = Target::Wire(Wire {
            gate: buffer_index,
            input: CurveAddGate::<C, InnerC>::WIRE_GROUP_ACC_X,
        });
        let result_y = Target::Wire(Wire {
            gate: buffer_index,
            input: CurveAddGate::<C, InnerC>::WIRE_GROUP_ACC_Y,
        });
        AffinePointTarget {
            x: result_x,
            y: result_y,
        }
    }

    pub fn curve_double<InnerC: Curve<BaseField = C::ScalarField>>(
        &mut self,
        p: AffinePointTarget,
    ) -> AffinePointTarget {
        let idx_dbl = self.num_gates();
        self.add_gate_no_constants(CurveDblGate::<C, InnerC>::new(idx_dbl));
        self.copy(
            p.x,
            Target::Wire(Wire {
                gate: idx_dbl,
                input: CurveDblGate::<C, InnerC>::WIRE_X_OLD,
            }),
        );
        self.copy(
            p.y,
            Target::Wire(Wire {
                gate: idx_dbl,
                input: CurveDblGate::<C, InnerC>::WIRE_Y_OLD,
            }),
        );
        AffinePointTarget {
            x: Target::Wire(Wire {
                gate: idx_dbl,
                input: CurveDblGate::<C, InnerC>::WIRE_X_NEW,
            }),
            y: Target::Wire(Wire {
                gate: idx_dbl,
                input: CurveDblGate::<C, InnerC>::WIRE_Y_NEW,
            }),
        }
    }

    pub fn curve_sub<InnerC: Curve<BaseField = C::ScalarField>>(
        &mut self,
        p_1: AffinePointTarget,
        p_2: AffinePointTarget,
    ) -> AffinePointTarget {
        let neg_p_2 = self.curve_neg::<InnerC>(p_2);
        self.curve_add::<InnerC>(p_1, neg_p_2)
    }

    pub fn curve_mul<InnerC: Curve<BaseField = C::ScalarField>>(
        &mut self,
        mul: CurveMulOp,
    ) -> AffinePointTarget {
        self.curve_msm::<InnerC>(&[mul])
    }

    pub fn curve_mul_endo<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &mut self,
        mul: CurveMulOp,
    ) -> CurveMulEndoResult {
        let result = self.curve_msm_endo::<InnerC>(&[mul]);
        CurveMulEndoResult {
            mul_result: result.msm_result,
            actual_scalar: result.actual_scalars[0],
        }
    }

    /// Note: This assumes the most significant bit of each scalar is unset. This occurs with high
    /// probability if the field size is slightly larger than a power of two and the scalars are
    /// uniformly random.
    pub fn curve_msm<InnerC: Curve<BaseField = C::ScalarField>>(
        &mut self,
        parts: &[CurveMulOp],
    ) -> AffinePointTarget {
        // We assume each most significant bit is unset; see the note in the method doc.
        let f_bits = C::ScalarField::BITS - 1;

        let all_bits: Vec<Vec<Target>> = parts
            .iter()
            .map(|part| self.split_binary(part.scalar, f_bits))
            .collect();

        // Normally we would start with zero, but to avoid exceptional cases, we start with some
        // random nonzero point and subtract it later. This avoids exceptional cases with high
        // probability. A malicious prover may be able to craft an input which leads to an
        // exceptional case, but this isn't a problem as our curve gates will be unsatisfiable in
        // exceptional cases.
        let mut filler = blake_hash_base_field_to_curve::<InnerC>(InnerC::BaseField::ZERO);
        let mut acc = self.constant_affine_point(filler);
        let mut scalar_accs = vec![self.zero_wire(); parts.len()];

        for i in (0..f_bits).rev() {
            // Route the accumulator to the first curve addition gate's inputs.
            self.copy(
                acc.x,
                Target::Wire(Wire {
                    gate: self.num_gates(),
                    input: CurveAddGate::<C, InnerC>::WIRE_GROUP_ACC_X,
                }),
            );
            self.copy(
                acc.y,
                Target::Wire(Wire {
                    gate: self.num_gates(),
                    input: CurveAddGate::<C, InnerC>::WIRE_GROUP_ACC_Y,
                }),
            );

            for (j, part) in parts.iter().enumerate() {
                let bit = all_bits[j][i];

                let idx_add = self.num_gates();
                self.add_gate_no_constants(CurveAddGate::<C, InnerC>::new(idx_add));
                self.copy(
                    scalar_accs[j],
                    Target::Wire(Wire {
                        gate: idx_add,
                        input: CurveAddGate::<C, InnerC>::WIRE_SCALAR_ACC_OLD,
                    }),
                );
                scalar_accs[j] = Target::Wire(Wire {
                    gate: idx_add,
                    input: CurveAddGate::<C, InnerC>::WIRE_SCALAR_ACC_NEW,
                });
                self.copy(
                    part.point.x,
                    Target::Wire(Wire {
                        gate: idx_add,
                        input: CurveAddGate::<C, InnerC>::WIRE_ADDEND_X,
                    }),
                );
                self.copy(
                    part.point.y,
                    Target::Wire(Wire {
                        gate: idx_add,
                        input: CurveAddGate::<C, InnerC>::WIRE_ADDEND_Y,
                    }),
                );
                self.copy(
                    bit,
                    Target::Wire(Wire {
                        gate: idx_add,
                        input: CurveAddGate::<C, InnerC>::WIRE_SCALAR_BIT,
                    }),
                );
            }

            // Double the accumulator.
            let idx_dbl = self.num_gates();
            self.add_gate_no_constants(CurveDblGate::<C, InnerC>::new(idx_dbl));
            // No need to route the double gate's inputs, because the last add gate would have
            // constrained them.
            // Normally, we will take the double gate's outputs as the new accumulator. If we just
            // completed the last iteration of the MSM though, then we don't want to perform a final
            // doubling, so we will take its inputs as the result instead.
            if i == f_bits - 1 {
                acc = AffinePointTarget {
                    x: Target::Wire(Wire {
                        gate: idx_dbl,
                        input: CurveDblGate::<C, InnerC>::WIRE_X_OLD,
                    }),
                    y: Target::Wire(Wire {
                        gate: idx_dbl,
                        input: CurveDblGate::<C, InnerC>::WIRE_Y_OLD,
                    }),
                };
            } else {
                acc = AffinePointTarget {
                    x: Target::Wire(Wire {
                        gate: idx_dbl,
                        input: CurveDblGate::<C, InnerC>::WIRE_X_NEW,
                    }),
                    y: Target::Wire(Wire {
                        gate: idx_dbl,
                        input: CurveDblGate::<C, InnerC>::WIRE_Y_NEW,
                    }),
                };
            }

            // Also double the filler, so we can subtract out a rescaled version later.
            filler = filler.double();
        }

        // Subtract (a rescaled version of) the arbitrary nonzero value that we started with.
        let filler_target = self.constant_affine_point(filler);
        acc = self.curve_sub::<InnerC>(acc, filler_target);

        // Assert that each accumulation of scalar bits matches the original scalar.
        for (j, part) in parts.iter().enumerate() {
            self.copy(scalar_accs[j], part.scalar);
        }

        acc
    }

    /// Like `curve_msm`, but uses the endomorphism described in the Halo paper.
    pub fn curve_msm_endo<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &mut self,
        parts: &[CurveMulOp],
    ) -> CurveMsmEndoResult {
        let zero = self.zero_wire();

        // We assume each most significant bit is unset; see the note in curve_msm's method doc.
        let f_bits = C::ScalarField::BITS - 1;
        let scalar_bits = self.security_bits;
        let scalar_dibits = (f_bits - scalar_bits) / 2;

        // To keep things simple for now, we only handle the case of |F| ~= 2^254 and lambda = 128.
        assert_eq!(f_bits, 254);
        assert_eq!(scalar_bits, 128);
        assert_eq!(scalar_dibits, 63);

        // We split each scalar into 128 bits and 63 dibits. The bits are used in the MSM, while the
        // dibits are ignored, except that we need to include them in our sum-of-limbs computation
        // in order to verify that the decomposition was correct.
        let (all_bits, all_dibits): (Vec<Vec<Target>>, Vec<Vec<Target>>) = parts
            .iter()
            .map(|part| self.split_binary_and_base_4(part.scalar, scalar_bits, scalar_dibits))
            .unzip();

        // Normally we would start with zero, but to avoid exceptional cases, we start with some
        // random nonzero point and subtract it later. This avoids exceptional cases with high
        // probability. A malicious prover may be able to craft an input which leads to an
        // exceptional case, but this isn't a problem as our curve gates will be unsatisfiable in
        // exceptional cases.
        let mut filler = blake_hash_base_field_to_curve::<InnerC>(InnerC::BaseField::ZERO);
        let mut acc = self.constant_affine_point(filler);

        // For each scalar, we maintain two accumulators. The unsigned one is for computing a
        // weighted sum of bits and dibits in the usual manner, so that we can later check that this
        // sum equals the original scalar. The signed one is for computing n(s) for each scalar s.
        // This is the "actual" scalar by which the associated point was multiplied, accounting for
        // the endomorphism.
        let mut scalar_acc_unsigned = Vec::new();
        let mut scalar_acc_signed = Vec::new();

        // As in the Halo paper, we process two scalar bits at a time.
        for i in (0..scalar_bits).step_by(2).rev() {
            // Route the point accumulator to the first gate's inputs.
            self.copy(
                acc.x,
                Target::Wire(Wire {
                    gate: self.num_gates(),
                    input: CurveEndoGate::<C, InnerC>::WIRE_GROUP_ACC_X,
                }),
            );
            self.copy(
                acc.y,
                Target::Wire(Wire {
                    gate: self.num_gates(),
                    input: CurveEndoGate::<C, InnerC>::WIRE_GROUP_ACC_Y,
                }),
            );

            for (j, part) in parts.iter().enumerate() {
                let bit_0 = all_bits[j][i];
                let bit_1 = all_bits[j][i + 1];

                let gate = self.num_gates();
                self.add_gate_no_constants(CurveEndoGate::<C, InnerC>::new(gate));

                self.copy(
                    part.point.x,
                    Target::Wire(Wire {
                        gate,
                        input: CurveEndoGate::<C, InnerC>::WIRE_ADDEND_X,
                    }),
                );
                self.copy(
                    part.point.y,
                    Target::Wire(Wire {
                        gate,
                        input: CurveEndoGate::<C, InnerC>::WIRE_ADDEND_Y,
                    }),
                );
                self.copy(
                    bit_0,
                    Target::Wire(Wire {
                        gate,
                        input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_BIT_0,
                    }),
                );
                self.copy(
                    bit_1,
                    Target::Wire(Wire {
                        gate,
                        input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_BIT_1,
                    }),
                );

                // If this is the first pair of scalar bits being processed, route 0 to the scalar accumulators.
                if i == scalar_bits - 2 {
                    self.copy(
                        zero,
                        Target::Wire(Wire {
                            gate,
                            input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_ACC_UNSIGNED,
                        }),
                    );
                    self.copy(
                        zero,
                        Target::Wire(Wire {
                            gate,
                            input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_ACC_SIGNED,
                        }),
                    );
                }

                // If this is the last pair of scalar bits being processed, save the final
                // accumulator states.
                // Since CurveEndoGate will store these in the "next" gate, but this is the last
                // CurveEndoGate for this scalar, we need to add an extra BufferGate to receive them.
                if i == 0 {
                    let gate = self.num_gates();
                    self.add_gate_no_constants(BufferGate::new(gate));
                    scalar_acc_unsigned.push(Target::Wire(Wire {
                        gate,
                        input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_ACC_UNSIGNED,
                    }));
                    scalar_acc_signed.push(Target::Wire(Wire {
                        gate,
                        input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_ACC_SIGNED,
                    }));
                }
            }

            // Double the accumulator.
            let gate = self.num_gates();
            self.add_gate_no_constants(CurveDblGate::<C, InnerC>::new(gate));
            // No need to route the double gate's inputs, because the last endo gate would have
            // constrained them.
            // Normally, we will take the double gate's outputs as the new accumulator. If we just
            // completed the last iteration of the MSM though, then we don't want to perform a final
            // doubling, so we will take its inputs as the result instead.
            if i == scalar_bits - 1 {
                acc = AffinePointTarget {
                    x: Target::Wire(Wire {
                        gate,
                        input: CurveDblGate::<C, InnerC>::WIRE_X_OLD,
                    }),
                    y: Target::Wire(Wire {
                        gate,
                        input: CurveDblGate::<C, InnerC>::WIRE_Y_OLD,
                    }),
                };
            } else {
                acc = AffinePointTarget {
                    x: Target::Wire(Wire {
                        gate,
                        input: CurveDblGate::<C, InnerC>::WIRE_X_NEW,
                    }),
                    y: Target::Wire(Wire {
                        gate,
                        input: CurveDblGate::<C, InnerC>::WIRE_Y_NEW,
                    }),
                };
            }

            // Also double the filler, so we can subtract out a rescaled version later.
            filler = filler.double();
        }

        // Subtract (a rescaled version of) the arbitrary nonzero value that we started with.
        let filler_target = self.constant_affine_point(filler);
        acc = self.curve_sub::<InnerC>(acc, filler_target);

        // By now we've accumulated all the bits of each scalar, but we also need to accumulate the dibits.
        for j in 0..parts.len() {
            for dibits_chunk in all_dibits[j].chunks(Base4SumGate::<C>::NUM_LIMBS) {
                assert_eq!(dibits_chunk.len(), Base4SumGate::<C>::NUM_LIMBS);

                let gate = self.num_gates();
                self.add_gate_no_constants(Base4SumGate::new(gate));
                self.copy(
                    scalar_acc_unsigned[j],
                    Target::Wire(Wire {
                        gate,
                        input: Base4SumGate::<C>::WIRE_ACC_OLD,
                    }),
                );
                scalar_acc_unsigned[j] = Target::Wire(Wire {
                    gate,
                    input: Base4SumGate::<C>::WIRE_ACC_NEW,
                });

                for (i, &dibit) in dibits_chunk.iter().enumerate() {
                    self.copy(
                        dibit,
                        Target::Wire(Wire {
                            gate,
                            input: Base4SumGate::<C>::WIRE_LIMB_0 + i,
                        }),
                    );
                }
            }
        }

        // Finally, assert that each unsigned accumulator matches the original scalar.
        for (j, part) in parts.iter().enumerate() {
            self.copy(scalar_acc_unsigned[j], part.scalar);
        }

        CurveMsmEndoResult {
            msm_result: acc,
            actual_scalars: scalar_acc_signed,
        }
    }

    /// Adds a gate to the circuit, without doing any routing.
    pub fn add_gate_no_constants<G: Gate<C>>(&mut self, gate: G) {
        self.add_gate(gate, Vec::new());
    }

    /// Adds a gate to the circuit, without doing any routing.
    pub fn add_gate<G: Gate<C>>(&mut self, gate: G, gate_constants: Vec<C::ScalarField>) {
        debug_assert!(G::PREFIX.len() + gate_constants.len() <= NUM_CONSTANTS);

        // Merge the gate type's prefix bits with the given gate config constants.
        let mut all_constants = Vec::new();
        for &prefix_bit in G::PREFIX {
            all_constants.push(C::ScalarField::from_canonical_bool(prefix_bit));
        }
        all_constants.extend(gate_constants);

        // Pad if not all constants were used.
        while all_constants.len() < NUM_CONSTANTS {
            all_constants.push(C::ScalarField::ZERO);
        }

        self.gate_constants.push(all_constants);
        self.add_generator(gate);
        *self.gate_counts.entry(G::NAME).or_insert(0) += 1;
    }

    pub fn add_generator<G: WitnessGenerator<C::ScalarField>>(&mut self, gate: G) {
        self.generators.push(Box::new(gate));
    }

    pub fn num_gates(&self) -> usize {
        self.gate_constants.len()
    }

    /// Add a copy constraint between two routing targets.
    pub fn copy(&mut self, target_1: Target, target_2: Target) {
        self.copy_constraints.push((target_1, target_2));
    }

    /// Add a copy constraint between two affine targets.
    pub fn copy_curve(
        &mut self,
        affine_target_1: AffinePointTarget,
        affine_target_2: AffinePointTarget,
    ) {
        self.copy_constraints
            .push((affine_target_1.x, affine_target_2.x));
        self.copy_constraints
            .push((affine_target_1.y, affine_target_2.y));
    }

    /// Adds a gate with random wire values. By adding `k` of these gates, we can ensure that
    /// nothing is learned by opening the wire polynomials at `k` points outside of H.
    fn add_blinding_gate(&mut self) {
        let gate = self.num_gates();
        self.add_gate_no_constants(BufferGate::new(gate));
        for input in 0..NUM_WIRES {
            self.add_generator(RandomGenerator {
                target: Target::Wire(Wire { gate, input }),
            });
        }

        struct RandomGenerator {
            target: Target,
        }

        impl<F: Field> WitnessGenerator<F> for RandomGenerator {
            fn dependencies(&self) -> Vec<Target> {
                vec![]
            }

            fn generate(
                &self,
                _constants: &Vec<Vec<F>>,
                _witness: &PartialWitness<F>,
            ) -> PartialWitness<F> {
                let mut result = PartialWitness::new();
                result.set_target(self.target, F::rand());
                result
            }
        }
    }

    pub fn build(mut self) -> Circuit<C> {
        // Since we will open each polynomial at three points outside of H, we need three random
        // values to ensure nothing is learned from the out-of-H openings.
        for _i in 0..3 {
            self.add_blinding_gate();
        }

        // Print gate counts.
        println!("Gate counts:");
        for (gate, count) in &self.gate_counts {
            println!("{}: {}", gate, count);
        }
        println!();

        // Pad to a power of two.
        println!("Total gates before padding: {}", self.num_gates());
        while !self.num_gates().is_power_of_two() {
            // Add an empty gate.
            self.add_gate_no_constants(BufferGate::new(self.num_gates()));
        }
        println!("Total gates after padding: {}", self.num_gates());

        let degree = self.num_gates();
        let degree_pow = log2_strict(degree);
        let routing_target_partitions = self.get_routing_partitions();
        let wire_partitions = routing_target_partitions.to_wire_partitions();
        let sigma = wire_partitions.to_sigma();

        let CircuitBuilder {
            security_bits,
            public_input_index: num_public_inputs,
            gate_constants,
            generators,
            ..
        } = self;

        let fft_precomputation_n = fft_precompute(degree);
        let fft_precomputation_8n = fft_precompute(degree * 8);

        let subgroup_generator_n = C::ScalarField::primitive_root_of_unity(degree_pow);
        let subgroup_generator_8n = C::ScalarField::primitive_root_of_unity(degree_pow + 3);
        let subgroup_n = C::ScalarField::cyclic_subgroup_known_order(subgroup_generator_n, degree);
        let subgroup_8n =
            C::ScalarField::cyclic_subgroup_known_order(subgroup_generator_8n, 8 * degree);

        // TODO: Shouldn't this be random?
        let pedersen_g: Vec<_> = (0..degree)
            .map(|i| blake_hash_usize_to_curve::<C>(i))
            .collect();
        let pedersen_h = blake_hash_usize_to_curve::<C>(degree);
        let u = blake_hash_usize_to_curve::<C>(degree + 1);

        let w = 8; // TODO: Should really be set dynamically based on MSM size.
        let pedersen_g_msm_precomputation =
            msm_precompute(&AffinePoint::batch_to_projective(&pedersen_g), w);

        // While gate_constants is indexed by gate index first, this is indexed by wire index first.
        let wire_constants = transpose::<C::ScalarField>(&gate_constants);

        let constants_coeffs: Vec<Vec<C::ScalarField>> =
            values_to_coeffs(&wire_constants, &fft_precomputation_n);
        let constants_8n = coeffs_to_values_padded(&constants_coeffs, &fft_precomputation_8n);
        let c_constants = PolynomialCommitment::coeffs_vec_to_commitments(
            &constants_coeffs,
            &pedersen_g_msm_precomputation,
            pedersen_h,
            false, // Circuit blinding is not necessary here.
        );

        // Convert sigma's values to scalar field elements and split it into degree-n chunks.
        let sigma_chunks = sigma_polynomials(sigma, degree, subgroup_generator_n);

        // Compute S_sigma, then a commitment to it.
        let s_sigma_coeffs = values_to_coeffs(&sigma_chunks, &fft_precomputation_n);
        let s_sigma_values_8n = coeffs_to_values_padded(&s_sigma_coeffs, &fft_precomputation_8n);
        let c_s_sigmas = PolynomialCommitment::coeffs_vec_to_commitments(
            &s_sigma_coeffs,
            &pedersen_g_msm_precomputation,
            pedersen_h,
            false, // Circuit blinding is not necessary here.
        );

        Circuit {
            security_bits,
            num_public_inputs,
            gate_constants,
            routing_target_partitions,
            generators,
            subgroup_generator_n,
            subgroup_generator_8n,
            subgroup_n,
            subgroup_8n,
            pedersen_g,
            pedersen_h,
            u,
            constants_coeffs,
            constants_8n,
            c_constants,
            s_sigma_coeffs,
            s_sigma_values_8n,
            c_s_sigmas,
            pedersen_g_msm_precomputation,
            fft_precomputation_n,
            fft_precomputation_8n,
        }
    }

    fn get_routing_partitions(&self) -> TargetPartitions {
        let mut partitions = TargetPartitions::new();

        for i in 0..self.virtual_target_index {
            partitions.add_partition(Target::VirtualTarget(VirtualTarget { index: i }));
        }

        for gate in 0..self.num_gates() {
            for input in 0..NUM_WIRES {
                partitions.add_partition(Target::Wire(Wire { gate, input }));
            }
        }

        for &(a, b) in &self.copy_constraints {
            partitions.merge(a, b);
        }

        partitions
    }
}
