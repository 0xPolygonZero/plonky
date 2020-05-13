use std::collections::HashMap;
use std::time::Instant;

use crate::{AffinePoint, Curve, Field, generate_rescue_constants, HaloEndomorphismCurve};
use crate::plonk_gates::{ArithmeticGate, BufferGate, CurveAddGate, CurveDblGate, Gate, PublicInputGate, RescueStepAGate, RescueStepBGate, Base4SumGate};
use crate::util::ceil_div_usize;

pub(crate) const NUM_WIRES: usize = 9;
pub(crate) const NUM_ROUTED_WIRES: usize = 6;
pub(crate) const NUM_ADVICE_WIRES: usize = NUM_WIRES - NUM_ROUTED_WIRES;
pub(crate) const NUM_CONSTANTS: usize = 5;
pub(crate) const GRID_WIDTH: usize = 65;
// This is currently dominated by Base4SumGate. It has degree-4n constraints, and its prefix is 4
// bits long, so its filtered constraints are degree-8n. Dividing by Z_H makes t degree-7n.
pub(crate) const QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER: usize = 7;

pub struct PartialWitness<F: Field> {
    wire_values: HashMap<Target, F>,
}

impl<F: Field> PartialWitness<F> {
    pub fn new() -> Self {
        PartialWitness { wire_values: HashMap::new() }
    }

    pub fn get_target(&self, target: Target) -> F {
        self.wire_values[&target]
    }

    pub fn set_target(&mut self, target: Target, value: F) {
        let old_value = self.wire_values.insert(target, value);
        debug_assert!(old_value.is_none(), "Target was set twice");
    }

    pub fn set_wire(&mut self, wire: Wire, value: F) {
        self.set_target(Target::Wire(wire), value);
    }

    pub fn get_wire(&self, wire: Wire) -> F {
        self.get_target(Target::Wire(wire))
    }
}

pub struct Witness<F: Field> {
    wire_values: Vec<Vec<F>>,
}

pub trait WitnessGenerator<F: Field>: 'static {
    fn dependencies(&self) -> Vec<Target>;

    /// Given a partial witness, return any newly generated values. The caller will merge them in.
    fn generate(&self, circuit: Circuit<F>, witness: &PartialWitness<F>) -> PartialWitness<F>;
}

pub struct Circuit<F: Field> {
    pub gate_constants: Vec<Vec<F>>,
    pub routing_target_partitions: RoutingTargetPartitions,
    pub generators: Vec<Box<dyn WitnessGenerator<F>>>,
}

impl<F: Field> Circuit<F> {
    pub fn num_gates(&self) -> usize {
        self.gate_constants.len()
    }

    pub fn generate_witness(&self) {
        let start = Instant::now();
        println!("Witness generation took {}s", start.elapsed().as_secs_f32());
        todo!()
    }
}

/// A sort of proxy wire, in the context of routing and witness generation. It is not an actual
/// witness element (i.e. wire) itself, but it can be copy-constrained to wires, listed as a
/// dependency in generators, etc.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct VirtualTarget {
    pub index: usize,
}

/// Represents a wire in the circuit.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct Wire {
    /// The index of the associated gate.
    pub gate: usize,
    /// The index of the gate input wherein this wire is inserted.
    pub input: usize,
}

/// A routing target.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum Target {
    VirtualTarget(VirtualTarget),
    Wire(Wire),
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct PublicInput {
    pub index: usize,
}

/// See `PublicInputGate` for an explanation of how we make public inputs routable.
impl PublicInput {
    fn original_wire(&self) -> Wire {
        let gate = self.index / NUM_WIRES * 2;
        let input = self.index % NUM_WIRES;
        Wire { gate, input }
    }

    pub fn routable_target(&self) -> Target {
        let Wire { mut gate, mut input } = self.original_wire();
        if input > NUM_ROUTED_WIRES {
            gate += 1;
            input -= NUM_ROUTED_WIRES;
        }
        Target::Wire(Wire { gate, input })
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct AffinePointTarget {
    x: Target,
    y: Target,
}

impl AffinePointTarget {
    pub fn to_vec(&self) -> Vec<Target> {
        vec![self.x, self.y]
    }
}

pub struct CircuitBuilder<F: Field> {
    pub(crate) security_bits: usize,
    public_input_index: usize,
    virtual_target_index: usize,
    gate_counts: HashMap<&'static str, usize>,
    gate_constants: Vec<Vec<F>>,
    copy_constraints: Vec<(Target, Target)>,
    generators: Vec<Box<dyn WitnessGenerator<F>>>,
    constant_wires: HashMap<F, Target>,
}

/// Represents a scalar * point multiplication operation.
pub struct CurveMulOp {
    pub scalar: Target,
    pub point: AffinePointTarget,
}

pub struct CurveMulEndoResult {
    pub mul_result: AffinePointTarget,
    pub actual_scalar: Target,
}

pub struct CurveMsmEndoResult {
    pub msm_result: AffinePointTarget,
    /// While `msm` computes a sum of `[s] P` terms, `msm_endo` computes a sum of `[n(s)] P` terms
    /// for some injective `n`. Here we return each `n(s)`, i.e., the scalar by which the point was
    /// actually multiplied.
    pub actual_scalars: Vec<Target>,
}

impl<F: Field> CircuitBuilder<F> {
    pub fn new(security_bits: usize) -> Self {
        CircuitBuilder {
            security_bits,
            public_input_index: 0,
            virtual_target_index: 0,
            gate_counts: HashMap::new(),
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
        (0..n).map(|i| self.stage_public_input()).collect()
    }

    /// Add `PublicInputGate`s which enable public inputs to be routed. Should be called after all
    /// `stage_public_input[s]` calls, but before any gates are added.
    pub fn route_public_inputs(&mut self) {
        debug_assert_eq!(self.num_gates(), 0, "Must be called before any gates are added");
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
        (0..n).map(|i| self.add_virtual_target()).collect()
    }

    pub fn add_virtual_point_target(&mut self) -> AffinePointTarget {
        let x = self.add_virtual_target();
        let y = self.add_virtual_target();
        AffinePointTarget { x, y }
    }

    pub fn add_virtual_point_targets(&mut self, n: usize) -> Vec<AffinePointTarget> {
        (0..n).map(|i| self.add_virtual_point_target()).collect()
    }

    pub fn zero_wire(&mut self) -> Target {
        self.constant_wire(F::ZERO)
    }

    pub fn one_wire(&mut self) -> Target {
        self.constant_wire(F::ONE)
    }

    pub fn neg_one_wire(&mut self) -> Target {
        self.constant_wire(F::NEG_ONE)
    }

    pub fn constant_wire(&mut self, c: F) -> Target {
        if self.constant_wires.contains_key(&c) {
            self.constant_wires[&c]
        } else {
            let result = self.create_constant_wire(c);
            self.constant_wires.insert(c, result);
            result
        }
    }

    pub fn constant_wire_u32(&mut self, c: u32) -> Target {
        self.constant_wire(F::from_canonical_u32(c))
    }

    fn create_constant_wire(&mut self, c: F) -> Target {
        let index = self.num_gates();
        self.add_gate(BufferGate::new(index), vec![c]);
        Target::Wire(Wire { gate: index, input: BufferGate::WIRE_BUFFER_CONST })
    }

    pub fn constant_affine_point<C: Curve<BaseField=F>>(&mut self, point: AffinePoint<C>) -> AffinePointTarget {
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
        let zero = self.zero_wire();
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
        self.add_gate(ArithmeticGate::new(index), vec![F::ONE, F::ONE, F::ZERO]);
        self.copy(x, Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_MULTIPLICAND_0 }));
        self.copy(one, Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_MULTIPLICAND_1 }));
        self.copy(y, Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_ADDEND }));
        Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_OUTPUT })
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
        self.add_gate(ArithmeticGate::new(index), vec![F::ONE, F::NEG_ONE, F::ZERO]);
        self.copy(x, Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_MULTIPLICAND_0 }));
        self.copy(one, Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_MULTIPLICAND_1 }));
        self.copy(y, Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_ADDEND }));
        Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_OUTPUT })
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
        self.add_gate(ArithmeticGate::new(index), vec![F::ONE, F::ZERO, F::ZERO]);
        self.copy(x, Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_MULTIPLICAND_0 }));
        self.copy(y, Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_MULTIPLICAND_1 }));
        self.copy(zero, Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_ADDEND }));
        Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_OUTPUT })
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

            fn generate(&self, circuit: Circuit<F>, witness: &PartialWitness<F>) -> PartialWitness<F> {
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
        let f_bits = F::BITS - 1;
        assert_eq!(f_bits, 254, "We currently only handle fields of size 2^254 + epsilon");
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
        for chunk in dibits.chunks(Base4SumGate::NUM_LIMBS) {
            assert_eq!(chunk.len(), Base4SumGate::NUM_LIMBS, "Should not have a partial chunk");

            let index = self.num_gates();
            self.add_gate_no_constants(Base4SumGate::new(index));
            for i in 0..chunk.len() {
                self.copy(sum, Target::Wire(Wire { gate: index, input: Base4SumGate::WIRE_ACC_OLD }));
                self.copy(chunk[i], Target::Wire(Wire { gate: index, input: Base4SumGate::WIRE_LIMB_0 + i }));
                sum = Target::Wire(Wire { gate: index, input: Base4SumGate::WIRE_ACC_NEW })
            }
        }
        self.copy(sum, x);

        x_sqrt
    }

    /// Compute `x^power`, where `power` is a constant.
    pub fn exp_constant(&mut self, x: Target, power: F) -> Target {
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
        self.exp_constant(x, F::from_canonical_usize(power))
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

            fn generate(&self, circuit: Circuit<F>, witness: &PartialWitness<F>) -> PartialWitness<F> {
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
        self.add_gate(ArithmeticGate::new(index), vec![F::ONE, F::ONE, F::ZERO]);
        self.copy(x, Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_MULTIPLICAND_0 }));
        self.copy(y, Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_MULTIPLICAND_1 }));
        self.copy(z, Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_ADDEND }));
        Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_OUTPUT })
    }

    /// Multiply and subtract; i.e. computes `x * y - z`.
    pub fn mul_sub(&mut self, x: Target, y: Target, z: Target) -> Target {
        let index = self.num_gates();
        self.add_gate(ArithmeticGate::new(index), vec![F::ONE, F::NEG_ONE, F::ZERO]);
        self.copy(x, Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_MULTIPLICAND_0 }));
        self.copy(y, Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_MULTIPLICAND_1 }));
        self.copy(z, Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_ADDEND }));
        Target::Wire(Wire { gate: index, input: ArithmeticGate::<F>::WIRE_OUTPUT })
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

            fn generate(&self, circuit: Circuit<F>, witness: &PartialWitness<F>) -> PartialWitness<F> {
                let x = witness.wire_values[&self.x];
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
        let generator = SplitGenerator { x, bits: bits.clone(), dibits: dibits.clone() };
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

    pub fn rescue_sponge(
        &mut self,
        inputs: &[Target],
        num_outputs: usize,
    ) -> Vec<Target> {
        // This is a r=2, c=1 sponge function with a single absorption and a single squeeze.
        let zero = self.zero_wire();
        let mut state = [zero, zero, zero];

        // Absorb all input chunks.
        for input_chunk in inputs.chunks(2) {
            for i in 0..input_chunk.len() {
                state[i] = self.add(state[i], input_chunk[i]);
            }
            state = self.rescue_permutation_3x3(state);
        }

        // Squeeze until we have the desired number of outputs.
        let mut outputs = Vec::new();
        while outputs.len() < num_outputs {
            outputs.push(state[0]);
            if outputs.len() < num_outputs {
                outputs.push(state[1]);
            }
            if outputs.len() < num_outputs {
                state = self.rescue_permutation_3x3(state);
            }
        }

        outputs
    }

    pub fn rescue_permutation_3x3(&mut self, inputs: [Target; 3]) -> [Target; 3] {
        let all_constants = generate_rescue_constants(3, self.security_bits);
        let mut state = inputs;
        for (a_constants, b_constants) in all_constants.into_iter() {
            state = self.rescue_round(state, a_constants, b_constants);
        }
        state
    }

    fn rescue_round(
        &mut self,
        inputs: [Target; 3],
        a_constants: Vec<F>,
        b_constants: Vec<F>,
    ) -> [Target; 3] {
        let a_index = self.num_gates();
        let a_gate = RescueStepAGate::new(a_index);
        self.add_gate(a_gate, a_constants);

        let b_index = self.num_gates();
        let b_gate = RescueStepBGate::new(b_index);
        self.add_gate(b_gate, b_constants);

        let a_in_0_target = Target::Wire(Wire { gate: a_index, input: RescueStepAGate::<F>::WIRE_INPUT_0 });
        let a_in_1_target = Target::Wire(Wire { gate: a_index, input: RescueStepAGate::<F>::WIRE_INPUT_1 });
        let a_in_2_target = Target::Wire(Wire { gate: a_index, input: RescueStepAGate::<F>::WIRE_INPUT_2 });
        let a_out_0_target = Target::Wire(Wire { gate: a_index, input: RescueStepAGate::<F>::WIRE_OUTPUT_0 });
        let a_out_1_target = Target::Wire(Wire { gate: a_index, input: RescueStepAGate::<F>::WIRE_OUTPUT_1 });
        let a_out_2_target = Target::Wire(Wire { gate: a_index, input: RescueStepAGate::<F>::WIRE_OUTPUT_2 });

        let b_in_0_target = Target::Wire(Wire { gate: b_index, input: RescueStepBGate::<F>::WIRE_INPUT_0 });
        let b_in_1_target = Target::Wire(Wire { gate: b_index, input: RescueStepBGate::<F>::WIRE_INPUT_1 });
        let b_in_2_target = Target::Wire(Wire { gate: b_index, input: RescueStepBGate::<F>::WIRE_INPUT_2 });
        let b_out_0_target = Target::Wire(Wire { gate: b_index, input: RescueStepBGate::<F>::WIRE_OUTPUT_0 });
        let b_out_1_target = Target::Wire(Wire { gate: b_index, input: RescueStepBGate::<F>::WIRE_OUTPUT_1 });
        let b_out_2_target = Target::Wire(Wire { gate: b_index, input: RescueStepBGate::<F>::WIRE_OUTPUT_2 });

        self.copy(inputs[0], a_in_0_target);
        self.copy(inputs[1], a_in_1_target);
        self.copy(inputs[2], a_in_2_target);
        self.copy(a_out_0_target, b_in_0_target);
        self.copy(a_out_1_target, b_in_1_target);
        self.copy(a_out_2_target, b_in_2_target);

        [b_out_0_target, b_out_1_target, b_out_2_target]
    }

    /// Assert that a given coordinate pair is on the curve `C`.
    pub fn curve_assert_valid<C: Curve<BaseField=F>>(&mut self, p: AffinePointTarget) {
        // TODO
    }

    pub fn curve_neg<C: Curve<BaseField=F>>(
        &mut self,
        p: AffinePointTarget,
    ) -> AffinePointTarget {
        let neg_y = self.neg(p.y);
        AffinePointTarget { x: p.x, y: neg_y }
    }

    pub fn curve_add<C: Curve<BaseField=F>>(
        &mut self,
        p_1: AffinePointTarget,
        p_2: AffinePointTarget,
    ) -> AffinePointTarget {
        let add_index = self.num_gates();
        self.add_gate_no_constants(CurveAddGate::<C>::new(add_index));
        let buffer_index = self.num_gates();
        self.add_gate_no_constants(BufferGate::new(buffer_index));

        // TODO: Wiring.

        let result_x = Target::Wire(Wire { gate: buffer_index, input: CurveAddGate::<C>::WIRE_GROUP_ACC_X });
        let result_y = Target::Wire(Wire { gate: buffer_index, input: CurveAddGate::<C>::WIRE_GROUP_ACC_Y });
        AffinePointTarget { x: result_x, y: result_y }
    }


    pub fn curve_double<C: Curve<BaseField=F>>(
        &mut self,
        p: AffinePointTarget,
    ) -> AffinePointTarget {
        let idx_dbl = self.num_gates();
        self.add_gate_no_constants(CurveDblGate::<C>::new(idx_dbl));
        self.copy(p.x, Target::Wire(Wire { gate: idx_dbl, input: CurveDblGate::<C>::WIRE_X_OLD }));
        self.copy(p.y, Target::Wire(Wire { gate: idx_dbl, input: CurveDblGate::<C>::WIRE_Y_OLD }));
        AffinePointTarget {
            x: Target::Wire(Wire { gate: idx_dbl, input: CurveDblGate::<C>::WIRE_X_NEW }),
            y: Target::Wire(Wire { gate: idx_dbl, input: CurveDblGate::<C>::WIRE_Y_NEW }),
        }
    }

    pub fn curve_sub<C: Curve<BaseField=F>>(
        &mut self,
        p_1: AffinePointTarget,
        p_2: AffinePointTarget,
    ) -> AffinePointTarget {
        let neg_p_2 = self.curve_neg::<C>(p_2);
        self.curve_add::<C>(p_1, neg_p_2)
    }

    pub fn curve_mul<C: Curve<BaseField=F>>(&mut self, mul: CurveMulOp) -> AffinePointTarget {
        self.curve_msm::<C>(&[mul])
    }

    pub fn curve_mul_endo<C: HaloEndomorphismCurve<BaseField=F>>(
        &mut self,
        mul: CurveMulOp,
    ) -> CurveMulEndoResult {
        let result = self.curve_msm_endo::<C>(&[mul]);
        CurveMulEndoResult {
            mul_result: result.msm_result,
            actual_scalar: result.actual_scalars[0],
        }
    }

    /// Note: This assumes the most significant bit of each scalar is unset. This occurs with high
    /// probability if the field size is slightly larger than a power of two and the scalars are
    /// uniformly random.
    pub fn curve_msm<C: Curve<BaseField=F>>(
        &mut self,
        parts: &[CurveMulOp],
    ) -> AffinePointTarget {
        // We assume each most significant bit is unset; see the note in the method doc.
        let f_bits = F::BITS - 1;

        let all_bits: Vec<Vec<Target>> = parts.iter()
            .map(|part| self.split_binary(part.scalar, f_bits))
            .collect();

        // Normally we would start with zero, but to avoid exceptional cases, we start with some
        // arbitrary nonzero point and subtract it later. This avoids exception with high
        // probability provided that the scalars and points are random. (We don't worry about
        // non-random inputs from malicious provers, since our curve gates will be unsatisfiable in
        // exceptional cases.)
        let mut filler = C::GENERATOR_AFFINE;
        let mut acc = self.constant_affine_point(filler);
        let mut scalar_accs = vec![self.zero_wire(); parts.len()];

        for i in (0..f_bits).rev() {
            // Route the accumulator to the first curve addition gate's inputs.
            self.copy(acc.x, Target::Wire(
                Wire { gate: self.num_gates(), input: CurveAddGate::<C>::WIRE_GROUP_ACC_X }));
            self.copy(acc.y, Target::Wire(
                Wire { gate: self.num_gates(), input: CurveAddGate::<C>::WIRE_GROUP_ACC_Y }));

            for (j, part) in parts.iter().enumerate() {
                let bit = all_bits[j][i];

                let idx_add = self.num_gates();
                self.add_gate_no_constants(CurveAddGate::<C>::new(idx_add));
                self.copy(scalar_accs[j], Target::Wire(
                    Wire { gate: idx_add, input: CurveAddGate::<C>::WIRE_SCALAR_ACC_OLD }));
                scalar_accs[j] = Target::Wire(
                    Wire { gate: idx_add, input: CurveAddGate::<C>::WIRE_SCALAR_ACC_NEW });
                self.copy(part.point.x, Target::Wire(
                    Wire { gate: idx_add, input: CurveAddGate::<C>::WIRE_ADDEND_X }));
                self.copy(part.point.y, Target::Wire(
                    Wire { gate: idx_add, input: CurveAddGate::<C>::WIRE_ADDEND_Y }));
                self.copy(bit, Target::Wire(
                    Wire { gate: idx_add, input: CurveAddGate::<C>::WIRE_SCALAR_BIT }));
            }

            // Double the accumulator.
            let idx_dbl = self.num_gates();
            self.add_gate_no_constants(CurveDblGate::<C>::new(idx_dbl));
            // No need to route the double gate's inputs, because the last add gate would have
            // constrained them.
            // Normally, we will take the double gate's outputs as the new accumulator. If we just
            // completed the last iteration of the MSM though, then we don't want to perform a final
            // doubling, so we will take its inputs as the result instead.
            if i == f_bits - 1 {
                acc = AffinePointTarget {
                    x: Target::Wire(Wire { gate: idx_dbl, input: CurveDblGate::<C>::WIRE_X_OLD }),
                    y: Target::Wire(Wire { gate: idx_dbl, input: CurveDblGate::<C>::WIRE_Y_OLD }),
                };
            } else {
                acc = AffinePointTarget {
                    x: Target::Wire(Wire { gate: idx_dbl, input: CurveDblGate::<C>::WIRE_X_NEW }),
                    y: Target::Wire(Wire { gate: idx_dbl, input: CurveDblGate::<C>::WIRE_Y_NEW }),
                };
            }

            // Also double the filler, so we can subtract out a rescaled version later.
            filler = filler.double();
        }

        // Subtract (a rescaled version of) the arbitrary nonzero value that we started with.
        let filler_target = self.constant_affine_point(filler);
        acc = self.curve_sub::<C>(acc, filler_target);

        // Assert that each accumulation of scalar bits matches the original scalar.
        for (j, part) in parts.iter().enumerate() {
            self.copy(scalar_accs[j], part.scalar);
        }

        acc
    }

    /// Like `curve_msm`, but uses the endomorphism described in the Halo paper.
    pub fn curve_msm_endo<C: HaloEndomorphismCurve<BaseField=F>>(
        &mut self,
        parts: &[CurveMulOp],
    ) -> CurveMsmEndoResult {
        todo!()
    }

    /// Adds a gate to the circuit, without doing any routing.
    fn add_gate_no_constants<G: Gate<F>>(&mut self, gate: G) {
        self.add_gate(gate, Vec::new());
    }

    /// Adds a gate to the circuit, without doing any routing.
    pub fn add_gate<G: Gate<F>>(&mut self, gate: G, gate_constants: Vec<F>) {
        // Merge the gate type's prefix bits with the given gate config constants.
        debug_assert!(G::PREFIX.len() + gate_constants.len() <= NUM_CONSTANTS);
        let mut all_constants = Vec::new();
        for &prefix_bit in G::PREFIX {
            all_constants.push(F::from_canonical_bool(prefix_bit));
        }
        all_constants.extend(gate_constants);
        self.gate_constants.push(all_constants);
        self.add_generator(gate);
        *self.gate_counts.entry(G::NAME).or_insert(0) += 1;
    }

    pub fn add_generator<G: WitnessGenerator<F>>(&mut self, gate: G) {
        self.generators.push(Box::new(gate));
    }

    pub fn num_gates(&self) -> usize {
        self.gate_constants.len()
    }

    /// Add a copy constraint between two routing targets.
    pub fn copy(&mut self, target_1: Target, target_2: Target) {
        self.copy_constraints.push((target_1, target_2));
    }

    pub fn build(self) -> Circuit<F> {
        let routing_target_partitions = self.get_routing_partitions();
        let CircuitBuilder { gate_counts, gate_constants, generators, .. } = self;

        // Print gate counts.
        println!("Gate counts:");
        for (gate, count) in gate_counts {
            println!("{}: {}", gate, count);
        }
        println!();

        Circuit { gate_constants, routing_target_partitions, generators }
    }

    fn get_routing_partitions(&self) -> RoutingTargetPartitions {
        let mut partitions = RoutingTargetPartitions::new();

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

pub struct RoutingTargetPartitions {
    partitions: Vec<Vec<Target>>,
    indices: HashMap<Target, usize>,
}

impl RoutingTargetPartitions {
    fn new() -> Self {
        Self { partitions: Vec::new(), indices: HashMap::new() }
    }

    /// Add a new partition with a single member.
    fn add_partition(&mut self, target: Target) {
        let index = self.partitions.len();
        self.partitions.push(vec![target]);
        self.indices.insert(target, index);
    }

    /// Merge the two partitions containing the two given targets. Does nothing if the targets are
    /// already members of the same partition.
    fn merge(&mut self, a: Target, b: Target) {
        let a_index = self.indices[&a];
        let b_index = self.indices[&b];
        if a_index != b_index {
            // Merge a's partition into b's partition, leaving a's partition empty.
            // We have to clone because Rust's borrow checker doesn't know that
            // self.partitions[b_index] and self.partitions[b_index] are disjoint.
            let mut a_partition = self.partitions[a_index].clone();
            let b_partition = &mut self.partitions[b_index];
            for a_sibling in &a_partition {
                *self.indices.get_mut(a_sibling).unwrap() = b_index;
            }
            b_partition.append(&mut a_partition);
        }
    }

    fn to_gate_inputs(&self) -> GateInputPartitions {
        // Here we just drop all CircuitInputs, leaving all GateInputs.
        let mut partitions = Vec::new();
        let mut indices = HashMap::new();

        for old_partition in &self.partitions {
            let mut new_partition = Vec::new();
            for target in old_partition {
                if let Target::Wire(gi) = *target {
                    new_partition.push(gi);
                }
            }
            partitions.push(new_partition);
        }

        for (&target, &index) in &self.indices {
            if let Target::Wire(gi) = target {
                indices.insert(gi, index);
            }
        }

        GateInputPartitions { partitions, indices }
    }
}

struct GateInputPartitions {
    partitions: Vec<Vec<Wire>>,
    indices: HashMap<Wire, usize>,
}
