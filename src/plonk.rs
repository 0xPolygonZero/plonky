use std::collections::HashMap;
use std::time::Instant;

use crate::{Field, generate_rescue_constants, Curve, HaloEndomorphismCurve};
use crate::plonk_gates::{BufferGate, Gate, RescueStepAGate, CurveAddGate, MaddGate};
use std::marker::PhantomData;

pub(crate) const NUM_WIRES: usize = 9;
pub(crate) const NUM_CONSTANTS: usize = 5;
pub(crate) const GRID_WIDTH: usize = 65;
pub(crate) const QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER: usize = 7;

pub struct PartialWitness<F: Field> {
    pub wire_values: HashMap<GateInput, F>,
}

impl<F: Field> PartialWitness<F> {
    pub fn new() -> Self {
        PartialWitness { wire_values: HashMap::new() }
    }
}

pub struct Witness<F: Field> {
    wire_values: Vec<Vec<F>>,
}

pub trait WitnessGenerator<F: Field> {
    fn dependencies(&self) -> Vec<GateInput>;

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

#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct CircuitInput {
    pub index: usize,
}

#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct GateInput {
    pub gate: usize,
    pub input: usize,
}

#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum RoutingTarget {
    CircuitInput(CircuitInput),
    GateInput(GateInput),
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct AffinePointTarget {
    x: RoutingTarget,
    y: RoutingTarget,
}

pub struct CircuitBuilder<F: Field> {
    circuit_input_index: usize,
    gate_constants: Vec<Vec<F>>,
    copy_constraints: Vec<(RoutingTarget, RoutingTarget)>,
    generators: Vec<Box<dyn WitnessGenerator<F>>>,
    constant_wires: HashMap<F, RoutingTarget>,
}

/// A component of an MSM, or in other words, an individual scalar-group multiplication.
pub struct MsmPart {
    pub scalar_bits: Vec<RoutingTarget>,
    pub x: RoutingTarget,
    pub y: RoutingTarget,
}

pub struct MsmResult {
    pub x: RoutingTarget,
    pub y: RoutingTarget,
    /// For each part, we return the weighted sum of the given scalar bits.
    pub scalars: Vec<RoutingTarget>,
}

pub struct MsmEndoResult {
    pub msm_result: MsmResult,
    /// While `msm` computes a sum of `[s] P` terms, `msm_end` computes a sum of `[n(s)] P` terms
    /// for some injective `n`. Here we return each `n(s)`, i.e., the scalar by which the point was
    /// actually multiplied.
    pub actual_scalars: Vec<RoutingTarget>,
}

impl<F: Field> CircuitBuilder<F> {
    pub fn new() -> Self {
        CircuitBuilder {
            circuit_input_index: 0,
            gate_constants: Vec::new(),
            copy_constraints: Vec::new(),
            generators: Vec::new(),
            constant_wires: HashMap::new(),
        }
    }

    pub fn add_public_input(&mut self) -> GateInput {
        let index = self.num_gates();
        self.add_gate_no_constants(BufferGate::new(index));
        GateInput { gate: index, input: BufferGate::WIRE_BUFFER_PI }
    }

    pub fn add_public_inputs(&mut self, n: usize) -> Vec<GateInput> {
        (0..n).map(|i| self.add_public_input()).collect()
    }

    pub fn add_circuit_input(&mut self) -> CircuitInput {
        let index = self.circuit_input_index;
        self.circuit_input_index += 1;
        CircuitInput { index }
    }

    pub fn zero_wire(&mut self) -> RoutingTarget {
        self.constant_wire(F::ZERO)
    }

    pub fn one_wire(&mut self) -> RoutingTarget {
        self.constant_wire(F::ONE)
    }

    pub fn neg_one_wire(&mut self) -> RoutingTarget {
        self.constant_wire(F::NEG_ONE)
    }

    pub fn constant_wire(&mut self, c: F) -> RoutingTarget {
        if self.constant_wires.contains_key(&c) {
            self.constant_wires[&c]
        } else {
            let result = self.create_constant_wire(c);
            self.constant_wires.insert(c, result);
            result
        }
    }

    pub fn constant_wire_u32(&mut self, c: u32) -> RoutingTarget {
        self.constant_wire(F::from_canonical_u32(c))
    }

    fn create_constant_wire(&mut self, c: F) -> RoutingTarget {
        let index = self.num_gates();
        self.add_gate(BufferGate::new(index), vec![c]);
        RoutingTarget::GateInput(GateInput { gate: index, input: BufferGate::WIRE_BUFFER_CONST })
    }

    pub fn add(&mut self, x: RoutingTarget, y: RoutingTarget) -> RoutingTarget {
        let one = self.one_wire();
        let index = self.num_gates();
        self.add_gate(MaddGate::new(index), vec![F::ONE, F::ONE, F::ZERO]);
        self.copy(x, RoutingTarget::GateInput(GateInput { gate: index, input: MaddGate::<F>::WIRE_MULTIPLICAND_0 }));
        self.copy(one, RoutingTarget::GateInput(GateInput { gate: index, input: MaddGate::<F>::WIRE_MULTIPLICAND_1 }));
        self.copy(y, RoutingTarget::GateInput(GateInput { gate: index, input: MaddGate::<F>::WIRE_ADDEND }));
        RoutingTarget::GateInput(GateInput { gate: index, input: MaddGate::<F>::WIRE_OUTPUT })
    }

    pub fn mul(&mut self, x: RoutingTarget, y: RoutingTarget) -> RoutingTarget {
        let zero = self.zero_wire();
        let index = self.num_gates();
        self.add_gate(MaddGate::new(index), vec![F::ONE, F::ZERO, F::ZERO]);
        self.copy(x, RoutingTarget::GateInput(GateInput { gate: index, input: MaddGate::<F>::WIRE_MULTIPLICAND_0 }));
        self.copy(y, RoutingTarget::GateInput(GateInput { gate: index, input: MaddGate::<F>::WIRE_MULTIPLICAND_1 }));
        self.copy(zero, RoutingTarget::GateInput(GateInput { gate: index, input: MaddGate::<F>::WIRE_ADDEND }));
        RoutingTarget::GateInput(GateInput { gate: index, input: MaddGate::<F>::WIRE_OUTPUT })
    }

    pub fn neg(&mut self, x: RoutingTarget) -> RoutingTarget {
        let neg_one = self.neg_one_wire();
        self.mul(x, neg_one)
    }

    pub fn rescue_hash_2x1(&mut self, inputs: [RoutingTarget; 2]) -> RoutingTarget {
        self.rescue_hash_2x2(inputs)[0]
    }

    /// TODO: Instead define a variable-output hash.
    pub fn rescue_hash_2x2(&mut self, inputs: [RoutingTarget; 2]) -> [RoutingTarget; 2] {
        // This is a width-3 sponge function with a single absorption and a single squeeze.
        let zero = self.zero_wire();
        let outputs = self.rescue_permutation_3x3([inputs[0], inputs[1], zero]);
        [outputs[0], outputs[1]]
    }

    pub fn rescue_permutation_3x3(&mut self, inputs: [RoutingTarget; 3]) -> [RoutingTarget; 3] {
        let all_constants = generate_rescue_constants(3);

        let first_gate_index = self.num_gates();
        for constants in all_constants.into_iter() {
            let gate = RescueStepAGate::new(self.num_gates());
            self.add_gate(gate, constants);
        }
        let last_gate_index = self.num_gates() - 1;

        let in_0_target = RoutingTarget::GateInput(GateInput { gate: first_gate_index, input: RescueStepAGate::<F>::WIRE_INPUT_0 });
        let in_1_target = RoutingTarget::GateInput(GateInput { gate: first_gate_index, input: RescueStepAGate::<F>::WIRE_INPUT_1 });
        let in_2_target = RoutingTarget::GateInput(GateInput { gate: first_gate_index, input: RescueStepAGate::<F>::WIRE_INPUT_1 });

        let out_0_target = RoutingTarget::GateInput(GateInput { gate: last_gate_index, input: RescueStepAGate::<F>::WIRE_OUTPUT_0 });
        let out_1_target = RoutingTarget::GateInput(GateInput { gate: last_gate_index, input: RescueStepAGate::<F>::WIRE_OUTPUT_1 });
        let out_2_target = RoutingTarget::GateInput(GateInput { gate: last_gate_index, input: RescueStepAGate::<F>::WIRE_OUTPUT_1 });

        self.copy(inputs[0], in_0_target);
        self.copy(inputs[1], in_1_target);
        self.copy(inputs[2], in_2_target);

        [out_0_target, out_1_target, out_2_target]
    }

    pub fn curve_neg<C: Curve<BaseField = F>>(
        &mut self,
        p: AffinePointTarget,
    ) -> AffinePointTarget {
        let neg_y = self.neg(p.y);
        AffinePointTarget { x: p.x, y: neg_y }
    }

    pub fn curve_add<C: Curve<BaseField = F>>(
        &mut self,
        p_1: AffinePointTarget,
        p_2: AffinePointTarget,
    ) -> AffinePointTarget {
        let add_index = self.num_gates();
        self.add_gate_no_constants(CurveAddGate::<C>::new(add_index));
        let buffer_index = self.num_gates();
        self.add_gate_no_constants(BufferGate::new(buffer_index));

        // TODO: Wiring.

        let result_x = RoutingTarget::GateInput(GateInput { gate: buffer_index, input: CurveAddGate::<C>::WIRE_GROUP_ACC_X });
        let result_y = RoutingTarget::GateInput(GateInput { gate: buffer_index, input: CurveAddGate::<C>::WIRE_GROUP_ACC_Y });
        AffinePointTarget { x: result_x, y: result_y }
    }

    pub fn curve_sub<C: Curve<BaseField = F>>(
        &mut self,
        p_1: AffinePointTarget,
        p_2: AffinePointTarget,
    ) -> AffinePointTarget {
        let neg_p_2 = self.curve_neg::<C>(p_2);
        self.curve_add::<C>(p_1, neg_p_2)
    }

    pub fn curve_msm<C: Curve<BaseField = F>>(&mut self, parts: &[MsmPart]) -> MsmResult {
        // Normally we would start with zero, but to avoid exceptional cases, we start with some
        // arbitrary nonzero point and subtract it later. This avoids exception with high
        // probability provided that the scalars and points are random. (We don't worry about
        // non-random inputs from malicious provers, since our curve gates will be unsatisfiable in
        // exceptional cases.)
        let mut acc_x = self.constant_wire(C::GENERATOR_AFFINE.x);
        let mut acc_y = self.constant_wire(C::GENERATOR_AFFINE.y);
        let mut scalars = vec![self.zero_wire(); parts.len()];

        let max_bits = parts.iter().map(|p| p.scalar_bits.len()).max().expect("Empty MSM");
        for i in 0..max_bits {
            for part in parts {
                if part.scalar_bits.len() > i {
                    // TODO: Wiring.
                    self.add_gate_no_constants(CurveAddGate::<C>::new(self.num_gates()));
                }
            }
        }

        // Subtract the arbitrary nonzero value that we started with.
        todo!();

        MsmResult {
            x: acc_x,
            y: acc_y,
            scalars,
        }
    }

    /// Like `add_msm`, but uses the endomorphism described in the Halo paper.
    pub fn curve_msm_endo<C: HaloEndomorphismCurve<BaseField = F>>(&mut self, parts: &[MsmPart]) -> MsmEndoResult {
        todo!()
    }

    /// Adds a gate to the circuit, without doing any routing.
    fn add_gate_no_constants<G: Gate<F>>(&mut self, gate: G) {
        self.add_gate(gate, Vec::new());
    }

    /// Adds a gate to the circuit, without doing any routing.
    fn add_gate<G: Gate<F>>(&mut self, gate: G, gate_constants: Vec<F>) {
        // Merge the gate type's prefix bits with the given gate config constants.
        debug_assert!(G::PREFIX.len() + gate_constants.len() <= NUM_CONSTANTS);
        let mut all_constants = Vec::new();
        for &prefix_bit in G::PREFIX {
            all_constants.push(if prefix_bit { F::ONE } else { F::ZERO });
        }
        all_constants.extend(gate_constants);
        self.gate_constants.push(all_constants);
        self.generators.push(Box::new(gate));
    }

    fn num_gates(&self) -> usize {
        self.gate_constants.len()
    }

    /// Add a copy constraint between two routing targets.
    pub fn copy(&mut self, target_1: RoutingTarget, target_2: RoutingTarget) {
        self.copy_constraints.push((target_1, target_2));
    }

    pub fn build(self) -> Circuit<F> {
        let routing_target_partitions = self.get_routing_partitions();
        let CircuitBuilder { circuit_input_index, gate_constants, copy_constraints, generators, constant_wires } = self;
        Circuit { gate_constants, routing_target_partitions, generators }
    }

    fn get_routing_partitions(&self) -> RoutingTargetPartitions {
        let mut partitions = RoutingTargetPartitions::new();

        for i in 0..self.circuit_input_index {
            partitions.add_partition(RoutingTarget::CircuitInput(CircuitInput { index: i }));
        }

        for gate in 0..self.num_gates() {
            for input in 0..NUM_WIRES {
                partitions.add_partition(RoutingTarget::GateInput(GateInput { gate, input }));
            }
        }

        for &(a, b) in &self.copy_constraints {
            partitions.merge(a, b);
        }

        partitions
    }
}

pub struct RoutingTargetPartitions {
    partitions: Vec<Vec<RoutingTarget>>,
    indices: HashMap<RoutingTarget, usize>,
}

impl RoutingTargetPartitions {
    fn new() -> Self {
        Self { partitions: Vec::new(), indices: HashMap::new() }
    }

    /// Add a new partition with a single member.
    fn add_partition(&mut self, target: RoutingTarget) {
        let index = self.partitions.len();
        self.partitions.push(vec![target]);
        self.indices.insert(target, index);
    }

    /// Merge the two partitions containing the two given targets. Does nothing if the targets are
    /// already members of the same partition.
    fn merge(&mut self, a: RoutingTarget, b: RoutingTarget) {
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
                if let &RoutingTarget::GateInput(gi) = target {
                    new_partition.push(gi);
                }
            }
            partitions.push(new_partition);
        }

        for (&target, &index) in &self.indices {
            if let RoutingTarget::GateInput(gi) = target {
                indices.insert(gi, index);
            }
        }

        GateInputPartitions { partitions, indices }
    }
}

struct GateInputPartitions {
    partitions: Vec<Vec<GateInput>>,
    indices: HashMap<GateInput, usize>,
}
