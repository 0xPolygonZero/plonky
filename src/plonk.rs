use std::collections::HashMap;
use std::time::Instant;

use crate::Field;
use crate::plonk_gates::{BufferGate, Gate};

pub(crate) const NUM_WIRES: usize = 9;
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
    fn generate(&self, witness: &PartialWitness<F>) -> PartialWitness<F>;
}

pub struct Circuit<F: Field> {
    gate_ids: Vec<usize>,
    routing_target_partitions: RoutingTargetPartitions,
    generators: Vec<Box<dyn WitnessGenerator<F>>>,
}

impl<F: Field> Circuit<F> {
    pub fn num_gates(&self) -> usize {
        self.gate_ids.len()
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

pub struct CircuitBuilder<F: Field> {
    circuit_input_index: usize,
    gate_ids: Vec<usize>,
    copy_constraints: Vec<(RoutingTarget, RoutingTarget)>,
    generators: Vec<Box<dyn WitnessGenerator<F>>>,
}

pub struct MsmEndoPart {
    scalar: RoutingTarget,
    x: RoutingTarget,
    y: RoutingTarget,
    truncate_to_128: bool,
}

impl<F: Field> CircuitBuilder<F> {
    pub fn new() -> Self {
        CircuitBuilder {
            circuit_input_index: 0,
            gate_ids: Vec::new(),
            copy_constraints: Vec::new(),
            generators: Vec::new(),
        }
    }

    pub fn add_public_input(&mut self) -> GateInput {
        let index = self.gate_ids.len();
        self.add_gate(BufferGate { index });
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

    pub fn add_rescue(&mut self) {
        for r in 0..16 {
            // TODO
        }
    }

    pub fn add_msm_endo(&mut self, parts: &[MsmEndoPart]) {
        todo!();
    }

    /// Adds a gate to the circuit, without doing any routing.
    fn add_gate<G: Gate<F>>(&mut self, gate: G) {
        self.gate_ids.push(G::ID);
        self.generators.push(Box::new(gate));
    }

    /// Add a copy constraint between two routing targets.
    pub fn copy(&mut self, target_1: RoutingTarget, target_2: RoutingTarget) {
        self.copy_constraints.push((target_1, target_2));
    }

    pub fn build(self) -> Circuit<F> {
        let routing_target_partitions = self.get_routing_partitions();
        let CircuitBuilder { circuit_input_index, gate_ids, copy_constraints, generators } = self;
        Circuit { gate_ids, routing_target_partitions, generators }
    }

    fn get_routing_partitions(&self) -> RoutingTargetPartitions {
        let mut partitions = RoutingTargetPartitions::new();

        for i in 0..self.circuit_input_index {
            partitions.add_partition(RoutingTarget::CircuitInput(CircuitInput { index: i }));
        }

        for gate in 0..self.gate_ids.len() {
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

struct RoutingTargetPartitions {
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
