use crate::Field;
use std::time::Instant;
use std::collections::HashMap;
use crate::plonk_gates::Gate;

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
    fn generate_witness(&self) {
        let start = Instant::now();
        println!("Witness generation took {}s", start.elapsed().as_secs_f32());
        todo!()
    }
}

#[derive(Eq, PartialEq, Hash, Debug)]
pub struct CircuitInput {
    pub index: usize,
}

#[derive(Eq, PartialEq, Hash, Debug)]
pub struct GateInput {
    pub gate: usize,
    pub input: usize,
}

#[derive(Eq, PartialEq, Hash, Debug)]
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
        todo!()
    }

    pub fn add_public_inputs(&mut self, n: usize) -> Vec<GateInput> {
        todo!()
    }

    pub fn add_circuit_input(&mut self) -> CircuitInput {
        let index = self.circuit_input_index;
        self.circuit_input_index += 1;
        CircuitInput { index }
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

    pub fn build(&self) -> Circuit<F> {
        // TODO: Shift indices
        // TODO: Add dummy gates through which public inputs can be routed to other gates.
        todo!()
    }

    fn get_routing_partitions(&self) -> RoutingTargetPartitions {
        todo!()
    }
}

struct RoutingTargetPartitions {
    partitions: Vec<Vec<RoutingTarget>>,
    indices: HashMap<RoutingTarget, usize>,
}

impl RoutingTargetPartitions {
    fn to_gate_inputs(&self) -> GateInputPartitions {
        todo!()
    }
}

struct GateInputPartitions {
    partitions: Vec<Vec<GateInput>>,
    indices: HashMap<GateInput, usize>,
}
