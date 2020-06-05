use std::time::Instant;

use plonky::{BufferGate, Circuit, CircuitBuilder, PartialWitness, recursive_verification_circuit, Tweedledee, Tweedledum};

const INNER_PROOF_DEGREE_POW: usize = 13;
const INNER_PROOF_DEGREE: usize = 1 << INNER_PROOF_DEGREE_POW;
const SECURITY_BITS: usize = 128;

fn main() {
    println!("Generating inner circuit");
    let start = Instant::now();
    let inner_circuit = generate_trivial_circuit();
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

    println!("Generating inner witness");
    let start = Instant::now();
    let inner_witness = inner_circuit.generate_witness(PartialWitness::new());
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

    println!("Generating inner proof");
    let start = Instant::now();
    let inner_proof = inner_circuit.generate_proof::<Tweedledee>(inner_witness).unwrap();
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

    println!("Generating recursion circuit...");
    let start = Instant::now();
    let recursion_circuit = recursive_verification_circuit::<Tweedledee, Tweedledum>(
        INNER_PROOF_DEGREE_POW, SECURITY_BITS);
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!("Gate count: {}", recursion_circuit.circuit.degree());
    println!();

    // Populate inputs.
    let mut recursion_inputs = PartialWitness::new();
    if let Err(e) = recursion_circuit.proof.populate_witness(&mut recursion_inputs, inner_proof) {
        panic!("Failed to populate inputs: {:?}", e);
    }

    println!("Generating recursion witness...");
    let start = Instant::now();
    let recursion_witness = recursion_circuit.circuit.generate_witness(recursion_inputs);
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

    println!("Generating recursion proof...");
    let start = Instant::now();
    recursion_circuit.circuit.generate_proof::<Tweedledum>(recursion_witness).unwrap();
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

    println!("Verifying proof...");
    let start = Instant::now();
    // TODO
    println!("Finished in {}s", start.elapsed().as_secs_f64());
}

fn generate_trivial_circuit() -> Circuit<Tweedledum> {
    let mut builder = CircuitBuilder::new(SECURITY_BITS);
    builder.route_public_inputs();
    while builder.num_gates() < INNER_PROOF_DEGREE {
        builder.add_gate_no_constants(BufferGate::new(builder.num_gates()));
    }
    builder.build()
}
