use std::time::Instant;

use plonky::{recursive_verification_circuit, Tweedledum, Tweedledee, PartialWitness, Curve};

fn main() {
    let degree_pow = 13;

    println!("Generating circuit...");
    let start = Instant::now();
    let recursive_circuit = recursive_verification_circuit::<Tweedledee, Tweedledum>(degree_pow, 128);
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!("Gate count: {}", recursive_circuit.circuit.degree());
    println!();

    // Dummy inputs.
    let mut inputs = PartialWitness::new();
    inputs.set_point_target(recursive_circuit.proof.c_plonk_z, Tweedledum::GENERATOR_AFFINE);
    inputs.set_point_target(recursive_circuit.proof.halo_g, Tweedledum::GENERATOR_AFFINE);
    // TODO: Populate other targets with dummy data.

    println!("Generating witness...");
    let start = Instant::now();
    let witness = recursive_circuit.circuit.generate_witness(inputs);
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

    println!("Generating proof...");
    let start = Instant::now();
    recursive_circuit.circuit.generate_proof::<Tweedledum>(witness);
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

    println!("Verifying proof...");
    let start = Instant::now();
    // TODO
    println!("Finished in {}s", start.elapsed().as_secs_f64());
}
