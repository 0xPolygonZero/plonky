use std::time::Instant;

use plonky::{recursive_verification_circuit, Tweedledum};

fn main() {
    let degree_pow = 13;

    println!("Generating circuit...");
    let start = Instant::now();
    let recursive_circuit = recursive_verification_circuit::<Tweedledum>(degree_pow);
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!("Gate count: {}", recursive_circuit.circuit.num_gates());
    println!();

    println!("Generating witness...");
    let start = Instant::now();
    // TODO
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

    println!("Generating proof...");
    let start = Instant::now();
    // TODO
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

    println!("Verifying proof...");
    let start = Instant::now();
    // TODO
    println!("Finished in {}s", start.elapsed().as_secs_f64());
}
