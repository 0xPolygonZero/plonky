use std::time::Instant;

use plonky::{Curve, Field, recursive_verification_circuit, Tweedledee, Tweedledum};

const DEGREE: usize = 1 << 12;
const MSM_COUNT: usize = 10 + 1 + 7; // 10 wires, Z, and 7 components of t

type C = Tweedledee;
type SF = <C as Curve>::ScalarField;

fn main() {
    let degree_pow = 13;

    println!("Generating circuit...");
    let start = Instant::now();
    let recursive_circuit = recursive_verification_circuit::<Tweedledum>(degree_pow);
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!("Gate count: {}", recursive_circuit.circuit.num_gates());

    println!("Generating witness...");
    let start = Instant::now();
    todo!();
}
