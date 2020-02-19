use std::time::Instant;

use plonky::{Bls12Scalar, G1_GENERATOR, msm_execute, msm_precompute};

const DEGREE: usize = 1 << 13;

fn main() {
    let w = 15;
    let mut generators = Vec::with_capacity(DEGREE);
    let mut scalars = Vec::with_capacity(DEGREE);
    for _i in 0..DEGREE {
        generators.push(G1_GENERATOR);
        scalars.push(Bls12Scalar::rand());
    }

    let start = Instant::now();
    println!("Precomputing...");
    let precomputation = msm_precompute(&generators, w);
    println!("Finished in {}s", start.elapsed().as_secs_f64());

    let start = Instant::now();
    println!("Computing MSM...");
    let result = msm_execute(&precomputation, &scalars, w);
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!("Result: {:?}", result);
}
