//! Runs MSMs with various thread pool sizes, and outputs timing data.

use std::time::Instant;

use plonky::{Curve, Field, msm_execute_parallel, msm_precompute, Tweedledum, ProjectivePoint};

fn main() {
    type C = Tweedledum;
    type SF = <C as Curve>::ScalarField;
    let size_log = 14;
    let size = 1 << size_log;
    let window_size = 12;

    let mut generators: Vec<_> = (0..size)
        .map(|_| C::convert(SF::rand()) * C::GENERATOR_PROJECTIVE)
        .collect();

    let scalars: Vec<_> = (0..size)
        .map(|_| SF::rand())
        .collect();

    let precomputation = msm_precompute(&generators, window_size);

    for threads in 1..=100 {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
        let start = Instant::now();
        pool.install(|| msm_execute_parallel(&precomputation, &scalars));
        println!("MSM with size=2^{}, threads={} took {}s", size_log, threads, start.elapsed().as_secs_f64());
    }
}