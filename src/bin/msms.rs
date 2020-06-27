//! Runs MSMs with various thread pool sizes, and outputs timing data.

use std::time::{Instant, Duration};

use plonky::{Curve, Field, msm_execute_parallel, msm_precompute, Tweedledum, ProjectivePoint, MsmPrecomputation};

type C = Tweedledum;
type SF = <C as Curve>::ScalarField;

fn main() {
    let min_window_size = 11;
    let max_window_size = 12;
    let min_threads = 1;
    let max_threads = 100;
    let repetitions = 3;
    let terms_log = 14;
    let terms = 1 << terms_log;

    let mut generators: Vec<_> = (0..terms)
        .map(|_| C::convert(SF::rand()) * C::GENERATOR_PROJECTIVE)
        .collect();

    for window_size in min_window_size..=max_window_size {
        let precomputation = msm_precompute(&generators, window_size);

        for threads in min_threads..=max_threads {
            let mut durations = Vec::new();
            for _ in 0..repetitions {
                durations.push(time_msm(terms, window_size, threads, &precomputation));
            }
            durations.sort();
            let average_duration = durations.iter().sum::<Duration>() / durations.len() as u32;
            let median_duration = durations[durations.len() / 2];
            println!("MSMs with terms=2^{}, threads={}, window_size={}: average={:.4}s, median={:.4}s",
                     terms_log, threads, window_size,
                     average_duration.as_secs_f64(), median_duration.as_secs_f64());
        }
    }
}

fn time_msm(terms: usize, window_size: usize, threads: usize, precomputation: &MsmPrecomputation<C>) -> Duration {
    let scalars: Vec<_> = (0..terms)
        .map(|_| SF::rand())
        .collect();

    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();

    pool.install(|| {
        let start = Instant::now();
        msm_execute_parallel(&precomputation, &scalars);
        start.elapsed()
    })
}