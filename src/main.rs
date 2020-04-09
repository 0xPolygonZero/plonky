use std::time::Instant;

use plonky::{Curve, fft_precompute, fft_with_precomputation, FftPrecomputation, Field, msm_execute_parallel, msm_precompute, MsmPrecomputation, ProjectivePoint, recursive_verification_circuit, Tweedledee, TWEEDLEDEE_GENERATOR_PROJECTIVE, Tweedledum};

const DEGREE: usize = 1 << 12;
const MSM_COUNT: usize = 10 + 1 + 7; // 10 wires, Z, and 7 components of t

type C = Tweedledee;
type SF = <C as Curve>::ScalarField;

fn main() {
    // Configure the main thread pool size.
    // rayon::ThreadPoolBuilder::new().num_threads(8).build_global().unwrap();

    let degree_pow = 13;
    recursive_verification_circuit::<Tweedledum>(degree_pow);

    // run_all_ffts();
    //
    // let mut generators = Vec::with_capacity(DEGREE);
    // let mut scalars = Vec::with_capacity(DEGREE);
    // for _i in 0..DEGREE {
    //     generators.push(TWEEDLEDEE_GENERATOR_PROJECTIVE);
    //     scalars.push(SF::rand());
    // }

    // Here's a quick Python snippet to calculate optimal window sizes:
    //     degree = 2**12
    //     parallelism = 4
    //     field_bits = 254
    //     group_ops = lambda w: 2**w + degree * ceil(field_bits / w) / parallelism
    //     min(range(1, 50), key=group_ops)
    // This is oversimplified though, as it doesn't account for the different summation methods we
    // use for different problem sizes. So some trial and error is needed to find the best config.

    // for w in 11..=11 {
    //     println!();
    //     println!("MSM WITH WINDOW SIZE {}", w);
    //     run_msms(w, &generators, &scalars);
    // }
}

fn run_all_ffts() {
    // As per the paper, we do 8 FFTs of size 4n, 5 FFTs of size 2n and 12 FFTs of size n.
    let start = Instant::now();

    let fft_8n_precomputation = fft_precompute(8 * DEGREE);
    for _i in 0..11 {
        run_fft(8 * DEGREE, &fft_8n_precomputation);
    }

    let fft_n_precomputation = fft_precompute(DEGREE);
    for _i in 0..7 {
        run_fft(DEGREE, &fft_n_precomputation);
    }

    println!("All FFTs took {}s", start.elapsed().as_secs_f64());
}

fn run_fft(size: usize, precomputation: &FftPrecomputation<SF>) {
    let mut coefficients = Vec::new();
    for i in 0..size {
        coefficients.push(SF::from_canonical_usize(i));
    }

    println!("Running FFT of size {}...", size);
    let start = Instant::now();
    let _result = fft_with_precomputation(&coefficients, precomputation);
    println!("Finished in {}s", start.elapsed().as_secs_f64());
}

fn run_msms(w: usize, generators: &[ProjectivePoint<C>], scalars: &[SF]) {
    println!("Precomputing...");
    let start = Instant::now();
    let precomputation = msm_precompute(generators, w);
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

    let start = Instant::now();
    for _i in 0..MSM_COUNT {
        run_msm(&precomputation, w, scalars);
    }
    println!("All MSMs took {}s", start.elapsed().as_secs_f64());
}

fn run_msm(precomputation: &MsmPrecomputation<C>, w: usize, scalars: &[SF]) {
    println!("Computing MSM in parallel...");
    let start = Instant::now();
    msm_execute_parallel(&precomputation, scalars, w);
    println!("Finished in {}s", start.elapsed().as_secs_f64());
}
