use std::time::Instant;

use plonky::{Bls12Scalar, msm_execute, msm_precompute, msm_execute_parallel, G1_GENERATOR_PROJECTIVE, G1ProjectivePoint, fft, fft_precompute, FftPrecomputation, fft_with_precomputation, MsmPrecomputation};

const DEGREE: usize = 1 << 17;

fn main() {
    // Configure the main thread pool size.
    rayon::ThreadPoolBuilder::new().num_threads(12).build_global().unwrap();

    run_all_ffts();

    let mut generators = Vec::with_capacity(DEGREE);
    let mut scalars = Vec::with_capacity(DEGREE);
    for _i in 0..DEGREE {
        generators.push(G1_GENERATOR_PROJECTIVE);
        scalars.push(Bls12Scalar::rand());
    }

    // Here's a quick Python snippet to calculate optimal window sizes:
    //     degree = 2**17
    //     parallelism = 8
    //     field_bits = 253
    //     group_ops = lambda w: 2**w + degree * ceil(field_bits / w) / parallelism
    //     min(range(1, 50), key=group_ops)
    // This is oversimplified though, as it doesn't account for the different summation methods we
    // use for different problem sizes. So some trial and error is needed to find the best config.

    for w in 14..=14 {
        println!();
        println!("MSM WITH WINDOW SIZE {}", w);
        run_msms(w, &generators, &scalars);
    }
}

fn run_all_ffts() {
    // As per the paper, we do 8 FFTs of size 4n, 5 FFTs of size 2n and 12 FFTs of size n.
    let start = Instant::now();

    let fft_4n_precomputation = fft_precompute(4 * DEGREE);
    for _i in 0..8 {
        run_fft(4 * DEGREE, &fft_4n_precomputation);
    }

    let fft_2n_precomputation = fft_precompute(2 * DEGREE);
    for _i in 0..5 {
        run_fft(2 * DEGREE, &fft_2n_precomputation);
    }

    let fft_n_precomputation = fft_precompute(DEGREE);
    for _i in 0..12 {
        run_fft(DEGREE, &fft_n_precomputation);
    }

    println!("All FFTs took {}s", start.elapsed().as_secs_f64());
}

fn run_fft(size: usize, precomputation: &FftPrecomputation) {
    let mut coefficients = Vec::new();
    for i in 0..size {
        coefficients.push(Bls12Scalar::from_canonical_usize(i));
    }

    println!("Running FFT of size {}...", size);
    let start = Instant::now();
    let _result = fft_with_precomputation(&coefficients, precomputation);
    println!("Finished in {}s", start.elapsed().as_secs_f64());
}

fn run_msms(w: usize, generators: &[G1ProjectivePoint], scalars: &[Bls12Scalar]) {
    println!("Precomputing...");
    let start = Instant::now();
    let precomputation = msm_precompute(generators, w);
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

//    println!("Computing MSM with one thread...");
//    let start = Instant::now();
//    let result = msm_execute(&precomputation, scalars, w);
//    println!("Finished in {}s", start.elapsed().as_secs_f64());
//    println!("Result: {:?}", result.to_affine());
//    println!();

    let start = Instant::now();
    for _i in 0..9 {
        run_msm(&precomputation, w, generators, scalars);
    }
    println!("All MSMs took {}s", start.elapsed().as_secs_f64());
}

fn run_msm(precomputation: &MsmPrecomputation, w: usize, generators: &[G1ProjectivePoint], scalars: &[Bls12Scalar]) {
    println!("Computing MSM in parallel...");
    let start = Instant::now();
    let result = msm_execute_parallel(&precomputation, scalars, w);
    println!("Finished in {}s", start.elapsed().as_secs_f64());
//    println!("Result: {:?}", result.to_affine());
}
