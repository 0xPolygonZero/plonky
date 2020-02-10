use std::time::Instant;

use executors::Executor;
use executors::threadpool_executor::ThreadPoolExecutor;

use plonky::{Bls12Scalar, G1_GENERATOR};

const DEGREE: usize = 1 << 17;

fn main() {
    let executor = ThreadPoolExecutor::new(num_cpus::get() * 2);

    for _i in 0..DEGREE {
        executor.execute(|| {
            let scalar = Bls12Scalar::rand();
            scalar * G1_GENERATOR;
        });
    }

    let now = Instant::now();
    println!("Executing {} multiplications...", DEGREE);
    executor.shutdown().expect("pool to shut down");
    let duration = now.elapsed();
    println!("Completed in {:.2}s ({:.2}ms per mul)",
             duration.as_secs_f64(),
             duration.div_f64(DEGREE as f64).as_secs_f64() * 1000.0);
}
