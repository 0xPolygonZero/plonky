use criterion::{black_box, Criterion};
use criterion::criterion_group;
use criterion::{BenchmarkId, criterion_main, SamplingMode};

use plonky::{Field, TweedledeeBase, fft_precompute, ifft_with_precomputation_power_of_2, fft_with_precomputation_power_of_2};
use std::time::Duration;

type F = TweedledeeBase;

const DEGREE_LOG_MIN: usize = 1;
const DEGREE_LOG_MAX: usize = 18;

fn degree_logs() -> Vec<usize> {
    (DEGREE_LOG_MIN..=DEGREE_LOG_MAX).collect()
}

fn fft(c: &mut Criterion) {
    let mut group = c.benchmark_group("fft");

    for degree_log in degree_logs() {
        let degree = 1 << degree_log;
        let precomputation = fft_precompute(degree);
        let coeffs: Vec<F> = (0..degree).map(|_| F::rand()).collect();
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("2_exp_{}", degree_log)),
            &degree_log, 
            |b, &_degree_log| {
                b.iter(|| {
                    fft_with_precomputation_power_of_2(black_box(&coeffs), &precomputation);
                });
            }
        );
    }
}

fn ifft(c: &mut Criterion) {
    let mut group = c.benchmark_group("ifft");

    for degree_log in degree_logs() {
        let degree = 1 << degree_log;
        let precomputation = fft_precompute(degree);
        let points: Vec<F> = (0..degree).map(|_| F::rand()).collect();
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("2_exp_{}", degree_log)),
            &degree_log, 
            |b, &_degree_log| {
                b.iter(|| {
                    ifft_with_precomputation_power_of_2(black_box(&points), &precomputation);
                });
            }
        );
    }
}

criterion_group!(
    name = benches;
    config = Criterion::default()
        .sample_size(10)
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(1));
    targets = fft, ifft
);

criterion_main!(benches);
