use criterion::{black_box, Criterion};
use criterion::criterion_group;
use criterion::criterion_main;

use plonky::{Bls12377Scalar, Field, rescue_hash_1_to_1};

fn criterion_benchmark(c: &mut Criterion) {
    let x = Bls12377Scalar::from_canonical([11111111, 22222222, 33333333, 44444444]);
    let y = Bls12377Scalar::from_canonical([44444444, 55555555, 66666666, 77777777]);

    c.bench_function("Bls12Scalar field addition", move |b| b.iter(|| {
        black_box(y) + black_box(x)
    }));

    c.bench_function("Bls12Scalar field subtraction (no underflow)", move |b| b.iter(|| {
        black_box(y) - black_box(x)
    }));

    c.bench_function("Bls12Scalar field subtraction (underflow)", move |b| b.iter(|| {
        black_box(x) - black_box(y)
    }));

    c.bench_function("Bls12Scalar field multiplication", move |b| b.iter(|| {
        black_box(x) * black_box(y)
    }));

    c.bench_function("Bls12Scalar field squaring", move |b| b.iter(|| {
        black_box(x).square()
    }));

    c.bench_function("Bls12Scalar field inversion", move |b| b.iter(|| {
        black_box(x).multiplicative_inverse()
    }));

    c.bench_function("Bls12Scalar Rescue hash", move |b| b.iter(|| {
        rescue_hash_1_to_1(black_box(x), 128)
    }));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
