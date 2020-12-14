use criterion::{black_box, Criterion};
use criterion::criterion_group;
use criterion::criterion_main;

use plonky::{Bls12377Base, Field};

fn criterion_benchmark(c: &mut Criterion) {
    let x = Bls12377Base::from_canonical([11111111, 22222222, 33333333, 44444444, 55555555, 66666666]);
    let y = Bls12377Base::from_canonical([44444444, 55555555, 66666666, 77777777, 88888888, 99999999]);

    c.bench_function("Bls12Base field addition", move |b| b.iter(|| {
        black_box(y) + black_box(x)
    }));

    c.bench_function("Bls12Base field subtraction (no underflow)", move |b| b.iter(|| {
        black_box(y) - black_box(x)
    }));

    c.bench_function("Bls12Base field subtraction (underflow)", move |b| b.iter(|| {
        black_box(x) - black_box(y)
    }));

    c.bench_function("Bls12Base field multiplication", move |b| b.iter(|| {
        black_box(x) * black_box(y)
    }));

    c.bench_function("Bls12Base field squaring", move |b| b.iter(|| {
        black_box(x).square()
    }));

    c.bench_function("Bls12Base field inversion", move |b| b.iter(|| {
        black_box(x).multiplicative_inverse()
    }));

    c.bench_function("Bls12Base field exp", move |b| b.iter(|| {
        black_box(x).exp(black_box(y))
    }));

}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
