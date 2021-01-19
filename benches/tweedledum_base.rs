use criterion::{black_box, Criterion};
use criterion::criterion_group;
use criterion::criterion_main;

use plonky::{Field, TweedledumBase};

fn criterion_benchmark(c: &mut Criterion) {
    let x = TweedledumBase::from_canonical([11111111, 22222222, 33333333, 44444444]);
    let y = TweedledumBase::from_canonical([44444444, 55555555, 66666666, 77777777]);

    c.bench_function("TweedledumBase field addition", move |b| b.iter(|| {
        black_box(y) + black_box(x)
    }));

    c.bench_function("TweedledumBase field subtraction (no underflow)", move |b| b.iter(|| {
        black_box(y) - black_box(x)
    }));

    c.bench_function("TweedledumBase field subtraction (underflow)", move |b| b.iter(|| {
        black_box(x) - black_box(y)
    }));

    c.bench_function("TweedledumBase field multiplication", move |b| b.iter(|| {
        black_box(x) * black_box(y)
    }));

    c.bench_function("TweedledumBase field squaring", move |b| b.iter(|| {
        black_box(x).square()
    }));

    c.bench_function("TweedledumBase field inversion", move |b| b.iter(|| {
        black_box(x).multiplicative_inverse()
    }));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
