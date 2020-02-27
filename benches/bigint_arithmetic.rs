use criterion::{black_box, Criterion};
use criterion::criterion_group;
use criterion::criterion_main;

use plonky::{cmp_6_6, mul_6_6};

fn criterion_benchmark(c: &mut Criterion) {
    let x = [11111111, 22222222, 33333333, 44444444, 55555555, 66666666];
    let y = [44444444, 55555555, 66666666, 77777777, 88888888, 99999999];

    c.bench_function("[u64; 6] widening multiplication", move |b| b.iter(|| {
        mul_6_6(black_box(x), black_box(y))
    }));

    c.bench_function("[u64; 6] comparison (lhs == rhs)", move |b| b.iter(|| {
        cmp_6_6(black_box(x), black_box(x))
    }));

    c.bench_function("[u64; 6] comparison (lhs < rhs)", move |b| b.iter(|| {
        cmp_6_6(black_box(x), black_box(y))
    }));

    c.bench_function("[u64; 6] comparison (lhs > rhs)", move |b| b.iter(|| {
        cmp_6_6(black_box(y), black_box(x))
    }));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
