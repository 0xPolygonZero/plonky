use criterion::{black_box, Criterion};
use criterion::criterion_group;
use criterion::criterion_main;

use plonky::{Bls12Base, cmp_6_6, mul_6_6};

fn criterion_benchmark(c: &mut Criterion) {
    let x = [11111111, 22222222, 33333333, 44444444, 55555555, 66666666];
    let y = [44444444, 55555555, 66666666, 77777777, 88888888, 99999999];

    let x_bls12base = Bls12Base::from_canonical(x);
    let y_bls12base = Bls12Base::from_canonical(y);

    c.bench_function("[u64; 6] widening multiplication", move |b| b.iter(|| {
        mul_6_6(black_box(x), black_box(y));
    }));

    c.bench_function("Bls12Base field addition", move |b| b.iter(|| {
        black_box(y_bls12base) + black_box(x_bls12base);
    }));

    c.bench_function("[u64; 6] comparison (lhs == rhs)", move |b| b.iter(|| {
        cmp_6_6(black_box(x), black_box(x));
    }));

    c.bench_function("[u64; 6] comparison (lhs < rhs)", move |b| b.iter(|| {
        cmp_6_6(black_box(x), black_box(y));
    }));

    c.bench_function("[u64; 6] comparison (lhs > rhs)", move |b| b.iter(|| {
        cmp_6_6(black_box(y), black_box(x));
    }));

    c.bench_function("Bls12Base field subtraction (no underflow)", move |b| b.iter(|| {
        black_box(y_bls12base) - black_box(x_bls12base);
    }));

    c.bench_function("Bls12Base field subtraction (underflow)", move |b| b.iter(|| {
        black_box(x_bls12base) - black_box(y_bls12base);
    }));

    c.bench_function("Bls12Base field multiplication", move |b| b.iter(|| {
        black_box(x_bls12base) * black_box(y_bls12base);
    }));

    c.bench_function("Bls12Base field squaring", move |b| b.iter(|| {
        black_box(x_bls12base).square()
    }));

    c.bench_function("Bls12Base field inversion", move |b| b.iter(|| {
        black_box(x_bls12base).multiplicative_inverse()
    }));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
