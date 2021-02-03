use criterion::{black_box, Criterion};
use criterion::criterion_group;
use criterion::criterion_main;

use plonky::{Field, TweedledeeBase, DairaRepr};

fn criterion_benchmark(c: &mut Criterion) {
    let x = TweedledeeBase::from_canonical([11111111, 22222222, 33333333, 44444444]);
    let y = TweedledeeBase::from_canonical([44444444, 55555555, 66666666, 77777777]);

    let d = [44444444u64, 55555555, 66666666, 77777777, 88888888, 99999999, 11111111, 22222222];

    c.bench_function("TweedledeeBase field addition", move |b| b.iter(|| {
        black_box(y) + black_box(x)
    }));

    c.bench_function("TweedledeeBase field subtraction (no underflow)", move |b| b.iter(|| {
        black_box(y) - black_box(x)
    }));

    c.bench_function("TweedledeeBase field subtraction (underflow)", move |b| b.iter(|| {
        black_box(x) - black_box(y)
    }));

    c.bench_function("TweedledeeBase field multiplication", move |b| b.iter(|| {
        black_box(x) * black_box(y)
    }));

    c.bench_function("TweedledeeBase field squaring", move |b| b.iter(|| {
        black_box(x).square()
    }));

    c.bench_function("TweedledeeBase Daira reduction", move |b| b.iter(|| {
        TweedledeeBase::_reduce(black_box(d))
    }));

    c.bench_function("TweedledeeBase field inversion", move |b| b.iter(|| {
        black_box(x).multiplicative_inverse()
    }));

    c.bench_function("TweedledeeBase field exp", move |b| b.iter(|| {
        black_box(x).exp(black_box(y))
    }));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
