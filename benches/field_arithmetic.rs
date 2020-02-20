use criterion::{black_box, Criterion};
use criterion::criterion_group;
use criterion::criterion_main;
use num::BigUint;

use plonky::{Bls12Base, mul_6_6};

fn criterion_benchmark(c: &mut Criterion) {
    let x = [11111111, 22222222, 33333333, 44444444, 55555555, 66666666];
    let y = [44444444, 55555555, 66666666, 77777777, 88888888, 99999999];

    let x_bls12base = Bls12Base::from_canonical(x);
    let y_bls12base = Bls12Base::from_canonical(y);

    let x_biguint = u64_slice_to_biguint(&x);
    let y_biguint = u64_slice_to_biguint(&y);
    let m_biguint = u64_slice_to_biguint(&Bls12Base::ORDER);

    {
        let x_biguint = x_biguint.clone();
        let y_biguint = y_biguint.clone();
        c.bench_function("BigUint widening multiplication", move |b| b.iter(|| {
            let x_biguint = x_biguint.clone();
            let y_biguint = y_biguint.clone();
            black_box(x_biguint) * black_box(y_biguint);
        }));
    }

    c.bench_function("[u64] widening multiplication", move |b| b.iter(|| {
        mul_6_6(black_box(x), black_box(y));
    }));

    {
        let x_biguint = x_biguint.clone();
        let y_biguint = y_biguint.clone();
        let m_biguint = m_biguint.clone();
        c.bench_function("BigUint field multiplication", move |b| b.iter(|| {
            let x_biguint = x_biguint.clone();
            let y_biguint = y_biguint.clone();
            let m_biguint = m_biguint.clone();
            black_box(x_biguint) * black_box(y_biguint) % black_box(m_biguint);
        }));
    }

    c.bench_function("Bls12Base field addition", move |b| b.iter(|| {
        black_box(y_bls12base) + black_box(x_bls12base);
    }));

    c.bench_function("Bls12Base field subtraction", move |b| b.iter(|| {
        black_box(y_bls12base) - black_box(x_bls12base);
    }));

    c.bench_function("Bls12Base field multiplication", move |b| b.iter(|| {
        black_box(x_bls12base) * black_box(y_bls12base);
    }));

    c.bench_function("Bls12Base field squaring", move |b| b.iter(|| {
        black_box(x_bls12base).square()
    }));
}

fn u64_slice_to_biguint(n: &[u64]) -> BigUint {
    let mut bytes_le = Vec::new();
    for n_i in n {
        for j in 0..8 {
            bytes_le.push((n_i >> j * 8) as u8);
        }
    }
    BigUint::from_bytes_le(&bytes_le)
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
