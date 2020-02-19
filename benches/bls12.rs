use criterion::{black_box, Criterion};
use criterion::criterion_group;
use criterion::criterion_main;
use num::BigUint;

use plonky::{Bls12Base, Bls12Scalar, G1ProjectivePoint, mul_6_6};

fn criterion_benchmark(c: &mut Criterion) {
    // We want a scalar with a Hamming weight of 0.5, to simulate the "average case".
    let s_part = 0b1010101010101010101010101010101010101010101010101010101010101010u64;
    let s = Bls12Scalar { limbs: [s_part, s_part, s_part, s_part] };

    let g_x = Bls12Base { limbs: [11111111, 22222222, 33333333, 44444444, 55555555, 66666666] };
    let g_y = Bls12Base { limbs: [44444444, 55555555, 66666666, 77777777, 88888888, 99999999] };
    let g = G1ProjectivePoint { x: g_x, y: g_y, z: Bls12Base::ONE };

    c.bench_function("BLS12 G1 addition", move |b| b.iter(|| {
        black_box(g) + black_box(g);
    }));

    c.bench_function("BLS12 G1 doubling", move |b| b.iter(|| {
        black_box(g).double();
    }));

    c.bench_function("BLS12 G1 multiplication", move |b| b.iter(|| {
        black_box(s) * black_box(g);
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
