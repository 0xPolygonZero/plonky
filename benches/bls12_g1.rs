use criterion::{black_box, Criterion};
use criterion::criterion_group;
use criterion::criterion_main;
use num::BigUint;

use plonky::{Bls12Base, Bls12Scalar, G1ProjectivePoint, mul_6_6, G1AffinePoint, G1_GENERATOR_AFFINE};

fn criterion_benchmark(c: &mut Criterion) {
    // We want a scalar with a Hamming weight of 0.5, to simulate the "average case".
    let s_part = 0b1010101010101010101010101010101010101010101010101010101010101010u64;
    let s = Bls12Scalar { limbs: [s_part, s_part, s_part, s_part] };

    let p1_affine = G1_GENERATOR_AFFINE;
    let p2_affine = (p1_affine + p1_affine).to_affine();

    // Use a non-zero z to make sure we don't hit any fast paths checking for z=1.
    let p1_projective = G1ProjectivePoint { x: p1_affine.x * 2, y: p1_affine.y * 2, z: Bls12Base::TWO };
    let p2_projective = p1_projective + p1_projective;

    c.bench_function("BLS12 G1 affine + affine = projective addition", move |b| b.iter(|| {
        let _result: G1ProjectivePoint = black_box(p1_affine) + black_box(p2_affine);
    }));

    c.bench_function("BLS12 G1 projective + affine = projective addition", move |b| b.iter(|| {
        let _result: G1ProjectivePoint = black_box(p1_projective) + black_box(p2_affine);
    }));

    c.bench_function("BLS12 G1 projective + projective = projective addition", move |b| b.iter(|| {
        let _result: G1ProjectivePoint = black_box(p1_projective) + black_box(p2_projective);
    }));

    c.bench_function("BLS12 G1 projective doubling", move |b| b.iter(|| {
        black_box(p1_projective).double();
    }));

    c.bench_function("BLS12 G1 projective multiplication", move |b| b.iter(|| {
        black_box(s) * black_box(p1_projective);
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
