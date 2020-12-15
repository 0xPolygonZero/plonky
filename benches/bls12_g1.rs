use criterion::{black_box, Criterion};
use criterion::criterion_group;
use criterion::criterion_main;

use plonky::{Bls12377, Bls12377Base, Bls12377Scalar, Curve, Field, ProjectivePoint};

fn criterion_benchmark(c: &mut Criterion) {
    // We want a scalar with a Hamming weight of 0.5, to simulate the "average case".
    let s_part = 0b1010101010101010101010101010101010101010101010101010101010101010u64;
    let s = Bls12377Scalar { limbs: [s_part, s_part, s_part, s_part] };

    let p1_affine = Bls12377::GENERATOR_AFFINE;
    let p2_affine = (p1_affine + p1_affine).to_affine();

    // Use a non-zero z to make sure we don't hit any fast paths checking for z=1.
    let p1_projective = ProjectivePoint::nonzero(p1_affine.x * 2, p1_affine.y * 2, Bls12377Base::TWO);
    let p2_projective = p1_projective + p1_projective;

    c.bench_function("BLS12 G1 affine + affine = projective addition", move |b| b.iter(|| {
        black_box(p1_affine) + black_box(p2_affine)
    }));

    c.bench_function("BLS12 G1 projective + affine = projective addition", move |b| b.iter(|| {
        black_box(p1_projective) + black_box(p2_affine)
    }));

    c.bench_function("BLS12 G1 projective + projective = projective addition", move |b| b.iter(|| {
        black_box(p1_projective) + black_box(p2_projective)
    }));

    c.bench_function("BLS12 G1 projective doubling", move |b| b.iter(|| {
        black_box(p1_projective).double()
    }));

    c.bench_function("BLS12 G1 affine doubling", move |b| b.iter(|| {
        black_box(p1_affine).double()
    }));

    c.bench_function("BLS12 G1 projective multiplication", move |b| b.iter(|| {
        black_box(Bls12377::convert(s)) * black_box(p1_projective)
    }));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
