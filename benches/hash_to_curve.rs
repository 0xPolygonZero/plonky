use criterion::{black_box, Criterion};
use criterion::criterion_group;
use criterion::criterion_main;

use plonky::{Field, TweedledeeBase, hash_base_field_to_curve, blake_hash_base_field_to_curve, Tweedledee};

fn criterion_benchmark(c: &mut Criterion) {
    let x = TweedledeeBase::rand();

    c.bench_function("Hash using Rescue", move |b| b.iter(|| {
        hash_base_field_to_curve::<Tweedledee>(black_box(x), 128)
    }));

    c.bench_function("Hash using Blake", move |b| b.iter(|| {
        blake_hash_base_field_to_curve::<Tweedledee>(black_box(x))
    }));

}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
