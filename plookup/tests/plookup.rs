use anyhow::Result;
use num::{BigUint, Integer};
use plonky::{biguint_to_field, field_to_biguint, Curve, Field, Tweedledee};
use plookup::plookup::prove;
use plookup::table::Table;
use plookup::verifier::verify;
use rand::seq::SliceRandom;
use rand::{thread_rng, Rng};

fn add<F: Field>(a: [F; 2]) -> F {
    a[0] + a[1]
}

fn add_4_bits<F: Field>(a: [F; 2]) -> F {
    biguint_to_field(
        (field_to_biguint(a[0]) + field_to_biguint(a[1])).mod_floor(&BigUint::from(1usize << 4)),
    )
}

fn xor_4_bits<F: Field>(a: [F; 2]) -> F {
    biguint_to_field(
        (field_to_biguint(a[0]) ^ field_to_biguint(a[1])).mod_floor(&BigUint::from(1usize << 4)),
    )
}

fn mul_3_times_3_bits<F: Field>(a: [F; 3]) -> F {
    biguint_to_field(
        (field_to_biguint(a[0]) * field_to_biguint(a[1]) * field_to_biguint(a[2]))
            .mod_floor(&BigUint::from(1usize << 3)),
    )
}

#[test]
fn test_plookup() -> Result<()> {
    type C = Tweedledee;
    type SF = <C as Curve>::ScalarField;
    let mut rng = thread_rng();
    let n = 15;
    let t = (0..rng.gen_range(1, n / 2))
        .map(|_| SF::rand())
        .collect::<Vec<_>>();
    let f = (0..n)
        .map(|_| *t.choose(&mut rng).unwrap())
        .collect::<Vec<_>>();
    let proof = prove::<C>(&f, &t)?;
    verify(&t, &proof)?;
    Ok(())
}

#[test]
fn test_plookup_table() -> Result<()> {
    type C = Tweedledee;
    type SF = <C as Curve>::ScalarField;
    let mut rng = thread_rng();
    let n = 15;
    const WIDTH: usize = 10;
    let t = Table(
        (0..rng.gen_range(1, n / 2))
            .map(|_| {
                let mut arr = [SF::ZERO; WIDTH];
                (0..WIDTH).for_each(|i| arr[i] = SF::rand());
                arr
            })
            .collect::<Vec<_>>(),
    );
    let f = Table(
        (0..n)
            .map(|_| *t.0.choose(&mut rng).unwrap())
            .collect::<Vec<_>>(),
    );
    let proof = t.prove_row::<C>(&f)?;
    t.verify(&proof)?;
    Ok(())
}
#[test]
fn test_plookup_function() -> Result<()> {
    const WIDTH: usize = 3;
    type C = Tweedledee;
    type SF = <C as Curve>::ScalarField;
    let mut rng = thread_rng();
    let n: usize = 15;
    let domain = (0..rng.gen_range(1, n / 2))
        .map(|_| [SF::rand(), SF::rand()])
        .collect::<Vec<_>>();
    let t = Table::<SF, WIDTH>::from_function(&add, &domain);
    let f = Table(
        (0..n)
            .map(|_| *t.0.choose(&mut rng).unwrap())
            .collect::<Vec<_>>(),
    );
    let proof = t.prove_row::<C>(&f)?;
    t.verify(&proof)?;
    Ok(())
}

#[test]
fn test_add_4_bits() -> Result<()> {
    const WIDTH: usize = 3;
    const BITS: usize = 4;
    type C = Tweedledee;
    type SF = <C as Curve>::ScalarField;
    let mut rng = thread_rng();
    let n: usize = 15;
    let domain = (0..(1 << BITS))
        .map(|i| SF::from_canonical_usize(i))
        .collect::<Vec<_>>();
    let t = Table::<SF, WIDTH>::from_function_cartesian(&add_4_bits::<SF>, &domain);
    let f = Table(
        (0..n)
            .map(|_| *t.0.choose(&mut rng).unwrap())
            .collect::<Vec<_>>(),
    );
    let proof = t.prove_row::<C>(&f)?;
    t.verify(&proof)?;
    Ok(())
}

#[test]
fn test_xor_4_bits() -> Result<()> {
    const WIDTH: usize = 3;
    const BITS: usize = 4;
    type C = Tweedledee;
    type SF = <C as Curve>::ScalarField;
    let mut rng = thread_rng();
    let n: usize = 15;
    let domain = (0..(1 << BITS))
        .map(|i| SF::from_canonical_usize(i))
        .collect::<Vec<_>>();
    let t = Table::<SF, WIDTH>::from_function_cartesian(&xor_4_bits::<SF>, &domain);
    let f = Table(
        (0..n)
            .map(|_| *t.0.choose(&mut rng).unwrap())
            .collect::<Vec<_>>(),
    );
    let proof = t.prove_row::<C>(&f)?;
    t.verify(&proof)?;
    Ok(())
}

#[test]
fn test_mul_3_bits() -> Result<()> {
    const WIDTH: usize = 4;
    const BITS: usize = 3;
    type C = Tweedledee;
    type SF = <C as Curve>::ScalarField;
    let mut rng = thread_rng();
    let n: usize = 15;
    let domain = (0..(1 << BITS))
        .map(|i| SF::from_canonical_usize(i))
        .collect::<Vec<_>>();
    let t = Table::<SF, WIDTH>::from_function_cartesian(&mul_3_times_3_bits::<SF>, &domain);
    let f = Table(
        (0..n)
            .map(|_| *t.0.choose(&mut rng).unwrap())
            .collect::<Vec<_>>(),
    );
    let proof = t.prove_row::<C>(&f)?;
    t.verify(&proof)?;
    Ok(())
}
