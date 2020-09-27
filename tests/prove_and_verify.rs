use anyhow::Result;
use plonky::{
    blake_hash_base_field_to_curve, msm_parallel, rescue_hash_1_to_1, verify_proof, AffinePoint,
    Base4SumGate, Circuit, CircuitBuilder, Curve, CurveMulOp, Field, HaloCurve, PartialWitness,
    Target, Tweedledee, Tweedledum, Wire, Witness,
};
use rand::{thread_rng, Rng};
use std::time::Instant;

fn get_trivial_circuit<C: HaloCurve>(x: C::ScalarField) -> (Circuit<C>, Witness<C::ScalarField>) {
    let mut builder = CircuitBuilder::<C>::new(128);
    let t = builder.constant_wire(x);
    builder.assert_zero(t);
    let mut partial_witness = PartialWitness::new();
    partial_witness.set_target(t, x);
    let circuit = builder.build();
    let witness = circuit.generate_witness(partial_witness);
    (circuit, witness)
}

#[test]
fn test_proof_trivial() -> Result<()> {
    let (circuit, witness) = get_trivial_circuit(<Tweedledee as Curve>::ScalarField::ZERO);
    let proof = circuit
        .generate_proof::<Tweedledum>(&witness, &[], true)
        .unwrap();
    verify_proof::<Tweedledee, Tweedledum>(&[], &proof, &[], &circuit.into(), true)?;

    Ok(())
}

#[test]
#[allow(clippy::same_item_push)]
fn test_proof_trivial_circuit_many_proofs() -> Result<()> {
    let mut old_proofs = Vec::new();
    for _ in 0..10 {
        let (circuit, witness) = get_trivial_circuit(<Tweedledee as Curve>::ScalarField::ZERO);
        let proof = circuit
            .generate_proof::<Tweedledum>(&witness, &[], true)
            .unwrap();
        let vk = circuit.to_vk();
        let old_proof = verify_proof::<Tweedledee, Tweedledum>(&[], &proof, &[], &vk, false)
            .expect("Invalid proof")
            .unwrap();
        old_proofs.push(old_proof);
    }
    let (circuit, witness) = get_trivial_circuit(<Tweedledee as Curve>::ScalarField::ZERO);
    let proof = circuit
        .generate_proof::<Tweedledum>(&witness, &old_proofs, true)
        .unwrap();
    let vk = circuit.to_vk();
    verify_proof::<Tweedledee, Tweedledum>(&[], &proof, &old_proofs, &vk, true)?;

    Ok(())
}

#[test]
fn test_proof_sum() -> Result<()> {
    let mut builder = CircuitBuilder::<Tweedledee>::new(128);
    let t1 = builder.constant_wire(<Tweedledee as Curve>::ScalarField::ZERO);
    let t2 = builder.constant_wire(<Tweedledee as Curve>::ScalarField::ZERO);
    let s = builder.add(t1, t2);
    builder.assert_zero(s);
    let mut partial_witness = PartialWitness::new();
    partial_witness.set_target(t1, <Tweedledee as Curve>::ScalarField::ZERO);
    partial_witness.set_target(t2, <Tweedledee as Curve>::ScalarField::ZERO);
    let circuit = builder.build();
    let witness = circuit.generate_witness(partial_witness);
    let proof = circuit
        .generate_proof::<Tweedledum>(&witness, &[], true)
        .unwrap();
    let vk = circuit.to_vk();
    verify_proof::<Tweedledee, Tweedledum>(&[], &proof, &[], &vk, true)?;

    Ok(())
}

#[test]
#[ignore]
fn test_proof_sum_big() -> Result<()> {
    let now = Instant::now();
    let mut builder = CircuitBuilder::<Tweedledee>::new(128);
    let ts = (0..10_000)
        .map(|i| builder.constant_wire(<Tweedledee as Curve>::ScalarField::from_canonical_usize(i)))
        .collect::<Vec<_>>();
    let s = builder.add_many(&ts);
    let x = builder.add_virtual_target();
    let z = builder.sub(s, x);
    builder.assert_zero(z);
    let mut partial_witness = PartialWitness::new();
    partial_witness.set_target(
        x,
        <Tweedledee as Curve>::ScalarField::from_canonical_usize((10_000 * 9_999) / 2),
    );
    dbg!(now.elapsed());
    let circuit = builder.build();
    dbg!(now.elapsed());
    let witness = circuit.generate_witness(partial_witness);
    dbg!(now.elapsed());
    let proof = circuit
        .generate_proof::<Tweedledum>(&witness, &[], true)
        .unwrap();
    dbg!(now.elapsed());
    let vk = circuit.to_vk();
    dbg!(now.elapsed());
    verify_proof::<Tweedledee, Tweedledum>(&[], &proof, &[], &vk, true)?;
    dbg!(now.elapsed());

    Ok(())
}

#[test]
fn test_proof_quadratic() -> Result<()> {
    let mut builder = CircuitBuilder::<Tweedledee>::new(128);
    let one = builder.one_wire();
    let t = builder.add_virtual_target();
    let t_sq = builder.square(t);
    let quad = builder.add_many(&[one, t, t_sq]);
    let seven = builder.constant_wire(<Tweedledee as Curve>::ScalarField::from_canonical_usize(7));
    let res = builder.sub(quad, seven);
    builder.assert_zero(res);
    let mut partial_witness = PartialWitness::new();
    partial_witness.set_target(t, <Tweedledee as Curve>::ScalarField::TWO);
    let circuit = builder.build();
    let witness = circuit.generate_witness(partial_witness);
    let proof = circuit
        .generate_proof::<Tweedledum>(&witness, &[], true)
        .unwrap();
    let vk = circuit.to_vk();
    verify_proof::<Tweedledee, Tweedledum>(&[], &proof, &[], &vk, true)?;

    Ok(())
}

#[test]
fn test_proof_quadratic_public_input() -> Result<()> {
    type F = <Tweedledee as Curve>::ScalarField;
    let mut builder = CircuitBuilder::<Tweedledee>::new(128);
    let seven_pi = builder.stage_public_input();
    builder.route_public_inputs();
    let one = builder.one_wire();
    let t = builder.add_virtual_target();
    let t_sq = builder.square(t);
    let quad = builder.add_many(&[one, t, t_sq]);
    let seven = seven_pi.routable_target();
    let res = builder.sub(quad, seven);
    builder.assert_zero(res);
    let mut partial_witness = PartialWitness::new();
    partial_witness.set_target(seven, F::from_canonical_usize(7));
    partial_witness.set_target(t, F::TWO);
    let circuit = builder.build();
    let witness = circuit.generate_witness(partial_witness);
    let proof = circuit
        .generate_proof::<Tweedledum>(&witness, &[], true)
        .unwrap();
    let vk = circuit.to_vk();
    verify_proof::<Tweedledee, Tweedledum>(&[F::from_canonical_usize(7)], &proof, &[], &vk, false)?;

    Ok(())
}

#[test]
fn test_proof_sum_public_input() -> Result<()> {
    type F = <Tweedledee as Curve>::ScalarField;
    let mut builder = CircuitBuilder::<Tweedledee>::new(128);
    let n = 100;
    let factorial_usize = (1..=n).sum();
    let factors_pis = builder.stage_public_inputs(n);
    builder.route_public_inputs();
    let res = builder.constant_wire(F::from_canonical_usize(factorial_usize));
    let mut factorial = builder.zero_wire();
    factors_pis.iter().for_each(|pi| {
        factorial = builder.add(factorial, pi.routable_target());
    });
    builder.copy(factorial, res);
    let mut partial_witness = PartialWitness::new();
    (0..n).for_each(|i| {
        partial_witness.set_public_input(factors_pis[i], F::from_canonical_usize(i + 1));
    });
    let circuit = builder.build();
    let witness = circuit.generate_witness(partial_witness);
    let now = Instant::now();
    let proof = circuit
        .generate_proof::<Tweedledum>(&witness, &[], true)
        .unwrap();
    dbg!(now.elapsed());
    let vk = circuit.to_vk();
    verify_proof::<Tweedledee, Tweedledum>(
        &(1..=n).map(F::from_canonical_usize).collect::<Vec<_>>(),
        &proof,
        &[],
        &vk,
        false,
    )?;
    dbg!(now.elapsed());

    Ok(())
}

#[test]
fn test_proof_factorial_public_input() -> Result<()> {
    type F = <Tweedledee as Curve>::ScalarField;
    let mut builder = CircuitBuilder::<Tweedledee>::new(128);
    let n = 20;
    let factorial_usize = (1..=n).product();
    let factors_pis = builder.stage_public_inputs(n);
    builder.route_public_inputs();
    let res = builder.constant_wire(F::from_canonical_usize(factorial_usize));
    let mut factorial = builder.one_wire();
    factors_pis.iter().for_each(|pi| {
        factorial = builder.mul(factorial, pi.routable_target());
    });
    builder.copy(factorial, res);
    let mut partial_witness = PartialWitness::new();
    (0..n).for_each(|i| {
        partial_witness.set_public_input(factors_pis[i], F::from_canonical_usize(i + 1));
    });
    let circuit = builder.build();
    let witness = circuit.generate_witness(partial_witness);
    let proof = circuit
        .generate_proof::<Tweedledum>(&witness, &[], true)
        .unwrap();
    let vk = circuit.to_vk();
    verify_proof::<Tweedledee, Tweedledum>(
        &(1..=n).map(F::from_canonical_usize).collect::<Vec<_>>(),
        &proof,
        &[],
        &vk,
        true,
    )?;

    Ok(())
}

#[test]
fn test_proof_public_input() -> Result<()> {
    // Set many random public inputs
    let n = thread_rng().gen_range(1, 200);
    let values = (0..n)
        .map(|_| <Tweedledee as Curve>::ScalarField::rand())
        .collect::<Vec<_>>();
    let mut builder = CircuitBuilder::<Tweedledee>::new(128);
    let pis = (0..n)
        .map(|_| builder.stage_public_input())
        .collect::<Vec<_>>();
    builder.route_public_inputs();
    let mut partial_witness = PartialWitness::new();
    (0..n).for_each(|i| partial_witness.set_public_input(pis[i], values[i]));
    let circuit = builder.build();
    let witness = circuit.generate_witness(partial_witness);
    let proof = circuit
        .generate_proof::<Tweedledum>(&witness, &[], true)
        .unwrap();
    let pis = circuit.get_public_inputs(&witness);
    // Check that the public inputs are set correctly in the proof.
    assert_eq!(values, pis);
    let vk = circuit.to_vk();
    verify_proof::<Tweedledee, Tweedledum>(&values, &proof, &[], &vk, true)?;

    Ok(())
}

#[test]
fn test_rescue_hash() -> Result<()> {
    type F = <Tweedledee as Curve>::ScalarField;
    let x = F::rand();
    let h = rescue_hash_1_to_1(x, 128);
    let mut builder = CircuitBuilder::<Tweedledee>::new(128);
    let t = builder.add_virtual_target();
    let h_pur = builder.rescue_hash_n_to_1(&[t]);
    let c = builder.constant_wire(h);
    let should_be_zero = builder.sub(h_pur, c);
    builder.assert_zero(should_be_zero);
    let mut partial_witness = PartialWitness::new();
    partial_witness.set_target(t, x);
    let circuit = builder.build();
    let witness = circuit.generate_witness(partial_witness);
    let proof = circuit
        .generate_proof::<Tweedledum>(&witness, &[], true)
        .unwrap();
    let vk = circuit.to_vk();
    verify_proof::<Tweedledee, Tweedledum>(&[], &proof, &[], &vk, true)?;

    Ok(())
}

#[test]
fn test_curve_add() -> Result<()> {
    type F = <Tweedledee as Curve>::ScalarField;

    let a = blake_hash_base_field_to_curve::<Tweedledum>(F::rand());
    let b = blake_hash_base_field_to_curve::<Tweedledum>(F::rand());
    let sum = (a + b).to_affine();

    let mut builder = CircuitBuilder::<Tweedledee>::new(128);

    let ta = builder.add_virtual_point_target();
    let tb = builder.add_virtual_point_target();
    let tsum_purported = builder.curve_add::<Tweedledum>(ta, tb);
    let tsum_true = builder.constant_affine_point(sum);
    builder.copy_curve(tsum_purported, tsum_true);

    let mut partial_witness = PartialWitness::new();
    partial_witness.set_point_target(ta, a);
    partial_witness.set_point_target(tb, b);

    let circuit = builder.build();
    let witness = circuit.generate_witness(partial_witness);

    let proof = circuit
        .generate_proof::<Tweedledum>(&witness, &[], true)
        .unwrap();

    let vk = circuit.to_vk();
    verify_proof::<Tweedledee, Tweedledum>(&[], &proof, &[], &vk, true)?;

    Ok(())
}

#[test]
fn test_curve_msm() -> Result<()> {
    type SF = <Tweedledee as Curve>::ScalarField;
    type BF = <Tweedledee as Curve>::BaseField;
    let n = 10;
    let xs = (0..n).map(|_| SF::rand()).collect::<Vec<_>>();
    let ps = (0..n)
        .map(|_| blake_hash_base_field_to_curve::<Tweedledee>(BF::rand()))
        .collect::<Vec<_>>();
    let res = msm_parallel(&xs, &AffinePoint::batch_to_projective(&ps), 8);
    let mut builder = CircuitBuilder::<Tweedledum>::new(128);
    let txs = builder.add_virtual_targets(n);
    let tps = builder.add_virtual_point_targets(n);
    let tres_purported = builder.curve_msm::<Tweedledee>(
        &(0..n)
            .map(|i| CurveMulOp {
                scalar: txs[i],
                point: tps[i],
            })
            .collect::<Vec<_>>(),
    );
    let tres_true = builder.constant_affine_point(res.to_affine());
    builder.copy_curve(tres_purported, tres_true);

    let mut partial_witness = PartialWitness::new();
    partial_witness.set_targets(
        &txs,
        &xs.into_iter()
            .map(|x| x.try_convert().unwrap())
            .collect::<Vec<_>>(),
    );
    partial_witness.set_point_targets(&tps, &ps);
    let circuit = builder.build();
    let witness = circuit.generate_witness(partial_witness);
    let proof = circuit
        .generate_proof::<Tweedledee>(&witness, &[], true)
        .unwrap();

    let vk = circuit.to_vk();
    verify_proof::<Tweedledum, Tweedledee>(&[], &proof, &[], &vk, true)?;

    Ok(())
}

#[test]
fn test_base_4_sum() -> Result<()> {
    type C = Tweedledee;
    type InnerC = Tweedledum;
    type SF = <C as Curve>::ScalarField;
    type B4 = Base4SumGate<C>;

    let mut builder = CircuitBuilder::<C>::new(128);

    let mut rng = thread_rng();
    let limbs = (0..B4::NUM_LIMBS) //(0..B4::NUM_LIMBS)
        .map(|_| SF::from_canonical_usize(rng.gen_range(0, 4)))
        .collect::<Vec<_>>();
    let x = limbs.iter().fold(SF::ZERO, |acc, &l| acc.quadruple() + l);
    let y = SF::rand();

    let t_y = builder.constant_wire(y);
    let t_x_y = builder.constant_wire(
        x + {
            let mut tmp_y = y;
            (0..B4::NUM_LIMBS).for_each(|_| tmp_y = tmp_y.quadruple());
            tmp_y
        },
    );

    let index = builder.num_gates();
    builder.add_gate_no_constants(B4::new(index));

    builder.copy(
        t_y,
        Target::Wire(Wire {
            gate: index,
            input: B4::WIRE_ACC_OLD,
        }),
    );
    builder.copy(
        t_x_y,
        Target::Wire(Wire {
            gate: index,
            input: B4::WIRE_ACC_NEW,
        }),
    );

    let t_limbs = (0..B4::NUM_LIMBS)
        .map(|i| {
            Target::Wire(Wire {
                gate: index,
                input: B4::wire_limb(i),
            })
        })
        .collect::<Vec<_>>();

    let circuit = builder.build();
    let mut partial_witness = PartialWitness::new();
    partial_witness.set_targets(&t_limbs, &limbs);
    let witness = circuit.generate_witness(partial_witness);
    let proof = circuit.generate_proof::<Tweedledum>(&witness, &[], true)?;
    verify_proof::<C, InnerC>(&[], &proof, &[], &circuit.into(), true)?;

    Ok(())
}

#[test]
fn test_curve_double_gate() -> Result<()> {
    type C = Tweedledee;
    type InnerC = Tweedledum;
    type SF = <C as Curve>::ScalarField;

    let mut builder = CircuitBuilder::<C>::new(128);

    let p = blake_hash_base_field_to_curve::<InnerC>(SF::rand());
    let t_p = builder.constant_affine_point(p);
    let t_double_p_true = builder.constant_affine_point(p.double());

    let t_double_p_purported = builder.curve_double::<InnerC>(t_p);

    builder.copy_curve(t_double_p_purported, t_double_p_true);

    let circuit = builder.build();
    let witness = circuit.generate_witness(PartialWitness::new());
    let proof = circuit.generate_proof::<Tweedledum>(&witness, &[], true)?;
    verify_proof::<C, InnerC>(&[], &proof, &[], &circuit.into(), true)?;

    Ok(())
}
