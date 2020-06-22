use plonky::{
    verify_proof_circuit, verify_proof_vk, Circuit, CircuitBuilder, Curve, Field, HaloCurve,
    PartialWitness, Tweedledee, Tweedledum, Witness,
};
use std::time::Instant;

// Make sure it's the same as in `plonk.rs`.
const NUM_WIRES: usize = 9;

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
fn test_proof_trivial_vk() {
    let (circuit, witness) = get_trivial_circuit(<Tweedledee as Curve>::ScalarField::ZERO);
    let proof = circuit
        .generate_proof::<Tweedledum>(witness, &[], true)
        .unwrap();
    let vk = circuit.to_vk();
    assert!(verify_proof_vk::<Tweedledee, Tweedledum>(&[], &proof, &[], &vk,).is_ok());
}

#[test]
fn test_proof_trivial_circuit_many_proofs_vk() {
    let mut old_proofs = Vec::new();
    for _ in 0..10 {
        let (circuit, witness) = get_trivial_circuit(<Tweedledee as Curve>::ScalarField::ZERO);
        let proof = circuit
            .generate_proof::<Tweedledum>(witness.clone(), &[], true)
            .unwrap();
        let vk = circuit.to_vk();
        assert!(verify_proof_vk::<Tweedledee, Tweedledum>(&[], &proof, &[], &vk,).is_ok());
        old_proofs.push(proof.into());
    }
    let (circuit, witness) = get_trivial_circuit(<Tweedledee as Curve>::ScalarField::ZERO);
    let proof = circuit
        .generate_proof::<Tweedledum>(witness.clone(), &old_proofs, true)
        .unwrap();
    let vk = circuit.to_vk();
    assert!(verify_proof_vk::<Tweedledee, Tweedledum>(&[], &proof, &old_proofs, &vk,).is_ok());
}

#[test]
fn test_proof_sum_vk() {
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
        .generate_proof::<Tweedledum>(witness, &[], true)
        .unwrap();
    let vk = circuit.to_vk();
    assert!(verify_proof_vk::<Tweedledee, Tweedledum>(&[], &proof, &[], &vk,).is_ok());
}

#[test]
#[ignore]
fn test_proof_sum_big_vk() {
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
        .generate_proof::<Tweedledum>(witness, &[], true)
        .unwrap();
    dbg!(now.elapsed());
    let vk = circuit.to_vk();
    dbg!(now.elapsed());
    assert!(verify_proof_vk::<Tweedledee, Tweedledum>(&[], &proof, &[], &vk,).is_ok());
    dbg!(now.elapsed());
}

#[test]
fn test_proof_quadratic_vk() {
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
        .generate_proof::<Tweedledum>(witness, &[], true)
        .unwrap();
    let vk = circuit.to_vk();
    assert!(verify_proof_vk::<Tweedledee, Tweedledum>(&[], &proof, &[], &vk,).is_ok());
}

#[test]
fn test_proof_public_input1_vk() {
    // Set public inputs pi1 = 2 and check that pi1 - 2 == 0.
    let mut builder = CircuitBuilder::<Tweedledee>::new(128);
    let pi = builder.stage_public_input();
    builder.route_public_inputs();
    let t1 = pi.routable_target();
    let t2 = builder.constant_wire(<Tweedledee as Curve>::ScalarField::TWO);
    let t3 = builder.sub(t1, t2);
    builder.assert_zero(t3);
    let mut partial_witness = PartialWitness::new();
    partial_witness.set_target(t1, <Tweedledee as Curve>::ScalarField::TWO);
    let circuit = builder.build();
    let witness = circuit.generate_witness(partial_witness);
    let proof = circuit
        .generate_proof::<Tweedledum>(witness, &[], true)
        .unwrap();
    // Check that the public input is set correctly in the proof.
    assert_eq!(
        proof.o_public_inputs[0].o_wires[0],
        <Tweedledee as Curve>::ScalarField::TWO
    );
    let vk = circuit.to_vk();
    assert!(verify_proof_vk::<Tweedledee, Tweedledum>(
        &[<Tweedledee as Curve>::ScalarField::TWO],
        &proof,
        &[],
        &vk,
    )
    .is_ok());
}

#[test]
#[ignore]
fn test_proof_public_input2_vk() {
    // Set many random public inputs
    let mut builder = CircuitBuilder::<Tweedledee>::new(128);
    let pis = (0..200)
        .map(|_| builder.stage_public_input())
        .collect::<Vec<_>>();
    builder.route_public_inputs();
    let tis = pis.iter().map(|p| p.original_wire()).collect::<Vec<_>>();
    let mut partial_witness = PartialWitness::new();
    let values = tis
        .iter()
        .map(|&_t| <Tweedledee as Curve>::ScalarField::rand())
        .collect::<Vec<_>>();
    tis.iter().zip(values.iter()).for_each(|(&t, &v)| {
        partial_witness.set_wire(t, v);
    });
    let circuit = builder.build();
    let witness = circuit.generate_witness(partial_witness);
    let proof = circuit
        .generate_proof::<Tweedledum>(witness, &[], true)
        .unwrap();
    // Check that the public inputs are set correctly in the proof.
    values.iter().enumerate().for_each(|(i, &v)| {
        assert_eq!(
            v,
            proof.o_public_inputs[i / NUM_WIRES].o_wires[i % NUM_WIRES],
        )
    });
    let vk = circuit.to_vk();
    assert!(verify_proof_vk::<Tweedledee, Tweedledum>(&values, &proof, &[], &vk,).is_ok());
}
