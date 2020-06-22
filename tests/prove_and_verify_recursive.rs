use plonky::{recursive_verification_circuit, verify_proof_circuit, CircuitBuilder, Curve, Field, PartialWitness, Tweedledee, Tweedledum};
use std::time::Instant;

#[test]
fn test_proof_trivial_recursive() {
    let mut builder = CircuitBuilder::<Tweedledee>::new(128);
    let t = builder.constant_wire(<Tweedledee as Curve>::ScalarField::ZERO);
    // builder.assert_zero(t);
    let mut partial_witness = PartialWitness::new();
    partial_witness.set_target(t, <Tweedledee as Curve>::ScalarField::ZERO);
    let circuit = builder.build();
    let witness = circuit.generate_witness(partial_witness);
    let proof = circuit.generate_proof::<Tweedledum>(witness, &[], true).unwrap();
    assert!(verify_proof_circuit::<Tweedledee, Tweedledum>(&[], &proof, &[], &circuit, true).is_ok());

    let recursion_circuit =
        recursive_verification_circuit::<Tweedledum, Tweedledee>(circuit.degree_pow(), 128);
    let mut recursion_inputs = PartialWitness::new();
    if let Err(e) = recursion_circuit
        .proof
        .populate_witness(&mut recursion_inputs, proof)
    {
        panic!("Failed to populate inputs: {:?}", e);
    }
    let recursion_witness = recursion_circuit.circuit.generate_witness(recursion_inputs);
    let proof = recursion_circuit
        .circuit
        .generate_proof::<Tweedledee>(recursion_witness, &[], true)
        .unwrap();
    assert!(verify_proof_circuit::<Tweedledum, Tweedledee>(
        &[],
        &proof,
        &[],
        &recursion_circuit.circuit,
        true
    )
    .is_ok());
}
