use anyhow::Result;
use plonky::{recursive_verification_circuit, verify_proof, CircuitBuilder, Curve, Field, PartialWitness, Tweedledee, Tweedledum};

#[test]
fn test_proof_trivial_recursive() -> Result<()> {
    let mut builder = CircuitBuilder::<Tweedledee>::new(128);
    let t = builder.constant_wire(<Tweedledee as Curve>::ScalarField::ZERO);
    let mut partial_witness = PartialWitness::new();
    partial_witness.set_target(t, <Tweedledee as Curve>::ScalarField::ZERO);
    let inner_circuit = builder.build();
    let witness = inner_circuit.generate_witness(partial_witness);
    let inner_proof = inner_circuit
        .generate_proof::<Tweedledum>(witness, &[], true)
        .unwrap();
    let inner_vk = inner_circuit.to_vk();
    verify_proof::<Tweedledee, Tweedledum>(&[], &inner_proof, &[], &inner_vk, true)?;

    let recursion_circuit = recursive_verification_circuit::<Tweedledum, Tweedledee>(
        inner_circuit.degree_pow(),
        128,
        0,
        0,
    );
    let mut recursion_inputs = PartialWitness::new();
    if let Err(e) = recursion_circuit
        .proof
        .populate_witness(&mut recursion_inputs, inner_proof)
    {
        panic!("Failed to populate inputs: {:?}", e);
    }
    let recursion_witness = recursion_circuit.circuit.generate_witness(recursion_inputs);
    let proof = recursion_circuit
        .circuit
        .generate_proof::<Tweedledee>(recursion_witness, &[], true)
        .unwrap();
    let vk = recursion_circuit.circuit.to_vk();
    verify_proof::<Tweedledum, Tweedledee>(&[], &proof, &[], &vk, true)?;

    Ok(())
}
