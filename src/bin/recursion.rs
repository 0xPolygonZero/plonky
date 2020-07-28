use anyhow::Result;
use std::time::Instant;

use plonky::{recursive_verification_circuit, verify_proof, BufferGate, Circuit, CircuitBuilder, PartialWitness, Tweedledee, Tweedledum, OldProof, TweedledumBase};

const INNER_PROOF_DEGREE_POW: usize = 14;
const INNER_PROOF_DEGREE: usize = 1 << INNER_PROOF_DEGREE_POW;
const SECURITY_BITS: usize = 128;

fn main() -> Result<()> {
    println!("Generating inner circuit");
    let start = Instant::now();
    let inner_circuit: Circuit<Tweedledum> = generate_trivial_circuit();
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

    println!("Generating inner witness");
    let start = Instant::now();
    let inner_witness = inner_circuit.generate_witness(PartialWitness::new());
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

    println!("Generating inner proof");
    let start = Instant::now();
    let inner_old_proofs: Vec<OldProof<Tweedledum>> = vec![];
    let inner_proof = inner_circuit
        .generate_proof::<Tweedledee>(inner_witness, &inner_old_proofs, true, true)
        .unwrap();
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

    println!("Verifying inner proof");
    let start = Instant::now();
    let inner_vk = inner_circuit.to_vk();
    verify_proof::<Tweedledum, Tweedledee>(&[], &inner_proof, &inner_old_proofs, &inner_vk, true)?;
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

    println!("Generating recursion circuit...");
    let recursion_old_proofs: Vec<OldProof<Tweedledee>> = vec![];
    let start = Instant::now();
    let recursion_circuit = recursive_verification_circuit::<Tweedledee, Tweedledum>(
        // INNER_PROOF_DEGREE_POW,
        inner_proof.halo_l.len(),
        SECURITY_BITS,
        0,
        recursion_old_proofs.len(),
    );
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!("Gate count: {}", recursion_circuit.circuit.degree());
    println!();

    // Populate inputs.
    let mut recursion_inputs: PartialWitness<TweedledumBase> = PartialWitness::new();
    if let Err(e) = recursion_circuit
        .proof
        .populate_witness(&mut recursion_inputs, inner_proof.clone())
    {
        panic!("Failed to populate inputs: {:?}", e);
    }

    // Populate old proofs inputs.
    recursion_old_proofs
        .iter()
        .zip(recursion_circuit.old_proofs.iter())
        .for_each(|(p, pt)| {
            pt.populate_witness(&mut recursion_inputs, p).unwrap();
        });

    println!("Generating recursion witness...");
    let start = Instant::now();
    let recursion_witness = recursion_circuit.circuit.generate_witness(recursion_inputs);
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

    println!("Generating recursion proof...");
    let start = Instant::now();
    let proof = recursion_circuit
        .circuit
        // .generate_proof::<Tweedledum>(recursion_witness, &[inner_proof.into()], true)
        .generate_proof::<Tweedledum>(recursion_witness, &[], true, true)
        .unwrap();
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    println!();

    println!("Verifying proof...");
    let start = Instant::now();
    let pis = proof.get_public_inputs(recursion_circuit.circuit.num_public_inputs);
    println!(
        "Number of public inputs: {}",
        recursion_circuit.circuit.num_public_inputs
    );
    let vk = recursion_circuit.circuit.to_vk();
    verify_proof::<Tweedledee, Tweedledum>(&pis.unwrap(), &proof, &[], &vk, true)?;
    println!("Finished in {}s", start.elapsed().as_secs_f64());
    Ok(())
}

fn generate_trivial_circuit() -> Circuit<Tweedledum> {
    let mut builder = CircuitBuilder::new(SECURITY_BITS);
    builder.route_public_inputs();
    while builder.num_gates() < INNER_PROOF_DEGREE - 3 {
        builder.add_gate_no_constants(BufferGate::new(builder.num_gates()));
    }
    builder.build()
}
