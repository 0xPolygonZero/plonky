use crate::{Circuit, CircuitBuilder, Curve, Field, PartialWitness, Tweedledee, Tweedledum};

fn lol() {
    let mut builder = CircuitBuilder::<Tweedledee>::new(128);
    let t1 = builder.add_virtual_target();
    let t2 = builder.add_virtual_target();
    dbg!(builder.num_gates());
    builder.route_public_inputs();
    let s = builder.add(t1, t2);
    builder.assert_zero(s);
    let mut partial_witness = PartialWitness::new();
    let zero = <Tweedledee as Curve>::ScalarField::ZERO;
    partial_witness.set_target(t1, zero);
    partial_witness.set_target(t2, zero);
    dbg!("DONE");
    let circuit = builder.build();
    dbg!("DONE");
    let witness = circuit.generate_witness(partial_witness);
    dbg!("DONE");
    let proof = circuit.generate_proof::<Tweedledum>(witness);
}

fn lal() {
    let mut builder = CircuitBuilder::<Tweedledee>::new(128);
    let pi = builder.stage_public_input();
    let t1 = pi.routable_target();
    dbg!(builder.num_gates());
    builder.route_public_inputs();
    builder.assert_zero(t1);
    let mut partial_witness = PartialWitness::new();
    let zero = <Tweedledee as Curve>::ScalarField::ZERO;
    partial_witness.set_target(t1, zero);
    dbg!("DONE");
    let circuit = builder.build();
    dbg!("DONE");
    let witness = circuit.generate_witness(partial_witness);
    dbg!("DONE");
    let proof = circuit.generate_proof::<Tweedledum>(witness);
}

fn simple() {
let mut builder = CircuitBuilder::<Tweedledee>::new(128);
// let t = builder.add_virtual_target();
// builder.assert_zero(t);
let mut partial_witness = PartialWitness::new();
// partial_witness.set_target(t, <Tweedledee as Curve>::ScalarField::ZERO);
let circuit = builder.build();
let witness = circuit.generate_witness(partial_witness);
let proof = circuit.generate_proof::<Tweedledum>(witness);
}

#[test]
fn test_lol() {
    simple();
}
