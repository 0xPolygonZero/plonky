use crate::{Circuit, CircuitBuilder, Field, HaloEndomorphismCurve, NUM_CONSTANTS, NUM_WIRES, QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER, Target};

/// Wraps a `Circuit` for recursive verification, with an extra field for each input target.
pub struct RecursiveCircuit<F: Field> {
    /// A commitment to each wire polynomial.
    c_wires: Vec<Target>,
    /// A commitment to Z, in the context of the permutation argument.
    c_z: Target,
    /// A commitment to the quotient polynomial.
    c_t: Vec<Target>,

    /// The purported opening of each constant polynomial.
    o_constants: Vec<Target>,
    /// The purported opening of each wire polynomial.
    o_wires: Vec<Target>,
    /// The purported opening of Z, in the context of the permutation argument.
    o_z: Target,
    /// The purported opening of the quotient polynomial.
    o_t: Vec<Target>,

    /// L_i in the Halo reduction.
    l_i: Vec<Target>,
    /// R_i in the Halo reduction.
    r_i: Vec<Target>,

    pub circuit: Circuit<F>,
}

pub fn recursive_verification_circuit<C: HaloEndomorphismCurve>(degree_pow: usize) -> RecursiveCircuit<C::BaseField> {
    let mut builder = CircuitBuilder::<C::BaseField>::new();

    // TODO: Is this actually needed to avoid cyclic dependencies?
    // let inner_c_constants = builder.add_public_inputs(NUM_CONSTANTS);

    let inner_o_constants = builder.add_public_inputs(NUM_CONSTANTS);
    let inner_o_wires = builder.add_public_inputs(NUM_WIRES);
    let inner_o_z = builder.add_public_input();
    let inner_o_t = builder.add_public_inputs(QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER);
    let inner_u = builder.add_public_inputs(degree_pow);
    let inner_pi_hash = builder.add_public_input();

    let c_wires = builder.add_virtual_targets(NUM_WIRES);
    let c_z = builder.add_virtual_target();
    let c_t = builder.add_virtual_targets(QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER);

    let o_constants = builder.add_virtual_targets(NUM_CONSTANTS);
    let o_wires = builder.add_virtual_targets(NUM_WIRES);
    let o_z = builder.add_virtual_target();
    let o_t = builder.add_virtual_targets(QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER);

    // TODO: Verify that each prover polynomial commitment is on the curve.


    // Compute random challenges.
    let (beta, gamma) = builder.rescue_hash_n_to_2(&c_wires);
    let alpha = builder.rescue_hash_n_to_1(&vec![beta, c_z]);
    let zeta = builder.rescue_hash_n_to_1(&[vec![alpha], c_t.clone()].concat());
    let (v, u) = builder.rescue_hash_n_to_2(
        &[vec![zeta], o_constants.clone(), o_wires.clone(), vec![o_z], o_t.clone()].concat());

    let mut l_i = Vec::with_capacity(degree_pow);
    let mut r_i = Vec::with_capacity(degree_pow);
    for _i in 0..degree_pow {
        l_i.push(builder.add_virtual_target());
        r_i.push(builder.add_virtual_target());
    }

    let circuit = builder.build();
    RecursiveCircuit {
        c_wires,
        c_z,
        c_t,
        o_constants,
        o_wires,
        o_z,
        o_t,
        l_i,
        r_i,
        circuit,
    }
}
