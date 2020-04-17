use crate::{Circuit, CircuitBuilder, Field, HaloEndomorphismCurve, NUM_CONSTANTS, NUM_WIRES, QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER, Target};

pub struct RecursiveCircuit<F: Field> {
    /// A commitment to each wire polynomial.
    c_wires: Vec<Target>,
    /// A commitment to Z in the context of the permutation argument.
    c_z: Target,
    /// A commitment to the quotient polynomial.
    c_t: Vec<Target>,
    /// A commitment to the aggregate polynomial.
    c_agg: Target,

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

    // A commitment to each wire polynomial.
    let mut c_wires = Vec::with_capacity(NUM_WIRES);
    for _i in 0..NUM_WIRES {
        c_wires.push(builder.add_virtual_target());
    }

    // A commitment to Z, the polynomial used in the permutation argument.
    let c_z = builder.add_virtual_target();

    // A commitment to t, the quotient polynomial, split into several degree-n polynomials.
    let mut c_t = Vec::with_capacity(QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER);
    for _i in 0..QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER {
        c_t.push(builder.add_virtual_target());
    }

    // A commitment to F, a random linear combination of the other polynomials.
    let c_agg = builder.add_virtual_target();

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
        c_agg,
        l_i,
        r_i,
        circuit,
    }
}
