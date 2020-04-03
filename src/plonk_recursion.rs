use crate::{Field, CircuitInput, Circuit, CircuitBuilder, NUM_WIRES, QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER};

struct RecursiveCircuit<F: Field> {
    /// A commitment to each wire polynomial.
    c_wires: Vec<CircuitInput>,
    /// A commitment to Z in the context of the permutation argument.
    c_z: CircuitInput,
    /// A commitment to the quotient polynomial.
    c_t: Vec<CircuitInput>,

    /// L_i in the Halo reduction.
    l_i: Vec<CircuitInput>,
    /// R_i in the Halo reduction.
    r_i: Vec<CircuitInput>,

    circuit: Circuit<F>,
}

fn recursive_verification_circuit<F: Field>(degree_pow: usize) -> RecursiveCircuit<F> {
    let mut builder = CircuitBuilder::new();

    let inner_c_is_noop = builder.add_public_input();
    let inner_c_is_msm = builder.add_public_input();
    let inner_c_is_rescue = builder.add_public_input();
    let inner_c_is_base4sum = builder.add_public_input();
    let inner_c_is_madd = builder.add_public_input();
    let inner_c_const = builder.add_public_input();

    let inner_o_is_noop = builder.add_public_input();
    let inner_o_is_msm = builder.add_public_input();
    let inner_o_is_rescue = builder.add_public_input();
    let inner_o_is_base4sum = builder.add_public_input();
    let inner_o_is_madd = builder.add_public_input();
    let inner_o_t = builder.add_public_inputs(QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER);

    // A commitment to each wire polynomial.
    let mut c_wires = Vec::with_capacity(NUM_WIRES);
    for _i in 0..NUM_WIRES {
        c_wires.push(builder.add_circuit_input());
    }

    // A commitment to Z, the polynomial used in the permutation argument.
    let c_z = builder.add_circuit_input();

    // A commitment to t, the quotient polynomial, split into several degree-n polynomials.
    let mut c_t = Vec::with_capacity(QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER);
    for _i in 0..QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER {
        c_t.push(builder.add_circuit_input());
    }

    let mut l_i = Vec::with_capacity(degree_pow);
    let mut r_i = Vec::with_capacity(degree_pow);
    for _i in 0..degree_pow {
        l_i.push(builder.add_circuit_input());
        r_i.push(builder.add_circuit_input());
    }

    let circuit = builder.build();
    RecursiveCircuit {
        c_wires,
        c_z,
        c_t,
        l_i,
        r_i,
        circuit,
    }
}
