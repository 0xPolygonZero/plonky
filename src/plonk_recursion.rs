use crate::{Circuit, CircuitBuilder, Field, HaloEndomorphismCurve, NUM_CONSTANTS, NUM_WIRES, PublicInput, QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER, Target, AffinePointTarget, MsmPart};
use crate::plonk_gates::evaluate_all_constraints_recursively;

/// Wraps a `Circuit` for recursive verification with inputs for the proof data.
pub struct RecursiveCircuit<F: Field> {
    pub circuit: Circuit<F>,
    pub proof: ProofTarget,
}

pub struct ProofTarget {
    /// A commitment to each wire polynomial.
    c_wires: Vec<AffinePointTarget>,
    /// A commitment to Z, in the context of the permutation argument.
    c_plonk_z: AffinePointTarget,
    /// A commitment to the quotient polynomial.
    c_plonk_t: Vec<AffinePointTarget>,

    /// The opening of each polynomial at each `PublicInputGate` index.
    o_public_inputs: Vec<OpeningSetTarget>,
    /// The opening of each polynomial at `zeta`.
    o_local: OpeningSetTarget,
    /// The opening of each polynomial at `g * zeta`.
    o_right: OpeningSetTarget,
    /// The opening of each polynomial at `g^65 * zeta`.
    o_below: OpeningSetTarget,

    // Data for the previous proof in the recursive chain, which hasn't been fully verified.
    inner_beta: PublicInput,
    inner_gamma: PublicInput,
    inner_alpha: PublicInput,
    inner_zeta: PublicInput,
    inner_o_constants: Vec<PublicInput>,
    inner_o_plonk_sigmas: Vec<PublicInput>,
    inner_o_local_wires: Vec<PublicInput>,
    inner_o_right_wires: Vec<PublicInput>,
    inner_o_below_wires: Vec<PublicInput>,
    inner_o_plonk_z_local: PublicInput,
    inner_o_plonk_z_right: PublicInput,
    inner_o_plonk_t: Vec<PublicInput>,
    inner_halo_us: Vec<PublicInput>,

    /// L_i in the Halo reduction.
    halo_l_i: Vec<Target>,
    /// R_i in the Halo reduction.
    halo_r_i: Vec<Target>,
    /// The purported value of G, i.e. <s, G>, in the context of Halo.
    halo_g: Target,
}

impl ProofTarget {
    fn all_opening_sets(&self) -> Vec<OpeningSetTarget> {
        [
            self.o_public_inputs.as_slice(),
            &[self.o_local.clone(), self.o_right.clone(), self.o_below.clone()],
        ].concat()
    }

    fn all_opening_targets(&self) -> Vec<Target> {
        let targets_2d: Vec<Vec<Target>> = self.all_opening_sets().into_iter()
            .map(|set| set.to_vec())
            .collect();
        targets_2d.concat()
    }
}

#[derive(Clone)]
struct OpeningSetTarget {
    /// The purported opening of each constant polynomial.
    o_constants: Vec<Target>,
    /// The purported opening of each S_sigma polynomial in the context of Plonk's permutation argument.
    o_plonk_sigmas: Vec<Target>,
    /// The purported opening of each wire polynomial.
    o_wires: Vec<Target>,
    /// The purported opening of `Z`.
    o_plonk_z: Target,
    /// The purported opening of `t`.
    o_plonk_t: Vec<Target>,
}

impl OpeningSetTarget {
    fn to_vec(&self) -> Vec<Target> {
        [
            self.o_constants.as_slice(),
            self.o_plonk_sigmas.as_slice(),
            self.o_wires.as_slice(),
            &[self.o_plonk_z],
            self.o_plonk_t.as_slice(),
        ].concat()
    }
}

pub fn recursive_verification_circuit<C: HaloEndomorphismCurve>(
    degree_pow: usize,
) -> RecursiveCircuit<C::BaseField> {
    let mut builder = CircuitBuilder::<C::BaseField>::new();
    let proof = ProofTarget {
        c_wires: builder.add_virtual_point_targets(NUM_WIRES),
        c_plonk_z: builder.add_virtual_point_target(),
        c_plonk_t: builder.add_virtual_point_targets(QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER),
        o_public_inputs: todo!(),
        o_local: make_opening_set(&mut builder),
        o_right: make_opening_set(&mut builder),
        o_below: make_opening_set(&mut builder),
        inner_beta: builder.stage_public_input(),
        inner_gamma: builder.stage_public_input(),
        inner_alpha: builder.stage_public_input(),
        inner_zeta: builder.stage_public_input(),
        inner_o_constants: builder.stage_public_inputs(NUM_CONSTANTS),
        inner_o_plonk_sigmas: builder.stage_public_inputs(NUM_WIRES),
        inner_o_local_wires: builder.stage_public_inputs(NUM_WIRES),
        inner_o_right_wires: builder.stage_public_inputs(NUM_WIRES),
        inner_o_below_wires: builder.stage_public_inputs(NUM_WIRES),
        inner_o_plonk_z_local: builder.stage_public_input(),
        inner_o_plonk_z_right: builder.stage_public_input(),
        inner_o_plonk_t: builder.stage_public_inputs(QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER),
        inner_halo_us: builder.stage_public_inputs(degree_pow),
        halo_l_i: builder.add_virtual_targets(degree_pow),
        halo_r_i: builder.add_virtual_targets(degree_pow),
        halo_g: builder.add_virtual_target(),
    };
    builder.route_public_inputs();

    verify_assumptions::<C>(&mut builder, degree_pow, &proof);

    // TODO: Verify that each prover polynomial commitment is on the curve.

    // Compute random challenges.
    let (beta, gamma) = builder.rescue_hash_n_to_2(&flatten_points(&proof.c_wires));
    let alpha = builder.rescue_hash_n_to_1(&[vec![beta], proof.c_plonk_z.to_vec()].concat());
    let zeta = builder.rescue_hash_n_to_1(&[vec![alpha], flatten_points(&proof.c_plonk_t)].concat());
    let (v, u) = builder.rescue_hash_n_to_2(&[
        vec![zeta],
        proof.all_opening_targets(),
    ].concat());

    verify_ipas::<C>(&mut builder, &proof, u, v);

    let circuit = builder.build();
    RecursiveCircuit { circuit, proof }
}

fn flatten_points(points: &[AffinePointTarget]) -> Vec<Target> {
    let coordinate_pairs: Vec<Vec<Target>> = points.iter()
        .map(|p| p.to_vec())
        .collect();
    coordinate_pairs.concat()
}

fn verify_ipas<C: HaloEndomorphismCurve>(
    builder: &mut CircuitBuilder<C::BaseField>,
    proof: &ProofTarget,
    u: Target,
    v: Target,
) {
    // Reduce all polynomial commitments to a single one.
    let c_constants = todo!();
    let c_plonk_sigmas = todo!();
    let c_all: Vec<AffinePointTarget> = [
        c_constants,
        c_plonk_sigmas,
        proof.c_wires.clone(),
        vec![proof.c_plonk_z],
        proof.c_plonk_t.clone(),
    ].concat();
    let mut c_reduction_parts = Vec::new();
    let powers_of_u = powers(builder, u, c_all.len());
    for (&c, &power) in c_all.iter().zip(powers_of_u.iter()) {
        // TODO: Need to split into a mix of base 2 and 4, and verify the weighted sum.
        let scalar_bits = builder.split_binary(power, 128);
        c_reduction_parts.push(MsmPart { scalar_bits, addend: c });
    }
    let c_reduction_msm_result = builder.curve_msm_endo::<C>(&c_reduction_parts);
    let actual_scalars = c_reduction_msm_result.actual_scalars;
    let c_reduction = c_reduction_msm_result.msm_result;

    // For each opening location, we do a similar reduction, using the actual scalars above.
    let o_all = proof.all_opening_sets();
    let o_reductions: Vec<Target> = o_all.iter()
        .map(|opening_set| reduce_with_coefficients(
            builder, &opening_set.to_vec(), &actual_scalars))
        .collect();

    // TODO: Verify reduced IPA.
}

fn make_opening_set<F: Field>(builder: &mut CircuitBuilder<F>) -> OpeningSetTarget {
    OpeningSetTarget {
        o_constants: builder.add_virtual_targets(NUM_CONSTANTS),
        o_plonk_sigmas: builder.add_virtual_targets(NUM_WIRES),
        o_wires: builder.add_virtual_targets(NUM_WIRES),
        o_plonk_z: builder.add_virtual_target(),
        o_plonk_t: builder.add_virtual_targets(QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER),
    }
}

/// In our recursion scheme, to avoid non-native field arithmetic, each proof in a recursive chain
/// only partially verifies its inner proof. It outputs various challenges and openings, and the
/// following proof is expected to verify constraints upon that data. This function performs those
/// final verification steps.
fn verify_assumptions<C: HaloEndomorphismCurve>(
    builder: &mut CircuitBuilder<C::BaseField>,
    degree_pow: usize,
    proof: &ProofTarget,
) {
    let degree = 1 << degree_pow;
    let degree_f = C::BaseField::from_canonical_usize(degree);
    let degree_wire = builder.constant_wire(degree_f);

    let one = builder.one_wire();

    // Convert inner proof data from `PublicInput`s to `Target`s.
    let o_constants: Vec<Target> = proof.inner_o_constants.iter().map(PublicInput::routable_target).collect();
    let o_sigmas: Vec<Target> = proof.inner_o_plonk_sigmas.iter().map(PublicInput::routable_target).collect();
    let o_local_wires: Vec<Target> = proof.inner_o_local_wires.iter().map(PublicInput::routable_target).collect();
    let o_right_wires: Vec<Target> = proof.inner_o_right_wires.iter().map(PublicInput::routable_target).collect();
    let o_below_wires: Vec<Target> = proof.inner_o_below_wires.iter().map(PublicInput::routable_target).collect();
    let beta = proof.inner_beta.routable_target();
    let gamma = proof.inner_gamma.routable_target();
    let alpha = proof.inner_alpha.routable_target();
    let zeta = proof.inner_zeta.routable_target();
    let o_z_local = proof.inner_o_plonk_z_local.routable_target();
    let o_z_right = proof.inner_o_plonk_z_right.routable_target();

    // Evaluate zeta^degree.
    let mut zeta_power_d = zeta;
    for _i in 0..degree_pow {
        zeta_power_d = builder.double(zeta_power_d);
    }

    // Evaluate Z_H(zeta) = zeta^degree - 1.
    let zero_eval = builder.sub(zeta_power_d, one);

    // Evaluate L_1(zeta) = (zeta^degree - 1) / (degree * (zeta - 1)).
    let zeta_minus_one = builder.sub(zeta, one);
    let lagrange_1_eval_denominator = builder.mul(degree_wire, zeta_minus_one);
    let lagrange_1_eval = builder.div(zero_eval, lagrange_1_eval_denominator);

    // Compute Z(zeta) f'(zeta) - Z(g * zeta) g'(zeta), which should vanish on H.
    let mut f_prime = one;
    let mut g_prime = one;
    let quadratic_nonresidues = C::BaseField::generate_quadratic_nonresidues(NUM_WIRES - 1);
    for i in 0..NUM_WIRES {
        let s_id = if i == 0 {
            zeta
        } else {
            let k_i = builder.constant_wire(quadratic_nonresidues[i - 1]);
            builder.mul(k_i, zeta)
        };
        let beta_s_id = builder.mul(beta, s_id);
        let beta_s_sigma = builder.mul(beta, o_sigmas[i]);
        let f_prime_part = builder.add_many(&[o_local_wires[i], beta_s_id, gamma]);
        let g_prime_part = builder.add_many(&[o_local_wires[i], beta_s_sigma, gamma]);
        f_prime = builder.mul(f_prime, f_prime_part);
        g_prime = builder.mul(g_prime, g_prime_part);
    }
    let z_f_prime = builder.mul(o_z_local, f_prime);
    let z_shifted_g_prime = builder.mul(o_z_right, g_prime);
    let vanishing_v_shift_term = builder.sub(z_f_prime, z_shifted_g_prime);

    // Evaluate the function which is supposed to vanish on H. It is a sum of several terms which
    // should vanish, each weighted by a different power of alpha.
    let o_z_minus_1 = builder.sub(o_z_local, one);
    let vanishing_z_1_term = builder.mul(o_z_minus_1, lagrange_1_eval);
    let constraint_terms = evaluate_all_constraints_recursively::<C>(
        builder, &o_constants, &o_local_wires, &o_right_wires, &o_below_wires);
    let vanishing_terms = [
        vec![vanishing_z_1_term],
        vec![vanishing_v_shift_term],
        constraint_terms
    ].concat();
    let vanishing_eval = reduce_with_powers(builder, &vanishing_terms, alpha);

    // Evaluate the quotient polynomial, and assert that it matches the prover's opening.
    let quotient_eval = builder.div(vanishing_eval, zero_eval);
    let inner_o_plonk_t_targets: Vec<Target> =
        proof.inner_o_plonk_t.iter()
            .map(|pi| pi.routable_target())
            .collect();
    let inner_o_plonk_t_eval = eval_composite_poly(builder, &inner_o_plonk_t_targets, zeta_power_d);
    builder.copy(quotient_eval, inner_o_plonk_t_eval);
}

/// Computes a sum of terms weighted by the given coefficients.
fn reduce_with_coefficients<F: Field>(
    builder: &mut CircuitBuilder<F>,
    terms: &[Target],
    coefficients: &[Target],
) -> Target {
    let mut reduction = builder.zero_wire();
    for (i, &term) in terms.iter().enumerate() {
        reduction = builder.mul_add(coefficients[i], term, reduction);
    }
    reduction
}

/// Computes a sum of terms weighted by powers of alpha.
fn reduce_with_powers<F: Field>(
    builder: &mut CircuitBuilder<F>,
    terms: &[Target],
    alpha: Target,
) -> Target {
    let mut sum = builder.zero_wire();
    for &term in terms.iter().rev() {
        sum = builder.mul_add(sum, alpha, term);
    }
    sum
}

/// Compute `[x^0, x^1, ..., x^(n - 1)]`.
fn powers<F: Field>(builder: &mut CircuitBuilder<F>, x: Target, n: usize) -> Vec<Target> {
    let mut powers = Vec::new();
    let mut current = builder.one_wire();
    for i in 0..n {
        if i != 0 {
            current = builder.mul(current, x);
        }
        powers.push(current);
    }
    powers
}

/// In Plonk, some polynomials are broken up into degree-d components. Given an evaluation of each
/// component at some point zeta, this function evaluates the composite polynomial at zeta.
fn eval_composite_poly<F: Field>(
    builder: &mut CircuitBuilder<F>,
    component_evals: &[Target],
    zeta_power_d: Target,
) -> Target {
    reduce_with_powers(builder, component_evals, zeta_power_d)
}

/// Evaluate `g(X, {u_i})` as defined in the Halo paper.
fn halo_g<F: Field>(builder: &mut CircuitBuilder<F>, x: Target, us: &[Target]) -> Target {
    let mut product = builder.one_wire();
    let mut x_power = x;
    for &u_i in us {
        let u_i_inv = builder.inv(u_i);
        let term = builder.mul_add(u_i_inv, x_power, u_i);
        product = builder.mul(product, term);
        x_power = builder.double(x_power);
    }
    product
}
