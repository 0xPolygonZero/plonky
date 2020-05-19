use crate::{AffinePointTarget, Circuit, CircuitBuilder, CurveMulOp, Field, HaloEndomorphismCurve, NUM_CONSTANTS, NUM_ROUTED_WIRES, NUM_WIRES, OpeningSetTarget, ProofTarget, PublicInput, QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER, Target, Curve};
use crate::plonk_gates::evaluate_all_constraints_recursively;
use crate::util::ceil_div_usize;

/// Wraps a `Circuit` for recursive verification with inputs for the proof data.
pub struct RecursiveCircuit<C: Curve> {
    pub circuit: Circuit<C>,
    pub proof: ProofTarget,
}

/// Public inputs of the recursive circuit. This contains data for the inner proof which is needed
/// to complete verification of it.
struct RecursionPublicInputs {
    beta: PublicInput,
    gamma: PublicInput,
    alpha: PublicInput,
    zeta: PublicInput,
    o_constants: Vec<PublicInput>,
    o_plonk_sigmas: Vec<PublicInput>,
    o_local_wires: Vec<PublicInput>,
    o_right_wires: Vec<PublicInput>,
    o_below_wires: Vec<PublicInput>,
    o_plonk_z_local: PublicInput,
    o_plonk_z_right: PublicInput,
    o_plonk_t: Vec<PublicInput>,
    halo_u_l: Vec<PublicInput>,
    halo_u_r: Vec<PublicInput>,
}

impl RecursionPublicInputs {
    fn to_vec(&self) -> Vec<PublicInput> {
        [
            &[self.beta, self.gamma, self.alpha, self.zeta],
            self.o_constants.as_slice(),
            self.o_plonk_sigmas.as_slice(),
            self.o_local_wires.as_slice(),
            self.o_right_wires.as_slice(),
            self.o_below_wires.as_slice(),
            &[self.o_plonk_z_local, self.o_plonk_z_right],
            self.o_plonk_t.as_slice(),
            self.halo_u_l.as_slice(),
            self.halo_u_r.as_slice(),
        ].concat()
    }
}

/// The number of `PublicInputGate`s needed to route the given number of public inputs.
fn num_public_input_gates(num_public_inputs: usize) -> usize {
    ceil_div_usize(num_public_inputs, NUM_WIRES)
}

pub fn recursive_verification_circuit<C: Curve, InnerC: HaloEndomorphismCurve<BaseField=C::ScalarField>>(
    degree_pow: usize,
    security_bits: usize,
) -> RecursiveCircuit<C> {
    let mut builder = CircuitBuilder::<C>::new(security_bits);
    let public_inputs = RecursionPublicInputs {
        beta: builder.stage_public_input(),
        gamma: builder.stage_public_input(),
        alpha: builder.stage_public_input(),
        zeta: builder.stage_public_input(),
        o_constants: builder.stage_public_inputs(NUM_CONSTANTS),
        o_plonk_sigmas: builder.stage_public_inputs(NUM_ROUTED_WIRES),
        o_local_wires: builder.stage_public_inputs(NUM_WIRES),
        o_right_wires: builder.stage_public_inputs(NUM_WIRES),
        o_below_wires: builder.stage_public_inputs(NUM_WIRES),
        o_plonk_z_local: builder.stage_public_input(),
        o_plonk_z_right: builder.stage_public_input(),
        o_plonk_t: builder.stage_public_inputs(QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER),
        halo_u_l: builder.stage_public_inputs(degree_pow),
        halo_u_r: builder.stage_public_inputs(degree_pow),
    };

    let num_public_inputs = public_inputs.to_vec().len();
    let num_public_input_gates = num_public_input_gates(num_public_inputs);

    let proof = ProofTarget {
        c_wires: builder.add_virtual_point_targets(NUM_WIRES),
        c_plonk_z: builder.add_virtual_point_target(),
        c_plonk_t: builder.add_virtual_point_targets(QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER),
        o_public_inputs: make_opening_sets(&mut builder, num_public_input_gates),
        o_local: make_opening_set(&mut builder),
        o_right: make_opening_set(&mut builder),
        o_below: make_opening_set(&mut builder),
        halo_l_i: builder.add_virtual_point_targets(degree_pow),
        halo_r_i: builder.add_virtual_point_targets(degree_pow),
        halo_g: builder.add_virtual_point_target(),
    };
    builder.route_public_inputs();

    // Flatten the list of public input openings.
    let o_public_inputs: Vec<Target> = proof.o_public_inputs.iter().cloned()
        .flat_map(|opening_set| opening_set.o_wires)
        .collect();

    verify_assumptions::<C, InnerC>(&mut builder, degree_pow, &public_inputs, &o_public_inputs);

    // TODO: Verify that each prover polynomial commitment is on the curve.
    // Can call curve_assert_valid.

    // Compute random challenges.
    let (beta, gamma) = builder.rescue_hash_n_to_2(&flatten_points(&proof.c_wires));
    let alpha = builder.rescue_hash_n_to_1(&[vec![beta], proof.c_plonk_z.to_vec()].concat());
    let zeta = builder.rescue_hash_n_to_1(&[vec![alpha], flatten_points(&proof.c_plonk_t)].concat());
    let (v, u, x) = builder.rescue_hash_n_to_3(&[
        vec![zeta],
        proof.all_opening_targets(),
    ].concat());

    // Compute IPA challenges.
    let mut transcript_state = v;
    let mut ipa_challenges = Vec::new();
    for i in 0..degree_pow {
        let u_i = builder.rescue_hash_n_to_1(&[
            vec![transcript_state],
            proof.halo_l_i[i].to_vec(),
            proof.halo_r_i[i].to_vec(),
        ].concat());
        ipa_challenges.push(u_i);
        transcript_state = u_i;
    }

    let (u_l, u_r) = verify_all_ipas::<C, InnerC>(&mut builder, &proof, u, v, x, ipa_challenges);

    // "Outputs" data relating to assumption which still need to be verified by the next proof.
    builder.copy(public_inputs.beta.routable_target(), beta);
    builder.copy(public_inputs.gamma.routable_target(), gamma);
    builder.copy(public_inputs.alpha.routable_target(), alpha);
    builder.copy(public_inputs.zeta.routable_target(), zeta);
    for i in 0..NUM_CONSTANTS {
        builder.copy(public_inputs.o_constants[i].routable_target(), proof.o_local.o_constants[i]);
    }
    for i in 0..NUM_ROUTED_WIRES {
        builder.copy(public_inputs.o_plonk_sigmas[i].routable_target(), proof.o_local.o_plonk_sigmas[i]);
    }
    for i in 0..NUM_WIRES {
        builder.copy(public_inputs.o_local_wires[i].routable_target(), proof.o_local.o_wires[i]);
        builder.copy(public_inputs.o_right_wires[i].routable_target(), proof.o_right.o_wires[i]);
        builder.copy(public_inputs.o_below_wires[i].routable_target(), proof.o_below.o_wires[i]);
    }
    builder.copy(public_inputs.o_plonk_z_local.routable_target(), proof.o_local.o_plonk_z);
    builder.copy(public_inputs.o_plonk_z_right.routable_target(), proof.o_right.o_plonk_z);
    for i in 0..degree_pow {
        builder.copy(public_inputs.halo_u_l[i].routable_target(), u_l[i]);
        builder.copy(public_inputs.halo_u_r[i].routable_target(), u_r[i]);
    }

    let circuit = builder.build();
    RecursiveCircuit { circuit, proof }
}

fn flatten_points(points: &[AffinePointTarget]) -> Vec<Target> {
    let coordinate_pairs: Vec<Vec<Target>> = points.iter()
        .map(|p| p.to_vec())
        .collect();
    coordinate_pairs.concat()
}

/// Verify all IPAs in the given proof. Return `(u_l, u_r)`, which roughly correspond to `u` and
/// `u^{-1}` in the Halo paper, respectively.
fn verify_all_ipas<C: Curve, InnerC: HaloEndomorphismCurve<BaseField=C::ScalarField>>(
    builder: &mut CircuitBuilder<C>,
    proof: &ProofTarget,
    u: Target,
    v: Target,
    x: Target,
    ipa_challenges: Vec<Target>,
) -> (Vec<Target>, Vec<Target>) {
    // Reduce all polynomial commitments to a single one, i.e. a random combination of them.
    // TODO: Configure the actual constants and permutations of whatever circuit we wish to verify.
    // For now, we use a dummy point for each of those polynomial commitments.
    let dummy_point = builder.constant_affine_point(InnerC::GENERATOR_AFFINE);
    let c_constants = vec![dummy_point; NUM_CONSTANTS];
    let c_plonk_sigmas = vec![dummy_point; NUM_ROUTED_WIRES];
    let c_all: Vec<AffinePointTarget> = [
        c_constants,
        c_plonk_sigmas,
        proof.c_wires.clone(),
        vec![proof.c_plonk_z],
        proof.c_plonk_t.clone(),
    ].concat();
    let mut c_reduction_muls = Vec::new();
    let powers_of_u = powers(builder, u, c_all.len());
    for (&c, &power) in c_all.iter().zip(powers_of_u.iter()) {
        let mul = CurveMulOp { scalar: power, point: c };
        c_reduction_muls.push(mul);
    }
    let c_reduction_msm_result = builder.curve_msm_endo::<InnerC>(&c_reduction_muls);
    let actual_scalars = c_reduction_msm_result.actual_scalars;
    let c_reduction = c_reduction_msm_result.msm_result;

    // For each opening set, we do a similar reduction, using the actual scalars above.
    let opening_set_reductions: Vec<Target> = proof.all_opening_sets().iter()
        .map(|opening_set| reduce_with_coefficients(
            builder, &opening_set.to_vec(), &actual_scalars))
        .collect();

    // Then, we reduce the above opening set reductions to a single value.
    let reduced_opening = reduce_with_powers(builder, &opening_set_reductions, v);

    verify_ipa::<C, InnerC>(builder, proof, c_reduction, reduced_opening, x, ipa_challenges)
}

/// Verify the final IPA. Return `(u_l, u_r)`, which roughly correspond to `u` and `u^{-1}` in the
/// Halo paper, respectively.
fn verify_ipa<C: Curve, InnerC: HaloEndomorphismCurve<BaseField=C::ScalarField>>(
    builder: &mut CircuitBuilder<C>,
    proof: &ProofTarget,
    p: AffinePointTarget,
    c: Target,
    x: Target,
    ipa_challenges: Vec<Target>,
) -> (Vec<Target>, Vec<Target>) {
    // Now we begin IPA verification by computing P' and u' as in Protocol 1 of Bulletproofs.
    // In Protocol 1 we compute u' = [x] u, but we leverage to endomorphism, instead computing
    // u' = [n(x)] u.
    let u = builder.constant_affine_point(InnerC::GENERATOR_AFFINE);
    let u_prime = builder.curve_mul_endo::<InnerC>(CurveMulOp { scalar: x, point: u }).mul_result;

    // Compute [c] [n(x)] u = [c] u'.
    let u_n_x_c = builder.curve_mul::<InnerC>(CurveMulOp { scalar: c, point: u_prime });
    let p_prime = builder.curve_add::<InnerC>(p, u_n_x_c);

    // Compute Q as defined in the Halo paper.
    let mut q_muls = Vec::new();
    for i in 0..proof.halo_l_i.len() {
        let l_i = proof.halo_l_i[i];
        let r_i = proof.halo_r_i[i];
        let l_challenge = ipa_challenges[i];
        let r_challenge = builder.inv(l_challenge);
        q_muls.push(CurveMulOp { scalar: l_challenge, point: l_i });
        q_muls.push(CurveMulOp { scalar: r_challenge, point: r_i });
    }
    let q_msm_result = builder.curve_msm_endo::<InnerC>(&q_muls);
    let q = builder.curve_add::<InnerC>(p_prime, q_msm_result.msm_result);

    // TODO: Perform ZK opening protocol.

    // Take the square roots of the actual scalars as u_l and u_r, which will be used elsewhere
    // in the protocol.
    let mut u_l = Vec::new();
    let mut u_r = Vec::new();
    for (i, &scalar) in q_msm_result.actual_scalars.iter().enumerate() {
        let sqrt = builder.deterministic_square_root(scalar);
        if i % 2 == 0 {
            u_l.push(sqrt);
        } else {
            u_r.push(sqrt);
        }
    }

    (u_l, u_r)
}

fn make_opening_set<C: Curve>(builder: &mut CircuitBuilder<C>) -> OpeningSetTarget {
    OpeningSetTarget {
        o_constants: builder.add_virtual_targets(NUM_CONSTANTS),
        o_plonk_sigmas: builder.add_virtual_targets(NUM_ROUTED_WIRES),
        o_wires: builder.add_virtual_targets(NUM_WIRES),
        o_plonk_z: builder.add_virtual_target(),
        o_plonk_t: builder.add_virtual_targets(QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER),
    }
}

fn make_opening_sets<C: Curve>(builder: &mut CircuitBuilder<C>, n: usize) -> Vec<OpeningSetTarget> {
    (0..n).map(|i| make_opening_set(builder)).collect()
}

/// In our recursion scheme, to avoid non-native field arithmetic, each proof in a recursive chain
/// only partially verifies its inner proof. It outputs various challenges and openings, and the
/// following proof is expected to verify constraints upon that data. This function performs those
/// final verification steps.
fn verify_assumptions<C: Curve, InnerC: HaloEndomorphismCurve<BaseField=C::ScalarField>>(
    builder: &mut CircuitBuilder<C>,
    degree_pow: usize,
    public_inputs: &RecursionPublicInputs,
    o_public_inputs: &[Target],
) {
    // TODO: Still need to verify evaluation of g().

    let degree = 1 << degree_pow;
    let degree_f = C::ScalarField::from_canonical_usize(degree);
    let degree_wire = builder.constant_wire(degree_f);

    let one = builder.one_wire();

    // Convert inner proof data from `PublicInput`s to `Target`s.
    let o_constants: Vec<Target> = public_inputs.o_constants.iter().map(|pi| o_public_inputs[pi.index]).collect();
    let o_sigmas: Vec<Target> = public_inputs.o_plonk_sigmas.iter().map(|pi| o_public_inputs[pi.index]).collect();
    let o_local_wires: Vec<Target> = public_inputs.o_local_wires.iter().map(|pi| o_public_inputs[pi.index]).collect();
    let o_right_wires: Vec<Target> = public_inputs.o_right_wires.iter().map(|pi| o_public_inputs[pi.index]).collect();
    let o_below_wires: Vec<Target> = public_inputs.o_below_wires.iter().map(|pi| o_public_inputs[pi.index]).collect();
    let beta = o_public_inputs[public_inputs.beta.index];
    let gamma = o_public_inputs[public_inputs.gamma.index];
    let alpha = o_public_inputs[public_inputs.alpha.index];
    let zeta = o_public_inputs[public_inputs.zeta.index];
    let o_z_local = o_public_inputs[public_inputs.o_plonk_z_local.index];
    let o_z_right = o_public_inputs[public_inputs.o_plonk_z_right.index];

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
    let quadratic_nonresidues = C::ScalarField::generate_quadratic_nonresidues(NUM_ROUTED_WIRES - 1);
    for i in 0..NUM_ROUTED_WIRES {
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
    let constraint_terms = evaluate_all_constraints_recursively::<C, InnerC>(
        builder, &o_constants, &o_local_wires, &o_right_wires, &o_below_wires);
    let vanishing_terms = [
        vec![vanishing_z_1_term],
        vec![vanishing_v_shift_term],
        constraint_terms
    ].concat();
    let vanishing_eval = reduce_with_powers(builder, &vanishing_terms, alpha);

    // Evaluate the quotient polynomial, and assert that it matches the prover's opening.
    let quotient_eval = builder.div(vanishing_eval, zero_eval);
    let t_components: Vec<Target> =
        public_inputs.o_plonk_t.iter()
            .map(|pi| o_public_inputs[pi.index])
            .collect();
    let o_plonk_t_eval = eval_composite_poly(builder, &t_components, zeta_power_d);
    builder.copy(quotient_eval, o_plonk_t_eval);
}

/// Computes a sum of terms weighted by the given coefficients.
fn reduce_with_coefficients<C: Curve>(
    builder: &mut CircuitBuilder<C>,
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
fn reduce_with_powers<C: Curve>(
    builder: &mut CircuitBuilder<C>,
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
fn powers<C: Curve>(builder: &mut CircuitBuilder<C>, x: Target, n: usize) -> Vec<Target> {
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
fn eval_composite_poly<C: Curve>(
    builder: &mut CircuitBuilder<C>,
    component_evals: &[Target],
    zeta_power_d: Target,
) -> Target {
    reduce_with_powers(builder, component_evals, zeta_power_d)
}

/// Evaluate `g(X, {u_i})` as defined in the Halo paper.
fn halo_g<C: Curve>(builder: &mut CircuitBuilder<C>, x: Target, us: &[Target]) -> Target {
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
