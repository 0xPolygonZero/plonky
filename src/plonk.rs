use std::borrow::Borrow;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::time::Instant;

use anyhow::Result;
use rand_chacha::ChaCha8Rng;
use rand_chacha::rand_core::SeedableRng;

use crate::{AffinePoint, blake_hash_usize_to_curve, Curve, fft_precompute, fft_with_precomputation_power_of_2, FftPrecomputation, Field, generate_rescue_constants, ifft_with_precomputation_power_of_2, msm_execute, msm_parallel, MsmPrecomputation, OpeningSet, ProjectivePoint, Proof, rescue_hash_n_to_1, rescue_hash_n_to_2, rescue_hash_n_to_3, evaluate_all_constraints, divide_by_z_h, HaloCurve, msm_precompute};
use crate::plonk_challenger::Challenger;
use crate::plonk_gates::{ArithmeticGate, Base4SumGate, BufferGate, CurveAddGate, CurveDblGate, CurveEndoGate, Gate, PublicInputGate, RescueStepAGate, RescueStepBGate};
use crate::plonk_util::{eval_l_1, eval_zero_poly, halo_n, powers, reduce_with_powers};
use crate::util::{ceil_div_usize, log2_strict, transpose};

pub(crate) const NUM_WIRES: usize = 9;
pub(crate) const NUM_ROUTED_WIRES: usize = 6;
pub(crate) const NUM_ADVICE_WIRES: usize = NUM_WIRES - NUM_ROUTED_WIRES;
pub(crate) const NUM_CONSTANTS: usize = 5;
pub(crate) const GRID_WIDTH: usize = 65;
// This is currently dominated by Base4SumGate. It has degree-4n constraints, and its prefix is 4
// bits long, so its filtered constraints are degree-8n. Dividing by Z_H makes t degree-7n.
pub(crate) const QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER: usize = 7;

pub struct PartialWitness<F: Field> {
    wire_values: HashMap<Target, F>,
}

impl<F: Field> PartialWitness<F> {
    pub fn new() -> Self {
        PartialWitness { wire_values: HashMap::new() }
    }

    pub fn is_empty(&self) -> bool {
        self.wire_values.is_empty()
    }

    pub fn contains_target(&self, target: Target) -> bool {
        self.wire_values.contains_key(&target)
    }

    pub fn contains_wire(&self, wire: Wire) -> bool {
        self.contains_target(Target::Wire(wire))
    }

    pub fn contains_all_targets(&self, targets: &[Target]) -> bool {
        targets.iter().all(|&t| self.contains_target(t))
    }

    pub fn all_populated_targets(&self) -> Vec<Target> {
        self.wire_values.keys().cloned().collect()
    }

    pub fn get_target(&self, target: Target) -> F {
        self.wire_values[&target]
    }

    pub fn get_wire(&self, wire: Wire) -> F {
        self.get_target(Target::Wire(wire))
    }

    pub fn set_target(&mut self, target: Target, value: F) {
        let opt_old_value = self.wire_values.insert(target, value);
        if let Some(old_value) = opt_old_value {
            debug_assert_eq!(old_value, value, "Target was set twice with different values");
        }
    }

    pub fn set_point_target<InnerC: Curve<BaseField=F>>(
        &mut self,
        point_target: AffinePointTarget,
        point: AffinePoint<InnerC>,
    ) {
        self.set_target(point_target.x, point.x);
        self.set_target(point_target.y, point.y);
    }

    pub fn set_wire(&mut self, wire: Wire, value: F) {
        self.set_target(Target::Wire(wire), value);
    }

    pub fn extend(&mut self, other: PartialWitness<F>) {
        for (target, value) in other.wire_values {
            self.set_target(target, value);
        }
    }
}

pub struct Witness<F: Field> {
    wire_values: Vec<Vec<F>>,
}

impl<F: Field> Witness<F> {
    pub fn get(&self, wire: Wire) -> F {
        self.wire_values[wire.gate][wire.input]
    }
}

pub trait WitnessGenerator<F: Field>: 'static {
    fn dependencies(&self) -> Vec<Target>;

    /// Given a partial witness, return any newly generated values. The caller will merge them in.
    fn generate(&self, constants: &Vec<Vec<F>>, witness: &PartialWitness<F>) -> PartialWitness<F>;
}

/// Contains all data needed to generate and/or verify proofs.
pub struct Circuit<C: HaloCurve> {
    pub security_bits: usize,
    pub num_public_inputs: usize,
    pub gate_constants: Vec<Vec<C::ScalarField>>,
    pub routing_target_partitions: TargetPartitions,
    pub generators: Vec<Box<dyn WitnessGenerator<C::ScalarField>>>,
    /// A generator of `subgroup_n`.
    pub subgroup_generator_n: C::ScalarField,
    /// A generator of `subgroup_8n`.
    pub subgroup_generator_8n: C::ScalarField,
    /// A multiplicative subgroup of order n, where n is our number of gates.
    pub subgroup_n: Vec<C::ScalarField>,
    /// A multiplicative subgroup of order 8n, where n is our number of gates.
    pub subgroup_8n: Vec<C::ScalarField>,
    /// The generators used for binding elements in a Pedersen commitment.
    pub pedersen_g: Vec<AffinePoint<C>>,
    /// The generator used for blinding Pedersen commitments.
    pub pedersen_h: AffinePoint<C>,
    /// The generator U used in Halo.
    pub u: AffinePoint<C>,
    /// Each constant polynomial, in coefficient form.
    pub constants_coeffs: Vec<Vec<C::ScalarField>>,
    /// Each constant polynomial, in point-value form, low-degree extended to be degree 8n.
    pub constants_8n: Vec<Vec<C::ScalarField>>,
    /// A commitment to each constant polynomial.
    pub c_constants: Vec<AffinePoint<C>>,
    /// Each permutation polynomial, in coefficient form.
    pub s_sigma_coeffs: Vec<Vec<C::ScalarField>>,
    /// Each permutation polynomial, low-degree extended to be degree 8n.
    pub s_sigma_values_8n: Vec<Vec<C::ScalarField>>,
    /// A commitment to each permutation polynomial.
    pub c_s_sigmas: Vec<AffinePoint<C>>,
    /// A precomputation used for MSMs involving `generators`.
    pub pedersen_g_msm_precomputation: MsmPrecomputation<C>,
    /// A precomputation used for FFTs of degree n, where n is the number of gates.
    pub fft_precomputation_n: FftPrecomputation<C::ScalarField>,
    /// A precomputation used for FFTs of degree 8n, where n is the number of gates.
    pub fft_precomputation_8n: FftPrecomputation<C::ScalarField>,
}

impl<C: HaloCurve> Circuit<C> {
    pub fn degree(&self) -> usize {
        self.gate_constants.len()
    }

    pub fn degree_pow(&self) -> usize {
        log2_strict(self.degree())
    }

    // TODO: For now we assume that there's exactly one embedded curve, InnerC.
    // Ideally it should be possible to use any number of embedded curves (including zero),
    // and we should add a set of curve gates for each embedded curve.
    pub fn generate_proof<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &self,
        witness: Witness<C::ScalarField>,
    ) -> Result<Proof<C>> {
        let mut challenger = Challenger::new(self.security_bits);

        // Convert the witness both to coefficient form and a degree-8n LDE.
        let wires_coeffs: Vec<Vec<C::ScalarField>> = witness.wire_values.iter()
            .map(|values| ifft_with_precomputation_power_of_2(values, &self.fft_precomputation_n))
            .collect();
        let wires_coeffs_8n: Vec<_> = wires_coeffs.iter()
            .map(|coeffs| self.pad_to_8n(coeffs))
            .collect();
        let wire_values_8n: Vec<_> = wires_coeffs_8n.iter()
            .map(|coeffs| fft_with_precomputation_power_of_2(coeffs, &self.fft_precomputation_8n))
            .collect();

        // Commit to the wire polynomials.
        let c_wires: Vec<ProjectivePoint<C>> = wires_coeffs.iter()
            .map(|coeffs| pedersen_hash(coeffs, &self.pedersen_g_msm_precomputation))
            .collect();
        let c_wires = ProjectivePoint::batch_to_affine(&c_wires);

        // Generate a random beta and gamma from the transcript.
        challenger.observe_affine_points(&c_wires);
        let (beta_bf, gamma_bf) = challenger.get_2_challenges();
        let beta_sf = beta_bf.try_convert::<C::ScalarField>()?;
        let gamma_sf = gamma_bf.try_convert::<C::ScalarField>()?;

        // Generate Z, which is used in Plonk's permutation argument.
        let mut plonk_z_points_n = vec![C::ScalarField::ONE];
        for i in 1..self.degree() {
            let x = self.subgroup_n[i];
            let mut numerator = C::ScalarField::ONE;
            let mut denominator = C::ScalarField::ONE;
            for j in 0..NUM_ROUTED_WIRES {
                let wire_value = witness.wire_values[j][i - 1];
                let k_i = get_subgroup_shift::<C::ScalarField>(j);
                let s_id = k_i * x;
                let s_sigma = self.s_sigma_values_8n[j][8 * i];
                numerator = numerator * (wire_value + beta_sf * s_id + gamma_sf);
                denominator = denominator * (wire_value + beta_sf * s_sigma + gamma_sf);
            }
            let last = *plonk_z_points_n.last().unwrap();
            plonk_z_points_n.push(last * numerator / denominator);
        }

        // Commit to Z.
        let plonk_z_coeffs = ifft_with_precomputation_power_of_2(&plonk_z_points_n, &self.fft_precomputation_n);
        let c_plonk_z = pedersen_hash::<C>(&plonk_z_coeffs, &self.pedersen_g_msm_precomputation).to_affine();

        // Generate a random alpha from the transcript.
        challenger.observe_affine_point(c_plonk_z);
        let alpha_bf = challenger.get_challenge();
        let alpha_sf = alpha_bf.try_convert::<C::ScalarField>()?;

        // Generate the vanishing polynomial.
        let vanishing_coeffs = self.vanishing_poly_coeffs::<InnerC>(
            &wire_values_8n, alpha_sf, beta_sf, gamma_sf, &plonk_z_coeffs);

        // Compute the quotient polynomial, t(x) = vanishing(x) / Z_H(x).
        let plonk_t_coeffs: Vec<C::ScalarField> = divide_by_z_h(&vanishing_coeffs, self.degree());

        // Split t into degree-n chunks.
        let mut plonk_t_coeff_chunks: Vec<Vec<C::ScalarField>> = plonk_t_coeffs.chunks(self.degree())
            .map(|chunk| chunk.to_vec()).collect();

        // Commit to the quotient polynomial.
        let c_plonk_t: Vec<ProjectivePoint<C>> = plonk_t_coeff_chunks.iter()
            .map(|coeffs| pedersen_hash::<C>(&plonk_t_coeffs, &self.pedersen_g_msm_precomputation))
            .collect();
        let c_plonk_t = ProjectivePoint::batch_to_affine(&c_plonk_t);

        // Generate a random zeta from the transcript.
        challenger.observe_affine_points(&c_plonk_t);
        let zeta_bf = challenger.get_challenge();
        let zeta_sf = C::try_convert_b2s(zeta_bf).expect("should fit in both fields with high probability");

        // Open all polynomials at each PublicInputGate index.
        let num_public_input_gates = ceil_div_usize(self.num_public_inputs, NUM_WIRES);
        let o_public_inputs: Vec<OpeningSet<C::ScalarField>> = (0..num_public_input_gates)
            // We place PublicInputGates at indices 0, 2, 4, ...
            .map(|i| i * 2)
            .map(|i| self.open_all_polynomials(
                &wires_coeffs, &plonk_z_coeffs, &plonk_t_coeff_chunks,
                C::ScalarField::from_canonical_usize(i)))
            .collect();

        // Open all polynomials at zeta, zeta * g, and zeta * g^65.
        let o_local = self.open_all_polynomials(&wires_coeffs, &plonk_z_coeffs, &plonk_t_coeff_chunks,
                                                zeta_sf);
        let o_right = self.open_all_polynomials(&wires_coeffs, &plonk_z_coeffs, &plonk_t_coeff_chunks,
                                                zeta_sf * self.subgroup_generator_n);
        let o_below = self.open_all_polynomials(&wires_coeffs, &plonk_z_coeffs, &plonk_t_coeff_chunks,
                                                zeta_sf * self.subgroup_generator_n.exp_usize(GRID_WIDTH));

        // Get a list of all opened values, to append to the transcript.
        let all_opening_sets: Vec<OpeningSet<C::ScalarField>> = [
            o_public_inputs.clone(),
            vec![o_local.clone(), o_right.clone(), o_below.clone()],
        ].concat();
        let all_opened_values_sf: Vec<C::ScalarField> = all_opening_sets.iter()
            .map(|os| os.to_vec()).collect::<Vec<_>>()
            .concat();
        let all_opened_values_bf: Vec<_> = all_opened_values_sf.into_iter()
            .map(|f| C::try_convert_s2b(f)
                .expect("For now, we assume that all opened values fit in both fields"))
            .collect();

        // Generate random v, u, and x from the transcript.
        challenger.observe_elements(&all_opened_values_bf);
        let (v_bf, u_bf, x) = challenger.get_3_challenges();
        let v_sf = v_bf.try_convert::<C::ScalarField>()?;
        let u_sf = u_bf.try_convert::<C::ScalarField>()?;

        // Make a list of all polynomials' commitments and coefficients, to be reduced later.
        // This must match the order of OpeningSet::to_vec.
        let all_commits = [
            self.c_constants.clone(),
            self.c_s_sigmas.clone(),
            c_wires.clone(),
            vec![c_plonk_z],
            c_plonk_t.clone(),
        ].concat();
        let all_coeffs = [
            self.constants_coeffs.clone(),
            self.s_sigma_coeffs.clone(),
            wires_coeffs,
            vec![plonk_z_coeffs],
            plonk_t_coeff_chunks,
        ].concat();

        // Normally we would reduce these lists using powers of u, but for the sake of efficiency
        // (particularly in the recursive verifier) we instead use n(u^i) for each u^i, where n is
        // the injective function related to the Halo endomorphism. Here we compute n(u^i).
        let actual_scalars: Vec<C::ScalarField> = powers(u_sf, all_commits.len()).iter()
            .map(|u_power| halo_n::<C>(&u_power.to_canonical_bool_vec()[..self.security_bits]))
            .collect();

        // Reduce the commitment list to a single commitment.
        let reduced_commit = msm_parallel(
            &actual_scalars,
            &all_commits.iter().map(AffinePoint::to_projective).collect::<Vec<_>>(),
            2);

        // Reduce the coefficient list to a single set of polynomial coefficients.
        let mut reduced_coeffs = vec![C::ScalarField::ZERO; self.degree()];
        for (i, coeffs) in all_coeffs.iter().enumerate() {
            for (j, &c) in coeffs.iter().enumerate() {
                reduced_coeffs[j] = reduced_coeffs[j] + actual_scalars[i] * c;
            }
        }

        // Reduce each opening set to a single point.
        let opening_set_reductions: Vec<C::ScalarField> = all_opening_sets.iter()
            .map(|opening_set| C::ScalarField::inner_product(&actual_scalars, &opening_set.to_vec()))
            .collect();

        // Then, we reduce the above opening set reductions to a single value.
        let reduced_opening = reduce_with_powers(&opening_set_reductions, v_sf);

        // Final IPA proof.
        let mut halo_a = reduced_coeffs;
        let mut halo_b = powers(zeta_sf, self.degree());
        let mut halo_g = AffinePoint::batch_to_projective(&self.pedersen_g);
        let mut halo_l = Vec::new();
        let mut halo_r = Vec::new();
        for j in (1..self.degree_pow()).rev() {
            let n = 1 << j;
            let middle = n / 2;

            assert_eq!(halo_a.len(), n);
            assert_eq!(halo_b.len(), n);
            assert_eq!(halo_g.len(), n);

            let a_lo = &halo_a[..middle];
            let a_hi = &halo_a[middle..];
            let b_lo = &halo_b[..middle];
            let b_hi = &halo_b[middle..];
            let g_lo = &halo_g[..middle];
            let g_hi = &halo_g[middle..];

            let l_j_blinding_factor = C::ScalarField::rand();
            let r_j_blinding_factor = C::ScalarField::rand();

            let window_size = 8;

            // L_i = <a_lo, G_hi> + [l_j] H + [<a_lo, b_hi>] U.
            let halo_l_j = msm_parallel(a_lo, g_hi, window_size)
                + C::convert(l_j_blinding_factor) * self.pedersen_h.to_projective()
                + C::convert(C::ScalarField::inner_product(a_lo, b_hi)) * self.u.to_projective();
            halo_l.push(halo_l_j);
            // R_i = <a_hi, G_lo> + [r_j] H + [<a_hi, b_lo>] U.
            let halo_r_j = msm_parallel(a_hi, g_lo, window_size)
                + C::convert(r_j_blinding_factor) * self.pedersen_h.to_projective()
                + C::convert(C::ScalarField::inner_product(a_hi, b_lo)) * self.u.to_projective();
            halo_r.push(halo_r_j);

            challenger.observe_proj_points(&[halo_l_j, halo_r_j]);
            let l_challenge_bf = challenger.get_challenge();
            let r_challenge_bf = l_challenge_bf.multiplicative_inverse().expect("This is improbable!");
            let l_challenge_sf = l_challenge_bf.try_convert::<C::ScalarField>()?;
            let r_challenge_sf = r_challenge_bf.try_convert::<C::ScalarField>()?;

            halo_a = C::ScalarField::add_slices(
                &l_challenge_sf.scale_slice(a_lo),
                &r_challenge_sf.scale_slice(a_hi));
            halo_b = C::ScalarField::add_slices(
                &l_challenge_sf.scale_slice(b_lo),
                &r_challenge_sf.scale_slice(b_hi));
            halo_g = ProjectivePoint::<C>::add_slices(
                &l_challenge_sf.scale_proj_point_slice(g_lo),
                &r_challenge_sf.scale_proj_point_slice(g_hi));
        }

        debug_assert_eq!(halo_g.len(), 1);
        let halo_g = halo_g[0].to_affine();

        Ok(Proof {
            c_wires,
            c_plonk_z,
            c_plonk_t,
            o_public_inputs,
            o_local,
            o_right,
            o_below,
            halo_l: ProjectivePoint::batch_to_affine(&halo_l),
            halo_r: ProjectivePoint::batch_to_affine(&halo_r),
            halo_g,
        })
    }

    fn vanishing_poly_coeffs<InnerC: HaloCurve<BaseField=C::ScalarField>>(
        &self,
        wire_values_8n: &[Vec<C::ScalarField>],
        alpha_sf: C::ScalarField,
        beta_sf: C::ScalarField,
        gamma_sf: C::ScalarField,
        plonk_z_coeffs: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        // Low degree extend Z.
        let plonk_z_points_8n = fft_with_precomputation_power_of_2(
            &self.pad_to_8n(&plonk_z_coeffs),
            &self.fft_precomputation_8n);

        // We will evaluate the vanishing polynomial at 8n points, then interpolate.
        let mut vanishing_points: Vec<C::ScalarField> = Vec::new();
        for (i, &x) in self.subgroup_8n.iter().enumerate() {
            // Load the constant polynomials' values at x.
            let mut local_constant_values = Vec::new();
            for j in 0..NUM_CONSTANTS {
                local_constant_values.push(self.constants_8n[j][i]);
            }

            // Load the wire polynomials' values at x, g x (the "right" position), and g^WIDTH x
            // (the "below" position). Note that a shift of 1 in the degree-n subgroup corresponds
            // to a shift of 8 in the degree-8n subgroup.
            let i_right = (i + 8) % (8 * self.degree());
            let i_below = (i + 8 * GRID_WIDTH) % (8 * self.degree());
            let mut local_wire_values = Vec::new();
            let mut right_wire_values = Vec::new();
            let mut below_wire_values = Vec::new();
            for j in 0..NUM_WIRES {
                local_wire_values.push(wire_values_8n[j][i]);
                right_wire_values.push(wire_values_8n[j][i_right]);
                below_wire_values.push(wire_values_8n[j][i_below]);
            }

            let constraint_terms = evaluate_all_constraints::<C, InnerC>(
                &local_constant_values, &local_wire_values, &right_wire_values, &below_wire_values);

            // Evaluate the L_1(x) (Z(x) - 1) vanishing term.
            let z_x = plonk_z_points_8n[i];
            let z_gz = plonk_z_points_8n[i_right];
            let vanishing_z_1_term = eval_l_1(self.degree(), x) * (z_x - C::ScalarField::ONE);

            // Evaluate the Z(x) f'(x) - g'(x) Z(g x) term.
            let mut f_prime = C::ScalarField::ONE;
            let mut g_prime = C::ScalarField::ONE;
            for j in 0..NUM_ROUTED_WIRES {
                let wire_value = wire_values_8n[i][j];
                let k_i = get_subgroup_shift::<C::ScalarField>(j);
                let s_id = k_i * x;
                let s_sigma = self.s_sigma_values_8n[j][i];
                f_prime = f_prime * (wire_value + beta_sf * s_id + gamma_sf);
                g_prime = g_prime * (wire_value + beta_sf * s_sigma + gamma_sf);
            }
            let vanishing_v_shift_term = f_prime * z_x - g_prime * z_gz;

            let vanishing_terms = [
                vec![vanishing_z_1_term],
                vec![vanishing_v_shift_term],
                constraint_terms
            ].concat();

            vanishing_points.push(reduce_with_powers(&vanishing_terms, alpha_sf));
        }

        ifft_with_precomputation_power_of_2(&vanishing_points, &self.fft_precomputation_8n)
    }

    /// Zero-pad a list of polynomial coefficients to a length of 8n, which is the degree at which
    /// we do most polynomial arithmetic.
    fn pad_to_8n(&self, coeffs: &[C::ScalarField]) -> Vec<C::ScalarField> {
        let eight_n = 8 * self.degree();
        let mut result = coeffs.to_vec();
        while result.len() < eight_n {
            result.push(C::ScalarField::ZERO);
        }
        result
    }

    /// Open each polynomial at the given point, `zeta`.
    fn open_all_polynomials(
        &self,
        wire_coeffs: &Vec<Vec<C::ScalarField>>,
        plonk_z_coeffs: &Vec<C::ScalarField>,
        plonk_t_coeffs: &Vec<Vec<C::ScalarField>>,
        zeta: C::ScalarField,
    ) -> OpeningSet<C::ScalarField> {
        let mut powers_of_zeta = vec![C::ScalarField::ONE];
        for _i in 1..self.degree() {
            powers_of_zeta.push(*powers_of_zeta.last().unwrap() * zeta);
        }

        OpeningSet {
            o_constants: self.constants_coeffs.iter().map(|c| C::ScalarField::inner_product(c, &powers_of_zeta)).collect(),
            o_plonk_sigmas: self.s_sigma_coeffs.iter().map(|c| C::ScalarField::inner_product(c, &powers_of_zeta)).collect(),
            o_wires: wire_coeffs.iter().map(|c| C::ScalarField::inner_product(c, &powers_of_zeta)).collect(),
            o_plonk_z: C::ScalarField::inner_product(&plonk_z_coeffs, &powers_of_zeta),
            o_plonk_t: plonk_t_coeffs.iter().map(|c| C::ScalarField::inner_product(c, &powers_of_zeta)).collect(),
        }
    }

    pub fn generate_witness(&self, inputs: PartialWitness<C::ScalarField>) -> Witness<C::ScalarField> {
        let start = Instant::now();

        // Index generator indices by their dependencies.
        let mut generator_indices_by_deps: HashMap<Target, Vec<usize>> = HashMap::new();
        for (i, generator) in self.generators.iter().enumerate() {
            for dep in generator.dependencies() {
                generator_indices_by_deps.entry(dep).or_insert_with(|| Vec::new()).push(i);
            }
        }

        // We start with the inputs as our witness, and execute any copy constraints.
        let mut witness = inputs;
        witness.extend(self.generate_copies(&witness, &witness.all_populated_targets()));

        // Build a list of "pending" generators which are ready to run.
        let mut pending_generator_indices = HashSet::new();
        for (i, generator) in self.generators.iter().enumerate() {
            let generator: &dyn WitnessGenerator<C::ScalarField> = generator.borrow();
            if witness.contains_all_targets(&generator.dependencies()) {
                pending_generator_indices.insert(i);
            }
        }

        // We will also keep track of which generators have already run.
        let mut completed_generator_indices = HashSet::new();

        // Now we repeat the following:
        // - Run all pending generators, keeping track of any targets that were just populated.
        // - For any newly-set targets, execute any relevant copy constraints, again tracking any
        //   newly-populated targets.
        // - Generate a new set of pending generators based on the newly-populated targets.
        while !pending_generator_indices.is_empty() {
            let mut populated_targets: Vec<Target> = Vec::new();

            for &generator_idx in &pending_generator_indices {
                let generator: &dyn WitnessGenerator<C::ScalarField> = self.generators[generator_idx].borrow();
                let result = generator.generate(&self.gate_constants, &witness);
                populated_targets.extend(result.all_populated_targets());
                witness.extend(result);
                completed_generator_indices.insert(generator_idx);
            }

            let copy_result = self.generate_copies(&witness, &populated_targets);
            populated_targets.extend(copy_result.all_populated_targets());
            witness.extend(copy_result);

            // Refresh the set of pending generators.
            pending_generator_indices.clear();
            for target in populated_targets {
                let no_indices = Vec::new();
                let affected_generator_indices = generator_indices_by_deps.get(&target).unwrap_or(&no_indices);

                for &generator_idx in affected_generator_indices {
                    // If this generator is not already pending or completed, and its dependencies
                    // are all satisfied, then add it as a pending generator.
                    let generator: &dyn WitnessGenerator<C::ScalarField> = self.generators[generator_idx].borrow();
                    if !pending_generator_indices.contains(&generator_idx)
                        && !completed_generator_indices.contains(&generator_idx)
                        && witness.contains_all_targets(&generator.dependencies()) {
                        pending_generator_indices.insert(generator_idx);
                    }
                }
            }
        }

        debug_assert_eq!(completed_generator_indices.len(), self.generators.len(),
                         "Only {} of {} generators could be run",
                         completed_generator_indices.len(), self.generators.len());

        let mut wire_values: Vec<Vec<C::ScalarField>> = Vec::new();
        for i in 0..self.degree() {
            let mut gate_i_wires = Vec::new();
            for j in 0..NUM_WIRES {
                let wire = Wire { gate: i, input: j };
                let value = if witness.contains_wire(wire) {
                    witness.get_wire(wire)
                } else {
                    // In our circuit model, a lot of wires are unused. We just set them to zero.
                    C::ScalarField::ZERO
                };
                gate_i_wires.push(value);
            }
            wire_values.push(gate_i_wires);
        }

        println!("Witness generation took {}s", start.elapsed().as_secs_f32());
        Witness { wire_values }
    }

    /// For the given set of targets, find any copy constraints involving those targets and populate
    /// the witness with copies as needed.
    fn generate_copies(&self, witness: &PartialWitness<C::ScalarField>, targets: &[Target]) -> PartialWitness<C::ScalarField> {
        let mut result = PartialWitness::new();

        for &target in targets {
            let value = witness.get_target(target);
            let partition = self.routing_target_partitions.get_partition(target);

            for &sibling in partition {
                if witness.contains_target(sibling) {
                    // This sibling's value was already set; make sure it has the same value.
                    debug_assert_eq!(witness.get_target(sibling), value);
                } else {
                    result.set_target(sibling, value);
                }
            }
        }

        result
    }
}

/// Like `pedersen_commit`, but with no blinding factor.
fn pedersen_hash<C: Curve>(
    xs: &[C::ScalarField],
    pedersen_g_msm_precomputation: &MsmPrecomputation<C>,
) -> ProjectivePoint<C> {
    msm_execute(pedersen_g_msm_precomputation, xs)
}

fn pedersen_commit<C: Curve>(
    xs: &[C::ScalarField],
    opening: C::ScalarField,
    h: AffinePoint<C>,
    pedersen_g_msm_precomputation: &MsmPrecomputation<C>,
) -> ProjectivePoint<C> {
    // TODO: Couldn't get this working with *.
    let h = h.to_projective();
    let mul_precomputation = h.mul_precompute();
    let blinding_term = h.mul_with_precomputation(opening, mul_precomputation);

    msm_execute(pedersen_g_msm_precomputation, xs) + blinding_term
}

/// A sort of proxy wire, in the context of routing and witness generation. It is not an actual
/// witness element (i.e. wire) itself, but it can be copy-constrained to wires, listed as a
/// dependency in generators, etc.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct VirtualTarget {
    pub index: usize,
}

/// Represents a wire in the circuit.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct Wire {
    /// The index of the associated gate.
    pub gate: usize,
    /// The index of the gate input wherein this wire is inserted.
    pub input: usize,
}

/// A routing target.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum Target {
    VirtualTarget(VirtualTarget),
    Wire(Wire),
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct PublicInput {
    pub index: usize,
}

/// See `PublicInputGate` for an explanation of how we make public inputs routable.
impl PublicInput {
    fn original_wire(&self) -> Wire {
        let gate = self.index / NUM_WIRES * 2;
        let input = self.index % NUM_WIRES;
        Wire { gate, input }
    }

    pub fn routable_target(&self) -> Target {
        let Wire { mut gate, mut input } = self.original_wire();
        if input > NUM_ROUTED_WIRES {
            gate += 1;
            input -= NUM_ROUTED_WIRES;
        }
        Target::Wire(Wire { gate, input })
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct AffinePointTarget {
    pub x: Target,
    pub y: Target,
}

impl AffinePointTarget {
    pub fn to_vec(&self) -> Vec<Target> {
        vec![self.x, self.y]
    }
}

/// Represents a scalar * point multiplication operation.
pub struct CurveMulOp {
    pub scalar: Target,
    pub point: AffinePointTarget,
}

pub struct CurveMulEndoResult {
    pub mul_result: AffinePointTarget,
    pub actual_scalar: Target,
}

pub struct CurveMsmEndoResult {
    pub msm_result: AffinePointTarget,
    /// While `msm` computes a sum of `[s] P` terms, `msm_endo` computes a sum of `[n(s)] P` terms
    /// for some injective `n`. Here we return each `n(s)`, i.e., the scalar by which the point was
    /// actually multiplied.
    pub actual_scalars: Vec<Target>,
}

pub struct CircuitBuilder<C: HaloCurve> {
    pub(crate) security_bits: usize,
    public_input_index: usize,
    virtual_target_index: usize,
    gate_counts: BTreeMap<&'static str, usize>,
    gate_constants: Vec<Vec<C::ScalarField>>,
    copy_constraints: Vec<(Target, Target)>,
    generators: Vec<Box<dyn WitnessGenerator<C::ScalarField>>>,
    constant_wires: HashMap<C::ScalarField, Target>,
}

impl<C: HaloCurve> CircuitBuilder<C> {
    pub fn new(security_bits: usize) -> Self {
        CircuitBuilder {
            security_bits,
            public_input_index: 0,
            virtual_target_index: 0,
            gate_counts: BTreeMap::new(),
            gate_constants: Vec::new(),
            copy_constraints: Vec::new(),
            generators: Vec::new(),
            constant_wires: HashMap::new(),
        }
    }

    pub fn stage_public_input(&mut self) -> PublicInput {
        let index = self.public_input_index;
        self.public_input_index += 1;
        PublicInput { index }
    }

    pub fn stage_public_inputs(&mut self, n: usize) -> Vec<PublicInput> {
        (0..n).map(|i| self.stage_public_input()).collect()
    }

    /// Add `PublicInputGate`s which enable public inputs to be routed. Should be called after all
    /// `stage_public_input[s]` calls, but before any gates are added.
    pub fn route_public_inputs(&mut self) {
        debug_assert_eq!(self.num_gates(), 0, "Must be called before any gates are added");
        let num_pi_gates = ceil_div_usize(self.public_input_index, NUM_WIRES);
        for i in 0..num_pi_gates {
            self.add_gate_no_constants(PublicInputGate::new(i * 2));
            self.add_gate_no_constants(BufferGate::new(i * 2 + 1));
        }
    }

    pub fn add_virtual_target(&mut self) -> Target {
        let index = self.virtual_target_index;
        self.virtual_target_index += 1;
        Target::VirtualTarget(VirtualTarget { index })
    }

    pub fn add_virtual_targets(&mut self, n: usize) -> Vec<Target> {
        (0..n).map(|i| self.add_virtual_target()).collect()
    }

    pub fn add_virtual_point_target(&mut self) -> AffinePointTarget {
        let x = self.add_virtual_target();
        let y = self.add_virtual_target();
        AffinePointTarget { x, y }
    }

    pub fn add_virtual_point_targets(&mut self, n: usize) -> Vec<AffinePointTarget> {
        (0..n).map(|i| self.add_virtual_point_target()).collect()
    }

    pub fn zero_wire(&mut self) -> Target {
        self.constant_wire(C::ScalarField::ZERO)
    }

    pub fn one_wire(&mut self) -> Target {
        self.constant_wire(C::ScalarField::ONE)
    }

    pub fn two_wire(&mut self) -> Target {
        self.constant_wire(C::ScalarField::TWO)
    }

    pub fn neg_one_wire(&mut self) -> Target {
        self.constant_wire(C::ScalarField::NEG_ONE)
    }

    pub fn constant_wire(&mut self, c: C::ScalarField) -> Target {
        if self.constant_wires.contains_key(&c) {
            self.constant_wires[&c]
        } else {
            let result = self.create_constant_wire(c);
            self.constant_wires.insert(c, result);
            result
        }
    }

    pub fn constant_wire_u32(&mut self, c: u32) -> Target {
        self.constant_wire(C::ScalarField::from_canonical_u32(c))
    }

    fn create_constant_wire(&mut self, c: C::ScalarField) -> Target {
        let index = self.num_gates();
        self.add_gate(BufferGate::new(index), vec![c]);
        Target::Wire(Wire { gate: index, input: BufferGate::<C>::WIRE_BUFFER_CONST })
    }

    pub fn constant_affine_point<InnerC: Curve<BaseField=C::ScalarField>>(
        &mut self,
        point: AffinePoint<InnerC>,
    ) -> AffinePointTarget {
        assert!(!point.zero);
        AffinePointTarget {
            x: self.constant_wire(point.x),
            y: self.constant_wire(point.y),
        }
    }

    pub fn assert_zero(&mut self, x: Target) {
        let zero = self.zero_wire();
        self.copy(x, zero);
    }

    pub fn assert_binary(&mut self, x: Target) {
        let zero = self.zero_wire();
        let one = self.one_wire();

        let x_minus_1 = self.sub(x, one);
        let product = self.mul(x, x_minus_1);
        self.assert_zero(product);
    }

    pub fn add(&mut self, x: Target, y: Target) -> Target {
        let zero = self.zero_wire();
        if x == zero {
            return y;
        }
        if y == zero {
            return x;
        }

        let one = self.one_wire();
        let index = self.num_gates();
        self.add_gate(ArithmeticGate::new(index), vec![C::ScalarField::ONE, C::ScalarField::ONE, C::ScalarField::ZERO]);
        self.copy(x, Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_0 }));
        self.copy(one, Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_1 }));
        self.copy(y, Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_ADDEND }));
        Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_OUTPUT })
    }

    pub fn add_many(&mut self, terms: &[Target]) -> Target {
        let mut sum = self.zero_wire();
        for term in terms {
            sum = self.add(sum, *term);
        }
        sum
    }

    pub fn double(&mut self, x: Target) -> Target {
        self.add(x, x)
    }

    pub fn sub(&mut self, x: Target, y: Target) -> Target {
        let zero = self.zero_wire();
        if y == zero {
            return x;
        }

        let one = self.one_wire();
        let index = self.num_gates();
        self.add_gate(ArithmeticGate::new(index), vec![C::ScalarField::ONE, C::ScalarField::NEG_ONE, C::ScalarField::ZERO]);
        self.copy(x, Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_0 }));
        self.copy(one, Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_1 }));
        self.copy(y, Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_ADDEND }));
        Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_OUTPUT })
    }

    pub fn mul(&mut self, x: Target, y: Target) -> Target {
        let one = self.one_wire();
        if x == one {
            return y;
        }
        if y == one {
            return x;
        }

        let zero = self.zero_wire();
        let index = self.num_gates();
        self.add_gate(ArithmeticGate::new(index), vec![C::ScalarField::ONE, C::ScalarField::ZERO, C::ScalarField::ZERO]);
        self.copy(x, Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_0 }));
        self.copy(y, Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_1 }));
        self.copy(zero, Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_ADDEND }));
        Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_OUTPUT })
    }

    pub fn mul_many(&mut self, terms: &[Target]) -> Target {
        let mut product = self.one_wire();
        for term in terms {
            product = self.mul(product, *term);
        }
        product
    }

    pub fn square(&mut self, x: Target) -> Target {
        self.mul(x, x)
    }

    /// Note: This assumes the most significant bit of each scalar is unset. This occurs with high
    /// probability if the field size is slightly larger than a power of two and the inputs are
    /// uniformly random.
    pub fn deterministic_square_root(&mut self, x: Target) -> Target {
        // Assume x != 0. Let y, z be the square roots of x. Since y + z = |F|, and |F| is odd, the
        // parity of y and z must differ, so we can enforce determinism by checking for a certain
        // parity bit state. We chose a parity bit of 0 since this also works for the x = 0 case.

        struct SqrtGenerator {
            x: Target,
            x_sqrt: Target,
        }

        impl<F: Field> WitnessGenerator<F> for SqrtGenerator {
            fn dependencies(&self) -> Vec<Target> {
                vec![self.x]
            }

            fn generate(&self, constants: &Vec<Vec<F>>, witness: &PartialWitness<F>) -> PartialWitness<F> {
                let x_value = witness.get_target(self.x);
                let mut x_sqrt_value = x_value.square_root().expect("Not square");

                if x_sqrt_value.to_canonical_bool_vec()[0] {
                    // The parity bit is 1; we want the other square root.
                    x_sqrt_value = -x_sqrt_value;
                    debug_assert!(!x_sqrt_value.to_canonical_bool_vec()[0]);
                }

                let mut result = PartialWitness::new();
                result.set_target(self.x_sqrt, x_sqrt_value);
                result
            }
        }

        let x_sqrt = self.add_virtual_target();
        self.add_generator(SqrtGenerator { x, x_sqrt });

        // We assume each most significant bit is unset; see the note in the method doc.
        let f_bits = C::ScalarField::BITS - 1;
        assert_eq!(f_bits, 254, "We currently only handle fields of size 2^254 + epsilon");
        let (bits, dibits) = self.split_binary_and_base_4(x_sqrt, 2, 126);

        // Verify that x_sqrt * x_sqrt = x.
        let x_sqrt_squared = self.square(x_sqrt);
        self.copy(x_sqrt_squared, x);

        // Verify that the parity bit is 0, and the other bit is binary.
        self.assert_zero(bits[0]);
        self.assert_binary(bits[1]);

        // Verify the decomposition by taking a weighted sum of the limbs. Since the first bit is
        // always zero, we start with the second bit (scaled by its weight of 2) and then add the
        // base 4 limbs.
        let mut sum = self.double(bits[1]);
        for chunk in dibits.chunks(Base4SumGate::<C>::NUM_LIMBS) {
            assert_eq!(chunk.len(), Base4SumGate::<C>::NUM_LIMBS, "Should not have a partial chunk");

            let index = self.num_gates();
            self.add_gate_no_constants(Base4SumGate::new(index));
            for i in 0..chunk.len() {
                self.copy(sum, Target::Wire(Wire { gate: index, input: Base4SumGate::<C>::WIRE_ACC_OLD }));
                self.copy(chunk[i], Target::Wire(Wire { gate: index, input: Base4SumGate::<C>::WIRE_LIMB_0 + i }));
                sum = Target::Wire(Wire { gate: index, input: Base4SumGate::<C>::WIRE_ACC_NEW })
            }
        }
        self.copy(sum, x);

        x_sqrt
    }

    /// Compute `x^power`, where `power` is a constant.
    pub fn exp_constant(&mut self, x: Target, power: C::ScalarField) -> Target {
        let power_bits = power.num_bits();
        let mut current = x;
        let mut product = self.one_wire();

        for (i, limb) in power.to_canonical_u64_vec().iter().enumerate() {
            for j in 0..64 {
                // If we've gone through all the 1 bits already, no need to keep squaring.
                let bit_index = i * 64 + j;
                if bit_index == power_bits {
                    return product;
                }

                if (limb >> j & 1) != 0 {
                    product = self.mul(product, current);
                }
                current = self.square(current);
            }
        }

        product
    }

    /// Compute `x^power`, where `power` is a constant `usize`.
    pub fn exp_constant_usize(&mut self, x: Target, power: usize) -> Target {
        self.exp_constant(x, C::ScalarField::from_canonical_usize(power))
    }

    pub fn inv(&mut self, x: Target) -> Target {
        struct InverseGenerator {
            x: Target,
            x_inv: Target,
        }

        impl<F: Field> WitnessGenerator<F> for InverseGenerator {
            fn dependencies(&self) -> Vec<Target> {
                vec![self.x]
            }

            fn generate(&self, constants: &Vec<Vec<F>>, witness: &PartialWitness<F>) -> PartialWitness<F> {
                let x_value = witness.get_target(self.x);
                let x_inv_value = x_value.multiplicative_inverse().expect("x = 0");

                let mut result = PartialWitness::new();
                result.set_target(self.x_inv, x_inv_value);
                result
            }
        }

        let x_inv = self.add_virtual_target();
        self.add_generator(InverseGenerator { x, x_inv });

        // Enforce that x * x_inv = 1.
        let product = self.mul(x, x_inv);
        let one = self.one_wire();
        self.copy(product, one);

        x_inv
    }

    pub fn div(&mut self, x: Target, y: Target) -> Target {
        let y_inv = self.inv(y);
        self.mul(x, y_inv)
    }

    /// Multiply and add; i.e. computes `x * y + z`.
    pub fn mul_add(&mut self, x: Target, y: Target, z: Target) -> Target {
        let index = self.num_gates();
        self.add_gate(ArithmeticGate::new(index), vec![C::ScalarField::ONE, C::ScalarField::ONE, C::ScalarField::ZERO]);
        self.copy(x, Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_0 }));
        self.copy(y, Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_1 }));
        self.copy(z, Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_ADDEND }));
        Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_OUTPUT })
    }

    /// Multiply and subtract; i.e. computes `x * y - z`.
    pub fn mul_sub(&mut self, x: Target, y: Target, z: Target) -> Target {
        let index = self.num_gates();
        self.add_gate(ArithmeticGate::new(index), vec![C::ScalarField::ONE, C::ScalarField::NEG_ONE, C::ScalarField::ZERO]);
        self.copy(x, Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_0 }));
        self.copy(y, Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_1 }));
        self.copy(z, Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_ADDEND }));
        Target::Wire(Wire { gate: index, input: ArithmeticGate::<C>::WIRE_OUTPUT })
    }

    /// Computes `-x`.
    pub fn neg(&mut self, x: Target) -> Target {
        let neg_one = self.neg_one_wire();
        self.mul(x, neg_one)
    }

    /// Splits `x` into its binary representation. Note that this method merely adds a generator to
    /// populate the bit wires; it does not enforce constraints to verify the decomposition.
    fn split_binary(&mut self, x: Target, num_bits: usize) -> Vec<Target> {
        let (bits, _dibits) = self.split_binary_and_base_4(x, num_bits, 0);
        bits
    }

    /// Splits `x` into a combination of binary and base 4 limbs. Note that this method merely adds
    /// a generator to populate the limb wires; it does not enforce constraints to verify the
    /// decomposition.
    fn split_binary_and_base_4(
        &mut self,
        x: Target,
        num_bits: usize,
        num_dibits: usize,
    ) -> (Vec<Target>, Vec<Target>) {
        struct SplitGenerator {
            x: Target,
            bits: Vec<Target>,
            dibits: Vec<Target>,
        }

        impl<F: Field> WitnessGenerator<F> for SplitGenerator {
            fn dependencies(&self) -> Vec<Target> {
                vec![self.x]
            }

            fn generate(&self, constants: &Vec<Vec<F>>, witness: &PartialWitness<F>) -> PartialWitness<F> {
                let x = witness.wire_values[&self.x];
                let x_bits = x.to_canonical_bool_vec();

                let mut result = PartialWitness::new();
                for i in 0..self.bits.len() {
                    result.set_target(self.bits[i], F::from_canonical_bool(x_bits[i]));
                }
                for i in 0..self.dibits.len() {
                    let bit_1 = x_bits[self.bits.len() + i * 2];
                    let bit_2 = x_bits[self.bits.len() + i * 2 + 1];
                    let dibit = if bit_1 { 1 } else { 0 } + if bit_2 { 2 } else { 0 };
                    result.set_target(self.dibits[i], F::from_canonical_u32(dibit));
                }
                result
            }
        }

        let bits = self.add_virtual_targets(num_bits);
        let dibits = self.add_virtual_targets(num_dibits);
        let generator = SplitGenerator { x, bits: bits.clone(), dibits: dibits.clone() };
        self.add_generator(generator);
        (bits, dibits)
    }

    pub fn rescue_hash_n_to_1(&mut self, inputs: &[Target]) -> Target {
        self.rescue_sponge(inputs, 1)[0]
    }

    pub fn rescue_hash_n_to_2(&mut self, inputs: &[Target]) -> (Target, Target) {
        let outputs = self.rescue_sponge(inputs, 2);
        (outputs[0], outputs[1])
    }

    pub fn rescue_hash_n_to_3(&mut self, inputs: &[Target]) -> (Target, Target, Target) {
        let outputs = self.rescue_sponge(inputs, 3);
        (outputs[0], outputs[1], outputs[2])
    }

    pub fn rescue_sponge(
        &mut self,
        inputs: &[Target],
        num_outputs: usize,
    ) -> Vec<Target> {
        // This is a r=2, c=1 sponge function with a single absorption and a single squeeze.
        let zero = self.zero_wire();
        let mut state = [zero, zero, zero];

        // Absorb all input chunks.
        for input_chunk in inputs.chunks(2) {
            for i in 0..input_chunk.len() {
                state[i] = self.add(state[i], input_chunk[i]);
            }
            state = self.rescue_permutation_3x3(state);
        }

        // Squeeze until we have the desired number of outputs.
        let mut outputs = Vec::new();
        while outputs.len() < num_outputs {
            outputs.push(state[0]);
            if outputs.len() < num_outputs {
                outputs.push(state[1]);
            }
            if outputs.len() < num_outputs {
                state = self.rescue_permutation_3x3(state);
            }
        }

        outputs
    }

    pub fn rescue_permutation_3x3(&mut self, inputs: [Target; 3]) -> [Target; 3] {
        let all_constants = generate_rescue_constants(3, self.security_bits);
        let mut state = inputs;
        for (a_constants, b_constants) in all_constants.into_iter() {
            state = self.rescue_round(state, a_constants, b_constants);
        }
        state
    }

    fn rescue_round(
        &mut self,
        inputs: [Target; 3],
        a_constants: Vec<C::ScalarField>,
        b_constants: Vec<C::ScalarField>,
    ) -> [Target; 3] {
        let a_index = self.num_gates();
        let a_gate = RescueStepAGate::new(a_index);
        self.add_gate(a_gate, a_constants);

        let b_index = self.num_gates();
        let b_gate = RescueStepBGate::new(b_index);
        self.add_gate(b_gate, b_constants);

        let a_in_0_target = Target::Wire(Wire { gate: a_index, input: RescueStepAGate::<C>::WIRE_INPUT_0 });
        let a_in_1_target = Target::Wire(Wire { gate: a_index, input: RescueStepAGate::<C>::WIRE_INPUT_1 });
        let a_in_2_target = Target::Wire(Wire { gate: a_index, input: RescueStepAGate::<C>::WIRE_INPUT_2 });
        let a_out_0_target = Target::Wire(Wire { gate: a_index, input: RescueStepAGate::<C>::WIRE_OUTPUT_0 });
        let a_out_1_target = Target::Wire(Wire { gate: a_index, input: RescueStepAGate::<C>::WIRE_OUTPUT_1 });
        let a_out_2_target = Target::Wire(Wire { gate: a_index, input: RescueStepAGate::<C>::WIRE_OUTPUT_2 });

        let b_in_0_target = Target::Wire(Wire { gate: b_index, input: RescueStepBGate::<C>::WIRE_INPUT_0 });
        let b_in_1_target = Target::Wire(Wire { gate: b_index, input: RescueStepBGate::<C>::WIRE_INPUT_1 });
        let b_in_2_target = Target::Wire(Wire { gate: b_index, input: RescueStepBGate::<C>::WIRE_INPUT_2 });
        let b_out_0_target = Target::Wire(Wire { gate: b_index, input: RescueStepBGate::<C>::WIRE_OUTPUT_0 });
        let b_out_1_target = Target::Wire(Wire { gate: b_index, input: RescueStepBGate::<C>::WIRE_OUTPUT_1 });
        let b_out_2_target = Target::Wire(Wire { gate: b_index, input: RescueStepBGate::<C>::WIRE_OUTPUT_2 });

        self.copy(inputs[0], a_in_0_target);
        self.copy(inputs[1], a_in_1_target);
        self.copy(inputs[2], a_in_2_target);
        self.copy(a_out_0_target, b_in_0_target);
        self.copy(a_out_1_target, b_in_1_target);
        self.copy(a_out_2_target, b_in_2_target);

        [b_out_0_target, b_out_1_target, b_out_2_target]
    }

    /// Assert that a given coordinate pair is on the curve `C`.
    pub fn curve_assert_valid<InnerC: Curve<BaseField=C::ScalarField>>(&mut self, p: AffinePointTarget) {
        // Recall the short Weierstrass equation: y^2 = x^3 + a*x + b.
        let a = self.constant_wire(InnerC::A);
        let b = self.constant_wire(InnerC::B);
        let y_squared = self.square(p.y);
        let x_cubed = self.exp_constant_usize(p.x, 3);
        let a_x_plus_b = self.mul_add(a, p.x, b);
        let rhs = self.add(x_cubed, a_x_plus_b);
        self.copy(y_squared, rhs);
    }

    pub fn curve_neg<InnerC: Curve<BaseField=C::ScalarField>>(
        &mut self,
        p: AffinePointTarget,
    ) -> AffinePointTarget {
        let neg_y = self.neg(p.y);
        AffinePointTarget { x: p.x, y: neg_y }
    }

    pub fn curve_add<InnerC: Curve<BaseField=C::ScalarField>>(
        &mut self,
        p_1: AffinePointTarget,
        p_2: AffinePointTarget,
    ) -> AffinePointTarget {
        let add_index = self.num_gates();
        self.add_gate_no_constants(CurveAddGate::<C, InnerC>::new(add_index));
        let buffer_index = self.num_gates();
        self.add_gate_no_constants(BufferGate::new(buffer_index));

        // TODO: Wiring.

        let result_x = Target::Wire(Wire { gate: buffer_index, input: CurveAddGate::<C, InnerC>::WIRE_GROUP_ACC_X });
        let result_y = Target::Wire(Wire { gate: buffer_index, input: CurveAddGate::<C, InnerC>::WIRE_GROUP_ACC_Y });
        AffinePointTarget { x: result_x, y: result_y }
    }


    pub fn curve_double<InnerC: Curve<BaseField=C::ScalarField>>(
        &mut self,
        p: AffinePointTarget,
    ) -> AffinePointTarget {
        let idx_dbl = self.num_gates();
        self.add_gate_no_constants(CurveDblGate::<C, InnerC>::new(idx_dbl));
        self.copy(p.x, Target::Wire(Wire { gate: idx_dbl, input: CurveDblGate::<C, InnerC>::WIRE_X_OLD }));
        self.copy(p.y, Target::Wire(Wire { gate: idx_dbl, input: CurveDblGate::<C, InnerC>::WIRE_Y_OLD }));
        AffinePointTarget {
            x: Target::Wire(Wire { gate: idx_dbl, input: CurveDblGate::<C, InnerC>::WIRE_X_NEW }),
            y: Target::Wire(Wire { gate: idx_dbl, input: CurveDblGate::<C, InnerC>::WIRE_Y_NEW }),
        }
    }

    pub fn curve_sub<InnerC: Curve<BaseField=C::ScalarField>>(
        &mut self,
        p_1: AffinePointTarget,
        p_2: AffinePointTarget,
    ) -> AffinePointTarget {
        let neg_p_2 = self.curve_neg::<InnerC>(p_2);
        self.curve_add::<InnerC>(p_1, neg_p_2)
    }

    pub fn curve_mul<InnerC: Curve<BaseField=C::ScalarField>>(&mut self, mul: CurveMulOp) -> AffinePointTarget {
        self.curve_msm::<InnerC>(&[mul])
    }

    pub fn curve_mul_endo<InnerC: HaloCurve<BaseField=C::ScalarField>>(
        &mut self,
        mul: CurveMulOp,
    ) -> CurveMulEndoResult {
        let result = self.curve_msm_endo::<InnerC>(&[mul]);
        CurveMulEndoResult {
            mul_result: result.msm_result,
            actual_scalar: result.actual_scalars[0],
        }
    }

    /// Note: This assumes the most significant bit of each scalar is unset. This occurs with high
    /// probability if the field size is slightly larger than a power of two and the scalars are
    /// uniformly random.
    pub fn curve_msm<InnerC: Curve<BaseField=C::ScalarField>>(
        &mut self,
        parts: &[CurveMulOp],
    ) -> AffinePointTarget {
        // We assume each most significant bit is unset; see the note in the method doc.
        let f_bits = C::ScalarField::BITS - 1;

        let all_bits: Vec<Vec<Target>> = parts.iter()
            .map(|part| self.split_binary(part.scalar, f_bits))
            .collect();

        // Normally we would start with zero, but to avoid exceptional cases, we start with some
        // arbitrary nonzero point and subtract it later. This avoids exception with high
        // probability provided that the scalars and points are random. (We don't worry about
        // non-random inputs from malicious provers, since our curve gates will be unsatisfiable in
        // exceptional cases.)
        let mut filler = InnerC::GENERATOR_AFFINE;
        let mut acc = self.constant_affine_point(filler);
        let mut scalar_accs = vec![self.zero_wire(); parts.len()];

        for i in (0..f_bits).rev() {
            // Route the accumulator to the first curve addition gate's inputs.
            self.copy(acc.x, Target::Wire(
                Wire { gate: self.num_gates(), input: CurveAddGate::<C, InnerC>::WIRE_GROUP_ACC_X }));
            self.copy(acc.y, Target::Wire(
                Wire { gate: self.num_gates(), input: CurveAddGate::<C, InnerC>::WIRE_GROUP_ACC_Y }));

            for (j, part) in parts.iter().enumerate() {
                let bit = all_bits[j][i];

                let idx_add = self.num_gates();
                self.add_gate_no_constants(CurveAddGate::<C, InnerC>::new(idx_add));
                self.copy(scalar_accs[j], Target::Wire(
                    Wire { gate: idx_add, input: CurveAddGate::<C, InnerC>::WIRE_SCALAR_ACC_OLD }));
                scalar_accs[j] = Target::Wire(
                    Wire { gate: idx_add, input: CurveAddGate::<C, InnerC>::WIRE_SCALAR_ACC_NEW });
                self.copy(part.point.x, Target::Wire(
                    Wire { gate: idx_add, input: CurveAddGate::<C, InnerC>::WIRE_ADDEND_X }));
                self.copy(part.point.y, Target::Wire(
                    Wire { gate: idx_add, input: CurveAddGate::<C, InnerC>::WIRE_ADDEND_Y }));
                self.copy(bit, Target::Wire(
                    Wire { gate: idx_add, input: CurveAddGate::<C, InnerC>::WIRE_SCALAR_BIT }));
            }

            // Double the accumulator.
            let idx_dbl = self.num_gates();
            self.add_gate_no_constants(CurveDblGate::<C, InnerC>::new(idx_dbl));
            // No need to route the double gate's inputs, because the last add gate would have
            // constrained them.
            // Normally, we will take the double gate's outputs as the new accumulator. If we just
            // completed the last iteration of the MSM though, then we don't want to perform a final
            // doubling, so we will take its inputs as the result instead.
            if i == f_bits - 1 {
                acc = AffinePointTarget {
                    x: Target::Wire(Wire { gate: idx_dbl, input: CurveDblGate::<C, InnerC>::WIRE_X_OLD }),
                    y: Target::Wire(Wire { gate: idx_dbl, input: CurveDblGate::<C, InnerC>::WIRE_Y_OLD }),
                };
            } else {
                acc = AffinePointTarget {
                    x: Target::Wire(Wire { gate: idx_dbl, input: CurveDblGate::<C, InnerC>::WIRE_X_NEW }),
                    y: Target::Wire(Wire { gate: idx_dbl, input: CurveDblGate::<C, InnerC>::WIRE_Y_NEW }),
                };
            }

            // Also double the filler, so we can subtract out a rescaled version later.
            filler = filler.double();
        }

        // Subtract (a rescaled version of) the arbitrary nonzero value that we started with.
        let filler_target = self.constant_affine_point(filler);
        acc = self.curve_sub::<InnerC>(acc, filler_target);

        // Assert that each accumulation of scalar bits matches the original scalar.
        for (j, part) in parts.iter().enumerate() {
            self.copy(scalar_accs[j], part.scalar);
        }

        acc
    }

    /// Like `curve_msm`, but uses the endomorphism described in the Halo paper.
    pub fn curve_msm_endo<InnerC: HaloCurve<BaseField=C::ScalarField>>(
        &mut self,
        parts: &[CurveMulOp],
    ) -> CurveMsmEndoResult {
        let zero = self.zero_wire();

        // We assume each most significant bit is unset; see the note in curve_msm's method doc.
        let f_bits = C::ScalarField::BITS - 1;
        let scalar_bits = self.security_bits;
        let scalar_dibits = (f_bits - scalar_bits) / 2;

        // To keep things simple for now, we only handle the case of |F| ~= 2^254 and lambda = 128.
        assert_eq!(f_bits, 254);
        assert_eq!(scalar_bits, 128);
        assert_eq!(scalar_dibits, 63);

        // We split each scalar into 128 bits and 63 dibits. The bits are used in the MSM, while the
        // dibits are ignored, except that we need to include them in our sum-of-limbs computation
        // in order to verify that the decomposition was correct.
        let (all_bits, all_dibits): (Vec<Vec<Target>>, Vec<Vec<Target>>) = parts.iter()
            .map(|part| self.split_binary_and_base_4(part.scalar, scalar_bits, scalar_dibits))
            .unzip();

        // Normally we would start with zero, but to avoid exceptional cases, we start with some
        // arbitrary nonzero point and subtract it later. This avoids exception with high
        // probability provided that the scalars and points are random. (We don't worry about
        // non-random inputs from malicious provers, since our curve gates will be unsatisfiable in
        // exceptional cases.)
        let mut filler = InnerC::GENERATOR_AFFINE;
        let mut acc = self.constant_affine_point(filler);

        // For each scalar, we maintain two accumulators. The unsigned one is for computing a
        // weighted sum of bits and dibits in the usual manner, so that we can later check that this
        // sum equals the original scalar. The signed one is for computing n(s) for each scalar s.
        // This is the "actual" scalar by which the associated point was multiplied, accounting for
        // the endomorphism.
        let mut scalar_acc_unsigned = Vec::new();
        let mut scalar_acc_signed = Vec::new();

        // As in the Halo paper, we process two scalar bits at a time.
        for i in (0..scalar_bits).step_by(2).rev() {
            // Route the point accumulator to the first gate's inputs.
            self.copy(acc.x, Target::Wire(Wire { gate: self.num_gates(), input: CurveEndoGate::<C, InnerC>::WIRE_GROUP_ACC_X }));
            self.copy(acc.y, Target::Wire(Wire { gate: self.num_gates(), input: CurveEndoGate::<C, InnerC>::WIRE_GROUP_ACC_Y }));

            for (j, part) in parts.iter().enumerate() {
                let bit_0 = all_bits[j][i];
                let bit_1 = all_bits[j][i + 1];

                let gate = self.num_gates();
                self.add_gate_no_constants(CurveEndoGate::<C, InnerC>::new(gate));

                self.copy(part.point.x, Target::Wire(Wire { gate, input: CurveEndoGate::<C, InnerC>::WIRE_ADDEND_X }));
                self.copy(part.point.y, Target::Wire(Wire { gate, input: CurveEndoGate::<C, InnerC>::WIRE_ADDEND_Y }));
                self.copy(bit_0, Target::Wire(Wire { gate, input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_BIT_0 }));
                self.copy(bit_1, Target::Wire(Wire { gate, input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_BIT_1 }));

                // If this is the first pair of scalar bits being processed, route 0 to the scalar accumulators.
                if i == scalar_bits - 2 {
                    self.copy(zero, Target::Wire(Wire { gate, input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_ACC_UNSIGNED }));
                    self.copy(zero, Target::Wire(Wire { gate, input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_ACC_SIGNED }));
                }

                // If this is the last pair of scalar bits being processed, save the final
                // accumulator states.
                // Since CurveEndoGate will store these in the "next" gate, but this is the last
                // CurveEndoGate for this scalar, we need to add an extra BufferGate to receive them.
                if i == 0 {
                    let gate = self.num_gates();
                    self.add_gate_no_constants(BufferGate::new(gate));
                    scalar_acc_unsigned.push(Target::Wire(Wire { gate, input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_ACC_UNSIGNED }));
                    scalar_acc_signed.push(Target::Wire(Wire { gate, input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_ACC_SIGNED }));
                }
            }

            // Double the accumulator.
            let gate = self.num_gates();
            self.add_gate_no_constants(CurveDblGate::<C, InnerC>::new(gate));
            // No need to route the double gate's inputs, because the last endo gate would have
            // constrained them.
            // Normally, we will take the double gate's outputs as the new accumulator. If we just
            // completed the last iteration of the MSM though, then we don't want to perform a final
            // doubling, so we will take its inputs as the result instead.
            if i == scalar_bits - 1 {
                acc = AffinePointTarget {
                    x: Target::Wire(Wire { gate, input: CurveDblGate::<C, InnerC>::WIRE_X_OLD }),
                    y: Target::Wire(Wire { gate, input: CurveDblGate::<C, InnerC>::WIRE_Y_OLD }),
                };
            } else {
                acc = AffinePointTarget {
                    x: Target::Wire(Wire { gate, input: CurveDblGate::<C, InnerC>::WIRE_X_NEW }),
                    y: Target::Wire(Wire { gate, input: CurveDblGate::<C, InnerC>::WIRE_Y_NEW }),
                };
            }

            // Also double the filler, so we can subtract out a rescaled version later.
            filler = filler.double();
        }

        // Subtract (a rescaled version of) the arbitrary nonzero value that we started with.
        let filler_target = self.constant_affine_point(filler);
        acc = self.curve_sub::<InnerC>(acc, filler_target);

        // By now we've accumulated all the bits of each scalar, but we also need to accumulate the dibits.
        for j in 0..parts.len() {
            for dibits_chunk in all_dibits[j].chunks(Base4SumGate::<C>::NUM_LIMBS) {
                assert_eq!(dibits_chunk.len(), Base4SumGate::<C>::NUM_LIMBS);

                let gate = self.num_gates();
                self.add_gate_no_constants(Base4SumGate::new(gate));
                self.copy(scalar_acc_unsigned[j], Target::Wire(Wire { gate, input: Base4SumGate::<C>::WIRE_ACC_OLD }));
                scalar_acc_unsigned[j] = Target::Wire(Wire { gate, input: Base4SumGate::<C>::WIRE_ACC_NEW });

                for (i, &dibit) in dibits_chunk.iter().enumerate() {
                    self.copy(dibit, Target::Wire(Wire { gate, input: Base4SumGate::<C>::WIRE_LIMB_0 + i }));
                }
            }
        }

        // Finally, assert that each unsigned accumulator matches the original scalar.
        for (j, part) in parts.iter().enumerate() {
            self.copy(scalar_acc_unsigned[j], part.scalar);
        }

        CurveMsmEndoResult { msm_result: acc, actual_scalars: scalar_acc_signed }
    }

    /// Adds a gate to the circuit, without doing any routing.
    fn add_gate_no_constants<G: Gate<C>>(&mut self, gate: G) {
        self.add_gate(gate, Vec::new());
    }

    /// Adds a gate to the circuit, without doing any routing.
    pub fn add_gate<G: Gate<C>>(&mut self, gate: G, gate_constants: Vec<C::ScalarField>) {
        debug_assert!(G::PREFIX.len() + gate_constants.len() <= NUM_CONSTANTS);

        // Merge the gate type's prefix bits with the given gate config constants.
        let mut all_constants = Vec::new();
        for &prefix_bit in G::PREFIX {
            all_constants.push(C::ScalarField::from_canonical_bool(prefix_bit));
        }
        all_constants.extend(gate_constants);

        // Pad if not all constants were used.
        while all_constants.len() < NUM_CONSTANTS {
            all_constants.push(C::ScalarField::ZERO);
        }

        self.gate_constants.push(all_constants);
        self.add_generator(gate);
        *self.gate_counts.entry(G::NAME).or_insert(0) += 1;
    }

    pub fn add_generator<G: WitnessGenerator<C::ScalarField>>(&mut self, gate: G) {
        self.generators.push(Box::new(gate));
    }

    pub fn num_gates(&self) -> usize {
        self.gate_constants.len()
    }

    /// Add a copy constraint between two routing targets.
    pub fn copy(&mut self, target_1: Target, target_2: Target) {
        self.copy_constraints.push((target_1, target_2));
    }

    /// Adds a gate with random wire values. By adding `k` of these gates, we can ensure that
    /// nothing is learned by opening the wire polynomials at `k` points outside of H.
    fn add_blinding_gate(&mut self) {
        let gate = self.num_gates();
        self.add_gate_no_constants(BufferGate::new(gate));
        for input in 0..NUM_WIRES {
            self.add_generator(RandomGenerator { target: Target::Wire(Wire { gate, input }) });
        }

        struct RandomGenerator {
            target: Target,
        }

        impl<F: Field> WitnessGenerator<F> for RandomGenerator {
            fn dependencies(&self) -> Vec<Target> {
                vec![]
            }

            fn generate(&self, constants: &Vec<Vec<F>>, witness: &PartialWitness<F>) -> PartialWitness<F> {
                let mut result = PartialWitness::new();
                result.set_target(self.target, F::rand());
                result
            }
        }
    }

    pub fn build(mut self) -> Circuit<C> {
        // Since we will open each polynomial at three points outside of H, we need three random
        // values to ensure nothing is learned from the out-of-H openings.
        for _i in 0..3 {
            self.add_blinding_gate();
        }

        // Print gate counts.
        println!("Gate counts:");
        for (gate, count) in &self.gate_counts {
            println!("{}: {}", gate, count);
        }
        println!();

        // Pad to a power of two.
        println!("Total gates before padding: {}", self.num_gates());
        while !self.num_gates().is_power_of_two() {
            // Add an empty gate.
            self.add_gate_no_constants(BufferGate::new(self.num_gates()));
        }
        println!("Total gates after padding: {}", self.num_gates());

        let degree = self.num_gates();
        let degree_pow = log2_strict(degree);
        let routing_target_partitions = self.get_routing_partitions();
        let wire_partitions = routing_target_partitions.to_wire_partitions();
        let sigma = wire_partitions.to_sigma();

        let CircuitBuilder {
            security_bits,
            public_input_index: num_public_inputs,
            gate_constants,
            generators,
            ..
        } = self;

        let fft_precomputation_n = fft_precompute(degree);
        let fft_precomputation_8n = fft_precompute(degree * 8);

        let subgroup_generator_n = C::ScalarField::primitive_root_of_unity(degree_pow);
        let subgroup_generator_8n = C::ScalarField::primitive_root_of_unity(degree_pow + 3);
        let subgroup_n = C::ScalarField::cyclic_subgroup_known_order(subgroup_generator_n, degree);
        let subgroup_8n = C::ScalarField::cyclic_subgroup_known_order(subgroup_generator_n, 8 * degree);

        let pedersen_g: Vec<_> = (0..degree).map(|i| blake_hash_usize_to_curve::<C>(i)).collect();
        let pedersen_h = blake_hash_usize_to_curve::<C>(degree);
        let u = blake_hash_usize_to_curve::<C>(degree + 1);

        let w = 8; // TODO: Should really be set dynamically based on MSM size.
        let pedersen_g_msm_precomputation = msm_precompute(&AffinePoint::batch_to_projective(&pedersen_g), w);

        // While gate_constants is indexed by gate index first, this is indexed by wire index first.
        let wire_constants = transpose::<C::ScalarField>(&gate_constants);

        let constants_coeffs: Vec<Vec<C::ScalarField>> = wire_constants.iter()
            .map(|values| ifft_with_precomputation_power_of_2(values, &fft_precomputation_n))
            .collect();
        let constants_8n = constants_coeffs.iter()
            .map(|coeffs| fft_with_precomputation_power_of_2(coeffs, &fft_precomputation_8n))
            .collect();
        let c_constants: Vec<ProjectivePoint<C>> = constants_coeffs.iter()
            .map(|coeffs| pedersen_hash(&coeffs, &pedersen_g_msm_precomputation))
            .collect();
        let c_constants = ProjectivePoint::batch_to_affine(&c_constants);

        // Convert sigma's values to scalar field elements and split it into degree-n chunks.
        let sigma_chunks: Vec<Vec<C::ScalarField>> = sigma.into_iter()
            .map(|x| C::ScalarField::from_canonical_usize(x))
            .collect::<Vec<_>>()
            .chunks(degree)
            .map(|chunk| chunk.to_vec())
            .collect();

        // Compute S_sigma, then a commitment to it.
        let s_sigma_coeffs: Vec<Vec<C::ScalarField>> = sigma_chunks.iter()
            .map(|sigma_chunk| ifft_with_precomputation_power_of_2(sigma_chunk, &fft_precomputation_n))
            .collect();
        let s_sigma_values_8n: Vec<Vec<C::ScalarField>> = s_sigma_coeffs.iter()
            .map(|coeffs| fft_with_precomputation_power_of_2(coeffs, &fft_precomputation_8n))
            .collect();
        let c_s_sigmas: Vec<ProjectivePoint<C>> = s_sigma_coeffs.iter()
            .map(|coeffs| pedersen_hash(&coeffs, &pedersen_g_msm_precomputation))
            .collect();
        let c_s_sigmas = ProjectivePoint::batch_to_affine(&c_s_sigmas);

        Circuit {
            security_bits,
            num_public_inputs,
            gate_constants,
            routing_target_partitions,
            generators,
            subgroup_generator_n,
            subgroup_generator_8n,
            subgroup_n,
            subgroup_8n,
            pedersen_g,
            pedersen_h,
            u,
            constants_coeffs,
            constants_8n,
            c_constants,
            s_sigma_coeffs,
            s_sigma_values_8n,
            c_s_sigmas,
            pedersen_g_msm_precomputation,
            fft_precomputation_n,
            fft_precomputation_8n,
        }
    }

    fn get_routing_partitions(&self) -> TargetPartitions {
        let mut partitions = TargetPartitions::new();

        for i in 0..self.virtual_target_index {
            partitions.add_partition(Target::VirtualTarget(VirtualTarget { index: i }));
        }

        for gate in 0..self.num_gates() {
            for input in 0..NUM_WIRES {
                partitions.add_partition(Target::Wire(Wire { gate, input }));
            }
        }

        for &(a, b) in &self.copy_constraints {
            partitions.merge(a, b);
        }

        partitions
    }
}

pub struct TargetPartitions {
    partitions: Vec<Vec<Target>>,
    indices: HashMap<Target, usize>,
}

impl TargetPartitions {
    fn new() -> Self {
        Self { partitions: Vec::new(), indices: HashMap::new() }
    }

    fn get_partition(&self, target: Target) -> &[Target] {
        &self.partitions[self.indices[&target]]
    }

    /// Add a new partition with a single member.
    fn add_partition(&mut self, target: Target) {
        let index = self.partitions.len();
        self.partitions.push(vec![target]);
        self.indices.insert(target, index);
    }

    /// Merge the two partitions containing the two given targets. Does nothing if the targets are
    /// already members of the same partition.
    fn merge(&mut self, a: Target, b: Target) {
        let a_index = self.indices[&a];
        let b_index = self.indices[&b];
        if a_index != b_index {
            // Merge a's partition into b's partition, leaving a's partition empty.
            // We have to clone because Rust's borrow checker doesn't know that
            // self.partitions[b_index] and self.partitions[b_index] are disjoint.
            let mut a_partition = self.partitions[a_index].clone();
            let b_partition = &mut self.partitions[b_index];
            for a_sibling in &a_partition {
                *self.indices.get_mut(a_sibling).unwrap() = b_index;
            }
            b_partition.append(&mut a_partition);
        }
    }

    fn to_wire_partitions(&self) -> WirePartitions {
        // Here we just drop all CircuitInputs, leaving all GateInputs.
        let mut partitions = Vec::new();
        let mut indices = HashMap::new();

        for old_partition in &self.partitions {
            let mut new_partition = Vec::new();
            for target in old_partition {
                if let Target::Wire(gi) = *target {
                    new_partition.push(gi);
                }
            }
            partitions.push(new_partition);
        }

        for (&target, &index) in &self.indices {
            if let Target::Wire(gi) = target {
                indices.insert(gi, index);
            }
        }

        WirePartitions { partitions, indices }
    }
}

struct WirePartitions {
    partitions: Vec<Vec<Wire>>,
    indices: HashMap<Wire, usize>,
}

impl WirePartitions {
    /// Find a wire's "neighbor" in the context of Plonk's "extended copy constraints" check. In
    /// other words, find the next wire in the given wire's partition. If the given wire is last in
    /// its partition, this will loop around. If the given wire has a partition all to itself, it
    /// is considered its own neighbor.
    fn get_neighbor(&self, wire: Wire) -> Wire {
        let partition = &self.partitions[self.indices[&wire]];
        let n = partition.len();
        for i in 0..n {
            if partition[i] == wire {
                let neighbor_index = (i + 1) % n;
                return partition[neighbor_index];
            }
        }
        panic!("Wire not found in the expected partition")
    }

    /// Generates sigma in the context of Plonk, which is a map from `[kn]` to `[kn]`, where `k` is
    /// the number of wires and `n` is the number of gates.
    fn to_sigma(&self) -> Vec<usize> {
        let kn = self.indices.len();
        let k = NUM_WIRES;
        let n = kn / k;
        debug_assert_eq!(k * n, kn);

        let mut sigma = Vec::new();
        for input in 0..k {
            for gate in 0..n {
                let wire = Wire { gate, input };
                let neighbor = self.get_neighbor(wire);
                sigma.push(neighbor.input * n + neighbor.gate);
            }
        }
        sigma
    }
}

/// Returns `k_i`, the multiplier used in `S_ID_i` in the context of Plonk's permutation argument.
pub(crate) fn get_subgroup_shift<F: Field>(i: usize) -> F {
    // The optimized variant of Plonk's permutation argument calls for NUM_ROUTED_WIRES shifts,
    // k_1, ..., k_n, which result in distinct cosets. The paper suggests a method which is
    // fairly straightforward when only three shifts are needed, but seems a bit complex and
    // expensive if more are needed.

    // We will "cheat" and just use random field elements. Since our subgroup has |F*|/degree
    // possible cosets, the probability of a collision is negligible for large fields.

    if i == 0 {
        // We use a trivial shift of 1 for k_1, as in the paper, to save a multiplication.
        return F::ONE;
    } else {
        let mut rng = ChaCha8Rng::seed_from_u64(i as u64);
        F::rand_from_rng(&mut rng)
    }
}
