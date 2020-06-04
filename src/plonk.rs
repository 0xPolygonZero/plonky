use std::borrow::Borrow;
use std::collections::{HashMap, HashSet};
use std::time::Instant;

use anyhow::Result;

use crate::partition::{get_subgroup_shift, TargetPartitions};
use crate::plonk_challenger::Challenger;
use crate::plonk_util::{
    coeffs_to_commitments, coeffs_to_values_padded, eval_coeffs, eval_l_1, halo_n,
    pad_to_8n, pedersen_hash, permutation_polynomial, powers, reduce_with_powers, values_to_coeffs,
};
use crate::target::{Target};
use crate::util::{ceil_div_usize, log2_strict};
use crate::witness::{PartialWitness, Witness, WitnessGenerator};
use crate::{
    divide_by_z_h, evaluate_all_constraints, fft_precompute,
    fft_with_precomputation_power_of_2, generate_rescue_constants,
    ifft_with_precomputation_power_of_2, msm_execute, msm_parallel, msm_precompute,
    rescue_hash_n_to_1, rescue_hash_n_to_2, rescue_hash_n_to_3, AffinePoint, CircuitBuilder,
    ConstantGate, Curve, FftPrecomputation, Field, HaloCurve, MsmPrecomputation, OpeningSet,
    ProjectivePoint, Proof,
};

pub(crate) const NUM_WIRES: usize = 9;
pub(crate) const NUM_ROUTED_WIRES: usize = 6;
pub(crate) const NUM_ADVICE_WIRES: usize = NUM_WIRES - NUM_ROUTED_WIRES;
pub(crate) const NUM_CONSTANTS: usize = 5;
pub(crate) const GRID_WIDTH: usize = 65;
// This is currently dominated by Base4SumGate. It has degree-4n constraints, and its prefix is 4
// bits long, so its filtered constraints are degree-8n. Dividing by Z_H makes t degree-7n.
pub(crate) const QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER: usize = 7;

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
        let wire_values_by_wire_index = &witness.transpose();
        let wires_coeffs = values_to_coeffs(&wire_values_by_wire_index, &self.fft_precomputation_n);
        let wire_values_8n = coeffs_to_values_padded(&wires_coeffs, &self.fft_precomputation_8n);

        // Commit to the wire polynomials.
        let c_wires = coeffs_to_commitments(&wires_coeffs, &self.pedersen_g_msm_precomputation);

        // Generate a random beta and gamma from the transcript.
        challenger.observe_affine_points(&c_wires);
        let (beta_bf, gamma_bf) = challenger.get_2_challenges();
        let beta_sf = beta_bf.try_convert::<C::ScalarField>()?;
        let gamma_sf = gamma_bf.try_convert::<C::ScalarField>()?;

        let plonk_z_points_n = permutation_polynomial(
            self.degree(),
            &self.subgroup_8n,
            &witness,
            &self.s_sigma_values_8n,
            beta_sf,
            gamma_sf,
        );
        // Commit to Z.
        let plonk_z_coeffs =
            ifft_with_precomputation_power_of_2(&plonk_z_points_n, &self.fft_precomputation_n);
        let c_plonk_z =
            pedersen_hash::<C>(&plonk_z_coeffs, &self.pedersen_g_msm_precomputation).to_affine();

        // Generate a random alpha from the transcript.
        challenger.observe_affine_point(c_plonk_z);
        let alpha_bf = challenger.get_challenge();
        let alpha_sf = alpha_bf.try_convert::<C::ScalarField>()?;

        // Generate the vanishing polynomial.
        let vanishing_coeffs = self.vanishing_poly_coeffs::<InnerC>(
            &wire_values_8n,
            alpha_sf,
            beta_sf,
            gamma_sf,
            &plonk_z_coeffs,
        );

        // Compute the quotient polynomial, t(x) = vanishing(x) / Z_H(x).
        let plonk_t_coeffs: Vec<C::ScalarField> = divide_by_z_h(&vanishing_coeffs, self.degree());

        // Split t into degree-n chunks.
        let plonk_t_coeff_chunks: Vec<Vec<C::ScalarField>> = plonk_t_coeffs
            .chunks(self.degree())
            .map(|chunk| chunk.to_vec())
            .collect();

        // Commit to the quotient polynomial.
        let c_plonk_t =
            coeffs_to_commitments(&plonk_t_coeff_chunks, &self.pedersen_g_msm_precomputation);

        // Generate a random zeta from the transcript.
        challenger.observe_affine_points(&c_plonk_t);
        let zeta_bf = challenger.get_challenge();
        let zeta_sf =
            C::try_convert_b2s(zeta_bf).expect("should fit in both fields with high probability");

        // Open all polynomials at each PublicInputGate index.
        let num_public_input_gates = ceil_div_usize(self.num_public_inputs, NUM_WIRES);
        let o_public_inputs: Vec<OpeningSet<C::ScalarField>> = (0..num_public_input_gates)
            // We place PublicInputGates at indices 0, 2, 4, ...
            .map(|i| i * 2)
            .map(|i| {
                self.open_all_polynomials(
                    &wires_coeffs,
                    &plonk_z_coeffs,
                    &plonk_t_coeff_chunks,
                    C::ScalarField::from_canonical_usize(i),
                )
            })
            .collect();

        // Open all polynomials at zeta, zeta * g, and zeta * g^65.
        let o_local = self.open_all_polynomials(
            &wires_coeffs,
            &plonk_z_coeffs,
            &plonk_t_coeff_chunks,
            zeta_sf,
        );
        let o_right = self.open_all_polynomials(
            &wires_coeffs,
            &plonk_z_coeffs,
            &plonk_t_coeff_chunks,
            zeta_sf * self.subgroup_generator_n,
        );
        let o_below = self.open_all_polynomials(
            &wires_coeffs,
            &plonk_z_coeffs,
            &plonk_t_coeff_chunks,
            zeta_sf * self.subgroup_generator_n.exp_usize(GRID_WIDTH),
        );

        // Get a list of all opened values, to append to the transcript.
        let all_opening_sets: Vec<OpeningSet<C::ScalarField>> = [
            o_public_inputs.clone(),
            vec![o_local.clone(), o_right.clone(), o_below.clone()],
        ]
        .concat();
        let all_opened_values_sf: Vec<C::ScalarField> = all_opening_sets
            .iter()
            .map(|os| os.to_vec())
            .collect::<Vec<_>>()
            .concat();
        let all_opened_values_bf: Vec<_> = all_opened_values_sf
            .into_iter()
            .map(|f| {
                C::try_convert_s2b(f)
                    .expect("For now, we assume that all opened values fit in both fields")
            })
            .collect();

        // Generate random v, u, and x from the transcript.
        challenger.observe_elements(&all_opened_values_bf);
        let (v_bf, u_bf, _x) = challenger.get_3_challenges();
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
        ]
        .concat();
        let all_coeffs = [
            self.constants_coeffs.clone(),
            self.s_sigma_coeffs.clone(),
            wires_coeffs,
            vec![plonk_z_coeffs],
            plonk_t_coeff_chunks,
        ]
        .concat();

        // Normally we would reduce these lists using powers of u, but for the sake of efficiency
        // (particularly in the recursive verifier) we instead use n(u^i) for each u^i, where n is
        // the injective function related to the Halo endomorphism. Here we compute n(u^i).
        let actual_scalars: Vec<C::ScalarField> = powers(u_sf, all_commits.len())
            .iter()
            .map(|u_power| halo_n::<C>(&u_power.to_canonical_bool_vec()[..self.security_bits]))
            .collect();

        // Reduce the commitment list to a single commitment.
        let _reduced_commit = msm_parallel(
            &actual_scalars,
            &all_commits
                .iter()
                .map(AffinePoint::to_projective)
                .collect::<Vec<_>>(),
            2,
        );

        // Reduce the coefficient list to a single set of polynomial coefficients.
        let mut reduced_coeffs = vec![C::ScalarField::ZERO; self.degree()];
        for (i, coeffs) in all_coeffs.iter().enumerate() {
            for (j, &c) in coeffs.iter().enumerate() {
                reduced_coeffs[j] = reduced_coeffs[j] + actual_scalars[i] * c;
            }
        }

        // Reduce each opening set to a single point.
        let opening_set_reductions: Vec<C::ScalarField> = all_opening_sets
            .iter()
            .map(|opening_set| {
                C::ScalarField::inner_product(&actual_scalars, &opening_set.to_vec())
            })
            .collect();

        // Then, we reduce the above opening set reductions to a single value.
        let _reduced_opening = reduce_with_powers(&opening_set_reductions, v_sf);

        // Final IPA proof.
        let mut halo_a = reduced_coeffs;
        let mut halo_b = powers(zeta_sf, self.degree());
        let mut halo_g = AffinePoint::batch_to_projective(&self.pedersen_g);
        let mut halo_l = Vec::new();
        let mut halo_r = Vec::new();
        for j in (1..=self.degree_pow()).rev() {
            let n = 1 << j;
            let middle = n / 2;

            debug_assert_eq!(halo_a.len(), n);
            debug_assert_eq!(halo_b.len(), n);
            debug_assert_eq!(halo_g.len(), n);

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
            let r_challenge_bf = l_challenge_bf
                .multiplicative_inverse()
                .expect("This is improbable!");
            let l_challenge_sf = l_challenge_bf.try_convert::<C::ScalarField>()?;
            let r_challenge_sf = r_challenge_bf.try_convert::<C::ScalarField>()?;

            halo_a = C::ScalarField::add_slices(
                &l_challenge_sf.scale_slice(a_lo),
                &r_challenge_sf.scale_slice(a_hi),
            );
            halo_b = C::ScalarField::add_slices(
                &l_challenge_sf.scale_slice(b_lo),
                &r_challenge_sf.scale_slice(b_hi),
            );
            halo_g = ProjectivePoint::<C>::add_slices(
                &l_challenge_sf.scale_proj_point_slice(g_lo),
                &r_challenge_sf.scale_proj_point_slice(g_hi),
            );
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

    fn vanishing_poly_coeffs<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &self,
        wire_values_8n: &[Vec<C::ScalarField>],
        alpha_sf: C::ScalarField,
        beta_sf: C::ScalarField,
        gamma_sf: C::ScalarField,
        plonk_z_coeffs: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        // Low degree extend Z.
        let plonk_z_points_8n = fft_with_precomputation_power_of_2(
            &pad_to_8n(&plonk_z_coeffs),
            &self.fft_precomputation_8n,
        );

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
                &local_constant_values,
                &local_wire_values,
                &right_wire_values,
                &below_wire_values,
            );

            // Evaluate the L_1(x) (Z(x) - 1) vanishing term.
            let z_x = plonk_z_points_8n[i];
            let z_gz = plonk_z_points_8n[i_right];
            let vanishing_z_1_term = eval_l_1(self.degree(), x) * (z_x - C::ScalarField::ONE);

            // Evaluate the Z(x) f'(x) - g'(x) Z(g x) term.
            let mut f_prime = C::ScalarField::ONE;
            let mut g_prime = C::ScalarField::ONE;
            for j in 0..NUM_ROUTED_WIRES {
                let wire_value = wire_values_8n[j][i];
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
                constraint_terms,
            ]
            .concat();

            vanishing_points.push(reduce_with_powers(&vanishing_terms, alpha_sf));
        }

        ifft_with_precomputation_power_of_2(&vanishing_points, &self.fft_precomputation_8n)
    }

    /// Open each polynomial at the given point, `zeta`.
    fn open_all_polynomials(
        &self,
        wire_coeffs: &Vec<Vec<C::ScalarField>>,
        plonk_z_coeffs: &Vec<C::ScalarField>,
        plonk_t_coeffs: &Vec<Vec<C::ScalarField>>,
        zeta: C::ScalarField,
    ) -> OpeningSet<C::ScalarField> {
        let powers_of_zeta = powers(zeta, self.degree());

        OpeningSet {
            o_constants: eval_coeffs(&self.constants_coeffs, &powers_of_zeta),
            o_plonk_sigmas: eval_coeffs(&self.s_sigma_coeffs, &powers_of_zeta),
            o_wires: eval_coeffs(&wire_coeffs, &powers_of_zeta),
            o_plonk_z: C::ScalarField::inner_product(&plonk_z_coeffs, &powers_of_zeta),
            o_plonk_t: eval_coeffs(&plonk_t_coeffs, &powers_of_zeta),
        }
    }

    pub fn generate_witness(
        &self,
        inputs: PartialWitness<C::ScalarField>,
    ) -> Witness<C::ScalarField> {
        let start = Instant::now();

        // Index generator indices by their dependencies.
        let mut generator_indices_by_deps: HashMap<Target, Vec<usize>> = HashMap::new();
        for (i, generator) in self.generators.iter().enumerate() {
            for dep in generator.dependencies() {
                generator_indices_by_deps
                    .entry(dep)
                    .or_insert_with(|| Vec::new())
                    .push(i);
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
                let generator: &dyn WitnessGenerator<C::ScalarField> =
                    self.generators[generator_idx].borrow();
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
                let affected_generator_indices = generator_indices_by_deps
                    .get(&target)
                    .unwrap_or(&no_indices);

                for &generator_idx in affected_generator_indices {
                    // If this generator is not already pending or completed, and its dependencies
                    // are all satisfied, then add it as a pending generator.
                    let generator: &dyn WitnessGenerator<C::ScalarField> =
                        self.generators[generator_idx].borrow();
                    if !pending_generator_indices.contains(&generator_idx)
                        && !completed_generator_indices.contains(&generator_idx)
                        && witness.contains_all_targets(&generator.dependencies())
                    {
                        pending_generator_indices.insert(generator_idx);
                    }
                }
            }
        }

        debug_assert_eq!(
            completed_generator_indices.len(),
            self.generators.len(),
            "Only {} of {} generators could be run",
            completed_generator_indices.len(),
            self.generators.len()
        );

        println!("Witness generation took {}s", start.elapsed().as_secs_f32());
        Witness::from_partial(&witness, self.degree())
    }

    /// For the given set of targets, find any copy constraints involving those targets and populate
    /// the witness with copies as needed.
    fn generate_copies(
        &self,
        witness: &PartialWitness<C::ScalarField>,
        targets: &[Target],
    ) -> PartialWitness<C::ScalarField> {
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

#[cfg(test)]
mod tests {
    use crate::{CircuitBuilder, Curve, Field, PartialWitness, Tweedledee, Tweedledum};

    #[test]
    fn test_generate_proof() {
        let mut builder = CircuitBuilder::<Tweedledee>::new(128);
        let t = builder.add_virtual_target();
        builder.assert_zero(t);
        let mut partial_witness = PartialWitness::new();
        partial_witness.set_target(t, <Tweedledee as Curve>::ScalarField::ZERO);
        let circuit = builder.build();
        let witness = circuit.generate_witness(partial_witness);
        let _proof = circuit.generate_proof::<Tweedledum>(witness);
    }

    #[test]
    fn test_generate_sum() {
        let mut builder = CircuitBuilder::<Tweedledee>::new(128);
        let t1 = builder.add_virtual_target();
        let t2 = builder.add_virtual_target();
        let s = builder.add(t1, t2);
        builder.assert_zero(s);
        let mut partial_witness = PartialWitness::new();
        partial_witness.set_target(t1, <Tweedledee as Curve>::ScalarField::ZERO);
        partial_witness.set_target(t2, <Tweedledee as Curve>::ScalarField::ZERO);
        let circuit = builder.build();
        let witness = circuit.generate_witness(partial_witness);
        let _proof = circuit.generate_proof::<Tweedledum>(witness);
    }
}
