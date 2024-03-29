use std::borrow::Borrow;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::time::Instant;

use anyhow::Result;
use rayon::prelude::*;

use crate::halo::batch_opening_proof;
use crate::partition::{get_subgroup_shift, TargetPartitions};
use crate::plonk_challenger::Challenger;
use crate::plonk_proof::{OldProof, Proof};
use crate::plonk_util::{commit_polynomials, eval_l_1, eval_polys, eval_zero_poly, pad_to_8n, permutation_polynomial, polynomials_to_values_padded, powers, reduce_with_powers, scale_polynomials, values_to_polynomials};
use crate::poly_commit::PolynomialCommitment;
use crate::polynomial::Polynomial;
use crate::target::Target;
use crate::util::{ceil_div_usize, log2_strict};
use crate::witness::{PartialWitness, Witness, WitnessGenerator};
use crate::{evaluate_all_constraints, fft_with_precomputation_power_of_2, AffinePoint, FftPrecomputation, Field, HaloCurve, MsmPrecomputation, OpeningSet, VerificationKey};

pub(crate) const NUM_WIRES: usize = 9;
pub(crate) const NUM_ROUTED_WIRES: usize = 6;
pub(crate) const NUM_ADVICE_WIRES: usize = NUM_WIRES - NUM_ROUTED_WIRES;
pub(crate) const NUM_CONSTANTS: usize = 6;
pub(crate) const GRID_WIDTH: usize = 65;
// This is currently dominated by Base4SumGate. It has degree-4n constraints, and its prefix is 4
// bits long, so its filtered constraints are degree-8n. Dividing by Z_H makes t degree-7n.
pub(crate) const QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER: usize = 7;

/// Contains all data needed to generate and/or verify proofs.
pub struct Circuit<C: HaloCurve> {
    pub security_bits: usize,
    pub num_public_inputs: usize,
    pub num_gates_without_pis: usize,
    pub gate_constants: Vec<Vec<C::ScalarField>>,
    pub routing_target_partitions: TargetPartitions<C::ScalarField>,
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
    pub constant_polynomials: Vec<Polynomial<C::ScalarField>>,
    /// Each constant polynomial, in point-value form, low-degree extended to be degree 8n.
    pub constants_8n: Vec<Vec<C::ScalarField>>,
    /// A commitment to each constant polynomial.
    pub c_constants: Vec<PolynomialCommitment<C>>,
    /// Each permutation polynomial, in coefficient form.
    pub s_sigma_polynomials: Vec<Polynomial<C::ScalarField>>,
    /// Each permutation polynomial, low-degree extended to be degree 8n.
    pub s_sigma_values_8n: Vec<Vec<C::ScalarField>>,
    /// A commitment to each permutation polynomial.
    pub c_s_sigmas: Vec<PolynomialCommitment<C>>,
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
        witness: &Witness<C::ScalarField>,
        old_proofs: &[OldProof<C>],
        blinding_commitments: bool,
    ) -> Result<Proof<C>> {
        let mut challenger = Challenger::new(self.security_bits);

        // Convert the witness both to coefficient form and a degree-8n LDE.
        let wire_values_by_wire_index = &witness.transpose();
        let wire_polynomials =
            values_to_polynomials(wire_values_by_wire_index, &self.fft_precomputation_n);
        let wire_values_8n =
            polynomials_to_values_padded(&wire_polynomials, &self.fft_precomputation_8n);

        // Commit to the wire polynomials.
        let c_wires = commit_polynomials(
            &wire_polynomials,
            &self.pedersen_g_msm_precomputation,
            self.pedersen_h,
            blinding_commitments,
        );

        let num_public_input_gates = ceil_div_usize(self.num_public_inputs, NUM_WIRES);
        // Compute the wire coefficients when the public input gates are set to zero.
        let wire_polynomials_no_pis = {
            let mut wire_values_by_wire_no_pis = wire_values_by_wire_index.clone();
            wire_values_by_wire_no_pis.iter_mut().for_each(|w| {
                for i in 0..num_public_input_gates {
                    // Set the wire value at the public input gate to zero.
                    w[self.num_gates_without_pis + 2 * i] = C::ScalarField::ZERO;
                }
            });
            values_to_polynomials(&wire_values_by_wire_no_pis, &self.fft_precomputation_n)
        };

        // Generate a random beta and gamma from the transcript.
        challenger
            .observe_affine_points(&PolynomialCommitment::commitments_to_affine_vec(&c_wires));
        let (beta_bf, gamma_bf) = challenger.get_2_challenges();
        let beta_sf = beta_bf.try_convert::<C::ScalarField>()?;
        let gamma_sf = gamma_bf.try_convert::<C::ScalarField>()?;

        let plonk_z_points_n = permutation_polynomial(
            self.degree(),
            &self.subgroup_n,
            witness,
            &self.s_sigma_values_8n,
            beta_sf,
            gamma_sf,
        );
        // Commit to Z.
        let plonk_z_polynomial =
            Polynomial::from_evaluations(&plonk_z_points_n, &self.fft_precomputation_n);
        let c_plonk_z = plonk_z_polynomial.commit(
            &self.pedersen_g_msm_precomputation,
            self.pedersen_h,
            blinding_commitments,
        );

        // Generate a random alpha from the transcript.
        challenger.observe_affine_point(c_plonk_z.to_affine());
        let alpha_bf = challenger.get_challenge();
        let alpha_sf = alpha_bf.try_convert::<C::ScalarField>()?;

        // Generate the vanishing polynomial.
        let vanishing_poly = self.vanishing_poly::<InnerC>(
            &wire_values_8n,
            alpha_sf,
            beta_sf,
            gamma_sf,
            &plonk_z_polynomial[..],
        );

        if cfg!(debug_assertions) {
            // Check that the vanishing polynomial indeed vanishes.
            self.subgroup_n.iter().enumerate().for_each(|(i, &x)| {
                assert!(
                    vanishing_poly.eval(x).is_zero(),
                    "{}-th gate constraints are not satisfied",
                    i
                );
            });
        }

        // Compute the quotient polynomial, t(x) = vanishing(x) / Z_H(x).
        let mut plonk_t_poly = vanishing_poly.divide_by_z_h(self.degree());

        if cfg!(debug_assertions) {
            // Check that division was performed correctly by evaluating at a random point.
            let x = C::ScalarField::rand();
            assert_eq!(
                plonk_t_poly.eval(x),
                vanishing_poly.eval(x) / eval_zero_poly(self.degree(), x),
            );
        }

        // Pad the coefficients to 7n.
        plonk_t_poly.pad(QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER * self.degree());

        // Split t into degree-n chunks.
        let plonk_t_poly_chunks = plonk_t_poly
            .coeffs()
            .chunks(self.degree())
            .map(|chunk| Polynomial::from(chunk.to_vec()))
            .collect::<Vec<_>>();

        // Commit to the quotient polynomial.
        let c_plonk_t = commit_polynomials(
            &plonk_t_poly_chunks,
            &self.pedersen_g_msm_precomputation,
            self.pedersen_h,
            blinding_commitments,
        );

        // Combine the coefficients in `wires_coeffs_no_pis` using a linear combination weighted by `alpha`.
        let vanishing_pis_poly =
            scale_polynomials(wire_polynomials_no_pis, alpha_sf, self.degree());
        // `vanishing_pis_poly` vanishes at the public input gates. It is thus divisible by the vanishing
        // polynomial at the public input gates. The quotient is computed here.
        let pis_quotient_poly = {
            // The vanishing polynomial of a set `S` is `prod_{s \in S} (X-s)`.
            // TODO: Faster implementation.
            let pis_quotient_denominator = (0..num_public_input_gates).fold(
                Polynomial::from(vec![C::ScalarField::ONE]),
                |acc, i| {
                    let mut ans =
                        acc.mul(&vec![-self.subgroup_n[self.num_gates_without_pis + 2 * i], C::ScalarField::ONE].into());
                    ans.trim();
                    ans
                },
            );
            let mut ans = vanishing_pis_poly
                .polynomial_division(&pis_quotient_denominator)
                .0;
            if cfg!(debug_assertions) {
                // Check that division was performed correctly by evaluating at a random point.
                let x = C::ScalarField::rand();
                assert_eq!(
                    ans.eval(x),
                    vanishing_pis_poly.eval(x) / pis_quotient_denominator.eval(x)
                );
            }
            ans.pad(self.degree());
            ans
        };
        // Commit to the public inputs quotient polynomial.
        let c_pis_quotient = pis_quotient_poly.commit(
            &self.pedersen_g_msm_precomputation,
            self.pedersen_h,
            blinding_commitments,
        );

        let public_inputs = (0..self.num_public_inputs)
            .map(|i| wire_values_by_wire_index[i % NUM_WIRES][self.num_gates_without_pis + 2 * (i / NUM_WIRES)])
            .collect::<Vec<_>>();

        // Observe the `t` polynomial commitment.
        challenger
            .observe_affine_points(&PolynomialCommitment::commitments_to_affine_vec(&c_plonk_t));
        // If the proof doesn't output public inputs, observe the `pis_quotient` polynomial commitment and the public inputs.
        challenger.observe_affine_point(c_pis_quotient.to_affine());
        // Observe the public inputs
        challenger.observe_elements(
            &C::ScalarField::try_convert_all(&public_inputs)
                .expect("Public inputs should fit in both fields"),
        );
        // Observe the old proofs' `G` points.
        old_proofs
            .iter()
            .for_each(|old_proof| challenger.observe_affine_point(old_proof.halo_g));
        // Generate a random `zeta` from the transcript.
        let zeta_bf = challenger.get_challenge();
        let zeta_sf =
            C::try_convert_b2s(zeta_bf).expect("should fit in both fields with high probability");

        // Open all polynomials at zeta, zeta * g, and zeta * g^65.
        let o_local = self.open_all_polynomials(
            &wire_polynomials,
            &plonk_z_polynomial,
            &plonk_t_poly_chunks,
            old_proofs,
            &pis_quotient_poly,
            zeta_sf,
        );
        let o_right = self.open_all_polynomials(
            &wire_polynomials,
            &plonk_z_polynomial,
            &plonk_t_poly_chunks,
            old_proofs,
            &pis_quotient_poly,
            zeta_sf * self.subgroup_generator_n,
        );
        let o_below = self.open_all_polynomials(
            &wire_polynomials,
            &plonk_z_polynomial,
            &plonk_t_poly_chunks,
            old_proofs,
            &pis_quotient_poly,
            zeta_sf * self.subgroup_generator_n.exp_usize(GRID_WIDTH),
        );

        // Get a list of all opened values, to append to the transcript.
        let all_opening_sets: Vec<OpeningSet<C::ScalarField>> =
            vec![o_local.clone(), o_right.clone(), o_below.clone()];
        let all_opened_values_sf: Vec<C::ScalarField> = all_opening_sets
            .iter()
            .map(|os| os.to_vec())
            .collect::<Vec<_>>()
            .concat();
        let all_opened_values_bf: Vec<_> = all_opened_values_sf
            .into_iter()
            .map(|f| {
                // TODO: Fix this, this can regularly fail if a public input is not in range for example.
                C::try_convert_s2b(f)
                    .expect("For now, we assume that all opened values fit in both fields")
            })
            .collect();

        // Generate random v, u, and x from the transcript.
        challenger.observe_elements(&all_opened_values_bf);
        let (v_bf, u_bf, u_scaling_bf) = challenger.get_3_challenges();
        let v_sf = v_bf.try_convert::<C::ScalarField>()?;
        let u_sf = u_bf.try_convert::<C::ScalarField>()?;
        let u_scaling_sf = u_scaling_bf.try_convert::<C::ScalarField>()?;

        let old_proofs_polys = old_proofs
            .iter()
            .map(|p| Polynomial::from(p.coeffs()))
            .collect::<Vec<_>>();
        let all_coeffs = [
            self.constant_polynomials.clone(),
            self.s_sigma_polynomials.clone(),
            wire_polynomials,
            vec![plonk_z_polynomial],
            plonk_t_poly_chunks,
            old_proofs_polys,
            vec![pis_quotient_poly],
        ]
        .concat();

        let commitments = [
            self.c_constants.clone(),
            self.c_s_sigmas.clone(),
            c_wires.clone(),
            vec![c_plonk_z],
            c_plonk_t.clone(),
            old_proofs
                .iter()
                .map(|old| old.halo_g.into())
                .collect::<Vec<_>>(),
            vec![c_pis_quotient],
        ]
        .concat();

        let opening_points = vec![
            zeta_sf,
            zeta_sf * self.subgroup_generator_n,
            zeta_sf * self.subgroup_generator_n.exp_usize(GRID_WIDTH),
        ];

        let halo_proof = batch_opening_proof(
            &all_coeffs.iter().map(|c| &c[..]).collect::<Vec<_>>(),
            &commitments,
            &opening_points,
            &self.pedersen_g,
            self.pedersen_h.to_projective(),
            self.u,
            u_sf,
            v_sf,
            u_scaling_sf,
            self.degree(),
            self.security_bits,
            &mut challenger,
        )?;

        Ok(Proof {
            c_wires: c_wires.iter().map(|c| c.to_affine()).collect(),
            c_plonk_z: c_plonk_z.to_affine(),
            c_plonk_t: c_plonk_t.iter().map(|c| c.to_affine()).collect(),
            c_pis_quotient: c_pis_quotient.to_affine(),
            o_local,
            o_right,
            o_below,
            halo_g: halo_proof.halo_g,
            halo_l: halo_proof.halo_l,
            halo_r: halo_proof.halo_r,
            schnorr_proof: halo_proof.schnorr_proof,
        })
    }

    fn vanishing_poly<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &self,
        wire_values_8n: &[Vec<C::ScalarField>],
        alpha_sf: C::ScalarField,
        beta_sf: C::ScalarField,
        gamma_sf: C::ScalarField,
        plonk_z_coeffs: &[C::ScalarField],
    ) -> Polynomial<C::ScalarField> {
        let degree = self.degree();
        let k_is = (0..NUM_ROUTED_WIRES)
            .map(get_subgroup_shift::<C::ScalarField>)
            .collect::<Vec<_>>();
        // Low degree extend Z.
        let plonk_z_points_8n = fft_with_precomputation_power_of_2(
            &pad_to_8n(plonk_z_coeffs),
            &self.fft_precomputation_8n,
        );

        // We will evaluate the vanishing polynomial at 8n points, then interpolate.
        let vanishing_points = self
            .subgroup_8n
            .par_iter()
            .enumerate()
            .map(|(i, &x)| {
                // Load the constant polynomials' values at x.
                let mut local_constant_values = Vec::new();
                for j in 0..NUM_CONSTANTS {
                    local_constant_values.push(self.constants_8n[j][i]);
                }

                // Load the wire polynomials' values at x, g x (the "right" position), and g^WIDTH x
                // (the "below" position). Note that a shift of 1 in the degree-n subgroup corresponds
                // to a shift of 8 in the degree-8n subgroup.
                let i_right = (i + 8) % (8 * degree);
                let i_below = (i + 8 * GRID_WIDTH) % (8 * degree);
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
                let vanishing_z_1_term = eval_l_1(degree, x) * (z_x - C::ScalarField::ONE);

                // Evaluate the Z(x) f'(x) - g'(x) Z(g x) term.
                let mut f_prime = C::ScalarField::ONE;
                let mut g_prime = C::ScalarField::ONE;
                for j in 0..NUM_ROUTED_WIRES {
                    let wire_value = wire_values_8n[j][i];
                    let k_i = k_is[j];
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

                reduce_with_powers(&vanishing_terms, alpha_sf)
            })
            .collect::<Vec<_>>();

        Polynomial::from_evaluations(&vanishing_points, &self.fft_precomputation_8n)
    }

    /// Open each polynomial at the given point, `zeta`.
    fn open_all_polynomials(
        &self,
        wire_poly: &[Polynomial<C::ScalarField>],
        plonk_z_poly: &Polynomial<C::ScalarField>,
        plonk_t_polys: &[Polynomial<C::ScalarField>],
        old_proofs: &[OldProof<C>],
        pi_quotient_poly: &Polynomial<C::ScalarField>,
        zeta: C::ScalarField,
    ) -> OpeningSet<C::ScalarField> {
        let powers_of_zeta = powers(zeta, self.degree());

        OpeningSet {
            o_constants: eval_polys(&self.constant_polynomials, &powers_of_zeta),
            o_plonk_sigmas: eval_polys(&self.s_sigma_polynomials, &powers_of_zeta),
            o_wires: eval_polys(wire_poly, &powers_of_zeta),
            o_plonk_z: plonk_z_poly.eval_from_power(&powers_of_zeta),
            o_plonk_t: eval_polys(plonk_t_polys, &powers_of_zeta),
            o_pi_quotient: pi_quotient_poly.eval_from_power(&powers_of_zeta),
            o_old_proofs: old_proofs
                .iter()
                .map(|p| p.evaluate_g(zeta))
                .collect::<Vec<_>>(),
        }
    }

    /// Generates a `PartialWitness`, which maps `Target`s to their values. Although
    /// `PartialWitness` is designed as a sparse representation, the result here should have an
    /// entry for every target in the circuit.
    pub fn generate_partial_witness(
        &self,
        inputs: PartialWitness<C::ScalarField>,
    ) -> PartialWitness<C::ScalarField> {
        let start = Instant::now();

        // Index generator indices by their dependencies.
        let mut generator_indices_by_deps: HashMap<Target<C::ScalarField>, Vec<usize>> =
            HashMap::new();
        for (i, generator) in self.generators.iter().enumerate() {
            for dep in generator.dependencies() {
                generator_indices_by_deps
                    .entry(dep)
                    .or_insert_with(Vec::new)
                    .push(i);
            }
        }

        // We start with the inputs as our witness, and execute any copy constraints.
        let mut witness = inputs;

        // Replace public inputs targets by their corresponding wires in the circuit.
        witness.replace_public_inputs(self.num_gates_without_pis);

        let mut copy_result = self.generate_copies(&witness, &witness.all_populated_targets());
        copy_result.copy_buffer_to_pi_gate(self.num_gates_without_pis);
        witness.extend(copy_result);

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
            let mut populated_targets: Vec<Target<C::ScalarField>> = Vec::new();

            for &generator_idx in &pending_generator_indices {
                let generator: &dyn WitnessGenerator<C::ScalarField> =
                    self.generators[generator_idx].borrow();
                let result = generator.generate(&self.gate_constants, &witness);
                populated_targets.extend(result.all_populated_targets());
                witness.extend(result);
                completed_generator_indices.insert(generator_idx);
            }

            let mut copy_result = self.generate_copies(&witness, &populated_targets);
            copy_result.copy_buffer_to_pi_gate(self.num_gates_without_pis);
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

        // TODO: Fix this.
        // debug_assert_eq!(
        //     completed_generator_indices.len(),
        //     self.generators.len(),
        //     "Only {} of {} generators could be run",
        //     completed_generator_indices.len(),
        //     self.generators.len()
        // );

        info!("Witness generation took {}s", start.elapsed().as_secs_f32());
        witness
    }

    pub fn generate_witness(
        &self,
        inputs: PartialWitness<C::ScalarField>,
    ) -> Witness<C::ScalarField> {
        let partial_witness = self.generate_partial_witness(inputs);
        Witness::from_partial(&partial_witness, self.degree())
    }

    /// For the given set of targets, find any copy constraints involving those targets and populate
    /// the witness with copies as needed.
    fn generate_copies(
        &self,
        witness: &PartialWitness<C::ScalarField>,
        targets: &[Target<C::ScalarField>],
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

    pub fn to_vk(&self) -> VerificationKey<C> {
        VerificationKey {
            c_constants: self
                .c_constants
                .iter()
                .map(|c| c.to_affine())
                .collect::<Vec<_>>(),
            c_s_sigmas: self
                .c_s_sigmas
                .iter()
                .map(|c| c.to_affine())
                .collect::<Vec<_>>(),
            degree: self.degree(),
            num_public_inputs: self.num_public_inputs,
            num_gates_without_pis: self.num_gates_without_pis,
            security_bits: self.security_bits,
            pedersen_g_msm_precomputation: Some(self.pedersen_g_msm_precomputation.clone()),
            fft_precomputation: Some(self.fft_precomputation_n.clone()),
        }
    }

    pub fn get_public_inputs(&self, witness: &Witness<C::ScalarField>) -> Vec<C::ScalarField> {
        (0..self.num_public_inputs)
            .map(|i| witness.get_indices(self.num_gates_without_pis + 2 * (i / NUM_WIRES), i % NUM_WIRES))
            .collect()
    }
}

impl<C: HaloCurve> Debug for Circuit<C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Circuit of size {}.", self.degree())
    }
}
