use anyhow::Result;

use crate::{AffinePoint, AffinePointTarget, Curve, Field, PartialWitness, Target};

pub struct Proof<C: Curve> {
    /// A commitment to each wire polynomial.
    pub c_wires: Vec<AffinePoint<C>>,
    /// A commitment to Z, in the context of the permutation argument.
    pub c_plonk_z: AffinePoint<C>,
    /// A commitment to the quotient polynomial.
    pub c_plonk_t: Vec<AffinePoint<C>>,

    /// The opening of each polynomial at each `PublicInputGate` index.
    pub o_public_inputs: Vec<OpeningSet<C::ScalarField>>,
    /// The opening of each polynomial at `zeta`.
    pub o_local: OpeningSet<C::ScalarField>,
    /// The opening of each polynomial at `g * zeta`.
    pub o_right: OpeningSet<C::ScalarField>,
    /// The opening of each polynomial at `g^65 * zeta`.
    pub o_below: OpeningSet<C::ScalarField>,

    /// L in the Halo reduction.
    pub halo_l: Vec<AffinePoint<C>>,
    /// R in the Halo reduction.
    pub halo_r: Vec<AffinePoint<C>>,
    /// The purported value of G, i.e. <s, G>, in the context of Halo.
    pub halo_g: AffinePoint<C>,
    /// The two scalar used in the final Schnorr protocol of the Halo opening proof.
    pub schnorr_proof: (C::ScalarField, C::ScalarField),
}

impl<C: Curve> Proof<C> {
    pub fn all_opening_sets(&self) -> Vec<OpeningSet<C::ScalarField>> {
        [
            self.o_public_inputs.as_slice(),
            &[
                self.o_local.clone(),
                self.o_right.clone(),
                self.o_below.clone(),
            ],
        ]
        .concat()
    }
}

pub struct ProofTarget {
    /// A commitment to each wire polynomial.
    pub c_wires: Vec<AffinePointTarget>,
    /// A commitment to Z, in the context of the permutation argument.
    pub c_plonk_z: AffinePointTarget,
    /// A commitment to the quotient polynomial.
    pub c_plonk_t: Vec<AffinePointTarget>,

    /// The opening of each polynomial at each `PublicInputGate` index.
    pub o_public_inputs: Vec<OpeningSetTarget>,
    /// The opening of each polynomial at `zeta`.
    pub o_local: OpeningSetTarget,
    /// The opening of each polynomial at `g * zeta`.
    pub o_right: OpeningSetTarget,
    /// The opening of each polynomial at `g^65 * zeta`.
    pub o_below: OpeningSetTarget,

    /// L_i in the Halo reduction.
    pub halo_l_i: Vec<AffinePointTarget>,
    /// R_i in the Halo reduction.
    pub halo_r_i: Vec<AffinePointTarget>,
    /// The purported value of G, i.e. <s, G>, in the context of Halo.
    pub halo_g: AffinePointTarget,
}

impl ProofTarget {
    /// `log_2(d)`, where `d` is the degree of the proof being verified.
    fn degree_pow(&self) -> usize {
        self.halo_l_i.len()
    }

    pub fn all_opening_sets(&self) -> Vec<OpeningSetTarget> {
        [
            self.o_public_inputs.as_slice(),
            &[
                self.o_local.clone(),
                self.o_right.clone(),
                self.o_below.clone(),
            ],
        ]
        .concat()
    }

    pub fn all_opening_targets(&self) -> Vec<Target> {
        let targets_2d: Vec<Vec<Target>> = self
            .all_opening_sets()
            .into_iter()
            .map(|set| set.to_vec())
            .collect();
        targets_2d.concat()
    }

    pub fn populate_witness<C: Curve>(
        &self,
        witness: &mut PartialWitness<C::BaseField>,
        values: Proof<C>,
    ) {
        witness.set_point_targets(&self.c_wires, &values.c_wires);
        witness.set_point_target(self.c_plonk_z, values.c_plonk_z);
        witness.set_point_targets(&self.c_plonk_t, &values.c_plonk_t);

        debug_assert_eq!(self.o_public_inputs.len(), values.o_public_inputs.len());
        for (o_pi_targets, o_pi_values) in self.o_public_inputs.iter().zip(values.o_public_inputs) {
            o_pi_targets.populate_witness(witness, o_pi_values);
        }

        self.o_local.populate_witness(witness, values.o_local);
        self.o_right.populate_witness(witness, values.o_right);
        self.o_below.populate_witness(witness, values.o_below);

        witness.set_point_targets(&self.halo_l_i, &values.halo_l);
        witness.set_point_targets(&self.halo_r_i, &values.halo_r);
        witness.set_point_target(self.halo_g, values.halo_g);
    }
}

/// The opening of each Plonk polynomial at a particular point.
#[derive(Clone, Debug)]
pub struct OpeningSet<F: Field> {
    /// The purported opening of each constant polynomial.
    pub o_constants: Vec<F>,
    /// The purported opening of each S_sigma polynomial in the context of Plonk's permutation argument.
    pub o_plonk_sigmas: Vec<F>,
    /// The purported opening of each wire polynomial.
    pub o_wires: Vec<F>,
    /// The purported opening of `Z`.
    pub o_plonk_z: F,
    /// The purported opening of `t`.
    pub o_plonk_t: Vec<F>,
}

impl<F: Field> OpeningSet<F> {
    pub fn to_vec(&self) -> Vec<F> {
        [
            self.o_constants.as_slice(),
            self.o_plonk_sigmas.as_slice(),
            self.o_wires.as_slice(),
            &[self.o_plonk_z],
            self.o_plonk_t.as_slice(),
        ]
        .concat()
    }
}

/// The opening of each Plonk polynomial at a particular point.
#[derive(Clone)]
pub struct OpeningSetTarget {
    /// The purported opening of each constant polynomial.
    pub o_constants: Vec<Target>,
    /// The purported opening of each S_sigma polynomial in the context of Plonk's permutation argument.
    pub o_plonk_sigmas: Vec<Target>,
    /// The purported opening of each wire polynomial.
    pub o_wires: Vec<Target>,
    /// The purported opening of `Z`.
    pub o_plonk_z: Target,
    /// The purported opening of `t`.
    pub o_plonk_t: Vec<Target>,
}

impl OpeningSetTarget {
    pub fn to_vec(&self) -> Vec<Target> {
        [
            self.o_constants.as_slice(),
            self.o_plonk_sigmas.as_slice(),
            self.o_wires.as_slice(),
            &[self.o_plonk_z],
            self.o_plonk_t.as_slice(),
        ]
        .concat()
    }

    pub fn populate_witness<InnerBF: Field, InnerSF: Field>(
        &self,
        witness: &mut PartialWitness<InnerBF>,
        values: OpeningSet<InnerSF>,
    ) -> Result<()> {
        // TODO: We temporarily assume that each opened value fits in both fields.
        witness.set_targets(
            &self.o_constants,
            &InnerSF::try_convert_all::<InnerBF>(&values.o_constants)?,
        );
        witness.set_targets(
            &self.o_plonk_sigmas,
            &InnerSF::try_convert_all::<InnerBF>(&values.o_plonk_sigmas)?,
        );
        witness.set_targets(
            &self.o_wires,
            &InnerSF::try_convert_all::<InnerBF>(&values.o_wires)?,
        );
        witness.set_target(self.o_plonk_z, values.o_plonk_z.try_convert::<InnerBF>()?);
        witness.set_targets(
            &self.o_plonk_t,
            &InnerSF::try_convert_all::<InnerBF>(&values.o_plonk_t)?,
        );
        Ok(())
    }
}
