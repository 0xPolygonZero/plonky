use crate::{Target, AffinePoint, Curve, AffinePointTarget, Field};

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
            &[self.o_local.clone(), self.o_right.clone(), self.o_below.clone()],
        ].concat()
    }

    pub fn all_opening_targets(&self) -> Vec<Target> {
        let targets_2d: Vec<Vec<Target>> = self.all_opening_sets().into_iter()
            .map(|set| set.to_vec())
            .collect();
        targets_2d.concat()
    }
}

/// The opening of each Plonk polynomial at a particular point.
#[derive(Clone)]
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
        ].concat()
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
        ].concat()
    }
}