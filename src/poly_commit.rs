use crate::MsmPrecomputation;
use crate::plonk_util::pedersen_hash;
use crate::{AffinePoint, Curve, Field, ProjectivePoint};

#[derive(Debug, Copy, Clone)]
/// Represents a curve point either in affine or projective coordinates.
enum CurvePoint<C: Curve> {
    Affine(AffinePoint<C>),
    Projective(ProjectivePoint<C>),
}

#[derive(Debug, Copy, Clone)]
/// A Bulletproof polynomial commitment.
/// `commitment` is the actual commitment, while randomness is the scalar factor used for blinding.
pub struct PolynomialCommitment<C: Curve> {
    commitment: CurvePoint<C>,
    pub randomness: C::ScalarField,
}

impl<C: Curve> PolynomialCommitment<C> {
    /// Creates a polynomial commitment from a vector of coefficients.
    /// If `blinding` is true, a random blinding factor is used. Otherwise, it is set to zero.
    pub(crate) fn coeffs_to_commitment(
        coeffs: &[C::ScalarField],
        msm_precomputation: &MsmPrecomputation<C>,
        blinding_point: AffinePoint<C>,
        blinding: bool,
    ) -> Self {
        let blinding_factor = if blinding {
            C::ScalarField::rand()
        } else {
            C::ScalarField::ZERO
        };
        let proj = pedersen_hash(coeffs, msm_precomputation)
            + C::convert(blinding_factor) * blinding_point.to_projective();
        Self {
            commitment: CurvePoint::Projective(proj),
            randomness: blinding_factor,
        }
    }

    /// Creates a list of polynomial commitments from a list of polynomials in coefficients vector form.
    pub(crate) fn coeffs_vec_to_commitments(
        coefficients_vec: &[Vec<C::ScalarField>],
        msm_precomputation: &MsmPrecomputation<C>,
        blinding_point: AffinePoint<C>,
        blinding: bool,
    ) -> Vec<Self> {
        let mut comms: Vec<_> = coefficients_vec
            .iter()
            .map(|coeffs| {
                Self::coeffs_to_commitment(coeffs, &msm_precomputation, blinding_point, blinding)
            })
            .collect();
        Self::batch_to_affine(&mut comms);
        comms
    }

    /// Changes all commitments in projective coordinates to affine coordinates.
    pub fn batch_to_affine(comms: &mut [Self]) {
        let proj_indices = comms
            .iter()
            .enumerate()
            .filter(|(_i, p)| matches!(p.commitment, CurvePoint::Projective(_)))
            .map(|(i, _p)| i)
            .collect::<Vec<_>>();
        let projs = proj_indices
            .iter()
            .map(|&i| {
                if let CurvePoint::Projective(x) = comms[i].commitment {
                    x
                } else {
                    unreachable!()
                }
            })
            .collect::<Vec<_>>();
        let affs = ProjectivePoint::batch_to_affine(&projs);
        for (&i, a) in proj_indices.iter().zip(affs.iter()) {
            comms[i].commitment = CurvePoint::Affine(*a);
        }
    }

    /// Returns the commitment point in affine coordinates.
    /// If this method is to be used on a list of commitments, some of which are in projective coordinates,
    /// `Self::batch_to_affine` should be run first for better performances.
    pub fn to_affine(&self) -> AffinePoint<C> {
        match self.commitment {
            CurvePoint::Affine(p) => p,
            CurvePoint::Projective(p) => p.to_affine(),
        }
    }
}
