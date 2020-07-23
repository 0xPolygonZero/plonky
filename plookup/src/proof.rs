use crate::openings::PlookupOpenings;
use crate::plookup::SECURITY_BITS;
use anyhow::{anyhow, Result};
use plonky::halo::OpeningProof;
use plonky::plonk_challenger::Challenger;
use plonky::plonk_util::halo_n;
use plonky::{AffinePoint, Curve, Field, HaloCurve, PolynomialCommitment};

// TODO: The verifier should somehow have access to `c_f`, so we can probably remove it from the proof or make it optional.
pub struct PlookupProof<C: HaloCurve> {
    pub c_f: AffinePoint<C>,
    pub c_t: AffinePoint<C>,
    pub c_h1: AffinePoint<C>,
    pub c_h2: AffinePoint<C>,
    pub c_z: AffinePoint<C>,
    pub c_quotient: AffinePoint<C>,
    pub openings: PlookupOpenings<C::ScalarField>,
    pub halo_proof: OpeningProof<C>,
    pub n: usize,
}

impl<C: HaloCurve>
    From<(
        Vec<PolynomialCommitment<C>>,
        PlookupOpenings<C::ScalarField>,
        OpeningProof<C>,
        usize,
    )> for PlookupProof<C>
{
    fn from(
        mut t: (
            Vec<PolynomialCommitment<C>>,
            PlookupOpenings<C::ScalarField>,
            OpeningProof<C>,
            usize,
        ),
    ) -> Self {
        PolynomialCommitment::batch_to_affine(&mut t.0);
        Self {
            c_f: t.0[0].to_affine(),
            c_t: t.0[1].to_affine(),
            c_h1: t.0[2].to_affine(),
            c_h2: t.0[3].to_affine(),
            c_z: t.0[4].to_affine(),
            c_quotient: t.0[5].to_affine(),
            openings: t.1,
            halo_proof: t.2,
            n: t.3,
        }
    }
}

pub struct PlookupProofChallenge<C: Curve> {
    pub beta: C::ScalarField,
    pub gamma: C::ScalarField,
    pub alpha: C::ScalarField,
    pub zeta: C::ScalarField,
    pub v: C::ScalarField,
    pub u: C::ScalarField,
    pub u_scaling: C::ScalarField,
    pub halo_us: Vec<C::ScalarField>,
    pub schnorr_challenge: C::ScalarField,
}

impl<C: HaloCurve> PlookupProof<C> {
    pub fn get_challenges(&self) -> Result<PlookupProofChallenge<C>> {
        let mut challenger = Challenger::new(SECURITY_BITS);
        let error_msg = "Conversion from base to scalar field failed.";
        challenger.observe_affine_points(&[self.c_f, self.c_t, self.c_h1, self.c_h2]);
        let (beta_bf, gamma_bf) = challenger.get_2_challenges();
        let beta = C::try_convert_b2s(beta_bf).map_err(|_| anyhow!(error_msg))?;
        let gamma = C::try_convert_b2s(gamma_bf).map_err(|_| anyhow!(error_msg))?;
        challenger.observe_affine_point(self.c_z);
        let alpha_bf = challenger.get_challenge();
        let alpha = C::try_convert_b2s(alpha_bf).map_err(|_| anyhow!(error_msg))?;
        challenger.observe_affine_point(self.c_quotient);
        let zeta_bf = challenger.get_challenge();
        let zeta = C::try_convert_b2s(zeta_bf).map_err(|_| anyhow!(error_msg))?;
        let openings_bf: Vec<_> = self
            .openings
            .to_vec()
            .into_iter()
            .map(|f| {
                C::try_convert_s2b(f)
                    .expect("For now, we assume that all opened values fit in both fields")
            })
            .collect();
        challenger.observe_elements(&openings_bf);
        let (v_bf, u_bf, u_scaling_bf) = challenger.get_3_challenges();
        let v = C::try_convert_b2s(v_bf).map_err(|_| anyhow!(error_msg))?;
        let u = C::try_convert_b2s(u_bf).map_err(|_| anyhow!(error_msg))?;
        let u_scaling = C::try_convert_b2s(u_scaling_bf).map_err(|_| anyhow!(error_msg))?;

        // Compute IPA challenges.
        let mut halo_us = Vec::new();
        for i in 0..self.halo_proof.halo_l.len() {
            challenger
                .observe_affine_points(&[self.halo_proof.halo_l[i], self.halo_proof.halo_r[i]]);
            let r_bf = challenger.get_challenge();
            let r_sf = r_bf.try_convert::<C::ScalarField>()?;
            let r_bits = &r_sf.to_canonical_bool_vec()[..SECURITY_BITS];
            let u_j_squared = halo_n::<C>(r_bits);
            let u_j = u_j_squared
                .square_root()
                .expect("Prover should have ensured that n(r) is square");
            halo_us.push(u_j);
        }

        // Compute challenge for Schnorr protocol.
        challenger.observe_affine_point(self.halo_proof.schnorr_proof.r);
        let schnorr_challenge_bf = challenger.get_challenge();
        let schnorr_challenge =
            C::try_convert_b2s(schnorr_challenge_bf).map_err(|_| anyhow!(error_msg))?;

        Ok(PlookupProofChallenge {
            beta,
            gamma,
            alpha,
            zeta,
            v,
            u,
            u_scaling,
            halo_us,
            schnorr_challenge,
        })
    }
}
