use crate::{blake_hash_usize_to_curve, Curve, Field, ProjectivePoint, plonk_challenger::Challenger, msm_execute, msm_precompute, util::log2_ceil};
use anyhow::{ensure, Result};
use rand::Rng;
use std::marker::PhantomData;

trait PolynomialCommitment {
    type Parameters;
    type Field: Field;
    type Commitment;
    type Opening;

    fn generate_parameters<R: Rng>(r: &mut R, max_degree: usize) -> Self::Parameters;

    fn commit(params: &Self::Parameters, coeffs: &[Self::Field], randomness: &Self::Field) -> Result<Self::Commitment>;

    fn open(
        params: &Self::Parameters,
        coeffs: &[Self::Field],
        point: &Self::Field,
        commitment: &Self::Commitment,
    ) -> (Self::Field, Self::Opening);

    fn verify(
        params: &Self::Parameters,
        commitment: &Self::Commitment,
        point: &Self::Field,
        proof: &Self::Opening,
    ) -> bool;
}

struct DLPC<C: Curve> {
    phantom_curve: PhantomData<C>,
}

struct DLPCOpening<C: Curve>(Vec<ProjectivePoint<C>>, Vec<ProjectivePoint<C>>, ProjectivePoint<C>, C::ScalarField);

impl<C: Curve> DLPC<C> {
    const SECURITY_BITS: usize = 128;
}

impl<C: Curve> PolynomialCommitment for DLPC<C> {

    type Field = C::ScalarField;
    type Parameters = (Vec<ProjectivePoint<C>>, ProjectivePoint<C>);
    type Commitment = ProjectivePoint<C>;
    type Opening = DLPCOpening<C>;

    fn generate_parameters<R: Rng>(r: &mut R, max_degree: usize) -> Self::Parameters {
        ((0..max_degree)
            .map(|_| blake_hash_usize_to_curve::<C>(r.next_u64() as usize).to_projective())
            .collect(), blake_hash_usize_to_curve::<C>(r.next_u64() as usize).to_projective())
    }

    fn commit(params: &Self::Parameters, coeffs: &[Self::Field], randomness: &Self::Field) -> Result<Self::Commitment> {
        ensure!(
            coeffs.len() > params.0.len(),
            "Polynomial degree is too large"
        );
        let gs = [&params.0[..], &[params.1]].concat();
        let precomputation = msm_precompute(&gs, 8);
        Ok(msm_execute(&precomputation, &[coeffs, &[*randomness]].concat()))
    }

    fn open(params: &Self::Parameters, coeffs: &[Self::Field], point: &Self::Field, commitment: &Self::Commitment) -> (Self::Field, Self::Opening) {
        let value = evaluate(coeffs, *point);
        let challenger = Challenger::<Self::Field>::new(Self::SECURITY_BITS);
        challenger.observe_proj_point_other_curve(*commitment);
        challenger.observe_elements(&[*point, value]);
        let chal = challenger.get_challenge();
        let h = C::convert(chal) * params.1;
        let degree = coeffs.len();
        assert!((degree+1).is_power_of_two());
        let z_pow_vec = vec![Self::Field::ONE]; 
        let cs = coeffs.to_vec();
        let gs = params.0.to_vec();
        let mut ls = Vec::new();
        let mut rs = Vec::new();
        (1..=degree).for_each(|_| {z_pow_vec.push(*z_pow_vec.last().unwrap()*(*point))});
        let mut n = degree + 1;
        for i in 0..log2_ceil(degree+1) {
            let mut s_l = &gs[..n/2].to_vec();
            let mut s_r = &gs[n/2..].to_vec();
            s_l.push(h);
            s_r.push(h);
            let l = Self::commit(&s_l, coeffs, randomness);
        }
        todo!()
    }
}


struct KZG10;

impl KZG10 {
    /// Taken from the scipr-lab implementation here: https://github.com/scipr-lab/poly-commit/blob/master/src/kzg10/mod.rs
    ///
    /// For now, only produces G1 points, for a polynomial of degree max_degree
    pub fn setup<C: Curve>(
        generator: ProjectivePoint<C>,
        max_degree: u64,
    ) -> Vec<ProjectivePoint<C>> {
        let alpha = C::ScalarField::rand();

        let mut powers_of_alpha = vec![C::ScalarField::ONE];
        let mut cur = alpha;
        for _ in 0..max_degree {
            powers_of_alpha.push(cur);
            cur = cur * alpha;
        }

        //TODO: implement fixed-base scalar multiplication

        let mut projective_params = Vec::new();

        projective_params.push(generator);

        for pow in powers_of_alpha {
            // TODO: Investigate error.
            // projective_params.push(pow * generator);
        }

        //TODO: normalize projective coordinates and pass them as affine points

        projective_params
    }
}

fn evaluate<F: Field>(coeffs: &[F], z: F) -> F {
    let mut z_pow = F::ONE;
    coeffs.iter().fold(F::ZERO, |acc, c| {
        let tmp = c * z_pow;
        z_pow = z_pow * z;
        tmp
    })
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_setup() {
        //let generator_x = Bls12Base { limbs: [] };
        //let generator = G1ProjectivePoint {x: }
    }
}