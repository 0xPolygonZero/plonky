use crate::{G1ProjectivePoint, Bls12Scalar};

pub struct KZG10 {

}

impl KZG10 {

    /// Taken from the scipr-lab implementation here: https://github.com/scipr-lab/poly-commit/blob/master/src/kzg10/mod.rs
    ///
    /// For now, only produces G1 points, for a polynomial of degree max_degree
    pub fn setup (generator: G1ProjectivePoint, max_degree: u64) -> Vec<G1ProjectivePoint> {

        let alpha = Bls12Scalar::rand();

        let mut powers_of_alpha = vec![Bls12Scalar::ONE];
        let mut cur = alpha;
        for _ in 0..max_degree {
            powers_of_alpha.push(cur);
            cur = cur * alpha;
        }

        //TODO: implement fixed-base scalar multiplication

        let mut projective_params = Vec::new();

        projective_params.push(generator);

        for pow in powers_of_alpha {
            projective_params.push(pow * generator);
        }

        //TODO: normalize projective coordinates and pass them as affine points

        projective_params
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use num::BigUint;

    use crate::field::{Bls12Base, Bls12Scalar};
    use crate::G1ProjectivePoint;

    #[test]
    fn test_setup() {
        //let generator_x = Bls12Base { limbs: [] };
        //let generator = G1ProjectivePoint {x: }
    }
}
