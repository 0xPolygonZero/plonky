use std::ops::{Add, Neg};

use crate::{Bls12Base, Bls12Scalar};
use std::iter::Sum;

// Parameters taken from the implementation of Bls12-377 in Zexe found here:
// https://github.com/scipr-lab/zexe/blob/master/algebra/src/curves/bls12_377/g1.rs

const A: Bls12Base = Bls12Base::ZERO;
const B: Bls12Base = Bls12Base::ONE;

const COFACTOR: &'static [u64] = &[0x0, 0x170b5d4430000000];

const COFACTOR_INV: Bls12Scalar = Bls12Scalar {
    limbs: [
        2013239619100046060,
        4201184776506987597,
        2526766393982337036,
        1114629510922847535,
    ]
};

/// 81937999373150964239938255573465948239988671502647976594219695644855304257327692006745978603320413799295628339695
const G1_GENERATOR_X: Bls12Base = Bls12Base {
    limbs: [2742467569752756724, 14217256487979144792, 6635299530028159197, 8509097278468658840,
        14518893593143693938, 46181716169194829]
};

/// 241266749859715473739788878240585681733927191168601896383759122102112907357779751001206799952863815012735208165030
const G1_GENERATOR_Y: Bls12Base = Bls12Base {
    limbs: [9336971515457667571, 28021381849722296, 18085035374859187530, 14013031479170682136,
        3369780711397861396, 35370409237953649]
};

pub const G1_GENERATOR_AFFINE: G1AffinePoint = G1AffinePoint {
    x: G1_GENERATOR_X,
    y: G1_GENERATOR_Y,
};

pub const G1_GENERATOR_PROJECTIVE: G1ProjectivePoint = G1ProjectivePoint {
    x: G1_GENERATOR_X,
    y: G1_GENERATOR_Y,
    z: Bls12Base::ONE,
};

#[derive(Hash, Copy, Clone, Debug, PartialEq, Eq)]
pub struct G1AffinePoint {
    pub x: Bls12Base,
    pub y: Bls12Base,
}

impl G1AffinePoint {
    pub const ZERO: G1AffinePoint = G1AffinePoint { x: Bls12Base::ZERO, y: Bls12Base::ZERO };

    pub fn is_zero(&self) -> bool {
        *self == G1AffinePoint::ZERO
    }

    pub fn to_projective(&self) -> G1ProjectivePoint {
        G1ProjectivePoint { x: self.x, y: self.y, z: Bls12Base::ONE }
    }
}

#[derive(Hash, Copy, Clone, Debug)]
pub struct G1ProjectivePoint {
    pub x: Bls12Base,
    pub y: Bls12Base,
    pub z: Bls12Base,
}

impl G1ProjectivePoint {
    pub const ZERO: G1ProjectivePoint = G1ProjectivePoint { x: Bls12Base::ZERO, y: Bls12Base::ONE, z: Bls12Base::ZERO };

    pub fn is_zero(&self) -> bool {
        self.z == Bls12Base::ZERO
    }

    pub fn double(&self) -> G1ProjectivePoint {
        let G1ProjectivePoint { x: x1, y: y1, z: z1 } = *self;
        let xx = x1.square();
        let zz = z1.square();
        let w = xx.triple();
        let s = y1 * z1.double();
        let ss = s.square();
        let sss = s * ss;
        let r = y1 * s;
        let rr = r.square();
        let b = (x1 + r).square() - xx - rr;
        let h = w.square() - b.double();
        let x3 = h * s;
        let y3 = w * (b - h) - rr.double();
        let z3 = sss;
        G1ProjectivePoint { x: x3, y: y3, z: z3 }
    }

    pub fn new(x: Bls12Base, y: Bls12Base, z: Bls12Base) -> G1ProjectivePoint {
        debug_assert!(G1ProjectivePoint::is_on_curve(x, y, z) /*&& is_in_subgroup(x, y, z)*/);
        G1ProjectivePoint { x, y, z }
    }

    fn is_on_curve(x: Bls12Base, y: Bls12Base, z: Bls12Base) -> bool {
        if z == Bls12Base::ZERO {
            true
        } else {
            let y = y / z;
            let x = x / z;
            y * y == x * x * x + B
        }
    }

    /*
    fn is_in_subgroup(x: Bls12Base, y: Bls12Base, z: Bls12Base) -> bool {

    }
    */

    pub fn to_affine(&self) -> G1AffinePoint {
        if self.is_zero() {
            G1AffinePoint::ZERO
        } else {
            G1AffinePoint { x: self.x / self.z, y: self.y / self.z }
        }
    }

    pub fn batch_to_affine(proj_points: &[Self]) -> Vec<G1AffinePoint> {
        let n = proj_points.len();
        let zs: Vec<Bls12Base> = proj_points.iter().map(|pp| pp.z).collect();
        let z_inv_opts = Bls12Base::batch_multiplicative_inverse_opt(&zs);

        let mut result = Vec::with_capacity(n);
        for i in 0..n {
            let pp = proj_points[i];
            result.push(match z_inv_opts[i] {
                Some(z_inv) => G1AffinePoint { x: pp.x * z_inv, y: pp.y * z_inv },
                None => G1AffinePoint::ZERO,
            });
        }
        result
    }
}

impl PartialEq for G1ProjectivePoint {
    fn eq(&self, other: &Self) -> bool {
        let self_zero = self.is_zero();
        let other_zero = other.is_zero();
        if self_zero || other_zero {
            return self_zero == other_zero;
        }

        // We want to compare (x1/z1, y1/z1) == (x2/z2, y2/z2).
        // But to avoid field division, it is better to compare (x1*z2, y1*z2) == (x2*z1, y2*z1).
        self.x * other.z == other.x * self.z
            && self.y * other.z == other.y * self.z
    }
}

impl Eq for G1ProjectivePoint {}

impl Neg for G1AffinePoint {
    type Output = G1AffinePoint;

    fn neg(self) -> Self::Output {
        G1AffinePoint { x: self.x, y: -self.y }
    }
}

impl Neg for G1ProjectivePoint {
    type Output = G1ProjectivePoint;

    fn neg(self) -> Self::Output {
        G1ProjectivePoint { x: self.x, y: -self.y, z: self.z }
    }
}

#[cfg(test)]
mod tests {
    use crate::{Bls12Base, Bls12Scalar, G1_GENERATOR_PROJECTIVE, G1ProjectivePoint};

    #[test]
    fn test_naive_multiplication() {
        let g = G1_GENERATOR_PROJECTIVE;
        let ten = Bls12Scalar::from_canonical_u64(10);
        let product = mul_naive(ten, g);
        let sum = g + g + g + g + g + g + g + g + g + g;
        assert_eq!(product, sum);
    }

    #[test]
    fn test_g1_multiplication() {
        let lhs = Bls12Scalar::from_canonical([11111111, 22222222, 33333333, 44444444]);
        assert_eq!(lhs * G1_GENERATOR_PROJECTIVE, mul_naive(lhs, G1_GENERATOR_PROJECTIVE));
    }

    /// A simple, somewhat inefficient implementation of multiplication which is used as a reference
    /// for correctness.
    fn mul_naive(lhs: Bls12Scalar, rhs: G1ProjectivePoint) -> G1ProjectivePoint {
        let mut g = rhs;
        let mut sum = G1ProjectivePoint::ZERO;
        for limb in lhs.to_canonical().iter() {
            for j in 0..64 {
                if (limb >> j & 1u64) != 0u64 {
                    sum = sum + g;
                }
                g = g.double();
            }
        }
        sum
    }
}
