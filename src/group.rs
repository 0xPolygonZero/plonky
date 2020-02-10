use std::ops::{Add, Mul};

use chashmap::CHashMap;

use lazy_static::lazy_static;

use crate::{Bls12Base, Bls12Scalar};

// Parameters taken from the implementation of Bls12-377 in Zexe found here:
// https://github.com/scipr-lab/zexe/blob/master/algebra/src/curves/bls12_377/g1.rs

const A: Bls12Base = Bls12Base { limbs: [0x0, 0x0, 0x0, 0x0, 0x0, 0x0] };

const B: Bls12Base = Bls12Base {
    limbs: [
        0x2cdffffffffff68,
        0x51409f837fffffb1,
        0x9f7db3a98a7d3ff2,
        0x7b4e97b76e7c6305,
        0x4cf495bf803c84e8,
        0x008d6661e2fdf49a,
    ]
};

const COFACTOR: &'static [u64] = &[0x0, 0x170b5d4430000000];

const COFACTOR_INV: Bls12Scalar = Bls12Scalar {
    limbs: [
        2013239619100046060,
        4201184776506987597,
        2526766393982337036,
        1114629510922847535,
    ]
};

// 81937999373150964239938255573465948239988671502647976594219695644855304257327692006745978603320413799295628339695
pub const G1_GENERATOR_X: Bls12Base = Bls12Base {
    limbs: [
        0x260f33b9772451f4,
        0xc54dd773169d5658,
        0x5c1551c469a510dd,
        0x761662e4425e1698,
        0xc97d78cc6f065272,
        0xa41206b361fd4d,
    ]
};

// 241266749859715473739788878240585681733927191168601896383759122102112907357779751001206799952863815012735208165030
pub const G1_GENERATOR_Y: Bls12Base = Bls12Base {
    limbs: [
        0x8193961fb8cb81f3,
        0x638d4c5f44adb8,
        0xfafaf3dad4daf54a,
        0xc27849e2d655cd18,
        0x2ec3ddb401d52814,
        0x7da93326303c71,
    ]
};

#[derive(Eq, PartialEq, Hash, Copy, Clone, Debug)]
pub struct G1ProjectivePoint {
    pub x: Bls12Base,
    pub y: Bls12Base,
    pub z: Bls12Base,
}

impl Add<G1ProjectivePoint> for G1ProjectivePoint {
    type Output = G1ProjectivePoint;

    /// Safe version of addition with non-zero checks
    /// From https://www.hyperelliptic.org/EFD/g1p/data/shortw/projective/addition/add-1998-cmo-2
    fn add(self, rhs: G1ProjectivePoint) -> Self::Output {
        if self.is_zero() {
            rhs
        } else if rhs.is_zero() {
            self
        } else if self.x == -rhs.x {
            G1ProjectivePoint::ZERO
        } else {
            let y1z2 = self.y * rhs.z;
            let x1z2 = self.x * rhs.z;
            let z1z2 = self.z * rhs.z;
            let u = rhs.y * self.z - y1z2;
            let uu = u * u;
            let v = rhs.x * self.z - x1z2;
            let vv = v * v;
            let vvv = v * vv;
            let r = vv * x1z2;
            let a = uu * z1z2 - vvv - (r + r);
            let x3 = v * a;
            let y3 = u * (r - a) - vvv * y1z2;
            let z3 = vvv * z1z2;
            G1ProjectivePoint { x: x3, y: y3, z: z3 }
        }
    }
}

const WINDOW_BITS: usize = 4;
const BASE: usize = 1 << WINDOW_BITS;
const DIGITS: usize = (Bls12Scalar::BITS + WINDOW_BITS - 1) / WINDOW_BITS;

/// Precomputed state used for Bls12Scalar x G1ProjectivePoint multiplications,
/// specific to a particular generator.
#[derive(Copy, Clone)]
struct G1GeneratorPrecomputations {
    /// [(2^w)^i] g for each i < DIGITS
    powers: [G1ProjectivePoint; DIGITS],
}

// TODO: Use compressed coordinates in the cache.
lazy_static! {
    static ref G1_MUL_PRECOMPUTATIONS: CHashMap<G1ProjectivePoint, G1GeneratorPrecomputations> = CHashMap::new();
}

fn get_precomputation(g: G1ProjectivePoint) -> G1GeneratorPrecomputations {
    match G1_MUL_PRECOMPUTATIONS.get(&g) {
        Some(x) => *x,
        None => {
            let precomputation = precompute(g);
            G1_MUL_PRECOMPUTATIONS.insert(g, precomputation);
            precomputation
        }
    }
}

fn precompute(g: G1ProjectivePoint) -> G1GeneratorPrecomputations {
    let mut powers = [G1ProjectivePoint::ZERO; DIGITS];
    powers[0] = g;
    for i in 1..DIGITS {
        powers[i] = powers[i - 1];
        for _j in 0..WINDOW_BITS {
            powers[i] = powers[i] + powers[i];
        }
    }
    G1GeneratorPrecomputations { powers }
}

impl Mul<G1ProjectivePoint> for Bls12Scalar {
    type Output = G1ProjectivePoint;

    fn mul(self, rhs: G1ProjectivePoint) -> Self::Output {
        // Yao's method; see https://koclab.cs.ucsb.edu/teaching/ecc/eccPapers/Doche-ch09.pdf
        let precomputed_powers = get_precomputation(rhs).powers;

        let digits = to_digits(&self);

        let mut y = G1ProjectivePoint::ZERO;
        let mut u = G1ProjectivePoint::ZERO;
        let mut ops = 0; // todo
        for j in (1..BASE).rev() {
            for (i, &digit) in digits.iter().enumerate() {
                if digit == j as u64 {
                    u = u + precomputed_powers[i];
                    ops += 1;
                }
            }
            y = y + u;
            ops += 1;
        }
//        println!("{} {} {}", ops, BASE, DIGITS);
        y
    }
}

fn to_digits(x: &Bls12Scalar) -> [u64; DIGITS] {
    debug_assert!(64 % WINDOW_BITS == 0,
                  "For simplicity, only power-of-two window sizes are handled for now");
    let digits_per_u64 = 64 / WINDOW_BITS;
    let mut digits = [0; DIGITS];
    for (i, limb) in x.limbs.iter().enumerate() {
        for j in 0..digits_per_u64 {
            digits[i*digits_per_u64 + j] = (limb >> (j * WINDOW_BITS)) % BASE as u64;
        }
    }

    // TODO
//    for digit in digits.iter() {
//        print!("{} ", digit);
//    }
//    println!();

    digits
}

impl G1ProjectivePoint {
    pub fn is_zero(&self) -> bool {
        self.z == Bls12Base::ZERO
    }

    const ZERO: G1ProjectivePoint = G1ProjectivePoint { x: Bls12Base::ZERO, y: Bls12Base::ZERO, z: Bls12Base::ZERO };

    /// Doubling of a G1 point
    /// From https://www.hyperelliptic.org/EFD/g1p/data/shortw/projective/doubling/dbl-2007-bl
    pub fn double(&self) -> G1ProjectivePoint {
        let w = self.x * self.x * 3u64;
        let s = self.y * self.z;
        let ss = s * s;
        let sss = s * ss;
        let r = self.y * s;
        let b = self.x * r;
        let h = w * w - b * 8u64;
        let x3 = h * s * 2u64;
        let y3 = w * (b * 4u64 - h) - r * r * 8u64;
        let z3 = sss * 8u64;
        G1ProjectivePoint { x: x3, y: y3, z: z3 }
    }

    pub fn new(x: Bls12Base, y: Bls12Base, z: Bls12Base) -> G1ProjectivePoint {
        assert!(G1ProjectivePoint::is_on_curve(x, y, z) /*&& is_in_subgroup(x, y, z)*/);
        G1ProjectivePoint { x: x, y: y, z: z }
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
}

#[cfg(test)]
mod tests {
    use crate::{Bls12Base, Bls12Scalar, G1ProjectivePoint};

    #[test]
    fn test_g1_multiplication() {
        let lhs = Bls12Scalar { limbs: [11111111, 22222222, 33333333, 44444444] };
        let s_part = 0b1010111010101010101000101010101010101011101010101010101010101010u64;
        let lhs = Bls12Scalar { limbs: [s_part, s_part, s_part, s_part] };
        // TODO: Not on the curve.
        let rhs = G1ProjectivePoint { x: Bls12Base::ZERO, y: Bls12Base::ONE, z: Bls12Base::ONE };
        assert_eq!(lhs * rhs, mul_naive(lhs, rhs));
    }

    /// A simple, somewhat inefficient implementation of multiplication which is used as a reference
    /// for correctness.
    fn mul_naive(lhs: Bls12Scalar, rhs: G1ProjectivePoint) -> G1ProjectivePoint {
        let mut g = rhs;
        let mut sum = G1ProjectivePoint::ZERO;
        for limb in lhs.limbs.iter() {
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
