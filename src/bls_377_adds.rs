use std::ops::Add;

use crate::{G1AffinePoint, G1ProjectivePoint, Field};

impl Add<G1ProjectivePoint> for G1ProjectivePoint {
    type Output = G1ProjectivePoint;

    fn add(self, rhs: G1ProjectivePoint) -> Self::Output {
        if self.is_zero() {
            return rhs;
        }
        if rhs.is_zero() {
            return self;
        }

        let G1ProjectivePoint { x: x1, y: y1, z: z1 } = self;
        let G1ProjectivePoint { x: x2, y: y2, z: z2 } = rhs;

        let x1z2 = x1 * z2;
        let y1z2 = y1 * z2;
        let x2z1 = x2 * z1;
        let y2z1 = y2 * z1;

        // Check if we're doubling or adding inverses.
        if x1z2 == x2z1 {
            if y1z2 == y2z1 {
                // TODO: inline to avoid redundant muls.
                return self.double();
            }
            if y1z2 == -y2z1 {
                return G1ProjectivePoint::ZERO;
            }
        }

        let z1z2 = z1 * z2;
        let u = y2z1 - y1z2;
        let uu = u.square();
        let v = x2z1 - x1z2;
        let vv = v.square();
        let vvv = v * vv;
        let r = vv * x1z2;
        let a = uu * z1z2 - vvv - r.double();
        let x3 = v * a;
        let y3 = u * (r - a) - vvv * y1z2;
        let z3 = vvv * z1z2;
        G1ProjectivePoint { x: x3, y: y3, z: z3 }
    }
}

impl Add<G1AffinePoint> for G1ProjectivePoint {
    type Output = G1ProjectivePoint;

    fn add(self, rhs: G1AffinePoint) -> Self::Output {
        if self.is_zero() {
            return rhs.to_projective();
        }
        if rhs.is_zero() {
            return self;
        }

        let G1ProjectivePoint { x: x1, y: y1, z: z1 } = self;
        let G1AffinePoint { x: x2, y: y2 } = rhs;

        let x2z1 = x2 * z1;
        let y2z1 = y2 * z1;

        // Check if we're doubling or adding inverses.
        if x1 == x2z1 {
            if y1 == y2z1 {
                // TODO: inline to avoid redundant muls.
                return self.double();
            }
            if y1 == -y2z1 {
                return G1ProjectivePoint::ZERO;
            }
        }

        let u = y2z1 - y1;
        let uu = u.square();
        let v = x2z1 - x1;
        let vv = v.square();
        let vvv = v * vv;
        let r = vv * x1;
        let a = uu * z1 - vvv - r.double();
        let x3 = v * a;
        let y3 = u * (r - a) - vvv * y1;
        let z3 = vvv * z1;
        G1ProjectivePoint { x: x3, y: y3, z: z3 }
    }
}

impl Add<G1AffinePoint> for G1AffinePoint {
    type Output = G1ProjectivePoint;

    fn add(self, rhs: G1AffinePoint) -> Self::Output {
        if self.is_zero() {
            return rhs.to_projective();
        }
        if rhs.is_zero() {
            return self.to_projective();
        }

        let G1AffinePoint { x: x1, y: y1 } = self;
        let G1AffinePoint { x: x2, y: y2 } = rhs;

        // Check if we're doubling or adding inverses.
        if x1 == x2 {
            if y1 == y2 {
                return self.to_projective().double();
            }
            if y1 == -y2 {
                return G1ProjectivePoint::ZERO;
            }
        }

        let u = y2 - y1;
        let uu = u.square();
        let v = x2 - x1;
        let vv = v.square();
        let vvv = v * vv;
        let r = vv * x1;
        let a = uu - vvv - r.double();
        let x3 = v * a;
        let y3 = u * (r - a) - vvv * y1;
        let z3 = vvv;
        G1ProjectivePoint { x: x3, y: y3, z: z3 }
    }
}
