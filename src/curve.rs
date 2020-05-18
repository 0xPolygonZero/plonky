use std::ops::Neg;

use crate::Field;

/// A short Weierstrass curve.
pub trait Curve: 'static + Sized + Copy {
    type BaseField: Field;
    type ScalarField: Field;

    const A: Self::BaseField;
    const B: Self::BaseField;

    const GENERATOR_AFFINE: AffinePoint<Self>;

    const GENERATOR_PROJECTIVE: ProjectivePoint<Self> = ProjectivePoint {
        x: Self::GENERATOR_AFFINE.x,
        y: Self::GENERATOR_AFFINE.y,
        z: Self::BaseField::ONE,
        zero: false,
    };
}

/// A curve with the endomorphism described in the Halo paper, i.e. `phi((x, y)) = (zeta_p x, y)`,
/// where `phi(P) = [zeta_q] P` for some `zeta_q` of multiplicative order 3.
pub trait HaloEndomorphismCurve: Curve {
    const ZETA: Self::BaseField;
    const ZETA_SCALAR: Self::ScalarField;
}

/// A point on a short Weierstrass curve, represented in affine coordinates.
#[derive(Copy, Clone, Debug)]
pub struct AffinePoint<C: Curve> {
    pub x: C::BaseField,
    pub y: C::BaseField,
    pub zero: bool,
}

impl<C: Curve> AffinePoint<C> {
    pub const ZERO: Self = Self {
        x: C::BaseField::ZERO,
        y: C::BaseField::ZERO,
        zero: true,
    };

    pub fn nonzero(x: C::BaseField, y: C::BaseField) -> Self {
        let point = Self { x, y, zero: false };
        debug_assert!(point.is_valid());
        point
    }

    pub fn is_valid(&self) -> bool {
        let Self { x, y, zero } = *self;
        zero || y.square() == x.cube() + C::A * x + C::B
    }

    pub fn to_projective(&self) -> ProjectivePoint<C> {
        let Self { x, y, zero } = *self;
        ProjectivePoint {
            x,
            y,
            z: C::BaseField::ONE,
            zero,
        }
    }

    pub fn double(&self) -> Self {
        // TODO: This is a very lazy implementation...
        self.to_projective().double().to_affine()
    }
}

impl<C: Curve> PartialEq for AffinePoint<C> {
    fn eq(&self, other: &Self) -> bool {
        let AffinePoint {
            x: x1,
            y: y1,
            zero: zero1,
        } = *self;
        let AffinePoint {
            x: x2,
            y: y2,
            zero: zero2,
        } = *other;
        if zero1 || zero2 {
            return zero1 == zero2;
        }
        x1 == x2 && y1 == y2
    }
}

impl<C: Curve> Eq for AffinePoint<C> {}

/// A point on a short Weierstrass curve, represented in projective coordinates.
#[derive(Copy, Clone, Debug)]
pub struct ProjectivePoint<C: Curve> {
    pub x: C::BaseField,
    pub y: C::BaseField,
    pub z: C::BaseField,
    pub zero: bool,
}

impl<C: Curve> ProjectivePoint<C> {
    pub const ZERO: Self = Self {
        x: C::BaseField::ZERO,
        y: C::BaseField::ZERO,
        z: C::BaseField::ZERO,
        zero: true,
    };

    pub fn nonzero(x: C::BaseField, y: C::BaseField, z: C::BaseField) -> Self {
        let point = Self {
            x,
            y,
            z,
            zero: false,
        };
        debug_assert!(point.is_valid());
        point
    }

    pub fn is_valid(&self) -> bool {
        self.to_affine().is_valid()
    }

    pub fn to_affine(&self) -> AffinePoint<C> {
        let Self { x, y, z, zero } = *self;
        if zero {
            AffinePoint::ZERO
        } else {
            let z_inv = z.multiplicative_inverse_assuming_nonzero();
            AffinePoint::nonzero(x * z_inv, y * z_inv)
        }
    }

    pub fn batch_to_affine(proj_points: &[Self]) -> Vec<AffinePoint<C>> {
        let n = proj_points.len();
        let zs: Vec<C::BaseField> = proj_points.iter().map(|pp| pp.z).collect();
        let z_inv_opts = C::BaseField::batch_multiplicative_inverse_opt(&zs);

        let mut result = Vec::with_capacity(n);
        for i in 0..n {
            let Self { x, y, z, zero } = proj_points[i];
            result.push(if zero {
                AffinePoint::ZERO
            } else {
                let z_inv = z_inv_opts[i].unwrap();
                AffinePoint::nonzero(x * z_inv, y * z_inv)
            });
        }
        result
    }

    pub fn double(&self) -> Self {
        let Self { x, y, z, zero } = *self;
        if zero {
            return ProjectivePoint::ZERO;
        }

        let xx = x.square();
        let zz = z.square();
        let mut w = xx.triple();
        if C::A.is_nonzero() {
            w = w + C::A * zz;
        }
        let s = y.double() * z;
        let r = y * s;
        let rr = r.square();
        let b = (x + r).square() - (xx + rr);
        let h = w.square() - b.double();
        let x3 = h * s;
        let y3 = w * (b - h) - rr.double();
        let z3 = s.cube();
        Self {
            x: x3,
            y: y3,
            z: z3,
            zero: false,
        }
    }
}

impl<C: Curve> PartialEq for ProjectivePoint<C> {
    fn eq(&self, other: &Self) -> bool {
        let ProjectivePoint {
            x: x1,
            y: y1,
            z: z1,
            zero: zero1,
        } = *self;
        let ProjectivePoint {
            x: x2,
            y: y2,
            z: z2,
            zero: zero2,
        } = *other;
        if zero1 || zero2 {
            return zero1 == zero2;
        }

        // We want to compare (x1/z1, y1/z1) == (x2/z2, y2/z2).
        // But to avoid field division, it is better to compare (x1*z2, y1*z2) == (x2*z1, y2*z1).
        x1 * z2 == x2 * z1 && y1 * z2 == y2 * z1
    }
}

impl<C: Curve> Eq for ProjectivePoint<C> {}

impl<C: Curve> Neg for AffinePoint<C> {
    type Output = AffinePoint<C>;

    fn neg(self) -> Self::Output {
        let AffinePoint { x, y, zero } = self;
        AffinePoint { x, y: -y, zero }
    }
}

impl<C: Curve> Neg for ProjectivePoint<C> {
    type Output = ProjectivePoint<C>;

    fn neg(self) -> Self::Output {
        let ProjectivePoint { x, y, z, zero } = self;
        ProjectivePoint { x, y: -y, z, zero }
    }
}
