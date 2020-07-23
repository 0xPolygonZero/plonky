use plonky::plonk_util::eval_poly;
use plonky::Field;

pub struct Opening<F: Field> {
    pub local: F,
    pub right: F,
}

impl<F: Field> Opening<F> {
    fn to_vec(&self) -> Vec<F> {
        vec![self.local, self.right]
    }
}

impl<F: Field> From<(F, F)> for Opening<F> {
    fn from(t: (F, F)) -> Self {
        Self {
            local: t.0,
            right: t.1,
        }
    }
}

/// Struct containing all polynomial openings used in the Plookup protocol.
pub struct PlookupOpenings<F: Field> {
    pub f: Opening<F>,
    pub t: Opening<F>,
    pub h1: Opening<F>,
    pub h2: Opening<F>,
    pub z: Opening<F>,
    pub quotient: Opening<F>,
}

impl<F: Field> PlookupOpenings<F> {
    pub fn to_vec(&self) -> Vec<F> {
        [
            self.f.to_vec(),
            self.t.to_vec(),
            self.h1.to_vec(),
            self.h2.to_vec(),
            self.z.to_vec(),
            self.quotient.to_vec(),
        ]
        .concat()
    }

    pub fn local(&self) -> Vec<F> {
        vec![
            self.f.local,
            self.t.local,
            self.h1.local,
            self.h2.local,
            self.z.local,
            self.quotient.local,
        ]
    }

    pub fn right(&self) -> Vec<F> {
        vec![
            self.f.right,
            self.t.right,
            self.h1.right,
            self.h2.right,
            self.z.right,
            self.quotient.right,
        ]
    }
}

pub fn open_all_polynomials<F: Field>(
    f_coeffs: &[F],
    t_coeffs: &[F],
    h1_coeffs: &[F],
    h2_coeffs: &[F],
    z_coeffs: &[F],
    quotient_coeffs: &[F],
    zeta: F,
    generator: F,
) -> PlookupOpenings<F> {
    let right = zeta * generator;
    let f = (eval_poly(f_coeffs, zeta), eval_poly(f_coeffs, right)).into();
    let t = (eval_poly(t_coeffs, zeta), eval_poly(t_coeffs, right)).into();
    let h1 = (eval_poly(h1_coeffs, zeta), eval_poly(h1_coeffs, right)).into();
    let h2 = (eval_poly(h2_coeffs, zeta), eval_poly(h2_coeffs, right)).into();
    let z = (eval_poly(z_coeffs, zeta), eval_poly(z_coeffs, right)).into();
    let quotient = (
        eval_poly(quotient_coeffs, zeta),
        eval_poly(quotient_coeffs, right),
    )
        .into();
    PlookupOpenings {
        f,
        t,
        h1,
        h2,
        z,
        quotient,
    }
}
