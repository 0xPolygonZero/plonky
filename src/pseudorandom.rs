use crate::Field;

pub trait PRG<F: Field> {
    fn next_bool(&mut self) -> bool;

    fn next_u32(&mut self) -> u32;

    fn next_field(&mut self) -> F;
}

pub struct PRFBasedPRG<F: Field, Prf: PRF<F>> {
    prf: Prf,
    state: F,
}

impl<F: Field, Prf: PRF<F>> PRFBasedPRG<F, Prf> {
    pub fn seeded(prf: Prf, seed: F) -> Self {
        PRFBasedPRG { prf, state: seed }
    }

    pub fn unseeded(prf: Prf) -> Self {
        Self::seeded(prf, F::ZERO)
    }
}

impl<F: Field, Prf: PRF<F>> PRG<F> for PRFBasedPRG<F, Prf> {
    fn next_bool(&mut self) -> bool {
        self.next_u32() & 1 != 0
    }

    fn next_u32(&mut self) -> u32 {
        self.next_field().to_canonical_u32_vec()[0]
    }

    fn next_field(&mut self) -> F {
        let f = self.prf.rand(self.state);
        self.state = f;
        f
    }
}

pub trait PRF<F: Field> {
    fn rand(&self, x: F) -> F;
}
