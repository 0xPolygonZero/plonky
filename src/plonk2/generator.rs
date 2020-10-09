use crate::{Field, Target2, PartialWitness2};

/// A generator participates in the generation of the witness.
pub trait WitnessGenerator2<F: Field> {
    /// Targets to be "watched" by this generator. Whenever a target in the watch list is populated,
    /// the generator will be queued to run.
    fn watch_list(&self) -> Vec<Target2<F>>;

    /// Run this generator, returning a `PartialWitness` containing any new witness elements, and a
    /// flag indicating whether the generator is finished. If the flag is true, the generator will
    /// never be run again, otherwise it will be queued for another run next time a target in its
    /// watch list is populated.
    fn run(&mut self, witness: &PartialWitness2<F>) -> (PartialWitness2<F>, bool);
}

/// A generator which runs once after a list of dependencies is present in the witness.
pub trait SimpleGenerator<F: Field> {
    fn dependencies(&self) -> Vec<Target2<F>>;

    fn run_once(&mut self, witness: &PartialWitness2<F>) -> PartialWitness2<F>;
}

impl<F: Field> WitnessGenerator2<F> for dyn SimpleGenerator<F> {
    fn watch_list(&self) -> Vec<Target2<F>> {
        self.dependencies()
    }

    fn run(&mut self, witness: &PartialWitness2<F>) -> (PartialWitness2<F>, bool) {
        if witness.contains_all(&self.dependencies()) {
            (self.run_once(witness), true)
        } else {
            (PartialWitness2::new(), false)
        }
    }
}
