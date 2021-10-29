use std::collections::HashMap;

use crate::{Field, Target2, Wire};

#[derive(Debug)]
pub struct PartialWitness2<F: Field> {
    target_values: HashMap<Target2<F>, F>,
}

impl<F: Field> PartialWitness2<F> {
    pub fn new() -> Self {
        PartialWitness2 {
            target_values: HashMap::new(),
        }
    }

    pub fn singleton(target: Target2<F>, value: F) -> Self {
        let mut witness = PartialWitness2::new();
        witness.set_target(target, value);
        witness
    }

    pub fn is_empty(&self) -> bool {
        self.target_values.is_empty()
    }

    pub fn get_target(&self, target: Target2<F>) -> F {
        self.target_values[&target]
    }

    pub fn try_get_target(&self, target: Target2<F>) -> Option<F> {
        self.target_values.get(&target).cloned()
    }

    pub fn get_wire(&self, wire: Wire) -> F {
        self.get_target(Target2::Wire(wire))
    }

    pub fn contains(&self, target: Target2<F>) -> bool {
        self.target_values.contains_key(&target)
    }

    pub fn contains_all(&self, targets: &[Target2<F>]) -> bool {
        targets.iter().all(|&t| self.contains(t))
    }

    pub fn set_target(&mut self, target: Target2<F>, value: F) {
        self.target_values.insert(target, value);
    }

    pub fn set_wire(&mut self, wire: Wire, value: F) {
        self.set_target(Target2::Wire(wire), value)
    }
}
