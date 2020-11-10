use std::collections::HashMap;

use crate::{Field, Target2};

#[derive(Debug)]
pub struct PartialWitness2<F: Field> {
    wire_values: HashMap<Target2<F>, F>,
}

impl<F: Field> PartialWitness2<F> {
    pub fn new() -> Self {
        PartialWitness2 {
            wire_values: HashMap::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.wire_values.is_empty()
    }

    pub fn try_get(&self, target: Target2<F>) -> Option<F> {
        self.wire_values.get(&target).cloned()
    }

    pub fn contains(&self, target: Target2<F>) -> bool {
        self.wire_values.contains_key(&target)
    }

    pub fn contains_all(&self, targets: &[Target2<F>]) -> bool {
        targets.iter().all(|&t| self.contains(t))
    }

    pub fn set(&mut self, target: Target2<F>, value: F) {
        self.wire_values.insert(target, value);
    }
}
