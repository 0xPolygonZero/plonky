use std::collections::HashMap;

use crate::{Field, Target2};

#[derive(Debug)]
pub(crate) struct Partitions2<F: Field> {
    partitions: Vec<Vec<Target2<F>>>,
    indices: HashMap<Target2<F>, usize>,
}

impl<F: Field> Partitions2<F> {
    pub fn new() -> Self {
        Self {
            partitions: Vec::new(),
            indices: HashMap::new(),
        }
    }

    /// Adds the targets as new singleton partitions if they are not already present, then merges
    /// their partitions if the targets are not already in the same partition.
    pub fn merge(&mut self, a: Target2<F>, b: Target2<F>) {
        let a_index = self.get_index(a);
        let b_index = self.get_index(b);

        if a_index != b_index {
            // Merge a's partition into b's partition, leaving a's partition empty.
            // We have to clone because Rust's borrow checker doesn't know that
            // self.partitions[b_index] and self.partitions[b_index] are disjoint.
            let mut a_partition = self.partitions[a_index].clone();
            let b_partition = &mut self.partitions[b_index];
            for a_sibling in &a_partition {
                *self.indices.get_mut(a_sibling).unwrap() = b_index;
            }
            b_partition.append(&mut a_partition);
        }
    }

    /// Gets the partition index of a given target. If the target is not present, adds it as a new
    /// singleton partition and returns the new partition's index.
    fn get_index(&mut self, target: Target2<F>) -> usize {
        if let Some(&index) = self.indices.get(&target) {
            index
        } else {
            let index = self.partitions.len();
            self.partitions.push(vec![target]);
            index
        }
    }
}
