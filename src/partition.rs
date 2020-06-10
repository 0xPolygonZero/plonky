use crate::{Field, Target, Wire, NUM_ROUTED_WIRES, NUM_WIRES};
use rand_chacha::rand_core::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct TargetPartitions {
    partitions: Vec<Vec<Target>>,
    indices: HashMap<Target, usize>,
}

impl TargetPartitions {
    pub fn new() -> Self {
        Self {
            partitions: Vec::new(),
            indices: HashMap::new(),
        }
    }

    pub fn get_partition(&self, target: Target) -> &[Target] {
        &self.partitions[self.indices[&target]]
    }

    /// Add a new partition with a single member.
    pub fn add_partition(&mut self, target: Target) {
        let index = self.partitions.len();
        self.partitions.push(vec![target]);
        self.indices.insert(target, index);
    }

    /// Merge the two partitions containing the two given targets. Does nothing if the targets are
    /// already members of the same partition.
    pub fn merge(&mut self, a: Target, b: Target) {
        let a_index = self.indices[&a];
        let b_index = self.indices[&b];
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

    pub fn to_wire_partitions(&self) -> WirePartitions {
        // Here we just drop all CircuitInputs, leaving all GateInputs.
        let mut partitions = Vec::new();
        let mut indices = HashMap::new();

        for old_partition in &self.partitions {
            let mut new_partition = Vec::new();
            for target in old_partition {
                if let Target::Wire(w) = *target {
                    new_partition.push(w);
                }
            }
            partitions.push(new_partition);
        }

        for (&target, &index) in &self.indices {
            if let Target::Wire(gi) = target {
                indices.insert(gi, index);
            }
        }

        let result = WirePartitions {
            partitions,
            indices,
        };
        result.assert_valid();
        result
    }
}

pub struct WirePartitions {
    partitions: Vec<Vec<Wire>>,
    indices: HashMap<Wire, usize>,
}

impl WirePartitions {
    fn assert_valid(&self) {
        for partition in &self.partitions {
            for wire in partition {
                if wire.input >= NUM_ROUTED_WIRES {
                    assert_eq!(
                        partition.len(),
                        1,
                        "Non-routed wires should not be in a partition containing other wires"
                    );
                }
            }
        }
    }

    /// Find a wire's "neighbor" in the context of Plonk's "extended copy constraints" check. In
    /// other words, find the next wire in the given wire's partition. If the given wire is last in
    /// its partition, this will loop around. If the given wire has a partition all to itself, it
    /// is considered its own neighbor.
    fn get_neighbor(&self, wire: Wire) -> Wire {
        let partition = &self.partitions[self.indices[&wire]];
        let n = partition.len();
        for i in 0..n {
            if partition[i] == wire {
                let neighbor_index = (i + 1) % n;
                return partition[neighbor_index];
            }
        }
        panic!("Wire not found in the expected partition")
    }

    /// Generates sigma in the context of Plonk, which is a map from `[kn]` to `[kn]`, where `k` is
    /// the number of routed wires and `n` is the number of gates.
    pub fn to_sigma(&self) -> Vec<usize> {
        debug_assert_eq!(self.indices.len() % NUM_WIRES, 0);
        let num_all_wires = self.indices.len();
        let num_gates = num_all_wires / NUM_WIRES;

        let mut sigma = Vec::new();
        for input in 0..NUM_ROUTED_WIRES {
            for gate in 0..num_gates {
                let wire = Wire { gate, input };
                let neighbor = self.get_neighbor(wire);
                sigma.push(neighbor.input * num_gates + neighbor.gate);
            }
        }
        sigma
    }
}

/// Returns `k_i`, the multiplier used in `S_ID_i` in the context of Plonk's permutation argument.
pub(crate) fn get_subgroup_shift<F: Field>(i: usize) -> F {
    // The optimized variant of Plonk's permutation argument calls for NUM_ROUTED_WIRES shifts,
    // k_1, ..., k_n, which result in distinct cosets. The paper suggests a method which is
    // fairly straightforward when only three shifts are needed, but seems a bit complex and
    // expensive if more are needed.

    // We will "cheat" and just use random field elements. Since our subgroup has |F*|/degree
    // possible cosets, the probability of a collision is negligible for large fields.

    if i == 0 {
        // We use a trivial shift of 1 for k_1, as in the paper, to save a multiplication.
        F::ONE
    } else {
        let mut rng = ChaCha8Rng::seed_from_u64(i as u64);
        F::rand_from_rng(&mut rng)
    }
}
