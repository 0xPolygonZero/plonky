use crate::custom_gates::tree_builder::gates_to_tree;
use crate::{Gate, HaloCurve};
use std::collections::HashMap;
use std::sync::Arc;

pub type GatePrefixes = HashMap<String, Vec<bool>>;

// impl GatePrefixes {
//     pub fn prefix<G: Gate<C>>(&self, gate: &G) -> Vec<bool> {
//         self.get(gate.name())
//             .expect(&format!("Gate {} not found", gate.name()))
//             .clone()
//     }
//     pub fn prefix_from_str(&self, gate_name: &str) -> Vec<bool> {
//         self.get(gate_name)
//             .expect(&format!("Gate {} not found", gate_name))
//             .clone()
//     }
// }

#[derive(Clone)]
pub struct GateCollection<C: HaloCurve> {
    pub gates: Vec<Arc<dyn Gate<C>>>,
    pub prefixes: GatePrefixes,
}

impl<C: HaloCurve> From<Vec<Arc<dyn Gate<C>>>> for GateCollection<C> {
    fn from(gates: Vec<Arc<dyn Gate<C>>>) -> Self {
        let tree = gates_to_tree(
            &gates
                .iter()
                .map(|g| (g.name(), g.num_constants(), g.degree()))
                .collect::<Vec<_>>(),
        );
        let prefixes = tree.prefixes();
        prefixes
            .iter()
            .for_each(|(s, p)| trace!("Prefix for gate {}: {:?}", s, p));
        GateCollection { gates, prefixes }
    }
}

impl<C: HaloCurve> GateCollection<C> {
    pub fn prefix<G: Gate<C>>(&self, gate: &G) -> Vec<bool> {
        self.prefixes
            .get(gate.name())
            .unwrap_or_else(|| panic!("Gate {} not found.", gate.name()))
            .clone()
    }
    pub fn prefix_from_str(&self, gate_name: &str) -> Vec<bool> {
        self.prefixes
            .get(gate_name)
            .unwrap_or_else(|| panic!("Gate {} not found.", gate_name))
            .clone()
    }
}
