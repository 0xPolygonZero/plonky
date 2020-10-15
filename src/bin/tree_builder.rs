use plonky::custom_gates::tree_builder::gates_to_tree;
use plonky::util::get_canonical_gates_containers;
use plonky::{Tweedledee, Tweedledum};

fn main() {
    type C = Tweedledum;
    type InnerC = Tweedledee;

    let gates = get_canonical_gates_containers::<C, InnerC>()
        .into_iter()
        .map(|g| (g.name(), g.num_constants(), g.degree()))
        .collect::<Vec<_>>();
    let tree = gates_to_tree(&gates);
    // tree.graph();
    dbg!(tree.prefixes());
}
