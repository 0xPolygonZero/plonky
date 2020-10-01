use plonky::{CurveAddGate, RescueStepBGate, RescueStepAGate, ArithmeticGate, ConstantGate, BufferGate, PublicInputGate, Base4SumGate, CurveEndoGate, CurveDblGate, Tweedledum, Tweedledee, Gate };
#[macro_use]
extern crate log;

pub(crate) const NUM_CONSTANTS: usize = 6;
pub(crate) const QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER: usize = 7;


#[derive(Debug, Clone)]
struct Tree {
    node: String,
    left: Option<Box<Tree>>,
    right: Option<Box<Tree>>,
}

fn label(s: &str, i: &mut usize) -> String {
    if s.is_empty() {
        *i += 1;
        format!("X{}", i)
    } else {
        s.to_string()
    }
}

impl Tree {
    fn new(node: String, left: Option<Tree>, right: Option<Tree>) -> Self {
        Tree {
            node,
            left: left.map(Box::new),
            right: right.map(Box::new),

        }
    }

    /// Saves the tree to a `graph.dot` in Graphviz format.
    fn graph(&self) {
        let mut i = 0;
        let mut res = vec![];
        let mut first = self.clone();
        first.node = label(&self.node, &mut i);
        let mut to_show = vec![first];
        while !to_show.is_empty() {
            let n = to_show.pop().unwrap();
            if let Some(l) = &n.left {
                let mut ll = *l.clone();
                ll.node = label(&l.node, &mut i);
                res.push(format!("{} -> {}", &n.node, &ll.node));
                to_show.insert(0, ll);
            }
            if let Some(r) = &n.right {
                let mut rr = *r.clone();
                rr.node = label(&r.node, &mut i);
                res.push(format!("{} -> {}", &n.node, &rr.node));
                to_show.insert(0, rr);
            }
        }
        let mut fin = (1..=i).map(|j| {
            format!("X{} [label=\"\"]", j)
        }).collect::<Vec<_>>();
        fin.extend(res);
        let fin: String = fin.join("\n");
        let graph = format!("digraph G {{\ngraph [ordering=\"out\"]\n{}\n}}", fin);
        std::fs::write("graph.dot", graph).expect("Cannot write to file graph.dot");
    }
}

type GateDesc = (String, usize, usize, usize);

fn construct_tree(gates: &[GateDesc], depth: usize, level: usize) -> Vec<(Tree, Vec<GateDesc>)> {
    if gates.is_empty() || gates[0].1 - 1 < level {
        return vec![];
    } else if depth == 0 {
        return vec![(Tree::new(gates[0].0.clone(), None, None), gates[1..].to_vec())];
    }

    let mut left = construct_tree(gates, depth - 1, level + 1);
    left.push((Tree::new(gates[0].0.clone(), None, None), gates[1..].to_vec()));
    let mut ans = vec![];
    for (l, gs) in left.into_iter() {
        if gs.is_empty() {
            ans.push((Tree::new(String::new(), Some(l.clone()), None), gs.clone()));
        }
        let mut right = construct_tree(&gs, depth - 1, level + 1);
        if !gs.is_empty() {
            right.push((Tree::new(gs[0].0.clone(), None, None), gs[1..].to_vec()));
        }
        for (r, gs) in right.into_iter() {
            ans.push((Tree::new(String::new(), Some(l.clone()), Some(r)), gs));
        }
    }
    ans
}

fn gates_to_tree(gates: &[(&str, usize, usize)]) -> Tree {
    let mut gates = gates.iter().map(|x| {
        let max_prefix_len = (NUM_CONSTANTS - x.1).min(QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER + 1 - x.2);
        (x.0.to_string(), max_prefix_len, x.1, x.2)
    }).collect::<Vec<_>>();
    gates.sort_by(|x, y| {
        x.1.cmp(&y.1)
    });

    for d in 1..100 {
        let trees = construct_tree(&gates, d, 0);
        let tree = trees.iter().find(|x| x.1.is_empty());
        if let Some(t) = tree {
            info!("Tree of depth {} found", d);
            return t.0.clone();
        }
    }
    panic!("No tree found.")
}

fn main() {
    pretty_env_logger::init();
    type C = Tweedledum;
    type InnerC = Tweedledee;

    let gates = [
        (CurveAddGate::<C, InnerC>::NAME, CurveAddGate::<C, InnerC>::NUM_CONSTANTS, CurveAddGate::<C, InnerC>::DEGREE),
        (CurveDblGate::<C, InnerC>::NAME, CurveDblGate::<C, InnerC>::NUM_CONSTANTS, CurveDblGate::<C, InnerC>::DEGREE),
        (CurveEndoGate::<C, InnerC>::NAME, CurveEndoGate::<C, InnerC>::NUM_CONSTANTS, CurveEndoGate::<C, InnerC>::DEGREE),
        (Base4SumGate::<C>::NAME, Base4SumGate::<C>::NUM_CONSTANTS, Base4SumGate::<C>::DEGREE),
        (PublicInputGate::<C>::NAME, PublicInputGate::<C>::NUM_CONSTANTS, PublicInputGate::<C>::DEGREE),
        (BufferGate::<C>::NAME, BufferGate::<C>::NUM_CONSTANTS, BufferGate::<C>::DEGREE),
        (ConstantGate::<C>::NAME, ConstantGate::<C>::NUM_CONSTANTS, ConstantGate::<C>::DEGREE),
        (ArithmeticGate::<C>::NAME, ArithmeticGate::<C>::NUM_CONSTANTS, ArithmeticGate::<C>::DEGREE),
        (RescueStepAGate::<C>::NAME, RescueStepAGate::<C>::NUM_CONSTANTS, RescueStepAGate::<C>::DEGREE),
        (RescueStepBGate::<C>::NAME, RescueStepBGate::<C>::NUM_CONSTANTS, RescueStepBGate::<C>::DEGREE),
    ];

    let tree = gates_to_tree(&gates);
    tree.graph();
}