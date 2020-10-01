use plonky::{CurveAddGate, RescueStepBGate, RescueStepAGate, ArithmeticGate, ConstantGate, BufferGate, PublicInputGate, Base4SumGate, CurveEndoGate, CurveDblGate, Tweedledum, Tweedledee, Gate, Curve, fft_with_precomputation_power_of_2};

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
    fn graph(&self) {
        let mut i = 0;
        let mut res = vec![];
        let mut first = self.clone();
        first.node = label(&self.node, &mut i);
        let mut to_show = vec![self.clone()];
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
        (1..=i).for_each(|j| {
            println!("X{} [label=\"\"]", j);
        });
        println!("{}", res.join("\n"));
    }

    fn show(&self) {
        let mut levels = vec![vec![Some(self.clone())]];
        loop {
            if levels.last().unwrap().iter().all(|x| x.is_none()) {
                levels.pop();
                break;
            }
            let mut new = vec![];
            for n in levels.last().unwrap() {
                if let Some(nn) = n {
                    if let Some(l) = &nn.left {
                        new.push(Some(*l.clone()));
                    } else {
                        new.push(None);
                    }
                    if let Some(r) = &nn.right {
                        new.push(Some(*r.clone()));
                    } else {
                        new.push(None);
                    }
                } else {
                    new.push(None);
                    new.push(None);
                }
            }
            levels.push(new);
        }
        let d = levels.len();
        for (i, l) in levels.into_iter().enumerate() {
            let sep1 = if i == d - 1 { "".to_string() } else {
                "\t\t\t".repeat(2usize.pow((d - 2 - i) as u32))
            };
            let sep2 = "\t\t\t".repeat(2usize.pow((d - 1 - i) as u32));
            println!(
                "{}{}", sep1,
                l.into_iter().map(|o| {
                    if let Some(n) = o {
                        if n.node.is_empty() {
                            "O".to_string()
                        } else {
                            n.node
                        }
                    } else {
                        "X".to_string()
                    }
                }).collect::<Vec<_>>().join(&sep2)
            );
        }

        // let mut to_show = vec![(self, 0, 'o')];
        // let mut cur = 0;
        // while !to_show.is_empty() {
        //     let n = to_show.pop();
        //     if let Some(t) = n {
        //         if t.1 > cur {
        //             println!();
        //             cur = t.1;
        //         }
        //         if t.0.node.is_empty() {
        //             print!("O\t");
        //         } else {
        //             print!("{}{}\t",t.2, t.0.node);
        //         }
        //         if let Some(l) = &t.0.left {
        //             to_show.insert(0, (l, t.1+1, 'l'));
        //         }
        //         if let Some(r) = &t.0.right {
        //             to_show.insert(0, (r, t.1+1, 'r'));
        //         }
        //     }
        // }
    }
}

type GateDesc = (String, usize, usize, usize);

fn construct_tree(gates: &[GateDesc], depth: usize, level: usize) -> Vec<(Tree, Vec<GateDesc>)> {
    // dbg!(depth, gates);
    if gates.is_empty() || gates[0].1-1 < level {
        return vec![];
    } else if depth == 0 {
        return vec![(Tree::new(gates[0].0.clone(), None, None), gates[1..].to_vec())];
    }

    // let decremented_gates = gates.iter().map(|g| (g.0.clone(), g.1 - 1, g.2, g.3)).collect::<Vec<_>>();
    let mut left = construct_tree(gates, depth - 1, level + 1);
    left.push((Tree::new(gates[0].0.clone(), None, None), gates[1..].to_vec()));
    let mut ans = vec![];
    for (l, gs) in left.into_iter() {
        if gs.is_empty() {
            ans.push((Tree::new(String::new(), Some(l.clone()), None), gs.clone()));
        }
        let mut right = construct_tree(&gs, depth - 1, level+1);
        if !gs.is_empty() {
            right.push((Tree::new(gs[0].0.clone(), None, None), gs[1..].to_vec()));
        }
        for (r, gs) in right.into_iter() {
// if gs.is_empty() {
//     ans.push((Tree::new(String::new(), Some(l.clone()), Some(r)), gs));
// }
            ans.push((Tree::new(String::new(), Some(l.clone()), Some(r)), gs));
        }
    }
    ans
}

fn main() {
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

    let bam = gates;
    let mut boom = bam.iter().map(|x| {
        let max_prefix_len = (NUM_CONSTANTS - x.1).min(QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER + 1 - x.2);
        (x.0.to_string(), max_prefix_len, x.1, x.2)
    }).collect::<Vec<_>>();
    boom.sort_by(|x, y| {
        x.1.cmp(&y.1)
    });
// dbg!(gates);
    dbg!(&boom);

// dbg!(construct_tree(&boom[..3].to_vec(), 2));
    let lol = construct_tree(&boom[..].to_vec(), 6, 0);
    dbg!(lol.iter().map(|x| x.1.len()).collect::<Vec<_>>());
    let tt = lol.iter().find(|x| x.1.is_empty());
    if let Some(t) = tt {
        dbg!(&t.0);
        dbg!(&t.1);
        t.0.show();
        t.0.graph();
    }
}