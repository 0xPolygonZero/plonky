use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::iter::{Product, Sum};
use std::ops::{Add, Mul, Neg, Sub};
use std::ptr;
use std::rc::Rc;

use crate::{Field, Wire};

pub(crate) struct EvaluationVars<'a, F: Field> {
    pub(crate) local_constants: &'a [F],
    pub(crate) next_constants: &'a [F],
    pub(crate) local_wire_values: &'a [F],
    pub(crate) next_wire_values: &'a [F],
}

/// A polynomial over all the variables that are subject to constraints (local constants, next
/// constants, local wire values, and next wire values). This representation does not require any
/// particular form; it permits arbitrary forms such as `(x + 1)^3 + y z`.
// Implementation note: This is a wrapper because we want to hide complexity behind
// `ConstraintPolynomialInner` and `ConstraintPolynomialRef`. In particular, the caller shouldn't
// need to know that we use reference counting internally, and shouldn't have to deal with wrapper
// types related to reference counting.
#[derive(Clone)]
pub struct ConstraintPolynomial<F: Field>(ConstraintPolynomialRef<F>);

impl<F: Field> ConstraintPolynomial<F> {
    pub fn constant(c: F) -> Self {
        Self::from_inner(ConstraintPolynomialInner::Constant(c))
    }

    pub fn constant_usize(c: usize) -> Self {
        Self::constant(F::from_canonical_usize(c))
    }

    pub fn zero() -> Self {
        Self::constant(F::ZERO)
    }

    pub fn one() -> Self {
        Self::constant(F::ONE)
    }

    pub fn local_constant(index: usize) -> Self {
        Self::from_inner(ConstraintPolynomialInner::LocalConstant(index))
    }

    pub fn next_constant(index: usize) -> Self {
        Self::from_inner(ConstraintPolynomialInner::NextConstant(index))
    }

    pub fn local_wire_value(index: usize) -> Self {
        Self::from_inner(ConstraintPolynomialInner::LocalWireValue(index))
    }

    pub fn next_wire_value(index: usize) -> Self {
        Self::from_inner(ConstraintPolynomialInner::NextWireValue(index))
    }

    // TODO: Have these take references?
    pub fn add(&self, rhs: &Self) -> Self {
        // TODO: Special case for either operand being 0.
        Self::from_inner(ConstraintPolynomialInner::Sum {
            lhs: self.0.clone(),
            rhs: rhs.0.clone(),
        })
    }

    pub fn sub(&self, rhs: &Self) -> Self {
        // TODO: Special case for either operand being 0.
        // TODO: Faster to have a dedicated ConstraintPolynomialInner::Difference?
        // TODO: `self + -rhs`?
        self.add(&rhs.neg())
    }

    pub fn double(&self) -> Self {
        self.clone().add(self)
    }

    pub fn triple(&self) -> Self {
        self * 3
    }

    pub fn mul(&self, rhs: &Self) -> Self {
        // TODO: Special case for either operand being 1.
        Self::from_inner(ConstraintPolynomialInner::Product {
            lhs: self.0.clone(),
            rhs: rhs.0.clone(),
        })
    }

    pub fn exp(&self, exponent: usize) -> Self {
        Self::from_inner(ConstraintPolynomialInner::Exponentiation {
            base: self.0.clone(),
            exponent,
        })
    }

    pub fn square(&self) -> Self {
        self * self
    }

    pub(crate) fn degree(&self) -> usize {
        (self.0).0.degree()
    }

    /// Returns the set of wires that this constraint would depend on if it were applied at a
    /// certain gate index.
    pub(crate) fn dependencies(&self, gate: usize) -> Vec<Wire> {
        let mut deps = HashSet::new();
        self.0.0.add_dependencies(gate, &mut deps);
        deps.into_iter().collect()
    }

    /// Find the largest input index among the wires this constraint depends on.
    pub(crate) fn max_wire_input_index(&self) -> Option<usize> {
        self.dependencies(0)
            .into_iter()
            .map(|wire| wire.input)
            .max()
    }

    pub(crate) fn max_constant_index(&self) -> Option<usize> {
        let mut indices = HashSet::new();
        self.0.0.add_constant_indices(&mut indices);
        indices.into_iter().max()
    }

    pub(crate) fn evaluate(&self, vars: EvaluationVars<F>) -> F {
        let results = Self::evaluate_all(&[self.clone()], vars);
        assert_eq!(results.len(), 1);
        results[0]
    }

    /// Evaluate multiple constraint polynomials simultaneously. This can be more efficient than
    /// evaluating them sequentially, since shared intermediate results will only be computed once.
    pub(crate) fn evaluate_all(
        polynomials: &[ConstraintPolynomial<F>],
        vars: EvaluationVars<F>,
    ) -> Vec<F> {
        let mut mem = HashMap::new();
        polynomials.iter()
            .map(|p| p.0.evaluate_memoized(&vars, &mut mem))
            .collect()
    }

    fn from_inner(inner: ConstraintPolynomialInner<F>) -> Self {
        Self(ConstraintPolynomialRef::new(inner))
    }
}

impl<F: Field> Neg for ConstraintPolynomial<F> {
    type Output = Self;

    fn neg(self) -> Self {
        // TODO: Faster to have a dedicated ConstraintPolynomialInner::Negation?
        self * ConstraintPolynomial::constant(F::NEG_ONE)
    }
}

impl<F: Field> Neg for &ConstraintPolynomial<F> {
    type Output = ConstraintPolynomial<F>;

    fn neg(self) -> ConstraintPolynomial<F> {
        self.clone().neg()
    }
}

/// Generates the following variants of a binary operation:
/// - `Self . Self`
/// - `&Self . Self`
/// - `Self . &Self`
/// - `&Self . &Self`
/// - `Self . F`
/// - `&Self . F`
/// - `Self . usize`
/// - `&Self . usize`
/// where `Self` is `ConstraintPolynomial<F>`.
///
/// Takes the following arguments:
/// - `$trait`: the name of the binary operation trait to implement
/// - `$method`: the name of the method in the trait. It is assumed that `ConstraintPolynomial`
///   contains a method with the same name, implementing the `Self . Self` variant.
macro_rules! binop_variants {
    ($trait:ident, $method:ident) => {
        impl<F: Field> $trait<Self> for ConstraintPolynomial<F> {
            type Output = Self;

            fn $method(self, rhs: Self) -> Self {
                ConstraintPolynomial::$method(&self, &rhs)
            }
        }

        impl<F: Field> $trait<&Self> for ConstraintPolynomial<F> {
            type Output = Self;

            fn $method(self, rhs: &Self) -> Self {
                ConstraintPolynomial::$method(&self, rhs)
            }
        }

        impl<F: Field> $trait<ConstraintPolynomial<F>> for &ConstraintPolynomial<F> {
            type Output = ConstraintPolynomial<F>;

            fn $method(self, rhs: ConstraintPolynomial<F>) -> Self::Output {
                ConstraintPolynomial::$method(self, &rhs)
            }
        }

        impl<F: Field> $trait for &ConstraintPolynomial<F> {
            type Output = ConstraintPolynomial<F>;

            fn $method(self, rhs: Self) -> Self::Output {
                ConstraintPolynomial::$method(self, rhs)
            }
        }

        impl<F: Field> $trait<F> for ConstraintPolynomial<F> {
            type Output = Self;

            fn $method(self, rhs: F) -> Self {
                ConstraintPolynomial::$method(&self, &ConstraintPolynomial::constant(rhs))
            }
        }

        impl<F: Field> $trait<F> for &ConstraintPolynomial<F> {
            type Output = ConstraintPolynomial<F>;

            fn $method(self, rhs: F) -> Self::Output {
                ConstraintPolynomial::$method(self, &ConstraintPolynomial::constant(rhs))
            }
        }

        impl<F: Field> $trait<usize> for ConstraintPolynomial<F> {
            type Output = Self;

            fn $method(self, rhs: usize) -> Self {
                ConstraintPolynomial::$method(&self, &ConstraintPolynomial::constant_usize(rhs))
            }
        }

        impl<F: Field> $trait<usize> for &ConstraintPolynomial<F> {
            type Output = ConstraintPolynomial<F>;

            fn $method(self, rhs: usize) -> Self::Output {
                ConstraintPolynomial::$method(self, &ConstraintPolynomial::constant_usize(rhs))
            }
        }
    };
}

binop_variants!(Add, add);
binop_variants!(Sub, sub);
binop_variants!(Mul, mul);

impl<F: Field> Sum for ConstraintPolynomial<F> {
    fn sum<I: Iterator<Item=Self>>(iter: I) -> Self {
        iter.fold(
            ConstraintPolynomial::zero(),
            |sum, x| sum + x)
    }
}

impl<F: Field> Product for ConstraintPolynomial<F> {
    fn product<I: Iterator<Item=Self>>(iter: I) -> Self {
        iter.fold(
            ConstraintPolynomial::one(),
            |product, x| product * x)
    }
}

enum ConstraintPolynomialInner<F: Field> {
    Constant(F),

    LocalConstant(usize),
    NextConstant(usize),
    LocalWireValue(usize),
    NextWireValue(usize),

    Sum {
        lhs: ConstraintPolynomialRef<F>,
        rhs: ConstraintPolynomialRef<F>,
    },
    Product {
        lhs: ConstraintPolynomialRef<F>,
        rhs: ConstraintPolynomialRef<F>,
    },
    Exponentiation {
        base: ConstraintPolynomialRef<F>,
        exponent: usize,
    },
}

impl<F: Field> ConstraintPolynomialInner<F> {
    fn add_dependencies(&self, gate: usize, deps: &mut HashSet<Wire>) {
        match self {
            ConstraintPolynomialInner::Constant(_) => (),
            ConstraintPolynomialInner::LocalConstant(_) => (),
            ConstraintPolynomialInner::NextConstant(_) => (),
            ConstraintPolynomialInner::LocalWireValue(i) =>
                { deps.insert(Wire { gate, input: *i }); },
            ConstraintPolynomialInner::NextWireValue(i) =>
                { deps.insert(Wire { gate: gate + 1, input: *i }); }
            ConstraintPolynomialInner::Sum { lhs, rhs } => {
                lhs.0.add_dependencies(gate, deps);
                rhs.0.add_dependencies(gate, deps);
            },
            ConstraintPolynomialInner::Product { lhs, rhs } => {
                lhs.0.add_dependencies(gate, deps);
                rhs.0.add_dependencies(gate, deps);
            },
            ConstraintPolynomialInner::Exponentiation { base, exponent: _ } => {
                base.0.add_dependencies(gate, deps);
            },
        }
    }

    fn add_constant_indices(&self, indices: &mut HashSet<usize>) {
        match self {
            ConstraintPolynomialInner::Constant(_) => (),
            ConstraintPolynomialInner::LocalConstant(i) => { indices.insert(*i); },
            ConstraintPolynomialInner::NextConstant(i) => { indices.insert(*i); },
            ConstraintPolynomialInner::LocalWireValue(_) => (),
            ConstraintPolynomialInner::NextWireValue(_) => (),
            ConstraintPolynomialInner::Sum { lhs, rhs } => {
                lhs.0.add_constant_indices(indices);
                rhs.0.add_constant_indices(indices);
            },
            ConstraintPolynomialInner::Product { lhs, rhs } => {
                lhs.0.add_constant_indices(indices);
                rhs.0.add_constant_indices(indices);
            },
            ConstraintPolynomialInner::Exponentiation { base, exponent: _ } => {
                base.0.add_constant_indices(indices);
            },
        }
    }

    fn evaluate(
        &self,
        vars: &EvaluationVars<F>,
        mem: &mut HashMap<ConstraintPolynomialRef<F>, F>,
    ) -> F {
        match self {
            ConstraintPolynomialInner::Constant(c) => *c,
            ConstraintPolynomialInner::LocalConstant(i) => vars.local_constants[*i],
            ConstraintPolynomialInner::NextConstant(i) => vars.next_constants[*i],
            ConstraintPolynomialInner::LocalWireValue(i) => vars.local_wire_values[*i],
            ConstraintPolynomialInner::NextWireValue(i) => vars.next_wire_values[*i],
            ConstraintPolynomialInner::Sum { lhs, rhs } => {
                let lhs = lhs.evaluate_memoized(vars, mem);
                let rhs = rhs.evaluate_memoized(vars, mem);
                lhs + rhs
            },
            ConstraintPolynomialInner::Product { lhs, rhs } => {
                let lhs = lhs.evaluate_memoized(vars, mem);
                let rhs = rhs.evaluate_memoized(vars, mem);
                lhs * rhs
            },
            ConstraintPolynomialInner::Exponentiation { base, exponent } => {
                let base = base.evaluate_memoized(vars, mem);
                base.exp_usize(*exponent)
            },
        }
    }

    fn degree(&self) -> usize {
        match self {
            ConstraintPolynomialInner::Constant(_) => 0,
            ConstraintPolynomialInner::LocalConstant(_) => 1,
            ConstraintPolynomialInner::NextConstant(_) => 1,
            ConstraintPolynomialInner::LocalWireValue(_) => 1,
            ConstraintPolynomialInner::NextWireValue(_) => 1,
            ConstraintPolynomialInner::Sum { lhs, rhs } => lhs.0.degree().max(rhs.0.degree()),
            ConstraintPolynomialInner::Product { lhs, rhs } => lhs.0.degree() + rhs.0.degree(),
            ConstraintPolynomialInner::Exponentiation { base, exponent } => base.0.degree() * exponent,
        }
    }
}

/// Wraps `Rc<ConstraintPolynomialRef>`, and implements `Hash` and `Eq` based on references rather
/// than content. This is useful when we want to use constraint polynomials as `HashMap` keys, but
/// we want address-based hashing for performance reasons.
#[derive(Clone)]
struct ConstraintPolynomialRef<F: Field>(Rc<ConstraintPolynomialInner<F>>);

impl<F: Field> ConstraintPolynomialRef<F> {
    fn new(inner: ConstraintPolynomialInner<F>) -> Self {
        Self(Rc::new(inner))
    }

    fn evaluate_memoized(
        &self,
        vars: &EvaluationVars<F>,
        mem: &mut HashMap<Self, F>,
    ) -> F {
        if let Some(&result) = mem.get(self) {
            result
        } else {
            let result = self.0.evaluate(vars, mem);
            mem.insert(self.clone(), result);
            result
        }
    }
}

impl<F: Field> PartialEq for ConstraintPolynomialRef<F> {
    fn eq(&self, other: &Self) -> bool {
        ptr::eq(&*self.0, &*other.0)
    }
}

impl<F: Field> Eq for ConstraintPolynomialRef<F> {}

impl<F: Field> Hash for ConstraintPolynomialRef<F> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        ptr::hash(&*self.0, state);
    }
}
