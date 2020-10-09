use crate::{ArithmeticGate, Base4SumGate, BufferGate, ConstantGate, CurveAddGate, CurveDblGate, CurveEndoGate, Gate, HaloCurve, PublicInputGate, RescueStepAGate, RescueStepBGate};

// TODO: Can this impl usize?
pub(crate) fn ceil_div_usize(a: usize, b: usize) -> usize {
    (a + b - 1) / b
}

pub(crate) fn pad_to_multiple_usize(a: usize, b: usize) -> usize {
    ceil_div_usize(a, b) * b
}

/// Computes `ceil(log_2(n))`.
pub(crate) fn log2_ceil(n: usize) -> usize {
    n.next_power_of_two().trailing_zeros() as usize
}

/// Computes `log_2(n)`, panicking if `n` is not a power of two.
pub fn log2_strict(n: usize) -> usize {
    assert!(n.is_power_of_two(), "Not a power of two");
    log2_ceil(n)
}

pub(crate) fn transpose<T: Clone>(matrix: &[Vec<T>]) -> Vec<Vec<T>> {
    let old_rows = matrix.len();
    let old_cols = matrix[0].len();
    let mut transposed = vec![Vec::with_capacity(old_rows); old_cols];
    for new_r in 0..old_cols {
        for new_c in 0..old_rows {
            transposed[new_r].push(matrix[new_c][new_r].clone());
        }
    }
    transposed
}

// Needed because of issues related to https://github.com/rust-lang/rust/issues/61083
pub struct GateContainer<C: HaloCurve> {
    pub gate: Box<dyn Gate<C>>,
}

impl<C: HaloCurve> GateContainer<C> {
    pub fn name(&self) -> &'static str {
        self.gate.name()
    }
    pub fn degree(&self) -> usize {
        self.gate.degree()
    }
    pub fn num_constants(&self) -> usize {
        self.gate.num_constants()
    }
    pub fn prefix(&self) -> &'static [bool] {
        self.gate.prefix()
    }
}

pub fn get_canonical_gates<C: HaloCurve, InnerC: HaloCurve<BaseField = C::ScalarField>>(
) -> Vec<Box<dyn Gate<C>>> {
    vec![
        Box::new(CurveAddGate::<C, InnerC>::new(0)),
        Box::new(CurveDblGate::<C, InnerC>::new(0)),
        Box::new(CurveEndoGate::<C, InnerC>::new(0)),
        Box::new(Base4SumGate::<C>::new(0)),
        Box::new(PublicInputGate::<C>::new(0)),
        Box::new(BufferGate::<C>::new(0)),
        Box::new(ConstantGate::<C>::new(0)),
        Box::new(ArithmeticGate::<C>::new(0)),
        Box::new(RescueStepAGate::<C>::new(0)),
        Box::new(RescueStepBGate::<C>::new(0)),
    ]
}

pub fn get_canonical_gates_containers<
    C: HaloCurve,
    InnerC: HaloCurve<BaseField = C::ScalarField>,
>() -> Vec<GateContainer<C>> {
    vec![
        GateContainer {
            gate: Box::new(CurveAddGate::<C, InnerC>::new(0)),
        },
        GateContainer {
            gate: Box::new(CurveDblGate::<C, InnerC>::new(0)),
        },
        GateContainer {
            gate: Box::new(CurveEndoGate::<C, InnerC>::new(0)),
        },
        GateContainer {
            gate: Box::new(Base4SumGate::<C>::new(0)),
        },
        GateContainer {
            gate: Box::new(PublicInputGate::<C>::new(0)),
        },
        GateContainer {
            gate: Box::new(BufferGate::<C>::new(0)),
        },
        GateContainer {
            gate: Box::new(ConstantGate::<C>::new(0)),
        },
        GateContainer {
            gate: Box::new(ArithmeticGate::<C>::new(0)),
        },
        GateContainer {
            gate: Box::new(RescueStepAGate::<C>::new(0)),
        },
        GateContainer {
            gate: Box::new(RescueStepBGate::<C>::new(0)),
        },
    ]
}
