use crate::{Field, ConstraintPolynomial};
use std::any::TypeId;
use once_cell::sync::Lazy;
use std::sync::Mutex;
use std::collections::HashMap;
use std::marker::PhantomData;

static CACHED_MDS_MATRICES: Lazy<Mutex<HashMap<MdsMatrixKey, UnparameterizedMdsMatrix>>> = Lazy::new(|| {
    Mutex::new(HashMap::new())
});

/// A key for looking up a cached MDS matrix.
#[derive(Eq, PartialEq, Hash)]
struct MdsMatrixKey {
    field_type_id: TypeId,
    size: usize,
}

impl MdsMatrixKey {
    fn new<F: Field>(size: usize) -> MdsMatrixKey {
        MdsMatrixKey { field_type_id: TypeId::of::<F>(), size }
    }
}

pub struct MdsMatrix<F: Field> {
    unparameterized: UnparameterizedMdsMatrix,
    _phantom: PhantomData<F>,
}

impl<F: Field> MdsMatrix<F> {
    /// Returns the width and height of this matrix.
    pub fn size(&self) -> usize {
        self.unparameterized.rows.len()
    }

    /// Returns the entry at row `r` and column `c`.
    pub fn get(&self, r: usize, c: usize) -> F {
        F::from_canonical_u64_vec(self.unparameterized.rows[r][c].clone())
    }
}

/// A representation of an MDS matrix which does not involve Field types.
#[derive(Clone)]
struct UnparameterizedMdsMatrix {
    rows: Vec<Vec<Vec<u64>>>,
}

/// Apply an MDS matrix to the given vector of field elements.
pub(crate) fn apply_mds<F: Field>(vec: Vec<F>) -> Vec<F> {
    let n = vec.len();
    let mds = mds_matrix::<F>(n);

    (0..n)
        .map(|r| (0..n)
            .map(|c| mds.get(r, c) * vec[c])
            .fold(F::ZERO, |acc, x| acc + x))
        .collect()
}

/// Applies an MDS matrix to the given vector of constraint polynomials.
pub(crate) fn apply_mds_constraint_polys<F: Field>(
    vec: Vec<ConstraintPolynomial<F>>,
) -> Vec<ConstraintPolynomial<F>> {
    let n = vec.len();
    let mds = mds_matrix::<F>(n);

    (0..n)
        .map(|r| (0..n)
            .map(|c| &vec[c] * mds.get(r, c))
            .sum())
        .collect()
}

/// Returns entry `(r, c)` of an `n` by `n` MDS matrix.
pub(crate) fn mds_matrix<F: Field>(n: usize) -> MdsMatrix<F> {
    let mut cached_matrices = CACHED_MDS_MATRICES.lock().unwrap();
    let key = MdsMatrixKey::new::<F>(n);
    let unparameterized = cached_matrices.entry(key).or_insert_with(|| generate_mds_matrix::<F>(n)).clone();
    MdsMatrix { unparameterized, _phantom: PhantomData }
}

fn generate_mds_matrix<F: Field>(n: usize) -> UnparameterizedMdsMatrix {
    let mut rows: Vec<Vec<Vec<u64>>> = Vec::new();
    for r in 0..n {
        let mut row = Vec::new();
        for c in 0..n {
            // We use a Cauchy matrix with x_r = n + r, y_c = c.
            let x = F::from_canonical_usize(n + r);
            let y = F::from_canonical_usize(c);
            let entry = (x - y).multiplicative_inverse().unwrap().to_canonical_u64_vec();
            row.push(entry);
        }
        rows.push(row);
    }
    UnparameterizedMdsMatrix { rows }
}
