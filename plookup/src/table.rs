use crate::plookup::{prove, SECURITY_BITS};
use crate::proof::PlookupProof;
use crate::verifier::verify;
use anyhow::Result;
use itertools::Itertools;
use plonky::plonk_challenger::Challenger;
use plonky::{Field, HaloCurve};

/// A `Table` is a list of rows of field elements, all with size `N`.
pub struct Table<F: Field, const N: usize>(pub Vec<[F; N]>);

impl<F: Field, const N: usize> Table<F, N> {
    /// Constructs a `Table` with width `N` from a function `F^(N-1) -> F` and an evaluation domain
    /// `domain \subseteq F^(N-1)`.
    pub fn from_function<const M: usize>(f: &dyn Fn([F; M]) -> F, domain: &[[F; M]]) -> Self {
        assert_eq!(M, N - 1);
        Self(
            domain
                .iter()
                .map(|&a| {
                    let mut arr = [F::ZERO; N];
                    for (x, &y) in arr.iter_mut().zip(a.iter()) {
                        *x = y;
                    }
                    arr[N - 1] = f(a);
                    arr
                })
                .collect::<Vec<_>>(),
        )
    }

    /// Constructs a `Table` with width `N` from a function `F^(N-1) -> F` and an evaluation domain
    /// `domain^(N-1) \subseteq F^(N-1)`.
    pub fn from_function_cartesian<const M: usize>(f: &dyn Fn([F; M]) -> F, domain: &[F]) -> Self {
        assert_eq!(M, N - 1);
        let cartesian_domain = (0..M).map(|_| domain.iter()).multi_cartesian_product();
        Self(
            cartesian_domain
                .map(|a| {
                    let mut arr = [F::ZERO; N];
                    for (x, &&y) in arr.iter_mut().zip(a.iter()) {
                        *x = y;
                    }
                    let mut src = [F::ZERO; M];
                    src.copy_from_slice(&a.iter().map(|x| **x).collect::<Vec<_>>());
                    arr[N - 1] = f(src);
                    arr
                })
                .collect::<Vec<_>>(),
        )
    }

    /// Get a verifier challenge from a `Table`. Used to reduce a `Table` to a vector.
    fn get_challenge(&self) -> F {
        let mut challenger = Challenger::new(SECURITY_BITS);
        for a in &self.0 {
            challenger.observe_elements(a);
        }
        challenger.get_challenge()
    }

    /// Reduces a `Table` to a vector using a random challenge.
    pub fn to_vec(&self) -> Vec<F> {
        let alpha = self.get_challenge();
        self.0
            .iter()
            .map(|a| a.iter().fold(F::ZERO, |acc, &x| alpha * acc + x))
            .collect()
    }

    /// Reduces a `Table` and a subset thereof to two vectors using a random challenge.
    fn to_vecs_with_row_witnesses(&self, ws: &Self) -> (Vec<F>, Vec<F>) {
        let alpha = self.get_challenge();
        (
            self.0
                .iter()
                .map(|a| a.iter().fold(F::ZERO, |acc, &x| alpha * acc + x))
                .collect(),
            ws.0.iter()
                .map(|a| a.iter().fold(F::ZERO, |acc, &x| alpha * acc + x))
                .collect(),
        )
    }

    /// Reduces a `Table` and vector of columns to two vectors using a random challenge.
    fn to_vecs_with_column_witnesses(&self, ws: &[Vec<F>; N]) -> (Vec<F>, Vec<F>) {
        let alpha = self.get_challenge();
        let h = ws[0].len();
        (
            self.0
                .iter()
                .map(|a| a.iter().fold(F::ZERO, |acc, &x| alpha * acc + x))
                .collect(),
            (0..h)
                .map(|i| {
                    (0..N)
                        .map(|j| ws[j][i])
                        .fold(F::ZERO, |acc, x| alpha * acc + x)
                })
                .collect(),
        )
    }

    /// Proves that `ws` is a subtable of `self` using the Plookup protocol.
    pub fn prove_row<C: HaloCurve<ScalarField = F>>(&self, ws: &Self) -> Result<PlookupProof<C>> {
        let (t, f) = self.to_vecs_with_row_witnesses(ws);
        prove(&f, &t)
    }

    /// Proves that the columns `ws` form a subtable of `self` using the Plookup protocol.
    pub fn prove_column<C: HaloCurve<ScalarField = F>>(
        &self,
        ws: &[Vec<F>; N],
    ) -> Result<PlookupProof<C>> {
        let (t, f) = self.to_vecs_with_column_witnesses(ws);
        prove(&f, &t)
    }

    /// Verifies that a proof is valid for a table `t`.
    pub fn verify<C: HaloCurve<ScalarField = F>>(&self, proof: &PlookupProof<C>) -> Result<()> {
        verify(&self.to_vec(), proof)
    }
}
