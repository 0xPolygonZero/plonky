use crate::Bls12377Scalar;

/// A Plonk constraint has the form `L a + R b + O c + M a b + C = 0`.
struct PlonkConstraint {
    /// The coefficient of a, the left input.
    l: Bls12377Scalar,
    /// The coefficient of b, the right input.
    r: Bls12377Scalar,
    /// The coefficient of c, the output.
    o: Bls12377Scalar,
    /// The coefficient of ab, the product.
    m: Bls12377Scalar,
    /// The constant which stands on its own.
    c: Bls12377Scalar,
}

struct Witness {
    l: Vec<Bls12377Scalar>,
    r: Vec<Bls12377Scalar>,
    o: Vec<Bls12377Scalar>,
}

impl Witness {
    pub fn num_gates(&self) -> usize {
        self.l.len()
    }

    /// Get the i'th wire within the concatenated list (a, b, c).
    fn get_wire(&self, mut i: usize) -> Bls12377Scalar {
        let n = self.num_gates();
        if i < n {
            return self.l[i];
        }
        i -= n;
        if i < n {
            return self.r[i];
        }
        i -= n;
        self.o[i]
    }
}
