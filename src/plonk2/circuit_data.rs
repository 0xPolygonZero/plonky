use crate::{FftPrecomputation, MsmPrecomputation, AffinePoint, Curve, Proof2};
use crate::plonk2::prover::prove2;
use crate::plonk2::verifier::verify2;

#[derive(Copy, Clone)]
pub struct CircuitConfig {
    pub num_wires: usize,
    pub num_routed_wires: usize,
    pub security_bits: usize,
}

impl CircuitConfig {
    pub fn advice_wires(&self) -> usize {
        self.num_wires - self.num_routed_wires
    }
}

/// Circuit data required by the prover or the verifier.
pub struct CircuitData<C: Curve> {
    prover_only: ProverOnlyCircuitData<C>,
    verifier_only: VerifierOnlyCircuitData,
    common: CommonCircuitData<C>,
}

impl<C: Curve> CircuitData<C> {
    pub fn prove2(&self) -> Proof2 {
        prove2(&self.prover_only, &self.common)
    }

    pub fn verify2(&self) {
        verify2(&self.verifier_only, &self.common)
    }
}

/// Circuit data required by the prover.
pub struct ProverCircuitData<C: Curve> {
    prover_only: ProverOnlyCircuitData<C>,
    common: CommonCircuitData<C>,
}

impl<C: Curve> ProverCircuitData<C> {
    pub fn prove2(&self) -> Proof2 {
        prove2(&self.prover_only, &self.common)
    }
}

/// Circuit data required by the prover.
pub struct VerifierCircuitData<C: Curve> {
    verifier_only: VerifierOnlyCircuitData,
    common: CommonCircuitData<C>,
}

impl<C: Curve> VerifierCircuitData<C> {
    pub fn verify2(&self) {
        verify2(&self.verifier_only, &self.common)
    }
}

/// Circuit data required by the prover, but not the verifier.
pub(crate) struct ProverOnlyCircuitData<C: Curve> {
    /// A precomputation used for FFTs of degree 8n, where n is the number of gates.
    pub fft_precomputation_8n: FftPrecomputation<C::ScalarField>,

    /// A precomputation used for MSMs involving `generators`.
    pub pedersen_g_msm_precomputation: MsmPrecomputation<C>,
}

/// Circuit data required by the verifier, but not the prover.
pub(crate) struct VerifierOnlyCircuitData {}

/// Circuit data required by both the prover and the verifier.
pub(crate) struct CommonCircuitData<C: Curve> {
    pub config: CircuitConfig,

    pub degree: usize,

    /// A commitment to each constant polynomial.
    pub c_constants: Vec<AffinePoint<C>>,

    /// A commitment to each permutation polynomial.
    pub c_s_sigmas: Vec<AffinePoint<C>>,

    /// A precomputation used for MSMs involving `generators`.
    pub pedersen_g_msm_precomputation: MsmPrecomputation<C>,

    /// A precomputation used for FFTs of degree n, where n is the number of gates.
    pub fft_precomputation_n: FftPrecomputation<C::ScalarField>,
}
