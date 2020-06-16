use std::marker::PhantomData;

use crate::{
    AffinePoint, AffinePointTarget, CircuitBuilder, Curve, Field, HaloCurve, ProjectivePoint,
    rescue_sponge, Target,
};

/// Observes prover messages, and generates challenges by hashing the transcript.
pub(crate) struct Challenger<F: Field> {
    transcript: Vec<F>,
    security_bits: usize,
}

impl<F: Field> Challenger<F> {
    pub(crate) fn new(security_bits: usize) -> Challenger<F> {
        Challenger {
            transcript: Vec::new(),
            security_bits,
        }
    }

    pub(crate) fn observe_element(&mut self, element: F) {
        self.transcript.push(element);
    }

    pub(crate) fn observe_elements(&mut self, elements: &[F]) {
        for &element in elements {
            self.observe_element(element);
        }
    }

    pub(crate) fn observe_affine_point<C: Curve<BaseField=F>>(&mut self, point: AffinePoint<C>) {
        // debug_assert!(!point.zero);
        self.observe_element(point.x);
        self.observe_element(point.y);
    }

    pub(crate) fn observe_affine_points<C: Curve<BaseField=F>>(
        &mut self,
        points: &[AffinePoint<C>],
    ) {
        for &point in points {
            self.observe_affine_point(point);
        }
    }

    pub(crate) fn observe_proj_point<C: Curve<BaseField=F>>(
        &mut self,
        point: ProjectivePoint<C>,
    ) {
        self.observe_affine_point(point.to_affine());
    }

    pub(crate) fn observe_proj_points<C: Curve<BaseField=F>>(
        &mut self,
        points: &[ProjectivePoint<C>],
    ) {
        self.observe_affine_points(&ProjectivePoint::batch_to_affine(points));
    }

    pub(crate) fn observe_proj_point_other_curve<C: Curve<ScalarField=F>>(
        &mut self,
        point: ProjectivePoint<C>,
    ) {
        let p = point.to_affine();
        let x = p.x.try_convert::<F>().expect("Element is too large");
        let y = p.x.try_convert::<F>().expect("Element is too large");
        self.observe_elements(&[x, y]);
    }

    pub(crate) fn get_challenge(&mut self) -> F {
        self.get_n_challenges(1)[0]
    }

    pub(crate) fn get_2_challenges(&mut self) -> (F, F) {
        let challenges = self.get_n_challenges(2);
        (challenges[0], challenges[1])
    }

    pub(crate) fn get_3_challenges(&mut self) -> (F, F, F) {
        let challenges = self.get_n_challenges(3);
        (challenges[0], challenges[1], challenges[2])
    }

    pub(crate) fn get_n_challenges(&mut self, n: usize) -> Vec<F> {
        let challenges = rescue_sponge(self.transcript.clone(), n, self.security_bits);
        self.transcript = vec![challenges[0]];
        challenges
    }
}

/// Observes prover messages, and generates challenges by hashing the transcript.
pub(crate) struct RecursiveChallenger<C: HaloCurve> {
    transcript: Vec<Target>,
    _phantom: PhantomData<C>,
}

impl<C: HaloCurve> RecursiveChallenger<C> {
    pub(crate) fn new() -> RecursiveChallenger<C> {
        RecursiveChallenger {
            transcript: Vec::new(),
            _phantom: PhantomData,
        }
    }

    pub(crate) fn observe_element(&mut self, target: Target) {
        self.transcript.push(target);
    }

    pub(crate) fn observe_elements(&mut self, targets: &[Target]) {
        for &target in targets {
            self.observe_element(target);
        }
    }

    pub(crate) fn observe_affine_point(&mut self, point: AffinePointTarget) {
        self.observe_element(point.x);
        self.observe_element(point.y);
    }

    pub(crate) fn observe_affine_points(&mut self, points: &[AffinePointTarget]) {
        for &point in points {
            self.observe_affine_point(point);
        }
    }

    pub(crate) fn get_challenge(&mut self, builder: &mut CircuitBuilder<C>) -> Target {
        self.get_n_challenges(builder, 1)[0]
    }

    pub(crate) fn get_2_challenges(&mut self, builder: &mut CircuitBuilder<C>) -> (Target, Target) {
        let challenges = self.get_n_challenges(builder, 2);
        (challenges[0], challenges[1])
    }

    pub(crate) fn get_3_challenges(
        &mut self,
        builder: &mut CircuitBuilder<C>,
    ) -> (Target, Target, Target) {
        let challenges = self.get_n_challenges(builder, 3);
        (challenges[0], challenges[1], challenges[2])
    }

    pub(crate) fn get_n_challenges(
        &mut self,
        builder: &mut CircuitBuilder<C>,
        n: usize,
    ) -> Vec<Target> {
        let challenges = builder.rescue_sponge(&self.transcript, n);
        self.transcript = vec![challenges[0]];
        challenges
    }

    // TODO: Implement recursive `get_affine_point`.
}
