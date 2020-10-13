#![allow(clippy::type_complexity)]
use crate::gates::gate_collection::GatePrefixes;
use crate::plonk_util::halo_n;
use crate::{blake_hash_base_field_to_curve, AffinePoint, Base4SumGate, BufferGate, CircuitBuilder, Curve, CurveAddGate, CurveDblGate, CurveEndoGate, Field, HaloCurve, PartialWitness, Target, Wire, WitnessGenerator};
use std::marker::PhantomData;

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct AffinePointTarget<C: Curve> {
    pub x: Target<C::BaseField>,
    pub y: Target<C::BaseField>,
}

impl<C: Curve> AffinePointTarget<C> {
    pub fn to_vec(&self) -> Vec<Target<C::BaseField>> {
        vec![self.x, self.y]
    }
}

/// Represents a scalar * point multiplication operation on `InnerC`.
/// `scalar` is modelled here in the "wrong" field `InnerC::BaseField = C::ScalarField` for coherence.
/// Thus, all scalar operations should be done preemptively in the correct field `InnerC::ScalarField = C::BaseField".
pub struct CurveMulOp<C: Curve, InnerC: Curve<BaseField = C::ScalarField>> {
    pub scalar: Target<C::ScalarField>,
    pub point: AffinePointTarget<InnerC>,
}

pub struct CurveMulEndoResult<C: Curve, InnerC: Curve<BaseField = C::ScalarField>> {
    pub mul_result: AffinePointTarget<InnerC>,
    pub actual_scalar: Target<C::ScalarField>,
}

pub struct CurveMsmEndoResult<C: Curve, InnerC: Curve<BaseField = C::ScalarField>> {
    pub msm_result: AffinePointTarget<InnerC>,
    /// While `msm` computes a sum of `[s] P` terms, `msm_endo` computes a sum of `[n(s)] P` terms
    /// for some injective `n`. Here we return each `n(s)`, i.e., the scalar by which the point was
    /// actually multiplied.
    pub actual_scalars: Vec<Target<C::ScalarField>>,
}

impl<C: HaloCurve> CircuitBuilder<C> {
    pub fn constant_affine_point<InnerC: Curve<BaseField = C::ScalarField>>(
        &mut self,
        point: AffinePoint<InnerC>,
    ) -> AffinePointTarget<InnerC> {
        assert!(!point.zero);
        AffinePointTarget {
            x: self.constant_wire(point.x),
            y: self.constant_wire(point.y),
        }
    }

    /// Add a copy constraint between two affine targets.
    pub fn copy_curve<InnerC: Curve<BaseField = C::ScalarField>>(
        &mut self,
        affine_target_1: AffinePointTarget<InnerC>,
        affine_target_2: AffinePointTarget<InnerC>,
    ) {
        self.copy(affine_target_1.x, affine_target_2.x);
        self.copy(affine_target_1.y, affine_target_2.y);
    }

    /// Assert that a given coordinate pair is on the curve `C`.
    pub fn curve_assert_valid<InnerC: Curve<BaseField = C::ScalarField>>(
        &mut self,
        p: AffinePointTarget<InnerC>,
    ) {
        // Recall the short Weierstrass equation: y^2 = x^3 + a*x + b.
        let a = self.constant_wire(InnerC::A);
        let b = self.constant_wire(InnerC::B);
        let y_squared = self.square(p.y);
        let x_cubed = self.exp_constant_usize(p.x, 3);
        let a_x_plus_b = self.mul_add(a, p.x, b);
        let rhs = self.add(x_cubed, a_x_plus_b);
        self.copy(y_squared, rhs);
    }

    pub fn curve_neg<InnerC: Curve<BaseField = C::ScalarField>>(
        &mut self,
        p: AffinePointTarget<InnerC>,
    ) -> AffinePointTarget<InnerC> {
        let neg_y = self.neg(p.y);
        AffinePointTarget { x: p.x, y: neg_y }
    }

    pub fn curve_add<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &mut self,
        p_1: AffinePointTarget<InnerC>,
        p_2: AffinePointTarget<InnerC>,
    ) -> AffinePointTarget<InnerC> {
        // Add a CurveAddGate, then add a BufferGate to receive the updated accumulator state.
        let add_index = self.num_gates();
        self.add_gate_no_constants(CurveAddGate::<C, InnerC>::new(add_index));
        let buffer_index = self.num_gates();
        self.add_gate_no_constants(BufferGate::new(buffer_index));

        let group_acc_x = Target::Wire(Wire {
            gate: add_index,
            input: CurveAddGate::<C, InnerC>::WIRE_GROUP_ACC_X,
        });
        let group_acc_y = Target::Wire(Wire {
            gate: add_index,
            input: CurveAddGate::<C, InnerC>::WIRE_GROUP_ACC_Y,
        });
        let scalar_acc_old = Target::Wire(Wire {
            gate: add_index,
            input: CurveAddGate::<C, InnerC>::WIRE_SCALAR_ACC_OLD,
        });
        let scalar_acc_new = Target::Wire(Wire {
            gate: add_index,
            input: CurveAddGate::<C, InnerC>::WIRE_SCALAR_ACC_NEW,
        });
        let addend_x = Target::Wire(Wire {
            gate: add_index,
            input: CurveAddGate::<C, InnerC>::WIRE_ADDEND_X,
        });
        let addend_y = Target::Wire(Wire {
            gate: add_index,
            input: CurveAddGate::<C, InnerC>::WIRE_ADDEND_Y,
        });
        let scalar_bit = Target::Wire(Wire {
            gate: add_index,
            input: CurveAddGate::<C, InnerC>::WIRE_SCALAR_BIT,
        });
        let result_x = Target::Wire(Wire {
            gate: buffer_index,
            input: CurveAddGate::<C, InnerC>::WIRE_GROUP_ACC_X,
        });
        let result_y = Target::Wire(Wire {
            gate: buffer_index,
            input: CurveAddGate::<C, InnerC>::WIRE_GROUP_ACC_Y,
        });

        // Wire inputs
        self.copy(group_acc_x, p_1.x);
        self.copy(group_acc_y, p_1.y);
        self.copy(addend_x, p_2.x);
        self.copy(addend_y, p_2.y);

        // The scalar bit should always be 1, since we always want to perform the add.
        self.generate_constant(scalar_bit, C::ScalarField::ONE);

        // The scalar accumulator should change from 0 to 1. This enforces that scalar_bit = 1.
        let zero = self.zero_wire();
        let one = self.one_wire();
        self.copy(scalar_acc_old, zero);
        self.copy(scalar_acc_new, one);

        AffinePointTarget {
            x: result_x,
            y: result_y,
        }
    }

    pub fn curve_double<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &mut self,
        p: AffinePointTarget<InnerC>,
    ) -> AffinePointTarget<InnerC> {
        let idx_dbl = self.num_gates();
        self.add_gate_no_constants(CurveDblGate::<C, InnerC>::new(idx_dbl));
        self.copy(
            p.x,
            Target::Wire(Wire {
                gate: idx_dbl,
                input: CurveDblGate::<C, InnerC>::WIRE_X_OLD,
            }),
        );
        self.copy(
            p.y,
            Target::Wire(Wire {
                gate: idx_dbl,
                input: CurveDblGate::<C, InnerC>::WIRE_Y_OLD,
            }),
        );
        AffinePointTarget {
            x: Target::Wire(Wire {
                gate: idx_dbl,
                input: CurveDblGate::<C, InnerC>::WIRE_X_NEW,
            }),
            y: Target::Wire(Wire {
                gate: idx_dbl,
                input: CurveDblGate::<C, InnerC>::WIRE_Y_NEW,
            }),
        }
    }

    pub fn curve_sub<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &mut self,
        p_1: AffinePointTarget<InnerC>,
        p_2: AffinePointTarget<InnerC>,
    ) -> AffinePointTarget<InnerC> {
        let neg_p_2 = self.curve_neg::<InnerC>(p_2);
        self.curve_add::<InnerC>(p_1, neg_p_2)
    }

    pub fn curve_mul<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &mut self,
        mul: CurveMulOp<C, InnerC>,
    ) -> AffinePointTarget<InnerC> {
        self.curve_msm::<InnerC>(&[mul])
    }

    /// Computes `[n(s)] p`.
    pub fn curve_mul_endo<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &mut self,
        mul: CurveMulOp<C, InnerC>,
    ) -> CurveMulEndoResult<C, InnerC> {
        let result = self.curve_msm_endo::<InnerC>(&[mul]);
        CurveMulEndoResult {
            mul_result: result.msm_result,
            actual_scalar: result.actual_scalars[0],
        }
    }

    /// Computes `[1 / n(s)] p`.
    pub fn curve_mul_inv_endo<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &mut self,
        mul: CurveMulOp<C, InnerC>,
    ) -> CurveMulEndoResult<C, InnerC> {
        // We witness r = [1 / n(s)] p, then verify that [n(s)] r = p, and return r.
        let CurveMulOp { scalar, point } = mul;

        struct ResultGenerator<C: HaloCurve, InnerC: HaloCurve<BaseField = C::ScalarField>> {
            mul: CurveMulOp<C, InnerC>,
            result: AffinePointTarget<InnerC>,
            security_bits: usize,
            _phantom: PhantomData<C>,
        }

        impl<C: HaloCurve, InnerC: HaloCurve<BaseField = C::ScalarField>>
            WitnessGenerator<InnerC::BaseField> for ResultGenerator<C, InnerC>
        {
            fn dependencies(&self) -> Vec<Target<InnerC::BaseField>> {
                vec![
                    self.mul.scalar.convert::<InnerC::BaseField>(),
                    self.mul.point.x,
                    self.mul.point.y,
                ]
            }

            fn generate(
                &self,
                prefixes: &GatePrefixes,
                _constants: &[Vec<InnerC::BaseField>],
                witness: &PartialWitness<InnerC::BaseField>,
            ) -> PartialWitness<InnerC::BaseField> {
                let scalar = witness
                    .get_target(self.mul.scalar.convert())
                    .try_convert::<InnerC::ScalarField>()
                    .expect("Improbable");
                let scalar_bits = &scalar.to_canonical_bool_vec()[..self.security_bits];
                let n_scalar = halo_n::<InnerC>(scalar_bits);
                let n_scalar_inv = n_scalar
                    .multiplicative_inverse()
                    .expect("Can't invert zero");

                let point = witness.get_point_target(self.mul.point);
                let result = (InnerC::convert(n_scalar_inv) * point.to_projective()).to_affine();

                let mut result_witness = PartialWitness::new();
                result_witness.set_point_target(self.result, result);
                result_witness
            }
        }

        let result = self.add_virtual_point_target();
        self.add_generator(ResultGenerator::<C, InnerC> {
            mul,
            result,
            security_bits: self.security_bits,
            _phantom: PhantomData,
        });

        // Compute [n(s)] r, and verify that it matches p.
        let mul_result = self.curve_mul_endo::<InnerC>(CurveMulOp {
            scalar,
            point: result,
        });
        self.copy_curve(mul_result.mul_result, point);

        CurveMulEndoResult {
            mul_result: result,
            actual_scalar: mul_result.actual_scalar,
        }
    }

    /// Note: This assumes the most significant bit of each scalar is unset. This occurs with high
    /// probability if the field size is slightly larger than a power of two and the scalars are
    /// uniformly random.
    pub fn curve_msm<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &mut self,
        parts: &[CurveMulOp<C, InnerC>],
    ) -> AffinePointTarget<InnerC> {
        // We assume each most significant bit is unset; see the note in the method doc.
        let f_bits = C::ScalarField::BITS - 1;

        let all_bits: Vec<Vec<Target<C::ScalarField>>> = parts
            .iter()
            .map(|part| self.split_binary(part.scalar.convert(), f_bits))
            .collect();

        // Normally we would start with zero, but to avoid exceptional cases, we start with some
        // random nonzero point and subtract it later. This avoids exceptional cases with high
        // probability. A malicious prover may be able to craft an input which leads to an
        // exceptional case, but this isn't a problem as our curve gates will be unsatisfiable in
        // exceptional cases.
        let mut filler = blake_hash_base_field_to_curve::<InnerC>(InnerC::BaseField::ZERO);
        let mut acc = self.constant_affine_point(filler);
        let mut scalar_accs = vec![self.zero_wire(); parts.len()];

        for i in (0..f_bits).rev() {
            // Route the accumulator to the first curve addition gate's inputs.
            self.copy(
                acc.x,
                Target::Wire(Wire {
                    gate: self.num_gates(),
                    input: CurveAddGate::<C, InnerC>::WIRE_GROUP_ACC_X,
                }),
            );
            self.copy(
                acc.y,
                Target::Wire(Wire {
                    gate: self.num_gates(),
                    input: CurveAddGate::<C, InnerC>::WIRE_GROUP_ACC_Y,
                }),
            );

            for (j, part) in parts.iter().enumerate() {
                let bit = all_bits[j][i];

                let idx_add = self.num_gates();
                self.add_gate_no_constants(CurveAddGate::<C, InnerC>::new(idx_add));
                self.copy(
                    scalar_accs[j],
                    Target::Wire(Wire {
                        gate: idx_add,
                        input: CurveAddGate::<C, InnerC>::WIRE_SCALAR_ACC_OLD,
                    }),
                );
                scalar_accs[j] = Target::Wire(Wire {
                    gate: idx_add,
                    input: CurveAddGate::<C, InnerC>::WIRE_SCALAR_ACC_NEW,
                });
                self.copy(
                    part.point.x,
                    Target::Wire(Wire {
                        gate: idx_add,
                        input: CurveAddGate::<C, InnerC>::WIRE_ADDEND_X,
                    }),
                );
                self.copy(
                    part.point.y,
                    Target::Wire(Wire {
                        gate: idx_add,
                        input: CurveAddGate::<C, InnerC>::WIRE_ADDEND_Y,
                    }),
                );
                self.copy(
                    bit,
                    Target::Wire(Wire {
                        gate: idx_add,
                        input: CurveAddGate::<C, InnerC>::WIRE_SCALAR_BIT,
                    }),
                );
            }

            // Double the accumulator.
            let idx_dbl = self.num_gates();
            self.add_gate_no_constants(CurveDblGate::<C, InnerC>::new(idx_dbl));
            // No need to route the double gate's inputs, because the last add gate would have
            // constrained them.
            // Normally, we will take the double gate's outputs as the new accumulator. If we just
            // completed the last iteration of the MSM though, then we don't want to perform a final
            // doubling, so we will take its inputs as the result instead.
            if i == 0 {
                acc = AffinePointTarget {
                    x: Target::Wire(Wire {
                        gate: idx_dbl,
                        input: CurveDblGate::<C, InnerC>::WIRE_X_OLD,
                    }),
                    y: Target::Wire(Wire {
                        gate: idx_dbl,
                        input: CurveDblGate::<C, InnerC>::WIRE_Y_OLD,
                    }),
                };
            } else {
                acc = AffinePointTarget {
                    x: Target::Wire(Wire {
                        gate: idx_dbl,
                        input: CurveDblGate::<C, InnerC>::WIRE_X_NEW,
                    }),
                    y: Target::Wire(Wire {
                        gate: idx_dbl,
                        input: CurveDblGate::<C, InnerC>::WIRE_Y_NEW,
                    }),
                };

                // Also double the filler, so we can subtract out a rescaled version later.
                filler = filler.double();
            }
        }

        // Subtract (a rescaled version of) the arbitrary nonzero value that we started with.
        let filler_target = self.constant_affine_point(filler);
        acc = self.curve_sub::<InnerC>(acc, filler_target);

        // Assert that each accumulation of scalar bits matches the original scalar.
        for (j, part) in parts.iter().enumerate() {
            self.copy(scalar_accs[j], part.scalar.convert());
        }

        acc
    }

    /// Like `curve_msm`, but uses the endomorphism described in the Halo paper.
    pub fn curve_msm_endo<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &mut self,
        parts: &[CurveMulOp<C, InnerC>],
    ) -> CurveMsmEndoResult<C, InnerC> {
        let zero = self.zero_wire();

        // We assume each most significant bit is unset; see the note in curve_msm's method doc.
        let f_bits = C::ScalarField::BITS - 1;
        let scalar_bits = self.security_bits;
        let scalar_dibits = (f_bits - scalar_bits) / 2;

        // To keep things simple for now, we only handle the case of |F| ~= 2^254 and lambda = 128.
        assert_eq!(f_bits, 254);
        assert_eq!(scalar_bits, 128);
        assert_eq!(scalar_dibits, 63);

        // We split each scalar into 128 bits and 63 dibits. The bits are used in the MSM, while the
        // dibits are ignored, except that we need to include them in our sum-of-limbs computation
        // in order to verify that the decomposition was correct.
        let (all_bits, all_dibits): (
            Vec<Vec<Target<C::ScalarField>>>,
            Vec<Vec<Target<C::ScalarField>>>,
        ) = parts
            .iter()
            .map(|part| {
                self.split_binary_and_base_4(part.scalar.convert(), scalar_bits, scalar_dibits)
            })
            .unzip();

        // Normally we would start with zero, but to avoid exceptional cases, we start with some
        // random nonzero point and subtract it later. This avoids exceptional cases with high
        // probability. A malicious prover may be able to craft an input which leads to an
        // exceptional case, but this isn't a problem as our curve gates will be unsatisfiable in
        // exceptional cases.
        let mut filler = blake_hash_base_field_to_curve::<InnerC>(InnerC::BaseField::ZERO);
        let mut acc = self.constant_affine_point(filler);

        // For each scalar, we maintain two accumulators. The unsigned one is for computing a
        // weighted sum of bits and dibits in the usual manner, so that we can later check that this
        // sum equals the original scalar. The signed one is for computing n(s) for each scalar s.
        // This is the "actual" scalar by which the associated point was multiplied, accounting for
        // the endomorphism.
        let mut scalar_acc_unsigned = Vec::new();
        let mut scalar_acc_signed = Vec::new();

        // As in the Halo paper, we process two scalar bits at a time.
        for i in (0..scalar_bits).step_by(2).rev() {
            // Route the point accumulator to the first gate's inputs.
            self.copy(
                acc.x,
                Target::Wire(Wire {
                    gate: self.num_gates(),
                    input: CurveEndoGate::<C, InnerC>::WIRE_GROUP_ACC_X,
                }),
            );
            self.copy(
                acc.y,
                Target::Wire(Wire {
                    gate: self.num_gates(),
                    input: CurveEndoGate::<C, InnerC>::WIRE_GROUP_ACC_Y,
                }),
            );

            for (j, part) in parts.iter().enumerate() {
                let bit_0 = all_bits[j][i];
                let bit_1 = all_bits[j][i + 1];

                let gate = self.num_gates();
                self.add_gate_no_constants(CurveEndoGate::<C, InnerC>::new(gate));

                self.copy(
                    part.point.x,
                    Target::Wire(Wire {
                        gate,
                        input: CurveEndoGate::<C, InnerC>::WIRE_ADDEND_X,
                    }),
                );
                self.copy(
                    part.point.y,
                    Target::Wire(Wire {
                        gate,
                        input: CurveEndoGate::<C, InnerC>::WIRE_ADDEND_Y,
                    }),
                );
                self.copy(
                    bit_0,
                    Target::Wire(Wire {
                        gate,
                        input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_BIT_0,
                    }),
                );
                self.copy(
                    bit_1,
                    Target::Wire(Wire {
                        gate,
                        input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_BIT_1,
                    }),
                );

                // If this is the first pair of scalar bits being processed, route 0 to the scalar accumulators.
                if i == scalar_bits - 2 {
                    self.copy(
                        zero,
                        Target::Wire(Wire {
                            gate,
                            input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_ACC_UNSIGNED,
                        }),
                    );
                    self.copy(
                        zero,
                        Target::Wire(Wire {
                            gate,
                            input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_ACC_SIGNED,
                        }),
                    );
                }

                // If this is the last pair of scalar bits being processed, save the final
                // accumulator states.
                // Since CurveEndoGate will store these in the "next" gate, but this is the last
                // CurveEndoGate for this scalar, we need to add an extra BufferGate to receive them.
                if i == 0 {
                    let gate = self.num_gates();
                    self.add_gate_no_constants(BufferGate::new(gate));
                    scalar_acc_unsigned.push(Target::Wire(Wire {
                        gate,
                        input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_ACC_UNSIGNED,
                    }));
                    scalar_acc_signed.push(Target::Wire(Wire {
                        gate,
                        input: CurveEndoGate::<C, InnerC>::WIRE_SCALAR_ACC_SIGNED,
                    }));
                }
            }

            // Double the accumulator.
            let gate = self.num_gates();
            self.add_gate_no_constants(CurveDblGate::<C, InnerC>::new(gate));
            // No need to route the double gate's inputs, because the last endo gate would have
            // constrained them.
            // Normally, we will take the double gate's outputs as the new accumulator. If we just
            // completed the last iteration of the MSM though, then we don't want to perform a final
            // doubling, so we will take its inputs as the result instead.
            if i == 0 {
                acc = AffinePointTarget {
                    x: Target::Wire(Wire {
                        gate,
                        input: CurveDblGate::<C, InnerC>::WIRE_X_OLD,
                    }),
                    y: Target::Wire(Wire {
                        gate,
                        input: CurveDblGate::<C, InnerC>::WIRE_Y_OLD,
                    }),
                };
            } else {
                acc = AffinePointTarget {
                    x: Target::Wire(Wire {
                        gate,
                        input: CurveDblGate::<C, InnerC>::WIRE_X_NEW,
                    }),
                    y: Target::Wire(Wire {
                        gate,
                        input: CurveDblGate::<C, InnerC>::WIRE_Y_NEW,
                    }),
                };

                // Also double the filler, so we can subtract out a rescaled version later.
                filler = filler.double();
            }
        }

        // Subtract (a rescaled version of) the arbitrary nonzero value that we started with.
        let filler_target = self.constant_affine_point(filler);
        acc = self.curve_sub::<InnerC>(acc, filler_target);

        // By now we've accumulated all the bits of each scalar, but we also need to accumulate the dibits.
        for j in 0..parts.len() {
            for dibits_chunk in all_dibits[j].chunks(Base4SumGate::<C>::NUM_LIMBS) {
                assert_eq!(dibits_chunk.len(), Base4SumGate::<C>::NUM_LIMBS);

                let gate = self.num_gates();
                self.add_gate_no_constants(Base4SumGate::new(gate));
                self.copy(
                    scalar_acc_unsigned[j],
                    Target::Wire(Wire {
                        gate,
                        input: Base4SumGate::<C>::WIRE_ACC_OLD,
                    }),
                );
                scalar_acc_unsigned[j] = Target::Wire(Wire {
                    gate,
                    input: Base4SumGate::<C>::WIRE_ACC_NEW,
                });

                for (i, &dibit) in dibits_chunk.iter().enumerate() {
                    self.copy(
                        dibit,
                        Target::Wire(Wire {
                            gate,
                            input: Base4SumGate::<C>::wire_limb(i),
                        }),
                    );
                }
            }
        }

        // Finally, assert that each unsigned accumulator matches the original scalar.
        for (j, part) in parts.iter().enumerate() {
            self.copy(scalar_acc_unsigned[j], part.scalar.convert());
        }

        CurveMsmEndoResult {
            msm_result: acc,
            actual_scalars: scalar_acc_signed,
        }
    }
}

#[cfg(test)]
mod tests {
    use anyhow::Result;

    use crate::util::get_canonical_gates;
    use crate::{verify_proof, CircuitBuilder, Curve, CurveMulOp, Field, PartialWitness, Tweedledee, Tweedledum};

    #[test]
    // TODO: This fails because curve_mul_endo has a flaw.
    #[ignore]
    fn test_curve_mul_inv_endo() -> Result<()> {
        type C = Tweedledee;
        type InnerC = Tweedledum;
        type SF = <C as Curve>::ScalarField;
        let mut builder = CircuitBuilder::<Tweedledee>::new::<InnerC>(128);

        // This is an arbitrary nonzero scalar.
        let scalar = builder.constant_wire(SF::FIVE);

        // Let p1 be an arbitrary point.
        let p1 = builder.constant_affine_point::<InnerC>(InnerC::GENERATOR_AFFINE);

        // Let p2 = [n(s)] p1.
        let mul_result = builder.curve_mul_endo::<InnerC>(CurveMulOp { scalar, point: p1 });
        let p2 = mul_result.mul_result;

        // Let p3 = [1 / n(s)] p2.
        let mul_inv_result = builder.curve_mul_inv_endo::<InnerC>(CurveMulOp { scalar, point: p2 });
        let p3 = mul_inv_result.mul_result;

        // Since the scalars cancel, p1 = p3.
        builder.copy_curve(p1, p3);

        let circuit = builder.build();
        let witness = circuit.generate_witness(PartialWitness::new());
        let proof = circuit.generate_proof::<InnerC>(&witness, &[], false)?;
        let vk = circuit.to_vk();
        verify_proof::<C, InnerC>(
            &[],
            &proof,
            &[],
            &vk,
            true,
            get_canonical_gates::<C, InnerC>().into(),
        )?;

        Ok(())
    }
}
