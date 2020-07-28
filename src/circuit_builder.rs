use std::collections::{BTreeMap, HashMap};

use crate::gates::*;
use crate::plonk_util::{commit_polynomials, halo_n, polynomials_to_values_padded, sigma_polynomials, values_to_polynomials};
use crate::util::{ceil_div_usize, log2_strict, transpose};
use crate::{blake_hash_base_field_to_curve, blake_hash_usize_to_curve, fft_precompute, generate_rescue_constants, msm_precompute, AffinePoint, AffinePointTarget, Circuit, CurveMsmEndoResult, CurveMulEndoResult, CurveMulOp, Field, HaloCurve, PartialWitness, PublicInput, Target, TargetPartitions, VirtualTarget, Wire, WitnessGenerator, NUM_CONSTANTS, NUM_WIRES, BoundedTarget};
use std::marker::PhantomData;
use num::{Zero, BigUint};

pub struct CircuitBuilder<C: HaloCurve> {
    pub(crate) security_bits: usize,
    public_input_index: usize,
    virtual_target_index: usize,
    gate_counts: BTreeMap<&'static str, usize>,
    gate_constants: Vec<Vec<C::ScalarField>>,
    copy_constraints: Vec<(Target, Target)>,
    generators: Vec<Box<dyn WitnessGenerator<C::ScalarField>>>,
    constant_wires: HashMap<C::ScalarField, Target>,
}

impl<C: HaloCurve> CircuitBuilder<C> {
    pub fn new(security_bits: usize) -> Self {
        CircuitBuilder {
            security_bits,
            public_input_index: 0,
            virtual_target_index: 0,
            gate_counts: BTreeMap::new(),
            gate_constants: Vec::new(),
            copy_constraints: Vec::new(),
            generators: Vec::new(),
            constant_wires: HashMap::new(),
        }
    }

    pub fn stage_public_input(&mut self) -> PublicInput {
        let index = self.public_input_index;
        self.public_input_index += 1;
        PublicInput { index }
    }

    pub fn stage_public_inputs(&mut self, n: usize) -> Vec<PublicInput> {
        (0..n).map(|_i| self.stage_public_input()).collect()
    }

    /// Add `PublicInputGate`s which enable public inputs to be routed. Should be called after all
    /// `stage_public_input[s]` calls, but before any gates are added.
    pub fn route_public_inputs(&mut self) {
        debug_assert_eq!(
            self.num_gates(),
            0,
            "Must be called before any gates are added"
        );
        let num_pi_gates = ceil_div_usize(self.public_input_index, NUM_WIRES);
        for i in 0..num_pi_gates {
            self.add_gate_no_constants(PublicInputGate::new(i * 2));
            self.add_gate_no_constants(BufferGate::new(i * 2 + 1));
        }
    }

    pub fn add_virtual_target(&mut self) -> Target {
        let index = self.virtual_target_index;
        self.virtual_target_index += 1;
        Target::VirtualTarget(VirtualTarget { index })
    }

    pub fn add_virtual_targets(&mut self, n: usize) -> Vec<Target> {
        (0..n).map(|_i| self.add_virtual_target()).collect()
    }

    pub fn add_virtual_point_target<InnerC: HaloCurve>(&mut self) -> AffinePointTarget<InnerC> {
        let x = self.add_virtual_target();
        let y = self.add_virtual_target();
        AffinePointTarget::new(x, y)
    }

    pub fn add_virtual_point_targets<InnerC: HaloCurve>(
        &mut self,
        n: usize,
    ) -> Vec<AffinePointTarget<InnerC>> {
        (0..n).map(|_i| self.add_virtual_point_target()).collect()
    }

    pub fn zero_wire(&mut self) -> Target {
        self.constant_wire(C::ScalarField::ZERO)
    }

    pub fn zero_bounded_target(&mut self) -> BoundedTarget {
        let zero = self.zero_wire();
        BoundedTarget { target: zero, max: BigUint::zero() }
    }

    pub fn one_wire(&mut self) -> Target {
        self.constant_wire(C::ScalarField::ONE)
    }

    pub fn two_wire(&mut self) -> Target {
        self.constant_wire(C::ScalarField::TWO)
    }

    pub fn neg_one_wire(&mut self) -> Target {
        self.constant_wire(C::ScalarField::NEG_ONE)
    }

    pub fn constant_wire(&mut self, c: C::ScalarField) -> Target {
        if self.constant_wires.contains_key(&c) {
            self.constant_wires[&c]
        } else {
            let result = self.create_constant_wire(c);
            self.constant_wires.insert(c, result);
            result
        }
    }

    pub fn constant_wires(&mut self, constants: &[C::ScalarField]) -> Vec<Target> {
        constants.iter().map(|&c| self.constant_wire(c)).collect()
    }

    pub fn constant_wire_u32(&mut self, c: u32) -> Target {
        self.constant_wire(C::ScalarField::from_canonical_u32(c))
    }

    fn create_constant_wire(&mut self, c: C::ScalarField) -> Target {
        // We will create a ConstantGate and pass c as its first (and only) constant, which will
        // cause it to populate its output wire with the same value c.
        let gate = self.num_gates();
        self.add_gate(ConstantGate::new(gate), vec![c]);
        Target::Wire(Wire {
            gate,
            input: ConstantGate::<C>::WIRE_OUTPUT,
        })
    }

    pub fn constant_affine_point<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &mut self,
        point: AffinePoint<InnerC>,
    ) -> AffinePointTarget<InnerC> {
        assert!(!point.zero);
        AffinePointTarget::new(
            self.constant_wire(point.x),
            self.constant_wire(point.y),
        )
    }

    /// Adds a generator to populate the given target with the given constant.
    pub fn generate_constant(&mut self, target: Target, c: C::ScalarField) {
        struct ConstantGenerator<F: Field> {
            target: Target,
            c: F,
        }

        impl<F: Field> WitnessGenerator<F> for ConstantGenerator<F> {
            fn dependencies(&self) -> Vec<Target> {
                Vec::new()
            }

            fn generate(
                &self,
                _constants: &Vec<Vec<F>>,
                _witness: &PartialWitness<F>,
            ) -> PartialWitness<F> {
                let mut result = PartialWitness::new();
                result.set_target(self.target, self.c);
                result
            }
        }

        self.add_generator(ConstantGenerator { target, c });
    }

    pub fn assert_zero(&mut self, x: Target) {
        let zero = self.zero_wire();
        self.copy(x, zero);
    }

    pub fn assert_one(&mut self, x: Target) {
        let one = self.one_wire();
        self.copy(x, one);
    }

    pub fn assert_binary(&mut self, x: Target) {
        // This is typically implemented with a constraint like x * (x - 1) = 0.
        // We rewrite this as x * x - x = 0, which requires just one gate in our model.
        let lhs = self.mul_sub(x, x, x);
        self.assert_zero(lhs);
    }

    /// Assert that each of the given targets is less than 4.
    pub fn assert_all_base_4(&mut self, limbs: &[Target]) {
        // We will leverage Base4SumGate, which checks that each of its limbs is base 4.
        for chunk in limbs.chunks(Base4SumGate::<C>::NUM_ROUTED_LIMBS) {
            let gate = self.num_gates();
            self.add_gate_no_constants(Base4SumGate::new(gate));

            // We don't care about Base4SumGate's accumulator wires, but we need to pass some
            // (arbitrary) value to the old accumulator wire in order for the generator to run.
            self.generate_constant(
                Target::Wire(Wire { gate, input: Base4SumGate::<C>::WIRE_ACC_OLD }),
                C::ScalarField::ZERO, // This value is arbitrary.
            );

            // Route each limb to one of Base4SumGate's routed limb wires.
            for (i, &limb) in chunk.iter().enumerate() {
                self.copy(
                    limb,
                    Target::Wire(Wire {
                        gate,
                        input: Base4SumGate::<C>::wire_limb(i),
                    })
                )
            }
        }
    }

    pub fn assert_nonzero(&mut self, x: Target) {
        // An element is nonzero iff it has an inverse.
        self.inv(x);
    }

    /// Returns `if x == 0 { 1 } else { 0 }`.
    pub fn is_zero(&mut self, x: Target) -> Target {
        // This is similar to the technique described in
        // https://github.com/mir-protocol/r1cs-workshop/blob/master/workshop.pdf

        let is_zero: Target = self.add_virtual_target();

        // m will hold if x != 0 { -1 / x } else { 1 }.
        let m = self.add_virtual_target();

        struct IsZeroGenerator {
            x: Target,
            m: Target,
            is_zero: Target,
        }

        impl<F: Field> WitnessGenerator<F> for IsZeroGenerator {
            fn dependencies(&self) -> Vec<Target> {
                vec![self.x]
            }

            fn generate(
                &self,
                _constants: &Vec<Vec<F>>,
                witness: &PartialWitness<F>,
            ) -> PartialWitness<F> {
                let x = witness.get_target(self.x);
                let (m, is_zero) = if x.is_zero() {
                    (F::ONE, F::ONE)
                } else {
                    (-x.multiplicative_inverse().unwrap(), F::ZERO)
                };

                let mut result = PartialWitness::new();
                result.set_target(self.m, m);
                result.set_target(self.is_zero, is_zero);
                result
            }
        }

        self.add_generator(IsZeroGenerator { x, m, is_zero });

        // Enforce that is_zero = x * m + 1.
        let one = self.one_wire();
        let x_m_plus_1 = self.mul_add(x, m, one);
        self.copy(is_zero, x_m_plus_1);

        // Enforce that is_zero * x = 0.
        let is_zero_x = self.mul(is_zero, x);
        self.assert_zero(is_zero_x);

        is_zero
    }

    /// Returns `if x != 0 { 1 } else { 0 }`.
    pub fn is_nonzero(&mut self, x: Target) -> Target {
        let one = self.one_wire();
        let is_zero = self.is_zero(x);
        self.sub(one, is_zero)
    }

    /// Returns `if x == y { 1 } else { 0 }`.
    pub fn is_equal(&mut self, x: Target, y: Target) -> Target {
        let diff = self.sub(x, y);
        self.is_zero(diff)
    }

    /// Returns `if x != y { 1 } else { 0 }`.
    pub fn is_not_equal(&mut self, x: Target, y: Target) -> Target {
        let diff = self.sub(x, y);
        self.is_nonzero(diff)
    }

    /// Selects `x` or `y` based on `b`, which is assumed to be binary.
    /// In particular, this returns `if b { x } else { y }`.
    pub fn select(&mut self, b: Target, x: Target, y: Target) -> Target {
        // This can be computed various ways, e.g.
        //     b x + (1 - b) y
        //     b x + y - b y
        //     y + b (x - y)
        // We will actually compute it as
        //     b x - (b y - y)
        // since this can be done with two mul_sub calls.

        let b_y_minus_y = self.mul_sub(b, y, y);
        self.mul_sub(b, x, b_y_minus_y)
    }

    /// Returns the negation of a bit `b`, which is assumed to be in `{0, 1}`.
    pub fn not(&mut self, b: Target) -> Target {
        let one = self.one_wire();
        self.sub(one, b)
    }

    pub fn add(&mut self, x: Target, y: Target) -> Target {
        let zero = self.zero_wire();
        if x == zero {
            return y;
        }
        if y == zero {
            return x;
        }

        let one = self.one_wire();
        let index = self.num_gates();
        self.add_gate(
            ArithmeticGate::new(index),
            vec![C::ScalarField::ONE, C::ScalarField::ONE],
        );
        self.copy(
            x,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_0,
            }),
        );
        self.copy(
            one,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_1,
            }),
        );
        self.copy(
            y,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_ADDEND,
            }),
        );
        Target::Wire(Wire {
            gate: index,
            input: ArithmeticGate::<C>::WIRE_OUTPUT,
        })
    }

    pub fn add_many(&mut self, terms: &[Target]) -> Target {
        let mut sum = self.zero_wire();
        for term in terms {
            sum = self.add(sum, *term);
        }
        sum
    }

    pub fn double(&mut self, x: Target) -> Target {
        self.add(x, x)
    }

    pub fn sub(&mut self, x: Target, y: Target) -> Target {
        let zero = self.zero_wire();
        if y == zero {
            return x;
        }

        let one = self.one_wire();
        let index = self.num_gates();
        self.add_gate(
            ArithmeticGate::new(index),
            vec![C::ScalarField::ONE, C::ScalarField::NEG_ONE],
        );
        self.copy(
            x,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_0,
            }),
        );
        self.copy(
            one,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_1,
            }),
        );
        self.copy(
            y,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_ADDEND,
            }),
        );
        Target::Wire(Wire {
            gate: index,
            input: ArithmeticGate::<C>::WIRE_OUTPUT,
        })
    }

    pub fn mul(&mut self, x: Target, y: Target) -> Target {
        let one = self.one_wire();
        if x == one {
            return y;
        }
        if y == one {
            return x;
        }

        let zero = self.zero_wire();
        let index = self.num_gates();
        self.add_gate(
            ArithmeticGate::new(index),
            vec![C::ScalarField::ONE, C::ScalarField::ZERO],
        );
        self.copy(
            x,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_0,
            }),
        );
        self.copy(
            y,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_1,
            }),
        );
        self.copy(
            zero,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_ADDEND,
            }),
        );
        Target::Wire(Wire {
            gate: index,
            input: ArithmeticGate::<C>::WIRE_OUTPUT,
        })
    }

    pub fn mul_many(&mut self, terms: &[Target]) -> Target {
        let mut product = self.one_wire();
        for term in terms {
            product = self.mul(product, *term);
        }
        product
    }

    pub fn square(&mut self, x: Target) -> Target {
        self.mul(x, x)
    }

    /// Note: This assumes the most significant bit of each scalar is unset. This occurs with high
    /// probability if the field size is slightly larger than a power of two and the inputs are
    /// uniformly random.
    pub fn deterministic_square_root(&mut self, x: Target) -> Target {
        // Assume x != 0. Let y, z be the square roots of x. Since y + z = |F|, and |F| is odd, the
        // parity of y and z must differ, so we can enforce determinism by checking for a certain
        // parity bit state. We chose a parity bit of 0 since this also works for the x = 0 case.

        struct SqrtGenerator {
            x: Target,
            x_sqrt: Target,
        }

        impl<F: Field> WitnessGenerator<F> for SqrtGenerator {
            fn dependencies(&self) -> Vec<Target> {
                vec![self.x]
            }

            fn generate(
                &self,
                _constants: &Vec<Vec<F>>,
                witness: &PartialWitness<F>,
            ) -> PartialWitness<F> {
                let x_value = witness.get_target(self.x);
                let mut x_sqrt_value = x_value.square_root().expect("Not square");

                if x_sqrt_value.to_canonical_bool_vec()[0] {
                    // The parity bit is 1; we want the other square root.
                    x_sqrt_value = -x_sqrt_value;
                    debug_assert!(!x_sqrt_value.to_canonical_bool_vec()[0]);
                }

                let mut result = PartialWitness::new();
                result.set_target(self.x_sqrt, x_sqrt_value);
                result
            }
        }

        let x_sqrt = self.add_virtual_target();
        self.add_generator(SqrtGenerator { x, x_sqrt });

        // We assume each most significant bit is unset; see the note in the method doc.
        let f_bits = C::ScalarField::BITS - 1;
        assert_eq!(
            f_bits, 254,
            "We currently only handle fields of size 2^254 + epsilon"
        );
        let (bits, dibits) = self.split_binary_and_base_4(x_sqrt, 2, 126);

        // Verify that x_sqrt * x_sqrt = x.
        let x_sqrt_squared = self.square(x_sqrt);
        self.copy(x_sqrt_squared, x);

        // Verify that the parity bit is 0, and the other bit is binary.
        self.assert_zero(bits[0]);
        self.assert_binary(bits[1]);

        // Verify the decomposition by taking a weighted sum of the limbs. Since the first bit is
        // always zero, we start with the second bit (scaled by its weight of 2) and then add the
        // base 4 limbs.
        let mut sum = self.double(bits[1]);
        for chunk in dibits.chunks(Base4SumGate::<C>::NUM_LIMBS) {
            assert_eq!(
                chunk.len(),
                Base4SumGate::<C>::NUM_LIMBS,
                "Should not have a partial chunk"
            );

            let index = self.num_gates();
            self.add_gate_no_constants(Base4SumGate::new(index));
            for i in 0..chunk.len() {
                self.copy(
                    sum,
                    Target::Wire(Wire {
                        gate: index,
                        input: Base4SumGate::<C>::WIRE_ACC_OLD,
                    }),
                );
                self.copy(
                    chunk[i],
                    Target::Wire(Wire {
                        gate: index,
                        input: Base4SumGate::<C>::wire_limb(i),
                    }),
                );
                sum = Target::Wire(Wire {
                    gate: index,
                    input: Base4SumGate::<C>::WIRE_ACC_NEW,
                })
            }
        }
        self.copy(sum, x);

        x_sqrt
    }

    /// Compute `x^power`, where `power` is a constant.
    pub fn exp_constant(&mut self, x: Target, power: C::ScalarField) -> Target {
        let power_bits = power.num_bits();
        let mut current = x;
        let mut product = self.one_wire();

        for (i, limb) in power.to_canonical_u64_vec().iter().enumerate() {
            for j in 0..64 {
                // If we've gone through all the 1 bits already, no need to keep squaring.
                let bit_index = i * 64 + j;
                if bit_index == power_bits {
                    return product;
                }

                if (limb >> j & 1) != 0 {
                    product = self.mul(product, current);
                }
                current = self.square(current);
            }
        }

        product
    }

    /// Compute `x^power`, where `power` is a constant `usize`.
    pub fn exp_constant_usize(&mut self, x: Target, power: usize) -> Target {
        self.exp_constant(x, C::ScalarField::from_canonical_usize(power))
    }

    pub fn inv(&mut self, x: Target) -> Target {
        struct InverseGenerator {
            x: Target,
            x_inv: Target,
        }

        impl<F: Field> WitnessGenerator<F> for InverseGenerator {
            fn dependencies(&self) -> Vec<Target> {
                vec![self.x]
            }

            fn generate(
                &self,
                _constants: &Vec<Vec<F>>,
                witness: &PartialWitness<F>,
            ) -> PartialWitness<F> {
                let x_value = witness.get_target(self.x);
                let x_inv_value = x_value.multiplicative_inverse().expect("x = 0");

                let mut result = PartialWitness::new();
                result.set_target(self.x_inv, x_inv_value);
                result
            }
        }

        let x_inv = self.add_virtual_target();
        self.add_generator(InverseGenerator { x, x_inv });

        // Enforce that x * x_inv = 1.
        let product = self.mul(x, x_inv);
        let one = self.one_wire();
        self.copy(product, one);

        x_inv
    }

    pub fn div(&mut self, x: Target, y: Target) -> Target {
        let y_inv = self.inv(y);
        self.mul(x, y_inv)
    }

    /// Multiply and add; i.e. computes `x * y + z`.
    pub fn mul_add(&mut self, x: Target, y: Target, z: Target) -> Target {
        let index = self.num_gates();
        self.add_gate(
            ArithmeticGate::new(index),
            vec![C::ScalarField::ONE, C::ScalarField::ONE],
        );
        self.copy(
            x,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_0,
            }),
        );
        self.copy(
            y,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_1,
            }),
        );
        self.copy(
            z,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_ADDEND,
            }),
        );
        Target::Wire(Wire {
            gate: index,
            input: ArithmeticGate::<C>::WIRE_OUTPUT,
        })
    }

    pub(crate) fn bounded_mul_add(
        &mut self,
        x: &BoundedTarget,
        y: &BoundedTarget,
        z: &BoundedTarget,
    ) -> BoundedTarget {
        let target = self.mul_add(x.target, y.target, z.target);
        let max = &x.max * &y.max + &z.max;
        BoundedTarget { target, max }
    }

    /// Multiply and subtract; i.e. computes `x * y - z`.
    pub fn mul_sub(&mut self, x: Target, y: Target, z: Target) -> Target {
        let index = self.num_gates();
        self.add_gate(
            ArithmeticGate::new(index),
            vec![C::ScalarField::ONE, C::ScalarField::NEG_ONE],
        );
        self.copy(
            x,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_0,
            }),
        );
        self.copy(
            y,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_MULTIPLICAND_1,
            }),
        );
        self.copy(
            z,
            Target::Wire(Wire {
                gate: index,
                input: ArithmeticGate::<C>::WIRE_ADDEND,
            }),
        );
        Target::Wire(Wire {
            gate: index,
            input: ArithmeticGate::<C>::WIRE_OUTPUT,
        })
    }

    /// Computes `-x`.
    pub fn neg(&mut self, x: Target) -> Target {
        let neg_one = self.neg_one_wire();
        self.mul(x, neg_one)
    }

    /// Splits `x` into its binary representation. Note that this method merely adds a generator to
    /// populate the bit wires; it does not enforce constraints to verify the decomposition.
    fn split_binary(&mut self, x: Target, num_bits: usize) -> Vec<Target> {
        let (bits, _dibits) = self.split_binary_and_base_4(x, num_bits, 0);
        bits
    }
    
    /// Splits `x` into its base 4 representation. Note that this method merely adds a generator to
    /// populate the bit wires; it does not enforce constraints to verify the decomposition.
    fn split_base_4(&mut self, x: Target, num_dibits: usize) -> Vec<Target> {
        let (_bits, dibits) = self.split_binary_and_base_4(x, 0, num_dibits);
        dibits
    }

    /// Splits `x` into a combination of binary and base 4 limbs. Note that this method merely adds
    /// a generator to populate the limb wires; it does not enforce constraints to verify the
    /// decomposition.
    fn split_binary_and_base_4(
        &mut self,
        x: Target,
        num_bits: usize,
        num_dibits: usize,
    ) -> (Vec<Target>, Vec<Target>) {
        struct SplitGenerator {
            x: Target,
            bits: Vec<Target>,
            dibits: Vec<Target>,
        }

        impl<F: Field> WitnessGenerator<F> for SplitGenerator {
            fn dependencies(&self) -> Vec<Target> {
                vec![self.x]
            }

            fn generate(
                &self,
                _constants: &Vec<Vec<F>>,
                witness: &PartialWitness<F>,
            ) -> PartialWitness<F> {
                let x = witness.get_target(self.x);
                let x_bits = x.to_canonical_bool_vec();

                let mut result = PartialWitness::new();
                for i in 0..self.bits.len() {
                    result.set_target(self.bits[i], F::from_canonical_bool(x_bits[i]));
                }
                for i in 0..self.dibits.len() {
                    let bit_1 = x_bits[self.bits.len() + i * 2];
                    let bit_2 = x_bits[self.bits.len() + i * 2 + 1];
                    let dibit = if bit_1 { 1 } else { 0 } + if bit_2 { 2 } else { 0 };
                    result.set_target(self.dibits[i], F::from_canonical_u32(dibit));
                }
                result
            }
        }

        let bits = self.add_virtual_targets(num_bits);
        let dibits = self.add_virtual_targets(num_dibits);
        let generator = SplitGenerator {
            x,
            bits: bits.clone(),
            dibits: dibits.clone(),
        };
        self.add_generator(generator);
        (bits, dibits)
    }

    /// Asserts that the given target's value is small enough to fit in the given number of dibits.
    ///
    /// Note: This is most efficient when `num_dibits` is a multiple of `Base4SumGate::NUM_LIMBS`.
    pub(crate) fn assert_dibit_length(&mut self, x: Target, num_dibits: usize) {
        // Get the purported base 4 decomposition of x.
        let dibits = self.split_base_4(x, num_dibits);

        // Accumulate each full chunk of NUM_LIMBS dibits using a Base4SumGate.
        let mut sum = self.zero_wire();
        let chunks = dibits.chunks_exact(Base4SumGate::<C>::NUM_LIMBS);
        for chunk in chunks.clone() {
            let gate = self.num_gates();
            self.add_gate_no_constants(Base4SumGate::new(gate));

            // Route sum into WIRE_ACC_OLD.
            self.copy(
                sum,
                Target::Wire(Wire {
                    gate,
                    input: Base4SumGate::<C>::WIRE_ACC_OLD,
                }),
            );

            for (i, &dibit) in chunk.iter().enumerate() {
                self.copy(
                    dibit,
                    Target::Wire(Wire {
                        gate,
                        input: Base4SumGate::<C>::wire_limb(i),
                    })
                );
            }

            // Take WIRE_ACC_NEW as our updated sum.
            sum = Target::Wire(Wire {
                gate,
                input: Base4SumGate::<C>::WIRE_ACC_NEW,
            });
        }

        // If there is a partial chunk of dibits, it would be difficult to accumulate it with
        // Base4SumGate, e.g. since it would always perform NUM_LIMBS doublings. So instead, we
        // will accumulate it with simple arithmetic operations.
        if !chunks.remainder().is_empty() {
            self.assert_all_base_4(chunks.remainder());
            let four = self.constant_wire_u32(4);
            for &dibit in chunks.remainder() {
                sum = self.mul_add(sum, four, dibit);
            }
        }

        // Ensure that the weighted sum we computed matches the original value, x.
        self.copy(sum, x);
    }

    pub fn rescue_hash_n_to_1(&mut self, inputs: &[Target]) -> Target {
        self.rescue_sponge(inputs, 1)[0]
    }

    pub fn rescue_hash_n_to_2(&mut self, inputs: &[Target]) -> (Target, Target) {
        let outputs = self.rescue_sponge(inputs, 2);
        (outputs[0], outputs[1])
    }

    pub fn rescue_hash_n_to_3(&mut self, inputs: &[Target]) -> (Target, Target, Target) {
        let outputs = self.rescue_sponge(inputs, 3);
        (outputs[0], outputs[1], outputs[2])
    }

    pub fn rescue_sponge(&mut self, inputs: &[Target], num_outputs: usize) -> Vec<Target> {
        let zero = self.zero_wire();
        let mut state = vec![zero; RESCUE_SPONGE_WIDTH];

        // Absorb all input chunks.
        for input_chunk in inputs.chunks(RESCUE_SPONGE_WIDTH - 1) {
            for i in 0..input_chunk.len() {
                state[i] = self.add(state[i], input_chunk[i]);
            }
            state = self.rescue_permutation(&state);
        }

        // Squeeze until we have the desired number of outputs.
        let mut outputs = Vec::new();
        loop {
            for i in 0..(RESCUE_SPONGE_WIDTH - 1) {
                outputs.push(state[i]);
                if outputs.len() == num_outputs {
                    return outputs;
                }
            }
            state = self.rescue_permutation(&state);
        }
    }

    pub fn rescue_permutation(&mut self, inputs: &[Target]) -> Vec<Target> {
        assert_eq!(inputs.len(), RESCUE_SPONGE_WIDTH);

        // Route the input wires.
        for i in 0..RESCUE_SPONGE_WIDTH {
            self.copy(
                inputs[i],
                Target::Wire(Wire {
                    gate: self.num_gates(),
                    input: RescueStepAGate::<C>::wire_acc(i),
                }),
            );
        }

        let all_constants = generate_rescue_constants(RESCUE_SPONGE_WIDTH, self.security_bits);
        for (a_constants, b_constants) in all_constants.into_iter() {
            let a_index = self.num_gates();
            let a_gate = RescueStepAGate::new(a_index);
            self.add_gate(a_gate, a_constants);

            let b_index = self.num_gates();
            let b_gate = RescueStepBGate::new(b_index);
            self.add_gate(b_gate, b_constants);
        }

        // Use a BufferGate to receive the final accumulator states.
        let gate = self.num_gates();
        self.add_gate_no_constants(BufferGate::new(gate));
        (0..RESCUE_SPONGE_WIDTH)
            .map(|i| {
                Target::Wire(Wire {
                    gate,
                    input: RescueStepBGate::<C>::wire_acc(i),
                })
            })
            .collect()
    }

    /// Assert that a given coordinate pair is on the curve `C`.
    pub fn curve_assert_valid<InnerC: HaloCurve<BaseField = C::ScalarField>>(
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

    pub fn curve_neg<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &mut self,
        p: AffinePointTarget<InnerC>,
    ) -> AffinePointTarget<InnerC> {
        let neg_y = self.neg(p.y);
        AffinePointTarget::new(p.x, neg_y)
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

        AffinePointTarget::new(result_x, result_y)
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
        AffinePointTarget::new(
            Target::Wire(Wire {
                gate: idx_dbl,
                input: CurveDblGate::<C, InnerC>::WIRE_X_NEW,
            }),
            Target::Wire(Wire {
                gate: idx_dbl,
                input: CurveDblGate::<C, InnerC>::WIRE_Y_NEW,
            })
        )
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
        mul: CurveMulOp<InnerC>,
    ) -> AffinePointTarget<InnerC> {
        self.curve_msm::<InnerC>(&[mul])
    }

    /// Computes `[n(s)] p`.
    pub fn curve_mul_endo<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &mut self,
        mul: CurveMulOp<InnerC>,
    ) -> CurveMulEndoResult<InnerC> {
        let result = self.curve_msm_endo::<InnerC>(&[mul]);
        CurveMulEndoResult {
            mul_result: result.msm_result,
            actual_scalar: result.actual_scalars[0],
            phantom: PhantomData,
        }
    }

    /// Computes `[1 / n(s)] p`.
    pub fn curve_mul_inv_endo<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &mut self,
        mul: CurveMulOp<InnerC>,
    ) -> CurveMulEndoResult<InnerC> {
        // We witness r = [1 / n(s)] p, then verify that [n(s)] r = p, and return r.
        let CurveMulOp { scalar, point, .. } = mul;

        struct ResultGenerator<InnerC: HaloCurve> {
            mul: CurveMulOp<InnerC>,
            result: AffinePointTarget<InnerC>,
            security_bits: usize,
            _phantom: PhantomData<InnerC>,
        }

        impl<InnerC: HaloCurve> WitnessGenerator<InnerC::BaseField> for ResultGenerator<InnerC> {
            fn dependencies(&self) -> Vec<Target> {
                vec![self.mul.scalar, self.mul.point.x, self.mul.point.y]
            }

            fn generate(
                &self,
                _constants: &Vec<Vec<InnerC::BaseField>>,
                witness: &PartialWitness<InnerC::BaseField>,
            ) -> PartialWitness<InnerC::BaseField> {
                let scalar = witness
                    .get_target(self.mul.scalar)
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
        self.add_generator(ResultGenerator::<InnerC> {
            mul,
            result,
            security_bits: self.security_bits,
            _phantom: PhantomData,
        });

        // Compute [n(s)] r, and verify that it matches p.
        let mul_result = self.curve_mul_endo::<InnerC>(
            CurveMulOp::new(scalar, result));
        self.copy_curve(mul_result.mul_result, point);

        CurveMulEndoResult {
            mul_result: result,
            actual_scalar: mul_result.actual_scalar,
            phantom: PhantomData,
        }
    }

    /// Note: This assumes the most significant bit of each scalar is unset. This occurs with high
    /// probability if the field size is slightly larger than a power of two and the scalars are
    /// uniformly random.
    pub fn curve_msm<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &mut self,
        parts: &[CurveMulOp<InnerC>],
    ) -> AffinePointTarget<InnerC> {
        // We assume each most significant bit is unset; see the note in the method doc.
        let f_bits = C::ScalarField::BITS - 1;

        let all_bits: Vec<Vec<Target>> = parts
            .iter()
            .map(|part| self.split_binary(part.scalar, f_bits))
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
                acc = AffinePointTarget::new(
                    Target::Wire(Wire {
                        gate: idx_dbl,
                        input: CurveDblGate::<C, InnerC>::WIRE_X_OLD,
                    }),
                    Target::Wire(Wire {
                        gate: idx_dbl,
                        input: CurveDblGate::<C, InnerC>::WIRE_Y_OLD,
                    }),
                );
            } else {
                acc = AffinePointTarget::new(
                    Target::Wire(Wire {
                        gate: idx_dbl,
                        input: CurveDblGate::<C, InnerC>::WIRE_X_NEW,
                    }),
                    Target::Wire(Wire {
                        gate: idx_dbl,
                        input: CurveDblGate::<C, InnerC>::WIRE_Y_NEW,
                    }),
                );

                // Also double the filler, so we can subtract out a rescaled version later.
                filler = filler.double();
            }
        }

        // Subtract (a rescaled version of) the arbitrary nonzero value that we started with.
        let filler_target = self.constant_affine_point(filler);
        acc = self.curve_sub::<InnerC>(acc, filler_target);

        // Assert that each accumulation of scalar bits matches the original scalar.
        for (j, part) in parts.iter().enumerate() {
            self.copy(scalar_accs[j], part.scalar);
        }

        acc
    }

    /// Like `curve_msm`, but uses the endomorphism described in the Halo paper.
    pub fn curve_msm_endo<InnerC: HaloCurve<BaseField = C::ScalarField>>(
        &mut self,
        parts: &[CurveMulOp<InnerC>],
    ) -> CurveMsmEndoResult<InnerC> {
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
        let (all_bits, all_dibits): (Vec<Vec<Target>>, Vec<Vec<Target>>) = parts
            .iter()
            .map(|part| self.split_binary_and_base_4(part.scalar, scalar_bits, scalar_dibits))
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
                acc = AffinePointTarget::new(
                    Target::Wire(Wire {
                        gate,
                        input: CurveDblGate::<C, InnerC>::WIRE_X_OLD,
                    }),
                    Target::Wire(Wire {
                        gate,
                        input: CurveDblGate::<C, InnerC>::WIRE_Y_OLD,
                    }),
                );
            } else {
                acc = AffinePointTarget::new(
                    Target::Wire(Wire {
                        gate,
                        input: CurveDblGate::<C, InnerC>::WIRE_X_NEW,
                    }),
                    Target::Wire(Wire {
                        gate,
                        input: CurveDblGate::<C, InnerC>::WIRE_Y_NEW,
                    }),
                );

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
            self.copy(scalar_acc_unsigned[j], part.scalar);
        }

        CurveMsmEndoResult {
            msm_result: acc,
            actual_scalars: scalar_acc_signed,
            phantom: PhantomData,
        }
    }

    /// Adds a gate to the circuit, without doing any routing.
    pub fn add_gate_no_constants<G: Gate<C>>(&mut self, gate: G) {
        self.add_gate(gate, Vec::new());
    }

    /// Adds a gate to the circuit, without doing any routing.
    pub fn add_gate<G: Gate<C>>(&mut self, gate: G, gate_constants: Vec<C::ScalarField>) {
        debug_assert!(G::PREFIX.len() + gate_constants.len() <= NUM_CONSTANTS);

        // Merge the gate type's prefix bits with the given gate config constants.
        let mut all_constants = Vec::new();
        for &prefix_bit in G::PREFIX {
            all_constants.push(C::ScalarField::from_canonical_bool(prefix_bit));
        }
        all_constants.extend(gate_constants);

        // Pad if not all constants were used.
        while all_constants.len() < NUM_CONSTANTS {
            all_constants.push(C::ScalarField::ZERO);
        }

        self.gate_constants.push(all_constants);
        self.add_generator(gate);
        *self.gate_counts.entry(G::NAME).or_insert(0) += 1;
    }

    pub fn add_generator<G: WitnessGenerator<C::ScalarField>>(&mut self, generator: G) {
        self.generators.push(Box::new(generator));
    }

    pub fn num_gates(&self) -> usize {
        self.gate_constants.len()
    }

    /// Add a copy constraint between two routing targets.
    pub fn copy(&mut self, target_1: Target, target_2: Target) {
        self.copy_constraints.push((target_1, target_2));
    }

    /// Add a copy constraint between two affine targets.
    pub fn copy_curve<InnerC: HaloCurve>(
        &mut self,
        affine_target_1: AffinePointTarget<InnerC>,
        affine_target_2: AffinePointTarget<InnerC>,
    ) {
        self.copy_constraints
            .push((affine_target_1.x, affine_target_2.x));
        self.copy_constraints
            .push((affine_target_1.y, affine_target_2.y));
    }

    /// Enforces a copy constraint between the two targets if the condition is non-zero.
    pub fn conditional_copy(&mut self, condition: Target, target_1: Target, target_2: Target) {
        let conditional_target_1 = self.mul(condition, target_1);
        let conditional_target_2 = self.mul(condition, target_2);
        self.copy(conditional_target_1, conditional_target_2);
    }

    /// Adds a gate with random wire values. By adding `k` of these gates, we can ensure that
    /// nothing is learned by opening the wire polynomials at `k` points outside of H.
    fn add_blinding_gate(&mut self) {
        let gate = self.num_gates();
        self.add_gate_no_constants(BufferGate::new(gate));
        for input in 0..NUM_WIRES {
            self.add_generator(RandomGenerator {
                target: Target::Wire(Wire { gate, input }),
            });
        }

        struct RandomGenerator {
            target: Target,
        }

        impl<F: Field> WitnessGenerator<F> for RandomGenerator {
            fn dependencies(&self) -> Vec<Target> {
                vec![]
            }

            fn generate(
                &self,
                _constants: &Vec<Vec<F>>,
                _witness: &PartialWitness<F>,
            ) -> PartialWitness<F> {
                let mut result = PartialWitness::new();
                result.set_target(self.target, F::rand());
                result
            }
        }
    }

    pub fn build(mut self) -> Circuit<C> {
        // Since we will open each polynomial at three points outside of H, we need three random
        // values to ensure nothing is learned from the out-of-H openings.
        for _i in 0..3 {
            self.add_blinding_gate();
        }

        // Print gate counts.
        println!("Gate counts:");
        for (gate, count) in &self.gate_counts {
            println!("{}: {}", gate, count);
        }
        println!();

        // Pad to a power of two.
        println!("Total gates before padding: {}", self.num_gates());
        while !self.num_gates().is_power_of_two() {
            // Add an empty gate.
            self.add_gate_no_constants(BufferGate::new(self.num_gates()));
        }
        println!("Total gates after padding: {}", self.num_gates());

        let degree = self.num_gates();
        let degree_pow = log2_strict(degree);
        let routing_target_partitions = self.get_routing_partitions();
        let wire_partitions = routing_target_partitions.to_wire_partitions();
        let sigma = wire_partitions.to_sigma();

        let CircuitBuilder {
            security_bits,
            public_input_index: num_public_inputs,
            gate_constants,
            generators,
            ..
        } = self;

        let fft_precomputation_n = fft_precompute(degree);
        let fft_precomputation_8n = fft_precompute(degree * 8);

        let subgroup_generator_n = C::ScalarField::primitive_root_of_unity(degree_pow);
        let subgroup_generator_8n = C::ScalarField::primitive_root_of_unity(degree_pow + 3);
        let subgroup_n = C::ScalarField::cyclic_subgroup_known_order(subgroup_generator_n, degree);
        let subgroup_8n =
            C::ScalarField::cyclic_subgroup_known_order(subgroup_generator_8n, 8 * degree);

        let pedersen_g: Vec<_> = (0..degree).map(blake_hash_usize_to_curve::<C>).collect();
        let pedersen_h = blake_hash_usize_to_curve::<C>(degree);
        let u = blake_hash_usize_to_curve::<C>(degree + 1);

        let w = 11; // TODO: Should really be set dynamically based on MSM size.
        let pedersen_g_msm_precomputation =
            msm_precompute(&AffinePoint::batch_to_projective(&pedersen_g), w);

        // While gate_constants is indexed by gate index first, this is indexed by wire index first.
        let wire_constants = transpose::<C::ScalarField>(&gate_constants);

        let constant_polynomials = values_to_polynomials(&wire_constants, &fft_precomputation_n);
        let constants_8n =
            polynomials_to_values_padded(&constant_polynomials, &fft_precomputation_8n);
        let c_constants = commit_polynomials(
            constant_polynomials.as_slice(),
            &pedersen_g_msm_precomputation,
            pedersen_h,
            false, // Circuit blinding is not necessary here.
        );

        // Convert sigma's values to scalar field elements and split it into degree-n chunks.
        let sigma_chunks = sigma_polynomials(sigma, degree, subgroup_generator_n);

        // Compute S_sigma, then a commitment to it.
        let s_sigma_polynomials = values_to_polynomials(&sigma_chunks, &fft_precomputation_n);
        let s_sigma_values_8n =
            polynomials_to_values_padded(&s_sigma_polynomials, &fft_precomputation_8n);
        let c_s_sigmas = commit_polynomials(
            s_sigma_polynomials.as_slice(),
            &pedersen_g_msm_precomputation,
            pedersen_h,
            false, // Circuit blinding is not necessary here.
        );

        Circuit {
            security_bits,
            num_public_inputs,
            gate_constants,
            routing_target_partitions,
            generators,
            subgroup_generator_n,
            subgroup_generator_8n,
            subgroup_n,
            subgroup_8n,
            pedersen_g,
            pedersen_h,
            u,
            constant_polynomials,
            constants_8n,
            c_constants,
            s_sigma_polynomials,
            s_sigma_values_8n,
            c_s_sigmas,
            pedersen_g_msm_precomputation,
            fft_precomputation_n,
            fft_precomputation_8n,
        }
    }

    fn get_routing_partitions(&self) -> TargetPartitions {
        let mut partitions = TargetPartitions::new();

        for i in 0..self.virtual_target_index {
            partitions.add_partition(Target::VirtualTarget(VirtualTarget { index: i }));
        }

        for gate in 0..self.num_gates() {
            for input in 0..NUM_WIRES {
                partitions.add_partition(Target::Wire(Wire { gate, input }));
            }
        }

        for &(a, b) in &self.copy_constraints {
            partitions.merge(a, b);
        }

        partitions
    }
}

#[cfg(test)]
mod tests {
    use anyhow::Result;

    use crate::{verify_proof, CircuitBuilder, Curve, CurveMulOp, Field, PartialWitness, Tweedledee, Tweedledum};

    #[test]
    // TODO: This fails because curve_mul_endo has a flaw.
    #[ignore]
    fn test_curve_mul_inv_endo() -> Result<()> {
        type C = Tweedledee;
        type InnerC = Tweedledum;
        type SF = <C as Curve>::ScalarField;
        let mut builder = CircuitBuilder::<Tweedledee>::new(128);

        // This is an arbitrary nonzero scalar.
        let scalar = builder.constant_wire(SF::FIVE);

        // Let p1 be an arbitrary point.
        let p1 = builder.constant_affine_point::<InnerC>(InnerC::GENERATOR_AFFINE);

        // Let p2 = [n(s)] p1.
        let mul_result = builder.curve_mul_endo::<InnerC>(CurveMulOp::new(scalar, p1));
        let p2 = mul_result.mul_result;

        // Let p3 = [1 / n(s)] p2.
        let mul_inv_result = builder.curve_mul_inv_endo::<InnerC>(CurveMulOp::new(scalar, p2));
        let p3 = mul_inv_result.mul_result;

        // Since the scalars cancel, p1 = p3.
        builder.copy_curve(p1, p3);

        let circuit = builder.build();
        let witness = circuit.generate_witness(PartialWitness::new());
        let proof = circuit.generate_proof::<InnerC>(witness, &[], false, true)?;
        let vk = circuit.to_vk();
        verify_proof::<C, InnerC>(&[], &proof, &[], &vk, true)?;

        Ok(())
    }
}
