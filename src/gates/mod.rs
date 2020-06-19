//! For reference, here is our gate prefix tree:
//!
//! ```text
//! 101001 PublicInputGate
//! 101000 CurveAddGate
//! 10111* CurveDblGate
//! 11**** CurveEndoGate
//! 1000** Base4SumGate
//! 101010 BufferGate
//! 10110* ConstantGate
//! 1001** ArithmeticGate
//! 00**** RescueStepAGate
//! 01**** RescueStepBGate
//! ```
//!
//! The `*`s above represent constants which are not used in the gate prefix, and are thus available
//! for gate configuration.

pub use arithmetic::*;
pub use base_4_sum::*;
pub use buffer::*;
pub use constant::*;
pub use curve_add::*;
pub use curve_dbl::*;
pub use curve_endo::*;
pub use public_input::*;
pub use rescue_a::*;
pub use rescue_b::*;

use crate::{CircuitBuilder, Field, HaloCurve, Target, WitnessGenerator};

mod arithmetic;
mod base_4_sum;
mod buffer;
mod constant;
mod curve_add;
mod curve_dbl;
mod curve_endo;
mod public_input;
mod rescue_a;
mod rescue_b;

pub const RESCUE_SPONGE_WIDTH: usize = 4;

pub fn evaluate_all_constraints<C: HaloCurve, InnerC: HaloCurve<BaseField = C::ScalarField>>(
    local_constant_values: &[C::ScalarField],
    local_wire_values: &[C::ScalarField],
    right_wire_values: &[C::ScalarField],
    below_wire_values: &[C::ScalarField],
) -> Vec<C::ScalarField> {
    let constraint_sets_per_gate = vec![
        CurveAddGate::<C, InnerC>::evaluate_filtered(
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        CurveDblGate::<C, InnerC>::evaluate_filtered(
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        CurveEndoGate::<C, InnerC>::evaluate_filtered(
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        Base4SumGate::<C>::evaluate_filtered(
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        PublicInputGate::<C>::evaluate_filtered(
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        BufferGate::<C>::evaluate_filtered(
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        ConstantGate::<C>::evaluate_filtered(
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        ArithmeticGate::<C>::evaluate_filtered(
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        RescueStepAGate::<C>::evaluate_filtered(
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        RescueStepBGate::<C>::evaluate_filtered(
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
    ];

    let mut unified_constraint_set = vec![];
    for constraint_sets in constraint_sets_per_gate {
        while unified_constraint_set.len() < constraint_sets.len() {
            unified_constraint_set.push(C::ScalarField::ZERO);
        }
        for i in 0..constraint_sets.len() {
            unified_constraint_set[i] = unified_constraint_set[i] + constraint_sets[i];
        }
    }
    unified_constraint_set
}

pub fn evaluate_all_constraints_recursively<
    C: HaloCurve,
    InnerC: HaloCurve<BaseField = C::ScalarField>,
>(
    builder: &mut CircuitBuilder<C>,
    local_constant_values: &[Target],
    local_wire_values: &[Target],
    right_wire_values: &[Target],
    below_wire_values: &[Target],
) -> Vec<Target> {
    let constraint_sets_per_gate = vec![
        CurveAddGate::<C, InnerC>::evaluate_filtered_recursively(
            builder,
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        CurveDblGate::<C, InnerC>::evaluate_filtered_recursively(
            builder,
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        CurveEndoGate::<C, InnerC>::evaluate_filtered_recursively(
            builder,
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        Base4SumGate::<C>::evaluate_filtered_recursively(
            builder,
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        PublicInputGate::<C>::evaluate_filtered_recursively(
            builder,
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        BufferGate::<C>::evaluate_filtered_recursively(
            builder,
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        ConstantGate::<C>::evaluate_filtered_recursively(
            builder,
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        ArithmeticGate::<C>::evaluate_filtered_recursively(
            builder,
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        RescueStepAGate::<C>::evaluate_filtered_recursively(
            builder,
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
        RescueStepBGate::<C>::evaluate_filtered_recursively(
            builder,
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        ),
    ];

    let mut unified_constraint_set = vec![];
    for constraint_set in constraint_sets_per_gate {
        while unified_constraint_set.len() < constraint_set.len() {
            unified_constraint_set.push(builder.zero_wire());
        }
        for i in 0..constraint_set.len() {
            unified_constraint_set[i] = builder.add(unified_constraint_set[i], constraint_set[i]);
        }
    }
    unified_constraint_set
}

/// Computes `x * (x - 1)`, which should vanish iff `x` is binary.
fn assert_binary_recursively<C: HaloCurve>(builder: &mut CircuitBuilder<C>, x: Target) -> Target {
    let one = builder.one_wire();
    let x_minus_one = builder.sub(x, one);
    builder.mul(x, x_minus_one)
}

/// Computes `x * y - 1`, which should vanish iff `x` and `y` are inverses.
fn assert_inverses_recursively<C: HaloCurve>(
    builder: &mut CircuitBuilder<C>,
    x: Target,
    y: Target,
) -> Target {
    let one = builder.one_wire();
    let x_y = builder.mul(x, y);
    builder.sub(x_y, one)
}

pub trait Gate<C: HaloCurve>: WitnessGenerator<C::ScalarField> {
    const NAME: &'static str;

    /// In order to combine the constraints of various gate types into a unified constraint set, we
    /// assign each gate type a binary prefix such that no two prefixes overlap.
    const PREFIX: &'static [bool];

    fn evaluate_filtered(
        local_constant_values: &[C::ScalarField],
        local_wire_values: &[C::ScalarField],
        right_wire_values: &[C::ScalarField],
        below_wire_values: &[C::ScalarField],
    ) -> Vec<C::ScalarField> {
        let filter = Self::evaluate_prefix_filter(local_constant_values);
        let unfiltered = Self::evaluate_unfiltered(
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        );
        unfiltered.into_iter().map(|u| filter * u).collect()
    }

    fn evaluate_filtered_recursively(
        builder: &mut CircuitBuilder<C>,
        local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        below_wire_values: &[Target],
    ) -> Vec<Target> {
        let filter = Self::evaluate_prefix_filter_recursively(builder, local_constant_values);
        let unfiltered = Self::evaluate_unfiltered_recursively(
            builder,
            local_constant_values,
            local_wire_values,
            right_wire_values,
            below_wire_values,
        );
        unfiltered
            .into_iter()
            .map(|u| builder.mul(filter, u))
            .collect()
    }

    fn evaluate_prefix_filter(local_constant_values: &[C::ScalarField]) -> C::ScalarField {
        let mut product = C::ScalarField::ONE;
        for (i, &bit) in Self::PREFIX.iter().enumerate() {
            let c = local_constant_values[i];
            if bit {
                product = product * c;
            } else {
                product = product * (C::ScalarField::ONE - c);
            }
        }
        product
    }

    fn evaluate_prefix_filter_recursively(
        builder: &mut CircuitBuilder<C>,
        local_constant_values: &[Target],
    ) -> Target {
        let one = builder.one_wire();
        let mut product = one;
        for (i, &bit) in Self::PREFIX.iter().enumerate() {
            let c = local_constant_values[i];
            let term = if bit { c } else { builder.sub(one, c) };
            product = builder.mul(product, term);
        }
        product
    }

    /// Evaluate the constraints implied by this gate at the given challenge point.
    ///
    /// For example, if the gate computes `c = a * b`, this should return `[c(x) - a(x) * b(x)]`,
    /// where `x` is the challenge point.
    fn evaluate_unfiltered(
        local_constant_values: &[C::ScalarField],
        local_wire_values: &[C::ScalarField],
        right_wire_values: &[C::ScalarField],
        below_wire_values: &[C::ScalarField],
    ) -> Vec<C::ScalarField>;

    /// Like the other `evaluate` method, but in the context of a recursive circuit.
    fn evaluate_unfiltered_recursively(
        builder: &mut CircuitBuilder<C>,
        local_constant_values: &[Target],
        local_wire_values: &[Target],
        right_wire_values: &[Target],
        below_wire_values: &[Target],
    ) -> Vec<Target>;
}

/// Test that a gate's constraints are within degree 8n, including the gate prefix filter.
#[macro_export]
macro_rules! test_gate_low_degree {
    ($method:ident, $curve:ty, $gate:ty) => {
        #[test]
        #[ignore] // Too slow to run regularly.
        #[allow(non_snake_case)]
        fn $method() {
            type C = $curve;
            type SF = <C as $crate::curve::Curve>::ScalarField;

            let n = 256;
            let fft_precomputation_n = $crate::fft::fft_precompute::<SF>(n);
            let fft_precomputation_8n = $crate::fft::fft_precompute::<SF>(8 * n);
            let fft_precomputation_16n = $crate::fft::fft_precompute::<SF>(16 * n);

            // Generate random constant and wire polynomials.
            let mut constant_values_n: Vec<Vec<SF>> =
                vec![Vec::new(); $crate::plonk::NUM_CONSTANTS];
            let mut wire_values_n: Vec<Vec<SF>> = vec![Vec::new(); $crate::plonk::NUM_WIRES];
            for i in 0..n {
                for points in constant_values_n.iter_mut() {
                    points.push(<SF as $crate::field::Field>::rand())
                }
                for points in wire_values_n.iter_mut() {
                    points.push(<SF as $crate::field::Field>::rand())
                }
            }

            // Low-degree extend them to 16n values.
            let mut constant_coeffs_16n =
                $crate::plonk_util::values_to_coeffs(&constant_values_n, &fft_precomputation_n);
            let mut wire_coeffs_16n =
                $crate::plonk_util::values_to_coeffs(&wire_values_n, &fft_precomputation_n);
            for coeffs in constant_coeffs_16n
                .iter_mut()
                .chain(wire_coeffs_16n.iter_mut())
            {
                while coeffs.len() < 16 * n {
                    coeffs.push(<SF as $crate::field::Field>::ZERO);
                }
            }
            let constant_values_16n: Vec<Vec<SF>> = constant_coeffs_16n
                .iter()
                .map(|coeffs| {
                    $crate::fft::fft_with_precomputation_power_of_2(coeffs, &fft_precomputation_16n)
                })
                .collect();
            let wire_values_16n: Vec<Vec<SF>> = wire_coeffs_16n
                .iter()
                .map(|coeffs| {
                    $crate::fft::fft_with_precomputation_power_of_2(coeffs, &fft_precomputation_16n)
                })
                .collect();

            // Make sure each extended polynomial is still degree n.
            for values_16n in constant_values_16n.iter().chain(wire_values_16n.iter()) {
                assert!(
                    $crate::plonk_util::polynomial_degree(values_16n, &fft_precomputation_16n) < n
                );
            }

            let constant_values_16n_t = $crate::util::transpose(&constant_values_16n);
            let wire_values_16n_t = $crate::util::transpose(&wire_values_16n);

            // Evaluate constraints at each of our 16n points.
            let mut constraint_values_16n: Vec<Vec<SF>> = Vec::new();
            for i in 0..16 * n {
                let constraints: Vec<SF> =
                    <$gate as $crate::gates::Gate<C>>::evaluate_filtered(
                        &constant_values_16n_t[i],
                        &wire_values_16n_t[i],
                        &wire_values_16n_t[(i + 16) % (16 * n)],
                        &wire_values_16n_t[(i + 16 * $crate::plonk::GRID_WIDTH) % (16 * n)],
                    );
                for (j, &c) in constraints.iter().enumerate() {
                    if constraint_values_16n.len() <= j {
                        constraint_values_16n.push(Vec::new());
                    }
                    constraint_values_16n[j].push(c);
                }
            }

            // Check that the degree of each constraint is within the limit.
            let constraint_degrees = constraint_values_16n
                .iter()
                .map(|c| $crate::plonk_util::polynomial_degree(c, &fft_precomputation_16n))
                .collect::<Vec<_>>();
            let max_degree_excl = (crate::plonk::QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER + 1) * n;
            for (i, &deg) in constraint_degrees.iter().enumerate() {
                assert!(
                    deg < max_degree_excl,
                    "Constraint at index {} has degree {}; should be less than {}n = {}",
                    i,
                    deg,
                    crate::plonk::QUOTIENT_POLYNOMIAL_DEGREE_MULTIPLIER + 1,
                    max_degree_excl
                );
            }
        }
    };
}
