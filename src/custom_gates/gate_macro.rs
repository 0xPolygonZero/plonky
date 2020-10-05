use crate::custom_gates::multivariate_polynomial::MultivariatePolynomial;
use crate::{CircuitBuilder, Curve, Field, Gate, HaloCurve, PartialWitness, Target, Wire, WitnessGenerator, NUM_WIRES};
use std::marker::PhantomData;

fn evaluate_polys<C: Curve>(
    polys: &[MultivariatePolynomial<C::ScalarField, NUM_WIRES>],
    local_constant_values: &[C::ScalarField],
    local_wire_values: &[C::ScalarField],
    right_wire_values: &[C::ScalarField],
    below_wire_values: &[C::ScalarField],
) -> Vec<C::ScalarField> {
    polys.iter().map(|p| {
        let mut xs = [C::ScalarField::ZERO; NUM_WIRES];
        xs.copy_from_slice(local_constant_values);
        p.eval(xs)
    }).collect()
}

macro_rules! gate {
    ($name: ident, $name_string: expr, $polynomials: expr) => {
        pub struct $name<C: HaloCurve, InnerC: Curve<BaseField = C::ScalarField>> {
            pub index: usize,
            _phantom_oc: PhantomData<C>,
            _phantom_ic: PhantomData<InnerC>,
        }
        impl<C: HaloCurve, InnerC: Curve<BaseField = C::ScalarField>> Gate<C> for $name<C, InnerC> {
            const NAME: &'static str = $name_string;
            const DEGREE: usize = 3;
            const NUM_CONSTANTS: usize = 2;

            const PREFIX: &'static [bool] = &[true, false, false, true];

            fn evaluate_unfiltered(
                local_constant_values: &[C::ScalarField],
                local_wire_values: &[C::ScalarField],
                right_wire_values: &[C::ScalarField],
                below_wire_values: &[C::ScalarField],
            ) -> Vec<C::ScalarField> {
                evaluate_polys(&$polynomials,
                local_constant_values,
                local_wire_values,
                right_wire_values,
                below_wire_values
                )
            }

            fn evaluate_unfiltered_recursively(
                builder: &mut CircuitBuilder<C>,
                _local_constant_values: &[Target<C::ScalarField>],
                local_wire_values: &[Target<C::ScalarField>],
                right_wire_values: &[Target<C::ScalarField>],
                _below_wire_values: &[Target<C::ScalarField>],
            ) -> Vec<Target<C::ScalarField>> {
                todo!()
            }
        }

        impl<C: HaloCurve, InnerC: Curve<BaseField = C::ScalarField>>
            WitnessGenerator<C::ScalarField> for $name<C, InnerC>
        {
            fn dependencies(&self) -> Vec<Target<C::ScalarField>> {
                todo!()
            }

            fn generate(
                &self,
                _constants: &[Vec<C::ScalarField>],
                witness: &PartialWitness<InnerC::BaseField>,
            ) -> PartialWitness<InnerC::BaseField> {
                todo!()
            }
        }
    };
}

gate!(BamGate, "BamGate", poly);

fn bam<C: Curve>() {
    // let polys = Vec::<MultivariatePolynomial<C::ScalarField, NUM_WIRES>>::new();
    let poly = MultivariatePolynomial::<C::ScalarField, NUM_WIRES>::zero();
    dbg!(&polys);
}

// fn polys_to_evaluation<C: Curve>(
//     polys: &[MultivariatePolynomial<C::ScalarField, NUM_WIRES>],
// ) -> Box<
//     dyn Fn(
//         &[C::ScalarField],
//         &[C::ScalarField],
//         &[C::ScalarField],
//         &[C::ScalarField],
//     ) -> Vec<C::ScalarField>,
// > {
//     let eval = |(
//          &[C::ScalarField],
//          &[C::ScalarField],
//          &[C::ScalarField],
//          &[C::ScalarField],
//     )| {
//         polys.iter().map(|p| {
//
//         }).collect::<Vec<_>>()
//     };
//     Box::new(eval)
// }

