use crate::{AffinePoint, Curve, Field};
use serde::de::Error as DeError;
use serde::de::Visitor;
use serde::ser::Error as SerdeError;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::fmt;
use std::io::{Error, ErrorKind, Read, Result, Write};

pub trait ToBytes {
    fn write<W: Write>(&self, writer: W) -> Result<()>;
}

pub trait FromBytes: Sized {
    fn read<R: Read>(reader: R) -> Result<Self>;
}

impl<F: Field> ToBytes for F {
    fn write<W: Write>(&self, mut writer: W) -> Result<()> {
        writer.write_all(&self.to_canonical_u8_vec())
    }
}

impl<F: Field> FromBytes for F {
    fn read<R: Read>(mut reader: R) -> Result<Self> {
        let mut buf = vec![0u8; Self::BYTES];
        reader.read_exact(&mut buf)?;
        Self::from_canonical_u8_vec(buf.to_vec())
            .map_err(|e| Error::new(ErrorKind::Other, e.to_string()))
    }
}

impl<C: Curve> ToBytes for AffinePoint<C> {
    fn write<W: Write>(&self, mut writer: W) -> Result<()> {
        let zero = if self.zero { 1 } else { 0 };
        let odd = if self.y.to_canonical_u64_vec()[0] % 2 == 1 {
            2
        } else {
            0
        };
        let mask: u8 = zero | odd;
        writer.write_all(&[mask])?;
        writer.write_all(&self.x.to_canonical_u8_vec())
    }
}

impl<C: Curve> FromBytes for AffinePoint<C> {
    fn read<R: Read>(mut reader: R) -> Result<Self> {
        let mut mask = vec![0u8];
        reader.read_exact(&mut mask)?;
        let mask = mask[0];
        if mask & 1 == 1 {
            return Ok(AffinePoint {
                x: C::BaseField::ZERO,
                y: C::BaseField::ZERO,
                zero: true,
            });
        }
        let mut buf = vec![0u8; C::BaseField::BYTES];
        reader.read_exact(&mut buf)?;
        let x = C::BaseField::from_canonical_u8_vec(buf.to_vec())
            .map_err(|e| Error::new(ErrorKind::Other, e.to_string()))?;
        let square_candidate = x.cube() + C::A * x + C::B;
        let y = square_candidate
            .square_root()
            .ok_or_else(|| Error::new(ErrorKind::Other, "Invalid x coordinate"))?;
        if (y.to_canonical_u64_vec()[0] % 2) as u8 == (mask & 2) >> 1 {
            Ok(AffinePoint::nonzero(x, y))
        } else {
            Ok(AffinePoint::nonzero(x, -y))
        }
    }
}

impl<C: Curve> Serialize for AffinePoint<C> {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut buf = vec![];
        self.write(&mut buf)
            .map_err(|e| S::Error::custom(format!("{}", e)))?;
        serializer.serialize_bytes(&buf)
    }
}

impl<'de, C: Curve> Deserialize<'de> for AffinePoint<C> {
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct AffinePointVisitor<C: Curve> {
            phantom: std::marker::PhantomData<C>,
        }

        impl<'de, C: Curve> Visitor<'de> for AffinePointVisitor<C> {
            type Value = AffinePoint<C>;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                write!(formatter, "An affine point.")
            }

            fn visit_bytes<E: DeError>(self, v: &[u8]) -> std::result::Result<Self::Value, E> {
                AffinePoint::<C>::read(v).map_err(|e| DeError::custom(format!("{}", e)))
            }
        }
        deserializer.deserialize_bytes(AffinePointVisitor {
            phantom: std::marker::PhantomData,
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{blake_hash_base_field_to_curve, Bls12377, Bls12377Base, Bls12377Scalar, CircuitBuilder, HaloCurve, PartialWitness, Proof, Tweedledee, TweedledeeBase, Tweedledum, TweedledumBase, VerificationKey};

    macro_rules! test_field_serialization {
        ($field:ty, $test_name:ident) => {
            #[test]
            fn $test_name() -> Result<()> {
                let x = <$field>::rand();
                let mut buf = [0u8; <$field>::BYTES];
                x.write(&mut buf[..])?;
                let y = <$field>::read(&buf[..])?;
                assert_eq!(x, y);

                // Serde (de)serialization
                let ser = bincode::serialize(&x).unwrap();
                let y = bincode::deserialize(&ser).unwrap();
                assert_eq!(x, y);

                Ok(())
            }
        };
    }

    macro_rules! test_curve_serialization {
        ($curve:ty, $basefield:ty, $test_name:ident) => {
            #[test]
            fn $test_name() -> Result<()> {
                // Random element
                let f = <$basefield>::rand();
                let p = blake_hash_base_field_to_curve::<$curve>(f);
                let mut buf = [0u8; <$basefield>::BYTES + 1];
                p.write(&mut buf[..])?;
                let q = AffinePoint::<$curve>::read(&buf[..])?;
                assert_eq!(p, q);
                // Point at infinity
                let zero = AffinePoint::<$curve> {
                    x: <$basefield>::ZERO,
                    y: <$basefield>::ZERO,
                    zero: true,
                };
                let mut buf = [0u8; <$basefield>::BYTES + 1];
                zero.write(&mut buf[..])?;
                let q = AffinePoint::<$curve>::read(&buf[..])?;
                assert_eq!(zero, q);

                // Serde (de)serialization
                let ser = bincode::serialize(&p).unwrap();
                let q = bincode::deserialize(&ser).unwrap();
                assert_eq!(p, q);
                let ser = bincode::serialize(&zero).unwrap();
                let q = bincode::deserialize(&ser).unwrap();
                assert_eq!(zero, q);

                Ok(())
            }
        };
    }

    test_field_serialization!(TweedledeeBase, test_tweedledee_base_serialization);
    test_field_serialization!(TweedledumBase, test_tweedledum_base_serialization);
    test_field_serialization!(Bls12377Base, test_bls_base_serialization);
    test_field_serialization!(Bls12377Scalar, test_bls_scalar_serialization);
    test_curve_serialization!(
        Tweedledee,
        <Tweedledee as Curve>::BaseField,
        test_tweedledee_curve_serialization
    );
    test_curve_serialization!(
        Tweedledum,
        <Tweedledum as Curve>::BaseField,
        test_tweedledum_curve_serialization
    );
    test_curve_serialization!(
        Bls12377,
        <Bls12377 as Curve>::BaseField,
        test_bls_curve_serialization
    );

    // Generate a proof and verification key for the factorial circuit.
    fn get_circuit_vk<C: HaloCurve, InnerC: HaloCurve<BaseField = C::ScalarField>>(
    ) -> (Proof<C>, VerificationKey<C>) {
        let mut builder = CircuitBuilder::<C>::new(128);
        let n = 10;
        let factorial_usize = (1..=n).product();
        let factors_pis = builder.stage_public_inputs(n);
        builder.route_public_inputs();
        let res = builder.constant_wire(C::ScalarField::from_canonical_usize(factorial_usize));
        let mut factorial = builder.one_wire();
        factors_pis.iter().for_each(|pi| {
            factorial = builder.mul(factorial, pi.routable_target());
        });
        builder.copy(factorial, res);
        let mut partial_witness = PartialWitness::new();
        (0..n).for_each(|i| {
            partial_witness
                .set_public_input(factors_pis[i], C::ScalarField::from_canonical_usize(i + 1));
        });
        let circuit = builder.build();
        let witness = circuit.generate_witness(partial_witness);
        let proof = circuit
            .generate_proof::<InnerC>(&witness, &[], true)
            .unwrap();
        let vk = circuit.to_vk();
        (proof, vk)
    }

    #[test]
    fn test_proof_vk_serialization_tweedledee() {
        let (proof, vk) = get_circuit_vk::<Tweedledee, Tweedledum>();
        let ser_proof = bincode::serialize(&proof).unwrap();
        let ser_vk = bincode::serialize(&vk).unwrap();

        let der_proof = bincode::deserialize(&ser_proof).unwrap();
        let der_vk: VerificationKey<Tweedledee> = bincode::deserialize(&ser_vk).unwrap();

        assert_eq!(proof, der_proof);
        assert_eq!(vk.c_constants, der_vk.c_constants);
        assert_eq!(vk.c_s_sigmas, der_vk.c_s_sigmas);
        assert_eq!(vk.degree, der_vk.degree);
        assert_eq!(vk.num_public_inputs, der_vk.num_public_inputs);
        assert_eq!(vk.security_bits, der_vk.security_bits);
    }

    #[test]
    fn test_proof_vk_serialization_tweedledum() {
        let (proof, vk) = get_circuit_vk::<Tweedledum, Tweedledee>();
        let ser_proof = bincode::serialize(&proof).unwrap();
        let ser_vk = bincode::serialize(&vk).unwrap();

        let der_proof = bincode::deserialize(&ser_proof).unwrap();
        let der_vk: VerificationKey<Tweedledum> = bincode::deserialize(&ser_vk).unwrap();

        assert_eq!(proof, der_proof);
        assert_eq!(vk.c_constants, der_vk.c_constants);
        assert_eq!(vk.c_s_sigmas, der_vk.c_s_sigmas);
        assert_eq!(vk.degree, der_vk.degree);
        assert_eq!(vk.num_public_inputs, der_vk.num_public_inputs);
        assert_eq!(vk.security_bits, der_vk.security_bits);
    }
}
