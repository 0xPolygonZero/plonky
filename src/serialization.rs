use crate::{
    AffinePoint, Bls12377, Bls12377Base, Bls12377Scalar, Curve, Field, Tweedledee, TweedledeeBase,
    Tweedledum, TweedledumBase,
};
use std::io::{Error, ErrorKind, Read, Result, Write};

pub trait ToBytes {
    fn write<W: Write>(&self, writer: W) -> Result<()>;
}

pub trait FromBytes: Sized {
    fn read<R: Read>(reader: R) -> Result<Self>;
}

macro_rules! impl_field {
    ($field:ty) => {
        impl ToBytes for $field {
            fn write<W: Write>(&self, mut writer: W) -> Result<()> {
                writer.write_all(&self.to_canonical_u8_vec())
            }
        }
        impl FromBytes for $field {
            fn read<R: Read>(mut reader: R) -> Result<Self> {
                let mut buf = [0u8; Self::BYTES];
                reader.read_exact(&mut buf)?;
                Self::from_canonical_u8_vec(buf.to_vec())
                    .map_err(|e| Error::new(ErrorKind::Other, e.to_string()))
            }
        }
    };
}

impl_field!(TweedledeeBase);
impl_field!(TweedledumBase);
impl_field!(Bls12377Base);
impl_field!(Bls12377Scalar);

macro_rules! impl_curve {
    ($curve:ty, $field:ty) => {
        impl ToBytes for AffinePoint<$curve> {
            fn write<W: Write>(&self, mut writer: W) -> Result<()> {
                let zero = if self.zero { 1 } else { 0 };
                let odd = if self.y.to_canonical_u64_vec()[0] % 2 == 1 {
                    2
                } else {
                    0
                };
                let mask: u8 = zero | odd;
                writer.write(&[mask])?;
                writer.write_all(&self.x.to_canonical_u8_vec())
            }
        }

        impl FromBytes for AffinePoint<$curve> {
            fn read<R: Read>(mut reader: R) -> Result<Self> {
                let mut mask = vec![0u8];
                reader.read_exact(&mut mask)?;
                let mask = mask[0];
                if mask & 1 == 1 {
                    return Ok(AffinePoint {
                        x: <$field>::ZERO,
                        y: <$field>::ZERO,
                        zero: true,
                    });
                }
                let mut buf = [0u8; <$field>::BYTES];
                reader.read_exact(&mut buf)?;
                let x = <$field>::from_canonical_u8_vec(buf.to_vec())
                    .map_err(|e| Error::new(ErrorKind::Other, e.to_string()))?;
                let square_candidate = x.cube() + <$curve>::A * x + <$curve>::B;
                let y = square_candidate
                    .square_root()
                    .ok_or(Error::new(ErrorKind::Other, "Invalid x coordinate"))?;
                if (y.to_canonical_u64_vec()[0] % 2) as u8 == (mask & 2) >> 1 {
                    return Ok(AffinePoint::nonzero(x, y));
                } else {
                    return Ok(AffinePoint::nonzero(x, -y));
                }
            }
        }
    };
}

impl_curve!(Tweedledee, <Tweedledee as Curve>::BaseField);
impl_curve!(Tweedledum, <Tweedledum as Curve>::BaseField);
impl_curve!(Bls12377, <Bls12377 as Curve>::BaseField);

#[cfg(test)]
mod test {
    use super::*;
    use crate::blake_hash_base_field_to_curve;

    macro_rules! test_field_serialization {
        ($field:ty, $test_name:ident) => {
            #[test]
            fn $test_name() -> Result<()> {
                let x = <$field>::rand();
                let mut buf = [0u8; <$field>::BYTES];
                x.write(&mut buf[..])?;
                let y = <$field>::read(&buf[..])?;
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
}
