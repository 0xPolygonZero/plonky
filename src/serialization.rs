use crate::{Field, TweedledeeBase, TweedledumBase, Bls12377Base, Bls12377Scalar};
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

macro_rules! test_field_serialization {
    ($field:ty, $test_name:ident) => {
        #[test]
        fn $test_name() -> Result<()> {
            let x = <$field>::rand();
            let mut buf: Vec<u8> = Vec::new();
            x.write(&mut buf)?;
            let mut arr = [0u8; <$field>::BYTES];
            arr.copy_from_slice(buf.as_slice());
            let y = <$field>::read(&arr[..])?;
            assert_eq!(x, y);
            Ok(())
        }
    };
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{TweedledeeBase, TweedledumBase, Bls12377Base, Bls12377Scalar};
    test_field_serialization!(TweedledeeBase, test_tweedledee_base_serialization);
    test_field_serialization!(TweedledumBase, test_tweedledum_base_serialization);
    test_field_serialization!(Bls12377Base, test_bls_base_serialization);
    test_field_serialization!(Bls12377Scalar, test_bls_scalar_serialization);
}

// #[test]
// fn lol() -> Result<()> {
//     let x = TweedledeeBase::rand();
//     let mut buf: Vec<u8> = Vec::new();
//     x.write(&mut buf)?;
//     let mut arr = [0u8; 32];
//     arr.copy_from_slice(buf.as_slice());
//     let y = TweedledeeBase::read(&arr[..])?;
//     assert_eq!(x, y);
//     Ok(())
// }
