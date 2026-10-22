//! A compact binary encoding of serde types, for the parts of a `.fieldcraft` file that hold
//! large numeric data with a structure of its own ([`crate::results`]).
//!
//! Little-endian, not self-describing: the type being read says what comes next.
//!
//! - `bool` one byte (0 or 1); integers and floats their fixed width (`usize` as `u64`); `char`
//!   as `u32`.
//! - A string or byte string: `u64` length, then the bytes (UTF-8 for a string).
//! - `Option`: one byte, 0 for none, 1 then the value.
//! - A sequence or map: `u64` count, then the elements (a map's as key, value).
//! - A tuple, tuple struct or array: its elements, no count.
//! - A struct: per field in declaration order a presence byte (1, then the value; 0 for a field
//!   its type leaves out, `skip_serializing_if`), so a field left out reads back as its default.
//! - An enum: `u32` variant index, then the variant's content as above.
//!
//! The encoding follows the types' field and variant order: a part encoded this way carries a
//! version of its own, which changes with any change to those types.
//!
//! Reading treats the bytes as hostile: every count is checked against the bytes left before
//! anything is read for it (each element takes at least one byte), so nothing is allocated
//! beyond a small multiple of the input; a malformed value, missing bytes or bytes left over
//! are errors.

use serde::de::{self, DeserializeSeed, IntoDeserializer, Visitor};
use serde::ser::{self, Serialize};

/// An encoding or decoding error.
#[derive(Debug, Clone, PartialEq)]
pub struct CodecError(pub String);

impl std::fmt::Display for CodecError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for CodecError {}

impl ser::Error for CodecError {
    fn custom<T: std::fmt::Display>(msg: T) -> Self {
        CodecError(msg.to_string())
    }
}

impl de::Error for CodecError {
    fn custom<T: std::fmt::Display>(msg: T) -> Self {
        CodecError(msg.to_string())
    }
}

type Result<T> = std::result::Result<T, CodecError>;

/// `value` encoded, appended to `out`.
pub fn encode_into<T: Serialize + ?Sized>(value: &T, out: &mut Vec<u8>) -> Result<()> {
    value.serialize(&mut Encoder { out })
}

/// `value` encoded.
pub fn encode<T: Serialize + ?Sized>(value: &T) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    encode_into(value, &mut out)?;
    Ok(out)
}

/// A `T` decoded from all of `bytes` (bytes left over are an error).
pub fn decode<'de, T: de::Deserialize<'de>>(bytes: &'de [u8]) -> Result<T> {
    let mut d = Decoder { input: bytes };
    let value = T::deserialize(&mut d)?;
    if !d.input.is_empty() {
        return Err(CodecError(format!("{} bytes left over after the value", d.input.len())));
    }
    Ok(value)
}

struct Encoder<'a> {
    out: &'a mut Vec<u8>,
}

impl Encoder<'_> {
    fn len(&mut self, n: usize) {
        self.out.extend_from_slice(&(n as u64).to_le_bytes());
    }
}

/// A sequence or map must say its length up front.
fn unknown_length() -> CodecError {
    CodecError("a sequence or map of unknown length cannot be encoded".into())
}

impl ser::Serializer for &mut Encoder<'_> {
    type Ok = ();
    type Error = CodecError;
    type SerializeSeq = Self;
    type SerializeTuple = Self;
    type SerializeTupleStruct = Self;
    type SerializeTupleVariant = Self;
    type SerializeMap = Self;
    type SerializeStruct = Self;
    type SerializeStructVariant = Self;

    fn serialize_bool(self, v: bool) -> Result<()> {
        self.out.push(u8::from(v));
        Ok(())
    }
    fn serialize_i8(self, v: i8) -> Result<()> {
        self.out.extend_from_slice(&v.to_le_bytes());
        Ok(())
    }
    fn serialize_i16(self, v: i16) -> Result<()> {
        self.out.extend_from_slice(&v.to_le_bytes());
        Ok(())
    }
    fn serialize_i32(self, v: i32) -> Result<()> {
        self.out.extend_from_slice(&v.to_le_bytes());
        Ok(())
    }
    fn serialize_i64(self, v: i64) -> Result<()> {
        self.out.extend_from_slice(&v.to_le_bytes());
        Ok(())
    }
    fn serialize_u8(self, v: u8) -> Result<()> {
        self.out.push(v);
        Ok(())
    }
    fn serialize_u16(self, v: u16) -> Result<()> {
        self.out.extend_from_slice(&v.to_le_bytes());
        Ok(())
    }
    fn serialize_u32(self, v: u32) -> Result<()> {
        self.out.extend_from_slice(&v.to_le_bytes());
        Ok(())
    }
    fn serialize_u64(self, v: u64) -> Result<()> {
        self.out.extend_from_slice(&v.to_le_bytes());
        Ok(())
    }
    fn serialize_f32(self, v: f32) -> Result<()> {
        self.out.extend_from_slice(&v.to_le_bytes());
        Ok(())
    }
    fn serialize_f64(self, v: f64) -> Result<()> {
        self.out.extend_from_slice(&v.to_le_bytes());
        Ok(())
    }
    fn serialize_char(self, v: char) -> Result<()> {
        self.serialize_u32(u32::from(v))
    }
    fn serialize_str(self, v: &str) -> Result<()> {
        self.serialize_bytes(v.as_bytes())
    }
    fn serialize_bytes(self, v: &[u8]) -> Result<()> {
        self.len(v.len());
        self.out.extend_from_slice(v);
        Ok(())
    }
    fn serialize_none(self) -> Result<()> {
        self.out.push(0);
        Ok(())
    }
    fn serialize_some<T: Serialize + ?Sized>(self, value: &T) -> Result<()> {
        self.out.push(1);
        value.serialize(self)
    }
    fn serialize_unit(self) -> Result<()> {
        Ok(())
    }
    fn serialize_unit_struct(self, _: &'static str) -> Result<()> {
        Ok(())
    }
    fn serialize_unit_variant(self, _: &'static str, index: u32, _: &'static str) -> Result<()> {
        self.serialize_u32(index)
    }
    fn serialize_newtype_struct<T: Serialize + ?Sized>(self, _: &'static str, value: &T) -> Result<()> {
        value.serialize(self)
    }
    fn serialize_newtype_variant<T: Serialize + ?Sized>(self, _: &'static str, index: u32, _: &'static str, value: &T) -> Result<()> {
        self.serialize_u32(index)?;
        value.serialize(self)
    }
    fn serialize_seq(self, len: Option<usize>) -> Result<Self> {
        self.len(len.ok_or_else(unknown_length)?);
        Ok(self)
    }
    fn serialize_tuple(self, _: usize) -> Result<Self> {
        Ok(self)
    }
    fn serialize_tuple_struct(self, _: &'static str, _: usize) -> Result<Self> {
        Ok(self)
    }
    fn serialize_tuple_variant(self, _: &'static str, index: u32, _: &'static str, _: usize) -> Result<Self> {
        self.serialize_u32(index)?;
        Ok(self)
    }
    fn serialize_map(self, len: Option<usize>) -> Result<Self> {
        self.len(len.ok_or_else(unknown_length)?);
        Ok(self)
    }
    fn serialize_struct(self, _: &'static str, _: usize) -> Result<Self> {
        Ok(self)
    }
    fn serialize_struct_variant(self, _: &'static str, index: u32, _: &'static str, _: usize) -> Result<Self> {
        self.serialize_u32(index)?;
        Ok(self)
    }
    fn is_human_readable(&self) -> bool {
        false
    }
}

impl ser::SerializeSeq for &mut Encoder<'_> {
    type Ok = ();
    type Error = CodecError;
    fn serialize_element<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<()> {
        value.serialize(&mut **self)
    }
    fn end(self) -> Result<()> {
        Ok(())
    }
}

impl ser::SerializeTuple for &mut Encoder<'_> {
    type Ok = ();
    type Error = CodecError;
    fn serialize_element<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<()> {
        value.serialize(&mut **self)
    }
    fn end(self) -> Result<()> {
        Ok(())
    }
}

impl ser::SerializeTupleStruct for &mut Encoder<'_> {
    type Ok = ();
    type Error = CodecError;
    fn serialize_field<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<()> {
        value.serialize(&mut **self)
    }
    fn end(self) -> Result<()> {
        Ok(())
    }
}

impl ser::SerializeTupleVariant for &mut Encoder<'_> {
    type Ok = ();
    type Error = CodecError;
    fn serialize_field<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<()> {
        value.serialize(&mut **self)
    }
    fn end(self) -> Result<()> {
        Ok(())
    }
}

impl ser::SerializeMap for &mut Encoder<'_> {
    type Ok = ();
    type Error = CodecError;
    fn serialize_key<T: Serialize + ?Sized>(&mut self, key: &T) -> Result<()> {
        key.serialize(&mut **self)
    }
    fn serialize_value<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<()> {
        value.serialize(&mut **self)
    }
    fn end(self) -> Result<()> {
        Ok(())
    }
}

impl ser::SerializeStruct for &mut Encoder<'_> {
    type Ok = ();
    type Error = CodecError;
    fn serialize_field<T: Serialize + ?Sized>(&mut self, _: &'static str, value: &T) -> Result<()> {
        self.out.push(1);
        value.serialize(&mut **self)
    }
    fn skip_field(&mut self, _: &'static str) -> Result<()> {
        self.out.push(0);
        Ok(())
    }
    fn end(self) -> Result<()> {
        Ok(())
    }
}

impl ser::SerializeStructVariant for &mut Encoder<'_> {
    type Ok = ();
    type Error = CodecError;
    fn serialize_field<T: Serialize + ?Sized>(&mut self, _: &'static str, value: &T) -> Result<()> {
        self.out.push(1);
        value.serialize(&mut **self)
    }
    fn skip_field(&mut self, _: &'static str) -> Result<()> {
        self.out.push(0);
        Ok(())
    }
    fn end(self) -> Result<()> {
        Ok(())
    }
}

struct Decoder<'de> {
    input: &'de [u8],
}

impl<'de> Decoder<'de> {
    fn take(&mut self, n: usize) -> Result<&'de [u8]> {
        if n > self.input.len() {
            return Err(CodecError(format!("the data ends early: {n} bytes wanted, {} left", self.input.len())));
        }
        let (head, rest) = self.input.split_at(n);
        self.input = rest;
        Ok(head)
    }

    fn array<const N: usize>(&mut self) -> Result<[u8; N]> {
        let mut a = [0u8; N];
        a.copy_from_slice(self.take(N)?);
        Ok(a)
    }

    fn u8(&mut self) -> Result<u8> {
        Ok(self.array::<1>()?[0])
    }

    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.array()?))
    }

    fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.array()?))
    }

    /// A count of elements that each take at least one byte: no more than the bytes left.
    fn count(&mut self) -> Result<usize> {
        let n = self.u64()?;
        match usize::try_from(n) {
            Ok(n) if n <= self.input.len() => Ok(n),
            _ => Err(CodecError(format!("a count of {n} exceeds the {} bytes left", self.input.len()))),
        }
    }

    fn bytes(&mut self) -> Result<&'de [u8]> {
        let n = self.count()?;
        self.take(n)
    }

    fn str(&mut self) -> Result<&'de str> {
        std::str::from_utf8(self.bytes()?).map_err(|e| CodecError(format!("a string is not UTF-8: {e}")))
    }
}

fn not_self_describing() -> CodecError {
    CodecError("this encoding is read by type only".into())
}

impl<'de> de::Deserializer<'de> for &mut Decoder<'de> {
    type Error = CodecError;

    fn deserialize_any<V: Visitor<'de>>(self, _: V) -> Result<V::Value> {
        Err(not_self_describing())
    }
    fn deserialize_bool<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        match self.u8()? {
            0 => v.visit_bool(false),
            1 => v.visit_bool(true),
            b => Err(CodecError(format!("{b} is not a boolean"))),
        }
    }
    fn deserialize_i8<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        v.visit_i8(i8::from_le_bytes(self.array()?))
    }
    fn deserialize_i16<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        v.visit_i16(i16::from_le_bytes(self.array()?))
    }
    fn deserialize_i32<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        v.visit_i32(i32::from_le_bytes(self.array()?))
    }
    fn deserialize_i64<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        v.visit_i64(i64::from_le_bytes(self.array()?))
    }
    fn deserialize_u8<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        v.visit_u8(self.u8()?)
    }
    fn deserialize_u16<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        v.visit_u16(u16::from_le_bytes(self.array()?))
    }
    fn deserialize_u32<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        v.visit_u32(self.u32()?)
    }
    fn deserialize_u64<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        v.visit_u64(self.u64()?)
    }
    fn deserialize_f32<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        v.visit_f32(f32::from_le_bytes(self.array()?))
    }
    fn deserialize_f64<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        v.visit_f64(f64::from_le_bytes(self.array()?))
    }
    fn deserialize_char<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        let n = self.u32()?;
        v.visit_char(char::from_u32(n).ok_or_else(|| CodecError(format!("{n} is not a character")))?)
    }
    fn deserialize_str<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        v.visit_borrowed_str(self.str()?)
    }
    fn deserialize_string<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        v.visit_borrowed_str(self.str()?)
    }
    fn deserialize_bytes<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        v.visit_borrowed_bytes(self.bytes()?)
    }
    fn deserialize_byte_buf<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        v.visit_borrowed_bytes(self.bytes()?)
    }
    fn deserialize_option<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        match self.u8()? {
            0 => v.visit_none(),
            1 => v.visit_some(self),
            b => Err(CodecError(format!("{b} is not an option tag"))),
        }
    }
    fn deserialize_unit<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        v.visit_unit()
    }
    fn deserialize_unit_struct<V: Visitor<'de>>(self, _: &'static str, v: V) -> Result<V::Value> {
        v.visit_unit()
    }
    fn deserialize_newtype_struct<V: Visitor<'de>>(self, _: &'static str, v: V) -> Result<V::Value> {
        v.visit_newtype_struct(self)
    }
    fn deserialize_seq<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        let left = self.count()?;
        v.visit_seq(Elements { d: self, left })
    }
    fn deserialize_tuple<V: Visitor<'de>>(self, len: usize, v: V) -> Result<V::Value> {
        v.visit_seq(Elements { d: self, left: len })
    }
    fn deserialize_tuple_struct<V: Visitor<'de>>(self, _: &'static str, len: usize, v: V) -> Result<V::Value> {
        v.visit_seq(Elements { d: self, left: len })
    }
    fn deserialize_map<V: Visitor<'de>>(self, v: V) -> Result<V::Value> {
        let left = self.count()?;
        v.visit_map(Entries { d: self, left })
    }
    fn deserialize_struct<V: Visitor<'de>>(self, _: &'static str, fields: &'static [&'static str], v: V) -> Result<V::Value> {
        v.visit_map(Fields { d: self, fields })
    }
    fn deserialize_enum<V: Visitor<'de>>(self, _: &'static str, _: &'static [&'static str], v: V) -> Result<V::Value> {
        v.visit_enum(Variant { d: self })
    }
    fn deserialize_identifier<V: Visitor<'de>>(self, _: V) -> Result<V::Value> {
        Err(not_self_describing())
    }
    fn deserialize_ignored_any<V: Visitor<'de>>(self, _: V) -> Result<V::Value> {
        Err(not_self_describing())
    }
    fn is_human_readable(&self) -> bool {
        false
    }
}

/// The elements of a sequence or tuple.
struct Elements<'a, 'de> {
    d: &'a mut Decoder<'de>,
    left: usize,
}

impl<'de> de::SeqAccess<'de> for Elements<'_, 'de> {
    type Error = CodecError;
    fn next_element_seed<T: DeserializeSeed<'de>>(&mut self, seed: T) -> Result<Option<T::Value>> {
        if self.left == 0 {
            return Ok(None);
        }
        self.left -= 1;
        seed.deserialize(&mut *self.d).map(Some)
    }
    fn size_hint(&self) -> Option<usize> {
        Some(self.left)
    }
}

/// The entries of a map.
struct Entries<'a, 'de> {
    d: &'a mut Decoder<'de>,
    left: usize,
}

impl<'de> de::MapAccess<'de> for Entries<'_, 'de> {
    type Error = CodecError;
    fn next_key_seed<K: DeserializeSeed<'de>>(&mut self, seed: K) -> Result<Option<K::Value>> {
        if self.left == 0 {
            return Ok(None);
        }
        self.left -= 1;
        seed.deserialize(&mut *self.d).map(Some)
    }
    fn next_value_seed<V: DeserializeSeed<'de>>(&mut self, seed: V) -> Result<V::Value> {
        seed.deserialize(&mut *self.d)
    }
    fn size_hint(&self) -> Option<usize> {
        Some(self.left)
    }
}

/// A struct's fields, in declaration order, each behind its presence byte.
struct Fields<'a, 'de> {
    d: &'a mut Decoder<'de>,
    fields: &'static [&'static str],
}

impl<'de> de::MapAccess<'de> for Fields<'_, 'de> {
    type Error = CodecError;
    fn next_key_seed<K: DeserializeSeed<'de>>(&mut self, seed: K) -> Result<Option<K::Value>> {
        while let Some((name, rest)) = self.fields.split_first() {
            self.fields = rest;
            match self.d.u8()? {
                0 => continue,
                1 => return seed.deserialize(name.into_deserializer()).map(Some),
                b => return Err(CodecError(format!("{b} is not a field tag (field `{name}`)"))),
            }
        }
        Ok(None)
    }
    fn next_value_seed<V: DeserializeSeed<'de>>(&mut self, seed: V) -> Result<V::Value> {
        seed.deserialize(&mut *self.d)
    }
}

/// An enum: its variant index, then its content.
struct Variant<'a, 'de> {
    d: &'a mut Decoder<'de>,
}

impl<'de> de::EnumAccess<'de> for Variant<'_, 'de> {
    type Error = CodecError;
    type Variant = Self;
    fn variant_seed<V: DeserializeSeed<'de>>(self, seed: V) -> Result<(V::Value, Self)> {
        let index = self.d.u32()?;
        let v = seed.deserialize(IntoDeserializer::<CodecError>::into_deserializer(index))?;
        Ok((v, self))
    }
}

impl<'de> de::VariantAccess<'de> for Variant<'_, 'de> {
    type Error = CodecError;
    fn unit_variant(self) -> Result<()> {
        Ok(())
    }
    fn newtype_variant_seed<T: DeserializeSeed<'de>>(self, seed: T) -> Result<T::Value> {
        seed.deserialize(&mut *self.d)
    }
    fn tuple_variant<V: Visitor<'de>>(self, len: usize, v: V) -> Result<V::Value> {
        v.visit_seq(Elements { d: self.d, left: len })
    }
    fn struct_variant<V: Visitor<'de>>(self, fields: &'static [&'static str], v: V) -> Result<V::Value> {
        v.visit_map(Fields { d: self.d, fields })
    }
}

#[cfg(test)]
mod tests {
    use serde::{Deserialize, Serialize};

    use super::*;

    #[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
    enum Shape {
        Point,
        Circle(f64),
        Box { min: [f64; 3], max: [f64; 3] },
        Pair(u8, i64),
    }

    #[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
    struct Sample {
        name: String,
        labels: Vec<u32>,
        values: Vec<[f32; 3]>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        note: Option<String>,
        #[serde(default, skip_serializing_if = "is_zero")]
        count: usize,
        flags: Vec<bool>,
        shapes: Vec<Shape>,
        nested: Vec<Vec<f64>>,
        pairs: Vec<(usize, [f64; 3])>,
        map: std::collections::BTreeMap<String, i32>,
        maybe: Option<Option<u16>>,
        letter: char,
        unit: (),
    }

    fn is_zero(n: &usize) -> bool {
        *n == 0
    }

    fn sample() -> Sample {
        Sample {
            name: "Flow — étape".into(),
            labels: vec![1, 5, u32::MAX],
            values: vec![[1.0, -0.0, f32::NAN], [f32::INFINITY, 2.5, -3.0]],
            note: Some("kept".into()),
            count: 7,
            flags: vec![true, false],
            shapes: vec![Shape::Point, Shape::Circle(2.0), Shape::Box { min: [0.0; 3], max: [1.0; 3] }, Shape::Pair(9, -4)],
            nested: vec![vec![], vec![1.0, 2.0]],
            pairs: vec![(3, [1.0, 2.0, 3.0])],
            map: [("a".to_string(), -1), ("b".to_string(), 2)].into_iter().collect(),
            maybe: Some(None),
            letter: 'λ',
            unit: (),
        }
    }

    fn same_bits(a: &Sample, b: &Sample) -> bool {
        let bits = |s: &Sample| s.values.iter().flatten().map(|v| v.to_bits()).collect::<Vec<_>>();
        bits(a) == bits(b) && Sample { values: vec![], ..a.clone() } == Sample { values: vec![], ..b.clone() }
    }

    #[test]
    fn values_round_trip_bit_for_bit() {
        let s = sample();
        let bytes = encode(&s).unwrap();
        let back: Sample = decode(&bytes).unwrap();
        assert!(same_bits(&s, &back), "{back:?}");
        assert_eq!(encode(&back).unwrap(), bytes, "the same value, the same bytes");
        // Fields left out read back as their defaults.
        let lean = Sample { note: None, count: 0, ..sample() };
        let back: Sample = decode(&encode(&lean).unwrap()).unwrap();
        assert!(same_bits(&lean, &back));
    }

    #[test]
    fn hostile_bytes_are_refused_without_large_allocations() {
        let bytes = encode(&sample()).unwrap();
        // Every truncation is an error, never a panic.
        for n in 0..bytes.len() {
            assert!(decode::<Sample>(&bytes[..n]).is_err(), "{n}");
        }
        // Bytes left over.
        let mut longer = bytes.clone();
        longer.push(0);
        assert!(decode::<Sample>(&longer).unwrap_err().0.contains("left over"));
        // A count larger than the bytes left is refused before anything is read for it.
        let mut huge = Vec::new();
        huge.extend_from_slice(&u64::MAX.to_le_bytes());
        assert!(decode::<Vec<[f64; 3]>>(&huge).unwrap_err().0.contains("exceeds"));
        let mut many = Vec::new();
        many.extend_from_slice(&(1u64 << 40).to_le_bytes());
        many.extend_from_slice(&[0; 64]);
        assert!(decode::<Vec<u8>>(&many).is_err());
        // Bad tags.
        assert!(decode::<bool>(&[2]).is_err());
        assert!(decode::<Option<u8>>(&[7, 1]).is_err());
        assert!(decode::<Shape>(&99u32.to_le_bytes()).is_err());
        assert!(decode::<char>(&0xD800u32.to_le_bytes()).is_err());
        let mut text = 2u64.to_le_bytes().to_vec();
        text.extend_from_slice(&[0xff, 0xfe]);
        assert!(decode::<String>(&text).is_err());
        // Field tags other than 0 and 1.
        let mut tagged = encode(&sample()).unwrap();
        tagged[0] = 2;
        assert!(decode::<Sample>(&tagged).is_err());
        // Random bytes never panic.
        let mut x: u64 = 0x9e37_79b9_7f4a_7c15;
        for _ in 0..2000 {
            let mut b = bytes.clone();
            for _ in 0..4 {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                let i = (x as usize) % b.len();
                b[i] = (x >> 32) as u8;
            }
            let _ = decode::<Sample>(&b);
        }
    }
}
