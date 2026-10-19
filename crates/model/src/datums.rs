//! Datum planes, axes and points, and the partitions that cut the geometry by datum planes.
//!
//! A datum is a construction reference in model coordinates: a plane (a point and its unit
//! normal), an axis (a point and its unit direction) or a point. Partitions, mirrors and symmetry
//! conditions name planes; points are inputs for planes and axes (their coordinates are taken
//! when those are made). A **partition** cuts bodies of the geometry by a datum plane into cells
//! (the geometry layer does the cutting); the model keeps the partitions in the order they were
//! made, because each cuts the cells the ones before it left.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::geometry::MAX_GEOMETRY_NAME_LEN;
use crate::{Error, Model, Result, check_name, context, insert_new, not_in_use, remove_existing, replace_existing, require};

/// Most partitions a model holds.
pub const MAX_PARTITIONS: usize = 256;
/// Most bodies one partition names.
pub const MAX_PARTITION_BODIES: usize = 1000;

/// A direction shorter than this has none; one this close to unit length is kept bit for bit.
const UNIT_ROUND_OFF: f64 = 1e-12; // unit-factors: the length of a direction vector

/// A datum plane, axis or point.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Datum {
    /// A point on the plane and its unit normal.
    Plane { point: [f64; 3], normal: [f64; 3] },
    /// A point on the axis and its unit direction.
    Axis { point: [f64; 3], direction: [f64; 3] },
    /// A point.
    Point { point: [f64; 3] },
}

impl Datum {
    /// The datum with its point checked and its direction made a unit vector.
    fn checked(self) -> Result<Datum> {
        let unit = |what: &str, v: [f64; 3]| -> Result<[f64; 3]> {
            let len = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            if v.iter().all(|x| x.is_finite()) && len > UNIT_ROUND_OFF && len.is_finite() {
                // A unit vector stays bit for bit, so a model read back equals the one written.
                Ok(if (len - 1.0).abs() <= UNIT_ROUND_OFF { v } else { v.map(|x| x / len) })
            } else {
                Err(Error::Invalid(format!("the {what} must be a finite, non-zero vector")))
            }
        };
        let point = |p: [f64; 3]| -> Result<[f64; 3]> {
            crate::check_coordinates(&p)?;
            Ok(p)
        };
        Ok(match self {
            Datum::Plane { point: p, normal } => Datum::Plane { point: point(p)?, normal: unit("normal", normal)? },
            Datum::Axis { point: p, direction } => Datum::Axis { point: point(p)?, direction: unit("direction", direction)? },
            Datum::Point { point: p } => Datum::Point { point: point(p)? },
        })
    }

    /// The plane's point and unit normal; `None` for an axis or a point.
    pub fn plane(&self) -> Option<([f64; 3], [f64; 3])> {
        match *self {
            Datum::Plane { point, normal } => Some((point, normal)),
            Datum::Axis { .. } | Datum::Point { .. } => None,
        }
    }

    /// The axis's point and unit direction; `None` for a plane or a point.
    pub fn axis(&self) -> Option<([f64; 3], [f64; 3])> {
        match *self {
            Datum::Axis { point, direction } => Some((point, direction)),
            Datum::Plane { .. } | Datum::Point { .. } => None,
        }
    }

    /// The point's coordinates; `None` for a plane or an axis.
    pub fn point(&self) -> Option<[f64; 3]> {
        match *self {
            Datum::Point { point } => Some(point),
            Datum::Plane { .. } | Datum::Axis { .. } => None,
        }
    }

    /// `plane`, `axis` or `point`.
    pub fn kind(&self) -> &'static str {
        match self {
            Datum::Plane { .. } => "plane",
            Datum::Axis { .. } => "axis",
            Datum::Point { .. } => "point",
        }
    }
}

/// A partition: the bodies cut by a datum plane.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Partition {
    pub name: String,
    /// The datum plane that cuts.
    pub datum: String,
    /// The bodies it cuts, as they are when it comes (cells of earlier partitions included);
    /// empty: every body the plane divides.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub bodies: Vec<String>,
    /// Kept in the model but left out of it: the design is cut without it.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub suppressed: bool,
}

impl Model {
    pub fn datums(&self) -> &BTreeMap<String, Datum> {
        &self.datums
    }

    pub fn datum(&self, name: &str) -> Option<&Datum> {
        self.datums.get(name)
    }

    pub fn add_datum(&mut self, name: &str, datum: Datum) -> Result<()> {
        let datum = datum.checked().map_err(|e| context("datum", name, e))?;
        insert_new(&mut self.datums, "datum", name, datum)
    }

    /// Replace a datum; a datum plane a partition uses stays a plane.
    pub fn replace_datum(&mut self, name: &str, datum: Datum) -> Result<()> {
        let datum = datum.checked().map_err(|e| context("datum", name, e))?;
        if datum.plane().is_none() {
            not_in_use("datum", name, self.datum_users(name))?;
        }
        replace_existing(&mut self.datums, "datum", name, datum)
    }

    pub fn remove_datum(&mut self, name: &str) -> Result<()> {
        not_in_use("datum", name, self.datum_users(name))?;
        remove_existing(&mut self.datums, "datum", name)?;
        self.suppressed_datums.remove(name);
        Ok(())
    }

    /// Whether datum `name` is suppressed: kept, but not usable as an input and not cutting.
    pub fn datum_suppressed(&self, name: &str) -> bool {
        self.suppressed_datums.contains(name)
    }

    /// Suppress datum `name` (or bring it back).
    pub fn set_datum_suppressed(&mut self, name: &str, suppressed: bool) -> Result<()> {
        require(&self.datums, "datum", name)?;
        if suppressed {
            self.suppressed_datums.insert(name.to_string());
        } else {
            self.suppressed_datums.remove(name);
        }
        Ok(())
    }

    /// Datum `name` as an input to another definition: refused while it is suppressed.
    pub fn usable_datum(&self, name: &str) -> Result<&Datum> {
        let d = self.datums.get(name).ok_or_else(|| Error::NotFound(format!("datum `{name}`")))?;
        if self.datum_suppressed(name) {
            return Err(Error::Invalid(format!("datum `{name}` is suppressed: unsuppress it to use it")));
        }
        Ok(d)
    }

    /// Suppress partition `name` (or bring it back): the design is cut without it.
    pub fn set_partition_suppressed(&mut self, name: &str, suppressed: bool) -> Result<()> {
        let p = self.partitions.iter_mut().find(|p| p.name == name).ok_or_else(|| Error::NotFound(format!("partition `{name}`")))?;
        p.suppressed = suppressed;
        Ok(())
    }

    /// Whether partition `p` cuts: neither it nor its datum is suppressed.
    pub fn partition_cuts(&self, p: &Partition) -> bool {
        !p.suppressed && !self.datum_suppressed(&p.datum)
    }

    /// The definitions that refer to a datum (the partitions that cut by it), as `kind \`name\``
    /// strings.
    pub fn datum_users(&self, name: &str) -> Vec<String> {
        let partitions = self.partitions.iter().filter(|p| p.datum == name).map(|p| format!("partition `{}`", p.name));
        partitions.chain(self.face_splits().iter().filter(|s| s.datum == name).map(|s| format!("face split `{}`", s.name))).collect()
    }

    /// The partitions, in the order they cut.
    pub fn partitions(&self) -> &[Partition] {
        &self.partitions
    }

    /// Add a partition after the others.
    pub fn add_partition(&mut self, partition: Partition) -> Result<()> {
        let name = partition.name.clone();
        check_name("partition", &name)?;
        if self.partitions.iter().any(|p| p.name == name) {
            return Err(Error::Invalid(format!("partition `{name}` already exists")));
        }
        if self.partitions.len() >= MAX_PARTITIONS {
            return Err(Error::Limit(format!("a model holds at most {MAX_PARTITIONS} partitions")));
        }
        let invalid = |m: String| context("partition", &name, Error::Invalid(m));
        match self.datums.get(&partition.datum) {
            None => return Err(context("partition", &name, Error::NotFound(format!("datum `{}`", partition.datum)))),
            Some(d) if d.plane().is_none() => {
                return Err(invalid(format!("datum `{}` is a datum {}; a partition cuts by a plane", partition.datum, d.kind())));
            }
            Some(_) => {}
        }
        if partition.bodies.len() > MAX_PARTITION_BODIES {
            return Err(invalid(format!("at most {MAX_PARTITION_BODIES} bodies")));
        }
        for (i, b) in partition.bodies.iter().enumerate() {
            if b.trim().is_empty() || b.len() > MAX_GEOMETRY_NAME_LEN {
                return Err(invalid(format!("body names are 1 to {MAX_GEOMETRY_NAME_LEN} bytes")));
            }
            if partition.bodies.get(..i).is_some_and(|before| before.contains(b)) {
                return Err(invalid(format!("body `{b}` is listed twice")));
            }
        }
        self.partitions.push(partition);
        Ok(())
    }

    pub fn remove_partition(&mut self, name: &str) -> Result<Partition> {
        let i = self.partitions.iter().position(|p| p.name == name).ok_or_else(|| Error::NotFound(format!("partition `{name}`")))?;
        Ok(self.partitions.remove(i))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plane(z: f64) -> Datum {
        Datum::Plane { point: [0.0, 0.0, z], normal: [0.0, 0.0, 2.0] }
    }

    /// A suppressed datum is kept but not usable, and stops its partitions cutting; a suppressed
    /// partition stops cutting; both are saved (#291).
    #[test]
    fn suppressed_datums_and_partitions_are_kept_but_left_out() {
        let mut m = Model::default();
        m.add_datum("Mid", plane(5.0)).unwrap();
        m.add_datum("Low", plane(1.0)).unwrap();
        let p = |name: &str, datum: &str| Partition { name: name.into(), datum: datum.into(), bodies: vec![], suppressed: false };
        m.add_partition(p("Half", "Mid")).unwrap();
        m.add_partition(p("Base", "Low")).unwrap();
        assert!(m.partitions().iter().all(|p| m.partition_cuts(p)));
        m.set_datum_suppressed("Mid", true).unwrap();
        m.set_partition_suppressed("Base", true).unwrap();
        assert!(m.partitions().iter().all(|p| !m.partition_cuts(p)));
        assert!(m.usable_datum("Mid").unwrap_err().to_string().contains("suppressed"));
        assert!(m.usable_datum("Low").is_ok() && m.datum("Mid").is_some());
        assert!(m.set_datum_suppressed("Nope", true).is_err() && m.set_partition_suppressed("Nope", true).is_err());
        let json = serde_json::to_string(&m).unwrap();
        let back: Model = serde_json::from_str(&json).unwrap();
        assert_eq!(back, m);
        assert!(back.datum_suppressed("Mid") && back.partitions().iter().any(|p| p.suppressed));
        m.set_datum_suppressed("Mid", false).unwrap();
        assert!(m.usable_datum("Mid").is_ok());
        let plain = serde_json::to_string(&Model::default()).unwrap();
        assert!(!plain.contains("suppressed"), "nothing suppressed, nothing written");
    }

    #[test]
    fn datums_are_normalised_and_partitions_refer_to_planes() {
        let mut m = Model::default();
        m.add_datum("Mid", plane(5.0)).unwrap();
        assert_eq!(m.datum("Mid"), Some(&Datum::Plane { point: [0.0, 0.0, 5.0], normal: [0.0, 0.0, 1.0] }));
        m.add_datum("Axis", Datum::Axis { point: [0.0; 3], direction: [1.0, 1.0, 0.0] }).unwrap();
        let Some(Datum::Axis { direction, .. }) = m.datum("Axis") else { panic!() };
        assert!((direction[0] - 0.5f64.sqrt()).abs() < 1e-15);

        let p = |name: &str, datum: &str| Partition { name: name.into(), datum: datum.into(), bodies: vec![], suppressed: false };
        m.add_partition(p("Half", "Mid")).unwrap();
        assert!(m.add_partition(p("Half", "Mid")).is_err());
        assert!(m.add_partition(p("Bad", "Axis")).is_err());
        m.add_datum("P", Datum::Point { point: [1.0, 2.0, 3.0] }).unwrap();
        assert_eq!((m.datum("P").and_then(Datum::point), m.datum("P").map(Datum::kind)), (Some([1.0, 2.0, 3.0]), Some("point")));
        assert!(m.add_partition(p("Bad", "P")).is_err());
        assert!(m.add_datum("Q", Datum::Point { point: [f64::NAN, 0.0, 0.0] }).is_err());
        assert_eq!(m.datum_users("Mid"), ["partition `Half`"]);
        assert!(m.add_partition(p("Bad", "Nope")).is_err());
        assert!(m.add_partition(Partition { bodies: vec!["A".into(), "A".into()], ..p("Bad", "Mid") }).is_err());
        assert!(matches!(m.remove_datum("Mid"), Err(Error::InUse(_))));
        assert!(matches!(m.replace_datum("Mid", Datum::Axis { point: [0.0; 3], direction: [1.0, 0.0, 0.0] }), Err(Error::InUse(_))));
        m.replace_datum("Mid", plane(3.0)).unwrap();
        assert_eq!(m.partitions().len(), 1);

        // The model round-trips through JSON, partitions in order.
        m.add_partition(Partition { bodies: vec!["Body1/a".into()], ..p("Quarter", "Mid") }).unwrap();
        let back: Model = serde_json::from_value(serde_json::to_value(&m).unwrap()).unwrap();
        assert_eq!(back, m);
        m.remove_partition("Half").unwrap();
        m.remove_partition("Quarter").unwrap();
        m.remove_datum("Mid").unwrap();
    }

    #[test]
    fn hostile_datums_are_an_error() {
        let mut m = Model::default();
        for d in [
            Datum::Plane { point: [0.0; 3], normal: [0.0; 3] },
            Datum::Plane { point: [f64::NAN, 0.0, 0.0], normal: [0.0, 0.0, 1.0] },
            Datum::Plane { point: [1e300, 0.0, 0.0], normal: [0.0, 0.0, 1.0] },
            Datum::Axis { point: [0.0; 3], direction: [f64::INFINITY, 0.0, 0.0] },
            Datum::Axis { point: [0.0; 3], direction: [1e-300, 0.0, 0.0] },
        ] {
            assert!(m.add_datum("D", d).is_err(), "{d:?}");
        }
        assert!(m.add_datum("1bad", plane(0.0)).is_err());
        assert!(m.datums().is_empty());
    }
}
