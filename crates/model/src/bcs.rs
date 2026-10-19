//! Boundary conditions: prescribed degrees of freedom on a set's nodes or a reference point.
//!
//! Degrees of freedom are numbered 1–3 for translations along X, Y, Z and 4–6 for rotations about
//! them, in the global axes. Which steps a boundary condition is active in is listed explicitly by
//! each step ([`crate::steps`]).

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::{Error, Model, Result, context, finite, insert_new, not_in_use, remove_existing, replace_existing, require};

/// A set of degrees of freedom (1–6), written as a sorted list such as `[1, 2, 3]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(try_from = "Vec<u8>", into = "Vec<u8>")]
pub struct Dofs(u8);

impl Dofs {
    pub fn all() -> Dofs {
        Dofs(0b111_111)
    }

    pub fn none() -> Dofs {
        Dofs(0)
    }

    pub fn contains(self, dof: u8) -> bool {
        (1..=6).contains(&dof) && self.0 & (1 << (dof - 1)) != 0
    }

    pub fn is_empty(self) -> bool {
        self.0 == 0
    }

    /// The degrees of freedom as bits: bit `d - 1` for degree of freedom `d`.
    pub fn bits(self) -> u8 {
        self.0
    }

    /// The degrees of freedom of bits (bits above the sixth are ignored).
    pub fn from_bits(bits: u8) -> Dofs {
        Dofs(bits & 0b111_111)
    }

    /// The degrees of freedom in order.
    pub fn list(self) -> Vec<u8> {
        (1..=6).filter(|d| self.contains(*d)).collect()
    }
}

impl TryFrom<Vec<u8>> for Dofs {
    type Error = Error;
    fn try_from(v: Vec<u8>) -> Result<Dofs> {
        let mut mask = 0u8;
        for d in v {
            if !(1..=6).contains(&d) {
                return Err(Error::Invalid(format!("degree of freedom {d} is not in 1..=6")));
            }
            mask |= 1 << (d - 1);
        }
        Ok(Dofs(mask))
    }
}

impl From<Dofs> for Vec<u8> {
    fn from(d: Dofs) -> Vec<u8> {
        d.list()
    }
}

/// What a boundary condition, force or moment acts on.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Target {
    /// The nodes of a set.
    Set(String),
    /// A reference point.
    ReferencePoint(String),
}

impl Target {
    pub fn set(&self) -> Option<&str> {
        match self {
            Target::Set(s) => Some(s),
            Target::ReferencePoint(_) => None,
        }
    }

    pub fn reference_point(&self) -> Option<&str> {
        match self {
            Target::ReferencePoint(r) => Some(r),
            Target::Set(_) => None,
        }
    }

    pub(crate) fn check_exists(&self, model: &Model) -> Result<()> {
        match self {
            Target::Set(s) => require(&model.sets, "set", s),
            Target::ReferencePoint(r) => require(&model.reference_points, "reference point", r),
        }
    }
}

impl std::fmt::Display for Target {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Target::Set(s) => write!(f, "set `{s}`"),
            Target::ReferencePoint(r) => write!(f, "reference point `{r}`"),
        }
    }
}

/// The plane of a symmetry boundary condition, named by its normal: an axis (`"x"`), or any
/// direction (`{"normal": [1, 1, 0]}`, not necessarily of unit length).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Symmetry {
    /// The plane normal to X: no X translation, no rotation about Y or Z.
    X,
    Y,
    Z,
    /// A plane normal to a direction: no translation along it, no rotation about the two axes
    /// in the plane (n · u = 0, and θ ∥ n at nodes with rotations).
    Normal([f64; 3]),
}

/// A normal whose component along an axis is within this of 1 (unit length) is that axis.
const ALONG_AXIS: f64 = 1e-12; // unit-factors: a component of a unit vector

impl Symmetry {
    /// The plane normal to `normal`: [`Symmetry::X`], `Y` or `Z` when it lies along an axis,
    /// otherwise [`Symmetry::Normal`] with the unit normal; `None` for a zero or non-finite one.
    pub fn from_normal(normal: [f64; 3]) -> Option<Symmetry> {
        let unit = unit(normal)?;
        Some(match (0..3).find(|&k| (unit[k].abs() - 1.0).abs() < ALONG_AXIS) {
            Some(0) => Symmetry::X,
            Some(1) => Symmetry::Y,
            Some(2) => Symmetry::Z,
            _ => Symmetry::Normal(unit),
        })
    }

    /// The unit normal (`None` for a zero or non-finite [`Symmetry::Normal`]).
    pub fn unit_normal(self) -> Option<[f64; 3]> {
        match self {
            Symmetry::X => Some([1.0, 0.0, 0.0]),
            Symmetry::Y => Some([0.0, 1.0, 0.0]),
            Symmetry::Z => Some([0.0, 0.0, 1.0]),
            Symmetry::Normal(n) => unit(n),
        }
    }

    /// The axis (0–2) the normal lies along, so the condition holds global degrees of freedom;
    /// `None` for an inclined plane, which is a linear constraint instead.
    pub fn axis(self) -> Option<usize> {
        let n = self.unit_normal()?;
        (0..3).find(|&k| (n[k].abs() - 1.0).abs() < ALONG_AXIS)
    }

    /// The held global degrees of freedom: those of an axis-aligned plane; none for an inclined
    /// one, which holds combinations of them ([`Symmetry::unit_normal`]).
    pub fn dofs(self) -> Dofs {
        match self.axis() {
            Some(0) => Dofs(0b110_001),
            Some(1) => Dofs(0b101_010),
            Some(2) => Dofs(0b011_100),
            _ => Dofs::none(),
        }
    }

    /// The inclined plane's unit normal; `None` for one normal to an axis.
    pub fn inclined(self) -> Option<[f64; 3]> {
        if self.axis().is_some() { None } else { self.unit_normal() }
    }

    /// Two unit vectors in the plane, perpendicular to each other and to the unit normal n:
    /// t₁ = n × e / |n × e| with e the axis n leans on least, and t₂ = n × t₁ (`None` for a
    /// zero or non-finite normal).
    pub fn in_plane(self) -> Option<[[f64; 3]; 2]> {
        let n = self.unit_normal()?;
        let least = (0..3).fold(0, |m, k| if n[k].abs() < n[m].abs() { k } else { m });
        let mut e = [0.0; 3];
        e[least] = 1.0;
        let cross = |a: [f64; 3], b: [f64; 3]| [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
        let t1 = unit(cross(n, e))?;
        Some([t1, cross(n, t1)])
    }

    fn validate(self) -> Result<()> {
        if let Symmetry::Normal(n) = self
            && unit(n).is_none()
        {
            return Err(Error::Invalid(format!("the symmetry plane's normal must be a finite, non-zero vector, got {n:?}")));
        }
        Ok(())
    }
}

/// `v` scaled to unit length; `None` for a zero or non-finite vector.
fn unit(v: [f64; 3]) -> Option<[f64; 3]> {
    let scale = v.iter().fold(0.0f64, |m, c| m.max(c.abs()));
    if !(scale > 0.0 && v.iter().all(|c| c.is_finite())) {
        return None;
    }
    let w = v.map(|c| c / scale);
    let length = w.iter().map(|c| c * c).sum::<f64>().sqrt();
    Some(w.map(|c| c / length))
}

/// A boundary condition.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum BoundaryCondition {
    /// Every degree of freedom held at zero.
    Fixed { target: Target },
    /// Prescribed values for some degrees of freedom: `values[i]` is degree of freedom `i + 1`
    /// (`null` leaves it free, `0` holds it).
    Displacement { target: Target, values: [Option<f64>; 6] },
    /// Mirror symmetry about a plane through the nodes.
    Symmetry { target: Target, plane: Symmetry },
    /// A prescribed temperature (degree of freedom [`TEMPERATURE_DOF`]) at the target's nodes,
    /// in heat transfer and coupled temperature–displacement steps.
    Temperature { target: Target, value: f64 },
}

/// The model's number of the temperature degree of freedom (the deck formats' convention).
pub const TEMPERATURE_DOF: u8 = 11;

impl BoundaryCondition {
    pub fn target(&self) -> &Target {
        match self {
            BoundaryCondition::Fixed { target }
            | BoundaryCondition::Displacement { target, .. }
            | BoundaryCondition::Symmetry { target, .. }
            | BoundaryCondition::Temperature { target, .. } => target,
        }
    }

    /// The constrained structural degrees of freedom (1–6); none for a temperature.
    pub fn dofs(&self) -> Dofs {
        match self {
            BoundaryCondition::Temperature { .. } => Dofs::none(),
            BoundaryCondition::Fixed { .. } => Dofs::all(),
            BoundaryCondition::Displacement { values, .. } => {
                Dofs(values.iter().enumerate().filter(|(_, v)| v.is_some()).fold(0, |m, (i, _)| m | (1 << i)))
            }
            BoundaryCondition::Symmetry { plane, .. } => plane.dofs(),
        }
    }

    /// Whether it drives the model on its own, as a load does: a non-zero prescribed
    /// displacement or a non-zero prescribed temperature (a support only holds).
    pub fn drives(&self) -> bool {
        match self {
            BoundaryCondition::Displacement { values, .. } => values.iter().flatten().any(|v| *v != 0.0),
            BoundaryCondition::Temperature { value, .. } => *value != 0.0,
            BoundaryCondition::Fixed { .. } | BoundaryCondition::Symmetry { .. } => false,
        }
    }

    /// Whether it prescribes temperatures (rather than motion).
    pub fn is_thermal(&self) -> bool {
        matches!(self, BoundaryCondition::Temperature { .. })
    }

    /// The prescribed value of a constrained degree of freedom (1–6, or [`TEMPERATURE_DOF`]);
    /// `None` when it is free.
    pub fn value(&self, dof: u8) -> Option<f64> {
        match self {
            BoundaryCondition::Temperature { value, .. } => (dof == TEMPERATURE_DOF).then_some(*value),
            BoundaryCondition::Displacement { values, .. } => values.get(usize::from(dof).checked_sub(1)?).copied().flatten(),
            _ => self.dofs().contains(dof).then_some(0.0),
        }
    }

    pub fn validate(&self) -> Result<()> {
        if let BoundaryCondition::Symmetry { plane, .. } = self {
            plane.validate()?;
        }
        if let BoundaryCondition::Temperature { target, value } = self {
            finite("temperature", *value)?;
            if let Target::ReferencePoint(r) = target {
                return Err(Error::Invalid(format!("reference points have no temperature (`{r}`); prescribe it on a set")));
            }
        }
        if let BoundaryCondition::Displacement { values, .. } = self {
            for (i, v) in values.iter().enumerate() {
                if let Some(v) = v {
                    finite(&format!("degree of freedom {}", i + 1), *v)?;
                }
            }
            if values.iter().all(Option::is_none) {
                return Err(Error::Invalid("a displacement boundary condition needs at least one value".into()));
            }
        }
        Ok(())
    }
}

impl Model {
    pub fn boundary_conditions(&self) -> &BTreeMap<String, BoundaryCondition> {
        &self.boundary_conditions
    }

    pub fn boundary_condition(&self, name: &str) -> Option<&BoundaryCondition> {
        self.boundary_conditions.get(name)
    }

    fn check_bc(&self, name: &str, bc: &BoundaryCondition) -> Result<()> {
        bc.validate().map_err(|e| context("boundary condition", name, e))?;
        bc.target().check_exists(self).map_err(|e| context("boundary condition", name, e))
    }

    pub fn add_boundary_condition(&mut self, name: &str, bc: BoundaryCondition) -> Result<()> {
        self.check_bc(name, &bc)?;
        insert_new(&mut self.boundary_conditions, "boundary condition", name, bc)
    }

    pub fn replace_boundary_condition(&mut self, name: &str, bc: BoundaryCondition) -> Result<()> {
        self.check_bc(name, &bc)?;
        replace_existing(&mut self.boundary_conditions, "boundary condition", name, bc)
    }

    /// The steps that list a boundary condition, as `step \`name\`` strings.
    pub fn boundary_condition_users(&self, name: &str) -> Vec<String> {
        self.steps.iter().filter(|s| s.boundary_conditions.contains(name)).map(|s| format!("step `{}`", s.name)).collect()
    }

    /// Remove a boundary condition that no step lists.
    pub fn remove_boundary_condition(&mut self, name: &str) -> Result<()> {
        not_in_use("boundary condition", name, self.boundary_condition_users(name))?;
        remove_existing(&mut self.boundary_conditions, "boundary condition", name).map(|_| ())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dofs_and_values() {
        let d: Dofs = serde_json::from_value(serde_json::json!([3, 1, 1])).unwrap();
        assert_eq!(d.list(), vec![1, 3]);
        assert_eq!(serde_json::to_value(d).unwrap(), serde_json::json!([1, 3]));
        assert!(serde_json::from_value::<Dofs>(serde_json::json!([0])).is_err());
        assert!(serde_json::from_value::<Dofs>(serde_json::json!([7])).is_err());
        assert_eq!(Symmetry::X.dofs().list(), vec![1, 5, 6]);
        assert_eq!(Symmetry::Y.dofs().list(), vec![2, 4, 6]);
        assert_eq!(Symmetry::Z.dofs().list(), vec![3, 4, 5]);

        let t = Target::Set("s".into());
        let bc = BoundaryCondition::Displacement { target: t.clone(), values: [Some(0.0), None, Some(1.5), None, None, None] };
        assert_eq!(bc.dofs().list(), vec![1, 3]);
        assert_eq!((bc.value(1), bc.value(2), bc.value(3), bc.value(0), bc.value(9)), (Some(0.0), None, Some(1.5), None, None));
        assert_eq!(BoundaryCondition::Fixed { target: t.clone() }.value(6), Some(0.0));
        assert!(BoundaryCondition::Displacement { target: t.clone(), values: [None; 6] }.validate().is_err());
        assert!(BoundaryCondition::Displacement { target: t, values: [Some(f64::NAN), None, None, None, None, None] }.validate().is_err());
        let hot = BoundaryCondition::Temperature { target: Target::Set("s".into()), value: 80.0 };
        assert_eq!((hot.value(TEMPERATURE_DOF), hot.value(1), hot.dofs().list()), (Some(80.0), None, vec![]));
        assert!(hot.is_thermal());
        assert!(BoundaryCondition::Temperature { target: Target::Set("s".into()), value: f64::INFINITY }.validate().is_err());
        assert!(BoundaryCondition::Temperature { target: Target::ReferencePoint("RP".into()), value: 1.0 }.validate().is_err());
        let json = serde_json::json!({"type": "symmetry", "target": {"reference_point": "RP"}, "plane": "z"});
        let parsed: BoundaryCondition = serde_json::from_value(json.clone()).unwrap();
        assert_eq!(serde_json::to_value(&parsed).unwrap(), json);
    }

    #[test]
    fn inclined_symmetry_planes() {
        let tilt = Symmetry::from_normal([1.0, 1.0, 0.0]).unwrap();
        let h = 1.0 / 2f64.sqrt();
        assert_eq!(tilt, Symmetry::Normal([h, h, 0.0]));
        assert_eq!((tilt.axis(), tilt.dofs(), tilt.inclined()), (None, Dofs::none(), Some([h, h, 0.0])));
        assert_eq!(Symmetry::from_normal([0.0, -2.0, 0.0]), Some(Symmetry::Y));
        assert_eq!(Symmetry::Normal([0.0, 0.0, 3.0]).dofs(), Symmetry::Z.dofs());
        assert_eq!(Symmetry::Normal([0.0, 0.0, 3.0]).inclined(), None);
        for bad in [[0.0; 3], [f64::NAN, 1.0, 0.0], [f64::INFINITY, 0.0, 0.0]] {
            assert_eq!(Symmetry::from_normal(bad), None);
            let bc = BoundaryCondition::Symmetry { target: Target::Set("s".into()), plane: Symmetry::Normal(bad) };
            assert!(bc.validate().is_err(), "{bad:?}");
        }
        // The tiniest and the largest finite normals still give a unit one.
        assert_eq!(Symmetry::from_normal([1e-300, 1e-300, 0.0]), Some(Symmetry::Normal([h, h, 0.0])));
        assert_eq!(Symmetry::from_normal([f64::MAX, f64::MAX, 0.0]), Some(Symmetry::Normal([h, h, 0.0])));
        for plane in [tilt, Symmetry::X, Symmetry::Normal([0.3, -2.0, 7.0])] {
            let n = plane.unit_normal().unwrap();
            let [a, b] = plane.in_plane().unwrap();
            let dot = |u: [f64; 3], v: [f64; 3]| (0..3).map(|k| u[k] * v[k]).sum::<f64>();
            for (u, v, want) in [(a, a, 1.0), (b, b, 1.0), (a, b, 0.0), (a, n, 0.0), (b, n, 0.0)] {
                assert!((dot(u, v) - want).abs() < 1e-15, "{plane:?}");
            }
        }
        let json = serde_json::json!({"type": "symmetry", "target": {"set": "s"}, "plane": {"normal": [1.0, 1.0, 0.0]}});
        let parsed: BoundaryCondition = serde_json::from_value(json.clone()).unwrap();
        assert_eq!(parsed, BoundaryCondition::Symmetry { target: Target::Set("s".into()), plane: Symmetry::Normal([1.0, 1.0, 0.0]) });
        assert_eq!(serde_json::to_value(&parsed).unwrap(), json);
    }
}
