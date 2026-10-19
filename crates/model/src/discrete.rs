//! Discrete elements: point masses with rotary inertia, springs and dashpots.
//!
//! They are definitions, not mesh elements: each acts on a target (the nodes of a set, or a
//! reference point) and survives remeshing like a load does. Every node of the target gets its
//! own copy, with the full value: a point mass of 2 on a set of five nodes adds 2 at each of them
//! (a mass shared over a face is a point mass on a reference point coupled to the face).
//!
//! - A **point mass** adds a translational mass, the same in every direction or one per global
//!   direction, and a rotary inertia tensor about the node (which needs rotational degrees of
//!   freedom: a reference point, or a shell or beam node).
//! - A **spring** resists the relative displacement of its two ends along a direction with a
//!   force `stiffness × relative displacement`; a **dashpot** resists the relative velocity with
//!   `coefficient × relative velocity`. With no second end, the spring or dashpot connects each
//!   node of the target to the ground. The direction is a global degree of freedom (1–3
//!   translations, 4–6 rotations: then a torsional spring or dashpot), a fixed vector, or the
//!   line through the two ends (two-point only, from the undeformed positions).
//!
//! The linear (small displacement) form is the textbook spring element: for a unit direction n,
//! the element matrix is k [n nᵀ, −n nᵀ; −n nᵀ, n nᵀ] on the two ends' translations (R. D. Cook,
//! D. S. Malkus, M. E. Plesha, R. J. Witt, *Concepts and Applications of Finite Element
//! Analysis*, 4th ed., Wiley 2002, §2.2 and §2.10), and k n nᵀ for a grounded one. Dashpots
//! have the same form on velocities; static and frequency steps ignore them.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::bcs::Target;
use crate::{Error, Model, Result, context, finite, insert_new, remove_existing, replace_existing};

/// Most copies of discrete elements a model resolves to (definitions × target nodes).
pub const MAX_DISCRETE: usize = 10_000_000;

/// A principal minor of an inertia above minus this fraction of its largest component is
/// non-negative (to round-off).
const INERTIA_ROUND_OFF: f64 = 1e-12; // unit-factors: relative to the largest component

/// A translational mass: the same in every direction, or one per global direction.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum Mass {
    Isotropic(f64),
    Directional([f64; 3]),
}

impl Mass {
    /// The mass in X, Y and Z.
    pub fn components(self) -> [f64; 3] {
        match self {
            Mass::Isotropic(m) => [m; 3],
            Mass::Directional(m) => m,
        }
    }
}

/// A point mass and rotary inertia at every node of a target.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PointMass {
    pub target: Target,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mass: Option<Mass>,
    /// The inertia tensor about the node in global axes, `[I11, I22, I33, I12, I13, I23]`: the
    /// components of the symmetric tensor J with angular momentum H = J ω (products of inertia
    /// carry their sign: I12 = −∫ x y dm).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inertia: Option<[f64; 6]>,
}

impl PointMass {
    pub fn validate(&self) -> Result<()> {
        if self.mass.is_none() && self.inertia.is_none() {
            return Err(Error::Invalid("a point mass needs a mass, an inertia or both".into()));
        }
        if let Some(m) = self.mass {
            for c in m.components() {
                if finite("mass", c)? < 0.0 {
                    return Err(Error::Invalid(format!("a mass is not negative, got {c}")));
                }
            }
            if m.components().iter().all(|c| *c == 0.0) && self.inertia.is_none() {
                return Err(Error::Invalid("a point mass of zero adds nothing".into()));
            }
        }
        if let Some(i) = self.inertia {
            for c in i {
                finite("inertia", c)?;
            }
            // Positive semidefinite: every principal minor is non-negative (to round-off).
            let j = inertia_matrix(i);
            let scale = i.iter().fold(0.0f64, |m, c| m.max(c.abs()));
            let tolerance = INERTIA_ROUND_OFF * scale;
            let minors = [
                j[0][0],
                j[1][1],
                j[2][2],
                j[0][0] * j[1][1] - j[0][1] * j[0][1],
                j[0][0] * j[2][2] - j[0][2] * j[0][2],
                j[1][1] * j[2][2] - j[1][2] * j[1][2],
                j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1]) - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
                    + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]),
            ];
            let powers = [1, 1, 1, 2, 2, 2, 3];
            if minors.iter().zip(powers).any(|(m, p)| *m < -tolerance * scale.powi(p - 1)) {
                return Err(Error::Invalid(
                    "an inertia tensor must be positive semidefinite (non-negative moments, products no larger than the moments allow)".into(),
                ));
            }
            if i.iter().all(|c| *c == 0.0) && self.mass.is_none_or(|m| m.components().iter().all(|c| *c == 0.0)) {
                return Err(Error::Invalid("a point mass of zero adds nothing".into()));
            }
        }
        Ok(())
    }
}

/// `[I11, I22, I33, I12, I13, I23]` as a symmetric matrix.
pub fn inertia_matrix(i: [f64; 6]) -> [[f64; 3]; 3] {
    let [i11, i22, i33, i12, i13, i23] = i;
    [[i11, i12, i13], [i12, i22, i23], [i13, i23, i33]]
}

/// The direction a spring or dashpot acts in.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Direction {
    /// A global degree of freedom: 1–3 translations along X, Y, Z, 4–6 rotations about them.
    Dof(u8),
    /// A fixed direction in global axes (translations; normalized when used).
    Vector([f64; 3]),
    /// The line through the two ends' undeformed positions (two-point only).
    Line,
}

impl Direction {
    /// Whether it acts on rotations.
    pub fn is_rotational(self) -> bool {
        matches!(self, Direction::Dof(4..=6))
    }

    fn validate(self, two_point: bool) -> Result<()> {
        match self {
            Direction::Dof(d) if !(1..=6).contains(&d) => Err(Error::Invalid(format!("degree of freedom {d} is not in 1..=6"))),
            Direction::Dof(_) => Ok(()),
            Direction::Vector(v) => {
                for c in v {
                    finite("direction", c)?;
                }
                if v.iter().all(|c| *c == 0.0) {
                    return Err(Error::Invalid("a direction vector must not be zero".into()));
                }
                Ok(())
            }
            Direction::Line if !two_point => Err(Error::Invalid("the direction `line` needs a second end (it runs between the two ends)".into())),
            Direction::Line => Ok(()),
        }
    }
}

/// A linear spring or dashpot between the nodes of `first` and the ground, or between `first`
/// and `second` (each one node: a vertex set or a reference point).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Spring {
    /// Force (or moment) = stiffness × relative displacement (or rotation).
    Spring {
        first: Target,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        second: Option<Target>,
        direction: Direction,
        stiffness: f64,
    },
    /// Force (or moment) = coefficient × relative velocity (or angular velocity); dynamic steps
    /// only.
    Dashpot {
        first: Target,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        second: Option<Target>,
        direction: Direction,
        coefficient: f64,
    },
}

impl Spring {
    pub fn first(&self) -> &Target {
        match self {
            Spring::Spring { first, .. } | Spring::Dashpot { first, .. } => first,
        }
    }

    pub fn second(&self) -> Option<&Target> {
        match self {
            Spring::Spring { second, .. } | Spring::Dashpot { second, .. } => second.as_ref(),
        }
    }

    pub fn direction(&self) -> Direction {
        match self {
            Spring::Spring { direction, .. } | Spring::Dashpot { direction, .. } => *direction,
        }
    }

    /// The stiffness or the damping coefficient.
    pub fn value(&self) -> f64 {
        match self {
            Spring::Spring { stiffness, .. } => *stiffness,
            Spring::Dashpot { coefficient, .. } => *coefficient,
        }
    }

    pub fn is_dashpot(&self) -> bool {
        matches!(self, Spring::Dashpot { .. })
    }

    /// `spring` or `dashpot`.
    pub fn kind(&self) -> &'static str {
        if self.is_dashpot() { "dashpot" } else { "spring" }
    }

    /// The targets it acts on.
    pub fn targets(&self) -> impl Iterator<Item = &Target> {
        std::iter::once(self.first()).chain(self.second())
    }

    pub fn validate(&self) -> Result<()> {
        let what = if self.is_dashpot() { "damping coefficient" } else { "stiffness" };
        let v = finite(what, self.value())?;
        if v <= 0.0 {
            return Err(Error::Invalid(format!("a {} needs a positive {what}, got {v}", self.kind())));
        }
        self.direction().validate(self.second().is_some())?;
        if self.second() == Some(self.first()) {
            return Err(Error::Invalid(format!("a {} cannot connect a target to itself", self.kind())));
        }
        Ok(())
    }
}

/// A resolved direction: a global degree of freedom (1–6), or a unit vector acting on the
/// translations (from a fixed vector, or along the line of the ends).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ResolvedDirection {
    Dof(u8),
    Vector([f64; 3]),
    Line([f64; 3]),
}

impl ResolvedDirection {
    /// The coefficients of the relative motion on each end's six degrees of freedom: the
    /// spring's elongation is Σ c (u_second − u_first).
    pub fn coefficients(self) -> [f64; 6] {
        let mut c = [0.0; 6];
        match self {
            ResolvedDirection::Dof(d) => {
                if let Some(x) = usize::from(d).checked_sub(1).and_then(|i| c.get_mut(i)) {
                    *x = 1.0;
                }
            }
            ResolvedDirection::Vector(n) | ResolvedDirection::Line(n) => c[..3].copy_from_slice(&n),
        }
        c
    }

    pub fn is_rotational(self) -> bool {
        matches!(self, ResolvedDirection::Dof(4..=6))
    }

    /// The unit vector it acts along (for a rotation, the axis it turns about).
    pub fn axis(self) -> [f64; 3] {
        let c = self.coefficients();
        if self.is_rotational() { [c[3], c[4], c[5]] } else { [c[0], c[1], c[2]] }
    }
}

/// One point mass at one node.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedMass {
    /// The definition it comes from.
    pub name: String,
    pub node: u32,
    /// Mass in X, Y, Z.
    pub mass: [f64; 3],
    /// `[I11, I22, I33, I12, I13, I23]` (zero when there is none).
    pub inertia: [f64; 6],
}

impl ResolvedSpring {
    /// `spring` or `dashpot`.
    pub fn kind(&self) -> &'static str {
        if self.dashpot { "dashpot" } else { "spring" }
    }
}

/// One spring or dashpot at one node (grounded) or between two nodes.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedSpring {
    pub name: String,
    pub dashpot: bool,
    pub first: u32,
    /// `None`: to the ground.
    pub second: Option<u32>,
    pub direction: ResolvedDirection,
    /// Stiffness, or damping coefficient.
    pub value: f64,
}

impl Model {
    pub fn masses(&self) -> &BTreeMap<String, PointMass> {
        &self.masses
    }

    pub fn mass(&self, name: &str) -> Option<&PointMass> {
        self.masses.get(name)
    }

    pub fn springs(&self) -> &BTreeMap<String, Spring> {
        &self.springs
    }

    pub fn spring(&self, name: &str) -> Option<&Spring> {
        self.springs.get(name)
    }

    fn check_mass(&self, name: &str, m: &PointMass) -> Result<()> {
        m.validate().map_err(|e| context("point mass", name, e))?;
        m.target.check_exists(self).map_err(|e| context("point mass", name, e))
    }

    pub fn add_mass(&mut self, name: &str, mass: PointMass) -> Result<()> {
        self.check_mass(name, &mass)?;
        insert_new(&mut self.masses, "point mass", name, mass)
    }

    pub fn replace_mass(&mut self, name: &str, mass: PointMass) -> Result<()> {
        self.check_mass(name, &mass)?;
        replace_existing(&mut self.masses, "point mass", name, mass)
    }

    pub fn remove_mass(&mut self, name: &str) -> Result<()> {
        remove_existing(&mut self.masses, "point mass", name).map(|_| ())
    }

    fn check_spring(&self, name: &str, s: &Spring) -> Result<()> {
        let what = s.kind();
        s.validate().map_err(|e| context(what, name, e))?;
        for t in s.targets() {
            t.check_exists(self).map_err(|e| context(what, name, e))?;
        }
        Ok(())
    }

    pub fn add_spring(&mut self, name: &str, spring: Spring) -> Result<()> {
        self.check_spring(name, &spring)?;
        insert_new(&mut self.springs, "spring", name, spring)
    }

    pub fn replace_spring(&mut self, name: &str, spring: Spring) -> Result<()> {
        self.check_spring(name, &spring)?;
        replace_existing(&mut self.springs, "spring", name, spring)
    }

    pub fn remove_spring(&mut self, name: &str) -> Result<()> {
        remove_existing(&mut self.springs, "spring", name).map(|_| ())
    }

    /// The definitions among point masses, springs and dashpots that act on a target.
    pub(crate) fn discrete_users(&self, uses: impl Fn(&Target) -> bool) -> Vec<String> {
        let mut users: Vec<String> = self.masses.iter().filter(|(_, m)| uses(&m.target)).map(|(n, _)| format!("point mass `{n}`")).collect();
        users.extend(self.springs.iter().filter(|(_, s)| s.targets().any(&uses)).map(|(n, s)| format!("{} `{n}`", s.kind())));
        users
    }

    /// The coordinates of a mesh node or a reference point's node label.
    pub fn label_xyz(&self, label: u32) -> Option<[f64; 3]> {
        self.nodes.get(&label).map(|n| n.xyz).or_else(|| self.reference_points.values().find(|p| p.label == label).map(|p| p.xyz))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Set;

    #[test]
    fn definitions_validate_and_round_trip() {
        let mut m = Model::default();
        m.add_set("tip", Set { nodes: vec![1], ..Set::default() }).unwrap();
        let json = serde_json::json!({"target": {"set": "tip"}, "mass": [1.0, 2.0, 0.0], "inertia": [1.0, 2.0, 3.0, 0.5, 0.0, 0.0]});
        let pm: PointMass = serde_json::from_value(json.clone()).unwrap();
        assert_eq!(pm.mass, Some(Mass::Directional([1.0, 2.0, 0.0])));
        assert_eq!(serde_json::to_value(&pm).unwrap(), json);
        m.add_mass("pm", pm).unwrap();
        let iso: PointMass = serde_json::from_value(serde_json::json!({"target": {"set": "tip"}, "mass": 3})).unwrap();
        assert_eq!(iso.mass.map(Mass::components), Some([3.0; 3]));
        for bad in [
            serde_json::json!({"target": {"set": "tip"}}),
            serde_json::json!({"target": {"set": "tip"}, "mass": -1}),
            serde_json::json!({"target": {"set": "tip"}, "mass": 0}),
            serde_json::json!({"target": {"set": "tip"}, "inertia": [1, 1, 1, 2, 0, 0]}),
            serde_json::json!({"target": {"set": "tip"}, "inertia": [-1, 1, 1, 0, 0, 0]}),
            serde_json::json!({"target": {"set": "tip"}, "inertia": [0, 0, 0, 0, 0, 0]}),
        ] {
            let pm: PointMass = serde_json::from_value(bad.clone()).unwrap();
            assert!(m.add_mass("x", pm).is_err(), "{bad}");
        }
        assert!(serde_json::from_value::<PointMass>(serde_json::json!({"target": {"set": "tip"}, "mass": 1, "masss": 2})).is_err());
        let missing: PointMass = serde_json::from_value(serde_json::json!({"target": {"set": "nope"}, "mass": 1})).unwrap();
        assert!(matches!(m.add_mass("x", missing), Err(Error::NotFound(_))));

        let s: Spring =
            serde_json::from_value(serde_json::json!({"type": "spring", "first": {"set": "tip"}, "direction": {"dof": 3}, "stiffness": 100}))
                .unwrap();
        m.add_spring("k", s).unwrap();
        let d: Spring = serde_json::from_value(
            serde_json::json!({"type": "dashpot", "first": {"set": "tip"}, "direction": {"vector": [1, 1, 0]}, "coefficient": 2}),
        )
        .unwrap();
        assert!(d.is_dashpot());
        m.add_spring("c", d).unwrap();
        for bad in [
            serde_json::json!({"type": "spring", "first": {"set": "tip"}, "direction": "line", "stiffness": 1}),
            serde_json::json!({"type": "spring", "first": {"set": "tip"}, "direction": {"dof": 7}, "stiffness": 1}),
            serde_json::json!({"type": "spring", "first": {"set": "tip"}, "direction": {"vector": [0, 0, 0]}, "stiffness": 1}),
            serde_json::json!({"type": "spring", "first": {"set": "tip"}, "direction": {"dof": 1}, "stiffness": 0}),
            serde_json::json!({"type": "spring", "first": {"set": "tip"}, "direction": {"dof": 1}, "stiffness": f64::NAN}),
            serde_json::json!({"type": "dashpot", "first": {"set": "tip"}, "second": {"set": "tip"}, "direction": {"dof": 1}, "coefficient": 1}),
        ] {
            let Ok(s) = serde_json::from_value::<Spring>(bad.clone()) else { continue };
            assert!(m.add_spring("x", s).is_err(), "{bad}");
        }
        assert!(
            serde_json::from_value::<Spring>(serde_json::json!({"type": "dashpot", "first": {"set": "tip"}, "direction": "line", "stiffness": 1}))
                .is_err()
        );
        // The set is in use.
        assert!(matches!(m.remove_set("tip"), Err(Error::InUse(_))));
        m.remove_mass("pm").unwrap();
        m.remove_spring("k").unwrap();
        m.remove_spring("c").unwrap();
        m.remove_set("tip").unwrap();
    }

    #[test]
    fn directions_give_unit_coefficients() {
        assert_eq!(ResolvedDirection::Dof(5).coefficients(), [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]);
        assert_eq!(ResolvedDirection::Line([0.6, 0.8, 0.0]).coefficients(), [0.6, 0.8, 0.0, 0.0, 0.0, 0.0]);
        assert!(ResolvedDirection::Dof(4).is_rotational() && !ResolvedDirection::Dof(3).is_rotational());
    }

    #[test]
    fn resolution_expands_and_checks_targets() {
        use crate::{Procedure, ReferencePoint, Step};
        let mut m = Model::default();
        m.set_node(1, [0.0; 3]).unwrap();
        m.set_node(2, [3.0, 4.0, 0.0]).unwrap();
        m.add_set("a", Set { nodes: vec![1], ..Set::default() }).unwrap();
        m.add_set("b", Set { nodes: vec![2], ..Set::default() }).unwrap();
        m.add_set("ab", Set { nodes: vec![1, 2], ..Set::default() }).unwrap();
        m.add_reference_point("RP", ReferencePoint { xyz: [0.0, 0.0, 1.0], label: 0 }).unwrap();
        m.add_step(Step::new("Modes", Procedure::Frequency { modes: 1, min_hz: None, max_hz: None }), None).unwrap();
        let spring = |first: &str, second: Option<&str>, direction| Spring::Spring {
            first: Target::Set(first.into()),
            second: second.map(|s| Target::Set(s.into())),
            direction,
            stiffness: 10.0,
        };
        m.add_spring("line", spring("a", Some("b"), Direction::Line)).unwrap();
        m.add_spring("ground", spring("ab", None, Direction::Vector([0.0, 0.0, 2.0]))).unwrap();
        m.add_mass("m", PointMass { target: Target::Set("ab".into()), mass: Some(Mass::Isotropic(2.0)), inertia: None }).unwrap();
        m.add_mass("j", PointMass { target: Target::ReferencePoint("RP".into()), mass: None, inertia: Some([1.0, 1.0, 1.0, 0.0, 0.0, 0.0]) })
            .unwrap();
        // No mesh elements: the discrete elements are the model.
        let r = m.resolve().unwrap();
        assert_eq!(r.springs.len(), 3);
        assert_eq!(r.springs[0].direction, ResolvedDirection::Vector([0.0, 0.0, 1.0]));
        assert_eq!((r.springs[2].first, r.springs[2].second), (1, Some(2)));
        assert_eq!(r.springs[2].direction, ResolvedDirection::Line([0.6, 0.8, 0.0]));
        assert_eq!(r.masses.len(), 3);
        assert_eq!(r.masses.iter().map(|x| x.node).collect::<Vec<_>>(), vec![99_999_999, 1, 2]);

        // Errors: two-point ends with several nodes, rotations at solid nodes, inertia at mesh nodes.
        m.add_spring("wide", spring("ab", Some("b"), Direction::Dof(1))).unwrap();
        m.add_spring("twist", spring("a", None, Direction::Dof(6))).unwrap();
        m.add_mass("bad", PointMass { target: Target::Set("a".into()), mass: None, inertia: Some([1.0, 1.0, 1.0, 0.0, 0.0, 0.0]) }).unwrap();
        let errors = m.resolve().unwrap_err().to_string();
        for want in ["spring `wide`: a two-point spring connects one node to one node", "spring `twist`: acts on a rotation", "point mass `bad`"] {
            assert!(errors.contains(want), "{errors}");
        }
        // A reference point used by a discrete element cannot be deleted.
        assert!(matches!(m.remove_reference_point("RP"), Err(Error::InUse(_))));
        // Torsion at a reference point is fine.
        m.remove_spring("twist").unwrap();
        m.remove_spring("wide").unwrap();
        m.remove_mass("bad").unwrap();
        m.add_spring(
            "twist",
            Spring::Spring { first: Target::ReferencePoint("RP".into()), second: None, direction: Direction::Dof(6), stiffness: 1.0 },
        )
        .unwrap();
        m.resolve().unwrap();
        // The definitions survive a file round trip.
        let back: Model = serde_json::from_str(&serde_json::to_string(&m).unwrap()).unwrap();
        assert_eq!(back, m);
    }
}
