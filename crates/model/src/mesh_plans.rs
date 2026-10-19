//! Mesh plans: how a body of an opened design is to be hex meshed, derived from the features
//! that built it (`mesh.plan`), kept so it can be inspected, edited and applied again
//! (`mesh.plan.apply`) after the design changes.
//!
//! A plan names the body, the element size and the sweep direction, lists the features it maps
//! and the template each one takes (a sweep of the profile, an O-grid round a hole, layers
//! through a shell, an imprint of a joined tool, a copy of a pattern instance), the partition cuts
//! its templates need (each with the feature it is for and why), and the features it does not
//! map yet with the reason. Plans are definitions of the model: journalled with the command that
//! made them and saved in the `.fieldcraft` file.

use serde::{Deserialize, Serialize};

use crate::{Error, Result, check_coordinates, check_name};

/// Most entries in each list of a plan.
pub const MAX_PLAN_ITEMS: usize = 10_000;

/// A direction or axis is a unit vector when its length is within this of 1.
const UNIT_LENGTH: f64 = 1e-6; // unit-factors: the length of a direction vector
/// A normal shorter than this is zero.
const ZERO_NORMAL: f64 = 1e-9; // unit-factors: the length of a direction vector

/// What a plan meshes with today.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PlanMesher {
    /// The whole body swept in one piece (`mesh.hex.generate`, swept methods).
    #[default]
    Sweep,
    /// The body cut by the plan's partitions, its cells swept together (`mesh.cells`).
    Cells,
    /// One block topology laid over the whole body (a grid with O-grids round its circles and a
    /// blend layer along its fillets, chamfers and shell walls) projected onto its faces
    /// (projection blocking, `mesh.hex.project`): no partitions, one conforming mesh.
    Project,
}

/// A feature the plan maps, and the template it takes.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PlanFeature {
    /// The feature's id (`F<id>` in persistent names).
    pub id: u64,
    pub name: String,
    /// The feature's type (`extrude`, `box`, `hole`, …).
    #[serde(rename = "type")]
    pub kind: String,
    /// `sweep`, `block`, `butterfly`, `o_grid`, `layers`, `imprint`, `partition` or `copy`.
    pub template: String,
    /// What the template does with this feature, in a sentence.
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub note: String,
    /// The body's faces the feature made (persistent names): the surfaces its blocks' boundary
    /// nodes project onto (the projection blocking mesher, #399).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub faces: Vec<String>,
}

/// The kind of ring a circle of the body takes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RingKind {
    /// Material inside the circle (a cylinder, a boss): a square core and four blocks round it.
    Butterfly,
    /// Material outside the circle (a hole, a circular cut): four blocks between the circle and
    /// a square frame round it.
    OGrid,
}

/// A circle of the body, square to the sweep, with the ring of blocks that goes round it: the
/// block topology the projection mesher (#399) builds, its corners at 45° on the circle and its
/// radial edges projected onto `faces`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PlanRing {
    /// The feature the circle comes from.
    pub feature: u64,
    pub kind: RingKind,
    /// The circle's centre (a point on its axis), mm.
    pub center: [f64; 3],
    /// Its axis (unit, along the sweep).
    pub axis: [f64; 3],
    pub radius: f64,
    /// Half the side of the square frame (an O-grid's outer boundary, a butterfly's core is
    /// `radius / 2`), mm.
    pub frame: f64,
    /// The cylindrical faces the ring's inner boundary projects onto (persistent names).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub faces: Vec<String>,
}

/// What a cut is for.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CutRole {
    /// Square to the sweep, where a joined or cut tool's section starts or ends.
    Level,
    /// Along the sweep: the frame of an O-grid or butterfly, or of a cut opening.
    #[default]
    Frame,
}

/// A partition cut the plan needs: a plane through `point` with normal `normal`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PlanCut {
    pub point: [f64; 3],
    pub normal: [f64; 3],
    #[serde(default)]
    pub role: CutRole,
    /// The feature whose template needs it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub feature: Option<u64>,
    pub reason: String,
}

/// A feature the plan does not map, and why.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PlanUnmapped {
    pub id: u64,
    pub name: String,
    #[serde(rename = "type")]
    pub kind: String,
    pub reason: String,
}

/// A mesh plan of one body.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MeshPlan {
    /// The body of the design as opened (not a partition cell).
    pub body: String,
    /// Element size, mm.
    pub size: f64,
    /// The chordal tolerance the body and its cells are tessellated with, mm (small enough for
    /// its smallest curve); `None`: the meshers' default.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tolerance: Option<f64>,
    /// Elements through the thinnest wall the features give (a shell's thickness, a plate's).
    pub layers: u32,
    #[serde(default)]
    pub mesher: PlanMesher,
    /// The sweep direction (unit), when the mapped features are prismatic along one.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub direction: Option<[f64; 3]>,
    /// Every direction all mapped features are prismatic along (`direction` first).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub directions: Vec<[f64; 3]>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub features: Vec<PlanFeature>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub cuts: Vec<PlanCut>,
    /// The rings round the body's circles (O-grids at holes, butterflies in bosses).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub rings: Vec<PlanRing>,
    /// The other bodies that are copies of this one (a pattern's or a mirror's), meshed with
    /// plans of their own when this plan is applied.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub copies: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub unmapped: Vec<PlanUnmapped>,
    /// What the planner decided and why (the size it took, a cut it left out).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub notes: Vec<String>,
}

fn unit_length(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

impl MeshPlan {
    pub fn check(&self) -> Result<()> {
        if self.body.trim().is_empty() || self.body.len() > 1024 {
            return Err(Error::Invalid("a mesh plan names its body".into()));
        }
        if !(self.size.is_finite() && self.size > 0.0) {
            return Err(Error::Invalid(format!("the element size must be positive, got {}", self.size)));
        }
        if self.tolerance.is_some_and(|t| !(t.is_finite() && t > 0.0)) {
            return Err(Error::Invalid("the tolerance must be positive".into()));
        }
        if !(1..=100).contains(&self.layers) {
            return Err(Error::Invalid(format!("the layers through a wall must be 1 to 100, got {}", self.layers)));
        }
        if self.features.len() > MAX_PLAN_ITEMS
            || self.cuts.len() > MAX_PLAN_ITEMS
            || self.unmapped.len() > MAX_PLAN_ITEMS
            || self.notes.len() > MAX_PLAN_ITEMS
            || self.rings.len() > MAX_PLAN_ITEMS
            || self.copies.len() > MAX_PLAN_ITEMS
            || self.features.iter().any(|f| f.faces.len() > MAX_PLAN_ITEMS)
            || self.directions.len() > 3
        {
            return Err(Error::Limit(format!("a mesh plan lists at most {MAX_PLAN_ITEMS} features, cuts and unmapped features, and 3 directions")));
        }
        for d in self.direction.iter().chain(&self.directions) {
            check_coordinates(d)?;
            if (unit_length(*d) - 1.0).abs() > UNIT_LENGTH {
                return Err(Error::Invalid(format!("a sweep direction is a unit vector, got {d:?}")));
            }
        }
        for c in &self.cuts {
            check_coordinates(&c.point)?;
            check_coordinates(&c.normal)?;
            if !(unit_length(c.normal) > ZERO_NORMAL) {
                return Err(Error::Invalid("a cut's normal must not be zero".into()));
            }
        }
        for r in &self.rings {
            check_coordinates(&r.center)?;
            check_coordinates(&r.axis)?;
            if (unit_length(r.axis) - 1.0).abs() > UNIT_LENGTH {
                return Err(Error::Invalid(format!("a ring's axis is a unit vector, got {:?}", r.axis)));
            }
            if !(r.radius.is_finite() && r.radius > 0.0 && r.frame.is_finite() && r.frame > 0.0) {
                return Err(Error::Invalid("a ring's radius and frame must be positive".into()));
            }
            if r.faces.len() > MAX_PLAN_ITEMS {
                return Err(Error::Limit(format!("a ring lists at most {MAX_PLAN_ITEMS} faces")));
            }
        }
        if self.copies.iter().any(|c| c.trim().is_empty() || c == &self.body) {
            return Err(Error::Invalid("a plan's copies are other bodies, by name".into()));
        }
        Ok(())
    }

    /// The name a plan of `body` gets by default: `Plan-<body>` with every character a
    /// definition name does not take as `_`.
    pub fn default_name(body: &str) -> String {
        let clean: String = body.chars().map(|c| if c.is_ascii_alphanumeric() || c == '_' || c == '-' { c } else { '_' }).collect();
        let name: String = format!("Plan-{clean}").chars().take(crate::MAX_NAME_LEN).collect();
        if check_name("mesh plan", &name).is_ok() { name } else { "Plan".into() }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plan() -> MeshPlan {
        MeshPlan {
            body: "Base".into(),
            size: 2.0,
            tolerance: Some(0.01),
            layers: 2,
            mesher: PlanMesher::Cells,
            direction: Some([0.0, 0.0, 1.0]),
            directions: vec![[0.0, 0.0, 1.0]],
            features: vec![PlanFeature {
                id: 1,
                name: "Box1".into(),
                kind: "box".into(),
                template: "block".into(),
                note: String::new(),
                faces: vec!["F1:+z".into()],
            }],
            cuts: vec![PlanCut { point: [1.0, 0.0, 0.0], normal: [1.0, 0.0, 0.0], role: CutRole::Frame, feature: Some(1), reason: "frame".into() }],
            rings: vec![PlanRing {
                feature: 1,
                kind: RingKind::OGrid,
                center: [0.0; 3],
                axis: [0.0, 0.0, 1.0],
                radius: 1.0,
                frame: 2.0,
                faces: Vec::new(),
            }],
            copies: vec!["Base.2".into()],
            unmapped: Vec::new(),
            notes: Vec::new(),
        }
    }

    #[test]
    fn a_plan_round_trips_and_checks() {
        let p = plan();
        p.check().unwrap();
        let back: MeshPlan = serde_json::from_value(serde_json::to_value(&p).unwrap()).unwrap();
        assert_eq!(back, p);
    }

    #[test]
    fn hostile_plans_are_refused() {
        let bad = [
            MeshPlan { size: f64::NAN, ..plan() },
            MeshPlan { size: 0.0, ..plan() },
            MeshPlan { layers: 0, ..plan() },
            MeshPlan { tolerance: Some(-1.0), ..plan() },
            MeshPlan { body: String::new(), ..plan() },
            MeshPlan { direction: Some([0.0, 0.0, 2.0]), ..plan() },
            MeshPlan { rings: vec![PlanRing { radius: -1.0, ..plan().rings[0].clone() }], ..plan() },
            MeshPlan { rings: vec![PlanRing { axis: [0.0, 0.0, 0.0], ..plan().rings[0].clone() }], ..plan() },
            MeshPlan { copies: vec!["Base".into()], ..plan() },
            MeshPlan {
                cuts: vec![PlanCut { point: [0.0; 3], normal: [0.0; 3], role: CutRole::Level, feature: None, reason: String::new() }],
                ..plan()
            },
            MeshPlan {
                cuts: vec![PlanCut {
                    point: [f64::INFINITY, 0.0, 0.0],
                    normal: [1.0, 0.0, 0.0],
                    role: CutRole::Frame,
                    feature: None,
                    reason: String::new(),
                }],
                ..plan()
            },
        ];
        for p in bad {
            assert!(p.check().is_err(), "{p:?}");
        }
    }

    #[test]
    fn default_names_are_definition_names() {
        assert_eq!(MeshPlan::default_name("Base"), "Plan-Base");
        assert_eq!(MeshPlan::default_name("Pins.3"), "Plan-Pins_3");
        assert!(check_name("mesh plan", &MeshPlan::default_name(&"x".repeat(500))).is_ok());
    }
}
