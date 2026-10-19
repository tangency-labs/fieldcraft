//! The feature graph of an opened design: for each body, the features of the SolidCraft timeline
//! that built it, with their parameters evaluated and the persistent face and edge names each
//! one left on the body.
//!
//! A design opened from `.solidcraft` is rebuilt from its timeline, so how every body was made is
//! known, not guessed: an extrude of a sketch profile along a direction, a box, a shell of a
//! thickness, a hole of a diameter. The provenance-driven mesh plan (`mesh.plan`) reads it.
//!
//! Which features built a body comes from two places: the persistent face names (each names the
//! features that made the face, `F8:side:l1`, `F3:inner:F1:+y`, `F82:-z@F84.2` for a pattern
//! copy), and the features that name the body (targets, tools, the body a fillet or shell
//! works on). A feature whose faces were all removed later and that names no body (a cut that a
//! later cut swallowed) is left out.

use std::collections::{BTreeMap, BTreeSet};

use serde::Serialize;
use solidcraft_doc::expr::Kind;
use solidcraft_doc::{Document, FeatureKind, ModelState, Operation, SolvedSketch, plain};
use solidcraft_geom::{Seg2, Vec3};
use solidcraft_kernel::guard;

use crate::{Design, Result, Shape};

/// Most features listed per body (a design's timeline is capped far above any real part).
pub const MAX_FEATURES_PER_BODY: usize = 10_000;

/// The timeline a design was rebuilt from, kept with the design for [`Design::features`].
pub(crate) struct History {
    pub(crate) doc: Document,
    pub(crate) state: std::sync::Arc<ModelState>,
    /// Error per feature id, from the rebuild.
    pub(crate) errors: BTreeMap<u64, String>,
}

impl std::fmt::Debug for History {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("History").field("features", &self.doc.features.len()).field("errors", &self.errors).finish()
    }
}

/// The features that built one body, in timeline order.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct BodyFeatures {
    pub body: String,
    pub features: Vec<FeatureNode>,
}

/// One feature of the timeline as it built a body.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct FeatureNode {
    /// The feature's id: `F<id>` in persistent names.
    pub id: u64,
    pub name: String,
    /// The definition's type (`extrude`, `box`, `cylinder`, `hole`, `shell`, `fillet`, …).
    #[serde(rename = "type")]
    pub kind: String,
    /// `new_body`, `join`, `cut` or `intersect`, for features that have one.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub operation: Option<String>,
    /// The numeric inputs evaluated: lengths in mm, angles in radians, counts as numbers, keyed
    /// by their snake-case label (`distance`, `radius`, `thickness`, …).
    pub values: BTreeMap<String, f64>,
    /// The keys of `values` that are lengths.
    #[serde(skip)]
    pub lengths: BTreeSet<String>,
    /// Where a primitive sits (a box's corner, a cylinder's base, a sphere's centre), mm.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub placement: Option<[f64; 3]>,
    /// The direction the feature builds along: an extrude's (the sketch normal, flipped for a
    /// negative extrude), a cylinder's axis, a hole's drilling direction.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub direction: Option<[f64; 3]>,
    /// An extrude in both directions (`symmetric` or a second distance).
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    pub two_sided: bool,
    /// The sketch the feature takes its profiles from, with the profiles it uses.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sketch: Option<SketchGeometry>,
    /// Pattern and mirror copies: each copy's 4×4 transform (row-major), the identity excluded.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub transforms: Vec<[[f64; 4]; 4]>,
    /// The bodies a join, cut or intersect changes (empty: all), or a combine's target.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub targets: Vec<String>,
    /// The tool bodies of a combine.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub tools: Vec<String>,
    /// The names the feature gives the bodies it makes (a combine's tool is found by them).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub creates: Vec<String>,
    /// The definition as stored in the design (expressions unevaluated).
    pub definition: serde_json::Value,
    /// The body's faces this feature made or shaped (their persistent names name it).
    pub faces: Vec<String>,
    /// The body's edges on those faces.
    pub edges: Vec<String>,
    /// Why the feature failed to rebuild, when it did.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

/// A sketch plane and the profiles a feature uses from it (sketch coordinates, mm).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SketchGeometry {
    pub id: u64,
    pub name: String,
    pub origin: [f64; 3],
    pub x_dir: [f64; 3],
    pub y_dir: [f64; 3],
    pub normal: [f64; 3],
    pub profiles: Vec<ProfileGeometry>,
}

/// A closed profile: its outer loop and holes.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ProfileGeometry {
    /// Index in the solved sketch's profile list.
    pub index: usize,
    pub area: f64,
    pub outer: Vec<Curve>,
    pub holes: Vec<Vec<Curve>>,
    /// The sketch curve ids of the outer loop and of each hole.
    pub outer_curves: Vec<String>,
    pub hole_curves: Vec<Vec<String>>,
}

/// A curve of a profile loop (sketch coordinates). A loop of arcs about one centre that closes is
/// one [`Curve::Circle`].
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Curve {
    Line {
        a: [f64; 2],
        b: [f64; 2],
    },
    Arc {
        center: [f64; 2],
        radius: f64,
        start: f64,
        sweep: f64,
    },
    Circle {
        center: [f64; 2],
        radius: f64,
    },
    /// A spline, conic or other curve the planner does not decompose.
    Other {
        kind: String,
    },
}

impl BodyFeatures {
    /// Every length (by `length`) and area (by `area`) of the graph converted: the graph is in
    /// millimetres, as the design.
    pub fn convert(&mut self, length: &dyn Fn(f64) -> f64, area: &dyn Fn(f64) -> f64) {
        for f in &mut self.features {
            for (k, v) in f.values.iter_mut() {
                if f.lengths.contains(k) {
                    *v = length(*v);
                }
            }
            if let Some(p) = &mut f.placement {
                *p = p.map(length);
            }
            for t in &mut f.transforms {
                for row in t.iter_mut().take(3) {
                    row[3] = length(row[3]);
                }
            }
            if let Some(sk) = &mut f.sketch {
                sk.origin = sk.origin.map(length);
                for p in &mut sk.profiles {
                    p.area = area(p.area);
                    for c in p.outer.iter_mut().chain(p.holes.iter_mut().flatten()) {
                        match c {
                            Curve::Line { a, b } => {
                                *a = a.map(length);
                                *b = b.map(length);
                            }
                            Curve::Arc { center, radius, start: _, sweep: _ } | Curve::Circle { center, radius } => {
                                *center = center.map(length);
                                *radius = length(*radius);
                            }
                            Curve::Other { kind: _ } => {}
                        }
                    }
                }
            }
        }
    }
}

fn xyz(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

/// `v` scaled to unit length, components within 1e-12 of zero made zero (a drilling direction
/// found from the geometry carries round-off).
fn unit(v: [f64; 3]) -> Option<[f64; 3]> {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    (n.is_finite() && n > 1e-12).then(|| v.map(|c| if (c / n).abs() < 1e-12 { 0.0 } else { c / n }))
}

/// The ids of the features a persistent face name mentions (`F3:inner:F1:+y` → 3, 1).
pub fn feature_ids(name: &str) -> impl Iterator<Item = u64> + '_ {
    let bytes = name.as_bytes();
    (0..bytes.len()).filter_map(move |i| {
        let starts = bytes.get(i) == Some(&b'F') && (i == 0 || !bytes.get(i - 1).is_some_and(|b| b.is_ascii_alphanumeric()));
        if !starts {
            return None;
        }
        let digits: &str = name.get(i + 1..)?.split(|c: char| !c.is_ascii_digit()).next()?;
        if digits.is_empty() { None } else { digits.parse().ok() }
    })
}

/// `Taper angle` → `taper_angle`, `Distance (0 to 1)` → `distance_0_to_1`.
fn snake(label: &str) -> String {
    let mut out = String::with_capacity(label.len());
    for c in label.chars() {
        if c.is_ascii_alphanumeric() {
            out.push(c.to_ascii_lowercase());
        } else if !out.ends_with('_') && !out.is_empty() {
            out.push('_');
        }
    }
    out.trim_end_matches('_').to_string()
}

fn operation_name(op: Operation) -> &'static str {
    match op {
        Operation::NewBody => "new_body",
        Operation::Join => "join",
        Operation::Cut => "cut",
        Operation::Intersect => "intersect",
    }
}

fn curves(segs: &[Seg2]) -> Vec<Curve> {
    // A loop of arcs about one centre and radius that turns a full circle is a circle.
    let arcs: Vec<(solidcraft_geom::Vec2, f64, f64)> =
        segs.iter().filter_map(|s| if let Seg2::Arc { center, radius, sweep, .. } = s { Some((*center, *radius, *sweep)) } else { None }).collect();
    if !arcs.is_empty() && arcs.len() == segs.len() {
        let (c0, r0, _) = arcs[0];
        let same = arcs.iter().all(|(c, r, _)| (c.x - c0.x).abs() < 1e-9 && (c.y - c0.y).abs() < 1e-9 && (r - r0).abs() < 1e-9);
        let turn: f64 = arcs.iter().map(|a| a.2).sum();
        if same && (turn.abs() - std::f64::consts::TAU).abs() < 1e-6 {
            return vec![Curve::Circle { center: [c0.x, c0.y], radius: r0 }];
        }
    }
    segs.iter()
        .map(|s| match s {
            Seg2::Line { a, b } => Curve::Line { a: [a.x, a.y], b: [b.x, b.y] },
            Seg2::Arc { center, radius, start, sweep } => Curve::Arc { center: [center.x, center.y], radius: *radius, start: *start, sweep: *sweep },
            other => Curve::Other { kind: format!("{other:?}").split(['{', ' ', '(']).next().unwrap_or("curve").to_ascii_lowercase() },
        })
        .collect()
}

fn sketch_geometry(ss: &SolvedSketch, sel: Option<&solidcraft_doc::ProfileSel>) -> SketchGeometry {
    let indices: Vec<usize> = match sel {
        Some(sel) => solidcraft_doc::profile_selection_indices(ss, sel).unwrap_or_else(|_| (0..ss.profiles.len()).collect()),
        None => (0..ss.profiles.len()).collect(),
    };
    let profiles = indices
        .into_iter()
        .filter_map(|i| {
            let p = ss.profiles.get(i)?;
            Some(ProfileGeometry {
                index: i,
                area: p.area,
                outer: curves(&p.region.outer.segs),
                holes: p.region.holes.iter().map(|h| curves(&h.segs)).collect(),
                outer_curves: p.outer_curves.clone(),
                hole_curves: p.hole_curves.clone(),
            })
        })
        .collect();
    SketchGeometry {
        id: ss.feature,
        name: ss.name.clone(),
        origin: xyz(ss.plane.origin),
        x_dir: xyz(ss.plane.x),
        y_dir: xyz(ss.plane.y),
        normal: xyz(ss.plane.normal()),
        profiles,
    }
}

impl History {
    /// The feature `id` of the timeline, described against `body` (its faces and edges).
    fn node(&self, f: &solidcraft_doc::Feature, faces: Vec<String>, edges: Vec<String>) -> FeatureNode {
        let (vals, _) = self.doc.param_values();
        let eval = |e: &str, k: Kind| Document::eval_in(&vals, e, k).ok().filter(|v| v.is_finite());
        let mut values = BTreeMap::new();
        let mut lengths = BTreeSet::new();
        for (label, e, k) in f.kind.inputs() {
            if let Some(v) = eval(&e, k) {
                let mut key = snake(label);
                let mut n = 2;
                while values.contains_key(&key) {
                    key = format!("{}_{n}", snake(label));
                    n += 1;
                }
                if k == Kind::Length {
                    lengths.insert(key.clone());
                }
                values.insert(key, v);
            }
        }
        let point = |p: &[solidcraft_doc::Expr; 3]| -> Option<[f64; 3]> {
            Some([eval(&p[0], Kind::Length)?, eval(&p[1], Kind::Length)?, eval(&p[2], Kind::Length)?])
        };
        let sketch_of = |id: solidcraft_doc::SketchId| self.state.sketch(*id);
        let mut node = FeatureNode {
            id: f.id,
            name: f.name.clone(),
            kind: String::new(),
            operation: None,
            values,
            lengths,
            placement: None,
            direction: None,
            two_sided: false,
            sketch: None,
            transforms: Vec::new(),
            targets: Vec::new(),
            tools: Vec::new(),
            creates: f.body_names.clone(),
            definition: serde_json::to_value(&f.kind).unwrap_or_default(),
            faces,
            edges,
            error: self.errors.get(&f.id).cloned(),
        };
        node.kind = node.definition.get("type").and_then(serde_json::Value::as_str).unwrap_or("feature").to_string();
        match &f.kind {
            FeatureKind::Extrude { sketch, profiles, extent, operation, targets, .. } => {
                node.operation = Some(operation_name(*operation).into());
                node.targets = plain(targets);
                if let Some(ss) = sketch_of(*sketch) {
                    let n = xyz(ss.plane.normal());
                    let negative = matches!(extent.direction, solidcraft_doc::Direction::Negative);
                    node.direction = unit(if negative { n.map(|c| -c) } else { n });
                    node.two_sided = matches!(extent.direction, solidcraft_doc::Direction::Symmetric) || extent.distance2.is_some();
                    node.sketch = Some(sketch_geometry(ss, Some(profiles)));
                }
            }
            FeatureKind::Revolve { sketch, profiles, operation, targets, .. } | FeatureKind::Sweep { sketch, profiles, operation, targets, .. } => {
                node.operation = Some(operation_name(*operation).into());
                node.targets = plain(targets);
                node.sketch = sketch_of(*sketch).map(|ss| sketch_geometry(ss, Some(profiles)));
            }
            FeatureKind::Loft { operation, targets, .. } => {
                node.operation = Some(operation_name(*operation).into());
                node.targets = plain(targets);
            }
            FeatureKind::Box { corner, operation, .. } => {
                node.operation = Some(operation_name(*operation).into());
                node.placement = point(corner);
                node.direction = Some([0.0, 0.0, 1.0]);
            }
            FeatureKind::Cylinder { base, axis, operation, .. } => {
                node.operation = Some(operation_name(*operation).into());
                node.placement = point(base);
                node.direction = unit(xyz(**axis));
            }
            FeatureKind::Sphere { center, operation, .. } | FeatureKind::Torus { center, operation, .. } => {
                node.operation = Some(operation_name(*operation).into());
                node.placement = point(center);
            }
            FeatureKind::Hole { position, direction, .. } => {
                node.operation = Some("cut".into());
                node.placement = Some(xyz(**position));
                node.direction = unit(xyz(**direction));
            }
            FeatureKind::Combine { target, tools, operation, .. } => {
                node.operation = Some(operation_name(*operation).into());
                node.targets = vec![target.to_string()];
                node.tools = plain(tools);
            }
            FeatureKind::Shell { body, .. } | FeatureKind::Fillet { body, .. } | FeatureKind::Chamfer { body, .. } => {
                node.targets = body.iter().map(ToString::to_string).collect();
            }
            FeatureKind::Split { body, .. } => node.targets = vec![body.to_string()],
            FeatureKind::Pattern { pattern, bodies, .. } => {
                node.targets = plain(bodies);
                if let Ok(m) = solidcraft_doc::pattern_matrices(&self.doc, &self.state, pattern) {
                    // SolidCraft's matrices are column-major (`m[3]` the translation); the graph's are row-major.
                    let rows = |m: [[f64; 4]; 4]| std::array::from_fn(|i| std::array::from_fn(|j| m[j][i]));
                    node.transforms = m.into_iter().filter(|m| *m != solidcraft_doc::IDENTITY).map(rows).collect();
                }
            }
            FeatureKind::Mirror { plane, bodies, .. } => {
                node.targets = plain(bodies);
                if let Ok(p) = self.doc.resolve_plane_in(&vals, plane, 0, Some(&self.state))
                    && let Some(n) = unit(xyz(p.normal()))
                {
                    // x' = x − 2 ((x − o)·n) n
                    let o = xyz(p.origin);
                    let d = 2.0 * (o[0] * n[0] + o[1] * n[1] + o[2] * n[2]);
                    let mut m = solidcraft_doc::IDENTITY;
                    for (i, row) in m.iter_mut().take(3).enumerate() {
                        for (j, c) in row.iter_mut().take(3).enumerate() {
                            *c -= 2.0 * n[i] * n[j];
                        }
                        row[3] = d * n[i];
                    }
                    node.transforms = vec![m];
                }
            }
            _ => {}
        }
        node
    }
}

impl Design {
    /// The features that built each body of the design as opened (`body`: only that one), in
    /// timeline order. Partition cells are described by the body they were cut from. An error
    /// for a design not opened from a timeline (shapes built in code) or a body it lacks.
    pub fn features(&self, body: Option<&str>) -> Result<Vec<BodyFeatures>> {
        let history = self.history.as_deref().ok_or_else(|| crate::Error::Invalid("the design has no feature timeline".into()))?;
        let shapes: Vec<&Shape> = match body {
            Some(name) => {
                let source = self.source_of(name).ok_or_else(|| crate::Error::Invalid(format!("no body `{name}`")))?;
                self.opened.iter().filter(|b| b.name == source).collect()
            }
            None => self.opened.iter().collect(),
        };
        shapes.into_iter().map(|shape| guard("feature graph", || Ok(body_features(history, shape))).map_err(crate::Error::from)?).collect()
    }
}

fn body_features(history: &History, shape: &Shape) -> Result<BodyFeatures> {
    // Faces by the features their names mention (a merged face: the names of its B-rep faces).
    let mut faces: BTreeMap<u64, Vec<String>> = BTreeMap::new();
    for f in shape.faces() {
        let ids: BTreeSet<u64> =
            if f.members.is_empty() { feature_ids(&f.name).collect() } else { f.members.iter().flat_map(|m| feature_ids(m)).collect() };
        for id in ids {
            faces.entry(id).or_default().push(f.name.clone());
        }
    }
    let edges = shape.edges()?;
    // The body as SolidCraft names it before split pieces were numbered (`Bracket.2` was
    // `Bracket (1)`): features name it by its base.
    let base = shape.name.rsplit_once('.').filter(|(_, n)| n.bytes().all(|b| b.is_ascii_digit())).map_or(shape.name.as_str(), |(b, _)| b);
    let names_body = |f: &solidcraft_doc::Feature| {
        let mentioned = |n: &str| n == shape.name || n == base || crate::design::split_base(n) == base;
        // The bodies the feature changes. A combine's tools are not among them: kept tools come
        // out unchanged, and the others are gone.
        let targets: Vec<&solidcraft_doc::BodyName> = match &f.kind {
            FeatureKind::Extrude { targets, .. }
            | FeatureKind::Revolve { targets, .. }
            | FeatureKind::Sweep { targets, .. }
            | FeatureKind::Loft { targets, .. } => targets.iter().collect(),
            FeatureKind::Combine { target, .. } => vec![target],
            FeatureKind::Shell { body, .. } | FeatureKind::Fillet { body, .. } | FeatureKind::Chamfer { body, .. } => body.iter().collect(),
            FeatureKind::Split { body, .. } => vec![body],
            // A body pattern or mirror: the original is the instance the copies repeat.
            FeatureKind::Pattern { bodies, .. } | FeatureKind::Mirror { bodies, .. } => bodies.iter().collect(),
            _ => Vec::new(),
        };
        targets.into_iter().map(|b| &**b).chain(f.body_names.iter().map(String::as_str)).any(mentioned)
    };
    let mut features = Vec::new();
    for f in history.doc.features.iter().filter(|f| !f.suppressed) {
        if matches!(
            f.kind,
            FeatureKind::Sketch { .. }
                | FeatureKind::ConstructionPlane { .. }
                | FeatureKind::ConstructionAxis { .. }
                | FeatureKind::ConstructionPoint { .. }
        ) {
            continue;
        }
        let own = faces.remove(&f.id).unwrap_or_default();
        if own.is_empty() && !names_body(f) {
            continue;
        }
        let tagged: BTreeSet<u32> = own.iter().filter_map(|n| shape.face_index(n)).collect();
        let on: Vec<String> = edges.iter().filter(|e| e.faces.iter().any(|i| tagged.contains(i))).map(|e| e.name.clone()).collect();
        features.push(history.node(f, own, on));
        if features.len() >= MAX_FEATURES_PER_BODY {
            break;
        }
    }
    Ok(BodyFeatures { body: shape.name.clone(), features })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn feature_ids_come_from_every_name_part() {
        let ids = |n: &str| feature_ids(n).collect::<Vec<_>>();
        assert_eq!(ids("F3:inner:F2:blend:F1:+y|F1:-x"), vec![3, 2, 1, 1]);
        assert_eq!(ids("F82:-z@F84.2"), vec![82, 84]);
        assert_eq!(ids("F44:side:l12"), vec![44]);
        assert_eq!(ids("partition:Half:a"), Vec::<u64>::new());
        assert_eq!(ids("OFF7"), Vec::<u64>::new());
        assert_eq!(ids("F"), Vec::<u64>::new());
    }

    #[test]
    fn the_plate_is_an_extrude_of_its_outline_and_a_cut_of_a_circle() {
        use crate::sample::{Plate, plate_design};
        let plate = Plate::default();
        let design = crate::open_design(&plate_design(&plate).unwrap()).unwrap();
        let graph = design.features(None).unwrap();
        assert_eq!(graph.len(), 1);
        let f = &graph[0].features;
        assert_eq!(
            f.iter().map(|n| (n.kind.as_str(), n.operation.as_deref())).collect::<Vec<_>>(),
            [("extrude", Some("new_body")), ("extrude", Some("cut"))]
        );
        let (outline, hole) = (&f[0], &f[1]);
        // Expressions come out evaluated: the thickness parameter, and `thickness + 1 mm`.
        assert_eq!(outline.values.get("distance"), Some(&plate.thickness));
        assert_eq!(hole.values.get("distance"), Some(&(plate.thickness + 1.0)));
        assert_eq!(outline.direction, Some([0.0, 0.0, 1.0]));
        let sketch = outline.sketch.as_ref().unwrap();
        assert_eq!(sketch.normal, [0.0, 0.0, 1.0]);
        assert_eq!(sketch.profiles.len(), 1);
        assert!(sketch.profiles[0].outer.iter().all(|c| matches!(c, Curve::Line { .. })) && sketch.profiles[0].outer.len() == 4);
        let circle = &hole.sketch.as_ref().unwrap().profiles[0].outer;
        let c = plate.depth / 2.0;
        assert_eq!(circle, &[Curve::Circle { center: [c, c], radius: plate.hole_diameter / 2.0 }]);
        // Each feature's faces: the outline's four sides and caps, the hole's two halves.
        assert_eq!(outline.faces.len(), 6, "{:?}", outline.faces);
        assert_eq!(hole.faces.len(), 2, "{:?}", hole.faces);
        assert!(!hole.edges.is_empty());
        assert!(outline.error.is_none() && hole.error.is_none());
        // Bodies the design lacks are an error, not an empty graph.
        assert!(design.features(Some("Nope")).is_err());
    }

    #[test]
    fn labels_become_snake_case_keys() {
        assert_eq!(snake("Taper angle"), "taper_angle");
        assert_eq!(snake("Distance (0 to 1)"), "distance_0_to_1");
        assert_eq!(snake("Radius"), "radius");
    }

    /// Primitives, a shell, a hole, a combine, a body pattern and a mirror: each body lists the
    /// features that built it, with their placements, directions, targets and transforms.
    #[cfg(feature = "scripts")]
    #[test]
    fn a_timeline_of_primitives_and_copies_is_listed_per_body() {
        let script = r#"[
            {"command": "solid.box", "params": {"length": 40, "width": 30, "height": 20, "body_name": "Case"}},
            {"command": "solid.shell", "params": {"faces": [[20, 15, 20]], "thickness": 2, "body": "Case"}},
            {"command": "solid.cylinder", "params": {"diameter": 6, "height": 5, "base": [10, 10, 2], "operation": "join", "targets": ["Case"]}},
            {"command": "solid.hole", "params": {"position": [10, 10, 7], "diameter": 2, "depth": 4}},
            {"command": "solid.box", "params": {"length": 4, "width": 4, "height": 20, "corner": [30, 20, 0], "body_name": "Slot"}},
            {"command": "solid.combine", "params": {"target": "Case", "tools": ["Slot"], "operation": "cut"}},
            {"command": "solid.box", "params": {"length": 1, "width": 1, "height": 3, "corner": [60, 0, 0], "body_name": "Pin"}},
            {"command": "solid.pattern.rectangular", "params": {"bodies": ["Pin"], "dir1": [1, 0, 0], "count1": 3, "spacing1": 2}},
            {"command": "solid.box", "params": {"length": 2, "width": 2, "height": 2, "corner": [0, 50, 0], "body_name": "Foot"}},
            {"command": "solid.mirror", "params": {"bodies": ["Foot"], "plane": {"origin": [10, 0, 0], "normal": [1, 0, 0]}}}
        ]"#;
        let design = crate::open_design(&crate::script_design(script).unwrap()).unwrap();
        let case = design.features(Some("Case")).unwrap();
        let kinds: Vec<&str> = case[0].features.iter().map(|f| f.kind.as_str()).collect();
        assert_eq!(kinds, ["box", "shell", "cylinder", "hole", "box", "combine"], "{case:#?}");
        let f = &case[0].features;
        assert_eq!(f[0].placement, Some([0.0, 0.0, 0.0]));
        assert_eq!((f[0].values["length"], f[0].values["width"], f[0].values["height"]), (40.0, 30.0, 20.0));
        assert_eq!(f[1].values["thickness"], 2.0);
        assert_eq!((f[2].values["radius"], f[2].direction, f[2].operation.as_deref()), (3.0, Some([0.0, 0.0, 1.0]), Some("join")));
        assert_eq!((f[3].values["diameter"], f[3].direction, f[3].placement), (2.0, Some([0.0, 0.0, -1.0]), Some([10.0, 10.0, 7.0])));
        assert_eq!((f[5].targets.as_slice(), f[5].tools.as_slice()), (&["Case".to_string()][..], &["Slot".to_string()][..]));
        // The cut tool's faces are the slot's walls in the case.
        assert!(!f[4].faces.is_empty());
        // The pattern's copies: each a body of its own with the pattern's transforms.
        let pins: Vec<&str> = design.bodies.iter().map(|b| b.name.as_str()).filter(|n| n.starts_with("Pin")).collect();
        assert_eq!(pins.len(), 3, "{pins:?}");
        for pin in pins {
            let g = design.features(Some(pin)).unwrap();
            let kinds: Vec<&str> = g[0].features.iter().map(|f| f.kind.as_str()).collect();
            assert_eq!(kinds, ["box", "pattern"], "{pin}");
            let t = &g[0].features[1].transforms;
            assert_eq!(t.len(), 2);
            assert!(t.iter().any(|m| (m[0][3] - 2.0).abs() < 1e-9) && t.iter().any(|m| (m[0][3] - 4.0).abs() < 1e-9), "{t:?}");
        }
        // The mirror: x' = 20 − x.
        let foot = design.bodies.iter().find(|b| b.name.starts_with("Foot") && b.name != "Foot").map(|b| b.name.clone()).unwrap_or_default();
        let g = design.features(Some(&foot)).unwrap();
        let mirror = g[0].features.iter().find(|f| f.kind == "mirror").unwrap_or_else(|| panic!("{foot}: {g:#?}"));
        assert_eq!(mirror.transforms, vec![[[-1.0, 0.0, 0.0, 20.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]]);
        // The cut tool is gone; the graph of every body comes in the design's order.
        let all = design.features(None).unwrap();
        assert_eq!(all.iter().map(|b| b.body.as_str()).collect::<Vec<_>>(), design.bodies.iter().map(|b| b.name.as_str()).collect::<Vec<_>>());
    }
}
