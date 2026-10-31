//! Provenance-driven mesh plans (#402): the features that built a body (`geometry.features`)
//! become its hex-meshing plan, and the plan is carried out with today's meshers.
//!
//! `mesh.plan` reads the body's feature graph and maps each feature to a template:
//! - an extrude of a sketch: the sweep, along its stored direction, of the profile's cross-section
//!   (circles butterflied, circular holes O-gridded);
//! - a box: a block (or a part of the swept cross-section); a cylinder: a butterfly, an O-grid
//!   when it cuts;
//! - a hole: an O-grid at its stored diameter and axis;
//! - a shell: layers through its stored thickness (the element size takes the layers through
//!   the thinnest wall);
//! - a join or cut (and a combine): the tool's cross-section imprinted on the swept one;
//! - a pattern or mirror: the instance meshed and repeated;
//! - a split body: its partition cuts.
//!
//! The sweep direction is the direction most mapped features are prismatic along (the first new
//! body's stored direction on a tie); a feature that is not prismatic along it is reported as not
//! mapped. Fillets, chamfers, revolves, sweeps and lofts are reported as not mapped (their blend,
//! revolve and sweep templates are step 3 of #402); nothing is faked.
//!
//! The butterflies and O-grids need a frame round each circle: the plan cuts the body by planes
//! square to the sweep on either side of every circle that is not the body's own outline
//! ([`frame_cuts`]), placed so the frame is square (as far from the circle as the body's nearest
//! side, when that is close) and never on or near a face of the body, through another circle, or
//! next to another cut. `mesh.plan.apply` adds those cuts as partitions of the body (only those),
//! sweeps its cells together (`mesh.cells`), or sweeps the whole body when the plan has no cuts
//! (`mesh.hex.generate`), and reports the quality and whether the cells joined node to node.
//!
//! The plan is a definition of the model ([`MeshPlan`]): saved with it, inspectable, editable
//! (apply an edited plan object) and applied again after the design changes. The projection
//! blocking mesher (#399) consumes the same plan: its direction, templates and cuts.

use std::collections::BTreeSet;

use fieldcraft_geometry::{Curve, Design, FaceSurface, FeatureNode, Shape};
use fieldcraft_mesh::project::{ProjectOptions, Topology};
use fieldcraft_model::geometry::Entity;
use fieldcraft_model::units::{Quantity, Typed};
use fieldcraft_model::{CutRole, Datum, MeshPlan, Partition, PlanCut, PlanFeature, PlanMesher, PlanRing, PlanUnmapped, RingKind};
use serde_json::{Value, json};

/// A vector shorter than this has no direction.
const MIN_NORM: f64 = 1e-12; // unit-factors: the length of a direction vector
/// Unit directions are parallel (or square) within this.
const SAME_DIRECTION: f64 = 1e-6; // unit-factors: a component of a unit vector
/// A unit direction's component taken as not zero.
const NON_ZERO_COMPONENT: f64 = 1e-9; // unit-factors: a component of a unit vector
/// Circles agree within this fraction of their radius (at least a millimetre).
const SAME_CIRCLE: f64 = 1e-6; // unit-factors: relative to the radius
/// Concentric circles of radii this close (relative) are one.
const CONCENTRIC: f64 = 1e-9; // unit-factors: relative to the radius
/// A circle's centre is on a cylinder's axis within this, relative to its radius.
const ON_AXIS: f64 = 1e-9; // unit-factors: relative to the radius
/// Sketch lines whose directions cross by less than this are parallel.
const SAME_LINE: f64 = 1e-9; // unit-factors: the sine of an angle
/// The meshers' default chordal tolerance, relative to the diagonal.
const DEFAULT_RELATIVE_TOLERANCE: f64 = 1e-3; // unit-factors: relative to the diagonal
/// Corner planes this close, relative to the diagonal, are one.
const CORNER_EPS: f64 = 1e-6; // unit-factors: relative to the diagonal
/// Scaled Jacobians this close are equal.
const SAME_QUALITY: f64 = 1e-9; // unit-factors: a difference of scaled Jacobians
/// How much worse than the whole cell a mesh round the corners may be.
const QUALITY_SLACK: f64 = 1e-3; // unit-factors: a difference of scaled Jacobians
/// Distances of the geometry this close are one (a cut on a face, a section's end).
const SAME_LENGTH: Typed = Typed::new(1e-6, "mm"); // unit-factors: an absolute amount, in the geometry frame's millimetres
/// Coordinates this close are one (a sketch's, a cell's bounds).
const SAME_COORDINATE: Typed = Typed::new(1e-9, "mm"); // unit-factors: an absolute amount, in the geometry frame's millimetres
/// A sketch line shorter than this has no direction.
const MIN_SEGMENT: Typed = Typed::new(1e-12, "mm"); // unit-factors: a floor against a zero length, in the geometry frame's millimetres
/// Axes whose squared distance apart is below this are one.
const SAME_AXIS: Typed = Typed::new(1e-12, "mm^2"); // unit-factors: an absolute amount, in the geometry frame's millimetres

mod blocks;
mod decompose;
mod pillow;

use crate::requirement::Requirement;
use crate::schema::ProblemKind;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

/// How a message shows a number of the geometry (millimetres): in the display system, with its
/// unit ([`crate::units::shown`]).
pub(crate) type Show<'a> = &'a dyn Fn(Quantity, f64) -> String;

/// An absolute length of the geometry frame ([`crate::units::geometry_amount`]).
fn length(t: Typed) -> f64 {
    crate::units::geometry_amount(Quantity::Length, t)
}

/// [`SAME_AXIS`] in the geometry frame.
fn axis_gap() -> f64 {
    crate::units::geometry_amount(Quantity::Area, SAME_AXIS)
}

/// [`Show`] for `s`'s display system.
fn show_for(s: &Session) -> impl Fn(Quantity, f64) -> String + use<> {
    let display = s.units().clone();
    move |q, v| crate::units::shown(&display, crate::units::geometry(), q, v)
}

/// Elements a plan aims at most for one body unless told otherwise.
pub const DEFAULT_MAX_ELEMENTS: usize = 300_000;
/// The most a plan may aim at.
const MAX_ELEMENTS: usize = 5_000_000;
/// Element size by default, as a fraction of the body's bounding-box diagonal (as the meshers).
const DEFAULT_RELATIVE_SIZE: f64 = 1.0 / 25.0;
/// The finest tessellation a plan asks for, as a fraction of the diagonal (as the meshers take).
const MIN_RELATIVE_TOLERANCE: f64 = 1e-5;
/// The share of the body's volume the chords of its curved faces may cost at the plan's size: a
/// linear element's face spans a chord of the surface, and the sag between them is volume lost
/// (or gained, in a hole).
const CHORD_VOLUME_BUDGET: f64 = 0.005;

/// The principal curvatures' magnitudes, larger first, of an analytic surface at `p` (none for a
/// plane or a free-form surface).
fn curvatures(s: &FaceSurface, p: V3) -> Option<(f64, f64)> {
    let off_axis = |origin: V3, axis: V3| {
        let v = std::array::from_fn(|k| p[k] - origin[k]);
        let a = unit(axis)?;
        let along = dot(v, a);
        let r: V3 = std::array::from_fn(|k| v[k] - along * a[k]);
        Some(dot(r, r).sqrt())
    };
    let k = match *s {
        FaceSurface::Plane { .. } => return None,
        FaceSurface::Cylinder { radius, .. } => (1.0 / radius, 0.0),
        FaceSurface::Sphere { radius, .. } => (1.0 / radius, 1.0 / radius),
        FaceSurface::Cone { apex, axis, half_angle } => (half_angle.cos() / off_axis(apex, axis)?, 0.0),
        // Round the tube 1/minor; round the axis cos φ / ρ, with ρ the distance from the axis.
        FaceSurface::Torus { center, axis, major, minor } => {
            let rho = off_axis(center, axis)?;
            (1.0 / minor, ((rho - major) / (minor * rho)).abs())
        }
    };
    let k = (k.0.abs(), k.1.abs());
    (k.0.is_finite() && k.1.is_finite()).then_some(if k.0 >= k.1 { k } else { (k.1, k.0) })
}

/// The largest element size whose chords cost each curved face at most [`CHORD_VOLUME_BUDGET`]
/// of the body's `volume`, with the curved face that needs the finest size and its radius.
///
/// A linear face of size h across a surface of principal curvatures κ1, κ2 sags h²(κ1 + κ2)/12
/// on average below it, so a face of area A costs A h²(κ1 + κ2)/12 of volume. A face's share is
/// bounded two ways, and the looser bound is enough: its cost of all the faces against the body's volume (a
/// small blend costs little whatever its curvature), and per face its cost against the volume it
/// bounds itself, the turning angle per element (a thin tube's wall, whose inner and outer chords nearly
/// cancel, needs no finer size than a solid bar of its radius).
fn curvature_size(surfaces: &[Option<FaceSurface>], mesh: &fieldcraft_geometry::SurfaceMesh, volume: f64) -> Option<(f64, usize, f64)> {
    let n = surfaces.len();
    // Per face: ∫(κ1 + κ2) dA, area, ∫ κ1 (κ1 + κ2) dA, ∫ κ1 dA.
    let mut sums = vec![[0.0f64; 4]; n];
    for (t, &f) in mesh.triangles.iter().zip(&mesh.tri_face) {
        let Some(Some(surface)) = surfaces.get(f as usize) else { continue };
        let [a, b, c] = t.map(|i| mesh.positions.get(i as usize).copied().unwrap_or_default());
        let centroid = std::array::from_fn(|k| (a[k] + b[k] + c[k]) / 3.0);
        let Some((k1, k2)) = curvatures(surface, centroid) else { continue };
        let area = 0.5 * {
            let n = cross(std::array::from_fn(|k| b[k] - a[k]), std::array::from_fn(|k| c[k] - a[k]));
            dot(n, n).sqrt()
        };
        if let Some(s) = sums.get_mut(f as usize) {
            s[0] += area * (k1 + k2);
            s[1] += area;
            s[2] += area * k1 * (k1 + k2);
            s[3] += area * k1;
        }
    }
    let e = CHORD_VOLUME_BUDGET;
    // The whole body's cost against its volume, then the face that needs the finest size by its
    // own turning angle, never finer than the body's bound.
    let total: f64 = sums.iter().map(|s| s[0]).sum();
    let global = (12.0 * e * volume.abs() / total).sqrt();
    let (local, face, radius) = sums
        .iter()
        .enumerate()
        .filter(|(_, s)| s[1] > 0.0 && s[2] > 0.0 && s[3] > 0.0)
        // A circle split into arcs of angle θ loses θ²/6 of its area, a sphere's θ²/3: h² κ1 (κ1 + κ2)/6.
        .map(|(i, s)| ((6.0 * e * s[1] / s[2]).sqrt(), i, s[1] / s[3]))
        .filter(|(h, _, _)| h.is_finite() && *h > 0.0)
        .min_by(|a, b| a.0.total_cmp(&b.0))?;
    let h = global.max(local);
    h.is_finite().then_some((h, face, radius))
}

pub const PLAN: CommandSpec = CommandSpec {
    id: "mesh.plan",
    label: "Mesh Plan",
    params: r#"{"body": "Base", "name"?: "Plan-Base", "size"?: length (largest element size; default 1/25 of the body's diagonal), "layers"?: 2 (elements through the thinnest wall: a shell's thickness, a plate's; 1..=10), "max_elements"?: 300000 (the size grows past the walls' when more would be needed)} — the body's hex-meshing plan derived from the features that built it (geometry.features), saved in the model under `name` (replacing a plan of that name); nothing is meshed (mesh.plan.apply). The plan: the sweep `direction` (and every direction the mapped features are prismatic along), the element `size`, `tolerance` and `layers`; per mapped feature its `template` (sweep of an extrude's profile along its stored direction, block for a box, butterfly for a cylinder or a circular boss, o_grid for a hole or a circular cut, layers through a shell's thickness, imprint of a join, cut or combine tool's cross-section, copy for a pattern or mirror, partition for a split body) and the `faces` it made (the projection targets); the `rings` round the body's circles (o_grid or butterfly: centre, axis, radius, half-side of the square `frame`, the cylindrical `faces` they project onto); the partition `cuts` the frames and imprints need (each with its feature and reason, never on a face of the body); the `mesher` (cells when there are cuts, else sweep); the `copies` (the other bodies a pattern or mirror made from it); the features it does not map with the reason (`unmapped`: fillets, chamfers, revolves, sweeps and lofts are step 3; a feature not prismatic along the sweep needs projection blocking, #399) and `notes`; with unmapped features that the sweep cannot mesh, a `decomposition` proposal (not saved): the partition `cuts` (point, normal) that box those features' faces into one cell, one element size clear of them, inside the body and clear of its faces (`open_sides` where no such cut fits), the boxed `tet_region`, an upper bound of its share of the volume, and the `commands` that make the cuts and mesh the cells with mesh.cells (the boxed cell tet, joined to the hex cells by pyramids)"#,
    form: Some(super::forms::MESH_PLAN),
    kind: CommandKind::Model,
    requires: &[Requirement::Design],
    run: plan,
};

pub const APPLY: CommandSpec = CommandSpec {
    id: "mesh.plan.apply",
    label: "Apply Mesh Plan",
    params: r#"{"plan": "Plan-Base" (a saved plan) | {plan object, as mesh.plan replies it; edited}, "name"?: "Plan-Base" (the name an object is saved under; default Plan-<body>), "symmetry"?: true (a projection plan on a body mirror symmetric across planes through its middle square to the sweep meshes one cell and mirrors it: false meshes the whole body), "corners"?: true (with symmetry: when the cell's two ends along another in-plane axis are alike — the planes of its block layout mirror one another from both ends inward, and its surface there too — the cell is cut midway into the next sheet of blocks and at its image, the low end meshed and mirrored onto the high one, and the middle meshed with both cut faces fixed to those meshes, so both ends mesh exactly alike; reported in `corners`; false meshes the cell whole)} — carries the plan out with today's meshers: replaces the plan's earlier partitions (`<plan>-cut<k>`) with its `cuts` (datum planes and partitions of the plan's body only; a cut the kernel cannot make, or that leaves a cell no mesher takes, is left out and reported in `skipped_cuts`), then meshes the body: its cells swept together with structured hex (mesh.cells, `default` hex) when there are cuts, else the body swept (mesh.hex.generate, swept), at the plan's size; each of its `copies` is meshed the same way by a plan of its own (`Plan-<copy>`, saved). No fallback: a cell that does not mesh all-hex is the error. Replies the cuts made and skipped, the cells, the element count, the smallest scaled Jacobian, `all_hex`, whether the cells joined node to node into one mesh (`conforming`), the mapped features, what today's meshers cannot build yet (`pending`: the O-grid and butterfly rings, which need the projection blocking mesher, #399), the `unmapped` features, the copies and the mesher's reply (`mesh`); `untagged_faces` names each face of the meshed bodies no element face is tagged with (a set by it finds nothing: only a face narrower than the elements beside it); `replaced` names the meshes it replaced. A plan object is checked and saved first"#,
    form: Some(super::forms::MESH_PLAN_APPLY),
    kind: CommandKind::Model,
    requires: &[Requirement::Design],
    run: apply,
};

type V3 = [f64; 3];

fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn unit(v: V3) -> Option<V3> {
    let n = dot(v, v).sqrt();
    (n.is_finite() && n > MIN_NORM).then(|| v.map(|c| c / n))
}

/// Parallel or opposite, within 1e-6 of a radian.
fn parallel(a: V3, b: V3) -> bool {
    match (unit(a), unit(b)) {
        (Some(a), Some(b)) => dot(a, b).abs() > 1.0 - SAME_DIRECTION,
        _ => false,
    }
}

fn perpendicular(a: V3, b: V3) -> bool {
    match (unit(a), unit(b)) {
        (Some(a), Some(b)) => dot(a, b).abs() < SAME_DIRECTION,
        _ => false,
    }
}

/// A direction up to its sign: the first non-zero component positive.
fn canonical(d: V3) -> Option<V3> {
    let d = unit(d)?;
    let first = d.iter().copied().find(|c| c.abs() > NON_ZERO_COMPONENT).unwrap_or(1.0);
    Some(if first < 0.0 { d.map(|c| if c == 0.0 { 0.0 } else { -c }) } else { d })
}

const AXES: [V3; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// Whether the faces `faces` (persistent names on `shape`) are all prismatic along `d`: planes
/// along or across it, cylinders about an axis along it; a merged face on no one surface, when
/// each of its B-rep faces is.
fn faces_prismatic(shape: &Shape, surfaces: &[Option<FaceSurface>], faces: &[String], d: V3) -> bool {
    let prismatic = |s: Option<FaceSurface>| match s {
        Some(FaceSurface::Plane { normal, .. }) => parallel(normal, d) || perpendicular(normal, d),
        Some(FaceSurface::Cylinder { axis, .. }) => parallel(axis, d),
        _ => false,
    };
    faces.iter().all(|name| {
        let Some(i) = shape.face_index(name) else { return true };
        match surfaces.get(i as usize).copied().flatten() {
            Some(s) => prismatic(Some(s)),
            None => shape.member_surfaces(i).is_ok_and(|m| m.len() > 1 && m.into_iter().all(prismatic)),
        }
    })
}

/// The kinds of feature the plan has a template for.
fn template_of(f: &FeatureNode) -> Option<&'static str> {
    let op = f.operation.as_deref().unwrap_or("new_body");
    Some(match (f.kind.as_str(), op) {
        ("extrude", "new_body") => "sweep",
        ("extrude", "join" | "cut" | "intersect") => "imprint",
        ("box", "new_body") => "block",
        ("box", _) => "imprint",
        ("cylinder", "cut") => "o_grid",
        ("cylinder", _) => "butterfly",
        ("hole", _) => "o_grid",
        ("shell", _) => "layers",
        ("combine", _) => "imprint",
        ("pattern" | "mirror", _) => "copy",
        ("split", _) => "partition",
        ("fillet" | "chamfer", _) => "blend",
        _ => return None,
    })
}

/// Why a feature with no template is not mapped.
fn unmapped_reason(f: &FeatureNode) -> String {
    match f.kind.as_str() {
        "revolve" => "revolving the profile mesh about the stored axis is step 3 of the provenance plan (#402); no template yet".into(),
        "sweep" | "loft" => {
            "sweeping the profile along the stored path or between the sections is step 3 of the provenance plan (#402); no template yet".into()
        }
        other => format!("no template for `{other}` features yet"),
    }
}

/// The directions a feature builds along, or is prismatic along by its nature.
fn feature_directions(f: &FeatureNode) -> Vec<V3> {
    match f.kind.as_str() {
        "box" => AXES.to_vec(),
        _ => f.direction.into_iter().collect(),
    }
}

/// Whether `f` keeps to a sweep along `d`.
fn compatible(f: &FeatureNode, d: V3, shape: &Shape, surfaces: &[Option<FaceSurface>]) -> bool {
    match f.kind.as_str() {
        "extrude" | "cylinder" | "hole" if f.direction.is_some_and(|x| parallel(x, d)) => true,
        "box" => AXES.iter().any(|a| parallel(*a, d)),
        // Copies and combines bring the features they repeat or join, which are judged themselves.
        "pattern" | "mirror" | "combine" | "split" => true,
        // A blend is a layer of blocks along the body's sides, whatever its edges' direction.
        "fillet" | "chamfer" => true,
        _ => faces_prismatic(shape, surfaces, &f.faces, d),
    }
}

/// A circle a butterfly or O-grid goes round: in the plane square to the sweep.
#[derive(Debug, Clone, Copy)]
struct Circle {
    center: V3,
    radius: f64,
    feature: u64,
    kind: RingKind,
    /// Whether the circle needs a frame cut round it (not the body's own outline).
    framed: bool,
}

/// The circles of `f` along `d`: a cylinder's or hole's, an extrude's profile circles (in world
/// coordinates). A circle that is the body's whole outline (a lone cylinder) needs no frame.
fn circles(f: &FeatureNode, d: V3) -> Vec<Circle> {
    let mut out = Vec::new();
    let cut = f.operation.as_deref() == Some("cut");
    match f.kind.as_str() {
        "cylinder" | "hole" if f.direction.is_some_and(|x| parallel(x, d)) => {
            let r = f.values.get("radius").copied().or_else(|| f.values.get("diameter").map(|v| v / 2.0));
            if let (Some(c), Some(r)) = (f.placement, r) {
                let kind = if cut { RingKind::OGrid } else { RingKind::Butterfly };
                out.push(Circle { center: c, radius: r, feature: f.id, kind, framed: f.operation.as_deref() != Some("new_body") });
            }
        }
        "extrude" if f.direction.is_some_and(|x| parallel(x, d)) => {
            let Some(sk) = &f.sketch else { return out };
            let world = |p: [f64; 2]| -> V3 { std::array::from_fn(|k| sk.origin[k] + p[0] * sk.x_dir[k] + p[1] * sk.y_dir[k]) };
            for p in &sk.profiles {
                // A profile whose outline is one circle is the body's own when the extrude makes
                // a new body; a joined boss's or a cut's otherwise.
                let outline_is_body = f.operation.as_deref() == Some("new_body");
                // An outline is material inside (unless the extrude cuts), a hole material outside.
                let loops = std::iter::once((&p.outer, outline_is_body, cut)).chain(p.holes.iter().map(|h| (h, false, !cut)));
                for (lp, own, hollow) in loops {
                    if let [Curve::Circle { center, radius }] = lp.as_slice() {
                        let kind = if hollow { RingKind::OGrid } else { RingKind::Butterfly };
                        out.push(Circle { center: world(*center), radius: *radius, feature: f.id, kind, framed: !own });
                    }
                }
            }
        }
        _ => {}
    }
    out
}

/// The circles of `f` where the body has them: a feature's stored placement is the original's,
/// and a pattern or mirror copy has its circles elsewhere. Each stored circle takes the centre of
/// a cylindrical face of the feature on this body of its radius, about an axis along `d`, each
/// face once; a stored circle with no such face keeps its place.
fn placed(stored: Vec<Circle>, f: &FeatureNode, d: V3, shape: &Shape, surfaces: &[Option<FaceSurface>]) -> Vec<Circle> {
    let mut axes: Vec<(V3, f64)> = Vec::new();
    for name in &f.faces {
        let Some(i) = shape.face_index(name) else { continue };
        if let Some(FaceSurface::Cylinder { origin, axis, radius }) = surfaces.get(i as usize).copied().flatten()
            && parallel(axis, d)
        {
            // The axis's point in the plane through the origin square to `d`.
            let along = dot(origin, d);
            let p: V3 = std::array::from_fn(|k| origin[k] - along * d[k]);
            let same = |(q, r): &(V3, f64)| {
                (r - radius).abs() < SAME_CIRCLE * radius.max(1.0) && (0..3).all(|k| (q[k] - p[k]).abs() < SAME_CIRCLE * radius.max(1.0))
            };
            if !axes.iter().any(same) {
                axes.push((p, radius));
            }
        }
    }
    let mut out = Vec::with_capacity(stored.len());
    for mut c in stored {
        let near = |(p, r): &(V3, f64)| {
            let off: V3 = std::array::from_fn(|k| c.center[k] - p[k]);
            let along = dot(off, d);
            ((r - c.radius).abs() < SAME_CIRCLE * c.radius.max(1.0), dot(off, off) - along * along)
        };
        let best = axes.iter().enumerate().filter(|(_, a)| near(a).0).min_by(|a, b| near(a.1).1.total_cmp(&near(b.1).1)).map(|(k, _)| k);
        if let Some(k) = best {
            let (p, _) = axes.swap_remove(k);
            // Keep the stored point's level along the sweep.
            let along = dot(c.center, d);
            c.center = std::array::from_fn(|i| p[i] + along * d[i]);
        }
        out.push(c);
    }
    out
}

/// What the cuts must keep clear of: the body's extent, its planar faces (normal, offset, name)
/// and the points of its curved faces, and the smallest distance a cut keeps.
struct Clearance<'a> {
    corners: Vec<V3>,
    planes: Vec<(V3, f64, &'a str)>,
    curved: Vec<Vec<V3>>,
    gap: f64,
}

impl Clearance<'_> {
    /// The body's extent along `n`.
    fn extent(&self, n: V3) -> (f64, f64) {
        self.corners.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), c| (l.min(dot(*c, n)), h.max(dot(*c, n))))
    }

    /// Whether the plane `n·x = t` lies inside the body, `gap` from its extent.
    fn inside(&self, n: V3, t: f64) -> bool {
        let (lo, hi) = self.extent(n);
        t > lo + self.gap && t < hi - self.gap
    }

    /// Whether the plane passes through a curved face (between its points' extremes along `n`,
    /// by more than `gap`): the partitions' tessellation does not close along such a cut.
    fn crosses_curved(&self, n: V3, t: f64) -> bool {
        let same = length(SAME_LENGTH);
        self.curved.iter().any(|pts| {
            let (lo, hi) = pts.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), p| (l.min(dot(*p, n)), h.max(dot(*p, n))));
            t > lo + same && t < hi - same
        })
    }

    /// The planar faces square to `n` within `gap` of the plane `n·x = t`, with their distance.
    fn faces_near(&self, n: V3, t: f64) -> impl Iterator<Item = (f64, &str)> + '_ {
        self.planes
            .iter()
            .filter(move |(m, _, _)| parallel(*m, n))
            .map(move |(m, off, name)| ((t - off * dot(*m, n).signum()).abs(), *name))
            .filter(|(x, _)| *x < self.gap)
    }
}

/// A joined or cutting tool's place along the sweep and, for a cut, its sides across it: the
/// levels and the frame of its imprint.
struct Imprint {
    feature: u64,
    cut: bool,
    /// Where its cross-section starts and ends along the sweep direction.
    extent: Option<(f64, f64)>,
    /// A cut's side planes along the sweep (normal, offset): a rectangle's four.
    sides: Vec<(V3, f64)>,
    /// The middle of the cut's rectangle: which way is out of the opening from each side.
    center: Option<V3>,
}

/// The imprint of `f` (a join or cut, `cut` saying which) along `d`.
fn imprint_of(f: &FeatureNode, d: V3, cut: bool) -> Option<Imprint> {
    let v = |k: &str| f.values.get(k).copied();
    let mut sides = Vec::new();
    let mut center = None;
    let extent = match f.kind.as_str() {
        "box" => {
            let p = f.placement?;
            let dims = [v("length")?, v("width")?, v("height")?];
            let k = AXES.iter().position(|a| parallel(*a, d))?;
            for j in (0..3).filter(|&j| j != k) {
                sides.push((AXES[j], p[j]));
                sides.push((AXES[j], p[j] + dims[j]));
            }
            center = Some(std::array::from_fn(|i| p[i] + 0.5 * dims[i]));
            Some((p[k], p[k] + dims[k]))
        }
        "cylinder" => {
            let (p, a) = (f.placement?, f.direction?);
            parallel(a, d).then(|| (dot(p, d), dot(p, d) + v("height").unwrap_or_default() * dot(a, d)))
        }
        "hole" => {
            let (p, a) = (f.placement?, f.direction?);
            v("depth").filter(|_| parallel(a, d)).map(|depth| (dot(p, d), dot(p, d) + depth * dot(a, d)))
        }
        "extrude" => {
            let (sk, a) = (f.sketch.as_ref()?, f.direction?);
            if !parallel(a, d) || f.two_sided {
                return None;
            }
            // A rectangle's sides (four lines, each along a sketch axis).
            let same = length(SAME_COORDINATE);
            for p in &sk.profiles {
                let lines: Vec<([f64; 2], [f64; 2])> =
                    p.outer.iter().filter_map(|c| if let Curve::Line { a, b } = c { Some((*a, *b)) } else { None }).collect();
                let square = lines.len() == 4
                    && p.outer.len() == 4
                    && p.holes.is_empty()
                    && lines.iter().all(|(a, b)| (a[0] - b[0]).abs() < same || (a[1] - b[1]).abs() < same);
                if square {
                    let world = |q: [f64; 2]| -> V3 { std::array::from_fn(|k| sk.origin[k] + q[0] * sk.x_dir[k] + q[1] * sk.y_dir[k]) };
                    let mid = lines.iter().fold([0.0; 2], |m, (a, _)| [m[0] + a[0] / 4.0, m[1] + a[1] / 4.0]);
                    center = Some(world(mid));
                    for (a, b) in lines {
                        let (wa, wb) = (world(a), world(b));
                        if let Some(n) = unit(cross(std::array::from_fn(|k| wb[k] - wa[k]), d)) {
                            sides.push((n, dot(n, wa)));
                        }
                    }
                }
            }
            let t0 = dot(sk.origin, d);
            Some((t0, t0 + v("distance")? * dot(a, d)))
        }
        _ => None,
    };
    Some(Imprint { feature: f.id, cut, extent, sides: if cut { sides } else { Vec::new() }, center })
}

/// The cuts that imprint the tools: a plane square to the sweep where each tool's cross-section
/// starts or ends inside the body (a level), and for a cutting tool with a rectangular section,
/// its four sides extended across the body (the frame of the opening). A plane goes in only
/// inside the body (`gap` from its extent), not through a curved face, not near but off a face, not on a face a feature
/// other than the tools made (a level exactly on a tool's own face is that face carried on: no
/// sliver), and not next to another cut. A frame side on a face of the body (an opening's end
/// wall) is moved one element (twice `gap`) out of the opening: the kernel cannot cut along a
/// face (#427, #455), and a strip one element wide is left between the cut and the face.
fn imprint_cuts(tools: &[Imprint], d: V3, body: &Clearance, cuts: &mut Vec<PlanCut>, show: Show) {
    let same = length(SAME_LENGTH);
    let tool_ids: Vec<u64> = tools.iter().map(|t| t.feature).collect();
    // `off`: the distance from `t` of the face a frame side was moved off, allowed exactly.
    let mut add = |n: V3, t: f64, reason: String, feature: u64, off: Option<f64>| {
        if !body.inside(n, t) || body.crosses_curved(n, t) {
            return;
        }
        // Near a face but off it, or on a face some other feature made (a floor a boss stands
        // on): a sliver, or a cut along a face that the kernel splits badly.
        let near_face = body.faces_near(n, t).any(|(x, name)| {
            if let Some(o) = off {
                return (x - o).abs() > same;
            }
            x > same || !fieldcraft_geometry::features::feature_ids(name).any(|id| tool_ids.contains(&id))
        });
        if near_face || cuts.iter().any(|k| parallel(k.normal, n) && (dot(k.point, n) - t).abs() < body.gap) {
            return;
        }
        let role = if parallel(n, d) { CutRole::Level } else { CutRole::Frame };
        cuts.push(PlanCut { point: n.map(|x| x * t), normal: n, role, feature: Some(feature), reason });
    };
    for tool in tools {
        let what = if tool.cut { "cut" } else { "joined" };
        if let Some((a, b)) = tool.extent {
            for (t, end) in [(a.min(b), "starts"), (a.max(b), "ends")] {
                add(d, t, format!("level: the {what} section of F{} {end} here", tool.feature), tool.feature, None);
            }
        }
        for (n, off) in &tool.sides {
            let why = format!("frame: a side of the opening F{} cuts, carried across the body", tool.feature);
            let on_face = body.faces_near(*n, *off).any(|(x, _)| x <= same);
            match tool.center.filter(|_| on_face) {
                Some(c) => {
                    let out = (off - dot(*n, c)).signum();
                    let shift = 2.0 * body.gap;
                    add(
                        *n,
                        off + out * shift,
                        format!("{why}, {} out of the opening off its end face", show(Quantity::Length, shift)),
                        tool.feature,
                        Some(shift),
                    );
                }
                None => add(*n, *off, why, tool.feature, None),
            }
        }
    }
}

/// Most profile segments [`profile_wall`] compares pairwise.
const MAX_WALL_SEGMENTS: usize = 2000;

/// The thinnest wall of a profile: the smallest distance between two parallel lines of its
/// loops with material between them (the point halfway between, over their overlap, inside the
/// outline and outside its holes). Arcs count only as the polygon their ends make.
fn profile_wall(p: &fieldcraft_geometry::ProfileGeometry) -> Option<f64> {
    let polygon = |lp: &[Curve]| -> Vec<[f64; 2]> {
        let mut pts = Vec::new();
        for c in lp {
            match c {
                Curve::Line { a, .. } => pts.push(*a),
                Curve::Arc { center, radius, start, sweep } => {
                    for k in 0..8 {
                        let t = start + sweep * f64::from(k) / 8.0;
                        pts.push([center[0] + radius * t.cos(), center[1] + radius * t.sin()]);
                    }
                }
                Curve::Circle { center, radius } => {
                    for k in 0..16 {
                        let t = std::f64::consts::TAU * f64::from(k) / 16.0;
                        pts.push([center[0] + radius * t.cos(), center[1] + radius * t.sin()]);
                    }
                }
                Curve::Other { .. } => {}
            }
        }
        pts
    };
    let inside = |poly: &[[f64; 2]], q: [f64; 2]| {
        let mut odd = false;
        for (i, a) in poly.iter().enumerate() {
            let b = poly[(i + 1) % poly.len()];
            if (a[1] > q[1]) != (b[1] > q[1]) && q[0] < a[0] + (q[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]) {
                odd = !odd;
            }
        }
        odd
    };
    let outer = polygon(&p.outer);
    let holes: Vec<Vec<[f64; 2]>> = p.holes.iter().map(|h| polygon(h)).collect();
    let material = |q: [f64; 2]| inside(&outer, q) && !holes.iter().any(|h| inside(h, q));
    let lines: Vec<([f64; 2], [f64; 2])> = std::iter::once(&p.outer)
        .chain(&p.holes)
        .flatten()
        .filter_map(|c| if let Curve::Line { a, b } = c { Some((*a, *b)) } else { None })
        .take(MAX_WALL_SEGMENTS)
        .collect();
    let (min_segment, same) = (length(MIN_SEGMENT), length(SAME_COORDINATE));
    let mut best: Option<f64> = None;
    for (i, (a, b)) in lines.iter().enumerate() {
        let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
        let len = (dx * dx + dy * dy).sqrt();
        if !(len > min_segment) {
            continue;
        }
        let (ux, uy) = (dx / len, dy / len);
        for (c, e) in lines.iter().skip(i + 1) {
            let (ex, ey) = (e[0] - c[0], e[1] - c[1]);
            if (ux * ey - uy * ex).abs() > SAME_LINE * (ex * ex + ey * ey).sqrt().max(min_segment) {
                continue;
            }
            let gap = ux * (c[1] - a[1]) - uy * (c[0] - a[0]);
            if gap.abs() < same {
                continue;
            }
            // The overlap of the two along the direction.
            let along = |q: [f64; 2]| (q[0] - a[0]) * ux + (q[1] - a[1]) * uy;
            let (lo, hi) = (along(*c).min(along(*e)).max(0.0), along(*c).max(along(*e)).min(len));
            if hi - lo <= same {
                continue;
            }
            let m = 0.5 * (lo + hi);
            let q = [a[0] + ux * m - uy * gap / 2.0, a[1] + uy * m + ux * gap / 2.0];
            if material(q) && best.is_none_or(|x| gap.abs() < x) {
                best = Some(gap.abs());
            }
        }
    }
    best
}

/// Two in-plane axes square to `d`: the other two principal axes when `d` is one.
fn in_plane_axes(d: V3) -> [V3; 2] {
    if let Some(k) = AXES.iter().position(|a| parallel(*a, d)) {
        let others: Vec<V3> = AXES.iter().enumerate().filter(|(i, _)| *i != k).map(|(_, a)| *a).collect();
        if let [u, v] = others.as_slice() {
            return [*u, *v];
        }
    }
    let least = AXES.iter().copied().min_by(|a, b| dot(*a, d).abs().total_cmp(&dot(*b, d).abs())).unwrap_or([1.0, 0.0, 0.0]);
    let u = unit(cross(d, least)).unwrap_or([1.0, 0.0, 0.0]);
    [u, unit(cross(d, u)).unwrap_or([0.0, 1.0, 0.0])]
}

/// The frame cuts round the circles: for each circle (concentric ones together, the largest
/// framing them) and each in-plane axis, a plane on either side, square to the axis, at a
/// distance `h` from the centre; none on a side where the body's side is within 4 radii (it is
/// the frame there). `h` is first the distance to the body's side on the other side (a square
/// frame, when that is 1.4 to 4 radii), else 2, 1.75, 2.5, 1.5 or 3 radii, the first
/// that keeps the plane inside the body, off its curved faces, at least `gap` off every face of
/// the body square to the axis, clear of the other circles (1.3 of their radius), and not next to another cut (one
/// within 0.6 of the radius serves instead). A side with no such place gets no cut, and a note.
fn frame_cuts(circles: &[Circle], d: V3, body: &Clearance, cuts: &mut Vec<PlanCut>, show: Show) -> Vec<String> {
    let axis_gap = axis_gap();
    let mut groups: Vec<Circle> = Vec::new();
    for c in circles.iter().filter(|c| c.framed) {
        let same_axis = |g: &Circle| {
            let off = std::array::from_fn::<f64, 3, _>(|k| c.center[k] - g.center[k]);
            let along = dot(off, d);
            dot(off, off) - along * along < axis_gap
        };
        match groups.iter_mut().find(|g| same_axis(g)) {
            Some(g) if c.radius > g.radius => *g = *c,
            Some(_) => {}
            None => groups.push(*c),
        }
    }
    let gap = body.gap;
    let mut notes = Vec::new();
    for g in &groups {
        let r = g.radius;
        for a in in_plane_axes(d) {
            let (amin, amax) = body.extent(a);
            let c = dot(g.center, a);
            for sign in [1.0, -1.0] {
                let opposite = if sign > 0.0 { c - amin } else { amax - c };
                let same = if sign > 0.0 { amax - c } else { c - amin };
                // The body's side close by on this side frames the circle itself.
                if same <= 4.0 * r {
                    continue;
                }
                let mut hs: Vec<f64> = [opposite].into_iter().filter(|h| (1.4 * r..=4.0 * r).contains(h)).collect();
                hs.extend([2.0, 1.75, 2.5, 1.5, 3.0].map(|k| k * r));
                let mut placed = false;
                for h in hs {
                    let t = c + sign * h;
                    if !body.inside(a, t) || body.crosses_curved(a, t) || body.faces_near(a, t).next().is_some() {
                        continue;
                    }
                    if groups.iter().any(|o| !std::ptr::eq(o, g) && (t - dot(o.center, a)).abs() < 1.3 * o.radius) {
                        continue;
                    }
                    if cuts.iter().any(|k| parallel(k.normal, a) && (dot(k.point, a) - t).abs() < (0.6 * r).max(gap)) {
                        placed = true;
                        break;
                    }
                    cuts.push(PlanCut {
                        point: a.map(|x| x * t),
                        normal: a,
                        role: CutRole::Frame,
                        feature: Some(g.feature),
                        reason: format!(
                            "frame of the O-grid round the circle of radius {} at F{}: {} from its centre",
                            show(Quantity::Length, r),
                            g.feature,
                            show(Quantity::Length, h)
                        ),
                    });
                    placed = true;
                    break;
                }
                if !placed {
                    notes.push(format!(
                        "F{}: no frame cut on the {} side along {:?} (the body's side or a face is the frame there)",
                        g.feature,
                        if sign > 0.0 { "positive" } else { "negative" },
                        a
                    ));
                }
            }
        }
    }
    notes
}

/// The rings round the circles: each circle once (concentric ones the largest), with the
/// half-side of its square frame (an O-grid's: the nearest frame cut or side of the body square
/// to an in-plane axis; a butterfly's core: half the radius) and the cylindrical faces of the
/// feature it projects onto.
fn plan_rings(
    circles: &[Circle],
    d: V3,
    body: &Clearance,
    cuts: &[PlanCut],
    shape: &Shape,
    surfaces: &[Option<FaceSurface>],
    features: &[FeatureNode],
) -> Vec<PlanRing> {
    let axis_gap = axis_gap();
    let mut out: Vec<PlanRing> = Vec::new();
    for c in circles {
        let same = |r: &PlanRing| {
            let off: V3 = std::array::from_fn(|k| c.center[k] - r.center[k]);
            let along = dot(off, d);
            dot(off, off) - along * along < axis_gap
        };
        // Concentric circles each get a ring (a boss and its hole: squares inside squares).
        if out.iter().any(|r| same(r) && (r.radius - c.radius).abs() < CONCENTRIC * c.radius.max(1.0)) {
            continue;
        }
        let frame = match c.kind {
            RingKind::Butterfly => 0.5 * c.radius,
            RingKind::OGrid => in_plane_axes(d)
                .iter()
                .flat_map(|&a| {
                    let (lo, hi) = body.extent(a);
                    let t = dot(c.center, a);
                    let planes = cuts.iter().filter(move |k| parallel(k.normal, a)).map(move |k| (dot(k.point, a) - t).abs());
                    [t - lo, hi - t].into_iter().chain(planes)
                })
                .filter(|x| *x > c.radius)
                .fold(f64::INFINITY, f64::min),
        };
        let faces = features
            .iter()
            .filter(|f| f.id == c.feature)
            .flat_map(|f| &f.faces)
            .filter(|name| {
                let Some(i) = shape.face_index(name) else { return false };
                match surfaces.get(i as usize).copied().flatten() {
                    Some(FaceSurface::Cylinder { origin, axis, radius }) => {
                        let off: V3 = std::array::from_fn(|k| c.center[k] - origin[k]);
                        let along = dot(off, d);
                        parallel(axis, d)
                            && (radius - c.radius).abs() < SAME_CIRCLE * c.radius.max(1.0)
                            && dot(off, off) - along * along < ON_AXIS * c.radius.max(1.0)
                    }
                    _ => false,
                }
            })
            .cloned()
            .collect();
        if let Some(axis) = canonical(d) {
            out.push(PlanRing {
                feature: c.feature,
                kind: c.kind,
                center: c.center,
                axis,
                radius: c.radius,
                frame: if frame.is_finite() { frame } else { 2.0 * c.radius },
                faces,
            });
        }
    }
    out
}

/// The other bodies of the design the plan's copy features (a pattern, a mirror) made from the
/// same timeline: each lists that feature too.
fn copies_of(design: &Design, body: &str, features: &[PlanFeature]) -> Result<Vec<String>> {
    let copy_ids: Vec<u64> = features.iter().filter(|f| f.template == "copy").map(|f| f.id).collect();
    if copy_ids.is_empty() {
        return Ok(Vec::new());
    }
    let mut out = Vec::new();
    for graph in design.features(None)? {
        if graph.body != body && graph.features.iter().any(|f| copy_ids.contains(&f.id)) {
            out.push(graph.body);
        }
    }
    Ok(out)
}

/// The shape of `body`: a body of the design as opened, or a cell partitions cut from one.
fn plan_shape<'a>(design: &'a Design, body: &str) -> Option<&'a Shape> {
    design.source_body(body).or_else(|| design.body(body))
}

/// Build the plan of `body`: a body of the design as opened, or a cell of one (planned from the
/// features that built the body it was cut from, with the rings the cell holds).
pub(crate) fn build(design: &Design, body: &str, size: Option<f64>, layers: u32, max_elements: usize, show: Show) -> Result<MeshPlan> {
    let length = |v: f64| show(Quantity::Length, v);
    let shape = plan_shape(design, body).ok_or_else(|| EngineError::Params(format!("no body `{body}` in the design")))?;
    let source = design.source_of(body).unwrap_or(body);
    let is_cell = source != body;
    let graph = design.features(Some(source))?;
    let features: Vec<FeatureNode> = graph.into_iter().next().map(|b| b.features).unwrap_or_default();
    let surfaces = shape.face_surfaces()?;
    let m = shape.measure()?;
    let (lo, hi) = (m.bbox_min, m.bbox_max);
    let diagonal = dot(std::array::from_fn(|k| hi[k] - lo[k]), std::array::from_fn(|k| hi[k] - lo[k])).sqrt();
    if !(diagonal.is_finite() && diagonal > 0.0) {
        return Err(EngineError::Params(format!("body `{body}` has no extent")));
    }

    // The sweep direction: the candidate most mapped features keep to, the first new body's own
    // on a tie (candidates in the order the features give them).
    let mapped: Vec<&FeatureNode> = features.iter().filter(|f| template_of(f).is_some() && f.error.is_none()).collect();
    let primary = mapped.iter().find(|f| f.operation.as_deref() == Some("new_body")).and_then(|f| feature_directions(f).first().copied());
    let mut candidates: Vec<V3> = Vec::new();
    for d in primary.into_iter().chain(mapped.iter().flat_map(|f| feature_directions(f))).filter_map(canonical) {
        if !candidates.iter().any(|c| parallel(*c, d)) {
            candidates.push(d);
        }
    }
    let score = |d: V3| mapped.iter().filter(|f| compatible(f, d, shape, &surfaces)).count();
    let best = candidates.iter().map(|d| score(*d)).max().unwrap_or(0);
    let directions: Vec<V3> = candidates.iter().copied().filter(|d| score(*d) == best && best > 0).collect();
    let direction = directions.first().copied();

    let mut plan_features = Vec::new();
    let mut unmapped = Vec::new();
    let mut notes = Vec::new();
    let mut circle_list = Vec::new();
    let mut imprints = Vec::new();
    // The thinnest wall the features give, with what gives it.
    let mut thin: Option<(f64, String)> = None;
    let mut thinner = |t: f64, why: String| {
        if t.is_finite() && t > 0.0 && thin.as_ref().is_none_or(|(x, _)| t < *x) {
            thin = Some((t, why));
        }
    };
    for f in &features {
        let unmap = |reason: String| PlanUnmapped { id: f.id, name: f.name.clone(), kind: f.kind.clone(), reason };
        if let Some(e) = &f.error {
            unmapped.push(unmap(format!("failed to rebuild: {e}")));
            continue;
        }
        let Some(template) = template_of(f) else {
            let mut reason = unmapped_reason(f);
            if let Some(d) = direction
                && !f.faces.is_empty()
                && faces_prismatic(shape, &surfaces, &f.faces, d)
            {
                reason.push_str("; its faces are prismatic along the sweep, so the sweep meshes them as part of the cross-section");
            }
            unmapped.push(unmap(reason));
            continue;
        };
        let Some(d) = direction else {
            unmapped.push(unmap("no sweep direction keeps to the body's features".into()));
            continue;
        };
        if !compatible(f, d, shape, &surfaces) {
            unmapped.push(unmap(format!("not prismatic along the sweep direction {d:?}: it needs blocks of its own (projection blocking, #399)")));
            continue;
        }
        let v = |k: &str| f.values.get(k).copied();
        let note = match (f.kind.as_str(), template) {
            ("extrude", "sweep") => {
                if let Some(t) = v("distance").filter(|_| f.direction.is_some_and(|x| parallel(x, d))) {
                    thinner(t, format!("{} extent", f.name));
                }
                for p in f.sketch.iter().flat_map(|sk| &sk.profiles) {
                    if let Some(t) = profile_wall(p) {
                        thinner(t, format!("{} profile's thinnest wall", f.name));
                    }
                }
                let along = if f.direction.is_some_and(|x| parallel(x, d)) {
                    "along its stored direction"
                } else {
                    "across its extrude direction (its profile is prismatic along the sweep)"
                };
                format!("the profile's cross-section swept {along}; circles butterflied, circular holes O-gridded")
            }
            ("extrude", _) => "its profile imprinted on the swept cross-section between the levels of its extent".into(),
            ("box", t) => {
                let dims: Vec<f64> = ["length", "width", "height"].iter().filter_map(|k| v(k)).collect();
                if let Some(min) = dims.iter().copied().reduce(f64::min) {
                    thinner(min, format!("{} thinnest side", f.name));
                }
                if t == "block" {
                    "a mapped block (part of the swept cross-section)".into()
                } else {
                    "its rectangle imprinted on the swept cross-section".into()
                }
            }
            ("cylinder", "butterfly") => format!("a butterfly round its circle (radius {})", length(v("radius").unwrap_or_default())),
            ("cylinder", _) => format!("an O-grid round its circle (radius {})", length(v("radius").unwrap_or_default())),
            ("hole", _) => format!("an O-grid round the hole (diameter {}) along its stored axis", length(v("diameter").unwrap_or_default())),
            ("shell", _) => {
                let t = v("thickness").unwrap_or_default();
                thinner(t, format!("{} thickness", f.name));
                format!("{layers} layers through the {} wall", length(t))
            }
            ("combine", _) => format!("the cross-sections of {} imprinted on the target", f.tools.join(", ")),
            ("pattern" | "mirror", _) => format!("the instance meshed and repeated by {} transforms", f.transforms.len()),
            ("split", _) => "the pieces meshed as cells of the split".into(),
            ("fillet", _) => {
                format!("a blend layer of blocks along its edges, {} deep, projected onto the fillet faces", length(v("radius").unwrap_or_default()))
            }
            ("chamfer", _) => {
                format!(
                    "a blend layer of blocks along its edges, {} deep, projected onto the chamfer faces",
                    length(v("distance").unwrap_or_default())
                )
            }
            _ => String::new(),
        };
        circle_list.extend(placed(circles(f, d), f, d, shape, &surfaces));
        // A tool: joined or cut by its own operation, or the tool body of a combine.
        let role = match f.operation.as_deref() {
            Some("join") => Some(false),
            Some("cut") => Some(true),
            _ => features.iter().find(|c| c.kind == "combine" && c.tools.iter().any(|t| f.creates.contains(t))).and_then(|c| {
                match c.operation.as_deref() {
                    Some("join") => Some(false),
                    Some("cut") => Some(true),
                    _ => None,
                }
            }),
        };
        let template = if role.is_some() && template == "block" { "imprint" } else { template };
        if let Some(cut) = role {
            imprints.extend(imprint_of(f, d, cut));
        }
        plan_features.push(PlanFeature {
            id: f.id,
            name: f.name.clone(),
            kind: f.kind.clone(),
            template: template.into(),
            note,
            faces: f.faces.clone(),
        });
    }

    // The element size: at most the one asked for (or the meshers' default), small enough for
    // `layers` through the thinnest wall, but not so small the body takes more than
    // `max_elements`.
    let asked = size.unwrap_or(DEFAULT_RELATIVE_SIZE * diagonal);
    let mut h = asked;
    if let Some((t, why)) = &thin {
        let through = t / f64::from(layers);
        if through < h {
            notes.push(format!("size {}: {layers} layers through the {} wall ({why})", length(through), length(*t)));
            h = through;
        }
    }
    // The curved faces: fine enough that the chords keep the volume.
    let chords = shape.tessellate(DEFAULT_RELATIVE_TOLERANCE * diagonal)?;
    if let Some((c, face, r)) = curvature_size(&surfaces, &chords, m.volume)
        && c < h
    {
        let name = shape.faces().get(face).map_or_else(|| format!("face {face}"), |f| f.name.clone());
        notes.push(format!(
            "size {}: the chords of {name} (radius ~{}) cost at most {:.1} % of the volume",
            length(c),
            length(r),
            100.0 * CHORD_VOLUME_BUDGET
        ));
        h = c;
    }
    let floor = (m.volume.abs() / max_elements as f64).cbrt();
    if h < floor {
        notes.push(format!(
            "size {} instead of {}: at most {max_elements} elements (fewer layers through the thinnest wall)",
            length(floor),
            length(h)
        ));
        h = floor;
    }
    let mut cuts = Vec::new();
    let mut rings = Vec::new();
    // Where the frame cuts could not go: the swept cells' business only.
    let mut frame_notes = Vec::new();
    if let Some(d) = direction {
        let planes = surfaces
            .iter()
            .zip(shape.faces())
            .filter_map(|(s, f)| match s {
                Some(FaceSurface::Plane { normal, point }) => Some((*normal, dot(*normal, *point), f.name.as_str())),
                _ => None,
            })
            .collect();
        // The curved faces' points, from a tessellation at the meshers' default tolerance.
        let mesh = shape.tessellate(DEFAULT_RELATIVE_TOLERANCE * diagonal)?;
        let mut curved: Vec<Vec<V3>> = vec![Vec::new(); surfaces.len()];
        for (t, &f) in mesh.triangles.iter().zip(&mesh.tri_face) {
            if !matches!(surfaces.get(f as usize), Some(Some(FaceSurface::Plane { .. })))
                && let Some(pts) = curved.get_mut(f as usize)
            {
                pts.extend(t.iter().filter_map(|&i| mesh.positions.get(i as usize).copied()));
            }
        }
        curved.retain(|p| !p.is_empty());
        let corners = (0..8)
            .map(|i| [if i & 1 == 0 { lo[0] } else { hi[0] }, if i & 2 == 0 { lo[1] } else { hi[1] }, if i & 4 == 0 { lo[2] } else { hi[2] }])
            .collect();
        let body = Clearance { corners, planes, curved, gap: 0.5 * h };
        imprint_cuts(&imprints, d, &body, &mut cuts, show);
        frame_notes = frame_cuts(&circle_list, d, &body, &mut cuts, show);
        rings = plan_rings(&circle_list, d, &body, &cuts, shape, &surfaces, &features);
        // A cell holds the rings whose axis runs through it.
        if is_cell {
            let same = crate::units::geometry_amount(Quantity::Length, SAME_COORDINATE);
            rings.retain(|r| {
                (0..3).all(|k| {
                    (r.center[k] - d[k] * dot(r.center, d)) >= lo[k] - d[k] * dot(lo, d) - same
                        && (r.center[k] - d[k] * dot(r.center, d)) <= hi[k] - d[k] * dot(hi, d) + same
                })
            });
        }
    }
    // The tessellation: the meshers' default, finer where the body has small curves (a hundredth
    // of the smallest radius), so short blend edges stay edges.
    let radius = surfaces
        .iter()
        .filter_map(|s| match s {
            Some(FaceSurface::Cylinder { radius, .. }) | Some(FaceSurface::Sphere { radius, .. }) => Some(*radius),
            Some(FaceSurface::Torus { minor, .. }) => Some(*minor),
            _ => None,
        })
        .fold(f64::INFINITY, f64::min);
    let tolerance = (DEFAULT_RELATIVE_TOLERANCE * diagonal).min(0.01 * radius).max(MIN_RELATIVE_TOLERANCE * diagonal);
    if tolerance < DEFAULT_RELATIVE_TOLERANCE * diagonal {
        notes.push(format!("tolerance {}: a hundredth of the smallest radius ({})", length(tolerance), length(radius)));
    }
    let mesher = if cuts.is_empty() { PlanMesher::Sweep } else { PlanMesher::Cells };
    // A cell's copies are its body's (meshed with it), not its own.
    let copies = if !is_cell { copies_of(design, body, &plan_features)? } else { Vec::new() };
    let mut plan = MeshPlan {
        body: body.into(),
        size: h,
        tolerance: Some(tolerance),
        layers,
        mesher,
        direction,
        directions,
        features: plan_features,
        cuts,
        rings,
        copies,
        unmapped,
        notes,
    };
    plan.check()?;
    // Rings and blends: one block topology over the whole body when it fits (no cuts needed);
    // else the frame cuts and swept cells, with the reason.
    match layout_of(design, &plan, show)? {
        Some(Ok(layout)) => {
            plan.mesher = PlanMesher::Project;
            plan.cuts.clear();
            // The rings' frames are the layout's.
            if let Some(d) = plan.direction {
                for r in &mut plan.rings {
                    if let Some(half) = layout.frame_of(r.center, d).filter(|_| r.kind == RingKind::OGrid) {
                        r.frame = half;
                    }
                }
            }
            plan.notes.push(layout.summary(show));
        }
        Some(Err(why)) => {
            plan.notes.extend(frame_notes);
            plan.notes.push(format!("no projection layout ({why}): the frame cuts box the rings and the cells are swept"));
        }
        None => plan.notes.extend(frame_notes),
    }
    plan.check()?;
    Ok(plan)
}

/// The block layout of `plan` for the projection mesher, when the plan has rings or blends
/// (`None` when it has neither: the sweep meshes it); `Err` inside says why the body does not
/// take one.
fn layout_of(design: &Design, plan: &MeshPlan, show: Show) -> Result<Option<std::result::Result<blocks::Layout, String>>> {
    let blends: Vec<&PlanFeature> = plan.features.iter().filter(|f| f.template == "blend").collect();
    if plan.rings.is_empty() && blends.is_empty() {
        return Ok(None);
    }
    let Some(direction) = plan.direction else { return Ok(Some(Err("no sweep direction".into()))) };
    let shape = plan_shape(design, &plan.body).ok_or_else(|| EngineError::Params(format!("no body `{}` in the design", plan.body)))?;
    let surfaces = shape.face_surfaces()?;
    let m = shape.measure()?;
    let diagonal = (0..3).map(|k| (m.bbox_max[k] - m.bbox_min[k]).powi(2)).sum::<f64>().sqrt();
    let mesh = shape.tessellate(plan.tolerance.unwrap_or(DEFAULT_RELATIVE_TOLERANCE * diagonal).max(1e-5 * diagonal))?;
    let blend_faces = blends.iter().flat_map(|f| &f.faces).filter_map(|n| shape.face_index(n)).collect();
    let shell = if plan.features.iter().any(|f| f.kind == "shell") {
        design
            .features(Some(design.source_of(&plan.body).unwrap_or(&plan.body)))?
            .into_iter()
            .flat_map(|b| b.features)
            .filter(|f| f.kind == "shell")
            .find_map(|f| f.values.get("thickness").copied())
    } else {
        None
    };
    let input =
        blocks::Input { shape, surfaces: &surfaces, mesh: &mesh, direction, rings: &plan.rings, blends: blend_faces, shell, size: plan.size, show };
    // The grid layout, else the ring-and-pillow layout (a section at an angle to the axes); the
    // grid's reason when neither fits.
    Ok(Some(blocks::layout(&input).or_else(|why| pillow::layout(&input).map_err(|other| format!("{why}; ring and pillow: {other}")))))
}

fn plan(s: &mut Session, p: &Value) -> Result<Value> {
    let body: String = super::req(p, "body")?;
    let size: Option<f64> = super::opt(p, "size")?;
    if size.is_some_and(|v| !(v.is_finite() && v > 0.0)) {
        return Err(EngineError::param("size", ProblemKind::NotPositive));
    }
    let layers: u32 = super::opt(p, "layers")?.unwrap_or(2);
    if !(1..=10).contains(&layers) {
        return Err(EngineError::Params("`layers` must be 1 to 10".into()));
    }
    let max_elements: usize = super::opt(p, "max_elements")?.unwrap_or(DEFAULT_MAX_ELEMENTS);
    if !(1..=MAX_ELEMENTS).contains(&max_elements) {
        return Err(EngineError::Params(format!("`max_elements` must be 1 to {MAX_ELEMENTS}")));
    }
    let name: String = super::opt(p, "name")?.unwrap_or_else(|| MeshPlan::default_name(&body));
    let show = show_for(s);
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    let plan = build(design, &body, size, layers, max_elements, &show)?;
    s.model.set_mesh_plan(&name, plan.clone())?;
    let decomposition = decompose::proposal(design, &plan)?;
    let mut reply = serde_json::to_value(&plan).map_err(|e| EngineError::Params(e.to_string()))?;
    if let Some(o) = reply.as_object_mut() {
        o.insert("name".into(), json!(name));
        if let Some(mut d) = decomposition {
            // Its commands are for the user to run: their numbers in the display system.
            for c in d.get_mut("commands").and_then(Value::as_array_mut).into_iter().flatten() {
                if let Some(id) = c.get("command").and_then(Value::as_str).map(str::to_string) {
                    let shown = crate::units::geometry_params_shown(s, &id, &c["params"])?;
                    c["params"] = shown;
                }
            }
            o.insert("decomposition".into(), d);
        }
    }
    Ok(reply)
}

/// The partitions (and their datums) an earlier apply of plan `name` made.
fn plan_cut_names(s: &Session, name: &str) -> Vec<String> {
    let prefix = format!("{name}-cut");
    s.model
        .partitions()
        .iter()
        .map(|p| p.name.clone())
        .filter(|n| n.strip_prefix(&prefix).is_some_and(|k| !k.is_empty() && k.bytes().all(|b| b.is_ascii_digit())))
        .collect()
}

fn apply(s: &mut Session, p: &Value) -> Result<Value> {
    let (name, plan) = match p.get("plan") {
        Some(Value::String(n)) => (n.clone(), s.model.mesh_plan(n).cloned().ok_or_else(|| EngineError::Params(format!("no mesh plan `{n}`")))?),
        Some(v @ Value::Object(_)) => {
            let plan: MeshPlan = serde_json::from_value(v.clone()).map_err(|e| EngineError::Params(format!("`plan`: {e}")))?;
            let name: String = super::opt(p, "name")?.unwrap_or_else(|| MeshPlan::default_name(&plan.body));
            s.model.set_mesh_plan(&name, plan.clone())?;
            (name, plan)
        }
        _ => return Err(EngineError::Params("`plan` must be a saved plan's name or a plan object".into())),
    };
    let symmetry: bool = super::opt(p, "symmetry")?.unwrap_or(true);
    let corners: bool = super::opt(p, "corners")?.unwrap_or(true);
    let symmetric = if symmetry { apply_symmetric(s, &name, &plan, corners)? } else { None };
    let mut reply = match symmetric {
        Some(r) => r,
        None => apply_one(s, &name, &plan)?,
    };
    // The copies: each meshed by a plan of its own at this plan's size, saved beside it.
    let mut copies = Vec::with_capacity(plan.copies.len());
    for copy in &plan.copies {
        fieldcraft_model::progress::stage(&format!("Plan: meshing the copy {copy}"), None, "");
        fieldcraft_model::progress::check()?;
        let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
        if design.source_body(copy).is_none() {
            return Err(EngineError::Params(format!("mesh plan `{name}`: no copy `{copy}` in the design")));
        }
        let show = show_for(s);
        let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
        let mut theirs = build(design, copy, Some(plan.size), plan.layers, MAX_ELEMENTS, &show)?;
        // Its copies are this body and the others: they are meshed here, not again.
        theirs.copies.clear();
        theirs.size = plan.size;
        let copy_name = MeshPlan::default_name(copy);
        s.model.set_mesh_plan(&copy_name, theirs.clone())?;
        let r = apply_one(s, &copy_name, &theirs)?;
        take_replaced(&mut reply, &r);
        copies.push(json!({"body": copy, "plan": copy_name, "cells": r["cells"], "elements": r["elements"], "min_scaled_jacobian": r["min_scaled_jacobian"], "conforming": r["conforming"]}));
    }
    let untagged = untagged_faces(s, &plan_bodies(s, p)?);
    let mins: Vec<f64> = std::iter::once(&reply).chain(&copies).filter_map(|r| r["min_scaled_jacobian"].as_f64()).collect();
    if let Some(o) = reply.as_object_mut() {
        if !copies.is_empty() {
            o.insert("copies".into(), json!(copies));
            o.insert("min_scaled_jacobian".into(), json!(mins.iter().copied().reduce(f64::min)));
        }
        o.insert("untagged_faces".into(), json!(untagged));
        let results: Vec<&str> = s.solved_steps().collect();
        if o.get("replaced").is_some_and(|r| r.as_array().is_some_and(|a| !a.is_empty())) && !results.is_empty() {
            o.insert("results_cleared".into(), json!(results));
        }
    }
    Ok(reply)
}

/// The faces of `bodies` (as the design has them) that no element face of the mesh is tagged
/// with: a set, load or relaxation by such a face finds nothing. Only a face narrower than the
/// elements beside it should be one.
pub fn untagged_faces(s: &Session, bodies: &[String]) -> Vec<Value> {
    let c = s.model.classification();
    let mut tagged: BTreeSet<u32> = BTreeSet::new();
    for t in c.elements().values() {
        tagged.extend(t.faces.values().copied());
    }
    let tagged: BTreeSet<(&str, &str)> = tagged
        .iter()
        .filter_map(|&t| c.entity(t))
        .filter_map(|g| match &g.entity {
            Entity::Face(name) => Some((g.body.as_str(), name.as_str())),
            _ => None,
        })
        .collect();
    let Some(design) = s.design.as_deref() else { return Vec::new() };
    bodies
        .iter()
        .filter_map(|b| design.body(b))
        .flat_map(|b| b.faces().iter().map(move |f| (b.name.as_str(), f.name.as_str())))
        .filter(|key| !tagged.contains(key))
        .map(|(body, face)| json!({"body": body, "face": face}))
        .collect()
}

/// The bodies `mesh.plan.apply` fills (replace.rs): the cells of the plan's body and of its
/// copies, as the design has them now (the body itself when it is not cut).
pub(crate) fn plan_bodies(s: &Session, p: &Value) -> Result<Vec<String>> {
    let plan = match p.get("plan") {
        Some(Value::String(n)) => s.model.mesh_plan(n).cloned().ok_or_else(|| EngineError::Params(format!("no mesh plan `{n}`")))?,
        Some(v @ Value::Object(_)) => serde_json::from_value(v.clone()).map_err(|e| EngineError::Params(format!("`plan`: {e}")))?,
        _ => return Err(EngineError::Params("`plan` must be a saved plan's name or a plan object".into())),
    };
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    let sources: Vec<&str> = std::iter::once(plan.body.as_str()).chain(plan.copies.iter().map(String::as_str)).collect();
    Ok(design.bodies.iter().filter(|b| design.source_of(&b.name).is_some_and(|x| sources.contains(&x))).map(|b| b.name.clone()).collect())
}

/// Move a mesher's `replaced` and `stale_sets` (replace.rs) into `into`, which reports them once
/// for the whole command.
fn take_replaced(into: &mut Value, from: &Value) {
    for key in ["replaced", "stale_sets"] {
        let Some(items) = from.get(key).and_then(Value::as_array).filter(|a| !a.is_empty()) else { continue };
        if let Some(o) = into.as_object_mut() {
            let list = o.entry(key).or_insert_with(|| json!([]));
            if let Some(list) = list.as_array_mut() {
                for item in items {
                    if !list.contains(item) {
                        list.push(item.clone());
                    }
                }
            }
        }
    }
}

/// Carry out one plan: its cuts, then the mesh of the cells (or of the body when it has none).
fn apply_one(s: &mut Session, name: &str, plan: &MeshPlan) -> Result<Value> {
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    let source = plan_shape(design, &plan.body)
        .ok_or_else(|| EngineError::Params(format!("the plan is of body `{}`, which the design does not have", plan.body)))?;
    // The chordal tolerance the cells are tessellated with, from the body's extent (as mesh.cells
    // takes it from the cells'), for the check of the cells and the mesher alike.
    let tolerance = match plan.tolerance {
        Some(t) => t,
        None => {
            let b = source.measure()?;
            DEFAULT_RELATIVE_TOLERANCE * (0..3).map(|k| (b.bbox_max[k] - b.bbox_min[k]).powi(2)).sum::<f64>().sqrt()
        }
    };
    if plan.mesher == PlanMesher::Project {
        return apply_projection(s, name, plan, tolerance);
    }
    let cuts: Vec<&PlanCut> = plan.cuts.iter().collect();
    let (cleared, made, skipped) = cut(s, name, &plan.body, &cuts, tolerance)?;
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    let cells: Vec<String> = design.bodies.iter().filter(|b| design.source_of(&b.name) == Some(plan.body.as_str())).map(|b| b.name.clone()).collect();
    // No fallback to a coarser plan: a cell the hex meshers cannot take is the error, with the
    // cuts that made it, so the plan is edited rather than quietly paved.
    let mesh = if cells.len() > 1 {
        (super::cells::CELLS.run)(s, &json!({"bodies": cells, "size": plan.size, "default": "hex", "tolerance": tolerance}))
    } else {
        let body = cells.first().cloned().unwrap_or_else(|| plan.body.clone());
        (super::hex::HEX_GENERATE.run)(s, &json!({"body": body, "size": plan.size, "method": "swept", "tolerance": tolerance}))
    }
    .map_err(|e| {
        if e.is_cancelled() {
            e
        } else {
            EngineError::Params(format!("mesh plan `{name}`: the cells {cells:?} (cuts {made:?}) do not mesh all-hex: {e}"))
        }
    })?;
    let summary = summary(&mesh, &cells);
    // What today's meshers cannot express yet: the rings round circles are built by the
    // projection mesher (#399); until then the frame cuts box the circle and the sweep meshes the
    // cell's cross-section.
    let pending: Vec<Value> = plan
        .rings
        .iter()
        .map(|r| {
            let what = match r.kind {
                RingKind::OGrid => "O-grid ring round the hole",
                RingKind::Butterfly => "butterfly round the circle",
            };
            json!({"feature": r.feature, "ring": r.kind, "radius": r.radius, "needs": format!("the {what} is built by the projection blocking mesher (#399); today the frame cuts box it and the swept cell's cross-section is meshed by the sweep's quad mesher")})
        })
        .collect();
    let mut reply = json!({
        "name": name,
        "body": plan.body,
        "cuts": made,
        "skipped_cuts": skipped.iter().map(|(n, why)| json!({"name": n, "reason": why})).collect::<Vec<_>>(),
        "cells": cells,
        "size": plan.size,
        "tolerance": tolerance,
        "elements": summary.elements,
        "min_scaled_jacobian": summary.min_scaled_jacobian,
        "all_hex": summary.all_hex,
        "conforming": summary.conforming,
        "mapped": plan.features.iter().map(|f| json!({"id": f.id, "name": f.name, "type": f.kind, "template": f.template})).collect::<Vec<_>>(),
        "pending": pending,
        "unmapped": plan.unmapped,
        "mesh": mesh,
    });
    reply = cleared.report_if_any(s, reply);
    if let Some(m) = reply.get_mut("mesh").and_then(Value::as_object_mut) {
        let theirs = json!({"replaced": m.remove("replaced"), "stale_sets": m.remove("stale_sets")});
        m.remove("results_cleared");
        take_replaced(&mut reply, &theirs);
    }
    if let Some(o) = reply.as_object_mut() {
        // Reported once for the whole command, by apply.
        o.remove("results_cleared");
    }
    Ok(reply)
}

/// Most surface points the symmetry test checks.
const SYMMETRY_SAMPLES: usize = 3000;

/// The planes the body of `plan` is mirror symmetric about, square to the in-plane principal axes
/// of the sweep `d` through the middle of its extent (normal, offset), none through a ring: every
/// point of its surface (the middles of its tessellation's triangles, a little inside and a little
/// outside along the triangle's normal) has its mirror image inside and outside alike (the
/// generalized winding number of the tessellation).
fn symmetry_planes(shape: &Shape, d: V3, rings: &[PlanRing]) -> Result<Vec<(V3, f64)>> {
    let m = shape.measure()?;
    let (lo, hi) = (m.bbox_min, m.bbox_max);
    let Some(mirror) = Mirror::of(shape)? else { return Ok(Vec::new()) };
    let mut out = Vec::new();
    for a in in_plane_axes(d).into_iter().filter(|a| AXES.iter().any(|x| parallel(*x, *a))) {
        let c = 0.5 * (dot(lo, a) + dot(hi, a));
        if rings.iter().any(|r| (dot(r.center, a) - c).abs() < 1.01 * r.radius) {
            continue;
        }
        if mirror.symmetric(a, c, |_| true) {
            out.push((a, c));
        }
    }
    Ok(out)
}

/// The mirror test of a body's surface: points of its faces (the middles of its tessellation's
/// triangles, a little inside and a little outside along the triangle's normal) whose mirror
/// images are inside and outside alike (the generalized winding number of the tessellation).
struct Mirror {
    mesh: fieldcraft_geometry::SurfaceMesh,
    /// Each triangle's middle and normal: a point of a face, away from its edges.
    samples: Vec<(V3, V3)>,
    delta: f64,
}

impl Mirror {
    /// `None` for a body with no extent.
    fn of(shape: &Shape) -> Result<Option<Mirror>> {
        let m = shape.measure()?;
        let diagonal = (0..3).map(|k| (m.bbox_max[k] - m.bbox_min[k]).powi(2)).sum::<f64>().sqrt();
        if !(diagonal.is_finite() && diagonal > 0.0) {
            return Ok(None);
        }
        // Fine enough that the chords lie well within `delta` of the surface.
        let mesh = shape.tessellate(2e-4 * diagonal)?;
        let samples: Vec<(V3, V3)> = mesh
            .triangles
            .iter()
            .filter_map(|t| {
                let [a, b, c] = t.map(|i| mesh.positions.get(i as usize).copied().unwrap_or_default());
                let n = unit(cross(std::array::from_fn(|k| b[k] - a[k]), std::array::from_fn(|k| c[k] - a[k])))?;
                Some((std::array::from_fn(|k| (a[k] + b[k] + c[k]) / 3.0), n))
            })
            .collect();
        Ok(Some(Mirror { mesh, samples, delta: 2e-3 * diagonal }))
    }

    /// Whether the samples `keep` takes (at most [`SYMMETRY_SAMPLES`] of them, evenly) have
    /// their mirror images across the plane `x·a = c` inside and outside alike.
    fn symmetric(&self, a: V3, c: f64, keep: impl Fn(V3) -> bool) -> bool {
        let kept: Vec<&(V3, V3)> = self.samples.iter().filter(|(p, _)| keep(*p)).collect();
        let step = (kept.len() / SYMMETRY_SAMPLES).max(1);
        let mirror = |x: V3| -> V3 { std::array::from_fn(|k| x[k] - 2.0 * (dot(x, a) - c) * a[k]) };
        let delta = self.delta;
        kept.iter().step_by(step).all(|(p, n)| {
            let inside = mirror(std::array::from_fn(|k| p[k] - delta * n[k]));
            let outside = mirror(std::array::from_fn(|k| p[k] + delta * n[k]));
            blocks::winding(&self.mesh, inside) > 0.5 && blocks::winding(&self.mesh, outside) < 0.5
        })
    }
}

/// Carry out a projection plan on a mirror-symmetric body by its symmetry: the body cut by its
/// planes of symmetry (`<plan>-cut<k>`), the cell on the negative side of every plane planned
/// and meshed (`Plan-<cell>`), and its mesh mirrored across each plane in turn onto the other
/// cells, joined node to node on the planes — so mirror-image regions mesh alike. `None` when
/// the body has no symmetry, or its cell does not take a projection layout (the whole body is
/// meshed then).
fn apply_symmetric(s: &mut Session, name: &str, plan: &MeshPlan, corners: bool) -> Result<Option<Value>> {
    if plan.mesher != PlanMesher::Project || !plan.copies.is_empty() {
        return Ok(None);
    }
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    let (Some(shape), Some(d)) = (design.source_body(&plan.body), plan.direction) else { return Ok(None) };
    // A body other partitions have cut already is its cells (the whole body's plan is refused).
    if plan_cut_names(s, name).is_empty() && design.body(&plan.body).is_none() {
        return Ok(None);
    }
    let planes = symmetry_planes(shape, d, &plan.rings)?;
    if planes.is_empty() {
        return Ok(None);
    }
    let tolerance = match plan.tolerance {
        Some(t) => t,
        None => {
            let b = shape.measure()?;
            DEFAULT_RELATIVE_TOLERANCE * (0..3).map(|k| (b.bbox_max[k] - b.bbox_min[k]).powi(2)).sum::<f64>().sqrt()
        }
    };
    let cuts: Vec<PlanCut> = planes
        .iter()
        .map(|(a, c)| PlanCut {
            point: a.map(|x| x * c),
            normal: *a,
            role: CutRole::Frame,
            feature: None,
            reason: "symmetry: the body is its own mirror image across this plane".into(),
        })
        .collect();
    let refs: Vec<&PlanCut> = cuts.iter().collect();
    let (cleared, made, skipped) = cut(s, name, &plan.body, &refs, tolerance)?;
    let undo = |s: &mut Session| -> Result<()> { cut(s, name, &plan.body, &[], tolerance).map(|_| ()) };
    if !skipped.is_empty() {
        undo(s)?;
        return Ok(None);
    }
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    let cells: Vec<(String, V3)> = design
        .bodies
        .iter()
        .filter(|b| design.source_of(&b.name) == Some(plan.body.as_str()))
        .filter_map(|b| b.measure().ok().map(|m| (b.name.clone(), std::array::from_fn(|k| 0.5 * (m.bbox_min[k] + m.bbox_max[k])))))
        .collect();
    let fundamental = cells.iter().find(|(_, mid)| planes.iter().all(|(a, c)| dot(*mid, *a) < *c)).map(|(n, _)| n.clone());
    let Some(cell) = fundamental else {
        undo(s)?;
        return Ok(None);
    };
    let mut cell_plan = build(design, &cell, Some(plan.size), plan.layers, MAX_ELEMENTS, &show_for(s))?;
    cell_plan.size = plan.size;
    cell_plan.tolerance = Some(tolerance);
    if cell_plan.mesher != PlanMesher::Project {
        undo(s)?;
        return Ok(None);
    }
    let copied = if corners { copy_corner(s, name, plan, &cuts, &cell, &cell_plan, tolerance)? } else { Err(None) };
    let (r, mut meshed, corner) = match copied {
        Ok(c) => (c.reply, c.cells, Some(c.summary)),
        Err(why) => {
            let corner = why.map(|w| json!({"skipped": w}));
            let cell_name = MeshPlan::default_name(&cell);
            s.model.set_mesh_plan(&cell_name, cell_plan.clone())?;
            (apply_projection(s, &cell_name, &cell_plan, tolerance)?, vec![cell.clone()], corner)
        }
    };
    // The cells now (the corner copy cuts the body further).
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    let cells: Vec<String> = design.bodies.iter().filter(|b| design.source_of(&b.name) == Some(plan.body.as_str())).map(|b| b.name.clone()).collect();
    let meshed_first = meshed.clone();
    for (a, c) in &planes {
        fieldcraft_model::progress::stage("Plan: mirroring the cell", None, "");
        fieldcraft_model::progress::check()?;
        (super::mirror::MIRROR.run)(s, &json!({"plane": {"point": a.map(|x| x * c), "normal": a}, "bodies": meshed}))?;
        meshed = cells.iter().filter(|n| !s.model.body_elements(n).is_empty()).cloned().collect();
    }
    if let Some(missing) = cells.iter().find(|n| s.model.body_elements(n).is_empty()) {
        return Err(EngineError::Params(format!("mesh plan `{name}`: the mirror images do not reach cell `{missing}`")));
    }
    let elements: usize = cells.iter().map(|n| s.model.body_elements(n).len()).sum();
    let mut reply = json!({
        "name": name,
        "body": plan.body,
        "mesher": "project",
        "symmetry": planes.iter().map(|(a, c)| json!({"normal": a, "offset": c})).collect::<Vec<_>>(),
        "layout": r["layout"],
        "blend_layers": r["blend_layers"],
        "cuts": made,
        "skipped_cuts": [],
        "cells": cells,
        "meshed": meshed_first,
        "corners": corner,
        "size": plan.size,
        "tolerance": tolerance,
        "elements": elements,
        "min_scaled_jacobian": r["min_scaled_jacobian"],
        "all_hex": true,
        // One cell's mesh and its mirror images, joined node to node on the planes.
        "conforming": true,
        "mapped": plan.features.iter().map(|f| json!({"id": f.id, "name": f.name, "type": f.kind, "template": f.template})).collect::<Vec<_>>(),
        "pending": [],
        "unmapped": plan.unmapped,
        "mesh": r["mesh"],
    });
    reply = cleared.report_if_any(s, reply);
    take_replaced(&mut reply, &r);
    if let Some(o) = reply.as_object_mut() {
        o.remove("results_cleared");
    }
    Ok(Some(reply))
}

/// The fundamental cell meshed by copying one end's corner onto the other end.
struct CornerCopy {
    /// The middle cell's mesh reply, with the smallest scaled Jacobian of both meshes.
    reply: Value,
    /// The cells meshed: the corner, the middle and the corner's image.
    cells: Vec<String>,
    summary: Value,
}

/// Where across the sheet past a corner the candidate cuts lie. Holding the cut's sheet planar
/// costs the corner some freedom to smooth, and how much depends on where it runs, so each
/// candidate is meshed and the best kept.
const CORNER_FRACTIONS: [f64; 4] = [0.5, 0.625, 0.75, 0.875];

/// A block's corners in parameter space (hex8 order).
const BLOCK_CORNERS: [[usize; 3]; 8] = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];

/// The corner planes of layout `t` along axis `a`, over a body from `lo` to `hi` along it: the
/// planes of the layout's vertices from the low end inward whose mirror images across the
/// middle `c` are planes too, up to the first that is not; the candidate cuts between the last
/// of them and that one ([`CORNER_FRACTIONS`] of the way), and `c`. `None` when there is none (or only the end), or the
/// planes mirror past the middle (the body's own symmetry meshes it then).
fn corner_planes(t: &Topology, a: V3, lo: f64, hi: f64, eps: f64) -> Option<(Vec<f64>, f64)> {
    let c = 0.5 * (lo + hi);
    let mut xs: Vec<f64> = t.vertices.iter().map(|v| dot(*v, a)).collect();
    xs.sort_by(f64::total_cmp);
    xs.dedup_by(|x, y| (*x - *y).abs() <= eps);
    let has = |x: f64| xs.iter().any(|y| (y - x).abs() <= eps);
    let k = xs.iter().take_while(|&&x| x < c - eps && has(2.0 * c - x)).count();
    if k < 2 {
        return None;
    }
    let (&p, &q) = (xs.get(k - 1)?, xs.get(k)?);
    if q >= c - eps {
        return None;
    }
    // Cuts across the sheet from midway to near its far plane; the image must cut no plane, and
    // the high end hold as many planes beyond it as the low end before the cut.
    let cuts: Vec<f64> = CORNER_FRACTIONS
        .iter()
        .map(|f| p + f * (q - p))
        .filter(|&x1| {
            let x2 = 2.0 * c - x1;
            xs.iter().filter(|&&x| x < x1).count() == xs.iter().filter(|&&x| x > x2).count() && !has(x2)
        })
        .collect();
    (!cuts.is_empty()).then_some((cuts, c))
}

/// The blocks of `t` that `keep` takes, over their own vertices.
fn sub_topology(t: &Topology, keep: impl Fn(&[u32; 8]) -> bool) -> Topology {
    let mut index: std::collections::BTreeMap<u32, u32> = std::collections::BTreeMap::new();
    let mut out = Topology::default();
    for block in t.blocks.iter().filter(|b| keep(b)) {
        let b = block.map(|v| {
            *index.entry(v).or_insert_with(|| {
                out.vertices.push(t.vertices.get(v as usize).copied().unwrap_or_default());
                out.vertices.len() as u32 - 1
            })
        });
        out.blocks.push(b);
    }
    out
}

/// Split the sheet of blocks of `t` the plane `p·a = x` runs through, there.
fn split_at(t: &mut Topology, a: V3, x: f64) -> Result<()> {
    let span = |t: &Topology, b: usize| -> Option<(f64, f64)> {
        let xs = t.corners(b).ok()?.map(|p| dot(p, a));
        Some((xs.iter().copied().fold(f64::INFINITY, f64::min), xs.iter().copied().fold(f64::NEG_INFINITY, f64::max)))
    };
    let Some(b) = (0..t.blocks.len()).find(|&b| span(t, b).is_some_and(|(lo, hi)| lo < x && x < hi)) else {
        return Err(EngineError::Params(format!("no block crosses the corner cut at {x}")));
    };
    let (d, _) = t.direction_along(b, a)?;
    let c = t.corners(b)?;
    // Where the block's side-0 and side-1 faces are along `a`, on average.
    let side = |s: usize| {
        let at: Vec<f64> = BLOCK_CORNERS.iter().zip(&c).filter(|(bits, _)| bits[d] == s).map(|(_, p)| dot(*p, a)).collect();
        at.iter().sum::<f64>() / at.len().max(1) as f64
    };
    let f = (x - side(0)) / (side(1) - side(0));
    if !(f > 0.0 && f < 1.0) {
        return Err(EngineError::Params(format!("the corner cut at {x} does not cross block {b}")));
    }
    t.split(b, d, &[f])?;
    Ok(())
}

/// Mesh the fundamental cell `cell` of a symmetric plan by copying one corner onto the other.
/// Along an in-plane axis the body is not symmetric about, when the planes of the cell's layout
/// mirror one another from both ends inward up to a pair (an end with the same features, such as
/// an enclosure's corners with their bosses), and the cell's surface beyond them is mirror
/// symmetric too, the cell is cut midway into the next sheet of blocks and at its image. The low
/// corner is meshed with its blocks of the layout and mirrored onto the high one; the middle is
/// meshed with its blocks, both cut faces fixed (taken from the corner meshes as they are). Both
/// corners then mesh exactly alike. `None`, with nothing changed, when there is no such corner
/// or its meshes fail.
fn copy_corner(
    s: &mut Session,
    name: &str,
    plan: &MeshPlan,
    sym: &[PlanCut],
    cell: &str,
    cell_plan: &MeshPlan,
    tolerance: f64,
) -> Result<std::result::Result<CornerCopy, Option<String>>> {
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    let (Some(shape), Some(d)) = (design.body(cell), cell_plan.direction) else { return Ok(Err(None)) };
    let layout = match layout_of(design, cell_plan, &show_for(s))? {
        Some(Ok(l)) => l,
        _ => return Ok(Err(None)),
    };
    let m = shape.measure()?;
    let (lo, hi) = (m.bbox_min, m.bbox_max);
    let diagonal = (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt();
    let eps = CORNER_EPS * diagonal.max(1.0);
    let Some(mirror) = Mirror::of(shape)? else { return Ok(Err(None)) };
    let mut found = None;
    for a in in_plane_axes(d).into_iter().filter(|a| AXES.iter().any(|x| parallel(*x, *a))) {
        if sym.iter().any(|c| parallel(c.normal, a)) {
            continue;
        }
        let Some((cuts, c)) = corner_planes(&layout.topology, a, dot(lo, a), dot(hi, a), eps) else { continue };
        let clear = |x: f64| cell_plan.rings.iter().all(|r| (dot(r.center, a) - x).abs() >= 1.01 * r.radius);
        let cuts: Vec<f64> = cuts.into_iter().filter(|&x| clear(x) && clear(2.0 * c - x)).collect();
        // The surface beyond the deepest cut mirrors, so beyond every cut it does.
        let Some(deepest) = cuts.iter().copied().reduce(f64::max) else { continue };
        if mirror.symmetric(a, c, |p| dot(p, a) < deepest || dot(p, a) > 2.0 * c - deepest) {
            found = Some((a, cuts, c));
            break;
        }
    }
    let Some((a, at, c)) = found else { return Ok(Err(None)) };
    let before = s.snapshot();
    match corner_meshes(s, name, plan, sym, cell, &layout.topology, a, &at, c, tolerance) {
        Ok(mut copy) => {
            // The layout is the cell's, split at the cuts.
            if let Some(o) = copy.reply.as_object_mut() {
                o.insert("layout".into(), json!(layout.summary(&show_for(s))));
                o.insert(
                    "blend_layers".into(),
                    json!(layout.layers.iter().map(|(side, depth, what)| json!({"side": side, "depth": depth, "for": what})).collect::<Vec<_>>()),
                );
            }
            Ok(Ok(copy))
        }
        Err(e) if e.is_cancelled() => Err(e),
        // A corner copy that does not mesh leaves the cell to be meshed whole.
        Err(e) => {
            s.restore(before);
            let shown = s.units().clone();
            Ok(Err(Some(e.shown(crate::units::geometry(), &shown).to_string())))
        }
    }
}

/// Carry out [`copy_corner`] along axis `a`: cuts at `at` (low, high), mirror plane `p·a = c`.
#[allow(clippy::too_many_arguments)]
fn corner_meshes(
    s: &mut Session,
    name: &str,
    plan: &MeshPlan,
    sym: &[PlanCut],
    cell: &str,
    t0: &Topology,
    a: V3,
    candidates: &[f64],
    c: f64,
    tolerance: f64,
) -> Result<CornerCopy> {
    let field = super::sizes::sizing(s)?.local_field(plan.size)?;
    let opts = fieldcraft_mesh::project::ProjectOptions { size: plan.size, field, ..Default::default() };
    let whole_body = super::project::open_body(s, cell, Some(tolerance), &json!({}))?;
    let merge = super::project::merge_distance(None, whole_body.diagonal())?;
    // What the cell meshed whole gets (no corner may be worse), and per candidate cut the whole
    // cell meshed with its layout split there and the sheets at the cut and its image held
    // planar, the corner taken from that mesh: independent meshes, made side by side. The best
    // corner wins.
    let whole = &whole_body;
    // The job's progress on the pool's threads too, so a cancel stops the candidates.
    let progress = fieldcraft_model::progress::current();
    let runs = fieldcraft_parallel::map(candidates.len() + 1, |k| {
        fieldcraft_model::progress::scoped(progress.clone(), || match k.checked_sub(1).and_then(|i| candidates.get(i)) {
            None => super::project::project_only(whole, cell, t0, &opts).map(|m| (m.report.min_scaled_jacobian, None)),
            Some(&x1) => corner_of_whole(whole, cell, t0, a, [x1, 2.0 * c - x1], &opts).map(|r| (r.sj, Some(([x1, 2.0 * c - x1], r)))),
        })
    });
    let mut reference: Option<f64> = None;
    let mut best: Option<(f64, [f64; 2], CornerOfWhole)> = None;
    for run in runs {
        match run {
            Ok((sj, None)) => reference = Some(sj),
            Ok((sj, Some((at, corner)))) => {
                if best.as_ref().is_none_or(|b| sj > b.0) {
                    best = Some((sj, at, corner));
                }
            }
            Err(e) if e.is_cancelled() => return Err(e),
            Err(_) => {}
        }
    }
    let reference = reference.ok_or_else(|| EngineError::Params("the cell does not mesh whole".into()))?;
    let Some((corner_sj, at, corner)) = best else {
        return Err(EngineError::Params("no corner cut meshes".into()));
    };
    let [x1, x2] = at;
    if corner_sj < reference - SAME_QUALITY {
        return Err(EngineError::Params(format!(
            "the best corner has a smallest scaled Jacobian of {corner_sj:.3}, below the {reference:.3} of the cell meshed whole"
        )));
    }
    let mut cuts: Vec<PlanCut> = sym.to_vec();
    for x in at {
        cuts.push(PlanCut {
            point: a.map(|k| k * x),
            normal: a,
            role: CutRole::Frame,
            feature: None,
            reason: "corner: the body's two ends mesh alike, the low one meshed and mirrored onto the high one".into(),
        });
    }
    let refs: Vec<&PlanCut> = cuts.iter().collect();
    let (_, _, skipped) = cut(s, name, &plan.body, &refs, tolerance)?;
    if !skipped.is_empty() {
        return Err(EngineError::Params("a corner cut was left out".into()));
    }
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    // The fundamental cells (on the negative side of every symmetry plane), by where they lie
    // along the axis.
    let fundamental: Vec<(String, f64)> = design
        .bodies
        .iter()
        .filter(|b| design.source_of(&b.name) == Some(plan.body.as_str()))
        .filter_map(|b| b.measure().ok().map(|m| (b.name.clone(), std::array::from_fn::<f64, 3, _>(|k| 0.5 * (m.bbox_min[k] + m.bbox_max[k])))))
        .filter(|(_, mid)| sym.iter().all(|cut| dot(*mid, cut.normal) < dot(cut.point, cut.normal)))
        .map(|(n, mid)| (n, dot(mid, a)))
        .collect();
    let pick = |test: &dyn Fn(f64) -> bool| -> Result<String> {
        let mut found = fundamental.iter().filter(|(_, x)| test(*x)).map(|(n, _)| n.clone());
        match (found.next(), found.next()) {
            (Some(n), None) => Ok(n),
            _ => Err(EngineError::Params("the corner cuts do not leave one cell at each end and one between".into())),
        }
    };
    let (low, middle, high) = (pick(&|x| x < x1)?, pick(&|x| x > x1 && x < x2)?, pick(&|x| x > x2)?);
    let centre = |b: &[u32; 8]| b.iter().map(|&v| corner.t.vertices.get(v as usize).map_or(0.0, |p| dot(*p, a))).sum::<f64>() / 8.0;
    let middle_blocks = sub_topology(&corner.t, |b| (x1..x2).contains(&centre(b)));
    let body = super::project::open_body(s, &low, Some(tolerance), &json!({}))?;
    super::project::add_hexes(s, &body, &corner.nodes, &corner.hexes, &corner.tags, merge)?;
    (super::mirror::MIRROR.run)(s, &json!({"plane": {"point": a.map(|k| k * c), "normal": a}, "bodies": [low]}))?;
    let body = super::project::open_body(s, &middle, Some(tolerance), &json!({}))?;
    let (fixed, fixed_faces) = cut_face_patches(s, &body, &middle, plan.size)?;
    if fixed_faces.len() < 2 {
        return Err(EngineError::Params("the middle cell does not meet both corners' meshes".into()));
    }
    let opts = fieldcraft_mesh::project::ProjectOptions { fixed, ..opts };
    let mut mesh = super::project::project_body(s, &middle, body, &middle_blocks, &opts, merge)?;
    let middle_sj = mesh["min_scaled_jacobian"].as_f64().unwrap_or(-1.0);
    if middle_sj < reference - QUALITY_SLACK {
        return Err(EngineError::Params(format!(
            "the middle meshed round the corners has a smallest scaled Jacobian of {middle_sj:.3}, below the {reference:.3} of the cell meshed whole"
        )));
    }
    let corner_sj = Some(corner_sj);
    let summary = json!({
        "axis": a,
        "cuts": at,
        "meshed": low,
        "image": high,
        "middle": middle,
        "fixed": fixed_faces,
        "min_scaled_jacobian": {"whole_cell": reference, "corner": corner_sj, "middle": mesh["min_scaled_jacobian"]},
    });
    let min = [corner_sj, mesh["min_scaled_jacobian"].as_f64()].into_iter().flatten().reduce(f64::min);
    if let Some(o) = mesh.as_object_mut() {
        o.remove("topology");
    }
    let reply = json!({"min_scaled_jacobian": min, "mesh": mesh});
    Ok(CornerCopy { reply, cells: vec![low, middle, high], summary })
}

/// The low corner (below `at[0]`) of the cell meshed whole with layout `t0` split at both cuts
/// `at`, the sheets there held planar: its smallest scaled Jacobian, the split layout, and its
/// hexes over their own nodes (those on the cut sheet put exactly on the plane).
#[allow(clippy::type_complexity)]
fn corner_of_whole(body: &super::project::Body, cell: &str, t0: &Topology, a: V3, at: [f64; 2], opts: &ProjectOptions) -> Result<CornerOfWhole> {
    let x1 = at[0];
    let mut t = t0.clone();
    for x in at {
        split_at(&mut t, a, x)?;
    }
    let sheets = ProjectOptions { planes: at.iter().map(|&x| (a.map(|k| k * x), a)).collect(), ..opts.clone() };
    let whole = super::project::project_only(body, cell, &t, &sheets)?;
    let along = |p: &[f64; 3]| dot(*p, a);
    let hex_centre = |h: &[u32; 8]| h.iter().map(|&v| whole.nodes.get(v as usize).map_or(0.0, along)).sum::<f64>() / 8.0;
    let beyond: std::collections::BTreeSet<u32> = whole.hexes.iter().filter(|h| hex_centre(h) >= x1).flat_map(|h| h.iter().copied()).collect();
    let mut index: std::collections::BTreeMap<u32, u32> = std::collections::BTreeMap::new();
    let mut nodes: Vec<[f64; 3]> = Vec::new();
    let mut hexes: Vec<[u32; 8]> = Vec::new();
    let mut element_of: std::collections::BTreeMap<u32, u32> = std::collections::BTreeMap::new();
    for (e, h) in whole.hexes.iter().enumerate().filter(|(_, h)| hex_centre(h) < x1) {
        element_of.insert(e as u32, hexes.len() as u32);
        hexes.push(h.map(|v| {
            *index.entry(v).or_insert_with(|| {
                let mut p = whole.nodes.get(v as usize).copied().unwrap_or_default();
                if beyond.contains(&v) {
                    let off = x1 - along(&p);
                    p = std::array::from_fn(|k| p[k] + off * a[k]);
                }
                nodes.push(p);
                nodes.len() as u32 - 1
            })
        }));
    }
    if hexes.is_empty() {
        return Err(EngineError::Params("the corner holds no elements".into()));
    }
    let sj = hexes
        .iter()
        .map(|h| {
            let x = h.map(|v| nodes.get(v as usize).copied().unwrap_or_default());
            fieldcraft_mesh::quality::element_quality(fieldcraft_model::ElementKind::Hex8, &x)
                .ok()
                .and_then(|q| q.get(fieldcraft_mesh::quality::Metric::ScaledJacobian))
                .unwrap_or(-1.0)
        })
        .fold(f64::INFINITY, f64::min);
    let names = body.face_names();
    let tags =
        whole.boundary.iter().filter_map(|f| Some(((*element_of.get(&f.element)?, f.face), names.get(f.brep_face as usize)?.clone()))).collect();
    Ok(CornerOfWhole { sj, t, nodes, hexes, tags })
}

/// The low corner [`corner_of_whole`] takes from the whole cell's mesh.
struct CornerOfWhole {
    sj: f64,
    /// The layout split at the cuts.
    t: Topology,
    nodes: Vec<[f64; 3]>,
    hexes: Vec<[u32; 8]>,
    /// The B-rep face (by name) the whole cell's mesh tagged each outer element face with.
    tags: std::collections::BTreeMap<(u32, u8), String>,
}

/// The fixed patches of body `name` (opened as `b`): each of its cut faces another cell's mesh
/// already covers, with the faces' names.
fn cut_face_patches(
    s: &Session,
    b: &super::project::Body,
    name: &str,
    size: f64,
) -> Result<(Vec<fieldcraft_mesh::project::FixedPatch>, Vec<String>)> {
    let mut fixed = Vec::new();
    let mut faces: Vec<String> = Vec::new();
    let cut_faces: Vec<String> = b.face_names().iter().filter(|n| n.starts_with("partition:")).cloned().collect();
    for face in cut_faces {
        let patch = super::project::fixed_patch(s, b, name, std::slice::from_ref(&face), size)?;
        if !patch.quads.is_empty() {
            fixed.push(patch);
            faces.push(face);
        }
    }
    Ok((fixed, faces))
}

/// Carry out a plan whose mesher is projection blocking: its earlier partitions removed, one block
/// layout over the whole body (rings, blend layer) projected and smoothed, the body's mesh
/// replaced. No fallback: a layout that does not fit, or a mesh below the projection's quality
/// floor, is the error.
fn apply_projection(s: &mut Session, name: &str, plan: &MeshPlan, tolerance: f64) -> Result<Value> {
    let (cleared, _, _) = cut(s, name, &plan.body, &[], tolerance)?;
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    let layout = match layout_of(design, plan, &show_for(s))? {
        Some(Ok(l)) => l,
        Some(Err(why)) => return Err(EngineError::Params(format!("mesh plan `{name}`: no projection layout for body `{}`: {why}", plan.body))),
        None => return Err(EngineError::Params(format!("mesh plan `{name}`: projection blocking needs rings or blends; set the mesher to sweep"))),
    };
    // Partitions of the user's that cut the body leave cells, each meshed on its own.
    if design.body(&plan.body).is_none() {
        let cells: Vec<&str> =
            design.bodies.iter().filter(|b| design.source_of(&b.name) == Some(plan.body.as_str())).map(|b| b.name.as_str()).collect();
        return Err(EngineError::Params(format!(
            "mesh plan `{name}`: body `{}` is cut into the cells {cells:?} by partitions; plan and apply a cell (mesh.plan {{\"body\": \"{}\"}}), or delete the partitions",
            plan.body,
            cells.first().copied().unwrap_or_default()
        )));
    }
    let body = super::project::open_body(s, &plan.body, Some(tolerance), &json!({}))?;
    let merge = super::project::merge_distance(None, body.diagonal())?;
    // Cut faces another cell's mesh already covers are taken as given, so the cells meet node
    // to node there.
    let (fixed, fixed_faces) = cut_face_patches(s, &body, &plan.body, plan.size)?;
    let field = super::sizes::sizing(s)?.local_field(plan.size)?;
    let opts = fieldcraft_mesh::project::ProjectOptions { size: plan.size, fixed, field, ..Default::default() };
    let mut mesh = super::project::project_body(s, &plan.body, body, &layout.topology, &opts, merge)
        .map_err(|e| EngineError::context(format!("mesh plan `{name}`: the projected layout does not mesh"), e))?;
    if let Some(m) = mesh.as_object_mut() {
        // The block listing is long; the plan reports the layout instead.
        m.remove("topology");
    }
    let mut reply = json!({
        "name": name,
        "body": plan.body,
        "mesher": "project",
        "layout": layout.summary(&show_for(s)),
        "blend_layers": layout.layers.iter().map(|(side, depth, what)| json!({"side": side, "depth": depth, "for": what})).collect::<Vec<_>>(),
        "cuts": [],
        "skipped_cuts": [],
        "cells": [plan.body],
        "size": plan.size,
        "tolerance": tolerance,
        "fixed": fixed_faces,
        "blocks": mesh["blocks"],
        "elements": mesh["elements"],
        "min_scaled_jacobian": mesh["min_scaled_jacobian"],
        "all_hex": true,
        // One block topology, one mesh: every block face shared node to node.
        "conforming": true,
        "mapped": plan.features.iter().map(|f| json!({"id": f.id, "name": f.name, "type": f.kind, "template": f.template})).collect::<Vec<_>>(),
        "pending": [],
        "unmapped": plan.unmapped,
        "mesh": mesh,
    });
    reply = cleared.report_if_any(s, reply);
    if let Some(m) = reply.get_mut("mesh").and_then(Value::as_object_mut) {
        let theirs = json!({"replaced": m.remove("replaced"), "stale_sets": m.remove("stale_sets")});
        m.remove("results_cleared");
        take_replaced(&mut reply, &theirs);
    }
    if let Some(o) = reply.as_object_mut() {
        o.remove("results_cleared");
    }
    Ok(reply)
}

/// Replace plan `name`'s partitions of `body` with `cuts` (named `<name>-cut<k>` after their
/// place in the plan). A cut the kernel cannot make (or that divides nothing left), or that
/// leaves a cell whose boundary no mesher takes at `tolerance` (a plane the tessellation does
/// not close along), is left out: the latest such cut first, one at a time. Returns the meshes
/// the first recut removed, the cuts made and those left out with why.
#[allow(clippy::type_complexity)]
fn cut(
    s: &mut Session,
    name: &str,
    body: &str,
    cuts: &[&PlanCut],
    tolerance: f64,
) -> Result<(super::replace::Cleared, Vec<String>, Vec<(String, String)>)> {
    let display = s.units().clone();
    for old in plan_cut_names(s, name) {
        s.model.remove_partition(&old)?;
        if s.model.datum_users(&old).is_empty() && s.model.datum(&old).is_some() {
            s.model.remove_datum(&old)?;
        }
    }
    let mut names = Vec::with_capacity(cuts.len());
    for (k, c) in cuts.iter().enumerate() {
        let n = format!("{name}-cut{}", k + 1);
        if s.model.datum(&n).is_some() {
            return Err(EngineError::Params(format!("datum `{n}` exists and is not the plan's: rename it or the plan")));
        }
        s.model.add_datum(&n, Datum::Plane { point: c.point, normal: c.normal })?;
        s.model.add_partition(Partition { name: n.clone(), datum: n.clone(), bodies: vec![body.to_string()], suppressed: false })?;
        names.push(n);
    }
    let cleared = super::datum::recut(s)?;
    let mut skipped: Vec<(String, String)> = Vec::new();
    for _ in 0..=names.len() {
        fieldcraft_model::progress::stage("Plan: checking the cells", None, "");
        fieldcraft_model::progress::check()?;
        let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
        let live: Vec<&String> = names.iter().filter(|n| !skipped.iter().any(|(k, _)| k == *n)).collect();
        let mut failed: Vec<(String, String)> = live
            .iter()
            .filter_map(|n| {
                let feature = format!("partition {n}");
                let messages: Vec<&str> = design.errors.iter().filter(|e| e.feature == feature).map(|e| e.message.as_str()).collect();
                (!messages.is_empty()).then(|| ((*n).clone(), messages.join("; ")))
            })
            .collect();
        if failed.is_empty() {
            let bad = design.bodies.iter().filter(|b| design.source_of(&b.name) == Some(body)).find_map(|cell| {
                let why = super::datum::unmeshable(cell, Some(tolerance), &display).err()?;
                // The latest of the plan's cuts that bound the cell.
                let latest = live.iter().rev().find(|n| cell.faces().iter().any(|f| f.name.starts_with(&format!("partition:{n}:"))))?;
                Some(((*latest).clone(), format!("cell `{}` it bounds has a boundary no mesher takes ({why})", cell.name)))
            });
            failed.extend(bad);
        }
        if failed.is_empty() {
            break;
        }
        for (n, _) in &failed {
            s.model.remove_partition(n)?;
            s.model.remove_datum(n)?;
        }
        skipped.extend(failed);
        // The meshes this recut removes are those of cells the first one made, which had none.
        super::datum::recut(s)?;
    }
    let made = names.into_iter().filter(|n| !skipped.iter().any(|(f, _)| f == n)).collect();
    Ok((cleared, made, skipped))
}

/// What a mesher's reply says of the mesh a plan made.
struct Summary {
    min_scaled_jacobian: Option<f64>,
    elements: usize,
    /// Every element a hexahedron.
    all_hex: bool,
    /// The cells joined node to node into one mesh: every interface merged (none tied), and the
    /// merged interfaces connect every cell.
    conforming: bool,
}

/// The [`Summary`] of a mesher's reply (`mesh.hex.generate`'s `bodies`, or `mesh.cells`'
/// `cells` and `interfaces`) over `cells`.
fn summary(mesh: &Value, cells: &[String]) -> Summary {
    let mut min: Option<f64> = None;
    let (mut elements, mut all_hex) = (0usize, true);
    let rows = mesh.get("bodies").or_else(|| mesh.get("cells")).and_then(Value::as_array).into_iter().flatten();
    for b in rows {
        if let Some(v) = b.get("quality").and_then(|q| q.get("min_scaled_jacobian")).and_then(Value::as_f64) {
            min = Some(min.map_or(v, |m: f64| m.min(v)));
        }
        elements += b.get("elements").and_then(Value::as_u64).map_or(0, |n| usize::try_from(n).unwrap_or(usize::MAX));
        all_hex &= b.get("kind").and_then(Value::as_str).is_some_and(|k| k.starts_with("hex"));
    }
    // The cells the merged interfaces join, grown from the first.
    let interfaces: Vec<&Value> = mesh.get("interfaces").and_then(Value::as_array).into_iter().flatten().collect();
    let merged = |i: &Value| i.get("joined").and_then(Value::as_str) == Some("merged");
    let mut joined: Vec<&str> = cells.first().map(String::as_str).into_iter().collect();
    loop {
        let before = joined.len();
        for i in interfaces.iter().copied().filter(|i| merged(i)) {
            let pair: Vec<&str> = i.get("cells").and_then(Value::as_array).into_iter().flatten().filter_map(Value::as_str).collect();
            if let [a, b] = pair.as_slice() {
                match (joined.contains(a), joined.contains(b)) {
                    (true, false) => joined.push(b),
                    (false, true) => joined.push(a),
                    _ => {}
                }
            }
        }
        if joined.len() == before {
            break;
        }
    }
    let conforming = interfaces.iter().all(|i| merged(i)) && cells.iter().all(|c| joined.contains(&c.as_str()));
    Summary { min_scaled_jacobian: min, elements, all_hex, conforming }
}

#[cfg(test)]
mod tests {
    use serde_json::{Value, json};

    use crate::Session;

    fn open_script(script: &str) -> (crate::tests::TempDir, Session) {
        let bytes = fieldcraft_geometry::script_design(script).unwrap();
        let (dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        (dir, s)
    }

    fn open(commands: Value) -> (crate::tests::TempDir, Session) {
        open_script(&json!({"commands": commands}).to_string())
    }

    fn templates(plan: &Value) -> Vec<(String, String)> {
        plan["features"].as_array().unwrap().iter().map(|f| (f["type"].as_str().unwrap().into(), f["template"].as_str().unwrap().into())).collect()
    }

    fn pairs(list: &[(&str, &str)]) -> Vec<(String, String)> {
        list.iter().map(|(a, b)| (a.to_string(), b.to_string())).collect()
    }

    /// Plan, apply (the whole body, not by its symmetry), and check the mesh is all-hex,
    /// conforming and fills the body.
    fn apply(s: &mut Session, plan: &Value, min_sj: f64) -> Value {
        let r = s.execute("mesh.plan.apply", &json!({"plan": plan["name"], "symmetry": false})).unwrap();
        assert_eq!((r["all_hex"].as_bool(), r["conforming"].as_bool()), (Some(true), Some(true)), "{r}");
        let sj = r["min_scaled_jacobian"].as_f64().unwrap();
        assert!(sj >= min_sj, "min scaled Jacobian {sj} < {min_sj}: {r}");
        r
    }

    /// The mesh's volume against the bodies'.
    fn volume_matches(s: &mut Session) {
        let q = s.execute("mesh.quality", &json!({})).unwrap();
        let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
        let geometry: f64 = bodies.as_array().unwrap().iter().map(|b| b["volume"].as_f64().unwrap()).sum();
        let mesh = q["volume"].as_f64().unwrap_or_else(|| q["metrics"]["volume"]["sum"].as_f64().unwrap_or_default());
        if mesh > 0.0 {
            assert!((mesh - geometry).abs() < 0.02 * geometry, "mesh volume {mesh} vs {geometry}");
        }
    }

    #[test]
    fn a_plate_with_holes_is_one_block_layout_with_an_o_grid_round_each_hole() {
        let (_d, mut s) = open(json!([
            {"command": "solid.box", "params": {"length": 60, "width": 40, "height": 6, "body_name": "Plate"}},
            {"command": "solid.hole", "params": {"position": [15, 20, 6], "diameter": 8, "depth": 6}},
            {"command": "solid.hole", "params": {"position": [45, 20, 6], "diameter": 8, "depth": 6}},
        ]));
        let plan = s.execute("mesh.plan", &json!({"body": "Plate"})).unwrap();
        assert_eq!(templates(&plan), pairs(&[("box", "block"), ("hole", "o_grid"), ("hole", "o_grid")]));
        assert_eq!(plan["direction"], json!([0.0, 0.0, 1.0]));
        // Two O-grid rings, each projecting onto its hole's walls, framed short of the other hole.
        let rings = plan["rings"].as_array().unwrap();
        assert_eq!(rings.len(), 2, "{plan}");
        for r in rings {
            assert_eq!((r["kind"].as_str(), r["radius"].as_f64()), (Some("o_grid"), Some(4.0)));
            assert!(!r["faces"].as_array().unwrap().is_empty(), "{r}");
            assert!(r["frame"].as_f64().unwrap() > 4.0);
        }
        // One block layout over the whole plate: no partition cuts.
        assert_eq!((plan["mesher"].as_str(), plan.get("cuts")), (Some("project"), None), "{plan}");
        let r = apply(&mut s, &plan, 0.5);
        assert_eq!(r["cells"], json!(["Plate"]), "{r}");
        assert_eq!(r["pending"], json!([]), "{r}");
        // The ring nodes lie on the hole walls: the holes are round.
        assert!(r["mesh"]["max_gap"].as_f64().unwrap() < 1e-6, "{r}");
        volume_matches(&mut s);
        // Saved: applied again by name it replaces its mesh.
        let again = s.execute("mesh.plan.apply", &json!({"plan": "Plan-Plate", "symmetry": false})).unwrap();
        assert_eq!(again["elements"], r["elements"]);
        assert!(s.model.partitions().is_empty());
        assert!(!again["replaced"].as_array().unwrap().is_empty(), "{again}");
    }

    #[test]
    fn a_sketched_plate_sweeps_its_profile_and_o_grids_its_circle() {
        let bytes = fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap();
        let (_d, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        let body = s.execute("geometry.bodies", &json!({})).unwrap()[0]["name"].as_str().unwrap().to_string();
        let plan = s.execute("mesh.plan", &json!({"body": body})).unwrap();
        assert_eq!(templates(&plan), pairs(&[("extrude", "sweep"), ("extrude", "imprint")]));
        assert_eq!(plan["rings"].as_array().unwrap().len(), 1, "{plan}");
        assert_eq!(plan["rings"][0]["kind"], "o_grid");
        apply(&mut s, &plan, 0.5);
        volume_matches(&mut s);
    }

    #[test]
    fn a_cylinder_is_a_butterfly_swept_along_its_axis() {
        let (_d, mut s) = open(json!([{"command": "solid.cylinder", "params": {"diameter": 20, "height": 30, "body_name": "Cyl"}}]));
        let plan = s.execute("mesh.plan", &json!({"body": "Cyl"})).unwrap();
        assert_eq!(templates(&plan), pairs(&[("cylinder", "butterfly")]));
        assert_eq!((plan["mesher"].as_str(), plan.get("cuts")), (Some("project"), None));
        assert_eq!(plan["rings"][0]["kind"], "butterfly");
        assert_eq!(plan["rings"][0]["frame"], 5.0);
        apply(&mut s, &plan, 0.7);
        volume_matches(&mut s);
    }

    #[test]
    fn a_joined_boss_is_a_butterfly_on_the_block() {
        let (_d, mut s) = open(json!([
            {"command": "solid.box", "params": {"length": 40, "width": 40, "height": 10, "body_name": "Base"}},
            {"command": "solid.cylinder", "params": {"diameter": 12, "height": 8, "base": [20, 20, 10], "operation": "join", "targets": ["Base"]}},
        ]));
        let plan = s.execute("mesh.plan", &json!({"body": "Base"})).unwrap();
        assert_eq!(templates(&plan), pairs(&[("box", "block"), ("cylinder", "butterfly")]));
        assert_eq!(plan["rings"][0]["kind"], "butterfly");
        apply(&mut s, &plan, 0.5);
        volume_matches(&mut s);
    }

    #[test]
    fn a_shelled_box_takes_layers_through_its_wall() {
        let (_d, mut s) = open(json!([
            {"command": "solid.box", "params": {"length": 40, "width": 30, "height": 20, "body_name": "Case"}},
            {"command": "solid.shell", "params": {"faces": [[20, 15, 20]], "thickness": 2, "body": "Case"}},
        ]));
        let plan = s.execute("mesh.plan", &json!({"body": "Case", "layers": 2})).unwrap();
        assert_eq!(templates(&plan), pairs(&[("box", "block"), ("shell", "layers")]));
        // Two elements through the 2 mm wall.
        assert_eq!(plan["size"], 1.0);
        apply(&mut s, &plan, 0.9);
        volume_matches(&mut s);
    }

    #[test]
    fn a_patterned_and_mirrored_part_meshes_every_copy() {
        let (_d, mut s) = open(json!([
            {"command": "solid.box", "params": {"length": 10, "width": 10, "height": 4, "body_name": "Tab"}},
            {"command": "solid.hole", "params": {"position": [5, 5, 4], "diameter": 4, "depth": 4}},
            {"command": "solid.pattern.rectangular", "params": {"bodies": ["Tab"], "dir1": [1, 0, 0], "count1": 3, "spacing1": 15}},
            {"command": "solid.mirror", "params": {"bodies": ["Tab"], "plane": {"origin": [0, 0, 0], "normal": [0, 1, 0]}}},
        ]));
        let plan = s.execute("mesh.plan", &json!({"body": "Tab.1"})).unwrap();
        assert_eq!(templates(&plan), pairs(&[("box", "block"), ("hole", "o_grid"), ("pattern", "copy"), ("mirror", "copy")]));
        let copies = plan["copies"].as_array().unwrap();
        assert_eq!(copies.len(), 3, "{plan}");
        let r = apply(&mut s, &plan, 0.5);
        assert_eq!(r["copies"].as_array().unwrap().len(), 3, "{r}");
        for c in r["copies"].as_array().unwrap() {
            assert_eq!(c["conforming"], true, "{c}");
            assert_eq!(c["elements"], r["elements"], "a copy meshes as the instance: {c}");
        }
        // Every body meshed.
        for b in s.execute("geometry.bodies", &json!({})).unwrap().as_array().unwrap() {
            assert!(!s.model.body_elements(b["name"].as_str().unwrap()).is_empty(), "{b}");
        }
    }

    /// A cut's rectangular opening whose ends are faces of the body (a pocket's end walls): the
    /// kernel cannot cut along a face (#427), so those frame sides move one element (twice the
    /// gap) out of the opening; the sides off every face stay where they are (#455).
    #[test]
    fn a_frame_side_on_an_opening_s_end_face_moves_one_element_out_of_it() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        let z = [0.0, 0.0, 1.0];
        let body = super::Clearance {
            corners: vec![[-20.0, -20.0, 0.0], [190.0, 20.0, 16.0]],
            // The pocket's end walls (x = 48 and x = 130), facing into it.
            planes: vec![(x, 48.0, "F6:side:l4"), (x.map(|c| -c), -130.0, "F6:side:l2")],
            curved: Vec::new(),
            gap: 1.5,
        };
        let pocket = super::Imprint {
            feature: 6,
            cut: true,
            extent: None,
            sides: vec![(x, 48.0), (x, 130.0), (y, -7.0), (y, 7.0)],
            center: Some([89.0, 0.0, 4.0]),
        };
        let mut cuts = Vec::new();
        super::imprint_cuts(&[pocket], z, &body, &mut cuts, &|q, v| crate::units::shown(crate::units::geometry(), crate::units::geometry(), q, v));
        let at = |n: [f64; 3]| -> Vec<f64> {
            let mut v: Vec<f64> = cuts.iter().filter(|c| super::parallel(c.normal, n)).map(|c| super::dot(c.point, n)).collect();
            v.sort_by(f64::total_cmp);
            v
        };
        assert_eq!(at(x), vec![45.0, 133.0], "{cuts:?}");
        assert_eq!(at(y), vec![-7.0, 7.0], "{cuts:?}");
    }

    #[test]
    fn a_fillet_is_a_blend_layer_and_hostile_plans_are_refused() {
        let (_d, mut s) = open(json!([
            {"command": "solid.box", "params": {"length": 20, "width": 20, "height": 20, "body_name": "B"}},
            {"command": "solid.fillet", "params": {"edges": [[20, 10, 20]], "radius": 2}},
        ]));
        let plan = s.execute("mesh.plan", &json!({"body": "B"})).unwrap();
        assert_eq!(templates(&plan), pairs(&[("box", "block"), ("fillet", "blend")]));
        assert_eq!(plan["mesher"], "project", "{plan}");
        let r = apply(&mut s, &plan, 0.5);
        // The layer along the filleted edge, on the two sides it runs between, at its radius.
        let layers = r["blend_layers"].as_array().unwrap();
        assert_eq!(layers.len(), 2, "{r}");
        assert!(layers.iter().all(|l| (l["depth"].as_f64().unwrap() - 2.0).abs() < 1e-6), "{r}");
        volume_matches(&mut s);
        let mut bad = plan.clone();
        bad["size"] = json!(-1.0);
        assert!(s.execute("mesh.plan.apply", &json!({"plan": bad})).is_err());
        let mut bad = plan.clone();
        bad["rings"] = json!([{"feature": 1, "kind": "o_grid", "center": [0, 0, 0], "axis": [0, 0, 0], "radius": 1, "frame": 2}]);
        assert!(s.execute("mesh.plan.apply", &json!({"plan": bad})).is_err());
        assert!(s.execute("mesh.plan", &json!({"body": "Nope"})).is_err());
        assert!(s.execute("mesh.plan", &json!({"body": "B", "layers": 0})).is_err());
        assert!(s.execute("mesh.plan", &json!({"body": "B", "size": f64::MAX})).is_ok());
        assert!(s.execute("mesh.plan.apply", &json!({"plan": "Nope"})).is_err());
    }

    /// A body partitions have cut is planned whole, and a cell of it on its own (from the features
    /// that built the body, with the rings the cell holds) (#473).
    #[test]
    fn a_partitioned_body_and_its_cells_can_be_planned() {
        let (_d, mut s) = open(json!([
            {"command": "solid.box", "params": {"length": 60, "width": 40, "height": 6, "body_name": "Plate"}},
            {"command": "solid.hole", "params": {"position": [15, 20, 6], "diameter": 8, "depth": 6}},
            {"command": "solid.hole", "params": {"position": [45, 20, 6], "diameter": 8, "depth": 6}},
        ]));
        s.execute("geometry.datum", &json!({"name": "Mid", "by": {"type": "point_normal", "point": [30, 0, 0], "normal": [1, 0, 0]}})).unwrap();
        s.execute("geometry.partition", &json!({"name": "Half", "datum": "Mid"})).unwrap();
        let whole = s.execute("mesh.plan", &json!({"body": "Plate"})).unwrap();
        assert_eq!(whole["rings"].as_array().map(Vec::len), Some(2), "{whole}");
        let cell = s.execute("mesh.plan", &json!({"body": "Plate/a"})).unwrap();
        assert_eq!(cell["rings"].as_array().map(Vec::len), Some(1), "{cell}");
        assert_eq!(cell["mesher"], "project", "{cell}");
        let r = s.execute("mesh.plan.apply", &json!({"plan": cell["name"]})).unwrap();
        assert_eq!(r["all_hex"], true, "{r}");
        assert!(r["min_scaled_jacobian"].as_f64().unwrap() >= 0.5, "{r}");
        // The whole body's plan cannot mesh a body partitions have cut: it says to plan a cell.
        let e = s.execute("mesh.plan.apply", &json!({"plan": whole["name"]})).unwrap_err().to_string();
        assert!(e.contains("cut into the cells") && e.contains("Plate/a"), "{e}");
    }

    /// The same from a saved model opened again (the design it refers to rebuilt on opening).
    #[test]
    fn a_saved_model_partitioned_after_opening_plans_its_body_and_cells() {
        let bytes = fieldcraft_geometry::script_design(
            &json!({"commands": [
                {"command": "solid.box", "params": {"length": 60, "width": 40, "height": 6, "body_name": "Plate"}},
                {"command": "solid.hole", "params": {"position": [15, 20, 6], "diameter": 8, "depth": 6}},
            ]})
            .to_string(),
        )
        .unwrap();
        let (_d, path) = crate::tests::design_file(&bytes);
        let model = std::path::Path::new(&path).with_file_name("model.fieldcraft").to_string_lossy().into_owned();
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        s.execute("file.save_as", &json!({"path": model})).unwrap();
        let mut s = Session::default();
        s.execute("file.open", &json!({"path": model})).unwrap();
        s.execute("geometry.datum", &json!({"name": "SymX", "by": {"type": "offset", "plane": "yz", "offset": 30}})).unwrap();
        s.execute("geometry.partition", &json!({"name": "HalfX", "datum": "SymX", "bodies": ["Plate"]})).unwrap();
        let cell = s.execute("mesh.plan", &json!({"body": "Plate/a"})).unwrap();
        assert_eq!(cell["rings"].as_array().map(Vec::len), Some(1), "{cell}");
        let whole = s.execute("mesh.plan", &json!({"body": "Plate"})).unwrap();
        assert_eq!(whole["body"], "Plate", "{whole}");
        // The plate is symmetric too, but the user's partition has cut it: refused, not mirrored.
        let e = s.execute("mesh.plan.apply", &json!({"plan": whole["name"]})).unwrap_err().to_string();
        assert!(e.contains("cut into the cells"), "{e}");
        s.execute("mesh.plan.apply", &json!({"plan": cell["name"]})).unwrap();
    }

    /// The feature graph of every corpus design and the crank arm: every body lists the features
    /// that built it, and a plan builds for every body (fillets, revolves and lofts unmapped).
    #[test]
    fn every_corpus_design_has_a_feature_graph_and_a_plan() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/corpus");
        let mut scripts: Vec<(String, String)> = std::fs::read_dir(&dir)
            .unwrap()
            .filter_map(|e| {
                let path = e.ok()?.path().join("recipe.json");
                Some((path.parent()?.file_name()?.to_string_lossy().into_owned(), std::fs::read_to_string(path).ok()?))
            })
            .collect();
        scripts.push(("crank-arm".into(), crate::sample::CRANK_ARM.into()));
        scripts.sort();
        assert!(scripts.len() > 10);
        for (name, script) in scripts {
            let (_d, mut s) = open_script(&script);
            let graph = s.execute("geometry.features", &json!({})).unwrap();
            for b in graph["bodies"].as_array().unwrap() {
                let features = b["features"].as_array().unwrap();
                assert!(!features.is_empty(), "{name}: {b}");
                // Every feature's faces are faces of the body.
                let body = b["body"].as_str().unwrap();
                let plan = s.execute("mesh.plan", &json!({"body": body}));
                let plan = plan.unwrap_or_else(|e| panic!("{name}/{body}: {e}"));
                let mapped = plan["features"].as_array().map_or(0, Vec::len);
                let unmapped = plan["unmapped"].as_array().map_or(0, Vec::len);
                assert_eq!(mapped + unmapped, features.len(), "{name}/{body}: {plan}");
            }
        }
    }
}

/// Acceptance of the projection layouts (#402 step 3): parts built the way enclosures are, each
/// meshed all-hex from its plan alone, conforming, holes round.
#[cfg(test)]
mod acceptance {
    use serde_json::{Value, json};

    use super::tests_support::*;
    use crate::Session;

    #[test]
    fn a_plate_with_four_holes_has_round_o_grid_rings() {
        let (_d, mut s) = open(json!([
            {"command": "solid.box", "params": {"length": 80, "width": 60, "height": 8, "body_name": "Plate"}},
            {"command": "solid.hole", "params": {"position": [15, 15, 8], "diameter": 8, "depth": 8}},
            {"command": "solid.hole", "params": {"position": [65, 15, 8], "diameter": 8, "depth": 8}},
            {"command": "solid.hole", "params": {"position": [15, 45, 8], "diameter": 8, "depth": 8}},
            {"command": "solid.hole", "params": {"position": [65, 45, 8], "diameter": 8, "depth": 8}},
        ]));
        let r = plan_and_apply(&mut s, "Plate");
        assert!(r["layout"].as_str().unwrap().contains("4 circle(s) ringed"), "{r}");
        accept(&mut s, &r, 0.5, 10.0);
    }

    /// An enclosure lid: a plate with rounded corners, a screw hole in each corner close to the
    /// fillet and a square opening. Each hole's ring reaches the plate's sides and takes the
    /// corner fillet in (no blend layer, no fan into the hole), its frame centred on the hole.
    #[test]
    fn a_lid_with_screw_holes_in_its_rounded_corners_rings_each_hole_to_the_fillet() {
        let (_d, mut s) = open(json!([
            {"command": "solid.box", "params": {"length": 128, "width": 98, "height": 4, "body_name": "Lid"}},
            {"command": "solid.fillet", "params": {"edges": [[0, 0, 2], [128, 0, 2], [128, 98, 2], [0, 98, 2]], "radius": 8}},
            {"command": "solid.hole", "params": {"position": [7, 7, 4], "diameter": 3.5, "depth": 4}},
            {"command": "solid.hole", "params": {"position": [121, 7, 4], "diameter": 3.5, "depth": 4}},
            {"command": "solid.hole", "params": {"position": [7, 91, 4], "diameter": 3.5, "depth": 4}},
            {"command": "solid.hole", "params": {"position": [121, 91, 4], "diameter": 3.5, "depth": 4}},
            {"command": "solid.box", "params": {"length": 22, "width": 22, "height": 4, "corner": [41, 40, 0], "operation": "cut"}},
        ]));
        let plan = s.execute("mesh.plan", &json!({"body": "Lid", "size": 1.5, "layers": 2})).unwrap();
        assert_eq!(plan["mesher"], "project", "{plan}");
        // Every frame is the same square, centred on its hole and reaching the plate's sides.
        let rings = plan["rings"].as_array().unwrap();
        assert_eq!(rings.len(), 4, "{plan}");
        for r in rings {
            assert!((r["frame"].as_f64().unwrap() - 7.0).abs() < 1e-6, "{r}");
        }
        let notes = plan["notes"].as_array().unwrap();
        assert!(notes.iter().all(|n| !n.as_str().unwrap().contains("frame cut")), "{plan}");
        let r = s.execute("mesh.plan.apply", &json!({"plan": plan["name"]})).unwrap();
        assert!(r["layout"].as_str().unwrap().contains("4 circle(s) ringed"), "{r}");
        assert_eq!(r["blend_layers"], json!([]), "the rings take the corner fillets in: {r}");
        accept(&mut s, &r, 0.6, 6.0);
    }

    /// The enclosure base of the test below without its standoffs: symmetric across x = 28 and
    /// y = 21, so the plan meshes one quarter and mirrors it: the four corners mesh exactly
    /// alike, the smallest scaled Jacobian round each corner boss well above 0.5.
    #[test]
    fn an_enclosure_base_meshes_its_four_corners_alike() {
        let mut commands = vec![
            json!({"command": "solid.box", "params": {"length": 56, "width": 42, "height": 16, "body_name": "Base"}}),
            json!({"command": "solid.fillet", "params": {"edges": [[0, 0, 8], [56, 0, 8], [56, 42, 8], [0, 42, 8]], "radius": 5}}),
            json!({"command": "solid.shell", "params": {"faces": [[28, 21, 16]], "thickness": 3, "body": "Base"}}),
        ];
        let corners = [(-1.0, -1.0), (57.0, -1.0), (-1.0, 43.0), (57.0, 43.0)];
        for (x, y) in corners {
            commands.push(json!({"command": "solid.cylinder", "params": {"diameter": 13, "height": 16, "base": [x, y, 0], "operation": "join", "targets": ["Base"]}}));
        }
        for (x, y) in corners {
            commands.push(json!({"command": "solid.hole", "params": {"position": [x, y, 16], "diameter": 3.3, "depth": 8}}));
        }
        let (_d, mut s) = open(json!(commands));
        let plan = s.execute("mesh.plan", &json!({"body": "Base", "size": 1.5, "layers": 2})).unwrap();
        let r = s.execute("mesh.plan.apply", &json!({"plan": plan["name"]})).unwrap();
        accept(&mut s, &r, 0.5, 10.0);
        // Every element's scaled Jacobian, by the corner (quadrant) its centroid lies in.
        let m = s.execute("mesh.quality.mark", &json!({"limits": {"scaled_jacobian": 1.01, "min_dihedral": null, "max_dihedral": null, "aspect_ratio": null, "min_angle": null, "max_angle": null}})).unwrap();
        let mut worst = [f64::INFINITY; 4];
        for e in m["elements"].as_array().unwrap() {
            let at = &e["at"];
            let q = usize::from(at[0].as_f64().unwrap() > 28.0) + 2 * usize::from(at[1].as_f64().unwrap() > 21.0);
            for f in e["failed"].as_array().unwrap() {
                if f["metric"] == "scaled_jacobian" {
                    worst[q] = worst[q].min(f["value"].as_f64().unwrap());
                }
            }
        }
        // Meshed by its symmetry: one quarter, mirrored across both planes.
        assert_eq!(r["symmetry"].as_array().map(Vec::len), Some(2), "{r}");
        assert_eq!(r["cells"].as_array().map(Vec::len), Some(4), "{r}");
        let (lo, hi) = worst.iter().fold((f64::INFINITY, 0.0f64), |(l, h), w| (l.min(*w), h.max(*w)));
        assert!(hi - lo < 1e-9, "per corner min scaled Jacobian {worst:?}: {r}");
        assert!(lo > 0.65, "per corner min scaled Jacobian {worst:?}");
    }

    /// The enclosure base meshed by its symmetry (one quarter, mirrored twice): each corner
    /// boss cuts the inner fillet into two pieces and a sliver narrower than the elements, and
    /// every face of every cell, the slivers and their mirror images among them, has element
    /// faces tagged with it, so a set by any face finds the mesh (#494).
    #[test]
    fn a_symmetric_filleted_base_tags_every_face_of_every_cell() {
        let mut commands = vec![
            json!({"command": "solid.box", "params": {"length": 56, "width": 42, "height": 16, "body_name": "Base"}}),
            json!({"command": "solid.fillet", "params": {"edges": [[0, 0, 8], [56, 0, 8], [56, 42, 8], [0, 42, 8]], "radius": 5}}),
            json!({"command": "solid.shell", "params": {"faces": [[28, 21, 16]], "thickness": 3, "body": "Base"}}),
        ];
        for (x, y) in [(-1.0, -1.0), (57.0, -1.0), (-1.0, 43.0), (57.0, 43.0)] {
            commands.push(json!({"command": "solid.cylinder", "params": {"diameter": 13, "height": 16, "base": [x, y, 0], "operation": "join", "targets": ["Base"]}}));
        }
        let (_d, mut s) = open(json!(commands));
        let plan = s.execute("mesh.plan", &json!({"body": "Base", "size": 1.5, "layers": 2})).unwrap();
        let r = s.execute("mesh.plan.apply", &json!({"plan": plan["name"]})).unwrap();
        assert_eq!(r["cells"].as_array().map(Vec::len), Some(4), "{r}");
        assert_eq!(r["untagged_faces"], json!([]), "{r}");
        let design = s.design.clone().unwrap();
        let mut blends = 0;
        for cell in r["cells"].as_array().unwrap() {
            let body = design.body(cell.as_str().unwrap()).unwrap();
            for (k, face) in body.faces().iter().enumerate() {
                blends += usize::from(face.name.contains("blend"));
                let set = format!("F{k}");
                s.execute("model.set.create", &json!({"name": set, "geometry": [{"body": body.name, "face": face.name}]})).unwrap();
                let faces = s.model.resolve_set(&set).map(|x| x.faces.len());
                assert!(faces.as_ref().is_ok_and(|n| *n > 0), "face {} of {}: {faces:?}", face.name, body.name);
                s.execute("model.set.delete", &json!({"name": set})).unwrap();
            }
        }
        assert!(blends >= 4, "the inner fillet's pieces");
    }

    /// An enclosure base: a box with rounded vertical edges, shelled open at the top, a screw
    /// boss joined at each corner (centred outside it, nearly touching the inner fillet) with a
    /// blind hole, and standoffs on the floor. Each boss's ring frames the whole circle,
    /// subdivided where the walls join it; the inner fillet is taken in by the ring's corner and
    /// the boss's near touch with it leaves no sliver.
    #[test]
    fn an_enclosure_base_rings_each_corner_boss_where_the_walls_join_it() {
        let mut commands = vec![
            json!({"command": "solid.box", "params": {"length": 56, "width": 42, "height": 16, "body_name": "Base"}}),
            json!({"command": "solid.fillet", "params": {"edges": [[0, 0, 8], [56, 0, 8], [56, 42, 8], [0, 42, 8]], "radius": 5}}),
            json!({"command": "solid.shell", "params": {"faces": [[28, 21, 16]], "thickness": 3, "body": "Base"}}),
        ];
        for (x, y) in [(-1.0, -1.0), (57.0, -1.0), (-1.0, 43.0), (57.0, 43.0)] {
            commands.push(json!({"command": "solid.cylinder", "params": {"diameter": 13, "height": 16, "base": [x, y, 0], "operation": "join", "targets": ["Base"]}}));
        }
        for (x, y) in [(-1.0, -1.0), (57.0, -1.0), (-1.0, 43.0), (57.0, 43.0)] {
            commands.push(json!({"command": "solid.hole", "params": {"position": [x, y, 16], "diameter": 3.3, "depth": 8}}));
        }
        for (x, y) in [(16.5, 16.5), (39.5, 16.5)] {
            commands.push(json!({"command": "solid.cylinder", "params": {"diameter": 6, "height": 4, "base": [x, y, 3], "operation": "join", "targets": ["Base"]}}));
        }
        let (_d, mut s) = open(json!(commands));
        let plan = s.execute("mesh.plan", &json!({"body": "Base", "size": 1.5, "layers": 2})).unwrap();
        assert_eq!(plan["mesher"], "project", "{plan}");
        let r = s.execute("mesh.plan.apply", &json!({"plan": plan["name"]})).unwrap();
        // Symmetric across x = 28 only (the standoffs stand on one side of y = 21): one half
        // meshed (two bosses, their holes and a standoff ringed) and mirrored.
        assert_eq!(r["symmetry"].as_array().map(Vec::len), Some(1), "{r}");
        assert!(r["layout"].as_str().unwrap().contains("5 circle(s) ringed"), "{r}");
        accept(&mut s, &r, 0.5, 10.0);
    }

    /// The enclosure base's bosses cut a 0.42 mm strip out of the middle of the inner fillet
    /// (#511, #513): geometry.check finds the slivers and geometry.face.absorb absorbs them.
    /// The absorbed faces stay absorbed in the cells the plan's symmetry cuts the body into, and
    /// the projection keeps to the exact surfaces of a merged face's pieces, so the mesh is no
    /// worse for them, whole or by symmetry.
    #[test]
    fn absorbed_slivers_stay_absorbed_in_the_plan_s_cells() {
        let mut commands = vec![
            json!({"command": "solid.box", "params": {"length": 56, "width": 42, "height": 16, "body_name": "Base"}}),
            json!({"command": "solid.fillet", "params": {"edges": [[0, 0, 8], [56, 0, 8], [56, 42, 8], [0, 42, 8]], "radius": 5}}),
            json!({"command": "solid.shell", "params": {"faces": [[28, 21, 16]], "thickness": 3, "body": "Base"}}),
        ];
        for (x, y) in [(-1.0, -1.0), (57.0, -1.0), (-1.0, 43.0), (57.0, 43.0)] {
            commands.push(json!({"command": "solid.cylinder", "params": {"diameter": 13, "height": 16, "base": [x, y, 0], "operation": "join", "targets": ["Base"]}}));
        }
        let (_d, mut s) = open(json!(commands));
        let check = s.execute("geometry.check", &json!({"body": "Base", "size": 1.5, "mark": true})).unwrap();
        assert_eq!(check["counts"]["sliver_face"], 4, "{check}");
        let absorbed = s.execute("geometry.face.absorb", &json!({"from_check": true})).unwrap();
        assert_eq!(absorbed["absorbed"].as_array().map(Vec::len), Some(4), "{absorbed}");
        let plan = s.execute("mesh.plan", &json!({"body": "Base", "size": 1.5, "layers": 2})).unwrap();
        assert_eq!(plan["mesher"], "project", "{plan}");
        let r = s.execute("mesh.plan.apply", &json!({"plan": plan["name"], "symmetry": false})).unwrap();
        accept(&mut s, &r, 0.6, 10.0);
        let r = s.execute("mesh.plan.apply", &json!({"plan": plan["name"]})).unwrap();
        assert_eq!(r["symmetry"].as_array().map(Vec::len), Some(2), "{r}");
        let design = s.design.as_deref().unwrap();
        assert!(design.errors.is_empty(), "{:?}", design.errors);
        // Each quarter holds one corner, its sliver absorbed.
        let cells: Vec<_> = design.bodies.iter().filter(|b| b.name.starts_with("Base/")).collect();
        assert_eq!(cells.len(), 4);
        assert!(cells.iter().all(|c| c.merged_faces().len() == 1), "{:?}", cells.iter().map(|c| c.merged_faces()).collect::<Vec<_>>());
        assert_eq!(design.source_body("Base").map(|b| b.merged_faces().len()), Some(4));
        accept(&mut s, &r, 0.65, 10.0);
    }

    /// The double 0.4 mm step's slivers merge into faces of other orientations (#546): before the
    /// fix the merged faces classified free-form and the plan lost its cells layout and fell back
    /// to sweep. An absorbed face keeps the classification of the face it was absorbed into, so
    /// the plan keeps its layout and the hex route survives the clean-up.
    #[test]
    fn absorbed_step_slivers_keep_the_cells_plan_and_the_hex_route() {
        let (_d, mut s) = open(json!([
            {"command": "solid.box", "params": {"corner": [0, 0, 0], "length": 60, "width": 40, "height": 10, "body_name": "Steps"}},
            {"command": "solid.box", "params": {"corner": [0, 0, 10], "length": 60, "width": 0.4, "height": 0.4, "operation": "join", "targets": ["Steps"]}},
            {"command": "solid.box", "params": {"corner": [0, 0, 10], "length": 0.4, "width": 40, "height": 0.4, "operation": "join", "targets": ["Steps"]}}
        ]));
        let plan = s.execute("mesh.plan", &json!({"body": "Steps", "size": 2.0, "layers": 2})).unwrap();
        assert_eq!(plan["mesher"], "cells", "{plan}");
        let check = s.execute("geometry.check", &json!({"body": "Steps", "size": 2.0, "mark": true})).unwrap();
        assert_eq!(check["counts"]["sliver_face"], 3, "{check}");
        let absorbed = s.execute("geometry.face.absorb", &json!({"from_check": true})).unwrap();
        assert_eq!(absorbed["absorbed"].as_array().map(Vec::len), Some(3), "{absorbed}");
        let plan = s.execute("mesh.plan", &json!({"body": "Steps", "size": 2.0, "layers": 2})).unwrap();
        assert_eq!(plan["mesher"], "cells", "{plan}");
        let hex = s.execute("mesh.hex.generate", &json!({"body": "Steps", "size": 2.0})).unwrap();
        let sj = hex.pointer("/bodies/0/quality/min_scaled_jacobian").and_then(|v| v.as_f64()).unwrap_or(0.0);
        assert!(sj >= 0.5, "{hex}");
    }

    #[test]
    fn curved_bodies_mesh_within_a_percent_of_their_volume_at_the_plan_s_size() {
        // The plan's size keeps the chords of a torus (swept) and a cylinder (butterfly) within
        // the volume budget; the butterfly's rim is divided by its length, not its class's mean.
        for cmds in [
            json!([{"command": "solid.torus", "params": {"major": 30.0, "minor": 8.0, "body_name": "B"}}]),
            json!([{"command": "solid.cylinder", "params": {"radius": 10.0, "height": 40.0, "body_name": "B"}}]),
        ] {
            let (_d, mut s) = open(cmds.clone());
            let body = s.execute("geometry.bodies", &json!({})).unwrap()[0]["volume"].as_f64().unwrap();
            let plan = s.execute("mesh.plan", &json!({"body": "B"})).unwrap();
            assert!(plan["notes"].to_string().contains("the chords of"), "{plan}");
            s.execute("mesh.plan.apply", &json!({"plan": plan["name"]})).unwrap();
            let mesh = s.execute("mesh.quality", &json!({})).unwrap()["volume"].as_f64().unwrap();
            assert!(((mesh - body) / body).abs() < 0.01, "{cmds}: mesh {mesh} against {body}");
        }
    }

    /// The crank arm sample, whole (#455, #190): its sides are at an angle to the axes and its
    /// root arcs meet the bosses at a kink, so it takes the ring-and-pillow layout: each boss an
    /// annulus from its bore or pedal hole to its circle through every level, the arm a carved
    /// grid with a pillow layer to the outline and the bosses' circles, a corner block at each
    /// root. One conforming all-hex mesh from the plan alone, at three sizes.
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn the_crank_arm_meshes_from_its_plan_with_rings_and_a_pillow() {
        for size in [2.5, 3.0, 4.0] {
            let mut s = Session::default();
            s.execute("file.sample", &json!({"name": "crank-arm"})).unwrap();
            s.execute("mesh.clear", &json!({})).unwrap();
            s.execute("geometry.partition.delete", &json!({"name": "Half"})).unwrap();
            let plan = s.execute("mesh.plan", &json!({"body": "CrankArm", "size": size})).unwrap();
            assert_eq!(plan["mesher"], "project", "size {size}: {plan}");
            assert!(plan["cuts"].as_array().is_none_or(Vec::is_empty), "{plan}");
            let r = s.execute("mesh.plan.apply", &json!({"plan": plan["name"]})).unwrap();
            let layout = r["layout"].as_str().unwrap();
            assert!(layout.contains("2 circle(s) ringed") && layout.contains("pillow layer"), "{r}");
            // The rims are divided by their own length (no chord longer than the size), which
            // thins the bore's ring elements a little at 2.5 mm (edge ratio 7.1).
            accept(&mut s, &r, 0.5, 7.5);
        }
    }

    #[test]
    fn a_filleted_box_meshes_structured_round_its_corners() {
        let (_d, mut s) = open(json!([
            {"command": "solid.box", "params": {"length": 40, "width": 30, "height": 20, "body_name": "Box"}},
            {"command": "solid.fillet", "params": {"edges": [[0, 0, 10], [40, 0, 10], [40, 30, 10], [0, 30, 10]], "radius": 4}},
        ]));
        let r = plan_and_apply(&mut s, "Box");
        assert_eq!(r["blend_layers"].as_array().unwrap().len(), 4, "{r}");
        accept(&mut s, &r, 0.5, 10.0);
    }

    #[test]
    fn a_box_filleted_on_every_edge_takes_a_blend_layer_on_every_side() {
        let recipe: Value = serde_json::from_str(
            &std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/corpus/rounded-box/recipe.json")).unwrap(),
        )
        .unwrap();
        let (_d, mut s) = open(recipe["commands"].clone());
        let r = plan_and_apply(&mut s, "Box");
        assert_eq!(r["blend_layers"].as_array().unwrap().len(), 6, "{r}");
        accept(&mut s, &r, 0.5, 10.0);
    }
}

/// What the plan tests share.
#[cfg(test)]
mod tests_support {
    use serde_json::{Value, json};

    use crate::Session;

    pub(super) fn open(commands: Value) -> (crate::tests::TempDir, Session) {
        let bytes = fieldcraft_geometry::script_design(&json!({"commands": commands}).to_string()).unwrap();
        let (dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        (dir, s)
    }

    /// Plan `body` (it must take the projection layout) and apply the plan to the whole body.
    pub(super) fn plan_and_apply(s: &mut Session, body: &str) -> Value {
        let plan = s.execute("mesh.plan", &json!({"body": body})).unwrap();
        assert_eq!(plan["mesher"], "project", "{plan}");
        assert!(plan["unmapped"].as_array().is_none_or(Vec::is_empty), "{plan}");
        s.execute("mesh.plan.apply", &json!({"plan": plan["name"], "symmetry": false})).unwrap()
    }

    /// All hex, one conforming mesh on the geometry, the volume the body's, no element below
    /// `min_sj` or above `max_edge_ratio`, no free faces inside, every face of the body tagged.
    pub(super) fn accept(s: &mut Session, r: &Value, min_sj: f64, max_edge_ratio: f64) {
        assert_eq!((r["all_hex"].as_bool(), r["conforming"].as_bool()), (Some(true), Some(true)), "{r}");
        assert_eq!(r["untagged_faces"], json!([]), "faces no element face is tagged with");
        let q = s.execute("mesh.quality", &json!({})).unwrap();
        let sj = q["metrics"]["scaled_jacobian"]["min"].as_f64().unwrap();
        let er = q["metrics"]["edge_ratio"]["max"].as_f64().unwrap();
        assert!(sj >= min_sj, "min scaled Jacobian {sj}: {q}");
        assert!(er <= max_edge_ratio, "max edge ratio {er}: {q}");
        assert!(r["mesh"]["max_gap"].as_f64().unwrap() < 0.05, "boundary nodes off the geometry: {r}");
        let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
        let geometry: f64 = bodies.as_array().unwrap().iter().map(|b| b["volume"].as_f64().unwrap()).sum();
        let mesh = q["volume"].as_f64().unwrap();
        assert!((mesh - geometry).abs() < 0.01 * geometry, "mesh volume {mesh} vs {geometry}");
        let dup = s.execute("mesh.check.duplicates", &json!({})).unwrap();
        assert_eq!((dup["counts"]["coincident"].as_u64(), dup["counts"]["cracks"].as_u64()), (Some(0), Some(0)), "{dup}");
    }
}
