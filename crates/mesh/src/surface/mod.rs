//! Boundary surface preparation: the step between a tessellated body and volume meshing.
//!
//! [`BoundarySurface::new`] turns a face-tagged triangle surface (SolidCraft's tessellation, or a
//! tagged 3MF) into one that every volume mesher can rely on, or fails with the list of
//! [`Defect`]s, each located by input triangle or vertex index and a position:
//!
//! 1. **Coordinates and indices** are checked: finite, within [`fieldcraft_model::MAX_COORD`],
//!    in range.
//! 2. **Weld.** Vertices closer than a distance (relative to the bounding-box diagonal, or
//!    absolute: [`SurfaceOptions`]) become one, in input order (the earliest vertex of a cluster
//!    is its representative and keeps its position), so the result does not depend on hashing.
//!    SolidCraft's tessellation is a triangle soup (three vertices per triangle), and two faces'
//!    copies of a shared boundary point agree only to within the kernel's tolerance, so this is
//!    where its faces join. Triangles that collapse are dropped; coordinates below
//!    [`crate::predicates::MIN_MAGNITUDE`] snap to zero so every point is in the exact
//!    predicates' domain.
//! 3. **Degenerate triangles** (exactly collinear corners) are rejected.
//! 4. **Topology.** Every edge has exactly two triangles that traverse it in opposite directions
//!    (closed, edge-manifold, consistently oriented), and the triangles around each vertex form a
//!    single fan (vertex-manifold).
//! 5. **Self-intersection.** A bounding volume hierarchy pairs triangles whose boxes touch; each
//!    pair is decided exactly ([`intersect`]): triangles may meet only in the vertices or the edge
//!    they share.
//! 6. **Shells and orientation.** Connected components, their genus from the Euler
//!    characteristic, their signed volume (divergence theorem), and their nesting from the
//!    generalized winding number: an outer shell must enclose positive volume, a void (a shell
//!    inside an odd number of others) negative.
//! 7. **Feature edges**: edges between different B-rep faces, and edges whose dihedral angle
//!    exceeds a threshold. The B-rep face of every triangle is kept, so loads and boundary
//!    conditions on CAD faces map to mesh faces.
//!
//! # Sources
//! - J. R. Shewchuk, *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric
//!   Predicates*, DCG 18:305–363, 1997: every intersection decision is a sign of `orient3d` or
//!   `orient2d` ([`crate::predicates`]).
//! - P. Guigue, O. Devillers, *Fast and Robust Triangle–Triangle Overlap Test Using Orientation
//!   Predicates*, J. Graphics Tools 8(1):25–32, 2003: triangle overlap decided from orientation
//!   signs alone (here through the edge–triangle reduction in [`intersect`]).
//! - C. Ericson, *Real-Time Collision Detection*, Morgan Kaufmann, 2005, ch. 6: bounding volume
//!   hierarchies, top-down median-split construction and the self-collision traversal.
//! - A. Jacobson, L. Kavan, O. Sorkine-Hornung, *Robust Inside-Outside Segmentation using
//!   Generalized Winding Numbers*, ACM TOG 32(4), 2013, with the triangle solid angle of A. van
//!   Oosterom, J. Strackee, *The Solid Angle of a Plane Triangle*, IEEE TBME 30(2):125–126, 1983:
//!   shell nesting.

mod bvh;
pub mod intersect;
mod topology;
pub(crate) mod weld;

#[cfg(test)]
mod tests;

use std::fmt;

use fieldcraft_model::{MAX_COORD, MAX_ELEMENTS, MAX_NODES};

use crate::message;
use crate::message::{Message, length, point, volume};
use crate::predicates::Point3;
use crate::vec3::{cross, norm, sub};
use crate::{Error, Result};

/// The default relative weld distance.
const DEFAULT_WELD: f64 = 1e-9; // unit-factors: relative to the diagonal

/// Defects listed in an error; the count covers the rest.
pub const MAX_LISTED_DEFECTS: usize = 100;
/// Cap on exact triangle-pair tests in the self-intersection check.
pub const MAX_PAIR_TESTS: u64 = 1_000_000_000;
/// Cap on triangle solid angles summed to nest shells (shells × triangles, after box culling).
pub const MAX_WINDING_WORK: u64 = 500_000_000;

/// The largest weld distance, as a fraction of the bounding-box diagonal.
pub const MAX_RELATIVE_WELD: f64 = 1e-3; // unit-factors: relative to the diagonal

/// How [`BoundarySurface::new`] welds and finds features.
///
/// Vertices closer than the larger of `weld_tolerance` × the bounding-box diagonal and
/// `weld_distance` are welded. Tessellations of B-rep faces agree on their shared boundary points
/// only to within the kernel's tolerance, so a caller with CAD input sets `weld_distance` from it
/// (twice the kernel's linear tolerance: each face's copy of a boundary point is within it of the
/// exact point).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceOptions {
    /// Relative weld distance, in `[0, 1e-3]`; nonzero values below 1e-15 are refused (finer than
    /// the coordinates themselves). 0 with `weld_distance` 0 welds identical coordinates only.
    pub weld_tolerance: f64,
    /// Absolute weld distance in model units, at least 0 (the combined distance may be at most
    /// [`MAX_RELATIVE_WELD`] of the bounding-box diagonal). A tessellation carries its own
    /// (`fieldcraft_geometry::SurfaceMesh::weld_distance`).
    pub weld_distance: f64,
    /// An edge inside one B-rep face is a feature when the normals of its two triangles differ
    /// by more than this angle, in degrees, in `(0, 180]` (180: face boundaries only).
    pub feature_angle: f64,
}

impl Default for SurfaceOptions {
    fn default() -> Self {
        SurfaceOptions { weld_tolerance: DEFAULT_WELD, weld_distance: 0.0, feature_angle: 45.0 }
    }
}

/// A defect that keeps a triangle surface from bounding a volume mesh. Triangle and vertex
/// numbers are indices into the input; `at` is a position near the defect.
#[derive(Debug, Clone, PartialEq)]
pub enum Defect {
    /// A coordinate is not a finite number within ±[`MAX_COORD`].
    BadCoordinate { vertex: u32 },
    /// A triangle refers to a vertex that does not exist.
    BadIndex { triangle: u32 },
    /// A triangle with exactly collinear corners.
    Degenerate { triangle: u32, at: [f64; 3] },
    /// An edge with only one triangle: the surface has a hole there.
    OpenEdge { vertices: [u32; 2], triangle: u32, at: [f64; 3] },
    /// An edge shared by more than two triangles.
    NonManifoldEdge { vertices: [u32; 2], triangles: u32, at: [f64; 3] },
    /// Two triangles traverse their shared edge in the same direction: one of them is flipped.
    Flipped { vertices: [u32; 2], triangles: [u32; 2], at: [f64; 3] },
    /// The triangles around a vertex form more than one fan (two sheets touch at a point).
    NonManifoldVertex { vertex: u32, at: [f64; 3] },
    /// Two triangles meet outside the vertices and edge they share.
    SelfIntersection { triangles: [u32; 2], at: [f64; 3] },
    /// A shell enclosing negative volume although it is not inside another shell (its triangles
    /// face inwards), or positive volume although it is a void.
    Inverted { shell: u32, triangle: u32, volume: f64 },
    /// A shell whose enclosed volume is too small to tell its orientation from rounding.
    Flat { shell: u32, triangle: u32 },
}

impl Defect {
    /// A short, stable name for the kind of defect (`open_edge`, `self_intersection`, …).
    pub fn kind(&self) -> &'static str {
        match self {
            Defect::BadCoordinate { .. } => "bad_coordinate",
            Defect::BadIndex { .. } => "bad_index",
            Defect::Degenerate { .. } => "degenerate",
            Defect::OpenEdge { .. } => "open_edge",
            Defect::NonManifoldEdge { .. } => "non_manifold_edge",
            Defect::Flipped { .. } => "flipped",
            Defect::NonManifoldVertex { .. } => "non_manifold_vertex",
            Defect::SelfIntersection { .. } => "self_intersection",
            Defect::Inverted { .. } => "inverted",
            Defect::Flat { .. } => "flat",
        }
    }

    /// A position near the defect, when it has one.
    pub fn at(&self) -> Option<[f64; 3]> {
        match self {
            Defect::BadCoordinate { .. } | Defect::BadIndex { .. } | Defect::Inverted { .. } | Defect::Flat { .. } => None,
            Defect::Degenerate { at, .. }
            | Defect::OpenEdge { at, .. }
            | Defect::NonManifoldEdge { at, .. }
            | Defect::Flipped { at, .. }
            | Defect::NonManifoldVertex { at, .. }
            | Defect::SelfIntersection { at, .. } => Some(*at),
        }
    }
}

impl Defect {
    /// What the defect is, its points as data.
    pub fn message(&self) -> Message {
        match *self {
            Defect::BadCoordinate { vertex } => {
                message![format!("vertex {vertex}: a coordinate is not a finite number within ±"), length(MAX_COORD)]
            }
            Defect::BadIndex { triangle } => message![format!("triangle {triangle} refers to a vertex that does not exist")],
            Defect::Degenerate { triangle, at } => message![format!("triangle {triangle} at "), point(at), " has collinear corners"],
            Defect::OpenEdge { vertices: [a, b], triangle, at } => {
                message![format!("open edge {a}–{b} at "), point(at), format!(" (only triangle {triangle} uses it)")]
            }
            Defect::NonManifoldEdge { vertices: [a, b], triangles, at } => {
                message![format!("edge {a}–{b} at "), point(at), format!(" is shared by {triangles} triangles")]
            }
            Defect::Flipped { vertices: [a, b], triangles: [s, t], at } => {
                message![format!("triangles {s} and {t} are oriented inconsistently across edge {a}–{b} at "), point(at)]
            }
            Defect::NonManifoldVertex { vertex, at } => message![format!("vertex {vertex} at "), point(at), " joins separate sheets"],
            Defect::SelfIntersection { triangles: [s, t], at } => message![format!("triangles {s} and {t} intersect near "), point(at)],
            Defect::Inverted { shell, triangle, volume: v } => message![
                format!("shell {shell} (from triangle {triangle}) is inside out: it encloses "),
                volume(v),
                " where the opposite sign is expected"
            ],
            Defect::Flat { shell, triangle } => message![format!("shell {shell} (from triangle {triangle}) encloses no measurable volume")],
        }
    }
}

impl fmt::Display for Defect {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.message().fmt(f)
    }
}

/// The defects found by the first check that failed: up to [`MAX_LISTED_DEFECTS`], and how many
/// there are in all (a lower bound for self-intersections, whose search stops at the cap).
#[derive(Debug, Clone, PartialEq)]
pub struct Defects {
    pub listed: Vec<Defect>,
    pub total: usize,
}

impl Defects {
    fn new() -> Self {
        Defects { listed: Vec::new(), total: 0 }
    }

    fn push(&mut self, d: Defect) {
        self.total += 1;
        if self.listed.len() < MAX_LISTED_DEFECTS {
            self.listed.push(d);
        }
    }

    fn is_full(&self) -> bool {
        self.listed.len() >= MAX_LISTED_DEFECTS
    }

    /// `Err` when there is at least one defect.
    fn into_result(self) -> Result<()> {
        if self.total == 0 { Ok(()) } else { Err(Error::Surface(self)) }
    }
}

impl Defects {
    /// The defects, their points as data.
    pub fn message(&self) -> Message {
        let s = if self.total == 1 { "" } else { "s" };
        let mut out = message![format!("{} defect{s}", self.total)];
        for d in self.listed.iter().take(3) {
            out = out.then(message!["; "]).then(d.message());
        }
        if self.total > 3 {
            out.0.push("; …".into());
        }
        out
    }

    /// The defects as text, their points shown in `units`.
    pub fn render(&self, units: &dyn crate::message::Units) -> String {
        self.message().render(units)
    }
}

impl fmt::Display for Defects {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.render(&crate::message::Bare))
    }
}

/// A connected component of the surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Shell {
    /// The shell's first triangle (output index); [`BoundarySurface::tri_shell`] lists the rest.
    pub first_triangle: u32,
    pub triangles: u32,
    pub vertices: u32,
    /// Handles (through holes): 0 for a block, 1 for a plate with a hole.
    pub genus: u32,
    /// Signed enclosed volume: positive for an outer shell, negative for a void.
    pub volume: f64,
    /// Inside an odd number of other shells: a cavity in the solid.
    pub void: bool,
}

/// Why an edge is a feature.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FeatureKind {
    /// Between two B-rep faces.
    FaceBoundary,
    /// Inside one face, with a dihedral angle above the threshold.
    Sharp,
}

/// An edge the volume mesh must keep: its vertices (ascending), its two triangles (the one that
/// traverses it from `vertices[0]` to `vertices[1]` first) and the angle between their normals
/// in degrees.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FeatureEdge {
    pub vertices: [u32; 2],
    pub triangles: [u32; 2],
    pub angle: f64,
    pub kind: FeatureKind,
}

/// A welded, closed, 2-manifold, outward-oriented, self-intersection-free triangle surface with
/// its B-rep face tags and feature edges: the boundary every volume mesher starts from.
///
/// Every point is admissible for [`crate::predicates`] ([`BoundarySurface::point3`]), triangles
/// are counter-clockwise seen from outside the solid (voids included: their normals point into
/// the cavity), and all output is a deterministic function of the input.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundarySurface {
    points: Vec<[f64; 3]>,
    triangles: Vec<[u32; 3]>,
    tri_face: Vec<u32>,
    neighbors: Vec<[u32; 3]>,
    tri_shell: Vec<u32>,
    shells: Vec<Shell>,
    feature_edges: Vec<FeatureEdge>,
    source_vertex: Vec<u32>,
    source_triangle: Vec<u32>,
}

impl BoundarySurface {
    /// Prepare `triangles` over `positions` (counter-clockwise seen from outside) with the B-rep
    /// face of each triangle (`tri_face`; empty means one face). Fails with
    /// [`Error::Surface`] listing the defects of the first check that fails, with
    /// [`Error::Invalid`] on unusable options or sizes, or with [`Error::Limit`] when a search cap
    /// is reached.
    pub fn new(positions: &[[f64; 3]], triangles: &[[u32; 3]], tri_face: &[u32], options: &SurfaceOptions) -> Result<Self> {
        check_options(options)?;
        if positions.len() > MAX_NODES || triangles.len() > MAX_ELEMENTS {
            return Err(Error::Invalid(format!(
                "{} vertices and {} triangles exceed the limits ({MAX_NODES}, {MAX_ELEMENTS})",
                positions.len(),
                triangles.len()
            )));
        }
        if !tri_face.is_empty() && tri_face.len() != triangles.len() {
            return Err(Error::Invalid(format!("{} face tags for {} triangles", tri_face.len(), triangles.len())));
        }
        if triangles.is_empty() {
            return Err(Error::Invalid("the surface has no triangles".into()));
        }
        check_input(positions, triangles)?;

        let welded = weld::weld(positions, triangles, options.weld_tolerance, options.weld_distance)?;
        if welded.triangles.is_empty() {
            return Err(Error::Invalid("every triangle collapses when the vertices are welded".into()));
        }
        let pts: Vec<Point3> = welded.points.iter().map(|&p| Point3::try_from(p)).collect::<Result<_>>()?;
        let tri_face: Vec<u32> = welded.source_triangle.iter().map(|&t| tri_face.get(t as usize).copied().unwrap_or(0)).collect();
        let surface =
            Surface { pts: &pts, triangles: &welded.triangles, source_vertex: &welded.source_vertex, source_triangle: &welded.source_triangle };

        let mut defects = Defects::new();
        for (t, tri) in surface.triangles.iter().enumerate() {
            if intersect::is_degenerate(surface.corners(tri)?) {
                defects.push(Defect::Degenerate { triangle: surface.source_tri(t), at: surface.centroid(t) });
            }
        }
        defects.into_result()?;

        let topo = topology::Topology::build(&surface)?;
        bvh::self_intersections(&surface)?.into_result()?;
        let shells = topology::shells(&surface, &topo)?;
        let feature_edges = topology::feature_edges(&surface, &topo, &tri_face, options.feature_angle);
        let weld::Welded { points, triangles, source_vertex, source_triangle } = welded;
        Ok(BoundarySurface {
            points,
            triangles,
            tri_face,
            neighbors: topo.neighbors,
            tri_shell: shells.tri_shell,
            shells: shells.shells,
            feature_edges,
            source_vertex,
            source_triangle,
        })
    }

    /// Welded vertex positions.
    pub fn points(&self) -> &[[f64; 3]] {
        &self.points
    }

    /// Vertex `v` as a predicate point.
    pub fn point3(&self, v: u32) -> Option<Point3> {
        self.points.get(v as usize).and_then(|&p| Point3::try_from(p).ok())
    }

    /// Triangles over [`Self::points`], counter-clockwise seen from outside the solid.
    pub fn triangles(&self) -> &[[u32; 3]] {
        &self.triangles
    }

    /// The B-rep face of each triangle (the input's tag, 0 when there was none).
    pub fn tri_face(&self) -> &[u32] {
        &self.tri_face
    }

    /// For each triangle `[a, b, c]`, the triangles across its edges `ab`, `bc` and `ca`.
    pub fn neighbors(&self) -> &[[u32; 3]] {
        &self.neighbors
    }

    /// The shell of each triangle (an index into [`Self::shells`]).
    pub fn tri_shell(&self) -> &[u32] {
        &self.tri_shell
    }

    /// Shells in order of their first triangle.
    pub fn shells(&self) -> &[Shell] {
        &self.shells
    }

    /// The enclosed volume: the sum of the shells' signed volumes.
    pub fn volume(&self) -> f64 {
        self.shells.iter().map(|s| s.volume).sum()
    }

    /// Feature edges, sorted by vertices.
    pub fn feature_edges(&self) -> &[FeatureEdge] {
        &self.feature_edges
    }

    /// Feature corners: vertices where the number of feature edges is neither 0 nor 2 (where
    /// three faces meet, or where a feature curve ends), ascending.
    pub fn feature_corners(&self) -> Vec<u32> {
        let mut count = vec![0u32; self.points.len()];
        for e in &self.feature_edges {
            for v in e.vertices {
                if let Some(c) = count.get_mut(v as usize) {
                    *c += 1;
                }
            }
        }
        count.iter().enumerate().filter(|(_, c)| !matches!(c, 0 | 2)).filter_map(|(v, _)| u32::try_from(v).ok()).collect()
    }

    /// The input vertex each point comes from (the first of the vertices welded into it).
    pub fn source_vertex(&self) -> &[u32] {
        &self.source_vertex
    }

    /// The input triangle each triangle comes from (collapsed ones are dropped).
    pub fn source_triangle(&self) -> &[u32] {
        &self.source_triangle
    }
}

fn check_options(o: &SurfaceOptions) -> Result<()> {
    let w = o.weld_tolerance;
    if !(w == 0.0 || (1e-15..=MAX_RELATIVE_WELD).contains(&w)) {
        return Err(Error::Invalid(format!(
            "weld tolerance must be 0 or in [1e-15, {MAX_RELATIVE_WELD}] (relative to the bounding-box diagonal), got {w}"
        )));
    }
    if !(o.weld_distance >= 0.0 && o.weld_distance <= MAX_COORD) {
        return Err(Error::invalid(crate::message!["weld distance must be a length of at least 0, got ", crate::message::length(o.weld_distance)]));
    }
    if !(o.feature_angle > 0.0 && o.feature_angle <= 180.0) {
        return Err(Error::Invalid(format!("feature angle must be in (0, 180] degrees, got {}", o.feature_angle)));
    }
    Ok(())
}

fn check_input(positions: &[[f64; 3]], triangles: &[[u32; 3]]) -> Result<()> {
    let mut defects = Defects::new();
    for (v, p) in positions.iter().enumerate() {
        if p.iter().any(|c| !(c.abs() <= MAX_COORD)) {
            defects.push(Defect::BadCoordinate { vertex: index(v) });
        }
    }
    for (t, tri) in triangles.iter().enumerate() {
        if tri.iter().any(|&v| v as usize >= positions.len()) {
            defects.push(Defect::BadIndex { triangle: index(t) });
        }
    }
    defects.into_result()
}

/// An index known to fit (sizes are capped well below `u32::MAX`).
pub(crate) fn index(i: usize) -> u32 {
    u32::try_from(i).unwrap_or(u32::MAX)
}

/// The welded surface the checks run on.
pub(crate) struct Surface<'a> {
    pts: &'a [Point3],
    triangles: &'a [[u32; 3]],
    source_vertex: &'a [u32],
    source_triangle: &'a [u32],
}

impl Surface<'_> {
    fn point(&self, v: u32) -> Result<&Point3> {
        self.pts.get(v as usize).ok_or_else(|| Error::Failed(format!("surface: no vertex {v}")))
    }

    fn corners(&self, t: &[u32; 3]) -> Result<[&Point3; 3]> {
        Ok([self.point(t[0])?, self.point(t[1])?, self.point(t[2])?])
    }

    fn xyz(&self, v: u32) -> [f64; 3] {
        self.pts.get(v as usize).map_or([0.0; 3], Point3::coords)
    }

    /// The corner positions of triangle `t`.
    fn tri_xyz(&self, t: usize) -> [[f64; 3]; 3] {
        self.triangles.get(t).copied().unwrap_or([0; 3]).map(|v| self.xyz(v))
    }

    fn centroid(&self, t: usize) -> [f64; 3] {
        let [a, b, c] = self.tri_xyz(t);
        [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, (a[2] + b[2] + c[2]) / 3.0]
    }

    fn midpoint(&self, a: u32, b: u32) -> [f64; 3] {
        let (a, b) = (self.xyz(a), self.xyz(b));
        [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, (a[2] + b[2]) / 2.0]
    }

    fn source_tri(&self, t: usize) -> u32 {
        self.source_triangle.get(t).copied().unwrap_or(u32::MAX)
    }

    fn source_vert(&self, v: u32) -> u32 {
        self.source_vertex.get(v as usize).copied().unwrap_or(u32::MAX)
    }

    /// The unnormalized normal of triangle `t` (floating point; for angles only).
    fn normal(&self, t: usize) -> [f64; 3] {
        let [a, b, c] = self.tri_xyz(t);
        cross(sub(b, a), sub(c, a))
    }
}

/// An axis-aligned box (empty until a point is added).
#[derive(Debug, Clone, Copy)]
pub(crate) struct Box3 {
    pub lo: [f64; 3],
    pub hi: [f64; 3],
}

impl Box3 {
    pub const EMPTY: Box3 = Box3 { lo: [f64::INFINITY; 3], hi: [f64::NEG_INFINITY; 3] };

    pub fn add_point(&mut self, p: [f64; 3]) {
        self.add(&Box3 { lo: p, hi: p });
    }

    pub fn add(&mut self, o: &Box3) {
        for ((lo, hi), (olo, ohi)) in self.lo.iter_mut().zip(&mut self.hi).zip(o.lo.iter().zip(&o.hi)) {
            *lo = lo.min(*olo);
            *hi = hi.max(*ohi);
        }
    }

    /// Closed overlap.
    pub fn touches(&self, o: &Box3) -> bool {
        (0..3).all(|k| self.lo[k] <= o.hi[k] && o.lo[k] <= self.hi[k])
    }

    pub fn contains(&self, p: [f64; 3]) -> bool {
        self.touches(&Box3 { lo: p, hi: p })
    }

    pub fn diagonal(&self) -> f64 {
        norm(sub(self.hi, self.lo))
    }
}
