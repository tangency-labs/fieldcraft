//! Tetrahedral mesh quality: Delaunay refinement to a size field, ODT smoothing and aggressive
//! improvement of a constrained Delaunay tetrahedralization ([`crate::cdt`]), then optionally
//! quadratic elements with their mid-edge nodes on the geometry.
//!
//! # Pipeline ([`quality_mesh`], `docs/meshing/literature.md` Part D Step 3)
//! 1. **Boundary entities.** The CDT's boundary faces carry their B-rep face; an edge between two
//!    faces of different B-rep faces, or whose faces turn by more than the feature angle within a
//!    face whose exact surface is not known, is a feature edge (curve). Vertices are corners
//!    (fixed), curve vertices (move along their curve), face vertices (move on their face) or
//!    interior.
//! 2. **Size field** ([`geometry::GradedSize`]): the global size, limited near curved boundary
//!    to `θ r` (θ radians per element along a curve of radius `r`, from the exact surface or
//!    estimated from the tessellation) and in thin walls to their thickness
//!    ([`geometry::Tessellation::thickness`]), graded at a bounded rate.
//! 3. **Coarsening** ([`refine::coarsen`]): edge collapses where the boundary recovery's
//!    Steiner points and the tessellation are denser than the size field.
//! 4. **Refinement interleaved with smoothing** (Tournois et al. 2009): rounds of constrained
//!    Delaunay refinement ([`refine`], Shewchuk 1998, Si 2010) and ODT smoothing ([`smooth`],
//!    Chen–Xu 2004, Alliez et al. 2005), boundary vertices kept on their faces by projection
//!    ([`geometry`]): face vertices on the tessellation, curve vertices on the exact curves from
//!    their insertion on (Boivin–Ollivier-Gooch 2002; Cheng–Dey–Shewchuk 2012).
//! 5. **Aggressive improvement** ([`improve`], Klingner–Shewchuk 2007): smoothing, edge
//!    removal, multi-face removal and vertex insertion, monotone in the worst dihedral sine.
//! 6. **Gate and report**: every tetrahedron admissible by FieldCraft's one criterion, the volume
//!    equal to what the mesh's own boundary encloses; minimum and maximum dihedral angles and
//!    radius ratios ([`TetQuality`]).
//! 7. **Quadratic elements** ([`quadratic`]): mid-edge nodes, projected onto the B-rep faces and
//!    curves on the boundary; an element that would not be valid keeps straight edges there.
//!
//! Every loop is bounded (insertions, tetrahedra, passes); a stage that stops early says so in
//! the report. Everything is deterministic: fixed orders, no hashing into output, and only
//! correctly rounded arithmetic in decisions.

// The tetrahedron kernels index their four corners and six edges by slot; ranges over them are
// clearer than iterator chains (as in the fem and solver crates).
#![allow(clippy::needless_range_loop)]

pub mod canonical;
pub mod complex;
pub mod geometry;
pub mod improve;
pub mod measure;
pub mod ops;
pub mod quadratic;
pub mod refine;
pub mod smooth;
mod winding;

#[cfg(test)]
mod tests;

use std::collections::BTreeMap;

use fieldcraft_model::{ElementKind, SeedConstraint};

use crate::cdt::TetMesh;
use crate::predicates::{Point3, snap_to_domain};
use crate::quality::{GateOptions, GateReport, Metric, Validity, element_quality, enclosed_volume, validity_gate_elements};
use crate::surface::BoundarySurface;
use crate::{Error, Result};

use complex::{IdMap, NONE, TetComplex, VertexKind, face_toward, sorted2, sorted3};
use geometry::{AnalyticSurface, ExactCurves, GradedSize, Tessellation};
use improve::{ImproveOptions, ImproveStats};
use refine::{Context, RefineOptions, RefineStats};

/// A node is on a seeded edge within this, relative to the size.
const ON_EDGE: f64 = 1e-6; // unit-factors: relative to the size
/// A tetrahedron of quality below this is a sliver.
const SLIVER: f64 = 1e-6; // unit-factors: a normalized quality

/// Element order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Order {
    Linear,
    Quadratic,
}

/// Settings of [`quality_mesh`].
#[derive(Debug, Clone, PartialEq)]
pub struct QualityOptions {
    /// Target edge length; `None`: a twentieth of the bounding-box diagonal.
    pub size: Option<f64>,
    /// Radians per element along curved boundary (0: no curvature refinement).
    pub curvature_angle: f64,
    /// Elements across a thin wall: the size is at most the local thickness over this (0: no
    /// proximity sizing).
    pub proximity: f64,
    /// Largest change of size per unit length.
    pub grading: f64,
    /// Smallest size, as a fraction of the global size.
    pub min_size: f64,
    /// Radius–edge bound of the refinement.
    pub radius_edge: f64,
    /// Interleaved refinement and ODT rounds (0: no refinement).
    pub rounds: usize,
    /// ODT passes per round.
    pub odt_passes: usize,
    /// Run the improvement schedule.
    pub improve: bool,
    pub improve_options: ImproveOptions,
    /// An edge whose two boundary faces turn by more than this many degrees is a feature edge.
    pub feature_angle: f64,
    pub order: Order,
    pub max_tets: usize,
    pub max_points: usize,
    /// Edge seeds: the size at each seeded node is the seed's spacing there, and the size field
    /// grows from it at the grading rate. A seed finer than the global size refines along its
    /// edge; one coarser leaves it to the global size.
    pub seeds: Vec<crate::seeds::EdgeSeed>,
    /// Local sizes (size controls): the size field is at most each one's size in its region and
    /// grows from it at the grading rate.
    pub locals: Vec<crate::size::LocalSize>,
    /// B-rep face tags whose boundary triangles are kept exactly: their vertices stay where they
    /// are and no triangle of theirs is split, flipped or removed. The faces a pyramid
    /// transition shares with the tetrahedra ([`crate::transition`]) are frozen, so the two
    /// meshes stay conforming.
    pub frozen_faces: Vec<u32>,
}

impl Default for QualityOptions {
    fn default() -> Self {
        QualityOptions {
            size: None,
            curvature_angle: std::f64::consts::PI / 8.0,
            proximity: 1.0,
            grading: 0.3,
            min_size: 0.02,
            radius_edge: 2.0,
            rounds: 3,
            odt_passes: 3,
            improve: true,
            improve_options: ImproveOptions::default(),
            feature_angle: 45.0,
            order: Order::Linear,
            max_tets: 5_000_000,
            max_points: 1_000_000,
            seeds: Vec::new(),
            locals: Vec::new(),
            frozen_faces: Vec::new(),
        }
    }
}

/// Quality statistics of a tetrahedral mesh (corner geometry).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct TetQuality {
    pub tets: usize,
    pub min_dihedral: f64,
    pub max_dihedral: f64,
    pub min_radius_ratio: f64,
    pub mean_radius_ratio: f64,
    /// The smallest scaled Jacobian at the worst corner (as `mesh.quality` reports it). The
    /// dihedral range does not bound it: a needle has every dihedral angle near 60° or 90° and a
    /// nearly flat corner (#547).
    pub min_scaled_jacobian: f64,
    /// The mean, over the tetrahedra, of the scaled Jacobian at each one's worst corner (the
    /// values [`Self::min_scaled_jacobian`] is the smallest of).
    pub mean_scaled_jacobian: f64,
    /// Tetrahedra with a dihedral angle below 10° or above 170°.
    pub slivers: usize,
    /// Tetrahedra not admissible by FieldCraft's validity criterion.
    pub inadmissible: usize,
}

/// The quality of tetrahedra given by corner coordinates, with the metrics `mesh.quality`
/// reports ([`crate::quality::element_quality`]).
pub fn tet_quality_stats(points: &[[f64; 3]], tets: &[[u32; 4]]) -> Result<TetQuality> {
    let mut s = TetQuality { min_dihedral: 180.0, min_radius_ratio: 1.0, min_scaled_jacobian: 1.0, ..TetQuality::default() };
    let (mut sum, mut sum_sj) = (0.0, 0.0);
    for t in tets {
        let p = |v: u32| points.get(v as usize).copied().ok_or_else(|| Error::Failed("quality: vertex id".into()));
        let q = element_quality(ElementKind::Tet4, &[p(t[0])?, p(t[1])?, p(t[2])?, p(t[3])?])?;
        let m = |metric: Metric| q.get(metric).unwrap_or(f64::NAN);
        let (lo, hi, rr) = (m(Metric::MinDihedral), m(Metric::MaxDihedral), m(Metric::RadiusRatio));
        s.min_dihedral = s.min_dihedral.min(lo);
        s.max_dihedral = s.max_dihedral.max(hi);
        s.min_radius_ratio = s.min_radius_ratio.min(rr);
        s.min_scaled_jacobian = s.min_scaled_jacobian.min(m(Metric::ScaledJacobian));
        sum += rr;
        sum_sj += m(Metric::ScaledJacobian);
        s.slivers += usize::from(!(lo >= 10.0 && hi <= 170.0));
        s.inadmissible += usize::from(q.validity() != Some(Validity::Valid));
    }
    s.tets = tets.len();
    if tets.is_empty() {
        s.min_dihedral = 0.0;
        s.min_radius_ratio = 0.0;
        s.min_scaled_jacobian = 0.0;
    } else {
        s.mean_radius_ratio = sum / tets.len() as f64;
        s.mean_scaled_jacobian = sum_sj / tets.len() as f64;
    }
    Ok(s)
}

/// A boundary face of the output: tetrahedron, the slot of its opposite vertex, B-rep face.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TaggedFace {
    pub tet: u32,
    pub opposite: u8,
    pub face: u32,
}

impl TaggedFace {
    /// The face's 1-based number in the model's tet4 (and tet10) face order.
    pub fn element_face(&self) -> Option<u8> {
        let k = usize::from(self.opposite);
        let n = ElementKind::Tet4.faces().iter().position(|f| !f.contains(&k))?;
        u8::try_from(n + 1).ok()
    }
}

/// What [`quality_mesh`] did.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct QualityReport {
    /// The constrained Delaunay tetrahedralization it started from.
    pub input: TetQuality,
    pub output: TetQuality,
    pub refine: Vec<RefineStats>,
    pub odt_moves: usize,
    /// Edges collapsed where the input was denser than the size field.
    pub collapsed: usize,
    /// Boundary vertices that could not be moved all the way onto their B-rep face or curve (the
    /// move would have inverted or spoilt a tetrahedron); they stay between the tessellation and
    /// the geometry.
    pub unsnapped: usize,
    /// The largest distance of an unsnapped vertex from its face or curve.
    pub snap_gap: f64,
    pub improve: ImproveStats,
    /// The global size used.
    pub size: f64,
    /// Volume of the mesh, and what its boundary encloses.
    pub volume: f64,
    pub boundary_volume: f64,
    /// Quadratic elements whose curved edges had to stay straight to be valid.
    pub straightened: usize,
    /// Why a stage stopped early, if one did, and the edge seeds whose nodes could not all be
    /// placed.
    pub stopped: Vec<String>,
    /// Per edge seed: its name, the count it asks for and the boundary edges the mesh has along
    /// it.
    pub seeds: Vec<(String, usize, usize)>,
}

/// The quality mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct QualityMesh {
    /// Corner vertices, then (quadratic) mid-edge nodes.
    pub points: Vec<[f64; 3]>,
    /// Elements in the model's node order (4 or 10 nodes).
    pub elements: Vec<Vec<u32>>,
    pub kind: ElementKind,
    pub boundary: Vec<TaggedFace>,
    pub report: QualityReport,
}

fn bounds(points: &[[f64; 3]]) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in points {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    (lo, hi)
}

/// The size field of a boundary: `global` limited by curvature, thin walls, edge seeds and local
/// sizes, graded.
pub fn size_field(surface: &BoundarySurface, geometry: &Tessellation, global: f64, o: &QualityOptions) -> GradedSize {
    let (lo, hi) = bounds(surface.points());
    let floor = o.min_size * global;
    let mut samples = Vec::new();
    if o.curvature_angle > 0.0 {
        let radii = geometry::curvature_radii(surface, geometry);
        for (p, r) in surface.points().iter().zip(radii) {
            let h = (o.curvature_angle * r).max(floor);
            if h < global {
                samples.push((*p, h));
            }
        }
    }
    if o.proximity > 0.0 {
        // Thin walls: elements no larger than the wall's thickness over `proximity`.
        for (p, t) in surface.points().iter().zip(geometry.thickness()) {
            let h = (t / o.proximity).max(floor);
            if h < global {
                samples.push((*p, h));
            }
        }
    }
    for seed in &o.seeds {
        for (p, h) in crate::seeds::node_sizes(seed) {
            let h = h.max(floor);
            if h < global {
                samples.push((p, h));
            }
        }
    }
    GradedSize::new(global, &samples, &o.locals, lo, hi, o.grading)
}

/// The seeds in an order of their own (by name, then polyline), so the nodes [`seed_curves`]
/// inserts, and with them the mesh, do not depend on the order they were given in (#894).
fn seed_order(seeds: &[crate::seeds::EdgeSeed]) -> Vec<&crate::seeds::EdgeSeed> {
    let bits = |s: &crate::seeds::EdgeSeed| s.points.iter().flatten().map(|x| x.to_bits()).collect::<Vec<_>>();
    let mut out: Vec<_> = seeds.iter().collect();
    out.sort_by(|a, b| a.name.cmp(&b.name).then_with(|| bits(a).cmp(&bits(b))));
    out
}

/// Put the seeds' nodes on their feature curves, pinned (corner vertices, which neither
/// smoothing nor coarsening moves), so the curve is divided where the seed says
/// ([`place_seed_node`]). An `exact` seed's curve then keeps only its nodes and is locked
/// ([`hold_exact`]); a `min` or `free` seed's edges may still be split. Returns a line for each
/// seed some of whose nodes could not be placed, so a count that is not held is said.
fn seed_curves(cx: &mut TetComplex, ctx: &Context, seeds: &[&crate::seeds::EdgeSeed], global: f64) -> Result<Vec<String>> {
    let mut short = Vec::new();
    for seed in seeds {
        let tol = 2.0 * seed.tolerance + ON_EDGE * global;
        let on = |x: [f64; 3]| crate::seeds::distance_to(&seed.points, x) <= tol;
        let nodes = crate::seeds::place(&seed.points, &seed.fractions(seed.count, false))?;
        // The curve edges along the seed, kept up to date as nodes split them.
        let mut along = edges_on(cx, &on)?;
        let mut missed = 0;
        for &p in nodes.iter().skip(1).take(nodes.len().saturating_sub(2)) {
            missed += usize::from(!place_seed_node(cx, ctx, p, &mut along, global)?);
        }
        if missed > 0 {
            short.push(format!("{missed} nodes of the edge seed {} could not be placed", seed.describe()));
        }
        if seed.constraint == SeedConstraint::Exact {
            hold_exact(cx, &on, &along)?;
        }
    }
    Ok(short)
}

/// Put seed node `p` on the curve edge of `along` (the seed's curve edges, updated) nearest it.
/// Within [`refine::MIN_GAP`] of the local size of a vertex of that edge it takes that vertex,
/// moved onto the node when it is not a corner already. Otherwise it is inserted with a cavity
/// that leaves out every tetrahedron whose new one would not be admissible
/// ([`ops::insert_delaunay_admissible`]): splitting a long edge beside a short one in line with
/// it, the plain Delaunay cavity makes a flat tetrahedron and the node was lost. Returns whether
/// the node is placed.
fn place_seed_node(cx: &mut TetComplex, ctx: &Context, p: [f64; 3], along: &mut Vec<((u32, u32), u32)>, global: f64) -> Result<bool> {
    let mut best: Option<(f64, ((u32, u32), u32))> = None;
    for &e in along.iter() {
        let d = crate::seeds::distance_to(&[cx.xyz(e.0.0)?, cx.xyz(e.0.1)?], p);
        if best.is_none_or(|x| d < x.0) {
            best = Some((d, e));
        }
    }
    let Some((_, ((a, b), c))) = best else { return Ok(false) };
    let gap = refine::MIN_GAP * ctx.h(p, global);
    let near = [a, b].into_iter().map(|v| Ok((measure::dist(cx.xyz(v)?, p), v))).collect::<Result<Vec<_>>>()?;
    if let Some(&(_, v)) = near.iter().filter(|x| x.0 <= gap).min_by(|x, y| x.0.total_cmp(&y.0)) {
        if cx.kind(v) != VertexKind::Corner {
            if let Some(q) = ctx.place(cx.kind(v), p) {
                cx.try_move(v, q)?;
            }
            cx.set_kind(v, VertexKind::Corner);
        }
        return Ok(true);
    }
    let Some(q) = ctx.place(VertexKind::Curve(c), p) else { return Ok(false) };
    let Some(v) = ops::insert_delaunay_admissible(cx, q, VertexKind::Curve(c), ops::Site::Edge(a, b), gap)? else {
        return Ok(false);
    };
    cx.set_kind(v, VertexKind::Corner);
    along.retain(|&(e, _)| e != (a, b));
    along.extend([(sorted2(a, v), c), (sorted2(v, b), c)]);
    Ok(true)
}

/// The feature edges (with their curve) with both ends where `on` holds.
fn edges_on(cx: &TetComplex, on: &impl Fn([f64; 3]) -> bool) -> Result<Vec<((u32, u32), u32)>> {
    let mut out = Vec::new();
    for (e, c) in cx.curve_edges() {
        if on(cx.xyz(e.0)?) && on(cx.xyz(e.1)?) {
            out.push((e, c));
        }
    }
    Ok(out)
}

/// Make an exact seed's nodes the only vertices of its curve (where `on` holds) and keep them so
/// (#894): every curve vertex there that is not one of them (not a corner) is collapsed onto a
/// neighbour along its curve, and the curve's edges are locked
/// ([`TetComplex::lock_curve_edges`]). Refinement and improvement would otherwise split them
/// where a finer seed, or the curvature, beside the seed grades the size field down.
fn hold_exact(cx: &mut TetComplex, on: &impl Fn([f64; 3]) -> bool, along: &[((u32, u32), u32)]) -> Result<()> {
    let mut extra: Vec<u32> = along.iter().flat_map(|&((a, b), _)| [a, b]).collect();
    extra.retain(|&v| matches!(cx.kind(v), VertexKind::Curve(_)));
    extra.sort_unstable();
    extra.dedup();
    let mut collapsed = false;
    for v in extra {
        // `neighbours` is sorted, so the choice is deterministic.
        for u in cx.neighbours(v)? {
            if cx.curve(v, u).is_some() && ops::collapse(cx, v, u, f64::INFINITY, SLIVER)? {
                collapsed = true;
                break;
            }
        }
    }
    // A collapse changes the curve's edges: find them again.
    let edges: Vec<(u32, u32)> = if collapsed { edges_on(cx, on)? } else { along.to_vec() }.into_iter().map(|(e, _)| e).collect();
    cx.lock_curve_edges(&edges);
    Ok(())
}

/// How many boundary edges of a mesh lie along each seed's edge: both ends on its polyline, and
/// the middle within a fifth of the edge's length of it (an element edge on a curved edge cuts
/// its chord).
fn seed_counts(points: &[[f64; 3]], tets: &[[u32; 4]], boundary: &[TaggedFace], seeds: &[crate::seeds::EdgeSeed]) -> Vec<(String, usize, usize)> {
    if seeds.is_empty() {
        return Vec::new();
    }
    let mut edges: Vec<(u32, u32)> = Vec::with_capacity(boundary.len() * 3);
    for b in boundary {
        let Some(t) = tets.get(b.tet as usize) else { continue };
        let f = face_toward(*t, usize::from(b.opposite));
        for (a, c) in [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])] {
            edges.push(sorted2(a, c));
        }
    }
    edges.sort_unstable();
    edges.dedup();
    let (lo, hi) = bounds(points);
    let eps = ON_EDGE * measure::dist(lo, hi);
    seeds
        .iter()
        .map(|seed| {
            let tol = 2.0 * seed.tolerance + eps;
            let on = |p: [f64; 3], slack: f64| crate::seeds::distance_to(&seed.points, p) <= tol + slack;
            let n = edges
                .iter()
                .filter(|(a, c)| {
                    let (Some(&pa), Some(&pc)) = (points.get(*a as usize), points.get(*c as usize)) else { return false };
                    let m = [(pa[0] + pc[0]) / 2.0, (pa[1] + pc[1]) / 2.0, (pa[2] + pc[2]) / 2.0];
                    on(pa, 0.0) && on(pc, 0.0) && on(m, 0.2 * measure::dist(pa, pc))
                })
                .count();
            (seed.name.clone(), seed.count, n)
        })
        .collect()
}

/// The CDT as a [`TetComplex`] with B-rep face tags, feature curves and vertex kinds; and the
/// feature curves' face pairs. A face with an exact surface in `analytic` (by tag) has no
/// feature curves inside it ([`refine::feature_curves`]).
pub fn complex_from_cdt(
    surface: &BoundarySurface,
    cdt: &TetMesh,
    analytic: &[Option<AnalyticSurface>],
    o: &QualityOptions,
) -> Result<(TetComplex, Vec<[u32; 2]>)> {
    let tri_face = surface.tri_face();
    let tags: Vec<(u32, usize, u32)> =
        cdt.boundary.iter().map(|b| (b.tet, usize::from(b.opposite), tri_face.get(b.triangle as usize).copied().unwrap_or(NONE))).collect();
    let kinds = vec![VertexKind::Interior; cdt.points.len()];
    let mut cx =
        TetComplex::with_invalid(&cdt.points, kinds, &cdt.tets, &tags, &[], o.max_points.max(cdt.points.len()), o.max_tets.max(cdt.tets.len()))?;
    let feature_cos = geometry::sin_cos(o.feature_angle.clamp(0.0, 180.0).to_radians()).1;
    let smooth: Vec<bool> = analytic.iter().map(Option::is_some).collect();
    let (table, edges) = refine::feature_curves(&cx, feature_cos, &smooth)?;
    for &((a, b), c) in &edges {
        cx.set_curve(a, b, Some(c));
    }
    // Vertex kinds from the boundary faces and curve edges at each vertex.
    let mut face_tag: Vec<Option<u32>> = vec![None; cdt.points.len()];
    let mut mixed = vec![false; cdt.points.len()];
    let mut frozen = vec![false; cdt.points.len()];
    cx.freeze_faces(&o.frozen_faces);
    for (t, vs) in cx.live_tets() {
        let (ns, tg) = (cx.adj(t)?, cx.tags(t)?);
        for k in 0..4 {
            if ns[k] != NONE {
                continue;
            }
            for v in face_toward(vs, k) {
                let i = v as usize;
                if o.frozen_faces.contains(&tg[k])
                    && let Some(f) = frozen.get_mut(i)
                {
                    *f = true;
                }
                match face_tag.get(i).copied().flatten() {
                    None => {
                        if let Some(s) = face_tag.get_mut(i) {
                            *s = Some(tg[k]);
                        }
                    }
                    Some(x) if x != tg[k] => {
                        if let Some(m) = mixed.get_mut(i) {
                            *m = true;
                        }
                    }
                    _ => {}
                }
            }
        }
    }
    let mut curve_count: IdMap<u32, Vec<u32>> = IdMap::default();
    for &((a, b), c) in &edges {
        curve_count.entry(a).or_default().push(c);
        curve_count.entry(b).or_default().push(c);
    }
    for v in 0..cdt.points.len() as u32 {
        let kind = match (face_tag.get(v as usize).copied().flatten(), curve_count.get(&v).map(Vec::as_slice)) {
            _ if frozen.get(v as usize).copied().unwrap_or(false) => VertexKind::Corner,
            (None, _) => VertexKind::Interior,
            (Some(tag), None) if !mixed.get(v as usize).copied().unwrap_or(true) => VertexKind::Face(tag),
            (Some(_), Some(&[c1, c2])) if c1 == c2 => VertexKind::Curve(c1),
            _ => VertexKind::Corner,
        };
        cx.set_kind(v, kind);
    }
    Ok((cx, table))
}

/// Most sweeps of [`snap_to_geometry`] over the boundary vertices.
const SNAP_SWEEPS: usize = 16;
/// The fractions of its way to the geometry a vertex tries when the whole way is refused.
const SNAP_STEPS: [f64; 3] = [0.5, 0.25, 0.125];
/// A vertex this close to its projection, relative to its coordinates, is on its geometry
/// (projecting again only moves it by rounding). Scaled by `1 + |x|`, as the projections' own
/// round-off test in `geometry.rs` is: relative away from the origin, with a floor of this much
/// model length near it, far below any part's tolerances.
const ON_GEOMETRY: f64 = 1e-12; // unit-factors: relative to the coordinates' size

/// What [`snap_to_geometry`] left off the geometry.
#[derive(Debug, Clone, Default, PartialEq)]
struct Snapped {
    /// Vertices still off their B-rep face or curve.
    off: usize,
    /// The largest distance of one from its face or curve.
    gap: f64,
    /// Why the sweeps stopped while vertices were still moving, if they did.
    stopped: Option<String>,
}

/// Move every face and curve vertex onto its B-rep face or curve where that keeps the
/// tetrahedra around it positive (and, when `guarded`, their worst quality at least the smaller
/// of what it was and sin 17.5°, so a finished mesh is not spoilt).
///
/// A vertex whose whole move is refused goes part of the way (half, a quarter, an eighth), and
/// the sweeps repeat while anything moves: where refinement made elements smaller than the
/// tessellation's chordal gap, the vertices of a facet are all off the surface by about as much,
/// and no one of them can move alone without spoiling its star, but the layer moves out
/// together, a step at a time (#674). Returns the vertices left off their geometry, and the
/// largest distance of one from it.
fn snap_to_geometry(cx: &mut TetComplex, ctx: &Context, guarded: bool) -> Result<Snapped> {
    // Only the vertices a sweep did not put on their geometry go into the next one.
    let mut off: Vec<u32> =
        (0..cx.point_count() as u32).filter(|&v| cx.is_used(v) && matches!(cx.kind(v), VertexKind::Face(_) | VertexKind::Curve(_))).collect();
    let mut capped = true;
    for sweep in 0..SNAP_SWEEPS {
        if snap_sweep(cx, ctx, guarded, &mut off, sweep == 0)? == 0 {
            capped = false;
            break;
        }
    }
    let mut left = Snapped::default();
    for &v in &off {
        if let Some((_, gap, on)) = off_by(cx, ctx, v)?
            && gap > on
        {
            left.off += 1;
            left.gap = left.gap.max(gap);
        }
    }
    if capped && left.off > 0 {
        left.stopped = Some(format!("snapping to the geometry stopped after {SNAP_SWEEPS} sweeps with {} vertices still off it", left.off));
    }
    Ok(left)
}

/// Vertex `v`'s point on its face or curve, its distance from it, and the distance within which
/// it is on it ([`ON_GEOMETRY`]).
fn off_by(cx: &TetComplex, ctx: &Context, v: u32) -> Result<Option<(Point3, f64, f64)>> {
    let x = cx.xyz(v)?;
    Ok(ctx.place(cx.kind(v), x).map(|p| (p, measure::dist(p.coords(), x), ON_GEOMETRY * (1.0 + measure::norm(x)))))
}

/// One sweep of [`snap_to_geometry`] over the vertices `off` its geometry: moves each that is
/// (the `first` sweep also by rounding) all the way or part of it, and keeps in `off` those it
/// did not put on it. Returns how many moved.
fn snap_sweep(cx: &mut TetComplex, ctx: &Context, guarded: bool, off: &mut Vec<u32>, first: bool) -> Result<usize> {
    let (mut moved, mut still) = (0, Vec::new());
    for &v in off.iter() {
        let Some((p, gap, on)) = off_by(cx, ctx, v)? else { continue };
        if gap <= if first { 0.0 } else { on } {
            continue;
        }
        let x = cx.xyz(v)?;
        let star = cx.star(v)?;
        let floor = if guarded { ops::worst(cx, &star)?.min(smooth::ODT_FLOOR) } else { f64::NEG_INFINITY };
        let partial = SNAP_STEPS.iter().filter_map(|&s| Point3::try_from(measure::lerp(x, p.coords(), s).map(snap_to_domain)).ok());
        let mut whole = false;
        for (i, q) in std::iter::once(p).chain(partial).enumerate() {
            if smooth::star_quality(cx, v, q.coords(), &star)? >= floor && cx.try_move(v, q)? {
                (moved, whole) = (moved + 1, i == 0);
                break;
            }
        }
        if !whole {
            still.push(v);
        }
    }
    *off = still;
    Ok(moved)
}

/// What refinement, smoothing and improvement place boundary vertices on (`plain`: the
/// tessellation alone; `exact`: with the faces' exact surfaces).
///
/// Curve vertices go on the exact curves from the moment they are inserted, as Delaunay
/// refinement on a curved boundary splits a boundary curve at a point of the curve itself
/// (Boivin and Ollivier-Gooch 2002; Cheng, Dey and Shewchuk 2012 on piecewise-smooth
/// complexes). Placed on the tessellation's chord instead, a curve vertex inserted below the
/// facet size sits up to the chord's sag off its curve, the faces beside it are refined around
/// it there, and the final snap has no room left to move it (#674).
///
/// Face vertices stay on the tessellation, a piecewise-flat surface whose creases they can mesh
/// around; the finished mesh then moves them onto the exact surfaces ([`snap_to_geometry`]),
/// and tet10 puts its mid-edge nodes there. Projected onto the exact surfaces while improving,
/// face vertices on both sides of a blend's tangent edge (where two faces meet without a
/// crease) come out too close to flat to separate.
fn refinement_geometry<'a>(plain: &'a Tessellation, exact: &'a Tessellation) -> ExactCurves<'a> {
    ExactCurves { faces: plain, curves: exact }
}

/// The volume enclosed by the boundary faces of a mesh ([`enclosed_volume`]).
fn boundary_volume(points: &[[f64; 3]], tets: &[[u32; 4]], boundary: &[TaggedFace]) -> Result<f64> {
    let mut triangles = Vec::with_capacity(boundary.len());
    for b in boundary {
        let t = tets.get(b.tet as usize).ok_or_else(|| Error::Failed("boundary face of no tetrahedron".into()))?;
        // face_toward points into the tetrahedron; outward is its reverse.
        let f = face_toward(*t, usize::from(b.opposite));
        triangles.push([f[0] as usize, f[2] as usize, f[1] as usize]);
    }
    enclosed_volume(points, &triangles)
}

use fieldcraft_model::progress;

/// Run the quality pipeline on a constrained Delaunay tetrahedralization of `surface` (see the
/// module docs). `analytic[face]` is the exact surface of a B-rep face when known.
pub fn quality_mesh(surface: &BoundarySurface, cdt: &TetMesh, analytic: Vec<Option<AnalyticSurface>>, o: &QualityOptions) -> Result<QualityMesh> {
    let mut report = QualityReport { input: tet_quality_stats(&cdt.points, &cdt.tets)?, ..QualityReport::default() };
    let (lo, hi) = bounds(surface.points());
    let global = match o.size {
        Some(h) if h > 0.0 && h.is_finite() => h,
        Some(h) => return Err(Error::invalid(crate::message!["size ", crate::message::length(h), " must be a positive number"])),
        None => measure::dist(lo, hi) / 20.0,
    };
    if !(global > 0.0 && global.is_finite()) {
        return Err(Error::Invalid("the boundary has no extent".into()));
    }
    report.size = global;
    let (mut cx, curve_faces) = complex_from_cdt(surface, cdt, &analytic, o)?;
    let geometry = Tessellation::new(surface, analytic);
    let size = size_field(surface, &geometry, global, o);
    let plain = Tessellation::new(surface, Vec::new());
    let refining = refinement_geometry(&plain, &geometry);
    let ctx = Context { size: &size, geometry: &refining, curve_faces: &curve_faces };
    let exact = Context { size: &size, geometry: &geometry, curve_faces: &curve_faces };
    // Remove what rounding left flat first, so refinement starts from a valid mesh.
    let flat =
        |cx: &TetComplex| -> Vec<(u32, [u32; 4])> { cx.live_tets().filter(|(_, vs)| ops::tet_quality(cx, *vs).is_ok_and(|q| q < SLIVER)).collect() };
    let bad: Vec<u32> = flat(&cx).into_iter().map(|(t, _)| t).collect();
    if !bad.is_empty() {
        improve::topological_pass(&mut cx, &bad)?;
        for (_, vs) in flat(&cx) {
            for v in vs {
                smooth::smooth_vertex(&mut cx, &ctx, v)?;
            }
        }
    }
    // From here on, no change or move lets the skin pass through itself.
    cx.index_skin(global)?;
    // Coarsen what the boundary recovery made denser than the size field asks for.
    report.collapsed = refine::coarsen(&mut cx, &ctx, 0.6, 1.4, 0.05, 6)?;
    report.stopped.extend(seed_curves(&mut cx, &ctx, &seed_order(&o.seeds), global)?);
    let ropt = RefineOptions { radius_edge: o.radius_edge, ..RefineOptions::default() };
    progress::stage("Refinement round", Some(o.rounds as f64), "");
    for round in 0..o.rounds {
        progress::count((round + 1) as f64);
        progress::check()?;
        let stats = refine::refine(&mut cx, &ctx, &ropt)?;
        let quiet = stats.interior + stats.boundary == 0;
        let stopped = stats.stopped.clone();
        report.refine.push(stats);
        for _ in 0..o.odt_passes {
            report.odt_moves += smooth::odt_pass(&mut cx, &ctx)?;
        }
        if let Some(s) = stopped {
            report.stopped.push(s);
            break;
        }
        if quiet {
            break;
        }
    }
    if o.improve {
        report.improve = improve::improve(&mut cx, &ctx, &o.improve_options)?;
        if let Some(s) = &report.improve.stopped {
            report.stopped.push(s.clone());
        }
    }
    // Put the boundary vertices on the exact surfaces (the tessellation is within its chordal
    // tolerance of them).
    progress::check()?;
    progress::stage("Snapping to the geometry", None, "");
    let snapped = snap_to_geometry(&mut cx, &exact, true)?;
    (report.unsnapped, report.snap_gap) = (snapped.off, snapped.gap);
    report.stopped.extend(snapped.stopped);
    let c = cx.compact()?;
    let boundary: Vec<TaggedFace> = c.boundary.iter().map(|&(tet, k, face)| TaggedFace { tet, opposite: k as u8, face }).collect();
    report.output = tet_quality_stats(&c.points, &c.tets)?;
    report.seeds = seed_counts(&c.points, &c.tets, &boundary, &o.seeds);
    report.boundary_volume = boundary_volume(&c.points, &c.tets, &boundary)?;
    let gate = linear_gate(&c.points, &c.tets, &boundary, report.boundary_volume)?;
    report.volume = gate.volume;
    gate.into_result()?;
    match o.order {
        Order::Linear => {
            Ok(QualityMesh { elements: c.tets.iter().map(|t| t.to_vec()).collect(), points: c.points, kind: ElementKind::Tet4, boundary, report })
        }
        Order::Quadratic => {
            let curve_of: IdMap<(u32, u32), u32> = c.curves.iter().map(|&(e, k)| (sorted2(e.0, e.1), k)).collect();
            let q = quadratic::tet10(&c.points, &c.tets, &boundary, &curve_of, &curve_faces, &geometry)?;
            report.straightened = q.straightened;
            Ok(QualityMesh { points: q.points, elements: q.elements, kind: ElementKind::Tet10, boundary, report })
        }
    }
}

/// How many tetrahedra hold each face.
fn face_use(tets: &[[u32; 4]]) -> BTreeMap<[u32; 3], u32> {
    let mut count = BTreeMap::new();
    for t in tets {
        for k in 0..4 {
            *count.entry(sorted3(face_toward(*t, k))).or_default() += 1;
        }
    }
    count
}

/// The tetrahedra (1-based) that share none of their faces with another one, in a mesh of more
/// than one: each is a separate piece of boundary, joined to the mesh by vertices or an edge.
fn detached_tets(tets: &[[u32; 4]], count: &BTreeMap<[u32; 3], u32>) -> Vec<u32> {
    if tets.len() < 2 {
        return Vec::new();
    }
    tets.iter()
        .enumerate()
        .filter(|(_, t)| (0..4).all(|k| count.get(&sorted3(face_toward(**t, k))).copied().unwrap_or(0) < 2))
        .map(|(i, _)| i as u32 + 1)
        .collect()
}

/// Where the mesh is joined to itself by less than a face (#671): the boundary edges that are not
/// in exactly two boundary faces, and the vertices round which the tetrahedra fall into groups
/// joined by no face holding the vertex, or the boundary faces into more than one fan. A group
/// of tetrahedra that share faces with each other but meet the rest only at an edge or a vertex
/// shows here, where [`detached_tets`] sees only a single tetrahedron. The boundary surface the
/// mesher starts from is a manifold ([`BoundarySurface`] refuses an edge or a vertex that joins
/// separate sheets), so a valid mesh has none.
fn pinches(points: usize, tets: &[[u32; 4]]) -> (Vec<[u32; 2]>, Vec<u32>) {
    // Every face with its tetrahedron, grouped by its sorted vertices.
    let mut faces: Vec<([u32; 3], u32)> = Vec::with_capacity(tets.len() * 4);
    for (t, vs) in tets.iter().enumerate() {
        for k in 0..4 {
            faces.push((sorted3(face_toward(*vs, k)), t as u32));
        }
    }
    faces.sort_unstable();
    // Union-find over (tetrahedron, corner) slots, t * 4 + j, joined through the faces that hold
    // the corner's vertex; and over (boundary face, corner) slots, joined through boundary edges.
    let slot = |t: u32, v: u32| tets.get(t as usize).and_then(|c| c.iter().position(|&x| x == v)).map(|j| t as usize * 4 + j);
    let mut star = Sets::new(tets.len() * 4);
    let mut skin: Vec<[u32; 3]> = Vec::new();
    for run in faces.chunk_by(|a, b| a.0 == b.0) {
        match run {
            [(f, _)] => skin.push(*f),
            [(f, t1), (_, t2)] => {
                for v in *f {
                    if let (Some(a), Some(b)) = (slot(*t1, v), slot(*t2, v)) {
                        star.join(a, b);
                    }
                }
            }
            // A face held by more than two: the fold check reports it.
            _ => {}
        }
    }
    let mut edges: Vec<([u32; 2], usize)> = Vec::with_capacity(skin.len() * 3);
    for (i, f) in skin.iter().enumerate() {
        for (a, b) in [(f[0], f[1]), (f[1], f[2]), (f[0], f[2])] {
            edges.push(([a, b], i));
        }
    }
    edges.sort_unstable();
    let mut fans = Sets::new(skin.len() * 3);
    let corner = |i: usize, v: u32| skin.get(i).and_then(|f| f.iter().position(|&x| x == v)).map(|j| i * 3 + j);
    let mut bad_edges = Vec::new();
    for run in edges.chunk_by(|a, b| a.0 == b.0) {
        match run {
            [(e, i), (_, j)] => {
                for v in *e {
                    if let (Some(a), Some(b)) = (corner(*i, v), corner(*j, v)) {
                        fans.join(a, b);
                    }
                }
            }
            [(e, _), ..] => bad_edges.push(*e),
            _ => {}
        }
    }
    // A vertex is pinched when its slots have more than one root, in either structure.
    let mut pinched = vec![false; points];
    let star_roots = tets.iter().enumerate().flat_map(|(t, vs)| vs.iter().enumerate().map(move |(j, &v)| (v, t * 4 + j)));
    mark_split(&mut pinched, star_roots.map(|(v, s)| (v, star.root(s))));
    let fan_roots = skin.iter().enumerate().flat_map(|(i, f)| f.iter().enumerate().map(move |(j, &v)| (v, i * 3 + j)));
    mark_split(&mut pinched, fan_roots.map(|(v, s)| (v, fans.root(s))));
    let bad_vertices = pinched.iter().enumerate().filter(|(_, x)| **x).map(|(v, _)| v as u32).collect();
    (bad_edges, bad_vertices)
}

/// Mark the vertices whose (vertex, root) pairs name more than one root.
fn mark_split(pinched: &mut [bool], roots: impl Iterator<Item = (u32, usize)>) {
    let mut first = vec![usize::MAX; pinched.len()];
    for (v, r) in roots {
        if let (Some(f), Some(x)) = (first.get_mut(v as usize), pinched.get_mut(v as usize)) {
            *x |= *f != usize::MAX && *f != r;
            *f = r;
        }
    }
}

/// A union-find over `n` items (path halving; the smaller root wins, so roots are deterministic).
struct Sets(Vec<usize>);

impl Sets {
    fn new(n: usize) -> Sets {
        Sets((0..n).collect())
    }

    fn root(&mut self, mut i: usize) -> usize {
        while let Some(&p) = self.0.get(i) {
            if p == i {
                break;
            }
            let grand = self.0.get(p).copied().unwrap_or(p);
            if let Some(slot) = self.0.get_mut(i) {
                *slot = grand;
            }
            i = grand;
        }
        i
    }

    fn join(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.root(a), self.root(b));
        if ra != rb
            && let Some(slot) = self.0.get_mut(ra.max(rb))
        {
            *slot = ra.min(rb);
        }
    }
}

/// The constrained Delaunay tetrahedralization itself as a [`QualityMesh`] (no quality stage),
/// through the same validity gate: every tetrahedron admissible, the volume equal to what its
/// boundary encloses.
pub fn raw_mesh(surface: &BoundarySurface, cdt: &TetMesh) -> Result<QualityMesh> {
    let tri_face = surface.tri_face();
    let boundary: Vec<TaggedFace> = cdt
        .boundary
        .iter()
        .map(|b| TaggedFace { tet: b.tet, opposite: b.opposite, face: tri_face.get(b.triangle as usize).copied().unwrap_or(NONE) })
        .collect();
    let stats = tet_quality_stats(&cdt.points, &cdt.tets)?;
    let boundary_volume = boundary_volume(&cdt.points, &cdt.tets, &boundary)?;
    let gate = linear_gate(&cdt.points, &cdt.tets, &boundary, boundary_volume)?.into_result()?;
    let report = QualityReport { input: stats, output: stats, volume: gate.volume, boundary_volume, ..QualityReport::default() };
    Ok(QualityMesh { points: cdt.points.clone(), elements: cdt.tets.iter().map(|t| t.to_vec()).collect(), kind: ElementKind::Tet4, boundary, report })
}

/// The validity gate of a linear mesh: every tetrahedron admissible, joined to the others by a
/// face (when there are several), no face held by more than two, no part of the mesh joined to
/// the rest by an edge or a vertex alone ([`pinches`]), the faces listed as `boundary`
/// exactly those held by one, and the volume equal to what the boundary encloses within 1e-9.
/// The tetrahedra round each interior edge wind about it exactly once, and less than once round
/// each boundary edge ([`winding`], exact): this catches a star that overlaps itself round an
/// edge, not every intersection between unrelated parts of a mesh.
pub fn linear_gate(points: &[[f64; 3]], tets: &[[u32; 4]], boundary: &[TaggedFace], boundary_volume: f64) -> Result<GateReport> {
    let count = face_use(tets);
    // A mesh folded over itself: the volume check misses it when the fold is consistently
    // oriented, its two layers adding up to what the boundary encloses.
    let folded: Vec<&[u32; 3]> = count.iter().filter(|(_, n)| **n > 2).map(|(f, _)| f).collect();
    if let Some(first) = folded.first() {
        let at = points.get(first[0] as usize).copied().unwrap_or_default();
        return Err(Error::failed(crate::message![
            format!("mesh validity gate: {} faces are held by more than two tetrahedra (the first near ", folded.len()),
            crate::message::point(at),
            "): the mesh overlaps itself there"
        ]));
    }
    // A boundary face that two tetrahedra hold: the mesh closed onto its own boundary there (two
    // sides of a narrow gap meeting in one face), which neither the volume nor the fold sees.
    let mut listed = Vec::with_capacity(boundary.len());
    for b in boundary {
        let t = tets.get(b.tet as usize).ok_or_else(|| Error::Failed("gate: boundary face of no tetrahedron".into()))?;
        listed.push(sorted3(face_toward(*t, usize::from(b.opposite))));
    }
    listed.sort_unstable();
    listed.dedup();
    let closed: Vec<&[u32; 3]> = listed.iter().filter(|f| count.get(*f).copied().unwrap_or(0) != 1).collect();
    if let Some(first) = closed.first() {
        let at = points.get(first[0] as usize).copied().unwrap_or_default();
        return Err(Error::failed(crate::message![
            format!("mesh validity gate: {} boundary faces are held by two tetrahedra (the first near ", closed.len()),
            crate::message::point(at),
            "): the mesh closes onto its own boundary there"
        ]));
    }
    let once = count.values().filter(|n| **n == 1).count();
    if listed.len() != boundary.len() || once != listed.len() {
        return Err(Error::Failed(format!(
            "mesh validity gate: the mesh lists {} boundary faces where {once} faces are held by one tetrahedron",
            boundary.len()
        )));
    }
    let detached = detached_tets(tets, &count);
    drop(count);
    if !detached.is_empty() {
        let list = detached.iter().take(10).map(u32::to_string).collect::<Vec<_>>().join(", ");
        return Err(Error::Failed(format!(
            "mesh validity gate: {} tetrahedra share no face with the rest of the mesh ({list}): held by vertices alone, they would make the stiffness singular",
            detached.len()
        )));
    }
    let (edges, vertices) = pinches(points.len(), tets);
    if let Some(&[a, b]) = edges.first() {
        let (pa, pb) = (points.get(a as usize).copied().unwrap_or_default(), points.get(b as usize).copied().unwrap_or_default());
        return Err(Error::failed(crate::message![
            format!("mesh validity gate: {} boundary edges are not in exactly two boundary faces (the first, {a}–{b}, near ", edges.len()),
            crate::message::point([(pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0, (pa[2] + pb[2]) / 2.0]),
            "): the mesh is joined to itself by an edge alone there, which would make the stiffness singular"
        ]));
    }
    if let Some(&v) = vertices.first() {
        return Err(Error::failed(crate::message![
            format!("mesh validity gate: the mesh is pinched at {} vertices (the first, {v}, near ", vertices.len()),
            crate::message::point(points.get(v as usize).copied().unwrap_or_default()),
            "): tetrahedra joined to the rest by a vertex or an edge alone would make the stiffness singular"
        ]));
    }
    let coords = |t: &[u32; 4]| -> Result<[[f64; 3]; 4]> {
        let p = |v: u32| points.get(v as usize).copied().ok_or_else(|| Error::Failed("gate: vertex id".into()));
        Ok([p(t[0])?, p(t[1])?, p(t[2])?, p(t[3])?])
    };
    let elements = tets.iter().enumerate().map(|(i, t)| Ok((i as u32 + 1, ElementKind::Tet4, coords(t)?))).collect::<Result<Vec<_>>>()?;
    let report = validity_gate_elements(
        elements,
        &GateOptions { expected_volume: Some(boundary_volume), volume_tolerance: crate::quality::VOLUME_TOLERANCE },
    )?;
    // The winding count needs every tetrahedron positive: a report that failed is returned as it
    // is, and the caller turns it into an error.
    if report.passed() {
        winding::check(points, tets)?;
    }
    Ok(report)
}
