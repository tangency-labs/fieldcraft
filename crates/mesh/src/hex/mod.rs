//! Structured all-hexahedral meshing of sweepable and mappable bodies.
//!
//! [`hex_mesh`] meshes one body's boundary surface (face-tagged, from surface preparation) with
//! hexahedra only, by the first method that fits:
//! 1. **Extrude** ([`layered`]): a straight multi-level sweep: every face is across the direction
//!    (a source or target level) or along it (a linking wall), possibly at several levels with
//!    regions of different heights (stepped blocks, ribs on plates, shelled boxes). The cross-
//!    section is imprinted, quad-meshed region by region and stacked layer by layer.
//! 2. **Sweep** ([`sweep`]): a one-to-one sweep between a planar source and a planar target face
//!    through rings of logically rectangular linking faces (bent tubes, lofts, partial
//!    revolutions): linking faces get transfinite grids projected onto their surfaces, and every
//!    layer's interior nodes come from the source mesh by the least-squares affine map of the
//!    boundary rings (Roca and Sarrate) with boundary residual correction (Knupp).
//! 3. **Revolve** ([`layered`] in a rotating frame): bodies of revolution whose cross-section does
//!    not reach the axis (rings, tori, grooved pulleys).
//! 4. **Mapped** ([`mapped`]): a logical cube (six four-sided faces, three edges at each of eight
//!    corners), meshed by transfinite interpolation of its face grids (Gordon and Thiel).
//!
//! A body none of these fits is [`Error::Unsupported`] with each method's reason: the input for
//! decomposition into sweepable and mappable regions.
//!
//! **Partitioned parts.** [`hex_mesh_cells`] ([`cells`]) meshes a part split into cells that
//! meet face to face (SolidCraft split bodies) as one conforming mesh: every cell a mapped block,
//! interval assignment over all cells together, each shared face meshed once and its nodes
//! merged by identity. [`hex_mesh_swept_cells`] sweeps such cells together instead: one
//! cross-section imprinted with every cell's walls, one interval assignment and one set of
//! layers, so cells that share a sweep conform across their cut faces.
//!
//! **Common steps.** Interval assignment ([`interval`]) picks the division counts of shared curves
//! together. Boundary nodes on analytic faces are placed on the exact surface; quadratic hexahedra
//! (hex20) get their mid-edge nodes projected onto the faces (and onto both faces of a B-rep
//! edge). Every boundary quadrilateral is tagged with the B-rep face it lies on, so model sets
//! named by face resolve. The output passes the certified validity gate (every element's
//! Jacobian determinant proven positive, [`crate::quality`]), a tiling check (the elements'
//! volume equals the volume their boundary encloses, so nothing overlaps), and a fidelity check
//! (the volume differs from the boundary surface's by no more than the gap between the mesh
//! boundary and the surface accounts for).

mod arrangement;
mod brep;
mod cells;
mod decompose;
mod extrude_to;
pub(crate) mod geom;
mod interval;
mod layered;
mod mapped;
mod region;
mod section;
mod sweep;

#[cfg(test)]
mod tests;

use std::collections::HashMap;

pub use cells::{Cell, CellsMesh, GROWTH, MAX_CELLS, hex_mesh_cells, hex_mesh_swept_cells};
pub use extrude_to::{ExtrudeBoundary, ExtrudeLayers, ExtrudePath, ExtrudeRole, ExtrudeSpec, Extruded, MAX_BIAS, MAX_LAYERS, extrude_to_surface};
pub use geom::Surface;

use crate::seeds::EdgeSeed;

use brep::Brep;
use geom::{TriangleTree, add, dist, scale, unit};

use crate::message;
use crate::message::Message;
use crate::quad::QuadOptions;
use crate::quality::{GateOptions, GateReport, Metric, element_quality, validity_gate_checks};
use crate::size::SizeField;
use crate::surface::BoundarySurface;
use crate::vec3::{P, cross, dot, sub};
use crate::{Error, Result};
use fieldcraft_model::ElementKind;

/// A projection step this small, relative to the tolerance, has converged.
const CONVERGED: f64 = 1e-3; // unit-factors: relative to the tolerance

/// Which structured method to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum HexMethod {
    /// The first method that fits, in the order extrude, sweep, revolve, mapped.
    Auto,
    /// The first sweeping method that fits: extrude, sweep, revolve.
    Swept,
    Extrude,
    Sweep,
    Revolve,
    Mapped,
}

impl HexMethod {
    pub fn name(self) -> &'static str {
        match self {
            HexMethod::Auto => "auto",
            HexMethod::Swept => "swept",
            HexMethod::Extrude => "extrude",
            HexMethod::Sweep => "sweep",
            HexMethod::Revolve => "revolve",
            HexMethod::Mapped => "mapped",
        }
    }

    pub fn from_name(name: &str) -> Option<HexMethod> {
        [HexMethod::Auto, HexMethod::Swept, HexMethod::Extrude, HexMethod::Sweep, HexMethod::Revolve, HexMethod::Mapped]
            .into_iter()
            .find(|m| m.name() == name)
    }
}

/// Options of [`hex_mesh`].
#[derive(Debug, Clone, PartialEq)]
pub struct HexOptions {
    pub method: HexMethod,
    /// Hex20 (mid-edge nodes on the geometry) instead of hex8.
    pub quadratic: bool,
    /// The chordal tolerance the boundary surface was tessellated with (0 when exact): copies of a
    /// curve at different levels closer than this are one curve, and the fidelity check allows
    /// it.
    pub chordal_tolerance: f64,
    /// The most hexahedra ([`Error::Limit`] beyond).
    pub max_elements: usize,
    /// Options of the face meshers.
    pub quad: QuadOptions,
    /// Edge seeds: the B-rep edges they lie on get their counts (exact seeds are held, and a
    /// method that cannot keep them names them) and their spacing.
    pub seeds: Vec<EdgeSeed>,
}

impl Default for HexOptions {
    fn default() -> Self {
        HexOptions {
            method: HexMethod::Auto,
            quadratic: false,
            chordal_tolerance: 0.0,
            max_elements: 5_000_000,
            quad: QuadOptions::default(),
            seeds: Vec::new(),
        }
    }
}

/// A boundary face of the mesh and the B-rep face it lies on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BoundaryFace {
    /// Element index (0-based, into [`HexMesh::element`]).
    pub element: u32,
    /// Face number (1-based, as in [`ElementKind::faces`]).
    pub face: u8,
    /// The B-rep face tag (`tri_face`).
    pub brep_face: u32,
}

/// How the mesh was made and how good it is.
#[derive(Debug, Clone, PartialEq)]
pub struct HexReport {
    pub method: HexMethod,
    /// What the method found (direction, levels, regions, layers, …).
    pub description: Message,
    /// Smallest and mean scaled Jacobian (Knupp) over the elements.
    pub min_scaled_jacobian: f64,
    pub mean_scaled_jacobian: f64,
    /// The elements' total volume, and the volume the boundary surface encloses.
    pub volume: f64,
    pub surface_volume: f64,
    /// The validity gate's report.
    pub gate: GateReport,
}

/// An all-hexahedral mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct HexMesh {
    pub nodes: Vec<[f64; 3]>,
    /// [`ElementKind::Hex8`] or [`ElementKind::Hex20`].
    pub kind: ElementKind,
    /// Element node indices, `kind.node_count()` per element, in the model's node order.
    pub connectivity: Vec<u32>,
    pub boundary: Vec<BoundaryFace>,
    pub report: HexReport,
}

impl HexMesh {
    pub fn element_count(&self) -> usize {
        self.connectivity.len() / self.kind.node_count()
    }

    /// The node indices of element `e`.
    pub fn element(&self, e: usize) -> Option<&[u32]> {
        let n = self.kind.node_count();
        self.connectivity.get(e * n..(e + 1) * n)
    }
}

/// What every method gets.
pub(crate) struct Ctx<'a> {
    pub brep: &'a Brep<'a>,
    pub tree: &'a TriangleTree<'a>,
    pub field: &'a dyn SizeField,
    pub opts: &'a HexOptions,
    /// Classification tolerance (model units).
    pub tol: f64,
    /// The seed on each B-rep edge (index into `opts.seeds`, and whether the edge runs from the
    /// seed's end to its start).
    pub edge_seeds: Vec<Option<(usize, bool)>>,
}

impl<'a> Ctx<'a> {
    pub(crate) fn new(brep: &'a Brep<'a>, tree: &'a TriangleTree<'a>, field: &'a dyn SizeField, opts: &'a HexOptions) -> Result<Ctx<'a>> {
        Ctx::with_tolerance(brep, tree, field, opts, RELATIVE_TOLERANCE * brep.diagonal)
    }

    /// A context with the classification tolerance `tol` (bodies meshed together share one).
    pub(crate) fn with_tolerance(
        brep: &'a Brep<'a>,
        tree: &'a TriangleTree<'a>,
        field: &'a dyn SizeField,
        opts: &'a HexOptions,
        tol: f64,
    ) -> Result<Ctx<'a>> {
        let curves: Vec<Vec<P>> = (0..brep.edges.len()).map(|e| brep.edge_points(e, true)).collect();
        let edge_seeds = crate::seeds::match_curves(&opts.seeds, &curves, seed_tolerance(opts, tol))?;
        Ok(Ctx { brep, tree, field, opts, tol, edge_seeds })
    }

    /// The seed on B-rep edge `e`, and whether it runs against the edge.
    pub(crate) fn seed_of(&self, e: usize) -> Option<(&EdgeSeed, bool)> {
        let (i, rev) = self.edge_seeds.get(e).copied().flatten()?;
        Some((self.opts.seeds.get(i)?, rev))
    }

    /// A division-count variable for B-rep edge `e` with goal `goal`, its seed applied.
    pub(crate) fn edge_var(&self, ia: &mut interval::Intervals, e: usize, goal: f64) -> usize {
        let v = ia.var(goal, 1);
        if let Some((seed, _)) = self.seed_of(e) {
            ia.seed(v, seed);
        }
        v
    }

    /// `n + 1` points along B-rep edge `e` (traversed forwards or not): at the seed's spacing
    /// when it has one, or the spacing `carried` to it from a biased seed on an edge it is tied
    /// to ([`Ctx::carried`]), else evenly in the size field.
    pub(crate) fn divide_edge(&self, e: usize, forward: bool, n: usize, carried: &[Option<(usize, bool)>]) -> Result<Vec<P>> {
        let pts = self.brep.edge_points(e, forward);
        let carried = carried.get(e).copied().flatten().and_then(|(i, rev)| Some((self.opts.seeds.get(i)?, rev)));
        match self.seed_of(e).or(carried) {
            Some((seed, rev)) => crate::seeds::place(&pts, &seed.fractions(n, rev == forward)),
            None => crate::size::divide_polyline(self.field, &pts, n),
        }
    }

    /// The seeds carried along `ties` (uses of B-rep edges divided alike) from biased seeds to
    /// unseeded edges ([`crate::seeds::carry`]), by index into the seeds.
    pub(crate) fn carried(&self, ties: &[[crate::seeds::Use; 2]]) -> Vec<Option<(usize, bool)>> {
        let mut own = self.edge_seeds.clone();
        own.resize(self.brep.edges.len(), None);
        crate::seeds::carry(&own, |&i| self.opts.seeds.get(i).is_some_and(|s| s.bias.is_some()), ties)
    }
}

/// How far apart a seed's polyline and a mesher's curve of the same edge may lie: both are
/// tessellations of it.
pub(crate) fn seed_tolerance(opts: &HexOptions, tol: f64) -> f64 {
    4.0 * opts.chordal_tolerance.max(0.0) + 16.0 * tol
}

/// What a method produces: linear hexahedra in the model's node order.
pub(crate) struct Raw {
    pub nodes: Vec<P>,
    pub hexes: Vec<[u32; 8]>,
    pub description: Message,
}

/// Refinements of the cross-section curves tried when a face mesher fails (each halves them).
const MAX_REFINE: u32 = 2;

/// Classification tolerance relative to the bounding-box diagonal.
pub(crate) const RELATIVE_TOLERANCE: f64 = 1e-6; // unit-factors: relative to the diagonal

/// Mesh the solid bounded by `surface` with hexahedra. `faces[tag]` is the analytic surface of
/// B-rep face `tag` when known (exact node placement; free-form faces use the triangles);
/// `field` sets the element size. A body no structured method fits is [`Error::Unsupported`]
/// with the reasons.
pub fn hex_mesh(surface: &BoundarySurface, faces: &[Option<Surface>], field: &dyn SizeField, opts: &HexOptions) -> Result<HexMesh> {
    if !(opts.chordal_tolerance >= 0.0 && opts.chordal_tolerance.is_finite()) {
        return Err(Error::invalid(crate::message![
            "the chordal tolerance must be a non-negative number, got ",
            crate::message::length(opts.chordal_tolerance)
        ]));
    }
    if surface.shells().len() != 1 {
        return Err(Error::Unsupported(format!(
            "the boundary has {} shells; structured meshing takes one solid without internal voids",
            surface.shells().len()
        )));
    }
    check_estimate(&[surface], field, opts)?;
    let brep = Brep::new(surface, faces)?;
    let tree = TriangleTree::new(surface)?;
    let ctx = Ctx::new(&brep, &tree, field, opts)?;
    let methods: &[HexMethod] = match opts.method {
        HexMethod::Auto => &[HexMethod::Extrude, HexMethod::Sweep, HexMethod::Revolve, HexMethod::Mapped],
        HexMethod::Swept => &[HexMethod::Extrude, HexMethod::Sweep, HexMethod::Revolve],
        HexMethod::Extrude => &[HexMethod::Extrude],
        HexMethod::Sweep => &[HexMethod::Sweep],
        HexMethod::Revolve => &[HexMethod::Revolve],
        HexMethod::Mapped => &[HexMethod::Mapped],
    };
    let mut reasons: Vec<Message> = Vec::new();
    for &method in methods {
        fieldcraft_model::progress::check()?;
        let mut last = None;
        match method {
            HexMethod::Extrude | HexMethod::Revolve => {
                let setups = match candidate_setups(std::slice::from_ref(&ctx), method) {
                    Ok(s) => s,
                    Err(why) => {
                        reasons.push(message![format!("{}: ", method.name())].then(why));
                        continue;
                    }
                };
                for setup in setups {
                    fieldcraft_model::progress::check()?;
                    // The cross-section is laid out once per frame; refinement retries reuse it.
                    let prepared = match layered::prepare(std::slice::from_ref(&ctx), setup) {
                        Ok(p) => p,
                        Err(e) if e.kind() == Some(crate::Kind::Limit) => return Err(e),
                        Err(e) => {
                            last = Some(e);
                            continue;
                        }
                    };
                    match with_refinement(&ctx, method, |refine| layered::mesh(std::slice::from_ref(&ctx), &prepared, refine).map(|m| m.0)) {
                        Ok(mesh) => return Ok(mesh),
                        Err(e) if e.kind() == Some(crate::Kind::Limit) => return Err(e),
                        Err(e) => last = Some(e),
                    }
                }
            }
            HexMethod::Sweep => match with_refinement(&ctx, method, |refine| sweep::mesh(&ctx, refine)) {
                Ok(mesh) => return Ok(mesh),
                Err(e) if e.kind() == Some(crate::Kind::Limit) => return Err(e),
                Err(e) => last = Some(e),
            },
            _ => match mapped::mesh(&ctx).and_then(|raw| finish(&Cells::one(&ctx), raw, method)) {
                Ok(mesh) => return Ok(mesh),
                Err(e) if e.kind() == Some(crate::Kind::Limit) => return Err(e),
                Err(e) => last = Some(e),
            },
        }
        if let Some(e) = last {
            reasons.push(message![format!("{}: ", method.name())].then(plain(&e)));
        }
    }
    Err(Error::Measured(
        crate::Kind::Unsupported,
        message!["no structured hex method fits this body ("].then(Message::join(reasons, "; ")).then(message![")"]),
    ))
}

/// Refuse a mesh far beyond the element limit before any work: the solids' volume over the cube
/// of the smallest size the field asks for on their boundaries.
fn check_estimate(surfaces: &[&BoundarySurface], field: &dyn SizeField, opts: &HexOptions) -> Result<()> {
    let (mut h_min, mut volume) = (f64::INFINITY, 0.0);
    for s in surfaces {
        for p in s.points() {
            h_min = h_min.min(crate::size::checked_size(field, *p)?);
        }
        volume += s.volume().abs();
    }
    let estimate = volume / h_min.powi(3);
    if !(estimate <= 8.0 * opts.max_elements as f64) {
        return Err(Error::Limit(format!("about {estimate:.3e} hexahedra at this size; the limit is {}", opts.max_elements)));
    }
    Ok(())
}

/// The frames a layered method can try, fewest levels first, or the first reason none fits.
pub(crate) fn candidate_setups(ctxs: &[Ctx], method: HexMethod) -> std::result::Result<Vec<layered::Setup>, Message> {
    let breps: Vec<&Brep> = ctxs.iter().map(|c| c.brep).collect();
    let frames = if method == HexMethod::Extrude { layered::translations(&breps) } else { layered::rotations(&breps) };
    let mut fitting: Vec<(usize, usize, layered::Setup)> = Vec::new();
    let mut first_reason = None;
    for (k, f) in frames.into_iter().enumerate() {
        match layered::setup(ctxs, &f) {
            Ok(s) => fitting.push((s.level_count(), k, s)),
            Err(e) => {
                first_reason.get_or_insert_with(|| plain(&e));
            }
        }
    }
    if fitting.is_empty() {
        return Err(first_reason.unwrap_or_else(|| message![if method == HexMethod::Extrude { "no direction to try" } else { "no axis to try" }]));
    }
    fitting.sort_by_key(|x| (x.0, x.1));
    Ok(fitting.into_iter().take(layered::MAX_CANDIDATES).map(|x| x.2).collect())
}

/// Run a method's mesher and [`finish`], and when that fails (a face mesher on a boundary too
/// coarse for it, or curved elements the gate rejects) again with the cross-section's curves
/// halved, up to [`MAX_REFINE`] times. Limits and unsupported shapes end it at once.
fn with_refinement(ctx: &Ctx, method: HexMethod, run: impl Fn(u32) -> Result<Raw>) -> Result<HexMesh> {
    let mut last = Error::Failed("no attempt".into());
    for refine in 0..=MAX_REFINE {
        fieldcraft_model::progress::check()?;
        match run(refine).and_then(|raw| finish(&Cells::one(ctx), raw, method)) {
            Ok(mesh) => return Ok(mesh),
            Err(e) if e.kind() == Some(crate::Kind::Failed) => last = e,
            Err(e) => return Err(e),
        }
    }
    Err(last)
}

/// An error's message without the "not supported yet" prefix (for the combined reasons).
fn plain(e: &Error) -> Message {
    match e {
        Error::Unsupported(m) => message![m.clone()],
        Error::Measured(crate::Kind::Unsupported, m) => m.clone(),
        other => other.message(),
    }
}

/// Hex face `k` (0-based) corner positions in the element.
fn hex_face(k: usize) -> [usize; 4] {
    ElementKind::Hex8.faces().get(k).and_then(|f| <[usize; 4]>::try_from(*f).ok()).unwrap_or([0, 1, 2, 3])
}

/// The cells a mesh spans, for [`finish`]: one context per cell and the cell of each element.
pub(crate) struct Cells<'c, 'a> {
    pub ctxs: &'c [Ctx<'a>],
    /// The cell of each hexahedron (index into `ctxs`); empty: all in the first.
    pub of: &'c [u32],
    /// The cells' names, for messages.
    pub names: &'c [&'c str],
}

impl<'c, 'a> Cells<'c, 'a> {
    pub(crate) fn one(ctx: &'c Ctx<'a>) -> Self {
        Cells { ctxs: std::slice::from_ref(ctx), of: &[], names: &[] }
    }

    fn name(&self, c: usize) -> String {
        self.names.get(c).map_or_else(|| c.to_string(), |n| format!("`{n}`"))
    }

    fn cell(&self, e: usize) -> usize {
        self.of.get(e).map_or(0, |&c| c as usize)
    }

    fn ctx(&self, e: usize) -> Result<&Ctx<'a>> {
        self.ctxs.get(self.cell(e)).ok_or_else(|| Error::Failed(format!("element {e} has no cell")))
    }
}

/// Buckets of a uniform grid over the cells' bounding boxes: which cells may be near a point.
struct CellGrid {
    lo: P,
    step: f64,
    dims: [usize; 3],
    buckets: Vec<Vec<usize>>,
}

/// Cosine of the largest angle between a boundary face's normal and the reversed normal of
/// another cell's triangle for the face to lie on that triangle (about 25°): a face on another
/// cell faces it, while faces that merely meet at an edge (a fillet beside a split) do not.
const FACING: f64 = 0.9;

/// Buckets per axis of a [`CellGrid`] at most.
const MAX_GRID_DIM: usize = 64;

impl CellGrid {
    /// The cells' boxes grown by `gap`.
    fn new(ctxs: &[Ctx], gap: f64) -> CellGrid {
        let boxes: Vec<(P, P)> = ctxs
            .iter()
            .map(|c| {
                c.brep.surface.points().iter().fold(([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]), |(mut lo, mut hi), p| {
                    for k in 0..3 {
                        lo[k] = lo[k].min(p[k] - gap);
                        hi[k] = hi[k].max(p[k] + gap);
                    }
                    (lo, hi)
                })
            })
            .collect();
        let (lo, hi) = boxes.iter().fold(([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]), |(mut lo, mut hi), (a, b)| {
            for k in 0..3 {
                lo[k] = lo[k].min(a[k]);
                hi[k] = hi[k].max(b[k]);
            }
            (lo, hi)
        });
        let per_axis = (crate::vec3::cbrt(ctxs.len() as f64).ceil() as usize * 2).clamp(1, MAX_GRID_DIM);
        let extent = (0..3).map(|k| hi[k] - lo[k]).fold(0.0, f64::max);
        let step = if extent > 0.0 && extent.is_finite() { extent / per_axis as f64 } else { 1.0 };
        let dim = |k: usize| (((hi[k] - lo[k]) / step).ceil() as usize).clamp(1, MAX_GRID_DIM);
        let dims = [dim(0), dim(1), dim(2)];
        let mut grid = CellGrid { lo, step, dims, buckets: vec![Vec::new(); dims[0] * dims[1] * dims[2]] };
        for (c, (a, b)) in boxes.iter().enumerate() {
            let (i0, i1) = (grid.index(*a), grid.index(*b));
            for i in i0[0]..=i1[0] {
                for j in i0[1]..=i1[1] {
                    for k in i0[2]..=i1[2] {
                        if let Some(bucket) = grid.buckets.get_mut((k * dims[1] + j) * dims[0] + i) {
                            bucket.push(c);
                        }
                    }
                }
            }
        }
        grid
    }

    fn index(&self, p: P) -> [usize; 3] {
        [0, 1, 2].map(|k| {
            let x = ((p[k] - self.lo[k]) / self.step).floor();
            if x.is_finite() && x > 0.0 { (x as usize).min(self.dims[k] - 1) } else { 0 }
        })
    }

    fn near(&self, p: P) -> &[usize] {
        let [i, j, k] = self.index(p);
        self.buckets.get((k * self.dims[1] + j) * self.dims[0] + i).map(Vec::as_slice).unwrap_or(&[])
    }
}

/// A boundary quadrilateral of a mesh being finished: the element face, its cell and the B-rep
/// face (index into the cell's faces) it lies on.
struct Skin {
    element: u32,
    face: u8,
    cell: usize,
    brep: usize,
}

/// Tag, check, upgrade to hex20 when asked, gate and measure. Besides the gate, the tiling and
/// the fidelity checks, no boundary face of a mesh of several cells may lie on another cell's
/// boundary (a free internal face: cells that do not conform there), and the mesh boundary must
/// be a closed 2-manifold of quadrilaterals (every edge on exactly two).
fn finish(cells: &Cells, raw: Raw, method: HexMethod) -> Result<HexMesh> {
    let Raw { mut nodes, hexes, description } = raw;
    let first = cells.ctxs.first().ok_or_else(|| Error::Failed("no cell to finish".into()))?;
    let opts = first.opts;
    if hexes.is_empty() {
        return Err(Error::Failed("the method produced no elements".into()));
    }
    if hexes.len() > opts.max_elements {
        return Err(Error::Limit(format!("the mesh has {} hexahedra; the limit is {}", hexes.len(), opts.max_elements)));
    }
    if !cells.of.is_empty() && cells.of.len() != hexes.len() {
        return Err(Error::Failed("a cell is needed for every element".into()));
    }
    fieldcraft_model::check_coordinates(nodes.as_flattened())?;
    let at = |nodes: &[P], i: u32| nodes.get(i as usize).copied().ok_or_else(|| Error::Failed(format!("node {i} out of range")));
    // Boundary faces: hex faces used once.
    let mut uses: HashMap<[u32; 4], u32> = HashMap::with_capacity(hexes.len() * 3);
    for h in &hexes {
        for k in 0..6 {
            let mut key = hex_face(k).map(|i| h.get(i).copied().unwrap_or(0));
            key.sort_unstable();
            *uses.entry(key).or_insert(0) += 1;
        }
    }
    if let Some(n) = uses.values().copied().max().filter(|&n| n > 2) {
        return Err(Error::Failed(format!("a hexahedron face is shared by {n} elements (the mesh is not a manifold)")));
    }
    // A boundary face this close to another cell's boundary lies on it.
    let gap = 2.0 * opts.chordal_tolerance + cells.ctxs.iter().map(|c| c.tol).fold(0.0, f64::max);
    let others = (cells.ctxs.len() > 1).then(|| CellGrid::new(cells.ctxs, gap));
    let mut skin: Vec<Skin> = Vec::new();
    let mut boundary = Vec::new();
    let mut fidelity = 0.0;
    // Boundary edges of the quadrilaterals: each must be on exactly two of them.
    let mut rim: HashMap<(u32, u32), u32> = HashMap::new();
    for (e, h) in hexes.iter().enumerate() {
        for k in 0..6 {
            let corners = hex_face(k).map(|i| h.get(i).copied().unwrap_or(0));
            let mut key = corners;
            key.sort_unstable();
            if uses.get(&key) != Some(&1) {
                continue;
            }
            for m in 0..4 {
                let (a, b) = (corners[m], corners[(m + 1) % 4]);
                *rim.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
            let ctx = cells.ctx(e)?;
            let [a, b, c, d] = [at(&nodes, corners[0])?, at(&nodes, corners[1])?, at(&nodes, corners[2])?, at(&nodes, corners[3])?];
            let centre = scale(add(add(a, b), add(c, d)), 0.25);
            let n = cross(sub(c, a), sub(d, b));
            let area = 0.5 * dot(n, n).sqrt();
            let Some(n) = unit(n) else { return Err(Error::Failed(format!("element {e} has a degenerate boundary face"))) };
            let hit = ctx
                .tree
                .closest(centre, &|t| dot(ctx.tree.normal(t), n) > 0.0)
                .ok_or_else(|| Error::Failed(format!("no boundary face faces the same way as element {e}'s face {}", k + 1)))?;
            let f = ctx.brep.face_of_triangle(hit.triangle).ok_or_else(|| Error::Failed("triangle without a face".into()))?;
            let cell = cells.cell(e);
            if let Some(grid) = &others {
                for &o in grid.near(centre) {
                    let Some(other) = cells.ctxs.get(o).filter(|_| o != cell) else { continue };
                    if other.tree.closest_within(centre, gap, &|t| dot(other.tree.normal(t), n) < -FACING).is_some() {
                        return Err(Error::Failed(format!(
                            "element {e}'s face {} lies on cell {}'s boundary without a neighbour there: cells {} and {} do not conform (split them so that they meet face to face)",
                            k + 1,
                            cells.name(o),
                            cells.name(cell),
                            cells.name(o)
                        )));
                    }
                }
            }
            // The geometry's own face tag (a member of a composite face keeps its own).
            let tag = ctx.brep.surface.tri_face().get(hit.triangle as usize).copied().unwrap_or_else(|| ctx.brep.face_tag(f));
            fidelity += area * (hit.dist + opts.chordal_tolerance);
            boundary.push(BoundaryFace { element: e as u32, face: k as u8 + 1, brep_face: tag });
            skin.push(Skin { element: e as u32, face: k as u8 + 1, cell, brep: f });
        }
    }
    if let Some((&(a, b), &n)) = rim.iter().filter(|(_, n)| **n != 2).min_by_key(|(k, _)| **k) {
        return Err(Error::Failed(format!("the mesh boundary is not closed: the edge between nodes {a} and {b} is on {n} boundary faces")));
    }
    let kind = if opts.quadratic { ElementKind::Hex20 } else { ElementKind::Hex8 };
    let mut connectivity: Vec<u32> = Vec::with_capacity(hexes.len() * kind.node_count());
    if opts.quadratic {
        let mids = mid_edge_nodes(cells, &mut nodes, &hexes, &skin)?;
        for (h, m) in hexes.iter().zip(&mids) {
            connectivity.extend_from_slice(h);
            connectivity.extend_from_slice(m);
        }
    } else {
        for h in &hexes {
            connectivity.extend_from_slice(h);
        }
    }
    let npe = kind.node_count();
    let coords = |e: usize| -> Result<Vec<P>> { connectivity.get(e * npe..(e + 1) * npe).unwrap_or(&[]).iter().map(|&i| at(&nodes, i)).collect() };
    // The volume the mesh's own boundary encloses: equal to the elements' total exactly when they
    // tile it without overlapping.
    let mut enclosed = 0.0;
    for b in &boundary {
        let x = coords(b.element as usize)?;
        let local = kind.faces().get(b.face as usize - 1).copied().unwrap_or(&[]);
        let face: Vec<P> = local.iter().filter_map(|&i| x.get(i).copied()).collect();
        enclosed += face_flux(&face);
    }
    // One certified pass: every element's quality, with its Jacobian check for the gate.
    let (mut min_sj, mut sum_sj) = (f64::INFINITY, 0.0);
    let mut checks = Vec::with_capacity(hexes.len());
    for e in 0..hexes.len() {
        let q = element_quality(kind, &coords(e)?)?;
        let sj = q.get(Metric::ScaledJacobian).unwrap_or(-1.0);
        min_sj = min_sj.min(sj);
        sum_sj += sj;
        checks.push((e as u32 + 1, q.jacobian.ok_or_else(|| Error::Failed("an element was not certified".into()))?));
    }
    let gate = validity_gate_checks(
        checks.iter().map(|(l, c)| (*l, c)),
        &GateOptions { expected_volume: Some(enclosed), volume_tolerance: crate::quality::VOLUME_TOLERANCE },
    )?
    .into_result()
    .map_err(|e| Error::failed(plain(&e).then(message![" ("]).then(description.clone()).then(message![")"])))?;
    let surface_volume: f64 = cells.ctxs.iter().map(|c| c.brep.surface.volume()).sum();
    let slack = 2.0 * fidelity + crate::quality::VOLUME_TOLERANCE * surface_volume.abs();
    if (gate.volume - surface_volume).abs() > slack {
        return Err(Error::failed(
            crate::message![
                "the mesh volume ",
                crate::message::volume(gate.volume),
                " differs from the boundary's ",
                crate::message::volume(surface_volume),
                " by more than the gap between them allows (",
                crate::message::volume(slack),
                "): the mesh does not fill the body ("
            ]
            .then(description)
            .then(message![")"]),
        ));
    }
    let report = HexReport {
        method,
        description,
        min_scaled_jacobian: min_sj,
        mean_scaled_jacobian: sum_sj / hexes.len().max(1) as f64,
        volume: gate.volume,
        surface_volume,
        gate,
    };
    Ok(HexMesh { nodes, kind, connectivity, boundary, report })
}

/// `∫ x n_x dA` over a quad4 or quad8 face (outward normal): its contribution to the enclosed
/// volume. The 3 × 3 Gauss rule integrates it exactly (degree at most 5 per direction).
pub(crate) fn face_flux(face: &[P]) -> f64 {
    const G: [(f64, f64); 3] = [(-0.774_596_669_241_483_4, 5.0 / 9.0), (0.0, 8.0 / 9.0), (0.774_596_669_241_483_4, 5.0 / 9.0)];
    let mut sum = 0.0;
    for (r, wr) in G {
        for (s, ws) in G {
            let [n, dr, ds] = quad_shape(face.len() == 8, r, s);
            let (mut x, mut xr, mut xs) = ([0.0; 3], [0.0; 3], [0.0; 3]);
            for (k, p) in face.iter().enumerate().take(8) {
                x = add(x, scale(*p, n[k]));
                xr = add(xr, scale(*p, dr[k]));
                xs = add(xs, scale(*p, ds[k]));
            }
            sum += wr * ws * x[0] * cross(xr, xs)[0];
        }
    }
    sum
}

/// Shape functions and their `r`, `s` derivatives of quad4 (the first four entries) or quad8 at
/// `(r, s)`.
fn quad_shape(quadratic: bool, r: f64, s: f64) -> [[f64; 8]; 3] {
    const C: [(f64, f64); 4] = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];
    let mut out = [[0.0; 8]; 3];
    for (k, (ri, si)) in C.iter().enumerate() {
        if quadratic {
            out[0][k] = 0.25 * (1.0 + ri * r) * (1.0 + si * s) * (ri * r + si * s - 1.0);
            out[1][k] = 0.25 * ri * (1.0 + si * s) * (2.0 * ri * r + si * s);
            out[2][k] = 0.25 * si * (1.0 + ri * r) * (ri * r + 2.0 * si * s);
        } else {
            out[0][k] = 0.25 * (1.0 + ri * r) * (1.0 + si * s);
            out[1][k] = 0.25 * ri * (1.0 + si * s);
            out[2][k] = 0.25 * si * (1.0 + ri * r);
        }
    }
    if quadratic {
        // Mid-side nodes on sides 0-1 (s = −1), 1-2 (r = 1), 2-3 (s = 1), 3-0 (r = −1).
        out[0][4..].copy_from_slice(&[
            0.5 * (1.0 - r * r) * (1.0 - s),
            0.5 * (1.0 + r) * (1.0 - s * s),
            0.5 * (1.0 - r * r) * (1.0 + s),
            0.5 * (1.0 - r) * (1.0 - s * s),
        ]);
        out[1][4..].copy_from_slice(&[-r * (1.0 - s), 0.5 * (1.0 - s * s), -r * (1.0 + s), -0.5 * (1.0 - s * s)]);
        out[2][4..].copy_from_slice(&[-0.5 * (1.0 - r * r), -s * (1.0 + r), 0.5 * (1.0 - r * r), -s * (1.0 - r)]);
    }
    out
}

/// Projection rounds for a point on several faces (alternating projections onto each).
const PROJECTION_ROUNDS: usize = 24;

/// The closest point to `p` on B-rep face index `f`: its analytic surface, or its triangles.
pub(crate) fn project_to_face(ctx: &Ctx, f: usize, p: P) -> Option<P> {
    let face = ctx.brep.faces.get(f)?;
    if let Some(q) = face.surface.and_then(|s| s.project(p)) {
        return Some(q);
    }
    let tag = face.id;
    ctx.tree.closest(p, &|t| ctx.brep.tag(t) == Some(tag)).map(|h| h.point)
}

/// The point on every listed face (their common curve or corner) nearest `p`, by alternating
/// projections.
pub(crate) fn project_to_faces(ctx: &Ctx, faces: &[usize], p: P) -> P {
    let pairs: Vec<(&Ctx, usize)> = faces.iter().map(|&f| (ctx, f)).collect();
    project_alternating(&pairs, ctx.tol, p)
}

/// The point on every listed face of the listed cells nearest `p`, by alternating projections
/// until a round moves it less than a thousandth of `tol`.
fn project_alternating(faces: &[(&Ctx, usize)], tol: f64, p: P) -> P {
    let mut x = p;
    for _ in 0..PROJECTION_ROUNDS {
        let before = x;
        for &(ctx, f) in faces {
            if let Some(q) = project_to_face(ctx, f, x) {
                x = q;
            }
        }
        if faces.len() < 2 || dist(before, x) <= tol * CONVERGED {
            break;
        }
    }
    x
}

/// Mid-edge nodes of every hexahedron (Hex20 order: positions 8..20), shared between elements
/// (of any cell); boundary edges get theirs on the B-rep faces of the boundary quads around them.
fn mid_edge_nodes(cells: &Cells, nodes: &mut Vec<P>, hexes: &[[u32; 8]], skin: &[Skin]) -> Result<Vec<[u32; 12]>> {
    // B-rep faces (cell, face) at each boundary edge.
    let mut edge_faces: HashMap<(u32, u32), Vec<(usize, usize)>> = HashMap::new();
    for b in skin {
        let Some(h) = hexes.get(b.element as usize) else { continue };
        let c = hex_face(b.face as usize - 1).map(|i| h.get(i).copied().unwrap_or(0));
        for k in 0..4 {
            let (a, z) = (c[k], c[(k + 1) % 4]);
            let list = edge_faces.entry((a.min(z), a.max(z))).or_default();
            if !list.contains(&(b.cell, b.brep)) {
                list.push((b.cell, b.brep));
            }
        }
    }
    let edges = ElementKind::Hex20.edges();
    let tol = cells.ctxs.iter().map(|c| c.tol).fold(f64::INFINITY, f64::min);
    let mut mid: HashMap<(u32, u32), u32> = HashMap::new();
    let mut out = Vec::with_capacity(hexes.len());
    for h in hexes {
        let mut m = [0u32; 12];
        for (k, e) in edges.iter().enumerate() {
            let (Some(&i), Some(&j)) = (e.first().and_then(|&x| h.get(x)), e.get(1).and_then(|&x| h.get(x))) else {
                return Err(Error::Failed("hex20 edge table".into()));
            };
            let key = (i.min(j), i.max(j));
            let id = match mid.get(&key) {
                Some(&id) => id,
                None => {
                    let (a, b) = (nodes.get(i as usize).copied().unwrap_or([0.0; 3]), nodes.get(j as usize).copied().unwrap_or([0.0; 3]));
                    let p = scale(add(a, b), 0.5);
                    let p = match edge_faces.get(&key) {
                        Some(fs) => {
                            let pairs: Vec<(&Ctx, usize)> = fs.iter().filter_map(|&(c, f)| cells.ctxs.get(c).map(|ctx| (ctx, f))).collect();
                            project_alternating(&pairs, tol, p)
                        }
                        None => p,
                    };
                    nodes.push(p);
                    let id = nodes.len() as u32 - 1;
                    mid.insert(key, id);
                    id
                }
            };
            if let Some(slot) = m.get_mut(k) {
                *slot = id;
            }
        }
        out.push(m);
    }
    fieldcraft_model::check_coordinates(nodes.as_flattened())?;
    Ok(out)
}
