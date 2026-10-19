//! Unstructured all-hexahedral meshing of any closed boundary: an adaptive grid cut into
//! conforming hexahedra, a core of grid elements inside the part, two boundary layers that carry
//! the part's faces, sharp curves and corners, and an optimisation that fits them to the
//! geometry while keeping every element valid.
//!
//! # Pipeline
//! 1. **Frame.** The grid is aligned with the part: the normals of its two largest planar
//!    patches when they are orthogonal and cover more area than the coordinate axes do, the
//!    coordinate axes otherwise (octree methods are orientation dependent; N. Pietroni et al.,
//!    *Hex-mesh generation and processing: a survey*, ACM TOG 42(2), 2022, §4.6).
//! 2. **Grid** ([`grid`]). Root cubes over the part's box, refined as a 27-tree where a B-rep
//!    face curves (its normals turn by more than [`OctreeOptions::curvature_angle`] inside a
//!    cell), where the part is thin (fewer than [`OctreeOptions::thickness_cells`] cells across a
//!    wall, measured along the inward normal of each boundary triangle to a face that faces it)
//!    and where a patch is narrow, then strongly balanced and closed under the transition
//!    templates. The levels the thinnest wall and narrowest patch need are found before the grid
//!    is built, and the root size is chosen with them (see [`octree_hex_mesh`]'s plan): a part
//!    that needs more levels or cells than allowed fails at once, with the numbers.
//! 3. **Conforming hexahedra** ([`templates`]): the 3-refinement node-marking templates
//!    (R. Schneiders, 1996; Y. Ito, A. M. Shih, B. K. Soni, IJNME 77(13), 2009), so every grid
//!    leaf becomes 1, 4, 11, 22 or 27 hexahedra and the grid is conforming by construction.
//! 4. **Core.** Elements whose eight corners are inside the part (exact ray parity on the
//!    boundary triangles) and at least [`OctreeOptions::core_offset`] of their cell size from it.
//!    Elements are dropped until the core's boundary is a 2-manifold and every boundary vertex can
//!    move outward (its quads' normals lie in one half-space: the buffer rule of H. Tong,
//!    E. Halilaj, Y. J. Zhang, *HybridOctree_Hex*, arXiv 2401.05984, 2024); its components and
//!    genus must match the part's shells.
//! 5. **Feature map** ([`map`]): the part's patches, sharp curves and corners mapped onto the
//!    core's boundary quads (X. Gao, H. Shen, D. Panozzo, *Feature preserving octree-based
//!    hexahedral meshing*, CGF 38(5), 2019, §3.3), patches first.
//! 6. **Padding and pillowing** ([`topo`]): one layer of hexahedra on the core's boundary (so no
//!    element keeps two faces on the boundary), then each patch's layer is pillowed, which inserts
//!    a ring of quads along every sharp curve on both sides: a curve then runs between quads that
//!    each have one edge on it, so straight and curved features are captured without degenerate
//!    quads (Gao et al. 2019, §3.4, Fig. 9).
//! 7. **Fit** ([`smooth`]): the foldover-free energy of V. Garanzha et al., *Foldover-free maps
//!    in 50 lines of code*, ACM TOG 40(4), 2021, plus a fitting term that draws every boundary
//!    vertex onto its corner, curve or patch with a growing weight (Gao et al. 2019, §3.5); the
//!    barrier keeps every element valid while the boundary moves. Each boundary vertex then goes
//!    exactly onto its feature where its elements stay valid; where they would not (a blade, a
//!    near-tangent contact) the mesh stays off by [`OctreeReport::boundary_distance`].
//! 8. **Gate.** The worst elements are improved by compass search, then every element is proven
//!    valid by the certified Jacobian check (Bernstein bounds, [`crate::quality`]; corner
//!    Jacobians alone do not prove a trilinear hexahedron valid, J. Marschner et al., CGF 39(5),
//!    2020); elements that are not are repaired locally first. A failure in 4–8 names where the
//!    grid is too coarse; the leaves there are split once more and the mesh is rebuilt, at most
//!    [`OctreeOptions::attempts`] times, and a rebuild may not grow the grid more than
//!    [`RETRY_GROWTH`] times (Gao et al. 2019, §3.3).
//!
//! # Why 3-refinement templates
//! The dual schemes of M. Livesu, L. Pitzalis, G. Cherchi, *Optimal dual schemes for adaptive
//! grid based hexmeshing* (ACM TOG 41(2), 2022) give a simpler singular structure, but their
//! installation depends on a paired grid whose transition vertices the reference pipeline
//! (L. Pitzalis et al., *Generalized adaptive refinement for grid-based hexahedral meshing*, ACM
//! TOG 40(6), 2021) selects with an integer linear program, and on heuristics that are retried
//! with other parameters when they fail. The node-marking templates need only the balance and
//! closure rules and are proven conforming by an exhaustive test. Their limit is the one
//! Schneiders states: they close convex refined regions only, so a leaf next to a concave corner
//! of a finer region is split, and the split spreads to the region's orthogonal hull (H. Tong,
//! Y. J. Zhang, *MCHex*, arXiv 2511.02064, 2025, meet the same limit with the same templates).
//! Curvature therefore refines only [`OctreeOptions::curvature_levels`] levels; deeper levels
//! follow thickness and narrow faces, whose refined regions are mostly convex. On curved parts
//! this closure is most of the grid: the corpus's rounded box goes from 8 266 refined cells to
//! 85 460 after closure, the loft from 6 650 to 80 438, the vase from 4 206 to 53 346 (thin walls
//! cost little: the I-beam's 469 304 cells are all refinement). Schneiders' modified edge rule
//! (an edge with a marked end split in three; R. Schneiders, VKI lecture series 2000-4, §5.3)
//! makes every concave face a full 3 × 3 face and the seven-marked cell the 27-split, but the
//! cell with six marks (a concave edge) then needs a template whose two faces at the unmarked
//! edge both fold their three columns into that edge, and no valid construction of it with
//! usable element quality is published; the dual schemes need a 2-refinement grid. Both remain
//! the upgrade path.
//!
//! # Determinism and bounds
//! The same boundary and options give the same mesh: grid cells are processed in sorted order,
//! maps are only looked up, the optimiser sweeps vertices in index order. Grid cells, elements,
//! refinement depth, closure rounds, optimiser sweeps and attempts are capped.

mod geom;
mod grid;
mod map;
mod smooth;
mod templates;
mod topo;

#[cfg(test)]
mod tests;

use std::cell::{Cell, RefCell};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use fieldcraft_model::ElementKind;

pub use crate::hex::BoundaryFace;
use crate::hex::face_flux;
use crate::message;
use crate::message::{Message, length};
use crate::quality::{GateOptions, Metric, Validity, element_quality, validity_gate_checks};
use crate::size::{SizeField, checked_size};
use crate::surface::BoundarySurface;
use crate::vec3::{P, cross, dot, norm, sub};
use crate::{Error, Result};

use geom::{Aabb, Geometry};
use grid::Tree;
use smooth::{Constraint, Optimizer};

/// Options and caps of [`octree_hex_mesh`].
#[derive(Debug, Clone, PartialEq)]
pub struct OctreeOptions {
    /// The most refinement levels below the root cells (each one a third of the size), at most
    /// 8. The levels the part's walls and narrow faces need are found before the grid is built
    /// (see [`octree_hex_mesh`]); a part that needs more fails at once.
    pub max_depth: u8,
    /// The most grid cells and the most output elements.
    pub max_elements: usize,
    /// A cell is split where the normals of one B-rep face inside it turn by more than this
    /// (degrees).
    pub curvature_angle: f64,
    /// Levels below the root that curvature may refine (see the module docs).
    pub curvature_levels: u8,
    /// A cell is split where the part is less than this many cells thick.
    pub thickness_cells: f64,
    /// Core elements keep this fraction of their cell size from the boundary.
    pub core_offset: f64,
    /// Optimiser sweeps over all vertices.
    pub max_sweeps: usize,
    /// Builds, counting the rebuilds after refining where a stage failed (1 = no retry).
    pub attempts: usize,
    /// Boundary edges whose faces meet at more than this angle (degrees between normals) are
    /// sharp features the mesh captures with its own edges; smoother face boundaries are not.
    pub sharp_angle: f64,
    /// Edges whose normals turn by more than this (degrees: thin blades, near-tangent contacts)
    /// are not captured; the mesh rounds them.
    pub blade_angle: f64,
}

impl Default for OctreeOptions {
    fn default() -> Self {
        OctreeOptions {
            max_depth: 4,
            max_elements: 5_000_000,
            curvature_angle: 25.0,
            curvature_levels: 1,
            thickness_cells: 3.0,
            core_offset: 0.3,
            max_sweeps: 300,
            attempts: 4,
            sharp_angle: 30.0,
            blade_angle: 160.0,
        }
    }
}

/// What [`octree_hex_mesh`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct OctreeReport {
    /// Builds run (1 + rebuilds after a failure).
    pub attempts: usize,
    /// Grid leaves of the last build.
    pub grid_cells: usize,
    /// The root cells' edge and the deepest level of the grid (see [`octree_hex_mesh`]).
    pub root_size: f64,
    pub levels: u8,
    /// The thinnest wall and the narrowest patch the grid was chosen for (`∞` when none).
    pub thinnest_wall: f64,
    pub narrowest_patch: f64,
    /// Grid leaves per level (root level first) of the last build: after the refinement
    /// criteria (size, curvature, thickness, narrow patches and retry hints), and after balance
    /// and closure.
    pub refined_cells: Vec<usize>,
    pub closed_cells: Vec<usize>,
    /// Core elements (from the grid) and boundary-layer elements.
    pub core_elements: usize,
    pub layer_elements: usize,
    /// The grid's axes.
    pub frame: [[f64; 3]; 3],
    /// Sharp corners, curves and patches of the boundary.
    pub corners: usize,
    pub curves: usize,
    pub patches: usize,
    /// Optimiser sweeps of the last build.
    pub sweeps: usize,
    /// Smallest and mean corner scaled Jacobian.
    pub min_scaled_jacobian: f64,
    pub mean_scaled_jacobian: f64,
    /// Exact volume of the elements and the boundary's enclosed volume.
    pub volume: f64,
    pub surface_volume: f64,
    /// The largest distance of a boundary vertex from its corner, curve or patch: zero when
    /// every boundary vertex lies on the geometry; where the geometry is sharper than valid
    /// elements can follow (blades, slits), the mesh stays this far off.
    pub boundary_distance: f64,
    /// Why earlier builds were refined.
    pub retries: Vec<Message>,
}

/// An all-hex8 mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct OctreeHexMesh {
    pub nodes: Vec<[f64; 3]>,
    /// Model node order; every element is proven valid.
    pub hexes: Vec<[u32; 8]>,
    pub boundary: Vec<BoundaryFace>,
    pub report: OctreeReport,
}

/// A rebuild after a failure may have at most this many times the grid cells of the build
/// before…
pub const RETRY_GROWTH: usize = 4;
/// … or this many cells, whichever is more.
const RETRY_MIN_CELLS: usize = 250_000;
/// The fitting weights of the successive optimisation stages.
const FIT_STAGES: &[f64] = &[1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0];
/// A boundary vertex goes exactly onto its feature only if its elements keep at least this
/// corner scaled Jacobian.
const MIN_SNAP_QUALITY: f64 = 0.05;
/// Elements below this corner scaled Jacobian are improved by compass search.
const IMPROVE_BELOW: f64 = 0.3;
/// The two faces of a wall open by at most this angle (degrees) from facing each other. Where a
/// convex edge is sharper (a wedge), the part is thinner at a distance `d` from the edge than
/// `2 d tan(30°) ≈ 1.15 d`, so the cells that resolve the edge do not resolve its thickness; at
/// blunter edges they do, and the edge is a feature, not a wall.
const WALL_OPENING: f64 = 60.0;
/// Closure rounds of the grid (each splits at least one leaf; the cell cap bounds them too).
const CLOSURE_ROUNDS: usize = 1000;
/// Rounds of the core clean-up (each drops at least one element).
const CORE_ROUNDS: usize = 1000;

/// A failure that refining at `at` may cure.
struct Retry {
    at: Vec<P>,
    why: Message,
}

enum Attempt {
    Done(Box<OctreeHexMesh>),
    Retry(Retry),
}

/// Mesh the solid bounded by `surface` with hexahedra of about `field`'s size (the root cell
/// size is the field's largest size ([`SizeField::range`]), else the field at the centre of the
/// part's box; cells are refined where the field is smaller than they are; finer cells follow the part's curvature,
/// thickness and narrow faces). See the module docs.
pub fn octree_hex_mesh(surface: &BoundarySurface, field: &dyn SizeField, options: &OctreeOptions) -> Result<OctreeHexMesh> {
    check_options(options)?;
    let g = Geometry::new(surface, options.sharp_angle, options.blade_angle)?;
    let thickness = triangle_thickness(&g);
    let ctx = plan(&g, field, options, &thickness)?;
    let mut hints: Vec<P> = Vec::new();
    let mut retries: Vec<Message> = Vec::new();
    for attempt in 0..options.attempts {
        fieldcraft_model::progress::check()?;
        let cap = if attempt == 0 {
            options.max_elements
        } else {
            options.max_elements.min(RETRY_GROWTH.saturating_mul(ctx.last_leaves.get()).max(RETRY_MIN_CELLS))
        };
        let built = match ctx.build(&hints, cap) {
            Err(e) if attempt > 0 && e.is_limit() => {
                retries.push(message![format!("refining for another build needs more than {cap} grid cells (")].then(e.detail()).then(message![")"]));
                break;
            }
            other => other?,
        };
        match built {
            Attempt::Done(mut m) => {
                m.report.attempts = attempt + 1;
                m.report.retries = retries;
                return Ok(*m);
            }
            Attempt::Retry(r) => {
                retries.push(r.why);
                hints.extend(r.at);
            }
        }
    }
    Err(Error::failed(message![format!("octree hex meshing failed after {} builds: ", retries.len())].then(Message::join(retries, "; "))))
}

/// What the part's walls and narrow patches ask of the grid, measured before it is built.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Needs {
    /// The thinnest wall (`∞` when no two faces face each other).
    wall: f64,
    /// The narrowest patch (`∞` when no patch is bounded by curves).
    patch: f64,
    /// The largest cell edge that puts [`OctreeOptions::thickness_cells`] cells across every wall
    /// and two across every patch.
    edge: f64,
}

impl Needs {
    fn of(g: &Geometry, thickness: &[f64], options: &OctreeOptions) -> Needs {
        let wall = thickness.iter().copied().filter(|t| t.is_finite() && *t > 0.0).fold(f64::INFINITY, f64::min);
        let patch = g.patches.iter().map(|p| p.width).filter(|w| w.is_finite() && *w > 0.0).fold(f64::INFINITY, f64::min);
        let by_wall = if options.thickness_cells > 0.0 { wall / options.thickness_cells } else { f64::INFINITY };
        Needs { wall, patch, edge: by_wall.min(0.5 * patch) }
    }

    /// Why the cells must be as small as [`Needs::edge`].
    fn reason(&self, options: &OctreeOptions) -> Message {
        if options.thickness_cells > 0.0 && self.wall / options.thickness_cells <= 0.5 * self.patch {
            message![
                "the part has walls ",
                length(self.wall),
                format!(" thick, which {} cells across need cells of ", options.thickness_cells),
                length(self.edge)
            ]
        } else {
            message!["the part has a face ", length(self.patch), " wide, which two cells across need cells of ", length(self.edge)]
        }
    }
}

/// The levels below a root cell of edge `h` down to cells no larger than `edge` (`None` past
/// [`grid::MAX_LEVEL`]).
fn levels_below(h: f64, edge: f64) -> Option<u8> {
    (0..=grid::MAX_LEVEL).find(|&d| h / 3f64.powi(i32::from(d)) <= edge * (1.0 + 1e-9))
}

/// Choose the grid before building it. An octree's cell sizes are the root size divided by
/// powers of 3, so the cells across the thinnest wall come out up to 3 times smaller than the
/// wall needs, and the cells along every wall are counted in cubes of that. The root cells at
/// the requested size, with the walls refined as deep as they need, come first; when that does
/// not fit the levels or cells allowed, the root cells are shrunk (by less than 3 times) so that
/// one level fewer reaches exactly the size the walls need, which puts far fewer cells along
/// them. When neither fits, meshing fails here, before anything is built, naming the wall, the
/// levels and the cells each would take. Rebuilds after a failure refine within the levels the
/// plan reached, and at least one level below [`OctreeOptions::curvature_levels`].
fn plan<'a>(g: &'a Geometry<'a>, field: &'a dyn SizeField, options: &'a OctreeOptions, thickness: &'a [f64]) -> Result<Context<'a>> {
    let frame = choose_frame(g);
    let (lo, hi) = frame_box(g, &frame);
    let centre = from_frame(&frame, [0, 1, 2].map(|k| 0.5 * (lo[k] + hi[k])));
    // The root cells at the largest size; the depth reaches the smallest (size controls).
    let h = match field.range() {
        Some((_, coarsest)) if coarsest > 0.0 && coarsest.is_finite() => coarsest,
        _ => checked_size(field, centre)?,
    };
    if let Some((finest, _)) = field.range()
        && finest > 0.0
        && finest < h
    {
        match levels_below(h, finest) {
            Some(d) if d <= options.max_depth => {}
            d => {
                let n = d.map_or_else(|| format!("more than {}", grid::MAX_LEVEL), |l| l.to_string());
                return Err(Error::limit(message![
                    "a size of ",
                    length(finest),
                    " below the root size ",
                    length(h),
                    format!(" needs {n} levels (each a third); `max_depth` is {}: raise it, or use a larger local size", options.max_depth)
                ]));
            }
        }
    }
    let needs = Needs::of(g, thickness, options);
    let deep = levels_below(h, needs.edge);
    let mut anchors: Vec<(f64, Option<u8>)> = vec![(h, deep)];
    if let Some(d) = deep.filter(|&d| d >= 1) {
        anchors.push((needs.edge * 3f64.powi(i32::from(d) - 1), Some(d - 1)));
    }
    let mut best: Option<(usize, Context<'a>)> = None;
    let mut refused: Vec<Message> = Vec::new();
    for (root, levels) in anchors {
        let Some(levels) = levels.filter(|&l| l <= options.max_depth) else {
            let n = levels.map_or_else(|| format!("more than {}", grid::MAX_LEVEL), |l| l.to_string());
            refused.push(message![format!("{n} levels below a size of "), length(root)]);
            continue;
        };
        let ctx = Context::new(g, field, options, thickness, frame, root, (lo, hi), needs)?;
        match ctx.refine(options.max_elements)? {
            Refined::Within(tree) => {
                let mut closed = tree.clone();
                match closed.close(CLOSURE_ROUNDS, &|k| ctx.outside(k)) {
                    Ok(()) => {
                        let n = closed.leaf_count();
                        // Rebuilds refine within the levels the plan reached, and one level
                        // below the curvature levels.
                        let deepest = closed.level_counts().iter().rposition(|&c| c > 0).and_then(|l| u8::try_from(l).ok()).unwrap_or(0);
                        let depth = deepest.max(options.curvature_levels.saturating_add(1)).min(options.max_depth);
                        let mut tree = tree;
                        tree.set_depth(depth)?;
                        closed.set_depth(depth)?;
                        *ctx.refined_cells.borrow_mut() = tree.level_counts();
                        *ctx.refined.borrow_mut() = Some(tree);
                        *ctx.closed.borrow_mut() = Some(closed);
                        best = Some((n, ctx));
                        break;
                    }
                    Err(e) if e.is_limit() => {
                        refused.push(message![
                            format!("{levels} levels below a size of "),
                            length(root),
                            format!(", {} grid cells refined and more than {} after balancing and closure", tree.leaf_count(), options.max_elements)
                        ]);
                    }
                    Err(e) => return Err(e),
                }
            }
            Refined::Over { level, at_least } => {
                refused.push(message![
                    format!("{levels} levels below a size of "),
                    length(root),
                    format!(", more than {at_least} grid cells by level {level}")
                ]);
            }
        }
    }
    match best {
        Some((_, ctx)) => Ok(ctx),
        None if needs.edge.is_finite() => Err(Error::limit(needs.reason(options).then(message![": "]).then(Message::join(refused, "; or ")).then(
            message![format!(
                "; the limits are {} levels (`max_depth`) and {} cells (`max_elements`); use a smaller size, more levels or cells, or a structured (swept) hex mesh",
                options.max_depth, options.max_elements
            )],
        ))),
        None => Err(Error::limit(message!["the grid at a size of ", length(h), " needs "].then(Message::join(refused, "; or ")).then(message![
            format!("; the limits are {} levels (`max_depth`) and {} cells (`max_elements`)", options.max_depth, options.max_elements)
        ]))),
    }
}

/// The part's box in the frame's coordinates.
fn frame_box(g: &Geometry, frame: &[[f64; 3]; 3]) -> (P, P) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for &p in g.surface.points() {
        let u = to_frame(frame, p);
        for k in 0..3 {
            lo[k] = lo[k].min(u[k]);
            hi[k] = hi[k].max(u[k]);
        }
    }
    (lo, hi)
}

/// A grid refined by the criteria, or how far the refinement got before the cell limit.
enum Refined {
    Within(Tree),
    Over { level: u8, at_least: usize },
}

fn check_options(o: &OctreeOptions) -> Result<()> {
    let ok = o.max_depth <= 8
        && o.curvature_levels <= o.max_depth
        && o.max_elements >= 1
        && o.curvature_angle.is_finite()
        && o.curvature_angle > 0.0
        && o.thickness_cells.is_finite()
        && o.thickness_cells >= 0.0
        && o.core_offset.is_finite()
        && (0.0..1.0).contains(&o.core_offset)
        && (1..=16).contains(&o.attempts)
        && o.sharp_angle.is_finite()
        && o.sharp_angle > 0.0
        && o.blade_angle.is_finite()
        && o.blade_angle > o.sharp_angle
        && o.blade_angle <= 180.0
        && o.max_sweeps <= 100_000;
    if ok { Ok(()) } else { Err(Error::Invalid(format!("octree options out of range: {o:?}"))) }
}

fn from_frame(f: &[[f64; 3]; 3], u: P) -> P {
    [0, 1, 2].map(|k| u[0] * f[0][k] + u[1] * f[1][k] + u[2] * f[2][k])
}

fn to_frame(f: &[[f64; 3]; 3], p: P) -> P {
    [dot(p, f[0]), dot(p, f[1]), dot(p, f[2])]
}

/// The grid axes (rows): the coordinate axes, or the normals of the two largest orthogonal planar
/// patches when more boundary area is aligned with them.
fn choose_frame(g: &Geometry) -> [[f64; 3]; 3] {
    let tol = libm::sin(1f64.to_radians());
    let normal = |t: u32| g.normals.get(t as usize).copied().unwrap_or([0.0; 3]);
    let mut planar: Vec<(f64, usize, P)> = Vec::new();
    for (i, p) in g.patches.iter().enumerate() {
        let Some(&t0) = p.triangles.first() else { continue };
        let n0 = normal(t0);
        if p.triangles.iter().all(|&t| norm(cross(normal(t), n0)) < tol) {
            planar.push((p.area, i, n0));
        }
    }
    planar.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let Some(&(_, _, a0)) = planar.first() else { return identity };
    let Some(&(_, _, a1)) = planar.iter().find(|x| dot(x.2, a0).abs() < tol) else { return identity };
    let a1 = {
        let v = sub(a1, a0.map(|x| x * dot(a1, a0)));
        let l = norm(v);
        v.map(|x| x / l)
    };
    let candidate = [a0, a1, cross(a0, a1)];
    let score = |f: &[[f64; 3]; 3]| -> f64 {
        (0..g.surface.triangles().len() as u32)
            .filter(|&t| f.iter().any(|a| norm(cross(normal(t), *a)) < tol))
            .map(|t| {
                let [a, b, c] = g.triangle(t);
                0.5 * norm(cross(sub(b, a), sub(c, a)))
            })
            .sum()
    };
    if score(&candidate) > score(&identity) * (1.0 + 1e-9) { candidate } else { identity }
}

/// The distance through the part from each boundary triangle: a ray from its centre along the
/// inward normal to the far side of the wall, a triangle whose normal opens from the opposite of
/// this one by at most [`WALL_OPENING`] (`∞` when the ray first leaves through another face). A
/// ray that starts near a convex edge leaves through the face beyond the edge after a distance
/// that only measures how near the edge it started; the faces of a wall face each other.
fn triangle_thickness(g: &Geometry) -> Vec<f64> {
    let far = 2.0 * norm(sub(g.bounds.hi, g.bounds.lo)) + 1.0;
    let opposed = -libm::cos(WALL_OPENING.to_radians());
    (0..g.surface.triangles().len() as u32)
        .map(|t| {
            let n = g.normals.get(t as usize).copied().unwrap_or([0.0; 3]);
            match g.ray(g.triangle_centre(t), n.map(|x| -x), far, t) {
                Some((hit, s)) if g.normals.get(hit as usize).is_some_and(|&m| dot(m, n) <= opposed) => s,
                _ => f64::INFINITY,
            }
        })
        .collect()
}

fn centroid(points: &[P], hx: &[u32; 8]) -> P {
    [0, 1, 2].map(|k| hx.iter().map(|&v| points.get(v as usize).map_or(0.0, |p| p[k])).sum::<f64>() / 8.0)
}

fn edge_mean(points: &[P], hx: &[u32; 8]) -> f64 {
    let x = hx.map(|v| points.get(v as usize).copied().unwrap_or([0.0; 3]));
    let edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)];
    let m = edges.iter().map(|&(a, b)| norm(sub(x[a], x[b]))).sum::<f64>() / 12.0;
    if m > 0.0 { m } else { 1.0 }
}

/// The part being meshed and the state carried from one build to the next.
struct Context<'a> {
    g: &'a Geometry<'a>,
    field: &'a dyn SizeField,
    options: &'a OctreeOptions,
    frame: [[f64; 3]; 3],
    origin: P,
    h: f64,
    roots: [i64; 3],
    thickness: &'a [f64],
    needs: Needs,
    outside_memo: RefCell<HashMap<grid::Key, bool>>,
    /// The grid refined by the criteria (before hints and closure), and closed, kept from the
    /// plan: the first build takes the closed grid, rebuilds refine the other near their hints.
    refined: RefCell<Option<Tree>>,
    closed: RefCell<Option<Tree>>,
    /// Grid leaves of the last build.
    last_leaves: Cell<usize>,
    /// Leaves per level of the last build after the refinement criteria.
    refined_cells: RefCell<Vec<usize>>,
}

/// The mesh while its boundary layers are built: points, elements, their target sizes and how
/// each point may move.
struct Work {
    points: Vec<P>,
    hexes: Vec<[u32; 8]>,
    size: Vec<f64>,
    role: Vec<Constraint>,
    core_count: usize,
}

impl<'a> Context<'a> {
    /// Root cells of edge `h` over the part's box `(lo, hi)` (in the frame's coordinates).
    #[allow(clippy::too_many_arguments)]
    fn new(
        g: &'a Geometry<'a>,
        field: &'a dyn SizeField,
        options: &'a OctreeOptions,
        thickness: &'a [f64],
        frame: [[f64; 3]; 3],
        h: f64,
        (lo, hi): (P, P),
        needs: Needs,
    ) -> Result<Context<'a>> {
        if !(h.is_finite() && h > 0.0) {
            return Err(Error::invalid(message!["a root cell size of ", length(h)]));
        }
        let centre_u = [0, 1, 2].map(|k| 0.5 * (lo[k] + hi[k]));
        let mut roots = [0i64; 3];
        for k in 0..3 {
            let n = ((hi[k] - lo[k]) / h + 1.0).ceil();
            if !(n.is_finite() && n >= 1.0 && n <= options.max_elements as f64) {
                return Err(Error::limit(message!["a size of ", length(h), format!(" needs {n} root cells along an axis")]));
            }
            roots[k] = n as i64;
        }
        let origin = [0, 1, 2].map(|k| centre_u[k] - 0.5 * roots[k] as f64 * h);
        Ok(Context {
            g,
            field,
            options,
            frame,
            origin,
            h,
            roots,
            thickness,
            needs,
            outside_memo: RefCell::new(HashMap::new()),
            refined: RefCell::new(None),
            closed: RefCell::new(None),
            last_leaves: Cell::new(0),
            refined_cells: RefCell::new(Vec::new()),
        })
    }
}

impl Context<'_> {
    fn edge(&self, level: u8) -> f64 {
        self.h / 3f64.powi(i32::from(level))
    }

    /// The world box (inflated by `pad` cell edges) and the centre of a cell.
    fn cell_box(&self, (level, i): grid::Key, pad: f64) -> (Aabb, P) {
        let e = self.edge(level);
        let mut b = Aabb::EMPTY;
        for c in 0..8 {
            b.add(from_frame(&self.frame, [0, 1, 2].map(|k| self.origin[k] + (i[k] as f64 + ((c >> k) & 1) as f64) * e)));
        }
        let centre = from_frame(&self.frame, [0, 1, 2].map(|k| self.origin[k] + (i[k] as f64 + 0.5) * e));
        let p = pad * e;
        (Aabb { lo: b.lo.map(|x| x - p), hi: b.hi.map(|x| x + p) }, centre)
    }

    fn needs_split(&self, key: grid::Key) -> Result<bool> {
        let e = self.edge(key.0);
        let (bx, centre) = self.cell_box(key, 0.05);
        if e > checked_size(self.field, centre)? * (1.0 + 1e-9) {
            return Ok(true);
        }
        let reach = 0.5 * 3f64.sqrt() * e * 1.05;
        let cos = libm::cos(self.options.curvature_angle.to_radians());
        let mut first: BTreeMap<u32, P> = BTreeMap::new();
        for t in self.g.overlapping(&bx) {
            let [a, b, c] = self.g.triangle(t);
            if geom::dist2(centre, geom::closest_on_triangle(centre, a, b, c)) > reach * reach {
                continue;
            }
            let n = self.g.normals.get(t as usize).copied().unwrap_or([0.0; 3]);
            let face = self.g.surface.tri_face().get(t as usize).copied().unwrap_or(0);
            let f = *first.entry(face).or_insert(n);
            if key.0 < self.options.curvature_levels && dot(f, n) < cos {
                return Ok(true);
            }
            // (With a relative margin, so cells planned at exactly the size a wall needs stay.)
            if self.thickness.get(t as usize).is_some_and(|&th| th * (1.0 + 1e-9) < self.options.thickness_cells * e) {
                return Ok(true);
            }
            // Narrow patches (a hole's wall, a thin rim) need cells across them.
            if self.g.tri_patch.get(t as usize).and_then(|&p| self.g.patches.get(p as usize)).is_some_and(|pt| pt.width * (1.0 + 1e-9) < 2.0 * e) {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// Whether a cell lies wholly outside the part (its centre outside and farther from the
    /// boundary than its half diagonal): no core element can come from it, so its transition
    /// pattern does not matter and it is not meshed. Memoised.
    fn outside(&self, key: grid::Key) -> bool {
        if let Some(&v) = self.outside_memo.borrow().get(&key) {
            return v;
        }
        let (_, centre) = self.cell_box(key, 0.0);
        let half = 0.5 * 3f64.sqrt() * self.edge(key.0);
        let far = self.g.nearest(centre, |_| true).is_some_and(|x| x.2 > half * half * (1.0 + 1e-9));
        let v = far && self.g.inside(centre) == Some(false);
        self.outside_memo.borrow_mut().insert(key, v);
        v
    }

    /// The level-`level` cells within one cell of a hint point.
    fn hint_cells(&self, hints: &[P], level: u8) -> HashSet<grid::Key> {
        let e = self.edge(level);
        let mut out = HashSet::new();
        for &p in hints {
            let u = to_frame(&self.frame, p);
            let i = [0, 1, 2].map(|k| ((u[k] - self.origin[k]) / e).floor());
            if !i.iter().all(|x| x.is_finite() && x.abs() < 1e15) {
                continue;
            }
            for d in 0..27i64 {
                out.insert((level, [i[0] as i64 + d % 3 - 1, i[1] as i64 + (d / 3) % 3 - 1, i[2] as i64 + d / 9 - 1]));
            }
        }
        out
    }

    /// The grid refined level by level where [`Context::needs_split`] asks, within `max_leaves`
    /// cells; past them, how far it got and a lower bound on the cells it needs (the cells so far
    /// and the splits the level still asks for).
    fn refine(&self, max_leaves: usize) -> Result<Refined> {
        let roots = self.roots.iter().try_fold(1usize, |a, &n| usize::try_from(n).ok().and_then(|n| a.checked_mul(n))).unwrap_or(usize::MAX);
        let mut tree = match Tree::new(self.roots, self.options.max_depth, max_leaves) {
            Err(ref e) if e.kind() == Some(crate::Kind::Limit) => return Ok(Refined::Over { level: 0, at_least: roots }),
            other => other?,
        };
        for level in 0..self.options.max_depth {
            fieldcraft_model::progress::check()?;
            let keys: Vec<grid::Key> = tree.sorted_leaves().into_iter().filter(|k| k.0 == level).collect();
            for (i, &key) in keys.iter().enumerate() {
                if !self.needs_split(key)? {
                    continue;
                }
                match tree.split(key) {
                    Ok(()) => {}
                    Err(e) if e.kind() == Some(crate::Kind::Limit) && tree.leaf_count() + 26 > max_leaves => {
                        let mut more = 1usize;
                        for &k in keys.get(i + 1..).unwrap_or(&[]) {
                            if self.needs_split(k)? {
                                more += 1;
                            }
                        }
                        return Ok(Refined::Over { level: level + 1, at_least: tree.leaf_count().saturating_add(26usize.saturating_mul(more)) });
                    }
                    Err(e) => return Err(e),
                }
            }
        }
        Ok(Refined::Within(tree))
    }

    /// The refined, balanced and closed grid; leaves near `hints` (where an earlier build
    /// failed) are split once more.
    fn tree(&self, hints: &[P], max_leaves: usize) -> Result<Tree> {
        if hints.is_empty()
            && let Some(t) = self.closed.borrow_mut().take()
            && t.leaf_count() <= max_leaves
        {
            return Ok(t);
        }
        let kept = self.refined.borrow().clone();
        let mut tree = match kept {
            Some(t) => t,
            None => match self.refine(max_leaves)? {
                Refined::Within(t) => t,
                Refined::Over { level, at_least } => {
                    return Err(Error::Limit(format!("the grid needs more than {at_least} cells by level {level}; the limit is {max_leaves}")));
                }
            },
        };
        tree.set_max_leaves(max_leaves)?;
        if !hints.is_empty() {
            let near: Vec<HashSet<grid::Key>> = (0..tree.depth).map(|l| self.hint_cells(hints, l)).collect();
            let split: Vec<grid::Key> =
                tree.sorted_leaves().into_iter().filter(|k| near.get(usize::from(k.0)).is_some_and(|s| s.contains(k))).collect();
            for key in split {
                tree.split(key)?;
            }
        }
        *self.refined_cells.borrow_mut() = tree.level_counts();
        tree.close(CLOSURE_ROUNDS, &|k| self.outside(k))?;
        Ok(tree)
    }

    fn build(&self, hints: &[P], max_leaves: usize) -> Result<Attempt> {
        let g = self.g;
        let tree = self.tree(hints, max_leaves)?;
        self.last_leaves.set(tree.leaf_count());
        let lattice = tree.mesh(&|k| self.outside(k))?;
        if lattice.hexes.len() > self.options.max_elements {
            return Err(Error::Limit(format!("the grid has {} elements; the limit is {}", lattice.hexes.len(), self.options.max_elements)));
        }
        // A lattice coordinate over the root span is an exact ratio, rounded once, so a point's
        // position does not depend on how deep the lattice is.
        let span = tree.span(0) as f64;
        let points: Vec<P> =
            lattice.points.iter().map(|q| from_frame(&self.frame, [0, 1, 2].map(|k| self.origin[k] + q[k] as f64 / span * self.h))).collect();

        // Core.
        let dist: Vec<f64> = points.iter().map(|&p| g.nearest(p, |_| true).map_or(0.0, |x| x.2.sqrt())).collect();
        let keep = self.core(&lattice, &points, &dist);
        let core: Vec<usize> = (0..lattice.hexes.len()).filter(|&i| keep[i]).collect();
        if core.is_empty() {
            let at = (0..g.surface.triangles().len() as u32).map(|t| g.triangle_centre(t)).collect();
            return Ok(Attempt::Retry(Retry { at, why: message!["no grid element lies inside the part"] }));
        }
        let mut w = Work {
            hexes: core.iter().map(|&i| lattice.hexes[i]).collect(),
            size: core.iter().map(|&i| self.edge(lattice.level[i])).collect(),
            role: Vec::new(),
            points,
            core_count: core.len(),
        };
        let qb = topo::boundary(&w.hexes);
        let q: Vec<[u32; 4]> = qb.iter().map(|x| x.1).collect();
        if let Some(r) = check_topology(g, &q) {
            return Ok(Attempt::Retry(r));
        }

        // Feature map, padding and pillows.
        let fmap = {
            let pts = &w.points;
            match map::map_features(g, &q, &|v| pts.get(v as usize).copied().unwrap_or([0.0; 3])) {
                Ok(m) => m,
                Err(f) => return Ok(Attempt::Retry(Retry { at: f.at, why: f.why })),
            }
        };
        self.layers(&mut w, &qb, &fmap, &dist)?;
        if w.hexes.len() > self.options.max_elements {
            return Err(Error::Limit(format!("the mesh needs {} elements; the limit is {}", w.hexes.len(), self.options.max_elements)));
        }

        // Fit, then improve and repair.
        let mut opt = Optimizer::new(g, std::mem::take(&mut w.points), &w.hexes, &w.size, w.role.clone());
        let mut ran = 0;
        let sweeps = (self.options.max_sweeps / FIT_STAGES.len()).max(1);
        for &weight in FIT_STAGES {
            opt.fit = weight;
            ran = opt.run(sweeps, 1e-4);
        }
        let (_, boundary_distance) = opt.snap(MIN_SNAP_QUALITY);
        let inverted: Vec<P> = w
            .hexes
            .iter()
            .filter(|hx| smooth::min_scaled_jacobian(&hx.map(|v| opt.points[v as usize])) <= 0.0)
            .map(|hx| centroid(&opt.points, hx))
            .collect();
        if !inverted.is_empty() {
            let why = format!("{} elements stay inverted after untangling", inverted.len());
            return Ok(Attempt::Retry(Retry { at: inverted, why: why.into() }));
        }
        opt.improve_worst(IMPROVE_BELOW, 4);
        let uncertified = opt.repair_certified(6);
        if !uncertified.is_empty() {
            let at: Vec<P> = uncertified.iter().filter_map(|&e| w.hexes.get(e as usize)).map(|hx| centroid(&opt.points, hx)).collect();
            let why = format!("{} elements fail the certified validity check", at.len());
            return Ok(Attempt::Retry(Retry { at, why: why.into() }));
        }
        w.points = opt.points;
        self.finish(&w, &tree, ran, boundary_distance)
    }

    /// The core: elements wholly inside the part, cleaned to a 2-manifold boundary whose
    /// vertices can all move outward.
    fn core(&self, lattice: &grid::LatticeMesh, points: &[P], dist: &[f64]) -> Vec<bool> {
        let inside: Vec<bool> = points.iter().map(|&p| self.g.inside(p) == Some(true)).collect();
        let mut keep: Vec<bool> = lattice
            .hexes
            .iter()
            .zip(&lattice.level)
            .map(|(hx, &l)| {
                let floor = self.options.core_offset * self.edge(l);
                hx.iter().all(|&v| inside.get(v as usize) == Some(&true) && dist.get(v as usize).is_some_and(|&d| d >= floor))
            })
            .collect();
        for _ in 0..CORE_ROUNDS {
            let mut drop = topo::nonmanifold(&lattice.hexes, &keep);
            if drop.is_empty() {
                drop = topo::hidden_vertices(&lattice.hexes, &keep, &|v| points.get(v as usize).copied().unwrap_or([0.0; 3]), 0.1);
            }
            if drop.is_empty() {
                break;
            }
            for d in drop {
                if let Some(k) = keep.get_mut(d as usize) {
                    *k = false;
                }
            }
        }
        keep
    }

    /// The boundary layers: padding on the core's boundary quads `qb`, whose outer copies take the
    /// map's roles and start halfway to the boundary along the core's normals, then a pillow
    /// under each patch, whose copies start a quarter of the way into their elements. Every
    /// element starts valid or nearly so; the fit untangles the rest.
    fn layers(&self, w: &mut Work, qb: &[(topo::FaceRef, [u32; 4])], fmap: &map::FeatureMap, dist: &[f64]) -> Result<()> {
        let q: Vec<[u32; 4]> = qb.iter().map(|x| x.1).collect();
        let core_count = w.core_count;
        let copy = topo::pad(&mut w.points, &mut w.hexes, &q)?;
        w.role = vec![Constraint::Free; w.points.len()];
        let mut normal: BTreeMap<u32, P> = BTreeMap::new();
        for f in &q {
            let p = f.map(|v| w.points[v as usize]);
            let n = cross(sub(p[2], p[0]), sub(p[3], p[1]));
            let l = norm(n);
            if l > 0.0 {
                for &v in f {
                    let e = normal.entry(v).or_insert([0.0; 3]);
                    for k in 0..3 {
                        e[k] += n[k] / l;
                    }
                }
            }
        }
        for (&v, &c) in &copy {
            if let Some(r) = fmap.role.get(&v) {
                w.role[c as usize] = *r;
            }
            let n = normal.get(&v).copied().unwrap_or([0.0; 3]);
            let l = norm(n);
            let d = 0.5 * dist.get(v as usize).copied().unwrap_or(0.0);
            if l > 0.0 {
                w.points[c as usize] = [0, 1, 2].map(|k| w.points[v as usize][k] + d * n[k] / l);
            }
        }
        for ((owner, _), _) in qb {
            w.size.push(0.5 * w.size.get(*owner as usize).copied().unwrap_or(self.h));
        }
        let patches: BTreeSet<u32> = fmap.quad_patch.iter().copied().collect();
        for p in patches {
            let inside: Vec<bool> =
                (0..w.hexes.len()).map(|i| i >= core_count && i < core_count + q.len() && fmap.quad_patch.get(i - core_count) == Some(&p)).collect();
            let before = w.hexes.len();
            let centroids: Vec<P> = w.hexes.iter().map(|hx| centroid(&w.points, hx)).collect();
            let copies = topo::pillow(&mut w.points, &mut w.hexes, &inside)?;
            w.role.resize(w.points.len(), Constraint::Free);
            for pc in &copies {
                let orig = w.points[pc.original as usize];
                let n = pc.elements.len().max(1) as f64;
                let mut c = [0.0; 3];
                for &e in &pc.elements {
                    let ce = centroids.get(e as usize).copied().unwrap_or(orig);
                    for k in 0..3 {
                        c[k] += ce[k] / n;
                    }
                }
                w.points[pc.copy as usize] = [0, 1, 2].map(|k| orig[k] + 0.25 * (c[k] - orig[k]));
                w.role[pc.copy as usize] = match w.role.get(pc.original as usize).copied().unwrap_or(Constraint::Free) {
                    Constraint::Free => Constraint::Free,
                    _ => Constraint::Patch(p),
                };
            }
            for hx in &w.hexes[before..] {
                w.size.push(edge_mean(&w.points, hx));
            }
        }
        Ok(())
    }

    /// The gate and the output: every element proven valid, nodes renumbered, boundary faces with
    /// their B-rep faces.
    fn finish(&self, w: &Work, tree: &Tree, sweeps: usize, boundary_distance: f64) -> Result<Attempt> {
        let g = self.g;
        // One certified pass: every element's quality (the same metric as mesh.quality and
        // mesh.hex.generate report) and its Jacobian check for the gate.
        let mut bad: Vec<P> = Vec::new();
        let mut checks = Vec::with_capacity(w.hexes.len());
        let (mut sj_min, mut sj_sum) = (f64::INFINITY, 0.0);
        for (e, hx) in w.hexes.iter().enumerate() {
            let q = element_quality(ElementKind::Hex8, &hx.map(|v| w.points[v as usize]))?;
            let check = q.jacobian.ok_or_else(|| Error::Failed("an element was not certified".into()))?;
            if check.validity != Validity::Valid {
                bad.push(centroid(&w.points, hx));
            }
            let sj = q.get(Metric::ScaledJacobian).unwrap_or(-1.0);
            sj_min = sj_min.min(sj);
            sj_sum += sj;
            checks.push((e as u32 + 1, check));
        }
        if !bad.is_empty() {
            let why = format!("{} elements fail the certified validity check", bad.len());
            return Ok(Attempt::Retry(Retry { at: bad, why: why.into() }));
        }
        // The volume the mesh's own boundary encloses equals the elements' total exactly when they
        // tile it without overlapping.
        let faces = topo::boundary(&w.hexes);
        let enclosed: f64 = faces.iter().map(|(_, quad)| face_flux(&quad.map(|v| w.points[v as usize]))).sum();
        let gate =
            validity_gate_checks(checks.iter().map(|(l, c)| (*l, c)), &GateOptions { expected_volume: Some(enclosed), volume_tolerance: 1e-9 })?
                .into_result()?;
        let volume = gate.volume;
        let mut new_id = vec![u32::MAX; w.points.len()];
        let mut nodes = Vec::new();
        let mut hexes = Vec::with_capacity(w.hexes.len());
        for hx in &w.hexes {
            let mut o = [0u32; 8];
            for (k, &v) in hx.iter().enumerate() {
                if new_id[v as usize] == u32::MAX {
                    new_id[v as usize] = u32::try_from(nodes.len()).map_err(|_| Error::Limit("too many nodes".into()))?;
                    nodes.push(w.points[v as usize]);
                }
                o[k] = new_id[v as usize];
            }
            hexes.push(o);
        }
        let mut boundary = Vec::new();
        for ((e, f), quad) in faces {
            let patch = quad.iter().find_map(|&v| match w.role.get(v as usize) {
                Some(Constraint::Patch(p)) => Some(*p),
                _ => None,
            });
            let c = [0, 1, 2].map(|k| quad.iter().map(|&v| w.points[v as usize][k]).sum::<f64>() / 4.0);
            let brep_face = match patch.and_then(|p| g.patches.get(p as usize)) {
                Some(pt) if pt.faces.len() == 1 => pt.faces[0],
                _ => {
                    let keep = |t: u32| patch.is_none_or(|p| g.tri_patch.get(t as usize) == Some(&p));
                    g.nearest(c, keep).and_then(|(t, _, _)| g.surface.tri_face().get(t as usize).copied()).unwrap_or(0)
                }
            };
            boundary.push(BoundaryFace { element: e, face: f + 1, brep_face });
        }
        let report = OctreeReport {
            attempts: 0,
            grid_cells: tree.leaf_count(),
            root_size: self.h,
            levels: tree.level_counts().iter().rposition(|&n| n > 0).and_then(|l| u8::try_from(l).ok()).unwrap_or(0),
            thinnest_wall: self.needs.wall,
            narrowest_patch: self.needs.patch,
            refined_cells: self.refined_cells.borrow().clone(),
            closed_cells: tree.level_counts(),
            core_elements: w.core_count,
            layer_elements: w.hexes.len() - w.core_count,
            frame: self.frame,
            corners: g.corners.len(),
            curves: g.curves.len(),
            patches: g.patches.len(),
            sweeps,
            min_scaled_jacobian: sj_min,
            mean_scaled_jacobian: sj_sum / w.hexes.len().max(1) as f64,
            volume,
            surface_volume: g.surface.volume(),
            boundary_distance,
            retries: Vec::new(),
        };
        Ok(Attempt::Done(Box::new(OctreeHexMesh { nodes, hexes, boundary, report })))
    }
}

/// The core boundary's components and genus against the part's shells; on a mismatch, a retry
/// that refines near the whole boundary.
fn check_topology(g: &Geometry, quads: &[[u32; 4]]) -> Option<Retry> {
    let edges = topo::quad_edges(quads);
    let mut comp: Vec<usize> = (0..quads.len()).collect();
    for qs in edges.values() {
        for w in qs.windows(2) {
            let (a, b) = (topo::find(&mut comp, w[0] as usize), topo::find(&mut comp, w[1] as usize));
            comp[a] = b;
        }
    }
    // Per component: vertices, edges, faces.
    let mut stats: BTreeMap<usize, (BTreeSet<u32>, i64, i64)> = BTreeMap::new();
    for (i, q) in quads.iter().enumerate() {
        let s = stats.entry(topo::find(&mut comp, i)).or_default();
        s.2 += 1;
        s.0.extend(q.iter().copied());
    }
    for qs in edges.values() {
        if let Some(&q0) = qs.first()
            && let Some(s) = stats.get_mut(&topo::find(&mut comp, q0 as usize))
        {
            s.1 += 1;
        }
    }
    let mut got: Vec<i64> = stats.values().map(|s| (2 - (s.0.len() as i64 - s.1 + s.2)) / 2).collect();
    let mut want: Vec<i64> = g.surface.shells().iter().map(|s| i64::from(s.genus)).collect();
    got.sort_unstable();
    want.sort_unstable();
    if got == want {
        return None;
    }
    let at = (0..g.surface.triangles().len() as u32).map(|t| g.triangle_centre(t)).collect();
    Some(Retry { at, why: message![format!("the core's boundary has components of genus {got:?}, the part {want:?}")] })
}
