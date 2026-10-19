//! Planar all-quadrilateral meshing: mapped (transfinite) meshes of four-sided regions and free
//! meshes of polygons with holes, on prescribed boundary nodes.
//!
//! These are the face meshers under the mapped and swept hexahedral meshers and shell meshing.
//! Both take the boundary nodes as given (shared edges must share nodes), put them first in the
//! output in input order, and return counterclockwise `Quad4` connectivity whose every element is
//! proven valid by the exact corner test ([`crate::quality::planar_quad_validity`]).
//!
//! # Mapped: [`mapped_quad_mesh`], [`transfinite_grid`]
//! Transfinite interpolation (Gordon and Hall 1973), see [`tfi`]. Opposite sides need equal
//! division counts; a region whose TFI folds (strongly concave sides) is an error, not a tangled
//! mesh.
//!
//! # Free: [`free_quad_mesh`]
//! **Method.** An indirect tri-to-quad method in three steps:
//! 1. Keep every other boundary point of each loop (the parity class with the sharper turns, so
//!    corners at even spacing survive) and build the constrained Delaunay triangulation of that
//!    half-resolution region, refined to twice the target size ([`cdt`]: Bowyer–Watson, Sloan
//!    1993 recovery, Ruppert 1995 refinement without segment splitting).
//! 2. Recombine adjacent triangles into strictly convex quadrilaterals, greedily, best first by
//!    the smallest corner scaled Jacobian (Lee and Lo, *A new scheme for the generation of a
//!    graded quadrilateral mesh*, Computers & Structures 52(5):847–857, 1994).
//! 3. Split every element once: a quadrilateral into four, a leftover triangle into three, through
//!    edge mid nodes and a centre node. On a boundary edge the mid node is the prescribed point the
//!    half-resolution loop skipped, so the final boundary nodes are exactly the input points. The
//!    split is the Catmull–Clark topology (Catmull and Clark, CAD 10(6):350–355, 1978) and appears
//!    as the all-quad completion of recombination in Lee and Lo 1994.
//!
//! Then constrained Laplacian smoothing ([`smooth`]) and the exact validity gate.
//!
//! **Why this method.** The alternatives considered, from their papers:
//! - *Q-Morph* (Owen, Staten, Canann, Saigal, *Q-Morph: an indirect approach to advancing front
//!   quad meshing*, IJNME 44(9):1317–1340, 1999) gives better-aligned rows along the boundary, but
//!   needs front-edge state classification, side-edge recovery by swaps, seam closing and special
//!   handling where fronts collide; the paper reports the cases that need it, and each is a place
//!   to fail on concave faces with many holes. An odd front remainder still ends in a triangle
//!   to be removed.
//! - *Paving* (Blacker and Stephenson, IJNME 32(4):811–847, 1991) has the same front-collision
//!   problem with no background triangulation to fall back on; it is the highest-risk choice.
//! - *Blossom-Quad* (Remacle, Lambrechts, Seny, Marchandise, Johnen, Geuzaine, *Blossom-Quad: a
//!   non-uniform quadrilateral mesh generator using a minimum-cost perfect-matching algorithm*,
//!   IJNME 89(9):1102–1119, 2012; paper only, Gmsh's code is GPL) gives the best matching, but a
//!   perfect matching of the triangle graph need not exist, so it still needs a fallback, and it
//!   needs Edmonds' weighted blossom algorithm (cubic in the worst case).
//!
//! The split method is **all-quad by construction**: any triangulation of the half-resolution
//! region becomes quadrilaterals only, whatever the matching found, so the only requirement is an
//! even number of boundary points on every loop. Robustness on concave regions with holes is
//! that of the constrained Delaunay triangulation, which exists for every valid input and is
//! built here with exact predicates. The implementation risk is low: one triangulator (which the
//! surface and tetrahedral work need anyway), a greedy matching and a local split. The cost is
//! quality: the split creates irregular nodes (valence 3 at triangle centres, 6 at some coarse
//! vertices) and triangle-born quads have a 120° corner; smoothing removes most of the
//! distortion. If that proves too weak, the greedy matching can be replaced by Blossom-Quad's
//! minimum-cost matching without touching the rest (fewer leftover triangles, same guarantee).
//!
//! **Requirements and failure modes.**
//! - Every loop has an even number of points, at least 6 when there is more than one loop (a
//!   single 4-point loop is meshed as one quadrilateral). An odd loop is an error, never repaired:
//!   the points are prescribed by neighbouring faces, so adding one would break conformity. Choose
//!   divisions with [`crate::size::loop_divisions`] or [`crate::size::even_divisions`].
//! - The chords joining every other boundary point must not cross (a boundary too coarse for a
//!   thin feature); both parity classes are tried before this is reported.
//! - Every output element is checked with exact predicates; a mesh that smoothing cannot make
//!   valid is an error.
//! - Hostile input (NaN, coordinates beyond the model limit, repeated points, zero-length edges,
//!   crossing or touching loops, holes outside the outer loop or inside another hole) is an
//!   [`Error::Invalid`]; size and work limits are [`Error::Limit`].

mod cdt;
mod domain;
mod free;
mod guard;
pub(crate) mod smooth;
pub mod tfi;

pub use domain::{MAX_BOUNDARY_POINTS, MAX_LOOPS};
pub use tfi::{MAX_GRID_NODES, lattice_index, transfinite_faces, transfinite_grid, transfinite_grid_blended, transfinite_volume};

use crate::predicates::{Point2, Sign, snap_to_domain};
use crate::quality::{planar_corner, quad_axes, turns_left};
use crate::size::SizeField;
use crate::{Error, Result};

/// Limits and smoothing of the quad meshers.
#[derive(Debug, Clone, PartialEq)]
pub struct QuadOptions {
    /// The most quadrilaterals of one mesh ([`Error::Limit`] beyond).
    pub max_elements: usize,
    /// Smoothing sweeps of the free mesher (0 turns smoothing off).
    pub smoothing_sweeps: usize,
}

impl Default for QuadOptions {
    fn default() -> Self {
        QuadOptions { max_elements: 4_000_000, smoothing_sweeps: 10 }
    }
}

/// What a quad mesher produced and how good it is.
#[derive(Debug, Clone, PartialEq)]
pub struct QuadStats {
    /// Smallest corner scaled Jacobian over the mesh (1 for squares).
    pub min_scaled_jacobian: f64,
    /// Largest skew (0 for rectangles).
    pub max_skew: f64,
    /// Largest aspect ratio (1 for squares).
    pub max_aspect_ratio: f64,
    /// Smoothing sweeps run, and whether smoothing converged before its sweep limit.
    pub smoothing_sweeps: usize,
    pub smoothing_converged: bool,
    /// Free mesher: triangle pairs recombined, and triangles split on their own.
    pub recombined: usize,
    pub leftover_triangles: usize,
}

impl Default for QuadStats {
    fn default() -> Self {
        QuadStats {
            min_scaled_jacobian: 1.0,
            max_skew: 0.0,
            max_aspect_ratio: 1.0,
            smoothing_sweeps: 0,
            smoothing_converged: true,
            recombined: 0,
            leftover_triangles: 0,
        }
    }
}

/// A planar all-quadrilateral mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct QuadMesh {
    /// Node coordinates. `nodes[..boundary_nodes]` are the input boundary points, unchanged and
    /// in input order (loop by loop), so node `k` of the mesh is boundary point `k`.
    pub nodes: Vec<[f64; 2]>,
    /// Counterclockwise corner node ids (the `Quad4` node order).
    pub quads: Vec<[u32; 4]>,
    pub boundary_nodes: usize,
    pub stats: QuadStats,
}

/// A mapped mesh with its structure: node `(i, j)`, `0 ≤ i ≤ nu`, `0 ≤ j ≤ nv`, is
/// `mesh.nodes[grid[j·(nu + 1) + i]]`, with `i` running from corner 0 to corner 1 and `j` from
/// corner 1 to corner 2 (for a sweep: the source face's structured layout).
#[derive(Debug, Clone, PartialEq)]
pub struct MappedMesh {
    pub mesh: QuadMesh,
    pub nu: usize,
    pub nv: usize,
    pub grid: Vec<u32>,
}

/// Free all-quadrilateral mesh of a planar region (see the module docs for the method).
///
/// `loops[0]` is the outer boundary, the others are holes; each loop lists its boundary nodes
/// once (not closed), in either orientation, with an even count. `size` sets the interior
/// density; the boundary density is the given nodes.
pub fn free_quad_mesh(loops: &[Vec<[f64; 2]>], size: &dyn SizeField, opts: &QuadOptions) -> Result<QuadMesh> {
    free::mesh(loops, size, opts)
}

/// A planar triangle mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct TriMesh {
    /// Node coordinates; `nodes[..boundary_nodes]` are the input boundary points in input order.
    pub nodes: Vec<[f64; 2]>,
    /// Counterclockwise corner node ids (the `Tri3` node order).
    pub triangles: Vec<[u32; 3]>,
    pub boundary_nodes: usize,
}

/// Free triangle mesh of a planar region on prescribed boundary nodes: the constrained Delaunay
/// triangulation of the loops ([`cdt`]: Bowyer–Watson, Sloan 1993 recovery), refined to the size
/// field by Ruppert-style circumcentre insertion that never splits a boundary edge (Steiner points
/// keep a clearance of half the local boundary spacing from the boundary), so the boundary nodes
/// are exactly the input points. Every triangle is checked counterclockwise with exact
/// predicates. Loops as for [`free_quad_mesh`], with any number of points (at least 3 each).
pub fn free_tri_mesh(loops: &[Vec<[f64; 2]>], size: &dyn SizeField, opts: &QuadOptions) -> Result<TriMesh> {
    let region = domain::region(loops)?;
    let mut chords = Vec::new();
    for l in &region.loops {
        for (i, &a) in l.iter().enumerate() {
            let b = l.get((i + 1) % l.len()).copied().ok_or_else(|| Error::Failed("loop index".into()))?;
            let (pa, pb) = (*region.point(a)?, *region.point(b)?);
            let (ca, cb) = (pa.coords(), pb.coords());
            let m = Point2::new(snap_to_domain(0.5 * (ca[0] + cb[0])), snap_to_domain(0.5 * (ca[1] + cb[1])))?;
            chords.push((pa, pb, m));
        }
    }
    let guard = guard::Guard::new(&chords)?;
    let mut cdt = cdt::Cdt::build(&region, opts.max_elements.saturating_add(4 * region.pts.len() + 16))?;
    cdt.refine(size, 1.0, &|p, h| guard.allows(p, cdt::CLEARANCE * h))?;
    let tri = cdt.into_triangulation();
    if tri.tris.len() > opts.max_elements {
        return Err(Error::Limit(format!("the mesh would have {} triangles; the limit is {}", tri.tris.len(), opts.max_elements)));
    }
    for t in &tri.tris {
        let [a, b, c] = t.map(|v| tri.p2.get(v as usize));
        match (a, b, c) {
            (Some(a), Some(b), Some(c)) if crate::predicates::orient2d(a, b, c) == Sign::Positive => {}
            _ => return Err(Error::Failed("a triangle is inverted or degenerate".into())),
        }
    }
    Ok(TriMesh { nodes: tri.pts, triangles: tri.tris, boundary_nodes: region.pts.len() })
}

/// Mapped mesh of a four-sided region: `boundary` lists the region's boundary nodes once (either
/// orientation), and `corners` are the positions in it of the four corners in boundary order.
/// Opposite sides need the same number of divisions. The grid is transfinite interpolation
/// ([`transfinite_grid`]); a grid that folds is an error.
pub fn mapped_quad_mesh(boundary: &[[f64; 2]], corners: [usize; 4], opts: &QuadOptions) -> Result<MappedMesh> {
    let region = domain::region(&[boundary.to_vec()])?;
    let n = boundary.len();
    if let Some(c) = corners.iter().find(|&&c| c >= n) {
        return Err(Error::Invalid(format!("corner position {c} is outside the boundary of {n} points")));
    }
    let c0 = corners[0];
    let d = corners.map(|c| (c + n - c0) % n);
    if !(d[0] == 0 && d[0] < d[1] && d[1] < d[2] && d[2] < d[3]) {
        return Err(Error::Invalid("the corners must be four distinct positions in boundary order".into()));
    }
    let (nu, nv) = (d[1], d[2] - d[1]);
    if d[3] - d[2] != nu || n - d[3] != nv {
        return Err(Error::Invalid(format!("opposite sides need the same number of divisions: {} and {}, {} and {}", nu, d[3] - d[2], nv, n - d[3])));
    }
    let at = |k: usize| region.pts.get((c0 + k) % n).copied().ok_or_else(|| Error::Failed("boundary index".into()));
    let side = |from: usize, len: usize, forward: bool| -> Result<Vec<[f64; 2]>> {
        (0..=len).map(|k| if forward { at(from + k) } else { at(from + n - k) }).collect()
    };
    // Boundary positions relative to corner 0: bottom 0..nu, right nu..nu+nv, top runs backwards
    // from d[3], left backwards from 0 (= n).
    let bottom = side(0, nu, true)?;
    let right = side(d[1], nv, true)?;
    let top = side(d[3], nu, false)?;
    let left = side(n, nv, false)?;
    let points = transfinite_grid(&bottom, &right, &top, &left)?;
    // Node ids: boundary nodes keep their input positions; interior nodes follow, row by row.
    let mut nodes = region.pts.clone();
    let mut grid = Vec::with_capacity(points.len());
    for j in 0..=nv {
        for i in 0..=nu {
            let pos = if j == 0 {
                Some(i)
            } else if i == nu {
                Some(nu + j)
            } else if j == nv {
                Some(d[3] - i)
            } else if i == 0 {
                Some((n - j) % n)
            } else {
                None
            };
            let id = match pos {
                Some(p) => ((c0 + p) % n) as u32,
                None => {
                    let x = points.get(j * (nu + 1) + i).copied().ok_or_else(|| Error::Failed("grid index".into()))?;
                    nodes.push(domain::admissible(x)?);
                    nodes.len() as u32 - 1
                }
            };
            grid.push(id);
        }
    }
    let ccw = domain::loop_orientation(&region, &(0..n as u32).collect::<Vec<_>>())? == Sign::Positive;
    let g = |i: usize, j: usize| grid.get(j * (nu + 1) + i).copied().ok_or_else(|| Error::Failed("grid index".into()));
    let mut quads = Vec::with_capacity(nu * nv);
    for j in 0..nv {
        for i in 0..nu {
            let (a, b, c, e) = (g(i, j)?, g(i + 1, j)?, g(i + 1, j + 1)?, g(i, j + 1)?);
            quads.push(if ccw { [a, b, c, e] } else { [a, e, c, b] });
        }
    }
    if quads.len() > opts.max_elements {
        return Err(Error::Limit(format!("the mesh would have {} quadrilaterals; the limit is {}", quads.len(), opts.max_elements)));
    }
    let mesh = finish(nodes, quads, n, &QuadOptions { smoothing_sweeps: 0, ..opts.clone() }, QuadStats::default())?;
    Ok(MappedMesh { mesh, nu, nv, grid })
}

/// Smooth (free nodes only), gate every element with the exact validity test and measure.
pub(crate) fn finish(
    mut nodes: Vec<[f64; 2]>,
    quads: Vec<[u32; 4]>,
    boundary_nodes: usize,
    opts: &QuadOptions,
    mut stats: QuadStats,
) -> Result<QuadMesh> {
    // Constructed nodes (mid nodes, centres) stay in the predicates' domain.
    for n in nodes.iter_mut() {
        *n = [snap_to_domain(n[0]), snap_to_domain(n[1])];
    }
    if opts.smoothing_sweeps > 0 {
        let (sweeps, converged) = smooth::smooth(&mut nodes, &quads, boundary_nodes, opts.smoothing_sweeps)?;
        stats.smoothing_sweeps = sweeps;
        stats.smoothing_converged = converged;
    }
    let mut invalid = 0usize;
    let mut first_invalid = None;
    let (mut sj, mut skew, mut aspect) = (f64::INFINITY, 0.0_f64, 1.0_f64);
    for (e, q) in quads.iter().enumerate() {
        let mut c = [[0.0; 2]; 4];
        for (k, &v) in q.iter().enumerate() {
            c[k] = nodes.get(v as usize).copied().ok_or_else(|| Error::Failed("node index".into()))?;
        }
        let [a, b, d, f] = c.map(|x| Point2::new(x[0], x[1]));
        if !turns_left([&a?, &b?, &d?, &f?]) {
            invalid += 1;
            first_invalid.get_or_insert(e);
        }
        sj = sj.min((0..4).fold(f64::INFINITY, |m, i| m.min(planar_corner(&c, i).1.unwrap_or(0.0))));
        let (r, k) = quad_axes(c.map(|x| [x[0], x[1], 0.0]));
        aspect = aspect.max(r);
        skew = skew.max(k);
    }
    if invalid > 0 {
        let at = first_invalid.and_then(|e| quads.get(e)).and_then(|q| nodes.get(q[0] as usize)).copied().unwrap_or_default();
        return Err(Error::failed(crate::message![
            format!("{invalid} of {} quadrilaterals are inverted or degenerate (first near ", quads.len()),
            crate::message::point(at),
            "); refine the boundary there or use a finer size"
        ]));
    }
    if !quads.is_empty() {
        stats.min_scaled_jacobian = sj;
        stats.max_skew = skew;
        stats.max_aspect_ratio = aspect;
    }
    Ok(QuadMesh { nodes, quads, boundary_nodes, stats })
}

#[cfg(test)]
mod tests;
