//! Mapped blocks: a structured mesh of one six-sided region given by its eight corners and,
//! optionally, curved edges (circular and elliptical arcs).
//!
//! The lattice follows the transfinite interpolation shared with the structured hex mesher
//! ([`crate::quad::transfinite_faces`], Gordon and Hall 1973, for the faces;
//! [`crate::quad::transfinite_volume`], Gordon and Thiel 1982, for the interior): the edge curves
//! at their (graded) parameters, each face the Coons patch of its four edges and the volume
//! blended from the faces, both by the same lattice parameters. It reproduces every edge exactly
//! and gives the trilinear map for straight edges. A 2D patch is the Coons patch of its four
//! edges.
//!
//! Elements: hex8 and hex20 per cell; tet4 and tet10 by Kuhn's split of each cell into six
//! tetrahedra along its main diagonal (conforming between the cells of a block); wedge6 by
//! splitting each cell's first-second-direction face along its diagonal and extruding along the
//! third. Quadratic elements take their mid-side nodes from the same map, so they lie on the
//! curved edges.
//!
//! Every element is checked with the certified Jacobian test ([`crate::quality`]); a block
//! whose corners are not in hex8 order or whose curves fold it is an error.

use fieldcraft_model::{ElementKind, MAX_COORD};

use crate::quad::{lattice_index, transfinite_faces, transfinite_grid_blended, transfinite_volume};
use crate::quality::{Validity, element_validity, planar_quad_validity};
use crate::vec3::{P, cross, dot, norm, sub};
use crate::{Error, Result};

/// Points on an arc within this, relative to its radius.
const ON_ARC: f64 = 1e-9; // unit-factors: relative to the radius
/// An arc's plane needs this much sine between its sides, relative to the radius squared.
const MIN_SINE: f64 = 1e-12; // unit-factors: relative to the radius squared
/// The least sweep of an arc.
const MIN_SWEEP: f64 = 1e-12; // unit-factors: an angle in radians
/// An arc short of a half turn by less than this is ambiguous.
const HALF_TURN_SLACK: f64 = 1e-9; // unit-factors: an angle in radians
/// The smallest grading ratio a block takes.
const MIN_BIAS: f64 = 1e-6; // unit-factors: a ratio of element sizes
/// The largest grading ratio a block takes.
const MAX_BIAS: f64 = 1e6; // unit-factors: a ratio of element sizes

/// Most elements one block makes.
pub const MAX_BLOCK_ELEMENTS: usize = 4_000_000;

/// Most divisions along one direction.
pub const MAX_DIVISIONS: usize = 100_000;

/// The corners of the block's twelve edges, grouped by direction: four along 1→2, four along
/// 1→4, four along 1→5, each from the corner where its direction's parameter is 0 to the one
/// where it is 1. Block edge numbers in [`BlockMesh::edges`] index this table.
pub const EDGES: [[usize; 2]; 12] = [[0, 1], [3, 2], [4, 5], [7, 6], [0, 3], [1, 2], [4, 7], [5, 6], [0, 4], [1, 5], [2, 6], [3, 7]];

/// Corner `c` of the block in parameter space (hex8 order).
const CORNERS: [[usize; 3]; 8] = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];

/// The shape of one block edge between its two corners.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Curve {
    Line,
    /// The shorter circular arc about `center` (both corners at the same distance from it).
    Arc {
        center: P,
    },
    /// The shorter arc of the ellipse `center + a cos t + b sin t` (`a` and `b` perpendicular)
    /// between the two corners, which must lie on it.
    Ellipse {
        center: P,
        axes: [P; 2],
    },
}

/// One block to mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct Block {
    pub kind: ElementKind,
    /// The corners in hex8 order: 1-2-3-4 counter-clockwise seen from 5-6-7-8's side.
    pub corners: [P; 8],
    /// Elements along the directions 1→2, 1→4 and 1→5.
    pub divisions: [usize; 3],
    /// Last element size over first along each direction (geometric grading); 1 is uniform.
    pub bias: [f64; 3],
    /// The curve of each edge, by corner pair (0-based corners).
    pub curves: Vec<([usize; 2], Curve)>,
}

/// The mesh of a block: points, elements as point indices in the model's node order, and where
/// element faces and edges lie on the block's faces and edges.
#[derive(Debug, Clone, PartialEq)]
pub struct BlockMesh {
    pub points: Vec<P>,
    pub elements: Vec<Vec<u32>>,
    /// `(element, face 1-based, block face 0..6)`: block faces in hex8 face order.
    pub faces: Vec<(u32, u8, u8)>,
    /// `(element, edge 1-based, block edge 0..12)`: block edges as in [`EDGES`].
    pub edges: Vec<(u32, u8, u8)>,
    /// The point at each block corner (a patch fills the first four).
    pub corners: [u32; 8],
}

impl Block {
    /// The block edge ([`EDGES`] index) joining two 0-based corners, either way round.
    pub fn edge_index(corners: [usize; 2]) -> Option<usize> {
        EDGES.iter().position(|e| *e == corners || *e == [corners[1], corners[0]])
    }
}

fn add_scaled(a: P, s: f64, b: P) -> P {
    [a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2]]
}

/// An edge curve resolved against its end points: `point(s)` for `s` in `[0, 1]`.
#[derive(Debug, Clone, Copy)]
enum Path {
    Line {
        from: P,
        to: P,
    },
    /// The ends are kept so the curve passes through the corners exactly (not to rounding).
    Ellipse {
        center: P,
        a: P,
        b: P,
        start: f64,
        sweep: f64,
        ends: [P; 2],
    },
}

impl Path {
    fn new(curve: Curve, from: P, to: P) -> Result<Path> {
        let (center, a, b) = match curve {
            Curve::Line => return Ok(Path::Line { from, to }),
            Curve::Arc { center } => {
                let (p, q) = (sub(from, center), sub(to, center));
                let r = norm(p);
                if !(r > 0.0) || (norm(q) - r).abs() > ON_ARC * r.max(norm(q)) {
                    return Err(Error::Invalid(format!(
                        "an arc's corners {from:?} and {to:?} must be at the same distance from its centre {center:?}"
                    )));
                }
                let n = cross(p, q);
                if !(norm(n) > MIN_SINE * r * r) {
                    return Err(Error::Invalid(format!(
                        "the arc from {from:?} to {to:?} about {center:?} has no plane (it spans 0 or 180 degrees); split it"
                    )));
                }
                let b = cross(n, p);
                (center, p, [b[0] * r / norm(b), b[1] * r / norm(b), b[2] * r / norm(b)])
            }
            Curve::Ellipse { center, axes: [a, b] } => (center, a, b),
        };
        let (la, lb) = (norm(a), norm(b));
        if !(la > 0.0 && lb > 0.0) || dot(a, b).abs() > ON_ARC * la * lb {
            return Err(Error::Invalid("an ellipse's two axes must be non-zero and perpendicular".into()));
        }
        let angle = |p: P| -> Result<f64> {
            let d = sub(p, center);
            let (c, s) = (dot(d, a) / (la * la), dot(d, b) / (lb * lb));
            let back = add_scaled(add_scaled(center, c, a), s, b);
            if norm(sub(back, p)) > ON_ARC * la.max(lb) || ((c * c + s * s).sqrt() - 1.0).abs() > ON_ARC {
                return Err(Error::Invalid(format!("corner {p:?} is not on the edge's ellipse or arc about {center:?}")));
            }
            Ok(crate::vec3::atan2(s, c))
        };
        let start = angle(from)?;
        let mut sweep = angle(to)? - start;
        if sweep > std::f64::consts::PI {
            sweep -= 2.0 * std::f64::consts::PI;
        } else if sweep <= -std::f64::consts::PI {
            sweep += 2.0 * std::f64::consts::PI;
        }
        if !(sweep.abs() > MIN_SWEEP) || std::f64::consts::PI - sweep.abs() < HALF_TURN_SLACK {
            return Err(Error::Invalid(format!("the arc from {from:?} to {to:?} spans 0 or 180 degrees; split it")));
        }
        Ok(Path::Ellipse { center, a, b, start, sweep, ends: [from, to] })
    }

    fn point(&self, s: f64) -> P {
        match *self {
            // The ends are the corners exactly (not to rounding): shared faces and the
            // transfinite grids, which need equal corners, rely on it.
            Path::Line { from, .. } | Path::Ellipse { ends: [from, _], .. } if s == 0.0 => from,
            Path::Line { to, .. } | Path::Ellipse { ends: [_, to], .. } if s == 1.0 => to,
            Path::Line { from, to } => add_scaled(from, s, sub(to, from)),
            Path::Ellipse { center, a, b, start, sweep, .. } => {
                let t = start + s * sweep;
                add_scaled(add_scaled(center, libm::cos(t), a), libm::sin(t), b)
            }
        }
    }
}

/// The parameters of `n` elements along one direction with last/first size ratio `bias`, at
/// `order` points per element (`order · n + 1` values from 0 to 1; mid-side values halfway).
fn parameters(n: usize, bias: f64, order: usize) -> Vec<f64> {
    let q = if n > 1 { libm::pow(bias, 1.0 / (n - 1) as f64) } else { 1.0 };
    let mut ends = Vec::with_capacity(n + 1);
    let (mut sum, mut size) = (0.0, 1.0);
    ends.push(0.0);
    for _ in 0..n {
        sum += size;
        ends.push(sum);
        size *= q;
    }
    let ends: Vec<f64> = ends.iter().map(|e| e / sum).collect();
    let mut out = Vec::with_capacity(order * n + 1);
    for w in ends.windows(2) {
        let (a, b) = (w[0], w[1]);
        out.push(a);
        if order == 2 {
            out.push(0.5 * (a + b));
        }
    }
    out.push(1.0);
    out
}

/// The block's lattice of `params.len()` points per direction (node `(i, j, k)` at
/// `(k·n1 + j)·n0 + i`): the twelve edge curves at their parameters, the six faces as transfinite
/// grids of their edges ([`transfinite_grid`]) and the interior by the trivariate interpolation of
/// the faces ([`transfinite_volume`]), blended by the lattice parameters.
fn lattice_points(paths: &[Path; 12], params: &[Vec<f64>; 3]) -> Result<Vec<P>> {
    let n = [0, 1, 2].map(|d| params[d].len() - 1);
    let mut x: Vec<Option<P>> = vec![None; params.iter().map(Vec::len).product()];
    for (k, [from, _]) in EDGES.iter().enumerate() {
        let along = k / 4;
        for (q, &t) in params[along].iter().enumerate() {
            let c = [0, 1, 2].map(|d| if d == along { q } else { CORNERS[*from][d] * n[d] });
            if let Some(slot) = x.get_mut(lattice_index(n, c)) {
                slot.get_or_insert(paths[k].point(t));
            }
        }
    }
    let blend = [params[0].as_slice(), params[1].as_slice(), params[2].as_slice()];
    transfinite_faces(&mut x, n, Some(blend), |_, _, _| Ok(()))?;
    transfinite_volume(&mut x, n, blend)?;
    x.into_iter().map(|p| p.ok_or_else(|| Error::Failed("a block lattice node is missing".into()))).collect()
}

/// Lattice offsets of each element of one cell, in units of the lattice (two per cell for
/// quadratic kinds), in the model's node order.
fn cell_elements(kind: ElementKind) -> Vec<Vec<[usize; 3]>> {
    // Corners in units of the lattice, then (quadratic kinds) a mid-side node on each edge in
    // the kind's edge order.
    let nodes = |corners: &[[usize; 3]]| -> Vec<[usize; 3]> {
        if !kind.is_quadratic() {
            return corners.to_vec();
        }
        let mut out: Vec<[usize; 3]> = corners.iter().map(|c| c.map(|x| 2 * x)).collect();
        for edge in kind.edges() {
            if let (Some(a), Some(b)) = (edge.first().and_then(|&i| corners.get(i)), edge.get(1).and_then(|&i| corners.get(i))) {
                out.push([a[0] + b[0], a[1] + b[1], a[2] + b[2]]);
            }
        }
        out
    };
    match kind {
        ElementKind::Hex8 | ElementKind::Hex20 => vec![nodes(&CORNERS)],
        ElementKind::Tet4 | ElementKind::Tet10 => {
            [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]
                .iter()
                .map(|perm| {
                    let mut v = [[0usize; 3]; 4];
                    for s in 0..3 {
                        v[s + 1] = v[s];
                        v[s + 1][perm[s]] = 1;
                    }
                    // Even permutations are positively oriented; swap two nodes of the others.
                    let odd = matches!(perm, [0, 2, 1] | [1, 0, 2] | [2, 1, 0]);
                    if odd {
                        v.swap(1, 2);
                    }
                    nodes(&v)
                })
                .collect()
        }
        ElementKind::Wedge6 => vec![
            vec![[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1]],
            vec![[0, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1], [0, 1, 1]],
        ],
        _ => Vec::new(),
    }
}

/// Mesh one block.
pub fn block(b: &Block) -> Result<BlockMesh> {
    let kind = b.kind;
    if !matches!(kind, ElementKind::Hex8 | ElementKind::Hex20 | ElementKind::Tet4 | ElementKind::Tet10 | ElementKind::Wedge6) {
        return Err(Error::Invalid(format!("a block is meshed with hex8, hex20, tet4, tet10 or wedge6 elements, not {}", kind.name())));
    }
    if b.corners.iter().flatten().any(|c| !c.is_finite() || c.abs() > MAX_COORD) {
        return Err(Error::Invalid(format!("block corners must be finite and within ±{MAX_COORD}")));
    }
    if b.divisions.iter().any(|&n| n == 0 || n > MAX_DIVISIONS) {
        return Err(Error::Invalid(format!("divisions must be 1 to {MAX_DIVISIONS} in each direction")));
    }
    if b.bias.iter().any(|&r| !(r.is_finite() && (MIN_BIAS..=MAX_BIAS).contains(&r))) {
        return Err(Error::Invalid("bias must be between 1e-6 and 1e6 (1 is uniform)".into()));
    }
    let per_cell = cell_elements(kind).len();
    let cells = b.divisions.iter().try_fold(1usize, |acc, &n| acc.checked_mul(n));
    let count = cells.and_then(|c| c.checked_mul(per_cell)).filter(|&n| n <= MAX_BLOCK_ELEMENTS).ok_or_else(|| {
        Error::Limit(format!("a block makes at most {MAX_BLOCK_ELEMENTS} elements; {:?} divisions of {} are more", b.divisions, kind.name()))
    })?;

    let mut curves: [Option<Curve>; 12] = [None; 12];
    for (pair, curve) in &b.curves {
        let k = Block::edge_index(*pair)
            .ok_or_else(|| Error::Invalid(format!("corners {} and {} are not joined by a block edge", pair[0] + 1, pair[1] + 1)))?;
        if curves[k].replace(*curve).is_some() {
            return Err(Error::Invalid(format!("the edge from corner {} to {} is given twice", pair[0] + 1, pair[1] + 1)));
        }
    }
    let mut paths = [Path::Line { from: [0.0; 3], to: [0.0; 3] }; 12];
    for (k, [from, to]) in EDGES.iter().enumerate() {
        paths[k] = Path::new(curves[k].unwrap_or(Curve::Line), b.corners[*from], b.corners[*to])?;
    }

    let order = if kind.is_quadratic() { 2 } else { 1 };
    let params: [Vec<f64>; 3] = [0, 1, 2].map(|d| parameters(b.divisions[d], b.bias[d], order));
    let lattice_xyz = lattice_points(&paths, &params)?;
    let size = [0, 1, 2].map(|d| params[d].len());
    let index = |c: [usize; 3]| c[0] + size[0] * (c[1] + size[1] * c[2]);

    // Elements as lattice coordinates; then number the lattice points they use, in lattice order.
    let pattern = cell_elements(kind);
    let mut lattice: Vec<Vec<[usize; 3]>> = Vec::with_capacity(count);
    for k in 0..b.divisions[2] {
        for j in 0..b.divisions[1] {
            for i in 0..b.divisions[0] {
                let base = [i * order, j * order, k * order];
                for e in &pattern {
                    lattice.push(e.iter().map(|o| [base[0] + o[0], base[1] + o[1], base[2] + o[2]]).collect());
                }
            }
        }
    }
    let mut used = vec![false; size.iter().product::<usize>()];
    for c in lattice.iter().flatten() {
        if let Some(slot) = used.get_mut(index(*c)) {
            *slot = true;
        }
    }
    let mut number: Vec<Option<u32>> = vec![None; used.len()];
    let mut points = Vec::new();
    for k in 0..size[2] {
        for j in 0..size[1] {
            for i in 0..size[0] {
                let at = index([i, j, k]);
                if used.get(at) == Some(&true)
                    && let Some(slot) = number.get_mut(at)
                {
                    *slot = Some(u32::try_from(points.len()).map_err(|_| Error::Limit("too many points in one block".into()))?);
                    points.push(lattice_xyz.get(at).copied().unwrap_or([f64::NAN; 3]));
                }
            }
        }
    }
    let point_of = |c: [usize; 3]| number.get(index(c)).copied().flatten();
    let last = [0, 1, 2].map(|d| size[d] - 1);
    let mut corners = [0u32; 8];
    for (slot, c) in corners.iter_mut().zip(CORNERS) {
        *slot = point_of([c[0] * last[0], c[1] * last[1], c[2] * last[2]]).ok_or_else(|| Error::Failed("a block corner has no point".into()))?;
    }

    // Block faces in hex8 face order: (direction, at the end?).
    const FACES: [(usize, bool); 6] = [(2, false), (2, true), (1, false), (0, true), (1, true), (0, false)];
    let on_face = |c: &[usize; 3], f: (usize, bool)| c[f.0] == if f.1 { last[f.0] } else { 0 };
    // Each block edge: the direction it runs along and the lattice coordinates fixed on it.
    let block_edges: [(usize, [usize; 3]); 12] = std::array::from_fn(|be| {
        let [from, _] = EDGES[be];
        (be / 4, [0, 1, 2].map(|d| CORNERS[from][d] * last[d]))
    });
    let corner_count = kind.corner_count();
    let mut elements = Vec::with_capacity(lattice.len());
    let mut faces = Vec::new();
    let mut edges = Vec::new();
    for (e, nodes) in lattice.iter().enumerate() {
        let label = u32::try_from(e).map_err(|_| Error::Limit("too many elements in one block".into()))?;
        let ids: Vec<u32> =
            nodes.iter().map(|c| point_of(*c)).collect::<Option<_>>().ok_or_else(|| Error::Failed("an element node has no point".into()))?;
        let xyz: Vec<P> = ids.iter().map(|&i| points.get(i as usize).copied().unwrap_or([f64::NAN; 3])).collect();
        let check = element_validity(kind, &xyz)?;
        if check.validity != Validity::Valid {
            let cell = nodes.first().map(|c| [c[0] / order, c[1] / order, c[2] / order]).unwrap_or_default();
            return Err(Error::Invalid(format!(
                "the block's {} element in cell {cell:?} is inverted or degenerate: give the corners in hex8 order (1-2-3-4 counter-clockwise seen from 5-6-7-8) and curves that do not fold the block",
                kind.name()
            )));
        }
        for (f, local) in kind.faces().iter().enumerate() {
            let cs: Vec<&[usize; 3]> = local.iter().filter(|&&i| i < corner_count).filter_map(|&i| nodes.get(i)).collect();
            for (bf, side) in FACES.iter().enumerate() {
                if cs.iter().all(|c| on_face(c, *side)) {
                    faces.push((label, f as u8 + 1, bf as u8));
                }
            }
        }
        for (k, local) in kind.edges().iter().enumerate() {
            let (Some(a), Some(c)) = (local.first().and_then(|&i| nodes.get(i)), local.get(1).and_then(|&i| nodes.get(i))) else { continue };
            for (be, (along, fixed)) in block_edges.iter().enumerate() {
                let on = |p: &[usize; 3]| (0..3).filter(|d| d != along).all(|d| p[d] == fixed[d]);
                if on(a) && on(c) {
                    edges.push((label, k as u8 + 1, be as u8));
                }
            }
        }
        elements.push(ids);
    }
    Ok(BlockMesh { points, elements, faces, edges, corners })
}

/// The corners of a patch's four edges, as [`EDGES`] for blocks: two along 1→2 (at the start
/// and the end of 1→4), then two along 1→4.
pub const PATCH_EDGES: [[usize; 2]; 4] = [[0, 1], [3, 2], [0, 3], [1, 2]];

/// One four-sided planar patch to mesh with 2D elements (plane stress, plane strain and
/// axisymmetric sections live in the XY plane, so its corners and curves are at z = 0).
#[derive(Debug, Clone, PartialEq)]
pub struct Patch {
    /// quad4, quad8, tri3 or tri6.
    pub kind: ElementKind,
    /// The corners counter-clockwise seen from +z.
    pub corners: [P; 4],
    /// Elements along 1→2 and 1→4.
    pub divisions: [usize; 2],
    pub bias: [f64; 2],
    pub curves: Vec<([usize; 2], Curve)>,
}

impl Patch {
    /// The patch edge ([`PATCH_EDGES`] index) joining two 0-based corners, either way round.
    pub fn edge_index(corners: [usize; 2]) -> Option<usize> {
        PATCH_EDGES.iter().position(|e| *e == corners || *e == [corners[1], corners[0]])
    }
}

/// Lattice offsets of each 2D element of one cell, in the model's node order.
fn cell_elements_2d(kind: ElementKind) -> Vec<Vec<[usize; 2]>> {
    let nodes = |corners: &[[usize; 2]]| -> Vec<[usize; 2]> {
        if !kind.is_quadratic() {
            return corners.to_vec();
        }
        let mut out: Vec<[usize; 2]> = corners.iter().map(|c| c.map(|x| 2 * x)).collect();
        for edge in kind.edges() {
            if let (Some(a), Some(b)) = (edge.first().and_then(|&i| corners.get(i)), edge.get(1).and_then(|&i| corners.get(i))) {
                out.push([a[0] + b[0], a[1] + b[1]]);
            }
        }
        out
    };
    match kind {
        ElementKind::Quad4 | ElementKind::Quad8 => vec![nodes(&[[0, 0], [1, 0], [1, 1], [0, 1]])],
        ElementKind::Tri3 | ElementKind::Tri6 => vec![nodes(&[[0, 0], [1, 0], [1, 1]]), nodes(&[[0, 0], [1, 1], [0, 1]])],
        _ => Vec::new(),
    }
}

/// Mesh one patch with 2D elements by the Coons (bilinearly blended) transfinite interpolation of
/// its four edges. Each element's corners must turn counter-clockwise (exact predicates); mid-side
/// nodes come from the map, so quadratic elements follow the curved edges. In the result,
/// `faces` is empty and `edges` lists element edges on the patch edges ([`PATCH_EDGES`]).
pub fn patch(b: &Patch) -> Result<BlockMesh> {
    let kind = b.kind;
    if !matches!(kind, ElementKind::Quad4 | ElementKind::Quad8 | ElementKind::Tri3 | ElementKind::Tri6) {
        return Err(Error::Invalid(format!("a patch is meshed with quad4, quad8, tri3 or tri6 elements, not {}", kind.name())));
    }
    if b.corners.iter().flatten().any(|c| !c.is_finite() || c.abs() > MAX_COORD) {
        return Err(Error::Invalid(format!("patch corners must be finite and within ±{MAX_COORD}")));
    }
    if b.corners.iter().any(|c| c[2] != 0.0) {
        return Err(Error::Invalid("2D elements lie in the XY plane: give the patch corners with z = 0".into()));
    }
    if b.divisions.iter().any(|&n| n == 0 || n > MAX_DIVISIONS) {
        return Err(Error::Invalid(format!("divisions must be 1 to {MAX_DIVISIONS} in each direction")));
    }
    if b.bias.iter().any(|&r| !(r.is_finite() && (MIN_BIAS..=MAX_BIAS).contains(&r))) {
        return Err(Error::Invalid("bias must be between 1e-6 and 1e6 (1 is uniform)".into()));
    }
    let pattern = cell_elements_2d(kind);
    b.divisions[0]
        .checked_mul(b.divisions[1])
        .and_then(|c| c.checked_mul(pattern.len()))
        .filter(|&n| n <= MAX_BLOCK_ELEMENTS)
        .ok_or_else(|| Error::Limit(format!("a patch makes at most {MAX_BLOCK_ELEMENTS} elements; {:?} divisions are more", b.divisions)))?;
    let mut curves: [Option<Curve>; 4] = [None; 4];
    for (pair, curve) in &b.curves {
        let k = Patch::edge_index(*pair)
            .ok_or_else(|| Error::Invalid(format!("corners {} and {} are not joined by a patch edge", pair[0] + 1, pair[1] + 1)))?;
        if let Curve::Arc { center } | Curve::Ellipse { center, .. } = curve
            && center[2] != 0.0
        {
            return Err(Error::Invalid("a patch's curves lie in the XY plane: give their centres with z = 0".into()));
        }
        if curves[k].replace(*curve).is_some() {
            return Err(Error::Invalid(format!("the edge from corner {} to {} is given twice", pair[0] + 1, pair[1] + 1)));
        }
    }
    let order = if kind.is_quadratic() { 2 } else { 1 };
    let params: [Vec<f64>; 2] = [0, 1].map(|d| parameters(b.divisions[d], b.bias[d], order));
    let mut sides: [Vec<P>; 4] = Default::default();
    for (k, [from, to]) in PATCH_EDGES.iter().enumerate() {
        let path = Path::new(curves[k].unwrap_or(Curve::Line), b.corners[*from], b.corners[*to])?;
        sides[k] = params[k / 2].iter().map(|&t| path.point(t)).collect();
    }
    // Bottom (1→2), right (2→3), top (4→3), left (1→4): the Coons patch of the four sides.
    let [bottom, top, left, right] = &sides;
    let grid = transfinite_grid_blended(bottom, right, top, left, &params[0], &params[1])?;
    let at = |i: usize, j: usize| -> P {
        let mut x = grid.get(j * params[0].len() + i).copied().unwrap_or([f64::NAN; 3]);
        x[2] = 0.0;
        x
    };
    let size = [params[0].len(), params[1].len()];
    let mut used = vec![false; size[0] * size[1]];
    let mut lattice: Vec<Vec<[usize; 2]>> = Vec::new();
    for j in 0..b.divisions[1] {
        for i in 0..b.divisions[0] {
            for e in &pattern {
                let nodes: Vec<[usize; 2]> = e.iter().map(|o| [i * order + o[0], j * order + o[1]]).collect();
                for c in &nodes {
                    if let Some(slot) = used.get_mut(c[0] + size[0] * c[1]) {
                        *slot = true;
                    }
                }
                lattice.push(nodes);
            }
        }
    }
    let mut number: Vec<Option<u32>> = vec![None; used.len()];
    let mut points = Vec::new();
    for j in 0..size[1] {
        for i in 0..size[0] {
            let k = i + size[0] * j;
            if used.get(k) == Some(&true)
                && let Some(slot) = number.get_mut(k)
            {
                *slot = Some(u32::try_from(points.len()).map_err(|_| Error::Limit("too many points in one patch".into()))?);
                points.push(at(i, j));
            }
        }
    }
    let point_of = |c: [usize; 2]| number.get(c[0] + size[0] * c[1]).copied().flatten();
    let last = [size[0] - 1, size[1] - 1];
    let mut corners = [0u32; 8];
    for (slot, c) in corners.iter_mut().zip([[0, 0], [last[0], 0], [last[0], last[1]], [0, last[1]]]) {
        *slot = point_of(c).ok_or_else(|| Error::Failed("a patch corner has no point".into()))?;
    }
    // Each patch edge: the direction it runs along and the lattice coordinate fixed on it.
    let patch_edges: [(usize, usize); 4] = [(0, 0), (0, last[1]), (1, 0), (1, last[0])];
    let mut elements = Vec::with_capacity(lattice.len());
    let mut on_edges = Vec::new();
    for (e, nodes) in lattice.iter().enumerate() {
        let label = u32::try_from(e).map_err(|_| Error::Limit("too many elements in one patch".into()))?;
        let ids: Vec<u32> =
            nodes.iter().map(|c| point_of(*c)).collect::<Option<_>>().ok_or_else(|| Error::Failed("an element node has no point".into()))?;
        let xy: Vec<[f64; 2]> =
            ids.iter().take(kind.corner_count()).map(|&i| points.get(i as usize).map_or([f64::NAN; 2], |p| [p[0], p[1]])).collect();
        let valid = match xy.as_slice() {
            [a, b, c, d] => planar_quad_validity(&[*a, *b, *c, *d])? == Validity::Valid,
            [a, b, c] => (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) > 0.0,
            _ => false,
        };
        if !valid {
            let cell = nodes.first().map(|c| [c[0] / order, c[1] / order]).unwrap_or_default();
            return Err(Error::Invalid(format!(
                "the patch's {} element in cell {cell:?} is inverted or degenerate: give the corners counter-clockwise seen from +z and curves that do not fold the patch",
                kind.name()
            )));
        }
        for (k, local) in kind.edges().iter().enumerate() {
            let (Some(a), Some(c)) = (local.first().and_then(|&i| nodes.get(i)), local.get(1).and_then(|&i| nodes.get(i))) else { continue };
            for (pe, (along, fixed)) in patch_edges.iter().enumerate() {
                let other = 1 - along;
                if a[other] == *fixed && c[other] == *fixed {
                    on_edges.push((label, k as u8 + 1, pe as u8));
                }
            }
        }
        elements.push(ids);
    }
    Ok(BlockMesh { points, elements, faces: Vec::new(), edges: on_edges, corners })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit(kind: ElementKind, divisions: [usize; 3]) -> Block {
        let corners = CORNERS.map(|c| [c[0] as f64, c[1] as f64, c[2] as f64]);
        Block { kind, corners, divisions, bias: [1.0; 3], curves: Vec::new() }
    }

    fn volume(kind: ElementKind, m: &BlockMesh) -> f64 {
        m.elements
            .iter()
            .map(|e| {
                let xyz: Vec<P> = e.iter().map(|&i| m.points[i as usize]).collect();
                element_validity(kind, &xyz).unwrap().volume
            })
            .sum()
    }

    #[test]
    fn every_kind_fills_a_cube_with_valid_elements_and_shared_nodes() {
        for (kind, per_cell, points) in [
            (ElementKind::Hex8, 1, 27),
            (ElementKind::Hex20, 1, 81),
            (ElementKind::Tet4, 6, 27),
            (ElementKind::Tet10, 6, 125),
            (ElementKind::Wedge6, 2, 27),
        ] {
            let m = block(&unit(kind, [2, 2, 2])).unwrap();
            assert_eq!(m.elements.len(), 8 * per_cell, "{kind:?}");
            assert_eq!(m.points.len(), points, "{kind:?}");
            assert!((volume(kind, &m) - 1.0).abs() < 1e-12, "{kind:?}");
            // Each block face is covered by element faces of total area 1: count them instead.
            let per_face = |bf: u8| m.faces.iter().filter(|f| f.2 == bf).count();
            let expect = if matches!(kind, ElementKind::Hex8 | ElementKind::Hex20) { 4 } else { 8 };
            for bf in 0..6 {
                let n = per_face(bf);
                // Wedges have quads on the sides and triangles top and bottom.
                let want = if kind == ElementKind::Wedge6 && bf >= 2 { 4 } else { expect };
                assert_eq!(n, want, "{kind:?} face {bf}");
            }
            for be in 0..12 {
                // Two segments per block edge (an edge segment may belong to several elements).
                let mut segments: Vec<[u32; 2]> = m
                    .edges
                    .iter()
                    .filter(|e| e.2 == be)
                    .map(|&(e, k, _)| {
                        let local = kind.edges()[k as usize - 1];
                        let el = &m.elements[e as usize];
                        let mut s = [el[local[0]], el[local[1]]];
                        s.sort();
                        s
                    })
                    .collect();
                segments.sort();
                segments.dedup();
                assert_eq!(segments.len(), 2, "{kind:?} edge {be}");
            }
            assert_eq!(m.points[m.corners[6] as usize], [1.0, 1.0, 1.0]);
        }
    }

    #[test]
    fn arcs_and_ellipses_put_every_edge_node_on_the_curve() {
        // A quarter annulus a = 1, b = 2, extruded 0.5, hex20 so mid-side nodes are on the arcs.
        let (a, b) = (1.0, 2.0);
        let corners = [[a, 0.0, 0.0], [b, 0.0, 0.0], [0.0, b, 0.0], [0.0, a, 0.0], [a, 0.0, 0.5], [b, 0.0, 0.5], [0.0, b, 0.5], [0.0, a, 0.5]];
        let o = |z: f64| [0.0, 0.0, z];
        let curves = vec![
            ([1, 2], Curve::Arc { center: o(0.0) }),
            ([0, 3], Curve::Arc { center: o(0.0) }),
            ([5, 6], Curve::Arc { center: o(0.5) }),
            ([4, 7], Curve::Ellipse { center: o(0.5), axes: [[a, 0.0, 0.0], [0.0, a, 0.0]] }),
        ];
        let m = block(&Block { kind: ElementKind::Hex20, corners, divisions: [3, 6, 2], bias: [1.0; 3], curves }).unwrap();
        for p in &m.points {
            let r = libm::hypot(p[0], p[1]);
            assert!(r > a - 1e-12 && r < b + 1e-12, "{p:?}");
        }
        let exact = std::f64::consts::PI / 4.0 * (b * b - a * a) * 0.5;
        let v = volume(ElementKind::Hex20, &m);
        assert!((v - exact).abs() < 1e-4 * exact, "{v} vs {exact}");
        // The inner face's nodes are on r = a.
        let elements = &m.elements;
        let inner: Vec<u32> = m
            .faces
            .iter()
            .filter(|f| f.2 == 5)
            .flat_map(|&(e, f, _)| ElementKind::Hex20.faces()[f as usize - 1].iter().map(move |&i| elements[e as usize][i]))
            .collect();
        assert!(!inner.is_empty());
        assert!(inner.iter().all(|&i| (m.points[i as usize][0].hypot(m.points[i as usize][1]) - a).abs() < 1e-12));
    }

    #[test]
    fn patches_mesh_a_quarter_annulus_with_every_2d_kind() {
        let (a, b) = (1.0, 2.0);
        let o = [0.0; 3];
        for (kind, per_cell) in [(ElementKind::Quad4, 1), (ElementKind::Quad8, 1), (ElementKind::Tri3, 2), (ElementKind::Tri6, 2)] {
            let p = Patch {
                kind,
                corners: [[a, 0.0, 0.0], [b, 0.0, 0.0], [0.0, b, 0.0], [0.0, a, 0.0]],
                divisions: [3, 6],
                bias: [1.0; 2],
                curves: vec![([1, 2], Curve::Arc { center: o }), ([0, 3], Curve::Arc { center: o })],
            };
            let m = patch(&p).unwrap();
            assert_eq!(m.elements.len(), 18 * per_cell, "{kind:?}");
            assert!(m.points.iter().all(|q| q[2] == 0.0 && q[0].hypot(q[1]) > a - 1e-12 && q[0].hypot(q[1]) < b + 1e-12));
            // Element edges on each patch edge: 3, 3, 6, 6 segments.
            for (pe, n) in [(0, 3), (1, 3), (2, 6), (3, 6)] {
                assert_eq!(m.edges.iter().filter(|e| e.2 == pe).count(), n, "{kind:?} edge {pe}");
            }
            if kind.is_quadratic() {
                // Mid-side nodes of the inner arc lie on it.
                let on_inner = m.edges.iter().filter(|e| e.2 == 2).all(|&(e, k, _)| {
                    let mid = kind.edges()[k as usize - 1][2];
                    let q = m.points[m.elements[e as usize][mid] as usize];
                    (q[0].hypot(q[1]) - a).abs() < 1e-12
                });
                assert!(on_inner, "{kind:?}");
            }
        }
        let mut clockwise = Patch {
            kind: ElementKind::Quad4,
            corners: [[0.0; 3], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0]],
            divisions: [1, 1],
            bias: [1.0; 2],
            curves: vec![],
        };
        assert!(patch(&clockwise).is_err());
        clockwise.corners.swap(1, 3);
        assert!(patch(&clockwise).is_ok());
        clockwise.corners[0][2] = 1.0;
        assert!(patch(&clockwise).is_err(), "off the XY plane");
        clockwise.corners[0][2] = 0.0;
        clockwise.kind = ElementKind::Hex8;
        assert!(patch(&clockwise).is_err());
    }

    #[test]
    fn grading_and_bad_blocks() {
        let mut g = unit(ElementKind::Hex8, [4, 1, 1]);
        g.bias = [8.0, 1.0, 1.0];
        let m = block(&g).unwrap();
        let xs: Vec<f64> = (0..5).map(|i| m.points[i][0]).collect();
        let (first, last) = (xs[1] - xs[0], xs[4] - xs[3]);
        assert!((last / first - 8.0).abs() < 1e-9, "{xs:?}");

        let mut inverted = unit(ElementKind::Hex8, [1, 1, 1]);
        inverted.corners.swap(1, 3);
        assert!(block(&inverted).is_err());
        let mut off = unit(ElementKind::Hex8, [1, 1, 1]);
        off.curves = vec![([0, 1], Curve::Arc { center: [0.3, 0.0, 0.0] })];
        assert!(block(&off).is_err(), "corners at different distances from the centre");
        let mut half = unit(ElementKind::Hex8, [1, 1, 1]);
        half.curves = vec![([0, 1], Curve::Arc { center: [0.5, 0.0, 0.0] })];
        assert!(block(&half).is_err(), "180 degrees has no plane");
        let mut diagonal = unit(ElementKind::Hex8, [1, 1, 1]);
        diagonal.curves = vec![([0, 6], Curve::Line)];
        assert!(block(&diagonal).is_err(), "not an edge");
        for divisions in [[0, 1, 1], [MAX_DIVISIONS + 1, 1, 1], [MAX_DIVISIONS, MAX_DIVISIONS, MAX_DIVISIONS]] {
            assert!(block(&unit(ElementKind::Hex8, divisions)).is_err(), "{divisions:?}");
        }
        assert!(block(&unit(ElementKind::Quad4, [1, 1, 1])).is_err());
        let mut nan = unit(ElementKind::Hex8, [1, 1, 1]);
        nan.corners[0][0] = f64::NAN;
        assert!(block(&nan).is_err());
        let mut bias = unit(ElementKind::Hex8, [2, 1, 1]);
        bias.bias[0] = 0.0;
        assert!(block(&bias).is_err());
    }
}
