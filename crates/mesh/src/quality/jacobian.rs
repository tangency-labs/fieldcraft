//! The Jacobian determinant of a solid element as an exact polynomial, its certified sign and
//! the element's exact volume.
//!
//! The geometry map of every solid element FieldCraft knows is a polynomial on its reference
//! domain, so its Jacobian determinant is one too: constant for tet4, cubic for tet10, degree
//! (2, 2) on triangle × segment for wedge6, degree (2, 2, 2) for hex8 and (5, 5, 5) for hex20. It
//! is built here in Bernstein form ([`super::bernstein`]) and its sign is settled by bounds and
//! subdivision, the method of Johnen, Remacle and Geuzaine (J. Comput. Phys. 233, 2013).
//!
//! Positive corner Jacobians do not make a trilinear hexahedron valid: its determinant has
//! degree 2 in each variable and can be negative inside while positive at all eight corners
//! (Marschner, Palmer, Zhang, Solomon, *Hexahedral mesh repair via sum-of-squares relaxation*,
//! CGF 39(5), 2020; Knupp, *Achieving finite element mesh quality via optimization of the
//! Jacobian matrix norm and associated quantities*, IJNME 48, 2000). Sampling at Gauss points is
//! no proof either. The check here is a proof: `Valid` means `det J > 0` on the whole element.
//!
//! Node orders (corners first, then mid-side nodes):
//! - tet4/tet10: mid-side nodes 5–10 on edges 1-2, 2-3, 3-1, 1-4, 2-4, 3-4; positive when node 4
//!   lies on the side of face 1-2-3 that its right-hand normal points to.
//! - wedge6: triangle 1-2-3, then 4-5-6 above it; positive when 4-5-6 lies on the side of 1-2-3's
//!   right-hand normal.
//! - pyramid5: base 1-2-3-4, then the apex 5 on the side of its right-hand normal. Its map is
//!   rational in the pyramid's coordinates (ξ, η, ζ), but its Jacobian determinant is not: with
//!   the base's bilinear map `B(s, t)` on the unit square (ξ = (2s − 1)(1 − ζ), likewise η),
//!   the columns are `∂x/∂ξ = ½ B_s`, `∂x/∂η = ½ B_t` and `∂x/∂ζ = x₅ − B + (s − ½) B_s +
//!   (t − ½) B_t`, all bilinear in (s, t) and constant along each line to the apex, so the
//!   determinant is a polynomial of degree (3, 3) on the square and the element's volume is
//!   4/3 of its integral there.
//! - hex8/hex20: face 1-2-3-4, then 5-6-7-8 above it; mid-side nodes 9–12 on the edges of 1-2-3-4,
//!   13–16 on those of 5-6-7-8, 17–20 on edges 1-5, 2-6, 3-7, 4-8.

use fieldcraft_model::ElementKind;

use super::bernstein::{Bernstein, Factor};
use crate::{Error, Result};

/// Reference coordinates of the hex corners on the unit cube, in node order.
pub(crate) const HEX_CORNERS: [[usize; 3]; 8] = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
/// Hex edges in the order of hex20's mid-side nodes.
pub(crate) const HEX_EDGES: [(usize, usize); 12] = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)];
/// Tet edges in the order of tet10's mid-side nodes.
pub const TET_EDGES: [(usize, usize); 6] = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];

/// A Jacobian determinant is positive when it exceeds this fraction of its largest magnitude.
/// Anything smaller is degenerate to rounding and counts as invalid.
pub const RELATIVE_TOLERANCE: f64 = 1e-12; // unit-factors: relative to the largest magnitude
/// Subdivision caps for the certificate: depth of one piece and pieces per element.
pub const MAX_DEPTH: u32 = 30;
pub const MAX_PIECES: usize = 4096;

/// Whether an element's Jacobian determinant is positive everywhere.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Validity {
    /// Proven positive on the whole element.
    Valid,
    /// Zero or negative somewhere (degenerate or inverted): a point where it is was found.
    Invalid,
    /// The subdivision caps were reached without a proof either way. Treat as invalid.
    Undetermined,
}

/// The certified Jacobian check of one solid element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JacobianCheck {
    pub validity: Validity,
    /// A lower bound of `min det J` divided by an upper bound of `max |det J|`, in `[-1, 1]`;
    /// 1 for an affine element (straight-sided tet, parallelepiped). For a valid element it comes
    /// from the subdivision that proved it, so it is conservative.
    pub ratio: f64,
    /// The exact signed volume (`∫ det J` over the reference domain).
    pub volume: f64,
}

/// The geometry map's Jacobian of one element: its columns (each three coordinate polynomials)
/// and its determinant.
pub(crate) struct Jacobian {
    pub columns: [[Bernstein; 3]; 3],
    pub det: Bernstein,
}

fn domain(kind: ElementKind) -> Result<(Vec<Factor>, [(usize, usize); 3])> {
    let seg = |deg| Factor { dim: 1, deg };
    Ok(match kind {
        ElementKind::Tet4 => (vec![Factor { dim: 3, deg: 1 }], [(0, 1), (0, 2), (0, 3)]),
        ElementKind::Tet10 => (vec![Factor { dim: 3, deg: 2 }], [(0, 1), (0, 2), (0, 3)]),
        ElementKind::Wedge6 => (vec![Factor { dim: 2, deg: 1 }, seg(1)], [(0, 1), (0, 2), (1, 1)]),
        ElementKind::Hex8 => (vec![seg(1); 3], [(0, 1), (1, 1), (2, 1)]),
        ElementKind::Hex20 => (vec![seg(2); 3], [(0, 1), (1, 1), (2, 1)]),
        other => return Err(Error::Invalid(format!("{} is not a solid element", other.name()))),
    })
}

/// Free components of a simplex vertex `v` scaled by `deg` (vertex 0 is all zeros).
fn simplex_vertex(v: usize, deg: usize) -> [usize; 3] {
    let mut m = [0; 3];
    if v > 0
        && let Some(x) = m.get_mut(v - 1)
    {
        *x = deg;
    }
    m
}

fn mid(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])]
}

/// The quadratic Bernstein control point between two end points for a curve through `m` at the
/// parameter midpoint: `2m − (a + b)/2`.
fn quadratic_control(a: [f64; 3], m: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [2.0 * m[0] - 0.5 * (a[0] + b[0]), 2.0 * m[1] - 0.5 * (a[1] + b[1]), 2.0 * m[2] - 0.5 * (a[2] + b[2])]
}

/// The 20-node serendipity hexahedron's shape functions at reference point `x ∈ [-1, 1]³`.
fn serendipity(x: [f64; 3], xyz: &[[f64; 3]]) -> [f64; 3] {
    let mut out = [0.0; 3];
    let mut add = |n: f64, p: Option<&[f64; 3]>| {
        if let Some(p) = p {
            for (o, c) in out.iter_mut().zip(p) {
                *o += n * c;
            }
        }
    };
    let corner = |k: usize| HEX_CORNERS.get(k).map(|c| c.map(|v| 2.0 * v as f64 - 1.0)).unwrap_or([0.0; 3]);
    for k in 0..8 {
        let c = corner(k);
        let n = 0.125 * (1.0 + x[0] * c[0]) * (1.0 + x[1] * c[1]) * (1.0 + x[2] * c[2]) * (x[0] * c[0] + x[1] * c[1] + x[2] * c[2] - 2.0);
        add(n, xyz.get(k));
    }
    for (e, &(a, b)) in HEX_EDGES.iter().enumerate() {
        let r = mid(corner(a), corner(b));
        let mut n = 0.25;
        for (xi, ri) in x.iter().zip(r) {
            n *= if ri == 0.0 { 1.0 - xi * xi } else { 1.0 + xi * ri };
        }
        add(n, xyz.get(8 + e));
    }
    out
}

/// The geometry map in Bernstein form: one polynomial per coordinate.
fn geometry(kind: ElementKind, xyz: &[[f64; 3]]) -> Result<[Bernstein; 3]> {
    let (factors, _) = domain(kind)?;
    let zero = Bernstein::zero(&factors)?;
    let mut g = [zero.clone(), zero.clone(), zero];
    let mut put = |m: [usize; 3], p: [f64; 3]| {
        for (poly, c) in g.iter_mut().zip(p) {
            poly.set(m, c);
        }
    };
    let node = |i: usize| xyz.get(i).copied().ok_or_else(|| Error::Invalid(format!("{} needs {} nodes", kind.name(), kind.node_count())));
    match kind {
        ElementKind::Tet4 => {
            for v in 0..4 {
                put(simplex_vertex(v, 1), node(v)?);
            }
        }
        ElementKind::Tet10 => {
            for v in 0..4 {
                put(simplex_vertex(v, 2), node(v)?);
            }
            for (e, &(a, b)) in TET_EDGES.iter().enumerate() {
                let (ma, mb) = (simplex_vertex(a, 1), simplex_vertex(b, 1));
                put([ma[0] + mb[0], ma[1] + mb[1], ma[2] + mb[2]], quadratic_control(node(a)?, node(4 + e)?, node(b)?));
            }
        }
        ElementKind::Wedge6 => {
            for (i, base) in [[0, 0], [1, 0], [0, 1]].into_iter().enumerate() {
                put([base[0], base[1], 0], node(i)?);
                put([base[0], base[1], 1], node(i + 3)?);
            }
        }
        ElementKind::Hex8 => {
            for (i, c) in HEX_CORNERS.iter().enumerate() {
                put(*c, node(i)?);
            }
        }
        ElementKind::Hex20 => {
            if xyz.len() < 20 {
                return Err(Error::Invalid("hex20 needs 20 nodes".into()));
            }
            // The serendipity map lies in the triquadratic space, so interpolating it at the 27
            // triquadratic nodes reproduces it; then convert each axis from Lagrange to Bernstein.
            let mut v = [[[[0.0; 3]; 3]; 3]; 3];
            for (a, plane) in v.iter_mut().enumerate() {
                for (b, row) in plane.iter_mut().enumerate() {
                    for (c, p) in row.iter_mut().enumerate() {
                        *p = serendipity([a as f64 - 1.0, b as f64 - 1.0, c as f64 - 1.0], xyz);
                    }
                }
            }
            // One full pass per axis (the conversion is a tensor product).
            for axis in 0..3 {
                for i in 0..3 {
                    for j in 0..3 {
                        let at = |k: usize| match axis {
                            0 => (k, i, j),
                            1 => (i, k, j),
                            _ => (i, j, k),
                        };
                        let r = [0, 1, 2].map(|k| {
                            let (a, b, c) = at(k);
                            v[a][b][c]
                        });
                        let (a, b, c) = at(1);
                        v[a][b][c] = quadratic_control(r[0], r[1], r[2]);
                    }
                }
            }
            for (a, plane) in v.iter().enumerate() {
                for (b, row) in plane.iter().enumerate() {
                    for (c, p) in row.iter().enumerate() {
                        put([a, b, c], *p);
                    }
                }
            }
        }
        other => return Err(Error::Invalid(format!("{} is not a solid element", other.name()))),
    }
    Ok(g)
}

/// The pyramid's Jacobian columns (see the module docs) as bilinear polynomials on the unit
/// square of its base.
fn pyramid_columns(xyz: &[[f64; 3]]) -> Result<[[Bernstein; 3]; 3]> {
    let node = |i: usize| xyz.get(i).copied().ok_or_else(|| Error::Invalid("pyramid5 needs 5 nodes".into()));
    let [p0, p1, p2, p3, apex] = [node(0)?, node(1)?, node(2)?, node(3)?, node(4)?];
    let square = [Factor { dim: 1, deg: 1 }, Factor { dim: 1, deg: 1 }];
    let zero = Bernstein::zero(&square)?;
    let mut columns = [[zero.clone(), zero.clone(), zero.clone()], [zero.clone(), zero.clone(), zero.clone()], [zero.clone(), zero.clone(), zero]];
    // A bilinear polynomial's coefficients are its values at the square's corners.
    for (s, t) in [(0usize, 0usize), (1, 0), (1, 1), (0, 1)] {
        let (sf, tf) = (s as f64, t as f64);
        let base = |k: usize| p0[k] * (1.0 - sf) * (1.0 - tf) + p1[k] * sf * (1.0 - tf) + p2[k] * sf * tf + p3[k] * (1.0 - sf) * tf;
        let b_s = |k: usize| (p1[k] - p0[k]) * (1.0 - tf) + (p2[k] - p3[k]) * tf;
        let b_t = |k: usize| (p3[k] - p0[k]) * (1.0 - sf) + (p2[k] - p1[k]) * sf;
        for (k, top) in apex.iter().enumerate() {
            let values = [0.5 * b_s(k), 0.5 * b_t(k), top - base(k) + (sf - 0.5) * b_s(k) + (tf - 0.5) * b_t(k)];
            for (column, v) in columns.iter_mut().zip(values) {
                if let Some(c) = column.get_mut(k) {
                    c.set([s, t, 0], v);
                }
            }
        }
    }
    Ok(columns)
}

/// The Jacobian of a solid element whose node coordinates are already checked.
pub(crate) fn jacobian(kind: ElementKind, xyz: &[[f64; 3]]) -> Result<Jacobian> {
    let columns = if kind == ElementKind::Pyramid5 {
        pyramid_columns(xyz)?
    } else {
        let (_, dirs) = domain(kind)?;
        let g = geometry(kind, xyz)?;
        let column =
            |(f, k): (usize, usize)| -> Result<[Bernstein; 3]> { Ok([g[0].derivative(f, k)?, g[1].derivative(f, k)?, g[2].derivative(f, k)?]) };
        [column(dirs[0])?, column(dirs[1])?, column(dirs[2])?]
    };
    let [a, b, c] = &columns;
    // det [a b c] = a · (b × c)
    let cross = |i: usize, j: usize| -> Result<Bernstein> { b[i].mul(&c[j])?.axpy(-1.0, &b[j].mul(&c[i])?) };
    let det = a[0].mul(&cross(1, 2)?)?.axpy(1.0, &a[1].mul(&cross(2, 0)?)?)?.axpy(1.0, &a[2].mul(&cross(0, 1)?)?)?;
    Ok(Jacobian { columns, det })
}

/// A piece of the reference domain during certification: the determinant on it and its vertices
/// in reference coordinates, per factor.
struct Piece {
    det: Bernstein,
    /// The determinant's smallest coefficient on this piece.
    min: f64,
    verts: Vec<Vec<[f64; 3]>>,
    depth: u32,
}

impl Piece {
    fn new(det: Bernstein, verts: Vec<Vec<[f64; 3]>>, depth: u32) -> Self {
        let min = det.bounds().0;
        Piece { det, min, verts, depth }
    }

    /// The longest edge over all factors (first on ties, so the choice is deterministic).
    fn longest_edge(&self) -> (usize, usize, usize) {
        let mut best = (0, 0, 1);
        let mut len = -1.0;
        for (f, vs) in self.verts.iter().enumerate() {
            for (i, a) in vs.iter().enumerate() {
                for (j, b) in vs.iter().enumerate().skip(i + 1) {
                    let l = (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2);
                    if l > len {
                        len = l;
                        best = (f, i, j);
                    }
                }
            }
        }
        best
    }

    fn split(&self) -> Result<(Piece, Piece)> {
        let (f, i, j) = self.longest_edge();
        let (a, b) = self.det.bisect(f, i, j)?;
        let vs = self.verts.get(f).ok_or_else(|| Error::Failed("certificate: lost a factor".into()))?;
        let (pi, pj) = match (vs.get(i), vs.get(j)) {
            (Some(p), Some(q)) => (*p, *q),
            _ => return Err(Error::Failed("certificate: lost a vertex".into())),
        };
        let m = mid(pi, pj);
        let child = |det: Bernstein, replace: usize| {
            let mut verts = self.verts.clone();
            if let Some(v) = verts.get_mut(f).and_then(|vs| vs.get_mut(replace)) {
                *v = m;
            }
            Piece::new(det, verts, self.depth + 1)
        };
        Ok((child(a, j), child(b, i)))
    }
}

/// An upper bound of `|∂x/∂ξ_a|` over the whole element for one Jacobian column (a vector
/// polynomial in Bernstein form): the longest control vector, since the polynomial stays in the
/// convex hull of its control points.
fn column_norm_bound(column: &[Bernstein; 3]) -> f64 {
    let [x, y, z] = column.each_ref().map(Bernstein::coefficients);
    x.iter().zip(&y).zip(&z).map(|((a, b), c)| (a * a + b * b + c * c).sqrt()).fold(0.0, f64::max)
}

/// The determinant's floor for [`certify`] that proves FieldCraft's validity criterion
/// ([`fieldcraft_model::admissible_jacobian`]) on the whole element: with `N` the product of the
/// columns' norm bounds, `det J(ξ) > 2 · MIN_NORMALIZED_JACOBIAN · N` everywhere gives
/// `det J(ξ) / Π|J_a(ξ)| > MIN_NORMALIZED_JACOBIAN` at every point, so the solvers, which check
/// the same criterion at their integration points, accept every element the gate accepts. The
/// factor 2 covers the rounding of the solvers' own Jacobian evaluation.
pub(crate) fn admissible_floor(columns: &[[Bernstein; 3]; 3]) -> f64 {
    let n: f64 = columns.iter().map(column_norm_bound).product();
    2.0 * fieldcraft_model::MIN_NORMALIZED_JACOBIAN * n
}

/// Decide whether a Jacobian determinant exceeds `floor` (and [`RELATIVE_TOLERANCE`] of its
/// largest magnitude) on the whole reference domain.
pub(crate) fn certify(det: Bernstein, floor: f64) -> Result<JacobianCheck> {
    let (lo, hi) = det.bounds();
    let volume = det.integral();
    let scale = lo.abs().max(hi.abs());
    if !scale.is_finite() || !volume.is_finite() || !floor.is_finite() {
        return Err(Error::Invalid("the element's Jacobian is not a finite number".into()));
    }
    let tol = (RELATIVE_TOLERANCE * scale).max(floor);
    let ratio_of = |v: f64| if scale > 0.0 { v / scale } else { 0.0 };
    let verts = det.factors().iter().map(|f| (0..=f.dim).map(|k| simplex_vertex(k, 1).map(|x| x as f64)).collect()).collect();
    let mut stack = vec![Piece::new(det, verts, 0)];
    let mut lower = f64::INFINITY;
    let mut undetermined = false;
    let mut pieces = 0;
    while let Some(piece) = stack.pop() {
        pieces += 1;
        if pieces > MAX_PIECES {
            undetermined = true;
            break;
        }
        if piece.det.vertex_values().iter().any(|v| !(*v > tol)) {
            return Ok(JacobianCheck { validity: Validity::Invalid, ratio: ratio_of(lo), volume });
        }
        if piece.min > tol {
            lower = lower.min(piece.min);
            continue;
        }
        if piece.depth >= MAX_DEPTH {
            undetermined = true;
            continue;
        }
        let (a, b) = piece.split()?;
        // Look at the more doubtful piece first: a negative value is found sooner.
        let (first, second) = if a.min <= b.min { (a, b) } else { (b, a) };
        stack.push(second);
        stack.push(first);
    }
    if undetermined {
        return Ok(JacobianCheck { validity: Validity::Undetermined, ratio: ratio_of(lo), volume });
    }
    Ok(JacobianCheck { validity: Validity::Valid, ratio: ratio_of(lower), volume })
}
