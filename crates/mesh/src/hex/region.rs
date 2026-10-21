//! Planar face meshing for the sweepers: the scheme a region gets, the interval constraints the
//! scheme needs, and the quadrilateral mesh on the assigned boundary nodes.
//!
//! - **Mapped**: one loop with four corners (vertex classification by interior angle, as in
//!   submapping: Whiteley, White, Benzley, Blacker, *Two and three-quarter dimensional meshing
//!   facilitators*, Engineering with Computers 12:144–154, 1996): opposite sides equal,
//!   transfinite interpolation ([`crate::quad::mapped_quad_mesh`]).
//! - **Annulus**: two smooth loops (no corners), one inside the other: the two loops get equal
//!   counts and the ring between them is a periodic structured grid, each radial line joining a
//!   node of the outer loop to its partner on the inner loop (the rotation that pairs them most
//!   closely) — the O-grid of a tube's cross-section.
//! - **Triangle**: one loop with three corners (a fillet's quarter disk, the corner of a
//!   pocket cut off by a partition): the triangle primitive of midpoint subdivision (Tam and
//!   Armstrong, *2D finite element mesh generation by medial axis subdivision*, Advances in
//!   Engineering Software 13(5–6):313–324, 1991; Li, Armstrong, McKeag, *Quad mesh generation for
//!   k-sided faces and hex mesh generation for trivalent polyhedra*, Finite Elements in Analysis
//!   and Design 26(4):279–301, 1997): lines from a centre to one node on each side split the
//!   region into three mapped quadrilateral blocks. Side `k` with `n_k` divisions is split into
//!   `m_{k-1}` and `m_{k+1}` divisions, the counts of the lines that end on the sides before
//!   and after it, so `m_k = (n_{k+1} + n_{k+2} − n_k) / 2`: the loop's count is even and each
//!   side shorter than the other two together.
//! - **Disk**: one smooth loop (no corners; a straight side running tangent into an arc is no
//!   corner): an O-grid ("butterfly"): a mapped rectangle in the middle on the region's principal
//!   axes, inset from the boundary by the same distance all round (on a circle, the square whose
//!   corners are halfway to the boundary; on a slot, a long thin block), and rings of
//!   quadrilaterals from it out to the loop, then smoothed. The loop's count is a multiple of
//!   four.
//! - **Free**: anything else: even loops for the free all-quad mesher
//!   ([`crate::quad::free_quad_mesh`]).
//!
//! A structured scheme whose grid would fold falls back to the free mesher (its constraints imply
//! even loops).

use std::f64::consts::PI;

use super::interval::Intervals;
use crate::quad::{QuadMesh, QuadOptions, QuadStats, free_quad_mesh, mapped_quad_mesh, transfinite_grid};
use crate::quality::{Validity, planar_quad_validity};
use crate::size::SizeField;
use crate::{Error, Result};

type Q = [f64; 2];

/// The most nodes the annulus pairing compares (its cost is quadratic in the loop length).
const MAX_ANNULUS_NODES: usize = 4096;

/// Smoothing sweeps over the interior of a triangle primitive.
const TRIANGLE_SWEEPS: usize = 20;

/// How close (radians) an interior angle must come to a classification threshold to count as on
/// it. Exact geometry lands on the thresholds (an O-grid's 45° cuts meet straight sides at
/// exactly 135°), and the computed angle is then rounded either way by the last bits of its
/// points: those differ between platforms' `sin`/`cos` and between tessellations, and the
/// scheme must not (#556). An angle within this of a threshold is taken as exactly on it, and
/// on a threshold is smooth, neither a corner nor reflex: that keeps the slot's intended O-grid.
pub(crate) const ANGLE_TOLERANCE: f64 = 1e-9; // unit-factors: an angle in radians

/// A corner: an interior angle clearly below 135° (exactly 135° is not one).
pub(crate) fn is_corner(angle: f64) -> bool {
    angle < 0.75 * PI - ANGLE_TOLERANCE
}

/// A reflex turn no mapped or triangle region can take: clearly above 225°.
pub(crate) fn is_reflex(angle: f64) -> bool {
    angle > 1.25 * PI + ANGLE_TOLERANCE
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Scheme {
    /// Corner positions (curve uses) in the single loop.
    Mapped([usize; 4]),
    /// Corner positions (curve uses) in the single loop.
    Triangle([usize; 3]),
    Annulus,
    Disk,
    Free,
}

/// The interior angle at every curve end of a loop, given each curve use's oriented polyline
/// (with the loop counter-clockwise when `ccw`).
fn angles(lp: &[Vec<Q>], ccw: bool) -> Option<Vec<f64>> {
    let n = lp.len();
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let prev = lp.get((k + n - 1) % n)?;
        let next = lp.get(k)?;
        let (a, o, b) = (*prev.get(prev.len().checked_sub(2)?)?, *next.first()?, *next.get(1)?);
        let (u, v) = ([o[0] - a[0], o[1] - a[1]], [b[0] - o[0], b[1] - o[1]]);
        let turn = crate::vec3::atan2(u[0] * v[1] - u[1] * v[0], u[0] * v[0] + u[1] * v[1]);
        out.push(PI - if ccw { turn } else { -turn });
    }
    Some(out)
}

fn signed_area(pts: &[Q]) -> f64 {
    let n = pts.len();
    0.5 * (0..n).map(|k| pts[k][0] * pts[(k + 1) % n][1] - pts[(k + 1) % n][0] * pts[k][1]).sum::<f64>()
}

/// Choose the scheme of a region from its loops: per loop, the oriented polyline of each curve
/// use (loop 0 outer).
pub(crate) fn classify(loops: &[Vec<Vec<Q>>]) -> Scheme {
    let ccw = |lp: &Vec<Vec<Q>>| {
        let pts: Vec<Q> = lp.iter().flat_map(|c| c.iter().take(c.len().saturating_sub(1)).copied()).collect();
        signed_area(&pts) > 0.0
    };
    match loops {
        [lp] => {
            let Some(a) = angles(lp, ccw(lp)) else { return Scheme::Free };
            if a.iter().any(|&x| is_reflex(x)) {
                return Scheme::Free;
            }
            let corners: Vec<usize> = (0..a.len()).filter(|&k| a.get(k).is_some_and(|&x| is_corner(x))).collect();
            if corners.is_empty() {
                return Scheme::Disk;
            }
            match corners[..] {
                [a, b, c] => Scheme::Triangle([a, b, c]),
                [a, b, c, d] => Scheme::Mapped([a, b, c, d]),
                _ => Scheme::Free,
            }
        }
        [outer, inner] => {
            let smooth = |lp: &Vec<Vec<Q>>| angles(lp, ccw(lp)).is_some_and(|a| a.iter().all(|&x| !is_corner(x) && !is_reflex(x)));
            if smooth(outer) && smooth(inner) { Scheme::Annulus } else { Scheme::Free }
        }
        _ => Scheme::Free,
    }
}

/// The curve variables of each side of a loop between its corners (positions in increasing
/// order); the last side wraps round to the first corner.
fn loop_sides<T: Copy, const N: usize>(lp: &[T], c: [usize; N]) -> [Vec<T>; N] {
    std::array::from_fn(|k| {
        let (a, b) = (c[k], c.get(k + 1).copied().unwrap_or(c[0] + lp.len()));
        (a..b).filter_map(|i| lp.get(i % lp.len()).copied()).collect()
    })
}

/// The opposite sides of a mapped loop that are one curve each, as uses traversed the same way
/// (the loop runs along one side and back along the other), for carrying a biased seed's
/// spacing across ([`crate::seeds::carry`]).
pub(crate) fn opposite_ties(lp: &[(usize, bool)], scheme: Scheme) -> Vec<[crate::seeds::Use; 2]> {
    let Scheme::Mapped(c) = scheme else { return Vec::new() };
    let [s0, s1, s2, s3] = loop_sides(lp, c);
    [(s0, s2), (s1, s3)]
        .into_iter()
        .filter_map(|(a, b)| match (a.as_slice(), b.as_slice()) {
            ([a], [(b, fwd)]) => Some([*a, (*b, !fwd)]),
            _ => None,
        })
        .collect()
}

/// Add the scheme's constraints on the loops' curve variables.
pub(crate) fn constrain(ia: &mut Intervals, loops: &[Vec<usize>], scheme: Scheme) {
    match scheme {
        Scheme::Mapped(c) => {
            let Some(lp) = loops.first() else { return };
            let [s0, s1, s2, s3] = loop_sides(lp, c);
            ia.equal.push((s0, s2));
            ia.equal.push((s1, s3));
        }
        Scheme::Triangle(c) => {
            // Even, and each side at least two divisions shorter than the other two together
            // (each line from the centre at least one division).
            let Some(lp) = loops.first() else { return };
            let sides = loop_sides(lp, c);
            for k in 0..3 {
                let others = sides[(k + 1) % 3].iter().chain(&sides[(k + 2) % 3]).copied().collect();
                ia.exceed.push((others, sides[k].clone(), 2));
            }
            ia.multiple.push((lp.clone(), 2));
        }
        Scheme::Annulus => {
            if let (Some(o), Some(i)) = (loops.first(), loops.get(1)) {
                ia.equal.push((o.clone(), i.clone()));
                ia.multiple.push((o.clone(), 2));
                ia.at_least.push((o.clone(), 6));
            }
        }
        Scheme::Disk => {
            if let Some(lp) = loops.first() {
                ia.multiple.push((lp.clone(), 4));
                ia.at_least.push((lp.clone(), 8));
            }
        }
        Scheme::Free => {
            let min = if loops.len() > 1 { 6 } else { 4 };
            for lp in loops {
                ia.multiple.push((lp.clone(), 2));
                ia.at_least.push((lp.clone(), min));
            }
        }
    }
}

/// Mesh a region on its boundary nodes: `loops[k]` are loop `k`'s nodes in order, `starts[k]` the
/// positions where its curve uses start. Returns the mesh (boundary nodes first, in the order of
/// the loops) and the scheme used.
pub(crate) fn mesh(loops: &[Vec<Q>], starts: &[Vec<usize>], scheme: Scheme, size: &dyn SizeField, opts: &QuadOptions) -> Result<(QuadMesh, Scheme)> {
    match mesh_structured(loops, starts, scheme, size, opts)? {
        Some(m) => Ok((m, scheme)),
        None => Ok((free_quad_mesh(loops, size, opts)?, Scheme::Free)),
    }
}

/// The structured mesh of a region with its scheme, or `None` when the scheme is free or its
/// grid would fold.
pub(crate) fn mesh_structured(
    loops: &[Vec<Q>],
    starts: &[Vec<usize>],
    scheme: Scheme,
    size: &dyn SizeField,
    opts: &QuadOptions,
) -> Result<Option<QuadMesh>> {
    Ok(match scheme {
        Scheme::Mapped(c) => {
            let s0 = starts.first();
            let pos = c.map(|k| s0.and_then(|s| s.get(k)).copied().unwrap_or(0));
            match loops.first().map(|b| mapped_quad_mesh(b, pos, opts)) {
                Some(Ok(m)) => Some(m.mesh),
                Some(Err(e)) if e.kind() == Some(crate::Kind::Limit) => return Err(e),
                _ => None,
            }
        }
        Scheme::Triangle(c) => {
            let s0 = starts.first();
            let pos = c.map(|k| s0.and_then(|s| s.get(k)).copied().unwrap_or(0));
            loops.first().map(|b| triangle(b, pos, opts)).transpose()?.flatten()
        }
        Scheme::Annulus => annulus(loops, size, opts)?,
        Scheme::Disk => disk(loops, size, opts)?,
        Scheme::Free => None,
    })
}

/// The triangle primitive on one loop with corners at node positions `corners` (in loop order),
/// or `None` when the side counts allow no split or a block would fold.
fn triangle(boundary: &[Q], corners: [usize; 3], opts: &QuadOptions) -> Result<Option<QuadMesh>> {
    let n = boundary.len();
    if !(6..=MAX_ANNULUS_NODES).contains(&n) || corners.iter().any(|&c| c >= n) {
        return Ok(None);
    }
    // Corners in loop order from the first; side k runs from corner k to corner k + 1.
    let mut c = corners;
    c.sort_unstable();
    let len = |k: usize| -> usize {
        let (a, b) = (c[k % 3], c[(k + 1) % 3]);
        (b + n - a) % n
    };
    let sides = [len(0), len(1), len(2)];
    let total: usize = sides.iter().sum();
    if total != n || !total.is_multiple_of(2) {
        return Ok(None);
    }
    // m[k]: divisions of the line from the centre to side k's split node.
    let mut m = [0usize; 3];
    for k in 0..3 {
        let other = sides[(k + 1) % 3] + sides[(k + 2) % 3];
        if other <= sides[k] {
            return Ok(None);
        }
        m[k] = (other - sides[k]) / 2;
    }
    let count = (0..3).map(|k| m[k] * m[(k + 1) % 3]).sum::<usize>();
    if count > opts.max_elements {
        return Err(Error::Limit(format!("the mesh would have {count} quadrilaterals; the limit is {}", opts.max_elements)));
    }
    let at = |i: usize| boundary.get(i % n).copied().unwrap_or([0.0; 2]);
    // Side k is split after m[k + 2] divisions (the line count of the block at corner k).
    let split: [usize; 3] = std::array::from_fn(|k| (c[k] + m[(k + 2) % 3]) % n);
    let sp: [Q; 3] = split.map(at);
    let centre = [(sp[0][0] + sp[1][0] + sp[2][0]) / 3.0, (sp[0][1] + sp[1][1] + sp[2][1]) / 3.0];
    let mut nodes: Vec<Q> = boundary.to_vec();
    nodes.push(centre);
    let centre_id = n as u32;
    // Line k from the split node of side k to the centre, its interior nodes new.
    let mut lines: Vec<Vec<u32>> = Vec::with_capacity(3);
    for k in 0..3 {
        let mut ids = vec![split[k] as u32];
        for j in 1..m[k] {
            let t = j as f64 / m[k] as f64;
            nodes.push([sp[k][0] + t * (centre[0] - sp[k][0]), sp[k][1] + t * (centre[1] - sp[k][1])]);
            ids.push(nodes.len() as u32 - 1);
        }
        ids.push(centre_id);
        lines.push(ids);
    }
    let ccw = signed_area(boundary) > 0.0;
    let mut quads = Vec::with_capacity(count);
    // Block k at corner k: bottom corner k → split k (side k), right line k, top line k − 1
    // (split k − 1 → centre), left corner k → split k − 1 (side k − 1, backwards).
    for k in 0..3 {
        let prev = (k + 2) % 3;
        let (nu, nv) = (m[prev], m[k]);
        let bottom: Vec<u32> = (0..=nu).map(|j| ((c[k] + j) % n) as u32).collect();
        let left: Vec<u32> = (0..=nv).map(|j| ((c[k] + n - j) % n) as u32).collect();
        let (Some(right), Some(top)) = (lines.get(k), lines.get(prev)) else { return Ok(None) };
        if !fill_block(&mut nodes, [&bottom, right, top, &left], !ccw, &mut quads)? {
            return Ok(None);
        }
    }
    // The centre and the lines from it are placed by construction only: smooth them, as the
    // midpoint-subdivision papers do, keeping every quadrilateral valid.
    crate::quad::smooth::smooth(&mut nodes, &quads, n, TRIANGLE_SWEEPS)?;
    Ok(Some(QuadMesh { nodes, quads, boundary_nodes: n, stats: QuadStats::default() }))
}

/// Fill the block bounded by node sequences `bottom` and `top` (left to right) and `left` and
/// `right` (bottom to top) with a transfinite grid, pushing its interior nodes and its
/// quadrilaterals (reversed when `flip`). `false` when a quadrilateral would not be valid.
fn fill_block(nodes: &mut Vec<Q>, [bottom, right, top, left]: [&[u32]; 4], flip: bool, quads: &mut Vec<[u32; 4]>) -> Result<bool> {
    let (nu, nv) = (bottom.len().saturating_sub(1), left.len().saturating_sub(1));
    let at = |nodes: &[Q], i: u32| nodes.get(i as usize).copied().unwrap_or([0.0; 2]);
    let pts = |ids: &[u32], nodes: &[Q]| ids.iter().map(|&i| at(nodes, i)).collect::<Vec<Q>>();
    let grid = transfinite_grid(&pts(bottom, nodes), &pts(right, nodes), &pts(top, nodes), &pts(left, nodes))?;
    let mut id = vec![0u32; (nu + 1) * (nv + 1)];
    for j in 0..=nv {
        for i in 0..=nu {
            let v = if j == 0 {
                bottom.get(i).copied()
            } else if j == nv {
                top.get(i).copied()
            } else if i == 0 {
                left.get(j).copied()
            } else if i == nu {
                right.get(j).copied()
            } else {
                nodes.push(grid.get(j * (nu + 1) + i).copied().unwrap_or([0.0; 2]));
                Some(nodes.len() as u32 - 1)
            };
            if let (Some(slot), Some(v)) = (id.get_mut(j * (nu + 1) + i), v) {
                *slot = v;
            }
        }
    }
    let g = |i: usize, j: usize| id.get(j * (nu + 1) + i).copied().unwrap_or(0);
    for j in 0..nv {
        for i in 0..nu {
            let q = [g(i, j), g(i + 1, j), g(i + 1, j + 1), g(i, j + 1)];
            let q = if flip { [q[0], q[3], q[2], q[1]] } else { q };
            if planar_quad_validity(&q.map(|v| at(nodes, v)))? != Validity::Valid {
                return Ok(false);
            }
            quads.push(q);
        }
    }
    Ok(true)
}

/// The O-grid between two loops of equal length, or `None` when it would fold.
fn annulus(loops: &[Vec<Q>], size: &dyn SizeField, opts: &QuadOptions) -> Result<Option<QuadMesh>> {
    let [outer, inner] = loops else { return Ok(None) };
    let n = outer.len();
    if n != inner.len() || !(3..=MAX_ANNULUS_NODES).contains(&n) {
        return Ok(None);
    }
    // Index maps that make both loops counter-clockwise.
    let ccw_order = |pts: &[Q]| -> Vec<usize> { if signed_area(pts) > 0.0 { (0..pts.len()).collect() } else { (0..pts.len()).rev().collect() } };
    let (oi, ii) = (ccw_order(outer), ccw_order(inner));
    let o = |k: usize| oi.get(k % n).and_then(|&i| outer.get(i)).copied().unwrap_or([0.0; 2]);
    let inn = |k: usize| ii.get(k % n).and_then(|&i| inner.get(i)).copied().unwrap_or([0.0; 2]);
    let d2 = |a: Q, b: Q| (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2);
    let shift = (0..n)
        .min_by(|&s, &t| {
            let cost = |s: usize| (0..n).map(|k| d2(o(k), inn(k + s))).sum::<f64>();
            cost(s).total_cmp(&cost(t)).then(s.cmp(&t))
        })
        .unwrap_or(0);
    let mean_gap = (0..n).map(|k| d2(o(k), inn(k + shift)).sqrt()).sum::<f64>() / n as f64;
    let mid = |k: usize| {
        let (a, b) = (o(k), inn(k + shift));
        [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5, 0.0]
    };
    let h = crate::size::checked_size(size, mid(0))?;
    let m = (mean_gap / h).round().max(1.0) as usize;
    let quads_count = n.saturating_mul(m);
    if quads_count > opts.max_elements {
        return Err(Error::Limit(format!("the mesh would have {quads_count} quadrilaterals; the limit is {}", opts.max_elements)));
    }
    // Node ids: outer loop 0..n, inner loop n..2n (input order), then the rings between.
    let outer_ids: Vec<u32> = (0..n).map(|k| oi.get(k % n).copied().unwrap_or(0) as u32).collect();
    let inner_ids: Vec<u32> = (0..n).map(|k| (n + ii.get((k + shift) % n).copied().unwrap_or(0)) as u32).collect();
    let mut nodes: Vec<Q> = outer.iter().chain(inner.iter()).copied().collect();
    let mut quads = Vec::with_capacity(quads_count);
    if !rings(&mut nodes, &outer_ids, &inner_ids, m, true, &mut quads)? {
        return Ok(None);
    }
    Ok(Some(QuadMesh { nodes, quads, boundary_nodes: 2 * n, stats: QuadStats::default() }))
}

/// Rings of quadrilaterals between a closed outer node sequence (counter-clockwise) and an inner
/// one of the same length, node `k` of each joined by a straight radial line divided into `m`.
/// `false` when a quadrilateral would not be strictly convex, unless `check` is off (a caller
/// that smooths and checks the mesh itself).
fn rings(nodes: &mut Vec<Q>, outer: &[u32], inner: &[u32], m: usize, check: bool, quads: &mut Vec<[u32; 4]>) -> Result<bool> {
    let n = outer.len();
    let at = |nodes: &[Q], i: u32| nodes.get(i as usize).copied().unwrap_or([0.0; 2]);
    let mut ring_ids: Vec<Vec<u32>> = vec![outer.to_vec()];
    for r in 1..m {
        let t = r as f64 / m as f64;
        let mut ids = Vec::with_capacity(n);
        for k in 0..n {
            let (a, b) = (at(nodes, outer.get(k).copied().unwrap_or(0)), at(nodes, inner.get(k).copied().unwrap_or(0)));
            nodes.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
            ids.push(nodes.len() as u32 - 1);
        }
        ring_ids.push(ids);
    }
    ring_ids.push(inner.to_vec());
    for w in ring_ids.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        for k in 0..n {
            let (Some(&p0), Some(&p1), Some(&q1), Some(&q0)) = (a.get(k), a.get((k + 1) % n), b.get((k + 1) % n), b.get(k)) else {
                return Ok(false);
            };
            // Outer ring counter-clockwise, inner ring inside it: (p0, p1, q1, q0) turns left.
            let q = [p0, p1, q1, q0];
            if check && planar_quad_validity(&q.map(|i| at(nodes, i)))? != Validity::Valid {
                return Ok(false);
            }
            quads.push(q);
        }
    }
    Ok(true)
}

/// The O-grid of one smooth loop whose length is a multiple of four (at least 8), or `None` when
/// it would fold.
fn disk(loops: &[Vec<Q>], size: &dyn SizeField, opts: &QuadOptions) -> Result<Option<QuadMesh>> {
    let [outer] = loops else { return Ok(None) };
    let n = outer.len();
    if n < 8 || !n.is_multiple_of(4) || n > MAX_ANNULUS_NODES {
        return Ok(None);
    }
    let ccw = signed_area(outer) > 0.0;
    let at = |nodes: &[Q], i: u32| nodes.get(i as usize).copied().unwrap_or([0.0; 2]);
    // The area centroid and second moments about it (Green's theorem over the polygon).
    let (mut a2, mut cx, mut cy) = (0.0, 0.0, 0.0);
    for k in 0..n {
        let (p, r) = (outer[k], outer[(k + 1) % n]);
        let c = p[0] * r[1] - r[0] * p[1];
        a2 += c;
        cx += (p[0] + r[0]) * c;
        cy += (p[1] + r[1]) * c;
    }
    if !(a2.abs() > 0.0) {
        return Ok(None);
    }
    let centre = [cx / (3.0 * a2), cy / (3.0 * a2)];
    let (mut ixx, mut iyy, mut ixy) = (0.0, 0.0, 0.0);
    for k in 0..n {
        let (p, r) = (outer[k], outer[(k + 1) % n]);
        let (p, r) = ([p[0] - centre[0], p[1] - centre[1]], [r[0] - centre[0], r[1] - centre[1]]);
        let c = p[0] * r[1] - r[0] * p[1];
        ixx += c * (p[1] * p[1] + p[1] * r[1] + r[1] * r[1]);
        iyy += c * (p[0] * p[0] + p[0] * r[0] + r[0] * r[0]);
        ixy += c * (p[0] * r[1] + 2.0 * p[0] * p[1] + 2.0 * r[0] * r[1] + r[0] * p[1]);
    }
    // ∫y², ∫x² (both × 12) and ∫xy (× 24), signed as the loop: the region's long direction is
    // the principal axis of the largest ∫ along², the angle of 0.5 · atan2(2 ∫xy, ∫x² − ∫y²).
    let sign = a2.signum();
    let theta = 0.5 * crate::vec3::atan2(sign * ixy, sign * (iyy - ixx));
    let (u, v) = {
        let (s, c) = (libm::sin(theta), libm::cos(theta));
        ([c, s], [-s, c])
    };
    let local = |p: Q| {
        let d = [p[0] - centre[0], p[1] - centre[1]];
        [d[0] * u[0] + d[1] * u[1], d[0] * v[0] + d[1] * v[1]]
    };
    // Half-widths of the region along its axes, from the centroid to the farther boundary.
    let (mut a, mut b) = (0.0f64, 0.0f64);
    for &p in outer.iter() {
        let l = local(p);
        a = a.max(l[0].abs());
        b = b.max(l[1].abs());
    }
    if !(a > 0.0 && b > 0.0) {
        return Ok(None);
    }
    // The middle block: a rectangle on the axes inset from the boundary by the same distance all
    // round, its half-widths (a − b + c, c) with c = b / (2√2) (on a circle, the square whose
    // corners are halfway to the boundary). A square middle block on an elongated region (a slot,
    // a rounded pocket) squashes the elements at its long sides (#198).
    let half_short = b / (2.0 * std::f64::consts::SQRT_2);
    let half_long = (a - b + half_short).max(half_short);
    let signs: [Q; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
    let world = |l: Q| [centre[0] + l[0] * u[0] + l[1] * v[0], centre[1] + l[0] * u[1] + l[1] * v[1]];
    let corners: [Q; 4] = signs.map(|s| world([s[0] * half_long, s[1] * half_short]));
    // Each corner's partner on the loop is a node near the point the inset carries the corner
    // out to, along its diagonal; the loop nodes between two partners pair with the side between
    // the corners, so the sides' divisions follow: q1 along the long axis, q2 across, their sum
    // half the loop's.
    let inset = b - half_short;
    let ccw_order: Vec<u32> = if ccw { (0..n as u32).collect() } else { (0..n as u32).rev().collect() };
    let targets: [Q; 4] = std::array::from_fn(|k| world([signs[k][0] * (half_long + inset), signs[k][1] * (half_short + inset)]));
    let dist = |i: usize, t: Q| {
        let p = at(outer, ccw_order.get(i % n).copied().unwrap_or(0));
        (p[0] - t[0]).powi(2) + (p[1] - t[1]).powi(2)
    };
    let nearest = |t: Q| (0..n).min_by(|&x, &y| dist(x, t).total_cmp(&dist(y, t))).unwrap_or(0);
    let p: [usize; 4] = targets.map(nearest);
    let run = |k: usize| (p[(k + 1) % 4] + n - p[k]) % n;
    let half = n / 2;
    // The four partners sit at `start`, `start + q1`, `start + half` and `start + half + q1`:
    // the start and q1 (near the mean of the two long runs) that put them nearest their targets.
    let mean = (run(0) + run(2)) / 2;
    let mut best = (f64::INFINITY, p[0], mean.clamp(1, half - 1));
    for q1 in mean.saturating_sub(2).max(1)..=(mean + 2).min(half - 1) {
        for shift in 0..5 {
            let start = (p[0] + n + shift - 2) % n;
            let offsets = [0, q1, half, half + q1];
            let cost: f64 = (0..4).map(|k| dist(start + offsets[k], targets[k])).sum();
            if cost < best.0 {
                best = (cost, start, q1);
            }
        }
    }
    let (_, start, q1) = best;
    let q2 = half - q1;
    let offsets = [0, q1, q1 + q2, 2 * q1 + q2];
    let lengths = [q1, q2, q1, q2];
    let order: Vec<u32> = (0..n).map(|k| ccw_order.get((start + k) % n).copied().unwrap_or(0)).collect();
    let mut nodes: Vec<Q> = outer.clone();
    let mut square: Vec<u32> = Vec::with_capacity(n);
    // Each side's nodes spaced as the loop nodes they pair with (by chord length along the
    // run), so the lines between the two run straight across and do not shear.
    for k in 0..4 {
        let (c0, c1) = (corners[k], corners[(k + 1) % 4]);
        let run: Vec<Q> = (0..=lengths[k]).map(|j| at(outer, order.get((offsets[k] + j) % n).copied().unwrap_or(0))).collect();
        let along: Vec<f64> = run.windows(2).map(|w| ((w[1][0] - w[0][0]).powi(2) + (w[1][1] - w[0][1]).powi(2)).sqrt()).collect();
        let total: f64 = along.iter().sum();
        let mut acc = 0.0;
        for j in 0..lengths[k] {
            let t = if total > 0.0 { acc / total } else { j as f64 / lengths[k] as f64 };
            acc += along.get(j).copied().unwrap_or(0.0);
            nodes.push([c0[0] + t * (c1[0] - c0[0]), c0[1] + t * (c1[1] - c0[1])]);
            square.push(nodes.len() as u32 - 1);
        }
    }
    let gap = (0..n).map(|k| {
        let (a, b) = (at(&nodes, order.get(k).copied().unwrap_or(0)), at(&nodes, square.get(k).copied().unwrap_or(0)));
        ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt()
    });
    let mean_gap = gap.sum::<f64>() / n as f64;
    let h = crate::size::checked_size(size, [centre[0], centre[1], 0.0])?;
    let m = (mean_gap / h).round().max(1.0) as usize;
    let count = n.saturating_mul(m).saturating_add(q1 * q2);
    if count > opts.max_elements {
        return Err(Error::Limit(format!("the mesh would have {count} quadrilaterals; the limit is {}", opts.max_elements)));
    }
    let mut quads = Vec::with_capacity(count);
    if !rings(&mut nodes, &order, &square, m, false, &mut quads)? {
        return Ok(None);
    }
    // The middle block: a q1 × q2 transfinite grid on its perimeter (counter-clockwise from
    // corner 0).
    let side = |k: usize, rev: bool| -> Vec<u32> {
        let mut v: Vec<u32> = (0..=lengths[k]).map(|j| square.get((offsets[k] + j) % n).copied().unwrap_or(0)).collect();
        if rev {
            v.reverse();
        }
        v
    };
    let (bottom, right, top, left) = (side(0, false), side(1, false), side(2, true), side(3, true));
    if !fill_block(&mut nodes, [&bottom, &right, &top, &left], false, &mut quads)? {
        return Ok(None);
    }
    // The rings and the middle block are placed by construction only; where the loop's nodes
    // crowd round the block's corners (a slot's round ends) the straight radial lines shear the
    // rings, even fold them on a coarse loop: smooth them (the smoother untangles first), then
    // every quadrilateral must be valid.
    crate::quad::smooth::smooth(&mut nodes, &quads, n, TRIANGLE_SWEEPS)?;
    for q in &quads {
        if planar_quad_validity(&q.map(|i| at(&nodes, i)))? != Validity::Valid {
            return Ok(None);
        }
    }
    Ok(Some(QuadMesh { nodes, quads, boundary_nodes: n, stats: QuadStats::default() }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::size::UniformSize;

    /// An angle of exactly 135°, the threshold, is not a corner however the last bits of its
    /// points round (Windows' `sin` against glibc's, or two tessellations): a trapezoid with one
    /// stays a triangle rather than flipping to mapped when that vertex moves by an ulp (#556).
    #[test]
    fn an_angle_on_the_threshold_is_not_a_corner_whichever_way_it_rounds() {
        // (0, 0) 90°, (3, 0) 45°, (2, 1) 135°, (0, 1) 90°: three corners and a 135° bend.
        let trapezoid = |p: Q| {
            let pts = [[0.0, 0.0], [3.0, 0.0], p, [0.0, 1.0]];
            let curves: Vec<Vec<Q>> = (0..4).map(|k| vec![pts[k], pts[(k + 1) % 4]]).collect();
            classify(&[curves])
        };
        let steps: [fn(f64) -> f64; 3] = [|v| v, f64::next_up, f64::next_down];
        for nx in steps {
            for ny in steps {
                let mut p = [nx(2.0), ny(1.0)];
                // A few ulps further each way as well.
                for _ in 0..3 {
                    assert_eq!(trapezoid(p), Scheme::Triangle([0, 1, 3]), "{p:?}");
                    p = [nx(p[0]), ny(p[1])];
                }
            }
        }
    }

    fn circle(r: f64, n: usize, ccw: bool) -> Vec<Q> {
        (0..n)
            .map(|k| {
                let t = 2.0 * PI * k as f64 / n as f64 * if ccw { 1.0 } else { -1.0 };
                [r * t.cos(), r * t.sin()]
            })
            .collect()
    }

    #[test]
    fn rings_get_an_o_grid() {
        let (outer, inner) = (circle(5.0, 16, true), circle(4.0, 16, false));
        let polyline = |pts: &Vec<Q>| {
            let mut p = pts.clone();
            p.push(pts[0]);
            vec![p]
        };
        assert_eq!(classify(&[polyline(&outer), polyline(&inner)]), Scheme::Annulus);
        let (m, s) = mesh(&[outer.clone(), inner.clone()], &[vec![0], vec![0]], Scheme::Annulus, &UniformSize(0.5), &QuadOptions::default()).unwrap();
        assert_eq!(s, Scheme::Annulus);
        assert_eq!(m.quads.len(), 16 * 2);
        assert_eq!(m.boundary_nodes, 32);
        assert_eq!(&m.nodes[..16], &outer[..]);
        // Unequal loops fall back to the free mesher.
        let (_, s) =
            mesh(&[circle(5.0, 16, true), circle(2.0, 8, true)], &[vec![0], vec![0]], Scheme::Annulus, &UniformSize(1.0), &QuadOptions::default())
                .unwrap();
        assert_eq!(s, Scheme::Free);
    }

    #[test]
    fn disks_get_an_o_grid() {
        let c = circle(5.0, 16, false);
        let mut closed = c.clone();
        closed.push(c[0]);
        assert_eq!(classify(&[vec![closed]]), Scheme::Disk);
        let (m, s) = mesh(std::slice::from_ref(&c), &[vec![0]], Scheme::Disk, &UniformSize(1.5), &QuadOptions::default()).unwrap();
        assert_eq!(s, Scheme::Disk);
        assert_eq!(&m.nodes[..16], &c[..]);
        // 16 · 2 ring quadrilaterals and the 4 × 4 square.
        assert_eq!(m.quads.len(), 16 * 2 + 16);
        let area: f64 = m.quads.iter().map(|q| signed_area(&q.map(|i| m.nodes[i as usize]))).sum();
        assert!((area - signed_area(&circle(5.0, 16, true))).abs() < 1e-9, "{area}");
    }

    #[test]
    fn squares_are_mapped() {
        let sq: Vec<Vec<Q>> =
            vec![vec![[0.0, 0.0], [1.0, 0.0]], vec![[1.0, 0.0], [1.0, 1.0]], vec![[1.0, 1.0], [0.0, 1.0]], vec![[0.0, 1.0], [0.0, 0.0]]];
        assert_eq!(classify(std::slice::from_ref(&sq)), Scheme::Mapped([0, 1, 2, 3]));
        let l: Vec<Vec<Q>> = vec![
            vec![[0.0, 0.0], [2.0, 0.0]],
            vec![[2.0, 0.0], [2.0, 1.0]],
            vec![[2.0, 1.0], [1.0, 1.0]],
            vec![[1.0, 1.0], [1.0, 2.0]],
            vec![[1.0, 2.0], [0.0, 2.0]],
            vec![[0.0, 2.0], [0.0, 0.0]],
        ];
        assert_eq!(classify(&[l]), Scheme::Free);
    }

    /// A fillet's quarter disk (two straight sides and an arc, three corners) gets the triangle
    /// primitive: three mapped blocks around a centre, every quadrilateral convex, the boundary
    /// nodes kept and the area filled.
    #[test]
    fn quarter_disks_get_the_triangle_primitive() {
        // Corners at the origin, (5, 0) and (0, 5); sides of 4, 6 and 4 divisions.
        let mut pts: Vec<Q> = (0..4).map(|k| [5.0 * k as f64 / 4.0, 0.0]).collect();
        pts.extend((0..6).map(|k| {
            let t = 0.5 * PI * k as f64 / 6.0;
            [5.0 * t.cos(), 5.0 * t.sin()]
        }));
        pts.extend((0..4).map(|k| [0.0, 5.0 - 5.0 * k as f64 / 4.0]));
        let polyline = |a: usize, b: usize| (a..=b).map(|k| pts[k % pts.len()]).collect::<Vec<Q>>();
        let curves = vec![polyline(0, 4), polyline(4, 10), polyline(10, 14)];
        assert_eq!(classify(&[curves]), Scheme::Triangle([0, 1, 2]));
        let (m, s) =
            mesh(std::slice::from_ref(&pts), &[vec![0, 4, 10]], Scheme::Triangle([0, 1, 2]), &UniformSize(1.0), &QuadOptions::default()).unwrap();
        assert_eq!(s, Scheme::Triangle([0, 1, 2]));
        assert_eq!(&m.nodes[..14], &pts[..]);
        // Lines of (6 + 4 − 4) / 2 = 3, (4 + 4 − 6) / 2 = 1 and 3 divisions: 3·1 + 1·3 + 3·3 quads.
        assert_eq!(m.quads.len(), 15);
        for q in &m.quads {
            assert_eq!(planar_quad_validity(&q.map(|i| m.nodes[i as usize])).unwrap(), Validity::Valid);
        }
        let area: f64 = m.quads.iter().map(|q| signed_area(&q.map(|i| m.nodes[i as usize]))).sum();
        assert!((area - signed_area(&pts)).abs() < 1e-9, "{area}");
        // A side as long as the other two together has no split: the free mesher takes it.
        let (_, s) =
            mesh(std::slice::from_ref(&pts), &[vec![0, 1, 2]], Scheme::Triangle([0, 1, 2]), &UniformSize(1.0), &QuadOptions::default()).unwrap();
        assert_eq!(s, Scheme::Free);
    }

    /// A rounded rectangle `w` × `h` with corner radius `r` (counter-clockwise from the bottom
    /// side's start), each side and each quarter arc one curve use divided about `size` apart (a
    /// side of no length left out): every junction is a straight side tangent to an arc (180°),
    /// so the loop has no corners.
    fn rounded_rectangle(w: f64, h: f64, r: f64, size: f64) -> (Vec<Q>, Vec<Vec<Q>>, Vec<usize>) {
        let (a, b) = (0.5 * w - r, 0.5 * h - r);
        let count = |len: f64| ((len / size).round() as usize).max(1);
        let mut curves: Vec<Vec<Q>> = Vec::new();
        for k in 0..4 {
            let (half, other) = if k % 2 == 0 { (a, b) } else { (b, a) };
            // Side k runs along direction d at distance other + r from the centre, then the arc
            // round the corner after it.
            let t = 0.5 * PI * k as f64;
            let (d, nrm) = ([t.cos(), t.sin()], [t.sin(), -t.cos()]);
            let at = |s: f64, off: f64| [d[0] * s + nrm[0] * off, d[1] * s + nrm[1] * off];
            let n = count(2.0 * half);
            if half > 0.0 {
                curves.push((0..=n).map(|j| at(-half + 2.0 * half * j as f64 / n as f64, other + r)).collect());
            }
            let centre = at(half, other);
            let m = count(0.5 * PI * r);
            curves.push(
                (0..=m)
                    .map(|j| {
                        let u = t - 0.5 * PI + 0.5 * PI * j as f64 / m as f64;
                        [centre[0] + r * u.cos(), centre[1] + r * u.sin()]
                    })
                    .collect(),
            );
        }
        let mut pts = Vec::new();
        let mut starts = Vec::new();
        for c in &curves {
            starts.push(pts.len());
            pts.extend(c.iter().take(c.len() - 1).copied());
        }
        (pts, curves, starts)
    }

    /// The smallest corner scaled Jacobian of a planar quad mesh.
    fn min_scaled_jacobian(m: &QuadMesh) -> f64 {
        m.quads
            .iter()
            .flat_map(|q| {
                let p = q.map(|i| m.nodes[i as usize]);
                (0..4).map(move |i| crate::quality::planar_corner(&p, i).1.unwrap_or(-1.0))
            })
            .fold(f64::INFINITY, f64::min)
    }

    /// A slot-like region whose straight sides run tangent into its corner arcs (#198): its
    /// tangent junctions are no corners, so it is a disk, and its O-grid's middle block follows
    /// its elongation, so no element is squashed (a square middle block on a 6:1 region gave a
    /// scaled Jacobian of 0.16).
    #[test]
    fn elongated_disks_keep_their_elements_square() {
        for (w, h, r, size) in [
            (82.0, 14.0, 3.0, 1.0),
            (82.0, 14.0, 5.0, 3.0),
            (82.0, 14.0, 5.0, 2.0),
            (82.0, 14.0, 7.0, 1.5),
            (40.0, 30.0, 5.0, 2.0),
            (10.0, 10.0, 5.0, 0.5),
        ] {
            let (mut pts, curves, starts) = rounded_rectangle(w, h, r, size);
            assert_eq!(classify(std::slice::from_ref(&curves)), Scheme::Disk);
            // Pad to a multiple of four nodes, as the interval assignment would.
            while pts.len() % 4 != 0 {
                let (p, q) = (pts[0], pts[1]);
                pts.insert(1, [0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1])]);
            }
            let (m, s) = mesh(std::slice::from_ref(&pts), &[starts], Scheme::Disk, &UniformSize(size), &QuadOptions::default()).unwrap();
            assert_eq!(s, Scheme::Disk, "{w} × {h}");
            assert_eq!(&m.nodes[..pts.len()], &pts[..]);
            let area: f64 = m.quads.iter().map(|q| signed_area(&q.map(|i| m.nodes[i as usize]))).sum();
            assert!((area - signed_area(&pts)).abs() < 1e-9 * w * h, "{w} × {h}: {area}");
            let min = min_scaled_jacobian(&m);
            assert!(min > 0.5, "{w} × {h} r {r}: min scaled Jacobian {min}");
        }
    }
}
