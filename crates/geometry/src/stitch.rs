//! Stitching: closing the gaps of a triangulated surface under a tolerance, for meshing.
//!
//! An imported surface (a 3MF from another tool, faces written apart) can be open along seams
//! where two sides were tessellated separately: their boundary points lie a little apart, or one
//! side has points where the other has a straight run (a T-junction). [`stitch`] closes such
//! gaps when they are under a tolerance:
//! 1. **Vertices.** Points on open edges closer than the tolerance become one point (the first in
//!    index order), and triangles that collapse go.
//! 2. **T-junctions.** A point on an open edge that lies within the tolerance of the inside of
//!    another open edge splits that edge's triangle in two through it.
//!
//! Only open edges are touched, so a surface without gaps is unchanged; tags stay with the
//! triangles (a split triangle's two halves keep its tag). The result is deterministic: no
//! hash-order iteration (a grid in a `BTreeMap`, points in index order).

use std::collections::BTreeMap;

use crate::{Error, Result, SurfaceMesh};

/// A point is inside an edge when its parameter is this far from the ends.
const INSIDE: f64 = 1e-9; // unit-factors: a fraction of the edge

/// Most vertices a surface stitched may have.
pub const MAX_VERTICES: usize = 20_000_000;
/// Most grid cells the T-junction search visits in one pass (a bound, far above what the grid
/// sized to the open edges needs).
const MAX_CELL_VISITS: usize = 200_000_000;
/// Most passes splitting T-junctions (each split can open another).
const MAX_PASSES: usize = 16;

/// What [`stitch`] did.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct StitchReport {
    /// Open edges before and after.
    pub open_before: usize,
    pub open_after: usize,
    /// Points made one with another.
    pub merged: usize,
    /// Triangles split at a T-junction.
    pub split: usize,
    /// Triangles that collapsed and went.
    pub removed: usize,
}

fn dist2(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|k| (a[k] - b[k]).powi(2)).sum()
}

/// For each point of `points` (indices into `positions`, in increasing order), the first of them
/// within `d` (in that order): a grid of cell `d`; identical coordinates only for `d` = 0.
fn cluster(positions: &[[f64; 3]], points: &[usize], d: f64) -> BTreeMap<usize, usize> {
    let mut rep = BTreeMap::new();
    if !(d > 0.0) {
        let mut first: BTreeMap<[u64; 3], usize> = BTreeMap::new();
        for &i in points {
            let Some(p) = positions.get(i) else { continue };
            rep.insert(i, *first.entry(p.map(f64::to_bits)).or_insert(i));
        }
        return rep;
    }
    let cell = |p: [f64; 3]| p.map(|c| (c / d).floor() as i64);
    let mut grid: BTreeMap<[i64; 3], Vec<usize>> = BTreeMap::new();
    for &i in points {
        let Some(&p) = positions.get(i) else { continue };
        let c = cell(p);
        let mut found: Option<usize> = None;
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    for &j in grid.get(&[c[0] + dx, c[1] + dy, c[2] + dz]).into_iter().flatten() {
                        if positions.get(j).is_some_and(|&q| dist2(p, q) <= d * d) && found.is_none_or(|f| j < f) {
                            found = Some(j);
                        }
                    }
                }
            }
        }
        match found {
            Some(j) => {
                rep.insert(i, j);
            }
            None => {
                grid.entry(c).or_default().push(i);
                rep.insert(i, i);
            }
        }
    }
    rep
}

/// The open edges of `triangles`: undirected sides on one triangle only, as (low, high) with
/// the triangle and its side.
fn open_edges(triangles: &[[u32; 3]]) -> Vec<((u32, u32), usize, usize)> {
    let mut sides: BTreeMap<(u32, u32), Vec<(usize, usize)>> = BTreeMap::new();
    for (t, tri) in triangles.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            sides.entry((a.min(b), a.max(b))).or_default().push((t, k));
        }
    }
    sides.into_iter().filter_map(|(e, on)| if let [(t, k)] = on.as_slice() { Some((e, *t, *k)) } else { None }).collect()
}

/// `mesh` with its gaps under `tolerance` (mm) closed (see the module docs).
pub fn stitch(mesh: &SurfaceMesh, tolerance: f64) -> Result<(SurfaceMesh, StitchReport)> {
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        // Its value is in the geometry's millimetres, which the caller may not show: not named.
        return Err(Error::Invalid("the stitch tolerance must be a finite, positive length".into()));
    }
    if mesh.positions.len() > MAX_VERTICES {
        return Err(Error::Invalid(format!("a surface stitched has at most {MAX_VERTICES} vertices")));
    }
    let n = u32::try_from(mesh.positions.len()).map_err(|_| Error::Invalid("too many vertices".into()))?;
    if mesh.triangles.iter().flatten().any(|&v| v >= n) || mesh.tri_face.len() != mesh.triangles.len() {
        return Err(Error::Invalid("the surface has a triangle with a vertex out of range, or triangles without tags".into()));
    }
    let mut out = mesh.clone();
    let mut report = StitchReport::default();
    // Copies of one point (within the surface's own weld distance) are one point first.
    let all: Vec<usize> = (0..mesh.positions.len()).collect();
    let same = cluster(&mesh.positions, &all, mesh.weld_distance.max(0.0));
    let canon = |v: u32, map: &BTreeMap<usize, usize>| map.get(&(v as usize)).map_or(v, |&r| r as u32);
    for t in &mut out.triangles {
        *t = t.map(|v| canon(v, &same));
    }
    let open = open_edges(&out.triangles);
    report.open_before = open.len();
    // 1. Open points within the tolerance are one point.
    let mut points: Vec<usize> = open.iter().flat_map(|((a, b), _, _)| [*a as usize, *b as usize]).collect();
    points.sort_unstable();
    points.dedup();
    let snap = cluster(&out.positions, &points, tolerance);
    report.merged = snap.iter().filter(|(i, r)| i != r).count();
    for t in &mut out.triangles {
        *t = t.map(|v| canon(v, &snap));
    }
    let before = out.triangles.len();
    let keep: Vec<bool> = out.triangles.iter().map(|t| t[0] != t[1] && t[1] != t[2] && t[0] != t[2]).collect();
    let mut k = keep.iter();
    out.triangles.retain(|_| k.next().copied().unwrap_or(false));
    let mut k = keep.iter();
    out.tri_face.retain(|_| k.next().copied().unwrap_or(false));
    report.removed = before - out.triangles.len();
    // 2. T-junctions: an open point on the inside of another open edge splits its triangle.
    for _ in 0..MAX_PASSES {
        let open = open_edges(&out.triangles);
        let mut points: Vec<u32> = open.iter().flat_map(|((a, b), _, _)| [*a, *b]).collect();
        points.sort_unstable();
        points.dedup();
        // The open points in a grid with cells as large as the open edges are long on average
        // (and no smaller than the tolerance), so an edge looks at the few cells around it.
        let lengths: f64 =
            open.iter().filter_map(|((a, b), _, _)| Some(dist2(*out.positions.get(*a as usize)?, *out.positions.get(*b as usize)?).sqrt())).sum();
        let cell = (lengths / open.len().max(1) as f64).max(tolerance);
        let key = |p: [f64; 3]| p.map(|c| (c / cell).floor() as i64);
        let mut grid: BTreeMap<[i64; 3], Vec<u32>> = BTreeMap::new();
        for &v in &points {
            if let Some(&p) = out.positions.get(v as usize) {
                grid.entry(key(p)).or_default().push(v);
            }
        }
        let mut visits = 0usize;
        // The first edge (in order) each triangle side is split at, by the point nearest its start.
        let mut splits: BTreeMap<usize, (usize, u32)> = BTreeMap::new();
        for ((a, b), t, side) in &open {
            let (Some(&pa), Some(&pb)) = (out.positions.get(*a as usize), out.positions.get(*b as usize)) else { continue };
            let ab = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
            let l2 = dist2(pa, pb);
            if !(l2 > 0.0) {
                continue;
            }
            let mut best: Option<(f64, u32)> = None;
            let lo = key([0, 1, 2].map(|k| pa[k].min(pb[k]) - tolerance));
            let hi = key([0, 1, 2].map(|k| pa[k].max(pb[k]) + tolerance));
            let span = (0..3).map(|k| usize::try_from(hi[k] - lo[k] + 1).unwrap_or(usize::MAX)).fold(1usize, usize::saturating_mul);
            visits = visits.saturating_add(span);
            if visits > MAX_CELL_VISITS {
                return Err(Error::Invalid(format!("the open edges are too many or too uneven to stitch ({} open edges)", open.len())));
            }
            let near = (lo[0]..=hi[0])
                .flat_map(|x| (lo[1]..=hi[1]).flat_map(move |y| (lo[2]..=hi[2]).map(move |z| [x, y, z])))
                .flat_map(|c| grid.get(&c).into_iter().flatten().copied());
            for v in near {
                if v == *a || v == *b {
                    continue;
                }
                let Some(&p) = out.positions.get(v as usize) else { continue };
                let s = ((p[0] - pa[0]) * ab[0] + (p[1] - pa[1]) * ab[1] + (p[2] - pa[2]) * ab[2]) / l2;
                if !(s > INSIDE && s < 1.0 - INSIDE) {
                    continue;
                }
                let q = [pa[0] + s * ab[0], pa[1] + s * ab[1], pa[2] + s * ab[2]];
                if dist2(p, q) <= tolerance * tolerance && best.is_none_or(|(bs, _)| s < bs) {
                    best = Some((s, v));
                }
            }
            if let Some((_, v)) = best {
                splits.entry(*t).or_insert((*side, v));
            }
        }
        if splits.is_empty() {
            break;
        }
        for (t, (side, v)) in splits {
            let (Some(&tri), Some(&tag)) = (out.triangles.get(t), out.tri_face.get(t)) else { continue };
            let (a, b, c) = (tri[side], tri[(side + 1) % 3], tri[(side + 2) % 3]);
            if let Some(slot) = out.triangles.get_mut(t) {
                *slot = [a, v, c];
            }
            out.triangles.push([v, b, c]);
            out.tri_face.push(tag);
            report.split += 1;
        }
    }
    report.open_after = open_edges(&out.triangles).len();
    Ok((out, report))
}

#[doc(hidden)]
/// The unit cube's surface, each face with its own four corners (a triangle soup by face),
/// tagged by face; the top face (+z) raised by `lift` and, with `t_junction`, its +x side
/// split at its middle on the side face only.
pub fn test_cube(lift: f64, t_junction: bool) -> SurfaceMesh {
    let mut m = SurfaceMesh::default();
    // Each face: four corners counter-clockwise from outside.
    let faces: [[[f64; 3]; 4]; 6] = [
        [[0., 0., 0.], [0., 1., 0.], [1., 1., 0.], [1., 0., 0.]],
        [[0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]],
        [[0., 0., 0.], [1., 0., 0.], [1., 0., 1.], [0., 0., 1.]],
        [[0., 1., 0.], [0., 1., 1.], [1., 1., 1.], [1., 1., 0.]],
        [[0., 0., 0.], [0., 0., 1.], [0., 1., 1.], [0., 1., 0.]],
        [[1., 0., 0.], [1., 1., 0.], [1., 1., 1.], [1., 0., 1.]],
    ];
    for (f, q) in faces.iter().enumerate() {
        let base = m.positions.len() as u32;
        let q = q.map(|p| if f == 1 { [p[0], p[1], p[2] + lift] } else { p });
        m.positions.extend(q);
        if t_junction && f == 5 {
            // The +x face's top side has a point at its middle: a fan from it.
            m.positions.push([1.0, 0.5, 1.0]);
            let mid = base + 4;
            for t in [[base, base + 1, mid], [base + 1, base + 2, mid], [base + 3, base, mid]] {
                m.triangles.push(t);
                m.tri_face.push(f as u32);
            }
            continue;
        }
        for t in [[base, base + 1, base + 2], [base, base + 2, base + 3]] {
            m.triangles.push(t);
            m.tri_face.push(f as u32);
        }
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_gap_under_the_tolerance_closes() {
        let (m, r) = stitch(&test_cube(0.01, false), 0.05).unwrap();
        assert!(r.open_before > 0 && r.open_after == 0, "{r:?}");
        assert_eq!(r.merged, 4, "the raised top's corners");
        assert!((m.enclosed_volume() - 1.0).abs() < 0.02, "{}", m.enclosed_volume());
        assert_eq!(m.triangles.len(), m.tri_face.len());
        // Over the tolerance: left open.
        let (_, r) = stitch(&test_cube(0.1, false), 0.05).unwrap();
        assert_eq!(r.open_after, r.open_before);
        // A closed surface is unchanged.
        let (m, r) = stitch(&test_cube(0.0, false), 0.05).unwrap();
        assert_eq!((r.open_before, r.open_after, r.merged, r.split, r.removed), (0, 0, 0, 0, 0));
        assert_eq!(m.triangles.len(), 12);
    }

    #[test]
    fn a_t_junction_splits_the_triangle_across() {
        let (m, r) = stitch(&test_cube(0.0, true), 1e-6).unwrap();
        assert!(r.open_before > 0, "{r:?}");
        assert_eq!((r.open_after, r.split), (0, 1), "{r:?}");
        assert_eq!(m.tri_face.iter().filter(|&&f| f == 1).count(), 3, "the top's triangle split in two");
        assert!((m.enclosed_volume() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn hostile_input_is_an_error() {
        let c = test_cube(0.0, false);
        for t in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(stitch(&c, t).is_err());
        }
        let mut bad = c.clone();
        bad.triangles.push([0, 1, 999]);
        bad.tri_face.push(0);
        assert!(stitch(&bad, 0.1).is_err());
        let mut untagged = c;
        untagged.tri_face.pop();
        assert!(stitch(&untagged, 0.1).is_err());
    }

    /// A strip of 2000 squares along x, tessellated apart from a strip beside it whose points
    /// sit half-way between the first's (T-junctions all along a 4000-edge seam), closes
    /// quickly through the grid.
    #[test]
    fn a_long_seam_of_t_junctions_closes() {
        let n = 2000u32;
        let mut m = SurfaceMesh::default();
        // Strip A: y in [0, 1], points at x = 0..=n; strip B: y in [-1, 0], its seam points at
        // half-steps, its far side at whole steps.
        for i in 0..=n {
            m.positions.push([f64::from(i), 0.0, 0.0]);
            m.positions.push([f64::from(i), 1.0, 0.0]);
        }
        for i in 0..n {
            let (a, b, c, d) = (2 * i, 2 * i + 2, 2 * i + 3, 2 * i + 1);
            m.triangles.extend([[a, b, c], [a, c, d]]);
            m.tri_face.extend([0, 0]);
        }
        let base = m.positions.len() as u32;
        m.positions.push([0.0, 0.0, 0.0]);
        for i in 0..n {
            m.positions.push([f64::from(i) + 0.5, 0.0, 0.0]);
        }
        m.positions.push([f64::from(n), 0.0, 0.0]);
        let far = m.positions.len() as u32;
        for i in 0..=n {
            m.positions.push([f64::from(i), -1.0, 0.0]);
        }
        // Seam points s_0..s_{n+1} (0, 0.5, 1.5, …, n), far points f_0..f_n: a fan per step.
        for i in 0..n {
            let (s0, s1, f0, f1) = (base + i, base + i + 1, far + i, far + i + 1);
            m.triangles.extend([[s0, f0, f1], [s0, f1, s1]]);
            m.tri_face.extend([1, 1]);
        }
        let (s0, s1, f) = (base + n, base + n + 1, far + n);
        m.triangles.push([s0, f, s1]);
        m.tri_face.push(1);
        let started = std::time::Instant::now();
        let (out, r) = stitch(&m, 1e-6).unwrap();
        assert!(r.open_before > 4000, "{r:?}");
        assert!(r.split >= n as usize, "{r:?}");
        assert!(out.triangles.len() > m.triangles.len());
        assert!(started.elapsed().as_secs_f64() < 5.0, "{:?}", started.elapsed());
    }
}
