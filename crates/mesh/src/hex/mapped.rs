//! Mapped meshing of a logical cube: six faces with one loop of four edges each, three edges at
//! each of eight corners.
//!
//! The corners get logical coordinates `(i, j, k) ∈ {0, 1}³` from the edge graph; the four edges
//! along each logical direction share one division count (interval assignment). Edges are divided
//! on the geometry (nodes on both faces of each edge), each face gets the transfinite grid of its
//! four edges projected onto its surface (Gordon and Hall 1973), and the volume is the trivariate
//! transfinite interpolation of the six face grids: faces blended linearly, minus edges blended
//! bilinearly, plus corners blended trilinearly (W. J. Gordon, L. C. Thiel, *Transfinite
//! mappings and their application to grid generation*, Applied Mathematics and Computation
//! 10–11:171–233, 1982).

use std::collections::BTreeMap;

use super::brep::Brep;
use super::interval::Intervals;
use super::{Ctx, Raw, project_to_face, project_to_faces};
use crate::quad::{lattice_index, transfinite_faces, transfinite_volume};
use crate::size::ideal_divisions;
use crate::vec3::{P, det, sub};
use crate::{Error, Result};

/// The most nodes of a mapped block.
pub(crate) const MAX_NODES: usize = 1 << 26;

/// A logical cube found in a B-rep: its corners, its edges by logical direction and its faces.
pub(crate) struct Cube {
    /// The surface vertex at each corner, `corner[i][j][k]`, right-handed.
    pub corner: [[[u32; 2]; 2]; 2],
    /// The four edges along each logical direction `d`, oriented from the low corner to the high
    /// one, at the other two coordinates `(s, t) ∈ {0, 1}²` in the order (0, 0), (0, 1), (1, 0),
    /// (1, 1).
    pub groups: [[(usize, bool); 4]; 3],
    /// The B-rep face (index into `faces`) on each side of each direction: `faces[d][side]`.
    pub faces: [[usize; 2]; 3],
}

impl Cube {
    /// The corner at logical coordinates `(i, j, k)` (each 0 or 1).
    pub(crate) fn at(&self, i: usize, j: usize, k: usize) -> u32 {
        self.corner[i & 1][j & 1][k & 1]
    }
}

/// Recognize a logical cube: six faces with one loop of four edges each, three edges at each of
/// eight corners. Anything else is [`Error::Unsupported`] with the reason.
pub(crate) fn recognize(brep: &Brep) -> Result<Cube> {
    if brep.faces.len() != 6 || brep.edges.len() != 12 {
        return Err(Error::Unsupported(format!("not a logical cube: {} faces and {} edges (6 and 12 needed)", brep.faces.len(), brep.edges.len())));
    }
    if let Some(f) = brep.faces.iter().find(|f| f.loops.len() != 1 || f.loops.first().map(Vec::len) != Some(4)) {
        return Err(Error::Unsupported(format!("not a logical cube: face {} is not bounded by one loop of four edges", f.id)));
    }
    // Corner graph.
    let mut adj: BTreeMap<u32, Vec<(u32, usize)>> = BTreeMap::new();
    for (e, edge) in brep.edges.iter().enumerate() {
        if edge.closed() {
            return Err(Error::Unsupported("not a logical cube: a closed edge".into()));
        }
        adj.entry(edge.first()).or_default().push((edge.last(), e));
        adj.entry(edge.last()).or_default().push((edge.first(), e));
    }
    if adj.len() != 8 || adj.values().any(|v| v.len() != 3) {
        return Err(Error::Unsupported("not a logical cube: the corners do not each have three edges".into()));
    }
    let neighbours = |v: u32| -> Vec<u32> { adj.get(&v).map(|l| l.iter().map(|x| x.0).collect()).unwrap_or_default() };
    let common = |a: u32, b: u32, not: u32| -> Option<u32> {
        let nb = neighbours(b);
        neighbours(a).into_iter().find(|x| *x != not && nb.contains(x))
    };
    let c000 = *adj.keys().next().ok_or_else(|| Error::Failed("no corner".into()))?;
    let n = neighbours(c000);
    let (Some(&a), Some(&b), Some(&c)) = (n.first(), n.get(1), n.get(2)) else { return Err(Error::Failed("corner".into())) };
    // Right-handed: (a, b, c) along (i, j, k) when the triple product is positive.
    let p = |v: u32| brep.point(v);
    let (a, b, c) = if det(sub(p(a), p(c000)), sub(p(b), p(c000)), sub(p(c), p(c000))) > 0.0 { (a, b, c) } else { (a, c, b) };
    let c110 = common(a, b, c000);
    let c101 = common(a, c, c000);
    let c011 = common(b, c, c000);
    let (Some(c110), Some(c101), Some(c011)) = (c110, c101, c011) else {
        return Err(Error::Unsupported("not a logical cube: the corner graph is not a cube".into()));
    };
    let c111 = common(c110, c101, a).filter(|&x| neighbours(c011).contains(&x));
    let Some(c111) = c111 else { return Err(Error::Unsupported("not a logical cube: the corner graph is not a cube".into())) };
    let corner = [[[c000, c], [b, c011]], [[a, c101], [c110, c111]]];
    let at = |i: usize, j: usize, k: usize| -> u32 { corner[i][j][k] };
    let mut distinct: Vec<u32> = corner.iter().flatten().flatten().copied().collect();
    distinct.sort_unstable();
    distinct.dedup();
    if distinct.len() != 8 {
        return Err(Error::Unsupported("not a logical cube: corners repeat".into()));
    }
    // The edge between two corners, oriented from the first.
    let edge_between = |u: u32, v: u32| -> Option<(usize, bool)> {
        adj.get(&u)?.iter().find(|x| x.0 == v).map(|x| (x.1, brep.edges.get(x.1).is_some_and(|e| e.first() == u)))
    };
    let mut groups = [[(0usize, true); 4]; 3];
    for (d, group) in groups.iter_mut().enumerate() {
        for (m, slot) in group.iter_mut().enumerate() {
            let (s, t) = (m / 2, m % 2);
            let (lo, hi) = match d {
                0 => (at(0, s, t), at(1, s, t)),
                1 => (at(s, 0, t), at(s, 1, t)),
                _ => (at(s, t, 0), at(s, t, 1)),
            };
            *slot = edge_between(lo, hi).ok_or_else(|| Error::Unsupported("not a logical cube: a corner pair has no edge".into()))?;
        }
    }
    // The face bounded by the four edges of each side.
    let mut faces = [[0usize; 2]; 3];
    for (d, row) in faces.iter_mut().enumerate() {
        for (side, slot) in row.iter_mut().enumerate() {
            let ring = face_ring(d, side).map(|[i, j, k]| at(i, j, k));
            let mut bounding: Vec<usize> = (0..4).filter_map(|m| edge_between(ring[m], ring[(m + 1) % 4]).map(|x| x.0)).collect();
            bounding.sort_unstable();
            *slot = brep
                .faces
                .iter()
                .position(|f| {
                    let mut es: Vec<usize> = f.loops.iter().flatten().map(|x| x.0).collect();
                    es.sort_unstable();
                    es == bounding
                })
                .ok_or_else(|| Error::Unsupported("not a logical cube: a side's four edges bound no face".into()))?;
        }
    }
    Ok(Cube { corner, groups, faces })
}

/// The corners of side `side` of direction `d` in the order (0, 0), (1, 0), (1, 1), (0, 1) of its
/// two other directions (lower first).
pub(crate) fn face_ring(d: usize, side: usize) -> [[usize; 3]; 4] {
    let (da, db) = other(d);
    [(0, 0), (1, 0), (1, 1), (0, 1)].map(|(a, b)| {
        let mut c = [0usize; 3];
        c[d] = side;
        c[da] = a;
        c[db] = b;
        c
    })
}

/// The two directions other than `d`, lower first.
pub(crate) fn other(d: usize) -> (usize, usize) {
    match d {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    }
}

pub(crate) fn mesh(ctx: &Ctx) -> Result<Raw> {
    let brep = ctx.brep;
    let cube = recognize(brep)?;
    let groups = cube.groups;
    let mut ia = Intervals::default();
    for e in 0..brep.edges.len() {
        ctx.edge_var(&mut ia, e, ideal_divisions(ctx.field, &brep.edge_points(e, true))?);
    }
    for g in &groups {
        let first = g.first().map(|x| x.0).unwrap_or(0);
        for &(e, _) in g.iter().skip(1) {
            ia.equal.push((vec![first], vec![e]));
        }
    }
    let counts = ia.solve()?;
    let count = |g: &[(usize, bool); 4]| g.first().and_then(|x| counts.get(x.0)).copied().unwrap_or(1);
    let (ni, nj, nk) = (count(&groups[0]), count(&groups[1]), count(&groups[2]));
    let total = (ni + 1).checked_mul(nj + 1).and_then(|x| x.checked_mul(nk + 1)).filter(|&x| x <= MAX_NODES);
    let total = total.ok_or_else(|| Error::Limit(format!("a {ni} × {nj} × {nk} block exceeds {MAX_NODES} nodes")))?;
    if ni * nj * nk > ctx.opts.max_elements {
        return Err(Error::Limit(format!("the mesh would have {} hexahedra; the limit is {}", ni * nj * nk, ctx.opts.max_elements)));
    }
    // A biased seed's spacing carries to the unseeded edges parallel to it: the four edges of a
    // direction, all low corner to high.
    let ties: Vec<[(usize, bool); 2]> = groups.iter().flat_map(|g| g.iter().skip(1).map(move |&u| [g[0], u])).collect();
    let carried = ctx.carried(&ties);
    // Edge node positions, low corner to high.
    let edge_pts = |(e, fwd): (usize, bool), n: usize| -> Result<Vec<P>> {
        let pts = ctx.divide_edge(e, fwd, n, &carried)?;
        let faces = brep.edges.get(e).map(|x| x.faces.to_vec()).unwrap_or_default();
        let last = pts.len().saturating_sub(1);
        Ok(pts.iter().enumerate().map(|(k, q)| if k == 0 || k == last { *q } else { project_to_faces(ctx, &faces, *q) }).collect())
    };
    let idx = |i: usize, j: usize, k: usize| lattice_index([ni, nj, nk], [i, j, k]);
    let mut x: Vec<Option<P>> = vec![None; total];
    // Edges: E[d][s][t] along direction d at the other two coordinates (s, t) ∈ {0, 1}².
    for (d, g) in groups.iter().enumerate() {
        let n = [ni, nj, nk][d];
        for (m, &use_) in g.iter().enumerate() {
            let (s, t) = (m / 2, m % 2);
            for (q, pt) in edge_pts(use_, n)?.into_iter().enumerate() {
                let (i, j, k) = match d {
                    0 => (q, s * nj, t * nk),
                    1 => (s * ni, q, t * nk),
                    _ => (s * ni, t * nj, q),
                };
                if let Some(slot) = x.get_mut(idx(i, j, k)) {
                    slot.get_or_insert(pt);
                }
            }
        }
    }
    // Faces: grids from their edges, projected onto the B-rep face their four edges bound.
    let n = [ni, nj, nk];
    transfinite_faces(&mut x, n, None, |d, side, grid| {
        let (da, db) = other(d);
        let face = cube.faces.get(d).and_then(|r| r.get(side)).copied();
        let (na, nb) = (n[da], n[db]);
        for b in 1..nb {
            for a in 1..na {
                if let Some(q) = grid.get_mut(b * (na + 1) + a) {
                    *q = face.and_then(|f| project_to_face(ctx, f, *q)).unwrap_or(*q);
                }
            }
        }
        Ok(())
    })?;
    // Volume: trivariate transfinite interpolation, blended by index.
    let uniform = |n: usize| -> Vec<f64> { (0..=n).map(|q| q as f64 / n as f64).collect() };
    let (bi, bj, bk) = (uniform(ni), uniform(nj), uniform(nk));
    transfinite_volume(&mut x, [ni, nj, nk], [&bi, &bj, &bk])?;
    let nodes: Vec<P> = x.into_iter().map(|p| p.ok_or_else(|| Error::Failed("mapped grid node missing".into()))).collect::<Result<_>>()?;
    let id = |i: usize, j: usize, k: usize| idx(i, j, k) as u32;
    let mut hexes = Vec::with_capacity(ni * nj * nk);
    for k in 0..nk {
        for j in 0..nj {
            for i in 0..ni {
                hexes.push([
                    id(i, j, k),
                    id(i + 1, j, k),
                    id(i + 1, j + 1, k),
                    id(i, j + 1, k),
                    id(i, j, k + 1),
                    id(i + 1, j, k + 1),
                    id(i + 1, j + 1, k + 1),
                    id(i, j + 1, k + 1),
                ]);
            }
        }
    }
    Ok(Raw { nodes, hexes, description: crate::message![format!("mapped block of {ni} × {nj} × {nk}")] })
}
