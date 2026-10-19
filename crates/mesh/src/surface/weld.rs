//! Vertex welding in input order: each vertex joins the earliest representative within the
//! tolerance, or becomes one. Lookups go through a hash map, but nothing is ever iterated out of
//! it, so the result does not depend on hashing.

use std::collections::HashMap;

use super::{Box3, MAX_RELATIVE_WELD, index};
use crate::predicates::MIN_MAGNITUDE;
use crate::vec3::{dot, sub};
use crate::{Error, Result};

pub(crate) struct Welded {
    pub points: Vec<[f64; 3]>,
    pub triangles: Vec<[u32; 3]>,
    pub source_vertex: Vec<u32>,
    pub source_triangle: Vec<u32>,
}

/// Weld the vertices `triangles` use (indices already checked) within the larger of `relative`
/// × the bounding-box diagonal and `absolute`, snap coordinates below [`MIN_MAGNITUDE`] to zero
/// and drop triangles that collapse. A distance above [`MAX_RELATIVE_WELD`] of the diagonal is an
/// error; one below 1e-15 of it welds identical coordinates only.
pub(crate) fn weld(positions: &[[f64; 3]], triangles: &[[u32; 3]], relative: f64, absolute: f64) -> Result<Welded> {
    let mut used = vec![false; positions.len()];
    for v in triangles.iter().flatten() {
        if let Some(u) = used.get_mut(*v as usize) {
            *u = true;
        }
    }
    let referenced = || positions.iter().zip(&used).enumerate().filter(|(_, (_, u))| **u).map(|(v, (p, _))| (v, *p));
    let mut bounds = Box3::EMPTY;
    for (_, p) in referenced() {
        bounds.add_point(p);
    }
    let diagonal = bounds.diagonal();
    let tol = (relative * diagonal).max(absolute);
    if tol > MAX_RELATIVE_WELD * diagonal {
        return Err(Error::invalid(crate::message![
            "weld distance ",
            crate::message::length(tol),
            format!(" is more than {MAX_RELATIVE_WELD} of the surface's size ("),
            crate::message::length(diagonal),
            "): it would merge distinct vertices"
        ]));
    }

    let mut lookup = if tol >= 1e-15 * diagonal && tol > 0.0 {
        Lookup::Grid { tol, lo: bounds.lo, heads: HashMap::new(), next: Vec::new() }
    } else {
        Lookup::Exact(HashMap::new())
    };
    let mut map = vec![u32::MAX; positions.len()];
    let mut points: Vec<[f64; 3]> = Vec::new();
    let mut source_vertex = Vec::new();
    for (v, p) in referenced() {
        let w = match lookup.find(p, &points) {
            Some(w) => w,
            None => {
                let w = index(points.len());
                points.push(p);
                source_vertex.push(index(v));
                lookup.insert(p, w);
                w
            }
        };
        if let Some(m) = map.get_mut(v) {
            *m = w;
        }
    }
    for x in points.iter_mut().flatten() {
        if x.abs() < MIN_MAGNITUDE {
            *x = 0.0;
        }
    }

    let mut out = Vec::with_capacity(triangles.len());
    let mut source_triangle = Vec::with_capacity(triangles.len());
    for (t, tri) in triangles.iter().enumerate() {
        let w = tri.map(|v| map.get(v as usize).copied().unwrap_or(u32::MAX));
        if w[0] != w[1] && w[1] != w[2] && w[2] != w[0] {
            out.push(w);
            source_triangle.push(index(t));
        }
    }
    if tol > 0.0 {
        split_t_junctions(&points, &mut out, &mut source_triangle, tol);
    }
    Ok(Welded { points, triangles: out, source_vertex, source_triangle })
}

/// Rounds of T-junction splitting, and the most open edges looked at.
const T_JUNCTION_ROUNDS: usize = 8;
const MAX_OPEN_EDGES: usize = 100_000;

/// Close T-junctions: where a vertex of one face's boundary lies on an edge of the neighbouring
/// face (within `tol`) that the tessellation did not divide there, split that edge's triangle at
/// it. Only edges no other triangle uses are looked at, so a closed surface is left as it is.
/// The standard repair of tessellations whose faces are triangulated independently (Barequet and
/// Kumar 1997, *Repairing CAD models*).
fn split_t_junctions(points: &[[f64; 3]], triangles: &mut Vec<[u32; 3]>, source: &mut Vec<u32>, tol: f64) {
    use std::collections::{BTreeMap, BTreeSet};
    for _ in 0..T_JUNCTION_ROUNDS {
        let mut used: BTreeMap<(u32, u32), u32> = BTreeMap::new();
        for tri in triangles.iter() {
            for k in 0..3 {
                *used.entry((tri[k], tri[(k + 1) % 3])).or_default() += 1;
            }
        }
        // Edges used once, in one direction only.
        let open: BTreeSet<(u32, u32)> = used.iter().filter(|(e, n)| **n == 1 && !used.contains_key(&(e.1, e.0))).map(|(e, _)| *e).collect();
        if open.is_empty() || open.len() > MAX_OPEN_EDGES {
            return;
        }
        let ends: BTreeSet<u32> = open.iter().flat_map(|&(a, b)| [a, b]).collect();
        let at = |v: u32| points.get(v as usize).copied();
        let mut changed = false;
        let mut next = Vec::with_capacity(triangles.len() + 8);
        let mut next_source = Vec::with_capacity(triangles.len() + 8);
        for (tri, &src) in triangles.iter().zip(source.iter()) {
            let mut split = None;
            for k in 0..3 {
                let (a, b, c) = (tri[k], tri[(k + 1) % 3], tri[(k + 2) % 3]);
                if !open.contains(&(a, b)) {
                    continue;
                }
                let (Some(pa), Some(pb)) = (at(a), at(b)) else { continue };
                let ab = sub(pb, pa);
                let len2 = dot(ab, ab);
                if !(len2 > 0.0) {
                    continue;
                }
                let mut on: Vec<(f64, u32)> = ends
                    .iter()
                    .filter(|&&v| v != a && v != b && v != c)
                    .filter_map(|&v| {
                        let p = at(v)?;
                        let t = dot(sub(p, pa), ab) / len2;
                        let q = [pa[0] + ab[0] * t, pa[1] + ab[1] * t, pa[2] + ab[2] * t];
                        let d = sub(p, q);
                        let inside = t * len2.sqrt() > tol && (1.0 - t) * len2.sqrt() > tol;
                        (inside && dot(d, d) <= tol * tol).then_some((t, v))
                    })
                    .collect();
                if !on.is_empty() {
                    on.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
                    split = Some((a, b, c, on));
                    break;
                }
            }
            match split {
                Some((a, b, c, on)) => {
                    let mut prev = a;
                    for (_, v) in on.iter().copied().chain(std::iter::once((1.0, b))) {
                        next.push([prev, v, c]);
                        next_source.push(src);
                        prev = v;
                    }
                    changed = true;
                }
                None => {
                    next.push(*tri);
                    next_source.push(src);
                }
            }
        }
        *triangles = next;
        *source = next_source;
        if !changed {
            return;
        }
    }
}

/// Representatives found so far.
enum Lookup {
    /// Identical coordinates only (keyed by bits; `+ 0.0` turns -0 into +0).
    Exact(HashMap<[u64; 3], u32>),
    /// Cells of side `tol`: a vertex's representatives are in its cell or the 26 around it.
    /// Coordinates over the cell size stay below 1e15 (the distance is at least 1e-15 of the
    /// diagonal), exact in i64. Each cell keeps the head of a chain through `next`.
    Grid { tol: f64, lo: [f64; 3], heads: HashMap<[i64; 3], u32>, next: Vec<u32> },
}

impl Lookup {
    fn cell(tol: f64, lo: [f64; 3], p: [f64; 3]) -> [i64; 3] {
        [0, 1, 2].map(|k| ((p[k] - lo[k]) / tol).floor() as i64)
    }

    /// The earliest representative within the distance of `p`.
    fn find(&self, p: [f64; 3], points: &[[f64; 3]]) -> Option<u32> {
        match self {
            Lookup::Exact(seen) => seen.get(&p.map(|x| (x + 0.0).to_bits())).copied(),
            Lookup::Grid { tol, lo, heads, next } => {
                let c = Self::cell(*tol, *lo, p);
                let mut best: Option<u32> = None;
                for dx in -1..=1_i64 {
                    for dy in -1..=1_i64 {
                        for dz in -1..=1_i64 {
                            let key = [c[0].saturating_add(dx), c[1].saturating_add(dy), c[2].saturating_add(dz)];
                            let mut w = heads.get(&key).copied().unwrap_or(u32::MAX);
                            while let Some(&q) = points.get(w as usize) {
                                let d = sub(p, q);
                                if dot(d, d) <= tol * tol && best.is_none_or(|b| w < b) {
                                    best = Some(w);
                                }
                                w = next.get(w as usize).copied().unwrap_or(u32::MAX);
                            }
                        }
                    }
                }
                best
            }
        }
    }

    /// Record `p` as representative `w` (the next index).
    fn insert(&mut self, p: [f64; 3], w: u32) {
        match self {
            Lookup::Exact(seen) => {
                seen.insert(p.map(|x| (x + 0.0).to_bits()), w);
            }
            Lookup::Grid { tol, lo, heads, next } => {
                let head = heads.insert(Self::cell(*tol, *lo, p), w);
                next.push(head.unwrap_or(u32::MAX));
            }
        }
    }
}
