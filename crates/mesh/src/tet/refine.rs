//! Constrained Delaunay refinement to a size field.
//!
//! A tetrahedron is bad when its radius–edge ratio exceeds a bound `B` (about 2) or its
//! circumradius exceeds the local target size; its circumcentre is inserted (Delaunay insertion,
//! [`ops::insert_delaunay`]). A circumcentre outside the domain, or inside the diametral sphere of
//! a boundary face (it *encroaches* the face), is not inserted; the encroached boundary face is
//! split instead, at its circumcentre when that lies inside it and otherwise at the midpoint of
//! its longest edge (which splits a feature edge when the edge is one), and every boundary point
//! is projected onto its B-rep face or feature curve ([`super::geometry`]). A feature curve is
//! split at a point of the exact curve where its faces' surfaces are known, not of the
//! tessellation's chord, so a curve refined below the facet size stays on the curve
//! (C. Boivin, C. Ollivier-Gooch, *Guaranteed-quality triangular mesh generation for domains
//! with curved boundaries*, IJNME 55(10):1185–1213, 2002; S.-W. Cheng, T. K. Dey, J. R.
//! Shewchuk, *Delaunay Mesh Generation*, CRC 2012, on piecewise-smooth complexes).
//!
//! Sources: J. R. Shewchuk, *Tetrahedral Mesh Generation by Delaunay Refinement*, SoCG 1998
//! (radius–edge refinement with encroachment rules, termination for B = 2 and input angles of 90°
//! or more); H. Si, *Constrained Delaunay tetrahedral mesh generation and refinement*, FEAD
//! 46(1–2), 2010 (refinement of a constrained Delaunay tetrahedralization with small input
//! angles). CAD boundaries have small angles, so termination is guaranteed here by bounds
//! instead: a boundary face smaller than a fraction of the local size is not split, a
//! tetrahedron whose point could not be inserted is not retried, and the number of insertions
//! and tetrahedra are capped (the stop is reported).

use crate::predicates::{Point3, snap_to_domain};
use crate::size::SizeField;
use crate::{Error, Result};

use super::complex::{IdMap, NONE, TetComplex, VertexKind, Walk, face_toward, sorted2, sorted4};
use super::geometry::BoundaryGeometry;
use super::measure::{self, P};
use super::ops::{self, Site};

/// What the boundary insertions project onto.
pub struct Context<'a> {
    pub size: &'a dyn SizeField,
    pub geometry: &'a dyn BoundaryGeometry,
    /// The B-rep faces of each feature curve (by curve id).
    pub curve_faces: &'a [[u32; 2]],
}

impl Context<'_> {
    /// The local target size, or `fallback` when the field misbehaves.
    pub fn h(&self, p: P, fallback: f64) -> f64 {
        let h = self.size.size(p);
        if h > 0.0 && h.is_finite() { h } else { fallback }
    }

    /// `p` moved onto the boundary entity of a vertex of `kind`, as an admissible point.
    pub fn place(&self, kind: VertexKind, p: P) -> Option<Point3> {
        let q = match kind {
            VertexKind::Interior => p,
            VertexKind::Face(f) => self.geometry.project_face(f, p)?,
            VertexKind::Curve(c) => self.geometry.project_curve(*self.curve_faces.get(c as usize)?, p)?,
            VertexKind::Corner => return None,
        };
        Point3::try_from(q.map(snap_to_domain)).ok()
    }
}

/// Refinement settings.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefineOptions {
    /// The radius–edge bound `B`.
    pub radius_edge: f64,
    /// A tetrahedron is too large when its circumradius exceeds this times the local size (the
    /// regular tetrahedron with edge `h` has circumradius 0.61 h).
    pub size_factor: f64,
    /// Boundary faces whose circumradius is below this times the local size are not split.
    pub min_face: f64,
    /// A tetrahedron is refined for its shape (radius–edge ratio) only when its circumradius is
    /// at least this times the local size; smaller ones are left to smoothing and improvement,
    /// so badly shaped input near the boundary does not cascade into refinement far below the
    /// target size.
    pub min_shape_radius: f64,
    /// Most vertices inserted.
    pub max_insertions: usize,
}

impl Default for RefineOptions {
    fn default() -> Self {
        RefineOptions { radius_edge: 2.0, size_factor: 0.95, min_face: 0.5, min_shape_radius: 0.5, max_insertions: 2_000_000 }
    }
}

/// What a refinement did.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RefineStats {
    pub interior: usize,
    pub boundary: usize,
    pub skipped: usize,
    /// Why it stopped before every tetrahedron was good, if it did.
    pub stopped: Option<String>,
}

/// Refinement never puts a vertex closer than this times the local size to another.
pub const MIN_GAP: f64 = 0.05;

/// Most tetrahedra a Delaunay cavity may hold.
pub const MAX_CAVITY: usize = 4096;
/// Most steps of a point location walk.
const MAX_WALK: usize = 1 << 16;

/// The boundary face of a tetrahedron whose diametral sphere holds `p`, if any (the largest).
fn encroached(cx: &TetComplex, tets: &[u32], p: P) -> Result<Option<(u32, usize)>> {
    let mut best: Option<(f64, u32, usize)> = None;
    for &t in tets {
        let vs = cx.tet(t)?;
        let ns = cx.adj(t)?;
        for k in 0..4 {
            if ns[k] != NONE {
                continue;
            }
            let f = face_toward(vs, k);
            let (a, b, c) = (cx.xyz(f[0])?, cx.xyz(f[1])?, cx.xyz(f[2])?);
            let Some(o) = measure::triangle_circumcenter(a, b, c) else { continue };
            let r = measure::dist(o, a);
            if measure::dist(o, p) < r && best.is_none_or(|(br, _, _)| r > br) {
                best = Some((r, t, k));
            }
        }
    }
    Ok(best.map(|(_, t, k)| (t, k)))
}

/// The tetrahedra whose circumsphere holds `p`, grown from `t` across interior faces.
fn conflicts(cx: &TetComplex, t: u32, p: &Point3) -> Result<Vec<u32>> {
    let mut out = vec![t];
    let mut seen: IdMap<u32, ()> = IdMap::default();
    seen.insert(t, ());
    let mut i = 0;
    while let Some(&u) = out.get(i) {
        i += 1;
        if out.len() > MAX_CAVITY {
            break;
        }
        for n in cx.adj(u)? {
            if n != NONE && !seen.contains_key(&n) {
                seen.insert(n, ());
                if cx.in_sphere(cx.tet(n)?, p)? {
                    out.push(n);
                }
            }
        }
    }
    Ok(out)
}

/// Split boundary face `k` of tetrahedron `t`: at its circumcentre when inside it, else at the
/// midpoint of its longest edge. Returns the new vertex.
pub fn split_face(cx: &mut TetComplex, ctx: &Context, t: u32, k: usize, min_face: f64) -> Result<Option<u32>> {
    let vs = cx.tet(t)?;
    if cx.adj(t)?.get(k) != Some(&NONE) {
        return Ok(None);
    }
    let tag = cx.tags(t)?.get(k).copied().unwrap_or(NONE);
    let f = face_toward(vs, k);
    let x = [cx.xyz(f[0])?, cx.xyz(f[1])?, cx.xyz(f[2])?];
    let Some(o) = measure::triangle_circumcenter(x[0], x[1], x[2]) else { return Ok(None) };
    let r = measure::dist(o, x[0]);
    let centroid = measure::scale(measure::add(measure::add(x[0], x[1]), x[2]), 1.0 / 3.0);
    if r < min_face * ctx.h(centroid, r) {
        return Ok(None);
    }
    // Barycentric coordinates of o (signed areas against the face normal).
    let n = measure::cross(measure::sub(x[1], x[0]), measure::sub(x[2], x[0]));
    let inside = (0..3).all(|i| {
        let (a, b) = (x[(i + 1) % 3], x[(i + 2) % 3]);
        measure::dot(measure::cross(measure::sub(b, a), measure::sub(o, a)), n) >= 0.0
    });
    if inside && tag != NONE {
        // A point inside the diametral sphere of one of the face's feature edges encroaches it:
        // split the edge instead (Shewchuk 1998).
        for i in 0..3 {
            let (u, w) = (f[i], f[(i + 1) % 3]);
            if cx.curve(u, w).is_some() {
                let (a, b) = (x[i], x[(i + 1) % 3]);
                if measure::dist(o, measure::mid(a, b)) < 0.5 * measure::dist(a, b) {
                    return split_edge(cx, ctx, u, w, tag, min_face);
                }
            }
        }
        let kind = VertexKind::Face(tag);
        let Some(p) = ctx.place(kind, o) else { return Ok(None) };
        // The projection may only move the point a little, and never next to a corner.
        let shortest = (0..3).map(|i| measure::dist(x[i], x[(i + 1) % 3])).fold(f64::INFINITY, f64::min);
        if measure::dist(p.coords(), o) > 0.5 * r || x.iter().any(|c| measure::dist(*c, p.coords()) < 0.25 * shortest) {
            return Ok(None);
        }
        return ops::insert_delaunay_apart(cx, p, kind, Site::Face(t, k), MAX_CAVITY, MIN_GAP * ctx.h(o, r));
    }
    // The longest edge.
    let (mut u, mut w, mut len) = (f[0], f[1], -1.0);
    for i in 0..3 {
        let l = measure::dist(x[i], x[(i + 1) % 3]);
        if l > len {
            (u, w, len) = (f[i], f[(i + 1) % 3], l);
        }
    }
    split_edge(cx, ctx, u, w, tag, min_face)
}

/// Split boundary edge `uw` at its midpoint, projected onto its feature curve or onto B-rep face
/// `tag`. An edge shorter than `2 · min_face` times the local size is not split.
pub fn split_edge(cx: &mut TetComplex, ctx: &Context, u: u32, w: u32, tag: u32, min_face: f64) -> Result<Option<u32>> {
    let m = measure::mid(cx.xyz(u)?, cx.xyz(w)?);
    let len = measure::dist(cx.xyz(u)?, cx.xyz(w)?);
    if len < 2.0 * min_face * ctx.h(m, len) {
        return Ok(None);
    }
    let kind = match cx.curve(u, w) {
        Some(c) => VertexKind::Curve(c),
        None if tag != NONE => VertexKind::Face(tag),
        None => return Ok(None),
    };
    // A midpoint inside the diametral sphere of a nearby feature edge encroaches it: split that
    // edge instead (Shewchuk 1998), or face vertices pile up against the curve, each split
    // halving their distance to it.
    if matches!(kind, VertexKind::Face(_)) {
        let mut near = vec![u, w];
        if let Some(ring) = cx.edge_ring(u, w)? {
            near.extend(ring.verts.first().copied());
            near.extend(ring.verts.last().copied());
        }
        if let Some((a, b)) = encroached_curve_edge(cx, &near, m)? {
            return split_edge(cx, ctx, a, b, tag, min_face);
        }
    }
    let Some(p) = ctx.place(kind, m) else { return Ok(None) };
    // The projection may only move the midpoint a little (a quarter of the edge).
    if measure::dist(p.coords(), m) > 0.25 * len {
        return Ok(None);
    }
    ops::insert_delaunay_apart(cx, p, kind, Site::Edge(u, w), MAX_CAVITY, MIN_GAP * ctx.h(m, len))
}

/// A feature-curve edge at one of the vertices `near` whose diametral sphere holds `p`.
fn encroached_curve_edge(cx: &mut TetComplex, near: &[u32], p: P) -> Result<Option<(u32, u32)>> {
    for &c in near {
        for n in cx.neighbours(c)? {
            if cx.curve(c, n).is_some() {
                let (a, b) = (cx.xyz(c)?, cx.xyz(n)?);
                if measure::dist(p, measure::mid(a, b)) < 0.5 * measure::dist(a, b) {
                    return Ok(Some((c, n)));
                }
            }
        }
    }
    Ok(None)
}

/// Whether a tetrahedron is bad, and its circumcentre.
fn badness(x: &[P; 4], ctx: &Context, o: &RefineOptions) -> Option<P> {
    let c = measure::circumcenter(x)?;
    let r = measure::dist(c, x[0]);
    let (lmin, lmax) = measure::edge_range(x);
    let h = ctx.h(measure::centroid(x), lmax);
    ((r > o.radius_edge * lmin && r > o.min_shape_radius * h) || r > o.size_factor * h).then_some(c)
}

/// Refine until no tetrahedron is bad (see the module docs).
pub fn refine(cx: &mut TetComplex, ctx: &Context, o: &RefineOptions) -> Result<RefineStats> {
    let mut stats = RefineStats::default();
    let mut queue: Vec<u32> = cx.live_tets().map(|(t, _)| t).collect();
    let mut given_up: IdMap<[u32; 4], ()> = IdMap::default();
    let mut head = 0usize;
    let mut inserted = 0usize;
    while head < queue.len() {
        let t = queue.get(head).copied().unwrap_or(NONE);
        head += 1;
        // A round can insert hundreds of thousands of vertices: stop between them when asked.
        if head.is_multiple_of(1024) {
            fieldcraft_model::progress::check()?;
        }
        // Keep the queue's memory bounded: drop the processed prefix now and then.
        if head > 1 << 20 {
            queue.drain(..head);
            head = 0;
        }
        if !cx.is_live(t) {
            continue;
        }
        let vs = cx.tet(t)?;
        let key = sorted4(vs);
        if given_up.contains_key(&key) {
            continue;
        }
        let x = cx.coords(vs)?;
        let Some(c) = badness(&x, ctx, o) else { continue };
        if inserted >= o.max_insertions {
            stats.stopped = Some(format!("refinement stopped at {} inserted vertices", o.max_insertions));
            break;
        }
        if cx.live_count() + 64 >= cx.max_tets() || cx.point_count() + 64 >= cx.max_points() {
            stats.stopped = Some(format!("refinement stopped at {} tetrahedra and {} vertices", cx.live_count(), cx.point_count()));
            break;
        }
        let result = match Point3::try_from(c.map(snap_to_domain)) {
            Err(_) => None,
            Ok(p) => {
                // The boundary face the circumcentre lies beyond or encroaches, if any.
                let (inside, face) = match cx.walk(t, &p, MAX_WALK)? {
                    None => (None, None),
                    Some(Walk::Outside(tb, k)) => (None, Some((tb, k))),
                    Some(Walk::Inside(ti)) => (Some(ti), encroached(cx, &conflicts(cx, ti, &p)?, c)?),
                };
                match (inside, face) {
                    (_, Some((tb, k))) => {
                        let v = split_face(cx, ctx, tb, k, o.min_face)?;
                        stats.boundary += usize::from(v.is_some());
                        v
                    }
                    (Some(ti), None) => {
                        let v = ops::insert_delaunay_apart(cx, p, VertexKind::Interior, Site::Interior(ti), MAX_CAVITY, MIN_GAP * ctx.h(c, 1.0))?;
                        stats.interior += usize::from(v.is_some());
                        v
                    }
                    (None, None) => None,
                }
            }
        };
        match result {
            Some(v) => {
                inserted += 1;
                queue.extend(cx.star(v)?);
                // The bad tetrahedron may survive a boundary split: look at it again.
                if cx.is_live(t) {
                    queue.push(t);
                }
            }
            None => {
                stats.skipped += 1;
                given_up.insert(key, ());
            }
        }
    }
    Ok(stats)
}

/// The feature curves of a boundary: an edge between two boundary faces is a feature edge when
/// their tags differ, or when their normals turn by more than `feature_cos` (cosine) and their
/// tag is not `smooth` (by tag; a tag past its end is not). A face whose exact surface is known
/// is smooth inside, so a turn between two of its triangles is the tessellation's own: a needle
/// triangle along a trim, its corners nearly in line, can even face the other way. Returns the
/// curve table (face pairs) and the edges with their curve.
pub fn feature_curves(cx: &TetComplex, feature_cos: f64, smooth: &[bool]) -> Result<(Vec<[u32; 2]>, Vec<((u32, u32), u32)>)> {
    // Boundary edge -> (tag, oriented normal) of each side.
    let mut sides: IdMap<(u32, u32), Vec<(u32, P)>> = IdMap::default();
    for (t, vs) in cx.live_tets() {
        let ns = cx.adj(t)?;
        let tags = cx.tags(t)?;
        for k in 0..4 {
            if ns[k] != NONE {
                continue;
            }
            let f = face_toward(vs, k);
            let n = measure::triangle_normal(cx.xyz(f[0])?, cx.xyz(f[1])?, cx.xyz(f[2])?);
            for i in 0..3 {
                sides.entry(sorted2(f[i], f[(i + 1) % 3])).or_default().push((tags[k], n));
            }
        }
    }
    let mut keys: Vec<(u32, u32)> = sides.keys().copied().collect();
    keys.sort_unstable();
    let mut table: Vec<[u32; 2]> = Vec::new();
    let mut index: IdMap<[u32; 2], u32> = IdMap::default();
    let mut edges = Vec::new();
    for e in keys {
        let s = sides.get(&e).map(Vec::as_slice).unwrap_or_default();
        let [(ta, na), (tb, nb)] = s else {
            return Err(Error::Failed(format!("boundary edge {e:?} has {} faces", s.len())));
        };
        if ta != tb || (measure::dot(*na, *nb) < feature_cos && !smooth.get(*ta as usize).copied().unwrap_or(false)) {
            let pair = [(*ta).min(*tb), (*ta).max(*tb)];
            let id = match index.get(&pair) {
                Some(&id) => id,
                None => {
                    let id = table.len() as u32;
                    table.push(pair);
                    index.insert(pair, id);
                    id
                }
            };
            edges.push((e, id));
        }
    }
    Ok((table, edges))
}

/// Coarsening: collapse edges shorter than `short` times the local size, shortest (relative)
/// first, each onto whichever end it may move to ([`ops::collapse`]); new edges stay below
/// `long` times the size, and the worst tetrahedron around does not drop below
/// `min(its old worst, floor)`. At most `passes` passes. Returns the collapses.
pub fn coarsen(cx: &mut TetComplex, ctx: &Context, short: f64, long: f64, floor: f64, passes: usize) -> Result<usize> {
    let mut done = 0;
    for _ in 0..passes {
        let mut edges: Vec<(f64, u32, u32)> = Vec::new();
        let mut seen: IdMap<(u32, u32), ()> = IdMap::default();
        for (_, vs) in cx.live_tets() {
            for &(i, j) in &super::complex::EDGES {
                let e = sorted2(vs[i], vs[j]);
                if seen.insert(e, ()).is_some() {
                    continue;
                }
                let (a, b) = (cx.xyz(e.0)?, cx.xyz(e.1)?);
                let l = measure::dist(a, b);
                let h = ctx.h(measure::mid(a, b), l);
                if l < short * h {
                    edges.push((l / h, e.0, e.1));
                }
            }
        }
        edges.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)).then(x.2.cmp(&y.2)));
        let mut changed = 0;
        for (_, a, b) in edges {
            if !cx.is_used(a) || !cx.is_used(b) {
                continue;
            }
            let (pa, pb) = (cx.xyz(a)?, cx.xyz(b)?);
            let max_len = long * ctx.h(measure::mid(pa, pb), measure::dist(pa, pb));
            if ops::collapse(cx, a, b, max_len, floor)? || ops::collapse(cx, b, a, max_len, floor)? {
                changed += 1;
            }
        }
        done += changed;
        if changed == 0 {
            break;
        }
    }
    Ok(done)
}
