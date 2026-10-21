//! Local operations on the [`TetComplex`]: edge removal, face removal (2-3 flip and multi-face
//! removal), Delaunay vertex insertion with a star-shaped cavity, and the optimal-cavity vertex
//! insertion of the improvement stage.
//!
//! Every operation is a [`TetComplex::replace`], so it either keeps a valid mesh or changes
//! nothing; the topological ones are hill-climbing steps that are applied only when the worst new
//! tetrahedron is better than the worst one they remove.
//!
//! Sources:
//! - J. R. Shewchuk, *Two Discrete Optimization Algorithms for the Topological Improvement of
//!   Tetrahedral Meshes*, 2002: edge removal by Klincsek's dynamic programming over the ring of
//!   the edge (§2, including boundary edges ending in a 2-2 flip), and multi-face removal (§3).
//! - B. M. Klingner, J. R. Shewchuk, *Aggressive Tetrahedral Mesh Improvement*, IMR 2007, §3.4:
//!   vertex insertion whose cavity maximizes the worst new tetrahedron (the cut of Figure 3).
//! - A. Bowyer 1981, D. F. Watson 1981; Shewchuk, *Tetrahedral Mesh Generation by Delaunay
//!   Refinement*, SoCG 1998: Delaunay insertion; the cavity is shrunk until it is star-shaped from
//!   the new vertex, so a mesh that is no longer exactly Delaunay (after smoothing) stays valid.

use crate::Result;
use crate::predicates::{Point3, Sign};

use super::complex::{Change, EDGES, IdMap, NONE, TetComplex, VertexKind, face_toward, sorted2, sorted3};
use super::measure::{self, P};

/// The quality objective of a tetrahedron given by vertex ids ([`measure::quality`]).
pub fn tet_quality(cx: &TetComplex, vs: [u32; 4]) -> Result<f64> {
    Ok(measure::quality(&cx.coords(vs)?))
}

/// The worst quality of a set of tetrahedra.
pub fn worst(cx: &TetComplex, tets: &[u32]) -> Result<f64> {
    let mut q = f64::INFINITY;
    for &t in tets {
        q = q.min(tet_quality(cx, cx.tet(t)?)?);
    }
    Ok(q)
}

/// The largest ring an edge removal considers (Shewchuk 2002 recommends 7; up to 10).
pub const MAX_RING: usize = 10;

/// The cosine of the largest angle (20°) between the normals of two boundary faces on one B-rep
/// face that a boundary edge removal (ending in a 2-2 flip) or the deletion of a tetrahedron
/// capping them may re-triangulate. The size field holds the faces of a curved B-rep face to
/// about 22.5° of turn per element, so faces of one B-rep face fold by about that much, and the
/// change moves the boundary by an element's chordal error.
pub const BOUNDARY_FLIP_COS: f64 = 0.939_692_620_785_908_4;

/// Klincsek's table: the best worst quality of a triangulation of the chain `v_i … v_j`
/// (with the closing edge `v_i v_j`), and the apex `k` that achieves it.
fn klincsek(m: usize, q: &mut dyn FnMut(usize, usize, usize) -> Result<f64>) -> Result<(Vec<f64>, Vec<usize>)> {
    let mut best = vec![f64::INFINITY; m * m];
    let mut arg = vec![0usize; m * m];
    for i in (0..m.saturating_sub(2)).rev() {
        for j in i + 2..m {
            let mut b = f64::NEG_INFINITY;
            let mut a = i + 1;
            for k in i + 1..j {
                let sub = best.get(i * m + k).copied().unwrap_or(f64::NEG_INFINITY).min(best.get(k * m + j).copied().unwrap_or(f64::NEG_INFINITY));
                if sub <= b {
                    continue;
                }
                let v = sub.min(q(i, k, j)?);
                if v > b {
                    b = v;
                    a = k;
                }
            }
            if let Some(s) = best.get_mut(i * m + j) {
                *s = b;
            }
            if let Some(s) = arg.get_mut(i * m + j) {
                *s = a;
            }
        }
    }
    Ok((best, arg))
}

/// The triangles `(i, k, j)` of the optimal triangulation of `v_i … v_j`.
fn triangles(arg: &[usize], m: usize, i: usize, j: usize, out: &mut Vec<(usize, usize, usize)>) {
    let mut stack = vec![(i, j)];
    while let Some((i, j)) = stack.pop() {
        if j < i + 2 {
            continue;
        }
        let k = arg.get(i * m + j).copied().unwrap_or(i + 1);
        out.push((i, k, j));
        stack.push((i, k));
        stack.push((k, j));
    }
}

/// The two tetrahedra a ring triangle `(v_i, v_k, v_j)` (in ring order) makes with the edge's
/// ends: positive when the ring runs as in [`TetComplex::edge_ring`].
fn ring_tets(a: u32, b: u32, vi: u32, vk: u32, vj: u32) -> [[u32; 4]; 2] {
    [[vi, vj, vk, a], [vi, vk, vj, b]]
}

/// Remove edge `ab` if that improves the worst tetrahedron around it (Shewchuk 2002 §2). A
/// boundary edge is removed only when its two boundary faces lie on the same B-rep face, are
/// nearly coplanar ([`BOUNDARY_FLIP_COS`]) and the edge is not a feature edge; the operation then
/// ends in a 2-2 flip of the boundary. Returns whether the mesh changed.
pub fn remove_edge(cx: &mut TetComplex, a: u32, b: u32) -> Result<bool> {
    if cx.curve(a, b).is_some() {
        return Ok(false);
    }
    let Some(t) = cx.edge_tet(a, b)? else { return Ok(false) };
    remove_edge_at(cx, t, a, b)
}

/// [`remove_edge`] from a tetrahedron `t` that holds the edge.
pub fn remove_edge_at(cx: &mut TetComplex, t: u32, a: u32, b: u32) -> Result<bool> {
    let vs = cx.tet(t)?;
    if cx.curve(a, b).is_some() || !vs.contains(&a) || !vs.contains(&b) {
        return Ok(false);
    }
    let Some(ring) = cx.edge_ring_at(t, a, b)? else { return Ok(false) };
    let m = ring.verts.len();
    if m < 3 || m > MAX_RING + usize::from(!ring.closed) {
        return Ok(false);
    }
    let v = ring.verts.clone();
    let x = |id: u32| cx.xyz(id);
    let (xa, xb) = (x(a)?, x(b)?);
    // An open ring runs from boundary face (a, b, v_0) to (a, b, v_{m−1}): both must lie on one
    // flat B-rep face, and they become (a, v_0, v_{m−1}) and (b, v_0, v_{m−1}) (a 2-2 flip).
    let mut new_boundary = Vec::new();
    if !ring.closed {
        let (first, last) = (ring.tets.first().copied().unwrap_or(NONE), ring.tets.last().copied().unwrap_or(NONE));
        let (v0, vl) = (v.first().copied().unwrap_or(NONE), v.last().copied().unwrap_or(NONE));
        let (Some(t0), Some(t1)) = (boundary_tag(cx, first, [a, b, v0])?, boundary_tag(cx, last, [a, b, vl])?) else { return Ok(false) };
        let n0 = measure::triangle_normal(xa, xb, x(v0)?);
        let n1 = measure::triangle_normal(xa, x(vl)?, xb);
        if t0 != t1 || !(measure::dot(n0, n1) >= BOUNDARY_FLIP_COS) {
            return Ok(false);
        }
        // Oriented into the mesh, as the faces of the new tetrahedra (v_0, v_l, v_k, a) and
        // (v_0, v_k, v_l, b) opposite v_k.
        new_boundary.push(([v0, a, vl], t0));
        new_boundary.push(([v0, vl, b], t0));
    }
    let before = worst(cx, &ring.tets)?;
    let xv = v.iter().map(|&id| x(id)).collect::<Result<Vec<P>>>()?;
    let mut q = |i: usize, k: usize, j: usize| -> Result<f64> {
        let (pi, pk, pj) = match (xv.get(i), xv.get(k), xv.get(j)) {
            (Some(p), Some(q), Some(r)) => (*p, *q, *r),
            _ => return Ok(f64::NEG_INFINITY),
        };
        Ok(measure::quality(&[pi, pj, pk, xa]).min(measure::quality(&[pi, pk, pj, xb])))
    };
    // Triangulate the polygon v_0 … v_{m−1}, closed by the edge v_{m−1} v_0.
    let (best, arg) = klincsek(m, &mut q)?;
    let total = best.get(m - 1).copied().unwrap_or(f64::NEG_INFINITY);
    if !(total > before) {
        return Ok(false);
    }
    let mut tris = Vec::new();
    triangles(&arg, m, 0, m - 1, &mut tris);
    let mut new = Vec::with_capacity(2 * tris.len());
    for (i, k, j) in tris {
        let (vi, vk, vj) = match (v.get(i), v.get(k), v.get(j)) {
            (Some(p), Some(q), Some(r)) => (*p, *q, *r),
            _ => return Ok(false),
        };
        new.extend(ring_tets(a, b, vi, vk, vj));
    }
    let change = Change { old: ring.tets.clone(), new, new_boundary, drop_curves: Vec::new() };
    Ok(cx.replace(&change)?.is_some())
}

/// The tag of the boundary face of tetrahedron `t` with vertices `face`, if it is one.
pub fn boundary_tag(cx: &TetComplex, t: u32, face: [u32; 3]) -> Result<Option<u32>> {
    if t == NONE {
        return Ok(None);
    }
    let vs = cx.tet(t)?;
    let key = sorted3(face);
    let ns = cx.adj(t)?;
    let tags = cx.tags(t)?;
    for k in 0..4 {
        if sorted3(face_toward(vs, k)) == key && ns[k] == NONE {
            return Ok(Some(tags[k]));
        }
    }
    Ok(None)
}

/// Remove the interior face opposite vertex `k` of tetrahedron `t` by the best multi-face removal
/// that includes it (Shewchuk 2002 §3), if that improves the worst tetrahedron. With a single
/// face this is the 2-3 flip. Returns whether the mesh changed.
pub fn remove_face(cx: &mut TetComplex, t: u32, k: usize) -> Result<bool> {
    let vs = cx.tet(t)?;
    let u = cx.adj(t)?.get(k).copied().unwrap_or(NONE);
    if u == NONE {
        return Ok(false);
    }
    let a = vs[k];
    let us = cx.tet(u)?;
    let face = face_toward(vs, k);
    let b = us.iter().copied().find(|v| !face.contains(v)).unwrap_or(NONE);
    if b == NONE || cx.curve(a, b).is_some() || cx.edge_tet(a, b)?.is_some() {
        return Ok(false);
    }
    // Sandwiched faces between a and b: faces f with tetrahedra (a ∪ f) and (b ∪ f). Grow the
    // set from the given face across shared edges while every new tetrahedron (a, b, e) on the
    // ring stays better than the worst removed one; take the best prefix found (a greedy form of
    // Shewchuk's optimal search, which visits the same candidates).
    let mut faces: Vec<[u32; 3]> = vec![face];
    let mut old: Vec<u32> = vec![t, u];
    let mut best: Option<(f64, Vec<u32>, Vec<[u32; 4]>)> = None;
    for _round in 0..8 {
        let removed_q = worst(cx, &old)?;
        // The boundary edges of the face set (edges in exactly one face), oriented as in the face
        // seen from a: each gives the new tetrahedron (a, b, e) positive.
        let mut count: IdMap<(u32, u32), (u32, u32, usize)> = IdMap::default();
        for f in &faces {
            for i in 0..3 {
                let (p, q) = (f[i], f[(i + 1) % 3]);
                let e = count.entry(sorted2(p, q)).or_insert((p, q, 0));
                e.2 += 1;
            }
        }
        let mut ring: Vec<(u32, u32)> = count.values().filter(|e| e.2 == 1).map(|e| (e.0, e.1)).collect();
        ring.sort_unstable();
        let mut new = Vec::with_capacity(ring.len());
        let mut q = f64::INFINITY;
        for &(p, r) in &ring {
            // face (p, r, x) has a on its positive side, so (p, r, a, b)… orient explicitly.
            let cand = [p, r, b, a];
            let cand = if cx.orient(cand)? == Sign::Positive { cand } else { [r, p, b, a] };
            q = q.min(tet_quality(cx, cand)?);
            new.push(cand);
        }
        if q > removed_q && best.as_ref().is_none_or(|(bq, _, _)| q > *bq) {
            best = Some((q, old.clone(), new));
        }
        // Grow: a sandwiched neighbour face across a ring edge.
        let mut grown = false;
        for &(p, r) in &ring {
            if let Some((f, ta, tb)) = sandwiched_across(cx, a, b, p, r, &faces)? {
                faces.push(f);
                old.push(ta);
                old.push(tb);
                grown = true;
                break;
            }
        }
        if !grown {
            break;
        }
    }
    let Some((_, old, new)) = best else { return Ok(false) };
    let change = Change { old, new, new_boundary: Vec::new(), drop_curves: Vec::new() };
    Ok(cx.replace(&change)?.is_some())
}

/// A face `(p, r, x)` sandwiched between `a` and `b` (tetrahedra `a ∪ f` and `b ∪ f` both exist)
/// next to edge `pr`, not already in `faces`.
fn sandwiched_across(cx: &mut TetComplex, a: u32, b: u32, p: u32, r: u32, faces: &[[u32; 3]]) -> Result<Option<([u32; 3], u32, u32)>> {
    let Some(ring) = cx.edge_ring(p, r)? else { return Ok(None) };
    // Tetrahedra around pr containing a, and containing b.
    for &ta in &ring.tets {
        let va = cx.tet(ta)?;
        if !va.contains(&a) {
            continue;
        }
        let x = va.iter().copied().find(|&v| v != a && v != p && v != r).unwrap_or(NONE);
        if x == NONE || x == b {
            continue;
        }
        let f = [p, r, x];
        if faces.iter().any(|g| sorted3(*g) == sorted3(f)) {
            continue;
        }
        let k = va.iter().position(|&v| v == a).unwrap_or(0);
        let tb = cx.adj(ta)?.get(k).copied().unwrap_or(NONE);
        if tb != NONE && cx.tet(tb)?.contains(&b) {
            return Ok(Some((f, ta, tb)));
        }
    }
    Ok(None)
}

/// How an inserted vertex is placed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Site {
    /// Inside the domain, in (or on) tetrahedron `t`.
    Interior(u32),
    /// On boundary face `k` of tetrahedron `t` (the vertex lies on that face's B-rep face).
    Face(u32, usize),
    /// On the boundary edge `ab` (a feature edge, which is split, or an edge inside one face).
    Edge(u32, u32),
}

/// Insert a vertex at `p` with the Delaunay cavity: the tetrahedra whose circumsphere contains
/// `p`, grown from the site's tetrahedra across interior faces, then shrunk until every face of
/// the cavity sees `p` (star-shaped). A boundary site also removes the boundary faces of the
/// cavity on the same B-rep face(s) whose circumcircle holds `p` and that connect to the site
/// without crossing a feature edge, and fans the new vertex to their outline. Returns the new
/// vertex, or `None` when no valid cavity exists (nothing changes then).
pub fn insert_delaunay(cx: &mut TetComplex, p: Point3, kind: VertexKind, site: Site, max_cavity: usize) -> Result<Option<u32>> {
    insert_delaunay_apart(cx, p, kind, site, max_cavity, 0.0)
}

/// [`insert_delaunay`], refused when a vertex of the cavity lies closer than `min_gap` to `p`
/// (the nearest vertex of a Delaunay mesh is a vertex of the cavity), so refinement cannot
/// pile vertices onto one spot.
pub fn insert_delaunay_apart(cx: &mut TetComplex, p: Point3, kind: VertexKind, site: Site, max_cavity: usize, min_gap: f64) -> Result<Option<u32>> {
    let (roots, root_faces) = site_roots(cx, site)?;
    if roots.is_empty() {
        return Ok(None);
    }
    // Grow the conflict region.
    let mut cavity: Vec<u32> = roots.clone();
    let mut inset: IdMap<u32, ()> = cavity.iter().map(|&t| (t, ())).collect();
    let mut i = 0;
    while let Some(&t) = cavity.get(i) {
        i += 1;
        if cavity.len() > max_cavity {
            return Ok(None);
        }
        for n in cx.adj(t)? {
            if n == NONE || inset.contains_key(&n) {
                continue;
            }
            if cx.in_sphere(cx.tet(n)?, &p)? {
                inset.insert(n, ());
                cavity.push(n);
            }
        }
    }
    if min_gap > 0.0 {
        let px = p.coords();
        for &t in &cavity {
            for v in cx.tet(t)? {
                if measure::dist(cx.xyz(v)?, px) < min_gap {
                    return Ok(None);
                }
            }
        }
    }
    let curve_split = match site {
        Site::Edge(a, b) => cx.curve(a, b).map(|c| (a, b, c)),
        _ => None,
    };
    // Shrink until star-shaped; boundary faces of the cavity on the site's face(s) are removed
    // when connected to the site and in conflict.
    for _ in 0..max_cavity.max(8) {
        let plan = plan_cavity(cx, &cavity, &p, site, &root_faces)?;
        match plan {
            Plan::Ready { new, new_boundary } => {
                let len = cx.journal.steps.len();
                let v = cx.add_point(p, kind)?;
                let new: Vec<[u32; 4]> = new.into_iter().map(|f| [f[0], f[1], f[2], v]).collect();
                let new_boundary: Vec<([u32; 3], u32)> = new_boundary.into_iter().map(|(e, tag)| ([e.0, e.1, v], tag)).collect();
                let change =
                    Change { old: cavity.clone(), new, new_boundary, drop_curves: curve_split.map(|(a, b, _)| (a, b)).into_iter().collect() };
                if cx.replace(&change)?.is_some() {
                    if let Some((a, b, c)) = curve_split {
                        cx.set_curve(a, b, None);
                        cx.set_curve(a, v, Some(c));
                        cx.set_curve(v, b, Some(c));
                    }
                    return Ok(Some(v));
                }
                // Undo the vertex.
                if cx.journal.recording {
                    cx.rollback(len)?;
                } else {
                    cx.drop_last_point(v)?;
                }
                return Ok(None);
            }
            Plan::Drop(t) => {
                if roots.contains(&t) {
                    return Ok(None);
                }
                cavity.retain(|&c| c != t);
                // Keep the cavity connected to the roots.
                cavity = connected(cx, &cavity, &roots)?;
            }
            Plan::Fail => return Ok(None),
        }
    }
    Ok(None)
}

/// The tetrahedra a site must remove, and its boundary faces (as sorted keys) it must remove.
fn site_roots(cx: &mut TetComplex, site: Site) -> Result<(Vec<u32>, Vec<[u32; 3]>)> {
    Ok(match site {
        Site::Interior(t) => (vec![t], Vec::new()),
        Site::Face(t, k) => {
            let vs = cx.tet(t)?;
            if cx.adj(t)?.get(k) != Some(&NONE) {
                return Ok((Vec::new(), Vec::new()));
            }
            (vec![t], vec![sorted3(face_toward(vs, k))])
        }
        Site::Edge(a, b) => {
            let Some(ring) = cx.edge_ring(a, b)? else { return Ok((Vec::new(), Vec::new())) };
            if ring.closed {
                return Ok((Vec::new(), Vec::new()));
            }
            let mut faces = Vec::new();
            if let (Some(&v0), Some(&vl)) = (ring.verts.first(), ring.verts.last()) {
                faces.push(sorted3([a, b, v0]));
                faces.push(sorted3([a, b, vl]));
            }
            (ring.tets, faces)
        }
    })
}

enum Plan {
    /// New tetrahedra as faces (to be joined to the new vertex) and new boundary edges.
    Ready {
        new: Vec<[u32; 3]>,
        new_boundary: Vec<((u32, u32), u32)>,
    },
    /// Drop this tetrahedron from the cavity and try again.
    Drop(u32),
    Fail,
}

/// Plan the fill of a cavity from `p`.
fn plan_cavity(cx: &TetComplex, cavity: &[u32], p: &Point3, site: Site, root_faces: &[[u32; 3]]) -> Result<Plan> {
    let inset: IdMap<u32, ()> = cavity.iter().map(|&t| (t, ())).collect();
    // Boundary faces of cavity tetrahedra: candidates for removal (only for boundary sites).
    let mut bfaces: Vec<([u32; 3], [u32; 3], u32, u32)> = Vec::new(); // (oriented, sorted, tag, tet)
    let mut outer: Vec<([u32; 3], u32)> = Vec::new(); // (oriented face, owning tet)
    for &t in cavity {
        let vs = cx.tet(t)?;
        let ns = cx.adj(t)?;
        let tags = cx.tags(t)?;
        for k in 0..4 {
            let n = ns[k];
            if n != NONE && inset.contains_key(&n) {
                continue;
            }
            let f = face_toward(vs, k);
            if n == NONE {
                bfaces.push((f, sorted3(f), tags[k], t));
            } else {
                outer.push((f, t));
            }
        }
    }
    // Removed boundary faces D: for a boundary site, the root faces plus faces connected to them
    // across non-feature edges, on the same tag, whose circumcircle contains p (approximately:
    // the diametral sphere of the triangle).
    let mut removed: Vec<usize> = Vec::new();
    if !matches!(site, Site::Interior(_)) {
        let mut queue: Vec<usize> = bfaces.iter().enumerate().filter(|(_, f)| root_faces.contains(&f.1)).map(|(i, _)| i).collect();
        if queue.len() != root_faces.len() {
            return Ok(Plan::Fail);
        }
        let px = p.coords();
        let mut taken = vec![false; bfaces.len()];
        for &i in &queue {
            if let Some(s) = taken.get_mut(i) {
                *s = true;
            }
        }
        let mut qi = 0;
        while let Some(&i) = queue.get(qi) {
            qi += 1;
            let Some(&(f, _, tag, _)) = bfaces.get(i) else { continue };
            for e in 0..3 {
                let (u, w) = (f[e], f[(e + 1) % 3]);
                if cx.curve(u, w).is_some() {
                    continue;
                }
                for (j, g) in bfaces.iter().enumerate() {
                    if taken.get(j).copied().unwrap_or(true) || g.2 != tag || !(g.0.contains(&u) && g.0.contains(&w)) {
                        continue;
                    }
                    let [ga, gb, gc] = [cx.xyz(g.0[0])?, cx.xyz(g.0[1])?, cx.xyz(g.0[2])?];
                    let Some(o) = measure::triangle_circumcenter(ga, gb, gc) else { continue };
                    let r = measure::dist(o, ga);
                    if measure::dist(o, px) < r {
                        if let Some(s) = taken.get_mut(j) {
                            *s = true;
                        }
                        queue.push(j);
                    }
                }
            }
        }
        removed = queue;
        // Vertices strictly inside D would be deleted: only allowed for face vertices of the
        // same tag that no kept face uses.
        let mut in_kept: IdMap<u32, ()> = IdMap::default();
        for (j, g) in bfaces.iter().enumerate() {
            if !removed.contains(&j) {
                for v in g.0 {
                    in_kept.insert(v, ());
                }
            }
        }
        for &(f, _) in &outer {
            for v in f {
                in_kept.insert(v, ());
            }
        }
        for &i in &removed {
            if let Some(g) = bfaces.get(i) {
                for v in g.0 {
                    if !in_kept.contains_key(&v) && !matches!(cx.kind(v), VertexKind::Face(_) | VertexKind::Interior) {
                        return Ok(Plan::Fail);
                    }
                }
            }
        }
    }
    // Every kept face must see p strictly.
    let mut new = Vec::with_capacity(outer.len() + bfaces.len());
    for &(f, t) in &outer {
        if cx.orient_with(f, p)? != Sign::Positive {
            return Ok(Plan::Drop(t));
        }
        new.push(f);
    }
    for (j, &(f, _, _, t)) in bfaces.iter().enumerate() {
        if removed.contains(&j) {
            continue;
        }
        if cx.orient_with(f, p)? != Sign::Positive {
            return Ok(Plan::Drop(t));
        }
        new.push(f);
    }
    // The outline of D (edges in exactly one removed face) becomes new boundary, fanned to p.
    let mut new_boundary = Vec::new();
    if !removed.is_empty() {
        let mut count: IdMap<(u32, u32), ((u32, u32), u32, usize)> = IdMap::default();
        for &i in &removed {
            let Some(&(f, _, tag, _)) = bfaces.get(i) else { continue };
            for e in 0..3 {
                let (u, w) = (f[e], f[(e + 1) % 3]);
                let c = count.entry(sorted2(u, w)).or_insert(((u, w), tag, 0));
                c.2 += 1;
            }
        }
        let split = match site {
            Site::Edge(a, b) => Some(sorted2(a, b)),
            _ => None,
        };
        // D must be a disk, its outline one loop. Grown face by face, it can close around a
        // kept boundary face whose circumcircle misses p: the fan over that inner loop and the
        // kept face would enclose tetrahedra joined to the mesh by vertices alone. D is
        // connected by construction; with Euler characteristic 1 and every outline vertex on two
        // outline edges, it is a disk.
        let mut vertices: Vec<u32> = removed.iter().filter_map(|&i| bfaces.get(i)).flat_map(|g| g.0).collect();
        vertices.sort_unstable();
        vertices.dedup();
        let euler = vertices.len() as i64 - count.len() as i64 + removed.len() as i64;
        let mut degree: IdMap<u32, usize> = IdMap::default();
        for c in count.values().filter(|c| c.2 == 1) {
            *degree.entry(c.0.0).or_default() += 1;
            *degree.entry(c.0.1).or_default() += 1;
        }
        if euler != 1 || degree.values().any(|&d| d != 2) {
            return Ok(Plan::Fail);
        }
        let mut edges: Vec<((u32, u32), u32)> = count.into_iter().filter(|(k, c)| c.2 == 1 && Some(*k) != split).map(|(_, c)| (c.0, c.1)).collect();
        edges.sort_unstable();
        new_boundary = edges;
    }
    // A cavity whose tetrahedra all go must still leave something to fill.
    if new.is_empty() {
        return Ok(Plan::Fail);
    }
    Ok(Plan::Ready { new, new_boundary })
}

/// The tetrahedra of `cavity` connected to the roots through faces inside the cavity.
fn connected(cx: &TetComplex, cavity: &[u32], roots: &[u32]) -> Result<Vec<u32>> {
    let inset: IdMap<u32, ()> = cavity.iter().map(|&t| (t, ())).collect();
    let mut seen: IdMap<u32, ()> = IdMap::default();
    let mut out: Vec<u32> = Vec::new();
    for &r in roots {
        if inset.contains_key(&r) && seen.insert(r, ()).is_none() {
            out.push(r);
        }
    }
    let mut i = 0;
    while let Some(&t) = out.get(i) {
        i += 1;
        for n in cx.adj(t)? {
            if n != NONE && inset.contains_key(&n) && seen.insert(n, ()).is_none() {
                out.push(n);
            }
        }
    }
    Ok(out)
}

/// The edges of a set of tetrahedra, sorted and unique.
pub fn edges_of(cx: &TetComplex, tets: &[u32]) -> Result<Vec<(u32, u32)>> {
    let mut out = Vec::with_capacity(tets.len() * 6);
    for &t in tets {
        let vs = cx.tet(t)?;
        out.extend(EDGES.iter().map(|&(i, j)| sorted2(vs[i], vs[j])));
    }
    out.sort_unstable();
    out.dedup();
    Ok(out)
}

/// Collapse vertex `a` onto its neighbour `b` (edge contraction): every tetrahedron around `a`
/// takes `b` in its place and those holding the edge disappear. Allowed when `a` may move to `b`
/// without leaving its boundary entity — an interior vertex anywhere, a face vertex along an edge
/// inside its B-rep face, a curve vertex along its curve, never a corner — when the boundary
/// faces it changes keep their orientation within 30°, every new edge is at most `max_len`, and
/// the worst new tetrahedron is at least `min(worst old, floor)`. Returns whether it collapsed.
///
/// Edge contraction is the coarsening operation of mesh adaptation by local modification (X. Li,
/// M. S. Shephard, M. W. Beall, *3D anisotropic mesh adaptation by mesh modification*, CMAME
/// 194, 2005); here it removes the boundary recovery's closely spaced Steiner points and dense
/// tessellation vertices where the size field asks for larger elements.
pub fn collapse(cx: &mut TetComplex, a: u32, b: u32, max_len: f64, floor: f64) -> Result<bool> {
    let ka = cx.kind(a);
    let allowed = match ka {
        VertexKind::Corner => false,
        VertexKind::Interior => true,
        VertexKind::Curve(c) => cx.curve(a, b) == Some(c),
        VertexKind::Face(tag) => {
            if cx.curve(a, b).is_some() {
                false
            } else {
                match cx.kind(b) {
                    VertexKind::Interior => false,
                    _ => {
                        // The edge must be a boundary edge inside face `tag`.
                        let Some(ring) = cx.edge_ring(a, b)? else { return Ok(false) };
                        if ring.closed {
                            false
                        } else {
                            let ends = [ring.tets.first().copied(), ring.tets.last().copied()];
                            let (v0, vl) = (ring.verts.first().copied().unwrap_or(NONE), ring.verts.last().copied().unwrap_or(NONE));
                            boundary_tag(cx, ends[0].unwrap_or(NONE), [a, b, v0])? == Some(tag)
                                && boundary_tag(cx, ends[1].unwrap_or(NONE), [a, b, vl])? == Some(tag)
                        }
                    }
                }
            }
        }
    };
    if !allowed || cx.edge_tet(a, b)?.is_none() {
        return Ok(false);
    }
    let star = cx.star(a)?;
    let before = worst(cx, &star)?;
    let xb = cx.xyz(b)?;
    let mut new = Vec::with_capacity(star.len());
    let mut new_boundary = Vec::new();
    for &t in &star {
        let vs = cx.tet(t)?;
        if vs.contains(&b) {
            continue;
        }
        let moved = vs.map(|v| if v == a { b } else { v });
        for v in moved {
            if v != b && measure::dist(cx.xyz(v)?, xb) > max_len {
                return Ok(false);
            }
        }
        let (ns, tags) = (cx.adj(t)?, cx.tags(t)?);
        for k in 0..4 {
            if ns[k] != NONE || vs[k] == a {
                continue;
            }
            let f = face_toward(vs, k);
            let g = f.map(|v| if v == a { b } else { v });
            let n0 = measure::triangle_normal(cx.xyz(f[0])?, cx.xyz(f[1])?, cx.xyz(f[2])?);
            let n1 = measure::triangle_normal(cx.xyz(g[0])?, cx.xyz(g[1])?, cx.xyz(g[2])?);
            if !(measure::dot(n0, n1) >= COLLAPSE_NORMAL_COS) {
                return Ok(false);
            }
            new_boundary.push((g, tags[k]));
        }
        new.push(moved);
    }
    if new.is_empty() {
        return Ok(false);
    }
    let mut after = f64::INFINITY;
    for vs in &new {
        after = after.min(tet_quality(cx, *vs)?);
    }
    if !(after > 0.0) || after < before.min(floor) {
        return Ok(false);
    }
    // Curve edges at a move to b.
    let mut drop = Vec::new();
    let mut add = Vec::new();
    if let VertexKind::Curve(c) = ka {
        for u in cx.neighbours(a)? {
            if cx.curve(a, u) == Some(c) {
                drop.push((a, u));
                if u != b {
                    add.push((b, u, c));
                }
            }
        }
    }
    let change = Change { old: star, new, new_boundary, drop_curves: drop.clone() };
    if cx.replace(&change)?.is_none() {
        return Ok(false);
    }
    for (x, y) in drop {
        cx.set_curve(x, y, None);
    }
    for (x, y, c) in add {
        cx.set_curve(x, y, Some(c));
    }
    Ok(true)
}

/// cos 30°: the largest turn of a boundary face's normal a collapse may cause.
const COLLAPSE_NORMAL_COS: f64 = 0.866_025_403_784_438_6;
