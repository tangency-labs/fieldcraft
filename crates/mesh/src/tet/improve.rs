//! Aggressive mesh improvement: the schedule of Klingner and Shewchuk (IMR 2007, §4.2,
//! Figure 4) over smoothing, topological transformations and vertex insertion, monotone in the
//! worst quality (the minimum dihedral sine).
//!
//! Each pass smooths; if that is not enough progress, a topological pass (edge removal,
//! multi-face removal, and deletion of flat tetrahedra lying on one B-rep face) follows; if still
//! not enough, an insertion pass tries new vertices in the worst tetrahedra (a boundary face's
//! barycentre first, then the tetrahedron's, then the midpoints of its boundary edges). An
//! insertion is a composite operation: the vertex goes in with a Delaunay cavity, is smoothed,
//! and topological passes and smoothing run around it; if the worst of the tetrahedra it made
//! or changed is not better than the worst of those it removed, everything is rolled back
//! ([`super::complex::Journal`]). Progress is judged on the worst quality and seven thresholded
//! means (§4.2); the schedule stops after three passes without progress, or at its pass and
//! insertion caps (reported).
//!
//! The paper chooses the inserted vertex's cavity by an optimal graph cut (§3.4); this
//! implementation uses the Delaunay cavity followed by the same smoothing and topological
//! repair, which the paper's own composite operation relies on to make insertion pay off.

use crate::Result;
use crate::predicates::Point3;

use super::complex::{EDGES, IdMap, NONE, Step, TetComplex, VertexKind, Walk, face_toward, sorted4};
use super::measure::{self, P};
use super::ops::{self, Site, remove_edge_at, remove_face, tet_quality};
use super::refine::{Context, MIN_GAP};
use super::smooth::smooth_vertex;

/// Settings of the schedule.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImproveOptions {
    pub max_passes: usize,
    pub max_insertions: usize,
    /// Tetrahedra at or above this quality are left alone by the topological and smoothing
    /// passes (default sin 40°, the paper's desperation threshold: a mesh whose tetrahedra all
    /// have dihedral angles between 40° and 140° needs nothing more).
    pub good: f64,
    /// Insertion runs every pass while the worst quality is below this (sin 30°: dihedral
    /// angles between 30° and 150°).
    pub target: f64,
}

impl Default for ImproveOptions {
    fn default() -> Self {
        ImproveOptions { max_passes: 10, max_insertions: 50_000, good: SIN_40, target: 0.5 }
    }
}

/// What the schedule did.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ImproveStats {
    pub passes: usize,
    pub smoothed: usize,
    pub topological: usize,
    pub inserted: usize,
    /// The worst quality (smallest dihedral sine) after each pass.
    pub worst: Vec<f64>,
    pub stopped: Option<String>,
}

/// sin 1°, 5°, 10°, 15°, 25°, 35°, 45° (the thresholds of §4.2), as literals so they are the same
/// everywhere.
pub const SIN_DEG: [f64; 7] = [
    0.017_452_406_437_283_51,
    0.087_155_742_747_658_17,
    0.173_648_177_666_930_33,
    0.258_819_045_102_520_74,
    0.422_618_261_740_699_44,
    0.573_576_436_351_046,
    std::f64::consts::FRAC_1_SQRT_2,
];

/// Most tetrahedra one insertion pass tries (the worst first).
pub const MAX_INSERTION_LIST: usize = 2000;

/// sin 40°.
pub const SIN_40: f64 = 0.642_787_609_686_539_3;

/// The worst quality and the seven thresholded means.
fn indicators(cx: &TetComplex) -> Result<[f64; 8]> {
    let mut out = [0.0; 8];
    out[0] = f64::INFINITY;
    let mut n = 0usize;
    for (_, vs) in cx.live_tets() {
        let q = tet_quality(cx, vs)?;
        out[0] = out[0].min(q);
        for (i, s) in SIN_DEG.iter().enumerate() {
            out[i + 1] += q.min(*s);
        }
        n += 1;
    }
    for m in out.iter_mut().skip(1) {
        *m /= n.max(1) as f64;
    }
    Ok(out)
}

fn progressed(now: &[f64; 8], before: &[f64; 8]) -> bool {
    now[0] > before[0] || now.iter().zip(before).skip(1).any(|(a, b)| *a >= b + 1e-4)
}

/// Smooth every movable vertex of a tetrahedron below `good`. Returns the moves.
fn smoothing_pass(cx: &mut TetComplex, ctx: &Context, good: f64) -> Result<usize> {
    let mut want = vec![false; cx.point_count()];
    for (_, vs) in cx.live_tets() {
        if tet_quality(cx, vs)? < good {
            for v in vs {
                if let Some(w) = want.get_mut(v as usize) {
                    *w = true;
                }
            }
        }
    }
    let mut moved = 0;
    for (v, w) in want.iter().enumerate() {
        if *w && smooth_vertex(cx, ctx, v as u32)? {
            moved += 1;
        }
    }
    Ok(moved)
}

/// Delete a tetrahedron capping two boundary faces of one B-rep face that fold by at most 20°
/// ([`ops::BOUNDARY_FLIP_COS`]): its dihedral angle at their shared edge is at least 160°, and
/// its two other faces cover the same piece of surface.
fn delete_flat(cx: &mut TetComplex, t: u32) -> Result<bool> {
    let vs = cx.tet(t)?;
    let (ns, tags) = (cx.adj(t)?, cx.tags(t)?);
    let bnd: Vec<usize> = (0..4).filter(|&k| ns[k] == NONE).collect();
    let &[i, j] = bnd.as_slice() else { return Ok(false) };
    if tags[i] != tags[j] {
        return Ok(false);
    }
    let normal = |k: usize| -> Result<P> {
        let f = face_toward(vs, k);
        Ok(measure::triangle_normal(cx.xyz(f[0])?, cx.xyz(f[1])?, cx.xyz(f[2])?))
    };
    // The inward normals of the two faces of such a cap point nearly the same way.
    let (a, b) = (normal(i)?, normal(j)?);
    if !(measure::dot(a, b) >= ops::BOUNDARY_FLIP_COS) {
        return Ok(false);
    }
    cx.delete_boundary_tet(t)
}

/// A flat tetrahedron spanning two boundary faces of different B-rep faces that meet at a
/// nearly flat crease (two faces tangent along their curve): its two apexes lie almost on the
/// crease, so no flip or move fixes it. Collapse an apex face vertex onto an end of the crease
/// edge ([`ops::collapse`]), which takes the tetrahedron away.
fn collapse_into_crease(cx: &mut TetComplex, t: u32) -> Result<bool> {
    let vs = cx.tet(t)?;
    if tet_quality(cx, vs)? > 0.1 {
        return Ok(false);
    }
    let (ns, tags) = (cx.adj(t)?, cx.tags(t)?);
    let bnd: Vec<usize> = (0..4).filter(|&k| ns[k] == NONE).collect();
    let &[i, j] = bnd.as_slice() else { return Ok(false) };
    if tags[i] == tags[j] {
        return Ok(false);
    }
    // The faces opposite i and j share the edge of the two other vertices; vs[i] and vs[j] are
    // the apexes (vs[i] lies on face j, vs[j] on face i).
    let edge: Vec<u32> = (0..4).filter(|&k| k != i && k != j).map(|k| vs[k]).collect();
    let (ex, longest) = {
        let x = cx.coords(vs)?;
        (edge.clone(), measure::edge_range(&x).1)
    };
    for apex in [vs[i], vs[j]] {
        if !matches!(cx.kind(apex), VertexKind::Face(_)) {
            continue;
        }
        for &e in &ex {
            if ops::collapse(cx, apex, e, 2.0 * longest, 0.0)? {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// Try to remove each tetrahedron of `list`: its edges, then its faces. Returns the changes.
pub fn topological_pass(cx: &mut TetComplex, list: &[u32]) -> Result<usize> {
    let mut changes = 0;
    for &t in list {
        if !cx.is_live(t) {
            continue;
        }
        let vs = cx.tet(t)?;
        if delete_flat(cx, t)? || collapse_into_crease(cx, t)? {
            changes += 1;
            continue;
        }
        let mut gone = false;
        for &(i, j) in &EDGES {
            if remove_edge_at(cx, t, vs[i], vs[j])? {
                gone = true;
                break;
            }
        }
        if gone {
            changes += 1;
            continue;
        }
        for k in 0..4 {
            if remove_face(cx, t, k)? {
                changes += 1;
                break;
            }
        }
    }
    Ok(changes)
}

fn below(cx: &TetComplex, good: f64) -> Result<Vec<u32>> {
    let mut out = Vec::new();
    for (t, vs) in cx.live_tets() {
        if tet_quality(cx, vs)? < good {
            out.push(t);
        }
    }
    Ok(out)
}

/// The quality changes since journal step `len`: the worst quality of the tetrahedra removed
/// (at their coordinates before) and of those made or moved (now).
fn change_quality(cx: &mut TetComplex, len: usize) -> Result<(f64, f64)> {
    let steps: Vec<Step> = cx.journal.steps.get(len..).unwrap_or_default().to_vec();
    // Earliest position of each moved vertex.
    let mut first: IdMap<u32, Point3> = IdMap::default();
    for s in &steps {
        if let Step::Moved { v, from } = s {
            first.entry(*v).or_insert(*from);
        }
    }
    let old_xyz = |cx: &TetComplex, v: u32| -> Result<P> { first.get(&v).map_or_else(|| cx.xyz(v), |p| Ok(p.coords())) };
    // Net removed and added tetrahedra by sorted vertices.
    let mut net: IdMap<[u32; 4], (i32, [u32; 4])> = IdMap::default();
    let key = sorted4;
    for s in &steps {
        match s {
            Step::Replaced { removed, added, .. } => {
                for r in removed {
                    net.entry(key(r.vs)).or_insert((0, r.vs)).0 -= 1;
                }
                for &vs in added {
                    net.entry(key(vs)).or_insert((0, vs)).0 += 1;
                }
            }
            Step::Deleted { removed, .. } => net.entry(key(removed.vs)).or_insert((0, removed.vs)).0 -= 1,
            _ => {}
        }
    }
    let (mut before, mut after) = (f64::INFINITY, f64::INFINITY);
    let mut keys: Vec<[u32; 4]> = net.keys().copied().collect();
    keys.sort_unstable();
    for k in &keys {
        let Some(&(c, vs)) = net.get(k) else { continue };
        if c < 0 {
            let x = [old_xyz(cx, vs[0])?, old_xyz(cx, vs[1])?, old_xyz(cx, vs[2])?, old_xyz(cx, vs[3])?];
            before = before.min(measure::quality(&x));
        }
    }
    // The tetrahedra made, and those around moved vertices.
    for k in &keys {
        if let Some(&(c, vs)) = net.get(k)
            && c > 0
        {
            after = after.min(tet_quality(cx, vs)?);
        }
    }
    let mut moved: Vec<u32> = first.keys().copied().collect();
    moved.sort_unstable();
    for v in moved {
        if !cx.is_used(v) {
            continue;
        }
        for t in cx.star(v)? {
            let vs = cx.tet(t)?;
            if net.get(&key(vs)).is_some_and(|e| e.0 > 0) {
                continue;
            }
            after = after.min(tet_quality(cx, vs)?);
            let x = [old_xyz(cx, vs[0])?, old_xyz(cx, vs[1])?, old_xyz(cx, vs[2])?, old_xyz(cx, vs[3])?];
            before = before.min(measure::quality(&x));
        }
    }
    Ok((before, after))
}

/// The composite insertion of §3.5 (see the module docs). Returns whether it was kept.
fn attempt_insert(cx: &mut TetComplex, ctx: &Context, p: Point3, kind: VertexKind, site: Site) -> Result<bool> {
    let was = cx.journal.recording;
    cx.journal.recording = true;
    let len = cx.journal.steps.len();
    let kept = (|| -> Result<bool> {
        // Refused within MIN_GAP of the local size of an existing vertex, as in refinement: a
        // tetrahedron that stays below the target would otherwise draw an insertion every pass
        // and pile vertices around one of its own, whose star then grows without bound (#311).
        let gap = MIN_GAP * ctx.h(p.coords(), 0.0);
        let Some(v) = ops::insert_delaunay_apart(cx, p, kind, site, gap)? else { return Ok(false) };
        for _ in 0..4 {
            smooth_vertex(cx, ctx, v)?;
            if !cx.is_used(v) {
                break;
            }
            let star = cx.star(v)?;
            if topological_pass(cx, &star)? == 0 {
                break;
            }
        }
        if cx.is_used(v) {
            let mut around: Vec<u32> = cx.neighbours(v)?;
            around.push(v);
            for u in around {
                smooth_vertex(cx, ctx, u)?;
            }
        }
        let (before, after) = change_quality(cx, len)?;
        Ok(after > before && after.is_finite())
    })();
    let keep = match kept {
        Ok(k) => k,
        Err(e) => {
            cx.rollback(len)?;
            cx.journal.recording = was;
            return Err(e);
        }
    };
    if !keep {
        cx.rollback(len)?;
    }
    if !was {
        cx.journal.steps.truncate(len);
    }
    cx.journal.recording = was;
    Ok(keep)
}

/// Try inserting a vertex in each tetrahedron of `list` (§4.2). Returns the insertions kept.
fn insertion_pass(cx: &mut TetComplex, ctx: &Context, list: &[u32], budget: usize, good: f64) -> Result<usize> {
    let mut kept = 0;
    for &t in list {
        if kept >= budget {
            break;
        }
        if !cx.is_live(t) {
            continue;
        }
        let vs = cx.tet(t)?;
        if tet_quality(cx, vs)? >= good {
            continue;
        }
        let ns = cx.adj(t)?;
        let tags = cx.tags(t)?;
        let mut done = false;
        // A tetrahedron lying on two boundary faces that meet at a nearly flat edge (two faces
        // that are tangent along their curve) is flat however its vertices move: the edge needs
        // a second tetrahedron. A vertex inside, off the edge's midpoint along the inward
        // normal, splits the edge's dihedral angle.
        let bnd: Vec<usize> = (0..4).filter(|&k| ns[k] == NONE).collect();
        if let &[i, j] = bnd.as_slice() {
            let edge: Vec<u32> = (0..4).filter(|&k| k != i && k != j).map(|k| vs[k]).collect();
            if let &[a, b] = edge.as_slice() {
                let n = |k: usize| -> Result<P> {
                    let f = face_toward(vs, k);
                    Ok(measure::triangle_normal(cx.xyz(f[0])?, cx.xyz(f[1])?, cx.xyz(f[2])?))
                };
                let inward = measure::add(n(i)?, n(j)?);
                let len = measure::norm(inward);
                let (pa, pb) = (cx.xyz(a)?, cx.xyz(b)?);
                let l = measure::dist(pa, pb);
                // Unit inward normals within 25° of each other: |n_i + n_j| > 2 cos 12.5°.
                if len > 1.952 {
                    let target = measure::add(measure::mid(pa, pb), measure::scale(inward, 0.4 * l / len));
                    if let Ok(p) = Point3::try_from(target)
                        && let Some(Walk::Inside(at)) = cx.walk(t, &p, 1 << 12)?
                        && attempt_insert(cx, ctx, p, VertexKind::Interior, Site::Interior(at))?
                    {
                        done = true;
                    }
                }
            }
        }
        // Boundary faces' barycentres.
        for k in 0..4 {
            if done || ns[k] != NONE || tags[k] == NONE {
                continue;
            }
            let f = face_toward(vs, k);
            let g = measure::scale(measure::add(measure::add(cx.xyz(f[0])?, cx.xyz(f[1])?), cx.xyz(f[2])?), 1.0 / 3.0);
            let kind = VertexKind::Face(tags[k]);
            if let Some(p) = ctx.place(kind, g)
                && cx.is_live(t)
                && attempt_insert(cx, ctx, p, kind, Site::Face(t, k))?
            {
                done = true;
            }
        }
        // The barycentre.
        if !done && cx.is_live(t) && cx.tet(t)? == vs {
            let g = measure::centroid(&cx.coords(vs)?);
            if let Ok(p) = Point3::try_from(g)
                && attempt_insert(cx, ctx, p, VertexKind::Interior, Site::Interior(t))?
            {
                done = true;
            }
        }
        // Boundary edges' midpoints.
        if !done && cx.is_live(t) && cx.tet(t)? == vs {
            for &(i, j) in &EDGES {
                let (a, b) = (vs[i], vs[j]);
                let Some(ring) = cx.edge_ring(a, b)? else { continue };
                if ring.closed {
                    continue;
                }
                let tag = ring
                    .tets
                    .first()
                    .and_then(|&r| {
                        let rv = cx.tet(r).ok()?;
                        let rt = cx.tags(r).ok()?;
                        let rn = cx.adj(r).ok()?;
                        (0..4).find(|&k| rn[k] == NONE && face_toward(rv, k).contains(&a) && face_toward(rv, k).contains(&b)).map(|k| rt[k])
                    })
                    .unwrap_or(NONE);
                let kind = match cx.curve(a, b) {
                    Some(c) => VertexKind::Curve(c),
                    None if tag != NONE => VertexKind::Face(tag),
                    None => continue,
                };
                let m = measure::mid(cx.xyz(a)?, cx.xyz(b)?);
                if let Some(p) = ctx.place(kind, m)
                    && attempt_insert(cx, ctx, p, kind, Site::Edge(a, b))?
                {
                    done = true;
                    break;
                }
                if !cx.is_live(t) {
                    break;
                }
            }
        }
        kept += usize::from(done);
    }
    Ok(kept)
}

/// The improvement schedule (see the module docs). Each pass smooths and runs a topological
/// pass over the tetrahedra below the threshold, then an insertion pass when the worst
/// tetrahedron did not improve; it stops after three passes without progress (§4.2).
pub fn improve(cx: &mut TetComplex, ctx: &Context, o: &ImproveOptions) -> Result<ImproveStats> {
    let mut stats = ImproveStats::default();
    stats.smoothed += smoothing_pass(cx, ctx, o.good)?;
    let list = below(cx, o.good)?;
    stats.topological += topological_pass(cx, &list)?;
    let mut failed = 0;
    // A pass ends where the next begins: its indicators are the next one's starting point.
    let mut carried: Option<[f64; 8]> = None;
    while failed < 3 {
        if stats.passes >= o.max_passes {
            stats.stopped = Some(format!("improvement stopped after {} passes", o.max_passes));
            break;
        }
        stats.passes += 1;
        fieldcraft_model::progress::stage("Improvement pass", Some(o.max_passes as f64), "");
        fieldcraft_model::progress::count(stats.passes as f64);
        fieldcraft_model::progress::check()?;
        let before = match carried.take() {
            Some(b) => b,
            None => indicators(cx)?,
        };
        stats.smoothed += smoothing_pass(cx, ctx, o.good)?;
        let list = below(cx, o.good)?;
        stats.topological += topological_pass(cx, &list)?;
        // Insertion when smoothing and flips did not make enough progress (§4.2), or left a
        // tetrahedron below the target.
        let mid = indicators(cx)?;
        if (!progressed(&mid, &before) || mid[0] < o.target) && stats.inserted < o.max_insertions {
            // The worst 3.5% of the tetrahedra, or (desperation pass, §4.2) every one with a
            // dihedral angle below 40° or above 140° (sine below sin 40°).
            let mut scored: Vec<(f64, u32)> = Vec::new();
            for (t, vs) in cx.live_tets() {
                scored.push((tet_quality(cx, vs)?, t));
            }
            scored.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
            let list: Vec<u32> = if failed == 1 {
                scored.iter().filter(|s| s.0 < SIN_40).take(MAX_INSERTION_LIST).map(|s| s.1).collect()
            } else {
                let n = (scored.len() * 35).div_ceil(1000).clamp(1, MAX_INSERTION_LIST);
                scored.iter().take(n).map(|s| s.1).collect()
            };
            let budget = o.max_insertions - stats.inserted;
            stats.inserted += insertion_pass(cx, ctx, &list, budget, o.good)?;
        }
        let now = indicators(cx)?;
        carried = Some(now);
        stats.worst.push(now[0]);
        if progressed(&now, &before) {
            failed = 0;
        } else {
            failed += 1;
        }
    }
    if stats.stopped.is_none() && stats.inserted >= o.max_insertions {
        stats.stopped = Some(format!("improvement stopped at {} inserted vertices", o.max_insertions));
    }
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cdt::tests::{surface, uv_sphere};
    use crate::cdt::{CdtOptions, tetrahedralize};
    use crate::size::UniformSize;
    use crate::surface::{BoundarySurface, SurfaceOptions};
    use crate::tet::geometry::Tessellation;
    use crate::tet::{QualityOptions, complex_from_cdt};

    /// An insertion closer than `MIN_GAP` × the local size to an existing vertex is refused,
    /// whatever it would do to the quality (#311): a stubborn tetrahedron must not draw vertex
    /// after vertex onto one of its own.
    #[test]
    fn an_insertion_beside_an_existing_vertex_is_refused() {
        const H: f64 = 100.0;
        let (p, t) = uv_sphere(6, 10);
        let p: Vec<[f64; 3]> = p.iter().map(|c| c.map(|x| 5.0 * x)).collect();
        let s = surface(&p, &t, [0.0; 3]);
        let faces = vec![0u32; s.triangles().len()];
        let s = BoundarySurface::new(s.points(), s.triangles(), &faces, &SurfaceOptions::default()).unwrap();
        let cdt = tetrahedralize(&s, &CdtOptions::default()).unwrap();
        let (mut cx, curve_faces) = complex_from_cdt(&s, &cdt, &[], &QualityOptions::default()).unwrap();
        let (geometry, size) = (Tessellation::new(&s, Vec::new()), UniformSize(H));
        let ctx = Context { size: &size, geometry: &geometry, curve_faces: &curve_faces };
        let (points, live) = (cx.point_count(), cx.live_count());
        let mut tried = 0;
        for (t, vs) in cx.live_tets().collect::<Vec<_>>() {
            if !cx.is_live(t) {
                continue;
            }
            // A point half the gap from each vertex toward the centroid: unguarded, the coarse
            // sphere keeps some of these.
            let g = measure::centroid(&cx.coords(vs).unwrap());
            for v in vs {
                let x = cx.xyz(v).unwrap();
                let d = measure::sub(g, x);
                let q = measure::add(x, measure::scale(d, 0.5 * MIN_GAP * H / measure::norm(d)));
                let q = Point3::try_from(q).unwrap();
                if let Some(Walk::Inside(at)) = cx.walk(t, &q, 1 << 12).unwrap() {
                    tried += 1;
                    assert!(!attempt_insert(&mut cx, &ctx, q, VertexKind::Interior, Site::Interior(at)).unwrap());
                }
            }
        }
        assert!(tried > 0);
        assert_eq!((cx.point_count(), cx.live_count()), (points, live));
        cx.validate().unwrap();
    }
}
