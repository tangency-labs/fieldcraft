//! The mutable tetrahedral mesh of the quality stage: explicit vertices (doubles), tetrahedra in
//! tet4 order with neighbours, the boundary faces with their B-rep face tags, and the feature
//! edges (subsegments) with their curve.
//!
//! Every change goes through [`TetComplex::replace`]: a set of tetrahedra is replaced by another
//! set filling the same region. The replacement is checked before anything is changed — every new
//! tetrahedron positive by the exact `orient3d` and admissible by FieldCraft's validity criterion, every face of the region's outside matched once,
//! every new face either shared by two new tetrahedra, matched outside, or declared a new
//! boundary face — and refused (no change) otherwise. So every operation built on it (flips, edge
//! removal, insertion) either keeps a valid, conforming mesh or does nothing, and every operation
//! can be undone by the inverse replacement ([`Journal`]).
//!
//! Layout: per tetrahedron four vertex ids in tet4 order (`det[b − a, c − a, d − a] > 0`), four
//! neighbour ids (neighbour `k` across the face opposite vertex `k`, [`NONE`] on the boundary) and
//! four face tags (the B-rep face of a boundary face, [`NONE`] inside). Slots of removed
//! tetrahedra are recycled in a fixed order, so the same input gives the same mesh.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

use fieldcraft_model::{normalized_jacobian, tet4_jacobian};

use crate::predicates::{Point3, Sign, insphere, tet_volume_sign};
use crate::{Error, Result};

/// Every tetrahedron a change creates has a normalized Jacobian above this: FieldCraft's validity
/// criterion ([`fieldcraft_model::admissible_jacobian`]) with a margin of a thousand, so the
/// finished mesh passes the gate and the solvers whatever the operations did.
pub const ADMISSIBLE_MARGIN: f64 = 1e3 * fieldcraft_model::MIN_NORMALIZED_JACOBIAN; // unit-factors: a multiple of a normalized Jacobian

/// No tetrahedron, no vertex, no tag.
pub const NONE: u32 = u32::MAX;
/// Marks a recycled tetrahedron slot.
const DEAD: u32 = u32::MAX - 1;

/// A multiplicative hasher for small integer keys. The maps it serves are only looked up or
/// drained in a sorted order, never iterated into output.
#[derive(Default)]
pub struct IdHasher(u64);

impl Hasher for IdHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.write_u64(u64::from(b));
        }
    }

    fn write_u32(&mut self, x: u32) {
        self.write_u64(u64::from(x));
    }

    fn write_u64(&mut self, x: u64) {
        self.0 = (self.0.rotate_left(5) ^ x).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
    }
}

/// A hash map keyed by vertex ids.
pub type IdMap<K, V> = HashMap<K, V, BuildHasherDefault<IdHasher>>;

pub(crate) fn corrupt(what: &str) -> Error {
    Error::Failed(format!("tetrahedral mesh is inconsistent: {what}"))
}

/// The face of `vs` opposite slot `k`, ordered so that vertex `k` is on its positive side (the
/// face's right-hand normal points into the tetrahedron).
pub fn face_toward(vs: [u32; 4], k: usize) -> [u32; 3] {
    let [a, b, c, d] = vs;
    match k {
        0 => [b, d, c],
        1 => [a, c, d],
        2 => [a, d, b],
        _ => [a, b, c],
    }
}

/// A triangle's vertices sorted: the key of an unoriented face.
pub fn sorted3([a, b, c]: [u32; 3]) -> [u32; 3] {
    let mut f = [a, b, c];
    f.sort_unstable();
    f
}

/// A tetrahedron's vertices sorted: the key of an unoriented tetrahedron.
pub fn sorted4(vs: [u32; 4]) -> [u32; 4] {
    let mut k = vs;
    k.sort_unstable();
    k
}

/// An edge's vertices sorted.
pub fn sorted2(a: u32, b: u32) -> (u32, u32) {
    if a < b { (a, b) } else { (b, a) }
}

/// The six edges of a tetrahedron, as slot pairs (the tet10 mid-side node order).
pub use crate::quality::jacobian::TET_EDGES as EDGES;

/// What a vertex may do when the mesh is improved.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VertexKind {
    /// Fixed: a corner of the boundary (where feature curves meet, or an isolated input vertex
    /// that must stay).
    Corner,
    /// On a feature curve (between B-rep faces): moves along it.
    Curve(u32),
    /// On a B-rep face: moves on it.
    Face(u32),
    /// Inside: moves freely.
    Interior,
}

/// A replacement of tetrahedra ([`TetComplex::replace`]).
#[derive(Debug, Clone, Default)]
pub struct Change {
    /// The tetrahedra removed.
    pub old: Vec<u32>,
    /// The tetrahedra added, tet4 order.
    pub new: Vec<[u32; 4]>,
    /// Boundary faces the change creates, oriented into the mesh (as [`face_toward`] gives them
    /// for the new tetrahedron), with their B-rep face tag. A change
    /// that lists some may also remove boundary faces of the old tetrahedra; one that lists none
    /// keeps the boundary exactly.
    pub new_boundary: Vec<([u32; 3], u32)>,
    /// Feature edges the change removes (split or collapsed; the caller registers what replaces
    /// them).
    pub drop_curves: Vec<(u32, u32)>,
}

/// A removed tetrahedron: its vertices and face tags, enough to put it back.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Removed {
    pub vs: [u32; 4],
    pub tags: [u32; 4],
}

/// One undoable step.
#[derive(Debug, Clone)]
pub enum Step {
    /// The tetrahedra removed, and those added (by vertices: slots are recycled).
    Replaced {
        removed: Vec<Removed>,
        added: Vec<[u32; 4]>,
        raw: Box<Raw>,
    },
    Deleted {
        removed: Removed,
        raw: Box<Raw>,
    },
    Moved {
        v: u32,
        from: Point3,
    },
    Added {
        v: u32,
    },
    Curve {
        edge: (u32, u32),
        curve: Option<u32>,
    },
}

/// The data a change overwrites, to put it back exactly: the slots it touches (vertices,
/// neighbours, face tags), the slot count, the free list, the vertices' incidence lists and the
/// live count.
#[derive(Debug, Clone, Default)]
pub struct Raw {
    slots: Vec<(u32, [u32; 4], [u32; 4], [u32; 4])>,
    tets_len: usize,
    free: Vec<u32>,
    vtets: Vec<(u32, Vec<u32>)>,
    live: usize,
}

/// The steps since a checkpoint, newest last.
#[derive(Debug, Clone, Default)]
pub struct Journal {
    pub steps: Vec<Step>,
    pub recording: bool,
}

/// The mutable tetrahedral mesh (see the module docs).
#[derive(Debug, Clone)]
pub struct TetComplex {
    pts: Vec<Point3>,
    kinds: Vec<VertexKind>,
    tets: Vec<[u32; 4]>,
    adj: Vec<[u32; 4]>,
    tags: Vec<[u32; 4]>,
    free: Vec<u32>,
    /// The tetrahedra holding each vertex (complete even where a vertex's tetrahedra touch
    /// only at it, so stars never miss a piece).
    vtets: Vec<Vec<u32>>,
    curves: IdMap<(u32, u32), u32>,
    /// Boundary face tags no change may remove or add a face of ([`Self::freeze_faces`]).
    frozen: Vec<u32>,
    mark: Vec<u32>,
    epoch: u32,
    live: usize,
    max_tets: usize,
    max_points: usize,
    pub journal: Journal,
}

/// Where a walk toward a point ended ([`TetComplex::walk`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Walk {
    /// In this (closed) tetrahedron.
    Inside(u32),
    /// The point is beyond boundary face `k` of this tetrahedron.
    Outside(u32, usize),
}

/// A boundary face given to [`TetComplex::new`]: tetrahedron, the slot of the vertex opposite the
/// face, and the B-rep face tag.
pub type BoundaryTag = (u32, usize, u32);

impl TetComplex {
    /// A complex from vertices, tetrahedra (tet4 order, positive), boundary faces and feature
    /// edges. Faces shared by two tetrahedra become neighbours; every other face must be listed
    /// in `boundary`. Caps: `max_points` vertices and `max_tets` tetrahedra at any time.
    pub fn new(
        points: &[[f64; 3]],
        kinds: Vec<VertexKind>,
        tets: &[[u32; 4]],
        boundary: &[BoundaryTag],
        curves: &[((u32, u32), u32)],
        max_points: usize,
        max_tets: usize,
    ) -> Result<Self> {
        Self::build(points, kinds, tets, boundary, curves, max_points, max_tets, true)
    }

    /// [`Self::new`] for a mesh that may have flat or inverted tetrahedra (a tetrahedralization
    /// rounded to doubles, to be repaired by the operations, which only create positive ones).
    pub fn with_invalid(
        points: &[[f64; 3]],
        kinds: Vec<VertexKind>,
        tets: &[[u32; 4]],
        boundary: &[BoundaryTag],
        curves: &[((u32, u32), u32)],
        max_points: usize,
        max_tets: usize,
    ) -> Result<Self> {
        Self::build(points, kinds, tets, boundary, curves, max_points, max_tets, false)
    }

    #[allow(clippy::too_many_arguments)]
    fn build(
        points: &[[f64; 3]],
        kinds: Vec<VertexKind>,
        tets: &[[u32; 4]],
        boundary: &[BoundaryTag],
        curves: &[((u32, u32), u32)],
        max_points: usize,
        max_tets: usize,
        strict: bool,
    ) -> Result<Self> {
        if points.len() > max_points || points.len() >= NONE as usize - 1 {
            return Err(Error::Limit(format!("{} vertices; the limit is {max_points}", points.len())));
        }
        if tets.len() > max_tets || tets.len() >= DEAD as usize {
            return Err(Error::Limit(format!("{} tetrahedra; the limit is {max_tets}", tets.len())));
        }
        if kinds.len() != points.len() {
            return Err(Error::Invalid("one vertex kind per point".into()));
        }
        let pts = points.iter().map(|p| Point3::try_from(*p)).collect::<Result<Vec<_>>>()?;
        let mut cx = TetComplex {
            vtets: vec![Vec::new(); pts.len()],
            pts,
            kinds,
            tets: Vec::with_capacity(tets.len()),
            adj: vec![[NONE; 4]; tets.len()],
            tags: vec![[NONE; 4]; tets.len()],
            free: Vec::new(),
            curves: IdMap::default(),
            frozen: Vec::new(),
            mark: vec![0; tets.len()],
            epoch: 0,
            live: tets.len(),
            max_tets,
            max_points,
            journal: Journal::default(),
        };
        let n = cx.pts.len() as u32;
        let mut faces: IdMap<[u32; 3], (u32, usize)> = IdMap::default();
        for (t, vs) in tets.iter().enumerate() {
            if vs.iter().any(|&v| v >= n) || (0..4).any(|i| (i + 1..4).any(|j| vs[i] == vs[j])) {
                return Err(Error::Invalid(format!("tetrahedron {t} has a bad vertex list {vs:?}")));
            }
            if strict && cx.orient(*vs)? != Sign::Positive {
                return Err(Error::Invalid(format!("tetrahedron {t} is not positive")));
            }
            cx.tets.push(*vs);
            cx.link(t as u32, *vs);
            for k in 0..4 {
                let key = sorted3(face_toward(*vs, k));
                match faces.remove(&key) {
                    Some((u, j)) => {
                        cx.set_adj(t as u32, k, u)?;
                        cx.set_adj(u, j, t as u32)?;
                    }
                    None => {
                        faces.insert(key, (t as u32, k));
                    }
                }
            }
        }
        for &(t, k, tag) in boundary {
            let vs = cx.tet(t)?;
            if k > 3 || cx.adj(t)?.get(k).copied() != Some(NONE) {
                return Err(Error::Invalid(format!("boundary face {k} of tetrahedron {t} is not on the boundary")));
            }
            faces.remove(&sorted3(face_toward(vs, k)));
            if let Some(s) = cx.tags.get_mut(t as usize).and_then(|f| f.get_mut(k)) {
                *s = tag;
            }
        }
        if !faces.is_empty() {
            return Err(Error::Invalid(format!("{} faces have one tetrahedron and no boundary tag", faces.len())));
        }
        for &(e, c) in curves {
            cx.curves.insert(sorted2(e.0, e.1), c);
        }
        Ok(cx)
    }

    // ---- access ----

    pub fn point_count(&self) -> usize {
        self.pts.len()
    }

    pub fn point(&self, v: u32) -> Result<Point3> {
        self.pts.get(v as usize).copied().ok_or_else(|| corrupt("vertex id"))
    }

    pub fn xyz(&self, v: u32) -> Result<[f64; 3]> {
        self.point(v).map(|p| p.coords())
    }

    pub fn coords(&self, vs: [u32; 4]) -> Result<[[f64; 3]; 4]> {
        Ok([self.xyz(vs[0])?, self.xyz(vs[1])?, self.xyz(vs[2])?, self.xyz(vs[3])?])
    }

    pub fn kind(&self, v: u32) -> VertexKind {
        self.kinds.get(v as usize).copied().unwrap_or(VertexKind::Corner)
    }

    pub fn set_kind(&mut self, v: u32, k: VertexKind) {
        if let Some(s) = self.kinds.get_mut(v as usize) {
            *s = k;
        }
    }

    pub fn slots(&self) -> usize {
        self.tets.len()
    }

    pub fn live_count(&self) -> usize {
        self.live
    }

    pub fn max_tets(&self) -> usize {
        self.max_tets
    }

    pub fn max_points(&self) -> usize {
        self.max_points
    }

    pub fn is_live(&self, t: u32) -> bool {
        self.tets.get(t as usize).is_some_and(|vs| vs[0] != DEAD)
    }

    pub fn tet(&self, t: u32) -> Result<[u32; 4]> {
        self.tets.get(t as usize).copied().filter(|vs| vs[0] != DEAD).ok_or_else(|| corrupt("tetrahedron id"))
    }

    pub fn adj(&self, t: u32) -> Result<[u32; 4]> {
        self.adj.get(t as usize).copied().ok_or_else(|| corrupt("tetrahedron id"))
    }

    pub fn tags(&self, t: u32) -> Result<[u32; 4]> {
        self.tags.get(t as usize).copied().ok_or_else(|| corrupt("tetrahedron id"))
    }

    pub fn live_tets(&self) -> impl Iterator<Item = (u32, [u32; 4])> + '_ {
        self.tets.iter().enumerate().filter(|(_, vs)| vs[0] != DEAD).map(|(t, vs)| (t as u32, *vs))
    }

    /// Whether a vertex is still used by a tetrahedron.
    pub fn is_used(&self, v: u32) -> bool {
        self.vtets.get(v as usize).is_some_and(|l| !l.is_empty())
    }

    /// The curve of a feature edge.
    pub fn curve(&self, a: u32, b: u32) -> Option<u32> {
        self.curves.get(&sorted2(a, b)).copied()
    }

    pub fn curve_edges(&self) -> Vec<((u32, u32), u32)> {
        let mut e: Vec<_> = self.curves.iter().map(|(k, v)| (*k, *v)).collect();
        e.sort_unstable();
        e
    }

    /// Keep the boundary faces tagged with `tags` exactly as they are: [`Self::replace`] refuses
    /// every change that would remove one of them (a split, flip or collapse) or add one. Their
    /// vertices must also be [`VertexKind::Corner`]s, so nothing moves them.
    pub fn freeze_faces(&mut self, tags: &[u32]) {
        self.frozen = tags.to_vec();
        self.frozen.sort_unstable();
        self.frozen.dedup();
    }

    fn is_frozen(&self, tag: u32) -> bool {
        self.frozen.binary_search(&tag).is_ok()
    }

    pub fn set_curve(&mut self, a: u32, b: u32, curve: Option<u32>) {
        let key = sorted2(a, b);
        let before = match curve {
            Some(c) => self.curves.insert(key, c),
            None => self.curves.remove(&key),
        };
        if self.journal.recording {
            self.journal.steps.push(Step::Curve { edge: key, curve: before });
        }
    }

    // ---- predicates ----

    /// The exact orientation of four vertices (positive: tet4 order).
    pub fn orient(&self, [a, b, c, d]: [u32; 4]) -> Result<Sign> {
        Ok(tet_volume_sign(&self.point(a)?, &self.point(b)?, &self.point(c)?, &self.point(d)?))
    }

    /// The exact orientation of three vertices and a point.
    pub fn orient_with(&self, [a, b, c]: [u32; 3], p: &Point3) -> Result<Sign> {
        Ok(tet_volume_sign(&self.point(a)?, &self.point(b)?, &self.point(c)?, p))
    }

    /// Whether `p` lies strictly inside the circumsphere of the positive tetrahedron `vs`.
    pub fn in_sphere(&self, vs: [u32; 4], p: &Point3) -> Result<bool> {
        let [a, b, c, d] = [self.point(vs[0])?, self.point(vs[1])?, self.point(vs[2])?, self.point(vs[3])?];
        // Shewchuk's insphere is positive inside for orient3d(a, b, c, d) > 0, which is the
        // negative tet4 orientation.
        Ok(insphere(&a, &b, &c, &d, p) == Sign::Negative)
    }

    // ---- traversal ----

    fn next_epoch(&mut self) {
        self.epoch = self.epoch.wrapping_add(1);
        if self.epoch == 0 {
            self.mark.iter_mut().for_each(|m| *m = 0);
            self.epoch = 1;
        }
    }

    fn visited(&self, t: u32) -> bool {
        self.mark.get(t as usize).is_some_and(|&m| m == self.epoch)
    }

    fn visit(&mut self, t: u32) {
        if let Some(m) = self.mark.get_mut(t as usize) {
            *m = self.epoch;
        }
    }

    /// A live tetrahedron containing vertex `v`.
    pub fn incident(&self, v: u32) -> Result<u32> {
        self.vtets.get(v as usize).and_then(|l| l.first()).copied().ok_or_else(|| corrupt("vertex without a tetrahedron"))
    }

    /// The tetrahedra containing vertex `v`, by slot.
    pub fn star(&mut self, v: u32) -> Result<Vec<u32>> {
        let mut out = self.vtets.get(v as usize).cloned().unwrap_or_default();
        if out.is_empty() {
            return Err(corrupt("vertex without a tetrahedron"));
        }
        out.sort_unstable();
        Ok(out)
    }

    fn link(&mut self, t: u32, vs: [u32; 4]) {
        for v in vs {
            if let Some(l) = self.vtets.get_mut(v as usize) {
                l.push(t);
            }
        }
    }

    fn unlink(&mut self, t: u32, vs: [u32; 4]) {
        for v in vs {
            if let Some(l) = self.vtets.get_mut(v as usize)
                && let Some(i) = l.iter().position(|&u| u == t)
            {
                l.swap_remove(i);
            }
        }
    }

    /// The vertices joined to `v` by an edge, sorted.
    pub fn neighbours(&mut self, v: u32) -> Result<Vec<u32>> {
        let mut out = Vec::new();
        for t in self.star(v)? {
            out.extend(self.tet(t)?.into_iter().filter(|&u| u != v));
        }
        out.sort_unstable();
        out.dedup();
        Ok(out)
    }

    /// A tetrahedron containing edge `ab`, if the edge exists.
    pub fn edge_tet(&mut self, a: u32, b: u32) -> Result<Option<u32>> {
        let t = self.incident(a)?;
        if self.tet(t)?.contains(&b) {
            return Ok(Some(t));
        }
        Ok(self.star(a)?.into_iter().find(|&t| self.tet(t).is_ok_and(|vs| vs.contains(&b))))
    }

    /// The ring of tetrahedra around edge `ab`: the link vertices `v_0, …, v_{m−1}` and the
    /// tetrahedra `t_i = (a, b, v_i, v_{i+1})` (positive in that order). For an interior edge the
    /// ring is closed (`t_{m−1}` joins `v_{m−1}` to `v_0`, `closed` is true); for a boundary edge
    /// it is open: `m` link vertices and `m − 1` tetrahedra, from boundary face to boundary face.
    pub fn edge_ring(&mut self, a: u32, b: u32) -> Result<Option<EdgeRing>> {
        let Some(t0) = self.edge_tet(a, b)? else { return Ok(None) };
        self.edge_ring_at(t0, a, b)
    }

    /// [`Self::edge_ring`] from a tetrahedron `t0` known to hold the edge. `None` when the
    /// tetrahedra around the edge do not form one fan (an edge pinched between two groups), where
    /// no operation around the edge applies.
    pub fn edge_ring_at(&mut self, t0: u32, a: u32, b: u32) -> Result<Option<EdgeRing>> {
        // Orient the first tetrahedron as (a, b, c, d) positive.
        let (c, d) = {
            let vs = self.tet(t0)?;
            let others: Vec<u32> = vs.iter().copied().filter(|&v| v != a && v != b).collect();
            let (&c, &d) = match (others.first(), others.get(1)) {
                (Some(c), Some(d)) => (c, d),
                _ => return Ok(None),
            };
            if self.orient([a, b, c, d])? == Sign::Positive { (c, d) } else { (d, c) }
        };
        let mut verts = vec![c, d];
        let mut tets = vec![t0];
        let cap = self.max_tets.min(1 << 16);
        // Forward: across the face (a, b, v_{i+1}) opposite v_i.
        let mut closed = false;
        let (mut t, mut prev, mut cur) = (t0, c, d);
        loop {
            let vs = self.tet(t)?;
            let Some(k) = vs.iter().position(|&v| v == prev) else { return Ok(None) };
            let n = self.adj(t)?.get(k).copied().unwrap_or(NONE);
            if n == NONE {
                break;
            }
            if n == t0 {
                closed = true;
                break;
            }
            let nv = self.tet(n)?;
            let Some(next) = nv.iter().copied().find(|&v| v != a && v != b && v != cur) else { return Ok(None) };
            tets.push(n);
            if next != c {
                verts.push(next);
            }
            prev = cur;
            cur = next;
            t = n;
            if tets.len() > cap {
                return Ok(None);
            }
        }
        if closed {
            // The last tetrahedron joins the last link vertex back to v_0.
            if verts.len() != tets.len() {
                return Ok(None);
            }
            let mut ring = EdgeRing { a, b, verts, tets, closed: true };
            self.orient_ring(&mut ring)?;
            return Ok(Some(ring));
        }
        // Backward from t0: across the face (a, b, v_0) opposite v_1.
        let (mut t, mut prev, mut cur) = (t0, d, c);
        loop {
            let vs = self.tet(t)?;
            let Some(k) = vs.iter().position(|&v| v == prev) else { return Ok(None) };
            let n = self.adj(t)?.get(k).copied().unwrap_or(NONE);
            if n == NONE {
                break;
            }
            let nv = self.tet(n)?;
            let Some(next) = nv.iter().copied().find(|&v| v != a && v != b && v != cur) else { return Ok(None) };
            tets.insert(0, n);
            verts.insert(0, next);
            prev = cur;
            cur = next;
            t = n;
            if tets.len() > cap {
                return Ok(None);
            }
        }
        let mut ring = EdgeRing { a, b, verts, tets, closed: false };
        self.orient_ring(&mut ring)?;
        Ok(Some(ring))
    }

    /// Make the ring run so that `(a, b, v_i, v_{i+1})` is positive, judged by its first
    /// tetrahedron that is not flat (the first one may be flat or inverted in a mesh being
    /// repaired).
    fn orient_ring(&self, ring: &mut EdgeRing) -> Result<()> {
        let m = ring.verts.len();
        for i in 0..ring.tets.len() {
            let (Some(&p), Some(&q)) = (ring.verts.get(i), ring.verts.get((i + 1) % m)) else { break };
            match self.orient([ring.a, ring.b, p, q])? {
                Sign::Positive => return Ok(()),
                Sign::Negative => {
                    ring.verts.reverse();
                    if ring.closed {
                        // t_i joined v_i and v_{i+1}; reversed, the wrap tetrahedron stays last.
                        let last = ring.tets.pop();
                        ring.tets.reverse();
                        ring.tets.extend(last);
                    } else {
                        ring.tets.reverse();
                    }
                    return Ok(());
                }
                Sign::Zero => {}
            }
        }
        Ok(())
    }

    /// Walk from tetrahedron `start` toward `p` (visibility walk, faces tried in a fixed
    /// pseudo-random order so it cannot cycle forever), at most `cap` steps.
    pub fn walk(&self, start: u32, p: &Point3, cap: usize) -> Result<Option<Walk>> {
        let mut t = start;
        let mut rng: u32 = 0x9E37_79B9;
        for _ in 0..cap {
            let vs = self.tet(t)?;
            rng ^= rng << 13;
            rng ^= rng >> 17;
            rng ^= rng << 5;
            let first = (rng % 4) as usize;
            let mut moved = false;
            for i in 0..4 {
                let k = (first + i) % 4;
                if self.orient_with(face_toward(vs, k), p)? == Sign::Negative {
                    let n = self.adj(t)?.get(k).copied().unwrap_or(NONE);
                    if n == NONE {
                        return Ok(Some(Walk::Outside(t, k)));
                    }
                    t = n;
                    moved = true;
                    break;
                }
            }
            if !moved {
                return Ok(Some(Walk::Inside(t)));
            }
        }
        Ok(None)
    }

    // ---- mutation ----

    fn set_adj(&mut self, t: u32, k: usize, n: u32) -> Result<()> {
        let slot = self.adj.get_mut(t as usize).and_then(|a| a.get_mut(k)).ok_or_else(|| corrupt("neighbour slot"))?;
        *slot = n;
        Ok(())
    }

    /// Add a vertex (no tetrahedron uses it yet).
    pub fn add_point(&mut self, p: Point3, kind: VertexKind) -> Result<u32> {
        if self.pts.len() >= self.max_points {
            return Err(Error::Limit(format!("{} vertices; the limit is {}", self.pts.len(), self.max_points)));
        }
        let v = self.pts.len() as u32;
        self.pts.push(p);
        self.kinds.push(kind);
        self.vtets.push(Vec::new());
        if self.journal.recording {
            self.journal.steps.push(Step::Added { v });
        }
        Ok(v)
    }

    /// The data a change of `slots` and of the incidence of `verts` would overwrite.
    fn snapshot(&self, slots: &[u32], verts: &[u32]) -> Raw {
        let mut ids: Vec<u32> = slots.to_vec();
        ids.sort_unstable();
        ids.dedup();
        let mut vs: Vec<u32> = verts.to_vec();
        vs.sort_unstable();
        vs.dedup();
        Raw {
            slots: ids
                .into_iter()
                .filter_map(|t| {
                    let i = t as usize;
                    Some((t, *self.tets.get(i)?, *self.adj.get(i)?, *self.tags.get(i)?))
                })
                .collect(),
            tets_len: self.tets.len(),
            free: self.free.clone(),
            vtets: vs.into_iter().filter_map(|v| Some((v, self.vtets.get(v as usize)?.clone()))).collect(),
            live: self.live,
        }
    }

    /// Put back what [`Self::snapshot`] saved.
    fn restore(&mut self, raw: Raw) {
        for (v, l) in raw.vtets {
            if let Some(s) = self.vtets.get_mut(v as usize) {
                *s = l;
            }
        }
        for (t, vs, adj, tags) in raw.slots {
            let i = t as usize;
            if let (Some(a), Some(b), Some(c)) = (self.tets.get_mut(i), self.adj.get_mut(i), self.tags.get_mut(i)) {
                (*a, *b, *c) = (vs, adj, tags);
            }
        }
        self.tets.truncate(raw.tets_len);
        self.adj.truncate(raw.tets_len);
        self.tags.truncate(raw.tets_len);
        self.mark.truncate(raw.tets_len);
        self.free = raw.free;
        self.live = raw.live;
    }

    /// Remove the last vertex added, when no tetrahedron uses it (an insertion that failed).
    pub fn drop_last_point(&mut self, v: u32) -> Result<()> {
        if v as usize + 1 == self.pts.len() && !self.is_used(v) {
            self.pts.pop();
            self.kinds.pop();
            self.vtets.pop();
            return Ok(());
        }
        Err(corrupt("dropping a vertex in use"))
    }

    /// Move a vertex. The caller has checked that every tetrahedron around it stays positive
    /// ([`Self::move_is_valid`]).
    pub fn move_point(&mut self, v: u32, p: Point3) -> Result<()> {
        let slot = self.pts.get_mut(v as usize).ok_or_else(|| corrupt("vertex id"))?;
        let from = *slot;
        *slot = p;
        if self.journal.recording {
            self.journal.steps.push(Step::Moved { v, from });
        }
        Ok(())
    }

    /// Whether every tetrahedron of `star` (around `v`) stays positive and admissible with `v`
    /// at `p`.
    pub fn move_is_valid(&self, v: u32, p: &Point3, star: &[u32]) -> Result<bool> {
        for &t in star {
            let vs = self.tet(t)?;
            let k = vs.iter().position(|&u| u == v).ok_or_else(|| corrupt("star"))?;
            if self.orient_with(face_toward(vs, k), p)? != Sign::Positive {
                return Ok(false);
            }
            let mut x = self.coords(vs)?;
            if let Some(c) = x.get_mut(k) {
                *c = p.coords();
            }
            if !(normalized_jacobian(&tet4_jacobian(&x)) > ADMISSIBLE_MARGIN) {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn alloc(&mut self, vs: [u32; 4]) -> Result<u32> {
        if let Some(t) = self.free.pop() {
            let slot = self.tets.get_mut(t as usize).ok_or_else(|| corrupt("free slot"))?;
            *slot = vs;
            if let Some(a) = self.adj.get_mut(t as usize) {
                *a = [NONE; 4];
            }
            if let Some(a) = self.tags.get_mut(t as usize) {
                *a = [NONE; 4];
            }
            self.link(t, vs);
            return Ok(t);
        }
        if self.tets.len() >= self.max_tets || self.tets.len() >= DEAD as usize - 1 {
            return Err(Error::Limit(format!("more than {} tetrahedra", self.max_tets)));
        }
        self.tets.push(vs);
        self.adj.push([NONE; 4]);
        self.tags.push([NONE; 4]);
        self.mark.push(0);
        let t = (self.tets.len() - 1) as u32;
        self.link(t, vs);
        Ok(t)
    }

    /// Replace tetrahedra (see the module docs). `Ok(None)` when the change is not a valid
    /// replacement; nothing is changed then. `Ok(Some(ids))`: the new tetrahedra's ids, in the
    /// order of `change.new`.
    pub fn replace(&mut self, change: &Change) -> Result<Option<Vec<u32>>> {
        let old = &change.old;
        if old.is_empty() || change.new.is_empty() {
            return Ok(None);
        }
        let n_pts = self.pts.len() as u32;
        // Old region: its outside faces, keyed by sorted vertices.
        self.next_epoch();
        for &t in old {
            if !self.is_live(t) || self.visited(t) {
                return Ok(None);
            }
            self.visit(t);
        }
        // key -> (neighbour outside or NONE, the neighbour's slot pointing back, boundary tag,
        // the face as the old tetrahedron saw it)
        let mut outside: IdMap<[u32; 3], (u32, usize, u32, [u32; 3])> = IdMap::default();
        let mut old_edges: Vec<(u32, u32)> = Vec::new();
        for &t in old {
            let vs = self.tet(t)?;
            let ns = self.adj(t)?;
            let tg = self.tags(t)?;
            for &(i, j) in &EDGES {
                old_edges.push(sorted2(vs[i], vs[j]));
            }
            for k in 0..4 {
                let n = ns[k];
                if n != NONE && self.visited(n) {
                    continue;
                }
                let key = sorted3(face_toward(vs, k));
                let back = if n == NONE {
                    0
                } else {
                    let nv = self.tet(n)?;
                    let nn = self.adj(n)?;
                    (0..4).find(|&j| nn[j] == t && sorted3(face_toward(nv, j)) == key).ok_or_else(|| corrupt("neighbour does not point back"))?
                };
                if outside.insert(key, (n, back, tg[k], face_toward(vs, k))).is_some() {
                    return Ok(None);
                }
            }
        }
        let neighbours: Vec<u32> = outside.values().map(|e| e.0).filter(|&n| n != NONE).collect();
        // New tetrahedra: positive, distinct vertices.
        for vs in &change.new {
            if vs.iter().any(|&v| v >= n_pts) || (0..4).any(|i| (i + 1..4).any(|j| vs[i] == vs[j])) {
                return Ok(None);
            }
            if self.orient(*vs)? != Sign::Positive || !(normalized_jacobian(&tet4_jacobian(&self.coords(*vs)?)) > ADMISSIBLE_MARGIN) {
                return Ok(None);
            }
        }
        if self.live - old.len() + change.new.len() > self.max_tets {
            return Err(Error::Limit(format!("more than {} tetrahedra", self.max_tets)));
        }
        // Match the new faces.
        let mut inner: IdMap<[u32; 3], (usize, usize)> = IdMap::default();
        let mut links: Vec<(usize, usize, Link)> = Vec::with_capacity(change.new.len() * 4);
        let declared: IdMap<[u32; 3], ([u32; 3], u32)> = change.new_boundary.iter().map(|(f, tag)| (sorted3(*f), (*f, *tag))).collect();
        if declared.len() != change.new_boundary.len() {
            return Ok(None);
        }
        let mut used_declared = 0usize;
        for (i, vs) in change.new.iter().enumerate() {
            for k in 0..4 {
                let key = sorted3(face_toward(*vs, k));
                if let Some((j, l)) = inner.remove(&key) {
                    // Two tetrahedra on opposite sides of their shared face.
                    let other = change.new.get(j).map(|u| face_toward(*u, l)).unwrap_or([NONE; 3]);
                    if !same_cycle(face_toward(*vs, k), [other[0], other[2], other[1]]) {
                        return Ok(None);
                    }
                    links.push((i, k, Link::New(j)));
                    links.push((j, l, Link::New(i)));
                    continue;
                }
                inner.insert(key, (i, k));
            }
        }
        let mut pending: Vec<([u32; 3], (usize, usize))> = inner.into_iter().collect();
        pending.sort_unstable();
        for (key, (i, k)) in pending {
            let face = change.new.get(i).map(|vs| face_toward(*vs, k)).unwrap_or([NONE; 3]);
            let declared_face = declared.get(&key).copied();
            used_declared += usize::from(declared_face.is_some());
            if let Some((n, back, tag, seen)) = outside.remove(&key) {
                // The new tetrahedron lies on the side of the face the old one did.
                if !same_cycle(face, seen) || (n == NONE && self.is_frozen(tag) && declared_face.is_some_and(|d| d.1 != tag)) {
                    return Ok(None);
                }
                links.push((i, k, if n == NONE { Link::Boundary(declared_face.map_or(tag, |d| d.1)) } else { Link::Old(n, back) }));
            } else if let Some((oriented, tag)) = declared_face {
                // A new boundary face, oriented as declared (into the mesh).
                if !same_cycle(face, oriented) || self.is_frozen(tag) {
                    return Ok(None);
                }
                links.push((i, k, Link::Boundary(tag)));
            } else {
                return Ok(None);
            }
        }
        if used_declared != declared.len() {
            return Ok(None);
        }
        // What is left outside must be boundary faces the change removes, and none frozen.
        if outside.values().any(|(n, _, tag, _)| *n != NONE || self.is_frozen(*tag)) || (!outside.is_empty() && change.new_boundary.is_empty()) {
            return Ok(None);
        }
        // Feature edges strictly inside the old region must survive (or be the split one).
        if !self.curves.is_empty() {
            let mut kept: Vec<(u32, u32)> = change.new.iter().flat_map(|vs| EDGES.iter().map(move |&(i, j)| sorted2(vs[i], vs[j]))).collect();
            kept.sort_unstable();
            for e in &old_edges {
                if self.curves.contains_key(e) && !change.drop_curves.iter().any(|&(a, b)| sorted2(a, b) == *e) && kept.binary_search(e).is_err() {
                    return Ok(None);
                }
            }
        }
        // A vertex the change removes (in no new tetrahedron) must have every tetrahedron around
        // it in the old region; otherwise a tetrahedron outside would keep a vertex whose faces
        // were taken away.
        let mut kept_vertices: Vec<u32> = change.new.iter().flatten().copied().collect();
        kept_vertices.sort_unstable();
        kept_vertices.dedup();
        let mut gone: Vec<u32> = Vec::new();
        for &t in old {
            for v in self.tet(t)? {
                if kept_vertices.binary_search(&v).is_err() && !gone.contains(&v) {
                    gone.push(v);
                }
            }
        }
        if !gone.is_empty() {
            let region: IdMap<u32, ()> = old.iter().map(|&t| (t, ())).collect();
            for v in gone {
                for t in self.star(v)? {
                    if !region.contains_key(&t) {
                        return Ok(None);
                    }
                }
            }
        }
        // Commit (saving what it overwrites when the journal records).
        let raw = if self.journal.recording {
            let mut slots: Vec<u32> = old.clone();
            slots.extend(neighbours.iter().copied());
            // Slots the allocations below may take from the free list.
            slots.extend(self.free.iter().rev().take(change.new.len()).copied());
            let mut verts: Vec<u32> = change.new.iter().flatten().copied().collect();
            for &t in old {
                verts.extend(self.tet(t)?);
            }
            Some(Box::new(self.snapshot(&slots, &verts)))
        } else {
            None
        };
        let mut removed = Vec::with_capacity(old.len());
        for &t in old {
            let vs = self.tet(t)?;
            removed.push(Removed { vs, tags: self.tags(t)? });
            self.unlink(t, vs);
            if let Some(s) = self.tets.get_mut(t as usize) {
                *s = [DEAD; 4];
            }
        }
        // Recycle in a fixed order: the lowest slots are handed out first.
        let mut freed: Vec<u32> = old.clone();
        freed.sort_unstable_by(|a, b| b.cmp(a));
        self.free.extend(freed);
        let mut ids = Vec::with_capacity(change.new.len());
        for vs in &change.new {
            ids.push(self.alloc(*vs)?);
        }
        for (i, k, link) in links {
            let t = *ids.get(i).ok_or_else(|| corrupt("new tetrahedron"))?;
            match link {
                Link::New(j) => {
                    let u = *ids.get(j).ok_or_else(|| corrupt("new tetrahedron"))?;
                    self.set_adj(t, k, u)?;
                }
                Link::Old(n, back) => {
                    self.set_adj(t, k, n)?;
                    self.set_adj(n, back, t)?;
                }
                Link::Boundary(tag) => {
                    self.set_adj(t, k, NONE)?;
                    if let Some(s) = self.tags.get_mut(t as usize).and_then(|f| f.get_mut(k)) {
                        *s = tag;
                    }
                }
            }
        }
        self.live = self.live - old.len() + ids.len();
        if self.journal.recording {
            self.journal.steps.push(Step::Replaced { removed, added: change.new.clone(), raw: raw.unwrap_or_default() });
        }
        Ok(Some(ids))
    }

    /// Delete a tetrahedron that has exactly two boundary faces with the same tag meeting at an
    /// edge that is not a feature edge: its two other faces become boundary faces with that tag
    /// (the surface is retriangulated across the other diagonal of the quadrilateral). Only for a
    /// flat tetrahedron, whose removal does not change the domain; the caller decides that.
    /// Returns whether it was deleted.
    pub fn delete_boundary_tet(&mut self, t: u32) -> Result<bool> {
        if !self.is_live(t) {
            return Ok(false);
        }
        let vs = self.tet(t)?;
        let ns = self.adj(t)?;
        let tags = self.tags(t)?;
        let bnd: Vec<usize> = (0..4).filter(|&k| ns[k] == NONE).collect();
        let &[i, j] = bnd.as_slice() else { return Ok(false) };
        let tag = tags[i];
        if tags[j] != tag || tag == NONE {
            return Ok(false);
        }
        // The two boundary faces share the edge of the two other vertices.
        let inner: Vec<usize> = (0..4).filter(|&k| k != i && k != j).collect();
        let &[k, l] = inner.as_slice() else { return Ok(false) };
        if self.curve(vs[k], vs[l]).is_some() {
            return Ok(false);
        }
        // Each neighbour gains a boundary face (both, when they are one tetrahedron); one left
        // with every face on the boundary would float free of the mesh, joined to it by its
        // vertices alone (a singular stiffness).
        for m in [k, l] {
            let n = ns[m];
            let gained = if ns[k] == ns[l] { 2 } else { 1 };
            let open = self.adj(n)?.iter().filter(|&&x| x == NONE).count();
            if open + gained >= 4 {
                return Ok(false);
            }
        }
        let raw = self.journal.recording.then(|| Box::new(self.snapshot(&[t, ns[k], ns[l]], &vs)));
        for m in [k, l] {
            let n = ns[m];
            let nn = self.adj(n)?;
            let back = (0..4).find(|&x| nn[x] == t).ok_or_else(|| corrupt("neighbour does not point back"))?;
            self.set_adj(n, back, NONE)?;
            if let Some(s) = self.tags.get_mut(n as usize).and_then(|f| f.get_mut(back)) {
                *s = tag;
            }
        }
        self.unlink(t, vs);
        if let Some(s) = self.tets.get_mut(t as usize) {
            *s = [DEAD; 4];
        }
        self.free.push(t);
        self.live -= 1;
        if self.journal.recording {
            self.journal.steps.push(Step::Deleted { removed: Removed { vs, tags }, raw: raw.unwrap_or_default() });
        }
        Ok(true)
    }

    /// Check the data structure: neighbours point back across the same face (seen from opposite
    /// sides), boundary faces have no neighbour, every used vertex's tetrahedron holds it, the
    /// live count is right. For tests and diagnostics.
    pub fn validate(&self) -> Result<()> {
        let mut live = 0;
        for (t, vs) in self.live_tets() {
            live += 1;
            if self.orient(vs)? != Sign::Positive {
                return Err(corrupt(&format!("tetrahedron {t} is not positive")));
            }
            let ns = self.adj(t)?;
            for k in 0..4 {
                let n = ns[k];
                if n == NONE {
                    continue;
                }
                let nv = self.tet(n).map_err(|_| corrupt(&format!("tetrahedron {t} has a dead neighbour {n}")))?;
                let nn = self.adj(n)?;
                let f = face_toward(vs, k);
                let Some(j) = (0..4).find(|&j| nn[j] == t) else { return Err(corrupt(&format!("neighbour {n} of {t} does not point back"))) };
                let g = face_toward(nv, j);
                if !same_cycle(f, [g[0], g[2], g[1]]) {
                    return Err(corrupt(&format!("tetrahedra {t} and {n} do not share face {f:?} oppositely ({g:?})")));
                }
            }
        }
        if live != self.live {
            return Err(corrupt(&format!("live count {} but {live} live tetrahedra", self.live)));
        }
        let mut held = 0usize;
        for (v, l) in self.vtets.iter().enumerate() {
            for &t in l {
                if !self.tet(t).is_ok_and(|vs| vs.contains(&(v as u32))) {
                    return Err(corrupt(&format!("vertex {v} lists tetrahedron {t} that does not hold it")));
                }
            }
            held += l.len();
        }
        if held != 4 * live {
            return Err(corrupt("vertex incidence does not match the tetrahedra"));
        }
        Ok(())
    }

    /// Undo every step recorded since the journal had `len` steps.
    pub fn rollback(&mut self, len: usize) -> Result<()> {
        let recording = self.journal.recording;
        self.journal.recording = false;
        while self.journal.steps.len() > len {
            let Some(step) = self.journal.steps.pop() else { break };
            match step {
                Step::Replaced { raw, .. } | Step::Deleted { raw, .. } => self.restore(*raw),
                Step::Moved { v, from } => {
                    if let Some(s) = self.pts.get_mut(v as usize) {
                        *s = from;
                    }
                }
                Step::Added { v } => {
                    if v as usize + 1 == self.pts.len() {
                        self.pts.pop();
                        self.kinds.pop();
                        self.vtets.pop();
                    }
                }
                Step::Curve { edge, curve } => {
                    match curve {
                        Some(c) => self.curves.insert(edge, c),
                        None => self.curves.remove(&edge),
                    };
                }
            }
        }
        self.journal.recording = recording;
        Ok(())
    }

    /// The live mesh, compacted: vertices still used (renumbered in order), tetrahedra in slot
    /// order, boundary faces `(tet, opposite slot, tag)`, and the old id of each vertex.
    pub fn compact(&self) -> Result<Compacted> {
        let mut id = vec![NONE; self.pts.len()];
        let mut points = Vec::new();
        let mut source = Vec::new();
        for (v, p) in self.pts.iter().enumerate() {
            if self.is_used(v as u32) {
                if let Some(s) = id.get_mut(v) {
                    *s = points.len() as u32;
                }
                points.push(p.coords());
                source.push(v as u32);
            }
        }
        let mut tets = Vec::with_capacity(self.live);
        let mut boundary = Vec::new();
        for (t, vs) in self.live_tets() {
            let mapped = vs.map(|v| id.get(v as usize).copied().unwrap_or(NONE));
            if mapped.contains(&NONE) {
                return Err(corrupt("a tetrahedron uses a removed vertex"));
            }
            let ns = self.adj(t)?;
            let tg = self.tags(t)?;
            for k in 0..4 {
                if ns[k] == NONE {
                    boundary.push((tets.len() as u32, k, tg[k]));
                }
            }
            tets.push(mapped);
        }
        let mut curves = Vec::new();
        for ((a, b), c) in self.curve_edges() {
            let (ia, ib) = (id.get(a as usize).copied().unwrap_or(NONE), id.get(b as usize).copied().unwrap_or(NONE));
            if ia != NONE && ib != NONE {
                curves.push(((ia, ib), c));
            }
        }
        Ok(Compacted { points, tets, boundary, curves, source })
    }
}

/// Whether `g` is `f` with its vertices rotated (the same oriented triangle).
fn same_cycle(f: [u32; 3], g: [u32; 3]) -> bool {
    g == f || g == [f[1], f[2], f[0]] || g == [f[2], f[0], f[1]]
}

/// How a new tetrahedron's face is matched in [`TetComplex::replace`].
#[derive(Debug, Clone, Copy)]
enum Link {
    New(usize),
    Old(u32, usize),
    Boundary(u32),
}

/// The ring around an edge ([`TetComplex::edge_ring`]).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EdgeRing {
    pub a: u32,
    pub b: u32,
    pub verts: Vec<u32>,
    pub tets: Vec<u32>,
    pub closed: bool,
}

/// A compacted mesh ([`TetComplex::compact`]).
#[derive(Debug, Clone, PartialEq)]
pub struct Compacted {
    pub points: Vec<[f64; 3]>,
    pub tets: Vec<[u32; 4]>,
    pub boundary: Vec<(u32, usize, u32)>,
    pub curves: Vec<((u32, u32), u32)>,
    pub source: Vec<u32>,
}
