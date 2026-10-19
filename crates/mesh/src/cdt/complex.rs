//! The working tetrahedralization of the boundary recovery: tetrahedra over explicit and implicit
//! ([`Lnc`](crate::predicates::Lnc)) vertices, with the operations the recovery needs: point
//! location, Bowyer–Watson insertion of Steiner points, vertex stars, edge lookup, and the
//! replacement of a set of tetrahedra by another filling the same region.
//!
//! Layout: per tetrahedron four vertex ids in tet4 order (`det[b − a, c − a, d − a] > 0`, decided
//! exactly) and four neighbour ids, neighbour `i` across the face opposite vertex `i`; [`NONE`] on
//! the hull. The hull is the bounding box the recovery adds around the surface, so every Steiner
//! point is strictly inside it and no ghost tetrahedra are needed.

pub(super) use crate::tet::complex::{IdMap, NONE, face_toward, sorted3};

use crate::delaunay::Triangulation;
use crate::predicates::{GenericPoint, Sign, generic_insphere_sos, generic_lex, generic_orient3d};
use crate::{Error, Result};

/// Marks a recycled tetrahedron slot.
const DEAD: u32 = u32::MAX - 1;

pub(super) fn corrupt(what: &str) -> Error {
    Error::Failed(format!("constrained tetrahedralization is inconsistent: {what}"))
}

/// The working tetrahedralization.
#[derive(Debug, Clone)]
pub(super) struct Complex {
    points: Vec<GenericPoint>,
    /// A box certainly containing each point (implicit points are not exactly representable).
    lo: Vec<[f64; 3]>,
    hi: Vec<[f64; 3]>,
    /// Rounded coordinates (floating-point estimates only, never decisions).
    approx: Vec<[f64; 3]>,
    tets: Vec<[u32; 4]>,
    adj: Vec<[u32; 4]>,
    free: Vec<u32>,
    vertex_tet: Vec<u32>,
    /// Traversal marks: a tetrahedron is visited in the current traversal when its mark equals
    /// `epoch`.
    mark: Vec<u32>,
    epoch: u32,
    rng: u64,
    /// Most tetrahedron slots allowed.
    max_tets: usize,
    live: usize,
}

/// Where a point lies ([`Complex::locate`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Located {
    /// In the closed tetrahedron and not one of its vertices.
    Inside(u32),
    /// Equal to this vertex.
    Vertex(u32),
}

impl Complex {
    /// The finite tetrahedra of a Delaunay triangulation of `points` (in the same vertex order),
    /// whose hull faces get no neighbour.
    pub(super) fn from_delaunay(tri: &Triangulation, points: Vec<GenericPoint>, max_tets: usize) -> Result<Self> {
        let mut id = vec![NONE; tri.slots()];
        let mut tets = Vec::new();
        for (t, vs) in tri.finite_tets() {
            if let Some(slot) = id.get_mut(t as usize) {
                *slot = tets.len() as u32;
            }
            tets.push(vs);
        }
        if tets.len() > max_tets {
            return Err(Error::Limit(format!("{} tetrahedra; the limit is {max_tets}", tets.len())));
        }
        let mut adj = Vec::with_capacity(tets.len());
        for (t, _) in tri.finite_tets() {
            let ns = tri.neighbors(t).ok_or_else(|| corrupt("Delaunay neighbours"))?;
            adj.push(ns.map(|n| id.get(n as usize).copied().unwrap_or(NONE)));
        }
        let (lo, hi) = points.iter().map(GenericPoint::bounds).unzip();
        let approx = points.iter().map(GenericPoint::approx).collect();
        let mut cx = Complex {
            vertex_tet: vec![NONE; points.len()],
            points,
            lo,
            hi,
            approx,
            mark: vec![0; tets.len()],
            live: tets.len(),
            tets,
            adj,
            free: Vec::new(),
            epoch: 0,
            rng: 0x9E37_79B9_7F4A_7C15,
            max_tets,
        };
        for (t, vs) in cx.tets.iter().enumerate() {
            for &v in vs {
                if let Some(s) = cx.vertex_tet.get_mut(v as usize) {
                    *s = t as u32;
                }
            }
        }
        Ok(cx)
    }

    // ---- read access ----

    pub(super) fn points(&self) -> &[GenericPoint] {
        &self.points
    }

    pub(super) fn pt(&self, v: u32) -> Result<&GenericPoint> {
        self.points.get(v as usize).ok_or_else(|| corrupt("vertex id"))
    }

    /// Vertex `v`'s coordinates rounded to doubles: for estimates and messages only.
    pub(super) fn approx(&self, v: u32) -> Result<[f64; 3]> {
        self.approx.get(v as usize).copied().ok_or_else(|| corrupt("vertex id"))
    }

    /// The box certainly containing vertex `v`.
    pub(super) fn bounds(&self, v: u32) -> ([f64; 3], [f64; 3]) {
        match (self.lo.get(v as usize), self.hi.get(v as usize)) {
            (Some(&l), Some(&h)) => (l, h),
            _ => ([f64::NEG_INFINITY; 3], [f64::INFINITY; 3]),
        }
    }

    /// The box certainly containing all of `vs`.
    pub(super) fn bounds_of(&self, vs: &[u32]) -> ([f64; 3], [f64; 3]) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for &v in vs {
            let (l, h) = self.bounds(v);
            for k in 0..3 {
                lo[k] = lo[k].min(l[k]);
                hi[k] = hi[k].max(h[k]);
            }
        }
        (lo, hi)
    }

    pub(super) fn slots(&self) -> usize {
        self.tets.len()
    }

    pub(super) fn live_count(&self) -> usize {
        self.live
    }

    /// The vertices of a live tetrahedron.
    pub(super) fn tet(&self, t: u32) -> Result<[u32; 4]> {
        self.tets.get(t as usize).copied().filter(|v| v[0] != DEAD).ok_or_else(|| corrupt("tetrahedron id"))
    }

    pub(super) fn is_live(&self, t: u32) -> bool {
        self.tets.get(t as usize).is_some_and(|v| v[0] != DEAD)
    }

    pub(super) fn adj(&self, t: u32) -> Result<[u32; 4]> {
        self.tet(t)?;
        self.adj.get(t as usize).copied().ok_or_else(|| corrupt("tetrahedron id"))
    }

    /// The live tetrahedra in slot order.
    pub(super) fn live_tets(&self) -> impl Iterator<Item = (u32, [u32; 4])> + '_ {
        self.tets.iter().enumerate().filter(|(_, v)| v[0] != DEAD).map(|(t, v)| (t as u32, *v))
    }

    /// The sign of `det[b − a, c − a, d − a]` (tet4 orientation), exact.
    pub(super) fn volume(&self, [a, b, c, d]: [u32; 4]) -> Result<Sign> {
        Ok(-generic_orient3d(self.pt(a)?, self.pt(b)?, self.pt(c)?, self.pt(d)?))
    }

    /// Whether `x` is inside the perturbed circumsphere of the positive tetrahedron `vs`.
    pub(super) fn in_sphere(&self, [a, b, c, d]: [u32; 4], x: u32) -> Result<bool> {
        // `insphere` wants Shewchuk-positive orientation (tet4-negative): swap the first two.
        Ok(generic_insphere_sos(self.pt(b)?, self.pt(a)?, self.pt(c)?, self.pt(d)?, self.pt(x)?) == Sign::Positive)
    }

    /// Whether two vertices are the same point (exact).
    fn same_point(&self, u: u32, v: u32) -> Result<bool> {
        Ok(generic_lex(self.pt(u)?, self.pt(v)?).is_eq())
    }

    // ---- traversal marks ----

    /// A new traversal: marks `epoch` (visited) and `epoch + 1` (rejected) are fresh.
    fn next_epoch(&mut self) {
        if self.epoch >= u32::MAX - 4 {
            self.mark.fill(0);
            self.epoch = 0;
        }
        self.epoch += 2;
    }

    fn visited(&self, t: u32) -> bool {
        self.mark.get(t as usize) == Some(&self.epoch)
    }

    fn visit(&mut self, t: u32) {
        if let Some(m) = self.mark.get_mut(t as usize) {
            *m = self.epoch;
        }
    }

    fn rejected(&self, t: u32) -> bool {
        self.mark.get(t as usize) == Some(&(self.epoch + 1))
    }

    fn reject(&mut self, t: u32) {
        if let Some(m) = self.mark.get_mut(t as usize) {
            *m = self.epoch + 1;
        }
    }

    // ---- stars and edges ----

    /// A live tetrahedron incident to `v`.
    fn incident(&self, v: u32) -> Result<u32> {
        let t = self.vertex_tet.get(v as usize).copied().unwrap_or(NONE);
        if self.tets.get(t as usize).is_some_and(|vs| vs[0] != DEAD && vs.contains(&v)) {
            return Ok(t);
        }
        // The hint is maintained on every allocation, so this scan is a safety net only.
        self.live_tets().find(|(_, vs)| vs.contains(&v)).map(|(t, _)| t).ok_or_else(|| corrupt("a vertex is in no tetrahedron"))
    }

    /// The tetrahedra incident to `v` (its star), into `out`.
    pub(super) fn star(&mut self, v: u32, out: &mut Vec<u32>) -> Result<()> {
        out.clear();
        let start = self.incident(v)?;
        self.next_epoch();
        self.visit(start);
        out.push(start);
        let mut i = 0;
        while let Some(&t) = out.get(i) {
            i += 1;
            let vs = self.tet(t)?;
            let ns = self.adj(t)?;
            for (&u, &n) in vs.iter().zip(&ns) {
                if u == v || n == NONE || self.visited(n) {
                    continue;
                }
                self.visit(n);
                out.push(n);
            }
            if out.len() > self.live {
                return Err(corrupt("a vertex star larger than the mesh"));
            }
        }
        Ok(())
    }

    /// Whether `u v` is an edge.
    pub(super) fn has_edge(&mut self, u: u32, v: u32, scratch: &mut Vec<u32>) -> Result<bool> {
        self.star(u, scratch)?;
        for &t in scratch.iter() {
            if self.tet(t)?.contains(&v) {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// The tetrahedra whose (conservative) box meets the box `lo..hi`, connected to the star of
    /// `v` through such tetrahedra; at most `cap` of them.
    pub(super) fn tets_near(&mut self, v: u32, lo: [f64; 3], hi: [f64; 3], cap: usize) -> Result<Vec<u32>> {
        let mut out = Vec::new();
        self.star(v, &mut out)?;
        self.next_epoch();
        for &t in &out {
            self.visit(t);
        }
        let mut i = 0;
        while let Some(&t) = out.get(i) {
            i += 1;
            if out.len() >= cap {
                break;
            }
            for n in self.adj(t)? {
                if n == NONE || self.visited(n) {
                    continue;
                }
                self.visit(n);
                let (l, h) = self.bounds_of(&self.tet(n)?);
                if (0..3).all(|k| l[k] <= hi[k] && h[k] >= lo[k]) {
                    out.push(n);
                }
            }
        }
        Ok(out)
    }

    // ---- point location ----

    /// The visibility walk towards `p` from `start` (Devillers, Pion, Teillaud 2002), capped at
    /// the number of tetrahedra, then a scan.
    pub(super) fn locate(&mut self, p: u32, start: u32) -> Result<Located> {
        let mut t = if self.is_live(start) { start } else { self.live_tets().next().map(|(t, _)| t).ok_or_else(|| corrupt("empty mesh"))? };
        let mut prev = NONE;
        for _ in 0..self.live.saturating_add(16) {
            let vs = self.tet(t)?;
            let ns = self.adj(t)?;
            self.rng ^= self.rng << 13;
            self.rng ^= self.rng >> 7;
            self.rng ^= self.rng << 17;
            let first = (self.rng >> 61) as usize;
            let mut next = None;
            for i in 0..4 {
                let k = (first + i) % 4;
                let n = ns.get(k).copied().unwrap_or(NONE);
                if (n != prev || n == NONE) && self.beyond(vs, k, p)? {
                    next = Some(n);
                    break;
                }
            }
            match next {
                None => return self.inside(vs, t, p),
                Some(NONE) => return Err(Error::Failed("a Steiner point lies outside the bounding box".into())),
                Some(n) => {
                    prev = t;
                    t = n;
                }
            }
        }
        // The walk may cycle in a non-Delaunay mesh: scan.
        let live: Vec<(u32, [u32; 4])> = self.live_tets().collect();
        for (t, vs) in live {
            let mut inside = true;
            for k in 0..4 {
                if self.beyond(vs, k, p)? {
                    inside = false;
                    break;
                }
            }
            if inside {
                return self.inside(vs, t, p);
            }
        }
        Err(Error::Failed("a Steiner point is in no tetrahedron".into()))
    }

    /// Whether `p` is strictly beyond the face of `vs` opposite slot `k`.
    fn beyond(&self, mut vs: [u32; 4], k: usize, p: u32) -> Result<bool> {
        if let Some(s) = vs.get_mut(k) {
            *s = p;
        }
        Ok(self.volume(vs)? == Sign::Negative)
    }

    fn inside(&self, vs: [u32; 4], t: u32, p: u32) -> Result<Located> {
        for v in vs {
            if self.same_point(v, p)? {
                return Ok(Located::Vertex(v));
            }
        }
        Ok(Located::Inside(t))
    }

    // ---- modification ----

    /// Add a point (not yet in any tetrahedron); its vertex id.
    pub(super) fn push_point(&mut self, p: GenericPoint) -> Result<u32> {
        if self.points.len() >= (NONE - 2) as usize {
            return Err(Error::Limit("too many vertices".into()));
        }
        let (l, h) = p.bounds();
        self.approx.push(p.approx());
        self.points.push(p);
        self.lo.push(l);
        self.hi.push(h);
        self.vertex_tet.push(NONE);
        Ok((self.points.len() - 1) as u32)
    }

    /// Insert vertex `v` (from [`Complex::push_point`]) by Bowyer–Watson with the perturbed
    /// in-sphere test: the mesh must be Delaunay for that predicate, which it is until the facet
    /// recovery starts. The walk starts at a tetrahedron incident to `near`. Returns the vertices
    /// of the deleted tetrahedra (the only ones whose edges can have changed).
    pub(super) fn insert(&mut self, v: u32, near: u32) -> Result<Vec<u32>> {
        let start = self.incident(near)?;
        let seed = match self.locate(v, start)? {
            Located::Vertex(u) => {
                let at = self.approx(u)?;
                return Err(Error::Measured(
                    crate::Kind::Failed,
                    crate::message![format!("a Steiner point coincides with vertex {u} at "), crate::message::point(at)],
                ));
            }
            Located::Inside(t) => t,
        };
        // The cavity: every tetrahedron whose perturbed circumsphere contains v, grown from the
        // one containing it (which always conflicts: v is in the closed tetrahedron, so strictly
        // inside its ball).
        self.next_epoch();
        let mut cavity = vec![seed];
        self.visit(seed);
        let mut boundary: Vec<(u32, usize)> = Vec::new();
        let mut i = 0;
        while let Some(&t) = cavity.get(i) {
            i += 1;
            for (j, n) in self.adj(t)?.into_iter().enumerate() {
                if n == NONE || self.rejected(n) {
                    boundary.push((t, j));
                    continue;
                }
                if self.visited(n) {
                    continue;
                }
                if self.in_sphere(self.tet(n)?, v)? {
                    self.visit(n);
                    cavity.push(n);
                } else {
                    self.reject(n);
                    boundary.push((t, j));
                }
            }
            if cavity.len() > self.live {
                return Err(corrupt("a cavity larger than the mesh"));
            }
        }
        let mut new = Vec::with_capacity(boundary.len());
        for (t, j) in boundary {
            let mut vs = self.tet(t)?;
            if let Some(s) = vs.get_mut(j) {
                *s = v;
            }
            if self.volume(vs)? != Sign::Positive {
                return Err(corrupt("a Bowyer–Watson cavity is not star-shaped from the new point"));
            }
            new.push(vs);
        }
        let mut touched: Vec<u32> = Vec::with_capacity(4 * cavity.len());
        for &t in &cavity {
            touched.extend(self.tet(t)?);
        }
        self.replace(&cavity, &new)?;
        touched.sort_unstable();
        touched.dedup();
        Ok(touched)
    }

    fn alloc(&mut self, vs: [u32; 4]) -> Result<u32> {
        let t = match self.free.pop() {
            Some(t) => {
                *self.tets.get_mut(t as usize).ok_or_else(|| corrupt("free slot"))? = vs;
                *self.adj.get_mut(t as usize).ok_or_else(|| corrupt("free slot"))? = [NONE; 4];
                t
            }
            None => {
                if self.tets.len() >= self.max_tets {
                    return Err(Error::Limit(format!("more than {} tetrahedra", self.max_tets)));
                }
                self.tets.push(vs);
                self.adj.push([NONE; 4]);
                self.mark.push(0);
                (self.tets.len() - 1) as u32
            }
        };
        for v in vs {
            if let Some(s) = self.vertex_tet.get_mut(v as usize) {
                *s = t;
            }
        }
        self.live += 1;
        Ok(t)
    }

    fn kill(&mut self, t: u32) -> Result<()> {
        *self.tets.get_mut(t as usize).ok_or_else(|| corrupt("tetrahedron id"))? = [DEAD; 4];
        self.free.push(t);
        self.live = self.live.saturating_sub(1);
        Ok(())
    }

    fn set_adj(&mut self, t: u32, i: usize, n: u32) -> Result<()> {
        *self.adj.get_mut(t as usize).and_then(|a| a.get_mut(i)).ok_or_else(|| corrupt("tetrahedron id"))? = n;
        Ok(())
    }

    /// Replace the tetrahedra `old` by `new` (positive tet4 tuples) filling the same region: each
    /// face of a new tetrahedron is either shared with exactly one other new tetrahedron or is a
    /// face of the region's boundary, where it takes over the outside neighbour. Anything else is
    /// reported, never patched. The new slot ids.
    pub(super) fn replace(&mut self, old: &[u32], new: &[[u32; 4]]) -> Result<Vec<u32>> {
        // Read the outside links before any slot is recycled.
        // Boundary face → (outside neighbour, its slot back, covered yet).
        let mut external: IdMap<[u32; 3], (u32, usize, bool)> = IdMap::with_capacity_and_hasher(old.len() * 2, Default::default());
        self.next_epoch();
        for &t in old {
            self.visit(t);
        }
        for &t in old {
            let vs = self.tet(t)?;
            for (j, n) in self.adj(t)?.into_iter().enumerate() {
                if n != NONE && self.visited(n) {
                    continue;
                }
                let back = if n == NONE { 0 } else { self.adj(n)?.iter().position(|&m| m == t).ok_or_else(|| corrupt("one-sided adjacency"))? };
                if external.insert(sorted3(face_toward(vs, j)), (n, back, false)).is_some() {
                    return Err(corrupt("a region boundary face is listed twice"));
                }
            }
        }
        for &t in old {
            self.kill(t)?;
        }
        let mut ids = Vec::with_capacity(new.len());
        let mut open: IdMap<[u32; 3], (u32, usize)> = IdMap::with_capacity_and_hasher(new.len() * 2, Default::default());
        let mut matched_external = 0usize;
        for &vs in new {
            let t = self.alloc(vs)?;
            ids.push(t);
            for k in 0..4 {
                let key = sorted3(face_toward(vs, k));
                if let Some(entry) = external.get_mut(&key) {
                    let (n, back, covered) = *entry;
                    if covered {
                        return Err(corrupt("a region boundary face is covered twice"));
                    }
                    entry.2 = true;
                    self.set_adj(t, k, n)?;
                    if n != NONE {
                        self.set_adj(n, back, t)?;
                    }
                    matched_external += 1;
                } else if let Some((u, l)) = open.remove(&key) {
                    self.set_adj(t, k, u)?;
                    self.set_adj(u, l, t)?;
                } else {
                    open.insert(key, (t, k));
                }
            }
        }
        if !open.is_empty() || matched_external != external.len() {
            return Err(corrupt(&format!(
                "a retriangulated region does not close: {} unmatched inner faces, {} of {} boundary faces covered",
                open.len(),
                matched_external,
                external.len()
            )));
        }
        Ok(ids)
    }
}
