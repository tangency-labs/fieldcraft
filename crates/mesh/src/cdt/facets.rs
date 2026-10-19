//! Facet recovery (Diazzi et al. 2023, §3.4 and §4.4). Every input triangle is a facet; once all
//! segments are edges, a facet is *missing* when an edge of the tetrahedralization pierces its
//! interior. The tetrahedra with such edges form face-connected *cavities* (each meets the
//! facet's plane inside the facet only, because the facet's edges are mesh edges). Each cavity is
//! deleted and split by the facet's plane into two half-cavities, which are refilled one after
//! the other by the modified gift-wrapping of §4.4 (after Shewchuk 2000, *Sweep algorithms for
//! constructing higher-dimensional constrained Delaunay triangulations*): repeatedly take an
//! open triangle `σ` of the half-cavity's front and join it to the vertex `w` of the half-cavity
//! such that
//! - (i) `w` is on the open side of `σ` (positive volume),
//! - (ii) the tetrahedron meets every cavity, front and constraint triangle properly
//!   ([`super::meet`]), and contains no other cavity vertex,
//! - (iii) among those, the perturbed circumsphere of `σ w` contains no other candidate: the
//!   candidate that is first in the pencil of spheres through `σ`, with the same symbolic
//!   perturbation as the Delaunay tetrahedralization, so both half-cavities induce the same
//!   triangulation on the facet.
//!
//! The first half-cavity leaves its triangles on the facet's plane open; the second closes them.
//! Constraint triangles already recovered that lie inside a cavity stay constraints: both of
//! their sides enter the front.
//!
//! Si's cavity expansion with local Delaunay tetrahedralizations (§3.4) is a speed-up the paper
//! uses before falling back to gift-wrapping; gift-wrapping alone is the always-correct path
//! (§4.4: it builds the constrained Delaunay tetrahedralization of any half-cavity that admits
//! one), so it is the one implemented here.

use std::collections::{HashMap, HashSet};

use crate::predicates::{GenericPoint, Sign, generic_orient3d};
use crate::{Error, Result};

use super::Recovery;
use super::complex::{NONE, corrupt, face_toward, sorted3};
use super::meet::{in_closed_tet, meet_properly};

/// Most rounds of finding and refilling the cavities of one facet (one suffices in exact
/// arithmetic).
const MAX_ROUNDS_PER_FACET: usize = 16;

/// A box: `(lo, hi)`.
type Bounds = ([f64; 3], [f64; 3]);

/// The open front of a half-cavity: faces, their boxes, and whether each is still open.
type Front<'a> = (&'a [[u32; 3]], &'a [Bounds], &'a [bool]);

/// An oriented triangle rotated so its smallest id comes first: the key of an oriented face.
fn canon([a, b, c]: [u32; 3]) -> [u32; 3] {
    if a <= b && a <= c {
        [a, b, c]
    } else if b <= a && b <= c {
        [b, c, a]
    } else {
        [c, a, b]
    }
}

fn boxes_meet((l1, h1): Bounds, (l2, h2): Bounds) -> bool {
    (0..3).all(|k| l1[k] <= h2[k] && l2[k] <= h1[k])
}

/// The side of the facet's plane each vertex is on, cached for one facet.
struct Plane {
    corners: [GenericPoint; 3],
    side: HashMap<u32, Sign>,
}

impl Recovery {
    fn plane(&self, f: usize) -> Result<Plane> {
        let tri = *self.triangles.get(f).ok_or_else(|| corrupt("facet id"))?;
        let corner = |v: u32| self.cx.pt(v).cloned();
        Ok(Plane { corners: [corner(tri[0])?, corner(tri[1])?, corner(tri[2])?], side: HashMap::new() })
    }

    fn facet_side(&self, plane: &mut Plane, v: u32) -> Result<Sign> {
        if let Some(&s) = plane.side.get(&v) {
            return Ok(s);
        }
        let [a, b, c] = &plane.corners;
        // tet4 volume sign of (a, b, c, v): positive on the side the facet's normal points to.
        let s = -generic_orient3d(a, b, c, self.cx.pt(v)?);
        plane.side.insert(v, s);
        Ok(s)
    }

    /// Whether edge `p q` crosses the facet's open interior.
    fn pierces(&self, plane: &mut Plane, p: u32, q: u32) -> Result<bool> {
        if self.facet_side(plane, p)? * self.facet_side(plane, q)? != Sign::Negative {
            return Ok(false);
        }
        let [a, b, c] = &plane.corners;
        let (pp, pq) = (self.cx.pt(p)?, self.cx.pt(q)?);
        let s0 = generic_orient3d(pp, pq, a, b);
        if s0.is_zero() {
            return Ok(false);
        }
        Ok(generic_orient3d(pp, pq, b, c) == s0 && generic_orient3d(pp, pq, c, a) == s0)
    }

    fn crosses(&self, plane: &mut Plane, t: u32) -> Result<bool> {
        let vs = self.cx.tet(t)?;
        for (i, j) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
            if let (Some(&p), Some(&q)) = (vs.get(i), vs.get(j))
                && self.pierces(plane, p, q)?
            {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// The vertices on facet `f`'s boundary: its corners and the Steiner points on its edges.
    pub(super) fn facet_ring(&self, f: usize) -> Result<Vec<u32>> {
        let tri = *self.triangles.get(f).ok_or_else(|| corrupt("facet id"))?;
        let mut ring = Vec::new();
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            let s = *self.edge_segment.get(&(a.min(b), a.max(b))).ok_or_else(|| corrupt("facet edge without a segment"))?;
            let chain = self.segment(s)?.chain_from(a);
            ring.extend(chain.iter().take(chain.len().saturating_sub(1)));
        }
        Ok(ring)
    }

    /// The face-connected sets of tetrahedra crossing facet `f`'s interior. Each contains a
    /// tetrahedron incident to a vertex of the facet's boundary `ring`.
    fn crossing_cavities(&mut self, plane: &mut Plane, ring: &[u32]) -> Result<Vec<Vec<u32>>> {
        let mut cavities = Vec::new();
        let mut seen: HashSet<u32> = HashSet::new();
        let mut star = Vec::new();
        for &v in ring {
            self.cx.star(v, &mut star)?;
            for &t in &star {
                if seen.contains(&t) || !self.crosses(plane, t)? {
                    continue;
                }
                seen.insert(t);
                let mut cavity = vec![t];
                let mut i = 0;
                while let Some(&c) = cavity.get(i) {
                    i += 1;
                    for n in self.cx.adj(c)? {
                        if n != NONE && !seen.contains(&n) && self.crosses(plane, n)? {
                            seen.insert(n);
                            cavity.push(n);
                        }
                    }
                }
                cavities.push(cavity);
            }
        }
        Ok(cavities)
    }

    /// Whether facet `f` is a union of mesh triangles: no tetrahedron crosses it.
    pub(super) fn facet_present(&mut self, f: usize) -> Result<bool> {
        let ring = self.facet_ring(f)?;
        let mut plane = self.plane(f)?;
        Ok(self.crossing_cavities(&mut plane, &ring)?.is_empty())
    }

    /// Recover facet `f`, then register its triangles as constraints. Filling one cavity leaves
    /// the others untouched (they share no tetrahedron), so one scan finds them all; the scan is
    /// repeated until it finds none.
    pub(super) fn recover_facet(&mut self, f: usize) -> Result<()> {
        let ring = self.facet_ring(f)?;
        let mut plane = self.plane(f)?;
        for _ in 0..MAX_ROUNDS_PER_FACET {
            let cavities = self.crossing_cavities(&mut plane, &ring)?;
            if cavities.is_empty() {
                return self.register_facet(f, &ring);
            }
            for cavity in cavities {
                self.stats.cavities += 1;
                self.fill_cavity(f, &mut plane, &cavity)?;
            }
        }
        Err(Error::Limit(format!("facet {f} was still crossed after {MAX_ROUNDS_PER_FACET} rounds of cavity filling")))
    }

    /// Record the mesh triangles that make up facet `f` (all three vertices on its boundary).
    fn register_facet(&mut self, f: usize, ring: &[u32]) -> Result<()> {
        let on: HashSet<u32> = ring.iter().copied().collect();
        let mut star = Vec::new();
        for &v in ring {
            self.cx.star(v, &mut star)?;
            for &t in &star {
                let vs = self.cx.tet(t)?;
                for k in 0..4 {
                    let face = face_toward(vs, k);
                    if face.contains(&v) && face.iter().all(|u| on.contains(u)) {
                        let key = sorted3(face);
                        match self.constraint.insert(key, f as u32) {
                            Some(g) if g != f as u32 => {
                                return Err(Error::Failed(format!("input triangles {g} and {f} overlap")));
                            }
                            _ => {}
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// Delete the cavity and refill its two halves by gift-wrapping.
    fn fill_cavity(&mut self, f: usize, plane: &mut Plane, cavity: &[u32]) -> Result<()> {
        let inside: HashSet<u32> = cavity.iter().copied().collect();
        // Per half: the boundary faces (cavity on their positive side) and the constraint faces
        // inside it, and the other sides of those constraints.
        let (mut fixed_up, mut fixed_down) = (Vec::new(), Vec::new());
        let (mut back_up, mut back_down) = (Vec::new(), Vec::new());
        let mut vertices: Vec<u32> = Vec::new();
        for &t in cavity {
            let vs = self.cx.tet(t)?;
            vertices.extend(vs);
            for (j, n) in self.cx.adj(t)?.into_iter().enumerate() {
                let face = face_toward(vs, j);
                let interior = n != NONE && inside.contains(&n);
                if interior && (n < t || !self.constraint.contains_key(&sorted3(face))) {
                    continue;
                }
                let mut pos = false;
                let mut neg = false;
                for v in face {
                    match self.facet_side(plane, v)? {
                        Sign::Positive => pos = true,
                        Sign::Negative => neg = true,
                        Sign::Zero => {}
                    }
                }
                let (fixed, back) = match (pos, neg) {
                    (true, false) => (&mut fixed_up, &mut back_up),
                    (false, true) => (&mut fixed_down, &mut back_down),
                    _ => {
                        let at = self.cx.approx(face[0])?;
                        return Err(Error::Measured(
                            crate::Kind::Failed,
                            crate::message![
                                format!("facet {f}: a cavity face near "),
                                crate::message::point(at),
                                " is not on one side of the facet's plane"
                            ],
                        ));
                    }
                };
                fixed.push(face);
                if interior {
                    // A recovered constraint inside the cavity: both sides are open.
                    let [a, b, c] = face;
                    back.push([a, c, b]);
                }
            }
        }
        vertices.sort_unstable();
        vertices.dedup();
        let mut sides = Vec::with_capacity(vertices.len());
        for &v in &vertices {
            sides.push(self.facet_side(plane, v)?);
        }
        let above: Vec<u32> = vertices.iter().zip(&sides).filter(|(_, s)| **s != Sign::Negative).map(|(v, _)| *v).collect();
        let below: Vec<u32> = vertices.iter().zip(&sides).filter(|(_, s)| **s != Sign::Positive).map(|(v, _)| *v).collect();
        let on_plane = |v: u32| plane.side.get(&v) == Some(&Sign::Zero);

        let front_up = fixed_up.iter().chain(&back_up).copied().collect();
        let (mut new, open) = self.gift_wrap(f, front_up, &fixed_up, &above, Some(&on_plane))?;
        if open.iter().any(|face| !face.iter().all(|&v| on_plane(v))) {
            return Err(Error::Failed(format!("facet {f}: the upper half-cavity was left open off the facet")));
        }
        // The facet's triangles made by the upper half bound the lower one.
        fixed_down.extend(open.iter().copied());
        let front_down = fixed_down.iter().chain(&back_down).copied().collect();
        let (lower, rest) = self.gift_wrap(f, front_down, &fixed_down, &below, None)?;
        if !rest.is_empty() {
            return Err(Error::Failed(format!("facet {f}: the lower half-cavity did not close")));
        }
        new.extend(lower);
        let used: HashSet<u32> = new.iter().flatten().copied().collect();
        if let Some(&lost) = vertices.iter().find(|v| !used.contains(v)) {
            let at = self.cx.approx(lost)?;
            return Err(Error::Measured(
                crate::Kind::Failed,
                crate::message![format!("facet {f}: refilling a cavity lost vertex {lost} at "), crate::message::point(at)],
            ));
        }
        self.stats.gift_wrapped_tets += new.len();
        self.cx.replace(cavity, &new)?;
        Ok(())
    }

    /// Fill a half-cavity from its open `front` faces (each with the unfilled region on its
    /// positive side), using only `vertices`. `fixed` are the triangles bounding it, which no new
    /// tetrahedron may cross. Faces on the facet's plane are left open when `skip` says so. The new
    /// tetrahedra and the faces left open.
    fn gift_wrap(
        &mut self,
        f: usize,
        front: Vec<[u32; 3]>,
        fixed: &[[u32; 3]],
        vertices: &[u32],
        skip: Option<&dyn Fn(u32) -> bool>,
    ) -> Result<(Vec<[u32; 4]>, Vec<[u32; 3]>)> {
        let mut faces: Vec<[u32; 3]> = Vec::with_capacity(front.len() * 2);
        let mut boxes: Vec<([f64; 3], [f64; 3])> = Vec::with_capacity(front.len() * 2);
        let mut alive: Vec<bool> = Vec::new();
        let mut index: HashMap<[u32; 3], usize> = HashMap::new();
        let mut queue: Vec<usize> = Vec::new();
        for face in front {
            if index.insert(canon(face), faces.len()).is_some() {
                return Err(corrupt("a cavity face is listed twice"));
            }
            queue.push(faces.len());
            faces.push(face);
            boxes.push(self.cx.bounds_of(&face));
            alive.push(true);
        }
        queue.reverse();
        let fixed_boxes: Vec<([u32; 3], Bounds)> = fixed.iter().map(|&t| (t, self.cx.bounds_of(&t))).collect();
        let cap = 64 * (faces.len() + vertices.len()) + 64;
        let mut tets: Vec<[u32; 4]> = Vec::new();
        while let Some(i) = queue.pop() {
            if !alive.get(i).copied().unwrap_or(false) {
                continue;
            }
            let sigma = *faces.get(i).ok_or_else(|| corrupt("front index"))?;
            if let Some(on) = skip
                && sigma.iter().all(|&v| on(v))
            {
                continue;
            }
            let w = self.apex(f, sigma, (&faces, &boxes, &alive), &fixed_boxes, vertices)?;
            let tet = [sigma[0], sigma[1], sigma[2], w];
            tets.push(tet);
            if tets.len() > cap {
                return Err(Error::Limit(format!("facet {f}: gift-wrapping a cavity made more than {cap} tetrahedra")));
            }
            if let Some(a) = alive.get_mut(i) {
                *a = false;
            }
            index.remove(&canon(sigma));
            for k in 0..3 {
                let g = face_toward(tet, k);
                if let Some(j) = index.remove(&canon(g)) {
                    if let Some(a) = alive.get_mut(j) {
                        *a = false;
                    }
                    continue;
                }
                let r = [g[0], g[2], g[1]];
                if index.insert(canon(r), faces.len()).is_some() {
                    let at = self.cx.approx(r[0])?;
                    return Err(Error::Measured(
                        crate::Kind::Failed,
                        crate::message![format!("facet {f}: gift-wrapping overlapped itself near "), crate::message::point(at)],
                    ));
                }
                queue.push(faces.len());
                faces.push(r);
                boxes.push(self.cx.bounds_of(&r));
                alive.push(true);
            }
        }
        let open: Vec<[u32; 3]> = faces.iter().zip(&alive).filter(|(_, a)| **a).map(|(f, _)| *f).collect();
        Ok((tets, open))
    }

    /// The apex for front face `sigma`: conditions (i)–(iii) of the module docs.
    fn apex(&mut self, f: usize, sigma: [u32; 3], front: Front<'_>, fixed: &[([u32; 3], Bounds)], vertices: &[u32]) -> Result<u32> {
        let mut candidates = Vec::new();
        for &w in vertices {
            if !sigma.contains(&w) && self.cx.volume([sigma[0], sigma[1], sigma[2], w])? == Sign::Positive {
                candidates.push(w);
            }
        }
        let mut excluded = vec![false; candidates.len()];
        for _ in 0..candidates.len() {
            let mut best: Option<(usize, u32)> = None;
            for (k, &w) in candidates.iter().enumerate() {
                if excluded.get(k).copied().unwrap_or(true) {
                    continue;
                }
                best = match best {
                    None => Some((k, w)),
                    Some((bk, b)) => {
                        if self.cx.in_sphere([sigma[0], sigma[1], sigma[2], b], w)? {
                            Some((k, w))
                        } else {
                            Some((bk, b))
                        }
                    }
                };
            }
            let Some((k, w)) = best else { break };
            if self.valid(sigma, w, front, fixed, vertices)? {
                return Ok(w);
            }
            if let Some(e) = excluded.get_mut(k) {
                *e = true;
            }
        }
        let at = self.cx.approx(sigma[0])?;
        Err(Error::Measured(
            crate::Kind::Failed,
            crate::message![format!("facet {f}: gift-wrapping found no valid apex for a cavity face near "), crate::message::point(at)],
        ))
    }

    /// Condition (ii): the tetrahedron `sigma w` meets every front and fixed triangle properly
    /// and contains no other cavity vertex.
    fn valid(&self, sigma: [u32; 3], w: u32, (faces, boxes, alive): Front<'_>, fixed: &[([u32; 3], Bounds)], vertices: &[u32]) -> Result<bool> {
        let tet = [sigma[0], sigma[1], sigma[2], w];
        let bx = self.cx.bounds_of(&tet);
        let points = self.cx.points();
        for &u in vertices {
            if !tet.contains(&u) && boxes_meet(bx, self.cx.bounds(u)) && in_closed_tet(points, tet, u)? {
                return Ok(false);
            }
        }
        for ((face, fb), _) in faces.iter().zip(boxes).zip(alive).filter(|(_, a)| **a) {
            if boxes_meet(bx, *fb) && !meet_properly(points, &tet, face)? {
                return Ok(false);
            }
        }
        for (face, fb) in fixed {
            if boxes_meet(bx, *fb) && !meet_properly(points, &tet, face)? {
                return Ok(false);
            }
        }
        Ok(true)
    }
}
