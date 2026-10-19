//! Independent checks of a triangulation: they re-derive every property from the arrays with
//! the exact predicates; nothing here relies on how the insertion built them.

use super::{DEAD, GHOST, NONE, Triangulation, volume_sign};
use crate::predicates::{Point3, Sign, insphere, insphere_sos};
use crate::quality::shape::{det, sub};
use crate::{Error, Result};

fn failed<T>(what: String) -> Result<T> {
    Err(Error::Failed(format!("Delaunay verification: {what}")))
}

/// `det[b − a, c − a, d − a] / 6` in floating point (a measure for the volume balance only).
pub(super) fn volume(p: [&Point3; 4]) -> f64 {
    let [a, b, c, d] = p.map(Point3::coords);
    det(sub(b, a), sub(c, a), sub(d, a)) / 6.0
}

/// Whether `e` is outside the circumsphere of the positive-volume tetrahedron `[a, b, c, d]`:
/// never inside exactly, and outside with the perturbation.
fn outside_sphere([a, b, c, d]: [&Point3; 4], e: &Point3) -> bool {
    // `insphere` wants Shewchuk-positive orientation: swap the first two vertices.
    match insphere(b, a, c, d, e) {
        Sign::Positive => false,
        Sign::Negative => true,
        Sign::Zero => insphere_sos(b, a, c, d, e) == Sign::Negative,
    }
}

/// The face of `vs` opposite slot `k`, sorted.
pub(super) fn face(vs: [u32; 4], k: usize) -> [u32; 3] {
    let mut f = [NONE; 3];
    for (o, (_, &v)) in f.iter_mut().zip(vs.iter().enumerate().filter(|&(i, _)| i != k)) {
        *o = v;
    }
    f.sort_unstable();
    f
}

impl Triangulation {
    fn point_of(&self, v: u32) -> Result<&Point3> {
        match self.points.get(v as usize) {
            Some(p) => Ok(p),
            None => failed(format!("vertex {v} does not exist")),
        }
    }

    fn points_of(&self, [a, b, c, d]: [u32; 4]) -> Result<[&Point3; 4]> {
        Ok([self.point_of(a)?, self.point_of(b)?, self.point_of(c)?, self.point_of(d)?])
    }

    fn points_with<'a>(&'a self, vs: [u32; 4], k: usize, p: &'a Point3) -> Result<[&'a Point3; 4]> {
        let mut out = [p; 4];
        for (i, o) in out.iter_mut().enumerate() {
            if i != k {
                *o = self.point_of(vs.get(i).copied().unwrap_or(NONE))?;
            }
        }
        Ok(out)
    }

    /// The vertex of `o` across its shared face with `t`, and the slot of that face in `o`.
    fn across(&self, t: u32, o: u32) -> Result<([u32; 4], usize)> {
        let Some(ws) = self.tet(o) else { return failed(format!("tetrahedron {t} has a deleted neighbour {o}")) };
        let adj = self.adj.get(o as usize).copied().unwrap_or([NONE; 4]);
        match (adj.iter().position(|&n| n == t), adj.iter().rposition(|&n| n == t)) {
            (Some(k), Some(l)) if k == l => Ok((ws, k)),
            _ => failed(format!("tetrahedra {t} and {o} are not mutual neighbours exactly once")),
        }
    }

    /// Check the triangulation: adjacency is mutual and matches shared faces, every finite
    /// tetrahedron has positive volume, the hull is convex, every interior face is locally
    /// Delaunay (the opposite vertex is not inside the circumsphere, exactly, and is outside
    /// with the perturbation), and every vertex is in a tetrahedron. As a floating-point sanity
    /// check against a hull covered twice, the finite volume must match the volume the hull faces
    /// enclose to 1e-9. By Delaunay's lemma, local Delaunay faces in a
    /// triangulation of the hull make every circumsphere empty. `O(tetrahedra)`.
    pub fn verify(&self) -> Result<()> {
        let n = self.points.len();
        let mut seen = vec![false; n];
        let (mut finite_volume, mut hull_volume, mut scale) = (0.0f64, 0.0f64, 0.0f64);
        let reference = self.tets.iter().flatten().copied().find(|&v| v < DEAD);
        for (t, &vs) in self.tets.iter().enumerate() {
            if vs[0] == DEAD {
                continue;
            }
            let t = t as u32;
            let adj = self.adj.get(t as usize).copied().unwrap_or([NONE; 4]);
            for (i, &v) in vs.iter().enumerate() {
                if vs.iter().skip(i + 1).any(|&w| w == v) {
                    return failed(format!("tetrahedron {t} repeats vertex {v}"));
                }
                if v == GHOST {
                    continue;
                }
                if self.duplicate_of(v).is_some() {
                    return failed(format!("duplicate point {v} is a vertex"));
                }
                match seen.get_mut(v as usize) {
                    Some(s) => *s = true,
                    None => return failed(format!("tetrahedron {t} has vertex {v}, which does not exist")),
                }
            }
            // Per face: the neighbour's vertices and the slot of the shared face in it.
            let mut across = [([NONE; 4], 0); 4];
            for (j, (&o, slot)) in adj.iter().zip(&mut across).enumerate() {
                if o == t {
                    return failed(format!("tetrahedron {t} is its own neighbour"));
                }
                let (ws, k) = self.across(t, o)?;
                if face(vs, j) != face(ws, k) {
                    return failed(format!("tetrahedra {t} and {o} are neighbours without sharing a face"));
                }
                *slot = (ws, k);
            }
            let apex = |j: usize| across.get(j).and_then(|&(ws, k)| ws.get(k).copied()).unwrap_or(NONE);
            match vs.iter().position(|&v| v == GHOST) {
                None => {
                    let p = self.points_of(vs)?;
                    if volume_sign(p) != Sign::Positive {
                        return failed(format!("tetrahedron {t} is flat or inverted"));
                    }
                    let vol = volume(p);
                    finite_volume += vol;
                    scale += vol.abs();
                    // The test is symmetric in the two tetrahedra: check each face once.
                    for (j, &o) in adj.iter().enumerate() {
                        let u = apex(j);
                        if u != GHOST && t < o && !outside_sphere(p, self.point_of(u)?) {
                            return failed(format!("vertex {u} is inside the circumsphere of tetrahedron {t}"));
                        }
                    }
                }
                Some(g) => {
                    let fs = across.get(g).map(|&(ws, _)| ws).unwrap_or([NONE; 4]);
                    if fs.contains(&GHOST) {
                        return failed(format!("hull face of {t} borders another ghost"));
                    }
                    let inner = self.point_of(apex(g))?;
                    if volume_sign(self.points_with(vs, g, inner)?) != Sign::Negative {
                        return failed(format!("ghost tetrahedron {t} faces into the hull"));
                    }
                    for j in (0..4).filter(|&j| j != g) {
                        let apex = self.point_of(apex(j))?;
                        match volume_sign(self.points_with(vs, g, apex)?) {
                            Sign::Positive => return failed(format!("the hull is not convex at ghost tetrahedron {t}")),
                            Sign::Negative => {}
                            Sign::Zero => {
                                // Coplanar hull faces: the perturbed decision of the finite
                                // tetrahedron under this face.
                                if !outside_sphere(self.points_of(fs)?, apex) {
                                    return failed(format!("coplanar hull faces at ghost tetrahedron {t} are not perturbed-Delaunay"));
                                }
                            }
                        }
                    }
                    if let Some(r) = reference {
                        hull_volume -= volume(self.points_with(vs, g, self.point_of(r)?)?);
                    }
                }
            }
        }
        if self.dimension() == 3 {
            for (v, s) in seen.iter().enumerate() {
                let v = v as u32;
                let triangulated = self.duplicate_of(v).is_none();
                if *s != triangulated {
                    return failed(format!("vertex {v} is {}in a tetrahedron", if *s { "" } else { "not " }));
                }
                if triangulated && !self.incident_tet(v).and_then(|t| self.tet(t)).is_some_and(|vs| vs.contains(&v)) {
                    return failed(format!("vertex {v} has no valid incident tetrahedron"));
                }
            }
        } else if self.tets.iter().any(|vs| vs[0] != DEAD) {
            return failed("a lower-dimensional triangulation has tetrahedra".into());
        }
        if (finite_volume - hull_volume).abs() > 1e-9 * scale {
            return Err(Error::failed(crate::message![
                "Delaunay verification: the tetrahedra fill ",
                crate::message::volume(finite_volume),
                ", the hull encloses ",
                crate::message::volume(hull_volume)
            ]));
        }
        Ok(())
    }

    /// Brute force: no triangulated vertex is inside any circumsphere (exactly) or the perturbed
    /// circumsphere, and none is beyond a hull face. `O(points × tetrahedra)`; for tests.
    pub fn verify_empty_spheres(&self) -> Result<()> {
        let vertices: Vec<u32> = (0..self.points.len() as u32).filter(|&v| self.duplicate_of(v).is_none()).collect();
        for (t, &vs) in self.tets.iter().enumerate() {
            if vs[0] == DEAD {
                continue;
            }
            for &u in &vertices {
                if vs.contains(&u) {
                    continue;
                }
                let e = self.point_of(u)?;
                match vs.iter().position(|&v| v == GHOST) {
                    None => {
                        if !outside_sphere(self.points_of(vs)?, e) {
                            return failed(format!("vertex {u} is inside the circumsphere of tetrahedron {t}"));
                        }
                    }
                    Some(g) => {
                        if volume_sign(self.points_with(vs, g, e)?) == Sign::Positive {
                            return failed(format!("vertex {u} is beyond the hull face of ghost tetrahedron {t}"));
                        }
                    }
                }
            }
        }
        Ok(())
    }
}
