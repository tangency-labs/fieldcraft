//! Virtual topology: a face split in two along its curve through a plane, for meshing only.
//!
//! A partition cuts whole bodies into cells with the kernel; a [`FaceSplit`] splits one face of
//! a body along the curve where a plane crosses it, in what the meshers take from the body (the
//! B-rep is unchanged): the face's triangles are cut along the plane, the pieces on the plane's
//! negative side are `split:<name>:a` and the others `split:<name>:b` (the partition's
//! convention), both on the face's surface, and the cut is a new edge between them, so the
//! meshers put a row of nodes along it. Edges of the face that the plane crosses are cut there
//! too. Splits apply after the merges ([`crate::virtual_faces`]), by persistent names, whenever
//! the design is opened; a merge cannot name a split face.

use std::collections::BTreeMap;

use crate::edges::EDGE_TOLERANCE;
use crate::{Design, Error, Face, FeatureError, Result, Shape, SurfaceMesh};

/// The largest coordinate a split plane may have.
const MAX_COORD: f64 = 1e9; // unit-factors: a cap on hostile input
/// A split plane's normal must be longer than this.
const MIN_NORMAL: f64 = 1e-12; // unit-factors: the length of a direction vector
/// A point within this of the plane, relative to the coordinates' size, is on it.
const ON_PLANE: f64 = 1e-12; // unit-factors: relative to the coordinates' size

/// Most splits a design takes.
pub const MAX_SPLITS: usize = 1024;

/// A face of a body split by a plane.
#[derive(Debug, Clone, PartialEq)]
pub struct FaceSplit {
    /// Names the two faces (`split:<name>:a|b`).
    pub name: String,
    pub body: String,
    pub face: String,
    /// A point on the plane and its normal (any length but zero).
    pub point: [f64; 3],
    pub normal: [f64; 3],
}

/// One split as a body applies it: the face it cuts (kept as the `a` side), the `b` face's
/// index, and the plane (unit normal).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Cut {
    a: u32,
    b: u32,
    point: [f64; 3],
    normal: [f64; 3],
}

impl Cut {
    pub(crate) fn faces(&self) -> (u32, u32) {
        (self.a, self.b)
    }
}

fn signed(c: &Cut, p: [f64; 3]) -> f64 {
    (0..3).map(|k| (p[k] - c.point[k]) * c.normal[k]).sum()
}

impl Design {
    /// The design with `splits` applied in order (after its partitions and merges). A split that
    /// fails (a body or face that is not there, a plane that does not cross the face) is listed
    /// in [`Design::errors`] (feature `split <name>`) and changes nothing.
    pub fn split(mut self, splits: &[FaceSplit]) -> Result<Design> {
        if splits.len() > MAX_SPLITS {
            return Err(Error::Invalid(format!("a design takes at most {MAX_SPLITS} face splits")));
        }
        for sp in splits {
            let feature = format!("split {}", sp.name);
            let Some(body) = self.bodies.iter_mut().find(|b| b.name == sp.body) else {
                self.errors.push(FeatureError { feature, message: format!("no body `{}`", sp.body), bodies: vec![sp.body.clone()] });
                continue;
            };
            match body.split_face(&sp.name, &sp.face, sp.point, sp.normal) {
                Ok(split) => *body = split,
                Err(e) => self.errors.push(FeatureError { feature, message: e.to_string(), bodies: vec![sp.body.clone()] }),
            }
        }
        if self.bodies.iter().any(|b| !b.cuts.is_empty()) {
            let opened = self.source.iter().map(|s| self.bodies.iter().find(|b| b.name == s.name).unwrap_or(s).clone()).collect();
            self.opened = std::sync::Arc::new(opened);
        }
        Ok(self)
    }
}

impl Shape {
    /// This shape with face `face` split along the plane through `point` normal to `normal`:
    /// `split:<name>:a` on its negative side, `split:<name>:b` on the other.
    pub fn split_face(&self, name: &str, face: &str, point: [f64; 3], normal: [f64; 3]) -> Result<Shape> {
        let len = normal.iter().map(|v| v * v).sum::<f64>().sqrt();
        if !(point.iter().chain(&normal).all(|v| v.is_finite() && v.abs() <= MAX_COORD) && len > MIN_NORMAL) {
            return Err(Error::Invalid("a split needs a plane with a finite point and a non-zero normal".into()));
        }
        let named = self.faces_named(face);
        let [a] = named.as_slice() else {
            return Err(Error::Invalid(format!("body `{}` has no face `{face}` (or several pieces of it: name one)", self.name)));
        };
        let (an, bn) = (format!("split:{name}:a"), format!("split:{name}:b"));
        if self.faces.iter().any(|f| f.name == an || f.name == bn) {
            return Err(Error::Invalid(format!("body `{}` already has the faces of a split `{name}`", self.name)));
        }
        let b = u32::try_from(self.faces.len()).map_err(|_| Error::Invalid("too many faces".into()))?;
        let mut shape = self.clone();
        let original = shape.faces.get(*a as usize).cloned().ok_or_else(|| Error::Invalid("no such face".into()))?;
        let members = if original.members.is_empty() { vec![original.name.clone()] } else { original.members.clone() };
        if let Some(f) = shape.faces.get_mut(*a as usize) {
            *f = Face { name: an, exact: original.exact, members: members.clone() };
        }
        shape.faces.push(Face { name: bn, exact: original.exact, members });
        shape.cuts.push(Cut { a: *a, b, point, normal: normal.map(|v| v / len) });
        // The plane must cross the face: both pieces have triangles.
        let mesh = shape.tessellate((self.size()? * EDGE_TOLERANCE).max(1e-4))?;
        for side in [*a, b] {
            if !mesh.tri_face.contains(&side) {
                return Err(Error::Invalid(format!("the plane does not cross face `{face}` of body `{}`", self.name)));
            }
        }
        Ok(shape)
    }
}

/// `mesh` with each of `cuts` applied in order (see the module docs).
pub(crate) fn apply(mesh: &mut SurfaceMesh, cuts: &[Cut]) {
    for c in cuts {
        cut(mesh, c);
    }
}

fn lerp(p: [f64; 3], q: [f64; 3], t: f64) -> [f64; 3] {
    [p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1]), p[2] + t * (q[2] - p[2])]
}

fn cut(mesh: &mut SurfaceMesh, c: &Cut) {
    let scale = mesh.positions.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs())).max(1.0);
    let eps = ON_PLANE * scale;
    let side = |d: f64| {
        if d > eps {
            1
        } else if d < -eps {
            -1
        } else {
            0
        }
    };
    // A new point on the side (u, v) where the plane crosses it, once per side.
    let mut made: BTreeMap<(u32, u32), u32> = BTreeMap::new();
    let mut segments: Vec<(u32, u32)> = Vec::new();
    let n = mesh.triangles.len();
    for t in 0..n {
        let (Some(&tri), Some(&tag)) = (mesh.triangles.get(t), mesh.tri_face.get(t)) else { continue };
        if tag != c.a {
            continue;
        }
        let pos = |i: u32, m: &SurfaceMesh| m.positions.get(i as usize).copied().unwrap_or([0.0; 3]);
        let s = tri.map(|i| side(signed(c, pos(i, mesh))));
        let tag_of = |sg: i32| if sg > 0 { c.b } else { c.a };
        if s.iter().all(|&x| x >= 0) || s.iter().all(|&x| x <= 0) {
            let sg = s.iter().copied().find(|&x| x != 0).unwrap_or(-1);
            if let Some(f) = mesh.tri_face.get_mut(t) {
                *f = tag_of(sg);
            }
            // A side on the plane between pieces of both sides is the cut too (picked up when
            // chaining from the other triangles' points; left to the edges otherwise).
            continue;
        }
        let mut point_on = |u: u32, v: u32, m: &mut SurfaceMesh| -> u32 {
            let key = (u.min(v), u.max(v));
            if let Some(&i) = made.get(&key) {
                return i;
            }
            let (p, q) = (pos(key.0, m), pos(key.1, m));
            let (dp, dq) = (signed(c, p), signed(c, q));
            let x = lerp(p, q, dp / (dp - dq));
            let i = u32::try_from(m.positions.len()).unwrap_or(u32::MAX);
            m.positions.push(x);
            made.insert(key, i);
            i
        };
        let rot = |k: usize| [tri[k], tri[(k + 1) % 3], tri[(k + 2) % 3]];
        let rs = |k: usize| [s[k], s[(k + 1) % 3], s[(k + 2) % 3]];
        let mut out: Vec<([u32; 3], i32)> = Vec::with_capacity(3);
        if let Some(k) = (0..3).find(|&k| s[k] == 0) {
            // A corner on the plane, the other two on either side.
            let ([z, u, v], [_, su, sv]) = (rot(k), rs(k));
            let m = point_on(u, v, mesh);
            out.push(([z, u, m], su));
            out.push(([z, m, v], sv));
            segments.push((z, m));
        } else {
            // One corner alone on its side.
            let k = (0..3).find(|&k| s[k] != s[(k + 1) % 3] && s[k] != s[(k + 2) % 3]).unwrap_or(0);
            let ([l, u, v], [sl, su, _]) = (rot(k), rs(k));
            let (pu, pv) = (point_on(l, u, mesh), point_on(l, v, mesh));
            out.push(([l, pu, pv], sl));
            out.push(([pu, u, v], su));
            out.push(([pu, v, pv], su));
            segments.push((pu, pv));
        }
        let mut first = true;
        for (tri, sg) in out {
            if first {
                if let (Some(slot), Some(f)) = (mesh.triangles.get_mut(t), mesh.tri_face.get_mut(t)) {
                    *slot = tri;
                    *f = tag_of(sg);
                }
                first = false;
            } else {
                mesh.triangles.push(tri);
                mesh.tri_face.push(tag_of(sg));
            }
        }
    }
    // The face's edges: cut where the plane crosses them, each piece on its side's face.
    let mut edges = Vec::with_capacity(mesh.edges.len() + 1);
    let mut edge_faces = Vec::with_capacity(mesh.edges.len() + 1);
    for (points, faces) in std::mem::take(&mut mesh.edges).into_iter().zip(std::mem::take(&mut mesh.edge_faces)) {
        if !faces.contains(&c.a) {
            edges.push(points);
            edge_faces.push(faces);
            continue;
        }
        let mut piece: Vec<[f64; 3]> = Vec::new();
        let mut piece_side = 0;
        let mut flush = |piece: &mut Vec<[f64; 3]>, sg: i32| {
            if piece.len() >= 2 {
                let on = if sg > 0 { c.b } else { c.a };
                edges.push(std::mem::take(piece));
                edge_faces.push(faces.iter().map(|&f| if f == c.a { on } else { f }).collect());
            }
            piece.clear();
        };
        for (k, &p) in points.iter().enumerate() {
            let sp = side(signed(c, p));
            if k > 0
                && let Some(&q) = points.get(k - 1)
            {
                let sq = side(signed(c, q));
                if sp != 0 && sq != 0 && sp != sq {
                    let (dq, dp) = (signed(c, q), signed(c, p));
                    let x = lerp(q, p, dq / (dq - dp));
                    piece.push(x);
                    flush(&mut piece, piece_side);
                    piece.push(x);
                }
            }
            if sp != 0 {
                piece_side = sp;
            }
            piece.push(p);
        }
        flush(&mut piece, piece_side);
    }
    // The cut: its segments chained into polylines, between the two faces. Copies of one point
    // (a tessellation may repeat points) are one point for the chaining.
    let close = (2.0 * mesh.weld_distance).max(eps);
    let mut canon: BTreeMap<u32, u32> = BTreeMap::new();
    let mut reps: Vec<(u32, [f64; 3])> = Vec::new();
    for &i in segments.iter().flat_map(|(u, v)| [u, v]) {
        if canon.contains_key(&i) {
            continue;
        }
        let p = mesh.positions.get(i as usize).copied().unwrap_or([0.0; 3]);
        let r = match reps.iter().find(|(_, q)| (0..3).map(|k| (p[k] - q[k]).powi(2)).sum::<f64>().sqrt() <= close) {
            Some((r, _)) => *r,
            None => {
                reps.push((i, p));
                i
            }
        };
        canon.insert(i, r);
    }
    let segments: Vec<(u32, u32)> =
        segments.iter().map(|(u, v)| (canon.get(u).copied().unwrap_or(*u), canon.get(v).copied().unwrap_or(*v))).filter(|(u, v)| u != v).collect();
    for line in chain(&segments) {
        let points: Vec<[f64; 3]> = line.iter().filter_map(|&i| mesh.positions.get(i as usize).copied()).collect();
        if points.len() >= 2 {
            edges.push(points);
            edge_faces.push(vec![c.a, c.b]);
        }
    }
    mesh.edges = edges;
    mesh.edge_faces = edge_faces;
}

/// Segments (point index pairs) chained end to end into polylines, in a fixed order.
fn chain(segments: &[(u32, u32)]) -> Vec<Vec<u32>> {
    let mut at: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for (i, &(u, v)) in segments.iter().enumerate() {
        at.entry(u).or_default().push(i);
        at.entry(v).or_default().push(i);
    }
    let mut used = vec![false; segments.len()];
    let mut lines = Vec::new();
    // Start at ends (points on one segment) first, then what is left (closed loops).
    let ends: Vec<usize> = at.values().filter(|s| s.len() == 1).flatten().copied().collect();
    let starts: Vec<usize> = ends.into_iter().chain(0..segments.len()).collect();
    for s in starts {
        if used.get(s).copied().unwrap_or(true) {
            continue;
        }
        let Some(&(u, v)) = segments.get(s) else { continue };
        // From an end segment, start at its free point.
        let (mut line, mut last) = if at.get(&v).is_some_and(|x| x.len() == 1) { (vec![v, u], u) } else { (vec![u, v], v) };
        if let Some(x) = used.get_mut(s) {
            *x = true;
        }
        while let Some(&next) = at.get(&last).and_then(|x| x.iter().find(|&&i| !used.get(i).copied().unwrap_or(true))) {
            if let Some(x) = used.get_mut(next) {
                *x = true;
            }
            let Some(&(a, b)) = segments.get(next) else { break };
            last = if a == last { b } else { a };
            line.push(last);
        }
        lines.push(line);
    }
    lines
}

#[cfg(test)]
mod tests {
    use crate::open_design;
    use crate::sample::{self, Plate};

    use super::*;

    /// The plate's top (z = 10, 40 × 30 with the Ø10 hole at (15, 15)) split at x = 30: two
    /// faces on the top's plane, a straight cut edge between them along x = 30, the top's rim
    /// edges at y = 0 and y = 30 cut there; triangles all on one side of the plane each.
    #[test]
    fn the_plates_top_splits_along_a_plane() {
        let d = open_design(&sample::plate_design(&Plate::default()).unwrap()).unwrap();
        let body = d.bodies[0].name.clone();
        let b = &d.bodies[0];
        let surfaces = b.face_surfaces().unwrap();
        let top = b
            .faces()
            .iter()
            .zip(&surfaces)
            .find(|(_, s)| matches!(s, Some(crate::FaceSurface::Plane { normal, point }) if normal[2] > 0.5 && point[2] > 9.0))
            .unwrap()
            .0
            .name
            .clone();
        let split = |x: f64| FaceSplit { name: "S".into(), body: body.clone(), face: top.clone(), point: [x, 0.0, 0.0], normal: [1.0, 0.0, 0.0] };
        let d = d.split(&[split(30.0)]).unwrap();
        assert!(d.errors.is_empty(), "{:?}", d.errors);
        let b = &d.bodies[0];
        assert_eq!(b.faces().len(), 9);
        let (fa, fb) = (b.face_index("split:S:a").unwrap(), b.face_index("split:S:b").unwrap());
        let s = b.face_surfaces().unwrap();
        assert_eq!(s[fa as usize], s[fb as usize]);
        let mesh = b.tessellate(0.01).unwrap();
        for (t, &f) in mesh.triangles.iter().zip(&mesh.tri_face) {
            let x = t.map(|i| mesh.positions[i as usize][0]);
            if f == fa {
                assert!(x.iter().all(|&v| v <= 30.0 + 1e-9), "{x:?}");
            } else if f == fb {
                assert!(x.iter().all(|&v| v >= 30.0 - 1e-9), "{x:?}");
            }
        }
        // The area is kept: 40 × 30 less the hole, split 30 × 30 − π25 and 10 × 30.
        let area = |f: u32| -> f64 {
            mesh.triangles
                .iter()
                .zip(&mesh.tri_face)
                .filter(|(_, g)| **g == f)
                .map(|(t, _)| {
                    let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
                    let (u, v) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
                    0.5 * ((u[1] * v[2] - u[2] * v[1]).powi(2) + (u[2] * v[0] - u[0] * v[2]).powi(2) + (u[0] * v[1] - u[1] * v[0]).powi(2)).sqrt()
                })
                .sum()
        };
        assert!((area(fb) - 300.0).abs() < 1e-6, "{}", area(fb));
        assert!((area(fa) - (900.0 - std::f64::consts::PI * 25.0)).abs() < 0.5, "{}", area(fa));
        let edges = b.edges().unwrap();
        let cut: Vec<_> = edges.iter().filter(|e| e.faces.contains(&fa) && e.faces.contains(&fb)).collect();
        assert_eq!(cut.len(), 1, "{:?}", edges.iter().map(|e| &e.name).collect::<Vec<_>>());
        assert!((cut[0].length() - 30.0).abs() < 1e-9 && cut[0].points.iter().all(|p| (p[0] - 30.0).abs() < 1e-9));
        // The top's rim edges along x are cut at x = 30: four pieces where there were two.
        let rim = edges.iter().filter(|e| e.faces.iter().any(|&f| f == fa || f == fb) && !(e.faces.contains(&fa) && e.faces.contains(&fb))).count();
        assert_eq!(rim, 2 + 2 + 2 + 2, "four sides, the two along x in two pieces, and the hole's two rim arcs");
    }

    #[test]
    fn bad_splits_are_listed() {
        let d = open_design(&sample::plate_design(&Plate::default()).unwrap()).unwrap();
        let body = d.bodies[0].name.clone();
        let face = d.bodies[0].faces()[0].name.clone();
        let sp = |n: &str, b: &str, f: &str, x: f64, normal: [f64; 3]| FaceSplit {
            name: n.into(),
            body: b.into(),
            face: f.into(),
            point: [x, 0.0, 0.0],
            normal,
        };
        let d = d
            .split(&[
                sp("Far", &body, &face, 1000.0, [1.0, 0.0, 0.0]),
                sp("NoFace", &body, "nope", 20.0, [1.0, 0.0, 0.0]),
                sp("NoBody", "x", &face, 20.0, [1.0, 0.0, 0.0]),
                sp("Zero", &body, &face, 20.0, [0.0; 3]),
            ])
            .unwrap();
        assert_eq!(d.errors.len(), 4, "{:?}", d.errors);
        assert_eq!(d.bodies[0].faces().len(), 8);
    }
}
