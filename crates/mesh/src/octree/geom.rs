//! The boundary as the hex mesher sees it: a triangle tree for nearest-point, ray and box
//! queries, the exact inside test, and the feature graph (corners, curves, patches).
//!
//! # Feature graph
//! Every feature edge of the [`BoundarySurface`] (between two B-rep faces, or a crease inside a
//! face) whose dihedral angle reaches the sharp angle is kept; smooth face boundaries (a fillet
//! meeting its faces tangentially) are not features, the faces on either side share a patch.
//! Neither are blades and slits, edges whose normals turn by more than the blade angle (a 3°
//! wedge, a sphere touching a plane): no valid hexahedra can follow them, so the mesh rounds
//! them within the boundary tolerance. *Corners* are vertices where the number of feature edges is not two, or where a
//! curve turns by more than [`CORNER_TURN`]; *curves* are the chains of feature edges between
//! corners, or closed loops; *patches* are the connected sets of triangles between feature edges,
//! so every patch lies in one B-rep face and is bounded by curves. This is the partition the hex
//! mesh's boundary must reproduce (Gao, Shen, Panozzo 2019, §3.3).

use std::collections::BTreeMap;

pub(crate) use crate::hex::geom::closest_on_triangle;
use crate::predicates::{Axis, Point3, Sign, orient2d_projected, orient3d, snap_to_domain};
use crate::surface::BoundarySurface;
use crate::vec3::{P, cross, dot, norm, sub};
use crate::{Error, Result};

/// A box is padded by this, relative to its size.
const PAD: f64 = 1e-9; // unit-factors: relative to the box's size
/// Round-off allowed, relative to the numbers compared (or of a barycentric coordinate).
const ROUND_OFF: f64 = 1e-12; // unit-factors: relative to the terms compared

/// A degree-2 feature vertex whose curve turns by more than this (degrees) is a corner.
pub(crate) const CORNER_TURN: f64 = 30.0;
/// Triangles per tree leaf.
const LEAF: usize = 4;

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Aabb {
    pub lo: P,
    pub hi: P,
}

impl Aabb {
    pub const EMPTY: Aabb = Aabb { lo: [f64::INFINITY; 3], hi: [f64::NEG_INFINITY; 3] };
    pub fn add(&mut self, p: P) {
        for ((lo, hi), x) in self.lo.iter_mut().zip(&mut self.hi).zip(p) {
            *lo = lo.min(x);
            *hi = hi.max(x);
        }
    }
    fn join(&mut self, o: &Aabb) {
        self.add(o.lo);
        self.add(o.hi);
    }
    pub fn overlaps(&self, o: &Aabb) -> bool {
        (0..3).all(|k| self.lo[k] <= o.hi[k] && o.lo[k] <= self.hi[k])
    }
    /// Squared distance from a point to the box (0 inside).
    fn dist2(&self, p: P) -> f64 {
        (0..3).map(|k| (self.lo[k] - p[k]).max(0.0).max(p[k] - self.hi[k])).map(|d| d * d).sum()
    }
    /// Whether the segment `a + t (b − a)`, `t ∈ [0, 1]`, may meet the box (slab test, widened by
    /// a relative margin so rounding never drops a hit).
    fn meets_segment(&self, a: P, b: P) -> bool {
        let (mut t0, mut t1) = (0.0f64, 1.0f64);
        for k in 0..3 {
            let pad = PAD * (self.hi[k] - self.lo[k]).abs().max(1e-300) + ROUND_OFF * (a[k].abs() + b[k].abs());
            let (lo, hi) = (self.lo[k] - pad, self.hi[k] + pad);
            let d = b[k] - a[k];
            if d == 0.0 {
                if a[k] < lo || a[k] > hi {
                    return false;
                }
                continue;
            }
            let (mut u, mut v) = ((lo - a[k]) / d, (hi - a[k]) / d);
            if u > v {
                std::mem::swap(&mut u, &mut v);
            }
            t0 = t0.max(u);
            t1 = t1.min(v);
            if t0 > t1 {
                return false;
            }
        }
        true
    }
}

#[derive(Debug, Clone, Copy)]
struct Node {
    bounds: Aabb,
    children: Option<usize>,
    start: usize,
    end: usize,
}

/// A median-split tree over triangle boxes.
#[derive(Debug, Clone)]
struct Tree {
    nodes: Vec<Node>,
    order: Vec<u32>,
}

impl Tree {
    fn build(boxes: &[Aabb]) -> Tree {
        let mut order: Vec<u32> = (0..boxes.len() as u32).collect();
        let empty = Node { bounds: Aabb::EMPTY, children: None, start: 0, end: 0 };
        let mut nodes = vec![empty];
        let mut work = vec![(0usize, 0usize, boxes.len())];
        while let Some((ni, start, end)) = work.pop() {
            let Some(slice) = order.get_mut(start..end) else { continue };
            let (mut bounds, mut centres) = (Aabb::EMPTY, Aabb::EMPTY);
            for &t in slice.iter() {
                if let Some(b) = boxes.get(t as usize) {
                    bounds.join(b);
                    centres.add([0, 1, 2].map(|k| b.lo[k] + b.hi[k]));
                }
            }
            let mut node = Node { bounds, children: None, start, end };
            if end - start > LEAF {
                let ext = [0, 1, 2].map(|k| centres.hi[k] - centres.lo[k]);
                let axis = if ext[0] >= ext[1] && ext[0] >= ext[2] {
                    0
                } else if ext[1] >= ext[2] {
                    1
                } else {
                    2
                };
                let key = |t: &u32| boxes.get(*t as usize).map_or(0.0, |b| b.lo[axis] + b.hi[axis]);
                let mid = (end - start) / 2;
                slice.select_nth_unstable_by(mid, |a, b| key(a).total_cmp(&key(b)).then(a.cmp(b)));
                let l = nodes.len();
                nodes.push(empty);
                nodes.push(empty);
                node.children = Some(l);
                work.push((l, start, start + mid));
                work.push((l + 1, start + mid, end));
            }
            if let Some(n) = nodes.get_mut(ni) {
                *n = node;
            }
        }
        Tree { nodes, order }
    }

    fn leaf(&self, n: &Node) -> &[u32] {
        self.order.get(n.start..n.end).unwrap_or(&[])
    }
}

pub(crate) fn lerp(a: P, b: P, t: f64) -> P {
    [0, 1, 2].map(|k| a[k] + (b[k] - a[k]) * t)
}

/// The closest point of segment `ab` to `p` and its parameter.
pub(crate) fn closest_on_segment(p: P, a: P, b: P) -> (P, f64) {
    let ab = sub(b, a);
    let l2 = dot(ab, ab);
    let t = if l2 > 0.0 { (dot(sub(p, a), ab) / l2).clamp(0.0, 1.0) } else { 0.0 };
    (lerp(a, b, t), t)
}

pub(crate) fn dist2(a: P, b: P) -> f64 {
    let d = sub(a, b);
    dot(d, d)
}

/// A chain of feature edges.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Curve {
    /// Surface vertices in order; a closed loop repeats its first vertex at the end.
    pub vertices: Vec<u32>,
    pub closed: bool,
    /// The patches on either side (sorted, deduplicated: one entry for a crease seam).
    pub patches: Vec<u32>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Patch {
    pub triangles: Vec<u32>,
    /// The B-rep faces of its triangles (sorted, deduplicated; one for a well-formed boundary).
    pub faces: Vec<u32>,
    pub area: f64,
    /// Twice the area over the length of the bounding curves: the width of a strip, the
    /// diameter of a disc, `∞` for a closed surface without curves.
    pub width: f64,
}

/// The boundary with its queries and feature graph.
pub(crate) struct Geometry<'a> {
    pub surface: &'a BoundarySurface,
    tree: Tree,
    pub boxes: Vec<Aabb>,
    pub normals: Vec<P>,
    pub tri_patch: Vec<u32>,
    pub patches: Vec<Patch>,
    pub curves: Vec<Curve>,
    /// Corner vertices, ascending.
    pub corners: Vec<u32>,
    pub bounds: Aabb,
    /// Far enough to leave the boundary from any point of its box.
    far: f64,
}

impl<'a> Geometry<'a> {
    pub fn new(surface: &'a BoundarySurface, sharp_angle: f64, blade_angle: f64) -> Result<Geometry<'a>> {
        let pts = surface.points();
        let tris = surface.triangles();
        let xyz = |v: u32| pts.get(v as usize).copied().unwrap_or([0.0; 3]);
        let mut bounds = Aabb::EMPTY;
        let mut boxes = Vec::with_capacity(tris.len());
        let mut normals = Vec::with_capacity(tris.len());
        for t in tris {
            let mut b = Aabb::EMPTY;
            for &v in t {
                b.add(xyz(v));
            }
            bounds.join(&b);
            boxes.push(b);
            let n = cross(sub(xyz(t[1]), xyz(t[0])), sub(xyz(t[2]), xyz(t[0])));
            let l = norm(n);
            normals.push(if l > 0.0 { n.map(|c| c / l) } else { [0.0; 3] });
        }
        if tris.is_empty() || !(0..3).all(|k| bounds.hi[k] >= bounds.lo[k]) {
            return Err(Error::Invalid("the boundary has no triangles".into()));
        }
        let tree = Tree::build(&boxes);
        let diag = norm(sub(bounds.hi, bounds.lo));
        let far = 4.0 * diag.max(fieldcraft_model::units::min_geometry_length())
            + (0..3).map(|k| bounds.lo[k].abs().max(bounds.hi[k].abs())).fold(0.0, f64::max);

        // Feature edges and corners.
        // A boundary between two B-rep faces is kept or dropped as a whole, by the median angle
        // along it, so tessellation noise never breaks a curve into pieces; creases inside a
        // face go edge by edge.
        let face_of = |t: u32| surface.tri_face().get(t as usize).copied().unwrap_or(0);
        let mut pair_angles: BTreeMap<(u32, u32), Vec<f64>> = BTreeMap::new();
        for e in surface.feature_edges() {
            let (a, b) = (face_of(e.triangles[0]), face_of(e.triangles[1]));
            if a != b {
                pair_angles.entry((a.min(b), a.max(b))).or_default().push(e.angle);
            }
        }
        let pair_kept: BTreeMap<(u32, u32), bool> = pair_angles
            .into_iter()
            .map(|(k, mut v)| {
                v.sort_by(f64::total_cmp);
                let m = v.get(v.len() / 2).copied().unwrap_or(0.0);
                (k, m >= sharp_angle && m <= blade_angle)
            })
            .collect();
        let kept = |e: &&crate::surface::FeatureEdge| {
            let (a, b) = (face_of(e.triangles[0]), face_of(e.triangles[1]));
            if a != b { pair_kept.get(&(a.min(b), a.max(b))).copied().unwrap_or(false) } else { e.angle >= sharp_angle && e.angle <= blade_angle }
        };
        let mut edge_set: BTreeMap<(u32, u32), [u32; 2]> = BTreeMap::new();
        for e in surface.feature_edges().iter().filter(kept) {
            edge_set.insert((e.vertices[0], e.vertices[1]), e.triangles);
        }
        // Patches: flood fill across non-feature edges.
        let neighbors = surface.neighbors();
        let mut tri_patch = vec![u32::MAX; tris.len()];
        let mut patches: Vec<Patch> = Vec::new();
        for seed in 0..tris.len() {
            if tri_patch[seed] != u32::MAX {
                continue;
            }
            let id = patches.len() as u32;
            let mut patch = Patch { triangles: Vec::new(), faces: Vec::new(), area: 0.0, width: f64::INFINITY };
            tri_patch[seed] = id;
            let mut stack = vec![seed as u32];
            while let Some(t) = stack.pop() {
                patch.triangles.push(t);
                let tri = tris.get(t as usize).copied().unwrap_or([0; 3]);
                let n = neighbors.get(t as usize).copied().unwrap_or([u32::MAX; 3]);
                for k in 0..3 {
                    let (a, b) = (tri[k], tri[(k + 1) % 3]);
                    if edge_set.contains_key(&(a.min(b), a.max(b))) {
                        continue;
                    }
                    if let Some(slot) = tri_patch.get_mut(n[k] as usize)
                        && *slot == u32::MAX
                    {
                        *slot = id;
                        stack.push(n[k]);
                    }
                }
            }
            patch.triangles.sort_unstable();
            for &t in &patch.triangles {
                let tri = tris.get(t as usize).copied().unwrap_or([0; 3]);
                patch.area += 0.5 * norm(cross(sub(xyz(tri[1]), xyz(tri[0])), sub(xyz(tri[2]), xyz(tri[0]))));
                patch.faces.push(surface.tri_face().get(t as usize).copied().unwrap_or(0));
            }
            patch.faces.sort_unstable();
            patch.faces.dedup();
            patches.push(patch);
        }

        // A crease with the same patch on both sides (one inside a face that does not run from
        // feature to feature, such as a tessellation crease across a fillet) separates nothing:
        // no region border can follow it and its open ends could never be mapped, so it is not a
        // feature.
        edge_set.retain(|_, t| t.map(|t| tri_patch.get(t as usize).copied()).iter().collect::<std::collections::BTreeSet<_>>().len() == 2);
        // The widths count only the curves that bound the patches.
        for patch in &mut patches {
            let mut perimeter = 0.0;
            for &t in &patch.triangles {
                let tri = tris.get(t as usize).copied().unwrap_or([0; 3]);
                for k in 0..3 {
                    let (a, b) = (tri[k], tri[(k + 1) % 3]);
                    if edge_set.contains_key(&(a.min(b), a.max(b))) {
                        perimeter += norm(sub(xyz(b), xyz(a)));
                    }
                }
            }
            patch.width = if perimeter > 0.0 { 2.0 * patch.area / perimeter } else { f64::INFINITY };
        }
        let mut adj: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
        for &(a, b) in edge_set.keys() {
            adj.entry(a).or_default().push(b);
            adj.entry(b).or_default().push(a);
        }
        let mut corners: Vec<u32> = adj
            .iter()
            .filter(|(v, n)| {
                if n.len() != 2 {
                    return true;
                }
                let (p, a, b) = (xyz(**v), xyz(n[0]), xyz(n[1]));
                let (u, w) = (sub(p, a), sub(b, p));
                crate::vec3::angle_deg(u, w) > CORNER_TURN
            })
            .map(|(v, _)| *v)
            .collect();
        corners.sort_unstable();

        // Curves: chains from corners, then closed loops.
        let mut used: BTreeMap<(u32, u32), bool> = edge_set.keys().map(|k| (*k, false)).collect();
        let is_corner = |v: u32| corners.binary_search(&v).is_ok();
        let mut curves: Vec<Curve> = Vec::new();
        let trace = |start: u32, next: u32, used: &mut BTreeMap<(u32, u32), bool>, curves: &mut Vec<Curve>| {
            let mut verts = vec![start];
            let (mut prev, mut cur) = (start, next);
            loop {
                let key = (prev.min(cur), prev.max(cur));
                if let Some(u) = used.get_mut(&key) {
                    *u = true;
                }
                verts.push(cur);
                if cur == start || is_corner(cur) || verts.len() > edge_set.len() + 2 {
                    break;
                }
                let nxt = adj.get(&cur).and_then(|n| n.iter().copied().find(|&w| w != prev && used.get(&(cur.min(w), cur.max(w))) == Some(&false)));
                match nxt {
                    Some(w) => {
                        prev = cur;
                        cur = w;
                    }
                    None => break,
                }
            }
            let closed = verts.first() == verts.last() && verts.len() > 2;
            let mut ps: Vec<u32> = Vec::new();
            for w in verts.windows(2) {
                if let Some(t) = edge_set.get(&(w[0].min(w[1]), w[0].max(w[1]))) {
                    ps.extend(t.iter().filter_map(|&t| tri_patch.get(t as usize).copied()));
                }
            }
            ps.sort_unstable();
            ps.dedup();
            curves.push(Curve { vertices: verts, closed, patches: ps });
        };
        for &c in &corners {
            let ns = adj.get(&c).cloned().unwrap_or_default();
            for n in ns {
                if used.get(&(c.min(n), c.max(n))) == Some(&false) {
                    trace(c, n, &mut used, &mut curves);
                }
            }
        }
        let keys: Vec<(u32, u32)> = used.keys().copied().collect();
        for k in keys {
            if used.get(&k) == Some(&false) {
                trace(k.0, k.1, &mut used, &mut curves);
            }
        }
        Ok(Geometry { surface, tree, boxes, normals, tri_patch, patches, curves, corners, bounds, far })
    }

    pub fn xyz(&self, v: u32) -> P {
        self.surface.points().get(v as usize).copied().unwrap_or([0.0; 3])
    }

    /// The centroid of triangle `t`.
    pub fn triangle_centre(&self, t: u32) -> P {
        let [a, b, c] = self.triangle(t);
        [0, 1, 2].map(|k| (a[k] + b[k] + c[k]) / 3.0)
    }

    pub fn triangle(&self, t: u32) -> [P; 3] {
        let tri = self.surface.triangles().get(t as usize).copied().unwrap_or([0; 3]);
        tri.map(|v| self.xyz(v))
    }

    /// The closest boundary point to `p` among triangles that pass `keep`: (triangle, point,
    /// squared distance).
    pub fn nearest(&self, p: P, keep: impl Fn(u32) -> bool) -> Option<(u32, P, f64)> {
        let mut best: Option<(u32, P, f64)> = None;
        let mut stack = vec![0usize];
        while let Some(n) = stack.pop() {
            let Some(node) = self.tree.nodes.get(n) else { continue };
            let bd = node.bounds.dist2(p);
            if best.is_some_and(|b| bd >= b.2) {
                continue;
            }
            match node.children {
                Some(l) => {
                    let (dl, dr) = (
                        self.tree.nodes.get(l).map_or(f64::INFINITY, |x| x.bounds.dist2(p)),
                        self.tree.nodes.get(l + 1).map_or(f64::INFINITY, |x| x.bounds.dist2(p)),
                    );
                    if dl <= dr {
                        stack.extend([l + 1, l]);
                    } else {
                        stack.extend([l, l + 1]);
                    }
                }
                None => {
                    for &t in self.tree.leaf(node) {
                        if !keep(t) {
                            continue;
                        }
                        let [a, b, c] = self.triangle(t);
                        let q = closest_on_triangle(p, a, b, c);
                        let d = dist2(p, q);
                        if best.is_none_or(|x| d < x.2 || (d == x.2 && t < x.0)) {
                            best = Some((t, q, d));
                        }
                    }
                }
            }
        }
        best
    }

    /// The triangles whose boxes overlap `b`, ascending.
    pub fn overlapping(&self, b: &Aabb) -> Vec<u32> {
        let mut out = Vec::new();
        let mut stack = vec![0usize];
        while let Some(n) = stack.pop() {
            let Some(node) = self.tree.nodes.get(n) else { continue };
            if !node.bounds.overlaps(b) {
                continue;
            }
            match node.children {
                Some(l) => stack.extend([l, l + 1]),
                None => out.extend(self.tree.leaf(node).iter().copied().filter(|&t| self.boxes.get(t as usize).is_some_and(|x| x.overlaps(b)))),
            }
        }
        out.sort_unstable();
        out
    }

    /// The first triangle hit by the ray from `p` along `d` (unit) within `max`, with its
    /// distance (floating point; used to measure, never to decide).
    pub fn ray(&self, p: P, d: P, max: f64, skip: u32) -> Option<(u32, f64)> {
        let q = [0, 1, 2].map(|k| p[k] + d[k] * max);
        let mut best: Option<(u32, f64)> = None;
        let mut stack = vec![0usize];
        while let Some(n) = stack.pop() {
            let Some(node) = self.tree.nodes.get(n) else { continue };
            if !node.bounds.meets_segment(p, q) {
                continue;
            }
            match node.children {
                Some(l) => stack.extend([l, l + 1]),
                None => {
                    for &t in self.tree.leaf(node) {
                        if t == skip {
                            continue;
                        }
                        let [a, b, c] = self.triangle(t);
                        if let Some(s) = ray_triangle(p, d, a, b, c)
                            && s > 0.0
                            && s <= max
                            && best.is_none_or(|x| s < x.1)
                        {
                            best = Some((t, s));
                        }
                    }
                }
            }
        }
        best
    }

    /// Whether `p` is strictly inside the solid, decided exactly: the parity of the boundary
    /// crossings of a segment from `p` to a far point, tried along fixed directions until one
    /// meets no edge or vertex and does not start on the boundary. `None` when every direction is
    /// degenerate (p lies on the boundary) or `p` is not admissible.
    pub fn inside(&self, p: P) -> Option<bool> {
        const DIRS: [P; 7] = [
            [0.5773, 0.6123, 0.5401],
            [-0.6301, 0.4372, 0.6417],
            [0.3989, -0.7068, 0.5839],
            [0.7211, 0.4109, -0.5577],
            [-0.4402, -0.6005, -0.6675],
            [0.1351, 0.9013, -0.4116],
            [-0.8803, 0.1237, 0.4583],
        ];
        let pp = Point3::new(snap_to_domain(p[0]), snap_to_domain(p[1]), snap_to_domain(p[2])).ok()?;
        'dirs: for d in DIRS {
            let q = [0, 1, 2].map(|k| snap_to_domain(p[k] + d[k] * self.far));
            let Ok(qq) = Point3::new(q[0], q[1], q[2]) else { continue };
            let mut parity = false;
            let mut stack = vec![0usize];
            while let Some(n) = stack.pop() {
                let Some(node) = self.tree.nodes.get(n) else { continue };
                if !node.bounds.meets_segment(p, q) {
                    continue;
                }
                match node.children {
                    Some(l) => stack.extend([l, l + 1]),
                    None => {
                        for &t in self.tree.leaf(node) {
                            match self.crossing(&pp, &qq, t) {
                                Some(true) => parity = !parity,
                                Some(false) => {}
                                None => continue 'dirs,
                            }
                        }
                    }
                }
            }
            return Some(parity);
        }
        None
    }

    /// Whether segment `pq` crosses triangle `t` properly (`Some(true)`), misses it
    /// (`Some(false)`), or touches it degenerately (`None`): exact.
    fn crossing(&self, p: &Point3, q: &Point3, t: u32) -> Option<bool> {
        let tri = self.surface.triangles().get(t as usize)?;
        let a = self.surface.point3(tri[0])?;
        let b = self.surface.point3(tri[1])?;
        let c = self.surface.point3(tri[2])?;
        let (sp, sq) = (orient3d(&a, &b, &c, p), orient3d(&a, &b, &c, q));
        if sq.is_zero() {
            // The far end is outside the boundary's box, so never on a triangle; a segment in
            // the plane is ambiguous.
            return None;
        }
        if sp.is_zero() {
            // p in the triangle's plane: the segment leaves the plane at p, so it touches the
            // triangle only if p lies in it.
            let axis = [Axis::Z, Axis::X, Axis::Y].into_iter().find(|&ax| !orient2d_projected(&a, &b, &c, ax).is_zero())?;
            let o = orient2d_projected(&a, &b, &c, axis);
            let inside = [(&a, &b), (&b, &c), (&c, &a)].into_iter().all(|(u, v)| orient2d_projected(u, v, p, axis) != -o);
            return if inside { None } else { Some(false) };
        }
        if sp == sq {
            return Some(false);
        }
        let s = [orient3d(p, q, &a, &b), orient3d(p, q, &b, &c), orient3d(p, q, &c, &a)];
        if s.iter().any(|x| x.is_zero()) {
            let nz: Vec<Sign> = s.iter().copied().filter(|x| !x.is_zero()).collect();
            return if nz.windows(2).all(|w| w[0] == w[1]) { None } else { Some(false) };
        }
        Some(s[0] == s[1] && s[1] == s[2])
    }

    /// The closest point of curve `c` to `p`: (point, segment index, parameter on it).
    pub fn nearest_on_curve(&self, c: u32, p: P) -> Option<(P, usize, f64)> {
        let curve = self.curves.get(c as usize)?;
        let mut best: Option<(P, usize, f64, f64)> = None;
        for (i, w) in curve.vertices.windows(2).enumerate() {
            let (q, t) = closest_on_segment(p, self.xyz(w[0]), self.xyz(w[1]));
            let d = dist2(p, q);
            if best.is_none_or(|b| d < b.3) {
                best = Some((q, i, t, d));
            }
        }
        best.map(|(q, i, t, _)| (q, i, t))
    }
}

/// Möller–Trumbore ray–triangle distance (floating point).
fn ray_triangle(o: P, d: P, a: P, b: P, c: P) -> Option<f64> {
    let (e1, e2) = (sub(b, a), sub(c, a));
    let h = cross(d, e2);
    let det = dot(e1, h);
    if det.abs() < 1e-300 {
        return None;
    }
    let inv = 1.0 / det;
    let s = sub(o, a);
    let u = inv * dot(s, h);
    if !(-ROUND_OFF..=1.0 + ROUND_OFF).contains(&u) {
        return None;
    }
    let q = cross(s, e1);
    let v = inv * dot(d, q);
    if v < -ROUND_OFF || u + v > 1.0 + ROUND_OFF {
        return None;
    }
    Some(inv * dot(e2, q))
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::surface::SurfaceOptions;

    /// A closed axis-aligned box surface, one B-rep face per side.
    pub(crate) fn box_surface(lo: P, hi: P) -> BoundarySurface {
        let c = |i: usize| [if i & 1 == 0 { lo[0] } else { hi[0] }, if i & 2 == 0 { lo[1] } else { hi[1] }, if i & 4 == 0 { lo[2] } else { hi[2] }];
        let pts: Vec<P> = (0..8).map(c).collect();
        let quads = [[0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [2, 6, 7, 3], [0, 4, 6, 2], [1, 3, 7, 5]];
        let mut tris = Vec::new();
        let mut face = Vec::new();
        for (f, q) in quads.iter().enumerate() {
            tris.push([q[0], q[1], q[2]]);
            tris.push([q[0], q[2], q[3]]);
            face.extend([f as u32; 2]);
        }
        BoundarySurface::new(&pts, &tris, &face, &SurfaceOptions::default()).unwrap()
    }

    #[test]
    fn box_feature_graph_inside_and_nearest() {
        let s = box_surface([0.0; 3], [2.0, 3.0, 4.0]);
        let g = Geometry::new(&s, 30.0, 160.0).unwrap();
        assert_eq!(g.corners.len(), 8);
        assert_eq!(g.curves.len(), 12);
        assert_eq!(g.patches.len(), 6);
        assert!(g.curves.iter().all(|c| c.vertices.len() == 2 && c.patches.len() == 2 && !c.closed));
        assert_eq!(g.inside([1.0, 1.5, 2.0]), Some(true));
        assert_eq!(g.inside([3.0, 1.5, 2.0]), Some(false));
        assert_eq!(g.inside([1.0, 1.5, 0.0]), None, "on the boundary");
        assert_eq!(g.inside([0.0, 0.0, 0.0]), None, "a corner");
        let (_, q, d) = g.nearest([1.0, 1.0, -1.0], |_| true).unwrap();
        assert!((d - 1.0).abs() < 1e-12 && q[2] == 0.0);
        let top = g.tri_patch[2];
        let (_, q, _) = g.nearest([1.0, 1.0, -1.0], |t| g.tri_patch[t as usize] == top).unwrap();
        assert_eq!(q[2], 4.0, "restricted to the top patch");
        let (t, s) = g.ray([1.0, 1.5, 2.0], [0.0, 0.0, 1.0], 10.0, u32::MAX).unwrap();
        assert!((s - 2.0).abs() < 1e-12 && g.tri_patch[t as usize] == top);
    }
}
