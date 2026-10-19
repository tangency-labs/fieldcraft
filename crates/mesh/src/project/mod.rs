//! Projection blocking: a block topology laid over the part and projected onto its B-rep.
//!
//! Instead of cutting the geometry into sweepable cells, the analyst lays down a clean block
//! topology — a grid with blocks carved away, butterfly (O-grid) templates for round sections,
//! the seven-block sphere — whose outer faces stand for the part's surface. Every block face on
//! the outside of the topology is projected onto a target: one B-rep face, or a union of faces
//! (by default the smooth group of faces nearest it). Where two outer block faces with different
//! targets meet, their edge goes onto the curve where the targets meet, and a vertex on three
//! targets onto their common point ([`target`]). Nodes are then placed outside in: block
//! vertices, edge nodes (evenly by arc length along their curve), face nodes by transfinite
//! interpolation of their edges projected onto the face's target, and interior nodes by the
//! transfinite interpolation of the block's faces ([`crate::quad::transfinite_volume`], Gordon
//! and Hall 1973, Gordon and Thiel 1982). Elliptic smoothing ([`smooth`]) finishes the grid with
//! the boundary held on the geometry.
//!
//! The method follows the public descriptions of block-structured projection meshing (blocks in
//! index space, carved, projected face → surface, edge → curve, vertex → point, then smoothed);
//! the algorithms are the papers cited in each module.

pub mod device;
pub mod fixed;
mod lattice;
pub mod smooth;
pub mod target;
pub mod topology;

#[cfg(test)]
mod bench;
#[cfg(test)]
mod tests;

use std::collections::{BTreeMap, BTreeSet, HashMap};

use fieldcraft_model::ElementKind;

pub use fixed::FixedPatch;
pub use smooth::{QUALITY_FLOOR, SmoothReport, Smoothing};
pub use target::Geometry;
pub use topology::Topology;

use crate::hex::BoundaryFace;
use crate::message::Message;
use crate::quality::{Metric, element_quality};
use crate::surface::{BoundarySurface, FeatureKind};
use crate::vec3::{P, add, dist, lerp, scale};
use crate::{Error, Result};
use lattice::{Entity, build, edge_classes, face_key, face_uses};
use target::{Constraint, on_target, place};
use topology::{CORNERS, FACE_SIDES, corner_at, trilinear};

/// Faces meeting at less than this angle (degrees) along their shared edge belong to one smooth
/// group, the default target of a block face.
pub const SMOOTH_ANGLE: f64 = 20.0;

/// The default [`ProjectOptions::min_quality`]: a corner flattened toward 180° (a block corner
/// projected onto a smooth surface) scores far below it.
pub const DEFAULT_MIN_QUALITY: f64 = 0.2;

/// Sweeps of the inner vertices' displacement blend.
const VERTEX_SWEEPS: usize = 500;

/// Rounds of arc-length redistribution of an edge's nodes.
const EDGE_ROUNDS: usize = 4;

/// The default largest turn of the surface across one element: 30°.
pub const DEFAULT_MAX_ARC_ANGLE: f64 = std::f64::consts::PI / 6.0;

/// Samples per side of an outer block face when measuring how much its target turns.
const CURVATURE_SAMPLES: usize = 32;

/// The longest an element edge may be, in element sizes, when a class of parallel block edges
/// takes its divisions from the mean of its lengths.
const MAX_STRETCH: f64 = 2.0;

/// How to project.
#[derive(Debug, Clone, PartialEq)]
pub struct ProjectOptions {
    /// The element size: each class of parallel block edges gets `round(length / size)`
    /// divisions (at least one), `length / size` the geometric mean of its projected edges'
    /// lengths in sizes (measured by `field` when there is one), and enough that no element is
    /// longer than twice the size.
    pub size: f64,
    /// The largest angle (radians) the surface may turn through across one element on a curved
    /// target face: a class takes enough divisions that its elements are no longer than the
    /// face's radius of curvature times this (0 for no limit).
    pub max_arc_angle: f64,
    /// Division counts set by hand: `(block, direction 0..3, count)`; the count applies to the
    /// whole class of parallel edges.
    pub divisions: Vec<(usize, usize, usize)>,
    /// Targets set by hand: `(block, hex8 face, B-rep faces)`. Other outer faces take the smooth
    /// group of B-rep faces nearest them.
    pub faces: Vec<(usize, u8, Vec<u32>)>,
    pub smoothing: Smoothing,
    pub iterations: usize,
    pub max_elements: usize,
    /// The smallest scaled Jacobian accepted: a mesh with a worse element is an error that
    /// says where a template is missing (0 accepts any valid mesh).
    pub min_quality: f64,
    /// Already-meshed patches of the boundary taken as given ([`fixed`]): the block faces on
    /// them take their divisions and nodes from the patch, held through smoothing.
    pub fixed: Vec<FixedPatch>,
    /// Planes `(point, unit normal)` inner block faces are held on: a sheet of block faces whose
    /// corners lie on one (within a tenth of the size) keeps its nodes on it through smoothing, so
    /// the mesh can be cut there into two meshes with planar faces.
    pub planes: Vec<(P, P)>,
    /// A graded size field (size controls) in place of the uniform `size` when the divisions are
    /// chosen: an edge's length is then measured in units of the field ([`crate::size::ideal_divisions`]).
    /// Its global size is `size`.
    pub field: Option<crate::size::RegionSize>,
    /// Where the smoothing sweeps run.
    pub backend: Backend,
}

/// Where a projection's smoothing sweeps run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Backend {
    /// On the CPU, node by node: the reference, the same mesh on every machine.
    #[default]
    Cpu,
    /// On the CPU, colour by colour (no two nodes of a colour share an element), each colour in
    /// parallel: the same mesh on every machine and thread count, a slightly different one from
    /// [`Backend::Cpu`]'s; the reference of [`Backend::Device`].
    Coloured,
    /// Colour by colour with the free nodes on the registered sweep device ([`device`]): single
    /// precision there, deterministic per device, the same verdicts as [`Backend::Coloured`] and a
    /// smallest scaled Jacobian within 0.02 of it; an error when none is registered.
    Device,
}

impl Default for ProjectOptions {
    fn default() -> Self {
        ProjectOptions {
            size: 1.0,
            max_arc_angle: DEFAULT_MAX_ARC_ANGLE,
            divisions: Vec::new(),
            faces: Vec::new(),
            smoothing: Smoothing::Elliptic,
            iterations: 100,
            max_elements: 2_000_000,
            min_quality: DEFAULT_MIN_QUALITY,
            fixed: Vec::new(),
            planes: Vec::new(),
            field: None,
            backend: Backend::Cpu,
        }
    }
}

/// How it went.
#[derive(Debug, Clone, PartialEq)]
pub struct ProjectReport {
    pub blocks: usize,
    pub elements: usize,
    /// The smallest scaled Jacobian (corners) of the transfinite grid, before smoothing.
    pub initial_min_scaled_jacobian: f64,
    /// Certified smallest and mean scaled Jacobian of the result.
    pub min_scaled_jacobian: f64,
    pub mean_scaled_jacobian: f64,
    pub smoothing: SmoothReport,
    /// Farthest a boundary node lies from its target (surface, curve or point).
    pub max_gap: f64,
    pub volume: f64,
    pub surface_volume: f64,
    /// The B-rep faces each target stands for.
    pub targets: Vec<Vec<u32>>,
    /// Divisions per block.
    pub divisions: Vec<[usize; 3]>,
    /// Nodes taken from the fixed patches, held where they are.
    pub fixed_nodes: usize,
}

/// An all-hex mesh of a projected topology.
#[derive(Debug, Clone, PartialEq)]
pub struct ProjectedMesh {
    pub nodes: Vec<P>,
    pub hexes: Vec<[u32; 8]>,
    pub boundary: Vec<BoundaryFace>,
    pub report: ProjectReport,
}

/// Give each target face no boundary quad's centre landed on (a face narrower than the elements:
/// a fillet's sliver between a cut and the walls) one row of the quads it runs through (of a
/// target holding the face, centre within half the quad's shortest side of it): one quad of
/// those beside each other across the face (projecting onto it at nearly the same point), all
/// along it. A quad keeps its face when it is that face's
/// last, so no face loses its only quad to a neighbour.
fn claim_narrow_faces(boundary: &mut [BoundaryFace], centres: &[(P, f64, usize)], targets: &[Vec<u32>], g: &dyn Geometry) {
    let mut count: BTreeMap<u32, usize> = BTreeMap::new();
    for b in boundary.iter() {
        *count.entry(b.brep_face).or_insert(0) += 1;
    }
    let targeted: BTreeSet<u32> = targets.iter().flatten().copied().collect();
    for f in targeted {
        if count.contains_key(&f) {
            continue;
        }
        // (distance to the face, quad, its centre, the centre's foot on the face, shortest side)
        let mut near: Vec<(f64, usize, P, P, f64)> = centres
            .iter()
            .enumerate()
            .filter(|(_, (_, _, tg))| targets.get(*tg).is_some_and(|set| set.contains(&f)))
            .filter_map(|(i, (centre, shortest, _))| {
                let foot = g.project(f, *centre)?;
                let d = dist(foot, *centre);
                (d <= 0.5 * shortest).then_some((d, i, *centre, foot, *shortest))
            })
            .collect();
        near.sort_by(|a, b| a.0.total_cmp(&b.0));
        // The rows across the face: quads beside each other that project onto it at nearly the
        // same point.
        let mut rows: Vec<Vec<(usize, P, P, f64)>> = Vec::new();
        for (_, i, centre, foot, shortest) in near {
            let across = |(_, c, q, s): &(usize, P, P, f64)| {
                let apart = dist(*c, centre);
                apart < 1.5 * s.max(shortest) && dist(*q, foot) < 0.5 * apart
            };
            match rows.iter_mut().find(|row| row.iter().any(across)) {
                Some(row) => row.push((i, centre, foot, shortest)),
                None => rows.push(vec![(i, centre, foot, shortest)]),
            }
        }
        // One quad per row: from the neighbour with the most quads left (the nearest of equals),
        // so neighbours on both sides give alike.
        for row in rows {
            let left = |i: usize| boundary.get(i).and_then(|b| count.get(&b.brep_face)).copied().unwrap_or(0);
            let Some(&(i, ..)) = row.iter().max_by_key(|(i, ..)| left(*i)) else { continue };
            let Some(b) = boundary.get_mut(i) else { continue };
            let Some(left) = count.get_mut(&b.brep_face) else { continue };
            if *left > 1 {
                *left -= 1;
                b.brep_face = f;
                *count.entry(f).or_insert(0) += 1;
            }
        }
    }
}

/// Most planes inner faces are held on.
pub const MAX_PLANES: usize = 64;

/// The geometry with planes added as faces `first`, `first + 1`, ….
struct Planed<'a> {
    g: &'a dyn Geometry,
    first: u32,
    planes: &'a [(P, P)],
}

impl Geometry for Planed<'_> {
    fn project(&self, face: u32, p: P) -> Option<P> {
        match face.checked_sub(self.first) {
            Some(k) => {
                let (o, n) = self.planes.get(k as usize)?;
                Some(crate::vec3::sub(p, scale(*n, crate::vec3::dot(crate::vec3::sub(p, *o), *n))))
            }
            None => self.g.project(face, p),
        }
    }
    fn normal(&self, face: u32, p: P) -> Option<P> {
        match face.checked_sub(self.first) {
            Some(k) => self.planes.get(k as usize).map(|x| x.1),
            None => self.g.normal(face, p),
        }
    }
}

/// The face groups joined across smooth B-rep edges.
fn smooth_groups(surface: &BoundarySurface, faces: usize) -> Vec<usize> {
    let mut parent: Vec<usize> = (0..faces).collect();
    fn root(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    let tf = surface.tri_face();
    for e in surface.feature_edges() {
        if e.kind != FeatureKind::FaceBoundary || e.angle >= SMOOTH_ANGLE {
            continue;
        }
        let (Some(&a), Some(&b)) = (tf.get(e.triangles[0] as usize), tf.get(e.triangles[1] as usize)) else { continue };
        let (a, b) = (a as usize, b as usize);
        if a < faces && b < faces {
            let (ra, rb) = (root(&mut parent, a), root(&mut parent, b));
            if ra != rb {
                parent[ra.max(rb)] = ra.min(rb);
            }
        }
    }
    (0..faces).map(|f| root(&mut parent, f)).collect()
}

/// The B-rep face a sample `u` of block face `(d, s)` stands for: the nearest face that faces the
/// same way (normals within 60°), else the nearest face.
fn facing(g: &dyn Geometry, faces: &[u32], c: &[P; 8], d: usize, s: usize, u: [f64; 3]) -> Option<u32> {
    let p = trilinear(c, u);
    let (e1, e2) = other(d);
    let h = 1e-3;
    let diff = |k: usize| {
        let (mut a, mut b) = (u, u);
        a[k] -= h;
        b[k] += h;
        crate::vec3::sub(trilinear(c, b), trilinear(c, a))
    };
    let mut n = crate::vec3::cross(diff(e1), diff(e2));
    // Outward: along +d on side 1, against it on side 0.
    let inward = {
        let mut q = u;
        q[d] = if s == 1 { 1.0 - 0.1 } else { 0.1 };
        crate::vec3::sub(trilinear(c, q), p)
    };
    if crate::vec3::dot(n, inward) > 0.0 {
        n = scale(n, -1.0);
    }
    let ln = crate::vec3::norm(n);
    let mut best: Option<(f64, u32)> = None;
    let mut nearest: Option<(f64, u32)> = None;
    for &f in faces {
        let Some(q) = g.project(f, p) else { continue };
        let dq = dist(p, q);
        if nearest.is_none_or(|b| dq < b.0) {
            nearest = Some((dq, f));
        }
        let agrees = ln > 0.0 && g.normal(f, q).is_some_and(|m| crate::vec3::dot(m, n) > 0.5 * ln);
        if agrees && best.is_none_or(|b| dq < b.0) {
            best = Some((dq, f));
        }
    }
    best.or(nearest).map(|b| b.1)
}

/// The corners of block face `(d, s)` in cyclic order.
fn face_cycle(block: &[u32; 8], d: usize, s: usize) -> [u32; 4] {
    let (e1, e2) = other(d);
    [[0, 0], [1, 0], [1, 1], [0, 1]].map(|[a, b]| {
        let mut bits = [0; 3];
        bits[d] = s;
        bits[e1] = a;
        bits[e2] = b;
        block[corner_at(bits)]
    })
}

pub(crate) fn other(d: usize) -> (usize, usize) {
    match d {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    }
}

fn edge_key(a: u32, b: u32) -> [u32; 2] {
    [a.min(b), a.max(b)]
}

/// Mesh the topology `t` projected onto the body bounded by `surface`, through `g`.
pub fn project_mesh(t: &Topology, surface: &BoundarySurface, g: &dyn Geometry, opts: &ProjectOptions) -> Result<ProjectedMesh> {
    let order = match opts.backend {
        Backend::Cpu => smooth::Order::Sequential,
        Backend::Coloured => smooth::Order::Coloured,
        Backend::Device => smooth::Order::Device,
    };
    project_ordered(t, surface, g, opts, order)
}

/// [`project_mesh`] with the smoothing sweeps in `order`.
pub(crate) fn project_ordered(
    t: &Topology,
    surface: &BoundarySurface,
    g: &dyn Geometry,
    opts: &ProjectOptions,
    order: smooth::Order,
) -> Result<ProjectedMesh> {
    t.check()?;
    if !(opts.size > 0.0 && opts.size.is_finite()) {
        return Err(Error::invalid(crate::message!["the element size must be a positive number, got ", crate::message::length(opts.size)]));
    }
    let uses = face_uses(t)?;
    let face_count = surface.tri_face().iter().max().map_or(1, |m| *m as usize + 1);
    let groups = smooth_groups(surface, face_count);
    let diag = t.diagonal().max(opts.size);
    let tol = 1e-9 * diag;
    if opts.planes.len() > MAX_PLANES
        || opts.planes.iter().any(|(p, n)| !(p.iter().chain(n).all(|x| x.is_finite()) && (crate::vec3::norm(*n) - 1.0).abs() < 1e-6))
    {
        return Err(Error::Invalid(format!("at most {MAX_PLANES} planes, each a finite point and a unit normal")));
    }
    let planed = Planed { g, first: face_count as u32, planes: &opts.planes };
    let g: &dyn Geometry = &planed;

    // Targets of the outer block faces.
    let explicit: BTreeMap<(usize, u8), Vec<u32>> = opts.faces.iter().map(|(b, f, s)| ((*b, *f), s.clone())).collect();
    for ((b, f), set) in &explicit {
        let (d, s) = FACE_SIDES
            .get(usize::from(*f).wrapping_sub(1))
            .copied()
            .ok_or_else(|| Error::Invalid(format!("face {f} is not a hex8 face (1 to 6)")))?;
        let block = t.blocks.get(*b).ok_or_else(|| Error::Invalid(format!("there is no block {b}")))?;
        if uses.get(&face_key(block, d, s)) != Some(&1) {
            return Err(Error::Invalid(format!("face {f} of block {b} is inside the topology; only outer faces are projected")));
        }
        if set.is_empty() || set.iter().any(|&x| x as usize >= face_count) {
            return Err(Error::Invalid(format!("face {f} of block {b}: the target must name faces of the body")));
        }
    }
    let mut targets: Vec<Vec<u32>> = Vec::new();
    let mut target_of = |set: Vec<u32>| -> usize {
        match targets.iter().position(|t| *t == set) {
            Some(i) => i,
            None => {
                targets.push(set);
                targets.len() - 1
            }
        }
    };
    let all_faces: Vec<u32> = (0..face_count as u32).collect();
    let mut face_target: HashMap<[u32; 4], usize> = HashMap::new();
    let mut outer: Vec<(usize, usize, usize)> = Vec::new();
    for (b, block) in t.blocks.iter().enumerate() {
        let c = t.corners(b)?;
        for (fi, (d, s)) in FACE_SIDES.iter().copied().enumerate() {
            let key = face_key(block, d, s);
            if uses.get(&key) != Some(&1) {
                continue;
            }
            outer.push((b, d, s));
            let set = match explicit.get(&(b, fi as u8 + 1)) {
                Some(set) => {
                    let mut set = set.clone();
                    set.sort_unstable();
                    set.dedup();
                    set
                }
                None => {
                    // The smooth group most samples of the face lie nearest.
                    let (e1, e2) = other(d);
                    let mut votes: BTreeMap<usize, usize> = BTreeMap::new();
                    for (a, bb) in [(0.5, 0.5), (0.25, 0.25), (0.75, 0.25), (0.25, 0.75), (0.75, 0.75)] {
                        let mut u = [0.0; 3];
                        u[d] = s as f64;
                        u[e1] = a;
                        u[e2] = bb;
                        if let Some(f) = facing(g, &all_faces, &c, d, s, u) {
                            *votes.entry(groups.get(f as usize).copied().unwrap_or(0)).or_insert(0) += 1;
                        }
                    }
                    let group = votes
                        .iter()
                        .max_by_key(|(gr, n)| (**n, std::cmp::Reverse(**gr)))
                        .map(|(gr, _)| *gr)
                        .ok_or_else(|| Error::Failed(format!("face {} of block {b} has no B-rep face to project onto", fi + 1)))?;
                    (0..face_count as u32).filter(|&f| groups.get(f as usize) == Some(&group)).collect()
                }
            };
            face_target.insert(key, target_of(set));
        }
    }
    // Inner block faces on a plane: a target of their own.
    let mut sheets: Vec<(usize, usize, usize)> = Vec::new();
    let mut sheet_targets: BTreeSet<usize> = BTreeSet::new();
    for (b, block) in t.blocks.iter().enumerate() {
        for (d, s) in FACE_SIDES {
            let key = face_key(block, d, s);
            if uses.get(&key) != Some(&2) || face_target.contains_key(&key) {
                continue;
            }
            let on = opts.planes.iter().position(|(p, n)| {
                key.iter()
                    .all(|&v| t.vertices.get(v as usize).is_some_and(|x| crate::vec3::dot(crate::vec3::sub(*x, *p), *n).abs() <= 0.1 * opts.size))
            });
            if let Some(k) = on {
                let tg = target_of(vec![face_count as u32 + k as u32]);
                face_target.insert(key, tg);
                sheet_targets.insert(tg);
                sheets.push((b, d, s));
            }
        }
    }
    // Targets of the vertices and edges on the outside.
    let mut vertex_on: BTreeMap<u32, BTreeSet<usize>> = BTreeMap::new();
    let mut edge_on: BTreeMap<[u32; 2], BTreeSet<usize>> = BTreeMap::new();
    for &(b, d, s) in outer.iter().chain(&sheets) {
        let block = t.blocks.get(b).copied().unwrap_or([0; 8]);
        let tg = face_target.get(&face_key(&block, d, s)).copied().unwrap_or(0);
        let cyc = face_cycle(&block, d, s);
        for i in 0..4 {
            vertex_on.entry(cyc[i]).or_default().insert(tg);
            edge_on.entry(edge_key(cyc[i], cyc[(i + 1) % 4])).or_default().insert(tg);
        }
    }
    let as_vec = |s: Option<&BTreeSet<usize>>| s.map(|s| s.iter().copied().collect::<Vec<_>>()).unwrap_or_default();
    let vertex_c = |v: u32| Constraint::of(&as_vec(vertex_on.get(&v)));
    let edge_c = |e: [u32; 2]| Constraint::of(&as_vec(edge_on.get(&e)));

    // 1. Vertices.
    fieldcraft_model::progress::stage("Projection: block vertices", None, "");
    fieldcraft_model::progress::check()?;
    let mut vpos: Vec<P> = t.vertices.clone();
    for (v, p) in vpos.iter_mut().enumerate() {
        *p = place(g, &targets, vertex_c(v as u32), *p, tol);
    }
    // Vertices on fixed patches go to their patch nodes.
    let fixed = if opts.fixed.is_empty() {
        fixed::Fixed::default()
    } else {
        let target_of = |k: &[u32; 4]| face_target.get(k).and_then(|&i| targets.get(i)).cloned();
        fixed::match_patches(t, &outer, &target_of, &opts.fixed, &vpos)?
    };
    for (&v, &p) in &fixed.vertex {
        if let Some(slot) = vpos.get_mut(v as usize) {
            *slot = p;
        }
    }
    // Inner vertices follow: the projected vertices' moves blended through the block edges
    // (a bounded Jacobi solve of Laplace's equation for the displacement), so a block the
    // projection stretches does not leave its inner corners behind.
    let mut nbr: Vec<BTreeSet<u32>> = vec![BTreeSet::new(); vpos.len()];
    for block in &t.blocks {
        for (k, bits) in CORNERS.iter().enumerate() {
            for d in 0..3 {
                let mut e = *bits;
                e[d] = 1 - e[d];
                if let Some(set) = nbr.get_mut(block[k] as usize) {
                    set.insert(block[corner_at(e)]);
                }
            }
        }
    }
    let free: Vec<bool> = (0..vpos.len()).map(|v| vertex_c(v as u32) == Constraint::Free).collect();
    if free.iter().any(|f| *f) {
        let mut disp: Vec<P> = vpos.iter().zip(&t.vertices).map(|(a, b)| crate::vec3::sub(*a, *b)).collect();
        for _ in 0..VERTEX_SWEEPS {
            let mut change = 0.0f64;
            for v in 0..disp.len() {
                if !free[v] {
                    continue;
                }
                let Some(ns) = nbr.get(v).filter(|n| !n.is_empty()) else { continue };
                let mean = scale(ns.iter().fold([0.0; 3], |a, &w| add(a, disp[w as usize])), 1.0 / ns.len() as f64);
                change = change.max(dist(mean, disp[v]));
                disp[v] = mean;
            }
            if change <= tol {
                break;
            }
        }
        for v in 0..vpos.len() {
            if free[v] {
                vpos[v] = add(t.vertices[v], disp[v]);
            }
        }
    }
    // 2. Divisions per class of parallel edges: the count whose element length is nearest the
    // size, in ratio, over all the class's edges (the geometric mean of their projected
    // lengths). A class runs through every block its edges are parallel in, so its edges can
    // differ a lot in length (a cup's wall, 4 mm thick, and its floor's radial edges, twice
    // that); the longest edge alone would cut the thin blocks into flat slivers. No element is
    // longer than [`MAX_STRETCH`] sizes, so a coarse size still follows the geometry. An edge
    // that bends with a curved target is the exception: its elements' sides are chords of the
    // surface, and one longer than the size cuts off more of the part than the size allows (a
    // butterfly's rim shares its class with the shorter sides of its middle square), so a class
    // takes at least the geometric mean of its bent edges' lengths in sizes (rounded,
    // as the mean; an annulus's rim and bore balance each other as its straight edges do).
    let (classes, nclass) = edge_classes(t);
    let mut class_log = vec![(0.0f64, 0usize); nclass];
    // Per class, its longest edge in element sizes.
    let mut class_max = vec![0.0f64; nclass];
    let mut class_bent = vec![(0.0f64, 0usize); nclass];
    let mut seen: BTreeSet<[u32; 2]> = BTreeSet::new();
    for (block, cls) in t.blocks.iter().zip(&classes) {
        for d in 0..3 {
            for (k, bits) in CORNERS.iter().enumerate() {
                if bits[d] != 0 {
                    continue;
                }
                let mut e = *bits;
                e[d] = 1;
                let (a, b) = (block[k], block[corner_at(e)]);
                if !seen.insert(edge_key(a, b)) {
                    continue;
                }
                let (pa, pb) = (vpos.get(a as usize).copied().unwrap_or_default(), vpos.get(b as usize).copied().unwrap_or_default());
                let c = edge_c(edge_key(a, b));
                let mut len = 0.0;
                let mut line = Vec::with_capacity(17);
                line.push(pa);
                for q in 1..=16 {
                    let x = if q == 16 { pb } else { place(g, &targets, c, lerp(pa, pb, q as f64 / 16.0), tol) };
                    len += dist(*line.last().unwrap_or(&pa), x);
                    line.push(x);
                }
                // The edge's length in element sizes.
                let measured = match &opts.field {
                    Some(f) if len > tol => crate::size::ideal_divisions(f, &line)?,
                    _ => len / opts.size,
                };
                if let Some(m) = class_max.get_mut(cls[d]) {
                    *m = m.max(measured);
                }
                if len - dist(pa, pb) > tol
                    && let Some(l) = class_bent.get_mut(cls[d])
                {
                    l.0 += libm::log(measured);
                    l.1 += 1;
                }
                if let Some(l) = class_log.get_mut(cls[d]).filter(|_| len > tol) {
                    l.0 += libm::log(measured);
                    l.1 += 1;
                }
            }
        }
    }
    // Curved targets: across each outer block face, along each of its two directions, the
    // largest normal curvature of the target in that direction (at a grid of points projected
    // onto it, from the turn of the B-rep face's normal over a short step on that face) times
    // the face's length there, over the largest turn an element may span, is the least count
    // of that direction's class.
    let probe = 1e-3 * opts.size;
    let mut class_arc = vec![0usize; nclass];
    if opts.max_arc_angle > 0.0 {
        for &(b, d, s) in &outer {
            let Some(block) = t.blocks.get(b) else { continue };
            let Some(&tg) = face_target.get(&face_key(block, d, s)) else { continue };
            let Some(set) = targets.get(tg) else { continue };
            let c: [P; 8] = std::array::from_fn(|k| vpos.get(block[k] as usize).copied().unwrap_or_default());
            let (e1, e2) = other(d);
            let m = CURVATURE_SAMPLES;
            let mut grid: Vec<Option<(P, P, u32)>> = Vec::with_capacity((m + 1) * (m + 1));
            for j in 0..=m {
                for i in 0..=m {
                    let mut u = [0.0; 3];
                    u[d] = s as f64;
                    u[e1] = i as f64 / m as f64;
                    u[e2] = j as f64 / m as f64;
                    grid.push(target::on_target(g, set, trilinear(&c, u)));
                }
            }
            let at = |i: usize, j: usize| grid.get(j * (m + 1) + i).copied().flatten();
            for (dir, axis) in [(0, e1), (1, e2)] {
                let (mut kappa, mut length) = (0.0f64, 0.0f64);
                for line in 0..=m {
                    let mut run = 0.0;
                    for k in 0..m {
                        let (p, q) = if dir == 0 { (at(k, line), at(k + 1, line)) } else { (at(line, k), at(line, k + 1)) };
                        let (Some((pa, na, fa)), Some((pb, _, _))) = (p, q) else { continue };
                        let step = dist(pa, pb);
                        run += step;
                        if step > tol {
                            let x = add(pa, scale(crate::vec3::sub(pb, pa), probe / step));
                            let Some(qx) = g.project(fa, x) else { continue };
                            let Some(nx) = g.normal(fa, qx) else { continue };
                            let h = dist(pa, qx);
                            if h > 0.5 * probe {
                                kappa = kappa.max(libm::acos(crate::vec3::dot(na, nx).clamp(-1.0, 1.0)) / h);
                            }
                        }
                    }
                    length = length.max(run);
                }
                let need = (kappa * length / opts.max_arc_angle).ceil();
                if let Some(slot) = classes.get(b).and_then(|cl| cl.get(axis)).and_then(|&k| class_arc.get_mut(k)) {
                    *slot = (*slot).max(need.min(crate::block::MAX_DIVISIONS as f64) as usize);
                }
            }
        }
    }
    let mut class_n: Vec<usize> = class_log
        .iter()
        .zip(&class_max)
        .zip(&class_arc)
        .zip(&class_bent)
        .map(|(((&(sum, k), &longest), &arc), &(bent_sum, bent_k))| {
            let mean = if k == 0 { 1.0 } else { libm::exp(sum / k as f64).round() };
            let bent = if bent_k == 0 { 0 } else { libm::exp(bent_sum / bent_k as f64).round().min(crate::block::MAX_DIVISIONS as f64) as usize };
            (mean.max((longest / MAX_STRETCH).ceil()) as usize).max(arc).max(bent).clamp(1, crate::block::MAX_DIVISIONS)
        })
        .collect();
    for &(b, d, n) in &opts.divisions {
        let cls = classes.get(b).and_then(|c| c.get(d)).ok_or_else(|| Error::Invalid(format!("no direction {d} of block {b}")))?;
        if !(1..=crate::block::MAX_DIVISIONS).contains(&n) {
            return Err(Error::Invalid(format!("{n} divisions are out of range")));
        }
        if let Some(slot) = class_n.get_mut(*cls) {
            *slot = n;
        }
    }
    for (class, n) in fixed::class_divisions(t, &classes, &fixed)? {
        if let Some(slot) = class_n.get_mut(class) {
            if opts.divisions.iter().any(|&(b, d, _)| classes.get(b).and_then(|c| c.get(d)) == Some(&class)) && *slot != n {
                return Err(Error::Invalid(format!("divisions set by hand ({}) differ from the fixed patch's ({n})", *slot)));
            }
            *slot = n;
        }
    }
    let divisions: Vec<[usize; 3]> = classes.iter().map(|c| c.map(|k| class_n.get(k).copied().unwrap_or(1))).collect();
    let lat = build(t, &divisions, opts.max_elements)?;

    // Each node's constraint.
    let constraint: Vec<Constraint> = lat
        .entity
        .iter()
        .map(|e| match e {
            Entity::Vertex(v) => vertex_c(*v),
            Entity::Edge(k) => edge_c(*k),
            Entity::Face(k) => face_target.get(k).map_or(Constraint::Free, |&tg| Constraint::Surface(tg)),
            Entity::Inside => Constraint::Free,
        })
        .collect();
    let n = lat.entity.len();
    let mut nodes = vec![[0.0; 3]; n];
    let mut placed = vec![false; n];
    // Nodes on fixed patches: placed and held.
    let mut held = vec![false; n];
    for ((b, d, s), grid) in &fixed.faces {
        let (Some(nd), Some(bn)) = (lat.divisions.get(*b), lat.block_nodes.get(*b)) else { continue };
        let (e1, e2) = other(*d);
        for j in 0..=nd[e2] {
            for i in 0..=nd[e1] {
                let mut at = [0; 3];
                at[*d] = s * nd[*d];
                at[e1] = i;
                at[e2] = j;
                let (Some(&id), Some(&p)) = (bn.get(crate::quad::lattice_index(*nd, at)), grid.get(i + (nd[e1] + 1) * j)) else {
                    return Err(Error::Failed("a fixed patch's grid does not fit its block face".into()));
                };
                if let (Some(slot), Some(h), Some(pl)) = (nodes.get_mut(id as usize), held.get_mut(id as usize), placed.get_mut(id as usize)) {
                    *slot = p;
                    *h = true;
                    *pl = true;
                }
            }
        }
    }
    for (i, e) in lat.entity.iter().enumerate() {
        if let (Entity::Vertex(v), false) = (e, held[i]) {
            nodes[i] = vpos.get(*v as usize).copied().unwrap_or_default();
            placed[i] = true;
        }
    }
    // 3. Edge nodes: even by arc length along their curve or surface.
    fieldcraft_model::progress::stage("Projection: edge nodes", None, "");
    fieldcraft_model::progress::check()?;
    for (b, block) in t.blocks.iter().enumerate() {
        let nd = lat.divisions[b];
        let bn = &lat.block_nodes[b];
        for d in 0..3 {
            for (k, bits) in CORNERS.iter().enumerate() {
                if bits[d] != 0 {
                    continue;
                }
                let mut e = *bits;
                e[d] = 1;
                let c = edge_c(edge_key(block[k], block[corner_at(e)]));
                let ids: Vec<u32> = (0..=nd[d])
                    .map(|q| {
                        let mut at = [bits[0] * nd[0], bits[1] * nd[1], bits[2] * nd[2]];
                        at[d] = q;
                        bn.get(crate::quad::lattice_index(nd, at)).copied().unwrap_or(0)
                    })
                    .collect();
                if ids.iter().skip(1).take(nd[d].saturating_sub(1)).all(|&i| placed[i as usize]) {
                    continue;
                }
                let (pa, pb) = (nodes[ids[0] as usize], nodes[ids[nd[d]] as usize]);
                let m = nd[d];
                let mut line: Vec<P> = (0..=m).map(|q| place(g, &targets, c, lerp(pa, pb, q as f64 / m as f64), tol)).collect();
                line[0] = pa;
                line[m] = pb;
                if c != Constraint::Free {
                    for _ in 0..EDGE_ROUNDS {
                        line = even(&line).into_iter().map(|x| place(g, &targets, c, x, tol)).collect();
                        line[0] = pa;
                        line[m] = pb;
                    }
                }
                for (q, &id) in ids.iter().enumerate() {
                    if !placed[id as usize] {
                        nodes[id as usize] = line[q];
                        placed[id as usize] = true;
                    }
                }
            }
        }
    }
    // 4. Faces (projected when outside), then 5. the interior.
    fieldcraft_model::progress::stage("Projection: faces and interior", None, "");
    fieldcraft_model::progress::check()?;
    for (b, block) in t.blocks.iter().enumerate() {
        let nd = lat.divisions[b];
        let bn = &lat.block_nodes[b];
        let mut x: Vec<Option<P>> = bn.iter().map(|&i| placed[i as usize].then(|| nodes[i as usize])).collect();
        crate::quad::transfinite_faces(&mut x, nd, None, |d, side, grid| {
            if let Some(&tg) = face_target.get(&face_key(block, d, side)) {
                for p in grid.iter_mut() {
                    *p = place(g, &targets, Constraint::Surface(tg), *p, tol);
                }
            }
            Ok(())
        })?;
        for (slot, &i) in x.iter_mut().zip(bn) {
            if placed[i as usize] {
                *slot = Some(nodes[i as usize]);
            } else if let Some(p) = *slot {
                nodes[i as usize] = p;
                placed[i as usize] = true;
            }
        }
        let blend = nd.map(|m| (0..=m).map(|q| q as f64 / m as f64).collect::<Vec<f64>>());
        crate::quad::transfinite_volume(&mut x, nd, [&blend[0], &blend[1], &blend[2]])?;
        for (slot, &i) in x.iter().zip(bn) {
            if let (false, Some(p)) = (placed[i as usize], slot) {
                nodes[i as usize] = *p;
                placed[i as usize] = true;
            }
        }
    }
    if placed.iter().any(|p| !p) {
        return Err(Error::Failed("a node of the block lattice was not placed".into()));
    }
    let corners_of = |nodes: &[P], h: &[u32; 8]| h.map(|v| nodes.get(v as usize).copied().unwrap_or_default());
    let initial_min = smooth::scaled_jacobians(&nodes, &lat.hexes).into_iter().fold(f64::INFINITY, f64::min);

    // The outer quads (outward, hex8 face order) and the target each lies on.
    let hex8_faces = ElementKind::Hex8.faces();
    let mut quads: Vec<([u32; 4], u32, u8, usize)> = Vec::new();
    let mut sheet_quads: Vec<[u32; 4]> = Vec::new();
    for (e, (hex, (b, cell))) in lat.hexes.iter().zip(&lat.hex_cell).enumerate() {
        let block = t.blocks.get(*b as usize).copied().unwrap_or([0; 8]);
        let nd = lat.divisions.get(*b as usize).copied().unwrap_or([1; 3]);
        for (fi, (d, s)) in FACE_SIDES.iter().copied().enumerate() {
            if cell[d] != if s == 0 { 0 } else { nd[d] - 1 } {
                continue;
            }
            let Some(&tg) = face_target.get(&face_key(&block, d, s)) else { continue };
            let local = hex8_faces.get(fi).copied().unwrap_or(&[]);
            let q = [0, 1, 2, 3].map(|i| local.get(i).and_then(|&l| hex.get(l)).copied().unwrap_or(0));
            if sheet_targets.contains(&tg) {
                sheet_quads.push(q);
            } else {
                quads.push((q, e as u32, fi as u8 + 1, tg));
            }
        }
    }
    // Nodes on a sheet slide along it among their neighbours on it, as boundary nodes do.
    let boundary_quads: Vec<[u32; 4]> = quads.iter().map(|q| q.0).chain(sheet_quads).collect();
    fieldcraft_model::progress::stage("Projection: smoothing", None, "");
    fieldcraft_model::progress::check()?;
    let smoothing = smooth::smooth(
        smooth::Mesh { nodes: &mut nodes, hexes: &lat.hexes, constraint: &constraint, held: &held, boundary: &boundary_quads },
        g,
        &targets,
        smooth::Run { method: opts.smoothing, iterations: opts.iterations, tol: 1e-4 * opts.size, order },
    )?;

    fieldcraft_model::progress::check()?;
    // Certified quality, the gap to the geometry and the tags.
    let (mut min_sj, mut sum_sj, mut volume) = (f64::INFINITY, 0.0, 0.0);
    let mut invalid: Vec<usize> = Vec::new();
    let qualities = fieldcraft_parallel::map_grained(lat.hexes.len(), smooth::HEX_GRAIN, |e| {
        lat.hexes.get(e).map_or_else(
            || Err(Error::Failed("an element of the lattice is missing".into())),
            |h| element_quality(ElementKind::Hex8, &corners_of(&nodes, h)),
        )
    });
    for (e, q) in qualities.into_iter().enumerate() {
        let q = q?;
        let sj = q.get(Metric::ScaledJacobian).unwrap_or(-1.0);
        if q.validity().is_none_or(|v| v != crate::quality::Validity::Valid) {
            invalid.push(e);
        }
        min_sj = min_sj.min(sj);
        sum_sj += sj;
        volume += q.size;
    }
    let mut max_gap = 0.0f64;
    for (i, c) in constraint.iter().enumerate() {
        if *c != Constraint::Free {
            max_gap = max_gap.max(dist(nodes[i], place(g, &targets, *c, nodes[i], tol)));
        }
    }
    let mut boundary = Vec::with_capacity(quads.len());
    let mut centres = Vec::with_capacity(quads.len());
    let mut enclosed = 0.0;
    for (q, e, f, tg) in &quads {
        let pts = q.map(|v| nodes.get(v as usize).copied().unwrap_or_default());
        enclosed += crate::hex::face_flux(&pts);
        let centre = scale(pts.iter().fold([0.0; 3], |a, p| add(a, *p)), 0.25);
        let set = targets.get(*tg).map_or(&[][..], Vec::as_slice);
        let brep_face = on_target(g, set, centre).map_or_else(|| set.first().copied().unwrap_or(0), |h| h.2);
        boundary.push(BoundaryFace { element: *e, face: *f, brep_face });
        let shortest = (0..4).map(|i| dist(pts[i], pts[(i + 1) % 4])).fold(f64::INFINITY, f64::min);
        centres.push((centre, shortest, *tg));
    }
    claim_narrow_faces(&mut boundary, &centres, &targets, g);
    let surface_volume = surface.volume();
    let report = ProjectReport {
        blocks: t.blocks.len(),
        elements: lat.hexes.len(),
        initial_min_scaled_jacobian: initial_min,
        min_scaled_jacobian: min_sj,
        mean_scaled_jacobian: sum_sj / lat.hexes.len().max(1) as f64,
        smoothing,
        max_gap,
        volume,
        surface_volume,
        targets: targets.clone(),
        divisions: lat.divisions.clone(),
        fixed_nodes: held.iter().filter(|h| **h).count(),
    };
    // Where the worst element is, for the error messages.
    let worst = || {
        lat.hexes
            .iter()
            .zip(&lat.hex_cell)
            .min_by(|a, b| smooth::hex_scaled_jacobian(&corners_of(&nodes, a.0)).total_cmp(&smooth::hex_scaled_jacobian(&corners_of(&nodes, b.0))))
            .map_or_else(Message::default, |(h, c)| {
                let x = scale(corners_of(&nodes, h).iter().fold([0.0; 3], |a, p| add(a, *p)), 0.125);
                crate::message![format!("worst in block {} near ", c.0), crate::message::point(x)]
            })
    };
    if !invalid.is_empty() {
        return Err(Error::failed(
            crate::message![format!(
                "{} of {} elements are inverted, twisted or degenerate after projection (certified Jacobian; smallest corner scaled Jacobian {min_sj:.3}): ",
                invalid.len(),
                lat.hexes.len()
            )]
            .then(where_invalid(t, &lat, &nodes, &invalid))
            .then(crate::message![
                ". The blocks do not fit the part there — add an O-grid or a template where a block corner lands on a curved face, or split the blocks so the faces meeting at that junction are not squeezed"
            ]),
        ));
    }
    if min_sj < opts.min_quality {
        return Err(Error::failed(
            crate::message![format!("the smallest scaled Jacobian is {min_sj:.3}, below {} (", opts.min_quality)].then(worst()).then(
                crate::message![
                    "): a block corner flattens on a curved face — use an O-grid or a butterfly/sphere template there, or lower `min_quality`"
                ],
            ),
        ));
    }
    if (enclosed - volume).abs() > 1e-6 * volume.abs().max(1e-300) {
        return Err(Error::failed(crate::message![
            "the elements overlap or leave a gap inside (their volume ",
            crate::message::volume(volume),
            " differs from the ",
            crate::message::volume(enclosed),
            " their outer faces enclose): check that block faces meet face to face"
        ]));
    }
    if (volume - surface_volume).abs() > 0.2 * surface_volume.abs() {
        return Err(Error::failed(crate::message![
            "the mesh volume ",
            crate::message::volume(volume),
            " is far from the body's ",
            crate::message::volume(surface_volume),
            ": the blocks do not cover the part (every outer block face must stand for part of its surface)"
        ]));
    }
    Ok(ProjectedMesh { nodes, hexes: lat.hexes, boundary, report })
}

/// Most invalid elements an error names.
const NAMED_INVALID: usize = 4;

/// Where the invalid elements `invalid` are: per element its block, cell and centre, and when it
/// touches a block vertex (a junction), the vertex and the blocks that meet there.
fn where_invalid(t: &Topology, lat: &lattice::Lattice, nodes: &[P], invalid: &[usize]) -> Message {
    let mut parts: Vec<Message> = Vec::new();
    for &e in invalid.iter().take(NAMED_INVALID) {
        let (Some(h), Some((b, cell))) = (lat.hexes.get(e), lat.hex_cell.get(e)) else { continue };
        let x = scale(h.iter().fold([0.0; 3], |a, &v| add(a, nodes.get(v as usize).copied().unwrap_or_default())), 0.125);
        let mut s = crate::message![format!("block {b} cell {cell:?} near "), crate::message::point(x)];
        let vertices: BTreeSet<u32> = h
            .iter()
            .filter_map(|&v| match lat.entity.get(v as usize) {
                Some(Entity::Vertex(k)) => Some(*k),
                _ => None,
            })
            .collect();
        for k in vertices {
            let around: Vec<String> = t.blocks.iter().enumerate().filter(|(_, blk)| blk.contains(&k)).map(|(i, _)| i.to_string()).collect();
            s.0.push(format!(" at block vertex {k} where blocks {} meet", around.join(", ")).into());
        }
        parts.push(s);
    }
    if invalid.len() > NAMED_INVALID {
        parts.push(crate::message![format!("and {} more", invalid.len() - NAMED_INVALID)]);
    }
    Message::join(parts, "; ")
}

/// The polyline's points redistributed evenly by arc length (ends kept).
fn even(line: &[P]) -> Vec<P> {
    let m = line.len().saturating_sub(1);
    let mut acc = vec![0.0; line.len()];
    for q in 1..line.len() {
        acc[q] = acc[q - 1] + dist(line[q - 1], line[q]);
    }
    let total = acc.last().copied().unwrap_or(0.0);
    if !(total > 0.0) || m == 0 {
        return line.to_vec();
    }
    let mut out = Vec::with_capacity(line.len());
    let mut seg = 0;
    for q in 0..=m {
        let s = total * q as f64 / m as f64;
        while seg + 1 < m && acc[seg + 1] < s {
            seg += 1;
        }
        let l = acc[seg + 1] - acc[seg];
        let t = if l > 0.0 { ((s - acc[seg]) / l).clamp(0.0, 1.0) } else { 0.0 };
        out.push(lerp(line[seg], line[seg + 1], t));
    }
    out
}
