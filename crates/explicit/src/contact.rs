//! Penalty contact between element faces, with Coulomb friction.
//!
//! **Kinematics: segment to segment, both ways.** The normal contact pressure is integrated over
//! the overlap of each pair of facing faces: a secondary face and a main face it faces are
//! projected onto the secondary face's plane, each split into triangles, clipped against each
//! other (Sutherland–Hodgman), and the overlap is integrated with a three-point rule per
//! sub-triangle; at each point the gap g = (x_s − x_m) · n_m between the two faces (n_m the main
//! face's outward normal) gives a penalty pressure p = ε (−g) where they overlap, which acts on
//! the secondary face's nodes with weights N_i^s and on the main face's with −N_j^m. This is the
//! segmentation and integration of mortar contact (M. A. Puso, T. A. Laursen, *A mortar
//! segment-to-segment contact method for large deformation solid mechanics*, CMAME
//! 193:601–629, 2004) with the penalty applied at the integration points (K. A. Fischer, P.
//! Wriggers, *Frictionless 2D contact formulations for finite deformations based on the mortar
//! method*, Comput. Mech. 36:226–244, 2005): for a uniform pressure across flat faces both sides
//! receive exactly the consistent nodal forces ∫ p N dA, whatever the meshes, so the contact
//! patch test is passed — node-to-segment contact fails it on non-matching meshes (R. L. Taylor,
//! P. Papadopoulos, *On a patch test for contact problems in two dimensions*, 1991). Every
//! contact runs both ways — a pair with each surface as the secondary in turn, single-surface
//! (general and self) contact with every face against every other one that shares no node with
//! it — at half the penalty each, so neither side is privileged and no corner of either surface
//! sinks into the other. Two faces are in contact where they face each other (their normals more
//! than [`FACING`] from parallel, which keeps the faces of one body that meet at an edge apart)
//! and the gap is negative but no deeper than half the main face's element is thick (deeper, the
//! surfaces have passed through each other). Faces of eroded elements leave the contact.
//!
//! **Search.** Every increment: each surface's faces go into a bucket grid with their boxes grown
//! by their depth limit ([`fieldcraft_fem::contact::Grid`]: D. J. Benson, J. O. Hallquist, *A
//! single surface contact algorithm for the post-buckling analysis of shell structures*, CMAME
//! 78:141–163, 1990); each secondary face takes the main faces listed in the cells its box
//! overlaps. The faces are integrated in parallel, each on its own, and their forces added face by
//! face in a fixed order, so the result does not depend on how the faces are split among threads.
//!
//! **Penalty.** The penalty per area comes from the masses and the stable increment Δt:
//! ε = s γ min_i (m_i / A_i) / Δt², with m_i the lumped mass of node i and A_i the area that
//! bounds its contact stiffness, over the nodes of the contact. The stiffness entries of node i
//! in direction d sum to at most 2 ε Σ_f a_f |n_f,d| Σ_d′ |n_f,d′| over its faces f (a_f its
//! share of the face, a quarter of a quadrilateral or a third of a triangle, n_f the unit normal;
//! its weights over a face sum to one on both sides), so A_i = κ_i max_d Σ_f a_f |n_f,d| with
//! κ_i = max_f Σ_d |n_f,d| ∈ [1, √3]: the node's area on a flat surface in any orientation, and
//! no more than the largest of its areas facing one way at an edge or a corner. By Gershgorin's
//! theorem the contact then adds at most 2 s γ / Δt² to the square of any frequency (Weyl's
//! inequality adds it to the elements'); with γ = [`GAMMA`] and the
//! default increment scale (0.9 of the critical increment, ω Δt ≤ 1.8) the step stays under the
//! stability limit ω Δt = 2 (1.8² + 2γ ≤ 4) and contact never shortens the increment. The penalty
//! scale s ≤ 1 only softens it. The pressure stiffness this gives is about a fifth of the
//! smallest element's M/h (M the constrained modulus): stiffer than the element-based penalty
//! f K A²/V with its usual f = 0.1 (J. O. Hallquist, G. L. Goudreau, D. J. Benson, *Sliding
//! interfaces with contact-impact in large-scale Lagrangian computations*, CMAME 51:107–137,
//! 1985), and a uniform ε keeps a uniform pressure uniform.
//!
//! **Friction** at the secondary nodes by a regularised return map (P. Wriggers, *Computational
//! Contact Mechanics*, 2nd ed., Springer 2006, §5.2 and §10.3; J. C. Simo, T. A. Laursen, *An
//! augmented Lagrangian treatment of contact problems involving friction*, Comput. Struct.
//! 42(1):97–116, 1992): a node's normal force f_n is its share of the pressure, its tangential
//! force an elastic predictor, the previous force moved onto the current tangent plane less
//! k_t times the tangential motion of the node relative to the closest point of the main face that
//! pressed it most over the increment, t* = P t_n − k_t Δu_T with k_t the node's normal stiffness ε A_i / 2; if
//! |t*| ≤ μ f_n the node sticks (t = t*), otherwise it slips and the force returns to the cone,
//! t = μ f_n t*/|t*|, dissipating μ f_n (|t*| − μ f_n)/k_t. The reaction goes to the main face's
//! nodes with −N_j. The coefficient μ = μ_k + (μ_s − μ_k) e^{−d γ̇} follows the slip rate
//! γ̇ = |Δu_T|/Δt ([`fieldcraft_model::Friction`]). A node that leaves contact forgets its
//! tangential force.
//!
//! **Energy.** The contact forces enter the equations of motion like external forces; their work
//! (the penalty energy stored, plus what friction dissipates) is integrated with the trapezoidal
//! rule like the internal forces' (see [`crate::integrate`]), and the friction dissipation is
//! reported too.
//!
//! **Initial penetration.** Before the analysis each node of a secondary surface is projected
//! onto the faces of the main surface ([`fieldcraft_fem::contact::closest_point`]); a node behind
//! a face it lies over (within a fiftieth of the face's size), within the depth limit and with its
//! own faces facing it, is moved onto it (strain-free: the reference coordinates change before the
//! elements are built) or the analysis is refused with the node, the face and the depth, as the
//! contact's `initial_penetration` says.
//!
//! **Erosion.** The faces of eroded elements leave the contact; the faces their removal exposes
//! do not join it (eroding contact is later work).

use fieldcraft_fem::contact::{Aabb, Grid, Projection, bounds, closest_point};
use fieldcraft_fem::frame::{cross, dot, norm, sub};
use fieldcraft_model::{Contact, ContactSurfaces, ElementFace, ElementKind, Friction, InitialPenetration, Model, Resolved};
use serde::Serialize;

use crate::elements::{State, add_at};
use crate::system::NODE_DOFS;
use crate::{Error, Result};

/// The penalty's share of the stability budget: ε = s γ min(m/A) / Δt² (see the module docs).
pub const GAMMA: f64 = 0.3;

/// A node lies over a face when its lateral offset is at most this fraction of the face's size
/// (the initial-penetration check).
pub const EDGE_TOLERANCE: f64 = 0.02;

/// A point deeper than this fraction of the thickness of a face's element has passed through it.
pub const DEPTH_FRACTION: f64 = 0.5;

/// Facing faces: the dot product of their unit normals is below minus this.
pub const FACING: f64 = 0.1;

/// A node behind a face by less than this fraction of the coordinates' size is touching it:
/// coordinates are exact to a few units in the last place.
const ROUND_OFF: f64 = 1e-9; // unit-factors: relative to the coordinates' size

/// Faces per parallel task.
const BLOCK: usize = 64;

/// The three-point rule on a triangle, exact for quadratics: barycentric points (of the second
/// and third corners), weights summing to one.
const TRIANGLE_RULE: [([f64; 2], f64); 3] =
    [([1.0 / 6.0, 1.0 / 6.0], 1.0 / 3.0), ([2.0 / 3.0, 1.0 / 6.0], 1.0 / 3.0), ([1.0 / 6.0, 2.0 / 3.0], 1.0 / 3.0)];

type Vec3 = [f64; 3];

/// A face of a contact surface.
#[derive(Debug, Clone)]
pub struct Face {
    /// The corner nodes (system indices), in the face's order; the last unused on a triangle.
    pub nodes: [usize; 4],
    pub count: usize,
    /// The element's label, its position in [`crate::system::System::elements`] and the face
    /// number.
    pub element: u32,
    pub index: usize,
    pub face: u8,
    /// The nodes of the face's element (a node of it is never in contact with the face).
    pub owner: Vec<usize>,
    /// The penetration beyond which a point has passed through: [`DEPTH_FRACTION`] of the
    /// element's thickness under the face.
    pub depth: f64,
}

impl Face {
    fn corners(&self) -> &[usize] {
        self.nodes.get(..self.count).unwrap_or_default()
    }
}

/// The faces and nodes of one surface.
#[derive(Debug, Clone, Default)]
pub struct Surface {
    pub faces: Vec<Face>,
    /// Its nodes, sorted (system indices).
    pub nodes: Vec<usize>,
    /// Per node (in `nodes` order), the faces that hold it.
    node_faces: Vec<Vec<usize>>,
    /// Per node, the area of the surface that bounds its contact stiffness (undeformed; see
    /// [`Interface::forces`]).
    area: Vec<f64>,
}

/// A node's friction state in one pass.
#[derive(Debug, Clone, Copy, Default)]
struct Slip {
    /// The tangential force on the node.
    force: Vec3,
}

/// One pass: the faces of a surface (secondary) against the faces of a surface (main).
#[derive(Debug, Clone)]
struct Pass {
    secondary: usize,
    main: usize,
    /// Per node of the secondary surface.
    slip: Vec<Slip>,
}

/// A contact as the solver runs it.
#[derive(Debug, Clone)]
pub struct Interface {
    pub name: String,
    pub friction: Friction,
    pub scale: f64,
    pub surfaces: Vec<Surface>,
    passes: Vec<Pass>,
    /// Nodes moved onto the other surface before the analysis.
    pub adjusted: Vec<Adjustment>,
    /// min over the contact's nodes of mass / area, once the masses are known.
    mass_per_area: Option<f64>,
}

/// A node moved onto a face before the analysis.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, serde::Deserialize)]
pub struct Adjustment {
    pub node: u32,
    pub element: u32,
    pub face: u8,
    /// How far it started inside.
    pub depth: f64,
}

/// What one evaluation of the contact forces found, per contact.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Evaluation {
    /// The normal force across the contact (the integrated pressure, both ways).
    pub normal_force: f64,
    /// The deepest penetration.
    pub penetration: f64,
    /// Integration points in contact.
    pub active: usize,
    /// Friction dissipated over the increment that led here.
    pub dissipation: f64,
}

fn invalid(m: impl Into<String>) -> Error {
    Error::Invalid(m.into())
}

/// The volume of a solid element of corners `x` (hex8 or tet4), for the thickness under a face.
fn volume(kind: ElementKind, x: &[Vec3]) -> Option<f64> {
    match kind {
        ElementKind::Hex8 => <[Vec3; 8]>::try_from(x).ok().and_then(|a| fieldcraft_fem::uniform::Hex8::hex8(&a).ok()).map(|k| k.volume),
        ElementKind::Tet4 => <[Vec3; 4]>::try_from(x).ok().and_then(|a| fieldcraft_fem::uniform::Tet4::tet4(&a).ok()).map(|k| k.volume),
        _ => None,
    }
}

/// A face's area vector (its outward normal times its area, flat or bilinear).
fn area_vector(x: &[Vec3]) -> Vec3 {
    match x {
        [a, b, c] => cross(sub(*b, *a), sub(*c, *a)).map(|v| 0.5 * v),
        [a, b, c, d] => cross(sub(*c, *a), sub(*d, *b)).map(|v| 0.5 * v),
        _ => [0.0; 3],
    }
}

/// The largest side of a face's bounding box.
fn size(x: &[Vec3]) -> f64 {
    let b = bounds(x, 0.0);
    (0..3).map(|k| b[1][k] - b[0][k]).fold(0.0, f64::max)
}

/// Where a node is: its rest position plus its displacement.
fn position(xyz: &[Vec3], u: &[f64], node: usize) -> Vec3 {
    let x = xyz.get(node).copied().unwrap_or_default();
    [0, 1, 2].map(|k| x[k] + u.get(NODE_DOFS * node + k).copied().unwrap_or(0.0))
}

/// A node's displacement increment.
fn shift(du: &[f64], node: usize) -> Vec3 {
    [0, 1, 2].map(|k| du.get(NODE_DOFS * node + k).copied().unwrap_or(0.0))
}

impl Surface {
    /// The surface of element faces `faces`, its nodes, their faces and areas. `node` and
    /// `element` give a label's index among the analysed nodes and elements.
    fn new(
        model: &Model,
        faces: &[ElementFace],
        xyz: &[Vec3],
        node: &impl Fn(u32) -> Option<usize>,
        element: &impl Fn(u32) -> Option<usize>,
    ) -> Result<Surface> {
        let mut out = Surface::default();
        for f in faces {
            let at = |what: &str| invalid(format!("element {} face {}: {what}", f.element, f.face));
            let el = model.element(f.element).ok_or_else(|| at("the element is missing"))?;
            if !matches!(el.kind, ElementKind::Hex8 | ElementKind::Tet4) {
                return Err(Error::NotSupported(format!(
                    "contact on the faces of {} elements (element {}); explicit contact acts on hex8 and tet4 faces",
                    el.kind.name(),
                    f.element
                )));
            }
            let labels = model.face_nodes(*f).ok_or_else(|| at("no such face"))?;
            let corners: Vec<usize> = labels.iter().map(|l| node(*l).ok_or_else(|| at("a node is not analysed"))).collect::<Result<_>>()?;
            let owner: Vec<usize> = el.nodes.iter().filter_map(|l| node(*l)).collect();
            let index = element(f.element).ok_or_else(|| at("the element is not analysed"))?;
            let x: Vec<Vec3> = corners.iter().map(|n| xyz.get(*n).copied().unwrap_or_default()).collect();
            let ex: Vec<Vec3> = owner.iter().map(|n| xyz.get(*n).copied().unwrap_or_default()).collect();
            let area = norm(area_vector(&x));
            // Thickness V/A; a face of no area or an element of no volume falls back on its size.
            let thickness = volume(el.kind, &ex).filter(|v| *v > 0.0 && area > 0.0).map_or(size(&x), |v| v / area);
            let mut nodes = [usize::MAX; 4];
            for (slot, n) in nodes.iter_mut().zip(&corners) {
                *slot = *n;
            }
            out.faces.push(Face { nodes, count: corners.len(), element: f.element, index, face: f.face, owner, depth: DEPTH_FRACTION * thickness });
        }
        out.nodes = out.faces.iter().flat_map(|f| f.corners().iter().copied()).collect();
        out.nodes.sort_unstable();
        out.nodes.dedup();
        out.node_faces = vec![Vec::new(); out.nodes.len()];
        // Per node and direction, Σ_f A_f |n_f,d| / count, and the largest Σ_d |n_d| of its faces.
        let mut directional = vec![[0.0; 3]; out.nodes.len()];
        let mut spread = vec![1.0f64; out.nodes.len()];
        for (i, f) in out.faces.iter().enumerate() {
            let x: Vec<Vec3> = f.corners().iter().map(|n| xyz.get(*n).copied().unwrap_or_default()).collect();
            let a = area_vector(&x);
            let l = norm(a);
            let share = a.map(|c| c.abs() / f.count as f64);
            let kappa = if l > 0.0 { a.iter().map(|c| c.abs()).sum::<f64>() / l } else { 1.0 };
            for n in f.corners() {
                if let Ok(k) = out.nodes.binary_search(n) {
                    if let Some(list) = out.node_faces.get_mut(k) {
                        list.push(i);
                    }
                    if let (Some(d), Some(s)) = (directional.get_mut(k), spread.get_mut(k)) {
                        for c in 0..3 {
                            d[c] += share[c];
                        }
                        *s = s.max(kappa);
                    }
                }
            }
        }
        out.area = directional.iter().zip(&spread).map(|(d, s)| s * d.iter().copied().fold(0.0, f64::max)).collect();
        Ok(out)
    }

    /// The faces' corner positions at displacements `u`.
    fn positions(&self, xyz: &[Vec3], u: &[f64]) -> Vec<Vec<Vec3>> {
        self.faces.iter().map(|f| f.corners().iter().map(|n| position(xyz, u, *n)).collect()).collect()
    }

    /// The nodes' outward normals (the sum of their faces' area vectors; not normalised).
    fn node_normals(&self, faces: &[Vec<Vec3>]) -> Vec<Vec3> {
        self.node_faces
            .iter()
            .map(|list| {
                let mut n = [0.0; 3];
                for a in list.iter().filter_map(|i| faces.get(*i)).map(|x| area_vector(x)) {
                    for k in 0..3 {
                        n[k] += a[k];
                    }
                }
                n
            })
            .collect()
    }
}

/// The faces' grid for the global search at positions `positions`.
fn grid(surface: &Surface, positions: &[Vec<Vec3>]) -> Grid {
    let boxes: Vec<Aabb> = surface.faces.iter().zip(positions).map(|(f, x)| bounds(x, f.depth)).collect();
    Grid::new(&boxes, 4 * boxes.len() + 64)
}

/// The face a node at `x` (index `node`, outward normal `normal`) starts inside, among the grid's
/// candidates: over it, behind it within its depth limit, facing it; the nearest.
fn inside(surface: &Surface, positions: &[Vec<Vec3>], grid: &Grid, node: usize, x: Vec3, normal: Vec3) -> Option<(usize, Projection)> {
    let mut best: Option<(usize, Projection)> = None;
    for &c in grid.candidates(x) {
        let c = c as usize;
        let (Some(face), Some(fx)) = (surface.faces.get(c), positions.get(c)) else { continue };
        if face.owner.contains(&node) {
            continue;
        }
        let Some(q) = closest_point(fx, x) else { continue };
        let over = q.lateral <= EDGE_TOLERANCE * size(fx);
        let n = norm(normal);
        let facing = n == 0.0 || dot(normal, q.normal) < -FACING * n;
        if q.gap < 0.0 && -q.gap <= face.depth && over && facing && best.is_none_or(|b| q.distance < b.1.distance) {
            best = Some((c, q));
        }
    }
    best
}

/// A point of a plane in its own 2D coordinates.
type Point = [f64; 2];

fn cross2(o: Point, a: Point, b: Point) -> f64 {
    (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
}

/// The signed area of a polygon (counter-clockwise positive).
fn signed_area(p: &[Point]) -> f64 {
    let mut s = 0.0;
    for (i, a) in p.iter().enumerate() {
        let b = p.get((i + 1) % p.len()).copied().unwrap_or(*a);
        s += a[0] * b[1] - a[1] * b[0];
    }
    0.5 * s
}

/// `subject` clipped by the convex counter-clockwise polygon `clip` (Sutherland–Hodgman).
fn clip(subject: &[Point], clip: &[Point]) -> Vec<Point> {
    let mut out: Vec<Point> = subject.to_vec();
    for (i, a) in clip.iter().enumerate() {
        let Some(b) = clip.get((i + 1) % clip.len()).copied() else { break };
        let input = std::mem::take(&mut out);
        let side = |p: Point| cross2(*a, b, p);
        for (j, p) in input.iter().enumerate() {
            let Some(q) = input.get((j + 1) % input.len()).copied() else { break };
            let (sp, sq) = (side(*p), side(q));
            if sp >= 0.0 {
                out.push(*p);
            }
            if (sp >= 0.0) != (sq >= 0.0) {
                let t = sp / (sp - sq);
                out.push([p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])]);
            }
        }
        if out.is_empty() {
            break;
        }
    }
    out
}

/// The corners `t` of `points`, counter-clockwise.
fn counter_clockwise(t: &[usize; 3], points: &[Point]) -> Vec<Point> {
    let mut out: Vec<Point> = t.iter().filter_map(|i| points.get(*i).copied()).collect();
    if signed_area(&out) < 0.0 {
        out.reverse();
    }
    out
}

/// The triangles a face splits into (corner positions in its order).
fn triangles(count: usize) -> &'static [[usize; 3]] {
    if count == 3 { &[[0, 1, 2]] } else { &[[0, 1, 2], [0, 2, 3]] }
}

/// What the integration needs of each main face, worked out once per increment.
#[derive(Debug, Clone, Copy)]
struct FaceInfo {
    /// Area vector and its length.
    area: Vec3,
    length: f64,
    /// The lateral offset beyond which a point is past the face's edge.
    edge: f64,
    alive: bool,
}

fn face_info(surface: &Surface, positions: &[Vec<Vec3>], alive: &dyn Fn(&Face) -> bool) -> Vec<FaceInfo> {
    surface
        .faces
        .iter()
        .zip(positions)
        .map(|(f, x)| {
            let area = area_vector(x);
            FaceInfo { area, length: norm(area), edge: EDGE_TOLERANCE * size(x), alive: alive(f) }
        })
        .collect()
}

/// One secondary face against the main surface: its contributions to the nodal forces and to
/// each secondary node's normal force.
#[derive(Debug, Default)]
struct FaceForces {
    /// (node, force) contributions, in order.
    forces: Vec<(usize, Vec3)>,
    /// (secondary node, normal force share, main face) of each integration point.
    normal: Vec<(usize, f64, usize)>,
    penetration: f64,
    points: usize,
    total: f64,
}

/// The penalty pressure integrated over the overlaps of secondary face `f` (corners `fx`) with
/// the main faces facing it, at penalty `eps` per area.
#[allow(clippy::too_many_arguments)]
fn integrate_face(
    f: &Face,
    fx: &[Vec3],
    alive: bool,
    main: &Surface,
    mx: &[Vec<Vec3>],
    info: &[FaceInfo],
    grid: &Grid,
    eps: f64,
    found: &mut Vec<u32>,
) -> FaceForces {
    let mut out = FaceForces::default();
    if !alive {
        return out;
    }
    let area = area_vector(fx);
    let a = norm(area);
    if !(a > 0.0 && a.is_finite()) {
        return out;
    }
    let normal = area.map(|c| c / a);
    let n = fx.len() as f64;
    let center: Vec3 = [0, 1, 2].map(|k| fx.iter().map(|p| p[k]).sum::<f64>() / n);
    let Some(first) = fx.first().copied() else { return out };
    // A basis of the face's plane, (e1, e2, normal) right-handed.
    let r = sub(first, center);
    let r = sub(r, normal.map(|c| c * dot(r, normal)));
    let rl = norm(r);
    if !(rl > 0.0) {
        return out;
    }
    let e1 = r.map(|c| c / rl);
    let e2 = cross(normal, e1);
    let flat = |x: Vec3| -> Point {
        let d = sub(x, center);
        [dot(d, e1), dot(d, e2)]
    };
    let lift = |p: Point| -> Vec3 { [0, 1, 2].map(|k| center[k] + p[0] * e1[k] + p[1] * e2[k]) };
    let own: Vec<Point> = fx.iter().map(|x| flat(*x)).collect();
    let own_triangles: Vec<Vec<Point>> = triangles(f.count).iter().map(|t| counter_clockwise(t, &own)).collect();
    grid.candidates_in(&bounds(fx, 0.0), found);
    for &c in found.iter() {
        let c = c as usize;
        let (Some(m), Some(x), Some(mi)) = (main.faces.get(c), mx.get(c), info.get(c)) else { continue };
        if !mi.alive || m.corners().iter().any(|k| f.corners().contains(k)) {
            continue;
        }
        if !(mi.length > 0.0) || dot(mi.area, normal) >= -FACING * mi.length {
            continue;
        }
        let theirs: Vec<Point> = x.iter().map(|p| flat(*p)).collect();
        for s in &own_triangles {
            for tm in triangles(m.count) {
                let t = counter_clockwise(tm, &theirs);
                let overlap = clip(&t, s);
                if overlap.len() < 3 {
                    continue;
                }
                let Some(o) = overlap.first().copied() else { continue };
                for w in overlap.windows(2).skip(1) {
                    let (p, q) = (w[0], w[1]);
                    let tri_area = 0.5 * cross2(o, p, q);
                    if !(tri_area > 0.0) {
                        continue;
                    }
                    for (bary, weight) in TRIANGLE_RULE {
                        let y = [0, 1].map(|k| o[k] + bary[0] * (p[k] - o[k]) + bary[1] * (q[k] - o[k]));
                        let point = lift(y);
                        let (Some(qs), Some(qm)) = (closest_point(fx, point), closest_point(x, point)) else { continue };
                        let gap = dot(sub(qs.point, qm.point), qm.normal);
                        // Over the main face (its foot not clamped to an edge), behind it, within
                        // its depth.
                        if !(gap < 0.0 && -gap <= m.depth && qm.lateral <= mi.edge) {
                            continue;
                        }
                        let push = eps * (-gap) * weight * tri_area;
                        if !push.is_finite() {
                            continue;
                        }
                        let force = qm.normal.map(|k| push * k);
                        for (node, w) in f.corners().iter().zip(qs.weights) {
                            out.forces.push((*node, force.map(|k| w * k)));
                            out.normal.push((*node, w * push, c));
                        }
                        for (node, w) in m.corners().iter().zip(qm.weights) {
                            out.forces.push((*node, force.map(|k| -w * k)));
                        }
                        out.penetration = out.penetration.max(-gap);
                        out.points += 1;
                        out.total += push;
                    }
                }
            }
        }
    }
    out
}

impl Interface {
    /// The contact `name` of the model on the analysed nodes and elements: its surfaces, and the
    /// initial penetrations adjusted in `xyz` (or refused).
    pub fn new(
        model: &Model,
        resolved: &Resolved,
        name: &str,
        contact: &Contact,
        xyz: &mut [Vec3],
        labels: &[u32],
        elements: &[u32],
    ) -> Result<Interface> {
        let context = |e: Error| match e {
            Error::Invalid(m) => invalid(format!("contact `{name}`: {m}")),
            other => other,
        };
        let node = |l: u32| labels.binary_search(&l).ok();
        let element = |l: u32| elements.binary_search(&l).ok();
        let faces_of = |set: &str| resolved.sets.get(set).map(|s| s.faces.clone()).unwrap_or_default();
        let (surfaces, passes) = match &contact.surfaces {
            ContactSurfaces::Pair { main, secondary } => {
                let main = Surface::new(model, &faces_of(main), xyz, &node, &element).map_err(context)?;
                let secondary = Surface::new(model, &faces_of(secondary), xyz, &node, &element).map_err(context)?;
                (vec![main, secondary], vec![(1, 0), (0, 1)])
            }
            ContactSurfaces::General { surface } => {
                let faces = match surface {
                    Some(s) => faces_of(s),
                    None => model.exterior_faces(),
                };
                (vec![Surface::new(model, &faces, xyz, &node, &element).map_err(context)?], vec![(0, 0)])
            }
        };
        let passes = passes
            .into_iter()
            .map(|(secondary, main)| Pass { secondary, main, slip: vec![Slip::default(); surfaces.get(secondary).map_or(0, |s| s.nodes.len())] })
            .collect();
        let mut interface = Interface {
            name: name.to_string(),
            friction: contact.friction,
            scale: contact.penalty_scale,
            surfaces,
            passes,
            adjusted: Vec::new(),
            mass_per_area: None,
        };
        interface.adjust(xyz, labels, contact.initial_penetration)?;
        Ok(interface)
    }

    /// Move the nodes that start inside a face onto it, or refuse them.
    fn adjust(&mut self, xyz: &mut [Vec3], labels: &[u32], mode: InitialPenetration) -> Result<()> {
        let zero: Vec<f64> = Vec::new();
        // The model's extent: the round-off of a coordinate is relative to the larger of it and
        // this (so the test reads the same in any unit of length).
        let extent = xyz.iter().flatten().fold(0.0f64, |m, c| m.max(c.abs()));
        for pass in &self.passes {
            let (Some(nodes), Some(faces)) = (self.surfaces.get(pass.secondary), self.surfaces.get(pass.main)) else { continue };
            let positions = faces.positions(xyz, &zero);
            let grid = grid(faces, &positions);
            let own = nodes.positions(xyz, &zero);
            let normals = nodes.node_normals(&own);
            // Each move is decided on the positions before any move of this pass.
            let mut moves: Vec<(usize, Vec3, Adjustment)> = Vec::new();
            for (i, &n) in nodes.nodes.iter().enumerate() {
                let x = xyz.get(n).copied().unwrap_or_default();
                let normal = normals.get(i).copied().unwrap_or_default();
                let Some((c, q)) = inside(faces, &positions, &grid, n, x, normal) else { continue };
                let Some(face) = faces.faces.get(c) else { continue };
                // Coordinates are exact to a few units in the last place: a node a hair behind a
                // face it was meshed on is touching, not inside.
                let scale = x.iter().fold(extent, |m, c| m.max(c.abs()));
                if -q.gap <= ROUND_OFF * scale {
                    continue;
                }
                let label = labels.get(n).copied().unwrap_or(0);
                let adjustment = Adjustment { node: label, element: face.element, face: face.face, depth: -q.gap };
                if mode == InitialPenetration::Refuse {
                    return Err(invalid(format!(
                        "contact `{}`: node {label} starts {:.4e} inside face {} of element {} (move it out, or let the contact adjust it: initial_penetration \"adjust\")",
                        self.name, -q.gap, face.face, face.element
                    )));
                }
                moves.push((n, q.point, adjustment));
            }
            for (n, point, a) in moves {
                if let Some(x) = xyz.get_mut(n) {
                    *x = point;
                }
                self.adjusted.push(a);
            }
        }
        Ok(())
    }

    /// The number of faces.
    pub fn faces(&self) -> usize {
        self.surfaces.iter().map(|s| s.faces.len()).sum()
    }

    /// Add the contact forces at displacements `u` into `force` (per degree of freedom), with
    /// `du` the displacement increment that led to `u` over `dt_prev`, the stable increment
    /// `dt` and the lumped masses per degree of freedom `mass`; update the friction state.
    #[allow(clippy::too_many_arguments)]
    pub fn forces(
        &mut self,
        xyz: &[Vec3],
        u: &[f64],
        du: &[f64],
        dt_prev: f64,
        dt: f64,
        mass: &[f64],
        state: &State,
        force: &mut [f64],
    ) -> Evaluation {
        let mut e = Evaluation::default();
        let node_mass = |n: usize| (0..3).map(|k| mass.get(NODE_DOFS * n + k).copied().unwrap_or(0.0)).fold(f64::INFINITY, f64::min);
        let mass_per_area = *self.mass_per_area.get_or_insert_with(|| {
            self.surfaces
                .iter()
                .flat_map(|s| s.nodes.iter().zip(&s.area))
                .filter(|(_, a)| **a > 0.0)
                .map(|(n, a)| node_mass(*n) / a)
                .fold(f64::INFINITY, f64::min)
        });
        // Every contact runs both ways, at half the penalty each.
        let eps = 0.5 * self.scale * GAMMA * mass_per_area / (dt * dt);
        if !(eps > 0.0 && eps.is_finite()) {
            for pass in &mut self.passes {
                pass.slip.fill(Slip::default());
            }
            return e;
        }
        let alive = |f: &Face| state.elements.get(f.index).is_none_or(|s| s.failed.is_none());
        let positions: Vec<Vec<Vec<Vec3>>> = self.surfaces.iter().map(|s| s.positions(xyz, u)).collect();
        let grids: Vec<Grid> = self.surfaces.iter().zip(&positions).map(|(s, p)| grid(s, p)).collect();
        let friction = self.friction;
        for pass in &mut self.passes {
            let (Some(secondary), Some(main), Some(sx), Some(mx), Some(grid)) = (
                self.surfaces.get(pass.secondary),
                self.surfaces.get(pass.main),
                positions.get(pass.secondary),
                positions.get(pass.main),
                grids.get(pass.main),
            ) else {
                continue;
            };
            let info = face_info(main, mx, &alive);
            let blocks = fieldcraft_parallel::map_grained(secondary.faces.len().div_ceil(BLOCK), 1, |b| {
                let mut found = Vec::new();
                (b * BLOCK..((b + 1) * BLOCK).min(secondary.faces.len()))
                    .map(|i| match (secondary.faces.get(i), sx.get(i)) {
                        (Some(f), Some(fx)) => integrate_face(f, fx, alive(f), main, mx, &info, grid, eps, &mut found),
                        _ => FaceForces::default(),
                    })
                    .collect::<Vec<_>>()
            });
            // Per secondary node: its normal force, and the main face that gave it the most.
            let mut normal = vec![0.0; secondary.nodes.len()];
            let mut partner: Vec<Option<(usize, f64)>> = vec![None; secondary.nodes.len()];
            for face in blocks.iter().flatten() {
                for (node, f) in &face.forces {
                    add_at(force, *node, f);
                }
                for (node, share, c) in &face.normal {
                    if let Ok(k) = secondary.nodes.binary_search(node)
                        && let (Some(slot), Some(best)) = (normal.get_mut(k), partner.get_mut(k))
                    {
                        *slot += share;
                        if best.is_none_or(|(_, s)| *share > s) {
                            *best = Some((*c, *share));
                        }
                    }
                }
                e.normal_force += face.total;
                e.penetration = e.penetration.max(face.penetration);
                e.active += face.points;
            }
            if friction.is_frictionless() {
                continue;
            }
            // Friction at the secondary nodes, in node order, against the point under the node
            // of the main face that pressed it most.
            for (k, (&n, slip)) in secondary.nodes.iter().zip(pass.slip.iter_mut()).enumerate() {
                let fn_ = normal.get(k).copied().unwrap_or(0.0);
                let kt = eps * secondary.area.get(k).copied().unwrap_or(0.0);
                let hit = partner
                    .get(k)
                    .copied()
                    .flatten()
                    .filter(|_| fn_ > 0.0 && kt > 0.0)
                    .and_then(|(c, _)| mx.get(c).and_then(|x| closest_point(x, position(xyz, u, n))).map(|q| (c, q)));
                let Some((c, q)) = hit else {
                    *slip = Slip::default();
                    continue;
                };
                let Some(face) = main.faces.get(c) else { continue };
                let mut rel = shift(du, n);
                for (m, w) in face.corners().iter().zip(q.weights) {
                    let d = shift(du, *m);
                    for c in 0..3 {
                        rel[c] -= w * d[c];
                    }
                }
                let normal = q.normal;
                let rn = dot(rel, normal);
                let slide = [0, 1, 2].map(|c| rel[c] - rn * normal[c]);
                let old = slip.force;
                let on = dot(old, normal);
                let trial: Vec3 = std::array::from_fn(|c| old[c] - on * normal[c] - kt * slide[c]);
                let size = norm(trial);
                let rate = if dt_prev > 0.0 { norm(slide) / dt_prev } else { 0.0 };
                let limit = friction.coefficient(rate) * fn_;
                let tangential = if size <= limit {
                    trial
                } else {
                    e.dissipation += limit * (size - limit) / kt;
                    trial.map(|c| c * limit / size)
                };
                slip.force = tangential;
                add_at(force, n, &tangential);
                for (m, w) in face.corners().iter().zip(q.weights) {
                    add_at(force, *m, &tangential.map(|c| -w * c));
                }
            }
        }
        e
    }
}

/// Every contact of the model, built on the analysed nodes and elements (with the initial
/// adjustments made in `xyz`), and no more than `max_faces` faces in all.
pub fn interfaces(
    model: &Model,
    resolved: &Resolved,
    xyz: &mut [Vec3],
    labels: &[u32],
    elements: &[u32],
    max_faces: usize,
) -> Result<Vec<Interface>> {
    let mut out = Vec::with_capacity(model.contacts().len());
    let mut faces = 0usize;
    for (name, c) in model.contacts() {
        let i = Interface::new(model, resolved, name, c, xyz, labels, elements)?;
        faces = faces.saturating_add(i.faces());
        if faces > max_faces {
            return Err(Error::Limit(format!(
                "the contacts hold more than {max_faces} faces (contact `{name}` passes the cap); name smaller surfaces"
            )));
        }
        out.push(i);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clipping_two_triangles_gives_their_overlap() {
        let a = [[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]];
        let b = [[1.0, -1.0], [1.0, 3.0], [-1.0, 1.0]];
        // b clockwise: reversed it is counter-clockwise.
        let mut b2 = b.to_vec();
        if signed_area(&b2) < 0.0 {
            b2.reverse();
        }
        let o = clip(&a, &b2);
        // The overlap: the part of a with x ≤ 1 and y ≤ x + 2… here the triangle's left half,
        // area 2 − ½ = 1.5.
        assert!((signed_area(&o) - 1.5).abs() < 1e-12, "{o:?} {}", signed_area(&o));
        assert!(clip(&a, &[[5.0, 5.0], [6.0, 5.0], [5.0, 6.0]]).is_empty());
        let same = clip(&a, &a);
        assert!((signed_area(&same) - 2.0).abs() < 1e-12);
    }
}
