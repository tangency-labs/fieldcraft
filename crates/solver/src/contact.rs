//! Contact in nonlinear static steps: mortar (segment-to-segment) normal contact enforced by an
//! augmented Lagrangian, with its consistent tangent.
//!
//! **Discretisation: mortar, one pass.** A contact pair's `secondary` surface is the non-mortar
//! side, its `main` surface the mortar side (M. A. Puso, T. A. Laursen, *A mortar
//! segment-to-segment contact method for large deformation solid mechanics*, CMAME 193:601–629,
//! 2004). Each secondary node A carries one normal pressure p_A and one weighted gap
//!
//!   g̃_A = Σ_segments ∫ Φ_A (x_m − x_s) · n dA,
//!
//! integrated over the overlaps of its faces with the main faces they face
//! ([`fieldcraft_fem::mortar`]: projection along the secondary face's normal n, clipping,
//! seven-point rule), with Φ_A the secondary face's own shape functions. The pressure field is
//! interpolated like the displacements, so a uniform pressure across non-matching meshes is
//! transmitted exactly: both sides receive the consistent nodal forces ∫ p N dA of it, which is
//! the contact patch test that node-to-surface contact fails (R. L. Taylor, P. Papadopoulos,
//! *On a patch test for contact problems in two dimensions*, 1991; T. W. McDevitt, T. A.
//! Laursen, *A mortar-finite element formulation for frictional contact problems*, IJNME
//! 48:1525–1547, 2000). A node's gap is normalised by its reference area A_A = ∫ Φ_A dA over its
//! own faces: ĝ_A = g̃_A / A_A, a length (negative: penetration).
//!
//! **Enforcement: augmented Lagrangian** (J. C. Simo, T. A. Laursen, *An augmented Lagrangian
//! treatment of contact problems involving friction*, Comput. Struct. 42(1):97–116, 1992;
//! P. Alart, A. Curnier, *A mixed formulation for frictional contact problems prone to Newton
//! like solution methods*, CMAME 92:353–375, 1991). With the node's multiplier λ_A held fixed,
//! the pressure is p_A = ⟨λ_A − ε_A ĝ_A⟩ (⟨x⟩ = max(x, 0)), the derivative of the potential
//! Π = Σ_A A_A/(2 ε_A) ⟨λ_A − ε_A ĝ_A⟩², whose forces −p_A ∂g̃_A/∂x enter the internal forces and
//! whose Hessian
//!
//!   K_c = Σ_{A closed} (ε_A / A_A) ∂g̃_A ⊗ ∂g̃_A − p_A ∂²g̃_A
//!
//! the tangent: the consistent linearisation of the segmentation, the normal and the projection
//! (the second term is the geometric stiffness of the contact, which matters on curved and
//! sliding surfaces), both from the automatic differentiation of the mortar integrals
//! ([`fieldcraft_fem::dual`]). The penalty only regularises: once Newton has converged an
//! increment, the multipliers take the pressures (Uzawa's update λ_A ← p_A) and the increment is
//! solved again from its solution, until every closed node's gap is within [`GAP_TOLERANCE`] of
//! the thickness of the element under it and no open node keeps a multiplier
//! ([`Interface::augmentation_error`]; at most [`MAX_AUGMENTATIONS`] times, then with a warning).
//! The penalty per length is ε_A = s · 10 · 3Ḡ / h
//! over the secondary faces at the node (Ḡ the mean shear modulus of the element under the face,
//! h its thickness V/A, s the contact's `penalty_scale`): about ten times the stiffness of the
//! element under the node, stiff enough that a few augmentations reach the tolerance and soft
//! enough to leave the tangent well conditioned (3Ḡ rather than the constrained modulus keeps it
//! so for nearly incompressible rubber).
//!
//! **Active set.** The pressure's max() is the active set: each evaluation decides which nodes
//! are closed (p_A > 0) from the current gaps and multipliers, and Newton runs on the
//! semismooth system (Alart & Curnier 1991; A. Popp, M. W. Gee, W. A. Wall, *A finite
//! deformation mortar contact formulation using a primal–dual active set strategy*, IJNME
//! 79:1354–1391, 2009). The nonlinear solver cuts the increment back when the set alternates
//! between two states ([`CHATTER_LIMIT`] times: chattering).
//!
//! **Small and finite sliding.** Without nonlinear geometry the face pairs (the main faces within
//! a main element's thickness of each secondary face) and their mortar integrals
//! D_ai = ∫ Φ_a N_i and M_aj = ∫ Φ_a N_j are found once in the reference configuration, so that
//! g̃_a = n · (Σ_j M_aj x_j − Σ_i D_ai x_i) is linear in the positions (small sliding: exact, and
//! consistent with the small-strain elements); with nonlinear geometry they are searched and
//! integrated on the current configuration at every evaluation (finite sliding).
//!
//! **Search.** Every evaluation of finite sliding: the main faces' boxes, grown by their depth limit, go into a
//! bucket grid ([`fieldcraft_fem::contact::Grid`]); each secondary face takes the main faces in
//! the cells its box overlaps that face it (unit normals more than [`FACING`] from parallel) and
//! share no node with its element. A pair whose overlap reaches deeper than [`DEPTH_FRACTION`]
//! of the main element's thickness has passed through and is left out. The pairs are integrated
//! in parallel and summed in face order, so results do not depend on the thread count.
//!
//! **Initial overclosure.** Before the step, nodes whose gap starts negative are reported:
//! `adjust` (the default) takes the initial overclosure out of their gap (an offset of the
//! weighted gap: strain-free, the nodes are not moved), `refuse` refuses the step naming the
//! first node and its depth.
//!
//! **Output.** Per secondary node the pressure as the node's normal contact force over the area
//! it carries (∫ Φ_A over the overlaps: smoother than the multipliers, which are the coefficients
//! of the pressure field and ring at the contact's edge), the multiplier, the gap and the status;
//! per contact the force on the secondary surface and its normal part.
//!
//! **Single-surface (general) contact** is the same mortar contact with one face list as both
//! sides: a set's faces, or every exterior face of the solids. Each face is a secondary face
//! against every other one that faces it, that shares no node with its element (its neighbours)
//! and that it does not reach through more than half an element (the far side of a thin wall),
//! so two faces pressing together are paired both ways, symmetrically: self-contact of folds.
//!
//! **2D models.** In plane strain, plane stress and axisymmetric models the surfaces are element
//! edges (sets' edges; general contact defaults to every exterior edge of the 2D elements),
//! integrated by [`fieldcraft_fem::mortar::integrate_edges`] with the out-of-plane measure
//! (thickness, or 2π r) in the weights, so A_A, g̃_A and the forces are those of the whole ring
//! or slab. Each edge is oriented so that its normal points out of its element; the element's
//! thickness under it is its area over the edge's length. Points keep z = 0 and the 2D nodes have
//! no Z translation, so the tangent blocks keep their in-plane rows.
//!
//! **Scope.** Linear solid faces (tri3 and quad4 faces: tet4, wedge6 and hex8 elements) in 3D
//! models and linear edges (tri3 and quad4 elements) in 2D ones, with Coulomb friction
//! ([`Interface::friction`]). Quadratic faces and edges and arc-length steps are refused by name.

use std::collections::BTreeMap;

use fieldcraft_fem::contact::{Aabb, Grid, bounds};
use fieldcraft_fem::dual::{Dual, Hyper, Real, VARS};
use fieldcraft_fem::frame::{Vec3, dot, norm, unit};
use fieldcraft_fem::mortar::{Sample, area_vector, edge_area_vector, integrate, integrate_edges};
use fieldcraft_model::{Contact, ContactSurfaces, ElementEdge, ElementFace, ElementKind, Extent, Friction, InitialPenetration, Model, Resolved};
use fieldcraft_parallel as parallel;
use serde::Serialize;

use crate::system::{Mesh, translations};
use crate::{Error, Result};

/// Facing faces: the dot product of their unit normals is below minus this.
pub const FACING: f64 = 0.1;

/// A pair deeper than this fraction of the main element's thickness has passed through.
pub const DEPTH_FRACTION: f64 = 0.5;

/// ε = s · PENALTY_FACTOR · 3Ḡ / h (see the module docs).
pub const PENALTY_FACTOR: f64 = 10.0;

/// Augmentations stop when every closed node's gap is within this fraction of the thickness of
/// the element under it (and no open node keeps a multiplier).
pub const GAP_TOLERANCE: f64 = 1e-5;

/// Once an increment's active set has chattered, a node open by less than this fraction of its
/// element's thickness keeps its penalty in the tangent for the rest of the step (its
/// pressure, and so the residual, is still zero): a node that touches with no pressure, the edge
/// of a contact that lifts, otherwise flips between open and closed from one Newton iteration to
/// the next, each state's correction pushing it into the other, and a smaller increment does not
/// help under proportional loading.
pub const ENGAGED_BAND: f64 = 1e-3; // unit-factors: relative to the element's thickness

/// Most augmentations of one increment.
pub const MAX_AUGMENTATIONS: usize = 25;

/// A node behind a face by less than this fraction of the model's extent is touching it:
/// coordinates are exact to a few units in the last place.
const ROUND_OFF: f64 = 1e-9; // unit-factors: relative to the model's extent

/// An increment is cut back when the active set alternates between two states this many times
/// (and the step takes [`ENGAGED_BAND`] from then on).
pub const CHATTER_LIMIT: usize = 3;

/// Most faces over the contacts of a model.
pub const MAX_CONTACT_FACES: usize = 2_000_000;

/// Secondary faces per parallel task.
const GRAIN: usize = 16;

/// A face of a contact surface: its corners (mesh slots) and their translations.
#[derive(Debug, Clone)]
pub(crate) struct Face {
    pub(crate) nodes: Vec<usize>,
    /// X, Y and Z translations per corner (a 2D model's nodes have no Z).
    pub(crate) dofs: Vec<[Option<usize>; 3]>,
    /// The nodes of its element (never in contact with the face).
    owner: Vec<usize>,
    /// The penetration beyond which a point has passed through, and its element's thickness.
    depth: f64,
    pub(crate) thickness: f64,
    /// The penalty per length its element sets (before the penalty scale).
    stiffness: f64,
}

/// A contact pair as the solver runs it.
#[derive(Debug, Clone)]
pub struct Interface {
    pub name: String,
    secondary: Vec<Face>,
    main: Vec<Face>,
    /// The secondary nodes (mesh slots), sorted, and per node: its reference area A_A, its
    /// penalty ε_A, the thinnest element under it and the offset taken out of its weighted gap
    /// (initial overclosure).
    pub nodes: Vec<usize>,
    area: Vec<f64>,
    penalty: Vec<f64>,
    thickness: Vec<f64>,
    offset: Vec<f64>,
    /// Small sliding (steps without nonlinear geometry): the face pairs and their integrals,
    /// fixed in the reference configuration; `None`: finite sliding, searched and integrated
    /// on the current configuration at every evaluation.
    fixed: Option<Vec<Vec<Fixed>>>,
    /// Coulomb friction (frictionless: none).
    friction: Option<Friction>,
    /// A 2D model's out-of-plane measure (its surfaces are edges); `None`: solid faces.
    extent: Option<Extent>,
}

/// A contact's state carried from one converged increment to the next: per secondary node the
/// augmented Lagrange multiplier, and with friction the tangential traction and whether the node
/// was slipping.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Committed {
    pub multipliers: Vec<f64>,
    pub traction: Vec<Vec3>,
    pub slipping: Vec<bool>,
}

impl Committed {
    /// Take a converged state's friction (the multipliers are the augmentation's).
    pub fn commit(&mut self, state: &State) {
        self.traction.clone_from(&state.traction);
        self.slipping.clone_from(&state.slipping);
    }
}

/// Whether a secondary node presses on the main surface (and with friction, sticks or slips).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Status {
    Open,
    Closed,
    Stick,
    Slip,
}

impl Status {
    pub fn is_closed(self) -> bool {
        self != Status::Open
    }
}

impl ContactResult {
    /// The closed nodes.
    pub fn closed(&self) -> usize {
        self.status.iter().filter(|s| s.is_closed()).count()
    }
}

/// One evaluation of a contact: per secondary node (in [`Interface::nodes`] order).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct State {
    /// ĝ_A, `None` where no main face lies over the node.
    pub gap: Vec<Option<f64>>,
    pub pressure: Vec<f64>,
    /// ∫ Φ_A over the overlaps (the area the pressure acts on).
    pub support: Vec<f64>,
    /// The force the contact puts on each secondary node, and the node's unit outward normal
    /// (its faces' area vectors summed).
    pub nodal_force: Vec<Vec3>,
    pub normal: Vec<Vec3>,
    /// Reporting only: the unit normal of the overlapping face pairs, weighted by their
    /// nodal mortar supports. Uncontacted faces do not contribute. The friction law keeps
    /// its existing `normal`; changing a report must not change the assembled problem.
    pub report_normal: Vec<Vec3>,
    /// Friction: the tangential traction on each secondary node (on it, against its slip), its
    /// slip over the increment (relative to the main surface, normalised like the gap) and
    /// whether it slips; empty without friction.
    pub traction: Vec<Vec3>,
    pub slip: Vec<Vec3>,
    pub slipping: Vec<bool>,
}

impl State {
    /// A hash of which nodes are closed (for the chattering check).
    pub fn active_set(&self) -> u64 {
        // FNV-1a over the closed nodes' indices: deterministic on every platform.
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        for (i, _) in self.pressure.iter().enumerate().filter(|(_, p)| **p > 0.0) {
            for b in (i as u64).to_le_bytes() {
                h ^= u64::from(b);
                h = h.wrapping_mul(0x0100_0000_01b3);
            }
        }
        h
    }
}

/// A contact's result at the end of a step, per secondary node.
#[derive(Debug, Clone, PartialEq, Serialize, serde::Deserialize)]
pub struct ContactResult {
    pub name: String,
    /// The secondary surface's node labels, sorted.
    pub nodes: Vec<u32>,
    /// The contact pressure: the normal contact force on the node over the area it carries
    /// (∫ Φ_A over the overlaps), so that Σ p_A ∫Φ_A is the force across the contact.
    pub pressure: Vec<f64>,
    /// The mortar multipliers (the coefficients of the pressure field, Φ_A-interpolated).
    pub multiplier: Vec<f64>,
    /// Friction: the shear traction (the tangential contact force over the area, on the node)
    /// and the slip over the last increment, per node; empty without friction.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub shear: Vec<Vec3>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub slip: Vec<Vec3>,
    /// The normalised gap ĝ (negative: penetration), `None` where no main face is over the node.
    pub gap: Vec<Option<f64>>,
    pub status: Vec<Status>,
    /// The force the contact puts on the secondary surface, and its normal part Σ p_A ∫Φ_A.
    pub force: Vec3,
    pub normal_force: f64,
    /// The deepest penetration of a closed node.
    pub penetration: f64,
    /// Augmentations over the step.
    pub augmentations: usize,
}

/// The contribution of one face pair to the forces and the tangent.
struct PairTerms {
    /// Mesh slots of the pair's corners (secondary first) and their degrees of freedom.
    nodes: Vec<usize>,
    dofs: Vec<[Option<usize>; 3]>,
    /// Secondary node indices (in [`Interface::nodes`]) of the secondary corners.
    corners: Vec<usize>,
    /// ∂g̃_a/∂x per secondary corner (24 variables; corner c of the pair at 3c, main corners from 12).
    gradient: Vec<[f64; VARS]>,
    /// Σ_a p_a ∂²g̃_a, packed (tangent evaluations only).
    hessian: Option<Vec<f64>>,
}

/// The tangent blocks of a contact: per block its mesh slots and the dense matrix over their
/// translations (row-major, node by node).
pub type Blocks = Vec<(Vec<usize>, Vec<f64>)>;

fn position(mesh: &Mesh, u: &[f64], node: usize, dofs: &[Option<usize>; 3]) -> Vec3 {
    let x = mesh.xyz.get(node).copied().unwrap_or_default();
    [0, 1, 2].map(|k| x[k] + dofs[k].and_then(|d| u.get(d)).copied().unwrap_or(0.0))
}

/// The element's volume (hex8, tet4), or a 2D element's area (tri3, quad4), for the thickness
/// under a face or edge.
fn volume(kind: ElementKind, x: &[Vec3]) -> Option<f64> {
    match kind {
        ElementKind::Hex8 => <[Vec3; 8]>::try_from(x).ok().and_then(|a| fieldcraft_fem::uniform::Hex8::hex8(&a).ok()).map(|k| k.volume),
        ElementKind::Tet4 => <[Vec3; 4]>::try_from(x).ok().and_then(|a| fieldcraft_fem::uniform::Tet4::tet4(&a).ok()).map(|k| k.volume),
        ElementKind::Tri3 | ElementKind::Quad4 => {
            // The shoelace formula over the corners in the XY plane.
            let doubled: f64 = x.iter().zip(x.iter().cycle().skip(1)).map(|(a, b)| a[0] * b[1] - a[1] * b[0]).sum();
            Some(0.5 * doubled.abs())
        }
        _ => None,
    }
}

/// The depth of element `kind` (corners `element`) under its face or edge `x` of measure `area`:
/// the element's volume over the area, or the face's size where the volume is unknown.
pub(crate) fn thickness(kind: ElementKind, element: &[Vec3], x: &[Vec3], area: f64) -> f64 {
    volume(kind, element).filter(|v| *v > 0.0 && area > 0.0).map_or(size(x), |v| v / area)
}

/// The overlap integration of a face pair (`extent`: `None`) or an edge pair of a 2D model.
fn overlap<T: Real>(extent: Option<Extent>, xs: &[[T; 3]], xm: &[[T; 3]], f: impl FnMut(&Sample<T>)) -> Option<()> {
    match extent {
        None => integrate(xs, xm, f),
        Some(e) => integrate_edges(xs, xm, e, f),
    }
}

/// A face's area vector, or an edge's ([`edge_area_vector`]) in a 2D model.
fn area_of(extent: Option<Extent>, x: &[Vec3]) -> Option<Vec3> {
    match extent {
        None => area_vector(x),
        Some(e) => edge_area_vector(x, e),
    }
}

fn size(x: &[Vec3]) -> f64 {
    let b = bounds(x, 0.0);
    (0..3).map(|k| b[1][k] - b[0][k]).fold(0.0, f64::max)
}

/// The faces of set `set` for definition `owner` (``contact `name` ``, ``tie `name` ``), refused
/// unless they are linear faces of analysed solid elements.
pub(crate) fn faces(model: &Model, resolved: &Resolved, mesh: &Mesh, owner: &str, set: &str) -> Result<Vec<Face>> {
    face_list(model, mesh, owner, resolved.sets.get(set).map(|s| s.faces.as_slice()).unwrap_or_default())
}

/// The edges of set `set` of a 2D model for definition `owner`, refused unless they are linear
/// edges of analysed 2D elements.
fn edges(model: &Model, resolved: &Resolved, mesh: &Mesh, owner: &str, set: &str) -> Result<Vec<Face>> {
    edge_list(model, mesh, owner, resolved.sets.get(set).map(|s| s.edges.as_slice()).unwrap_or_default())
}

/// The faces `list` for definition `owner`, as [`faces`].
fn face_list(model: &Model, mesh: &Mesh, owner: &str, list: &[ElementFace]) -> Result<Vec<Face>> {
    list.iter()
        .map(|f| {
            let at = |what: String| Error::Invalid(format!("{owner}: element {} face {}: {what}", f.element, f.face));
            let labels = model.face_nodes(*f).ok_or_else(|| at("no such face".into()))?;
            surface_face(mesh, owner, f.element, labels, false, at)
        })
        .collect()
}

/// The edges `list` of a 2D model for definition `owner`, as [`edges`].
fn edge_list(model: &Model, mesh: &Mesh, owner: &str, list: &[ElementEdge]) -> Result<Vec<Face>> {
    list.iter()
        .map(|e| {
            let at = |what: String| Error::Invalid(format!("{owner}: element {} edge {}: {what}", e.element, e.edge));
            let labels = model.edge_nodes(*e).ok_or_else(|| at("no such edge".into()))?;
            surface_face(mesh, owner, e.element, labels, true, at)
        })
        .collect()
}

/// One face (or, `edge`, one edge of a 2D element) of element `element` with nodes `labels`.
fn surface_face(mesh: &Mesh, owner: &str, element: u32, labels: Vec<u32>, edge: bool, at: impl Fn(String) -> Error) -> Result<Face> {
    let el = mesh.element(element).ok_or_else(|| at("the element is not a solid or 2D element of the analysis".into()))?;
    let (linear, takes) = if edge {
        (labels.len() == 2 && el.kind.dimension() == 2, "the linear edges of tri3 and quad4 elements")
    } else {
        (matches!(labels.len(), 3 | 4) && el.kind.dimension() == 3, "the linear faces of tet4, wedge6 and hex8 elements")
    };
    if !linear {
        let what = if edge { "edges" } else { "faces" };
        return Err(Error::NotSupported(format!(
            "{owner} on the {what} of {} elements (element {element}) in implicit analyses; it takes {takes}",
            el.kind.name()
        )));
    }
    let mut nodes: Vec<usize> = labels.iter().map(|l| mesh.node(*l).ok_or_else(|| at(format!("node {l} is not analysed")))).collect::<Result<_>>()?;
    let mut x: Vec<Vec3> = nodes.iter().map(|n| mesh.xyz.get(*n).copied().unwrap_or_default()).collect();
    if edge {
        // Outward: the edge's normal (t_y, −t_x) away from its element's centroid.
        let count = el.xyz.len().max(1) as f64;
        let centroid: Vec3 = std::array::from_fn(|k| el.xyz.iter().map(|p| p[k]).sum::<f64>() / count);
        let (Some(a), Some(b)) = (x.first(), x.get(1)) else { return Err(at("the edge has no size".into())) };
        let mid = [0.5 * (a[0] + b[0]) - centroid[0], 0.5 * (a[1] + b[1]) - centroid[1]];
        if (b[1] - a[1]) * mid[0] - (b[0] - a[0]) * mid[1] < 0.0 {
            nodes.reverse();
            x.reverse();
        }
    }
    let dofs = nodes
        .iter()
        .map(|n| {
            let d = translations(mesh, *n)?;
            match d.as_slice() {
                [a, b, c] => Ok([Some(*a), Some(*b), Some(*c)]),
                [a, b] if edge => Ok([Some(*a), Some(*b), None]),
                _ => Err(Error::Internal("node translations".into())),
            }
        })
        .collect::<Result<Vec<_>>>()?;
    // The face's area, or the edge's length.
    let area = match x.as_slice() {
        [a, b] if edge => norm([b[0] - a[0], b[1] - a[1], 0.0]),
        _ => area_vector(&x).map_or(0.0, norm),
    };
    let thickness = thickness(el.kind, &el.xyz, &x, area);
    if !(thickness > 0.0 && thickness.is_finite()) {
        return Err(at(format!("the {} has no size", if edge { "edge" } else { "face" })));
    }
    let shear = (el.d[3][3] + el.d[4][4] + el.d[5][5]) / 3.0;
    let mut owner = el.nodes.clone();
    owner.sort_unstable();
    Ok(Face { nodes, dofs, owner, depth: DEPTH_FRACTION * thickness, thickness, stiffness: PENALTY_FACTOR * 3.0 * shear / thickness })
}

/// A face pair of small-sliding contact: its mortar integrals in the reference configuration,
/// D_ai = ∫ Φ_a N_i (secondary corners i) and M_aj = ∫ Φ_a N_j (main corners j), and the
/// secondary face's normal there, so that g̃_a = n · (Σ_j M_aj x_j − Σ_i D_ai x_i).
#[derive(Debug, Clone)]
struct Fixed {
    main: usize,
    d: [[f64; 4]; 4],
    m: [[f64; 4]; 4],
    normal: Vec3,
    support: [f64; 4],
}

impl Fixed {
    /// B_a, the slip operator ([`slip_operator`]) of the reference configuration.
    fn slip_jacobian(&self, a: usize) -> [[f64; VARS]; 3] {
        let mut j = [[0.0; VARS]; 3];
        let (Some(d), Some(m)) = (self.d.get(a), self.m.get(a)) else { return j };
        let n = self.normal;
        for k in 0..3 {
            for l in 0..3 {
                let proj = if k == l { 1.0 } else { 0.0 } - n[k] * n[l];
                for i in 0..4 {
                    j[k][3 * i + l] = d[i] * proj;
                    j[k][12 + 3 * i + l] = -m[i] * proj;
                }
            }
        }
        j
    }

    /// ∂g̃_a/∂x over the pair's variables (secondary corners from 0, main corners from 12).
    fn gradient(&self, a: usize) -> [f64; VARS] {
        let mut g = [0.0; VARS];
        let (Some(d), Some(m)) = (self.d.get(a), self.m.get(a)) else { return g };
        for k in 0..3 {
            for i in 0..4 {
                g[3 * i + k] = -d[i] * self.normal[k];
                g[12 + 3 * i + k] = m[i] * self.normal[k];
            }
        }
        g
    }
}

/// What a secondary face finds against the main surface: per main face, ∫ Φ_a g and ∫ Φ_a for
/// each corner a.
struct Found {
    main: usize,
    gap: [f64; 4],
    support: [f64; 4],
    /// The secondary projection normal actually used by this overlap's integrals.
    normal: Vec3,
}

/// A reporting direction from only the overlaps carrying a node's contact. Opposing supports
/// can cancel; a zero sum has no direction, so it reports no normal component instead of
/// inventing one. This does not participate in the contact or friction law.
fn report_normals(pairs: &[Vec<Found>], corner_index: &[Vec<usize>], support: &[f64]) -> Vec<Vec3> {
    let mut normals = vec![[0.0; 3]; support.len()];
    for (found, corners) in pairs.iter().zip(corner_index) {
        for pair in found {
            for (a, i) in corners.iter().enumerate() {
                let (Some(normal), Some(&total), Some(&area)) = (normals.get_mut(*i), support.get(*i), pair.support.get(a)) else { continue };
                if total > 0.0 {
                    // Scale before adding so large finite overlap areas cannot overflow
                    // the normal's length. This is a direction, not another area measure.
                    let weight = area / total;
                    for (component, direction) in normal.iter_mut().zip(pair.normal) {
                        *component += weight * direction;
                    }
                }
            }
        }
    }
    for normal in &mut normals {
        *normal = unit(*normal).unwrap_or_default();
    }
    normals
}

/// The weighted gaps g̃_a of one face pair (corners `xs`, `xm`) over any scalar.
fn weighted_gaps<T: Real>(extent: Option<Extent>, xs: &[[T; 3]], xm: &[[T; 3]]) -> Option<[T; 4]> {
    let mut g = [T::cst(0.0); 4];
    overlap(extent, xs, xm, |p| {
        let w = p.weight * p.gap();
        for a in 0..4 {
            g[a] = g[a] + w * p.secondary[a];
        }
    })?;
    Some(g)
}

/// The slip operator of one face pair (corners `xs`, `xm`): per secondary corner a the matrix
/// B_a with w_a = B_a Δx, the weighted tangential slip ∫ Φ_a (I − n n)(Δx_s − Δx_m) dA of the
/// pair's corners' motions Δx (rows: its components; columns: the pair's variables, secondary
/// corners from 0, main corners from 12), with the integrals on the current configuration: the
/// rate form of the mortar slip, B_a = [∫ Φ_a N_i P, −∫ Φ_a N_j P] (McDevitt & Laursen 2000; Puso &
/// Laursen 2004), P = I − n n.
fn slip_operator(extent: Option<Extent>, xs: &[Vec3], xm: &[Vec3]) -> Option<[[[f64; VARS]; 3]; 4]> {
    let mut b = [[[0.0; VARS]; 3]; 4];
    overlap(extent, xs, xm, |p| {
        let n = p.normal;
        for a in 0..4 {
            let wa = p.weight * p.secondary[a];
            if wa == 0.0 {
                continue;
            }
            for k in 0..3 {
                for l in 0..3 {
                    let proj = wa * (if k == l { 1.0 } else { 0.0 } - n[k] * n[l]);
                    for c in 0..4 {
                        b[a][k][3 * c + l] += proj * p.secondary[c];
                        b[a][k][12 + 3 * c + l] -= proj * p.main[c];
                    }
                }
            }
        }
    })?;
    Some(b)
}

/// The variables of a face pair: secondary corners at 0, main corners from 12.
fn variables<T: Real>(xs: &[Vec3], xm: &[Vec3], var: impl Fn(f64, usize) -> T) -> (Vec<[T; 3]>, Vec<[T; 3]>) {
    let s = xs.iter().enumerate().map(|(i, x)| [0, 1, 2].map(|k| var(x[k], 3 * i + k))).collect();
    let m = xm.iter().enumerate().map(|(j, x)| [0, 1, 2].map(|k| var(x[k], 12 + 3 * j + k))).collect();
    (s, m)
}

impl Interface {
    /// The contact pairs of a model, ready for the nonlinear solver (initial overclosures
    /// adjusted or refused, and reported in `warnings`).
    pub fn all(model: &Model, resolved: &Resolved, mesh: &Mesh, finite: bool, warnings: &mut Vec<String>) -> Result<Vec<Interface>> {
        let mut out = Vec::new();
        let mut total = 0usize;
        for (name, contact) in model.contacts() {
            let mut interface = Interface::new(model, resolved, mesh, name, contact)?;
            if !finite {
                interface.fixed = Some(interface.fixed_pairs(mesh)?);
            }
            total = total.saturating_add(interface.secondary.len() + interface.main.len());
            if total > MAX_CONTACT_FACES {
                return Err(Error::Limit(format!("the contacts have more than {MAX_CONTACT_FACES} faces")));
            }
            out.push(interface.overclosure(mesh, contact.initial_penetration, warnings)?);
        }
        Ok(out)
    }

    fn new(model: &Model, resolved: &Resolved, mesh: &Mesh, name: &str, contact: &Contact) -> Result<Interface> {
        let owner = format!("contact `{name}`");
        let extent = mesh.theory.map(|t| t.extent());
        // A set's faces or (2D) edges, or every exterior one of the model.
        let surface = |set: Option<&String>| match (set, extent) {
            (Some(s), None) => faces(model, resolved, mesh, &owner, s),
            (Some(s), Some(_)) => edges(model, resolved, mesh, &owner, s),
            (None, None) => face_list(model, mesh, &owner, &model.exterior_faces()),
            (None, Some(_)) => edge_list(model, mesh, &owner, &model.exterior_edges()),
        };
        let (secondary, main) = match &contact.surfaces {
            ContactSurfaces::Pair { main, secondary } => (surface(Some(secondary))?, surface(Some(main))?),
            // Single surface: every face against every other one that faces it and shares no node
            // with its element (its neighbours), both ways.
            ContactSurfaces::General { surface: set } => {
                let list = surface(set.as_ref())?;
                (list.clone(), list)
            }
        };
        let mut nodes: Vec<usize> = secondary.iter().flat_map(|f| f.nodes.iter().copied()).collect();
        nodes.sort_unstable();
        nodes.dedup();
        let (mut area, mut penalty, mut thickness) = (vec![0.0; nodes.len()], vec![0.0f64; nodes.len()], vec![f64::INFINITY; nodes.len()]);
        for f in &secondary {
            let x: Vec<Vec3> = f.nodes.iter().map(|n| mesh.xyz.get(*n).copied().unwrap_or_default()).collect();
            // ∫ Φ_a over the face itself (the face against its own mirror image).
            let mirror: Vec<Vec3> = x.iter().rev().copied().collect();
            let mut own = [0.0; 4];
            overlap(extent, &x, &mirror, |p| {
                for a in 0..4 {
                    own[a] += p.weight * p.secondary[a];
                }
            });
            for (n, share) in f.nodes.iter().zip(own) {
                if let Ok(i) = nodes.binary_search(n) {
                    if let Some(a) = area.get_mut(i) {
                        *a += share;
                    }
                    if let Some(e) = penalty.get_mut(i) {
                        *e = e.max(contact.penalty_scale * f.stiffness);
                    }
                    if let Some(h) = thickness.get_mut(i) {
                        *h = h.min(f.thickness);
                    }
                }
            }
        }
        if let Some(i) = area.iter().position(|a| !(*a > 0.0 && a.is_finite())) {
            let label = nodes.get(i).and_then(|n| mesh.labels.get(*n)).copied().unwrap_or(0);
            return Err(Error::Invalid(format!("contact `{name}`: the secondary faces at node {label} have no area")));
        }
        let offset = vec![0.0; nodes.len()];
        Ok(Interface {
            name: name.to_string(),
            secondary,
            main,
            nodes,
            area,
            penalty,
            thickness,
            offset,
            fixed: None,
            friction: (!contact.friction.is_frictionless()).then_some(contact.friction),
            extent,
        })
    }

    /// Adjust (offset) or refuse the nodes that start behind the main surface.
    fn overclosure(mut self, mesh: &Mesh, mode: InitialPenetration, warnings: &mut Vec<String>) -> Result<Interface> {
        let zero = vec![0.0; mesh.dofs()];
        let (gaps, _) = self.gaps(mesh, &zero)?;
        // The model's extent (relative, so the test reads the same in any unit of length).
        let scale = mesh.xyz.iter().flatten().fold(0.0f64, |m, c| m.max(c.abs()));
        let (mut count, mut deepest) = (0usize, (0.0f64, 0u32));
        for (i, g) in gaps.iter().enumerate() {
            let (Some(g), Some(a)) = (g, self.area.get(i)) else { continue };
            let depth = -g / a;
            // Coordinates are exact to a few units in the last place: a node a hair behind a
            // face it was meshed on is touching.
            if depth <= ROUND_OFF * scale {
                continue;
            }
            let label = self.nodes.get(i).and_then(|n| mesh.labels.get(*n)).copied().unwrap_or(0);
            if mode == InitialPenetration::Refuse {
                return Err(Error::Invalid(format!(
                    "contact `{}`: node {label} starts {depth:.4e} inside the main surface (move it out, or let the contact adjust it: initial_penetration \"adjust\")",
                    self.name
                )));
            }
            if let Some(o) = self.offset.get_mut(i) {
                *o = *g;
            }
            count += 1;
            if depth > deepest.0 {
                deepest = (depth, label);
            }
        }
        if count > 0 {
            warnings.push(format!(
                "contact `{}`: {count} node(s) start inside the main surface (deepest node {}, {:.4e}); their initial overclosure is taken out of the gap",
                self.name, deepest.1, deepest.0
            ));
        }
        Ok(self)
    }

    fn corners(mesh: &Mesh, u: &[f64], f: &Face) -> Vec<Vec3> {
        f.nodes.iter().zip(&f.dofs).map(|(n, d)| position(mesh, u, *n, d)).collect()
    }

    /// The face pairs: per secondary face, the main faces it faces and their integrals.
    fn pairs(&self, mesh: &Mesh, u: &[f64]) -> Result<Vec<Vec<Found>>> {
        let Some(fixed) = &self.fixed else { return self.search(mesh, u, 1.0, |_, _, found| found) };
        Ok(parallel::map_grained(self.secondary.len(), GRAIN, |s| {
            let (Some(f), Some(list)) = (self.secondary.get(s), fixed.get(s)) else { return Vec::new() };
            let xs = Interface::corners(mesh, u, f);
            list.iter()
                .map(|p| {
                    let xm = self.main.get(p.main).map(|m| Interface::corners(mesh, u, m)).unwrap_or_default();
                    let gap = std::array::from_fn(|a| {
                        let (Some(d), Some(m)) = (p.d.get(a), p.m.get(a)) else { return 0.0 };
                        let mut v = [0.0; 3];
                        for (x, c) in xm.iter().zip(m) {
                            (0..3).for_each(|k| v[k] += c * x[k]);
                        }
                        for (x, c) in xs.iter().zip(d) {
                            (0..3).for_each(|k| v[k] -= c * x[k]);
                        }
                        dot(v, p.normal)
                    });
                    Found { main: p.main, gap, support: p.support, normal: p.normal }
                })
                .collect()
        }))
    }

    /// The small-sliding pairs: found and integrated once, in the reference configuration.
    fn fixed_pairs(&self, mesh: &Mesh) -> Result<Vec<Vec<Fixed>>> {
        let zero = vec![0.0; mesh.dofs()];
        // Faces within a main element's thickness may touch later.
        self.search(mesh, &zero, 2.0, |xs, xm, found| {
            found
                .into_iter()
                .filter_map(|f| {
                    let xm = xm.get(f.main)?;
                    let (mut d, mut m, mut normal) = ([[0.0; 4]; 4], [[0.0; 4]; 4], [0.0; 3]);
                    overlap(self.extent, xs, xm, |p| {
                        normal = p.normal;
                        for a in 0..4 {
                            for i in 0..4 {
                                d[a][i] += p.weight * p.secondary[a] * p.secondary[i];
                                m[a][i] += p.weight * p.secondary[a] * p.main[i];
                            }
                        }
                    })?;
                    Some(Fixed { main: f.main, d, m, normal, support: f.support })
                })
                .collect()
        })
    }

    /// The global and local search at displacements `u`: per secondary face (corners `xs`), the
    /// main faces (corners `xm`, all of them) within `reach` times their depth limit that it
    /// faces, with their integrals, passed through `keep`.
    fn search<T: Send>(&self, mesh: &Mesh, u: &[f64], reach: f64, keep: impl Fn(&[Vec3], &[Vec<Vec3>], Vec<Found>) -> T + Sync) -> Result<Vec<T>> {
        let main_x: Vec<Vec<Vec3>> = self.main.iter().map(|f| Interface::corners(mesh, u, f)).collect();
        let boxes: Vec<Aabb> = self.main.iter().zip(&main_x).map(|(f, x)| bounds(x, reach * f.depth)).collect();
        let grid = Grid::new(&boxes, 4 * boxes.len() + 64);
        let main_normals: Vec<Option<Vec3>> = main_x.iter().map(|x| area_of(self.extent, x).filter(|a| norm(*a) > 0.0)).collect();
        Ok(parallel::map_grained(self.secondary.len(), GRAIN, |i| {
            let mut found = Vec::new();
            let Some(f) = self.secondary.get(i) else { return keep(&[], &main_x, found) };
            let xs = Interface::corners(mesh, u, f);
            let Some(ns) = area_of(self.extent, &xs).filter(|a| norm(*a) > 0.0) else { return keep(&xs, &main_x, found) };
            let mut candidates = Vec::new();
            grid.candidates_in(&bounds(&xs, 0.0), &mut candidates);
            for &c in &candidates {
                let c = c as usize;
                let (Some(m), Some(xm), Some(Some(nm))) = (self.main.get(c), main_x.get(c), main_normals.get(c)) else { continue };
                if m.nodes.iter().any(|n| f.owner.binary_search(n).is_ok()) || f.nodes.iter().any(|n| m.owner.binary_search(n).is_ok()) {
                    continue;
                }
                if dot(ns, *nm) >= -FACING * norm(ns) * norm(*nm) {
                    continue;
                }
                let (mut gap, mut support, mut normal, mut deepest, mut any) = ([0.0; 4], [0.0; 4], [0.0; 3], f64::INFINITY, false);
                overlap(self.extent, &xs, xm, |p| {
                    any = true;
                    // The projection normal is constant over a face-pair overlap.
                    normal = p.normal;
                    let g = p.gap();
                    deepest = deepest.min(g);
                    for a in 0..4 {
                        gap[a] += p.weight * g * p.secondary[a];
                        support[a] += p.weight * p.secondary[a];
                    }
                });
                if any && deepest >= -m.depth && gap.iter().chain(&support).all(|v| v.is_finite()) {
                    found.push(Found { main: c, gap, support, normal });
                }
            }
            keep(&xs, &main_x, found)
        }))
    }

    /// The weighted gaps g̃_A less the offsets (`None` where no main face lies over the node), and
    /// the supports ∫ Φ_A, at displacements `u`.
    fn gaps(&self, mesh: &Mesh, u: &[f64]) -> Result<(Vec<Option<f64>>, Vec<f64>)> {
        let pairs = self.pairs(mesh, u)?;
        let (gaps, support, _) = self.sums(&pairs);
        Ok((gaps, support))
    }

    /// Node sums of the pairs' integrals, in face order.
    fn sums(&self, pairs: &[Vec<Found>]) -> (Vec<Option<f64>>, Vec<f64>, Vec<Vec<usize>>) {
        let n = self.nodes.len();
        let (mut gap, mut support) = (vec![None; n], vec![0.0; n]);
        let mut index = Vec::with_capacity(self.secondary.len());
        for (f, found) in self.secondary.iter().zip(pairs) {
            let corners: Vec<usize> = f.nodes.iter().filter_map(|x| self.nodes.binary_search(x).ok()).collect();
            for p in found {
                for (a, i) in corners.iter().enumerate() {
                    if let (Some(g), Some(s)) = (gap.get_mut(*i), support.get_mut(*i)) {
                        *g = Some(g.unwrap_or(0.0) + p.gap.get(a).copied().unwrap_or(0.0));
                        *s += p.support.get(a).copied().unwrap_or(0.0);
                    }
                }
            }
            index.push(corners);
        }
        for (g, o) in gap.iter_mut().zip(&self.offset) {
            if let Some(g) = g {
                *g -= o;
            }
        }
        (gap, support, index)
    }

    /// The state a step starts from: no multipliers, no friction traction.
    pub fn committed(&self) -> Committed {
        let n = self.nodes.len();
        let friction = self.friction.is_some();
        Committed {
            multipliers: vec![0.0; n],
            traction: if friction { vec![[0.0; 3]; n] } else { Vec::new() },
            slipping: vec![false; if friction { n } else { 0 }],
        }
    }

    /// The contact's forces at displacements `u` (added to `internal`), its tangent blocks when
    /// `tangent`, and its state, from the last converged increment's state `committed` and
    /// displacements `reference`; nodes open by less than `band` of their element's thickness
    /// keep their penalty in the tangent ([`ENGAGED_BAND`]; 0: only closed nodes).
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate(
        &self,
        mesh: &Mesh,
        u: &[f64],
        reference: &[f64],
        committed: &Committed,
        band: f64,
        tangent: bool,
        internal: &mut [f64],
        blocks: &mut Blocks,
    ) -> Result<State> {
        let multipliers = committed.multipliers.as_slice();
        let pairs = self.pairs(mesh, u)?;
        let (weighted, support, corner_index) = self.sums(&pairs);
        let report_normal = report_normals(&pairs, &corner_index, &support);
        let n = self.nodes.len();
        let mut state = State {
            gap: vec![None; n],
            pressure: vec![0.0; n],
            support,
            nodal_force: vec![[0.0; 3]; n],
            normal: vec![[0.0; 3]; n],
            report_normal,
            traction: Vec::new(),
            slip: Vec::new(),
            slipping: Vec::new(),
        };
        for (f, corners) in self.secondary.iter().zip(&corner_index) {
            let Some(a) = area_of(self.extent, &Interface::corners(mesh, u, f)) else { continue };
            for i in corners {
                if let Some(v) = state.normal.get_mut(*i) {
                    (0..3).for_each(|k| v[k] += a[k]);
                }
            }
        }
        for v in &mut state.normal {
            *v = unit(*v).unwrap_or(*v);
        }
        // Nodes closed, or open by less than `band` of their element's thickness: their penalty
        // stays in the tangent.
        let mut engaged = vec![false; n];
        for i in 0..self.nodes.len() {
            let (Some(Some(g)), Some(a), Some(e)) = (weighted.get(i), self.area.get(i), self.penalty.get(i)) else { continue };
            let gap = g / a;
            let lambda = multipliers.get(i).copied().unwrap_or(0.0);
            if let (Some(slot), Some(p), Some(h)) = (state.gap.get_mut(i), state.pressure.get_mut(i), self.thickness.get(i)) {
                *slot = Some(gap);
                *p = (lambda - e * gap).max(0.0);
                engaged[i] = lambda - e * gap > -e * band * h;
            }
        }
        // The pairs with a closed secondary corner, linearised.
        let mut active: Vec<(usize, usize)> = Vec::new();
        for (s, found) in pairs.iter().enumerate() {
            let closed = corner_index.get(s).is_some_and(|c| c.iter().any(|i| engaged.get(*i).copied().unwrap_or(false)));
            if closed {
                active.extend((0..found.len()).map(|k| (s, k)));
            }
        }
        let terms = parallel::map_grained(active.len(), GRAIN, |t| -> Result<Option<PairTerms>> {
            let Some(&(s, k)) = active.get(t) else { return Ok(None) };
            let (Some(f), Some(found), Some(corners)) = (self.secondary.get(s), pairs.get(s).and_then(|p| p.get(k)), corner_index.get(s)) else {
                return Ok(None);
            };
            let Some(m) = self.main.get(found.main) else { return Ok(None) };
            let (xs, xm) = (Interface::corners(mesh, u, f), Interface::corners(mesh, u, m));
            let pressure: Vec<f64> = corners.iter().map(|i| state.pressure.get(*i).copied().unwrap_or(0.0)).collect();
            let fixed = self.fixed.as_ref().and_then(|f| f.get(s)).and_then(|l| l.get(k));
            let (gradient, hessian) = if let Some(fixed) = fixed {
                // Small sliding: linear in the positions.
                ((0..corners.len()).map(|a| fixed.gradient(a)).collect(), None)
            } else if tangent {
                let (vs, vm) = variables(&xs, &xm, Hyper::var);
                let Some(g) = weighted_gaps(self.extent, &vs, &vm) else { return Ok(None) };
                let mut h = vec![0.0; fieldcraft_fem::dual::PACKED];
                for (ga, p) in g.iter().zip(&pressure) {
                    h.iter_mut().zip(&ga.h).for_each(|(x, y)| *x += p * y);
                }
                (g.iter().take(corners.len()).map(|ga| ga.g).collect::<Vec<_>>(), Some(h))
            } else {
                let (vs, vm) = variables(&xs, &xm, Dual::var);
                let Some(g) = weighted_gaps(self.extent, &vs, &vm) else { return Ok(None) };
                (g.iter().take(corners.len()).map(|ga| ga.g).collect(), None)
            };
            let nodes = f.nodes.iter().chain(&m.nodes).copied().collect();
            let dofs = f.dofs.iter().chain(&m.dofs).copied().collect();
            Ok(Some(PairTerms { nodes, dofs, corners: corners.clone(), gradient, hessian }))
        });
        // Forces −p_A ∂g̃_A/∂x, and per closed node its gradient over the nodes it reaches.
        let mut node_gradients: BTreeMap<usize, BTreeMap<usize, Vec3>> = BTreeMap::new();
        for t in terms {
            let Some(t) = t? else { continue };
            for (a, i) in t.corners.iter().enumerate() {
                let p = state.pressure.get(*i).copied().unwrap_or(0.0);
                let Some(g) = t.gradient.get(a) else { continue };
                for (c, (node, dofs)) in t.nodes.iter().zip(&t.dofs).enumerate() {
                    let var = if c < t.corners.len() { 3 * c } else { 12 + 3 * (c - t.corners.len()) };
                    let gv = [0, 1, 2].map(|k| g.get(var + k).copied().unwrap_or(0.0));
                    let own = t.corners.get(c).and_then(|j| state.nodal_force.get_mut(*j));
                    let mut on = [0.0; 3];
                    for k in 0..3 {
                        let f = -p * gv[k];
                        if let Some(x) = dofs[k].and_then(|d| internal.get_mut(d)) {
                            *x += f;
                        }
                        on[k] = -f;
                    }
                    if let Some(own) = own {
                        (0..3).for_each(|k| own[k] += on[k]);
                    }
                    if tangent && engaged.get(*i).copied().unwrap_or(false) {
                        let entry = node_gradients.entry(*i).or_default().entry(*node).or_insert([0.0; 3]);
                        for k in 0..3 {
                            entry[k] += gv[k];
                        }
                    }
                }
            }
            if let Some(h) = &t.hessian {
                // −Σ_a p_a ∂²g̃_a over the pair's nodes.
                let count = t.nodes.len();
                let s = t.corners.len();
                let var = |c: usize, k: usize| if c < s { 3 * c + k } else { 12 + 3 * (c - s) + k };
                let width = 3 * count;
                let mut data = vec![0.0; width * width];
                for ci in 0..count {
                    for ki in 0..3 {
                        for cj in 0..count {
                            for kj in 0..3 {
                                let at = fieldcraft_fem::dual::packed(var(ci, ki), var(cj, kj));
                                if let Some(x) = data.get_mut((3 * ci + ki) * width + 3 * cj + kj) {
                                    *x = -h.get(at).copied().unwrap_or(0.0);
                                }
                            }
                        }
                    }
                }
                blocks.push((t.nodes.clone(), data));
            }
        }
        // (ε_A / A_A) ∂g̃_A ⊗ ∂g̃_A per closed node.
        for (i, gradient) in node_gradients {
            let (Some(e), Some(a)) = (self.penalty.get(i), self.area.get(i)) else { continue };
            let c = e / a;
            let nodes: Vec<usize> = gradient.keys().copied().collect();
            let g: Vec<f64> = gradient.values().flatten().copied().collect();
            let width = g.len();
            let mut data = vec![0.0; width * width];
            for (r, gr) in g.iter().enumerate() {
                for (s, gs) in g.iter().enumerate() {
                    if let Some(x) = data.get_mut(r * width + s) {
                        *x = c * gr * gs;
                    }
                }
            }
            blocks.push((nodes, data));
        }
        if let Some(law) = self.friction {
            self.friction(mesh, u, reference, committed, (&pairs, &corner_index, &active), law, tangent, &mut state, internal, blocks)?;
        }
        Ok(state)
    }

    /// Coulomb friction at the closed secondary nodes by a return map (Wriggers 2006, §5.2;
    /// Simo & Laursen 1992): the weighted tangential slip of each node since the last converged
    /// increment, w_A = Σ ∫ Φ_A (I − n n)(Δx_s − Δx_m) dA = B_A Δx ([`slip_operator`]) and
    /// ĝ_T = w_A / A_A; the trial traction
    /// t* = P t_n − ε_A ĝ_T (P projects the committed traction on the current tangent plane);
    /// stick, t = t*, when |t*| ≤ μ p_A, else slip back to the cone, t = μ p_A t*/|t*|, μ the
    /// static coefficient for a node that stuck in the last increment and the kinetic one for
    /// one that slipped. The traction does the virtual work t · B_A δx, so the forces −B_Aᵀ t enter
    /// the internal forces (B_A at the current configuration: the friction force is exactly μ
    /// times the normal one where nodes slip). The tangent keeps the symmetric stick/slip
    /// stiffness (ε_A/A_A) Bᵀ B, or (ε_A μ p_A/(A_A |t*|)) Bᵀ (I − t̂ t̂) B for slip, and leaves out
    /// the coupling to the pressure (μ t̂ ⊗ ∂p_A/∂x, not symmetric) and the change of B with the
    /// geometry: the factorization stays symmetric and the residual exact, so only the rate of
    /// convergence can suffer while nodes slip.
    #[allow(clippy::too_many_arguments)]
    fn friction(
        &self,
        mesh: &Mesh,
        u: &[f64],
        reference: &[f64],
        committed: &Committed,
        (pairs, corner_index, active): (&[Vec<Found>], &[Vec<usize>], &[(usize, usize)]),
        law: Friction,
        tangent: bool,
        state: &mut State,
        internal: &mut [f64],
        blocks: &mut Blocks,
    ) -> Result<()> {
        let n = self.nodes.len();
        state.traction = vec![[0.0; 3]; n];
        state.slip = vec![[0.0; 3]; n];
        state.slipping = vec![false; n];
        // Per active pair: its corners' weighted slips and their Jacobians.
        let terms = parallel::map_grained(active.len(), GRAIN, |t| {
            let &(s, k) = active.get(t)?;
            let (f, corners) = (self.secondary.get(s)?, corner_index.get(s)?);
            let m = self.main.get(pairs.get(s)?.get(k)?.main)?;
            let (xs, xm) = (Interface::corners(mesh, u, f), Interface::corners(mesh, u, m));
            let (rs, rm) = (Interface::corners(mesh, reference, f), Interface::corners(mesh, reference, m));
            let operator: Vec<[[f64; VARS]; 3]> = match self.fixed.as_ref().and_then(|f| f.get(s)).and_then(|l| l.get(k)) {
                Some(fixed) => (0..corners.len()).map(|a| fixed.slip_jacobian(a)).collect(),
                None => slip_operator(self.extent, &xs, &xm)?.into_iter().take(corners.len()).collect(),
            };
            let mut delta = [0.0; VARS];
            for (i, (x, r)) in xs.iter().zip(&rs).enumerate() {
                (0..3).for_each(|c| delta[3 * i + c] = x[c] - r[c]);
            }
            for (j, (x, r)) in xm.iter().zip(&rm).enumerate() {
                (0..3).for_each(|c| delta[12 + 3 * j + c] = x[c] - r[c]);
            }
            let slips: Vec<(Vec3, [[f64; VARS]; 3])> =
                operator.into_iter().map(|b| (b.map(|row| row.iter().zip(&delta).map(|(x, y)| x * y).sum()), b)).collect();
            let nodes: Vec<usize> = f.nodes.iter().chain(&m.nodes).copied().collect();
            let dofs: Vec<[Option<usize>; 3]> = f.dofs.iter().chain(&m.dofs).copied().collect();
            Some((nodes, dofs, corners.clone(), slips))
        });
        let terms: Vec<_> = terms.into_iter().flatten().collect();
        let var = |count: usize, c: usize, l: usize| if c < count { 3 * c + l } else { 12 + 3 * (c - count) + l };
        // Node sums: the slip, and with the tangent J_A over the nodes the node reaches.
        let mut total = vec![[0.0; 3]; n];
        let mut jacobians: BTreeMap<usize, BTreeMap<usize, [[f64; 3]; 3]>> = BTreeMap::new();
        for (nodes, _, corners, slips) in &terms {
            for (a, i) in corners.iter().enumerate() {
                let (Some((w, jac)), Some(sum)) = (slips.get(a), total.get_mut(*i)) else { continue };
                (0..3).for_each(|k| sum[k] += w[k]);
                if tangent && state.pressure.get(*i).is_some_and(|p| *p > 0.0) {
                    for (c, node) in nodes.iter().enumerate() {
                        let block = jacobians.entry(*i).or_default().entry(*node).or_insert([[0.0; 3]; 3]);
                        for k in 0..3 {
                            for l in 0..3 {
                                block[k][l] += jac[k].get(var(corners.len(), c, l)).copied().unwrap_or(0.0);
                            }
                        }
                    }
                }
            }
        }
        // The return map, and per node the tangent's middle matrix.
        let mut middle: Vec<Option<[[f64; 3]; 3]>> = vec![None; n];
        for i in 0..n {
            let (Some(&p), Some(&area), Some(&eps), Some(&normal)) =
                (state.pressure.get(i), self.area.get(i), self.penalty.get(i), state.normal.get(i))
            else {
                continue;
            };
            if !(p > 0.0) {
                continue;
            }
            let old = committed.traction.get(i).copied().unwrap_or_default();
            let along = dot(old, normal);
            let slip = total.get(i).copied().unwrap_or_default().map(|w| w / area);
            let trial: Vec3 = std::array::from_fn(|k| old[k] - along * normal[k] - eps * slip[k]);
            let mu = if committed.slipping.get(i).copied().unwrap_or(false) { law.kinetic_coefficient() } else { law.static_coefficient };
            let (limit, size) = (mu * p, norm(trial));
            let (traction, slipping, m) = if size <= limit {
                (trial, false, std::array::from_fn(|k| std::array::from_fn(|l| if k == l { eps / area } else { 0.0 })))
            } else {
                let dir = trial.map(|c| c / size);
                let c = eps * limit / (size * area);
                (
                    dir.map(|d| limit * d),
                    true,
                    std::array::from_fn(|k| std::array::from_fn(|l| c * (if k == l { 1.0 } else { 0.0 } - dir[k] * dir[l]))),
                )
            };
            if let (Some(t), Some(g), Some(s), Some(mm)) =
                (state.traction.get_mut(i), state.slip.get_mut(i), state.slipping.get_mut(i), middle.get_mut(i))
            {
                (*t, *g, *s, *mm) = (traction, slip, slipping, Some(m));
            }
        }
        // Forces −t · ∂w/∂x.
        for (_, dofs, corners, slips) in &terms {
            for (a, i) in corners.iter().enumerate() {
                let (Some(t), Some((_, jac))) = (state.traction.get(*i).copied(), slips.get(a)) else { continue };
                if t == [0.0; 3] {
                    continue;
                }
                for (c, d) in dofs.iter().enumerate() {
                    let on: Vec3 = std::array::from_fn(|l| (0..3).map(|k| t[k] * jac[k].get(var(corners.len(), c, l)).copied().unwrap_or(0.0)).sum());
                    for l in 0..3 {
                        if let Some(x) = d[l].and_then(|d| internal.get_mut(d)) {
                            *x -= on[l];
                        }
                    }
                    if let Some(own) = corners.get(c).and_then(|j| state.nodal_force.get_mut(*j)) {
                        (0..3).for_each(|l| own[l] += on[l]);
                    }
                }
            }
        }
        // Jᵀ M J per closed node.
        for (i, reach) in jacobians {
            let Some(Some(m)) = middle.get(i) else { continue };
            let nodes: Vec<usize> = reach.keys().copied().collect();
            let blocks_of: Vec<&[[f64; 3]; 3]> = reach.values().collect();
            let width = 3 * nodes.len();
            let mut data = vec![0.0; width * width];
            for (a, ja) in blocks_of.iter().enumerate() {
                for (b, jb) in blocks_of.iter().enumerate() {
                    for l in 0..3 {
                        for r in 0..3 {
                            let mut v = 0.0;
                            for k in 0..3 {
                                for q in 0..3 {
                                    v += ja[k][l] * m[k][q] * jb[q][r];
                                }
                            }
                            if let Some(x) = data.get_mut((3 * a + l) * width + 3 * b + r) {
                                *x = v;
                            }
                        }
                    }
                }
            }
            blocks.push((nodes, data));
        }
        Ok(())
    }

    /// The result of a converged state.
    pub fn result(&self, mesh: &Mesh, state: &State, augmentations: usize) -> ContactResult {
        let status: Vec<Status> = state
            .pressure
            .iter()
            .enumerate()
            .map(|(i, p)| match (*p > 0.0, self.friction.is_some(), state.slipping.get(i).copied().unwrap_or(false)) {
                (false, ..) => Status::Open,
                (true, false, _) => Status::Closed,
                (true, true, false) => Status::Stick,
                (true, true, true) => Status::Slip,
            })
            .collect();
        let penetration = state.gap.iter().zip(&status).filter(|(_, s)| s.is_closed()).filter_map(|(g, _)| *g).fold(0.0f64, |m, g| m.max(-g));
        // The force on the secondary nodes pushes them back along −n; its tangential part is the
        // friction's.
        let normal: Vec<f64> = state.nodal_force.iter().zip(&state.report_normal).map(|(f, n)| -dot(*f, *n)).collect();
        let per_area = |v: f64, a: f64| if a > 0.0 { v / a } else { 0.0 };
        let pressure = normal.iter().zip(&state.support).map(|(f, a)| per_area(*f, *a)).collect();
        let shear = if self.friction.is_some() {
            state
                .nodal_force
                .iter()
                .zip(&state.report_normal)
                .zip(&state.support)
                .map(|((f, n), a)| {
                    let along = dot(*f, *n);
                    std::array::from_fn(|k| per_area(f[k] - along * n[k], *a))
                })
                .collect()
        } else {
            Vec::new()
        };
        let force = state.nodal_force.iter().fold([0.0; 3], |s, f| [s[0] + f[0], s[1] + f[1], s[2] + f[2]]);
        ContactResult {
            name: self.name.clone(),
            nodes: self.nodes.iter().map(|n| mesh.labels.get(*n).copied().unwrap_or(0)).collect(),
            pressure,
            multiplier: state.pressure.clone(),
            shear,
            slip: state.slip.clone(),
            gap: state.gap.clone(),
            status,
            force,
            normal_force: normal.iter().sum(),
            penetration,
            augmentations,
        }
    }

    /// How far the augmentation is from converged at `state`, computed with `multipliers`: the
    /// largest gap of a closed node, or the gap equivalent λ_A / ε_A of the multiplier an open
    /// node keeps, over the thickness of the element under the node.
    pub fn augmentation_error(&self, multipliers: &[f64], state: &State) -> f64 {
        let mut worst = 0.0f64;
        for i in 0..self.nodes.len() {
            let (Some(p), Some(e), Some(h)) = (state.pressure.get(i), self.penalty.get(i), self.thickness.get(i)) else { continue };
            let lambda = multipliers.get(i).copied().unwrap_or(0.0);
            let gap = if *p > 0.0 { state.gap.get(i).copied().flatten().unwrap_or(0.0).abs() } else { lambda.max(0.0) / e };
            worst = worst.max(gap / h);
        }
        worst
    }
}

#[cfg(test)]
mod report_tests;
