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
//! own faces that have a main-face pair: ĝ_A = g̃_A / A_A, a length (negative: penetration).
//! Each paired secondary face contributes its full reference area once, even if it overlaps
//! several main faces. Unpaired faces contribute nothing. These reference contributions are
//! constant within a fixed pairing, preserving the linearisation below; finite sliding selects
//! the contributing faces anew at each evaluation, including faces that first pair later.
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
//! over the paired secondary faces at the node (Ḡ the mean shear modulus of the element under the face,
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
//! **Node normal.** A secondary node's normal is the sum of its overlaps' projection normals
//! weighted by its support, Σ ∫ Φ_A n dA, normalised (in small sliding the reference pairs'
//! normals): only the faces that carry the node's contact turn it, so the uncontacted sides of a
//! surface's corner leave it alone, and a quadratic pair's facets each count with their own
//! normal. The friction return map's tangent plane and the reported pressure and shear use it.
//!
//! **Output.** Per secondary node the pressure as the node's normal contact force over the area
//! it carries (∫ Φ_A over the overlaps: smoother than the multipliers, which are the coefficients
//! of the pressure field and ring at the contact's edge), the multiplier, the gap and the status;
//! per contact the force on the secondary surface and the sum of the local normal components.
//! On a curved surface this scalar sum is not a component of the resultant along one fixed
//! direction: compare `force` with global equilibrium and `normal_force` with integrated pressure.
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
//! **Quadratic faces** (tri6 and quad8: tet10 and hex20 elements), on either side and against
//! linear ones (M. A. Puso, T. A. Laursen, J. Solberg, *A segment-to-segment mortar contact
//! method for quadratic elements and large deformations*, CMAME 197:555–566, 2008). Each face is
//! segmented as four flat facets through its nodes ([`fieldcraft_fem::mortar::integrate_faces`]),
//! with the quadratic shape functions evaluated at the facets' points, so a pair of faces is
//! integrated facet pair by facet pair as linear faces are; the points' positions are the
//! quadratic interpolation of the nodes. The secondary side's multiplier basis Φ is the
//! transformed one of [`fieldcraft_fem::mortar::multiplier`]: the standard quadratic functions
//! integrate to nothing at a tri6's corners and to a negative area at a quad8's, which would
//! leave those nodes no area to carry a pressure; Φ spans the same space, is a partition of
//! unity (so the patch test holds) and is non-negative, so every node has a positive area A_A and
//! its weighted gap is a positive average of the gap over its faces. The linearisation is still
//! exact: the automatic differentiation runs per facet pair over the facets' 24 vertex
//! coordinates, chained to the nodes ([`quadratic_terms`]). The reported pressure of a node of a
//! quadratic face is its nodal force in the multiplier basis over its area
//! ([`Interface::in_multiplier_basis`]).
//!
//! **Scope.** Solid faces (tri3, quad4, tri6 and quad8 faces: tet4, tet10, wedge6, hex8 and hex20
//! elements) in 3D models and linear edges (tri3 and quad4 elements) in 2D ones, with Coulomb
//! friction ([`Interface::friction`]). Quadratic edges and arc-length steps are refused by name,
//! as are quadratic faces in ties ([`crate::tie`]) and a mid-side node shared by a tri6 and a
//! quad8 face of a secondary surface (the two transform the basis differently there).

use std::collections::BTreeMap;

use fieldcraft_fem::contact::{Aabb, Grid, bounds};
use fieldcraft_fem::dual::{Dual, Hyper, PACKED, Real, VARS, packed};
use fieldcraft_fem::frame::{Vec3, dot, norm, unit};
use fieldcraft_fem::mortar::{
    FACE_NODES, FaceSample, Sample, alpha, dot as dot_t, edge_area_vector, face_area_vector, facet_points, facet_vertices, facets, integrate,
    integrate_edges, integrate_faces, integrate_facet, interpolate_fixed, reversed, sub, vertex_weights,
};
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
/// element's thickness keeps its penalty in the tangent through that increment's retries and
/// augmentations (its pressure, and so the residual, is still zero): a node that touches with no pressure, the edge
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
/// (and its retries/augmentations take [`ENGAGED_BAND`] until the increment is accepted).
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
    /// The secondary nodes (mesh slots), sorted.
    pub nodes: Vec<usize>,
    /// ∫ Φ_a over each whole secondary face in the reference configuration.
    reference_area: Vec<[f64; FACE_NODES]>,
    penalty_scale: f64,
    /// The scales over every secondary face, used only to release multipliers at unpaired nodes.
    all_faces: Normalization,
    /// Offset taken out of the weighted gap (initial overclosure).
    offset: Vec<f64>,
    /// Small sliding (steps without nonlinear geometry): the face pairs and their integrals,
    /// fixed in the reference configuration; `None`: finite sliding, searched and integrated
    /// on the current configuration at every evaluation.
    fixed: Option<Vec<Vec<Fixed>>>,
    /// Coulomb friction (frictionless: none).
    friction: Option<Friction>,
    /// A 2D model's out-of-plane measure (its surfaces are edges); `None`: solid faces.
    extent: Option<Extent>,
    /// The mid-side nodes of the quadratic secondary faces, the corners at their ends (indices in
    /// [`Interface::nodes`]) and their faces' α, for the nodal forces in the multiplier basis.
    mid_sides: Vec<((usize, usize, usize), f64)>,
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
    normalization: Normalization,
    /// ĝ_A, `None` where no main face lies over the node.
    pub gap: Vec<Option<f64>>,
    pub pressure: Vec<f64>,
    /// ∫ Φ_A over the overlaps (the area the pressure acts on).
    pub support: Vec<f64>,
    /// The force the contact puts on each secondary node, and the node's unit outward normal:
    /// its overlaps' projection normals weighted by its support, Σ ∫ Φ_A n normalised (the
    /// reference pairs' in small sliding), zero where nothing overlaps the node or the overlaps
    /// cancel. The friction law's tangent plane and the reported pressure and shear share it.
    pub nodal_force: Vec<Vec3>,
    pub normal: Vec<Vec3>,
    /// Friction: the tangential traction on each secondary node (on it, against its slip), its
    /// slip over the increment (relative to the main surface, normalised like the gap) and
    /// whether it slips; empty without friction.
    pub traction: Vec<Vec3>,
    pub slip: Vec<Vec3>,
    pub slipping: Vec<bool>,
}

/// Reference scales summed over the secondary faces that currently have main-face pairs.
#[derive(Debug, Clone, PartialEq, Default)]
struct Normalization {
    area: Vec<f64>,
    penalty: Vec<f64>,
    thickness: Vec<f64>,
}

impl Normalization {
    /// The sums over the secondary faces given with their reference integrals and node indices.
    fn over<'a>(n: usize, penalty_scale: f64, faces: impl Iterator<Item = (&'a Face, &'a [f64; FACE_NODES], &'a [usize])>) -> Normalization {
        let mut out = Normalization { area: vec![0.0; n], penalty: vec![0.0; n], thickness: vec![f64::INFINITY; n] };
        for (face, own, corners) in faces {
            for (i, share) in corners.iter().zip(own) {
                if let (Some(a), Some(e), Some(h)) = (out.area.get_mut(*i), out.penalty.get_mut(*i), out.thickness.get_mut(*i)) {
                    *a += share;
                    *e = e.max(penalty_scale * face.stiffness);
                    *h = h.min(face.thickness);
                }
            }
        }
        out
    }

    /// Node `i`'s (A_A, ε_A, h_A), `None` where none of the faces summed is at it.
    fn at(&self, i: usize) -> Option<(f64, f64, f64)> {
        let a = *self.area.get(i)?;
        (a > 0.0).then_some((a, *self.penalty.get(i)?, *self.thickness.get(i)?))
    }
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
    /// (∫ Φ_A over the overlaps), so that Σ p_A ∫Φ_A is `normal_force`.
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
    /// The vector sum of the forces the contact puts on the secondary surface.
    pub force: Vec3,
    /// Sum of the local normal components, Σ p_A ∫Φ_A. On a curved contact this need not
    /// equal the resultant's magnitude or its component in a fixed load direction.
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
    /// ∂g̃_a/∂x per secondary node a over the pair's nodes (3 per node, secondary nodes first).
    gradient: Vec<Vec<f64>>,
    /// −Σ_a p_a ∂²g̃_a over the pair's nodes, dense (tangent evaluations only).
    hessian: Option<Vec<f64>>,
}

/// The variable of a linear pair's AD ([`variables`]) for coordinate `k` of its node `c` (`s`
/// secondary corners, then the main ones).
fn var(s: usize, c: usize, k: usize) -> usize {
    if c < s { 3 * c + k } else { 12 + 3 * (c - s) + k }
}

/// A linear pair's gradient over its AD variables, over its nodes instead (`s` secondary
/// corners, `count` nodes).
fn by_node(g: &[f64; VARS], s: usize, count: usize) -> Vec<f64> {
    (0..3 * count).map(|i| g.get(var(s, i / 3, i % 3)).copied().unwrap_or(0.0)).collect()
}

/// The tangent blocks of a contact: per block its mesh slots and the dense matrix over their
/// translations (row-major, node by node).
pub type Blocks = Vec<(Vec<usize>, Vec<f64>)>;

fn position(mesh: &Mesh, u: &[f64], node: usize, dofs: &[Option<usize>; 3]) -> Vec3 {
    let x = mesh.xyz.get(node).copied().unwrap_or_default();
    [0, 1, 2].map(|k| x[k] + dofs[k].and_then(|d| u.get(d)).copied().unwrap_or(0.0))
}

/// The element's volume (hex8, tet4, and hex20, tet10 by their corners), or a 2D element's area (tri3, quad4), for the thickness
/// under a face or edge.
fn volume(kind: ElementKind, x: &[Vec3]) -> Option<f64> {
    match kind {
        ElementKind::Hex8 => <[Vec3; 8]>::try_from(x).ok().and_then(|a| fieldcraft_fem::uniform::Hex8::hex8(&a).ok()).map(|k| k.volume),
        ElementKind::Tet4 => <[Vec3; 4]>::try_from(x).ok().and_then(|a| fieldcraft_fem::uniform::Tet4::tet4(&a).ok()).map(|k| k.volume),
        // Quadratic elements: their corners' (straight-sided) volume.
        ElementKind::Hex20 => x.get(..8).and_then(|c| volume(ElementKind::Hex8, c)),
        ElementKind::Tet10 => x.get(..4).and_then(|c| volume(ElementKind::Tet4, c)),
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

/// The overlap integration of a linear face pair (`extent`: `None`) or an edge pair of a 2D model.
fn overlap<T: Real>(extent: Option<Extent>, xs: &[[T; 3]], xm: &[[T; 3]], f: impl FnMut(&Sample<T>)) -> Option<()> {
    match extent {
        None => integrate(xs, xm, f),
        Some(e) => integrate_edges(xs, xm, e, f),
    }
}

/// The overlap integration of a face pair of any order ([`integrate_faces`]: a linear pair's
/// points are [`overlap`]'s) or of an edge pair of a 2D model.
fn integrate_pair<T: Real>(extent: Option<Extent>, xs: &[[T; 3]], xm: &[[T; 3]], mut f: impl FnMut(&FaceSample<T>)) -> Option<()> {
    match extent {
        None => integrate_faces(xs, xm, f),
        Some(e) => integrate_edges(xs, xm, e, |p| f(&FaceSample::linear(p))),
    }
}

/// A face's area vector (a quadratic face's facets' sum), or an edge's ([`edge_area_vector`]) in
/// a 2D model.
fn area_of(extent: Option<Extent>, x: &[Vec3]) -> Option<Vec3> {
    match extent {
        None => face_area_vector(x),
        Some(e) => edge_area_vector(x, e),
    }
}

fn size(x: &[Vec3]) -> f64 {
    let b = bounds(x, 0.0);
    (0..3).map(|k| b[1][k] - b[0][k]).fold(0.0, f64::max)
}

/// The faces of set `set` for definition `owner` (``contact `name` ``, ``tie `name` ``), refused
/// unless they are faces of analysed solid elements, linear unless `quadratic`.
pub(crate) fn faces(model: &Model, resolved: &Resolved, mesh: &Mesh, owner: &str, set: &str, quadratic: bool) -> Result<Vec<Face>> {
    face_list(model, mesh, owner, resolved.sets.get(set).map(|s| s.faces.as_slice()).unwrap_or_default(), quadratic)
}

/// The edges of set `set` of a 2D model for definition `owner`, refused unless they are linear
/// edges of analysed 2D elements.
fn edges(model: &Model, resolved: &Resolved, mesh: &Mesh, owner: &str, set: &str) -> Result<Vec<Face>> {
    edge_list(model, mesh, owner, resolved.sets.get(set).map(|s| s.edges.as_slice()).unwrap_or_default())
}

/// The faces `list` for definition `owner`, as [`faces`].
fn face_list(model: &Model, mesh: &Mesh, owner: &str, list: &[ElementFace], quadratic: bool) -> Result<Vec<Face>> {
    list.iter()
        .map(|f| {
            let at = |what: String| Error::Invalid(format!("{owner}: element {} face {}: {what}", f.element, f.face));
            let labels = model.face_nodes(*f).ok_or_else(|| at("no such face".into()))?;
            surface_face(mesh, owner, f.element, labels, Surface::Faces { quadratic }, at)
        })
        .collect()
}

/// The edges `list` of a 2D model for definition `owner`, as [`edges`].
fn edge_list(model: &Model, mesh: &Mesh, owner: &str, list: &[ElementEdge]) -> Result<Vec<Face>> {
    list.iter()
        .map(|e| {
            let at = |what: String| Error::Invalid(format!("{owner}: element {} edge {}: {what}", e.element, e.edge));
            let labels = model.edge_nodes(*e).ok_or_else(|| at("no such edge".into()))?;
            surface_face(mesh, owner, e.element, labels, Surface::Edges, at)
        })
        .collect()
}

/// What a surface is made of: faces of solid elements (quadratic ones too, or only linear
/// ones), or linear edges of 2D elements.
#[derive(Debug, Clone, Copy)]
enum Surface {
    Faces { quadratic: bool },
    Edges,
}

/// One face (or one edge of a 2D element) of element `element` with nodes `labels`.
fn surface_face(mesh: &Mesh, owner: &str, element: u32, labels: Vec<u32>, surface: Surface, at: impl Fn(String) -> Error) -> Result<Face> {
    let el = mesh.element(element).ok_or_else(|| at("the element is not a solid or 2D element of the analysis".into()))?;
    let edge = matches!(surface, Surface::Edges);
    let (taken, takes) = match surface {
        Surface::Edges => (labels.len() == 2 && el.kind.dimension() == 2, "the linear edges of tri3 and quad4 elements"),
        Surface::Faces { quadratic: false } => {
            (matches!(labels.len(), 3 | 4) && el.kind.dimension() == 3, "the linear faces of tet4, wedge6 and hex8 elements")
        }
        Surface::Faces { quadratic: true } => {
            (matches!(labels.len(), 3 | 4 | 6 | 8) && el.kind.dimension() == 3, "the faces of tet4, tet10, wedge6, hex8 and hex20 elements")
        }
    };
    if !taken {
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
        _ => face_area_vector(&x).map_or(0.0, norm),
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
/// D_ai = ∫ Φ_a N_i (secondary nodes i) and M_aj = ∫ Φ_a N_j (main nodes j), and the
/// secondary face's normal there, so that g̃_a = n · (Σ_j M_aj x_j − Σ_i D_ai x_i).
#[derive(Debug, Clone)]
struct Fixed {
    main: usize,
    d: [[f64; FACE_NODES]; FACE_NODES],
    m: [[f64; FACE_NODES]; FACE_NODES],
    normal: Vec3,
    support: [f64; FACE_NODES],
}

impl Fixed {
    /// B_a, the slip operator ([`slip_operator`]) of the reference configuration, for a pair of
    /// `s` secondary and `t` main nodes.
    fn slip_jacobian(&self, a: usize, s: usize, t: usize) -> [Vec<f64>; 3] {
        let mut j = [0, 1, 2].map(|_| vec![0.0; 3 * (s + t)]);
        let (Some(d), Some(m)) = (self.d.get(a), self.m.get(a)) else { return j };
        let n = self.normal;
        for (k, row) in j.iter_mut().enumerate() {
            for l in 0..3 {
                let proj = if k == l { 1.0 } else { 0.0 } - n[k] * n[l];
                for (i, d) in d.iter().take(s).enumerate() {
                    row[3 * i + l] = d * proj;
                }
                for (i, m) in m.iter().take(t).enumerate() {
                    row[3 * (s + i) + l] = -m * proj;
                }
            }
        }
        j
    }

    /// ∂g̃_a/∂x over the pair's nodes (`s` secondary, then `t` main).
    fn gradient(&self, a: usize, s: usize, t: usize) -> Vec<f64> {
        let mut g = vec![0.0; 3 * (s + t)];
        let (Some(d), Some(m)) = (self.d.get(a), self.m.get(a)) else { return g };
        for k in 0..3 {
            for (i, d) in d.iter().take(s).enumerate() {
                g[3 * i + k] = -d * self.normal[k];
            }
            for (i, m) in m.iter().take(t).enumerate() {
                g[3 * (s + i) + k] = m * self.normal[k];
            }
        }
        g
    }
}

/// What a secondary face finds against the main surface: per main face, ∫ Φ_a g and ∫ Φ_a for
/// each node a.
struct Found {
    main: usize,
    gap: [f64; FACE_NODES],
    support: [f64; FACE_NODES],
    /// ∫ Φ_a n per node a: the support-weighted projection normal of the overlap (a quadratic
    /// pair's facets each with their own).
    normal: [Vec3; FACE_NODES],
}

/// Each secondary node's unit normal: the sum of its overlaps' support-weighted normals
/// Σ ∫ Φ_A n, so that only the faces carrying the node's contact turn it (the uncontacted sides
/// at a corner do not). Opposing supports can cancel: a zero sum has no direction, and the node
/// gets none (zero) rather than an invented one.
fn normals(pairs: &[Vec<Found>], corner_index: &[Vec<usize>], support: &[f64]) -> Vec<Vec3> {
    let mut normals = vec![[0.0; 3]; support.len()];
    for (found, corners) in pairs.iter().zip(corner_index) {
        for pair in found {
            for (i, weighted) in corners.iter().zip(&pair.normal) {
                let (Some(normal), Some(&total)) = (normals.get_mut(*i), support.get(*i)) else { continue };
                if total > 0.0 {
                    // Over the node's support before adding, so that large finite areas cannot
                    // overflow the sum: it is a direction.
                    (0..3).for_each(|k| normal[k] += weighted[k] / total);
                }
            }
        }
    }
    for normal in &mut normals {
        *normal = unit(*normal).unwrap_or_default();
    }
    normals
}

/// The weighted gaps g̃_a of one linear face pair (corners `xs`, `xm`) over any scalar (quadratic
/// pairs: [`quadratic_terms`]).
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

/// The slip operator of one face pair (nodes `xs`, `xm`): per secondary node a the matrix
/// B_a with w_a = B_a Δx, the weighted tangential slip ∫ Φ_a (I − n n)(Δx_s − Δx_m) dA of the
/// pair's nodes' motions Δx (rows: its components; columns: the pair's nodes' translations,
/// secondary nodes first), with the integrals on the current configuration: the rate form of
/// the mortar slip, B_a = [∫ Φ_a N_i P, −∫ Φ_a N_j P] (McDevitt & Laursen 2000; Puso & Laursen
/// 2004), P = I − n n.
fn slip_operator(extent: Option<Extent>, xs: &[Vec3], xm: &[Vec3]) -> Option<Vec<[Vec<f64>; 3]>> {
    let (s, t) = (xs.len(), xm.len());
    let mut b: Vec<[Vec<f64>; 3]> = (0..s).map(|_| [0, 1, 2].map(|_| vec![0.0; 3 * (s + t)])).collect();
    integrate_pair(extent, xs, xm, |p| {
        let n = p.normal;
        for (ba, phi) in b.iter_mut().zip(p.multiplier) {
            let wa = p.weight * phi;
            if wa == 0.0 {
                continue;
            }
            for (k, row) in ba.iter_mut().enumerate() {
                for l in 0..3 {
                    let proj = wa * (if k == l { 1.0 } else { 0.0 } - n[k] * n[l]);
                    for (c, nc) in p.secondary.iter().take(s).enumerate() {
                        row[3 * c + l] += proj * nc;
                    }
                    for (c, nc) in p.main.iter().take(t).enumerate() {
                        row[3 * (s + c) + l] -= proj * nc;
                    }
                }
            }
        }
    })?;
    Some(b)
}

/// What one facet pair's pass computes ([`facet_terms`]).
#[derive(Clone, Copy)]
enum Pass<'a> {
    /// The weighted gaps and their direct derivatives.
    Gaps,
    /// Those, and T for the nodes' pressures.
    GapsAndT(&'a [f64]),
    /// Ψ for the nodes' pressures.
    Psi(&'a [f64]),
}

/// One facet pair's integrals over any scalar, with the pair's nodes `xs`, `xm` held as constants
/// (the facets' vertices, in `vs`, `vm`, are the variables): per secondary node a the weighted gap
/// G_a = ∫ Φ_a (x_m − x_s) · n with x = Σ N_i x_i, and the direct derivatives of the G_a with
/// respect to the nodes, ∂G_a/∂x_i = ∓∫ Φ_a N_i n (values, nodes `xs` then `xm`); with the
/// pressures p, T_i = Σ_a p_a ∂G_a/∂x_i, or Ψ = Σ_a p_a G_a alone ([`Pass`]). `any`: whether the
/// facets overlap.
struct FacetTerms<T> {
    gaps: [T; FACE_NODES],
    direct: Vec<Vec<f64>>,
    psi: T,
    t: Vec<T>,
    any: bool,
}

fn facet_terms<T: Real>(
    (ns, fs, xs, vs): (usize, usize, &[Vec3], &[[T; 3]]),
    (nm, fm, xm, vm): (usize, usize, &[Vec3], &[[T; 3]]),
    pass: Pass,
) -> Option<FacetTerms<T>> {
    let zero = T::cst(0.0);
    let width = 3 * (ns + nm);
    let (rows, t) = match pass {
        Pass::Gaps => (ns, 0),
        Pass::GapsAndT(_) => (ns, width),
        Pass::Psi(_) => (0, 0),
    };
    let mut out = FacetTerms { gaps: [zero; FACE_NODES], direct: vec![vec![0.0; width]; rows], psi: zero, t: vec![zero; t], any: false };
    integrate_facet((ns, fs), vs, (nm, fm), vm, |p| {
        out.any = true;
        let gap = dot_t(sub(interpolate_fixed(xm, &p.main), interpolate_fixed(xs, &p.secondary)), p.normal);
        // ∓N_i per node, secondary then main.
        let shapes = || p.secondary.iter().take(ns).map(|n| -*n).chain(p.main.iter().take(nm).copied());
        let pressure_weight = |pressure: &[f64]| p.weight * p.multiplier.iter().zip(pressure).fold(zero, |s, (m, pa)| s + *m * *pa);
        match pass {
            Pass::Psi(pressure) => out.psi = out.psi + pressure_weight(pressure) * gap,
            Pass::Gaps | Pass::GapsAndT(_) => {
                for ((g, direct), phi) in out.gaps.iter_mut().zip(out.direct.iter_mut()).zip(p.multiplier) {
                    let wa = p.weight * phi;
                    *g = *g + wa * gap;
                    for (c, n) in shapes().enumerate() {
                        let v = wa.val() * n.val();
                        for k in 0..3 {
                            if let Some(d) = direct.get_mut(3 * c + k) {
                                *d += v * p.normal[k].val();
                            }
                        }
                    }
                }
                if let Pass::GapsAndT(pressure) = pass {
                    let wp = pressure_weight(pressure);
                    for (c, n) in shapes().enumerate() {
                        let wn = wp * n;
                        for k in 0..3 {
                            if let Some(t) = out.t.get_mut(3 * c + k) {
                                *t = *t + wn * p.normal[k];
                            }
                        }
                    }
                }
            }
        }
    })?;
    Some(out)
}

/// ∂g̃_a/∂x per secondary node a over the pair's nodes (3 per node, secondary first) and, with the
/// nodes' `pressure`, the block −Σ_a p_a ∂²g̃_a, of a pair with a quadratic face (nodes `xs`,
/// `xm`). Per facet pair the automatic differentiation runs over the facets' vertices V (24
/// variables, as a linear pair's corners), which are fixed combinations V = W x of the nodes
/// ([`vertex_weights`]), with the nodes x in the quadratic interpolation of the points held
/// apart: g̃_a = G_a(V(x), x), linear in its second argument, so that
///
///   ∂g̃_a/∂x = ∂G_a/∂x + Wᵀ ∂G_a/∂V,
///   ∂²Ψ/∂x² = Wᵀ Ψ_VV W + Wᵀ ∂T/∂V + (∂T/∂V)ᵀ W   (Ψ = Σ_a p_a g̃_a, T = ∂Ψ/∂x at fixed V),
///
/// exact, at a few times a linear pair's cost per facet pair.
#[allow(clippy::type_complexity)]
fn quadratic_terms(xs: &[Vec3], xm: &[Vec3], pressure: Option<&[f64]>) -> Option<(Vec<Vec<f64>>, Option<Vec<f64>>)> {
    let (ns, nm) = (xs.len(), xm.len());
    let width = 3 * (ns + nm);
    let (ps, pm) = (facet_points(xs)?, facet_points(xm)?);
    let mains: Vec<(&[usize], Vec<Vec3>)> = facets(nm)?.iter().map(|m| (*m, facet_vertices(&pm, m))).collect();
    let mut gradient = vec![vec![0.0; width]; ns];
    let mut hessian = pressure.map(|_| vec![0.0; width * width]);
    // No pressure, no second derivatives.
    let pressed = pressure.filter(|p| p.iter().any(|v| *v != 0.0));
    let mut row = vec![0.0; width];
    for (i, fs) in facets(ns)?.iter().enumerate() {
        let vs = facet_vertices(&ps, fs);
        for (j, (fm, vm)) in mains.iter().enumerate() {
            // Per facet vertex: its first variable and the pair's nodes it combines.
            let slots: Vec<(usize, &[(usize, f64)], usize)> = fs
                .iter()
                .enumerate()
                .map(|(k, v)| (3 * k, vertex_weights(ns, *v), 0))
                .chain(fm.iter().enumerate().map(|(k, v)| (12 + 3 * k, vertex_weights(nm, *v), ns)))
                .collect();
            // Wᵀ applied to a gradient over the facets' variables, added to `out` (over the nodes).
            let chain = |out: &mut [f64], g: &[f64; VARS], scale: f64| {
                for &(first, weights, offset) in &slots {
                    for &(node, w) in weights {
                        for k in 0..3 {
                            if let (Some(o), Some(v)) = (out.get_mut(3 * (offset + node) + k), g.get(first + k)) {
                                *o += scale * w * v;
                            }
                        }
                    }
                }
            };
            let (a, b) = variables(&vs, vm, Dual::var);
            let pass = pressed.map_or(Pass::Gaps, Pass::GapsAndT);
            let terms = facet_terms((ns, i, xs, &a), (nm, j, xm, &b), pass)?;
            if !terms.any {
                continue;
            }
            for ((out, g), direct) in gradient.iter_mut().zip(&terms.gaps).zip(&terms.direct) {
                out.iter_mut().zip(direct).for_each(|(o, d)| *o += d);
                chain(out, &g.g, 1.0);
            }
            let (Some(h), Some(p)) = (hessian.as_mut(), pressed) else { continue };
            // −Wᵀ ∂T/∂V and its transpose.
            for (r, t) in terms.t.iter().enumerate() {
                row.fill(0.0);
                chain(&mut row, &t.g, -1.0);
                for (c, v) in row.iter().enumerate().filter(|(_, v)| **v != 0.0) {
                    if let Some(x) = h.get_mut(r * width + c) {
                        *x += v;
                    }
                    if let Some(x) = h.get_mut(c * width + r) {
                        *x += v;
                    }
                }
            }
            // −Wᵀ Ψ_VV W.
            let (a, b) = variables(&vs, vm, Hyper::var);
            let psi = facet_terms((ns, i, xs, &a), (nm, j, xm, &b), Pass::Psi(p))?.psi;
            for &(f1, w1, o1) in &slots {
                for &(f2, w2, o2) in &slots {
                    for k1 in 0..3 {
                        for k2 in 0..3 {
                            let v = psi.h.get(packed(f1 + k1, f2 + k2)).copied().unwrap_or(0.0);
                            if v == 0.0 {
                                continue;
                            }
                            for &(n1, a1) in w1 {
                                for &(n2, a2) in w2 {
                                    if let Some(x) = h.get_mut((3 * (o1 + n1) + k1) * width + 3 * (o2 + n2) + k2) {
                                        *x -= a1 * a2 * v;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    Some((gradient, hessian))
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
            (Some(s), None) => faces(model, resolved, mesh, &owner, s, true),
            (Some(s), Some(_)) => edges(model, resolved, mesh, &owner, s),
            (None, None) => face_list(model, mesh, &owner, &model.exterior_faces(), true),
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
        let mut reference_area = Vec::with_capacity(secondary.len());
        for f in &secondary {
            let x: Vec<Vec3> = f.nodes.iter().map(|n| mesh.xyz.get(*n).copied().unwrap_or_default()).collect();
            // ∫ Φ_a over the face itself (the face against its own mirror image).
            let mut own = [0.0; FACE_NODES];
            integrate_pair(extent, &x, &reversed(&x), |p| {
                for (o, phi) in own.iter_mut().zip(p.multiplier) {
                    *o += p.weight * phi;
                }
            });
            reference_area.push(own);
        }
        let corners: Vec<Vec<usize>> = secondary.iter().map(|f| f.nodes.iter().filter_map(|x| nodes.binary_search(x).ok()).collect()).collect();
        let faces = secondary.iter().zip(&reference_area).zip(&corners).map(|((f, own), c)| (f, own, c.as_slice()));
        let all_faces = Normalization::over(nodes.len(), contact.penalty_scale, faces);
        let label_of = |i: usize| nodes.get(i).and_then(|n| mesh.labels.get(*n)).copied().unwrap_or(0);
        if let Some(i) = all_faces.area.iter().position(|a| !(*a > 0.0 && a.is_finite())) {
            return Err(Error::Invalid(format!("contact `{name}`: the secondary faces at node {} have no area", label_of(i))));
        }
        let offset = vec![0.0; nodes.len()];
        // Each mid-side node with the node count of the faces it is on.
        let mut mid_sides: BTreeMap<(usize, usize, usize), usize> = BTreeMap::new();
        for f in secondary.iter().filter(|f| f.nodes.len() > 4) {
            let index = |k: usize| f.nodes.get(k).and_then(|n| nodes.binary_search(n).ok());
            let c = f.nodes.len() / 2;
            for k in 0..c {
                // Each edge once: two faces see it in opposite directions.
                if let (Some(m), Some(a), Some(b)) = (index(c + k), index(k), index((k + 1) % c)) {
                    let count = *mid_sides.entry((m, a.min(b), a.max(b))).or_insert(f.nodes.len());
                    // A tri6 and a quad8 transform the basis by different α at the node: no
                    // single basis holds there, so its pressure would be wrong.
                    if count != f.nodes.len() {
                        return Err(Error::NotSupported(format!(
                            "contact `{name}`: mid-side node {} is shared by a tri6 face and a quad8 face of the secondary surface (their multiplier bases differ there); mesh that surface with one kind of quadratic face, or make the other surface secondary",
                            label_of(m)
                        )));
                    }
                }
            }
        }
        Ok(Interface {
            name: name.to_string(),
            secondary,
            main,
            nodes,
            reference_area,
            penalty_scale: contact.penalty_scale,
            all_faces,
            offset,
            fixed: None,
            friction: (!contact.friction.is_frictionless()).then_some(contact.friction),
            extent,
            mid_sides: mid_sides.into_iter().map(|(edge, count)| (edge, alpha(count))).collect(),
        })
    }

    /// Adjust (offset) or refuse the nodes that start behind the main surface.
    fn overclosure(mut self, mesh: &Mesh, mode: InitialPenetration, warnings: &mut Vec<String>) -> Result<Interface> {
        let zero = vec![0.0; mesh.dofs()];
        let pairs = self.pairs(mesh, &zero)?;
        let (gaps, _, corner_index) = self.sums(&pairs);
        let normalization = self.normalization(&pairs, &corner_index);
        // The model's extent (relative, so the test reads the same in any unit of length).
        let scale = mesh.xyz.iter().flatten().fold(0.0f64, |m, c| m.max(c.abs()));
        let (mut count, mut deepest) = (0usize, (0.0f64, 0u32));
        for (i, g) in gaps.iter().enumerate() {
            let (Some(g), Some((a, _, _))) = (g, normalization.at(i)) else { continue };
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
                    Found { main: p.main, gap, support: p.support, normal: p.support.map(|a| p.normal.map(|n| a * n)) }
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
                    let (mut d, mut m, mut normal) = ([[0.0; FACE_NODES]; FACE_NODES], [[0.0; FACE_NODES]; FACE_NODES], [0.0; 3]);
                    let mut mean = [0.0; 3];
                    integrate_pair(self.extent, xs, xm, |p| {
                        normal = p.normal;
                        (0..3).for_each(|k| mean[k] += p.weight * p.normal[k]);
                        for a in 0..FACE_NODES {
                            for i in 0..FACE_NODES {
                                d[a][i] += p.weight * p.multiplier[a] * p.secondary[i];
                                m[a][i] += p.weight * p.multiplier[a] * p.main[i];
                            }
                        }
                    })?;
                    // A quadratic pair's facets may face apart (a curved face): their mean normal.
                    if xs.len() > 4 || xm.len() > 4 {
                        normal = unit(mean).unwrap_or(normal);
                    }
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
                let (mut gap, mut support, mut normal) = ([0.0; FACE_NODES], [0.0; FACE_NODES], [[0.0; 3]; FACE_NODES]);
                let (mut deepest, mut any) = (f64::INFINITY, false);
                integrate_pair(self.extent, &xs, xm, |p| {
                    any = true;
                    let g = p.gap();
                    deepest = deepest.min(g);
                    for a in 0..FACE_NODES {
                        gap[a] += p.weight * g * p.multiplier[a];
                        let w = p.weight * p.multiplier[a];
                        support[a] += w;
                        (0..3).for_each(|k| normal[a][k] += w * p.normal[k]);
                    }
                });
                if any && deepest >= -m.depth && gap.iter().chain(&support).chain(normal.iter().flatten()).all(|v| v.is_finite()) {
                    found.push(Found { main: c, gap, support, normal });
                }
            }
            keep(&xs, &main_x, found)
        }))
    }

    /// Each paired secondary face contributes once, independently of how many main faces
    /// overlap it. Using its reference integral keeps these scales constant within a pairing.
    fn normalization(&self, pairs: &[Vec<Found>], corner_index: &[Vec<usize>]) -> Normalization {
        let faces = self.secondary.iter().zip(&self.reference_area).zip(corner_index).zip(pairs);
        let paired = faces.filter(|(_, found)| !found.is_empty()).map(|(((f, own), c), _)| (f, own, c.as_slice()));
        Normalization::over(self.nodes.len(), self.penalty_scale, paired)
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
        let normal = normals(&pairs, &corner_index, &support);
        let n = self.nodes.len();
        let mut state = State {
            normalization: self.normalization(&pairs, &corner_index),
            gap: vec![None; n],
            pressure: vec![0.0; n],
            support,
            nodal_force: vec![[0.0; 3]; n],
            normal,
            traction: Vec::new(),
            slip: Vec::new(),
            slipping: Vec::new(),
        };
        // Nodes closed, or open by less than `band` of their element's thickness: their penalty
        // stays in the tangent.
        let mut engaged = vec![false; n];
        for i in 0..self.nodes.len() {
            let (Some(Some(g)), Some((a, e, h))) = (weighted.get(i), state.normalization.at(i)) else { continue };
            let gap = g / a;
            let lambda = multipliers.get(i).copied().unwrap_or(0.0);
            if let (Some(slot), Some(p)) = (state.gap.get_mut(i), state.pressure.get_mut(i)) {
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
            let (ns, count) = (corners.len(), xs.len() + xm.len());
            let (gradient, hessian) = if let Some(fixed) = fixed {
                // Small sliding: linear in the positions.
                ((0..ns).map(|a| fixed.gradient(a, ns, xm.len())).collect(), None)
            } else if xs.len() > 4 || xm.len() > 4 {
                let Some(terms) = quadratic_terms(&xs, &xm, tangent.then_some(pressure.as_slice())) else { return Ok(None) };
                terms
            } else if tangent {
                let (vs, vm) = variables(&xs, &xm, Hyper::var);
                let Some(g) = weighted_gaps(self.extent, &vs, &vm) else { return Ok(None) };
                let mut h = vec![0.0; PACKED];
                for (ga, p) in g.iter().zip(&pressure) {
                    h.iter_mut().zip(&ga.h).for_each(|(x, y)| *x += p * y);
                }
                // −Σ_a p_a ∂²g̃_a over the pair's nodes.
                let width = 3 * count;
                let mut data = vec![0.0; width * width];
                for ci in 0..count {
                    for ki in 0..3 {
                        for cj in 0..count {
                            for kj in 0..3 {
                                let at = packed(var(ns, ci, ki), var(ns, cj, kj));
                                if let Some(x) = data.get_mut((3 * ci + ki) * width + 3 * cj + kj) {
                                    *x = -h.get(at).copied().unwrap_or(0.0);
                                }
                            }
                        }
                    }
                }
                (g.iter().take(ns).map(|ga| by_node(&ga.g, ns, count)).collect::<Vec<_>>(), Some(data))
            } else {
                let (vs, vm) = variables(&xs, &xm, Dual::var);
                let Some(g) = weighted_gaps(self.extent, &vs, &vm) else { return Ok(None) };
                (g.iter().take(ns).map(|ga| by_node(&ga.g, ns, count)).collect(), None)
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
                    let gv = [0, 1, 2].map(|k| g.get(3 * c + k).copied().unwrap_or(0.0));
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
            if let Some(data) = t.hessian {
                blocks.push((t.nodes, data));
            }
        }
        // (ε_A / A_A) ∂g̃_A ⊗ ∂g̃_A per closed node.
        for (i, gradient) in node_gradients {
            let Some((a, e, _)) = state.normalization.at(i) else { continue };
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
    /// ĝ_T = w_A / A_A; the trial traction t* = P t_n − ε_A ĝ_T (P projects the committed
    /// traction on the tangent plane of the node normal);
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
            let operator: Vec<[Vec<f64>; 3]> = match self.fixed.as_ref().and_then(|f| f.get(s)).and_then(|l| l.get(k)) {
                Some(fixed) => (0..corners.len()).map(|a| fixed.slip_jacobian(a, xs.len(), xm.len())).collect(),
                None => slip_operator(self.extent, &xs, &xm)?.into_iter().take(corners.len()).collect(),
            };
            let delta: Vec<f64> = xs.iter().zip(&rs).chain(xm.iter().zip(&rm)).flat_map(|(x, r)| [0, 1, 2].map(|c| x[c] - r[c])).collect();
            let slips: Vec<(Vec3, [Vec<f64>; 3])> =
                operator.into_iter().map(|b| ([0, 1, 2].map(|k| b[k].iter().zip(&delta).map(|(x, y)| x * y).sum()), b)).collect();
            let nodes: Vec<usize> = f.nodes.iter().chain(&m.nodes).copied().collect();
            let dofs: Vec<[Option<usize>; 3]> = f.dofs.iter().chain(&m.dofs).copied().collect();
            Some((nodes, dofs, corners.clone(), slips))
        });
        let terms: Vec<_> = terms.into_iter().flatten().collect();
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
                                block[k][l] += jac[k].get(3 * c + l).copied().unwrap_or(0.0);
                            }
                        }
                    }
                }
            }
        }
        // The return map, and per node the tangent's middle matrix.
        let mut middle: Vec<Option<[[f64; 3]; 3]>> = vec![None; n];
        for i in 0..n {
            let (Some(&p), Some((area, eps, _)), Some(&normal)) = (state.pressure.get(i), state.normalization.at(i), state.normal.get(i)) else {
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
                    let on: Vec3 = std::array::from_fn(|l| (0..3).map(|k| t[k] * jac[k].get(3 * c + l).copied().unwrap_or(0.0)).sum());
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
        let nodal = self.in_multiplier_basis(&state.nodal_force);
        let normal: Vec<f64> = nodal.iter().zip(&state.normal).map(|(f, n)| -dot(*f, *n)).collect();
        let per_area = |v: f64, a: f64| if a > 0.0 { v / a } else { 0.0 };
        let pressure = normal.iter().zip(&state.support).map(|(f, a)| per_area(*f, *a)).collect();
        let shear = if self.friction.is_some() {
            nodal
                .iter()
                .zip(&state.normal)
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

    /// Nodal forces `forces` of the secondary nodes in the multiplier basis: F̃ = T F for Φ = T N
    /// ([`fieldcraft_fem::mortar::multiplier`]), so that a pressure field p gives F̃_A = ∫ p Φ_A n
    /// and F̃_A / ∫ Φ_A is a pressure at every node. A quadratic face's own nodal forces are not:
    /// a uniform pressure gives a tri6's corners none and a quad8's a pull. Linear faces' forces
    /// are their own (T = I).
    fn in_multiplier_basis(&self, forces: &[Vec3]) -> Vec<Vec3> {
        let mut out = forces.to_vec();
        for &((m, a, b), alpha) in &self.mid_sides {
            let f = forces.get(m).copied().unwrap_or_default();
            for c in [a, b] {
                if let Some(o) = out.get_mut(c) {
                    (0..3).for_each(|k| o[k] += alpha * f[k]);
                }
            }
            if let Some(o) = out.get_mut(m) {
                *o = f.map(|v| v * (1.0 - 2.0 * alpha));
            }
        }
        out
    }

    /// How far the augmentation is from converged at `state`, computed with `multipliers`: the
    /// largest gap of a closed node, or the gap equivalent λ_A / ε_A of the multiplier an open
    /// node keeps, over the thickness of the element under the node.
    pub fn augmentation_error(&self, multipliers: &[f64], state: &State) -> f64 {
        let mut worst = 0.0f64;
        for i in 0..self.nodes.len() {
            let (Some(p), Some((_, e, h))) = (state.pressure.get(i), state.normalization.at(i).or_else(|| self.all_faces.at(i))) else { continue };
            let lambda = multipliers.get(i).copied().unwrap_or(0.0);
            let gap = if *p > 0.0 { state.gap.get(i).copied().flatten().unwrap_or(0.0).abs() } else { lambda.max(0.0) / e };
            worst = worst.max(gap / h);
        }
        worst
    }
}

#[cfg(test)]
mod normal_tests;

#[cfg(test)]
mod tests {
    use super::*;

    /// The weighted gaps of a pair of any order, by value.
    fn gaps(xs: &[Vec3], xm: &[Vec3]) -> [f64; FACE_NODES] {
        let mut g = [0.0; FACE_NODES];
        integrate_faces(xs, xm, |p| {
            for (g, phi) in g.iter_mut().zip(p.multiplier) {
                *g += p.weight * p.gap() * phi;
            }
        })
        .unwrap();
        g
    }

    fn quadratic_face(corners: &[Vec3], bulge: f64) -> Vec<Vec3> {
        let n = corners.len();
        let mids = (0..n).map(|k| {
            let m: Vec3 = std::array::from_fn(|i| 0.5 * (corners[k][i] + corners[(k + 1) % n][i]));
            [m[0], m[1], m[2] + bulge * (k as f64 + 1.0)]
        });
        corners.iter().copied().chain(mids).collect()
    }

    /// The chain rule of a quadratic pair's gradient and Hessian through its facets, against
    /// central differences of the weighted gaps, for curved, tilted faces of each order.
    #[test]
    fn quadratic_pairs_differentiate_exactly() {
        let quad = quadratic_face(&[[0.0, 0.0, 0.0], [1.0, 0.05, 0.02], [1.1, 0.9, -0.03], [-0.05, 1.0, 0.01]], 0.01);
        let tri = quadratic_face(&[[0.0, 0.0, 0.0], [1.1, 0.1, 0.01], [0.2, 1.0, -0.02]], -0.01);
        let down = |x: &[Vec3]| -> Vec<Vec3> { reversed(x).iter().map(|p| [p[0] + 0.13, p[1] - 0.07, 0.08 - 0.5 * p[2]]).collect() };
        let linear: Vec<Vec3> = vec![[-0.2, -0.3, 0.12], [-0.25, 1.15, 0.08], [1.2, 1.1, 0.1], [0.9, -0.2, 0.15]];
        for (xs, xm) in [
            (&quad, down(&quad)),
            (&tri, down(&quad)),
            (&quad, down(&tri)),
            (&quad, linear.clone()),
            (&quad[..4].to_vec(), down(&quad)),
            (&tri[..3].to_vec(), down(&tri)),
        ] {
            let (ns, nm) = (xs.len(), xm.len());
            let p: Vec<f64> = (0..ns).map(|a| 1.0 + 0.1 * a as f64).collect();
            let (gradient, hessian) = quadratic_terms(xs, &xm, Some(&p)).unwrap();
            let (plain, none) = quadratic_terms(xs, &xm, None).unwrap();
            assert!(none.is_none() && plain.iter().flatten().zip(gradient.iter().flatten()).all(|(a, b)| (a - b).abs() < 1e-14));
            let hessian = hessian.unwrap();
            let width = 3 * (ns + nm);
            let shifted = |v: usize, by: f64| {
                let (mut s, mut m) = (xs.clone(), xm.clone());
                let (node, k) = (v / 3, v % 3);
                if node < ns {
                    s[node][k] += by
                } else {
                    m[node - ns][k] += by
                }
                (s, m)
            };
            let h = 1e-6;
            for v in 0..width {
                let ((sp, mp), (sm, mm)) = (shifted(v, h), shifted(v, -h));
                let (gp, gm) = (gaps(&sp, &mp), gaps(&sm, &mm));
                for a in 0..ns {
                    let fd = (gp[a] - gm[a]) / (2.0 * h);
                    assert!((gradient[a][v] - fd).abs() < 1e-7, "{ns}/{nm}: ∂g_{a}/∂x_{v}: {} vs {fd}", gradient[a][v]);
                }
                // The Hessian's row v against differences of the gradient.
                let (up, _) = quadratic_terms(&sp, &mp, None).unwrap();
                let (dn, _) = quadratic_terms(&sm, &mm, None).unwrap();
                for w in 0..width {
                    let fd: f64 = (0..ns).map(|a| p[a] * (up[a][w] - dn[a][w]) / (2.0 * h)).sum();
                    assert!((hessian[v * width + w] + fd).abs() < 1e-6, "{ns}/{nm}: H({v}, {w}): {} vs {}", -hessian[v * width + w], fd);
                }
            }
        }
    }
}
