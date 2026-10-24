//! Nonlinear static analysis: equilibrium followed in load increments with Newton–Raphson
//! iterations.
//!
//! **Load path.** The step's load factor λ runs from 0 to 1. The dead loads (forces, tractions,
//! gravity, body forces; pressure too without nonlinear geometry) are λ times their values, a
//! follower pressure λ p on the deformed faces, and every prescribed displacement λ times its
//! value (u = T_f q + λ T_p u_p, the constraint elimination of the linear solver,
//! [`crate::constraints`]).
//!
//! **Newton–Raphson** with the full consistent tangent, refactored every iteration (Bathe
//! 2014, §8.4; M. A. Crisfield, *Non-linear Finite Element Analysis of Solids and Structures*,
//! vol. 1, Wiley 1991, ch. 9): K_T δq = r with r = T_fᵀ (λ f_ext(u) − f_int(u)), then, when the
//! step asks for it (it is off by default: on large elastic rotations the full Newton step
//! converges in fewer iterations), a line search on the work of the residual along δq, G(s) = δq · r(q + s δq): s = 1 is taken when
//! |G(1)| ≤ 0.8 |G(0)| (always near the solution, so quadratic convergence is kept), otherwise
//! s is moved to the root of the line through (0, G(0)) and (s, G(s)), within [0.1, 1], at most
//! five times (Crisfield 1991, §9.3). An increment has converged when, at once,
//! - **force**: the largest residual force is at most `force` times the force scale, the largest
//!   force the step has seen (applied, internal, reaction), never less than a round-off floor of
//!   10⁻⁸ of a typical diagonal stiffness times the model's size; and the largest residual
//!   moment is at most `force` times the moment scale, the largest moment the step has seen, never
//!   less than the force scale times the model's size L (the diagonal of its bounding box);
//! - **displacement**: the largest translation of the correction is at most `displacement` times
//!   the largest translation of the increment (at least 10⁻¹² L), and its largest rotation at
//!   most `displacement` times the largest rotation of the increment, never less than that
//!   translation over L (the rotation that moves a point at the model's size by it);
//! - **energy**: the work of the correction against the residual, s δq · r, is at most `energy`
//!   times that of the first iteration (K.-J. Bathe, A. P. Cimento, *Some practical procedures
//!   for the solution of nonlinear finite element equations*, CMAME 22:59–85, 1980);
//!
//! or when the force residual is a thousandth of its tolerance (round-off: no correction can
//! improve it). With path-dependent materials, an increment that starts within the force
//! tolerance (a predictor that already solved it) waives the energy criterion, whose reference
//! would be round-off.
//!
//! **Norms by kind of freedom.** Rotations (today those of coupled reference points; shells and
//! beams are refused) are radians and their residuals moments, so neither shares a norm with the
//! translations (metres, newtons): a mixed norm would change with the unit system and with the
//! part's size, since a moment scales as a force times a length and a rotation not at all. Each
//! kind is measured against its own scale, as the field-by-field checks of nonlinear finite
//! element practice do, and the measures recorded per iteration are the larger of the two. The
//! work of the energy criterion is consistent across kinds (N·m either way) and stays one sum.
//! A model without rotations converges exactly as with a single norm.
//!
//! **Automatic incrementation.** Each increment starts from the linear extrapolation of the
//! last converged one (the free motion scaled by the ratio of their sizes), which keeps the
//! inner nodes in step with the prescribed ones; from the last converged state when the
//! extrapolation inverts an element. The first increment is the step's `initial_increment`; two
//! increments in a row that converge in at most four iterations let the next grow by half (up
//! to `max_increment`), one that needs more than eight shrinks it by a quarter. An increment is
//! cut back to a quarter of its size when it fails: an element inverts, the residual is not
//! finite, the residual grows in two successive iterations, the tangent is singular, or it does
//! not converge in `max_iterations`. The
//! step fails cleanly, naming the load factor and the last reason, when an increment falls below
//! `min_increment`, when it has taken `max_increments` increments, or when the Newton iterations
//! of the whole step reach [`crate::Limits::max_newton_iterations`].
//!
//! **Indefinite tangents.** A tangent that is not positive definite — a strained mixed element
//! of a nearly incompressible material, an unstable material, a state past a limit or
//! bifurcation point — is factored by the symmetric indefinite LDLᵀ ([`crate::ldlt`]) and Newton
//! goes on: the residual decides convergence, and every iteration records the tangent's negative
//! eigenvalues (its inertia), which the step summary reports when a converged state keeps some
//! (an unstable equilibrium reached under load control). Only the first tangent of a model that
//! can move without straining is an error (the undeformed model's Cholesky factor names where);
//! a singular tangent later is a cutback.
//!
//! **Follower pressure** and the Cholesky factorization: the load stiffness of a pressure on the
//! deformed surface is not symmetric (Schweizerhof & Ramm 1984); its skew part is left out of the
//! tangent (the symmetric part is kept), so the factorization stays the supernodal Cholesky of the
//! linear solver. The residual is always exact, so the converged solution is unchanged; only the
//! rate of convergence can drop where the skew part matters. Over a closed surface, or one whose
//! boundary is held, the skew parts of neighbouring faces cancel and nothing is left out; the
//! result reports the skew part's size relative to the load stiffness when it is not small.
//!
//! **Elements.** Solids (tet4, tet10, wedge6, hex8, hex20 full and reduced, hybrid hex8, hex20
//! and tet10) and 2D plane strain and axisymmetric elements, total Lagrangian
//! ([`fieldcraft_fem::finite`]); without nonlinear geometry an element of a plastic material
//! takes the same loops with small strain, and an elastic one its linear stiffness (K u).
//! Springs act along their fixed directions. Trusses, shells and beams (large rotations), plane
//! stress, incompatible modes (with nonlinear geometry, or with an inelastic material) and
//! centrifugal loads are refused by name.
//!
//! **Thermal strain.** The step's temperatures are its temperature loads, or, when it has
//! temperature boundary conditions or heat loads, their steady conduction solved first (one-way
//! coupling, as the linear coupled temperature–displacement step). They ramp with the load
//! factor from each material's reference temperature, ΔT(λ) = λ (T − T_ref). With nonlinear
//! geometry the deformation gradient splits multiplicatively, F = F_m F_θ with
//! F_θ = diag(1 + αᵢ ΔT), the material sees F_m and the stress is J_θ τ_m per reference volume
//! (fem's `Temperature`); with small strain the split is additive, ε_m = ε − α ΔT; elements on
//! their linear stiffness take the equivalent nodal forces, a dead load. A Johnson–Cook material
//! softens at the point's temperature.
//!
//! **Plasticity** (`Material::plastic`, J2 with the radial return explicit dynamics uses,
//! [`fieldcraft_fem::elastoplastic`]): Simo's exponential return map at finite strain, the
//! additive radial return at small strain, both with their consistent tangents, so Newton keeps
//! its quadratic rate. The plastic history is committed when an increment converges and
//! discarded when it is cut back. A static step has no time, so the Johnson–Cook strain-rate term
//! is refused (C = 0 is the rate-independent law), and kinematic hardening at finite strain is
//! refused (its back stress breaks the principal-frame return); damage, erosion and adiabatic
//! heating act in explicit dynamics only, with a warning. Results add the equivalent plastic
//! strain at the integration points.
//! Couplings stay the linear (small-rotation) constraints of the linear solver, with a warning.
//!
//! **Contact** ([`crate::contact`]): the model's contact pairs add their mortar forces to the
//! internal forces and their consistent tangent to the stiffness, whose pattern takes the blocks
//! of the closed nodes' face pairs at every evaluation. Each converged increment is solved again
//! with updated multipliers until the gaps are tight (the augmented Lagrangian's Uzawa loop, the
//! increment's `augmentations`, each iteration tagged with its `augmentation`); an increment whose
//! active set chatters is cut back. With contact the force scale of the convergence test is that
//! of the converged states, the loads and the current state, not the largest force of every
//! trial: a trial state can press a body deep into another before the active set settles.
//!
//! Element forces and tangents are computed in parallel and added in label order, so results are
//! the same bit for bit on any number of threads.

use fieldcraft_fem::constitutive::MaterialModel;
use fieldcraft_fem::dense::Matrix;
use fieldcraft_fem::element::{Extrapolation, extrapolation};
use fieldcraft_fem::finite::{Kinematics, Temperature};
use fieldcraft_fem::plane::Theory;
use fieldcraft_fem::quadrature::Rule;
use fieldcraft_fem::{DofKind, DofSet, ElementType, PointStress, face, von_mises};
use fieldcraft_model::{
    ArcConstraint, ArcLength, BoundaryCondition, Extent, Hardening, Imperfection, Load, Model, NonlinearStatic, Plastic, Procedure, SectionKind, Step,
};
use fieldcraft_parallel as parallel;
use serde::Serialize;

use crate::analysis::factor_symmetric;
use crate::analysis::{
    ElementStress, GRAIN, PointResult, ROUND, Setup, StaticResult, Stats, budget, factor_error, factor_raw, ignored_dashpots, pattern, reduction,
    split,
};
use crate::cholesky::FactorError;
use crate::constraints::Reduction;
use crate::contact::{self, ContactResult, Interface};
use crate::ldlt::{Kind, Symmetric};
use crate::sparse::{Csr, NodeGraph, add_block, add_entry};
use crate::system::{Mesh, add_element_thermal_loads, loads, temperature_loads, translation_block, translations};
use crate::{Error, Limits, Result};

#[cfg(test)]
mod band_tests;

/// One Newton iteration of an increment.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, serde::Deserialize)]
pub struct Iteration {
    /// Largest residual force over the force scale, or largest residual moment over the moment
    /// scale, whichever is larger.
    pub force: f64,
    /// Largest translation of the correction over the increment's translation reference, or
    /// largest rotation over its rotation reference, whichever is larger.
    pub displacement: f64,
    /// Work of the correction against the residual over the first iteration's.
    pub energy: f64,
    /// The line search's step along the correction (1: the full Newton step).
    pub line_search: f64,
    /// Negative eigenvalues of the tangent the correction was solved with (an indefinite
    /// tangent is factored by LDLᵀ, [`crate::ldlt`]).
    #[serde(default, skip_serializing_if = "is_zero")]
    pub negative_pivots: usize,
    /// Contact: the augmentation it belongs to (0: the increment's first solution).
    #[serde(default, skip_serializing_if = "is_zero")]
    pub augmentation: usize,
}

fn is_zero(n: &usize) -> bool {
    *n == 0
}

/// One increment (converged, or cut back).
#[derive(Debug, Clone, PartialEq, Serialize, serde::Deserialize)]
pub struct Increment {
    /// From 1, attempts included.
    pub number: usize,
    /// The load factor reached (when converged) or attempted.
    pub load_factor: f64,
    /// Its size, as a fraction of the step.
    pub size: f64,
    pub iterations: Vec<Iteration>,
    pub converged: bool,
    /// Why it was cut back.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cutback: Option<String>,
    /// Negative eigenvalues of the tangent at the converged state (arc-length steps, which
    /// factor it for the next predictor).
    #[serde(default, skip_serializing_if = "is_zero")]
    pub negative_pivots: usize,
    /// Contact: the times the increment was solved again with updated multipliers.
    #[serde(default, skip_serializing_if = "is_zero")]
    pub augmentations: usize,
}

/// Where the tangent's inertia changes along an arc-length path.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CriticalKind {
    /// The load factor turns (snap-through).
    Limit,
    /// The load factor goes on: another branch crosses this one.
    Bifurcation,
}

/// A critical point of an arc-length path: the tangent's count of negative eigenvalues changes
/// between two converged states.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, serde::Deserialize)]
pub struct CriticalPoint {
    pub kind: CriticalKind,
    /// The increment whose converged state has the new count.
    pub increment: usize,
    /// For a limit point the extreme load factor of the converged states around it, for a
    /// bifurcation the last one before it.
    pub load_factor: f64,
    /// Negative eigenvalues before and after.
    pub negative_before: usize,
    pub negative_after: usize,
}

/// The displacements and reactions at the end of a converged increment.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Frame {
    pub load_factor: f64,
    /// Per node, in [`StaticResult::nodes`] order.
    pub displacement: Vec<[f64; 3]>,
    /// The reaction forces at the nodes that carry one, as (index in [`StaticResult::nodes`],
    /// force), by index; nodes left out carry none. Sparse: only supported or prescribed nodes
    /// appear, so a frame costs little beyond its displacements.
    pub reaction: Vec<(usize, [f64; 3])>,
}

impl Frame {
    /// The reaction force at the node of index `i` (zero when it carries none).
    pub fn reaction_at(&self, i: usize) -> [f64; 3] {
        self.reaction.binary_search_by_key(&i, |(k, _)| *k).ok().and_then(|k| self.reaction.get(k)).map_or([0.0; 3], |(_, f)| *f)
    }
}

/// The increments of a nonlinear static step and the displacements after each converged one.
#[derive(Debug, Clone, PartialEq, Default, serde::Serialize, serde::Deserialize)]
pub struct History {
    pub increments: Vec<Increment>,
    /// The undeformed state (load factor 0), then one per converged increment while they fit
    /// [`MAX_FRAME_VALUES`].
    pub frames: Vec<Frame>,
    /// Newton iterations over the step.
    pub iterations: usize,
    /// Arc-length steps: where the tangent's inertia changes along the path.
    pub critical_points: Vec<CriticalPoint>,
}

impl History {
    /// The load factor of every frame.
    pub fn load_factors(&self) -> Vec<f64> {
        self.frames.iter().map(|f| f.load_factor).collect()
    }
}

/// Most values the frames of one step keep (3 per node per frame for the displacements, 4 per
/// reacted node for the reactions).
pub const MAX_FRAME_VALUES: usize = 30_000_000;

/// Line search: accept a step whose residual work is at most this fraction of the first.
const LINE_SEARCH_ACCEPT: f64 = 0.8;
/// Line search: most trial steps, and the range of the step.
const LINE_SEARCH_TRIALS: usize = 5;
const LINE_SEARCH_MIN: f64 = 0.1;
/// Cutback factor of a failed increment, growth after fast ones and shrinking after slow ones.
const CUTBACK: f64 = 0.25;
const GROWTH: f64 = 1.5;
const SHRINK: f64 = 0.75;
const FAST_ITERATIONS: usize = 4;
const SLOW_ITERATIONS: usize = 8;
/// A force residual this fraction of its tolerance is round-off.
const ROUND_OFF: f64 = 1e-3; // unit-factors: relative to the force tolerance
/// Increments and times within this fraction of each other are the same: what is left of summing
/// increments.
const SAME_TIME: f64 = 1e-9; // unit-factors: relative to the increment
/// The smallest displacement scale, as a fraction of the model's size.
const LENGTH_FLOOR: f64 = 1e-12; // unit-factors: relative to the model's size
/// A follower load's skew stiffness above this fraction of the whole is worth a note.
const NOTABLE_ASYMMETRY: f64 = 1e-6; // unit-factors: a ratio of norms

/// A pressure on element faces (or 2D edges) that follows the deformation.
struct Follower {
    element: u32,
    /// Layout slots of the face's nodes, and their reference positions.
    nodes: Vec<usize>,
    xyz: Vec<[f64; 3]>,
    magnitude: f64,
    /// The out-of-plane extent of a 2D edge; `None` for a solid face.
    extent: Option<Extent>,
    /// The translations of its nodes (x, y, z, or x, y per node).
    dofs: Vec<usize>,
}

/// Everything an evaluation of the out-of-balance forces needs.
struct Problem<'a> {
    mesh: &'a Mesh,
    red: &'a Reduction,
    /// Dead loads at λ = 1.
    dead: Vec<f64>,
    followers: Vec<Follower>,
    /// Springs (linear along their fixed directions), on the pattern.
    springs: Option<Csr>,
    /// Large-deformation kinematics.
    finite: bool,
    /// Each element's global degrees of freedom, in its matrices' order, and the kinds per node.
    element_dofs: Vec<Vec<usize>>,
    element_kinds: Vec<DofSet>,
    /// Each element's material model.
    materials: Vec<MaterialModel>,
    /// The step's temperature field at λ = 1, per mesh node (`None` where the field leaves a
    /// node at its material's reference temperature); `None` without one.
    field: Option<Vec<Option<f64>>>,
    /// The tangent's pattern, zero, and the mesh's node graph it comes from (contact adds the
    /// blocks of its face pairs to it).
    pattern: Csr,
    graph: NodeGraph,
    /// The contact pairs ([`crate::contact`]).
    contacts: Vec<Interface>,
}

/// The internal and external forces at a state, and the tangent when asked for.
struct Evaluation {
    internal: Vec<f64>,
    external: Vec<f64>,
    tangent: Option<Csr>,
    /// Every element's material history at this state (committed when the increment
    /// converges).
    history: Vec<Vec<f64>>,
    /// Each contact's state.
    contact: Vec<contact::State>,
}

/// Why an attempt at an increment failed, when that calls for a cutback rather than an error.
enum Failure {
    Cutback(String),
    Fatal(Error),
}

impl From<Error> for Failure {
    fn from(e: Error) -> Failure {
        match e {
            // An inverted element or a material that cannot take the deformation: a smaller
            // increment may avoid it.
            Error::Element { label, source: fieldcraft_fem::Error::Degenerate(m) } => Failure::Cutback(format!("element {label}: {m}")),
            Error::NotConverged(m) => Failure::Cutback(m),
            e => Failure::Fatal(e),
        }
    }
}

/// One hash of the contacts' active sets.
fn active_sets(e: &Evaluation) -> u64 {
    e.contact.iter().fold(0u64, |h, s| h.rotate_left(17) ^ s.active_set())
}

fn norm_inf(v: &[f64]) -> f64 {
    v.iter().fold(0.0, |m, x| m.max(x.abs()))
}

/// The largest magnitudes of `v`'s translational and rotational entries (`rotation` marks the
/// rotations; entries past its end are translations): forces and moments of a force vector,
/// translations and rotations of a motion.
fn norm_inf_split(v: &[f64], rotation: &[bool]) -> (f64, f64) {
    v.iter().enumerate().fold(
        (0.0f64, 0.0f64),
        |(t, r), (i, x)| {
            if rotation.get(i).copied().unwrap_or(false) { (t, r.max(x.abs())) } else { (t.max(x.abs()), r) }
        },
    )
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// ½(K + sign Kᵀ) of a row-major m × m matrix: its symmetric (sign 1) or skew (−1) part.
fn half(k: &[f64], m: usize, sign: f64) -> Vec<f64> {
    (0..m * m).map(|i| 0.5 * (k.get(i).copied().unwrap_or(0.0) + sign * k.get((i % m) * m + i / m).copied().unwrap_or(0.0))).collect()
}

/// T_fᵀ K T_p (s u_p): the free equations' forces of the prescribed motions scaled by `s`.
fn prescribed_force(red: &Reduction, k: &Csr, s: f64) -> Vec<f64> {
    let jump = red.expand_scaled(&vec![0.0; red.neq()], s);
    let mut kj = vec![0.0; jump.len()];
    k.mul(&jump, &mut kj);
    red.restrict(&kj)
}

/// The out-of-balance forces over the free equations, T_fᵀ (f_ext − f_int).
fn free_residual(red: &Reduction, e: &Evaluation) -> Vec<f64> {
    let r: Vec<f64> = e.external.iter().zip(&e.internal).map(|(a, b)| a - b).collect();
    red.restrict(&r)
}

impl Problem<'_> {
    /// Nodal displacements of an element (3 per node for solids, 2 for 2D elements).
    fn element_vector(&self, e: usize, u: &[f64]) -> Vec<f64> {
        self.element_dofs.get(e).map(|dofs| dofs.iter().map(|d| u.get(*d).copied().unwrap_or(0.0)).collect()).unwrap_or_default()
    }

    /// How element `e` is evaluated: total Lagrangian with nonlinear geometry, with small
    /// strain when its material is inelastic, and by its linear stiffness otherwise (`None`).
    fn kinematics(&self, e: usize) -> Option<Kinematics> {
        if self.finite { Some(Kinematics::Finite) } else { self.materials.get(e).filter(|m| m.is_path_dependent()).map(|_| Kinematics::Small) }
    }

    /// The temperatures of element `e` at load factor λ, which ramps them from the reference
    /// temperature to the field's (J2 Johnson–Cook softening and thermal strain): the
    /// material's expansion and reference temperature, or no expansion from the Johnson–Cook
    /// room temperature; `None` without a field or for a material neither affects.
    fn temperature(&self, e: usize, lambda: f64) -> Option<Temperature> {
        let field = self.field.as_ref()?;
        let el = self.mesh.elements.get(e)?;
        let softens = el.plastic.as_ref().is_some_and(|j2| j2.is_temperature_dependent());
        let (alpha, reference) = match el.expansion {
            Some(x) => x,
            None if softens => ([0.0; 3], el.plastic.as_ref().map_or(0.0, |j2| j2.initial_temperature())),
            None => return None,
        };
        let nodal: Vec<f64> =
            el.nodes.iter().map(|n| field.get(*n).copied().flatten().map_or(reference, |t| reference + lambda * (t - reference))).collect();
        (softens || nodal.iter().any(|t| *t != reference)).then_some(Temperature { alpha, reference, nodal })
    }

    /// Internal forces and the tangent of element `e` at displacements `u` (all degrees of
    /// freedom) and load factor λ.
    fn element(&self, e: usize, u: &[f64], lambda: f64, tangent: bool, committed: &[f64]) -> Result<(Vec<f64>, Option<Matrix>, Vec<f64>)> {
        let el = self.mesh.elements.get(e).ok_or_else(|| Error::Internal("element index".into()))?;
        let at = |err| Error::element(el.label, err);
        let ue = self.element_vector(e, u);
        let kernel = el.kernel()?;
        let Some(kin) = self.kinematics(e) else {
            let k = kernel.stiffness(&el.d).map_err(at)?;
            let f = k.mul_vec(&ue);
            return Ok((f, tangent.then_some(k), committed.to_vec()));
        };
        let material = self.materials.get(e).ok_or_else(|| Error::Internal("element material".into()))?;
        let mut updated = committed.to_vec();
        let r = kernel.nonlinear(material, &ue, kin, self.temperature(e, lambda).as_ref(), committed, &mut updated, tangent).map_err(at)?;
        Ok((r.force, r.tangent, updated))
    }

    /// The forces of a follower pressure at λ (on [`Follower::dofs`]) and its load stiffness
    /// (full, not symmetric).
    fn follower(&self, f: &Follower, u: &[f64], lambda: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        let at = |err| Error::element(f.element, err);
        let ue: Vec<f64> = f.dofs.iter().map(|d| u.get(*d).copied().unwrap_or(0.0)).collect();
        let p = lambda * f.magnitude;
        let (force, k) = match f.extent {
            None => {
                let (force, k) = face::follower_pressure(&f.xyz, ue.as_chunks::<3>().0, p).map_err(at)?;
                (force.into_iter().flatten().collect(), k)
            }
            Some(extent) => {
                let (force, k) = face::follower_edge_pressure(&f.xyz, ue.as_chunks::<2>().0, extent, p).map_err(at)?;
                (force.into_iter().flatten().collect(), k)
            }
        };
        Ok((force, k))
    }

    /// The forces at displacements `u` and load factor λ, and the tangent (symmetric part of the
    /// load stiffness included) when `tangent`.
    /// `contact` is each contact's state at the last converged increment (augmented Lagrange
    /// multipliers, friction tractions; none: zero) and `reference` the displacements there (the
    /// start of the friction's slip increment; none: zero).
    fn evaluate(
        &self,
        u: &[f64],
        lambda: f64,
        tangent: bool,
        committed: &[Vec<f64>],
        contact: &[contact::Committed],
        reference: &[f64],
        band: f64,
    ) -> Result<Evaluation> {
        let layout = &self.mesh.layout;
        let n = u.len();
        let mut internal = vec![0.0; n];
        // Contact first: its face pairs extend the tangent's pattern.
        let mut blocks = contact::Blocks::new();
        let mut states = Vec::with_capacity(self.contacts.len());
        for (i, c) in self.contacts.iter().enumerate() {
            let empty = contact::Committed::default();
            states.push(c.evaluate(self.mesh, u, reference, contact.get(i).unwrap_or(&empty), band, tangent, &mut internal, &mut blocks)?);
        }
        let mut k = tangent.then(|| {
            if blocks.is_empty() {
                self.pattern.clone()
            } else {
                self.graph.with_blocks(blocks.iter().map(|(nodes, _)| nodes.as_slice())).zero_matrix(layout)
            }
        });
        if let Some(k) = k.as_mut() {
            // Contact blocks are over X, Y and Z per node; a 2D model keeps the in-plane rows.
            let kinds = self.mesh.translation_kinds();
            for (nodes, data) in &blocks {
                add_block(k, layout, nodes, kinds, &translation_block(data, nodes.len(), kinds))?;
            }
        }
        let mut history = Vec::with_capacity(self.mesh.elements.len());
        for round in parallel::blocks(self.mesh.elements.len(), ROUND) {
            let start = round.start;
            let computed = parallel::map_grained(round.len(), GRAIN, |i| {
                self.element(start + i, u, lambda, tangent, committed.get(start + i).map(Vec::as_slice).unwrap_or_default())
            });
            for (i, c) in computed.into_iter().enumerate() {
                let (f, ke, updated) = c?;
                history.push(updated);
                let e = start + i;
                let el = self.mesh.elements.get(e).ok_or_else(|| Error::Internal("element index".into()))?;
                for (d, v) in self.element_dofs.get(e).into_iter().flatten().zip(&f) {
                    if let Some(x) = internal.get_mut(*d) {
                        *x += v;
                    }
                }
                if let (Some(k), Some(ke)) = (k.as_mut(), ke) {
                    let kinds = self.element_kinds.get(e).copied().ok_or_else(|| Error::Internal("element kinds".into()))?;
                    add_block(k, layout, &el.nodes, kinds, &ke.data)?;
                }
            }
        }
        if let Some(s) = &self.springs {
            let mut ku = vec![0.0; s.n];
            s.mul(u.get(..s.n).unwrap_or_default(), &mut ku);
            for (x, v) in internal.iter_mut().zip(&ku) {
                *x += v;
            }
            // On the mesh's pattern, which contact may have extended.
            if let Some(k) = k.as_mut() {
                if k.ptr == s.ptr && k.idx == s.idx {
                    k.val.iter_mut().zip(&s.val).for_each(|(x, v)| *x += v);
                } else {
                    for r in 0..s.n {
                        let (cols, vals) = s.row(r);
                        for (c, v) in cols.iter().zip(vals).filter(|(_, v)| **v != 0.0) {
                            add_entry(k, r, *c, *v)?;
                        }
                    }
                }
            }
        }
        let mut external: Vec<f64> = self.dead.iter().map(|f| lambda * f).collect();
        for f in &self.followers {
            let (force, kl) = self.follower(f, u, lambda)?;
            for (d, v) in f.dofs.iter().zip(&force) {
                if let Some(x) = external.get_mut(*d) {
                    *x += v;
                }
            }
            if let Some(k) = k.as_mut() {
                add_block(k, layout, &f.nodes, self.mesh.translation_kinds(), &half(&kl, f.dofs.len(), 1.0))?;
            }
        }
        if internal.iter().chain(&external).any(|v| !v.is_finite()) {
            return Err(Error::NotConverged("the forces are not finite".into()));
        }
        Ok(Evaluation { internal, external, tangent: k, history, contact: states })
    }

    /// The loads at load factor 1 at displacements `u` (dead loads and follower pressures on
    /// the current faces): ∂f_ext/∂λ.
    fn reference_load(&self, u: &[f64]) -> Result<Vec<f64>> {
        let mut external = self.dead.clone();
        for f in &self.followers {
            let (force, _) = self.follower(f, u, 1.0)?;
            for (d, v) in f.dofs.iter().zip(&force) {
                if let Some(x) = external.get_mut(*d) {
                    *x += v;
                }
            }
        }
        Ok(external)
    }

    /// The size of the follower load stiffness's skew part relative to the whole, over the free
    /// degrees of freedom (Frobenius norms), at `u` and λ; `None` without followers.
    fn follower_asymmetry(&self, u: &[f64], lambda: f64) -> Result<Option<f64>> {
        if self.followers.is_empty() {
            return Ok(None);
        }
        let layout = &self.mesh.layout;
        let (mut skew, mut whole) = (self.pattern.clone(), self.pattern.clone());
        for f in &self.followers {
            let (_, kl) = self.follower(f, u, lambda)?;
            add_block(&mut skew, layout, &f.nodes, self.mesh.translation_kinds(), &half(&kl, f.dofs.len(), -1.0))?;
            add_block(&mut whole, layout, &f.nodes, self.mesh.translation_kinds(), &kl)?;
        }
        let mut free = vec![false; skew.n];
        for d in &self.red.eq_dof {
            if let Some(x) = free.get_mut(*d) {
                *x = true;
            }
        }
        let norm = |a: &Csr| -> f64 {
            let mut sum = 0.0;
            for r in (0..a.n).filter(|r| free.get(*r).copied().unwrap_or(false)) {
                let (cols, vals) = a.row(r);
                for (c, v) in cols.iter().zip(vals) {
                    if free.get(*c).copied().unwrap_or(false) {
                        sum += v * v;
                    }
                }
            }
            sum.sqrt()
        };
        let total = norm(&whole);
        Ok((total > 0.0).then(|| norm(&skew) / total))
    }
}

/// Mutable state of the Newton loop over the step.
struct Run<'a> {
    problem: &'a Problem<'a>,
    settings: &'a NonlinearStatic,
    limits: &'a Limits,
    /// The largest force seen so far (∞-norm over the translations), and the round-off floor of
    /// the force scale.
    force_scale: f64,
    force_floor: Option<f64>,
    /// The largest moment seen so far (∞-norm over the rotations); the moment scale is at least
    /// the force scale times [`Run::length`] ([`Run::moment_scale`]).
    moment_scale: f64,
    /// The model's size (the diagonal of its bounding box): the characteristic length that
    /// relates moments to forces and rotations to translations.
    length: f64,
    /// Which degrees of freedom are rotations (reference points' today; their forces are
    /// moments), over every degree of freedom and over the free equations: the convergence
    /// tests measure them apart from the translations.
    rotation: Vec<bool>,
    rotation_eq: Vec<bool>,
    /// The largest prescribed translation and rotation at load factor 1.
    prescribed: (f64, f64),
    iterations: usize,
    factor_entries: usize,
    /// Every element's material history at the end of the last converged increment
    /// (integration points × the material's history length).
    committed: Vec<Vec<f64>>,
    /// Whether a material keeps history (the tangent predictor and its convergence rules).
    path_dependent: bool,
    /// How the undeformed model's tangent factors (Cholesky), once asked.
    undeformed: Option<std::result::Result<(), FactorError>>,
    /// The contacts' multipliers per secondary node (the last converged increment's, or the
    /// current augmentation's), and the augmentations over the step.
    contact: Vec<contact::Committed>,
    /// The displacements of the last converged increment (all degrees of freedom).
    reference: Vec<f64>,
    /// The band of the engaged contact nodes ([`contact::ENGAGED_BAND`] once an increment's
    /// active set has chattered, through that increment's retries and augmentations; else 0).
    band: f64,
    augmentations: usize,
    /// The augmentation being solved in the current increment.
    augmentation: usize,
    /// The largest force and moment of the converged states (contact's scales,
    /// [`Run::residual`]).
    converged_scale: f64,
    converged_moment: f64,
}

/// The outcome of a converged increment.
struct Converged {
    q: Vec<f64>,
    evaluation: Evaluation,
    /// The work of the increment's first correction (the energy criterion's reference).
    energy: Option<f64>,
}

impl Run<'_> {
    /// Commit a fully converged load-control increment, after its contact augmentations.
    fn commit_increment(&mut self, q: &[f64], lambda: f64, evaluation: &Evaluation) {
        self.committed.clone_from(&evaluation.history);
        for (l, s) in self.contact.iter_mut().zip(&evaluation.contact) {
            l.commit(s);
        }
        self.reference = self.problem.red.expand_scaled(q, lambda);
        let (external, external_moment) = norm_inf_split(&evaluation.external, &self.rotation);
        let (internal, internal_moment) = norm_inf_split(&evaluation.internal, &self.rotation);
        self.converged_scale = self.converged_scale.max(external).max(internal);
        self.converged_moment = self.converged_moment.max(external_moment).max(internal_moment);
        // A near-open node has no force but keeps a penalty tangent while the band is active.
        // That temporary anti-chatter aid must not poison a later increment's predictor/Newton
        // tangent. Only clear it after this increment and every augmentation have committed.
        self.band = 0.0;
    }

    /// The residual over the free equations, T_fᵀ (f_ext − f_int), and the force scale.
    ///
    /// With contact a trial state can be far from equilibrium (a body pressed deep into another
    /// before the active set settles, its contact forces huge), so the scale is not the largest
    /// force of every trial but that of the converged states, the loads and the current state.
    fn residual(&mut self, e: &Evaluation) -> Vec<f64> {
        let (external, external_moment) = norm_inf_split(&e.external, &self.rotation);
        let (internal, internal_moment) = norm_inf_split(&e.internal, &self.rotation);
        let (current, moment) = (external.max(internal), external_moment.max(internal_moment));
        if self.problem.contacts.is_empty() {
            self.force_scale = self.force_scale.max(current);
            self.moment_scale = self.moment_scale.max(moment);
        } else {
            self.force_scale = self.converged_scale.max(current);
            self.moment_scale = self.converged_moment.max(moment);
        }
        free_residual(self.problem.red, e)
    }

    /// The round-off floor of the force scale, from the first tangent's mean diagonal: a
    /// translation's stiffness as it is, a rotation's over the length squared (the stiffness of
    /// the translation a rotation gives at the model's size).
    fn seed_force_floor(&mut self, kff: &Csr) {
        if self.force_floor.is_none() {
            let diagonal = kff.diagonal();
            let lever = self.length * self.length;
            let sum: f64 = diagonal
                .iter()
                .enumerate()
                .map(|(i, d)| if self.rotation_eq.get(i).copied().unwrap_or(false) { d.abs() / lever } else { d.abs() })
                .sum();
            let mean = sum / diagonal.len().max(1) as f64;
            self.force_floor = Some(1e-8 * mean * self.length);
        }
    }

    /// The force scale: the largest force seen, at least the round-off floor.
    fn scale(&self) -> f64 {
        self.force_scale.max(self.force_floor.unwrap_or(0.0)).max(f64::MIN_POSITIVE)
    }

    /// The moment scale: the largest moment seen, at least the force scale times the model's
    /// size (the moment of the largest force about a point of the model).
    fn moment_scale(&self) -> f64 {
        self.moment_scale.max(self.scale() * self.length).max(f64::MIN_POSITIVE)
    }

    /// The force criterion's measure of a residual over the free equations: its largest force
    /// over the force scale or its largest moment over the moment scale, whichever is larger.
    fn force_measure(&self, r: &[f64]) -> f64 {
        let (force, moment) = norm_inf_split(r, &self.rotation_eq);
        let measure = force / self.scale();
        if moment > 0.0 { measure.max(moment / self.moment_scale()) } else { measure }
    }

    /// Whether a residual's forces are within `fraction` of the force scale and its moments
    /// within `fraction` of the moment scale.
    fn force_within(&self, r: &[f64], fraction: f64) -> bool {
        let (force, moment) = norm_inf_split(r, &self.rotation_eq);
        force <= fraction * self.scale() && moment <= fraction * self.moment_scale()
    }

    /// The displacement criterion's references from the largest translation and rotation of
    /// the increment so far: a translation's, at least 10⁻¹² of the model's size; a rotation's,
    /// at least the rotation that moves a point at the model's size by the translation's
    /// reference (so a rotation the increment barely changes is measured by what its error
    /// moves, not against round-off).
    fn motion_references(&self, translation: f64, rotation: f64) -> (f64, f64) {
        let t = translation.max(LENGTH_FLOOR * self.length);
        (t, rotation.max(t / self.length))
    }

    /// The displacement criterion's measure of the correction `s · dq`: its largest translation
    /// over the translation reference or its largest rotation over the rotation reference,
    /// whichever is larger.
    fn correction_measure(&self, dq: &[f64], s: f64, (translation, rotation): (f64, f64)) -> f64 {
        let (t, r) = norm_inf_split(dq, &self.rotation_eq);
        let measure = s * t / translation;
        if r > 0.0 { measure.max(s * r / rotation) } else { measure }
    }

    fn state(&self, q: &[f64], lambda: f64, tangent: bool) -> std::result::Result<Evaluation, Failure> {
        Ok(self.problem.evaluate(
            &self.problem.red.expand_scaled(q, lambda),
            lambda,
            tangent,
            &self.committed,
            &self.contact,
            &self.reference,
            self.band,
        )?)
    }

    /// How the tangent of the undeformed, unloaded model factors (Cholesky): a first tangent
    /// that is not positive definite, at a state the first increment has already strained, is
    /// then the material or the geometry, not a model that can move freely (this error).
    fn undeformed(&mut self) -> std::result::Result<std::result::Result<(), FactorError>, Failure> {
        if let Some(u) = &self.undeformed {
            return Ok(u.clone());
        }
        let mut e = self.problem.evaluate(&vec![0.0; self.problem.mesh.dofs()], 0.0, true, &self.committed, &[], &[], 0.0)?;
        let k = e.tangent.take().ok_or_else(|| Failure::Fatal(Error::Internal("tangent".into())))?;
        let [kff] = self.problem.red.reduce([&k])?;
        let u = factor_raw(&kff, self.limits).map(|_| ());
        self.undeformed = Some(u.clone());
        Ok(u)
    }

    /// Factor a tangent: Cholesky, or LDLᵀ when it is not positive definite (an indefinite
    /// tangent: a strained mixed element, an unstable material, a state past a limit or
    /// bifurcation point; Newton proceeds on it and the result counts its negative pivots). A
    /// model that can move freely is an error at its first tangent; a singular tangent later
    /// is a cutback.
    fn factor(&mut self, kff: &Csr, first: bool) -> std::result::Result<Symmetric, Failure> {
        let free_motion = |run: &mut Self| -> std::result::Result<Option<Error>, Failure> {
            if !first {
                return Ok(None);
            }
            Ok(match run.undeformed()? {
                Err(e @ (FactorError::NotPositiveDefinite { .. } | FactorError::Singular { .. })) => {
                    Some(factor_error(e, run.problem.mesh, run.problem.red))
                }
                _ => None,
            })
        };
        match factor_symmetric(kff, self.limits, Kind::Auto) {
            Ok(f @ Symmetric::Cholesky(_)) => Ok(f),
            Ok(f) => match free_motion(self)? {
                Some(e) => Err(Failure::Fatal(e)),
                None => Ok(f),
            },
            Err(FactorError::NotPositiveDefinite { .. } | FactorError::Singular { .. }) => match free_motion(self)? {
                Some(e) => Err(Failure::Fatal(e)),
                None => Err(Failure::Cutback("the tangent stiffness is singular (a limit point or bifurcation, or an unstable material)".into())),
            },
            Err(e) => Err(Failure::Fatal(factor_error(e, self.problem.mesh, self.problem.red))),
        }
    }

    /// The tangent predictor of an increment from the converged `q0` at load factor `from` to
    /// `to`, and its work δq · r (the increment's reference for the energy criterion, since the
    /// Newton iterations after it only correct it): one linear solve with the tangent of the converged state for the change of the
    /// loads and of the prescribed displacements, K δq = T_fᵀ (f_ext(to) − f_int(u₀) − K T_p Δλ u_p).
    /// It spreads a prescribed motion into the interior before the material sees it, which a
    /// path-dependent material needs: started from the boundary's jump alone, the elements next
    /// to it would flow plastically far beyond the solution. `None` when the tangent cannot be
    /// factored (the increment then starts from `q0`).
    fn tangent_predictor(&mut self, q0: &[f64], from: f64, to: f64) -> std::result::Result<Option<(Vec<f64>, f64)>, Failure> {
        let red = self.problem.red;
        let mut e = self.problem.evaluate(&red.expand_scaled(q0, from), to, true, &self.committed, &self.contact, &self.reference, self.band)?;
        let k = e.tangent.take().ok_or_else(|| Failure::Fatal(Error::Internal("tangent".into())))?;
        let rf: Vec<f64> = free_residual(red, &e).iter().zip(prescribed_force(red, &k, to - from)).map(|(a, b)| a - b).collect();
        let mut dq = rf.clone();
        let [kff] = red.reduce([&k])?;
        drop(k);
        let chol = match factor_raw(&kff, self.limits) {
            Ok(chol) => chol,
            Err(FactorError::Cancelled) => return Err(Failure::Fatal(fieldcraft_model::Error::Cancelled.into())),
            Err(_) => return Ok(None),
        };
        chol.solve(&mut dq)?;
        let energy = dot(&dq, &rf).abs();
        Ok(dq.iter().all(|v| v.is_finite()).then(|| (q0.iter().zip(&dq).map(|(a, d)| a + d).collect(), energy)))
    }

    /// Newton iterations from the converged `q0` (at the previous load factor) to load factor
    /// `lambda`.
    fn increment(
        &mut self,
        q0: &[f64],
        predictor: Option<&[f64]>,
        predicted_energy: Option<f64>,
        lambda: f64,
        size: f64,
        iterations: &mut Vec<Iteration>,
    ) -> std::result::Result<Converged, Failure> {
        let red = self.problem.red;
        let tol = self.settings.tolerances;
        // Start from the extrapolated predictor when it is a valid state, else from q0.
        let start = predictor.map(|p| self.state(p, lambda, true).map(|e| (p.to_vec(), e)));
        let (mut q, mut eval) = match start {
            Some(Ok(s)) => s,
            Some(Err(Failure::Fatal(e))) => return Err(Failure::Fatal(e)),
            Some(Err(Failure::Cutback(_))) | None => (q0.to_vec(), self.state(q0, lambda, true)?),
        };
        let mut r = self.residual(&eval);
        // With path-dependent materials, an increment that starts within the force tolerance (a
        // predictor that already solved it) has only round-off left for the energy criterion to
        // measure against: it is waived.
        let started_converged = self.path_dependent && self.force_floor.is_some() && self.force_within(&r, tol.force);
        let prescribed = (size * self.prescribed.0, size * self.prescribed.1);
        // The energy criterion's reference: the first correction's work, or the tangent
        // predictor's when the increment started from one that is still in use.
        let mut first_energy: Option<f64> = predicted_energy.filter(|_| q.as_slice() != q0);
        let mut history: Vec<f64> = Vec::new();
        // Contact: the active sets of the iterations, and how often they went back to the one
        // before the last (chattering).
        let mut sets = (None, active_sets(&eval));
        let mut flips = 0usize;
        let max_iterations = usize::try_from(self.settings.max_iterations).unwrap_or(usize::MAX);
        for _ in 0..max_iterations {
            fieldcraft_model::progress::check().map_err(Error::from)?;
            // The round-off rule: nothing left to correct.
            if self.force_within(&r, ROUND_OFF * tol.force) && self.force_floor.is_some() {
                return Ok(Converged { q, evaluation: eval, energy: first_energy });
            }
            if self.iterations >= self.limits.max_newton_iterations {
                return Err(Failure::Fatal(Error::Limit(format!(
                    "the step took {} Newton iterations, the limit (raise max_newton_iterations or loosen the tolerances)",
                    self.iterations
                ))));
            }
            self.iterations += 1;
            let k = eval.tangent.take().ok_or_else(|| Failure::Fatal(Error::Internal("tangent".into())))?;
            let [kff] = red.reduce([&k])?;
            // The first tangent of the step: a singular one is a modelling error, not a cutback.
            let first = self.force_floor.is_none();
            if first {
                self.seed_force_floor(&kff);
            }
            drop(k);
            let factor = self.factor(&kff, first)?;
            self.factor_entries = factor.entries();
            drop(kff);
            let mut dq = r.clone();
            factor.solve(&mut dq)?;
            if dq.iter().any(|v| !v.is_finite()) {
                return Err(Failure::Cutback("the Newton correction is not finite".into()));
            }
            let g0 = dot(&dq, &r);
            let (s, next) = self.line_search(&q, &dq, lambda, g0)?;
            for (x, d) in q.iter_mut().zip(&dq) {
                *x += s * d;
            }
            eval = next;
            if !self.problem.contacts.is_empty() {
                let set = active_sets(&eval);
                if sets.0 == Some(set) && set != sets.1 {
                    flips += 1;
                    if flips >= contact::CHATTER_LIMIT {
                        // Through this increment's retries/augmentations the nodes about to
                        // touch stay in the tangent; the accepted commit clears the band.
                        self.band = contact::ENGAGED_BAND;
                        return Err(Failure::Cutback("the contact status alternates between two states (chattering)".into()));
                    }
                }
                sets = (Some(sets.1), set);
            }
            let r_old = std::mem::replace(&mut r, self.residual(&eval));
            let energy = (s * dot(&dq, &r_old)).abs();
            let first = *first_energy.get_or_insert(energy);
            let moved: Vec<f64> = q.iter().zip(q0).map(|(a, b)| a - b).collect();
            let (translated, rotated) = norm_inf_split(&moved, &self.rotation_eq);
            let reference = self.motion_references(translated.max(prescribed.0), rotated.max(prescribed.1));
            let it = Iteration {
                force: self.force_measure(&r),
                displacement: self.correction_measure(&dq, s, reference),
                energy: if first > 0.0 { energy / first } else { 0.0 },
                line_search: s,
                negative_pivots: factor.inertia().negative,
                augmentation: self.augmentation,
            };
            iterations.push(it);
            if !it.force.is_finite() {
                return Err(Failure::Cutback("the residual is not finite".into()));
            }
            let energy_met = it.energy <= tol.energy || started_converged;
            if (it.force <= tol.force && it.displacement <= tol.displacement && energy_met) || it.force <= ROUND_OFF * tol.force {
                return Ok(Converged { q, evaluation: eval, energy: first_energy });
            }
            history.push(it.force);
            if let [.., a, b, c] = history[..]
                && c > b
                && b > a
            {
                return Err(Failure::Cutback(format!("the residual grows ({a:.2e}, {b:.2e}, {c:.2e} of the force scale)")));
            }
        }
        Err(Failure::Cutback(format!("no convergence in {} iterations", self.settings.max_iterations)))
    }

    /// The augmented Lagrangian loop of a converged increment (Uzawa: the contacts' multipliers
    /// take the converged pressures, and the increment is solved again from its solution) until
    /// the gaps are within [`contact::GAP_TOLERANCE`] ([`Interface::augmentation_error`]), at
    /// most [`contact::MAX_AUGMENTATIONS`] times (then with a warning). Without contact it
    /// returns the increment as it is.
    fn augment(
        &mut self,
        mut c: Converged,
        q0: &[f64],
        lambda: f64,
        size: f64,
        iterations: &mut Vec<Iteration>,
        count: &mut usize,
        warnings: &mut Vec<String>,
    ) -> std::result::Result<Converged, Failure> {
        if self.problem.contacts.is_empty() {
            return Ok(c);
        }
        loop {
            let error = self
                .problem
                .contacts
                .iter()
                .zip(&self.contact)
                .zip(&c.evaluation.contact)
                .fold(0.0f64, |m, ((i, l), s)| m.max(i.augmentation_error(&l.multipliers, s)));
            if error <= contact::GAP_TOLERANCE {
                return Ok(c);
            }
            if *count >= contact::MAX_AUGMENTATIONS {
                warnings.push(format!(
                    "contact: at load factor {lambda:.6} a gap is still {error:.1e} of its element's thickness after {} augmentations (the penalty's regularisation)",
                    contact::MAX_AUGMENTATIONS
                ));
                return Ok(c);
            }
            *count += 1;
            self.augmentations += 1;
            self.augmentation = *count;
            for (l, s) in self.contact.iter_mut().zip(&c.evaluation.contact) {
                l.multipliers.clone_from(&s.pressure);
            }
            // From the converged state, measured against the increment's start and its first
            // correction like the increment itself.
            let start = std::mem::take(&mut c.q);
            c = self.increment(q0, Some(&start), c.energy, lambda, size, iterations)?;
        }
    }

    /// The line search along `dq` from `q`: the accepted step and the state there (with the
    /// tangent).
    fn line_search(&mut self, q: &[f64], dq: &[f64], lambda: f64, g0: f64) -> std::result::Result<(f64, Evaluation), Failure> {
        let at = |s: f64| -> Vec<f64> { q.iter().zip(dq).map(|(a, b)| a + s * b).collect() };
        // The full step, with its tangent (taken in most iterations).
        let full = self.state(&at(1.0), lambda, true);
        if !self.settings.line_search {
            return full.map(|e| (1.0, e));
        }
        let red = self.problem.red;
        let work = |e: &Evaluation| dot(dq, &free_residual(red, e));
        let mut s = 1.0;
        let mut g = match &full {
            Ok(e) => {
                let g1 = work(e);
                if g1.abs() <= LINE_SEARCH_ACCEPT * g0.abs() {
                    return full.map(|e| (1.0, e));
                }
                Some(g1)
            }
            Err(Failure::Cutback(_)) => None,
            Err(Failure::Fatal(_)) => return full.map(|e| (1.0, e)),
        };
        let mut best: Option<(f64, f64)> = g.map(|g| (1.0, g.abs()));
        for _ in 0..LINE_SEARCH_TRIALS {
            // Root of the line through (0, g0) and (s, g), or half the step after a failure.
            let next = match g {
                Some(gs) if (g0 - gs).abs() > 0.0 => (s * g0 / (g0 - gs)).clamp(LINE_SEARCH_MIN, 1.0),
                Some(_) => break,
                None => (0.5 * s).max(LINE_SEARCH_MIN),
            };
            if next == s {
                break;
            }
            s = next;
            g = match self.state(&at(s), lambda, false) {
                Ok(e) => Some(work(&e)),
                Err(Failure::Cutback(_)) => None,
                Err(f) => return Err(f),
            };
            if let Some(gs) = g {
                if best.is_none_or(|(_, b)| gs.abs() < b) {
                    best = Some((s, gs.abs()));
                }
                if gs.abs() <= LINE_SEARCH_ACCEPT * g0.abs() {
                    break;
                }
            }
        }
        let Some((s, _)) = best else {
            return Err(Failure::Cutback("every step along the Newton correction inverts an element".into()));
        };
        if s == 1.0 {
            return full.map(|e| (1.0, e));
        }
        self.state(&at(s), lambda, true).map(|e| (s, e))
    }
}

/// Refuse what nonlinear steps cannot take yet, and warn about what they simplify.
fn check_supported(model: &Model, mesh: &Mesh, step: &Step, settings: &NonlinearStatic, warnings: &mut Vec<String>) -> Result<()> {
    let at = |what: String| Error::NotSupported(format!("{what} in nonlinear steps (step `{}`)", step.name));
    if let Some(t) = mesh.trusses.first() {
        return Err(at(format!("trusses (element {})", t.label)));
    }
    if let Some(s) = mesh.structural.first() {
        return Err(at(format!("shells and beams, whose large rotations are not implemented yet (element {})", s.label)));
    }
    for e in &mesh.elements {
        if !settings.nonlinear_geometry && e.hyperelastic.is_some() {
            return Err(at(format!("a hyperelastic material without nonlinear geometry (element {})", e.label)));
        }
    }
    if let Some(name) = step.loads.iter().find(|n| matches!(model.load(n), Some(Load::Centrifugal { .. }))) {
        return Err(at(format!("centrifugal loads (load `{name}`), which depend on the deformed position")));
    }
    plasticity(model, step, settings, warnings)?;
    if let Some(name) = model.contacts().keys().next()
        && settings.arc_length.is_some()
    {
        return Err(at(format!("contact (`{name}`) with arc-length control")));
    }
    if settings.nonlinear_geometry
        && let Some(name) = model.couplings().keys().next()
    {
        warnings.push(format!(
            "coupling `{name}` and the others are linear constraints (small rotations of their reference points) in this nonlinear step"
        ));
    }
    Ok(())
}

/// Refuse the plasticity a static step cannot take, and warn about what it leaves out, for the
/// materials of the solid and 2D sections.
fn plasticity(model: &Model, step: &Step, settings: &NonlinearStatic, warnings: &mut Vec<String>) -> Result<()> {
    let at = |what: String| Error::NotSupported(format!("{what} in nonlinear steps (step `{}`)", step.name));
    let used: std::collections::BTreeSet<&str> = model
        .sections()
        .values()
        .filter(|s| matches!(s.kind, SectionKind::Solid {}) || Theory::from_section(&s.kind).is_some())
        .map(|s| s.material.as_str())
        .collect();
    for name in used {
        let Some(m) = model.material(name) else { continue };
        match &m.plastic {
            Some(Plastic::JohnsonCook { c, specific_heat, .. }) => {
                if *c != 0.0 {
                    return Err(at(format!(
                        "the Johnson–Cook strain-rate term of material `{name}` (C = {c}): a static step has no time, so no strain rate; set C to 0 or use explicit dynamics"
                    )));
                }
                if specific_heat.is_some() {
                    warnings.push(format!(
                        "material `{name}`: a static step is isothermal, so plastic work does not heat it (adiabatic heating acts in explicit dynamics)"
                    ));
                }
            }
            Some(Plastic::Bilinear { hardening: Hardening::Kinematic, .. }) if settings.nonlinear_geometry => {
                return Err(at(format!(
                    "the kinematic hardening of material `{name}` at finite strain (its back stress is not coaxial with the elastic stretch, which the exponential return map needs); use nonlinear_geometry false or isotropic hardening"
                )));
            }
            _ => {}
        }
        if m.plastic.is_some() && (m.damage.is_some() || m.erosion.is_some()) {
            warnings.push(format!("material `{name}`: damage and erosion act in explicit dynamics; this static step leaves them out"));
        }
    }
    Ok(())
}

/// The step's temperature field per mesh node: the steady heat solution of its temperature
/// boundary conditions and heat loads when it has some (one-way coupling, as the coupled
/// temperature–displacement step of the linear solver), or its temperature loads.
fn temperatures(model: &Model, resolved: &fieldcraft_model::Resolved, mesh: &Mesh, step: &Step, limits: &Limits) -> Result<Option<Vec<Option<f64>>>> {
    let conducts = step.boundary_conditions.iter().any(|b| model.boundary_conditions().get(b).is_some_and(BoundaryCondition::is_thermal))
        || step.loads.iter().any(|l| model.load(l).is_some_and(Load::is_thermal));
    let loads = temperature_loads(model, resolved, mesh, step)?;
    match (conducts, loads) {
        (true, Some(_)) => Err(Error::Invalid(format!(
            "step `{}` gives temperatures both by temperature loads and by heat conditions (temperature boundary conditions, heat loads); use one",
            step.name
        ))),
        (true, None) => Ok(Some(crate::heat::temperatures_at(model, resolved, step, limits, &mesh.labels)?)),
        (false, field) => Ok(field),
    }
}

/// The follower pressures of a step (with nonlinear geometry), and the step without them for
/// the dead loads.
fn followers(model: &Model, resolved: &fieldcraft_model::Resolved, mesh: &Mesh, step: &Step) -> Result<(Vec<Follower>, Step)> {
    let mut out = Vec::new();
    let mut dead = step.clone();
    for name in &step.loads {
        let Some(Load::Pressure { set, magnitude }) = model.load(name) else { continue };
        dead.loads.remove(name);
        let slots = |labels: Vec<u32>| -> Result<(Vec<usize>, Vec<[f64; 3]>, Vec<usize>)> {
            let idx: Vec<usize> = labels
                .iter()
                .map(|l| mesh.node(*l).ok_or_else(|| Error::Invalid(format!("load `{name}`: node {l} belongs to no solid element"))))
                .collect::<Result<_>>()?;
            let xyz = idx.iter().filter_map(|i| mesh.xyz.get(*i).copied()).collect();
            let mut dofs = Vec::with_capacity(3 * idx.len());
            for n in &idx {
                dofs.extend(translations(mesh, *n)?);
            }
            Ok((idx, xyz, dofs))
        };
        if let Some(theory) = mesh.theory {
            if resolved.sets.get(set).is_some_and(|s| !s.faces.is_empty()) {
                return Err(Error::Invalid(format!("load `{name}`: a 2D model's pressures act on element edges, and set `{set}` has element faces")));
            }
            for g in resolved.planar_edges(set) {
                let labels =
                    model.edge_nodes(g).ok_or_else(|| Error::Invalid(format!("load `{name}`: element {} edge {} is missing", g.element, g.edge)))?;
                let (nodes, xyz, dofs) = slots(labels)?;
                out.push(Follower { element: g.element, nodes, xyz, magnitude: *magnitude, extent: Some(theory.extent()), dofs });
            }
        } else {
            for fc in resolved.sets.get(set).map(|s| s.faces.as_slice()).unwrap_or_default() {
                let labels = model
                    .face_nodes(*fc)
                    .ok_or_else(|| Error::Invalid(format!("load `{name}`: element {} face {} is missing", fc.element, fc.face)))?;
                let (nodes, xyz, dofs) = slots(labels)?;
                out.push(Follower { element: fc.element, nodes, xyz, magnitude: *magnitude, extent: None, dofs });
            }
        }
    }
    Ok((out, dead))
}

/// Solve a nonlinear static step: the step named `step`, or the model's only nonlinear static
/// step.
pub fn solve_nonlinear(model: &Model, step: Option<&str>, limits: &Limits) -> Result<StaticResult> {
    // A job cancelled before it starts stops here, before any setup work.
    fieldcraft_model::progress::check()?;
    let chosen = model.select_step(step, |p| matches!(p, Procedure::NonlinearStatic(_)), "nonlinear static")?;
    let Procedure::NonlinearStatic(settings) = &chosen.procedure else {
        return Err(Error::Internal("nonlinear static step".into()));
    };
    settings.validate()?;
    let Some(imperfection) = &settings.imperfection else { return solve_on(model, chosen.name.as_str(), limits, Vec::new()) };
    let (perturbed, note) = imperfect(model, imperfection, limits)?;
    solve_on(&perturbed, chosen.name.as_str(), limits, vec![note])
}

/// The model with a buckling mode, scaled to a largest component of `scale`, added to its node
/// coordinates (the buckling step solved first), and a note saying so.
fn imperfect(model: &Model, imperfection: &Imperfection, limits: &Limits) -> Result<(Model, String)> {
    let Imperfection { step, mode, scale } = imperfection;
    if !matches!(model.step(step).map(|s| &s.procedure), Some(Procedure::Buckle { .. })) {
        return Err(Error::Invalid(format!("the imperfection names step `{step}`, which is not a buckling step")));
    }
    let buckled = crate::buckle::solve_buckle(model, Some(step), limits)?;
    let shape =
        usize::try_from(*mode).ok().and_then(|m| m.checked_sub(1)).and_then(|m| buckled.modes.get(m)).ok_or_else(|| {
            Error::Invalid(format!("the imperfection takes mode {mode} of step `{step}`, which found {} modes", buckled.modes.len()))
        })?;
    let mut out = model.clone();
    for (label, u) in buckled.nodes.iter().zip(&shape.shape) {
        let node = model.node(*label).ok_or_else(|| Error::Internal(format!("node {label}")))?;
        let xyz: [f64; 3] = std::array::from_fn(|k| node.xyz[k] + scale * u[k]);
        out.set_node(*label, xyz)?;
    }
    let note = format!(
        "solved on an imperfect geometry: mode {mode} of buckling step `{step}` (load factor {:.6e}) scaled to {scale} added to the node coordinates",
        shape.eigenvalue
    );
    Ok((out, note))
}

fn solve_on(model: &Model, step: &str, limits: &Limits, notes: Vec<String>) -> Result<StaticResult> {
    let resolved = model.resolve()?;
    let step = model.step(step).ok_or_else(|| Error::Internal("nonlinear static step".into()))?;
    let Procedure::NonlinearStatic(settings) = &step.procedure else {
        return Err(Error::Internal("nonlinear static step".into()));
    };
    let mesh = Mesh::new(model, &resolved, limits.max_dofs)?;
    let mut warnings: Vec<String> =
        notes.into_iter().chain(resolved.warnings.iter().map(ToString::to_string)).chain(mesh.warnings.iter().cloned()).collect();
    check_supported(model, &mesh, step, settings, &mut warnings)?;
    let (graph, order) = pattern(&mesh);
    ignored_dashpots(&mesh, &mut warnings);
    let setup = Setup { model, resolved: &resolved, mesh: &mesh, step };
    let red = reduction(&setup, true, &order, budget(limits, &mesh, &graph, 1), &mut warnings)?;
    if red.neq() == 0 {
        return Err(Error::Invalid("every degree of freedom is prescribed; there is nothing to solve".into()));
    }
    let (follower_loads, dead_step) =
        if settings.nonlinear_geometry { followers(model, &resolved, &mesh, step)? } else { (Vec::new(), step.clone()) };
    let dead = loads(model, &resolved, &mesh, &dead_step, &red.active)?;
    let field = temperatures(model, &resolved, &mesh, step, limits)?;
    let springs = if mesh.discrete.springs.iter().any(|s| !s.dashpot) {
        let mut s = graph.zero_matrix(&mesh.layout);
        mesh.discrete.add_stiffness(&mut s)?;
        Some(s)
    } else {
        None
    };
    let contacts = Interface::all(model, &resolved, &mesh, settings.nonlinear_geometry, &mut warnings)?;
    let element_kinds = mesh.elements.iter().map(|e| Ok(e.kernel()?.dofs().node)).collect::<Result<Vec<_>>>()?;
    let element_dofs = mesh
        .elements
        .iter()
        .zip(&element_kinds)
        .map(|(e, kinds)| {
            mesh.layout.element_dofs(&e.nodes, *kinds).ok_or_else(|| Error::Internal(format!("element {} degrees of freedom", e.label)))
        })
        .collect::<Result<Vec<_>>>()?;
    let mut problem = Problem {
        mesh: &mesh,
        red: &red,
        dead,
        followers: follower_loads,
        springs,
        finite: settings.nonlinear_geometry,
        element_dofs,
        element_kinds,
        materials: mesh
            .elements
            .iter()
            .map(|e| match &e.plastic {
                Some(j2) => MaterialModel::Plastic(j2.clone()),
                None => MaterialModel::from_elasticity(&e.d, e.hyperelastic.as_ref()),
            })
            .collect(),
        field,
        pattern: graph.zero_matrix(&mesh.layout),
        graph,
        contacts,
    };
    // Elements on their linear stiffness take the thermal strain as the equivalent nodal
    // forces, a dead load proportional to the load factor like the temperatures.
    if let Some(field) = problem.field.take() {
        let linear: Vec<bool> = (0..mesh.elements.len()).map(|i| problem.kinematics(i).is_none()).collect();
        add_element_thermal_loads(&mut problem.dead, &mesh, &field, |i| linear.get(i).copied().unwrap_or(false))?;
        problem.field = Some(field);
    }
    let (lo, hi) = mesh.xyz.iter().fold(([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]), |(lo, hi), p| {
        (std::array::from_fn(|k| lo[k].min(p[k])), std::array::from_fn(|k| hi[k].max(p[k])))
    });
    let length = (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt().max(f64::MIN_POSITIVE);
    let committed = mesh
        .elements
        .iter()
        .zip(&problem.materials)
        .map(|(e, m)| Ok(vec![0.0; e.kernel()?.rule().points().len() * m.history_len()]))
        .collect::<Result<Vec<_>>>()?;
    // The undeformed state first: an element that is bad before anything moves is an error in
    // the model, not a reason to cut increments back.
    problem.evaluate(&vec![0.0; mesh.dofs()], 0.0, false, &committed, &[], &[], 0.0)?;
    let path_dependent = problem.materials.iter().any(MaterialModel::is_path_dependent);
    // Rotations (radians, their forces moments) are measured apart from translations.
    let rotation: Vec<bool> = (0..mesh.dofs()).map(|d| matches!(mesh.layout.locate(d), Some((_, DofKind::Rx | DofKind::Ry | DofKind::Rz)))).collect();
    let is_rotation = |d: usize| rotation.get(d).copied().unwrap_or(false);
    let mut run = Run {
        problem: &problem,
        settings,
        limits,
        force_scale: 0.0,
        force_floor: None,
        moment_scale: 0.0,
        length,
        prescribed: (red.largest_prescribed_where(|d| !is_rotation(d)), red.largest_prescribed_where(is_rotation)),
        rotation_eq: red.eq_dof.iter().map(|d| is_rotation(*d)).collect(),
        rotation,
        iterations: 0,
        factor_entries: 0,
        committed,
        path_dependent,
        undeformed: None,
        contact: problem.contacts.iter().map(Interface::committed).collect(),
        reference: vec![0.0; mesh.dofs()],
        band: 0.0,
        augmentations: 0,
        augmentation: 0,
        converged_scale: 0.0,
        converged_moment: 0.0,
    };

    let mut history = History::default();
    let (mut frame_values, mut frames_full) = (0usize, false);
    let values = 3 * mesh.labels.len();
    let mut keep_frame = |history: &mut History, lambda: f64, u: &[f64], at: Option<&Evaluation>, warnings: &mut Vec<String>| {
        if frames_full {
            return;
        }
        let reaction: Vec<(usize, [f64; 3])> = at
            .map(|e| {
                let out_of_balance: Vec<f64> = e.internal.iter().zip(&e.external).map(|(a, b)| a - b).collect();
                split(&red.reactions_of(&out_of_balance), &mesh).0.into_iter().enumerate().filter(|(_, f)| f.iter().any(|c| *c != 0.0)).collect()
            })
            .unwrap_or_default();
        let these = values.saturating_add(reaction.len().saturating_mul(4));
        if frame_values.saturating_add(these) > MAX_FRAME_VALUES {
            warnings.push(format!("the frames stop at load factor {lambda}: the step keeps at most {MAX_FRAME_VALUES} values"));
            frames_full = true;
            return;
        }
        frame_values += these;
        history.frames.push(Frame { load_factor: lambda, displacement: split(u, &mesh).0, reaction });
    };
    keep_frame(&mut history, 0.0, &vec![0.0; mesh.dofs()], None, &mut warnings);
    let (q, lambda, last) = if let Some(arc) = &settings.arc_length {
        if problem.field.is_some() {
            return Err(Error::NotSupported(format!("temperatures in arc-length steps (step `{}`)", step.name)));
        }
        arc_length(&mut run, arc, step, &mut history, &mut warnings, &mut keep_frame)?
    } else {
        let mut q = vec![0.0; red.neq()];
        // The step fraction t and the load factor λ(t) of the step's load path (λ = t without one).
        let (mut t, mut lambda) = (0.0f64, 0.0f64);
        let mut size = settings.initial_increment;
        let mut fast = 0usize;
        let mut last: Option<Evaluation> = None;
        // The last converged increment's free motion and the change of the load factor over it.
        let mut previous: Option<(Vec<f64>, f64)> = None;
        let max_increments = usize::try_from(settings.max_increments).unwrap_or(usize::MAX);
        fieldcraft_model::progress::stage("Load factor", Some(1.0), "");
        while t < 1.0 {
            fieldcraft_model::progress::count(lambda);
            fieldcraft_model::progress::check()?;
            if history.increments.len() >= max_increments {
                return Err(Error::Limit(format!(
                    "step `{}`: {} increments reached the step's cap at load factor {lambda:.6} (raise max_increments)",
                    step.name, settings.max_increments
                )));
            }
            // An increment ends at the next point of the load path at the latest.
            let turn = settings.next_turn(t);
            let remaining = turn - t;
            let (target, this) = if size >= remaining * (1.0 - SAME_TIME) { (turn, remaining) } else { (t + size, size) };
            let target_lambda = settings.load_factor(target);
            let change = if settings.load_path.is_some() { target_lambda - lambda } else { this };
            let mut iterations = Vec::new();
            let number = history.increments.len() + 1;
            // Linear extrapolation of the last converged increment's free motion, in proportion to
            // the change of the load factor (so it turns with the load path).
            let mut predictor: Option<Vec<f64>> = previous
                .as_ref()
                .filter(|(_, prev)| *prev != 0.0)
                .map(|(dq, prev): &(Vec<f64>, f64)| q.iter().zip(dq).map(|(a, d)| a + change / prev * d).collect());
            // Path-dependent materials start without one from the tangent predictor.
            let mut predicted_energy = None;
            if predictor.is_none() && path_dependent {
                match run.tangent_predictor(&q, lambda, target_lambda) {
                    Ok(Some((p, energy))) => (predictor, predicted_energy) = (Some(p), Some(energy)),
                    Ok(None) | Err(Failure::Cutback(_)) => {}
                    Err(Failure::Fatal(e)) => return Err(e),
                }
            }
            let saved = run.contact.clone();
            let mut augmentations = 0;
            run.augmentation = 0;
            let solved = run
                .increment(&q, predictor.as_deref(), predicted_energy, target_lambda, change.abs(), &mut iterations)
                .and_then(|c| run.augment(c, &q, target_lambda, change.abs(), &mut iterations, &mut augmentations, &mut warnings));
            if solved.is_err() {
                run.contact = saved;
            }
            match solved {
                Ok(c) => {
                    previous = Some((c.q.iter().zip(&q).map(|(a, b)| a - b).collect(), change));
                    // The increment's own Newton iterations size the next one (augmentations
                    // solve the same increment again).
                    let count = iterations.iter().filter(|i| i.augmentation == 0).count();
                    history.increments.push(Increment {
                        number,
                        load_factor: target_lambda,
                        size: this,
                        iterations,
                        converged: true,
                        cutback: None,
                        negative_pivots: 0,
                        augmentations,
                    });
                    (t, lambda) = (target, target_lambda);
                    q = c.q;
                    keep_frame(&mut history, lambda, &red.expand_scaled(&q, lambda), Some(&c.evaluation), &mut warnings);
                    run.commit_increment(&q, lambda, &c.evaluation);
                    last = Some(c.evaluation);
                    fast = if count <= FAST_ITERATIONS { fast + 1 } else { 0 };
                    if fast >= 2 {
                        size *= GROWTH;
                    } else if count > SLOW_ITERATIONS {
                        size *= SHRINK;
                    }
                    size = size.clamp(settings.min_increment, settings.max_increment);
                }
                Err(Failure::Fatal(e)) => return Err(e),
                Err(Failure::Cutback(reason)) => {
                    history.increments.push(Increment {
                        number,
                        load_factor: target_lambda,
                        size: this,
                        iterations,
                        converged: false,
                        cutback: Some(reason.clone()),
                        negative_pivots: 0,
                        augmentations: 0,
                    });
                    fast = 0;
                    size = this * CUTBACK;
                    if size < settings.min_increment * (1.0 - SAME_TIME) {
                        return Err(Error::NotConverged(format!(
                            "step `{}` stopped at load factor {lambda:.6} (step fraction {t:.6}): the increment to {target_lambda:.6} failed ({reason}) and a smaller one would be below the minimum increment {}",
                            step.name, settings.min_increment
                        )));
                    }
                }
            }
        }
        (q, lambda, last)
    };
    history.iterations = run.iterations;
    let u = red.expand_scaled(&q, lambda);
    let evaluation = match last {
        Some(e) => e,
        None => problem.evaluate(&u, lambda, false, &run.committed, &run.contact, &run.reference, 0.0)?,
    };
    if let Some((at, negative)) = history
        .increments
        .iter()
        .filter(|_| settings.arc_length.is_none())
        .filter(|i| i.converged)
        .filter_map(|i| i.iterations.last().map(|it| (i.load_factor, it.negative_pivots)))
        .find(|(_, n)| *n > 0)
    {
        warnings.push(format!(
            "the converged tangent had {negative} negative eigenvalue(s) from load factor {at:.6} on: an unstable equilibrium under load control (past a limit or bifurcation point, or an unstable material)"
        ));
    }
    if let Some(a) = problem.follower_asymmetry(&u, lambda)?
        && a > NOTABLE_ASYMMETRY
    {
        warnings.push(format!(
            "follower pressure: the skew part of its load stiffness ({a:.1e} of it, over the free degrees of freedom) is left out of the tangent; the solution is exact, Newton may converge more slowly"
        ));
    }
    let contacts = problem.contacts.iter().zip(&evaluation.contact).map(|(c, s)| c.result(&mesh, s, run.augmentations)).collect();
    finish(&mesh, step, &red, &problem, u, lambda, evaluation, &run.committed, history, run.factor_entries, contacts, warnings)
}

/// The final state as a static result: displacements, reactions, Cauchy stresses and
/// logarithmic strains.
#[allow(clippy::too_many_arguments)]
fn finish(
    mesh: &Mesh,
    step: &Step,
    red: &Reduction,
    problem: &Problem,
    u: Vec<f64>,
    lambda: f64,
    e: Evaluation,
    committed: &[Vec<f64>],
    history: History,
    factor_entries: usize,
    contacts: Vec<ContactResult>,
    mut warnings: Vec<String>,
) -> Result<StaticResult> {
    let (displacement, rotation, point_motion) = split(&u, mesh);
    let out_of_balance: Vec<f64> = e.internal.iter().zip(&e.external).map(|(a, b)| a - b).collect();
    let reactions = red.reactions_of(&out_of_balance);
    let (reaction, reaction_moment, point_reaction) = split(&reactions, mesh);
    let (mut applied_total, mut reaction_total) = ([0.0; 3], [0.0; 3]);
    for (dof, (fv, rv)) in e.external.iter().zip(&reactions).enumerate() {
        let Some(d) = mesh.kind(dof).and_then(DofKind::direction) else { continue };
        if let (Some(a), Some(r)) = (applied_total.get_mut(d), reaction_total.get_mut(d)) {
            *a += fv;
            *r += rv;
        }
    }
    let reference_points = mesh
        .points
        .iter()
        .zip(point_motion.iter().zip(&point_reaction))
        .map(|(p, (d, r))| PointResult { name: p.name.clone(), label: p.label, displacement: *d, reaction: *r })
        .collect();
    #[allow(clippy::type_complexity)]
    let all_points = parallel::map_grained(mesh.elements.len(), GRAIN, |i| -> Result<(Vec<PointStress>, Vec<f64>, &'static ElementType, Rule)> {
        let el = mesh.elements.get(i).ok_or_else(|| Error::Internal("element index".into()))?;
        let at = |err| Error::element(el.label, err);
        let ue = problem.element_vector(i, &u);
        let kernel = el.kernel()?;
        let history = committed.get(i).map(Vec::as_slice).unwrap_or_default();
        let material = problem.materials.get(i).ok_or_else(|| Error::Internal("element material".into()))?;
        let points = match problem.kinematics(i) {
            Some(kin) => kernel.nonlinear_stresses(material, &ue, kin, problem.temperature(i, lambda).as_ref(), history),
            None => {
                let per_node: Vec<[f64; 3]> = el.nodes.iter().map(|n| displacement.get(*n).copied().unwrap_or_default()).collect();
                let thermal = problem.field.as_deref().and_then(|f| el.thermal_strain(f)).map(|mut t| {
                    t.delta.iter_mut().for_each(|d| *d *= lambda);
                    t
                });
                kernel.stresses_with(&el.d, &per_node, thermal.as_ref())
            }
        }
        .map_err(at)?;
        let len = material.history_len();
        let plastic: Vec<f64> = if len > 0 { history.chunks(len).filter_map(|h| material.plastic_strain(h)).collect() } else { Vec::new() };
        Ok((points, plastic, kernel.ty(), kernel.rule()))
    });
    let mut extrapolations: Vec<(&'static ElementType, Rule, Extrapolation)> = Vec::new();
    let mut elements = Vec::with_capacity(mesh.elements.len());
    let mut sum = vec![[0.0; 6]; mesh.labels.len()];
    let mut count = vec![0u32; mesh.labels.len()];
    for (el, computed) in mesh.elements.iter().zip(all_points) {
        let (points, plastic_strain, ty, rule) = computed?;
        let values: Vec<[f64; 6]> = points.iter().map(|p| p.stress).collect();
        let index = match extrapolations.iter().position(|(t, r, _)| std::ptr::eq(*t, ty) && *r == rule) {
            Some(i) => i,
            None => {
                extrapolations.push((ty, rule, extrapolation(ty, rule).map_err(|err| Error::element(el.label, err))?));
                extrapolations.len() - 1
            }
        };
        let nodal = extrapolations.get(index).map(|(_, _, x)| x.apply(&values)).unwrap_or_default();
        for (n, s) in el.nodes.iter().zip(&nodal) {
            if let (Some(acc), Some(c)) = (sum.get_mut(*n), count.get_mut(*n)) {
                for k in 0..6 {
                    acc[k] += s[k];
                }
                *c += 1;
            }
        }
        elements.push(ElementStress { label: el.label, kind: el.kind, points, nodal, plastic_strain });
    }
    let stress: Vec<[f64; 6]> = sum.iter().zip(&count).map(|(s, c)| s.map(|v| if *c > 0 { v / f64::from(*c) } else { 0.0 })).collect();
    let von_mises = stress.iter().map(von_mises).collect();
    let scale = crate::system::balance_scale(&applied_total, &reaction_total, &reaction);
    let ground = mesh.discrete.ground_force(&u, mesh);
    let imbalance = (0..3).map(|c| (applied_total[c] + reaction_total[c] + ground[c]).abs()).fold(0.0, f64::max);
    if scale > 0.0 && imbalance > crate::system::BALANCE_TOLERANCE * scale {
        warnings.push(format!("loads and reactions differ by {imbalance:.3e} (relative {:.1e})", imbalance / scale));
    }
    Ok(StaticResult {
        step: step.name.clone(),
        nodes: mesh.labels.clone(),
        displacement,
        rotation,
        reaction,
        reaction_moment,
        stress,
        von_mises,
        elements,
        shells: Vec::new(),
        beams: Vec::new(),
        reference_points,
        applied_total,
        reaction_total,
        stats: Stats { dofs: mesh.dofs(), equations: red.neq(), factor_entries, lanczos_vectors: 0 },
        warnings,
        temperature: None,
        nonlinear: Some(history),
        contacts,
    })
}

/// What one arc-length increment starts from: the tangent at the converged state, factored, and
/// the tangent displacement of the unit load.
struct ArcStart {
    factor: Symmetric,
    load_displacement: Vec<f64>,
}

impl Run<'_> {
    /// The arc length's inner product of two motions over the free equations: translations as
    /// they are, rotations times L² (the translation a rotation gives at the model's size), so
    /// the arc is a length² in any unit system and at any size of the part. Without rotations it
    /// is the plain dot product, term for term.
    fn arc_dot(&self, a: &[f64], b: &[f64]) -> f64 {
        let lever = self.length * self.length;
        a.iter().zip(b).enumerate().map(|(i, (x, y))| if self.rotation_eq.get(i).copied().unwrap_or(false) { x * y * lever } else { x * y }).sum()
    }

    /// The tangent at (q, λ) factored (Cholesky, or LDLᵀ past a critical point), the unit load
    /// ∂r/∂λ at that state over the free equations (the loads at factor 1, less K T_p u_p for
    /// prescribed motions), and the evaluation's residual.
    fn arc_tangent(&mut self, e: &mut Evaluation, u: &[f64], first: bool) -> std::result::Result<(Symmetric, Vec<f64>), Failure> {
        let red = self.problem.red;
        let k = e.tangent.take().ok_or_else(|| Failure::Fatal(Error::Internal("tangent".into())))?;
        let mut g = red.restrict(&self.problem.reference_load(u)?);
        if red.largest_prescribed() > 0.0 {
            g.iter_mut().zip(prescribed_force(red, &k, 1.0)).for_each(|(a, b)| *a -= b);
        }
        let [kff] = red.reduce([&k])?;
        drop(k);
        if first {
            self.seed_force_floor(&kff);
        }
        let factor = self.factor(&kff, first)?;
        self.factor_entries = factor.entries();
        Ok((factor, g))
    }

    /// The converged state's tangent (from its evaluation `e`, which keeps its forces) and the
    /// unit load's tangent displacement.
    fn arc_start(&mut self, e: &mut Evaluation, q: &[f64], lambda: f64, first: bool) -> std::result::Result<ArcStart, Failure> {
        let (factor, mut g) = self.arc_tangent(e, &self.problem.red.expand_scaled(q, lambda), first)?;
        factor.solve(&mut g)?;
        if g.iter().any(|v| !v.is_finite()) || norm_inf(&g) == 0.0 {
            return Err(Failure::Fatal(Error::Invalid("the step's loads do not move the model, so an arc length cannot be measured".into())));
        }
        Ok(ArcStart { factor, load_displacement: g })
    }

    /// One arc-length increment from the converged (q0, λ0) (M. A. Crisfield, *A fast
    /// incremental/iterative solution procedure that handles "snap-through"*, Computers &
    /// Structures 13:55–62, 1981): the constraint Δqᵀ W Δq + w Δλ² = Δl² (W the identity on
    /// translations and L² on rotations, [`Run::arc_dot`]; w = 0 cylindrical, w = ‖q_f‖²_W
    /// spherical) with the load factor an unknown. Predictor along the tangent, its sign
    /// by the angle with the last increment (M. Feng, D. Perić, D. R. J. Owen, *A new criterion for
    /// determination of initial loading parameter in arc-length methods*, Computers & Structures
    /// 58(3):479–485, 1996); each correction δq = δq_r + δλ δq_g solves the quadratic for δλ and
    /// takes the root whose increment is closest in direction to the previous iterate's (the
    /// angle criterion, Crisfield 1981).
    #[allow(clippy::too_many_arguments)]
    fn arc_increment(
        &mut self,
        q0: &[f64],
        lambda0: f64,
        start: &ArcStart,
        arc: f64,
        weight: f64,
        previous: Option<&(Vec<f64>, f64)>,
        iterations: &mut Vec<Iteration>,
    ) -> std::result::Result<(Vec<f64>, f64, Evaluation), Failure> {
        let tol = self.settings.tolerances;
        let qg = &start.load_displacement;
        let sign = match previous {
            Some((dq, dl)) if self.arc_dot(dq, qg) + weight * dl < 0.0 => -1.0,
            _ => 1.0,
        };
        let mut dlam = sign * arc / (self.arc_dot(qg, qg) + weight).sqrt();
        let mut dq: Vec<f64> = qg.iter().map(|g| dlam * g).collect();
        let mut correction = f64::INFINITY;
        let mut first_energy: Option<f64> = None;
        let max_iterations = usize::try_from(self.settings.max_iterations).unwrap_or(usize::MAX);
        let at = |dq: &[f64]| -> Vec<f64> { q0.iter().zip(dq).map(|(a, b)| a + b).collect() };
        for it in 0..=max_iterations {
            let q = at(&dq);
            let mut e = self.state(&q, lambda0 + dlam, true)?;
            let r = self.residual(&e);
            let force = self.force_measure(&r);
            if !force.is_finite() {
                return Err(Failure::Cutback("the residual is not finite".into()));
            }
            if it > 0 && ((force <= tol.force && correction <= tol.displacement) || force <= ROUND_OFF * tol.force) {
                return Ok((q, lambda0 + dlam, e));
            }
            if it == max_iterations {
                break;
            }
            if self.iterations >= self.limits.max_newton_iterations {
                return Err(Failure::Fatal(Error::Limit(format!(
                    "the step took {} Newton iterations, the limit (raise max_newton_iterations or loosen the tolerances)",
                    self.iterations
                ))));
            }
            self.iterations += 1;
            let u = self.problem.red.expand_scaled(&q, lambda0 + dlam);
            let (factor, mut dqg) = self.arc_tangent(&mut e, &u, false)?;
            let mut dqr = r.clone();
            factor.solve(&mut dqr)?;
            factor.solve(&mut dqg)?;
            // (Δq + δq_r + δλ δq_g)² + w (Δλ + δλ)² = Δl².
            let base: Vec<f64> = dq.iter().zip(&dqr).map(|(a, b)| a + b).collect();
            let a = self.arc_dot(&dqg, &dqg) + weight;
            let b = 2.0 * (self.arc_dot(&dqg, &base) + weight * dlam);
            let c = self.arc_dot(&base, &base) + weight * dlam * dlam - arc * arc;
            let disc = b * b - 4.0 * a * c;
            if !(disc >= 0.0) || !(a > 0.0) {
                return Err(Failure::Cutback("the arc-length constraint has no real root (the arc is too long here)".into()));
            }
            let root = disc.sqrt();
            // Stable roots of a x² + b x + c (no cancellation).
            let qr = -0.5 * (b + b.signum() * root);
            let roots = if qr == 0.0 { [0.0, 0.0] } else { [qr / a, c / qr] };
            let cosine = |dl: f64| self.arc_dot(&dq, &base) + dl * self.arc_dot(&dq, &dqg) + weight * dlam * (dlam + dl);
            let dl = if cosine(roots[0]) >= cosine(roots[1]) { roots[0] } else { roots[1] };
            let step: Vec<f64> = dqr.iter().zip(&dqg).map(|(r, g)| r + dl * g).collect();
            if step.iter().any(|v| !v.is_finite()) {
                return Err(Failure::Cutback("the arc-length correction is not finite".into()));
            }
            dq.iter_mut().zip(&step).for_each(|(a, s)| *a += s);
            dlam += dl;
            let (translated, rotated) = norm_inf_split(&dq, &self.rotation_eq);
            correction = self.correction_measure(&step, 1.0, self.motion_references(translated, rotated));
            let energy = dot(&step, &r).abs();
            let first = *first_energy.get_or_insert(energy);
            iterations.push(Iteration {
                force,
                displacement: correction,
                energy: if first > 0.0 { energy / first } else { 0.0 },
                line_search: 1.0,
                negative_pivots: factor.inertia().negative,
                augmentation: 0,
            });
        }
        Err(Failure::Cutback(format!("no convergence in {} iterations", self.settings.max_iterations)))
    }
}

/// Follow a step's equilibrium path by arc length (Riks 1979, Crisfield 1981): the load factor
/// is an unknown, so the path goes through limit points. The first increment's load factor
/// change is the step's initial increment, which sets the arc length; increments grow and
/// shrink as under load control (the arc length in proportion, within the minimum and maximum
/// increments' arcs). The tangent at every converged state is factored for the next predictor
/// (Cholesky, LDLᵀ when indefinite); its negative eigenvalues are recorded per increment, and
/// where their count changes a critical point is reported: a limit point when the load factor
/// turns there, a bifurcation otherwise. The step ends at the maximum load factor, at the
/// displacement limit, or at its increment cap (with a warning).
#[allow(clippy::type_complexity)]
fn arc_length(
    run: &mut Run,
    arc: &ArcLength,
    step: &Step,
    history: &mut History,
    warnings: &mut Vec<String>,
    keep_frame: &mut dyn FnMut(&mut History, f64, &[f64], Option<&Evaluation>, &mut Vec<String>),
) -> Result<(Vec<f64>, f64, Option<Evaluation>)> {
    let settings = run.settings;
    let red = run.problem.red;
    let mesh = run.problem.mesh;
    let limit = match &arc.displacement_limit {
        Some(l) => {
            let slot =
                mesh.node(l.node).ok_or_else(|| Error::Invalid(format!("the displacement limit's node {} has no degrees of freedom", l.node)))?;
            let dofs = translations(mesh, slot)?;
            let dof = dofs
                .get(usize::from(l.direction).saturating_sub(1))
                .copied()
                .ok_or_else(|| Error::Invalid(format!("the displacement limit's node {} has no direction {}", l.node, l.direction)))?;
            Some((dof, l.displacement.abs()))
        }
        None => None,
    };
    let fatal = |f: Failure| match f {
        Failure::Fatal(e) => e,
        Failure::Cutback(m) => Error::NotConverged(format!("step `{}`: {m}", step.name)),
    };
    let mut q = vec![0.0; red.neq()];
    let mut lambda = 0.0f64;
    let mut initial = run.state(&q, lambda, true).map_err(fatal)?;
    let mut start = run.arc_start(&mut initial, &q, lambda, true).map_err(fatal)?;
    // States: load factor and negative eigenvalues, from the undeformed one.
    let mut states: Vec<(usize, f64, usize)> = vec![(0, 0.0, start.factor.inertia().negative)];
    let scale = (run.arc_dot(&start.load_displacement, &start.load_displacement)).sqrt();
    let weight = if arc.constraint == ArcConstraint::Spherical { scale * scale } else { 0.0 };
    let unit = (scale * scale + weight).sqrt();
    let mut size = settings.initial_increment * unit;
    let (smallest, largest) = (settings.min_increment * unit, settings.max_increment * unit);
    let mut previous: Option<(Vec<f64>, f64)> = None;
    let mut fast = 0usize;
    let mut last: Option<Evaluation> = None;
    let max_increments = usize::try_from(settings.max_increments).unwrap_or(usize::MAX);
    loop {
        fieldcraft_model::progress::check()?;
        if history.increments.len() >= max_increments {
            warnings.push(format!(
                "the arc-length step `{}` ended at its cap of {} increments, at load factor {lambda:.6}",
                step.name, settings.max_increments
            ));
            break;
        }
        let number = history.increments.len() + 1;
        let mut iterations = Vec::new();
        match run.arc_increment(&q, lambda, &start, size, weight, previous.as_ref(), &mut iterations) {
            Ok((q1, lambda1, mut evaluation)) => {
                let change = lambda1 - lambda;
                previous = Some((q1.iter().zip(&q).map(|(a, b)| a - b).collect(), change));
                let count = iterations.len();
                (q, lambda) = (q1, lambda1);
                run.committed.clone_from(&evaluation.history);
                let u = red.expand_scaled(&q, lambda);
                keep_frame(history, lambda, &u, Some(&evaluation), warnings);
                start = run.arc_start(&mut evaluation, &q, lambda, false).map_err(fatal)?;
                let negative = start.factor.inertia().negative;
                history.increments.push(Increment {
                    number,
                    load_factor: lambda,
                    size: change,
                    iterations,
                    converged: true,
                    cutback: None,
                    negative_pivots: negative,
                    augmentations: 0,
                });
                states.push((number, lambda, negative));
                last = Some(evaluation);
                fast = if count <= FAST_ITERATIONS { fast + 1 } else { 0 };
                if fast >= 2 {
                    size *= GROWTH;
                } else if count > SLOW_ITERATIONS {
                    size *= SHRINK;
                }
                size = size.clamp(smallest, largest);
                let reached = limit.is_some_and(|(dof, value)| u.get(dof).is_some_and(|x| x.abs() >= value));
                if lambda >= arc.max_load_factor || reached {
                    break;
                }
            }
            Err(Failure::Fatal(e)) => return Err(e),
            Err(Failure::Cutback(reason)) => {
                history.increments.push(Increment {
                    number,
                    load_factor: lambda,
                    size: 0.0,
                    iterations,
                    converged: false,
                    cutback: Some(reason.clone()),
                    negative_pivots: 0,
                    augmentations: 0,
                });
                fast = 0;
                size *= CUTBACK;
                if size < smallest * (1.0 - SAME_TIME) {
                    return Err(Error::NotConverged(format!(
                        "arc-length step `{}` stopped at load factor {lambda:.6}: an increment failed ({reason}) and a shorter arc would be below the minimum increment",
                        step.name
                    )));
                }
            }
        }
    }
    history.critical_points = critical_points(&states);
    Ok((q, lambda, last))
}

/// Critical points from the converged states' (increment, load factor, negative eigenvalues):
/// where the count changes, a limit point if the load factor turns within a state of it, a
/// bifurcation otherwise.
fn critical_points(states: &[(usize, f64, usize)]) -> Vec<CriticalPoint> {
    let mut out = Vec::new();
    for k in 1..states.len() {
        let (Some(&(_, l0, n0)), Some(&(increment, _, n1))) = (states.get(k - 1), states.get(k)) else { continue };
        if n0 == n1 {
            continue;
        }
        let around: Vec<f64> = states.get(k.saturating_sub(2)..(k + 2).min(states.len())).unwrap_or_default().iter().map(|s| s.1).collect();
        let turns = around.windows(3).any(|w| (w[1] - w[0]) * (w[2] - w[1]) < 0.0);
        let (kind, load_factor) = if turns {
            // The extreme load factor next to the change.
            let rising = around.windows(2).next().is_some_and(|w| w[1] >= w[0]);
            let extreme =
                if rising { around.iter().copied().fold(f64::NEG_INFINITY, f64::max) } else { around.iter().copied().fold(f64::INFINITY, f64::min) };
            (CriticalKind::Limit, extreme)
        } else {
            (CriticalKind::Bifurcation, l0)
        };
        out.push(CriticalPoint { kind, increment, load_factor, negative_before: n0, negative_after: n1 });
    }
    out
}
