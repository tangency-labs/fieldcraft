//! The resonant modes of an electromagnetic cavity (docs/rfc/0001-electromagnetics.md §3.3):
//! the frequencies at which a field E exists with
//!
//! curl (1/μᵣ) curl E = k₀² εᵣ E,   E × n = 0 on the perfect conductors,
//!
//! k₀ = 2π f / c₀, on Nédélec edge elements (`fieldcraft_fem::nedelec`). The weak form is
//! K x = k₀² M x, K the curl–curl and M the mass matrix of the free edge (and face) unknowns.
//!
//! **The null space.** Every gradient of a potential that vanishes on the conductors solves it
//! with k₀ = 0, a seventh of the unknowns of a first-order mesh. Shift-invert Lanczos would fill
//! its basis with them, so each iterate is kept M-orthogonal to them:
//! x ← x − G (Gᵀ M G)⁻¹ Gᵀ M x, G the discrete gradient (P. Arbenz, R. Geus, *Multilevel
//! preconditioned iterative eigensolvers for Maxwell eigenvalue problems*, Appl. Numer. Math.
//! 54(2):107–121, 2005; RFC §8.4). Its columns are the potentials of the nodes off the
//! conductors, at order 2 the free edges' functions too, and one potential for each separate
//! conductor beyond the first (constant on that conductor: a floating conductor carries one).
//! A constant potential over a connected part of the mesh has no gradient, so each part not
//! touching the first conductor leaves out one potential more (its smallest node's, or a
//! conductor's). Gᵀ M G is then positive definite and factored once.
//!
//! **The shift.** Below every mode: σ = −(π/D)² · min(1/μᵣ) / max(εᵣ), D the diagonal of the
//! model's bounding box, so K − σM is positive definite (Cholesky) and the wanted modes are the
//! nearest above σ; or the `min_hz` asked for, when it is above zero (LDLᵀ). (π/D)² is the scale
//! of an empty cavity's k₀²; the materials scale each mode's Rayleigh quotient
//! ∫ ν |curl E|² / ∫ ε |E|² by at least min ν / max ε, so a filled cavity's modes and its shift
//! stay as far apart relatively, and Lanczos converges as fast as in an empty one.
//!
//! **No mode missed.** Once the modes are found, K − σ′M is factored with σ′ just below the
//! highest one: its negative eigenvalues (the LDLᵀ inertia, Sylvester's law) less the null
//! space's dimension are the modes below σ′, which a step from zero must all have returned.
//! Lanczos accepts θ = 1/(k₀² − σ) within [`TOLERANCE`] · θ of an eigenvalue, so the highest
//! k₀² found is within `TOLERANCE` · (k₀² − σ) of its discrete value, and so is each copy of a
//! degenerate one; σ′ is twice that below it, so the rest of a degenerate cluster that was not
//! asked for is not counted. With a `min_hz`, the negative eigenvalues of the shifted matrix
//! (K − σM at that bound) take the null space's place, and the difference is the modes between
//! the bound and σ′. The null space counted is G's columns: it is all of ker K when every
//! closed curve in the cavity bounds a surface in it, or the walls are all conductors. With
//! magnetic walls (no conductor on part of the boundary) on a domain with a handle, such as a
//! ring, ker K also holds harmonic fields that are no gradient; they are counted as modes below
//! σ′ that were not found, and the warning says so.
//!
//! **Fields.** Each mode's E at the nodes: the elements' fields at their corners (and a tet10's
//! mid-side nodes, halfway along its straight edges), averaged. E is
//! continuous along a conductor-free interface only tangentially, so the average smooths the
//! normal jump at a change of material (RFC §8.5; element-nodal results come later). A mode is
//! scaled so that its largest nodal |E| is 1, with its largest component positive.
//!
//! **Conductor losses.** A wall of a good conductor (a [`BoundaryCondition::Conductor`]) is held
//! as a perfect conductor, and each mode's quality factor is the energy it stores over the power
//! the walls absorb per radian (D. M. Pozar, *Microwave Engineering*, 4th ed., Wiley 2012, §6.7,
//! the perturbation method): Q = ω W / P_c, with W = (ε₀/2) ∫ εᵣ |E|² dV the stored energy and
//! P_c = (Rₛ/2) ∮ |H_t|² dS, Rₛ = √(ω μ₀ μ_c / (2σ)) (eq. 1.131) and H = j curl E / (ω μ₀ μᵣ)
//! (Faraday's law for e^{jωt}), so
//!
//! Q = ω³ ε₀ μ₀² ∫ εᵣ |E|² dV / Σ_walls Rₛ ∮ |n × curl E|² / μᵣ² dS.
//!
//! The integral over the volume is xᵀ M x; over each wall face it takes the curl of the element
//! behind it at the face's edge midpoints, exact for the quadratic |curl E|² of order 2. Its
//! error falls as h² at both orders (the trace of curl E on the walls). The walls' small
//! reactance Xₛ = Rₛ lowers each frequency by f / (2Q) as well; it is not applied.
//!
//! The elements are tetrahedra (`tet4`, and `tet10` whose edges are taken straight); every
//! material in use has electromagnetic properties ([`fieldcraft_model::Electromagnetic`]).

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_fem::dof::{DofLayout, Entity};
use fieldcraft_fem::nedelec::{Barycentric, EdgeTet, Order, Potential, declare, gradient, slots};
use fieldcraft_model::{BoundaryCondition, ElementKind, Model, Procedure, Resolved, SectionKind};

use crate::analysis::{BLOCK, Stats, TOLERANCE, factor_symmetric, lanczos_basis};
use crate::cholesky::{Factor, FactorError, UpperColumns};
use crate::eigen::{Operator, Wanted, lanczos};
use crate::ldlt::{Kind, Symmetric};
use crate::ordering::fill_reducing_order;
use crate::sparse::{Csr, NodeGraph, add_block};
use crate::{Error, Limits, Result};

/// The speed of light in vacuum, m/s (exact by the SI's definition of the metre, 2019).
pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;

/// The magnetic constant μ₀, H/m: measured since the SI of 2019, 4π × 10⁻⁷ to 1e-9.
pub const VACUUM_PERMEABILITY: f64 = 1.256_637_062_12e-6;

/// One resonant mode.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CavityMode {
    /// Frequency, Hz.
    pub frequency: f64,
    /// The electric field at each node of [`CavityResult::nodes`], scaled to a largest |E| of 1;
    /// all zero (and a warning) for a mode with no field at the nodes.
    pub field: Vec<[f64; 3]>,
    /// The quality factor of the conductors' losses, when a conductor is active. Not saved with
    /// the results yet (it would change their encoding): a model read back has none until solved.
    #[serde(skip)]
    pub quality: Option<f64>,
}

/// What a cavity step computes.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CavityResult {
    pub step: String,
    /// The node labels the fields are given at.
    pub nodes: Vec<u32>,
    /// Ascending frequency.
    pub modes: Vec<CavityMode>,
    /// The order of the edge elements.
    pub order: u8,
    /// The dimension of the null space projected out (gradients of potentials).
    pub null_space: usize,
    pub stats: Stats,
    pub warnings: Vec<String>,
}

/// One tetrahedron: label, corner labels, a tet10's mid-side nodes (each with its edge's
/// corners), 1/μᵣ and εᵣ.
struct Tet {
    label: u32,
    corners: [u32; 4],
    midsides: Vec<(u32, [u32; 2])>,
    nu: f64,
    eps: f64,
}

impl Tet {
    /// The nodes the field is given at, with their barycentric coordinates in the element on
    /// its sorted corners: the corners, and a tet10's mid-side nodes halfway along their (straight)
    /// edges.
    fn samples(&self) -> Vec<(u32, Barycentric)> {
        let mut sorted = self.corners;
        sorted.sort_unstable();
        let at = |c: u32| sorted.iter().position(|s| *s == c);
        let unit = |i: usize| -> Barycentric { std::array::from_fn(|j| if j == i { 1.0 } else { 0.0 }) };
        let corners = sorted.iter().enumerate().map(|(i, c)| (*c, unit(i)));
        let midsides = self.midsides.iter().filter_map(|(m, [a, b])| {
            let (i, j) = (at(*a)?, at(*b)?);
            Some((*m, std::array::from_fn(|k| if k == i || k == j { 0.5 } else { 0.0 })))
        });
        corners.chain(midsides).collect()
    }
}

/// The tetrahedra of every section, by label.
fn tetrahedra(model: &Model, resolved: &Resolved, warnings: &mut Vec<String>) -> Result<Vec<Tet>> {
    let mut out = Vec::new();
    let mut curved = false;
    for (name, section) in model.sections() {
        if !matches!(section.kind, SectionKind::Solid {}) {
            return Err(Error::NotSupported(format!("electromagnetics in {} sections (section `{name}`); it has tetrahedra", section.kind.name())));
        }
        let material = model
            .material(&section.material)
            .ok_or_else(|| Error::Invalid(format!("section `{name}`: material `{}` is missing", section.material)))?;
        let em = material.electromagnetic.ok_or_else(|| {
            Error::Invalid(format!("material `{}` has no electromagnetic properties (relative permittivity and permeability)", section.material))
        })?;
        let set = resolved.sets.get(&section.set).ok_or_else(|| Error::Invalid(format!("set `{}` is not resolved", section.set)))?;
        for &label in &set.elements {
            let e = model.element(label).ok_or_else(|| Error::Invalid(format!("element {label} is missing")))?;
            curved |= e.kind == ElementKind::Tet10;
            let corners = match (e.kind, e.nodes.get(..4)) {
                (ElementKind::Tet4 | ElementKind::Tet10, Some(&[a, b, c, d])) => [a, b, c, d],
                (kind, _) => {
                    return Err(Error::NotSupported(format!("electromagnetics on {} elements (element {label}); it has tetrahedra", kind.name())));
                }
            };
            let midsides = match e.kind {
                ElementKind::Tet10 => ElementKind::Tet10
                    .edges()
                    .iter()
                    .filter_map(|edge| match **edge {
                        [i, j, m] => Some((*e.nodes.get(m)?, [*e.nodes.get(i)?, *e.nodes.get(j)?])),
                        _ => None,
                    })
                    .collect(),
                _ => Vec::new(),
            };
            out.push(Tet { label, corners, midsides, nu: 1.0 / em.permeability, eps: em.permittivity });
        }
    }
    if curved {
        warnings.push("tet10 elements are taken with straight edges (curved edge elements are not supported yet)".into());
    }
    out.sort_by_key(|t| t.label);
    Ok(out)
}

/// The edges and faces the active conductors hold, the nodes on them grouped into connected
/// conductors (by the held edges), and the faces that absorb power.
struct Walls {
    edges: BTreeSet<[u32; 2]>,
    faces: BTreeSet<[u32; 3]>,
    conductors: Vec<Vec<u32>>,
    lossy: Vec<LossyFace>,
}

/// A face of a good conductor: the tetrahedron behind it, its corners and the metal's σ and μᵣ.
struct LossyFace {
    element: u32,
    nodes: [u32; 3],
    conductivity: f64,
    permeability: f64,
}

fn walls(model: &Model, resolved: &Resolved, bcs: &BTreeSet<String>, corners: &BTreeMap<u32, [u32; 4]>) -> Result<Walls> {
    let (mut edges, mut faces, mut lossy) = (BTreeSet::new(), BTreeSet::new(), Vec::new());
    for name in bcs {
        let (target, metal) = match model.boundary_condition(name) {
            Some(BoundaryCondition::PerfectConductor { target }) => (target, None),
            Some(BoundaryCondition::Conductor { target, conductivity, permeability }) => (target, Some((*conductivity, *permeability))),
            _ => continue,
        };
        let set = target
            .set()
            .and_then(|s| resolved.sets.get(s))
            .ok_or_else(|| Error::Invalid(format!("boundary condition `{name}`: its set is not resolved")))?;
        for f in &set.faces {
            let nodes = model
                .face_nodes(*f)
                .ok_or_else(|| Error::Invalid(format!("boundary condition `{name}`: element {} has no face {}", f.element, f.face)))?;
            // A face left out would be a magnetic wall, and the resonances wrong with no word.
            let c = corners.get(&f.element).ok_or_else(|| {
                Error::Invalid(format!("boundary condition `{name}`: element {} face {} is on no tetrahedron of a section", f.element, f.face))
            })?;
            let mut key: Vec<u32> = nodes.into_iter().filter(|n| c.contains(n)).collect();
            key.sort_unstable();
            let [a, b, d] = key[..] else {
                return Err(Error::Invalid(format!("boundary condition `{name}`: element {} face {} is not a triangle", f.element, f.face)));
            };
            faces.insert([a, b, d]);
            edges.extend([[a, b], [a, d], [b, d]]);
            if let Some((conductivity, permeability)) = metal {
                lossy.push(LossyFace { element: f.element, nodes: [a, b, d], conductivity, permeability });
            }
        }
    }
    // Conductors: the nodes of the held edges, joined by them.
    let conductors = connected(edges.iter().map(|[a, b]| (*a, *b))).into_values().collect();
    Ok(Walls { edges, faces, conductors, lossy })
}

/// The connected parts of a graph given by its links, by their smallest member, each part's
/// members ascending (union–find, the smaller root kept).
fn connected<T: Ord + Copy>(links: impl IntoIterator<Item = (T, T)>) -> BTreeMap<T, Vec<T>> {
    fn find<T: Ord + Copy>(root: &mut BTreeMap<T, T>, n: T) -> T {
        let mut r = n;
        while let Some(&p) = root.get(&r).filter(|p| **p != r) {
            r = p;
        }
        root.insert(n, r);
        r
    }
    let mut root = BTreeMap::new();
    for (a, b) in links {
        let (ra, rb) = (find(&mut root, a), find(&mut root, b));
        root.insert(ra.max(rb), ra.min(rb));
    }
    let mut parts: BTreeMap<T, Vec<T>> = BTreeMap::new();
    for n in root.keys().copied().collect::<Vec<_>>() {
        let r = find(&mut root, n);
        parts.entry(r).or_default().push(n);
    }
    parts
}

/// A column of the discrete gradient: a potential whose gradient is a null-space field.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Column {
    /// A node off the conductors.
    Node(u32),
    /// A conductor beyond the first, by its place in [`Walls::conductors`].
    Conductor(usize),
    /// A free edge's quadratic function (order 2), by its layout slot.
    Edge(usize),
}

/// A sparse matrix with `cols` columns, by rows.
struct Rows {
    cols: usize,
    ptr: Vec<usize>,
    idx: Vec<usize>,
    val: Vec<f64>,
}

impl Rows {
    fn row(&self, r: usize) -> (&[usize], &[f64]) {
        let (a, b) = (self.ptr.get(r).copied().unwrap_or(0), self.ptr.get(r + 1).copied().unwrap_or(0));
        (self.idx.get(a..b).unwrap_or_default(), self.val.get(a..b).unwrap_or_default())
    }

    fn rows(&self) -> usize {
        self.ptr.len().saturating_sub(1)
    }

    /// Aᵀ.
    fn transpose(&self) -> Rows {
        let mut count = vec![0usize; self.cols + 1];
        for c in &self.idx {
            if let Some(x) = count.get_mut(c + 1) {
                *x += 1;
            }
        }
        for c in 0..self.cols {
            let before = count.get(c).copied().unwrap_or(0);
            if let Some(x) = count.get_mut(c + 1) {
                *x += before;
            }
        }
        let mut next = count.clone();
        let (mut idx, mut val) = (vec![0; self.idx.len()], vec![0.0; self.val.len()]);
        for r in 0..self.rows() {
            let (cols, vals) = self.row(r);
            for (c, v) in cols.iter().zip(vals) {
                if let Some(at) = next.get_mut(*c) {
                    if let (Some(i), Some(x)) = (idx.get_mut(*at), val.get_mut(*at)) {
                        *i = r;
                        *x = *v;
                    }
                    *at += 1;
                }
            }
        }
        Rows { cols: self.rows(), ptr: count, idx, val }
    }

    /// A B, rows in order and each row's sums in the order of A's columns (Gustavson).
    fn times(&self, b: &Rows) -> Rows {
        let (mut ptr, mut idx, mut val) = (vec![0], Vec::new(), Vec::new());
        let mut acc = vec![0.0; b.cols];
        let mut seen = vec![false; b.cols];
        let mut touched = Vec::new();
        for r in 0..self.rows() {
            let (cols, vals) = self.row(r);
            for (k, a) in cols.iter().zip(vals) {
                let (bc, bv) = b.row(*k);
                for (c, v) in bc.iter().zip(bv) {
                    if let (Some(s), Some(x)) = (seen.get_mut(*c), acc.get_mut(*c)) {
                        if !*s {
                            *s = true;
                            touched.push(*c);
                        }
                        *x += a * v;
                    }
                }
            }
            touched.sort_unstable();
            for c in touched.drain(..) {
                if let (Some(s), Some(x)) = (seen.get_mut(c), acc.get_mut(c)) {
                    idx.push(c);
                    val.push(*x);
                    *s = false;
                    *x = 0.0;
                }
            }
            ptr.push(idx.len());
        }
        Rows { cols: b.cols, ptr, idx, val }
    }

    fn of(a: &Csr) -> Rows {
        Rows { cols: a.n, ptr: a.ptr.clone(), idx: a.idx.clone(), val: a.val.clone() }
    }

    /// y = A x.
    fn mul(&self, x: &[f64]) -> Vec<f64> {
        (0..self.rows())
            .map(|r| {
                let (cols, vals) = self.row(r);
                cols.iter().zip(vals).map(|(c, v)| v * x.get(*c).copied().unwrap_or(0.0)).sum()
            })
            .collect()
    }
}

/// The rows and columns `keep` of `a`, in that order.
fn restrict(a: &Csr, keep: &[usize]) -> Csr {
    let mut new = vec![usize::MAX; a.n];
    for (k, d) in keep.iter().enumerate() {
        if let Some(x) = new.get_mut(*d) {
            *x = k;
        }
    }
    let (mut ptr, mut idx, mut val) = (vec![0], Vec::new(), Vec::new());
    for d in keep {
        let (cols, vals) = a.row(*d);
        let mut row: Vec<(usize, f64)> =
            cols.iter().zip(vals).filter_map(|(c, v)| new.get(*c).copied().filter(|n| *n != usize::MAX).map(|n| (n, *v))).collect();
        row.sort_by_key(|e| e.0);
        idx.extend(row.iter().map(|e| e.0));
        val.extend(row.iter().map(|e| e.1));
        ptr.push(idx.len());
    }
    Csr { n: keep.len(), ptr, idx, val }
}

/// A positive definite matrix factored in a fill-reducing order of its own.
struct Ordered {
    order: Vec<usize>,
    factor: Factor,
}

impl Ordered {
    fn new(a: &Csr, points: &[[f64; 3]], limits: &Limits) -> Result<Ordered> {
        let adjacency: Vec<Vec<usize>> = (0..a.n).map(|r| a.row(r).0.to_vec()).collect();
        let order = fill_reducing_order(points, &adjacency);
        let factor = Factor::new(&UpperColumns::from_symmetric(&restrict(a, &order)), limits.max_factor_entries).map_err(factor_error)?;
        Ok(Ordered { order, factor })
    }

    fn solve(&self, b: &[f64]) -> Result<Vec<f64>> {
        let mut x: Vec<f64> = self.order.iter().map(|i| b.get(*i).copied().unwrap_or(0.0)).collect();
        self.factor.solve(&mut x)?;
        let mut out = vec![0.0; b.len()];
        for (i, v) in self.order.iter().zip(x) {
            if let Some(o) = out.get_mut(*i) {
                *o = v;
            }
        }
        Ok(out)
    }
}

fn factor_error(e: FactorError) -> Error {
    match e {
        FactorError::TooLarge { entries, limit } => Error::factor_limit(entries, limit),
        FactorError::OutOfMemory { entries } => Error::Limit(format!("not enough memory for the factorization ({entries} entries)")),
        FactorError::Cancelled => Error::Model(fieldcraft_model::Error::Cancelled),
        FactorError::Internal(what) => Error::Internal(what),
        FactorError::NotPositiveDefinite { .. } | FactorError::Singular { .. } => {
            Error::Singular("the cavity's matrices are singular: a region no element reaches, or elements that are flat to working precision".into())
        }
    }
}

/// Shift-invert with the null space projected out.
struct Projected<'a> {
    m: &'a Csr,
    shifted: Symmetric,
    g: Rows,
    gt: Rows,
    gmg: Ordered,
}

impl Projected<'_> {
    /// x ← x − G (Gᵀ M G)⁻¹ Gᵀ M x.
    fn project(&self, x: &mut [f64]) -> crate::Result<()> {
        let mut mx = vec![0.0; x.len()];
        self.m.mul(x, &mut mx);
        let z = self.gmg.solve(&self.gt.mul(&mx))?;
        for (xi, gz) in x.iter_mut().zip(self.g.mul(&z)) {
            *xi -= gz;
        }
        Ok(())
    }
}

impl Operator for Projected<'_> {
    fn n(&self) -> usize {
        self.m.n
    }

    fn mass(&self, x: &[f64], y: &mut [f64]) {
        self.m.mul(x, y);
    }

    fn solve(&self, x: &mut [f64]) -> crate::Result<()> {
        self.shifted.solve(x)?;
        self.project(x)
    }

    /// A block's solves together (the factor read once), then each projected.
    fn solve_many(&self, xs: &mut [Vec<f64>]) -> crate::Result<()> {
        let mut views: Vec<&mut [f64]> = xs.iter_mut().map(Vec::as_mut_slice).collect();
        self.shifted.solve_many(&mut views)?;
        xs.iter_mut().try_for_each(|x| self.project(x))
    }
}

/// The resonant modes of the model's cavity step (`step`, or its only one).
pub fn solve_cavity(model: &Model, step: Option<&str>, limits: &Limits) -> Result<CavityResult> {
    fieldcraft_model::progress::check()?;
    let resolved = model.resolve()?;
    let step = model.select_step(step, |p| matches!(p, Procedure::Cavity { .. }), "cavity")?;
    let Procedure::Cavity { modes, min_hz, max_hz, order } = step.procedure else {
        return Err(Error::Invalid(format!("step `{}` is a {} step, not a cavity", step.name, step.procedure.name())));
    };
    let degree = order;
    let order = match degree {
        1 => Order::One,
        2 => Order::Two,
        d => return Err(Error::Invalid(format!("edge elements are of order 1 or 2, got {d}"))),
    };
    let mut warnings: Vec<String> = resolved.warnings.iter().map(ToString::to_string).collect();
    let tets = tetrahedra(model, &resolved, &mut warnings)?;
    if tets.is_empty() {
        return Err(Error::Invalid("the model has no tetrahedra for the field".into()));
    }
    let corners: BTreeMap<u32, [u32; 4]> = tets.iter().map(|t| (t.label, t.corners)).collect();
    let walls = walls(model, &resolved, &step.boundary_conditions, &corners)?;
    if walls.faces.is_empty() {
        warnings.push("no conductor is active: every boundary is a magnetic wall (n × H = 0)".into());
    }

    // The layout: every tetrahedron's edges (and faces at order 2).
    let mut builder = DofLayout::builder();
    for t in &tets {
        declare(&mut builder, t.corners, order).map_err(|e| Error::Element { label: t.label, source: e })?;
    }
    let layout = builder.build();
    if layout.dofs() > limits.max_dofs {
        return Err(Error::Limit(format!("the cavity has {} unknowns; the limit is {}", layout.dofs(), limits.max_dofs)));
    }
    let xyz = |l: u32| model.nodes().get(&l).map(|n| n.xyz).unwrap_or_default();
    let element_slots: Vec<Vec<usize>> = tets
        .iter()
        .map(|t| slots(&layout, t.corners, order).ok_or_else(|| Error::Internal(format!("element {}: an edge the layout lacks", t.label))))
        .collect::<Result<_>>()?;
    let graph = NodeGraph::new(layout.slots(), element_slots.iter());
    let (mut k, mut m) = (graph.zero_matrix(&layout), graph.zero_matrix(&layout));
    let element =
        |t: &Tet| EdgeTet::of_labels(order, t.corners.map(xyz), t.corners).map(|(e, _)| e).map_err(|e| Error::Element { label: t.label, source: e });
    for (t, s) in tets.iter().zip(&element_slots) {
        fieldcraft_model::progress::check()?;
        let e = element(t)?;
        let element_error = |source| Error::Element { label: t.label, source };
        add_block(&mut k, &layout, s, order.kinds(), &e.curl_curl([t.nu; 3]).map_err(element_error)?.data)?;
        add_block(&mut m, &layout, s, order.kinds(), &e.mass([t.eps; 3]).map_err(element_error)?.data)?;
    }

    // The free unknowns, in a fill-reducing order of the slots by their points.
    let centre = |ls: &[u32]| {
        let w = 1.0 / ls.len() as f64;
        ls.iter().fold([0.0; 3], |p, l| [0, 1, 2].map(|a| p[a] + w * xyz(*l)[a]))
    };
    let points: Vec<[f64; 3]> = (0..layout.slots())
        .map(|s| match layout.entity(s) {
            Some(Entity::Edge(e)) => centre(&e),
            Some(Entity::Face(f)) => centre(&f),
            Some(Entity::Node(l) | Entity::Point(l)) => xyz(l),
            None => [0.0; 3],
        })
        .collect();
    let held = |s: usize| match layout.entity(s) {
        Some(Entity::Edge(e)) => walls.edges.contains(&e),
        Some(Entity::Face(f)) => walls.faces.contains(&f),
        _ => true,
    };
    let free: Vec<usize> = fill_reducing_order(&points, &graph.adjacency).into_iter().filter(|s| !held(*s)).flat_map(|s| layout.range(s)).collect();
    if free.is_empty() {
        return Err(Error::Invalid("every edge lies on a perfect conductor; there is no field to resonate".into()));
    }
    // One full matrix at a time beside the restricted ones.
    let kf = restrict(&k, &free);
    drop(k);
    let mf = restrict(&m, &free);
    drop(m);

    // The discrete gradient on the free unknowns, a column per potential: the free nodes, a
    // conductor beyond the first (the sum of its nodes' potentials), the free edges (order 2).
    let mut row_of = vec![None; layout.dofs()];
    for (r, d) in free.iter().enumerate() {
        if let Some(x) = row_of.get_mut(*d) {
            *x = Some(r);
        }
    }
    let row_of = |d: usize| row_of.get(d).copied().flatten();
    let conductor_of: BTreeMap<u32, usize> = walls.conductors.iter().enumerate().flat_map(|(c, ns)| ns.iter().map(move |n| (*n, c))).collect();
    let all_nodes: BTreeSet<u32> = tets.iter().flat_map(|t| t.corners).collect();
    // The nodal potentials (a conductor's nodes as one) joined by the free edges: a connected
    // part's constant potential has no gradient on the free unknowns, so each part gives up one
    // potential, the first conductor's (zero) where it holds it, else its smallest.
    let vertex = |n: u32| conductor_of.get(&n).map_or(Column::Node(n), |c| Column::Conductor(*c));
    let free_edges = (0..layout.slots()).filter(|s| !held(*s)).filter_map(|s| match layout.entity(s) {
        Some(Entity::Edge([a, b])) => Some((vertex(a), vertex(b))),
        _ => None,
    });
    let pinned: BTreeSet<Column> = connected(free_edges)
        .into_values()
        .filter_map(|part| if part.contains(&Column::Conductor(0)) { Some(Column::Conductor(0)) } else { part.first().copied() })
        .chain(walls.conductors.first().map(|_| Column::Conductor(0)))
        .collect();
    let column = |p: Potential| -> Option<Column> {
        match p {
            Potential::Node(n) => Some(vertex(n)).filter(|v| !pinned.contains(v)),
            Potential::Edge(s) => (!held(s)).then_some(Column::Edge(s)),
        }
    };
    // Summed by (row, column): two nodes of one conductor fall into one column.
    let mut entries: BTreeMap<(usize, Column), f64> = BTreeMap::new();
    for (d, p, v) in gradient(&layout) {
        if let (Some(r), Some(c)) = (row_of(d), column(p)) {
            *entries.entry((r, c)).or_default() += v;
        }
    }
    let potentials: Vec<Column> = entries.keys().map(|k| k.1).collect::<BTreeSet<_>>().into_iter().collect();
    let index: BTreeMap<Column, usize> = potentials.iter().enumerate().map(|(i, p)| (*p, i)).collect();
    let mut g = Rows { cols: potentials.len(), ptr: vec![0], idx: Vec::with_capacity(entries.len()), val: Vec::with_capacity(entries.len()) };
    let mut by_row: Vec<Vec<(usize, f64)>> = vec![Vec::new(); free.len()];
    for ((r, c), v) in entries {
        if let (Some(row), Some(i)) = (by_row.get_mut(r), index.get(&c)) {
            row.push((*i, v));
        }
    }
    for mut row in by_row {
        row.sort_by_key(|e| e.0);
        g.idx.extend(row.iter().map(|e| e.0));
        g.val.extend(row.iter().map(|e| e.1));
        g.ptr.push(g.idx.len());
    }
    let null_space = potentials.len();
    let gt = g.transpose();
    let gmg = Rows::of(&mf).times(&g);
    let gmg = gt.times(&gmg);
    let gmg = Csr { n: gmg.rows(), ptr: gmg.ptr, idx: gmg.idx, val: gmg.val };
    let potential_points: Vec<[f64; 3]> = potentials
        .iter()
        .map(|p| match *p {
            Column::Node(n) => xyz(n),
            Column::Edge(s) => points.get(s).copied().unwrap_or_default(),
            Column::Conductor(c) => walls.conductors.get(c).and_then(|ns| ns.first()).map(|n| xyz(*n)).unwrap_or_default(),
        })
        .collect();
    let gmg = Ordered::new(&gmg, &potential_points, limits)?;

    // The shift, and the modes.
    let (lo, hi) = all_nodes.iter().copied().fold(([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]), |(lo, hi), n| {
        let p = xyz(n);
        ([0, 1, 2].map(|a| lo[a].min(p[a])), [0, 1, 2].map(|a| hi[a].max(p[a])))
    });
    let diagonal2: f64 = (0..3).map(|a| (hi[a] - lo[a]) * (hi[a] - lo[a])).sum();
    let k2 = |hz: f64| {
        let k = 2.0 * std::f64::consts::PI * hz / SPEED_OF_LIGHT;
        k * k
    };
    let hz = |k2: f64| k2.sqrt() * SPEED_OF_LIGHT / (2.0 * std::f64::consts::PI);
    // A lower bound of zero is no bound: σ = 0 would leave K − σM singular on the gradients.
    let lower = min_hz.filter(|f| *f > 0.0);
    let (nu_min, eps_max) = tets.iter().fold((f64::INFINITY, 0.0f64), |(n, e), t| (n.min(t.nu), e.max(t.eps)));
    let sigma = lower.map_or(-std::f64::consts::PI * std::f64::consts::PI / diagonal2 * nu_min / eps_max, k2);
    if !sigma.is_finite() {
        return Err(Error::Invalid(format!("the shift {sigma} of the cavity's eigenproblem is not finite")));
    }
    let shifted = factor_symmetric(&kf.add_scaled(-sigma, &mf)?, limits, Kind::Auto).map_err(factor_error)?;
    // Below the lower bound: the null space and every mode under it (the check's baseline).
    let below_lower = lower.map(|_| shifted.inertia().negative);
    let factor_entries = shifted.entries();
    let n = mf.n;
    let wanted =
        Wanted { count: usize::try_from(modes).unwrap_or(usize::MAX).min(n), min: min_hz.map_or(0.0, k2), max: max_hz.map_or(f64::INFINITY, k2) };
    let op = Projected { m: &mf, shifted, g, gt, gmg };
    let pairs = lanczos(&op, sigma, wanted, BLOCK, lanczos_basis(limits, n), TOLERANCE)?;
    // The modes are positive; a value at zero to the eigensolver's tolerance, relative to the
    // largest found, is a gradient the projection let through.
    let largest = pairs.values.iter().copied().fold(0.0, f64::max);
    let found: Vec<(f64, &Vec<f64>)> = pairs.values.iter().copied().zip(&pairs.vectors).filter(|(v, _)| *v > TOLERANCE * largest).collect();
    // Its factor is not needed beside the check's.
    drop(op);
    if found.is_empty() {
        let range = match (min_hz, max_hz) {
            (None, None) => String::new(),
            (lo, hi) => format!(" between {} and {} Hz", lo.unwrap_or(0.0), hi.map_or("∞".into(), |h| h.to_string())),
        };
        return Err(Error::Invalid(format!("no resonant mode was found{range}")));
    }

    // No mode between the lower bound (or zero) and the highest one missed: the inertia of K − σ′M.
    if let Some(s) = found.last().map(|(top, _)| top - 2.0 * TOLERANCE * (top - sigma)).filter(|s| *s > sigma) {
        let check = factor_symmetric(&kf.add_scaled(-s, &mf)?, limits, Kind::Indefinite).map_err(factor_error)?;
        let below = check.inertia().negative.saturating_sub(below_lower.unwrap_or(null_space));
        let found_below = found.iter().filter(|(v, _)| *v < s).count();
        if below != found_below {
            let from = lower.map_or(String::new(), |l| format!(" between {l} Hz and"));
            let to = if lower.is_some() { "" } else { " below" };
            // A boundary face of one tetrahedron that no conductor holds is a magnetic wall.
            let mut faces: BTreeMap<[u32; 3], u32> = BTreeMap::new();
            for [a, b, c, d] in tets.iter().map(|t| t.corners) {
                for mut f in [[a, b, c], [a, b, d], [a, c, d], [b, c, d]] {
                    f.sort_unstable();
                    *faces.entry(f).or_default() += 1;
                }
            }
            let magnetic = faces.iter().any(|(f, n)| *n == 1 && !walls.faces.contains(f));
            let why = if magnetic && below > found_below {
                "; with magnetic walls on a cavity with a handle (a ring), fields that are no gradient also lie at zero and are counted"
            } else {
                ""
            };
            warnings.push(format!("{below} modes lie{from}{to} {:.6e} Hz, and {found_below} were found{why}", hz(s)));
        }
    }

    // Fields at the nodes.
    let samples: Vec<Vec<(u32, Barycentric)>> = tets.iter().map(Tet::samples).collect();
    let nodes: Vec<u32> = samples.iter().flatten().map(|(n, _)| *n).collect::<BTreeSet<_>>().into_iter().collect();
    let at: BTreeMap<u32, usize> = nodes.iter().enumerate().map(|(i, n)| (*n, i)).collect();
    let mut count = vec![0u32; nodes.len()];
    for (n, _) in samples.iter().flatten() {
        if let Some(x) = at.get(n).and_then(|i| count.get_mut(*i)) {
            *x += 1;
        }
    }
    // Element by element: its nodes' field in every mode at once (the sums at each node in
    // element order, as one mode at a time would take them).
    let mut fields = vec![vec![[0.0; 3]; nodes.len()]; found.len()];
    for ((t, s), points) in tets.iter().zip(&element_slots).zip(&samples) {
        fieldcraft_model::progress::check()?;
        let e = element(t)?;
        let dofs = layout.element_dofs(s, order.kinds()).unwrap_or_default();
        let rows: Vec<Option<usize>> = dofs.iter().map(|d| row_of(*d)).collect();
        for (c, l) in points {
            let values = e.values(*l);
            let Some(node) = at.get(c).copied() else { continue };
            for (field, (_, x)) in fields.iter_mut().zip(&found) {
                let v = rows.iter().zip(&values).fold([0.0; 3], |a, (r, v)| {
                    let c = r.and_then(|r| x.get(r)).copied().unwrap_or(0.0);
                    [0, 1, 2].map(|k| a[k] + c * v[k])
                });
                if let Some(f) = field.get_mut(node) {
                    *f = [0, 1, 2].map(|k| f[k] + v[k]);
                }
            }
        }
    }
    let quality = conductor_quality(&walls.lossy, &tets, &element_slots, &layout, order, &mf, &row_of, &element, &found)?;
    let result_modes: Vec<CavityMode> = fields
        .into_iter()
        .zip(&found)
        .zip(quality)
        .map(|((field, (value, _)), quality)| {
            let field: Vec<[f64; 3]> = field.into_iter().zip(&count).map(|(f, c)| f.map(|v| v / f64::from((*c).max(1)))).collect();
            // Largest |E| 1, its largest component positive (the first such, in node and axis order).
            let size = |f: &[f64; 3]| f.iter().map(|v| v * v).sum::<f64>();
            let big = field.iter().fold(None::<&[f64; 3]>, |best, f| if best.is_none_or(|b| size(f) > size(b)) { Some(f) } else { best });
            let scale = big.filter(|b| size(b) > 0.0).map_or(1.0, |b| {
                let c = b.iter().copied().fold(0.0f64, |a, v| if v.abs() > a.abs() { v } else { a });
                c.signum() / size(b).sqrt()
            });
            CavityMode { frequency: hz(*value), field: field.into_iter().map(|f| f.map(|v| v * scale)).collect(), quality }
        })
        .collect();
    for (k, mode) in result_modes.iter().enumerate() {
        if mode.field.iter().all(|f| f.iter().all(|v| *v == 0.0)) {
            warnings.push(format!("mode {} has no field at the nodes (it lives on the faces' functions); its field is shown as zero", k + 1));
        }
    }
    // Against the count asked for, not the one held to the free degrees of freedom.
    if result_modes.len() < usize::try_from(modes).unwrap_or(usize::MAX) {
        warnings.push(format!("{} of the {} modes asked for were found", result_modes.len(), modes));
    }
    fieldcraft_model::progress::check()?;
    let stats = Stats { dofs: layout.dofs(), equations: n, factor_entries, lanczos_vectors: pairs.basis };
    Ok(CavityResult { step: step.name.clone(), nodes, modes: result_modes, order: degree, null_space, stats, warnings })
}

/// Each mode's quality factor from the conductors' losses (the module's **Conductor losses**);
/// `None` for every mode when no conductor absorbs power.
#[allow(clippy::too_many_arguments)]
fn conductor_quality(
    lossy: &[LossyFace],
    tets: &[Tet],
    element_slots: &[Vec<usize>],
    layout: &DofLayout,
    order: Order,
    mass: &Csr,
    row_of: &dyn Fn(usize) -> Option<usize>,
    element: &dyn Fn(&Tet) -> Result<EdgeTet>,
    found: &[(f64, &Vec<f64>)],
) -> Result<Vec<Option<f64>>> {
    if lossy.is_empty() {
        return Ok(vec![None; found.len()]);
    }
    let index: BTreeMap<u32, usize> = tets.iter().enumerate().map(|(i, t)| (t.label, i)).collect();
    let omega = |k2: f64| k2.sqrt() * SPEED_OF_LIGHT;
    // Σ Rₛ ∮ |n × curl E|² / μᵣ² dS per mode, face by face in the walls' order.
    let mut loss = vec![0.0; found.len()];
    for face in lossy {
        fieldcraft_model::progress::check()?;
        let i = index
            .get(&face.element)
            .copied()
            .ok_or_else(|| Error::Internal(format!("conductor face on element {}, which is no tetrahedron", face.element)))?;
        let (Some(t), Some(s)) = (tets.get(i), element_slots.get(i)) else {
            return Err(Error::Internal(format!("element {}: no slots", face.element)));
        };
        let e = element(t)?;
        let rows: Vec<Option<usize>> = layout.element_dofs(s, order.kinds()).unwrap_or_default().into_iter().map(row_of).collect();
        // The face's corners in the element on its sorted corners; the opposite corner's coordinate is 0.
        let mut sorted = t.corners;
        sorted.sort_unstable();
        let at = |n: u32| sorted.iter().position(|c| *c == n);
        let [Some(a), Some(b), Some(c)] = face.nodes.map(at) else {
            return Err(Error::Internal(format!("element {}: a conductor face is not one of its faces", face.element)));
        };
        let corners = e.corners();
        let (Some(pa), Some(pb), Some(pc)) = (corners.get(a).copied(), corners.get(b).copied(), corners.get(c).copied()) else {
            return Err(Error::Internal(format!("element {}: a corner out of range", face.element)));
        };
        let (u, v) = ([0, 1, 2].map(|k| pb[k] - pa[k]), [0, 1, 2].map(|k| pc[k] - pa[k]));
        let cross = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        let twice_area = cross.iter().map(|x| x * x).sum::<f64>().sqrt();
        if twice_area.is_nan() || twice_area <= 0.0 {
            return Err(Error::Invalid(format!("element {}: a conductor face has no area", face.element)));
        }
        let n = cross.map(|x| x / twice_area);
        // The edge midpoints, each a third of the area: exact for a quadratic integrand.
        let points: [Barycentric; 3] = [(a, b), (b, c), (c, a)].map(|(i, j)| std::array::from_fn(|k| if k == i || k == j { 0.5 } else { 0.0 }));
        let curls: Vec<Vec<[f64; 3]>> = points.iter().map(|l| e.curl_values(*l)).collect();
        let weight = twice_area / 6.0 * t.nu * t.nu;
        for (sum, (k2, x)) in loss.iter_mut().zip(found) {
            let w = omega(*k2);
            let rs = (w * VACUUM_PERMEABILITY * face.permeability / (2.0 * face.conductivity)).sqrt();
            let integral: f64 = curls
                .iter()
                .map(|curl| {
                    let h = rows.iter().zip(curl).fold([0.0; 3], |h, (r, c)| {
                        let x = r.and_then(|r| x.get(r)).copied().unwrap_or(0.0);
                        [0, 1, 2].map(|k| h[k] + x * c[k])
                    });
                    let t = [n[1] * h[2] - n[2] * h[1], n[2] * h[0] - n[0] * h[2], n[0] * h[1] - n[1] * h[0]];
                    t.iter().map(|v| v * v).sum::<f64>()
                })
                .sum();
            *sum += rs * weight * integral;
        }
    }
    let mut mx = vec![0.0; mass.n];
    Ok(found
        .iter()
        .zip(loss)
        .map(|((k2, x), loss)| {
            mass.mul(x, &mut mx);
            let energy: f64 = x.iter().zip(&mx).map(|(a, b)| a * b).sum();
            let w = omega(*k2);
            // ω³ ε₀ μ₀² = ω³ μ₀ / c².
            let q = w * w * w * VACUUM_PERMEABILITY / (SPEED_OF_LIGHT * SPEED_OF_LIGHT) * energy / loss;
            (q.is_finite() && q > 0.0).then_some(q)
        })
        .collect())
}
