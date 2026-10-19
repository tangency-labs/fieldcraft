//! Heat transfer: steady and transient conduction with prescribed temperatures, surface fluxes,
//! convection and internal heat sources.
//!
//! Every node of the solid elements carries one temperature ([`fieldcraft_fem::heat::HEAT_DOFS`]).
//! The conductivity K (with the convection films) and capacity C are assembled on the node
//! graph in element label order ([`fieldcraft_fem::heat`]); prescribed temperatures are
//! eliminated with the ties (each secondary node's temperature follows the main nodes' with the
//! displacement tie's mortar weights, [`crate::tie`]), the free ones numbered by nested
//! dissection and solved by sparse Cholesky.
//!
//! - **Steady**: K T = f.
//! - **Transient**, the generalized trapezoidal (θ) method with n equal increments Δt over the
//!   step (T. J. R. Hughes, *The Finite Element Method*, Dover 2000, §8.2; unconditionally stable
//!   for θ ≥ ½, θ = 1 backward Euler, θ = ½ Crank–Nicolson):
//!   (C + θΔt K) Tₙ₊₁ = (C − (1 − θ)Δt K) Tₙ + Δt f, from a uniform initial temperature, with
//!   the prescribed temperatures held from the first increment on. The matrix is factored once.
//!
//! The heat flux −k∇T is reported at the integration points and averaged at the nodes; the heat
//! flow at prescribed temperatures is the reaction K T − f (C Ṫ included in transient steps).

use fieldcraft_fem::element::{Extrapolation, extrapolation};
use fieldcraft_fem::frame::norm;
use fieldcraft_fem::heat::{Conduction, PointFlux, face_film, face_weights};
use fieldcraft_fem::mortar::area_vector;
use fieldcraft_fem::{DofKind, DofLayout, DofSet};
use fieldcraft_model::{ElementKind, Load, Model, Procedure, Resolved, SectionKind, Step, TEMPERATURE_DOF};

use crate::cholesky::{Factor, FactorError, UpperColumns};
use crate::constraints::{Budget, Reduction, Row};
use crate::contact::thickness;
use crate::ordering::fill_reducing_order;
use crate::sparse::{Csr, NodeGraph, add_block};
use crate::system::hold;
use crate::tie::{TieFace, TieNodes, rows_of};
use crate::{Error, Limits, Result, Stats};

/// Most temperature values a transient result keeps (frames × nodes); frames are thinned
/// evenly to stay below it, the last increment always kept.
pub const MAX_TEMPERATURE_VALUES: usize = 50_000_000;

/// A solid element conducting heat.
#[derive(Debug, Clone)]
pub struct HeatElem {
    pub label: u32,
    pub kind: ElementKind,
    /// Indices into [`HeatMesh::labels`].
    pub nodes: Vec<usize>,
    pub xyz: Vec<[f64; 3]>,
    /// Conductivity along X, Y, Z.
    pub conductivity: [f64; 3],
    /// ρc, when the material has a density and a specific heat.
    pub capacity: Option<f64>,
}

impl HeatElem {
    fn kernel(&self) -> Result<Conduction<'_>> {
        Conduction::new(self.kind, &self.xyz).map_err(|e| Error::element(self.label, e))
    }
}

/// The conducting mesh: the nodes of the solid elements by label, and the elements by label.
#[derive(Debug, Clone)]
pub struct HeatMesh {
    pub labels: Vec<u32>,
    pub xyz: Vec<[f64; 3]>,
    pub elements: Vec<HeatElem>,
    pub layout: DofLayout,
}

impl HeatMesh {
    pub fn new(model: &Model, resolved: &Resolved, max_dofs: usize) -> Result<HeatMesh> {
        let mut elements = Vec::new();
        for (name, section) in model.sections() {
            if !matches!(section.kind, SectionKind::Solid {}) {
                return Err(Error::NotSupported(format!(
                    "heat transfer in {} sections (section `{name}`); it has solid elements so far",
                    section.kind.name()
                )));
            }
            let material = model
                .material(&section.material)
                .ok_or_else(|| Error::Invalid(format!("section `{name}`: material `{}` is missing", section.material)))?;
            let thermal = material.thermal.unwrap_or_default();
            let conductivity =
                thermal.conductivity.ok_or_else(|| Error::Invalid(format!("material `{}` has no conductivity", section.material)))?.axes();
            let capacity = material.density.zip(thermal.specific_heat).map(|(rho, c)| rho * c);
            let set = resolved.sets.get(&section.set).ok_or_else(|| Error::Invalid(format!("set `{}` is not resolved", section.set)))?;
            for &label in &set.elements {
                let el = model.element(label).ok_or_else(|| Error::Invalid(format!("element {label} is missing")))?;
                elements.push((label, el.kind, el.nodes.clone(), conductivity, capacity));
            }
        }
        elements.sort_by_key(|e| e.0);
        let layout = DofLayout::uniform(elements.iter().flat_map(|e| e.2.iter().copied()), DofSet::TEMPERATURE);
        if layout.dofs() > max_dofs {
            return Err(Error::Limit(format!("the model has {} temperatures; the limit is {max_dofs}", layout.dofs())));
        }
        let labels: Vec<u32> = (0..layout.nodes()).filter_map(|s| layout.label(s)).collect();
        let xyz: Vec<[f64; 3]> = labels
            .iter()
            .map(|l| model.node(*l).map(|n| n.xyz).ok_or_else(|| Error::Invalid(format!("node {l} is missing"))))
            .collect::<Result<_>>()?;
        let mut mesh = HeatMesh { labels, xyz, elements: Vec::with_capacity(elements.len()), layout };
        for (label, kind, nodes, conductivity, capacity) in elements {
            let idx: Vec<usize> = nodes.iter().filter_map(|n| mesh.layout.node(*n)).collect();
            let xyz = idx.iter().filter_map(|i| mesh.xyz.get(*i).copied()).collect();
            mesh.elements.push(HeatElem { label, kind, nodes: idx, xyz, conductivity, capacity });
        }
        Ok(mesh)
    }

    fn node_of(&self, what: &str, label: u32) -> Result<usize> {
        self.layout.node(label).ok_or_else(|| Error::Invalid(format!("{what}: node {label} belongs to no solid element")))
    }

    fn element(&self, label: u32) -> Option<&HeatElem> {
        self.elements.binary_search_by_key(&label, |e| e.label).ok().and_then(|i| self.elements.get(i))
    }
}

/// The heat flux in one element.
#[derive(Debug, Clone, PartialEq)]
pub struct ElementFlux {
    pub label: u32,
    pub kind: ElementKind,
    pub points: Vec<PointFlux>,
}

/// The heat flow through one heat load.
#[derive(Debug, Clone, PartialEq)]
pub struct LoadFlow {
    pub name: String,
    /// `surface_flux`, `convection` or `heat_source`.
    pub kind: &'static str,
    /// Heat per time into the body (negative: out of it).
    pub heat: f64,
}

/// The result of a heat transfer step.
#[derive(Debug, Clone, PartialEq)]
pub struct HeatResult {
    pub step: String,
    /// Node labels with a temperature, sorted; per-node vectors follow this order.
    pub nodes: Vec<u32>,
    /// The time of each frame (one frame at time 0 for a steady step).
    pub times: Vec<f64>,
    /// Nodal temperatures per frame; the last frame is the end of the step.
    pub frames: Vec<Vec<f64>>,
    /// Heat flux averaged at the nodes, at the end of the step.
    pub flux: Vec<[f64; 3]>,
    /// Heat flux at the integration points, at the end of the step.
    pub elements: Vec<ElementFlux>,
    /// The heat each prescribed temperature supplies to the body at the end of the step
    /// (negative where it takes heat away; zero at free nodes).
    pub reaction: Vec<f64>,
    /// The heat each heat load carries into the body at the end of the step, in the step's
    /// order (negative where it takes heat away, as a convection to a colder sink does).
    pub loads: Vec<LoadFlow>,
    /// Heat supplied by fluxes, convection and sources (at the end temperatures), and the total
    /// reaction heat flow; in a steady state they cancel.
    pub applied_total: f64,
    pub reaction_total: f64,
    pub stats: Stats,
    pub warnings: Vec<String>,
}

impl HeatResult {
    /// The temperatures at the end of the step.
    pub fn temperature(&self) -> &[f64] {
        self.frames.last().map(Vec::as_slice).unwrap_or_default()
    }
}

/// The assembled heat system of a step: conductivity (with films), capacity when transient,
/// and the heat supplied.
struct HeatSystem {
    k: Csr,
    c: Option<Csr>,
    f: Vec<f64>,
}

fn assemble(model: &Model, resolved: &Resolved, mesh: &HeatMesh, graph: &NodeGraph, step: &Step, transient: bool) -> Result<HeatSystem> {
    let layout = &mesh.layout;
    let mut k = graph.zero_matrix(layout);
    let mut c = transient.then(|| graph.zero_matrix(layout));
    for e in &mesh.elements {
        let kernel = e.kernel()?;
        let at = |err| Error::element(e.label, err);
        add_block(&mut k, layout, &e.nodes, DofSet::TEMPERATURE, &kernel.conductivity(e.conductivity).map_err(at)?.data)?;
        if let Some(c) = c.as_mut() {
            let rho_c = e.capacity.ok_or_else(|| Error::Invalid(format!("element {}: its material needs a density and a specific heat", e.label)))?;
            add_block(c, layout, &e.nodes, DofSet::TEMPERATURE, &kernel.capacity(rho_c).map_err(at)?.data)?;
        }
    }
    let mut f = vec![0.0; layout.dofs()];
    let add = |f: &mut [f64], nodes: &[usize], v: &[f64]| {
        for (n, x) in nodes.iter().zip(v) {
            if let Some(y) = f.get_mut(*n) {
                *y += x;
            }
        }
    };
    for name in &step.loads {
        let load = model.load(name).ok_or_else(|| Error::Invalid(format!("load `{name}` is missing")))?;
        let faces = || load_faces(model, resolved, mesh, name, load);
        match load {
            Load::SurfaceFlux { flux, .. } => {
                for (label, idx, xyz) in faces()? {
                    let w = face_weights(&xyz).map_err(|e| Error::element(label, e))?;
                    add(&mut f, &idx, &w.iter().map(|w| w * flux).collect::<Vec<_>>());
                }
            }
            Load::Convection { film, sink, .. } => {
                for (label, idx, xyz) in faces()? {
                    let at = |e| Error::element(label, e);
                    add_block(&mut k, layout, &idx, DofSet::TEMPERATURE, &face_film(&xyz, *film).map_err(at)?.data)?;
                    let w = face_weights(&xyz).map_err(at)?;
                    add(&mut f, &idx, &w.iter().map(|w| w * film * sink).collect::<Vec<_>>());
                }
            }
            Load::HeatSource { set, power } => {
                for label in resolved.sets.get(set).map(|s| s.elements.as_slice()).unwrap_or_default() {
                    let e = mesh.element(*label).ok_or_else(|| Error::Invalid(format!("load `{name}`: element {label} is not a solid element")))?;
                    let q = e.kernel()?.source(*power).map_err(|err| Error::element(e.label, err))?;
                    add(&mut f, &e.nodes, &q);
                }
            }
            // Mechanical loads and temperature fields do not act on heat conduction.
            _ => {}
        }
    }
    Ok(HeatSystem { k, c, f })
}

/// The element faces of heat load `name`'s set: each face's element, node indices and
/// coordinates.
fn load_faces(model: &Model, resolved: &Resolved, mesh: &HeatMesh, name: &str, load: &Load) -> Result<Vec<(u32, Vec<usize>, Vec<[f64; 3]>)>> {
    let set = load.set().unwrap_or_default();
    resolved
        .sets
        .get(set)
        .map(|s| s.faces.as_slice())
        .unwrap_or_default()
        .iter()
        .map(|fc| {
            let labels =
                model.face_nodes(*fc).ok_or_else(|| Error::Invalid(format!("load `{name}`: element {} face {} is missing", fc.element, fc.face)))?;
            let idx: Vec<usize> = labels.iter().map(|l| mesh.node_of(&format!("load `{name}`"), *l)).collect::<Result<_>>()?;
            let xyz = idx.iter().filter_map(|i| mesh.xyz.get(*i).copied()).collect();
            Ok((fc.element, idx, xyz))
        })
        .collect()
}

/// The heat each of the step's heat loads carries into the body at temperatures `t`: a flux
/// ∫q dA, a source ∫Q dV, a convection h ∫(T∞ − T) dA (negative where it carries heat away).
fn load_flows(model: &Model, resolved: &Resolved, mesh: &HeatMesh, step: &Step, t: &[f64]) -> Result<Vec<LoadFlow>> {
    let mut out = Vec::new();
    for name in &step.loads {
        let Some(load) = model.load(name).filter(|l| l.is_thermal()) else { continue };
        let mut heat = 0.0;
        match load {
            Load::SurfaceFlux { flux, .. } => {
                for (label, _, xyz) in load_faces(model, resolved, mesh, name, load)? {
                    heat += flux * face_weights(&xyz).map_err(|e| Error::element(label, e))?.iter().sum::<f64>();
                }
            }
            // h Σ w_j (T∞ − T_j): the film matrix's columns sum to h w_j (Σ N = 1).
            Load::Convection { film, sink, .. } => {
                for (label, idx, xyz) in load_faces(model, resolved, mesh, name, load)? {
                    let w = face_weights(&xyz).map_err(|e| Error::element(label, e))?;
                    heat += film * w.iter().zip(&idx).map(|(w, i)| w * (sink - t.get(*i).copied().unwrap_or(*sink))).sum::<f64>();
                }
            }
            Load::HeatSource { set, power } => {
                for label in resolved.sets.get(set).map(|s| s.elements.as_slice()).unwrap_or_default() {
                    if let Some(e) = mesh.element(*label) {
                        heat += e.kernel()?.source(*power).map_err(|err| Error::element(e.label, err))?.iter().sum::<f64>();
                    }
                }
            }
            _ => continue,
        }
        out.push(LoadFlow { name: name.clone(), kind: load.kind(), heat });
    }
    Ok(out)
}

/// The prescribed temperature of every node (`None` = free) under the step's boundary
/// conditions.
fn prescribed(model: &Model, resolved: &Resolved, mesh: &HeatMesh, step: &Step) -> Result<Vec<Option<f64>>> {
    let mut out = vec![None; mesh.layout.dofs()];
    for name in &step.boundary_conditions {
        let bc = model.boundary_condition(name).ok_or_else(|| Error::Invalid(format!("boundary condition `{name}` is missing")))?;
        let Some(v) = bc.value(TEMPERATURE_DOF) else { continue };
        for &label in resolved.target_nodes(bc.target()) {
            let Some(n) = mesh.layout.node(label) else { continue };
            hold(out.get_mut(n), v, |old| format!("node {label} is held at both {old} and {v} (boundary condition `{name}` and another)"))?;
        }
    }
    Ok(out)
}

fn factor(a: &Csr, limits: &Limits, mesh: &HeatMesh, eqs: &Reduction) -> Result<Factor> {
    Factor::new(&UpperColumns::from_symmetric(a), limits.max_factor_entries).map_err(|e| match e {
        FactorError::TooLarge { entries, limit } => Error::factor_limit(entries, limit),
        FactorError::OutOfMemory { entries } => Error::Limit(format!("not enough memory for the factorization ({entries} entries)")),
        FactorError::Internal(what) => Error::Internal(what),
        FactorError::Cancelled => Error::Model(fieldcraft_model::Error::Cancelled),
        FactorError::NotPositiveDefinite { column } | FactorError::Singular { column } => {
            let node = eqs.eq_dof.get(column).and_then(|d| mesh.layout.locate(*d)).map(|(slot, _)| slot);
            let label = node.and_then(|n| mesh.labels.get(n)).copied().unwrap_or(0);
            Error::Singular(format!(
                "the temperature is undetermined at node {label}: no temperature boundary condition or convection reaches its part of the model"
            ))
        }
    })
}

/// y = A x.
fn mul(a: &Csr, x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; a.n];
    a.mul(x, &mut y);
    y
}

/// The temperature history of a heat transfer or coupled step and what produced it.
struct Solution {
    mesh: HeatMesh,
    eqs: Reduction,
    k: Csr,
    c: Option<Csr>,
    f: Vec<f64>,
    times: Vec<f64>,
    frames: Vec<Vec<f64>>,
    /// Ṫ over the last increment (zero when steady).
    rate: Vec<f64>,
    factor_entries: usize,
    warnings: Vec<String>,
}

impl Solution {
    fn temperature(&self) -> &[f64] {
        self.frames.last().map(Vec::as_slice).unwrap_or_default()
    }
}

/// The faces of set `set` for tie `owner`, refused unless they are linear faces of conducting
/// solid elements.
fn tie_faces(model: &Model, resolved: &Resolved, mesh: &HeatMesh, owner: &str, set: &str) -> Result<Vec<TieFace>> {
    let list = resolved.sets.get(set).map(|s| s.faces.as_slice()).unwrap_or_default();
    list.iter()
        .map(|fc| {
            let at = |what: String| Error::Invalid(format!("{owner}: element {} face {}: {what}", fc.element, fc.face));
            let labels = model.face_nodes(*fc).ok_or_else(|| at("no such face".into()))?;
            let el = mesh.element(fc.element).ok_or_else(|| at("the element is not a solid element of the analysis".into()))?;
            if !(matches!(labels.len(), 3 | 4) && el.kind.dimension() == 3) {
                return Err(Error::NotSupported(format!(
                    "{owner} on the faces of {} elements (element {}) in heat transfer; it takes the linear faces of tet4, wedge6 and hex8 elements",
                    el.kind.name(),
                    fc.element
                )));
            }
            let nodes: Vec<usize> = labels.iter().map(|l| mesh.node_of(owner, *l)).collect::<Result<_>>()?;
            let x: Vec<[f64; 3]> = nodes.iter().filter_map(|n| mesh.xyz.get(*n).copied()).collect();
            let thickness = thickness(el.kind, &el.xyz, &x, area_vector(&x).map_or(0.0, norm));
            if !(thickness > 0.0 && thickness.is_finite()) {
                return Err(at("the face has no size".into()));
            }
            Ok(TieFace { nodes, thickness })
        })
        .collect()
}

/// The temperature constraint rows of the model's ties: each secondary node's temperature
/// follows the main nodes' with the displacement tie's mortar weights ([`crate::tie`]).
fn tie_rows<'m>(model: &'m Model, resolved: &Resolved, mesh: &HeatMesh, warnings: &mut Vec<String>) -> Result<Vec<Row<'m>>> {
    if model.ties().is_empty() {
        return Ok(Vec::new());
    }
    let nodes = TieNodes { xyz: &mesh.xyz, labels: &mesh.labels, faces: |owner: &str, set: &str| tie_faces(model, resolved, mesh, owner, set) };
    let dof = |slot: usize, kind: DofKind| mesh.layout.dof(slot, kind).ok_or_else(|| Error::Internal("node temperature".into()));
    rows_of(model, &nodes, &[DofKind::T], dof, warnings)
}

/// Solve the step's temperatures (steady unless the step is transient).
///
/// The prescribed temperatures and the ties are eliminated by transformation
/// ([`crate::constraints::Reduction`]): T = T_f q + T_p T̄ (T̄ the prescribed values), with the free temperatures q numbered
/// by nested dissection. Without ties it is the plain elimination of the prescribed values.
fn run(model: &Model, resolved: &Resolved, step: &Step, limits: &Limits) -> Result<Solution> {
    let transient = match step.procedure {
        Procedure::HeatTransfer { transient } => transient,
        _ => None,
    };
    let mesh = HeatMesh::new(model, resolved, limits.max_dofs)?;
    let graph = NodeGraph::new(mesh.labels.len(), mesh.elements.iter().map(|e| e.nodes.as_slice()));
    let mut warnings = Vec::new();
    let rows = tie_rows(model, resolved, &mesh, &mut warnings)?;
    let fixed = prescribed(model, resolved, &mesh, step)?;
    let order = fill_reducing_order(&mesh.xyz, &graph.adjacency);
    let budget = Budget { max_bytes: limits.max_factor_entries.saturating_mul(8), pattern_entries: graph.entries(&mesh.layout), matrices: 1 };
    let eqs = Reduction::new(&mesh.layout, std::iter::empty(), &order, &rows, |_| Ok(fixed), budget)?;
    if eqs.neq() == 0 {
        return Err(Error::Invalid("every temperature is prescribed; there is nothing to solve".into()));
    }
    let HeatSystem { k, c, f } = assemble(model, resolved, &mesh, &graph, step, transient.is_some())?;
    let n = mesh.labels.len();
    let mut q = Vec::with_capacity(eqs.neq());
    // The temperatures from the free ones q, checked finite.
    let expand = |q: &[f64]| -> Result<Vec<f64>> {
        let t = eqs.expand(q, true);
        if t.iter().any(|x| !x.is_finite()) {
            return Err(Error::Singular("the temperatures are not finite; the model is ill-conditioned".into()));
        }
        Ok(t)
    };
    let (times, frames, rate, factor_entries) = match transient.zip(c.as_ref()) {
        Some((tr, c)) => {
            let steps = tr.increments();
            let dt = tr.end_time / steps as f64;
            // (C + θΔt K) Tₙ₊₁ = (C − (1 − θ)Δt K) Tₙ + Δt f, the prescribed values held from the
            // first increment on (their part moves to the right-hand side).
            let a = c.add_scaled(tr.theta * dt, &k)?;
            let b = c.add_scaled(-(1.0 - tr.theta) * dt, &k)?;
            let [reduced] = eqs.reduce([&a])?;
            let chol = factor(&reduced, limits, &mesh, &eqs)?;
            let every = (steps.saturating_mul(n) / MAX_TEMPERATURE_VALUES.max(1)).saturating_add(1);
            let mut t = vec![tr.initial_temperature; n];
            let mut previous = t.clone();
            let (mut frames, mut times) = (vec![t.clone()], vec![0.0]);
            let mut rhs = vec![0.0; n];
            fieldcraft_model::progress::stage("Time increment", Some(steps as f64), "");
            for i in 1..=steps {
                fieldcraft_model::progress::count((i - 1) as f64);
                fieldcraft_model::progress::check()?;
                b.mul(&t, &mut rhs);
                for (r, fj) in rhs.iter_mut().zip(&f) {
                    *r += dt * fj;
                }
                q.clear();
                q.extend(eqs.rhs(&rhs, &a));
                chol.solve(&mut q)?;
                previous = std::mem::replace(&mut t, expand(&q)?);
                if i % every == 0 || i == steps {
                    frames.push(t.clone());
                    times.push(if i == steps { tr.end_time } else { i as f64 * dt });
                }
            }
            let rate = t.iter().zip(&previous).map(|(a, b)| (a - b) / dt).collect();
            (times, frames, rate, chol.entries())
        }
        None => {
            let [reduced] = eqs.reduce([&k])?;
            let chol = factor(&reduced, limits, &mesh, &eqs)?;
            q.extend(eqs.rhs(&f, &k));
            chol.solve(&mut q)?;
            (vec![0.0], vec![expand(&q)?], vec![0.0; n], chol.entries())
        }
    };
    Ok(Solution { mesh, eqs, k, c, f, times, frames, rate, factor_entries, warnings })
}

/// The steady temperatures of a coupled temperature–displacement step, by node label (sorted).
pub(crate) fn temperatures(model: &Model, resolved: &Resolved, step: &Step, limits: &Limits) -> Result<(Vec<u32>, Vec<f64>)> {
    let mut s = run(model, resolved, step, limits)?;
    let t = s.frames.pop().unwrap_or_default();
    Ok((s.mesh.labels, t))
}

/// The steady temperatures of a step's thermal conditions at the nodes `labels` (sorted), `None`
/// at the nodes conduction does not reach.
pub(crate) fn temperatures_at(model: &Model, resolved: &Resolved, step: &Step, limits: &Limits, labels: &[u32]) -> Result<Vec<Option<f64>>> {
    let (solved, t) = temperatures(model, resolved, step, limits)?;
    Ok(labels.iter().map(|l| solved.binary_search(l).ok().and_then(|i| t.get(i)).copied()).collect())
}

/// Solve a heat transfer step: the step named `step`, or the model's only one.
pub fn solve_heat(model: &Model, step: Option<&str>, limits: &Limits) -> Result<HeatResult> {
    // A job cancelled before it starts stops here, before any setup work.
    fieldcraft_model::progress::check()?;
    let resolved = model.resolve()?;
    let step = model.select_step(step, |p| matches!(p, Procedure::HeatTransfer { .. }), "heat transfer")?;
    let mut warnings: Vec<String> = resolved.warnings.iter().map(ToString::to_string).collect();
    let mut s = run(model, &resolved, step, limits)?;
    warnings.append(&mut s.warnings);
    let t = s.temperature();
    // Reactions: the heat each prescribed node supplies, T_pᵀ (K T + C Ṫ − f) (a tied node held
    // through its main nodes passes its share to them).
    let kt = mul(&s.k, t);
    let ct = s.c.as_ref().map(|c| mul(c, &s.rate));
    let residual: Vec<f64> =
        kt.iter().zip(&s.f).enumerate().map(|(i, (kt, f))| kt + ct.as_ref().and_then(|c| c.get(i)).copied().unwrap_or(0.0) - f).collect();
    let reaction = s.eqs.reactions_of(&residual);
    // Heat supplied: fluxes and sources, and convection at the end temperatures, Σ (f − K_film T)
    // = Σ (f − K T) since conduction alone has zero row sums.
    let applied_total: f64 = kt.iter().zip(&s.f).map(|(a, b)| b - a).sum();
    let reaction_total: f64 = reaction.iter().sum();
    let (flux, elements) = fluxes(&s.mesh, t)?;
    let loads = load_flows(model, &resolved, &s.mesh, step, t)?;
    // Measured against the gross heat flows, which cancel in the net sums.
    let scale = reaction.iter().chain(&s.f).map(|v| v.abs()).sum::<f64>();
    if s.c.is_none() && scale > 0.0 && (applied_total + reaction_total).abs() > crate::system::BALANCE_TOLERANCE * scale {
        warnings.push(format!("heat supplied and removed differ by {:.3e}", (applied_total + reaction_total).abs()));
    }
    let stats = Stats { dofs: s.mesh.layout.dofs(), equations: s.eqs.neq(), factor_entries: s.factor_entries, lanczos_vectors: 0 };
    Ok(HeatResult {
        step: step.name.clone(),
        nodes: s.mesh.labels,
        times: s.times,
        frames: s.frames,
        flux,
        elements,
        reaction,
        loads,
        applied_total,
        reaction_total,
        stats,
        warnings,
    })
}

/// The heat flux at the integration points and averaged at the nodes (one extrapolation per
/// element kind).
fn fluxes(mesh: &HeatMesh, t: &[f64]) -> Result<(Vec<[f64; 3]>, Vec<ElementFlux>)> {
    let mut sum = vec![[0.0; 3]; mesh.labels.len()];
    let mut count = vec![0u32; mesh.labels.len()];
    let mut out = Vec::with_capacity(mesh.elements.len());
    let mut extrapolations: Vec<(ElementKind, Extrapolation)> = Vec::new();
    for e in &mesh.elements {
        let kernel = e.kernel()?;
        let at = |err| Error::element(e.label, err);
        let te: Vec<f64> = e.nodes.iter().map(|n| t.get(*n).copied().unwrap_or(0.0)).collect();
        let points = kernel.flux(e.conductivity, &te).map_err(at)?;
        let padded: Vec<[f64; 6]> = points.iter().map(|p| [p.flux[0], p.flux[1], p.flux[2], 0.0, 0.0, 0.0]).collect();
        let index = match extrapolations.iter().position(|(k, _)| *k == e.kind) {
            Some(i) => i,
            None => {
                extrapolations.push((e.kind, extrapolation(kernel.ty, kernel.ty.full).map_err(at)?));
                extrapolations.len() - 1
            }
        };
        let nodal = extrapolations.get(index).map(|(_, x)| x.apply(&padded)).unwrap_or_default();
        for (n, v) in e.nodes.iter().zip(&nodal) {
            if let (Some(s), Some(c)) = (sum.get_mut(*n), count.get_mut(*n)) {
                for k in 0..3 {
                    s[k] += v[k];
                }
                *c += 1;
            }
        }
        out.push(ElementFlux { label: e.label, kind: e.kind, points });
    }
    let flux = sum.iter().zip(&count).map(|(s, c)| s.map(|v| if *c > 0 { v / f64::from(*c) } else { 0.0 })).collect();
    Ok((flux, out))
}
