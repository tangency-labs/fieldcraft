//! Heat transfer: steady and transient conduction with prescribed temperatures, surface fluxes,
//! convection, radiation to an ambient and internal heat sources.
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
//! - **Radiation** to an ambient, a heat loss r(T) = ∫ N ε σ (T⁴ − T_amb⁴) dA from each node of
//!   the faces, is nonlinear: each solve (the steady one, or each increment, with r weighted by
//!   θ as K is) runs Newton's method on it, r(T) ≈ r(Tᵏ) + J (T − Tᵏ) with the tangent
//!   J = ∫ N Nᵀ 4 ε σ T³ dA, refactoring K (+ C) + J each iteration, until the heat out of balance
//!   in every free equation is at most [`RADIATION_TOLERANCE`] of the largest heat flow, above
//!   the equation's own round-off (at
//!   most [`MAX_RADIATION_ITERATIONS`] per solve, and at most [`Limits::max_newton_iterations`]:
//!   every increment has its own count, so a long transient is bounded by its increments). A
//!   steady solve starts from the lumped balance: the uniform temperature at which the faces
//!   radiate the heat supplied (or the highest prescribed temperature, if higher), and is refused
//!   when that is 0 K (no heat, no held temperature, surroundings at 0 K). Temperatures are
//!   absolute (the model is in kelvin); faces see only the ambient, not each other. Without
//!   radiation the solve is linear and factored once, as above.
//!
//! The heat flux −k∇T is reported at the integration points and averaged at the nodes; the heat
//! flow at prescribed temperatures is the reaction K T − f (C Ṫ included in transient steps).
//!
//! **Temperature-dependent properties** ([`fieldcraft_model::curves`]): a material with
//! conductivity, specific heat or density tables makes the problem nonlinear, K(T) T = f. It is
//! solved by Picard (successive substitution) iteration: each element takes its properties at the
//! mean of its nodal temperatures, the system is assembled and solved again, until no temperature
//! moves by more than [`PICARD_TOLERANCE`] of the largest; the steady step starts from the
//! constants, each transient increment from the previous one's temperatures and properties
//! (fully implicit in the properties). At most [`MAX_PICARD_ITERATIONS`] solves per step or
//! increment (without radiation); a run that stops there says so in its warnings. With radiation
//! as well, the two share one loop ([`Problem::settle`]): each Newton iteration first takes the
//! properties at the last temperatures, and the solve ends when the heat balances and the
//! properties have settled, within [`MAX_RADIATION_ITERATIONS`]. With a conductivity linear in the
//! temperature the Kirchhoff transform θ = ∫k dT is linear in x across a slab, and the element
//! mean is exact for linear elements, so the nodal temperatures are exact (H. S. Carslaw, J. C.
//! Jaeger, *Conduction of Heat in Solids*, 2nd ed., Oxford 1959, §1.7).

use std::borrow::Cow;

use fieldcraft_fem::element::{Extrapolation, extrapolation};
use fieldcraft_fem::frame::norm;
use fieldcraft_fem::heat::{Conduction, PointFlux, face_film, face_radiation, face_weights};
use fieldcraft_fem::mortar::area_vector;
use fieldcraft_fem::{DofKind, DofLayout, DofSet};
use fieldcraft_model::{
    ElementKind, Load, Model, Procedure, Resolved, STEFAN_BOLTZMANN, SectionKind, Step, TEMPERATURE_DOF, TableRange, TemperatureDependent,
    TemperatureTable,
};

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

/// Most Newton iterations of one solve with radiation (the steady solve or one increment), the
/// temperature-dependent properties' iteration included.
pub const MAX_RADIATION_ITERATIONS: usize = 50;

/// A radiation solve has converged when the heat out of balance in each free equation is at most
/// this share of the largest net heat flow at a node (supplied, conducted or stored, and
/// radiated), above the equation's round-off ([`ROUND_OFF`]).
pub const RADIATION_TOLERANCE: f64 = 1e-10;

/// The round-off allowed in a free equation's heat balance, in units of ε times the sizes of its
/// terms (|g| + |A| |T| + s |r|): where large conductances meet nearly equal temperatures, the
/// balance cannot be computed more closely than that, however well the temperatures are solved.
const ROUND_OFF: f64 = 16.0 * f64::EPSILON;

/// Most solves of a step (steady) or an increment (transient) with temperature-dependent
/// properties and no radiation.
pub const MAX_PICARD_ITERATIONS: usize = 100;

/// The Picard iteration has converged when no temperature moves by more than this fraction of
/// the largest absolute temperature.
pub const PICARD_TOLERANCE: f64 = 1e-10;

/// The temperature tables of one material, beside the constants they replace.
#[derive(Debug, Clone)]
pub struct Variable {
    pub material: String,
    /// Its conductivity, specific heat and density tables (no viscosity).
    tables: TemperatureDependent,
    /// The constants: conductivity along X, Y, Z, specific heat and density.
    k: [f64; 3],
    c: Option<f64>,
    rho: Option<f64>,
}

impl Variable {
    /// The conductivity and the capacity ρc at `kelvin`.
    fn at(&self, kelvin: f64) -> ([f64; 3], Option<f64>) {
        let value = |t: &Option<TemperatureTable>, constant: Option<f64>| t.as_ref().map(|t| t.value(kelvin)).or(constant);
        let d = &self.tables;
        // A conductivity table is isotropic; without one the constant keeps its axes.
        let k = d.conductivity.as_ref().map_or(self.k, |t| [t.value(kelvin); 3]);
        (k, value(&d.density, self.rho).zip(value(&d.specific_heat, self.c)).map(|(rho, c)| rho * c))
    }
}

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
    /// The index of its material's temperature tables in [`HeatMesh::variables`].
    pub variable: Option<usize>,
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
    /// The materials with temperature tables (heat transfer is then nonlinear).
    pub variables: Vec<Variable>,
}

impl HeatMesh {
    pub fn new(model: &Model, resolved: &Resolved, max_dofs: usize) -> Result<HeatMesh> {
        let mut elements = Vec::new();
        let mut variables: Vec<Variable> = Vec::new();
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
            let variable = material.temperature_dependent.as_ref().filter(|d| d.conducts()).map(|d| {
                match variables.iter().position(|v| v.material == section.material) {
                    Some(i) => i,
                    None => {
                        variables.push(Variable {
                            material: section.material.clone(),
                            tables: TemperatureDependent { viscosity: None, ..d.clone() },
                            k: conductivity,
                            c: thermal.specific_heat,
                            rho: material.density,
                        });
                        variables.len() - 1
                    }
                }
            });
            let set = resolved.sets.get(&section.set).ok_or_else(|| Error::Invalid(format!("set `{}` is not resolved", section.set)))?;
            for &label in &set.elements {
                let el = model.element(label).ok_or_else(|| Error::Invalid(format!("element {label} is missing")))?;
                elements.push((label, el.kind, el.nodes.clone(), conductivity, capacity, variable));
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
        let mut mesh = HeatMesh { labels, xyz, elements: Vec::with_capacity(elements.len()), layout, variables };
        for (label, kind, nodes, conductivity, capacity, variable) in elements {
            let idx: Vec<usize> = nodes.iter().filter_map(|n| mesh.layout.node(*n)).collect();
            let xyz = idx.iter().filter_map(|i| mesh.xyz.get(*i).copied()).collect();
            mesh.elements.push(HeatElem { label, kind, nodes: idx, xyz, conductivity, capacity, variable });
        }
        Ok(mesh)
    }

    /// The mean of an element's nodal temperatures `t`.
    fn mean(e: &HeatElem, t: &[f64]) -> f64 {
        e.nodes.iter().filter_map(|n| t.get(*n)).sum::<f64>() / e.nodes.len().max(1) as f64
    }

    /// Take each temperature-dependent element's properties at its mean temperature in `t`.
    fn update(&mut self, t: &[f64]) {
        for e in &mut self.elements {
            if let Some(v) = e.variable.and_then(|i| self.variables.get(i)) {
                (e.conductivity, e.capacity) = v.at(HeatMesh::mean(e, t));
            }
        }
    }

    /// The materials whose element temperatures in any of `frames` leave the range of a table
    /// the step reads (a steady one the conductivity alone).
    fn ranges(&self, frames: &[Vec<f64>], transient: bool) -> Vec<TableRange> {
        let mut span = vec![[f64::INFINITY, f64::NEG_INFINITY]; self.variables.len()];
        for t in frames {
            for e in &self.elements {
                if let Some(s) = e.variable.and_then(|i| span.get_mut(i)) {
                    let mean = HeatMesh::mean(e, t);
                    *s = [s[0].min(mean), s[1].max(mean)];
                }
            }
        }
        self.variables
            .iter()
            .zip(span)
            .filter_map(|(v, [low, high])| {
                if transient {
                    v.tables.range(&v.material, low, high)
                } else {
                    let conduction = TemperatureDependent { conductivity: v.tables.conductivity.clone(), ..TemperatureDependent::default() };
                    conduction.range(&v.material, low, high)
                }
            })
            .collect()
    }

    fn node_of(&self, what: &str, label: u32) -> Result<usize> {
        self.layout.node(label).ok_or_else(|| Error::Invalid(format!("{what}: node {label} belongs to no solid element")))
    }

    fn element(&self, label: u32) -> Option<&HeatElem> {
        self.elements.binary_search_by_key(&label, |e| e.label).ok().and_then(|i| self.elements.get(i))
    }
}

/// Whether the Picard iterate `next` has settled from `previous`.
fn settled(previous: &[f64], next: &[f64]) -> bool {
    let scale = next.iter().fold(0.0f64, |m, t| m.max(t.abs()));
    previous.iter().zip(next).all(|(a, b)| (a - b).abs() <= PICARD_TOLERANCE * scale)
}

/// The warning of a Picard iteration that stopped at its cap of `solves`; `what` names the step
/// or increment.
fn unsettled(what: &str, solves: usize) -> String {
    format!(
        "the temperature-dependent properties did not settle in {solves} solves ({what}): the temperatures are those of the last solve; a property that changes steeply with the temperature, or shorter increments, need more"
    )
}

/// The heat flux in one element.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ElementFlux {
    pub label: u32,
    pub kind: ElementKind,
    pub points: Vec<PointFlux>,
}

/// The heat flow through one heat load.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct LoadFlow {
    pub name: String,
    /// `surface_flux`, `convection`, `radiation` or `heat_source`.
    pub kind: String,
    /// Heat per time into the body (negative: out of it).
    pub heat: f64,
}

/// The result of a heat transfer step.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
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
    /// order (negative where it takes heat away, as a convection to a colder sink or a
    /// radiation to colder surroundings does).
    pub loads: Vec<LoadFlow>,
    /// Heat supplied by fluxes, convection, radiation and sources (at the end temperatures), and the total
    /// reaction heat flow; in a steady state they cancel.
    pub applied_total: f64,
    pub reaction_total: f64,
    pub stats: Stats,
    pub warnings: Vec<String>,
    /// The temperature tables read beyond their ranges, in kelvin: warnings the caller words
    /// in its display system. Reported by the solve, not saved: the results encoding skips it
    /// (none when read back), so results saved before it read unchanged.
    #[serde(skip)]
    pub table_ranges: Vec<TableRange>,
}

impl HeatResult {
    /// The temperatures at the end of the step.
    pub fn temperature(&self) -> &[f64] {
        self.frames.last().map(Vec::as_slice).unwrap_or_default()
    }

    /// Nodal heat flux in `frame`, in the order of [`Self::nodes`].
    ///
    /// The model must be unchanged since the solve. Earlier frames recover −k∇T from
    /// their temperatures with the same extrapolation and averaging as the final frame;
    /// the final frame borrows the stored flux. No additional frame history is kept.
    pub fn flux_at(&self, model: &Model, frame: usize) -> Result<Cow<'_, [[f64; 3]]>> {
        let t = self.frames.get(frame).ok_or_else(|| Error::Invalid(format!("step `{}` has no frame {frame}", self.step)))?;
        if frame + 1 == self.frames.len() {
            return Ok(Cow::Borrowed(&self.flux));
        }
        if t.len() != self.nodes.len() || t.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid(format!("step `{}` frame {frame} has invalid temperatures", self.step)));
        }
        let mut mesh = HeatMesh::new(model, &model.resolve()?, self.nodes.len())?;
        if mesh.labels != self.nodes {
            return Err(Error::Invalid(format!("step `{}` temperatures do not match the conducting mesh", self.step)));
        }
        // The conductivity of the frame's own temperatures, where the material follows a table.
        mesh.update(t);
        Ok(Cow::Owned(fluxes(&mesh, t)?.0))
    }
}

/// The assembled heat system of a step: conductivity (with films), capacity when transient,
/// the heat supplied, and the faces that radiate.
struct HeatSystem {
    k: Csr,
    c: Option<Csr>,
    f: Vec<f64>,
    radiation: Vec<RadiatingFace>,
}

/// A face radiating to an ambient: its element, node indices and coordinates, ε σ, the
/// ambient temperature and its area.
struct RadiatingFace {
    label: u32,
    nodes: Vec<usize>,
    xyz: Vec<[f64; 3]>,
    es: f64,
    ambient: f64,
    area: f64,
}

/// The faces of heat load `name` when it is a radiation (none otherwise).
fn radiating_faces(model: &Model, resolved: &Resolved, mesh: &HeatMesh, name: &str, load: &Load) -> Result<Vec<RadiatingFace>> {
    let Load::Radiation { emissivity, ambient, .. } = load else { return Ok(Vec::new()) };
    load_faces(model, resolved, mesh, name, load)?
        .into_iter()
        .map(|(label, nodes, xyz)| {
            let area = face_weights(&xyz).map_err(|e| Error::element(label, e))?.iter().sum();
            Ok(RadiatingFace { label, nodes, xyz, es: emissivity * STEFAN_BOLTZMANN, ambient: *ambient, area })
        })
        .collect()
}

/// The heat each node loses by radiation at temperatures `t`; with `tangent` = (M, s), also
/// s J(T) added into M (the conduction's pattern) and J(T) T returned.
fn radiate(faces: &[RadiatingFace], layout: &DofLayout, t: &[f64], mut tangent: Option<(&mut Csr, f64)>) -> Result<(Vec<f64>, Vec<f64>)> {
    let mut loss = vec![0.0; t.len()];
    let mut jt = vec![0.0; if tangent.is_some() { t.len() } else { 0 }];
    for face in faces {
        let te: Vec<f64> = face.nodes.iter().map(|n| t.get(*n).copied().unwrap_or(0.0)).collect();
        let (l, mut m) = face_radiation(&face.xyz, &te, face.es, face.ambient).map_err(|e| Error::element(face.label, e))?;
        let mt = m.mul_vec(&te);
        for (k, n) in face.nodes.iter().enumerate() {
            if let Some(x) = loss.get_mut(*n) {
                *x += l.get(k).copied().unwrap_or(0.0);
            }
            if let Some(x) = jt.get_mut(*n) {
                *x += mt.get(k).copied().unwrap_or(0.0);
            }
        }
        if let Some((matrix, scale)) = tangent.as_mut() {
            m.data.iter_mut().for_each(|v| *v *= *scale);
            add_block(matrix, layout, &face.nodes, DofSet::TEMPERATURE, &m.data)?;
        }
    }
    Ok((loss, jt))
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
    let mut radiation = Vec::new();
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
            Load::Radiation { .. } => radiation.extend(radiating_faces(model, resolved, mesh, name, load)?),
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
    Ok(HeatSystem { k, c, f, radiation })
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
/// ∫q dA, a source ∫Q dV, a convection h ∫(T∞ − T) dA, a radiation ε σ ∫(T_amb⁴ − T⁴) dA
/// (negative where it carries heat away).
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
            Load::Radiation { .. } => {
                let faces = radiating_faces(model, resolved, mesh, name, load)?;
                heat -= radiate(&faces, &mesh.layout, t, None)?.0.iter().sum::<f64>();
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
        out.push(LoadFlow { name: name.clone(), kind: load.kind().to_string(), heat });
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
                "the temperature is undetermined at node {label}: no temperature boundary condition, convection or radiation reaches its part of the model"
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
    /// The heat each node loses by radiation at the end temperatures (zero without radiation).
    radiated: Vec<f64>,
    factor_entries: usize,
    warnings: Vec<String>,
    ranges: Vec<TableRange>,
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
    let mut mesh = HeatMesh::new(model, resolved, limits.max_dofs)?;
    let graph = NodeGraph::new(mesh.labels.len(), mesh.elements.iter().map(|e| e.nodes.as_slice()));
    let mut warnings = Vec::new();
    let rows = tie_rows(model, resolved, &mesh, &mut warnings)?;
    let fixed = prescribed(model, resolved, &mesh, step)?;
    let held = fixed.iter().flatten().fold(0.0_f64, |m, v| m.max(*v));
    let order = fill_reducing_order(&mesh.xyz, &graph.adjacency);
    let budget = Budget { max_bytes: limits.max_factor_entries.saturating_mul(8), pattern_entries: graph.entries(&mesh.layout), matrices: 1 };
    let eqs = Reduction::new(&mesh.layout, std::iter::empty(), &order, &rows, |_| Ok(fixed), budget)?;
    if eqs.neq() == 0 {
        return Err(Error::Invalid("every temperature is prescribed; there is nothing to solve".into()));
    }
    // A transient step's properties start at the initial temperature, a steady one's at the
    // constants.
    let taken = transient.filter(|_| !mesh.variables.is_empty()).map(|tr| vec![tr.initial_temperature; mesh.labels.len()]);
    if let Some(t) = &taken {
        mesh.update(t);
    }
    let HeatSystem { k, c, f, radiation } = assemble(model, resolved, &mesh, &graph, step, transient.is_some())?;
    let n = mesh.labels.len();
    let mut p = Problem { model, resolved, step, limits, graph, mesh, eqs, radiation, k, c, f, taken, factored: None, entries: 0 };
    let (times, frames, rate, radiated) = match transient.filter(|_| p.c.is_some()) {
        Some(tr) => {
            let steps = tr.increments();
            let dt = tr.end_time / steps as f64;
            let every = (steps.saturating_mul(n) / MAX_TEMPERATURE_VALUES.max(1)).saturating_add(1);
            let mut t = vec![tr.initial_temperature; n];
            let mut radiated = radiate(&p.radiation, &p.mesh.layout, &t, None)?.0;
            let mut previous = t.clone();
            let (mut frames, mut times) = (vec![t.clone()], vec![0.0]);
            let mut stopped = None;
            // A = C + θΔt K and B = C − (1 − θ)Δt K, built once when nothing makes them follow
            // the temperature.
            let scaled = |p: &Problem| -> Result<(Csr, Csr)> {
                let c = p.c.as_ref().ok_or_else(|| Error::Internal("transient heat capacity".into()))?;
                Ok((c.add_scaled(tr.theta * dt, &p.k)?, c.add_scaled(-(1.0 - tr.theta) * dt, &p.k)?))
            };
            let fixed = if p.mesh.variables.is_empty() && p.radiation.is_empty() { Some(scaled(&p)?) } else { None };
            fieldcraft_model::progress::stage("Time increment", Some(steps as f64), "");
            for i in 1..=steps {
                fieldcraft_model::progress::count((i - 1) as f64);
                // (C + θΔt K) Tₙ₊₁ = (C − (1 − θ)Δt K) Tₙ + Δt f, the prescribed values held from the
                // first increment on (their part moves to the right-hand side), the radiation
                // weighted alike, θ r(Tₙ₊₁) + (1 − θ) r(Tₙ). C and K are those of the properties
                // last taken (fully implicit in the properties).
                let system = |p: &Problem| {
                    let (a, mut rhs) = match &fixed {
                        Some((a, b)) => (Cow::Borrowed(a), mul(b, &t)),
                        None => {
                            let (a, b) = scaled(p)?;
                            (Cow::Owned(a), mul(&b, &t))
                        }
                    };
                    for ((r, fj), lost) in rhs.iter_mut().zip(&p.f).zip(&radiated) {
                        *r += dt * fj - (1.0 - tr.theta) * dt * lost;
                    }
                    Ok::<_, Error>((a, rhs))
                };
                let (next, lost, settled) = p.settle(t.clone(), tr.theta * dt, system).map_err(|e| at_increment(e, i))?;
                if !settled {
                    stopped.get_or_insert(i);
                }
                radiated = lost;
                previous = std::mem::replace(&mut t, next);
                if i % every == 0 || i == steps {
                    frames.push(t.clone());
                    times.push(if i == steps { tr.end_time } else { i as f64 * dt });
                }
            }
            if let Some(i) = stopped {
                warnings.push(unsettled(&format!("first in increment {i}"), p.cap()));
            }
            let rate = t.iter().zip(&previous).map(|(a, b)| (a - b) / dt).collect();
            (times, frames, rate, radiated)
        }
        None => {
            let start = if p.radiation.is_empty() {
                vec![0.0; n]
            } else {
                let start = steady_start(&p.radiation, &p.f, held);
                if start <= 0.0 {
                    return Err(Error::Invalid(
                        "the radiation's surroundings are at 0 K and nothing supplies heat or holds a temperature: the steady temperature is 0 K everywhere, where the radiation fixes nothing; add heat, hold a temperature or raise the ambient".into(),
                    ));
                }
                vec![start; n]
            };
            let (t, radiated, settled) = p.settle(start, 1.0, |p: &Problem| Ok((Cow::Owned(p.k.clone()), p.f.clone())))?;
            if !settled {
                warnings.push(unsettled("the steady step", p.cap()));
            }
            (vec![0.0], vec![t], vec![0.0; n], radiated)
        }
    };
    let ranges = p.mesh.ranges(&frames, transient.is_some());
    let Problem { mesh, eqs, k, c, f, entries: factor_entries, .. } = p;
    Ok(Solution { mesh, eqs, k, c, f, times, frames, rate, radiated, factor_entries, warnings, ranges })
}

/// The increment an error of a transient step happened in, when it is a radiation solve that
/// did not converge.
fn at_increment(e: Error, i: usize) -> Error {
    match e {
        Error::NotConverged(what) => Error::NotConverged(format!("increment {i}: {what}")),
        e => e,
    }
}

/// Where a steady radiation solve starts: the uniform temperature at which the faces radiate
/// the heat supplied, (Σf / Σ ε σ A + mean T_amb⁴)^¼, or the highest prescribed temperature
/// when that is higher.
fn steady_start(faces: &[RadiatingFace], f: &[f64], held: f64) -> f64 {
    let es_area: f64 = faces.iter().map(|face| face.es * face.area).sum();
    let ambient4: f64 = faces.iter().map(|face| face.es * face.area * face.ambient.powi(4)).sum();
    let supplied = f.iter().sum::<f64>().max(0.0);
    let lumped = if es_area > 0.0 { ((supplied + ambient4) / es_area).powf(0.25) } else { 0.0 };
    if lumped.is_finite() { lumped.max(held) } else { held }
}

/// The heat problem of one step as it is iterated: the mesh with its properties as last taken,
/// the matrices assembled from them, and the radiating faces.
struct Problem<'a> {
    model: &'a Model,
    resolved: &'a Resolved,
    step: &'a Step,
    limits: &'a Limits,
    graph: NodeGraph,
    mesh: HeatMesh,
    eqs: Reduction,
    radiation: Vec<RadiatingFace>,
    k: Csr,
    c: Option<Csr>,
    f: Vec<f64>,
    /// The temperatures the temperature-dependent properties were last taken at (`None`: the
    /// constants).
    taken: Option<Vec<f64>>,
    /// The factor of the last linear system, kept while its matrix stays the same.
    factored: Option<Factor>,
    /// The size of the last factor.
    entries: usize,
}

impl Problem<'_> {
    /// Most solves of one [`Problem::settle`]: Newton's with radiation, Picard's without.
    fn cap(&self) -> usize {
        if self.radiation.is_empty() { MAX_PICARD_ITERATIONS } else { MAX_RADIATION_ITERATIONS.min(self.limits.max_newton_iterations) }
    }

    /// Reduce, factor and solve M T = rhs (the prescribed values applied), the temperatures
    /// checked finite.
    fn solve(&mut self, m: &Csr, rhs: &[f64], keep: bool) -> Result<Vec<f64>> {
        let chol = match self.factored.take() {
            Some(chol) => chol,
            None => {
                let [reduced] = self.eqs.reduce([m])?;
                factor(&reduced, self.limits, &self.mesh, &self.eqs)?
            }
        };
        let mut q = self.eqs.rhs(rhs, m);
        chol.solve(&mut q)?;
        self.entries = chol.entries();
        if keep {
            self.factored = Some(chol);
        }
        let t = self.eqs.expand(&q, true);
        if t.iter().any(|x| !x.is_finite()) {
            return Err(Error::Singular("the temperatures are not finite; the model is ill-conditioned".into()));
        }
        Ok(t)
    }

    /// Solve A(T) T + s r(T) = g(T) from `start`, `system` giving A and g of the properties last
    /// taken, r the radiation. Both nonlinearities share one loop: each iteration takes the
    /// temperature-dependent properties at the last temperatures (unless they settled there), and
    /// Newton's method linearizes the radiation there, (A + s J(Tᵏ)) Tᵏ⁺¹ = g − s (r(Tᵏ) − J(Tᵏ) Tᵏ),
    /// until every free equation's residual g − A T − s r(T) is at most [`RADIATION_TOLERANCE`] of
    /// the largest net heat flow at a node, |g| + |A T| + s |r|, plus its own round-off,
    /// [`ROUND_OFF`] (|g| + |A| |T| + s |r|), and the last solve moved no temperature by more
    /// than [`PICARD_TOLERANCE`] from those its properties were taken at. At most
    /// [`Problem::cap`] solves; radiation still out of balance then is an error, properties that
    /// have not settled are reported by the `false` returned. Returns T, r(T) and whether the
    /// properties settled.
    fn settle<'f>(
        &mut self,
        start: Vec<f64>,
        s: f64,
        system: impl Fn(&Problem) -> Result<(Cow<'f, Csr>, Vec<f64>)>,
    ) -> Result<(Vec<f64>, Vec<f64>, bool)> {
        let radiating = !self.radiation.is_empty();
        let variable = !self.mesh.variables.is_empty();
        let cap = self.cap();
        let mut t = start;
        let (mut a, mut g) = system(self)?;
        let mut lost = vec![0.0; t.len()];
        let (mut solved, mut steady, mut balanced) = (false, !variable, !radiating);
        let (mut residual, mut flows) = (f64::INFINITY, 0.0);
        for k in 0..=cap {
            fieldcraft_model::progress::check()?;
            if variable && solved {
                steady = self.taken.as_ref().is_some_and(|p| settled(p, &t));
                // (Not at the last check: nothing would solve with them.)
                if !steady && k < cap {
                    self.mesh.update(&t);
                    let HeatSystem { k, c, f, .. } = assemble(self.model, self.resolved, &self.mesh, &self.graph, self.step, self.c.is_some())?;
                    (self.k, self.c, self.f, self.taken, self.factored) = (k, c, f, Some(t.clone()), None);
                    (a, g) = system(self)?;
                }
            }
            if !radiating {
                if solved && steady {
                    return Ok((t, lost, true));
                }
                if k == cap {
                    break;
                }
                t = self.solve(&a, &g, true)?;
                solved = true;
                continue;
            }
            let mut m = a.as_ref().clone();
            let jt;
            (lost, jt) = radiate(&self.radiation, &self.mesh.layout, &t, Some((&mut m, s)))?;
            let at = mul(&a, &t);
            // |A|, entry by entry: the size of the conducted (and stored) heat.
            let mut size = a.as_ref().clone();
            size.val.iter_mut().for_each(|v| *v = v.abs());
            let sizes = mul(&size, &t.iter().map(|v| v.abs()).collect::<Vec<_>>());
            let out: Vec<f64> = g.iter().zip(&at).zip(&lost).map(|((g, a), r)| g - a - s * r).collect();
            flows = g.iter().zip(&at).zip(&lost).fold(0.0_f64, |m, ((g, a), r)| m.max(g.abs() + a.abs() + s * r.abs()));
            let terms: Vec<f64> = g.iter().zip(&sizes).zip(&lost).map(|((g, c), r)| g.abs() + c + s * r.abs()).collect();
            // (A tie's weights may be negative, so a restricted size can only come out smaller:
            // stricter, never looser.)
            let (out, terms) = (self.eqs.restrict(&out), self.eqs.restrict(&terms));
            residual = out.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            balanced = out.iter().zip(&terms).all(|(r, size)| r.abs() <= RADIATION_TOLERANCE * flows + ROUND_OFF * size.abs());
            if balanced && steady {
                return Ok((t, lost, true));
            }
            if k == cap {
                break;
            }
            let rhs: Vec<f64> = g.iter().zip(lost.iter().zip(&jt)).map(|(g, (r, jt))| g - s * (r - jt)).collect();
            t = self.solve(&m, &rhs, false)?;
            solved = true;
        }
        if !balanced {
            return Err(Error::NotConverged(format!(
                "the radiation did not settle in {cap} Newton iterations (heat out of balance {residual:.3e} against flows of {flows:.3e})"
            )));
        }
        Ok((t, lost, false))
    }
}

/// The steady temperatures of a step's thermal conditions at the nodes `labels` (sorted), `None`
/// at the nodes conduction does not reach; the heat solution's own warnings (properties that did
/// not settle) go to `warnings`, the tables it read beyond their ranges to `ranges`.
pub(crate) fn temperatures_at(
    model: &Model,
    resolved: &Resolved,
    step: &Step,
    limits: &Limits,
    labels: &[u32],
    warnings: &mut Vec<String>,
    ranges: &mut Vec<TableRange>,
) -> Result<Vec<Option<f64>>> {
    let mut s = run(model, resolved, step, limits)?;
    warnings.extend(s.warnings.into_iter().map(|w| format!("the step's heat conduction: {w}")));
    ranges.append(&mut s.ranges);
    let t = s.frames.pop().unwrap_or_default();
    Ok(labels.iter().map(|l| s.mesh.labels.binary_search(l).ok().and_then(|i| t.get(i)).copied()).collect())
}

/// Solve a heat transfer step: the step named `step`, or the model's only one.
pub fn solve_heat(model: &Model, step: Option<&str>, limits: &Limits) -> Result<HeatResult> {
    // A job cancelled before it starts stops here, before any setup work.
    fieldcraft_model::progress::check()?;
    let resolved = model.resolve()?;
    let step = model.select_step(step, |p| matches!(p, Procedure::HeatTransfer { .. }), "heat transfer")?;
    let mut warnings: Vec<String> = resolved.warnings.iter().map(ToString::to_string).collect();
    warnings.extend(fieldcraft_model::initial::left_out(model, step));
    let mut s = run(model, &resolved, step, limits)?;
    warnings.append(&mut s.warnings);
    let t = s.temperature();
    // Reactions: the heat each prescribed node supplies, T_pᵀ (K T + r + C Ṫ − f) (a tied node
    // held through its main nodes passes its share to them).
    let kt = mul(&s.k, t);
    // K T plus what radiates away.
    let out: Vec<f64> = kt.iter().zip(&s.radiated).map(|(a, b)| a + b).collect();
    let ct = s.c.as_ref().map(|c| mul(c, &s.rate));
    let residual: Vec<f64> =
        out.iter().zip(&s.f).enumerate().map(|(i, (o, f))| o + ct.as_ref().and_then(|c| c.get(i)).copied().unwrap_or(0.0) - f).collect();
    let reaction = s.eqs.reactions_of(&residual);
    // Heat supplied: fluxes and sources, and convection and radiation at the end temperatures,
    // Σ (f − K_film T − r) = Σ (f − K T − r) since conduction alone has zero row sums.
    let applied_total: f64 = out.iter().zip(&s.f).map(|(a, b)| b - a).sum();
    let reaction_total: f64 = reaction.iter().sum();
    let (flux, elements) = fluxes(&s.mesh, t)?;
    let loads = load_flows(model, &resolved, &s.mesh, step, t)?;
    // Measured against the gross heat flows, which cancel in the net sums.
    let scale = reaction.iter().chain(&s.f).chain(&s.radiated).map(|v| v.abs()).sum::<f64>();
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
        table_ranges: s.ranges,
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
