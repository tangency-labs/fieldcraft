//! Solver decks: the solver-neutral model written as Abaqus input, Nastran bulk data and LS-DYNA
//! keyword files. Each format is its own module (`abaqus`, `nastran`, `lsdyna`), written from the
//! vendor's public reference manual in our own words; examples and expected files are our own.
//! Writers read [`fieldcraft_model::Model`] and its `resolve()` result, never another format.
//!
//! The model stores SI; a deck is written in one consistent unit system ([`Preset`]), stated in
//! its header: [`express`] converts the model into it and the writers write its numbers as they
//! are (Nastran's `PARAM,WTMASS,1.0` and LS-DYNA's `*CONTROL_UNITS` say so too).
//!
//! Writers produce text; saving it is the caller's business (the engine's `deck.write`), so the
//! crate stays free of the file system. [`FORMATS`] lists every writer by the name `deck.write`
//! takes.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]

pub mod abaqus;
pub mod lsdyna;
pub mod nastran;

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_model::nodal::edge_pressure;
use fieldcraft_model::units::{DisplaySystem, Preset};
use fieldcraft_model::{BoundaryCondition, Diagnostic, Load, Model, Resolved, ResolvedDirection};

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum Error {
    /// The model has something this format cannot express (named, with the format).
    #[error("not supported by {format}: {what}")]
    Unsupported { format: &'static str, what: String },
    /// The model is not ready to write (its check failed).
    #[error(transparent)]
    Model(#[from] fieldcraft_model::Error),
    #[error("limit reached: {0}")]
    Limit(String),
    /// Invalid or unsupported input data, at its original one-based physical line.
    #[error("Abaqus input line {line}: {message}")]
    Parse { line: usize, message: String },
}

pub type Result<T> = std::result::Result<T, Error>;

/// A deck format.
pub struct Format {
    /// The name commands use (`abaqus`, …).
    pub name: &'static str,
    /// The usual file extension, without the dot.
    pub extension: &'static str,
    writer: fn(&Model, &Resolved, Preset) -> Result<String>,
}

/// A written deck: its text and the model check's warnings.
#[derive(Debug, Clone, PartialEq)]
pub struct Deck {
    pub text: String,
    pub warnings: Vec<Diagnostic>,
}

impl Format {
    /// Write the model (SI) in this format, in the consistent system `units`.
    pub fn write(&self, model: &Model, units: Preset) -> Result<Deck> {
        // The check's warnings come from the model as stored (its magnitudes are judged in SI);
        // the deck from the model expressed in `units`.
        let warnings = model.resolve()?.warnings;
        let expressed = express(model, units)?;
        let resolved = expressed.resolve()?;
        let text = (self.writer)(&expressed, &resolved, units)?;
        Ok(Deck { text, warnings })
    }

    /// Write only the named steps (in the model's order), or every step when `steps` is `None`;
    /// for runs that hold some steps only, such as a Nastran run (static or frequency steps).
    pub fn write_steps(&self, model: &Model, steps: Option<&[String]>, units: Preset) -> Result<Deck> {
        match steps {
            None => self.write(model, units),
            Some(names) => self.write(&with_steps(model, names)?, units),
        }
    }
}

/// The model, stored in SI, with every number expressed in `units`: what a deck in that system
/// writes as it is.
pub fn express(model: &Model, units: Preset) -> Result<Model> {
    let system = DisplaySystem::new(units);
    Ok(model.converted(&|q, v| system.express(q, v).ok())?)
}

/// What a preset's units are, in words, for a deck's header.
pub(crate) fn units_in_words(units: Preset) -> &'static str {
    match units {
        Preset::MKgS => "metre, kilogram, second, newton, pascal, kelvin",
        Preset::MmTonneS => "millimetre, tonne, second, newton, megapascal, degree Celsius",
        Preset::MmKgMs => "millimetre, kilogram, millisecond, kilonewton, gigapascal, degree Celsius",
        Preset::MmGMs => "millimetre, gram, millisecond, newton, megapascal, degree Celsius",
        Preset::InLbfS => "inch, lbf*s^2/in, second, pound-force, psi, degree Fahrenheit",
    }
}

/// The model with only the named steps, in its own order. Definitions the dropped steps used
/// stay (the check warns about the ones no step uses).
pub fn with_steps(model: &Model, names: &[String]) -> Result<Model> {
    if names.is_empty() {
        return Err(Error::Model(fieldcraft_model::Error::Invalid("choose at least one step to write".into())));
    }
    if let Some(n) = names.iter().find(|n| model.step(n).is_none()) {
        let shown: String = n.chars().take(fieldcraft_model::MAX_NAME_LEN).collect();
        return Err(Error::Model(fieldcraft_model::Error::NotFound(format!("step `{shown}`"))));
    }
    let mut out = model.clone();
    let dropped: Vec<String> = model.steps().iter().map(|s| s.name.clone()).filter(|s| !names.contains(s)).collect();
    for name in dropped {
        out.remove_step(&name)?;
    }
    Ok(out)
}

/// A fluid step: solved by FieldCraft's lattice Boltzmann solver, never written to a deck.
pub(crate) fn fluid_step(format: &'static str, step: &str) -> Error {
    Error::Unsupported { format, what: format!("fluid step `{step}` (solve it with solve.fluid; write the structural steps only)") }
}

/// Every writer's first check: the model has no fluid step.
pub(crate) fn refuse_fluid(format: &'static str, model: &Model) -> Result<()> {
    match model.steps().iter().find(|s| s.procedure.fluid_domain().is_some()) {
        Some(step) => Err(fluid_step(format, &step.name)),
        None => Ok(()),
    }
}

/// Every deck format, one row per writer.
pub const FORMATS: &[Format] = &[
    Format { name: "abaqus", extension: "inp", writer: abaqus::write_resolved },
    Format { name: "nastran", extension: "bdf", writer: nastran::write_resolved },
    Format { name: "lsdyna", extension: "k", writer: lsdyna::write_resolved },
];

/// The format called `name`.
pub fn format(name: &str) -> Option<&'static Format> {
    FORMATS.iter().find(|f| f.name == name)
}

/// Fail with [`Error::Unsupported`] naming every item, or succeed when there are none. Writers
/// collect everything their format cannot express and report it at once.
pub(crate) fn unsupported(format: &'static str, items: BTreeSet<String>) -> Result<()> {
    if items.is_empty() { Ok(()) } else { Err(Error::Unsupported { format, what: items.into_iter().collect::<Vec<_>>().join("; ") }) }
}

/// Heat transfer in a model, which no writer expresses yet: heat transfer and coupled
/// temperature–displacement steps, temperature boundary conditions, heat loads, temperature
/// fields and initial (temperature or stress) fields. Thermal material properties alone are left out of the decks (nothing uses them).
pub(crate) fn thermal_items(model: &Model) -> BTreeSet<String> {
    let mut out = BTreeSet::new();
    for step in model.steps().iter().filter(|s| s.procedure.is_thermal()) {
        out.insert(format!("step `{}`: {} steps (not written to decks yet)", step.name, step.procedure.name()));
    }
    for name in model.boundary_conditions().iter().filter(|(_, b)| b.is_thermal()).map(|(n, _)| n) {
        out.insert(format!("boundary condition `{name}`: prescribed temperatures (not written to decks yet)"));
    }
    for (name, load) in model.loads().iter().filter(|(_, l)| l.is_thermal() || matches!(l, Load::Temperature { .. })) {
        out.insert(format!("load `{name}`: {} loads (not written to decks yet)", load.kind()));
    }
    for (name, field) in model.initial_fields() {
        out.insert(format!("initial field `{name}`: initial {} fields (not written to decks yet)", field.kind()));
    }
    out
}

/// The concentrated loads of a set of active loads: the resolved nodal values of every force,
/// traction (degrees of freedom 1–3) and moment (4–6), summed per node and degree of freedom so
/// loads sharing a node add up; zero sums are left out.
pub fn concentrated_loads<'a>(model: &Model, r: &Resolved, active: impl IntoIterator<Item = &'a String>) -> BTreeMap<(u32, u8), f64> {
    let mut out: BTreeMap<(u32, u8), f64> = BTreeMap::new();
    for name in active {
        let first_dof = match model.load(name) {
            Some(Load::Moment { .. }) => 4,
            Some(Load::Force { .. } | Load::Traction { .. }) => 1,
            _ => continue,
        };
        for (node, v) in r.nodal_loads.get(name).map(Vec::as_slice).unwrap_or_default() {
            for (dof, c) in (first_dof..).zip(v) {
                *out.entry((*node, dof)).or_default() += c;
            }
        }
    }
    out.retain(|_, v| *v != 0.0);
    out
}

/// The consistent nodal forces (degrees of freedom 1 and 2) of the active pressures on a 2D
/// model's element edges, for formats that write them as concentrated loads: the model's edge
/// integration ([`fieldcraft_model::nodal::edge_pressure`]), a positive pressure pushing
/// into the element, summed per node and direction; zero sums are left out.
pub fn edge_pressure_loads<'a>(model: &Model, r: &Resolved, active: impl IntoIterator<Item = &'a str>) -> Result<BTreeMap<(u32, u8), f64>> {
    let mut out: BTreeMap<(u32, u8), f64> = BTreeMap::new();
    for name in active {
        let Some(Load::Pressure { set, magnitude }) = model.load(name) else { continue };
        for g in r.planar_edges(set) {
            let Some(extent) = r.planar.get(&g.element).copied() else { continue };
            let labels = model.edge_nodes(g).unwrap_or_default();
            let xyz: Vec<[f64; 3]> = labels.iter().filter_map(|n| model.node(*n)).map(|n| n.xyz).collect();
            for (node, f) in labels.iter().zip(edge_pressure(&xyz, extent, *magnitude)?) {
                for (dof, c) in (1u8..).zip(&f[..2]) {
                    *out.entry((*node, dof)).or_default() += c;
                }
            }
        }
    }
    out.retain(|_, v| *v != 0.0);
    Ok(out)
}

/// How a degree of freedom is constrained apart from the inclined symmetry planes, for choosing
/// the dependent one of a plane's equation (formats eliminate it, so it may carry no other
/// constraint).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Held {
    Free,
    /// Held at this value wherever the equation applies: it cannot depend on the others, and an
    /// equation of held ones alone is redundant when the values satisfy it.
    Always(f64),
    /// Held in some steps only, or made dependent by a coupling or tie: it cannot be dependent,
    /// nor can the equation be dropped.
    Sometimes,
}

/// A linear equation Σ cⱼ u(node, dofⱼ) = 0 of an inclined symmetry plane, dependent degree of
/// freedom first.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Equation {
    pub node: u32,
    pub terms: Vec<(u8, f64)>,
}

/// A coefficient below this fraction of its row's largest is round-off of the elimination.
const NEGLIGIBLE: f64 = 1e-12; // unit-factors: relative to the row's largest coefficient
/// A held node is on its inclined plane when the equation's residual is below this fraction of
/// the sum of its terms' magnitudes.
const ON_PLANE: f64 = 1e-9; // unit-factors: relative to the equation's terms

/// The equations of the inclined symmetry conditions `names` (others are skipped), node by node
/// in label order: n · u = 0 and, at a node in `rotating`, t₁ · θ = 0 and t₂ · θ = 0 for the
/// plane's in-plane axes ([`fieldcraft_model::Symmetry`]). The rows of one node (several planes
/// meet at an edge) are reduced against each other, so each has a dependent degree of freedom
/// of its own: the largest coefficient among those `held` says are free. A row of held degrees
/// of freedom only is left out when their values satisfy it; otherwise it is a problem, as is a
/// row with none free that cannot be left out.
pub(crate) fn symmetry_equations<'a>(
    model: &Model,
    r: &Resolved,
    names: impl IntoIterator<Item = &'a String>,
    rotating: &BTreeSet<u32>,
    held: impl Fn(u32, u8) -> Held,
    problems: &mut BTreeSet<String>,
) -> Vec<Equation> {
    let mut rows: BTreeMap<u32, Vec<(&str, [f64; 6])>> = BTreeMap::new();
    let mut seen = BTreeSet::new();
    for name in names {
        if !seen.insert(name.as_str()) {
            continue;
        }
        let Some(BoundaryCondition::Symmetry { target, plane }) = model.boundary_condition(name) else { continue };
        let (Some(n), Some([t1, t2])) = (plane.inclined(), plane.in_plane()) else { continue };
        for &node in r.target_nodes(target) {
            let at = rows.entry(node).or_default();
            at.push((name, [n[0], n[1], n[2], 0.0, 0.0, 0.0]));
            if rotating.contains(&node) || target.reference_point().is_some() {
                at.push((name, [0.0, 0.0, 0.0, t1[0], t1[1], t1[2]]));
                at.push((name, [0.0, 0.0, 0.0, t2[0], t2[1], t2[2]]));
            }
        }
    }
    let mut out = Vec::new();
    for (node, rows) in rows {
        // The equations of this node so far: (dependent index, coefficients).
        let mut done: Vec<(usize, [f64; 6])> = Vec::new();
        for (name, mut c) in rows {
            let scale = c.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            for (d, e) in &done {
                let f = c[*d] / e[*d];
                for k in 0..6 {
                    c[k] -= f * e[k];
                }
                c[*d] = 0.0;
            }
            for v in &mut c {
                if v.abs() <= NEGLIGIBLE * scale {
                    *v = 0.0;
                }
            }
            if c.iter().all(|v| *v == 0.0) {
                continue;
            }
            let dof = |k: usize| k as u8 + 1;
            let pivot = (0..6).filter(|&k| c[k] != 0.0 && held(node, dof(k)) == Held::Free).fold(None, |best: Option<usize>, k| match best {
                Some(b) if c[b].abs() >= c[k].abs() => Some(b),
                _ => Some(k),
            });
            let Some(pivot) = pivot else {
                let mut residual = 0.0;
                let mut size = 0.0;
                let mut always = true;
                for k in (0..6).filter(|&k| c[k] != 0.0) {
                    match held(node, dof(k)) {
                        Held::Always(v) => {
                            residual += c[k] * v;
                            size += (c[k] * v).abs();
                        }
                        _ => always = false,
                    }
                }
                if !always {
                    problems.insert(format!(
                        "boundary condition `{name}`: at node {node} every degree of freedom the inclined symmetry plane relates is held in some steps only or follows a coupling or tie, so none can be eliminated by its equation"
                    ));
                } else if residual.abs() > ON_PLANE * size {
                    problems.insert(format!("boundary condition `{name}`: node {node} is held off the inclined symmetry plane"));
                }
                continue;
            };
            let mut terms = vec![(dof(pivot), c[pivot])];
            terms.extend((0..6).filter(|&k| k != pivot && c[k] != 0.0).map(|k| (dof(k), c[k])));
            out.push(Equation { node, terms });
            done.push((pivot, c));
        }
    }
    out
}

/// The global degree of freedom (1–6) a spring or dashpot direction acts along, when it is one
/// (a vector or line along an axis counts: the sign does not matter to a spring).
pub(crate) fn axis_dof(d: ResolvedDirection) -> Option<u8> {
    match d {
        ResolvedDirection::Dof(k) => Some(k),
        ResolvedDirection::Vector(n) | ResolvedDirection::Line(n) => {
            let axis = n.iter().position(|c| c.abs() == 1.0)?;
            (n.iter().filter(|c| **c != 0.0).count() == 1).then_some(axis as u8 + 1)
        }
    }
}

/// `base`, or `base-2`, `base-3`, … : the first that is not in `taken` (compared without case),
/// which it is then added to.
pub(crate) fn unique_name(base: &str, taken: &mut BTreeSet<String>) -> String {
    let mut candidate = base.to_string();
    let mut k = 2u32;
    while taken.contains(&candidate.to_ascii_uppercase()) {
        candidate = format!("{base}-{k}");
        k = k.saturating_add(1);
    }
    taken.insert(candidate.to_ascii_uppercase());
    candidate
}
