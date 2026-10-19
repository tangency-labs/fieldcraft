//! LS-DYNA keyword decks for implicit linear static and eigenvalue analyses.
//!
//! The deck is written from the public LS-DYNA keyword manual (R14), in our own words and with
//! our own examples. What it writes:
//! - `*KEYWORD long=y`: every field 20 columns wide (see [`card`] for why), then `*NODE` for the
//!   mesh and the reference points;
//! - one `*PART` per section and element kind, with `*SECTION_SOLID`, `*SECTION_SHELL` or
//!   `*SECTION_BEAM` (rods as trusses) and the element formulation mapped from the section's
//!   ([`elements`]); `*HOURGLASS` where a reduced-integration formulation needs it;
//! - `*MAT_ELASTIC` for isotropic and `*MAT_ORTHOTROPIC_ELASTIC` (global material axes) for
//!   orthotropic materials; von Mises plasticity as `*MAT_PLASTIC_KINEMATIC` (bilinear: the
//!   tangent modulus E·H/(E + H) of the uniaxial curve, β = 1 isotropic or 0 kinematic, the
//!   erosion plastic strain as its failure strain) or `*MAT_PIECEWISE_LINEAR_PLASTICITY` (a
//!   hardening table as a `*DEFINE_CURVE` of yield stress against plastic strain, the erosion
//!   plastic strain as its failure strain). Johnson–Cook plasticity and damage are reported: their
//!   LS-DYNA material needs an equation of state and belongs to explicit analyses, which this
//!   writer does not write yet;
//! - `*ELEMENT_SOLID` (tetrahedra and wedges as degenerate hexahedra, 10-node tetrahedra),
//!   `*ELEMENT_SOLID_H20`, `*ELEMENT_SHELL` and `*ELEMENT_BEAM`;
//! - `*SET_NODE_LIST` and `*SET_SEGMENT` for the sets boundary conditions, loads and couplings
//!   use (with the set's name as title);
//! - couplings: kinematic as `*CONSTRAINED_NODAL_RIGID_BODY` (the reference point first in its
//!   node set), distributing as `*CONSTRAINED_INTERPOLATION` (the reference point dependent on
//!   the set's translations, equal weights);
//! - contact: a pair as `*CONTACT_AUTOMATIC_SURFACE_TO_SURFACE` between segment sets (the
//!   secondary surface first; a tie as `*CONTACT_TIED_SURFACE_TO_SURFACE` likewise), general
//!   contact as `*CONTACT_AUTOMATIC_SINGLE_SURFACE` on a segment
//!   set or on every part (`surfatyp` 5); the static and kinetic coefficients and the decay
//!   coefficient (`fs`, `fd`, `dc`: the same exponential decay law) and the penalty scale on both
//!   sides (`sfsa`, `sfsb`, a factor on LS-DYNA's own default stiffness). LS-DYNA moves nodes
//!   that start inside the other surface out, as `adjust` does; a refusal is noted in a comment;
//! - per step: `*CONTROL_IMPLICIT_GENERAL`, `*CONTROL_IMPLICIT_SOLUTION` (linear) or
//!   `*CONTROL_IMPLICIT_EIGENVALUE`, `*CONTROL_TERMINATION` (a static step is one implicit step
//!   from time 0 to 1 with every load ramped by one curve, `*DEFINE_CURVE`); a nonlinear static
//!   step runs over the same time with the same ramp (or its load path as the curve), its first
//!   increment as `dt0`, the nonlinear
//!   solver with a stiffness reformation every iteration (`ilimit` 1, full Newton) and its
//!   iteration cap as `maxref`, and `*CONTROL_IMPLICIT_AUTO` with its minimum and maximum
//!   increments (its increment cap and tolerances have no LS-DYNA counterpart and are left to
//!   LS-DYNA's defaults; a step without nonlinear geometry and hyperelastic materials in nonlinear
//!   steps are refused); boundary conditions
//!   as `*BOUNDARY_SPC_*_ID` (named) and `*BOUNDARY_PRESCRIBED_MOTION_*_ID` for non-zero
//!   displacements, a symmetry plane not normal to X, Y or Z as `*BOUNDARY_SPC_*_ID` in a
//!   `*DEFINE_COORDINATE_SYSTEM` whose x axis is the plane's normal; forces, moments and tractions as `*LOAD_NODE_SET`/`_POINT` from the
//!   resolved nodal values; pressure as `*LOAD_SEGMENT_SET`; gravity as `*LOAD_BODY_X/Y/Z`
//!   (with `*LOAD_BODY_PARTS` when it acts on a set); output as `*DATABASE_BINARY_D3PLOT`,
//!   `*DATABASE_SPCFORC` and `*DATABASE_EXTENT_BINARY`.
//!
//! Point masses are `*ELEMENT_MASS` (one mass for every direction) or `*ELEMENT_INERTIA` (the
//! inertia tensor's components, with the mass); springs and dashpots are `*ELEMENT_DISCRETE` with
//! one `*PART`, `*SECTION_DISCRETE` (torsional for a rotational degree of freedom) and
//! `*MAT_SPRING_ELASTIC` or `*MAT_DAMPER_VISCOUS` per definition, acting along the line of their
//! nodes or along a `*DEFINE_SD_ORIENTATION` vector. A spring or dashpot to the ground connects
//! to an extra node at the same place, held by `*BOUNDARY_SPC_NODE` (labels after the mesh's).
//! Element, part, section and material ids continue after the mesh's. A mass that differs by
//! direction and the three-node truss (rod3) are reported.
//!
//! Several steps are independent analyses in the model (nothing carries over), which is what
//! LS-DYNA's `*CASE` runs: each step is a case between `*CASE_BEGIN_n` and `*CASE_END_n`, run in
//! order when LS-DYNA is started with `CASE` on its command line. A single step is written
//! without cases.
//!
//! 2D models (plane stress and plane strain in the XY plane) are `*SECTION_SHELL` parts with
//! ELFORM 12 and 13 ([`elements`]); plane strain is per unit thickness in LS-DYNA, so only a
//! thickness of 1 is written; pressures on 2D element edges are the consistent nodal forces
//! (`*LOAD_NODE_POINT`).
//!
//! What LS-DYNA cannot take is an [`Error::Unsupported`] naming the item: beam sections (not yet
//! written), body forces per volume, kinematic couplings of some degrees of freedom only, boundary
//! conditions on a coupled reference point, gravity on different sets, hybrid formulations other
//! than on tet4, reduced integration of quadratic elements, materials without a density, more
//! than 10 shell integration points.

mod card;
mod elements;

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_model::{
    BoundaryCondition, ContactSurfaces, CouplingKind, Diagnostic, Dofs, Elastic, Element, ElementKind, Hardening, InitialPenetration, Load, Material,
    Model, OutputVariable, Plastic, Procedure, Resolved, ResolvedDirection, Section, SectionKind, Step, Target,
};

use crate::{Error, Result};
use card::{Deck, Field};
use card::{WIDTH, ints};
use elements::{Hourglass, Keyword};
use fieldcraft_model::units::Preset;

/// The format's name in [`Error::Unsupported`].
pub(crate) const FORMAT: &str = "LS-DYNA";

/// The positions of an element's model-order nodes in the order LS-DYNA takes them (repeats
/// for degenerate forms), for readers and checks of written decks.
pub fn node_order(kind: ElementKind) -> &'static [usize] {
    elements::order(kind)
}

/// A written deck.
#[derive(Debug, Clone, PartialEq)]
pub struct Written {
    /// The deck, ASCII.
    pub text: String,
    /// Reals that did not fit a 20-column field exactly and were rounded (to at least 13
    /// significant digits).
    pub rounded_values: usize,
    /// Number of `*CASE` analyses (0 when the model has one step and no cases are written).
    pub cases: usize,
    /// The model check's warnings (a written model has no errors).
    pub warnings: Vec<Diagnostic>,
}

/// Largest deck written, in bytes (a deck is built in memory).
pub const MAX_DECK_BYTES: usize = 1 << 31;

/// Most shell integration points LS-DYNA's Gauss rule takes.
const MAX_SHELL_POINTS: u32 = 10;

/// A resolved model that breaks its own invariants (unreachable after `resolve`).
fn broken(what: String) -> Error {
    Error::Model(fieldcraft_model::Error::Invalid(format!("inconsistent resolved model: {what}")))
}

/// One `*PART`: a section's elements of one kind.
struct Part<'a> {
    id: u32,
    name: &'a str,
    section: &'a Section,
    kind: ElementKind,
    elform: i32,
    material: u32,
    hourglass: Option<u32>,
    elements: Vec<(u32, &'a Element)>,
}

/// Every id the deck uses, decided before anything is written.
struct Plan<'a> {
    parts: Vec<Part<'a>>,
    /// Material name → MID, for the materials sections use.
    materials: BTreeMap<&'a str, u32>,
    hourglasses: BTreeMap<Hourglass, u32>,
    /// Set name → SID (node and segment sets have separate id spaces; both use this number).
    sets: BTreeMap<&'a str, u32>,
    /// Sets written as node sets, and as segment sets.
    node_sets: BTreeSet<&'a str>,
    segment_sets: BTreeSet<&'a str>,
    /// Forces and moments written on a node set: load → (value per node, SID).
    uniform: BTreeMap<&'a str, ([f64; 3], u32)>,
    /// The part set gravity acts on, when it acts on a set.
    gravity_parts: Option<Vec<u32>>,
    /// Material name → the curve id of its hardening table.
    curves: BTreeMap<&'a str, u32>,
}

/// Write the model as an LS-DYNA keyword deck.
pub fn write(model: &Model, units: Preset) -> Result<Written> {
    let expressed = crate::express(model, units)?;
    let r = expressed.resolve()?;
    write_with(&expressed, &r, units)
}

/// The writer of [`crate::FORMATS`]: the deck's text from a model expressed in `units` and
/// resolved.
pub(crate) fn write_resolved(model: &Model, r: &Resolved, units: Preset) -> Result<String> {
    write_with(model, r, units).map(|w| w.text)
}

/// `*CONTROL_UNITS`' length, time, mass and temperature words for a preset.
fn control_units(units: Preset) -> [&'static str; 4] {
    match units {
        Preset::MKgS => ["m", "sec", "kg", "k"],
        Preset::MmTonneS => ["mm", "sec", "mtrc_ton", "c"],
        Preset::MmKgMs => ["mm", "ms", "kg", "c"],
        Preset::MmGMs => ["mm", "ms", "g", "c"],
        Preset::InLbfS => ["in", "sec", "lbf*s2/in", "f"],
    }
}

fn write_with(model: &Model, r: &Resolved, units: Preset) -> Result<Written> {
    crate::refuse_fluid(FORMAT, model)?;
    let estimate = (model.nodes().len() + model.reference_points().len()) * 130 + model.elements().len() * 260;
    if estimate > MAX_DECK_BYTES {
        return Err(Error::Limit(format!("the deck would be about {estimate} bytes; at most {MAX_DECK_BYTES} are written")));
    }
    let plan = plan(model, r)?;
    let mut d = Deck { text: String::with_capacity(estimate), rounded: 0 };
    let steps = model.steps();
    let cases = if steps.len() > 1 { steps.len() } else { 0 };

    d.line("*KEYWORD long=y");
    d.comment("Written by FieldCraft: implicit linear analysis, long format (20-column fields).");
    d.comment(&format!("Units: {} ({}, consistent).", units.label(), crate::units_in_words(units)));
    d.line("*CONTROL_UNITS");
    d.header(&["length", "time", "mass", "temp"]);
    d.line(&control_units(units).iter().map(|w| format!("{w:>WIDTH$}")).collect::<String>());
    if cases > 0 {
        d.comment(&format!("{cases} steps are written as *CASE analyses: run LS-DYNA with CASE on its command line."));
    }
    nodes(&mut d, model);
    materials(&mut d, model, &plan)?;
    parts(&mut d, &plan)?;
    element_cards(&mut d, &plan)?;
    sets(&mut d, model, r, &plan)?;
    couplings(&mut d, model, r, &plan)?;
    contacts(&mut d, model, &plan)?;
    discrete(&mut d, model, r, &plan)?;
    if let Some(parts) = &plan.gravity_parts {
        d.line("*SET_PART_LIST");
        d.header(&["sid"]);
        d.card([Field::from(1u32)]);
        d.list(parts);
        d.line("*LOAD_BODY_PARTS");
        d.header(&["psid"]);
        d.card([Field::from(1u32)]);
    }
    let mut constraint_id = 0u32;
    for (i, step) in steps.iter().enumerate() {
        let n = u32::try_from(i + 1).map_err(|_| Error::Limit("too many steps".into()))?;
        if cases > 0 {
            d.line(&format!("*CASE_BEGIN_{n}"));
        }
        write_step(&mut d, model, r, &plan, step, n, &mut constraint_id)?;
        if cases > 0 {
            d.line(&format!("*CASE_END_{n}"));
        }
    }
    d.line("*END");
    if d.text.len() > MAX_DECK_BYTES {
        return Err(Error::Limit(format!("the deck is {} bytes; at most {MAX_DECK_BYTES} are written", d.text.len())));
    }
    Ok(Written { text: d.text, rounded_values: d.rounded, cases, warnings: r.warnings.clone() })
}

fn plan<'a>(model: &'a Model, r: &Resolved) -> Result<Plan<'a>> {
    let used: BTreeSet<&str> = model.sections().values().map(|s| s.material.as_str()).collect();
    let materials: BTreeMap<&str, u32> = used.into_iter().enumerate().map(|(i, m)| Ok((m, id(i)?))).collect::<Result<_>>()?;
    // Everything LS-DYNA cannot express, reported at once.
    let mut refused = crate::thermal_items(model);
    let mut parts = Vec::new();
    let mut hourglasses: BTreeMap<Hourglass, u32> = BTreeMap::new();
    for (name, sec) in model.sections() {
        let mut by_kind: BTreeMap<ElementKind, Vec<(u32, &Element)>> = BTreeMap::new();
        for e in r.sets.get(&sec.set).map(|s| s.elements.as_slice()).unwrap_or_default() {
            let el = model.element(*e).ok_or_else(|| broken(format!("element {e} of section `{name}`")))?;
            by_kind.entry(el.kind).or_default().push((*e, el));
        }
        if let SectionKind::PlaneStrain { thickness } = sec.kind
            && thickness != 1.0
        {
            refused.insert(format!("section `{name}`: a plane strain thickness of {thickness} (LS-DYNA's plane strain is per unit thickness)"));
        }
        if let SectionKind::Shell { points, .. } = sec.kind
            && points > MAX_SHELL_POINTS
        {
            refused.insert(format!(
                "section `{name}`: {points} integration points through the thickness (LS-DYNA's Gauss rule takes at most {MAX_SHELL_POINTS})"
            ));
        }
        for (kind, list) in by_kind {
            let f = match elements::elform(kind, &sec.kind, sec.formulation, name) {
                Ok(f) => f,
                Err(why) => {
                    refused.insert(why);
                    continue;
                }
            };
            let next = id(hourglasses.len())?;
            let hourglass = f.hourglass.map(|h| *hourglasses.entry(h).or_insert(next));
            parts.push(Part {
                id: id(parts.len())?,
                name,
                section: sec,
                kind,
                elform: f.elform,
                material: materials.get(sec.material.as_str()).copied().unwrap_or_default(),
                hourglass,
                elements: list,
            });
        }
    }
    let sets: BTreeMap<&str, u32> = model.sets().keys().enumerate().map(|(i, n)| Ok((n.as_str(), id(i)?))).collect::<Result<_>>()?;

    // Which sets are written, and what the steps need that LS-DYNA may not express.
    let coupled: BTreeSet<&str> = model.couplings().values().map(|c| c.reference_point.as_str()).collect();
    let mut node_sets = BTreeSet::new();
    let mut segment_sets = BTreeSet::new();
    let mut uniform = BTreeMap::new();
    for (name, c) in model.couplings() {
        if c.kind == CouplingKind::Kinematic && c.dofs != Dofs::all() {
            refused.insert(format!(
                "coupling `{name}`: a kinematic coupling of some degrees of freedom only (LS-DYNA's nodal rigid bodies release degrees of freedom reliably only for two nodes)"
            ));
        }
        if c.kind == CouplingKind::Distributing {
            node_sets.insert(c.set.as_str());
        }
    }
    // Contact surfaces are segment sets.
    for c in model.contacts().values() {
        segment_sets.extend(c.surfaces.sets());
    }
    for t in model.ties().values() {
        segment_sets.extend([t.main.as_str(), t.secondary.as_str()]);
    }
    let mut gravity_sets: BTreeSet<Option<&str>> = BTreeSet::new();
    for step in model.steps() {
        if matches!(step.procedure, Procedure::ExplicitDynamic(_)) {
            refused.insert(format!("step `{}`: explicit dynamic steps (this writer writes implicit analyses so far)", step.name));
        }
        if let Procedure::NonlinearStatic(n) = &step.procedure {
            if !n.nonlinear_geometry {
                refused.insert(format!(
                    "step `{}`: a nonlinear step without nonlinear geometry (LS-DYNA's nonlinear implicit analysis always follows the deformed shape)",
                    step.name
                ));
            }
            for material in model.hyperelastic_in_use() {
                refused.insert(format!(
                    "step `{}`: hyperelastic material `{material}` in a nonlinear step (LS-DYNA's rubber materials define their volumetric energy differently; not written yet)",
                    step.name
                ));
            }
        }
        let is_static = matches!(step.procedure, Procedure::Static {} | Procedure::NonlinearStatic(_));
        for b in &step.boundary_conditions {
            let Some(bc) = model.boundary_conditions().get(b) else { continue };
            match bc.target() {
                Target::Set(s) => {
                    node_sets.insert(s.as_str());
                }
                Target::ReferencePoint(p) if coupled.contains(p.as_str()) => {
                    refused.insert(format!(
                        "boundary condition `{b}`: reference point `{p}` is coupled, and LS-DYNA takes no constraint on a node of a nodal rigid body or on an interpolation's dependent node; constrain the coupled set instead"
                    ));
                }
                Target::ReferencePoint(_) => {}
            }
        }
        if !is_static {
            continue;
        }
        for l in &step.loads {
            let Some(load) = model.loads().get(l) else { continue };
            match load {
                Load::Force { target, .. } | Load::Moment { target, .. } => {
                    if let Target::Set(s) = target
                        && let Some(v) = uniform_on_set(r, l, s)
                        && let Some(sid) = sets.get(s.as_str())
                    {
                        node_sets.insert(s.as_str());
                        uniform.insert(l.as_str(), (v, *sid));
                    }
                }
                Load::Traction { .. } => {}
                Load::Pressure { .. } if r.is_planar() => {}
                Load::Pressure { set, .. } => {
                    segment_sets.insert(set.as_str());
                }
                Load::Centrifugal { .. } => {
                    refused.insert(format!("load `{l}`: a centrifugal load (not written to LS-DYNA decks yet)"));
                }
                Load::Gravity { set, .. } => {
                    gravity_sets.insert(set.as_deref());
                }
                Load::BodyForce { .. } => {
                    refused.insert(format!("load `{l}`: a body force per unit volume (LS-DYNA's body loads are accelerations; use gravity)"));
                }
                // In thermal_items.
                Load::SurfaceFlux { .. } | Load::Convection { .. } | Load::HeatSource { .. } | Load::Temperature { .. } => {}
            }
        }
    }
    let gravity_parts = match gravity_sets.into_iter().collect::<Vec<_>>().as_slice() {
        [] | [None] => None,
        [Some(set)] => parts_of_set(r, &parts, set).map_err(|why| refused.insert(why)).ok(),
        _ => {
            refused.insert("gravity on different sets (LS-DYNA applies every body load to one part subset per deck)".into());
            None
        }
    };
    if r.is_planar() {
        for name in model.sections().iter().filter(|(_, s)| matches!(s.kind, SectionKind::Rod { .. })).map(|(n, _)| n) {
            refused.insert(format!("section `{name}`: trusses in a 2D model (LS-DYNA's 2D elements take no truss)"));
        }
    }
    for m in &r.masses {
        if m.mass[0] != m.mass[1] || m.mass[0] != m.mass[2] {
            refused.insert(format!("point mass `{}`: a mass that differs by direction (LS-DYNA's mass elements take one mass)", m.name));
        }
    }
    let mut curves = BTreeMap::new();
    for name in materials.keys() {
        let Some(m) = model.material(name) else { continue };
        if m.density.is_none() {
            refused.insert(format!("material `{name}` has no density, and LS-DYNA requires a mass density on every material"));
        }
        if matches!(m.plastic, Some(Plastic::JohnsonCook { .. })) || m.damage.is_some() {
            refused.insert(format!(
                "material `{name}`: Johnson–Cook plasticity and damage (LS-DYNA's Johnson–Cook material needs an equation of state and is for explicit analyses, which this writer does not write yet)"
            ));
        }
        if matches!(m.plastic, Some(Plastic::Table { .. })) {
            // Curve ids after the steps' load ramps (ids 1 to the number of steps).
            let lcid = u32::try_from(model.steps().len() + 1 + curves.len()).map_err(|_| Error::Limit("too many curves".into()))?;
            curves.insert(*name, lcid);
        }
    }
    crate::unsupported(FORMAT, refused)?;
    Ok(Plan { parts, materials, hourglasses, sets, node_sets, segment_sets, uniform, gravity_parts, curves })
}

/// A 1-based id from a 0-based index.
fn id(i: usize) -> Result<u32> {
    u32::try_from(i + 1).map_err(|_| Error::Limit("more than 4294967295 items of one kind".into()))
}

/// The parts that make up a set's elements exactly, or why not.
fn parts_of_set(r: &Resolved, parts: &[Part], set: &str) -> std::result::Result<Vec<u32>, String> {
    let elements: BTreeSet<u32> = r.sets.get(set).map(|s| s.elements.iter().copied().collect()).unwrap_or_default();
    let mut ids = Vec::new();
    let mut covered = 0;
    for p in parts {
        let inside = p.elements.iter().filter(|(e, _)| elements.contains(e)).count();
        if inside == p.elements.len() {
            ids.push(p.id);
            covered += inside;
        } else if inside > 0 {
            return Err(format!(
                "gravity on set `{set}`, which holds part of section `{}`'s elements (LS-DYNA applies body loads to whole parts)",
                p.name
            ));
        }
    }
    if covered != elements.len() || ids.is_empty() {
        return Err(format!("gravity on set `{set}`, whose elements are not whole sections"));
    }
    Ok(ids)
}

/// The value every node of a force or moment gets, when it is one value on exactly the set's
/// nodes (then the load is written on the node set).
fn uniform_on_set(r: &Resolved, load: &str, set: &str) -> Option<[f64; 3]> {
    let values = r.nodal_loads.get(load)?;
    let first = values.first()?.1;
    let nodes = &r.sets.get(set)?.nodes;
    (values.len() == nodes.len() && values.iter().zip(nodes).all(|((n, v), m)| n == m && *v == first)).then_some(first)
}

fn nodes(d: &mut Deck, model: &Model) {
    // The mesh nodes are sorted; the few reference points are merged in by label.
    let mut points: Vec<(u32, [f64; 3])> = model.reference_points().values().map(|p| (p.label, p.xyz)).collect();
    points.sort_unstable_by_key(|p| p.0);
    let mut points = points.into_iter().peekable();
    d.line("*NODE");
    d.header(&["nid", "x", "y", "z"]);
    let card = |d: &mut Deck, label: u32, [x, y, z]: [f64; 3]| d.card([Field::from(label), Field::from(x), Field::from(y), Field::from(z)]);
    for (label, n) in model.nodes() {
        while let Some((l, xyz)) = points.next_if(|p| p.0 < *label) {
            card(d, l, xyz);
        }
        card(d, *label, n.xyz);
    }
    for (l, xyz) in points {
        card(d, l, xyz);
    }
}

fn materials(d: &mut Deck, model: &Model, plan: &Plan) -> Result<()> {
    for (name, &mid) in &plan.materials {
        let m = model.material(name).ok_or_else(|| broken(format!("material `{name}`")))?;
        let rho = m.density.unwrap_or_default();
        match m.elastic {
            Some(Elastic::Isotropic { young, poisson }) if m.plastic.is_some() => plastic_material(d, name, mid, rho, young, poisson, m, plan)?,
            Some(Elastic::Isotropic { young, poisson }) => {
                d.line("*MAT_ELASTIC_TITLE");
                d.line(name);
                d.header(&["mid", "ro", "e", "pr"]);
                d.card([Field::from(mid), Field::from(rho), Field::from(young), Field::from(poisson)]);
            }
            Some(Elastic::EngineeringConstants { young: [e1, e2, e3], poisson: [n12, n13, n23], shear: [g12, g13, g23] }) => {
                // LS-DYNA takes the minor ratios: ν_ba = ν_ab·E_b/E_a.
                d.line("*MAT_ORTHOTROPIC_ELASTIC_TITLE");
                d.line(name);
                d.header(&["mid", "ro", "ea", "eb", "ec", "prba", "prca", "prcb"]);
                let mut card = vec![Field::from(mid)];
                card.extend([rho, e1, e2, e3, n12 * e2 / e1, n13 * e3 / e1, n23 * e3 / e2].map(Field::from));
                d.card(card);
                d.header(&["gab", "gbc", "gca", "aopt"]);
                d.card([g12, g23, g13, 2.0].map(Field::from));
                d.header(&["xp", "yp", "zp", "a1", "a2", "a3"]);
                d.card([Field::Blank, Field::Blank, Field::Blank, Field::from(1.0), Field::from(0.0), Field::from(0.0)]);
                d.header(&["v1", "v2", "v3", "d1", "d2", "d3"]);
                d.card([Field::Blank, Field::Blank, Field::Blank, Field::from(0.0), Field::from(1.0), Field::from(0.0)]);
            }
            None => return Err(broken(format!("material `{name}` has no elastic constants"))),
        }
    }
    Ok(())
}

/// The plastic materials: `*MAT_PLASTIC_KINEMATIC` or `*MAT_PIECEWISE_LINEAR_PLASTICITY` with its
/// curve, in place of `*MAT_ELASTIC`.
fn plastic_material(d: &mut Deck, name: &str, mid: u32, rho: f64, young: f64, poisson: f64, m: &Material, plan: &Plan) -> Result<()> {
    let fail = m.erosion.map_or(Field::Blank, |e| Field::from(e.plastic_strain));
    match &m.plastic {
        Some(Plastic::Bilinear { yield_stress, hardening_modulus, hardening }) => {
            d.line("*MAT_PLASTIC_KINEMATIC_TITLE");
            d.line(name);
            d.header(&["mid", "ro", "e", "pr", "sigy", "etan", "beta"]);
            let tangent = young * hardening_modulus / (young + hardening_modulus);
            let beta = if *hardening == Hardening::Kinematic { 0.0 } else { 1.0 };
            d.card([
                Field::from(mid),
                Field::from(rho),
                Field::from(young),
                Field::from(poisson),
                Field::from(*yield_stress),
                Field::from(tangent),
                Field::from(beta),
            ]);
            d.header(&["src", "srp", "fs"]);
            d.card([Field::Blank, Field::Blank, fail]);
        }
        Some(Plastic::Table { points }) => {
            let lcid = plan.curves.get(name).copied().ok_or_else(|| broken(format!("material `{name}` has no curve")))?;
            d.line("*MAT_PIECEWISE_LINEAR_PLASTICITY_TITLE");
            d.line(name);
            d.header(&["mid", "ro", "e", "pr", "sigy", "etan", "fail"]);
            let first = m.plastic.as_ref().map_or(0.0, Plastic::initial_yield);
            d.card([Field::from(mid), Field::from(rho), Field::from(young), Field::from(poisson), Field::from(first), Field::Blank, fail]);
            d.header(&["c", "p", "lcss"]);
            d.card([Field::Blank, Field::Blank, Field::from(lcid)]);
            d.header(&["eps1"]);
            d.card([Field::Blank]);
            d.header(&["es1"]);
            d.card([Field::Blank]);
            d.line("*DEFINE_CURVE_TITLE");
            d.line(&format!("{name} hardening"));
            d.header(&["lcid"]);
            d.card([Field::from(lcid)]);
            d.header(&["a1", "o1"]);
            for p in points {
                d.card([Field::from(p.strain), Field::from(p.stress)]);
            }
        }
        _ => return Err(broken(format!("material `{name}`: a plasticity LS-DYNA cannot take"))),
    }
    Ok(())
}

fn parts(d: &mut Deck, plan: &Plan) -> Result<()> {
    for p in &plan.parts {
        match p.section.kind {
            SectionKind::Solid {} => {
                d.line("*SECTION_SOLID");
                d.header(&["secid", "elform"]);
                d.card([Field::from(p.id), Field::from(p.elform)]);
            }
            SectionKind::Shell { thickness, offset, points } => {
                d.line("*SECTION_SHELL");
                d.header(&["secid", "elform", "shrf", "nip"]);
                d.card([Field::from(p.id), Field::from(p.elform), Field::from(5.0 / 6.0), Field::from(points)]);
                d.header(&["t1", "t2", "t3", "t4", "nloc"]);
                // The model's offset is the reference (node) surface's distance from the
                // mid-surface along the normal, as a fraction of the thickness; NLOC puts the
                // mid-surface at −NLOC·t/2 from the nodes.
                let nloc = if offset == 0.0 { Field::Blank } else { Field::from(2.0 * offset) };
                d.card([Field::from(thickness), Field::Blank, Field::Blank, Field::Blank, nloc]);
            }
            SectionKind::Rod { area } => {
                d.line("*SECTION_BEAM");
                d.header(&["secid", "elform"]);
                d.card([Field::from(p.id), Field::from(p.elform)]);
                d.header(&["a"]);
                d.card([Field::from(area)]);
            }
            SectionKind::Beam { .. } => return Err(broken(format!("section `{}`: beams have no LS-DYNA formulation", p.name))),
            SectionKind::PlaneStress { thickness } | SectionKind::PlaneStrain { thickness } => {
                d.line("*SECTION_SHELL");
                d.header(&["secid", "elform"]);
                d.card([Field::from(p.id), Field::from(p.elform)]);
                d.header(&["t1", "t2", "t3", "t4"]);
                d.card([Field::from(thickness), Field::from(thickness), Field::from(thickness), Field::from(thickness)]);
            }
            SectionKind::Axisymmetric {} => return Err(broken(format!("section `{}`: axisymmetric sections have no LS-DYNA formulation", p.name))),
        }
    }
    for (h, hgid) in &plan.hourglasses {
        d.line("*HOURGLASS");
        d.header(&["hgid", "ihq", "qm"]);
        match h {
            Hourglass::Solid => d.card([Field::from(*hgid), Field::from(6), Field::from(1.0)]),
            Hourglass::Shell => d.card([Field::from(*hgid), Field::from(4)]),
        }
    }
    let several: BTreeMap<&str, usize> = plan.parts.iter().fold(BTreeMap::new(), |mut m, p| {
        *m.entry(p.name).or_default() += 1;
        m
    });
    for p in &plan.parts {
        d.line("*PART");
        if several.get(p.name).copied().unwrap_or_default() > 1 {
            d.line(&format!("{} {}", p.name, p.kind.name()));
        } else {
            d.line(p.name);
        }
        d.header(&["pid", "secid", "mid", "eosid", "hgid"]);
        let hg = p.hourglass.map_or(Field::Blank, Field::from);
        d.card([Field::from(p.id), Field::from(p.id), Field::from(p.material), Field::Blank, hg]);
    }
    Ok(())
}

fn element_cards(d: &mut Deck, plan: &Plan) -> Result<()> {
    let mut by_keyword: BTreeMap<Keyword, Vec<(u32, u32, &Element)>> = BTreeMap::new();
    for p in &plan.parts {
        by_keyword.entry(Keyword::of(p.kind)).or_default().extend(p.elements.iter().map(|(label, el)| (*label, p.id, *el)));
    }
    for (keyword, mut elements) in by_keyword {
        elements.sort_unstable_by_key(|e| e.0);
        d.line(keyword.name());
        match keyword {
            Keyword::Shell | Keyword::Beam => d.header(&["eid", "pid", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8"]),
            Keyword::Solid | Keyword::SolidH20 => {
                d.header(&["eid", "pid"]);
                d.header(&["n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8", "n9", "n10"]);
            }
        }
        for (label, pid, el) in elements {
            let nodes = elements::connectivity(el.kind, &el.nodes).ok_or_else(|| broken(format!("element {label} has {} nodes", el.nodes.len())))?;
            match keyword {
                Keyword::Shell | Keyword::Beam => d.card([Field::from(label), Field::from(pid)].into_iter().chain(ints(&nodes))),
                Keyword::Solid | Keyword::SolidH20 => {
                    d.card([Field::from(label), Field::from(pid)]);
                    for chunk in nodes.chunks(10) {
                        d.card(ints(chunk));
                    }
                }
            }
        }
    }
    Ok(())
}

fn sets(d: &mut Deck, model: &Model, r: &Resolved, plan: &Plan) -> Result<()> {
    for name in &plan.node_sets {
        let (Some(sid), Some(set)) = (plan.sets.get(name), r.sets.get(*name)) else { continue };
        d.line("*SET_NODE_LIST_TITLE");
        d.line(name);
        d.header(&["sid"]);
        d.card([Field::from(*sid)]);
        d.list(&set.nodes);
    }
    for name in &plan.segment_sets {
        let (Some(sid), Some(set)) = (plan.sets.get(name), r.sets.get(*name)) else { continue };
        d.line("*SET_SEGMENT_TITLE");
        d.line(name);
        d.header(&["sid"]);
        d.card([Field::from(*sid)]);
        d.header(&["n1", "n2", "n3", "n4"]);
        for f in &set.faces {
            let nodes = model.face_nodes(*f).ok_or_else(|| broken(format!("element {} face {}", f.element, f.face)))?;
            // Corner nodes only: LS-DYNA recognises the quadratic face from its corners and
            // spreads the pressure over the mid-side nodes too.
            let corners = match nodes.len() {
                3 | 6 => nodes.get(..3).zip(nodes.get(2..3)),
                _ => nodes.get(..4).zip(Some(&[][..])),
            };
            let (first, repeat) = corners.ok_or_else(|| broken(format!("element {} face {} has {} nodes", f.element, f.face, nodes.len())))?;
            d.card(ints(first).chain(ints(repeat)));
        }
    }
    Ok(())
}

fn couplings(d: &mut Deck, model: &Model, r: &Resolved, plan: &Plan) -> Result<()> {
    let mut next_set = id(model.sets().len())?;
    let mut next_part = id(plan.parts.len())?;
    let mut next_interpolation = 1u32;
    for (name, c) in model.couplings() {
        let rp =
            r.reference_point_nodes.get(&c.reference_point).copied().ok_or_else(|| broken(format!("coupling `{name}`: no reference point node")))?;
        let nodes = r.sets.get(&c.set).map(|s| s.nodes.as_slice()).unwrap_or_default();
        match c.kind {
            CouplingKind::Kinematic => {
                // The reference point leads: it is the first node of the rigid body's set.
                d.line("*SET_NODE_LIST_TITLE");
                d.line(&format!("{name} (rigid body)"));
                d.header(&["sid"]);
                d.card([Field::from(next_set)]);
                let members: Vec<u32> = std::iter::once(rp).chain(nodes.iter().copied().filter(|n| *n != rp)).collect();
                d.list(&members);
                d.line("*CONSTRAINED_NODAL_RIGID_BODY_TITLE");
                d.line(name);
                d.header(&["pid", "cid", "nsid"]);
                d.card([Field::from(next_part), Field::Blank, Field::from(next_set)]);
                next_set = next_set.checked_add(1).ok_or_else(|| Error::Limit("too many node sets".into()))?;
                next_part = next_part.checked_add(1).ok_or_else(|| Error::Limit("too many parts".into()))?;
            }
            CouplingKind::Distributing => {
                let sid = plan.sets.get(c.set.as_str()).copied().ok_or_else(|| broken(format!("coupling `{name}`: set `{}`", c.set)))?;
                let ddof: i64 = c.dofs.list().iter().fold(0, |acc, d| acc * 10 + i64::from(*d));
                d.comment(&format!("coupling `{name}`"));
                d.line("*CONSTRAINED_INTERPOLATION");
                d.header(&["icid", "dnid", "ddof", "cidd", "ityp", "idnsw", "fgm"]);
                // FGM = 1: the dependent node is not part of the mesh.
                d.card([
                    Field::from(next_interpolation),
                    Field::from(rp),
                    Field::Int(ddof),
                    Field::Blank,
                    Field::from(1),
                    Field::Blank,
                    Field::from(1),
                ]);
                d.header(&["inid", "idof"]);
                d.card([Field::from(sid), Field::from(123)]);
                next_interpolation += 1;
            }
        }
    }
    Ok(())
}

/// Contacts (see the module docs): a pair as an automatic surface-to-surface contact between the
/// segment sets, general contact as an automatic single-surface contact.
fn contacts(d: &mut Deck, model: &Model, plan: &Plan) -> Result<()> {
    for (name, c) in model.contacts() {
        let sid = |set: &str| plan.sets.get(set).copied().ok_or_else(|| broken(format!("contact `{name}`: set `{set}`")));
        d.comment(&format!("contact `{name}`"));
        if c.initial_penetration == InitialPenetration::Refuse {
            d.comment("FieldCraft refuses nodes that start inside the other surface; LS-DYNA moves them out (as `adjust` does)");
        }
        match &c.surfaces {
            ContactSurfaces::Pair { main, secondary } => {
                d.line("*CONTACT_AUTOMATIC_SURFACE_TO_SURFACE");
                d.header(&["surfa", "surfb", "surfatyp", "surfbtyp"]);
                d.card([Field::from(sid(secondary)?), Field::from(sid(main)?), Field::from(0), Field::from(0)]);
            }
            ContactSurfaces::General { surface: Some(s) } => {
                d.line("*CONTACT_AUTOMATIC_SINGLE_SURFACE");
                d.header(&["surfa", "surfb", "surfatyp"]);
                d.card([Field::from(sid(s)?), Field::Blank, Field::from(0)]);
            }
            ContactSurfaces::General { surface: None } => {
                // SURFATYP 5: every part.
                d.line("*CONTACT_AUTOMATIC_SINGLE_SURFACE");
                d.header(&["surfa", "surfb", "surfatyp"]);
                d.card([Field::from(0), Field::Blank, Field::from(5)]);
            }
        }
        let f = c.friction;
        d.header(&["fs", "fd", "dc"]);
        d.card([Field::from(f.static_coefficient), Field::from(f.kinetic_coefficient()), Field::from(f.decay)]);
        d.header(&["sfsa", "sfsb"]);
        d.card([Field::from(c.penalty_scale), Field::from(c.penalty_scale)]);
    }
    // Ties: a tied surface-to-surface contact between the segment sets, secondary first.
    for (name, t) in model.ties() {
        let sid = |set: &str| plan.sets.get(set).copied().ok_or_else(|| broken(format!("tie `{name}`: set `{set}`")));
        d.comment(&format!("tie `{name}`"));
        if t.position_tolerance.is_some() {
            d.comment("FieldCraft's position tolerance: LS-DYNA ties the nodes within its own search distance");
        }
        d.line("*CONTACT_TIED_SURFACE_TO_SURFACE");
        d.header(&["surfa", "surfb", "surfatyp", "surfbtyp"]);
        d.card([Field::from(sid(&t.secondary)?), Field::from(sid(&t.main)?), Field::from(0), Field::from(0)]);
        d.header(&["fs", "fd", "dc"]);
        d.card([Field::from(0.0), Field::from(0.0), Field::from(0.0)]);
        d.header(&["sfsa", "sfsb"]);
        d.card([Field::from(1.0), Field::from(1.0)]);
    }
    Ok(())
}

/// Point masses, springs and dashpots (see the module docs).
fn discrete(d: &mut Deck, model: &Model, r: &Resolved, plan: &Plan) -> Result<()> {
    if r.masses.is_empty() && r.springs.is_empty() {
        return Ok(());
    }
    let limit = |what: &str| Error::Limit(format!("no {what} id left for point masses, springs and dashpots"));
    let mut eid = model.elements().keys().next_back().copied().unwrap_or(0);
    let mut next_eid = || -> Result<u32> {
        eid = eid.checked_add(1).filter(|e| *e <= fieldcraft_model::MAX_LABEL).ok_or_else(|| limit("element"))?;
        Ok(eid)
    };
    d.comment("Point masses");
    for m in &r.masses {
        let mass = m.mass[0];
        if m.inertia.iter().all(|c| *c == 0.0) {
            d.line("*ELEMENT_MASS");
            d.header(&["eid", "nid", "mass"]);
            d.card([Field::from(next_eid()?), Field::from(m.node), Field::from(mass)]);
        } else {
            let [i11, i22, i33, i12, i13, i23] = m.inertia;
            d.line("*ELEMENT_INERTIA");
            d.header(&["eid", "nid", "csid"]);
            d.card([Field::from(next_eid()?), Field::from(m.node)]);
            d.header(&["ixx", "ixy", "ixz", "iyy", "iyz", "izz", "mass"]);
            d.card([i11, i12, i13, i22, i23, i33, mass].map(Field::from));
        }
    }
    if r.springs.is_empty() {
        return Ok(());
    }
    // Ids after the mesh's parts, the rigid bodies' parts and the materials.
    let kinematic = model.couplings().values().filter(|c| c.kind == CouplingKind::Kinematic).count();
    let mut pid = u32::try_from(plan.parts.len() + kinematic).map_err(|_| limit("part"))?;
    let mut mid = u32::try_from(plan.materials.len()).map_err(|_| limit("material"))?;
    let mut vid = 0u32;
    // Ground nodes: labels after the mesh's, skipping the reference points'.
    let taken: BTreeSet<u32> = model.reference_points().values().map(|p| p.label).collect();
    let mut node = model.nodes().keys().next_back().copied().unwrap_or(0);
    let mut ground: Vec<(u32, [f64; 3])> = Vec::new();
    let mut elements: Vec<Vec<Field>> = Vec::new();
    let mut current: Option<&str> = None;
    let mut vid_of_part = 0u32;
    for s in &r.springs {
        let rotational = s.direction.is_rotational();
        if current != Some(s.name.as_str()) {
            current = Some(&s.name);
            pid = pid.checked_add(1).ok_or_else(|| limit("part"))?;
            mid = mid.checked_add(1).ok_or_else(|| limit("material"))?;
            d.comment(&format!("{} `{}`", s.kind(), s.name));
            if s.dashpot {
                d.line("*MAT_DAMPER_VISCOUS_TITLE");
                d.line(&s.name);
                d.header(&["mid", "dc"]);
            } else {
                d.line("*MAT_SPRING_ELASTIC_TITLE");
                d.line(&s.name);
                d.header(&["mid", "k"]);
            }
            d.card([Field::from(mid), Field::from(s.value)]);
            d.line("*SECTION_DISCRETE");
            d.header(&["secid", "dro"]);
            d.card([Field::from(pid), Field::from(i32::from(rotational))]);
            d.line("*PART");
            d.line(&s.name);
            d.header(&["pid", "secid", "mid"]);
            d.card([Field::from(pid), Field::from(pid), Field::from(mid)]);
            let axis = (!matches!(s.direction, ResolvedDirection::Line(_))).then(|| s.direction.axis());
            let v = match axis {
                Some([x, y, z]) => {
                    vid += 1;
                    d.line("*DEFINE_SD_ORIENTATION");
                    d.header(&["vid", "iop", "xt", "yt", "zt"]);
                    d.card([Field::from(vid), Field::from(0), Field::from(x), Field::from(y), Field::from(z)]);
                    vid
                }
                None => 0,
            };
            vid_of_part = v;
        }
        let second = match s.second {
            Some(b) => b,
            None => {
                node = node.checked_add(1).ok_or_else(|| limit("node"))?;
                while taken.contains(&node) {
                    node = node.checked_add(1).ok_or_else(|| limit("node"))?;
                }
                if node > fieldcraft_model::MAX_LABEL {
                    return Err(limit("node"));
                }
                ground.push((node, model.label_xyz(s.first).ok_or_else(|| broken(format!("node {}", s.first)))?));
                node
            }
        };
        elements.push(vec![Field::from(next_eid()?), Field::from(pid), Field::from(s.first), Field::from(second), Field::from(vid_of_part)]);
    }
    d.line("*ELEMENT_DISCRETE");
    d.header(&["eid", "pid", "n1", "n2", "vid"]);
    for e in elements {
        d.card(e);
    }
    if !ground.is_empty() {
        d.comment("The ground ends of springs and dashpots, held in every direction");
        d.line("*NODE");
        d.header(&["nid", "x", "y", "z"]);
        for (l, [x, y, z]) in &ground {
            d.card([Field::from(*l), Field::from(*x), Field::from(*y), Field::from(*z)]);
        }
        d.line("*BOUNDARY_SPC_NODE");
        d.header(&["nid", "cid", "dofx", "dofy", "dofz", "dofrx", "dofry", "dofrz"]);
        for (l, _) in &ground {
            d.card([Field::from(*l), Field::Blank].into_iter().chain([1; 6].map(Field::from)));
        }
    }
    Ok(())
}

/// LS-DYNA's degree-of-freedom code for prescribed motion and nodal loads: 1–3 translations,
/// 5–7 rotations.
fn dof_code(dof: u8) -> i32 {
    if dof <= 3 { i32::from(dof) } else { i32::from(dof) + 1 }
}

fn write_step(d: &mut Deck, model: &Model, r: &Resolved, plan: &Plan, step: &Step, n: u32, constraint_id: &mut u32) -> Result<()> {
    d.line("*TITLE");
    d.line(&step.name);
    let is_static = matches!(step.procedure, Procedure::Static {} | Procedure::NonlinearStatic(_));
    let first_increment = if let Procedure::NonlinearStatic(n) = &step.procedure { n.initial_increment } else { 1.0 };
    d.line("*CONTROL_IMPLICIT_GENERAL");
    d.header(&["imflag", "dt0"]);
    d.card([Field::from(1), Field::from(first_increment)]);
    match &step.procedure {
        Procedure::Fluid { .. } => return Err(crate::fluid_step(FORMAT, &step.name)),
        Procedure::Static {} => {
            d.line("*CONTROL_IMPLICIT_SOLUTION");
            d.header(&["nsolvr"]);
            d.card([Field::from(1)]);
        }
        Procedure::NonlinearStatic(n) => {
            if n.arc_length.is_some() {
                return Err(Error::Unsupported { format: FORMAT, what: format!("step `{}`: arc-length control", step.name) });
            }
            if n.imperfection.is_some() {
                return Err(Error::Unsupported {
                    format: FORMAT,
                    what: format!("step `{}`: an imperfection from a buckling mode (FieldCraft computes the mode when it solves)", step.name),
                });
            }
            // Nonlinear with a stiffness reformation every iteration (full Newton), at most
            // `max_iterations` per increment; LS-DYNA's own convergence norms and line search.
            d.line("*CONTROL_IMPLICIT_SOLUTION");
            d.header(&["nsolvr", "ilimit", "maxref"]);
            d.card([Field::from(2), Field::from(1), Field::from(n.max_iterations)]);
            d.line("*CONTROL_IMPLICIT_AUTO");
            d.header(&["iauto", "iteopt", "itewin", "dtmin", "dtmax"]);
            d.card([Field::from(1), Field::Blank, Field::Blank, Field::from(n.min_increment), Field::from(n.max_increment)]);
        }
        &Procedure::Frequency { modes, min_hz, max_hz } => {
            d.line("*CONTROL_IMPLICIT_EIGENVALUE");
            d.header(&["neig", "center", "lflag", "lftend", "rflag", "rhtend"]);
            let bound = |b: Option<f64>| match b {
                Some(v) => [Field::from(1), Field::from(v)],
                None => [Field::Blank, Field::Blank],
            };
            let mut card = vec![Field::from(modes), Field::Blank];
            card.extend(bound(min_hz));
            card.extend(bound(max_hz));
            d.card(card);
            if step.output.contains(&OutputVariable::Stress) || step.output.contains(&OutputVariable::Strain) {
                d.header(&["isolid", "ibeam", "ishell", "itshell", "mstres"]);
                d.card([Field::Blank, Field::Blank, Field::Blank, Field::Blank, Field::from(1)]);
            }
        }
        Procedure::ExplicitDynamic(_) | Procedure::HeatTransfer { .. } | Procedure::CoupledTemperatureDisplacement {} | Procedure::Buckle { .. } => {
            return Err(Error::Unsupported { format: FORMAT, what: format!("step `{}`: {} steps", step.name, step.procedure.name()) });
        }
    }
    d.line("*CONTROL_TERMINATION");
    d.header(&["endtim"]);
    d.card([Field::from(1.0)]);

    // Every load and prescribed motion of a static step follows one ramp from 0 to 1.
    let prescribes = step
        .boundary_conditions
        .iter()
        .filter_map(|b| model.boundary_conditions().get(b))
        .any(|b| (1..=6).any(|dof| b.value(dof).is_some_and(|v| v != 0.0)));
    let loads: Vec<(&String, &Load)> = step.loads.iter().filter_map(|l| model.loads().get_key_value(l)).collect();
    if is_static && (prescribes || !loads.is_empty()) {
        d.line("*DEFINE_CURVE_TITLE");
        d.line(&format!("{} ramp", step.name));
        d.header(&["lcid"]);
        d.card([Field::from(n)]);
        d.header(&["a1", "o1"]);
        let path = match &step.procedure {
            Procedure::NonlinearStatic(s) => s.load_path.clone(),
            _ => None,
        };
        for [t, f] in path.unwrap_or_else(|| vec![[0.0, 0.0], [1.0, 1.0]]) {
            d.card([Field::from(t), Field::from(f)]);
        }
    }

    for b in &step.boundary_conditions {
        let Some(bc) = model.boundary_conditions().get(b) else { continue };
        boundary_condition(d, r, plan, b, bc, is_static, n, constraint_id)?;
    }
    if is_static {
        let mut gravity = [0.0f64; 3];
        let mut pointwise = Vec::new();
        let mut edge_pressures = Vec::new();
        for (name, load) in &loads {
            match load {
                Load::Force { .. } | Load::Moment { .. } if plan.uniform.contains_key(name.as_str()) => {
                    on_node_set(d, plan, name, matches!(load, Load::Moment { .. }), n);
                }
                Load::Force { .. } | Load::Moment { .. } | Load::Traction { .. } => pointwise.push(*name),
                Load::Pressure { .. } if r.is_planar() => edge_pressures.push(name.as_str()),
                Load::Pressure { set, magnitude } => {
                    let sid = plan.sets.get(set.as_str()).copied().unwrap_or_default();
                    d.comment(&format!("load `{name}`"));
                    d.line("*LOAD_SEGMENT_SET");
                    d.header(&["ssid", "lcid", "sf"]);
                    d.card([Field::from(sid), Field::from(n), Field::from(*magnitude)]);
                }
                Load::Gravity { acceleration, .. } => {
                    for (g, a) in gravity.iter_mut().zip(acceleration) {
                        *g += a;
                    }
                }
                Load::BodyForce { .. } => return Err(broken(format!("load `{name}`: a body force"))),
                Load::Centrifugal { .. } => return Err(broken(format!("load `{name}`: a centrifugal load"))),
                Load::SurfaceFlux { .. } | Load::Convection { .. } | Load::HeatSource { .. } | Load::Temperature { .. } => {
                    return Err(broken(format!("load `{name}`: a {} load", load.kind())));
                }
            }
        }
        // The other concentrated loads node by node, summed where they share a node.
        let mut point = crate::concentrated_loads(model, r, pointwise);
        for (key, v) in crate::edge_pressure_loads(model, r, edge_pressures)? {
            *point.entry(key).or_default() += v;
        }
        point.retain(|_, v| *v != 0.0);
        if !point.is_empty() {
            d.line("*LOAD_NODE_POINT");
            d.header(&["nid", "dof", "lcid", "sf"]);
            for ((node, dof), v) in point {
                d.card([Field::from(node), Field::from(dof_code(dof)), Field::from(n), Field::from(v)]);
            }
        }
        // A base acceleration loads the model against its direction, so gravity g is written
        // with scale factor −g.
        for (axis, g) in ["X", "Y", "Z"].iter().zip(gravity) {
            if g != 0.0 {
                d.line(&format!("*LOAD_BODY_{axis}"));
                d.header(&["lcid", "sf"]);
                d.card([Field::from(n), Field::from(-g)]);
            }
        }
        d.line("*DATABASE_BINARY_D3PLOT");
        d.header(&["dt"]);
        d.card([Field::from(1.0)]);
        if step.output.contains(&OutputVariable::ReactionForce) {
            d.line("*DATABASE_SPCFORC");
            d.header(&["dt"]);
            d.card([Field::from(1.0)]);
        }
        if step.output.contains(&OutputVariable::Strain) {
            d.line("*DATABASE_EXTENT_BINARY");
            d.header(&["neiph", "neips", "maxint", "strflg"]);
            d.card([Field::Blank, Field::Blank, Field::Blank, Field::from(1)]);
        }
    } else if !loads.is_empty() {
        d.comment("Loads have no effect on natural frequencies and are not written in this step.");
    }
    Ok(())
}

fn boundary_condition(
    d: &mut Deck,
    r: &Resolved,
    plan: &Plan,
    name: &str,
    bc: &BoundaryCondition,
    is_static: bool,
    lcid: u32,
    constraint_id: &mut u32,
) -> Result<()> {
    let (entity, suffix) = match bc.target() {
        Target::Set(s) => (plan.sets.get(s.as_str()).copied().ok_or_else(|| broken(format!("boundary condition `{name}`: set `{s}`")))?, "SET"),
        Target::ReferencePoint(p) => {
            (r.reference_point_nodes.get(p).copied().ok_or_else(|| broken(format!("boundary condition `{name}`: reference point `{p}`")))?, "NODE")
        }
    };
    // A symmetry plane not normal to an axis: held in a coordinate system whose x axis is its
    // normal and whose xy plane holds its first in-plane axis (translation x, rotations y and z).
    if let BoundaryCondition::Symmetry { plane, .. } = bc
        && let (Some(n), Some([t1, _])) = (plane.inclined(), plane.in_plane())
    {
        *constraint_id = constraint_id.saturating_add(1);
        let cid = *constraint_id;
        d.line("*DEFINE_COORDINATE_SYSTEM");
        d.header(&["cid", "xo", "yo", "zo", "xl", "yl", "zl", "cidl"]);
        d.card([Field::from(cid), Field::from(0.0), Field::from(0.0), Field::from(0.0), Field::from(n[0]), Field::from(n[1]), Field::from(n[2])]);
        d.header(&["xp", "yp", "zp"]);
        d.card([Field::from(t1[0]), Field::from(t1[1]), Field::from(t1[2])]);
        d.line(&format!("*BOUNDARY_SPC_{suffix}_ID"));
        id_card(d, constraint_id, name);
        d.header(&["nsid", "cid", "dofx", "dofy", "dofz", "dofrx", "dofry", "dofrz"]);
        d.card([Field::from(entity), Field::from(cid)].into_iter().chain([1, 0, 0, 0, 1, 1].map(|v: i32| Field::from(v))));
        return Ok(());
    }
    let dofs = bc.dofs();
    // Held degrees of freedom (in a frequency step every constrained one is held).
    let held: Vec<u8> = dofs.list().into_iter().filter(|dof| !is_static || bc.value(*dof).is_none_or(|v| v == 0.0)).collect();
    if !held.is_empty() {
        d.line(&format!("*BOUNDARY_SPC_{suffix}_ID"));
        id_card(d, constraint_id, name);
        d.header(&["nsid", "cid", "dofx", "dofy", "dofz", "dofrx", "dofry", "dofrz"]);
        d.card([Field::from(entity), Field::Blank].into_iter().chain((1..=6u8).map(|dof| Field::from(i32::from(held.contains(&dof))))));
    }
    if is_static {
        for dof in dofs.list() {
            let Some(v) = bc.value(dof).filter(|v| *v != 0.0) else { continue };
            d.line(&format!("*BOUNDARY_PRESCRIBED_MOTION_{suffix}_ID"));
            id_card(d, constraint_id, name);
            d.header(&["nsid", "dof", "vad", "lcid", "sf"]);
            d.card([Field::from(entity), Field::from(dof_code(dof)), Field::from(2), Field::from(lcid), Field::from(v)]);
        }
    }
    Ok(())
}

/// The ID card of an `_ID` keyword: the next id and a heading.
fn id_card(d: &mut Deck, id: &mut u32, heading: &str) {
    *id = id.saturating_add(1);
    d.line(&format!("{id:>WIDTH$}{heading}"));
}

/// A force or moment with one value on every node of its set ([`Plan::uniform`]), on the node set.
fn on_node_set(d: &mut Deck, plan: &Plan, name: &str, moment: bool, lcid: u32) {
    let Some(&(v, sid)) = plan.uniform.get(name) else { return };
    d.comment(&format!("load `{name}`"));
    d.line("*LOAD_NODE_SET");
    d.header(&["nsid", "dof", "lcid", "sf"]);
    let first = if moment { 4 } else { 1 };
    for (dof, c) in (first..).zip(v) {
        if c != 0.0 {
            d.card([Field::from(sid), Field::from(dof_code(dof)), Field::from(lcid), Field::from(c)]);
        }
    }
}
