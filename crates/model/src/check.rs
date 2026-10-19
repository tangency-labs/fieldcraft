//! The model check and the resolved view.
//!
//! [`Model::check`] lists everything that keeps the model from being analysed (errors) or looks
//! wrong (warnings), each naming the item it is about. [`Model::resolve`] gives the view the deck
//! writers and the solver read: every set resolved against the current mesh, the section of every
//! element and the node label of every reference point; it fails when the check finds an error.

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};

use crate::loads::Distribution;
use crate::{BoundaryCondition, Elastic, ElementKind, Error, Extent, Load, Model, Procedure, Result, SectionKind, Target};

/// How serious a [`Diagnostic`] is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Severity {
    /// The model cannot be analysed until it is fixed.
    Error,
    /// The model can be analysed but something looks wrong.
    Warning,
}

/// One finding of the check.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct Diagnostic {
    pub severity: Severity,
    /// The item it is about, such as ``section `web` `` or `model`.
    pub subject: String,
    pub message: String,
}

impl std::fmt::Display for Diagnostic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self.severity {
            Severity::Error => "error",
            Severity::Warning => "warning",
        };
        write!(f, "{s}: {}: {}", self.subject, self.message)
    }
}

/// A set resolved against the mesh. Every list is sorted and has no repeats.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct ResolvedSet {
    pub nodes: Vec<u32>,
    pub elements: Vec<u32>,
    pub faces: Vec<crate::ElementFace>,
    /// Element edges (every element kind; 2D models load the edges of their 2D elements).
    #[serde(default)]
    pub edges: Vec<crate::ElementEdge>,
}

impl ResolvedSet {
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty() && self.elements.is_empty() && self.faces.is_empty() && self.edges.is_empty()
    }
}

/// The model resolved against its mesh: what writers and the solver read besides the
/// definitions themselves. When resolution succeeds, every element is in exactly one section's
/// set (`sets[&section.set].elements`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Resolved {
    /// Every set, by name.
    pub sets: BTreeMap<String, ResolvedSet>,
    /// The node label of every reference point ([`Model::reference_point_nodes`]).
    pub reference_point_nodes: BTreeMap<String, u32>,
    /// Forces, moments and tractions as nodal values, sorted by node: what writers and the solver
    /// apply, so every one of them distributes a load the same way. A force or traction gives
    /// force components, a moment moment components.
    pub nodal_loads: BTreeMap<String, Vec<(u32, [f64; 3])>>,
    /// Every point mass at every node of its target, by definition name, then node.
    pub masses: Vec<crate::ResolvedMass>,
    /// Every spring and dashpot at every node of its target (or between its two ends), by
    /// definition name, then node.
    pub springs: Vec<crate::ResolvedSpring>,
    /// The check's warnings (a resolved model has no errors).
    pub warnings: Vec<Diagnostic>,
    /// How each element of a 2D section (plane stress, plane strain, axisymmetric) extends out
    /// of its plane; empty for a model without 2D sections. A model is 2D or not as a whole.
    pub planar: BTreeMap<u32, Extent>,
}

impl Resolved {
    /// Whether the model is a 2D model (its sections are plane stress, plane strain or
    /// axisymmetric).
    pub fn is_planar(&self) -> bool {
        !self.planar.is_empty()
    }

    /// How a 2D model extends out of its plane: revolved when any of its sections is
    /// axisymmetric (the check refuses a mix), else the first element's thickness; `None` for a
    /// model without 2D sections.
    pub fn extent(&self) -> Option<Extent> {
        if self.planar.values().any(|e| *e == Extent::Revolved) { Some(Extent::Revolved) } else { self.planar.values().next().copied() }
    }

    /// The edges of a set that belong to 2D elements: where a 2D model's surface loads act.
    pub fn planar_edges(&self, set: &str) -> Vec<crate::ElementEdge> {
        self.sets.get(set).map(|s| s.edges.iter().filter(|e| self.planar.contains_key(&e.element)).copied().collect()).unwrap_or_default()
    }

    /// The node labels point masses, springs and dashpots act on (repeats possible).
    pub fn discrete_labels(&self) -> impl Iterator<Item = u32> + '_ {
        self.masses.iter().map(|m| m.node).chain(self.springs.iter().flat_map(|s| std::iter::once(s.first).chain(s.second)))
    }

    /// The nodes a boundary condition, force or moment acts on.
    pub fn target_nodes(&self, target: &Target) -> &[u32] {
        match target {
            Target::Set(s) => self.sets.get(s).map(|r| r.nodes.as_slice()).unwrap_or_default(),
            Target::ReferencePoint(r) => self.reference_point_nodes.get(r).map(std::slice::from_ref).unwrap_or_default(),
        }
    }
}

/// Labels listed in a message before it says how many more there are.
const LISTED: usize = 5;

fn first_labels(labels: &[u32]) -> String {
    let shown: Vec<String> = labels.iter().take(LISTED).map(u32::to_string).collect();
    match labels.len().saturating_sub(LISTED) {
        0 => shown.join(", "),
        more => format!("{}, and {more} more", shown.join(", ")),
    }
}

#[derive(Default)]
struct Diagnostics(Vec<Diagnostic>);

impl Diagnostics {
    fn error(&mut self, subject: impl Into<String>, message: impl Into<String>) {
        self.0.push(Diagnostic { severity: Severity::Error, subject: subject.into(), message: message.into() });
    }

    fn warning(&mut self, subject: impl Into<String>, message: impl Into<String>) {
        self.0.push(Diagnostic { severity: Severity::Warning, subject: subject.into(), message: message.into() });
    }

    /// A warning when a magnitude check found one.
    fn warning_if(&mut self, subject: &str, message: Option<String>) {
        if let Some(m) = message {
            self.warning(subject, m);
        }
    }
}

impl Model {
    /// Everything that keeps the model from being analysed or looks wrong, errors first, each
    /// naming its item.
    pub fn check(&self) -> Vec<Diagnostic> {
        self.analyze().1
    }

    /// Resolve the model for writing or solving. Fails with every error the check finds.
    pub fn resolve(&self) -> Result<Resolved> {
        let (mut resolved, diagnostics) = self.analyze();
        let errors: Vec<String> = diagnostics.iter().filter(|d| d.severity == Severity::Error).map(ToString::to_string).collect();
        if !errors.is_empty() {
            return Err(Error::Invalid(format!("the model is not ready for analysis: {}", errors.join("; "))));
        }
        resolved.warnings = diagnostics;
        Ok(resolved)
    }

    /// The mesh nodes with rotational degrees of freedom: those of shells (2D elements that are
    /// not `planar`: plane stress, plane strain or axisymmetric) and of beams. Reference points
    /// rotate too.
    pub fn rotational_nodes(&self, planar: impl Fn(u32) -> bool) -> BTreeSet<u32> {
        self.elements
            .iter()
            .filter(|(l, e)| (e.kind.dimension() == 2 && !planar(**l)) || e.kind == ElementKind::Beam2)
            .flat_map(|(_, e)| e.nodes.iter().copied())
            .collect()
    }

    /// The mesh nodes with rotational degrees of freedom under the model's own sections: those of
    /// shells and beams, not of plane stress, plane strain or axisymmetric elements.
    pub fn rotating_nodes(&self) -> BTreeSet<u32> {
        let planar: BTreeSet<u32> = self
            .sections
            .values()
            .filter(|s| s.kind.extent().is_some())
            .filter_map(|s| self.resolve_set_partial(&s.set))
            .flat_map(|r| r.elements)
            .collect();
        self.rotational_nodes(|l| planar.contains(&l))
    }

    fn analyze(&self) -> (Resolved, Vec<Diagnostic>) {
        let mut d = Diagnostics::default();
        let mut resolved = Resolved::default();
        let shown = self.display_units.clone().unwrap_or_default();
        let modulus =
            |what: &str, v: f64| crate::units::check_magnitude(what, crate::units::Quantity::Pressure, v, crate::units::MODULUS_RANGE, &shown);
        // Whether a step needs elasticity (every step but pure heat transfer).
        let structural = self.steps.is_empty() || self.steps.iter().any(|s| s.procedure.is_structural());

        // A model of fluid steps alone needs no mesh: its domains carry their own grid; one of
        // point masses, springs and dashpots alone is its own mesh.
        let fluid_only = !self.steps.is_empty() && self.steps.iter().all(|s| s.procedure.fluid_domain().is_some());
        if self.elements.is_empty() && !fluid_only && self.masses.is_empty() && self.springs.is_empty() {
            d.error("model", "the mesh has no elements");
        }
        // Sets.
        for (name, set) in &self.sets {
            let (r, missing) = self.resolve_definition(set);
            for m in missing {
                d.error(format!("set `{name}`"), m);
            }
            if r.is_empty() {
                if self.set_users(name).is_empty() {
                    d.warning(format!("set `{name}`"), "resolves to nothing on the current mesh");
                } else {
                    d.error(format!("set `{name}`"), "resolves to nothing on the current mesh");
                }
            }
            resolved.sets.insert(name.clone(), r);
        }
        // The elements of 2D sections.
        for sec in self.sections.values() {
            if let (Some(extent), Some(r)) = (sec.kind.extent(), resolved.sets.get(&sec.set)) {
                for &e in &r.elements {
                    resolved.planar.insert(e, extent);
                }
            }
        }
        let empty = crate::ResolvedSet::default();
        let set_of = |name: &str| resolved.sets.get(name).unwrap_or(&empty);

        // Materials: magnitudes for the unit system.
        for (name, m) in &self.materials {
            let subject = format!("material `{name}`");
            match m.elastic {
                Some(Elastic::Isotropic { young, .. }) => d.warning_if(&subject, modulus("Young's modulus", young)),
                Some(Elastic::EngineeringConstants { young, shear, .. }) => {
                    for (i, e) in young.iter().enumerate() {
                        d.warning_if(&subject, modulus(&format!("E{}", i + 1), *e));
                    }
                    for (n, g) in crate::materials::SHEAR_MODULI.iter().zip(shear) {
                        d.warning_if(&subject, modulus(n, g));
                    }
                }
                None => {}
            }
            if let Some(rho) = m.density {
                let check = crate::units::check_magnitude("density", crate::units::Quantity::Density, rho, crate::units::DENSITY_RANGE, &shown);
                d.warning_if(&subject, check);
            }
        }

        // Sections: element kinds, formulation, one section per element.
        let mut owner: BTreeMap<u32, &str> = BTreeMap::new();
        for (name, sec) in &self.sections {
            let subject = format!("section `{name}`");
            let r = set_of(&sec.set);
            if r.elements.is_empty() {
                d.error(&subject, format!("set `{}` has no elements", sec.set));
            }
            let fits = |k: ElementKind| match sec.kind {
                SectionKind::Solid {} => k.dimension() == 3,
                SectionKind::Shell { .. } | SectionKind::PlaneStress { .. } | SectionKind::PlaneStrain { .. } | SectionKind::Axisymmetric {} => {
                    k.dimension() == 2
                }
                SectionKind::Beam { .. } => k == ElementKind::Beam2,
                SectionKind::Rod { .. } => matches!(k, ElementKind::Rod2 | ElementKind::Rod3),
            };
            // The linear element with an incompatible-modes formulation for this kind of section.
            let incompatible_kind = if sec.kind.is_planar() { ElementKind::Quad4 } else { ElementKind::Hex8 };
            let mut wrong = Vec::new();
            let mut not_hex8 = Vec::new();
            let mut along = Vec::new();
            let mut taken: BTreeMap<&str, Vec<u32>> = BTreeMap::new();
            for &e in &r.elements {
                let Some(el) = self.elements.get(&e) else { continue };
                if !fits(el.kind) {
                    wrong.push(e);
                }
                if sec.formulation.incompatible_modes && el.kind != incompatible_kind {
                    not_hex8.push(e);
                }
                if let SectionKind::Beam { orientation, .. } = sec.kind
                    && self.parallel_to_axis(&el.nodes, orientation)
                {
                    along.push(e);
                }
                match owner.get(&e) {
                    Some(other) => taken.entry(other).or_default().push(e),
                    None => {
                        owner.insert(e, name);
                    }
                }
            }
            if !wrong.is_empty() {
                d.error(&subject, format!("a {} section cannot hold elements {}", sec.kind.name(), first_labels(&wrong)));
            }
            if !not_hex8.is_empty() {
                d.error(
                    &subject,
                    format!(
                        "incompatible modes apply to {} elements only, not elements {} (in a {} section)",
                        incompatible_kind.name(),
                        first_labels(&not_hex8),
                        sec.kind.name()
                    ),
                );
            }

            if !along.is_empty() {
                d.error(&subject, format!("the orientation vector is parallel to elements {}", first_labels(&along)));
            }
            for (other, labels) in taken {
                d.error(&subject, format!("elements {} are also in section `{other}`", first_labels(&labels)));
            }
            if structural && self.materials.get(&sec.material).is_some_and(|m| m.elastic.is_none()) {
                d.error(&subject, format!("material `{}` has no elastic constants", sec.material));
            }
        }
        let unsectioned: Vec<u32> = self.elements.keys().filter(|e| !owner.contains_key(e)).copied().collect();
        if !unsectioned.is_empty() {
            d.error("model", format!("elements {} have no section", first_labels(&unsectioned)));
        }
        self.check_planar(&resolved, &mut d);

        // Mass: frequency steps and gravity need a density in every material in use.
        let needs_mass = self
            .steps
            .iter()
            .find_map(|s| {
                let why = match s.procedure {
                    Procedure::Frequency { .. } => "extracts frequencies",
                    Procedure::ExplicitDynamic(_) => "is dynamic",
                    Procedure::HeatTransfer { transient: Some(_) } => "is transient heat transfer",
                    Procedure::Static {}
                    | Procedure::Buckle { .. }
                    | Procedure::NonlinearStatic(_)
                    | Procedure::Fluid { .. }
                    | Procedure::HeatTransfer { transient: None }
                    | Procedure::CoupledTemperatureDisplacement {} => {
                        return None;
                    }
                };
                Some(format!("step `{}` {why}", s.name))
            })
            .or_else(|| {
                self.loads.iter().find_map(|(n, l)| match l {
                    Load::Gravity { .. } => Some(format!("load `{n}` is gravity")),
                    Load::Centrifugal { .. } => Some(format!("load `{n}` is centrifugal")),
                    _ => None,
                })
            });
        if let Some(why) = needs_mass {
            let used: BTreeSet<&str> = self.sections.values().map(|s| s.material.as_str()).collect();
            for m in used {
                if self.materials.get(m).is_some_and(|m| m.density.is_none()) {
                    d.error(format!("material `{m}`"), format!("has no density, and {why}"));
                }
            }
        }

        // Heat transfer needs conductivities (and heat capacities when transient); thermal strain
        // needs expansion coefficients.
        let used: BTreeSet<&str> = self.sections.values().map(|s| s.material.as_str()).collect();
        let thermal = |m: &str| self.materials.get(m).and_then(|m| m.thermal);
        for s in self.steps.iter().filter(|s| s.procedure.is_thermal()) {
            for &m in &used {
                let t = thermal(m);
                if t.and_then(|t| t.conductivity).is_none() {
                    d.error(format!("material `{m}`"), format!("has no conductivity, and step `{}` conducts heat", s.name));
                }
                if matches!(s.procedure, Procedure::HeatTransfer { transient: Some(_) }) && t.and_then(|t| t.specific_heat).is_none() {
                    d.error(format!("material `{m}`"), format!("has no specific heat, and step `{}` is transient heat transfer", s.name));
                }
            }
        }
        let strains_thermally = self.steps.iter().any(|s| {
            matches!(s.procedure, Procedure::CoupledTemperatureDisplacement {})
                || s.loads.iter().any(|l| matches!(self.loads.get(l), Some(Load::Temperature { .. })))
        });
        if strains_thermally {
            for &m in &used {
                if thermal(m).and_then(|t| t.expansion).is_none() {
                    d.warning(format!("material `{m}`"), "has no expansion coefficient, so temperatures do not strain it");
                }
            }
        }

        // Reference points and couplings.
        resolved.reference_point_nodes = self.reference_point_nodes();
        for (name, label) in &resolved.reference_point_nodes {
            if self.nodes.contains_key(label) {
                d.error(format!("reference point `{name}`"), format!("its node label {label} is also a mesh node"));
            }
        }
        for (name, c) in &self.couplings {
            if set_of(&c.set).nodes.is_empty() {
                d.error(format!("coupling `{name}`"), format!("set `{}` has no nodes", c.set));
            }
        }

        // Contacts: their surfaces are element faces of solids, or edges of 2D elements.
        for (name, c) in &self.contacts {
            for set in c.surfaces.sets() {
                let faces = &set_of(set).faces;
                let solid = |f: &crate::ElementFace| self.elements.get(&f.element).is_some_and(|e| e.kind.dimension() == 3);
                if faces.is_empty() && resolved.planar_edges(set).is_empty() {
                    d.error(
                        format!("contact `{name}`"),
                        format!("set `{set}` has no element faces (a contact surface is a set of faces, or of 2D element edges)"),
                    );
                } else if !faces.iter().all(solid) {
                    d.error(
                        format!("contact `{name}`"),
                        format!("set `{set}` has faces of elements that are not solids (contact acts on the faces of solid elements so far)"),
                    );
                }
            }
        }

        // Ties: their surfaces are element faces of solids too.
        for (name, t) in &self.ties {
            for set in [&t.main, &t.secondary] {
                let faces = &set_of(set).faces;
                let solid = |f: &crate::ElementFace| self.elements.get(&f.element).is_some_and(|e| e.kind.dimension() == 3);
                if faces.is_empty() {
                    d.error(format!("tie `{name}`"), format!("set `{set}` has no element faces (a tie surface is a set of faces)"));
                } else if !faces.iter().all(solid) {
                    d.error(
                        format!("tie `{name}`"),
                        format!("set `{set}` has faces of elements that are not solids (ties bond the faces of solid elements so far)"),
                    );
                }
            }
        }

        // Nodes with rotational degrees of freedom: those of shells and beams, and reference points
        // (worked out only when a boundary condition or load needs it).
        let rotational = std::cell::OnceCell::new();
        let has_rotations = |t: &Target| match t {
            Target::ReferencePoint(_) => true,
            Target::Set(s) => {
                let rotational = rotational.get_or_init(|| self.rotational_nodes(|l| resolved.planar.contains_key(&l)));
                set_of(s).nodes.iter().any(|n| rotational.contains(n))
            }
        };

        // Boundary conditions.
        for (name, bc) in &self.boundary_conditions {
            let subject = format!("boundary condition `{name}`");
            if let Target::Set(s) = bc.target()
                && set_of(s).nodes.is_empty()
            {
                d.error(&subject, format!("set `{s}` has no nodes"));
            }
            if let BoundaryCondition::Displacement { values, target } = bc
                && values.get(3..).is_some_and(|r| r.iter().any(Option::is_some))
                && !has_rotations(target)
            {
                d.warning(&subject, format!("{target} has only solid-element nodes, which have no rotations to prescribe"));
            }
        }

        // Point masses, springs and dashpots.
        let (masses, springs) = self.resolve_discrete(&resolved, &mut d, &has_rotations);

        // Loads.
        let mut nodal_loads = BTreeMap::new();
        for (name, load) in &self.loads {
            let subject = format!("load `{name}`");
            match self.nodal_values(load, &resolved) {
                Ok(Some(values)) => {
                    nodal_loads.insert(name.clone(), values);
                }
                Ok(None) => {}
                Err(m) => d.error(&subject, m),
            }
            match load {
                Load::Moment { target, .. } if !has_rotations(target) => {
                    d.error(&subject, format!("{target} has no nodes with rotations; apply the moment to a reference point coupled to it"))
                }
                Load::Pressure { set, .. } if set_of(set).faces.is_empty() && resolved.planar_edges(set).is_empty() => {
                    d.error(&subject, format!("set `{set}` has no element faces (or 2D element edges) to press on"))
                }
                Load::BodyForce { set, .. } if set_of(set).elements.is_empty() => d.error(&subject, format!("set `{set}` has no elements")),
                Load::Gravity { acceleration, set } => {
                    let g = acceleration.iter().map(|c| c * c).sum::<f64>().sqrt();
                    let check =
                        crate::units::check_magnitude("gravity", crate::units::Quantity::Acceleration, g, crate::units::GRAVITY_RANGE, &shown);
                    d.warning_if(&subject, check);
                    if let Some(s) = set
                        && set_of(s).elements.is_empty()
                    {
                        d.error(&subject, format!("set `{s}` has no elements"));
                    }
                }
                Load::Centrifugal { set: Some(s), .. } if set_of(s).elements.is_empty() => d.error(&subject, format!("set `{s}` has no elements")),
                _ => {}
            }
            if resolved.is_planar() {
                Self::check_planar_load(load, resolved.extent(), &subject, &mut d);
            }
        }

        // Steps and the activity table.
        if self.steps.is_empty() {
            d.error("model", "there is no analysis step");
        }
        for s in &self.steps {
            let subject = format!("step `{}`", s.name);
            let prescribes_motion = s
                .boundary_conditions
                .iter()
                .filter_map(|b| self.boundary_conditions.get(b))
                .any(|b| (1..=6).any(|dof| b.value(dof).is_some_and(|v| v != 0.0)));
            for set in s.procedure.sets() {
                if set_of(set).nodes.is_empty() {
                    d.error(&subject, format!("set `{set}` has no nodes"));
                }
            }
            match &s.procedure {
                Procedure::Static {} | Procedure::NonlinearStatic(_) | Procedure::Buckle { .. } => {
                    if s.loads.is_empty() && !prescribes_motion {
                        d.warning(&subject, "nothing loads the model in this step");
                    }
                }
                Procedure::Frequency { .. } => {
                    if !s.loads.is_empty() {
                        d.warning(&subject, "loads have no effect on natural frequencies");
                    }
                }
                Procedure::ExplicitDynamic(e) => {
                    let moves = e.initial_velocities.iter().any(|v| v.velocity.iter().any(|c| *c != 0.0));
                    if s.loads.is_empty() && !prescribes_motion && !moves {
                        d.warning(&subject, "nothing loads or moves the model in this step");
                    }
                }
                // The fluid domain holds its own conditions; nothing structural applies.
                Procedure::Fluid { .. } => continue,
                Procedure::HeatTransfer { .. } | Procedure::CoupledTemperatureDisplacement {} => {
                    let bcs = s.boundary_conditions.iter().filter_map(|b| self.boundary_conditions.get(b));
                    let loads: Vec<&Load> = s.loads.iter().filter_map(|l| self.loads.get(l)).collect();
                    let held = bcs.clone().any(BoundaryCondition::is_thermal) || loads.iter().any(|l| matches!(l, Load::Convection { .. }));
                    if !held && matches!(s.procedure, Procedure::HeatTransfer { transient: None } | Procedure::CoupledTemperatureDisplacement {}) {
                        d.warning(&subject, "no temperature boundary condition or convection is active, so the temperature level is undetermined");
                    }
                }
            }
            for l in s.loads.iter().filter(|l| self.loads.get(*l).is_some_and(|l| l.is_thermal() && !s.procedure.is_thermal())) {
                d.warning(&subject, format!("heat load `{l}` has no effect in a {} step", s.procedure.name()));
            }
            // A free body is a natural case of transient dynamics (flight, impact); heat transfer
            // has no rigid motion.
            if s.procedure.is_structural()
                && s.boundary_conditions.iter().filter_map(|b| self.boundary_conditions.get(b)).all(BoundaryCondition::is_thermal)
                && !s.procedure.is_dynamic()
            {
                d.warning(&subject, "no boundary condition is active, so the model is free to move as a rigid body");
            }
        }
        let active_bcs: BTreeSet<&String> = self.steps.iter().flat_map(|s| &s.boundary_conditions).collect();
        let active_loads: BTreeSet<&String> = self.steps.iter().flat_map(|s| &s.loads).collect();
        for n in self.boundary_conditions.keys().filter(|n| !active_bcs.contains(n)) {
            d.warning(format!("boundary condition `{n}`"), "is not active in any step");
        }
        for n in self.loads.keys().filter(|n| !active_loads.contains(n)) {
            d.warning(format!("load `{n}`"), "is not active in any step");
        }

        resolved.nodal_loads = nodal_loads;
        resolved.masses = masses;
        resolved.springs = springs;
        let mut out = d.0;
        out.sort_by_key(|a| a.severity);
        (resolved, out)
    }

    /// Point masses, springs and dashpots expanded to their nodes, reporting what cannot be
    /// expanded.
    fn resolve_discrete(
        &self,
        resolved: &Resolved,
        d: &mut Diagnostics,
        has_rotations: &dyn Fn(&Target) -> bool,
    ) -> (Vec<crate::ResolvedMass>, Vec<crate::ResolvedSpring>) {
        use crate::{Direction, ResolvedDirection, ResolvedMass, ResolvedSpring};
        let mut masses = Vec::new();
        let mut springs = Vec::new();
        let mut count = 0usize;
        let mut over = false;
        let mut take = |n: usize| {
            count = count.saturating_add(n);
            over |= count > crate::discrete::MAX_DISCRETE;
            !over
        };
        for (name, m) in &self.masses {
            let subject = format!("point mass `{name}`");
            let nodes = resolved.target_nodes(&m.target);
            if nodes.is_empty() {
                d.error(&subject, format!("{} has no nodes", m.target));
                continue;
            }
            if m.inertia.is_some_and(|i| i.iter().any(|c| *c != 0.0)) && !has_rotations(&m.target) {
                d.error(&subject, format!("{} has no nodes with rotations for its inertia; put it on a reference point", m.target));
            }
            if !take(nodes.len()) {
                break;
            }
            let mass = m.mass.map_or([0.0; 3], crate::Mass::components);
            let inertia = m.inertia.unwrap_or_default();
            masses.extend(nodes.iter().map(|n| ResolvedMass { name: name.clone(), node: *n, mass, inertia }));
        }
        for (name, s) in &self.springs {
            let subject = format!("{} `{name}`", s.kind());
            let first = resolved.target_nodes(s.first());
            let rotational = s.direction().is_rotational();
            if first.is_empty() {
                d.error(&subject, format!("{} has no nodes", s.first()));
                continue;
            }
            if rotational && s.targets().any(|t| !has_rotations(t)) {
                d.error(&subject, "acts on a rotation, and its nodes have none; connect reference points (or shell and beam nodes)");
            }
            let unit = |v: [f64; 3]| {
                let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                (len > 0.0 && len.is_finite()).then(|| v.map(|c| c / len))
            };
            let pairs: Vec<(u32, Option<u32>)> = match s.second() {
                None => first.iter().map(|n| (*n, None)).collect(),
                Some(t) => {
                    let second = resolved.target_nodes(t);
                    if first.len() != 1 || second.len() != 1 {
                        d.error(
                            &subject,
                            format!(
                                "a two-point {} connects one node to one node; {} has {} and {t} has {} (use a vertex set or a reference point at each end)",
                                s.kind(),
                                s.first(),
                                first.len(),
                                second.len()
                            ),
                        );
                        continue;
                    }
                    match (first.first(), second.first()) {
                        (Some(a), Some(b)) if a == b => {
                            d.error(&subject, format!("both ends are node {a}"));
                            continue;
                        }
                        (Some(a), Some(b)) => vec![(*a, Some(*b))],
                        _ => continue,
                    }
                }
            };
            if !take(pairs.len()) {
                break;
            }
            for (a, b) in pairs {
                let direction = match s.direction() {
                    Direction::Dof(k) => ResolvedDirection::Dof(k),
                    Direction::Vector(v) => match unit(v) {
                        Some(n) => ResolvedDirection::Vector(n),
                        None => {
                            d.error(&subject, "its direction vector has no length");
                            break;
                        }
                    },
                    Direction::Line => {
                        let ends = b.and_then(|b| Some((self.label_xyz(a)?, self.label_xyz(b)?)));
                        match ends.and_then(|(x, y)| unit([y[0] - x[0], y[1] - x[1], y[2] - x[2]])) {
                            Some(n) => ResolvedDirection::Line(n),
                            None => {
                                d.error(&subject, "its two ends are at the same point, so it has no line to act along");
                                break;
                            }
                        }
                    }
                };
                springs.push(ResolvedSpring { name: name.clone(), dashpot: s.is_dashpot(), first: a, second: b, direction, value: s.value() });
            }
        }
        if over {
            d.error("model", format!("point masses, springs and dashpots resolve to more than {} copies", crate::discrete::MAX_DISCRETE));
        }
        (masses, springs)
    }

    /// A force, moment or traction as nodal values (`None` for other loads), or why it cannot be
    /// distributed.
    fn nodal_values(&self, load: &Load, r: &Resolved) -> std::result::Result<Option<Vec<(u32, [f64; 3])>>, String> {
        let scale = |v: [f64; 3], k: f64| v.map(|c| c * k);
        let (target, vector, distribution) = match load {
            Load::Force { target, vector, distribution } | Load::Moment { target, vector, distribution } => (target, *vector, *distribution),
            Load::Traction { set, vector } => {
                let faces = r.sets.get(set).map(|s| s.faces.as_slice()).unwrap_or_default();
                let edges = r.planar_edges(set);
                if faces.is_empty() && edges.is_empty() {
                    return Err(format!("set `{set}` has no element faces (or 2D element edges) to load"));
                }
                let w = if faces.is_empty() { self.edge_weights(&edges, &r.planar) } else { self.surface_weights(faces) };
                return Ok(Some(w.map_err(|e| e.to_string())?.into_iter().map(|(n, a)| (n, scale(*vector, a))).collect()));
            }
            _ => return Ok(None),
        };
        let empty = ResolvedSet::default();
        let set = match target {
            Target::ReferencePoint(rp) => {
                return r
                    .reference_point_nodes
                    .get(rp)
                    .map(|n| Some(vec![(*n, vector)]))
                    .ok_or_else(|| format!("reference point `{rp}` has no node"));
            }
            Target::Set(s) => r.sets.get(s).unwrap_or(&empty),
        };
        let name = target.set().unwrap_or_default();
        if set.nodes.is_empty() {
            return Err(format!("set `{name}` has no nodes"));
        }
        let edges = if set.faces.is_empty() { r.planar_edges(name) } else { Vec::new() };
        Ok(Some(match distribution {
            Distribution::PerNode => set.nodes.iter().map(|n| (*n, vector)).collect(),
            Distribution::Total if set.faces.is_empty() && edges.is_empty() => {
                let share = 1.0 / set.nodes.len() as f64;
                set.nodes.iter().map(|n| (*n, scale(vector, share))).collect()
            }
            Distribution::Total => {
                let w = if edges.is_empty() { self.surface_weights(&set.faces) } else { self.edge_weights(&edges, &r.planar) };
                let w = w.map_err(|e| e.to_string())?;
                let area: f64 = w.values().sum();
                if !(area > 0.0 && area.is_finite()) {
                    return Err(format!("the faces of set `{name}` have no area to share the load by"));
                }
                w.into_iter().map(|(n, a)| (n, scale(vector, a / area))).collect()
            }
        }))
    }

    /// The rules of a 2D model: one kind of 2D model (plane or axisymmetric) and nothing else but
    /// trusses (plane models only), in the XY plane (x ≥ 0 when axisymmetric), elements numbered
    /// counter-clockwise.
    fn check_planar(&self, r: &Resolved, d: &mut Diagnostics) {
        if !r.is_planar() {
            return;
        }
        let axisymmetric = r.extent() == Some(Extent::Revolved);
        // Trusses join a plane model with its two in-plane unknowns per node.
        let truss = |s: &crate::Section| matches!(s.kind, SectionKind::Rod { .. }) && !axisymmetric;
        let others: Vec<&str> = self.sections.iter().filter(|(_, s)| !s.kind.is_planar() && !truss(s)).map(|(n, _)| n.as_str()).collect();
        if !others.is_empty() {
            d.error(
                "model",
                format!(
                    "2D sections (plane stress, plane strain, axisymmetric) and other sections ({}) cannot be mixed in one model yet",
                    others.join(", ")
                ),
            );
        }
        if axisymmetric && r.planar.values().any(|e| *e != Extent::Revolved) {
            d.error("model", "axisymmetric and plane (stress or strain) sections cannot be mixed in one model");
        }
        let (mut off_plane, mut negative_radius, mut clockwise) = (Vec::new(), Vec::new(), Vec::new());
        for &label in r.planar.keys() {
            let Some(e) = self.elements.get(&label) else { continue };
            let xyz: Vec<[f64; 3]> = e.nodes.iter().filter_map(|n| self.nodes.get(n)).map(|n| n.xyz).collect();
            if xyz.iter().any(|p| p[2] != 0.0) {
                off_plane.push(label);
            }
            if axisymmetric && xyz.iter().any(|p| p[0] < 0.0) {
                negative_radius.push(label);
            }
            // Twice the signed area of the corner polygon (shoelace).
            let corners = xyz.get(..e.kind.corner_count()).unwrap_or_default();
            let area2: f64 = corners.iter().zip(corners.iter().cycle().skip(1)).map(|(a, b)| a[0] * b[1] - b[0] * a[1]).sum();
            if !(area2 > 0.0) {
                clockwise.push(label);
            }
        }
        for sec in self.sections.values().filter(|s| truss(s)) {
            for label in r.sets.get(&sec.set).map(|s| s.elements.as_slice()).unwrap_or_default() {
                let Some(e) = self.elements.get(label) else { continue };
                if e.nodes.iter().filter_map(|n| self.nodes.get(n)).any(|n| n.xyz[2] != 0.0) {
                    off_plane.push(*label);
                }
            }
        }
        if !off_plane.is_empty() {
            d.error("model", format!("2D elements {} have nodes off the XY plane (z ≠ 0); 2D models lie in the XY plane", first_labels(&off_plane)));
        }
        if !negative_radius.is_empty() {
            d.error(
                "model",
                format!(
                    "axisymmetric elements {} have nodes at negative x; x is the radius, so the model lies in x ≥ 0",
                    first_labels(&negative_radius)
                ),
            );
        }
        if !clockwise.is_empty() {
            d.error(
                "model",
                format!(
                    "2D elements {} are clockwise or have no area in the XY plane (2D elements number their nodes counter-clockwise)",
                    first_labels(&clockwise)
                ),
            );
        }
    }

    /// A load in a 2D model: in the XY plane, and a rotation about the axis the model can carry.
    fn check_planar_load(load: &Load, extent: Option<Extent>, subject: &str, d: &mut Diagnostics) {
        let in_plane = match load {
            Load::Force { vector, .. } | Load::Traction { vector, .. } | Load::BodyForce { vector, .. } => vector[2] == 0.0,
            Load::Gravity { acceleration, .. } => acceleration[2] == 0.0,
            Load::Moment { .. } | Load::Pressure { .. } => true,
            // Scalars.
            Load::SurfaceFlux { .. } | Load::Convection { .. } | Load::HeatSource { .. } | Load::Temperature { .. } => true,
            Load::Centrifugal { point, axis, .. } => {
                let ok = if extent == Some(Extent::Revolved) {
                    axis[0] == 0.0 && axis[2] == 0.0 && point[0] == 0.0 && point[2] == 0.0
                } else {
                    axis[0] == 0.0 && axis[1] == 0.0
                };
                if !ok {
                    d.error(
                        subject,
                        "in a 2D model a rotation is about the Y axis (axisymmetric: through x = 0, z = 0) or about an axis along Z (plane models)",
                    );
                }
                true
            }
        };
        if !in_plane {
            d.error(subject, "a 2D model is loaded in the XY plane; this load has a Z component");
        }
    }

    /// Whether a vector is parallel to the axis of a two-node element (or the element has no
    /// length).
    fn parallel_to_axis(&self, nodes: &[u32], v: [f64; 3]) -> bool {
        let [Some(a), Some(b)] = [nodes.first(), nodes.get(1)].map(|n| n.and_then(|n| self.nodes.get(n))) else { return false };
        let t = [b.xyz[0] - a.xyz[0], b.xyz[1] - a.xyz[1], b.xyz[2] - a.xyz[2]];
        let c = [t[1] * v[2] - t[2] * v[1], t[2] * v[0] - t[0] * v[2], t[0] * v[1] - t[1] * v[0]];
        let norm = |x: [f64; 3]| (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        !(norm(c) > 1e-9 * norm(t) * norm(v))
    }
}
