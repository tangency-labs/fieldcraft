//! The solver-neutral finite element model.
//!
//! A [`Model`] holds what every solver deck describes, in no solver's vocabulary:
//! - the mesh: nodes and elements addressed by the positive integer labels the decks use
//!   ([`mesh`]), and how the mesh lies on the geometry ([`geometry`]: each element, element face,
//!   element edge and node tagged with the persistent name of the body, face, edge or vertex it
//!   was meshed on);
//! - named sets ([`sets`]) defined by geometry names or explicit labels, resolved to nodes,
//!   elements and element faces against the current mesh, so they survive remeshing;
//! - materials ([`materials`]), sections with their element
//!   formulation ([`sections`]), reference points and couplings ([`constraints`]), point masses,
//!   springs and dashpots ([`discrete`]);
//! - boundary conditions ([`bcs`]), loads ([`loads`]) and steps with an explicit table of which
//!   boundary conditions and loads are active in each step ([`steps`]);
//! - the check that names every missing or bad item, and the resolved view the deck writers and
//!   the solver read ([`check`]).
//!
//! Solver formats (Abaqus input, Nastran bulk data, LS-DYNA keyword) translate to and from this
//! model; nothing here knows any of them.
//!
//! Every number is stored in SI (m, kg, s, N, Pa, K, rad; [`units`]), whatever the model is
//! shown in: the display system (`display_units`, or the user's when the model sets none) is
//! applied at the boundaries (commands, forms, replies, files of other formats), never here.
//!
//! Every mutation validates its input, because labels, names and numbers come from files, scripts
//! and the control channel: a bad value is an [`Error`], never a panic. Definitions that refer to
//! other definitions (a section to its set and material, a load to its set) must refer to ones
//! that exist, and a definition still referred to cannot be deleted. What depends on the mesh (a
//! set's geometry names and labels) is checked when the model is resolved, because remeshing
//! replaces the mesh under the definitions.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]

pub mod bcs;
pub mod check;
pub mod constraints;
pub mod contact;
mod convert;
pub mod datums;
pub mod discrete;
pub mod fluid_zones;
pub mod fluids;
pub mod geometry;
pub mod geometry_names;
pub mod loads;
pub mod materials;
mod mesh;
pub mod mesh_plans;
pub mod nodal;
pub mod progress;
pub mod relaxations;
pub mod rename;
pub mod sections;
pub mod seeds;
pub mod sets;
pub mod sizes;
pub mod steps;
pub mod tie;
pub mod units;
pub mod virtual_faces;

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};

pub use bcs::{BoundaryCondition, Dofs, Symmetry, TEMPERATURE_DOF, Target};
pub use check::{Diagnostic, Resolved, ResolvedSet, Severity};
pub use constraints::{Coupling, CouplingKind, ReferencePoint};
pub use contact::{Contact, ContactSurfaces, Friction, InitialPenetration};
pub use datums::{Datum, Partition};
pub use discrete::{Direction, Mass, PointMass, ResolvedDirection, ResolvedMass, ResolvedSpring, Spring};
pub use fluid_zones::{Axis, Fan, FanDirection, PorousZone, Resistance};
pub use fluids::{
    BoxFace, BoxFaces, BoxThermalFaces, ConductingSolid, FaceCondition, FluidBoundary, FluidDomain, FluidThermal, Outflow, Profile, SurfaceBoundary,
    SurfaceCondition, SurfaceThermalBoundary, ThermalCondition, ThermalContact,
};
pub use geometry::{BodyMesh, Classification, Entity, GeometryRef};
pub use geometry_names::NameUse;
pub use loads::{Distribution, Load, centrifugal_acceleration};
pub use materials::{Damage, Directional, Elastic, Erosion, Fluid, Hardening, HardeningPoint, Hyperelastic, Material, Plastic, Thermal};
pub use mesh::{
    Element, ElementEdge, ElementFace, ElementKind, MIN_NORMALIZED_JACOBIAN, Node, admissible_jacobian, check_coordinates, normalized_jacobian,
    tet4_jacobian,
};
pub use mesh_plans::{CutRole, MeshPlan, PlanCut, PlanFeature, PlanMesher, PlanRing, PlanUnmapped, RingKind};
pub use relaxations::{Relaxation, RelaxationCell};
pub use sections::{BeamProfile, Extent, Formulation, Integration, Section, SectionKind};
pub use seeds::{Bias, BiasToward, Divisions, Seed, SeedConstraint};
pub use sets::Set;
pub use sizes::{SizeControl, SizeRegion};
pub use steps::{
    ArcConstraint, ArcLength, Averaging, BulkViscosity, Collision, DisplacementLimit, ExplicitDynamic, FlowMonitor, Hourglass, HourglassKind,
    Imperfection, InitialVelocity, MAX_LOAD_PATH_POINTS, MAX_NONLINEAR_INCREMENTS, MAX_NONLINEAR_ITERATIONS, MassScaling, NonlinearStatic,
    OutputVariable, Procedure, RigidWall, Steady, Step, StepInputs, StressRate, Tolerances, Transient, Turbulence,
};
pub use tie::Tie;
pub use units::DisplaySystem;
pub use virtual_faces::{Defeature, FaceMerge, FaceSplit, Stitch};

/// The largest node or element label. Nastran's GRID and element ids stop at 99 999 999 (eight
/// digits in a small field), the tightest of the three formats.
pub const MAX_LABEL: u32 = 99_999_999;

/// Caps on model size, so hostile input cannot exhaust memory.
pub const MAX_NODES: usize = 50_000_000;
pub const MAX_ELEMENTS: usize = 50_000_000;

/// Coordinates beyond this magnitude are rejected (metres).
pub const MAX_COORD: f64 = 1e9; // unit-factors: a cap on hostile input, in metres (SI)

/// Most definitions of one kind (sets, materials, sections, …) a model holds.
pub const MAX_DEFINITIONS: usize = 100_000;

/// Longest name of a definition. Names are ASCII: a letter, then letters, digits, `_` or `-`, so
/// every deck format can carry them unquoted.
pub const MAX_NAME_LEN: usize = 64;

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum Error {
    #[error("invalid input: {0}")]
    Invalid(String),
    #[error("not found: {0}")]
    NotFound(String),
    #[error("limit reached: {0}")]
    Limit(String),
    /// A definition cannot be deleted while another refers to it.
    #[error("in use: {0}")]
    InUse(String),
    /// The command was cancelled while it ran ([`progress::Progress::cancel`]); it left
    /// nothing behind.
    #[error("cancelled")]
    Cancelled,
}

pub type Result<T> = std::result::Result<T, Error>;

/// The finite element model.
///
/// It deserializes by replaying what it reads through the same checked methods commands use, so
/// a model read from a file holds the same invariants as one built by commands.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "ModelData")]
pub struct Model {
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    nodes: BTreeMap<u32, Node>,
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    elements: BTreeMap<u32, Element>,
    #[serde(skip_serializing_if = "Classification::is_unused")]
    classification: Classification,
    /// The unit system the model is shown and typed in, when it sets its own (a model of format
    /// 2 or older shows the system it declared); `None` follows the user's.
    #[serde(skip_serializing_if = "Option::is_none")]
    display_units: Option<DisplaySystem>,
    sets: BTreeMap<String, Set>,
    materials: BTreeMap<String, Material>,
    sections: BTreeMap<String, Section>,
    reference_points: BTreeMap<String, ReferencePoint>,
    couplings: BTreeMap<String, Coupling>,
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    contacts: BTreeMap<String, Contact>,
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    ties: BTreeMap<String, Tie>,
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    masses: BTreeMap<String, PointMass>,
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    springs: BTreeMap<String, Spring>,
    boundary_conditions: BTreeMap<String, BoundaryCondition>,
    loads: BTreeMap<String, Load>,
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    fluid_domains: BTreeMap<String, FluidDomain>,
    /// In analysis order.
    steps: Vec<Step>,
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    datums: BTreeMap<String, Datum>,
    /// In the order they cut.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    partitions: Vec<Partition>,
    /// Datums kept but left out of the model ([`Model::datum_suppressed`]).
    #[serde(skip_serializing_if = "BTreeSet::is_empty")]
    suppressed_datums: BTreeSet<String>,
    /// Edge seeds, sorted by edge.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    seeds: Vec<Seed>,
    /// Size controls (`mesh.size.set`), by name.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    size_controls: BTreeMap<String, SizeControl>,
    /// The size field's growth rate (`mesh.size.growth`); `None`: the default.
    #[serde(skip_serializing_if = "Option::is_none")]
    size_growth: Option<f64>,
    /// Relaxation fields (`mesh.node.relax`), by name.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    relaxations: BTreeMap<String, Relaxation>,
    /// Mesh plans (`mesh.plan`), by name.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    mesh_plans: BTreeMap<String, MeshPlan>,
    /// Face merges (virtual topology), in the order they apply.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    face_merges: Vec<FaceMerge>,
    /// Stitches of imported surfaces, in the order they apply.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    stitches: Vec<Stitch>,
    /// Face splits along datum planes, in the order they apply (after the merges).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    face_splits: Vec<FaceSplit>,
    /// Features of the design taken out for meshing, in the order they were made.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    defeatures: Vec<Defeature>,
}

/// A [`Model`] as stored, before it is checked.
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ModelData {
    #[serde(default)]
    nodes: BTreeMap<u32, Node>,
    #[serde(default)]
    elements: BTreeMap<u32, Element>,
    #[serde(default)]
    classification: geometry::ClassificationData,
    #[serde(default)]
    display_units: Option<DisplaySystem>,
    #[serde(default)]
    sets: BTreeMap<String, Set>,
    #[serde(default)]
    materials: BTreeMap<String, Material>,
    #[serde(default)]
    sections: BTreeMap<String, Section>,
    #[serde(default)]
    reference_points: BTreeMap<String, ReferencePoint>,
    #[serde(default)]
    couplings: BTreeMap<String, Coupling>,
    #[serde(default)]
    contacts: BTreeMap<String, Contact>,
    #[serde(default)]
    ties: BTreeMap<String, Tie>,
    #[serde(default)]
    masses: BTreeMap<String, PointMass>,
    #[serde(default)]
    springs: BTreeMap<String, Spring>,
    #[serde(default)]
    boundary_conditions: BTreeMap<String, BoundaryCondition>,
    #[serde(default)]
    loads: BTreeMap<String, Load>,
    #[serde(default)]
    fluid_domains: BTreeMap<String, FluidDomain>,
    #[serde(default)]
    steps: Vec<Step>,
    #[serde(default)]
    datums: BTreeMap<String, Datum>,
    #[serde(default)]
    partitions: Vec<Partition>,
    #[serde(default)]
    suppressed_datums: BTreeSet<String>,
    #[serde(default)]
    seeds: Vec<Seed>,
    #[serde(default)]
    size_controls: BTreeMap<String, SizeControl>,
    #[serde(default)]
    size_growth: Option<f64>,
    #[serde(default)]
    relaxations: BTreeMap<String, Relaxation>,
    #[serde(default)]
    mesh_plans: BTreeMap<String, MeshPlan>,
    #[serde(default)]
    face_merges: Vec<FaceMerge>,
    #[serde(default)]
    stitches: Vec<Stitch>,
    #[serde(default)]
    face_splits: Vec<FaceSplit>,
    #[serde(default)]
    defeatures: Vec<Defeature>,
}

impl TryFrom<ModelData> for Model {
    type Error = Error;

    fn try_from(d: ModelData) -> Result<Model> {
        let mut m = Model { display_units: d.display_units, ..Model::default() };
        for (label, node) in d.nodes {
            m.set_node(label, node.xyz)?;
        }
        for (label, e) in d.elements {
            m.set_element(label, e.kind, e.nodes)?;
        }
        for (i, entity) in d.classification.entities.iter().enumerate() {
            if usize::try_from(m.geometry_tag(entity)?).ok() != Some(i) {
                return Err(Error::Invalid(format!("{entity} is listed twice among the geometry tags")));
            }
        }
        for (element, tags) in d.classification.elements {
            if let Some(t) = tags.region {
                m.tag_element(element, t)?;
            }
            for (face, t) in tags.faces {
                m.tag_element_face(element, face, t)?;
            }
            for (edge, t) in tags.edges {
                m.tag_element_edge(element, edge, t)?;
            }
        }
        for (node, t) in d.classification.nodes {
            m.tag_node(node, t)?;
        }
        for (n, v) in d.sets {
            m.add_set(&n, v)?;
        }
        for (n, v) in d.materials {
            m.add_material(&n, v)?;
        }
        for (n, v) in d.sections {
            m.add_section(&n, v)?;
        }
        for (n, v) in d.reference_points {
            m.add_reference_point(&n, v)?;
        }
        for (n, v) in d.couplings {
            m.add_coupling(&n, v)?;
        }
        for (n, v) in d.contacts {
            m.add_contact(&n, v)?;
        }
        for (n, v) in d.ties {
            m.add_tie(&n, v)?;
        }
        for (n, v) in d.masses {
            m.add_mass(&n, v)?;
        }
        for (n, v) in d.springs {
            m.add_spring(&n, v)?;
        }
        for (n, v) in d.boundary_conditions {
            m.add_boundary_condition(&n, v)?;
        }
        for (n, v) in d.loads {
            m.add_load(&n, v)?;
        }
        for (n, v) in d.fluid_domains {
            m.add_fluid_domain(&n, v)?;
        }
        for s in d.steps {
            m.add_step(s, None)?;
        }
        for (n, v) in d.datums {
            m.add_datum(&n, v)?;
        }
        for n in d.suppressed_datums {
            m.set_datum_suppressed(&n, true)?;
        }
        for p in d.partitions {
            m.add_partition(p)?;
        }
        for s in d.seeds {
            m.set_seed(s)?;
        }
        for (n, v) in d.size_controls {
            m.set_size_control(&n, v)?;
        }
        m.set_size_growth(d.size_growth)?;
        for (n, v) in d.relaxations {
            m.set_relaxation(&n, v)?;
        }
        for (n, v) in d.mesh_plans {
            m.set_mesh_plan(&n, v)?;
        }
        for v in d.face_merges {
            m.add_face_merge(v)?;
        }
        for v in d.stitches {
            m.add_stitch(v)?;
        }
        for v in d.face_splits {
            m.add_face_split(v)?;
        }
        for v in d.defeatures {
            m.set_defeature(v)?;
        }
        Ok(m)
    }
}

/// Check a definition's name: 1..=[`MAX_NAME_LEN`] ASCII characters, a letter first, then
/// letters, digits, `_` or `-`.
pub fn check_name(what: &str, name: &str) -> Result<()> {
    let mut chars = name.chars();
    let ok = name.len() <= MAX_NAME_LEN
        && chars.next().is_some_and(|c| c.is_ascii_alphabetic())
        && chars.all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-');
    if ok {
        Ok(())
    } else {
        let shown: String = name.chars().take(MAX_NAME_LEN).collect();
        let rule = format!("{what} name `{shown}` must be 1 to {MAX_NAME_LEN} characters: a letter, then letters, digits, `_` or `-`");
        Err(Error::Invalid(match suggest_name(name) {
            Some(fixed) => format!("{rule} (try `{fixed}`)"),
            None => rule,
        }))
    }
}

/// The valid name closest to `name`: spaces and other characters become `_`, a leading
/// non-letter gets an `N` before it, and it is cut to [`MAX_NAME_LEN`]. `None` when nothing of it
/// is usable (empty, or no letter or digit at all).
pub fn suggest_name(name: &str) -> Option<String> {
    let trimmed = name.trim();
    if !trimmed.chars().any(|c| c.is_ascii_alphanumeric()) {
        return None;
    }
    let mut fixed: String = trimmed.chars().map(|c| if c.is_ascii_alphanumeric() || c == '-' { c } else { '_' }).collect();
    if !fixed.starts_with(|c: char| c.is_ascii_alphabetic()) {
        fixed.insert(0, 'N');
    }
    fixed.truncate(MAX_NAME_LEN);
    Some(fixed)
}

/// Check that a value is a finite number, naming it in the error.
pub(crate) fn finite(what: &str, v: f64) -> Result<f64> {
    if v.is_finite() { Ok(v) } else { Err(Error::Invalid(format!("{what} must be a finite number, got {v}"))) }
}

/// Check that a value is a finite number above zero, naming it in the error.
pub(crate) fn positive(what: &str, v: f64) -> Result<f64> {
    if v > 0.0 && v.is_finite() { Ok(v) } else { Err(Error::Invalid(format!("{what} must be a positive number, got {v}"))) }
}

/// Add a new definition under a fresh, valid name.
fn insert_new<T>(map: &mut BTreeMap<String, T>, what: &str, name: &str, value: T) -> Result<()> {
    check_name(what, name)?;
    if map.contains_key(name) {
        return Err(Error::Invalid(format!("{what} `{name}` already exists")));
    }
    if map.len() >= MAX_DEFINITIONS {
        return Err(Error::Limit(format!("a model holds at most {MAX_DEFINITIONS} of each kind of definition")));
    }
    map.insert(name.to_string(), value);
    Ok(())
}

/// Replace an existing definition.
fn replace_existing<T>(map: &mut BTreeMap<String, T>, what: &str, name: &str, value: T) -> Result<()> {
    let slot = map.get_mut(name).ok_or_else(|| Error::NotFound(format!("{what} `{name}`")))?;
    *slot = value;
    Ok(())
}

/// Prefix an error with the definition it is about.
pub(crate) fn context(what: &str, name: &str, e: Error) -> Error {
    let shown: String = name.chars().take(MAX_NAME_LEN).collect();
    match e {
        Error::Invalid(m) => Error::Invalid(format!("{what} `{shown}`: {m}")),
        Error::NotFound(m) => Error::NotFound(format!("{what} `{shown}`: {m}")),
        Error::Limit(m) => Error::Limit(format!("{what} `{shown}`: {m}")),
        Error::InUse(m) => Error::InUse(format!("{what} `{shown}`: {m}")),
        Error::Cancelled => Error::Cancelled,
    }
}

/// Remove an existing definition.
fn remove_existing<T>(map: &mut BTreeMap<String, T>, what: &str, name: &str) -> Result<T> {
    map.remove(name).ok_or_else(|| Error::NotFound(format!("{what} `{name}`")))
}

/// Fail with [`Error::InUse`] when any referrer is listed.
fn not_in_use(what: &str, name: &str, referrers: Vec<String>) -> Result<()> {
    if referrers.is_empty() { Ok(()) } else { Err(Error::InUse(format!("{what} `{name}` is used by {}", referrers.join(", ")))) }
}

fn require<T>(map: &BTreeMap<String, T>, what: &str, name: &str) -> Result<()> {
    if map.contains_key(name) { Ok(()) } else { Err(Error::NotFound(format!("{what} `{name}`"))) }
}

impl Model {
    /// The unit system the model is shown in when it sets its own; `None` follows the user's.
    pub fn display_units(&self) -> Option<&DisplaySystem> {
        self.display_units.as_ref()
    }

    /// Show the model in `units` (`None`: the user's system). Stored values do not change.
    pub fn set_display_units(&mut self, units: Option<DisplaySystem>) {
        self.display_units = units;
    }

    /// A copy with the display system and every definition but no mesh: no nodes, elements or
    /// geometry tags. A `.fieldcraft` file stores this part as text and the mesh in binary.
    pub fn without_mesh(&self) -> Model {
        Model {
            nodes: BTreeMap::new(),
            elements: BTreeMap::new(),
            classification: Classification::default(),
            display_units: self.display_units.clone(),
            sets: self.sets.clone(),
            materials: self.materials.clone(),
            sections: self.sections.clone(),
            reference_points: self.reference_points.clone(),
            couplings: self.couplings.clone(),
            contacts: self.contacts.clone(),
            ties: self.ties.clone(),
            masses: self.masses.clone(),
            springs: self.springs.clone(),
            boundary_conditions: self.boundary_conditions.clone(),
            loads: self.loads.clone(),
            fluid_domains: self.fluid_domains.clone(),
            steps: self.steps.clone(),
            datums: self.datums.clone(),
            partitions: self.partitions.clone(),
            suppressed_datums: self.suppressed_datums.clone(),
            seeds: self.seeds.clone(),
            size_controls: self.size_controls.clone(),
            size_growth: self.size_growth,
            relaxations: self.relaxations.clone(),
            mesh_plans: self.mesh_plans.clone(),
            face_merges: self.face_merges.clone(),
            stitches: self.stitches.clone(),
            face_splits: self.face_splits.clone(),
            defeatures: self.defeatures.clone(),
        }
    }

    /// The relaxation fields, by name.
    pub fn relaxations(&self) -> &BTreeMap<String, Relaxation> {
        &self.relaxations
    }

    pub fn relaxation(&self, name: &str) -> Option<&Relaxation> {
        self.relaxations.get(name)
    }

    /// Add a relaxation field, or replace the one of that name.
    pub fn set_relaxation(&mut self, name: &str, field: Relaxation) -> Result<()> {
        field.check().map_err(|e| context("relaxation", name, e))?;
        if self.relaxations.contains_key(name) {
            replace_existing(&mut self.relaxations, "relaxation", name, field)
        } else {
            insert_new(&mut self.relaxations, "relaxation", name, field)
        }
    }

    pub fn remove_relaxation(&mut self, name: &str) -> Result<Relaxation> {
        remove_existing(&mut self.relaxations, "relaxation", name)
    }

    /// The mesh plans, by name.
    pub fn mesh_plans(&self) -> &BTreeMap<String, MeshPlan> {
        &self.mesh_plans
    }

    pub fn mesh_plan(&self, name: &str) -> Option<&MeshPlan> {
        self.mesh_plans.get(name)
    }

    /// Add a mesh plan, or replace the one of that name.
    pub fn set_mesh_plan(&mut self, name: &str, plan: MeshPlan) -> Result<()> {
        plan.check().map_err(|e| context("mesh plan", name, e))?;
        if self.mesh_plans.contains_key(name) {
            replace_existing(&mut self.mesh_plans, "mesh plan", name, plan)
        } else {
            insert_new(&mut self.mesh_plans, "mesh plan", name, plan)
        }
    }

    pub fn remove_mesh_plan(&mut self, name: &str) -> Result<MeshPlan> {
        remove_existing(&mut self.mesh_plans, "mesh plan", name)
    }

    /// Counts of each kind of definition, for summaries.
    pub fn definition_counts(&self) -> [(&'static str, usize); 21] {
        [
            ("sets", self.sets.len()),
            ("materials", self.materials.len()),
            ("sections", self.sections.len()),
            ("reference_points", self.reference_points.len()),
            ("couplings", self.couplings.len()),
            ("contacts", self.contacts.len()),
            ("ties", self.ties.len()),
            ("masses", self.masses.len()),
            ("springs", self.springs.len()),
            ("boundary_conditions", self.boundary_conditions.len()),
            ("loads", self.loads.len()),
            ("fluid_domains", self.fluid_domains.len()),
            ("steps", self.steps.len()),
            ("datums", self.datums.len()),
            ("partitions", self.partitions.len()),
            ("seeds", self.seeds.len()),
            ("size_controls", self.size_controls.len()),
            ("face_merges", self.face_merges.len()),
            ("stitches", self.stitches.len()),
            ("face_splits", self.face_splits.len()),
            ("defeatures", self.defeatures.len()),
        ]
    }
}

#[cfg(test)]
mod name_tests {
    use super::*;

    #[test]
    fn a_bad_name_suggests_a_fixed_one() {
        let e = check_name("load", "Hot air").unwrap_err().to_string();
        assert!(e.contains("try `Hot_air`"), "{e}");
        assert_eq!(suggest_name("1st load"), Some("N1st_load".into()));
        assert_eq!(suggest_name("  "), None);
        assert!(check_name("load", &suggest_name(&"x ".repeat(80)).unwrap()).is_ok());
    }
}
