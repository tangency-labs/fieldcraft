//! The browser (left): the model as a tree. Geometry (the design's bodies, datums and
//! partitions) and Mesh with eye toggles, then the model's
//! definitions (sets, materials, sections, steps, loads and boundary conditions, interactions)
//! with counts (Loads & BCs and Interactions with eye toggles for their viewport symbols), and
//! the solved steps' results. Clicking a definition, datum or partition selects it (its
//! properties show on the right), a double click opens its edit form, and clicking a result step
//! shows it. Each item the viewport draws on its own (bodies, surfaces, datums, boundary
//! conditions, loads, reference points, couplings) has its own eye, and a right click on an item
//! opens its menu ([`crate::menu`]: edit, hide, isolate, active in step, delete).
//!
//! The tree also shows where the model is in the analysis workflow
//! ([`fieldcraft_engine::progress`]): each top-level node that is a setup stage carries a mark in
//! one column at the tree's left edge, so the state of every stage reads down one line (✓ done,
//! → the next stage, ⚠ when its command cannot run yet), with the reason and the place in its hover, and a node that still lacks what its
//! stage needs shows a ghost row ("＋ Add material…") that starts the providing command. Optional
//! stages (seeds and partitions) carry no mark.
//!
//! Rows are the kit's tree rows (`solidcraft_ui_kit::tree_row`), as SolidCraft's browser draws
//! them; the setup marks' column at their left, the motion arrow on a moving boundary condition's
//! icon and the ghost rows are FieldCraft's own (the kit's rows have no workflow stages).

use egui::{Color32, FontId, Pos2, Rect, Sense, Stroke, Vec2, pos2, vec2};
use fieldcraft_engine::progress::{Stage, StageState, Status};
use fieldcraft_engine::schema::Collection;
use fieldcraft_model::Model;
use serde::Serialize;
use serde_json::{Value, json};

use crate::ModelApp;
use crate::menu::{self, Action, Target};
use crate::theme::Tokens;
use solidcraft_ui_kit as kit;
use solidcraft_ui_kit::BackingTooltipExt as _;

/// The browser's width.
pub const WIDTH: f32 = 250.0;
/// The kit's indent per tree level.
const INDENT: f32 = 14.0;
/// The state marks' column at the tree's left edge (every row keeps it, so the tree lines up).
const MARK_W: f32 = 14.0;

/// The browser's top-level nodes below the model, in order, with the setup stages each shows
/// (`ui.inspect`'s `browser`). Seeds and partitions ([`Stage::Prepare`]) are optional and sit
/// under Geometry and Mesh without a mark.
pub const NODES: [(&str, &[Stage]); 10] = [
    ("geometry", &[Stage::Geometry]),
    ("mesh", &[Stage::Mesh]),
    ("sets", &[Stage::Sets]),
    ("materials", &[Stage::Material]),
    ("sections", &[Stage::Section]),
    ("steps", &[Stage::Step]),
    ("loads", &[Stage::BoundaryConditions, Stage::Loads]),
    ("interactions", &[]),
    ("fluids", &[]),
    ("results", &[Stage::Solve, Stage::Results]),
];

fn node_stages(node: &str) -> &'static [Stage] {
    NODES.iter().find(|(n, _)| *n == node).map_or(&[], |(_, st)| st)
}

/// A stage's hover: what it needs and where; `None` once done or when optional.
fn stage_tip(s: &StageState) -> Option<String> {
    match s.status {
        Status::Next => Some(format!("Next: {} — {}", s.label, s.place)),
        Status::Blocked | Status::Outdated => Some(s.detail.clone().unwrap_or_else(|| s.place.clone())),
        Status::Todo => Some(format!("{} still to do — {}", s.label, s.place)),
        Status::Done | Status::Optional => None,
    }
}

/// A node's state and hover from its stages: the most pressing of them (next, blocked, to do),
/// else done; `None` for nodes that are no setup stage.
fn node_state(stages: &[StageState], node: &str) -> Option<(Status, String)> {
    let own: Vec<&StageState> = stages.iter().filter(|s| node_stages(node).contains(&s.stage)).collect();
    let status = [Status::Outdated, Status::Next, Status::Blocked, Status::Todo, Status::Done, Status::Optional]
        .into_iter()
        .find(|st| own.iter().any(|s| s.status == *st))?;
    let tip: Vec<String> = own.iter().filter_map(|s| stage_tip(s)).collect();
    let tip = if tip.is_empty() { "Done".to_string() } else { tip.join("\n") };
    Some((status, tip))
}

/// A ghost row's text: what its stage's command adds.
fn ghost_label(stage: Stage) -> &'static str {
    match stage {
        Stage::Geometry => "Open design…",
        Stage::Prepare => "Seed edges…",
        Stage::Mesh => "Generate mesh…",
        Stage::Material => "Add material…",
        Stage::Section => "Add section…",
        Stage::Sets => "Add set…",
        Stage::BoundaryConditions => "Add boundary condition…",
        Stage::Loads => "Add load…",
        Stage::Step => "Add step…",
        Stage::Solve => "Solve…",
        Stage::Results => "Show results…",
    }
}

/// Where a stage's ghost row was drawn (`ui.inspect`'s `layout`).
pub fn ghost_region(stage: Stage) -> &'static str {
    match stage {
        Stage::Geometry => "browser/add/geometry",
        Stage::Prepare => "browser/add/prepare",
        Stage::Mesh => "browser/add/mesh",
        Stage::Material => "browser/add/material",
        Stage::Section => "browser/add/section",
        Stage::Sets => "browser/add/sets",
        Stage::BoundaryConditions => "browser/add/boundary_conditions",
        Stage::Loads => "browser/add/loads",
        Stage::Step => "browser/add/step",
        Stage::Solve => "browser/add/solve",
        Stage::Results => "browser/add/results",
    }
}

/// The stages of `node` still to do: each shows a ghost row, except a blocked one behind an
/// earlier stage of the same node still to do (Results offers Solve, not Show results, until a
/// solve).
fn ghosts<'s>(stages: &'s [StageState], node: &str) -> Vec<&'s StageState> {
    let mut pending = false;
    let mut out = Vec::new();
    for s in stages.iter().filter(|s| node_stages(node).contains(&s.stage)) {
        let todo = matches!(s.status, Status::Next | Status::Todo | Status::Blocked);
        if todo && !(pending && s.status == Status::Blocked) {
            out.push(s);
        }
        pending |= todo;
    }
    out
}

/// Each top-level node's setup state for `ui.inspect`: `[{node, state, tip, add: [{stage,
/// label, command, status}]}]`; `state` is null for nodes that are no setup stage.
pub fn describe(app: &mut ModelApp) -> Value {
    let stages = app.setup_stages();
    Value::Array(
        NODES
            .iter()
            .map(|(node, _)| {
                let state = node_state(&stages, node);
                let add: Vec<Value> = ghosts(&stages, node)
                    .into_iter()
                    .map(|s| json!({"stage": s.stage, "label": ghost_label(s.stage), "command": s.command, "status": s.status}))
                    .collect();
                json!({"node": node, "state": state.as_ref().map(|(st, _)| st), "tip": state.map(|(_, tip)| tip), "add": add})
            })
            .collect(),
    )
}

/// A kind of model definition.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Def {
    Set,
    Material,
    Section,
    Step,
    BoundaryCondition,
    Load,
    ReferencePoint,
    Coupling,
    Contact,
    /// A tie between two surfaces (`interaction.tie`).
    Tie,
    FluidDomain,
    /// A datum plane, axis or point.
    Datum,
    Partition,
    /// A point mass (`model.mass.create`).
    Mass,
    /// A spring or dashpot (`model.spring.create`).
    Spring,
    /// A local element size (`mesh.size.set`).
    SizeControl,
    /// Faces meshed as one (`geometry.faces.merge`), or a face absorbed into a neighbour
    /// (`geometry.face.absorb`, a merge that keeps its neighbour's name).
    FaceMerge,
    /// A face split along a datum plane (`geometry.face.split`).
    FaceSplit,
    /// An imported surface's gaps closed (`geometry.stitch`).
    Stitch,
    /// Features of the design taken out for meshing (`geometry.defeature`).
    Defeature,
    /// The temperature or stress the nodes start from (`model.initial_field.map`).
    InitialField,
    /// A body's saved hex-meshing plan (`mesh.plan`), applied by `mesh.plan.apply`.
    MeshPlan,
    /// A relaxation field kept by `mesh.node.relax`, applied by `mesh.node.relax.apply`.
    Relaxation,
}

impl Def {
    pub const ALL: [Def; 23] = [
        Def::Set,
        Def::Material,
        Def::Section,
        Def::Step,
        Def::BoundaryCondition,
        Def::Load,
        Def::ReferencePoint,
        Def::Coupling,
        Def::Contact,
        Def::Tie,
        Def::FluidDomain,
        Def::Datum,
        Def::Partition,
        Def::Mass,
        Def::Spring,
        Def::SizeControl,
        Def::FaceMerge,
        Def::FaceSplit,
        Def::Stitch,
        Def::Defeature,
        Def::InitialField,
        Def::MeshPlan,
        Def::Relaxation,
    ];

    /// The kind `model.rename` takes for it; `None` for a partition (its cut faces are named
    /// after it) and the meshing helpers `model.rename` does not take.
    pub fn rename_kind(self) -> Option<&'static str> {
        Some(match self {
            Def::Set => "set",
            Def::Material => "material",
            Def::Section => "section",
            Def::Step => "step",
            Def::BoundaryCondition => "boundary_condition",
            Def::Load => "load",
            Def::ReferencePoint => "reference_point",
            Def::Coupling => "coupling",
            Def::Contact => "contact",
            Def::Tie => "tie",
            Def::FluidDomain => "fluid_domain",
            Def::Datum => "datum",
            Def::Mass => "point_mass",
            Def::Spring => "spring",
            Def::Partition
            | Def::SizeControl
            | Def::FaceMerge
            | Def::FaceSplit
            | Def::Stitch
            | Def::Defeature
            | Def::InitialField
            | Def::MeshPlan
            | Def::Relaxation => return None,
        })
    }

    /// The kind of setup symbols a definition of this kind draws, if any.
    pub fn symbols(self) -> Option<fieldcraft_render::glyphs::Definition> {
        use fieldcraft_render::glyphs::Definition;
        Some(match self {
            Def::BoundaryCondition => Definition::BoundaryCondition,
            Def::Load => Definition::Load,
            Def::ReferencePoint => Definition::ReferencePoint,
            Def::Coupling => Definition::Coupling,
            Def::Contact => Definition::Contact,
            Def::Tie => Definition::Tie,
            Def::Step => Definition::Step,
            Def::Set
            | Def::Material
            | Def::Section
            | Def::FluidDomain
            | Def::Datum
            | Def::Partition
            | Def::Mass
            | Def::Spring
            | Def::SizeControl
            | Def::FaceMerge
            | Def::FaceSplit
            | Def::Stitch
            | Def::Defeature
            | Def::InitialField
            | Def::MeshPlan
            | Def::Relaxation => return None,
        })
    }

    /// The reference fields that name one of these.
    pub fn collection(self) -> Collection {
        match self {
            Def::Set => Collection::Set,
            Def::Material => Collection::Material,
            Def::Section => Collection::Section,
            Def::Step => Collection::Step,
            Def::BoundaryCondition => Collection::BoundaryCondition,
            Def::Load => Collection::Load,
            Def::ReferencePoint => Collection::ReferencePoint,
            Def::Coupling => Collection::Coupling,
            Def::Contact => Collection::Contact,
            Def::Tie => Collection::Tie,
            Def::FluidDomain => Collection::FluidDomain,
            Def::Datum => Collection::Datum,
            Def::Partition => Collection::Partition,
            Def::Mass => Collection::Mass,
            Def::Spring => Collection::Spring,
            Def::SizeControl => Collection::SizeControl,
            Def::FaceMerge => Collection::FaceMerge,
            Def::FaceSplit => Collection::FaceSplit,
            Def::Stitch => Collection::Stitch,
            Def::Defeature => Collection::Defeature,
            Def::InitialField => Collection::InitialField,
            Def::MeshPlan => Collection::MeshPlan,
            Def::Relaxation => Collection::Relaxation,
        }
    }

    /// Its short name: `model.<stem>.create|edit|delete` for the model's definitions, and the
    /// row ids' kind.
    pub fn stem(self) -> &'static str {
        match self {
            Def::Set => "set",
            Def::Material => "material",
            Def::Section => "section",
            Def::Step => "step",
            Def::BoundaryCondition => "boundary_condition",
            Def::Load => "load",
            Def::ReferencePoint => "reference_point",
            Def::Coupling => "coupling",
            Def::Contact => "contact",
            Def::Tie => "tie",
            Def::FluidDomain => "fluid_domain",
            Def::Datum => "datum",
            Def::Partition => "partition",
            Def::Mass => "mass",
            Def::Spring => "spring",
            Def::SizeControl => "size_control",
            Def::FaceMerge => "face_merge",
            Def::FaceSplit => "face_split",
            Def::Stitch => "stitch",
            Def::Defeature => "defeature",
            Def::InitialField => "initial_field",
            Def::MeshPlan => "mesh_plan",
            Def::Relaxation => "relaxation",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Def::Set => "Set",
            Def::Material => "Material",
            Def::Section => "Section",
            Def::Step => "Step",
            Def::BoundaryCondition => "Boundary Condition",
            Def::Load => "Load",
            Def::ReferencePoint => "Reference Point",
            Def::Coupling => "Coupling",
            Def::Contact => "Contact",
            Def::Tie => "Tie",
            Def::FluidDomain => "Fluid Domain",
            Def::Datum => "Datum",
            Def::Partition => "Partition",
            Def::Mass => "Point Mass",
            Def::Spring => "Spring or Dashpot",
            Def::SizeControl => "Size Control",
            Def::FaceMerge => "Face Merge",
            Def::FaceSplit => "Face Split",
            Def::Stitch => "Stitch",
            Def::Defeature => "Defeature",
            Def::InitialField => "Initial Field",
            Def::MeshPlan => "Mesh Plan",
            Def::Relaxation => "Relaxation Field",
        }
    }

    /// Its icon: the catalogue entry that creates it.
    pub fn icon(self) -> &'static str {
        match self {
            Def::Set => "model.set.create",
            Def::Material => "model.material.create",
            Def::Section => "model.section.create",
            Def::Step => "model.step.create",
            Def::BoundaryCondition => "model.boundary_condition.create",
            Def::Load => "model.load.create",
            Def::ReferencePoint => "model.reference_point.create",
            Def::Coupling => "model.coupling.create",
            Def::Contact => "interaction.contact.pair",
            Def::Tie => "interaction.tie",
            Def::FluidDomain => "model.fluid_domain.create",
            Def::Datum => "geometry.datum",
            Def::Partition => "geometry.partition",
            Def::Mass => "model.mass.create",
            Def::Spring => "model.spring.create",
            Def::SizeControl => "mesh.size.set",
            Def::FaceMerge => "geometry.faces.merge",
            Def::FaceSplit => "geometry.face.split",
            Def::Stitch => "geometry.stitch",
            Def::Defeature => "geometry.defeature",
            // Mapped from the Loads entry (the catalogue has no entry of its own for it).
            Def::InitialField => "model.load.create",
            Def::MeshPlan => "mesh.plan",
            Def::Relaxation => "mesh.node.relax",
        }
    }

    /// The command of this kind with `verb` (`create`, `edit`, `delete`); contacts are
    /// interactions, created as pairs or general contact; datums and partitions are geometry
    /// (a datum's edit command depends on its kind: [`Def::edit`]). A size control and a
    /// defeature are edited by making them again under their name, and size controls are
    /// deleted by `mesh.size.clear` ([`Def::delete_params`]). Face merges, splits, stitches,
    /// initial fields, mesh plans and relaxation fields cannot be edited ([`Def::edit_command`]);
    /// plans and fields are applied instead ([`Def::apply`]).
    pub fn command(self, verb: &str) -> String {
        match (self, verb) {
            (Def::SizeControl, "create" | "edit") => "mesh.size.set".into(),
            (Def::SizeControl, "delete") => "mesh.size.clear".into(),
            (Def::Defeature, "create" | "edit") => "geometry.defeature".into(),
            (Def::FaceMerge | Def::FaceSplit | Def::Stitch | Def::MeshPlan | Def::Relaxation, "create") => self.icon().into(),
            (Def::FaceMerge | Def::FaceSplit | Def::Stitch | Def::Defeature | Def::SizeControl | Def::MeshPlan | Def::Relaxation, _) => {
                format!("{}.{verb}", self.icon())
            }
            (Def::InitialField, "create") => "model.initial_field.map".into(),
            (Def::Contact, "create") => "interaction.contact.pair".into(),
            (Def::Contact, _) => format!("interaction.contact.{verb}"),
            (Def::Tie, "create") => "interaction.tie".into(),
            (Def::Tie, _) => format!("interaction.tie.{verb}"),
            (Def::Datum | Def::Partition, "create") => self.icon().into(),
            (Def::Datum | Def::Partition, _) => format!("geometry.{}.{verb}", self.stem()),
            _ => format!("model.{}.{verb}", self.stem()),
        }
    }

    /// The delete command's parameters for definition `name`: `{"name": name}`, or
    /// `{"names": [name]}` for a size control's `mesh.size.clear`.
    pub fn delete_params(self, name: &str) -> Value {
        match self {
            Def::SizeControl => json!({"names": [name]}),
            _ => json!({"name": name}),
        }
    }

    /// The command that applies definition `name` to the mesh, filled in (`mesh.plan.apply`,
    /// `mesh.node.relax.apply`); `None` for kinds that are not applied, or when it is gone.
    pub fn apply(self, m: &Model, name: &str) -> Option<(&'static str, Value)> {
        match self {
            Def::MeshPlan => m.mesh_plan(name).map(|_| ("mesh.plan.apply", json!({"plan": name}))),
            Def::Relaxation => m.relaxation(name).map(|_| ("mesh.node.relax.apply", json!({"field": name}))),
            _ => None,
        }
    }

    /// Whether definitions of this kind can be edited at all (by an edit command, or by making
    /// them again under their name).
    fn editable(self) -> bool {
        !matches!(self, Def::Partition | Def::FaceMerge | Def::FaceSplit | Def::Stitch | Def::InitialField | Def::MeshPlan | Def::Relaxation)
    }

    /// The command that edits definition `name` (without building its parameters); `None` when
    /// it is gone or cannot be edited.
    pub fn edit_command(self, m: &Model, name: &str) -> Option<String> {
        match self {
            _ if !self.editable() => None,
            Def::Datum => m.datum(name).map(|d| fieldcraft_engine::cmd::datum_edit_params(name, d).0.to_string()),
            _ => self.exists(m, name).then(|| self.command("edit")),
        }
    }

    /// The command that edits definition `name` and its parameters as it is, in `units` (its
    /// edit form, filled in); `None` when it is gone or cannot be edited (a partition).
    pub fn edit(self, m: &Model, name: &str, units: &fieldcraft_model::DisplaySystem) -> Option<(String, Value)> {
        let (command, params) = match self {
            _ if !self.editable() => return None,
            Def::Datum => m.datum(name).map(|d| {
                let (command, params) = fieldcraft_engine::cmd::datum_edit_params(name, d);
                (command.to_string(), params)
            })?,
            _ => (self.command("edit"), definition(m, self, name)?),
        };
        // The model is in SI, the form in the display units: an edit that changes nothing must
        // give back the same values.
        let shown = fieldcraft_engine::units::shown_params(&command, &params, units).ok()?;
        Some((command, shown))
    }

    /// The viewport item it is, for the kinds the viewport draws one by one (and can hide).
    pub fn item_kind(self) -> Option<fieldcraft_engine::ItemKind> {
        use fieldcraft_engine::ItemKind as K;
        match self {
            Def::BoundaryCondition => Some(K::BoundaryCondition),
            Def::Load => Some(K::Load),
            Def::ReferencePoint => Some(K::ReferencePoint),
            Def::Coupling => Some(K::Coupling),
            Def::Datum => Some(K::Datum),
            _ => None,
        }
    }

    /// The definitions that refer to `name`, as `kind \`name\`` strings: what refuses its
    /// deletion (a datum's partitions are deleted with it instead, when confirmed).
    pub fn users(self, m: &Model, name: &str) -> Vec<String> {
        match self {
            Def::Set => m.set_users(name),
            Def::Material => m.material_users(name),
            Def::BoundaryCondition => m.boundary_condition_users(name),
            Def::Load => m.load_users(name),
            Def::ReferencePoint => m.reference_point_users(name),
            Def::FluidDomain => m.fluid_domain_users(name),
            Def::Datum => m.datum_users(name),
            Def::Section
            | Def::Step
            | Def::Coupling
            | Def::Contact
            | Def::Tie
            | Def::Partition
            | Def::Mass
            | Def::Spring
            | Def::SizeControl
            | Def::FaceMerge
            | Def::FaceSplit
            | Def::Stitch
            | Def::Defeature
            | Def::InitialField
            | Def::MeshPlan
            | Def::Relaxation => Vec::new(),
        }
    }

    /// How many definitions of this kind the model has.
    pub fn count(self, m: &Model) -> usize {
        match self {
            Def::Set => m.sets().len(),
            Def::Material => m.materials().len(),
            Def::Section => m.sections().len(),
            Def::Step => m.steps().len(),
            Def::BoundaryCondition => m.boundary_conditions().len(),
            Def::Load => m.loads().len(),
            Def::ReferencePoint => m.reference_points().len(),
            Def::Coupling => m.couplings().len(),
            Def::Contact => m.contacts().len(),
            Def::Tie => m.ties().len(),
            Def::FluidDomain => m.fluid_domains().len(),
            Def::Datum => m.datums().len(),
            Def::Partition => m.partitions().len(),
            Def::Mass => m.masses().len(),
            Def::Spring => m.springs().len(),
            Def::SizeControl => m.size_controls().len(),
            Def::FaceMerge => m.face_merges().len(),
            Def::FaceSplit => m.face_splits().len(),
            Def::Stitch => m.stitches().len(),
            Def::Defeature => m.defeatures().len(),
            Def::InitialField => m.initial_fields().len(),
            Def::MeshPlan => m.mesh_plans().len(),
            Def::Relaxation => m.relaxations().len(),
        }
    }

    /// Whether the model has a definition of this kind called `name`.
    pub fn exists(self, m: &Model, name: &str) -> bool {
        match self {
            Def::Set => m.sets().contains_key(name),
            Def::Material => m.materials().contains_key(name),
            Def::Section => m.sections().contains_key(name),
            Def::Step => m.step(name).is_some(),
            Def::BoundaryCondition => m.boundary_conditions().contains_key(name),
            Def::Load => m.loads().contains_key(name),
            Def::ReferencePoint => m.reference_points().contains_key(name),
            Def::Coupling => m.couplings().contains_key(name),
            Def::Contact => m.contacts().contains_key(name),
            Def::Tie => m.tie(name).is_some(),
            Def::FluidDomain => m.fluid_domain(name).is_some(),
            Def::Datum => m.datum(name).is_some(),
            Def::Partition => m.partitions().iter().any(|p| p.name == name),
            Def::Mass => m.mass(name).is_some(),
            Def::Spring => m.spring(name).is_some(),
            Def::SizeControl => m.size_controls().contains_key(name),
            Def::FaceMerge => m.face_merges().iter().any(|f| f.name == name),
            Def::FaceSplit => m.face_splits().iter().any(|f| f.name == name),
            Def::Stitch => m.stitches().iter().any(|f| f.name == name),
            Def::Defeature => m.defeatures().iter().any(|f| f.name == name),
            Def::InitialField => m.initial_field(name).is_some(),
            Def::MeshPlan => m.mesh_plan(name).is_some(),
            Def::Relaxation => m.relaxation(name).is_some(),
        }
    }

    /// The names of this kind's definitions, in the model's order.
    pub fn names(self, m: &Model) -> Vec<String> {
        match self {
            Def::Set => m.sets().keys().cloned().collect(),
            Def::Material => m.materials().keys().cloned().collect(),
            Def::Section => m.sections().keys().cloned().collect(),
            Def::Step => m.steps().iter().map(|s| s.name.clone()).collect(),
            Def::BoundaryCondition => m.boundary_conditions().keys().cloned().collect(),
            Def::Load => m.loads().keys().cloned().collect(),
            Def::ReferencePoint => m.reference_points().keys().cloned().collect(),
            Def::Coupling => m.couplings().keys().cloned().collect(),
            Def::Contact => m.contacts().keys().cloned().collect(),
            Def::Tie => m.ties().keys().cloned().collect(),
            Def::FluidDomain => m.fluid_domains().keys().cloned().collect(),
            Def::Datum => m.datums().keys().cloned().collect(),
            Def::Partition => m.partitions().iter().map(|p| p.name.clone()).collect(),
            Def::Mass => m.masses().keys().cloned().collect(),
            Def::Spring => m.springs().keys().cloned().collect(),
            Def::SizeControl => m.size_controls().keys().cloned().collect(),
            Def::FaceMerge => m.face_merges().iter().map(|f| f.name.clone()).collect(),
            Def::FaceSplit => m.face_splits().iter().map(|f| f.name.clone()).collect(),
            Def::Stitch => m.stitches().iter().map(|f| f.name.clone()).collect(),
            Def::Defeature => m.defeatures().iter().map(|f| f.name.clone()).collect(),
            Def::InitialField => m.initial_fields().keys().cloned().collect(),
            Def::MeshPlan => m.mesh_plans().keys().cloned().collect(),
            Def::Relaxation => m.relaxations().keys().cloned().collect(),
        }
    }
}

/// A definition as JSON (its command parameters' shape, with its name); `None` when it is gone.
pub fn definition(m: &Model, kind: Def, name: &str) -> Option<Value> {
    let v = match kind {
        // Kept in the model's own shape; as their commands take them.
        Def::SizeControl => Ok(size_params(m.size_controls().get(name)?)),
        Def::FaceMerge => Ok(merge_params(m.face_merges().iter().find(|f| f.name == name)?)),
        Def::FaceSplit => {
            let f = m.face_splits().iter().find(|f| f.name == name)?;
            Ok(json!({"face": {"body": f.body, "face": f.face}, "datum": f.datum}))
        }
        Def::Stitch => serde_json::to_value(m.stitches().iter().find(|f| f.name == name)?),
        Def::Defeature => serde_json::to_value(m.defeatures().iter().find(|f| f.name == name)?),
        // What they are, not their values (a value per node, a cell per relaxed element): their
        // making commands' parameters that they keep.
        Def::InitialField => {
            let f = m.initial_field(name)?;
            Ok(json!({"quantity": f.kind(), "nodes": f.len(), "source": f.source()}))
        }
        Def::MeshPlan => {
            let p = m.mesh_plan(name)?;
            Ok(json!({"body": p.body, "layers": p.layers}))
        }
        Def::Relaxation => {
            let r = m.relaxation(name)?;
            Ok(json!({"bodies": r.bodies, "reapply": r.reapply}))
        }
        Def::Mass => serde_json::to_value(m.mass(name)?),
        Def::Spring => serde_json::to_value(m.spring(name)?),
        Def::Set => serde_json::to_value(m.sets().get(name)?),
        Def::Material => serde_json::to_value(m.materials().get(name)?),
        Def::Section => serde_json::to_value(m.sections().get(name)?),
        Def::Step => serde_json::to_value(m.step(name)?),
        Def::BoundaryCondition => serde_json::to_value(m.boundary_conditions().get(name)?),
        Def::Load => serde_json::to_value(m.loads().get(name)?),
        Def::ReferencePoint => serde_json::to_value(m.reference_points().get(name)?),
        Def::Coupling => serde_json::to_value(m.couplings().get(name)?),
        Def::Contact => serde_json::to_value(m.contacts().get(name)?),
        Def::Tie => serde_json::to_value(m.tie(name)?),
        Def::FluidDomain => serde_json::to_value(m.fluid_domain(name)?),
        Def::Datum => serde_json::to_value(m.datum(name)?),
        Def::Partition => serde_json::to_value(m.partitions().iter().find(|p| p.name == name)?),
    }
    .ok()?;
    let mut v = if v.is_object() { v } else { json!({"value": v}) };
    if let Some(o) = v.as_object_mut() {
        o.insert("name".into(), json!(name));
    }
    Some(v)
}

/// `mesh.size.set`'s parameters for a size control (without its name).
fn size_params(c: &fieldcraft_model::SizeControl) -> Value {
    use fieldcraft_model::{Entity, SizeRegion};
    let region = match &c.region {
        // A control is on faces or on bodies (`mesh.size.set` takes one or the other).
        SizeRegion::Geometry(refs) if refs.iter().all(|r| matches!(r.entity, Entity::Body)) => {
            json!({"bodies": refs.iter().map(|r| r.body.clone()).collect::<Vec<_>>()})
        }
        SizeRegion::Geometry(refs) => json!({"faces": refs}),
        SizeRegion::Box { min, max } => json!({"box": {"min": min, "max": max}}),
        SizeRegion::Sphere { center, radius } => json!({"sphere": {"center": center, "radius": radius}}),
    };
    let mut v = region;
    v["size"] = json!(c.size);
    v
}

/// A face merge's parameters: `geometry.faces.merge`'s, or `geometry.face.absorb`'s for an absorb
/// (the face, and the neighbour it went into).
fn merge_params(f: &fieldcraft_model::FaceMerge) -> Value {
    match &f.into {
        Some(into) => {
            let face = f.faces.iter().find(|n| *n != into).cloned().unwrap_or_default();
            json!({"face": {"body": f.body, "face": face}, "into": into})
        }
        None => json!({"faces": f.faces.iter().map(|n| json!({"body": f.body, "face": n})).collect::<Vec<_>>()}),
    }
}

/// One row's look.
#[derive(Default)]
struct Row<'a> {
    depth: usize,
    /// Some(open) for rows that fold.
    fold: Option<bool>,
    /// Some(visible) for rows with an eye.
    eye: Option<bool>,
    icon: Option<&'static str>,
    /// A motion arrow over the icon's corner: a boundary condition that moves the part.
    moves: bool,
    /// A colour swatch instead of an icon (element kinds).
    swatch: Option<Color32>,
    label: &'a str,
    /// The label is a name someone gave (shown as typed, never translated).
    name: bool,
    /// Dim text after the label (a count).
    note: Option<String>,
    selected: bool,
    dim: bool,
    /// The setup state mark in the left column (top-level nodes).
    mark: Option<Status>,
    /// The hover text (a setup state's reason and place).
    tip: Option<String>,
}

/// What happened on a row.
#[derive(Default)]
struct RowResp {
    clicked: bool,
    double: bool,
    fold: bool,
    eye: bool,
    /// The row's response, for its context menu (`None` when the row was not drawn).
    response: Option<egui::Response>,
}

/// The setup marks' column, then the kit's tree row.
fn draw_row(app: &mut ModelApp, ui: &mut egui::Ui, id: egui::Id, row: &Row) -> RowResp {
    let t = Tokens::get(app.ui.dark);
    let kit_row = kit::Row {
        depth: row.depth,
        fold: row.fold,
        eye: row.eye,
        icon: row.icon.unwrap_or(""),
        label: row.label,
        name: row.name,
        selected: row.selected,
        dim: row.dim,
        tip: row.tip.as_deref(),
        swatch: if row.icon.is_none() { row.swatch } else { None },
        count: row.note.as_deref().unwrap_or(""),
        ..Default::default()
    };
    let r = ui
        .horizontal(|ui| {
            ui.spacing_mut().item_spacing.x = 0.0;
            let (column, _) = ui.allocate_exact_size(vec2(MARK_W, kit::tree_row_rect(ui, &kit_row).height()), Sense::hover());
            if let Some(status) = row.mark
                && ui.is_rect_visible(column)
            {
                paint_mark(ui.painter(), Rect::from_center_size(column.center() + vec2(1.0, 0.0), Vec2::splat(9.0)), status, &t);
            }
            kit::tree_row(ui, id, &kit_row)
        })
        .inner;
    let Some(resp) = r.resp else { return RowResp::default() };
    if row.moves
        && let (Some(rect), Some(label_left)) = (r.rect, r.label_left)
    {
        // An arrow out of the icon's lower right corner, on a disc of the row's background. The
        // kit gives where the label starts; its 16 pt icon ends 5 pt before it (`tree_row`).
        let icon = Rect::from_min_size(pos2(label_left - 21.0, rect.center().y - 8.0), Vec2::splat(16.0));
        let c = icon.right_bottom() - vec2(3.0, 3.0);
        let k = solidcraft_theme::Tokens::get();
        let disc = if row.selected {
            solidcraft_theme::TabAccent::current(ui.ctx()).soft
        } else if r.hovered {
            k.hover
        } else {
            t.panel
        };
        ui.painter().circle_filled(c, 4.5, disc);
        let s = Stroke::new(1.4, t.accent);
        ui.painter().line_segment([c + vec2(-3.0, -3.0), c + vec2(2.5, 2.5)], s);
        ui.painter().line_segment([c + vec2(2.5, 2.5), c + vec2(-0.5, 2.5)], s);
        ui.painter().line_segment([c + vec2(2.5, 2.5), c + vec2(2.5, -0.5)], s);
    }
    resp.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, true, row.label));
    RowResp { clicked: r.clicked, double: r.double, fold: r.fold, eye: r.eye, response: Some(resp) }
}

/// A row for one item: its eye when the viewport draws it on its own, selection by a click and
/// its edit form by a double click (definitions, datums, partitions), and its menu by a right
/// click.
fn item_row(app: &mut ModelApp, ui: &mut egui::Ui, id: egui::Id, mut row: Row, target: Target, acts: &mut Vec<Act>) {
    row.name = true;
    if let Some(item) = target.item() {
        row.eye = Some(app.session.shows(item.kind, &item.name));
    }
    match &target {
        Target::Def(kind, name) => row.selected = app.selected.as_ref().is_some_and(|(k, n)| k == kind && n == name),
        Target::Seed(edge) => row.selected = app.selected_seed.as_ref() == Some(edge),
        Target::View(name) => row.selected = app.session.view.current.as_ref() == Some(name),
        Target::Body(_) | Target::Surface(_) | Target::Render(_) | Target::RigidWall { .. } => {}
    }
    let r = draw_row(app, ui, id, &row);
    // A right click selects the row its menu acts on (another item's highlight would say the
    // menu is about that one).
    if r.response.as_ref().is_some_and(egui::Response::secondary_clicked) {
        acts.push(Act::Focus(target.clone()));
    }
    if let Some(resp) = &r.response {
        crate::frame::context_menu(resp, |ui| {
            if let Some(a) = menu::show(app, ui, &target) {
                acts.push(Act::Menu(target.clone(), a));
            }
        });
    }
    if r.eye {
        acts.push(Act::Menu(target, if row.eye == Some(true) { Action::Hide } else { Action::Show }));
    } else if r.double && matches!(target, Target::View(_)) {
        acts.push(Act::Menu(target, Action::Restore));
    } else if r.double && matches!(target, Target::Def(..) | Target::Seed(_) | Target::Render(_) | Target::RigidWall { .. }) {
        acts.push(Act::Menu(target, Action::Edit));
    } else if r.clicked {
        match target {
            Target::Def(kind, name) => acts.push(Act::Select(kind, name)),
            Target::Seed(edge) => acts.push(Act::SelectSeed(edge)),
            Target::Body(_) | Target::Surface(_) | Target::View(_) | Target::Render(_) | Target::RigidWall { .. } => {}
        }
    }
}

/// ✓ done (success), → next (accent), ⚠ blocked (muted), ⚠ out of date (warm); nothing for the rest.
fn paint_mark(p: &egui::Painter, r: Rect, status: Status, t: &Tokens) {
    let at = |x: f32, y: f32| Pos2::new(r.left() + x * r.width(), r.top() + y * r.height());
    match status {
        Status::Done => {
            let s = Stroke::new(1.6, t.ok);
            p.line_segment([at(0.05, 0.55), at(0.4, 0.88)], s);
            p.line_segment([at(0.4, 0.88), at(0.95, 0.15)], s);
        }
        Status::Next => {
            let s = Stroke::new(1.6, t.accent);
            p.line_segment([at(0.0, 0.5), at(0.95, 0.5)], s);
            p.line_segment([at(0.55, 0.15), at(0.95, 0.5)], s);
            p.line_segment([at(0.55, 0.85), at(0.95, 0.5)], s);
        }
        Status::Blocked => {
            let s = Stroke::new(1.2, t.text_dim);
            p.add(egui::Shape::closed_line(vec![at(0.5, 0.0), at(1.05, 1.0), at(-0.05, 1.0)], s));
            p.line_segment([at(0.5, 0.35), at(0.5, 0.65)], s);
            p.circle_filled(at(0.5, 0.84), 0.7, t.text_dim);
        }
        Status::Outdated => {
            let s = Stroke::new(1.3, t.warm);
            p.add(egui::Shape::closed_line(vec![at(0.5, 0.0), at(1.05, 1.0), at(-0.05, 1.0)], s));
            p.line_segment([at(0.5, 0.35), at(0.5, 0.65)], s);
            p.circle_filled(at(0.5, 0.84), 0.8, t.warm);
        }
        Status::Todo | Status::Optional => {}
    }
}

/// A ghost row under a node that lacks what its stage needs: "＋ Add material…", muted, the next
/// stage's in the accent colour. Clicking it (or Enter/Space when it has the keyboard focus)
/// starts the stage's command; true when activated.
fn draw_ghost(app: &mut ModelApp, ui: &mut egui::Ui, s: &StageState) -> bool {
    let t = Tokens::get(app.ui.dark);
    let (r, _) = ui.allocate_exact_size(vec2(ui.available_width(), kit::TREE_ROW_HEIGHT), Sense::hover());
    app.place(ghost_region(s.stage), r);
    if !ui.is_rect_visible(r) {
        return false;
    }
    let resp = ui.interact(r, egui::Id::unique(("fc_row_add", ghost_region(s.stage))), Sense::click());
    if resp.hovered() {
        ui.painter().rect_filled(r, 3.0, t.hover);
    }
    if resp.has_focus() {
        ui.painter().rect_stroke(r.shrink(0.5), 3.0, Stroke::new(1.0, t.accent), egui::StrokeKind::Inside);
    }
    let col = if s.status == Status::Next { t.accent_ink } else { t.text_dim };
    // Aligned with the node's children: past the fold arrow and the eye's gap (as the kit's rows
    // place their icons), the plus in the icon's place.
    let x = r.left() + MARK_W + 2.0 + 2.0 * INDENT + 14.0 + solidcraft_theme::Tokens::get().space_10;
    let c = pos2(x + 8.0, r.center().y);
    let plus = Stroke::new(1.3, col);
    ui.painter().line_segment([c - vec2(4.0, 0.0), c + vec2(4.0, 0.0)], plus);
    ui.painter().line_segment([c - vec2(0.0, 4.0), c + vec2(0.0, 4.0)], plus);
    let label = ghost_label(s.stage);
    let font = FontId::proportional(solidcraft_theme::Tokens::get().font_label);
    ui.painter().with_clip_rect(r).text(pos2(x + 21.0, r.center().y), egui::Align2::LEFT_CENTER, label, font, col);
    let resp = resp.backing_on_hover_text(stage_tip(s).unwrap_or_default());
    resp.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, true, label));
    resp.clicked()
}

/// What the browser asks for after drawing (it cannot run commands while it borrows the model).
enum Act {
    Fold(String),
    Select(Def, String),
    /// Select (or deselect) the seed on this edge.
    SelectSeed(fieldcraft_model::GeometryRef),
    /// Select the item a right click opened its menu on (nothing for an item that cannot be).
    Focus(Target),
    ShowResult(String),
    /// A ghost row's command.
    Start(&'static str),
    /// An item menu's entry, or an item's eye.
    Menu(Target, Action),
}

/// A definition's icon in the tree, by what it is: thermal boundary conditions and loads carry
/// the heat icon (a temperature is no support), and a boundary condition that moves the part
/// (a non-zero prescribed displacement) its kind's icon with a motion arrow (`true`).
fn item_icon(m: &Model, kind: Def, name: &str) -> (&'static str, bool) {
    const HEAT: &str = "model.step.heat";
    match kind {
        Def::BoundaryCondition => match m.boundary_condition(name) {
            Some(bc) if bc.is_thermal() => (HEAT, false),
            Some(bc) => (kind.icon(), bc.drives()),
            None => (kind.icon(), false),
        },
        Def::Load => match m.loads().get(name) {
            Some(l) if l.is_thermal() || matches!(l, fieldcraft_model::Load::Temperature { .. }) => (HEAT, false),
            _ => (kind.icon(), false),
        },
        _ => (kind.icon(), false),
    }
}

/// A seed in a few words: `8 exact`, `size 2.5 min, bias 4 end`, `face, 8 exact`.
fn seed_note(s: &fieldcraft_model::Seed) -> String {
    let what = match s.divisions {
        fieldcraft_model::Divisions::Count(n) => format!("{n}"),
        fieldcraft_model::Divisions::Size(h) => format!("size {h}"),
    };
    let bias = s.bias.map(|b| format!(", bias {} {}", b.ratio, json_name(&b.toward))).unwrap_or_default();
    let face = if s.on_face() { "face, " } else { "" };
    format!("{face}{what} {}{bias}", s.constraint.name())
}

fn json_name<T: serde::Serialize>(v: &T) -> String {
    serde_json::to_value(v).ok().and_then(|v| v.as_str().map(str::to_string)).unwrap_or_default()
}

pub fn browser(app: &mut ModelApp, ui: &mut egui::Ui) {
    let t = Tokens::get(app.ui.dark);
    let mut acts: Vec<Act> = Vec::new();
    egui::Panel::left("fc_browser")
        .exact_size(WIDTH)
        .resizable(false)
        .show_separator_line(true)
        .frame(egui::Frame::NONE.fill(t.panel).inner_margin(egui::Margin::symmetric(6, 6)))
        .show(ui, |ui| {
            app.place("browser", ui.max_rect());
            // The title band as the kit's tree cards draw theirs (`solidcraft_ui_kit::tree_header`).
            let k = solidcraft_theme::Tokens::get();
            let (band, _) = ui.allocate_exact_size(vec2(ui.available_width(), k.form_header_height), Sense::hover());
            let title = Rect::from_min_max(band.min + vec2(k.space_12, 0.0), band.max);
            kit::paint_text(ui, title, "Browser", kit::type_face::semibold(k.font_heading), k.text, false);
            k.secondary_hline(ui.painter(), band.x_range(), band.bottom(), k.dialog_rule);
            ui.add_space(4.0);
            egui::ScrollArea::vertical().id_salt("fc_browser_scroll").auto_shrink([false, false]).show(ui, |ui| {
                ui.spacing_mut().item_spacing.y = 0.0;
                rows(app, ui, &mut acts);
            });
        });
    for a in acts {
        match a {
            Act::Fold(key) => {
                if !app.tree_open.remove(&key) {
                    app.tree_open.insert(key);
                }
            }
            Act::Select(kind, name) => {
                let again = app.selected.as_ref() == Some(&(kind, name.clone()));
                app.select_definition((!again).then_some((kind, name)));
            }
            Act::SelectSeed(edge) => {
                if app.selected_seed.as_ref() == Some(&edge) {
                    app.selected_seed = None;
                } else {
                    app.select_seed(edge);
                }
            }
            Act::Focus(target) => match target {
                Target::Def(kind, name) => app.select_definition(Some((kind, name))),
                Target::Seed(edge) => app.select_seed(edge),
                Target::Body(_) | Target::Surface(_) | Target::View(_) | Target::Render(_) | Target::RigidWall { .. } => app.select_definition(None),
            },
            Act::Menu(target, action) => menu::apply(app, target, action),
            Act::ShowResult(step) => {
                // Errors show in the status bar.
                let _ = app.run("results.deformed", json!({"step": step}));
            }
            Act::Start(command) => {
                // A command that cannot run yet says why in the status bar.
                let _ = app.start(command);
            }
        }
    }
}

/// One item row of the tree: what it is (its menu's target) and how it reads.
#[derive(Clone, Debug, PartialEq)]
struct Item {
    target: Target,
    icon: &'static str,
    label: String,
    /// Dim text after the label.
    note: Option<String>,
    /// Shown dim (suppressed).
    dim: bool,
    /// A boundary condition that moves the part ([`Row::moves`]).
    moves: bool,
}

impl Item {
    fn new(target: Target, icon: &'static str, label: impl Into<String>, note: Option<String>) -> Item {
        Item { target, icon, label: label.into(), note, dim: false, moves: false }
    }

    /// A definition of the model, with its kind's icon.
    fn def(m: &Model, kind: Def, name: &str, note: Option<String>) -> Item {
        let (icon, moves) = item_icon(m, kind, name);
        Item { moves, ..Item::new(Target::Def(kind, name.to_string()), icon, name, note) }
    }

    fn draw(self, app: &mut ModelApp, ui: &mut egui::Ui, depth: usize, acts: &mut Vec<Act>) {
        let id = egui::Id::unique(("fc_row_item", &self.target));
        let row = Row { depth, icon: Some(self.icon), moves: self.moves, label: &self.label, note: self.note, dim: self.dim, ..Default::default() };
        item_row(app, ui, id, row, self.target, acts);
    }
}

/// What Geometry lists under the design: datum planes, axes and points, the partitions in the
/// order they cut, then what is kept for meshing only: face merges and absorbs, face splits,
/// stitches and defeatures.
fn geometry_items(m: &Model) -> Vec<Item> {
    let mut out = Vec::new();
    for (n, d) in m.datums() {
        let kind = d.kind();
        let icon = match kind {
            "axis" => "geometry.datum.axis",
            "point" => "geometry.datum.point",
            _ => "geometry.datum",
        };
        let off = m.datum_suppressed(n);
        let note = if off { format!("{kind} · suppressed") } else { kind.into() };
        out.push(Item { dim: off, ..Item::new(Target::Def(Def::Datum, n.clone()), icon, n, Some(note)) });
    }
    for p in m.partitions() {
        let off = !m.partition_cuts(p);
        let note = if off { format!("by {} · suppressed", p.datum) } else { format!("by {}", p.datum) };
        out.push(Item { dim: off, ..Item::new(Target::Def(Def::Partition, p.name.clone()), "geometry.partition", &p.name, Some(note)) });
    }
    for f in m.face_merges() {
        let note = match &f.into {
            Some(into) => format!("{} · absorbed into {into}", f.body),
            None => format!("{} · {} faces", f.body, f.faces.len()),
        };
        out.push(Item::def(m, Def::FaceMerge, &f.name, Some(note)));
    }
    for f in m.face_splits() {
        out.push(Item::def(m, Def::FaceSplit, &f.name, Some(format!("{} · by {}", f.body, f.datum))));
    }
    for f in m.stitches() {
        out.push(Item::def(m, Def::Stitch, &f.name, Some(format!("on {}", f.surface))));
    }
    for f in m.defeatures() {
        let n = f.features.len();
        out.push(Item::def(m, Def::Defeature, &f.name, Some(format!("{n} feature{}", if n == 1 { "" } else { "s" }))));
    }
    out
}

/// What Mesh lists after the seeds: the local size controls, by where they apply, the saved
/// mesh plans, by body, and the relaxation fields.
fn mesh_items(m: &Model) -> Vec<Item> {
    use fieldcraft_model::{Entity, SizeRegion};
    let mut out: Vec<Item> = m
        .size_controls()
        .iter()
        .map(|(n, c)| {
            let on = match &c.region {
                SizeRegion::Geometry(refs) if refs.iter().all(|r| matches!(r.entity, Entity::Body)) => "bodies",
                SizeRegion::Geometry(_) => "faces",
                SizeRegion::Box { .. } => "box",
                SizeRegion::Sphere { .. } => "sphere",
            };
            Item::def(m, Def::SizeControl, n, Some(on.into()))
        })
        .collect();
    out.extend(m.mesh_plans().iter().map(|(n, p)| Item::def(m, Def::MeshPlan, n, Some(format!("for {}", p.body)))));
    out.extend(m.relaxations().iter().map(|(n, r)| {
        let note = if r.bodies.is_empty() { None } else { Some(format!("on {}", r.bodies.join(", "))) };
        Item::def(m, Def::Relaxation, n, note)
    }));
    out
}

/// The definition groups below Mesh: key, label and the kinds each lists.
const GROUPS: [(&str, &str, &[Def]); 7] = [
    ("sets", "Sets", &[Def::Set]),
    ("materials", "Materials", &[Def::Material]),
    ("sections", "Sections", &[Def::Section]),
    ("steps", "Steps", &[Def::Step]),
    ("loads", "Loads & BCs", &[Def::BoundaryCondition, Def::Load, Def::InitialField]),
    ("interactions", "Interactions", &[Def::ReferencePoint, Def::Coupling, Def::Contact, Def::Tie, Def::Spring, Def::Mass]),
    ("fluids", "Fluid Domains", &[Def::FluidDomain]),
];

/// What a definition group lists: its kinds' definitions, and under Interactions the rigid
/// walls of the explicit steps (`Rigid wall on front (Impact)`), each its own item.
fn group_items(m: &Model, key: &str, kinds: &[Def]) -> Vec<Item> {
    let mut out: Vec<Item> = Vec::new();
    for kind in kinds {
        for name in kind.names(m) {
            let note = match kind {
                Def::Spring => m.spring(&name).map(|s| s.kind().to_string()),
                Def::InitialField => m.initial_field(&name).map(|f| format!("initial {}", f.kind())),
                _ => None,
            };
            out.push(Item::def(m, *kind, &name, note));
        }
    }
    if key == "interactions" {
        for s in m.steps() {
            if let fieldcraft_model::Procedure::ExplicitDynamic(e) = &s.procedure {
                for (index, w) in e.rigid_walls.iter().enumerate() {
                    let target = Target::RigidWall { step: s.name.clone(), index };
                    out.push(Item::new(target, "interaction.rigid_wall", format!("Rigid wall on {} ({})", w.set, s.name), None));
                }
            }
        }
    }
    out
}

/// What Results lists after the solved steps: the saved views and the renders made from them.
fn saved_items(app: &ModelApp) -> Vec<Item> {
    let view = &app.session.view;
    let views = view.views.iter().map(|v| {
        let note = (view.current.as_deref() == Some(v.name.as_str())).then(|| "shown".to_string());
        Item::new(Target::View(v.name.clone()), "ui.views", &v.name, note)
    });
    let renders =
        view.renders.iter().map(|r| Item::new(Target::Render(r.name.clone()), "results.export.vtu", &r.name, Some(format!("of {}", r.view))));
    views.chain(renders).collect()
}

fn eye_flag<'u>(state: &'u mut crate::UiState, group: &str) -> Option<&'u mut bool> {
    match group {
        "loads" => Some(&mut state.show_loads),
        "interactions" => Some(&mut state.show_interactions),
        _ => None,
    }
}

/// The ghost rows of `node`'s stages still to do.
fn ghost_rows(app: &mut ModelApp, ui: &mut egui::Ui, stages: &[StageState], node: &str, acts: &mut Vec<Act>) {
    for s in ghosts(stages, node) {
        if draw_ghost(app, ui, s) {
            acts.push(Act::Start(s.command));
        }
    }
}

fn rows(app: &mut ModelApp, ui: &mut egui::Ui, acts: &mut Vec<Act>) {
    let t = Tokens::get(app.ui.dark);
    let stages = app.setup_stages();
    let state = |node: &str| node_state(&stages, node).map_or((None, None), |(st, tip)| (Some(st), Some(tip)));
    let open = |app: &ModelApp, key: &str| app.tree_open.contains(key);
    let name = app.model_name();
    let dirty = app.session.is_dirty();
    draw_row(
        app,
        ui,
        egui::Id::unique("fc_row_root"),
        &Row { icon: Some("ui.views"), label: &name, note: dirty.then(|| "unsaved".into()), ..Default::default() },
    );
    // Geometry.
    let design = app.session.design.clone();
    let surfaces: Vec<String> = app.session.surfaces.iter().map(|s| s.name.clone()).collect();
    let bodies = design.as_ref().map_or(0, |d| d.bodies.len());
    let (mark, tip) = state("geometry");
    let r = draw_row(
        app,
        ui,
        egui::Id::unique("fc_row_geometry"),
        &Row {
            depth: 1,
            fold: Some(open(app, "geometry")),
            eye: Some(app.ui.show_geometry),
            icon: Some("geometry.open"),
            label: "Geometry",
            note: Some(format!("{}", bodies + surfaces.len())),
            mark,
            tip,
            ..Default::default()
        },
    );
    if r.eye {
        app.ui.show_geometry = !app.ui.show_geometry;
    } else if r.fold || r.clicked {
        acts.push(Act::Fold("geometry".into()));
    }
    if open(app, "geometry") {
        if let Some(d) = &design {
            let label = d.name.clone();
            draw_row(
                app,
                ui,
                egui::Id::unique("fc_row_design"),
                &Row { depth: 2, icon: Some("geometry.open"), label: &label, note: Some(format!("{bodies} bodies")), ..Default::default() },
            );
            for (i, b) in d.bodies.iter().enumerate() {
                let row = Row { depth: 3, icon: Some("geometry.bodies"), label: &b.name, ..Default::default() };
                item_row(app, ui, egui::Id::unique(("fc_row_body", i)), row, Target::Body(b.name.clone()), acts);
            }
        }
        for (i, s) in surfaces.iter().enumerate() {
            let row = Row { depth: 2, icon: Some("geometry.import.mesh"), label: s, ..Default::default() };
            item_row(app, ui, egui::Id::unique(("fc_row_surface", i)), row, Target::Surface(s.clone()), acts);
        }
        for item in geometry_items(&app.session.model) {
            item.draw(app, ui, 2, acts);
        }
    }
    ghost_rows(app, ui, &stages, "geometry", acts);
    // Mesh.
    let (nodes, elements) = (app.session.model.nodes().len(), app.session.model.elements().len());
    let (mark, tip) = state("mesh");
    let r = draw_row(
        app,
        ui,
        egui::Id::unique("fc_row_mesh"),
        &Row {
            depth: 1,
            fold: Some(open(app, "mesh")),
            eye: Some(app.ui.show_mesh),
            icon: Some("mesh.tet.generate"),
            label: "Mesh",
            note: Some(format!("{elements}")),
            mark,
            tip,
            ..Default::default()
        },
    );
    if r.eye {
        app.ui.show_mesh = !app.ui.show_mesh;
    } else if r.fold {
        acts.push(Act::Fold("mesh".into()));
    } else if r.clicked {
        // A click shows how the mesh was made (#343) and opens the node; again, deselects it.
        if app.selected_mesh {
            app.selected_mesh = false;
        } else {
            app.select_mesh();
        }
        if !open(app, "mesh") {
            acts.push(Act::Fold("mesh".into()));
        }
    }
    if open(app, "mesh") {
        let r = draw_row(
            app,
            ui,
            egui::Id::unique("fc_row_nodes"),
            &Row {
                depth: 2,
                eye: Some(app.ui.show_nodes),
                icon: Some("mesh.node.add"),
                label: "Nodes",
                note: Some(format!("{nodes}")),
                ..Default::default()
            },
        );
        if r.eye {
            app.ui.show_nodes = !app.ui.show_nodes;
        }
        let revision = app.session.revision();
        if app.element_counts.as_ref().is_none_or(|(r, _)| *r != revision) {
            app.element_counts = Some((revision, app.session.model.element_counts()));
        }
        for (kind, n) in app.element_counts.as_ref().map(|(_, c)| c.clone()).unwrap_or_default() {
            draw_row(
                app,
                ui,
                egui::Id::unique(("fc_row_kind", kind.name())),
                &Row { depth: 2, swatch: Some(t.element(kind)), label: kind.name(), note: Some(format!("{n}")), ..Default::default() },
            );
        }
        // Seeds on edges and faces: a click selects the seed (its properties on the right) and its
        // edge or face in the viewport; a double click or the menu edits or deletes it.
        let seeds: Vec<(fieldcraft_model::GeometryRef, String)> = app.session.model.seeds().iter().map(|s| (s.on.clone(), seed_note(s))).collect();
        if !seeds.is_empty() {
            let r = draw_row(
                app,
                ui,
                egui::Id::unique("fc_row_seeds"),
                &Row {
                    depth: 2,
                    fold: Some(open(app, "seeds")),
                    icon: Some("mesh.seed.edge"),
                    label: "Seeds",
                    note: Some(format!("{}", seeds.len())),
                    ..Default::default()
                },
            );
            if r.fold || r.clicked {
                acts.push(Act::Fold("seeds".into()));
            }
            if open(app, "seeds") {
                for (i, (edge, note)) in seeds.into_iter().enumerate() {
                    let label = format!("{} · {}", edge.body, edge.entity.name().unwrap_or_default());
                    let row = Row { depth: 3, icon: Some("mesh.seed.edge"), label: &label, note: Some(note), ..Default::default() };
                    item_row(app, ui, egui::Id::unique(("fc_row_seed", i)), row, Target::Seed(edge), acts);
                }
            }
        }
        for item in mesh_items(&app.session.model) {
            item.draw(app, ui, 2, acts);
        }
    }
    ghost_rows(app, ui, &stages, "mesh", acts);
    // Definitions.
    for (key, label, kinds) in GROUPS {
        let items = group_items(&app.session.model, key, kinds);
        let icon = kinds.first().map(|k| k.icon());
        let eye = eye_flag(&mut app.ui, key).map(|v| *v);
        let (mark, tip) = state(key);
        let r = draw_row(
            app,
            ui,
            egui::Id::unique(("fc_row_group", key)),
            &Row {
                depth: 1,
                fold: Some(open(app, key)),
                eye,
                icon,
                label,
                note: Some(format!("{}", items.len())),
                dim: items.is_empty(),
                mark,
                tip,
                ..Default::default()
            },
        );
        if r.eye
            && let Some(v) = eye_flag(&mut app.ui, key)
        {
            *v = !*v;
        } else if r.fold || r.clicked {
            acts.push(Act::Fold(key.into()));
        }
        if open(app, key) {
            for item in items {
                item.draw(app, ui, 2, acts);
            }
        }
        ghost_rows(app, ui, &stages, key, acts);
    }
    // Results.
    let steps: Vec<(String, bool)> = app.session.result_steps().map(|(n, c)| (n.to_string(), c)).collect();
    let shown = app.session.results_view.as_ref().map(|v| v.step.clone());
    let (mark, tip) = state("results");
    let r = draw_row(
        app,
        ui,
        egui::Id::unique("fc_row_results"),
        &Row {
            depth: 1,
            fold: Some(open(app, "results")),
            icon: Some("results.contour"),
            label: "Results",
            note: Some(format!("{}", steps.len())),
            dim: steps.is_empty(),
            mark,
            tip,
            ..Default::default()
        },
    );
    if r.fold || r.clicked {
        acts.push(Act::Fold("results".into()));
    }
    if open(app, "results") {
        for (step, current) in steps {
            let label = if current { step.clone() } else { format!("{step} (stale)") };
            let selected = shown.as_deref() == Some(step.as_str());
            let r = draw_row(
                app,
                ui,
                egui::Id::unique(("fc_row_result", &step)),
                &Row { depth: 2, icon: Some("results.deformed"), label: &label, selected, dim: !current, ..Default::default() },
            );
            if r.clicked {
                acts.push(Act::ShowResult(step));
            }
        }
        for item in saved_items(app) {
            item.draw(app, ui, 2, acts);
        }
    }
    ghost_rows(app, ui, &stages, "results", acts);
}

#[cfg(test)]
pub(crate) mod tests {
    use fieldcraft_engine::Session;

    use super::*;
    use Status::*;

    /// The impact sample (an explicit step with a rigid wall) holding one of every other kind
    /// the browser lists: a point mass, a spring, a size control, a face merge and an absorb, a
    /// face split, a stitch, a defeature, a saved view and a render of it.
    pub(crate) fn every_kind() -> ModelApp {
        let mut app = ModelApp::new(Session::default());
        app.run_script(&fieldcraft_engine::sample::impact()).unwrap();
        app.run("model.mass.create", json!({"name": "motor", "target": {"set": "tip"}, "mass": 2.5})).unwrap();
        app.run(
            "model.spring.create",
            json!({"name": "mount", "type": "spring", "stiffness": 1000, "first": {"set": "tip"}, "direction": {"dof": 1}}),
        )
        .unwrap();
        app.run("mesh.size.set", json!({"name": "Fine", "size": 1, "box": {"min": [0, 0, 0], "max": [10, 10, 10]}})).unwrap();
        app.run("geometry.datum", json!({"name": "Mid", "by": {"type": "offset", "plane": "yz", "offset": 15}})).unwrap();
        // Merges, splits, stitches and defeatures need a design or a surface to be made by their
        // commands; the model keeps them by name, as a model file read back would.
        let m = &mut app.session.model;
        let merge = |name: &str, faces: [&str; 2], into: Option<&str>| fieldcraft_model::FaceMerge {
            name: name.into(),
            body: "Base".into(),
            faces: faces.map(str::to_string).to_vec(),
            into: into.map(str::to_string),
        };
        m.add_face_merge(merge("Fillet", ["F3:blend#0", "F3:blend#1"], None)).unwrap();
        m.add_face_merge(merge("Sliver", ["F4:side#2", "F1:+z"], Some("F1:+z"))).unwrap();
        m.add_face_split(fieldcraft_model::FaceSplit { name: "Band".into(), body: "Base".into(), face: "F1:+z".into(), datum: "Mid".into() })
            .unwrap();
        m.add_stitch(fieldcraft_model::Stitch { name: "Gaps".into(), surface: "Mesh1".into(), tolerance: 0.05 }).unwrap();
        m.set_defeature(fieldcraft_model::Defeature { name: "Small".into(), features: vec![3, 7] }).unwrap();
        let warm = fieldcraft_model::InitialField::Temperature { values: [(1, 300.0)].into(), source: "moulding.fieldcraft".into() };
        m.set_initial_field("Warm", warm).unwrap();
        m.set_relaxation("Relax1", fieldcraft_model::Relaxation { bodies: vec!["Base".into()], cells: Vec::new(), reapply: false }).unwrap();
        let plan = json!({"body": "Base", "size": 2.0, "layers": 2, "mesher": "sweep"});
        m.set_mesh_plan("Plan-Base", serde_json::from_value(plan).unwrap()).unwrap();
        app.run("view.save", json!({"name": "Front", "camera": {"yaw": 0, "pitch": 0}, "results": false})).unwrap();
        app.run("results.render.save", json!({"name": "Clip", "view": "Front", "step": "Impact"})).unwrap();
        app
    }

    /// Every item the tree lists, under each node.
    pub(crate) fn listed(app: &ModelApp) -> Vec<Target> {
        let m = &app.session.model;
        let mut out: Vec<Item> = geometry_items(m);
        out.extend(mesh_items(m));
        for (key, _, kinds) in GROUPS {
            out.extend(group_items(m, key, kinds));
        }
        out.extend(saved_items(app));
        out.into_iter().map(|i| i.target).collect()
    }

    /// Every kind of definition the model keeps has its place in the tree (#A1-08): each is
    /// listed once, so it can be selected, edited and deleted from its row; rigid walls are
    /// items of their own, and saved views and renders are listed under Results.
    #[test]
    fn every_kind_the_model_keeps_is_listed() {
        let app = every_kind();
        let listed = listed(&app);
        for kind in Def::ALL {
            let names = kind.names(&app.session.model);
            let rows: Vec<&String> = listed
                .iter()
                .filter_map(|t| match t {
                    Target::Def(k, n) if *k == kind => Some(n),
                    _ => None,
                })
                .collect();
            assert_eq!(rows, names.iter().collect::<Vec<_>>(), "{kind:?}");
        }
        for kind in [
            Def::Mass,
            Def::Spring,
            Def::SizeControl,
            Def::FaceMerge,
            Def::FaceSplit,
            Def::Stitch,
            Def::Defeature,
            Def::Datum,
            Def::InitialField,
            Def::MeshPlan,
            Def::Relaxation,
        ] {
            assert!(!kind.names(&app.session.model).is_empty(), "the fixture has a {kind:?}");
        }
        assert!(listed.contains(&Target::RigidWall { step: "Impact".into(), index: 0 }), "{listed:?}");
        assert!(listed.contains(&Target::View("Front".into())) && listed.contains(&Target::Render("Clip".into())));
        // Each kind sits in one place.
        let mut placed: Vec<Def> = GROUPS.iter().flat_map(|(_, _, k)| k.iter().copied()).collect();
        placed.extend([Def::Datum, Def::Partition, Def::FaceMerge, Def::FaceSplit, Def::Stitch, Def::Defeature]);
        placed.extend([Def::SizeControl, Def::MeshPlan, Def::Relaxation]);
        placed.sort_by_key(|k| format!("{k:?}"));
        let mut all = Def::ALL.to_vec();
        all.sort_by_key(|k| format!("{k:?}"));
        assert_eq!(placed, all);
        // An absorb reads as one, with the face it went into.
        let geometry = geometry_items(&app.session.model);
        let sliver = geometry.iter().find(|i| i.label == "Sliver").and_then(|i| i.note.clone());
        assert_eq!(sliver.as_deref(), Some("Base · absorbed into F1:+z"));
    }

    /// The new kinds read back as their commands' parameters (in the model's SI): a size
    /// control as mesh.size.set takes it, an absorb as geometry.face.absorb, a split as
    /// geometry.face.split.
    #[test]
    fn meshing_helpers_read_back_as_their_parameters() {
        let app = every_kind();
        let m = &app.session.model;
        let size = definition(m, Def::SizeControl, "Fine").unwrap();
        assert_eq!((size["name"].clone(), size["box"]["max"].clone()), (json!("Fine"), json!([0.01, 0.01, 0.01])), "{size}");
        let sliver = definition(m, Def::FaceMerge, "Sliver").unwrap();
        assert_eq!(sliver, json!({"name": "Sliver", "face": {"body": "Base", "face": "F4:side#2"}, "into": "F1:+z"}));
        let fillet = definition(m, Def::FaceMerge, "Fillet").unwrap();
        assert_eq!(fillet["faces"][1], json!({"body": "Base", "face": "F3:blend#1"}));
        let band = definition(m, Def::FaceSplit, "Band").unwrap();
        assert_eq!(band, json!({"name": "Band", "face": {"body": "Base", "face": "F1:+z"}, "datum": "Mid"}));
        assert_eq!(definition(m, Def::Defeature, "Small").unwrap(), json!({"name": "Small", "features": [3, 7]}));
    }

    /// Each top-level node's state, and the commands of its ghost rows.
    fn marks(app: &mut ModelApp) -> Vec<(&'static str, Option<Status>, Vec<&'static str>)> {
        let stages = app.setup_stages();
        NODES
            .iter()
            .map(|(n, _)| (*n, node_state(&stages, n).map(|(st, _)| st), ghosts(&stages, n).into_iter().map(|s| s.command).collect()))
            .collect()
    }

    fn mark(app: &mut ModelApp, node: &str) -> Option<Status> {
        marks(app).into_iter().find(|(n, ..)| *n == node).and_then(|(_, st, _)| st)
    }

    /// The setup's progress on the tree from an empty model to a solved one: one next stage at a
    /// time, a warning where a command cannot run yet, ghost rows on nodes still lacking theirs,
    /// and no mark on nodes that are no setup stage.
    #[test]
    fn nodes_mark_the_setup_progress() {
        let mut app = ModelApp::new(Session::default());
        let m = marks(&mut app);
        assert_eq!(m[0], ("geometry", Some(Next), vec!["geometry.open"]));
        assert_eq!(m[1], ("mesh", Some(Blocked), vec!["mesh.tet.generate"]));
        assert_eq!(mark(&mut app, "materials"), Some(Todo));
        assert_eq!(mark(&mut app, "results"), Some(Blocked));
        assert_eq!(mark(&mut app, "interactions"), None);
        assert_eq!(mark(&mut app, "fluids"), None);
        let tip = node_state(&app.setup_stages(), "mesh").map(|(_, tip)| tip).unwrap_or_default();
        assert!(tip.starts_with("Needs geometry") && tip.ends_with("— Geometry › Open Design"), "{tip}");
        let next = |app: &mut ModelApp| marks(app).into_iter().filter(|(_, st, _)| *st == Some(Next)).map(|(n, ..)| n).collect::<Vec<_>>();
        assert_eq!(next(&mut app), ["geometry"]);

        // A design: meshing is next.
        let dir = std::env::temp_dir().join(format!("fieldcraft-browser-marks-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let design = dir.join("plate.solidcraft");
        std::fs::write(&design, fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap()).unwrap();
        app.run("geometry.open", json!({"path": design.to_string_lossy()})).unwrap();
        assert_eq!(marks(&mut app)[0], ("geometry", Some(Done), vec![]));
        assert_eq!(next(&mut app), ["mesh"]);

        // Meshed: the material is next; the section waits for it.
        let body = app.session.design.as_ref().unwrap().bodies[0].name.clone();
        app.run("mesh.tet.generate", json!({"body": body, "tolerance": 2.0})).unwrap();
        assert_eq!(mark(&mut app, "mesh"), Some(Done));
        assert_eq!(next(&mut app), ["materials"]);
        assert_eq!(mark(&mut app, "sections"), Some(Blocked));
        let _ = std::fs::remove_dir_all(&dir);

        // Completely defined: solving is next, under Results.
        let mut app = ModelApp::new(Session::default());
        app.run_script(&fieldcraft_engine::sample::cantilever()).unwrap();
        for (node, st, add) in marks(&mut app) {
            match node {
                "results" => assert_eq!((st, add), (Some(Next), vec!["solve.run"])),
                "interactions" | "fluids" => assert_eq!(st, None),
                _ => assert_eq!((node, st, add), (node, Some(Done), vec![])),
            }
        }
        // Solved: every stage done, no ghost rows left.
        app.run("solve.run", json!({})).unwrap();
        assert!(marks(&mut app).iter().all(|(_, st, add)| matches!(st, None | Some(Done)) && add.is_empty()), "{:?}", marks(&mut app));
        let v = describe(&mut app);
        assert_eq!(v[2], json!({"node": "sets", "state": "done", "tip": "Done", "add": []}));
    }

    /// A prescribed displacement and a temperature read apart from a support in the tree (#354,
    /// #355).
    #[test]
    fn driving_and_thermal_conditions_have_their_own_icons() {
        let mut s = Session::default();
        s.run_script(&fieldcraft_engine::sample::rubber_seal()).unwrap();
        assert_eq!(item_icon(&s.model, Def::BoundaryCondition, "Press"), ("model.boundary_condition.create", true));
        assert_eq!(item_icon(&s.model, Def::BoundaryCondition, "Base"), ("model.boundary_condition.create", false));
        let mut s = Session::default();
        s.execute("file.sample", &json!({"name": "heated-bracket"})).unwrap();
        assert_eq!(item_icon(&s.model, Def::BoundaryCondition, "Heater"), ("model.step.heat", false));
        assert_eq!(item_icon(&s.model, Def::Load, "Air").0, "model.step.heat");
    }

    /// Ties are listed under Interactions, read back and edited like the others (#366).
    #[test]
    fn ties_are_interactions() {
        let mut app = ModelApp::new(Session::default());
        app.run_script(&fieldcraft_engine::sample::cantilever()).unwrap();
        let sets = Def::Set.names(&app.session.model);
        app.run("interaction.tie", json!({"name": "glue", "main": sets[0], "secondary": sets[1]})).unwrap();
        assert_eq!(Def::Tie.names(&app.session.model), ["glue"]);
        assert_eq!(Def::Tie.edit(&app.session.model, "glue", app.session.units()).map(|(c, _)| c).as_deref(), Some("interaction.tie.edit"));
        assert!(crate::describe::rows(&app.session.model, Def::Tie, "glue", app.session.units()).is_some_and(|r| !r.is_empty()));
    }

    #[test]
    fn definitions_read_back_as_their_parameters() {
        let mut app = ModelApp::new(Session::default());
        app.run_script(&fieldcraft_engine::sample::cantilever()).unwrap();
        let m = &app.session.model;
        assert_eq!(Def::Step.names(m), ["Load", "Modes"]);
        let steel = definition(m, Def::Material, "Steel").unwrap();
        assert_eq!(steel["name"], "Steel");
        assert!(steel.get("elastic").is_some(), "{steel}");
        assert!(definition(m, Def::Load, "nope").is_none());
        for kind in Def::ALL {
            let no_edit = [Def::Partition, Def::FaceMerge, Def::FaceSplit, Def::Stitch, Def::InitialField, Def::MeshPlan, Def::Relaxation];
            let verbs: &[&str] = if no_edit.contains(&kind) { &["create", "delete"] } else { &["create", "edit", "delete"] };
            for verb in verbs {
                assert!(fieldcraft_engine::find_command(&kind.command(verb)).is_some(), "{kind:?} {verb}");
            }
            assert!(fieldcraft_engine::catalogue::find(kind.icon()).is_some());
        }
    }
}
