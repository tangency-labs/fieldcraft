//! The capability catalogue: everything FieldCraft does or will do, in one table.
//!
//! Each [`Capability`] is a toolbar entry: its tab and panel, its label and one-line tip, and
//! either the command of the same id it runs ([`Status::Ready`]) or the milestone that brings it
//! ([`Status::Planned`], shown greyed as "Coming soon"). Planned entries are not commands: a
//! command always runs. The icon of an entry is `assets/icons/<id>.svg`.
//!
//! Every registered command is claimed here, as an entry's command or as one of its related
//! commands ([`Capability::more`], the entry's drop-down); the tests check it. A new command
//! claims its place by turning its planned entry into a ready one (or by adding an entry with
//! an icon), see `docs/features.md`.

use std::fmt::Write;

use serde::{Deserialize, Serialize};

/// The toolbar tabs, in the order of an analyst's workflow.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Tab {
    #[default]
    Geometry,
    Mesh,
    Properties,
    Assembly,
    Steps,
    Loads,
    Interactions,
    Solve,
    Results,
    Tools,
}

impl Tab {
    pub const ALL: [Tab; 10] =
        [Tab::Geometry, Tab::Mesh, Tab::Properties, Tab::Assembly, Tab::Steps, Tab::Loads, Tab::Interactions, Tab::Solve, Tab::Results, Tab::Tools];

    /// Its id in JSON (`loads`).
    pub fn id(self) -> &'static str {
        match self {
            Tab::Geometry => "geometry",
            Tab::Mesh => "mesh",
            Tab::Properties => "properties",
            Tab::Assembly => "assembly",
            Tab::Steps => "steps",
            Tab::Loads => "loads",
            Tab::Interactions => "interactions",
            Tab::Solve => "solve",
            Tab::Results => "results",
            Tab::Tools => "tools",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Tab::Geometry => "Geometry",
            Tab::Mesh => "Mesh",
            Tab::Properties => "Properties",
            Tab::Assembly => "Assembly",
            Tab::Steps => "Steps",
            Tab::Loads => "Loads & BCs",
            Tab::Interactions => "Interactions",
            Tab::Solve => "Solve",
            Tab::Results => "Results",
            Tab::Tools => "Tools",
        }
    }
}

/// The workspaces: one per stage of an analysis, each a set of toolbar tabs. Tools is on every
/// workspace.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Workspace {
    /// Geometry and meshing.
    #[default]
    Prepare,
    /// Materials, sections, sets, steps, loads, boundary conditions and interactions.
    Model,
    /// Checking, writing decks and solving.
    Analysis,
    /// Looking at results.
    Results,
}

impl Workspace {
    pub const ALL: [Workspace; 4] = [Workspace::Prepare, Workspace::Model, Workspace::Analysis, Workspace::Results];

    pub fn label(self) -> &'static str {
        match self {
            Workspace::Prepare => "Prepare",
            Workspace::Model => "Model",
            Workspace::Analysis => "Analysis",
            Workspace::Results => "Results",
        }
    }

    /// One line: what it is for.
    pub fn about(self) -> &'static str {
        match self {
            Workspace::Prepare => "Open geometry and mesh it",
            Workspace::Model => "Materials, sections, steps, loads and boundary conditions",
            Workspace::Analysis => "Check the model, write decks and solve",
            Workspace::Results => "Deformed shapes, contours, modes and probes",
        }
    }

    /// Its tabs, in order.
    pub fn tabs(self) -> &'static [Tab] {
        match self {
            Workspace::Prepare => &[Tab::Geometry, Tab::Mesh, Tab::Tools],
            Workspace::Model => &[Tab::Properties, Tab::Assembly, Tab::Steps, Tab::Loads, Tab::Interactions, Tab::Tools],
            Workspace::Analysis => &[Tab::Solve, Tab::Tools],
            Workspace::Results => &[Tab::Results, Tab::Tools],
        }
    }

    /// The workspace that shows `tab`: `current` when it does, else the first that does.
    pub fn showing(tab: Tab, current: Workspace) -> Workspace {
        if current.tabs().contains(&tab) { current } else { Workspace::ALL.into_iter().find(|w| w.tabs().contains(&tab)).unwrap_or(current) }
    }

    /// The catalogue entry whose icon stands for the workspace.
    pub fn icon(self) -> &'static str {
        match self {
            Workspace::Prepare => "mesh.tet.generate",
            Workspace::Model => "model.material.create",
            Workspace::Analysis => "solve.run",
            Workspace::Results => "results.contour",
        }
    }
}

/// Whether an entry runs today.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    /// Runs the command whose id is the entry's.
    Ready,
    /// Coming in this roadmap milestone (`M2`…).
    Planned { milestone: &'static str },
    /// Runs `deck.write` with this `format` once `fieldcraft_decks::FORMATS` has a writer for it;
    /// until then coming in `milestone`.
    Deck { format: &'static str, milestone: &'static str },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Capability {
    /// Our own dotted id; also the icon's file name.
    pub id: &'static str,
    pub label: &'static str,
    pub tab: Tab,
    pub panel: &'static str,
    /// One line: what it does.
    pub tip: &'static str,
    pub status: Status,
    /// Related commands in the entry's drop-down (edit, delete…).
    pub more: &'static [&'static str],
    /// Commands the entry's own form covers (an option of it does the same), or that act on what
    /// it makes from the browser's item menu (edit, delete): claimed by the entry for the command
    /// search and the ribbon's stage, but not listed in the ribbon.
    pub covers: &'static [&'static str],
}

impl Capability {
    /// The command it runs; `None` while it is coming soon.
    pub fn command(&self) -> Option<&'static str> {
        match self.status {
            Status::Ready => Some(self.id),
            Status::Planned { .. } => None,
            Status::Deck { format, .. } => fieldcraft_decks::format(format).map(|_| crate::cmd::deck::WRITE.id),
        }
    }

    /// The parameters its command starts with (an object; empty for most entries).
    pub fn preset(&self) -> serde_json::Value {
        match self.status {
            Status::Deck { format, .. } => serde_json::json!({ "format": format }),
            Status::Ready | Status::Planned { .. } => serde_json::json!({}),
        }
    }

    /// The milestone that brings it, while it is coming soon.
    pub fn milestone(&self) -> Option<&'static str> {
        match self.status {
            Status::Planned { milestone } => Some(milestone),
            Status::Deck { milestone, .. } if self.command().is_none() => Some(milestone),
            Status::Ready | Status::Deck { .. } => None,
        }
    }

    /// The hover text: label, tip and command, or when it comes.
    pub fn tooltip(&self) -> String {
        match (self.command(), self.milestone()) {
            (Some(command), _) => format!("{}\n{}\n{command}", self.label, self.tip),
            (None, m) => format!("{} — coming soon ({})\n{}", self.label, m.unwrap_or("later"), self.tip),
        }
    }
}

const fn ready(id: &'static str, label: &'static str, tab: Tab, panel: &'static str, tip: &'static str) -> Capability {
    Capability { id, label, tab, panel, tip, status: Status::Ready, more: &[], covers: &[] }
}

const fn ready_more(
    id: &'static str,
    label: &'static str,
    tab: Tab,
    panel: &'static str,
    tip: &'static str,
    more: &'static [&'static str],
) -> Capability {
    Capability { more, ..ready(id, label, tab, panel, tip) }
}

const fn covering(
    id: &'static str,
    label: &'static str,
    tab: Tab,
    panel: &'static str,
    tip: &'static str,
    covers: &'static [&'static str],
) -> Capability {
    Capability { covers, ..ready(id, label, tab, panel, tip) }
}

const fn soon(id: &'static str, label: &'static str, tab: Tab, panel: &'static str, milestone: &'static str, tip: &'static str) -> Capability {
    Capability { id, label, tab, panel, tip, status: Status::Planned { milestone }, more: &[], covers: &[] }
}

const fn deck(id: &'static str, label: &'static str, format: &'static str, tip: &'static str) -> Capability {
    Capability { id, label, tab: Tab::Solve, panel: "Decks", tip, status: Status::Deck { format, milestone: "M3" }, more: &[], covers: &[] }
}

use Tab::*;

/// Every capability in toolbar order: grouped by tab (in [`Tab::ALL`] order), then by panel.
pub static CATALOGUE: &[Capability] = &[
    // ---- Geometry ----
    ready("geometry.open", "Open Design", Geometry, "Import", "Open a SolidCraft design: its bodies with persistent face names"),
    ready("geometry.import.mesh", "Tagged Mesh", Geometry, "Import", "Read a face-tagged surface mesh (3MF)"),
    ready_more(
        "geometry.bodies",
        "Bodies",
        Geometry,
        "Inspect",
        "List the open design's bodies, volumes and named faces, edges and vertices",
        &["geometry.centre", "geometry.features"],
    ),
    ready("geometry.tessellate", "Tessellate", Geometry, "Inspect", "Tessellate a body into a face-tagged boundary surface"),
    ready(
        "geometry.check",
        "Check Geometry",
        Geometry,
        "Inspect",
        "Short edges, small and sliver faces and near faces against the element size, located",
    ),
    ready(
        "mesh.surface.check",
        "Check Surface",
        Geometry,
        "Inspect",
        "Weld, closed, manifold, orientation and self-intersection checks; feature edges",
    ),
    // Edit and delete of datums and partitions are on the browser's item menu, not the ribbon.
    covering(
        "geometry.datum",
        "Datum Plane",
        Geometry,
        "Prepare",
        "A datum plane: principal plane and offset, three points or datum points, point and normal, normal to an axis",
        &["geometry.datum.edit", "geometry.datum.delete", "geometry.datum.suppress"],
    ),
    covering(
        "geometry.datum.axis",
        "Datum Axis",
        Geometry,
        "Prepare",
        "A datum axis: principal axis, two points or datum points, point and direction",
        &["geometry.datum.axis.edit"],
    ),
    covering(
        "geometry.datum.point",
        "Datum Point",
        Geometry,
        "Prepare",
        "A datum point: coordinates, a mesh vertex, along an edge, a face's centre, an axis through a plane",
        &["geometry.datum.point.edit"],
    ),
    covering(
        "geometry.partition",
        "Partition",
        Geometry,
        "Prepare",
        "Cut bodies into cells by a datum plane; the two faces of each cut are paired by name",
        &["geometry.partition.delete", "geometry.partition.suppress"],
    ),
    covering(
        "geometry.faces.merge",
        "Merge Faces",
        Geometry,
        "Prepare",
        "Mesh adjacent faces as one (virtual topology): a fillet in pieces, a face split by an imprint; the B-rep is unchanged",
        &["geometry.faces.merge.delete"],
    ),
    ready(
        "geometry.face.absorb",
        "Absorb Faces",
        Geometry,
        "Prepare",
        "Absorb tiny and sliver faces into a neighbour for meshing, one picked or every one Check Geometry found",
    ),
    covering(
        "geometry.face.split",
        "Split Face",
        Geometry,
        "Prepare",
        "Split a face along a datum plane for meshing only: two faces on its surface, the cut an edge between them",
        &["geometry.face.split.delete"],
    ),
    covering(
        "geometry.defeature",
        "Defeature",
        Geometry,
        "Prepare",
        "Take small holes, fillets, chamfers and bosses out for meshing, from a proposal list with each one on or off",
        &["geometry.defeature.propose", "geometry.defeature.delete"],
    ),
    covering(
        "geometry.stitch",
        "Stitch",
        Geometry,
        "Prepare",
        "Close an imported surface's gaps under a tolerance: open points welded, T-junctions split",
        &["geometry.stitch.delete"],
    ),
    soon("geometry.midsurface", "Midsurface", Geometry, "Prepare", "M8", "Extract midsurfaces of thin parts for shell meshes"),
    // ---- Mesh ----
    soon("mesh.size.global", "Global Size", Mesh, "Size", "M2", "One element size for the whole part"),
    ready_more(
        "mesh.seed.edge",
        "Seed Edges",
        Mesh,
        "Size",
        "Element counts or sizes per design edge, with bias and an exact, minimum or free constraint; the meshers honour them",
        &["mesh.seed.clear", "mesh.seed.list"],
    ),
    ready_more(
        "mesh.size.set",
        "Local Size",
        Mesh,
        "Size",
        "Element size on faces or bodies, in a box or in a sphere, growing away at a set rate; every mesher honours it",
        &["mesh.size.growth", "mesh.size.clear", "mesh.size.list"],
    ),
    soon("mesh.size.adapt", "Curvature", Mesh, "Size", "M2", "Refine by curvature and proximity"),
    ready(
        "mesh.tet.generate",
        "Tet Mesh",
        Mesh,
        "Tet",
        "Quality tet4 or tet10 mesh of a body to a size: refined, smoothed and improved, boundary on the B-rep faces, faces tagged",
    ),
    soon("mesh.order.quadratic", "Quadratic", Mesh, "Tet", "M2", "Second-order elements with mid-side nodes on the geometry"),
    ready("mesh.delaunay.points", "Delaunay Points", Mesh, "Tet", "Delaunay tetrahedralization of a point list (debugging aid)"),
    ready(
        "mesh.planar",
        "2D Mesh",
        Mesh,
        "2D",
        "Quad or triangle mesh of a planar face or polygon in the XY plane for plane stress, plane strain and axisymmetric models; edges tagged",
    ),
    soon("mesh.quad.faces", "Quad Faces", Mesh, "Structured Hex", "M5", "All-quad meshes of faces on shared edge nodes"),
    covering(
        "mesh.hex.generate",
        "Hex Mesh",
        Mesh,
        "Structured Hex",
        "Structured hex mesh of a body: automatic, or extruded, swept, revolved or mapped (the Method drop-down)",
        &["mesh.hex.sweep", "mesh.hex.mapped"],
    ),
    covering(
        "mesh.cells",
        "Partitioned Mesh",
        Mesh,
        "Structured Hex",
        "One mesh of a partitioned part, a method per cell (hex, tet, or all hex as conforming blocks); cells joined across partitions, tet cells to hex cells through conforming pyramids",
        &["mesh.hex.blocks"],
    ),
    ready(
        "mesh.hex.extrude",
        "Extrude to Surface",
        Mesh,
        "Structured Hex",
        "Hex layers from a meshed or planar section up to a curved target surface, following the side faces",
    ),
    ready_more(
        "mesh.plan",
        "Mesh Plan",
        Mesh,
        "Structured Hex",
        "A hex-meshing plan from the features that built a body (sweeps, blocks, O-grids, layers, imprints, copies); saved, editable, applied with the structured meshers",
        &["mesh.plan.apply"],
    ),
    ready_more(
        "mesh.block",
        "Mesh Block",
        Mesh,
        "Structured Hex",
        "A mapped block of hex, wedge or quad elements between corner points; projection blocking (a block topology projected onto the part)",
        &["mesh.hex.project"],
    ),
    soon(
        "mesh.hex.decompose",
        "Decompose",
        Mesh,
        "Structured Hex",
        "M5",
        "Find sweepable and mappable regions of 3D parts (extruded cross-sections are already cut into blocks by Hex Mesh)",
    ),
    ready(
        "mesh.hex.octree",
        "Octree Hex",
        Mesh,
        "Automatic Hex",
        "Unstructured all-hex mesh of any part: adaptive grid, boundary layers on every sharp edge, untangled and proven valid",
    ),
    soon("mesh.hex.dominant", "Hex-Dominant", Mesh, "Automatic Hex", "M6", "Hex-dominant mesh by recombining tets"),
    soon("mesh.hex.polycube", "Polycube", Mesh, "Automatic Hex", "M6", "Polycube hex mesh of blocky parts"),
    soon("mesh.smooth", "Smooth", Mesh, "Improve", "M6", "Smooth node positions, keeping every element valid"),
    soon("mesh.untangle", "Untangle", Mesh, "Improve", "M6", "Repair inverted elements"),
    soon("mesh.layers", "Layers", Mesh, "Improve", "M6", "Boundary padding and layers through thin walls"),
    ready_more(
        "mesh.node.relax",
        "Relax Nodes",
        Mesh,
        "Improve",
        "Relax the nodes inside a sphere to fix spikes and slivers; surface nodes slide on the geometry; kept as a field to apply again after a remesh",
        &["mesh.node.relax.apply"],
    ),
    ready(
        "mesh.node.move",
        "Move Nodes",
        Mesh,
        "Improve",
        "Drag nodes of the mesh, the nodes round them following by a falloff or harmonically, on the geometry; the original drawn as a ghost and the elements coloured by how their quality changes",
    ),
    ready_more("mesh.node.add", "Node", Mesh, "Edit", "Add a node", &["mesh.node.remove"]),
    ready_more("mesh.element.add", "Element", Mesh, "Edit", "Add an element on existing nodes", &["mesh.element.remove", "mesh.clear"]),
    ready("mesh.mirror", "Mirror Mesh", Mesh, "Edit", "Mirror a body's mesh across a plane, joined on the plane (mesh half a symmetric part)"),
    ready_more(
        "mesh.node.merge",
        "Merge Nodes",
        Mesh,
        "Edit",
        "Merge coincident nodes, keeping sets, loads and boundary conditions on the survivors; check for duplicates and cracks, renumber",
        &["mesh.check.duplicates", "mesh.renumber"],
    ),
    soon("mesh.region.remesh", "Remesh Region", Mesh, "Edit", "M6", "Remesh a region in place; split and collapse elements"),
    soon("mesh.import.orphan", "Orphan Mesh", Mesh, "Edit", "M6", "Import a mesh without geometry"),
    ready("mesh.quality", "Quality", Mesh, "Verify", "Quality metrics and the certified validity check, failing elements listed"),
    ready(
        "mesh.quality.mark",
        "Mark Bad Elements",
        Mesh,
        "Verify",
        "Mark the elements below quality limits as a set, located and drawn in the marker colour",
    ),
    soon("mesh.verify", "Verify Set", Mesh, "Verify", "M2", "Quality and validity per set"),
    ready("mesh.export.stl", "Export STL", Mesh, "Verify", "Write the mesh's outer surface (or a set's faces) as binary STL"),
    ready_more("mesh.tag", "Tag Geometry", Mesh, "Verify", "Tag elements, faces, edges and nodes with the geometry they mesh", &["mesh.tag.clear"]),
    // ---- Properties ----
    ready_more(
        "model.material.create",
        "Material",
        Properties,
        "Materials",
        "Elastic (isotropic, orthotropic) material with density; neo-Hookean hyperelasticity for nonlinear steps",
        &["model.material.edit", "model.material.delete"],
    ),
    ready("model.material.plasticity", "Plasticity", Properties, "Materials", "Von Mises plasticity: bilinear, hardening table or Johnson–Cook"),
    soon("model.material.hyperelastic", "Hyperelastic", Properties, "Materials", "M8", "Hyperelastic material models"),
    ready("model.material.damage", "Damage", Properties, "Materials", "Johnson–Cook damage and element erosion"),
    soon("model.material.thermal", "Thermal", Properties, "Materials", "M8", "Conductivity, specific heat and expansion"),
    soon(
        "model.material.eos",
        "Equation of State",
        Properties,
        "Materials",
        "M11",
        "Polynomial, ideal gas, Mie–Grüneisen, JWL and tabulated equations of state",
    ),
    ready_more(
        "model.section.create",
        "Section",
        Properties,
        "Sections",
        "Solid, shell, beam, rod or 2D (plane stress, plane strain, axisymmetric) section assigned to a set, with its element formulation",
        &["model.section.edit", "model.section.delete"],
    ),
    soon("model.section.composite", "Composite", Properties, "Sections", "M8", "Composite layups and continuum shells"),
    soon("model.orientation", "Orientation", Properties, "Sections", "M8", "Material orientations"),
    ready_more(
        "model.set.create",
        "Set",
        Properties,
        "Sets",
        "A set of nodes, elements and faces, or of named geometry that survives remeshing",
        &["model.set.edit", "model.set.delete", "model.set.resolve"],
    ),
    ready("model.units", "Units", Properties, "Model", "The units the model is shown and typed in"),
    // ---- Assembly ----
    soon("assembly.instance", "Instances", Assembly, "Instances", "M7", "Part instances placed from SolidCraft occurrences"),
    // ---- Steps ----
    ready_more(
        "model.step.create",
        "Step",
        Steps,
        "Steps",
        "A linear or nonlinear static, frequency, heat or explicit step with its output",
        &["model.step.edit", "model.step.delete"],
    ),
    ready("model.step.activate", "Activate", Steps, "Steps", "Choose which loads and boundary conditions act in a step"),
    ready("model.activity", "Step Table", Steps, "Steps", "The step by load and boundary condition table"),
    soon("model.amplitude", "Amplitude", Steps, "Procedures", "M7", "Time and frequency amplitudes"),
    soon("model.step.nonlinear", "Nonlinear", Steps, "Procedures", "M7", "Nonlinear static steps"),
    soon("model.step.explicit", "Explicit", Steps, "Procedures", "M8", "Explicit dynamic steps"),
    soon("model.step.heat", "Heat Transfer", Steps, "Procedures", "M8", "Heat transfer steps"),
    // ---- Loads & BCs ----
    ready_more(
        "model.boundary_condition.create",
        "Boundary Condition",
        Loads,
        "Boundary Conditions",
        "Fixed, displacement or symmetry condition on a set or reference point",
        &["model.boundary_condition.edit", "model.boundary_condition.delete", "model.boundary_condition.symmetry"],
    ),
    ready_more(
        "model.load.create",
        "Load",
        Loads,
        "Loads",
        "Force, moment, pressure, traction, gravity or body force",
        &["model.load.edit", "model.load.delete"],
    ),
    ready_more(
        "model.fluid_domain.create",
        "Fluid Domain",
        Loads,
        "Flow",
        "The box of fluid cells, its fluid, and the obstacle bodies voxelised from their face-tagged surfaces",
        &["model.fluid_domain.edit", "model.fluid_domain.delete"],
    ),
    ready_more(
        "model.fluid_domain.boundary",
        "Inlet & Outlet",
        Loads,
        "Flow",
        "Velocity inlets, pressure outlets, walls, moving and slip walls, periodic faces; held temperatures, heat fluxes and adiabatic faces; conducting solids' temperature limits",
        &["model.fluid_domain.thermal_boundary", "model.fluid_domain.limit"],
    ),
    // ---- Interactions ----
    ready_more(
        "model.reference_point.create",
        "Reference Point",
        Interactions,
        "Constraints",
        "A reference point for couplings, loads and boundary conditions",
        &["model.reference_point.edit", "model.reference_point.delete"],
    ),
    ready_more(
        "model.coupling.create",
        "Coupling",
        Interactions,
        "Constraints",
        "Kinematic or distributing coupling of a set to a reference point",
        &["model.coupling.edit", "model.coupling.delete"],
    ),
    ready_more(
        "interaction.tie",
        "Tie",
        Interactions,
        "Constraints",
        "Tie non-matching meshes together (dual mortar constraints)",
        &["interaction.tie.edit", "interaction.tie.delete"],
    ),
    ready_more(
        "interaction.contact.general",
        "General Contact",
        Interactions,
        "Contact",
        "Contact between everything that may touch, or self contact of a surface, with friction (explicit dynamics and nonlinear static steps)",
        &["interaction.contact.edit", "interaction.contact.delete"],
    ),
    // Its edit and delete rows are General Contact's: they act on contacts of either kind.
    ready(
        "interaction.contact.pair",
        "Contact Pair",
        Interactions,
        "Contact",
        "Contact between two surfaces: with friction in explicit dynamics, mortar contact with friction in nonlinear static steps",
    ),
    ready(
        "interaction.contact.find",
        "Find Contacts",
        Interactions,
        "Contact",
        "Find the parts of the mesh that touch, and make contact pairs of them",
    ),
    ready_more(
        "interaction.rigid_wall",
        "Rigid Wall",
        Interactions,
        "Contact",
        "A fixed rigid plane a set's nodes cannot pass, with friction, in an explicit dynamics step",
        &["interaction.rigid_wall.delete"],
    ),
    soon("interaction.bolt", "Bolt Load", Interactions, "Connectors", "M7", "Bolt pretension"),
    ready_more(
        "model.spring.create",
        "Spring",
        Interactions,
        "Connectors",
        "Spring or dashpot to the ground or between two points, along a direction or their line",
        &["model.spring.edit", "model.spring.delete"],
    ),
    ready_more(
        "model.mass.create",
        "Point Mass",
        Interactions,
        "Connectors",
        "Point mass and rotary inertia at nodes or a reference point",
        &["model.mass.edit", "model.mass.delete"],
    ),
    soon("interaction.connector", "Connector", Interactions, "Connectors", "M7", "Connector elements between points"),
    soon("interaction.fsi", "FSI Interface", Interactions, "Coupled", "M13", "The fluid–structure interface"),
    // ---- Solve ----
    ready("solve.run", "Solve", Solve, "Run", "Solve every step in order and show the results"),
    ready("model.check", "Check Model", Solve, "Check", "Name every missing or bad definition before writing or solving"),
    deck("deck.write.abaqus", "Abaqus Input", "abaqus", "Write an Abaqus input file (.inp)"),
    deck("deck.write.nastran", "Nastran Deck", "nastran", "Write Nastran bulk data (SOL 101, 103)"),
    deck("deck.write.lsdyna", "LS-DYNA Deck", "lsdyna", "Write an LS-DYNA keyword file"),
    soon("deck.read", "Read Deck", Solve, "Decks", "M9", "Read a deck, keeping keywords FieldCraft does not know"),
    ready_more(
        "solve.static",
        "Static",
        Solve,
        "Implicit",
        "Linear static (or coupled temperature–displacement) solve of a step with FieldCraft's own solver; nonlinear static, buckling and heat transfer in the drop-down",
        &["solve.nonlinear", "solve.buckle", "solve.heat"],
    ),
    ready("solve.frequency", "Frequency", Solve, "Implicit", "Natural frequencies and mode shapes of a frequency step"),
    ready(
        "solve.explicit",
        "Explicit",
        Solve,
        "Explicit",
        "Transient dynamics of an explicit step by central differences (contact and erosion coming in M10)",
    ),
    soon("solve.blast", "Blast & Impact", Solve, "Explicit", "M14", "EOS-driven fluids with fluid–structure contact (ALE or SPH)"),
    ready_more(
        "solve.fluid",
        "Lattice Boltzmann",
        Solve,
        "Fluids",
        "Lattice Boltzmann flow of a fluid step on the CPU or a GPU, with wall forces, and the heat it carries",
        &["solve.settings"],
    ),
    soon("solve.fsi", "Fluid–Structure", Solve, "Fluids", "M13", "Coupled fluid and structure solve"),
    // ---- Results ----
    ready_more("results.deformed", "Deformed", Results, "Display", "The deformed shape of a solved step, to scale", &["results.clear"]),
    ready("results.contour", "Contours", Results, "Display", "Contours of displacement and stress, with a legend"),
    ready("results.modes", "Mode Shapes", Results, "Display", "Mode shapes and frequencies of a frequency step"),
    ready_more(
        "results.animate",
        "Animate",
        Results,
        "Display",
        "Animate the deformed shape, a mode or an explicit step's frames",
        &["results.frame"],
    ),
    ready_more(
        "solve.probe",
        "Probe",
        Results,
        "Query",
        "Displacement, rotation, reaction and stress at nodes, element section forces, a mode shape, or an explicit field frame (or click the model)",
        &["solve.history"],
    ),
    ready_more("results.reactions", "Reactions", Results, "Query", "Reaction force and moment sums per set", &["results.export.csv"]),
    soon("results.free_body", "Free Body", Results, "Query", "M9", "Free-body forces across a cut"),
    ready("results.xy", "XY Plot", Results, "Plots", "Plot a quantity against the load factor, time or mode"),
    ready("results.path", "Path Plot", Results, "Plots", "Interpolate a solid nodal field along a reference-coordinate path, with CSV export"),
    soon("results.cut", "Cut", Results, "Plots", "M9", "Section cuts through the model"),
    soon("results.display_group", "Display Group", Results, "Plots", "M9", "Show and hide parts of the model"),
    ready_more(
        "results.streamlines",
        "Streamlines",
        Results,
        "Plots",
        "A fluid step's velocity, pressure or temperature on slices and obstacles, speed or temperature along streamlines, and signed wall heat per area",
        &["results.slice", "results.surface"],
    ),
    ready_more(
        "results.export.vtu",
        "Export VTU",
        Results,
        "Files",
        "Write a solved step as VTK unstructured grids for ParaView",
        &["results.export.animation", "results.render.save", "results.render.delete"],
    ),
    soon("results.read", "Read Results", Results, "Files", "M9", "Read results in public formats"),
    soon("results.convergence", "Convergence", Results, "Files", "M9", "Mesh convergence studies"),
    // ---- Tools ----
    // Renaming starts from a definition's item menu (Rename…); the search finds it here too.
    ready_more("model.info", "Model Info", Tools, "Model", "Node and element counts and the model summary", &["model.rename"]),
    // New starts from the File menu; the search finds it here too.
    ready_more("file.open", "Open Model", Tools, "File", "Open a .fieldcraft model", &["file.new"]),
    ready("file.sample", "Samples", Tools, "File", "Open a sample model, ready to solve"),
    ready("file.save", "Save", Tools, "File", "Save the model to its .fieldcraft file"),
    ready("file.save_as", "Save As", Tools, "File", "Save the model to a new .fieldcraft file"),
    ready_more(
        "file.recover",
        "Recover",
        Tools,
        "File",
        "Open a model autosaved by an app that crashed or closed with unsaved changes",
        &["file.recovery.list", "file.recovery.discard"],
    ),
    // Hide, show and isolate act on the browser's items, from their eyes and item menus; saved
    // views are saved and shown from the Results panel's Views row.
    covering(
        "view.show_all",
        "Show All",
        Tools,
        "View",
        "Show every item hidden in the viewport (hide, show and isolate items from the browser)",
        &["view.hide", "view.show", "view.isolate", "view.save", "view.restore", "view.delete", "view.list"],
    ),
    ready("edit.undo", "Undo", Tools, "History", "Undo the last command that changed the model"),
    ready("edit.redo", "Redo", Tools, "History", "Redo the command undone last"),
    ready_more(
        "diagnostics.bug_report",
        "Report a Bug",
        Tools,
        "Help",
        "Draft a GitHub bug report with the version, the OS and the last error; you review it, and nothing is sent until you submit it",
        &["diagnostics.crash_reports", "diagnostics.crash_report", "diagnostics.crash_reports.delete", "diagnostics.crash_test", "app.update.check"],
    ),
];

/// The catalogue entry with this id.
pub fn find(id: &str) -> Option<&'static Capability> {
    CATALOGUE.iter().find(|c| c.id == id)
}

/// The entry that claims a command: the one running it, else the one listing it among its
/// related commands or covering it. Every registered command has one (a test checks it).
pub fn entry_of(command: &str) -> Option<&'static Capability> {
    CATALOGUE
        .iter()
        .find(|c| c.command() == Some(command))
        .or_else(|| CATALOGUE.iter().find(|c| c.more.contains(&command) || c.covers.contains(&command)))
}

/// The tab the ribbon follows to when `command` has run from anywhere but the ribbon (a script,
/// the control channel, MCP, the command search, the browser), as a person's click on that tab
/// would: its entry's tab. `None` for commands on the Tools tab (files, history, the model
/// summary: on every workspace and no stage of an analysis) and for unknown commands. A ribbon
/// click runs a command of the tab shown, so it never moves the ribbon (a test checks it).
pub fn stage_tab(command: &str) -> Option<Tab> {
    entry_of(command).map(|c| c.tab).filter(|t| *t != Tab::Tools)
}

/// The entries of a tab (one run of the catalogue).
pub fn tab_entries(tab: Tab) -> &'static [Capability] {
    let start = CATALOGUE.iter().position(|c| c.tab == tab).unwrap_or(CATALOGUE.len());
    let rest = CATALOGUE.get(start..).unwrap_or_default();
    rest.get(..rest.iter().take_while(|c| c.tab == tab).count()).unwrap_or_default()
}

/// The panels of a tab in order, each the run of its entries (never empty).
pub fn panels(tab: Tab) -> impl Iterator<Item = &'static [Capability]> {
    tab_entries(tab).chunk_by(|a, b| a.panel == b.panel)
}

/// `docs/features.md`: every capability by tab and panel, ready or coming soon.
pub fn features_page() -> String {
    let ready = CATALOGUE.iter().filter(|c| c.command().is_some()).count();
    let mut md = String::from(
        "# Features\n\n\
         Generated from the capability catalogue (`crates/engine/src/catalogue.rs`) by\n\
         `cargo run -p fieldcraft-cli -- features > docs/features.md`; a test fails when it is out of date.\n\
         Every entry is a toolbar button. **Ready** entries run their command; **coming soon** entries are\n\
         shown greyed with the roadmap milestone that brings them (ROADMAP.md). Scripts, the command line,\n\
         the control channel and AI agents over MCP (docs/mcp.md) run the same commands.\n\n",
    );
    let _ = writeln!(md, "{ready} ready, {} coming soon.", CATALOGUE.len() - ready);
    for tab in Tab::ALL {
        let _ = write!(md, "\n## {}\n\n| Panel | Capability | Status | Command |\n|---|---|---|---|\n", tab.label());
        for c in tab_entries(tab) {
            let status = c.milestone().map_or_else(|| "ready".to_string(), |m| format!("coming soon ({m})"));
            let preset = c.preset();
            let main =
                c.command().map(|id| if preset.as_object().is_some_and(|o| !o.is_empty()) { format!("`{id} {preset}`") } else { format!("`{id}`") });
            let commands: Vec<String> =
                main.into_iter().chain(c.more.iter().map(|id| format!("`{id}`"))).chain(c.covers.iter().map(|id| format!("(`{id}`)"))).collect();
            let _ = writeln!(md, "| {} | {} — {} | {status} | {} |", c.panel, c.label, c.tip, commands.join(", "));
        }
    }
    md
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{commands, find_command};

    #[test]
    fn every_tab_is_on_a_workspace_and_every_command_has_an_entry() {
        for tab in Tab::ALL {
            assert!(Workspace::ALL.iter().any(|w| w.tabs().contains(&tab)), "{tab:?} is on no workspace");
        }
        assert!(Workspace::ALL.iter().all(|w| w.tabs().contains(&Tab::Tools)));
        assert!(Workspace::ALL.iter().all(|w| find(w.icon()).is_some()));
        assert_eq!(Workspace::showing(Tab::Results, Workspace::Prepare), Workspace::Results);
        assert_eq!(Workspace::showing(Tab::Tools, Workspace::Model), Workspace::Model);
        for c in commands() {
            assert!(entry_of(c.id).is_some(), "{} has no catalogue entry", c.id);
        }
        assert_eq!(entry_of("model.material.edit").map(|c| c.id), Some("model.material.create"));
        for tab in Tab::ALL {
            assert_eq!(serde_json::to_value(tab).ok().and_then(|v| v.as_str().map(str::to_string)).as_deref(), Some(tab.id()));
        }
    }

    /// Each command has one place in the ribbon: no two menu rows (an entry or a drop-down row)
    /// run the same command with the same preset, and a command an entry's form covers is no row.
    #[test]
    fn no_command_appears_twice_across_ribbon_menus() {
        let mut rows: Vec<String> = Vec::new();
        for c in CATALOGUE {
            let entry = c.command().map(|id| format!("{id} {}", c.preset()));
            for row in entry.into_iter().chain(c.more.iter().map(|id| format!("{id} {{}}"))) {
                assert!(!rows.contains(&row), "{}: `{row}` is in the ribbon twice", c.id);
                rows.push(row);
            }
        }
        for c in CATALOGUE {
            for m in c.covers {
                assert!(!rows.iter().any(|r| r.split(' ').next() == Some(m)), "{}: `{m}` is covered by its form and also a ribbon row", c.id);
            }
        }
    }

    #[test]
    fn the_ribbon_follows_each_command_to_its_stage() {
        assert_eq!(stage_tab("mesh.hex.sweep"), Some(Tab::Mesh));
        assert_eq!(stage_tab("model.material.create"), Some(Tab::Properties));
        assert_eq!(stage_tab("model.material.edit"), Some(Tab::Properties));
        assert_eq!(stage_tab("model.load.create"), Some(Tab::Loads));
        assert_eq!(stage_tab("solve.static"), Some(Tab::Solve));
        assert_eq!(stage_tab("results.contour"), Some(Tab::Results));
        assert_eq!(stage_tab("edit.undo"), None);
        assert_eq!(stage_tab("file.save"), None);
        assert_eq!(stage_tab("model.info"), None);
        assert_eq!(stage_tab("no.such"), None);
        // A ribbon click (an entry's command, or one of its drop-down's) stays on the tab shown.
        for c in CATALOGUE {
            for command in c.command().into_iter().chain(c.more.iter().copied()).chain(c.covers.iter().copied()) {
                assert!(stage_tab(command).is_none_or(|t| t == c.tab), "{}: `{command}` would move the ribbon from {:?}", c.id, c.tab);
            }
        }
    }

    #[test]
    fn every_command_is_claimed_and_every_ready_entry_runs_a_command() {
        let mut claimed: Vec<&str> = Vec::new();
        for c in CATALOGUE {
            if let Some(command) = c.command() {
                assert!(find_command(command).is_some(), "{}: `{command}` is not a registered command", c.id);
                claimed.push(command);
            } else {
                assert!(c.more.is_empty() && c.covers.is_empty(), "{}: a planned entry has no commands", c.id);
            }
            for m in c.more.iter().chain(c.covers) {
                assert!(find_command(m).is_some(), "{}: `{m}` is not a registered command", c.id);
                claimed.push(m);
            }
        }
        for spec in commands() {
            let n = claimed.iter().filter(|c| **c == spec.id).count();
            assert!(n >= 1, "command `{}` has no catalogue entry (a ready entry or an entry's `more`)", spec.id);
        }
    }

    /// Deck entries follow the writers: one entry per format, ready exactly when its writer is in
    /// `fieldcraft_decks::FORMATS`, running `deck.write` with its format.
    #[test]
    fn deck_entries_follow_the_writers() {
        let decks: Vec<(&Capability, &str)> =
            CATALOGUE.iter().filter_map(|c| if let Status::Deck { format, .. } = c.status { Some((c, format)) } else { None }).collect();
        for f in fieldcraft_decks::FORMATS {
            assert!(decks.iter().any(|(_, name)| *name == f.name), "deck format `{}` has no catalogue entry", f.name);
        }
        for (c, format) in decks {
            let written = fieldcraft_decks::format(format).is_some();
            assert_eq!(c.command(), written.then_some("deck.write"), "{}", c.id);
            assert_eq!(c.milestone().is_none(), written, "{}", c.id);
            assert_eq!(c.preset(), serde_json::json!({"format": format}));
        }
        assert_eq!(find("deck.write.abaqus").and_then(Capability::command), Some("deck.write"));
    }

    #[test]
    fn entries_are_well_formed() {
        let mut ids: Vec<&str> = CATALOGUE.iter().map(|c| c.id).collect();
        for c in CATALOGUE {
            assert!(c.id.split('.').count() >= 2 && c.id.chars().all(|ch| ch.is_ascii_lowercase() || ch == '.' || ch == '_'), "{}", c.id);
            assert!(!c.label.is_empty() && !c.panel.is_empty() && !c.tip.is_empty(), "{}", c.id);
            if let Some(m) = c.milestone() {
                assert!(m.strip_prefix('M').is_some_and(|n| n.parse::<u32>().is_ok_and(|n| (1..=14).contains(&n))), "{}: {m}", c.id);
            }
        }
        ids.sort_unstable();
        let n = ids.len();
        ids.dedup();
        assert_eq!(ids.len(), n, "catalogue ids are unique");
        // Tabs are runs in Tab::ALL order, and each panel is one run within its tab.
        let total: usize = Tab::ALL.iter().map(|t| tab_entries(*t).len()).sum();
        assert_eq!(total, CATALOGUE.len(), "every tab's entries are together, in Tab::ALL order");
        let order: Vec<Tab> =
            CATALOGUE.iter().map(|c| c.tab).collect::<Vec<_>>().chunk_by(|a, b| a == b).filter_map(|r| r.first().copied()).collect();
        assert_eq!(order, Tab::ALL);
        for tab in Tab::ALL {
            let names: Vec<&str> = panels(tab).filter_map(|p| p.first().map(|c| c.panel)).collect();
            assert!(!names.is_empty(), "{tab:?} has no entries");
            assert!(names.iter().enumerate().all(|(i, n)| !names[..i].contains(n)), "{tab:?}: a panel's entries are together: {names:?}");
        }
    }

    #[test]
    fn features_page_is_current() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/features.md");
        let committed = std::fs::read_to_string(&path).unwrap_or_default().replace("\r\n", "\n");
        assert!(committed == features_page(), "docs/features.md is out of date: run `cargo run -p fieldcraft-cli -- features > docs/features.md`");
    }
}
