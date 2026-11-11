//! What a command needs in the session before it can run: its [`Requirement`]s.
//!
//! Each [`crate::CommandSpec`] lists its requirements (`requires`) next to its id. A requirement
//! is something the command cannot run without — an open design to mesh, a material to give
//! plasticity to, a result to contour — never a workflow order: a script that builds a model step
//! by step meets each command's requirements as it goes.
//!
//! [`Session::execute`] checks them before running a command, so a command whose requirement is
//! missing fails at once, with no partial work, naming each missing piece and the toolbar entry
//! that provides it ([`Unavailable`]). The desktop UI reads the same check
//! ([`Session::availability`]) to grey the ribbon buttons, the command search and Solve, with the
//! same text as their tooltip. Every check is a count or a few comparisons, cheap enough to run
//! for every button on every frame.

use serde::Serialize;

use crate::{CommandSpec, Session, catalogue};

/// Something a command needs in the session.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Requirement {
    /// An open SolidCraft design (`geometry.open`).
    Design,
    /// Geometry to mesh: an open design or an imported surface.
    Geometry,
    /// A datum plane, axis or point.
    Datum,
    /// A datum plane (what a partition cuts by).
    DatumPlane,
    DatumAxis,
    DatumPoint,
    /// A partition of the design.
    Partition,
    /// A face merge.
    FaceMerge,
    /// A stitch.
    Stitch,
    /// A face split.
    FaceSplit,
    /// A defeature.
    Defeature,
    /// An edge seed.
    Seed,
    /// A size control.
    SizeControl,
    /// At least one mesh node.
    Nodes,
    /// At least one element.
    Mesh,
    /// Something to analyse: elements, point masses or springs, or only fluid steps (whose
    /// domains carry their own grid) — the model check's first rule.
    Analysable,
    /// A material for the elements: met when there is a material or no element.
    ElementMaterial,
    /// A section for the elements: met when there is a section or no element.
    ElementSection,
    Set,
    /// A set or a reference point (what boundary conditions, masses and springs act on).
    Target,
    Material,
    Section,
    ReferencePoint,
    Coupling,
    Contact,
    Tie,
    Mass,
    Spring,
    BoundaryCondition,
    Load,
    FluidDomain,
    Step,
    /// A solved step's result (current or stale).
    Result,
    /// Results on the display.
    ResultsShown,
}

impl Requirement {
    /// Whether the session meets it.
    pub fn met(self, s: &Session) -> bool {
        let m = &s.model;
        let has_elements = !m.elements().is_empty();
        match self {
            Requirement::Design => s.design.is_some(),
            Requirement::Geometry => s.design.is_some() || !s.surfaces.is_empty(),
            Requirement::Datum => !m.datums().is_empty(),
            Requirement::DatumPlane => m.datums().values().any(|d| d.plane().is_some()),
            Requirement::DatumAxis => m.datums().values().any(|d| d.axis().is_some()),
            Requirement::DatumPoint => m.datums().values().any(|d| d.point().is_some()),
            Requirement::Partition => !m.partitions().is_empty(),
            Requirement::FaceMerge => !m.face_merges().is_empty(),
            Requirement::Stitch => !m.stitches().is_empty(),
            Requirement::FaceSplit => !m.face_splits().is_empty(),
            Requirement::Defeature => !m.defeatures().is_empty(),
            Requirement::Seed => !m.seeds().is_empty(),
            Requirement::SizeControl => !m.size_controls().is_empty(),
            Requirement::Nodes => !m.nodes().is_empty(),
            Requirement::Mesh => has_elements,
            Requirement::Analysable => {
                has_elements
                    || !m.masses().is_empty()
                    || !m.springs().is_empty()
                    || (!m.steps().is_empty() && m.steps().iter().all(|st| st.procedure.fluid_domain().is_some()))
            }
            Requirement::ElementMaterial => !has_elements || !m.materials().is_empty(),
            Requirement::ElementSection => !has_elements || !m.sections().is_empty(),
            Requirement::Set => !m.sets().is_empty(),
            Requirement::Target => !m.sets().is_empty() || !m.reference_points().is_empty(),
            Requirement::Material => !m.materials().is_empty(),
            Requirement::Section => !m.sections().is_empty(),
            Requirement::ReferencePoint => !m.reference_points().is_empty(),
            Requirement::Coupling => !m.couplings().is_empty(),
            Requirement::Contact => !m.contacts().is_empty(),
            Requirement::Tie => !m.ties().is_empty(),
            Requirement::Mass => !m.masses().is_empty(),
            Requirement::Spring => !m.springs().is_empty(),
            Requirement::BoundaryCondition => !m.boundary_conditions().is_empty(),
            Requirement::Load => !m.loads().is_empty(),
            Requirement::FluidDomain => !m.fluid_domains().is_empty(),
            Requirement::Step => !m.steps().is_empty(),
            Requirement::Result => s.result_steps().next().is_some(),
            Requirement::ResultsShown => s.results_view.is_some(),
        }
    }

    /// What is missing, as a noun phrase ("a material").
    pub fn what(self) -> &'static str {
        match self {
            Requirement::Design => "an open design",
            Requirement::Geometry => "geometry (an open design or an imported surface)",
            Requirement::Datum => "a datum",
            Requirement::DatumPlane => "a datum plane",
            Requirement::DatumAxis => "a datum axis",
            Requirement::DatumPoint => "a datum point",
            Requirement::Partition => "a partition",
            Requirement::FaceMerge => "a face merge",
            Requirement::Stitch => "a stitch",
            Requirement::FaceSplit => "a face split",
            Requirement::Defeature => "a defeature",
            Requirement::Seed => "a seed",
            Requirement::SizeControl => "a local size",
            Requirement::Nodes => "mesh nodes",
            Requirement::Mesh | Requirement::Analysable => "a mesh",
            Requirement::ElementMaterial | Requirement::Material => "a material",
            Requirement::ElementSection | Requirement::Section => "a section",
            Requirement::Set => "a set",
            Requirement::Target => "a set or a reference point",
            Requirement::ReferencePoint => "a reference point",
            Requirement::Coupling => "a coupling",
            Requirement::Contact => "a contact",
            Requirement::Tie => "a tie",
            Requirement::Mass => "a point mass",
            Requirement::Spring => "a spring or dashpot",
            Requirement::BoundaryCondition => "a boundary condition",
            Requirement::Load => "a load",
            Requirement::FluidDomain => "a fluid domain",
            Requirement::Step => "a step",
            Requirement::Result => "a solved step",
            Requirement::ResultsShown => "results on the display",
        }
    }

    /// The command that provides it.
    pub fn provider(self) -> &'static str {
        match self {
            Requirement::Design | Requirement::Geometry => "geometry.open",
            Requirement::Datum | Requirement::DatumPlane => "geometry.datum",
            Requirement::DatumAxis => "geometry.datum.axis",
            Requirement::DatumPoint => "geometry.datum.point",
            Requirement::Partition => "geometry.partition",
            Requirement::FaceMerge => "geometry.faces.merge",
            Requirement::Stitch => "geometry.stitch",
            Requirement::FaceSplit => "geometry.face.split",
            Requirement::Defeature => "geometry.defeature",
            Requirement::Seed => "mesh.seed.edge",
            Requirement::SizeControl => "mesh.size.set",
            Requirement::Nodes => "mesh.node.add",
            Requirement::Mesh | Requirement::Analysable => "mesh.tet.generate",
            Requirement::ElementMaterial | Requirement::Material => "model.material.create",
            Requirement::ElementSection | Requirement::Section => "model.section.create",
            Requirement::Set | Requirement::Target => "model.set.create",
            Requirement::ReferencePoint => "model.reference_point.create",
            Requirement::Coupling => "model.coupling.create",
            Requirement::Contact => "interaction.contact.pair",
            Requirement::Tie => "interaction.tie",
            Requirement::Mass => "model.mass.create",
            Requirement::Spring => "model.spring.create",
            Requirement::BoundaryCondition => "model.boundary_condition.create",
            Requirement::Load => "model.load.create",
            Requirement::FluidDomain => "model.fluid_domain.create",
            Requirement::Step => "model.step.create",
            Requirement::Result => "solve.run",
            Requirement::ResultsShown => "results.contour",
        }
    }

    /// Where the provider is on the toolbar: "Properties › Material".
    pub fn place(self) -> String {
        let id = self.provider();
        match catalogue::entry_of(id) {
            Some(c) => format!("{} › {}", c.tab.label(), c.label),
            None => id.to_string(),
        }
    }
}

impl crate::schema::Collection {
    /// What an item of the collection is as a requirement: its [`Requirement::provider`] is the
    /// command that makes one, which a form offers when a reference field has nothing to pick.
    pub fn requirement(self) -> Requirement {
        use crate::schema::Collection as C;
        match self {
            C::Set => Requirement::Set,
            C::Material => Requirement::Material,
            C::Section => Requirement::Section,
            C::ReferencePoint => Requirement::ReferencePoint,
            C::Coupling => Requirement::Coupling,
            C::Contact => Requirement::Contact,
            C::Tie => Requirement::Tie,
            C::Mass => Requirement::Mass,
            C::Spring => Requirement::Spring,
            C::BoundaryCondition => Requirement::BoundaryCondition,
            C::Load => Requirement::Load,
            C::FluidDomain => Requirement::FluidDomain,
            C::Step => Requirement::Step,
            C::SolvedStep => Requirement::Result,
            C::Body => Requirement::Design,
            C::Datum => Requirement::Datum,
            C::DatumPlane => Requirement::DatumPlane,
            C::DatumAxis => Requirement::DatumAxis,
            C::DatumPoint => Requirement::DatumPoint,
            C::Partition => Requirement::Partition,
            C::FaceMerge => Requirement::FaceMerge,
            C::Stitch => Requirement::Stitch,
            C::FaceSplit => Requirement::FaceSplit,
            C::SizeControl => Requirement::SizeControl,
            C::Defeature => Requirement::Defeature,
        }
    }
}

/// A requirement the session does not meet.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Missing(pub Requirement);

impl std::fmt::Display for Missing {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Needs {} — {}", self.0.what(), self.0.place())
    }
}

/// Why a command cannot run now: every requirement it misses, in its spec's order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Unavailable {
    pub command: &'static str,
    pub missing: Vec<Missing>,
}

impl std::fmt::Display for Unavailable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} cannot run yet. {}", self.command, self.joined("; "))
    }
}

impl Unavailable {
    fn joined(&self, separator: &str) -> String {
        self.missing.iter().map(ToString::to_string).collect::<Vec<_>>().join(separator)
    }

    /// The missing pieces, one per line.
    pub fn lines(&self) -> String {
        self.joined("\n")
    }

    /// A hover text: `head`, then the missing pieces when there are any.
    pub fn tip(head: &str, missing: Option<&Unavailable>) -> String {
        match missing {
            Some(u) => format!("{head}\n\n{}", u.lines()),
            None => head.to_string(),
        }
    }

    /// As JSON for the control channel and listings: `{"reason", "missing": [{"requirement",
    /// "what", "provider", "place"}]}`.
    pub fn to_json(&self) -> serde_json::Value {
        let missing: Vec<serde_json::Value> = self
            .missing
            .iter()
            .map(|Missing(r)| serde_json::json!({"requirement": r, "what": r.what(), "provider": r.provider(), "place": r.place()}))
            .collect();
        serde_json::json!({"reason": self.lines(), "missing": missing})
    }
}

impl Session {
    /// Whether `spec` can run now: `Err` names every requirement it misses.
    pub fn availability(&self, spec: &CommandSpec) -> std::result::Result<(), Unavailable> {
        let missing: Vec<Missing> = spec.requires.iter().filter(|r| !r.met(self)).map(|&r| Missing(r)).collect();
        if missing.is_empty() { Ok(()) } else { Err(Unavailable { command: spec.id, missing }) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{EngineError, commands};
    use serde_json::json;

    fn missing(s: &Session, id: &str) -> Vec<Requirement> {
        let spec = crate::find_command(id).unwrap();
        s.availability(spec).err().map(|u| u.missing.into_iter().map(|Missing(r)| r).collect()).unwrap_or_default()
    }

    /// Every command in an empty session is either available or names what it misses, with the
    /// toolbar entry that provides it; and running one that misses something fails at once with
    /// that same message, leaving the session as it was.
    #[test]
    fn every_command_in_an_empty_session_runs_or_names_what_it_misses() {
        for spec in commands() {
            let mut s = Session::default();
            match s.availability(spec) {
                Ok(()) => {}
                Err(u) => {
                    assert!(!u.missing.is_empty(), "{}", spec.id);
                    for Missing(r) in &u.missing {
                        assert!(crate::find_command(r.provider()).is_some(), "{}: {r:?}", spec.id);
                        assert!(catalogue::entry_of(r.provider()).is_some(), "{}: {r:?} has no toolbar entry", spec.id);
                        assert!(!r.met(&s));
                    }
                    let journal = s.journal.len();
                    let e = s.execute(spec.id, &json!({})).unwrap_err();
                    assert_eq!(e, EngineError::Unavailable(u.clone()), "{}", spec.id);
                    assert!(e.to_string().contains("Needs "), "{e}");
                    assert_eq!(s.journal.len(), journal);
                }
            }
        }
    }

    /// Each provider meets what it provides on its own: a requirement is never one its provider
    /// itself cannot meet from an empty model (no circular greying).
    #[test]
    fn providers_are_not_greyed_by_what_they_provide() {
        for spec in commands() {
            for r in spec.requires {
                let provider = crate::find_command(r.provider()).unwrap();
                assert!(!provider.requires.contains(r), "{} provides {r:?} but requires it", provider.id);
            }
        }
    }

    /// A collection's requirement is met exactly when it has a name to pick, and its provider is
    /// a command on the toolbar.
    #[test]
    fn collections_require_what_they_list() {
        use crate::schema::Collection as C;
        let all = [
            C::Set,
            C::Material,
            C::Section,
            C::ReferencePoint,
            C::Coupling,
            C::Contact,
            C::Tie,
            C::Mass,
            C::Spring,
            C::BoundaryCondition,
            C::Load,
            C::FluidDomain,
            C::Step,
            C::SolvedStep,
            C::Body,
            C::Datum,
            C::DatumPlane,
            C::DatumAxis,
            C::DatumPoint,
            C::Partition,
            C::FaceMerge,
            C::Stitch,
            C::FaceSplit,
            C::Defeature,
        ];
        let empty = Session::default();
        for c in all {
            let r = c.requirement();
            assert!(empty.names(c).is_empty() && !r.met(&empty), "{c:?}");
            assert!(catalogue::entry_of(r.provider()).is_some(), "{c:?}");
        }
        let mut s = Session::default();
        s.execute("model.material.create", &json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}})).unwrap();
        assert!(C::Material.requirement().met(&s) && s.names(C::Material) == ["Steel"]);
    }

    #[test]
    fn targeted_requirements() {
        let empty = Session::default();
        assert_eq!(missing(&empty, "model.material.plasticity"), [Requirement::Material]);
        assert_eq!(missing(&empty, "results.contour"), [Requirement::Result]);
        assert_eq!(missing(&empty, "solve.static"), [Requirement::Analysable, Requirement::Step]);
        assert_eq!(missing(&empty, "mesh.tet.generate"), [Requirement::Geometry]);
        assert!(missing(&empty, "model.material.create").is_empty());
        assert!(missing(&empty, "model.step.create").is_empty());

        // Meshed, nothing defined: solving needs a material, a section and a step.
        let mut meshed = Session::default();
        let cube = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
        meshed.execute("mesh.block", &json!({"body": "B", "kind": "hex8", "corners": cube, "divisions": [1, 1, 1]})).unwrap();
        assert_eq!(missing(&meshed, "solve.static"), [Requirement::ElementMaterial, Requirement::ElementSection, Requirement::Step]);
        let e = meshed.execute("solve.run", &json!({})).unwrap_err().to_string();
        assert_eq!(
            e,
            "solve.run cannot run yet. Needs a material — Properties › Material; Needs a section — Properties › Section; Needs a step — Steps › Step"
        );

        // A defined, solved cantilever: everything a result needs is there.
        let mut solved = Session::default();
        crate::cmd::solve::tests::cantilever(&mut solved);
        assert!(missing(&solved, "solve.run").is_empty());
        solved.execute("solve.run", &json!({})).unwrap();
        for id in ["results.contour", "results.deformed", "solve.probe", "results.clear", "model.material.plasticity"] {
            assert!(missing(&solved, id).is_empty(), "{id}");
        }
    }
}
