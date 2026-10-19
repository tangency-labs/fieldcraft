//! Renaming a definition and every reference to it, so a name can change while other
//! definitions use it (deleting and recreating is refused then).
//!
//! What refers to what: sets are named by sections, couplings, contacts, ties, point masses and
//! springs, boundary conditions, loads and explicit steps (history nodes, initial velocities,
//! rigid walls, mass scaling); materials by sections and fluid domains; reference points by
//! couplings, boundary conditions, loads, point masses and springs; boundary conditions and
//! loads by the steps they are active in; steps by the nonlinear steps whose imperfection is one
//! of their buckling modes; fluid domains by fluid steps; datum planes by the partitions that
//! cut by them. Partitions are not renamed: their cut faces are named after them.

use std::collections::BTreeMap;

use crate::{Error, Model, Procedure, Result, Target, check_name};

/// A kind of definition [`Model::rename`] renames.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Set,
    Material,
    Section,
    Step,
    BoundaryCondition,
    Load,
    ReferencePoint,
    Coupling,
    Contact,
    Tie,
    PointMass,
    Spring,
    FluidDomain,
    Datum,
}

impl Kind {
    /// The kind in words, as errors name it.
    pub fn word(self) -> &'static str {
        match self {
            Kind::Set => "set",
            Kind::Material => "material",
            Kind::Section => "section",
            Kind::Step => "step",
            Kind::BoundaryCondition => "boundary condition",
            Kind::Load => "load",
            Kind::ReferencePoint => "reference point",
            Kind::Coupling => "coupling",
            Kind::Contact => "contact",
            Kind::Tie => "tie",
            Kind::PointMass => "point mass",
            Kind::Spring => "spring",
            Kind::FluidDomain => "fluid domain",
            Kind::Datum => "datum",
        }
    }
}

/// Move `map`'s entry `from` to `to`: `from` must exist and `to` must be free.
fn rekey<T>(map: &mut BTreeMap<String, T>, what: &str, from: &str, to: &str) -> Result<()> {
    if map.contains_key(to) {
        return Err(Error::Invalid(format!("{what} `{to}` already exists")));
    }
    let value = map.remove(from).ok_or_else(|| Error::NotFound(format!("{what} `{from}`")))?;
    map.insert(to.to_string(), value);
    Ok(())
}

/// `name` becomes `to` where it is `from`.
fn swap(name: &mut String, from: &str, to: &str) {
    if name == from {
        *name = to.to_string();
    }
}

fn swap_opt(name: &mut Option<String>, from: &str, to: &str) {
    if let Some(n) = name {
        swap(n, from, to);
    }
}

fn swap_target(target: &mut Target, kind: Kind, from: &str, to: &str) {
    match (target, kind) {
        (Target::Set(s), Kind::Set) | (Target::ReferencePoint(s), Kind::ReferencePoint) => swap(s, from, to),
        _ => {}
    }
}

impl Model {
    /// Rename definition `from` of `kind` to `to`, rewriting every reference to it. A name in use
    /// by another definition of the kind, or an invalid one, is refused and nothing changes;
    /// renaming to the same name does nothing.
    pub fn rename(&mut self, kind: Kind, from: &str, to: &str) -> Result<()> {
        let what = kind.word();
        check_name(what, to)?;
        if from == to {
            return self.exists(kind, from);
        }
        // The definition moves first: it fails before anything else changes.
        match kind {
            Kind::Set => rekey(&mut self.sets, what, from, to)?,
            Kind::Material => rekey(&mut self.materials, what, from, to)?,
            Kind::Section => rekey(&mut self.sections, what, from, to)?,
            Kind::BoundaryCondition => rekey(&mut self.boundary_conditions, what, from, to)?,
            Kind::Load => rekey(&mut self.loads, what, from, to)?,
            Kind::ReferencePoint => rekey(&mut self.reference_points, what, from, to)?,
            Kind::Coupling => rekey(&mut self.couplings, what, from, to)?,
            Kind::Contact => rekey(&mut self.contacts, what, from, to)?,
            Kind::Tie => rekey(&mut self.ties, what, from, to)?,
            Kind::PointMass => rekey(&mut self.masses, what, from, to)?,
            Kind::Spring => rekey(&mut self.springs, what, from, to)?,
            Kind::FluidDomain => rekey(&mut self.fluid_domains, what, from, to)?,
            Kind::Datum => rekey(&mut self.datums, what, from, to)?,
            Kind::Step => {
                if self.steps.iter().any(|s| s.name == to) {
                    return Err(Error::Invalid(format!("step `{to}` already exists")));
                }
                let step = self.steps.iter_mut().find(|s| s.name == from).ok_or_else(|| Error::NotFound(format!("step `{from}`")))?;
                step.name = to.to_string();
            }
        }
        self.rewrite_references(kind, from, to);
        Ok(())
    }

    fn exists(&self, kind: Kind, name: &str) -> Result<()> {
        let found = match kind {
            Kind::Set => self.sets.contains_key(name),
            Kind::Material => self.materials.contains_key(name),
            Kind::Section => self.sections.contains_key(name),
            Kind::BoundaryCondition => self.boundary_conditions.contains_key(name),
            Kind::Load => self.loads.contains_key(name),
            Kind::ReferencePoint => self.reference_points.contains_key(name),
            Kind::Coupling => self.couplings.contains_key(name),
            Kind::Contact => self.contacts.contains_key(name),
            Kind::Tie => self.ties.contains_key(name),
            Kind::PointMass => self.masses.contains_key(name),
            Kind::Spring => self.springs.contains_key(name),
            Kind::FluidDomain => self.fluid_domains.contains_key(name),
            Kind::Datum => self.datums.contains_key(name),
            Kind::Step => self.steps.iter().any(|s| s.name == name),
        };
        if found { Ok(()) } else { Err(Error::NotFound(format!("{} `{name}`", kind.word()))) }
    }

    /// Every reference to `from` of `kind` now names `to`.
    fn rewrite_references(&mut self, kind: Kind, from: &str, to: &str) {
        match kind {
            Kind::Set | Kind::ReferencePoint => {
                let set = kind == Kind::Set;
                for c in self.couplings.values_mut() {
                    swap(if set { &mut c.set } else { &mut c.reference_point }, from, to);
                }
                for b in self.boundary_conditions.values_mut() {
                    use crate::BoundaryCondition as B;
                    let (B::Fixed { target } | B::Displacement { target, .. } | B::Symmetry { target, .. } | B::Temperature { target, .. }) = b;
                    swap_target(target, kind, from, to);
                }
                for l in self.loads.values_mut() {
                    use crate::Load as L;
                    match l {
                        L::Force { target, .. } | L::Moment { target, .. } | L::Temperature { target, .. } => swap_target(target, kind, from, to),
                        L::Traction { set: s, .. }
                        | L::Pressure { set: s, .. }
                        | L::BodyForce { set: s, .. }
                        | L::SurfaceFlux { set: s, .. }
                        | L::Convection { set: s, .. }
                        | L::HeatSource { set: s, .. } => {
                            if set {
                                swap(s, from, to);
                            }
                        }
                        L::Gravity { set: s, .. } | L::Centrifugal { set: s, .. } => {
                            if set {
                                swap_opt(s, from, to);
                            }
                        }
                    }
                }
                for m in self.masses.values_mut() {
                    swap_target(&mut m.target, kind, from, to);
                }
                for s in self.springs.values_mut() {
                    use crate::Spring as S;
                    let (S::Spring { first, second, .. } | S::Dashpot { first, second, .. }) = s;
                    swap_target(first, kind, from, to);
                    if let Some(t) = second {
                        swap_target(t, kind, from, to);
                    }
                }
                if set {
                    self.rewrite_set_only(from, to);
                }
            }
            Kind::Material => {
                for s in self.sections.values_mut() {
                    swap(&mut s.material, from, to);
                }
                for d in self.fluid_domains.values_mut() {
                    swap(&mut d.material, from, to);
                    for c in d.thermal.iter_mut().flat_map(|t| t.solids.iter_mut()) {
                        swap(&mut c.material, from, to);
                    }
                }
            }
            Kind::BoundaryCondition | Kind::Load => {
                for step in &mut self.steps {
                    let active = if kind == Kind::Load { &mut step.loads } else { &mut step.boundary_conditions };
                    if active.remove(from) {
                        active.insert(to.to_string());
                    }
                }
            }
            Kind::Step => {
                for step in &mut self.steps {
                    if let Procedure::NonlinearStatic(n) = &mut step.procedure
                        && let Some(i) = &mut n.imperfection
                    {
                        swap(&mut i.step, from, to);
                    }
                }
            }
            Kind::FluidDomain => {
                for step in &mut self.steps {
                    if let Procedure::Fluid { domain, .. } = &mut step.procedure {
                        swap(domain, from, to);
                    }
                }
            }
            Kind::Datum => {
                for p in &mut self.partitions {
                    swap(&mut p.datum, from, to);
                }
                if self.suppressed_datums.remove(from) {
                    self.suppressed_datums.insert(to.to_string());
                }
            }
            Kind::Section | Kind::Coupling | Kind::Contact | Kind::Tie | Kind::PointMass | Kind::Spring => {}
        }
    }

    /// The references only a set has: sections, contacts, ties and explicit steps.
    fn rewrite_set_only(&mut self, from: &str, to: &str) {
        for s in self.sections.values_mut() {
            swap(&mut s.set, from, to);
        }
        for c in self.contacts.values_mut() {
            match &mut c.surfaces {
                crate::ContactSurfaces::Pair { main, secondary } => {
                    swap(main, from, to);
                    swap(secondary, from, to);
                }
                crate::ContactSurfaces::General { surface } => swap_opt(surface, from, to),
            }
        }
        for t in self.ties.values_mut() {
            swap(&mut t.main, from, to);
            swap(&mut t.secondary, from, to);
        }
        for step in &mut self.steps {
            if let Procedure::ExplicitDynamic(e) = &mut step.procedure {
                swap_opt(&mut e.history_nodes, from, to);
                for v in &mut e.initial_velocities {
                    swap(&mut v.set, from, to);
                }
                for w in &mut e.rigid_walls {
                    swap(&mut w.set, from, to);
                }
                if let Some(m) = &mut e.mass_scaling {
                    swap_opt(&mut m.set, from, to);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BoundaryCondition, Coupling, CouplingKind, Datum, Dofs, Load, Material, Partition, ReferencePoint, Section, Set, Step, Tie};

    /// A model whose definitions refer to each other: set `s` in a section, a coupling, a tie, a
    /// boundary condition and a load; material `m` in the section; reference point `rp` in the
    /// coupling and a load; both loads and the condition active in step `go`; datum `d` cut by.
    fn model() -> Model {
        let mut m = Model::default();
        for s in ["s", "t"] {
            m.add_set(s, Set { nodes: vec![1], ..Set::default() }).unwrap();
        }
        m.add_material("m", Material { elastic: Some(crate::Elastic::Isotropic { young: 1.0, poisson: 0.3 }), ..Material::default() }).unwrap();
        m.add_section("sec", Section { set: "s".into(), material: "m".into(), kind: crate::SectionKind::Solid {}, formulation: Default::default() })
            .unwrap();
        m.add_reference_point("rp", ReferencePoint { xyz: [0.0; 3], label: 0 }).unwrap();
        m.add_coupling("c", Coupling { reference_point: "rp".into(), set: "s".into(), kind: CouplingKind::Kinematic, dofs: Dofs::all() }).unwrap();
        m.add_tie("tie", Tie::new("s", "t")).unwrap();
        m.add_boundary_condition("fix", BoundaryCondition::Fixed { target: Target::Set("s".into()) }).unwrap();
        m.add_load("p", Load::Pressure { set: "s".into(), magnitude: 1.0 }).unwrap();
        m.add_load("f", Load::Force { target: Target::ReferencePoint("rp".into()), vector: [1.0, 0.0, 0.0], distribution: Default::default() })
            .unwrap();
        m.add_step(Step::new("go", Procedure::Static {}), None).unwrap();
        m.set_boundary_condition_active("go", "fix", true).unwrap();
        m.set_load_active("go", "p", true).unwrap();
        m.set_load_active("go", "f", true).unwrap();
        m.add_datum("d", Datum::Plane { point: [0.0; 3], normal: [0.0, 0.0, 1.0] }).unwrap();
        m.add_partition(Partition { name: "half".into(), datum: "d".into(), bodies: Vec::new(), suppressed: false }).unwrap();
        m
    }

    /// Whatever is renamed, every reference follows: the old name has no users left, the new
    /// one has all of them, and the model reads back from its file.
    #[test]
    fn renaming_rewrites_every_reference() {
        let mut m = model();
        let users = m.set_users("s");
        assert_eq!(users.len(), 5, "{users:?}");
        m.rename(Kind::Set, "s", "S2").unwrap();
        assert!(m.set_users("s").is_empty() && m.sets().get("s").is_none());
        assert_eq!(m.set_users("S2"), users);

        let users = m.material_users("m");
        m.rename(Kind::Material, "m", "steel").unwrap();
        assert_eq!(m.material_users("steel"), users);
        assert_eq!(m.section("sec").map(|s| s.material.as_str()), Some("steel"));

        let users = m.reference_point_users("rp");
        m.rename(Kind::ReferencePoint, "rp", "RP-1").unwrap();
        assert!(m.reference_point_users("rp").is_empty());
        assert_eq!(m.reference_point_users("RP-1"), users);

        m.rename(Kind::BoundaryCondition, "fix", "clamp").unwrap();
        m.rename(Kind::Load, "f", "push").unwrap();
        m.rename(Kind::Step, "go", "Load").unwrap();
        let step = m.steps().iter().find(|s| s.name == "Load").unwrap();
        assert!(step.boundary_conditions.contains("clamp") && step.loads.contains("push") && step.loads.contains("p"));

        m.rename(Kind::Datum, "d", "Mid").unwrap();
        assert_eq!(m.datum_users("Mid"), ["partition `half`"]);

        let json = serde_json::to_value(&m).unwrap();
        let back: Model = serde_json::from_value(json).unwrap();
        assert_eq!(back, m);
    }

    /// A taken, invalid or missing name changes nothing; the same name is a no-op.
    #[test]
    fn a_rename_that_cannot_be_done_changes_nothing() {
        let mut m = model();
        let before = m.clone();
        assert!(m.rename(Kind::Set, "s", "t").unwrap_err().to_string().contains("already exists"));
        assert!(m.rename(Kind::Set, "nope", "x").is_err());
        assert!(m.rename(Kind::Set, "s", "").is_err());
        assert!(m.rename(Kind::Step, "go", "go").is_ok());
        assert!(m.rename(Kind::Step, "gone", "gone").is_err());
        assert_eq!(m, before);
    }
}
