//! Where the model names the design's geometry, and renaming a face, edge or vertex everywhere.
//!
//! Definitions refer to the design's entities by persistent name (`F2:top`, `F2:side:l1|F2:top`).
//! When the design is opened again, a name the design no longer has leaves the definition
//! naming nothing; [`Model::name_uses`] lists every such use so the opener can say which, and
//! [`Model::rename_entity`] moves every use of a name, and the mesh's tags, to the entity's new
//! name. The uses: sets, edge seeds, size controls, mesh plans, face merges and splits, fluid
//! domains' obstacle faces and fluid steps' monitors; boundary conditions, loads, contacts and
//! ties reach geometry through their sets. Defeatures name features by id, which do not change.

use crate::geometry::{Entity, GeometryRef};
use crate::{Error, FluidBoundary, Model, Procedure, Result, SizeRegion};

/// One place the model names a design entity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NameUse {
    /// What names it, for messages: ``set `cut` ``, ``edge seed``, ``size control `Hole` ``, …
    pub owner: String,
    pub reference: GeometryRef,
}

/// The face references a fluid boundary list holds.
fn boundary_faces(faces: &[FluidBoundary]) -> impl Iterator<Item = &GeometryRef> {
    faces.iter().filter_map(|f| match f {
        FluidBoundary::Surface(r) => Some(r),
        FluidBoundary::Box(_) => None,
    })
}

/// `name` as a reference to a face of `body` (invalid names are skipped: nothing can match them).
fn face_use(owner: &str, body: &str, name: &str, out: &mut Vec<NameUse>) {
    if let Ok(reference) = GeometryRef::face(body, name) {
        out.push(NameUse { owner: owner.to_string(), reference });
    }
}

impl Model {
    /// Every place a definition names a design entity, in a fixed order (the mesh's own tags are
    /// [`crate::Classification::entities`]).
    pub fn name_uses(&self) -> Vec<NameUse> {
        let mut out = Vec::new();
        let mut push = |owner: String, r: &GeometryRef| out.push(NameUse { owner, reference: r.clone() });
        for (name, set) in &self.sets {
            for r in &set.geometry {
                push(format!("set `{name}`"), r);
            }
        }
        for seed in &self.seeds {
            push("edge seed".into(), &seed.edge);
        }
        for (name, c) in &self.size_controls {
            if let SizeRegion::Geometry(refs) = &c.region {
                for r in refs {
                    push(format!("size control `{name}`"), r);
                }
            }
        }
        for (name, d) in &self.fluid_domains {
            let owner = format!("fluid domain `{name}`");
            for s in &d.surfaces {
                push(owner.clone(), &s.face);
            }
            for s in d.thermal.iter().flat_map(|t| &t.surfaces) {
                push(owner.clone(), &s.face);
            }
        }
        for step in &self.steps {
            if let Procedure::Fluid { steady: Some(steady), .. } = &step.procedure {
                for m in &steady.monitors {
                    if let crate::FlowMonitor::Force { faces } | crate::FlowMonitor::Heat { faces } = m {
                        for r in boundary_faces(faces) {
                            push(format!("step `{}` monitor", step.name), r);
                        }
                    }
                }
            }
        }
        for (name, plan) in &self.mesh_plans {
            let owner = format!("mesh plan `{name}`");
            for f in plan.features.iter().flat_map(|f| &f.faces).chain(plan.rings.iter().flat_map(|r| &r.faces)) {
                face_use(&owner, &plan.body, f, &mut out);
            }
        }
        for m in &self.face_merges {
            let owner = format!("face merge `{}`", m.name);
            for f in m.faces.iter().chain(&m.into) {
                face_use(&owner, &m.body, f, &mut out);
            }
        }
        for s in &self.face_splits {
            face_use(&format!("face split `{}`", s.name), &s.body, &s.face, &mut out);
        }
        out
    }

    /// Rename entity `from` to `to` (a face, edge or vertex of the same body, as the design names
    /// it now) in every use [`Model::name_uses`] lists and in the mesh's tags, so the mesh and
    /// every definition follow it. Only exact uses of `from` move (not the names of its pieces).
    /// Refused, changing nothing, when `to` is already a tag of the mesh or an edge seed's edge:
    /// two names would then become one, which the caller must decide. Returns the uses moved.
    pub fn rename_entity(&mut self, from: &GeometryRef, to: &GeometryRef) -> Result<usize> {
        let (Some(old), Some(new)) = (from.entity.name(), to.entity.name()) else {
            return Err(Error::Invalid("only faces, edges and vertices are renamed".into()));
        };
        if from.body != to.body || from.entity.kind() != to.entity.kind() {
            return Err(Error::Invalid(format!("{from} cannot become {to}: a rename keeps the body and the kind")));
        }
        if from == to {
            return Ok(0);
        }
        let tagged = self.classification.entities().contains(from);
        if tagged && self.classification.entities().contains(to) {
            return Err(Error::Invalid(format!("the mesh is tagged with both {from} and {to}")));
        }
        if self.seeds.iter().any(|s| &s.edge == from) && self.seeds.iter().any(|s| &s.edge == to) {
            return Err(Error::Invalid(format!("both {from} and {to} are seeded")));
        }
        let uses = self.name_uses();
        if let Some(u) = uses.iter().find(|u| &u.reference == from && uses.iter().any(|v| v.owner == u.owner && &v.reference == to)) {
            return Err(Error::Invalid(format!("{} names both {from} and {to}", u.owner)));
        }
        // The definitions as renamed must read back from a file, as `fieldcraft_io::read` replays
        // them through the same checks; tried on a copy without the mesh, which they don't need.
        let mut trial = self.definitions();
        trial.rename_in_definitions(from, to, old, new);
        let readable = serde_json::to_value(&trial)
            .map_err(|e| Error::Invalid(e.to_string()))
            .and_then(|v| <Model as serde::Deserialize>::deserialize(v).map(|_| ()).map_err(|e| Error::Invalid(e.to_string())));
        if let Err(e) = readable {
            return Err(Error::Invalid(format!("renaming {from} to {to} would leave a model that cannot be read back ({e})")));
        }
        let moved = self.rename_in_definitions(from, to, old, new);
        if tagged {
            self.classification.rename_entity(from, to);
        }
        Ok(moved)
    }

    /// The model without its mesh (nodes, elements and their tags): its definitions alone.
    fn definitions(&mut self) -> Model {
        let nodes = std::mem::take(&mut self.nodes);
        let elements = std::mem::take(&mut self.elements);
        let classification = std::mem::take(&mut self.classification);
        let copy = self.clone();
        self.nodes = nodes;
        self.elements = elements;
        self.classification = classification;
        copy
    }

    /// [`Model::rename_entity`] in the definitions (not the mesh's tags); the uses moved.
    fn rename_in_definitions(&mut self, from: &GeometryRef, to: &GeometryRef, old: &str, new: &str) -> usize {
        let mut moved = 0;
        let mut swap = |r: &mut GeometryRef| {
            if r == from {
                *r = to.clone();
                moved += 1;
            }
        };
        for set in self.sets.values_mut() {
            set.geometry.iter_mut().for_each(&mut swap);
        }
        for seed in &mut self.seeds {
            swap(&mut seed.edge);
        }
        self.seeds.sort_by(|a, b| a.edge.cmp(&b.edge));
        for c in self.size_controls.values_mut() {
            if let SizeRegion::Geometry(refs) = &mut c.region {
                refs.iter_mut().for_each(&mut swap);
            }
        }
        for d in self.fluid_domains.values_mut() {
            d.surfaces.iter_mut().for_each(|s| swap(&mut s.face));
            for s in d.thermal.iter_mut().flat_map(|t| &mut t.surfaces) {
                swap(&mut s.face);
            }
        }
        for step in &mut self.steps {
            if let Procedure::Fluid { steady: Some(steady), .. } = &mut step.procedure {
                for m in &mut steady.monitors {
                    if let crate::FlowMonitor::Force { faces } | crate::FlowMonitor::Heat { faces } = m {
                        for f in faces {
                            if let FluidBoundary::Surface(r) = f {
                                swap(r);
                            }
                        }
                    }
                }
            }
        }
        // Faces held as plain names beside their body.
        if matches!(from.entity, Entity::Face(_)) {
            let mut rename = |body: &str, name: &mut String| {
                if body == from.body && name == old {
                    *name = new.to_string();
                    moved += 1;
                }
            };
            for plan in self.mesh_plans.values_mut() {
                let body = plan.body.clone();
                for f in plan.features.iter_mut().flat_map(|f| &mut f.faces).chain(plan.rings.iter_mut().flat_map(|r| &mut r.faces)) {
                    rename(&body, f);
                }
            }
            for m in &mut self.face_merges {
                for f in m.faces.iter_mut().chain(&mut m.into) {
                    rename(&m.body, f);
                }
            }
            for s in &mut self.face_splits {
                rename(&s.body, &mut s.face);
            }
        }
        moved
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{FaceMerge, FaceSplit, Set, SizeControl};

    fn face(n: &str) -> GeometryRef {
        GeometryRef::face("B", n).unwrap()
    }

    fn edge(n: &str) -> GeometryRef {
        GeometryRef::new("B", Entity::Edge(n.into())).unwrap()
    }

    #[test]
    fn a_rename_moves_every_use_and_the_mesh_tag() {
        let mut m = Model::default();
        m.add_set("cut", Set::of_geometry(vec![face("F9:face:4"), face("F1:top")])).unwrap();
        m.set_size_control("Fine", SizeControl { region: SizeRegion::Geometry(vec![face("F9:face:4")]), size: 1.0 }).unwrap();
        m.add_face_merge(FaceMerge { name: "M".into(), body: "B".into(), faces: vec!["F9:face:4".into(), "F1:top".into()], into: None }).ok();
        m.add_datum("D", crate::Datum::Plane { point: [0.0; 3], normal: [0.0, 0.0, 1.0] }).unwrap();
        m.add_face_split(FaceSplit { name: "S".into(), body: "B".into(), face: "F9:face:4".into(), datum: "D".into() }).unwrap();
        let tag = m.geometry_tag(&face("F9:face:4")).unwrap();
        let before = m.name_uses();
        assert!(before.iter().any(|u| u.owner == "set `cut`" && u.reference == face("F9:face:4")));
        let moved = m.rename_entity(&face("F9:face:4"), &face("F9:section")).unwrap();
        let after = m.name_uses();
        // Every use of the old name moved, nothing else did.
        assert_eq!(moved, before.iter().filter(|u| u.reference == face("F9:face:4")).count(), "{before:?}");
        assert!(after.iter().all(|u| u.reference != face("F9:face:4")), "{after:?}");
        assert_eq!(after.len(), before.len());
        assert_eq!(m.classification().entity(tag), Some(&face("F9:section")));
        assert_eq!(m.geometry_tag(&face("F9:section")).unwrap(), tag);
    }

    /// A merge of `A` and `C` would list `C` twice after `A` became `C`: the file saves, and then
    /// never reads back ("face `C` is listed twice"). Refused, changing nothing.
    #[test]
    fn a_rename_onto_a_name_its_owner_lists_already_is_refused() {
        let mut m = Model::default();
        m.add_face_merge(FaceMerge { name: "M".into(), body: "B".into(), faces: vec!["A".into(), "C".into()], into: None }).unwrap();
        let before = m.clone();
        let e = m.rename_entity(&face("A"), &face("C")).unwrap_err().to_string();
        assert!(e.contains("face merge `M` names both"), "{e}");
        assert_eq!(m, before);
        // The read-back check stands behind it: a merge whose list holds the new name twice.
        let mut definitions = m.definitions();
        definitions.rename_in_definitions(&face("A"), &face("C"), "A", "C");
        let v = serde_json::to_value(&definitions).unwrap();
        assert!(<Model as serde::Deserialize>::deserialize(v).is_err());
    }

    #[test]
    fn a_rename_that_would_join_two_names_is_refused() {
        let mut m = Model::default();
        m.geometry_tag(&face("A")).unwrap();
        m.geometry_tag(&face("B")).unwrap();
        assert!(m.rename_entity(&face("A"), &face("B")).is_err());
        assert!(m.rename_entity(&face("A"), &edge("A|C")).is_err());
        assert_eq!(m.classification().entities(), &[face("A"), face("B")]);
    }
}
