//! Named sets: regions of the model that sections, couplings, boundary conditions and loads
//! apply to.
//!
//! A set is defined by geometry names, explicit labels or both; it is resolved against the
//! current mesh into nodes, elements and element faces ([`ResolvedSet`]). Defining sets by
//! geometry is the norm: the definition outlives the mesh, so a remesh (or an edit of the design
//! followed by a remesh) keeps every set, section and load in place. Explicit labels are for
//! meshes with no geometry, such as one read from a deck.
//!
//! What each member contributes:
//! - a body: the elements filling it, their nodes, and the element faces on the body's faces;
//! - a face: the element faces on it, the shell elements filling it, and their nodes;
//! - an edge: the element edges on it (what 2D models load), their nodes, the line elements
//!   filling it, and their nodes;
//! - a vertex: the node on it;
//! - node labels: those nodes; element labels: those elements and their nodes; element faces:
//!   those faces and their nodes; element edges: those edges and their nodes.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::check::ResolvedSet;
use crate::mesh::check_label;
use crate::{
    ElementEdge, ElementFace, Error, GeometryRef, MAX_ELEMENTS, MAX_NODES, Model, Result, context, insert_new, not_in_use, remove_existing,
    replace_existing,
};

/// Most geometry references in one set.
pub const MAX_SET_REFERENCES: usize = 100_000;

/// A set's definition. Members are kept sorted and without repeats.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Set {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub geometry: Vec<GeometryRef>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub nodes: Vec<u32>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub elements: Vec<u32>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub faces: Vec<ElementFace>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub edges: Vec<ElementEdge>,
}

impl Set {
    /// A set of geometry entities.
    pub fn of_geometry(geometry: Vec<GeometryRef>) -> Set {
        Set { geometry, ..Set::default() }
    }

    pub fn is_empty(&self) -> bool {
        self.geometry.is_empty() && self.nodes.is_empty() && self.elements.is_empty() && self.faces.is_empty() && self.edges.is_empty()
    }

    /// Check the definition's size and labels, and sort its members.
    pub fn normalized(mut self) -> Result<Set> {
        if self.is_empty() {
            return Err(Error::Invalid("a set needs geometry, nodes, elements, element faces or element edges".into()));
        }
        if self.geometry.len() > MAX_SET_REFERENCES {
            return Err(Error::Limit(format!("a set names at most {MAX_SET_REFERENCES} geometry entities")));
        }
        if self.nodes.len() > MAX_NODES || self.elements.len() > MAX_ELEMENTS || self.faces.len() > MAX_ELEMENTS || self.edges.len() > MAX_ELEMENTS {
            return Err(Error::Limit("a set lists at most as many members as a model holds".into()));
        }
        for &l in &self.nodes {
            check_label("node", l)?;
        }
        for &l in self.elements.iter().chain(self.faces.iter().map(|f| &f.element)).chain(self.edges.iter().map(|e| &e.element)) {
            check_label("element", l)?;
        }
        if let Some(f) = self.faces.iter().find(|f| f.face == 0) {
            return Err(Error::Invalid(format!("element {} face 0: faces are numbered from 1", f.element)));
        }
        if let Some(e) = self.edges.iter().find(|e| e.edge == 0) {
            return Err(Error::Invalid(format!("element {} edge 0: edges are numbered from 1", e.element)));
        }
        for v in [&mut self.nodes, &mut self.elements] {
            v.sort_unstable();
            v.dedup();
        }
        self.faces.sort_unstable();
        self.faces.dedup();
        self.edges.sort_unstable();
        self.edges.dedup();
        self.geometry.sort();
        self.geometry.dedup();
        Ok(self)
    }
}

impl Model {
    pub fn sets(&self) -> &BTreeMap<String, Set> {
        &self.sets
    }

    pub fn set(&self, name: &str) -> Option<&Set> {
        self.sets.get(name)
    }

    pub fn add_set(&mut self, name: &str, set: Set) -> Result<()> {
        let set = set.normalized().map_err(|e| context("set", name, e))?;
        insert_new(&mut self.sets, "set", name, set)
    }

    pub fn replace_set(&mut self, name: &str, set: Set) -> Result<()> {
        let set = set.normalized().map_err(|e| context("set", name, e))?;
        replace_existing(&mut self.sets, "set", name, set)
    }

    /// The sets that list, by label, a node or an element of `removed` (after a remesh they
    /// name labels that are gone or now mean other mesh entities; sets by geometry re-resolve).
    pub fn sets_using(&self, removed: &crate::BodyMesh) -> Vec<String> {
        let element = |e: &u32| removed.elements.contains(e);
        self.sets
            .iter()
            .filter(|(_, s)| {
                s.nodes.iter().any(|n| removed.nodes.contains(n))
                    || s.elements.iter().any(element)
                    || s.faces.iter().any(|f| element(&f.element))
                    || s.edges.iter().any(|e| element(&e.element))
            })
            .map(|(name, _)| name.clone())
            .collect()
    }

    pub fn remove_set(&mut self, name: &str) -> Result<()> {
        not_in_use("set", name, self.set_users(name))?;
        remove_existing(&mut self.sets, "set", name).map(|_| ())
    }

    /// The definitions that refer to a set, as `kind \`name\`` strings.
    pub fn set_users(&self, set: &str) -> Vec<String> {
        let mut users = self.set_model_users(set);
        users
            .extend(self.boundary_conditions.iter().filter(|(_, b)| b.target().set() == Some(set)).map(|(n, _)| format!("boundary condition `{n}`")));
        users.extend(self.loads.iter().filter(|(_, l)| l.set() == Some(set)).map(|(n, _)| format!("load `{n}`")));
        users.extend(self.steps.iter().filter(|s| s.procedure.sets().any(|x| x == set)).map(|s| format!("step `{}`", s.name)));
        users
    }

    /// The users of a set that act in every step (sections, couplings, contacts, ties, point
    /// masses and springs): all but boundary conditions, loads and steps, which act only in the
    /// steps that list them ([`Model::step_inputs`]).
    fn set_model_users(&self, set: &str) -> Vec<String> {
        let mut users = Vec::new();
        users.extend(self.sections.iter().filter(|(_, s)| s.set == set).map(|(n, _)| format!("section `{n}`")));
        users.extend(self.couplings.iter().filter(|(_, c)| c.set == set).map(|(n, _)| format!("coupling `{n}`")));
        users.extend(self.contact_users(set));
        users.extend(self.tie_users(set));
        users.extend(self.foundations.iter().filter(|(_, f)| f.set == set).map(|(n, _)| format!("foundation `{n}`")));
        users.extend(self.discrete_users(|t| t.set() == Some(set)));
        users
    }

    /// Whether solving `step` reads set `set`: a definition acting in every step refers to it,
    /// or one of the step's own boundary conditions and loads, or its procedure.
    pub(crate) fn step_uses_set(&self, step: &crate::Step, set: &str) -> bool {
        !self.set_model_users(set).is_empty()
            || step.boundary_conditions.iter().any(|n| self.boundary_conditions.get(n).is_some_and(|b| b.target().set() == Some(set)))
            || step.loads.iter().any(|n| self.loads.get(n).is_some_and(|l| l.set() == Some(set)))
            || step.procedure.sets().any(|x| x == set)
    }

    /// Resolve one set against the current mesh. The error names every member that is not on
    /// the mesh.
    pub fn resolve_set(&self, name: &str) -> Result<ResolvedSet> {
        let set = self.sets.get(name).ok_or_else(|| Error::NotFound(format!("set `{name}`")))?;
        let (resolved, missing) = self.resolve_definition(set);
        if missing.is_empty() { Ok(resolved) } else { Err(Error::NotFound(format!("set `{name}`: {}", missing.join("; ")))) }
    }

    /// Resolve one set against the current mesh, leaving out the members that are not on it
    /// (for display: a set with a missing member still shows what it has). `None` when there is
    /// no set called `name`.
    pub fn resolve_set_partial(&self, name: &str) -> Option<ResolvedSet> {
        self.sets.get(name).map(|set| self.resolve_definition(set).0)
    }

    /// Resolve a definition, and list the members that are not on the mesh.
    pub(crate) fn resolve_definition(&self, set: &Set) -> (ResolvedSet, Vec<String>) {
        let c = &self.classification;
        let mut nodes = Vec::new();
        let mut elements = Vec::new();
        let mut faces = Vec::new();
        let mut edges = Vec::new();
        let mut missing = Vec::new();

        // Geometry: mark the selected tags by role, then make one pass over the tagged elements.
        let n = c.entities().len();
        let (mut regions, mut face_tags, mut edge_tags, mut vertex_tags) = (vec![false; n], vec![false; n], vec![false; n], vec![false; n]);
        let mark = |marks: &mut Vec<bool>, tag: u32| {
            if let Some(m) = marks.get_mut(tag as usize) {
                *m = true;
            }
        };
        for r in &set.geometry {
            let tags = c.tags_selected_by(r);
            if tags.is_empty() {
                missing.push(if c.is_empty() {
                    format!("{r} is not on the mesh (the mesh carries no geometry tags)")
                } else {
                    format!("{r} is not on the mesh")
                });
                continue;
            }
            for &t in &tags {
                match r.entity.dimension() {
                    0 => mark(&mut vertex_tags, t),
                    1 => {
                        mark(&mut regions, t);
                        mark(&mut edge_tags, t);
                    }
                    2 => {
                        mark(&mut regions, t);
                        mark(&mut face_tags, t);
                    }
                    _ => mark(&mut regions, t),
                }
            }
            if r.entity.dimension() == 3 {
                // A body's skin: every tagged face of the body.
                for (i, e) in c.entities().iter().enumerate() {
                    if e.body == r.body && e.entity.dimension() == 2 {
                        mark(&mut face_tags, u32::try_from(i).unwrap_or(u32::MAX));
                    }
                }
            }
        }
        let on = |marks: &[bool], tag: u32| marks.get(tag as usize).copied().unwrap_or(false);
        if !set.geometry.is_empty() {
            for (&label, tags) in c.elements() {
                let Some(element) = self.elements.get(&label) else { continue };
                if tags.region.is_some_and(|t| on(&regions, t)) {
                    elements.push(label);
                    nodes.extend_from_slice(&element.nodes);
                }
                for (&face, &t) in &tags.faces {
                    if on(&face_tags, t) {
                        faces.push(ElementFace { element: label, face });
                    }
                }
                for (&edge, &t) in &tags.edges {
                    if on(&edge_tags, t) {
                        edges.push(ElementEdge { element: label, edge });
                    }
                }
            }
            nodes.extend(c.nodes().iter().filter(|(_, t)| on(&vertex_tags, **t)).map(|(n, _)| n));
        }

        // Explicit labels.
        for &n in &set.nodes {
            if self.nodes.contains_key(&n) {
                nodes.push(n);
            } else {
                missing.push(format!("node {n} is not in the mesh"));
            }
        }
        for &l in &set.elements {
            match self.elements.get(&l) {
                Some(e) => {
                    elements.push(l);
                    nodes.extend_from_slice(&e.nodes);
                }
                None => missing.push(format!("element {l} is not in the mesh")),
            }
        }
        for &f in &set.faces {
            match self.elements.get(&f.element) {
                Some(e) if usize::from(f.face) <= e.kind.faces().len() => faces.push(f),
                Some(e) => missing.push(format!("element {} ({}) has no face {}", f.element, e.kind.name(), f.face)),
                None => missing.push(format!("element {} is not in the mesh", f.element)),
            }
        }
        for &g in &set.edges {
            match self.elements.get(&g.element) {
                Some(e) if usize::from(g.edge) <= e.kind.edges().len() => edges.push(g),
                Some(e) => missing.push(format!("element {} ({}) has no edge {}", g.element, e.kind.name(), g.edge)),
                None => missing.push(format!("element {} is not in the mesh", g.element)),
            }
        }
        elements.sort_unstable();
        elements.dedup();
        faces.sort_unstable();
        faces.dedup();
        edges.sort_unstable();
        edges.dedup();
        for f in &faces {
            if let Some(e) = self.elements.get(&f.element) {
                let local = usize::from(f.face).checked_sub(1).and_then(|i| e.kind.faces().get(i)).copied().unwrap_or_default();
                nodes.extend(local.iter().filter_map(|&i| e.nodes.get(i)));
            }
        }
        for g in &edges {
            nodes.extend(self.edge_nodes(*g).unwrap_or_default());
        }
        nodes.sort_unstable();
        nodes.dedup();
        // Long lists of missing labels are summarised.
        if missing.len() > 20 {
            let more = missing.len() - 20;
            missing.truncate(20);
            missing.push(format!("and {more} more"));
        }
        (ResolvedSet { nodes, elements, faces, edges }, missing)
    }

    /// The node labels of an element edge, in the edge's order ([`crate::ElementKind::edges`]:
    /// its ends, then the mid-side node of a quadratic element).
    pub fn edge_nodes(&self, edge: ElementEdge) -> Option<Vec<u32>> {
        let e = self.elements.get(&edge.element)?;
        let local = e.kind.edges().get(usize::from(edge.edge).checked_sub(1)?)?;
        local.iter().map(|&i| e.nodes.get(i).copied()).collect()
    }

    /// The node labels of an element face, in the face's order ([`crate::ElementKind::faces`]).
    pub fn face_nodes(&self, face: ElementFace) -> Option<Vec<u32>> {
        let e = self.elements.get(&face.element)?;
        let local = e.kind.faces().get(usize::from(face.face).checked_sub(1)?)?;
        local.iter().map(|&i| e.nodes.get(i).copied()).collect()
    }
}
