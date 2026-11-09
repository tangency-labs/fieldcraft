//! Reference points and couplings.
//!
//! A reference point is a single node off the mesh with three translations and three rotations.
//! A coupling ties it to a set's nodes, so a force, moment or boundary condition on the point acts
//! on a whole face (a bolt hole, a mounting pad) without hand-picked nodes.
//!
//! Each reference point keeps the node label it is written under in every deck. Unless one is
//! given, it gets the highest free label at the top of the label range when it is created, far
//! above the labels meshers use (they count up from 1), so remeshing, adding nodes and adding or
//! removing other reference points never renumber it. A mesh node on the same label is an error
//! of the check.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::bcs::Dofs;
use crate::{Error, MAX_LABEL, Model, Result, check_coordinates, context, insert_new, not_in_use, remove_existing, replace_existing, require};

/// A reference point.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferencePoint {
    pub xyz: [f64; 3],
    /// Its node label in decks; 0 when creating means "choose one" (see the module docs).
    #[serde(default)]
    pub label: u32,
}

/// How a coupling ties the reference point to the set's nodes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CouplingKind {
    /// The nodes follow the reference point rigidly in the coupled degrees of freedom.
    Kinematic,
    /// The reference point's load is spread over the nodes (weighted equally), and its motion is
    /// their weighted average; the region stays free to deform.
    Distributing,
}

/// A coupling between a reference point and a set's nodes.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Coupling {
    pub reference_point: String,
    pub set: String,
    pub kind: CouplingKind,
    /// The coupled degrees of freedom (all six by default).
    #[serde(default = "Dofs::all")]
    pub dofs: Dofs,
}

impl Model {
    pub fn reference_points(&self) -> &BTreeMap<String, ReferencePoint> {
        &self.reference_points
    }

    pub fn reference_point(&self, name: &str) -> Option<&ReferencePoint> {
        self.reference_points.get(name)
    }

    /// Check a reference point's coordinates and label, choosing a label when it has none.
    fn place_reference_point(&self, name: &str, mut point: ReferencePoint, current: Option<u32>) -> Result<ReferencePoint> {
        check_coordinates(&point.xyz).map_err(|e| context("reference point", name, e))?;
        let taken = |l: u32| self.reference_points.iter().any(|(n, p)| n != name && p.label == l);
        if point.label == 0 {
            point.label = match current {
                Some(l) => l,
                None => (1..=MAX_LABEL)
                    .rev()
                    .take(crate::MAX_DEFINITIONS + 1)
                    .find(|l| !taken(*l))
                    .ok_or_else(|| Error::Limit(format!("no node label left for reference point `{name}`")))?,
            };
        }
        if point.label > MAX_LABEL {
            return Err(context("reference point", name, Error::Invalid(format!("label {} is outside 1..={MAX_LABEL}", point.label))));
        }
        if let Some((other, _)) = self.reference_points.iter().find(|(n, p)| *n != name && p.label == point.label) {
            return Err(context("reference point", name, Error::Invalid(format!("label {} is taken by reference point `{other}`", point.label))));
        }
        Ok(point)
    }

    pub fn add_reference_point(&mut self, name: &str, point: ReferencePoint) -> Result<()> {
        let point = self.place_reference_point(name, point, None)?;
        insert_new(&mut self.reference_points, "reference point", name, point)
    }

    /// Replace a reference point; a label of 0 keeps its current one.
    pub fn replace_reference_point(&mut self, name: &str, point: ReferencePoint) -> Result<()> {
        let current = self.reference_points.get(name).map(|p| p.label).ok_or_else(|| Error::NotFound(format!("reference point `{name}`")))?;
        let point = self.place_reference_point(name, point, Some(current))?;
        replace_existing(&mut self.reference_points, "reference point", name, point)
    }

    pub fn remove_reference_point(&mut self, name: &str) -> Result<()> {
        not_in_use("reference point", name, self.reference_point_users(name))?;
        remove_existing(&mut self.reference_points, "reference point", name).map(|_| ())
    }

    /// The definitions that refer to a reference point, as `kind \`name\`` strings.
    pub fn reference_point_users(&self, name: &str) -> Vec<String> {
        let mut users: Vec<String> =
            self.couplings.iter().filter(|(_, c)| c.reference_point == name).map(|(n, _)| format!("coupling `{n}`")).collect();
        users.extend(
            self.boundary_conditions
                .iter()
                .filter(|(_, b)| b.target().reference_point() == Some(name))
                .map(|(n, _)| format!("boundary condition `{n}`")),
        );
        users.extend(self.loads.iter().filter(|(_, l)| l.reference_point() == Some(name)).map(|(n, _)| format!("load `{n}`")));
        users.extend(self.pretensions.iter().filter(|(_, p)| p.control == name).map(|(n, _)| format!("pretension `{n}`")));
        users.extend(self.discrete_users(|t| t.reference_point() == Some(name)));
        users
    }

    /// The node label of each reference point in decks.
    pub fn reference_point_nodes(&self) -> BTreeMap<String, u32> {
        self.reference_points.iter().map(|(n, p)| (n.clone(), p.label)).collect()
    }

    pub fn couplings(&self) -> &BTreeMap<String, Coupling> {
        &self.couplings
    }

    pub fn coupling(&self, name: &str) -> Option<&Coupling> {
        self.couplings.get(name)
    }

    fn check_coupling(&self, name: &str, c: &Coupling) -> Result<()> {
        if c.dofs.is_empty() {
            return Err(context("coupling", name, Error::Invalid("couples no degree of freedom".into())));
        }
        require(&self.reference_points, "reference point", &c.reference_point).map_err(|e| context("coupling", name, e))?;
        require(&self.sets, "set", &c.set).map_err(|e| context("coupling", name, e))
    }

    pub fn add_coupling(&mut self, name: &str, coupling: Coupling) -> Result<()> {
        self.check_coupling(name, &coupling)?;
        insert_new(&mut self.couplings, "coupling", name, coupling)
    }

    pub fn replace_coupling(&mut self, name: &str, coupling: Coupling) -> Result<()> {
        self.check_coupling(name, &coupling)?;
        replace_existing(&mut self.couplings, "coupling", name, coupling)
    }

    pub fn remove_coupling(&mut self, name: &str) -> Result<()> {
        remove_existing(&mut self.couplings, "coupling", name).map(|_| ())
    }
}
