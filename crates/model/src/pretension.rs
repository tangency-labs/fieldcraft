//! A cut with a constant axial displacement jump, controlled by a scalar reference DOF.
//! Virtual work: q = n·(u_first-u_second), so a positive conjugate force tightens
//! the two bolt halves. A prescribed q locks the adjustment for a later equilibrium.
use crate::{Error, Model, Result, insert_new};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PretensionSection {
    /// Conforming, duplicate-node cut surfaces, resolved as node sets.
    pub first: String,
    pub second: String,
    /// Unit direction from the first bolt half to the second, fixed in reference space.
    pub normal: [f64; 3],
    /// Dedicated reference point; Ux represents the scalar cut adjustment, not global X.
    pub control: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum PretensionControl {
    Force {
        value: f64,
    },
    Adjustment {
        value: f64,
    },
    /// Hold the adjustment reached by the preceding equilibrium state.
    Fixed,
}
impl PretensionControl {
    pub fn validate(&self) -> Result<()> {
        let value = match self {
            Self::Force { value } | Self::Adjustment { value } => *value,
            Self::Fixed => return Ok(()),
        };
        if !value.is_finite() {
            return Err(Error::Invalid("pretension value must be finite".into()));
        }
        Ok(())
    }
}
impl Model {
    pub fn pretensions(&self) -> &BTreeMap<String, PretensionSection> {
        &self.pretensions
    }
    pub fn add_pretension(&mut self, name: &str, section: PretensionSection) -> Result<()> {
        if self.pretensions.len() >= 1000 {
            return Err(Error::Limit("at most 1000 pretension sections".into()));
        }
        if section.first == section.second || self.set(&section.first).is_none() || self.set(&section.second).is_none() {
            return Err(Error::Invalid("pretension requires two different existing cut sets".into()));
        }
        let norm = section.normal.iter().map(|v| v * v).sum::<f64>();
        if !norm.is_finite() || (norm - 1.0).abs() > 1e-10 {
            return Err(Error::Invalid("pretension normal must be a finite unit vector".into()));
        }
        if self.reference_point(&section.control).is_none() || !self.reference_point_users(&section.control).is_empty() {
            return Err(Error::Invalid("pretension control must be an existing unused reference point".into()));
        }
        insert_new(&mut self.pretensions, "pretension section", name, section)
    }
    pub fn remove_pretension(&mut self, name: &str) -> Result<()> {
        if self.steps.iter().any(|s| s.pretension.contains_key(name)) {
            return Err(Error::Invalid(format!("pretension section `{name}` is controlled in a step")));
        }
        self.pretensions.remove(name).map(|_| ()).ok_or_else(|| Error::NotFound(format!("pretension section `{name}`")))
    }
}
