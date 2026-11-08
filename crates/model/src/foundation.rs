//! Linear surface foundations on the reference surface, in SI units.
use crate::{Error, Model, Result, insert_new, remove_existing, replace_existing};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Foundation {
    pub set: String,
    /// Normal stiffness per area (N/m³).
    pub normal: f64,
    /// Isotropic in-plane stiffness per area (N/m³).
    #[serde(default)]
    pub tangential: f64,
    /// Normal and tangential viscous coefficients per area (N s/m³).
    #[serde(default)]
    pub normal_damping: f64,
    #[serde(default)]
    pub tangential_damping: f64,
    /// Structural damping H = loss_factor K_foundation (dimensionless).
    #[serde(default)]
    pub loss_factor: f64,
}
impl Foundation {
    pub fn validate(&self) -> Result<()> {
        if self.set.trim().is_empty()
            || [self.normal, self.tangential, self.normal_damping, self.tangential_damping, self.loss_factor]
                .iter()
                .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(Error::Invalid("foundation needs a set and finite non-negative stiffness/damping coefficients".into()));
        }
        Ok(())
    }
}
impl Model {
    pub fn foundations(&self) -> &BTreeMap<String, Foundation> {
        &self.foundations
    }
    pub fn foundation(&self, name: &str) -> Option<&Foundation> {
        self.foundations.get(name)
    }
    fn check_foundation(&self, foundation: &Foundation) -> Result<()> {
        foundation.validate()?;
        if self.set(&foundation.set).is_none() {
            return Err(Error::NotFound(format!("foundation set `{}`", foundation.set)));
        }
        Ok(())
    }
    pub fn add_foundation(&mut self, name: &str, foundation: Foundation) -> Result<()> {
        self.check_foundation(&foundation)?;
        insert_new(&mut self.foundations, "foundation", name, foundation)
    }
    pub fn replace_foundation(&mut self, name: &str, foundation: Foundation) -> Result<()> {
        self.check_foundation(&foundation)?;
        replace_existing(&mut self.foundations, "foundation", name, foundation)
    }
    pub fn remove_foundation(&mut self, name: &str) -> Result<()> {
        remove_existing(&mut self.foundations, "foundation", name).map(|_| ())
    }
}
