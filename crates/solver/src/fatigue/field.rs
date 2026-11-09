//! Stored nodal fatigue field. Missing coverage is explicit; no interpolation of
//! uncomputed damage. Life is repetition of the counted record, in SI seconds.
use crate::{Error, Result};
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct FatigueField {
    pub damage: Vec<Option<f64>>,
    pub life_seconds: Vec<Option<f64>>,
    /// Exact command inputs, retained for reproducibility, not executed on load.
    pub inputs_json: String,
    /// Definition hashes of every reference step (including multichannel sources).
    pub dependencies: Vec<(String, [u8; 32])>,
}
impl FatigueField {
    pub fn validate(&self, nodes: usize) -> Result<()> {
        if self.damage.len() != nodes
            || self.life_seconds.len() != nodes
            || self.damage.iter().flatten().any(|v| !v.is_finite() || *v < 0.)
            || self.life_seconds.iter().flatten().any(|v| !v.is_finite() || *v <= 0.)
            || self.damage.iter().zip(&self.life_seconds).any(|(d, l)| l.is_some() && d.is_none_or(|d| d <= 0.))
            || self.inputs_json.len() > 16 * 1024 * 1024
            || self.dependencies.is_empty()
            || self.dependencies.len() > 64
        {
            return Err(Error::Invalid("invalid stored fatigue field dimensions, values or provenance".into()));
        }
        Ok(())
    }
}
