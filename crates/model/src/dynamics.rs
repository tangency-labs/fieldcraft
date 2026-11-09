//! Linear transient and pseudo-acceleration response-spectrum contracts in SI.
use crate::{Error, Result, amplitude::Amplitude, harmonic::Damping};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Integration {
    #[default]
    AverageAcceleration,
    LinearAcceleration,
    Hht {
        alpha: f64,
    },
}
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Method {
    #[default]
    Direct,
    Modal {
        modes: u32,
        #[serde(default)]
        damping_ratios: Vec<f64>,
    },
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BaseAcceleration {
    pub direction: [f64; 3],
    pub magnitude: f64,
    pub amplitude: Amplitude,
}
/// Prescribed uniform base translation. Initial base velocity/acceleration are
/// explicit because a displacement table alone does not specify those states.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BaseDisplacement {
    pub direction: [f64; 3],
    pub magnitude: f64,
    pub amplitude: Amplitude,
    pub initial_velocity: f64,
    pub initial_acceleration: f64,
}
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ErrorEstimator {
    #[default]
    ZienkiewiczXie,
    QuadraticAcceleration,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Adaptive {
    pub min_step: f64,
    pub max_step: f64,
    pub absolute_tolerance: f64,
    #[serde(default)]
    pub relative_tolerance: f64,
    pub max_attempts: u32,
    #[serde(default)]
    pub estimator: ErrorEstimator,
}
impl Adaptive {
    pub fn validate(&self, interval: f64) -> Result<()> {
        if !self.min_step.is_finite()
            || !self.max_step.is_finite()
            || self.min_step <= 0.0
            || self.max_step < self.min_step
            || self.max_step > interval
            || !self.absolute_tolerance.is_finite()
            || self.absolute_tolerance <= 0.0
            || !self.relative_tolerance.is_finite()
            || !(0.0..=1.0).contains(&self.relative_tolerance)
            || !(1..=1_000_000).contains(&self.max_attempts)
        {
            return Err(Error::Invalid("adaptive transient requires 0 < min_step <= max_step <= time_step, positive absolute_tolerance, relative_tolerance in [0,1], and 1..1000000 max_attempts".into()));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Transient {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub adaptive: Option<Adaptive>,
    pub time_step: f64,
    pub increments: u32,
    #[serde(default = "one")]
    pub output_every: u32,
    #[serde(default)]
    pub integration: Integration,
    #[serde(default)]
    pub method: Method,
    #[serde(default)]
    pub damping: Damping,
    #[serde(default)]
    pub load_amplitudes: BTreeMap<String, Amplitude>,
    /// Node label to XYZ translation, in m or m/s. Rotational initial states unsupported.
    #[serde(default, deserialize_with = "node_vectors")]
    pub initial_displacement: BTreeMap<u32, [f64; 3]>,
    #[serde(default, deserialize_with = "node_vectors")]
    pub initial_velocity: BTreeMap<u32, [f64; 3]>,
    #[serde(default)]
    pub base_acceleration: Option<BaseAcceleration>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub base_displacement: Option<BaseDisplacement>,
}
// Tagged procedure deserialization buffers map keys as strings. Decode labels
// explicitly so the command path and direct JSON loading have the same contract.
fn node_vectors<'de, D: serde::Deserializer<'de>>(deserializer: D) -> std::result::Result<BTreeMap<u32, [f64; 3]>, D::Error> {
    let raw = BTreeMap::<String, [f64; 3]>::deserialize(deserializer)?;
    let mut result = BTreeMap::new();
    for (key, value) in raw {
        let label = key.parse::<u32>().map_err(serde::de::Error::custom)?;
        if key != label.to_string() {
            return Err(serde::de::Error::custom("node labels must be canonical unsigned integers"));
        }
        result.insert(label, value);
    }
    Ok(result)
}
fn one() -> u32 {
    1
}
pub fn direction(d: [f64; 3]) -> Result<()> {
    let norm: f64 = d.iter().map(|v| v * v).sum();
    if d.iter().any(|v| !v.is_finite()) || (norm - 1.0).abs() > 1e-10 {
        return Err(Error::Invalid("base direction must be a finite unit vector".into()));
    }
    Ok(())
}
pub fn ratios(modes: u32, values: &[f64]) -> Result<()> {
    if !(1..=4096).contains(&modes)
        || ![0, 1, modes as usize].contains(&values.len())
        || values.iter().any(|x| !x.is_finite() || *x < 0.0 || *x >= 1.0)
    {
        return Err(Error::Invalid("modal method needs 1..4096 modes and zero, one or modes damping ratios in [0,1)".into()));
    }
    Ok(())
}
impl Transient {
    pub fn validate(&self) -> Result<()> {
        if !self.time_step.is_finite()
            || self.time_step <= 0.0
            || !(1..=1_000_000).contains(&self.increments)
            || !(self.time_step * self.increments as f64).is_finite()
            || self.output_every == 0
            || self.output_every > self.increments
        {
            return Err(Error::Invalid(
                "transient requires positive finite time_step, 1..1000000 increments and output_every in 1..increments".into(),
            ));
        }
        if let Integration::Hht { alpha } = self.integration
            && (!alpha.is_finite() || !(-1.0 / 3.0..=0.0).contains(&alpha))
        {
            return Err(Error::Invalid("HHT alpha in [-1/3,0]".into()));
        }
        if let Method::Modal { modes, damping_ratios } = &self.method {
            ratios(*modes, damping_ratios)?;
        }
        if let Some(a) = &self.adaptive {
            a.validate(self.time_step)?;
            if !matches!(self.integration, Integration::AverageAcceleration)
                || !matches!(self.method, Method::Direct)
                || self.base_displacement.is_some()
            {
                return Err(Error::Invalid("adaptive transient currently supports direct average_acceleration without base_displacement".into()));
            }
        }
        self.damping.validate()?;
        if self.damping.loss_factor != 0.0 {
            return Err(Error::Invalid("transient structural loss factor is unsupported; use viscous damping".into()));
        }
        if self.load_amplitudes.len() > 10000 || self.initial_displacement.len() > 10000 || self.initial_velocity.len() > 10000 {
            return Err(Error::Invalid("transient input lists exceed 10000 entries".into()));
        }
        for a in self.load_amplitudes.values() {
            a.validate()?;
        }
        for (label, v) in self.initial_displacement.iter().chain(&self.initial_velocity) {
            if *label == 0 || v.iter().any(|x| !x.is_finite()) {
                return Err(Error::Invalid("invalid transient initial node or state".into()));
            }
        }
        if let Some(base) = &self.base_acceleration {
            direction(base.direction)?;
            if !base.magnitude.is_finite() {
                return Err(Error::Invalid("base acceleration magnitude must be finite".into()));
            }
            base.amplitude.validate()?;
        }
        if self.base_acceleration.is_some() && self.base_displacement.is_some() {
            return Err(Error::Invalid("choose base_acceleration or base_displacement, not both".into()));
        }
        if let Some(base) = &self.base_displacement {
            direction(base.direction)?;
            if [base.magnitude, base.initial_velocity, base.initial_acceleration].iter().any(|x| !x.is_finite()) {
                return Err(Error::Invalid("base displacement and initial derivatives must be finite".into()));
            }
            base.amplitude.validate()?;
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Combination {
    Srss,
    Cqc,
}
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DirectionCombination {
    #[default]
    Srss,
    Percent1004040,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SpectrumTable {
    pub damping_ratio: f64,
    pub points: Vec<[f64; 2]>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SpectrumDirection {
    pub direction: [f64; 3],
    pub scale: f64,
}
fn x_direction() -> [f64; 3] {
    [1.0, 0.0, 0.0]
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Spectrum {
    pub modes: u32,
    /// Used for every mode unless modal_damping_ratios overrides it.
    pub damping_ratio: f64,
    #[serde(default)]
    pub modal_damping_ratios: Vec<f64>,
    pub combination: Combination,
    #[serde(default = "x_direction")]
    pub direction: [f64; 3],
    /// Compatibility single table, at damping_ratio. Mutually exclusive with tables.
    #[serde(default)]
    pub pseudo_acceleration: Vec<[f64; 2]>,
    #[serde(default)]
    pub period_tables: Vec<SpectrumTable>,
    /// Ascending damping ratios; each table has ascending frequency in Hz.
    #[serde(default)]
    pub tables: Vec<SpectrumTable>,
    /// Empty uses direction at unit scale. Otherwise 1..3 orthogonal directions.
    #[serde(default)]
    pub directions: Vec<SpectrumDirection>,
    #[serde(default)]
    pub direction_combination: DirectionCombination,
}
fn table_valid(p: &[[f64; 2]]) -> bool {
    (2..=100000).contains(&p.len())
        && p.iter().all(|p| p[0].is_finite() && p[0] > 0.0 && p[1].is_finite() && p[1] >= 0.0)
        && p.windows(2).all(|w| w[0][0] < w[1][0] && (w[1][0] - w[0][0]).is_finite())
}
fn interpolate(points: &[[f64; 2]], x: f64) -> Result<f64> {
    let first = points.first().ok_or_else(|| Error::Invalid("empty spectrum table".into()))?;
    let last = points.last().ok_or_else(|| Error::Invalid("empty spectrum table".into()))?;
    if !x.is_finite() || x < first[0] || x > last[0] {
        return Err(Error::Invalid(format!("spectrum axis does not cover {x}")));
    }
    if x == last[0] {
        return Ok(last[1]);
    }
    let j = points.partition_point(|p| p[0] <= x);
    let a = points.get(j.saturating_sub(1)).ok_or_else(|| Error::Invalid("spectrum interval".into()))?;
    let b = points.get(j).ok_or_else(|| Error::Invalid("spectrum interval".into()))?;
    let w = (x - a[0]) / (b[0] - a[0]);
    Ok((1.0 - w) * a[1] + w * b[1])
}
impl Spectrum {
    fn active_tables(&self) -> &[SpectrumTable] {
        if self.period_tables.is_empty() { &self.tables } else { &self.period_tables }
    }

    pub fn mode_damping(&self, mode: usize) -> f64 {
        self.modal_damping_ratios.get(if self.modal_damping_ratios.len() == 1 { 0 } else { mode }).copied().unwrap_or(self.damping_ratio)
    }
    pub fn excitations(&self) -> Vec<SpectrumDirection> {
        if self.directions.is_empty() { vec![SpectrumDirection { direction: self.direction, scale: 1.0 }] } else { self.directions.clone() }
    }
    pub fn validate(&self) -> Result<()> {
        if !self.period_tables.is_empty() && (!self.tables.is_empty() || !self.pseudo_acceleration.is_empty()) {
            return Err(Error::Invalid("choose frequency tables or period_tables, not both".into()));
        }
        let tables = self.active_tables();
        ratios(self.modes, &[self.damping_ratio])?;
        ratios(self.modes, &self.modal_damping_ratios)?;
        direction(self.direction)?;
        if self.directions.len() > 3 {
            return Err(Error::Invalid("spectrum supports at most three directions".into()));
        }
        for (i, d) in self.directions.iter().enumerate() {
            direction(d.direction)?;
            if !d.scale.is_finite() || d.scale < 0.0 {
                return Err(Error::Invalid("direction scale must be finite and nonnegative".into()));
            }
            for prior in self.directions.iter().take(i) {
                let dot: f64 = d.direction.iter().zip(prior.direction).map(|(a, b)| a * b).sum();
                if dot.abs() > 1e-10 {
                    return Err(Error::Invalid("spectrum directions must be orthogonal".into()));
                }
            }
        }
        if !self.directions.is_empty() && self.direction != x_direction() {
            return Err(Error::Invalid("do not supply both direction and directions".into()));
        }
        if tables.is_empty() {
            if !table_valid(&self.pseudo_acceleration) {
                return Err(Error::Invalid("spectrum needs a valid single table or damping tables".into()));
            }
        } else if !self.pseudo_acceleration.is_empty()
            || tables.len() > 100
            || tables.iter().any(|t| !t.damping_ratio.is_finite() || !(0.0..1.0).contains(&t.damping_ratio) || !table_valid(&t.points))
            || tables.windows(2).any(|w| w[0].damping_ratio >= w[1].damping_ratio)
        {
            return Err(Error::Invalid("spectrum damping tables must be ascending, finite and exclusive of the single table (at most 100)".into()));
        }
        for i in 0..self.modes as usize {
            let z = self.mode_damping(i);
            if matches!(self.combination, Combination::Cqc) && z == 0.0 {
                return Err(Error::Invalid("CQC requires positive damping".into()));
            }
            if tables.is_empty() {
                if z != self.damping_ratio {
                    return Err(Error::Invalid("single spectrum table cannot interpolate modal damping".into()));
                }
            } else if tables.first().is_none_or(|t| z < t.damping_ratio) || tables.last().is_none_or(|t| z > t.damping_ratio) {
                return Err(Error::Invalid("spectrum damping tables do not bracket every modal damping ratio".into()));
            }
        }
        Ok(())
    }
    /// Linear in the declared axis, then linear in damping ratio. No extrapolation.
    /// Call validate before modal iteration.
    pub fn ordinate(&self, frequency: f64, damping: f64) -> Result<f64> {
        if !frequency.is_finite() || frequency <= 0.0 || !damping.is_finite() {
            return Err(Error::Invalid("invalid spectrum query".into()));
        }
        let tables = self.active_tables();
        let x = if self.period_tables.is_empty() { frequency } else { 1.0 / frequency };
        if tables.is_empty() {
            if damping != self.damping_ratio {
                return Err(Error::Invalid("single spectrum damping mismatch".into()));
            }
            return interpolate(&self.pseudo_acceleration, x);
        }
        if let Some(t) = tables.iter().find(|t| t.damping_ratio == damping) {
            return interpolate(&t.points, x);
        }
        let i = tables.partition_point(|t| t.damping_ratio < damping);
        if i == 0 {
            return Err(Error::Invalid("damping below spectrum tables".into()));
        }
        let lo = tables.get(i - 1).ok_or_else(|| Error::Invalid("damping interval".into()))?;
        let hi = tables.get(i).ok_or_else(|| Error::Invalid("damping above spectrum tables".into()))?;
        let w = (damping - lo.damping_ratio) / (hi.damping_ratio - lo.damping_ratio);
        let value = (1.0 - w) * interpolate(&lo.points, x)? + w * interpolate(&hi.points, x)?;
        if !value.is_finite() {
            return Err(Error::Invalid("spectrum interpolation overflow".into()));
        }
        Ok(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn spectrum() -> Spectrum {
        Spectrum {
            modes: 1,
            damping_ratio: 0.05,
            modal_damping_ratios: vec![],
            combination: Combination::Cqc,
            direction: [1.0, 0.0, 0.0],
            pseudo_acceleration: vec![],
            period_tables: vec![],
            directions: vec![],
            direction_combination: DirectionCombination::Srss,
            tables: vec![
                SpectrumTable { damping_ratio: 0.02, points: vec![[1.0, 2.0], [3.0, 6.0]] },
                SpectrumTable { damping_ratio: 0.08, points: vec![[1.0, 4.0], [3.0, 12.0]] },
            ],
        }
    }
    #[test]
    fn two_axis_interpolation_and_no_extrapolation() -> Result<()> {
        let s = spectrum();
        s.validate()?;
        assert!((s.ordinate(2.0, 0.05)? - 6.0).abs() < 1e-14);
        assert_eq!(s.ordinate(1.0, 0.02)?, 2.0);
        assert_eq!(s.ordinate(3.0, 0.08)?, 12.0);
        assert!(s.ordinate(0.9, 0.05).is_err());
        assert!(s.ordinate(2.0, 0.1).is_err());
        assert!(s.ordinate(f64::NAN, 0.05).is_err());
        Ok(())
    }
    #[test]
    fn absent_base_displacement_preserves_existing_step_serialization() {
        // Definition hashes serialize the step: a new null field would discard
        // previously saved histories even when the input physics is unchanged.
        let t: Transient = serde_json::from_str(r#"{"time_step":0.001,"increments":1000}"#).unwrap();
        let encoded = serde_json::to_value(&t).unwrap();
        assert!(encoded.get("base_displacement").is_none());
        assert_eq!(serde_json::from_value::<Transient>(encoded).unwrap(), t);
    }
}
