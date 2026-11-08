//! Real measured-mode observations in a common sensor basis. Reference-DOF
//! normalization differentiates the quotient of Nelson's mass-normalized vectors;
//! its arbitrary modal scale and sign therefore cancel. Nelson (1976),
//! doi:10.2514/3.7211; MAC assignment references are in `pairing`.
use super::{FrequencySensitivity, ModePairing, ModeShapeDerivative, pair_modes};
use crate::{Error, Result};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ModalSensor {
    pub node: u32,
    pub dof: u8,
}
#[derive(Clone, Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ShapeNormalization {
    Reference { reference_sensor: usize },
}
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ModalObservationGroup {
    pub id: String,
    pub step: String,
    pub sensors: Vec<ModalSensor>,
    pub measured_modes: Vec<Vec<f64>>,
    pub min_mac: f64,
    pub min_assignment_margin: f64,
    pub normalization: ShapeNormalization,
}
#[derive(Debug, Serialize)]
pub struct ModalPairingReport {
    pub group: String,
    pub pairing: ModePairing,
}

pub(super) fn sample(v: &ModeShapeDerivative, s: &ModalSensor) -> Result<(f64, f64)> {
    let i = v.nodes.binary_search(&s.node).map_err(|_| Error::Invalid(format!("modal sensor node {} missing", s.node)))?;
    let (a, b, j) = match s.dof {
        1..=3 => (&v.shape_si, &v.derivative_si, usize::from(s.dof - 1)),
        4..=6 => (&v.rotation_si, &v.rotation_derivative_si, usize::from(s.dof - 4)),
        _ => return Err(Error::Invalid("modal sensor dof must be 1..=6".into())),
    };
    Ok((
        *a.get(i).and_then(|x| x.get(j)).ok_or_else(|| Error::Invalid("modal sensor component missing".into()))?,
        *b.get(i).and_then(|x| x.get(j)).ok_or_else(|| Error::Invalid("modal sensor derivative missing".into()))?,
    ))
}

pub(super) fn validate(groups: &[ModalObservationGroup]) -> Result<()> {
    if groups.len() > 64 {
        return Err(Error::Invalid("at most 64 modal observation groups".into()));
    }
    for (i, g) in groups.iter().enumerate() {
        if g.id.is_empty() || groups[..i].iter().any(|x| x.id == g.id) || g.sensors.is_empty() || g.sensors.len() > 4096 {
            return Err(Error::Invalid("modal groups require unique nonempty ids and 1..=4096 sensors".into()));
        }
        let ShapeNormalization::Reference { reference_sensor } = g.normalization;
        if reference_sensor >= g.sensors.len() {
            return Err(Error::Invalid("reference sensor out of range".into()));
        }
        // Do not mix translations and rotations in an unweighted Euclidean MAC.
        let rotational = g.sensors[0].dof > 3;
        for (j, s) in g.sensors.iter().enumerate() {
            if !(1..=6).contains(&s.dof) || (s.dof > 3) != rotational || g.sensors[..j].iter().any(|x| x.node == s.node && x.dof == s.dof) {
                return Err(Error::Invalid("modal sensors must be unique, valid, and share translation/rotation dimensions".into()));
            }
        }
        if g.measured_modes.is_empty() || g.measured_modes.len() > 64 || g.measured_modes.iter().any(|x| x.len() != g.sensors.len()) {
            return Err(Error::Invalid("measured modes require 1..=64 vectors in the sensor basis".into()));
        }
        for row in &g.measured_modes {
            normalized(row, &vec![0.0; row.len()], reference_sensor, reference_sensor)?;
        }
        if !g.min_mac.is_finite() || !(0.0..=1.0).contains(&g.min_mac) || !g.min_assignment_margin.is_finite() || g.min_assignment_margin < 0.0 {
            return Err(Error::Invalid("invalid MAC thresholds".into()));
        }
    }
    Ok(())
}

pub(super) fn pairing(g: &ModalObservationGroup, result: &FrequencySensitivity) -> Result<ModePairing> {
    let candidates = result
        .modes
        .iter()
        .map(|m| {
            let v = m.vector.as_ref().ok_or_else(|| Error::Internal("missing Nelson vector".into()))?;
            g.sensors.iter().map(|s| sample(v, s).map(|x| x.0)).collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    pair_modes(&g.measured_modes, &candidates, g.min_mac, g.min_assignment_margin)
}

pub(super) fn normalized(values: &[f64], derivatives: &[f64], reference: usize, sensor: usize) -> Result<(f64, f64)> {
    if values.len() != derivatives.len() || values.iter().chain(derivatives).any(|x| !x.is_finite()) {
        return Err(Error::Invalid("nonfinite modal normalization inputs".into()));
    }
    let get = |x: &[f64], i| x.get(i).copied().ok_or_else(|| Error::Invalid("modal sensor index out of range".into()));
    let r = get(values, reference)?;
    let dr = get(derivatives, reference)?;
    let scale = values.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
    if scale == 0.0 || r.abs() <= 1e-8 * scale {
        return Err(Error::Invalid("reference modal component too small; choose another sensor".into()));
    }
    let value = get(values, sensor)? / r;
    let derivative = (get(derivatives, sensor)? - value * dr) / r;
    if !value.is_finite() || !derivative.is_finite() {
        return Err(Error::Invalid("modal normalization overflow".into()));
    }
    Ok((value, derivative))
}

pub(super) fn component(
    g: &ModalObservationGroup,
    result: &FrequencySensitivity,
    pairing: &ModePairing,
    measured_mode: usize,
    sensor: usize,
) -> Result<(f64, f64)> {
    let ShapeNormalization::Reference { reference_sensor } = g.normalization;
    if sensor == reference_sensor {
        return Err(Error::Invalid("reference component is constant and supplies no independent information".into()));
    }
    let index = pairing.pairs.get(measured_mode).ok_or_else(|| Error::Invalid("measured mode index out of range".into()))?.computed_mode;
    let v = result.modes.get(index).and_then(|x| x.vector.as_ref()).ok_or_else(|| Error::Internal("paired vector missing".into()))?;
    let samples = g.sensors.iter().map(|s| sample(v, s)).collect::<Result<Vec<_>>>()?;
    let (values, derivatives): (Vec<_>, Vec<_>) = samples.into_iter().unzip();
    normalized(&values, &derivatives, reference_sensor, sensor)
}
