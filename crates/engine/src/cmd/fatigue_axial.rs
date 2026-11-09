//! Explicit nominal axial channel, not a nodal tensor or native result type.
use crate::{EngineError, Result};
use fieldcraft_solver::fatigue::{
    MAX_HISTORY_POINTS,
    curve::{SnCurve, StressLife},
    life::StressMean,
    miner, rainflow,
};
use serde::Deserialize;
use serde_json::{Value, json};

const MEASURE: &str = "nominal_axial_reference_area";
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Input {
    axial_history: History,
    source_id: String,
    step: String,
    max_time_step_seconds: f64,
    curve: SnCurve,
    correction: StressMean,
    curve_stress_measure: String,
    calibration_source: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct History {
    step: String,
    element: u32,
    location: String,
    stress_measure: String,
    axis: String,
    reference_material_axis: [f64; 3],
    reference_area_m2: f64,
    stress_unit: String,
    time_unit: String,
    samples: Vec<Sample>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Sample {
    time: f64,
    stress: f64,
}
fn bad(message: &str) -> EngineError {
    EngineError::invalid("axial_history", message)
}
pub(super) fn run(p: &Value) -> Result<Value> {
    let p: Input = serde_json::from_value(p.clone()).map_err(|e| bad(&e.to_string()))?;
    let h = &p.axial_history;
    if h.stress_measure != MEASURE || p.curve_stress_measure != MEASURE {
        return Err(bad(
            "only nominal_axial_reference_area with a matching S-N and correction-strength calibration is supported; PK2 and Cauchy are not interchangeable",
        ));
    }
    if h.axis != "material_following_axial" || h.location != "constant_rod2_midpoint" {
        return Err(bad(
            "only a material_following_axial constant_rod2_midpoint scalar is supported; no global tensor or nonproportional equivalent stress",
        ));
    }
    let norm: f64 = h.reference_material_axis.iter().map(|x| x * x).sum();
    if !norm.is_finite() || (norm - 1.).abs() > 1e-10 || !h.reference_area_m2.is_finite() || h.reference_area_m2 <= 0. {
        return Err(bad("need a unit reference material axis and positive reference area"));
    }
    if h.element == 0 || p.step.is_empty() || p.step != h.step {
        return Err(bad("need a positive source element and matching nonempty step"));
    }
    for value in [&p.source_id, &p.calibration_source] {
        if value.trim().is_empty() || value.len() > 1024 {
            return Err(bad("source_id and calibration_source must identify the raw source and nominal-stress calibration (1..1024 bytes)"));
        }
    }
    if h.stress_unit != "Pa" || h.time_unit != "s" {
        return Err(bad("nominal stress and time must explicitly use Pa and s"));
    }
    if !p.max_time_step_seconds.is_finite() || p.max_time_step_seconds <= 0. || !(2..=MAX_HISTORY_POINTS).contains(&h.samples.len()) {
        return Err(bad("need 2..1000000 samples and a positive finite maximum sample gap"));
    }
    for (i, sample) in h.samples.iter().enumerate() {
        if !sample.time.is_finite() || !sample.stress.is_finite() {
            return Err(bad("nominal samples must be finite"));
        }
        if i > 0 {
            let dt = sample.time - h.samples[i - 1].time;
            if !dt.is_finite() || dt <= 0. || dt > p.max_time_step_seconds * (1. + 1e-10) {
                return Err(bad("sample times must increase within the declared maximum gap"));
            }
        }
    }
    // The length guard above ensures both endpoints exist without indexing panics.
    let duration = h.samples.last().zip(h.samples.first()).map(|(b, a)| b.time - a.time).ok_or_else(|| bad("missing time span"))?;
    if !duration.is_finite() {
        return Err(bad("duration overflow"));
    }
    let samples: Vec<_> = h.samples.iter().map(|s| s.stress).collect();
    p.curve.validate()?;
    p.correction.corrected_amplitude(0., 0.)?;
    let cycles = rainflow(&samples)?;
    let damage = miner(&cycles, &p.curve, p.correction)?;
    let mut cumulative = 0.;
    let rows = cycles.iter().map(|c| -> Result<Value> {
        let increment = miner(std::slice::from_ref(c), &p.curve, p.correction)?;
        cumulative += increment;
        Ok(json!({"range_pa":c.range,"amplitude_pa":c.range/2.,"mean_pa":c.mean,"count":c.count,"damage":increment,"cumulative_damage":cumulative}))
    }).collect::<Result<Vec<_>>>()?;
    let life = (damage > 0.).then(|| duration / damage).filter(|x| x.is_finite());
    Ok(json!({"source_id":p.source_id,"source_step":p.step,"element":h.element,"location":h.location,
        "stress_measure":h.stress_measure,"curve_stress_measure":p.curve_stress_measure,"calibration_source":p.calibration_source,
        "axis":h.axis,"reference_material_axis":h.reference_material_axis,"reference_area_m2":h.reference_area_m2,
        "units":{"stress":"Pa","time":"s"},"samples":samples.len(),"cycles":rows,"cycle_count":cycles.iter().map(|c| c.count).sum::<f64>(),
        "damage":damage,"duration_seconds":duration,"repeated_record_life_seconds":life,"persisted":false,
        "method":"nominal_axial_rainflow_miner","limitations":"caller-verified nominal-stress calibration including mean-correction strengths; no PK2/Cauchy conversion, nodal tensor, general nonproportional criterion or native persistence; finite-record residual half cycles; no periodic closure or sampling certification"}))
}
