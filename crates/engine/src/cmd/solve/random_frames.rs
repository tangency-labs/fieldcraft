//! One-sided PSD projection of solved unit-load harmonic phasors.
//! Relative coordinates: M u'' + C u' + K u = -M r a_g;
//! absolute acceleration transfer r - omega² H (Newland, Random Vibrations, 3rd ed.).
use crate::{EngineError, Result, Session, StepResult};
use fieldcraft_solver::{HarmonicResult, RandomInput, StepResult as ImplicitResult};
use serde_json::{Value, json};

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct StressProjection {
    weights: [f64; 6],
    input_model: InputModel,
    narrow_band: Option<NarrowBandInput>,
    fatigue: Option<NarrowBandInput>,
}

#[derive(Clone, Copy, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
enum InputModel {
    CommonScalarUnitPattern,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct NarrowBandInput {
    #[serde(default)]
    method: fieldcraft_solver::fatigue::spectral::Method,
    curve: fieldcraft_solver::fatigue::Basquin,
    duration_seconds: f64,
    mean_stress_pa: f64,
    correction: fieldcraft_solver::fatigue::life::StressMean,
    stationary_gaussian: bool,
}

pub(super) fn run(s: &Session, p: &Value, step: &str, table: &Value) -> Result<Value> {
    let result = s.result(step).ok_or_else(|| EngineError::invalid("step", "no stored result for this step; run solve.harmonic first"))?;
    let StepResult::Implicit(ImplicitResult::Harmonic(h)) = result.as_ref() else {
        return Err(EngineError::invalid("step", "the stored result is not harmonic"));
    };
    let table: Vec<[f64; 2]> = serde_json::from_value(table.clone()).map_err(|e| EngineError::invalid("input_psd", e.to_string()))?;
    if table.len() < 2
        || table.len() > 10000
        || table.iter().any(|r| r.iter().any(|v| !v.is_finite() || *v < 0.))
        || table.windows(2).any(|w| w[1][0] <= w[0][0])
    {
        return Err(EngineError::invalid("input_psd", "need 2..10000 finite nonnegative [Hz, PSD] rows with strictly increasing frequencies"));
    }
    let direction: Option<[f64; 3]> = p
        .get("base_direction")
        .map(|v| serde_json::from_value(v.clone()).map_err(|e| EngineError::invalid("base_direction", e.to_string())))
        .transpose()?;
    if direction.is_some_and(|v| v.iter().any(|x| !x.is_finite()) || (v.iter().map(|x| x * x).sum::<f64>() - 1.).abs() > 1e-10) {
        return Err(EngineError::invalid(
            "base_direction",
            "expected a finite unit translation direction; stored frames must be from the matching unit inertial load -M r",
        ));
    }
    let nodes: Vec<u32> = match p.get("nodes") {
        Some(v) => serde_json::from_value(v.clone()).map_err(|e| EngineError::invalid("nodes", e.to_string()))?,
        None => h.nodes.clone(),
    };
    if nodes.is_empty() || nodes.len() > h.nodes.len() || nodes.iter().collect::<std::collections::BTreeSet<_>>().len() != nodes.len() {
        return Err(EngineError::invalid("nodes", "expected nonempty unique solved node labels"));
    }
    let indices = nodes
        .iter()
        .map(|n| h.nodes.binary_search(n).map_err(|_| EngineError::invalid("nodes", format!("node {n} has no harmonic result"))))
        .collect::<Result<Vec<_>>>()?;
    if h.frames.len() < 2 || indices.len().checked_mul(h.frames.len()).is_none_or(|n| n > 4_000_000) {
        return Err(EngineError::invalid("nodes", "PSD projection needs at least two frames and at most 4000000 node-frame samples; select nodes"));
    }
    let psd = h
        .frames
        .iter()
        .map(|frame| {
            let f = frame.frequency;
            let first = table.first().ok_or_else(|| EngineError::invalid("input_psd", "empty table"))?;
            let last = table.last().ok_or_else(|| EngineError::invalid("input_psd", "empty table"))?;
            if !f.is_finite() || f < first[0] || f > last[0] {
                return Err(EngineError::invalid("input_psd", format!("table does not cover {f} Hz")));
            }
            let i = table.partition_point(|r| r[0] <= f).saturating_sub(1).min(table.len() - 2);
            let a = table.get(i).ok_or_else(|| EngineError::invalid("input_psd", "invalid interval"))?;
            let b = table.get(i + 1).ok_or_else(|| EngineError::invalid("input_psd", "invalid interval"))?;
            Ok(a[1] + (b[1] - a[1]) * (f - a[0]) / (b[0] - a[0]))
        })
        .collect::<Result<Vec<_>>>()?;
    let projection: Option<StressProjection> = p
        .get("stress_projection")
        .map(|v| serde_json::from_value(v.clone()).map_err(|e| EngineError::invalid("stress_projection", e.to_string())))
        .transpose()?;
    let projected = projection
        .as_ref()
        .map(|config| -> Result<Vec<Value>> {
            let InputModel::CommonScalarUnitPattern = config.input_model;
            if config.narrow_band.is_some() && config.fatigue.is_some() {
                return Err(EngineError::invalid("stress_projection", "use fatigue or legacy narrow_band, not both"));
            }
            if config.narrow_band.as_ref().is_some_and(|f| f.method != fieldcraft_solver::fatigue::spectral::Method::NarrowBand) {
                return Err(EngineError::invalid("narrow_band", "use fatigue for broadband methods"));
            }
            let fatigue_input = config.fatigue.as_ref().or(config.narrow_band.as_ref());
            if fatigue_input.is_some_and(|f| !f.stationary_gaussian) {
                return Err(EngineError::invalid("stress_projection", "narrow-band damage requires explicit stationary Gaussian applicability"));
            }
            nodes
                .iter()
                .map(|&node| -> Result<Value> {
                    use fieldcraft_solver::fatigue::spectral;
                    let samples = spectral::from_harmonic(h, node, &config.weights, &psd)?;
                    let moments = spectral::moments(&samples)?;
                    let fatigue = fatigue_input
                        .map(|f| -> Result<Value> {
                            let r = spectral::damage(&samples, &f.curve, f.duration_seconds, f.mean_stress_pa, f.correction, f.method)?;
                            Ok(json!({"method":f.method, "damage":r.damage,
                    "narrow_band_limit":r.moments_hz[0]>0.0 && 1.0-r.irregularity<=1e-8,
                    "life_seconds":if r.damage>0.0 { let life=f.duration_seconds/r.damage; life.is_finite().then_some(life) } else { None },
                    "zero_crossing_rate_hz":r.zero_crossing_rate, "irregularity":r.irregularity,
                    "duration_seconds":f.duration_seconds, "mean_stress_pa":f.mean_stress_pa,
                    "applicability":"stationary Gaussian scalar linear stress; spectral approximation; bandwidth not certified"}))
                        })
                        .transpose()?;
                    Ok(json!({"node":node,"psd":samples.iter().map(|s| s.psd).collect::<Vec<_>>(),
                "moments_hz":moments,"rms_pa":moments[0].sqrt(),"fatigue":fatigue}))
                })
                .collect()
        })
        .transpose()?;
    let mut u = Vec::new();
    let mut acc = Vec::new();
    let mut stress = Vec::new();
    let with_stress = h.frames.iter().all(|f| f.stress_real.is_some() && f.stress_imaginary.is_some());
    for (frame, input) in h.frames.iter().zip(&psd) {
        let mut ur = Vec::new();
        let mut ar = Vec::new();
        let mut sr = Vec::new();
        let w2 = (2. * std::f64::consts::PI * frame.frequency).powi(2);
        for &k in &indices {
            let re = frame.real.get(k).ok_or_else(|| EngineError::invalid("step", "missing real displacement"))?;
            let im = frame.imaginary.get(k).ok_or_else(|| EngineError::invalid("step", "missing imaginary displacement"))?;
            let up: [f64; 3] = std::array::from_fn(|j| (re[j] * re[j] + im[j] * im[j]) * input);
            let ap: [f64; 3] = std::array::from_fn(|j| ((direction.map_or(0., |d| d[j]) - w2 * re[j]).powi(2) + (w2 * im[j]).powi(2)) * input);
            ur.push(up);
            ar.push(ap);
            if with_stress {
                let r = frame.stress_real.as_ref().and_then(|v| v.get(k)).ok_or_else(|| EngineError::invalid("step", "missing stress"))?;
                let i = frame.stress_imaginary.as_ref().and_then(|v| v.get(k)).ok_or_else(|| EngineError::invalid("step", "missing stress"))?;
                sr.push(std::array::from_fn::<_, 6, _>(|j| (r[j] * r[j] + i[j] * i[j]) * input));
            }
        }
        u.push(ur);
        acc.push(ar);
        if with_stress {
            stress.push(sr);
        }
    }
    let (up, uc, ut) = summaries(h, &u)?;
    let (ap, ac, at) = summaries(h, &acc)?;
    let (sp, sc, st) = if with_stress { summaries(h, &stress)? } else { (Vec::new(), Vec::new(), Vec::new()) };
    let probes = nodes.iter().enumerate().map(|(k,n)| {
        let samples = h.frames.iter().enumerate().map(|(j,f)| json!({"frequency":f.frequency,"displacement_psd":up.get(j).and_then(|r|r.get(k)),"acceleration_psd":ap.get(j).and_then(|r|r.get(k)),"displacement_components_psd":u.get(j).and_then(|r|r.get(k)),"acceleration_components_psd":acc.get(j).and_then(|r|r.get(k)),"stress_components_psd":stress.get(j).and_then(|r|r.get(k))})).collect::<Vec<_>>();
        json!({"node":n,"displacement_rms":ut.get(k),"acceleration_rms":at.get(k),"displacement_components_rms":uc.get(k),"acceleration_components_rms":ac.get(k),"stress_components_rms":sc.get(k),"samples":samples})
    }).collect::<Vec<_>>();
    Ok(
        json!({"step":step,"method":"harmonic_frame_psd","nodes":nodes,"frequencies":h.frames.iter().map(|f|f.frequency).collect::<Vec<_>>(),"input_psd":psd,"node_displacement_psd":up,"node_acceleration_psd":ap,"node_displacement_rms":ut,"node_acceleration_rms":at,"stress_psd":sp,"stress_rms":st,"projected_stress":projected,"projection_weights":projection.as_ref().map(|v|v.weights),"projection_psd_unit":"Pa²/Hz","probes":probes,"acceleration_reference":if direction.is_some(){"absolute_base_motion"}else{"stationary_support"},"displacement_reference":if direction.is_some(){"relative_to_base"}else{"stationary_support"},"units":{"frequency":"Hz","displacement_rms":"m","acceleration_rms":"m/s²","stress_rms":"Pa","psd":"squared SI response unit per Hz"}}),
    )
}

type Summary<const N: usize> = (Vec<Vec<f64>>, Vec<[f64; N]>, Vec<f64>);
fn summaries<const N: usize>(h: &HarmonicResult, field: &[Vec<[f64; N]>]) -> Result<Summary<N>> {
    let count = field.first().map_or(0, Vec::len);
    let totals = field.iter().map(|r| r.iter().map(|v| v.iter().sum()).collect()).collect::<Vec<Vec<f64>>>();
    let mut components = Vec::new();
    let mut rms = Vec::new();
    for k in 0..count {
        let mut value = [0.; N];
        for (c, dest) in value.iter_mut().enumerate() {
            let samples = h
                .frames
                .iter()
                .zip(field)
                .map(|(f, row)| {
                    let v = row.get(k).ok_or_else(|| EngineError::invalid("step", "inconsistent frame dimensions"))?;
                    Ok(RandomInput { frequency_hz: f.frequency, psd: v[c] })
                })
                .collect::<Result<Vec<_>>>()?;
            *dest = fieldcraft_solver::integrate_psd(&samples)?.sqrt();
        }
        rms.push(value.iter().map(|x| x * x).sum::<f64>().sqrt());
        components.push(value);
    }
    Ok((totals, components, rms))
}

#[cfg(test)]
mod fatigue_tests {
    use super::*;
    pub(super) fn stored_session(second: f64) -> Session {
        let frames = [1.0, 3.0]
            .into_iter()
            .map(|frequency| fieldcraft_solver::harmonic::Frame {
                frequency,
                real: vec![[0.0; 3]],
                imaginary: vec![[0.0; 3]],
                rotation_real: vec![],
                rotation_imaginary: vec![],
                points_real: vec![],
                points_imaginary: vec![],
                stress_real: Some(vec![[2.0, second, 0.0, 0.0, 0.0, 0.0]]),
                stress_imaginary: Some(vec![[0.0; 6]]),
                relative_residual: 0.0,
                residual_norm_si: 0.0,
                projected_residual: None,
                input_power: 0.0,
                dissipated_power: 0.0,
            })
            .collect();
        let h = HarmonicResult {
            step: "synthetic".into(),
            nodes: vec![7],
            reference_points: vec![],
            frames,
            natural_frequencies: vec![],
            residual_vectors: 0,
            stats: Default::default(),
            warnings: vec![],
        };
        let mut s = Session::default();
        s.store_result(StepResult::Implicit(ImplicitResult::Harmonic(h)));
        s
    }
    pub(super) fn params() -> Value {
        json!({"step":"synthetic","nodes":[7],"input_psd":[[1,1],[3,1]],
            "stress_projection":{"weights":[0.5,0.5,0,0,0,0],"input_model":"common_scalar_unit_pattern",
                "narrow_band":{"curve":{"amplitude_ref":2,"cycles_ref":100,"exponent":2,"endurance_limit":0},
                "duration_seconds":10,"mean_stress_pa":0,"correction":{"type":"none"},"stationary_gaussian":true}}})
    }
    #[test]
    fn stored_phasors_to_projected_damage_and_phase_cancellation() {
        let p = params();
        let mut s = stored_session(2.0);
        let r = s.execute("solve.random", &p).unwrap();
        assert_eq!(r["projected_stress"][0]["psd"], json!([4.0, 4.0]));
        let d = r["projected_stress"][0]["fatigue"]["damage"].as_f64().unwrap();
        assert!((d / (0.4 * (13.0_f64 / 3.0).sqrt()) - 1.0).abs() < 1e-12);
        let mut opposite = stored_session(-2.0);
        let cancelled = opposite.execute("solve.random", &p).unwrap();
        assert_eq!(r["stress_psd"], cancelled["stress_psd"]);
        assert_eq!(cancelled["projected_stress"][0]["psd"], json!([0.0, 0.0]));
        assert_eq!(cancelled["projected_stress"][0]["fatigue"]["damage"], 0.0);
        // Source fields remain SI when display units change; changing model units
        // may invalidate stored results, so re-store the same SI test fixture.
        let source = s.result("synthetic").unwrap().as_ref().clone();
        s.execute("model.units", &json!({"system":"mm-t-s-N-MPa"})).unwrap();
        s.store_result(source);
        assert_eq!(s.execute("solve.random", &p).unwrap()["projected_stress"], r["projected_stress"]);
    }
    #[test]
    fn unavailable_coherence_and_stress_are_not_silently_substituted() {
        let mut s = stored_session(2.0);
        let mut p = params();
        p["stress_projection"]["input_model"] = json!("independent");
        assert!(s.execute("solve.random", &p).is_err());
        p = params();
        p["stress_projection"]["narrow_band"]["stationary_gaussian"] = json!(false);
        assert!(s.execute("solve.random", &p).is_err());
        p = params();
        p.as_object_mut().unwrap().remove("input_psd");
        assert!(s.execute("solve.random", &p).is_err());
        p = params();
        let mut source = s.result("synthetic").unwrap().as_ref().clone();
        if let StepResult::Implicit(ImplicitResult::Harmonic(h)) = &mut source {
            h.frames[0].stress_imaginary = None;
        }
        s.store_result(source);
        assert!(s.execute("solve.random", &p).is_err());
    }
}

#[cfg(test)]
mod broadband_command_tests {
    use super::*;
    #[test]
    fn broadband_method_dispatch_on_stored_frames() {
        let mut s = super::fatigue_tests::stored_session(2.0);
        for method in ["dirlik", "benasciutti_tovo"] {
            let mut p = super::fatigue_tests::params();
            let mut config = p["stress_projection"]["narrow_band"].take();
            p["stress_projection"].as_object_mut().unwrap().remove("narrow_band");
            config["method"] = json!(method);
            p["stress_projection"]["fatigue"] = config;
            let out = s.execute("solve.random", &p).unwrap();
            assert_eq!(out["projected_stress"][0]["fatigue"]["method"], method);
            assert!(out["projected_stress"][0]["fatigue"]["damage"].as_f64().unwrap() > 0.0);
            p["stress_projection"]["fatigue"]["mean_stress_pa"] = json!(1.0);
            assert!(s.execute("solve.random", &p).is_err());
        }
    }
}
