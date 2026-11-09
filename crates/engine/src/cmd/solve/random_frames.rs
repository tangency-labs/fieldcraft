//! One-sided PSD projection of solved unit-load harmonic phasors.
//! Relative coordinates: M u'' + C u' + K u = -M r a_g;
//! absolute acceleration transfer r - omega² H (Newland, Random Vibrations, 3rd ed.).
use crate::{EngineError, Result, Session, StepResult};
use fieldcraft_solver::{HarmonicResult, RandomInput, StepResult as ImplicitResult};
use serde_json::{Value, json};

pub(super) fn run(s: &Session, p: &Value, step: &str, table: &Value) -> Result<Value> {
    let result = s.result(step).ok_or_else(|| EngineError::invalid("step", "no stored result for this step; run solve.harmonic first"))?;
    let StepResult::Implicit(ImplicitResult::Harmonic(h)) = result.as_ref() else {
        return Err(EngineError::invalid("step", "the stored result is not harmonic"));
    };
    let mut out = evaluate(h, p, step, table)?;
    element_output(&s.model, h, p, &mut out)?;
    Ok(out)
}

pub(super) fn evaluate(h: &HarmonicResult, p: &Value, step: &str, table: &Value) -> Result<Value> {
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
        json!({"step":step,"method":"harmonic_frame_psd","nodes":nodes,"frequencies":h.frames.iter().map(|f|f.frequency).collect::<Vec<_>>(),"input_psd":psd,"node_displacement_psd":up,"node_acceleration_psd":ap,"node_displacement_rms":ut,"node_acceleration_rms":at,"stress_psd":sp,"stress_rms":st,"probes":probes,"acceleration_reference":if direction.is_some(){"absolute_base_motion"}else{"stationary_support"},"displacement_reference":if direction.is_some(){"relative_to_base"}else{"stationary_support"},"units":{"frequency":"Hz","displacement_rms":"m","acceleration_rms":"m/s²","stress_rms":"Pa","psd":"squared SI response unit per Hz"}}),
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

/// Recover each solid contribution before squaring it. In particular, this is
/// not the square of a stress averaged across neighbouring elements.
pub(super) fn element_output(model: &fieldcraft_model::Model, h: &HarmonicResult, p: &Value, out: &mut Value) -> Result<()> {
    let Some(elements) = p.get("elements") else {
        if p.get("element_phasors").is_some() {
            return Err(EngineError::invalid("elements", "element_phasors requires element labels"));
        }
        return Ok(());
    };
    let labels: Vec<u32> = serde_json::from_value(elements.clone()).map_err(|e| EngineError::invalid("elements", e.to_string()))?;
    let phasors: bool = super::opt(p, "element_phasors")?.unwrap_or(false);
    let recovery = fieldcraft_solver::harmonic::ElementStressRecovery::new(model, h, &labels, &super::limits(p)?)?;
    let input: Vec<f64> = serde_json::from_value(out["input_psd"].clone()).map_err(|e| EngineError::invalid("input_psd", e.to_string()))?;
    if input.len() != h.frames.len() {
        return Err(EngineError::invalid("input_psd", "inconsistent frame count"));
    }
    let mut spectra = Vec::new();
    let mut real = Vec::new();
    let mut imaginary = Vec::new();
    for (frame, psd) in h.frames.iter().zip(input) {
        let re = recovery.recover(&frame.real)?;
        let im = recovery.recover(&frame.imaginary)?;
        spectra.push(re.iter().zip(&im).map(|(r, i)| std::array::from_fn::<_, 6, _>(|c| (r[c] * r[c] + i[c] * i[c]) * psd)).collect());
        if phasors {
            real.push(re);
            imaginary.push(im);
        }
    }
    let (_, rms, _) = summaries(h, &spectra)?;
    out["element_stress"] = json!({
        "locations":recovery.locations,"components":["xx","yy","zz","xy","xz","yz"],
        "frequencies":h.frames.iter().map(|f| f.frequency).collect::<Vec<_>>(),
        "components_psd":spectra,"components_rms":rms,
        "stress_real":if phasors { Some(real) } else { None },
        "stress_imaginary":if phasors { Some(imaginary) } else { None },
        "units":{"xyz":"m","stress":"Pa","psd":"Pa²/Hz","rms":"Pa"},
        "recovery":"global Cartesian elastic stresses; integration-rule order (one-based), element-node extrapolation without inter-element averaging; component RMS, not von Mises RMS"
    });
    Ok(())
}
