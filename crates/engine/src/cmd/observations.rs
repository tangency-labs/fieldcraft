//! Read-only import at the units/solver boundary. Experimental labels and result indices
//! are different identities; the caller must state their correspondence explicitly.
use std::collections::BTreeSet;

use fieldcraft_io::{observations, roots::Access};
use serde_json::{Value, json};

use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const IMPORT: CommandSpec = CommandSpec {
    id: "observations.import",
    label: "Import Measured Modal Frequencies",
    params: r#"{"path":"modes.csv","step":"Modes","mode_indices":[6,7,8],"sigma_hz":0.1} — numeric CSV columns mode,frequency_hz and optional damping_ratio,sigma_hz. Returns observations for solve.inverse.joint plus original labels, source SHA256, and optional sensor mappings. mode_indices explicitly maps rows to zero-based solver indices; never assumes a rigid-body offset. Positive sigma_hz is required in the file or as fallback parameter; no measurement uncertainty is invented. Optional sensors:[{id,xyz_m:[x,y,z],direction_global:[dx,dy,dz]}] require max_distance_m; SI coordinates and signed global unit directions. No model mutation or solve. Optional format=uff reads ASCII 15/55/58 and imports normal-mode frequencies; raw_uff retains nodes, complex shapes and FRFs in file units without converting them to solver responses. Binary UFF is refused."#,
    form: Some(super::forms::OBSERVATIONS_IMPORT),
    kind: CommandKind::Session,
    requires: &[],
    run: import,
};

pub const READ: CommandSpec = CommandSpec {
    id: "observations.read",
    label: "Read ASCII UFF Measurements",
    params: r#"{"path":"measurements.uff"} — read ASCII UFF 15 nodes, 55 normal-mode vectors and 58 complex FRFs. Returns raw data, skipped datasets and source SHA256. Values retain file units and coordinate IDs; no unit, coordinate, Fourier-sign or modal-scale convention is inferred. Does not map sensors or create solver observations. Unsupported variants and damaged files are refused."#,
    form: Some(super::forms::OBSERVATIONS_READ),
    kind: CommandKind::Session,
    requires: &[],
    run: read,
};

/// Measured accelerance from the narrow CSV dialect of `fieldcraft_io::observations`.
fn accelerance(text: &str) -> Result<Value> {
    let samples = observations::read_accelerance_csv(text)?;
    Ok(json!({"quantity": "accelerance", "unit": "mm/s^2/N", "samples": samples}))
}

fn read(s: &mut Session, p: &Value) -> Result<Value> {
    let path = super::file_param(s, p, "path", Access::Read)?;
    let bytes = super::file::read_file(&path, fieldcraft_io::uff::MAX_BYTES)?;
    let text = std::str::from_utf8(&bytes).map_err(|_| EngineError::Params("UFF must be ASCII text".into()))?;
    let format = super::opt::<String>(p, "format")?.unwrap_or_else(|| "uff".into());
    if format == "accelerance_csv" {
        return Ok(json!({"frequency_response": accelerance(text)?, "source": {"path": path, "sha256": fieldcraft_io::sha256_hex(&bytes)},
            "warnings": ["Accelerance and its complex parts are kept as measured; no displacement conversion or solver observation is made."]}));
    }
    if format != "uff" {
        return Err(EngineError::Params("format must be uff or accelerance_csv".into()));
    }
    let data = fieldcraft_io::uff::read(text)?;
    Ok(json!({"data":data,"source":{"path":path,"sha256":fieldcraft_io::sha256_hex(&bytes)},
        "warnings":["File units, coordinate systems, Fourier convention and modal normalization require explicit interpretation before solver use."]}))
}

fn import(s: &mut Session, p: &Value) -> Result<Value> {
    let path = super::file_param(s, p, "path", Access::Read)?;
    let step: String = super::req(p, "step")?;
    if step.is_empty() || step.len() > 256 || step.chars().any(char::is_control) {
        return Err(EngineError::Params("step must be nonempty, at most 256 bytes without control characters".into()));
    }
    let bytes = super::file::read_file(&path, observations::MAX_BYTES)?;
    let text = std::str::from_utf8(&bytes).map_err(|_| EngineError::Params("modal CSV must be UTF-8".into()))?;
    let format = super::opt::<String>(p, "format")?.unwrap_or_else(|| "csv".into());
    if format == "accelerance_csv" {
        return Ok(json!({"frequency_response": accelerance(text)?, "observations": [], "sensor_mappings": [],
            "source": {"path": path, "sha256": fieldcraft_io::sha256_hex(&bytes)},
            "warnings": ["The quantity is accelerance; frequency-response observations need displacement, so none were generated."]}));
    }
    let raw_uff = if format == "uff" {
        Some(fieldcraft_io::uff::read(text)?)
    } else if format == "csv" {
        None
    } else {
        return Err(EngineError::Params("format must be csv or uff".into()));
    };
    let modes = if let Some(raw) = &raw_uff {
        if raw.modes.is_empty() {
            return Err(EngineError::Params("UFF contains no supported modal frequency observations".into()));
        }
        raw.modes
            .iter()
            .map(|m| observations::Mode {
                mode: m.label,
                frequency_hz: m.frequency_hz,
                damping_ratio: Some(m.viscous_damping_ratio),
                sigma_hz: None,
                source_line: m.source_line,
            })
            .collect()
    } else {
        observations::read_modal_csv(text)?
    };
    let indices: Vec<usize> = super::req(p, "mode_indices")?;
    if indices.len() != modes.len()
        || indices.iter().copied().collect::<BTreeSet<_>>().len() != indices.len()
        || indices.iter().any(|i| *i >= 100_000)
    {
        return Err(EngineError::Params("mode_indices must give one unique zero-based index below 100000 per measured row".into()));
    }
    let sigma: Option<f64> = super::opt(p, "sigma_hz")?;
    if sigma.is_some_and(|v| !v.is_finite() || v <= 0.0) {
        return Err(EngineError::Params("sigma_hz must be finite and positive".into()));
    }
    let mut output = Vec::new();
    let mut sigma_sources = Vec::new();
    for (mode, index) in modes.iter().zip(indices) {
        let sd = mode.sigma_hz.or(sigma).ok_or_else(|| {
            EngineError::Params(format!("mode {} has no sigma_hz; supply measurement uncertainty or an explicit weighting assumption", mode.mode))
        })?;
        output.push(json!({"response":{"kind":"frequency","step":step,"mode":index},"value_si":mode.frequency_hz,"standard_deviation_si":sd}));
        sigma_sources.push(if mode.sigma_hz.is_some() { "file" } else { "parameter" });
    }
    let sensors: Vec<observations::Sensor> = super::opt(p, "sensors")?.unwrap_or_default();
    let mappings = if sensors.is_empty() { Vec::new() } else { observations::map_sensors(&s.model, &sensors, super::req(p, "max_distance_m")?)? };
    Ok(json!({"observations":output,"measured_modes":modes,"sigma_sources":sigma_sources,
        "raw_uff":raw_uff,"sensor_mappings":mappings,"source":{"path":path,"sha256":fieldcraft_io::sha256_hex(&bytes)},
        "warnings":["Frequency pairing is caller-supplied, not MAC-verified.","Damping is retained as metadata, not fitted by this import."]}))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accelerance_csv_preserves_quantity_complex_parts_and_hash_without_displacement_conversion() {
        let path = std::env::temp_dir().join(format!("fieldcraft-accel-{}.csv", std::process::id()));
        std::fs::write(&path, "frequency_hz,real_accel_mm_s2_per_N,imag_accel_mm_s2_per_N\n1,-2,3\n2,4,-5\n").unwrap();
        let mut session = Session::default();
        let result = import(&mut session, &json!({"path":path,"step":"H","format":"accelerance_csv"})).unwrap();
        assert_eq!(result["frequency_response"]["quantity"], "accelerance");
        assert_eq!(result["frequency_response"]["unit"], "mm/s^2/N");
        assert_eq!(result["frequency_response"]["samples"][0]["imaginary"], 3.0);
        assert!(result["observations"].as_array().unwrap().is_empty());
        assert_eq!(result["source"]["sha256"].as_str().unwrap().len(), 64);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn read_accelerance_csv_preserves_quantity_and_complex_parts() {
        let path = std::env::temp_dir().join(format!("fieldcraft-read-accel-{}.csv", std::process::id()));
        std::fs::write(&path, "frequency_hz,real_accel_mm_s2_per_N,imag_accel_mm_s2_per_N\n1,-2,3\n").unwrap();
        let mut session = Session::default();
        let result = read(&mut session, &json!({"path":path,"format":"accelerance_csv"})).unwrap();
        assert_eq!(result["frequency_response"]["quantity"], "accelerance");
        assert_eq!(result["frequency_response"]["samples"][0]["imaginary"], 3.0);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn measured_table_generates_joint_observations_without_guessing_sigma_or_pairing() {
        let path = std::env::temp_dir().join(format!("fieldcraft-observation-{}.csv", std::process::id()));
        std::fs::write(&path, "mode,frequency_hz,damping_ratio\n1,6.376,0.013\n2,16.104,0.013\n").unwrap();
        let mut session = Session::default();
        let mut params = json!({"path":path,"step":"Modes","mode_indices":[6,7]});
        assert!(import(&mut session, &params).is_err());
        params["sigma_hz"] = json!(0.1);
        let result = import(&mut session, &params).unwrap();
        assert_eq!(result["observations"][0]["response"]["mode"], 6);
        assert_eq!(result["observations"][0]["value_si"], 6.376);
        assert_eq!(result["sigma_sources"][0], "parameter");
        let _: Vec<fieldcraft_solver::sensitivity::JointObservation> = serde_json::from_value(result["observations"].clone()).unwrap();
        assert!(session.model.nodes().is_empty());
        params["mode_indices"] = json!([0, 0]);
        assert!(import(&mut session, &params).is_err());
        std::fs::remove_file(path).unwrap();
    }
}
