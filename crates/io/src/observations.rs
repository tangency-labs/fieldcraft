//! Measured modal tables and sensor locations. CSV is a deliberately narrow numeric dialect:
//! comma separators, a named header, no quoting, finite numbers, bounded input. Unsupported
//! columns are rejected rather than interpreted as measurement units. No solver is run here.
use std::collections::BTreeSet;

use fieldcraft_model::Model;
use serde::{Deserialize, Serialize};

use crate::{Error, Result};

pub const MAX_BYTES: u64 = 16 << 20;
const MAX_ROWS: usize = 100_000;
const MAX_LINE: usize = 4096;

#[derive(Clone, Debug, Serialize)]
pub struct Mode {
    /// Experimental label, never implicitly a solver result index.
    pub mode: u32,
    pub frequency_hz: f64,
    pub damping_ratio: Option<f64>,
    pub sigma_hz: Option<f64>,
    pub source_line: usize,
}

fn error(line: usize, message: impl std::fmt::Display) -> Error {
    Error::Invalid(format!("observations line {line}: {message}"))
}

fn number(text: &str, line: usize) -> Result<f64> {
    let x = text.parse::<f64>().map_err(|_| error(line, "expected a number"))?;
    if !x.is_finite() {
        return Err(error(line, "number must be finite"));
    }
    Ok(x)
}

/// CSV columns: mode,frequency_hz[,damping_ratio][,sigma_hz], in any order.
/// Mode labels are positive and unique; frequency and sigma are positive. No noise level is
/// invented when absent. UTF-8 BOM and CRLF are allowed; comma-decimal locale is not.
pub fn read_modal_csv(text: &str) -> Result<Vec<Mode>> {
    if text.len() as u64 > MAX_BYTES {
        return Err(Error::Limit("observation file exceeds 16 MiB".into()));
    }
    let mut lines = text.trim_start_matches('\u{feff}').lines();
    let header = lines.next().ok_or_else(|| error(1, "missing header"))?;
    if header.len() > MAX_LINE {
        return Err(error(1, "header too long"));
    }
    let columns: Vec<_> = header.split(',').map(str::trim).collect();
    let mut unique = BTreeSet::new();
    for c in &columns {
        if !matches!(*c, "mode" | "frequency_hz" | "damping_ratio" | "sigma_hz") || !unique.insert(*c) {
            return Err(error(1, format!("unknown or duplicate column `{c}`")));
        }
    }
    let mode_col = columns.iter().position(|c| *c == "mode").ok_or_else(|| error(1, "missing mode"))?;
    let frequency_col = columns.iter().position(|c| *c == "frequency_hz").ok_or_else(|| error(1, "missing frequency_hz"))?;
    let damping_col = columns.iter().position(|c| *c == "damping_ratio");
    let sigma_col = columns.iter().position(|c| *c == "sigma_hz");
    let mut labels = BTreeSet::new();
    let mut result = Vec::new();
    for (i, line) in lines.enumerate() {
        let line_no = i + 2;
        if line.len() > MAX_LINE || result.len() >= MAX_ROWS {
            return Err(error(line_no, "row or file limit exceeded"));
        }
        if line.trim().is_empty() {
            return Err(error(line_no, "empty data row"));
        }
        let fields: Vec<_> = line.split(',').map(str::trim).collect();
        if fields.len() != columns.len() {
            return Err(error(line_no, "column count does not match header"));
        }
        let field = |col| fields.get(col).copied().ok_or_else(|| error(line_no, "missing field"));
        let mode = field(mode_col)?.parse::<u32>().map_err(|_| error(line_no, "mode must be a positive integer"))?;
        if mode == 0 || !labels.insert(mode) {
            return Err(error(line_no, "mode must be positive and unique"));
        }
        let frequency_hz = number(field(frequency_col)?, line_no)?;
        if frequency_hz <= 0.0 {
            return Err(error(line_no, "frequency_hz must be positive"));
        }
        let optional = |col: Option<usize>| -> Result<Option<f64>> { col.map(|c| number(field(c)?, line_no)).transpose() };
        let damping_ratio = optional(damping_col)?;
        if damping_ratio.is_some_and(|d| !(0.0..1.0).contains(&d)) {
            return Err(error(line_no, "damping_ratio must be in [0,1)"));
        }
        let sigma_hz = optional(sigma_col)?;
        if sigma_hz.is_some_and(|s| s <= 0.0) {
            return Err(error(line_no, "sigma_hz must be positive"));
        }
        result.push(Mode { mode, frequency_hz, damping_ratio, sigma_hz, source_line: line_no });
    }
    if result.is_empty() {
        return Err(error(2, "no observations"));
    }
    Ok(result)
}

/// Header of the measured accelerance CSV dialect; the column names carry the unit.
pub const ACCELERANCE_HEADER: &str = "frequency_hz,real_accel_mm_s2_per_N,imag_accel_mm_s2_per_N";

/// One measured complex accelerance sample in mm/s² per N, exactly as written in the file.
#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
pub struct AccelerancePoint {
    pub frequency_hz: f64,
    pub real: f64,
    pub imaginary: f64,
}

/// CSV with the header [`ACCELERANCE_HEADER`] and three finite numbers per row. Frequencies are
/// nonnegative and strictly increasing. Nothing is converted: the Fourier sign convention and the
/// quantity (accelerance, not receptance) are left for the caller to state explicitly.
pub fn read_accelerance_csv(text: &str) -> Result<Vec<AccelerancePoint>> {
    if text.len() as u64 > MAX_BYTES {
        return Err(Error::Limit("observation file exceeds 16 MiB".into()));
    }
    let mut lines = text.trim_start_matches('\u{feff}').lines();
    let header = lines.next().ok_or_else(|| error(1, "missing header"))?;
    if header.trim() != ACCELERANCE_HEADER {
        return Err(error(1, format!("accelerance header must be `{ACCELERANCE_HEADER}`")));
    }
    let mut samples: Vec<AccelerancePoint> = Vec::new();
    for (i, line) in lines.enumerate() {
        let line_no = i + 2;
        if line.len() > MAX_LINE || samples.len() >= MAX_ROWS {
            return Err(error(line_no, "row or file limit exceeded"));
        }
        let fields: Vec<_> = line.split(',').map(str::trim).collect();
        let [frequency, real, imaginary] = fields.as_slice() else {
            return Err(error(line_no, "expected frequency_hz, real and imaginary columns"));
        };
        let point =
            AccelerancePoint { frequency_hz: number(frequency, line_no)?, real: number(real, line_no)?, imaginary: number(imaginary, line_no)? };
        if point.frequency_hz < 0.0 || samples.last().is_some_and(|last| last.frequency_hz >= point.frequency_hz) {
            return Err(error(line_no, "frequency_hz must be nonnegative and strictly increasing"));
        }
        samples.push(point);
    }
    if samples.is_empty() {
        return Err(error(2, "no samples"));
    }
    Ok(samples)
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Sensor {
    pub id: String,
    pub xyz_m: [f64; 3],
    /// Unit vector in the same global Cartesian frame as the model, retaining sign.
    pub direction_global: [f64; 3],
}

#[derive(Debug, Serialize)]
pub struct Mapping {
    pub sensor: String,
    pub node: u32,
    pub distance_m: f64,
    pub direction_global: [f64; 3],
}

/// Nearest-node mapping in SI. Exact ties use the smaller label (model nodes are ordered).
/// Direction is retained as a projection, never rounded to the dominant Cartesian axis.
/// A work cap prevents hostile sensor × node inputs from taking unbounded time.
pub fn map_sensors(model: &Model, sensors: &[Sensor], max_distance_m: f64) -> Result<Vec<Mapping>> {
    if !max_distance_m.is_finite() || max_distance_m < 0.0 {
        return Err(error(0, "max_distance_m must be finite and nonnegative"));
    }
    if sensors.len() > MAX_ROWS || sensors.len().saturating_mul(model.nodes().len()) > 10_000_000 {
        return Err(Error::Limit("sensor mapping exceeds 100000 sensors or 10000000 distance evaluations".into()));
    }
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for sensor in sensors {
        if sensor.id.is_empty() || sensor.id.len() > 128 || sensor.id.chars().any(char::is_control) || !seen.insert(&sensor.id) {
            return Err(error(0, "sensor IDs must be nonempty, unique, at most 128 bytes without control characters"));
        }
        if sensor.xyz_m.iter().chain(sensor.direction_global.iter()).any(|v| !v.is_finite()) {
            return Err(error(0, "sensor position and direction must be finite"));
        }
        let [dx, dy, dz] = sensor.direction_global;
        let norm = dx.hypot(dy).hypot(dz);
        if (norm - 1.0).abs() > 1e-10 {
            return Err(error(0, "sensor direction must be a global unit vector"));
        }
        let [sx, sy, sz] = sensor.xyz_m;
        let mut best: Option<(u32, f64)> = None;
        for (&label, node) in model.nodes() {
            let [x, y, z] = node.xyz;
            let d = (sx - x).hypot(sy - y).hypot(sz - z);
            if best.is_none_or(|(_, old)| d < old) {
                best = Some((label, d));
            }
        }
        let (node, distance_m) = best.ok_or_else(|| error(0, "no model nodes for sensor mapping"))?;
        if !distance_m.is_finite() || distance_m > max_distance_m {
            return Err(error(0, format!("sensor {} nearest node {node} distance {distance_m} m exceeds {max_distance_m} m", sensor.id)));
        }
        out.push(Mapping { sensor: sensor.id.clone(), node, distance_m, direction_global: sensor.direction_global });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn measured_frequency_table_retains_labels_damping_and_missing_sigma() {
        let modes = read_modal_csv("mode,frequency_hz,damping_ratio\r\n1,6.376,0.013\r\n3,33.125,0.0083\r\n").unwrap();
        assert_eq!(modes[1].mode, 3);
        assert_eq!(modes[1].frequency_hz, 33.125);
        assert_eq!(modes[0].damping_ratio, Some(0.013));
        assert_eq!(modes[0].sigma_hz, None);
    }

    #[test]
    fn malformed_and_ambiguous_tables_are_refused() {
        for text in [
            "",
            "mode,hz\n1,3",
            "mode,frequency_hz\n1,NaN",
            "mode,frequency_hz\n1,inf",
            "mode,frequency_hz\n0,3",
            "mode,frequency_hz\n1,3\n1,4",
            "mode,frequency_hz\n1,3,4",
            "mode,frequency_hz,sigma_hz\n1,3,0",
            "mode,frequency_hz,damping_ratio\n1,3,1.3",
            "mode,frequency_hz\n1,3\n\n",
        ] {
            assert!(read_modal_csv(text).is_err(), "{text}");
        }
    }

    #[test]
    fn accelerance_table_keeps_values_and_refuses_damaged_rows() {
        let text = format!("{ACCELERANCE_HEADER}\r\n1,-2,3\r\n2,4,-5\r\n");
        let samples = read_accelerance_csv(&text).unwrap();
        assert_eq!(samples[1], AccelerancePoint { frequency_hz: 2.0, real: 4.0, imaginary: -5.0 });
        for bad in [
            String::new(),
            "frequency_hz,real,imag\n1,2,3".into(),
            format!("{ACCELERANCE_HEADER}\n"),
            format!("{ACCELERANCE_HEADER}\n1,2"),
            format!("{ACCELERANCE_HEADER}\n1,2,NaN"),
            format!("{ACCELERANCE_HEADER}\n2,1,1\n1,1,1"),
            format!("{ACCELERANCE_HEADER}\n-1,1,1"),
            format!("{ACCELERANCE_HEADER}\n1,1,1\n\n"),
        ] {
            assert!(read_accelerance_csv(&bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn mapping_preserves_sign_and_ties_and_rejects_distance() {
        let mut model = Model::default();
        model.set_node(7, [1.0, 0.0, 0.0]).unwrap();
        model.set_node(3, [-1.0, 0.0, 0.0]).unwrap();
        let sensors = [Sensor { id: "tip-z".into(), xyz_m: [0.0; 3], direction_global: [0.0, 0.0, -1.0] }];
        let out = map_sensors(&model, &sensors, 1.0).unwrap();
        assert_eq!(out[0].node, 3);
        assert_eq!(out[0].distance_m, 1.0);
        assert_eq!(out[0].direction_global, [0.0, 0.0, -1.0]);
        assert!(map_sensors(&model, &sensors, 0.9).is_err());
        assert!(map_sensors(&model, &sensors, f64::NAN).is_err());
    }
}
