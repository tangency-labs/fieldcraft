//! Direct harmonic points and constant damping from the public Abaqus keyword reference:
//! https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-steadystatedynamics.htm
//! and simakey-r-damping.htm. Frequency ranges with bias are deliberately not approximated.
use super::{Block, canonical, error, number, quantity, row_fields};
use crate::{Error, Result};
use fieldcraft_model::DisplaySystem;
use fieldcraft_model::harmonic::{Damping, Frequencies, Harmonic, MAX_FREQUENCIES, Method};
use fieldcraft_model::units::Quantity;
use std::collections::BTreeMap;

fn unsupported(b: &Block<'_>, what: &str) -> Error {
    Error::Unsupported { format: "Abaqus", what: format!("line {}: {what}", b.line) }
}

pub(super) fn damping(b: &Block<'_>, units: &DisplaySystem) -> Result<Damping> {
    b.only(&["ALPHA", "BETA", "STRUCTURAL"])?;
    if b.params.values().any(Option::is_none) {
        return Err(error(b.line, "damping parameters need numeric values"));
    }
    if !b.rows.is_empty() {
        return Err(unsupported(b, "frequency/temperature dependent damping"));
    }
    let value = |key, q| -> Result<f64> {
        let x = b.value(key).map(|s| number(s, b.line)).transpose()?.unwrap_or(0.0);
        quantity(units, q, x, b.line)
    };
    let result = Damping {
        alpha: value("ALPHA", Quantity::Frequency)?,
        beta: value("BETA", Quantity::Time)?,
        loss_factor: b.value("STRUCTURAL").map(|s| number(s, b.line)).transpose()?.unwrap_or(0.0),
    };
    result.validate().map_err(|e| error(b.line, e.to_string()))?;
    Ok(result)
}

pub(super) fn harmonic(b: &Block<'_>, units: &DisplaySystem) -> Result<Harmonic> {
    if b.params.keys().any(|k| !["DIRECT", "INTERVAL", "FREQUENCYSCALE"].contains(&k.as_str())) {
        return Err(unsupported(b, "additional steady-state dynamics parameters"));
    }
    if ["INTERVAL", "FREQUENCYSCALE"].iter().any(|k| matches!(b.params.get(*k), Some(None))) {
        return Err(error(b.line, "harmonic interval and scale need values"));
    }
    if !matches!(b.params.get("DIRECT"), Some(None)) {
        return Err(unsupported(b, "only DIRECT steady-state dynamics is mapped"));
    }
    if canonical(b.value("INTERVAL").unwrap_or("RANGE")) != "RANGE" {
        return Err(unsupported(b, "harmonic intervals based on eigenfrequencies"));
    }
    if !matches!(canonical(b.value("FREQUENCYSCALE").unwrap_or("LOGARITHMIC")).as_str(), "LINEAR" | "LOGARITHMIC") {
        return Err(error(b.line, "invalid harmonic frequency scale"));
    }
    if b.rows.len() > MAX_FREQUENCIES {
        return Err(error(b.line, "too many harmonic frequency points"));
    }
    let mut points = Vec::new();
    for row in &b.rows {
        let f = row_fields(row.text, row.line)?;
        let low = number(f.first().map(String::as_str).unwrap_or_default(), row.line)?;
        let high = f.get(1).filter(|s| !s.is_empty()).map(|s| number(s, row.line)).transpose()?.unwrap_or(0.0);
        if high != 0.0 {
            return Err(unsupported(b, "harmonic frequency ranges/bias; provide explicit single-frequency rows"));
        }
        // With upper=0 the public syntax ignores the remaining entries on that row.
        points.push(quantity(units, Quantity::Frequency, low, row.line)?);
    }
    let result =
        Harmonic { frequencies: Frequencies::List(points), method: Method::Direct, damping: Damping::default(), load_factors: BTreeMap::new() };
    result.validate().map_err(|e| error(b.line, e.to_string()))?;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::abaqus::read::{Limits, read_analysis_limited};
    use fieldcraft_model::{Procedure, units::Preset};

    const MODEL: &str = "*NODE\n1,0,0,0\n2,1,0,0\n*ELEMENT,TYPE=B31,ELSET=E\n1,1,2\n*MATERIAL,NAME=M\n*ELASTIC\n200e9,0.3\n*DENSITY\n7800\n*DAMPING,ALPHA=2,BETA=0.001,STRUCTURAL=0.02\n*BEAM SECTION,ELSET=E,MATERIAL=M,SECTION=RECT\n0.1,0.1\n0,0,-1\n*BOUNDARY\n1,ENCASTRE\n";

    #[test]
    fn direct_points_and_material_and_global_damping_are_native() {
        let text = format!(
            "{MODEL}*STEP,NAME=FRF\n*STEADY STATE DYNAMICS,DIRECT,FREQUENCY SCALE=LINEAR\n0,0\n10,0\n*GLOBAL DAMPING,ALPHA=3\n*CLOAD,REAL\n2,1,100\n*END STEP\n"
        );
        let got = read_analysis_limited(&text, Preset::MKgS, Limits::default()).unwrap();
        assert!(got.analysis_complete);
        assert_eq!(got.model.material("M").unwrap().damping, Some(Damping { alpha: 2.0, beta: 0.001, loss_factor: 0.02 }));
        let Procedure::Harmonic(h) = &got.model.step("FRF").unwrap().procedure else { panic!("harmonic") };
        assert_eq!(h.frequencies.points().unwrap(), [0.0, 10.0]);
        assert_eq!(h.damping.alpha, 3.0);
        let out = crate::abaqus::write(&got.model, Preset::MKgS).unwrap();
        let again = read_analysis_limited(&out, Preset::MKgS, Limits::default()).unwrap();
        assert!(again.analysis_complete, "{:?}", again.warnings);
        assert_eq!(again.model.step("FRF").unwrap().procedure, got.model.step("FRF").unwrap().procedure);
    }

    #[test]
    fn unknown_ranges_and_imaginary_loads_are_not_silently_solved() {
        for body in ["*STEADY STATE DYNAMICS,DIRECT\n1,10,20\n", "*STEADY STATE DYNAMICS,DIRECT\n1,0\n*CLOAD,IMAGINARY\n2,1,1\n"] {
            let got = read_analysis_limited(&format!("{MODEL}*STEP\n{body}*END STEP\n"), Preset::MKgS, Limits::default()).unwrap();
            assert!(!got.analysis_complete);
            assert!(got.model.steps().is_empty());
        }
    }
}
