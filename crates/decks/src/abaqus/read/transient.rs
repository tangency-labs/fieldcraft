//! Conservative DYNAMIC/DIRECT mapping from the public keyword reference:
//! https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-dynamic.htm
//! Omitted HHT parameters are not assumed to mean Newmark average acceleration.
use super::{Block, Curve, Family, State, allow, canonical, error, label, number, row_fields, unsupported, value};
use crate::Result;
use fieldcraft_model::{DisplaySystem, Model, Procedure, Step, amplitude::Amplitude, units::Quantity};
use std::collections::BTreeMap;

pub(super) fn procedure(b: &Block<'_>, step: &Block<'_>, model: &Model, units: &DisplaySystem) -> Result<(Procedure, f64)> {
    allow(b, &["DIRECT", "APPLICATION", "ALPHA", "BETA", "GAMMA"])?;
    if !matches!(b.params.get("DIRECT"), Some(None)) {
        return Err(unsupported(b.line, "DYNAMIC requires DIRECT without NO STOP; adaptive increments are not mapped"));
    }
    if canonical(value(b, "APPLICATION", "TRANSIENTFIDELITY")?) != "TRANSIENTFIDELITY" {
        return Err(unsupported(b.line, "DYNAMIC APPLICATION other than TRANSIENT FIDELITY"));
    }
    for (key, expected) in [("ALPHA", 0.0), ("BETA", 0.25), ("GAMMA", 0.5)] {
        let Some(Some(text)) = b.params.get(key) else {
            return Err(unsupported(b.line, format!("DYNAMIC requires explicit {key}; numerical damping defaults are not guessed")));
        };
        if number(text, b.line)? != expected {
            return Err(unsupported(b.line, "DYNAMIC currently maps only ALPHA=0, BETA=.25, GAMMA=.5 (average acceleration)"));
        }
    }
    if !model.steps().is_empty() || step.params.contains_key("PERTURBATION") {
        return Err(unsupported(b.line, "DYNAMIC base state or continuation from preceding steps is not mapped"));
    }
    if let Some((id, element)) = model.elements().iter().next() {
        return Err(unsupported(
            b.line,
            format!("DYNAMIC mass policy for element {id} {} is not mapped; current scope is explicit MASS/SPRING1", element.kind.name()),
        ));
    }
    if !model.ties().is_empty()
        || !model.foundations().is_empty()
        || !model.contacts().is_empty()
        || !model.couplings().is_empty()
        || !model.reference_points().is_empty()
        || !model.initial_fields().is_empty()
    {
        return Err(unsupported(b.line, "DYNAMIC interaction/reference-point/initial-field combination"));
    }
    if b.rows.len() != 1 {
        return Err(error(b.line, "DYNAMIC requires one increment/period row"));
    }
    let row = b.rows.first().ok_or_else(|| error(b.line, "missing DYNAMIC row"))?;
    let fields = row_fields(row.text, row.line)?;
    if fields.len() != 2 {
        return Err(unsupported(row.line, "DYNAMIC requires exactly time increment and period"));
    }
    let dt = number(fields.first().map(String::as_str).unwrap_or_default(), row.line)?;
    let period = number(fields.get(1).map(String::as_str).unwrap_or_default(), row.line)?;
    let count = period / dt;
    if dt <= 0.0 || period <= 0.0 || !count.is_finite() || !(1.0..=1_000_000.0).contains(&count) || (count - count.round()).abs() > 1e-10 * count {
        return Err(unsupported(row.line, "DYNAMIC period must be 1..1000000 whole fixed increments"));
    }
    let increments = count.round() as u32;
    if increments > label(value(step, "INC", "100")?, step.line)? {
        return Err(unsupported(step.line, "DYNAMIC requested duration exceeds STEP INC; no silent shortening"));
    }
    let t = fieldcraft_model::dynamics::Transient {
        time_step: super::super::quantity(units, Quantity::Time, dt, row.line)?,
        increments,
        output_every: 1,
        adaptive: None,
        integration: Default::default(),
        method: Default::default(),
        damping: Default::default(),
        load_amplitudes: BTreeMap::new(),
        initial_displacement: BTreeMap::new(),
        initial_velocity: BTreeMap::new(),
        base_acceleration: None,
        base_displacement: None,
    };
    t.validate().map_err(|e| error(b.line, e.to_string()))?;
    Ok((Procedure::Transient(t), period))
}

fn common(state: &State, line: usize) -> Result<&str> {
    let mut name: Option<&str> = None;
    for ((family, _, _, dof), driven) in state {
        if *family == Family::Boundary {
            if driven.value != 0.0 || driven.amplitude.is_some() {
                return Err(unsupported(driven.line, "DYNAMIC nonzero or amplitude-driven BOUNDARY is not mapped"));
            }
            continue;
        }
        if *family != Family::Force || *dof > 3 {
            return Err(unsupported(driven.line, "DYNAMIC supports translational CLOAD only"));
        }
        let curve = driven.amplitude.as_deref().ok_or_else(|| unsupported(driven.line, "DYNAMIC requires an explicit load AMPLITUDE"))?;
        if name.is_some_and(|n| n != curve) {
            return Err(unsupported(driven.line, "DYNAMIC separate load amplitudes are not mapped"));
        }
        name = Some(curve);
    }
    name.ok_or_else(|| unsupported(line, "DYNAMIC requires a driven CLOAD"))
}

pub(super) fn validate_history(b: &Block<'_>, state: &State, curves: &BTreeMap<String, Curve>) -> Result<()> {
    let name = common(state, b.line)?;
    let curve = curves.get(name).ok_or_else(|| error(b.line, "missing DYNAMIC amplitude"))?;
    if curve.points.len() < 2 {
        return Err(unsupported(b.line, "DYNAMIC amplitude requires at least two samples"));
    }
    Ok(())
}

pub(super) fn attach(step: &mut Step, state: &State, curves: &BTreeMap<String, Curve>, units: &DisplaySystem) -> Result<()> {
    let Procedure::Transient(t) = &mut step.procedure else { return Ok(()) };
    let name = common(state, 1)?;
    let curve = curves.get(name).ok_or_else(|| error(1, "missing DYNAMIC amplitude"))?;
    // This independent first step starts at total time zero, so both time scopes coincide.
    let points = curve
        .points
        .iter()
        .map(|(time, value)| Ok([super::super::quantity(units, Quantity::Time, *time, 1)?, *value]))
        .collect::<Result<Vec<_>>>()?;
    let amplitude = if curve.smooth { Amplitude::SmoothStep { points } } else { Amplitude::Tabular { points } };
    for name in &step.loads {
        t.load_amplitudes.insert(name.clone(), amplitude.clone());
    }
    t.validate().map_err(|e| error(1, e.to_string()))
}

#[cfg(test)]
mod tests {
    use crate::abaqus::read::{Limits, read_analysis_limited};
    use fieldcraft_model::{Procedure, units::Preset};
    const INPUT: &str = "*NODE\n7,0,0,0\n*ELEMENT,TYPE=MASS,ELSET=M\n10,7\n*ELEMENT,TYPE=SPRING1,ELSET=K\n11,7\n*MASS,ELSET=M\n2\n*SPRING,ELSET=K\n1\n200\n*BOUNDARY\n7,2,3\n*AMPLITUDE,NAME=Pulse,TIME=TOTAL TIME\n0,0,.1,1,.2,0,.4,0\n*STEP,NAME=Pulse,NLGEOM=NO,INC=400\n*DYNAMIC,DIRECT,APPLICATION=TRANSIENT FIDELITY,ALPHA=0,BETA=.25,GAMMA=.5\n.001,.4\n*CLOAD,AMPLITUDE=Pulse\n7,1,1\n*END STEP\n";
    #[test]
    fn direct_transient_preserves_curve_time_and_reference_load() {
        let r = read_analysis_limited(INPUT, Preset::MKgS, Limits::default()).unwrap();
        assert!(r.analysis_complete, "{:?}", r.warnings);
        let step = &r.model.steps()[0];
        let Procedure::Transient(t) = &step.procedure else { panic!("wrong procedure") };
        assert_eq!(t.increments, 400);
        assert_eq!(t.time_step, 0.001);
        assert_eq!(t.load_amplitudes.len(), 1);
        let fieldcraft_model::Load::Force { vector, .. } = r.model.load(step.loads.iter().next().unwrap()).unwrap() else { panic!("wrong load") };
        assert_eq!(*vector, [1.0, 0.0, 0.0]);
        let step_time = read_analysis_limited(&INPUT.replace("TOTAL TIME", "STEP TIME"), Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(step.procedure, step_time.model.steps()[0].procedure);
        let text = crate::abaqus::write(&r.model, Preset::MKgS).unwrap();
        let back = read_analysis_limited(&text, Preset::MKgS, Limits::default()).unwrap();
        assert!(back.analysis_complete, "{:?}", back.warnings);
        assert!(text.contains("*DYNAMIC"));
    }
    #[test]
    fn direct_transient_refuses_unmapped_numerics_and_history() {
        for (from, to) in [
            ("ALPHA=0,", ""),
            ("ALPHA=0", "ALPHA=-.05"),
            ("DYNAMIC,DIRECT,", "DYNAMIC,DIRECT=NO STOP,"),
            (".001,.4", ".003,.4"),
            ("INC=400", "INC=100"),
            ("*CLOAD,AMPLITUDE=Pulse", "*BOUNDARY\n7,2,3,.1\n*CLOAD,AMPLITUDE=Pulse"),
            ("CLOAD,AMPLITUDE=Pulse", "CLOAD"),
            ("*DYNAMIC,DIRECT", "*MODAL DYNAMIC,DIRECT"),
        ] {
            let source = INPUT.replace(from, to);
            let r = read_analysis_limited(&source, Preset::MKgS, Limits::default()).unwrap();
            assert!(!r.analysis_complete, "accepted {from} -> {to}");
            assert!(r.model.steps().is_empty());
        }
    }
}
