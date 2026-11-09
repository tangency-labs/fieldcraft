//! Restricted modal histories from the public MODAL DYNAMIC, MODAL DAMPING,
//! SPECTRUM and RESPONSE SPECTRUM keyword references (abaqus.uclouvain.be).
//! The preceding unwindowed basis is compiled into the native re-extraction
//! contract; it is not a reference to a cached commercial modal result.
use super::{Block, Family, State, allow, canonical, error, label, name, number, row_fields, unsupported, value};
use crate::Result;
use fieldcraft_model::{
    DisplaySystem, Model, Procedure,
    dynamics::{Combination, DirectionCombination, Method, Spectrum, SpectrumDirection, Transient},
    units::Quantity,
};
use std::collections::BTreeMap;

const MAX_MODAL_INCREMENTS: f64 = 1_000_000.0; // unit-factors: bounded increment count, not a unit conversion

pub(super) struct Table {
    points: Vec<[f64; 2]>,
    damping: f64,
}
pub(super) fn tables(blocks: &[Block<'_>], units: &DisplaySystem) -> Result<BTreeMap<String, Table>> {
    let mut tables = BTreeMap::new();
    for b in blocks.iter().filter(|b| b.keyword == "SPECTRUM") {
        allow(b, &["NAME", "TYPE"])?;
        if canonical(value(b, "TYPE", "ACCELERATION")?) != "ACCELERATION" {
            return Err(unsupported(b.line, "SPECTRUM TYPE: only ACCELERATION is mapped"));
        }
        let mut points = Vec::new();
        let mut damping = None;
        let mut ordinate = None;
        for row in &b.rows {
            if points.len() >= 100000 {
                return Err(error(row.line, "SPECTRUM row limit"));
            }
            let f = row_fields(row.text, row.line)?;
            let [a, hz, z] = f.as_slice() else {
                return Err(error(row.line, "SPECTRUM requires magnitude,frequency,damping"));
            };
            let a = super::super::quantity(units, Quantity::Acceleration, number(a, row.line)?, row.line)?;
            let hz = super::super::quantity(units, Quantity::Frequency, number(hz, row.line)?, row.line)?;
            let z = number(z, row.line)?;
            if a <= 0.0 || hz <= 0.0 || !(0.0..1.0).contains(&z) || ordinate.is_some_and(|v| v != a) || damping.is_some_and(|v| v != z) {
                return Err(unsupported(
                    row.line,
                    "SPECTRUM requires one damping ratio and a positive constant acceleration table; logarithmic interpolation is not replaced",
                ));
            }
            if points.last().is_some_and(|p: &[f64; 2]| p[0] >= hz) {
                return Err(error(row.line, "SPECTRUM frequencies must ascend"));
            }
            ordinate = Some(a);
            damping = Some(z);
            points.push([hz, a]);
        }
        if points.len() < 2 {
            return Err(error(b.line, "SPECTRUM requires two frequency endpoints"));
        }
        let table = Table { points, damping: damping.ok_or_else(|| error(b.line, "empty SPECTRUM"))? };
        if tables.insert(name(b.required("NAME")?, b.line)?, table).is_some() {
            return Err(error(b.line, "duplicate SPECTRUM name"));
        }
    }
    Ok(tables)
}
fn damping(body: &[Block<'_>], modes: u32) -> Result<Vec<f64>> {
    let mut values = vec![0.0; modes as usize];
    let mut assigned = vec![false; modes as usize];
    for b in body.iter().filter(|b| b.keyword == "MODALDAMPING") {
        allow(b, &["DEFINITION", "VISCOUS"])?;
        if canonical(value(b, "DEFINITION", "MODE NUMBERS")?) != "MODENUMBERS"
            || canonical(value(b, "VISCOUS", "FRACTION OF CRITICAL DAMPING")?) != "FRACTIONOFCRITICALDAMPING"
        {
            return Err(unsupported(b.line, "MODAL DAMPING requires mode-number fractions of critical damping"));
        }
        for row in &b.rows {
            let f = row_fields(row.text, row.line)?;
            let [lo, hi, z] = f.as_slice() else { return Err(error(row.line, "MODAL DAMPING needs first,last,ratio")) };
            let lo = label(lo, row.line)?;
            let hi = label(hi, row.line)?;
            let z = number(z, row.line)?;
            if lo > hi || hi > modes || !(0.0..1.0).contains(&z) {
                return Err(unsupported(row.line, "MODAL DAMPING mode range or ratio outside extracted basis"));
            }
            for i in lo..=hi {
                let slot = assigned.get_mut((i - 1) as usize).ok_or_else(|| error(row.line, "modal index"))?;
                if *slot {
                    return Err(error(row.line, "overlapping MODAL DAMPING"));
                }
                *slot = true;
                if let Some(v) = values.get_mut((i - 1) as usize) {
                    *v = z;
                }
            }
        }
    }
    Ok(values)
}
pub(super) fn procedure(body: &[Block<'_>], model: &Model, tables: &BTreeMap<String, Table>, units: &DisplaySystem) -> Result<(Procedure, f64)> {
    let b = body
        .iter()
        .find(|b| matches!(b.keyword.as_str(), "MODALDYNAMIC" | "RESPONSESPECTRUM"))
        .ok_or_else(|| error(1, "missing projected procedure"))?;
    if body.iter().filter(|b| matches!(b.keyword.as_str(), "MODALDYNAMIC" | "RESPONSESPECTRUM")).count() != 1 {
        return Err(error(b.line, "multiple projected procedures"));
    }
    let source = model.steps().last().ok_or_else(|| unsupported(b.line, "modal procedure requires an immediately preceding FREQUENCY step"))?;
    let Procedure::Frequency { modes, min_hz: None, max_hz: None } = source.procedure else {
        return Err(unsupported(b.line, "modal procedure requires an immediately preceding unwindowed FREQUENCY step"));
    };
    if model.steps().len() != 1 || source.base_state.is_some() {
        return Err(unsupported(b.line, "modal FREQUENCY base state or preceding general history is not mapped"));
    }
    if let Some((id, e)) = model.elements().iter().next() {
        return Err(unsupported(b.line, format!("modal mass mapping for element {id} {}: only explicit MASS/SPRING1", e.kind.name())));
    }
    if !model.ties().is_empty()
        || !model.contacts().is_empty()
        || !model.couplings().is_empty()
        || !model.foundations().is_empty()
        || !model.reference_points().is_empty()
        || !model.initial_fields().is_empty()
    {
        return Err(unsupported(b.line, "modal TIE/CONTACT/COUPLING/FOUNDATION/reference point/INITIAL CONDITIONS combination"));
    }
    if body.iter().any(|b| b.keyword == "BOUNDARY") {
        return Err(unsupported(b.line, "modal BOUNDARY changes from the FREQUENCY basis"));
    }
    let ratios = damping(body, modes)?;
    if b.keyword == "MODALDYNAMIC" {
        allow(b, &["CONTINUE"])?;
        if canonical(value(b, "CONTINUE", "NO")?) != "NO" {
            return Err(unsupported(b.line, "MODAL DYNAMIC CONTINUE=YES"));
        }
        if b.rows.len() != 1 {
            return Err(error(b.line, "MODAL DYNAMIC requires one time row"));
        }
        let row = b.rows.first().ok_or_else(|| error(b.line, "missing time row"))?;
        let f = row_fields(row.text, row.line)?;
        let [dt, period] = f.as_slice() else { return Err(unsupported(row.line, "MODAL DYNAMIC requires increment,period")) };
        let dt = number(dt, row.line)?;
        let period = number(period, row.line)?;
        let n = period / dt;
        if dt <= 0.0 || !n.is_finite() || !(1.0..=MAX_MODAL_INCREMENTS).contains(&n) || (n - n.round()).abs() > 1e-10 * n {
            return Err(unsupported(row.line, "MODAL DYNAMIC period must contain whole fixed increments"));
        }
        let t = Transient {
            time_step: super::super::quantity(units, Quantity::Time, dt, row.line)?,
            increments: n.round() as u32,
            output_every: 1,
            adaptive: None,
            integration: Default::default(),
            method: Method::Modal { modes, damping_ratios: ratios },
            damping: Default::default(),
            load_amplitudes: BTreeMap::new(),
            initial_displacement: BTreeMap::new(),
            initial_velocity: BTreeMap::new(),
            base_acceleration: None,
            base_displacement: None,
        };
        t.validate().map_err(|e| error(b.line, e.to_string()))?;
        return Ok((Procedure::Transient(t), period));
    }
    allow(b, &["SUM", "COMP"])?;
    if canonical(value(b, "SUM", "ABS")?) != "SRSS" {
        return Err(unsupported(b.line, "RESPONSE SPECTRUM SUM: only explicit SRSS is mapped"));
    }
    let comp = match canonical(value(b, "COMP", "SRSS")?).as_str() {
        "SRSS" => DirectionCombination::Srss,
        "R40" => DirectionCombination::Percent1004040,
        _ => return Err(unsupported(b.line, "RESPONSE SPECTRUM COMP: only SRSS or R40")),
    };
    if !(1..=3).contains(&b.rows.len()) {
        return Err(unsupported(b.line, "RESPONSE SPECTRUM requires 1..3 directions"));
    }
    let mut key = None;
    let mut directions = Vec::new();
    for row in &b.rows {
        let f = row_fields(row.text, row.line)?;
        if !(4..=5).contains(&f.len()) {
            return Err(unsupported(
                row.line,
                "RESPONSE SPECTRUM row requires name,global direction cosines,optional scale; rigid response is not mapped",
            ));
        }
        let k = name(f.first().map(String::as_str).unwrap_or_default(), row.line)?;
        if key.as_ref().is_some_and(|v| v != &k) {
            return Err(unsupported(row.line, "RESPONSE SPECTRUM different tables by direction"));
        }
        key = Some(k);
        let get = |i| number(f.get(i).map(String::as_str).unwrap_or_default(), row.line);
        let scale = f.get(4).filter(|s| !s.is_empty()).map_or(Ok(1.0), |s| number(s, row.line))?;
        directions.push(SpectrumDirection { direction: [get(1)?, get(2)?, get(3)?], scale });
    }
    let key = key.ok_or_else(|| error(b.line, "missing spectrum reference"))?;
    let table = tables.get(&key).ok_or_else(|| unsupported(b.line, format!("RESPONSE SPECTRUM unknown spectrum {key}")))?;
    if ratios.iter().any(|v| *v != table.damping) {
        return Err(unsupported(b.line, "MODAL DAMPING does not match SPECTRUM damping table"));
    }
    let s = Spectrum {
        modes,
        damping_ratio: table.damping,
        modal_damping_ratios: Vec::new(),
        combination: Combination::Srss,
        direction: [1.0, 0.0, 0.0],
        pseudo_acceleration: table.points.clone(),
        period_tables: Vec::new(),
        tables: Vec::new(),
        directions,
        direction_combination: comp,
    };
    s.validate().map_err(|e| unsupported(b.line, e.to_string()))?;
    Ok((Procedure::ResponseSpectrum(s), 0.0))
}
pub(super) fn validate_state(state: &State) -> Result<()> {
    for ((family, _, _, _), v) in state {
        if *family != Family::Boundary || v.value != 0.0 || v.amplitude.is_some() {
            return Err(unsupported(v.line, "RESPONSE SPECTRUM requires zero fixed boundaries and no spatial loads"));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::abaqus::read::{Limits, read_analysis_limited};
    use fieldcraft_model::{Procedure, units::Preset};
    const MODEL: &str = "*NODE\n7,0,0,0\n*ELEMENT,TYPE=MASS,ELSET=M\n10,7\n*ELEMENT,TYPE=SPRING1,ELSET=X\n11,7\n*ELEMENT,TYPE=SPRING1,ELSET=Y\n12,7\n*MASS,ELSET=M\n2\n*SPRING,ELSET=X\n1\n200\n*SPRING,ELSET=Y\n2\n800\n*BOUNDARY\n7,3,3\n";
    const BASIS: &str = "*STEP,NAME=Modes,PERTURBATION\n*FREQUENCY\n2\n*END STEP\n";
    const MODAL: &str = "*STEP,NAME=Pulse,PERTURBATION\n*MODAL DYNAMIC,CONTINUE=NO\n.0005,.2\n*MODAL DAMPING\n1,1,.02\n2,2,.05\n*CLOAD,AMPLITUDE=A\n7,1,1\n7,2,2\n*END STEP\n";
    const AMP: &str = "*AMPLITUDE,NAME=A\n0,1,.2,1\n";
    const TABLE: &str = "*SPECTRUM,NAME=S,TYPE=ACCELERATION\n3,.1,.05\n3,100,.05\n";
    const RESPONSE: &str =
        "*STEP,NAME=Peak,PERTURBATION\n*RESPONSE SPECTRUM,SUM=SRSS,COMP=R40\nS,.6,.8,0,1\nS,-.8,.6,0,2\n*MODAL DAMPING\n1,2,.05\n*END STEP\n";
    #[test]
    fn modal_import_binds_all_preceding_modes_and_damping() {
        let r = read_analysis_limited(&format!("{MODEL}{AMP}{BASIS}{MODAL}"), Preset::MKgS, Limits::default()).unwrap();
        assert!(r.analysis_complete, "{:?}", r.warnings);
        let Procedure::Transient(t) = &r.model.steps()[1].procedure else { panic!("not transient") };
        assert_eq!(t.increments, 400);
        assert_eq!(t.method, fieldcraft_model::dynamics::Method::Modal { modes: 2, damping_ratios: vec![0.02, 0.05] });
        assert_eq!(t.load_amplitudes.len(), 2);
        for source in [
            format!("{MODEL}{AMP}{MODAL}"),
            format!("{MODEL}{AMP}{BASIS}{MODAL}").replace("CONTINUE=NO", "CONTINUE=YES"),
            format!("{MODEL}{AMP}{BASIS}{MODAL}").replace("*MODAL DYNAMIC", "*SELECT EIGENMODES\n1\n*MODAL DYNAMIC"),
            format!("{MODEL}{AMP}{BASIS}{MODAL}").replace("2,2,.05", "3,3,.05"),
        ] {
            let r = read_analysis_limited(&source, Preset::MKgS, Limits::default()).unwrap();
            assert!(!r.analysis_complete);
            assert!(r.model.steps().len() <= 1);
        }
    }
    #[test]
    fn competing_random_and_projected_procedures_cannot_overwrite_each_other() {
        let input = format!("{MODEL}{AMP}{BASIS}{MODAL}").replace("*MODAL DYNAMIC", "*RANDOM RESPONSE\n1,2\n*MODAL DYNAMIC");
        let r = read_analysis_limited(&input, Preset::MKgS, Limits::default()).unwrap();
        assert!(!r.analysis_complete);
        assert_eq!(r.model.steps().len(), 1);
        assert!(r.warnings.iter().any(|w| w.contains("multiple procedures")));
    }
    #[test]
    fn spectrum_units_directions_and_named_refusals() {
        let input = format!("{MODEL}{TABLE}{BASIS}{RESPONSE}");
        let r = read_analysis_limited(&input, Preset::MKgS, Limits::default()).unwrap();
        assert!(r.analysis_complete, "{:?}", r.warnings);
        let Procedure::ResponseSpectrum(s) = &r.model.steps()[1].procedure else { panic!("not spectrum") };
        assert_eq!(s.modes, 2);
        assert_eq!(s.directions[1].scale, 2.0);
        assert_eq!(s.pseudo_acceleration, vec![[0.1, 3.0], [100.0, 3.0]]);
        for (a, b) in [
            ("SUM=SRSS", "SUM=ABS"),
            ("COMP=R40", "COMP=ALGEBRAIC"),
            ("TYPE=ACCELERATION", "TYPE=G"),
            ("3,100,.05", "4,100,.05"),
            ("S,.6", "Missing,.6"),
            ("1,2,.05", "1,2,.02"),
            ("S,-.8,.6", "S,.6,.8"),
        ] {
            let r = read_analysis_limited(&input.replace(a, b), Preset::MKgS, Limits::default()).unwrap();
            assert!(!r.analysis_complete, "accepted {a} -> {b}");
        }
    }
}
