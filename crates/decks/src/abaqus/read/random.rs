//! Bounded random history mapping from the public *RANDOM RESPONSE, *PSD-DEFINITION,
//! *CORRELATION and *MODAL DAMPING references (abaqus.uclouvain.be/English/SIMACAEKEYRefMap/).
//! Initial scope: one fully correlated load case, real constant multiplier PSD,
//! all previously extracted modes, mode-number damping, no prestress or base motion.
use super::history::{allow, name, unsupported, value};
use super::{Block, canonical, error, label, number, quantity, row_fields};
use crate::Result;
use fieldcraft_model::random::{Excitation, ModalBand, ModalRandom, RandomVibration, SpectrumPoint};
use fieldcraft_model::units::Quantity;
use fieldcraft_model::{DisplaySystem, Model, Procedure};
use std::collections::BTreeMap;

pub(super) fn spectra(blocks: &[Block<'_>], units: &DisplaySystem) -> Result<BTreeMap<String, Vec<SpectrumPoint>>> {
    let mut out = BTreeMap::new();
    for b in blocks.iter().filter(|b| b.keyword == "PSD-DEFINITION") {
        allow(b, &["NAME", "TYPE"])?;
        if canonical(value(b, "TYPE", "FORCE")?) != "FORCE" {
            return Err(unsupported(b.line, "random import supports TYPE=FORCE only"));
        }
        let mut points = Vec::new();
        for row in &b.rows {
            if points.len() >= 10000 {
                return Err(error(row.line, "PSD point limit reached"));
            }
            let f = row_fields(row.text, row.line)?;
            let [real, imaginary, frequency] = f.as_slice() else {
                return Err(error(row.line, "PSD needs real, imaginary, frequency"));
            };
            let real = number(real, row.line)?;
            if real < 0.0 || number(imaginary, row.line)? != 0.0 {
                return Err(unsupported(row.line, "PSD must be real and nonnegative"));
            }
            points.push(SpectrumPoint {
                frequency: quantity(units, Quantity::Frequency, number(frequency, row.line)?, row.line)?,
                psd: quantity(units, Quantity::Time, real, row.line)?,
            });
        }
        if points.len() < 2 || points.windows(2).any(|w| w[0].frequency >= w[1].frequency) || points.iter().any(|p| p.frequency < 0.0) {
            return Err(error(b.line, "PSD requires at least two ascending nonnegative frequencies"));
        }
        if points.windows(2).any(|w| w[0].psd != w[1].psd) {
            return Err(unsupported(b.line, "nonconstant PSD requires explicit logarithmic interpolation mapping"));
        }
        if out.insert(name(b.required("NAME")?, b.line)?, points).is_some() {
            return Err(error(b.line, "duplicate PSD name"));
        }
    }
    Ok(out)
}

pub(super) fn procedure(
    body: &[Block<'_>],
    model: &Model,
    spectra: &BTreeMap<String, Vec<SpectrumPoint>>,
    units: &DisplaySystem,
) -> Result<Procedure> {
    let random = body.iter().find(|b| b.keyword == "RANDOMRESPONSE").ok_or_else(|| error(1, "missing random procedure"))?;
    allow(random, &[])?;
    if body.iter().filter(|b| b.keyword == "RANDOMRESPONSE").count() != 1 {
        return Err(error(random.line, "multiple random procedures in one step"));
    }
    let source = model.steps().last().ok_or_else(|| unsupported(random.line, "random response needs a preceding frequency step"))?;
    let Procedure::Frequency { modes, min_hz: None, max_hz: None } = source.procedure else {
        return Err(unsupported(random.line, "random response needs an immediately preceding unwindowed frequency extraction"));
    };
    if random.rows.len() != 1 {
        return Err(unsupported(random.line, "random import currently maps one frequency band"));
    }
    let row = random.rows.first().ok_or_else(|| error(random.line, "missing random frequency range"))?;
    let f = row_fields(row.text, row.line)?;
    if !(2..=5).contains(&f.len()) {
        return Err(error(row.line, "invalid random frequency range"));
    }
    let get = |i, default| f.get(i).filter(|v: &&String| !v.is_empty()).map_or(Ok(default), |s| number(s, row.line));
    let count = get(2, 20.0)?;
    let count = if count < 2.0 { 20.0 } else { count };
    if count.fract() != 0.0 || count > fieldcraft_model::harmonic::MAX_FREQUENCIES as f64 {
        return Err(error(row.line, "invalid random sample count"));
    }
    let scale = get(4, 0.0)?;
    if scale != 0.0 && scale != 1.0 {
        return Err(error(row.line, "frequency scale must be 0 or 1"));
    }
    let band = ModalBand {
        lower: quantity(units, Quantity::Frequency, get(0, 0.0)?, row.line)?,
        upper: quantity(units, Quantity::Frequency, get(1, 0.0)?, row.line)?,
        points: count as usize,
        bias: get(3, 3.0)?,
        logarithmic: scale == 0.0,
    };
    let correlations: Vec<_> = body.iter().filter(|b| b.keyword == "CORRELATION").collect();
    let [correlation] = correlations.as_slice() else {
        return Err(unsupported(random.line, "exactly one correlation definition is supported"));
    };
    allow(correlation, &["TYPE", "PSD", "COMPLEX"])?;
    if canonical(value(correlation, "TYPE", "CORRELATED")?) != "CORRELATED"
        || canonical(value(correlation, "COMPLEX", "NO")?) != "NO"
        || correlation.rows.len() != 1
    {
        return Err(unsupported(correlation.line, "only one real fully correlated load case is mapped"));
    }
    let row = correlation.rows.first().ok_or_else(|| error(correlation.line, "missing correlation row"))?;
    let f = row_fields(row.text, row.line)?;
    let [case, coefficient] = f.as_slice() else {
        return Err(error(row.line, "correlation needs load case, real scaling"));
    };
    let case = label(case, row.line)?;
    if number(coefficient, row.line)? != 1.0 {
        return Err(unsupported(row.line, "non-unit correlation scaling is not mapped yet"));
    }
    let spectrum =
        spectra.get(&name(correlation.required("PSD")?, correlation.line)?).ok_or_else(|| error(correlation.line, "missing PSD definition"))?.clone();
    if spectrum.first().is_none_or(|p| p.frequency > band.lower) || spectrum.last().is_none_or(|p| p.frequency < band.upper) {
        return Err(unsupported(random.line, "random band outside PSD definition"));
    }
    let mut damping = vec![0.0; modes as usize];
    let mut assigned = vec![false; modes as usize];
    for b in body {
        match b.keyword.as_str() {
            "MODALDAMPING" => {
                allow(b, &["DEFINITION"])?;
                if canonical(value(b, "DEFINITION", "MODENUMBERS")?) != "MODENUMBERS" {
                    return Err(unsupported(b.line, "modal damping by frequency is not mapped"));
                }
                for row in &b.rows {
                    let f = row_fields(row.text, row.line)?;
                    let [first, last, ratio] = f.as_slice() else {
                        return Err(error(row.line, "modal damping needs first,last,ratio"));
                    };
                    let (first, last) = (label(first, row.line)?, label(last, row.line)?);
                    let ratio = number(ratio, row.line)?;
                    if first > last || last > modes || ratio < 0.0 {
                        return Err(error(row.line, "invalid modal damping range"));
                    }
                    for mode in first..=last {
                        let i = (mode - 1) as usize;
                        let slot = assigned.get_mut(i).ok_or_else(|| error(row.line, "modal index out of bounds"))?;
                        if *slot {
                            return Err(error(row.line, "overlapping modal damping definitions"));
                        }
                        *slot = true;
                        if let Some(v) = damping.get_mut(i) {
                            *v = ratio;
                        }
                    }
                }
            }
            "CLOAD" => {
                if case == 1 || label(b.required("LOADCASE")?, b.line)? != case {
                    return Err(unsupported(b.line, "concentrated random loads require a single explicit case other than distributed-load case 1"));
                }
            }
            "DLOAD" | "DSLOAD" if case != 1 => return Err(unsupported(b.line, "distributed random loads belong to case 1")),
            _ => {}
        }
    }
    let r = RandomVibration {
        frequencies: Vec::new(),
        modal: Some(ModalRandom { source_step: source.name.clone(), modes, damping_ratios: damping, bands: vec![band] }),
        excitations: vec![Excitation::CorrelatedLoads { spectrum }],
        damping: Default::default(),
    };
    r.validate().map_err(|e| error(random.line, e.to_string()))?;
    Ok(Procedure::RandomVibration(r))
}
