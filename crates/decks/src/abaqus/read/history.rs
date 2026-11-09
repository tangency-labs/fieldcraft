//! Conservative analysis mapping from public Abaqus keyword semantics, not solver observations.
//! Sources: .../SIMACAEKEYRefMap/simakey-r-{step,static,frequency,buckle,boundary,
//! cload,dload,dsload,amplitude}.htm at https://abaqus.uclouvain.be/English/ .
//! Static amplitudes are evaluated at equilibrium end time; no transient integration is implied.
use super::{Block, Limits, Read, canonical, check_work, error, label, lex, number, row_fields};
use crate::{Error, Result};
use fieldcraft_model::units::{Preset, Quantity};
use fieldcraft_model::{BoundaryCondition, DisplaySystem, Distribution, ElementFace, Load, Model, Procedure, Set, Step, Target};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone)]
struct Curve {
    points: Vec<(f64, f64)>,
    smooth: bool,
    total: bool,
}
impl Curve {
    fn at(&self, x: f64) -> f64 {
        let Some(first) = self.points.first() else {
            return 0.0;
        };
        if x <= first.0 {
            return first.1;
        }
        for pair in self.points.windows(2) {
            let (a, b) = (pair[0], pair[1]);
            if x <= b.0 {
                let t = (x - a.0) / (b.0 - a.0);
                let w = if self.smooth { t * t * t * (10.0 + t * (-15.0 + 6.0 * t)) } else { t };
                return (1.0 - w) * a.1 + w * b.1;
            }
        }
        self.points.last().map_or(0.0, |p| p.1)
    }
}
#[derive(Clone)]
struct Driven {
    value: f64,
    amplitude: Option<String>,
    line: usize,
    multiple: bool,
}
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Family {
    Boundary,
    Force,
    PressureD,
    PressureS,
    Gravity,
}
type State = BTreeMap<(Family, String, u32, u8), Driven>;
type Touched = BTreeSet<(Family, String, u32, u8)>;

pub(super) fn unsupported(line: usize, message: impl Into<String>) -> Error {
    Error::Unsupported { format: "Abaqus analysis import", what: format!("line {line}: {}", message.into()) }
}
pub(super) fn allow(b: &Block<'_>, keys: &[&str]) -> Result<()> {
    if let Some(key) = b.params.keys().find(|k| !keys.contains(&k.as_str())) {
        return Err(unsupported(b.line, format!("*{} parameter {key}", b.keyword)));
    }
    Ok(())
}
pub(super) fn value<'a>(b: &'a Block<'_>, key: &str, default: &'a str) -> Result<&'a str> {
    match b.params.get(key) {
        None => Ok(default),
        Some(Some(v)) if !v.is_empty() => Ok(v),
        _ => Err(error(b.line, format!("{key} requires a value"))),
    }
}
pub(super) fn name(s: &str, line: usize) -> Result<String> {
    Ok(lex::unquote(s, line)?.to_ascii_uppercase())
}
fn curves(blocks: &[Block<'_>]) -> Result<BTreeMap<String, Curve>> {
    let mut out = BTreeMap::new();
    for b in blocks.iter().filter(|b| b.keyword == "AMPLITUDE") {
        allow(b, &["NAME", "DEFINITION", "TIME", "VALUE"])?;
        let definition = canonical(value(b, "DEFINITION", "TABULAR")?);
        if !matches!(definition.as_str(), "TABULAR" | "SMOOTHSTEP") || canonical(value(b, "VALUE", "RELATIVE")?) != "RELATIVE" {
            return Err(unsupported(b.line, "only relative TABULAR/SMOOTH STEP amplitudes are mapped"));
        }
        let time = canonical(value(b, "TIME", "STEPTIME")?);
        if !matches!(time.as_str(), "STEPTIME" | "TOTALTIME") {
            return Err(error(b.line, "invalid amplitude TIME"));
        }
        let mut points = Vec::new();
        for row in &b.rows {
            let f = row_fields(row.text, row.line)?;
            if f.len() % 2 != 0 {
                return Err(error(row.line, "amplitude requires time/value pairs"));
            }
            for pair in f.as_chunks::<2>().0 {
                if points.len() >= 100_000 {
                    return Err(error(row.line, "amplitude point limit reached"));
                }
                let x = number(&pair[0], row.line)?;
                let y = number(&pair[1], row.line)?;
                if x < 0.0 || points.last().is_some_and(|p: &(f64, f64)| p.0 >= x) {
                    return Err(error(row.line, "amplitude abscissae must be nonnegative and increasing"));
                }
                points.push((x, y));
            }
        }
        if points.is_empty() {
            return Err(error(b.line, "empty amplitude"));
        }
        let key = name(b.required("NAME")?, b.line)?;
        if out.insert(key, Curve { points, smooth: definition == "SMOOTHSTEP", total: time == "TOTALTIME" }).is_some() {
            return Err(error(b.line, "duplicate amplitude"));
        }
    }
    Ok(out)
}

/// The import is conservative: a retained unknown physical option cannot silently disappear
/// from a runnable step. Unsupported history stops this and subsequent step mappings.
pub(super) fn import(text: &str, read: &mut Read, preset: Preset, limits: Limits) -> Result<()> {
    let blocks = lex::blocks(text, super::MAX_RECORDS)?;
    let curves = match curves(&blocks) {
        Ok(v) => v,
        Err(e @ Error::Unsupported { .. }) => {
            read.warnings.push(e.to_string());
            return Ok(());
        }
        Err(e) => return Err(e),
    };
    let units = DisplaySystem::new(preset);
    let spectra = match super::random::spectra(&blocks, &units) {
        Ok(v) => v,
        Err(e @ Error::Unsupported { .. }) => {
            read.warnings.push(e.to_string());
            return Ok(());
        }
        Err(e) => return Err(e),
    };
    let mut in_step = false;
    let unknown: BTreeSet<usize> = read.unknown.iter().map(|u| u.line).collect();
    for b in &blocks {
        if b.keyword == "STEP" {
            in_step = true;
        }
        if !in_step && unknown.contains(&b.line) && !matches!(b.keyword.as_str(), "BOUNDARY" | "AMPLITUDE" | "PREPRINT" | "PSD-DEFINITION") {
            read.warnings.push(format!("analysis not imported: unhandled model option *{} at line {}; original retained", b.keyword, b.line));
            return Ok(());
        }
        if b.keyword == "ENDSTEP" {
            in_step = false;
        }
    }
    let mut state = State::new();
    let mut consumed = BTreeSet::new();
    let mut total = 0.0;
    let mut sequence = 0usize;
    let mut base_state: Option<String> = None;
    let mut index = 0;
    while index < blocks.len() {
        let Some(b) = blocks.get(index) else {
            break;
        };
        if matches!(b.keyword.as_str(), "AMPLITUDE" | "PSD-DEFINITION") {
            consumed.insert(b.line);
            index += 1;
            continue;
        }
        if b.keyword == "BOUNDARY" {
            apply(b, &mut state, &read.model, &curves, &mut BTreeMap::new(), &mut Touched::new(), limits)?;
            if state.values().any(|v| v.value != 0.0 || v.amplitude.is_some()) {
                return Err(unsupported(b.line, "initial boundaries must be zero and have no amplitude"));
            }
            consumed.insert(b.line);
            index += 1;
            continue;
        }
        if b.keyword != "STEP" {
            index += 1;
            continue;
        }
        let end = (index + 1..blocks.len())
            .find(|j| blocks.get(*j).is_some_and(|b| b.keyword == "ENDSTEP"))
            .ok_or_else(|| error(b.line, "STEP has no END STEP"))?;
        let body = blocks.get(index + 1..end).ok_or_else(|| error(b.line, "invalid step range"))?;
        let prepared = prepare(b, body, &state, &read.model, &HistoryData { curves: &curves, spectra: &spectra }, &units, limits);
        let (procedure, duration, mut next) = match prepared {
            Ok(v) => v,
            Err(e @ Error::Unsupported { .. }) => {
                read.warnings.push(format!("{e}; this and subsequent steps retained but not made runnable"));
                break;
            }
            Err(e) => return Err(e),
        };
        let perturbation = b.params.contains_key("PERTURBATION")
            || matches!(
                procedure,
                Procedure::Frequency { .. }
                    | Procedure::Buckle { .. }
                    | Procedure::Harmonic(_)
                    | Procedure::RandomVibration(_)
                    | Procedure::MatrixGeneration { .. }
            );
        if total > 0.0
            && (matches!(
                procedure,
                Procedure::Buckle { .. } | Procedure::Harmonic(_) | Procedure::RandomVibration(_) | Procedure::MatrixGeneration { .. }
            ) || (b.params.contains_key("PERTURBATION") && matches!(procedure, Procedure::Static {})))
        {
            read.warnings.push(format!(
                "line {}: base-state perturbation requires the prestress interface; this and subsequent steps are retained only",
                b.line
            ));
            break;
        }
        if matches!(procedure, Procedure::Frequency { .. }) {
            next.retain(|(family, _, _, _), _| *family == Family::Boundary);
            for driven in next.values_mut() {
                driven.value = 0.0;
                driven.amplitude = None;
            }
        }
        let step_name = lex::unquote(value(b, "NAME", "")?, b.line)?;
        let step_name = if step_name.is_empty() { format!("ImportedStep{}", read.model.steps().len() + 1) } else { step_name };
        let mut step = Step::new(step_name, procedure);
        if matches!(step.procedure, Procedure::Frequency { .. }) {
            step.base_state = base_state.clone();
        }
        if matches!(step.procedure, Procedure::Static {} | Procedure::NonlinearStatic(_)) {
            base_state = Some(step.name.clone());
        }
        materialize(
            &mut read.model,
            &mut step,
            &mut next,
            &Evaluation { curves: &curves, time: duration, total: total + duration, units: &units, limits },
            &mut sequence,
        )?;
        read.model.add_step(step, None).map_err(|e| error(b.line, e.to_string()))?;
        for part in blocks.get(index..=end).unwrap_or_default() {
            if !is_output(&part.keyword) {
                consumed.insert(part.line);
            }
        }
        if !perturbation {
            state = next;
            total += duration;
            if !total.is_finite() {
                return Err(error(b.line, "total analysis time overflow"));
            }
        }
        index = end + 1;
    }
    read.unknown.retain(|u| !consumed.contains(&u.line));
    read.unknown_count = read.unknown.len();
    read.analysis_complete = index == blocks.len();
    if !curves.is_empty() {
        read.warnings.push("static amplitude curves are evaluated at the equilibrium step end; the original curves remain in the saved source, not a transient solution".into());
    }
    if read.unknown.iter().any(|u| is_output(&u.keyword)) {
        read.warnings.push("output requests retained; native procedure default outputs are used".into());
    }
    Ok(())
}
fn is_output(k: &str) -> bool {
    matches!(k, "OUTPUT" | "NODEOUTPUT" | "ELEMENTOUTPUT" | "NODEPRINT" | "ELPRINT" | "ENERGYOUTPUT" | "PREPRINT")
}

struct HistoryData<'a> {
    curves: &'a BTreeMap<String, Curve>,
    spectra: &'a BTreeMap<String, Vec<fieldcraft_model::random::SpectrumPoint>>,
}

fn prepare(
    b: &Block<'_>,
    body: &[Block<'_>],
    state: &State,
    model: &Model,
    data: &HistoryData<'_>,
    units: &DisplaySystem,
    limits: Limits,
) -> Result<(Procedure, f64, State)> {
    let curves = data.curves;
    let spectra = data.spectra;
    allow(b, &["NAME", "NLGEOM", "PERTURBATION", "AMPLITUDE", "INC"])?;
    let nonlinear = match b.params.get("NLGEOM") {
        Some(None) => true,
        _ => match canonical(value(b, "NLGEOM", "NO")?).as_str() {
            "YES" => true,
            "NO" => false,
            _ => return Err(error(b.line, "NLGEOM expects YES or NO")),
        },
    };
    if nonlinear && (!model.steps().is_empty() || b.params.contains_key("PERTURBATION")) {
        return Err(unsupported(b.line, "NLGEOM: sequential equilibrium or PERTURBATION is not mapped"));
    }
    if nonlinear && canonical(value(b, "AMPLITUDE", "RAMP")?) != "RAMP" {
        return Err(unsupported(b.line, "NLGEOM with STEP AMPLITUDE other than RAMP"));
    }
    if b.params.get("PERTURBATION").is_some_and(Option::is_some) {
        return Err(error(b.line, "PERTURBATION is a flag"));
    }
    if !matches!(canonical(value(b, "AMPLITUDE", "RAMP")?).as_str(), "RAMP" | "STEP") {
        return Err(unsupported(b.line, "STEP amplitude setting"));
    }
    if !b.rows.is_empty() {
        return Err(unsupported(b.line, "step description rows are not mapped"));
    }
    let random = body.iter().any(|b| b.keyword == "RANDOMRESPONSE");
    let mut procedure = if random { Some(super::random::procedure(body, model, spectra, units)?) } else { None };
    let mut duration = 1.0;
    let mut global_damping = None;
    let mut next = state.clone();
    let mut ops = BTreeMap::new();
    let mut touched = Touched::new();
    for part in body {
        match part.keyword.as_str() {
            "MATRIXGENERATE" => {
                if procedure.is_some() {
                    return Err(error(part.line, "multiple analysis procedures"));
                }
                procedure = Some(super::matrices::procedure(part, body)?);
            }
            "MATRIXOUTPUT" if body.iter().any(|b| b.keyword == "MATRIXGENERATE") => {}
            "STATIC" => {
                allow(part, &[])?;
                if procedure.is_some() {
                    return Err(error(part.line, "multiple analysis procedures in one STEP"));
                }
                if part.rows.len() > 1 {
                    return Err(error(part.line, "STATIC has at most one data row"));
                }
                if let Some(row) = part.rows.first() {
                    let f = row_fields(row.text, row.line)?;
                    if f.len() > 4 {
                        return Err(error(row.line, "STATIC expects at most four increment/time entries"));
                    }
                    for token in f.iter().filter(|s| !s.is_empty()) {
                        if number(token, row.line)? < 0.0 {
                            return Err(error(row.line, "STATIC increments/time must be nonnegative"));
                        }
                    }
                    duration = f.get(1).filter(|s| !s.is_empty()).map(|s| number(s, row.line)).transpose()?.filter(|v| *v > 0.0).unwrap_or(1.0);
                }
                procedure = Some(if nonlinear { super::nonlinear::procedure(b, part, model)? } else { Procedure::Static {} });
            }
            "STEADYSTATEDYNAMICS" => {
                if procedure.is_some() {
                    return Err(error(part.line, "multiple analysis procedures in one STEP"));
                }
                procedure = Some(Procedure::Harmonic(super::dynamic::harmonic(part, units)?));
            }
            "GLOBALDAMPING" => {
                if global_damping.replace(super::dynamic::damping(part, units)?).is_some() {
                    return Err(error(part.line, "duplicate GLOBAL DAMPING"));
                }
            }
            "FREQUENCY" | "BUCKLE" => {
                allow(part, &["EIGENSOLVER", "NORMALIZATION"])?;
                if canonical(value(part, "EIGENSOLVER", "LANCZOS")?) != "LANCZOS" {
                    return Err(unsupported(part.line, "non-Lanczos eigen extraction"));
                }
                if !matches!(canonical(value(part, "NORMALIZATION", "MASS")?).as_str(), "MASS" | "DISPLACEMENT") {
                    return Err(error(part.line, "invalid modal normalization"));
                }
                if procedure.is_some() {
                    return Err(error(part.line, "multiple procedures"));
                }
                let row = part.rows.first().ok_or_else(|| error(part.line, "eigen extraction needs a mode count"))?;
                if part.rows.len() != 1 {
                    return Err(error(part.line, "eigen extraction needs one data row"));
                }
                let f = row_fields(row.text, row.line)?;
                let modes = label(f.first().map(String::as_str).unwrap_or_default(), row.line)?;
                if modes > fieldcraft_model::steps::MAX_MODES {
                    return Err(error(row.line, "mode count exceeds limit"));
                }
                let bound = |i: usize| f.get(i).filter(|s| !s.is_empty()).map(|s| number(s, row.line)).transpose();
                procedure = Some(if part.keyword == "FREQUENCY" {
                    if f.len() > 3 {
                        return Err(unsupported(row.line, "additional eigen extraction controls"));
                    }
                    Procedure::Frequency { modes, min_hz: bound(1)?, max_hz: bound(2)? }
                } else {
                    if f.len() > 1 {
                        return Err(unsupported(row.line, "buckling shift/extraction controls require explicit mapping"));
                    }
                    Procedure::Buckle { modes, shift: None }
                });
            }
            "RANDOMRESPONSE" | "MODALDAMPING" | "CORRELATION" if random => {}
            "CLOAD" if random && part.params.contains_key("LOADCASE") => {
                let mut params = part.params.clone();
                params.remove("LOADCASE");
                let mapped = Block {
                    line: part.line,
                    keyword: part.keyword.clone(),
                    params,
                    rows: part.rows.iter().map(|r| lex::Row { line: r.line, text: r.text }).collect(),
                };
                apply(&mapped, &mut next, model, curves, &mut ops, &mut touched, limits)?;
            }
            "BOUNDARY" | "CLOAD" | "DLOAD" | "DSLOAD" => apply(part, &mut next, model, curves, &mut ops, &mut touched, limits)?,
            k if is_output(k) => {}
            _ => return Err(unsupported(part.line, format!("*{} in STEP", part.keyword))),
        }
    }
    let mut procedure = procedure.ok_or_else(|| unsupported(b.line, "STEP has no supported procedure"))?;
    if let Some(damping) = global_damping {
        let Procedure::Harmonic(h) = &mut procedure else {
            return Err(unsupported(b.line, "GLOBAL DAMPING outside harmonic analysis"));
        };
        h.damping = damping;
    }
    if matches!(procedure, Procedure::Harmonic(_)) && next.iter().any(|((family, _, _, _), v)| *family == Family::Boundary && v.value != 0.0) {
        return Err(unsupported(b.line, "nonzero harmonic prescribed motion"));
    }
    if !matches!(procedure, Procedure::Static {}) && next.values().any(|v| v.amplitude.is_some()) {
        return Err(unsupported(b.line, "amplitude on non-static step"));
    }
    if matches!(procedure, Procedure::MatrixGeneration { .. })
        && next.iter().any(|((family, _, _, _), v)| *family != Family::Boundary || v.value != 0.0)
    {
        return Err(unsupported(b.line, "matrix export requires homogeneous boundaries and no active loads"));
    }
    if nonlinear && !matches!(procedure, Procedure::NonlinearStatic(_)) {
        return Err(unsupported(b.line, "NLGEOM=YES requires a supported *STATIC step"));
    }
    if nonlinear && next.values().any(|v| v.amplitude.is_some()) {
        return Err(unsupported(b.line, "NLGEOM with named load/boundary AMPLITUDE"));
    }
    if !nonlinear
        && model.steps().iter().any(|s| matches!(s.procedure, Procedure::NonlinearStatic(_)))
        && !matches!(procedure, Procedure::Frequency { .. })
    {
        return Err(unsupported(b.line, "NLGEOM continuation: subsequent general steps are not mapped"));
    }
    super::super::limited::enhanced(model, Some(&procedure))?;
    Ok((procedure, duration, next))
}

fn members(model: &Model, target: &str, nodes: bool, line: usize) -> Result<Vec<u32>> {
    let target = lex::unquote(target, line)?;
    if let Ok(id) = target.parse::<u32>() {
        if if nodes { model.node(id).is_some() } else { model.elements().contains_key(&id) } {
            return Ok(vec![id]);
        }
        return Err(error(line, format!("missing target label {id}")));
    }
    let set = model
        .sets()
        .iter()
        .find(|(key, _)| key.eq_ignore_ascii_case(&target))
        .map(|(_, set)| set)
        .ok_or_else(|| error(line, format!("missing target set {target}")))?;
    let labels = if nodes { &set.nodes } else { &set.elements };
    if labels.is_empty() {
        return Err(error(line, format!("target {target} has no {}", if nodes { "nodes" } else { "elements" })));
    }
    Ok(labels.clone())
}
fn apply(
    b: &Block<'_>,
    state: &mut State,
    model: &Model,
    curves: &BTreeMap<String, Curve>,
    ops: &mut BTreeMap<String, bool>,
    touched: &mut Touched,
    limits: Limits,
) -> Result<()> {
    allow(b, &["OP", "AMPLITUDE", "TYPE", "REAL"])?;
    if b.params.contains_key("REAL") && (b.keyword == "BOUNDARY" || b.params.get("REAL").is_some_and(Option::is_some)) {
        return Err(unsupported(b.line, "REAL must be a load flag"));
    }
    if b.keyword == "BOUNDARY" && canonical(value(b, "TYPE", "DISPLACEMENT")?) != "DISPLACEMENT" {
        return Err(unsupported(b.line, "boundary velocity/acceleration"));
    }
    if b.keyword != "BOUNDARY" && b.params.contains_key("TYPE") {
        return Err(unsupported(b.line, "load TYPE parameter"));
    }
    let op = canonical(value(b, "OP", "MOD")?);
    if !matches!(op.as_str(), "MOD" | "NEW") {
        return Err(error(b.line, "OP must be NEW or MOD"));
    }
    let new = op == "NEW";
    if let Some(old) = ops.get(&b.keyword) {
        if *old != new {
            return Err(unsupported(b.line, "mixed OP settings in one load family"));
        }
    } else {
        ops.insert(b.keyword.clone(), new);
        if new {
            state.retain(|(family, _, _, _), _| match b.keyword.as_str() {
                "BOUNDARY" => *family != Family::Boundary,
                "CLOAD" => *family != Family::Force,
                "DSLOAD" => *family != Family::PressureS,
                _ => !matches!(family, Family::PressureD | Family::Gravity),
            });
        }
    }
    let amplitude = b.value("AMPLITUDE").map(|s| name(s, b.line)).transpose()?;
    if b.params.contains_key("AMPLITUDE") && amplitude.is_none() {
        return Err(error(b.line, "AMPLITUDE requires a name"));
    }
    if amplitude.as_ref().is_some_and(|s| !curves.contains_key(s)) {
        return Err(error(b.line, "missing amplitude"));
    }
    for row in &b.rows {
        let f = row_fields(row.text, row.line)?;
        let target = f.first().map(String::as_str).unwrap_or_default();
        let mut insert = |family, node, dof, mut value| -> Result<()> {
            let identity = if family == Family::Boundary { String::new() } else { name(target, row.line)? };
            let key = (family, identity, node, dof);
            if state.len() >= limits.members && !state.contains_key(&key) {
                return Err(error(row.line, "analysis membership limit reached"));
            }
            let mut multiple = false;
            if family != Family::Boundary
                && let Some(old) = state.get(&key)
            {
                if touched.contains(&key) {
                    if old.amplitude != amplitude {
                        return Err(unsupported(row.line, "duplicate load rows with different amplitudes"));
                    }
                    value += old.value;
                    multiple = true;
                } else if old.multiple {
                    return Err(unsupported(row.line, "multiple identical load definitions must be removed with OP=NEW before modification"));
                }
            }
            touched.insert(key.clone());
            state.insert(key, Driven { value, amplitude: amplitude.clone(), line: row.line, multiple });
            Ok(())
        };
        match b.keyword.as_str() {
            "BOUNDARY" => {
                if !(2..=4).contains(&f.len()) {
                    return Err(error(row.line, "BOUNDARY expects target, dofs/type and optional magnitude"));
                }
                let kind = canonical(f.get(1).map(String::as_str).unwrap_or_default());
                let (dofs, magnitude) = if let Ok(first) = kind.parse::<u8>() {
                    let last = f
                        .get(2)
                        .filter(|s| !s.is_empty())
                        .map(|s| s.parse::<u8>().map_err(|_| error(row.line, "invalid boundary last DOF")))
                        .transpose()?
                        .unwrap_or(first);
                    if !(1..=6).contains(&first) || last < first || last > 6 {
                        return Err(unsupported(row.line, "only structural DOFs 1..6 are mapped"));
                    }
                    let magnitude = f.get(3).filter(|s| !s.is_empty()).map(|s| number(s, row.line)).transpose()?.unwrap_or(0.0);
                    ((first..=last).collect::<Vec<_>>(), magnitude)
                } else {
                    if f.len() != 2 {
                        return Err(error(row.line, "named boundary has no magnitude"));
                    }
                    let dofs = match kind.as_str() {
                        "ENCASTRE" => vec![1, 2, 3, 4, 5, 6],
                        "PINNED" => vec![1, 2, 3],
                        "XSYMM" => vec![1, 5, 6],
                        "YSYMM" => vec![2, 4, 6],
                        "ZSYMM" => vec![3, 4, 5],
                        "XASYMM" => vec![2, 3, 4],
                        "YASYMM" => vec![1, 3, 5],
                        "ZASYMM" => vec![1, 2, 6],
                        _ => return Err(unsupported(row.line, format!("boundary type {kind}"))),
                    };
                    (dofs, 0.0)
                };
                for node in members(model, target, true, row.line)? {
                    for dof in &dofs {
                        insert(Family::Boundary, node, *dof, magnitude)?;
                    }
                }
            }
            "CLOAD" => {
                if f.len() != 3 {
                    return Err(error(row.line, "CLOAD expects node/set, DOF, magnitude"));
                }
                let dof = f
                    .get(1)
                    .and_then(|s| s.parse::<u8>().ok())
                    .filter(|d| (1..=6).contains(d))
                    .ok_or_else(|| error(row.line, "CLOAD DOF must be 1..6"))?;
                let magnitude = number(f.get(2).map(String::as_str).unwrap_or_default(), row.line)?;
                for node in members(model, target, true, row.line)? {
                    insert(Family::Force, node, dof, magnitude)?;
                }
            }
            "DSLOAD" => {
                if f.len() != 3 || f.get(1).is_none_or(|v| canonical(v) != "P") {
                    return Err(unsupported(row.line, "only DSLOAD pressure P is mapped"));
                }
                let key = name(target, row.line)?;
                let set = model
                    .sets()
                    .iter()
                    .find(|(n, _)| n.eq_ignore_ascii_case(&key))
                    .map(|(_, s)| s)
                    .ok_or_else(|| error(row.line, "missing pressure surface"))?;
                if set.faces.is_empty() || !set.edges.is_empty() {
                    return Err(unsupported(row.line, "pressure requires element faces"));
                }
                let magnitude = number(f.get(2).map(String::as_str).unwrap_or_default(), row.line)?;
                for face in &set.faces {
                    insert(Family::PressureS, face.element, face.face, magnitude)?;
                }
            }
            _ => {
                let kind = canonical(f.get(1).map(String::as_str).unwrap_or_default());
                let elements = if target.is_empty() { model.elements().keys().copied().collect() } else { members(model, target, false, row.line)? };
                if kind == "GRAV" {
                    if !model.masses().is_empty() {
                        return Err(unsupported(row.line, "GRAV with imported point masses needs discrete load membership mapping"));
                    }
                    if f.len() != 6 {
                        return Err(error(row.line, "GRAV needs magnitude and three direction components"));
                    }
                    let magnitude = number(f.get(2).map(String::as_str).unwrap_or_default(), row.line)?;
                    let mut direction = [0.0; 3];
                    for (i, v) in direction.iter_mut().enumerate() {
                        *v = number(f.get(i + 3).map(String::as_str).unwrap_or_default(), row.line)?;
                    }
                    let norm = direction.iter().fold(0.0_f64, |v, x| v.hypot(*x));
                    if norm == 0.0 || !norm.is_finite() {
                        return Err(error(row.line, "GRAV direction must be finite and nonzero"));
                    }
                    for element in elements {
                        for (i, v) in direction.iter().enumerate() {
                            insert(Family::Gravity, element, (i + 1) as u8, magnitude * v / norm)?;
                        }
                    }
                } else if let Some(face) = kind.strip_prefix('P').and_then(|s| s.parse::<u8>().ok()) {
                    if f.len() != 3 {
                        return Err(error(row.line, "element pressure needs three entries"));
                    }
                    let magnitude = number(f.get(2).map(String::as_str).unwrap_or_default(), row.line)?;
                    for element in elements {
                        let e = model.elements().get(&element).ok_or_else(|| error(row.line, "missing pressure element"))?;
                        if super::super::face_id(e.kind, face).is_none() {
                            return Err(error(row.line, "invalid element pressure face"));
                        }
                        insert(Family::PressureD, element, face, magnitude)?;
                    }
                } else {
                    return Err(unsupported(row.line, format!("DLOAD {kind}")));
                }
            }
        }
    }
    Ok(())
}

struct Evaluation<'a> {
    curves: &'a BTreeMap<String, Curve>,
    time: f64,
    total: f64,
    units: &'a DisplaySystem,
    limits: Limits,
}
fn materialize(model: &mut Model, step: &mut Step, state: &mut State, evaluation: &Evaluation<'_>, sequence: &mut usize) -> Result<()> {
    let curves = evaluation.curves;
    let units = evaluation.units;
    let limits = evaluation.limits;
    let time = evaluation.time;
    let total = evaluation.total;
    let mut groups: BTreeMap<(Family, u8, u64, Option<String>), Vec<u32>> = BTreeMap::new();
    let mut accumulated: BTreeMap<(Family, u32, u8), f64> = BTreeMap::new();
    for ((family, _, id, dof), driven) in state.iter_mut() {
        check_work(groups.len())?;
        let mut v = driven.value;
        if let Some(name) = &driven.amplitude {
            let curve = curves.get(name).ok_or_else(|| error(driven.line, "missing amplitude"))?;
            v *= curve.at(if curve.total { total } else { time });
            if !curve.total {
                driven.value = v;
                driven.amplitude = None;
            }
        }
        let q = match family {
            Family::Boundary => {
                if *dof <= 3 {
                    Quantity::Length
                } else {
                    Quantity::Angle
                }
            }
            Family::Force => {
                if *dof <= 3 {
                    Quantity::Force
                } else {
                    Quantity::Moment
                }
            }
            Family::PressureD | Family::PressureS => Quantity::Pressure,
            _ => Quantity::Acceleration,
        };
        v = super::quantity(units, q, v, driven.line)?;
        *accumulated.entry((*family, *id, *dof)).or_default() += v;
    }
    for ((family, id, dof), v) in accumulated {
        if !v.is_finite() {
            return Err(error(1, "accumulated load is non-finite"));
        }
        let point = if matches!(family, Family::Boundary | Family::Force) {
            model.reference_points().iter().find(|(_, p)| p.label == id).map(|(n, _)| n.clone())
        } else {
            None
        };
        groups.entry((family, dof, v.to_bits(), point)).or_default().push(id);
    }
    for ((family, dof, bits, point), ids) in groups {
        if ids.len() > limits.members {
            return Err(error(1, "analysis member limit reached"));
        }
        let mut available = None;
        for _ in 0..=fieldcraft_model::MAX_DEFINITIONS {
            *sequence = sequence.saturating_add(1);
            let candidate = format!("InpAnalysis{sequence}");
            if model.set(&candidate).is_none() && !model.loads().contains_key(&candidate) && !model.boundary_conditions().contains_key(&candidate) {
                available = Some(candidate);
                break;
            }
        }
        let name = available.ok_or_else(|| error(1, "analysis name allocation limit reached"))?;
        let value = f64::from_bits(bits);
        let mut set = Set::default();
        match family {
            Family::Boundary | Family::Force => set.nodes = ids,
            Family::PressureD | Family::PressureS => set.faces = ids.into_iter().map(|element| ElementFace { element, face: dof }).collect(),
            _ => set.elements = ids,
        }
        let target = if let Some(point) = point {
            Target::ReferencePoint(point)
        } else {
            model.add_set(&name, set).map_err(|e| error(1, e.to_string()))?;
            Target::Set(name.clone())
        };
        if family == Family::Boundary {
            let mut values = [None; 6];
            if let Some(slot) = values.get_mut(usize::from(dof) - 1) {
                *slot = Some(value);
            }
            model.add_boundary_condition(&name, BoundaryCondition::Displacement { target, values }).map_err(|e| error(1, e.to_string()))?;
            step.boundary_conditions.insert(name);
        } else {
            let mut vector = [0.0; 3];
            if let Some(slot) = vector.get_mut((usize::from(dof) - 1) % 3) {
                *slot = value;
            }
            let load = match family {
                Family::Force if dof <= 3 => Load::Force { target, vector, distribution: Distribution::PerNode },
                Family::Force => Load::Moment { target, vector, distribution: Distribution::PerNode },
                Family::PressureD | Family::PressureS => Load::Pressure { set: name.clone(), magnitude: value },
                _ => Load::Gravity { acceleration: vector, set: Some(name.clone()) },
            };
            model.add_load(&name, load).map_err(|e| error(1, e.to_string()))?;
            step.loads.insert(name);
        }
    }
    // Eigen frequencies are converted only at the model boundary, like all other quantities.
    if let Procedure::Frequency { min_hz, max_hz, .. } = &mut step.procedure {
        for value in [min_hz, max_hz].into_iter().flatten() {
            *value = super::quantity(units, Quantity::Frequency, *value, 1)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    const MODEL: &str = "*NODE,NSET=Tip\n1,0,0,0\n2,1,0,0\n*BOUNDARY\n1,ENCASTRE\n";
    fn read(text: &str) -> Read {
        super::super::read_analysis_limited(text, Preset::MKgS, Limits::default()).unwrap()
    }
    fn force(read: &Read, step: &str) -> f64 {
        read.model
            .step(step)
            .unwrap()
            .loads
            .iter()
            .map(|n| match read.model.loads().get(n).unwrap() {
                Load::Force { vector, target: Target::Set(set), .. } => vector[0] * read.model.set(set).unwrap().nodes.len() as f64,
                _ => 0.0,
            })
            .sum()
    }
    #[test]
    fn amplitudes_step_and_total_time_propagate_differently() {
        for (time, expected) in [("STEP TIME", 500.0), ("TOTAL TIME", 1000.0)] {
            let input = format!(
                "{MODEL}*AMPLITUDE,NAME=A,TIME={time}\n0,0,2,2\n*STEP,NAME=One\n*STATIC\n1,1\n*CLOAD,AMPLITUDE=A\n2,1,500\n*END STEP\n*STEP,NAME=Two\n*STATIC\n1,1\n*END STEP\n"
            );
            let got = read(&input);
            assert!(got.analysis_complete);
            assert_eq!(force(&got, "One"), 500.0);
            assert_eq!(force(&got, "Two"), expected);
        }
    }
    #[test]
    fn load_identity_new_and_mod_follow_original_target_not_expanded_nodes() {
        let input = format!(
            "{MODEL}*STEP,NAME=One\n*STATIC\n*CLOAD\nTip,1,10\n*END STEP\n*STEP,NAME=Two\n*STATIC\n*CLOAD\n1,1,5\n*END STEP\n*STEP,NAME=Three\n*STATIC\n*CLOAD,OP=NEW\n2,1,7\n*END STEP\n"
        );
        let got = read(&input);
        assert_eq!(force(&got, "One"), 20.0);
        assert_eq!(force(&got, "Two"), 25.0);
        assert_eq!(force(&got, "Three"), 7.0);
    }
    #[test]
    fn unknown_physics_does_not_create_a_silently_incomplete_runnable_step() {
        let got = read(&format!("{MODEL}*STEP,NAME=Unsafe\n*STATIC\n*FUTURE PHYSICS\n10\n*END STEP\n"));
        assert!(!got.analysis_complete);
        assert!(got.model.steps().is_empty());
        assert!(got.unknown.iter().any(|u| u.keyword == "FUTUREPHYSICS"));
    }
    #[test]
    fn malformed_amplitude_and_missing_boundary_target_are_errors() {
        for text in [format!("{MODEL}*AMPLITUDE,NAME=A\n1,1,0,2\n"), format!("{MODEL}*STEP\n*STATIC\n*BOUNDARY\nMissing,1,3\n*END STEP\n")] {
            assert!(super::super::read_analysis_limited(&text, Preset::MKgS, Limits::default()).is_err());
        }
    }
    #[test]
    fn smooth_step_endpoint_and_frequency_become_native_definitions() {
        let got = read(&format!(
            "{MODEL}*AMPLITUDE,NAME=A,DEFINITION=SMOOTH STEP\n0,0,2,2\n*STEP,NAME=Pull\n*STATIC\n1,1\n*CLOAD,AMPLITUDE=A\n2,1,100\n*END STEP\n*STEP,NAME=Modes,PERTURBATION\n*FREQUENCY,EIGENSOLVER=LANCZOS\n3\n*END STEP\n"
        ));
        assert_eq!(force(&got, "Pull"), 100.0);
        assert!(matches!(got.model.step("Modes").unwrap().procedure, Procedure::Frequency { modes: 3, .. }));
        assert!(got.model.step("Modes").unwrap().loads.is_empty());
    }
}
