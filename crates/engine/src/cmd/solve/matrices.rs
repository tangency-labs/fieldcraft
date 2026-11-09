//! Native operator artifacts, independently of response results.
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};
use fieldcraft_io::roots::Access;
use fieldcraft_model::steps::MatrixKind;
use serde_json::{Value, json};
use std::fmt::Write;

pub const EXPORT: CommandSpec = CommandSpec {
    id: "solve.matrices",
    label: "Export Structural Matrices",
    params: r#"{"step":"Matrices", "path":"operators.json"} — assembles the requested unloaded linear structural operators with homogeneous constraints eliminated. Optional mass: consistent (default) or lumped selects inertia and material damping. Step Rayleigh damping is excluded. Optional path writes an SI manifest and sibling MATRIX INPUT five-column .mtx files (node, DOF, node, DOF, value; lower triangle by node/DOF labels). Rotations use radians. Without path returns metadata only. No factorization, response solve or result frames; each output file is atomic, the collection is not a transaction. Limits/threads as solve.static."#,
    form: Some(super::super::forms::SOLVE_MATRICES),
    kind: CommandKind::Model,
    requires: &[Requirement::Analysable, Requirement::Step],
    run,
};
fn key(k: MatrixKind) -> &'static str {
    match k {
        MatrixKind::Stiffness => "stiffness",
        MatrixKind::Mass => "mass",
        MatrixKind::ViscousDamping => "viscous_damping",
        MatrixKind::StructuralDamping => "structural_damping",
    }
}
pub(super) fn run(s: &mut Session, p: &Value) -> Result<Value> {
    let path = if p.get("path").is_some() { Some(super::super::file_param(s, p, "path", Access::Write)?) } else { None };
    let limits = super::limits(p)?;
    let threads = super::threads(s, p)?;
    let mass = super::opt::<fieldcraft_solver::MassKind>(p, "mass")?.unwrap_or_default();
    let r = super::with_threads(threads, || {
        fieldcraft_solver::export_matrices_with_mass(&s.model, p.get("step").and_then(Value::as_str), mass, &limits)
    })?;
    let mut out = json!({"step":r.step,"mass":mass,"equations":r.degrees.len(),"degrees":r.degrees,"units":"SI (m, kg, s, N, rad)","triangle":"lower by native node label and DOF", "constraint_reduction":"T^T A T; homogeneous prescribed and dependent coordinates eliminated", "convention":"K - omega^2 M + i(omega C + G)","warnings":r.warnings,"exported":path.is_some(),"matrices":{}});
    let mut destinations = Vec::new();
    if let Some(path) = &path {
        let manifest = std::path::PathBuf::from(super::super::file::absolute(path));
        let mut unique = std::collections::BTreeSet::from([s.file_access.check(&manifest, Access::Write)?]);
        let stem = manifest.file_stem().and_then(|s| s.to_str()).ok_or_else(|| EngineError::invalid("path", "manifest needs UTF-8 filename"))?;
        for &kind in r.matrices.keys() {
            let target = manifest.with_file_name(format!("{stem}-{}.mtx", key(kind)));
            let target = s.file_access.check(&target, Access::Write)?;
            if !unique.insert(target.clone()) {
                return Err(EngineError::invalid("path", "matrix output destinations alias"));
            }
            if target.exists() && !target.is_file() {
                return Err(EngineError::invalid("path", "matrix target is not a regular file"));
            }
            destinations.push((kind, target));
        }
    }
    for (&kind, a) in &r.matrices {
        let count = (0..a.n)
            .map(|i| {
                let (cols, vals) = a.row(i);
                cols.iter().zip(vals).filter(|(j, v)| r.degrees.get(**j) <= r.degrees.get(i) && (**v != 0.0 || **j == i)).count()
            })
            .sum::<usize>();
        out["matrices"][key(kind)] = json!({"entries":count});
    }
    for (kind, path) in destinations {
        let a = r.matrices.get(&kind).ok_or_else(|| EngineError::invalid("step", "missing operator"))?;
        let mut text = format!("** FieldCraft {} ; SI ; symmetric lower triangle by native node label and DOF\n", key(kind));
        for i in 0..a.n {
            let (node, dof) = r.degrees.get(i).ok_or_else(|| EngineError::invalid("step", "missing matrix row"))?;
            let (cols, vals) = a.row(i);
            for (&j, &v) in cols.iter().zip(vals) {
                let (other, axis) = r.degrees.get(j).ok_or_else(|| EngineError::invalid("step", "missing matrix column"))?;
                if (other, axis) > (node, dof) || (v == 0.0 && j != i) {
                    continue;
                }
                writeln!(text, "{node}, {dof}, {other}, {axis}, {v:.17e}").map_err(|e| EngineError::invalid("path", e.to_string()))?;
            }
            fieldcraft_model::progress::check()?;
        }
        fieldcraft_io::fs::write_output(&path, text.as_bytes())?;
        out["matrices"][key(kind)]["path"] = json!(path);
    }
    if let Some(path) = path {
        fieldcraft_io::fs::write_output(
            std::path::Path::new(&path),
            serde_json::to_string_pretty(&out).map_err(|e| EngineError::invalid("path", e.to_string()))?.as_bytes(),
        )?;
    }
    Ok(out)
}
