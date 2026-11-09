//! Inspect and rewrite a native result archive without solving or guessing its producer.
//! Usage: results_roundtrip input.fieldcraft output-directory
use fieldcraft_io::{
    Extras,
    results::{Solved, StoredRef},
};
use fieldcraft_solver::StepResult;
use serde_json::{Value, json};
use std::{error::Error, fs, path::PathBuf};

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args_os().skip(1);
    let input = PathBuf::from(args.next().ok_or("input archive required")?);
    let out = PathBuf::from(args.next().ok_or("new output directory required")?);
    if args.next().is_some() {
        return Err("unexpected arguments".into());
    }
    fs::create_dir(&out)?;
    let bytes = fs::read(&input)?;
    let f = fieldcraft_io::read(&bytes)?;
    if f.results_note.is_some() || f.results.is_empty() {
        return Err("no fully decoded results".into());
    }
    let mut fields = Vec::new();
    for stored in &f.results {
        let (kind, payload, fatigue): (&str, Value, Value) = match &stored.result {
            Solved::Implicit(StepResult::Static(r)) => ("static", serde_json::to_value(r)?, serde_json::to_value(&r.fatigue)?),
            Solved::Implicit(StepResult::Harmonic(r)) => ("harmonic", serde_json::to_value(r)?, serde_json::to_value(&r.fatigue)?),
            Solved::Implicit(StepResult::Dynamics(r)) => ("dynamics", serde_json::to_value(r)?, Value::Null),
            Solved::Implicit(StepResult::Frequency(r)) => ("frequency", serde_json::to_value(r)?, Value::Null),
            Solved::Implicit(StepResult::Heat(r)) => ("heat", serde_json::to_value(r)?, Value::Null),
            Solved::Explicit(r) => ("explicit", serde_json::to_value(r)?, Value::Null),
            Solved::Nonlinear(r) => ("nonlinear_rod", serde_json::to_value(r)?, Value::Null),
            Solved::Fluid(r) => ("fluid", serde_json::to_value(r)?, Value::Null),
        };
        fields.push(json!({"kind":kind,"definition":stored.definition,"payload":payload,"fatigue":fatigue}));
    }
    let extras = Extras {
        deck_source: f.deck_source.as_ref(),
        results: f.results.iter().map(|r| (r.result.step(), StoredRef { definition: r.definition, result: r.result.as_ref() })).collect(),
        embed: f.embedded.iter().map(|(k, v)| (k.clone(), v.as_slice())).collect(),
        run: f.run.clone(),
    };
    let rewritten = fieldcraft_io::write_with(&f.model, &f.geometry, &f.journal, &f.journal_header, &f.view, &extras)?;
    let reopened = fieldcraft_io::read(&rewritten)?;
    if reopened.results != f.results {
        return Err("complete result fields changed after native save/reopen".into());
    }
    fs::write(out.join("roundtrip.fieldcraft"), &rewritten)?;
    fs::write(out.join("complete-fields.json"), serde_json::to_vec_pretty(&fields)?)?;
    fs::write(
        out.join("report.json"),
        serde_json::to_vec_pretty(&json!({
            "input_sha256":fieldcraft_io::sha256_hex(&bytes), "output_sha256":fieldcraft_io::sha256_hex(&rewritten),
            "complete_results_equal":true, "results":fields.len(), "solves":0
        }))?,
    )?;
    Ok(())
}
