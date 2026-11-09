//! Bounded, non-GUI command replay with one native reply per command.
//! Usage: fatigue_replay recipe.json NEW_OUTPUT_DIRECTORY
use fieldcraft_engine::Session;
use serde_json::{Value, json};
use std::{error::Error, fs, path::PathBuf, time::Instant};
fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1);
    let input = PathBuf::from(args.next().ok_or("missing recipe")?);
    let output = PathBuf::from(args.next().ok_or("missing new output directory")?);
    if args.next().is_some() {
        return Err("unexpected argument".into());
    }
    let bytes = fs::read(&input)?;
    if bytes.len() > 16 * 1024 * 1024 {
        return Err("recipe exceeds 16 MiB".into());
    }
    let recipe: Value = serde_json::from_slice(&bytes)?;
    let rows = recipe.as_array().ok_or("recipe must be array")?;
    if rows.len() > 1000 {
        return Err("recipe exceeds 1000 commands".into());
    }
    fs::create_dir(&output)?;
    fs::write(output.join("recipe.json"), &bytes)?;
    let start = Instant::now();
    let mut s = Session::default();
    let mut probe_nodes: Option<(String, Vec<u32>)> = None;
    for (index, row) in rows.iter().enumerate() {
        let command = row.get("command").and_then(Value::as_str).ok_or("missing command")?;
        let mut params = row.get("params").ok_or("missing params")?.clone();
        if row.get("nodes_from_last_probe").and_then(Value::as_bool) == Some(true) {
            let (step, nodes) = probe_nodes.as_ref().ok_or("no previous resolved probe")?;
            if command != "solve.random" || params.get("step").and_then(Value::as_str) != Some(step.as_str()) {
                return Err("probe binding must target same harmonic step".into());
            }
            params.as_object_mut().ok_or("params not object")?.insert("nodes".into(), json!(nodes));
        }
        fs::write(output.join(format!("{index:03}-resolved-input.json")), serde_json::to_vec_pretty(&params)?)?;
        match s.execute(command, &params) {
            Ok(value) => {
                fs::write(output.join(format!("{index:03}-{command}.json")), serde_json::to_vec_pretty(&value)?)?;
                if command == "solve.probe" {
                    let step = params.get("step").and_then(Value::as_str).ok_or("probe needs step")?.to_string();
                    let rows = value.get("nodes").and_then(Value::as_array).ok_or("probe missing nodes")?;
                    let tolerance = row.get("max_probe_distance_si").and_then(Value::as_f64).ok_or("probe distance tolerance missing")?;
                    if !tolerance.is_finite() || tolerance <= 0.0 {
                        return Err("invalid probe tolerance".into());
                    }
                    let mut labels = Vec::new();
                    for node in rows {
                        let distance = node.get("distance").and_then(Value::as_f64).ok_or("probe distance missing")?;
                        if !distance.is_finite() || distance > tolerance {
                            return Err("probe outside spatial tolerance".into());
                        }
                        let label = node.get("node").and_then(Value::as_u64).ok_or("probe label missing")?;
                        labels.push(u32::try_from(label)?);
                    }
                    probe_nodes = Some((step, labels));
                }
                if command == "model.check" && !model_check_clear(&value) {
                    return Err("model.check reported errors".into());
                }
            }
            Err(e) => {
                fs::write(
                    output.join("failure.json"),
                    serde_json::to_vec_pretty(
                        &json!({"index":index,"command":command,"error":e.to_string(),"seconds":start.elapsed().as_secs_f64()}),
                    )?,
                )?;
                return Err(e.into());
            }
        }
    }
    fs::write(
        output.join("completed.json"),
        serde_json::to_vec_pretty(
            &json!({"commands":rows.len(),"seconds":start.elapsed().as_secs_f64(),"status":"commands_completed_not_engineering_acceptance"}),
        )?,
    )?;
    Ok(())
}

fn model_check_clear(value: &Value) -> bool {
    match value.get("errors") {
        Some(Value::Number(n)) => n.as_u64() == Some(0),
        Some(Value::Array(errors)) => errors.is_empty(),
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn model_check_counts_and_legacy_lists_are_fail_closed() {
        assert!(model_check_clear(&json!({"errors":0})));
        assert!(model_check_clear(&json!({"errors":[]})));
        for value in [json!({"errors":1}), json!({"errors":["bad model"]}), json!({}), json!({"errors":null})] {
            assert!(!model_check_clear(&value));
        }
    }
}
