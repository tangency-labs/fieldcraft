//! `fieldcraft-cli crashes`: the crash reports the desktop app wrote on this computer
//! (docs/crash-reports.md): list them, show one, replay one, delete them.

use std::io::Write;
use std::path::Path;

use fieldcraft_engine::crash;
use serde_json::json;

use crate::Error;

pub const USAGE: &str = "usage: fieldcraft-cli crashes [list [--json] | show <id> | replay <id> | dir | clear [<id>]]";

pub fn run(args: &[String], out: &mut impl Write) -> Result<(), Error> {
    let dir = crash::default_dir().ok_or("no crash folder on this system (set FIELDCRAFT_CRASH_DIR)")?;
    let id = || args.get(1).map(String::as_str).ok_or_else(|| USAGE.to_string());
    match args.first().map(String::as_str) {
        None | Some("list") => {
            let reports = crash::list(&dir);
            if args.iter().any(|a| a == "--json") {
                return print_json(out, &json!({ "dir": dir.to_string_lossy(), "reports": reports }));
            }
            writeln!(out, "crash reports in {}", dir.display())?;
            if reports.is_empty() {
                writeln!(out, "  (none)")?;
            }
            for r in reports {
                let times = if r.count > 1 { format!(" ×{}", r.count) } else { String::new() };
                let cmd = r.command.map(|c| format!(" in {c}")).unwrap_or_default();
                let replay = if r.replayable { "  [replayable]" } else { "" };
                writeln!(out, "  {}  {}{times}{cmd}  v{}{replay}\n      {}", r.id, r.kind, r.version, r.message)?;
            }
            Ok(())
        }
        Some("show") => print_json(out, &crash::read(&dir, id()?)?),
        Some("replay") => replay(&dir, id()?, out),
        Some("dir") => {
            writeln!(out, "{}", dir.display())?;
            Ok(())
        }
        Some("clear") => {
            writeln!(out, "deleted {} report(s)", crash::clear(&dir, args.get(1).map(String::as_str)))?;
            Ok(())
        }
        _ => Err(USAGE.into()),
    }
}

fn print_json(out: &mut impl Write, v: &serde_json::Value) -> Result<(), Error> {
    writeln!(out, "{}", serde_json::to_string_pretty(v).map_err(|e| e.to_string())?)?;
    Ok(())
}

/// Run a report's commands on its model and say whether the crash happens again. Files the
/// commands write go to a scratch folder, never over the user's files.
fn replay(dir: &Path, id: &str, out: &mut impl Write) -> Result<(), Error> {
    let report = crash::read(dir, id)?;
    if report["kind"] == "unclean_exit" {
        writeln!(out, "this crash ended the app without a report: if this process ends abruptly, it reproduced")?;
    }
    writeln!(out, "replaying on {}", report["model_path"].as_str().unwrap_or("a new model"))?;
    let scratch = std::env::temp_dir().join(format!("fieldcraft-replay-{}", std::process::id()));
    std::fs::create_dir_all(&scratch).map_err(|e| format!("{}: {e}", scratch.display()))?;
    let r = crash::replay(&report, &scratch);
    let _ = std::fs::remove_dir_all(&scratch);
    let r = r?;
    for (i, (cmd, outcome)) in r.steps.iter().enumerate() {
        match outcome {
            Ok(()) => writeln!(out, "  {:>3}. {cmd}: ok", i + 1)?,
            Err(e) => writeln!(out, "  {:>3}. {cmd}: {e}", i + 1)?,
        }
    }
    let said = |k: &str| report[k].as_str().unwrap_or_default().to_string();
    if r.reproduced {
        writeln!(out, "reproduced: the same internal error ({}) at {}", said("message"), said("location"))?;
    } else if let Some((l, m)) = r.panics.first() {
        writeln!(out, "an internal error happened, but not the reported one: {m} at {l}")?;
        writeln!(out, "(the report says: {} at {})", said("message"), said("location"))?;
    } else {
        writeln!(out, "did not reproduce: no internal error (fixed, or it depends on something the replay does not have)")?;
    }
    Ok(())
}
