//! Crash reports and bug report drafts (docs/crash-reports.md): list and read the reports the
//! desktop app wrote, draft a bug report, and (with `FIELDCRAFT_CRASH_TEST=1`) crash on purpose
//! to check that reporting works.

use std::path::PathBuf;

use fieldcraft_io::roots::Access;
use serde_json::{Value, json};

use super::{file_param, opt, req};
use crate::{CommandKind, CommandSpec, EngineError, Result, Session, bug_report, crash};

pub const CRASH_REPORTS: CommandSpec = CommandSpec {
    id: "diagnostics.crash_reports",
    label: "Crash Reports",
    params: r#"{"limit"?: 20, "dir"?: "…"} — local crash reports and ordinary solver diagnostics, newest first: id, kind (command_panic, recovered_panic, panic, unclean_exit, solver_error), time, message, command, count, version, path and whether it can be replayed (default folder: FIELDCRAFT_CRASH_DIR, a temporary folder in runs driven by automation, else the per-user data folder). Solver diagnostics contain no copied model or replay and do not announce an application crash"#,
    form: Some(super::forms::CRASH_REPORTS),
    kind: CommandKind::Session,
    requires: &[],
    run: crash_reports,
};

pub const CRASH_REPORT: CommandSpec = CommandSpec {
    id: "diagnostics.crash_report",
    label: "Crash Report",
    params: r#"{"id": "20261014-101500-…", "dir"?: "…"} — one local report: a crash's backtrace, command, recent commands and attachment paths, or a solver diagnostic's error, command parameters, model size, observed progress and backend context; unavailable fields remain null"#,
    form: Some(super::forms::CRASH_REPORT),
    kind: CommandKind::Session,
    requires: &[],
    run: crash_report,
};

pub const CRASH_REPORTS_DELETE: CommandSpec = CommandSpec {
    id: "diagnostics.crash_reports.delete",
    label: "Delete Crash Reports",
    params: r#"{"id": "20261014-101500-…"} or {"all": true}; "dir"?: "…" — deletes one crash report, or all of them, from this computer: how many were deleted"#,
    form: Some(super::forms::CRASH_REPORTS_DELETE),
    kind: CommandKind::Session,
    requires: &[],
    run: crash_reports_delete,
};

pub const BUG_REPORT: CommandSpec = CommandSpec {
    id: "diagnostics.bug_report",
    label: "Report a Bug",
    params: r#"{"title"?: "…", "observed"?: "what happened", "crash"?: "<crash report id>", "share_session"?: false, "include_model"?: false} — drafts a GitHub bug report: the issue link (the bug form filled with the version, the OS, what happened and, for a crash, where), and the files to attach. Only with share_session does it carry the model's commands, and only with include_model as well the model file. Nothing is sent: the desktop app's Help › Report a Bug shows the draft and opens the link when the user asks"#,
    form: Some(super::forms::BUG_REPORT),
    kind: CommandKind::Session,
    requires: &[],
    run: bug_report_draft,
};

pub const CRASH_TEST: CommandSpec = CommandSpec {
    id: "diagnostics.crash_test",
    label: "Crash Test",
    params: r#"{"mode": "panic|recovered|abort"} — crashes on purpose to check that crash reports work (only with FIELDCRAFT_CRASH_TEST=1): panic, an internal error in this command; recovered, one caught inside it; abort, the app ends at once (reported at its next start)"#,
    form: Some(super::forms::CRASH_TEST),
    kind: CommandKind::Session,
    requires: &[],
    run: crash_test,
};

fn crash_dir(s: &Session, p: &Value) -> Result<PathBuf> {
    crash_dir_for(s, p, Access::Read)
}

fn crash_dir_for(s: &Session, p: &Value, access: Access) -> Result<PathBuf> {
    if p.get("dir").is_some() {
        return Ok(PathBuf::from(file_param(s, p, "dir", access)?));
    }
    crash::default_dir().ok_or_else(|| EngineError::Params("no crash folder: set FIELDCRAFT_CRASH_DIR or give `dir`".into()))
}

fn crash_reports(s: &mut Session, p: &Value) -> Result<Value> {
    let dir = crash_dir(s, p)?;
    let limit = opt::<usize>(p, "limit")?.unwrap_or(crash::MAX_REPORTS);
    let mut reports = crash::list(&dir);
    reports.truncate(limit);
    Ok(json!({"dir": dir.to_string_lossy(), "reports": reports}))
}

fn crash_report(s: &mut Session, p: &Value) -> Result<Value> {
    let dir = crash_dir(s, p)?;
    let id = req::<String>(p, "id")?;
    crash::read(&dir, &id).map_err(EngineError::Params)
}

fn crash_reports_delete(s: &mut Session, p: &Value) -> Result<Value> {
    let id = opt::<String>(p, "id")?;
    let all = opt::<bool>(p, "all")?.unwrap_or(false);
    let dir = crash_dir_for(s, p, Access::Write)?;
    let deleted = match (id.as_deref(), all) {
        (Some(id), false) => match crash::clear(&dir, Some(id)) {
            0 => return Err(EngineError::Params(format!("no crash report `{}`", crash::limited(id, 64)))),
            n => n,
        },
        (None, true) => crash::clear(&dir, None),
        _ => return Err(EngineError::Params("give one of `id` (a crash report) or `all: true`".into())),
    };
    Ok(json!({"deleted": deleted}))
}

/// The draft input for `p` from the session: its commands and model (when shared), or a crash
/// report's replay and model.
pub fn bug_report_input(s: &Session, p: &Value) -> Result<bug_report::Input> {
    let share = opt::<bool>(p, "share_session")?.unwrap_or(false);
    let include_model = opt::<bool>(p, "include_model")?.unwrap_or(false);
    let mut input = bug_report::Input {
        title: opt::<String>(p, "title")?.unwrap_or_default(),
        observed: opt::<String>(p, "observed")?.unwrap_or_default(),
        share_session: share,
        include_model,
        ..Default::default()
    };
    match opt::<String>(p, "crash")? {
        Some(id) => {
            let report = crash::read(&crash_dir(s, p)?, &id).map_err(EngineError::Params)?;
            if share {
                let read = |k: &str| {
                    report
                        .get(k)
                        .and_then(Value::as_str)
                        .and_then(|f| fieldcraft_io::fs::read_capped(std::path::Path::new(f), crash::MAX_MODEL_BYTES as u64).ok())
                };
                let replay: Value = read("replay_path").and_then(|b| serde_json::from_slice(&b).ok()).unwrap_or(Value::Null);
                input.commands = replay.get("commands").and_then(Value::as_array).cloned().unwrap_or_default();
                input.model = if include_model { read("model_path") } else { None };
            }
            input.crash = Some(report);
        }
        None if share => {
            input.commands = s.journal.iter().map(|(command, params)| json!({"command": command, "params": params})).collect();
            input.model = if include_model { Some(s.file_bytes(None)?) } else { None };
        }
        None => {}
    }
    Ok(input)
}

fn bug_report_draft(s: &mut Session, p: &Value) -> Result<Value> {
    let input = bug_report_input(s, p)?;
    Ok(bug_report::draft(&input, &crash::Build::current().version_line()).to_json())
}

// The one deliberate panic and abort: they prove that crash reports work, and run only when
// asked for with FIELDCRAFT_CRASH_TEST=1.
#[allow(clippy::panic)]
fn crash_test(_: &mut Session, p: &Value) -> Result<Value> {
    let mode = req::<String>(p, "mode")?;
    if std::env::var_os("FIELDCRAFT_CRASH_TEST").is_none_or(|v| v != "1") {
        return Err(EngineError::Params("the crash test runs only with FIELDCRAFT_CRASH_TEST=1".into()));
    }
    match mode.as_str() {
        "panic" => std::panic::panic_any(String::from("crash test: an internal error in a command")),
        "recovered" => {
            let caught = std::panic::catch_unwind(|| std::panic::panic_any(String::from("crash test: an internal error caught inside a command")));
            Ok(json!({"caught": caught.is_err()}))
        }
        "abort" => std::process::abort(),
        _ => Err(EngineError::invalid("mode", "must be panic, recovered or abort")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `diagnostics.crash_reports.delete` deletes one report by id, then the rest; an unknown or
    /// hostile id and a missing choice are errors that delete nothing.
    #[test]
    fn crash_reports_are_deleted_one_or_all() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-crash-delete-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let ids = ["20261014-101500-0000000000000001", "20261014-101600-0000000000000002", "20261014-101700-0000000000000003"];
        for id in ids {
            std::fs::create_dir_all(dir.join(id)).unwrap();
            std::fs::write(dir.join(id).join("report.json"), r#"{"kind": "panic", "message": "x"}"#).unwrap();
        }
        let mut s = Session::default();
        let d = dir.to_string_lossy().to_string();
        let listed = |s: &mut Session| s.execute("diagnostics.crash_reports", &json!({"dir": d})).unwrap()["reports"].as_array().map_or(0, Vec::len);
        assert_eq!(listed(&mut s), 3);
        for bad in [
            json!({"dir": d}),
            json!({"dir": d, "id": "../x"}),
            json!({"dir": d, "id": "20261014-101500-00000000000000ff"}),
            json!({"dir": d, "id": ids[0], "all": true}),
        ] {
            assert!(s.execute("diagnostics.crash_reports.delete", &bad).is_err(), "{bad}");
        }
        assert_eq!(listed(&mut s), 3);
        assert_eq!(s.execute("diagnostics.crash_reports.delete", &json!({"dir": d, "id": ids[1]})).unwrap()["deleted"], 1);
        assert_eq!(listed(&mut s), 2);
        assert_eq!(s.execute("diagnostics.crash_reports.delete", &json!({"dir": d, "all": true})).unwrap()["deleted"], 2);
        assert_eq!(listed(&mut s), 0);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
