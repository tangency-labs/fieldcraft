use super::*;
use crate::Session;

fn fresh_dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("fieldcraft-crash-test-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    d
}

fn build() -> Build {
    Build { version: "9.9.9", commit: "abc1234" }
}

fn node(label: u32) -> Value {
    json!({"label": label, "xyz": [f64::from(label), 0.0, 0.0]})
}

/// Original 8 x 8 x 1 box: test stop/report semantics without a large sample or CAD.
fn diagnostic_flow(s: &mut Session) {
    s.execute("file.new", &json!({})).unwrap();
    s.execute("model.units", &json!({"system": "m-kg-s-N-Pa"})).unwrap();
    s.execute("model.material.create", &json!({"name": "Fluid", "density": 1, "fluid": {"type": "newtonian", "viscosity": 0.1}})).unwrap();
    s.execute("model.fluid_domain.create", &json!({"name": "Wind", "material": "Fluid", "min": [0,0,0], "max": [1,1,0.125], "spacing": 0.125,
        "faces": {"x_min": {"type": "velocity", "velocity": [0.1,0,0]}, "x_max": {"type": "pressure", "pressure": 0}, "z_min": {"type": "periodic"}, "z_max": {"type": "periodic"}}})).unwrap();
    s.execute("model.step.create", &json!({"name": "Flow", "procedure": {"type": "fluid", "domain": "Wind", "end_time": 0.01, "time_step": 0.001, "frames": 1, "steady": {"window": 1, "min_time": 20}}})).unwrap();
}

/// The child half of the process tests: does nothing unless a parent asks, then installs the
/// hook and crashes the way `FIELDCRAFT_CRASH_CHILD` says.
#[test]
fn child() {
    let (Some(dir), Ok(mode)) = (std::env::var_os("FIELDCRAFT_CRASH_CHILD_DIR"), std::env::var("FIELDCRAFT_CRASH_CHILD")) else {
        return;
    };
    let dir = PathBuf::from(dir);
    if mode == "solver_write_failure" {
        let mut s = Session::default();
        s.execute("file.sample", &json!({"name": "cantilever"})).unwrap();
        let e = s.execute("solve.static", &json!({"step": "Load", "max_dofs": 1})).unwrap_err();
        assert!(matches!(e, EngineError::Solver(fieldcraft_solver::Error::Limit(_))), "{e}");
        return;
    }
    if mode == "solver_error" {
        solver_error(&dir, "solve.fluid", &json!({"max_cells": 1}), &json!({"cells": 8}), "limit reached: 8 cells").unwrap();
    }
    if mode == "replay" {
        // Replay the one report (in a process of its own: the replay replaces the panic hook).
        let id = list(&dir)[0].id.clone();
        let r = replay(&read(&dir, &id).unwrap(), &std::env::temp_dir()).unwrap();
        let steps: Vec<_> = r.steps.iter().map(|(c, r)| json!({"command": c, "ok": r.is_ok()})).collect();
        println!("REPLAY {}", json!({"reproduced": r.reproduced, "nodes": r.session.model.nodes().len(), "steps": steps}));
        return;
    }
    let fresh = install(dir.clone(), build());
    println!("FRESH {}", serde_json::to_string(&fresh).unwrap());
    if mode == "idle" || mode == "solver_error" {
        end_session();
        return;
    }
    // Dying while starting up, before any command: in a graphics driver, with or without what
    // a native crash handler caught.
    if mode == "startup" || mode == "native" {
        phase("window");
        phase("adapter request");
        set_environment("Renderer: Vulkan Test GPU (driver 1.2.3)");
        phase("device creation");
        if mode == "native" {
            let (info, dump) = native_paths().unwrap();
            fs::write(
                info,
                r#"{"code":"0xC0000005","name":"EXCEPTION_ACCESS_VIOLATION","detail":"read of 0x0000000000000010","address":"0x00007FFB12345678","module":"C:\\Windows\\System32\\DriverStore\\igvk64.dll","offset":"0x64EA72","thread":4242}"#,
            )
            .unwrap();
            fs::write(dump, b"MDMP fake").unwrap();
        }
        std::process::abort();
    }
    let mut s = Session { crash_journal: true, ..Default::default() };
    s.execute("mesh.node.add", &node(1)).unwrap();
    s.execute("mesh.node.add", &node(2)).unwrap();
    // Commands that only read are not kept.
    for _ in 0..50 {
        s.execute("model.info", &json!({})).unwrap();
    }
    match mode.as_str() {
        "solver_multiple_steps" => {
            diagnostic_flow(&mut s);
            s.execute(
                "model.step.create",
                &json!({"name": "Second", "procedure": {"type": "fluid", "domain": "Wind", "end_time": 0.02, "time_step": 0.001, "frames": 1, "steady": {"window": 1, "min_time": 20}}}),
            )
            .unwrap();
            s.execute(
                "model.step.create",
                &json!({"name": "Physical", "procedure": {"type": "fluid", "domain": "Wind", "end_time": 0.03, "time_step": 0.001, "frames": 1}}),
            )
            .unwrap();
            let result = s.execute("solve.run", &json!({})).unwrap();
            assert!(result["steps"][0]["summary"]["diagnostic_report"].is_string());
            assert!(result["steps"][1]["summary"]["diagnostic_report"].is_string());
            assert!(result["steps"][2]["summary"].get("diagnostic_report").is_none());
            println!("MULTISTEP {}", result);
        }
        "solver_steady_cap" => {
            diagnostic_flow(&mut s);
            let capped = s.execute("solve.fluid", &json!({"step": "Flow", "backend": "cpu"})).unwrap();
            assert_eq!(capped["stop"]["reason"], "max_time");
            assert!(capped["diagnostic_report"].is_string());
            s.execute(
                "model.step.edit",
                &json!({"name": "Flow", "procedure": {"type": "fluid", "domain": "Wind", "end_time": 0.01, "time_step": 0.001, "frames": 1}}),
            )
            .unwrap();
            let physical = s.execute("solve.fluid", &json!({"step": "Flow", "backend": "cpu"})).unwrap();
            assert!(physical.get("diagnostic_report").is_none());
            println!("TRANSIENT {}", json!({"steady": capped["stop"]["reason"], "cells": capped["cells"], "physical_end": physical["end_time"]}));
        }
        "solver_cap" => {
            s.execute("file.sample", &json!({"name": "cantilever"})).unwrap();
            let e = s.execute("solve.static", &json!({"step": "Load", "max_dofs": 1})).unwrap_err();
            println!("SOLVE {}", json!({"error": e.to_string(), "nodes": s.model.nodes().len()}));
            assert!(matches!(e, EngineError::Context { .. }));
        }
        "solver_not_failure" => {
            s.execute("file.sample", &json!({"name": "cantilever"})).unwrap();
            let e = s.execute("solve.static", &json!({"max_dofs": 0})).unwrap_err();
            assert!(!e.to_string().contains("solver diagnostic"));
            let progress = fieldcraft_model::progress::Progress::new();
            progress.cancel();
            let e = fieldcraft_model::progress::scoped(Some(progress), || s.execute("solve.static", &json!({"step": "Load"}))).unwrap_err();
            assert!(e.is_cancelled());
        }
        "command" => {
            for _ in 0..3 {
                let r = s.execute("diagnostics.crash_test", &json!({"mode": "panic"}));
                assert!(matches!(r, Err(crate::EngineError::Internal { .. })), "{r:?}");
            }
            assert_eq!(s.model.nodes().len(), 2, "the model is unchanged");
        }
        "recovered" => {
            assert_eq!(s.execute("diagnostics.crash_test", &json!({"mode": "recovered"})).unwrap()["caught"], true);
        }
        "thread" => {
            let r = std::thread::Builder::new()
                .name("worker".into())
                .spawn(|| std::panic::panic_any(String::from("worker thread 42 failed")))
                .unwrap()
                .join();
            assert!(r.is_err());
        }
        "abort" => {
            activity("control ui.capture");
            let _ = s.execute("diagnostics.crash_test", &json!({"mode": "abort"}));
        }
        "autosave" => {
            // The journal starts again from each autosave.
            let mut a = crate::autosave::Autosaver::new(&dir.join("recovery")).unwrap();
            assert!(a.save(&s).unwrap());
            s.execute("mesh.node.add", &node(3)).unwrap();
            activity("control ui.capture");
            std::process::abort();
        }
        _ => {}
    }
    end_session();
}

/// Run [`child`] in a new process; returns its stdout and whether it succeeded.
fn run_child(dir: &Path, mode: &str) -> (String, bool) {
    let out = std::process::Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "crash::tests::child", "--test-threads=1", "--nocapture"])
        .env("FIELDCRAFT_CRASH_CHILD_DIR", dir)
        .env("FIELDCRAFT_CRASH_CHILD", mode)
        .env("FIELDCRAFT_CRASH_DIR", dir)
        .env("FIELDCRAFT_CRASH_TEST", "1")
        .env("RUST_BACKTRACE", "0")
        .output()
        .unwrap();
    (String::from_utf8_lossy(&out.stdout).into_owned(), out.status.success())
}

fn printed(out: &str, tag: &str) -> Value {
    let line = out.lines().find_map(|l| l.split_once(tag).map(|(_, rest)| rest)).unwrap_or_else(|| panic!("no {tag} in {out}"));
    serde_json::from_str(line.trim()).unwrap()
}

fn only_report(dir: &Path) -> (Summary, Value) {
    let all = list(dir);
    assert_eq!(all.len(), 1, "{all:?}");
    let s = all.into_iter().next().unwrap();
    let v = read(dir, &s.id).unwrap();
    (s, v)
}

fn nodes_in(model: &str) -> usize {
    let mut s = Session::default();
    s.execute("file.open", &json!({"path": model})).unwrap();
    s.model.nodes().len()
}

#[test]
fn ordinary_solver_records_are_readable_without_a_startup_crash_prompt() {
    let d = fresh_dir("solver_error");
    let (out, ok) = run_child(&d, "solver_error");
    assert!(ok, "{out}");
    assert_eq!(printed(&out, "FRESH "), json!([]));
    let (summary, report) = only_report(&d);
    assert_eq!(summary.kind, "solver_error");
    assert!(!summary.replayable);
    assert_eq!(report["solve"]["cells"], 8);
    assert_eq!(report["command"]["params"]["max_cells"], 1);
    assert!(!d.join(&summary.id).join(MODEL_FILE).exists());
    let path = solver_error(&d, "solve.fluid", &json!({"max_cells": 1}), &json!({"cells": 16}), "limit reached: 16 cells").unwrap();
    assert_eq!(path, d.join(&summary.id).join("report.json"));
    assert_eq!(read(&d, &summary.id).unwrap()["count"], 2);
    assert_eq!(read(&d, &summary.id).unwrap()["solve"]["cells"], 16);
    let other = solver_error(&d, "solve.fluid", &json!({"step": "Other", "max_cells": 1}), &json!({"cells": 16}), "limit reached: 16 cells").unwrap();
    assert_ne!(other, path, "a different step is not a repeat of this failed solve");
    let _ = fs::remove_dir_all(d);
}

#[test]
fn solver_record_redaction_preserves_json_escapes_and_bounds_large_fields() {
    let d = fresh_dir("solver_redaction");
    let home = std::env::var("HOME").or_else(|_| std::env::var("USERPROFILE")).unwrap_or_else(|_| "/example/home".into());
    let input = format!("{home}\\study\\input.json");
    let path = solver_error(&d, "solve.static", &json!({"input": input}), &json!({}), "solver redaction test").unwrap();
    let report: Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
    let expected = if home.len() > 1 && (std::env::var_os("HOME").is_some() || std::env::var_os("USERPROFILE").is_some()) {
        "~\\study\\input.json".to_string()
    } else {
        input
    };
    assert_eq!(report["command"]["params"]["input"], expected);
    let path = solver_error(
        &d,
        "solve.static",
        &json!({"large": "x".repeat(20 << 10)}),
        &json!({"large": "x".repeat(70 << 10)}),
        "oversized solver record test",
    )
    .unwrap();
    let report: Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
    assert!(report["command"]["params"].is_null());
    assert!(report["command"]["params_omitted"].as_str().unwrap().contains("16 KiB"));
    assert!(report["solve"].is_null());
    assert!(report["solve_omitted"].as_str().unwrap().contains("64 KiB"));
    let _ = fs::remove_dir_all(d);
}

#[test]
fn an_actual_solver_cap_records_its_command_and_preserves_the_model() {
    let d = fresh_dir("solver_cap");
    let (out, ok) = run_child(&d, "solver_cap");
    assert!(ok, "{out}");
    let outcome = printed(&out, "SOLVE ");
    assert!(outcome["error"].as_str().unwrap().contains("solver diagnostic saved to"));
    assert_eq!(outcome["nodes"], 315);
    let (_, report) = only_report(&d);
    assert_eq!(report["kind"], "solver_error");
    assert_eq!(report["command"]["command"], "solve.static");
    assert_eq!(report["command"]["params"]["max_dofs"], 1);
    assert_eq!(report["solve"]["requested_step"], "Load");
    assert_eq!(report["solve"]["backend"]["actual"], "cpu");
    assert!(report["solve"]["checkpoint_path"].is_null());
    let _ = fs::remove_dir_all(d);
}

#[test]
fn malformed_parameters_and_cancellation_do_not_create_solver_records() {
    let d = fresh_dir("solver_not_failure");
    let (out, ok) = run_child(&d, "solver_not_failure");
    assert!(ok, "{out}");
    assert!(list(&d).is_empty());
    let _ = fs::remove_dir_all(d);
}

#[test]
fn diagnostic_write_failure_does_not_replace_the_original_solver_error() {
    let d = fresh_dir("solver_write_failure");
    fs::write(&d, b"not a directory").unwrap();
    let (out, ok) = run_child(&d, "solver_write_failure");
    assert!(ok, "{out}");
    assert!(d.is_file());
    let _ = fs::remove_file(d);
}

#[test]
fn a_steady_cap_is_diagnostic_but_a_requested_physical_end_time_is_success() {
    let d = fresh_dir("solver_steady_cap");
    let (out, ok) = run_child(&d, "solver_steady_cap");
    assert!(ok, "{out}");
    assert_eq!(printed(&out, "TRANSIENT ")["steady"], "max_time");
    let (_, report) = only_report(&d);
    assert_eq!(report["solve"]["model"]["cells"], printed(&out, "TRANSIENT ")["cells"]);
    let _ = fs::remove_dir_all(d);
}

#[test]
fn an_earlier_steady_cap_is_recorded_at_its_own_child_boundary() {
    let d = fresh_dir("solver_multiple_steps");
    let (out, ok) = run_child(&d, "solver_multiple_steps");
    assert!(ok, "{out}");
    let result = printed(&out, "MULTISTEP ");
    let reports = list(&d);
    assert_eq!(reports.len(), 2);
    for (index, step) in ["Flow", "Second"].into_iter().enumerate() {
        let report = reports.iter().map(|r| read(&d, &r.id).unwrap()).find(|r| r["solve"]["requested_step"] == step).unwrap();
        assert_eq!(report["command"]["command"], "solve.fluid");
        assert_eq!(report["solve"]["parent_command"], "solve.run");
        assert_eq!(report["solve"]["summary"]["end_time"], result["steps"][index]["summary"]["end_time"]);
        assert_ne!(report["solve"]["summary"]["end_time"], result["steps"][2]["summary"]["end_time"]);
    }
    let _ = fs::remove_dir_all(d);
}

#[test]
fn a_command_panic_is_reported_once_with_the_model_and_the_command_and_replays() {
    let d = fresh_dir("command");
    let (_, ok) = run_child(&d, "command");
    assert!(ok, "the panics were caught");
    let (s, v) = only_report(&d);
    assert_eq!(s.kind, "command_panic");
    assert_eq!(s.count, 3, "repeats are counted, not written again");
    assert_eq!(s.command.as_deref(), Some("diagnostics.crash_test"));
    assert!(s.replayable);
    assert_eq!(v["format"], FORMAT);
    assert_eq!(v["app"]["version"], "9.9.9");
    assert_eq!(v["system"]["os"], std::env::consts::OS);
    assert!(v["message"].as_str().unwrap().contains("an internal error in a command"), "{v}");
    assert!(v["location"].as_str().unwrap().contains("diagnostics.rs"), "{v}");
    let bt = v["backtrace"].as_array().unwrap();
    assert!(bt.iter().any(|f| f.as_str().unwrap().contains("crash_test")), "{bt:?}");
    assert!(!bt.iter().any(|f| f.as_str().unwrap().contains("crash::hook")), "the hook's own frames are dropped: {bt:?}");
    assert_eq!(v["recent_commands"][0]["command"], "mesh.node.add");
    assert_eq!(v["recent_commands"].as_array().unwrap().len(), 3, "read-only commands are not kept: {v}");
    assert_eq!(v["recent_commands"][2]["ok"], false);
    // The model from before the command, and the command: together they reproduce the crash.
    assert_eq!(nodes_in(v["model_path"].as_str().unwrap()), 2);
    let replay: Value = serde_json::from_slice(&fs::read(v["replay_path"].as_str().unwrap()).unwrap()).unwrap();
    assert_eq!(replay["commands"], json!([{"command": "diagnostics.crash_test", "params": {"mode": "panic"}}]));
    assert_eq!(v["reproduce"], format!("fieldcraft-cli crashes replay {}", s.id));
    assert!(fs::read_dir(d.join("sessions")).unwrap().next().is_none(), "a clean exit leaves no journal");
    // The replay reaches the same internal error, from the same model.
    let (out, ok) = run_child(&d, "replay");
    assert!(ok, "{out}");
    let r = printed(&out, "REPLAY ");
    assert_eq!(r["reproduced"], true, "{r}");
    assert_eq!(r["nodes"], 2, "the replay ends in the state the crash left: {r}");
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn a_panic_caught_inside_a_command_is_reported_as_recovered() {
    let d = fresh_dir("recovered");
    let (_, ok) = run_child(&d, "recovered");
    assert!(ok);
    let (s, v) = only_report(&d);
    assert_eq!(s.kind, "recovered_panic");
    assert!(!s.replayable, "the command did not fail");
    assert!(v["message"].as_str().unwrap().contains("caught inside"), "{v}");
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn a_panic_outside_commands_is_reported_from_the_journal() {
    let d = fresh_dir("thread");
    let (_, ok) = run_child(&d, "thread");
    assert!(ok);
    let (s, v) = only_report(&d);
    assert_eq!(s.kind, "panic");
    assert_eq!(v["thread"], "worker");
    assert!(s.message.contains("worker thread 42 failed"));
    // The journal's commands from a new model rebuild it.
    let replay: Value = serde_json::from_slice(&fs::read(v["replay_path"].as_str().unwrap()).unwrap()).unwrap();
    assert_eq!(replay["commands"].as_array().unwrap().len(), 2, "{replay}");
    // A panic the hook reported is no unclean exit at the next start.
    let (out, _) = run_child(&d, "abort");
    let fresh = printed(&out, "FRESH ");
    assert_eq!(fresh.as_array().unwrap().len(), 1, "the earlier report, seen once: {fresh}");
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn an_unclean_exit_is_found_at_the_next_start_with_its_commands() {
    let d = fresh_dir("abort");
    let (_, ok) = run_child(&d, "abort");
    assert!(!ok, "the child aborted");
    assert!(list(&d).is_empty(), "nothing wrote a report while it died");
    assert_eq!(fs::read_dir(d.join("sessions")).unwrap().count(), 1, "its journal is left behind");
    let (out, ok) = run_child(&d, "idle");
    assert!(ok);
    let fresh = printed(&out, "FRESH ");
    assert_eq!(fresh[0]["kind"], "unclean_exit", "{fresh}");
    let (s, v) = only_report(&d);
    assert_eq!(s.command.as_deref(), Some("diagnostics.crash_test"), "the command in flight");
    assert!(s.message.contains("ended without a clean exit while running `diagnostics.crash_test`"), "{}", s.message);
    assert_eq!(v["last_activity"], "control ui.capture");
    let replay: Value = serde_json::from_slice(&fs::read(v["replay_path"].as_str().unwrap()).unwrap()).unwrap();
    let ids: Vec<&str> = replay["commands"].as_array().unwrap().iter().map(|c| c["command"].as_str().unwrap()).collect();
    assert_eq!(ids, ["mesh.node.add", "mesh.node.add", "diagnostics.crash_test"]);
    assert!(fs::read_dir(d.join("sessions")).unwrap().next().is_none(), "the dead journal is gone");
    // Told once: the next start offers nothing new.
    let (out, _) = run_child(&d, "idle");
    assert_eq!(printed(&out, "FRESH "), json!([]));
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn an_unclean_exit_after_an_autosave_starts_from_the_autosave() {
    let d = fresh_dir("autosave");
    let (_, ok) = run_child(&d, "autosave");
    assert!(!ok);
    run_child(&d, "idle");
    let (s, v) = only_report(&d);
    assert_eq!(s.kind, "unclean_exit");
    assert!(s.message.contains("between commands (last: control ui.capture)"), "{}", s.message);
    assert_eq!(nodes_in(v["model_path"].as_str().unwrap()), 2, "the autosaved model");
    let replay: Value = serde_json::from_slice(&fs::read(v["replay_path"].as_str().unwrap()).unwrap()).unwrap();
    assert_eq!(replay["commands"], json!([{"command": "mesh.node.add", "params": node(3)}]), "only what came after the autosave");
    let (out, ok) = run_child(&d, "replay");
    assert!(ok, "{out}");
    assert_eq!(printed(&out, "REPLAY ")["nodes"], 3, "the replay rebuilds the model the app had");
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn journals_parse_and_hostile_ids_are_refused() {
    let text = [
        r#"{"start":"x","pid":1,"version":"1"}"#,
        r#"{"base":"new"}"#,
        r#"{"run":"a","params":{}}"#,
        r#"{"ok":true}"#,
        r#"{"activity":"control ui.capture"}"#,
        r#"{"run":"b","params":{"x":1}}"#,
        r#"{"ok":false,"error":"no"}"#,
        "not json",
        r#"{"run":"c","params":{}}"#,
    ]
    .join("\n");
    let j = parse_journal(&text);
    assert_eq!(j.replay, vec![json!({"command": "a", "params": {}})]);
    assert_eq!(j.recent.len(), 2);
    assert_eq!(j.in_flight, Some(json!({"command": "c", "params": {}})));
    assert_eq!(j.activity.as_deref(), Some("control ui.capture"));
    assert!(j.complete && !j.reported);
    let d = fresh_dir("hostile");
    fs::create_dir_all(&d).unwrap();
    for id in ["", "../x", "a/b", &"f".repeat(65)] {
        assert!(read(&d, id).is_err(), "{id}");
        assert_eq!(clear(&d, Some(id)), 0);
    }
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn signatures_ignore_numbers() {
    assert_eq!(signature("a.rs:1", "node 12 of 400"), signature("a.rs:1", "node 7 of 3"));
    assert_ne!(signature("a.rs:1", "node 12"), signature("a.rs:2", "node 12"));
    assert_eq!(utc(0), "1970-01-01T00:00:00Z");
    assert_eq!(stamp(1_760_436_000), "20251014-100000");
}

#[test]
fn reports_are_pruned_to_the_limits() {
    let d = fresh_dir("prune");
    for i in 0..(MAX_REPORTS + 5) {
        let id = format!("20261014-1000{:02}-{:016x}", i % 60, i + 1);
        fs::create_dir_all(d.join(&id)).unwrap();
        let kind = if i == MAX_REPORTS + 4 { "recovered_panic" } else { "panic" };
        fs::write(d.join(&id).join("report.json"), json!({"kind": kind, "message": "m"}).to_string()).unwrap();
    }
    prune(&d, now());
    let left = list(&d);
    assert_eq!(left.len(), MAX_REPORTS);
    assert!(left.iter().all(|r| r.kind == "panic"), "recovered ones go first");
    assert_eq!(clear(&d, None), MAX_REPORTS);
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn a_crash_while_starting_up_names_the_step_and_the_renderer() {
    let d = fresh_dir("startup");
    let (_, ok) = run_child(&d, "startup");
    assert!(!ok, "the child aborted");
    let (_, ok) = run_child(&d, "idle");
    assert!(ok);
    let (s, v) = only_report(&d);
    assert_eq!(s.kind, "unclean_exit");
    assert!(s.message.contains("during startup (device creation"), "{}", s.message);
    assert_eq!(v["location"], "unclean exit during startup: device creation");
    assert_eq!(v["startup_phase"]["step"], "device creation");
    assert_eq!(v["system"]["environment"], "Renderer: Vulkan Test GPU (driver 1.2.3)", "{v}");
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn a_native_crash_is_reported_with_the_module_the_step_and_the_minidump() {
    let d = fresh_dir("native");
    let (_, ok) = run_child(&d, "native");
    assert!(!ok, "the child aborted");
    let (_, ok) = run_child(&d, "idle");
    assert!(ok);
    let (s, v) = only_report(&d);
    assert_eq!(s.kind, "unclean_exit");
    assert!(s.message.contains("EXCEPTION_ACCESS_VIOLATION (0xC0000005) in igvk64.dll+0x64EA72 during startup (device creation"), "{}", s.message);
    assert_eq!(v["location"], "native exception in igvk64.dll+0x64EA72");
    assert_eq!(v["native"]["thread"], 4242);
    let bt: Vec<&str> = v["backtrace"].as_array().unwrap().iter().map(|f| f.as_str().unwrap()).collect();
    assert_eq!(bt[0], "0xC0000005 EXCEPTION_ACCESS_VIOLATION (read of 0x0000000000000010)");
    assert!(bt[1].contains("igvk64.dll+0x64EA72"), "{bt:?}");
    assert_eq!(fs::read(v["minidump_path"].as_str().unwrap()).unwrap(), b"MDMP fake");
    assert_eq!(v["files"]["minidump"], MINIDUMP_FILE);
    // Nothing is left behind next to the journals.
    let left: Vec<_> = fs::read_dir(d.join("sessions")).unwrap().flatten().map(|e| e.file_name()).collect();
    assert!(left.iter().all(|n| n.to_string_lossy().ends_with(".jsonl")), "{left:?}");
    // The same crash again is counted, not written again.
    let _ = run_child(&d, "native");
    let _ = run_child(&d, "idle");
    let (s, _) = only_report(&d);
    assert_eq!(s.count, 2);
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn native_crash_files_without_their_journal_are_removed() {
    let d = fresh_dir("stray");
    fs::create_dir_all(d.join("sessions")).unwrap();
    fs::write(d.join("sessions/20261015-000000-1-0.native.json"), "{}").unwrap();
    fs::write(d.join("sessions/20261015-000000-1-0.dmp"), "MDMP").unwrap();
    let (_, ok) = run_child(&d, "idle");
    assert!(ok);
    assert!(fs::read_dir(d.join("sessions")).unwrap().next().is_none());
    assert!(list(&d).is_empty());
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn a_hostile_native_file_still_makes_a_bounded_report() {
    // Whatever the filter left (cut short, wrong types, huge strings) is data, never trusted.
    for native in [json!({}), json!({"code": 5, "module": ["x"], "thread": "t"}), json!({"code": "x".repeat(100_000), "module": "y".repeat(100_000)})]
    {
        let p = native_panic(&native, "1.0", "then", "between commands", "idle");
        assert!(p.message.len() < 2_000 && p.backtrace.iter().all(|f| f.len() < 1_000), "{}", p.message.len());
        assert!(p.location.starts_with("native exception in "));
    }
}

#[test]
fn startup_steps_survive_a_rebase_and_running_ends_startup() {
    let d = fresh_dir("rebase");
    fs::create_dir_all(d.join("sessions")).unwrap();
    let mut j = Journal::open(&d, &build(), now());
    j.phase = Some(json!({"phase": "first frame", "ms": 5}).to_string());
    j.environment = Some(json!({"environment": "Renderer: X"}).to_string());
    j.rebase(&json!("new"), 1);
    // Read through the journal's own handle, as `read_live_journal` does (Windows).
    let r = read_journal_file(j.file.as_ref().unwrap().try_clone().unwrap());
    assert_eq!(r.phase, Some(("first frame".to_string(), 5)));
    assert_eq!(r.environment.as_deref(), Some("Renderer: X"));
    j.write(&json!({"phase": RUNNING, "ms": 9}).to_string());
    let r = read_journal_file(j.file.as_ref().unwrap().try_clone().unwrap());
    assert_eq!(r.phase, Some((RUNNING.to_string(), 9)));
    assert!(from_journal(r, None).fields.get("startup_phase").is_none(), "a running app is past startup");
    drop(j);
    let _ = fs::remove_dir_all(&d);
}
