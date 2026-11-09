//! A cancelled job ([`fieldcraft_model::progress::Progress::cancel`]) stops every solver and
//! mesher promptly with `cancelled` and leaves the model as it was.

use std::time::{Duration, Instant};

use fieldcraft_geometry::sample::Plate;
use fieldcraft_model::progress::{self, Progress};
use serde_json::{Value, json};

use crate::Session;
use crate::sample::impact;

/// How long a cancelled command may keep running: one increment, time step or pass, far
/// below the minutes the runs below would take uncancelled.
const PROMPT: Duration = Duration::from_secs(5);

/// Run `command` on `s` the way a job does and cancel it once it reports the stage `stage`
/// (mid-run), or before it starts (`None`): it must fail with `cancelled`, soon after, leaving
/// the model and the journal as they were.
fn cancel(s: &mut Session, command: &str, params: Value, stage: Option<&str>) {
    cancel_when(s, command, params, stage.map(|stage| move |now: &progress::State| now.stage.starts_with(stage) && now.done > 0.0));
}

/// [`cancel`] once the reported state passes `when` (before it starts with `None`).
fn cancel_when(s: &mut Session, command: &str, params: Value, when: Option<impl Fn(&progress::State) -> bool>) {
    let (model, journal) = (s.model.clone(), s.journal.len());
    let design = s.design.as_ref().map(|d| d.bodies.iter().map(|b| b.name.clone()).collect::<Vec<_>>());
    let p = Progress::new();
    if when.is_none() {
        p.cancel();
    }
    let (r, waited) = std::thread::scope(|scope| {
        let worker = {
            let p = p.clone();
            let s = &mut *s;
            let params = params.clone();
            scope.spawn(move || progress::scoped(Some(p), || s.execute(command, &params)))
        };
        let mut at = Instant::now();
        if let Some(when) = &when {
            let deadline = Instant::now() + Duration::from_secs(600);
            while !worker.is_finished() && Instant::now() < deadline {
                let now = p.snapshot();
                if when(&now) {
                    at = Instant::now();
                    p.cancel();
                    eprintln!("{command}: cancellation requested at {now:?}");
                    break;
                }
                std::thread::sleep(Duration::from_millis(1));
            }
            assert!(p.is_cancelled(), "{command}: never reached the stage ({:?})", p.snapshot());
        }
        let r = worker.join().unwrap();
        (r, at.elapsed())
    });
    eprintln!("{command}: worker joined {waited:?} after cancellation; final progress {:?}", p.snapshot());
    let e = r.expect_err(&format!("{command}: ran to the end although cancelled"));
    assert!(e.is_cancelled(), "{command}: {e}");
    assert_eq!(e.to_string(), "cancelled", "{command}");
    assert!(waited < PROMPT, "{command}: stopped {waited:?} after the cancel");
    assert!(s.model == model, "{command}: the model changed");
    assert_eq!(s.journal.len(), journal, "{command}: journalled");
    assert_eq!(s.design.as_ref().map(|d| d.bodies.iter().map(|b| b.name.clone()).collect::<Vec<_>>()), design, "{command}: the design changed");
}

/// A sample's script run on a fresh session.
fn session(script: &Value) -> Session {
    let mut s = Session::default();
    s.run_script(script).unwrap();
    s
}

fn plate() -> (crate::tests::TempDir, Session) {
    let (dir, path) = crate::tests::sample_design_file(&Plate::default());
    let mut s = Session::default();
    s.execute("geometry.open", &json!({"path": path})).unwrap();
    (dir, s)
}

#[test]
fn explicit_stops_between_time_steps() {
    // A hundred times as long as the sample: hundreds of thousands of increments.
    let mut script = impact();
    for c in script.as_array_mut().unwrap() {
        if c["command"] == "model.step.create" {
            let p = &mut c["params"]["procedure"];
            (p["end_time"], p["output_interval"], p["history_interval"]) = (json!(4e-3), json!(4e-4), json!(2e-5));
        }
    }
    let mut s = session(&script);
    cancel(&mut s, "solve.explicit", json!({}), Some("Time"));
}

#[test]
fn fluid_stops_between_lattice_steps() {
    let mut s = Session::default();
    s.execute("file.sample", &json!({"name": "flow-past"})).unwrap();
    cancel(&mut s, "solve.fluid", json!({}), Some("Fluid steps"));
}

#[test]
fn nonlinear_stops_between_increments() {
    let mut s = Session::default();
    s.execute("file.sample", &json!({"name": "rubber-seal"})).unwrap();
    cancel(&mut s, "solve.nonlinear", json!({}), Some("Load factor"));
}

#[test]
fn transient_heat_stops_between_increments() {
    let mut s = Session::default();
    s.execute("file.sample", &json!({"name": "heated-bracket"})).unwrap();
    s.execute(
        "model.step.create",
        &json!({"name": "Warm", "procedure": {"type": "heat_transfer", "transient": {"end_time": 1e5, "increment": 1}},
            "boundary_conditions": ["Heater"], "loads": ["Air"]}),
    )
    .unwrap();
    cancel(&mut s, "solve.heat", json!({"step": "Warm"}), Some("Time increment"));
}

#[test]
fn static_frequency_and_solve_all_stop_at_the_factorization() {
    let mut s = Session::default();
    s.execute("file.sample", &json!({"name": "cantilever"})).unwrap();
    cancel(&mut s, "solve.static", json!({}), None);
    cancel(&mut s, "solve.frequency", json!({}), None);
    cancel(&mut s, "solve.run", json!({}), None);
}

#[test]
fn meshers_stop_and_leave_no_mesh() {
    let (_dir, mut s) = plate();
    cancel(&mut s, "mesh.tet.generate", json!({"body": "Body1", "size": 0.4}), Some("Refinement round"));
    cancel(&mut s, "mesh.tet.generate", json!({"body": "Body1", "size": 2.5}), None);
    cancel(&mut s, "mesh.hex.generate", json!({"body": "Body1", "size": 2.5}), None);
    cancel(&mut s, "mesh.hex.octree", json!({"body": "Body1", "size": 5.0, "max_depth": 2}), None);
    assert!(s.model.nodes().is_empty());
}

/// The partition, the cells mesher and the plan's apply (#545) stop at their stages: between
/// the kernel's cuts, between cells, between the projection's phases.
#[test]
fn partition_cells_and_plan_apply_stop_and_change_nothing() {
    let (_dir, mut s) = plate();
    s.execute("geometry.datum", &json!({"name": "X30", "by": {"type": "offset", "plane": "yz", "offset": 30}})).unwrap();
    s.execute("geometry.datum", &json!({"name": "Y5", "by": {"type": "offset", "plane": "zx", "offset": 5}})).unwrap();
    let plan = s.execute("mesh.plan", &json!({"body": "Body1"})).unwrap();
    cancel(&mut s, "mesh.plan.apply", json!({"plan": plan["name"]}), None);
    cancel_when(&mut s, "mesh.plan.apply", json!({"plan": plan["name"]}), Some(|now: &progress::State| now.stage.starts_with("Projection: ")));
    cancel(&mut s, "geometry.partition", json!({"name": "Half", "datum": "X30"}), None);
    s.execute("geometry.partition", &json!({"name": "Half", "datum": "X30"})).unwrap();
    // The second partition cuts again by both: cancelled once the first cut is made.
    cancel(&mut s, "geometry.partition", json!({"name": "Strip", "datum": "Y5"}), Some("Partition cuts"));
    let cells = json!({"bodies": ["Body1/a", "Body1/b"], "default": "tet", "size": 0.5});
    cancel(&mut s, "mesh.cells", cells.clone(), None);
    // Inside a cell's tet mesher: the cancel reaches the command, never a fallback.
    cancel(&mut s, "mesh.cells", cells, Some("Refinement round"));
    assert!(s.model.nodes().is_empty());
}
