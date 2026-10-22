//! Saved views (`view.save`, `view.restore`, `view.delete`, `view.list`): a named camera and
//! what the results display shows, kept in the file's view state ([`fieldcraft_io::views`]) and
//! shown again when the model opens. View state: saved with the model and marking it unsaved,
//! but not a definition, so neither journalled nor undone, and outside the results' definition
//! hash.
//!
//! A view is captured in the parameter shapes of the commands that set it (`ui.camera`'s,
//! `results.contour`'s…), with the unit system its bare numbers are in; showing it runs those
//! commands in turn with their numbers read in that system, and asks the hosts to aim their
//! camera ([`Session::request_camera`]).

use fieldcraft_io::views::{Params, ResultsSetup, SavedView};
use fieldcraft_model::DisplaySystem;
use fieldcraft_model::units::Quantity;
use serde_json::{Value, json};

use crate::results::{Field, Kind, ResultsView};
use crate::view::{Aim, Target};
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};
pub const SAVE: CommandSpec = CommandSpec {
    id: "view.save",
    label: "Save View",
    params: r#"{"name": "Cutaway", "camera"?: {"yaw": -70, "pitch": 18, "distance"?: 300, "target"?: [x, y, z] | {"center_of": "mass" | "bbox", "bodies"?}, "fit"?: true, "fov"?: 0}, "results"?: true} — keeps the camera (ui.camera's parameters; by default the one last aimed in this session: the app's viewport, or the last MCP picture's) and what the results display shows (`results: false` the camera only) as a named view, saved with the model and shown when it opens (it becomes the current view; a view of that name is replaced). Each part is stored in the parameter shape of the command that sets it (results.deformed, results.modes, results.frame, results.contour, results.slice, results.surface, results.streamlines), with `units` naming the system its bare numbers are in; view state, outside the definition and its results' hash"#,
    form: Some(super::forms::VIEW_SAVE),
    kind: CommandKind::Session,
    requires: &[],
    run: save,
};

pub const RESTORE: CommandSpec = CommandSpec {
    id: "view.restore",
    label: "Show Saved View",
    params: r#"{"name": "Cutaway"} — shows a saved view: runs its results commands in turn and aims the camera (the reply's `camera`, in the display units; the app's viewport and MCP pictures follow it); it becomes the current view. A view of a step not solved yet is shown once it is (`pending`); what no longer applies (a field the step has not) comes back in `warnings`"#,
    form: Some(super::forms::VIEW_NAME),
    kind: CommandKind::Session,
    requires: &[],
    run: restore,
};

pub const DELETE: CommandSpec = CommandSpec {
    id: "view.delete",
    label: "Delete Saved View",
    params: r#"{"name": "Cutaway"} — deletes a saved view (refused while a render starts from it, results.render.delete)"#,
    form: Some(super::forms::VIEW_NAME),
    kind: CommandKind::Session,
    requires: &[],
    run: delete,
};

pub const LIST: CommandSpec = CommandSpec {
    id: "view.list",
    label: "Saved Views",
    params: r#"{} — the saved views (each as stored), the current one and the renders"#,
    form: Some(&[]),
    kind: CommandKind::Session,
    requires: &[],
    run: list,
};

/// The `results.*` commands a saved view's results setup runs, in this order.
fn setup_commands(r: &ResultsSetup) -> Vec<(&'static str, Value)> {
    let with_step = |p: &Params| {
        let mut p = p.clone();
        p.insert("step".into(), json!(r.step));
        Value::Object(p)
    };
    let mut out = Vec::new();
    if let Some(p) = &r.deformed {
        out.push(("results.deformed", with_step(p)));
    }
    if let Some(m) = r.mode {
        out.push(("results.modes", json!({"step": r.step, "mode": m})));
    }
    if let Some(f) = r.frame {
        out.push(("results.frame", json!({"step": r.step, "frame": f})));
    }
    for (id, p) in
        [("results.contour", &r.contour), ("results.slice", &r.slice), ("results.surface", &r.surface), ("results.streamlines", &r.streamlines)]
    {
        if let Some(p) = p {
            out.push((id, with_step(p)));
        }
    }
    out
}

/// Run command `id` with `params` whose bare numbers are in `units`, as a script with that
/// header would (the results commands are session commands: nothing is journalled).
pub(crate) fn run_in(s: &mut Session, id: &str, params: &Value, units: &DisplaySystem) -> Result<Value> {
    let spec = crate::find_command(id).ok_or_else(|| EngineError::UnknownCommand(id.into()))?;
    let pinned = s.pinned_input.replace(units.clone());
    let result = (|| {
        let converted = match spec.form {
            Some(form) => crate::units::convert_params(form, params, units, crate::units::si())?,
            None => params.clone(),
        };
        let shown = s.units().clone();
        let mut reply = (spec.run)(s, &converted)?;
        crate::units::convert_reply(id, &mut reply, crate::units::si(), &shown)?;
        Ok(reply)
    })();
    s.pinned_input = pinned;
    result
}

/// What the results display shows now, as the commands that set it take it, in `display`.
fn capture_results(s: &Session, v: &ResultsView) -> Option<ResultsSetup> {
    let (r, _) = s.latest_result(&v.step)?;
    let display = s.units();
    let object = |v: Value| v.as_object().cloned();
    let fluid = matches!(r.kind(), Kind::Fluid(_));
    let mut out = ResultsSetup { step: v.step.clone(), frame: v.frame.map(|f| f as u64), ..ResultsSetup::default() };
    if !fluid {
        out.deformed = object(json!({"show": v.deformed, "scale": v.scale.map_or(json!("auto"), |x| json!(x))}));
    }
    if matches!(r.kind(), Kind::Frequency(_)) {
        out.mode = Some(v.mode as u64);
    }
    let mut contour = json!({"field": v.contour.map_or("none", Field::id), "bands": v.bands});
    if let Some(f) = v.contour
        && (v.range_min.is_some() || v.range_max.is_some())
    {
        let end = |x: Option<f64>| x.map_or(json!("auto"), |x| json!(f.shown(r, x, display)));
        contour["range"] = json!({"min": end(v.range_min), "max": end(v.range_max)});
    }
    out.contour = object(contour);
    if let Kind::Fluid(f) = r.kind() {
        let flow = crate::flow::FlowView::of(v, f);
        let length = |m: f64| crate::results::shown(Quantity::Length, m, display);
        out.slice = object(match &flow.slice {
            Some(sv) => json!({"axis": sv.axis.id(), "position": sv.position.map_or(json!("auto"), |p| json!(length(p)))}),
            None => json!({"axis": "none"}),
        });
        let surface = v.contour.filter(|f| *f == Field::WallHeatFlux).map(|_| json!({"show": flow.surface, "field": "wall_heat_flux"}));
        out.surface = object(surface.unwrap_or_else(|| json!({"show": flow.surface})));
        out.streamlines = object(match &flow.streamlines {
            Some(sl) => {
                let mut p = json!({"show": true, "count": sl.count, "max_steps": sl.max_steps});
                if let Some(seeds) = &sl.seeds {
                    p["seeds"] = json!(seeds.iter().map(|q| q.map(length)).collect::<Vec<_>>());
                }
                p
            }
            None => json!({"show": false}),
        });
    }
    Some(out)
}

/// `aim`'s lengths (distance, a target point) read in `from` and given in `to`.
fn aim_in(aim: &Aim, from: &DisplaySystem, to: &DisplaySystem) -> Aim {
    let length = |x: f64| from.to_si(Quantity::Length, x).and_then(|m| to.express(Quantity::Length, m)).unwrap_or(f64::NAN);
    let mut out = aim.clone();
    out.at.distance = aim.at.distance.map(length).filter(|d| d.is_finite() && *d > 0.0);
    if let Some(Target::Point(p)) = &aim.target {
        out.target = Some(Target::Point(fieldcraft_render::Vec3::new(length(p.x), length(p.y), length(p.z))));
    }
    out
}

fn view_name(p: &Value) -> Result<String> {
    let name = super::name(p)?.to_string();
    fieldcraft_io::views::check_view_name("view", &name).map_err(|e| EngineError::invalid("name", e.to_string()))?;
    Ok(name)
}

fn save(s: &mut Session, p: &Value) -> Result<Value> {
    let name = view_name(p)?;
    let camera = match p.get("camera").filter(|v| !v.is_null()) {
        Some(c) => Some(crate::view::camera(c)?),
        None => s.camera.clone(),
    };
    let with_results = super::opt::<bool>(p, "results")?.unwrap_or(true);
    let results = if with_results { s.results_view.as_ref().and_then(|v| capture_results(s, v)) } else { None };
    if camera.is_none() && results.is_none() {
        return Err(EngineError::Params("nothing to keep: give a `camera`, or show a result first".into()));
    }
    let view = SavedView { name: name.clone(), units: s.units().clone(), camera: camera.and_then(|a| a.json().as_object().cloned()), results };
    let replaced = match s.view.views.iter_mut().find(|v| v.name == name) {
        Some(v) => {
            *v = view.clone();
            true
        }
        None => {
            if s.view.views.len() >= fieldcraft_io::views::MAX_VIEWS {
                return Err(EngineError::Params(format!("a model keeps at most {} views", fieldcraft_io::views::MAX_VIEWS)));
            }
            s.view.views.push(view.clone());
            false
        }
    };
    s.view.current = Some(name);
    s.pending_view = None;
    Ok(json!({"view": view, "replaced": replaced, "views": s.view.views.iter().map(|v| &v.name).collect::<Vec<_>>()}))
}

/// Show saved view `name`: its results commands in turn (once its step is solved), and its
/// camera asked of the hosts. What no longer applies comes back as warnings; `pending` when its
/// step has no result yet.
pub(crate) fn show(s: &mut Session, name: &str) -> Result<Value> {
    let view = s.view.saved(name).cloned().ok_or_else(|| {
        let names: Vec<&str> = s.view.views.iter().map(|v| v.name.as_str()).collect();
        EngineError::Params(format!(
            "no saved view `{}` (views: {})",
            name.chars().take(80).collect::<String>(),
            if names.is_empty() { "none".into() } else { names.join(", ") }
        ))
    })?;
    let mut warnings = Vec::new();
    let mut pending = false;
    if let Some(r) = &view.results {
        if s.latest_result(&r.step).is_some() {
            for (id, params) in setup_commands(r) {
                if let Err(e) = run_in(s, id, &params, &view.units) {
                    warnings.push(format!("view `{name}`: {id}: {e}"));
                }
            }
            s.pending_view = None;
        } else {
            pending = true;
            s.pending_view = Some(name.to_string());
        }
    }
    let mut camera = Value::Null;
    if let Some(c) = &view.camera {
        match crate::view::camera(&Value::Object(c.clone())) {
            Ok(aim) => {
                let aim = aim_in(&aim, &view.units, s.units());
                camera = aim.json();
                s.request_camera(aim);
            }
            Err(e) => warnings.push(format!("view `{name}`: {e}")),
        }
    }
    s.view.current = Some(name.to_string());
    Ok(json!({"view": name, "camera": camera, "pending": pending, "warnings": warnings, "results": s.results_report()}))
}

/// After a step's result is stored: the view waiting for that step is shown.
pub(crate) fn show_pending(s: &mut Session, step: &str) {
    let Some(name) = s.pending_view.clone() else { return };
    if s.view.saved(&name).and_then(|v| v.results.as_ref()).is_some_and(|r| r.step == step) {
        let _ = show(s, &name);
    }
}

fn restore(s: &mut Session, p: &Value) -> Result<Value> {
    let name = super::name(p)?.to_string();
    show(s, &name)
}

fn delete(s: &mut Session, p: &Value) -> Result<Value> {
    let name = super::name(p)?.to_string();
    if let Some(r) = s.view.renders.iter().find(|r| r.view == name) {
        return Err(EngineError::Params(format!("render `{}` starts from view `{name}`: delete it first (results.render.delete)", r.name)));
    }
    let before = s.view.views.len();
    s.view.views.retain(|v| v.name != name);
    if s.view.views.len() == before {
        return Err(EngineError::Params(format!("no saved view `{}`", name.chars().take(80).collect::<String>())));
    }
    if s.view.current.as_deref() == Some(name.as_str()) {
        s.view.current = None;
    }
    if s.pending_view.as_deref() == Some(name.as_str()) {
        s.pending_view = None;
    }
    list(s, &Value::Null)
}

fn list(s: &mut Session, _: &Value) -> Result<Value> {
    Ok(json!({"views": s.view.views, "current": s.view.current, "renders": s.view.renders}))
}

#[cfg(test)]
pub(crate) mod tests {
    use serde_json::json;

    use super::*;

    fn dir(tag: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("fieldcraft-views-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    /// The heated Couette gap solved and shown as a cutaway: temperature clamped, a slice normal
    /// to z, frame 2, a few streamlines.
    pub(crate) fn cutaway() -> Session {
        let mut s = Session::default();
        crate::cmd::solve::tests::heated_couette(&mut s);
        s.execute("solve.fluid", &json!({})).unwrap();
        for (id, p) in [
            ("results.contour", json!({"field": "temperature", "bands": 6, "range": {"min": 0.2, "max": 0.8}})),
            ("results.slice", json!({"axis": "z", "position": 0.06})),
            ("results.frame", json!({"frame": 2})),
            ("results.streamlines", json!({"count": 4, "line": {"from": [0.06, 0.1, 0.06], "to": [0.06, 0.9, 0.06]}})),
        ] {
            s.execute(id, &p).unwrap();
        }
        s
    }

    /// #645: a saved view comes back when the model opens, the same display and the camera asked
    /// of the hosts, without typing commands; it is view state, not an edit.
    #[test]
    fn a_saved_view_comes_back_when_the_model_opens() {
        let d = dir("open");
        let mut s = cutaway();
        let (journal, revision) = (s.journal.len(), s.revision());
        let r = s.execute("view.save", &json!({"name": "Cutaway", "camera": {"yaw": -70, "pitch": 18, "distance": 300}})).unwrap();
        assert_eq!(r["replaced"], false);
        let stored = &r["view"];
        assert_eq!(stored["units"], "mm-t-s-N-MPa");
        assert_eq!(stored["camera"], json!({"yaw": -70.0, "pitch": 18.0, "distance": 300.0}));
        assert_eq!(stored["results"]["step"], "Shear");
        assert_eq!(stored["results"]["frame"], 2);
        assert_eq!(stored["results"]["contour"], json!({"field": "temperature", "bands": 6, "range": {"min": 0.2, "max": 0.8}}));
        assert_eq!(stored["results"]["slice"]["axis"], "z");
        assert_eq!(stored["results"]["streamlines"]["seeds"].as_array().map(Vec::len), Some(4));
        assert_eq!((s.journal.len(), s.revision()), (journal, revision), "not an edit");
        assert!(s.is_dirty(), "view state is saved with the model");
        assert_eq!(s.result_steps().collect::<Vec<_>>(), vec![("Shear", true)], "outside the results' definition");
        let path = d.join("couette.fieldcraft");
        s.execute("file.save_as", &json!({"path": path})).unwrap();

        // A fresh session (shown in SI): the same display, the camera in its units.
        let mut t = Session { default_units: crate::units::si().clone(), ..Session::default() };
        let asked = t.camera_requests();
        let r = t.execute("file.open", &json!({"path": path})).unwrap();
        assert_eq!(r["view"]["name"], "Cutaway", "{r}");
        assert_eq!(r["view"]["pending"], false);
        assert_eq!(r["warnings"], json!([]), "{r}");
        assert_eq!(t.results_view, s.results_view, "the same display");
        assert_eq!(t.camera_requests(), asked + 1);
        let aim = t.camera.clone().unwrap();
        assert!((aim.at.distance.unwrap() - 0.3).abs() < 1e-12 && aim.at.yaw == -70.0, "{aim:?}");
        assert!(!t.is_dirty(), "showing the saved view changes nothing to save");

        // Another view, then back; replace one; delete one.
        t.execute("results.contour", &json!({"field": "velocity"})).unwrap();
        t.execute("view.save", &json!({"name": "Speed"})).unwrap();
        assert_eq!(t.view.current.as_deref(), Some("Speed"));
        t.execute("view.restore", &json!({"name": "Cutaway"})).unwrap();
        assert_eq!(t.results_view, s.results_view);
        assert_eq!(t.execute("view.save", &json!({"name": "Speed", "results": false})).unwrap()["replaced"], true);
        let listed = t.execute("view.list", &json!({})).unwrap();
        assert_eq!(listed["views"].as_array().map(Vec::len), Some(2));
        t.execute("view.delete", &json!({"name": "Speed"})).unwrap();
        assert!(t.view.current.is_none(), "the current view deleted");
        assert!(t.execute("view.restore", &json!({"name": "Speed"})).unwrap_err().to_string().contains("views: Cutaway"));
        let _ = std::fs::remove_dir_all(&d);
    }

    /// A view of a step that has no result yet is shown once it is solved.
    #[test]
    fn a_view_waits_for_its_step_to_be_solved() {
        let d = dir("pending");
        let mut s = cutaway();
        s.execute("view.save", &json!({"name": "Cutaway"})).unwrap();
        let path = d.join("lean.fieldcraft");
        s.execute("file.save_as", &json!({"path": path, "results": false})).unwrap();
        let mut t = Session::default();
        let r = t.execute("file.open", &json!({"path": path})).unwrap();
        assert_eq!(r["view"]["pending"], true, "{r}");
        assert!(t.results_view.is_none());
        t.execute("solve.fluid", &json!({})).unwrap();
        assert_eq!(t.results_view, s.results_view, "shown once solved");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn hostile_view_params_change_nothing() {
        let mut s = cutaway();
        let empty = Session::default();
        let mut e = empty.clone();
        assert!(e.execute("view.save", &json!({"name": "V"})).unwrap_err().to_string().contains("nothing to keep"));
        for p in [
            json!({}),
            json!({"name": ""}),
            json!({"name": " "}),
            json!({"name": "x".repeat(100)}),
            json!({"name": "a\u{0}"}),
            json!({"name": "V", "camera": {"yaw": "left"}}),
            json!({"name": "V", "camera": {"yaw": 0, "pitch": 0, "zoom": 2}}),
            json!({"name": "V", "camera": [1, 2]}),
            json!({"name": "V", "results": "yes"}),
        ] {
            assert!(s.execute("view.save", &p).is_err(), "{p}");
        }
        assert!(s.view.views.is_empty());
        for p in [json!({}), json!({"name": 7}), json!({"name": "Nope"})] {
            assert!(s.execute("view.restore", &p).is_err(), "{p}");
            assert!(s.execute("view.delete", &p).is_err(), "{p}");
        }
        // A stored view whose parameters no longer apply comes back with warnings, not an error.
        s.execute("view.save", &json!({"name": "V"})).unwrap();
        if let Some(r) = s.view.views[0].results.as_mut() {
            r.contour = json!({"field": "von_mises"}).as_object().cloned();
            r.frame = Some(99);
        }
        let r = s.execute("view.restore", &json!({"name": "V"})).unwrap();
        assert_eq!(r["warnings"].as_array().map(Vec::len), Some(2), "{r}");
    }
}
