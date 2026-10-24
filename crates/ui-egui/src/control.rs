//! Programmatic control of the running app (agents, tests). The desktop app reads JSON lines
//! from its socket and hands each request here; replies are `{"ok": true, "result": …}` or
//! `{"ok": false, "error": "…"}`. docs/control-protocol.md lists the methods.
//!
//! Requests are untrusted input: numbers are checked and clamped, and nothing here panics.
//! Adapted from SolidCraft's `crates/ui-egui/src/control.rs` (MIT OR Apache-2.0), see
//! ATTRIBUTION.md.

use std::sync::mpsc::Sender;

use fieldcraft_render::StandardView;
use serde_json::{Value, json};

use crate::capture::CaptureSize;
use crate::{ModelApp, UiState};

/// `ui.render`'s image size when none is given.
const RENDER_DEFAULT: [usize; 2] = [1280, 800];
/// Most files one `ui.drop` drops, and the longest path it takes.
const MAX_DROP: usize = 16;
const MAX_DROP_PATH: usize = 4096;

pub struct ControlRequest {
    pub method: String,
    pub params: Value,
    pub reply: Sender<Value>,
}

impl ControlRequest {
    /// A request whose reply goes to `reply`. The channel is the caller's: the control server's
    /// thread makes it and waits on it, so the window never holds a receiver it could block on.
    pub fn new(method: impl Into<String>, params: Value, reply: Sender<Value>) -> Self {
        Self { method: method.into(), params, reply }
    }
}

/// The key of the interim reply to a request that waits for a background job
/// ([`Outcome::Pending`]): the server then waits for the final reply however long the job takes,
/// since the job always answers (its result, its error, or that it was cancelled).
pub const PENDING: &str = "pending";

pub enum Outcome {
    Done(Value),
    /// The reply comes when the background job finishes: the request's sender waits in
    /// [`crate::files::Loading::waiters`].
    Pending,
    /// A window screenshot, answered once a frame is presented.
    Screenshot {
        path: Option<String>,
    },
}

fn ok(v: Value) -> Outcome {
    Outcome::Done(json!({"ok": true, "result": v}))
}

fn err(e: impl std::fmt::Display) -> Outcome {
    Outcome::Done(json!({"ok": false, "error": e.to_string()}))
}

fn wrap(r: Result<Value, String>) -> Outcome {
    match r {
        Ok(v) => ok(v),
        Err(e) => err(e),
    }
}

fn button(name: Option<&str>) -> egui::PointerButton {
    match name {
        Some("right") => egui::PointerButton::Secondary,
        Some("middle") => egui::PointerButton::Middle,
        _ => egui::PointerButton::Primary,
    }
}

fn key_from(name: &str) -> Option<egui::Key> {
    egui::Key::from_name(name).or(match name.to_ascii_lowercase().as_str() {
        "enter" | "return" => Some(egui::Key::Enter),
        "esc" | "escape" => Some(egui::Key::Escape),
        "delete" => Some(egui::Key::Delete),
        "space" => Some(egui::Key::Space),
        "tab" => Some(egui::Key::Tab),
        _ => None,
    })
}

/// An image side for `ui.render`: absent gives `default`; anything but a whole number of pixels
/// in 1..=MAX_DIM is an error (the renderer also caps the total pixel count).
fn render_dim(p: &Value, key: &str, default: usize) -> Result<usize, String> {
    match p.get(key) {
        None | Some(Value::Null) => Ok(default),
        Some(v) => v
            .as_f64()
            .filter(|x| x.fract() == 0.0 && (1.0..=fieldcraft_render::MAX_DIM as f64).contains(x))
            .map(|x| x as usize)
            .ok_or_else(|| format!("`{key}` must be a whole number of pixels in 1..={}", fieldcraft_render::MAX_DIM)),
    }
}

/// The keys that make `ui.select` change the viewport's selection (without them it selects a
/// definition in the browser, `{kind, name}`).
pub const VIEWPORT_SELECT: [&str; 10] = ["filter", "box_mode", "mode", "clear", "at", "box", "geometry", "nodes", "elements", "faces"];

/// The viewport's selection: the filter, the items, their set members and the hover.
pub fn selection_json(app: &ModelApp) -> Value {
    let s = &app.pick.selection;
    json!({
        "filter": app.pick.filter,
        "box_mode": app.pick.box_mode,
        "count": s.len(),
        "capped": s.capped(),
        "summary": s.summary(),
        "items": s.items().iter().map(crate::selection::Item::to_json).collect::<Vec<_>>(),
        "params": s.params(),
        "hover": app.pick.hover.as_ref().map(crate::selection::Item::to_json),
    })
}

/// The largest window coordinate a request may give.
const MAX_POINT: f64 = 1e6; // unit-factors: a cap in window points

/// The window points in `v`, a list of `n` numbers (two per point).
fn points<const N: usize>(v: &Value, key: &str, what: &str) -> Result<[f32; N], String> {
    let bad = || format!("`{key}` must be {what} in window points");
    let list = v.as_array().filter(|a| a.len() == N).ok_or_else(bad)?;
    let mut out = [0.0; N];
    for (slot, x) in out.iter_mut().zip(list) {
        *slot = x.as_f64().filter(|x| x.is_finite() && x.abs() <= MAX_POINT).ok_or_else(bad)? as f32;
    }
    Ok(out)
}

/// `ui.select` on the viewport: set the filter and box mode, then pick by name and label
/// (`geometry`, `nodes`, `elements`, `faces`), at a point (`at`) or in a box (`box`), changing
/// the selection by `mode`; `clear` empties it first.
fn select(app: &mut ModelApp, p: &Value) -> Result<Value, String> {
    use crate::selection::{Filter, Mode, items_from_params};
    if let Some(f) = p.get("filter") {
        app.pick.filter = serde_json::from_value::<Filter>(f.clone())
            .map_err(|_| format!("`filter` must be one of {}", Filter::ALL.map(|f| json!(f).as_str().unwrap_or_default().to_string()).join(", ")))?;
    }
    if let Some(on) = p.get("box_mode") {
        app.pick.box_mode = on.as_bool().ok_or("`box_mode` must be true or false")?;
    }
    let mode = match p.get("mode") {
        None => Mode::Replace,
        Some(m) => serde_json::from_value::<Mode>(m.clone()).map_err(|_| "`mode` must be replace, add, remove or toggle".to_string())?,
    };
    let at = p.get("at").map(|v| points::<2>(v, "at", "a point [x, y]")).transpose()?;
    let corners = p.get("box").map(|v| points::<4>(v, "box", "[x0, y0, x1, y1]")).transpose()?;
    if corners.is_some() && !app.pick.filter.boxes() {
        return Err(format!("box selection picks nodes, element faces or elements, not {}", app.pick.filter.label().to_lowercase()));
    }
    if (at.is_some() || corners.is_some()) && app.viewport.rect.is_none() {
        return Err("the viewport has not been drawn yet".into());
    }
    let items = items_from_params(p)?;
    if p.get("clear").and_then(Value::as_bool) == Some(true) {
        app.clear_selection();
    }
    if !items.is_empty() {
        app.select(items, mode);
    }
    if let Some([x, y]) = at {
        app.select_at(egui::pos2(x, y), mode);
    }
    if let Some([x0, y0, x1, y1]) = corners {
        app.select_box(egui::pos2(x0, y0), egui::pos2(x1, y1), mode);
    }
    Ok(selection_json(app))
}

/// The path editor's real surface-pick action, also reached from viewport pointer clicks.
/// This changes only the local draft; sampling/export still runs the typed results.path.
fn path_action(app: &mut ModelApp, p: &Value) -> Result<Value, String> {
    let object = p.as_object().ok_or("ui.path parameters must be an object")?;
    if object.keys().any(|k| !["picking", "clear", "at"].contains(&k.as_str())) {
        return Err("ui.path takes picking, clear and at".into());
    }
    let flag = |key: &str| p.get(key).map(|v| v.as_bool().ok_or_else(|| format!("{key} must be true or false"))).transpose();
    let picking = flag("picking")?;
    let clear = flag("clear")?.unwrap_or(false);
    let at = p.get("at").map(|v| points::<2>(v, "at", "a point [x,y]")).transpose()?;
    if at.is_some() && app.viewport.rect.is_none() {
        return Err("the viewport has not been drawn yet".into());
    }
    if app.busy().is_some() {
        return Err("wait for the current job or cancel it before editing the path".into());
    }
    if app.path_editor.is_none() {
        app.open_path_plot();
    }
    let editor = app.path_editor.as_mut().ok_or("a current solid result is needed for a path")?;
    if at.is_some() && !picking.unwrap_or(editor.picking) {
        return Err("enable picking before adding a surface point".into());
    }
    if let Some(picking) = picking {
        editor.picking = picking;
    }
    if clear {
        editor.points.clear();
    }
    let before = editor.points.len();
    if let Some([x, y]) = at {
        app.path_pick_at(egui::pos2(x, y));
    }
    let editor = app.path_editor.as_ref().ok_or("path editor closed")?;
    Ok(json!({"params":editor.params(&app.session)?,"picking":editor.picking,"points":editor.points.len(),"added":editor.points.len()>before}))
}

/// The privacy choices, for `ui.inspect` and `ui.prompt`.
fn privacy_json(app: &ModelApp) -> Value {
    let p = app.prefs();
    json!({
        "available": app.usage_available(),
        "usage_sharing": p.usage_sharing,
        "usage_automation": p.usage_automation,
        "recording": fieldcraft_engine::usage::enabled(),
        "queued": fieldcraft_engine::usage::queued(),
    })
}

/// Answer `req` when the running job finishes; with none running, at once.
fn wait_job(app: &mut ModelApp, req: &ControlRequest) -> Outcome {
    match app.loading.as_mut() {
        Some(l) => {
            l.waiters.push(req.reply.clone());
            Outcome::Pending
        }
        None => ok(json!({"job": null})),
    }
}

pub fn handle(app: &mut ModelApp, ctx: &egui::Context, req: &ControlRequest) -> Outcome {
    let p = &req.params;
    let s = |k: &str| p.get(k).and_then(Value::as_str);
    let f = |k: &str| p.get(k).and_then(Value::as_f64).filter(|x| x.is_finite());
    let b = |k: &str| p.get(k).and_then(Value::as_bool).unwrap_or(false);
    let mods = || egui::Modifiers { shift: b("shift"), ctrl: b("ctrl"), command: b("ctrl"), alt: b("alt"), ..Default::default() };
    // A crash nothing catches (#248) still says which request came last.
    fieldcraft_engine::crash::activity(&format!("control {}", req.method));
    match req.method.as_str() {
        "engine.execute" => {
            let Some(id) = s("command") else { return err("missing `command`") };
            let params = p.get("params").cloned().filter(|v| !v.is_null()).unwrap_or(json!({}));
            // A solve or a mesher runs as a job: the window keeps drawing its progress, and the
            // reply comes when it is done however long that takes (or at once with `"wait":
            // false`, then `job.status` / `job.wait` / `job.cancel`).
            if crate::files::is_long(id) {
                return match app.start_job(id, params) {
                    Err(e) => err(e),
                    Ok(job) if p.get("wait").and_then(Value::as_bool) == Some(false) => ok(json!({"job": job})),
                    Ok(_) => wait_job(app, req),
                };
            }
            wrap(app.run(id, params))
        }
        "job.status" => ok(json!({"job": app.loading.as_ref().map(|l| l.to_json())})),
        "job.wait" => wait_job(app, req),
        "job.cancel" => {
            let job = app.loading.as_ref().map(|l| l.to_json());
            app.cancel_loading();
            ok(json!({"cancelled": job}))
        }
        "engine.script" => {
            let Some(steps) = p.get("commands") else { return err("missing `commands` (a list of {command, params})") };
            wrap(app.run_script(steps).map(|()| json!({"journal": app.session.journal.len()})))
        }
        "engine.commands" => ok(Value::Array(fieldcraft_engine::commands().iter().map(fieldcraft_engine::CommandSpec::listing).collect())),
        "ui.inspect" => ok(json!({
            "ui": app.ui,
            "camera": fieldcraft_engine::view::camera_shown(&app.cam, &app.session),
            "viewport": app.viewport.rect.map(|r| json!([r.left(), r.top(), r.width(), r.height()])),
            "renderer": if app.viewport.gpu.is_some() { "gpu" } else { "cpu" },
            "gpu": app.gpu_info,
            "dialog": app.dialog.as_ref().map(|d| json!({"command": d.command, "text": d.text, "error": d.error})),
            "status": app.status,
            "model": {"nodes": app.session.model.nodes().len(), "elements": app.session.model.elements().len()},
            "results": app.session.results_report(),
            "probe": app.probe,
            "busy": app.busy(),
            "job": app.loading.as_ref().map(|l| l.to_json()),
            "solve": solve_json(app),
            "selected": app.selected.as_ref().map(|(k, n)| json!({"kind": k, "name": n})).or_else(|| app.selected_seed.as_ref().map(|e| json!({"kind": "seed", "edge": e}))),
            "browser": crate::browser::describe(app),
            "nav": app.nav,
            "layout": app.layout.iter().map(|(k, r)| (k.to_string(), json!([r.left(), r.top(), r.width(), r.height()]))).collect::<serde_json::Map<String, Value>>(),
            "prompt": app.prompt.as_ref().map(|p| match p {
                crate::files::Prompt::Unsaved(next) => json!({"unsaved": match next {
                    crate::files::Replace::File(p) => p.clone(),
                    crate::files::Replace::Sample(n) => format!("sample {n}"),
                    crate::files::Replace::Update => "update".to_string(),
                    crate::files::Replace::New => "new".to_string(),
                }}),
                crate::files::Prompt::Recover => json!("recover"),
                crate::files::Prompt::Samples => json!("samples"),
                crate::files::Prompt::Design(p) => json!({"design": p}),
                crate::files::Prompt::ReplaceMesh(r) => json!({"replace_mesh": r.question, "warnings": r.warnings}),
                crate::files::Prompt::Delete(d) => json!({"delete": d.question(), "kind": d.kind, "name": d.name, "users": d.users, "allowed": d.allowed()}),
                crate::files::Prompt::Crashed(r) => json!({"crashed": r}),
                crate::files::Prompt::CrashReports(r) => json!({"crash_reports": r}),
                crate::files::Prompt::ReportBug => json!({"report_bug": app.report.as_ref().and_then(|f| f.drafted_json())}),
                crate::files::Prompt::Consent => json!({"consent": crate::privacy::QUESTION}),
                crate::files::Prompt::Settings => json!("settings"),
                crate::files::Prompt::Update => json!({"update": app.update_json()}),
            }),
            "privacy": privacy_json(app),
            "updates": app.update_json(),
            "recent": app.recent.paths,
            "version": crate::VERSION,
            "path_plot":app.path_editor.as_ref().map(|e|json!({"params":e.params(&app.session).ok(),"picking":e.picking,"selected":e.selected})),
        })),
        "ui.toolbar" => ok(crate::toolbar::describe(app)),
        "ui.path" => wrap(path_action(app, p)),
        "ui.settings" => {
            app.open_settings();
            ok(privacy_json(app))
        }
        "ui.select" if p.get("kind").is_none() && VIEWPORT_SELECT.iter().any(|k| p.get(*k).is_some()) => wrap(select(app, p)),
        "ui.selection" => ok(selection_json(app)),
        "ui.select" => {
            // Select a definition as a browser click does (its properties show on the right).
            let kind = p.get("kind").cloned().and_then(|k| serde_json::from_value::<crate::browser::Def>(k).ok());
            match (kind, s("name")) {
                (Some(kind), Some(name)) if crate::browser::definition(&app.session.model, kind, name).is_some() => {
                    app.selected = Some((kind, name.to_string()));
                    ok(Value::Null)
                }
                (Some(_), Some(name)) => err(format!("no such definition `{name}`")),
                _ => {
                    app.selected = None;
                    ok(Value::Null)
                }
            }
        }
        "ui.set" => {
            let mut cur = serde_json::to_value(&app.ui).unwrap_or(json!({}));
            if let (Some(o), Some(src)) = (cur.as_object_mut(), p.as_object()) {
                for (k, v) in src {
                    if !o.contains_key(k) {
                        return err(format!("unknown setting `{k}`"));
                    }
                    o.insert(k.clone(), v.clone());
                }
            }
            match serde_json::from_value::<UiState>(cur) {
                Ok(mut u) => {
                    // The workspace and its tab stay together: a tab brings its workspace, a
                    // workspace alone its first tab.
                    if p.get("tab").is_some() {
                        u.select_tab(u.tab);
                    } else if p.get("workspace").is_some() {
                        u.select_workspace(u.workspace);
                    }
                    app.ui = u;
                    ok(json!(app.ui))
                }
                Err(e) => err(e),
            }
        }
        "ui.view" => {
            match s("view").unwrap_or("iso") {
                "fit" => app.fit_view(),
                "home" => {
                    app.home_view();
                }
                v => match StandardView::parse(v) {
                    Some(sv) => app.set_view(sv),
                    None => return err(format!("unknown view `{v}` (front, back, left, right, top, bottom, iso, home, fit)")),
                },
            }
            ok(fieldcraft_engine::view::camera_shown(&app.cam, &app.session))
        }
        "ui.camera" => match fieldcraft_engine::view::camera(p).map_err(|e| e.to_string()).and_then(|a| app.aim_camera(&a)) {
            Ok(()) => ok(fieldcraft_engine::view::camera_shown(&app.cam, &app.session)),
            Err(e) => err(e),
        },
        "ui.start" => {
            let Some(id) = s("command") else { return err("missing `command`") };
            match app.start_programmatic(id) {
                Ok(()) => ok(json!({"dialog": app.dialog.as_ref().map(|d| d.command)})),
                Err(e) => err(e),
            }
        }
        "ui.dialog" => {
            if app.dialog.is_none() {
                return err("no command form is open (ui.start opens one)");
            }
            if b("cancel") {
                app.close_dialog();
                return ok(Value::Null);
            }
            if let (Some(d), Some(params)) = (app.dialog.as_mut(), p.get("params")) {
                d.text = serde_json::to_string_pretty(params).unwrap_or_default();
            }
            // The form's Target row: the viewport's selection into the parameters.
            if b("use_selection")
                && let Err(e) = app.use_selection_in_dialog()
            {
                return err(e);
            }
            if b("run") {
                return wrap(app.run_dialog());
            }
            // The form's OK button, as a click runs it: a mesher that would replace a mesh asks
            // first (`ui.inspect`'s `prompt`, answered with `ui.prompt`).
            if b("submit") {
                app.submit_dialog();
                return ok(json!({"prompt": app.prompt.is_some(), "busy": app.busy()}));
            }
            ok(json!(app.dialog.as_ref().map(|d| json!({"command": d.command, "text": d.text, "error": d.error}))))
        }
        "ui.item" => {
            // An item menu's entry, as a right click on the item's browser row and a click on the
            // entry do (crate::menu).
            let kind = s("kind");
            let name = s("name").unwrap_or_default().to_string();
            let target = match kind {
                Some("body") => crate::menu::Target::Body(name),
                Some("surface") => crate::menu::Target::Surface(name),
                Some(k) => match serde_json::from_value::<crate::browser::Def>(json!(k)) {
                    Ok(def) => crate::menu::Target::Def(def, name),
                    Err(_) => return err(format!("unknown item kind `{}`", k.chars().take(40).collect::<String>())),
                },
                None => return err("missing `kind`"),
            };
            let action = match (s("action"), s("step")) {
                (Some("edit"), _) => crate::menu::Action::Edit,
                (Some("rename"), _) => crate::menu::Action::Rename,
                (Some("suppress"), _) => crate::menu::Action::Suppress(true),
                (Some("unsuppress"), _) => crate::menu::Action::Suppress(false),
                (Some("hide"), _) => crate::menu::Action::Hide,
                (Some("show"), _) => crate::menu::Action::Show,
                (Some("isolate"), _) => crate::menu::Action::Isolate,
                (Some("show_all"), _) => crate::menu::Action::ShowAll,
                (Some("delete"), _) => crate::menu::Action::Delete,
                (Some(a @ ("activate" | "deactivate")), Some(step)) => {
                    crate::menu::Action::Active { step: step.to_string(), active: a == "activate" }
                }
                (a, _) => {
                    return err(format!(
                        "`action` {a:?}: edit, rename, suppress, unsuppress, hide, show, isolate, show_all, delete, or activate/deactivate with `step`"
                    ));
                }
            };
            crate::menu::apply(app, target, action);
            ok(json!({"status": app.status, "prompt": app.prompt.is_some(), "dialog": app.dialog.as_ref().map(|d| d.command)}))
        }
        "ui.prompt" => {
            let Some(prompt) = app.prompt.clone() else { return err("no prompt is open") };
            match (prompt, s("answer")) {
                (crate::files::Prompt::Unsaved(_), Some(a @ ("save" | "discard" | "cancel"))) => {
                    let answer = match a {
                        "save" => crate::files::Answer::Save,
                        "discard" => crate::files::Answer::Discard,
                        _ => crate::files::Answer::Cancel,
                    };
                    app.answer_unsaved(answer);
                    ok(json!({"busy": app.busy()}))
                }
                (crate::files::Prompt::Design(_), Some(a @ ("new" | "keep" | "cancel"))) => {
                    let answer = match a {
                        "new" => crate::files::DesignAnswer::New,
                        "keep" => crate::files::DesignAnswer::Keep,
                        _ => crate::files::DesignAnswer::Cancel,
                    };
                    app.answer_design(answer);
                    ok(json!({"busy": app.busy()}))
                }
                (crate::files::Prompt::ReplaceMesh(_), Some(a @ ("replace" | "cancel"))) => {
                    app.answer_replace_mesh(a == "replace");
                    ok(json!({"busy": app.busy()}))
                }
                (crate::files::Prompt::Delete(_), Some(a @ ("delete" | "cancel"))) => {
                    app.answer_delete(a == "delete");
                    ok(json!({"status": app.status}))
                }
                (crate::files::Prompt::Samples, Some("open")) => match s("name") {
                    Some(name) => {
                        app.prompt = None;
                        app.open_sample(name);
                        ok(json!({"prompt": app.prompt.is_some()}))
                    }
                    None => err("missing `name` (a sample)"),
                },
                (crate::files::Prompt::Recover | crate::files::Prompt::Samples, Some("close" | "cancel")) => {
                    app.prompt = None;
                    ok(Value::Null)
                }
                (crate::files::Prompt::Crashed(r), Some("report")) => {
                    app.report_bug(r.first().map(|c| c.id.clone()));
                    ok(json!({"prompt": "report_bug"}))
                }
                (crate::files::Prompt::CrashReports(_), Some("report")) => match s("id") {
                    Some(id) => {
                        app.report_bug(Some(id.to_string()));
                        ok(json!({"prompt": "report_bug"}))
                    }
                    None => err("missing `id` (a crash report)"),
                },
                (crate::files::Prompt::CrashReports(_), Some("delete")) => match s("id") {
                    Some(id) => wrap(app.delete_crash_reports(Some(id))),
                    None => err("missing `id` (a crash report)"),
                },
                (crate::files::Prompt::CrashReports(_), Some("delete_all")) => wrap(app.delete_crash_reports(None)),
                (crate::files::Prompt::Crashed(_) | crate::files::Prompt::CrashReports(_), Some("close" | "cancel")) => {
                    app.prompt = None;
                    ok(Value::Null)
                }
                // The form's fields; sending opens the browser, which only a click does.
                (crate::files::Prompt::ReportBug, Some("set")) => wrap(app.set_report(p)),
                (crate::files::Prompt::Consent, Some(a @ ("share" | "decline"))) => {
                    app.answer_consent(a == "share");
                    ok(privacy_json(app))
                }
                (crate::files::Prompt::Settings, Some("set")) if p.get("units").is_some() => {
                    match p.get("units").cloned().map(serde_json::from_value::<fieldcraft_model::DisplaySystem>) {
                        Some(Ok(units)) => {
                            app.set_default_units(units);
                            ok(json!({"units": app.prefs().units}))
                        }
                        _ => err("`units` must be a unit system (mm-t-s-N-MPa, m-kg-s-N-Pa, mm-kg-ms-kN-GPa, mm-g-ms-N-MPa, in-lbf-s-psi)"),
                    }
                }
                (crate::files::Prompt::Settings, Some("set")) if p.get("update_check").is_some() || p.get("update_skip").is_some() => {
                    let automatic = p.get("update_check").and_then(Value::as_bool).unwrap_or(app.prefs().update_check);
                    let skip = p.get("update_skip").and_then(Value::as_str).map(str::to_string);
                    app.set_update_prefs(automatic, skip);
                    ok(app.update_json())
                }
                (crate::files::Prompt::Update, Some(a @ ("restart" | "later"))) => {
                    app.answer_update(a == "restart");
                    ok(app.update_json())
                }
                (crate::files::Prompt::Settings, Some("set")) => {
                    if !app.usage_available() {
                        return err("this run cannot share usage (no analytics key, or FIELDCRAFT_TELEMETRY=0, CI or automation)");
                    }
                    let sharing = p.get("usage_sharing").and_then(Value::as_bool).unwrap_or(app.prefs().usage_sharing == Some(true));
                    let automation = p.get("usage_automation").and_then(Value::as_bool).unwrap_or(app.prefs().usage_automation);
                    app.set_usage(sharing, automation);
                    ok(privacy_json(app))
                }
                (crate::files::Prompt::Settings, Some("close" | "cancel")) => {
                    app.prompt = None;
                    ok(Value::Null)
                }
                (crate::files::Prompt::ReportBug, Some("close" | "cancel")) => {
                    app.prompt = None;
                    app.report = None;
                    ok(Value::Null)
                }
                (crate::files::Prompt::Recover, Some("recover")) => match s("id") {
                    Some(id) => {
                        app.recover(id);
                        ok(json!({"busy": app.busy()}))
                    }
                    None => err("missing `id` (from file.recovery.list)"),
                },
                (_, a) => err(format!(
                    "`answer` {a:?} does not answer this prompt (unsaved: save, discard, cancel; design: new, keep, cancel; replace_mesh: replace, cancel; delete: delete, cancel; recover: recover with id, close; samples: open with name, close; crashed: report, close; report_bug: set with title/observed/share_session/include_model, close; consent: share, decline; settings: set with usage_sharing/usage_automation or update_check/update_skip, close; update: restart, later)"
                )),
            }
        }
        "ui.move" | "ui.click" => {
            let (Some(x), Some(y)) = (f("x"), f("y")) else { return err("missing x/y") };
            let pos = egui::pos2(x as f32, y as f32);
            app.synthetic.push(egui::Event::PointerMoved(pos));
            if req.method == "ui.click" {
                let (button, modifiers) = (button(s("button")), mods());
                // `double: true` sends two clicks in one go, which egui reads as a double click.
                for _ in 0..if b("double") { 2 } else { 1 } {
                    app.synthetic.push(egui::Event::PointerButton { pos, button, pressed: true, modifiers });
                    app.synthetic.push(egui::Event::PointerButton { pos, button, pressed: false, modifiers });
                }
            }
            ok(Value::Null)
        }
        "ui.drag" => {
            let (Some(x0), Some(y0), Some(x1), Some(y1)) = (f("x0"), f("y0"), f("x1"), f("y1")) else { return err("missing x0/y0/x1/y1") };
            let (button, modifiers) = (button(s("button")), mods());
            let (a, z) = (egui::pos2(x0 as f32, y0 as f32), egui::pos2(x1 as f32, y1 as f32));
            app.synthetic.push(egui::Event::PointerMoved(a));
            app.synthetic.push(egui::Event::PointerButton { pos: a, button, pressed: true, modifiers });
            let steps = f("steps").unwrap_or(10.0).clamp(1.0, 200.0) as usize;
            for k in 1..=steps {
                app.synthetic.push(egui::Event::PointerMoved(a + (z - a) * (k as f32 / steps as f32)));
            }
            app.synthetic.push(egui::Event::PointerButton { pos: z, button, pressed: false, modifiers });
            ok(Value::Null)
        }
        "ui.scroll" => {
            let (Some(x), Some(y)) = (f("x"), f("y")) else { return err("missing x/y") };
            app.synthetic.push(egui::Event::PointerMoved(egui::pos2(x as f32, y as f32)));
            // Wheel lines; negative zooms out, positive zooms in.
            let lines = f("delta").unwrap_or(1.0).clamp(-100.0, 100.0);
            app.synthetic.push(egui::Event::MouseWheel {
                unit: egui::MouseWheelUnit::Line,
                delta: egui::vec2(0.0, lines as f32),
                modifiers: Default::default(),
                phase: egui::TouchPhase::Move,
            });
            ok(Value::Null)
        }
        "ui.key" => {
            let Some(key) = s("key").and_then(key_from) else { return err("unknown key") };
            let modifiers = egui::Modifiers { shift: b("shift"), command: b("cmd"), ctrl: b("cmd") || b("ctrl"), ..Default::default() };
            app.synthetic.push(egui::Event::Key { key, physical_key: None, pressed: true, repeat: false, modifiers });
            app.synthetic.push(egui::Event::Key { key, physical_key: None, pressed: false, repeat: false, modifiers });
            ok(Value::Null)
        }
        "ui.text" => {
            app.synthetic.push(egui::Event::Text(s("text").unwrap_or("").to_string()));
            ok(Value::Null)
        }
        "ui.drop" => {
            // Files dropped on the window: the next frame takes them with the window system's
            // drops, and they open as a picked file opens (the native picker cannot be driven).
            let paths: Option<Vec<String>> =
                p.get("paths").and_then(Value::as_array).and_then(|a| a.iter().map(|v| v.as_str().map(str::to_string)).collect());
            match paths {
                Some(paths) if !paths.is_empty() && paths.len() <= MAX_DROP && paths.iter().all(|p| !p.is_empty() && p.len() <= MAX_DROP_PATH) => {
                    let n = paths.len();
                    app.synthetic_drops.extend(paths);
                    ok(json!({"dropped": n}))
                }
                _ => err(format!("`paths` must be a list of 1 to {MAX_DROP} file paths")),
            }
        }
        "ui.worldToScreen" => {
            let pts: Vec<Value> = p
                .get("points")
                .and_then(Value::as_array)
                .into_iter()
                .flatten()
                .take(1000)
                .map(|v| {
                    let c = |i: usize| v.get(i).and_then(Value::as_f64).unwrap_or(f64::NAN);
                    // Given in the display units; the scene is in metres.
                    crate::viewport::to_screen(app, app.scene_point([c(0), c(1), c(2)])).map_or(Value::Null, |q| json!([q.x, q.y]))
                })
                .collect();
            ok(json!(pts))
        }
        "ui.resize" => {
            let (Some(w), Some(h)) = (f("width"), f("height")) else { return err("missing width/height") };
            // A size is a restored window's: one that started maximized on a small monitor too.
            ctx.send_viewport_cmd(egui::ViewportCommand::Maximized(false));
            ctx.send_viewport_cmd(egui::ViewportCommand::InnerSize(egui::vec2(w.clamp(320.0, 8192.0) as f32, h.clamp(240.0, 8192.0) as f32)));
            ok(Value::Null)
        }
        "ui.fullscreen" => {
            // Into or out of full screen (`on`, default true); the title bar follows next frame.
            let on = p.get("on").and_then(Value::as_bool).unwrap_or(true);
            ctx.send_viewport_cmd(egui::ViewportCommand::Fullscreen(on));
            ok(json!({"fullscreen": on}))
        }
        "ui.screenshot" => Outcome::Screenshot { path: s("path").map(str::to_string) },
        "ui.capture" => {
            let live = window_size(ctx);
            let dim = |k: &str, default: f32| match p.get(k) {
                None | Some(Value::Null) => Ok(default),
                Some(v) => v.as_f64().filter(|x| x.is_finite()).map(|x| x as f32).ok_or_else(|| format!("`{k}` must be a number")),
            };
            let size = match (dim("width", live.width), dim("height", live.height), dim("scale", live.scale)) {
                (Ok(width), Ok(height), Ok(scale)) => CaptureSize { width, height, scale },
                (Err(e), ..) | (_, Err(e), _) | (.., Err(e)) => return err(e),
            };
            Outcome::Done(capture_to(app, size, s("path")))
        }
        "ui.render" => {
            let Some(path) = s("path") else { return err("missing `path`") };
            let (w, h) = match (render_dim(p, "width", RENDER_DEFAULT[0]), render_dim(p, "height", RENDER_DEFAULT[1])) {
                (Ok(w), Ok(h)) => (w, h),
                (Err(e), _) | (_, Err(e)) => return err(e),
            };
            // Sizes on screen are in points: at the window's pixels per point unless given, so a
            // render at device resolution matches the window (#285).
            let scale =
                match p.get("scale") {
                    None | Some(Value::Null) => Ok(window_size(ctx).scale),
                    Some(v) => v.as_f64().map(|x| x as f32).filter(|x| crate::capture::SCALE_RANGE.contains(x)).ok_or_else(|| {
                        format!("`scale` must be from {} to {}", crate::capture::SCALE_RANGE.start(), crate::capture::SCALE_RANGE.end())
                    }),
                };
            let scale = match scale {
                Ok(s) => s,
                Err(e) => return err(e),
            };
            let (png, renderer) = match app.render_png(w, h, f64::from(scale)) {
                Ok(drawn) => drawn,
                Err(e) => return err(e),
            };
            match fieldcraft_engine::write_output(std::path::Path::new(path), &png) {
                Ok(()) => ok(json!({"path": path, "width": w, "height": h, "scale": scale, "bytes": png.len(), "renderer": renderer.report()})),
                Err(e) => err(format!("{path}: {e}")),
            }
        }
        "app.quit" => {
            app.quit_requested = true;
            ctx.send_viewport_cmd(egui::ViewportCommand::Close);
            ok(Value::Null)
        }
        other => err(format!("unknown method `{other}`")),
    }
}

/// The window's size in points and its pixels per point (1400 × 900 at 1 before a window was
/// laid out).
pub fn window_size(ctx: &egui::Context) -> CaptureSize {
    let r = ctx.content_rect();
    let fits = |v: f32| v.is_finite() && (1.0..=crate::capture::MAX_SIDE as f32).contains(&v);
    let scale = Some(ctx.pixels_per_point()).filter(|s| crate::capture::SCALE_RANGE.contains(s)).unwrap_or(1.0);
    if fits(r.width()) && fits(r.height()) {
        CaptureSize { width: r.width(), height: r.height(), scale }
    } else {
        CaptureSize { width: 1400.0, height: 900.0, scale }
    }
}

/// Capture the whole window offscreen ([`ModelApp::capture`]) and save it as PNG (a reply).
pub fn capture_to(app: &mut ModelApp, size: CaptureSize, path: Option<&str>) -> Value {
    match app.capture(size) {
        Ok(image) => {
            let mut reply = save_screenshot(&image, path);
            if let Some(r) = reply.get_mut("result").and_then(Value::as_object_mut) {
                r.insert("offscreen".into(), json!(true));
                r.insert("scale".into(), json!(size.scale));
            }
            reply
        }
        Err(e) => json!({"ok": false, "error": e}),
    }
}

/// Save a screenshot as PNG: at `path`, or `fieldcraft-shot.png` in the temporary folder.
pub fn save_screenshot(image: &egui::ColorImage, path: Option<&str>) -> Value {
    let [w, h] = image.size;
    let rgba: Vec<u8> = image.pixels.iter().flat_map(|c| c.to_array()).collect();
    let Some(png) = fieldcraft_render::png(w, h, &rgba) else { return json!({"ok": false, "error": "could not encode the screenshot"}) };
    let path = path.map(str::to_string).unwrap_or_else(|| std::env::temp_dir().join("fieldcraft-shot.png").to_string_lossy().to_string());
    match fieldcraft_engine::write_output(std::path::Path::new(&path), &png) {
        Ok(()) => json!({"ok": true, "result": {"path": path, "width": w, "height": h}}),
        Err(e) => json!({"ok": false, "error": format!("{path}: {e}")}),
    }
}

/// The title bar's Solve as `ui.inspect` reports it.
fn solve_json(app: &ModelApp) -> Value {
    let missing = app.unavailable("solve.run");
    json!({"enabled": app.busy().is_none() && missing.is_none(), "unavailable": missing.map(|u| u.to_json())})
}

#[cfg(test)]
mod tests {
    use fieldcraft_engine::Session;

    use super::*;

    fn call(app: &mut ModelApp, method: &str, params: Value) -> Value {
        let (tx, rx) = std::sync::mpsc::channel();
        let req = ControlRequest::new(method, params, tx);
        match handle(app, &egui::Context::default(), &req) {
            Outcome::Done(v) => v,
            Outcome::Screenshot { .. } => json!({"screenshot": true}),
            // A job's reply comes when it ends.
            Outcome::Pending => {
                app.wait_loaded();
                rx.try_recv().unwrap_or(Value::Null)
            }
        }
    }

    #[test]
    fn engine_view_and_settings_methods() {
        let mut app = ModelApp::new(Session::default());
        let r = call(&mut app, "engine.script", json!({"commands": fieldcraft_engine::sample::script()}));
        assert_eq!(r["ok"], true, "{r}");
        let r = call(&mut app, "engine.execute", json!({"command": "model.info"}));
        assert_eq!(r["result"]["elements"]["hex8"], 64);
        assert_eq!(call(&mut app, "engine.execute", json!({"command": "no.such"}))["ok"], false);
        assert!(call(&mut app, "engine.commands", json!({}))["result"].as_array().is_some_and(|a| a.len() == fieldcraft_engine::commands().len()));
        let r = call(&mut app, "ui.set", json!({"style": "wireframe", "show_nodes": true}));
        assert_eq!(r["result"]["style"], "wireframe");
        assert_eq!(app.ui.style, crate::Style::Wireframe);
        assert_eq!(call(&mut app, "ui.set", json!({"style": "cartoon"}))["ok"], false);
        assert_eq!(call(&mut app, "ui.set", json!({"nope": 1}))["ok"], false);
        assert_eq!(call(&mut app, "ui.view", json!({"view": "top"}))["ok"], true);
        assert!((app.cam.pitch - std::f64::consts::FRAC_PI_2).abs() < 1e-3);
        assert_eq!(call(&mut app, "ui.view", json!({"view": "sideways"}))["ok"], false);
        // Any orientation in degrees; hostile numbers are refused or clamped and leave a valid camera.
        let r = call(&mut app, "ui.camera", json!({"yaw": 90, "pitch": 30, "distance": 250}));
        assert!((app.cam.yaw - std::f64::consts::FRAC_PI_2).abs() < 1e-9 && (app.cam.pitch - 30f64.to_radians()).abs() < 1e-9, "{r}");
        assert_eq!(r["result"]["distance"], 250.0);
        for hostile in [
            json!({"yaw": 1}),
            json!({"yaw": "x", "pitch": 0}),
            json!({"yaw": 0, "pitch": 0, "distance": 0}),
            json!({"yaw": 0, "pitch": 0, "spin": 1}),
        ] {
            assert_eq!(call(&mut app, "ui.camera", hostile.clone())["ok"], false, "{hostile}");
        }
        // Given in the display units (mm), held in metres.
        assert!((app.cam.distance - 0.25).abs() < 1e-15);
        assert_eq!(call(&mut app, "ui.camera", json!({"yaw": -1e12, "pitch": 1e300, "distance": 1e300}))["ok"], true);
        // Aimed at a point it looks at it; at the model's centre of mass, framed (#467).
        let r = call(&mut app, "ui.camera", json!({"yaw": 30, "pitch": 20, "target": [1, 2, 3], "distance": 40}));
        assert!(app.cam.target.dist(fieldcraft_render::Vec3::new(0.001, 0.002, 0.003)) < 1e-15 && (app.cam.distance - 0.04).abs() < 1e-15, "{r}");
        assert_eq!(r["result"]["target"], json!({"x": 1.0, "y": 2.0, "z": 3.0}), "{r}");
        let r = call(&mut app, "ui.camera", json!({"yaw": 30, "pitch": 20, "target": {"center_of": "mass"}, "fit": true}));
        assert_eq!(r["ok"], true, "{r}");
        assert!(app.cam.target.dist(fieldcraft_render::Vec3::new(0.001, 0.002, 0.003)) > 1e-3 && app.cam.distance < 1e5, "{r}");
        assert_eq!(call(&mut app, "ui.camera", json!({"yaw": 0, "pitch": 0, "target": {"bodies": ["Nope"]}}))["ok"], false);
        assert!(app.cam.is_valid() && app.cam.pitch < std::f64::consts::FRAC_PI_2);
        assert_eq!(call(&mut app, "ui.inspect", json!({}))["result"]["model"]["elements"], 64 + 162);
        assert_eq!(call(&mut app, "ui.screenshot", json!({})), json!({"screenshot": true}));
        assert_eq!(call(&mut app, "no.method", json!({}))["ok"], false);
    }

    #[test]
    fn path_editor_control_validates_before_mutating_and_guards_pending_jobs() {
        let mut session = Session::default();
        session.run_script(&fieldcraft_engine::sample::cantilever()).unwrap();
        session.execute("solve.static", &json!({"step":"Load"})).unwrap();
        let mut app = ModelApp::new(session);
        for params in [json!([]), json!({"other":1}), json!({"picking":1}), json!({"clear":null}), json!({"at":[1]})] {
            assert_eq!(call(&mut app, "ui.path", params)["ok"], false);
            assert!(app.path_editor.is_none());
        }
        assert_eq!(call(&mut app, "ui.path", json!({"picking":true,"at":[1,2]}))["ok"], false);
        assert!(app.path_editor.is_none());
        assert_eq!(call(&mut app, "ui.path", json!({"picking":true}))["ok"], true);
        app.path_editor.as_mut().unwrap().points = vec![[0.02, 0.005, 0.02], [0.08, 0.005, 0.02]];
        app.sample_path().unwrap();
        assert_eq!(call(&mut app, "ui.path", json!({"clear":true,"picking":false}))["ok"], false);
        assert_eq!(app.path_editor.as_ref().unwrap().points.len(), 2);
        assert!(app.path_editor.as_ref().unwrap().picking);
        app.wait_loaded();
        assert_eq!(call(&mut app, "ui.path", json!({"clear":true,"picking":false}))["result"]["points"], 0);
        assert!(!app.path_editor.as_ref().unwrap().picking);
    }

    /// One frame of the whole window (1400 × 900 points), with the control channel's synthetic
    /// input, as the desktop app runs it.
    fn frame(app: &mut ModelApp, ctx: &egui::Context) {
        let mut raw =
            egui::RawInput { screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1400.0, 900.0))), ..Default::default() };
        app.raw_input_hook(&mut raw);
        let mut out = ctx.run_ui(raw, |ui| {
            app.logic(ui.ctx());
            app.ui(ui);
        });
        out.textures_delta.clear();
    }

    fn call_in(app: &mut ModelApp, ctx: &egui::Context, method: &str, params: Value) -> Value {
        let (tx, rx) = std::sync::mpsc::channel();
        let req = ControlRequest::new(method, params, tx);
        match handle(app, ctx, &req) {
            Outcome::Done(v) => v,
            Outcome::Screenshot { .. } => json!({"screenshot": true}),
            // A job's reply comes when it ends.
            Outcome::Pending => {
                app.wait_loaded();
                rx.try_recv().unwrap_or(Value::Null)
            }
        }
    }

    /// A meshed model's Materials node shows a highlighted ghost row; a click on it opens the
    /// material form, and `ui.inspect` reports each node's state.
    #[test]
    fn a_ghost_row_click_starts_its_command() {
        let mut app = ModelApp::new(Session::default());
        let ctx = egui::Context::default();
        let cube = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
        app.run("mesh.block", json!({"body": "B", "kind": "hex8", "corners": cube, "divisions": [1, 1, 1]})).unwrap();
        frame(&mut app, &ctx);
        let browser = call_in(&mut app, &ctx, "ui.inspect", json!({}))["result"]["browser"].clone();
        let materials = browser.as_array().and_then(|b| b.iter().find(|n| n["node"] == "materials")).cloned().unwrap_or_default();
        assert_eq!(materials["state"], "next", "{browser}");
        assert_eq!(materials["add"][0]["command"], "model.material.create");
        let at = app.layout[crate::browser::ghost_region(fieldcraft_engine::progress::Stage::Material)].center();
        call_in(&mut app, &ctx, "ui.click", json!({"x": at.x, "y": at.y}));
        for _ in 0..4 {
            frame(&mut app, &ctx);
        }
        assert_eq!(app.dialog.as_ref().map(|d| d.command), Some("model.material.create"));
    }

    /// The cantilever's free end (x = 100), looked at from +x: picked by a simulated click and
    /// by `ui.select`, made into a set through the set form's Target row, and resolved.
    #[test]
    fn the_cantilevers_end_is_picked_and_made_a_set_over_the_control_channel() {
        let mut app = ModelApp::new(Session::default());
        app.open_sample("cantilever");
        app.wait_loaded();
        let ctx = egui::Context::default();
        assert_eq!(call_in(&mut app, &ctx, "ui.view", json!({"view": "right"}))["ok"], true);
        app.fit_view();
        frame(&mut app, &ctx);
        let viewport = app.viewport.rect.unwrap();
        // A click on the end cell at y ∈ [0, 5], z ∈ [5, 10]: element 1 + 19 + 0·20 + 1·40, face 4.
        let at = crate::viewport::to_screen(&app, crate::results::tests::mm([100.0, 2.5, 7.5])).unwrap();
        assert_eq!(call_in(&mut app, &ctx, "ui.select", json!({"filter": "element_faces"}))["result"]["count"], 0);
        call_in(&mut app, &ctx, "ui.click", json!({"x": at.x, "y": at.y}));
        for _ in 0..4 {
            frame(&mut app, &ctx);
        }
        let r = call_in(&mut app, &ctx, "ui.selection", json!({}));
        assert_eq!(r["result"]["items"], json!([{"element": 60, "face": 4}]), "{r}");
        // What is under the pointer is the hover.
        assert_eq!(r["result"]["hover"], json!({"element": 60, "face": 4}));
        // The same pick by `ui.select {at}`, added to with Cmd (toggle) removes it again.
        let r = call_in(&mut app, &ctx, "ui.select", json!({"at": [at.x, at.y], "mode": "toggle"}));
        assert_eq!(r["result"]["count"], 0, "{r}");
        // A box over the whole viewport holds the end's eight element faces (the rest are hidden).
        let corners = [viewport.left(), viewport.top(), viewport.right(), viewport.bottom()];
        let r = call_in(&mut app, &ctx, "ui.select", json!({"box": corners}));
        assert_eq!(r["result"]["count"], 8, "{r}");
        assert!(r["result"]["items"].as_array().unwrap().iter().all(|f| f["face"] == 4));
        // The set form's Target row takes them; the set resolves to the end's faces and nodes.
        call_in(&mut app, &ctx, "ui.start", json!({"command": "model.set.create"}));
        let r = call_in(&mut app, &ctx, "ui.dialog", json!({"params": {"name": "end"}, "use_selection": true, "run": true}));
        assert_eq!(r["ok"], true, "{r}");
        let r = call_in(&mut app, &ctx, "engine.execute", json!({"command": "model.set.resolve", "params": {"name": "end"}}));
        assert_eq!((r["result"]["faces"]["count"].as_u64(), r["result"]["nodes"]["count"].as_u64()), (Some(8), Some(15)), "{r}");
        // A boundary condition's Target row names the set that has exactly the selection.
        call_in(&mut app, &ctx, "ui.start", json!({"command": "model.boundary_condition.create"}));
        let r = call_in(&mut app, &ctx, "ui.dialog", json!({"params": {"name": "Hold", "type": "fixed"}, "use_selection": true, "run": true}));
        assert_eq!(r["ok"], true, "{r}");
        let bc = serde_json::to_value(app.session.model.boundary_conditions().get("Hold").unwrap()).unwrap();
        assert_eq!(bc["target"], json!({"set": "end"}), "{bc}");
        // Nodes in the box: the end's 3 × 5, none of those behind it.
        let r = call_in(&mut app, &ctx, "ui.select", json!({"filter": "nodes", "box": corners}));
        assert_eq!(r["result"]["count"], 15, "{r}");
        // A selection no set has: Use selection names a set made only on OK, so Cancel leaves
        // the model as it was (#337).
        let r = call_in(&mut app, &ctx, "ui.select", json!({"nodes": [1, 2]}));
        assert_eq!(r["result"]["count"], 2, "{r}");
        let sets = app.session.model.sets().len();
        call_in(&mut app, &ctx, "ui.start", json!({"command": "model.boundary_condition.create"}));
        let r = call_in(&mut app, &ctx, "ui.dialog", json!({"params": {"name": "Pin", "type": "fixed"}, "use_selection": true}));
        assert_eq!(r["ok"], true, "{r}");
        assert_eq!(app.session.model.sets().len(), sets);
        app.close_dialog();
        assert_eq!(app.session.model.sets().len(), sets, "Cancel leaves no set behind");
        call_in(&mut app, &ctx, "ui.start", json!({"command": "model.boundary_condition.create"}));
        let r = call_in(&mut app, &ctx, "ui.dialog", json!({"params": {"name": "Pin", "type": "fixed"}, "use_selection": true, "run": true}));
        assert_eq!(r["ok"], true, "{r}");
        let bc = serde_json::to_value(app.session.model.boundary_conditions().get("Pin").unwrap()).unwrap();
        assert_eq!(bc["target"], json!({"set": "Pin-set"}), "{bc}");
        assert_eq!(app.session.model.resolve_set("Pin-set").unwrap().nodes, [1, 2]);
        assert_eq!(app.selected, Some((crate::browser::Def::BoundaryCondition, "Pin".into())), "what OK made is selected");
        // Escape clears the selection.
        call_in(&mut app, &ctx, "ui.key", json!({"key": "Escape"}));
        frame(&mut app, &ctx);
        assert_eq!(call_in(&mut app, &ctx, "ui.selection", json!({}))["result"]["count"], 0);
    }

    /// Text typed into one form never lands in the next one (#219): Young's modulus is focused
    /// and typed into, Plasticity opens straight over it, and its yield stress stays empty, both
    /// when the text arrives in the same frame as the switch and when it is typed afterwards.
    #[test]
    fn typing_never_crosses_into_the_next_form() {
        let ctx = egui::Context::default();
        for same_frame in [true, false] {
            let mut app = ModelApp::new(Session::default());
            // Plasticity needs a material to act on.
            app.run("model.material.create", json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 1, "poisson": 0.3}})).unwrap();
            call_in(&mut app, &ctx, "ui.start", json!({"command": "model.material.create"}));
            for _ in 0..2 {
                frame(&mut app, &ctx);
            }
            let opening = app.dialog.as_ref().unwrap().opening;
            let form = egui::Id::new(("fc_form", "model.material.create", opening));
            let young = ctx.read_response(egui::Id::new(("fc_form", form, "elastic/young"))).expect("the Young's modulus field").rect.center();
            call_in(&mut app, &ctx, "ui.click", json!({"x": young.x, "y": young.y}));
            for _ in 0..3 {
                frame(&mut app, &ctx);
            }
            assert!(ctx.memory(|m| m.focused()).is_some(), "the click focuses Young's modulus");
            if !same_frame {
                call_in(&mut app, &ctx, "ui.start", json!({"command": "model.material.plasticity"}));
            }
            call_in(&mut app, &ctx, "ui.text", json!({"text": "210000"}));
            if same_frame {
                call_in(&mut app, &ctx, "ui.start", json!({"command": "model.material.plasticity"}));
            }
            for _ in 0..3 {
                frame(&mut app, &ctx);
            }
            let d = app.dialog.as_ref().unwrap();
            assert_eq!(d.command, "model.material.plasticity");
            assert!(!d.text.contains("210000"), "same frame {same_frame}: {}", d.text);
            assert!(ctx.memory(|m| m.focused()).is_none(), "the new form does not inherit the focus");
        }
    }

    #[test]
    fn hostile_selections_are_refused_or_capped() {
        let mut app = ModelApp::new(Session::default());
        let ctx = egui::Context::default();
        // Before the viewport is drawn, a point has nothing to pick on.
        assert_eq!(call_in(&mut app, &ctx, "ui.select", json!({"at": [10, 10]}))["ok"], false);
        // An empty scene picks nothing anywhere.
        frame(&mut app, &ctx);
        for filter in ["faces", "edges", "vertices", "bodies", "nodes", "element_faces", "elements"] {
            let r = call_in(&mut app, &ctx, "ui.select", json!({"filter": filter, "at": [700, 450]}));
            assert_eq!(r["result"]["count"], 0, "{filter}: {r}");
        }
        for bad in [
            json!({"filter": "atoms"}),
            json!({"mode": "xor", "nodes": [1]}),
            json!({"at": [1e300, 5]}),
            json!({"at": "here"}),
            json!({"at": [1, 2, 3]}),
            json!({"box": [0, 0, 1]}),
            json!({"box": [0, f64::MAX, 1, 1]}),
            json!({"filter": "faces", "box": [0, 0, 100, 100]}),
            json!({"nodes": [0]}),
            json!({"faces": [{"element": 1}]}),
            json!({"geometry": [{"body": "B", "face": ""}]}),
            json!({"elements": vec![1; crate::selection::MAX_SELECTION + 1]}),
        ] {
            let r = call_in(&mut app, &ctx, "ui.select", bad.clone());
            assert_eq!(r["ok"], false, "{bad}: {r}");
        }
        // A huge selection is capped and says so.
        let all: Vec<u32> = (1..=crate::selection::MAX_SELECTION as u32).collect();
        call_in(&mut app, &ctx, "ui.select", json!({"nodes": all}));
        let r = call_in(&mut app, &ctx, "ui.select", json!({"elements": [1, 2], "mode": "add"}));
        assert_eq!((r["result"]["count"].as_u64(), r["result"]["capped"].as_bool()), (Some(crate::selection::MAX_SELECTION as u64), Some(true)));
        // The form's Target row: no form, a command without a target.
        assert_eq!(call_in(&mut app, &ctx, "ui.dialog", json!({"use_selection": true}))["ok"], false);
        call_in(&mut app, &ctx, "ui.start", json!({"command": "mesh.node.add"}));
        assert_eq!(call_in(&mut app, &ctx, "ui.dialog", json!({"use_selection": true}))["ok"], false);
        // `ui.select` with a definition still selects in the browser; `clear` empties the viewport's.
        assert_eq!(call_in(&mut app, &ctx, "ui.select", json!({"clear": true}))["result"]["count"], 0);
        frame(&mut app, &ctx);
    }

    #[test]
    fn the_command_form_is_driven_after_ui_start() {
        let mut app = ModelApp::new(Session::default());
        assert_eq!(call(&mut app, "ui.dialog", json!({"run": true}))["ok"], false);
        assert_eq!(call(&mut app, "ui.start", json!({"command": "mesh.node.add"}))["result"]["dialog"], "mesh.node.add");
        let r = call(&mut app, "ui.dialog", json!({"params": {"label": 3, "xyz": [1, 2, 3]}, "run": true}));
        assert_eq!(r["ok"], true, "{r}");
        assert!(app.dialog.is_none() && app.session.model.node(3).is_some());
    }

    /// A sample opens as a job however it is asked for (the Samples window, `engine.execute`):
    /// the window and the control channel answer while it is built and meshed (#552).
    #[test]
    fn a_sample_opens_while_the_window_answers() {
        let mut app = ModelApp::new(Session::default());
        assert_eq!(call(&mut app, "ui.start", json!({"command": "file.sample"}))["ok"], true);
        let r = call(&mut app, "ui.prompt", json!({"answer": "open", "name": "cantilever"}));
        assert_eq!(r["ok"], true, "{r}");
        let r = call(&mut app, "ui.inspect", json!({}));
        assert_eq!((r["result"]["busy"].as_str(), r["result"]["model"]["nodes"].as_u64()), (Some("Opening the Cantilever sample"), Some(0)), "{r}");
        // Commands and scripts wait for it: the opened session replaces the window's.
        let script = json!({"commands": [{"command": "mesh.node.add", "params": {"label": 1, "xyz": [0, 0, 0]}}]});
        let r = call(&mut app, "engine.script", script);
        assert!(r["error"].as_str().is_some_and(|e| e.contains("Opening the Cantilever sample")), "{r}");
        app.wait_loaded();
        assert!(!app.session.model.nodes().is_empty());
        let r = call(&mut app, "engine.execute", json!({"command": "file.sample", "params": {"name": "impact"}, "wait": false}));
        assert_eq!((r["ok"].as_bool(), app.busy()), (Some(true), Some("Opening the Bar Impact sample")), "{r}");
        let r = call(&mut app, "job.wait", json!({}));
        assert_eq!(r["result"]["sample"], "impact", "{r}");
        assert!(app.session.model.step("Impact").is_some() && !app.session.is_dirty());
    }

    /// `ui.drop` stands in for the window system: the next frame opens the files as a drop does.
    #[test]
    fn dropped_files_open_as_a_drop_does() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-drop-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let design = dir.join("plate.solidcraft");
        std::fs::write(&design, fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap()).unwrap();
        let mut app = ModelApp::new(Session::default());
        let r = call(&mut app, "ui.drop", json!({"paths": [design.to_string_lossy()]}));
        assert_eq!(r["result"]["dropped"], 1, "{r}");
        assert!(app.session.design.is_none(), "nothing opens before the next frame");
        app.logic(&egui::Context::default());
        assert_eq!(app.busy(), Some("Opening plate.solidcraft"));
        app.wait_loaded();
        assert_eq!(app.session.design.as_ref().map(|d| d.bodies.len()), Some(1));
        assert_eq!(app.status.as_ref().map(|s| s.text.as_str()), Some("Opened plate.solidcraft"));
        for bad in [json!({}), json!({"paths": []}), json!({"paths": [""]}), json!({"paths": [3]}), json!({"paths": vec!["a.solidcraft"; 17]})] {
            assert_eq!(call(&mut app, "ui.drop", bad.clone())["ok"], false, "{bad}");
        }
        assert!(app.synthetic_drops.is_empty());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn ui_render_writes_a_png_of_the_model() {
        let mut app = ModelApp::new(Session::default());
        app.run_script(&fieldcraft_engine::sample::script()).unwrap();
        let path = std::env::temp_dir().join(format!("fieldcraft-render-{}.png", std::process::id()));
        let r = call(&mut app, "ui.render", json!({"path": path.to_string_lossy(), "width": 120, "height": 80}));
        assert_eq!(r["ok"], true, "{r}");
        assert_eq!(r["result"]["renderer"]["kind"], "cpu", "{r}");
        let bytes = std::fs::read(&path).unwrap();
        assert!(bytes.starts_with(&[0x89, b'P', b'N', b'G']));
        let _ = std::fs::remove_file(&path);
        assert_eq!(call(&mut app, "ui.render", json!({}))["ok"], false);
    }

    #[test]
    fn ui_render_refuses_hostile_sizes() {
        let mut app = ModelApp::new(Session::default());
        app.run_script(&fieldcraft_engine::sample::script()).unwrap();
        let path = std::env::temp_dir().join(format!("fieldcraft-render-hostile-{}.png", std::process::id()));
        let path = path.to_string_lossy();
        for (w, h) in [
            (json!(100_000), json!(100_000)),
            (json!(8192), json!(8192)),
            (json!(1e12), json!(10)),
            (json!(-5), json!(10)),
            (json!(0), json!(10)),
            (json!(10.5), json!(10)),
            (json!("big"), json!(10)),
            (json!(f64::MAX), json!(f64::MIN)),
        ] {
            let r = call(&mut app, "ui.render", json!({"path": path, "width": w, "height": h}));
            assert_eq!(r["ok"], false, "{w} x {h}: {r}");
        }
        assert!(std::fs::metadata(&*path).is_err(), "nothing is written for a refused size");
    }

    #[test]
    fn hostile_requests_never_panic() {
        let mut app = ModelApp::new(Session::default());
        let hostile = [
            json!(null),
            json!([]),
            json!({"x": 1e300, "y": f64::MAX, "x0": -1e308, "y0": 0, "x1": 1e308, "y1": 1, "steps": 1e12, "delta": 1e300}),
            json!({"width": 1e12, "height": -5, "path": ""}),
            json!({"command": 5, "params": "x", "commands": {"a": 1}}),
            json!({"key": "NoSuchKey", "text": 7, "view": null, "points": [[1], "x", null]}),
            json!({"tab": "nope", "command": "mesh.tet.generate"}),
            json!({"tab": 3}),
        ];
        for method in [
            "engine.execute",
            "engine.script",
            "ui.set",
            "ui.view",
            "ui.camera",
            "ui.start",
            "ui.dialog",
            "ui.move",
            "ui.click",
            "ui.drag",
            "ui.scroll",
            "ui.key",
            "ui.text",
            "ui.worldToScreen",
            "ui.resize",
            "ui.render",
            "ui.capture",
            "ui.toolbar",
            "ui.prompt",
            "ui.item",
            "ui.select",
            "ui.path",
            "ui.selection",
            "ui.fullscreen",
            "ui.drop",
        ] {
            for p in &hostile {
                let _ = call(&mut app, method, p.clone());
            }
        }
        assert!(app.synthetic.len() < 1000);
    }

    #[test]
    fn the_toolbar_lists_every_entry_and_switches_tabs() {
        let mut app = ModelApp::new(Session::default());
        assert_eq!(call(&mut app, "ui.set", json!({"tab": "mesh"}))["ok"], true);
        let r = call(&mut app, "ui.toolbar", json!({}));
        assert_eq!(r["result"]["tab"], "mesh", "{r}");
        let entries = r["result"]["entries"].as_array().unwrap();
        assert_eq!(entries.len(), fieldcraft_engine::catalogue::CATALOGUE.len());
        let soon = entries.iter().find(|e| e["id"] == "mesh.order.quadratic").unwrap();
        assert_eq!((soon["coming_soon"].as_str(), soon["command"].is_null()), (Some("M2"), true));
        assert_eq!(call(&mut app, "ui.set", json!({"tab": "nope"}))["ok"], false);
        // A workspace brings its first tab; a tab brings its workspace.
        let r = call(&mut app, "ui.set", json!({"workspace": "model"}));
        assert_eq!((r["result"]["workspace"].as_str(), r["result"]["tab"].as_str()), (Some("model"), Some("properties")));
        let r = call(&mut app, "ui.set", json!({"tab": "solve"}));
        assert_eq!(r["result"]["workspace"], "analysis");
        let r = call(&mut app, "ui.start", json!({"command": "mesh.order.quadratic"}));
        assert_eq!((r["ok"].as_bool(), r["error"].as_str()), (Some(false), Some("Quadratic is coming soon (M2)")));
    }

    #[test]
    fn ui_capture_saves_the_whole_window() {
        let mut app = ModelApp::new(Session::default());
        let path = std::env::temp_dir().join(format!("fieldcraft-capture-{}.png", std::process::id()));
        let r = call(&mut app, "ui.capture", json!({"path": path.to_string_lossy(), "width": 640, "height": 400, "scale": 1}));
        assert_eq!(
            (r["ok"].as_bool(), r["result"]["width"].as_u64(), r["result"]["offscreen"].as_bool()),
            (Some(true), Some(640), Some(true)),
            "{r}"
        );
        assert!(std::fs::read(&path).unwrap().starts_with(&[0x89, b'P', b'N', b'G']));
        let _ = std::fs::remove_file(&path);
        for bad in [json!({"width": "x"}), json!({"width": 1e9}), json!({"scale": 0}), json!({"height": -1})] {
            assert_eq!(call(&mut app, "ui.capture", bad.clone())["ok"], false, "{bad}");
        }
    }

    #[test]
    fn screenshots_are_saved_as_png() {
        let image = egui::ColorImage::filled([4, 3], egui::Color32::RED);
        let path = std::env::temp_dir().join(format!("fieldcraft-shot-{}.png", std::process::id()));
        let r = save_screenshot(&image, Some(&path.to_string_lossy()));
        assert_eq!(r["ok"], true, "{r}");
        assert!(std::fs::read(&path).unwrap().starts_with(&[0x89, b'P', b'N', b'G']));
        let _ = std::fs::remove_file(&path);
    }
}
