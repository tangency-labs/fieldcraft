//! Help › Report a Bug, and the question the window asks after a crash ("FieldCraft closed
//! unexpectedly — Report it").
//!
//! The form shows exactly what the GitHub issue form will be filled with
//! ([`fieldcraft_engine::bug_report`]): the version, the OS, what happened (the last error by
//! default) and, after a crash, where it happened. **Share my session** adds the commands that
//! built the model, **Include the model file** the model as well. Nothing leaves the computer
//! until the user clicks **Open GitHub Issue**, which opens the filled form in their browser (they
//! submit it there, or not); files too big for the link are saved to a folder to attach from, and
//! a long session is also put on the clipboard.

use egui::RichText;
use fieldcraft_engine::bug_report::{self, Draft, Input};
use fieldcraft_engine::crash::{self, Summary};
use serde_json::{Value, json};

use crate::ModelApp;
use crate::files::Prompt;

/// The Report a Bug form.
#[derive(Debug, Default)]
pub struct ReportForm {
    pub title: String,
    pub observed: String,
    /// The crash report being reported (its id).
    pub crash: Option<String>,
    pub share_session: bool,
    pub include_model: bool,
    /// What the draft is built from (the session's commands and model, or the crash report's),
    /// read when the form opens and when the model is first included.
    input: Input,
    model_read: bool,
    /// The draft, and the form state it was made for.
    draft: Option<Draft>,
    drafted: Option<Key>,
    /// Why the input could not be read.
    pub problem: Option<String>,
}

type Key = (String, String, bool, bool);

impl ReportForm {
    /// The form and its draft as last made (for `ui.inspect`).
    pub fn drafted_json(&self) -> Option<Value> {
        let mut v = self.draft.as_ref()?.to_json();
        v["share_session"] = json!(self.share_session);
        v["include_model"] = json!(self.include_model);
        v["crash"] = json!(self.crash);
        Some(v)
    }

    /// The draft for the form as it is now (made again only when something changed).
    pub fn draft(&mut self) -> &Draft {
        let key = (self.title.clone(), self.observed.clone(), self.share_session, self.include_model);
        if self.drafted.as_ref() != Some(&key) {
            self.draft = None;
            self.drafted = Some(key);
        }
        self.draft.get_or_insert_with(|| {
            let mut input = self.input.clone();
            input.title = self.title.clone();
            input.observed = self.observed.clone();
            input.share_session = self.share_session;
            input.include_model = self.include_model;
            bug_report::draft(&input, &crash::Build::current().version_line())
        })
    }
}

/// The input for `form` from the session (or its crash report).
fn read_input(app: &ModelApp, form: &ReportForm) -> Result<Input, String> {
    let mut p = json!({"share_session": true, "include_model": form.include_model});
    if let Some(id) = &form.crash {
        p["crash"] = json!(id);
    }
    fieldcraft_engine::cmd::diagnostics::bug_report_input(&app.session, &p).map_err(|e| e.to_string())
}

impl ModelApp {
    /// Open Help › Report a Bug, for `crash` (a report's id) or what happened last.
    pub fn report_bug(&mut self, crash: Option<String>) {
        let mut form = ReportForm {
            observed: if crash.is_some() { String::new() } else { self.last_error.clone().unwrap_or_default() },
            crash,
            ..Default::default()
        };
        match read_input(self, &form) {
            Ok(input) => form.input = input,
            Err(e) => form.problem = Some(e),
        }
        self.report = Some(form);
        self.prompt = Some(Prompt::ReportBug);
        self.ui.menu = None;
    }

    /// Help › Crash Reports: the saved reports (`diagnostics.crash_reports`), newest first.
    pub fn open_crash_reports(&mut self) {
        // Read directly: a listing changes nothing, and the window shows it instead of the reply.
        #[allow(clippy::disallowed_methods, reason = "lists the crash report folder: a few small files")]
        let reports = match self.session.execute("diagnostics.crash_reports", &json!({})) {
            Ok(r) => r.get("reports").and_then(Value::as_array).map_or_else(Vec::new, |a| a.iter().filter_map(summary_of).collect()),
            Err(e) => {
                self.set_status(format!("Crash reports: {e}"), true);
                return;
            }
        };
        self.prompt = Some(Prompt::CrashReports(reports));
        self.ui.menu = None;
    }

    /// Delete one crash report, or all (`diagnostics.crash_reports.delete`), and list them again.
    pub fn delete_crash_reports(&mut self, id: Option<&str>) -> Result<Value, String> {
        let params = match id {
            Some(id) => json!({"id": id}),
            None => json!({"all": true}),
        };
        let r = self.run("diagnostics.crash_reports.delete", params)?;
        let n = r.get("deleted").and_then(Value::as_u64).unwrap_or(0);
        self.set_status(format!("Deleted {n} crash report{}", if n == 1 { "" } else { "s" }), false);
        self.open_crash_reports();
        Ok(r)
    }

    /// The host found these crash reports at startup: ask about the real crashes.
    pub fn announce_crashes(&mut self, fresh: &[Summary]) {
        let crashes: Vec<Summary> = fresh.iter().filter(|r| !matches!(r.kind.as_str(), "recovered_panic" | "solver_error")).cloned().collect();
        if crashes.is_empty() {
            return;
        }
        log::warn!(
            "{} crash report(s) since the last start: {}",
            crashes.len(),
            crashes.iter().map(|c| c.id.as_str()).collect::<Vec<_>>().join(", ")
        );
        if self.prompt.is_none() {
            self.prompt = Some(Prompt::Crashed(crashes));
        }
    }

    /// Change the open form (`ui.prompt` `set`): any of title, observed, share_session,
    /// include_model.
    pub fn set_report(&mut self, p: &Value) -> Result<Value, String> {
        let mut form = self.report.take().ok_or("Report a Bug is not open")?;
        if let Some(t) = p.get("title").and_then(Value::as_str) {
            form.title = t.chars().take(500).collect();
        }
        if let Some(t) = p.get("observed").and_then(Value::as_str) {
            form.observed = t.chars().take(20_000).collect();
        }
        if let Some(b) = p.get("share_session").and_then(Value::as_bool) {
            form.share_session = b;
        }
        if let Some(b) = p.get("include_model").and_then(Value::as_bool) {
            form.include_model = b;
        }
        self.include_model_if_asked(&mut form);
        let draft = form.draft().to_json();
        self.report = Some(form);
        Ok(draft)
    }

    /// Read the model once it is first included.
    fn include_model_if_asked(&mut self, form: &mut ReportForm) {
        if form.include_model && !form.model_read {
            form.model_read = true;
            match read_input(self, form) {
                Ok(input) => form.input = input,
                Err(e) => form.problem = Some(e),
            }
        }
    }

    /// Open the draft in the browser: attachments saved to a folder of their own, the session on
    /// the clipboard when it is too long for the link.
    fn send_report(&mut self, ctx: &egui::Context, draft: &Draft) {
        let mut saved = None;
        if !draft.attachments.is_empty() {
            match bug_report::save_attachments(draft) {
                Ok(dir) => saved = Some(dir),
                Err(e) => {
                    self.set_status(format!("Report a Bug: could not save the files to attach: {e}"), true);
                    return;
                }
            }
        }
        if let Some(text) = &draft.clipboard {
            ctx.copy_text(text.clone());
        }
        match &self.services.open_url {
            Some(open) => {
                if let Err(e) = open(&draft.url) {
                    self.set_status(format!("Report a Bug: could not open the browser: {e}"), true);
                    return;
                }
            }
            None => ctx.open_url(egui::OpenUrl::new_tab(draft.url.clone())),
        }
        log::info!("bug report opened in the browser ({} bytes, {} attachment(s))", draft.url.len(), draft.attachments.len());
        let text = match saved {
            Some(dir) => format!("The issue form opened in your browser: attach the files in {} and submit it there.", dir.display()),
            None => "The issue form opened in your browser: review it and submit it there.".to_string(),
        };
        self.set_status(text, false);
        self.prompt = None;
        self.report = None;
    }
}

/// A crash report as `diagnostics.crash_reports` lists it.
fn summary_of(v: &Value) -> Option<Summary> {
    let s = |k: &str| v.get(k).and_then(Value::as_str).unwrap_or_default().to_string();
    Some(Summary {
        id: v.get("id").and_then(Value::as_str)?.to_string(),
        kind: s("kind"),
        time: s("time"),
        message: s("message"),
        command: v.get("command").and_then(Value::as_str).map(str::to_string),
        count: v.get("count").and_then(Value::as_u64).unwrap_or(1),
        version: s("version"),
        path: s("path"),
        replayable: v.get("replayable").and_then(Value::as_bool).unwrap_or(false),
    })
}

fn size(bytes: usize) -> String {
    if bytes < 1024 { format!("{bytes} B") } else { format!("{} KB", bytes.div_ceil(1024)) }
}

fn kind_label(kind: &str) -> &'static str {
    match kind {
        "unclean_exit" => "FieldCraft closed unexpectedly",
        "command_panic" => "A command hit an internal error",
        "solver_error" => "A solve stopped",
        _ => "FieldCraft hit an internal error",
    }
}

/// The crash question and the Report a Bug form, when open.
pub fn show(app: &mut ModelApp, ctx: &egui::Context) {
    let t = crate::theme::Tokens::get(app.ui.dark);
    match app.prompt.clone() {
        Some(Prompt::Crashed(reports)) => {
            let Some(first) = reports.first() else {
                app.prompt = None;
                return;
            };
            let (mut report, mut close) = (false, false);
            crate::frame::window(ctx, kind_label(&first.kind), crate::frame::Width::Normal).id(egui::Id::new("fc_crashed")).show(ctx, |ui| {
                ui.label(match first.kind.as_str() {
                    "unclean_exit" => "FieldCraft closed unexpectedly the last time it ran.",
                    _ => "FieldCraft hit an internal error the last time it ran.",
                });
                ui.label(RichText::new(&first.message).color(t.text_dim));
                if reports.len() > 1 {
                    ui.label(RichText::new(format!("And {} more since.", reports.len() - 1)).color(t.text_dim));
                }
                ui.add_space(4.0);
                ui.label(
                    RichText::new(
                        "A crash report was saved on this computer. Reporting it opens a GitHub issue you review first; nothing is sent until you submit it.",
                    )
                    .color(t.text_dim),
                );
                ui.add_space(6.0);
                ui.horizontal(|ui| {
                    report = ui.button("Report It…").clicked();
                    close = ui.button("Not Now").clicked();
                });
            });
            if report {
                app.report_bug(Some(first.id.clone()));
            } else if close || ctx.input(|i| i.key_pressed(egui::Key::Escape)) {
                app.prompt = None;
            }
        }
        Some(Prompt::CrashReports(reports)) => {
            enum Act {
                Report(String),
                Copy(String),
                Delete(String),
                DeleteAll,
                Close,
            }
            let mut act = None;
            crate::frame::window(ctx, "Crash Reports", crate::frame::Width::Wide).id(egui::Id::new("fc_crash_reports")).show(ctx, |ui| {
                if reports.is_empty() {
                    ui.label("No crash reports are saved on this computer.");
                } else {
                    ui.label(
                        RichText::new("Saved on this computer, newest first. Reporting one opens a GitHub issue you review first; nothing is sent until you submit it.")
                            .color(t.text_dim),
                    );
                    ui.add_space(4.0);
                    egui::ScrollArea::vertical().max_height((ctx.content_rect().height() * 0.6).max(160.0)).show(ui, |ui| {
                        for r in &reports {
                            ui.separator();
                            ui.horizontal(|ui| {
                                ui.label(RichText::new(kind_label(&r.kind)).strong());
                                let when = if r.count > 1 { format!("{} · {} times", r.time, r.count) } else { r.time.clone() };
                                ui.label(RichText::new(when).color(t.text_dim));
                            });
                            ui.label(RichText::new(&r.message).color(t.text_dim).size(11.0));
                            ui.horizontal(|ui| {
                                if ui.button("Report…").clicked() {
                                    act = Some(Act::Report(r.id.clone()));
                                }
                                if ui.button("Copy Path").on_hover_text(&r.path).clicked() {
                                    act = Some(Act::Copy(r.path.clone()));
                                }
                                if ui.button("Delete").clicked() {
                                    act = Some(Act::Delete(r.id.clone()));
                                }
                            });
                        }
                    });
                }
                ui.add_space(6.0);
                ui.horizontal(|ui| {
                    if ui.add_enabled(!reports.is_empty(), egui::Button::new("Delete All")).clicked() {
                        act = Some(Act::DeleteAll);
                    }
                    if ui.button("Close").clicked() {
                        act = Some(Act::Close);
                    }
                });
            });
            match act {
                Some(Act::Report(id)) => app.report_bug(Some(id)),
                Some(Act::Copy(path)) => {
                    ctx.copy_text(path);
                    app.set_status("Copied the crash report's folder", false);
                }
                Some(Act::Delete(id)) => {
                    let _ = app.delete_crash_reports(Some(&id));
                }
                Some(Act::DeleteAll) => {
                    let _ = app.delete_crash_reports(None);
                }
                Some(Act::Close) => app.prompt = None,
                None if ctx.input(|i| i.key_pressed(egui::Key::Escape)) => app.prompt = None,
                None => {}
            }
        }
        Some(Prompt::ReportBug) => {
            let Some(mut form) = app.report.take() else {
                app.prompt = None;
                return;
            };
            let (mut send, mut close) = (false, false);
            let commands = form.input.commands.len();
            crate::frame::window(ctx, "Report a Bug", crate::frame::Width::Wide).id(egui::Id::new("fc_report_bug")).show(ctx, |ui| {
                ui.label(
                    RichText::new(
                        "This opens a GitHub issue in your browser, filled in with what is shown below. Nothing is sent until you submit it there.",
                    )
                    .color(t.text_dim),
                );
                if let Some(p) = &form.problem {
                    ui.label(RichText::new(p).color(t.error));
                }
                ui.add_space(4.0);
                ui.label("Title");
                let hint = form.draft().title.clone();
                ui.add(egui::TextEdit::singleline(&mut form.title).hint_text(hint).desired_width(f32::INFINITY));
                ui.label("What happened");
                ui.add(
                    egui::TextEdit::multiline(&mut form.observed)
                        .hint_text("What you did, and what went wrong")
                        .desired_rows(3)
                        .desired_width(f32::INFINITY),
                );
                ui.add_space(4.0);
                ui.checkbox(&mut form.share_session, "Share my session")
                    .on_hover_text("Attach the commands that built the model, to replay it. They can name your files.");
                ui.label(
                    RichText::new(format!("The {commands} command(s) that built the model, the version and the OS; they can name your files."))
                        .color(t.text_dim)
                        .size(11.0),
                );
                ui.add_enabled_ui(form.share_session, |ui| {
                    ui.checkbox(&mut form.include_model, "Include the model file")
                        .on_hover_text("Attach the model too (you attach the file yourself).");
                });
                if !form.share_session {
                    form.include_model = false;
                }
                app.include_model_if_asked(&mut form);
                let draft = form.draft().clone();
                ui.add_space(6.0);
                ui.horizontal(|ui| {
                    send = ui.button("Open GitHub Issue").clicked();
                    close = ui.button("Cancel").clicked();
                });
                ui.add_space(4.0);
                egui::CollapsingHeader::new("What is sent").id_salt("fc_report_sent").default_open(true).show(ui, |ui| {
                    egui::ScrollArea::vertical().max_height(180.0).show(ui, |ui| {
                        ui.label(RichText::new(format!("Title: {}", draft.title)).monospace().size(11.0));
                        for (k, v) in &draft.fields {
                            if !v.is_empty() {
                                ui.label(RichText::new(format!("{k}:")).strong().size(11.0));
                                ui.label(RichText::new(v).monospace().size(11.0));
                            }
                        }
                        for a in &draft.attachments {
                            ui.label(RichText::new(format!("Attach {} ({}): {}", a.name, size(a.bytes.len()), a.what)).size(11.0).color(t.accent));
                        }
                        if draft.clipboard.is_some() {
                            ui.label(RichText::new("The session goes on the clipboard too, to paste.").size(11.0).color(t.text_dim));
                        }
                    });
                });
            });
            let draft = form.draft().clone();
            app.report = Some(form);
            if send {
                app.send_report(ctx, &draft);
            } else if close || ctx.input(|i| i.key_pressed(egui::Key::Escape)) {
                app.prompt = None;
                app.report = None;
            }
        }
        _ => {}
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_engine::Session;

    #[test]
    fn report_a_bug_shows_the_last_error_and_shares_only_when_ticked() {
        let mut app = ModelApp::new(Session::default());
        app.run("mesh.node.add", json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
        let _ = app.run("mesh.node.remove", json!({"label": 99}));
        let error = app.last_error.clone().expect("the failure is the last error");
        app.report_bug(None);
        assert_eq!(app.prompt, Some(Prompt::ReportBug));
        let form = app.report.as_mut().unwrap();
        assert_eq!(form.observed, error);
        let d = form.draft().clone();
        assert!(d.attachments.is_empty() && !bug_report_text(&d).contains("mesh.node.add"), "no session unless shared");
        let d = app.set_report(&json!({"share_session": true})).unwrap();
        assert!(d["url"].as_str().unwrap().contains("mesh.node.add"), "{d}");
        let d = app.set_report(&json!({"include_model": true})).unwrap();
        assert_eq!(d["attachments"][0]["name"], "model.fieldcraft", "{d}");
        // The model goes only with the session.
        let d = app.set_report(&json!({"share_session": false})).unwrap();
        assert_eq!(d["attachments"], json!([]), "{d}");
    }

    /// Help › Crash Reports reads each report as `diagnostics.crash_reports` lists it.
    #[test]
    fn crash_reports_read_the_listing() {
        let r = Summary {
            id: "20261014-101500-0000000000000001".into(),
            kind: "unclean_exit".into(),
            time: "2026-10-14T10:15:00Z".into(),
            message: "ended".into(),
            command: Some("solve.run".into()),
            count: 2,
            version: "0.5.9".into(),
            path: "/x".into(),
            replayable: true,
        };
        assert_eq!(summary_of(&serde_json::to_value(&r).unwrap()), Some(r));
        assert_eq!(summary_of(&json!({"kind": "panic"})), None, "a report without an id is left out");
    }

    fn bug_report_text(d: &Draft) -> String {
        d.fields.iter().map(|(_, v)| v.as_str()).collect::<Vec<_>>().join("\n")
    }

    #[test]
    fn a_crash_at_startup_asks_and_report_it_opens_the_form() {
        let mut app = ModelApp::new(Session::default());
        let quiet = Summary {
            id: "20261014-101500-0000000000000001".into(),
            kind: "recovered_panic".into(),
            time: String::new(),
            message: String::new(),
            command: None,
            count: 1,
            version: String::new(),
            path: String::new(),
            replayable: false,
        };
        app.announce_crashes(std::slice::from_ref(&quiet));
        assert_eq!(app.prompt, None, "recovered panics do not interrupt");
        let solver = Summary { kind: "solver_error".into(), ..quiet.clone() };
        app.announce_crashes(&[solver]);
        assert_eq!(app.prompt, None, "ordinary solver records do not announce an application crash");
        let crashed = Summary { kind: "unclean_exit".into(), message: "FieldCraft ended without a clean exit".into(), ..quiet };
        app.announce_crashes(std::slice::from_ref(&crashed));
        assert!(matches!(&app.prompt, Some(Prompt::Crashed(r)) if r.len() == 1));
        app.report_bug(Some(crashed.id.clone()));
        let form = app.report.as_ref().unwrap();
        assert_eq!(form.crash.as_deref(), Some(crashed.id.as_str()));
        // The report is not in this test's crash folder: the form says so and still drafts.
        assert!(form.problem.as_deref().is_some_and(|p| p.contains("no crash report")), "{:?}", form.problem);
    }
}
