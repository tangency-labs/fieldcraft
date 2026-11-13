//! Anonymous usage sharing (docs/privacy.md): the question the first run asks once, and Help ›
//! Settings, where it is changed. The choice is kept in the preferences ([`crate::files::Prefs`]);
//! the events are [`fieldcraft_engine::usage`]'s, and the host sends them.
//!
//! The window records only when the host says this run may ([`ModelApp::with_usage`]: the build
//! has an analytics key, and neither `FIELDCRAFT_TELEMETRY=0`, CI nor automation turned it off)
//! and the user said yes.

use fieldcraft_engine::usage::{self, Event, Identity, Origin};

use crate::ModelApp;
use crate::files::{Prompt, dismissed, message, note};
use solidcraft_ui_kit as kit;
use solidcraft_ui_kit::BackingTooltipExt as _;

/// What is sent, event by event.
pub const PRIVACY_URL: &str = "https://github.com/tangency-labs/fieldcraft/blob/main/docs/privacy.md";
/// The question, as the first run asks it.
pub const QUESTION: &str = "Help improve FieldCraft by sharing anonymous usage: which tools you use, nothing from your models.";
const DETAIL: &str = "Counted: the commands you run (their names, how long they took, whether they failed), the steps from a design to results, \
                      the samples you open, and the version, OS and a random id for this install. Never sent: file names or paths, geometry, model \
                      contents, parameter values or messages.";

/// What the host provides when this run may share usage.
pub struct UsageHost {
    /// `desktop` or `web`.
    pub platform: &'static str,
    pub os: &'static str,
    pub arch: &'static str,
    /// A new random install id; `None` when there is no randomness to make one.
    pub new_id: Box<dyn Fn() -> Option<String>>,
}

impl ModelApp {
    /// This run may share usage (the host checked the build's key and the environment): ask the
    /// user once, or record if they already agreed.
    pub fn with_usage(mut self, host: UsageHost) -> Self {
        self.usage = Some(host);
        self.apply_usage();
        self
    }

    /// The user's preferences as last saved.
    pub fn prefs(&self) -> &crate::files::Prefs {
        &self.prefs
    }

    /// Whether this run can share usage at all.
    pub fn usage_available(&self) -> bool {
        self.usage.is_some()
    }

    /// Start or stop recording as the preferences say.
    fn apply_usage(&mut self) {
        let Some(host) = &self.usage else { return };
        if self.prefs.usage_sharing != Some(true) {
            usage::disable();
            return;
        }
        if self.prefs.install_id.is_empty() {
            let Some(id) = (host.new_id)() else {
                log::warn!("usage sharing: no random install id could be made");
                return;
            };
            self.prefs.install_id = id;
        }
        let identity = Identity {
            install_id: self.prefs.install_id.clone(),
            platform: host.platform,
            os: host.os,
            arch: host.arch,
            app_version: crate::VERSION,
        };
        usage::enable(identity, self.prefs.usage_automation);
    }

    /// Change the choices (the consent question, Help › Settings, `ui.prompt`) and keep them.
    pub fn set_usage(&mut self, sharing: bool, automation: bool) {
        self.prefs.usage_sharing = Some(sharing);
        self.prefs.usage_automation = automation;
        // Agreeing makes the install id, which is kept with the choice.
        self.apply_usage();
        self.prefs.save(self.services.prefs.as_ref());
    }

    /// Answer the first run's question.
    pub fn answer_consent(&mut self, share: bool) {
        if self.prompt == Some(Prompt::Consent) {
            self.prompt = None;
        }
        self.set_usage(share, self.prefs.usage_automation);
    }

    /// Open Help › Settings.
    pub fn open_settings(&mut self) {
        self.prompt = Some(Prompt::Settings);
        self.ui.menu = None;
    }

    /// The window is closing: the session's length.
    pub fn usage_closing(&self) {
        if self.usage.is_some() {
            usage::closing();
        }
    }

    /// Ask the question once, when nothing else is being asked.
    pub(crate) fn ask_consent(&mut self) {
        if self.usage.is_some() && self.prefs.usage_sharing.is_none() && self.prompt.is_none() && !self.consent_asked {
            self.consent_asked = true;
            self.prompt = Some(Prompt::Consent);
        }
    }

    /// Record something the user did in the window (not a script or an agent driving it).
    pub(crate) fn record_ui(&self, event: Event) {
        if !self.session.origin.is_some_and(Origin::automated) {
            usage::record(event);
        }
    }

    /// Run `f` with the session's commands counted as coming from `origin`.
    pub fn with_origin<R>(&mut self, origin: Origin, f: impl FnOnce(&mut Self) -> R) -> R {
        let before = self.session.origin.replace(origin);
        let r = f(self);
        self.session.origin = before;
        r
    }

    pub(crate) fn open_link(&mut self, ctx: &egui::Context, url: &str) {
        match &self.services.open_url {
            Some(open) => {
                if let Err(e) = open(url) {
                    self.set_status(format!("Could not open the browser: {e}"), true);
                }
            }
            None => ctx.open_url(egui::OpenUrl::new_tab(url)),
        }
    }
}

/// The Settings window's icon: the gear, as on the title bar's Settings button.
const SETTINGS_ICON: &str = "settings";

/// A window's section heading after its first: the kit's heading with a section's gap above it,
/// so the group starts apart from the paragraph before it.
fn section(ui: &mut egui::Ui, title: &str) {
    ui.add_space(solidcraft_theme::Tokens::get().space_12);
    kit::section_heading(ui, title);
}

/// The consent question and Help › Settings, when open.
pub fn show(app: &mut ModelApp, ctx: &egui::Context) {
    match app.prompt {
        Some(Prompt::Consent) => {
            let (mut link, mut updates) = (false, None);
            let (closed, action) =
                crate::frame::window(egui::Id::unique("fc_consent"), "ui.help", "Help Improve FieldCraft", crate::frame::Width::Normal)
                    .actions("Share Anonymous Usage", "Don't Share", None)
                    .show(ctx, true, |ui| {
                        message(ui, QUESTION);
                        note(ui, DETAIL);
                        note(ui, "You can change this any time in Help › Settings.");
                        link = ui.link("What is sent (docs/privacy.md)").backing_on_hover_text(PRIVACY_URL).clicked();
                        if app.updates_available() {
                            updates = update_check(ui, app.prefs.update_check);
                        }
                    });
            if link {
                app.open_link(ctx, PRIVACY_URL);
            }
            if let Some(on) = updates {
                app.set_update_prefs(on, None);
            }
            if action == kit::Action::Ok {
                app.answer_consent(true);
            } else if dismissed(closed, action) {
                app.answer_consent(false);
            }
        }
        Some(Prompt::Settings) => {
            let available = app.usage_available();
            let (mut sharing, mut automation) = (app.prefs.usage_sharing == Some(true), app.prefs.usage_automation);
            let (mut link, mut updates, mut unskip, mut check) = (false, None, false, false);
            let mut units = None;
            let (closed, action) = crate::frame::window(egui::Id::unique("fc_settings"), SETTINGS_ICON, "Settings", crate::frame::Width::Normal)
                .actions("", "Close", None)
                .show(ctx, true, |ui| {
                    kit::section_heading(ui, "Units");
                    units = units_choice(ui, app.prefs.units.preset());
                    section(ui, "Interface");
                    kit::checkbox(ui, "Translucent browser", &mut app.ui.browser_float)
                        .backing_on_hover_text("Show the view through the browser card");
                    if app.updates_available() {
                        section(ui, "Updates");
                        updates = update_check(ui, app.prefs.update_check);
                        let (now, again) = update_actions(ui, &app.prefs.update_skip);
                        check = now.clicked();
                        unskip = again.is_some_and(|r| r.clicked());
                    }
                    section(ui, "Privacy");
                    ui.add_enabled_ui(available, |ui| {
                        kit::checkbox(ui, "Share anonymous usage", &mut sharing).backing_on_hover_text(QUESTION);
                    });
                    ui.add_enabled_ui(available && sharing, |ui| {
                        ui.indent("fc_usage_automation", |ui| {
                            kit::checkbox(ui, "Include scripts and AI agents (MCP)", &mut automation)
                                .backing_on_hover_text("Count the commands scripts, the control channel and AI agents run, apart from yours");
                        });
                    });
                    if !available {
                        note(
                            ui,
                            "This build or run sends nothing: it has no analytics key, or FIELDCRAFT_TELEMETRY=0, CI or automation turned it off.",
                        );
                    }
                    note(ui, DETAIL);
                    note(ui, "Crash reports and the log stay on this computer; Help › Report a Bug sends only what you submit.");
                    link = ui.link("What is sent (docs/privacy.md)").backing_on_hover_text(PRIVACY_URL).clicked();
                });
            if available && (sharing != (app.prefs.usage_sharing == Some(true)) || automation != app.prefs.usage_automation) {
                app.set_usage(sharing, automation);
            }
            if updates.is_some() || unskip {
                app.set_update_prefs(updates.unwrap_or(app.prefs.update_check), unskip.then(String::new));
            }
            if let Some(preset) = units {
                app.set_default_units(fieldcraft_model::DisplaySystem::new(preset));
            }
            if check {
                app.check_for_updates();
            }
            if link {
                app.open_link(ctx, PRIVACY_URL);
            }
            if dismissed(closed, action) {
                app.prompt = None;
            }
        }
        _ => {}
    }
}

/// The display units' drop-down (Settings › Units); the preset chosen when it changed.
fn units_choice(ui: &mut egui::Ui, current: fieldcraft_model::units::Preset) -> Option<fieldcraft_model::units::Preset> {
    use fieldcraft_model::units::Preset;
    let labels: Vec<&str> = Preset::ALL.iter().map(|p| p.label()).collect();
    let mut selected = Preset::ALL.iter().position(|p| *p == current).unwrap_or(0);
    kit::choice_row(ui, "Display units", "fc_settings_units", &mut selected, &labels);
    note(
        ui,
        "Numbers are shown and typed in these units; a model can keep its own (Model › Units). Models are stored in SI, so nothing in them changes.",
    );
    Preset::ALL.get(selected).copied().filter(|p| *p != current)
}

/// The automatic update check's checkbox; the new value when it was clicked.
/// The Updates section's actions: Check Now, and Offer It Again while a version is skipped
/// (`skip`). The window's ordinary buttons, as its Close is.
fn update_actions(ui: &mut egui::Ui, skip: &str) -> (egui::Response, Option<egui::Response>) {
    ui.horizontal(|ui| {
        let check = kit::button(ui, "Check Now", kit::Button::Normal).backing_on_hover_text("Help › Check for Updates");
        let again = (!skip.is_empty())
            .then(|| kit::button(ui, "Offer It Again", kit::Button::Normal).backing_on_hover_text(format!("Skipping FieldCraft {skip}")));
        (check, again)
    })
    .inner
}

fn update_check(ui: &mut egui::Ui, on: bool) -> Option<bool> {
    let mut checked = on;
    let changed = kit::checkbox(ui, "Check for updates automatically", &mut checked)
        .backing_on_hover_text(
            "At start-up, at most once a day, FieldCraft asks GitHub Releases for the latest version. Nothing installs until you choose Install.",
        )
        .changed();
    note(ui, "The check sends only the request (your IP address and the FieldCraft version reach GitHub); off means no check.");
    changed.then_some(checked)
}

#[cfg(test)]
mod tests {
    use fieldcraft_engine::Session;
    use fieldcraft_engine::usage::TEST_LOCK;
    use serde_json::json;

    use super::*;
    use crate::files::{PrefsStore, Services};

    fn host() -> UsageHost {
        UsageHost { platform: "desktop", os: "test", arch: "test", new_id: Box::new(|| Some("11111111-2222-4333-8444-555555555555".into())) }
    }

    fn app(file: &std::path::Path) -> ModelApp {
        ModelApp::new(Session::default()).with_services(Services { prefs: Some(PrefsStore::file(file.to_path_buf())), ..Default::default() })
    }

    fn temp(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("fieldcraft-privacy-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir.join("prefs.json")
    }

    /// A later section starts a section's gap below the paragraph before it, not at the row pitch.
    #[test]
    fn a_later_section_stands_apart_from_the_text_before_it() {
        let ctx = egui::Context::default();
        crate::theme::apply(&ctx, true);
        let mut above = egui::Rect::NOTHING;
        let mut out = ctx.run_ui(egui::RawInput::default(), |ui| {
            above = ui.label("The paragraph before the section.").rect;
            section(ui, "Spacing check");
        });
        out.textures_delta.clear();
        let heading = solidcraft_ui_registry::rect_of("dialog:heading:Spacing check").expect("the kit publishes its heading");
        let gap = heading.top() - above.bottom();
        assert!(gap >= solidcraft_theme::Tokens::get().space_12, "the heading is {gap} pt below the paragraph");
    }

    /// Check Now and Offer It Again are the window's ordinary buttons, as tall as its Close.
    #[test]
    fn update_actions_match_the_close_button() {
        let ctx = egui::Context::default();
        crate::theme::apply(&ctx, true);
        let (mut close, mut actions) = (egui::Rect::NOTHING, Vec::new());
        let mut out = ctx.run_ui(egui::RawInput::default(), |ui| {
            let (now, again) = update_actions(ui, "0.9.0");
            actions = std::iter::once(now).chain(again).map(|r| r.rect).collect();
            close = kit::button(ui, "Close", kit::Button::Normal).rect;
        });
        out.textures_delta.clear();
        assert_eq!(actions.len(), 2);
        for r in actions {
            assert!((r.height() - close.height()).abs() < 0.5, "{r:?} against Close {close:?}");
        }
    }

    /// The Settings window shows the gear the title bar's Settings button shows, not the
    /// style cube.
    #[test]
    fn settings_shows_the_gear() {
        let icon = crate::icons::SET.resolve(SETTINGS_ICON);
        assert!(!icon.placeholder, "{SETTINGS_ICON} has no artwork");
        assert_eq!(icon.master, "settings");
    }

    /// The first run asks once; nothing is recorded before the answer or after a no, and the
    /// answer is kept: the next run neither asks again nor forgets it.
    #[test]
    fn consent_is_asked_once_and_kept() {
        let _g = TEST_LOCK.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        usage::disable();
        let file = temp("consent");
        // Without the host's go-ahead (no key, automation, FIELDCRAFT_TELEMETRY=0) nothing asks.
        let mut a = app(&file);
        a.ask_consent();
        assert_eq!(a.prompt, None);

        let mut a = app(&file).with_usage(host());
        a.ask_consent();
        assert_eq!(a.prompt, Some(Prompt::Consent));
        let _ = a.run("model.info", json!({}));
        assert_eq!(usage::queued(), 0, "nothing before the answer");
        a.answer_consent(false);
        assert_eq!(a.prompt, None);
        let _ = a.run("model.info", json!({}));
        assert!(!usage::enabled() && usage::queued() == 0);
        a.ask_consent();
        assert_eq!(a.prompt, None, "asked once");

        let mut again = app(&file).with_usage(host());
        again.ask_consent();
        assert_eq!((again.prompt.clone(), again.prefs.usage_sharing), (None, Some(false)));

        // Settings › Privacy turns it on: the install id is made and kept, and commands count.
        again.set_usage(true, true);
        assert!(usage::enabled());
        let _ = again.run("model.info", json!({}));
        let events: Vec<Event> = std::iter::from_fn(usage::take_batch).flat_map(|b| b.events).collect();
        assert!(events.contains(&Event::AppStarted), "{events:?}");
        assert!(events.iter().any(|e| matches!(e, Event::CommandRun { command: "model.info", from: Origin::Ui, .. })), "{events:?}");
        let kept = crate::files::Prefs::load(Some(&PrefsStore::file(file.clone())));
        assert_eq!((kept.usage_sharing, kept.install_id.as_str()), (Some(true), "11111111-2222-4333-8444-555555555555"));
        usage::disable();
        let _ = std::fs::remove_dir_all(file.parent().unwrap());
    }

    /// A cancelled form and a command clicked before the model is ready are counted, by command
    /// id and the requirements' names; the control channel's are not the user's.
    #[test]
    fn cancelled_forms_and_unavailable_clicks_are_counted() {
        let _g = TEST_LOCK.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        usage::disable();
        let file = temp("events");
        let mut a = app(&file).with_usage(host());
        a.set_usage(true, true);
        let _ = std::iter::from_fn(usage::take_batch).count();
        let _ = a.start("solve.run");
        a.start("model.material.create").unwrap();
        a.close_dialog();
        a.with_origin(Origin::Script, |a| {
            let _ = a.start("solve.run");
        });
        let events: Vec<Event> = std::iter::from_fn(usage::take_batch).flat_map(|b| b.events).collect();
        let unavailable: Vec<_> = events.iter().filter(|e| matches!(e, Event::UnavailableClicked { command: "solve.run", .. })).collect();
        assert_eq!(unavailable.len(), 1, "{events:?}");
        assert!(events.contains(&Event::FormCancelled { command: "model.material.create" }), "{events:?}");
        assert_eq!(a.session.origin, Some(Origin::Ui));
        usage::disable();
        let _ = std::fs::remove_dir_all(file.parent().unwrap());
    }
}
