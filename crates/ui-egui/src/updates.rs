//! In-app updates (docs/updates.md): the check at start-up (at most once a day, when the user
//! leaves it on) and from Help › Check for Updates, the quiet notice in the status bar
//! ("FieldCraft 0.5.9 is available — What's new · Install"), Help › Reinstall FieldCraft, and
//! the question that finishes an install. The host does the network and file work
//! ([`UpdateHost`], the desktop app's `update.rs`) on its own threads; the window only asks and
//! shows. Without a host (the browser, tests) nothing here runs and no network call is made.
//!
//! Nothing installs without the user's click: Install downloads and verifies the release, then
//! asks to restart. Restarting saves or asks about unsaved changes first, as opening another
//! model does; the swap happens after the window has closed.

use std::path::Path;

use egui::RichText;
use fieldcraft_update::{Check, Finish, Release, Staged};

use crate::ModelApp;
use crate::files::{Prompt, Replace};

/// What the host is asked to do, on its own thread.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Job {
    /// Ask for the latest release (`manual`: the user asked, so the answer is shown even when
    /// there is nothing new or the check fails).
    Check { manual: bool },
    /// Download, verify and stage this release.
    Stage(Release),
    /// Download, verify and stage the running version again (Help › Reinstall FieldCraft).
    Reinstall,
}

/// What the host reports back.
#[derive(Clone, Debug)]
pub enum Event {
    Checked {
        manual: bool,
        result: Result<Check, String>,
    },
    /// Bytes downloaded so far, and the total (0 when not known).
    Progress(u64, u64),
    Staged(Result<Staged, String>),
}

/// What the desktop app provides for updates.
pub struct UpdateHost {
    /// Why this copy cannot install updates itself (a development build, an unknown layout);
    /// `None` when it can. It still checks.
    pub unsupported: Option<String>,
    /// Automatic checks may run in this run (not under CI or automation, unless a test server
    /// is set). Help › Check for Updates works either way.
    pub automatic: bool,
    /// Start a job; its events come back on `events`.
    pub start: Box<dyn Fn(Job)>,
    pub events: crate::inbox::Inbox<Event>,
    /// Seconds since the Unix epoch.
    pub now: Box<dyn Fn() -> u64>,
    /// Open a downloaded package with the system's installer.
    pub open_file: Box<dyn Fn(&Path) -> Result<(), String>>,
}

/// The window's update state.
pub(crate) struct State {
    host: UpdateHost,
    /// A newer release the notice offers (not skipped, or asked for by hand).
    available: Option<Release>,
    /// What runs in the background, in words ("Checking for updates…").
    busy: Option<String>,
    /// The verified release, ready to finish, and how.
    staged: Option<(Staged, Finish)>,
    /// The automatic check of this run has been decided.
    auto_done: bool,
    /// Close the window to finish the install (sent with the next frame's context).
    close: bool,
}

impl ModelApp {
    /// Updates are available in this run (the desktop app).
    pub fn with_updates(mut self, host: UpdateHost) -> Self {
        self.updates = Some(State { host, available: None, busy: None, staged: None, auto_done: false, close: false });
        self
    }

    /// Whether this run checks for updates at all.
    pub fn updates_available(&self) -> bool {
        self.updates.is_some()
    }

    /// The automatic check, once per run, when it is on, due (a day since the last), allowed in
    /// this run, and the first run's question has been answered. Off means no network call.
    fn auto_check(&mut self) {
        let Some(u) = &mut self.updates else { return };
        if u.auto_done || self.prompt == Some(Prompt::Consent) {
            return;
        }
        u.auto_done = true;
        let now = (u.host.now)();
        if !u.host.automatic || !fieldcraft_update::due(self.prefs.update_check, self.prefs.update_checked_at, now) {
            return;
        }
        (u.host.start)(Job::Check { manual: false });
        self.prefs.update_checked_at = now;
        self.prefs.save(self.services.prefs.as_ref());
    }

    /// Help › Check for Updates.
    pub fn check_for_updates(&mut self) {
        self.ui.menu = None;
        let Some(u) = &mut self.updates else {
            self.set_status("This build does not check for updates", true);
            return;
        };
        u.busy = Some("Checking for updates…".into());
        (u.host.start)(Job::Check { manual: true });
        self.set_status("Checking for updates…", false);
    }

    /// Help › Reinstall FieldCraft: download, verify and swap in the running version again.
    pub fn reinstall(&mut self) {
        self.ui.menu = None;
        let Some(u) = &mut self.updates else { return };
        if let Some(why) = u.host.unsupported.clone() {
            self.set_status(format!("Cannot reinstall: {why}"), true);
            return;
        }
        u.busy = Some(format!("Downloading FieldCraft {}…", crate::VERSION));
        (u.host.start)(Job::Reinstall);
        self.set_status(format!("Downloading FieldCraft {} to reinstall it…", crate::VERSION), false);
    }

    /// The notice's Install: download and verify, then ask to restart. A copy that cannot
    /// install itself, or a release that is not signed, opens the release page instead.
    pub fn install_update(&mut self, ctx: &egui::Context) {
        let Some(u) = &mut self.updates else { return };
        let Some(release) = u.available.clone() else { return };
        let blocked = u.host.unsupported.clone().or_else(|| (!release.signed()).then(|| format!("FieldCraft {} is not signed", release.version)));
        if let Some(why) = blocked {
            self.set_status(format!("{why}: opening the download page"), true);
            self.open_link(ctx, &release.url);
            return;
        }
        u.busy = Some(format!("Downloading FieldCraft {}…", release.version));
        (u.host.start)(Job::Stage(release.clone()));
        self.set_status(format!("Downloading FieldCraft {}…", release.version), false);
    }

    /// Skip the offered version: the automatic check no longer offers it.
    pub fn skip_update(&mut self) {
        let Some(u) = &mut self.updates else { return };
        if let Some(r) = u.available.take() {
            self.prefs.update_skip = r.version.to_string();
            self.prefs.save(self.services.prefs.as_ref());
        }
    }

    /// Settings › Updates: automatic checks on or off, and the skipped version.
    pub fn set_update_prefs(&mut self, automatic: bool, skip: Option<String>) {
        self.prefs.update_check = automatic;
        if let Some(skip) = skip {
            self.prefs.update_skip = skip;
        }
        self.prefs.save(self.services.prefs.as_ref());
    }

    /// The release the notice offers.
    pub fn update_offered(&self) -> Option<&Release> {
        self.updates.as_ref()?.available.as_ref()
    }

    /// Restart to finish the install (or, `false`, finish it when the window closes). Unsaved
    /// changes are saved or asked about first, as when opening another model.
    pub fn answer_update(&mut self, restart: bool) {
        if self.prompt == Some(Prompt::Update) {
            self.prompt = None;
        }
        if !restart {
            if let Some(u) = &self.updates
                && let Some((s, Finish::Swap(_))) = &u.staged
            {
                let v = s.version.to_string();
                self.set_status(format!("FieldCraft {v} installs when you quit"), false);
            }
            return;
        }
        if self.has_unsaved_work() {
            self.prompt = Some(Prompt::Unsaved(Replace::Update));
            return;
        }
        self.close_to_update();
    }

    /// Close the window; the host swaps the staged release in after it has closed.
    pub(crate) fn close_to_update(&mut self) {
        if let Some(u) = &mut self.updates {
            u.close = true;
        }
    }

    /// The verified release to swap in once the window has closed (the desktop app's exit).
    pub fn take_staged_update(&mut self) -> Option<Staged> {
        let u = self.updates.as_mut()?;
        match &u.staged {
            Some((_, Finish::Swap(_))) => u.staged.take().map(|(s, _)| s),
            _ => None,
        }
    }

    /// Take the host's reports, run the automatic check when due, and close the window when
    /// an install is waiting on it.
    pub(crate) fn poll_updates(&mut self, ctx: &egui::Context) {
        self.auto_check();
        let Some(u) = &mut self.updates else { return };
        if std::mem::take(&mut u.close) {
            self.quit_requested = true;
            ctx.send_viewport_cmd(egui::ViewportCommand::Close);
            return;
        }
        let events: Vec<Event> = u.host.events.try_iter().collect();
        if u.busy.is_some() {
            ctx.request_repaint_after(std::time::Duration::from_millis(200));
        }
        for e in events {
            self.update_event(e);
        }
    }

    fn update_event(&mut self, e: Event) {
        let Some(u) = &mut self.updates else { return };
        match e {
            Event::Checked { manual, result } => {
                u.busy = None;
                match result {
                    Ok(check) if check.newer() && (manual || self.prefs.update_skip != check.latest.version.to_string()) => {
                        let v = check.latest.version.to_string();
                        u.available = Some(check.latest);
                        if manual {
                            self.set_status(format!("FieldCraft {v} is available"), false);
                        }
                    }
                    Ok(check) => {
                        // Nothing newer (or skipped): a notice from an earlier check goes.
                        u.available = None;
                        if manual {
                            self.set_status(format!("FieldCraft {} is up to date", check.current), false);
                        }
                    }
                    Err(e) if manual => self.set_status(format!("Could not check for updates: {e}"), true),
                    // An automatic check fails silently (logged).
                    Err(e) => log::info!("update check: {e}"),
                }
            }
            Event::Progress(done, total) => {
                if let Some(b) = &mut u.busy {
                    let head = b.split(" (").next().unwrap_or_default().to_string();
                    *b = match done.saturating_mul(100).checked_div(total) {
                        Some(percent) => format!("{head} ({percent}%)"),
                        None => format!("{head} ({} MB)", done / (1024 * 1024)),
                    };
                }
            }
            Event::Staged(Ok(staged)) => {
                u.busy = None;
                match fieldcraft_update::install::plan(&staged, true) {
                    Ok(finish) => {
                        u.staged = Some((staged, finish));
                        self.status = None;
                        self.prompt = Some(Prompt::Update);
                    }
                    Err(e) => self.set_status(format!("Cannot install: {e}"), true),
                }
            }
            Event::Staged(Err(e)) => {
                u.busy = None;
                self.set_status(format!("Update not installed: {e}"), true);
            }
        }
    }

    /// The web app: a new deployment is out; the status bar offers to reload into it.
    pub fn offer_reload(&mut self, reload: Box<dyn Fn()>) {
        self.reload = Some(reload);
    }

    /// Whether the status bar offers to reload into a new deployment.
    pub fn reload_offered(&self) -> bool {
        self.reload.is_some()
    }

    /// For the control channel's `ui.inspect`.
    pub fn update_json(&self) -> serde_json::Value {
        let Some(u) = &self.updates else { return serde_json::Value::Null };
        serde_json::json!({
            "automatic": self.prefs.update_check,
            "checked_at": self.prefs.update_checked_at,
            "skip": self.prefs.update_skip,
            "available": u.available.as_ref().map(|r| r.version.to_string()),
            "busy": u.busy,
            "staged": u.staged.as_ref().map(|(s, _)| s.version.to_string()),
            "unsupported": u.host.unsupported,
        })
    }
}

/// The status bar's notice: an update offered, or one downloading; in the browser, a new
/// deployment to reload into.
pub fn notice(app: &mut ModelApp, ui: &mut egui::Ui) {
    let t = crate::theme::Tokens::get(app.ui.dark);
    if let Some(reload) = &app.reload {
        let tip = if app.has_unsaved_work() { "Reloads the page; you are asked about the unsaved changes first" } else { "Reloads the page" };
        if ui.add(egui::Link::new(RichText::new("Reload to update").size(11.5))).on_hover_text(tip).clicked() {
            reload();
        }
        ui.label(RichText::new("A new version of FieldCraft is available —").size(11.5).color(t.text));
        ui.separator();
        return;
    }
    let Some(u) = &app.updates else { return };
    if let Some(busy) = &u.busy {
        ui.label(RichText::new(busy).size(11.5).color(t.text_dim));
        ui.separator();
        return;
    }
    let Some(release) = u.available.clone() else { return };
    // Right to left: the links, then the sentence.
    let small = |s: &str| RichText::new(s).size(11.5);
    let skip = ui.add(egui::Link::new(small("Skip"))).on_hover_text("Don't offer this version again (Help › Settings)").clicked();
    ui.label(small("·").color(t.text_dim));
    let install = ui.add(egui::Link::new(small("Install"))).on_hover_text("Download, verify and install it; you choose when to restart").clicked();
    ui.label(small("·").color(t.text_dim));
    let notes = ui.add(egui::Link::new(small("What's new"))).on_hover_text(release.url.as_str()).clicked();
    ui.label(small(&format!("FieldCraft {} is available —", release.version)).color(t.text));
    ui.separator();
    if notes {
        app.open_link(ui.ctx(), &release.url);
    }
    if install {
        app.install_update(ui.ctx());
    }
    if skip {
        app.skip_update();
    }
}

/// The question that finishes an install ([`Prompt::Update`]).
pub fn show(app: &mut ModelApp, ctx: &egui::Context) {
    if app.prompt != Some(Prompt::Update) {
        return;
    }
    let t = crate::theme::Tokens::get(app.ui.dark);
    let Some((staged, finish)) = app.updates.as_ref().and_then(|u| u.staged.clone()) else {
        app.prompt = None;
        return;
    };
    let reinstall = staged.version.to_string() == crate::VERSION;
    let (mut restart, mut later, mut open) = (false, false, false);
    crate::frame::window(ctx, "Update FieldCraft", crate::frame::Width::Normal).id(egui::Id::new("fc_update")).show(ctx, |ui| match &finish {
        Finish::Swap(_) => {
            ui.label(if reinstall {
                format!("FieldCraft {} was downloaded again and verified. Restart to reinstall it?", staged.version)
            } else {
                format!("FieldCraft {} is downloaded and verified. Restart now to finish installing it?", staged.version)
            });
            ui.label(
                RichText::new(format!("It replaces {}. Unsaved changes are saved or asked about first.", staged.install.describe()))
                    .color(t.text_dim)
                    .size(11.0),
            );
            ui.add_space(6.0);
            ui.horizontal(|ui| {
                restart = ui.button("Restart Now").clicked();
                later = ui.button("Later").on_hover_text("It installs when you quit FieldCraft").clicked();
            });
        }
        Finish::Package { file, command } => {
            ui.label(format!("FieldCraft {} is downloaded and verified:", staged.version));
            ui.label(RichText::new(file.display().to_string()).monospace().size(11.0));
            ui.label("Install it with your system's package manager:");
            let mut text = command.clone();
            ui.add(egui::TextEdit::singleline(&mut text).font(egui::TextStyle::Monospace).desired_width(f32::INFINITY));
            ui.add_space(6.0);
            ui.horizontal(|ui| {
                open = ui.button("Open with Installer").clicked();
                later = ui.button("Close").clicked();
            });
        }
    });
    if ctx.input(|i| i.key_pressed(egui::Key::Escape)) {
        later = true;
    }
    if open && let Finish::Package { file, .. } = &finish {
        let opened = app.updates.as_ref().map(|u| (u.host.open_file)(file));
        if let Some(Err(e)) = opened {
            app.set_status(format!("Could not open the installer: {e}"), true);
        }
        app.prompt = None;
    }
    if restart {
        app.answer_update(true);
    } else if later {
        app.answer_update(false);
    }
}

#[cfg(test)]
mod tests {
    use std::cell::RefCell;
    use std::rc::Rc;
    use std::sync::mpsc;

    use fieldcraft_engine::Session;
    use fieldcraft_update::{Install, Version};
    use serde_json::json;

    use super::*;
    use crate::files::{PrefsStore, Services};

    const DAY: u64 = fieldcraft_update::CHECK_INTERVAL;

    fn release(v: &str) -> Release {
        Release::from_github(&json!({"tag_name": format!("v{v}"), "html_url": "https://example.invalid/r", "assets": [
            {"name": "SHA256SUMS.txt", "browser_download_url": "https://example.invalid/s"},
            {"name": "SHA256SUMS.txt.minisig", "browser_download_url": "https://example.invalid/m"},
        ]}))
        .unwrap()
    }

    /// An app with a host that records the jobs it is given; the sender plays the host's thread.
    fn app(prefs: &Path, now: u64, automatic: bool) -> (ModelApp, Rc<RefCell<Vec<Job>>>, mpsc::Sender<Event>) {
        let jobs = Rc::new(RefCell::new(Vec::new()));
        let (tx, rx) = mpsc::channel();
        let j = jobs.clone();
        let host = UpdateHost {
            unsupported: None,
            automatic,
            start: Box::new(move |job| j.borrow_mut().push(job)),
            events: rx.into(),
            now: Box::new(move || now),
            open_file: Box::new(|_| Ok(())),
        };
        let a = ModelApp::new(Session::default())
            .with_services(Services { prefs: Some(PrefsStore::file(prefs.to_path_buf())), ..Default::default() })
            .with_updates(host);
        (a, jobs, tx)
    }

    fn temp(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("fieldcraft-updates-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir.join("prefs.json")
    }

    /// At most one automatic check a day, none when the user turned them off or the run does
    /// not allow them, and none while the first run's question is open.
    #[test]
    fn automatic_checks_are_daily_and_none_when_off() {
        let ctx = egui::Context::default();
        let file = temp("auto");
        let (mut a, jobs, _tx) = app(&file, 10 * DAY, true);
        a.poll_updates(&ctx);
        a.poll_updates(&ctx);
        assert_eq!(*jobs.borrow(), [Job::Check { manual: false }], "once per run");
        // The next run on the same day does not check; a day later it does.
        let (mut a, jobs, _tx) = app(&file, 10 * DAY + 3600, true);
        a.poll_updates(&ctx);
        assert!(jobs.borrow().is_empty());
        let (mut a, jobs, _tx) = app(&file, 11 * DAY, true);
        a.poll_updates(&ctx);
        assert_eq!(jobs.borrow().len(), 1);

        // Off in Settings: no network call, ever; Check for Updates still asks.
        a.set_update_prefs(false, None);
        let (mut a, jobs, _tx) = app(&file, 30 * DAY, true);
        a.poll_updates(&ctx);
        assert!(jobs.borrow().is_empty(), "off means no check");
        a.check_for_updates();
        assert_eq!(*jobs.borrow(), [Job::Check { manual: true }]);

        // A run that may not check automatically (CI, automation) does not.
        let (mut a, jobs, _tx) = app(&temp("ci"), 30 * DAY, false);
        a.poll_updates(&ctx);
        assert!(jobs.borrow().is_empty());

        // Not while the first run's question is open; after it is answered.
        let (mut a, jobs, _tx) = app(&temp("consent"), 30 * DAY, true);
        a.prompt = Some(Prompt::Consent);
        a.poll_updates(&ctx);
        assert!(jobs.borrow().is_empty());
        a.prompt = None;
        a.poll_updates(&ctx);
        assert_eq!(jobs.borrow().len(), 1);
        // Without a host (the browser, tests) there is nothing to call.
        let mut plain = ModelApp::new(Session::default());
        plain.poll_updates(&ctx);
        assert!(!plain.updates_available());
    }

    /// A newer release is offered (unless skipped, when the check was automatic); a manual
    /// check says when FieldCraft is up to date or why the check failed.
    #[test]
    fn checks_offer_newer_releases_and_skip_skipped_ones() {
        let ctx = egui::Context::default();
        let (mut a, _jobs, tx) = app(&temp("offer"), DAY, false);
        let current = Version::parse(crate::VERSION).unwrap();
        tx.send(Event::Checked { manual: false, result: Ok(Check { current: current.clone(), latest: release("99.0.0") }) }).unwrap();
        a.poll_updates(&ctx);
        assert_eq!(a.update_offered().map(|r| r.version.to_string()), Some("99.0.0".into()));
        a.skip_update();
        assert_eq!((a.update_offered(), a.prefs.update_skip.as_str()), (None, "99.0.0"));
        tx.send(Event::Checked { manual: false, result: Ok(Check { current: current.clone(), latest: release("99.0.0") }) }).unwrap();
        a.poll_updates(&ctx);
        assert_eq!(a.update_offered(), None, "skipped");
        tx.send(Event::Checked { manual: true, result: Ok(Check { current: current.clone(), latest: release("99.0.0") }) }).unwrap();
        a.poll_updates(&ctx);
        assert!(a.update_offered().is_some(), "asked for by hand");

        tx.send(Event::Checked { manual: true, result: Ok(Check { current: current.clone(), latest: release(crate::VERSION) }) }).unwrap();
        a.poll_updates(&ctx);
        assert_eq!(a.status.as_ref().map(|s| s.text.clone()), Some(format!("FieldCraft {} is up to date", crate::VERSION)));
        assert_eq!(a.update_offered(), None, "the earlier offer goes");
        tx.send(Event::Checked { manual: false, result: Err("offline".into()) }).unwrap();
        a.status = None;
        a.poll_updates(&ctx);
        assert!(a.status.is_none(), "an automatic check fails silently");
        tx.send(Event::Checked { manual: true, result: Err("offline".into()) }).unwrap();
        a.poll_updates(&ctx);
        assert!(a.status.as_ref().is_some_and(|s| s.error && s.text.contains("offline")));
    }

    /// Install downloads first; restarting asks about unsaved changes, then closes the window
    /// and hands the staged release to the host's exit.
    #[test]
    fn installing_asks_to_restart_and_saves_first() {
        let ctx = egui::Context::default();
        let dir = temp("install");
        let (mut a, jobs, tx) = app(&dir, DAY, false);
        let current = Version::parse(crate::VERSION).unwrap();
        tx.send(Event::Checked { manual: true, result: Ok(Check { current, latest: release("99.0.0") }) }).unwrap();
        a.poll_updates(&ctx);
        a.install_update(&ctx);
        assert_eq!(jobs.borrow().last(), Some(&Job::Stage(release("99.0.0"))));
        let staged = Staged {
            install: Install::TarGz { root: "/opt/fieldcraft".into() },
            version: Version::parse("99.0.0").unwrap(),
            dir: "/opt/.fieldcraft-update-99.0.0".into(),
            payload: "/opt/.fieldcraft-update-99.0.0/new/x".into(),
            entries: vec![],
        };
        tx.send(Event::Staged(Ok(staged.clone()))).unwrap();
        a.poll_updates(&ctx);
        assert_eq!(a.prompt, Some(Prompt::Update));

        // Unsaved work: asked first; Cancel keeps the window open.
        let _ = a.run("model.material.create", json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}}));
        a.answer_update(true);
        assert_eq!(a.prompt, Some(Prompt::Unsaved(Replace::Update)));
        a.answer_unsaved(crate::files::Answer::Cancel);
        a.poll_updates(&ctx);
        assert!(!a.quit_requested);
        a.prompt = Some(Prompt::Update);
        a.answer_update(true);
        a.answer_unsaved(crate::files::Answer::Discard);
        a.poll_updates(&ctx);
        assert!(a.quit_requested, "the window closes");
        assert_eq!(a.take_staged_update(), Some(staged));

        // A copy that cannot install itself opens the release page instead of downloading.
        let (mut a, jobs, tx) = app(&temp("dev"), DAY, false);
        a.updates.as_mut().unwrap().host.unsupported = Some("a development build".into());
        tx.send(Event::Checked { manual: true, result: Ok(Check { current: Version::parse(crate::VERSION).unwrap(), latest: release("99.0.0") }) })
            .unwrap();
        a.poll_updates(&ctx);
        a.install_update(&ctx);
        assert!(jobs.borrow().iter().all(|j| !matches!(j, Job::Stage(_))));
        assert!(a.status.as_ref().is_some_and(|s| s.text.contains("development build")));
        a.reinstall();
        assert!(jobs.borrow().iter().all(|j| *j != Job::Reinstall));
    }
}
