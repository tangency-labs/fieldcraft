//! The desktop app's side of updates (docs/updates.md): it runs the window's update jobs
//! ([`fieldcraft_ui_egui::updates`]) on background threads with [`fieldcraft_update`], and swaps
//! a staged release in after the window has closed ([`finish`]).
//!
//! Automatic checks are off in runs driven by automation (`--control`) and under CI, so test
//! and agent runs never reach GitHub on their own; `FIELDCRAFT_UPDATE_HOST` (a local test
//! server) turns them back on. Help › Check for Updates works in every run.

use std::cell::Cell;
use std::path::Path;
use std::sync::mpsc;

use fieldcraft_ui_egui::updates::{Event, Job, UpdateHost};
use fieldcraft_update::install::{self, Finish, Install};
use fieldcraft_update::net::{Client, Which};
use fieldcraft_update::stage::{self, Options};
use fieldcraft_update::{Release, Staged};

/// What the window needs for updates.
pub fn host(ctx: egui::Context, automation: bool) -> UpdateHost {
    let here = install::here();
    log::info!("updates: {}", here.describe());
    stage::clean_stale(&here);
    let env = |k: &str| std::env::var_os(k).is_some_and(|v| !v.is_empty());
    let automatic = env(fieldcraft_update::HOST_ENV) || !(automation || env("CI"));
    let (tx, events) = mpsc::channel::<Event>();
    let unsupported = here.unsupported();
    let start = move |job: Job| {
        let (tx, ctx, here) = (tx.clone(), ctx.clone(), here.clone());
        let spawned = std::thread::Builder::new().name("fieldcraft update".into()).spawn(move || {
            let send = |e: Event| {
                let _ = tx.send(e);
                ctx.request_repaint();
            };
            let client = Client::from_env();
            match job {
                Job::Check { manual } => send(Event::Checked { manual, result: client.check().map_err(|e| e.to_string()) }),
                Job::Stage(release) => send(Event::Staged(stage_release(&client, &release, &here, &send))),
                Job::Reinstall => {
                    let staged = client
                        .release(&Which::Version(fieldcraft_update::current()))
                        .map_err(|e| e.to_string())
                        .and_then(|r| stage_release(&client, &r, &here, &send));
                    send(Event::Staged(staged));
                }
            }
        });
        if let Err(e) = spawned {
            log::warn!("updates: no thread: {e}");
        }
    };
    UpdateHost { unsupported, automatic, start: Box::new(start), events: events.into(), now: Box::new(now), open_file: Box::new(open_file) }
}

/// Download, verify and stage `release`, reporting each whole percent.
fn stage_release(client: &Client, release: &Release, here: &Install, send: &dyn Fn(Event)) -> Result<Staged, String> {
    let downloads = fieldcraft_io::recovery::data_dir().map(|d| d.join("updates")).unwrap_or_else(|| std::env::temp_dir().join("fieldcraft-updates"));
    let last = Cell::new(u64::MAX);
    let progress = |done: u64, total: u64| {
        let step = done.saturating_mul(100).checked_div(total).unwrap_or(done / (1024 * 1024));
        if last.replace(step) != step {
            send(Event::Progress(done, total));
        }
    };
    let options = Options { key: fieldcraft_update::PUBLIC_KEY, downloads: &downloads, arch: std::env::consts::ARCH, progress: &progress };
    let staged = stage::stage(client, release, here, &options).map_err(|e| e.to_string());
    match &staged {
        Ok(s) => log::info!("updates: FieldCraft {} verified and staged at {}", s.version, s.payload.display()),
        Err(e) => log::warn!("updates: FieldCraft {} not staged: {e}", release.version),
    }
    staged
}

/// Swap `staged` in and start it: the window has closed, its autosave is closed.
pub fn finish(staged: &Staged) {
    match install::plan(staged, true) {
        Ok(Finish::Swap(plan)) => match stage::apply(&plan) {
            Ok(()) => log::info!("updates: FieldCraft {} installed ({})", staged.version, staged.install.describe()),
            Err(e) => log::error!("updates: FieldCraft {} not installed, the previous version is kept: {e}", staged.version),
        },
        Ok(Finish::Package { .. }) => {}
        Err(e) => log::error!("updates: {e}"),
    }
}

/// Seconds since the Unix epoch.
fn now() -> u64 {
    std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0)
}

/// Open a downloaded package with the system's handler for it (Linux: the software installer).
fn open_file(path: &Path) -> Result<(), String> {
    let opener = if cfg!(target_os = "macos") {
        "open"
    } else if cfg!(windows) {
        "explorer"
    } else {
        "xdg-open"
    };
    let mut child = std::process::Command::new(opener).arg(path).spawn().map_err(|e| format!("{opener}: {e}"))?;
    std::thread::spawn(move || {
        let _ = child.wait();
    });
    Ok(())
}
