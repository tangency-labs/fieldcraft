//! FieldCraft desktop app.
//!
//! Usage: `fieldcraft [--sample [NAME]] [--control PORT] [FILE…]` (files open as if dropped on the
//! window; NAME is a sample of `fieldcraft_engine::sample`, `parts` by default)
//!
//! `--sample` opens a sample model (`fieldcraft_engine::sample`). `--control PORT` (or
//! `FIELDCRAFT_CONTROL_PORT`) starts a localhost JSON-lines control server:
//! `{"id":1,"method":"engine.execute","params":{"command":"model.info","params":{}}}`. See
//! `fieldcraft_ui_egui::control` and docs/control-protocol.md.
#![cfg_attr(all(target_os = "windows", not(debug_assertions)), windows_subsystem = "windows")]
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
// Unsafe code is denied (Cargo.toml): only the Windows exception filter allows it (native_crash.rs).

mod control_server;
#[cfg(windows)]
mod native_crash;
mod terminate;
mod update;
mod usage;

use fieldcraft_engine::autosave::{self, Autosaver};
use fieldcraft_engine::{Session, crash, logfile};
use fieldcraft_ui_egui::{Filter, ModelApp, Services};

const USAGE: &str =
    "usage: fieldcraft [--sample [NAME]] [--control PORT] [FILE…]  (NAME: a sample, docs/examples.md; FILE: .fieldcraft, .solidcraft or .3mf)";
const DEPTH_BITS: u8 = 24;
/// The window draws its own caption buttons where the OS ones would sit apart from the app bar.
const CUSTOM_TITLEBAR: bool = !cfg!(target_os = "macos");
const MSAA: u16 = 4;
/// The reverse-DNS application id: the Wayland app_id / X11 class, matched to the packaged
/// `.desktop` entry (packaging/linux), and the macOS bundle id (packaging/macos/Info.plist.in).
const APP_ID: &str = "io.github.jyviko.fieldcraft";

/// The window and Dock icon: the 256 px export of the application icon (assets/app-icon).
const ICON_PNG: &[u8] = include_bytes!("../../../assets/app-icon/icon.iconset/icon_128x128@2x.png");

/// The application icon for the window; eframe's default icon if the PNG cannot be decoded.
fn app_icon() -> egui::IconData {
    eframe::icon_data::from_png_bytes(ICON_PNG).unwrap_or_else(|e| {
        log::warn!("app icon: {e}");
        egui::IconData::default()
    })
}

/// How often unsaved changes are autosaved for recovery (seconds).
const AUTOSAVE_SECONDS: f64 = 30.0;

/// The first window's inner size and the smallest the window may be (points).
const WINDOW_SIZE: egui::Vec2 = egui::vec2(1400.0, 900.0);
const MIN_WINDOW_SIZE: egui::Vec2 = egui::vec2(960.0, 600.0);
/// What the taskbar, dock or panels and the title bar may take off each side of the monitor: winit
/// reports the monitor's size but not its work area (points).
const WORK_AREA_RESERVE: f32 = 80.0;

/// Where the first window goes on a monitor of `monitor` points when the default size does not fit
/// in its work area (#601): as large as the monitor leaves after the reserve (never below the
/// minimum), centred across the monitor and at its top, so the system's cascade cannot push its
/// bottom under the taskbar. `None` when the default fits or the monitor's size is unknown. The
/// position is the monitor's own (winit gives no monitor origin): the primary monitor's, where a
/// first window opens.
fn fitted_window(monitor: egui::Vec2) -> Option<egui::Rect> {
    if !(monitor.is_finite() && monitor.min_elem() > 0.0) {
        return None;
    }
    let room = monitor - egui::Vec2::splat(WORK_AREA_RESERVE);
    let size = WINDOW_SIZE.min(room).max(MIN_WINDOW_SIZE);
    let at = egui::pos2(((monitor.x - size.x) * 0.5).max(0.0), 0.0);
    (size != WINDOW_SIZE).then(|| egui::Rect::from_min_size(at, size))
}

/// On a monitor too small for the default size, the first window starts maximized (the system
/// knows the work area exactly), restoring to a size and place that fit.
fn fit_first_window(ctx: &egui::Context) {
    let Some(fit) = ctx.input(|i| i.viewport().monitor_size).and_then(fitted_window) else { return };
    log::info!("window: {} × {} points do not fit the monitor; maximized, restoring to {fit:?}", WINDOW_SIZE.x, WINDOW_SIZE.y);
    ctx.send_viewport_cmd(egui::ViewportCommand::OuterPosition(fit.min));
    ctx.send_viewport_cmd(egui::ViewportCommand::InnerSize(fit.size()));
    ctx.send_viewport_cmd(egui::ViewportCommand::Maximized(true));
}

struct App {
    app: ModelApp,
    /// Autosave for crash recovery (`None` when the recovery folder is unavailable).
    autosaver: Option<Autosaver>,
    /// App time of the last autosave check.
    autosaved_at: f64,
    /// Frames drawn so far, up to 2 (the second means the first was presented: the app is up).
    frames: u8,
    /// Sends the anonymous usage events, when this run may (docs/privacy.md).
    usage: Option<usage::Sender>,
}

impl App {
    fn autosave(&mut self) {
        if let Some(a) = &mut self.autosaver
            && let Err(e) = a.save(&self.app.session)
        {
            log::warn!("autosave: {e}");
        }
    }
}

impl eframe::App for App {
    fn logic(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        if self.frames < 2 {
            if self.frames == 0 {
                fit_first_window(ctx);
            }
            self.frames += 1;
            // A second frame right away, which means the first was presented: the app is up.
            ctx.request_repaint();
            if self.frames == 2 {
                crash::phase(crash::RUNNING);
            }
        }
        self.app.logic(ctx);
        if self.app.now - self.autosaved_at >= AUTOSAVE_SECONDS {
            self.autosaved_at = self.app.now;
            self.autosave();
        }
        if self.autosaver.is_some() {
            ctx.request_repaint_after(std::time::Duration::from_secs_f64(AUTOSAVE_SECONDS));
        }
    }
    fn on_exit(&mut self) {
        if let Some(a) = self.autosaver.take()
            && let Err(e) = a.close(&self.app.session)
        {
            log::warn!("autosave: {e}");
        }
        self.app.usage_closing();
        if let Some(u) = self.usage.take() {
            u.finish();
        }
        // An update the user chose to install: swapped in now that the window has closed.
        if let Some(staged) = self.app.take_staged_update() {
            update::finish(&staged);
        }
        log::info!("FieldCraft closed");
        crash::end_session();
    }
    fn raw_input_hook(&mut self, _ctx: &egui::Context, raw: &mut egui::RawInput) {
        self.app.raw_input_hook(raw);
    }
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        self.app.ui(ui);
    }
}

/// Native file pickers and the recent-files file, for the UI.
fn services() -> Services {
    let dialog = |filters: &[Filter]| filters.iter().fold(rfd::FileDialog::new(), |d, f| d.add_filter(f.name, f.extensions));
    Services {
        pick_open: Some(Box::new(move |filters: &[Filter]| dialog(filters).pick_file().map(|p| p.to_string_lossy().to_string()))),
        pick_save: Some(Box::new(move |name: &str, filters: &[Filter]| {
            dialog(filters).set_file_name(name).save_file().map(|p| p.to_string_lossy().to_string())
        })),
        // Beside the recovery folder in the per-user data folder.
        recent_path: fieldcraft_io::recovery::default_dir().and_then(|d| d.parent().map(|p| p.join("recent.json"))),
        prefs: fieldcraft_io::recovery::default_dir()
            .and_then(|d| d.parent().map(|p| fieldcraft_ui_egui::files::PrefsStore::file(p.join("prefs.json")))),
        open_url: Some(Box::new(open_url)),
    }
}

/// Open a web link in the default browser, with the system's own opener: `open` (macOS),
/// `xdg-open` (Linux and the BSDs), `explorer` (Windows). The link is one argument, never
/// passed through a shell. Only https links are opened.
fn open_url(url: &str) -> Result<(), String> {
    if !url.starts_with("https://") {
        return Err(format!("not a web link: {}", url.chars().take(80).collect::<String>()));
    }
    let opener = if cfg!(target_os = "macos") {
        "open"
    } else if cfg!(windows) {
        "explorer"
    } else {
        "xdg-open"
    };
    let mut child = std::process::Command::new(opener).arg(url).spawn().map_err(|e| format!("{opener}: {e}"))?;
    // Reaped off the window thread; its exit status says nothing (`explorer` exits with 1 even
    // when it opened the link).
    std::thread::spawn(move || {
        let _ = child.wait();
    });
    Ok(())
}

/// Crash reports and the log file (docs/crash-reports.md): reports of earlier runs that crashed
/// come back, for the window to offer. Both live in the app's data folder (a temporary one in
/// runs driven by automation); `FIELDCRAFT_CRASH_DIR` moves the crash folder.
fn start_reporting() -> Vec<crash::Summary> {
    // First the crash reports: an unclean exit's report takes the dead run's lines from the log.
    let crashed = crash::default_dir().map(|dir| crash::install(dir, crash::Build::current())).unwrap_or_default();
    match fieldcraft_io::recovery::data_dir().map(|d| logfile::install(&d.join("logs"), log::LevelFilter::Info)) {
        Some(Ok(path)) => eprintln!("fieldcraft: log {}", path.display()),
        Some(Err(e)) => eprintln!("fieldcraft: no log file: {e}"),
        None => {}
    }
    let b = crash::Build::current();
    log::info!(
        "FieldCraft {} started ({}), crash reports in {}",
        b.version,
        b.version_line(),
        crash::default_dir().map(|d| d.display().to_string()).unwrap_or_default()
    );
    crashed
}

/// egui-wgpu is about to start a step of getting the graphics device and the window's surface
/// up: the crash journal keeps it, and the adapter once it is picked (#299).
fn startup_step(step: &'static str, detail: &str) {
    if step == eframe::egui_wgpu::STEP_ADAPTER_CHOSEN {
        crash::set_environment(&format!("Renderer: {detail}"));
    }
    crash::phase(step);
}

/// Command-line options.
#[derive(Debug, Default, PartialEq)]
struct Args {
    /// The sample to open (`--sample [NAME]`).
    sample: Option<String>,
    control: Option<u16>,
    /// Files to open (models, designs, tagged meshes).
    files: Vec<String>,
}

fn parse_args(args: impl IntoIterator<Item = String>, env_port: Option<String>) -> Result<Option<Args>, String> {
    let mut out = Args { control: env_port.and_then(|p| p.trim().parse().ok()), ..Default::default() };
    let mut args = args.into_iter().peekable();
    while let Some(a) = args.next() {
        match a.as_str() {
            "--sample" => {
                // An optional sample name follows.
                let name = args.peek().filter(|a| fieldcraft_engine::sample::find(a).is_some()).cloned();
                if name.is_some() {
                    args.next();
                }
                out.sample = Some(name.unwrap_or_else(|| "parts".into()));
            }
            "--control" => {
                let port = args.next().ok_or("--control needs a port")?;
                out.control = Some(port.parse().map_err(|_| format!("--control: `{port}` is not a port number"))?);
            }
            "--version" | "-V" => {
                println!("fieldcraft {}", fieldcraft_ui_egui::VERSION);
                return Ok(None);
            }
            "--help" | "-h" => {
                println!("{USAGE}");
                return Ok(None);
            }
            other if other.starts_with('-') => return Err(format!("unknown argument `{other}`\n{USAGE}")),
            file => out.files.push(file.to_string()),
        }
    }
    Ok(Some(out))
}

fn main() -> eframe::Result {
    // The lattice Boltzmann GPU backend, for `solve.fluid` with `backend: "gpu"`.
    fieldcraft_lbm_gpu::register();
    fieldcraft_mesh_gpu::register();
    // Pictures commands draw themselves (results.export.animation), with the viewport's pipeline.
    fieldcraft_ui_egui::register_picture_drawer();
    let args = match parse_args(std::env::args().skip(1), std::env::var("FIELDCRAFT_CONTROL_PORT").ok()) {
        Ok(Some(a)) => a,
        Ok(None) => return Ok(()),
        Err(e) => {
            eprintln!("fieldcraft: {e}");
            std::process::exit(2);
        }
    };
    if args.control.is_some() {
        // A run driven by automation keeps its autosaves (and recent files and preferences, which
        // live beside them), its crash reports and its log in a temporary folder, never the user's.
        fieldcraft_io::recovery::use_automation_dir();
    }
    let crashed = start_reporting();
    // A crash in native code (a graphics driver) says where it happened: the exception and the
    // module on Windows, and on every system the startup step the app was in.
    #[cfg(windows)]
    native_crash::install();
    eframe::egui_wgpu::set_startup_hook(startup_step);
    crash::phase("window");
    let mut options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("FieldCraft")
            .with_app_id(APP_ID)
            .with_icon(app_icon())
            .with_inner_size(WINDOW_SIZE)
            .with_min_inner_size(MIN_WINDOW_SIZE)
            .with_drag_and_drop(true)
            // The app bar is the title bar (as in the other Craft apps): on macOS under the
            // traffic lights, elsewhere with its own caption buttons and resize zones.
            .with_decorations(!CUSTOM_TITLEBAR)
            .with_fullsize_content_view(true)
            .with_titlebar_shown(false)
            .with_title_shown(false),
        depth_buffer: DEPTH_BITS,
        multisampling: MSAA,
        ..Default::default()
    };
    // Match headless rendering and compute: environment overrides cannot enable another API.
    if let eframe::egui_wgpu::WgpuSetup::CreateNew(create) = &mut options.wgpu_options.wgpu_setup {
        create.instance_descriptor.backends = fieldcraft_gpu::select::backends().map_err(|e| eframe::Error::AppCreation(Box::new(e)))?;
        create.instance_descriptor.backend_options = Default::default();
    }
    eframe::run_native(
        "FieldCraft",
        options,
        Box::new(move |cc| {
            crash::phase("app setup");
            let mut session = Session::default();
            // The window's model: its commands go into the crash journal.
            session.crash_journal = crash::enabled();
            let mut app = ModelApp::new(session).with_services(services());
            // Anonymous usage: asked once, recorded only when the user agreed (docs/privacy.md).
            let usage = match usage::host(args.control.is_some()) {
                Some(host) => {
                    app = app.with_usage(host);
                    usage::start()
                }
                None => None,
            };
            // Updates (docs/updates.md): checked at most daily when on, installed only on a click.
            fieldcraft_engine::enable_update_checks(true);
            app = app.with_updates(update::host(cc.egui_ctx.clone(), args.control.is_some()));
            app.custom_titlebar = CUSTOM_TITLEBAR;
            app.integrated_titlebar = cfg!(target_os = "macos");
            if let Some(rs) = &cc.wgpu_render_state {
                let info = rs.adapter.get_info();
                log::info!("graphics {:?} {} ({}), MSAA {}x (requested {MSAA}x)", info.backend, info.name, info.driver, rs.msaa_samples);
                crash::set_environment(&format!("Renderer: {:?} {} ({}), MSAA {}x", info.backend, info.name, info.driver, rs.msaa_samples));
                crash::phase("viewport pipelines (shader compile)");
                app.set_wgpu(rs, DEPTH_BITS, rs.msaa_samples, &cc.egui_ctx);
                // The graphics driver is loaded now; it may have put its own exception filter
                // in front of ours.
                #[cfg(windows)]
                native_crash::reinstall();
            }
            // `kill` or Ctrl-C closes the window as Close does: autosaved, no crash (#248).
            let ctx = cc.egui_ctx.clone();
            if let Err(e) = terminate::watch(move || {
                ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                ctx.request_repaint();
            }) {
                log::warn!("termination signals: {e}");
            }
            if let Some(port) = args.control {
                app = app.with_control(control_server::start(port, cc.egui_ctx.clone()));
            }
            if let Some(name) = &args.sample {
                app.open_sample(name);
            }
            for f in &args.files {
                app.open_path(f);
            }
            let autosaver = autosave::default_dir().and_then(|dir| {
                let listed = app.session.execute("file.recovery.list", &serde_json::json!({"dir": dir}));
                let recoverable = listed.ok().and_then(|v| v.as_array().map(Vec::len)).unwrap_or(0);
                if recoverable > 0 {
                    app.set_status(format!("{recoverable} autosaved model(s) can be recovered: File › Recover…"), false);
                }
                Autosaver::new(&dir).map_err(|e| log::warn!("autosave unavailable: {e}")).ok()
            });
            app.announce_crashes(&crashed);
            crash::phase("first frame");
            Ok(Box::new(App { app, autosaver, autosaved_at: 0.0, frames: 0, usage }))
        }),
    )
}

#[cfg(test)]
mod tests {
    #[test]
    fn the_app_icon_decodes() {
        let icon = super::app_icon();
        assert_eq!((icon.width, icon.height), (256, 256));
        assert_eq!(icon.rgba.len(), 256 * 256 * 4);
    }

    use super::*;

    #[test]
    fn the_first_window_fits_the_work_area() {
        let fit = |w, h| fitted_window(egui::vec2(w, h)).map(|r| (r.min.x, r.min.y, r.width(), r.height()));
        // Large monitors keep the default size and the system's placement; an unknown monitor too.
        assert_eq!(fit(2560.0, 1440.0), None);
        assert_eq!(fit(1536.0, 1080.0), None);
        assert_eq!(fit(0.0, 0.0), None);
        assert_eq!(fit(f32::NAN, 900.0), None);
        // 3440 × 1440 at 150 % (#601): 2293 × 960 points, a 912-point work area. Restored at the
        // top, its bottom (880) stays above the taskbar wherever the cascade would have put it.
        assert_eq!(fit(2293.0, 960.0), Some((446.5, 0.0, 1400.0, 880.0)));
        // 1920 × 1080 at 150 % and at 125 %.
        assert_eq!(fit(1280.0, 720.0), Some((40.0, 0.0, 1200.0, 640.0)));
        assert_eq!(fit(1536.0, 864.0), Some((68.0, 0.0, 1400.0, 784.0)));
        // Too narrow only: still at the top.
        assert_eq!(fit(1440.0, 1100.0), Some((40.0, 0.0, 1360.0, 900.0)));
        // Never below the minimum size.
        assert_eq!(fit(800.0, 600.0), Some((0.0, 0.0, 960.0, 600.0)));
    }

    fn parse(a: &[&str], env: Option<&str>) -> Result<Option<Args>, String> {
        parse_args(a.iter().map(|s| s.to_string()), env.map(str::to_string))
    }

    #[test]
    fn arguments() {
        assert_eq!(parse(&[], None), Ok(Some(Args::default())));
        assert_eq!(
            parse(&["--sample", "--control", "7070"], None),
            Ok(Some(Args { sample: Some("parts".into()), control: Some(7070), files: vec![] }))
        );
        assert_eq!(parse(&[], Some("9000")), Ok(Some(Args { control: Some(9000), ..Default::default() })));
        assert_eq!(parse(&[], Some("nope")), Ok(Some(Args::default())));
        assert!(parse(&["--control"], None).is_err());
        assert!(parse(&["--control", "99999"], None).is_err());
        assert!(parse(&["--frobnicate"], None).is_err());
        assert_eq!(parse(&["--sample", "impact"], None).map(|a| a.and_then(|a| a.sample)), Ok(Some("impact".into())));
        let a = parse(&["--sample", "x.fieldcraft"], None).unwrap().unwrap();
        assert_eq!((a.sample.as_deref(), a.files), (Some("parts"), vec!["x.fieldcraft".to_string()]));
        assert_eq!(parse(&["a.fieldcraft", "--sample"], None).map(|a| a.map(|a| a.files)), Ok(Some(vec!["a.fieldcraft".to_string()])));
        assert_eq!(parse(&["--help"], None), Ok(None));
    }
}
