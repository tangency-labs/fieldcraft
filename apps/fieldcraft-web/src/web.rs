//! The browser shell: the web services of the UI (file pickers, downloads), dropped files, the
//! autosave, and eframe's web runner.
//!
//! Structure follows VectorCraft's `apps/vectorcraft-web/src/web.rs` and SolidCraft's
//! `apps/solidcraft-web/src/main.rs` (MIT OR Apache-2.0), see ATTRIBUTION.md.

use std::cell::{Cell, RefCell};
use std::path::Path;
use std::rc::Rc;

use fieldcraft_engine::Session;
use fieldcraft_engine::autosave::Autosaver;
use fieldcraft_ui_egui::files::Kind;
use fieldcraft_ui_egui::{Filter, ModelApp, Services};
use wasm_bindgen::JsCast as _;

use crate::storage::{self, Store};

const CANVAS_ID: &str = "fieldcraft_canvas";
const LOADING_ID: &str = "fieldcraft_loading";
const DEPTH_BITS: u8 = 24;

/// The page's folder for files from the computer and the files commands write.
const FILES_DIR: &str = "/files";

/// How often unsaved changes are autosaved for recovery (seconds), as on the desktop.
const AUTOSAVE_SECONDS: f64 = 30.0;

/// Files read from the computer (picked together or dropped together), waiting for the next
/// frame: (name, bytes or why they could not be read) each.
type Batch = Vec<(String, Result<Vec<u8>, String>)>;
type Inbox = Rc<RefCell<Vec<Batch>>>;

/// Largest file read from the computer: what the page can hold.
const MAX_READ_BYTES: f64 = fieldcraft_io::memory::MAX_TOTAL_BYTES as f64;

/// Why a file of `size` bytes is not read, if it is too large.
fn too_large(size: f64) -> Option<String> {
    (size > MAX_READ_BYTES)
        .then(|| format!("{} MB is more than this page holds ({} MB)", (size / 1048576.0).round(), fieldcraft_io::memory::MAX_TOTAL_BYTES >> 20))
}

/// Hand files read from the computer to the next frame.
fn deliver(inbox: &Inbox, ctx: &egui::Context, batch: Batch) {
    inbox.borrow_mut().push(batch);
    ctx.request_repaint();
}

/// Where a file named `name` is kept in the page.
fn page_path(name: &str) -> String {
    format!("{FILES_DIR}/{}", crate::safe_name(name))
}

thread_local! {
    /// The model has unsaved changes: leaving the page asks first.
    static UNSAVED: Cell<bool> = const { Cell::new(false) };
}

pub fn js_err(e: wasm_bindgen::JsValue) -> String {
    e.as_string().unwrap_or_else(|| format!("{e:?}"))
}

fn query() -> String {
    web_sys::window().and_then(|w| w.location().search().ok()).unwrap_or_default()
}

fn element(id: &str) -> Option<web_sys::Element> {
    web_sys::window()?.document()?.get_element_by_id(id)
}

pub fn start() {
    eframe::WebLogger::init(log::LevelFilter::Info).ok();
    // The browser's clock for the recovery entries (the plain wasm target has none).
    fieldcraft_io::recovery::set_clock(|| js_sys::Date::now().max(0.0) as u64);
    let Some(window) = web_sys::window() else { return };
    ask_before_unload(&window);
    crate::usage::send_when_hidden(&window);
    let Some(canvas) = element(CANVAS_ID).and_then(|e| e.dyn_into::<web_sys::HtmlCanvasElement>().ok()) else {
        log::error!("missing <canvas id=\"{CANVAS_ID}\">");
        return;
    };
    wasm_bindgen_futures::spawn_local(async move {
        let result = eframe::WebRunner::new()
            .start(
                canvas,
                web_options(),
                Box::new(|cc| {
                    let inbox = Inbox::default();
                    let mut app = ModelApp::new(Session::default()).with_services(services(inbox.clone(), cc.egui_ctx.clone()));
                    // Anonymous usage: asked once, recorded only when the user agreed (docs/privacy.md).
                    if let Some(host) = crate::usage::host() {
                        app = app.with_usage(host);
                    }
                    if let Some(rs) = &cc.wgpu_render_state {
                        let info = rs.adapter.get_info();
                        log::info!("fieldcraft-web: graphics {:?} {}, MSAA {}x", info.backend, info.name, rs.msaa_samples);
                        app.set_wgpu(rs, DEPTH_BITS, rs.msaa_samples, &cc.egui_ctx);
                    }
                    // The autosaves a closed or crashed page left, back in the recovery folder.
                    let store = Store::load();
                    let listed = app.session.execute("file.recovery.list", &serde_json::json!({}));
                    let recoverable = listed.ok().and_then(|v| v.as_array().map(Vec::len)).unwrap_or(0);
                    let autosaver = Autosaver::new(storage::dir()).map_err(|e| log::warn!("autosave unavailable: {e}")).ok();
                    if let Some(name) = crate::sample_in_query(&query()) {
                        app.open_sample(&name);
                    }
                    if recoverable > 0 {
                        app.set_status(format!("{recoverable} autosaved model(s) can be recovered: File › Recover…"), false);
                    }
                    let deploy = crate::deploy::Watch::start();
                    Ok(Box::new(Shell { app, inbox, store, autosaver, autosaved_at: 0.0, storage_error: false, sent_at: 0.0, deploy }))
                }),
            )
            .await;
        match result {
            Ok(()) => {
                if let Some(el) = element(LOADING_ID) {
                    el.remove();
                }
            }
            Err(e) => failed(&format!("FieldCraft could not start: {}", js_err(e))),
        }
    });
}

/// The runner's options: wgpu (WebGPU where available, else WebGL2; `?webgl` forces WebGL2).
fn web_options() -> eframe::WebOptions {
    let mut options = eframe::WebOptions { depth_buffer: DEPTH_BITS, ..Default::default() };
    if query().split(['?', '&']).any(|kv| kv == "webgl")
        && let eframe::egui_wgpu::WgpuSetup::CreateNew(create) = &mut options.wgpu_options.wgpu_setup
    {
        create.instance_descriptor.backends = eframe::wgpu::Backends::GL;
    }
    options
}

/// Replace the loading screen with why the app could not start.
fn failed(why: &str) {
    log::error!("{why}");
    let Some(el) = element(LOADING_ID) else { return };
    let Some(document) = web_sys::window().and_then(|w| w.document()) else { return };
    el.set_text_content(None);
    // Best effort: the console has the message too.
    let _ = el.class_list().add_1("failed");
    for line in
        [why, "FieldCraft needs a browser with WebGPU or WebGL2 (a current Chrome, Edge, Firefox or Safari) and graphics acceleration turned on."]
    {
        if let Ok(p) = document.create_element("p") {
            p.set_text_content(Some(line));
            let _ = el.append_child(&p);
        }
    }
}

/// The browser's "leave site? changes you made may not be saved" prompt while the model has
/// unsaved changes (a page that is closing can't show its own dialog).
fn ask_before_unload(window: &web_sys::Window) {
    let on_unload = wasm_bindgen::closure::Closure::<dyn FnMut(web_sys::BeforeUnloadEvent)>::new(|e: web_sys::BeforeUnloadEvent| {
        if UNSAVED.with(Cell::get) {
            e.prevent_default();
            e.set_return_value("unsaved");
        }
    });
    if let Err(e) = window.add_event_listener_with_callback("beforeunload", on_unload.as_ref().unchecked_ref()) {
        log::warn!("couldn't ask before leaving: {}", js_err(e));
    }
    // The listener lives as long as the page.
    on_unload.forget();
}

/// The UI's host services in the browser.
fn services(inbox: Inbox, ctx: egui::Context) -> Services {
    Services {
        // The browser's file picker: the files arrive in a later frame (the picker and the reads
        // are asynchronous), so the call itself picks nothing.
        pick_open: Some(Box::new(move |filters: &[Filter]| {
            let (inbox, ctx) = (inbox.clone(), ctx.clone());
            let dialog = filters.iter().fold(rfd::AsyncFileDialog::new(), |d, f| d.add_filter(f.name, f.extensions));
            wasm_bindgen_futures::spawn_local(async move {
                let Some(files) = dialog.pick_files().await else { return };
                let mut batch = Vec::with_capacity(files.len());
                for f in files {
                    let read = match too_large(f.inner().size()) {
                        Some(why) => Err(why),
                        None => Ok(f.read().await),
                    };
                    batch.push((f.file_name(), read));
                }
                deliver(&inbox, &ctx, batch);
            });
            None
        })),
        // Saving writes into the page; the shell hands each written file out as a download.
        pick_save: Some(Box::new(|name: &str, _: &[Filter]| Some(page_path(name)))),
        // The recent list holds the page's files, which go with the page.
        recent_path: None,
        // Preferences are kept in the browser's local storage.
        prefs: Some(crate::usage::prefs_store()),
        // egui opens links itself in the browser.
        open_url: None,
        ..Default::default()
    }
}

/// The app with the browser's files around it.
struct Shell {
    app: ModelApp,
    inbox: Inbox,
    store: Store,
    /// `None` when autosave could not start.
    autosaver: Option<Autosaver>,
    /// App time of the last autosave check.
    autosaved_at: f64,
    /// The last autosave did not fit in the browser's storage (said once until one fits).
    storage_error: bool,
    /// App time usage events were last sent.
    sent_at: f64,
    /// Looks for a new deployment of the app.
    deploy: crate::deploy::Watch,
}

impl Shell {
    /// Files read from the computer go into the page's folder, then open: a model with the
    /// design and meshes it refers to picked or dropped with it, else each file in turn.
    fn open_arrived(&mut self) {
        let batches: Vec<_> = self.inbox.borrow_mut().drain(..).collect();
        for batch in batches {
            let mut paths = Vec::with_capacity(batch.len());
            for (name, read) in batch {
                let path = page_path(&name);
                match read.and_then(|bytes| fieldcraft_io::fs::insert(Path::new(&path), bytes).map_err(|e| e.to_string())) {
                    Ok(()) => paths.push(path),
                    Err(e) => self.app.set_status(format!("{name}: {e}"), true),
                }
            }
            let model = paths.iter().find(|p| Kind::of(p) == Some(Kind::Model)).cloned();
            match model {
                // Its geometry files open with it, from the same folder.
                Some(model) => self.app.open_path(&model),
                None => {
                    for p in &paths {
                        self.app.open_path(p);
                    }
                }
            }
        }
    }

    /// What commands wrote: autosaves to the browser's storage, everything else downloaded.
    fn hand_out_writes(&mut self) {
        for (path, bytes) in fieldcraft_io::fs::take_writes() {
            if path.parent() == Some(storage::dir()) {
                match self.store.put(&path, &bytes) {
                    Ok(()) => self.storage_error = false,
                    Err(e) if !self.storage_error => {
                        self.storage_error = true;
                        self.app.set_status(format!("Autosave: {e}. Save the model to keep it (File › Save downloads it)"), true);
                    }
                    Err(e) => log::warn!("autosave: {e}"),
                }
            } else if let Err(e) = download(&path, &bytes) {
                self.app.set_status(format!("{}: the download did not start ({e})", path.display()), true);
            }
        }
        self.store.sync_removals();
    }

    fn autosave(&mut self, ctx: &egui::Context) {
        let Some(a) = &mut self.autosaver else { return };
        if self.app.now - self.autosaved_at >= AUTOSAVE_SECONDS {
            self.autosaved_at = self.app.now;
            if let Err(e) = a.save(&self.app.session) {
                log::warn!("autosave: {e}");
            }
        }
        ctx.request_repaint_after(std::time::Duration::from_secs_f64(AUTOSAVE_SECONDS));
    }
}

impl eframe::App for Shell {
    fn logic(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.open_arrived();
        self.app.logic(ctx);
        self.autosave(ctx);
        if self.app.now - self.sent_at >= crate::usage::INTERVAL {
            self.sent_at = self.app.now;
            crate::usage::send();
        }
        UNSAVED.with(|u| u.set(self.app.session.is_dirty() && !self.app.session.journal.is_empty()));
        if self.app.now - self.deploy.looked_at >= crate::deploy::INTERVAL || self.app.now < self.deploy.looked_at {
            self.deploy.looked_at = self.app.now;
            self.deploy.look();
        }
        if self.deploy.found() && !self.app.reload_offered() {
            self.app.offer_reload(Box::new(crate::deploy::reload));
        }
        ctx.request_repaint_after(std::time::Duration::from_secs_f64(crate::deploy::INTERVAL));
    }

    /// Dropped files are read here (browsers read them asynchronously) and arrive through the
    /// inbox; the app never sees their names as paths.
    fn raw_input_hook(&mut self, ctx: &egui::Context, raw: &mut egui::RawInput) {
        let dropped = std::mem::take(&mut raw.dropped_files);
        if !dropped.is_empty() {
            let (inbox, ctx) = (self.inbox.clone(), ctx.clone());
            wasm_bindgen_futures::spawn_local(async move {
                let mut batch = Vec::with_capacity(dropped.len());
                for f in dropped {
                    let name = f.path().to_string_lossy().to_string();
                    let read = match f.web_file().and_then(|w| too_large(w.size())) {
                        Some(why) => Err(why),
                        None => f.bytes_async().await,
                    };
                    batch.push((name, read));
                }
                deliver(&inbox, &ctx, batch);
            });
        }
        self.app.raw_input_hook(raw);
    }

    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        self.app.ui(ui);
        // Once a frame, after everything that writes: autosaves, saves and decks.
        self.hand_out_writes();
    }
}

/// Hand `bytes` out as a download named after the last part of `path`.
fn download(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("fieldcraft").to_string();
    let window = web_sys::window().ok_or("no window")?;
    let document = window.document().ok_or("no document")?;
    let parts = js_sys::Array::of1(&js_sys::Uint8Array::from(bytes));
    let options = web_sys::BlobPropertyBag::new();
    options.set_type("application/octet-stream");
    let blob = web_sys::Blob::new_with_u8_array_sequence_and_options(&parts, &options).map_err(js_err)?;
    let url = web_sys::Url::create_object_url_with_blob(&blob).map_err(js_err)?;
    let a: web_sys::HtmlAnchorElement = document.create_element("a").map_err(js_err)?.dyn_into().map_err(|_| "not a link")?;
    a.set_href(&url);
    a.set_download(&name);
    a.style().set_property("display", "none").map_err(js_err)?;
    document.body().ok_or("no body")?.append_child(&a).map_err(js_err)?;
    a.click();
    a.remove();
    // Revoked once the click has been dispatched; the download keeps its own reference.
    let revoke = wasm_bindgen::closure::Closure::once_into_js(move || {
        let _ = web_sys::Url::revoke_object_url(&url);
    });
    // The download has started: a timer that can't be set leaves the URL allocated, no more.
    let _ = window.set_timeout_with_callback_and_timeout_and_arguments_0(revoke.unchecked_ref(), 10_000);
    Ok(())
}
