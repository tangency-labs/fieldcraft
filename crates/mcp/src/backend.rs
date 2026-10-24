//! Where tool calls end up: an in-process headless session, or a running app's control channel
//! (docs/control-protocol.md: `engine.execute`, `engine.commands`, `ui.view`, `ui.camera`,
//! `ui.render`, `ui.capture`, `ui.select`, `ui.selection`).

use std::io::{BufRead, BufReader, Read, Write};
use std::net::{TcpStream, ToSocketAddrs};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Duration;

use fieldcraft_engine::Session;
use fieldcraft_engine::view::{Aim, Palette, Section, Shown, Target};
use fieldcraft_gpu::Renderer;
use fieldcraft_io::roots::Roots;
use fieldcraft_render::{Orientation, StandardView};
use serde_json::{Value, json};

/// A framed picture as PNG: the viewport's GPU pipeline on the process's shared device, or the
/// CPU rasterizer (`FIELDCRAFT_GPU`).
#[cfg(not(target_arch = "wasm32"))]
fn draw(f: &fieldcraft_engine::view::Framed, w: u32, h: u32) -> Result<(Vec<u8>, Renderer), String> {
    use fieldcraft_gpu::capture;
    let policy = fieldcraft_gpu::select::Policy::from_env().map_err(|e| e.to_string())?;
    capture::png(capture::shared(), policy, &f.scene, &[&f.symbols], &f.cam, w, h).map_err(|e| e.to_string())
}

/// A framed picture as PNG, drawn by the CPU rasterizer (the browser has no headless device).
#[cfg(target_arch = "wasm32")]
fn draw(f: &fieldcraft_engine::view::Framed, w: u32, h: u32) -> Result<(Vec<u8>, Renderer), String> {
    let png = fieldcraft_render::render_layers_png(&f.scene, &[&f.symbols], &f.cam, w as usize, h as usize).map_err(|e| format!("render: {e}"))?;
    Ok((png, Renderer::Cpu { why: None }))
}

/// Register the engine's picture drawer ([`fieldcraft_engine::view::register_drawer`]): the
/// viewport's GPU pipeline on the process's shared device under `FIELDCRAFT_GPU`, as MCP pictures
/// are drawn, for the pictures commands draw themselves (`results.export.animation`).
#[cfg(not(target_arch = "wasm32"))]
pub fn register_drawer() {
    fieldcraft_engine::view::register_drawer(Box::new(|f, w, h| fieldcraft_gpu::capture::rgba_shared(&f.scene, &[&f.symbols], &f.cam, w, h)));
}

/// What drew a picture and what it shows, as JSON: `renderer` (`kind` hardware, software or cpu;
/// the adapter, or why the CPU drew it), `results` (step, scale, contour), `undrawn_bodies` and
/// `sampled_symbols`. The MCP `screenshot` reply and `fieldcraft-cli snapshot` print it.
pub fn picture_report(r: &Renderer, shown: &Shown) -> Value {
    let renderer = r.report();
    let contour = shown
        .contour
        .as_ref()
        .map(|(field, min, max)| json!({"field": field, "min": {"node": min.0, "value": min.1}, "max": {"node": max.0, "value": max.1}}));
    json!({
        "renderer": renderer,
        "results": shown.step.as_ref().map(|step| json!({"step": step, "scale": shown.scale, "contour": contour})),
        "undrawn_bodies": shown.failed,
        "sampled_symbols": shown.sampled,
    })
}

/// What a picture shows.
#[derive(Debug, Clone, PartialEq)]
pub struct PictureRequest {
    /// A standard view to turn to first (`iso`, `front`…; `home` and `fit` frame the model).
    pub view: Option<String>,
    /// A camera to aim first ([`fieldcraft_engine::view::camera`]; `view` and `camera` are not
    /// given together). Without a target it looks at the model's centre of mass.
    pub camera: Option<Aim>,
    /// A section through the mesh, as `fieldcraft-cli snapshot --section` (headless only).
    pub section: Option<Section>,
    /// The whole app window (`true`) or the model only.
    pub window: bool,
    /// Image size in pixels (a window capture: points); `None` keeps the default.
    pub width: Option<u32>,
    pub height: Option<u32>,
}

/// A PNG and what it shows.
#[derive(Debug, Clone, PartialEq)]
pub struct Picture {
    pub png: Vec<u8>,
    pub info: Value,
}

/// Something that runs FieldCraft commands and draws pictures.
pub trait Backend {
    /// Run a command (`engine.execute`).
    fn execute(&mut self, command: &str, params: Value) -> Result<Value, String>;
    /// Every command with its parameter docs (`engine.commands`).
    fn commands(&mut self) -> Result<Value, String>;
    /// A PNG of the model, or of the app's window when one is attached.
    fn picture(&mut self, request: &PictureRequest) -> Result<Picture, String>;
    /// A request to the app's window (`ui.select`, `ui.selection`): an error without one.
    fn window(&mut self, method: &str, params: Value) -> Result<Value, String>;
    /// The folders file access is confined to (checked by the server for the files it writes
    /// itself).
    fn roots(&self) -> &Roots;
    /// For the client's instructions.
    fn describe(&self) -> String;
}

/// Default size of a model picture.
pub(crate) const DEFAULT_SIZE: (u32, u32) = (960, 600);
/// Largest picture side the server asks for, in pixels.
pub(crate) const MAX_SIDE: u32 = 4096;
/// Largest image read back from the app.
const MAX_IMAGE_BYTES: u64 = 64 << 20;
/// Longest reply line accepted from the app.
const MAX_REPLY: u64 = 256 << 20;

fn standard_view(v: &str) -> Result<StandardView, String> {
    match v {
        "home" | "fit" => Ok(StandardView::Iso),
        _ => StandardView::parse(v).ok_or_else(|| format!("unknown view `{}`", v.chars().take(40).collect::<String>())),
    }
}

/// An in-process headless session; pictures come from the CPU rasterizer.
pub struct Headless {
    pub session: Session,
    /// The camera of pictures: a standard view frames everything drawn (unless it has a
    /// distance); a camera given aims at its target ([`fieldcraft_engine::view::frame_at`]).
    pub view: Aim,
    roots: Roots,
    /// The session's camera requests seen (a saved view shown aims the next picture).
    camera_seen: u64,
}

impl Default for Headless {
    fn default() -> Self {
        Headless::new(Roots::default())
    }
}

impl Headless {
    /// An empty session whose commands read and write files only as `roots` allow.
    pub fn new(roots: Roots) -> Self {
        fieldcraft_io::recovery::use_automation_dir();
        let mut session = Session::default();
        session.solver_diagnostics = true;
        session.file_access = roots.clone();
        Headless { session, view: Aim::new(Orientation::of(StandardView::Iso)), roots, camera_seen: 0 }
    }
}

impl Backend for Headless {
    fn execute(&mut self, command: &str, params: Value) -> Result<Value, String> {
        self.session.execute(command, &params).map_err(|e| e.to_string())
    }

    fn commands(&mut self) -> Result<Value, String> {
        Ok(fieldcraft_engine::commands().iter().map(fieldcraft_engine::CommandSpec::listing).collect())
    }

    fn picture(&mut self, r: &PictureRequest) -> Result<Picture, String> {
        // A saved view shown since the last picture aims this one (`view.restore`, a model opened).
        if self.session.camera_requests() != self.camera_seen {
            self.camera_seen = self.session.camera_requests();
            if let Some(a) = &self.session.camera {
                self.view = a.targeted();
            }
        }
        if let Some(v) = &r.view {
            self.view = Aim::new(Orientation::of(standard_view(v)?));
        }
        if let Some(a) = &r.camera {
            self.view = Aim { target: Some(a.target.clone().unwrap_or(Target::MODEL)), ..a.clone() };
        }
        let w = r.width.unwrap_or(DEFAULT_SIZE.0).clamp(16, MAX_SIDE);
        let h = r.height.unwrap_or(DEFAULT_SIZE.1).clamp(16, MAX_SIDE);
        let f = if self.view.target.is_some() {
            fieldcraft_engine::view::frame_at(&self.session, &self.view, w as usize, h as usize, &Palette::LIGHT, r.section)
                .map_err(|e| e.to_string())?
        } else {
            fieldcraft_engine::view::frame(&self.session, &self.view.at, w as usize, h as usize, &Palette::LIGHT, r.section)
        };
        let (png, renderer) = draw(&f, w, h)?;
        // The camera drawn with is the session's, what `view.save` keeps.
        self.session.camera = Some(fieldcraft_engine::view::aim_of(&f.cam, &self.session));
        let mut info = picture_report(&renderer, &f.shown);
        if let Some(o) = info.as_object_mut() {
            o.insert("source".into(), "model".into());
            o.insert("view".into(), self.view.json());
            o.insert("camera".into(), fieldcraft_engine::view::camera_shown(&f.cam, &self.session));
            o.insert("width".into(), w.into());
            o.insert("height".into(), h.into());
        }
        Ok(Picture { png, info })
    }

    fn window(&mut self, method: &str, _params: Value) -> Result<Value, String> {
        Err(format!("{method} needs the app's window: start the server with --connect to a running FieldCraft"))
    }

    fn roots(&self) -> &Roots {
        &self.roots
    }

    fn describe(&self) -> String {
        "a headless in-process FieldCraft session (no window; screenshots drawn offscreen with the viewport's GPU pipeline, or on the CPU without a GPU)".into()
    }
}

/// A running FieldCraft app reached through its loopback control channel
/// (`fieldcraft --control PORT`).
pub struct Remote {
    addr: String,
    conn: Option<(BufReader<TcpStream>, TcpStream)>,
    next_id: u64,
    roots: Roots,
}

static SHOT: AtomicU64 = AtomicU64::new(0);

impl Remote {
    pub fn connect(addr: &str) -> std::io::Result<Self> {
        let mut r = Self { addr: addr.to_string(), conn: None, next_id: 1, roots: Roots::default() };
        r.reconnect()?;
        Ok(r)
    }

    fn reconnect(&mut self) -> std::io::Result<()> {
        self.conn = None;
        let mut last = std::io::Error::new(std::io::ErrorKind::NotFound, format!("cannot resolve {}", self.addr));
        for sa in self.addr.to_socket_addrs()? {
            match TcpStream::connect_timeout(&sa, Duration::from_millis(800)) {
                Ok(s) => {
                    s.set_nodelay(true).ok();
                    // Meshing and solving big models take as long as they take: the app answers
                    // a solve or a mesher when its job ends (docs/control-protocol.md, Jobs), and
                    // a closed app closes the connection.
                    s.set_read_timeout(None).ok();
                    let read = s.try_clone()?;
                    self.conn = Some((BufReader::new(read), s));
                    return Ok(());
                }
                Err(e) => last = e,
            }
        }
        Err(last)
    }

    fn roundtrip(&mut self, line: &str) -> std::io::Result<String> {
        if self.conn.is_none() {
            self.reconnect()?;
        }
        let Some((reader, writer)) = self.conn.as_mut() else { return Err(std::io::Error::other("not connected")) };
        let sent = writer.write_all(line.as_bytes()).and_then(|()| writer.write_all(b"\n")).and_then(|()| writer.flush());
        if let Err(e) = sent {
            self.conn = None;
            return Err(e);
        }
        let mut reply = String::new();
        let n = match reader.by_ref().take(MAX_REPLY).read_line(&mut reply) {
            Ok(n) => n,
            Err(e) => {
                // A timed-out read leaves a late reply on the stream: start afresh next time.
                self.conn = None;
                return Err(e);
            }
        };
        if n == 0 || !reply.ends_with('\n') {
            self.conn = None;
            return Err(std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "connection closed (or reply too large)"));
        }
        Ok(reply)
    }

    /// One control-channel request.
    pub fn call(&mut self, method: &str, params: Value) -> Result<Value, String> {
        let id = self.next_id;
        self.next_id += 1;
        let line = json!({"id": id, "method": method, "params": params}).to_string();
        let reply = match self.roundtrip(&line) {
            Ok(r) => r,
            // One retry on a fresh connection (the app may have restarted).
            Err(_) => self.roundtrip(&line).map_err(|e| format!("FieldCraft at {}: {e}", self.addr))?,
        };
        let v: Value = serde_json::from_str(&reply).map_err(|e| format!("bad reply from the app: {e}"))?;
        if v.get("ok").and_then(Value::as_bool) == Some(true) {
            Ok(v.get("result").cloned().unwrap_or(Value::Null))
        } else {
            Err(v.get("error").and_then(Value::as_str).unwrap_or("error").to_string())
        }
    }
}

/// A temporary PNG path the app writes and the server reads back (both on this computer).
fn temp_png() -> PathBuf {
    let n = SHOT.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir().join(format!("fieldcraft-mcp-{}-{n}.png", std::process::id()))
}

impl Backend for Remote {
    fn execute(&mut self, command: &str, params: Value) -> Result<Value, String> {
        // `client` lets the app count an agent's commands apart from a script's (docs/privacy.md).
        self.call("engine.execute", json!({"command": command, "params": params, "client": "mcp"}))
    }

    fn commands(&mut self) -> Result<Value, String> {
        self.call("engine.commands", json!({}))
    }

    fn picture(&mut self, r: &PictureRequest) -> Result<Picture, String> {
        if let Some(v) = &r.view {
            standard_view(v)?;
            self.call("ui.view", json!({"view": v}))?;
        }
        if r.section.is_some() {
            return Err("screenshot: `section` draws headless only (the app's window has no section view)".into());
        }
        if let Some(a) = &r.camera {
            // Aimed as a headless picture is: at the model's centre of mass unless told otherwise.
            let a = Aim { target: Some(a.target.clone().unwrap_or(Target::MODEL)), ..a.clone() };
            self.call("ui.camera", a.json())?;
        }
        let path = temp_png();
        let shown = path.to_string_lossy().into_owned();
        let mut params = json!({"path": shown});
        let method = if r.window {
            for (k, v) in [("width", r.width), ("height", r.height)] {
                if let Some(v) = v {
                    params[k] = json!(v.clamp(16, MAX_SIDE));
                }
            }
            "ui.capture"
        } else {
            params["width"] = json!(r.width.unwrap_or(DEFAULT_SIZE.0).clamp(16, MAX_SIDE));
            params["height"] = json!(r.height.unwrap_or(DEFAULT_SIZE.1).clamp(16, MAX_SIDE));
            "ui.render"
        };
        let result = self.call(method, params);
        let read = result.and_then(|info| {
            let len = std::fs::metadata(&path).map_err(|e| format!("{shown}: {e} (is the app on another computer?)"))?.len();
            if len > MAX_IMAGE_BYTES {
                return Err(format!("{shown}: image too large ({len} bytes)"));
            }
            let png = std::fs::read(&path).map_err(|e| format!("{shown}: {e}"))?;
            Ok(Picture { png, info: json!({"source": if r.window { "window" } else { "model" }, "app": info}) })
        });
        let _ = std::fs::remove_file(&path);
        read
    }

    fn window(&mut self, method: &str, params: Value) -> Result<Value, String> {
        self.call(method, params)
    }

    fn roots(&self) -> &Roots {
        &self.roots
    }

    fn describe(&self) -> String {
        format!("connected to the FieldCraft app at {} (the user watches its window)", self.addr)
    }
}
