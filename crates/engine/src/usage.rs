//! Anonymous usage analytics, off until the user agrees (docs/privacy.md).
//!
//! The app records [`Event`]s here; a host (the desktop app's sender thread, the web page's
//! `fetch`) takes them in batches ([`take_batch`]) and posts [`Batch::body`] to PostHog's batch
//! API. Nothing here touches the network.
//!
//! What an event can hold is fixed by its type: command ids, sample names and requirements come
//! from the registries (`&'static str`), never from parameters, paths, file names or messages.
//! Recording does nothing until [`enable`] (the user agreed, and the host may send); [`disable`]
//! drops what is queued. The queue is bounded ([`MAX_QUEUED`]): when full, new events are dropped.
//!
//! [`Session::execute`](crate::Session::execute) counts every top-level command of a session
//! whose [`Session::origin`](crate::Session::origin) is set (the app's own model): its id, kind,
//! duration, whether it failed and how, and where it came from ([`Origin`]); the funnel's
//! [`Milestone`]s and [`Event::SampleOpened`] follow from the same commands.

use std::cell::Cell;
use std::collections::VecDeque;
use std::sync::Mutex;

use serde_json::{Map, Value, json};

use crate::requirement::Requirement;
use crate::{CommandKind, EngineError};

/// PostHog Cloud's EU ingestion host.
pub const DEFAULT_HOST: &str = "https://eu.i.posthog.com";
/// Events kept waiting to be sent; more are dropped.
pub const MAX_QUEUED: usize = 1000;
/// Events in one request.
pub const MAX_BATCH: usize = 100;

/// Where a command came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Origin {
    /// A click: the ribbon, a menu, a form, the browser.
    Ui,
    /// The command search or the command line.
    Palette,
    /// A script or the control channel.
    Script,
    /// An AI agent through the MCP server.
    Mcp,
}

impl Origin {
    pub fn as_str(self) -> &'static str {
        match self {
            Origin::Ui => "ui",
            Origin::Palette => "palette",
            Origin::Script => "script",
            Origin::Mcp => "mcp",
        }
    }

    /// Scripts and agents, which the user can leave out on their own.
    pub fn automated(self) -> bool {
        matches!(self, Origin::Script | Origin::Mcp)
    }
}

/// A step of the way from a design to results, counted once per run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Milestone {
    DesignOpened,
    Meshed,
    Solved,
    ResultsViewed,
}

impl Milestone {
    pub fn as_str(self) -> &'static str {
        match self {
            Milestone::DesignOpened => "design_opened",
            Milestone::Meshed => "meshed",
            Milestone::Solved => "solved",
            Milestone::ResultsViewed => "results_viewed",
        }
    }

    fn bit(self) -> u8 {
        1 << self as u8
    }
}

/// Something that happened, with only what docs/privacy.md lists.
#[derive(Clone, Debug, PartialEq)]
pub enum Event {
    AppStarted,
    /// The window closed after `seconds`.
    AppClosed {
        seconds: u64,
    },
    CommandRun {
        command: &'static str,
        kind: CommandKind,
        duration_ms: u64,
        error: Option<&'static str>,
        from: Origin,
    },
    Milestone(Milestone),
    SampleOpened {
        name: &'static str,
    },
    /// A command's form closed without running it.
    FormCancelled {
        command: &'static str,
    },
    /// A command clicked before the model had what it needs.
    UnavailableClicked {
        command: &'static str,
        missing: Vec<Requirement>,
    },
}

impl Event {
    pub fn name(&self) -> &'static str {
        match self {
            Event::AppStarted => "app_started",
            Event::AppClosed { .. } => "app_closed",
            Event::CommandRun { .. } => "command_run",
            Event::Milestone(_) => "milestone",
            Event::SampleOpened { .. } => "sample_opened",
            Event::FormCancelled { .. } => "form_cancelled",
            Event::UnavailableClicked { .. } => "unavailable_clicked",
        }
    }

    /// The event's own properties.
    pub fn properties(&self) -> Map<String, Value> {
        let v = match self {
            Event::AppStarted => json!({}),
            Event::AppClosed { seconds } => json!({"duration_s": seconds}),
            Event::CommandRun { command, kind, duration_ms, error, from } => {
                let kind = match kind {
                    CommandKind::Model => "model",
                    CommandKind::Session => "session",
                };
                let mut v = json!({"command": command, "kind": kind, "duration_ms": duration_ms, "ok": error.is_none(), "from": from.as_str()});
                if let Some(e) = error {
                    v["error"] = json!(e);
                }
                v
            }
            Event::Milestone(m) => json!({"name": m.as_str()}),
            Event::SampleOpened { name } => json!({"name": name}),
            Event::FormCancelled { command } => json!({"command": command}),
            Event::UnavailableClicked { command, missing } => json!({"command": command, "missing": missing}),
        };
        match v {
            Value::Object(m) => m,
            _ => Map::new(),
        }
    }
}

/// The class of an error, without its message (which can name the user's files and values).
pub fn error_class(e: &EngineError) -> &'static str {
    match e {
        EngineError::UnknownCommand(_) => "unknown_command",
        EngineError::Params(_) | EngineError::Param(_) => "params",
        EngineError::Model(_) => "model",
        EngineError::Geometry(_) => "geometry",
        EngineError::Mesh(_) => "mesh",
        EngineError::File(_) => "file",
        EngineError::Deck(_) => "deck",
        EngineError::Solver(_) => "solver",
        EngineError::Explicit(_) => "explicit",
        EngineError::Fluid(_) => "fluid",
        EngineError::Unavailable(_) => "unavailable",
        EngineError::Internal { .. } => "internal",
    }
}

/// Who is sending: the anonymous install and the build. The same for every event of a run.
#[derive(Clone, Debug, PartialEq)]
pub struct Identity {
    /// A random id made for this install (not derived from the machine or the user).
    pub install_id: String,
    /// `desktop` or `web`.
    pub platform: &'static str,
    /// The operating system (`macos`, `linux`, `windows`; `web` in a browser).
    pub os: &'static str,
    pub arch: &'static str,
    pub app_version: &'static str,
}

/// What [`take_batch`] hands the sender.
#[derive(Clone, Debug, PartialEq)]
pub struct Batch {
    pub identity: Identity,
    pub events: Vec<Event>,
}

impl Batch {
    /// The body of a POST to PostHog's `/batch/` endpoint with the project's API key: each event
    /// with the anonymous install id, the build, no person profile and no IP-based location.
    pub fn body(&self, api_key: &str) -> String {
        let id = &self.identity;
        let batch: Vec<Value> = self
            .events
            .iter()
            .map(|e| {
                let mut p = e.properties();
                p.insert("app_version".into(), json!(id.app_version));
                p.insert("platform".into(), json!(id.platform));
                p.insert("os".into(), json!(id.os));
                p.insert("arch".into(), json!(id.arch));
                p.insert("$process_person_profile".into(), json!(false));
                p.insert("$geoip_disable".into(), json!(true));
                json!({"event": e.name(), "distinct_id": id.install_id, "properties": p})
            })
            .collect();
        json!({"api_key": api_key, "batch": batch}).to_string()
    }
}

/// Whether this run may record at all, before asking the user: the build has a project key,
/// `FIELDCRAFT_TELEMETRY` (`telemetry`) is not `0`, and a run driven by automation (`--control`)
/// or CI records only when `FIELDCRAFT_TELEMETRY=1` says so.
pub fn allowed(has_key: bool, telemetry: Option<&str>, automation: bool, ci: bool) -> bool {
    let set = telemetry.map(|v| v.trim().to_ascii_lowercase());
    match set.as_deref() {
        Some("0" | "false" | "off" | "no") => false,
        Some("1" | "true" | "on" | "yes") => has_key,
        _ => has_key && !automation && !ci,
    }
}

struct State {
    identity: Identity,
    /// Count scripts and MCP runs too.
    automation: bool,
    queue: VecDeque<Event>,
    /// Milestones already sent this run.
    reached: u8,
    /// When recording started ([`now_ms`]), for the session's length.
    started: Option<f64>,
}

static STATE: Mutex<Option<State>> = Mutex::new(None);

fn with_state<R>(f: impl FnOnce(&mut Option<State>) -> R) -> R {
    // A poisoned lock still holds a usable queue.
    let mut g = STATE.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    f(&mut g)
}

/// Start recording, as `identity`, with [`Event::AppStarted`]; `automation` counts scripts and
/// MCP runs too. Already recording, only the settings change.
pub fn enable(identity: Identity, automation: bool) {
    let now = now_ms();
    with_state(|s| match s {
        Some(s) => {
            s.identity = identity;
            s.automation = automation;
        }
        None => *s = Some(State { identity, automation, queue: VecDeque::from([Event::AppStarted]), reached: 0, started: now }),
    });
}

/// The app is closing: [`Event::AppClosed`] with how long it recorded.
pub fn closing() {
    let now = now_ms();
    let started = with_state(|s| s.as_ref().map(|s| s.started));
    if let Some(started) = started {
        let seconds = match (started, now) {
            (Some(a), Some(b)) if b >= a => ((b - a) / MS_PER_SECOND).round() as u64,
            _ => 0,
        };
        record(Event::AppClosed { seconds });
    }
}

/// Stop recording and drop everything queued.
pub fn disable() {
    with_state(|s| *s = None);
}

pub fn enabled() -> bool {
    with_state(|s| s.is_some())
}

/// Queue an event, when recording (and when it is no script or MCP run the user left out).
pub fn record(event: Event) {
    with_state(|s| {
        let Some(s) = s else { return };
        match &event {
            Event::CommandRun { from, .. } if from.automated() && !s.automation => return,
            Event::Milestone(m) if s.reached & m.bit() != 0 => return,
            Event::Milestone(m) => s.reached |= m.bit(),
            _ => {}
        }
        if s.queue.len() < MAX_QUEUED {
            s.queue.push_back(event);
        }
    });
}

/// Up to [`MAX_BATCH`] queued events to send, oldest first; `None` when there are none.
pub fn take_batch() -> Option<Batch> {
    with_state(|s| {
        let s = s.as_mut()?;
        if s.queue.is_empty() {
            return None;
        }
        let n = s.queue.len().min(MAX_BATCH);
        Some(Batch { identity: s.identity.clone(), events: s.queue.drain(..n).collect() })
    })
}

/// Events waiting to be sent.
pub fn queued() -> usize {
    with_state(|s| s.as_ref().map_or(0, |s| s.queue.len()))
}

/// Milliseconds in a second, for the session's length and the native clock.
const MS_PER_SECOND: f64 = 1000.0; // unit-factors: wall-clock time of the usage counts, not a model quantity

/// A clock in milliseconds for command durations. The plain wasm target has none: the web page
/// sets one ([`set_clock`]).
static CLOCK: Mutex<Option<fn() -> f64>> = Mutex::new(None);

/// Use `clock` (milliseconds, any origin) for command durations.
pub fn set_clock(clock: fn() -> f64) {
    *CLOCK.lock().unwrap_or_else(std::sync::PoisonError::into_inner) = Some(clock);
}

fn now_ms() -> Option<f64> {
    if let Some(c) = *CLOCK.lock().unwrap_or_else(std::sync::PoisonError::into_inner) {
        return Some(c());
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        static START: std::sync::OnceLock<std::time::Instant> = std::sync::OnceLock::new();
        Some(START.get_or_init(std::time::Instant::now).elapsed().as_secs_f64() * MS_PER_SECOND)
    }
    #[cfg(target_arch = "wasm32")]
    None
}

thread_local! {
    /// Commands running on this thread: only the outermost is counted.
    static DEPTH: Cell<u32> = const { Cell::new(0) };
}

/// A command being counted, from [`Session::execute`](crate::Session::execute): `None` when it
/// is not (not recording, no origin, or a command inside another).
pub(crate) struct Timer {
    started: Option<f64>,
    from: Option<Origin>,
}

impl Timer {
    pub(crate) fn start(from: Option<Origin>) -> Timer {
        let outer = DEPTH.with(|d| {
            d.set(d.get() + 1);
            d.get() == 1
        });
        let from = from.filter(|_| outer && enabled());
        Timer { started: from.and_then(|_| now_ms()), from }
    }

    /// Whether this command is counted.
    pub(crate) fn counting(&self) -> bool {
        self.from.is_some()
    }

    /// The command ended: count it, and the milestone it reached.
    pub(crate) fn done(
        &self,
        command: &'static str,
        kind: CommandKind,
        error: Option<&EngineError>,
        milestone: Option<Milestone>,
        sample: Option<&'static str>,
    ) {
        let Some(from) = self.from else { return };
        let duration_ms = match (self.started, now_ms()) {
            (Some(a), Some(b)) if b >= a => (b - a).round() as u64,
            _ => 0,
        };
        record(Event::CommandRun { command, kind, duration_ms, error: error.map(error_class), from });
        if error.is_none() {
            if let Some(m) = milestone {
                record(Event::Milestone(m));
            }
            if let Some(name) = sample {
                record(Event::SampleOpened { name });
            }
        }
    }
}

impl Drop for Timer {
    fn drop(&mut self) {
        DEPTH.with(|d| d.set(d.get().saturating_sub(1)));
    }
}

/// The milestone a successful `command` reaches; `elements` is the element count before and
/// after it.
pub(crate) fn milestone(command: &str, elements: (usize, usize)) -> Option<Milestone> {
    match command {
        "geometry.open" | "geometry.import.mesh" | "file.open" => Some(Milestone::DesignOpened),
        c if c.starts_with("mesh.") && elements.1 > elements.0 => Some(Milestone::Meshed),
        // Solvers, not the queries on a solved step.
        c if c.starts_with("solve.") && !matches!(c, "solve.probe" | "solve.history") => Some(Milestone::Solved),
        c if c.starts_with("results.") && !matches!(c, "results.reactions" | "results.clear") => Some(Milestone::ResultsViewed),
        _ => None,
    }
}

/// Serialises the tests (and other code) that turn the process-wide recording on.
#[doc(hidden)]
pub static TEST_LOCK: Mutex<()> = Mutex::new(());

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Session;

    fn identity() -> Identity {
        Identity { install_id: "00000000-0000-4000-8000-000000000000".into(), platform: "desktop", os: "test", arch: "test", app_version: "0.0.0" }
    }

    fn take_all() -> Vec<Event> {
        std::iter::from_fn(take_batch).flat_map(|b| b.events).collect()
    }

    /// Every string in a payload is a key, a registry name or the identity, and every number a
    /// duration: no path, file name, parameter value or message can be in it.
    pub(crate) fn assert_clean(body: &str) {
        let v: Value = serde_json::from_str(body).unwrap();
        let mut allowed: Vec<String> = crate::commands().iter().map(|c| c.id.to_string()).collect();
        allowed.extend(crate::sample::SAMPLES.iter().map(|s| s.name.to_string()));
        for e in [
            "app_started",
            "app_closed",
            "command_run",
            "milestone",
            "sample_opened",
            "form_cancelled",
            "unavailable_clicked",
            "model",
            "session",
            "ui",
            "palette",
            "script",
            "mcp",
            "design_opened",
            "meshed",
            "solved",
            "results_viewed",
            "unknown_command",
            "params",
            "geometry",
            "mesh",
            "file",
            "deck",
            "solver",
            "explicit",
            "fluid",
            "unavailable",
            "internal",
            "KEY",
        ] {
            allowed.push(e.into());
        }
        let id = identity();
        allowed.extend([id.install_id, id.platform.into(), id.os.into(), id.arch.into(), id.app_version.into()]);
        fn walk(key: &str, v: &Value, allowed: &[String]) {
            match v {
                Value::String(s) => {
                    // A requirement's snake_case name.
                    let requirement = key == "missing" && s.chars().all(|c| c.is_ascii_lowercase() || c == '_');
                    assert!(allowed.contains(s) || requirement, "`{key}`: `{s}` is not a name FieldCraft defines");
                }
                Value::Number(_) => assert!(["duration_ms", "duration_s"].contains(&key), "`{key}`: a number that is not a duration"),
                Value::Array(a) => a.iter().for_each(|x| walk(key, x, allowed)),
                Value::Object(m) => m.iter().for_each(|(k, x)| walk(k, x, allowed)),
                Value::Bool(_) | Value::Null => {}
            }
        }
        walk("", &v, &allowed);
    }

    #[test]
    fn nothing_is_queued_without_consent() {
        let _g = TEST_LOCK.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        disable();
        let mut s = Session { origin: Some(Origin::Ui), ..Default::default() };
        s.run_script(&crate::sample::cantilever()).unwrap();
        record(Event::AppStarted);
        assert_eq!(queued(), 0);
        assert!(take_batch().is_none());
    }

    #[test]
    fn a_scripted_session_sends_only_names_and_durations() {
        let _g = TEST_LOCK.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        disable();
        enable(identity(), true);
        let dir = std::env::temp_dir().join(format!("fieldcraft-usage-secret-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let design = dir.join("my-secret-bracket.solidcraft");
        std::fs::write(&design, fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap()).unwrap();
        let mut s = Session { origin: Some(Origin::Script), ..Default::default() };
        s.execute("geometry.open", &json!({"path": design.to_string_lossy()})).unwrap();
        s.run_script(&crate::sample::cantilever()).unwrap();
        s.execute("model.material.create", &json!({"name": "Secret alloy 7", "elastic": {"young": 123456.75, "poisson": 0.3125}})).ok();
        s.execute("solve.run", &json!({})).unwrap();
        s.execute("results.contour", &json!({"field": "von_mises"})).unwrap();
        // A failure whose message names a path, and a sample.
        assert!(s.execute("file.open", &json!({"path": dir.join("missing-secret.fieldcraft").to_string_lossy()})).is_err());
        s.execute("file.sample", &json!({"name": "cantilever"})).unwrap();
        // An unknown id is user text: not counted.
        assert!(s.execute("secret.command", &json!({})).is_err());
        let events = take_all();
        let body = Batch { identity: identity(), events: events.clone() }.body("KEY");
        assert_clean(&body);
        for secret in ["secret", "Secret", "123456", "0.3125", "1000", dir.to_string_lossy().as_ref()] {
            assert!(!body.contains(secret), "{secret} in {body}");
        }
        let names: Vec<&str> = events.iter().map(Event::name).collect();
        assert!(names.contains(&"command_run") && names.contains(&"sample_opened"), "{names:?}");
        // Inner commands (the script's steps run by execute, a sample's own script) count once.
        let runs = events.iter().filter(|e| matches!(e, Event::CommandRun { command: "file.sample", .. })).count();
        assert_eq!(runs, 1);
        let reached: Vec<_> = events.iter().filter_map(|e| if let Event::Milestone(m) = e { Some(*m) } else { None }).collect();
        for m in [Milestone::DesignOpened, Milestone::Meshed, Milestone::Solved, Milestone::ResultsViewed] {
            assert_eq!(reached.iter().filter(|r| **r == m).count(), 1, "{m:?} once: {reached:?}");
        }
        assert!(
            events.iter().any(|e| matches!(e, Event::CommandRun { command: "file.open", error: Some("file"), from: Origin::Script, .. })),
            "{events:?}"
        );
        disable();
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn scripts_and_agents_can_be_left_out_and_the_queue_is_bounded() {
        let _g = TEST_LOCK.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        disable();
        enable(identity(), false);
        let run = |from| Event::CommandRun { command: "model.info", kind: CommandKind::Model, duration_ms: 1, error: None, from };
        record(run(Origin::Mcp));
        record(run(Origin::Script));
        record(run(Origin::Palette));
        assert_eq!(take_all(), [Event::AppStarted, run(Origin::Palette)]);
        for _ in 0..MAX_QUEUED + 50 {
            record(Event::AppStarted);
        }
        assert_eq!(queued(), MAX_QUEUED);
        let first = take_batch().unwrap();
        assert_eq!(first.events.len(), MAX_BATCH);
        assert_eq!(queued(), MAX_QUEUED - MAX_BATCH);
        // Turning it off drops what waits.
        disable();
        assert_eq!(queued(), 0);
        enable(identity(), true);
        closing();
        assert!(matches!(take_all()[..], [Event::AppStarted, Event::AppClosed { .. }]));
        disable();
    }

    #[test]
    fn the_body_is_a_posthog_batch() {
        let b = Batch {
            identity: identity(),
            events: vec![Event::UnavailableClicked { command: "solve.run", missing: vec![Requirement::Step, Requirement::Material] }],
        };
        let v: Value = serde_json::from_str(&b.body("phc_test")).unwrap();
        assert_eq!(v["api_key"], "phc_test");
        let e = &v["batch"][0];
        assert_eq!((e["event"].as_str(), e["distinct_id"].as_str()), (Some("unavailable_clicked"), Some(identity().install_id.as_str())));
        assert_eq!(e["properties"]["missing"], json!(["step", "material"]));
        assert_eq!(e["properties"]["$process_person_profile"], false);
        assert_eq!(e["properties"]["$geoip_disable"], true);
        assert_clean(&b.body("KEY"));
    }

    #[test]
    fn telemetry_is_off_without_a_key_in_automation_and_ci_and_when_turned_off() {
        assert!(allowed(true, None, false, false));
        assert!(!allowed(false, None, false, false));
        assert!(!allowed(false, Some("1"), false, false));
        assert!(!allowed(true, Some("0"), false, false));
        assert!(!allowed(true, Some(" off "), false, false));
        assert!(!allowed(true, None, true, false));
        assert!(!allowed(true, None, false, true));
        assert!(!allowed(true, Some(""), true, false));
        assert!(allowed(true, Some("1"), true, true));
    }
}
