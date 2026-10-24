//! Crash reports: enough to reproduce a crash and fix it, without slowing the app down or
//! filling the disk (docs/crash-reports.md).
//!
//! - **Panics** are caught by a hook. One inside a command becomes a report with the model as it
//!   was before the command and the command itself, so `fieldcraft-cli crashes replay <id>`
//!   reproduces it. A panic something inside the command caught and turned into an error is
//!   reported the same way, marked recovered.
//! - **Crashes no hook sees** (an abort, a stack overflow, a segfault, a crash inside a graphics
//!   driver, the process being killed) are found at the next start: every running app keeps a
//!   small, locked command journal; a journal nobody holds any more belongs to an app that died.
//!   Its report carries the model the journal starts from (the last autosave, the model's file,
//!   or a copy taken when the model was opened) and the commands run since, ending with the one
//!   in flight, plus the app's last control request and the dead run's lines of the log file
//!   ([`crate::logfile`]).
//! - **Solver errors**: an opted-in host writes a small local diagnostic at a failed solve's
//!   boundary, without copying the model or presenting a startup crash prompt.
//! - **Cost**: the crash journal does nothing until [`install`] is called; after that, two
//!   small unbuffered appends per top-level command of the app's own model
//!   ([`Session::crash_journal`]). Reports are assembled only after a panic, recovered crash,
//!   or an opted-in solver error, never during solver iterations.
//! - **Bounded**: the same crash (its signature: place + message, digits ignored) is counted, not
//!   written again; at most [`MAX_NEW_PER_SESSION`] new reports per run, [`MAX_REPORTS`]
//!   reports, [`MAX_TOTAL_BYTES`] on disk, nothing older than [`MAX_AGE_DAYS`] days.
//!
//! - **Startup steps** ([`phase`]): the journal keeps the step the app is in while it gets its
//!   window and graphics device up (window, graphics instance, adapter, device, shader
//!   compiles, first frame, [`RUNNING`]), so a crash inside a graphics driver still says where
//!   it happened, and the renderer ([`set_environment`]) once the adapter is picked.
//! - **Native crashes** (Windows): the app's exception filter writes what it caught (exception
//!   code, faulting module and offset) and a minidump next to the journal ([`native_paths`]);
//!   the next start puts both into the unclean exit's report.
//!
//! Reports are folders `crashes/<time>-<signature>/` holding `report.json`, and when known
//! `model.fieldcraft`, `replay.json`, `log.txt` and `minidump.dmp`. Nothing is ever sent anywhere.
//!
//! Adapted from SolidCraft's `crates/engine/src/crash.rs` (MIT OR Apache-2.0), see ATTRIBUTION.md.

use std::backtrace::Backtrace;
use std::cell::{Cell, RefCell};
use std::collections::HashMap;
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Mutex, OnceLock};

use serde_json::{Value, json};

use crate::{EngineError, Session};

/// Reports kept (the oldest, and recovered panics, go first).
pub const MAX_REPORTS: usize = 20;
/// Bytes all reports may take together.
pub const MAX_TOTAL_BYTES: u64 = 100 << 20;
/// Reports older than this are deleted.
pub const MAX_AGE_DAYS: u64 = 90;
/// New reports one run of the app may write (repeats are only counted).
pub const MAX_NEW_PER_SESSION: usize = 5;
/// A model bigger than this is not copied into a report.
pub const MAX_MODEL_BYTES: usize = 20 << 20;
/// A journal stops recording commands at this size (until the next autosave restarts it).
const MAX_JOURNAL_BYTES: u64 = 512 << 10;
/// Command parameters recorded per command.
const MAX_PARAMS_BYTES: usize = 16 << 10;
/// Recent commands shown in a report.
const RECENT: usize = 32;
const MAX_FRAMES: usize = 60;
const MAX_MESSAGE: usize = 2000;
/// The log lines copied into a report.
const LOG_TAIL_BYTES: usize = 64 << 10;
/// The report format.
pub const FORMAT: &str = "fieldcraft-crash/1";
/// The model file in a report folder.
pub const MODEL_FILE: &str = "model.fieldcraft";
/// The minidump of a native crash in a report folder (Windows).
pub const MINIDUMP_FILE: &str = "minidump.dmp";

/// What the app is (shown in every report).
#[derive(Clone, Debug, Default)]
pub struct Build {
    pub version: &'static str,
    /// The commit it was built from, when known.
    pub commit: &'static str,
}

impl Build {
    /// This build: the workspace version, and the commit when the build named one
    /// (`FIELDCRAFT_BUILD_SHA`, packaging/env.sh).
    pub fn current() -> Build {
        let sha = option_env!("FIELDCRAFT_BUILD_SHA").unwrap_or("");
        Build { version: env!("CARGO_PKG_VERSION"), commit: sha.get(..10).unwrap_or(sha) }
    }

    /// `v0.5.8 (abc1234567), macos aarch64`.
    pub fn version_line(&self) -> String {
        crate::bug_report::version_line(self.version, self.commit)
    }
}

/// One report, as listed.
#[derive(Clone, Debug, PartialEq, serde::Serialize)]
pub struct Summary {
    pub id: String,
    /// `command_panic`, `panic`, `recovered_panic`, `unclean_exit` or `solver_error`.
    pub kind: String,
    /// When it first happened (UTC).
    pub time: String,
    pub message: String,
    /// The command running when it happened, if one was.
    pub command: Option<String>,
    /// Times it happened.
    pub count: u64,
    pub version: String,
    /// The report's folder.
    pub path: String,
    /// `replay.json` is there: the crash can be replayed.
    pub replayable: bool,
}

struct State {
    dir: PathBuf,
    build: Build,
    journal: Mutex<Journal>,
    /// Signatures already reported (this run's count; the report's folder name).
    known: Mutex<HashMap<u64, (u64, String)>>,
    written: AtomicUsize,
    started: u64,
    /// When the app started, in milliseconds (for how far into the start a step came).
    started_ms: u64,
}

/// The startup step that ends startup: the app is up and drawing.
pub const RUNNING: &str = "running";

static STATE: OnceLock<State> = OnceLock::new();
/// The hook is writing a report (a panic inside it must not write another).
static IN_HOOK: AtomicBool = AtomicBool::new(false);
static SEQ: AtomicU64 = AtomicU64::new(0);

thread_local! {
    /// Commands running on this thread (their panics are reported when the outermost ends).
    static IN_COMMAND: Cell<u32> = const { Cell::new(0) };
    /// Panics caught inside the running command, waiting for it to end.
    static PENDING: RefCell<Vec<Panic>> = const { RefCell::new(Vec::new()) };
}

/// Crash reporting is on in this process.
pub fn enabled() -> bool {
    STATE.get().is_some()
}

/// The crash folder: `FIELDCRAFT_CRASH_DIR`, else `crashes` in the app's data folder (a
/// temporary folder in runs driven by automation, [`fieldcraft_io::recovery::data_dir`]).
pub fn default_dir() -> Option<PathBuf> {
    if let Some(s) = STATE.get() {
        return Some(s.dir.clone());
    }
    if let Some(d) = std::env::var_os("FIELDCRAFT_CRASH_DIR").filter(|d| !d.is_empty()) {
        return Some(PathBuf::from(d));
    }
    fieldcraft_io::recovery::data_dir().map(|d| d.join("crashes"))
}

/// A small local diagnostic for an ordinary solver failure, without a model, replay,
/// backtrace or startup crash announcement. Parameters and context are bounded independently;
/// repeated failures update one record. Failure to persist never changes the solver outcome.
pub(crate) fn solver_error(dir: &Path, command: &str, params: &Value, context: &Value, error: &str) -> Option<PathBuf> {
    static WRITING: Mutex<()> = Mutex::new(());
    static WRITTEN: AtomicUsize = AtomicUsize::new(0);
    let _writing = WRITING.lock().unwrap_or_else(|e| e.into_inner());
    let t = now();
    prune(dir, t);
    // A different step/request/backend is not a repeat of the same failed solve. Keep exact
    // parameters in the location key: crash signatures otherwise intentionally ignore digits.
    let key_params = bounded_json(params, MAX_PARAMS_BYTES).unwrap_or_else(|| "<oversized>".into());
    let key_backend = json!([
        context.pointer("/backend/actual"),
        context.pointer("/backend/observed_context/backend"),
        context.pointer("/backend/observed_context/adapter")
    ]);
    let sig = signature(&format!("solver_error:{command}:{key_params}:{key_backend}"), error);
    let previous = list(dir).into_iter().find(|s| s.kind == "solver_error" && sig_of(&s.id) == Some(sig));
    if previous.is_none() && WRITTEN.load(Ordering::Relaxed) >= MAX_NEW_PER_SESSION {
        return None;
    }
    let id = previous.as_ref().map(|s| s.id.clone()).unwrap_or_else(|| format!("{}-{sig:016x}", stamp(t)));
    let folder = dir.join(&id);
    let bounded = |v: &Value, limit| bounded_json(v, limit).and_then(|s| serde_json::from_str::<Value>(&s).ok());
    let parameters = bounded(params, MAX_PARAMS_BYTES);
    let solve = bounded(context, 64 << 10);
    let environment = STATE
        .get()
        .and_then(|s| s.journal.lock().ok().and_then(|j| j.environment.clone()))
        .and_then(|line| serde_json::from_str::<Value>(&line).ok())
        .and_then(|v| v.get("environment").cloned());
    let mut report = json!({
        "format": "fieldcraft-diagnostic/1",
        "id": id,
        "kind": "solver_error",
        "message": limited(error, MAX_MESSAGE),
        "command": {"command": limited(command, 128), "params": parameters, "params_omitted": parameters.is_none().then_some("parameters exceed the 16 KiB record bound")},
        "solve": solve,
        "solve_omitted": solve.is_none().then_some("context exceeds the 64 KiB record bound"),
        "app": {"version": Build::current().version_line()},
        "system": {"os": std::env::consts::OS, "arch": std::env::consts::ARCH, "environment": environment},
        "first_seen": previous.as_ref().map(|s| s.time.clone()).unwrap_or_else(|| utc(t)),
        "last_seen": utc(t),
        "count": previous.as_ref().map_or(1, |s| s.count.saturating_add(1)),
    });
    // Redact values, never serialized JSON: a Windows home path must not damage its escapes.
    fn redact(v: &mut Value) {
        match v {
            Value::String(s) => *s = without_home(s),
            Value::Array(a) => a.iter_mut().for_each(redact),
            Value::Object(o) => o.values_mut().for_each(redact),
            _ => {}
        }
    }
    redact(&mut report);
    let bytes = serde_json::to_vec_pretty(&report).ok()?;
    fs::create_dir_all(&folder).ok()?;
    let path = folder.join("report.json");
    fieldcraft_io::fs::write_atomic(&path, &bytes, false).ok()?;
    if previous.is_none() {
        WRITTEN.fetch_add(1, Ordering::Relaxed);
    }
    prune(dir, t);
    Some(path)
}

fn now() -> u64 {
    std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_secs())
}

fn now_ms() -> u64 {
    std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| u64::try_from(d.as_millis()).unwrap_or(u64::MAX))
}

/// `secs` since 1970 as (year, month, day, hour, minute, second) in UTC.
fn civil(secs: u64) -> (i64, i64, i64, u64, u64, u64) {
    let days = i64::try_from(secs / 86_400).unwrap_or(0);
    let rem = secs % 86_400;
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    (yoe + era * 400 + i64::from(m <= 2), m, d, rem / 3600, rem % 3600 / 60, rem % 60)
}

/// `2026-10-12T14:03:22Z`.
pub fn utc(secs: u64) -> String {
    let (y, m, d, h, mi, s) = civil(secs);
    format!("{y:04}-{m:02}-{d:02}T{h:02}:{mi:02}:{s:02}Z")
}

/// `20261012-140322`, for file and folder names.
pub fn stamp(secs: u64) -> String {
    let (y, m, d, h, mi, s) = civil(secs);
    format!("{y:04}{m:02}{d:02}-{h:02}{mi:02}{s:02}")
}

/// The current time as [`stamp`] gives it.
pub fn stamp_now() -> String {
    stamp(now())
}

/// A crash's identity: where it happened and what it said, numbers ignored (ids and sizes vary).
pub fn signature(location: &str, message: &str) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    let mut eat = |b: u8| {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0100_0000_01b3);
    };
    for b in location.bytes() {
        eat(b);
    }
    eat(b'|');
    let mut digit = false;
    for b in message.bytes().take(300) {
        if b.is_ascii_digit() {
            if !digit {
                eat(b'#');
            }
            digit = true;
        } else {
            digit = false;
            eat(b);
        }
    }
    h
}

fn valid_id(id: &str) -> bool {
    !id.is_empty() && id.len() <= 64 && id.chars().all(|c| c.is_ascii_alphanumeric() || c == '-')
}

/// `s` cut to at most `n` bytes (on a character boundary), saying how long it was.
pub fn limited(s: &str, n: usize) -> String {
    if s.len() <= n {
        return s.to_string();
    }
    let mut end = n;
    while !s.is_char_boundary(end) {
        end -= 1;
    }
    format!("{}… ({} bytes)", s.get(..end).unwrap_or_default(), s.len())
}

/// Serialise `v`, giving up past `limit` bytes.
fn bounded_json(v: &Value, limit: usize) -> Option<String> {
    struct W(Vec<u8>, usize);
    impl Write for W {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            if self.0.len() + b.len() > self.1 {
                return Err(std::io::Error::other("too long"));
            }
            self.0.extend_from_slice(b);
            Ok(b.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut w = W(Vec::new(), limit);
    serde_json::to_writer(&mut w, v).ok()?;
    String::from_utf8(w.0).ok()
}

/// The user's home folder in `s` shown as `~` (reports and drafts name no user folders).
pub fn without_home(s: &str) -> String {
    let home = std::env::var("HOME").or_else(|_| std::env::var("USERPROFILE")).unwrap_or_default();
    if home.len() > 1 { s.replace(&home, "~") } else { s.to_string() }
}

// ---------------------------------------------------------------------------------------------
// Panics

/// What a panic said.
#[derive(Clone, Debug)]
struct Panic {
    message: String,
    location: String,
    thread: String,
    backtrace: Vec<String>,
    signature: u64,
}

/// A backtrace as frames `symbol (file:line)`, without the panic machinery, the user's home
/// folder or the toolchain's hash.
fn frames(bt: &str) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for line in bt.lines() {
        let t = line.trim();
        if let Some(at) = t.strip_prefix("at ") {
            let mut at = without_home(at);
            if let Some(i) = at.find("/rustc/")
                && let Some(j) = at.get(i + 7..).and_then(|r| r.find('/'))
            {
                at = format!("rustc{}", at.get(i + 7 + j..).unwrap_or_default());
            }
            if let Some(last) = out.last_mut()
                && !last.ends_with(')')
            {
                last.push_str(&format!(" ({at})"));
            }
        } else if let Some((n, sym)) = t.split_once(": ")
            && n.chars().all(|c| c.is_ascii_digit())
        {
            out.push(sym.to_string());
        }
    }
    // Drop the hook and the panic machinery on top (everything up to where the panic began),
    // everything below the thread's start, and the standard library's own frames between.
    let entry = ["panicking::begin_panic", "panic::panic_any", "rust_begin_unwind", "core::panicking::", "__rust_end_short_backtrace"];
    let skip = out.iter().take(40).rposition(|f| entry.iter().any(|e| f.contains(e))).map_or(0, |i| i + 1);
    let mut v: Vec<String> = out.into_iter().skip(skip).collect();
    if let Some(end) = v.iter().position(|f| f.contains("__rust_begin_short_backtrace")) {
        v.truncate(end);
    }
    v.retain(|f| !f.contains("(rustc/"));
    v.truncate(MAX_FRAMES);
    v
}

/// Install the panic hook and start this run's journal; reports the runs that crashed without a
/// report (unclean exits). Returns the reports the user has not been told about yet, newest
/// first, for the app to offer. Call once, early (before the log file is opened: an unclean
/// exit's report takes the dead run's lines from it).
pub fn install(dir: PathBuf, build: Build) -> Vec<Summary> {
    let _ = fs::create_dir_all(dir.join("sessions"));
    let started = now();
    let started_ms = now_ms();
    for path in dead_journals(&dir) {
        unclean_exit(&dir, &build, &path);
        remove_journal(&path);
    }
    remove_stray_native_files(&dir);
    prune(&dir, started);
    // Reports the user has been told about (by id).
    let seen_file = dir.join(".seen");
    let mut seen = String::new();
    if let Ok(f) = File::open(&seen_file) {
        let _ = f.take(64 << 10).read_to_string(&mut seen);
    }
    let all = list(&dir);
    let fresh: Vec<Summary> = all.iter().filter(|s| s.kind != "solver_error" && !seen.lines().any(|l| l == s.id)).cloned().collect();
    let _ = fs::write(&seen_file, all.iter().map(|s| s.id.as_str()).collect::<Vec<_>>().join("\n"));
    let known = all.into_iter().filter_map(|s| sig_of(&s.id).map(|sig| (sig, (0, s.id)))).collect();
    let journal = Journal::open(&dir, &build, started);
    if STATE
        .set(State { dir, build, journal: Mutex::new(journal), known: Mutex::new(known), written: AtomicUsize::new(0), started, started_ms })
        .is_ok()
    {
        let previous = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            if !hook(info) {
                previous(info);
            }
        }));
    }
    fresh
}

/// What the app is doing outside commands (`"startup"`, `"running"`, a control request such as
/// `"control ui.capture"`): the journal keeps the latest, so a crash nothing catches still says
/// what came last.
pub fn activity(what: &str) {
    let Some(s) = STATE.get() else { return };
    if let Ok(mut j) = s.journal.lock() {
        // Bounded like the commands: a full journal keeps its last word.
        if j.activity_text.as_deref() == Some(what) || j.bytes > MAX_JOURNAL_BYTES {
            return;
        }
        let line = json!({ "activity": limited(what, 200), "at": utc(now()) }).to_string();
        j.write(&line);
        j.activity = Some(line);
        j.activity_text = Some(what.to_string());
    }
}

/// The app is about to start a step of getting up (`"window"`, `"device creation"`, …;
/// [`RUNNING`] once it draws): the journal keeps the latest, so a crash nothing catches (inside a
/// graphics driver, say) is still reported with the step it happened in.
pub fn phase(step: &str) {
    let Some(s) = STATE.get() else { return };
    if let Ok(mut j) = s.journal.lock() {
        if j.phase_step.as_deref() == Some(step) {
            return;
        }
        let line = json!({ "phase": limited(step, 200), "ms": now_ms().saturating_sub(s.started_ms) }).to_string();
        j.write(&line);
        j.phase = Some(line);
        j.phase_step = Some(step.to_string());
    }
}

/// The renderer, GPU and driver, once known (one short line): reports show it as
/// `system.environment`.
pub fn set_environment(text: &str) {
    let Some(s) = STATE.get() else { return };
    if let Ok(mut j) = s.journal.lock() {
        let line = json!({ "environment": limited(text, 500) }).to_string();
        j.write(&line);
        j.environment = Some(line);
    }
}

/// Where a native crash handler writes what it caught for this run: a small JSON file
/// (`code`, `name`, `address`, `module`, `offset`, `detail`, `thread`) and a minidump. The next
/// start reports both with this run's unclean exit. `None` before [`install`] or without a journal.
pub fn native_paths() -> Option<(PathBuf, PathBuf)> {
    let s = STATE.get()?;
    let j = s.journal.lock().ok()?;
    j.path.as_ref().map(|p| native_paths_for(p))
}

fn native_paths_for(journal: &Path) -> (PathBuf, PathBuf) {
    (journal.with_extension("native.json"), journal.with_extension("dmp"))
}

/// Native crash files whose journal is gone (their report was written, or never will be).
fn remove_stray_native_files(dir: &Path) {
    for path in fs::read_dir(dir.join("sessions")).into_iter().flatten().flatten().map(|e| e.path()).take(200) {
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or_default();
        let stem = name.strip_suffix(".native.json").or_else(|| name.strip_suffix(".dmp"));
        if let Some(stem) = stem
            && !path.with_file_name(format!("{stem}.jsonl")).exists()
        {
            let _ = fs::remove_file(&path);
        }
    }
}

/// The app is closing normally: its journal goes, and repeat counts are written.
pub fn end_session() {
    let Some(s) = STATE.get() else { return };
    if let Ok(mut j) = s.journal.lock()
        && let Some(path) = j.path.take()
    {
        j.file = None;
        remove_journal(&path);
    }
    let counts: Vec<(u64, String)> = s.known.lock().map(|k| k.values().filter(|(n, _)| *n > 0).cloned().collect()).unwrap_or_default();
    for (n, id) in counts {
        bump(&s.dir, &id, n);
    }
}

/// A journal and the model copy it starts from.
fn remove_journal(path: &Path) {
    let _ = fs::remove_file(path);
    let _ = fs::remove_file(base_copy(path));
}

/// Where a journal keeps the model it starts from, when that is no file of its own.
fn base_copy(journal: &Path) -> PathBuf {
    journal.with_extension("base.fieldcraft")
}

/// The panic hook. Returns true when the panic was recorded for a command (which will catch
/// it), so the default "thread panicked" message is not printed twice.
fn hook(info: &std::panic::PanicHookInfo<'_>) -> bool {
    let Some(state) = STATE.get() else { return false };
    if IN_HOOK.swap(true, Ordering::SeqCst) {
        return false;
    }
    let message = limited(info.payload_as_str().unwrap_or("(not a text message)"), MAX_MESSAGE);
    let location = info.location().map_or_else(|| "(unknown)".to_string(), |l| format!("{}:{}", l.file(), l.line()));
    log::error!("panic at {location}: {message}");
    let sig = signature(&location, &message);
    let in_command = IN_COMMAND.with(Cell::get) > 0;
    // Already reported (this run or before): count it, nothing else.
    let repeat = match state.known.lock() {
        Ok(mut k) => match k.get_mut(&sig) {
            Some(e) => {
                e.0 += 1;
                true
            }
            None => false,
        },
        Err(_) => true,
    };
    if repeat {
        IN_HOOK.store(false, Ordering::SeqCst);
        return in_command;
    }
    let thread = std::thread::current().name().unwrap_or("(unnamed)").to_string();
    let backtrace = frames(&Backtrace::force_capture().to_string());
    let p = Panic { message, location, thread, backtrace, signature: sig };
    if in_command {
        // Reported with the model and the command when the outermost command ends.
        PENDING.with(|q| {
            let mut q = q.borrow_mut();
            if q.len() < 8 {
                q.push(p);
            }
        });
    } else {
        // Outside any command: what the journal knows has to do.
        let extra = from_journal(read_live_journal(state), None);
        if let Some(id) = write_report(state, "panic", Some(&p), extra) {
            eprintln!("FieldCraft: internal error; crash report {}", state.dir.join(&id).display());
            if let Ok(mut j) = state.journal.try_lock() {
                j.record(&json!({ "reported": id }));
            }
        }
    }
    IN_HOOK.store(false, Ordering::SeqCst);
    in_command
}

/// Marks a command running on this thread; panics inside it wait for [`command_done`].
pub(crate) struct Scope(bool);

impl Scope {
    pub(crate) fn enter() -> Scope {
        if !enabled() {
            return Scope(false);
        }
        IN_COMMAND.with(|c| c.set(c.get() + 1));
        Scope(true)
    }

    /// The outermost command on this thread.
    pub(crate) fn top_level(&self) -> bool {
        self.0 && IN_COMMAND.with(Cell::get) == 1
    }
}

impl Drop for Scope {
    fn drop(&mut self) {
        if self.0 {
            IN_COMMAND.with(|c| c.set(c.get().saturating_sub(1)));
        }
    }
}

/// A top-level command of the app's own model is about to run: the journal records it. Returns
/// whether it did.
pub(crate) fn command_begin(scope: &Scope, s: &Session, id: &str, params: &Value) -> bool {
    let Some(state) = STATE.get() else { return false };
    if !s.crash_journal || !scope.top_level() {
        return false;
    }
    let Ok(mut j) = state.journal.lock() else { return false };
    if j.document != Some(s.opened()) {
        // The model on screen changed outside a recorded command (the app's first command).
        rebase_on_session(&mut j, s);
    }
    j.run(id, params);
    true
}

/// A top-level command ended (`scope` already dropped): the journal records how; when it
/// opened another model the journal starts again from it. Panics caught inside it are reported
/// with the session (restored to how it was before the command after an internal error) and
/// the command, which together reproduce them.
pub(crate) fn command_done(s: &Session, journaled: bool, changed: bool, id: &str, params: &Value, result: &crate::Result<Value>) {
    let Some(state) = STATE.get() else { return };
    if journaled && let Ok(mut j) = state.journal.lock() {
        j.done(result.as_ref().map(|_| ()).map_err(ToString::to_string), changed);
        if j.document != Some(s.opened()) {
            rebase_on_session(&mut j, s);
        }
    }
    if IN_COMMAND.with(Cell::get) > 0 {
        return;
    }
    let pending: Vec<Panic> = PENDING.with(|q| std::mem::take(&mut *q.borrow_mut()));
    if pending.is_empty() {
        return;
    }
    // The command failed with the panic: the session is as it was before it. A panic caught
    // inside the command (and turned into an error or worked around) left the command's own
    // outcome; it is reported, but as recovered, with the model as it is now.
    let crashed = matches!(result, Err(EngineError::Internal { .. }));
    let kind = if crashed { "command_panic" } else { "recovered_panic" };
    let mut notes = Vec::new();
    let model = match s.file_bytes(None) {
        Ok(b) if b.len() <= MAX_MODEL_BYTES => Some(b),
        Ok(b) => {
            notes.push(format!("the model ({} MB) is larger than {} MB and was not copied", b.len() >> 20, MAX_MODEL_BYTES >> 20));
            None
        }
        Err(e) => {
            notes.push(format!("the model could not be written: {e}"));
            None
        }
    };
    let recent = if journaled { read_live_journal(state).recent } else { Vec::new() };
    for p in pending {
        let extra = Extra {
            command: Some(json!({ "command": id, "params": params })),
            model: model.clone(),
            replay: (crashed && model.is_some()).then(|| json!({ "commands": [{ "command": id, "params": params }] })),
            recent: recent.clone(),
            base: json!(if crashed { "the model as it was when the command started" } else { "the model after the command" }),
            notes: notes.clone(),
            fields: serde_json::Map::new(),
            attach: Vec::new(),
            log_pid: Some(std::process::id()),
        };
        write_report(state, kind, Some(&p), extra);
    }
}

/// An autosave just wrote the session to `file`: replays start from there.
pub(crate) fn autosaved(s: &Session, file: &Path) {
    let Some(state) = STATE.get() else { return };
    if !s.crash_journal {
        return;
    }
    let Ok(mut j) = state.journal.lock() else { return };
    j.rebase(&json!({ "autosave": file.to_string_lossy() }), s.opened());
}

/// Start the journal again from the session as it is now: its file when that holds it, a new
/// model when it is empty, else a copy beside the journal.
fn rebase_on_session(j: &mut Journal, s: &Session) {
    let opened = s.opened();
    if !s.is_dirty()
        && let Some(p) = &s.path
    {
        j.rebase(&json!({ "file": p }), opened);
        return;
    }
    if s.journal.is_empty() && s.model.nodes().is_empty() && s.model.elements().is_empty() && s.design.is_none() && s.surfaces.is_empty() {
        j.rebase(&json!("new"), opened);
        return;
    }
    let base = j.path.as_deref().map(base_copy);
    let copied = base.as_ref().is_some_and(|b| match s.file_bytes(None) {
        Ok(bytes) if bytes.len() <= MAX_MODEL_BYTES => fieldcraft_io::fs::write_atomic(b, &bytes, false).is_ok(),
        _ => false,
    });
    match base.filter(|_| copied) {
        Some(b) => j.rebase(&json!({ "copy": b.to_string_lossy() }), opened),
        None => j.rebase(&json!("unknown"), opened),
    }
}

// ---------------------------------------------------------------------------------------------
// The journal

struct Journal {
    path: Option<PathBuf>,
    file: Option<File>,
    bytes: u64,
    full: bool,
    /// The document ([`Session::opened`]) the journal follows.
    document: Option<u64>,
    header: String,
    /// Where the running command's line starts (to take it back when it changed nothing).
    running_at: Option<u64>,
    /// The latest activity's line and text (kept across rebases).
    activity: Option<String>,
    activity_text: Option<String>,
    /// The latest startup step's line and name, and the renderer's line (kept across rebases).
    phase: Option<String>,
    phase_step: Option<String>,
    environment: Option<String>,
}

impl Journal {
    fn open(dir: &Path, build: &Build, started: u64) -> Journal {
        let seq = SEQ.fetch_add(1, Ordering::Relaxed);
        let path = dir.join("sessions").join(format!("{}-{}-{seq}.jsonl", stamp(started), std::process::id()));
        // Append-only handles cannot truncate on Windows; seek explicitly for each write.
        let file = OpenOptions::new().create(true).truncate(false).read(true).write(true).open(&path).ok().filter(|f| f.try_lock().is_ok());
        let header = json!({ "start": utc(started), "pid": std::process::id(), "version": build.version, "commit": build.commit }).to_string();
        let mut j = Journal {
            path: file.as_ref().map(|_| path),
            file,
            bytes: 0,
            full: false,
            document: None,
            header,
            running_at: None,
            activity: None,
            activity_text: None,
            phase: None,
            phase_step: None,
            environment: None,
        };
        j.rebase(&json!("new"), 0);
        j.document = None;
        j
    }

    fn write(&mut self, line: &str) {
        let Some(f) = &mut self.file else { return };
        if f.seek(SeekFrom::End(0)).and_then(|_| f.write_all(line.as_bytes())).and_then(|()| f.write_all(b"\n")).is_err() {
            self.file = None;
        }
        self.bytes += line.len() as u64 + 1;
    }

    fn record(&mut self, v: &Value) {
        self.write(&v.to_string());
    }

    fn rebase(&mut self, base: &Value, document: u64) {
        if let Some(f) = &self.file
            && f.set_len(0).is_err()
        {
            self.file = None;
        }
        self.bytes = 0;
        self.full = false;
        self.running_at = None;
        self.document = Some(document);
        let h = self.header.clone();
        self.write(&h);
        self.record(&json!({ "base": base }));
        for line in [self.activity.clone(), self.phase.clone(), self.environment.clone()].into_iter().flatten() {
            self.write(&line);
        }
    }

    fn run(&mut self, id: &str, params: &Value) {
        if self.full {
            return;
        }
        if self.bytes > MAX_JOURNAL_BYTES {
            self.full = true;
            self.record(&json!({ "full": true }));
            return;
        }
        self.running_at = Some(self.bytes);
        let id = Value::String(id.to_string()).to_string();
        let line = match bounded_json(params, MAX_PARAMS_BYTES) {
            Some(p) => format!("{{\"run\":{id},\"params\":{p}}}"),
            None => format!("{{\"run\":{id},\"params_omitted\":true}}"),
        };
        self.write(&line);
    }

    fn done(&mut self, outcome: std::result::Result<(), String>, changed: bool) {
        let at = self.running_at.take();
        if self.full {
            return;
        }
        match outcome {
            // Nothing changed: take the command's line back (one truncate, no write).
            Ok(()) if !changed => {
                if let (Some(at), Some(f)) = (at, &self.file) {
                    if f.set_len(at).is_ok() {
                        self.bytes = at;
                    } else {
                        self.file = None;
                    }
                }
            }
            Ok(()) => self.write("{\"ok\":true}"),
            Err(e) => self.record(&json!({ "ok": false, "error": limited(&e, 300) })),
        }
    }
}

/// What a journal says.
#[derive(Default)]
struct JournalRead {
    header: Value,
    base: Value,
    /// Commands since the base that changed something, in order.
    replay: Vec<Value>,
    /// The command that was running, if one was.
    in_flight: Option<Value>,
    recent: Vec<Value>,
    complete: bool,
    reported: bool,
    /// The latest activity outside commands.
    activity: Option<String>,
    /// The latest startup step and how many milliseconds after the start it began.
    phase: Option<(String, u64)>,
    /// The renderer, once known.
    environment: Option<String>,
}

// Read through the lock-owning file object: on Windows an independently opened handle cannot
// read a live journal's exclusively locked byte range.
fn read_live_journal(state: &State) -> JournalRead {
    state.journal.try_lock().ok().and_then(|j| j.file.as_ref()?.try_clone().ok().map(read_journal_file)).unwrap_or_default()
}

fn read_journal(path: &Path) -> JournalRead {
    File::open(path).ok().map(read_journal_file).unwrap_or_default()
}

fn read_journal_file(mut file: File) -> JournalRead {
    let mut text = String::new();
    if file.seek(SeekFrom::Start(0)).is_ok() {
        let _ = file.take(MAX_JOURNAL_BYTES + (64 << 10)).read_to_string(&mut text);
    }
    parse_journal(&text)
}

fn parse_journal(text: &str) -> JournalRead {
    let mut r = JournalRead { complete: true, ..Default::default() };
    let mut running: Option<Value> = None;
    for (i, line) in text.lines().enumerate() {
        let Ok(v) = serde_json::from_str::<Value>(line) else { continue };
        if i == 0 {
            r.header = v;
            continue;
        }
        if let Some(a) = v.get("activity").and_then(Value::as_str) {
            r.activity = Some(a.to_string());
            continue;
        }
        if let Some(p) = v.get("phase").and_then(Value::as_str) {
            r.phase = Some((p.to_string(), v.get("ms").and_then(Value::as_u64).unwrap_or(0)));
            continue;
        }
        if let Some(e) = v.get("environment").and_then(Value::as_str) {
            r.environment = Some(e.to_string());
            continue;
        }
        r.reported = v.get("reported").is_some();
        if let Some(b) = v.get("base") {
            r.base = b.clone();
            r.replay.clear();
            r.recent.clear();
            r.complete = true;
            running = None;
        } else if let Some(id) = v.get("run") {
            let c = match v.get("params") {
                Some(p) => json!({ "command": id, "params": p }),
                None => {
                    r.complete = false;
                    json!({ "command": id, "params_omitted": true })
                }
            };
            running = Some(c);
        } else if let Some(ok) = v.get("ok").and_then(Value::as_bool) {
            if let Some(c) = running.take() {
                let mut shown = c.clone();
                shown["ok"] = json!(ok);
                if let Some(e) = v.get("error") {
                    shown["error"] = e.clone();
                }
                r.recent.push(shown);
                if ok {
                    r.replay.push(c);
                }
            }
        } else if v.get("full").is_some() {
            r.complete = false;
        }
    }
    r.in_flight = running;
    let n = r.recent.len();
    r.recent.drain(..n.saturating_sub(RECENT));
    r
}

/// Journals of apps that are not running any more (nobody holds their lock).
fn dead_journals(dir: &Path) -> Vec<PathBuf> {
    fs::read_dir(dir.join("sessions"))
        .into_iter()
        .flatten()
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.extension().and_then(|x| x.to_str()) == Some("jsonl"))
        .take(100)
        .filter(|p| File::open(p).is_ok_and(|f| f.try_lock().is_ok()))
        .collect()
}

fn unclean_exit(dir: &Path, build: &Build, journal: &Path) {
    let j = read_journal(journal);
    let (native_file, dump_file) = native_paths_for(journal);
    let native = fs::read(&native_file).ok().and_then(|b| serde_json::from_slice::<Value>(&b).ok()).filter(Value::is_object);
    if j.reported {
        return;
    }
    let was = j.header.get("version").and_then(Value::as_str).unwrap_or("?").to_string();
    let when = j.header.get("start").and_then(Value::as_str).unwrap_or("?").to_string();
    let pid = j.header.get("pid").and_then(Value::as_u64).and_then(|p| u32::try_from(p).ok());
    let doing = j.in_flight.as_ref().and_then(|c| c.get("command")).and_then(Value::as_str).map(str::to_string);
    // A run that never got to drawing died while starting up: say at which step.
    let starting = j.phase.as_ref().filter(|(p, _)| p != RUNNING).cloned();
    let during = match (&starting, &doing, &j.activity) {
        (Some((p, ms)), ..) => format!("during startup ({p}, {:.1} s after launch)", std::time::Duration::from_millis(*ms).as_secs_f64()),
        (None, Some(c), _) => format!("while running `{c}`"),
        (None, None, Some(a)) => format!("between commands (last: {a})"),
        (None, None, None) => "between commands".to_string(),
    };
    // What the run was doing, for the signature: the step, the command or the last activity.
    let what = starting
        .as_ref()
        .map(|(p, _)| format!("startup: {p}"))
        .or_else(|| doing.clone())
        .or_else(|| j.activity.clone())
        .unwrap_or_else(|| "idle".into());
    let p = match &native {
        Some(n) => native_panic(n, &was, &when, &during, &what),
        None => Panic {
            message: format!("FieldCraft {was} (started {when}) ended without a clean exit {during}"),
            location: match (&starting, &doing) {
                (Some((p, _)), _) => format!("unclean exit during startup: {p}"),
                (None, c) => format!("unclean exit in {}", c.as_deref().unwrap_or("(between commands)")),
            },
            thread: "(unknown)".into(),
            backtrace: Vec::new(),
            signature: signature("unclean exit", &what),
        },
    };
    // A temporary state just for this report (the app's own is not set up yet).
    let state = State {
        dir: dir.to_path_buf(),
        build: build.clone(),
        journal: Mutex::new(Journal {
            path: None,
            file: None,
            bytes: 0,
            full: true,
            document: None,
            header: String::new(),
            running_at: None,
            activity: None,
            activity_text: None,
            phase: None,
            phase_step: None,
            environment: None,
        }),
        known: Mutex::new(list(dir).into_iter().filter_map(|s| sig_of(&s.id).map(|sig| (sig, (0, s.id)))).collect()),
        written: AtomicUsize::new(0),
        started: now(),
        started_ms: now_ms(),
    };
    if let Ok(mut k) = state.known.lock()
        && let Some(e) = k.get_mut(&p.signature)
    {
        e.0 += 1;
        let (n, id) = e.clone();
        drop(k);
        bump(dir, &id, n);
        return;
    }
    let mut extra = from_journal(j, pid);
    if let Some(n) = native {
        extra.fields.insert("native".into(), n);
        if dump_file.is_file() {
            extra.attach.push(("minidump".into(), dump_file, MINIDUMP_FILE.into()));
        }
    }
    write_report(&state, "unclean_exit", Some(&p), extra);
}

/// The report of a crash a native exception filter caught (`native`: what it wrote, untrusted:
/// every field is cut to length and may be missing).
fn native_panic(native: &Value, was: &str, when: &str, during: &str, what: &str) -> Panic {
    let field = |k: &str| native.get(k).and_then(Value::as_str).map(|v| limited(v, 300));
    let code = field("code").unwrap_or_else(|| "?".into());
    let name = field("name").unwrap_or_else(|| "exception".into());
    let module = field("module").unwrap_or_else(|| "(unknown module)".into());
    let short = module.rsplit(['\\', '/']).next().unwrap_or(&module).to_string();
    let offset = field("offset").unwrap_or_else(|| "?".into());
    let at = format!("{short}+{offset}");
    let mut backtrace = vec![format!("{code} {name}{}", field("detail").map(|d| format!(" ({d})")).unwrap_or_default())];
    backtrace.push(format!("at {module}+{offset} (address {})", field("address").unwrap_or_else(|| "?".into())));
    if let Some(t) = native.get("thread").and_then(Value::as_u64) {
        backtrace.push(format!("on thread {t}"));
    }
    Panic {
        message: format!("FieldCraft {was} (started {when}) crashed with {name} ({code}) in {at} {during}"),
        location: format!("native exception in {at}"),
        thread: native.get("thread").and_then(Value::as_u64).map_or_else(|| "(unknown)".into(), |t| format!("native thread {t}")),
        backtrace,
        // The place: the module and the offset into it (signature() ignores the digits).
        signature: signature("native exception", &format!("{code} {short} {what}")),
    }
}

// ---------------------------------------------------------------------------------------------
// Reports

struct Extra {
    command: Option<Value>,
    model: Option<Vec<u8>>,
    replay: Option<Value>,
    recent: Vec<Value>,
    base: Value,
    notes: Vec<String>,
    /// More top-level fields for the report (`last_activity`, `startup_phase`, `native`).
    fields: serde_json::Map<String, Value>,
    /// Files moved into the report's folder: (key in `files`, source, name in the folder).
    attach: Vec<(String, PathBuf, String)>,
    /// The process whose log lines go into the report.
    log_pid: Option<u32>,
}

/// Model and replay from a journal: the base model (autosave, the model's file or the copy
/// taken when it was opened) and the commands since, ending with the one in flight.
fn from_journal(j: JournalRead, log_pid: Option<u32>) -> Extra {
    let mut notes = Vec::new();
    let file = ["autosave", "file", "copy"].iter().find_map(|k| j.base.get(*k)).and_then(Value::as_str).map(PathBuf::from);
    let model = match &file {
        Some(f) => match fs::metadata(f) {
            Ok(m) if m.len() <= MAX_MODEL_BYTES as u64 => fs::read(f).ok(),
            Ok(_) => {
                notes.push(format!(
                    "the model ({}) is larger than {} MB and was not copied",
                    without_home(&f.to_string_lossy()),
                    MAX_MODEL_BYTES >> 20
                ));
                None
            }
            Err(_) => {
                notes.push(format!("the model the commands start from ({}) is gone", without_home(&f.to_string_lossy())));
                None
            }
        },
        None => None,
    };
    let fresh = j.base.as_str() == Some("new");
    let mut commands = j.replay;
    if let Some(c) = &j.in_flight {
        commands.push(c.clone());
    }
    if !j.complete {
        notes.push("the journal was cut short (too many or too large commands): the replay is incomplete".into());
    }
    if j.base.as_str() == Some("unknown") {
        notes.push("the model could not be copied when it was opened; it cannot be rebuilt from the journal".into());
    }
    if commands.iter().any(|c| matches!(c.get("command").and_then(Value::as_str), Some("edit.undo" | "edit.redo"))) {
        notes.push("the replay starts without undo history: edit.undo and edit.redo may not do what they did".into());
    }
    let replay = (model.is_some() || fresh) && j.complete && !commands.is_empty();
    let mut fields = serde_json::Map::new();
    if let Some(a) = j.activity {
        fields.insert("last_activity".into(), json!(a));
    }
    if let Some((p, ms)) = j.phase.filter(|(p, _)| p != RUNNING) {
        fields.insert("startup_phase".into(), json!({ "step": p, "ms": ms }));
    }
    if let Some(e) = j.environment {
        fields.insert("environment".into(), json!(e));
    }
    Extra {
        command: j.in_flight,
        model,
        replay: replay.then(|| json!({ "commands": commands })),
        recent: j.recent,
        base: match file {
            Some(f) => json!({ "file": without_home(&f.to_string_lossy()) }),
            None => j.base,
        },
        notes,
        fields,
        attach: Vec::new(),
        log_pid,
    }
}

fn write_report(state: &State, kind: &str, p: Option<&Panic>, extra: Extra) -> Option<String> {
    let sig = p.map_or(0, |p| p.signature);
    let mut k = state.known.lock().ok()?;
    if let Some(e) = k.get_mut(&sig) {
        e.0 += 1;
        return None;
    }
    if state.written.load(Ordering::SeqCst) >= MAX_NEW_PER_SESSION {
        return None;
    }
    state.written.fetch_add(1, Ordering::SeqCst);
    let t = now();
    let id = format!("{}-{sig:016x}", stamp(t));
    k.insert(sig, (0, id.clone()));
    drop(k);
    let folder = state.dir.join(&id);
    fs::create_dir_all(&folder).ok()?;
    let mut files = serde_json::Map::new();
    let mut notes = extra.notes;
    if let Some(m) = &extra.model {
        if m.len() <= MAX_MODEL_BYTES {
            if fieldcraft_io::fs::write_atomic(&folder.join(MODEL_FILE), m, false).is_ok() {
                files.insert("model".into(), json!(MODEL_FILE));
            }
        } else {
            notes.push(format!("the model is larger than {} MB and was not copied", MAX_MODEL_BYTES >> 20));
        }
    }
    if let Some(r) = &extra.replay
        && let Ok(b) = serde_json::to_vec_pretty(r)
        && fieldcraft_io::fs::write_atomic(&folder.join("replay.json"), &b, false).is_ok()
    {
        files.insert("replay".into(), json!("replay.json"));
    }
    if let Some(log) = crate::logfile::tail(extra.log_pid, LOG_TAIL_BYTES).filter(|l| !l.is_empty())
        && fieldcraft_io::fs::write_atomic(&folder.join("log.txt"), log.as_bytes(), false).is_ok()
    {
        files.insert("log".into(), json!("log.txt"));
    }
    for (key, from, name) in &extra.attach {
        let to = folder.join(name);
        if fs::rename(from, &to).is_ok() || fs::copy(from, &to).is_ok() {
            files.insert(key.clone(), json!(name));
        }
    }
    let reproduce = files.contains_key("replay").then(|| format!("fieldcraft-cli crashes replay {id}"));
    let mut report = json!({
        "format": FORMAT,
        "id": id,
        "kind": kind,
        "signature": format!("{sig:016x}"),
        "count": 1,
        "first_seen": utc(t),
        "last_seen": utc(t),
        "app": { "version": state.build.version, "commit": state.build.commit },
        "system": { "os": std::env::consts::OS, "arch": std::env::consts::ARCH },
        "session": { "uptime_s": t.saturating_sub(state.started) },
        "message": p.map(|p| without_home(&p.message)),
        "location": p.map(|p| without_home(&p.location)),
        "thread": p.map(|p| p.thread.clone()),
        "backtrace": p.map(|p| p.backtrace.clone()).unwrap_or_default(),
        "command": extra.command,
        "base": extra.base,
        "recent_commands": extra.recent,
        "files": files,
        "reproduce": reproduce,
        "notes": notes,
    });
    if let Some(r) = report.as_object_mut() {
        for (k, v) in extra.fields {
            // The renderer goes with the system it ran on.
            if k == "environment"
                && let Some(system) = r.get_mut("system").and_then(Value::as_object_mut)
            {
                system.insert(k, v);
                continue;
            }
            r.entry(k).or_insert(v);
        }
    }
    let bytes = serde_json::to_vec_pretty(&report).ok()?;
    fieldcraft_io::fs::write_atomic(&folder.join("report.json"), &bytes, false).ok()?;
    prune(&state.dir, t);
    Some(id)
}

/// Add `n` repeats to a report.
fn bump(dir: &Path, id: &str, n: u64) {
    if !valid_id(id) || n == 0 {
        return;
    }
    let path = dir.join(id).join("report.json");
    let Some(mut v) = read_report_file(&path) else { return };
    let count = v.get("count").and_then(Value::as_u64).unwrap_or(1).saturating_add(n);
    v["count"] = json!(count);
    v["last_seen"] = json!(utc(now()));
    if let Ok(b) = serde_json::to_vec_pretty(&v) {
        let _ = fieldcraft_io::fs::write_atomic(&path, &b, false);
    }
}

fn read_report_file(path: &Path) -> Option<Value> {
    let f = File::open(path).ok()?;
    let mut s = String::new();
    f.take(4 << 20).read_to_string(&mut s).ok()?;
    serde_json::from_str(&s).ok().filter(Value::is_object)
}

fn sig_of(id: &str) -> Option<u64> {
    id.rsplit('-').next().and_then(|h| u64::from_str_radix(h, 16).ok())
}

/// Seconds since 1970 at which report `id` was made (from its file's time).
fn report_time(dir: &Path, id: &str) -> u64 {
    fs::metadata(dir.join(id).join("report.json"))
        .and_then(|m| m.modified())
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map_or(0, |d| d.as_secs())
}

fn folder_bytes(p: &Path) -> u64 {
    fs::read_dir(p).into_iter().flatten().flatten().filter_map(|e| e.metadata().ok()).map(|m| m.len()).sum()
}

/// Report folders in `dir`, oldest first.
fn report_ids(dir: &Path) -> Vec<String> {
    let mut v: Vec<String> = fs::read_dir(dir)
        .into_iter()
        .flatten()
        .flatten()
        .filter(|e| e.file_type().is_ok_and(|t| t.is_dir()))
        .filter_map(|e| e.file_name().to_str().map(str::to_string))
        .filter(|n| valid_id(n) && n.len() == 32 && sig_of(n).is_some())
        .take(1000)
        .collect();
    v.sort();
    v
}

/// Keep the reports within the limits: drop expired ones, then recovered panics, then the
/// oldest, until at most [`MAX_REPORTS`] and [`MAX_TOTAL_BYTES`] remain.
pub fn prune(dir: &Path, now_secs: u64) {
    let mut keep: Vec<(String, bool, u64)> = Vec::new();
    for id in report_ids(dir) {
        let folder = dir.join(&id);
        let t = report_time(dir, &id);
        if (t > 0 && now_secs.saturating_sub(t) > MAX_AGE_DAYS * 86_400) || !folder.join("report.json").is_file() {
            let _ = fs::remove_dir_all(&folder);
            continue;
        }
        let recovered =
            read_report_file(&folder.join("report.json")).is_some_and(|v| v.get("kind").and_then(Value::as_str) == Some("recovered_panic"));
        keep.push((id, recovered, folder_bytes(&folder)));
    }
    let over = |k: &[(String, bool, u64)]| k.len() > MAX_REPORTS || k.iter().map(|e| e.2).sum::<u64>() > MAX_TOTAL_BYTES;
    while over(&keep) {
        let i = keep.iter().position(|e| e.1).unwrap_or(0);
        if i >= keep.len() {
            break;
        }
        let (id, _, _) = keep.remove(i);
        let _ = fs::remove_dir_all(dir.join(id));
    }
}

/// The reports in `dir`, newest first.
pub fn list(dir: &Path) -> Vec<Summary> {
    let mut out: Vec<Summary> = report_ids(dir)
        .into_iter()
        .rev()
        .filter_map(|id| {
            let folder = dir.join(&id);
            let v = read_report_file(&folder.join("report.json"))?;
            let s = |k: &str| v.get(k).and_then(Value::as_str).unwrap_or_default().to_string();
            Some(Summary {
                kind: s("kind"),
                time: s("first_seen"),
                message: limited(&s("message"), 300),
                command: v.get("command").and_then(|c| c.get("command")).and_then(Value::as_str).map(str::to_string),
                count: v.get("count").and_then(Value::as_u64).unwrap_or(1),
                version: v.get("app").and_then(|a| a.get("version")).and_then(Value::as_str).unwrap_or_default().to_string(),
                path: folder.to_string_lossy().into_owned(),
                replayable: folder.join("replay.json").is_file(),
                id,
            })
        })
        .collect();
    out.truncate(MAX_REPORTS * 2);
    out
}

/// One report's `report.json`, with the paths of its files.
pub fn read(dir: &Path, id: &str) -> std::result::Result<Value, String> {
    if !valid_id(id) {
        return Err(format!("no crash report `{}`", limited(id, 64)));
    }
    let folder = dir.join(id);
    let mut v = read_report_file(&folder.join("report.json")).ok_or_else(|| format!("no crash report `{id}`"))?;
    v["path"] = json!(folder.to_string_lossy());
    for (k, f) in [("model_path", MODEL_FILE), ("replay_path", "replay.json"), ("log_path", "log.txt"), ("minidump_path", MINIDUMP_FILE)] {
        if folder.join(f).is_file() {
            v[k] = json!(folder.join(f).to_string_lossy());
        }
    }
    Ok(v)
}

/// Delete one report, or all of them.
pub fn clear(dir: &Path, id: Option<&str>) -> usize {
    let ids = match id {
        Some(i) if valid_id(i) && report_ids(dir).iter().any(|r| r == i) => vec![i.to_string()],
        Some(_) => Vec::new(),
        None => report_ids(dir),
    };
    ids.iter().filter(|i| fs::remove_dir_all(dir.join(i)).is_ok()).count()
}

// ---------------------------------------------------------------------------------------------
// Replay

/// What replaying a report did.
#[derive(Debug)]
pub struct Replay {
    /// Each command and how it ended (`Ok` with nothing, or the error).
    pub steps: Vec<(String, std::result::Result<(), String>)>,
    /// Internal errors seen: (location, message).
    pub panics: Vec<(String, String)>,
    /// One of them is the reported one (the same signature).
    pub reproduced: bool,
    /// The session the replay ends with.
    pub session: Session,
}

static REPLAY_PANICS: Mutex<Vec<(String, String)>> = Mutex::new(Vec::new());

/// Replay a report ([`read`]): open its model (or start from an empty one), run its commands,
/// and say whether the same internal error happens again. Files the commands write can go only
/// into `scratch` (reads go anywhere: the model's geometry is where it was). The panic hook is
/// replaced while it runs: for the CLI and tests, never the app.
pub fn replay(report: &Value, scratch: &Path) -> std::result::Result<Replay, String> {
    let id = report.get("id").and_then(Value::as_str).unwrap_or_default();
    let replay_path = report.get("replay_path").and_then(Value::as_str).ok_or_else(|| format!("report {id} cannot be replayed (see its notes)"))?;
    let script: Value = fieldcraft_io::fs::read_capped(Path::new(replay_path), MAX_JOURNAL_BYTES * 4)
        .ok()
        .and_then(|b| serde_json::from_slice(&b).ok())
        .ok_or_else(|| format!("{replay_path}: not a replay"))?;
    let commands = script.get("commands").and_then(Value::as_array).cloned().unwrap_or_default();
    let mut s = Session::default();
    if let Some(m) = report.get("model_path").and_then(Value::as_str) {
        s.execute("file.open", &json!({ "path": m })).map_err(|e| format!("{m}: {e}"))?;
    }
    s.file_access = fieldcraft_io::roots::Roots::writes_within(scratch).map_err(|e| e.to_string())?;
    if let Ok(mut p) = REPLAY_PANICS.lock() {
        p.clear();
    }
    let previous = std::panic::take_hook();
    std::panic::set_hook(Box::new(|info| {
        let loc = info.location().map_or_else(|| "(unknown)".to_string(), |l| format!("{}:{}", l.file(), l.line()));
        let msg = info.payload_as_str().unwrap_or("(not a text message)").to_string();
        if let Ok(mut p) = REPLAY_PANICS.lock() {
            p.push((loc, limited(&msg, MAX_MESSAGE)));
        }
    }));
    let mut steps = Vec::new();
    for c in commands.iter().take(100_000) {
        let cmd = c.get("command").and_then(Value::as_str).unwrap_or_default().to_string();
        let r = s.execute(&cmd, c.get("params").unwrap_or(&Value::Null)).map(|_| ()).map_err(|e| e.to_string());
        steps.push((cmd, r));
    }
    std::panic::set_hook(previous);
    let panics = REPLAY_PANICS.lock().map(|p| p.clone()).unwrap_or_default();
    let want = report.get("signature").and_then(Value::as_str).unwrap_or_default();
    let reproduced = panics.iter().any(|(l, m)| format!("{:016x}", signature(l, m)) == want);
    Ok(Replay { steps, panics, reproduced, session: s })
}

#[cfg(test)]
#[path = "crash_tests.rs"]
mod tests;
