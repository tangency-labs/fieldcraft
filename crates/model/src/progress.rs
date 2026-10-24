//! Progress of a long-running command (a solve, a mesher), reported from deep inside it.
//!
//! A host that runs a command in the background (the desktop's jobs, the control channel)
//! makes a [`Progress`], installs it on the worker thread with [`scoped`] and reads
//! [`Progress::snapshot`] from another thread while the command runs. The command's layers
//! report where they are without taking a parameter for it: the engine names the step
//! ([`step`]), solvers and meshers their stage and counts ([`stage`], [`count`]). With no
//! [`Progress`] installed (scripts, tests, the CLI) every report is a no-op.
//!
//! The reporter is per thread: work handed to a thread pool installs the caller's
//! ([`current`]) inside it if it reports from there. Reporting never changes a result.
//!
//! The same [`Progress`] carries the host's request to stop ([`Progress::cancel`]): the long
//! loops of the solvers and meshers ask [`check`] at their natural points (an increment, a time
//! step, a refinement round, a supernode) and return [`crate::Error::Cancelled`] once it is
//! set. Cancellation is cooperative: the work stops at the next such point, and the caller
//! discards what it was building, so a cancelled command changes nothing.

use std::cell::RefCell;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

use serde::Serialize;

/// Where a command is.
#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct State {
    /// The step being solved: its name, index from 1 and the number of steps.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub step: Option<StepAt>,
    /// What it is doing now ("Delaunay", "Increment", "Improvement pass").
    pub stage: String,
    /// How far into the stage (increments, passes, time steps, tetrahedra made).
    pub done: f64,
    /// Where the stage ends, when that is known (`done / total` is then a fraction).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub total: Option<f64>,
    /// What `done` counts ("increments", "tets", "s").
    #[serde(skip_serializing_if = "str::is_empty")]
    pub unit: &'static str,
}

/// A step of a multi-step solve.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct StepAt {
    pub name: String,
    pub index: usize,
    pub count: usize,
}

/// Platform-neutral metadata observed while preparing a solver backend. None means that
/// setup has not observed the value; the report never opens a device to fill missing data.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub struct BackendDiagnostic {
    pub phase: String,
    pub backend: Option<String>,
    pub adapter: Option<String>,
    pub driver: Option<String>,
    pub driver_info: Option<String>,
    pub adapter_limits: Option<BackendLimits>,
    pub device_limits: Option<BackendLimits>,
    pub devices_opened: Option<u64>,
    pub instances_opened: Option<u64>,
    pub counter_scope: String,
}

/// The limits relevant to compute-buffer and device-allocation failures.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub struct BackendLimits {
    pub max_storage_buffers_per_shader_stage: u32,
    pub max_buffer_size: u64,
    pub max_storage_buffer_binding_size: u64,
}

impl State {
    /// The bounded part done, 0 to 1, when the stage has a known end.
    pub fn fraction(&self) -> Option<f64> {
        let total = self.total.filter(|t| *t > 0.0 && t.is_finite())?;
        Some((self.done / total).clamp(0.0, 1.0))
    }

    /// One line: `Step 1/2 Load · Increment 3 of 10`.
    pub fn line(&self) -> String {
        let mut parts = Vec::new();
        if let Some(s) = &self.step {
            parts.push(format!("Step {}/{} {}", s.index, s.count, s.name));
        }
        if !self.stage.is_empty() {
            let count = match self.total {
                Some(t) => format!(" {} of {}", number(self.done), number(t)),
                None if self.done > 0.0 => format!(" {}", number(self.done)),
                None => String::new(),
            };
            let unit = if self.unit.is_empty() || count.is_empty() { String::new() } else { format!(" {}", self.unit) };
            parts.push(format!("{}{count}{unit}", self.stage));
        }
        parts.join(" · ")
    }
}

/// A count as text: whole numbers without a fraction, others to 3 significant digits.
fn number(x: f64) -> String {
    if x.fract() == 0.0 && x.abs() < 1e15 { format!("{}", x as i64) } else { format!("{x:.3}") }
}

/// A command's progress, shared between the thread that runs it and those that show it.
#[derive(Debug, Default)]
pub struct Progress {
    state: Mutex<State>,
    cancelled: AtomicBool,
    backend: Mutex<Option<BackendDiagnostic>>,
    observations: AtomicU64,
}

impl Progress {
    pub fn new() -> Arc<Progress> {
        Arc::new(Progress::default())
    }

    /// The state now (the default if a reporter panicked while holding it).
    pub fn snapshot(&self) -> State {
        self.state.lock().map(|s| s.clone()).unwrap_or_default()
    }

    /// Reports observed, even when two reports have identical final values.
    pub fn observations(&self) -> u64 {
        self.observations.load(Ordering::Relaxed)
    }

    /// A new command keeps the host's cancellation token but no previous progress/device.
    pub fn begin_command(&self) {
        self.update(|s| *s = State::default());
        self.set_backend_diagnostic(None);
    }

    /// This reporter's backend context, independent of every other command's reporter.
    pub fn backend_diagnostic(&self) -> Option<BackendDiagnostic> {
        self.backend.lock().ok().and_then(|v| v.clone())
    }

    /// Replace setup metadata, or clear it at the beginning of another solve. This does not
    /// reset cancellation or progress and never runs in an iteration/cell loop.
    pub fn set_backend_diagnostic(&self, value: Option<BackendDiagnostic>) {
        if let Ok(mut context) = self.backend.lock() {
            *context = value;
        }
    }

    /// Ask the command to stop at its next check ([`check`]); from any thread.
    pub fn cancel(&self) {
        self.cancelled.store(true, Ordering::Relaxed);
    }

    /// Whether [`Progress::cancel`] was called.
    pub fn is_cancelled(&self) -> bool {
        self.cancelled.load(Ordering::Relaxed)
    }

    /// `Err(Cancelled)` once cancelled: for work on a thread that does not report (a thread
    /// pool's), which holds the caller's [`Progress`] from [`current`].
    pub fn check(&self) -> crate::Result<()> {
        if self.is_cancelled() { Err(crate::Error::Cancelled) } else { Ok(()) }
    }

    fn update(&self, f: impl FnOnce(&mut State)) {
        if let Ok(mut s) = self.state.lock() {
            f(&mut s);
            self.observations.fetch_add(1, Ordering::Relaxed);
        }
    }
}

thread_local! {
    static CURRENT: RefCell<Option<Arc<Progress>>> = const { RefCell::new(None) };
}

/// Run `f` with `progress` as this thread's reporter (none: reports are dropped); the previous
/// one comes back after, even if `f` panics.
pub fn scoped<R>(progress: Option<Arc<Progress>>, f: impl FnOnce() -> R) -> R {
    struct Restore(Option<Arc<Progress>>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0.take();
            CURRENT.with(|c| *c.borrow_mut() = previous);
        }
    }
    let previous = CURRENT.with(|c| c.replace(progress));
    let _restore = Restore(previous);
    f()
}

/// This thread's reporter, to install on another thread that reports for it.
pub fn current() -> Option<Arc<Progress>> {
    CURRENT.with(|c| c.borrow().clone())
}

fn with(f: impl FnOnce(&Progress)) {
    if let Some(p) = current() {
        f(&p);
    }
}

/// Record backend setup on this command's reporter, including on a cached-device reuse.
pub fn backend_diagnostic(context: BackendDiagnostic) {
    with(|p| p.set_backend_diagnostic(Some(context)));
}

/// Step `index` (from 1) of `count`, named `name`, starts.
pub fn step(name: &str, index: usize, count: usize) {
    with(|p| {
        p.set_backend_diagnostic(None);
        p.update(|s| {
            *s = State { step: Some(StepAt { name: name.to_string(), index, count }), ..State::default() };
        })
    });
}

/// A stage starts: its count is reset, its end and unit are `total` and `unit`.
pub fn stage(name: &str, total: Option<f64>, unit: &'static str) {
    with(|p| {
        p.update(|s| {
            s.stage.clear();
            s.stage.push_str(name);
            (s.done, s.total, s.unit) = (0.0, total, unit);
        })
    });
}

/// The current stage has `done` of its work done.
pub fn count(done: f64) {
    with(|p| p.update(|s| s.done = done));
}

/// `Err(Cancelled)` when this thread's command was cancelled; `Ok` with no [`Progress`]
/// installed.
pub fn check() -> crate::Result<()> {
    CURRENT.with(|c| c.borrow().as_ref().map_or(Ok(()), |p| p.check()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn backend_context_stays_with_its_reporter_when_commands_overlap() {
        let barrier = Arc::new(std::sync::Barrier::new(2));
        let reporters = [Progress::new(), Progress::new()];
        std::thread::scope(|scope| {
            for (index, reporter) in reporters.iter().enumerate() {
                let barrier = barrier.clone();
                scope.spawn(move || {
                    scoped(Some(reporter.clone()), || {
                        backend_diagnostic(BackendDiagnostic { adapter: Some(format!("adapter {index}")), ..Default::default() });
                        barrier.wait();
                        assert_eq!(reporter.backend_diagnostic().unwrap().adapter, Some(format!("adapter {index}")));
                    })
                });
            }
        });
        reporters[0].cancel();
        scoped(Some(reporters[0].clone()), || step("CPU after GPU", 2, 2));
        assert!(reporters[0].is_cancelled());
        assert!(reporters[0].backend_diagnostic().is_none());
        assert_eq!(reporters[1].backend_diagnostic().unwrap().adapter.as_deref(), Some("adapter 1"));
    }

    #[test]
    fn identical_final_progress_still_counts_as_new_observations() {
        let p = Progress::new();
        scoped(Some(p.clone()), || {
            stage("Iteration", Some(10.0), "increments");
            count(3.0);
            let before = p.observations();
            let state = p.snapshot();
            stage("Iteration", Some(10.0), "increments");
            count(3.0);
            assert_eq!(p.snapshot(), state);
            assert!(p.observations() > before);
        });
        p.cancel();
        p.begin_command();
        assert!(p.is_cancelled());
        assert_eq!(p.snapshot(), State::default());
    }

    #[test]
    fn reports_reach_the_installed_progress_only() {
        stage("ignored", None, "");
        let p = Progress::new();
        scoped(Some(p.clone()), || {
            step("Load", 1, 2);
            stage("Increment", Some(10.0), "increments");
            count(3.0);
        });
        let s = p.snapshot();
        assert_eq!(s.line(), "Step 1/2 Load · Increment 3 of 10 increments");
        assert_eq!(s.fraction(), Some(0.3));
        assert!(current().is_none());
        // Nothing installed: nothing reported, nothing panics.
        count(5.0);
        assert_eq!(p.snapshot().done, 3.0);
    }

    #[test]
    fn check_fails_once_cancelled() {
        assert_eq!(check(), Ok(()));
        let p = Progress::new();
        scoped(Some(p.clone()), || {
            assert_eq!(check(), Ok(()));
            p.cancel();
            assert_eq!(check(), Err(crate::Error::Cancelled));
        });
        // Not this thread's command any more.
        assert_eq!(check(), Ok(()));
    }
}
