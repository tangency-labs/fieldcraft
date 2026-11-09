//! Opt-in diagnostic wall times. No solver operands or scheduling depend on these clocks.
//! Records are inclusive (in particular, parallel solve spans must not be summed as wall time).
use std::io::Write;
use std::sync::OnceLock;
use std::time::Instant;

pub(crate) struct Span {
    name: &'static str,
    start: Option<Instant>,
}

impl Span {
    pub(crate) fn new(name: &'static str) -> Self {
        static ENABLED: OnceLock<bool> = OnceLock::new();
        let enabled = *ENABLED.get_or_init(|| std::env::var("FIELDCRAFT_SOLVER_TIMING").is_ok_and(|v| v == "1"));
        Self { name, start: enabled.then(Instant::now) }
    }
}

impl Drop for Span {
    fn drop(&mut self) {
        if let Some(start) = self.start {
            let seconds = start.elapsed().as_secs_f64();
            // Diagnostics must never turn an otherwise successful solve into an I/O error.
            let _ = writeln!(std::io::stderr().lock(), "FIELDCRAFT_TIMING\t{}\t{seconds:.9}", self.name);
        }
    }
}
