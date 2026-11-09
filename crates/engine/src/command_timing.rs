//! Opt-in complete solve-command wall time, including in-session result storage.
use std::io::Write;
use std::sync::OnceLock;
use std::time::Instant;

pub(crate) struct Span {
    command: &'static str,
    start: Instant,
}

impl Span {
    pub(crate) fn for_command(id: &str) -> Option<Self> {
        let command = match id {
            "solve.static" => "solve.static",
            "solve.frequency" => "solve.frequency",
            "solve.harmonic" => "solve.harmonic",
            _ => return None,
        };
        static ENABLED: OnceLock<bool> = OnceLock::new();
        let enabled = *ENABLED.get_or_init(|| std::env::var("FIELDCRAFT_SOLVER_TIMING").is_ok_and(|v| v == "1"));
        enabled.then(|| Self { command, start: Instant::now() })
    }
}

impl Drop for Span {
    fn drop(&mut self) {
        let seconds = self.start.elapsed().as_secs_f64();
        let _ = writeln!(std::io::stderr().lock(), "FIELDCRAFT_COMMAND_TIMING\t{}\t{seconds:.9}", self.command);
    }
}
