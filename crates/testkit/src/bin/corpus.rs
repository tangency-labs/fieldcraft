//! `corpus [--part NAME]... [--bless NAME|all]... [--jobs N] [--no-write]`: see
//! [`fieldcraft_testkit::runner`]; run it as `cargo xtask corpus`.
//!
//! `corpus --part NAME --design PATH` writes that part's `.solidcraft` design to `PATH` instead
//! (its recipe run, or its shape generated), for opening the part in the app or over MCP.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]

use std::path::Path;
use std::process::ExitCode;

/// Most parts meshed at once.
const MAX_JOBS: usize = 64;
const USAGE: &str = "usage: corpus [--part NAME]... [--bless NAME|all]... [--jobs N] [--no-write] | corpus --part NAME --design PATH";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut options = fieldcraft_testkit::runner::Options { write: true, ..Default::default() };
    let mut design = None;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let value = if matches!(a.as_str(), "--part" | "--bless" | "--design" | "--jobs") { it.next() } else { None };
        match (a.as_str(), value) {
            ("--part", Some(name)) => options.parts.push(name.clone()),
            ("--bless", Some(name)) => options.bless.push(name.clone()),
            ("--design", Some(path)) => design = Some(path.as_str()),
            ("--jobs", Some(n)) => match n.parse::<usize>() {
                Ok(n) if (1..=MAX_JOBS).contains(&n) => options.jobs = n,
                _ => {
                    eprintln!("error: --jobs takes a whole number from 1 to {MAX_JOBS}\n{USAGE}");
                    return ExitCode::FAILURE;
                }
            },
            ("--no-write", None) => options.write = false,
            (other, _) => {
                eprintln!("error: unknown argument `{other}` or a missing value\n{USAGE}");
                return ExitCode::FAILURE;
            }
        }
    }
    let result = match (design, options.parts.as_slice()) {
        (Some(path), [name]) => write_design(name, Path::new(path)),
        (Some(_), _) => Err(format!("--design needs one --part\n{USAGE}")),
        (None, _) => fieldcraft_testkit::runner::run(&fieldcraft_testkit::corpus::workspace_root(), &options),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

/// Write the design of the corpus part `name` to `path`.
fn write_design(name: &str, path: &Path) -> Result<(), String> {
    let parts = fieldcraft_testkit::corpus::parts().map_err(|e| e.to_string())?;
    let part = parts.iter().find(|p| p.name == name).ok_or_else(|| format!("no corpus part named `{name}`"))?;
    let bytes = fieldcraft_testkit::corpus::design_bytes(part).map_err(|e| e.to_string())?;
    std::fs::write(path, bytes).map_err(|e| format!("{}: {e}", path.display()))?;
    eprintln!("corpus: wrote {} ({name})", path.display());
    Ok(())
}
