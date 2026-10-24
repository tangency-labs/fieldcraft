//! fieldcraft-cli commands
//! fieldcraft-cli run <script.json>
//! fieldcraft-cli snapshot <script.json> --out shot.png [--view V | --camera JSON] [--size WxH] [--section x|y|z[=V]] [--field F]
//! fieldcraft-cli features
//! fieldcraft-cli mcp [--connect ADDR | --headless] [--automation-read-root DIR] [--automation-write-root DIR]
//! fieldcraft-cli crashes [list [--json] | show <id> | replay <id> | dir | clear [<id>]]
//! fieldcraft-cli update --check | --install
//! fieldcraft-cli --version
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]

mod crashes;
mod snapshot;
mod update;

use std::path::Path;
use std::process::ExitCode;

use fieldcraft_engine::Session;
use serde_json::json;

const USAGE: &str = "\
usage:
  fieldcraft-cli commands              list every command with its parameters
  fieldcraft-cli commands --json       the same as JSON, each command's parameters as JSON Schema
  fieldcraft-cli run <script.json>     run a script (a JSON array of {\"command\", \"params\"})
                                       and print the model summary
  fieldcraft-cli snapshot <script.json> --out FILE.png [--view V | --camera JSON]
                      [--size WxH] [--section x|y|z[=VALUE]] [--field FIELD]
                                       run a script and draw the model to a PNG with no
                                       window: the viewport's GPU pipeline on a hardware or
                                       software GPU, else the CPU rasterizer (FIELDCRAFT_GPU=
                                       auto|hardware|software|cpu); V: iso (default), x, y, z,
                                       front, back, left, right, top, bottom; size 1280x800;
                                       --section hides the elements beyond a plane (at the
                                       mesh's middle, or VALUE); --field contours a nodal
                                       field of the solved step (results.contour);
                                       --camera '{\"yaw\":30,\"pitch\":20,\"target\":...,
                                       \"fit\":true,\"fov\":35}' aims as MCP screenshot's camera
                                       (default target: the model's centre of mass)
  fieldcraft-cli features              print docs/features.md (every capability, ready or
                                       coming soon) from the capability catalogue
  fieldcraft-cli mcp [options]         serve the Model Context Protocol on stdin/stdout
                                       (docs/mcp.md):
      --headless                       an in-process session (the default)
      --connect ADDR                   drive the app at ADDR (`fieldcraft --control PORT`);
                                       headless when it cannot be reached
      --automation-read-root DIR       files may be read only inside DIR
      --automation-write-root DIR      files may be written only inside DIR (with either root,
                                       the access a root left out is not granted)
  fieldcraft-cli crashes [list|show <id>|replay <id>|dir|clear [<id>]]
                                       local crash reports and solver diagnostics
                                       (docs/crash-reports.md): list them (--json), show one,
                                       replay a crash (does it happen again?),
                                       print their folder, delete them
  fieldcraft-cli update --check        is a newer FieldCraft published? (app.update.check; JSON)
  fieldcraft-cli update --install      download the latest release, verify its signature and
                                       checksum, and install it over this copy (docs/updates.md)
  fieldcraft-cli --version             print the version
";

/// Largest script accepted (bytes).
const MAX_SCRIPT_BYTES: u64 = 256 * 1024 * 1024;

fn main() -> ExitCode {
    // The lattice Boltzmann GPU backend, for `solve.fluid` with `backend: "gpu"` (scripts and MCP).
    fieldcraft_lbm_gpu::register();
    fieldcraft_mesh_gpu::register();
    // Pictures commands draw themselves (results.export.animation), with the viewport's pipeline.
    fieldcraft_mcp::register_drawer();
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.first().map(String::as_str) {
        Some("commands") if args.get(1).map(String::as_str) == Some("--json") => {
            let all: Vec<_> = fieldcraft_engine::commands().iter().map(fieldcraft_engine::CommandSpec::listing).collect();
            serde_json::to_string_pretty(&all).map(|text| println!("{text}")).map_err(|e| e.to_string())
        }
        Some("commands") => {
            for c in fieldcraft_engine::commands() {
                println!("{:<24} {:<18} {}", c.id, c.label, c.params);
            }
            Ok(())
        }
        Some("features") => {
            print!("{}", fieldcraft_engine::catalogue::features_page());
            Ok(())
        }
        Some("mcp") => mcp(args.get(1..).unwrap_or(&[])),
        Some("crashes") => crashes::run(args.get(1..).unwrap_or(&[])),
        Some("update") => update::run(args.get(1..).unwrap_or(&[])),
        Some("snapshot") => snapshot::run(args.get(1..).unwrap_or(&[])),
        Some("run") => args.get(1).ok_or_else(|| format!("run needs a script\n\n{USAGE}")).and_then(|p| run(p)),
        Some("-V" | "--version") => {
            println!("fieldcraft-cli {}", env!("CARGO_PKG_VERSION"));
            Ok(())
        }
        Some("-h" | "--help" | "help") | None => {
            print!("{USAGE}");
            Ok(())
        }
        Some(other) => Err(format!("unknown command `{other}`\n\n{USAGE}")),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

/// The `mcp` options: each given once, as `--flag value` or `--flag=value`.
fn mcp_options(args: &[String]) -> Result<fieldcraft_mcp::Options, String> {
    let mut connect = None;
    let mut headless = false;
    let mut read = None;
    let mut write = None;
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        let (flag, inline) = match arg.split_once('=') {
            Some((f, v)) if f.starts_with("--") => (f, Some(v.to_string())),
            _ => (arg.as_str(), None),
        };
        let mut value = || inline.clone().or_else(|| it.next().cloned()).filter(|v| !v.is_empty()).ok_or_else(|| format!("{flag} needs a value"));
        let slot = match flag {
            "--headless" if inline.is_none() => {
                if std::mem::replace(&mut headless, true) {
                    return Err("--headless is given twice".into());
                }
                continue;
            }
            "--connect" => &mut connect,
            "--automation-read-root" => &mut read,
            "--automation-write-root" => &mut write,
            other => return Err(format!("mcp: unknown option `{other}`\n\n{USAGE}")),
        };
        if slot.replace(value()?).is_some() {
            return Err(format!("{flag} is given twice"));
        }
    }
    if headless && connect.is_some() {
        return Err("--connect and --headless exclude each other".into());
    }
    let roots = match (&read, &write) {
        (None, None) => fieldcraft_io::roots::Roots::default(),
        _ => fieldcraft_io::roots::Roots::new(read.as_deref().map(Path::new), write.as_deref().map(Path::new)).map_err(|e| e.to_string())?,
    };
    Ok(fieldcraft_mcp::Options { connect, roots })
}

fn mcp(args: &[String]) -> Result<(), String> {
    let (backend, note) = fieldcraft_mcp::backend(mcp_options(args)?)?;
    if let Some(note) = note {
        eprintln!("fieldcraft-cli mcp: {note}");
    }
    let mut server = fieldcraft_mcp::Server::new(backend);
    server.serve(std::io::stdin().lock(), std::io::stdout().lock()).map_err(|e| format!("mcp: {e}"))
}

/// The JSON report of a picture: what drew it and what it shows (the MCP `screenshot` reply's
/// fields).
pub fn picture_report(r: &fieldcraft_gpu::Renderer, shown: &fieldcraft_engine::view::Shown) -> serde_json::Value {
    fieldcraft_mcp::picture_report(r, shown)
}

/// The session after running the script at `path`.
fn session_of(path: &str) -> Result<Session, String> {
    fieldcraft_io::recovery::use_automation_dir();
    let size = std::fs::metadata(path).map_err(|e| format!("{path}: {e}"))?.len();
    if size > MAX_SCRIPT_BYTES {
        return Err(format!("{path}: larger than {} MB", MAX_SCRIPT_BYTES / (1024 * 1024)));
    }
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let script: serde_json::Value = serde_json::from_str(&text).map_err(|e| format!("{path}: {e}"))?;
    let mut session = Session::default();
    session.solver_diagnostics = true;
    session.run_script(&script).map_err(|e| e.to_string())?;
    Ok(session)
}

fn run(path: &str) -> Result<(), String> {
    let mut session = session_of(path)?;
    let info = session.execute("model.info", &json!({})).map_err(|e| e.to_string())?;
    println!("{}", serde_json::to_string_pretty(&info).map_err(|e| e.to_string())?);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn opts(args: &[&str]) -> Result<fieldcraft_mcp::Options, String> {
        mcp_options(&args.iter().map(|s| s.to_string()).collect::<Vec<_>>())
    }

    #[test]
    fn mcp_options_are_checked() {
        let dir = std::env::temp_dir();
        let d = dir.to_string_lossy();
        assert!(opts(&[]).unwrap().connect.is_none());
        assert_eq!(opts(&["--connect", "127.0.0.1:7070"]).unwrap().connect.as_deref(), Some("127.0.0.1:7070"));
        assert_eq!(opts(&["--connect=127.0.0.1:7070"]).unwrap().connect.as_deref(), Some("127.0.0.1:7070"));
        assert!(opts(&["--automation-read-root", &d]).unwrap().roots.is_confined());
        assert!(opts(&[&format!("--automation-write-root={d}")]).unwrap().roots.is_confined());
        for bad in [
            vec!["--connect"],
            vec!["--connect", "a", "--connect", "b"],
            vec!["--headless", "--headless"],
            vec!["--headless", "--connect", "a"],
            vec!["--automation-read-root", "/no/such/folder"],
            vec!["--automation-read-root="],
            vec!["--nope"],
            vec!["--headless=yes"],
        ] {
            assert!(opts(&bad).is_err(), "{bad:?}");
        }
    }
}
