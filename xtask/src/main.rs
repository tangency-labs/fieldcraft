//! Workspace tooling: `cargo xtask <command>`.
//!
//! Pure Rust (std + serde_json). External tools (`cargo`) are invoked through
//! `std::process::Command`. Adapted from SolidCraft's xtask (MIT OR Apache-2.0), see
//! ATTRIBUTION.md.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]

mod affected;
mod app_icon;
mod assets;
mod book;
mod icons;
mod layers;
mod licences;
mod skill;
mod test;
mod unit_factors;

use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};

const USAGE: &str = "\
usage: cargo xtask <command>

commands:
  layers          enforce the crate dependency layering (docs/architecture.md)
  assets          check that every non-code asset file is attributed in ATTRIBUTION.md
  icons [ID…]     draw the toolbar icons from their recipes (assets/icons/src/*.json) with
                  VectorCraft (VECTORCRAFT_CLI, else vectorcraft-cli on the PATH): writes each
                  icon's .vectorcraft source and .svg export
  app-icon        draw the application icon family (assets/app-icon) with VectorCraft, then
                  write the iconset, icns, ico and Icon Composer document (macOS: sips,
                  iconutil, magick; actool for Assets.car when available)
  book [--no-build]
                  build the documentation site (docs/book.toml, docs/SUMMARY.md) with mdbook
                  into target/book; links that leave docs/ point at the repository on GitHub;
                  --no-build stages target/book-src only
  licences        write THIRD-PARTY-LICENSES.txt (the licences of every shipped crate; the
                  release workflow runs it before packaging, the file is not committed)
  test [--heavy] [CARGO ARGS] [-- TEST ARGS]
                  build the tests with the ci profile, then run every test binary at once;
                  --heavy runs only the ignored tests instead (heavy verification and corpus
                  tests, and the benchmarks) and the doc tests alongside: the CI heavy job
  corpus [--part NAME]... [--bless NAME|all]... [--jobs N] [--no-write]
                  build every meshing test corpus part (tests/corpus, analytic shapes), tessellate
                  and mesh it, check each part against its fixture.json, and write the results
                  to docs/meshing/corpus.md (release build);
                  --bless rewrites the named parts' fixtures (docs/meshing/corpus.md says when);
                  --jobs N meshes N parts at once (the table's times are then not a part's alone);
                  --part NAME --design PATH writes that part's .solidcraft design instead
  skill-results [--part NAME] [--no-write]
                  mesh every corpus part and sample by following the meshing skill
                  (skills/fieldcraft-meshing) and write docs/meshing/skill-results.md and each
                  part's section picture (docs/meshing/skill/*.png; release build)
  oracle [--case NAME]... [--bless NAME|all]... [--quick] [--nightly] [--no-write]
                  run the verification and validation cases (tests/oracle, docs/verification.md)
                  against their references and fixtures; write docs/oracle.md after a full run;
                  --bless freezes a passing case's results in its fixture.json; --quick skips
                  cases marked slow; --nightly also runs the cases marked nightly, which are
                  otherwise checked by their fixtures (release build)
  affected --base REV [--head REV] [--github-output]
                  the workspace crates changed between REV and HEAD-REV (the working tree) and every crate that depends on them,
                  and the flags CI plans its jobs with: everything (manifest, lock file,
                  toolchain, xtask or workflows changed), solver (the heavy tests or the oracle
                  are affected), wasm (a crate the wasm check builds is affected)
  unit-factors    check for literal unit-conversion factors (1e-3, 1e6, 273.15, 25.4, …) written
                  outside the unit table, crates/model/src/units (docs/units.md); the geometry
                  frame is reported apart until its change lands
  wasm            cargo check every library crate and the web app for wasm32-unknown-unknown
  ci [--affected REV]
                  fmt --check, assets, layers, unit-factors, clippy -D warnings, then the tests
                  (not the heavy ones) with the wasm check alongside (stops at the first failing
                  step); --affected REV lints and tests only the crates `affected --base REV`
                  lists, and checks wasm only when it flags wasm
";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let rest: Vec<&str> = args.iter().skip(1).map(String::as_str).collect();
    let result = match args.first().map(String::as_str) {
        Some("layers") => cmd_layers(),
        Some("assets") => assets::run(&root()),
        Some("icons") => icons::run(&root(), &rest),
        Some("app-icon") => app_icon::run(&root()),
        Some("book") => book::run(&root(), !rest.contains(&"--no-build")),
        Some("licences") => licences::run(&root()),
        Some("test") => {
            let (heavy, rest) = match rest.split_first() {
                Some((&"--heavy", after)) => (true, after),
                _ => (false, &rest[..]),
            };
            let (cargo_args, test_args) = rest.iter().position(|a| *a == "--").map_or((rest, &[][..]), |i| (&rest[..i], &rest[i + 1..]));
            if heavy {
                let ignored: Vec<&str> = ["--ignored"].into_iter().chain(test_args.iter().copied()).collect();
                test::run(&[], cargo_args, &ignored, vec![test::doc_tests(cargo_args, test_args)])
            } else {
                test::run(&[], cargo_args, test_args, vec![])
            }
        }
        Some("affected") => affected::run(&rest),
        Some("ci") => match rest.as_slice() {
            [] => cmd_ci(),
            ["--affected", base] => cmd_ci_affected(base),
            _ => Err(format!("ci: expected no arguments or `--affected REV`\n\n{USAGE}")),
        },
        Some("corpus") => {
            // A binary of the testkit, so xtask itself stays free of SolidCraft.
            let mut c = cargo();
            c.args(["run", "--release", "-p", "fieldcraft-testkit", "--bin", "corpus", "--"]).args(&rest);
            run(c, "corpus")
        }
        Some("skill-results") => skill::run(&root(), &rest),
        Some("oracle") => {
            // A binary of the testkit, which has the engine and SolidCraft scripts for geometry.
            let mut c = cargo();
            c.args(["run", "--release", "-p", affected::ORACLE_PACKAGE, "--bin", "oracle", "--"]).args(&rest);
            run(c, "oracle")
        }
        Some("unit-factors") => unit_factors::run(&root()),
        Some("wasm") => cmd_wasm(),
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

/// Workspace root (parent of the xtask crate).
pub fn root() -> PathBuf {
    let here = Path::new(env!("CARGO_MANIFEST_DIR"));
    here.parent().unwrap_or(here).to_path_buf()
}

pub fn cargo() -> Command {
    let mut c = Command::new(std::env::var("CARGO").unwrap_or_else(|_| "cargo".into()));
    c.current_dir(root());
    c
}

fn run(mut cmd: Command, what: &str) -> Result<(), String> {
    eprintln!("$ {what}");
    let status = cmd.status().map_err(|e| format!("{what}: failed to spawn: {e}"))?;
    if status.success() { Ok(()) } else { Err(format!("{what}: exited with {status}")) }
}

pub fn metadata() -> Result<serde_json::Value, String> {
    let out = cargo().args(["metadata", "--format-version", "1", "--no-deps"]).output().map_err(|e| format!("cargo metadata: {e}"))?;
    if !out.status.success() {
        return Err(format!("cargo metadata failed:\n{}", String::from_utf8_lossy(&out.stderr)));
    }
    serde_json::from_slice(&out.stdout).map_err(|e| format!("cargo metadata: bad JSON: {e}"))
}

fn cmd_layers() -> Result<(), String> {
    let crates = layers::from_metadata(&metadata()?)?;
    println!("Dependency layering (docs/architecture.md)\n");
    println!("{:<28} {:<14} workspace deps", "crate", "layer");
    for c in &crates {
        let ws: Vec<String> = c
            .deps
            .iter()
            .filter(|d| d.workspace)
            .map(|d| {
                let k = match d.kind {
                    layers::DepKind::Normal => "",
                    layers::DepKind::Dev => " (dev)",
                    layers::DepKind::Build => " (build)",
                };
                format!("{}{k}", layers::short_name(&d.name))
            })
            .collect();
        println!("{:<28} {:<14} {}", c.name, layers::describe(layers::classify(&c.name)), ws.join(", "));
    }
    let violations = layers::check(&crates);
    println!();
    if violations.is_empty() {
        println!("OK: {} crates, no layering violations.", crates.len());
        Ok(())
    } else {
        println!("{} violation(s):", violations.len());
        for v in &violations {
            println!("  - {v}");
        }
        Err(format!("{} layering violation(s)", violations.len()))
    }
}

/// Never break wasm: every layered crate (everything below the apps) and the web app must build
/// for the browser.
fn cmd_wasm() -> Result<(), String> {
    let (c, n) = wasm_check()?;
    run(c, &format!("cargo check --target wasm32-unknown-unknown ({n} crates)"))?;
    println!("wasm: {n} crates build for wasm32-unknown-unknown");
    Ok(())
}

/// The app that runs in the browser (apps/fieldcraft-web), checked with the layered crates.
const WEB_APP: &str = "fieldcraft-web";

/// Whether the wasm check builds `pkg`: a layered crate or the web app.
pub fn wasm_crate(pkg: &str) -> bool {
    pkg == WEB_APP || matches!(layers::classify(pkg), Some(layers::Class::Layer(_)))
}

/// The wasm check and how many crates it covers. It uses the ci profile so the host build
/// scripts and proc macros built for clippy and the tests are reused.
fn wasm_check() -> Result<(Command, usize), String> {
    let crates = layers::from_metadata(&metadata()?)?;
    let set: Vec<String> = crates.iter().filter(|c| wasm_crate(&c.name)).map(|c| c.name.clone()).collect();
    if !set.iter().any(|c| c == WEB_APP) {
        return Err(format!("wasm: the web app `{WEB_APP}` is not in the workspace"));
    }
    let mut c = cargo();
    c.args(["check", "--target", "wasm32-unknown-unknown", "--profile", "ci"]);
    for p in &set {
        c.args(["-p", p]);
    }
    Ok((c, set.len()))
}

/// The quick gates first, so they fail fast; the wasm check runs while the test binaries do,
/// which leaves cargo's build lock free. The heavy tests and the doc tests are the CI heavy job's
/// (`cargo xtask test --heavy`).
fn cmd_ci() -> Result<(), String> {
    quick_gates()?;
    let mut c = cargo();
    c.args(["clippy", "--workspace", "--all-targets", "--profile", "ci", "--", "-D", "warnings"]);
    run(c, "cargo clippy -D warnings")?;
    let (wasm, n) = wasm_check()?;
    test::run(&[], &[], &[], vec![(format!("wasm32-unknown-unknown check ({n} crates)"), wasm)])?;
    eprintln!("ci: all gates passed");
    Ok(())
}

/// The gates that look at the whole tree and take seconds: they run on every change.
fn quick_gates() -> Result<(), String> {
    let mut c = cargo();
    c.args(["fmt", "--all", "--", "--check"]);
    run(c, "cargo fmt --check")?;
    assets::run(&root())?;
    cmd_layers()?;
    unit_factors::run(&root())
}

/// `ci --affected REV`: the quick gates, then clippy and the tests of the affected crates only, and
/// the wasm check when an affected crate is built for wasm. A change that runs everything is the
/// full gate.
fn cmd_ci_affected(base: &str) -> Result<(), String> {
    let (plan, files) = affected::plan_since(base, None)?;
    affected::print(&plan, base, files);
    if plan.everything {
        return cmd_ci();
    }
    quick_gates()?;
    if plan.crates.is_empty() {
        eprintln!("ci: no crate affected; the quick gates passed");
        return Ok(());
    }
    let packages: Vec<&str> = plan.crates.iter().map(String::as_str).collect();
    let mut c = cargo();
    c.args(["clippy", "--all-targets", "--profile", "ci"]);
    for p in &packages {
        c.args(["-p", p]);
    }
    c.args(["--", "-D", "warnings"]);
    run(c, &format!("cargo clippy -D warnings ({} crates)", packages.len()))?;
    let side = if plan.wasm {
        let (wasm, n) = wasm_check()?;
        vec![(format!("wasm32-unknown-unknown check ({n} crates)"), wasm)]
    } else {
        vec![]
    };
    test::run(&packages, &[], &[], side)?;
    eprintln!("ci: all gates passed for {} affected crate(s)", packages.len());
    Ok(())
}
