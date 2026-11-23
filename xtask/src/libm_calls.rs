//! The libm check: std transcendental methods (`exp`, `ln`, `sin`, `powf`, …) and their `f64::`
//! path forms in the CPU result crates, a check of `cargo xtask ci` (#1038).
//!
//! `f64::exp` and friends call the operating system's maths library (Apple's libm, glibc, the
//! CRT), which is not correctly rounded and differs in the last bits between platforms. A result
//! computed with one is not the same bits on every platform, which AGENTS.md "Determinism"
//! forbids. Those calls go through `libm::` (a port of musl's, pure Rust, the workspace
//! dependency). `sqrt` and `abs` are correctly rounded by IEEE 754 and are not flagged.
//!
//! A call that provably feeds no result (a diagnostic, a log line) is marked on its line with
//! `// libm: not a result: <why>`. Every CPU result crate is converted, so the check has no
//! list of files still to do.

use std::path::Path;

use crate::rust_source;

/// The CPU result crates: a result they compute is the same bits on every platform.
const CPU_CRATES: &[&str] = &["decks", "engine", "explicit", "fem", "lbm", "mesh", "model", "solver"];

/// The per-line exemption for a call that feeds no result.
const MARKER: &str = "// libm: not a result:";

pub fn run(root: &Path) -> Result<(), String> {
    let dirs: Vec<String> = CPU_CRATES.iter().map(|c| format!("crates/{c}")).collect();
    let dir_refs: Vec<&str> = dirs.iter().map(String::as_str).collect();
    let (sources, mut errors) = rust_source::read_tree(root, &dir_refs)?;

    let mut total = 0usize;
    for source in &sources {
        let lines: Vec<&str> = source.text.lines().collect();
        let marked = |line: usize| {
            let text = lines.get(line.wrapping_sub(1)).copied().unwrap_or("");
            rust_source::marker_reason(text, MARKER)
        };
        for (line, method) in &source.calls {
            match marked(*line) {
                Some("") => errors.push(format!("{}:{line}: `{MARKER}` needs a reason", source.rel)),
                Some(_) => {}
                None => {
                    total += 1;
                    errors.push(format!("{}:{line}: `{method}` is a platform maths call: use libm:: (#1038)", source.rel));
                }
            }
        }
        let called: Vec<usize> = source.calls.iter().map(|(line, _)| *line).collect();
        for (n, text) in lines.iter().enumerate() {
            if rust_source::marker_reason(text, MARKER).is_some() && !called.contains(&(n + 1)) {
                errors.push(format!("{}:{}: `{MARKER}` on a line without a transcendental call; remove it", source.rel, n + 1));
            }
        }
    }

    for e in &errors {
        println!("{e}");
    }
    println!("libm: {total} platform maths call(s) in the CPU crates");
    if errors.is_empty() {
        Ok(())
    } else {
        Err(format!("libm: {} problem(s): a transcendental in a result crate goes through libm:: (#1038)", errors.len()))
    }
}
