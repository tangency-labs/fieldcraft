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
//! `// libm: not a result: <why>`. A file that still has calls to convert is on [`KNOWN`] with
//! its count: the count may only go down (the list is a ratchet), and a file not on the list
//! must have none. The list empties as the crates are converted.

use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use crate::rust_source;

/// The CPU result crates: a result they compute is the same bits on every platform.
const CPU_CRATES: &[&str] = &["decks", "engine", "explicit", "fem", "lbm", "mesh", "model", "solver"];

/// The per-line exemption for a call that feeds no result.
const MARKER: &str = "// libm: not a result:";

/// Files with std calls still to convert, and how many. Removed as each crate converts.
const KNOWN: &[(&str, usize)] = &[
    ("crates/engine/src/cmd/render.rs", 1),
    ("crates/engine/src/path.rs", 2),
    ("crates/engine/src/sample.rs", 6),
    ("crates/engine/src/view.rs", 3),
    ("crates/fem/src/constitutive.rs", 3),
    ("crates/fem/src/contact.rs", 1),
    ("crates/fem/src/elastoplastic.rs", 5),
    ("crates/fem/src/elastoplastic/logarithmic.rs", 3),
    ("crates/fem/src/material.rs", 3),
    ("crates/fem/src/plasticity.rs", 6),
    ("crates/fem/src/uniform.rs", 2),
    ("crates/model/src/contact.rs", 1),
    ("crates/model/src/materials.rs", 5),
    ("crates/model/src/pvt.rs", 3),
    ("crates/model/src/sections.rs", 1),
    ("crates/model/src/units/quantity.rs", 1),
];

pub fn run(root: &Path) -> Result<(), String> {
    let dirs: Vec<String> = CPU_CRATES.iter().map(|c| format!("crates/{c}")).collect();
    let dir_refs: Vec<&str> = dirs.iter().map(String::as_str).collect();
    let (sources, mut errors) = rust_source::read_tree(root, &dir_refs)?;

    // The calls per file that carry no exemption: (line, method).
    let mut calls: BTreeMap<&str, Vec<(usize, &str)>> = BTreeMap::new();
    for source in &sources {
        let lines: Vec<&str> = source.text.lines().collect();
        let mut exempt = BTreeSet::new();
        for (line, _) in &source.calls {
            let text = lines.get(line.wrapping_sub(1)).copied().unwrap_or("");
            match rust_source::marker_reason(text, MARKER) {
                Some("") => errors.push(format!("{}:{line}: `{MARKER}` needs a reason", source.rel)),
                Some(_) => {
                    exempt.insert(*line);
                }
                None => {}
            }
        }
        for (line, method) in &source.calls {
            if !exempt.contains(line) {
                calls.entry(source.rel.as_str()).or_default().push((*line, method.as_str()));
            }
        }
        let called: BTreeSet<usize> = source.calls.iter().map(|(line, _)| *line).collect();
        for (n, text) in lines.iter().enumerate() {
            if rust_source::marker_reason(text, MARKER).is_some() && !called.contains(&(n + 1)) {
                errors.push(format!("{}:{}: `{MARKER}` on a line without a transcendental call; remove it", source.rel, n + 1));
            }
        }
    }

    let known: BTreeMap<&str, usize> = KNOWN.iter().copied().collect();
    let files: BTreeSet<&str> = calls.keys().copied().chain(known.keys().copied()).collect();
    let mut total = 0usize;
    for file in files {
        let found = calls.get(file).map_or(&[][..], Vec::as_slice);
        let expected = known.get(file).copied().unwrap_or(0);
        total += found.len();
        if found.len() > expected {
            let detail: Vec<String> = found.iter().map(|(line, m)| format!("{line} `{m}`")).collect();
            errors.push(format!(
                "{file}: {} std transcendental call(s), the known list allows {expected}: {} (use libm::, or mark `{MARKER} <why>`)",
                found.len(),
                detail.join(", ")
            ));
        } else if found.len() < expected {
            errors.push(format!("{file}: {} call(s) left, the known list says {expected}: lower its entry", found.len()));
        }
    }

    for e in &errors {
        println!("{e}");
    }
    println!("libm: {total} std transcendental call(s) in the CPU crates all on the known list");
    if errors.is_empty() {
        Ok(())
    } else {
        Err(format!("libm: {} problem(s): a transcendental in a result crate goes through libm:: (#1038)", errors.len()))
    }
}
