//! `cargo xtask unit-factors`: literal unit-conversion factors (1e-3, 1e6, 273.15, 25.4, …) written
//! outside the unit table (`crates/model/src/units/`), a check of `cargo xtask ci` (docs/units.md).
//!
//! Every conversion goes through the table (#579), so a float literal equal to a common factor is
//! taken for a hand-written conversion. A number that only looks like one, a tolerance or a cap,
//! is a named constant whose line carries the marker `// unit-factors: <why it is not a
//! conversion>`. The sources are read by [`crate::rust_source`]: comments, strings (help texts,
//! WGSL) and test code are not code that converts.
//!
//! One part of the tree is not counted, listed here with its reason: the frozen migration step
//! ([`FROZEN`]).

use std::collections::BTreeSet;
use std::path::Path;

use crate::rust_source;

/// The factors flagged, with what each usually converts.
const FACTORS: &[(f64, &str)] = &[
    (1e-3, "milli (mm↔m, ms↔s, g↔kg)"),
    (1e3, "kilo (t↔kg, m↔mm)"),
    (1e6, "mega (MPa↔Pa, mm²↔m²)"),
    (1e9, "giga (GPa↔Pa, mm³↔m³)"),
    (1e-6, "micro (Pa↔MPa)"),
    (1e-9, "nano (Pa↔GPa, m³↔mm³)"),
    (1e12, "t/mm³↔kg/m³"),
    (1e-12, "kg/m³↔t/mm³"),
    (273.15, "°C↔K"),
    (459.67, "°F↔°R"),
    (25.4, "in↔mm"),
    (0.0254, "in↔m"),
    (9.80665, "standard gravity"),
    (9.81, "gravity"),
    (9810.0, "gravity in mm/s²"),
];

/// The one place factors belong.
const TABLE: &str = "crates/model/src/units/";

/// Files that keep their literal factors, and why. Nothing else is exempt as a whole.
const FROZEN: &[(&str, &str)] = &[(
    "crates/io/src/migrate/v2_to_v3.rs",
    "the format 2 → 3 migration step is frozen: it converts old files with its own copy of the \
     factors, so a later change to the table cannot change how an old file migrates",
)];

/// The per-line exemption: `// unit-factors: <reason>` on the line of the literal.
const MARKER: &str = "// unit-factors:";

pub fn run(root: &Path) -> Result<(), String> {
    let (sources, mut errors) = rust_source::read_tree(root, &["crates", "apps"])?;

    let mut checked = Vec::new();
    let mut frozen = 0usize;
    for source in &sources {
        let rel = source.rel.as_str();
        if rel.starts_with(TABLE) {
            continue;
        }
        if FROZEN.iter().any(|(f, _)| *f == rel) {
            frozen += source.literals.iter().filter(|(_, v)| factor(*v).is_some()).count();
            continue;
        }
        let lines: Vec<&str> = source.text.lines().collect();
        let mut flagged = BTreeSet::new();
        for (line, value) in &source.literals {
            let Some(what) = factor(*value) else { continue };
            flagged.insert(*line);
            let text = lines.get(line.wrapping_sub(1)).copied().unwrap_or("");
            match rust_source::marker_reason(text, MARKER) {
                Some("") => errors.push(format!("{rel}:{line}: `{MARKER}` needs a reason")),
                Some(_) => {}
                None => checked.push(format!("{rel}:{line}: {value:e} ({what}): {}", text.trim())),
            }
        }
        for (n, text) in lines.iter().enumerate() {
            if rust_source::marker_reason(text, MARKER).is_some() && !flagged.contains(&(n + 1)) {
                errors.push(format!("{rel}:{}: `{MARKER}` on a line without a factor; remove it", n + 1));
            }
        }
    }

    for h in &checked {
        println!("{h}");
    }
    for e in &errors {
        println!("{e}");
    }
    for (file, why) in FROZEN {
        println!("unit-factors: {file} keeps its factors: {why}");
    }
    println!("unit-factors: {} literal factor(s) outside {TABLE}; {frozen} in the frozen step", checked.len());
    if checked.is_empty() && errors.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "unit-factors: {} literal factor(s) and {} marker error(s): convert through {TABLE}, or name the constant and \
             mark its line `{MARKER} <why it is not a conversion>` (docs/units.md)",
            checked.len(),
            errors.len()
        ))
    }
}

fn factor(value: f64) -> Option<&'static str> {
    FACTORS.iter().find(|(f, _)| *f == value).map(|(_, what)| *what)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn markers() {
        let marker = |line| rust_source::marker_reason(line, MARKER);
        assert_eq!(marker("const X: f64 = 1e-9; // unit-factors: relative to the size"), Some("relative to the size"));
        assert_eq!(marker("let x = 1e-9; // unit-factors:"), Some(""));
        assert_eq!(marker("let x = 1e-9;"), None);
    }
}
