//! `cargo xtask unit-factors`: report literal unit-conversion factors (1e-3, 1e6, 273.15, 25.4, …)
//! written outside the unit table (`crates/model/src/units/`).
//!
//! Issue #579 routes every conversion through the table. Until the later steps have removed the
//! existing factors this is a report and never fails; once the count is zero it becomes a CI check.
//! Many hits are tolerances or algorithm constants rather than conversions: the report lists them
//! all so each can be classified.

use std::path::Path;

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
const SKIP_DIRS: &[&str] = &["target", ".git", "plan", "dist", "xtask", "tests", "node_modules"];

pub fn run(root: &Path) -> Result<(), String> {
    let mut files = Vec::new();
    for dir in ["crates", "apps"] {
        walk(&root.join(dir), &mut files);
    }
    files.sort();
    let mut hits = Vec::new();
    for path in &files {
        let rel = path.strip_prefix(root).unwrap_or(path).to_string_lossy().replace('\\', "/");
        if rel.starts_with(TABLE) {
            continue;
        }
        let Ok(text) = std::fs::read_to_string(path) else { continue };
        for (n, line) in text.lines().enumerate() {
            // Unit tests sit at the end of a file behind `#[cfg(test)]`; their values are data.
            if line.trim_start().starts_with("#[cfg(test)]") {
                break;
            }
            let code = line.split("//").next().unwrap_or("");
            for value in literals(code) {
                if let Some((_, what)) = FACTORS.iter().find(|(f, _)| *f == value) {
                    hits.push(format!("{rel}:{}: {value:e} ({what}): {}", n + 1, line.trim()));
                }
            }
        }
    }
    for h in &hits {
        println!("{h}");
    }
    println!("unit-factors: {} literal factor(s) outside {TABLE} (report only; #579 drives this to zero)", hits.len());
    Ok(())
}

fn walk(dir: &Path, out: &mut Vec<std::path::PathBuf>) {
    let Ok(rd) = std::fs::read_dir(dir) else { return };
    for e in rd.flatten() {
        let p = e.path();
        let name = e.file_name().to_string_lossy().to_string();
        if p.is_dir() {
            if !SKIP_DIRS.contains(&name.as_str()) && !name.starts_with('.') {
                walk(&p, out);
            }
        } else if p.extension().is_some_and(|x| x == "rs") {
            out.push(p);
        }
    }
}

/// The floating-point literals on a line of code, parsed (`1_000.0`, `1e-3`,
/// `273.15_f64`). Digits inside identifiers (`vec3`) are skipped.
fn literals(code: &str) -> Vec<f64> {
    let chars: Vec<char> = code.chars().collect();
    let mut out = Vec::new();
    let mut i = 0;
    while let Some(&c) = chars.get(i) {
        let prev_ident = i > 0 && chars.get(i - 1).is_some_and(|p| p.is_alphanumeric() || *p == '_' || *p == '.');
        if !c.is_ascii_digit() || prev_ident {
            i += 1;
            continue;
        }
        let start = i;
        while let Some(&d) = chars.get(i) {
            let exp_sign = matches!(d, '+' | '-') && matches!(chars.get(i.wrapping_sub(1)), Some('e' | 'E'));
            if d.is_ascii_digit() || d == '_' || d == '.' || d == 'e' || d == 'E' || exp_sign {
                // `0..n` is a range, `1.max(x)` a method call: stop at a dot not followed by a digit.
                if d == '.' && !chars.get(i + 1).is_some_and(char::is_ascii_digit) {
                    break;
                }
                i += 1;
            } else {
                break;
            }
        }
        let text: String = chars.get(start..i).unwrap_or_default().iter().filter(|c| **c != '_').collect();
        let suffix_start = i;
        // Skip a type suffix (`f64`, `_f32`, `u8`).
        while chars.get(i).is_some_and(|c| c.is_alphanumeric() || *c == '_') {
            i += 1;
        }
        let suffix: String = chars.get(suffix_start..i).unwrap_or_default().iter().collect();
        // Conversions are floating-point arithmetic; integer literals are counts and caps.
        let float = text.contains(['.', 'e', 'E']) || suffix.trim_start_matches('_').starts_with('f');
        if float && let Ok(v) = text.parse::<f64>() {
            out.push(v);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::literals;

    #[test]
    fn literals_on_a_line() {
        assert_eq!(literals("let x = 1e-3 * y + 273.15_f64;"), vec![1e-3, 273.15]);
        assert_eq!(literals("for i in 0..1000 { vec3(1_000.0) }"), vec![1000.0]);
        assert_eq!(literals("a.max(1.0e6) - 2.max(b)"), vec![1e6]);
        assert_eq!(literals("t.0 + x1 + 25.4f32 + 3f64"), vec![25.4, 3.0]);
    }
}
