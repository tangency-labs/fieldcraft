//! Embeds the toolbar icons: `assets/icons/*.svg` become `ICONS`, a sorted `(id, svg)` table, and
//! `KIT_MASTERS`, the same artwork as the masters of the UI kit's icon set (`src/icons.rs`).

use std::fmt::Write;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let manifest = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR")?);
    let dir = manifest.join("../../assets/icons");
    println!("cargo:rerun-if-changed={}", dir.display());
    let mut ids: Vec<String> =
        std::fs::read_dir(&dir)?.flatten().filter_map(|e| e.file_name().to_str().and_then(|n| n.strip_suffix(".svg")).map(str::to_string)).collect();
    ids.sort();
    let mut out = String::from("/// Every toolbar icon by catalogue id, sorted by id.\npub static ICONS: &[(&str, &str)] = &[\n");
    let mut kit =
        String::from("/// The icons as the kit's masters: (id, family, SVG), sorted by id.\npub static KIT_MASTERS: &[(&str, &str, &str)] = &[\n");
    for id in &ids {
        let path = dir.join(format!("{id}.svg"));
        println!("cargo:rerun-if-changed={}", path.display());
        let _ = writeln!(out, "    ({id:?}, include_str!({:?})),", path.display().to_string());
        let master = KIT_ATTRS.iter().fold(std::fs::read_to_string(&path)?, |svg, attr| scaled(&svg, attr, KIT_SCALE));
        let _ = writeln!(kit, "    ({id:?}, \"\", {master:?}),");
    }
    out.push_str("];\n");
    kit.push_str("];\n");
    out.push_str(&kit);
    std::fs::write(PathBuf::from(std::env::var("OUT_DIR")?).join("icons.rs"), out)?;
    Ok(())
}

/// The kit draws a master on a 64-unit square; FieldCraft's artwork is drawn on 24 (`cargo xtask
/// icons`), so its coordinates and stroke weights are scaled by 64 / 24 for the kit.
const KIT_SCALE: f64 = 64.0 / 24.0;

/// The attributes that hold lengths in the artwork (paths, stroke weights and dashes).
const KIT_ATTRS: [&str; 3] = ["d", "stroke-width", "stroke-dasharray"];

/// Scaled coordinates keep three decimals.
const ROUND: f64 = 1000.0; // unit-factors: rounding to three decimals, not a unit conversion

/// `svg` with every number in the values of attribute `attr` multiplied by `factor`; path
/// commands, separators and other attributes are kept as they are.
fn scaled(svg: &str, attr: &str, factor: f64) -> String {
    let key = format!(" {attr}=\"");
    let mut out = String::with_capacity(svg.len());
    let mut rest = svg;
    while let Some(at) = rest.find(&key) {
        let (head, tail) = rest.split_at(at + key.len());
        out.push_str(head);
        let end = tail.find('"').unwrap_or(tail.len());
        let (value, after) = tail.split_at(end);
        let mut number = String::new();
        for c in value.chars().chain(std::iter::once(' ')) {
            // A minus sign starts a number, except an exponent's (`1e-3`).
            let sign = c == '-' && (number.is_empty() || number.ends_with(['e', 'E']));
            if c.is_ascii_digit() || c == '.' || sign || ((c == 'e' || c == 'E') && !number.is_empty()) {
                number.push(c);
                continue;
            }
            if !number.is_empty() {
                match number.parse::<f64>() {
                    Ok(x) => {
                        let _ = write!(out, "{}", (x * factor * ROUND).round() / ROUND);
                    }
                    Err(_) => out.push_str(&number),
                }
                number.clear();
            }
            out.push(c);
        }
        // The space chained after the value.
        out.pop();
        rest = after;
    }
    out.push_str(rest);
    out
}
