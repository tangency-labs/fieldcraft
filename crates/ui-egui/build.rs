//! Embeds FieldCraft's icons (`assets/icons/*.svg`) as masters of the UI kit's icon set,
//! `MASTERS`, sorted by id (`src/icons.rs`): in SolidCraft's colour roles and on the kit's square.

use std::fmt::Write;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let manifest = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR")?);
    let dir = manifest.join("../../assets/icons");
    println!("cargo:rerun-if-changed={}", dir.display());
    let mut ids: Vec<String> =
        std::fs::read_dir(&dir)?.flatten().filter_map(|e| e.file_name().to_str().and_then(|n| n.strip_suffix(".svg")).map(str::to_string)).collect();
    ids.sort();
    let mut out = String::from(
        "/// FieldCraft's icons as the kit's masters: (id, family, SVG), sorted by id.\npub static MASTERS: &[(&str, &str, &str)] = &[\n",
    );
    for id in &ids {
        let path = dir.join(format!("{id}.svg"));
        println!("cargo:rerun-if-changed={}", path.display());
        let svg =
            ROLES.iter().fold(std::fs::read_to_string(&path)?, |svg, (ours, theirs)| svg.replace(&format!("\"{ours}\""), &format!("\"{theirs}\"")));
        let master = KIT_ATTRS.iter().fold(svg, |svg, attr| scaled(&svg, attr, KIT_SCALE));
        let _ = writeln!(out, "    ({id:?}, \"\", {master:?}),");
    }
    out.push_str("];\n");
    std::fs::write(PathBuf::from(std::env::var("OUT_DIR")?).join("icons.rs"), out)?;
    Ok(())
}

/// FieldCraft's artwork keys its four colour roles (`PALETTE` in `xtask/src/icons.rs`); each
/// becomes the SolidCraft role it stands for (its `palette.toml` source hex): line art on the
/// background (`line`, as `currentColor`), the accent (`active-3`), a soft face (`passive-1`) and
/// warm marks (`component-orange-2`).
const ROLES: [(&str, &str); 4] = [("#dde2ea", "currentColor"), ("#4f9ce8", "#236ead"), ("#2d4e70", "#e2e8ef"), ("#f0b44a", "#c06f00")];

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
