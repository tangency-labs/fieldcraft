//! `cargo xtask icons [ID…]`: draw the toolbar icons with VectorCraft.
//!
//! Each icon is a recipe, `assets/icons/src/<id>.json`: a list of shapes on a 24 × 24 pixel
//! artboard, each a path (`d`), `circle` `[cx, cy, r]`, `rect` `[x, y, w, h]` (optional
//! `radius`), `line` `[x1, y1, x2, y2]` or `star` `[cx, cy, outer, inner, points]`, painted with
//! palette roles (`fill`, `stroke`: `line`, `accent`, `soft`, `warm`; a shape with neither is
//! stroked with `line`), an optional stroke `width` (default 1.5; round caps and joins) and
//! `dash` (`[dash, gap]`). The task turns every recipe into VectorCraft
//! commands, runs them headlessly with `vectorcraft-cli run` (`VECTORCRAFT_CLI`, else
//! `vectorcraft-cli` on the PATH) and writes `assets/icons/<id>.vectorcraft` (the editable
//! source) and `assets/icons/<id>.svg` (what the app draws). The roles' fixed colours are
//! keys: `crates/ui-egui/build.rs` maps each to the SolidCraft colour role it stands for, and the
//! UI kit draws those in the theme's colours (`crates/ui-egui/src/icons.rs`).

use std::path::Path;
use std::process::Command;

use serde_json::{Value, json};

/// Palette roles and the fixed colour each is drawn in: keys the app maps to SolidCraft's colour
/// roles. Keep in step with `ROLES` in `crates/ui-egui/build.rs`.
pub const PALETTE: &[(&str, &str)] = &[("line", "#dde2ea"), ("accent", "#4f9ce8"), ("soft", "#2d4e70"), ("warm", "#f0b44a")];

/// The artboard side (pixels) and the default stroke weight.
const SIZE: f64 = 24.0;
const STROKE: f64 = 1.5;

pub fn run(root: &Path, only: &[&str]) -> Result<(), String> {
    let dir = root.join("assets/icons");
    let src = dir.join("src");
    let mut recipes: Vec<String> = std::fs::read_dir(&src)
        .map_err(|e| format!("{}: {e}", src.display()))?
        .flatten()
        .filter_map(|e| e.file_name().to_str().and_then(|n| n.strip_suffix(".json")).map(str::to_string))
        .filter(|id| only.is_empty() || only.contains(&id.as_str()))
        .collect();
    recipes.sort();
    if recipes.is_empty() {
        return Err("no icon recipes to draw".into());
    }
    let cli = std::env::var("VECTORCRAFT_CLI").unwrap_or_else(|_| "vectorcraft-cli".into());
    for id in &recipes {
        let path = src.join(format!("{id}.json"));
        let text = std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
        let recipe: Value = serde_json::from_str(&text).map_err(|e| format!("{}: {e}", path.display()))?;
        let steps = commands(id, &recipe, &dir.join(format!("{id}.vectorcraft"))).map_err(|e| format!("{id}: {e}"))?;
        let mut c = Command::new(&cli);
        c.arg("run");
        for (cmd, params) in &steps {
            c.args(["--cmd", cmd, "--params", &params.to_string()]);
        }
        c.arg("--export").arg(dir.join(format!("{id}.svg")));
        let out = c.output().map_err(|e| format!("{cli}: {e} (set VECTORCRAFT_CLI to the vectorcraft-cli binary)"))?;
        if !out.status.success() {
            return Err(format!("{id}: vectorcraft-cli failed:\n{}{}", String::from_utf8_lossy(&out.stdout), String::from_utf8_lossy(&out.stderr)));
        }
    }
    println!("icons: drew {} icon(s) into assets/icons", recipes.len());
    Ok(())
}

fn colour(role: &Value) -> Result<Value, String> {
    match role {
        Value::Null => Ok(json!({"none": true})),
        Value::String(r) => {
            PALETTE.iter().find(|(name, _)| name == r).map(|(_, hex)| json!({"color": hex})).ok_or_else(|| format!("unknown role `{r}`"))
        }
        other => Err(format!("a role is a name, not {other}")),
    }
}

fn numbers<const N: usize>(v: &Value) -> Result<[f64; N], String> {
    let a = v.as_array().filter(|a| a.len() == N).ok_or_else(|| format!("expected {N} numbers, got {v}"))?;
    let mut out = [0.0; N];
    for (o, x) in out.iter_mut().zip(a) {
        *o = x.as_f64().filter(|x| x.is_finite()).ok_or_else(|| format!("expected numbers, got {v}"))?;
    }
    Ok(out)
}

/// The VectorCraft commands that draw one recipe and save it to `source`.
fn commands(id: &str, recipe: &Value, source: &Path) -> Result<Vec<(String, Value)>, String> {
    let shapes = recipe.get("shapes").and_then(Value::as_array).ok_or("a recipe is {\"shapes\": [...]}")?;
    let mut out = vec![("file.new".to_string(), json!({"width": SIZE, "height": SIZE, "units": "Pixels", "name": id, "created": null}))];
    for s in shapes {
        let create = if let Some(d) = s.get("d").and_then(Value::as_str) {
            ("path.create", json!({"d": d}))
        } else if let Some(c) = s.get("circle") {
            let [cx, cy, r] = numbers::<3>(c)?;
            ("shape.ellipse", json!({"x": cx - r, "y": cy - r, "width": 2.0 * r, "height": 2.0 * r}))
        } else if let Some(r) = s.get("rect") {
            let [x, y, width, height] = numbers::<4>(r)?;
            let radius = s.get("radius").and_then(Value::as_f64).unwrap_or(0.0);
            ("shape.rectangle", json!({"x": x, "y": y, "width": width, "height": height, "radius": radius}))
        } else if let Some(l) = s.get("line") {
            // A path, so the exported SVG keeps to the elements the UI kit's icon reader takes
            // (paths, circles, groups): VectorCraft exports a line shape as `<line>`.
            let [x1, y1, x2, y2] = numbers::<4>(l)?;
            ("path.create", json!({"d": format!("M{x1} {y1} L{x2} {y2}")}))
        } else if let Some(st) = s.get("star") {
            let [cx, cy, radius1, radius2, points] = numbers::<5>(st)?;
            ("shape.star", json!({"cx": cx, "cy": cy, "radius1": radius1, "radius2": radius2, "points": points}))
        } else {
            return Err(format!("a shape needs `d`, `circle`, `rect`, `line` or `star`: {s}"));
        };
        out.push((create.0.to_string(), create.1));
        let (fill, stroke) = match (s.get("fill"), s.get("stroke")) {
            (None, None) => (Value::Null, json!("line")),
            (f, st) => (f.cloned().unwrap_or(Value::Null), st.cloned().unwrap_or(Value::Null)),
        };
        out.push(("paint.setFill".into(), colour(&fill)?));
        out.push(("paint.setStroke".into(), colour(&stroke)?));
        if !stroke.is_null() {
            let width = s.get("width").and_then(Value::as_f64).unwrap_or(STROKE);
            let dash = s.get("dash").cloned().unwrap_or(Value::Null);
            out.push(("stroke.set".into(), json!({"weight": width, "cap": "round", "join": "round", "dash": dash})));
        }
    }
    out.push(("document.save".into(), json!({"path": source, "compress": true, "modified": null})));
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_recipe_becomes_vectorcraft_commands() {
        let recipe = json!({"shapes": [
            {"d": "M4 4 L20 20"},
            {"circle": [12, 12, 3], "fill": "accent"},
            {"rect": [2, 2, 4, 4], "radius": 1, "fill": "soft", "stroke": "line", "width": 1},
        ]});
        let steps = commands("x.y", &recipe, Path::new("/tmp/x.y.vectorcraft")).unwrap();
        let ids: Vec<&str> = steps.iter().map(|(c, _)| c.as_str()).collect();
        assert_eq!(
            ids,
            [
                "file.new",
                "path.create",
                "paint.setFill",
                "paint.setStroke",
                "stroke.set",
                "shape.ellipse",
                "paint.setFill",
                "paint.setStroke",
                "shape.rectangle",
                "paint.setFill",
                "paint.setStroke",
                "stroke.set",
                "document.save"
            ]
        );
        assert_eq!(steps[5].1, json!({"x": 9.0, "y": 9.0, "width": 6.0, "height": 6.0}));
        assert_eq!(steps[3].1, json!({"color": "#dde2ea"}));
        assert!(commands("x", &json!({"shapes": [{"d": "M0 0", "fill": "pink"}]}), Path::new("x")).is_err());
        assert!(commands("x", &json!({"shapes": [{"circle": [1, 2]}]}), Path::new("x")).is_err());
    }
}
