//! The `geometry.features` reply in words: each feature of a body's SolidCraft timeline as one
//! line a person reads ("Extrude 16 mm along +Z", "Hole Ø 14 through", "Fillet R 1.5 on 4
//! edges") and its remaining inputs as labelled details. Pure functions of the reply's JSON; the
//! context panel ([`crate::reply`]) draws them.

use serde_json::Value;

use crate::describe::quantity;

/// What a feature does, for its glyph.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Extrude,
    Cut,
    Revolve,
    Primitive,
    Hole,
    Fillet,
    Chamfer,
    Shell,
    Pattern,
    Mirror,
    Boolean,
    Other,
}

/// One feature made readable.
#[derive(Clone, Debug, PartialEq)]
pub struct Feature {
    pub id: u64,
    pub kind: Kind,
    pub name: String,
    /// What it does with its key inputs, in plain words.
    pub summary: String,
    /// Every other input and fact the reply has: (label, value).
    pub details: Vec<(String, String)>,
    /// The persistent names of the body's faces it made or shaped.
    pub faces: Vec<String>,
    pub error: Option<String>,
}

/// `v` to 3 significant figures, as the panel writes numbers.
fn num(v: f64) -> String {
    if v == 0.0 || !v.is_finite() {
        return quantity(v);
    }
    let scale = 10f64.powi(2 - v.abs().log10().floor() as i32);
    quantity((v * scale).round() / scale)
}

fn mm(v: f64) -> String {
    format!("{} mm", num(v))
}

fn deg(radians: f64) -> String {
    format!("{}°", num(radians.to_degrees()))
}

/// A direction as its axis (`+Z`) when it is one, else its components.
pub fn axis(d: [f64; 3]) -> String {
    const TOL: f64 = 1e-6;
    for (i, name) in ["X", "Y", "Z"].iter().enumerate() {
        let c = d.get(i).copied().unwrap_or(0.0);
        let rest: f64 = d.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, v)| v.abs()).sum();
        if (c.abs() - 1.0).abs() < TOL && rest < TOL {
            return format!("{}{name}", if c > 0.0 { "+" } else { "−" });
        }
    }
    format!("({}, {}, {})", num(d[0]), num(d[1]), num(d[2]))
}

/// A snake-case input key as a label: `taper_angle` → "Taper angle".
fn label(key: &str) -> String {
    let words = key.replace('_', " ");
    let mut c = words.chars();
    c.next().map(|f| f.to_uppercase().collect::<String>() + c.as_str()).unwrap_or_default()
}

/// An input's value with its unit, from its key: angles in degrees, counts and scales bare,
/// everything else a length in mm.
fn input(key: &str, v: f64) -> String {
    if key.contains("angle") {
        deg(v)
    } else if key.contains("quantity") || key.contains("scale") || key.contains("factor") || key == "ribs" {
        num(v)
    } else {
        mm(v)
    }
}

fn plural(n: usize, one: &str, many: &str) -> String {
    format!("{n} {}", if n == 1 { one } else { many })
}

fn vec3(v: &Value) -> Option<[f64; 3]> {
    serde_json::from_value(v.clone()).ok()
}

fn strings(v: &Value) -> Vec<String> {
    v.as_array().into_iter().flatten().filter_map(Value::as_str).map(str::to_string).collect()
}

/// The feature node `f` of a `geometry.features` reply made readable.
pub fn describe(f: &Value) -> Feature {
    let kind_name = f["type"].as_str().unwrap_or_default();
    let op = f["operation"].as_str();
    let def = &f["definition"];
    let get = |k: &str| f["values"][k].as_f64();
    let mut used: Vec<&str> = Vec::new();
    let mut take = |k: &'static str| {
        used.push(k);
        get(k)
    };
    let dir = vec3(&f["direction"]).map(axis);
    let along = dir.as_ref().map(|d| format!(" along {d}")).unwrap_or_default();
    let joined = match op {
        Some("join") => ", joined",
        Some("cut") => ", cut",
        Some("intersect") => ", intersected",
        _ => "",
    };
    let edges_n = def["edges"].as_array().map(Vec::len);
    let on_edges = edges_n.map(|n| format!(" on {}", plural(n, "edge", "edges"))).unwrap_or_default();
    let copies = f["transforms"].as_array().map_or(0, Vec::len);
    let (kind, summary) = match kind_name {
        "extrude" | "revolve" => {
            let (kind, verb) = match op {
                Some("cut") => (Kind::Cut, "Cut"),
                Some("intersect") => (Kind::Cut, "Intersect"),
                Some("join") if kind_name == "extrude" => (Kind::Extrude, "Join"),
                _ if kind_name == "revolve" => (Kind::Revolve, "Revolve"),
                _ => (Kind::Extrude, "Extrude"),
            };
            let mut s = if kind_name == "revolve" {
                format!("{verb} {}", take("angle").map(deg).unwrap_or_default())
            } else {
                let d = take("distance");
                if def["extent"]["through_all"] == true {
                    format!("{verb} through all{along}")
                } else if f["two_sided"] == true {
                    // A symmetric extrude's distance is the total in designs saved since
                    // SolidCraft 0.16 (`convention`), per side before.
                    let total = def["extent"]["direction"] == "symmetric" && def["extent"]["convention"] == "total_distance";
                    let per_side = if total { d.map(|d| d / 2.0) } else { d };
                    format!("{verb} {} both sides", per_side.map(mm).unwrap_or_default())
                } else {
                    format!("{verb} {}{along}", d.map(mm).unwrap_or_default())
                }
            };
            if let Some(a) = take("taper_angle").filter(|a| *a != 0.0) {
                s += &format!(", taper {}", deg(a));
            }
            if let Some(t) = take("thickness") {
                s += &format!(", thin {}", mm(t));
            }
            (kind, s.trim_end().to_string())
        }
        "box" => {
            let (l, w, h) = (take("length"), take("width"), take("height"));
            let dims: Vec<String> = [l, w, h].into_iter().flatten().map(num).collect();
            (Kind::Primitive, format!("Box {} mm{joined}", dims.join(" × ")))
        }
        "cylinder" => {
            let (r, h) = (take("radius"), take("height"));
            let s = format!("Cylinder Ø {} × {} mm{along}{joined}", r.map(|r| num(2.0 * r)).unwrap_or_default(), h.map(num).unwrap_or_default());
            (Kind::Primitive, s)
        }
        "sphere" => (Kind::Primitive, format!("Sphere Ø {}{joined}", take("radius").map(|r| mm(2.0 * r)).unwrap_or_default())),
        "torus" => {
            let (a, b) = (take("major_radius"), take("minor_radius"));
            (Kind::Primitive, format!("Torus R {} × r {}{joined}", a.map(num).unwrap_or_default(), b.map(mm).unwrap_or_default()))
        }
        "hole" => {
            let d = take("diameter").map(num).unwrap_or_default();
            let mut s = match take("depth") {
                Some(depth) => format!("Hole Ø {d} × {} deep", mm(depth)),
                None => format!("Hole Ø {d} through"),
            };
            if let (Some(a), Some(b)) = (take("counterbore_diameter"), take("counterbore_depth")) {
                s += &format!(", counterbore Ø {} × {}", num(a), mm(b));
            }
            if let (Some(a), Some(b)) = (take("countersink_diameter"), take("countersink_angle")) {
                s += &format!(", countersink Ø {} × {}", num(a), deg(b));
            }
            (Kind::Hole, s)
        }
        "fillet" => {
            let r = match (take("radius"), take("start_radius"), take("end_radius")) {
                (Some(r), _, _) => format!("R {}", num(r)),
                (None, Some(a), Some(b)) => format!("R {}–{}", num(a), num(b)),
                _ => String::new(),
            };
            (Kind::Fillet, format!("Fillet {r}{on_edges}").replace("  ", " "))
        }
        "chamfer" => {
            let d = take("distance").map(mm).unwrap_or_default();
            let a = take("angle").map(|a| format!(" × {}", deg(a))).unwrap_or_default();
            (Kind::Chamfer, format!("Chamfer {d}{a}{on_edges}"))
        }
        "shell" => {
            let faces = def["faces"].as_array().map(Vec::len).unwrap_or(0);
            let open = if faces > 0 { format!(", {} open", plural(faces, "face", "faces")) } else { String::new() };
            (Kind::Shell, format!("Shell, {} wall{open}", take("thickness").map(mm).unwrap_or_default()))
        }
        "pattern" | "linear_pattern" | "circular_pattern" => {
            let n = take("quantity");
            let s = match (n, take("spacing").or_else(|| take("distance")), take("total_angle")) {
                (Some(n), _, Some(a)) => format!("Pattern {} × over {}", num(n), deg(a)),
                (Some(n), Some(d), _) => format!("Pattern {} × {} apart", num(n), mm(d)),
                _ => format!("Pattern of {}", plural(copies + 1, "instance", "instances")),
            };
            (Kind::Pattern, s)
        }
        "mirror" => (Kind::Mirror, format!("Mirror, {}", plural(copies.max(1), "copy", "copies"))),
        "combine" => {
            let verb = match op {
                Some("cut") => "cut",
                Some("intersect") => "intersect",
                _ => "join",
            };
            let tools = strings(&f["tools"]);
            let targets = strings(&f["targets"]);
            let s = match (tools.is_empty(), targets.is_empty()) {
                (false, false) => format!("Combine: {verb} {} with {}", targets.join(", "), tools.join(", ")),
                _ => format!("Combine: {verb}"),
            };
            (Kind::Boolean, s)
        }
        other => {
            let mut s = label(other);
            let first: Vec<String> = f["values"]
                .as_object()
                .into_iter()
                .flatten()
                .take(2)
                .filter_map(|(k, v)| v.as_f64().map(|v| format!("{} {}", label(k).to_lowercase(), input(k, v))))
                .collect();
            if !first.is_empty() {
                s += &format!(", {}", first.join(", "));
            }
            (Kind::Other, format!("{s}{joined}"))
        }
    };
    let mut details: Vec<(String, String)> = Vec::new();
    for (k, v) in f["values"].as_object().into_iter().flatten() {
        if let Some(v) = v.as_f64()
            && !used.contains(&k.as_str())
        {
            details.push((label(k), input(k, v)));
        }
    }
    if let Some(sk) = f["sketch"].as_object() {
        let n = sk.get("profiles").and_then(Value::as_array).map_or(0, Vec::len);
        let name = sk.get("name").and_then(Value::as_str).unwrap_or("sketch");
        details.push(("Sketch".into(), format!("{name} ({})", plural(n, "profile", "profiles"))));
    }
    if let Some(d) = dir {
        details.push(("Direction".into(), d));
    }
    if let Some(p) = vec3(&f["placement"]) {
        details.push(("At".into(), format!("({}, {}, {}) mm", num(p[0]), num(p[1]), num(p[2]))));
    }
    if let Some(op) = op {
        details.push(("Operation".into(), op.replace('_', " ")));
    }
    for (key, name) in [("targets", "Changes"), ("tools", "Tools"), ("creates", "Creates")] {
        let v = strings(&f[key]);
        if !v.is_empty() {
            details.push((name.into(), v.join(", ")));
        }
    }
    let faces = strings(&f["faces"]);
    let edges = f["edges"].as_array().map_or(0, Vec::len);
    details.push(("Shapes".into(), format!("{}, {}", plural(faces.len(), "face", "faces"), plural(edges, "edge", "edges"))));
    Feature {
        id: f["id"].as_u64().unwrap_or(0),
        kind,
        name: f["name"].as_str().unwrap_or(kind_name).to_string(),
        summary,
        details,
        faces,
        error: f["error"].as_str().map(str::to_string),
    }
}

/// The feature's glyph, drawn in `rect` with lines in `line` and the part it adds or removes in
/// `accent` (no image assets: every glyph is a few strokes).
pub fn paint_glyph(painter: &egui::Painter, rect: egui::Rect, kind: Kind, line: egui::Color32, accent: egui::Color32) {
    use egui::{Shape, Stroke, pos2, vec2};
    let r = rect.shrink(rect.width() * 0.12);
    let (w, h) = (r.width(), r.height());
    let at = |x: f32, y: f32| pos2(r.left() + x * w, r.top() + y * h);
    let s = Stroke::new(1.3, line);
    let a = Stroke::new(1.5, accent);
    let arrow = |from: egui::Pos2, to: egui::Pos2| {
        painter.line_segment([from, to], a);
        let d = (to - from).normalized() * (w * 0.2);
        let n = vec2(-d.y, d.x);
        painter.add(Shape::convex_polygon(vec![to, to - d + n * 0.8, to - d - n * 0.8], accent, Stroke::NONE));
    };
    match kind {
        Kind::Extrude | Kind::Cut => {
            painter.rect_stroke(egui::Rect::from_min_max(at(0.0, 0.62), at(1.0, 1.0)), 1.0, s, egui::StrokeKind::Inside);
            if kind == Kind::Extrude { arrow(at(0.5, 0.62), at(0.5, 0.0)) } else { arrow(at(0.5, 0.0), at(0.5, 0.6)) }
        }
        Kind::Revolve => {
            painter.circle_stroke(r.center(), w * 0.45, s);
            arrow(at(0.15, 0.25), at(0.45, 0.05));
        }
        Kind::Primitive => {
            painter.rect_stroke(egui::Rect::from_min_max(at(0.0, 0.25), at(0.75, 1.0)), 1.0, s, egui::StrokeKind::Inside);
            painter.add(Shape::line(vec![at(0.0, 0.25), at(0.25, 0.0), at(1.0, 0.0), at(1.0, 0.75), at(0.75, 1.0)], s));
            painter.line_segment([at(0.75, 0.25), at(1.0, 0.0)], s);
        }
        Kind::Hole => {
            painter.rect_stroke(r, 1.0, s, egui::StrokeKind::Inside);
            painter.circle_stroke(r.center(), w * 0.24, a);
        }
        Kind::Fillet => {
            painter.add(Shape::line(vec![at(0.0, 1.0), at(0.0, 0.5)], s));
            painter.add(Shape::line(vec![at(0.5, 0.0), at(1.0, 0.0)], s));
            let arc: Vec<egui::Pos2> = (0..=8)
                .map(|i| {
                    let t = std::f32::consts::FRAC_PI_2 * i as f32 / 8.0;
                    at(0.5 - 0.5 * t.cos(), 0.5 - 0.5 * t.sin())
                })
                .collect();
            painter.add(Shape::line(arc, a));
        }
        Kind::Chamfer => {
            painter.add(Shape::line(vec![at(0.0, 1.0), at(0.0, 0.45)], s));
            painter.add(Shape::line(vec![at(0.55, 0.0), at(1.0, 0.0)], s));
            painter.line_segment([at(0.0, 0.45), at(0.55, 0.0)], a);
        }
        Kind::Shell => {
            painter.add(Shape::line(vec![at(0.0, 0.0), at(0.0, 1.0), at(1.0, 1.0), at(1.0, 0.0)], s));
            painter.add(Shape::line(vec![at(0.25, 0.0), at(0.25, 0.75), at(0.75, 0.75), at(0.75, 0.0)], a));
        }
        Kind::Pattern => {
            for (x, y) in [(0.2, 0.2), (0.8, 0.2), (0.2, 0.8), (0.8, 0.8)] {
                painter.circle_stroke(at(x, y), w * 0.16, if (x, y) == (0.2, 0.2) { s } else { a });
            }
        }
        Kind::Mirror => {
            let mut y = 0.0;
            while y < 1.0 {
                painter.line_segment([at(0.5, y), at(0.5, (y + 0.12).min(1.0))], s);
                y += 0.22;
            }
            painter.add(Shape::convex_polygon(vec![at(0.38, 0.15), at(0.38, 0.85), at(0.0, 0.85)], egui::Color32::TRANSPARENT, s));
            painter.add(Shape::convex_polygon(vec![at(0.62, 0.15), at(0.62, 0.85), at(1.0, 0.85)], egui::Color32::TRANSPARENT, a));
        }
        Kind::Boolean => {
            painter.circle_stroke(at(0.36, 0.5), w * 0.34, s);
            painter.circle_stroke(at(0.64, 0.5), w * 0.34, a);
        }
        Kind::Other => {
            painter.circle_filled(r.center(), w * 0.18, accent);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn summary(f: Value) -> String {
        describe(&f).summary
    }

    #[test]
    fn extrudes_and_cuts_say_how_far_and_which_way() {
        let f = json!({"type": "extrude", "operation": "new_body", "values": {"distance": 16.0}, "direction": [0, 0, 1], "definition": {}});
        assert_eq!(summary(f), "Extrude 16 mm along +Z");
        let f = json!({"type": "extrude", "operation": "cut", "values": {"distance": 8.0}, "direction": [0, 0, -1], "definition": {}});
        assert_eq!(summary(f), "Cut 8 mm along −Z");
        let f = json!({"type": "extrude", "operation": "cut", "values": {"distance": 0.0}, "direction": [0, 0, 1], "definition": {"extent": {"through_all": true}}});
        assert_eq!(summary(f), "Cut through all along +Z");
        let f = json!({"type": "extrude", "operation": "join", "two_sided": true, "values": {"distance": 10.0}, "definition": {}});
        assert_eq!(summary(f), "Join 10 mm both sides");
        let total = json!({"extent": {"direction": "symmetric", "convention": "total_distance"}});
        let f = json!({"type": "extrude", "operation": "join", "two_sided": true, "values": {"distance": 20.0}, "definition": total});
        assert_eq!(summary(f), "Join 10 mm both sides");
        let f = json!({"type": "extrude", "values": {"distance": 5.0, "taper_angle": 0.0523599}, "direction": [0.6, 0.8, 0], "definition": {}});
        assert_eq!(describe(&f).kind, Kind::Extrude);
        assert_eq!(summary(f), "Extrude 5 mm along (0.6, 0.8, 0), taper 3°");
    }

    #[test]
    fn primitives_give_their_size() {
        let f =
            json!({"type": "cylinder", "operation": "join", "values": {"radius": 23.0, "height": 22.0}, "direction": [0, 0, 1], "definition": {}});
        assert_eq!(summary(f), "Cylinder Ø 46 × 22 mm along +Z, joined");
        let f = json!({"type": "box", "values": {"length": 100.0, "width": 40.0, "height": 5.0}, "definition": {}});
        assert_eq!(summary(f), "Box 100 × 40 × 5 mm");
        let f = json!({"type": "sphere", "operation": "cut", "values": {"radius": 2.5}, "definition": {}});
        assert_eq!(summary(f), "Sphere Ø 5 mm, cut");
    }

    #[test]
    fn holes_say_through_or_how_deep() {
        let f = json!({"type": "hole", "operation": "cut", "values": {"diameter": 3.2}, "definition": {}});
        assert_eq!(describe(&f).kind, Kind::Hole);
        assert_eq!(summary(f), "Hole Ø 3.2 through");
        let f = json!({"type": "hole", "values": {"diameter": 6.0, "depth": 10.0, "counterbore_diameter": 11.0, "counterbore_depth": 6.5}, "definition": {}});
        assert_eq!(summary(f), "Hole Ø 6 × 10 mm deep, counterbore Ø 11 × 6.5 mm");
    }

    #[test]
    fn edge_features_count_their_edges() {
        let f = json!({"type": "fillet", "values": {"radius": 2.0}, "definition": {"edges": [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]]}});
        assert_eq!(describe(&f).kind, Kind::Fillet);
        assert_eq!(summary(f), "Fillet R 2 on 4 edges");
        let f = json!({"type": "fillet", "values": {"start_radius": 1.0, "end_radius": 3.0}, "definition": {"edges": [[0, 0, 0]]}});
        assert_eq!(summary(f), "Fillet R 1–3 on 1 edge");
        let f = json!({"type": "chamfer", "values": {"distance": 0.5}, "definition": {"edges": [[0, 0, 0], [1, 1, 1]]}});
        assert_eq!(describe(&f).kind, Kind::Chamfer);
        assert_eq!(summary(f), "Chamfer 0.5 mm on 2 edges");
    }

    #[test]
    fn shells_patterns_mirrors_and_combines() {
        let f = json!({"type": "shell", "values": {"thickness": 2.0}, "definition": {"faces": ["top"]}});
        assert_eq!(summary(f), "Shell, 2 mm wall, 1 face open");
        let f = json!({"type": "pattern", "values": {"quantity": 4.0, "spacing": 30.0}, "definition": {}});
        assert_eq!(describe(&f).kind, Kind::Pattern);
        assert_eq!(summary(f), "Pattern 4 × 30 mm apart");
        let f = json!({"type": "pattern", "values": {"quantity": 6.0, "total_angle": std::f64::consts::TAU}, "definition": {}});
        assert_eq!(summary(f), "Pattern 6 × over 360°");
        let f = json!({"type": "mirror", "values": {}, "transforms": [[[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]], "definition": {}});
        assert_eq!(describe(&f).kind, Kind::Mirror);
        assert_eq!(summary(f), "Mirror, 1 copy");
        let f = json!({"type": "combine", "operation": "cut", "targets": ["Base"], "tools": ["Tool"], "values": {}, "definition": {}});
        assert_eq!(describe(&f).kind, Kind::Boolean);
        assert_eq!(summary(f), "Combine: cut Base with Tool");
    }

    #[test]
    fn other_kinds_name_their_first_inputs_and_details_hold_the_rest() {
        let f = json!({"type": "offset_face", "values": {"distance": 1.5}, "definition": {}});
        assert_eq!(summary(f), "Offset face, distance 1.5 mm");
        let f = json!({
            "id": 2, "name": "Arm", "type": "extrude", "operation": "new_body", "values": {"distance": 16.0, "start_offset": 1.0},
            "direction": [0, 0, 1], "creates": ["CrankArm"], "sketch": {"name": "ArmProfile", "profiles": [{}]},
            "faces": ["F2:top", "F2:start"], "edges": ["e1"], "definition": {}, "error": "no profile"
        });
        let d = describe(&f);
        assert_eq!((d.id, d.name.as_str(), d.faces.len(), d.error.as_deref()), (2, "Arm", 2, Some("no profile")));
        let get = |k: &str| d.details.iter().find(|(l, _)| l == k).map(|(_, v)| v.as_str());
        assert_eq!(get("Start offset"), Some("1 mm"));
        assert_eq!(get("Distance"), None, "the summary already says it");
        assert_eq!(get("Sketch"), Some("ArmProfile (1 profile)"));
        assert_eq!(get("Creates"), Some("CrankArm"));
        assert_eq!(get("Shapes"), Some("2 faces, 1 edge"));
    }

    #[test]
    fn hostile_nodes_never_panic() {
        for f in
            [json!(null), json!([]), json!({"type": 3, "values": "x", "direction": [1]}), json!({"type": "hole", "values": {"diameter": f64::MAX}})]
        {
            let _ = describe(&f);
        }
        assert_eq!(axis([0.0, 0.0, 0.0]), "(0, 0, 0)");
    }
}
