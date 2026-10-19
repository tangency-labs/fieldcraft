//! The tools: their definitions (JSON schemas) and what they run. Every tool is a thin wrapper
//! over engine commands; the same list serves the headless and the connected backend.

use std::path::Path;

use fieldcraft_io::roots::Access;
use serde_json::{Map, Value, json};

use crate::backend::{Backend, PictureRequest};

/// Most commands in one `batch` call.
pub const MAX_STEPS: usize = 10_000;

/// The meshers the `mesh` tool runs, by method; a method whose command this build lacks is
/// left out.
const MESHERS: &[(&str, &str, &str)] = &[
    ("tet", "mesh.tet.generate", "constrained Delaunay tets of every body of the open design, one {body}, or an imported surface {surface}"),
    ("hex", "mesh.hex.generate", "all-hex mesh of a body"),
    ("hex_sweep", "mesh.hex.sweep", "swept all-hex mesh"),
    ("hex_mapped", "mesh.hex.mapped", "mapped all-hex mesh"),
    ("hex_blocks", "mesh.hex.blocks", "one conforming all-hex mesh of a design split into six-sided blocks {bodies, size, order}"),
    (
        "cells",
        "mesh.cells",
        "one mesh of a partitioned part, each cell with its own method {methods: {cell: method}, default, size, order}; cells joined across partitions (merged or tied)",
    ),
    (
        "hex_extrude",
        "mesh.hex.extrude",
        "hex layers from a section (faces, a neighbour's mesh or a 2D mesh) up to a target surface {source | source_set, to_surface, sides, direction | normal, layers | size, bias, body}",
    ),
    ("mirror", "mesh.mirror", "reflect meshed bodies across a datum plane and join them on it {plane, bodies}"),
    ("block", "mesh.block", "structured hex block"),
    (
        "hex_project",
        "mesh.hex.project",
        "projection blocking: a block topology (grid with carve, butterfly or 7-block sphere template, O-grids, extruded and connected faces) projected onto the body and elliptically smoothed {body, template, operations, size}",
    ),
    ("planar", "mesh.planar", "2D elements of a polygon or a planar face"),
];

/// The analyses the `solve` tool runs; one whose command this build lacks is left out.
const ANALYSES: &[(&str, &str)] = &[
    ("static", "solve.static"),
    ("nonlinear", "solve.nonlinear"),
    ("frequency", "solve.frequency"),
    ("buckle", "solve.buckle"),
    ("explicit", "solve.explicit"),
    ("fluid", "solve.fluid"),
    ("heat", "solve.heat"),
];

const VIEWS: &[&str] = &["iso", "home", "fit", "front", "back", "top", "bottom", "left", "right"];

fn available<T: Copy>(table: &'static [T], command: fn(T) -> &'static str) -> impl Iterator<Item = T> {
    table.iter().copied().filter(move |t| fieldcraft_engine::find_command(command(*t)).is_some())
}

fn obj(props: Value, required: &[&str]) -> Value {
    let mut o = json!({"type": "object", "properties": props, "additionalProperties": false});
    if !required.is_empty() {
        o["required"] = json!(required);
    }
    o
}

fn tool(name: &str, title: &str, desc: &str, schema: Value, read_only: bool, destructive: bool) -> Value {
    json!({
        "name": name,
        "title": title,
        "description": desc,
        "inputSchema": schema,
        "annotations": {"title": title, "readOnlyHint": read_only, "destructiveHint": destructive, "idempotentHint": read_only, "openWorldHint": false},
    })
}

fn string(desc: &str) -> Value {
    json!({"type": "string", "description": desc})
}

fn integer(desc: &str) -> Value {
    json!({"type": "integer", "description": desc})
}

/// Every tool, with its JSON schema.
pub fn tool_definitions() -> Vec<Value> {
    let methods: Vec<&str> = available(MESHERS, |m| m.1).map(|m| m.0).collect();
    let method_docs: Vec<String> = available(MESHERS, |m| m.1).map(|(m, c, d)| format!("{m} ({c}: {d})")).collect();
    let analyses: Vec<&str> = available(ANALYSES, |a| a.1).map(|a| a.0).collect();
    vec![
        tool(
            "list_commands",
            "List commands",
            "Every command with its label and parameter documentation. Ids are dotted and lower-case (geometry.open, \
             mesh.tet.generate, model.material.create, model.step.create, solve.static, results.contour). Optional \
             case-insensitive substring `filter` over id and label.",
            obj(json!({"filter": string("Substring to match, e.g. \"model.load\", \"mesh\", \"results\"")}), &[]),
            true,
            false,
        ),
        tool(
            "execute",
            "Execute a command",
            "Run one command with JSON parameters, exactly as the UI, scripts and the CLI do; it never opens a dialog. \
             Examples: model.units {system: \"mm-t-s-N-MPa\"}; model.set.create {name: \"root\", geometry: [{body: \"Body1\", \
             face: \"F1:start\"}]}; model.material.create {name: \"Steel\", elastic: {type: \"isotropic\", young: 210000, \
             poisson: 0.3}, density: 7.85e-9}; model.section.create {name: \"solid\", set: \"all\", material: \"Steel\", \
             kind: {type: \"solid\"}}; model.boundary_condition.create {name: \"fix\", type: \"fixed\", target: {set: \
             \"root\"}}; model.load.create {name: \"push\", type: \"force\", target: {set: \"tip\"}, vector: [0, 0, -100]}; \
             model.step.create {name: \"Load\", procedure: {type: \"static\"}, boundary_conditions: [\"fix\"], loads: \
             [\"push\"]}. Returns the command's result; a failure leaves the model unchanged and says why.",
            obj(
                json!({
                    "command": string("Command id (see list_commands)"),
                    "params": {"type": "object", "description": "The command's parameters"},
                }),
                &["command"],
            ),
            false,
            false,
        ),
        tool(
            "batch",
            "Run several commands",
            "Run commands in order, stopping at the first error. Each step is {command, params}. The reply has every \
             result and, on an error, the failing step (index from 0) and its error; the steps before it stay applied \
             (`undo` takes them back one by one).",
            obj(
                json!({
                    "commands": {
                        "type": "array",
                        "description": "Steps: [{\"command\": id, \"params\": {...}}, ...]",
                        "items": {"type": "object", "properties": {"command": {"type": "string"}, "params": {"type": "object"}}, "required": ["command"]},
                    },
                }),
                &["commands"],
            ),
            false,
            false,
        ),
        tool(
            "inspect_model",
            "Inspect the model",
            "The model: node and element counts, the unit system and the definitions by kind (model.info); the step \
             table (which boundary conditions and loads each step has, model.activity); the open design's bodies with \
             their measures, persistent face names (geometry.bodies), which sets name as {body, face}, and persistent \
             edge names `<face>|<face>` with lengths, which sets name as {body, edge}.",
            obj(json!({}), &[]),
            true,
            false,
        ),
        tool(
            "check_model",
            "Check the model",
            "The model check (model.check): every error and warning, each naming its item (a set member not on the \
             mesh, an element with no section, a missing density for a frequency step, implausible units…). A model \
             with errors does not solve or write. `progress` is the setup checklist: `next` names the next stage and \
             the command that does it, `stages` every stage (geometry … results) as done, next, todo, blocked or \
             optional.",
            obj(json!({}), &[]),
            true,
            false,
        ),
        tool(
            "open",
            "Open a file",
            "Open by extension: a .fieldcraft model (file.open: replaces the model; its design and meshes open again), a \
             SolidCraft design .solidcraft (geometry.open: its bodies become the geometry; `model` \"keep\", the default, \
             keeps the mesh and definitions for a reopen after an edit, \"new\" starts a new model; the reply's `related` \
             and `stale` say what no longer resolves on another part) or a face-tagged .3mf mesh (geometry.import.mesh: \
             added as a surface).",
            obj(
                json!({
                    "path": string("Path to a .fieldcraft, .solidcraft or .3mf file"),
                    "model": {"type": "string", "enum": ["keep", "new"], "description": "For a .solidcraft design: keep the current model (default) or start a new one"},
                }),
                &["path"],
            ),
            false,
            true,
        ),
        tool(
            "save",
            "Save the model",
            "Save the model as a .fieldcraft file (file.save): to `path`, or to the file it came from.",
            obj(json!({"path": string("Path to save to (default: the model's file)")}), &[]),
            false,
            false,
        ),
        tool(
            "write_deck",
            "Write a solver deck",
            "Write the model as a solver input deck (deck.write): Abaqus input (abaqus), Nastran bulk data (nastran) or \
             an LS-DYNA keyword deck (lsdyna). Fails naming every model check error and everything the format cannot \
             express.",
            obj(
                json!({
                    "format": {"type": "string", "enum": ["abaqus", "nastran", "lsdyna"]},
                    "path": string("Output file"),
                    "steps": {"type": "array", "items": {"type": "string"}, "description": "Only these steps (a Nastran run is static or frequency steps)"},
                }),
                &["format", "path"],
            ),
            false,
            false,
        ),
        tool(
            "mesh",
            "Mesh",
            &format!(
                "Mesh the geometry with the chosen method and parameters (list_commands shows each mesher's options), then \
                 report the mesh, the model's counts and the mesh quality (mesh.quality: certified validity, volume, metric \
                 extremes). Methods: {}.",
                method_docs.join("; ")
            ),
            obj(
                json!({
                    "method": {"type": "string", "enum": methods},
                    "params": {"type": "object", "description": "The mesher's parameters, e.g. {\"body\": \"Body1\", \"tolerance\": 0.5}"},
                }),
                &["method"],
            ),
            false,
            false,
        ),
        tool(
            "solve",
            "Solve a step",
            "Solve a step of the model and keep its result: static (linear static), frequency (modes), explicit \
             (explicit dynamics), fluid (lattice Boltzmann), heat when available. `step` may be left out when the model \
             has one step of that kind; `options` takes the solver's caps and settings (list_commands). The results \
             display then shows the step (show_results changes it).",
            obj(
                json!({
                    "analysis": {"type": "string", "enum": analyses},
                    "step": string("Step name"),
                    "options": {"type": "object", "description": "More solve parameters, e.g. {\"max_dofs\": 1000000} or {\"mass\": \"lumped\"}"},
                }),
                &["analysis"],
            ),
            false,
            false,
        ),
        tool(
            "show_results",
            "Show results",
            "Choose what the results display shows and read it back: the step, the deformed shape (`deformed`, `scale`: \
             a factor, 0 for the undeformed shape, or \"auto\"), contours of a nodal field (`field`: displacement, displacement_x/y/z, von_mises, \
             max/mid/min_principal, stress_xx/yy/zz/xy/xz/yz, or none; `bands` 2–24; `range` {min?, max?} or \"auto\" clamps the legend, values above max in the over colour), a frequency step's mode or an \
             explicit step's frame. The reply has the contour's range: its smallest and largest value and their nodes \
             (field displacement gives the largest displacement). Only `step`: show that step.",
            obj(
                json!({
                    "step": string("Solved step (default: the one shown)"),
                    "field": string("Contour field, or \"none\""),
                    "bands": integer("Contour bands, 2 to 24"),
                    "range": {"type": ["object", "string"], "description": "Legend range {\"min\"?, \"max\"?} (each a number or \"auto\"), or \"auto\"; goes with `field`"},
                    "deformed": {"type": "boolean", "description": "Draw the deformed shape"},
                    "scale": {"type": ["number", "string"], "description": "Deformation scale factor (0: undeformed) or \"auto\""},
                    "mode": integer("Mode of a frequency step, from 1"),
                    "frame": integer("Field frame of an explicit step, from 0"),
                    "time": {"type": "number", "description": "Show the explicit step's frame nearest this time"},
                }),
                &[],
            ),
            false,
            false,
        ),
        tool(
            "probe",
            "Probe results",
            "Stored results at nodes, or at the nodes nearest to points `at` (solve.probe): displacement, reaction, stress and von Mises (static), the shape of \
             mode `mode` (frequency), displacement and velocity in field frame `frame` (explicit) and, with `elements`, \
             explicit element values; for a fluid step, the flow at `points` and the force histories of `forces`.",
            obj(
                json!({
                    "step": string("Solved step"),
                    "nodes": {"type": "array", "items": {"type": "integer"}, "description": "Node labels"},
                    "elements": {"type": "array", "items": {"type": "integer"}, "description": "Element labels (explicit)"},
                    "at": {"type": "array", "items": {"type": "array"}, "description": "Points [x, y, z]; the nearest nodes are probed (structural steps)"},
                    "points": {"type": "array", "items": {"type": "array"}, "description": "Points [x, y, z] (fluid)"},
                    "forces": {"type": "array", "items": {"type": "string"}, "description": "Walls or faces whose force histories to return (fluid)"},
                    "mode": integer("Mode, from 1 (frequency)"),
                    "frame": integer("Field frame, from 0 (explicit, fluid)"),
                }),
                &[],
            ),
            true,
            false,
        ),
        tool(
            "reactions",
            "Reaction sums",
            "Reaction force and moment (about `about`) summed over each node set and over the whole model \
             (results.reactions); `sets` defaults to the sets the step's boundary conditions hold.",
            obj(
                json!({
                    "step": string("Solved static step"),
                    "sets": {"type": "array", "items": {"type": "string"}, "description": "Node sets"},
                    "about": {"type": "array", "items": {"type": "number"}, "description": "Moment point [x, y, z] (default origin)"},
                }),
                &[],
            ),
            true,
            false,
        ),
        tool(
            "screenshot",
            "Look at the model",
            "A PNG of the model and the results display. Headless: drawn offscreen framing everything drawn. Connected \
             to the app: the whole window drawn offscreen (source \"window\", the default; ui.capture) or the model with \
             the app's camera (source \"model\"; ui.render). `view` turns the camera to a standard view first, or `camera` \
             aims it: {yaw, pitch} in degrees (yaw 0 looks along +Y, -45 is the iso view's; pitch is the elevation, \
             clamped to ±90) at `target` (default: the centre of mass of the model drawn, so a sequence of yaws orbits \
             it), from `distance` (model units), or with `fit` the targeted bodies filling the frame; without either, \
             from the distance that keeps the whole model in view from every side. A camera given moves the app's too. \
             `section` cuts the mesh (headless).",
            obj(
                json!({
                    "view": {"type": "string", "enum": VIEWS, "description": "Standard view (default: keep the current one; headless starts at iso)"},
                    "camera": {
                        "type": "object",
                        "description": "Camera aim; kept for the next screenshot",
                        "properties": {
                            "yaw": {"type": "number", "description": "Degrees about Z; 0 looks along +Y"},
                            "pitch": {"type": "number", "description": "Degrees of elevation, -90 to 90"},
                            "distance": {"type": "number", "description": "From the target, in model units (not with fit)"},
                            "target": {
                                "description": "What to look at: [x, y, z], or {center_of: mass | bbox, bodies?: [design body names]} (default {center_of: mass}, the model drawn)",
                                "oneOf": [
                                    {"type": "array", "items": {"type": "number"}, "minItems": 3, "maxItems": 3},
                                    {"type": "object", "properties": {
                                        "center_of": {"type": "string", "enum": ["mass", "bbox"]},
                                        "bodies": {"type": "array", "items": {"type": "string"}},
                                    }},
                                ],
                            },
                            "fit": {"type": "boolean", "description": "Choose the distance so the targeted bodies (or the model) fill the frame from this direction"},
                            "fov": {"type": "number", "description": "Vertical field of view in degrees, 0 (orthographic, the default headless) to 150"},
                        },
                        "required": ["yaw", "pitch"],
                    },
                    "section": string("Cut the mesh normal to an axis, as fieldcraft-cli snapshot --section: x, y, z (at the middle) or AXIS=VALUE (y=2.5); headless only"),
                    "source": {"type": "string", "enum": ["window", "model"], "description": "window = the app's window, model = a model render"},
                    "width": integer("Width in pixels (window: points), 16 to 4096"),
                    "height": integer("Height in pixels (window: points), 16 to 4096"),
                    "path": string("Also save the PNG here"),
                }),
                &[],
            ),
            true,
            false,
        ),
        tool(
            "select",
            "Select in the viewport",
            "Change the app window's viewport selection (needs --connect; ui.select): `filter` (faces, edges, vertices, \
             bodies, nodes, element_faces, elements), `mode` (replace, add, remove, toggle), `clear`, items by name and \
             label as model.set.create takes them (`geometry`, `nodes`, `elements`, `faces`), a pick at a window point \
             `at` [x, y] or the visible items in a box `box` [x0, y0, x1, y1]. Returns the selection; its `params` are \
             a set's members.",
            obj(
                json!({
                    "filter": {"type": "string", "enum": ["faces", "edges", "vertices", "bodies", "nodes", "element_faces", "elements"]},
                    "box_mode": {"type": "boolean", "description": "Box selection on or off in the window"},
                    "mode": {"type": "string", "enum": ["replace", "add", "remove", "toggle"]},
                    "clear": {"type": "boolean", "description": "Empty the selection first"},
                    "geometry": {"type": "array", "items": {"type": "object"}, "description": "[{body, face|edge|vertex}]"},
                    "nodes": {"type": "array", "items": {"type": "integer"}},
                    "elements": {"type": "array", "items": {"type": "integer"}},
                    "faces": {"type": "array", "items": {"type": "object"}, "description": "[{element, face}]"},
                    "at": {"type": "array", "items": {"type": "number"}, "description": "Window point [x, y] to pick at"},
                    "box": {"type": "array", "items": {"type": "number"}, "description": "Window box [x0, y0, x1, y1]"},
                }),
                &[],
            ),
            false,
            false,
        ),
        tool(
            "selection",
            "Read the viewport selection",
            "The app window's viewport selection (needs --connect; ui.selection): filter, items, `params` (the \
             selection as a set's members, ready for model.set.create) and what is under the pointer.",
            obj(json!({}), &[]),
            true,
            false,
        ),
        tool("undo", "Undo", "Take back the last command that changed the model.", obj(json!({}), &[]), false, false),
        tool("redo", "Redo", "Run the last undone change again.", obj(json!({}), &[]), false, false),
    ]
}

fn type_matches(t: &str, v: &Value) -> bool {
    match t {
        "string" => v.is_string(),
        "number" => v.is_number(),
        "integer" => v.is_i64() || v.is_u64(),
        "boolean" => v.is_boolean(),
        "array" => v.is_array(),
        "object" => v.is_object(),
        "null" => v.is_null(),
        _ => true,
    }
}

/// Checks `args` against the tool's schema (top level): unknown tool, missing required
/// arguments, unknown arguments, type and enum mismatches.
pub(crate) fn check_arguments(name: &str, args: &Value) -> Result<(), String> {
    let defs = tool_definitions();
    let shown: String = name.chars().take(100).collect();
    let def = defs.iter().find(|d| d["name"] == name).ok_or_else(|| format!("unknown tool `{shown}`"))?;
    let schema = &def["inputSchema"];
    let Some(a) = args.as_object() else { return Err(format!("{name}: arguments must be an object")) };
    for r in schema["required"].as_array().into_iter().flatten().filter_map(Value::as_str) {
        if a.get(r).is_none_or(Value::is_null) {
            return Err(format!("{name}: missing required argument `{r}`"));
        }
    }
    let props = schema["properties"].as_object();
    for (k, v) in a {
        let Some(p) = props.and_then(|p| p.get(k)) else {
            return Err(format!("{name}: unknown argument `{}`", k.chars().take(100).collect::<String>()));
        };
        if v.is_null() {
            continue;
        }
        let ok = match &p["type"] {
            Value::String(t) => type_matches(t, v),
            Value::Array(ts) => ts.iter().filter_map(Value::as_str).any(|t| type_matches(t, v)),
            _ => true,
        };
        if !ok {
            return Err(format!("{name}: argument `{k}` must be of type {}", p["type"]));
        }
        if let Some(allowed) = p["enum"].as_array()
            && !allowed.contains(v)
        {
            return Err(format!("{name}: `{k}` must be one of {}", Value::Array(allowed.clone())));
        }
    }
    Ok(())
}

/// A tool result: one text block and the error flag.
fn text_result(v: &Value, is_error: bool) -> Value {
    let text = match v {
        Value::String(s) => s.clone(),
        v => serde_json::to_string_pretty(v).unwrap_or_default(),
    };
    let mut r = json!({"content": [{"type": "text", "text": text}]});
    if is_error {
        r["isError"] = json!(true);
    }
    r
}

/// Run a tool. A failure is a tool result with `isError: true` (the model reads the message and
/// corrects); nothing here panics on hostile arguments.
pub fn call_tool(b: &mut dyn Backend, name: &str, args: &Value) -> Value {
    let empty = Map::new();
    let a = args.as_object().unwrap_or(&empty);
    let r = match name {
        "screenshot" => return screenshot(b, a),
        "batch" => return batch(b, a),
        _ => dispatch(b, name, a),
    };
    match r {
        Ok(v) => text_result(&v, false),
        Err(e) => text_result(&Value::String(e), true),
    }
}

fn pick(a: &Map<String, Value>, keys: &[&str]) -> Value {
    Value::Object(keys.iter().filter_map(|k| a.get(*k).filter(|v| !v.is_null()).map(|v| (k.to_string(), v.clone()))).collect())
}

fn str_arg<'a>(a: &'a Map<String, Value>, k: &str) -> Option<&'a str> {
    a.get(k).and_then(Value::as_str)
}

/// An object argument (absent: empty).
fn object_arg(a: &Map<String, Value>, k: &str) -> Value {
    a.get(k).cloned().filter(Value::is_object).unwrap_or_else(|| json!({}))
}

fn dispatch(b: &mut dyn Backend, name: &str, a: &Map<String, Value>) -> Result<Value, String> {
    match name {
        "list_commands" => {
            let all = b.commands()?;
            let Some(f) = str_arg(a, "filter").map(str::to_lowercase).filter(|f| !f.is_empty()) else { return Ok(all) };
            let hit = |c: &Value| ["id", "label"].iter().any(|k| c.get(*k).and_then(Value::as_str).is_some_and(|s| s.to_lowercase().contains(&f)));
            Ok(Value::Array(all.as_array().map(|v| v.iter().filter(|c| hit(c)).cloned().collect()).unwrap_or_default()))
        }
        "execute" => {
            let id = str_arg(a, "command").ok_or("`command` is required")?;
            b.execute(id, object_arg(a, "params"))
        }
        "inspect_model" => Ok(json!({
            "model": b.execute("model.info", json!({}))?,
            "steps": b.execute("model.activity", json!({}))?,
            "bodies": b.execute("geometry.bodies", json!({}))?,
        })),
        "check_model" => b.execute("model.check", json!({})),
        "open" => {
            let path = str_arg(a, "path").ok_or("`path` is required")?;
            let ext = Path::new(path).extension().and_then(|e| e.to_str()).map(str::to_ascii_lowercase);
            let command = match ext.as_deref() {
                Some("fieldcraft") => "file.open",
                Some("solidcraft") => "geometry.open",
                Some("3mf") => "geometry.import.mesh",
                _ => return Err("open: a .fieldcraft model, a .solidcraft design or a .3mf mesh".into()),
            };
            let mut params = json!({"path": path});
            if let Some(model) = a.get("model").filter(|m| !m.is_null()) {
                if command != "geometry.open" {
                    return Err("open: `model` applies to a .solidcraft design only".into());
                }
                params["model"] = model.clone();
            }
            let r = b.execute(command, params)?;
            Ok(json!({"command": command, "result": r}))
        }
        "save" => b.execute("file.save", pick(a, &["path"])),
        "write_deck" => b.execute("deck.write", pick(a, &["format", "path", "steps"])),
        "mesh" => {
            let method = str_arg(a, "method").ok_or("`method` is required")?;
            let (_, command, _) = available(MESHERS, |m| m.1).find(|m| m.0 == method).ok_or_else(|| format!("unknown mesh method `{method}`"))?;
            let mesh = b.execute(command, object_arg(a, "params"))?;
            let model = b.execute("model.info", json!({}))?;
            let quality = b.execute("mesh.quality", json!({})).unwrap_or_else(|e| json!({"error": e}));
            Ok(json!({"command": command, "mesh": mesh, "model": model, "quality": quality}))
        }
        "solve" => {
            let analysis = str_arg(a, "analysis").ok_or("`analysis` is required")?;
            let (_, command) = available(ANALYSES, |x| x.1).find(|x| x.0 == analysis).ok_or_else(|| format!("unknown analysis `{analysis}`"))?;
            let mut params = object_arg(a, "options");
            if let Some(step) = a.get("step").filter(|s| !s.is_null()) {
                params["step"] = step.clone();
            }
            b.execute(command, params)
        }
        "show_results" => show_results(b, a),
        "probe" => b.execute("solve.probe", pick(a, &["step", "nodes", "at", "elements", "points", "forces", "mode", "frame"])),
        "reactions" => b.execute("results.reactions", pick(a, &["step", "sets", "about"])),
        "select" => {
            let p = pick(a, &["filter", "box_mode", "mode", "clear", "geometry", "nodes", "elements", "faces", "at", "box"]);
            if p.as_object().is_none_or(|o| o.is_empty()) {
                return Err("select: give `filter`, `mode`, `clear`, items, `at` or `box`".into());
            }
            b.window("ui.select", p)
        }
        "selection" => b.window("ui.selection", json!({})),
        "undo" => b.execute("edit.undo", json!({})),
        "redo" => b.execute("edit.redo", json!({})),
        other => Err(format!("unknown tool `{}`", other.chars().take(100).collect::<String>())),
    }
}

/// `show_results`: the `results.*` commands its arguments ask for, in a fixed order; the reply of
/// the last one (each replies with the whole display).
fn show_results(b: &mut dyn Backend, a: &Map<String, Value>) -> Result<Value, String> {
    let step = pick(a, &["step"]);
    let with_step = |mut p: Value| {
        if let (Some(o), Some(s)) = (p.as_object_mut(), step.get("step")) {
            o.insert("step".into(), s.clone());
        }
        p
    };
    let mut last = None;
    if let Some(mode) = a.get("mode").filter(|v| !v.is_null()) {
        last = Some(b.execute("results.modes", with_step(json!({"mode": mode})))?);
    }
    if a.get("frame").is_some_and(|v| !v.is_null()) || a.get("time").is_some_and(|v| !v.is_null()) {
        last = Some(b.execute("results.frame", with_step(pick(a, &["frame", "time"])))?);
    }
    if a.get("deformed").is_some_and(|v| !v.is_null()) || a.get("scale").is_some_and(|v| !v.is_null()) {
        let mut p = pick(a, &["scale"]);
        if let Some(show) = a.get("deformed").filter(|v| !v.is_null()) {
            p["show"] = show.clone();
        }
        last = Some(b.execute("results.deformed", with_step(p))?);
    }
    if a.get("field").is_some_and(|v| !v.is_null()) {
        last = Some(b.execute("results.contour", with_step(pick(a, &["field", "bands", "range"])))?);
    } else if a.get("bands").is_some_and(|v| !v.is_null()) || a.get("range").is_some_and(|v| !v.is_null()) {
        return Err("show_results: `bands` and `range` go with `field`".into());
    }
    match last {
        Some(v) => Ok(v),
        // Only a step: show it.
        None if step.get("step").is_some() => b.execute("results.deformed", step),
        None => Err("show_results: give `step`, `field`, `deformed`, `scale`, `mode`, `frame` or `time`".into()),
    }
}

/// `batch`: run steps in order, stopping at the first failure.
fn batch(b: &mut dyn Backend, a: &Map<String, Value>) -> Value {
    let Some(steps) = a.get("commands").and_then(Value::as_array) else {
        return text_result(&json!("`commands` must be an array of {command, params}"), true);
    };
    if steps.len() > MAX_STEPS {
        return text_result(&json!(format!("batch too long ({} steps, at most {MAX_STEPS})", steps.len())), true);
    }
    let mut results = Vec::with_capacity(steps.len());
    for (i, step) in steps.iter().enumerate() {
        let command = step.get("command").and_then(Value::as_str);
        let r = match command {
            Some(id) => b.execute(id, step.get("params").cloned().filter(Value::is_object).unwrap_or_else(|| json!({}))),
            None => Err("missing `command`".to_string()),
        };
        match r {
            Ok(v) => results.push(v),
            Err(e) => {
                return text_result(
                    &json!({"ok": false, "completed": i, "failed": {"index": i, "command": command, "error": e}, "results": results}),
                    true,
                );
            }
        }
    }
    text_result(&json!({"ok": true, "completed": results.len(), "results": results}), false)
}

/// `screenshot`: a PNG as an MCP image content block, with what it shows as text.
fn screenshot(b: &mut dyn Backend, a: &Map<String, Value>) -> Value {
    match capture(b, a) {
        Ok((png, info)) => json!({"content": [
            {"type": "image", "data": base64(&png), "mimeType": "image/png"},
            {"type": "text", "text": info.to_string()},
        ]}),
        Err(e) => text_result(&Value::String(e), true),
    }
}

fn capture(b: &mut dyn Backend, a: &Map<String, Value>) -> Result<(Vec<u8>, Value), String> {
    // The file is checked before anything is drawn.
    let keep = match str_arg(a, "path").filter(|p| !p.trim().is_empty()) {
        Some(p) => Some(b.roots().check(Path::new(p), Access::Write).map_err(|e| e.to_string())?),
        None => None,
    };
    let size = |k: &str| a.get(k).and_then(Value::as_u64).map(|v| u32::try_from(v).unwrap_or(u32::MAX));
    let camera = a.get("camera").filter(|v| !v.is_null()).map(fieldcraft_engine::view::camera).transpose().map_err(|e| e.to_string())?;
    let section = match a.get("section").filter(|v| !v.is_null()) {
        None => None,
        Some(v) => Some(
            v.as_str()
                .and_then(fieldcraft_engine::view::Section::parse)
                .ok_or("screenshot: `section` is x, y, z or AXIS=VALUE (x=12.5): the mesh cut normal to that axis")?,
        ),
    };
    if camera.is_some() && a.get("view").is_some_and(|v| !v.is_null()) {
        return Err("screenshot: give `view` or `camera`, not both".into());
    }
    let request = PictureRequest {
        view: str_arg(a, "view").map(str::to_string),
        camera,
        section,
        window: str_arg(a, "source") != Some("model"),
        width: size("width"),
        height: size("height"),
    };
    let picture = b.picture(&request)?;
    let mut info = picture.info;
    if let Some(path) = &keep {
        fieldcraft_io::fs::write_output(path, &picture.png).map_err(|e| e.to_string())?;
        info["path"] = json!(path.to_string_lossy());
    }
    info["bytes"] = json!(picture.png.len());
    Ok((picture.png, info))
}

/// Standard base64 with padding (RFC 4648).
pub(crate) fn base64(data: &[u8]) -> String {
    const T: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(data.len().div_ceil(3) * 4);
    for chunk in data.chunks(3) {
        let b = [chunk.first().copied().unwrap_or(0), chunk.get(1).copied().unwrap_or(0), chunk.get(2).copied().unwrap_or(0)];
        let n = (u32::from(b[0]) << 16) | (u32::from(b[1]) << 8) | u32::from(b[2]);
        for (i, shift) in [18u32, 12, 6, 0].into_iter().enumerate() {
            let digit = T.get(((n >> shift) & 63) as usize).copied().unwrap_or(b'=');
            out.push(char::from(if i <= chunk.len() { digit } else { b'=' }));
        }
    }
    out
}
