//! The engine's unit boundary (docs/units.md).
//!
//! The model stores SI and the solvers work in it. Commands take numbers in the session's
//! display system ([`crate::Session::units`]: the model's own, or the session's default): a bare
//! number is in the display system's unit for its quantity, `"210 GPa"` and
//! `{"value": 210, "unit": "GPa"}` carry their own. Before a command runs, its parameters are
//! converted, field by field from its form ([`crate::schema`]), into the system the command
//! works in ([`frame`]): SI for the model and the solvers, millimetres for the geometry and the
//! meshers (SolidCraft geometry is in millimetres and the meshers' tolerances are tuned for it;
//! the model is expressed in it around such a command, [`crate::Session`]'s frame). Its reply is
//! converted back into the display system ([`REPLIES`]), with the units it is in. The frame is
//! the command's own business: what goes in and comes out is in the display system, like every
//! other command's, and a message that names a number of the frame shows it through [`shown`].
//!
//! The numbers inside a JSON-valued parameter are described by its field's
//! [`Field::numbers`](crate::schema::Field::numbers) and converted alike.

use std::sync::OnceLock;

use fieldcraft_model::units::{DisplaySystem, Preset, Quantity, Unit, UnitError, ValueInput, WithUnit};
use serde_json::{Map, Value, json};

use crate::schema::{Dim, Field, Kind, MAX_DEPTH, Numbers, join};
use crate::{EngineError, Result};

/// SI: what the model stores and the solvers work in.
pub fn si() -> &'static DisplaySystem {
    static SI: OnceLock<DisplaySystem> = OnceLock::new();
    SI.get_or_init(|| DisplaySystem::new(Preset::MKgS))
}

/// What the geometry and the meshers work in: millimetres, in mm-t-s-N-MPa.
pub fn geometry() -> &'static DisplaySystem {
    fieldcraft_model::units::geometry_frame()
}

/// A length of the geometry (millimetres) in SI.
pub fn geometry_length(mm: f64) -> f64 {
    geometry().to_si(Quantity::Length, mm).unwrap_or(f64::NAN)
}

/// A point of the geometry (millimetres) in SI.
pub fn geometry_point(p: [f64; 3]) -> [f64; 3] {
    p.map(geometry_length)
}

/// A length in SI as the geometry (millimetres) has it.
pub fn to_geometry_length(m: f64) -> f64 {
    geometry().express(Quantity::Length, m).unwrap_or(f64::NAN)
}

/// A point in SI as the geometry (millimetres) has it.
pub fn to_geometry_point(p: [f64; 3]) -> [f64; 3] {
    p.map(to_geometry_length)
}

/// A session to preview a geometry command on (the interactive tools): the model expressed in
/// millimetres, the design and the imported surfaces shared. Nothing else of `s` is kept.
pub(crate) fn geometry_session(s: &crate::Session) -> Result<crate::Session> {
    let mm = geometry();
    let model = s.framed_model(mm)?;
    Ok(crate::Session { model, design: s.design.clone(), surfaces: s.surfaces.clone(), ..crate::Session::default() })
}

/// The parameters of command `id`, typed in `s`'s input system, in millimetres: what a preview of
/// a geometry command reads.
pub(crate) fn geometry_params(s: &crate::Session, id: &str, params: &Value) -> Result<Value> {
    match crate::find_command(id).and_then(|c| c.form) {
        Some(form) => convert_params(form, params, s.input_units(), geometry()),
        None => Ok(params.clone()),
    }
}

/// A number of `q` in `s`'s frame `from` (the geometry's millimetres, say) as a message shows it:
/// in the display system, with its unit (`0.19685 in`).
pub fn shown(display: &DisplaySystem, from: &DisplaySystem, q: Quantity, v: f64) -> String {
    from.to_si(q, v).and_then(|si| q.format(si, display)).unwrap_or_else(|_| fieldcraft_model::units::format_number(v))
}

/// `what` (a parameter of a command of the geometry frame) below `fraction` of the bounding-box
/// diagonal `diagonal` of `of` ("body `Base`'s"): the error names the fraction and the length it
/// is, which the command's boundary shows in the display units ([`crate::EngineError::shown`]).
pub(crate) fn below_diagonal(what: &str, fraction: f64, diagonal: f64, of: &str) -> EngineError {
    EngineError::Mesh(fieldcraft_mesh::Error::invalid(fieldcraft_mesh::message![
        format!("`{what}` must be at least {fraction:e} of {of} bounding-box diagonal ("),
        fieldcraft_mesh::message::length(fraction * diagonal),
        ")"
    ]))
}

/// A length of the geometry (millimetres) as a message to the user of `s` shows it.
pub fn shown_length(s: &crate::Session, mm: f64) -> String {
    shown(s.units(), geometry(), Quantity::Length, mm)
}

/// An area of the geometry (mm²) as a message to the user of `s` shows it.
pub fn shown_area(s: &crate::Session, mm2: f64) -> String {
    shown(s.units(), geometry(), Quantity::Area, mm2)
}

/// The parameters of command `id` written in millimetres (a command a reply proposes to run
/// next), in `s`'s display system: what running it reads.
pub(crate) fn geometry_params_shown(s: &crate::Session, id: &str, params: &Value) -> Result<Value> {
    match crate::find_command(id).and_then(|c| c.form) {
        Some(form) => convert_params(form, params, geometry(), s.units()),
        None => Ok(params.clone()),
    }
}

/// How the meshers' messages show their numbers ([`fieldcraft_mesh::message`]): in `from` (the
/// command's frame), shown in `to` (the display system).
pub struct MeshUnits<'a> {
    pub from: &'a DisplaySystem,
    pub to: &'a DisplaySystem,
}

impl fieldcraft_mesh::message::Units for MeshUnits<'_> {
    fn value(&self, q: Quantity, v: f64) -> f64 {
        self.from.to_si(q, v).and_then(|si| self.to.express(q, si)).unwrap_or(v)
    }
    fn symbol(&self, q: Quantity) -> String {
        let s = self.to.symbol(q);
        if s == "1" { String::new() } else { s.to_string() }
    }
}

/// A mesher's error (its numbers in the geometry's millimetres) as text in `display`.
pub fn mesh_error_text(e: &fieldcraft_mesh::Error, display: &DisplaySystem) -> String {
    e.clone().shown(&MeshUnits { from: geometry(), to: display }).to_string()
}

/// A mesher's message (its numbers in the geometry's millimetres) as `s`'s user reads it.
pub fn mesh_text(s: &crate::Session, m: &fieldcraft_mesh::message::Message) -> String {
    m.render(&MeshUnits { from: geometry(), to: s.units() })
}

pub use fieldcraft_model::units::{MIN_GEOMETRY_LENGTH as MIN_LENGTH, geometry_amount, min_geometry_length as min_length};

/// The system a command works in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Frame {
    /// SI: the model as stored (definitions, solves, results).
    Si,
    /// Millimetres: commands on the geometry, and the meshers.
    Geometry,
    /// The system its numbers were typed in (the display system): `mesh.block`, which meshes from
    /// typed corners and so gives the same mesh as before SI storage, to the last bit, and
    /// `mesh.export.stl`, whose file (which has no unit of its own) is in the display length.
    Input,
}

/// The frame of command `id`.
pub fn frame(id: &str) -> Frame {
    const INPUT: &[&str] = &["mesh.block", "mesh.export.stl", "mesh.import.orphan"];
    const SI_MESH: &[&str] = &[
        "mesh.node.add",
        "mesh.node.remove",
        "mesh.element.add",
        "mesh.element.remove",
        "mesh.clear",
        "mesh.renumber",
        "mesh.tag",
        "mesh.tag.clear",
    ];
    if INPUT.contains(&id) {
        Frame::Input
    } else if SI_MESH.contains(&id) {
        Frame::Si
    } else if id.starts_with("geometry.") || id.starts_with("mesh.") || id == "model.boundary_condition.symmetry" {
        Frame::Geometry
    } else {
        Frame::Si
    }
}

/// Whether `a` and `b` read a bare number of `q` alike (so passing it on changes nothing).
fn same_unit(a: &Unit, b: &Unit) -> bool {
    a == b || (a.dimension() == b.dimension() && a.factor() == b.factor() && a.is_si() == b.is_si() && a.symbol() == b.symbol())
}

fn unit_error(path: &str, e: UnitError) -> EngineError {
    EngineError::invalid(path, e.to_string())
}

/// The SI value of a command's value of quantity `q`: a bare number in `input`, or a value with
/// its unit; `None` when it is not one.
pub fn read(v: &Value, q: Quantity, input: &DisplaySystem) -> Option<std::result::Result<f64, UnitError>> {
    let input_value: ValueInput = match v {
        Value::Number(_) | Value::String(_) | Value::Object(_) => serde_json::from_value(v.clone()).ok()?,
        _ => return None,
    };
    Some(input_value.to_si(q, input))
}

/// What to do with each number a form describes.
trait Leaf {
    /// `v` at `path`, of quantity `q` (`None`: no unit), a bare number in `fixed` when set.
    fn number(&mut self, v: &mut Value, q: Option<Quantity>, fixed: Option<&'static str>, path: &str) -> Result<()>;
}

/// Converts every number to `target`: bare ones are in `input` (or their fixed unit).
struct ToTarget<'a> {
    input: &'a DisplaySystem,
    target: &'a DisplaySystem,
}

impl ToTarget<'_> {
    /// The value in the target system, from a bare number in `from`.
    fn bare(&self, x: f64, q: Quantity, from: &Unit, path: &str) -> Result<f64> {
        let to = self.target.unit_ref(q);
        if same_unit(from, to) {
            return if x.is_finite() { Ok(x) } else { Err(EngineError::param(path, crate::schema::ProblemKind::NotANumber)) };
        }
        let si = q.to_si(x, from).map_err(|e| unit_error(path, e))?;
        q.express(si, to).map_err(|e| unit_error(path, e))
    }
}

impl Leaf for ToTarget<'_> {
    fn number(&mut self, v: &mut Value, q: Option<Quantity>, fixed: Option<&'static str>, path: &str) -> Result<()> {
        // A field of a fixed unit is given to its command in that unit; the others in the target.
        let fixed_unit = fixed.map(Unit::parse).transpose().map_err(|e| unit_error(path, e))?;
        let out = match (&*v, q) {
            (Value::Null, _) => return Ok(()),
            // A number without a unit passes as it is, and so does a bare number of a fixed unit.
            (Value::Number(_), None) => return Ok(()),
            (Value::Number(_), Some(_)) if fixed_unit.is_some() => return Ok(()),
            (Value::Number(n), Some(q)) => {
                let x = n.as_f64().unwrap_or(f64::NAN);
                self.bare(x, q, self.input.unit_ref(q), path)?
            }
            // Text and objects of a field without a unit are the command's to read.
            (Value::String(_) | Value::Object(_), None) => return Ok(()),
            (Value::String(_) | Value::Object(_), Some(q)) => {
                let input: ValueInput = serde_json::from_value(v.clone())
                    .map_err(|_| EngineError::invalid(path, "must be a number, a value with its unit (\"210 GPa\") or {\"value\", \"unit\"}"))?;
                let si = match (&input, &fixed_unit) {
                    // Bare text in a fixed-unit field is in its fixed unit.
                    (ValueInput::Text(t), Some(unit)) if t.trim().parse::<f64>().is_ok() => {
                        q.to_si(t.trim().parse::<f64>().unwrap_or(f64::NAN), unit)
                    }
                    _ => input.to_si(q, self.input),
                }
                .map_err(|e| unit_error(path, e))?;
                let to = fixed_unit.as_ref().unwrap_or_else(|| self.target.unit_ref(q));
                q.express(si, to).map_err(|e| unit_error(path, e))?
            }
            // Not a number: left for the command to refuse.
            _ => return Ok(()),
        };
        *v = serde_json::Number::from_f64(out).map(Value::Number).ok_or_else(|| EngineError::invalid(path, "is out of range"))?;
        Ok(())
    }
}

/// Gives every bare number that follows the display system its unit, `{"value", "unit"}`: a
/// journal entry typed in another system than its journal's.
struct Annotate<'a> {
    input: &'a DisplaySystem,
}

impl Leaf for Annotate<'_> {
    fn number(&mut self, v: &mut Value, q: Option<Quantity>, fixed: Option<&'static str>, _: &str) -> Result<()> {
        if let (Value::Number(n), Some(q), None) = (&*v, q, fixed)
            && let Some(x) = n.as_f64()
        {
            let unit = self.input.symbol(q);
            if unit != "1" {
                *v = json!(WithUnit { value: x, unit: unit.to_string() });
            }
        }
        Ok(())
    }
}

/// Walk `params` by `fields`, calling `leaf` on every number.
fn walk(fields: &[Field], obj: &mut Map<String, Value>, prefix: &str, leaf: &mut dyn Leaf, depth: usize) -> Result<()> {
    if depth > MAX_DEPTH {
        return Ok(());
    }
    let parent = Value::Object(obj.clone());
    for f in fields {
        if let Kind::Variants { tag, .. } = f.kind
            && f.key.is_empty()
        {
            if let Some(v) = obj.get(tag).and_then(Value::as_str).and_then(|t| f.kind.variant(t)) {
                walk(v.fields, obj, prefix, leaf, depth + 1)?;
            }
            continue;
        }
        if let Some(v) = obj.get_mut(f.key) {
            field(f, v, &parent, &join(prefix, f.key), leaf, depth)?;
        }
    }
    Ok(())
}

fn field(f: &Field, v: &mut Value, parent: &Value, path: &str, leaf: &mut dyn Leaf, depth: usize) -> Result<()> {
    if depth > MAX_DEPTH {
        return Ok(());
    }
    if !f.numbers.is_empty() && (v.is_object() || v.is_array()) {
        inner_numbers(f.numbers, v, None, &mut Vec::new(), path, leaf, 0)?;
    }
    let parent = Some(parent);
    match f.kind {
        Kind::Number { dim, .. } => leaf.number(v, dim.resolve(parent, 0), dim.fixed(), path),
        Kind::NumberOrAuto { dim, .. } if v.as_str() != Some("auto") => leaf.number(v, dim.resolve(parent, 0), dim.fixed(), path),
        Kind::PerAxis { dim, .. } | Kind::Vector { dim, .. } => match v {
            Value::Array(items) => {
                for (i, x) in items.iter_mut().enumerate() {
                    leaf.number(x, dim.resolve(parent, i), dim.fixed(), path)?;
                }
                Ok(())
            }
            _ => leaf.number(v, dim.resolve(parent, 0), dim.fixed(), path),
        },
        Kind::Points { dim, .. } => {
            if let Value::Array(rows) = v {
                for row in rows.iter_mut() {
                    if let Value::Array(cols) = row {
                        for (c, x) in cols.iter_mut().enumerate() {
                            leaf.number(x, dim.resolve(parent, c), dim.fixed(), path)?;
                        }
                    }
                }
            }
            Ok(())
        }
        Kind::Group(fields) | Kind::NameOr { fields, .. } => match v {
            Value::Object(o) => walk(fields, o, path, leaf, depth + 1),
            _ => Ok(()),
        },
        Kind::Rows { fields, .. } => {
            if let Value::Array(rows) = v {
                for (i, row) in rows.iter_mut().enumerate() {
                    if let Value::Object(o) = row {
                        walk(fields, o, &format!("{path}/{i}"), leaf, depth + 1)?;
                    }
                }
            }
            Ok(())
        }
        Kind::Variants { tag, .. } => match v {
            Value::Object(o) => match o.get(tag).and_then(Value::as_str).and_then(|t| f.kind.variant(t)) {
                Some(var) => walk(var.fields, o, path, leaf, depth + 1),
                None => Ok(()),
            },
            _ => Ok(()),
        },
        Kind::Either(options) => {
            let holder = v.clone();
            if let Value::Object(o) = v {
                for (k, inner) in o.iter_mut() {
                    if let Some(opt) = options.iter().find(|opt| opt.key == k) {
                        field(opt, inner, &holder, &join(path, k), leaf, depth + 1)?;
                    }
                }
            }
            Ok(())
        }
        _ => Ok(()),
    }
}

/// The numbers of a JSON value `v` (at `at` inside it, under `parent`) by `rules`: every number of
/// a matched value goes to `leaf`.
fn inner_numbers(
    rules: Numbers,
    v: &mut Value,
    parent: Option<&Value>,
    at: &mut Vec<String>,
    path: &str,
    leaf: &mut dyn Leaf,
    depth: usize,
) -> Result<()> {
    if depth > 2 * MAX_DEPTH {
        return Ok(());
    }
    if let Some(dim) = rule_for(rules, at) {
        return each_number(v, dim, parent, 0, path, leaf, depth);
    }
    let holder = if v.is_object() && rules.iter().any(|(_, d)| depends(*d)) { Some(v.clone()) } else { None };
    match v {
        Value::Object(o) => {
            for (k, x) in o.iter_mut() {
                at.push(k.clone());
                inner_numbers(rules, x, holder.as_ref(), at, path, leaf, depth + 1)?;
                at.pop();
            }
        }
        Value::Array(a) => {
            for (i, x) in a.iter_mut().enumerate() {
                at.push(i.to_string());
                inner_numbers(rules, x, None, at, path, leaf, depth + 1)?;
                at.pop();
            }
        }
        _ => {}
    }
    Ok(())
}

/// Whether a dimension's quantity depends on the values beside it.
fn depends(d: Dim) -> bool {
    d != Dim::None && d.quantity().is_none()
}

/// Every number of `v` (a number, a value with its unit, or a list of lists of them) to `leaf`.
fn each_number(v: &mut Value, dim: Dim, parent: Option<&Value>, index: usize, path: &str, leaf: &mut dyn Leaf, depth: usize) -> Result<()> {
    match v {
        Value::Array(items) if depth <= 2 * MAX_DEPTH => {
            for (i, x) in items.iter_mut().enumerate() {
                each_number(x, dim, parent, i, path, leaf, depth + 1)?;
            }
            Ok(())
        }
        _ => leaf.number(v, dim.resolve(parent, index), dim.fixed(), path),
    }
}

/// `params` with every number its form describes converted into `target`: bare numbers are in
/// `input` (or their field's fixed unit), others carry their unit. A value of the wrong
/// dimension, an unknown unit or a number out of range is refused, naming the field.
pub fn convert_params(form: &[Field], params: &Value, input: &DisplaySystem, target: &DisplaySystem) -> Result<Value> {
    let mut out = params.clone();
    if let Value::Object(o) = &mut out {
        walk(form, o, "", &mut ToTarget { input, target }, 0)?;
    }
    Ok(out)
}

/// A definition of the model (stored in SI) as the parameters of `command` that make it, in
/// `display`: what an edit form is filled with and what a panel shows. Fields of a fixed unit
/// stay in it. A command without a form is refused, since its numbers cannot be told apart.
pub fn shown_params(command: &str, params: &Value, display: &DisplaySystem) -> Result<Value> {
    let form = crate::find_command(command)
        .and_then(|spec| spec.form)
        .ok_or_else(|| EngineError::Params(format!("{command} has no form to show its numbers by")))?;
    convert_params(form, params, si(), display)
}

/// `params` as a journal records them when they were typed in `input` but the journal's bare
/// numbers are in another system: every number that follows the display system written with its
/// unit, so replaying reads it as typed.
pub fn with_units(form: &[Field], params: &Value, input: &DisplaySystem) -> Value {
    let mut out = params.clone();
    if let Value::Object(o) = &mut out {
        // Annotating never fails.
        let _ = walk(form, o, "", &mut Annotate { input }, 0);
    }
    out
}

/// The numbers of a command's reply that carry a unit, by command: (command id, (JSON pointer
/// pattern, dimension) list). A `*` segment matches any key or index; a value is a number or a
/// list (of lists) of them. Replies are converted from the command's frame into the display
/// system by these, and get a `units` object naming the unit of each quantity they hold.
pub const REPLIES: &[(&str, &[(&str, Dim)])] = &[
    ("solve.static", STATIC),
    ("solve.nonlinear", NONLINEAR),
    (
        "solve.heat",
        &[
            ("/end_time", Dim::Time),
            ("/min_temperature/value", Dim::Temperature),
            ("/max_temperature/value", Dim::Temperature),
            ("/heat_supplied", Dim::Power),
            ("/heat_removed", Dim::Power),
            ("/heat_flows/loads/*/heat", Dim::Power),
            ("/heat_flows/reactions", Dim::Power),
            ("/heat_flows/imbalance", Dim::Power),
        ],
    ),
    (
        "solve.frequency",
        &[
            ("/total_mass", Dim::Mass),
            ("/modes/*/frequency", Dim::Hertz),
            ("/modes/*/eigenvalue", Dim::Eigenvalue),
            ("/modes/*/effective_mass", Dim::Mass),
        ],
    ),
    (
        "solve.explicit",
        &[
            ("/time", Dim::Time),
            ("/initial_increment", Dim::Time),
            ("/controlling_elements/*/increment", Dim::Time),
            ("/total_mass", Dim::Mass),
            ("/energy/*", Dim::Of(Quantity::Energy)),
            ("/eroded/*/time", Dim::Time),
            ("/eroded/*/internal_energy", Dim::Of(Quantity::Energy)),
            ("/contacts/*/adjusted/*/depth", Dim::Length),
            ("/contacts/*/peak_force", Dim::Force),
            ("/contacts/*/peak_time", Dim::Time),
            ("/contacts/*/first_contact", Dim::Time),
            ("/contacts/*/last_contact", Dim::Time),
            ("/contacts/*/max_penetration", Dim::Length),
            ("/contacts/*/friction", Dim::Of(Quantity::Energy)),
            ("/mass_scaling/target_increment", Dim::Time),
            ("/mass_scaling/added_mass", Dim::Mass),
            ("/mass_scaling/sections/*/1", Dim::Mass),
        ],
    ),
    (
        "solve.fluid",
        &[
            ("/time_step", Dim::Time),
            ("/end_time", Dim::Time),
            ("/velocity_scale", Dim::Velocity),
            ("/frames", Dim::Time),
            ("/max_speed", Dim::Velocity),
            ("/forces/*/final", Dim::Force),
            ("/forces/*/mean", Dim::Force),
            ("/temperature/min", Dim::Temperature),
            ("/temperature/max", Dim::Temperature),
            ("/heat/*/final", Dim::Power),
            ("/heat/*/mean", Dim::Power),
            ("/solids/*/power", Dim::Power),
            ("/solids/*/heat_to_fluid/*", Dim::Power),
            ("/solids/*/temperature/*", Dim::Temperature),
            ("/solids/*/peak", Dim::Temperature),
            ("/solids/*/peak_time", Dim::Time),
            ("/solids/*/limit", Dim::Temperature),
            ("/solids/*/margin", Dim::TemperatureDifference),
            ("/zones/*/flow/*", Dim::VolumeFlow),
            ("/zones/*/pressure_rise/*", Dim::Stress),
            ("/zones/*/measured_pressure_rise/*", Dim::Stress),
            ("/zones/*/pressure_drop/*", Dim::Stress),
            ("/stop/time", Dim::Time),
            ("/stop/window", Dim::Time),
            ("/stop/history/*/time", Dim::Time),
            ("/mass/first", Dim::Of(Quantity::Mass)),
            ("/mass/last", Dim::Of(Quantity::Mass)),
        ],
    ),
    (
        "solve.history",
        &[
            ("/time", Dim::Time),
            ("/increment", Dim::Time),
            ("/energies/*/*", Dim::Of(Quantity::Energy)),
            ("/momentum", Dim::Momentum),
            ("/constraint_force", Dim::Force),
            ("/contact_force", Dim::Force),
            ("/nodes/*/displacement", Dim::Length),
            ("/nodes/*/velocity", Dim::Velocity),
            ("/sets/*/displacement", Dim::Length),
            ("/sets/*/reaction", Dim::Force),
            ("/node/displacement", Dim::Length),
        ],
    ),
    (
        "results.reactions",
        &[
            ("/about", Dim::Length),
            ("/sets/*/force", Dim::Force),
            ("/sets/*/moment", Dim::Moment),
            ("/total/force", Dim::Force),
            ("/total/moment", Dim::Moment),
            ("/applied", Dim::Force),
        ],
    ),
];

/// The lengths of a mesh plan ([`fieldcraft_model::MeshPlan`]) as JSON: `mesh.plan`'s reply and
/// `mesh.plan.apply`'s `plan` object.
pub const MESH_PLAN: Numbers = &[
    ("/size", Dim::Length),
    ("/tolerance", Dim::Length),
    ("/cuts/*/point", Dim::Length),
    ("/rings/*/center", Dim::Length),
    ("/rings/*/radius", Dim::Length),
    ("/rings/*/frame", Dim::Length),
];

/// A datum as the datum commands reply it.
const DATUM: Numbers = &[("/point", Dim::Length)];

/// A mesher's volumes and size, at the top of its reply and per body.
const MESHED: Numbers = &[
    ("/size", Dim::Length),
    ("/volume", Dim::Of(Quantity::Volume)),
    ("/boundary_volume", Dim::Of(Quantity::Volume)),
    ("/surface_volume", Dim::Of(Quantity::Volume)),
    ("/body_volume", Dim::Of(Quantity::Volume)),
    ("/max_gap", Dim::Length),
    ("/bodies/*/size", Dim::Length),
    ("/bodies/*/volume", Dim::Of(Quantity::Volume)),
    ("/bodies/*/boundary_volume", Dim::Of(Quantity::Volume)),
    ("/bodies/*/surface_volume", Dim::Of(Quantity::Volume)),
    ("/bodies/*/grid/root_size", Dim::Length),
    ("/bodies/*/grid/thinnest_wall", Dim::Length),
    ("/bodies/*/grid/narrowest_patch", Dim::Length),
    ("/cells/*/quality/volume", Dim::Of(Quantity::Volume)),
    ("/topology/*/corners", Dim::Length),
];

/// The meshes `mesh.plan.apply` made, under `mesh`, and what it decided.
const PLAN_APPLIED: Numbers = &[
    ("/size", Dim::Length),
    ("/tolerance", Dim::Length),
    ("/pending/*/radius", Dim::Length),
    ("/blend_layers/*/depth", Dim::Length),
    ("/symmetry/*/offset", Dim::Length),
    ("/corners/cuts/*", Dim::Length),
    ("/mesh/size", Dim::Length),
    ("/mesh/volume", Dim::Of(Quantity::Volume)),
    ("/mesh/surface_volume", Dim::Of(Quantity::Volume)),
    ("/mesh/body_volume", Dim::Of(Quantity::Volume)),
    ("/mesh/max_gap", Dim::Length),
    ("/mesh/bodies/*/volume", Dim::Of(Quantity::Volume)),
    ("/mesh/bodies/*/surface_volume", Dim::Of(Quantity::Volume)),
    ("/mesh/cells/*/quality/volume", Dim::Of(Quantity::Volume)),
];

/// A size control's region as JSON, under `prefix`.
const SIZE_CONTROL: Numbers = &[
    ("/size", Dim::Length),
    ("/region/box/min", Dim::Length),
    ("/region/box/max", Dim::Length),
    ("/region/sphere/center", Dim::Length),
    ("/region/sphere/radius", Dim::Length),
    ("/controls/*/size", Dim::Length),
    ("/controls/*/region/box/min", Dim::Length),
    ("/controls/*/region/box/max", Dim::Length),
    ("/controls/*/region/sphere/center", Dim::Length),
    ("/controls/*/region/sphere/radius", Dim::Length),
];

/// The replies of the commands that work in the geometry's millimetres ([`Frame::Geometry`]).
const GEOMETRY_REPLIES: &[(&str, Numbers)] = &[
    (
        "geometry.bodies",
        &[
            ("/*/volume", Dim::Of(Quantity::Volume)),
            ("/*/area", Dim::Area),
            ("/*/centroid", Dim::Length),
            ("/*/bbox", Dim::Length),
            ("/*/edges/*/length", Dim::Length),
            ("/*/vertices/*/point", Dim::Length),
        ],
    ),
    ("geometry.centre", &[("/point", Dim::Length)]),
    (
        "geometry.check",
        &[("/size", Dim::Length), ("/findings/*/value", Dim::Finding), ("/findings/*/limit", Dim::Finding), ("/findings/*/at", Dim::Length)],
    ),
    ("geometry.datum", DATUM),
    ("geometry.datum.edit", DATUM),
    ("geometry.datum.axis", DATUM),
    ("geometry.datum.axis.edit", DATUM),
    ("geometry.datum.point", DATUM),
    ("geometry.datum.point.edit", DATUM),
    ("geometry.defeature.propose", &[("/candidates/*/size", Dim::Length)]),
    ("mesh.quality", &[("/volume", Dim::Of(Quantity::Volume)), ("/area", Dim::Area)]),
    ("mesh.quality.mark", &[("/elements/*/at", Dim::Length)]),
    (
        "mesh.surface.check",
        &[("/volume", Dim::Of(Quantity::Volume)), ("/shells/*/volume", Dim::Of(Quantity::Volume)), ("/defects/listed/*/at", Dim::Length)],
    ),
    ("mesh.tet.generate", MESHED),
    ("mesh.hex.generate", MESHED),
    ("mesh.hex.sweep", MESHED),
    ("mesh.hex.mapped", MESHED),
    ("mesh.hex.blocks", MESHED),
    ("mesh.hex.octree", MESHED),
    ("mesh.hex.project", MESHED),
    ("mesh.cells", MESHED),
    ("mesh.hex.extrude", &[("/volume", Dim::Of(Quantity::Volume)), ("/max_target_gap", Dim::Length)]),
    ("mesh.node.merge", &[("/tolerance", Dim::Length)]),
    ("mesh.check.duplicates", &[("/tolerance", Dim::Length), ("/coincident/*/2", Dim::Length)]),
    (
        "mesh.node.move",
        &[("/max_displacement", Dim::Length), ("/displacement", Dim::Length), ("/sphere/center", Dim::Length), ("/sphere/radius", Dim::Length)],
    ),
    (
        "mesh.node.relax",
        &[("/max_displacement", Dim::Length), ("/sphere/center", Dim::Length), ("/sphere/radius", Dim::Length), ("/rounds/*/reach", Dim::Length)],
    ),
    ("mesh.node.relax.apply", &[("/max_displacement", Dim::Length)]),
    ("mesh.seed.edge", &[("/seeds/*/size", Dim::Length), ("/seeds/*/on/*/length", Dim::Length)]),
    ("mesh.seed.list", &[("/*/size", Dim::Length), ("/*/on/*/length", Dim::Length)]),
    ("mesh.size.set", SIZE_CONTROL),
    ("mesh.size.list", SIZE_CONTROL),
    ("mesh.plan", MESH_PLAN),
    (
        "mesh.plan",
        &[
            ("/decomposition/cuts/*/point", Dim::Length),
            ("/decomposition/tet_region/min", Dim::Length),
            ("/decomposition/tet_region/max", Dim::Length),
        ],
    ),
    ("mesh.plan.apply", PLAN_APPLIED),
];

/// A static step's summary.
const STATIC: &[(&str, Dim)] = &[
    ("/max_displacement/value", Dim::Length),
    ("/max_von_mises/value", Dim::Stress),
    ("/applied_total", Dim::Force),
    ("/reaction_total", Dim::Force),
    ("/reference_points/*/displacement", Dim::Motion),
    ("/reference_points/*/reaction", Dim::Reaction),
    ("/max_temperature/value", Dim::Temperature),
];

/// A nonlinear static step's summary.
const NONLINEAR: &[(&str, Dim)] = &[
    ("/max_displacement/value", Dim::Length),
    ("/max_von_mises/value", Dim::Stress),
    ("/applied_total", Dim::Force),
    ("/reaction_total", Dim::Force),
    ("/reference_points/*/displacement", Dim::Motion),
    ("/reference_points/*/reaction", Dim::Reaction),
    ("/contacts/*/normal_force", Dim::Force),
    ("/contacts/*/force", Dim::Force),
    ("/contacts/*/max_pressure/value", Dim::Stress),
    ("/contacts/*/penetration", Dim::Length),
];

/// Whether JSON pointer `pattern` (`*` for any segment, empty for the root) names `path`.
fn matches(pattern: &str, path: &[String]) -> bool {
    let mut segments = pattern.split('/').filter(|s| !s.is_empty());
    let mut count = 0;
    // The path first: a pattern segment is taken only when a path segment is there for it.
    let all = path.iter().zip(segments.by_ref()).all(|(p, s)| {
        count += 1;
        s == "*" || s == p
    });
    all && count == path.len() && segments.next().is_none()
}

/// The dimension of the value at `path` (pointer segments) by `rules`.
fn rule_for(rules: &[(&'static str, Dim)], path: &[String]) -> Option<Dim> {
    rules.iter().find_map(|(pattern, dim)| matches(pattern, path).then_some(*dim))
}

/// Convert `reply` of command `id` from `from` (the command's frame) to `to` (the display
/// system), and name the units it is in.
pub fn convert_reply(id: &str, reply: &mut Value, from: &DisplaySystem, to: &DisplaySystem) -> Result<()> {
    let rules: Vec<(&'static str, Dim)> =
        REPLIES.iter().chain(GEOMETRY_REPLIES).filter(|(c, _)| *c == id).flat_map(|(_, rules)| rules.iter().copied()).collect();
    if rules.is_empty() {
        return Ok(());
    }
    let mut units = Map::new();
    let mut path = Vec::new();
    reply_values(reply, None, &mut path, &rules, from, to, &mut units)?;
    if let Value::Object(o) = reply
        && !units.is_empty()
        && !o.contains_key("units")
    {
        o.insert("units".into(), Value::Object(units));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn reply_values(
    v: &mut Value,
    parent: Option<&Value>,
    path: &mut Vec<String>,
    rules: &[(&'static str, Dim)],
    from: &DisplaySystem,
    to: &DisplaySystem,
    units: &mut Map<String, Value>,
) -> Result<()> {
    if path.len() > 2 * MAX_DEPTH {
        return Ok(());
    }
    if let Some(dim) = rule_for(rules, path) {
        return express_all(v, dim, parent, 0, from, to, units);
    }
    // The holder of the values below, for a dimension that depends on its siblings.
    let holder = if v.is_object() && rules.iter().any(|(_, d)| depends(*d)) { Some(v.clone()) } else { None };
    match v {
        Value::Object(o) => {
            for (k, x) in o.iter_mut() {
                path.push(k.clone());
                reply_values(x, holder.as_ref(), path, rules, from, to, units)?;
                path.pop();
            }
        }
        Value::Array(a) => {
            for (i, x) in a.iter_mut().enumerate() {
                path.push(i.to_string());
                reply_values(x, None, path, rules, from, to, units)?;
                path.pop();
            }
        }
        _ => {}
    }
    Ok(())
}

/// A reply number (or list) of dimension `dim`, from `from` into `to`.
fn express_all(
    v: &mut Value,
    dim: Dim,
    parent: Option<&Value>,
    index: usize,
    from: &DisplaySystem,
    to: &DisplaySystem,
    units: &mut Map<String, Value>,
) -> Result<()> {
    match v {
        Value::Array(a) => {
            for (i, x) in a.iter_mut().enumerate() {
                express_all(x, dim, parent, i, from, to, units)?;
            }
        }
        Value::Number(n) => {
            let Some(q) = dim.resolve(parent, index) else { return Ok(()) };
            let x = n.as_f64().unwrap_or(f64::NAN);
            let (a, b) = (from.unit_ref(q), to.unit_ref(q));
            let y =
                if same_unit(a, b) { x } else { q.to_si(x, a).and_then(|si| q.express(si, b)).map_err(|e| EngineError::Params(e.to_string()))? };
            if let Some(n) = serde_json::Number::from_f64(y) {
                *v = Value::Number(n);
            }
            if b.symbol() != "1" {
                units.insert(q.name().to_string(), json!(b.symbol()));
            }
        }
        _ => {}
    }
    Ok(())
}

/// How [`show_by_key`] shows a value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Show {
    /// An SI value of a quantity, in the display system's unit.
    Quantity(Quantity),
    /// A value without a unit, times a factor (a mass-normalised mode shape).
    Factor(f64),
}

/// Show every number of `v` (SI) in `display` by the key it is under: `rule(key, len, index)`
/// for entry `index` of a vector of `len` (a lone number is entry 0 of 1; a list of vectors, such
/// as a history, takes its vectors' length). Records each key's unit in `units`.
pub fn show_by_key(v: &mut Value, display: &DisplaySystem, rule: &dyn Fn(&str, usize, usize) -> Option<Show>, units: &mut Map<String, Value>) {
    fn number(x: &mut Value, show: Show, display: &DisplaySystem) {
        if let Some(n) = x.as_f64() {
            let y = match show {
                Show::Quantity(q) => display.express(q, n).unwrap_or(f64::NAN),
                Show::Factor(k) => n * k,
            };
            if let Some(y) = serde_json::Number::from_f64(y) {
                *x = Value::Number(y);
            }
        }
    }
    fn values(
        key: &str,
        x: &mut Value,
        display: &DisplaySystem,
        rule: &dyn Fn(&str, usize, usize) -> Option<Show>,
        units: &mut Map<String, Value>,
        depth: usize,
    ) {
        if depth > 2 * MAX_DEPTH {
            return;
        }
        match x {
            Value::Number(_) => {
                if let Some(show) = rule(key, 1, 0) {
                    note(key, show, display, units);
                    number(x, show, display);
                }
            }
            Value::Array(items) if items.iter().all(|i| i.is_number() || i.is_null()) => {
                let len = items.len();
                for (i, item) in items.iter_mut().enumerate() {
                    if let Some(show) = rule(key, len, i) {
                        note(key, show, display, units);
                        number(item, show, display);
                    }
                }
            }
            Value::Array(items) => {
                for item in items {
                    if item.is_object() {
                        walk_keys(item, display, rule, units, depth + 1);
                    } else {
                        values(key, item, display, rule, units, depth + 1);
                    }
                }
            }
            Value::Object(_) => walk_keys(x, display, rule, units, depth + 1),
            _ => {}
        }
    }
    fn note(key: &str, show: Show, display: &DisplaySystem, units: &mut Map<String, Value>) {
        let unit = match show {
            Show::Quantity(q) => Some(display.symbol(q)).filter(|u| *u != "1"),
            Show::Factor(_) => None,
        };
        units.entry(key.to_string()).or_insert_with(|| unit.map_or(Value::Null, |u| json!(u)));
    }
    fn walk_keys(
        v: &mut Value,
        display: &DisplaySystem,
        rule: &dyn Fn(&str, usize, usize) -> Option<Show>,
        units: &mut Map<String, Value>,
        depth: usize,
    ) {
        if let Value::Object(o) = v {
            for (k, x) in o.iter_mut() {
                values(k, x, display, rule, units, depth);
            }
        }
    }
    walk_keys(v, display, rule, units, 0);
}

#[cfg(test)]
mod geometry_tests;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::schema::*;

    const FORM: &[Field] = &[
        num("young", "E", Dim::Stress),
        num("angle", "Angle", Dim::Degree),
        free_vector("values", "Values", Dim::Motion, &["u1", "u2", "u3", "ur1", "ur2", "ur3"]),
        num("t", "T", Dim::Temperature),
        num("ratio", "Ratio", Dim::None),
    ];

    #[test]
    fn bare_numbers_are_in_the_input_system_and_units_override() {
        let mm = DisplaySystem::default();
        let p = json!({"young": 210000, "angle": 5, "values": [1, null, 2, 0.5, null, null], "t": "20 °C", "ratio": 0.3});
        let si_params = convert_params(FORM, &p, &mm, si()).unwrap();
        assert_eq!(si_params, json!({"young": 2.1e11, "angle": 5, "values": [0.001, null, 0.002, 0.5, null, null], "t": 293.15, "ratio": 0.3}));
        // Into millimetres a millimetre stays as typed.
        let framed =
            convert_params(FORM, &json!({"values": [3.97, 0, 0, 0, 0, 0], "young": {"value": 210, "unit": "GPa"}}), &mm, geometry()).unwrap();
        assert_eq!(framed, json!({"values": [3.97, 0.0, 0.0, 0.0, 0.0, 0.0], "young": 210000.0}));
        // A metre display: bare lengths are metres; the fixed degrees stay degrees, and a
        // fixed-unit field reaches its command in its unit whatever it is given in.
        let m = DisplaySystem::new(Preset::MKgS);
        assert_eq!(
            convert_params(FORM, &json!({"angle": 5, "values": [1, 0, 0, 0, 0, 0]}), &m, si()).unwrap(),
            json!({"angle": 5, "values": [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]})
        );
        assert_eq!(convert_params(FORM, &json!({"angle": "90 deg"}), &m, si()).unwrap(), json!({"angle": 90.0}));
    }

    /// The geometry frame's absolute amounts are written in units the table reads, and are the
    /// same number in the frame.
    #[test]
    fn geometry_amounts_read_in_millimetres() {
        for (q, t) in [
            (Quantity::Length, MIN_LENGTH),
            (Quantity::Area, fieldcraft_model::units::Typed::new(1e-12, "mm^2")),
            (Quantity::Volume, fieldcraft_model::units::Typed::new(1e-12, "mm^3")),
        ] {
            assert!(t.express(q, geometry()).is_ok(), "{t:?}");
            assert_eq!(geometry_amount(q, t).to_bits(), t.value.to_bits(), "{t:?}");
        }
    }

    #[test]
    fn hostile_units_are_refused_by_field() {
        let mm = DisplaySystem::default();
        for (p, field) in [
            (json!({"young": "210 mm"}), "young"),
            (json!({"young": "210 parsecs"}), "young"),
            (json!({"young": {"value": 1, "unit": "Pa", "x": 1}}), "young"),
            (json!({"young": "1e400 Pa"}), "young"),
            (json!({"young": 1e308}), "young"),
            (json!({"t": "-300 °C"}), "t"),
            (json!({"t": -300}), "t"),
            (json!({"values": ["5 MPa", 0, 0, 0, 0, 0]}), "values"),
            (json!({"values": [0, 0, 0, "1 m", 0, 0]}), "values"),
            (json!({"angle": "x".repeat(10_000)}), "angle"),
        ] {
            let e = convert_params(FORM, &p, &mm, si()).unwrap_err();
            assert!(matches!(&e, EngineError::Param(problem) if problem.path == field), "{p}: {e:?}");
        }
        // Not a number at all is left for the command.
        assert_eq!(convert_params(FORM, &json!({"young": true}), &mm, si()).unwrap(), json!({"young": true}));
    }

    #[test]
    fn a_journal_entry_typed_in_another_system_keeps_its_units() {
        let m = DisplaySystem::new(Preset::MKgS);
        let p = json!({"young": 2.1e11, "angle": 5, "values": [0.001, null, 0, 0.5, null, null], "ratio": 0.3});
        let annotated = with_units(FORM, &p, &m);
        assert_eq!(
            annotated,
            json!({"young": {"value": 2.1e11, "unit": "Pa"}, "angle": 5, "values": [{"value": 0.001, "unit": "m"}, null, {"value": 0.0, "unit": "m"}, {"value": 0.5, "unit": "rad"}, null, null], "ratio": 0.3})
        );
        // Read in any other system, it is what was typed.
        assert_eq!(convert_params(FORM, &annotated, &DisplaySystem::default(), si()).unwrap(), convert_params(FORM, &p, &m, si()).unwrap());
    }
}
