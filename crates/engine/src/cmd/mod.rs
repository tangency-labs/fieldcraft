//! The command registry. Each area has its own module; [`ALL`] lists every command.

mod block;
mod cells;
mod contact;
pub(crate) mod datum;
pub(crate) mod deck;
mod defeature;
pub(crate) mod define;
pub mod diagnostics;
mod edit;
mod extrude;
pub(crate) mod file;
pub(crate) mod forms;
mod geometry;
mod hex;
mod initial;
pub(crate) mod mesh;

pub use datum::edit_params as datum_edit_params;
pub use datum::preview as datum_preview;
pub use fieldcraft_mesh::relax::{ElementChange, Moved as NodeMoved};
pub use mesh::analytic_surfaces;
pub use node_move::{
    Capped as NodeMoveCapped, Prepared as NodeMovePrepared, nearest_node as node_move_nearest, prepare as node_move_prepare,
    prepare_capped as node_move_prepare_capped,
};
pub use nodes::merge_preview as node_merge_preview;
pub use plan::untagged_faces;
pub use relax::{Preview as RelaxPreview, preview as relax_preview, sphere as relax_sphere};
pub use solve::{FluidEstimate, fluid_estimate};
mod mirror;
mod model;
pub(crate) mod names;
pub(crate) mod node_move;
pub(crate) mod nodes;
pub(crate) mod octree;
mod orphan;
mod path;
pub(crate) mod plan;
mod project;
pub(crate) mod relax;
mod render;
pub mod replace;
mod results;
pub(crate) mod seeds;
pub(crate) mod sizes;
pub(crate) mod solve;
pub(crate) mod stitch;
mod symmetry;
mod tags;
pub mod update;
pub(crate) mod view;
pub(crate) mod views;
mod virtual_faces;
pub mod xy;

use crate::schema::ProblemKind;
use fieldcraft_model::MAX_LABEL;
use serde::Deserialize;
use serde_json::{Value, json};

use fieldcraft_io::roots::Access;

use crate::{CommandSpec, EngineError, Result, Session};

pub static ALL: &[CommandSpec] = &[
    model::INFO,
    geometry::OPEN,
    geometry::BODIES,
    geometry::FEATURES,
    geometry::CENTRE,
    geometry::CHECK,
    geometry::TESSELLATE,
    geometry::IMPORT_MESH,
    datum::DATUM_PLANE,
    datum::DATUM_AXIS,
    datum::DATUM_POINT,
    datum::DATUM_EDIT,
    datum::DATUM_AXIS_EDIT,
    datum::DATUM_POINT_EDIT,
    datum::DATUM_DELETE,
    view::HIDE,
    view::SHOW,
    view::ISOLATE,
    view::SHOW_ALL,
    views::SAVE,
    views::RESTORE,
    views::DELETE,
    views::LIST,
    datum::PARTITION,
    datum::DATUM_SUPPRESS,
    datum::PARTITION_SUPPRESS,
    datum::PARTITION_DELETE,
    virtual_faces::FACES_MERGE,
    virtual_faces::FACES_MERGE_DELETE,
    virtual_faces::FACE_ABSORB,
    virtual_faces::FACE_SPLIT,
    virtual_faces::FACE_SPLIT_DELETE,
    defeature::DEFEATURE,
    defeature::DEFEATURE_DELETE,
    defeature::DEFEATURE_PROPOSE,
    stitch::STITCH,
    stitch::STITCH_DELETE,
    mesh::NODE_ADD,
    mesh::IMPORT_ORPHAN,
    mesh::NODE_REMOVE,
    mesh::ELEMENT_ADD,
    mesh::ELEMENT_REMOVE,
    replace::CLEAR,
    mesh::QUALITY,
    mesh::QUALITY_MARK,
    mesh::EXPORT_STL,
    mesh::DELAUNAY_POINTS,
    mesh::SURFACE_CHECK,
    mesh::PLANAR,
    mesh::TET_GENERATE,
    octree::HEX_OCTREE,
    block::BLOCK,
    hex::HEX_GENERATE,
    hex::HEX_SWEEP,
    hex::HEX_MAPPED,
    hex::HEX_BLOCKS,
    project::HEX_PROJECT,
    extrude::HEX_EXTRUDE,
    cells::CELLS,
    plan::PLAN,
    plan::APPLY,
    seeds::SEED_EDGE,
    seeds::SEED_CLEAR,
    seeds::SEED_LIST,
    sizes::SIZE_SET,
    sizes::SIZE_GROWTH,
    sizes::SIZE_CLEAR,
    sizes::SIZE_LIST,
    mirror::MIRROR,
    nodes::MERGE,
    nodes::DUPLICATES,
    nodes::RENUMBER,
    relax::RELAX,
    relax::APPLY,
    node_move::MOVE,
    tags::TAG,
    tags::TAG_CLEAR,
    model::UNITS,
    model::CHECK,
    model::SET_RESOLVE,
    define::SET_CREATE,
    define::SET_EDIT,
    define::SET_DELETE,
    define::MATERIAL_CREATE,
    define::MATERIAL_EDIT,
    define::MATERIAL_DELETE,
    define::MATERIAL_PLASTICITY,
    define::MATERIAL_DAMAGE,
    define::MATERIAL_PVT,
    define::MATERIAL_SPECIFIC_VOLUME,
    define::MATERIAL_PHASE_CHANGE,
    define::SECTION_CREATE,
    define::SECTION_EDIT,
    define::SECTION_DELETE,
    define::REFERENCE_POINT_CREATE,
    define::REFERENCE_POINT_EDIT,
    define::REFERENCE_POINT_DELETE,
    define::COUPLING_CREATE,
    define::COUPLING_EDIT,
    define::COUPLING_DELETE,
    contact::PAIR,
    contact::GENERAL,
    contact::EDIT,
    contact::DELETE,
    contact::FIND,
    contact::TIE,
    contact::TIE_EDIT,
    contact::TIE_DELETE,
    contact::RIGID_WALL,
    contact::RIGID_WALL_DELETE,
    define::MASS_CREATE,
    define::MASS_EDIT,
    define::MASS_DELETE,
    define::SPRING_CREATE,
    define::SPRING_EDIT,
    define::SPRING_DELETE,
    define::BOUNDARY_CONDITION_CREATE,
    define::BOUNDARY_CONDITION_EDIT,
    define::BOUNDARY_CONDITION_DELETE,
    symmetry::SYMMETRY,
    define::LOAD_CREATE,
    define::LOAD_EDIT,
    define::LOAD_DELETE,
    initial::MAP,
    initial::DELETE,
    define::FLUID_DOMAIN_CREATE,
    define::FLUID_DOMAIN_EDIT,
    define::FLUID_DOMAIN_DELETE,
    define::FLUID_BOUNDARY,
    define::FLUID_THERMAL_BOUNDARY,
    define::FLUID_COMPONENT_LIMIT,
    define::STEP_CREATE,
    define::STEP_EDIT,
    define::STEP_DELETE,
    define::STEP_ACTIVATE,
    define::RENAME,
    define::ACTIVITY,
    file::NEW,
    file::OPEN,
    file::SAVE,
    file::SAVE_AS,
    file::RECOVERY_LIST,
    file::RECOVER,
    file::RECOVERY_DISCARD,
    file::SAMPLE,
    edit::UNDO,
    edit::REDO,
    deck::WRITE,
    deck::READ,
    solve::RUN,
    solve::STATIC,
    solve::NONLINEAR,
    solve::FREQUENCY,
    solve::BUCKLE,
    solve::HEAT,
    solve::EXPLICIT,
    solve::HISTORY,
    solve::SETTINGS,
    solve::FLUID,
    solve::PROBE,
    results::DEFORMED,
    results::CONTOUR,
    results::MODES,
    results::ANIMATE,
    results::FRAME,
    results::SLICE,
    results::SURFACE,
    results::STREAMLINES,
    results::REACTIONS,
    results::EXPORT_CSV,
    results::EXPORT_VTU,
    render::SAVE,
    render::DELETE,
    render::EXPORT,
    results::CLEAR,
    xy::XY,
    path::PATH,
    diagnostics::CRASH_REPORTS,
    diagnostics::CRASH_REPORT,
    diagnostics::CRASH_REPORTS_DELETE,
    diagnostics::BUG_REPORT,
    diagnostics::CRASH_TEST,
    update::CHECK,
];

/// The problem with required parameter `key`: [`ProblemKind::Required`] when it is left out,
/// else `kind`.
pub(crate) fn wrong(p: &Value, key: &str, kind: ProblemKind) -> EngineError {
    EngineError::param(key, if p.get(key).is_none_or(Value::is_null) { ProblemKind::Required } else { kind })
}

/// A required label parameter (1..=u32::MAX; the model checks its own range).
pub(crate) fn label(p: &Value, key: &str) -> Result<u32> {
    p.get(key)
        .and_then(Value::as_u64)
        .and_then(|v| u32::try_from(v).ok())
        .ok_or_else(|| wrong(p, key, ProblemKind::NotAWholeNumber { min: 1, max: i64::from(u32::MAX) }))
}

/// The first node and element labels after the model's largest, for `nodes` new nodes and
/// `elements` new elements, all within [`MAX_LABEL`].
pub(crate) fn fresh_labels(model: &fieldcraft_model::Model, nodes: usize, elements: usize) -> Result<(u32, u32)> {
    let first = |last: Option<&u32>, count: usize| {
        let first = last.map_or(Some(1), |l| l.checked_add(1))?;
        (u64::from(first) + count as u64 <= u64::from(MAX_LABEL) + 1).then_some(first)
    };
    match (first(model.nodes().keys().next_back(), nodes), first(model.elements().keys().next_back(), elements)) {
        (Some(n), Some(e)) => Ok((n, e)),
        _ => Err(EngineError::Params(format!("the mesh needs labels beyond {MAX_LABEL}"))),
    }
}

/// A required `[x, y, z]` parameter.
pub(crate) fn xyz(p: &Value, key: &str) -> Result<[f64; 3]> {
    p.get(key).and_then(xyz_value).ok_or_else(|| wrong(p, key, ProblemKind::NotAVector { len: 3 }))
}

/// A JSON `[x, y, z]`.
fn xyz_value(v: &Value) -> Option<[f64; 3]> {
    match v.as_array()?.as_slice() {
        [x, y, z] => Some([x.as_f64()?, y.as_f64()?, z.as_f64()?]),
        _ => None,
    }
}

/// A required list of `[x, y, z]` points, at most `max_len` long.
pub(crate) fn points(p: &Value, key: &str, max_len: usize) -> Result<Vec<[f64; 3]>> {
    let raw = p.get(key).and_then(Value::as_array).filter(|a| a.len() <= max_len);
    let raw = raw.ok_or_else(|| wrong(p, key, ProblemKind::NotAList { max: max_len, what: "[x, y, z]" }))?;
    raw.iter()
        .enumerate()
        .map(|(i, v)| xyz_value(v).ok_or_else(|| EngineError::invalid(key, format!("must be a list of [x, y, z] (entry {i} is not)"))))
        .collect()
}

/// A required list of labels (each in 1..=MAX_LABEL), at most `max_len` long.
pub(crate) fn labels(p: &Value, key: &str, max_len: usize) -> Result<Vec<u32>> {
    let raw = p
        .get(key)
        .and_then(Value::as_array)
        .filter(|a| a.len() <= max_len)
        .ok_or_else(|| wrong(p, key, ProblemKind::NotAList { max: max_len, what: "labels" }))?;
    raw.iter()
        .map(|v| v.as_u64().and_then(|n| u32::try_from(n).ok()).filter(|n| (1..=MAX_LABEL).contains(n)))
        .collect::<Option<Vec<u32>>>()
        .ok_or_else(|| EngineError::invalid(key, format!("must be labels in 1..={MAX_LABEL}")))
}

/// A required path parameter, as given.
fn path_param<'a>(p: &'a Value, key: &str) -> Result<&'a str> {
    p.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty() && s.len() <= fieldcraft_io::MAX_PATH_LEN && !s.contains('\0'))
        .ok_or_else(|| wrong(p, key, ProblemKind::NotAPath))
}

/// A required path parameter naming a file (or folder) the command reads or writes: the path
/// to use, checked against the session's [`Session::file_access`] (as given while unconfined;
/// resolved, links followed, when confined). Every command that touches a file named by a
/// parameter takes its path from here.
pub(crate) fn file_param(s: &Session, p: &Value, key: &str, access: Access) -> Result<String> {
    s.checked_path(path_param(p, key)?, access)
}

/// An optional number parameter.
pub(crate) fn optional_number(p: &Value, key: &str, default: f64) -> Result<f64> {
    match p.get(key) {
        None => Ok(default),
        Some(v) => v.as_f64().ok_or_else(|| EngineError::param(key, ProblemKind::NotANumber)),
    }
}

/// The `name` parameter.
pub(crate) fn name(p: &Value) -> Result<&str> {
    p.get("name").and_then(Value::as_str).ok_or_else(|| wrong(p, "name", ProblemKind::NotAName))
}

/// An optional parameter read as `T`'s JSON form; `None` when it is absent.
pub(crate) fn opt<'a, T: Deserialize<'a>>(p: &'a Value, key: &str) -> Result<Option<T>> {
    p.get(key).map(|v| T::deserialize(v).map_err(|e| EngineError::invalid(key, format!("is not valid: {e}")))).transpose()
}

/// A required parameter read as `T`'s JSON form.
pub(crate) fn req<'a, T: Deserialize<'a>>(p: &'a Value, key: &str) -> Result<T> {
    opt(p, key)?.ok_or_else(|| EngineError::param(key, ProblemKind::Required))
}

/// Entries listed in a report: the first ones, and how many there are.
pub(crate) const LISTED: usize = 100;

pub(crate) fn listed<T: serde::Serialize>(v: &[T]) -> Value {
    json!({"count": v.len(), "labels": v.get(..LISTED).unwrap_or(v)})
}
