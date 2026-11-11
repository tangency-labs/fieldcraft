//! Linear history import for the procedures emitted by the Abaqus writer.
//!
//! This module mirrors only the writer's current static, frequency and buckling records. It does
//! not interpret amplitudes, distributed loads, constraints, discrete elements, contacts or raw
//! source text; those records remain in `Read::unknown` for their own changes.

use super::{Block, Limits, Read, canonical, check_work, error, label, lex, number, quantity, row_fields};
use crate::{Error, Result};
use fieldcraft_model::units::{DisplaySystem, Preset, Quantity};
use fieldcraft_model::{BoundaryCondition, Distribution, Load, Model, Procedure, Set, Step, Symmetry, Target};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Family {
    Boundary,
    Force,
}

#[derive(Clone, Copy, PartialEq)]
enum BoundaryKind {
    Displacement,
    /// The input's `PINNED` shorthand.  The model representation is still a zero
    /// displacement on DOFs 1--3, but history must keep this spelling distinct
    /// from a direct DOF row so that changes require `OP=NEW`.
    Pinned,
    Fixed,
    Symmetry(Symmetry),
}

#[derive(Clone)]
struct Driven {
    value: f64,
    line: usize,
    boundary: Option<BoundaryKind>,
    members: usize,
    initial: Option<usize>,
}

type StateKey = (Family, String, u8);
type State = BTreeMap<StateKey, Driven>;

#[derive(Clone, Copy, Default)]
struct BoundaryOverlap {
    target: u32,
    nonzero: bool,
    multiple_zero: bool,
}

type BoundaryNodes = BTreeMap<(u32, u8), BoundaryOverlap>;

#[derive(Default)]
struct TargetIds {
    names: BTreeMap<String, u32>,
}

impl TargetIds {
    fn id(&mut self, target: &str) -> u32 {
        let next = self.names.len() as u32;
        *self.names.entry(target.to_string()).or_insert(next)
    }
}

#[derive(Clone, Default)]
struct HistoryState {
    values: State,
    /// Active boundary overlap information by resolved node/DOF. The target ID is assigned once
    /// per spelling; only the conflict facts needed by the next block are retained here.
    boundary_nodes: Arc<BoundaryNodes>,
    members: usize,
}

impl HistoryState {
    fn boundary_only(&self) -> Self {
        let values: State =
            self.values.iter().filter(|((family, _, _), _)| *family == Family::Boundary).map(|(key, driven)| (key.clone(), driven.clone())).collect();
        let members = values.values().map(|driven| driven.members).sum();
        Self { values, boundary_nodes: Arc::clone(&self.boundary_nodes), members }
    }
}

fn unsupported(line: usize, message: impl Into<String>) -> Error {
    Error::Unsupported { format: "Abaqus analysis import", what: format!("line {line}: {}", message.into()) }
}

fn allow(block: &Block<'_>, keys: &[&str]) -> Result<()> {
    if let Some(key) = block.params.keys().find(|key| !keys.contains(&key.as_str())) {
        return Err(unsupported(block.line, format!("*{} parameter {key}", block.keyword)));
    }
    Ok(())
}

fn value<'a>(block: &'a Block<'_>, key: &str, default: &'a str) -> Result<&'a str> {
    match block.params.get(key) {
        None => Ok(default),
        Some(Some(value)) if !value.is_empty() => Ok(value),
        _ => Err(error(block.line, format!("{key} requires a value"))),
    }
}

fn target_name(value: &str, line: usize) -> Result<String> {
    let value = lex::unquote(value, line)?;
    if let Ok(label) = value.parse::<u32>() {
        return Ok(label.to_string());
    }
    Ok(value.to_ascii_uppercase())
}

fn actual_set_name(model: &Model, input_sets: &BTreeSet<String>, target: &str, line: usize) -> Result<String> {
    let identity = target.to_ascii_uppercase();
    if !input_sets.contains(&identity) {
        return Err(error(line, format!("missing node set {target}")));
    }
    model.sets().keys().find(|current| current.eq_ignore_ascii_case(target)).cloned().ok_or_else(|| error(line, format!("missing node set {target}")))
}

fn members(model: &Model, input_sets: &BTreeSet<String>, target: &str, line: usize) -> Result<Vec<u32>> {
    let target = lex::unquote(target, line)?;
    if let Ok(label) = target.parse::<u32>() {
        if model.node(label).is_some() {
            return Ok(vec![label]);
        }
        return Err(error(line, format!("missing node label {label}")));
    }
    let set_name = actual_set_name(model, input_sets, &target, line)?;
    let mut nodes = model.set(&set_name).map(|set| set.nodes.clone()).ok_or_else(|| error(line, format!("missing node set {target}")))?;
    nodes.sort_unstable();
    nodes.dedup();
    if nodes.is_empty() {
        return Err(error(line, format!("target {target} has no nodes")));
    }
    Ok(nodes)
}

fn remove_family(state: &mut HistoryState, family: Family) {
    let removed: usize = state.values.iter().filter(|((current, _, _), _)| *current == family).map(|(_, driven)| driven.members).sum();
    if family == Family::Boundary {
        Arc::make_mut(&mut state.boundary_nodes).clear();
    }
    state.values.retain(|(current, _, _), _| *current != family);
    state.members = state.members.saturating_sub(removed);
}

fn insert(state: &mut HistoryState, key: StateKey, driven: Driven, target_nodes: &[u32], limits: Limits, line: usize) -> Result<bool> {
    let new = !state.values.contains_key(&key);
    if new && state.members.saturating_add(target_nodes.len()) > limits.members {
        return Err(error(line, "analysis membership limit reached"));
    }
    if new {
        state.members = state.members.saturating_add(target_nodes.len());
    }
    state.values.insert(key, driven);
    Ok(new)
}

fn apply(
    block: &Block<'_>,
    state: &mut HistoryState,
    model: &Model,
    input_sets: &BTreeSet<String>,
    limits: Limits,
    seen: &mut BTreeSet<StateKey>,
    target_ids: &mut TargetIds,
    clear_family: bool,
    initial: Option<usize>,
    frequency: bool,
    transformed_nodes: &BTreeSet<u32>,
) -> Result<()> {
    let boundary = block.keyword == "BOUNDARY";
    allow(block, if boundary { &["OP", "TYPE"][..] } else { &["OP"][..] })?;
    if boundary && canonical(value(block, "TYPE", "DISPLACEMENT")?) != "DISPLACEMENT" {
        return Err(unsupported(block.line, "boundary velocity and acceleration are not mapped"));
    }
    let op = canonical(value(block, "OP", "MOD")?);
    if !matches!(op.as_str(), "MOD" | "NEW") {
        return Err(error(block.line, "OP must be NEW or MOD"));
    }
    if clear_family {
        remove_family(state, if boundary { Family::Boundary } else { Family::Force });
    }
    let mut work = 0usize;
    for row in &block.rows {
        check_work(work)?;
        work = work.saturating_add(1);
        let fields = row_fields(row.text, row.line)?;
        if boundary {
            let target_text = fields.first().map(String::as_str).unwrap_or_default();
            if !(2..=4).contains(&fields.len()) {
                return Err(error(row.line, "BOUNDARY expects target, DOF/type and optional magnitude"));
            }
            let kind = canonical(fields.get(1).map(String::as_str).unwrap_or_default());
            let (dofs, value, boundary_kind) = if let Ok(first) = kind.parse::<u8>() {
                let last = fields
                    .get(2)
                    .filter(|value| !value.is_empty())
                    .map(|value| value.parse::<u8>().map_err(|_| error(row.line, "invalid boundary last DOF")))
                    .transpose()?
                    .unwrap_or(first);
                if !(1..=6).contains(&first) || last < first || last > 6 {
                    return Err(unsupported(row.line, "only structural DOFs 1..6 are mapped"));
                }
                let value = fields.get(3).filter(|value| !value.is_empty()).map(|value| number(value, row.line)).transpose()?.unwrap_or(0.0);
                ((first..=last).collect::<Vec<_>>(), value, BoundaryKind::Displacement)
            } else {
                if fields.len() != 2 {
                    return Err(error(row.line, "a named boundary has no magnitude"));
                }
                let boundary_kind = match kind.as_str() {
                    "ENCASTRE" => BoundaryKind::Fixed,
                    "XSYMM" => BoundaryKind::Symmetry(Symmetry::X),
                    "YSYMM" => BoundaryKind::Symmetry(Symmetry::Y),
                    "ZSYMM" => BoundaryKind::Symmetry(Symmetry::Z),
                    "PINNED" => BoundaryKind::Pinned,
                    _ => return Err(unsupported(row.line, format!("boundary type {kind}"))),
                };
                let dofs = match boundary_kind {
                    BoundaryKind::Fixed => vec![1, 2, 3, 4, 5, 6],
                    BoundaryKind::Symmetry(plane) => plane.dofs().list(),
                    BoundaryKind::Displacement => vec![1, 2, 3],
                    BoundaryKind::Pinned => vec![1, 2, 3],
                };
                (dofs, 0.0, boundary_kind)
            };
            let target = target_name(target_text, row.line)?;
            let target_id = target_ids.id(&target);
            let target_nodes = members(model, input_sets, target_text, row.line)?;
            for dof in dofs {
                let key = (Family::Boundary, target.clone(), dof);
                if !seen.insert(key.clone()) {
                    return Err(unsupported(row.line, "duplicate boundary definition for the same target and DOF in one step"));
                }
                for node in &target_nodes {
                    if let Some(old) = state.boundary_nodes.get(&(*node, dof)) {
                        let overlap_value = if frequency { 0.0 } else { value };
                        let conflict = old.multiple_zero && overlap_value != 0.0
                            || old.target != target_id && (old.nonzero && !frequency || overlap_value != 0.0);
                        if conflict {
                            return Err(unsupported(row.line, "boundary targets overlap with a nonzero value"));
                        }
                    }
                }
                if let Some(previous) = state.values.get(&key)
                    && previous.boundary != Some(boundary_kind)
                {
                    let message = match (previous.boundary, Some(boundary_kind)) {
                        (Some(BoundaryKind::Pinned), Some(BoundaryKind::Displacement))
                        | (Some(BoundaryKind::Displacement), Some(BoundaryKind::Pinned)) => {
                            "type-format and direct-format boundary changes require OP=NEW"
                        }
                        _ => "boundary type changes require OP=NEW",
                    };
                    return Err(unsupported(row.line, message));
                }
                let driven = Driven { value, line: row.line, boundary: Some(boundary_kind), members: target_nodes.len(), initial };
                insert(state, key, driven, &target_nodes, limits, row.line)?;
                if target_nodes.iter().any(|node| transformed_nodes.contains(node)) {
                    return Err(unsupported(row.line, "history target is covered by a *TRANSFORM local coordinate system"));
                }
                for node in &target_nodes {
                    check_work(work)?;
                    work = work.saturating_add(1);
                    let nodes = Arc::make_mut(&mut state.boundary_nodes);
                    if let Some(overlap) = nodes.get_mut(&(*node, dof)) {
                        if overlap.target == target_id {
                            overlap.nonzero = if frequency { false } else { value != 0.0 };
                        } else {
                            overlap.multiple_zero = true;
                            overlap.target = target_id;
                        }
                    } else {
                        nodes.insert(
                            (*node, dof),
                            BoundaryOverlap { target: target_id, nonzero: if frequency { false } else { value != 0.0 }, multiple_zero: false },
                        );
                    }
                }
            }
        } else {
            if fields.len() != 3 {
                return Err(error(row.line, "CLOAD expects node/set, DOF and magnitude"));
            }
            let dof = fields
                .get(1)
                .and_then(|value| value.parse::<u8>().ok())
                .filter(|dof| (1..=6).contains(dof))
                .ok_or_else(|| error(row.line, "CLOAD DOF must be 1..6"))?;
            let magnitude = number(fields.get(2).map(String::as_str).unwrap_or_default(), row.line)?;
            let target_text = fields.first().map(String::as_str).unwrap_or_default();
            let target = target_name(target_text, row.line)?;
            let target_nodes = members(model, input_sets, target_text, row.line)?;
            let key = (Family::Force, target, dof);
            if !seen.insert(key.clone()) {
                return Err(unsupported(row.line, "duplicate CLOAD definition for the same target and DOF in one step"));
            }
            let driven = Driven { value: magnitude, line: row.line, boundary: None, members: target_nodes.len(), initial };
            insert(state, key, driven, &target_nodes, limits, row.line)?;
            if target_nodes.iter().any(|node| transformed_nodes.contains(node)) {
                return Err(unsupported(row.line, "history target is covered by a *TRANSFORM local coordinate system"));
            }
        }
    }
    Ok(())
}

fn family_mode(body: &[Block<'_>], keyword: &str) -> Result<Option<bool>> {
    let mut found = false;
    let mut has_new = false;
    let mut all_new = true;
    for block in body.iter().filter(|block| block.keyword == keyword) {
        let op = canonical(value(block, "OP", "MOD")?);
        if !matches!(op.as_str(), "MOD" | "NEW") {
            return Err(error(block.line, "OP must be NEW or MOD"));
        }
        found = true;
        has_new |= op == "NEW";
        all_new &= op == "NEW";
    }
    if has_new && !all_new {
        return Err(unsupported(
            body.iter().find(|block| block.keyword == keyword).map(|block| block.line).unwrap_or_default(),
            format!("*{keyword} mixes OP=NEW and OP=MOD in one step"),
        ));
    }
    Ok(found.then_some(has_new))
}

fn prepare(
    block: &Block<'_>,
    body: &[Block<'_>],
    state: &HistoryState,
    model: &Model,
    input_sets: &BTreeSet<String>,
    limits: Limits,
    target_ids: &mut TargetIds,
    transformed_nodes: &BTreeSet<u32>,
) -> Result<(Procedure, HistoryState)> {
    allow(block, &["NAME", "NLGEOM"])?;
    if canonical(value(block, "NLGEOM", "NO")?) != "NO" {
        return Err(unsupported(block.line, "nonlinear geometry"));
    }
    if !block.rows.is_empty() {
        return Err(unsupported(block.line, "step description rows are not mapped"));
    }
    let boundary_mode = family_mode(body, "BOUNDARY")?;
    let force_mode = family_mode(body, "CLOAD")?;
    let mut procedure = None;
    let mut has_force_block = false;
    for part in body {
        match part.keyword.as_str() {
            "STATIC" => {
                allow(part, &[])?;
                if procedure.is_some() {
                    return Err(error(part.line, "multiple analysis procedures in one STEP"));
                }
                if part.rows.len() > 1 {
                    return Err(error(part.line, "STATIC has at most one data row"));
                }
                if let Some(row) = part.rows.first() {
                    let values = row_fields(row.text, row.line)?;
                    if values.len() > 4 {
                        return Err(error(row.line, "STATIC expects at most four increment/time entries"));
                    }
                    for value in values.iter().filter(|value| !value.is_empty()) {
                        if number(value, row.line)? <= 0.0 {
                            return Err(error(row.line, "STATIC increments/time must be positive"));
                        }
                    }
                }
                procedure = Some(Procedure::Static {});
            }
            "FREQUENCY" | "BUCKLE" => {
                allow(part, &["EIGENSOLVER"])?;
                if canonical(value(part, "EIGENSOLVER", "LANCZOS")?) != "LANCZOS" {
                    return Err(unsupported(part.line, "non-Lanczos eigen extraction"));
                }
                if procedure.is_some() {
                    return Err(error(part.line, "multiple analysis procedures"));
                }
                let row = part.rows.first().ok_or_else(|| error(part.line, "eigen extraction needs a mode count"))?;
                if part.rows.len() != 1 {
                    return Err(error(row.line, "eigen extraction needs one data row"));
                }
                let values = row_fields(row.text, row.line)?;
                let modes = label(values.first().map(String::as_str).unwrap_or_default(), row.line)?;
                if modes > fieldcraft_model::steps::MAX_MODES {
                    return Err(error(row.line, "mode count exceeds limit"));
                }
                let bound = |index: usize| values.get(index).filter(|value| !value.is_empty()).map(|value| number(value, row.line)).transpose();
                procedure = Some(if part.keyword == "FREQUENCY" {
                    if values.len() > 3 {
                        return Err(unsupported(row.line, "additional frequency controls"));
                    }
                    Procedure::Frequency { modes, min_hz: bound(1)?, max_hz: bound(2)? }
                } else {
                    if values.len() > 1 {
                        return Err(unsupported(row.line, "buckling shift and extraction controls"));
                    }
                    Procedure::Buckle { modes, shift: None }
                });
            }
            "BOUNDARY" => {}
            "CLOAD" => has_force_block = true,
            "OUTPUT" | "NODEOUTPUT" | "ELEMENTOUTPUT" | "NODEPRINT" | "ELPRINT" | "ENERGYOUTPUT" | "PREPRINT" => {}
            _ => return Err(unsupported(part.line, format!("*{} in STEP", part.keyword))),
        }
    }
    let procedure = procedure.ok_or_else(|| unsupported(block.line, "STEP has no supported procedure"))?;
    if matches!(procedure, Procedure::Frequency { .. }) && has_force_block {
        return Err(unsupported(block.line, "CLOAD has no effect in a frequency step"));
    }
    let frequency = matches!(procedure, Procedure::Frequency { .. });
    let mut next = state.clone();
    let mut seen = BTreeSet::new();
    let mut boundary_cleared = false;
    let mut force_cleared = false;
    for part in body {
        match part.keyword.as_str() {
            "BOUNDARY" => {
                let clear = boundary_mode == Some(true) && !boundary_cleared;
                apply(part, &mut next, model, input_sets, limits, &mut seen, target_ids, clear, None, frequency, transformed_nodes)?;
                boundary_cleared |= clear;
            }
            "CLOAD" => {
                let clear = force_mode == Some(true) && !force_cleared;
                apply(part, &mut next, model, input_sets, limits, &mut seen, target_ids, clear, None, false, transformed_nodes)?;
                force_cleared |= clear;
            }
            _ => {}
        }
    }
    Ok((procedure, next))
}

fn reserve_model_names(model: &Model) -> BTreeSet<String> {
    model.sets().keys().chain(model.boundary_conditions().keys()).chain(model.loads().keys()).map(|name| name.to_ascii_uppercase()).collect()
}

fn allocate_name(reserved: &mut BTreeSet<String>, sequence: &mut usize, line: usize) -> Result<String> {
    for _ in 0..=fieldcraft_model::MAX_DEFINITIONS {
        *sequence = sequence.saturating_add(1);
        let name = format!("InpAnalysis{sequence}");
        if reserved.insert(name.to_ascii_uppercase()) {
            return Ok(name);
        }
    }
    Err(error(line, "analysis name allocation limit reached"))
}

#[derive(Clone)]
struct BoundaryNode {
    kind: BoundaryKind,
    values: [Option<f64>; 6],
    targets: BTreeSet<String>,
    line: usize,
}

#[derive(Clone)]
struct BoundaryGroup {
    kind: BoundaryKind,
    values: [Option<f64>; 6],
    nodes: BTreeSet<u32>,
    targets: BTreeSet<String>,
    line: usize,
}

#[derive(Clone)]
struct LoadNode {
    vector: [f64; 3],
    targets: BTreeSet<String>,
    line: usize,
}

#[derive(Clone)]
struct LoadGroup {
    moment: bool,
    vector: [f64; 3],
    nodes: BTreeSet<u32>,
    targets: BTreeSet<String>,
    line: usize,
}

fn boundary_key(kind: BoundaryKind) -> u8 {
    match kind {
        BoundaryKind::Displacement => 0,
        BoundaryKind::Pinned => 1,
        BoundaryKind::Fixed => 2,
        BoundaryKind::Symmetry(Symmetry::X) => 3,
        BoundaryKind::Symmetry(Symmetry::Y) => 4,
        BoundaryKind::Symmetry(Symmetry::Z) => 5,
        BoundaryKind::Symmetry(Symmetry::Normal(_)) => 6,
    }
}

fn optional_bits(value: Option<f64>) -> u64 {
    value.map(|value| value.to_bits().wrapping_add(1)).unwrap_or(0)
}

fn target_key(target: &Target) -> String {
    match target {
        Target::Set(name) => format!("set:{name}"),
        Target::ReferencePoint(name) => format!("reference:{name}"),
    }
}

fn boundary_identity(boundary: &BoundaryCondition) -> String {
    match boundary {
        BoundaryCondition::Fixed { target } => format!("fixed:{}", target_key(target)),
        BoundaryCondition::Symmetry { target, plane } => {
            let plane = match plane {
                Symmetry::X => "x".to_string(),
                Symmetry::Y => "y".to_string(),
                Symmetry::Z => "z".to_string(),
                Symmetry::Normal([x, y, z]) => format!("normal:{:016x}:{:016x}:{:016x}", x.to_bits(), y.to_bits(), z.to_bits()),
            };
            format!("symmetry:{plane}:{}", target_key(target))
        }
        BoundaryCondition::Displacement { target, values } => {
            let values = values.iter().map(|value| format!("{:016x}", optional_bits(*value))).collect::<Vec<_>>().join(":");
            format!("displacement:{values}:{}", target_key(target))
        }
        BoundaryCondition::Temperature { target, value } => format!("temperature:{:016x}:{}", value.to_bits(), target_key(target)),
    }
}

fn load_identity(load: &Load) -> Option<String> {
    match load {
        Load::Force { target, vector: [x, y, z], distribution } => {
            Some(format!("force:{}:{:016x}:{:016x}:{:016x}:{distribution:?}", target_key(target), x.to_bits(), y.to_bits(), z.to_bits()))
        }
        Load::Moment { target, vector: [x, y, z], distribution } => {
            Some(format!("moment:{}:{:016x}:{:016x}:{:016x}:{distribution:?}", target_key(target), x.to_bits(), y.to_bits(), z.to_bits()))
        }
        _ => None,
    }
}

#[derive(Clone, Copy)]
enum LoadSource {
    Set,
    Node,
}

fn load_identity_for(load: &Load, source: LoadSource) -> Option<String> {
    // The model target may be the same set after materialization, but the
    // input target spelling is part of CLOAD identity: a set row and a node
    // row add instead of collapsing into one definition.
    let source = match source {
        LoadSource::Set => "set",
        LoadSource::Node => "node",
    };
    Some(format!("{source}:{}", load_identity(load)?))
}

struct SetIndex {
    by_name: BTreeMap<String, (String, Vec<u32>)>,
    by_nodes: BTreeMap<Vec<u32>, String>,
}

impl SetIndex {
    fn from_model(model: &Model) -> Self {
        let mut index = Self { by_name: BTreeMap::new(), by_nodes: BTreeMap::new() };
        for name in model.sets().keys() {
            let nodes = model.resolve_set_partial(name).map(|set| set.nodes).unwrap_or_default();
            index.insert(name, &nodes);
        }
        index
    }

    fn insert(&mut self, name: &str, nodes: &[u32]) {
        self.by_name.insert(name.to_ascii_uppercase(), (name.to_string(), nodes.to_vec()));
        self.by_nodes.entry(nodes.to_vec()).or_insert_with(|| name.to_string());
    }

    fn named(&self, identity: &str, nodes: &[u32]) -> Option<String> {
        self.by_name.get(identity).filter(|(_, current)| current == nodes).map(|(name, _)| name.clone())
    }

    fn alias_named(&mut self, identity: &str, name: &str, nodes: &[u32]) {
        self.by_name.insert(identity.to_string(), (name.to_string(), nodes.to_vec()));
    }
}

struct DefinitionIndex {
    boundaries: BTreeMap<String, String>,
    loads: BTreeMap<String, String>,
}

impl DefinitionIndex {
    fn from_model(model: &Model) -> Self {
        let boundaries = model.boundary_conditions().iter().map(|(name, boundary)| (boundary_identity(boundary), name.clone())).collect();
        let loads = model.loads().iter().filter_map(|(name, load)| load_identity_for(load, LoadSource::Set).map(|key| (key, name.clone()))).collect();
        Self { boundaries, loads }
    }
}

struct Materialized {
    sets: Vec<(String, Set)>,
    boundaries: Vec<(String, BoundaryCondition)>,
    loads: Vec<(String, Load)>,
    step: Step,
    sequence: usize,
}

fn restore_target(
    index: &mut SetIndex,
    identities: &BTreeSet<String>,
    nodes: &[u32],
    line: usize,
    sequence: &mut usize,
    reserved: &mut BTreeSet<String>,
    staged: &mut Vec<(String, Set)>,
) -> Result<String> {
    let named_target = identities.iter().any(|identity| identity.parse::<u32>().is_err());
    for identity in identities {
        if identity.parse::<u32>().is_ok() {
            continue;
        }
        if let Some(name) = index.named(identity, nodes) {
            return Ok(name);
        }
    }
    if !named_target && let Some(name) = index.by_nodes.get(nodes) {
        return Ok(name.clone());
    }
    let name = allocate_name(reserved, sequence, line)?;
    let set = Set { nodes: nodes.to_vec(), ..Set::default() };
    index.insert(&name, nodes);
    for identity in identities.iter().filter(|identity| identity.parse::<u32>().is_err()) {
        index.alias_named(identity, &name, nodes);
    }
    staged.push((name.clone(), set));
    Ok(name)
}

fn cached_members(
    model: &Model,
    input_sets: &BTreeSet<String>,
    cache: &mut BTreeMap<String, Vec<u32>>,
    target: &str,
    line: usize,
) -> Result<Vec<u32>> {
    if let Some(nodes) = cache.get(target) {
        return Ok(nodes.clone());
    }
    let nodes = members(model, input_sets, target, line)?;
    cache.insert(target.to_string(), nodes.clone());
    Ok(nodes)
}

fn materialize(
    model: &Model,
    input_sets: &BTreeSet<String>,
    mut step: Step,
    state: &HistoryState,
    units: &DisplaySystem,
    mut sequence: usize,
    frequency: bool,
    set_index: &mut SetIndex,
    definitions: &mut DefinitionIndex,
    reserved: &mut BTreeSet<String>,
    target_cache: &mut BTreeMap<String, Vec<u32>>,
) -> Result<Materialized> {
    let mut boundary_nodes: BTreeMap<(u8, u32), BoundaryNode> = BTreeMap::new();
    let mut boundary_targets: BTreeMap<(String, u8), BoundaryGroup> = BTreeMap::new();
    let mut load_nodes: BTreeMap<(bool, u32), LoadNode> = BTreeMap::new();
    let mut load_targets: BTreeMap<(String, bool), LoadGroup> = BTreeMap::new();
    let mut work = 0usize;
    for ((family, target, dof), driven) in &state.values {
        check_work(work)?;
        work = work.saturating_add(1);
        let nodes = cached_members(model, input_sets, target_cache, target, driven.line)?;
        match family {
            Family::Boundary => {
                let Some(kind) = driven.boundary else { continue };
                for node in nodes {
                    if target.parse::<u32>().is_ok() {
                        let key = (boundary_key(kind), node);
                        let entry = boundary_nodes.entry(key).or_insert_with(|| BoundaryNode {
                            kind,
                            values: [None; 6],
                            targets: BTreeSet::new(),
                            line: driven.line,
                        });
                        entry.targets.insert(target.clone());
                        if matches!(kind, BoundaryKind::Displacement) {
                            let kind_quantity = if *dof <= 3 { Quantity::Length } else { Quantity::Angle };
                            let value = quantity(units, kind_quantity, driven.value, driven.line)?;
                            let slot = entry
                                .values
                                .get_mut(usize::from(*dof - 1))
                                .ok_or_else(|| error(driven.line, "boundary DOF is outside the stored values"))?;
                            *slot = Some(value);
                        }
                    } else {
                        let key = (target.clone(), boundary_key(kind));
                        let entry = boundary_targets.entry(key).or_insert_with(|| BoundaryGroup {
                            kind,
                            values: [None; 6],
                            nodes: BTreeSet::new(),
                            targets: BTreeSet::new(),
                            line: driven.line,
                        });
                        entry.nodes.insert(node);
                        entry.targets.insert(target.clone());
                        if matches!(kind, BoundaryKind::Displacement) {
                            let kind_quantity = if *dof <= 3 { Quantity::Length } else { Quantity::Angle };
                            let value = quantity(units, kind_quantity, driven.value, driven.line)?;
                            let slot = entry
                                .values
                                .get_mut(usize::from(*dof - 1))
                                .ok_or_else(|| error(driven.line, "boundary DOF is outside the stored values"))?;
                            *slot = Some(value);
                        }
                    }
                }
            }
            Family::Force if !frequency => {
                let moment = *dof > 3;
                let kind_quantity = if moment { Quantity::Moment } else { Quantity::Force };
                let component = usize::from((*dof - 1) % 3);
                let value = quantity(units, kind_quantity, driven.value, driven.line)?;
                if target.parse::<u32>().is_ok() {
                    for node in nodes {
                        let entry = load_nodes.entry((moment, node)).or_insert_with(|| LoadNode {
                            vector: [0.0; 3],
                            targets: BTreeSet::new(),
                            line: driven.line,
                        });
                        entry.targets.insert(target.clone());
                        let slot = entry.vector.get_mut(component).ok_or_else(|| error(driven.line, "load DOF is outside the stored vector"))?;
                        *slot += value;
                    }
                } else {
                    let entry = load_targets.entry((target.clone(), moment)).or_insert_with(|| LoadGroup {
                        moment,
                        vector: [0.0; 3],
                        nodes: BTreeSet::new(),
                        targets: BTreeSet::new(),
                        line: driven.line,
                    });
                    entry.nodes.extend(nodes);
                    entry.targets.insert(target.clone());
                    let slot = entry.vector.get_mut(component).ok_or_else(|| error(driven.line, "load DOF is outside the stored vector"))?;
                    *slot += value;
                }
            }
            Family::Force => {}
        }
    }

    let mut boundary_groups: BTreeMap<(u8, [u64; 6]), BoundaryGroup> = BTreeMap::new();
    for ((_, node), entry) in boundary_nodes {
        let key = (boundary_key(entry.kind), entry.values.map(optional_bits));
        let group = boundary_groups.entry(key).or_insert_with(|| BoundaryGroup {
            kind: entry.kind,
            values: entry.values,
            nodes: BTreeSet::new(),
            targets: BTreeSet::new(),
            line: entry.line,
        });
        group.nodes.insert(node);
        group.targets.extend(entry.targets);
    }
    let mut load_groups: BTreeMap<(bool, [u64; 3]), LoadGroup> = BTreeMap::new();
    for ((moment, node), entry) in load_nodes {
        let key = (moment, entry.vector.map(f64::to_bits));
        let group = load_groups.entry(key).or_insert_with(|| LoadGroup {
            moment,
            vector: entry.vector,
            nodes: BTreeSet::new(),
            targets: BTreeSet::new(),
            line: entry.line,
        });
        group.nodes.insert(node);
        group.targets.extend(entry.targets);
    }

    let mut sets = Vec::new();
    let mut boundaries = Vec::new();
    let mut loads = Vec::new();
    let mut boundary_groups: Vec<BoundaryGroup> = boundary_groups.into_values().collect();
    boundary_groups.extend(boundary_targets.into_values());
    for group in boundary_groups {
        let nodes: Vec<u32> = group.nodes.into_iter().collect();
        let set_name = restore_target(set_index, &group.targets, &nodes, group.line, &mut sequence, reserved, &mut sets)?;
        let target = Target::Set(set_name);
        let boundary = match group.kind {
            BoundaryKind::Fixed => BoundaryCondition::Fixed { target },
            BoundaryKind::Symmetry(plane) => BoundaryCondition::Symmetry { target, plane },
            BoundaryKind::Displacement => BoundaryCondition::Displacement { target, values: group.values },
            BoundaryKind::Pinned => BoundaryCondition::Displacement { target, values: [Some(0.0), Some(0.0), Some(0.0), None, None, None] },
        };
        let key = boundary_identity(&boundary);
        let name = if let Some(name) = definitions.boundaries.get(&key) {
            name.clone()
        } else {
            let name = allocate_name(reserved, &mut sequence, group.line)?;
            definitions.boundaries.insert(key, name.clone());
            boundaries.push((name.clone(), boundary));
            name
        };
        step.boundary_conditions.insert(name);
    }
    let load_groups: Vec<(LoadSource, LoadGroup)> = load_targets
        .into_values()
        .map(|group| (LoadSource::Set, group))
        .chain(load_groups.into_values().map(|group| (LoadSource::Node, group)))
        .collect();
    for (source, group) in load_groups {
        if group.vector.iter().all(|value| *value == 0.0) {
            continue;
        }
        let nodes: Vec<u32> = group.nodes.into_iter().collect();
        let set_name = restore_target(set_index, &group.targets, &nodes, group.line, &mut sequence, reserved, &mut sets)?;
        let target = Target::Set(set_name);
        let load = if group.moment {
            Load::Moment { target, vector: group.vector, distribution: Distribution::PerNode }
        } else {
            Load::Force { target, vector: group.vector, distribution: Distribution::PerNode }
        };
        let Some(key) = load_identity_for(&load, source) else {
            return Err(error(group.line, "internal load mapping did not produce a force or moment"));
        };
        let name = if let Some(name) = definitions.loads.get(&key) {
            name.clone()
        } else {
            let name = allocate_name(reserved, &mut sequence, group.line)?;
            definitions.loads.insert(key, name.clone());
            loads.push((name.clone(), load));
            name
        };
        step.loads.insert(name);
    }
    Ok(Materialized { sets, boundaries, loads, step, sequence })
}

fn rollback(model: &mut Model, loads: &[String], boundaries: &[String], sets: &[String]) {
    for name in loads.iter().rev() {
        let _ = model.remove_load(name);
    }
    for name in boundaries.iter().rev() {
        let _ = model.remove_boundary_condition(name);
    }
    for name in sets.iter().rev() {
        let _ = model.remove_set(name);
    }
}

fn commit(model: &mut Model, plan: Materialized, line: usize) -> Result<()> {
    let Materialized { sets, boundaries, loads, step, .. } = plan;
    let mut added_sets = Vec::new();
    let mut added_boundaries = Vec::new();
    let mut added_loads = Vec::new();
    for (name, set) in sets {
        if let Err(add_error) = model.add_set(&name, set) {
            rollback(model, &added_loads, &added_boundaries, &added_sets);
            return Err(error(line, add_error.to_string()));
        }
        added_sets.push(name);
    }
    for (name, boundary) in boundaries {
        if let Err(add_error) = model.add_boundary_condition(&name, boundary) {
            rollback(model, &added_loads, &added_boundaries, &added_sets);
            return Err(error(line, add_error.to_string()));
        }
        added_boundaries.push(name);
    }
    for (name, load) in loads {
        if let Err(add_error) = model.add_load(&name, load) {
            rollback(model, &added_loads, &added_boundaries, &added_sets);
            return Err(error(line, add_error.to_string()));
        }
        added_loads.push(name);
    }
    if let Err(add_error) = model.add_step(step, None) {
        rollback(model, &added_loads, &added_boundaries, &added_sets);
        return Err(error(line, add_error.to_string()));
    }
    Ok(())
}

fn cancelled(error: &Error) -> bool {
    matches!(error, Error::Model(fieldcraft_model::Error::Cancelled))
}

fn with_line(error: Error, line: usize) -> Error {
    match error {
        Error::Unsupported { .. } | Error::Parse { .. } => error,
        other => unsupported(line, other.to_string()),
    }
}

fn report_unread(read: &mut Read, error: Error, line: usize, suffix: &str) -> Result<()> {
    if cancelled(&error) {
        return Err(error);
    }
    read.warnings.push(format!("{}; {suffix}", with_line(error, line)));
    Ok(())
}

fn rebuild_unknown(read: &mut Read, consumed: &BTreeSet<usize>) {
    let remaining: Vec<_> = read.unknown_all.iter().filter(|unknown| !consumed.contains(&unknown.line)).cloned().collect();
    read.unknown_count = remaining.len();
    read.unknown = remaining.iter().take(super::MAX_UNKNOWN_SHOWN).cloned().collect();
    read.unknown_all = remaining;
}

/// Import the supported linear history blocks and remove only those blocks from the unread report.
pub(super) fn import(text: &str, read: &mut Read, preset: Preset, limits: Limits) -> Result<()> {
    let blocks = lex::blocks(text, super::MAX_RECORDS)?;
    let input_sets: BTreeSet<String> = read.model.sets().keys().map(|name| name.to_ascii_uppercase()).collect();
    let unresolved_transform_line = read.transforms.iter().find(|transform| transform.target.is_none()).map(|transform| transform.line);
    if let Some(line) = unresolved_transform_line {
        report_unread(
            read,
            unsupported(line, "*TRANSFORM NSET cannot be resolved; local coordinate coverage is unknown"),
            line,
            "all history remains unread",
        )?;
        return Ok(());
    }
    let transform_targets: BTreeMap<String, usize> =
        read.transforms.iter().filter_map(|transform| transform.target.clone().map(|target| (target, transform.line))).collect();
    let mut transformed_nodes = BTreeSet::new();
    for (target, line) in transform_targets {
        if !input_sets.contains(&target.to_ascii_uppercase()) {
            report_unread(
                read,
                unsupported(line, "*TRANSFORM NSET cannot be resolved; local coordinate coverage is unknown"),
                line,
                "all history remains unread",
            )?;
            return Ok(());
        }
        let Some(name) = read.model.sets().keys().find(|name| name.eq_ignore_ascii_case(&target)) else {
            report_unread(
                read,
                unsupported(line, "*TRANSFORM NSET cannot be resolved; local coordinate coverage is unknown"),
                line,
                "all history remains unread",
            )?;
            return Ok(());
        };
        let Some(set) = read.model.set(name) else {
            report_unread(
                read,
                unsupported(line, "*TRANSFORM NSET cannot be resolved; local coordinate coverage is unknown"),
                line,
                "all history remains unread",
            )?;
            return Ok(());
        };
        if set.nodes.is_empty() {
            report_unread(
                read,
                unsupported(line, "*TRANSFORM NSET cannot be resolved; local coordinate coverage is unknown"),
                line,
                "all history remains unread",
            )?;
            return Ok(());
        }
        transformed_nodes.extend(&set.nodes);
    }
    let mut consumed = BTreeSet::new();
    let mut state = HistoryState::default();
    let mut target_ids = TargetIds::default();
    let mut sequence = 0usize;
    let mut history_members = 0usize;
    let mut set_index = SetIndex::from_model(&read.model);
    let mut definitions = DefinitionIndex::from_model(&read.model);
    let mut reserved = reserve_model_names(&read.model);
    let mut target_cache = BTreeMap::new();
    let units = DisplaySystem::new(preset);
    let mut initial_seen = BTreeSet::new();
    let mut initial_new_cleared = false;
    let mut seen_step = false;
    let reserved_step_names: BTreeSet<String> = blocks
        .iter()
        .filter(|candidate| candidate.keyword == "STEP")
        .filter_map(|candidate| candidate.params.get("NAME").and_then(|value| value.as_deref()))
        .filter_map(|raw| lex::unquote(raw, 0).ok())
        .filter(|name| !name.is_empty())
        .map(|name| name.to_ascii_uppercase())
        .collect();
    let mut imported_step_names: BTreeSet<String> = read.model.steps().iter().map(|step| step.name.to_ascii_uppercase()).collect();
    let mut generated_step_sequence = read.model.steps().len() + 1;
    let mut seen_general_step = false;
    let mut index = 0usize;
    while index < blocks.len() {
        let Some(block) = blocks.get(index) else { break };
        if block.keyword == "BOUNDARY" {
            if seen_step {
                report_unread(
                    read,
                    unsupported(block.line, "model-level *BOUNDARY appears after a *STEP"),
                    block.line,
                    "this and subsequent steps remain unread",
                )?;
                break;
            }
            let op = match value(block, "OP", "MOD") {
                Ok(value) => canonical(value),
                Err(error) => {
                    report_unread(read, error, block.line, "this and subsequent steps remain unread")?;
                    return Ok(());
                }
            };
            let clear = op == "NEW" && !initial_new_cleared;
            if let Err(error) = apply(
                block,
                &mut state,
                &read.model,
                &input_sets,
                limits,
                &mut initial_seen,
                &mut target_ids,
                clear,
                Some(block.line),
                false,
                &transformed_nodes,
            ) {
                report_unread(read, error, block.line, "all history remains unread")?;
                return Ok(());
            }
            initial_new_cleared |= clear;
            if state
                .values
                .iter()
                .any(|((family, _, _), driven)| *family == Family::Boundary && driven.initial == Some(block.line) && driven.value != 0.0)
            {
                report_unread(
                    read,
                    unsupported(block.line, "initial boundaries must be zero"),
                    block.line,
                    "this and subsequent steps remain unread",
                )?;
                return Ok(());
            }
            index += 1;
            continue;
        }
        if block.keyword != "STEP" {
            index += 1;
            continue;
        }
        let Some(end) = (index.saturating_add(1)..blocks.len()).find(|at| blocks.get(*at).is_some_and(|candidate| candidate.keyword == "ENDSTEP"))
        else {
            report_unread(read, error(block.line, "STEP has no END STEP"), block.line, "this and subsequent steps remain unread")?;
            break;
        };
        let Some(body) = blocks.get(index.saturating_add(1)..end) else {
            report_unread(read, error(block.line, "STEP scope is invalid"), block.line, "this and subsequent steps remain unread")?;
            break;
        };
        let raw_name = match value(block, "NAME", "") {
            Ok(value) => value,
            Err(error) => {
                report_unread(read, error, block.line, "this and subsequent steps remain unread")?;
                break;
            }
        };
        let explicit_name = match lex::unquote(raw_name, block.line) {
            Ok(name) => name,
            Err(error) => {
                report_unread(read, error, block.line, "this and subsequent steps remain unread")?;
                break;
            }
        };
        let step_name = if explicit_name.is_empty() {
            loop {
                let candidate = format!("ImportedStep{generated_step_sequence}");
                generated_step_sequence = generated_step_sequence.saturating_add(1);
                let key = candidate.to_ascii_uppercase();
                if !reserved_step_names.contains(&key) && !imported_step_names.contains(&key) {
                    break candidate;
                }
            }
        } else {
            explicit_name
        };
        let step_key = step_name.to_ascii_uppercase();
        if imported_step_names.contains(&step_key) {
            report_unread(
                read,
                unsupported(block.line, "step name already exists (names differing only by case are the same)"),
                block.line,
                "this and subsequent steps remain unread",
            )?;
            break;
        }
        let (procedure, next) = match prepare(block, body, &state, &read.model, &input_sets, limits, &mut target_ids, &transformed_nodes) {
            Ok(value) => value,
            Err(error) => {
                report_unread(read, error, block.line, "this and subsequent steps remain unread")?;
                break;
            }
        };
        if matches!(procedure, Procedure::Buckle { .. }) && seen_general_step {
            report_unread(
                read,
                unsupported(block.line, "buckling after a general analysis step is not mapped"),
                block.line,
                "this and subsequent steps remain unread",
            )?;
            break;
        }
        let mut step = Step::new(step_name, procedure.clone());
        let visible = if matches!(procedure, Procedure::Frequency { .. }) { next.boundary_only() } else { next.clone() };
        history_members = history_members.saturating_add(visible.members);
        if history_members > limits.members {
            report_unread(read, error(block.line, "history member limit reached"), block.line, "this and subsequent steps remain unread")?;
            break;
        }
        let attempt = (|| {
            if let Procedure::Frequency { min_hz, max_hz, .. } = &mut step.procedure {
                for bound in [min_hz.as_mut(), max_hz.as_mut()].into_iter().flatten() {
                    let current = *bound;
                    *bound = quantity(&units, Quantity::Frequency, current, block.line)?;
                }
            }
            let plan = materialize(
                &read.model,
                &input_sets,
                step,
                &visible,
                &units,
                sequence,
                matches!(procedure, Procedure::Frequency { .. }),
                &mut set_index,
                &mut definitions,
                &mut reserved,
                &mut target_cache,
            )?;
            let next_sequence = plan.sequence;
            commit(&mut read.model, plan, block.line)?;
            Ok::<usize, Error>(next_sequence)
        })();
        let next_sequence = match attempt {
            Ok(next_sequence) => next_sequence,
            Err(error) => {
                report_unread(read, error, block.line, "this and subsequent steps remain unread")?;
                break;
            }
        };
        sequence = next_sequence;
        imported_step_names.insert(step_key);
        let initial_used: BTreeSet<usize> = visible.values.values().filter_map(|driven| driven.initial).collect();
        consumed.extend(initial_used);
        if matches!(procedure, Procedure::Static {}) {
            state = next;
            seen_general_step = true;
        }
        seen_step = true;
        consumed.insert(block.line);
        if let Some(end_block) = blocks.get(end) {
            consumed.insert(end_block.line);
        }
        consumed.extend(
            body.iter()
                .filter(|part| matches!(part.keyword.as_str(), "STATIC" | "FREQUENCY" | "BUCKLE" | "BOUNDARY" | "CLOAD"))
                .map(|part| part.line),
        );
        index = end.saturating_add(1);
    }
    rebuild_unknown(read, &consumed);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::{Limits, Read, read_analysis_limited, read_in};
    use super::*;
    use fieldcraft_model::{BeamProfile, Elastic, ElementKind, Material, Model, Section, SectionKind};
    use std::collections::{BTreeMap, BTreeSet};
    use std::panic::{AssertUnwindSafe, catch_unwind};

    const MODEL: &str = "*NODE\n1,0,0,0\n2,1,0,0\n*ELEMENT,TYPE=T3D2,ELSET=Bar\n1,1,2\n*NSET,NSET=Bar\n1,2\n*NSET,NSET=Root\n1\n*NSET,NSET=Tip\n2\n*MATERIAL,NAME=Steel\n*ELASTIC\n210000,0.3\n*DENSITY\n1\n*SOLID SECTION,ELSET=Bar,MATERIAL=Steel\n1\n";
    const MULTI_MODEL: &str = "*NODE\n1,0,0,0\n2,1,0,0\n3,2,0,0\n4,3,0,0\n*ELEMENT,TYPE=T3D2,ELSET=Bar\n1,1,2\n2,2,3\n3,3,4\n*NSET,NSET=Bar\n1,2,3,4\n*NSET,NSET=Root\n1\n*NSET,NSET=Tip\n3,4\n*NSET,NSET=Mid\n2,3\n*MATERIAL,NAME=Steel\n*ELASTIC\n210000,0.3\n*DENSITY\n1\n*SOLID SECTION,ELSET=Bar,MATERIAL=Steel\n1\n";
    const ROUNDTRIP_MODEL: &str = "*NODE\n1,0,0,0\n2,1,0,0\n3,2,0,0\n4,3,0,0\n*ELEMENT,TYPE=T3D2,ELSET=Bar\n1,1,2\n2,2,3\n3,3,4\n*MATERIAL,NAME=Steel\n*ELASTIC\n210000,0.3\n*DENSITY\n1\n*SOLID SECTION,ELSET=Bar,MATERIAL=Steel\n1\n";
    const INVARIANT_MODEL: &str = "*NODE\n1,0,0,0\n2,1,0,0\n3,2,0,0\n4,3,0,0\n*ELEMENT,TYPE=T3D2,ELSET=Bar\n1,1,2\n2,2,3\n3,3,4\n*NSET,NSET=Root\n1\n*MATERIAL,NAME=Steel\n*ELASTIC\n210000,0.3\n*DENSITY\n1\n*SOLID SECTION,ELSET=Bar,MATERIAL=Steel\n1\n";
    const INVARIANT_BEAM_MODEL: &str = "*NODE\n1,0,0,0\n2,1,0,0\n3,2,0,0\n4,3,0,0\n*ELEMENT,TYPE=B31,ELSET=Bar\n1,1,2\n2,2,3\n3,3,4\n*NSET,NSET=Root\n1\n*MATERIAL,NAME=Steel\n*ELASTIC\n210000,0.3\n*DENSITY\n1\n*BEAM SECTION,ELSET=Bar,MATERIAL=Steel,SECTION=CIRC\n0.1\n0,1,0\n";

    fn assert_model_data_unchanged(deck: &str) {
        let base = read_in(deck, Preset::MKgS).unwrap();
        let got = read_analysis_limited(deck, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.nodes(), base.model.nodes());
        assert_eq!(got.model.elements(), base.model.elements());
        assert_eq!(got.model.materials(), base.model.materials());
        assert_eq!(got.model.sections(), base.model.sections());
        assert_eq!(got.transforms, base.transforms);
        for (name, set) in base.model.sets() {
            assert_eq!(got.model.set(name), Some(set), "source set {name}");
        }
    }

    fn roundtrip_base() -> Model {
        let mut model = read_in(ROUNDTRIP_MODEL, Preset::MKgS).unwrap().model;
        for (name, nodes) in
            [("Left", vec![1]), ("Right", vec![4]), ("SymX", vec![1, 2]), ("SymY", vec![2, 3]), ("Tip", vec![3, 4]), ("Mid", vec![2])]
        {
            model.add_set(name, Set { nodes, ..Set::default() }).unwrap();
        }
        model
    }

    fn roundtrip_set(name: &str) -> Target {
        Target::Set(name.into())
    }

    fn roundtrip_force(target: &str, x: f64) -> Load {
        Load::Force { target: roundtrip_set(target), vector: [x, 0.0, 0.0], distribution: Distribution::PerNode }
    }

    fn roundtrip_step(model: &mut Model, name: &str, boundaries: &[&str], loads: &[&str]) {
        roundtrip_step_with(model, name, Procedure::Static {}, boundaries, loads);
    }

    fn roundtrip_step_with(model: &mut Model, name: &str, procedure: Procedure, boundaries: &[&str], loads: &[&str]) {
        let mut step = Step::new(name, procedure);
        for boundary in boundaries {
            step.boundary_conditions.insert((*boundary).into());
        }
        for load in loads {
            step.loads.insert((*load).into());
        }
        model.add_step(step, None).unwrap();
    }

    fn roundtrip_moment_model() -> Model {
        let mut model = Model::default();
        model.set_node(1, [0.0, 0.0, 0.0]).unwrap();
        model.set_node(2, [1.0, 0.0, 0.0]).unwrap();
        model.set_element(1, ElementKind::Beam2, vec![1, 2]).unwrap();
        model.add_set("Bar", Set { elements: vec![1], ..Set::default() }).unwrap();
        model.add_set("Root", Set { nodes: vec![1], ..Set::default() }).unwrap();
        model.add_set("Tip", Set { nodes: vec![2], ..Set::default() }).unwrap();
        model
            .add_material(
                "Steel",
                Material { elastic: Some(Elastic::Isotropic { young: 210000.0, poisson: 0.3 }), density: Some(1.0), ..Material::default() },
            )
            .unwrap();
        model
            .add_section(
                "Section",
                Section {
                    set: "Bar".into(),
                    material: "Steel".into(),
                    kind: SectionKind::Beam { profile: BeamProfile::Circle { radius: 0.5 }, orientation: [0.0, 1.0, 0.0] },
                    formulation: Default::default(),
                },
            )
            .unwrap();
        model.add_boundary_condition("Fix", BoundaryCondition::Fixed { target: roundtrip_set("Root") }).unwrap();
        model
            .add_load("Moment", Load::Moment { target: roundtrip_set("Tip"), vector: [0.0, 5.0, 0.0], distribution: Distribution::PerNode })
            .unwrap();
        roundtrip_step(&mut model, "Moment", &["Fix"], &["Moment"]);
        model
    }

    fn roundtrip_target_nodes(model: &Model, target: &Target) -> Vec<u32> {
        match target {
            Target::Set(name) => model.resolve_set_partial(name).map(|set| set.nodes).unwrap_or_default(),
            Target::ReferencePoint(_) => Vec::new(),
        }
    }

    fn roundtrip_boundary_signatures(model: &Model, names: &BTreeSet<String>) -> BTreeSet<String> {
        names
            .iter()
            .filter_map(|name| model.boundary_condition(name))
            .map(|boundary| match boundary {
                BoundaryCondition::Fixed { target } => format!("fixed:{:?}", roundtrip_target_nodes(model, target)),
                BoundaryCondition::Symmetry { target, plane } => {
                    format!("symmetry:{plane:?}:{:?}", roundtrip_target_nodes(model, target))
                }
                BoundaryCondition::Displacement { target, values } => {
                    format!("displacement:{:?}:{values:?}", roundtrip_target_nodes(model, target))
                }
                BoundaryCondition::Temperature { target, value } => {
                    format!("temperature:{:?}:{value:?}", roundtrip_target_nodes(model, target))
                }
            })
            .collect()
    }

    fn roundtrip_load_signatures(model: &Model, names: &BTreeSet<String>) -> BTreeSet<String> {
        let mut totals: BTreeMap<(bool, u32), [f64; 3]> = BTreeMap::new();
        for load in names.iter().filter_map(|name| model.load(name)) {
            let (moment, target, vector, distribution) = match load {
                Load::Force { target, vector, distribution } => (false, target, vector, distribution),
                Load::Moment { target, vector, distribution } => (true, target, vector, distribution),
                _ => continue,
            };
            let nodes = roundtrip_target_nodes(model, target);
            for node in &nodes {
                let entry = totals.entry((moment, *node)).or_default();
                match distribution {
                    Distribution::PerNode => entry.iter_mut().zip(vector).for_each(|(out, value)| *out += value),
                    Distribution::Total => {
                        let share = 1.0 / nodes.len() as f64;
                        entry.iter_mut().zip(vector).for_each(|(out, value)| *out += value * share);
                    }
                }
            }
        }
        totals.into_iter().map(|((moment, node), vector)| format!("{moment}:{node}:{vector:?}")).collect()
    }

    fn roundtrip_history_lines(text: &str) -> Vec<String> {
        let lines: Vec<&str> = text.lines().collect();
        let start = lines.iter().position(|line| line.starts_with("*BOUNDARY") || line.starts_with("*STEP")).unwrap_or(lines.len());
        lines[start..].iter().filter(|line| !line.trim_start().starts_with("**")).map(|line| (*line).to_string()).collect()
    }

    fn assert_writer_history_roundtrip(label: &str, model: &Model) -> Read {
        let written = crate::abaqus::write(model, Preset::MKgS).unwrap_or_else(|error| panic!("{label}: writer refused input: {error}"));
        let imported = read_analysis_limited(&written, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(imported.model.steps().len(), model.steps().len(), "{label}: writer steps must all be imported");
        let history_unread: Vec<_> = imported
            .unknown_all
            .iter()
            .filter(|unknown| {
                !matches!(unknown.keyword.as_str(), "OUTPUT" | "NODEOUTPUT" | "ELEMENTOUTPUT" | "NODEPRINT" | "ELPRINT" | "ENERGYOUTPUT" | "PREPRINT")
            })
            .map(|unknown| format!("{}@{}", unknown.keyword, unknown.line))
            .collect();
        assert!(history_unread.is_empty(), "{label}: supported writer history stayed unread: {history_unread:?}");
        assert_eq!(imported.model.nodes(), model.nodes(), "{label}: nodes");
        assert_eq!(imported.model.elements(), model.elements(), "{label}: elements");
        assert_eq!(imported.model.materials(), model.materials(), "{label}: materials");
        assert_eq!(imported.model.sections().len(), model.sections().len(), "{label}: section count");
        for (original, imported_step) in model.steps().iter().zip(imported.model.steps()) {
            assert_eq!(original.name, imported_step.name, "{label}: step name");
            assert_eq!(format!("{:?}", original.procedure), format!("{:?}", imported_step.procedure), "{label}: procedure");
            assert_eq!(
                roundtrip_boundary_signatures(model, &original.boundary_conditions),
                roundtrip_boundary_signatures(&imported.model, &imported_step.boundary_conditions),
                "{label}: boundary signatures",
            );
            assert_eq!(
                roundtrip_load_signatures(model, &original.loads),
                roundtrip_load_signatures(&imported.model, &imported_step.loads),
                "{label}: load signatures",
            );
        }
        let rewritten = crate::abaqus::write(&imported.model, Preset::MKgS).unwrap_or_else(|error| panic!("{label}: rewrite refused: {error}"));
        assert_eq!(roundtrip_history_lines(&rewritten), roundtrip_history_lines(&written), "{label}: history text changed after read/write");
        imported
    }

    #[test]
    fn writer_history_roundtrip_matrix_preserves_each_targeted_object() {
        let mut two_fixed = roundtrip_base();
        two_fixed.add_boundary_condition("ClampLeft", BoundaryCondition::Fixed { target: roundtrip_set("Left") }).unwrap();
        two_fixed.add_boundary_condition("ClampRight", BoundaryCondition::Fixed { target: roundtrip_set("Right") }).unwrap();
        two_fixed.add_load("Push", roundtrip_force("Mid", 5.0)).unwrap();
        roundtrip_step(&mut two_fixed, "S1", &["ClampLeft", "ClampRight"], &["Push"]);

        let mut clamp_and_symmetry = roundtrip_base();
        clamp_and_symmetry.add_boundary_condition("Clamp", BoundaryCondition::Fixed { target: roundtrip_set("Left") }).unwrap();
        clamp_and_symmetry.add_boundary_condition("Sym", BoundaryCondition::Symmetry { target: roundtrip_set("SymX"), plane: Symmetry::X }).unwrap();
        clamp_and_symmetry.add_load("Push", roundtrip_force("Right", 5.0)).unwrap();
        roundtrip_step(&mut clamp_and_symmetry, "S1", &["Clamp", "Sym"], &["Push"]);

        let mut two_symmetries = roundtrip_base();
        two_symmetries.add_boundary_condition("Clamp", BoundaryCondition::Fixed { target: roundtrip_set("Right") }).unwrap();
        two_symmetries.add_boundary_condition("SymA", BoundaryCondition::Symmetry { target: roundtrip_set("SymX"), plane: Symmetry::X }).unwrap();
        two_symmetries.add_boundary_condition("SymB", BoundaryCondition::Symmetry { target: roundtrip_set("SymY"), plane: Symmetry::Y }).unwrap();
        two_symmetries.add_load("Push", roundtrip_force("Mid", 5.0)).unwrap();
        roundtrip_step(&mut two_symmetries, "S1", &["Clamp", "SymA", "SymB"], &["Push"]);

        let mut added_displacement = roundtrip_base();
        added_displacement.add_boundary_condition("Clamp", BoundaryCondition::Fixed { target: roundtrip_set("Left") }).unwrap();
        added_displacement
            .add_boundary_condition(
                "Pull",
                BoundaryCondition::Displacement { target: roundtrip_set("Right"), values: [Some(0.001), None, None, None, None, None] },
            )
            .unwrap();
        added_displacement.add_load("Push", roundtrip_force("Mid", 5.0)).unwrap();
        roundtrip_step(&mut added_displacement, "S1", &["Clamp"], &["Push"]);
        roundtrip_step(&mut added_displacement, "S2", &["Clamp", "Pull"], &[]);

        let mut released_fixed = roundtrip_base();
        released_fixed.add_boundary_condition("Clamp", BoundaryCondition::Fixed { target: roundtrip_set("Left") }).unwrap();
        released_fixed.add_boundary_condition("Hold", BoundaryCondition::Fixed { target: roundtrip_set("Right") }).unwrap();
        released_fixed.add_load("Push", roundtrip_force("Mid", 5.0)).unwrap();
        roundtrip_step(&mut released_fixed, "S1", &["Clamp", "Hold"], &["Push"]);
        roundtrip_step(&mut released_fixed, "S2", &["Clamp"], &["Push"]);

        let mut added_symmetry = roundtrip_base();
        added_symmetry.add_boundary_condition("Clamp", BoundaryCondition::Fixed { target: roundtrip_set("Right") }).unwrap();
        added_symmetry.add_boundary_condition("SymA", BoundaryCondition::Symmetry { target: roundtrip_set("SymX"), plane: Symmetry::X }).unwrap();
        added_symmetry.add_boundary_condition("SymB", BoundaryCondition::Symmetry { target: roundtrip_set("SymY"), plane: Symmetry::Y }).unwrap();
        added_symmetry.add_load("Push", roundtrip_force("Mid", 5.0)).unwrap();
        roundtrip_step(&mut added_symmetry, "S1", &["Clamp", "SymA"], &["Push"]);
        roundtrip_step(&mut added_symmetry, "S2", &["Clamp", "SymA", "SymB"], &["Push"]);

        let mut shared_load_nodes = roundtrip_base();
        shared_load_nodes.add_boundary_condition("Clamp", BoundaryCondition::Fixed { target: roundtrip_set("Left") }).unwrap();
        shared_load_nodes.add_load("A", roundtrip_force("Tip", 10.0)).unwrap();
        shared_load_nodes.add_load("B", roundtrip_force("Right", 5.0)).unwrap();
        roundtrip_step(&mut shared_load_nodes, "S1", &["Clamp"], &["A", "B"]);

        let mut frequency_boundary = roundtrip_base();
        frequency_boundary.add_boundary_condition("Clamp", BoundaryCondition::Fixed { target: roundtrip_set("Left") }).unwrap();
        frequency_boundary.add_boundary_condition("Modal", BoundaryCondition::Fixed { target: roundtrip_set("Right") }).unwrap();
        roundtrip_step(&mut frequency_boundary, "Base", &["Clamp"], &[]);
        roundtrip_step_with(
            &mut frequency_boundary,
            "Modes",
            Procedure::Frequency { modes: 3, min_hz: Some(1.0), max_hz: Some(4.0) },
            &["Clamp", "Modal"],
            &[],
        );

        let mut first_buckling = roundtrip_base();
        first_buckling.add_boundary_condition("Clamp", BoundaryCondition::Fixed { target: roundtrip_set("Left") }).unwrap();
        first_buckling.add_load("Push", roundtrip_force("Mid", 5.0)).unwrap();
        roundtrip_step_with(&mut first_buckling, "Buckling", Procedure::Buckle { modes: 2, shift: None }, &["Clamp"], &["Push"]);

        for (label, model) in [
            ("two fixed targets", two_fixed),
            ("clamp and symmetry", clamp_and_symmetry),
            ("two shared-node symmetries", two_symmetries),
            ("added displacement", added_displacement),
            ("released fixed", released_fixed),
            ("added symmetry", added_symmetry),
            ("shared-node loads", shared_load_nodes),
            ("frequency boundary", frequency_boundary),
            ("first buckling", first_buckling),
            ("beam moment", roundtrip_moment_model()),
        ] {
            let imported = assert_writer_history_roundtrip(label, &model);
            if label == "two fixed targets" {
                let targets: BTreeSet<Vec<u32>> = imported.model.steps()[0]
                    .boundary_conditions
                    .iter()
                    .filter_map(|name| imported.model.boundary_condition(name))
                    .filter_map(|boundary| match boundary {
                        BoundaryCondition::Fixed { target } => Some(roundtrip_target_nodes(&imported.model, target)),
                        _ => None,
                    })
                    .collect();
                assert_eq!(targets, [vec![1], vec![4]].into_iter().collect());
            }
        }
    }

    #[test]
    fn reads_writer_style_static_frequency_history_and_rejects_loaded_buckling() {
        let deck = format!(
            "{MODEL}*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=Static, NLGEOM=NO\n*STATIC\n1.0, 1.0\n*BOUNDARY\nRoot, ENCASTRE\n*CLOAD, OP=NEW\n2, 1, 10\n*END STEP\n*STEP, NAME=Modes\n*FREQUENCY, EIGENSOLVER=LANCZOS\n3\n*BOUNDARY\nRoot, ENCASTRE\n*END STEP\n*STEP, NAME=Buckle\n*BUCKLE, EIGENSOLVER=LANCZOS\n2\n*BOUNDARY\nRoot, ENCASTRE\n*CLOAD\n2, 1, 10\n*END STEP\n"
        );
        let got = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 2);
        assert!(matches!(got.model.steps()[0].procedure, Procedure::Static {}));
        assert!(matches!(got.model.steps()[1].procedure, Procedure::Frequency { modes: 3, .. }));
        assert!(got.model.steps()[0].loads.len() == 1 && got.model.steps()[1].loads.is_empty());
        assert!(got.warnings.iter().any(|warning| warning.contains("buckling after a general analysis step")), "{:?}", got.warnings);
    }

    #[test]
    fn step_names_are_case_insensitive_and_generated_names_avoid_future_names() {
        let duplicate =
            format!("{MODEL}*STEP, NAME=Run\n*STATIC\n*END STEP\n*STEP, NAME=run\n*STATIC\n*END STEP\n*STEP, NAME=After\n*STATIC\n*END STEP\n");
        let got = read_analysis_limited(&duplicate, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().iter().map(|step| step.name.as_str()).collect::<Vec<_>>(), ["Run"]);
        assert!(got.warnings.iter().any(|warning| warning.contains("already exists") && warning.contains("case")), "{:?}", got.warnings);

        let generated = format!("{MODEL}*STEP\n*STATIC\n*END STEP\n*STEP, NAME=importedstep2\n*STATIC\n*END STEP\n");
        let got = read_analysis_limited(&generated, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 2);
        assert_ne!(got.model.steps()[0].name.to_ascii_uppercase(), "IMPORTEDSTEP2");
        assert_eq!(got.model.steps()[1].name, "importedstep2");
        assert!(crate::abaqus::write(&got.model, Preset::MKgS).is_ok());
    }

    #[test]
    fn buckling_after_an_unloaded_general_step_stays_unread() {
        let deck = format!(
            "{MODEL}*STEP, NAME=Base\n*STATIC\n*END STEP\n*STEP, NAME=Buckle\n*BUCKLE\n2\n*END STEP\n*STEP, NAME=After\n*STATIC\n*END STEP\n"
        );
        let got = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 1);
        assert!(got.warnings.iter().any(|warning| warning.contains("buckling after a general analysis step")), "{:?}", got.warnings);
    }

    #[test]
    fn static_writer_roundtrip_keeps_history_lines() {
        let mut read = read_in(MODEL, Preset::MKgS).unwrap().model;
        read.add_set("Clamp", Set { nodes: vec![1], ..Set::default() }).unwrap();
        read.add_boundary_condition("InpAnalysis1", BoundaryCondition::Fixed { target: Target::Set("Clamp".into()) }).unwrap();
        read.add_load(
            "InpAnalysis2",
            Load::Force { target: Target::Set("Tip".into()), vector: [10.0, 0.0, 0.0], distribution: Distribution::PerNode },
        )
        .unwrap();
        let mut step = Step::new("Static", Procedure::Static {});
        step.boundary_conditions.insert("InpAnalysis1".into());
        step.loads.insert("InpAnalysis2".into());
        read.add_step(step, None).unwrap();
        let written = crate::abaqus::write(&read, Preset::MKgS).unwrap();
        let imported = read_analysis_limited(&written, Preset::MKgS, Limits::default()).unwrap();
        let rewritten = crate::abaqus::write(&imported.model, Preset::MKgS).unwrap();
        assert_eq!(written, rewritten);
        assert_eq!(imported.model.elements(), read.elements());
        assert_eq!(imported.model.steps().len(), read.steps().len());
        assert_eq!(imported.model, read);
    }

    #[test]
    fn frequency_writer_roundtrip_keeps_the_modal_step() {
        let mut read = read_in(MODEL, Preset::MKgS).unwrap().model;
        read.add_boundary_condition("InpAnalysis1", BoundaryCondition::Fixed { target: Target::Set("Root".into()) }).unwrap();
        let mut step = Step::new("Modes", Procedure::Frequency { modes: 3, min_hz: Some(1.0), max_hz: Some(4.0) });
        step.boundary_conditions.insert("InpAnalysis1".into());
        read.add_step(step, None).unwrap();
        let written = crate::abaqus::write(&read, Preset::MKgS).unwrap();
        let imported = read_analysis_limited(&written, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(crate::abaqus::write(&imported.model, Preset::MKgS).unwrap(), written);
        assert_eq!(imported.model, read);
    }

    #[test]
    fn buckle_writer_roundtrip_keeps_the_concentrated_load() {
        let mut read = read_in(MODEL, Preset::MKgS).unwrap().model;
        read.add_boundary_condition("InpAnalysis1", BoundaryCondition::Fixed { target: Target::Set("Root".into()) }).unwrap();
        read.add_load(
            "InpAnalysis2",
            Load::Force { target: Target::Set("Tip".into()), vector: [10.0, 0.0, 0.0], distribution: Distribution::PerNode },
        )
        .unwrap();
        let mut step = Step::new("Buckling", Procedure::Buckle { modes: 2, shift: None });
        step.boundary_conditions.insert("InpAnalysis1".into());
        step.loads.insert("InpAnalysis2".into());
        read.add_step(step, None).unwrap();
        let written = crate::abaqus::write(&read, Preset::MKgS).unwrap();
        let imported = read_analysis_limited(&written, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(crate::abaqus::write(&imported.model, Preset::MKgS).unwrap(), written);
        assert_eq!(imported.model, read);
    }

    #[test]
    fn multi_node_targets_keep_per_node_values_without_multiplication() {
        let deck =
            format!("{MULTI_MODEL}*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=Loads\n*STATIC\n*CLOAD\nTip, 1, 10\nTip, 2, 20\nTip, 5, 3\n*END STEP\n");
        let got = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.boundary_conditions().len(), 1);
        assert_eq!(got.model.loads().len(), 2);
        assert!(got.model.loads().values().any(|load| matches!(
            load,
            Load::Force { target: Target::Set(set), vector, distribution: Distribution::PerNode }
                if set == "Tip" && *vector == [10.0, 20.0, 0.0]
        )));
        assert!(got.model.loads().values().any(|load| matches!(
            load,
            Load::Moment { target: Target::Set(set), vector, distribution: Distribution::PerNode }
                if set == "Tip" && *vector == [0.0, 3.0, 0.0]
        )));
    }

    #[test]
    fn state_is_per_target_and_dof_across_mod_and_new_steps() {
        let deck = format!(
            "{MULTI_MODEL}*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=S1\n*STATIC\n*CLOAD\nTip, 1, 10\n*END STEP\n*STEP, NAME=S2\n*STATIC\n*CLOAD\n4, 1, 5\n*END STEP\n*STEP, NAME=S3\n*STATIC\n*CLOAD, OP=NEW\n3, 2, 7\n*END STEP\n"
        );
        let got = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 3);
        assert_eq!(got.model.steps()[0].loads.len(), 1);
        assert_eq!(got.model.steps()[1].loads.len(), 2);
        assert_eq!(got.model.steps()[2].loads.len(), 1);
        let s2: Vec<_> = got.model.steps()[1].loads.iter().filter_map(|name| got.model.load(name)).collect();
        assert!(
            s2.iter().any(
                |load| matches!(load, Load::Force { target: Target::Set(target), vector, .. } if target == "Tip" && *vector == [10.0, 0.0, 0.0])
            )
        );
        assert!(s2.iter().any(|load| matches!(load, Load::Force { vector, .. } if *vector == [5.0, 0.0, 0.0])));
        let s3: Vec<_> = got.model.steps()[2].loads.iter().filter_map(|name| got.model.load(name)).collect();
        assert!(s3.iter().all(|load| matches!(load, Load::Force { vector, .. } if *vector == [0.0, 7.0, 0.0])));
    }

    #[test]
    fn linear_perturbation_changes_do_not_leak_to_the_next_static_step() {
        let deck = format!(
            "{MULTI_MODEL}*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=Base\n*STATIC\n*END STEP\n*STEP, NAME=Modes\n*FREQUENCY\n2\n*BOUNDARY\nTip, 1, 1, 0.5\n*END STEP\n*STEP, NAME=After\n*STATIC\n*END STEP\n"
        );
        let got = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 3);
        assert_eq!(got.model.steps()[0].boundary_conditions, got.model.steps()[2].boundary_conditions);
        assert_ne!(got.model.steps()[1].boundary_conditions, got.model.steps()[2].boundary_conditions);
    }

    #[test]
    fn buckling_load_changes_do_not_replace_the_general_step_state() {
        let deck = format!(
            "{MULTI_MODEL}*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=Base\n*STATIC\n*CLOAD\n4, 1, 10\n*END STEP\n*STEP, NAME=Buckle\n*BUCKLE\n2\n*CLOAD\n4, 1, -1\n*END STEP\n*STEP, NAME=After\n*STATIC\n*END STEP\n"
        );
        let got = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 1);
        assert!(got.warnings.iter().any(|warning| warning.contains("buckling after a general analysis step")), "{:?}", got.warnings);
    }

    fn scale_deck(nodes: usize, displacement: bool) -> String {
        let mut deck = String::from("*NODE\n");
        for node in 1..=nodes {
            deck.push_str(&format!("{node},{node},0,0\n"));
        }
        deck.push_str("*ELEMENT,TYPE=T3D2,ELSET=Bar\n");
        for element in 1..nodes {
            deck.push_str(&format!("{element},{element},{}\n", element + 1));
        }
        deck.push_str(&format!("*NSET,NSET=All,GENERATE\n1,{nodes},1\n"));
        deck.push_str("*MATERIAL,NAME=Steel\n*ELASTIC\n210000,0.3\n*DENSITY\n1\n*SOLID SECTION,ELSET=Bar,MATERIAL=Steel\n1\n*BOUNDARY\n1,1,3\n*STEP,NAME=Scale\n*STATIC\n");
        if displacement {
            deck.push_str("*BOUNDARY\n");
            for node in 2..=nodes {
                deck.push_str(&format!("{node},2,2,{}\n", 1e-6 * node as f64));
            }
        } else {
            deck.push_str("*CLOAD\n");
            for node in 2..=nodes {
                deck.push_str(&format!("{node},1,{}\n", 0.001 * node as f64));
            }
        }
        deck.push_str("*END STEP\n");
        deck
    }

    #[test]
    fn distinct_node_values_scale_without_quadratic_grouping() {
        for displacement in [false, true] {
            let got = read_analysis_limited(&scale_deck(2_048, displacement), Preset::MKgS, Limits::default()).unwrap();
            if displacement {
                assert_eq!(got.model.boundary_conditions().len(), 2_048);
            } else {
                assert_eq!(got.model.loads().len(), 2_047);
            }
        }
    }

    #[test]
    #[ignore = "heavy: exercises distinct-value grouping at production scale"]
    fn distinct_node_values_large_scale_regression() {
        for displacement in [false, true] {
            let got = read_analysis_limited(&scale_deck(32_000, displacement), Preset::MKgS, Limits::default()).unwrap();
            assert_eq!(got.model.steps().len(), 1);
        }
    }

    #[test]
    fn history_members_share_the_global_limit_across_steps() {
        let mut deck = MODEL.to_string();
        for dof in 1..=6 {
            deck.push_str(&format!("*STEP, NAME=S{dof}\n*STATIC\n*CLOAD\n2, {dof}, {dof}\n*END STEP\n"));
        }
        let limits = Limits { members: 20, ..Limits::default() };
        let got = read_analysis_limited(&deck, Preset::MKgS, limits).unwrap();
        assert_eq!(got.model.steps().len(), 5);
        assert!(got.warnings.iter().any(|warning| warning.contains("history member limit")), "{:?}", got.warnings);
        assert!(got.unknown.iter().any(|unknown| unknown.keyword == "STEP"));
    }

    #[test]
    fn duplicate_targets_are_reported_without_guessing() {
        let additive = format!("{MULTI_MODEL}*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=Same\n*STATIC\n*CLOAD\nTip, 1, 10\n4, 1, 5\n*END STEP\n");
        let got = read_analysis_limited(&additive, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 1);
        assert!(got.model.loads().values().any(|load| matches!(load, Load::Force { vector, .. } if *vector == [10.0, 0.0, 0.0])));
        assert!(got.model.loads().values().any(|load| matches!(load, Load::Force { vector, .. } if *vector == [5.0, 0.0, 0.0])));

        for history in [
            "*STEP, NAME=Same\n*STATIC\n*CLOAD\nTip, 1, 10\nTip, 1, 5\n*END STEP\n",
            "*STEP, NAME=Same\n*STATIC\n*BOUNDARY\nTip, 1, 1, 0.001\nMid, 1, 1, 0.002\n*END STEP\n",
        ] {
            let deck = format!("{MULTI_MODEL}*BOUNDARY\nRoot, ENCASTRE\n{history}");
            let got = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
            assert!(got.model.steps().is_empty());
            assert!(got.warnings.iter().any(|warning| warning.contains("duplicate") || warning.contains("overlap")), "{:?}", got.warnings);
            assert_model_data_unchanged(&deck);
        }
    }

    #[test]
    fn identical_conditions_and_loads_are_shared_between_steps() {
        let deck = format!(
            "{MODEL}*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=S1\n*STATIC\n*CLOAD\n2, 1, 10\n*END STEP\n*STEP, NAME=F\n*FREQUENCY\n2\n*END STEP\n*STEP, NAME=S2\n*STATIC\n*END STEP\n"
        );
        let got = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 3);
        assert_eq!(got.model.boundary_conditions().len(), 1);
        assert_eq!(got.model.loads().len(), 1);
        assert_eq!(got.model.steps()[0].boundary_conditions, got.model.steps()[1].boundary_conditions);
        assert_eq!(got.model.steps()[1].boundary_conditions, got.model.steps()[2].boundary_conditions);
        assert_eq!(got.model.steps()[0].loads, got.model.steps()[2].loads);
    }

    #[test]
    fn history_mapping_errors_stay_unread_and_preserve_model_data() {
        for bad in [
            format!("{MODEL}*STEP,NAME=Bad\n*STATIC\n1,1\n*BOUNDARY\nMissing,1,1,0\n*END STEP\n"),
            format!("{MODEL}*STEP,NAME=Bad\n*FREQUENCY\n0\n*END STEP\n"),
            format!("{MODEL}*STEP,NAME=Bad\n*STATIC\n1e400,1\n*END STEP\n"),
            format!("{MODEL}*STEP,NAME=Bad\n*STATIC\n*CLOAD\n2,0,1\n*END STEP\n"),
            format!("{MODEL}*STEP,NAME=Bad\n*STATIC\n*CLOAD\n999999999,1,1\n*END STEP\n"),
            format!("{MODEL}*STEP,NAME=Bad\n*STATIC\n*CLOAD\n,\n*END STEP\n"),
        ] {
            let got = read_analysis_limited(&bad, Preset::MKgS, Limits::default()).unwrap();
            assert!(got.model.steps().is_empty());
            assert!(got.warnings.iter().all(|warning| warning.contains("line")), "{:?}", got.warnings);
            assert_model_data_unchanged(&bad);
        }
        let huge_modes = format!("{MODEL}*STEP,NAME=Bad\n*FREQUENCY\n100001\n*END STEP\n");
        let got = read_analysis_limited(&huge_modes, Preset::MKgS, Limits::default()).unwrap();
        assert!(got.model.steps().is_empty() && got.warnings.iter().any(|warning| warning.contains("mode count")), "{:?}", got.warnings);
        let duplicate = format!("{MODEL}*STEP,NAME=Same\n*STATIC\n*END STEP\n*STEP,NAME=Same\n*STATIC\n*END STEP\n");
        let got = read_analysis_limited(&duplicate, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 1);
        assert!(got.warnings.iter().any(|warning| warning.contains("already exists") && warning.contains("line")));
        assert_model_data_unchanged(&duplicate);
        let outside_unsupported = format!("{MODEL}*BOUNDARY\nRoot, XASYMM\n*STEP,NAME=Bad\n*STATIC\n*END STEP\n");
        let got = read_analysis_limited(&outside_unsupported, Preset::MKgS, Limits::default()).unwrap();
        assert!(got.model.steps().is_empty() && got.warnings.iter().any(|warning| warning.contains("XASYMM")));
        assert_model_data_unchanged(&outside_unsupported);
        let unsupported = format!("{MODEL}*STEP,NAME=Bad\n*STATIC\n*DYNAMIC\n*END STEP\n");
        assert_model_data_unchanged(&unsupported);
        let good = format!("{MODEL}*STEP,NAME=Good\n*STATIC\n1,1\n*CLOAD\n2,1,2\n*END STEP\n");
        for at in 0..=good.len() {
            let _ = read_analysis_limited(&good[..at], Preset::MKgS, Limits::default());
        }
    }

    #[test]
    fn target_identity_replaces_only_the_same_concentrated_load() {
        let deck = format!(
            "{MULTI_MODEL}*STEP, NAME=S1\n*STATIC\n*CLOAD\nTip, 1, 10\n*END STEP\n*STEP, NAME=S2\n*STATIC\n*CLOAD\nTip, 1, 5\n*END STEP\n*STEP, NAME=S3\n*STATIC\n*CLOAD\n4, 1, 5\n*END STEP\n"
        );
        let got = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
        let first: Vec<_> = got.model.steps()[0].loads.iter().filter_map(|name| got.model.load(name)).collect();
        let second: Vec<_> = got.model.steps()[1].loads.iter().filter_map(|name| got.model.load(name)).collect();
        let third: Vec<_> = got.model.steps()[2].loads.iter().filter_map(|name| got.model.load(name)).collect();
        assert!(first.iter().any(|load| matches!(load, Load::Force { vector, .. } if *vector == [10.0, 0.0, 0.0])));
        assert!(second.iter().any(|load| matches!(load, Load::Force { vector, .. } if *vector == [5.0, 0.0, 0.0])));
        assert!(third.iter().any(|load| matches!(load, Load::Force { vector, .. } if *vector == [5.0, 0.0, 0.0])));
        assert!(third.iter().filter(|load| matches!(load, Load::Force { vector, .. } if *vector == [5.0, 0.0, 0.0])).count() == 2);
    }

    #[test]
    fn zero_boundary_targets_can_overlap_and_named_targets_remain_distinct() {
        let deck = format!(
            "{MULTI_MODEL}*NSET, NSET=SymX\n1,2\n*NSET, NSET=SymY\n2,3\n*BOUNDARY\nRoot, ENCASTRE\nSymX, XSYMM\n*STEP, NAME=S1\n*STATIC\n*CLOAD\n4, 1, 1\n*END STEP\n*STEP, NAME=S2\n*STATIC\n*BOUNDARY\nSymY, YSYMM\n*END STEP\n*STEP, NAME=S3\n*STATIC\n*BOUNDARY\nSymX, XSYMM\n*END STEP\n"
        );
        let got = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 3);
        assert_eq!(got.model.steps()[0].boundary_conditions.len(), 2);
        assert_eq!(got.model.steps()[1].boundary_conditions.len(), 3);
        assert_eq!(got.model.steps()[2].boundary_conditions.len(), 3);
        assert!(
            got.model.boundary_conditions().values().any(|bc| matches!(bc, BoundaryCondition::Fixed { target: Target::Set(name) } if name == "Root"))
        );
        assert!(
            got.model
                .boundary_conditions()
                .values()
                .any(|bc| matches!(bc, BoundaryCondition::Symmetry { target: Target::Set(name), plane: Symmetry::X } if name == "SymX"))
        );
        assert!(
            got.model
                .boundary_conditions()
                .values()
                .any(|bc| matches!(bc, BoundaryCondition::Symmetry { target: Target::Set(name), plane: Symmetry::Y } if name == "SymY"))
        );
    }

    #[test]
    fn named_loads_stay_separate_and_same_name_element_members_are_not_history_nodes() {
        let deck = "*NODE\n1,0,0,0\n2,1,0,0\n3,2,0,0\n4,3,0,0\n*ELEMENT,TYPE=T3D2,ELSET=Bar\n1,1,2\n2,2,3\n3,3,4\n*NSET,NSET=Mix\n1\n*ELSET,ELSET=Mix\n3\n*NSET,NSET=Root\n2\n*MATERIAL,NAME=Steel\n*ELASTIC\n210000,0.3\n*DENSITY\n1\n*SOLID SECTION,ELSET=Bar,MATERIAL=Steel\n1\n*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=S1\n*STATIC\n*BOUNDARY\nMix, 2, 2\n*CLOAD\nMix, 1, 5\n*END STEP\n";
        let got = read_analysis_limited(deck, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 1);
        let load = got.model.loads().values().next().expect("one load");
        assert!(
            matches!(load, Load::Force { target: Target::Set(name), vector, .. } if name.starts_with("InpAnalysis") && *vector == [5.0, 0.0, 0.0])
        );
        let generated = got.model.sets().values().find(|set| set.nodes == [1]);
        assert!(generated.is_some(), "history target must contain only the node-set members");
        let generated_name = got.model.sets().iter().find(|(_, set)| set.nodes == [1]).map(|(name, _)| name).expect("generated history set");
        let boundary_target = got.model.boundary_conditions().values().find_map(|boundary| match boundary {
            BoundaryCondition::Displacement { target: Target::Set(name), .. } => Some(name),
            _ => None,
        });
        let load_target = got.model.loads().values().find_map(|load| match load {
            Load::Force { target: Target::Set(name), .. } => Some(name),
            _ => None,
        });
        assert_eq!(boundary_target, Some(generated_name));
        assert_eq!(load_target, Some(generated_name));
    }

    #[test]
    fn multiple_new_blocks_accumulate_and_mixed_modes_are_rejected() {
        let good = format!(
            "{MODEL}*STEP, NAME=Good\n*STATIC\n*BOUNDARY, OP=NEW\nRoot, ENCASTRE\n*BOUNDARY, OP=NEW\nTip, 1, 1, 0.2\n*CLOAD, OP=NEW\n2, 1, 3\n*CLOAD, OP=NEW\n1, 2, 4\n*END STEP\n"
        );
        let got = read_analysis_limited(&good, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 1);
        assert_eq!(got.model.steps()[0].loads.len(), 2);
        assert_eq!(got.model.steps()[0].boundary_conditions.len(), 2);

        let mixed = format!("{MODEL}*STEP, NAME=Bad\n*STATIC\n*BOUNDARY, OP=NEW\nRoot, ENCASTRE\n*BOUNDARY\nTip, 1, 1, 0.2\n*END STEP\n");
        let got = read_analysis_limited(&mixed, Preset::MKgS, Limits::default()).unwrap();
        assert!(got.model.steps().is_empty());
        assert!(got.warnings.iter().any(|warning| warning.contains("mixes OP=NEW and OP=MOD")), "{:?}", got.warnings);
    }

    #[test]
    fn history_members_are_limited_before_expanding_a_boundary_row() {
        let limits = Limits { members: 16, ..Limits::default() };
        let deck = "*NODE\n1,0,0,0\n2,1,0,0\n3,2,0,0\n4,3,0,0\n5,4,0,0\n6,5,0,0\n7,6,0,0\n8,7,0,0\n*NSET,NSET=Big\n1,2,3,4,5,6,7,8\n*STEP, NAME=TooLarge\n*STATIC\n*BOUNDARY\nBig, ENCASTRE\n*END STEP\n";
        let got = read_analysis_limited(deck, Preset::MKgS, limits).unwrap();
        assert!(got.model.steps().is_empty());
        assert!(got.warnings.iter().any(|warning| warning.contains("analysis membership limit")), "{:?}", got.warnings);
        assert_model_data_unchanged(deck);
    }

    #[test]
    fn generated_names_avoid_case_insensitive_source_names_and_resolved_sets() {
        let deck = format!(
            "{MULTI_MODEL}*NSET, NSET=inpanalysis1\n1\n*ELSET, ELSET=Mix\n1\n*NSET, NSET=Mix\n1\n*STEP, NAME=Names\n*STATIC\n*BOUNDARY\nMix, 1, 1, 0.1\n*CLOAD\n4, 1, 1\n*END STEP\n"
        );
        let got = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 1);
        let mix = got.model.set("Mix").expect("source set is retained");
        assert_eq!(mix.nodes, vec![1]);
        assert_eq!(got.model.boundary_conditions().len(), 1);
        let generated: Vec<_> = got.model.sets().keys().filter(|name| name.to_ascii_lowercase().starts_with("inpanalysis")).collect();
        assert!(generated.iter().any(|name| name.eq_ignore_ascii_case("inpanalysis2")));
        assert!(!got.model.sets().keys().any(|name| name == "InpAnalysis1"));
        assert!(
            matches!(got.model.boundary_conditions().get("InpAnalysis3"), Some(BoundaryCondition::Displacement { target: Target::Set(name), .. }) if name == "InpAnalysis2")
        );
    }

    #[test]
    fn unused_initial_boundaries_remain_unread_but_used_ones_are_consumed() {
        let no_step = format!("{MODEL}*BOUNDARY\nRoot, ENCASTRE\n");
        let got = read_analysis_limited(&no_step, Preset::MKgS, Limits::default()).unwrap();
        assert!(got.model.steps().is_empty());
        assert!(got.unknown.iter().any(|unknown| unknown.keyword == "BOUNDARY"));

        let unsupported = format!("{MODEL}*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=Bad\n*DYNAMIC\n*END STEP\n");
        let got = read_analysis_limited(&unsupported, Preset::MKgS, Limits::default()).unwrap();
        assert!(got.model.steps().is_empty());
        assert!(got.unknown.iter().any(|unknown| unknown.keyword == "BOUNDARY"));

        let used = format!("{MODEL}*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=Good\n*STATIC\n*END STEP\n");
        let got = read_analysis_limited(&used, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 1);
        assert!(!got.unknown.iter().any(|unknown| unknown.keyword == "BOUNDARY"));
    }

    #[test]
    fn unread_count_and_first_hundred_entries_are_rebuilt_after_history_consumption() {
        let mut deck = MODEL.to_string();
        for index in 0..50 {
            deck.push_str(&format!("*UNKNOWNBEFORE{index}\n"));
        }
        deck.push_str("*STEP, NAME=Good\n*STATIC\n*END STEP\n");
        for index in 0..75 {
            deck.push_str(&format!("*UNKNOWNAFTER{index}\n"));
        }
        let got = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 1);
        assert_eq!(got.unknown_count, 125);
        assert_eq!(got.unknown.len(), 100);
        assert!(got.unknown.iter().all(|unknown| unknown.keyword.starts_with("UNKNOWN")));
    }

    #[derive(Clone, Copy)]
    enum InvariantSet {
        SingleNode,
        TwoNodes,
    }

    #[derive(Clone, Copy)]
    enum InvariantPosition {
        Same,
        Later,
        New,
        Persistent,
    }

    #[derive(Clone)]
    struct ReplayLoadRow {
        target: String,
        nodes: Vec<u32>,
        dof: u8,
        value: f64,
    }

    fn invariant_case(
        set_kind: InvariantSet,
        set_first: bool,
        same_value: bool,
        position: InvariantPosition,
        moment: bool,
    ) -> (String, Vec<(bool, Vec<ReplayLoadRow>)>) {
        let set_nodes = match set_kind {
            InvariantSet::SingleNode => vec![4],
            InvariantSet::TwoNodes => vec![3, 4],
        };
        let set_block = format!("*NSET,NSET=Tip\n{}\n", set_nodes.iter().map(u32::to_string).collect::<Vec<_>>().join(","));
        let set_row = ReplayLoadRow { target: "Tip".into(), nodes: set_nodes.clone(), dof: if moment { 5 } else { 1 }, value: 5.0 };
        let node_rows: Vec<_> = if matches!(set_kind, InvariantSet::TwoNodes) {
            set_nodes
                .iter()
                .map(|node| ReplayLoadRow {
                    target: node.to_string(),
                    nodes: vec![*node],
                    dof: if moment { 5 } else { 1 },
                    value: if same_value { 5.0 } else { 7.0 },
                })
                .collect()
        } else {
            vec![ReplayLoadRow { target: "4".into(), nodes: vec![4], dof: if moment { 5 } else { 1 }, value: if same_value { 5.0 } else { 7.0 } }]
        };
        let (first_rows, second_rows) = if set_first { (vec![set_row.clone()], node_rows) } else { (node_rows, vec![set_row.clone()]) };
        let format_rows = |rows: &[ReplayLoadRow]| rows.iter().map(|row| format!("{}, {}, {}\n", row.target, row.dof, row.value)).collect::<String>();
        let mut history = Vec::new();
        let mut steps = Vec::new();
        match position {
            InvariantPosition::Same => {
                let rows = first_rows.iter().chain(second_rows.iter()).cloned().collect();
                history.push((false, rows));
                steps.push(format!(
                    "*STEP, NAME=S1\n*STATIC\n*CLOAD\n{}*END STEP\n",
                    format_rows(&first_rows.iter().chain(second_rows.iter()).cloned().collect::<Vec<_>>())
                ));
            }
            InvariantPosition::Later => {
                history.push((false, first_rows.clone()));
                history.push((false, second_rows.clone()));
                steps.push(format!("*STEP, NAME=S1\n*STATIC\n*CLOAD\n{}*END STEP\n", format_rows(&first_rows)));
                steps.push(format!("*STEP, NAME=S2\n*STATIC\n*CLOAD\n{}*END STEP\n", format_rows(&second_rows)));
            }
            InvariantPosition::New => {
                history.push((false, first_rows.clone()));
                let rows = first_rows.iter().chain(second_rows.iter()).cloned().collect::<Vec<_>>();
                history.push((true, rows.clone()));
                steps.push(format!("*STEP, NAME=S1\n*STATIC\n*CLOAD\n{}*END STEP\n", format_rows(&first_rows)));
                steps.push(format!("*STEP, NAME=S2\n*STATIC\n*CLOAD, OP=NEW\n{}*END STEP\n", format_rows(&rows)));
            }
            InvariantPosition::Persistent => {
                history.push((false, first_rows.clone()));
                history.push((false, second_rows.clone()));
                history.push((false, Vec::new()));
                steps.push(format!("*STEP, NAME=S1\n*STATIC\n*CLOAD\n{}*END STEP\n", format_rows(&first_rows)));
                steps.push(format!("*STEP, NAME=S2\n*STATIC\n*CLOAD\n{}*END STEP\n", format_rows(&second_rows)));
                steps.push("*STEP, NAME=S3\n*STATIC\n*END STEP\n".into());
            }
        }
        let model = if moment { INVARIANT_BEAM_MODEL } else { INVARIANT_MODEL };
        let deck = format!("{model}{set_block}*BOUNDARY\nRoot, ENCASTRE\n{}", steps.concat());
        (deck, history)
    }

    fn replay_load_totals(history: &[(bool, Vec<ReplayLoadRow>)]) -> Vec<BTreeMap<(bool, u32), [f64; 3]>> {
        let mut active: BTreeMap<(String, u8), (Vec<u32>, f64)> = BTreeMap::new();
        let mut totals = Vec::new();
        for (clear, rows) in history {
            if *clear {
                active.clear();
            }
            for row in rows {
                active.insert((row.target.clone(), row.dof), (row.nodes.clone(), row.value));
            }
            let mut step = BTreeMap::new();
            for ((_, dof), (nodes, value)) in &active {
                let moment = *dof > 3;
                let component = usize::from((*dof - 1) % 3);
                for node in nodes {
                    step.entry((moment, *node)).or_insert([0.0; 3])[component] += value;
                }
            }
            totals.push(step);
        }
        totals
    }

    fn imported_load_totals(model: &Model, step: usize) -> BTreeMap<(bool, u32), [f64; 3]> {
        let mut totals = BTreeMap::new();
        for name in &model.steps()[step].loads {
            let (moment, target, vector, distribution) = match model.load(name) {
                Some(Load::Force { target, vector, distribution }) => (false, target, vector, distribution),
                Some(Load::Moment { target, vector, distribution }) => (true, target, vector, distribution),
                _ => continue,
            };
            let nodes = roundtrip_target_nodes(model, target);
            for node in &nodes {
                let entry = totals.entry((moment, *node)).or_insert([0.0; 3]);
                match distribution {
                    Distribution::PerNode => entry.iter_mut().zip(vector).for_each(|(out, value)| *out += value),
                    Distribution::Total => {
                        let share = 1.0 / nodes.len() as f64;
                        entry.iter_mut().zip(vector).for_each(|(out, value)| *out += value * share);
                    }
                }
            }
        }
        totals
    }

    fn assert_totals_match(label: &str, expected: &BTreeMap<(bool, u32), [f64; 3]>, actual: &BTreeMap<(bool, u32), [f64; 3]>) {
        assert_eq!(actual.keys().collect::<Vec<_>>(), expected.keys().collect::<Vec<_>>(), "{label}: affected nodes");
        for (key, values) in expected {
            let got = actual.get(key).expect("key checked above");
            for (expected, actual) in values.iter().zip(got) {
                assert!((expected - actual).abs() < 1e-12, "{label}: {key:?}: expected {expected}, got {actual}");
            }
        }
    }

    fn assert_analysis_write_roundtrip(label: &str, first: &Read) {
        if first.model.steps().is_empty() {
            return;
        }
        let written =
            crate::abaqus::write(&first.model, Preset::MKgS).unwrap_or_else(|error| panic!("{label}: writer refused a readable history: {error}"));
        let reread = read_analysis_limited(&written, Preset::MKgS, Limits::default())
            .unwrap_or_else(|error| panic!("{label}: writer output was not readable: {error}"));
        assert_eq!(reread.model.steps().len(), first.model.steps().len(), "{label}: step count");
        for (index, (original, imported)) in first.model.steps().iter().zip(reread.model.steps()).enumerate() {
            assert_eq!(imported.name, original.name, "{label}: step {index} name");
            assert_eq!(format!("{:?}", imported.procedure), format!("{:?}", original.procedure), "{label}: step {index} procedure");
            assert_eq!(
                roundtrip_boundary_signatures(&first.model, &original.boundary_conditions),
                roundtrip_boundary_signatures(&reread.model, &imported.boundary_conditions),
                "{label}: step {index} boundary signatures",
            );
            assert_totals_match(
                &format!("{label}: step {index} load totals"),
                &imported_load_totals(&first.model, index),
                &imported_load_totals(&reread.model, index),
            );
        }
    }

    #[test]
    fn load_and_boundary_targets_preserve_per_node_sums_by_input_spelling() {
        for set_kind in [InvariantSet::SingleNode, InvariantSet::TwoNodes] {
            for set_first in [true, false] {
                for same_value in [true, false] {
                    for position in [InvariantPosition::Same, InvariantPosition::Later, InvariantPosition::New, InvariantPosition::Persistent] {
                        for moment in [false, true] {
                            let (deck, history) = invariant_case(set_kind, set_first, same_value, position, moment);
                            let expected = replay_load_totals(&history);
                            let imported = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
                            assert_eq!(imported.model.steps().len(), expected.len());
                            assert!(imported.unknown_all.is_empty(), "history was left unread: {:?}", imported.unknown_all);
                            for (step, expected) in expected.iter().enumerate() {
                                assert_totals_match("read", expected, &imported_load_totals(&imported.model, step));
                            }
                            let rewritten = crate::abaqus::write(&imported.model, Preset::MKgS).unwrap();
                            let reread = read_analysis_limited(&rewritten, Preset::MKgS, Limits::default()).unwrap();
                            assert_eq!(reread.model.steps().len(), expected.len());
                            for (step, expected) in expected.iter().enumerate() {
                                assert_totals_match("write/read", expected, &imported_load_totals(&reread.model, step));
                            }
                        }
                    }
                }
            }
        }

        let deck = format!(
            "{INVARIANT_MODEL}*NSET,NSET=Tip\n3,4\n*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=S1\n*STATIC\n*BOUNDARY\nTip, ENCASTRE\n4, 1, 1, 0\n*END STEP\n*STEP, NAME=S2\n*STATIC\n*BOUNDARY, OP=NEW\nTip, ENCASTRE\n4, 1, 1, 0\n*END STEP\n*STEP, NAME=S3\n*STATIC\n*END STEP\n"
        );
        let imported = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(imported.model.steps().len(), 3);
        let expected = |step: usize| {
            let mut values = BTreeMap::new();
            if step == 0 {
                for node in [1, 3, 4] {
                    for dof in 1..=6 {
                        values.insert((node, dof), 0.0);
                    }
                }
            } else {
                for node in [3, 4] {
                    for dof in 1..=6 {
                        values.insert((node, dof), 0.0);
                    }
                }
            }
            values
        };
        fn boundary_totals(model: &Model, step: usize) -> BTreeMap<(u32, u8), f64> {
            let mut totals = BTreeMap::new();
            for name in &model.steps()[step].boundary_conditions {
                let Some(boundary) = model.boundary_condition(name) else { continue };
                let (target, values) = match boundary {
                    BoundaryCondition::Fixed { target } => (target, [Some(0.0); 6]),
                    BoundaryCondition::Symmetry { target, plane } => {
                        let mut values = [None; 6];
                        for dof in plane.dofs().list() {
                            values[usize::from(dof - 1)] = Some(0.0);
                        }
                        (target, values)
                    }
                    BoundaryCondition::Displacement { target, values } => (target, *values),
                    BoundaryCondition::Temperature { target, value } => (target, [Some(*value), None, None, None, None, None]),
                };
                for node in roundtrip_target_nodes(model, target) {
                    for (index, value) in values.into_iter().enumerate() {
                        if let Some(value) = value {
                            totals.insert((node, (index + 1) as u8), value);
                        }
                    }
                }
            }
            totals
        }
        for step in 0..3 {
            assert_eq!(boundary_totals(&imported.model, step), expected(step));
        }

        let mixed = format!(
            "{MULTI_MODEL}*NSET,NSET=Mix\n4\n*ELSET,ELSET=Mix\n1\n*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=S1\n*STATIC\n*CLOAD\nMix, 1, 5\n*END STEP\n*STEP, NAME=S2\n*STATIC\n*CLOAD\n4, 1, 5\n*END STEP\n"
        );
        let mixed_read = read_analysis_limited(&mixed, Preset::MKgS, Limits::default()).unwrap();
        assert!(mixed_read.unknown_all.is_empty(), "{:?}", mixed_read.unknown_all);
        assert_totals_match(
            "same node and element set, step 1",
            &BTreeMap::from([((false, 4), [5.0, 0.0, 0.0])]),
            &imported_load_totals(&mixed_read.model, 0),
        );
        assert_totals_match(
            "same node and element set, step 2",
            &BTreeMap::from([((false, 4), [10.0, 0.0, 0.0])]),
            &imported_load_totals(&mixed_read.model, 1),
        );
        let mixed_reread =
            read_analysis_limited(&crate::abaqus::write(&mixed_read.model, Preset::MKgS).unwrap(), Preset::MKgS, Limits::default()).unwrap();
        assert_totals_match(
            "same node and element set, write/read step 2",
            &BTreeMap::from([((false, 4), [10.0, 0.0, 0.0])]),
            &imported_load_totals(&mixed_reread.model, 1),
        );

        let conflict = format!("{INVARIANT_MODEL}*NSET,NSET=Tip\n4\n*STEP, NAME=Bad\n*STATIC\n*BOUNDARY\nTip, 1, 1, 1\n4, 1, 2\n*END STEP\n");
        let rejected = read_analysis_limited(&conflict, Preset::MKgS, Limits::default()).unwrap();
        assert!(rejected.model.steps().is_empty());
        assert!(rejected.warnings.iter().any(|warning| warning.contains("overlap")), "{:?}", rejected.warnings);
    }

    #[test]
    fn transformed_nodes_keep_local_history_unread_without_changing_model_data() {
        let transform = "*TRANSFORM,NSET=Tip,TYPE=R\n0,1,0,-1,0,0\n";
        let history = |body: &str| format!("*STEP, NAME=S1\n*STATIC\n{body}*END STEP\n");
        let assert_unread = |deck: &str, marker: &str| {
            let base = read_in(deck, Preset::MKgS).unwrap();
            let got = read_analysis_limited(deck, Preset::MKgS, Limits::default()).unwrap();
            assert_eq!(got.model.nodes(), base.model.nodes());
            assert_eq!(got.model.elements(), base.model.elements());
            assert_eq!(got.model.materials(), base.model.materials());
            assert_eq!(got.model.sections(), base.model.sections());
            assert_eq!(got.model.sets(), base.model.sets());
            assert_eq!(got.transforms, base.transforms);
            assert!(got.model.steps().is_empty(), "{marker}: {:?}", got.model.steps());
            assert!(got.unknown.iter().any(|unknown| unknown.keyword == "TRANSFORM"));
            assert!(got.warnings.iter().any(|warning| warning.contains(marker)), "{marker}: {:?}", got.warnings);
        };

        assert_unread(&format!("{MODEL}{transform}{}", history("*CLOAD\nTip, 1, 10\n")), "local coordinate system");
        assert_unread(&format!("{MODEL}{transform}{}", history("*CLOAD\n2, 1, 10\n")), "local coordinate system");
        assert_unread(&format!("{MODEL}{transform}*BOUNDARY\nTip, ENCASTRE\n{}", history("")), "local coordinate system");
        assert_unread(&format!("{MODEL}*TRANSFORM,NSET=Missing,TYPE=R\n0,1,0,-1,0,0\n{}", history("*CLOAD\nRoot, 1, 10\n")), "cannot be resolved");

        let unaffected = format!("{MODEL}{transform}{}", history("*CLOAD\nRoot, 1, 10\n"));
        let got = read_analysis_limited(&unaffected, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 1);
        assert!(got.unknown.iter().any(|unknown| unknown.keyword == "TRANSFORM"));
        assert!(got.warnings.is_empty(), "{:?}", got.warnings);
    }

    #[test]
    fn frequency_boundary_overlap_uses_zero_only_for_conflicts_and_keeps_values() {
        let model = format!("{MODEL}*NSET,NSET=Other\n2\n");
        let before_and_after = format!(
            "{model}*STEP, NAME=Before\n*STATIC\n*BOUNDARY\nTip, 1, 1, 0.5\n*END STEP\n*STEP, NAME=Modes\n*FREQUENCY\n2\n*BOUNDARY\nOther, 1, 1, 0\n*END STEP\n*STEP, NAME=After\n*STATIC\n*END STEP\n"
        );
        let got = read_analysis_limited(&before_and_after, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 3);
        assert_eq!(got.model.steps()[0].boundary_conditions, got.model.steps()[2].boundary_conditions);
        let frequency_values: Vec<_> = got.model.steps()[1]
            .boundary_conditions
            .iter()
            .filter_map(|name| got.model.boundary_condition(name))
            .filter_map(|boundary| match boundary {
                BoundaryCondition::Displacement { values, .. } => Some(*values),
                _ => None,
            })
            .collect();
        assert!(frequency_values.iter().any(|values| values[0] == Some(0.5)), "{frequency_values:?}");
        assert!(frequency_values.iter().any(|values| values[0] == Some(0.0)), "{frequency_values:?}");

        let frequency_own = format!(
            "{model}*STEP, NAME=Modes\n*FREQUENCY\n2\n*BOUNDARY\nTip, 1, 1, 0.5\nOther, 1, 1, 0\n*END STEP\n*STEP, NAME=After\n*STATIC\n*END STEP\n"
        );
        let got = read_analysis_limited(&frequency_own, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 2);
        assert_eq!(got.model.steps()[0].boundary_conditions.len(), 2);
        assert!(
            got.model.steps()[0]
                .boundary_conditions
                .iter()
                .filter_map(|name| got.model.boundary_condition(name))
                .any(|boundary| matches!(boundary, BoundaryCondition::Displacement { values, .. } if values[0] == Some(0.5)))
        );
        assert!(
            got.model.steps()[0]
                .boundary_conditions
                .iter()
                .filter_map(|name| got.model.boundary_condition(name))
                .any(|boundary| matches!(boundary, BoundaryCondition::Displacement { values, .. } if values[0] == Some(0.0)))
        );
        assert!(got.model.steps()[1].boundary_conditions.is_empty());

        let static_control = format!("{model}*STEP, NAME=Bad\n*STATIC\n*BOUNDARY\nTip, 1, 1, 0.5\nOther, 1, 1, 0\n*END STEP\n");
        let rejected = read_analysis_limited(&static_control, Preset::MKgS, Limits::default()).unwrap();
        assert!(rejected.model.steps().is_empty());
        assert!(rejected.warnings.iter().any(|warning| warning.contains("overlap")), "{:?}", rejected.warnings);
    }

    #[test]
    fn pinned_and_direct_boundary_formats_have_distinct_history_state() {
        let assert_type_change_unread = |deck: &str, expected_steps, marker: &str| {
            let got = read_analysis_limited(deck, Preset::MKgS, Limits::default()).unwrap();
            assert_eq!(got.model.steps().len(), expected_steps, "{:?}", got.warnings);
            assert!(got.warnings.iter().any(|warning| warning.contains(marker)), "{:?}", got.warnings);
        };

        assert_type_change_unread(
            &format!(
                "{MODEL}*STEP, NAME=S1\n*STATIC\n*BOUNDARY\nTip, PINNED\n*END STEP\n*STEP, NAME=S2\n*STATIC\n*BOUNDARY\nTip, 1, 1, 0.5\n*END STEP\n"
            ),
            1,
            "type-format and direct-format",
        );
        assert_type_change_unread(
            &format!(
                "{MODEL}*STEP, NAME=S1\n*STATIC\n*BOUNDARY\nTip, 1, 1, 0.5\n*END STEP\n*STEP, NAME=S2\n*STATIC\n*BOUNDARY\nTip, PINNED\n*END STEP\n"
            ),
            1,
            "type-format and direct-format",
        );
        assert_type_change_unread(
            &format!("{MODEL}*BOUNDARY\nTip, PINNED\n*STEP, NAME=S1\n*STATIC\n*BOUNDARY\nTip, 1, 1, 0.5\n*END STEP\n"),
            0,
            "type-format and direct-format",
        );
        assert_type_change_unread(
            &format!(
                "{MODEL}*STEP, NAME=S1\n*STATIC\n*BOUNDARY\nTip, PINNED\n*END STEP\n*STEP, NAME=Modes\n*FREQUENCY\n2\n*BOUNDARY\nTip, 1, 1, 0.5\n*END STEP\n"
            ),
            1,
            "type-format and direct-format",
        );
        assert_type_change_unread(
            &format!(
                "{MODEL}*STEP, NAME=S1\n*STATIC\n*BOUNDARY\nTip, PINNED\n*END STEP\n*STEP, NAME=Buckle\n*BUCKLE\n2\n*BOUNDARY\nTip, 1, 1, 0.5\n*END STEP\n"
            ),
            1,
            "type-format and direct-format",
        );
        assert_type_change_unread(
            &format!(
                "{MODEL}*STEP, NAME=S1\n*STATIC\n*BOUNDARY\nTip, PINNED\n*END STEP\n*STEP, NAME=S2\n*STATIC\n*BOUNDARY\nTip, 1, 1, 0\n*END STEP\n"
            ),
            1,
            "type-format and direct-format",
        );
        assert_type_change_unread(
            &format!(
                "{MODEL}*STEP, NAME=S1\n*STATIC\n*BOUNDARY\nTip, 1, 1, 0\n*END STEP\n*STEP, NAME=S2\n*STATIC\n*BOUNDARY\nTip, PINNED\n*END STEP\n"
            ),
            1,
            "type-format and direct-format",
        );

        let op_new = format!(
            "{MODEL}*STEP, NAME=S1\n*STATIC\n*BOUNDARY\nTip, PINNED\n*END STEP\n*STEP, NAME=S2\n*STATIC\n*BOUNDARY, OP=NEW\nTip, 1, 1, 0.5\n*END STEP\n"
        );
        let got = read_analysis_limited(&op_new, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 2, "{:?}", got.warnings);

        let same_type = format!(
            "{MODEL}*STEP, NAME=S1\n*STATIC\n*BOUNDARY\nTip, PINNED\n*END STEP\n*STEP, NAME=S2\n*STATIC\n*BOUNDARY\nTip, PINNED\n*END STEP\n"
        );
        let got = read_analysis_limited(&same_type, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 2, "{:?}", got.warnings);

        let type_cases: [(&str, &[u8], Option<u8>); 5] = [
            ("PINNED", &[1, 2, 3], Some(5)),
            ("ENCASTRE", &[1, 2, 3, 4, 5, 6], None),
            ("XSYMM", &[1, 5, 6], Some(2)),
            ("YSYMM", &[2, 4, 6], Some(1)),
            ("ZSYMM", &[3, 4, 5], Some(1)),
        ];
        for (label, controlled, uncontrolled) in type_cases {
            let controlled = controlled[0];
            let deck = format!(
                "{MODEL}*STEP, NAME=S1\n*STATIC\n*BOUNDARY\nTip, {label}\n*END STEP\n*STEP, NAME=S2\n*STATIC\n*BOUNDARY\nTip, {controlled}, {controlled}, 0.1\n*END STEP\n"
            );
            assert_type_change_unread(&deck, 1, "require OP=NEW");
            if let Some(uncontrolled) = uncontrolled {
                let deck = format!(
                    "{MODEL}*STEP, NAME=S1\n*STATIC\n*BOUNDARY\nTip, {label}\n*END STEP\n*STEP, NAME=S2\n*STATIC\n*BOUNDARY\nTip, {uncontrolled}, {uncontrolled}, 0.1\n*END STEP\n"
                );
                let got = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
                assert_eq!(got.model.steps().len(), 2, "{label}: {:?}", got.warnings);
            }
        }
    }

    #[test]
    fn shared_node_zero_boundary_rows_roundtrip_for_type_and_pinned_forms() {
        let model = format!("{MULTI_MODEL}*NSET,NSET=End\n4\n");
        for (label, type_row, written_type, written_direct) in
            [("ENCASTRE", "Tip, ENCASTRE", "Tip, ENCASTRE", "End, 1, 1, 0.0"), ("PINNED", "Tip, PINNED", "Tip, 1, 1, 0.0", "End, 1, 1, 0.0")]
        {
            let deck = format!(
                "{model}*BOUNDARY\n{type_row}\n*STEP, NAME=S1\n*STATIC\n*END STEP\n*STEP, NAME=S2\n*STATIC\n*BOUNDARY\nEnd, 1, 1, 0\n*END STEP\n"
            );
            let first = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
            assert_eq!(first.model.steps().len(), 2, "{label}: {:?}", first.warnings);
            let written = crate::abaqus::write(&first.model, Preset::MKgS).unwrap();
            assert!(written.contains(written_type), "{label}: missing {written_type:?} in {written}");
            assert!(written.contains(written_direct), "{label}: missing {written_direct:?} in {written}");
            assert_analysis_write_roundtrip(label, &first);
        }
    }

    #[test]
    fn generated_history_sets_are_not_valid_input_targets() {
        let deck = format!(
            "{MULTI_MODEL}*STEP, NAME=First\n*STATIC\n*CLOAD\n4, 1, 5\n*END STEP\n*STEP, NAME=GeneratedTarget\n*STATIC\n*CLOAD\nInpAnalysis1, 1, 3\n*END STEP\n"
        );
        let got = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 1);
        assert!(got.model.sets().keys().any(|name| name.eq_ignore_ascii_case("InpAnalysis1")));
        assert!(got.warnings.iter().any(|warning| warning.contains("missing node set InpAnalysis1")), "{:?}", got.warnings);
    }

    #[test]
    fn model_boundary_after_first_step_stops_history_without_mutating_the_model() {
        let deck = format!("{MODEL}*STEP, NAME=First\n*STATIC\n*END STEP\n*BOUNDARY\nTip, ENCASTRE\n*STEP, NAME=Later\n*STATIC\n*END STEP\n");
        let got = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
        assert_eq!(got.model.steps().len(), 1);
        assert!(got.model.boundary_conditions().is_empty());
        assert!(got.warnings.iter().any(|warning| warning.contains("after a *STEP")), "{:?}", got.warnings);
    }

    #[test]
    fn repeated_transforms_keep_only_names_and_resolve_each_input_set_once() {
        let nodes = 2_048;
        let transforms = 6_000;
        let mut deck = String::from("*NODE\n");
        for node in 1..=nodes {
            deck.push_str(&format!("{node},{node},0,0\n"));
        }
        deck.push_str("*ELEMENT,TYPE=T3D2,ELSET=Bar\n");
        for element in 1..nodes {
            deck.push_str(&format!("{element},{element},{}\n", element + 1));
        }
        deck.push_str(&format!("*NSET,NSET=All,GENERATE\n1,{nodes},1\n*MATERIAL,NAME=Steel\n*ELASTIC\n210000,0.3\n*DENSITY\n1\n*SOLID SECTION,ELSET=Bar,MATERIAL=Steel\n1\n"));
        for _ in 0..transforms {
            deck.push_str("*TRANSFORM,NSET=All,TYPE=R\n0,1,0,-1,0,0\n");
        }
        let base = read_in(&deck, Preset::MKgS).unwrap();
        assert_eq!(base.transforms.len(), transforms);
        assert!(base.transforms.iter().all(|transform| transform.target.as_deref() == Some("ALL")));
        let with_history = format!("{deck}*STEP, NAME=Covered\n*STATIC\n*CLOAD\nAll, 1, 1\n*END STEP\n");
        let got = read_analysis_limited(&with_history, Preset::MKgS, Limits::default()).unwrap();
        assert!(got.model.steps().is_empty());
        assert!(got.warnings.iter().any(|warning| warning.contains("local coordinate system")), "{:?}", got.warnings);
    }

    fn next_mutation_value(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    fn toggle_pinned_row(lines: &mut [String]) -> bool {
        let candidate = lines.iter().position(|line| {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('*') {
                return false;
            }
            let fields: Vec<_> = trimmed.split(',').map(|field| field.trim().to_owned()).collect();
            if fields.len() == 2 && fields[1].eq_ignore_ascii_case("PINNED") {
                return true;
            }
            if !(3..=4).contains(&fields.len()) || fields[1] != "1" || fields[2] != "3" {
                return false;
            }
            fields.get(3).is_none_or(|value| value.is_empty() || value == "0" || value == "0.0")
        });
        let Some(index) = candidate else { return false };
        let fields: Vec<_> = lines[index].split(',').map(|field| field.trim().to_owned()).collect();
        if fields.get(1).is_some_and(|value| value.eq_ignore_ascii_case("PINNED")) {
            lines[index] = format!("{}, 1, 3", fields[0]);
        } else {
            lines[index] = format!("{}, PINNED", fields[0]);
        }
        true
    }

    fn mutate_deck(source: &str, state: &mut u64) -> String {
        let mut lines: Vec<String> = source.lines().map(str::to_owned).collect();
        let history_start = lines.iter().position(|line| line.trim_start().to_ascii_uppercase().starts_with("*STEP")).unwrap_or(lines.len());
        let history_len = lines.len().saturating_sub(history_start);
        match next_mutation_value(state) % 9 {
            0 if history_len > 1 => {
                let index = history_start + (next_mutation_value(state) as usize) % history_len;
                lines.remove(index);
            }
            1 if history_len > 0 => {
                let index = history_start + (next_mutation_value(state) as usize) % history_len;
                lines.insert(index, lines[index].clone());
            }
            2 if history_len > 1 => {
                let first = history_start + (next_mutation_value(state) as usize) % history_len;
                let second = history_start + (next_mutation_value(state) as usize) % history_len;
                lines.swap(first, second);
            }
            3 => {
                let keep = if history_len == 0 { lines.len() } else { history_start + (next_mutation_value(state) as usize) % (history_len + 1) };
                lines.truncate(keep);
            }
            4 if history_len > 0 => {
                let index = history_start + (next_mutation_value(state) as usize) % history_len;
                let mut fields: Vec<String> = lines[index].split(',').map(str::to_owned).collect();
                if fields.len() == 1 {
                    fields.push("BAD".into());
                }
                let field = (next_mutation_value(state) as usize) % fields.len();
                fields[field] = match next_mutation_value(state) % 6 {
                    0 => "0".into(),
                    1 => "1".into(),
                    2 => "BAD".into(),
                    3 => "Tip".into(),
                    4 => "Other".into(),
                    _ => "999999999".into(),
                };
                lines[index] = fields.join(",");
            }
            5 => {
                let candidate = lines.iter().position(|line| {
                    let upper = line.to_ascii_uppercase();
                    upper.starts_with("*STEP") && upper.contains("NAME=")
                });
                if let Some(index) = candidate {
                    let upper = lines[index].to_ascii_uppercase();
                    if let Some(name_start) = upper.find("NAME=") {
                        let start = name_start + "NAME=".len();
                        let end = upper[start..].find(',').map_or(lines[index].len(), |offset| start + offset);
                        let mut toggled = false;
                        let name: String = lines[index][start..end]
                            .chars()
                            .map(|character| {
                                if !toggled && character.is_ascii_alphabetic() {
                                    toggled = true;
                                    if character.is_ascii_lowercase() { character.to_ascii_uppercase() } else { character.to_ascii_lowercase() }
                                } else {
                                    character
                                }
                            })
                            .collect();
                        lines[index].replace_range(start..end, &name);
                    }
                }
            }
            6 => {
                let candidate = lines.iter().position(|line| {
                    let upper = line.to_ascii_uppercase();
                    upper.starts_with("*STEP") && upper.contains("NAME=")
                });
                if let Some(index) = candidate {
                    let fields: Vec<_> =
                        lines[index].split(',').filter(|field| !field.trim_start().to_ascii_uppercase().starts_with("NAME=")).collect();
                    lines[index] = fields.join(",");
                }
            }
            7 => {
                let candidate =
                    lines.iter().position(|line| matches!(line.trim().to_ascii_uppercase().as_str(), "*STATIC" | "*FREQUENCY" | "*BUCKLE"));
                if let Some(index) = candidate {
                    lines[index] = match lines[index].trim().to_ascii_uppercase().as_str() {
                        "*STATIC" => "*FREQUENCY".into(),
                        "*FREQUENCY" => "*BUCKLE".into(),
                        _ => "*STATIC".into(),
                    };
                }
            }
            8 if toggle_pinned_row(&mut lines) => {}
            8 => {
                let candidate = lines.iter().position(|line| {
                    let upper = line.trim_start().to_ascii_uppercase();
                    upper.starts_with("*BOUNDARY") || upper.starts_with("*CLOAD")
                });
                if let Some(index) = candidate {
                    let fields: Vec<_> = lines[index].split(',').map(str::to_owned).collect();
                    if fields.iter().any(|field| field.trim().eq_ignore_ascii_case("OP=NEW")) {
                        lines[index] = fields.into_iter().filter(|field| !field.trim().eq_ignore_ascii_case("OP=NEW")).collect::<Vec<_>>().join(",");
                    } else {
                        lines[index].push_str(", OP=NEW");
                    }
                }
            }
            _ => {}
        }
        lines.join("\n")
    }

    #[test]
    fn fixed_seed_mutations_never_panic_and_keep_model_data_when_both_reads_succeed() {
        let seeds = [
            format!(
                "{MODEL}*TRANSFORM,NSET=Tip,TYPE=R\n0,1,0,-1,0,0\n*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=S1\n*STATIC\n*CLOAD\n2, 1, 5\n*END STEP\n"
            ),
            format!(
                "{MULTI_MODEL}*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=S1\n*STATIC\n*BOUNDARY\nTip, 1, 1, 0.2\n*CLOAD\n4, 1, 3\n*END STEP\n*STEP, NAME=Modes\n*FREQUENCY\n2\n*BOUNDARY\nMid, 1, 1, 0\n*END STEP\n"
            ),
            format!(
                "{INVARIANT_MODEL}*NSET,NSET=LoadSet\n2,3\n*TRANSFORM,NSET=LoadSet,TYPE=R\n0,1,0,-1,0,0\n*STEP, NAME=S1\n*STATIC\n*BOUNDARY\nRoot, ENCASTRE\n*CLOAD\nLoadSet, 1, 7\n*END STEP\n"
            ),
            format!("{MODEL}*BOUNDARY\nRoot, ENCASTRE\n*STEP\n*STATIC\n*END STEP\n*STEP, NAME=importedstep2\n*STATIC\n*END STEP\n"),
            format!(
                "{MULTI_MODEL}*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=Base\n*STATIC\n*END STEP\n*STEP, NAME=Buckle\n*BUCKLE\n2\n*END STEP\n*STEP, NAME=After\n*STATIC\n*END STEP\n"
            ),
            format!("{MULTI_MODEL}*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=Modes\n*FREQUENCY\n2\n*END STEP\n*STEP, NAME=After\n*STATIC\n*END STEP\n"),
            format!(
                "{INVARIANT_BEAM_MODEL}*NSET,NSET=Tip\n4\n*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=Both\n*STATIC\n*CLOAD\nTip, 1, 5\nTip, 5, 2\n*END STEP\n"
            ),
            format!("{MULTI_MODEL}*BOUNDARY\nRoot, ENCASTRE\n*STEP, NAME=Pinned\n*STATIC\n*BOUNDARY\nTip, PINNED\n*CLOAD\n4, 1, 3\n*END STEP\n"),
        ];
        const MUTATIONS_PER_SEED: usize = 500;
        for (seed_index, seed) in seeds.iter().enumerate() {
            let mut state = 0x1400_0000_0000_0001_u64 + seed_index as u64;
            for mutation in 0..MUTATIONS_PER_SEED {
                let deck = mutate_deck(seed, &mut state);
                let pure = catch_unwind(AssertUnwindSafe(|| read_in(&deck, Preset::MKgS)));
                assert!(pure.is_ok(), "pure reader panicked for seed {seed_index}, mutation {mutation}");
                let pure = pure.unwrap();
                let analysis = catch_unwind(AssertUnwindSafe(|| read_analysis_limited(&deck, Preset::MKgS, Limits::default())));
                assert!(analysis.is_ok(), "analysis reader panicked for seed {seed_index}, mutation {mutation}");
                if let Ok(base) = pure {
                    let analysis = analysis.unwrap();
                    assert!(
                        analysis.is_ok(),
                        "analysis reader failed after pure read for seed {seed_index}, mutation {mutation}: {:?}",
                        analysis.err()
                    );
                    let got = analysis.unwrap();
                    assert_eq!(got.model.nodes(), base.model.nodes(), "nodes differ for seed {seed_index}, mutation {mutation}");
                    assert_eq!(got.model.elements(), base.model.elements(), "elements differ for seed {seed_index}, mutation {mutation}");
                    assert_eq!(got.model.materials(), base.model.materials(), "materials differ for seed {seed_index}, mutation {mutation}");
                    assert_eq!(got.model.sections(), base.model.sections(), "sections differ for seed {seed_index}, mutation {mutation}");
                    for (name, set) in base.model.sets() {
                        assert_eq!(got.model.set(name), Some(set), "input set {name} differs for seed {seed_index}, mutation {mutation}");
                    }
                    assert_analysis_write_roundtrip(&format!("seed {seed_index}, mutation {mutation}"), &got);
                }
            }
        }
    }

    #[test]
    fn unsupported_history_stays_reported_without_creating_a_partial_step() {
        let deck = format!("{MODEL}*STEP,NAME=Later\n*STATIC\n1,1\n*DLOAD\nBar, BX, 1\n*END STEP\n");
        let got = read_analysis_limited(&deck, Preset::MKgS, Limits::default()).unwrap();
        assert!(got.model.steps().is_empty());
        assert!(got.unknown.iter().any(|unknown| unknown.keyword == "DLOAD"));
    }

    #[test]
    fn source_fixture_has_the_expected_model_types() {
        let got = read_in(MODEL, Preset::MKgS).unwrap();
        assert_eq!(got.model.element(1).map(|e| e.kind), Some(ElementKind::Rod2));
        assert_eq!(got.model.material("Steel").map(|m| m.elastic), Some(Some(Elastic::Isotropic { young: 210000., poisson: 0.3 })));
        assert!(got.model.sections().values().any(|section| matches!(section.kind, SectionKind::Rod { .. })));
    }
}
