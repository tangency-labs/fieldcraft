//! Deterministic part-instance expansion before flat model construction.
//! Designed from the public Abaqus Assembly Definition and INSTANCE references:
//! https://abaqus.uclouvain.be/English/SIMACAEMODRefMap/simamod-c-partassy.htm
//! https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-instance.htm
//! Translations precede rotations; identifiers belong to their part instance.
use super::{Block, Limits, MAX_INPUT_BYTES, MAX_RECORDS, canonical, check_work, error, label, lex, number, row_fields, types::Type};
use crate::{Error, Result};
use fieldcraft_model::units::{Quantity, Unit};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InstanceMap {
    pub name: String,
    pub part: String,
    pub nodes: BTreeMap<u32, u32>,
    pub elements: BTreeMap<u32, u32>,
    /// Original local names to FieldCraft names (sets and orientations).
    pub names: BTreeMap<String, String>,
}

pub(super) struct Expanded {
    pub text: String,
    lines: Vec<usize>,
    pub instances: Vec<InstanceMap>,
    pub warnings: Vec<String>,
}
impl Expanded {
    fn push(&mut self, text: &str, line: usize) -> Result<()> {
        if self.text.len().saturating_add(text.len()).saturating_add(1) > MAX_INPUT_BYTES || self.lines.len() >= MAX_RECORDS {
            return Err(error(line, "expanded assembly stream limit reached"));
        }
        if text.len() > super::super::MAX_LINE {
            return Err(error(line, "expanded assembly record exceeds 256 bytes"));
        }
        self.text.push_str(text);
        self.text.push('\n');
        self.lines.push(line);
        Ok(())
    }
    fn header(&mut self, keyword: &str, params: &BTreeMap<String, Option<String>>, line: usize) -> Result<()> {
        let mut text = format!("*{keyword}");
        for (key, value) in params {
            let token = value.as_ref().map_or_else(|| key.clone(), |v| format!("{key}={v}"));
            if text.len() + token.len() + 2 > super::super::MAX_LINE - 1 {
                text.push(',');
                self.push(&text, line)?;
                text.clear();
            } else {
                text.push_str(", ");
            }
            text.push_str(&token);
        }
        self.push(&text, line)
    }
    pub fn original_line(&self, line: usize) -> usize {
        self.lines.get(line.saturating_sub(1)).copied().unwrap_or(line)
    }
    pub fn source_error(&self, e: Error) -> Error {
        match e {
            Error::Parse { line, message } => error(self.original_line(line), message),
            other => other,
        }
    }
}

#[derive(Clone, Default)]
struct Transform {
    shift: [f64; 3],
    origin: [f64; 3],
    axis: [f64; 3],
    angle: f64,
}
impl Transform {
    fn from_block(b: &Block<'_>) -> Result<Self> {
        if b.rows.len() > 2 {
            return Err(error(b.line, "INSTANCE accepts a translation and optional axis rotation"));
        }
        let mut out = Self::default();
        if let Some(row) = b.rows.first() {
            let f = row_fields(row.text, row.line)?;
            if f.len() != 3 {
                return Err(error(row.line, "INSTANCE translation needs three coordinates"));
            }
            for (i, x) in out.shift.iter_mut().enumerate() {
                *x = number(f.get(i).map(String::as_str).unwrap_or_default(), row.line)?;
            }
        }
        if let Some(row) = b.rows.get(1) {
            let f = row_fields(row.text, row.line)?;
            if f.len() != 7 {
                return Err(error(row.line, "INSTANCE rotation needs two axis points and an angle in degrees"));
            }
            let v: Vec<f64> = f.iter().map(|s| number(s, row.line)).collect::<Result<_>>()?;
            for i in 0..3 {
                if let (Some(o), Some(a), Some(x), Some(y)) = (out.origin.get_mut(i), out.axis.get_mut(i), v.get(i), v.get(i + 3)) {
                    *o = *x;
                    *a = y - x;
                }
            }
            let norm = out.axis.iter().fold(0.0_f64, |n, x| n.hypot(*x));
            if !norm.is_finite() || norm == 0.0 {
                return Err(error(row.line, "INSTANCE rotation axis must have finite nonzero length"));
            }
            for a in &mut out.axis {
                *a /= norm;
            }
            let unit = Unit::parse("deg").map_err(|e| error(row.line, e.to_string()))?;
            out.angle = Quantity::Angle.to_si(v.get(6).copied().unwrap_or_default(), &unit).map_err(|e| error(row.line, e.to_string()))?;
        }
        Ok(out)
    }
    fn vector(&self, v: [f64; 3]) -> [f64; 3] {
        if self.angle == 0.0 {
            return v;
        }
        let [x, y, z] = v;
        let [a, b, c] = self.axis;
        let dot = a * x + b * y + c * z;
        let cross = [b * z - c * y, c * x - a * z, a * y - b * x];
        let (sin, cos) = self.angle.sin_cos();
        std::array::from_fn(|i| v[i] * cos + cross[i] * sin + self.axis[i] * dot * (1.0 - cos))
    }
    fn point(&self, v: [f64; 3], line: usize) -> Result<[f64; 3]> {
        let moved = std::array::from_fn(|i| v[i] + self.shift[i] - self.origin[i]);
        let rotated = self.vector(moved);
        let out = std::array::from_fn(|i| rotated[i] + self.origin[i]);
        if out.iter().any(|x| !x.is_finite()) {
            return Err(error(line, "instance transform produced non-finite coordinates"));
        }
        Ok(out)
    }
}

struct Instance<'a> {
    name: String,
    part: String,
    line: usize,
    transform: Transform,
    blocks: Vec<&'a Block<'a>>,
}
fn empty(b: &Block<'_>) -> Result<()> {
    if !b.rows.is_empty() {
        return Err(error(b.line, format!("*{} has no data lines", b.keyword)));
    }
    Ok(())
}
fn id(s: &str, line: usize) -> Result<String> {
    let name = lex::unquote(s, line)?;
    if name.is_empty() || name.contains('.') {
        return Err(error(line, "assembly names must be nonempty and cannot contain a period"));
    }
    Ok(name.to_ascii_uppercase())
}

pub(super) fn expand(text: &str, limits: Limits) -> Result<Option<Expanded>> {
    let blocks = lex::blocks(text, MAX_RECORDS)?;
    if !blocks.iter().any(|b| matches!(b.keyword.as_str(), "PART" | "ENDPART" | "ASSEMBLY" | "ENDASSEMBLY" | "INSTANCE" | "ENDINSTANCE")) {
        return Ok(None);
    }
    let mut parts: BTreeMap<String, Vec<&Block<'_>>> = BTreeMap::new();
    let mut part: Option<String> = None;
    let mut assembly = None;
    let mut assembly_open = false;
    let mut instance: Option<Instance<'_>> = None;
    let mut instances = Vec::new();
    let mut global = Vec::new();
    for b in &blocks {
        check_work(b.line)?;
        match b.keyword.as_str() {
            "PART" => {
                b.only(&["NAME"])?;
                empty(b)?;
                if part.is_some() || assembly_open || instance.is_some() {
                    return Err(error(b.line, "nested PART or PART inside ASSEMBLY"));
                }
                let name = id(b.required("NAME")?, b.line)?;
                if parts.insert(name.clone(), Vec::new()).is_some() {
                    return Err(error(b.line, "duplicate PART"));
                }
                part = Some(name);
            }
            "ENDPART" => {
                b.only(&[])?;
                empty(b)?;
                if part.take().is_none() {
                    return Err(error(b.line, "END PART without PART"));
                }
            }
            "ASSEMBLY" => {
                b.only(&["NAME"])?;
                empty(b)?;
                if part.is_some() || assembly.is_some() || instance.is_some() {
                    return Err(error(b.line, "only one non-nested ASSEMBLY is supported"));
                }
                assembly = Some(id(b.required("NAME")?, b.line)?);
                assembly_open = true;
            }
            "ENDASSEMBLY" => {
                b.only(&[])?;
                empty(b)?;
                if !assembly_open || instance.is_some() {
                    return Err(error(b.line, "END ASSEMBLY without ASSEMBLY or before END INSTANCE"));
                }
                assembly_open = false;
            }
            "INSTANCE" => {
                b.only(&["NAME", "PART"])?;
                if !assembly_open || instance.is_some() || part.is_some() {
                    return Err(error(b.line, "INSTANCE must be directly inside ASSEMBLY"));
                }
                if instances.len() >= 10_000 {
                    return Err(error(b.line, "instance count limit reached"));
                }
                instance = Some(Instance {
                    name: id(b.required("NAME")?, b.line)?,
                    part: id(b.required("PART")?, b.line)?,
                    line: b.line,
                    transform: Transform::from_block(b)?,
                    blocks: Vec::new(),
                });
            }
            "ENDINSTANCE" => {
                b.only(&[])?;
                empty(b)?;
                instances.push(instance.take().ok_or_else(|| error(b.line, "END INSTANCE without INSTANCE"))?);
            }
            "STEP" if part.is_some() || assembly_open => return Err(error(b.line, "STEP must be outside PART and ASSEMBLY")),
            _ => {
                if let Some(i) = instance.as_mut() {
                    i.blocks.push(b);
                } else if let Some(p) = &part {
                    if let Some(rows) = parts.get_mut(p) {
                        rows.push(b);
                    }
                } else {
                    global.push((b, !assembly_open));
                }
            }
        }
    }
    if part.is_some() || assembly_open || instance.is_some() {
        return Err(error(blocks.last().map_or(1, |b| b.line), "unclosed PART, ASSEMBLY or INSTANCE"));
    }
    if !parts.is_empty() && instances.is_empty() {
        return Err(error(1, "parts have no instances; no assembly mesh to import"));
    }
    let mut next_node = 0;
    let mut next_element = 0;
    let global_blocks: Vec<&Block<'_>> = global.iter().map(|(b, _)| *b).collect();
    let (global_nodes, global_elements) = labels(&global_blocks)?;
    if let Some(v) = global_nodes.last() {
        next_node = *v;
    }
    if let Some(v) = global_elements.last() {
        next_element = *v;
    }
    let mut count_nodes = global_nodes.len();
    let mut count_elements = global_elements.len();
    let mut maps = BTreeMap::new();
    let mut used_names = BTreeSet::new();
    for b in &blocks {
        for key in ["NSET", "ELSET", "NAME"] {
            if let Some(v) = b.value(key) {
                used_names.insert(lex::unquote(v, b.line)?.to_ascii_uppercase());
            }
        }
    }
    let mut sequence = 0usize;
    let mut expanded_blocks = global.len();
    for i in &mut instances {
        let part = parts.get(&i.part).ok_or_else(|| error(i.line, format!("unknown part {}", i.part)))?;
        expanded_blocks = expanded_blocks.saturating_add(part.len()).saturating_add(i.blocks.len());
        if expanded_blocks > MAX_RECORDS {
            return Err(error(i.line, "expanded instance block limit reached"));
        }
        let mut local = part.clone();
        local.append(&mut i.blocks);
        i.blocks = local;
        let (nodes, elements) = labels(&i.blocks)?;
        count_nodes = count_nodes.saturating_add(nodes.len());
        count_elements = count_elements.saturating_add(elements.len());
        if count_nodes > limits.nodes || count_elements > limits.elements {
            return Err(error(i.line, "instance expansion exceeds node or element limit"));
        }
        let mut map =
            InstanceMap { name: i.name.clone(), part: i.part.clone(), nodes: BTreeMap::new(), elements: BTreeMap::new(), names: BTreeMap::new() };
        for (old, next, dest) in [(nodes, &mut next_node, &mut map.nodes), (elements, &mut next_element, &mut map.elements)] {
            for label in old {
                *next =
                    next.checked_add(1).filter(|v| *v <= fieldcraft_model::MAX_LABEL).ok_or_else(|| error(i.line, "expanded label limit reached"))?;
                dest.insert(label, *next);
            }
        }
        for b in &i.blocks {
            for key in ["NSET", "ELSET", "NAME"] {
                if let Some(value) = b.value(key) {
                    if key == "NAME" && !matches!(b.keyword.as_str(), "SURFACE" | "ORIENTATION") {
                        continue;
                    }
                    let key = id(value, b.line)?;
                    if let std::collections::btree_map::Entry::Vacant(entry) = map.names.entry(key) {
                        let mut available = None;
                        while sequence < MAX_RECORDS.saturating_mul(4) {
                            sequence += 1;
                            let name = format!("InpLocal{sequence}");
                            if used_names.insert(name.to_ascii_uppercase()) {
                                available = Some(name);
                                break;
                            }
                        }
                        entry.insert(available.ok_or_else(|| error(b.line, "instance name allocation limit reached"))?);
                    }
                }
            }
        }
        if maps.insert(i.name.clone(), map).is_some() {
            return Err(error(i.line, format!("duplicate INSTANCE {}", i.name)));
        }
    }
    let mut warnings = Vec::new();
    let instanced: BTreeSet<&str> = instances.iter().map(|i| i.part.as_str()).collect();
    let unused: Vec<&String> = parts.keys().filter(|part| !instanced.contains(part.as_str())).collect();
    if !unused.is_empty() {
        warnings.push(format!(
            "{} uninstanced parts retained only in original source (first 100): {:?}",
            unused.len(),
            unused.iter().take(100).collect::<Vec<_>>()
        ));
    }
    let mut orthotropic = BTreeSet::new();
    let mut material = None;
    for b in &blocks {
        if b.keyword == "MATERIAL" {
            material = Some(lex::unquote(b.required("NAME")?, b.line)?.to_ascii_uppercase());
        } else if b.keyword == "ELASTIC" && b.value("TYPE").is_some_and(|v| canonical(v) == "ENGINEERINGCONSTANTS") {
            if let Some(name) = &material {
                orthotropic.insert(name.clone());
            }
        } else if !super::material_behavior(&b.keyword) {
            material = None;
        }
    }
    let mut out = Expanded { text: String::new(), lines: Vec::new(), instances: maps.values().cloned().collect(), warnings };
    let assembly = assembly.as_deref().unwrap_or_default();
    let mut members = 0usize;
    for i in &instances {
        for b in &i.blocks {
            if i.transform.angle != 0.0 && b.keyword == "SOLIDSECTION" && !b.params.contains_key("ORIENTATION") {
                let material = lex::unquote(b.required("MATERIAL")?, b.line)?.to_ascii_uppercase();
                if orthotropic.contains(&material) {
                    return Err(error(b.line, "rotated orthotropic solid requires explicit part ORIENTATION; implicit axes are not inferred"));
                }
            }
            emit(&mut out, b, Some(i), &maps, (assembly, false), &mut members, limits.members)?;
        }
    }
    for (b, outside) in global {
        emit(&mut out, b, None, &maps, (assembly, outside), &mut members, limits.members)?;
    }
    Ok(Some(out))
}

fn labels(blocks: &[&Block<'_>]) -> Result<(BTreeSet<u32>, BTreeSet<u32>)> {
    let mut nodes = BTreeSet::new();
    let mut elements = BTreeSet::new();
    for b in blocks {
        if b.keyword == "NODE" {
            for row in &b.rows {
                let f = row_fields(row.text, row.line)?;
                let value = label(f.first().map(String::as_str).unwrap_or_default(), row.line)?;
                if !nodes.insert(value) {
                    return Err(error(row.line, "duplicate local node"));
                }
            }
        } else if b.keyword == "ELEMENT" {
            let count = Type::parse(b.required("TYPE")?, b.line)?.kind.node_count();
            let mut remaining = 0usize;
            for row in &b.rows {
                let f = row_fields(row.text, row.line)?;
                let mut entries = f.len();
                if remaining == 0 {
                    let value = label(f.first().map(String::as_str).unwrap_or_default(), row.line)?;
                    if !elements.insert(value) {
                        return Err(error(row.line, "duplicate local element"));
                    }
                    remaining = count;
                    entries = entries.saturating_sub(1);
                }
                remaining = remaining.checked_sub(entries).ok_or_else(|| error(row.line, "too many element nodes"))?;
                if remaining > 0 && !row.text.ends_with(',') {
                    return Err(error(row.line, "unfinished connectivity needs continuation"));
                }
            }
            if remaining > 0 {
                return Err(error(b.line, "incomplete element connectivity"));
            }
        }
    }
    Ok((nodes, elements))
}

fn reference(
    value: &str,
    local: Option<&str>,
    nodes: bool,
    maps: &BTreeMap<String, InstanceMap>,
    assembly: (&str, bool),
    line: usize,
) -> Result<String> {
    let value = lex::unquote(value, line)?.to_ascii_uppercase();
    let stripped = value.strip_prefix(&format!("{}.", assembly.0)).filter(|rest| assembly.1 || rest.contains('.')).unwrap_or(&value);
    let (scope, item) = match stripped.split_once('.') {
        Some((scope, item)) => (Some(scope), item),
        None => (local, stripped),
    };
    let Some(scope) = scope else {
        return Ok(item.to_string());
    };
    let map = maps.get(scope).ok_or_else(|| error(line, format!("unknown instance {scope}")))?;
    if let Ok(label) = item.parse::<u32>() {
        let ids = if nodes { &map.nodes } else { &map.elements };
        return ids
            .get(&label)
            .map(ToString::to_string)
            .ok_or_else(|| error(line, format!("unknown {} {scope}.{label}", if nodes { "node" } else { "element" })));
    }
    map.names.get(item).cloned().ok_or_else(|| error(line, format!("unknown local reference {scope}.{item}")))
}

fn emit(
    out: &mut Expanded,
    b: &Block<'_>,
    instance: Option<&Instance<'_>>,
    maps: &BTreeMap<String, InstanceMap>,
    assembly: (&str, bool),
    members: &mut usize,
    limit: usize,
) -> Result<()> {
    let scope = instance.map(|i| i.name.as_str());
    let mut params = b.params.clone();
    let set_scope = b.value("INSTANCE").map(|s| id(s, b.line)).transpose()?;
    if set_scope.is_some() && !matches!(b.keyword.as_str(), "NSET" | "ELSET") {
        return Err(error(b.line, "INSTANCE parameter only supported on NSET/ELSET"));
    }
    if let Some(name) = &set_scope {
        if scope.is_some() || !maps.contains_key(name) {
            return Err(error(b.line, "invalid set INSTANCE reference"));
        }
        params.remove("INSTANCE");
    }
    if b.params.contains_key("INSTANCE") && set_scope.is_none() {
        return Err(error(b.line, "INSTANCE requires a value"));
    }
    let target_scope = set_scope.as_deref().or(scope);
    // Only interpreted blocks are rewritten. Other text remains available in the source bundle.
    let mapped = matches!(
        b.keyword.as_str(),
        "NODE"
            | "ELEMENT"
            | "NSET"
            | "ELSET"
            | "SURFACE"
            | "SOLIDSECTION"
            | "SHELLSECTION"
            | "BEAMSECTION"
            | "ORIENTATION"
            | "BOUNDARY"
            | "CLOAD"
            | "DLOAD"
            | "DSLOAD"
    );
    if mapped {
        for key in ["NSET", "ELSET", "ORIENTATION"] {
            if let Some(Some(value)) = params.get_mut(key) {
                let reference_scope = if b.keyword == "NSET" && key == "ELSET" { target_scope } else { scope };
                *value = reference(value, reference_scope, key == "NSET", maps, assembly, b.line)?;
            }
        }
        if matches!(b.keyword.as_str(), "SURFACE" | "ORIENTATION")
            && let Some(Some(value)) = params.get_mut("NAME")
        {
            *value = reference(value, scope, false, maps, assembly, b.line)?;
        }
    }
    let generated = matches!(b.keyword.as_str(), "NSET" | "ELSET") && b.params.contains_key("GENERATE");
    if b.params.get("GENERATE").is_some_and(Option::is_some) {
        return Err(error(b.line, "GENERATE is a flag and does not take a value"));
    }
    if generated {
        params.remove("GENERATE");
    }
    if b.keyword == "BEAMSECTION" && instance.is_some() {
        let count = match canonical(b.required("SECTION")?).as_str() {
            "RECT" | "PIPE" => 2,
            "CIRC" => 1,
            "I" => 7,
            "BOX" => 6,
            _ => return Err(error(b.line, "unsupported beam section in instance")),
        };
        let row = b.rows.first().ok_or_else(|| error(b.line, "beam section has no dimensions"))?;
        if row_fields(row.text, row.line)?.len() != count || b.rows.len() > 2 {
            return Err(error(b.line, "instance beam section needs a dimension row and optional direction row"));
        }
    }
    out.header(&b.keyword, &params, b.line)?;
    let mut remaining = 0usize;
    for (index, row) in b.rows.iter().enumerate() {
        check_work(index)?;
        if !mapped {
            out.push(row.text, row.line)?;
            continue;
        }
        let mut f = row_fields(row.text, row.line)?;
        match b.keyword.as_str() {
            "NODE" if instance.is_some() => {
                if f.is_empty() || f.len() > 4 {
                    return Err(error(row.line, "NODE needs label and up to three coordinates"));
                }
                if b.value("SYSTEM").is_some_and(|v| canonical(v) != "R") {
                    return Err(error(b.line, "instance NODE requires Cartesian coordinates"));
                }
                let first = reference(f.first().map(String::as_str).unwrap_or_default(), scope, true, maps, assembly, row.line)?;
                let mut xyz = [0.0; 3];
                for (i, x) in xyz.iter_mut().enumerate() {
                    if let Some(v) = f.get(i + 1).filter(|v| !v.is_empty()) {
                        *x = number(v, row.line)?;
                    }
                }
                if let Some(i) = instance {
                    xyz = i.transform.point(xyz, row.line)?;
                }
                f = vec![first];
                f.extend(xyz.iter().map(ToString::to_string));
            }
            "ELEMENT" => {
                let count = Type::parse(b.required("TYPE")?, b.line)?.kind.node_count();
                for (j, value) in f.iter_mut().enumerate() {
                    if remaining == 0 && j == 0 {
                        *value = reference(value, scope, false, maps, assembly, row.line)?;
                        remaining = count;
                    } else {
                        *value = reference(value, scope, true, maps, assembly, row.line)?;
                        remaining = remaining.checked_sub(1).ok_or_else(|| error(row.line, "too many element nodes"))?;
                    }
                }
            }
            "NSET" | "ELSET" => {
                let node = b.keyword == "NSET";
                if generated {
                    if !(2..=3).contains(&f.len()) {
                        return Err(error(row.line, "GENERATE needs first,last,optional increment"));
                    }
                    let start = label(f.first().map(String::as_str).unwrap_or_default(), row.line)?;
                    let end = label(f.get(1).map(String::as_str).unwrap_or_default(), row.line)?;
                    let increment = f.get(2).map(|v| label(v, row.line)).transpose()?.unwrap_or(1);
                    if end < start || !(end - start).is_multiple_of(increment) {
                        return Err(error(row.line, "invalid GENERATE range"));
                    }
                    let count = ((end - start) / increment) as usize + 1;
                    *members = members.saturating_add(count);
                    if *members > limit {
                        return Err(error(row.line, "set expansion limit reached"));
                    }
                    for j in 0..count {
                        check_work(j)?;
                        let label = start + (j as u32) * increment;
                        out.push(&reference(&label.to_string(), target_scope, node, maps, assembly, row.line)?, row.line)?;
                    }
                    continue;
                }
                *members = members.saturating_add(f.len());
                if *members > limit {
                    return Err(error(row.line, "set member limit reached"));
                }
                for v in &mut f {
                    *v = reference(v, target_scope, node, maps, assembly, row.line)?;
                }
            }
            "BOUNDARY" | "CLOAD" | "DLOAD" | "DSLOAD" => {
                if let Some(v) = f.first_mut().filter(|v| !v.is_empty()) {
                    *v = reference(v, scope, matches!(b.keyword.as_str(), "BOUNDARY" | "CLOAD"), maps, assembly, row.line)?;
                }
            }
            "SURFACE" => {
                if let Some(v) = f.first_mut() {
                    *v = reference(v, scope, b.value("TYPE").is_some_and(|s| canonical(s) == "NODE"), maps, assembly, row.line)?;
                }
            }
            "ORIENTATION" if instance.is_some() && index == 0 => {
                if !matches!(f.len(), 6 | 9) {
                    return Err(error(row.line, "orientation requires 6 or 9 coordinates"));
                }
                if b.value("DEFINITION").is_some_and(|s| canonical(s) != "COORDINATES") {
                    return Err(error(b.line, "only coordinate orientations can be transformed"));
                }
                let rectangular = b.value("SYSTEM").is_none_or(|s| matches!(canonical(s).as_str(), "R" | "RECTANGULAR"));
                if rectangular && f.len() == 6 {
                    f.extend(["0".into(), "0".into(), "0".into()]);
                }
                for chunk in f.chunks_mut(3) {
                    let mut xyz = [0.0; 3];
                    for (x, v) in xyz.iter_mut().zip(chunk.iter()) {
                        *x = number(v, row.line)?;
                    }
                    if let Some(i) = instance {
                        xyz = i.transform.point(xyz, row.line)?;
                    }
                    for (v, x) in chunk.iter_mut().zip(xyz) {
                        *v = x.to_string();
                    }
                }
            }
            "BEAMSECTION" if instance.is_some() && index == 1 => {
                if f.len() != 3 {
                    return Err(error(row.line, "beam direction needs three coordinates"));
                }
                let mut xyz = [0.0; 3];
                for (x, v) in xyz.iter_mut().zip(&f) {
                    *x = number(v, row.line)?;
                }
                if let Some(i) = instance {
                    xyz = i.transform.vector(xyz);
                }
                f = xyz.iter().map(ToString::to_string).collect();
            }
            _ => {
                out.push(row.text, row.line)?;
                continue;
            }
        }
        let mut text = f.join(", ");
        if b.keyword == "ELEMENT" && remaining > 0 {
            text.push(',');
        }
        out.push(&text, row.line)?;
    }
    if b.keyword == "BEAMSECTION"
        && b.rows.len() == 1
        && let Some(instance) = instance
    {
        let direction = instance.transform.vector([0.0, 0.0, -1.0]);
        out.push(&direction.iter().map(ToString::to_string).collect::<Vec<_>>().join(", "), b.line)?;
    }
    Ok(())
}
