//! Flat model-data import from the public Abaqus keyword reference.
//!
//! Syntax: Dassault Systemes, Abaqus 2017, Input Syntax Rules and the *NODE, *ELEMENT,
//! *NSET, *ELSET, *SURFACE, *MATERIAL, *ELASTIC, *DENSITY and section keyword references:
//! <https://abaqus-docs.mit.edu/2017/English/SIMACAEMODRefMap/simamod-c-inputsyntax.htm>.
//! Only original handwritten or FieldCraft-writer-generated decks are used in tests.
//! No history, assembly, includes or raw-keyword passthrough is imported. Numbers in a deck
//! have no intrinsic units: callers choose a consistent preset; the model stores SI.
//! An unread option terminates a material-property block: later *ELASTIC/*DENSITY must have a
//! new *MATERIAL definition. This conservative boundary refuses ambiguous property association.

mod lex;
mod types;

use crate::Result;
use fieldcraft_model::progress;
use fieldcraft_model::units::{Preset, Quantity};
use fieldcraft_model::{BeamProfile, DisplaySystem, Elastic, ElementEdge, ElementFace, Material, Model, Section, SectionKind, Set};
use lex::{Block, canonical, error, fields, unquote};
use std::collections::{BTreeMap, BTreeSet};
use types::{Family, Type};

/// Maximum accepted source size and total expanded set membership. These are import bounds,
/// independent of the model's larger limits; a small file cannot request an enormous range.
pub const MAX_INPUT_BYTES: usize = 64 * 1024 * 1024;
pub const MAX_RECORDS: usize = 1_000_000;
pub const MAX_MEMBERS: usize = 2_000_000;
pub const MAX_UNKNOWN_SHOWN: usize = 100;
const MAX_SET_DEPTH: usize = 128;
/// Constant properties/sections need at most ten entries today; bound before collecting data.
const MAX_PROPERTY_VALUES: usize = 32;

/// An unread option; line is its original one-based physical keyword line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Unknown {
    pub keyword: String,
    pub line: usize,
}

/// A new model, without history. Unknown output is bounded; `unknown_count` counts all options.
#[derive(Debug, Clone)]
pub struct Read {
    pub model: Model,
    pub unknown: Vec<Unknown>,
    pub unknown_count: usize,
    pub warnings: Vec<String>,
}

/// Read a deck whose numbers are explicitly interpreted in SI.
pub fn read(text: &str) -> Result<Read> {
    read_in(text, Preset::MKgS)
}

/// Read a deck in the caller's chosen consistent units, converting each quantity at insertion.
pub fn read_in(text: &str, preset: Preset) -> Result<Read> {
    parse(text, preset, MAX_MEMBERS)
}

/// Cooperate with the host without a thread-local lookup for every small token.
fn check_work(index: usize) -> Result<()> {
    if index.is_multiple_of(256) {
        progress::check()?;
    }
    Ok(())
}

fn number(s: &str, line: usize) -> Result<f64> {
    let value = s.replace(['d', 'D'], "E").parse::<f64>().map_err(|_| error(line, format!("invalid number `{s}`")))?;
    if !value.is_finite() {
        return Err(error(line, format!("non-finite number `{s}`")));
    }
    Ok(value)
}
fn label(s: &str, line: usize) -> Result<u32> {
    let n = s.parse::<u32>().map_err(|_| error(line, format!("invalid label `{s}`")))?;
    if n == 0 || n > fieldcraft_model::MAX_LABEL {
        return Err(error(line, format!("label {n} is outside 1..={}", fieldcraft_model::MAX_LABEL)));
    }
    Ok(n)
}
fn named(value: &str, line: usize) -> Result<String> {
    let name = unquote(value, line)?;
    fieldcraft_model::check_name("imported", &name).map_err(|e| error(line, e.to_string()))?;
    Ok(name)
}
fn key(value: &str, line: usize) -> Result<String> {
    Ok(named(value, line)?.to_ascii_uppercase())
}
fn row_fields(text: &str, line: usize) -> Result<Vec<String>> {
    let mut values = fields(text, line)?;
    while values.last().is_some_and(|s| s.is_empty()) {
        values.pop();
    }
    Ok(values)
}
fn values(block: &Block<'_>) -> Result<Vec<f64>> {
    let mut out = Vec::new();
    for (index, row) in block.rows.iter().enumerate() {
        check_work(index)?;
        for token in row_fields(row.text, row.line)? {
            if out.len() >= MAX_PROPERTY_VALUES {
                return Err(error(row.line, "constant property/section exceeds 32 numeric entries"));
            }
            out.push(number(&token, row.line)?);
        }
    }
    Ok(out)
}
fn quantity(units: &DisplaySystem, q: Quantity, value: f64, line: usize) -> Result<f64> {
    units.to_si(q, value).map_err(|e| error(line, e.to_string()))
}
fn name_entry(names: &mut BTreeMap<String, String>, name: &str, line: usize) -> Result<String> {
    let canonical = name.to_ascii_uppercase();
    if let Some(old) = names.get(&canonical) {
        if old != name {
            return Err(error(line, format!("names `{old}` and `{name}` differ only in case")));
        }
    } else {
        if names.len() >= fieldcraft_model::MAX_DEFINITIONS {
            return Err(error(line, "definition limit reached"));
        }
        names.insert(canonical.clone(), name.to_string());
    }
    Ok(canonical)
}

#[derive(Debug, Clone)]
enum Member {
    Label(u32, usize),
    Ref(String, usize),
}
#[derive(Default)]
struct RawSet {
    members: Vec<Member>,
    line: usize,
}
struct ElementRecord {
    label: u32,
    typ: Type,
    nodes: Vec<u32>,
    line: usize,
}
struct Surface {
    name: String,
    line: usize,
    rows: Vec<(String, String, usize)>,
}
struct SectionRecord {
    set: String,
    material: String,
    keyword: String,
    params: BTreeMap<String, Option<String>>,
    values: Vec<f64>,
    line: usize,
}

fn append_member(set: &mut RawSet, member: Member, budget: &mut usize, limit: usize, line: usize) -> Result<()> {
    *budget = budget.saturating_add(1);
    check_work(*budget)?;
    if *budget > limit {
        return Err(error(line, "expanded set member limit reached"));
    }
    set.members.push(member);
    Ok(())
}

fn expand(
    name: &str,
    raw: &BTreeMap<String, RawSet>,
    cache: &mut BTreeMap<String, Vec<u32>>,
    active: &mut BTreeSet<String>,
    depth: usize,
    work: &mut usize,
    limit: usize,
    line: usize,
) -> Result<Vec<u32>> {
    progress::check()?;
    if depth > MAX_SET_DEPTH {
        return Err(error(line, "set nesting limit reached"));
    }
    if let Some(found) = cache.get(name) {
        *work = work.saturating_add(found.len());
        if *work > limit {
            return Err(error(line, "set expansion work limit reached"));
        }
        return Ok(found.clone());
    }
    let definition = raw.get(name).ok_or_else(|| error(line, format!("missing set `{name}`")))?;
    if !active.insert(name.to_string()) {
        return Err(error(line, format!("set `{name}` contains a reference cycle")));
    }
    let mut out = BTreeSet::new();
    for (index, member) in definition.members.iter().enumerate() {
        check_work(index)?;
        match member {
            Member::Label(value, at) => {
                *work = work.saturating_add(1);
                if *work > limit {
                    return Err(error(*at, "set expansion work limit reached"));
                }
                out.insert(*value);
            }
            Member::Ref(other, at) => {
                for (index, label) in expand(other, raw, cache, active, depth + 1, work, limit, *at)?.into_iter().enumerate() {
                    check_work(index)?;
                    out.insert(label);
                }
            }
        }
        if out.len() > limit {
            return Err(error(definition.line, "expanded set member limit reached"));
        }
    }
    active.remove(name);
    let mut labels = Vec::with_capacity(out.len());
    for (index, label) in out.into_iter().enumerate() {
        check_work(index)?;
        labels.push(label);
    }
    let out = labels;
    *work = work.saturating_add(out.len());
    if *work > limit {
        return Err(error(line, "set expansion work limit reached"));
    }
    cache.insert(name.to_string(), out.clone());
    Ok(out)
}

fn parse(text: &str, preset: Preset, member_limit: usize) -> Result<Read> {
    if text.len() > MAX_INPUT_BYTES {
        return Err(error(1, "input exceeds 64 MiB"));
    }
    let blocks = lex::blocks(text, MAX_RECORDS)?;
    let units = DisplaySystem::new(preset);
    let mut model = Model::default();
    model.set_display_units(Some(units.clone()));
    let mut names = BTreeMap::new();
    let mut mat_names = BTreeMap::new();
    let mut materials: BTreeMap<String, (Material, usize)> = BTreeMap::new();
    let (mut nsets, mut elsets): (BTreeMap<String, RawSet>, BTreeMap<String, RawSet>) = (BTreeMap::new(), BTreeMap::new());
    let mut node_from_element_sets = Vec::new();
    let mut elements = Vec::new();
    let mut element_labels = BTreeSet::new();
    let mut surfaces = Vec::new();
    let mut sections = Vec::new();
    let (mut unknown, mut unknown_count) = (Vec::new(), 0usize);
    let mut warnings = Vec::new();
    let (mut material, mut history) = (None, None);
    let mut budget = 0usize;
    for block in &blocks {
        progress::check()?;
        let k = block.keyword.as_str();
        let line = block.line;
        if ["INCLUDE", "SYSTEM", "PART", "ASSEMBLY", "INSTANCE", "ENDPART", "ENDASSEMBLY", "ENDINSTANCE"].contains(&k) {
            return Err(error(line, format!("*{k} is not supported by the flat model-data reader")));
        }
        if k == "STEP" {
            if history.is_some() {
                return Err(error(line, "nested *STEP is not valid; history is not imported"));
            }
            history = Some(line);
            material = None;
        }
        if k == "ENDSTEP" && history.is_none() {
            return Err(error(line, "*END STEP without an open *STEP"));
        }
        let ignored = history.is_some() || k == "ENDSTEP";
        if ignored {
            unknown_count += 1;
            if unknown.len() < MAX_UNKNOWN_SHOWN {
                unknown.push(Unknown { keyword: k.to_string(), line });
            }
            if k == "ENDSTEP" {
                history = None;
            }
            continue;
        }
        if !["MATERIAL", "ELASTIC", "DENSITY"].contains(&k) {
            material = None;
        }
        match k {
            "HEADING" => {
                block.only(&[])?;
            }
            "NODE" => {
                block.only(&["NSET", "SYSTEM"])?;
                if block.value("SYSTEM").is_some_and(|s| canonical(s) != "R") {
                    return Err(error(line, "*NODE coordinate-system transforms are unsupported; Cartesian SYSTEM=R only"));
                }
                let set = block.value("NSET").map(|s| named(s, line)).transpose()?.map(|name| name_entry(&mut names, &name, line)).transpose()?;
                for (index, row) in block.rows.iter().enumerate() {
                    check_work(index)?;
                    let f = row_fields(row.text, row.line)?;
                    if f.is_empty() || f.len() > 4 {
                        return Err(error(row.line, "*NODE needs a label and up to three coordinates; omitted coordinates are zero"));
                    }
                    let l = label(f.first().map(String::as_str).unwrap_or_default(), row.line)?;
                    if model.node(l).is_some() {
                        return Err(error(row.line, format!("duplicate node {l}")));
                    }
                    let mut xyz = [0.; 3];
                    for (i, value) in f.iter().skip(1).enumerate() {
                        if let Some(x) = xyz.get_mut(i) {
                            *x = quantity(&units, Quantity::Length, if value.is_empty() { 0. } else { number(value, row.line)? }, row.line)?;
                        }
                    }
                    model.set_node(l, xyz).map_err(|e| error(row.line, e.to_string()))?;
                    if let Some(name) = &set {
                        let raw = nsets.entry(name.clone()).or_default();
                        raw.line = line;
                        append_member(raw, Member::Label(l, row.line), &mut budget, member_limit, row.line)?;
                    }
                }
            }
            "ELEMENT" => {
                block.only(&["TYPE", "ELSET"])?;
                let typ = Type::parse(block.required("TYPE")?, line)?;
                let set = block.value("ELSET").map(|s| named(s, line)).transpose()?.map(|name| name_entry(&mut names, &name, line)).transpose()?;
                let mut pending: Option<ElementRecord> = None;
                for (index, row) in block.rows.iter().enumerate() {
                    check_work(index)?;
                    let mut f = row_fields(row.text, row.line)?.into_iter();
                    if pending.is_none() {
                        let l = label(&f.next().unwrap_or_default(), row.line)?;
                        if !element_labels.insert(l) {
                            return Err(error(row.line, format!("duplicate element {l}")));
                        }
                        pending = Some(ElementRecord { label: l, typ, nodes: Vec::new(), line: row.line });
                    }
                    if let Some(e) = pending.as_mut() {
                        for n in f {
                            e.nodes.push(label(&n, row.line)?);
                            if e.nodes.len() > typ.kind.node_count() {
                                return Err(error(row.line, format!("element {} has too many nodes", e.label)));
                            }
                        }
                    }
                    if pending.as_ref().is_some_and(|e| e.nodes.len() == typ.kind.node_count()) {
                        if let Some(e) = pending.take() {
                            if let Some(name) = &set {
                                let raw = elsets.entry(name.clone()).or_default();
                                raw.line = line;
                                append_member(raw, Member::Label(e.label, e.line), &mut budget, member_limit, e.line)?;
                            }
                            elements.push(e);
                        }
                    } else if !row.text.ends_with(',') {
                        return Err(error(row.line, "unfinished element connectivity needs a trailing-comma continuation"));
                    }
                }
                if let Some(e) = pending {
                    return Err(error(e.line, format!("element {} has incomplete connectivity", e.label)));
                }
                if block.rows.is_empty() {
                    return Err(error(line, "*ELEMENT has no data"));
                }
            }
            "NSET" | "ELSET" => {
                block.only(if k == "NSET" { &["NSET", "GENERATE", "ELSET"] } else { &["ELSET", "GENERATE"] })?;
                let name = named(block.required(k)?, line)?;
                let name = name_entry(&mut names, &name, line)?;
                if k == "NSET" && block.value("ELSET").is_some() {
                    if block.params.contains_key("GENERATE") || !block.rows.is_empty() {
                        return Err(error(line, "NSET from ELSET has no data or GENERATE"));
                    }
                    node_from_element_sets.push((name, key(block.required("ELSET")?, line)?, line));
                    continue;
                }
                let raw = if k == "NSET" { nsets.entry(name).or_default() } else { elsets.entry(name).or_default() };
                raw.line = line;
                for (index, row) in block.rows.iter().enumerate() {
                    check_work(index)?;
                    let f = row_fields(row.text, row.line)?;
                    if block.params.contains_key("GENERATE") {
                        if f.len() < 2 || f.len() > 3 {
                            return Err(error(row.line, "GENERATE needs first,last,optional increment"));
                        }
                        let start = label(f.first().map(String::as_str).unwrap_or_default(), row.line)?;
                        let end = label(f.get(1).map(String::as_str).unwrap_or_default(), row.line)?;
                        let increment = f.get(2).map(|s| label(s, row.line)).transpose()?.unwrap_or(1);
                        if end < start || (end - start) % increment != 0 {
                            return Err(error(row.line, "GENERATE needs an ascending exactly divisible range"));
                        }
                        let count = u64::from((end - start) / increment) + 1;
                        if count > member_limit.saturating_sub(budget) as u64 {
                            return Err(error(row.line, "GENERATE exceeds set member limit"));
                        }
                        for i in 0..count {
                            append_member(
                                raw,
                                Member::Label((u64::from(start) + i * u64::from(increment)) as u32, row.line),
                                &mut budget,
                                member_limit,
                                row.line,
                            )?;
                        }
                    } else {
                        for value in f {
                            let member = if value.bytes().next().is_some_and(|c| c.is_ascii_digit() || c == b'-' || c == b'+') {
                                Member::Label(label(&value, row.line)?, row.line)
                            } else {
                                Member::Ref(key(&value, row.line)?, row.line)
                            };
                            append_member(raw, member, &mut budget, member_limit, row.line)?;
                        }
                    }
                }
            }
            "SURFACE" => {
                block.only(&["TYPE", "NAME"])?;
                if block.value("TYPE").is_some_and(|t| canonical(t) != "ELEMENT") {
                    unknown_count += 1;
                    if unknown.len() < MAX_UNKNOWN_SHOWN {
                        unknown.push(Unknown { keyword: k.to_string(), line });
                    }
                    if warnings.len() < MAX_UNKNOWN_SHOWN {
                        warnings.push(format!(
                            "line {line}: *SURFACE TYPE={} is unread; only element-based surfaces are imported",
                            block.value("TYPE").unwrap_or_default()
                        ));
                    }
                    continue;
                }
                let name = named(block.required("NAME")?, line)?;
                name_entry(&mut names, &name, line)?;
                let mut rows = Vec::new();
                for (index, row) in block.rows.iter().enumerate() {
                    check_work(index)?;
                    let f = row_fields(row.text, row.line)?;
                    if f.len() != 2 {
                        return Err(error(row.line, "an element surface needs an element/set and face identifier"));
                    }
                    rows.push((f.first().cloned().unwrap_or_default(), canonical(f.get(1).map(String::as_str).unwrap_or_default()), row.line));
                }
                surfaces.push(Surface { name, line, rows });
            }
            "MATERIAL" => {
                block.only(&["NAME"])?;
                if !block.rows.is_empty() {
                    return Err(error(line, "*MATERIAL has no data lines"));
                }
                let name = named(block.required("NAME")?, line)?;
                let name = name_entry(&mut mat_names, &name, line)?;
                if materials.insert(name.clone(), (Material::default(), line)).is_some() {
                    return Err(error(line, "duplicate material"));
                }
                material = Some(name);
            }
            "ELASTIC" | "DENSITY" => {
                let name = material
                    .as_ref()
                    .ok_or_else(|| error(line, format!("*{k} outside a supported *MATERIAL block (an unread keyword ends material context)")))?;
                block.only(if k == "ELASTIC" { &["TYPE"] } else { &[] })?;
                let v = values(block)?;
                let m = materials.get_mut(name).ok_or_else(|| error(line, "missing active material"))?;
                if k == "DENSITY" {
                    if v.len() != 1 || m.0.density.is_some() {
                        return Err(error(line, "constant *DENSITY needs exactly one value and one definition"));
                    }
                    m.0.density = Some(quantity(&units, Quantity::Density, v.first().copied().unwrap_or_default(), line)?);
                } else {
                    if m.0.elastic.is_some() {
                        return Err(error(line, "duplicate *ELASTIC"));
                    }
                    let typ = canonical(block.value("TYPE").unwrap_or("ISOTROPIC"));
                    m.0.elastic = Some(match typ.as_str() {
                        "ISOTROPIC" if v.len() == 2 => Elastic::Isotropic {
                            young: quantity(&units, Quantity::Pressure, v.first().copied().unwrap_or_default(), line)?,
                            poisson: v.get(1).copied().unwrap_or_default(),
                        },
                        "ENGINEERINGCONSTANTS" if v.len() == 9 => {
                            let mut young = [0.; 3];
                            let mut poisson = [0.; 3];
                            let mut shear = [0.; 3];
                            for i in 0..3 {
                                if let Some(x) = young.get_mut(i) {
                                    *x = quantity(&units, Quantity::Pressure, v.get(i).copied().unwrap_or_default(), line)?;
                                }
                                if let Some(x) = poisson.get_mut(i) {
                                    *x = v.get(i + 3).copied().unwrap_or_default();
                                }
                                if let Some(x) = shear.get_mut(i) {
                                    *x = quantity(&units, Quantity::Pressure, v.get(i + 6).copied().unwrap_or_default(), line)?;
                                }
                            }
                            Elastic::EngineeringConstants { young, poisson, shear }
                        }
                        _ => return Err(error(line, format!("unsupported or incomplete constant *ELASTIC TYPE={typ}"))),
                    });
                }
            }
            "SOLIDSECTION" | "SHELLSECTION" | "BEAMSECTION" => {
                block.only(match k {
                    "SHELLSECTION" => &["ELSET", "MATERIAL", "OFFSET"],
                    "BEAMSECTION" => &["ELSET", "MATERIAL", "SECTION"],
                    _ => &["ELSET", "MATERIAL"],
                })?;
                let set = key(block.required("ELSET")?, line)?;
                let material = key(block.required("MATERIAL")?, line)?;
                sections.push(SectionRecord { set, material, keyword: k.to_string(), params: block.params.clone(), values: values(block)?, line });
            }
            _ => {
                unknown_count += 1;
                if unknown.len() < MAX_UNKNOWN_SHOWN {
                    unknown.push(Unknown { keyword: k.to_string(), line });
                }
            }
        }
    }
    if let Some(line) = history {
        return Err(error(line, "*STEP has no matching *END STEP; history is not imported"));
    }
    let mut types = BTreeMap::new();
    for (index, e) in elements.into_iter().enumerate() {
        check_work(index)?;
        model.set_element(e.label, e.typ.kind, e.nodes).map_err(|err| error(e.line, err.to_string()))?;
        types.insert(e.label, (e.typ, e.line));
    }
    for (name, (m, line)) in materials {
        progress::check()?;
        let spelling = mat_names.get(&name).ok_or_else(|| error(line, "missing material spelling"))?;
        model.add_material(spelling, m).map_err(|e| error(line, e.to_string()))?;
    }
    let mut work = 0usize;
    let mut expanded_elements = BTreeMap::new();
    let mut active = BTreeSet::new();
    for (name, raw) in &elsets {
        let labels = expand(name, &elsets, &mut expanded_elements, &mut active, 0, &mut work, member_limit, raw.line)?;
        for (index, l) in labels.iter().enumerate() {
            check_work(index)?;
            if model.element(*l).is_none() {
                return Err(error(raw.line, format!("element set `{name}` names missing element {l}")));
            }
        }
    }
    for (name, elset, line) in node_from_element_sets {
        progress::check()?;
        let elements = expanded_elements.get(&elset).ok_or_else(|| error(line, format!("missing element set `{elset}`")))?;
        let raw = nsets.entry(name).or_default();
        raw.line = line;
        for l in elements {
            if let Some(e) = model.element(*l) {
                for n in &e.nodes {
                    append_member(raw, Member::Label(*n, line), &mut budget, member_limit, line)?;
                }
            }
        }
    }
    let mut expanded_nodes = BTreeMap::new();
    for (name, raw) in &nsets {
        let labels = expand(name, &nsets, &mut expanded_nodes, &mut active, 0, &mut work, member_limit, raw.line)?;
        for (index, l) in labels.iter().enumerate() {
            check_work(index)?;
            if model.node(*l).is_none() {
                return Err(error(raw.line, format!("node set `{name}` names missing node {l}")));
            }
        }
    }
    let mut sets: BTreeMap<String, Set> = BTreeMap::new();
    for (name, labels) in &expanded_nodes {
        progress::check()?;
        sets.entry(name.clone()).or_default().nodes = labels.clone();
    }
    for (name, labels) in &expanded_elements {
        progress::check()?;
        sets.entry(name.clone()).or_default().elements = labels.clone();
    }
    for surface in surfaces {
        progress::check()?;
        let name = surface.name.to_ascii_uppercase();
        let set = sets.entry(name).or_default();
        for (target, face, line) in surface.rows {
            progress::check()?;
            let labels = if target.bytes().next().is_some_and(|c| c.is_ascii_digit()) {
                vec![label(&target, line)?]
            } else {
                expanded_elements
                    .get(&key(&target, line)?)
                    .cloned()
                    .ok_or_else(|| error(line, format!("surface `{}` names missing element set `{target}`", surface.name)))?
            };
            budget = budget.saturating_add(labels.len());
            if budget > member_limit {
                return Err(error(line, "surface expansion limit reached"));
            }
            for (index, l) in labels.into_iter().enumerate() {
                check_work(index)?;
                let (typ, _) = types.get(&l).ok_or_else(|| error(line, format!("surface `{}` names missing element {l}", surface.name)))?;
                let shell = typ.family == Family::Shell;
                let id = if shell {
                    match face.as_str() {
                        "SPOS" => 1,
                        "SNEG" => 2,
                        _ => return Err(error(line, "shell surfaces need SPOS or SNEG")),
                    }
                } else {
                    face.strip_prefix('S')
                        .and_then(|v| v.parse::<u8>().ok())
                        .filter(|v| *v > 0)
                        .ok_or_else(|| error(line, format!("invalid surface face `{face}`")))?
                };
                let valid = if matches!(typ.family, Family::Stress | Family::Strain | Family::Axis) {
                    if super::edge_id(typ.kind, id).is_some() {
                        set.edges.push(ElementEdge { element: l, edge: id });
                        true
                    } else {
                        false
                    }
                } else if super::face_id(typ.kind, id).is_some() {
                    set.faces.push(ElementFace { element: l, face: id });
                    true
                } else {
                    false
                };
                if !valid {
                    return Err(error(line, format!("element {l} has no face/edge `{face}`")));
                }
            }
        }
        if set.is_empty() {
            return Err(error(surface.line, format!("empty surface `{}`", surface.name)));
        }
    }
    for (name, set) in sets {
        progress::check()?;
        let spelling = names.get(&name).ok_or_else(|| error(1, "missing set spelling"))?;
        model.add_set(spelling, set).map_err(|e| error(1, e.to_string()))?;
    }
    let mut has_planar = false;
    let mut assigned = BTreeMap::new();
    for section in sections {
        progress::check()?;
        let line = section.line;
        let labels = expanded_elements
            .get(&section.set)
            .filter(|l| !l.is_empty())
            .ok_or_else(|| error(line, format!("section names missing/empty element set `{}`", section.set)))?;
        let mut typ = None;
        for (index, label) in labels.iter().enumerate() {
            check_work(index)?;
            if let Some(previous) = assigned.insert(*label, line) {
                return Err(error(line, format!("element {label} has overlapping section assignments (previous line {previous})")));
            }
            let (current, _) = types.get(label).ok_or_else(|| error(line, format!("section names missing element {label}")))?;
            if typ.is_some_and(|old: Type| old.family != current.family || old.formulation != current.formulation) {
                return Err(error(line, format!("section on `{}` mixes element types/formulations", section.set)));
            }
            typ = Some(*current);
        }
        let typ = typ.ok_or_else(|| error(line, "section has no elements"))?;
        if matches!(typ.family, Family::Stress | Family::Strain | Family::Axis) {
            has_planar = true;
        }
        let kind = section_kind(&section, typ, &units)?;
        let material = mat_names
            .get(&section.material)
            .filter(|name| model.material(name).is_some())
            .ok_or_else(|| error(line, format!("section names missing material `{}`", section.material)))?;
        let set = names.get(&section.set).ok_or_else(|| error(line, "missing section set spelling"))?;
        let name = format!("Section{}", model.sections().len() + 1);
        model
            .add_section(&name, Section { set: set.clone(), material: material.clone(), kind, formulation: typ.formulation })
            .map_err(|e| error(line, e.to_string()))?;
    }
    let mut unassigned = 0usize;
    for (index, (label, (typ, line))) in types.iter().enumerate() {
        check_work(index)?;
        if !has_planar && typ.family == Family::PlanarRod {
            return Err(error(
                *line,
                "T2D trusses without a planar continuum section cannot retain their 2D constraint in this model; import refused",
            ));
        }
        if has_planar && typ.family == Family::Rod {
            return Err(error(*line, "T3D trusses with a planar continuum would gain a different global constraint; import refused"));
        }
        unassigned += usize::from(!assigned.contains_key(label));
    }
    progress::check()?;
    if unassigned != 0 {
        warnings.push(format!(
            "{unassigned} mesh elements have no section: their input family/formulation is not stored on the mesh element. Assign compatible sections before analysis; reduced/hybrid/shell/plane context is not retained without them."
        ));
    }
    Ok(Read { model, unknown, unknown_count, warnings })
}

fn section_kind(raw: &SectionRecord, typ: Type, units: &DisplaySystem) -> Result<SectionKind> {
    let v = &raw.values;
    let line = raw.line;
    let len = |value: f64| quantity(units, Quantity::Length, value, line);
    if raw.keyword == "BEAMSECTION" {
        if typ.family != Family::Beam {
            return Err(error(line, "*BEAM SECTION needs B31 elements"));
        }
        let shape = canonical(raw.params.get("SECTION").and_then(|v| v.as_deref()).ok_or_else(|| error(line, "*BEAM SECTION needs SECTION="))?);
        let count = match shape.as_str() {
            "RECT" | "PIPE" => 2,
            "CIRC" => 1,
            "I" => 7,
            "BOX" => 6,
            _ => return Err(error(line, format!("unsupported beam section {shape}"))),
        };
        if v.len() != count && v.len() != count + 3 {
            return Err(error(line, "beam section needs geometric data and an optional three-component orientation"));
        }
        let mut dims = Vec::new();
        for x in v.iter().take(count) {
            dims.push(len(*x)?);
        }
        let at = |i| dims.get(i).copied().unwrap_or_default();
        let profile = match shape.as_str() {
            "RECT" => BeamProfile::Rectangle { width: at(0), height: at(1) },
            "PIPE" => BeamProfile::Pipe { radius: at(0), thickness: at(1) },
            "CIRC" => BeamProfile::Circle { radius: at(0) },
            "I" => {
                if v.first().copied().unwrap_or_default() != v.get(1).copied().unwrap_or_default() / 2.0 || at(2) != at(3) || at(4) != at(5) {
                    return Err(error(line, "only centred doubly symmetric I profiles are supported"));
                }
                BeamProfile::I { height: at(1), width: at(2), flange: at(4), web: at(6) }
            }
            _ => {
                if dims.iter().skip(3).any(|x| *x != at(2)) {
                    return Err(error(line, "only equal-thickness BOX profiles are supported"));
                }
                BeamProfile::Box { width: at(0), height: at(1), thickness: at(2) }
            }
        };
        let orientation = if v.len() == count {
            [0., 0., -1.]
        } else {
            [v.get(count).copied().unwrap_or_default(), v.get(count + 1).copied().unwrap_or_default(), v.get(count + 2).copied().unwrap_or_default()]
        };
        return Ok(SectionKind::Beam { profile, orientation });
    }
    if raw.keyword == "SHELLSECTION" {
        if typ.family != Family::Shell || v.is_empty() || v.len() > 2 {
            return Err(error(line, "*SHELL SECTION needs shell elements, thickness and optional point count"));
        }
        let points = v.get(1).copied().unwrap_or(5.0);
        if points.fract() != 0.0 || points < 1.0 || points > f64::from(fieldcraft_model::sections::MAX_SHELL_POINTS) {
            return Err(error(line, "invalid shell integration point count"));
        }
        let offset = raw.params.get("OFFSET").and_then(|x| x.as_deref()).map(|x| number(x, line)).transpose()?.unwrap_or_default();
        return Ok(SectionKind::Shell { thickness: len(v.first().copied().unwrap_or_default())?, offset, points: points as u32 });
    }
    if matches!(typ.family, Family::Beam | Family::Shell) {
        return Err(error(line, "*SOLID SECTION does not match these elements"));
    }
    if v.len() > 1 {
        return Err(error(line, "constant *SOLID SECTION accepts at most one attribute"));
    }
    let attribute = v.first().copied().unwrap_or(1.0);
    let amount = match typ.family {
        Family::Stress | Family::Strain => len(attribute)?,
        Family::Rod | Family::PlanarRod => quantity(units, Quantity::Area, attribute, line)?,
        _ => attribute,
    };
    Ok(typ.section(amount))
}

#[cfg(test)]
mod tests;
