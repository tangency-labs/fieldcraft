//! Clean-room ASCII MSH 4.1 import, from the Gmsh Reference Manual, sections 10.1 and 10.2:
//! <https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format>. No Gmsh implementation is used.
//!
//! Coordinates are unitless here. The caller supplies their length unit at its model boundary.
//! Boundary records select existing element sides; they never add an accidental shell layer.

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_model::{ElementEdge, ElementFace, ElementKind, MAX_COORD, MAX_DEFINITIONS, MAX_ELEMENTS, MAX_LABEL, MAX_NODES, Model, Set};

use crate::{Error, Result};

/// Largest text file accepted, and longest individual line (including ignored sections).
pub const MAX_BYTES: u64 = 256 << 20;
const MAX_LINE: usize = 64 << 10;

/// Work and storage bounds, also usable by callers with a smaller resource budget.
#[derive(Clone, Copy, Debug)]
pub struct Limits {
    pub nodes: usize,
    pub elements: usize,
    /// Global cap for physical memberships and for side matches, not a per-set allowance.
    pub members: usize,
}

impl Default for Limits {
    fn default() -> Self {
        Self { nodes: MAX_NODES, elements: MAX_ELEMENTS, members: MAX_ELEMENTS }
    }
}

/// How to treat surface records that do not belong to a volume element face.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Unmatched {
    #[default]
    Refuse,
    Shell,
}

/// A checked mesh and physical sets. Coordinates retain the file's unitless values.
#[derive(Debug)]
pub struct Imported {
    pub model: Model,
    pub skipped_sections: Vec<String>,
}

fn error(line: usize, message: impl std::fmt::Display) -> Error {
    Error::Invalid(format!("MSH line {line}: {message}"))
}

struct Lines<'a> {
    rest: &'a str,
    line: usize,
}

impl<'a> Lines<'a> {
    fn next(&mut self) -> Result<Option<(usize, &'a str)>> {
        if self.rest.is_empty() {
            return Ok(None);
        }
        let (line, rest) = self.rest.split_once('\n').unwrap_or((self.rest, ""));
        if line.len() > MAX_LINE {
            return Err(error(self.line, format!("line exceeds {MAX_LINE} bytes")));
        }
        let number = self.line;
        self.line += 1;
        self.rest = rest;
        Ok(Some((number, line.trim())))
    }

    fn section(&mut self) -> Result<Option<(&'a str, Tokens<'a>)>> {
        while let Some((line, header)) = self.next()? {
            if header.is_empty() {
                continue;
            }
            let name = header.strip_prefix('$').filter(|n| !n.is_empty() && n.len() <= 128 && !n.starts_with("End"));
            let Some(name) = name else { return Err(error(line, "expected a section header")) };
            let body = self.rest;
            let first = self.line;
            let end = format!("$End{name}");
            loop {
                let before = self.rest;
                let Some((_, text)) = self.next()? else { return Err(error(line, format!("missing {end}"))) };
                if text == end {
                    let bytes = body.len() - before.len();
                    let text = body.get(..bytes).ok_or_else(|| error(line, "invalid section boundary"))?;
                    return Ok(Some((name, Tokens::new(text, first))));
                }
            }
        }
        Ok(None)
    }
}

/// A bounded, streaming token cursor. Counting once per section avoids rescanning the suffix
/// for each hostile count. No allocation is proportional to an unchecked count.
struct Tokens<'a> {
    rest: &'a str,
    line: usize,
    remaining: usize,
}

impl<'a> Tokens<'a> {
    fn new(rest: &'a str, line: usize) -> Self {
        Self { rest, line, remaining: rest.split_whitespace().count() }
    }

    fn next(&mut self) -> Result<&'a str> {
        let trimmed = self.rest.trim_start();
        let whitespace = self.rest.get(..self.rest.len() - trimmed.len()).unwrap_or_default();
        self.line += whitespace.bytes().filter(|&b| b == b'\n').count();
        if trimmed.is_empty() {
            return Err(error(self.line, "truncated section"));
        }
        let end = trimmed.find(char::is_whitespace).unwrap_or(trimmed.len());
        let token = trimmed.get(..end).ok_or_else(|| error(self.line, "invalid token"))?;
        if token.len() > 256 {
            return Err(error(self.line, "token exceeds 256 bytes"));
        }
        self.rest = trimmed.get(end..).unwrap_or_default();
        self.remaining = self.remaining.saturating_sub(1);
        Ok(token)
    }

    fn number<T: std::str::FromStr>(&mut self, what: &str) -> Result<T> {
        let token = self.next()?;
        token.parse().map_err(|_| error(self.line, format!("invalid {what} `{token}`")))
    }

    fn count(&mut self, what: &str, cap: usize, tokens: usize) -> Result<usize> {
        let n: usize = self.number(what)?;
        if n > cap || n > self.remaining / tokens {
            return Err(error(self.line, format!("{what} {n} exceeds the limit or remaining input")));
        }
        Ok(n)
    }

    fn tag(&mut self, what: &str) -> Result<u32> {
        let tag: u32 = self.number(what)?;
        if !(1..=MAX_LABEL).contains(&tag) {
            return Err(error(self.line, format!("{what} {tag} is outside 1..={MAX_LABEL}")));
        }
        Ok(tag)
    }

    fn dimension(&mut self) -> Result<u8> {
        let d: u8 = self.number("entity dimension")?;
        if d > 3 {
            return Err(error(self.line, "entity dimension must be 0..=3"));
        }
        Ok(d)
    }

    fn coordinate(&mut self) -> Result<f64> {
        let n: f64 = self.number("coordinate")?;
        if !(n.abs() <= MAX_COORD) {
            return Err(error(self.line, format!("coordinate must be finite and within ±{MAX_COORD}")));
        }
        Ok(n)
    }

    fn name(&mut self) -> Result<&'a str> {
        let rest = self.rest.trim_start();
        self.line += self.rest.get(..self.rest.len() - rest.len()).unwrap_or_default().bytes().filter(|&b| b == b'\n').count();
        let quoted = rest.strip_prefix('"').ok_or_else(|| error(self.line, "expected a quoted physical name"))?;
        let end = quoted.find('"').filter(|&n| n <= 127).ok_or_else(|| error(self.line, "invalid physical name (at most 127 bytes)"))?;
        let name = quoted.get(..end).ok_or_else(|| error(self.line, "invalid physical name"))?;
        fieldcraft_model::check_name("physical group", name).map_err(|e| error(self.line, format!("physical name `{name}`: {e}")))?;
        self.rest = quoted.get(end + 1..).unwrap_or_default();
        // A valid model name contains no whitespace: its quotes form one token.
        self.remaining = self.remaining.saturating_sub(1);
        Ok(name)
    }

    fn done(&self) -> Result<()> {
        if self.rest.trim().is_empty() { Ok(()) } else { Err(error(self.line, "unexpected trailing section data")) }
    }
}

#[derive(Debug)]
struct Record {
    entity: (u8, u32),
    kind: Option<ElementKind>,
    nodes: Vec<u32>,
    line: usize,
}

struct Reader {
    model: Model,
    records: BTreeMap<u32, Record>,
    entities: BTreeMap<(u8, u32), Vec<u32>>,
    /// Validate after all sections: an Entities table may follow or extend earlier data.
    entity_refs: BTreeMap<(u8, u32), usize>,
    has_entities: bool,
    names: BTreeMap<(u8, u32), String>,
    limits: Limits,
    memberships: usize,
}

impl Reader {
    fn entities(&mut self, t: &mut Tokens<'_>) -> Result<()> {
        let mut counts = [0; 4];
        for count in &mut counts {
            *count = t.number::<usize>("entity count")?;
        }
        let total = counts.iter().try_fold(0usize, |sum, &n| sum.checked_add(n)).ok_or_else(|| error(t.line, "entity count overflow"))?;
        if total > self.limits.nodes.saturating_add(self.limits.elements).saturating_sub(self.entities.len()) || total > t.remaining / 5 {
            return Err(error(t.line, "entity counts exceed the limit or remaining input"));
        }
        for (d, &n) in counts.iter().enumerate() {
            for _ in 0..n {
                let tag = t.tag("entity tag")?;
                for _ in 0..if d == 0 { 3 } else { 6 } {
                    t.coordinate()?;
                }
                let n = t.count("physical tag count", MAX_DEFINITIONS, 1)?;
                self.memberships = budget(self.memberships, n, self.limits.members, t.line)?;
                let mut physical = Vec::new();
                for _ in 0..n {
                    physical.push(t.tag("physical tag")?);
                }
                physical.sort_unstable();
                physical.dedup();
                if d > 0 {
                    let n = t.count("bounding entity count", self.limits.elements, 1)?;
                    for _ in 0..n {
                        let signed: i64 = t.number("bounding entity tag")?;
                        if !(1..=u64::from(MAX_LABEL)).contains(&signed.unsigned_abs()) {
                            return Err(error(t.line, "invalid bounding entity tag"));
                        }
                    }
                }
                if self.entities.insert((d as u8, tag), physical).is_some() {
                    return Err(error(t.line, format!("duplicate entity tag {tag} in dimension {d}")));
                }
            }
        }
        t.done()
    }

    fn names(&mut self, t: &mut Tokens<'_>) -> Result<()> {
        let count = t.count("physical name count", MAX_DEFINITIONS.saturating_sub(self.names.len()), 3)?;
        for _ in 0..count {
            let d = t.dimension()?;
            let tag = t.tag("physical tag")?;
            let name = t.name()?.to_string();
            if self.names.insert((d, tag), name).is_some() {
                return Err(error(t.line, format!("duplicate physical name for ({d}, {tag})")));
            }
        }
        t.done()
    }

    fn nodes(&mut self, t: &mut Tokens<'_>) -> Result<()> {
        let blocks = t.count("node block count", self.limits.nodes, 4)?;
        let count = t.count("node count", self.limits.nodes.saturating_sub(self.model.nodes().len()), 4)?;
        let min: u32 = t.number("minimum node tag")?;
        let max: u32 = t.number("maximum node tag")?;
        let (mut found, mut low, mut high) = (0, u32::MAX, 0);
        for _ in 0..blocks {
            let dim = t.dimension()?;
            let entity = t.tag("node entity tag")?;
            self.note_entity(dim, entity, t.line)?;
            let parametric: u8 = t.number("parametric flag")?;
            if parametric > 1 {
                return Err(error(t.line, "parametric flag must be 0 or 1"));
            }
            let n = t.count("nodes in block", count.saturating_sub(found), 4 + usize::from(parametric * dim))?;
            let mut tags = Vec::new();
            for _ in 0..n {
                let tag = t.tag("node tag")?;
                low = low.min(tag);
                high = high.max(tag);
                tags.push(tag);
            }
            for tag in tags {
                let xyz = [t.coordinate()?, t.coordinate()?, t.coordinate()?];
                for _ in 0..parametric * dim {
                    t.coordinate()?;
                }
                if self.model.node(tag).is_some() {
                    return Err(error(t.line, format!("duplicate node tag {tag}")));
                }
                self.model.set_node(tag, xyz).map_err(|e| error(t.line, e))?;
            }
            found += n;
        }
        totals(t, "nodes", count, found, min, max, low, high)
    }

    fn note_entity(&mut self, d: u8, tag: u32, line: usize) -> Result<()> {
        if !self.entity_refs.contains_key(&(d, tag)) && self.entity_refs.len() >= self.limits.nodes.saturating_add(self.limits.elements) {
            return Err(error(line, "entity reference count exceeds limit"));
        }
        self.entity_refs.entry((d, tag)).or_insert(line);
        Ok(())
    }

    fn elements(&mut self, t: &mut Tokens<'_>) -> Result<()> {
        let blocks = t.count("element block count", self.limits.elements, 4)?;
        let count = t.count("element count", self.limits.elements.saturating_sub(self.records.len()), 2)?;
        let min: u32 = t.number("minimum element tag")?;
        let max: u32 = t.number("maximum element tag")?;
        let (mut found, mut low, mut high) = (0, u32::MAX, 0);
        for _ in 0..blocks {
            let dim = t.dimension()?;
            let entity = t.tag("element entity tag")?;
            self.note_entity(dim, entity, t.line)?;
            let code: u32 = t.number("element type")?;
            let kind = kind(code).map_err(|name| error(t.line, format!("unsupported element type {code} ({name})")))?;
            if kind.map_or(0, ElementKind::dimension) != dim {
                return Err(error(t.line, format!("element type {code} does not have entity dimension {dim}")));
            }
            let nodes = kind.map_or(1, ElementKind::node_count);
            let n = t.count("elements in block", count.saturating_sub(found), 1 + nodes)?;
            for _ in 0..n {
                let tag = t.tag("element tag")?;
                low = low.min(tag);
                high = high.max(tag);
                let mut connectivity = Vec::new();
                for _ in 0..nodes {
                    let node = t.tag("element node tag")?;
                    if self.model.node(node).is_none() {
                        return Err(error(t.line, format!("element {tag}: missing node {node}")));
                    }
                    if connectivity.contains(&node) {
                        return Err(error(t.line, format!("element {tag}: duplicate node {node}")));
                    }
                    connectivity.push(node);
                }
                // Derived from the manual's reference-edge figures: FC edge order, then Gmsh
                // midpoint index. Corners retain the same positive reference orientation.
                let permutation: &[usize] = match code {
                    11 => &[0, 1, 2, 3, 4, 5, 6, 7, 9, 8],
                    17 => &[0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 16, 18, 19, 17, 10, 12, 14, 15],
                    _ => &[],
                };
                if !permutation.is_empty() {
                    connectivity = permutation
                        .iter()
                        .map(|&i| connectivity.get(i).copied().ok_or_else(|| error(t.line, "short connectivity")))
                        .collect::<Result<_>>()?;
                }
                let record = Record { entity: (dim, entity), kind, nodes: connectivity, line: t.line };
                if self.records.insert(tag, record).is_some() {
                    return Err(error(t.line, format!("duplicate element tag {tag}")));
                }
            }
            found += n;
        }
        totals(t, "elements", count, found, min, max, low, high)
    }

    fn finish(mut self, unmatched: Unmatched, skipped_sections: Vec<String>) -> Result<Imported> {
        if self.has_entities {
            for (&(dim, tag), &line) in &self.entity_refs {
                if !self.entities.contains_key(&(dim, tag)) {
                    return Err(error(line, format!("missing entity ({dim}, {tag})")));
                }
            }
        }
        let dimension = self.records.values().map(|r| r.entity.0).max().unwrap_or(0);
        if dimension == 0 {
            return Err(error(1, "mesh needs line, surface or volume elements"));
        }
        // Index only the requested corner tuples, rather than materializing every volume side.
        let mut faces: BTreeMap<Vec<u32>, Vec<ElementFace>> =
            self.records.values().filter(|r| r.entity.0 == 2 && dimension == 3).map(|r| (corners(r), Vec::new())).collect();
        let mut edges: BTreeMap<Vec<u32>, Vec<ElementEdge>> =
            self.records.values().filter(|r| r.entity.0 == 1 && dimension > 1).map(|r| (corners(r), Vec::new())).collect();
        let mut side_matches = 0;
        for (&label, r) in &self.records {
            if r.entity.0 != dimension || dimension != 3 {
                continue;
            }
            if let Some(kind) = r.kind {
                for (i, face) in kind.faces().iter().enumerate() {
                    // Quadratic faces list their corners first, then the same number of mids.
                    let n = if kind.is_quadratic() { face.len() / 2 } else { face.len() };
                    let key = side_key(&r.nodes, face.get(..n).unwrap_or_default(), r.line)?;
                    if let Some(matches) = faces.get_mut(&key) {
                        side_matches = budget(side_matches, 1, self.limits.members, r.line)?;
                        matches.push(ElementFace { element: label, face: (i + 1) as u8 });
                    }
                }
            }
        }
        let mut kept = BTreeSet::new();
        for (&label, r) in &self.records {
            if r.entity.0 == dimension {
                kept.insert(label);
            } else if r.entity.0 == 2 && faces.get(&corners(r)).is_none_or(Vec::is_empty) {
                if unmatched != Unmatched::Shell {
                    return Err(error(
                        r.line,
                        format!("surface element {label} matches no volume face; use unmatched: shell only to keep an intentional shell"),
                    ));
                }
                kept.insert(label);
            }
        }
        for &label in &kept {
            let r = self.records.get(&label).ok_or_else(|| error(1, "missing retained element"))?;
            if let Some(kind) = r.kind {
                if kind.dimension() >= 2 {
                    for (i, edge) in kind.edges().iter().enumerate() {
                        let key = side_key(&r.nodes, edge.get(..2).unwrap_or_default(), r.line)?;
                        if let Some(matches) = edges.get_mut(&key) {
                            side_matches = budget(side_matches, 1, self.limits.members, r.line)?;
                            matches.push(ElementEdge { element: label, edge: (i + 1) as u8 });
                        }
                    }
                }
                self.model.set_element(label, kind, r.nodes.clone()).map_err(|e| error(r.line, e))?;
            }
        }
        let mut groups: BTreeMap<(u8, u32), Set> = BTreeMap::new();
        for (&label, r) in &self.records {
            let physical = self.entities.get(&r.entity).map(Vec::as_slice).unwrap_or_default();
            if r.kind.is_none() && physical.is_empty() {
                return Err(error(r.line, format!("point element {label} needs a 0-dimensional physical group")));
            }
            if r.entity.0 == 1 && dimension > 1 && edges.get(&corners(r)).is_none_or(Vec::is_empty) {
                return Err(error(r.line, format!("curve element {label} matches no element edge")));
            }
            for &tag in physical {
                if !groups.contains_key(&(r.entity.0, tag)) && groups.len() >= MAX_DEFINITIONS {
                    return Err(error(r.line, "too many physical groups"));
                }
                let set = groups.entry((r.entity.0, tag)).or_default();
                if kept.contains(&label) {
                    self.memberships = budget(self.memberships, 1, self.limits.members, r.line)?;
                    set.elements.push(label);
                } else if r.entity.0 == 2 {
                    let matches = faces.get(&corners(r)).map(Vec::as_slice).unwrap_or_default();
                    self.memberships = budget(self.memberships, matches.len(), self.limits.members, r.line)?;
                    set.faces.extend_from_slice(matches);
                } else if r.entity.0 == 1 {
                    let matches = edges.get(&corners(r)).map(Vec::as_slice).unwrap_or_default();
                    self.memberships = budget(self.memberships, matches.len(), self.limits.members, r.line)?;
                    set.edges.extend_from_slice(matches);
                } else {
                    self.memberships = budget(self.memberships, r.nodes.len(), self.limits.members, r.line)?;
                    set.nodes.extend_from_slice(&r.nodes);
                }
            }
        }
        for ((dim, tag), set) in groups {
            let name = self.names.get(&(dim, tag)).cloned().unwrap_or_else(|| format!("physical_{dim}_{tag}"));
            self.model.add_set(&name, set).map_err(|e| error(1, format!("physical group ({dim}, {tag}) `{name}`: {e}")))?;
        }
        Ok(Imported { model: self.model, skipped_sections })
    }
}

fn corners(r: &Record) -> Vec<u32> {
    let n = r.kind.map_or(1, ElementKind::corner_count);
    let mut key = r.nodes.get(..n).unwrap_or_default().to_vec();
    key.sort_unstable();
    key
}

fn side_key(nodes: &[u32], side: &[usize], line: usize) -> Result<Vec<u32>> {
    let mut key = side.iter().map(|&i| nodes.get(i).copied().ok_or_else(|| error(line, "short element side"))).collect::<Result<Vec<_>>>()?;
    key.sort_unstable();
    Ok(key)
}

fn budget(used: usize, n: usize, cap: usize, line: usize) -> Result<usize> {
    used.checked_add(n).filter(|&sum| sum <= cap).ok_or_else(|| error(line, "physical membership/side count exceeds limit"))
}

#[allow(clippy::too_many_arguments)]
fn totals(t: &Tokens<'_>, what: &str, count: usize, found: usize, min: u32, max: u32, low: u32, high: u32) -> Result<()> {
    if found != count || (count > 0 && (min != low || max != high)) || (count == 0 && (min != 0 || max != 0)) {
        return Err(error(t.line, format!("{what} count or min/max tags do not match the records")));
    }
    t.done()
}

fn kind(code: u32) -> std::result::Result<Option<ElementKind>, &'static str> {
    Ok(Some(match code {
        1 => ElementKind::Rod2,
        8 => ElementKind::Rod3,
        2 => ElementKind::Tri3,
        9 => ElementKind::Tri6,
        3 => ElementKind::Quad4,
        16 => ElementKind::Quad8,
        4 => ElementKind::Tet4,
        11 => ElementKind::Tet10,
        6 => ElementKind::Wedge6,
        5 => ElementKind::Hex8,
        17 => ElementKind::Hex20,
        15 => return Ok(None),
        7 => return Err("5-node pyramid"),
        10 => return Err("9-node quadrangle"),
        12 => return Err("27-node hexahedron"),
        13 => return Err("18-node prism"),
        14 => return Err("14-node pyramid"),
        18 => return Err("15-node prism"),
        19 => return Err("13-node pyramid"),
        20 => return Err("9-node third-order incomplete triangle"),
        21 => return Err("10-node third-order triangle"),
        22 => return Err("12-node fourth-order incomplete triangle"),
        23 => return Err("15-node fourth-order triangle"),
        24 => return Err("15-node fifth-order incomplete triangle"),
        25 => return Err("21-node fifth-order triangle"),
        26 => return Err("4-node third-order line"),
        27 => return Err("5-node fourth-order line"),
        28 => return Err("6-node fifth-order line"),
        29 => return Err("20-node third-order tetrahedron"),
        30 => return Err("35-node fourth-order tetrahedron"),
        31 => return Err("56-node fifth-order tetrahedron"),
        92 => return Err("64-node third-order hexahedron"),
        93 => return Err("125-node fourth-order hexahedron"),
        _ => return Err("unknown or higher-order element"),
    }))
}

/// Read only ASCII MSH 4.1. Every bound is intersected with the model's own caps. All validation
/// happens on a new model: refusal cannot mutate the caller's model. Unsupported metadata is
/// reported by section name; no periodic or partition semantics are guessed.
pub fn read_msh(text: &str, limits: Limits, unmatched: Unmatched) -> Result<Imported> {
    if text.len() as u64 > MAX_BYTES {
        return Err(Error::Limit("MSH text exceeds 256 MiB".into()));
    }
    let limits =
        Limits { nodes: limits.nodes.min(MAX_NODES), elements: limits.elements.min(MAX_ELEMENTS), members: limits.members.min(MAX_ELEMENTS) };
    let mut lines = Lines { rest: text, line: 1 };
    let Some(("MeshFormat", mut format)) = lines.section()? else { return Err(error(1, "missing $MeshFormat")) };
    let version = format.next()?;
    if version != "4.1" {
        return Err(error(format.line, format!("MSH version {version} not supported; expected 4.1")));
    }
    let mode: u8 = format.number("file type")?;
    if mode == 1 {
        return Err(error(format.line, "binary MSH not supported yet"));
    }
    if mode != 0 {
        return Err(error(format.line, "MSH file type must be 0 (ASCII)"));
    }
    let size: u8 = format.number("data size")?;
    if ![4, 8].contains(&size) {
        return Err(error(format.line, "MSH data size must be 4 or 8"));
    }
    format.done()?;
    let mut reader = Reader {
        model: Model::default(),
        records: BTreeMap::new(),
        entities: BTreeMap::new(),
        entity_refs: BTreeMap::new(),
        has_entities: false,
        names: BTreeMap::new(),
        limits,
        memberships: 0,
    };
    let mut seen = BTreeSet::new();
    let mut skipped = BTreeSet::new();
    while let Some((name, mut tokens)) = lines.section()? {
        seen.insert(name);
        match name {
            "MeshFormat" => return Err(error(tokens.line, "repeated $MeshFormat")),
            "Entities" => {
                reader.has_entities = true;
                reader.entities(&mut tokens)?;
            }
            "PhysicalNames" => reader.names(&mut tokens)?,
            "Nodes" => reader.nodes(&mut tokens)?,
            "Elements" => reader.elements(&mut tokens)?,
            "PartitionedEntities" => return Err(error(tokens.line, "partitioned MSH not supported yet")),
            _ => {
                if skipped.len() >= MAX_DEFINITIONS {
                    return Err(error(tokens.line, "too many ignored section names"));
                }
                skipped.insert(name.to_string());
            }
        }
    }
    if !seen.contains("Nodes") || !seen.contains("Elements") {
        return Err(error(1, "missing $Nodes or $Elements"));
    }
    reader.finish(unmatched, skipped.into_iter().collect())
}

#[cfg(test)]
mod tests;
