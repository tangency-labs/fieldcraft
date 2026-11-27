//! External result grids and time collections, from Kitware's
//! [VTK XML File Formats](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html).
//! A piece's point index is its identity. Coincident points can carry different field values.
//! In a PVD collection every frame owns its grid, including the coordinates at that time.

mod array;
#[cfg(test)]
mod tests;

use std::borrow::Cow;
use std::collections::BTreeMap;

use fieldcraft_model::units::{Quantity, Unit};
use quick_xml::{Reader, events::Event};
use roxmltree::{Document, Node, ParsingOptions};

use crate::external::{Cell, Field, Frame, Grid, MAX_CELLS, MAX_COMPONENTS, MAX_FIELDS, MAX_FRAME_BYTES, MAX_FRAMES, MAX_POINTS, cell_arity};
use crate::{Error, Result};

fn bad(message: impl std::fmt::Display) -> Error {
    Error::Invalid(format!("VTK XML: {message}"))
}

fn utf8(bytes: &[u8]) -> Result<&str> {
    std::str::from_utf8(bytes).map_err(|_| bad("invalid UTF-8 in XML metadata"))
}

fn document(xml: &str) -> Result<Document<'_>> {
    // Both the metadata tree and the result arrays have independent allocation limits.
    let options = ParsingOptions { nodes_limit: 64_000, ..Default::default() };
    Document::parse_with_options(xml, options).map_err(bad)
}

fn root<'a>(doc: &'a Document<'_>, kind: &str) -> Result<Node<'a, 'a>> {
    let root = doc.root_element();
    if !root.has_tag_name("VTKFile") || root.attribute("type") != Some(kind) {
        return Err(bad(format!("expected VTKFile of type {kind}")));
    }
    Ok(root)
}

fn one<'a, 'i>(parent: Node<'a, 'i>, tag: &str) -> Result<Node<'a, 'i>> {
    let mut items = parent.children().filter(|child| child.has_tag_name(tag));
    let node = items.next().ok_or_else(|| bad(format!("missing {tag}")))?;
    if items.next().is_some() {
        return Err(bad(format!("multiple {tag} elements are not supported")));
    }
    Ok(node)
}

fn attr<'a>(node: Node<'a, '_>, name: &str) -> Result<&'a str> {
    node.attribute(name).ok_or_else(|| bad(format!("missing {name} on {}", node.tag_name().name())))
}

fn count(node: Node<'_, '_>, name: &str, limit: usize) -> Result<usize> {
    let n = attr(node, name)?.parse::<usize>().map_err(bad)?;
    if n > limit {
        return Err(Error::Limit(format!("{name} exceeds {limit}")));
    }
    Ok(n)
}

/// Read XML tokens as bytes until the real AppendedData opening tag. Raw payload bytes
/// must never enter an XML parser. Comments, CDATA, and quoted attributes are distinct
/// token kinds, even if their text looks like a tag.
fn appended_start(bytes: &[u8]) -> Result<Option<usize>> {
    let mut reader = Reader::from_reader(bytes);
    let mut depth = 0usize;
    for _ in 0..128_000 {
        match reader.read_event().map_err(bad)? {
            Event::Start(e) if depth == 1 && e.name().as_ref() == "AppendedData" => {
                return Ok(Some(usize::try_from(reader.buffer_position()).map_err(bad)?));
            }
            Event::Start(_) => depth += 1,
            Event::End(_) => depth = depth.saturating_sub(1),
            Event::Eof => return Ok(None),
            _ => {}
        }
    }
    Err(Error::Limit("too many XML tokens before VTK appended data".into()))
}

/// AppendedData is the last child of VTKFile. Validate its trailer so even a comment
/// containing a literal closing tag after the payload cannot move an array boundary.
fn appended_end(body: &[u8]) -> Result<usize> {
    let close = b"</AppendedData";
    let mut allowance = body.len().saturating_mul(2);
    for (i, _) in body.windows(close.len()).enumerate().rev().filter(|(_, w)| *w == close) {
        let tail = body.get(i..).ok_or_else(|| bad("invalid AppendedData trailer"))?;
        let mut reader = Reader::from_reader(tail);
        reader.config_mut().allow_unmatched_ends = true;
        let mut ends = ["AppendedData", "VTKFile"].into_iter();
        loop {
            match reader.read_event() {
                Ok(Event::End(e)) if ends.next() == Some(e.name().as_ref()) => {}
                Ok(Event::Text(t)) if t.as_ref().bytes().all(|b| b.is_ascii_whitespace()) => {}
                Ok(Event::Comment(_) | Event::PI(_)) => {}
                Ok(Event::Eof) if ends.next().is_none() => return Ok(i),
                _ => break,
            }
        }
        allowance = allowance
            .checked_sub(usize::try_from(reader.buffer_position()).map_err(bad)?)
            .ok_or_else(|| Error::Limit("invalid or excessively complex VTU trailer".into()))?;
    }
    Err(bad("missing AppendedData closing tag"))
}

fn metadata(bytes: &[u8]) -> Result<(Cow<'_, str>, Option<&[u8]>)> {
    if bytes.len() > MAX_FRAME_BYTES {
        return Err(Error::Limit("VTU file too large".into()));
    }
    if !bytes.windows(b"<AppendedData".len()).any(|w| w == b"<AppendedData") {
        return Ok((Cow::Borrowed(utf8(bytes)?), None));
    }
    let Some(begin) = appended_start(bytes)? else {
        return Ok((Cow::Borrowed(utf8(bytes)?), None));
    };
    let remaining = bytes.get(begin..).ok_or_else(|| bad("invalid AppendedData bounds"))?;
    let end = begin + appended_end(remaining)?;
    let body = bytes.get(begin..end).ok_or_else(|| bad("invalid AppendedData bounds"))?;
    let sentinel = body.iter().position(|b| !b.is_ascii_whitespace()).ok_or_else(|| bad("missing AppendedData sentinel"))?;
    if body.get(sentinel) != Some(&b'_') {
        return Err(bad("AppendedData must start with '_'"));
    }
    let mut xml = utf8(bytes.get(..begin).ok_or_else(|| bad("invalid XML prefix"))?)?.to_owned();
    xml.push_str(utf8(bytes.get(end..).ok_or_else(|| bad("invalid XML suffix"))?)?);
    Ok((Cow::Owned(xml), body.get(sentinel + 1..)))
}

fn components(node: Node<'_, '_>) -> Result<usize> {
    let n = node.attribute("NumberOfComponents").unwrap_or("1").parse::<usize>().map_err(bad)?;
    if !(1..=MAX_COMPONENTS).contains(&n) {
        return Err(bad(format!("unsupported component count {n}")));
    }
    Ok(n)
}

fn cell_array<'a, 'i>(cells: Node<'a, 'i>, name: &str) -> Result<Node<'a, 'i>> {
    let mut arrays = cells.children().filter(|n| n.has_tag_name("DataArray") && n.attribute("Name") == Some(name));
    let node = arrays.next().ok_or_else(|| bad(format!("missing Cells array {name}")))?;
    if arrays.next().is_some() || components(node)? != 1 {
        return Err(bad(format!("expected one scalar Cells array named {name}")));
    }
    Ok(node)
}

/// The Lagrange vertices follow the corresponding linear cell's order, then the
/// edge-interior samples follow VTK's edge order. Quadratic curves, triangles, and
/// tetrahedra have exactly one point per edge and no face or volume points.
/// See Corona and Thompson, [Modeling Arbitrary-order Lagrange Finite Elements in the
/// Visualization Toolkit](https://www.kitware.com/modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/).
fn canonical_cell(vtk_type: u32, arity: usize) -> Result<u8> {
    let kind = match (vtk_type, arity) {
        (68, 2) => 3,
        (68, 3) => 21,
        (69, 3) => 5,
        (69, 6) => 22,
        (70, 4) => 9,
        (71, 4) => 10,
        (71, 10) => 24,
        (72, 8) => 12,
        (73, 6) => 13,
        (68..=73, _) => return Err(bad(format!("unsupported Lagrange cell type {vtk_type} with {arity} points"))),
        _ => {
            let kind = u8::try_from(vtk_type)
                .ok()
                .filter(|k| cell_arity(*k).is_some())
                .ok_or_else(|| bad(format!("unsupported VTK cell type {vtk_type}")))?;
            if cell_arity(kind) != Some(arity) {
                return Err(bad(format!("cell offsets give VTK cell type {kind} the wrong number of points")));
            }
            kind
        }
    };
    Ok(kind)
}

fn geometry(decoder: &array::Decoder<'_>, piece: Node<'_, '_>, np: usize, nc: usize, budget: &mut usize) -> Result<Grid> {
    crate::external::charge_storage(budget, crate::external::FRAME_BYTES)?;
    crate::external::charge_storage(budget, np.saturating_mul(24))?;
    crate::external::charge_storage(budget, nc.saturating_mul(crate::external::CELL_BYTES))?;
    let pos = one(one(piece, "Points")?, "DataArray")?;
    if components(pos)? != 3 {
        return Err(bad("Points must have 3 components"));
    }
    let coords = decoder.coordinates(pos, 3 * np)?;
    let points = coords.as_chunks::<3>().0.to_vec();
    drop(coords);
    let source = one(piece, "Cells")?;
    let types = decoder.indices(cell_array(source, "types")?, nc, 0)?;
    let offsets = decoder.indices(cell_array(source, "offsets")?, nc, types.capacity())?;
    let mut cells = Vec::with_capacity(nc);
    let mut entries = 0usize;
    for (number, offset) in types.into_iter().zip(offsets) {
        let end = offset as usize;
        let arity = end.checked_sub(entries).ok_or_else(|| bad("cell offsets must increase"))?;
        cells.push(Cell { vtk_type: canonical_cell(number, arity)?, points: Vec::new() });
        entries = end;
    }
    crate::external::charge_storage(budget, entries.saturating_mul(4))?;
    let indices = decoder.indices(cell_array(source, "connectivity")?, entries, 0)?;
    if indices.iter().any(|&index| index as usize >= np) {
        return Err(bad("cell references a missing point"));
    }
    let mut remaining = indices.as_slice();
    for cell in &mut cells {
        let arity = cell_arity(cell.vtk_type).ok_or_else(|| bad("invalid cell type"))?;
        let (points, tail) = remaining.split_at_checked(arity).ok_or_else(|| bad("invalid cell connectivity"))?;
        cell.points = points.to_vec();
        remaining = tail;
    }
    Ok(Grid { points, cells })
}

fn fields(
    decoder: &array::Decoder<'_>,
    piece: Node<'_, '_>,
    tag: &str,
    tuples: usize,
    remaining: usize,
    budget: &mut usize,
) -> Result<BTreeMap<String, Field>> {
    let mut result = BTreeMap::new();
    if !piece.children().any(|n| n.has_tag_name(tag)) {
        return Ok(result);
    }
    for node in one(piece, tag)?.children().filter(Node::is_element) {
        if !node.has_tag_name("DataArray") {
            return Err(bad(format!("unsupported {} in {tag}", node.tag_name().name())));
        }
        let name = attr(node, "Name")?;
        if !crate::external::valid_name(name) || result.contains_key(name) {
            return Err(bad(format!("empty, oversized, or duplicate array name in {tag}")));
        }
        if result.len() >= remaining {
            return Err(Error::Limit("too many VTU fields".into()));
        }
        let n = components(node)?;
        let entries = tuples.checked_mul(n).ok_or_else(|| bad("field is too large"))?;
        crate::external::charge_storage(budget, crate::external::FIELD_BYTES.saturating_add(name.len()).saturating_add(entries.saturating_mul(8)))?;
        let values = decoder.field(node, entries)?;
        result.insert(name.to_owned(), Field { components: n as u8, values });
    }
    Ok(result)
}

/// Read one serial unstructured grid. Time is supplied by the collection or its caller.
pub fn read(bytes: &[u8], time: f64) -> Result<Frame> {
    read_with_budget(bytes, time, MAX_FRAME_BYTES)
}

/// Read one frame within the storage still available to its series.
pub fn read_with_budget(bytes: &[u8], time: f64, remaining: usize) -> Result<Frame> {
    read_with_unit(bytes, time, remaining, None)
}

/// Read coordinates in the supplied length unit. VTU files carry no length unit.
/// Source field samples keep their values until explicitly interpreted as displacement.
pub fn read_with_unit(bytes: &[u8], time: f64, remaining: usize, unit: Option<&Unit>) -> Result<Frame> {
    if !time.is_finite() {
        return Err(bad("frame time must be finite"));
    }
    let (xml, payload) = metadata(bytes)?;
    let doc = document(&xml)?;
    let vtk = root(&doc, "UnstructuredGrid")?;
    decode_piece(vtk, payload, time, remaining, MAX_FRAME_BYTES, unit)
}

fn decode_piece(vtk: Node<'_, '_>, payload: Option<&[u8]>, time: f64, remaining: usize, frame_limit: usize, unit: Option<&Unit>) -> Result<Frame> {
    let mut budget = remaining.min(frame_limit);
    if let Some(unit) = unit {
        Quantity::Length.to_si(0.0, unit).map_err(bad)?;
        crate::external::charge_storage(&mut budget, unit.symbol().len())?;
    }
    let decoder = array::Decoder::new(vtk, payload, frame_limit)?;
    let piece = one(one(vtk, "UnstructuredGrid")?, "Piece")?;
    let np = count(piece, "NumberOfPoints", MAX_POINTS)?;
    let nc = count(piece, "NumberOfCells", MAX_CELLS)?;
    let mut grid = geometry(&decoder, piece, np, nc, &mut budget)?;
    if let Some(unit) = unit {
        for coordinate in grid.points.iter_mut().flatten() {
            *coordinate = Quantity::Length.to_si(*coordinate, unit).map_err(bad)?;
        }
    }
    let point_data = fields(&decoder, piece, "PointData", np, MAX_FIELDS, &mut budget)?;
    let cell_data = fields(&decoder, piece, "CellData", nc, MAX_FIELDS - point_data.len(), &mut budget)?;
    let frame = Frame { time, displacement_unit: unit.cloned(), grid, point_data, cell_data };
    frame.check()?;
    Ok(frame)
}

/// Indexes only hold metadata and paths; numeric arrays belong to the VTU piece.
pub const MAX_INDEX_BYTES: usize = 4 << 20;

fn index_document(bytes: &[u8]) -> Result<Document<'_>> {
    if bytes.len() > MAX_INDEX_BYTES {
        return Err(Error::Limit("VTK index file exceeds 4 MiB".into()));
    }
    document(utf8(bytes)?)
}

fn reference(node: Node<'_, '_>, attribute: &str, extensions: &[&str]) -> Result<String> {
    let path = attr(node, attribute)?.replace('\\', "/");
    let permitted = std::path::Path::new(&path).extension().is_some_and(|e| extensions.iter().any(|ext| e.eq_ignore_ascii_case(ext)));
    if !permitted
        || path.len() > crate::MAX_PATH_LEN
        || path.starts_with('/')
        || path.contains(':')
        || path.chars().any(char::is_control)
        || path.split('/').any(|part| part == "..")
    {
        return Err(bad(format!("{attribute} must name a relative {} file within the index directory", extensions.join(" or "))));
    }
    Ok(path)
}

/// Resolve the one serial VTU piece in a PUnstructuredGrid index. The caller confines the
/// returned path to the index directory and reads the actual topology and arrays from the VTU.
pub fn single_piece(bytes: &[u8]) -> Result<String> {
    let doc = index_document(bytes)?;
    let grid = one(root(&doc, "PUnstructuredGrid")?, "PUnstructuredGrid")?;
    let mut pieces = grid.children().filter(|n| n.has_tag_name("Piece"));
    let piece = pieces.next().ok_or_else(|| bad("PVTU has no Piece"))?;
    if pieces.next().is_some() {
        return Err(bad("multi-piece PVTU is not supported; export one serial grid per time step"));
    }
    reference(piece, "Source", &["vtu"])
}

/// Read a serial PVD time collection. Members are grids or single-piece indexes.
pub fn collection(bytes: &[u8]) -> Result<Vec<(f64, String)>> {
    let doc = index_document(bytes)?;
    let vtk = root(&doc, "Collection")?;
    let mut frames = Vec::new();
    let mut previous_group = None;
    for entry in one(vtk, "Collection")?.children().filter(Node::is_element) {
        if !entry.has_tag_name("DataSet") || entry.attribute("part").unwrap_or("0") != "0" {
            return Err(bad("expected one VTU DataSet per time step; multipart PVD is not supported"));
        }
        let group = entry.attribute("group").unwrap_or("");
        if previous_group.replace(group).is_some_and(|previous| previous != group) {
            return Err(bad("multiple PVD groups are not supported"));
        }
        let time = attr(entry, "timestep")?.parse::<f64>().map_err(bad)?;
        if !time.is_finite() || frames.last().is_some_and(|(prev, _)| time <= *prev) {
            return Err(bad("PVD time steps must be finite and strictly increasing"));
        }
        let path = reference(entry, "file", &["vtu", "pvtu"])?;
        if frames.len() >= MAX_FRAMES {
            return Err(Error::Limit("too many PVD frames".into()));
        }
        frames.push((time, path));
    }
    if frames.is_empty() {
        return Err(bad("empty PVD collection"));
    }
    Ok(frames)
}
