//! VTK XML unstructured grids (`.vtu`) and their time collections (`.pvd`), for ParaView, VisIt
//! and other tools that read the VTK formats.
//!
//! The layout follows Kitware's public "VTK File Formats" document (the XML formats and the cell
//! types): `VTKFile type="UnstructuredGrid"` holding one `Piece` with `PointData`, `CellData`,
//! `Points` and `Cells` (`connectivity`, `offsets`, `types`), and an optional `FieldData` of
//! whole-dataset values such as a mode's frequency. A `.pvd` file is a `Collection` of
//! `DataSet` entries, each a file and its time.
//!
//! Points are the model's nodes in label order and cells its elements in label order; the point
//! array `node_label` and the cell array `element_label` keep the labels. Element node order is
//! the model's ([`ElementKind`]), which is VTK's for every kind. That holds for the wedge too:
//! VTK's wedge has the parametric coordinates (0,0,0), (1,0,0), (0,1,0) for points 0-2 and the
//! same at ζ = 1 for 3-5, so the base 0-1-2's right-hand normal points to 3-4-5, as a FieldCraft
//! wedge's 1-2-3 points to 4-5-6 (VTK 9 computes a positive volume for it and lists its outward
//! base face as 0-2-1).
//!
//! Data arrays are written as text (`ascii`, numbers in the shortest form that reads back to the
//! same `f64`) or inline base64 (`binary`, each array prefixed with its byte count as a `UInt64`,
//! little endian). The same input always gives the same bytes. Reading and writing work on bytes
//! (nothing here touches the file system), so this builds for the browser.

use std::fmt::Write;

use fieldcraft_model::{ElementKind, Model};

use crate::{Error, Result};

/// The VTK cell type id of an element kind ("VTK File Formats", cell types).
pub fn cell_type(kind: ElementKind) -> u8 {
    match kind {
        ElementKind::Rod2 | ElementKind::Beam2 => 3, // VTK_LINE
        ElementKind::Rod3 => 21,                     // VTK_QUADRATIC_EDGE
        ElementKind::Tri3 => 5,                      // VTK_TRIANGLE
        ElementKind::Tri6 => 22,                     // VTK_QUADRATIC_TRIANGLE
        ElementKind::Quad4 => 9,                     // VTK_QUAD
        ElementKind::Quad8 => 23,                    // VTK_QUADRATIC_QUAD
        ElementKind::Tet4 => 10,                     // VTK_TETRA
        ElementKind::Tet10 => 24,                    // VTK_QUADRATIC_TETRA
        ElementKind::Wedge6 => 13,                   // VTK_WEDGE
        ElementKind::Pyramid5 => 14,                 // VTK_PYRAMID
        ElementKind::Hex8 => 12,                     // VTK_HEXAHEDRON
        ElementKind::Hex20 => 25,                    // VTK_QUADRATIC_HEXAHEDRON
    }
}

/// Where a field's values are: one tuple per point (node) or per cell (element).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Location {
    Point,
    Cell,
}

/// A data array: `components` values per point or cell, in point (node label) or cell (element
/// label) order. Integer fields are written as `Int32`; non-finite values are written as 0.
#[derive(Clone, Debug, PartialEq)]
pub struct Field {
    pub name: String,
    pub location: Location,
    pub components: usize,
    pub values: Vec<f64>,
    pub integer: bool,
}

/// How data arrays are written.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Encoding {
    #[default]
    Ascii,
    /// Inline base64 (`format="binary"`).
    Base64,
}

/// Most components a field may have (a 3×3 tensor).
pub const MAX_COMPONENTS: usize = 9;

/// The model's node labels in point order and element labels in cell order.
pub fn points_and_cells(model: &Model) -> (Vec<u32>, Vec<u32>) {
    (model.nodes().keys().copied().collect(), model.elements().keys().copied().collect())
}

/// The model and `fields` as a `.vtu` document; `field_data` are whole-dataset values (a mode's
/// frequency, a frame's time).
pub fn write_vtu(model: &Model, fields: &[Field], field_data: &[(&str, f64)], encoding: Encoding) -> Result<String> {
    write_vtu_in(model, fields, field_data, encoding, &[])
}

/// [`write_vtu`] stating the unit of each array (`Points`, the fields, the field data) in a
/// comment after the XML declaration: `<!-- units: Points mm; displacement mm; … -->`. The
/// numbers are written as given: the caller expresses them in those units.
pub fn write_vtu_in(model: &Model, fields: &[Field], field_data: &[(&str, f64)], encoding: Encoding, units: &[(String, String)]) -> Result<String> {
    let (nodes, elements) = points_and_cells(model);
    for f in fields {
        let count = match f.location {
            Location::Point => nodes.len(),
            Location::Cell => elements.len(),
        };
        if f.components == 0 || f.components > MAX_COMPONENTS {
            return Err(Error::Invalid(format!("field `{}` has {} components (1 to {MAX_COMPONENTS})", f.name, f.components)));
        }
        if Some(f.values.len()) != count.checked_mul(f.components) {
            return Err(Error::Invalid(format!("field `{}` has {} values for {count} tuples of {}", f.name, f.values.len(), f.components)));
        }
    }
    let index = |label: u32| nodes.binary_search(&label).map_err(|_| Error::Invalid(format!("an element uses missing node {label}")));
    let mut connectivity: Vec<i64> = Vec::new();
    let mut offsets: Vec<i64> = Vec::with_capacity(elements.len());
    let mut types: Vec<u8> = Vec::with_capacity(elements.len());
    for e in model.elements().values() {
        for &label in &e.nodes {
            connectivity.push(index(label)? as i64);
        }
        offsets.push(connectivity.len() as i64);
        types.push(cell_type(e.kind));
    }

    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\"?>\n");
    if !units.is_empty() {
        // A comment holds no `--`; unit symbols never do, names are kept to letters and `_`.
        let list: Vec<String> = units
            .iter()
            .map(|(name, unit)| {
                let clean = |t: &str| t.chars().filter(|c| c.is_alphanumeric() || "_·²³⁴⁻¹/()°%".contains(*c)).collect::<String>();
                format!("{} {}", clean(name), clean(unit))
            })
            .collect();
        let _ = writeln!(out, "<!-- units: {} -->", list.join("; "));
    }
    out.push_str("<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    out.push_str("  <UnstructuredGrid>\n");
    if !field_data.is_empty() {
        out.push_str("    <FieldData>\n");
        for (name, v) in field_data {
            let array = Array { name, kind: Kind::Float64, components: 1, tuples: Some(1), values: Values::F64(std::slice::from_ref(v)) };
            array.write(&mut out, encoding, 6);
        }
        out.push_str("    </FieldData>\n");
    }
    let _ = writeln!(out, "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">", nodes.len(), elements.len());
    let node_labels: Vec<f64> = nodes.iter().map(|&l| f64::from(l)).collect();
    let element_labels: Vec<f64> = elements.iter().map(|&l| f64::from(l)).collect();
    for (tag, location, label_name, labels) in
        [("PointData", Location::Point, "node_label", &node_labels), ("CellData", Location::Cell, "element_label", &element_labels)]
    {
        let _ = writeln!(out, "      <{tag}>");
        Array { name: label_name, kind: Kind::Int32, components: 1, tuples: None, values: Values::F64(labels) }.write(&mut out, encoding, 8);
        for f in fields.iter().filter(|f| f.location == location) {
            let kind = if f.integer { Kind::Int32 } else { Kind::Float64 };
            Array { name: &f.name, kind, components: f.components, tuples: None, values: Values::F64(&f.values) }.write(&mut out, encoding, 8);
        }
        let _ = writeln!(out, "      </{tag}>");
    }
    out.push_str("      <Points>\n");
    let xyz: Vec<f64> = model.nodes().values().flat_map(|n| n.xyz).collect();
    Array { name: "Points", kind: Kind::Float64, components: 3, tuples: None, values: Values::F64(&xyz) }.write(&mut out, encoding, 8);
    out.push_str("      </Points>\n");
    out.push_str("      <Cells>\n");
    Array { name: "connectivity", kind: Kind::Int64, components: 1, tuples: None, values: Values::I64(&connectivity) }.write(&mut out, encoding, 8);
    Array { name: "offsets", kind: Kind::Int64, components: 1, tuples: None, values: Values::I64(&offsets) }.write(&mut out, encoding, 8);
    Array { name: "types", kind: Kind::UInt8, components: 1, tuples: None, values: Values::U8(&types) }.write(&mut out, encoding, 8);
    out.push_str("      </Cells>\n");
    out.push_str("    </Piece>\n");
    out.push_str("  </UnstructuredGrid>\n");
    out.push_str("</VTKFile>\n");
    Ok(out)
}

/// A `.pvd` collection of `(time, file name)` entries, the names relative to the collection.
pub fn pvd(entries: &[(f64, &str)]) -> String {
    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\"?>\n");
    out.push_str("<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\">\n");
    out.push_str("  <Collection>\n");
    for (time, file) in entries {
        let _ = writeln!(out, "    <DataSet timestep=\"{}\" part=\"0\" file=\"{}\"/>", number(*time), escape(file));
    }
    out.push_str("  </Collection>\n");
    out.push_str("</VTKFile>\n");
    out
}

#[derive(Clone, Copy)]
enum Kind {
    Float64,
    Int32,
    Int64,
    UInt8,
}

impl Kind {
    fn name(self) -> &'static str {
        match self {
            Kind::Float64 => "Float64",
            Kind::Int32 => "Int32",
            Kind::Int64 => "Int64",
            Kind::UInt8 => "UInt8",
        }
    }
}

#[derive(Clone, Copy)]
enum Values<'a> {
    F64(&'a [f64]),
    I64(&'a [i64]),
    U8(&'a [u8]),
}

struct Array<'a> {
    name: &'a str,
    kind: Kind,
    components: usize,
    /// `NumberOfTuples` (field data needs it).
    tuples: Option<usize>,
    values: Values<'a>,
}

impl Array<'_> {
    fn write(&self, out: &mut String, encoding: Encoding, indent: usize) {
        let pad = " ".repeat(indent);
        let format = match encoding {
            Encoding::Ascii => "ascii",
            Encoding::Base64 => "binary",
        };
        let _ = write!(out, "{pad}<DataArray type=\"{}\" Name=\"{}\"", self.kind.name(), escape(self.name));
        if self.components != 1 {
            let _ = write!(out, " NumberOfComponents=\"{}\"", self.components);
        }
        if let Some(n) = self.tuples {
            let _ = write!(out, " NumberOfTuples=\"{n}\"");
        }
        let _ = writeln!(out, " format=\"{format}\">");
        out.push_str(&pad);
        out.push_str("  ");
        match encoding {
            Encoding::Ascii => self.ascii(out),
            Encoding::Base64 => base64(&self.bytes(), out),
        }
        out.push('\n');
        let _ = writeln!(out, "{pad}</DataArray>");
    }

    /// The values as text: a line per tuple would be long for big models, so one line, space
    /// separated.
    fn ascii(&self, out: &mut String) {
        let mut first = true;
        let mut sep = |out: &mut String| {
            if !first {
                out.push(' ');
            }
            first = false;
        };
        match (self.values, self.kind) {
            (Values::F64(v), Kind::Float64) => v.iter().for_each(|x| {
                sep(out);
                out.push_str(&number(*x));
            }),
            (Values::F64(v), _) => v.iter().for_each(|x| {
                sep(out);
                let _ = write!(out, "{}", int32(*x));
            }),
            (Values::I64(v), _) => v.iter().for_each(|x| {
                sep(out);
                let _ = write!(out, "{x}");
            }),
            (Values::U8(v), _) => v.iter().for_each(|x| {
                sep(out);
                let _ = write!(out, "{x}");
            }),
        }
    }

    /// The raw little-endian bytes, behind their byte count (`UInt64`).
    fn bytes(&self) -> Vec<u8> {
        let mut data = Vec::new();
        match (self.values, self.kind) {
            (Values::F64(v), Kind::Float64) => v.iter().for_each(|x| data.extend_from_slice(&finite(*x).to_le_bytes())),
            (Values::F64(v), _) => v.iter().for_each(|x| data.extend_from_slice(&int32(*x).to_le_bytes())),
            (Values::I64(v), _) => v.iter().for_each(|x| data.extend_from_slice(&x.to_le_bytes())),
            (Values::U8(v), _) => data.extend_from_slice(v),
        }
        let mut out = Vec::with_capacity(8 + data.len());
        out.extend_from_slice(&(data.len() as u64).to_le_bytes());
        out.extend_from_slice(&data);
        out
    }
}

fn finite(x: f64) -> f64 {
    if x.is_finite() { x } else { 0.0 }
}

fn int32(x: f64) -> i32 {
    let x = finite(x).round();
    if x >= f64::from(i32::MAX) {
        i32::MAX
    } else if x <= f64::from(i32::MIN) {
        i32::MIN
    } else {
        x as i32
    }
}

/// The shortest text that reads back to the same `f64` (Rust's `Debug` form: `1.0`, `0.1`,
/// `1e-7`), which C and C++ number parsers read; non-finite values as 0.
fn number(x: f64) -> String {
    format!("{:?}", finite(x))
}

fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            c if c.is_control() => out.push('_'),
            c => out.push(c),
        }
    }
    out
}

/// Standard base64 (RFC 4648, with padding).
fn base64(bytes: &[u8], out: &mut String) {
    const ALPHABET: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let sym = |i: u32| char::from(ALPHABET.get((i & 63) as usize).copied().unwrap_or(b'A'));
    for chunk in bytes.chunks(3) {
        let b = |k: usize| u32::from(chunk.get(k).copied().unwrap_or(0));
        let n = (b(0) << 16) | (b(1) << 8) | b(2);
        out.push(sym(n >> 18));
        out.push(sym(n >> 12));
        out.push(if chunk.len() > 1 { sym(n >> 6) } else { '=' });
        out.push(if chunk.len() > 2 { sym(n) } else { '=' });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A model of one element of `kind` with nodes at its reference positions (labels from 1).
    fn one(kind: ElementKind) -> Model {
        let p: &[[f64; 3]] = match kind {
            ElementKind::Rod2 | ElementKind::Beam2 => &[[0., 0., 0.], [1., 0., 0.]],
            ElementKind::Rod3 => &[[0., 0., 0.], [1., 0., 0.], [0.5, 0., 0.]],
            ElementKind::Tri3 => &[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]],
            ElementKind::Tri6 => &[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0.5, 0., 0.], [0.5, 0.5, 0.], [0., 0.5, 0.]],
            ElementKind::Quad4 => &[[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]],
            ElementKind::Quad8 => {
                &[[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0.5, 0., 0.], [1., 0.5, 0.], [0.5, 1., 0.], [0., 0.5, 0.]]
            }
            ElementKind::Tet4 => &[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            ElementKind::Tet10 => &[
                [0., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
                [0., 0., 1.],
                [0.5, 0., 0.],
                [0.5, 0.5, 0.],
                [0., 0.5, 0.],
                [0., 0., 0.5],
                [0.5, 0., 0.5],
                [0., 0.5, 0.5],
            ],
            ElementKind::Wedge6 => &[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [0., 1., 1.]],
            ElementKind::Pyramid5 => &[[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0.5, 0.5, 1.]],
            ElementKind::Hex8 => &[[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]],
            ElementKind::Hex20 => &[
                [0., 0., 0.],
                [1., 0., 0.],
                [1., 1., 0.],
                [0., 1., 0.],
                [0., 0., 1.],
                [1., 0., 1.],
                [1., 1., 1.],
                [0., 1., 1.],
                [0.5, 0., 0.],
                [1., 0.5, 0.],
                [0.5, 1., 0.],
                [0., 0.5, 0.],
                [0.5, 0., 1.],
                [1., 0.5, 1.],
                [0.5, 1., 1.],
                [0., 0.5, 1.],
                [0., 0., 0.5],
                [1., 0., 0.5],
                [1., 1., 0.5],
                [0., 1., 0.5],
            ],
        };
        let mut m = Model::default();
        // Labels in reverse, so point indices differ from label - 1.
        let n = p.len() as u32;
        for (i, xyz) in p.iter().enumerate() {
            m.set_node(10 * (n - i as u32), *xyz).unwrap();
        }
        m.set_element(7, kind, (0..n).map(|i| 10 * (n - i)).collect()).unwrap();
        m
    }

    /// A tiny XML reader for the tests: the text of the first `DataArray` called `name`, and
    /// the attributes of its tag.
    fn array<'a>(xml: &'a str, name: &str) -> (&'a str, &'a str) {
        let at = xml.find(&format!("Name=\"{name}\"")).unwrap_or_else(|| panic!("no array {name}"));
        let open = xml[..at].rfind("<DataArray").unwrap();
        let tag_end = at + xml[at..].find('>').unwrap();
        let close = tag_end + xml[tag_end..].find("</DataArray>").unwrap();
        (&xml[open..=tag_end], xml[tag_end + 1..close].trim())
    }

    fn numbers(text: &str) -> Vec<f64> {
        text.split_whitespace().map(|t| t.parse().unwrap()).collect()
    }

    /// Every tag closes in order (a well-formedness check of the element structure).
    fn balanced(xml: &str) {
        let mut stack: Vec<String> = Vec::new();
        let mut rest = xml;
        while let Some(i) = rest.find('<') {
            let end = i + rest[i..].find('>').unwrap();
            let tag = &rest[i + 1..end];
            rest = &rest[end + 1..];
            if tag.starts_with('?') || tag.ends_with('/') {
                continue;
            }
            let name = tag.split_whitespace().next().unwrap();
            if let Some(n) = name.strip_prefix('/') {
                assert_eq!(stack.pop().as_deref(), Some(n), "{xml}");
            } else {
                stack.push(name.to_string());
            }
        }
        assert!(stack.is_empty(), "{stack:?}");
    }

    #[test]
    fn every_kind_has_its_vtk_type_and_point_order() {
        for kind in ElementKind::ALL {
            let m = one(kind);
            let xml = write_vtu(&m, &[], &[], Encoding::Ascii).unwrap();
            balanced(&xml);
            assert!(xml.contains(&format!("NumberOfPoints=\"{}\" NumberOfCells=\"1\"", kind.node_count())));
            assert_eq!(numbers(array(&xml, "types").1), [f64::from(cell_type(kind))], "{kind:?}");
            assert_eq!(numbers(array(&xml, "offsets").1), [kind.node_count() as f64]);
            assert_eq!(numbers(array(&xml, "element_label").1), [7.0]);
            // The connectivity reaches the element's nodes through the point labels.
            let labels = numbers(array(&xml, "node_label").1);
            let written: Vec<u32> = numbers(array(&xml, "connectivity").1).iter().map(|&i| labels[i as usize] as u32).collect();
            let e = m.element(7).unwrap();
            assert_eq!(written, e.nodes, "{kind:?}");
        }
        let ids: Vec<u8> = ElementKind::ALL.iter().map(|&k| cell_type(k)).collect();
        // rod2, rod3, beam2, tri3, tri6, quad4, quad8, tet4, tet10, wedge6, pyramid5, hex8, hex20
        assert_eq!(ids, [3, 21, 3, 5, 22, 9, 23, 10, 24, 13, 14, 12, 25]);
    }

    /// The points of the one cell as written, in VTK's order.
    fn cell_points(xml: &str) -> Vec<[f64; 3]> {
        let p = numbers(array(xml, "Points").1);
        numbers(array(xml, "connectivity").1).iter().map(|&i| [p[3 * i as usize], p[3 * i as usize + 1], p[3 * i as usize + 2]]).collect()
    }

    fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    }
    fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
    }
    fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }
    fn centre(p: &[[f64; 3]]) -> [f64; 3] {
        let n = p.len() as f64;
        [0, 1, 2].map(|k| p.iter().map(|q| q[k]).sum::<f64>() / n)
    }

    /// VTK's convention for each 3D cell, from its parametric coordinates: the base face's
    /// right-hand normal points to the rest of the cell (tetrahedron 0-1-2 to 3, hexahedron
    /// 0-1-2-3 to 4-7, wedge 0-1-2 to 3-5, pyramid 0-1-2-3 to 4). The signed measure is positive when the written cell
    /// follows it; VTK 9.7's vtkCellSizeFilter and vtkMeshQuality agree (a unit wedge written
    /// this way measures +0.5, reordered 0-2-1-3-5-4 it measures -0.5).
    #[test]
    fn solid_cells_have_positive_volume_in_vtk_convention() {
        for (kind, base, rest) in [
            (ElementKind::Tet4, 3, 3..4),
            (ElementKind::Tet10, 3, 3..4),
            (ElementKind::Hex8, 4, 4..8),
            (ElementKind::Hex20, 4, 4..8),
            (ElementKind::Wedge6, 3, 3..6),
            (ElementKind::Pyramid5, 4, 4..5),
        ] {
            let p = cell_points(&write_vtu(&one(kind), &[], &[], Encoding::Ascii).unwrap());
            let b = &p[..base];
            let n = cross(sub(b[1], b[0]), sub(b[2], b[0]));
            assert!(dot(n, sub(centre(&p[rest]), centre(b))) > 0.0, "{kind:?}");
        }
    }

    #[test]
    fn quadratic_mid_nodes_sit_on_the_vtk_edges() {
        // VTK's quadratic cells list the mid-edge points after the corners, edge by edge.
        let edges: [(ElementKind, &[(usize, usize)]); 4] = [
            (ElementKind::Tri6, &[(0, 1), (1, 2), (2, 0)]),
            (ElementKind::Quad8, &[(0, 1), (1, 2), (2, 3), (3, 0)]),
            (ElementKind::Tet10, &[(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]),
            (ElementKind::Hex20, &[(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]),
        ];
        for (kind, list) in edges {
            let p = cell_points(&write_vtu(&one(kind), &[], &[], Encoding::Ascii).unwrap());
            let corners = kind.corner_count();
            for (k, (a, b)) in list.iter().enumerate() {
                assert_eq!(p[corners + k], centre(&[p[*a], p[*b]]), "{kind:?} edge {k}");
            }
        }
        let p = cell_points(&write_vtu(&one(ElementKind::Rod3), &[], &[], Encoding::Ascii).unwrap());
        assert_eq!(p[2], centre(&[p[0], p[1]]));
    }

    #[test]
    fn fields_are_written_by_location_and_checked() {
        let m = one(ElementKind::Tet4);
        let f = |name: &str, location, components, values: Vec<f64>| Field { name: name.into(), location, components, values, integer: false };
        let fields = [
            f("displacement", Location::Point, 3, (0..12).map(f64::from).collect()),
            f("von_mises", Location::Point, 1, vec![0.1, f64::NAN, 3e-12, 1e300]),
            f("damage", Location::Cell, 1, vec![0.5]),
            Field { integer: true, ..f("has_result", Location::Point, 1, vec![1.0, 0.0, 1.0, 1.0]) },
        ];
        let xml = write_vtu(&m, &fields, &[("frequency", 12.5)], Encoding::Ascii).unwrap();
        balanced(&xml);
        let (tag, text) = array(&xml, "displacement");
        assert!(tag.contains("type=\"Float64\"") && tag.contains("NumberOfComponents=\"3\""), "{tag}");
        assert_eq!(numbers(text).len(), 12);
        assert_eq!(array(&xml, "von_mises").1, "0.1 0.0 3e-12 1e300", "shortest round trip; NaN as 0");
        assert!(array(&xml, "has_result").0.contains("Int32"));
        let cell = &xml[xml.find("<CellData>").unwrap()..xml.find("</CellData>").unwrap()];
        assert!(cell.contains("Name=\"damage\"") && !cell.contains("displacement"));
        let fd = &xml[xml.find("<FieldData>").unwrap()..xml.find("</FieldData>").unwrap()];
        assert!(fd.contains("Name=\"frequency\" NumberOfTuples=\"1\"") && fd.contains("12.5"), "{fd}");
        assert_eq!(xml, write_vtu(&m, &fields, &[("frequency", 12.5)], Encoding::Ascii).unwrap(), "deterministic");

        assert!(write_vtu(&m, &[f("short", Location::Point, 3, vec![1.0])], &[], Encoding::Ascii).is_err());
        assert!(write_vtu(&m, &[f("none", Location::Cell, 0, vec![])], &[], Encoding::Ascii).is_err());
        assert!(write_vtu(&m, &[f("a\"<b", Location::Cell, 1, vec![1.0])], &[], Encoding::Ascii).unwrap().contains("Name=\"a&quot;&lt;b\""));
    }

    #[test]
    fn base64_arrays_hold_the_byte_count_and_little_endian_values() {
        let mut out = String::new();
        base64(b"Man", &mut out);
        base64(b"Ma", &mut out);
        base64(b"M", &mut out);
        assert_eq!(out, "TWFuTWE=TQ==", "RFC 4648 examples");

        let m = one(ElementKind::Tet4);
        let xml = write_vtu(&m, &[], &[], Encoding::Base64).unwrap();
        balanced(&xml);
        let (tag, text) = array(&xml, "types");
        assert!(tag.contains("format=\"binary\""));
        // 8-byte count (1) and one byte (10): 9 bytes, 12 base64 characters.
        let mut expected = String::new();
        base64(&[1, 0, 0, 0, 0, 0, 0, 0, 10], &mut expected);
        assert_eq!(text, expected);
        assert_eq!(array(&xml, "Points").1.len(), (8 + 4 * 3 * 8usize).div_ceil(3) * 4);
    }

    #[test]
    fn a_collection_lists_each_file_and_time() {
        let text = pvd(&[(0.0, "out_0000.vtu"), (0.25, "out_0001.vtu")]);
        balanced(&text);
        assert!(text.contains("<DataSet timestep=\"0.0\" part=\"0\" file=\"out_0000.vtu\"/>"), "{text}");
        assert!(text.contains("timestep=\"0.25\""));
    }
}
