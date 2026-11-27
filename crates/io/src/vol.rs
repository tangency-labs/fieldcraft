//! Netgen's ASCII mesh format. See `docs/formats/vol.md` for the supported records,
//! connectivity conventions, and sources. Coordinates are written in meters.

use std::{collections::BTreeMap, io::Write, path::Path};

use fieldcraft_mesh::predicates::{Sign, orient2d_finite, tet_volume_sign_finite};
use fieldcraft_mesh::quality::{Validity, element_validity};
use fieldcraft_model::{Element, ElementEdge, ElementFace, ElementKind, GeometryRef, Model};

use crate::{Error, Result, invalid, model_error};

/// A Hex8 mesh at this cap already needs about 275 MiB for temporary facet incidences
/// (six 24-byte entries per cell), before the node index and boundary tables.
/// Keep export memory below the cost of indexing a model at the model's much higher limits.
const MAX_CELLS: usize = 2_000_000;
const MAX_POINTS: usize = 2_000_000;

/// The name of a VOL domain and the model entity that supplied it.
#[derive(Debug, serde::Serialize)]
pub struct Domain {
    pub number: usize,
    pub name: String,
    pub source: DomainSource,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum DomainSource {
    ElementSet,
    SectionMaterial,
    CadBody,
    Unassigned,
}

#[derive(Debug)]
pub struct Summary {
    pub domains: Vec<Domain>,
    pub points: usize,
    pub cells: usize,
    pub facets: usize,
    pub dimension: u8,
}

/// Write a linear XY Tri3/Quad4 or 3D Tet4/Hex8 mesh. Explicit material sets partition the cells.
/// Boundary sets label exterior or shared edges in 2D and faces in 3D.
pub fn write(path: &Path, model: &Model, materials: &[String], boundaries: &[String]) -> Result<Summary> {
    #[cfg(not(target_arch = "wasm32"))]
    {
        crate::fs::write_output_stream(path, |out| serialize(out, model, materials, boundaries))
    }
    #[cfg(target_arch = "wasm32")]
    {
        let mut bytes = Vec::new();
        let summary = serialize(&mut bytes, model, materials, boundaries)?;
        crate::fs::write_output(path, &bytes)?;
        Ok(summary)
    }
}

/// One name per number, with a check that conversion to a VOL token does not conflate names.
#[derive(Default)]
struct Names<'a>(BTreeMap<String, Name<'a>>);

struct Name<'a> {
    number: usize,
    spelling: String,
    geometry: Option<&'a GeometryRef>,
}

fn file_name(name: &str) -> String {
    name.chars().map(|c| if c.is_whitespace() || c.is_control() || c == '#' { '_' } else { c }).collect()
}

impl<'a> Names<'a> {
    fn number(&mut self, name: &str, geometry: Option<&'a GeometryRef>) -> Result<usize> {
        let token = file_name(name);
        let next = self.0.len() + 1;
        let entry = self.0.entry(token).or_insert_with(|| Name { number: next, spelling: name.to_owned(), geometry });
        if entry.spelling != name {
            return invalid(format!("VOL names `{}` and `{name}` have the same file spelling", entry.spelling));
        }
        if entry.geometry != geometry {
            return match (entry.geometry, geometry) {
                (Some(first), Some(second)) => {
                    invalid(format!("CAD boundaries `{first}` and `{second}` both have VOL name `{name}`. Assign distinct boundary sets."))
                }
                _ => invalid(format!("A CAD boundary and a named group both have VOL name `{name}`. Assign distinct boundary sets.")),
            };
        }
        Ok(entry.number)
    }

    fn write(&self, out: &mut dyn Write, section: &str) -> std::io::Result<()> {
        writeln!(out, "\n{section}\n{}", self.0.len())?;
        for (token, name) in &self.0 {
            writeln!(out, "{} {token}", name.number)?;
        }
        Ok(())
    }

    fn write_planar_boundaries(&self, out: &mut dyn Write, descriptors: &BTreeMap<(usize, usize, usize), usize>) -> std::io::Result<()> {
        writeln!(out, "\nbcnames\n{}", descriptors.len())?;
        for (token, name) in &self.0 {
            let number = name.number;
            for (_, descriptor) in descriptors.range((number, 0, 0)..=(number, usize::MAX, usize::MAX)) {
                writeln!(out, "{descriptor} {token}")?;
            }
        }
        Ok(())
    }
}

fn material_assignments<'a>(model: &'a Model, sets: &'a [String]) -> Result<BTreeMap<u32, &'a str>> {
    let source: Vec<(&str, &str)> = if sets.is_empty() {
        model.sections().values().map(|section| (section.set.as_str(), section.material.as_str())).collect()
    } else {
        sets.iter().map(|name| (name.as_str(), name.as_str())).collect()
    };
    let mut assignments = BTreeMap::new();
    for (set, name) in source {
        let resolved = model.resolve_set(set).map_err(|e| model_error("VOL material", e))?;
        if resolved.elements.is_empty() {
            return invalid(format!("VOL material set `{set}` has no elements"));
        }
        for element in resolved.elements {
            if let Some(other) = assignments.insert(element, name)
                && other != name
            {
                return invalid(format!("element {element} has conflicting VOL materials `{other}` and `{name}`"));
            }
        }
    }
    Ok(assignments)
}

/// A cell's boundary edge in 2D or face in 3D.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Facet {
    element: u32,
    local: u8,
}

impl Facet {
    fn nodes(self, model: &Model, planar: bool) -> Result<Vec<u32>> {
        let nodes = if planar {
            model.edge_nodes(ElementEdge { element: self.element, edge: self.local })
        } else {
            model.face_nodes(ElementFace { element: self.element, face: self.local })
        };
        nodes.ok_or_else(|| Error::Invalid(format!("element {} has no boundary facet {}", self.element, self.local)))
    }
}

fn boundary_assignments<'a>(model: &Model, sets: &'a [String], planar: bool) -> Result<BTreeMap<Facet, &'a str>> {
    let mut assignments = BTreeMap::new();
    for name in sets {
        if matches!(name.as_str(), "default" | "interface") {
            return invalid(format!("VOL boundary set `{name}` uses a reserved fallback name; choose a different boundary set name"));
        }
        let resolved = model.resolve_set(name).map_err(|e| model_error("VOL boundary", e))?;
        let facets: Vec<_> = if planar {
            resolved.edges.into_iter().map(|edge| Facet { element: edge.element, local: edge.edge }).collect()
        } else {
            resolved.faces.into_iter().map(|face| Facet { element: face.element, local: face.face }).collect()
        };
        if facets.is_empty() {
            let kind = if planar { "element edges" } else { "element faces" };
            return invalid(format!("VOL boundary set `{name}` has no {kind}"));
        }
        for facet in facets {
            if let Some(other) = assignments.insert(facet, name.as_str())
                && other != name
            {
                return invalid(format!("element {} facet {} belongs to VOL boundary sets `{other}` and `{name}`", facet.element, facet.local));
            }
        }
    }
    Ok(assignments)
}

fn cad_facet_name(model: &Model, facet: Facet, planar: bool) -> Option<(String, &GeometryRef)> {
    let tags = model.classification();
    let element = tags.elements().get(&facet.element)?;
    let indices = if planar { &element.edges } else { &element.faces };
    let entity = tags.entity(*indices.get(&facet.local)?)?;
    Some((format!("{}.{}", entity.body, entity.entity.name()?), entity))
}

/// Planar cells are counterclockwise. Volume cells start with an outward-oriented base face.
fn cell_order(kind: ElementKind) -> Option<&'static [usize]> {
    match kind {
        ElementKind::Tri3 => Some(&[0, 1, 2]),
        ElementKind::Quad4 => Some(&[0, 1, 2, 3]),
        ElementKind::Tet4 => Some(&[0, 2, 1, 3]),
        ElementKind::Hex8 => Some(&[0, 3, 2, 1, 4, 7, 6, 5]),
        _ => None,
    }
}

fn validate_planar_cell(model: &Model, label: u32, element: &Element) -> Result<()> {
    let point = |index| -> Result<[f64; 2]> {
        let node = element.nodes.get(index).ok_or_else(|| Error::Invalid(format!("invalid connectivity for element {label}")))?;
        let [x, y, _] = model.node(*node).ok_or_else(|| Error::Invalid(format!("missing VOL node {node}")))?.xyz;
        Ok([x, y])
    };
    let positive = |a, b, c| {
        orient2d_finite(a, b, c).map(|sign| sign == Sign::Positive).map_err(|error| Error::Invalid(format!("VOL element {label}: {error}")))
    };
    let (a, b, c) = (point(0)?, point(1)?, point(2)?);
    let mut counterclockwise = positive(a, b, c)?;
    if counterclockwise && element.kind == ElementKind::Quad4 {
        let d = point(3)?;
        counterclockwise = positive(b, c, d)? && positive(c, d, a)? && positive(d, a, b)?;
    }
    if !counterclockwise {
        return invalid(format!("VOL planar element {label} must have counterclockwise, convex connectivity in the XY plane"));
    }
    Ok(())
}

/// A trilinear hex can fold between its corners. Check its local node ordering exactly
/// and require the mesh quality certificate for the entire reference cube.
fn validate_volume_cell(model: &Model, label: u32, element: &Element) -> Result<()> {
    let mut points = [[0.; 3]; 8];
    let xyz = points.get_mut(..element.nodes.len()).ok_or_else(|| Error::Invalid(format!("invalid connectivity for VOL element {label}")))?;
    for (point, node) in xyz.iter_mut().zip(&element.nodes) {
        *point = model.node(*node).ok_or_else(|| Error::Invalid(format!("missing VOL node {node}")))?.xyz;
    }
    let positive =
        |a, b, c, d| tet_volume_sign_finite(a, b, c, d).map(|s| s == Sign::Positive).map_err(|e| Error::Invalid(format!("VOL element {label}: {e}")));
    let valid = match *xyz {
        [a, b, c, d] if element.kind == ElementKind::Tet4 => positive(a, b, c, d)?,
        [a, b, c, d, e, f, g, h] if element.kind == ElementKind::Hex8 => {
            let corners = [[a, b, d, e], [b, c, a, f], [c, d, b, g], [d, a, c, h], [e, h, f, a], [f, e, g, b], [g, f, h, c], [h, g, e, d]];
            for [o, u, v, w] in corners {
                if !positive(o, u, v, w)? {
                    return invalid(format!("VOL Hex8 element {label} has an inverted or degenerate corner"));
                }
            }
            element_validity(ElementKind::Hex8, xyz).map_err(|e| Error::Invalid(format!("VOL Hex8 element {label}: {e}")))?.validity
                == Validity::Valid
        }
        _ => return invalid(format!("invalid connectivity for VOL element {label}")),
    };
    if !valid {
        return invalid(format!("VOL {} element {label} does not have a positive Jacobian throughout the element", element.kind.name()));
    }
    Ok(())
}

struct Boundary {
    facet: Facet,
    descriptor: usize,
    boundary: usize,
    inside: usize,
    outside: usize,
}

fn serialize(out: &mut dyn Write, model: &Model, material_sets: &[String], boundary_sets: &[String]) -> Result<Summary> {
    let cells = model.elements().len();
    let Some(first) = model.elements().values().next() else {
        return invalid("the mesh has no elements");
    };
    let dimension = first.kind.dimension();
    let planar = dimension == 2;
    if cells > MAX_CELLS || model.nodes().len() > MAX_POINTS {
        return Err(Error::Limit(format!("VOL export supports up to {MAX_CELLS} elements and {MAX_POINTS} nodes")));
    }
    if !model.ties().is_empty() {
        return invalid("VOL export requires a conforming mesh; the model has mesh ties");
    }
    for (&label, element) in model.elements() {
        if cell_order(element.kind).is_none() {
            return invalid(format!("VOL export supports Tri3, Quad4, Tet4, and Hex8, but element {label} is {}", element.kind.name()));
        }
        if element.kind.dimension() != dimension {
            return invalid("VOL export requires a single mesh dimension");
        }
        if planar {
            validate_planar_cell(model, label, element)?;
        } else {
            validate_volume_cell(model, label, element)?;
        }
    }

    let assigned_materials = material_assignments(model, material_sets)?;
    let assigned_boundaries = boundary_assignments(model, boundary_sets, planar)?;
    let tags = model.classification();
    let mut materials = Names::default();
    let mut domain_names = Vec::<Domain>::new();
    let mut domains = BTreeMap::new();
    let mut points = BTreeMap::new();
    let mut facets = Vec::new();
    facets.try_reserve_exact(cells * if planar { 4 } else { 6 }).map_err(|e| Error::Limit(format!("allocating VOL facet incidences: {e}")))?;
    for (&label, element) in model.elements() {
        let (material, source) = match assigned_materials.get(&label) {
            Some(&name) if material_sets.is_empty() => (name, DomainSource::SectionMaterial),
            Some(&name) => (name, DomainSource::ElementSet),
            None if material_sets.is_empty() => tags
                .elements()
                .get(&label)
                .and_then(|tag| tags.entity(tag.region?))
                .map_or(("default", DomainSource::Unassigned), |entity| (entity.body.as_str(), DomainSource::CadBody)),
            None => return invalid(format!("element {label} is missing from the VOL material sets")),
        };
        let number = materials.number(material, None)?;
        if let Some(existing) = domain_names.get(number - 1) {
            if existing.source != source {
                return invalid(format!("VOL domain `{material}` names different kinds of regions; supply distinct element sets"));
            }
        } else {
            domain_names.push(Domain { number, name: file_name(material), source });
        }
        domains.insert(label, number);
        points.extend(element.nodes.iter().map(|&node| (node, 0)));
        let sides = if planar { element.kind.edges() } else { element.kind.faces() };
        for (i, side) in sides.iter().enumerate() {
            let facet = Facet { element: label, local: i as u8 + 1 };
            // Node labels are positive. Zero padding distinguishes different facet types.
            let mut key = [0; 4];
            for (slot, &local) in key.iter_mut().zip(*side) {
                *slot = element.nodes.get(local).copied().ok_or_else(|| Error::Invalid(format!("invalid facet on element {label}")))?;
            }
            key.sort_unstable();
            facets.push((key, facet));
        }
    }
    for (i, (&label, index)) in points.iter_mut().enumerate() {
        if planar && model.node(label).is_some_and(|node| node.xyz[2] != 0.0) {
            return invalid(format!("VOL planar meshes must lie in the XY plane at z = 0, but node {label} does not"));
        }
        *index = i + 1;
    }
    let domain =
        |facet: Facet| domains.get(&facet.element).copied().ok_or_else(|| Error::Invalid(format!("no VOL domain for element {}", facet.element)));

    let mut boundaries = Names::default();
    let mut descriptors = BTreeMap::new();
    let mut boundary_elements = Vec::new();
    facets.sort_unstable();
    for group in facets.chunk_by(|a, b| a.0 == b.0) {
        let (facet, other) = match group {
            [(_, facet)] => (*facet, None),
            [(_, facet), (_, other)] => (*facet, Some(*other)),
            _ => return invalid("VOL export found a boundary facet shared by more than two elements"),
        };
        let mut inside = domain(facet)?;
        let mut outside = other.map(domain).transpose()?.unwrap_or(0);
        let explicit = assigned_boundaries.get(&facet).copied();
        let opposite = other.and_then(|side| assigned_boundaries.get(&side).copied());
        if let (Some(a), Some(b)) = (explicit, opposite)
            && a != b
        {
            return invalid(format!("shared facet belongs to VOL boundary sets `{a}` and `{b}`"));
        }
        // Keep one orientation for an interface even when the element labels interleave.
        let (facet, other) = match other {
            Some(other) if inside > outside => {
                std::mem::swap(&mut inside, &mut outside);
                (other, Some(facet))
            }
            _ => (facet, other),
        };
        let (name, geometry) = match explicit.or(opposite) {
            Some(name) => (name.to_owned(), None),
            None if inside != outside => cad_facet_name(model, facet, planar)
                .or_else(|| other.and_then(|side| cad_facet_name(model, side, planar)))
                .map_or_else(|| (if outside == 0 { "default".into() } else { "interface".into() }, None), |(name, geometry)| (name, Some(geometry))),
            None => continue,
        };
        let boundary = boundaries.number(&name, geometry)?;
        let next = descriptors.len() + 1;
        let descriptor = *descriptors.entry((boundary, inside, outside)).or_insert(next);
        boundary_elements.push(Boundary { facet, descriptor, boundary, inside, outside });
    }

    writeln!(out, "mesh3d\ndimension\n{dimension}\ngeomtype\n0").map_err(write_error)?;
    if planar {
        writeln!(out, "\nsurfaceelements\n{cells}").map_err(write_error)?;
    } else {
        writeln!(out, "\nsurfaceelements\n{}", boundary_elements.len()).map_err(write_error)?;
        for s in &boundary_elements {
            let nodes = s.facet.nodes(model, false)?;
            write!(out, "{} {} {} {} {}", s.descriptor, s.boundary, s.inside, s.outside, nodes.len()).map_err(write_error)?;
            write_nodes(out, &points, nodes)?;
        }
        writeln!(out, "\nvolumeelements\n{cells}").map_err(write_error)?;
    }
    for (&label, Element { kind, nodes }) in model.elements() {
        let order = cell_order(*kind).ok_or_else(|| Error::Invalid(format!("unsupported VOL element {label}")))?;
        let domain = domains.get(&label).ok_or_else(|| Error::Invalid(format!("missing VOL domain for element {label}")))?;
        if planar {
            write!(out, "{domain} {domain} 0 0 {}", order.len()).map_err(write_error)?;
        } else {
            write!(out, "{domain} {}", order.len()).map_err(write_error)?;
        }
        let nodes = order.iter().map(|&i| nodes.get(i).copied().ok_or_else(|| Error::Invalid(format!("invalid connectivity for element {label}"))));
        write_nodes(out, &points, nodes.collect::<Result<Vec<_>>>()?)?;
    }
    if planar {
        writeln!(out, "\nvolumeelements\n0\n\nedgesegmentsgi2\n{}", boundary_elements.len()).map_err(write_error)?;
        for s in &boundary_elements {
            let nodes = s.facet.nodes(model, true)?;
            let [start, end] = nodes.as_slice() else {
                return invalid("invalid VOL boundary segment");
            };
            let (start, end) = (point_index(&points, *start)?, point_index(&points, *end)?);
            let Boundary { inside, outside, descriptor, .. } = s;
            // In 2D the segment index identifies both its name and its domain pair.
            // The two endpoint records describe an uncurved segment without a CAD source.
            writeln!(out, "{descriptor} 0 {start} {end} -1 -1 {inside} {outside} {descriptor} 0 {descriptor} 1").map_err(write_error)?;
        }
    }
    writeln!(out, "\npoints\n{}", points.len()).map_err(write_error)?;
    for label in points.keys() {
        let node = model.node(*label).ok_or_else(|| Error::Invalid(format!("missing VOL node {label}")))?;
        let [x, y, z] = node.xyz;
        writeln!(out, "{x:e} {y:e} {z:e}").map_err(write_error)?;
    }
    materials.write(out, "materials").map_err(write_error)?;
    if planar {
        boundaries.write_planar_boundaries(out, &descriptors).map_err(write_error)?;
    } else {
        boundaries.write(out, "bcnames").map_err(write_error)?;
    }
    writeln!(out, "\nendmesh").map_err(write_error)?;
    Ok(Summary { domains: domain_names, points: points.len(), cells, facets: boundary_elements.len(), dimension })
}

fn write_error(error: std::io::Error) -> Error {
    Error::Io(format!("writing VOL mesh: {error}"))
}

fn write_nodes(out: &mut dyn Write, points: &BTreeMap<u32, usize>, nodes: impl IntoIterator<Item = u32>) -> Result<()> {
    for label in nodes {
        let index = point_index(points, label)?;
        write!(out, " {index}").map_err(write_error)?;
    }
    writeln!(out).map_err(write_error)
}

fn point_index(points: &BTreeMap<u32, usize>, label: u32) -> Result<usize> {
    points.get(&label).copied().ok_or_else(|| Error::Invalid(format!("missing VOL node {label}")))
}

#[cfg(test)]
mod format_tests;

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_model::{Entity, Set};

    fn two_tets() -> Model {
        let mut model = Model::default();
        for (label, xyz) in [
            (10, [0.0, 0.0, 0.0]),
            (20, [1.0, 0.0, 0.0]),
            (30, [0.0, 1.0, 0.0]),
            (40, [0.0, 0.0, 1.0]),
            (50, [0.0, 0.0, -1.0]),
            (60, [2.0, 2.0, 2.0]),
        ] {
            model.set_node(label, xyz).unwrap();
        }
        model.set_element(70, ElementKind::Tet4, vec![10, 20, 30, 40]).unwrap();
        model.set_element(80, ElementKind::Tet4, vec![10, 30, 20, 50]).unwrap();
        model
    }

    fn set(model: &mut Model, name: &str, elements: Vec<u32>, faces: Vec<ElementFace>) {
        model.add_set(name, Set { elements, faces, ..Set::default() }).unwrap();
    }

    #[test]
    fn writes_one_oriented_face_for_a_shared_named_interface() {
        let mut model = two_tets();
        set(&mut model, "steel", vec![70], vec![]);
        set(&mut model, "air", vec![80], vec![]);
        let shared = [ElementFace { element: 70, face: 1 }, ElementFace { element: 80, face: 1 }];
        set(&mut model, "joint", vec![], shared.to_vec());
        let cad_face = model.geometry_tag(&GeometryRef::face("Steel Plate", "Top").unwrap()).unwrap();
        model.tag_element_face(70, 2, cad_face).unwrap();

        let mut bytes = Vec::new();
        let summary = serialize(&mut bytes, &model, &["steel".into(), "air".into()], &["joint".into()]).unwrap();
        assert_eq!((summary.points, summary.cells, summary.facets), (5, 2, 7));
        let text = String::from_utf8(bytes).unwrap();
        assert!(text.contains("surfaceelements\n7\n1 1 1 2 3 1 3 2\n"));
        assert!(text.contains("volumeelements\n2\n1 4 1 3 2 4\n2 4 1 2 3 5\n"));
        assert!(text.contains("materials\n2\n2 air\n1 steel"));
        assert!(text.contains("1 joint\n"));
        assert!(text.contains(" Steel_Plate.Top\n"));
    }

    #[test]
    fn keeps_hex_quads_and_writes_the_vol_reference_order() {
        let mut model = Model::default();
        for (i, xyz) in
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]]
                .into_iter()
                .enumerate()
        {
            model.set_node(i as u32 + 1, xyz).unwrap();
        }
        model.set_element(10, ElementKind::Hex8, (1..=8).collect()).unwrap();
        let body = model.geometry_tag(&GeometryRef::body("Air Box").unwrap()).unwrap();
        model.tag_element(10, body).unwrap();
        let mut bytes = Vec::new();
        let summary = serialize(&mut bytes, &model, &[], &[]).unwrap();
        assert_eq!((summary.points, summary.cells, summary.facets), (8, 1, 6));
        let text = String::from_utf8(bytes).unwrap();
        assert!(text.contains("surfaceelements\n6\n1 1 1 0 4 1 4 3 2\n"));
        assert!(text.contains("volumeelements\n1\n1 8 1 4 3 2 5 8 7 6\n"));
        assert!(text.contains("materials\n1\n1 Air_Box\n"));
    }

    #[test]
    fn writes_planar_materials_and_one_segment_for_a_named_shared_edge() {
        let mut model = Model::default();
        for (i, xyz) in [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]].into_iter().enumerate() {
            model.set_node(i as u32 + 1, xyz).unwrap();
        }
        model.set_element(10, ElementKind::Tri3, vec![1, 2, 3]).unwrap();
        model.set_element(20, ElementKind::Tri3, vec![2, 4, 3]).unwrap();
        set(&mut model, "air", vec![10], vec![]);
        set(&mut model, "beam", vec![20], vec![]);
        let edges = vec![ElementEdge { element: 10, edge: 2 }, ElementEdge { element: 20, edge: 3 }];
        model.add_set("electrode", Set { edges, ..Set::default() }).unwrap();

        let mut bytes = Vec::new();
        let summary = serialize(&mut bytes, &model, &["air".into(), "beam".into()], &["electrode".into()]).unwrap();
        assert_eq!((summary.dimension, summary.points, summary.cells, summary.facets), (2, 4, 2, 5));
        let text = String::from_utf8(bytes).unwrap();
        super::format_tests::check_split_square(&text);
        assert!(summary.domains.iter().all(|domain| domain.source == DomainSource::ElementSet));

        let first = model.geometry_tag(&GeometryRef::new("a.b", Entity::Edge("c".into())).unwrap()).unwrap();
        let second = model.geometry_tag(&GeometryRef::new("a", Entity::Edge("b.c".into())).unwrap()).unwrap();
        model.tag_element_edge(20, 1, first).unwrap();
        model.tag_element_edge(20, 2, second).unwrap();
        let error = serialize(&mut Vec::new(), &model, &[], &[]).unwrap_err();
        assert!(error.to_string().contains("both have VOL name `a.b.c`"));

        model.set_node(4, [1.0, 1.0, 0.1]).unwrap();
        let error = serialize(&mut Vec::new(), &model, &[], &[]).unwrap_err();
        assert!(error.to_string().contains("node 4"));
    }

    #[test]
    fn writes_planar_cells_with_positive_jacobians() {
        let points = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0]];
        for kind in [ElementKind::Tri3, ElementKind::Quad4] {
            let mut model = Model::default();
            for (i, &point) in points.iter().take(kind.node_count()).enumerate() {
                model.set_node(i as u32 + 1, point).unwrap();
            }
            model.set_element(1, kind, (1..=kind.node_count() as u32).collect()).unwrap();

            let mut bytes = Vec::new();
            serialize(&mut bytes, &model, &[], &[]).unwrap();
            let text = String::from_utf8(bytes).unwrap();
            let connectivity = text.split_once("surfaceelements\n1\n").unwrap().1.lines().next().unwrap();
            let mut vertices = connectivity.split_whitespace().skip(5).map(|id| points[id.parse::<usize>().unwrap() - 1]);
            let [ax, ay, _] = vertices.next().unwrap();
            let [bx, by, _] = vertices.next().unwrap();
            let [cx, cy, _] = vertices.last().unwrap();
            let determinant = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
            assert!(determinant > 0.0, "{kind:?} has signed Jacobian {determinant}");

            model.set_element(1, kind, (1..=kind.node_count() as u32).rev().collect()).unwrap();
            let error = serialize(&mut Vec::new(), &model, &[], &[]).unwrap_err();
            assert!(matches!(error, Error::Invalid(message) if message.contains("counterclockwise")));
        }
    }

    #[test]
    fn preserves_nonzero_coordinates_below_the_expansion_domain() {
        let mut model = Model::default();
        let tiny = 1e-40;
        model.set_node(1, [tiny, 0., 0.]).unwrap();
        model.set_node(2, [1., 0., 0.]).unwrap();
        model.set_node(3, [0., 1., 0.]).unwrap();
        model.set_element(1, ElementKind::Tri3, vec![1, 2, 3]).unwrap();
        let mut bytes = Vec::new();
        serialize(&mut bytes, &model, &[], &[]).unwrap();
        let text = String::from_utf8(bytes).unwrap();
        assert!(text.contains("points\n3\n1e-40 0e0 0e0\n"));
    }

    #[test]
    fn refuses_ambiguous_faces_and_incomplete_materials() {
        let mut model = two_tets();
        set(&mut model, "steel", vec![70], vec![]);
        let error = serialize(&mut Vec::new(), &model, &["steel".into()], &[]).unwrap_err();
        assert!(error.to_string().contains("element 80 is missing"));

        set(&mut model, "first", vec![], vec![ElementFace { element: 70, face: 1 }]);
        set(&mut model, "second", vec![], vec![ElementFace { element: 80, face: 1 }]);
        let error = serialize(&mut Vec::new(), &model, &[], &["first".into(), "second".into()]).unwrap_err();
        assert!(error.to_string().contains("shared facet belongs to VOL boundary sets"));
    }

    #[test]
    fn refuses_distinct_cad_faces_with_the_same_qualified_name() {
        let mut model = two_tets();
        let first = model.geometry_tag(&GeometryRef::face("a.b", "c").unwrap()).unwrap();
        let second = model.geometry_tag(&GeometryRef::face("a", "b.c").unwrap()).unwrap();
        model.tag_element_face(70, 2, first).unwrap();
        model.tag_element_face(80, 2, first).unwrap();

        let mut bytes = Vec::new();
        serialize(&mut bytes, &model, &[], &[]).unwrap();
        assert_eq!(String::from_utf8(bytes).unwrap().matches(" a.b.c\n").count(), 1);

        model.tag_element_face(80, 2, second).unwrap();
        let error = serialize(&mut Vec::new(), &model, &[], &[]).unwrap_err().to_string();
        assert!(error.contains("face c of body a.b") && error.contains("face b.c of body a"));

        set(&mut model, "second", vec![], vec![ElementFace { element: 80, face: 2 }]);
        serialize(&mut Vec::new(), &model, &[], &["second".into()]).unwrap();
    }

    #[test]
    fn refuses_reserved_boundary_names_before_writing() {
        for name in ["default", "interface"] {
            let mut model = two_tets();
            set(&mut model, name, vec![], vec![ElementFace { element: 70, face: 2 }]);
            let mut bytes = Vec::new();
            let error = serialize(&mut bytes, &model, &[], &[name.into()]).unwrap_err();
            assert!(error.to_string().contains(&format!("`{name}` uses a reserved")));
            assert!(bytes.is_empty());
        }
    }

    #[test]
    fn rejects_inverted_or_degenerate_volume_cells() {
        let mut model = two_tets();
        for z in [-1., 0.] {
            model.set_node(40, [0., 0., z]).unwrap();
            let error = serialize(&mut Vec::new(), &model, &[], &[]).unwrap_err();
            assert!(error.to_string().contains("element 70"));
        }
        let mut model = Model::default();
        let bottom = [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]];
        for (i, [x, y]) in bottom.into_iter().enumerate() {
            model.set_node(i as u32 + 1, [x, y, 0.]).unwrap();
            model.set_node(i as u32 + 5, [x, y, -1.]).unwrap();
        }
        model.set_element(1, ElementKind::Hex8, (1..=8).collect()).unwrap();
        assert!(serialize(&mut Vec::new(), &model, &[], &[]).unwrap_err().to_string().contains("inverted or degenerate corner"));
        // All eight corners have positive Jacobians, but the trilinear map folds inside.
        for (i, [x, y]) in bottom.into_iter().enumerate() {
            model.set_node(i as u32 + 5, [-2. * x, -3. * y, 1.]).unwrap();
        }
        assert!(serialize(&mut Vec::new(), &model, &[], &[]).unwrap_err().to_string().contains("positive Jacobian throughout"));
    }

    #[test]
    fn domain_names_distinguish_assigned_materials_bodies_and_sets() {
        let mut model = two_tets();
        set(&mut model, "top", vec![70], vec![]);
        set(&mut model, "bottom", vec![80], vec![]);
        let tag = model.geometry_tag(&GeometryRef::body("Bottom body").unwrap()).unwrap();
        model.tag_element(80, tag).unwrap();
        model.add_material("steel", fieldcraft_model::Material { density: Some(7800.), ..Default::default() }).unwrap();
        model
            .add_section(
                "upper_material",
                fieldcraft_model::Section {
                    kind: fieldcraft_model::SectionKind::Solid {},
                    set: "top".into(),
                    material: "steel".into(),
                    formulation: Default::default(),
                },
            )
            .unwrap();
        let auto = serialize(&mut Vec::new(), &model, &[], &[]).unwrap();
        assert_eq!((auto.domains[0].name.as_str(), auto.domains[0].source), ("steel", DomainSource::SectionMaterial));
        assert_eq!((auto.domains[1].name.as_str(), auto.domains[1].source), ("Bottom_body", DomainSource::CadBody));
        let ambiguous = model.geometry_tag(&GeometryRef::body("steel").unwrap()).unwrap();
        model.tag_element(80, ambiguous).unwrap();
        assert!(serialize(&mut Vec::new(), &model, &[], &[]).unwrap_err().to_string().contains("names different kinds of regions"));
        let explicit = serialize(&mut Vec::new(), &model, &["top".into(), "bottom".into()], &[]).unwrap();
        assert_eq!(
            explicit.domains.iter().map(|d| (d.number, d.name.as_str(), d.source)).collect::<Vec<_>>(),
            [(1, "top", DomainSource::ElementSet), (2, "bottom", DomainSource::ElementSet)]
        );
    }

    #[test]
    fn refuses_equal_cad_and_explicit_names() {
        let geometry = GeometryRef::face("body", "face").unwrap();
        let mut names = Names::default();
        names.number("body.face", None).unwrap();
        let error = names.number("body.face", Some(&geometry)).unwrap_err();
        assert!(matches!(error, Error::Invalid(message) if message.contains("named group")));
    }
}
