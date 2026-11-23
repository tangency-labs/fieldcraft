//! Netgen's ASCII volume mesh format. Surface records carry a boundary number and two
//! domain numbers (zero means void). Coordinates are written in meters.
//! See the [VOL record layout](https://aboyle.ca/pubs/boyle2019d-workshop.pdf) and
//! [meshio's VOL connectivity table](https://github.com/nschloe/meshio/blob/main/src/meshio/netgen/_netgen.py).

use std::{collections::BTreeMap, io::Write, path::Path};

use fieldcraft_model::{Element, ElementFace, ElementKind, Model};

use crate::{Error, Result, invalid, model_error};

/// Bound the temporary face adjacency and node numbering tables before allocating them.
const MAX_CELLS: usize = 2_000_000;
const MAX_POINTS: usize = 2_000_000;

#[derive(Debug)]
pub struct Summary {
    pub points: usize,
    pub cells: usize,
    pub surfaces: usize,
}

/// Write a conforming Tet4/Hex8 mesh. Explicit material sets partition the cells. Explicit
/// boundary sets label exterior or shared faces. Section materials and CAD tags supply defaults.
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
struct Names(BTreeMap<String, (usize, String)>);

impl Names {
    fn number(&mut self, name: &str) -> Result<usize> {
        let token: String = name.chars().map(|c| if c.is_whitespace() || c.is_control() || c == '#' { '_' } else { c }).collect();
        let next = self.0.len() + 1;
        let (number, original) = self.0.entry(token).or_insert_with(|| (next, name.to_owned()));
        if original != name {
            return invalid(format!("VOL names `{original}` and `{name}` have the same file spelling"));
        }
        Ok(*number)
    }

    fn write(&self, out: &mut dyn Write, section: &str) -> std::io::Result<()> {
        writeln!(out, "\n{section}\n{}", self.0.len())?;
        for (token, (number, _)) in &self.0 {
            writeln!(out, "{number} {token}")?;
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
            return invalid(format!("VOL material set `{set}` has no volume elements"));
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

fn boundary_assignments<'a>(model: &Model, sets: &'a [String]) -> Result<BTreeMap<ElementFace, &'a str>> {
    let mut assignments = BTreeMap::new();
    for name in sets {
        let resolved = model.resolve_set(name).map_err(|e| model_error("VOL boundary", e))?;
        if resolved.faces.is_empty() {
            return invalid(format!("VOL boundary set `{name}` has no element faces"));
        }
        for face in resolved.faces {
            if let Some(other) = assignments.insert(face, name.as_str())
                && other != name
            {
                return invalid(format!("element {} face {} belongs to VOL boundary sets `{other}` and `{name}`", face.element, face.face));
            }
        }
    }
    Ok(assignments)
}

fn cad_face_name(model: &Model, face: ElementFace) -> Option<String> {
    let tags = model.classification();
    let tag = tags.elements().get(&face.element)?.faces.get(&face.face)?;
    let entity = tags.entity(*tag)?;
    Some(format!("{}.{}", entity.body, entity.entity.name()?))
}

/// Volume cells start with an outward-oriented base face.
fn cell_order(kind: ElementKind) -> Option<&'static [usize]> {
    match kind {
        ElementKind::Tet4 => Some(&[0, 2, 1, 3]),
        ElementKind::Hex8 => Some(&[0, 3, 2, 1, 4, 7, 6, 5]),
        _ => None,
    }
}

fn face_nodes(model: &Model, face: ElementFace) -> Result<Vec<u32>> {
    model.face_nodes(face).ok_or_else(|| Error::Invalid(format!("element {} has no face {}", face.element, face.face)))
}

struct Surface {
    face: ElementFace,
    descriptor: usize,
    boundary: usize,
    inside: usize,
    outside: usize,
}

fn serialize(out: &mut dyn Write, model: &Model, material_sets: &[String], boundary_sets: &[String]) -> Result<Summary> {
    let cells = model.elements().len();
    if cells == 0 {
        return invalid("the mesh has no volume elements");
    }
    if cells > MAX_CELLS || model.nodes().len() > MAX_POINTS {
        return Err(Error::Limit(format!("VOL export supports up to {MAX_CELLS} elements and {MAX_POINTS} nodes")));
    }
    if !model.ties().is_empty() {
        return invalid("VOL export requires a conforming mesh; the model has mesh ties");
    }
    for (&label, element) in model.elements() {
        if cell_order(element.kind).is_none() {
            return invalid(format!("VOL export supports Tet4 and Hex8, but element {label} is {}", element.kind.name()));
        }
    }

    let assigned_materials = material_assignments(model, material_sets)?;
    let assigned_boundaries = boundary_assignments(model, boundary_sets)?;
    let tags = model.classification();
    let mut materials = Names::default();
    let mut domains = BTreeMap::new();
    let mut points = BTreeMap::new();
    let mut faces = Vec::with_capacity(cells * 6);
    for (&label, element) in model.elements() {
        let material = match assigned_materials.get(&label) {
            Some(name) => name,
            None if material_sets.is_empty() => {
                tags.elements().get(&label).and_then(|tag| tags.entity(tag.region?)).map_or("default", |entity| entity.body.as_str())
            }
            None => return invalid(format!("element {label} is missing from the VOL material sets")),
        };
        domains.insert(label, materials.number(material)?);
        points.extend(element.nodes.iter().map(|&node| (node, 0)));
        for i in 0..element.kind.faces().len() {
            let face = ElementFace { element: label, face: i as u8 + 1 };
            // Node labels are positive. A zero distinguishes a triangle from a quadrilateral.
            let mut key = [0; 4];
            for (slot, node) in key.iter_mut().zip(face_nodes(model, face)?) {
                *slot = node;
            }
            key.sort_unstable();
            faces.push((key, face));
        }
    }
    for (i, index) in points.values_mut().enumerate() {
        *index = i + 1;
    }
    let domain =
        |face: ElementFace| domains.get(&face.element).copied().ok_or_else(|| Error::Invalid(format!("no VOL domain for element {}", face.element)));

    let mut boundaries = Names::default();
    let mut descriptors = BTreeMap::new();
    let mut surfaces = Vec::new();
    faces.sort_unstable();
    for group in faces.chunk_by(|a, b| a.0 == b.0) {
        let (face, other) = match group {
            [(_, face)] => (*face, None),
            [(_, face), (_, other)] => (*face, Some(*other)),
            _ => return invalid("VOL export found a face shared by more than two volume elements"),
        };
        let mut inside = domain(face)?;
        let mut outside = other.map(domain).transpose()?.unwrap_or(0);
        let explicit = assigned_boundaries.get(&face).copied();
        let opposite = other.and_then(|side| assigned_boundaries.get(&side).copied());
        if let (Some(a), Some(b)) = (explicit, opposite)
            && a != b
        {
            return invalid(format!("shared face belongs to VOL boundary sets `{a}` and `{b}`"));
        }
        // Keep one orientation for an interface even when the element labels interleave.
        let (face, other) = match other {
            Some(other) if inside > outside => {
                std::mem::swap(&mut inside, &mut outside);
                (other, Some(face))
            }
            _ => (face, other),
        };
        let name = match explicit.or(opposite) {
            Some(name) => name.to_owned(),
            None if inside != outside => cad_face_name(model, face)
                .or_else(|| other.and_then(|side| cad_face_name(model, side)))
                .unwrap_or_else(|| if outside == 0 { "default".into() } else { "interface".into() }),
            None => continue,
        };
        let boundary = boundaries.number(&name)?;
        let next = descriptors.len() + 1;
        let descriptor = *descriptors.entry((boundary, inside, outside)).or_insert(next);
        surfaces.push(Surface { face, descriptor, boundary, inside, outside });
    }

    let write = |out: &mut dyn Write| -> std::io::Result<()> {
        writeln!(out, "mesh3d\ndimension\n3\ngeomtype\n0\n\nsurfaceelements\n{}", surfaces.len())?;
        for s in &surfaces {
            let nodes = face_nodes(model, s.face).map_err(std::io::Error::other)?;
            write!(out, "{} {} {} {} {}", s.descriptor, s.boundary, s.inside, s.outside, nodes.len())?;
            write_nodes(out, &points, nodes)?;
        }
        writeln!(out, "\nvolumeelements\n{cells}")?;
        for (&label, Element { kind, nodes }) in model.elements() {
            let order = cell_order(*kind).ok_or_else(|| std::io::Error::other("unsupported VOL element"))?;
            let domain = domains.get(&label).ok_or_else(|| std::io::Error::other("missing VOL domain"))?;
            write!(out, "{domain} {}", order.len())?;
            let nodes = order
                .iter()
                .map(|&i| nodes.get(i).copied())
                .collect::<Option<Vec<_>>>()
                .ok_or_else(|| std::io::Error::other(format!("invalid connectivity for element {label}")))?;
            write_nodes(out, &points, nodes)?;
        }
        writeln!(out, "\npoints\n{}", points.len())?;
        for label in points.keys() {
            let node = model.node(*label).ok_or_else(|| std::io::Error::other(format!("missing VOL node {label}")))?;
            let [x, y, z] = node.xyz;
            writeln!(out, "{x:e} {y:e} {z:e}")?;
        }
        materials.write(out, "materials")?;
        boundaries.write(out, "bcnames")?;
        writeln!(out, "\nendmesh")
    };
    write(out).map_err(|e| Error::Io(format!("writing VOL mesh: {e}")))?;
    Ok(Summary { points: points.len(), cells, surfaces: surfaces.len() })
}

fn write_nodes(out: &mut dyn Write, points: &BTreeMap<u32, usize>, nodes: impl IntoIterator<Item = u32>) -> std::io::Result<()> {
    for label in nodes {
        let index = points.get(&label).ok_or_else(|| std::io::Error::other(format!("missing VOL node {label}")))?;
        write!(out, " {index}")?;
    }
    writeln!(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_model::{GeometryRef, Set};

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
        assert_eq!((summary.points, summary.cells, summary.surfaces), (5, 2, 7));
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
        assert_eq!((summary.points, summary.cells, summary.surfaces), (8, 1, 6));
        let text = String::from_utf8(bytes).unwrap();
        assert!(text.contains("surfaceelements\n6\n1 1 1 0 4 1 4 3 2\n"));
        assert!(text.contains("volumeelements\n1\n1 8 1 4 3 2 5 8 7 6\n"));
        assert!(text.contains("materials\n1\n1 Air_Box\n"));
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
        assert!(error.to_string().contains("shared face belongs to VOL boundary sets"));
    }
}
