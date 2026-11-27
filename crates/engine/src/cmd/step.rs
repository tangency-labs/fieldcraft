//! Export the placed CAD bodies, resolving geometry-bound set names before writing the file.

use std::collections::BTreeMap;
use std::path::Path;

use fieldcraft_geometry::{Design, Shape, StepBody, export_named_step, partition::base_name};
use fieldcraft_io::roots::Access;
use fieldcraft_model::{Entity, Model};
use serde_json::{Value, json};

use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const EXPORT: CommandSpec = CommandSpec {
    id: "geometry.export.step",
    label: "Export Named STEP",
    params: r#"{"path": "named.step"} — writes the current placed CAD bodies to a STEP file before or after meshing. Body and face references in model.set.create sets supply the solid and boundary names; every selected face receives the set name. Other bodies and faces keep their CAD names. Mesh, edge and vertex members do not supply names. A missing body or face or two different set names for the same entity stop the export. An existing CAD file is protected; choose a new STEP file name."#,
    form: Some(super::forms::EXPORT_STEP),
    kind: CommandKind::Session,
    requires: &[Requirement::Design],
    run: export,
};

#[derive(Default)]
struct SetNames<'a> {
    body: Option<&'a str>,
    // Each entry pairs a persistent face selector with the set that names it.
    faces: Vec<(&'a str, &'a str)>,
}

fn conflict(body: &str, face: Option<&str>, first: &str, second: &str) -> EngineError {
    let entity = face.map_or_else(|| format!("body `{body}`"), |f| format!("face `{f}` of body `{body}`"));
    EngineError::Params(format!("sets `{first}` and `{second}` both name {entity}; STEP can store one name per entity"))
}

fn face_names<'a>(shape: &Shape, selectors: &[(&str, &'a str)]) -> Result<BTreeMap<u32, &'a str>> {
    let mut names = BTreeMap::new();
    if selectors.is_empty() {
        return Ok(names);
    }
    // Index once for all sets. A set selects both its literal name and all numbered pieces
    // of that name (GeometryRef::selects), including pieces also named by a separate set.
    let mut indexed: BTreeMap<&str, Vec<u32>> = BTreeMap::new();
    for (i, face) in shape.faces().iter().enumerate() {
        let i = u32::try_from(i).map_err(|_| EngineError::Params(format!("too many faces in body `{}`", shape.name)))?;
        indexed.entry(&face.name).or_default().push(i);
        let base = base_name(&face.name);
        if base != face.name {
            indexed.entry(base).or_default().push(i);
        }
    }
    for &(selector, set) in selectors {
        let indices = indexed
            .get(selector)
            .ok_or_else(|| EngineError::Params(format!("set `{set}` refers to missing face `{selector}` of body `{}`", shape.name)))?;
        for &i in indices {
            if let Some(first) = names.insert(i, set)
                && first != set
            {
                return Err(conflict(&shape.name, Some(selector), first, set));
            }
        }
    }
    Ok(names)
}

fn named_bodies<'a>(design: &'a Design, model: &'a Model) -> Result<Vec<StepBody<'a>>> {
    let mut by_body: BTreeMap<&str, SetNames<'_>> = BTreeMap::new();
    for (set, members) in model.sets() {
        for reference in &members.geometry {
            match &reference.entity {
                Entity::Body => {
                    let names = by_body.entry(&reference.body).or_default();
                    if let Some(first) = names.body.replace(set) {
                        return Err(conflict(&reference.body, None, first, set));
                    }
                }
                Entity::Face(face) => by_body.entry(&reference.body).or_default().faces.push((face, set)),
                Entity::Edge(_) | Entity::Vertex(_) => {}
            }
        }
    }
    let mut bodies = Vec::with_capacity(design.bodies.len());
    for shape in &design.bodies {
        let names = by_body.remove(shape.name.as_str()).unwrap_or_default();
        bodies.push(StepBody { shape, name: names.body.unwrap_or(&shape.name), face_names: face_names(shape, &names.faces)? });
    }
    if let Some((body, names)) = by_body.first_key_value() {
        let set = names.body.or_else(|| names.faces.first().map(|&(_, set)| set)).unwrap_or_default();
        return Err(EngineError::Params(format!("set `{set}` refers to missing CAD body `{body}`")));
    }
    Ok(bodies)
}

fn export(s: &mut Session, p: &Value) -> Result<Value> {
    fieldcraft_model::progress::check()?;
    let path = super::file_param(s, p, "path", Access::Write)?;
    let target = Path::new(&path);
    if !target.extension().and_then(|e| e.to_str()).is_some_and(|e| e.eq_ignore_ascii_case("step") || e.eq_ignore_ascii_case("stp")) {
        return Err(EngineError::invalid("path", "choose a STEP file name ending in .step or .stp"));
    }
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a CAD design first".into()))?;
    let bodies = named_bodies(design, &s.model)?;
    let step = export_named_step(&bodies)?;
    fieldcraft_model::progress::check()?;
    fieldcraft_io::fs::write_output(target, step.as_bytes())?;
    Ok(json!({"path": path, "bodies": bodies.len(), "bytes": step.len()}))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn airbox() -> (crate::tests::TempDir, Session, std::path::PathBuf, Value) {
        let (dir, mut s) = super::super::geometry::tests::design(json!([
            {"command": "solid.box", "params": {"length": 5, "width": 4, "height": 3, "body_name": "Box"}}
        ]));
        let path = Path::new(s.sources.design.as_ref().unwrap().path()).with_file_name("named.step");
        let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
        (dir, s, path, bodies[0]["faces"].clone())
    }

    #[test]
    fn exports_body_and_shared_boundary_set_names_without_a_mesh() {
        let (_dir, mut s, path, faces) = airbox();
        let refs: Vec<_> = faces.as_array().unwrap().iter().take(2).map(|f| json!({"body":"Box","face":f["name"]})).collect();
        s.execute("model.set.create", &json!({"name":"outer","geometry":refs})).unwrap();
        s.execute("model.set.create", &json!({"name":"air","geometry":[{"body":"Box"}]})).unwrap();
        let journal = s.journal.len();
        let reply = s.execute("geometry.export.step", &json!({"path":path})).unwrap();
        let step = std::fs::read_to_string(&path).unwrap();
        assert_eq!(reply["bytes"], step.len());
        assert_eq!(reply["bodies"], 1);
        assert_eq!(step.matches("ADVANCED_FACE('outer',").count(), 2);
        assert!(step.contains("MANIFOLD_SOLID_BREP('air',"));
        assert!(step.contains(&format!("ADVANCED_FACE('{}',", faces[2]["name"].as_str().unwrap())));
        assert_eq!(s.journal.len(), journal);
        assert!(s.execute("geometry.export.step", &json!({"path":path})).is_err());
        assert_eq!(std::fs::read_to_string(&path).unwrap(), step);
    }

    #[test]
    fn missing_or_conflicting_geometry_set_names_refuse_the_export() {
        let (_dir, mut s, path, faces) = airbox();
        let face = json!({"body":"Box","face":faces[0]["name"]});
        for (first, second, expected) in [
            (face.clone(), face, "both name face"),
            (json!({"body":"Box"}), json!({"body":"Box"}), "both name body"),
            (json!({"body":"Box"}), json!({"body":"gone"}), "missing CAD body"),
            (json!({"body":"Box"}), json!({"body":"Box","face":"gone"}), "missing face"),
        ] {
            s.execute("model.set.create", &json!({"name":"first","geometry":[first]})).unwrap();
            s.execute("model.set.create", &json!({"name":"second","geometry":[second]})).unwrap();
            let error = s.execute("geometry.export.step", &json!({"path":path})).unwrap_err().to_string();
            assert!(error.contains(expected), "{error}");
            assert!(!path.exists());
            for name in ["first", "second"] {
                s.execute("model.set.delete", &json!({"name":name})).unwrap();
            }
        }
    }
}
