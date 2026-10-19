//! Stitching imported surfaces (`geometry.stitch`): their gaps under a tolerance closed for
//! meshing ([`fieldcraft_geometry::stitch`]).
//!
//! A stitch is part of the model (`fieldcraft_model::virtual_faces::Stitch`), naming the surface
//! by name. The surfaces keep what was read ([`TaggedSurface::original`]); [`restitch`] stitches
//! each again from it by the model's stitches, in order, whenever a surface is imported or the
//! model is opened, and after every stitch or delete, so the meshers always take the stitched
//! surface. Undo restores the surfaces with the model.

use std::sync::Arc;

use fieldcraft_geometry::tagged_mesh::TaggedSurface;
use fieldcraft_model::Stitch;
use serde_json::{Value, json};

use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const STITCH: CommandSpec = CommandSpec {
    id: "geometry.stitch",
    label: "Stitch Surface",
    params: r#"{"name": "Gaps", "surface": "Mesh1" (an imported surface, by name), "tolerance": 0.05 (mm)} — closes the gaps of an imported surface under the tolerance, for meshing: points on open edges closer than it become one point (triangles that collapse go), and a point on an open edge within it of the inside of another open edge splits that edge's triangle (a T-junction). Edges that are not open are never touched, so a closed surface is unchanged; triangles keep their face tags. Kept in the model and applied again, in order, whenever the surface is imported or the model opened (the file is not changed). Replies the open edges before and after (`open_before`, `open_after`: 0 means closed), the points `merged`, the triangles `split` and `removed`; the surface loses its mesh (`replaced`). A surface the session lacks is an error. One undo step"#,
    form: Some(super::forms::STITCH),
    kind: CommandKind::Model,
    requires: &[Requirement::Geometry],
    run: stitch,
};

pub const STITCH_DELETE: CommandSpec = CommandSpec {
    id: "geometry.stitch.delete",
    label: "Delete Stitch",
    params: r#"{"name": "Gaps"} — removes a stitch: the surface is stitched again by the ones left (as read when none); it loses its mesh (`replaced`)"#,
    form: Some(super::forms::STITCH_DELETE),
    kind: CommandKind::Model,
    requires: &[Requirement::Stitch],
    run: delete,
};

/// The name a surface goes by (and its mesh is tagged with): its own, or `Surface<index>`.
pub(crate) fn surface_name(index: usize, t: &TaggedSurface) -> String {
    if t.name.trim().is_empty() { format!("Surface{index}") } else { t.name.clone() }
}

/// Stitch the session's surfaces again from what was read, by the model's stitches (the model in
/// millimetres, as geometry commands see it). Replies each stitch's report by name, and the
/// stitches whose surface is not there or that failed (with why).
pub(crate) fn restitch(s: &mut Session) -> (Vec<(String, Value)>, Vec<(String, String)>) {
    let model = std::mem::take(&mut s.model);
    let out = restitch_surfaces(&mut s.surfaces, &model);
    s.model = model;
    out
}

/// [`restitch`] of `surfaces` by `model`'s stitches, its lengths in millimetres.
pub(crate) fn restitch_surfaces(
    surfaces: &mut [Arc<TaggedSurface>],
    model: &fieldcraft_model::Model,
) -> (Vec<(String, Value)>, Vec<(String, String)>) {
    let mut reports = Vec::new();
    let mut failed = Vec::new();
    let stitches = model.stitches().to_vec();
    for st in stitches.iter().filter(|st| !surfaces.iter().enumerate().any(|(i, t)| surface_name(i, t) == st.surface)) {
        failed.push((st.name.clone(), format!("no surface `{}`", st.surface)));
    }
    for (i, slot) in surfaces.iter_mut().enumerate() {
        let name = surface_name(i, slot);
        let mine: Vec<&Stitch> = stitches.iter().filter(|st| st.surface == name).collect();
        if mine.is_empty() && slot.original.is_none() {
            continue;
        }
        let raw = slot.original.clone().unwrap_or_else(|| Arc::new(slot.mesh.clone()));
        let mut mesh = raw.as_ref().clone();
        for st in mine {
            match fieldcraft_geometry::stitch::stitch(&mesh, st.tolerance) {
                Ok((m, r)) => {
                    mesh = m;
                    let report =
                        json!({"open_before": r.open_before, "open_after": r.open_after, "merged": r.merged, "split": r.split, "removed": r.removed});
                    reports.push((st.name.clone(), report));
                }
                Err(e) => failed.push((st.name.clone(), e.to_string())),
            }
        }
        let original = (!model.stitches().iter().all(|st| st.surface != name)).then_some(raw);
        *slot = Arc::new(TaggedSurface { name: slot.name.clone(), mesh, faces: slot.faces, original });
    }
    (reports, failed)
}

fn stitch(s: &mut Session, p: &Value) -> Result<Value> {
    let name = super::name(p)?.to_string();
    let surface: String = super::req(p, "surface")?;
    let tolerance: f64 = super::req(p, "tolerance")?;
    s.model.add_stitch(Stitch { name: name.clone(), surface: surface.clone(), tolerance })?;
    let (reports, failed) = restitch(s);
    if let Some((_, why)) = failed.iter().find(|(n, _)| *n == name) {
        return Err(EngineError::Params(format!("stitch `{name}`: {why}")));
    }
    let mut reply = reports.into_iter().find(|(n, _)| *n == name).map(|(_, r)| r).unwrap_or_else(|| json!({}));
    if let Some(o) = reply.as_object_mut() {
        o.insert("name".into(), json!(name));
        o.insert("surface".into(), json!(surface));
    }
    let cleared = super::replace::clear(&mut s.model, [surface.as_str()]);
    Ok(cleared.report_if_any(s, reply))
}

fn delete(s: &mut Session, p: &Value) -> Result<Value> {
    let name = super::name(p)?;
    let removed = s.model.remove_stitch(name)?;
    restitch(s);
    let cleared = super::replace::clear(&mut s.model, [removed.surface.as_str()]);
    Ok(cleared.report_if_any(s, json!({"deleted": name, "surface": removed.surface})))
}

#[cfg(test)]
mod tests {
    use serde_json::json;

    use crate::Session;

    /// A unit cube in a 3MF whose top face is raised 0.01 mm: open along the top's rim.
    fn gapped_cube(dir: &std::path::Path) -> String {
        let m = fieldcraft_geometry::stitch::test_cube(0.01, false);
        let vertices: String = m.positions.iter().map(|p| format!(r#"<vertex x="{}" y="{}" z="{}"/>"#, p[0], p[1], p[2])).collect();
        let triangles: String = m
            .triangles
            .iter()
            .zip(&m.tri_face)
            .map(|(t, f)| format!(r#"<triangle v1="{}" v2="{}" v3="{}" ft:face="{f}"/>"#, t[0], t[1], t[2]))
            .collect();
        let model = format!(
            r#"<model unit="millimeter" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02" xmlns:ft="{}"><resources><object id="1" name="Cube" ft:faces="6"><mesh><vertices>{vertices}</vertices><triangles>{triangles}</triangles></mesh></object></resources><build><item objectid="1"/></build></model>"#,
            fieldcraft_geometry::tagged_mesh::FACE_TAGS_NS
        );
        let path = dir.join("cube.3mf");
        std::fs::write(&path, super::super::geometry::tests::package(&model)).unwrap();
        path.to_str().unwrap().to_string()
    }

    /// The gapped cube does not mesh; stitched under a tolerance over the gap it does; the stitch
    /// is saved and applied again when the model is opened; undo and delete bring the gap back.
    #[test]
    fn a_stitched_surface_meshes_and_the_stitch_is_replayed() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-stitch-{}-{:?}", std::process::id(), std::thread::current().id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut s = Session::default();
        s.execute("geometry.import.mesh", &json!({"path": gapped_cube(&dir)})).unwrap();
        let tet = json!({"surface": 0, "size": 0.25});
        assert!(s.execute("mesh.tet.generate", &tet).is_err(), "open: no volume to fill");
        assert!(
            s.execute("geometry.stitch", &json!({"name": "Small", "surface": "Cube", "tolerance": 0.001})).unwrap()["open_after"].as_u64() > Some(0)
        );
        let r = s.execute("geometry.stitch", &json!({"name": "Gaps", "surface": "Cube", "tolerance": 0.05})).unwrap();
        assert_eq!((r["open_after"].as_u64(), r["merged"].as_u64()), (Some(0), Some(4)), "{r}");
        s.execute("mesh.tet.generate", &tet).unwrap();
        assert!(!s.model.elements().is_empty());
        let file = dir.join("stitched.fieldcraft");
        s.execute("file.save", &json!({"path": file.to_str().unwrap()})).unwrap();
        let mut t = Session::default();
        t.execute("file.open", &json!({"path": file.to_str().unwrap()})).unwrap();
        assert_eq!(t.surfaces[0].mesh, s.surfaces[0].mesh);
        assert!(t.surfaces[0].original.is_some());
        // Undo the mesh and the stitch: the gap is back.
        s.execute("edit.undo", &json!({})).unwrap();
        s.execute("edit.undo", &json!({})).unwrap();
        assert!(s.execute("mesh.tet.generate", &tet).is_err());
        s.execute("edit.redo", &json!({})).unwrap();
        s.execute("mesh.tet.generate", &tet).unwrap();
        let r = s.execute("geometry.stitch.delete", &json!({"name": "Gaps"})).unwrap();
        assert!(r["replaced"].is_array(), "{r}");
        assert!(s.execute("mesh.tet.generate", &tet).is_err());
        std::fs::remove_dir_all(dir).ok();
    }

    #[test]
    fn bad_stitches_change_nothing() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-stitch-bad-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut s = Session::default();
        s.execute("geometry.import.mesh", &json!({"path": gapped_cube(&dir)})).unwrap();
        let journal = s.journal.len();
        for p in [
            json!({}),
            json!({"name": "A", "surface": "Cube"}),
            json!({"name": "A", "surface": "Nope", "tolerance": 0.1}),
            json!({"name": "A", "surface": "Cube", "tolerance": 0}),
            json!({"name": "A", "surface": "Cube", "tolerance": "x"}),
            json!({"name": "1A", "surface": "Cube", "tolerance": 0.1}),
        ] {
            assert!(s.execute("geometry.stitch", &p).is_err(), "{p}");
        }
        assert!(s.execute("geometry.stitch.delete", &json!({"name": "A"})).is_err());
        assert_eq!(s.journal.len(), journal);
        assert!(s.model.stitches().is_empty() && s.surfaces[0].original.is_none());
        std::fs::remove_dir_all(dir).ok();
    }
}
