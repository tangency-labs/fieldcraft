//! Size controls: `mesh.size.set` asks for elements of a size on faces or bodies, in a box or in
//! a sphere; `mesh.size.growth` sets how fast the size grows away from them; `mesh.size.clear`
//! removes controls and `mesh.size.list` lists them with what each covers in the open design.
//! Controls live in the model (`fieldcraft_model::sizes`), are saved with it, and every mesher
//! reads them through [`Sizing`]: the tetrahedral mesher bakes them into its graded size field,
//! the structured hex meshers and the octree take the exact graded field
//! ([`fieldcraft_mesh::size::RegionSize`]), and the block projection (`mesh.hex.project`) divides
//! its block edges by it.

use std::collections::BTreeMap;

use fieldcraft_geometry::{Design, Shape};
use fieldcraft_mesh::size::{LocalSize, Region, RegionSize};
use fieldcraft_mesh::tet::QualityOptions;
use fieldcraft_model::{Entity, GeometryRef, SizeControl, SizeRegion};
use serde::Deserialize;
use serde_json::{Value, json};

use super::opt;
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

/// Most controls one `mesh.size.clear` names.
const MAX_NAMES: usize = 10_000;
/// The chordal tolerance of the tessellation a face control is resolved on, relative to the
/// body's diagonal (as the meshers' default).
const FACE_TOLERANCE: f64 = 1e-3;

pub const SIZE_SET: CommandSpec = CommandSpec {
    id: "mesh.size.set",
    label: "Local Size",
    params: r#"{"name": "Hole", "size": 0.5, and one of "faces": [{"body": "Body1", "face": "F3"}, …] (persistent face names, geometry.bodies; a name without its `#k` number takes every piece), "bodies": ["Body1", …], "box": {"min": [x, y, z], "max": [x, y, z]}, "sphere": {"center": [x, y, z], "radius": r}} — adds the size control `name`, or replaces it: elements of about `size` on those faces, in those bodies, in the box or in the sphere. Every mesher honours it: the size at a point is the smallest of the mesher's size and, over the controls, a control's size plus the growth rate (mesh.size.growth) times the distance to its region, so elements grow away from it gradually. A control finer than the mesher's size refines; a coarser one changes nothing, except on whole bodies, where it is the body's size in place of the mesher's `size` (meshers that mesh several bodies together, mesh.hex.blocks and mesh.cells, take it as a control on the body's surface). Saved with the model; mesh.size.list shows it"#,
    form: Some(super::forms::SIZE_SET),
    kind: CommandKind::Model,
    requires: &[],
    run: size_set,
};

pub const SIZE_GROWTH: CommandSpec = CommandSpec {
    id: "mesh.size.growth",
    label: "Size Growth",
    params: r#"{"rate": 0.3 (0.05..=2; omitted or null: back to the default, 0.3)} — the size field's growth rate: the largest change of element size per unit length away from a size control or an edge seed. mesh.tet.generate's `grading` overrides it for one run"#,
    form: Some(super::forms::SIZE_GROWTH),
    kind: CommandKind::Model,
    requires: &[],
    run: size_growth,
};

pub const SIZE_CLEAR: CommandSpec = CommandSpec {
    id: "mesh.size.clear",
    label: "Clear Local Sizes",
    params: r#"{"names"?: ["Hole", …]} — removes these size controls, or every one"#,
    form: Some(super::forms::SIZE_CLEAR),
    kind: CommandKind::Model,
    requires: &[Requirement::SizeControl],
    run: size_clear,
};

pub const SIZE_LIST: CommandSpec = CommandSpec {
    id: "mesh.size.list",
    label: "Local Sizes",
    params: "{} — the size controls by name, each with its region and size and, for faces and bodies, what it covers in the open design (`stale`: names the design no longer has), and the growth rate",
    form: Some(&[]),
    kind: CommandKind::Model,
    requires: &[],
    run: size_list,
};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BoxParams {
    min: [f64; 3],
    max: [f64; 3],
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct SphereParams {
    center: [f64; 3],
    radius: f64,
}

fn size_set(s: &mut Session, p: &Value) -> Result<Value> {
    let name: String = opt(p, "name")?.ok_or_else(|| EngineError::Params("`name` is required".into()))?;
    let size: f64 = opt(p, "size")?.ok_or_else(|| EngineError::Params("`size` is required: the element size".into()))?;
    let faces: Option<Vec<GeometryRef>> = opt(p, "faces")?;
    let bodies: Option<Vec<String>> = opt(p, "bodies")?;
    let boxed: Option<BoxParams> = opt(p, "box")?;
    let sphere: Option<SphereParams> = opt(p, "sphere")?;
    let given = [faces.is_some(), bodies.is_some(), boxed.is_some(), sphere.is_some()].iter().filter(|g| **g).count();
    if given != 1 {
        return Err(EngineError::Params("give one of `faces`, `bodies`, `box` or `sphere`".into()));
    }
    let region = if let Some(faces) = faces {
        if let Some(r) = faces.iter().find(|r| !matches!(r.entity, Entity::Face(_))) {
            return Err(EngineError::Params(format!("`faces` lists the {r}; name faces as {{body, face}}")));
        }
        SizeRegion::Geometry(faces)
    } else if let Some(bodies) = bodies {
        SizeRegion::Geometry(bodies.into_iter().map(|b| GeometryRef::new(b, Entity::Body)).collect::<fieldcraft_model::Result<_>>()?)
    } else if let Some(b) = boxed {
        SizeRegion::Box { min: b.min, max: b.max }
    } else if let Some(sp) = sphere {
        SizeRegion::Sphere { center: sp.center, radius: sp.radius }
    } else {
        return Err(EngineError::Params("give one of `faces`, `bodies`, `box` or `sphere`".into()));
    };
    let control = SizeControl { region, size };
    // Faces and bodies must name the open design's, when there is one.
    if let (Some(design), SizeRegion::Geometry(refs)) = (s.design.as_deref(), &control.region)
        && let Some(r) = refs.iter().find(|r| covered(design, r).is_ok_and(|c| c.is_empty()))
    {
        return Err(EngineError::Params(format!("the design has no {r} (geometry.bodies lists the faces)")));
    }
    s.model.set_size_control(&name, control.clone())?;
    Ok(report(s, &name, &control))
}

fn size_growth(s: &mut Session, p: &Value) -> Result<Value> {
    let rate: Option<f64> = opt(p, "rate")?;
    s.model.set_size_growth(rate)?;
    Ok(json!({"growth": s.model.size_growth(), "default": !s.model.size_growth_set()}))
}

fn size_clear(s: &mut Session, p: &Value) -> Result<Value> {
    let names: Option<Vec<String>> = opt(p, "names")?;
    let cleared = match names {
        None => {
            let all: Vec<String> = s.model.size_controls().keys().cloned().collect();
            if all.is_empty() {
                return Err(EngineError::Params("the model has no size controls".into()));
            }
            for n in &all {
                s.model.remove_size_control(n)?;
            }
            all.len()
        }
        Some(names) => {
            if names.is_empty() || names.len() > MAX_NAMES {
                return Err(EngineError::Params(format!("`names` must list 1 to {MAX_NAMES} size controls")));
            }
            for n in &names {
                s.model.remove_size_control(n)?;
            }
            names.len()
        }
    };
    Ok(json!({"cleared": cleared, "left": s.model.size_controls().len()}))
}

fn size_list(s: &mut Session, _: &Value) -> Result<Value> {
    let controls: serde_json::Map<String, Value> = s.model.size_controls().iter().map(|(n, c)| (n.clone(), report(s, n, c))).collect();
    Ok(json!({"controls": controls, "growth": s.model.size_growth()}))
}

/// A control's report: its definition and, for faces and bodies, the names each covers in the
/// open design.
fn report(s: &Session, name: &str, c: &SizeControl) -> Value {
    let mut v = json!({"name": name, "size": c.size, "region": c.region});
    if let (SizeRegion::Geometry(refs), Some(design)) = (&c.region, s.design.as_deref()) {
        let on: Vec<Value> = refs.iter().map(|r| json!({"target": r.to_string(), "covers": covered(design, r).unwrap_or_default()})).collect();
        v["stale"] = json!(on.iter().any(|o| o["covers"].as_array().is_none_or(Vec::is_empty)));
        v["on"] = json!(on);
    }
    v
}

/// The faces a face reference selects (its `#k` pieces too), or the body a body reference names.
fn covered(design: &Design, r: &GeometryRef) -> Result<Vec<String>> {
    let Some(body) = design.body(&r.body) else { return Ok(Vec::new()) };
    Ok(match r.entity {
        Entity::Body => vec![body.name.clone()],
        _ => body.faces().iter().filter(|f| selects(r, &body.name, &f.name)).map(|f| f.name.clone()).collect(),
    })
}

fn selects(r: &GeometryRef, body: &str, face: &str) -> bool {
    GeometryRef::new(body.to_string(), Entity::Face(face.to_string())).is_ok_and(|g| r.selects(&g))
}

/// The model's size controls resolved on the open design, as the meshers take them.
#[derive(Debug, Clone, Default)]
pub(crate) struct Sizing {
    /// Local sizes on faces, in boxes and in spheres.
    pub(crate) locals: Vec<LocalSize>,
    /// Per body, the size its body controls set.
    pub(crate) bodies: BTreeMap<String, f64>,
    /// Per body, its boundary surface's triangles (for the meshers that take a body control as
    /// a local size on its surface).
    surfaces: BTreeMap<String, Vec<[[f64; 3]; 3]>>,
    pub(crate) growth: f64,
}

/// The triangles of `shape`'s faces that `select` picks.
fn triangles(shape: &Shape, select: impl Fn(&str) -> bool) -> Result<Vec<[[f64; 3]; 3]>> {
    let tolerance = FACE_TOLERANCE * shape.size()?;
    let mesh = shape.tessellate(tolerance)?;
    let names = shape.faces();
    let point = |i: u32| mesh.positions.get(i as usize).copied();
    Ok(mesh
        .triangles
        .iter()
        .zip(&mesh.tri_face)
        .filter(|(_, f)| names.get(**f as usize).is_some_and(|n| select(&n.name)))
        .filter_map(|(t, _)| Some([point(t[0])?, point(t[1])?, point(t[2])?]))
        .collect())
}

/// Resolve the model's size controls on the open design. Face controls on faces the design does
/// not have, and body controls on bodies it does not have, are left out (mesh.size.list shows
/// them as stale).
pub(crate) fn sizing(s: &Session) -> Result<Sizing> {
    let mut out = Sizing { growth: s.model.size_growth(), ..Sizing::default() };
    for c in s.model.size_controls().values() {
        let region = match &c.region {
            SizeRegion::Box { min, max } => Region::Box { lo: *min, hi: *max },
            SizeRegion::Sphere { center, radius } => Region::Sphere { center: *center, radius: *radius },
            SizeRegion::Geometry(refs) => {
                let Some(design) = s.design.as_deref() else { continue };
                let mut tris = Vec::new();
                for r in refs {
                    let Some(shape) = design.body(&r.body) else { continue };
                    if r.entity == Entity::Body {
                        let size = out.bodies.entry(shape.name.clone()).or_insert(c.size);
                        *size = size.min(c.size);
                        if !out.surfaces.contains_key(&shape.name) {
                            out.surfaces.insert(shape.name.clone(), triangles(shape, |_| true)?);
                        }
                    } else {
                        tris.extend(triangles(shape, |f| selects(r, &shape.name, f))?);
                    }
                }
                if tris.is_empty() {
                    continue;
                }
                Region::Triangles(tris)
            }
        };
        out.locals.push(LocalSize { region, size: c.size });
    }
    Ok(out)
}

impl Sizing {
    /// The size of `body` meshed alone: its body control's, else `size`.
    pub(crate) fn body_size(&self, body: &str, size: f64) -> f64 {
        self.bodies.get(body).copied().unwrap_or(size)
    }

    /// The sizing of several bodies meshed together at one size, so the curves and faces they
    /// share are divided alike: the face, box and sphere controls, and each body control as a
    /// local size on that body's surface.
    pub(crate) fn joint(&self, bodies: &[&str]) -> Sizing {
        let mut locals = self.locals.clone();
        for b in bodies {
            if let (Some(&size), Some(tris)) = (self.bodies.get(*b), self.surfaces.get(*b)) {
                locals.push(LocalSize { region: Region::Triangles(tris.clone()), size });
            }
        }
        Sizing { locals, growth: self.growth, ..Sizing::default() }
    }

    /// The graded field of `body` meshed at `size` (its body control's size in place of it).
    pub(crate) fn field(&self, body: &str, size: f64) -> Result<RegionSize> {
        Ok(RegionSize::new(self.body_size(body, size), self.growth, self.locals.clone())?)
    }

    /// The graded field of the face, box and sphere controls about a size already chosen (body
    /// controls taken into it, or a mesh plan's size), for a mesher that takes a uniform size
    /// and an optional field; `None` when no control refines it.
    pub(crate) fn local_field(&self, size: f64) -> Result<Option<RegionSize>> {
        Ok(Some(RegionSize::new(size, self.growth, self.locals.clone())?).filter(|f| !f.is_uniform()))
    }

    /// The tetrahedral mesher's options for `body`: its body control's size in place of the
    /// global one, the local sizes, and the growth rate unless the command gives `grading`.
    pub(crate) fn tet_options(&self, body: &str, p: &Value, o: &QualityOptions) -> QualityOptions {
        QualityOptions {
            size: self.bodies.get(body).copied().or(o.size),
            locals: self.locals.clone(),
            grading: if p.get("grading").is_some() { o.grading } else { self.growth },
            ..o.clone()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plate() -> (Session, std::path::PathBuf) {
        super::super::datum::tests::plate()
    }

    /// The name of the plate's face that lies in the plane x = `x`.
    fn face_at_x(s: &Session, x: f64) -> String {
        let body = s.design.as_ref().unwrap().body("Body1").unwrap();
        let mesh = body.tessellate(0.1).unwrap();
        let faces = body.faces();
        (0..faces.len())
            .find(|&f| {
                let mut pts = mesh.triangles.iter().zip(&mesh.tri_face).filter(|(_, g)| **g as usize == f).flat_map(|(t, _)| t.iter());
                pts.clone().next().is_some() && pts.all(|&v| (mesh.positions[v as usize][0] - x).abs() < 1e-6)
            })
            .map(|f| faces[f].name.clone())
            .unwrap()
    }

    /// Every element edge of the model's tet4 or hex8 elements, once: its two end points, in
    /// millimetres as the design is.
    fn edges(s: &Session) -> Vec<([f64; 3], [f64; 3])> {
        const TET: &[[usize; 2]] = &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];
        const HEX: &[[usize; 2]] = &[[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]];
        let mut seen = std::collections::BTreeSet::new();
        let mut out = Vec::new();
        for e in s.model.elements().values() {
            let pairs = if e.nodes.len() == 4 { TET } else { HEX };
            for [a, b] in pairs {
                let (a, b) = (e.nodes[*a], e.nodes[*b]);
                if seen.insert((a.min(b), a.max(b))) {
                    let mm = crate::units::to_geometry_point;
                    out.push((mm(s.model.nodes()[&a].xyz), mm(s.model.nodes()[&b].xyz)));
                }
            }
        }
        out
    }

    fn length((a, b): &([f64; 3], [f64; 3])) -> f64 {
        ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
    }

    fn mean(v: &[f64]) -> f64 {
        assert!(!v.is_empty());
        v.iter().sum::<f64>() / v.len() as f64
    }

    #[test]
    fn controls_are_set_listed_cleared_saved_and_undone() {
        let (mut s, dir) = plate();
        let end = face_at_x(&s, 40.0);
        let r = s.execute("mesh.size.set", &json!({"name": "End", "size": 1.5, "faces": [{"body": "Body1", "face": end}]})).unwrap();
        assert_eq!(r["on"][0]["covers"], json!([end]), "{r}");
        assert_eq!(r["stale"], false);
        s.execute("mesh.size.set", &json!({"name": "Zone", "size": 2, "box": {"min": [0, 0, 0], "max": [5, 5, 10]}})).unwrap();
        s.execute("mesh.size.set", &json!({"name": "Zone", "size": 2, "sphere": {"center": [0, 0, 0], "radius": 3}})).unwrap();
        s.execute("mesh.size.set", &json!({"name": "Whole", "size": 4, "bodies": ["Body1"]})).unwrap();
        s.execute("mesh.size.growth", &json!({"rate": 0.5})).unwrap();
        let list = s.execute("mesh.size.list", &json!({})).unwrap();
        assert_eq!(list["controls"].as_object().unwrap().len(), 3, "{list}");
        assert_eq!(list["controls"]["Zone"]["region"]["sphere"]["radius"], 3.0);
        assert_eq!(list["growth"], 0.5);

        let file = dir.join("sized.fieldcraft");
        s.execute("file.save", &json!({"path": file.to_str().unwrap()})).unwrap();
        let mut back = Session::default();
        back.execute("file.open", &json!({"path": file.to_str().unwrap()})).unwrap();
        assert_eq!(back.model.size_controls(), s.model.size_controls());
        assert_eq!(back.model.size_growth(), 0.5);

        s.execute("mesh.size.clear", &json!({"names": ["Zone"]})).unwrap();
        assert_eq!(s.model.size_controls().len(), 2);
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(s.model.size_controls().len(), 3);
        assert!(s.execute("mesh.size.clear", &json!({"names": ["Nope"]})).is_err());
        assert_eq!(s.execute("mesh.size.clear", &json!({})).unwrap()["cleared"], 3);
        assert!(s.execute("mesh.size.clear", &json!({})).is_err(), "nothing to clear");
        s.execute("mesh.size.growth", &json!({})).unwrap();
        assert_eq!(s.model.size_growth(), fieldcraft_model::sizes::DEFAULT_GROWTH);
    }

    #[test]
    fn hostile_params_never_panic() {
        let (mut s, _) = plate();
        for p in [
            json!({}),
            json!({"name": "A"}),
            json!({"name": "A", "size": 1}),
            json!({"name": "A", "size": 0, "bodies": ["Body1"]}),
            json!({"name": "A", "size": -1, "bodies": ["Body1"]}),
            json!({"name": "A", "size": 1e300, "bodies": ["Body1"]}),
            json!({"name": "A", "size": "big", "bodies": ["Body1"]}),
            json!({"name": "A", "size": 1, "bodies": ["Nope"]}),
            json!({"name": "A", "size": 1, "bodies": []}),
            json!({"name": "A", "size": 1, "faces": [{"body": "Body1", "face": "Nope"}]}),
            json!({"name": "A", "size": 1, "faces": [{"body": "Body1", "edge": "E"}]}),
            json!({"name": "A", "size": 1, "bodies": ["Body1"], "box": {"min": [0, 0, 0], "max": [1, 1, 1]}}),
            json!({"name": "A", "size": 1, "box": {"min": [1, 0, 0], "max": [0, 1, 1]}}),
            json!({"name": "A", "size": 1, "box": {"min": [0, 0], "max": [1, 1, 1]}}),
            json!({"name": "A", "size": 1, "sphere": {"center": [0, 0, 0], "radius": 0}}),
            json!({"name": "A", "size": 1, "sphere": {"center": [0, 0, 0], "radius": 1, "x": 1}}),
            json!({"name": "1 bad", "size": 1, "bodies": ["Body1"]}),
        ] {
            assert!(s.execute("mesh.size.set", &p).is_err(), "{p}");
        }
        for p in [json!({"rate": 0}), json!({"rate": 5}), json!({"rate": "x"})] {
            assert!(s.execute("mesh.size.growth", &p).is_err(), "{p}");
        }
        assert!(s.execute("mesh.size.clear", &json!({"names": []})).is_err());
        assert!(s.model.size_controls().is_empty());
    }

    /// A size box halves the tetrahedra inside it; a face control gives its face the size it
    /// asks for, and the elements grow away from it no faster than the growth rate allows.
    #[test]
    fn tet_meshes_follow_size_controls() {
        let (mut s, _) = plate();
        s.execute("mesh.tet.generate", &json!({"body": "Body1", "size": 6, "tolerance": 0.2})).unwrap();
        let inside = |e: &&([f64; 3], [f64; 3])| [e.0, e.1].iter().all(|p| p[0] >= 30.0 && p[1] >= 0.0 && p[1] <= 10.0);
        let before = mean(&edges(&s).iter().filter(inside).map(length).collect::<Vec<_>>());
        s.execute("mesh.size.set", &json!({"name": "Zone", "size": 3, "box": {"min": [30, -1, -1], "max": [41, 10, 11]}})).unwrap();
        s.execute("mesh.tet.generate", &json!({"body": "Body1", "size": 6, "tolerance": 0.2})).unwrap();
        let after = mean(&edges(&s).iter().filter(inside).map(length).collect::<Vec<_>>());
        assert!((after / before - 0.5).abs() < 0.15, "box: {before} -> {after}");

        s.execute("mesh.size.clear", &json!({})).unwrap();
        let end = face_at_x(&s, 40.0);
        s.execute("mesh.size.set", &json!({"name": "End", "size": 1.5, "faces": [{"body": "Body1", "face": end}]})).unwrap();
        s.execute("mesh.size.growth", &json!({"rate": 0.5})).unwrap();
        s.execute("mesh.tet.generate", &json!({"body": "Body1", "size": 6, "tolerance": 0.2})).unwrap();
        let all = edges(&s);
        let on_face: Vec<f64> = all.iter().filter(|e| (e.0[0] - 40.0).abs() < 1e-6 && (e.1[0] - 40.0).abs() < 1e-6).map(length).collect();
        let m = mean(&on_face);
        assert!((1.2..1.95).contains(&m), "on the face: {m}");
        // Slabs away from the face: the mean edge stays near 1.5 + 0.5 d (and below 6).
        for d in [3.0, 6.0, 9.0] {
            let slab: Vec<f64> = all.iter().filter(|e| (40.0 - 0.5 * (e.0[0] + e.1[0]) - d).abs() < 1.0).map(length).collect();
            let want = (1.5f64 + 0.5 * d).min(6.0);
            let m = mean(&slab);
            assert!(m < want * 1.3 && m > want * 0.5, "{d} mm from the face: mean edge {m}, field {want}");
        }
    }

    /// A face control sets the hex divisions of the face's edges: the swept plate's end face gets
    /// elements of the size asked for.
    #[test]
    fn hex_meshes_follow_size_controls() {
        let (mut s, _) = plate();
        let end = face_at_x(&s, 40.0);
        s.execute("mesh.hex.generate", &json!({"body": "Body1", "size": 5})).unwrap();
        let on_end = |s: &Session| -> Vec<f64> {
            edges(s).iter().filter(|e| (e.0[0] - 40.0).abs() < 1e-6 && (e.1[0] - 40.0).abs() < 1e-6).map(length).collect()
        };
        let coarse = on_end(&s).into_iter().fold(0.0, f64::max);
        assert!(coarse > 3.0, "{coarse}");
        s.execute("mesh.size.set", &json!({"name": "End", "size": 1.25, "faces": [{"body": "Body1", "face": end}]})).unwrap();
        s.execute("mesh.hex.generate", &json!({"body": "Body1", "size": 5})).unwrap();
        let fine = on_end(&s);
        assert!(fine.iter().all(|&l| l < 1.25 * 1.3), "{fine:?}");
        assert!(mean(&fine) > 1.25 * 0.7, "{fine:?}");
    }

    /// The octree's cells on a face with a size control reach the level its size needs.
    #[test]
    fn octree_meshes_follow_size_controls() {
        let (mut s, _) = plate();
        let end = face_at_x(&s, 40.0);
        s.execute("mesh.size.set", &json!({"name": "End", "size": 2, "faces": [{"body": "Body1", "face": end}]})).unwrap();
        s.execute("mesh.hex.octree", &json!({"body": "Body1", "size": 6})).unwrap();
        let on_end: Vec<f64> = edges(&s).iter().filter(|e| (e.0[0] - 40.0).abs() < 1e-3 && (e.1[0] - 40.0).abs() < 1e-3).map(length).collect();
        assert!(!on_end.is_empty() && mean(&on_end) < 2.0 * 1.3, "{on_end:?}");
    }
}
