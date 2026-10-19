//! `model.boundary_condition.symmetry`: a symmetry condition on a datum plane, so half of a
//! symmetric part can be solved as a half.
//!
//! The condition holds the displacement normal to the plane, and the rotations about the two
//! axes in it, at zero (`type: "symmetry"`): on a plane normal to X, Y or Z those are global
//! degrees of freedom; on any other plane they are linear constraints along its normal. Its set is given, or made of the faces of the open
//! design's bodies that lie on the plane (the cut faces of a partition by that datum, `partition:
//! <name>:a|b`), so it follows the geometry through a remesh.

use fieldcraft_geometry::FaceSurface;
use fieldcraft_model::{BoundaryCondition, GeometryRef, Set, Symmetry, Target};
use serde_json::{Value, json};

use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

/// A face normal this close to the plane's, relative to its length, is parallel.
const PARALLEL: f64 = 1e-9; // unit-factors: relative to the normal's length

pub const SYMMETRY: CommandSpec = CommandSpec {
    id: "model.boundary_condition.symmetry",
    label: "Symmetry Plane",
    params: r#"{"name": "Sym", "datum": "Mid" (a datum plane, any orientation), "set"?: "cut" (default: a new set of the same name holding every face of the open design's bodies on the plane, e.g. a partition's cut faces), "bodies"?: ["Part/a"] (only their faces)} — a symmetry boundary condition on the plane: zero displacement normal to it and zero rotation about the axes in it, so half a symmetric part (loaded symmetrically, with half the load) is solved as a half. Add it to the steps (model.step.activate). Reports the plane as "x", "y", "z" or its unit normal"#,
    form: Some(super::forms::SYMMETRY),
    kind: CommandKind::Model,
    requires: &[Requirement::DatumPlane],
    run: symmetry,
};

/// A face lies on the plane within this fraction of its body's size.
const ON_PLANE: f64 = 1e-6; // unit-factors: relative to the body's size

fn symmetry(s: &mut Session, p: &Value) -> Result<Value> {
    let name = super::name(p)?.to_string();
    let datum: String = super::req(p, "datum")?;
    let (point, normal) = super::datum::datum_plane_of(&s.model, &datum)?;
    let unusable = || EngineError::Params(format!("datum `{datum}` has no usable normal ({normal:?})"));
    let plane = Symmetry::from_normal(normal).ok_or_else(unusable)?;
    let n = plane.unit_normal().ok_or_else(unusable)?;
    let along = |x: [f64; 3]| (0..3).map(|k| x[k] * n[k]).sum::<f64>();
    let set = match p.get("set") {
        Some(_) => super::req::<String>(p, "set")?,
        None => {
            let only: Option<Vec<String>> = super::opt(p, "bodies")?;
            let design =
                s.design.as_deref().ok_or_else(|| EngineError::Params("open a design to find the faces on the plane, or give `set`".into()))?;
            let mut faces = Vec::new();
            for b in &design.bodies {
                if only.as_ref().is_some_and(|o| !o.contains(&b.name)) {
                    continue;
                }
                let tol = ON_PLANE * b.size()?.max(1.0);
                for (f, surface) in b.faces().iter().zip(b.face_surfaces()?) {
                    if let Some(FaceSurface::Plane { normal: m, point: q }) = surface {
                        let length = m.iter().map(|c| c * c).sum::<f64>().sqrt();
                        let parallel = (along(m).abs() - length).abs() < PARALLEL * length;
                        if parallel && (along(q) - along(point)).abs() < tol {
                            faces.push(GeometryRef::face(b.name.as_str(), f.name.as_str())?);
                        }
                    }
                }
            }
            if faces.is_empty() {
                return Err(EngineError::Params(format!("no face of the design lies on datum `{datum}`: partition by it first, or give `set`")));
            }
            s.model.add_set(&name, Set::of_geometry(faces))?;
            name.clone()
        }
    };
    s.model.add_boundary_condition(&name, BoundaryCondition::Symmetry { target: Target::Set(set.clone()), plane })?;
    let faces = s.model.set(&set).map_or(0, |x| x.geometry.len());
    let reported = match plane.axis() {
        Some(k) => json!(["x", "y", "z"].get(k).copied().unwrap_or("x")),
        None => json!(n),
    };
    Ok(json!({"name": name, "set": set, "plane": reported, "faces": faces}))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A 40 × 10 × 10 bar along x (0..40, 0..10, 0..10), the datum `Mid` at y = 5.
    fn bar(half: bool) -> (crate::tests::TempDir, Session) {
        let commands = json!({"commands": [
            {"command": "solid.box", "params": {"length": 40, "width": 10, "height": 10, "body_name": "Bar"}},
        ]});
        let bytes = fieldcraft_geometry::script_design(&commands.to_string()).unwrap();
        let (dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        s.execute("geometry.datum", &json!({"name": "Mid", "by": {"type": "offset", "plane": "zx", "offset": 5}})).unwrap();
        if half {
            s.execute("geometry.partition", &json!({"name": "Half", "datum": "Mid"})).unwrap();
        }
        (dir, s)
    }

    fn face(s: &Session, body: &str, point: [f64; 3]) -> String {
        let shape = s.design.as_ref().unwrap().body(body).unwrap();
        let f = shape.nearest_faces(&[point], 0.01).unwrap()[0].unwrap();
        shape.faces()[f as usize].name.clone()
    }

    /// Clamp the x = 0 end and push the x = 40 end down (−z) with `force`, solve, and give the
    /// largest downward tip deflection and the clamp's reaction.
    fn bend(s: &mut Session, body: &str, force: f64, extra: &[&str]) -> (f64, f64) {
        let (root, tip) = (face(s, body, [0.0, 2.5, 5.0]), face(s, body, [40.0, 2.5, 5.0]));
        let mut bcs = vec!["Clamp".to_string()];
        bcs.extend(extra.iter().map(|x| x.to_string()));
        let script = json!([
            {"command": "model.set.create", "params": {"name": "all", "geometry": [{"body": body}]}},
            {"command": "model.set.create", "params": {"name": "root", "geometry": [{"body": body, "face": root}]}},
            {"command": "model.set.create", "params": {"name": "tip", "geometry": [{"body": body, "face": tip}]}},
            {"command": "model.material.create", "params": {"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}}},
            {"command": "model.section.create", "params": {"name": "Bar", "set": "all", "material": "Steel", "kind": {"type": "solid"}}},
            {"command": "model.boundary_condition.create", "params": {"name": "Clamp", "type": "fixed", "target": {"set": "root"}}},
            {"command": "model.load.create", "params": {"name": "Push", "type": "force", "target": {"set": "tip"}, "vector": [0, 0, -force]}},
            {"command": "model.step.create", "params": {"name": "Load", "procedure": {"type": "static"}, "boundary_conditions": bcs, "loads": ["Push"]}},
        ]);
        s.run_script(&script).unwrap();
        s.execute("solve.static", &json!({})).unwrap();
        let tip = s.model.resolve_set("tip").unwrap().nodes;
        let probe = s.execute("solve.probe", &json!({"step": "Load", "nodes": tip})).unwrap();
        let w = probe["nodes"].as_array().unwrap().iter().map(|n| -n["displacement"][2].as_f64().unwrap()).fold(f64::MIN, f64::max);
        let r = s.execute("results.reactions", &json!({"step": "Load", "sets": ["root"]})).unwrap();
        (w, r["sets"][0]["force"][2].as_f64().unwrap_or(f64::NAN))
    }

    /// Half the bar, cut on its symmetry plane y = 5, with the symmetry condition and half the
    /// load, bends as the whole bar does: the same deflection, half the reaction.
    #[test]
    fn a_half_model_with_symmetry_matches_the_whole() {
        let (_d1, mut whole) = bar(false);
        whole.execute("mesh.hex.mapped", &json!({"body": "Bar", "size": 2.5, "order": 2})).unwrap();
        let full = bend(&mut whole, "Bar", 1000.0, &[]);

        let (_d2, mut half) = bar(true);
        half.execute("mesh.hex.mapped", &json!({"body": "Bar/a", "size": 2.5, "order": 2})).unwrap();
        let r = half.execute("model.boundary_condition.symmetry", &json!({"name": "Sym", "datum": "Mid", "bodies": ["Bar/a"]})).unwrap();
        assert_eq!(r, json!({"name": "Sym", "set": "Sym", "plane": "y", "faces": 1}));
        assert_eq!(half.model.set("Sym").unwrap().geometry[0].entity, fieldcraft_model::geometry::Entity::Face("partition:Half:a".into()));
        let halved = bend(&mut half, "Bar/a", 500.0, &["Sym"]);
        assert!((halved.0 - full.0).abs() < 1e-9 * full.0, "deflection {} vs {}", halved.0, full.0);
        assert!((2.0 * halved.1 - full.1).abs() < 1e-9 * full.1.abs(), "reaction {} vs {}", halved.1, full.1);
    }

    /// A 20 × 20 × 2 plate (0..20 in x and y), cut on its diagonal plane x = y (datum `Diag`,
    /// normal (1, −1, 0)) into `Plate/a` and `Plate/b`, the half below the diagonal (y < x) meshed
    /// with quadratic tets.
    fn plate() -> (crate::tests::TempDir, Session, String) {
        let commands = json!({"commands": [
            {"command": "solid.box", "params": {"length": 20, "width": 20, "height": 2, "body_name": "Plate"}},
        ]});
        let bytes = fieldcraft_geometry::script_design(&commands.to_string()).unwrap();
        let (dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        s.execute("geometry.datum", &json!({"name": "Diag", "by": {"type": "point_normal", "point": [0, 0, 0], "normal": [1, -1, 0]}})).unwrap();
        s.execute("geometry.partition", &json!({"name": "Cut", "datum": "Diag"})).unwrap();
        let design = s.design.as_ref().unwrap();
        let below = design
            .bodies
            .iter()
            .find(|b| {
                let c = b.measure().unwrap().centroid;
                c[0] > c[1]
            })
            .unwrap()
            .name
            .clone();
        s.execute("mesh.tet.generate", &json!({"body": below, "size": 2, "order": "quadratic"})).unwrap();
        (dir, s, below)
    }

    /// Clamp the faces at y = 0 (and x = 0 when `whole`), press the top faces, solve; the
    /// displacement of every node of the half below the diagonal, by label.
    fn press(s: &mut Session, below: &str, whole: bool, extra: &[&str]) -> std::collections::BTreeMap<u64, [f64; 3]> {
        let bodies: Vec<String> = s.design.as_ref().unwrap().bodies.iter().map(|b| b.name.clone()).collect();
        let mut clamp = vec![json!({"body": below, "face": face(s, below, [10.0, 0.0, 1.0])})];
        let mut top = vec![json!({"body": below, "face": face(s, below, [15.0, 5.0, 2.0])})];
        if whole {
            let above = bodies.iter().find(|b| b.as_str() != below).unwrap();
            clamp.push(json!({"body": above, "face": face(s, above, [0.0, 10.0, 1.0])}));
            top.push(json!({"body": above, "face": face(s, above, [5.0, 15.0, 2.0])}));
        }
        let all: Vec<Value> = bodies.iter().filter(|b| whole || b.as_str() == below).map(|b| json!({"body": b})).collect();
        let mut bcs = vec!["Clamp".to_string()];
        bcs.extend(extra.iter().map(|x| x.to_string()));
        let script = json!([
            {"command": "model.set.create", "params": {"name": "all", "geometry": all}},
            {"command": "model.set.create", "params": {"name": "clamped", "geometry": clamp}},
            {"command": "model.set.create", "params": {"name": "top", "geometry": top}},
            {"command": "model.set.create", "params": {"name": "half", "geometry": [{"body": below}]}},
            {"command": "model.material.create", "params": {"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}}},
            {"command": "model.section.create", "params": {"name": "Plate", "set": "all", "material": "Steel", "kind": {"type": "solid"}}},
            {"command": "model.boundary_condition.create", "params": {"name": "Clamp", "type": "fixed", "target": {"set": "clamped"}}},
            {"command": "model.load.create", "params": {"name": "Press", "type": "pressure", "set": "top", "magnitude": 0.1}},
            {"command": "model.step.create", "params": {"name": "Load", "procedure": {"type": "static"}, "boundary_conditions": bcs, "loads": ["Press"]}},
        ]);
        s.run_script(&script).unwrap();
        s.execute("solve.static", &json!({})).unwrap();
        let nodes = s.model.resolve_set("half").unwrap().nodes;
        let probe = s.execute("solve.probe", &json!({"step": "Load", "nodes": nodes})).unwrap();
        probe["nodes"]
            .as_array()
            .unwrap()
            .iter()
            .map(|n| (n["node"].as_u64().unwrap(), [0, 1, 2].map(|k| n["displacement"][k].as_f64().unwrap())))
            .collect()
    }

    /// The half of the plate below its diagonal, with the inclined symmetry plane (found on the
    /// cut face) and its half of the pressure, deflects as the whole plate, meshed as that half
    /// mirrored across the diagonal (conforming there, so the whole is the same mesh twice).
    #[test]
    fn a_half_cut_on_an_inclined_plane_matches_the_whole() {
        let (_d1, mut whole, below) = plate();
        let r = whole.execute("mesh.mirror", &json!({"plane": "Diag", "bodies": [below]})).unwrap();
        let on_plane = r["on_plane"].as_u64().unwrap();
        let q = whole.execute("mesh.quality", &json!({})).unwrap();
        assert_eq!(q["valid"].as_u64().unwrap() as usize, whole.model.elements().len(), "{q}");
        assert!((q["volume"].as_f64().unwrap() - 800.0).abs() < 1e-9 * 800.0, "{q}");
        let full = press(&mut whole, &below, true, &[]);

        let (_d2, mut half, below) = plate();
        let nodes = half.model.nodes().len() as u64;
        assert_eq!(whole.model.nodes().len() as u64, 2 * nodes - on_plane, "the mirror shares the plane's nodes");
        let r = half.execute("model.boundary_condition.symmetry", &json!({"name": "Sym", "datum": "Diag", "bodies": [below]})).unwrap();
        let h = 1.0 / 2f64.sqrt();
        assert_eq!(r["faces"], json!(1), "{r}");
        let n: Vec<f64> = r["plane"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
        assert!((n[0] - h).abs() < 1e-15 && (n[1] + h).abs() < 1e-15 && n[2] == 0.0, "{r}");
        assert_eq!(half.model.resolve_set("Sym").unwrap().nodes.len() as u64, on_plane);
        let halved = press(&mut half, &below, false, &["Sym"]);

        assert_eq!(full.len(), halved.len());
        let size = full.values().fold(0.0f64, |m, u| m.max(u.iter().map(|c| c * c).sum::<f64>().sqrt()));
        let worst = halved.iter().map(|(l, u)| (0..3).map(|k| (u[k] - full[l][k]).powi(2)).sum::<f64>().sqrt()).fold(0.0f64, f64::max);
        eprintln!("half vs whole plate on the inclined plane: largest deflection {size:.6e}, difference {:.1e}", worst / size);
        assert!(size > 1e-4 && worst < 1e-9 * size, "{worst} of {size}");
        // On the plane the nodes move along it: n · u = 0.
        for l in half.model.resolve_set("Sym").unwrap().nodes {
            let u = halved[&u64::from(l)];
            assert!((u[0] - u[1]).abs() * h < 1e-12 * size, "node {l}: {u:?}");
        }
    }

    #[test]
    fn hostile_params_change_nothing() {
        let (_d, mut s) = bar(true);
        s.execute("geometry.datum", &json!({"name": "Out", "by": {"type": "offset", "plane": "xy", "offset": 100}})).unwrap();
        s.execute("geometry.datum", &json!({"name": "Tilt", "by": {"type": "point_normal", "point": [100, 0, 0], "normal": [1, 1, 0]}})).unwrap();
        s.execute("geometry.datum.axis", &json!({"name": "Ax", "by": {"type": "principal", "axis": "x"}})).unwrap();
        let model = s.model.clone();
        for p in [
            json!({"name": "S", "datum": "Out"}),
            json!({"name": "S", "datum": "Tilt"}),
            json!({"name": "S", "datum": "Ax"}),
            json!({"name": "S", "datum": "Nope"}),
            json!({"name": "S", "datum": "Mid", "set": "nope"}),
            json!({"name": "S", "datum": "Mid", "bodies": ["Nope"]}),
            json!({"name": "1S", "datum": "Mid"}),
            json!({"datum": "Mid"}),
        ] {
            assert!(s.execute("model.boundary_condition.symmetry", &p).is_err(), "{p}");
        }
        assert_eq!(s.model, model);
    }
}
