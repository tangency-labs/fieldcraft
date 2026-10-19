//! Decomposition proposals (#439): for a body with features the plan cannot map, the partition
//! cuts that box those features into a cell of their own, so the rest of the body is cells the
//! sweep meshes and only the boxed cell is left to tetrahedra (joined to the hex cells by
//! pyramids, `mesh.cells`).
//!
//! The region is the bounding box of the unmapped features' faces (those not already prismatic
//! along the sweep), in the frame of the sweep direction and two directions square to it. Each
//! side of the box becomes a cut one element size outside it, where that plane lies inside the
//! body, clear of its planar faces and of its curved faces (a cut through a curved face does not
//! close in the partitions' tessellation); a side with no such plane is left open and the tet
//! cell reaches the body's boundary there. Decomposing a part into sweepable volumes and a
//! remainder for a free mesher is the practice of T. D. Blacker, *The Cooper Tool*, Proc. 5th
//! International Meshing Roundtable, 1996 (docs/meshing/literature.md, B11). The proposal is a
//! reply, not a change: the agent applies the cuts (`commands`) and meshes the cells.

use fieldcraft_geometry::{Design, FaceSurface};
use fieldcraft_model::MeshPlan;
use serde_json::{Value, json};

use super::{AXES, Clearance, V3, cross, dot, faces_prismatic, parallel, plan_shape, unit};
use crate::Result;

/// The decomposition proposal for `plan`, or `None` when every feature is mapped (or the plan
/// has no sweep direction to make cells along).
pub(super) fn proposal(design: &Design, plan: &MeshPlan) -> Result<Option<Value>> {
    let Some(d) = plan.direction else { return Ok(None) };
    if plan.unmapped.is_empty() {
        return Ok(None);
    }
    let Some(shape) = plan_shape(design, &plan.body) else { return Ok(None) };
    let source = design.source_of(&plan.body).unwrap_or(&plan.body);
    let features: Vec<_> = design.features(Some(source))?.into_iter().next().map(|b| b.features).unwrap_or_default();
    let surfaces = shape.face_surfaces()?;
    // The unmapped features whose faces the sweep cannot mesh.
    let lumps: Vec<_> = features
        .iter()
        .filter(|f| plan.unmapped.iter().any(|u| u.id == f.id))
        .filter(|f| !f.faces.is_empty() && !faces_prismatic(shape, &surfaces, &f.faces, d))
        .collect();
    if lumps.is_empty() {
        return Ok(None);
    }
    let faces: Vec<u32> = lumps.iter().flat_map(|f| &f.faces).filter_map(|n| shape.face_index(n)).collect();
    let m = shape.measure()?;
    let diagonal = (0..3).map(|k| (m.bbox_max[k] - m.bbox_min[k]).powi(2)).sum::<f64>().sqrt();
    let mesh = shape.tessellate(plan.tolerance.unwrap_or(1e-3 * diagonal))?;
    // The frame: the sweep direction and two directions square to it (axes when it is one).
    let side = AXES.iter().copied().filter(|a| !parallel(*a, d)).min_by(|a, b| dot(*a, d).abs().total_cmp(&dot(*b, d).abs())).unwrap_or(AXES[0]);
    let e1 = unit(cross(d, side)).unwrap_or(side);
    let e2 = unit(cross(d, e1)).unwrap_or(side);
    let frame = [d, e1, e2].map(|a| AXES.iter().copied().find(|x| parallel(*x, a)).unwrap_or(a));
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    let mut curved: Vec<Vec<V3>> = vec![Vec::new(); surfaces.len()];
    for (t, &f) in mesh.triangles.iter().zip(&mesh.tri_face) {
        let pts: Vec<V3> = t.iter().filter_map(|&i| mesh.positions.get(i as usize).copied()).collect();
        if faces.contains(&f) {
            for p in &pts {
                for (k, a) in frame.iter().enumerate() {
                    lo[k] = lo[k].min(dot(*p, *a));
                    hi[k] = hi[k].max(dot(*p, *a));
                }
            }
        }
        if !matches!(surfaces.get(f as usize), Some(Some(FaceSurface::Plane { .. })))
            && let Some(c) = curved.get_mut(f as usize)
        {
            c.extend(pts);
        }
    }
    if lo.iter().chain(&hi).any(|v| !v.is_finite()) {
        return Ok(None);
    }
    curved.retain(|p| !p.is_empty());
    let planes = surfaces
        .iter()
        .zip(shape.faces())
        .filter_map(|(s, f)| match s {
            Some(FaceSurface::Plane { normal, point }) => Some((*normal, dot(*normal, *point), f.name.as_str())),
            _ => None,
        })
        .collect();
    let corners = (0..8)
        .map(|i| {
            [
                if i & 1 == 0 { m.bbox_min[0] } else { m.bbox_max[0] },
                if i & 2 == 0 { m.bbox_min[1] } else { m.bbox_max[1] },
                if i & 4 == 0 { m.bbox_min[2] } else { m.bbox_max[2] },
            ]
        })
        .collect();
    let body = Clearance { corners, planes, curved, gap: 0.5 * plan.size };
    let names: Vec<&str> = lumps.iter().map(|f| f.name.as_str()).collect();
    let mut cuts = Vec::new();
    let mut open = Vec::new();
    let mut region = [[0.0; 3]; 2];
    for (k, n) in frame.iter().enumerate() {
        let (blo, bhi) = body.extent(*n);
        region[0][k] = blo;
        region[1][k] = bhi;
        for (end, t) in [(0, lo[k] - plan.size), (1, hi[k] + plan.size)] {
            let which = if end == 0 { "below" } else { "above" };
            let why = if !body.inside(*n, t) {
                Some("the body ends there".to_string())
            } else if let Some((_, face)) = body.faces_near(*n, t).next() {
                Some(format!("it would lie on face `{face}`"))
            } else if body.crosses_curved(*n, t) {
                Some("it would cross a curved face".into())
            } else {
                None
            };
            match why {
                Some(why) => open.push(json!({"normal": n, "side": which, "reason": why})),
                None => {
                    region[end][k] = t;
                    let point = n.map(|c| c * t);
                    cuts.push(json!({"point": point, "normal": n, "reason": format!("boxes {} ({which} along {n:?})", names.join(", "))}));
                }
            }
        }
    }
    // The commands that make the cuts and mesh the cells (the boxed cell falls back to tet).
    let mut commands = Vec::new();
    for (i, c) in cuts.iter().enumerate() {
        let name = format!("{}-split{}", plan.body, i + 1);
        commands.push(json!({"command": "geometry.datum", "params": {"name": name, "by": {"type": "point_normal", "point": c["point"], "normal": c["normal"]}}}));
        commands.push(json!({"command": "geometry.partition", "params": {"name": name, "datum": name}}));
    }
    commands.push(json!({"command": "mesh.cells", "params": {"default": "auto", "size": plan.size, "transition": "pyramid"}}));
    let box_volume: f64 = (0..3).map(|k| region[1][k] - region[0][k]).product();
    let share = if m.volume.abs() > 0.0 { (box_volume / m.volume.abs()).min(1.0) } else { 1.0 };
    Ok(Some(json!({
        "features": names,
        "cuts": cuts,
        "open_sides": open,
        "tet_region": {"frame": frame, "min": region[0], "max": region[1]},
        "tet_share_at_most": share,
        "commands": commands,
        "note": "cut the body by these planes; the cell inside `tet_region` is meshed with tet, every other cell is a box or a slice of the sweep's cross-section the hex meshers take, and mesh.cells joins the tet cell to them through pyramids",
    })))
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use serde_json::{Value, json};

    use crate::Session;

    /// A 80 × 40 × 12 plate with a ball (radius 8) joined on its top: the sweep takes the plate,
    /// not the dome.
    fn plate_with_dome() -> (crate::tests::TempDir, Session) {
        let commands = json!({"commands": [
            {"command": "solid.box", "params": {"length": 80, "width": 40, "height": 12, "body_name": "Plate"}},
            {"command": "solid.sphere", "params": {"center": [60, 20, 12], "radius": 8, "body_name": "Ball", "operation": "new"}},
            {"command": "solid.combine", "params": {"target": "Plate", "tools": ["Ball"], "operation": "join"}},
        ]});
        let bytes = fieldcraft_geometry::script_design(&commands.to_string()).unwrap();
        let (dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        (dir, s)
    }

    /// Fix the end x = 0, pull the end x = 80 with `force` (N), solve, and return the mean
    /// displacement of the loaded end along `k`.
    fn pull(s: &mut Session, force: [f64; 3], k: usize) -> f64 {
        let m = &s.model;
        // The model is stored in metres; the plate is drawn in millimetres.
        let mm = crate::units::to_geometry_length;
        let at = |x: f64| -> Vec<u32> { m.nodes().iter().filter(|(_, n)| (mm(n.xyz[0]) - x).abs() < 1e-6).map(|(l, _)| *l).collect() };
        let (root, tip) = (at(0.0), at(80.0));
        // The loaded end's element faces, so the force is shared by consistent weights.
        let end: Vec<Value> = m
            .elements()
            .iter()
            .flat_map(|(&e, el)| (1..=el.kind.faces().len() as u8).map(move |f| (e, f)))
            .filter(|&(e, f)| {
                m.face_nodes(fieldcraft_model::ElementFace { element: e, face: f })
                    .is_some_and(|ns| ns.iter().all(|n| (mm(m.node(*n).unwrap().xyz[0]) - 80.0).abs() < 1e-6))
            })
            .map(|(e, f)| json!({"element": e, "face": f}))
            .collect();
        let script = json!([
            {"command": "model.set.create", "params": {"name": "all", "elements": m.elements().keys().collect::<Vec<_>>()}},
            {"command": "model.set.create", "params": {"name": "root", "nodes": root}},
            {"command": "model.set.create", "params": {"name": "tip", "faces": end}},
            {"command": "model.material.create", "params": {"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}}},
            {"command": "model.section.create", "params": {"name": "S", "set": "all", "material": "Steel", "kind": {"type": "solid"}}},
            {"command": "model.boundary_condition.create", "params": {"name": "Clamp", "type": "fixed", "target": {"set": "root"}}},
            {"command": "model.load.create", "params": {"name": "Pull", "type": "force", "target": {"set": "tip"}, "vector": force}},
            {"command": "model.step.create", "params": {"name": "Load", "procedure": {"type": "static"}, "boundary_conditions": ["Clamp"], "loads": ["Pull"]}},
        ]);
        s.run_script(&script).unwrap();
        s.execute("solve.static", &json!({})).unwrap();
        let probe = s.execute("solve.probe", &json!({"step": "Load", "nodes": tip})).unwrap();
        let w: Vec<f64> = probe["nodes"].as_array().unwrap().iter().map(|n| n["displacement"][k].as_f64().unwrap()).collect();
        w.iter().sum::<f64>() / w.len() as f64
    }

    /// #439's acceptance: the plan proposes cuts that box the dome; with them the plate is hex
    /// cells and the dome's cell tet, joined by pyramids into one valid, conforming mesh whose
    /// static answer matches the all-tet mesh's within 2 %.
    #[test]
    fn a_dome_on_a_plate_is_boxed_into_a_tet_cell_and_the_rest_is_hex() {
        let (_dir, mut s) = plate_with_dome();
        let plan = s.execute("mesh.plan", &json!({"body": "Plate", "size": 4})).unwrap();
        let d = &plan["decomposition"];
        assert!(d["features"].as_array().is_some_and(|f| !f.is_empty()), "{plan}");
        let cuts = d["cuts"].as_array().unwrap();
        assert!(cuts.len() >= 4, "{d}");
        let commands: Vec<Value> = d["commands"].as_array().unwrap().clone();
        let (cut_commands, mesh) = commands.split_at(commands.len() - 1);
        s.run_script(&json!(cut_commands)).unwrap();
        // The cell inside the region is tet; the others take the default (auto: hex).
        let (lo, hi) = (&d["tet_region"]["min"], &d["tet_region"]["max"]);
        let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
        // The region is given along its frame's axes (here coordinate axes).
        let axis = |k: usize| (0..3).find(|&i| d["tet_region"]["frame"][k][i].as_f64().unwrap().abs() > 0.5).unwrap();
        let inside = |b: &Value| {
            (0..3).all(|k| {
                let i = axis(k);
                b["bbox"][0][i].as_f64().unwrap() >= lo[k].as_f64().unwrap() - 1e-6
                    && b["bbox"][1][i].as_f64().unwrap() <= hi[k].as_f64().unwrap() + 1e-6
            })
        };
        let tet: Vec<String> = bodies.as_array().unwrap().iter().filter(|b| inside(b)).map(|b| b["name"].as_str().unwrap().to_string()).collect();
        assert_eq!(tet.len(), 1, "{bodies}");
        let mut params = mesh[0]["params"].clone();
        params["methods"] = json!({tet[0].clone(): "tet"});
        let r = s.execute("mesh.cells", &params).unwrap();
        let kinds: BTreeSet<&str> = s.model.elements().values().map(|e| e.kind.name()).collect();
        assert_eq!(kinds.into_iter().collect::<Vec<_>>(), ["hex8", "pyramid5", "tet4"], "{r}");
        assert!(r["interfaces"].as_array().unwrap().iter().all(|i| i["joined"] != "tied"), "{r}");
        let dup = s.execute("mesh.check.duplicates", &json!({})).unwrap();
        assert_eq!(dup["counts"], json!({"coincident": 0, "cracks": 0, "duplicate_elements": 0, "free_edges": 0}), "{dup}");
        let q = s.execute("mesh.quality", &json!({})).unwrap();
        assert_eq!(q["invalid"]["count"], 0, "{q}");
        // Tension along the plate, which the dome stiffens locally. (In bending, linear hexahedra
        // and linear tetrahedra lock differently at three elements through the thickness, so
        // the two meshes would differ by their discretisation, not by the transition.)
        let hex = s.model.elements().values().filter(|e| e.kind.name() == "hex8").count();
        let mixed = pull(&mut s, [8000.0, 0.0, 0.0], 0);
        let (_dir2, mut t) = plate_with_dome();
        t.execute("mesh.tet.generate", &json!({"body": "Plate", "size": 4})).unwrap();
        let all_tet = pull(&mut t, [8000.0, 0.0, 0.0], 0);
        eprintln!("mixed {mixed:.6e} ({hex} hex of {} elements), all tet {all_tet:.6e}", s.model.elements().len());
        assert!(((mixed - all_tet) / all_tet).abs() < 0.02, "{mixed} vs {all_tet}");
    }
}
