//! The filleted bracket meshed all-hex by partitioning (`tests/corpus/filleted-bracket-blocks`):
//! the upright's end clamped, the foot's end pushed down by 1 kN (steel). The conforming hex8 and
//! hex20 meshes converge: the tip deflection, the von Mises stress at the middle of the fillet and
//! the peak (at the upright's hole) agree between hex20 at 3 mm and 2 mm, and hex8 approaches them.
//! For reference, a 48 016-element tet10 mesh of the unpartitioned part (curvature-refined, from
//! the tet quality pipeline) gives a tip deflection of 0.5093 mm, 133.2 MPa at the fillet and a
//! peak of 291.8 MPa; hex20 at 2 mm gives 0.5095 mm, 132.4 MPa and 297.0 MPa.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use fieldcraft_engine::Session;
use fieldcraft_testkit::corpus;
use serde_json::{Value, json};

fn open(part: &str) -> (Session, std::path::PathBuf) {
    let parts = corpus::parts().unwrap();
    let p = parts.iter().find(|p| p.name == part).unwrap();
    let dir = std::env::temp_dir().join(format!("fieldcraft-bracket-{}-{part}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("part.solidcraft");
    std::fs::write(&path, corpus::design_bytes(p).unwrap()).unwrap();
    let mut s = Session::default();
    s.execute("geometry.open", &json!({"path": path})).unwrap();
    (s, dir)
}

/// Every (body, face) whose name is `face` (a face split into cells keeps its name in each).
fn faces_named(s: &mut Session, face: &str) -> Vec<Value> {
    let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
    let mut out = Vec::new();
    for b in bodies.as_array().unwrap() {
        for f in b["faces"].as_array().unwrap() {
            if f["name"] == face {
                out.push(json!({"body": b["name"], "face": face}));
            }
        }
    }
    out
}

/// Clamp the upright's top (y = 50), push the foot's end (x = 60) down by 1 kN; tip deflection,
/// von Mises at the middle of the fillet and the largest von Mises.
fn analyse(s: &mut Session) -> (f64, f64, f64) {
    let top = faces_named(s, "F2:side:l5");
    let tip = faces_named(s, "F2:side:l2");
    assert!(!top.is_empty() && !tip.is_empty());
    for (name, geometry) in [
        (
            "all",
            json!(
                s.execute("geometry.bodies", &json!({})).unwrap().as_array().unwrap().iter().map(|b| json!({"body": b["name"]})).collect::<Vec<_>>()
            ),
        ),
        ("top", json!(top)),
        ("tip", json!(tip)),
    ] {
        s.execute("model.set.create", &json!({"name": name, "geometry": geometry})).unwrap();
    }
    s.execute("model.material.create", &json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 210000.0, "poisson": 0.3}})).unwrap();
    s.execute("model.section.create", &json!({"name": "solid", "set": "all", "material": "Steel", "kind": {"type": "solid"}})).unwrap();
    s.execute("model.boundary_condition.create", &json!({"name": "clamp", "type": "fixed", "target": {"set": "top"}})).unwrap();
    s.execute("model.load.create", &json!({"name": "push", "type": "force", "target": {"set": "tip"}, "vector": [0.0, -1000.0, 0.0]})).unwrap();
    s.execute("model.step.create", &json!({"name": "Load", "procedure": {"type": "static"}, "boundary_conditions": ["clamp"], "loads": ["push"]}))
        .unwrap();
    let r = s.execute("solve.static", &json!({})).unwrap();
    let probe = s.execute("solve.probe", &json!({"step": "Load", "at": [[60.0, 4.0, 20.0], [9.757359312880714, 9.757359312880714, 20.0]]})).unwrap();
    let nodes = probe["nodes"].as_array().unwrap();
    let tip = nodes[0]["displacement"][1].as_f64().unwrap();
    let fillet = nodes[1]["von_mises"].as_f64().unwrap();
    (tip, fillet, r["max_von_mises"]["value"].as_f64().unwrap())
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn partitioned_bracket_converges_with_hex8_and_hex20() {
    let mut results = Vec::new();
    for (order, size) in [(1, 2.0), (2, 3.0), (2, 2.0)] {
        let (mut s, dir) = open("filleted-bracket-blocks");
        let m = s.execute("mesh.hex.blocks", &json!({"size": size, "order": order})).unwrap();
        assert_eq!(m["cells"], 30, "{m}");
        assert!(m["quality"]["min_scaled_jacobian"].as_f64().unwrap() > 0.6, "{m}");
        let r = analyse(&mut s);
        eprintln!("order {order}, size {size}: {} elements; tip {:.5} mm, fillet {:.2} MPa, peak {:.2} MPa", m["elements"], r.0, r.1, r.2);
        results.push(r);
        std::fs::remove_dir_all(dir).ok();
    }
    let close = |a: f64, b: f64, rel: f64| (a - b).abs() <= rel * b.abs();
    let [hex8, coarse, fine] = [results[0], results[1], results[2]];
    // The tip moves down; hex20 has converged on it, hex8 is within a few percent (stiffer).
    assert!(fine.0 < 0.0);
    assert!(close(coarse.0, fine.0, 0.005), "{coarse:?} {fine:?}");
    assert!(close(hex8.0, fine.0, 0.04) && hex8.0.abs() < fine.0.abs(), "{hex8:?} {fine:?}");
    // Stresses at the fillet and the peak.
    assert!(close(coarse.1, fine.1, 0.06), "{coarse:?} {fine:?}");
    assert!(close(coarse.2, fine.2, 0.02), "{coarse:?} {fine:?}");
    assert!(close(hex8.1, fine.1, 0.15) && close(hex8.2, fine.2, 0.10), "{hex8:?} {fine:?}");
}
