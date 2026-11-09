//! Classical perfect-cylinder bifurcation, independently generated quarter meshes.
//! Axial: NASA SP-8007 classical elastic result (knockdown factor = 1).
//! Dead radial traction: inextensional n=2 ring, p=4D/R^3; NOT follower pressure
//! (whose ring limit is 3D/R^3). The shell solver's load vector is fixed.
use fieldcraft_engine::Session;
use fieldcraft_fem::basis::Basis;
use serde_json::{Value, json};
use std::{collections::BTreeMap, f64::consts::PI};

fn model(tri: bool, pressure: bool, follower: bool, divisions: usize) -> (Value, f64) {
    let h = 0.02_f64;
    let young = 1e7;
    let nu = 0.0_f64;
    let length = if pressure { 0.2 } else { PI * h.sqrt() / 12.0_f64.powf(0.25) };
    let nz = if pressure { 2 } else { divisions };
    let nt = divisions;
    let mut nodes = BTreeMap::new();
    let mut xyz = vec![];
    let mut commands = vec![json!({"command":"model.units","params":{"system":"m-kg-s-N-Pa"}})];
    for j in 0..=2 * nz {
        for i in 0..=2 * nt {
            if !tri && i % 2 == 1 && j % 2 == 1 {
                continue;
            }
            let theta = PI / 2.0 * i as f64 / (2 * nt) as f64;
            let x = [if i == 2 * nt { 0.0 } else { theta.cos() }, theta.sin(), length * j as f64 / (2 * nz) as f64];
            let label = xyz.len() + 1;
            nodes.insert((i, j), label);
            xyz.push(x);
            commands.push(json!({"command":"mesh.node.add","params":{"label":label,"xyz":x}}));
        }
    }
    let mut elements = vec![];
    for j in 0..nz {
        for i in 0..nt {
            let a = 2 * i;
            let b = 2 * j;
            let n = |i, j| nodes[&(i, j)];
            if tri {
                elements.push(vec![n(a, b), n(a + 2, b), n(a + 2, b + 2), n(a + 1, b), n(a + 2, b + 1), n(a + 1, b + 1)]);
                elements.push(vec![n(a, b), n(a + 2, b + 2), n(a, b + 2), n(a + 1, b + 1), n(a + 1, b + 2), n(a, b + 1)]);
            } else {
                elements.push(vec![n(a, b), n(a + 2, b), n(a + 2, b + 2), n(a, b + 2), n(a + 1, b), n(a + 2, b + 1), n(a + 1, b + 2), n(a, b + 1)]);
            }
        }
    }
    for (i, n) in elements.iter().enumerate() {
        commands.push(json!({"command":"mesh.element.add","params":{"label":i+1,"kind":if tri {"tri6"} else {"quad8"},"nodes":n}}));
    }
    commands.extend([
        json!({"command":"model.set.create","params":{"name":"all","elements":(1..=elements.len()).collect::<Vec<_>>(),"nodes":(1..=xyz.len()).collect::<Vec<_>>()}}),
        json!({"command":"model.material.create","params":{"name":"E","elastic":{"type":"isotropic","young":young,"poisson":nu}}}),
        json!({"command":"model.section.create","params":{"name":"shell","set":"all","material":"E","kind":{"type":"shell","thickness":h}}})
    ]);
    let mut bcs = vec![];
    let mut bc = |name: &str, labels: Vec<usize>, values: Value| {
        commands.push(json!({"command":"model.set.create","params":{"name":name,"nodes":labels}}));
        commands.push(
            json!({"command":"model.boundary_condition.create","params":{"name":name,"type":"displacement","target":{"set":name},"values":values}}),
        );
        bcs.push(name.to_owned());
    };
    bc("theta0", nodes.iter().filter(|((i, _), _)| *i == 0).map(|(_, n)| *n).collect(), json!([null, 0, null, 0, null, 0]));
    bc("theta90", nodes.iter().filter(|((i, _), _)| *i == 2 * nt).map(|(_, n)| *n).collect(), json!([0, null, null, null, 0, 0]));
    if pressure {
        bc("plane_strain", nodes.values().copied().collect(), json!([null, null, 0, null, null, null]));
    } else {
        bc("bottom", nodes.iter().filter(|((_, j), _)| *j == 0).map(|(_, n)| *n).collect(), json!([0, 0, 0, null, null, null]));
        bc("top", nodes.iter().filter(|((_, j), _)| *j == 2 * nz).map(|(_, n)| *n).collect(), json!([0, 0, null, null, null, null]));
    }
    let mut loads = vec![[0.0; 3]; xyz.len()];
    if pressure {
        // Integrate fixed normal traction on the original isoparametric surface.
        let basis = if tri { Basis::Tri6 } else { Basis::Quad8 };
        let rule = if tri { fieldcraft_fem::quadrature::Rule::TriConical(3) } else { fieldcraft_fem::quadrature::Rule::Quad(3) };
        for conn in &elements {
            for q in rule.points() {
                let sh = basis.eval(q.xi);
                let mut g = [[0.0; 3]; 2];
                for (a, n) in conn.iter().enumerate() {
                    for (b, tangent) in g.iter_mut().enumerate() {
                        for (k, component) in tangent.iter_mut().enumerate() {
                            *component += sh.dn[a][b] * xyz[n - 1][k];
                        }
                    }
                }
                let normal = fieldcraft_fem::frame::cross(g[0], g[1]);
                for (a, n) in conn.iter().enumerate() {
                    for k in 0..3 {
                        loads[n - 1][k] -= q.weight * sh.n[a] * normal[k];
                    }
                }
            }
        }
    } else {
        for i in 0..=2 * nt {
            let w = if i == 0 || i == 2 * nt {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            };
            loads[nodes[&(i, 2 * nz)] - 1][2] = -w * PI / (12.0 * nt as f64);
        }
    }
    let mut names = vec![];
    for (i, f) in loads.iter().enumerate() {
        if !follower && f.iter().any(|v| *v != 0.0) {
            let name = format!("F{i}");
            names.push(name.clone());
            commands.push(json!({"command":"model.set.create","params":{"name":name,"nodes":[i+1]}}));
            commands.push(json!({"command":"model.load.create","params":{"name":name,"type":"force","target":{"set":name},"vector":f}}));
        }
    }
    if follower {
        commands.push(json!({"command":"model.set.create","params":{"name":"pressure_faces","faces":(1..=elements.len()).map(|i|json!({"element":i,"face":1})).collect::<Vec<_>>()}}));
        commands.push(json!({"command":"model.load.create","params":{"name":"pressure","type":"pressure","set":"pressure_faces","magnitude":1}}));
        names = vec!["pressure".to_owned()];
    }
    commands.push(json!({"command":"model.step.create","params":{"name":"Buckle","procedure":{"type":"buckle","modes":2,"follower_pressure":follower},"boundary_conditions":bcs,"loads":names}}));
    let reference = if pressure {
        (if follower { 3.0 } else { 4.0 }) * young * h.powi(3) / (12.0 * (1.0 - nu * nu))
    } else {
        young * h * h / (3.0 * (1.0 - nu * nu)).sqrt()
    };
    (json!(commands), reference)
}
fn run_case(tri: bool, pressure: bool, follower: bool) {
    let (recipe, reference) = model(tri, pressure, follower, 4);
    let mut s = Session::default();
    s.run_script(&recipe).unwrap();
    let result = s.execute("solve.buckle", &json!({"step":"Buckle"})).unwrap();
    let factor = result["modes"][0]["load_factor"].as_f64().unwrap();
    println!(
        "curved cylinder tri={tri} pressure={pressure} follower={follower} reference={reference} actual={factor} relative_error={} result={result}",
        factor / reference - 1.0
    );
    let dir = std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS").map(std::path::PathBuf::from).unwrap_or_else(std::env::temp_dir);
    let stem = format!("cylinder-tri{tri}-pressure{pressure}-follower{follower}");
    let archive = dir.join(format!("{stem}-{}.fieldcraft", std::process::id()));
    let p = json!({"step":"Buckle","nodes":[1,2,3],"mode":1});
    let before = s.execute("solve.probe", &p).unwrap();
    s.execute("file.save_as", &json!({"path":archive,"results":true})).unwrap();
    let mut open = Session::default();
    open.execute("file.open", &json!({"path":archive})).unwrap();
    assert_eq!(before, open.execute("solve.probe", &p).unwrap());
    if std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS").is_some() {
        std::fs::write(dir.join(format!("{stem}-recipe.json")), serde_json::to_vec_pretty(&recipe).unwrap()).unwrap();
        std::fs::write(dir.join(format!("{stem}-result.json")), serde_json::to_vec_pretty(&result).unwrap()).unwrap();
    } else {
        std::fs::remove_file(archive).unwrap();
    }
    // Fixed before runs: coarse quadratic mesh, thin-shell asymptotic reference.
    assert!((factor / reference - 1.0).abs() < 0.05, "{result}");
}
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn cylinders_through_buckling_commands_and_reopen() {
    for tri in [false, true] {
        for pressure in [false, true] {
            run_case(tri, pressure, false);
        }
    }
}
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn follower_cylinder_pressure_and_reopen() {
    for tri in [false, true] {
        run_case(tri, true, true);
    }
}
