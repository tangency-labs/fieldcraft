//! Independent LE7 shell model: ESRD NAFEMS guide (2018), p.23.
//! Mid-surface pressure and full-circumference consistent ring loads.
use fieldcraft_engine::Session;
use serde_json::{Value, json};
use std::f64::consts::PI;

fn recipe(per_m: usize) -> (Value, usize, usize) {
    let mut x = vec![[1.0, 0.0, 0.0]];
    let mut d = 0;
    for (a, b) in [(0.0_f64, 1.1136_f64), (1.1136, 1.4034), (1.4034, 1.5)] {
        let n = ((b - a) * per_m as f64).ceil() as usize;
        for i in 1..=n {
            x.push([1.0, a + (b - a) * i as f64 / n as f64, 0.0]);
        }
        if b == 1.4034 {
            d = x.len();
        }
    }
    let n = (PI / 2.0 * per_m as f64).ceil() as usize;
    for i in 1..=n {
        let a = PI / 2.0 * i as f64 / n as f64;
        x.push([if i == n { 0.0 } else { a.cos() }, 1.5 + a.sin(), 0.0]);
    }
    let mut c = vec![json!({"command":"model.units","params":{"system":"m-kg-s-N-Pa"}})];
    let mut f = vec![[0.0; 3]; x.len()];
    for (i, p) in x.iter().enumerate() {
        c.push(json!({"command":"mesh.node.add","params":{"label":i+1,"xyz":p}}));
    }
    for i in 0..x.len() - 1 {
        c.push(json!({"command":"mesh.element.add","params":{"label":i+1,"kind":"beam2","nodes":[i+1,i+2]}}));
        let dr = x[i + 1][0] - x[i][0];
        let dz = x[i + 1][1] - x[i][1];
        // 2*pi*p*integral N*r ds*n_out, n_out=(dz,-dr)/L.
        for (j, r) in [(i, (2.0 * x[i][0] + x[i + 1][0]) / 6.0), (i + 1, (x[i][0] + 2.0 * x[i + 1][0]) / 6.0)] {
            f[j][0] += 2.0 * PI * 1e6 * r * dz;
            f[j][1] -= 2.0 * PI * 1e6 * r * dr;
        }
    }
    c.extend([
        json!({"command":"model.set.create","params":{"name":"all","elements":(1..x.len()).collect::<Vec<_>>()}}),
        json!({"command":"model.set.create","params":{"name":"F","nodes":[1]}}),
        json!({"command":"model.material.create","params":{"name":"Steel","elastic":{"type":"isotropic","young":210e9,"poisson":0.3},"density":7800}}),
        json!({"command":"model.section.create","params":{"name":"shell","set":"all","material":"Steel","kind":{"type":"shell","thickness":0.025},"formulation":{"axisymmetric_shell":true}}}),
        json!({"command":"model.boundary_condition.create","params":{"name":"symmetry","type":"displacement","target":{"set":"F"},"values":[null,0,null,null,null,0]}})
    ]);
    for (i, f) in f.iter().enumerate() {
        c.push(json!({"command":"model.set.create","params":{"name":format!("N{i}"),"nodes":[i+1]}}));
        c.push(json!({"command":"model.load.create","params":{"name":format!("P{i}"),"type":"force","target":{"set":format!("N{i}")},"vector":f}}));
    }
    c.push(json!({"command":"model.step.create","params":{"name":"Pressure","procedure":{"type":"static"},"boundary_conditions":["symmetry"],"loads":(0..x.len()).map(|i|format!("P{i}")).collect::<Vec<_>>(),"output":["displacement","stress"]}}));
    c.push(json!({"command":"solve.static","params":{"step":"Pressure"}}));
    (json!(c), d, x.len())
}
#[test]
fn le7_shell_pole_commands_and_reopen() {
    for per_m in [100] {
        let (r, d, pole) = recipe(per_m);
        let mut s = Session::default();
        s.run_script(&r).unwrap();
        let p = json!({"step":"Pressure","nodes":[1,d,pole],"elements":[d-1,d]});
        let out = s.execute("solve.probe", &p).unwrap();
        // Positive shell normal is inward: bottom is the physical outer face.
        let left = out["elements"][0]["bottom"][1][1].as_f64().unwrap();
        let right = out["elements"][1]["bottom"][0][1].as_f64().unwrap();
        for value in [left, right] {
            assert!((value / 25.86e6 - 1.0).abs() < 0.01, "{out}");
        }
        assert_eq!(out["nodes"][2]["displacement"][0], 0.0);
        assert_eq!(out["nodes"][2]["rotation"][2], 0.0);
        assert!((out["nodes"][0]["reaction"][1].as_f64().unwrap() / (-PI * 1e6) - 1.0).abs() < 1e-8);
        println!(
            "LE7 per_m={per_m}, D={d}, pole={pole}, outer_left={left}, outer_right={right}, mean_relative_error={}",
            (0.5 * (left + right) / 25.86e6 - 1.0)
        );
        let dir = std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS").map(std::path::PathBuf::from).unwrap_or_else(std::env::temp_dir);
        let path = dir.join(format!("le7-shell-{per_m}-{}.fieldcraft", std::process::id()));
        s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
        let mut open = Session::default();
        open.execute("file.open", &json!({"path":path})).unwrap();
        assert_eq!(out, open.execute("solve.probe", &p).unwrap());
        s.execute(
            "model.boundary_condition.create",
            &json!({"name":"bad_pole","type":"displacement","target":{"set":format!("N{}",pole-1)},"values":[0.001,null,null,null,null,null]}),
        )
        .unwrap();
        let mut step = r.as_array().unwrap().iter().find(|c| c["command"] == "model.step.create").unwrap()["params"].clone();
        step["boundary_conditions"] = json!(["symmetry", "bad_pole"]);
        s.execute("model.step.edit", &step).unwrap();
        assert!(s.execute("solve.static", &json!({"step":"Pressure"})).unwrap_err().to_string().contains("pole"));
        if std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS").is_some() {
            std::fs::write(dir.join(format!("le7-shell-{per_m}-recipe.json")), serde_json::to_vec_pretty(&r).unwrap()).unwrap();
            std::fs::write(dir.join(format!("le7-shell-{per_m}-probe.json")), serde_json::to_vec_pretty(&out).unwrap()).unwrap();
        } else {
            std::fs::remove_file(path).unwrap();
        }
    }
}
