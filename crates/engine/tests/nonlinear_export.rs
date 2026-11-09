//! Native writer and postprocessing of retained real histories; never reruns time integration.
use fieldcraft_engine::Session;
use serde_json::{Value, json};
use std::{fs, path::Path};
fn model() -> Value {
    json!([
    {"command":"mesh.node.add","params":{"label":1,"xyz":[0,0,0]}},
    {"command":"mesh.node.add","params":{"label":2,"xyz":[0,-1,0]}},
    {"command":"mesh.element.add","params":{"label":1,"kind":"rod2","nodes":[1,2]}},
    {"command":"model.set.create","params":{"name":"rod","elements":[1]}},
    {"command":"model.set.create","params":{"name":"fixed","nodes":[1]}},
    {"command":"model.material.create","params":{"name":"elastic","elastic":{"type":"isotropic","young":1000,"poisson":0.25},"density":3}},
    {"command":"model.section.create","params":{"name":"rod","set":"rod","material":"elastic","kind":{"type":"rod","area":1}}},
    {"command":"model.boundary_condition.create","params":{"name":"fixed","type":"displacement","target":{"set":"fixed"},"values":[0,0,0,null,null,null]}},
    {"command":"model.step.create","params":{"name":"Motion","procedure":{"type":"transient","time_step":0.01,"increments":2,"initial_displacement":{"2":[0.8660254037844386,0.5,0]}},"boundary_conditions":["fixed"],"output":["displacement"]}}
    ])
}
#[test]
fn native_nonlinear_export_preserves_initial_state_and_refuses_lossy_options() {
    let dir = std::env::temp_dir().join(format!("fieldcraft-nonlinear-export-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let mut s = Session::default();
    s.execute("model.units", &json!({"system":"m-kg-s-N-Pa"})).unwrap();
    s.run_script(&model()).unwrap();
    let path = dir.join("new.inp");
    let reply = s.execute("deck.write.nonlinear", &json!({"path":path,"mass":"consistent"})).unwrap();
    let text = fs::read_to_string(&path).unwrap();
    assert!(text.contains("NLGEOM=YES") && text.contains("*DYNAMIC, DIRECT, NOHAF"));
    assert!(text.contains("NAME=FC_INITIALIZE") && text.contains("*BOUNDARY, OP=NEW"));
    assert!(text.contains("*USER ELEMENT") && !text.contains("TYPE=T3D2"));
    assert!(text.contains("*NODE OUTPUT, NSET=FC_NODES") && text.contains("U, V, A, RF"));
    assert_eq!(reply["output_precision"], "full");
    assert!(Path::new(reply["companion"].as_str().unwrap()).exists());
    let rejected = dir.join("rejected.inp");
    assert!(s.execute("deck.write.nonlinear", &json!({"path":rejected,"representation":"t3d2"})).unwrap_err().to_string().contains("StVK"));
    assert!(!rejected.exists() && !rejected.with_extension("c").exists());
    let mut bad = model();
    bad.as_array_mut().unwrap().last_mut().unwrap()["params"]["procedure"]["initial_velocity"] = json!({"2":[1,0,0]});
    let mut other = Session::default();
    other.execute("model.units", &json!({"system":"m-kg-s-N-Pa"})).unwrap();
    other.run_script(&bad).unwrap();
    assert!(other.execute("deck.write.nonlinear", &json!({"path":rejected})).unwrap_err().to_string().contains("initial_velocity"));
    fs::remove_dir_all(dir).unwrap();
}
#[test]
fn retained_histories_are_exported_and_recovered_without_solving() {
    let Ok(source) = std::env::var("FIELDCRAFT_NONLINEAR_RETAINED") else {
        return;
    };
    let output = std::env::var("FIELDCRAFT_NONLINEAR_EXPORT").unwrap();
    let source = Path::new(&source);
    let output = Path::new(&output);
    fs::create_dir_all(output).unwrap();
    let mut manifest = Vec::new();
    for shape in ["pendulum", "string"] {
        for mass in ["consistent", "lumped"] {
            for method in ["newmark", "hht"] {
                let name = format!("{shape}-{mass}-{method}-c0-dt0.0025");
                let original = source.join(&name);
                let out = output.join(&name);
                fs::create_dir_all(&out).unwrap();
                let input = fs::read(original.join("input.json")).unwrap();
                let actual = fs::read(original.join("actual.json")).unwrap();
                let commands: Value = serde_json::from_slice(&input).unwrap();
                let states: Value = serde_json::from_slice(&actual).unwrap();
                let mut s = Session::default();
                s.run_script(&commands).unwrap();
                let path = out.join("motion.inp");
                let reply = s.execute("deck.write.nonlinear", &json!({"path":path,"step":"Motion","mass":mass,"representation":"uel"})).unwrap();
                let resolved = s.model.resolve().unwrap();
                let labels: Vec<u32> = serde_json::from_value(states["nodes"].clone()).unwrap();
                let mut recovered = Vec::new();
                for frame in states["frames"].as_array().unwrap() {
                    let displacements: Vec<[f64; 3]> = serde_json::from_value(frame["displacement"].clone()).unwrap();
                    let mut forces = vec![[0.; 3]; labels.len()];
                    let mut elements = Vec::new();
                    let mut energy = 0.;
                    for (id, e) in s.model.elements() {
                        let x: Vec<_> = e.nodes.iter().map(|n| s.model.nodes().get(n).unwrap().xyz).collect();
                        let mut u = Vec::new();
                        let mut indices = Vec::new();
                        for n in &e.nodes {
                            let i = labels.iter().position(|v| v == n).unwrap();
                            indices.push(i);
                            u.extend(displacements[i]);
                        }
                        let section = s.model.sections().values().find(|v| resolved.sets[&v.set].elements.contains(id)).unwrap();
                        let fieldcraft_model::SectionKind::Rod { area } = section.kind else { panic!("rod") };
                        let Some(fieldcraft_model::Elastic::Isotropic { young, .. }) = s.model.material(&section.material).unwrap().elastic else {
                            panic!("isotropic")
                        };
                        let r = fieldcraft_solver::transient::nonlinear::recover_rod(&x, &u, young, area).unwrap();
                        for (local, index) in indices.iter().enumerate() {
                            for (d, value) in forces[*index].iter_mut().enumerate() {
                                *value += r.nodal_internal_force[3 * local + d];
                            }
                        }
                        energy += r.strain_energy;
                        elements.push(json!({"label":id,"nodes":e.nodes,"response":r}));
                    }
                    assert!((energy - frame["strain_energy"].as_f64().unwrap()).abs() < 1e-12);
                    recovered.push(json!({"time":frame["time"],"nodes":labels,"internal_force":forces,"elements":elements}));
                }
                fs::write(out.join("recovered-internal.json"), serde_json::to_vec_pretty(&recovered).unwrap()).unwrap();
                fs::write(out.join("native-reply.json"), serde_json::to_vec_pretty(&reply).unwrap()).unwrap();
                fs::write(
                    out.join("export-command.json"),
                    serde_json::to_vec_pretty(
                        &json!({"command":"deck.write.nonlinear","params":{"path":"motion.inp","step":"Motion","mass":mass,"representation":"uel"}}),
                    )
                    .unwrap(),
                )
                .unwrap();
                manifest.push(json!({"case":name,"input_sha256":fieldcraft_io::sha256_hex(&input),"actual_sha256":fieldcraft_io::sha256_hex(&actual),"time_integrations":0,"frames":recovered.len(),"force_provenance":"recovered with production rod kernel from retained actual displacement, not previously stored forces"}));
            }
        }
    }
    fs::write(output.join("manifest.json"), serde_json::to_vec_pretty(&manifest).unwrap()).unwrap();
    assert_eq!(manifest.len(), 8);
}
