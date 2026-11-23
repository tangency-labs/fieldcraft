//! The checked-in workshop scripts exercise the same commands as the CLI.
use fieldcraft_engine::Session;
use serde_json::Value;

fn check_recipe(script: &str, dimension: u8, boundary_names: &[&str]) {
    let mut commands: Value = serde_json::from_str(script).unwrap();
    let export = commands.as_array_mut().unwrap().last_mut().unwrap();
    assert_eq!(export["command"], "mesh.export.vol");
    let dir = std::env::temp_dir().join(format!("fieldcraft-vol-recipe-{}-{dimension}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(export["params"]["path"].as_str().unwrap());
    export["params"]["path"] = serde_json::to_value(&path).unwrap();
    let mut session = Session::default();
    let replies = session.run_script(&commands).unwrap();
    let summary = replies.last().unwrap();
    assert_eq!(summary["dimension"], dimension);
    assert_eq!(summary["cells"].as_u64().unwrap() as usize, session.model.elements().len());
    assert!(summary["domains"].as_array().unwrap().iter().all(|domain| domain["source"] == "element_set"));
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.starts_with(&format!("mesh3d\ndimension\n{dimension}\n")));
    let mut records = text.lines().skip_while(|&l| l != "bcnames").skip(1);
    let count: usize = records.next().unwrap().parse().unwrap();
    let names: Vec<&str> = records.take(count).map(|line| line.split_once(' ').unwrap().1).collect();
    for boundary in boundary_names {
        assert!(names.contains(boundary), "missing boundary `{boundary}`");
    }
    assert!(text.ends_with("endmesh\n"));
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn planar_mems_workshop_exports_its_named_edges() {
    check_recipe(
        include_str!("../../../docs/recipes/netgen-mems-pull-in.json"),
        2,
        &["ground", "outer", "beam_electrode", "tip", "beam_top", "clamp"],
    );
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn team13_workshop_exports_its_named_interfaces() {
    check_recipe(include_str!("../../../docs/recipes/netgen-team13.json"), 3, &["outer", "gamma_core", "gamma_coil", "winding_seams"]);
}
