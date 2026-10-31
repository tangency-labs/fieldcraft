//! The modelling skill's verification recipes (`skills/fieldcraft-modelling/references`) run
//! headless with today's commands and give the numbers `checks.md` states, so a command or
//! schema change cannot leave agents with stale recipes.

use serde_json::{Value, json};

use crate::Session;

const AXIAL_ROD: &str = include_str!("../../../skills/fieldcraft-modelling/references/axial-rod.json");
const STEADY_CONDUCTOR: &str = include_str!("../../../skills/fieldcraft-modelling/references/steady-conductor.json");

fn close(got: &Value, want: f64) -> bool {
    got.as_f64().is_some_and(|g| (g - want).abs() <= 1e-9 * want.abs().max(1.0))
}

/// Run a recipe in a fresh session and return its replies.
fn run(recipe: &Value) -> (Session, Vec<Value>) {
    let mut s = Session::default();
    let out = s.run_script(recipe).unwrap();
    (s, out)
}

#[test]
fn axial_rod_gives_the_closed_form_extension_stress_and_reaction() {
    let (mut s, out) = run(&serde_json::from_str(AXIAL_ROD).unwrap());
    let probe = out.last().unwrap();
    let (root, tip) = (&probe["nodes"][0], &probe["nodes"][1]);
    // u = F L / (E A) = 1000 · 0.3 / (70e9 · 1e-4), sigma = F / A = 10 MPa.
    assert!(close(&tip["displacement"][0], 4.285714285714286e-5), "{tip}");
    assert!(close(&root["reaction"][0], -1000.0), "{root}");
    assert!(close(&tip["stress"][0], 10e6), "{tip}");

    // Only the display changes: the same physical extension and stress, in mm and MPa.
    s.execute("model.units", &json!({"system": "mm-t-s-N-MPa"})).unwrap();
    let mm = s.execute("solve.probe", &json!({"step": "Axial", "nodes": [2]})).unwrap();
    assert!(close(&mm["nodes"][0]["displacement"][0], 0.04285714285714286), "{mm}");
    assert!(close(&mm["nodes"][0]["stress"][0], 10.0), "{mm}");
}

/// The steady conductor's midpoint temperature, axial flux, and the signed heat flow through
/// its cold and hot ends.
fn conductor(divisions: [u32; 3]) -> (f64, f64, f64, f64) {
    let mut recipe: Value = serde_json::from_str(STEADY_CONDUCTOR).unwrap();
    let block = recipe.as_array_mut().unwrap().iter_mut().find(|step| step["command"] == "mesh.block").unwrap();
    block["params"]["divisions"] = json!(divisions);
    let (_, out) = run(&recipe);
    let n = out.len();
    let mid = &out[n - 3]["nodes"][0];
    let total = |probe: &Value| probe["nodes"].as_array().unwrap().iter().map(|node| node["reaction"].as_f64().unwrap()).sum::<f64>();
    (mid["temperature"].as_f64().unwrap(), mid["heat_flux"][2].as_f64().unwrap(), total(&out[n - 2]), total(&out[n - 1]))
}

#[test]
fn steady_conductor_gives_the_linear_profile_and_opposite_end_rates() {
    for divisions in [[2, 3, 6], [4, 6, 12]] {
        let (mid, q_z, cold, hot) = conductor(divisions);
        // T linear from 293.15 K to 353.15 K over 0.12 m, k = 15 W/(m K), area 6e-4 m².
        assert!(close(&json!(mid), 323.15), "{divisions:?}: {mid}");
        assert!((q_z + 7500.0).abs() < 1e-6, "{divisions:?}: {q_z}");
        assert!((cold + 4.5).abs() < 1e-9 && (hot - 4.5).abs() < 1e-9, "{divisions:?}: {cold} {hot}");
    }
}
