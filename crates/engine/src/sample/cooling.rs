//! The electronics-cooling studies (docs/examples.md, discussion #584): a heat sink in forced
//! air and a chip on a PCB in a sealed box (natural convection). Each is thermal lattice Boltzmann flow past solids held at a temperature: there is no
//! conduction inside solids yet (#157), so a high-conductivity sink or radiator is an isothermal
//! surface and a component a held-temperature face.
//!
//! Every case is a function of its [`Run`] (cell size, time step, how long): the samples open
//! with the published settings ([`HEAT_SINK`], [`ENCLOSURE`]), the trunk tests solve a
//! coarse, short version and the heavy tests hold the published numbers.
//! `tools/cooling/make.py` opens the samples, solves them and renders the gallery's pictures.

use fieldcraft_geometry::Shape;
use serde_json::{Value, json};

/// The geometry: SolidCraft recipes, the sink's and the enclosure's in the corpus.
pub(super) const HEAT_SINK_DESIGN: &str = include_str!("../../../../tests/corpus/cooling-heat-sink/recipe.json");
pub(super) const ENCLOSURE_DESIGN: &str = include_str!("../../../../tests/corpus/cooling-enclosure/recipe.json");

/// How a case is solved.
#[derive(Clone, Copy, Debug)]
struct Run {
    /// Cell size, mm.
    spacing: f64,
    /// s; `None` lets the solver choose.
    time_step: Option<f64>,
    /// s: the end, or with `steady` the most it runs.
    end_time: f64,
    frames: u32,
    /// Stop once the flow, the temperature and the heat the hot body gives the air settle (#595).
    steady: bool,
}

/// Boussinesq expansion coefficient of air near room temperature, 1/K.
const AIR_EXPANSION: f64 = 3.4e-3;
/// Gravity along −z, mm/s².
const GRAVITY: [f64; 3] = [0.0, 0.0, -9810.0];

/// Heat sink: air speed (mm/s), sink and inlet temperatures (°C).
const SINK_SPEED: f64 = 100.0;
const SINK_TEMPERATURE: f64 = 50.0;
const SINK_INLET: f64 = 25.0;
/// The duct around the sink, mm.
const SINK_DUCT: ([f64; 3], [f64; 3]) = ([-30.0, -2.0, 0.0], [90.0, 42.0, 30.0]);
/// Fixed length: the wake behind the sink keeps shedding (the velocity never settles to the steady
/// rule's tolerance), while the heat holds within about 1 % from 2 s on.
const HEAT_SINK: Run = Run { spacing: 1.0, time_step: Some(0.00025), end_time: 3.0, frames: 30, steady: false };

/// Enclosure: chip and wall temperatures (°C).
const CHIP_TEMPERATURE: f64 = 40.0;
const ENCLOSURE_WALLS: f64 = 25.0;
/// To steady state (about 18 s of flow): the chip's heat and the walls' then agree.
const ENCLOSURE: Run = Run { spacing: 1.0, time_step: Some(0.0002), end_time: 40.0, frames: 40, steady: true };

/// The material every case's air is made of.
fn air() -> Value {
    let [density, viscosity] = super::AIR;
    json!({"name": "Air", "density": density, "fluid": {"type": "newtonian", "viscosity": viscosity},
        "thermal": {"conductivity": 0.026, "specific_heat": 1.005e9}})
}

/// Every face of `body`, as `{"body", "face"}`.
fn faces(bodies: &[Shape], body: &str) -> Vec<Value> {
    bodies.iter().filter(|b| b.name == body).flat_map(|b| b.faces().iter().map(move |f| json!({"body": b.name, "face": f.name}))).collect()
}

/// Every face of `body` held at `value` (an isothermal high-conductivity solid, or a component).
fn held(bodies: &[Shape], body: &str, value: f64) -> Vec<Value> {
    faces(bodies, body).into_iter().map(|face| json!({"face": face, "condition": {"type": "temperature", "value": value}})).collect()
}

/// The fluid step: `run`'s time and frames, and with `steady` a stop once the heat on `hot` settles.
fn step(name: &str, domain: &str, run: &Run, hot: Vec<Value>) -> Value {
    let mut procedure = json!({"type": "fluid", "domain": domain, "end_time": run.end_time, "frames": run.frames});
    if let Some(dt) = run.time_step {
        procedure["time_step"] = json!(dt);
    }
    if run.steady {
        procedure["steady"] = json!({"tolerance": 1e-3, "monitors": [{"type": "heat", "faces": hot}]});
    }
    json!({"command": "model.step.create", "params": {"name": name, "procedure": procedure}})
}

/// A case's commands: the units, the air, the fluid domain `domain` (its `params` but the name)
/// and the step `name` over it, with the heat of `hot` as the steady rule's monitor.
fn script(domain: &str, mut params: Value, name: &str, run: &Run, hot: Vec<Value>) -> Value {
    params["name"] = json!(domain);
    json!([
        {"command": "model.units", "params": {"system": "mm-t-s-N-MPa"}},
        {"command": "model.material.create", "params": air()},
        {"command": "model.fluid_domain.create", "params": params},
        step(name, domain, run, hot),
    ])
}

pub(super) fn heat_sink(bodies: &[Shape]) -> Value {
    heat_sink_with(bodies, &HEAT_SINK)
}

/// Mixed convection: air at [`SINK_INLET`] pushed at [`SINK_SPEED`] along a six-fin sink held at
/// [`SINK_TEMPERATURE`] in a closed duct, with buoyancy (at 0.1 m/s and 25 K, Ri ≈ 3).
fn heat_sink_with(bodies: &[Shape], run: &Run) -> Value {
    let (min, max) = SINK_DUCT;
    let v = [SINK_SPEED, 0.0, 0.0];
    script(
        "Duct",
        json!({
            "material": "Air", "min": min, "max": max, "spacing": run.spacing,
            "faces": {"x_min": {"type": "velocity", "velocity": v}, "x_max": {"type": "pressure", "pressure": 0},
                "y_min": {"type": "wall"}, "y_max": {"type": "wall"}, "z_min": {"type": "wall"}, "z_max": {"type": "wall"}},
            "obstacles": ["HeatSink"], "initial_velocity": v,
            "thermal": {"initial_temperature": SINK_INLET, "gravity": GRAVITY, "expansion": AIR_EXPANSION,
                "faces": {"x_min": {"type": "temperature", "value": SINK_INLET}},
                "surfaces": held(bodies, "HeatSink", SINK_TEMPERATURE)}
        }),
        "Cool",
        run,
        faces(bodies, "HeatSink"),
    )
}

pub(super) fn enclosure(bodies: &[Shape]) -> Value {
    enclosure_with(bodies, &ENCLOSURE)
}

/// Natural convection: a chip held at [`CHIP_TEMPERATURE`] on a vertical PCB in a sealed
/// 60 × 40 × 60 mm box whose walls are held at [`ENCLOSURE_WALLS`]; buoyancy drives the plume.
/// The buoyant start-up needs a small fixed time step. Convection only, no radiation.
fn enclosure_with(bodies: &[Shape], run: &Run) -> Value {
    let sides = ["x_min", "x_max", "y_min", "y_max", "z_min", "z_max"];
    let walls: serde_json::Map<String, Value> = sides.iter().map(|f| (f.to_string(), json!({"type": "wall"}))).collect();
    let held_walls: serde_json::Map<String, Value> =
        sides.iter().map(|f| (f.to_string(), json!({"type": "temperature", "value": ENCLOSURE_WALLS}))).collect();
    script(
        "Box",
        json!({
            "material": "Air", "min": [0, 0, 0], "max": [60, 40, 60], "spacing": run.spacing,
            "faces": walls, "obstacles": ["PCB", "Chip"], "initial_velocity": [0, 0, 0],
            "thermal": {"initial_temperature": ENCLOSURE_WALLS, "gravity": GRAVITY, "expansion": AIR_EXPANSION,
                "faces": held_walls, "surfaces": held(bodies, "Chip", CHIP_TEMPERATURE)}
        }),
        "Natural",
        run,
        faces(bodies, "Chip"),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Session;
    use std::sync::Arc;

    /// Air's volumetric heat capacity ρ c_p, mJ/(mm³ K).
    const AIR_RHO_CP: f64 = 1.2e-12 * 1.005e9;

    /// A session with `design` open and `script` (given its bodies) run.
    fn open(design: &str, script: impl Fn(&[Shape]) -> Value) -> Session {
        let design = fieldcraft_geometry::open_design(&fieldcraft_geometry::script_design(design).unwrap()).unwrap();
        let script = script(&design.bodies);
        let mut s = Session { design: Some(Arc::new(design)), ..Session::default() };
        s.run_script(&script).unwrap();
        s
    }

    /// The model checks clean: ready, no errors or warnings.
    fn assert_ready(s: &mut Session) {
        let check = s.execute("model.check", &json!({})).unwrap();
        assert_eq!((check["ready"].as_bool(), check["errors"].as_u64(), check["warnings"].as_u64()), (Some(true), Some(0), Some(0)), "{check}");
    }

    /// Solves on the CPU: the heat (W) `body` gives the air, and the report.
    fn solve(s: &mut Session, body: &str) -> (f64, Value) {
        let r = s.execute("solve.fluid", &json!({"backend": "cpu"})).unwrap();
        // The report's power is in the model's mW, positive into the wall.
        assert_eq!(r["units"]["power"], "mW", "{r}");
        let q = -r["heat"].as_array().unwrap().iter().filter(|h| h["body"] == body).map(|h| h["final"].as_f64().unwrap()).sum::<f64>() / 1000.0;
        (q, r)
    }

    /// The air stays between the coldest and the hottest held temperature, give or take the
    /// anti-bounce-back wall's undershoot at sharp convex corners.
    fn air_between(r: &Value, cold: f64, hot: f64, undershoot: f64) {
        let (min, max) = (r["temperature"]["min"].as_f64().unwrap(), r["temperature"]["max"].as_f64().unwrap());
        assert!(min > cold - undershoot && max < hot + 0.01, "air {min}..{max} °C: {r}");
    }

    /// A short run on the CPU.
    fn short(spacing: f64, time_step: f64, end_time: f64) -> Run {
        Run { spacing, time_step: Some(time_step), end_time, frames: 2, steady: false }
    }

    /// Within `tolerance` (relative) of `expected`: a band that holds the trunk tests' short runs
    /// (deterministic on the CPU) against an unnoticed change in the solver.
    fn near(value: f64, expected: f64, tolerance: f64) -> bool {
        (value - expected).abs() <= tolerance * expected.abs()
    }

    #[test]
    fn the_cooling_samples_open_ready_to_solve() {
        for name in ["cooling-heat-sink", "cooling-enclosure"] {
            let mut s = Session::default();
            s.execute("file.sample", &json!({"name": name})).unwrap();
            assert_ready(&mut s);
        }
    }

    /// The published cells, 100 steps (25 ms): the start-up's thin thermal layers on the fins
    /// take 6.3 W from the sink.
    #[test]
    fn a_short_heat_sink_run_warms_the_air() {
        let mut s = open(HEAT_SINK_DESIGN, |b| heat_sink_with(b, &short(1.0, 0.00025, 0.025)));
        let (q, r) = solve(&mut s, "HeatSink");
        assert!(near(q, 6.276, 0.05), "{q} W: {r}");
        air_between(&r, SINK_INLET, SINK_TEMPERATURE, 5.0);
    }

    /// 2 mm cells, a quarter second: the chip starts its plume, giving the air 69 mW.
    #[test]
    fn a_short_enclosure_run_takes_heat_from_the_chip() {
        let mut s = open(ENCLOSURE_DESIGN, |b| enclosure_with(b, &short(2.0, 0.0008, 0.25)));
        let (q, r) = solve(&mut s, "Chip");
        assert!(near(q, 0.06915, 0.05), "{q} W: {r}");
        air_between(&r, ENCLOSURE_WALLS, CHIP_TEMPERATURE, 0.5);
    }

    /// The published heat sink (docs/examples.md): 1.93 W at 3 s, a sink-to-inlet resistance of
    /// 13 K/W, and about half (ε 0.48) of the 3.98 W the duct's air could carry if it left at the
    /// sink's temperature. About four minutes on a CPU.
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn the_heat_sink_sheds_its_published_heat() {
        let mut s = open(HEAT_SINK_DESIGN, heat_sink);
        let (q, r) = solve(&mut s, "HeatSink");
        assert!(near(q, 1.925, 0.03), "{q} W: {r}");
        let (min, max) = SINK_DUCT;
        let limit = AIR_RHO_CP * SINK_SPEED * (max[1] - min[1]) * (max[2] - min[2]) * (SINK_TEMPERATURE - SINK_INLET) / 1000.0;
        assert!(near(q / limit, 0.484, 0.03), "ε {}: {r}", q / limit);
        air_between(&r, SINK_INLET, SINK_TEMPERATURE + 0.6, 5.0);
    }

    /// The enclosure to steady state on 2 mm cells: the chip's heat is what the walls take, and
    /// it is within 3 % of the published 1 mm figure (62.7 mW, h ≈ 7.9 W/m²K). The published run
    /// itself (90,000 steps on 144,000 cells) is too long for the heavy job on a CPU.
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn the_enclosure_settles_to_its_published_heat() {
        let run = Run { spacing: 2.0, time_step: Some(0.0008), frames: 2, ..ENCLOSURE };
        let mut s = open(ENCLOSURE_DESIGN, |b| enclosure_with(b, &run));
        let (q, r) = solve(&mut s, "Chip");
        assert_eq!(r["stop"]["reason"], "steady", "{r}");
        let walls = r["heat"].as_array().unwrap().iter().filter(|h| h["body"].is_null()).map(|h| h["final"].as_f64().unwrap()).sum::<f64>() / 1000.0;
        assert!(near(walls, q, 0.002), "walls {walls} W, chip {q} W");
        assert!(near(q, 0.06189, 0.01) && near(q, 0.0627, 0.03), "{q} W: {r}");
        air_between(&r, ENCLOSURE_WALLS, CHIP_TEMPERATURE, 0.5);
    }
}
