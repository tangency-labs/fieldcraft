//! A blower-cooled two-slot graphics card (docs/examples.md, #642): conjugate heat transfer from
//! the GPU die, the memory and the VRMs through a copper base, heat pipes and aluminium fins into
//! the air the blower pushes along the card and out through the I/O bracket.
//!
//! The geometry is `tools/cooling/gpu_card_recipe.py`'s recipe
//! (`crates/engine/samples/cooling-gpu-card.json`): a generic card of our own. Every powered or
//! heat-spreading part is a conducting solid of the fluid domain (#609); the shroud and the bracket
//! are walls.
//!
//! **The blower is a velocity inlet, not a fan curve.** The air box's top face lies just inside the
//! shroud's top, so its only fluid cells are the round intake over the impeller: the air comes in
//! there at the blower's design flow ([`Run::intake_speed`]) and at the case temperature, and
//! leaves through the bracket's slots (the box's x_min face, a convective outlet). The impeller is
//! drawn but not solved. A fan zone with a pressure–flow curve would set the flow from the fins'
//! resistance instead; fan zones are not on the GPU yet (#698).
//!
//! **A large-eddy simulation.** The jets between the VRMs and the fin base and through the fin
//! gaps reach cell Reynolds numbers far above what the lattice holds without a model (#642), so
//! the step runs Smagorinsky's subgrid eddy viscosity, with its eddy diffusivity of heat, and
//! stops on the running means of a statistically steady flow. The solver's Mach limit sets the
//! air's speed on a grid: the 0.36 mm cells the GPU holds today take 2.5 m/s at the intake (a
//! blower on a low setting; a full-speed one moves about twice that). Faster air needs cells of
//! 0.2–0.25 mm, 77–150 million, which the GPU backend can't hold yet (its populations sit in one
//! buffer, and the GB10 allows 2 GB per buffer; #699).

use fieldcraft_geometry::Shape;
use fieldcraft_model::units::{Quantity, Typed};
use serde_json::{Value, json};

use super::sample_value;

/// The geometry: the recipe `tools/cooling/gpu_card_recipe.py` writes.
pub(super) const DESIGN: &str = include_str!("../../samples/cooling-gpu-card.json");

/// How the card is solved.
#[derive(Clone, Copy, Debug)]
pub(super) struct Run {
    /// Cell size, mm.
    pub spacing: f64,
    /// s.
    pub time_step: f64,
    /// s: the most it runs (it stops once steady).
    pub end_time: f64,
    pub frames: u32,
    /// The air's speed into the intake, mm/s: the blower's flow over the intake's area.
    pub intake_speed: f64,
}

/// The published run, on the GPU: 0.36 mm cells (25 million: about the most whose 19 populations
/// fit one 2 GB device buffer, the GB10's limit), the intake at 2.5 m/s (7.1 l/s, 15 cfm) at
/// lattice Mach 0.095, so the jets in the gaps, two to three times faster, stay well inside the
/// lattice's speed. That puts the air's τ at 0.5028, which only the large-eddy simulation allows
/// ([`fieldcraft_lbm::setup::MIN_TAU_TURBULENT`]), and the heat's τ_T at 0.5052, just above its
/// floor; 0.5 s of flow time is about five passes of the air along the card.
pub(super) const CARD: Run = Run { spacing: 0.36, time_step: 7.92e-6, end_time: 0.5, frames: 1, intake_speed: 2500.0 };

/// The case air drawn in through the intake, °C.
const INTAKE: f64 = 38.0;
/// The steady stop rule's tolerance on the running means (relative change per window).
const STEADY_TOLERANCE: f64 = 1e-3; // unit-factors: relative
/// The share of the run discarded as start-up before the running means begin.
const START_UP: f64 = 0.2; // unit-factors: a fraction
/// Smagorinsky's constant: Lilly's 0.17 for isotropic turbulence. The gaps between fins, VRMs
/// and base plates are a few cells across, so the first cells off a wall are not in a log layer
/// and no wall model is used; Van Driest's damping is not on the GPU yet.
const SMAGORINSKY: f64 = 0.17;
/// The turbulent Prandtl number of the eddy diffusivity of heat.
const TURBULENT_PRANDTL: f64 = 0.9;

/// The air box, mm: from the bracket's outer face (the exhaust) to the shroud's inner end wall,
/// between its inner side walls (the box's walls stand in for them), from under the backplate to
/// just inside the shroud's top (37 mm against its 38). The top shell is solid there but for the
/// intake, the blower's inlet; the box stops inside the shroud's filleted top edges, which would
/// otherwise leave the inlet a slit along each edge, forced into a dead-end gap.
const BOX: ([f64; 3], [f64; 3]) = ([-2.0, 2.0, -4.0], [265.0, 109.0, 37.0]);

const DIE_POWER: Typed = Typed::new(275.0, "W");
const MEMORY_POWER: Typed = Typed::new(20.0, "W"); // eight chips, 2.5 W each
const VRM_POWER: Typed = Typed::new(25.0, "W");

/// A conducting material: density, conductivity (one value, or along x, y, z), specific heat.
struct Solid {
    name: &'static str,
    density: Typed,
    conductivity: [Typed; 3],
    specific_heat: Typed,
}

const fn iso(k: f64) -> [Typed; 3] {
    let t = Typed::new(k, "W/(m K)");
    [t, t, t]
}

const COPPER: Solid =
    Solid { name: "Copper", density: Typed::new(8960.0, "kg/m^3"), conductivity: iso(390.0), specific_heat: Typed::new(385.0, "J/(kg K)") };
const ALUMINIUM: Solid =
    Solid { name: "Aluminium", density: Typed::new(2700.0, "kg/m^3"), conductivity: iso(205.0), specific_heat: Typed::new(900.0, "J/(kg K)") };
const SILICON: Solid =
    Solid { name: "Silicon", density: Typed::new(2330.0, "kg/m^3"), conductivity: iso(130.0), specific_heat: Typed::new(700.0, "J/(kg K)") };
/// A multilayer board: copper planes make it conduct far better in-plane than through.
const FR4: Solid = Solid {
    name: "Board",
    density: Typed::new(1850.0, "kg/m^3"),
    conductivity: [Typed::new(25.0, "W/(m K)"), Typed::new(25.0, "W/(m K)"), Typed::new(0.5, "W/(m K)")],
    specific_heat: Typed::new(1100.0, "J/(kg K)"),
};
/// Packaged parts (moulding, substrate and silicon), as one effective material.
const PACKAGE: Solid =
    Solid { name: "Package", density: Typed::new(2000.0, "kg/m^3"), conductivity: iso(20.0), specific_heat: Typed::new(800.0, "J/(kg K)") };
/// A heat pipe as an orthotropic solid: very high along its axis (x), copper-wall across.
const HEAT_PIPE: Solid = Solid {
    name: "HeatPipe",
    density: Typed::new(8000.0, "kg/m^3"),
    conductivity: [Typed::new(5000.0, "W/(m K)"), Typed::new(50.0, "W/(m K)"), Typed::new(50.0, "W/(m K)")],
    specific_heat: Typed::new(385.0, "J/(kg K)"),
};

fn material(s: &Solid) -> Value {
    let k = s.conductivity.map(|k| sample_value(Quantity::ThermalConductivity, k));
    let conductivity = if k[0] == k[1] && k[1] == k[2] { json!(k[0]) } else { json!(k) };
    json!({"command": "model.material.create", "params": {"name": s.name, "density": sample_value(Quantity::Density, s.density),
        "thermal": {"conductivity": conductivity, "specific_heat": sample_value(Quantity::SpecificHeat, s.specific_heat)}}})
}

/// The conducting bodies: (body, material, power).
fn solids() -> Vec<(&'static str, &'static str, Option<Typed>)> {
    vec![
        ("Die", SILICON.name, Some(DIE_POWER)),
        ("Package", PACKAGE.name, None),
        ("Memory", PACKAGE.name, Some(MEMORY_POWER)),
        ("VRM", PACKAGE.name, Some(VRM_POWER)),
        ("Base", COPPER.name, None),
        ("HeatPipes", HEAT_PIPE.name, None),
        ("Fins", ALUMINIUM.name, None),
        ("Board", FR4.name, None),
        ("Backplate", ALUMINIUM.name, None),
    ]
}

pub(super) fn card(bodies: &[Shape]) -> Value {
    card_with(bodies, &CARD)
}

/// The card's commands for `run`: materials, the air box with the conducting solids and the
/// blower's intake, and the fluid step to steady state.
pub(super) fn card_with(_bodies: &[Shape], run: &Run) -> Value {
    let (lo, hi) = BOX;
    let hi = [0, 1, 2].map(|i| super::cooling::snap(lo[i], hi[i], run.spacing));
    let solids = solids();
    let obstacles: Vec<&str> = ["Shroud", "Bracket"].into_iter().chain(solids.iter().map(|s| s.0)).collect();
    let conducting: Vec<Value> = solids
        .iter()
        .map(|(body, material, power)| {
            let mut v = json!({"body": body, "material": material});
            if let Some(p) = power {
                v["power"] = json!(sample_value(Quantity::Power, *p));
            }
            v
        })
        .collect();
    let intake = sample_value(Quantity::Velocity, Typed::new(run.intake_speed, "mm/s"));
    let walls = |t: &str| json!({"type": t});
    let mut commands = vec![
        json!({"command": "model.units", "params": {"system": "mm-t-s-N-MPa"}}),
        json!({"command": "model.material.create", "params": super::cooling::air()}),
    ];
    commands.extend([&COPPER, &ALUMINIUM, &SILICON, &FR4, &PACKAGE, &HEAT_PIPE].into_iter().map(material));
    commands.push(json!({"command": "model.fluid_domain.create", "params": {
        "name": "Case", "material": "Air", "min": lo, "max": hi, "spacing": run.spacing,
        "faces": {"x_min": {"type": "pressure", "pressure": 0, "outflow": "convective"}, "x_max": walls("wall"),
            "y_min": walls("wall"), "y_max": walls("wall"), "z_min": walls("wall"),
            "z_max": {"type": "velocity", "velocity": [0.0, 0.0, -intake], "profile": "uniform"}},
        "obstacles": obstacles,
        "thermal": {"initial_temperature": INTAKE, "faces": {"z_max": {"type": "temperature", "value": INTAKE}}, "solids": conducting}
    }}));
    commands.push(json!({"command": "model.step.create", "params": {"name": "Cool", "procedure": {
        "type": "fluid", "domain": "Case", "end_time": run.end_time, "time_step": run.time_step, "frames": run.frames,
        "collision": "regularized",
        "turbulence": {"type": "smagorinsky", "constant": SMAGORINSKY, "prandtl": TURBULENT_PRANDTL},
        "steady": {"tolerance": STEADY_TOLERANCE, "monitors": [{"type": "mean_temperature"}],
            "averaging": {"start": run.end_time * START_UP}}
    }}}));
    Value::Array(commands)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Session;
    use std::sync::Arc;

    /// A session with the card's design open and `run`'s model.
    fn open(run: &Run) -> Session {
        let design = fieldcraft_geometry::open_design(&fieldcraft_geometry::script_design(DESIGN).unwrap()).unwrap();
        let script = card_with(&design.bodies, run);
        let mut s = Session { design: Some(Arc::new(design)), ..Session::default() };
        s.run_script(&script).unwrap();
        s
    }

    /// The card's model, on 4 mm cells to stay quick: the commands run and the model checks clean.
    #[test]
    fn the_card_builds_and_checks_clean() {
        let mut s = open(&Run { spacing: 4.0, ..CARD });
        let check = s.execute("model.check", &json!({})).unwrap();
        assert_eq!((check["ready"].as_bool(), check["errors"].as_u64()), (Some(true), Some(0)), "{check}");
    }

    /// Probe (not a check): a coarser card on the CPU with the LES on, `FIELDCRAFT_CARD_STEPS`
    /// steps (default 400) on `FIELDCRAFT_CARD_SPACING` mm cells with time step
    /// `FIELDCRAFT_CARD_DT` s and `FIELDCRAFT_CARD_SPEED` mm/s into the intake (default 0.5 mm,
    /// 1.46e-5 s and 1.9 m/s: the published run's τ and lattice Mach on the coarser cells),
    /// printing the speeds, the solids' temperatures and the time per step.
    #[test]
    #[ignore = "measurement"]
    fn probe_the_card_on_the_cpu() {
        let env = |k: &str, d: f64| std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d);
        let (steps, spacing) = (env("FIELDCRAFT_CARD_STEPS", 400.0), env("FIELDCRAFT_CARD_SPACING", 0.5));
        let (time_step, intake_speed) = (env("FIELDCRAFT_CARD_DT", 1.46e-5), env("FIELDCRAFT_CARD_SPEED", 1900.0));
        let run = Run { spacing, time_step, end_time: steps * time_step, frames: 1, intake_speed };
        let mut s = open(&run);
        let t0 = std::time::Instant::now();
        let r = s.execute("solve.fluid", &json!({"backend": "cpu"}));
        let dt = t0.elapsed().as_secs_f64();
        match r {
            Ok(r) => {
                println!(
                    "cells {} steps {} tau {} mach {} max_speed {} final {} ({:.3} s/step)",
                    r["cells"],
                    r["steps"],
                    r["tau"],
                    r["mach"],
                    r["max_speed"],
                    r["final_max_speed"],
                    dt / steps
                );
                for sr in r["solids"].as_array().into_iter().flatten() {
                    println!("{} {}", sr["body"], sr["temperature"]);
                }
                println!("warnings {}", r["warnings"]);
            }
            Err(e) => println!("after {dt:.0} s: {e}"),
        }
    }

    /// The published cells open and check clean.
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn the_card_opens_ready_to_solve() {
        let mut s = open(&CARD);
        let check = s.execute("model.check", &json!({})).unwrap();
        assert_eq!((check["ready"].as_bool(), check["errors"].as_u64()), (Some(true), Some(0)), "{check}");
        println!("check: {check}");
    }
}
