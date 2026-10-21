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
//! resistance instead; fan zones are not on the GPU yet, and the card's air speeds need it.
//!
//! **Real air speeds need fine cells.** The lattice needs its relaxation time above
//! [`fieldcraft_lbm::setup::MIN_TAU`] (the regularised collision above
//! [`fieldcraft_lbm::setup::REGULARIZED_TAU_LIMIT`] at the open faces), which caps the cell
//! Reynolds number at about 60–90, and the solver's Mach limit sets the rest: the 0.36 mm cells
//! the GPU holds today take 2.5 m/s at the intake (a blower on a low setting; a full-speed one
//! moves about twice that). Faster air needs cells of 0.2–0.25 mm, 77–150 million, which the GPU
//! backend can't hold yet (its populations sit in one buffer, and the GB10 allows 2 GB per buffer).

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

/// The published run, on the GPU: 0.36 mm cells (25 million: about the most whose 19 populations fit
/// one 2 GB device buffer, the GB10's limit), τ ≈ 0.508 for air (the regularised collision's
/// limit at the open faces), the intake at 2.5 m/s (7.1 l/s, 15 cfm; lattice Mach 0.28, near the
/// solver's 0.3), until steady.
pub(super) const CARD: Run = Run { spacing: 0.36, time_step: 2.32e-5, end_time: 1.0, frames: 1, intake_speed: 2500.0 };

/// The case air drawn in through the intake, °C.
const INTAKE: f64 = 38.0;
/// The steady stop rule's tolerance on the solids' heat and the flow (relative change per window).
const STEADY_TOLERANCE: f64 = 1e-3; // unit-factors: relative

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
        "steady": {"tolerance": STEADY_TOLERANCE, "monitors": [{"type": "mean_temperature"}]}
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
