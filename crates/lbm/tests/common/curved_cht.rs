//! Original model of Rahman et al., JNAME 5(2) (2008) 37–46,
//! doi:10.3329/jname.v5i2.2504, Fig. 1 and Table 1. The finite-conductivity cylinder is
//! passive: no prescribed surface temperature or internal source. The reference integrates
//! the heated RIGHT cavity wall (length 0.9 L), with theta = (T - T_in)/(T_hot - T_in).
//! The authors describe their own Galerkin finite-element/Newton solver, not a commercial tool.
//!
//! Box boundary conditions require short entrance/exit passages to represent the 0.1 L slots.
//! Their padding is an explicit model difference to measure separately, not curved-interface
//! error. All grids are multiples of ten, so the slots and cavity walls remain aligned.

use fieldcraft_lbm::{FluidResult, Obstacle, Setup};
use fieldcraft_model::{
    FaceCondition, FlowMonitor, FluidBoundary, FluidThermal, GeometryRef, Model, Outflow, Profile, Steady, SurfaceThermalBoundary, ThermalCondition,
};

use super::*;

const PERIODIC: FaceCondition = FaceCondition::Periodic {};

pub const RE: f64 = 100.0;
pub const PR: f64 = 0.71;
pub const RI: f64 = 1.0;
pub const DIAMETER: f64 = 0.2;
pub const CONDUCTIVITY_RATIO: f64 = 5.0;
/// The paper's finest (7512-element) value; its chosen 5918-element grid gives 4.832878.
pub const REFERENCE_NU: f64 = 4.842449;
pub const REFERENCE_MEAN_T: f64 = 0.197227;
/// Predeclared accuracy target, separate from the reference's approximately 0.2% grid spread.
pub const TARGET_RELATIVE: f64 = 0.01;

pub fn model_case(n: usize, padding: usize, segments: usize, end_time: f64, lattice_speed: f64) -> (Model, Vec<Obstacle>) {
    assert!([20, 40, 80, 160].contains(&n));
    assert!((1..=4).contains(&padding));
    assert!([256, 512, 1024].contains(&segments));
    let dx = 1.0 / n as f64;
    let extra = padding as f64 * dx;
    let nu = 1.0 / RE;
    let alpha = nu / PR;
    let mut d = domain([-extra, 0.0, 0.0], [1.0 + extra, 1.0, dx], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [1.0, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective };
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d.obstacles = vec!["LeftWall".into(), "RightWall".into(), "Cylinder".into()];
    let hot = GeometryRef::face("RightWall", "x_min").unwrap();
    let mut thermal =
        FluidThermal { initial_temperature: kelvin(0.0), gravity: [0.0, -1000.0 * RI, 0.0], expansion: Some(0.001), ..Default::default() };
    thermal.faces.x_min = ThermalCondition::Temperature { value: kelvin(0.0) };
    thermal.surfaces = vec![SurfaceThermalBoundary { face: hot.clone(), condition: ThermalCondition::Temperature { value: kelvin(1.0) } }];
    thermal.solids = vec![conducting("Cylinder", 0.0)];
    d.thermal = Some(thermal);
    let mut procedure = step(end_time, lattice_speed * dx, None);
    if let Procedure::Fluid { steady, .. } = &mut procedure {
        *steady = Some(Steady {
            tolerance: 1e-5,
            window: Some(5.0),
            min_time: Some(20.0),
            divergence: None,
            monitors: vec![FlowMonitor::Heat { faces: vec![FluidBoundary::Surface(hot)] }, FlowMonitor::MeanTemperature {}],
        });
    }
    let model = cht_model([1.0, nu, alpha], (1.0, [CONDUCTIVITY_RATIO * alpha; 3], 1.0), d, procedure);
    let obstacles = vec![
        obstacle("LeftWall", cuboid([-2.0, 0.1, -1.0], [0.0, 2.0, 1.0]), &CUBOID_FACES),
        obstacle("RightWall", cuboid([1.0, -1.0, -1.0], [2.0, 0.9, 1.0]), &CUBOID_FACES),
        obstacle("Cylinder", cylinder([0.5, 0.5], DIAMETER * 0.5, -1.0, 1.0, segments), &["side", "bottom", "top"]),
    ];
    (model, obstacles)
}

#[derive(Clone, Copy, Debug)]
pub struct Observation {
    pub nusselt: f64,
    pub mean_fluid: f64,
    pub solid_area: f64,
    pub heat_spread: f64,
}

pub fn observation(setup: &Setup, result: &FluidResult, padding: usize, segments: usize) -> Observation {
    let hot = result.heat.iter().find(|h| h.body.as_deref() == Some("RightWall") && h.name == "x_min").unwrap();
    let tail = &hot.power[hot.power.len() - (hot.power.len() / 10).max(1)..];
    let heat = -tail.iter().sum::<f64>() / tail.len() as f64;
    let alpha = 1.0 / (RE * PR);
    let dx = result.grid.spacing;
    let nusselt = heat / (alpha * 0.9 * dx);
    let frame = result.frames.last().unwrap();
    let mut fluid = Vec::new();
    let mut solid = Vec::new();
    for (c, &t) in frame.temperature.iter().enumerate() {
        let [x, y, _] = result.grid.centre(c);
        if x > 0.0 && x < 1.0 && y > 0.0 && y < 1.0 {
            if !result.solid[c] {
                fluid.push(t as f64);
            } else if (x - 0.5).powi(2) + (y - 0.5).powi(2) < (DIAMETER * 0.5).powi(2) {
                solid.push(t as f64);
            }
        }
    }
    let mean = fluid.iter().sum::<f64>() / fluid.len() as f64;
    let low = solid.iter().copied().fold(f64::INFINITY, f64::min);
    let high = solid.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let stopped = result.steady.as_ref().unwrap();
    let heat_spread = (tail.iter().copied().fold(f64::NEG_INFINITY, f64::max) - tail.iter().copied().fold(f64::INFINITY, f64::min)) / heat;
    println!(
        "curved_cht n={} padding={padding} segments={segments} Nu={nusselt:.9} reference={REFERENCE_NU:.6} error_percent={:.6} theta_fluid={mean:.9} theta_solid={low:.9}..{high:.9} cylinder_area={:.9} stop={} time={:.6} steps={} Mach={:.6} tau_T={:.6} backend={}",
        result.grid.cells[1],
        100.0 * (nusselt / REFERENCE_NU - 1.0),
        solid.len() as f64 * dx * dx,
        stopped.reason.name(),
        result.stats.end_time,
        result.stats.steps,
        result.stats.mach,
        setup.thermal.as_ref().unwrap().tau,
        result.stats.backend
    );
    assert!(nusselt.is_finite() && nusselt > 0.0 && mean > 0.0 && mean < 1.0);
    assert!(low >= -1e-5 && high <= 1.0 + 1e-5 && high > low, "a passive conducting solid has a solved, nonuniform temperature");
    Observation { nusselt, mean_fluid: mean, solid_area: solid.len() as f64 * dx * dx, heat_spread }
}
