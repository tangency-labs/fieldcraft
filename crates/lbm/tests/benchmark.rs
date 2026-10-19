//! Benchmark suite: comparisons with published numerical reference values (not analytic; the
//! analytic suite in `verification.rs` is the gate). Each case names its source. The cases are
//! long, so they are ignored by default:
//! `cargo test --release -p fieldcraft-lbm --test benchmark -- --ignored --nocapture`.

mod common;

use common::*;
use fieldcraft_lbm::{Limits, solve_fluid};
use fieldcraft_model::{FaceCondition, Procedure, Profile};

/// Benchmark "2D-1" (steady flow around a cylinder at Re = 20): M. Schäfer, S. Turek, *Benchmark
/// computations of laminar flow around a cylinder*, in: Flow Simulation with High-Performance
/// Computers II, NNFM 52:547–566, Vieweg 1996, DOI 10.1007/978-3-322-89849-4_39.
/// Channel 2.2 × 0.41, cylinder of diameter 0.1 centred at (0.2, 0.2), parabolic inlet with peak
/// 0.3 (mean 0.2), ν = 10⁻³, ρ = 1. Reference intervals: drag coefficient 5.57–5.59, lift
/// coefficient 0.0104–0.0110, pressure difference between the cylinder's front and back
/// 0.1172–0.1176. Run quasi-2D: one periodic cell across z.
/// Returns (c_D, c_L, Δp).
fn schafer_turek_2d1(dx: f64, time_step: Option<f64>) -> (f64, f64, f64) {
    let (h, l, d, um, nu) = (0.41, 2.2, 0.1, 0.3, 1e-3);
    let mut dom = domain([0.0; 3], [l, h, dx], dx);
    dom.faces.x_min = FaceCondition::Velocity { velocity: [um, 0.0, 0.0], profile: Profile::Parabolic };
    dom.faces.x_max = FaceCondition::Pressure { pressure: 0.0 };
    dom.faces.z_min = FaceCondition::Periodic {};
    dom.faces.z_max = FaceCondition::Periodic {};
    dom.obstacles = vec!["Cylinder".into()];
    let cyl = obstacle("Cylinder", cylinder([0.2, 0.2], d / 2.0, -0.1, 0.1, 1024), &["side", "bottom", "top"]);
    let procedure = Procedure::Fluid { domain: "Flow".into(), end_time: Some(10.0), steps: None, time_step, frames: 1, magic: None, steady: None };
    let m = model(1.0, nu, dom, procedure);
    let r = solve_fluid(&m, None, std::slice::from_ref(&cyl), &Limits::default()).unwrap();
    let side = r.forces.iter().find(|f| f.name == "side").unwrap();
    // The mean over the last tenth of the run, per unit depth.
    let n = side.force.len();
    let tail = &side.force[n - n / 10..];
    let mean = |a: usize| tail.iter().map(|v| v[a]).sum::<f64>() / tail.len() as f64 / dx;
    let mean_u = 2.0 / 3.0 * um;
    let coefficient = |f: f64| 2.0 * f / (mean_u * mean_u * d);
    // Pressure in the fluid cells just in front of and behind the cylinder on its axis
    // (the average of the two cell rows about y = 0.2).
    let p = |x: f64| {
        let at = |y: f64| f64::from(r.probe(1, [x, y, 0.5 * dx]).unwrap().2);
        0.5 * (at(0.2 - 0.5 * dx) + at(0.2 + 0.5 * dx))
    };
    let dp = p(0.15 - 0.5 * dx) - p(0.25 + 0.5 * dx);
    println!("Schäfer–Turek 2D-1, Δx = {dx}: {:?}; warnings {:?}", r.stats, r.warnings);
    (coefficient(mean(0)), coefficient(mean(1)), dp)
}

#[test]
#[ignore]
fn schafer_turek_2d1_cylinder_drag_and_lift() {
    // Lattice velocity 0.025 (Mach 0.043): the inlet's pressure, about 0.3, changes the density
    // by Δp / (ρ c_s²) in lattice units, so the time step shrinks with the spacing squared.
    for dx in [0.01, 0.005] {
        let dt = 0.025 * dx / 0.3;
        let (cd, cl, dp) = schafer_turek_2d1(dx, Some(dt * (dx / 0.01)));
        println!("Schäfer–Turek 2D-1, D/Δx = {}: c_D = {cd:.4} (5.57–5.59), c_L = {cl:.4} (0.0104–0.0110), Δp = {dp:.4} (0.1172–0.1176)", 0.1 / dx);
        if dx <= 0.005 {
            // A staircase cylinder of 20 cells (measured: c_D 5.773, +3.4%; c_L 0.0107, inside
            // the interval). Δp is read in the cells half a spacing off the staircase, where the
            // pressure changes fast, so it is reported, not asserted (measured 0.108).
            assert!((cd / 5.58 - 1.0).abs() < 0.05, "{cd}");
            assert!((cl / 0.0107 - 1.0).abs() < 0.2, "{cl}");
        }
    }
}
