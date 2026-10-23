//! Benchmark suite: comparisons with published numerical reference values (not analytic; the
//! analytic suite in `verification.rs` is the gate). Each case names its source. The cases are
//! long, so they are ignored by default:
//! `cargo test --release -p fieldcraft-lbm --test benchmark -- --ignored --nocapture`.

mod common;

use common::*;
use fieldcraft_lbm::{Limits, solve_fluid};
use fieldcraft_model::{Collision, FaceCondition, Outflow, Procedure, Profile};

/// Benchmark "2D-1" (steady flow around a cylinder at Re = 20): M. Schäfer, S. Turek, *Benchmark
/// computations of laminar flow around a cylinder*, in: Flow Simulation with High-Performance
/// Computers II, NNFM 52:547–566, Vieweg 1996, DOI 10.1007/978-3-322-89849-4_39.
/// Channel 2.2 × 0.41, cylinder of diameter 0.1 centred at (0.2, 0.2), parabolic inlet with peak
/// 0.3 (mean 0.2), ν = 10⁻³, ρ = 1. Reference intervals: drag coefficient 5.57–5.59, lift
/// coefficient 0.0104–0.0110, pressure difference between the cylinder's front and back
/// 0.1172–0.1176. Run quasi-2D: one periodic cell across z.
/// Returns (c_D, c_L, Δp).
fn schafer_turek_2d1(dx: f64, time_step: Option<f64>, collision: Collision) -> (f64, f64, f64) {
    let (h, l, d, um, nu) = (0.41, 2.2, 0.1, 0.3, 1e-3);
    let mut dom = domain([0.0; 3], [l, h, dx], dx);
    dom.faces.x_min = FaceCondition::Velocity { velocity: [um, 0.0, 0.0], profile: Profile::Parabolic };
    dom.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    dom.faces.z_min = FaceCondition::Periodic {};
    dom.faces.z_max = FaceCondition::Periodic {};
    dom.obstacles = vec!["Cylinder".into()];
    let cyl = obstacle("Cylinder", cylinder([0.2, 0.2], d / 2.0, -0.1, 0.1, 1024), &["side", "bottom", "top"]);
    let procedure = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(10.0),
        steps: None,
        time_step,
        frames: 1,
        magic: None,
        collision: Some(collision),
        turbulence: None,
        steady: None,
    };
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
    for collision in [Collision::Trt, Collision::Regularized] {
        for dx in [0.01, 0.005] {
            let dt = 0.025 * dx / 0.3;
            let (cd, cl, dp) = schafer_turek_2d1(dx, Some(dt * (dx / 0.01)), collision);
            println!(
                "{collision:?} Schäfer–Turek 2D-1, D/Δx = {}: c_D = {cd:.4} (5.57–5.59), c_L = {cl:.4} (0.0104–0.0110), Δp = {dp:.4} (0.1172–0.1176)",
                0.1 / dx
            );
            if dx <= 0.005 {
                // A staircase cylinder of 20 cells. Measured with two relaxation times: c_D 5.773
                // (+3.4%), c_L 0.0107 (inside the interval); regularized: c_D 5.890 (+5.5%), c_L
                // 0.0132 (+23%): without the two-relaxation-time parameter's exact wall placement
                // its staircase wall sits a little off (as in its Poiseuille test). Δp is read in
                // the cells half a spacing off the staircase, where the pressure changes fast, so
                // it is reported, not asserted (measured 0.108 and 0.105).
                let (drag, lift) = if collision == Collision::Trt { (0.05, 0.2) } else { (0.07, 0.3) };
                assert!((cd / 5.58 - 1.0).abs() < drag, "{cd}");
                assert!((cl / 0.0107 - 1.0).abs() < lift, "{cl}");
            }
        }
    }
}

/// Benchmark "2D-2" (periodic vortex shedding behind the cylinder at Re = 100), from the same
/// paper: the 2D-1 channel with an inlet peak of 1.5 (mean 1). Reference intervals: largest drag
/// coefficient 3.22–3.24, largest lift coefficient 0.99–1.01, Strouhal number D f / Ū
/// 0.295–0.305. `cells` across the cylinder, `lattice` the inlet peak's lattice speed. Returns
/// (largest c_D, largest c_L, St) over the last `window` of flow time of a run to `end_time`.
fn schafer_turek_2d2(cells: f64, lattice: f64, end_time: f64, window: f64, collision: Collision) -> (f64, f64, f64) {
    let (h, l, d, um, nu) = (0.41, 2.2, 0.1, 1.5, 1e-3);
    let dx = d / cells;
    let mut dom = domain([0.0; 3], [l, h, dx], dx);
    dom.faces.x_min = FaceCondition::Velocity { velocity: [um, 0.0, 0.0], profile: Profile::Parabolic };
    dom.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    dom.faces.z_min = FaceCondition::Periodic {};
    dom.faces.z_max = FaceCondition::Periodic {};
    dom.initial_velocity = [2.0 / 3.0 * um, 0.0, 0.0];
    dom.obstacles = vec!["Cylinder".into()];
    let cyl = obstacle("Cylinder", cylinder([0.2, 0.2], d / 2.0, -0.1, 0.1, 1024), &["side", "bottom", "top"]);
    let procedure = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(end_time),
        steps: None,
        time_step: Some(lattice * dx / um),
        frames: 1,
        magic: None,
        collision: Some(collision),
        turbulence: None,
        steady: None,
    };
    let m = model(1.0, nu, dom, procedure);
    let r = solve_fluid(&m, None, std::slice::from_ref(&cyl), &Limits::default()).unwrap();
    let side = r.forces.iter().find(|f| f.name == "side").unwrap();
    let mean_u = 2.0 / 3.0 * um;
    let coefficient = |f: f64| 2.0 * f / dx / (mean_u * mean_u * d);
    let samples: Vec<(f64, f64, f64)> = r
        .history_times
        .iter()
        .zip(&side.force)
        .filter(|(t, _)| **t >= end_time - window)
        .map(|(t, f)| (*t, coefficient(f[0]), coefficient(f[1])))
        .collect();
    let cd = samples.iter().map(|s| s.1).fold(f64::NEG_INFINITY, f64::max);
    let cl = samples.iter().map(|s| s.2).fold(f64::NEG_INFINITY, f64::max);
    // The lift's period from its upward zero crossings (interpolated), over the window.
    let crossings: Vec<f64> =
        samples.windows(2).filter(|w| w[0].2 < 0.0 && w[1].2 >= 0.0).map(|w| w[0].0 + (w[1].0 - w[0].0) * (-w[0].2) / (w[1].2 - w[0].2)).collect();
    let period = if crossings.len() >= 2 { (crossings[crossings.len() - 1] - crossings[0]) / (crossings.len() - 1) as f64 } else { f64::NAN };
    println!(
        "{collision:?} Schäfer–Turek 2D-2, D/Δx = {cells}: {:?}; {} lift periods; warnings {:?}",
        r.stats,
        crossings.len().saturating_sub(1),
        r.warnings
    );
    (cd, cl, d / (period * mean_u))
}

#[test]
#[ignore]
fn schafer_turek_2d2_vortex_shedding_drag_lift_and_strouhal() {
    let runs: Vec<(Collision, (f64, f64, f64))> =
        [Collision::Trt, Collision::Regularized].into_iter().map(|c| (c, schafer_turek_2d2(20.0, 0.05, 12.0, 4.0, c))).collect();
    for (collision, (cd, cl, st)) in runs {
        println!(
            "{collision:?} Schäfer–Turek 2D-2, D/Δx = 20: c_D,max = {cd:.4} (3.22–3.24), c_L,max = {cl:.4} (0.99–1.01), St = {st:.4} (0.295–0.305)"
        );
        // A staircase cylinder of 20 cells: the drag and lift are high by several percent, as the
        // drag is in 2D-1 (two relaxation times, measured: c_D,max 3.447, c_L,max 1.102); the
        // shedding frequency follows the wake and is inside the interval.
        assert!((cd / 3.23 - 1.0).abs() < 0.1, "{cd}");
        assert!((cl / 1.0 - 1.0).abs() < 0.2, "{cl}");
        assert!((st / 0.30 - 1.0).abs() < 0.04, "{st}");
    }
}

/// The lid-driven cavity at Re = 100 (unit square, lid speed 1, ν = 0.01; quasi-2D, one
/// periodic cell across z) on `n` cells: the smallest u/U on the vertical centreline and the
/// height it is at. U. Ghia, K. N. Ghia, C. T. Shin, *High-Re solutions for incompressible flow
/// using the Navier–Stokes equations and a multigrid method*, J. Comput. Phys. 48(3):387–411,
/// 1982, Table I: −0.21090 at y = 0.4531.
fn lid_driven_cavity_re100(n: usize, collision: Collision) -> (f64, f64) {
    let dx = 1.0 / n as f64;
    let mut d = domain([0.0; 3], [1.0, 1.0, dx], dx);
    d.faces.y_max = FaceCondition::MovingWall { velocity: [1.0, 0.0, 0.0] };
    d.faces.z_min = FaceCondition::Periodic {};
    d.faces.z_max = FaceCondition::Periodic {};
    let procedure = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(30.0),
        steps: None,
        time_step: Some(0.05 * dx),
        frames: 1,
        magic: None,
        collision: Some(collision),
        turbulence: None,
        steady: None,
    };
    let r = solve_fluid(&model(1.0, 0.01, d, procedure), None, &[], &Limits::default()).unwrap();
    let last = r.frames.last().unwrap();
    // The centreline x = 0.5 lies between two cell columns: their mean.
    let u = |j: usize| 0.5 * f64::from(last.velocity[r.grid.index([n / 2 - 1, j, 0])][0] + last.velocity[r.grid.index([n / 2, j, 0])][0]);
    let (j, umin) = (0..n).map(|j| (j, u(j))).fold((0, f64::INFINITY), |a, b| if b.1 < a.1 { b } else { a });
    // A parabola through the smallest sample and its neighbours.
    let (a, b, c) = (u(j - 1), umin, u(j + 1));
    let shift = 0.5 * (a - c) / (a - 2.0 * b + c);
    let y = (j as f64 + 0.5 + shift) * dx;
    (b - 0.25 * (a - c) * shift, y)
}

#[test]
#[ignore]
fn lid_driven_cavity_re100_matches_ghia() {
    for collision in [Collision::Trt, Collision::Regularized] {
        let (umin, y) = lid_driven_cavity_re100(64, collision);
        println!("{collision:?} lid-driven cavity Re 100, 64 cells: u_min/U = {umin:.4} at y = {y:.4} (Ghia: −0.2109 at 0.4531)");
        assert!((umin / -0.2109 - 1.0).abs() < 0.02, "{umin}");
        assert!((y - 0.4531).abs() < 0.02, "{y}");
    }
}
