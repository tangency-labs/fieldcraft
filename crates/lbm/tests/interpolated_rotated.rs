//! Interpolated walls on a channel turned against the grid: a plane Poiseuille channel `a:b`
//! against the lattice (26.57° is 2:1, 45° is 1:1), its walls moved by a fraction of a cell, driven
//! by a body force along the channel and periodic all round, no turbulence model. The exact answer
//! is the parabola between the wall planes. Each case runs 10 240 steps and reports the flow rate
//! against the exact one and the drift of the mean density per 10 000 steps; walls that conserve
//! mass hold the density, walls that do not make it drift.
//!
//! `ROTATED_RECORD=1` prints every case's numbers without asserting, to record new ones.
//!
//! The expected numbers below are what the code gives today, recorded as "known": the misses of
//! the linear rule stay visible as numbers, and a change that moves them says so here. The channel
//! of 8 cells across (τ = 0.65, lattice viscosity 0.05) is the coarse grid; a heavy test adds 16.
//! Theory: Bouzidi, Firdaouss & Lallemand, *Momentum transfer of a Boltzmann-lattice fluid with
//! boundaries*, Phys. Fluids 13:3452, 2001.

mod common;

use common::*;
use fieldcraft_lbm::Obstacle;
use fieldcraft_model::{FaceCondition, WallTreatment};

/// The lattice's mean speed the body force is set for.
const MEAN_SPEED: f64 = 0.02;
const STEPS: u64 = 10_240;

/// What a run measured.
#[derive(Debug)]
struct Run {
    /// Flow rate against the exact one, in percent.
    flow_rate: f64,
    /// Drift of the mean density over the fluid cells per 10 000 steps (lattice density 1).
    drift: f64,
    /// Wall links on which no crossing with the surface was found (kept halfway).
    no_crossing: u64,
}

/// The channel `cells` cells across, turned `a:b`, its walls moved `phase` cells, run [`STEPS`]
/// steps with the wall treatment `walls`.
fn rotated_channel((a, b): (usize, usize), cells: usize, phase: f64, walls: WallTreatment) -> Run {
    let nu_lat = 0.05;
    let dx = 1.0 / cells as f64;
    let dt = nu_lat * dx * dx;
    let width = cells as f64;
    let len = ((a * a + b * b) as f64).sqrt();
    let mut m = 1usize;
    while ((a * b * m) as f64 / len) < width + 4.0 {
        m += 1;
    }
    let (nx, ny, nz) = (a * m, b * m, 4usize);
    let pitch = (a * b * m) as f64 / len;
    let dir = [a as f64 / len, b as f64 / len];
    let n = [-dir[1], dir[0]];
    let (lx, ly, lz) = (nx as f64 * dx, ny as f64 * dx, nz as f64 * dx);
    let (p, ho, s0) = (pitch * dx, width * dx, phase * dx);
    let diag = (lx * lx + ly * ly).sqrt();
    // The wall strip and its copies one box apart that reach into the box.
    let mut obstacles: Vec<Obstacle> = Vec::new();
    for k in -4i32..=4 {
        let (lo, hi) = (ho + s0 + k as f64 * p, p + s0 + k as f64 * p);
        if hi < -p + 0.5 * dx || lo > p - 0.5 * dx {
            continue;
        }
        let t = 3.0 * diag;
        let origin = [n[0] * lo - dir[0] * t, n[1] * lo - dir[1] * t, -dx];
        let e = [[2.0 * t * dir[0], 2.0 * t * dir[1], 0.0], [n[0] * (hi - lo), n[1] * (hi - lo), 0.0], [0.0, 0.0, lz + 2.0 * dx]];
        obstacles.push(obstacle(&format!("Wall{}", k + 4), parallelepiped(origin, e), &CUBOID_FACES));
    }
    let mut d = domain([0.0; 3], [lx, ly, lz], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.obstacles = obstacles.iter().map(|o| o.body.clone()).collect();
    // g = 12 ν U / H² gives the exact mean speed (ν = 1 in the model, U in physical units).
    let g = 12.0 * MEAN_SPEED * dx / dt / (ho * ho);
    d.acceleration = [g * dir[0], g * dir[1], 0.0];
    let proc = with_walls(step(STEPS as f64 * dt, dt, None), walls);
    let (_, mut flow) = prepare(&model(1.0, 1.0, d, proc), &obstacles);
    let fluid: Vec<usize> = (0..flow.grid().len()).filter(|&c| flow.is_fluid(c)).collect();
    let mean_density = |flow: &fieldcraft_lbm::Flow| fluid.iter().map(|&c| flow.moments(c).map_or(0.0, |(r, _)| r)).sum::<f64>() / fluid.len() as f64;
    let start = mean_density(&flow);
    for _ in 0..STEPS {
        flow.step().unwrap();
    }
    let drift = (mean_density(&flow) - start) / STEPS as f64 * 10_000.0;
    let along: f64 = fluid.iter().map(|&c| flow.moments(c).map_or(0.0, |(_, u)| u[0] * dir[0] + u[1] * dir[1])).sum();
    let along_cells = ((nx * nx + ny * ny) as f64).sqrt();
    let flow_rate = along / along_cells / nz as f64 / (MEAN_SPEED * width) - 1.0;
    Run { flow_rate: 100.0 * flow_rate, drift, no_crossing: flow.wall_links().map_or(0, |w| w.no_crossing) }
}

/// One expected case: the angle, the grid, the wall offset in cells, the flow-rate error in
/// percent and the density drift per 10 000 steps.
struct Case {
    angle: (usize, usize),
    cells: usize,
    phase: f64,
    flow_rate: f64,
    drift: f64,
}

const fn case(angle: (usize, usize), cells: usize, phase: f64, flow_rate: f64, drift: f64) -> Case {
    Case { angle, cells, phase, flow_rate, drift }
}

/// Known: the linear rule on the turned channel. Flow rates miss by up to 1.6 % on 8 cells and the
/// mean density drifts by up to 2.7e-2 per 10 000 steps.
const LINEAR: [Case; 12] = [
    case((2, 1), 8, 0.0, -0.5419, -1.582e-02),
    case((2, 1), 8, 0.33, -0.7636, -2.669e-02),
    case((2, 1), 8, 0.67, -0.9289, 1.249e-02),
    case((1, 1), 8, 0.0, -1.6298, -6.991e-03),
    case((1, 1), 8, 0.25, -0.9163, -2.230e-02),
    case((1, 1), 8, 0.5, 0.6775, -6.162e-03),
    case((2, 1), 16, 0.0, -0.1126, -1.273e-03),
    case((2, 1), 16, 0.33, -0.1927, -3.099e-03),
    case((2, 1), 16, 0.67, -0.1485, 4.589e-04),
    case((1, 1), 16, 0.0, -0.5987, -1.061e-03),
    case((1, 1), 16, 0.25, -0.5841, -1.155e-03),
    case((1, 1), 16, 0.5, 0.1562, -1.544e-03),
];

fn check(cases: &[&Case], walls: WallTreatment) {
    for c in cases {
        let r = rotated_channel(c.angle, c.cells, c.phase, walls);
        println!(
            "{:?} {} cells, offset {}: flow rate {:+.4} % (expected {:+.4}), drift {:+.3e} (expected {:+.3e})",
            c.angle, c.cells, c.phase, r.flow_rate, c.flow_rate, r.drift, c.drift
        );
        if std::env::var("ROTATED_RECORD").is_ok() {
            // Print the numbers to record, assert nothing.
            continue;
        }
        assert_eq!(r.no_crossing, 0, "{:?} {} {}: links without a crossing", c.angle, c.cells, c.phase);
        assert!((r.flow_rate - c.flow_rate).abs() < 2e-3, "{:?} {} {}: flow rate {} against {}", c.angle, c.cells, c.phase, r.flow_rate, c.flow_rate);
        assert!(
            (r.drift - c.drift).abs() < 1e-5 + 0.01 * c.drift.abs(),
            "{:?} {} {}: drift {} against {}",
            c.angle,
            c.cells,
            c.phase,
            r.drift,
            c.drift
        );
    }
}

/// Two of the coarse cases, the linear rule's numbers (known misses).
#[test]
fn the_linear_rule_on_a_turned_channel_keeps_its_recorded_error() {
    check(&[&LINEAR[2], &LINEAR[3]], WallTreatment::Interpolated);
}

/// All twelve cases, the linear rule's numbers.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_linear_rule_on_a_turned_channel_keeps_its_recorded_error_everywhere() {
    check(&LINEAR.iter().collect::<Vec<_>>(), WallTreatment::Interpolated);
}
