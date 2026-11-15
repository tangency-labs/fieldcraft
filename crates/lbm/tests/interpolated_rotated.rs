//! Interpolated walls on a channel turned against the grid: a plane Poiseuille channel `a:b`
//! against the lattice (26.57° is 2:1, 45° is 1:1), its walls moved by a fraction of a cell, driven
//! by a body force along the channel and periodic all round, no turbulence model. The exact answer
//! is the parabola between the wall planes. Each case runs 10 240 steps and reports the flow rate
//! error (the sum of the velocity along the channel over the fluid cells' centres against the sum of
//! the exact parabola over the same centres, less 1, so the midpoint sum's own error is not in it)
//! and the drift of the mean density per 10 000 steps; walls that conserve mass hold the density,
//! walls that do not make it drift.
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
    // The exact parabola at a cell centre: its distance across the channel from the wall plane at s0, in cells, folded
    // into one period of the strips.
    let grid = *flow.grid();
    let exact: f64 = fluid
        .iter()
        .map(|&c| {
            let [x, y, _] = grid.centre(c);
            let s = ((n[0] * x + n[1] * y - s0) / dx).rem_euclid(pitch);
            if s < width { 6.0 * MEAN_SPEED * s * (width - s) / (width * width) } else { 0.0 }
        })
        .sum();
    let flow_rate = along / exact - 1.0;
    Run { flow_rate: 100.0 * flow_rate, drift }
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

/// Known: the linear rule on the turned channel. Flow rates miss by up to 1.9 % on 8 cells and the
/// mean density drifts by up to 3.3e-2 per 10 000 steps.
const LINEAR: [Case; 12] = [
    case((2, 1), 8, 0.0, -0.8242, -3.297e-2),
    case((2, 1), 8, 0.33, -1.0281, 1.252e-3),
    case((2, 1), 8, 0.67, -1.0976, 7.518e-3),
    case((1, 1), 8, 0.0, -1.9410, -7.201e-3),
    case((1, 1), 8, 0.25, -1.6530, -2.067e-2),
    case((1, 1), 8, 0.5, 1.8458, -3.923e-3),
    case((2, 1), 16, 0.0, -0.1491, -6.662e-3),
    case((2, 1), 16, 0.33, -0.2496, 1.422e-3),
    case((2, 1), 16, 0.67, -0.2119, 6.132e-4),
    case((1, 1), 16, 0.0, -0.8731, -1.176e-3),
    case((1, 1), 16, 0.25, -0.8478, -1.182e-3),
    case((1, 1), 16, 0.5, 0.2492, -1.434e-3),
];

/// Runs every case, prints what it measured, and fails at the end with every case that left its
/// recorded numbers (so one run gives all the numbers to record).
fn check(cases: &[&Case], walls: WallTreatment) {
    let mut failures = Vec::new();
    for c in cases {
        let r = rotated_channel(c.angle, c.cells, c.phase, walls);
        println!(
            "{:?} {} cells, offset {}: flow rate {:+.4} % (expected {:+.4}), drift {:+.3e} (expected {:+.3e})",
            c.angle, c.cells, c.phase, r.flow_rate, c.flow_rate, r.drift, c.drift
        );
        if (r.flow_rate - c.flow_rate).abs() >= 2e-3 {
            failures.push(format!("{:?} {} {}: flow rate {} against {}", c.angle, c.cells, c.phase, r.flow_rate, c.flow_rate));
        }
        if (r.drift - c.drift).abs() >= 1e-5 + 0.01 * c.drift.abs() {
            failures.push(format!("{:?} {} {}: drift {} against {}", c.angle, c.cells, c.phase, r.drift, c.drift));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
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
