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
use fieldcraft_lbm::setup::MIN_TAU;
use fieldcraft_lbm::{Obstacle, WallLinks};
use fieldcraft_model::{Collision, FaceCondition, WallTreatment};

/// The lattice's mean speed the body force is set for.
const MEAN_SPEED: f64 = 0.02;
const STEPS: u64 = 10_240;

/// What a run measured.
#[derive(Debug)]
struct Run {
    /// Flow rate against the exact one, in percent.
    flow_rate: f64,
    /// The relative L2 error of the velocity along the channel over the fluid cells' centres.
    l2: f64,
    /// Drift of the mean density over the fluid cells per 10 000 steps (lattice density 1).
    drift: f64,
    /// The largest speed of any fluid cell at the end of the run (lattice units).
    max_speed: f64,
    /// The largest exact parabola over the same fluid cells' centres (lattice units): the discrete peak.
    exact_peak: f64,
    /// How the wall links bounced back.
    links: WallLinks,
}

/// The collision and the lattice viscosity a channel runs with, and its mean speed (lattice units).
#[derive(Debug, Clone, Copy)]
struct Regime {
    collision: Collision,
    nu_lat: f64,
    mean_speed: f64,
}

/// The regime of the recorded cases: TRT at τ⁺ 0.65, the mean speed [`MEAN_SPEED`].
const BASE: Regime = Regime { collision: Collision::Trt, nu_lat: 0.05, mean_speed: MEAN_SPEED };

/// The channel `width` cells across, turned `a:b`, its walls moved `phase` cells, run `steps`
/// steps with the wall treatment `walls`.
fn rotated_channel((a, b): (usize, usize), width: f64, phase: f64, walls: WallTreatment, steps: u64) -> Run {
    run_channel((a, b), width, phase, walls, steps, BASE)
}

/// [`rotated_channel`] in the regime `regime`.
fn run_channel((a, b): (usize, usize), width: f64, phase: f64, walls: WallTreatment, steps: u64, regime: Regime) -> Run {
    let nu_lat = regime.nu_lat;
    let dx = 1.0 / width;
    let dt = nu_lat * dx * dx;
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
    let g = 12.0 * regime.mean_speed * dx / dt / (ho * ho);
    d.acceleration = [g * dir[0], g * dir[1], 0.0];
    let proc = with_walls(collide(steps as f64 * dt, dt, None, regime.collision), walls);
    let (_, mut flow) = prepare(&model(1.0, 1.0, d, proc), &obstacles);
    let fluid: Vec<usize> = (0..flow.grid().len()).filter(|&c| flow.is_fluid(c)).collect();
    let mean_density = |flow: &fieldcraft_lbm::Flow| fluid.iter().map(|&c| flow.moments(c).map_or(0.0, |(r, _)| r)).sum::<f64>() / fluid.len() as f64;
    let start = mean_density(&flow);
    for _ in 0..steps {
        flow.step().unwrap();
    }
    let drift = (mean_density(&flow) - start) / steps as f64 * 10_000.0;
    let max_speed = fluid.iter().map(|&c| flow.moments(c).map_or(0.0, |(_, u)| (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt())).fold(0.0, f64::max);
    // Each fluid cell's velocity along the channel and the exact parabola at its centre: the distance across the
    // channel from the wall plane at s0, in cells, folded into one period of the strips.
    let grid = *flow.grid();
    let pairs: Vec<(f64, f64)> = fluid
        .iter()
        .map(|&c| {
            let u = flow.moments(c).map_or(0.0, |(_, u)| u[0] * dir[0] + u[1] * dir[1]);
            let [x, y, _] = grid.centre(c);
            let s = ((n[0] * x + n[1] * y - s0) / dx).rem_euclid(pitch);
            (u, if s < width { 6.0 * regime.mean_speed * s * (width - s) / (width * width) } else { 0.0 })
        })
        .collect();
    let exact_peak = pairs.iter().map(|&(_, x)| x).fold(0.0, f64::max);
    let (along, exact) = pairs.iter().fold((0.0, 0.0), |(a, e), (u, x)| (a + u, e + x));
    let (e2, x2) = pairs.iter().fold((0.0, 0.0), |(e2, x2), (u, x)| (e2 + (u - x) * (u - x), x2 + x * x));
    Run { flow_rate: 100.0 * (along / exact - 1.0), l2: (e2 / x2).sqrt(), drift, max_speed, exact_peak, links: flow.wall_links().unwrap_or_default() }
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

/// Known: the linear rule on the turned channel. Flow rates miss by up to 1.4 % on 8 cells and the
/// mean density drifts by up to 2.7e-2 per 10 000 steps. The 2:1 rows at 0.67 (and 0.33 on 16 cells)
/// have links at q < 0.01. Their flow rates read 0.13, 0.06 and 0.006 points closer to exact while q
/// was held at 0.01 (before #1059): that floor widened the channel by a hundredth of a link, a
/// first-order error of the opposite sign to the linear rule's own. These are the rule's numbers.
const LINEAR: [Case; 12] = [
    case((2, 1), 8, 0.0, -0.6864, -1.582e-2),
    case((2, 1), 8, 0.33, -0.7247, -2.669e-2),
    case((2, 1), 8, 0.67, -0.8406, 1.197e-2),
    case((1, 1), 8, 0.0, -1.3624, -6.991e-3),
    case((1, 1), 8, 0.25, -1.1782, -2.230e-2),
    case((1, 1), 8, 0.5, 0.9205, -6.162e-3),
    case((2, 1), 16, 0.0, -0.1402, -1.273e-3),
    case((2, 1), 16, 0.33, -0.1664, -2.921e-3),
    case((2, 1), 16, 0.67, -0.1734, 3.859e-4),
    case((1, 1), 16, 0.0, -0.5402, -1.061e-3),
    case((1, 1), 16, 0.25, -0.5332, -1.155e-3),
    case((1, 1), 16, 0.5, 0.0983, -1.544e-3),
];

/// Recorded: the quadratic rule on the same channels. Its accuracy is checked by its order of convergence
/// ([`the_quadratic_rule_on_a_turned_channel_converges_at_second_order`]); these numbers only freeze
/// what it gives. The density still drifts: the walls do not conserve mass.
const QUADRATIC: [Case; 12] = [
    case((2, 1), 8, 0.0, 0.0545, -5.823e-3),
    case((2, 1), 8, 0.33, 0.3896, -1.104e-2),
    case((2, 1), 8, 0.67, 0.6444, -7.831e-3),
    case((1, 1), 8, 0.0, 0.5031, -8.965e-3),
    case((1, 1), 8, 0.25, -0.2425, -1.951e-2),
    case((1, 1), 8, 0.5, 0.5730, -3.644e-3),
    case((2, 1), 16, 0.0, 0.0200, 7.439e-4),
    case((2, 1), 16, 0.33, 0.1209, 2.029e-4),
    case((2, 1), 16, 0.67, 0.1333, -1.415e-3),
    case((1, 1), 16, 0.0, 0.1126, -9.363e-4),
    case((1, 1), 16, 0.25, 0.1057, -1.005e-3),
    case((1, 1), 16, 0.5, -0.0427, -4.774e-4),
];

/// Runs every case, prints what it measured, and fails at the end with every case that left its
/// recorded numbers (so one run gives all the numbers to record).
fn check(cases: &[&Case], walls: WallTreatment) {
    let mut failures = Vec::new();
    for c in cases {
        let r = rotated_channel(c.angle, c.cells as f64, c.phase, walls, STEPS);
        println!(
            "{:?} {} cells, offset {}: flow rate {:+.4} % (expected {:+.4}), drift {:+.3e} (expected {:+.3e})",
            c.angle, c.cells, c.phase, r.flow_rate, c.flow_rate, r.drift, c.drift
        );
        if r.links.no_crossing != 0 {
            failures.push(format!("{:?} {} {}: {} links without a crossing", c.angle, c.cells, c.phase, r.links.no_crossing));
        }
        if (r.flow_rate - c.flow_rate).abs() >= 2e-3 {
            failures.push(format!("{:?} {} {}: flow rate {} against {}", c.angle, c.cells, c.phase, r.flow_rate, c.flow_rate));
        }
        if (r.drift - c.drift).abs() >= 1e-5 + 0.01 * c.drift.abs() {
            failures.push(format!("{:?} {} {}: drift {} against {}", c.angle, c.cells, c.phase, r.drift, c.drift));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// A wall plane through cell centres on a spacing that is not a power of two (45°, offset 0, 12
/// cells): every wall link finds its crossing.
#[test]
fn a_wall_through_cell_centres_leaves_no_link_without_a_crossing() {
    let r = rotated_channel((1, 1), 12.0, 0.0, WallTreatment::Interpolated, 1);
    assert_eq!((r.links.no_crossing, r.links.interpolated), (0, 1472), "{:?}", r.links);
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

/// The quadratic rule, two coarse cases, its recorded numbers.
#[test]
fn the_quadratic_rule_on_a_turned_channel_keeps_its_recorded_error() {
    check(&[&QUADRATIC[2], &QUADRATIC[3]], WallTreatment::InterpolatedQuadratic);
}

/// All twelve cases, the quadratic rule's numbers.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_quadratic_rule_on_a_turned_channel_keeps_its_recorded_error_everywhere() {
    check(&QUADRATIC.iter().collect::<Vec<_>>(), WallTreatment::InterpolatedQuadratic);
}

/// Each run of the order study lasts this many diffusion times H²/ν: the slowest transient,
/// sin(πy/H), decays as exp(−π² ν t / H²), to e^(−4π²) ≈ 7e-18 of its start, far below the
/// discretisation error.
const DIFFUSION_TIMES: u64 = 4;
/// Bouzidi, Firdaouss & Lallemand (2001) show second order in space for their quadratic
/// interpolation. The fitted order must reach this bar, half a unit under the theoretical 2, to
/// allow for the fit's scatter over four grids.
const ORDER_BAR: f64 = 1.5;

/// One channel of the order study: its angle, the widths across it in lattice rows for each of the
/// four grids, and the wall offset in cells. Rows of direction `a:b` lie 1/√(a²+b²) cells apart, so
/// a width of `n` rows is `n / √(a²+b²)` cells. Each grid holds a whole number of rows, and the
/// offset is fixed in cells, so every grid cuts its links at the same fractions q and the fit
/// measures the rule's order, not how each grid happens to cut the wall.
struct OrderCase {
    angle: (usize, usize),
    rows: [usize; 4],
    phase: f64,
}

// 1:1 at offset 0 is left out: its wall runs through the centre row of cells, and rounding puts
// that row on a different side of the wall on each grid, so the fit measures that scatter, not the
// rule's order. The wall a millionth of a cell either side of that row is here instead: its normal
// links are cut at q = 1e-6/√2 ≈ 7.1e-7.
const ORDER_CASES: [OrderCase; 8] = [
    OrderCase { angle: (1, 1), rows: [11, 17, 23, 34], phase: 0.1 },
    OrderCase { angle: (1, 1), rows: [11, 17, 23, 34], phase: 0.25 },
    OrderCase { angle: (1, 1), rows: [11, 17, 23, 34], phase: 0.5 },
    OrderCase { angle: (1, 1), rows: [11, 17, 23, 34], phase: 1e-6 },
    OrderCase { angle: (1, 1), rows: [11, 17, 23, 34], phase: -1e-6 },
    OrderCase { angle: (2, 1), rows: [18, 27, 36, 54], phase: 0.1 },
    OrderCase { angle: (2, 1), rows: [18, 27, 36, 54], phase: 0.25 },
    OrderCase { angle: (2, 1), rows: [18, 27, 36, 54], phase: 0.67 },
];

/// The velocity's relative L2 error on each channel of [`ORDER_CASES`], on each of its four grids,
/// for the quadratic rule and for the linear rule beside it. The quadratic rule's fitted order must
/// reach [`ORDER_BAR`] on every channel.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_quadratic_rule_on_a_turned_channel_converges_at_second_order() {
    let nu_lat = 0.05;
    let mut failures = Vec::new();
    for c in &ORDER_CASES {
        let (a, b) = c.angle;
        let widths: [f64; 4] = c.rows.map(|n| n as f64 / ((a * a + b * b) as f64).sqrt());
        let study = |walls: WallTreatment| -> (Vec<f64>, Vec<f64>) {
            widths
                .iter()
                .map(|&width| {
                    let steps = (DIFFUSION_TIMES as f64 * width * width / nu_lat).ceil() as u64;
                    let r = rotated_channel(c.angle, width, c.phase, walls, steps);
                    println!("{:?} offset {} on {width:.3} cells, {}: links {:?}", c.angle, c.phase, walls.name(), r.links);
                    (r.l2, r.flow_rate)
                })
                .unzip()
        };
        let (quadratic, flow_q) = study(WallTreatment::InterpolatedQuadratic);
        let (linear, flow_l) = study(WallTreatment::Interpolated);
        let (oq, ol) = (fitted_order(&widths, &quadratic), fitted_order(&widths, &linear));
        println!(
            "{:?} offset {} on {widths:.3?} cells: L2 quadratic {} (order {oq:.2}), linear {} (order {ol:.2}); flow rate % quadratic {}, linear {}",
            c.angle,
            c.phase,
            sci(&quadratic),
            sci(&linear),
            sci(&flow_q),
            sci(&flow_l)
        );
        if oq <= ORDER_BAR {
            failures.push(format!("{:?} offset {}: quadratic order {oq:.2} is not above {ORDER_BAR} (L2 {})", c.angle, c.phase, sci(&quadratic)));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// The quadratic form is Lagrange interpolation: it keeps a constant, reproduces a parabola in
/// the populations for q < ½ and is halfway bounce-back at q = ½.
#[test]
fn the_quadratic_form_keeps_a_constant_reproduces_a_parabola_and_is_halfway_bounce_back_at_one_half() {
    use fieldcraft_lbm::flow::bouzidi_quadratic;
    for q in [0.05, 0.2, 0.4, 0.49, 0.5, 0.51, 0.6, 0.75, 0.95, 1.0f64] {
        let k = bouzidi_quadratic(q, 3.0, 3.0, 3.0, 3.0, 3.0);
        assert!((k - 3.0).abs() < 1e-12, "q {q}: {k}");
        if q < 0.5 {
            // The characteristic that starts 1 - 2q cells behind the cell, read from nodes 0, 1, 2.
            let p = |x: f64| 0.3 + 0.5 * x - 0.2 * x * x;
            let v = bouzidi_quadratic(q, p(0.0), p(1.0), p(2.0), 0.0, 0.0);
            assert!((v - p(1.0 - 2.0 * q)).abs() < 1e-12, "q {q}: {v}");
        }
    }
    assert!((bouzidi_quadratic(0.5, 0.7, 0.1, 0.2, 0.3, 0.4) - 0.7).abs() < 1e-12);
    assert!((bouzidi_quadratic(0.5 - 1e-12, 0.7, 0.1, 0.2, 0.3, 0.4) - 0.7).abs() < 1e-9);
}

/// The mean speed of the stability runs. The exact centreline speed is 1.5 times it, 0.0225, which
/// is under the stream-stability limit of TRT at τ⁺ 0.505, 5.6·(τ⁺ − ½) ≈ 0.028 (issue #1063, #1076;
/// not on develop yet).
const STABILITY_MEAN: f64 = 0.015;

/// The relative margin on the discrete exact peak: the flow-rate bar of [`check`], 2e-3 percent. check()
/// has no velocity bar of its own, so its flow-rate bar is the one tolerance this channel already uses.
const PEAK_MARGIN: f64 = 2e-3 / 100.0;

/// The density drift per 10 000 steps of the stability runs, at offsets +1e-6 and −1e-6, as the code
/// gives them on this branch. Measured: the floor code on develop 655d92fe gives TRT −3.119429e-2 and
/// −2.075217e-2, regularized −1.759055e-4 and −2.278871e-3, so removing the floor moves the drift by
/// under 1 % (TRT −1e-6: 0.07 %; regularized −1e-6: 0.9 %), not materially. The drift values are measured;
/// their cause is inferred: TRT at τ⁺ 0.505 with Λ = 3/16 has τ⁻ = 1/2 + (3/16)/(τ⁺ − ½) ≈ 38, so the
/// odd mode relaxes at 1/τ⁻ ≈ 0.026 per step against 1/τ⁺ ≈ 2.0 for the shear, and the interpolated wall's
/// mass defect is left unrelaxed at that rate (the cause is inferred, not checked against a second run).
const TRT_STABILITY_DRIFT: [f64; 2] = [-3.119429e-2, -2.076600e-2];
/// The regularized collision's drift in the same runs (its values are small; the same origin is not
/// inferred for them).
const REGULARIZED_STABILITY_DRIFT: [f64; 2] = [-1.759055e-4, -2.298570e-3];

/// Links with a fraction q below 0.01: the 45° wall a millionth of a cell off a row of cell centres
/// cuts its normal links at q ≈ 7.1e-7. Run at the viscosity floor of the setup, [`MIN_TAU`] (τ⁺ 0.505,
/// the lowest laminar τ for every collision), with TRT and the regularized collision, for
/// [`DIFFUSION_TIMES`] diffusion times H²/ν of the 11-cell channel (about 290 000 steps), so that a
/// slowly growing mode would show. The bound is the exact centreline speed of the parabola, 1.5 times
/// the mean, sampled at the fluid cells' centres, with the margin [`PEAK_MARGIN`] for a sound run's small
/// overshoot. A growing mode fails on the speed and on the drift. The drift is checked against its known value.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn links_with_a_tiny_fraction_keep_a_channel_bounded_at_the_viscosity_floor() {
    let nu_lat = (MIN_TAU - 0.5) / 3.0;
    let width = 11.0;
    let steps = (DIFFUSION_TIMES as f64 * width * width / nu_lat).ceil() as u64;
    let mut runs = Vec::new();
    for (collision, known) in [(Collision::Trt, TRT_STABILITY_DRIFT), (Collision::Regularized, REGULARIZED_STABILITY_DRIFT)] {
        for (k, phase) in [1e-6, -1e-6].into_iter().enumerate() {
            let regime = Regime { collision, nu_lat, mean_speed: STABILITY_MEAN };
            let r = run_channel((1, 1), width, phase, WallTreatment::Interpolated, steps, regime);
            println!(
                "{collision:?}, offset {phase:+e}, {steps} steps: max speed {:.6e} (discrete peak {:.6e}), drift {:+.6e} (known {:+.6e}), links {:?}",
                r.max_speed, r.exact_peak, r.drift, known[k], r.links
            );
            runs.push((collision, phase, known[k], r));
        }
    }
    for (collision, phase, known, r) in &runs {
        assert!(r.max_speed.is_finite() && r.drift.is_finite(), "{collision:?} {phase}: {r:?}");
        let bound = r.exact_peak * (1.0 + PEAK_MARGIN);
        assert!(
            r.max_speed <= bound,
            "{collision:?} {phase}: max speed {} above the discrete peak {} with margin {bound}",
            r.max_speed,
            r.exact_peak
        );
        assert!((r.drift - known).abs() < 1e-5 + 0.01 * known.abs(), "{collision:?} {phase}: drift {} against {known}", r.drift);
    }
}
