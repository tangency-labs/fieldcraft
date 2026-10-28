//! The stationary toy: the numerical linear analysis of the streamwise odd-even (k_x = π) mode of
//! the wall rows under the wall-function bounces, Han et al.'s on-node eq. 14
//! (`spalding_bounce`) and our halfway variant (`spalding_bounce_halfway`).
//!
//! In channel runs at τ̄ ≈ 0.50015 the halfway form grew an x-checkerboard of u_x in its wall
//! rows (about 3e-5 per step) until the run broke. Here a thin periodic box, 16 × (rows) × 2 cells
//! in lattice units, both walls with the wall model, carries a steady x-uniform flow: the laminar
//! profile between the walls whose first node reads the friction velocity the body force needs,
//! `g · (fluid rows) = 2 u_τ²`, and whose first-node speed is Spalding's (κ 0.41, B 5.5) at that
//! u_τ. The subgrid model is off (the setup's constant set to 0 after it is built: the wall
//! treatment alone). A checkerboard of 1e-6 in u_x is added to the first fluid row at each wall;
//! the base is uniform in x, so the k_x = π part of the field evolves on its own, linearly. Its
//! energy `E = Σ_rows Σ_components |û|²`, `û = Σ_{x,z} (−1)^x u / (n_x n_z)`, is recorded each step;
//! the growth rate per step is the least-squares slope of `½ ln E` over the second half of the run
//! (the decaying transients gone).
//!
//! The on-node box is 2N + 1 cells across (the outermost layers are the boundary nodes, 2N
//! between the walls, the first fluid node one cell from the wall), the halfway box 2N (the wall
//! half a cell outside, the first cell half a cell from it), N = 4.

mod common;

use common::*;
use fieldcraft_lbm::{Flow, Limits, Setup};
use fieldcraft_model::{Collision, FaceCondition, HrrSigma, Procedure, Turbulence, WallModel};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};

/// Half the channel in cells, the box's length and depth.
const N: usize = 4;
const NX: usize = 16;
const NZ: usize = 2;

/// The largest lattice speed the base may reach (at the centreline).
const MAX_SPEED: f64 = 0.1;

/// Spalding's y⁺(u⁺) with Han et al.'s constants.
fn spalding_y_plus(up: f64) -> f64 {
    let k = 0.41 * up;
    up + (-0.41f64 * 5.5).exp() * (k.exp() - 1.0 - k - k * k / 2.0 - k * k * k / 6.0)
}

/// One toy run's outcome.
#[derive(Debug, Clone, Copy)]
struct Toy {
    /// The growth rate per step of the k_x = π mode's amplitude.
    rate: f64,
    /// The first node's u⁺ and y⁺ of the base.
    u_plus: f64,
    y_plus: f64,
    /// The step the run broke at, if it did, and the step the mode fell below round-off, if it
    /// did (the rate is measured before it).
    broke: Option<u64>,
    decayed: Option<u64>,
}

/// The toy for wall model `wall`, collision `collision`, relaxation time `tau`, the first node's
/// u⁺ aimed at `u_plus` (lowered until the centreline is at most [`MAX_SPEED`]), run `steps`.
fn toy(wall: WallModel, collision: Collision, tau: f64, u_plus: f64, steps: u64) -> Toy {
    let on_node = wall == WallModel::SpaldingBounce;
    let nu = (tau - 0.5) / 3.0;
    // Rows across, fluid rows, the first node's distance from its wall, and each fluid row's
    // distance from the lower wall.
    let (ny, rows, d) = if on_node { (2 * N + 1, 2 * N - 1, 1.0) } else { (2 * N, 2 * N, 0.5) };
    let y_of = |j: usize| if on_node { j as f64 } else { j as f64 + 0.5 };
    let h = N as f64;
    // The base for a first-node u⁺: u_τ from y⁺ = u_τ d / ν, the body force from the balance, the
    // laminar profile above the first node.
    let base = |up: f64| {
        let ut = spalding_y_plus(up) * nu / d;
        let g = 2.0 * ut * ut / rows as f64;
        let u1 = ut * up;
        (ut, g, move |y: f64| u1 + g / (2.0 * nu) * ((h - d).powi(2) - (h - y).powi(2)))
    };
    let mut up = u_plus;
    while base(up).2(h) > MAX_SPEED && up > 0.05 {
        up *= 0.98;
    }
    let (ut, g, profile) = base(up);

    let mut dom = domain([0.0; 3], [NX as f64, ny as f64, NZ as f64], 1.0);
    for f in [&mut dom.faces.x_min, &mut dom.faces.x_max, &mut dom.faces.z_min, &mut dom.faces.z_max] {
        *f = PERIODIC;
    }
    dom.acceleration = [g, 0.0, 0.0];
    let mut p = collide(steps as f64, 1.0, None, collision);
    if let Procedure::Fluid { turbulence, sigma, .. } = &mut p {
        *turbulence = Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: Some(wall) });
        if collision == Collision::Hrr {
            *sigma = Some(HrrSigma::Value(0.99));
        }
    }
    let mut setup = Setup::new(&model(1.0, nu, dom, p), None, &[], &Limits::default()).unwrap();
    assert!((setup.tau_plus - tau).abs() < 1e-12, "τ {}", setup.tau_plus);
    // The wall treatment alone: no subgrid model.
    setup.smagorinsky = 0.0;
    let mut flow = Flow::new(&setup, &[]).unwrap();
    assert!(!flow.turbulent());
    let grid = *flow.grid();
    // The first fluid row at each wall.
    let (first, last) = if on_node { (1, ny - 2) } else { (0, ny - 1) };
    flow.set_equilibrium(|cell| {
        let [i, j, _] = grid.coords(cell);
        let y = y_of(j).min(2.0 * h - y_of(j));
        let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
        let seed = if j == first || j == last { 1e-6 * sign } else { 0.0 };
        (1.0, [profile(y) + seed, 0.0, 0.0])
    });
    let energy = |flow: &Flow| {
        let mut e = 0.0;
        for j in first..=last {
            let mut a = [0.0; 3];
            for k in 0..NZ {
                for i in 0..NX {
                    let (_, u) = flow.moments(grid.index([i, j, k])).unwrap();
                    let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
                    for c in 0..3 {
                        a[c] += sign * u[c];
                    }
                }
            }
            e += a.iter().map(|x| (x / (NX * NZ) as f64).powi(2)).sum::<f64>();
        }
        e
    };
    let mut log_e = Vec::with_capacity(steps as usize);
    let mut broke = None;
    for k in 0..steps {
        if flow.step().is_err() {
            broke = Some(k + 1);
            break;
        }
        let e = energy(&flow);
        if !(e.is_finite() && e > 0.0) {
            broke = Some(k + 1);
            break;
        }
        log_e.push(0.5 * e.ln());
    }
    // The least-squares slope over the second half of what ran while the mode stood above
    // round-off (an amplitude of 1e-14; the base is about 1e-2, its round-off 1e-18).
    let floor = log_e.iter().position(|a| *a < 1e-14f64.ln()).unwrap_or(log_e.len());
    let decayed = (floor < log_e.len()).then_some(floor as u64);
    let tail = &log_e[floor / 2..floor];
    let n = tail.len() as f64;
    let mean_t = (n - 1.0) / 2.0;
    let mean_y = tail.iter().sum::<f64>() / n;
    let (mut sxy, mut sxx) = (0.0, 0.0);
    for (t, y) in tail.iter().enumerate() {
        sxy += (t as f64 - mean_t) * (y - mean_y);
        sxx += (t as f64 - mean_t).powi(2);
    }
    let rate = if tail.len() < 2 { f64::NAN } else { sxy / sxx };
    Toy { rate, u_plus: up, y_plus: ut * d / nu, broke, decayed }
}

fn name(wall: WallModel) -> &'static str {
    match wall {
        WallModel::SpaldingBounce => "on-node (eq. 14)",
        _ => "halfway",
    }
}

/// The trunk's representative case: Han et al.'s on-node eq. 14 on the regularised collision at
/// τ̄ 0.5015, first node at u⁺ ≈ 5: the k_x = π mode does not grow.
#[test]
fn the_on_node_bounce_damps_the_odd_even_mode_of_its_wall_rows() {
    let t = toy(WallModel::SpaldingBounce, Collision::Regularized, 0.5015, 5.0, 4000);
    println!("on-node, regularised, τ̄ 0.5015, u⁺ {:.2} (y⁺ {:.2}): growth {:+.3e} per step", t.u_plus, t.y_plus, t.rate);
    assert!(t.broke.is_none() && t.rate <= 0.0, "{t:?}");
}

/// The full sweep: both forms, the regularised, hybrid (σ 0.99), MRT and TRT collisions,
/// τ̄ 0.50015 to 0.8, first node at u⁺ ≈ 5 and ≈ 11 (lowered where the centreline would pass
/// lattice speed 0.1); 20 000 steps each. Prints the table of growth rates per step.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_odd_even_mode_of_the_wall_rows_across_forms_collisions_and_tau() {
    let steps = 20_000;
    println!("| form | collision | τ̄ | u⁺ aimed | u⁺ | y⁺ | growth per step |");
    println!("|---|---|---|---|---|---|---|");
    for wall in [WallModel::SpaldingBounce, WallModel::SpaldingBounceHalfway] {
        for collision in [Collision::Regularized, Collision::Hrr, Collision::Mrt, Collision::Trt] {
            for tau in [0.50015, 0.5015, 0.505, 0.55, 0.8] {
                for aim in [5.0, 11.0] {
                    let t = toy(wall, collision, tau, aim, steps);
                    let growth = match (t.broke, t.decayed) {
                        (Some(k), _) => format!("{:+.3e} (broke at step {k})", t.rate),
                        (None, Some(k)) => format!("{:+.3e} (below 1e-14 at step {k})", t.rate),
                        (None, None) => format!("{:+.3e}", t.rate),
                    };
                    println!("| {} | {} | {tau} | {aim} | {:.2} | {:.2} | {growth} |", name(wall), collision.name(), t.u_plus, t.y_plus);
                }
            }
        }
    }
}
