//! The half-way wall (`spalding_bounce_surface`): Han, Ooka and Kikumoto's on-node bounce with the
//! wall half way between the first open cell and the node (the first solid layer), the law at
//! y = ½ cell.
//!
//! Laminar Poiseuille flow between two obstacle slabs, the first gate: the flow rate against the
//! analytic `g H³ / (12 ν)` for a width H equal to the open cells. Flow goes as the width cubed, so
//! it tests the effective wall position directly; with the law at the node (y = 1) the wall is
//! half a cell too far out, and the channel over-flows.

mod common;

use common::*;
use fieldcraft_lbm::Flow;
use fieldcraft_model::{Collision, FaceCondition, Procedure, Turbulence, WallModel};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};

/// Flow rate (physical, per unit span) of a Poiseuille flow of `open` cells between slabs two cells
/// thick, divided by the analytic one for the width of `open` cells; `distance` overrides the law's
/// y (cells) when given.
fn flow_rate_ratio(open: usize, model_kind: Option<WallModel>, distance: Option<f64>, steps: u64) -> f64 {
    let (nu, height) = (1.0, 1.0);
    let dx = height / open as f64;
    // Lattice viscosity 0.05 and a mean speed of 0.02 lattice units.
    let dt = 0.05 * dx * dx / nu;
    let u_mean = 0.02 * dx / dt;
    let g = 12.0 * nu * u_mean / (height * height);
    let span = 4.0 * dx;
    let mut d = domain([0.0, -2.0 * dx, 0.0], [span, height + 2.0 * dx, span], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    let obstacles = vec![
        obstacle("Floor", cuboid([-1.0, -2.0 * dx, -1.0], [span + 1.0, 0.0, span + 1.0]), &CUBOID_FACES),
        obstacle("Roof", cuboid([-1.0, height, -1.0], [span + 1.0, height + 2.0 * dx, span + 1.0]), &CUBOID_FACES),
    ];
    d.obstacles = obstacles.iter().map(|o| o.body.clone()).collect();
    d.acceleration = [g, 0.0, 0.0];
    let mut p = collide(steps as f64 * dt, dt, None, Collision::Trt);
    if let Procedure::Fluid { turbulence: t, .. } = &mut p {
        *t = Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: model_kind });
    }
    let (setup, mut flow) = prepare(&model(1.0, nu, d, p), &obstacles);
    if let Some(y) = distance {
        flow.set_wall_distance(y).unwrap();
    }
    for _ in 0..steps {
        flow.step().unwrap();
    }
    let _ = setup;
    flux(&flow, open) / (0.02 * open as f64)
}

/// Σ u over the open rows of the column (1, ·, 1) in lattice units: the flow rate in cells, to be
/// compared with the analytic mean speed (0.02) times the width in cells.
fn flux(flow: &Flow, open: usize) -> f64 {
    let grid = *flow.grid();
    let mut sum = 0.0;
    for j in 2..2 + open {
        let cell = grid.index([1, j, 1]);
        sum += flow.moments(cell).map_or(0.0, |m| m.1[0]);
    }
    sum
}

/// The laminar flow rates, measured: the law reads the second open cell (y = 1½), where a
/// parabola's slope is lower than at the wall by 1 − 1½/H, so a laminar passage over-flows (88 %
/// at 4 open cells, 22 % at 8, 5 % at 16); plain half-way bounce-back is 3 %, 0.6 % and 0.1 % over.
/// The law is for turbulent walls, whose second cell lies in the log layer. Pinned so a change that
/// moves them is seen, and with the law read at y = 2 (the node distance 1 plus the second cell's 1,
/// as on Han's wall) the channel over-flows far more.
#[test]
fn poiseuille_flow_rates_are_pinned() {
    let surface = Some(WallModel::SpaldingBounceSurface);
    let near = |got: f64, want: f64| (got - want).abs() < 2e-3;
    let mut failures = Vec::new();
    for (open, plain, half, far) in [(4usize, 1.0271, 1.8758, 2.6243), (8, 1.0056, 1.2187, 1.5935), (16, 1.0008, 1.0542, 1.2417)] {
        let steps = 40 * (open * open) as u64;
        let (p, h, o) =
            (flow_rate_ratio(open, None, None, steps), flow_rate_ratio(open, surface, None, steps), flow_rate_ratio(open, surface, Some(1.0), steps));
        println!("open {open}: plain {p:.4} y=1½ {h:.4} y=2 {o:.4}");
        if !(near(p, plain) && near(h, half) && near(o, far)) {
            failures.push(format!("open {open}: plain {p:.4} (want {plain}), y=1½ {h:.4} (want {half}), y=2 {o:.4} (want {far})"));
        }
        if o <= h {
            failures.push(format!("open {open}: the law read further out must over-flow more"));
        }
    }
    assert!(failures.is_empty(), "{failures:?}");
}

/// The wall layer: two nodes per slab side, none on a box face, and the node count in the report.
#[test]
fn the_slab_channel_has_one_node_row_per_wall() {
    let (report, nodes) = slab_channel(8);
    assert_eq!(report.nodes, 2 * 16, "a node row of 4 x 4 cells on each side of the channel");
    assert_eq!(report.plain_faces, 0);
    assert_eq!(nodes, 32);
    // 8 open cells between the walls is a narrow passage (fewer than 12); 16 is not.
    assert_eq!(report.narrow_nodes, report.nodes);
    assert_eq!(slab_channel(16).0.narrow_nodes, 0);
}

/// The report and node count of a channel `open` cells wide between two slabs two cells thick.
fn slab_channel(open: usize) -> (fieldcraft_lbm::wall_layer::WallLayerReport, usize) {
    let dx = 1.0 / open as f64;
    let mut d = domain([0.0, -2.0 * dx, 0.0], [4.0 * dx, 1.0 + 2.0 * dx, 4.0 * dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    let span = 4.0 * dx;
    let obstacles = vec![
        obstacle("Floor", cuboid([-1.0, -2.0 * dx, -1.0], [span + 1.0, 0.0, span + 1.0]), &CUBOID_FACES),
        obstacle("Roof", cuboid([-1.0, 1.0, -1.0], [span + 1.0, 1.0 + 2.0 * dx, span + 1.0]), &CUBOID_FACES),
    ];
    d.obstacles = obstacles.iter().map(|o| o.body.clone()).collect();
    let mut p = collide(10.0 * dx, dx, None, Collision::Trt);
    if let Procedure::Fluid { turbulence: t, .. } = &mut p {
        *t = Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: Some(WallModel::SpaldingBounceSurface) });
    }
    let (_, flow) = prepare(&model(1.0, 1.0, d, p), &obstacles);
    assert_eq!(flow.wall_nodes().map(|n| n.distance()), Some(0.5));
    (flow.wall_layer_report(), flow.wall_nodes().map_or(0, |n| n.nodes().len()))
}

/// A no-slip face of the box with open cells against it is refused: the half-way wall is for
/// obstacle walls (a slab in the box is the way to a wall).
#[test]
fn a_box_wall_face_is_refused() {
    let dx = 1.0 / 8.0;
    let mut d = domain([0.0; 3], [4.0 * dx, 1.0, 4.0 * dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    let mut p = collide(10.0 * dx, dx, None, Collision::Trt);
    if let Procedure::Fluid { turbulence: t, .. } = &mut p {
        *t = Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: Some(WallModel::SpaldingBounceSurface) });
    }
    let m = model(1.0, 1.0, d, p);
    let result = fieldcraft_lbm::Setup::new(&m, None, &[], &fieldcraft_lbm::Limits::default()).and_then(|s| Flow::new(&s, &[]));
    let message = result.err().map(|e| e.to_string()).unwrap_or_default();
    assert!(message.contains("not supported yet") && message.contains("spalding_bounce_surface"), "{message}");
}

/// A wall thinner than the spacing has no centre inside it: the leak is measured and, above 2 % of
/// the surface, refused; a thick enough slab reports none.
#[test]
fn a_wall_with_no_centre_inside_is_a_leak() {
    let dx = 1.0 / 8.0;
    let build = |thickness: f64| {
        let mut d = domain([0.0; 3], [4.0 * dx, 1.0, 4.0 * dx], dx);
        for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
            *f = PERIODIC;
        }
        // A plate across the channel; its centre line sits mid cell, so a thin one holds no centre.
        let lo = 0.5 - 0.5 * thickness;
        let obstacles = vec![obstacle("Plate", cuboid([-1.0, lo, -1.0], [1.0, lo + thickness, 1.0]), &CUBOID_FACES)];
        d.obstacles = obstacles.iter().map(|o| o.body.clone()).collect();
        let mut p = collide(10.0 * dx, dx, None, Collision::Trt);
        if let Procedure::Fluid { turbulence: t, .. } = &mut p {
            *t = Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: Some(WallModel::SpaldingBounceSurface) });
        }
        let m = model(1.0, 1.0, d, p);
        fieldcraft_lbm::Setup::new(&m, None, &obstacles, &fieldcraft_lbm::Limits::default()).and_then(|s| Flow::new(&s, &obstacles))
    };
    // 0.3 dx thick, centred on a cell face.
    let thin = build(0.3 * dx);
    let message = thin.err().map(|e| e.to_string()).unwrap_or_default();
    assert!(message.contains("no cell centre inside it"), "a thin plate on a cell face leaks: {message}");
    let thick = build(4.0 * dx).expect("a slab four cells thick is sealed");
    assert_eq!(thick.wall_layer_report().unsealed_share, 0.0);
    assert!(thick.wall_layer_report().nodes > 0);
}

/// A 2-cell checkerboard `(-1)^(i+k)` seeded in the first open row on each side of a slab channel,
/// `steps` steps later: the amplitudes (velocity along the wall, density) of that mode in the row,
/// divided by their initial values. `tau_minus_half` is τ - ½ of the symmetric part.
fn checkerboard(collision: Collision, tau_minus_half: f64, model_kind: Option<WallModel>, steps: u64) -> Result<(f64, f64), String> {
    let (nu, open) = (1.0, 8usize);
    let dx = 1.0 / open as f64;
    let nu_lattice = tau_minus_half / 3.0;
    let dt = nu_lattice * dx * dx / nu;
    let span = 4.0 * dx;
    let mut d = domain([0.0, -2.0 * dx, 0.0], [span, 1.0 + 2.0 * dx, span], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    let obstacles = vec![
        obstacle("Floor", cuboid([-1.0, -2.0 * dx, -1.0], [span + 1.0, 0.0, span + 1.0]), &CUBOID_FACES),
        obstacle("Roof", cuboid([-1.0, 1.0, -1.0], [span + 1.0, 1.0 + 2.0 * dx, span + 1.0]), &CUBOID_FACES),
    ];
    d.obstacles = obstacles.iter().map(|o| o.body.clone()).collect();
    let mut p = collide(steps as f64 * dt, dt, None, collision);
    if let Procedure::Fluid { turbulence: t, sigma, .. } = &mut p {
        *t = Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: model_kind });
        if collision == Collision::Hrr {
            *sigma = Some(fieldcraft_model::HrrSigma::Value(0.99));
        }
    }
    let m = model(1.0, nu, d, p);
    let setup = fieldcraft_lbm::Setup::new(&m, None, &obstacles, &fieldcraft_lbm::Limits::default()).map_err(|e| e.to_string())?;
    let mut flow = Flow::new(&setup, &obstacles).map_err(|e| e.to_string())?;
    let grid = *flow.grid();
    let (mean, a_u, a_rho) = (0.02, 0.004, 0.0004);
    let first = |j: usize| j == 2 || j == open + 1;
    flow.set_equilibrium(|cell| {
        let [i, j, k] = grid.coords(cell);
        let eta = (j as f64 - 1.5) / open as f64;
        let u = 6.0 * mean * eta * (1.0 - eta);
        let sign = if (i + k) % 2 == 0 { 1.0 } else { -1.0 };
        if first(j) { (1.0 + sign * a_rho, [u + sign * a_u, 0.0, 0.0]) } else { (1.0, [u, 0.0, 0.0]) }
    });
    let project = |flow: &Flow| {
        let (mut su, mut sr, mut count) = (0.0, 0.0, 0.0);
        for j in [2, open + 1] {
            for k in 0..4 {
                for i in 0..4 {
                    let sign = if (i + k) % 2 == 0 { 1.0 } else { -1.0 };
                    if let Some((rho, u)) = flow.moments(grid.index([i, j, k])) {
                        su += sign * u[0];
                        sr += sign * (rho - 1.0);
                    }
                    count += 1.0;
                }
            }
        }
        ((su / count).abs(), (sr / count).abs())
    };
    let (u0, r0) = project(&flow);
    for _ in 0..steps {
        flow.step().map_err(|e| e.to_string())?;
    }
    let (u1, r1) = project(&flow);
    Ok((u1 / u0, r1 / r0))
}

/// The checkerboard gate: a 2-cell mode seeded in the first open row decays below 1 % in 200 steps at
/// both ends of the viscosity range (the largest lattice viscosity, τ = 1, and τ - ½ = 2.4e-4) on every
/// collision that runs there. The two-relaxation-time collision does not damp the shortest waves as
/// τ → ½ whatever the walls do (plain half-way walls keep 38 % of the velocity mode, the law 44 %),
/// so there the law is held to plain walls' decay.
#[test]
fn the_checkerboard_mode_decays_with_the_law_at_both_ends_of_the_viscosity_range() {
    let law = Some(WallModel::SpaldingBounceSurface);
    for collision in [Collision::Trt, Collision::Regularized, Collision::Mrt, Collision::Hrr] {
        for tau in [0.5, 0.15, 2.4e-4] {
            let (u, r) = checkerboard(collision, tau, law, 200).expect("runs");
            if collision == Collision::Trt && tau < 1e-3 {
                let (plain_u, plain_r) = checkerboard(collision, tau, None, 200).expect("runs");
                assert!(
                    u < 1.25 * plain_u && r < 1.25 * plain_r,
                    "TRT at τ - ½ {tau}: law ({u:.3}, {r:.3}) against plain walls ({plain_u:.3}, {plain_r:.3})"
                );
            } else {
                assert!(u < 0.01 && r < 0.01, "{collision:?} at τ - ½ {tau}: the mode keeps u {u:.3e}, rho {r:.3e} of its amplitude");
            }
        }
    }
}

/// In a passage two cells wide each wall node's second open cell is the opposite wall's first, so no
/// node samples a second cell: every open cell's ν̃ is set once, from the law at the first cell's
/// distance (½ cell), and the same on every run.
#[test]
fn a_two_cell_passage_sets_each_nu_tilde_once() {
    let open = 2usize;
    let dx = 1.0 / open as f64;
    let span = 4.0 * dx;
    let mut d = domain([0.0, -2.0 * dx, 0.0], [span, 1.0 + 2.0 * dx, span], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    let obstacles = vec![
        obstacle("Floor", cuboid([-1.0, -2.0 * dx, -1.0], [span + 1.0, 0.0, span + 1.0]), &CUBOID_FACES),
        obstacle("Roof", cuboid([-1.0, 1.0, -1.0], [span + 1.0, 1.0 + 2.0 * dx, span + 1.0]), &CUBOID_FACES),
    ];
    d.obstacles = obstacles.iter().map(|o| o.body.clone()).collect();
    // Lattice viscosity 0.05 and a mean speed of 0.02 lattice units, as the flow-rate gate.
    let dt = 0.05 * dx * dx;
    d.acceleration = [12.0 * 0.02 * dx / dt, 0.0, 0.0];
    let mut p = collide(100.0 * dt, dt, None, Collision::Mrt);
    if let Procedure::Fluid { turbulence: t, .. } = &mut p {
        *t = Some(Turbulence::SpalartAllmaras { prandtl: 0.9, wall_model: Some(WallModel::SpalartAllmarasBounceSurface), inflow_ratio: 3.0 });
    }
    let (_, mut flow) = prepare(&model(1.0, 1.0, d, p), &obstacles);
    for _ in 0..100 {
        flow.step().unwrap();
    }
    let (law, samples) = flow.sa_wall_function();
    assert!(law.is_some() && !samples.is_empty());
    let nu_tilde = flow.spalart_allmaras().map(|f| f.nu_tilde().to_vec()).unwrap_or_default();
    for &(cell, y, u_tau) in samples {
        assert!(u_tau > 0.0, "cell {cell}: u_τ {u_tau}");
        assert_eq!(y, 0.5);
        let want = fieldcraft_lbm::sa::KAPPA * u_tau * y;
        let got = nu_tilde.get(cell).copied().unwrap_or(f64::NAN);
        assert!((got - want).abs() <= 1e-12 * want.max(1e-30), "cell {cell}: ν̃ {got} against κ u_τ y {want}");
    }
}
