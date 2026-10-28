//! Han, Ooka and Kikumoto's on-node wall-function bounce (`wall_model: "spalding_bounce"`, their
//! eq. 14; `crates/lbm/src/bounce.rs`): the assignments of one node, its mass, the wall force in a
//! coarse channel against the stresses, Spalding's law, the grid relation, the refusals and the
//! determinism.

mod common;

use common::*;
use fieldcraft_lbm::bounce::on_node;
use fieldcraft_lbm::lattice::{C, Q};
use fieldcraft_lbm::{Flow, Limits, Setup};
use fieldcraft_model::{Collision, FaceCondition, HrrSigma, Procedure, Turbulence, WallModel};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};

/// Our direction of velocity `v`.
fn dir(v: [i32; 3]) -> usize {
    C.iter().position(|c| *c == v).unwrap()
}

/// Populations that are all different: 0.05 + 0.001 q + small irrational parts.
fn populations() -> [f64; Q] {
    std::array::from_fn(|q| 0.05 + 0.001 * q as f64 + 1e-4 * ((q * q) as f64).sqrt().sin())
}

/// Eq. 14 written out in the paper's numbering (their Table A-1), mapped to ours by velocity
/// vector, for the wall at y = 0 with the fluid above: every assignment, exactly. The
/// populations from outside the box (their f3, f7, f11, f14, f17 before the collision) are
/// replaced by the images of those from the fluid, so the node's 19 populations keep their sum.
#[test]
fn eq_14_assigns_every_population_of_a_node_on_the_floor() {
    // The paper's D3Q19 (Table A-1).
    let paper: [[i32; 3]; Q] = [
        [0, 0, 0],
        [1, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
        [0, 0, 1],
        [0, 0, -1],
        [1, 1, 0],
        [-1, -1, 0],
        [1, 0, 1],
        [-1, 0, -1],
        [0, 1, 1],
        [0, -1, -1],
        [1, -1, 0],
        [-1, 1, 0],
        [1, 0, -1],
        [-1, 0, 1],
        [0, 1, -1],
        [0, -1, 1],
    ];
    // Their eq. 10b–d: the moments the table must give.
    let (mut jx, mut jy, mut jz) = ([0i32; Q], [0i32; Q], [0i32; Q]);
    for (q, c) in paper.iter().enumerate() {
        (jx[q], jy[q], jz[q]) = (c[0], c[1], c[2]);
    }
    assert_eq!(jx, [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0]);
    assert_eq!(jy, [0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1]);
    assert_eq!(jz, [0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1]);
    let ours: [usize; Q] = paper.map(dir);
    assert_eq!(ours, [0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 15, 16, 9, 10, 13, 14, 17, 18]);

    let f = populations();
    // The stress on the wall (momentum to it per step), along +x and −z: the paper's
    // τ_w,x = −s_x, τ_w,z = −s_z.
    let s = [3e-4, 0.0, -1e-4];
    let (tx, tz) = (-s[0], -s[2]);
    let mut out = f;
    let taken = on_node(&mut out, 1, 1, s);
    // In the paper's numbering: P(a) is our population of their index a.
    let p = |g: &[f64; Q], a: usize| g[ours[a]];
    // Before the collision, the ghost from outside: their 3, 7, 11, 14, 17 are the images of
    // 4, 13, 18, 8, 12.
    let ghost = |a: usize| match a {
        3 => p(&f, 4),
        7 => p(&f, 13),
        11 => p(&f, 18),
        14 => p(&f, 8),
        17 => p(&f, 12),
        _ => p(&f, a),
    };
    // Eq. 14a–b.
    let f7 = p(&f, 13) + 0.5 * tx;
    let f14 = p(&f, 8) - 0.5 * tx;
    let f11 = p(&f, 18) + 0.5 * tz;
    let f17 = p(&f, 12) - 0.5 * tz;
    assert_eq!(p(&out, 7), f7);
    assert_eq!(p(&out, 14), f14);
    assert_eq!(p(&out, 11), f11);
    assert_eq!(p(&out, 17), f17);
    // Eq. 14c, and the rest population kept.
    let want = [
        (0, ghost(0)),
        (1, ghost(2)),
        (2, ghost(1)),
        (3, ghost(4)),
        (4, ghost(3)),
        (5, ghost(6)),
        (6, ghost(5)),
        (8, f7),
        (9, ghost(10)),
        (10, ghost(9)),
        (12, f11),
        (13, f14),
        (15, ghost(16)),
        (16, ghost(15)),
        (18, f17),
    ];
    for (a, v) in want {
        assert_eq!(p(&out, a), v, "the paper's f{a}*");
    }
    // Mass: the 19 populations with the ghost, before and after.
    let before: f64 = (0..Q).map(ghost).sum();
    let after: f64 = out.iter().sum();
    assert!((after - before).abs() < 1e-15, "{after} vs {before}");
    // The momentum the node takes from the fluid: the stress along the wall, the bounce-back's
    // 2 Σ f of the populations arriving from the fluid along the normal.
    let down: f64 = [4, 8, 12, 13, 18].map(|a| p(&f, a)).iter().sum();
    assert!((taken[0] - s[0]).abs() < 1e-16 && (taken[2] - s[2]).abs() < 1e-16, "{taken:?}");
    assert!((taken[1] + 2.0 * down).abs() < 1e-15, "{taken:?} vs {down}");
}

/// On every one of the six faces, eq. 14 is the floor's turned: the populations into the fluid
/// carry what arrived from it (mass), less the stress along the wall (momentum), whatever the
/// stress; the node's 19 populations keep their sum with the ghost.
#[test]
fn eq_14_conserves_mass_and_carries_the_stress_on_every_face() {
    let f = populations();
    for axis in 0..3 {
        for inward in [1, -1] {
            let mut s = [2e-4, -3e-4, 5e-5];
            s[axis] = 0.0;
            let mut out = f;
            let taken = on_node(&mut out, axis, inward, s);
            let cn = |q: usize| C[q][axis] * inward;
            let into: f64 = (0..Q).filter(|q| cn(*q) > 0).map(|q| out[q]).sum();
            let from: f64 = (0..Q).filter(|q| cn(*q) < 0).map(|q| f[q]).sum();
            assert!((into - from).abs() < 1e-15, "axis {axis} {inward}: {into} vs {from}");
            let parallel_in: f64 = (0..Q).filter(|q| cn(*q) == 0).map(|q| f[q]).sum();
            let parallel_out: f64 = (0..Q).filter(|q| cn(*q) == 0).map(|q| out[q]).sum();
            assert!((parallel_in - parallel_out).abs() < 1e-15);
            // Into the wall: the values sent into the fluid, reversed.
            for q in (1..Q).filter(|q| cn(*q) < 0) {
                assert_eq!(out[q], out[fieldcraft_lbm::lattice::opposite(q)]);
            }
            let ghosted = 2.0 * from + parallel_in;
            let total: f64 = out.iter().sum();
            assert!((total - ghosted).abs() < 1e-15, "{total} vs {ghosted}");
            for b in (0..3).filter(|b| *b != axis) {
                assert!((taken[b] - s[b]).abs() < 1e-16, "axis {axis} {inward}: {taken:?} vs {s:?}");
            }
            assert!((taken[axis] + 2.0 * f64::from(inward) * from).abs() < 1e-15);
        }
    }
}

/// A turbulent step with Smagorinsky's model and the on-node bounce.
fn bounced(end_time: f64, time_step: f64, collision: Collision) -> Procedure {
    let mut p = collide(end_time, time_step, None, collision);
    if let Procedure::Fluid { turbulence, sigma, .. } = &mut p {
        *turbulence = Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: Some(WallModel::SpaldingBounce) });
        if collision == Collision::Hrr {
            *sigma = Some(HrrSigma::Value(0.99));
        }
    }
    p
}

/// A coarse channel (6 cells per half-height, Re_τ 180, periodic in x and z) whose walls are the
/// on-node bounce: the box is 2δ + Δx high (13 cells), the walls through the centres of its
/// outermost cells, started from a perturbed blunt profile, driven by a body force.
fn perturbed_on_node_channel(collision: Collision, steps: f64) -> (Setup, Flow) {
    use std::f64::consts::PI;
    let (n, re_tau, ub) = (6usize, 180.0, 15.7);
    let dx = 1.0 / n as f64;
    let (lx, lz) = ((2.0 * PI * n as f64).round() * dx, (PI * n as f64).round() * dx);
    let mut d = domain([0.0; 3], [lx, 2.0 + dx, lz], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.acceleration = [1.0, 0.0, 0.0];
    let dt = 0.1 * dx / ub;
    let (setup, mut flow) = prepare(&model(1.0, 1.0 / re_tau, d, bounced(steps * dt, dt, collision)), &[]);
    let grid = *flow.grid();
    assert_eq!(grid.cells[1], 2 * n + 1);
    let to_lattice = 1.0 / setup.velocity_scale();
    let (kx, kz) = (2.0 * PI / lx, 2.0 * PI / lz);
    flow.set_equilibrium(|cell| {
        let [x, y, z] = grid.centre(cell);
        // The walls at y = Δx/2 and 2 + Δx/2.
        let eta = y - 1.0 - 0.5 * dx;
        let wall = (1.0 - eta * eta).max(0.0);
        let u = 1.15 * ub * (1.0 - eta.abs().min(1.0).powi(8)) + 2.0 * (2.0 * kz * z).sin() * (kx * x).cos() * wall;
        let v = 1.5 * wall * (2.0 * kz * z).cos() * (3.0 * kx * x).sin();
        let w = 1.5 * wall * eta * (kx * x).sin() * (kz * z).cos();
        (1.0, [u * to_lattice, v * to_lattice, w * to_lattice])
    });
    (setup, flow)
}

/// The fluid's momentum and mass: its cells' populations, and the populations that have
/// arrived at the boundary nodes from the fluid and not yet gone back (eq. 14 sends them back in
/// the next step).
fn fluid_momentum_mass(flow: &Flow) -> ([f64; 3], f64) {
    let (mut p, mut m) = ([0.0; 3], 0.0);
    let mut add = |cell: usize, q: usize| {
        let f = flow.population(q, cell);
        m += f;
        for a in 0..3 {
            p[a] += f64::from(C[q][a]) * f;
        }
    };
    for cell in (0..flow.grid().len()).filter(|c| flow.is_fluid(*c)) {
        for q in 0..Q {
            add(cell, q);
        }
    }
    for w in flow.wall_nodes().unwrap().nodes() {
        for q in (0..Q).filter(|q| C[*q][usize::from(w.axis)] * i32::from(w.inward) < 0) {
            add(w.cell as usize, q);
        }
    }
    (p, m)
}

/// Spalding's y⁺(u⁺) with Han et al.'s constants.
fn spalding_y_plus(up: f64) -> f64 {
    let k = 0.41 * up;
    up + (-0.41f64 * 5.5).exp() * (k.exp() - 1.0 - k - k * k / 2.0 - k * k * k / 6.0)
}

/// In the perturbed coarse channel, step by step from the start, on every collision: the walls
/// take exactly the stresses the nodes sampled along them (x and z, to round-off: the stress is
/// the whole of the tangential momentum the nodes take), the fluid's momentum budget
/// P(t + 1) − P(t) = a M − F closes and its mass is conserved (P and M with the populations in
/// flight at the nodes), and each stress is ρ u_τ² along the interior node's velocity parallel to
/// the wall, u_τ the root of Spalding's law (κ 0.41, B 5.5) at one cell from the wall.
#[test]
fn the_on_node_bounce_gives_the_walls_exactly_the_law_s_stress() {
    for collision in [Collision::Trt, Collision::Regularized, Collision::Mrt, Collision::Hrr] {
        let (setup, mut flow) = perturbed_on_node_channel(collision, 40.0);
        assert!(flow.on_node_bounce() && !flow.wall_bounce() && !flow.wall_model());
        let a = flow.acceleration();
        let nu = (setup.tau_plus - 0.5) / 3.0;
        let (mut p, m0) = fluid_momentum_mass(&flow);
        let (mut budget, mut carried, mut law, mut mass, mut shear) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
        for _ in 0..setup.steps {
            let before: Vec<Option<(f64, [f64; 3])>> = (0..flow.grid().len()).map(|c| flow.moments(c)).collect();
            let bulk = flow.mass();
            flow.step().unwrap();
            let (next, m) = fluid_momentum_mass(&flow);
            let nodes = flow.wall_nodes().unwrap().nodes();
            let wall: [f64; 3] = [0, 1, 2].map(|k| flow.forces().iter().map(|f| f[k]).sum());
            let stress: [f64; 3] = [0, 1, 2].map(|k| nodes.iter().map(|w| w.stress[k]).sum());
            for k in 0..3 {
                budget = budget.max((next[k] - p[k] - (a[k] * bulk - wall[k])).abs() / next[0].abs());
            }
            for k in [0, 2] {
                carried = carried.max((wall[k] - stress[k]).abs() / stress[0].abs());
            }
            for w in nodes {
                let (rho, u) = before[w.interior as usize].unwrap();
                assert!(before[w.cell as usize].is_none(), "a boundary node is not a fluid cell");
                let speed = (u[0] * u[0] + u[2] * u[2]).sqrt();
                let s = (w.stress[0] * w.stress[0] + w.stress[2] * w.stress[2]).sqrt();
                let ut = (s / rho).sqrt();
                law = law.max((spalding_y_plus(speed / ut) / (ut / nu) - 1.0).abs());
                // Along the velocity parallel to the wall.
                shear = shear.max((w.stress[0] * u[2] - w.stress[2] * u[0]).abs() / (s * speed));
                assert_eq!(w.stress[1], 0.0);
            }
            mass = mass.max((m / m0 - 1.0).abs());
            p = next;
        }
        println!(
            "{} on-node bounce, {} steps: momentum budget {budget:.1e} of the momentum, carried vs sampled stress {carried:.1e}, Spalding's law {law:.1e}, direction {shear:.1e}, mass {mass:.1e}",
            collision.name(),
            setup.steps
        );
        assert!(budget < 1e-12 && carried < 1e-12 && law < 1e-10 && shear < 1e-12 && mass < 1e-13, "{budget} {carried} {law} {shear} {mass}");
    }
}

/// The grid relation: the outermost layers are the boundary nodes, the wall through their
/// centres, their interior nodes one cell further; they are no fluid cells (left out of the
/// cell count, the mass and the fields), and before the first step they hold the equilibrium.
#[test]
fn the_boundary_nodes_are_the_outermost_layers_and_not_fluid() {
    use fieldcraft_lbm::backend::Backend;
    let (_, mut flow) = perturbed_on_node_channel(Collision::Regularized, 1.0);
    let grid = *flow.grid();
    let [nx, ny, nz] = grid.cells;
    let nodes = flow.wall_nodes().unwrap().nodes().to_vec();
    assert_eq!(nodes.len(), 2 * nx * nz);
    assert_eq!(flow.fluid_cells(), nx * (ny - 2) * nz);
    for w in &nodes {
        let [_, j, _] = grid.coords(w.cell as usize);
        let [_, i, _] = grid.coords(w.interior as usize);
        assert!((j == 0 && i == 1 && w.inward == 1 && w.group == 2) || (j == ny - 1 && i == ny - 2 && w.inward == -1 && w.group == 3), "{w:?}");
        assert!(!flow.is_fluid(w.cell as usize) && flow.is_fluid(w.interior as usize));
    }
    assert!(nodes.windows(2).all(|p| p[0].cell < p[1].cell), "in cell order");
    flow.step().unwrap();
    let fields = flow.fields().unwrap();
    for w in &nodes {
        assert_eq!((fields.density[w.cell as usize], fields.velocity[w.cell as usize]), (0.0, [0.0; 3]));
    }
}

/// What the on-node bounce does not do yet is refused, naming it: a moving wall, an obstacle, a
/// wall that meets a face that is not periodic, heat, and a box too thin for an interior node.
#[test]
fn the_on_node_bounce_refuses_what_it_does_not_support() {
    let dx = 0.125;
    let channel = || {
        let mut d = domain([0.0; 3], [4.0 * dx, 1.0, 2.0 * dx], dx);
        for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
            *f = PERIODIC;
        }
        d
    };
    let refused = |d: fieldcraft_model::FluidDomain, obstacles: &[fieldcraft_lbm::Obstacle]| {
        let m = model(1.0, 0.01, d, bounced(10.0 * dx, dx, Collision::Mrt));
        match Setup::new(&m, None, obstacles, &Limits::default()) {
            Ok(setup) => Flow::new(&setup, obstacles).unwrap_err().to_string(),
            Err(e) => e.to_string(),
        }
    };
    let mut d = channel();
    d.faces.y_max = FaceCondition::MovingWall { velocity: [0.1, 0.0, 0.0] };
    let err = refused(d, &[]);
    assert!(err.contains("not supported yet") && err.contains("spalding_bounce") && err.contains("moving") && err.contains("y_max"), "{err}");
    let mut d = channel();
    d.faces.x_min = FaceCondition::Wall {};
    d.faces.x_max = FaceCondition::Wall {};
    let err = refused(d, &[]);
    assert!(err.contains("not supported yet") && err.contains("x_min") && err.contains("periodic"), "{err}");
    let mut d = channel();
    d.obstacles = vec!["Block".into()];
    let block = obstacle("Block", cuboid([0.1, 0.4, 0.0], [0.3, 0.6, 0.25]), &["side"]);
    let err = refused(d, &[block]);
    assert!(err.contains("not supported yet") && err.contains("obstacle"), "{err}");
    let mut d = domain([0.0; 3], [4.0 * dx, 2.0 * dx, 2.0 * dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    let err = refused(d, &[]);
    assert!(err.contains("not supported yet") && err.contains("at least 3"), "{err}");
    // Without a wall face it has nothing to do.
    let mut d = channel();
    d.faces.y_min = PERIODIC;
    d.faces.y_max = PERIODIC;
    let m = model(1.0, 0.01, d, bounced(10.0 * dx, dx, Collision::Mrt));
    let (_, flow) = prepare(&m, &[]);
    assert!(flow.wall_nodes().unwrap().nodes().is_empty());
}

/// With Van Driest's damping the interior nodes are the wall cells, one cell from the wall.
#[test]
fn damping_measures_the_distance_from_the_wall_through_the_nodes() {
    let dx = 0.125;
    let mut d = domain([0.0; 3], [4.0 * dx, 1.0 + dx, 2.0 * dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    let mut p = bounced(10.0 * dx, dx, Collision::Mrt);
    if let Procedure::Fluid { turbulence: Some(Turbulence::Smagorinsky { van_driest, .. }), .. } = &mut p {
        *van_driest = Some(26.0);
    }
    let (_, flow) = prepare(&model(1.0, 0.01, d, p), &[]);
    let grid = *flow.grid();
    let ny = grid.cells[1];
    assert!(flow.van_driest());
    for cell in 0..grid.len() {
        let j = grid.coords(cell)[1];
        let want = (j as f64).min((ny - 1 - j) as f64);
        match flow.wall_distance(cell) {
            Some(d) => assert!(flow.is_fluid(cell) && (d - want).abs() < 1e-12, "{cell}: {d} vs {want}"),
            None => assert!(!flow.is_fluid(cell)),
        }
    }
}

/// The same run on one thread and on several, bit for bit.
#[test]
fn on_node_results_are_bit_for_bit_the_same_on_one_and_many_threads() {
    let (n, re_tau, ub) = (6usize, 180.0, 15.7);
    let dx = 1.0 / n as f64;
    let mut d = domain([0.0; 3], [16.0, 2.0 + dx, 2.0], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.acceleration = [1.0, 0.0, 0.0];
    d.initial_velocity = [ub, 0.0, 0.0];
    let dt = 0.1 * dx / ub;
    for collision in [Collision::Mrt, Collision::Hrr] {
        let m = model(1.0, 1.0 / re_tau, d.clone(), bounced(30.0 * dt, dt, collision));
        let solve =
            |threads: usize| fieldcraft_parallel::with_threads(threads, || fieldcraft_lbm::solve_fluid(&m, None, &[], &Limits::default()).unwrap());
        let one = solve(1);
        assert!(one.grid.len() > 3 * 4096, "several parallel blocks: {}", one.grid.len());
        assert_eq!(one, solve(7));
        assert_eq!(one, solve(3));
    }
}
