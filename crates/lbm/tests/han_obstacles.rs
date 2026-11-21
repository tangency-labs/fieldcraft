//! Han, Ooka and Kikumoto's on-node wall-function bounce (`spalding_bounce`) on obstacle walls and
//! where a wall meets an inlet or outlet.

mod common;

use common::*;
use fieldcraft_lbm::Flow;
use fieldcraft_lbm::lattice::Q;
use fieldcraft_model::{Collision, FaceCondition, HrrSigma, Procedure, Turbulence, WallModel};

/// The force groups of every face of obstacle `body`.
fn body_groups(flow: &Flow, body: &str) -> Vec<usize> {
    flow.groups().iter().enumerate().filter(|(_, g)| g.body.as_deref() == Some(body)).map(|(k, _)| k).collect()
}

const PERIODIC: FaceCondition = FaceCondition::Periodic {};

/// A perturbed coarse channel (4 cells per half-height, Re_τ 180, periodic in x and z), `steps`
/// steps of `turbulence` on `collision`. Its walls are the y box faces, or with `slabs` two
/// obstacle slabs one cell thick in a box periodic in y too, one cell row more at each end: then
/// the boundary nodes are the fluid rows next to the slabs, at the same heights as the box faces'
/// node rows (the walls at y = Δx/2 and 2 + Δx/2 either way).
fn channel(collision: Collision, turbulence: Turbulence, steps: u64, slabs: bool) -> Flow {
    use std::f64::consts::PI;
    let (n, re_tau, ub) = (4usize, 180.0, 15.7);
    let dx = 1.0 / n as f64;
    let (lx, lz) = (12.0 * dx, 6.0 * dx);
    let mut d = if slabs { domain([0.0, -dx, 0.0], [lx, 2.0 + 2.0 * dx, lz], dx) } else { domain([0.0; 3], [lx, 2.0 + dx, lz], dx) };
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    let obstacles = if slabs {
        d.faces.y_min = PERIODIC;
        d.faces.y_max = PERIODIC;
        vec![
            obstacle("Floor", cuboid([-1.0, -dx, -1.0], [lx + 1.0, 0.0, lz + 1.0]), &CUBOID_FACES),
            obstacle("Roof", cuboid([-1.0, 2.0 + dx, -1.0], [lx + 1.0, 2.0 + 2.0 * dx, lz + 1.0]), &CUBOID_FACES),
        ]
    } else {
        Vec::new()
    };
    d.obstacles = obstacles.iter().map(|o| o.body.clone()).collect();
    d.acceleration = [1.0, 0.0, 0.0];
    let dt = 0.1 * dx / ub;
    let mut p = collide(steps as f64 * dt, dt, None, collision);
    if let Procedure::Fluid { turbulence: t, sigma, .. } = &mut p {
        *t = Some(turbulence);
        if collision == Collision::Hrr {
            *sigma = Some(HrrSigma::Value(0.99));
        }
    }
    let (setup, mut flow) = prepare(&model(1.0, 1.0 / re_tau, d, p), &obstacles);
    let grid = *flow.grid();
    let to_lattice = 1.0 / setup.velocity_scale();
    let (kx, kz) = (2.0 * PI / lx, 2.0 * PI / lz);
    flow.set_equilibrium(|cell| {
        let [x, y, z] = grid.centre(cell);
        let eta = y - 1.0 - 0.5 * dx;
        let wall = (1.0 - eta * eta).max(0.0);
        let u = 1.15 * ub * (1.0 - eta.abs().min(1.0).powi(8)) + 2.0 * (2.0 * kz * z).sin() * (kx * x).cos() * wall;
        let v = 1.5 * wall * (2.0 * kz * z).cos() * (3.0 * kx * x).sin();
        let w = 1.5 * wall * eta * (kx * x).sin() * (kz * z).cos();
        (1.0, [u * to_lattice, v * to_lattice, w * to_lattice])
    });
    for _ in 0..steps {
        flow.step().unwrap();
    }
    flow
}

fn smagorinsky() -> Turbulence {
    Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: Some(WallModel::SpaldingBounce) }
}

fn spalart_allmaras() -> Turbulence {
    Turbulence::SpalartAllmaras { prandtl: 0.9, wall_model: Some(WallModel::SpaldingBounce), inflow_ratio: 3.0 }
}

/// A channel whose walls are obstacle slabs is the box-face channel, bit for bit: every
/// population of every cell (the slabs' rows aside), each wall's force against its slab's, and
/// with Spalart–Allmaras the wall function's (y, u_τ) at every first fluid node.
#[test]
fn an_obstacle_channel_is_the_box_face_channel_bit_for_bit() {
    let cases: [(&str, Collision, Turbulence); 4] = [
        ("smagorinsky trt", Collision::Trt, smagorinsky()),
        ("smagorinsky mrt", Collision::Mrt, smagorinsky()),
        ("smagorinsky hrr", Collision::Hrr, smagorinsky()),
        ("spalart-allmaras mrt", Collision::Mrt, spalart_allmaras()),
    ];
    for (name, collision, turbulence) in cases {
        let (boxed, slabbed) = (channel(collision, turbulence, 200, false), channel(collision, turbulence, 200, true));
        let (gb, gs) = (*boxed.grid(), *slabbed.grid());
        assert_eq!(gs.cells[1], gb.cells[1] + 2, "{name}");
        let nodes = |f: &Flow| f.wall_nodes().unwrap().nodes().len();
        assert_eq!(nodes(&boxed), nodes(&slabbed), "{name}");
        assert!(slabbed.wall_nodes().unwrap().nodes().iter().all(|w| w.flat), "{name}: the slabs' nodes are flat");
        let shifted = |cell: usize| {
            let [i, j, k] = gb.coords(cell);
            gs.index([i, j + 1, k])
        };
        for cell in 0..gb.len() {
            for q in 0..Q {
                let (a, b) = (boxed.population(q, cell), slabbed.population(q, shifted(cell)));
                assert_eq!(a.to_bits(), b.to_bits(), "{name}: cell {cell} q {q}: {a} vs {b}");
            }
        }
        for (face, body) in [(2usize, "Floor"), (3, "Roof")] {
            let want = boxed.forces()[face];
            let groups = body_groups(&slabbed, body);
            let got: [f64; 3] = [0, 1, 2].map(|k| groups.iter().map(|g| slabbed.forces()[*g][k]).sum());
            assert_eq!(want.map(f64::to_bits), got.map(f64::to_bits), "{name}: {body} {want:?} vs {got:?}");
        }
        let (bw, sw) = (boxed.sa_wall_function().1, slabbed.sa_wall_function().1);
        assert_eq!(bw.len(), sw.len(), "{name}");
        for (a, b) in bw.iter().zip(sw) {
            assert_eq!((shifted(a.0), a.1.to_bits(), a.2.to_bits()), (b.0, b.1.to_bits(), b.2.to_bits()), "{name}");
        }
        if name.starts_with("spalart") {
            assert!(!bw.is_empty(), "the wall function reads the slabs' nodes");
        }
    }
}

/// A closed periodic box (2D, one cell in z) with a free plate one cell thick (two edges) and an
/// L of two blocks (convex and concave corners), driven both ways along x, on every collision:
/// the fluid's mass and momentum budget close (the populations in flight at the nodes counted),
/// and every node's Han pairs carry exactly the stress it applies, the rest of its tangential
/// momentum being its edge links'.
#[test]
fn plates_and_corners_keep_mass_and_carry_exactly_their_stress() {
    let dx = 1.0 / 16.0;
    for collision in [Collision::Trt, Collision::Regularized, Collision::Mrt, Collision::Hrr] {
        for sign in [1.0, -1.0] {
            let mut d = domain([0.0; 3], [48.0 * dx, 24.0 * dx, dx], dx);
            for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
                *f = PERIODIC;
            }
            let obstacles = vec![
                obstacle("Plate", cuboid([10.0 * dx, 6.0 * dx, -1.0], [30.0 * dx, 7.0 * dx, 1.0]), &CUBOID_FACES),
                obstacle("Base", cuboid([30.0 * dx, 14.0 * dx, -1.0], [42.0 * dx, 17.0 * dx, 1.0]), &CUBOID_FACES),
                obstacle("Post", cuboid([36.0 * dx, 17.0 * dx, -1.0], [42.0 * dx, 21.0 * dx, 1.0]), &CUBOID_FACES),
            ];
            d.obstacles = obstacles.iter().map(|o| o.body.clone()).collect();
            d.acceleration = [sign, 0.0, 0.0];
            d.initial_velocity = [sign, 0.0, 0.0];
            let dt = 0.05 * dx;
            let mut p = collide(300.0 * dt, dt, None, collision);
            if let Procedure::Fluid { turbulence: t, sigma, .. } = &mut p {
                *t = Some(smagorinsky());
                if collision == Collision::Hrr {
                    *sigma = Some(HrrSigma::Value(0.99));
                }
            }
            let (setup, mut flow) = prepare(&model(1.0, 0.02, d, p), &obstacles);
            let nodes = flow.wall_nodes().unwrap().nodes().to_vec();
            assert!(nodes.iter().any(|w| !w.flat) && nodes.iter().any(|w| w.flat), "edges and plane walls both");
            assert!(nodes.iter().any(|w| !w.sampled()), "a concave corner node");
            let a = flow.acceleration();
            let (mut p0, m0) = inventory(&flow);
            let (mut budget, mut mass, mut pairs, mut edge_share) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
            for _ in 0..setup.steps {
                let bulk = flow.mass();
                flow.step().unwrap();
                let (next, m) = inventory(&flow);
                let wall: [f64; 3] = [0, 1, 2].map(|k| flow.forces().iter().map(|f| f[k]).sum());
                for k in 0..3 {
                    budget = budget.max((next[k] - p0[k] - (a[k] * bulk - wall[k])).abs() / next[0].abs());
                }
                mass = mass.max((m / m0 - 1.0).abs());
                let (mut applied, mut edge) = (0.0f64, 0.0f64);
                for w in flow.wall_nodes().unwrap().nodes() {
                    for b in (0..3).filter(|b| *b != usize::from(w.axis)) {
                        pairs = pairs.max((w.taken[b] - w.edge[b] - w.stress[b]).abs());
                        applied += w.stress[b].abs();
                        edge += w.edge[b].abs();
                    }
                }
                edge_share = edge_share.max(edge / (applied + edge));
                p0 = next;
            }
            println!(
                "{} {sign:+}: momentum budget {budget:.1e} of the momentum, mass {mass:.1e}, pairs vs stress {pairs:.1e}, edge share of the tangential wall momentum ≤ {edge_share:.2}",
                collision.name()
            );
            assert!(budget < 1e-11 && mass < 1e-12 && pairs < 1e-15, "{budget} {mass} {pairs}");
        }
    }
}

/// The fluid's momentum and mass: its cells' populations, and the populations that have arrived
/// at each boundary node from the fluid along its Han and edge links (a step ends with the
/// streaming, so the node holds them; its rule sends their sum back in the next step).
fn inventory(flow: &Flow) -> ([f64; 3], f64) {
    use fieldcraft_lbm::bounce::Role;
    use fieldcraft_lbm::lattice::{C, opposite};
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
        for q in (1..Q).filter(|q| matches!(w.roles[*q], Role::Han | Role::Edge)) {
            add(w.cell as usize, opposite(q));
        }
    }
    (p, m)
}

/// Hostile geometry: gaps of one and two cells, a one-cell obstacle, an obstacle that leaves one
/// column of fluid, one that fills the box, one thinner than a cell, and an obstacle reaching an
/// inlet. Each is refused with a message or runs 20 steps (stopping with an error at worst), with
/// no panic; wherever it runs, the nodes keep the fluid's mass.
#[test]
fn hostile_obstacle_geometry_is_refused_or_runs_without_panic() {
    let dx = 1.0 / 8.0;
    let (lx, ly) = (16.0 * dx, 12.0 * dx);
    // (name, its blocks, an inlet and outlet, refused for leaving no fluid cell)
    let boxes: Vec<(&str, Vec<([f64; 3], [f64; 3])>, bool)> = vec![
        ("gap of one cell", vec![([-1.0, -1.0, -1.0], [lx + 1.0, 5.0 * dx, 1.0]), ([-1.0, 6.0 * dx, -1.0], [lx + 1.0, ly + 1.0, 1.0])], false),
        ("gap of two cells", vec![([-1.0, -1.0, -1.0], [lx + 1.0, 5.0 * dx, 1.0]), ([-1.0, 7.0 * dx, -1.0], [lx + 1.0, ly + 1.0, 1.0])], false),
        ("one cell", vec![([4.0 * dx, 4.0 * dx, -1.0], [5.0 * dx, 5.0 * dx, 1.0])], false),
        ("one column left", vec![([-1.0, -1.0, -1.0], [7.0 * dx, ly + 1.0, 1.0]), ([8.0 * dx, -1.0, -1.0], [lx + 1.0, ly + 1.0, 1.0])], false),
        ("fills the box", vec![([-1.0, -1.0, -1.0], [lx + 1.0, ly + 1.0, 1.0])], false),
        ("thinner than a cell", vec![([4.0 * dx, 4.1 * dx, -1.0], [9.0 * dx, 4.2 * dx, 1.0])], false),
        ("at an inlet", vec![([-1.0, -1.0, -1.0], [4.0 * dx, 3.0 * dx, 1.0])], true),
    ];
    for (name, cuboids, inlet) in boxes {
        let mut d = domain([0.0; 3], [lx, ly, dx], dx);
        for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
            *f = PERIODIC;
        }
        if inlet {
            d.faces.x_min = FaceCondition::Velocity { velocity: [1.0, 0.0, 0.0], profile: Default::default() };
            d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: fieldcraft_model::Outflow::Uniform };
        }
        let obstacles: Vec<fieldcraft_lbm::Obstacle> =
            cuboids.iter().enumerate().map(|(k, (lo, hi))| obstacle(&format!("Block{k}"), cuboid(*lo, *hi), &CUBOID_FACES)).collect();
        d.obstacles = obstacles.iter().map(|o| o.body.clone()).collect();
        d.acceleration = [1.0, 0.0, 0.0];
        d.initial_velocity = [1.0, 0.0, 0.0];
        let dt = 0.05 * dx;
        let mut p = collide(20.0 * dt, dt, None, Collision::Mrt);
        if let Procedure::Fluid { turbulence: t, .. } = &mut p {
            *t = Some(spalart_allmaras());
        }
        let m = model(1.0, 0.02, d, p);
        let built = fieldcraft_lbm::Setup::new(&m, None, &obstacles, &fieldcraft_lbm::Limits::default())
            .and_then(|setup| Flow::new(&setup, &obstacles).map(|flow| (setup, flow)));
        let (setup, mut flow) = match built {
            Err(e) => {
                let e = e.to_string();
                println!("{name}: refused: {e}");
                if inlet {
                    assert!(e.contains("not supported yet") && e.contains("x_min"), "{name}: {e}");
                } else {
                    assert!(e.contains("leaves no fluid cell") && e.contains("refine"), "{name}: {e}");
                    assert!(["gap of one cell", "gap of two cells", "one column left"].contains(&name), "{name}: {e}");
                }
                continue;
            }
            Ok(b) => b,
        };
        assert!(!inlet, "{name}: an obstacle node at an inlet runs");
        assert!(!["gap of one cell", "gap of two cells", "one column left"].contains(&name), "{name} runs with no fluid cell");
        let (_, m0) = inventory(&flow);
        let mut stopped = None;
        for _ in 0..setup.steps {
            if let Err(e) = flow.step() {
                stopped = Some(e.to_string());
                break;
            }
        }
        let (_, m1) = inventory(&flow);
        let nodes = flow.wall_nodes().map_or(0, |n| n.nodes().len());
        println!("{name}: {nodes} nodes, mass {:.1e}, {}", (m1 / m0 - 1.0).abs(), stopped.as_deref().unwrap_or("ran"));
        if stopped.is_none() && m0 > 0.0 {
            assert!((m1 / m0 - 1.0).abs() < 1e-12, "{name}: mass {m0} → {m1}");
        }
    }
}
