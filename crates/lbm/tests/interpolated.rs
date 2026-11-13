//! Interpolated walls (`walls: "interpolated"`, Bouzidi, Firdaouss & Lallemand 2001) against
//! halfway bounce-back: the order of accuracy of Poiseuille flow in a round pipe (theory: the
//! closed-form parabola), a plane channel whose walls sit at arbitrary fractions of a cell
//! (Bouzidi et al.'s own test), the Stokes drag of a periodic array of spheres at 3 cells per
//! radius (theory: Hasimoto 1959 with Sangani & Acrivos 1982, issue #798), the links that keep
//! halfway bounce-back, and the refusals.
//!
//! Each test prints its errors (`cargo test -p fieldcraft-lbm --test interpolated -- --nocapture`).

mod common;

use common::*;
use fieldcraft_lbm::{Error, Limits, Obstacle, Setup, WallLinks, solve_fluid};
use fieldcraft_mesh::tet::geometry::AnalyticSurface;
use fieldcraft_model::{BoxFace, FaceCondition, Procedure, WallTreatment};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};
const HALFWAY: WallTreatment = WallTreatment::Halfway;
const INTERPOLATED: WallTreatment = WallTreatment::Interpolated;

/// Poiseuille flow in a pipe of radius 1 along x (a tube obstacle with its exact inner cylinder),
/// driven by a body force, periodic along x and one cell deep; `k` cells per radius. Returns the
/// relative L2 error of the velocity over the box's cross-section, the momentum balance (the wall
/// force over the driving force, less 1), and the wall links.
fn pipe(k: usize, w: WallTreatment) -> (f64, f64, Option<WallLinks>) {
    let (r, nu, umax) = (1.0, 1.0 / 6.0, 0.1);
    let dx = r / k as f64;
    let half = r + 2.0 * dx;
    let mut d = domain([0.0, -half, -half], [dx, half, half], dx);
    d.faces.x_min = PERIODIC;
    d.faces.x_max = PERIODIC;
    let a = 4.0 * nu * umax / (r * r);
    d.acceleration = [a, 0.0, 0.0];
    d.obstacles = vec!["Tube".into()];
    let mut tube = obstacle("Tube", tube(-1.0, dx + 1.0, [0.0, 0.0], r, 3.0 * half, 512), &["inner", "outer", "start", "end"]);
    tube.exact = vec![Some(AnalyticSurface::Cylinder { origin: [0.0; 3], axis: [1.0, 0.0, 0.0], radius: r }), None, None, None];
    let m = model(1.0, nu, d, with_walls(step(5.0 * r * r / nu, dx * dx, None), w));
    let (setup, mut flow) = prepare(&m, std::slice::from_ref(&tube));
    run_all(&setup, &mut flow);
    let grid = *flow.grid();
    let (mut e2, mut u2) = (0.0, 0.0);
    for c in 0..grid.len() {
        let [_, y, z] = grid.centre(c);
        let rr = y * y + z * z;
        let exact = if rr < r * r { a / (4.0 * nu) * (r * r - rr) } else { 0.0 };
        let u = flow.moments(c).map_or(0.0, |(_, u)| u[0] * setup.velocity_scale());
        e2 += (u - exact) * (u - exact);
        u2 += exact * exact;
    }
    let wall: f64 = flow.groups().iter().zip(flow.forces()).filter(|(g, _)| g.body.is_some()).map(|(_, f)| f[0]).sum::<f64>() * setup.force_scale();
    // The body force acts on the fluid's mass, which interpolated walls do not conserve.
    let drive = a * flow.mass() * dx.powi(3);
    ((e2 / u2).sqrt(), wall / drive - 1.0, flow.wall_links())
}

/// The pipe on `k` cells per radius with both wall treatments, checked: every link into the
/// tube finds the cylinder and has a fluid cell behind it, the wall takes the driving force, and
/// the interpolation is several times closer. The L2 errors, halfway then interpolated.
fn pipe_pair(k: usize) -> (f64, f64) {
    let (h, balance_h, none) = pipe(k, HALFWAY);
    let (i, balance_i, links) = pipe(k, INTERPOLATED);
    assert!(none.is_none());
    let links = links.unwrap();
    println!("pipe, {k} cells per radius: L2 halfway {h:.3e}, interpolated {i:.3e}; balance {balance_h:.1e}, {balance_i:.1e}; links {links:?}");
    assert!(links.interpolated > 0 && links.halfway() == 0, "{links:?}");
    // At steady state the wall takes the driving force (with interpolated walls to how far the
    // mass has drifted since).
    assert!(balance_h.abs() < 1e-6 && balance_i.abs() < 1e-4, "{balance_h} {balance_i}");
    assert!(i < 0.3 * h, "{k}: {i} {h}");
    (h, i)
}

/// The pipe on each grid with both wall treatments ([`pipe_pair`]): the L2 errors, halfway then
/// interpolated, and their fitted orders.
fn pipe_study(cells: &[usize]) -> (Vec<f64>, Vec<f64>, f64, f64) {
    let (mut halfway, mut interpolated) = (Vec::new(), Vec::new());
    for &k in cells {
        let (h, i) = pipe_pair(k);
        halfway.push(h);
        interpolated.push(i);
    }
    let (oh, oi) = (fitted_order(cells, &halfway), fitted_order(cells, &interpolated));
    println!("pipe, {cells:?} cells per radius: halfway {} (order {oh:.2}), interpolated {} (order {oi:.2})", sci(&halfway), sci(&interpolated));
    (halfway, interpolated, oh, oi)
}

/// The round pipe on 8 cells per radius (τ = 1): interpolated walls put the wall on the cylinder,
/// several times closer than halfway bounce-back's staircase (5.8e-3 against 4.2e-2).
#[test]
fn a_coarse_pipe_is_closer_with_interpolated_walls() {
    pipe_pair(8);
}

/// The round pipe on 4 to 16 cells per radius (τ = 1): halfway bounce-back puts the wall on a
/// staircase (about first order), interpolated walls put it on the cylinder (about second).
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn pipe_poiseuille_converges_at_second_order_with_interpolated_walls() {
    let (_, _, halfway, interpolated) = pipe_study(&[4, 6, 8, 12, 16]);
    assert!(interpolated > 1.5, "interpolated order {interpolated}");
    assert!(halfway < 1.3, "halfway order {halfway}");
}

/// The same study on eight grids up to 32 cells per radius.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn pipe_poiseuille_orders_on_eight_grids() {
    let (_, _, halfway, interpolated) = pipe_study(&[4, 6, 8, 12, 16, 20, 24, 32]);
    assert!(interpolated > 1.7, "interpolated order {interpolated}");
    assert!(halfway < 1.3, "halfway order {halfway}");
}

/// A plane channel of height 1 between two slabs, its lower wall `shift` cells above a cell
/// face (the upper one too, the height being whole cells), `n` cells across, driven by a body
/// force, periodic in x and z (one cell each). Returns the largest error of the velocity over
/// its peak, at the fluid cells.
fn shifted_channel(n: usize, shift: f64, w: WallTreatment) -> f64 {
    let (h, nu, umax) = (1.0, 1.0 / 6.0, 0.1);
    let dx = h / n as f64;
    let a = 8.0 * nu * umax / (h * h);
    // Two cells of slab below the channel and above it.
    let (lo, hi) = ((2.0 + shift) * dx, (2.0 + shift) * dx + h);
    let top = (n + 5) as f64 * dx;
    let mut d = domain([0.0; 3], [dx, top, dx], dx);
    for f in [BoxFace::XMin, BoxFace::XMax, BoxFace::ZMin, BoxFace::ZMax] {
        *d.faces.get_mut(f) = PERIODIC;
    }
    d.acceleration = [a, 0.0, 0.0];
    d.obstacles = vec!["Floor".into(), "Ceiling".into()];
    let wide = |y0: f64, y1: f64| cuboid([-1.0, y0, -1.0], [1.0 + dx, y1, 1.0 + dx]);
    let obstacles = [obstacle("Floor", wide(-dx, lo), &CUBOID_FACES), obstacle("Ceiling", wide(hi, top + dx), &CUBOID_FACES)];
    let m = model(1.0, nu, d, with_walls(step(4.0 * h * h / nu, dx * dx, None), w));
    let (setup, mut flow) = prepare(&m, &obstacles);
    run_all(&setup, &mut flow);
    if w == INTERPOLATED {
        let links = flow.wall_links().unwrap();
        assert!(links.halfway() == 0, "{links:?}");
    }
    let mut err: f64 = 0.0;
    for j in 0..n + 5 {
        let y = (j as f64 + 0.5) * dx;
        let Some((_, u)) = flow.moments(j) else { continue };
        let exact = a / (2.0 * nu) * (y - lo) * (hi - y);
        err = err.max((u[0] * setup.velocity_scale() - exact).abs() / umax);
    }
    err
}

/// Bouzidi et al.'s test: wherever the walls sit within a cell, interpolated walls keep the
/// channel's parabola to second order, while halfway bounce-back moves each wall to the nearest
/// cell face (an error of up to half a cell, first order).
#[test]
fn a_channel_with_walls_anywhere_in_a_cell_keeps_its_parabola() {
    // 0: the walls on cell faces, where halfway bounce-back puts them; ½ would put them through
    // cell centres.
    let shifts = [0.0, 0.2, 0.4, 0.6, 0.8];
    let mut worst = Vec::new();
    for n in [8, 16, 32] {
        let halfway: Vec<f64> = shifts.iter().map(|&s| shifted_channel(n, s, HALFWAY)).collect();
        let interpolated: Vec<f64> = shifts.iter().map(|&s| shifted_channel(n, s, INTERPOLATED)).collect();
        println!("channel, {n} cells, shifts {shifts:?}: halfway {}, interpolated {}", sci(&halfway), sci(&interpolated));
        let max = |v: &[f64]| v.iter().copied().fold(0.0, f64::max);
        // Halfway is exact (Λ = 3/16) with the walls on cell faces, and so is the interpolation
        // there (q = ½ is halfway bounce-back).
        assert!(halfway[0] < 1e-9 && interpolated[0] < 1e-9, "{halfway:?} {interpolated:?}");
        assert!(max(&interpolated) < 0.25 * max(&halfway), "{interpolated:?} {halfway:?}");
        worst.push((max(&halfway), max(&interpolated)));
    }
    let h: Vec<f64> = worst.iter().map(|w| w.0).collect();
    let i: Vec<f64> = worst.iter().map(|w| w.1).collect();
    let (oh, oi) = (orders(&h, 2.0), orders(&i, 2.0));
    println!("channel, worst over shifts, 8 → 32 cells: halfway {} (orders {oh:.2?}), interpolated {} (orders {oi:.2?})", sci(&h), sci(&i));
    assert!(oi.iter().all(|o| *o > 1.8), "{oi:?}");
    assert!(oh.iter().all(|o| *o < 1.2), "{oh:?}");
}

/// A sphere of radius `r` cells (with its exact surface) at the centre of a periodic box of `l`
/// cells, driven by a small body force (scaled as R⁻³ from #798's 10⁻⁶ at R = 3) with
/// kinematic viscosity `nu` (lattice units; τ = 3ν + ½, Λ = 3/16), run `steps`: the measured K
/// over Hasimoto's, less 1.
fn sphere_drag(r: f64, l: usize, steps: u64, w: WallTreatment, nu: f64) -> f64 {
    let a = 1e-6 / (r / 3.0).powi(3);
    let side = l as f64;
    let mut d = domain([0.0; 3], [side; 3], 1.0);
    for f in BoxFace::ALL {
        *d.faces.get_mut(f) = PERIODIC;
    }
    d.acceleration = [a, 0.0, 0.0];
    d.obstacles = vec!["Ball".into()];
    let centre = [side / 2.0; 3];
    let mut ball = obstacle("Ball", sphere(centre, r, 4), &["skin"]);
    ball.exact = vec![Some(AnalyticSurface::Sphere { center: centre, radius: r })];
    let p = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: None,
        steps: Some(steps),
        time_step: Some(1.0),
        frames: 1,
        magic: None,
        sigma: None,
        collision: None,
        turbulence: None,
        steady: None,
        walls: w,
    };
    let m = model(1.0, nu, d, p);
    let (setup, mut flow) = prepare(&m, std::slice::from_ref(&ball));
    run_all(&setup, &mut flow);
    let cells = (l * l * l) as f64;
    let u: f64 = (0..l * l * l).filter_map(|c| flow.moments(c)).map(|(_, u)| u[0]).sum::<f64>() / cells;
    let skin = flow.groups().iter().position(|g| g.name == "skin").unwrap();
    let force = flow.forces()[skin][0];
    // Steady: the wall takes the driving force (on the fluid's mass).
    let drive = a * flow.mass();
    assert!((force - drive).abs() < 1e-3 * drive, "{force} vs {drive}");
    let k = force * cells / flow.fluid_cells() as f64 / (6.0 * std::f64::consts::PI * nu * r * u);
    let c = 4.0 / 3.0 * std::f64::consts::PI * r.powi(3) / cells;
    let error = k / hasimoto(c) - 1.0;
    println!(
        "sphere array, R = {r}, L = {l}, τ = {:.3}, {}: K = {k:.5}, Hasimoto {:.5} ({:+.2} %)",
        setup.tau_plus,
        w.name(),
        hasimoto(c),
        100.0 * error
    );
    error
}

/// #798's sphere array at 3 cells per radius, at τ = 1: halfway bounce-back's staircase reads the
/// drag 11.7 % high; interpolated walls put the sphere where it is.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn stokes_drag_of_a_coarse_sphere_array_with_interpolated_walls() {
    let halfway = sphere_drag(3.0, 16, 9000, HALFWAY, 1.0 / 6.0);
    let interpolated = sphere_drag(3.0, 16, 9000, INTERPOLATED, 1.0 / 6.0);
    assert!(halfway > 0.11, "{halfway}");
    assert!(interpolated.abs() < 0.02, "{interpolated}");
}

/// The sphere array at 3, 6 and 9 cells per radius at τ = 1, and at R = 3 over τ from 0.6 to 3
/// (#798's τ): halfway bounce-back does not depend on τ with Λ = 3/16; the linear interpolation's
/// wall moves with it, so its drag does.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn stokes_drag_of_a_sphere_array_with_interpolated_walls_over_radii_and_tau() {
    // Measured: halfway +11.7 %, +1.7 %, +2.3 % (its staircase error is not monotone in R, #798);
    // interpolated −1.64 %, −0.85 %, −0.32 %, closing in on the series every time.
    let mut last = f64::INFINITY;
    for (r, l) in [(3.0, 16), (6.0, 32), (9.0, 48)] {
        // Six viscous times of the box at τ = 1.
        let steps = (36 * l * l) as u64;
        let h = sphere_drag(r, l, steps, HALFWAY, 1.0 / 6.0);
        let i = sphere_drag(r, l, steps, INTERPOLATED, 1.0 / 6.0);
        assert!(i.abs() < 0.02 && i.abs() < last, "R = {r}: halfway {h}, interpolated {i}, previous {last}");
        last = i.abs();
    }
    for nu in [1.0 / 30.0, 1.0 / 12.0, 1.0 / 2.0, 5.0 / 6.0] {
        let steps = (1536.0 * 6.0 * (1.0 / 6.0) / nu) as u64;
        let h = sphere_drag(3.0, 16, steps, HALFWAY, nu);
        let i = sphere_drag(3.0, 16, steps, INTERPOLATED, nu);
        assert!((h - 0.1166).abs() < 1e-3, "halfway does not depend on τ: {h}");
        assert!(i.abs() < h.abs(), "ν = {nu}: halfway {h}, interpolated {i}");
    }
}

/// A channel one cell wide between walls nearer than half a link to its centres: its links
/// have no fluid cell behind them, keep halfway bounce-back, and the solve says how many.
#[test]
fn links_without_a_second_fluid_cell_keep_halfway_and_are_reported() {
    let dx = 0.1;
    let mut d = domain([0.0; 3], [dx, 1.0, dx], dx);
    for f in [BoxFace::XMin, BoxFace::XMax, BoxFace::ZMin, BoxFace::ZMax] {
        *d.faces.get_mut(f) = PERIODIC;
    }
    d.acceleration = [0.01, 0.0, 0.0];
    d.obstacles = vec!["Floor".into(), "Ceiling".into()];
    let wide = |y0: f64, y1: f64| cuboid([-1.0, y0, -1.0], [1.0 + dx, y1, 1.0 + dx]);
    // The cell centred at y = 0.45 is the only fluid one; the walls cross its links at 0.2 and
    // 0.3 of a cell.
    let obstacles = [obstacle("Floor", wide(-dx, 0.43), &CUBOID_FACES), obstacle("Ceiling", wide(0.48, 1.0 + dx), &CUBOID_FACES)];
    let m = model(1.0, 1.0 / 6.0, d, with_walls(step(1.0, dx * dx, None), INTERPOLATED));
    let (_, flow) = prepare(&m, &obstacles);
    let links = flow.wall_links().unwrap();
    // Up, down and the four diagonals up and down: ten links, every one below ½.
    assert_eq!(links, WallLinks { interpolated: 0, no_second_node: 10, no_crossing: 0 }, "{links:?}");
    let r = solve_fluid(&m, None, &obstacles, &Limits::default()).unwrap();
    assert_eq!(r.stats.wall_links, Some(links));
    assert!(r.warnings.iter().any(|w| w.contains("10 of 10 wall links kept halfway bounce-back")), "{:?}", r.warnings);
}

/// Interpolated walls refuse what would rebuild or reflect their links another way.
#[test]
fn interpolated_walls_refuse_wall_models_and_free_surfaces() {
    let base = || {
        let mut d = domain([0.0; 3], [1.0; 3], 0.1);
        d.obstacles = vec![];
        d
    };
    let mut p = with_walls(step(1.0, 0.01, None), INTERPOLATED);
    if let Procedure::Fluid { turbulence, .. } = &mut p {
        *turbulence = Some(fieldcraft_model::Turbulence::Smagorinsky {
            constant: 0.1,
            prandtl: 0.9,
            van_driest: None,
            wall_model: Some(fieldcraft_model::WallModel::Musker),
        });
    }
    let m = model(1.0, 0.01, base(), p);
    let e = Setup::new(&m, None, &[] as &[Obstacle], &Limits::default()).unwrap_err();
    assert!(matches!(&e, Error::Invalid(s) if s.contains("not supported yet: interpolated walls") && s.contains("wall model")), "{e}");
    let mut d = base();
    d.free_surface = Some(fieldcraft_model::FreeSurface { liquid: vec![fieldcraft_model::LiquidRegion { min: [0.0; 3], max: [1.0, 0.5, 1.0] }] });
    let m = model(1.0, 0.01, d, with_walls(step(1.0, 0.01, None), INTERPOLATED));
    let e = Setup::new(&m, None, &[] as &[Obstacle], &Limits::default()).unwrap_err();
    assert!(matches!(&e, Error::Invalid(s) if s.contains("not supported yet: interpolated walls") && s.contains("free surface")), "{e}");
}
