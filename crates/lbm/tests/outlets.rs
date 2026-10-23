//! Pressure outlets that a wake or vortices reach (#599): the convective outlet against the
//! uniform one (Zou–He) on the acceptance cases of #598, a vortex leaving the box, Schäfer–Turek's
//! cylinder at Re 100 with the outlet near and far, and a long run with backflow.
//!
//! Each test prints its numbers (`cargo test -p fieldcraft-lbm --test outlets -- --nocapture`);
//! the heavy ones run with `--include-ignored`.

mod common;

use common::*;
use fieldcraft_lbm::lattice::{Q, equilibrium};
use fieldcraft_lbm::{Flow, FluidResult, Limits, Obstacle, Setup, solve_fluid};
use fieldcraft_model::{BoxFace, Collision, FaceCondition, Model, Outflow, Procedure, Profile};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};
const BOTH: [Outflow; 2] = [Outflow::Uniform, Outflow::Convective];

fn fluid(end_time: f64, time_step: f64, frames: u32) -> Procedure {
    fluid_with(end_time, time_step, frames, None)
}

/// [`fluid`] with a collision operator.
fn fluid_with(end_time: f64, time_step: f64, frames: u32, collision: Option<Collision>) -> Procedure {
    Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(end_time),
        steps: None,
        time_step: Some(time_step),
        frames,
        magic: None,
        steady: None,
        collision,
        turbulence: None,
    }
}

/// Schäfer–Turek's channel (height 0.41, cylinder of diameter 0.1 at (0.2, 0.2), parabolic inlet
/// with peak `um`, ν = 10⁻³; quasi-2D: one periodic cell across z), the outlet at `l`, `segments`
/// facets round the cylinder: the model and the cylinder.
fn channel_model(l: f64, dx: f64, um: f64, outflow: Outflow, procedure: Procedure, segments: usize, initial: bool) -> (Model, Obstacle) {
    let (h, d, nu) = (0.41, 0.1, 1e-3);
    let mut dom = domain([0.0; 3], [l, h, dx], dx);
    dom.faces.x_min = FaceCondition::Velocity { velocity: [um, 0.0, 0.0], profile: Profile::Parabolic };
    dom.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow };
    dom.faces.z_min = PERIODIC;
    dom.faces.z_max = PERIODIC;
    if initial {
        // Started at the mean inlet speed, not from rest: the impulsive start's pressure waves
        // (ρ c u, far above the dynamic pressure) reflect between the inlet and outlet for a long
        // time in the weakly compressible scheme and would swamp the forces.
        dom.initial_velocity = [2.0 / 3.0 * um, 0.0, 0.0];
    }
    dom.obstacles = vec!["Cylinder".into()];
    let cyl = obstacle("Cylinder", cylinder([0.2, 0.2], d / 2.0, -0.1, 0.1, segments), &["side", "bottom", "top"]);
    (model(1.0, nu, dom, procedure), cyl)
}

/// [`channel_model`] solved.
fn channel(l: f64, dx: f64, um: f64, outflow: Outflow, procedure: Procedure, segments: usize, initial: bool) -> fieldcraft_lbm::Result<FluidResult> {
    let (m, cyl) = channel_model(l, dx, um, outflow, procedure, segments, initial);
    solve_fluid(&m, None, std::slice::from_ref(&cyl), &Limits::default())
}

/// The largest departure of u_x in layer `i` from the line through the two layers before it.
fn departure(r: &FluidResult, i: usize) -> f64 {
    let last = r.frames.last().unwrap();
    let ux = |i: usize, j: usize| f64::from(last.velocity[r.grid.index([i, j, 0])][0]);
    (0..r.grid.cells[1]).map(|j| (ux(i, j) - (2.0 * ux(i - 1, j) - ux(i - 2, j))).abs()).fold(0.0, f64::max)
}

/// #598's wake case: Schäfer–Turek 2D-1 (Re 20) cut 1.5 diameters behind the cylinder, where
/// the wake has not recovered. The uniform outlet puts a step into the last layer's velocity
/// (`a_wake_at_a_pressure_outlet_distorts_only_the_last_layer` in `verification.rs`: 9.2e-3
/// against 5.3e-3 inside); the convective one carries the wake out smoothly.
#[test]
fn a_wake_leaves_a_convective_outlet_without_a_step() {
    let (um, dx) = (0.3, 0.01);
    let r = channel(0.4, dx, um, Outflow::Convective, fluid(4.0, 0.025 * dx / um, 1), 256, false).unwrap();
    let nx = r.grid.cells[0];
    let (inside, outlet) = (departure(&r, nx - 2), departure(&r, nx - 1));
    println!("wake at a convective outlet: departure inside {inside:.2e}, outlet layer {outlet:.2e} (of {um})");
    assert!(outlet < 2e-3 && inside < 2e-3, "{inside} {outlet}");
}

/// #598's backflow case: Schäfer–Turek's channel at Re 100 cut 1 diameter behind the cylinder
/// (l = 0.3): the shedding wake reaches the outlet with reversed flow in some of its cells, which
/// hold the pressure (the backflow guard); the run stays stable and conserves mass, as the
/// uniform outlet does (`a_short_channel_at_re_100_survives_backflow_at_the_outlet`).
#[test]
fn a_convective_outlet_survives_backflow_at_re_100() {
    let (um, dx) = (1.5, 0.01);
    let r = channel(0.3, dx, um, Outflow::Convective, fluid(4.0, 0.05 * dx / um, 8), 256, false).unwrap();
    let [nx, ny, _] = r.grid.cells;
    let reversed = r.frames.iter().map(|f| (0..ny).filter(|j| f.velocity[r.grid.index([nx - 1, *j, 0])][0] < 0.0).count()).max().unwrap_or(0);
    let change = r.frames.last().unwrap().mass / r.frames[0].mass - 1.0;
    println!("Re 100 short channel, convective outlet: up to {reversed} reversed outlet cells, mass change {change:.2e}");
    assert!(reversed > 0, "no reversed outlet cells: the wake no longer reaches the outlet");
    assert!(change.abs() < 0.01, "{change}");
}

/// The same backflow case with the regularised collision. A convective outlet cell where the flow
/// comes back in takes the uniform outlet's closure ([`fieldcraft_lbm::outlet::backflow`]), which
/// under this collision holds the equilibrium of the face's density and the Zou–He normal
/// velocity; the cells being convected keep their non-equilibrium part (no override). Checked cell
/// by cell every 100 steps, and the convect path itself stays stable with the reversed cells
/// beside it.
#[test]
fn a_convective_outlet_survives_backflow_at_re_100_with_the_regularized_collision() {
    let (um, dx) = (1.5, 0.01);
    let (m, cyl) = channel_model(0.3, dx, um, Outflow::Convective, fluid_with(4.0, 0.05 * dx / um, 1, Some(Collision::Regularized)), 256, false);
    let (setup, mut flow) = prepare(&m, std::slice::from_ref(&cyl));
    let [nx, ny, _] = flow.grid().cells;
    let mass0 = flow.mass();
    let (mut held, mut convected, mut worst_held, mut least_convected) = (0, 0, 0.0f64, f64::MAX);
    for k in 1..=setup.steps {
        flow.step().unwrap();
        if k % 100 != 0 {
            continue;
        }
        for j in 0..ny {
            let cell = flow.grid().index([nx - 1, j, 0]);
            let Some((rho, u)) = flow.moments(cell) else { continue };
            // The largest departure of the populations from the equilibrium of the cell's state.
            let off = (0..Q).map(|q| (flow.population(q, cell) - equilibrium(q, rho, u)).abs()).fold(0.0, f64::max);
            if (rho - 1.0).abs() < 1e-12 && off < 1e-12 {
                held += 1;
                worst_held = worst_held.max((rho - 1.0).abs());
            } else {
                convected += 1;
                least_convected = least_convected.min(off);
            }
            if u[0] < 0.0 {
                // A reversed cell is held: the face's density, at equilibrium.
                assert!((rho - 1.0).abs() < 1e-12 && off < 1e-12, "step {k}, j {j}: ρ {rho}, off equilibrium by {off}");
            }
        }
    }
    let change = flow.mass() / mass0 - 1.0;
    println!(
        "Re 100 short channel, regularized collision, convective outlet: {held} held and {convected} convected outlet cell samples (the convected ones at least {least_convected:.1e} off equilibrium), mass change {change:.2e}"
    );
    assert!(held > 0 && convected > held, "{held} {convected}");
    assert!(least_convected > 1e-9, "{least_convected}");
    assert!(change.abs() < 0.01, "{change}");
}

/// Plane Poiseuille flow driven by the pressures alone (Δp = 12 μ Ū L / H²), both ends convective
/// pressure faces; plates 1 apart, 3 long, 16 cells across. The inlet end takes flow in
/// everywhere, so it is held like a uniform face; the outlet end holds its pressure cell by cell.
/// The profile is the parabola within 3 % and each end holds its pressure, as with uniform faces
/// (`pressure_driven_poiseuille_holds_its_end_pressures`).
#[test]
fn pressure_driven_poiseuille_holds_its_end_pressures_with_convective_faces() {
    let (h, l, nu, umax, n) = (1.0, 3.0, 0.05, 0.1, 16);
    let dx = h / n as f64;
    let dp = 12.0 * nu * (2.0 / 3.0 * umax) * l / (h * h);
    let mut d = domain([0.0; 3], [l, h, dx], dx);
    d.faces.x_min = FaceCondition::Pressure { pressure: dp, outflow: Outflow::Convective };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective };
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    let r = solve_fluid(&model(1.0, nu, d, fluid(4.0 * h * h / nu, 0.1 / 3.0 * dx * dx / nu, 1)), None, &[], &Limits::default()).unwrap();
    let last = r.frames.last().unwrap();
    let nx = r.grid.cells[0];
    let mut err: f64 = 0.0;
    for j in 0..n {
        let y = (j as f64 + 0.5) * dx;
        let exact = 4.0 * umax * y * (h - y) / (h * h);
        err = err.max((f64::from(last.velocity[r.grid.index([nx / 2, j, 0])][0]) - exact).abs() / umax);
    }
    let p = |i: usize| f64::from(last.pressure[r.grid.index([i, n / 2, 0])]);
    println!("pressure-driven Poiseuille, convective faces: profile error {err:.3e}, end pressures {:.5} / {:.5} of {dp:.5} / 0", p(0), p(nx - 1));
    assert!(err < 0.03, "{err}");
    assert!((p(0) / dp - 1.0).abs() < 1e-3 && (p(nx - 1) / dp).abs() < 1e-3, "{} {}", p(0), p(nx - 1));
}

/// A 16 × 8 box with an inlet at x_min, a wall at y_min and pressure faces at x_max and y_max
/// (both at gauge pressure 0) that meet in an edge. The edge cells belong to one face (x_max:
/// the first pressure face in `BoxFace::ALL` order), which holds them at the outlet density to
/// round-off whatever the two faces' outflow kinds; the corner with the inlet belongs to the inlet.
/// Two faces each applying their own rule to the edge cell left it 17 % off the outlet density.
#[test]
fn an_edge_of_two_pressure_faces_holds_the_outlet_density() {
    for (x_max, y_max) in
        [(Outflow::Convective, Outflow::Convective), (Outflow::Uniform, Outflow::Convective), (Outflow::Convective, Outflow::Uniform)]
    {
        let mut d = domain([0.0; 3], [16.0, 8.0, 1.0], 1.0);
        d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Uniform };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: x_max };
        d.faces.y_max = FaceCondition::Pressure { pressure: 0.0, outflow: y_max };
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        let m = model(1.0, 0.05, d, fluid(800.0, 1.0, 1));
        let (setup, mut flow) = prepare(&m, &[]);
        assert_eq!(setup.steps, 800);
        run_all(&setup, &mut flow);
        let edge = flow.grid().index([15, 7, 0]);
        let (rho, _) = flow.moments(edge).unwrap();
        let corner = flow.grid().index([0, 7, 0]);
        let (_, u) = flow.moments(corner).unwrap();
        println!("edge of x_max {x_max:?} and y_max {y_max:?}: ρ − ρ_out {:.2e}; inlet corner u_x {:.6}", rho - 1.0, u[0]);
        assert!((rho - 1.0).abs() < 1e-12, "{rho}");
        assert!((u[0] - 0.05).abs() < 1e-12, "{u:?}");
    }
}

/// A vortex (stream function `A exp(−r²/2R²)`, R = 0.12, peak swirl `swirl` U) carried by a
/// uniform stream U through a box `l` long with periodic sides, from x = 0.7 out through the
/// outlet: the setup and the flow at its initial state. A box 3 long, whose outlet the vortex
/// never reaches, is the reference: the largest difference from it over the short box is what the
/// outlet reflects.
fn vortex(l: f64, outflow: Outflow, swirl: f64) -> (Setup, Flow) {
    let (h, u0, nu, x0, r, n) = (1.0, 0.05, 1e-4, 0.7, 0.12, 50);
    let dx = h / n as f64;
    let mut d = domain([0.0; 3], [l, h, dx], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [u0, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow };
    for f in [BoxFace::YMin, BoxFace::YMax, BoxFace::ZMin, BoxFace::ZMax] {
        *d.faces.get_mut(f) = PERIODIC;
    }
    let m = model(1.0, nu, d, fluid(1.0 / u0, dx, 1));
    let (setup, mut flow) = prepare(&m, &[]);
    let grid = *flow.grid();
    // The swirl peaks at r = R: A / R e^{−1/2}.
    let a = swirl * u0 * r / (-0.5f64).exp();
    let to = 1.0 / setup.velocity_scale();
    flow.set_equilibrium(|c| {
        let [x, y, _] = grid.centre(c);
        let (dx, dy) = (x - x0, y - 0.5 * h);
        let e = (-(dx * dx + dy * dy) / (2.0 * r * r)).exp();
        // The radial balance dp/dr = ρ u_θ² / r of this swirl.
        let p = -a * a / (2.0 * r * r) * e * e;
        (setup.lattice_density(p), [(u0 - a * dy / (r * r) * e) * to, a * dx / (r * r) * e * to, 0.0])
    });
    (setup, flow)
}

/// Per outlet: the largest difference from the long box (of U) over the cells more than 5 layers
/// from the outlet, the same at the end, the reversed outlet cell samples, and the step it
/// diverged at.
fn vortex_runs(swirl: f64) -> Vec<(Outflow, f64, f64, usize, Option<u64>)> {
    let (setup, mut reference) = vortex(3.0, Outflow::Uniform, swirl);
    let mut runs: Vec<(Outflow, Flow, f64, f64, usize, Option<u64>)> =
        BOTH.into_iter().map(|o| (o, vortex(1.0, o, swirl).1, 0.0, 0.0, 0, None)).collect();
    let u0 = 0.05 / setup.velocity_scale();
    for k in 0..setup.steps {
        reference.step().unwrap();
        for (_, flow, worst, end, reversed, diverged) in runs.iter_mut().filter(|r| r.5.is_none()) {
            if flow.step().is_err() {
                *diverged = Some(k + 1);
                continue;
            }
            if k % 10 != 0 && k + 1 != setup.steps {
                continue;
            }
            let [nx, ny, _] = flow.grid().cells;
            let mut now: f64 = 0.0;
            for j in 0..ny {
                for i in 0..nx {
                    let (_, a) = flow.moments(flow.grid().index([i, j, 0])).unwrap();
                    let (_, b) = reference.moments(reference.grid().index([i, j, 0])).unwrap();
                    if i + 5 < nx {
                        now = now.max(((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt() / u0);
                    }
                    if i == nx - 1 && a[0] < 0.0 {
                        *reversed += 1;
                    }
                }
            }
            *worst = worst.max(now);
            *end = now;
        }
    }
    runs.into_iter()
        .map(|(outflow, _, worst, end, reversed, diverged)| {
            match diverged {
                Some(k) => println!("vortex of swirl {swirl} U through the outlet, {outflow:?}: diverged at step {k} of {}", setup.steps),
                None => println!(
                    "vortex of swirl {swirl} U through the outlet, {outflow:?}: reflected velocity up to {worst:.3e} U, {end:.3e} U at the end ({reversed} reversed outlet cell samples)"
                ),
            }
            (outflow, worst, end, reversed, diverged)
        })
        .collect()
}

/// With the two-relaxation-time collision (the vortex runs leave `collision` at its default): the
/// 0.8 ratio was measured with it. Under the regularised collision the uniform outlet holds the
/// equilibrium too, which changes what it reflects, so this ratio is not that collision's bound.
#[test]
fn a_vortex_leaves_a_convective_outlet_with_less_reflection() {
    let runs = vortex_runs(0.5);
    let [(_, uniform, _, _, d0), (_, convective, end, _, d1)] = runs[..] else { panic!() };
    assert_eq!((d0, d1), (None, None));
    assert!(convective < 0.8 * uniform && end < 0.01, "{convective} vs {uniform}, {end}");
}

/// A vortex whose swirl (1.8 U) reverses the flow where it crosses the outlet: the uniform outlet
/// diverges, the convective one holds the reversed cells and stays stable.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn a_strong_vortex_with_backflow_leaves_a_convective_outlet_stably() {
    let runs = vortex_runs(1.8);
    let (_, worst, end, reversed, diverged) = runs[1];
    assert_eq!(diverged, None);
    assert!(reversed > 0, "the vortex reverses the flow at the outlet");
    assert!(worst.is_finite() && end < 0.2, "{worst} {end}");
}

/// The steady wake of the first case against the same channel 1.6 long, whose outlet is 13
/// diameters behind the cylinder: the velocity in the short channel's last layers should be the
/// long one's there. Started at the mean inlet speed and run to t = 8 (the steady wake).
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn a_wake_at_a_convective_outlet_matches_a_long_channel() {
    let (um, dx) = (0.3, 0.01);
    let run = |l: f64, outflow| channel(l, dx, um, outflow, fluid(8.0, 0.025 * dx / um, 1), 256, true).unwrap();
    let long = run(1.6, Outflow::Uniform);
    let mut errors = Vec::new();
    for outflow in BOTH {
        let r = run(0.4, outflow);
        let (a, b) = (r.frames.last().unwrap(), long.frames.last().unwrap());
        let [nx, ny, _] = r.grid.cells;
        // The largest difference of the velocity from the long channel's in layer i.
        let error = |i: usize| {
            (0..ny)
                .map(|j| {
                    let (u, v) = (a.velocity[r.grid.index([i, j, 0])], b.velocity[long.grid.index([i, j, 0])]);
                    f64::from((u[0] - v[0]).powi(2) + (u[1] - v[1]).powi(2)).sqrt()
                })
                .fold(0.0, f64::max)
        };
        let e = [error(nx - 3), error(nx - 2), error(nx - 1)];
        println!("wake against a long channel, {outflow:?}: |u − u_long| in the last three layers {:.2e} {:.2e} {:.2e} (of {um})", e[0], e[1], e[2]);
        errors.push(e);
    }
    let (uniform, convective) = (errors[0], errors[1]);
    assert!(convective.iter().all(|e| *e < 0.01), "{convective:?}");
    assert!(convective[2] < 0.5 * uniform[2], "{convective:?} vs {uniform:?}");
}

/// Schäfer–Turek's channel at Re 100 cut 1 diameter behind the cylinder, run 20 s (60 000 steps,
/// about 60 shedding periods) with reversed flow at the outlet: stable, the mass bounded.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn a_convective_outlet_stays_stable_with_backflow_for_long() {
    let (um, dx) = (1.5, 0.01);
    let r = channel(0.3, dx, um, Outflow::Convective, fluid(20.0, 0.05 * dx / um, 80), 256, false).unwrap();
    let [nx, ny, _] = r.grid.cells;
    let reversed = r.frames.iter().filter(|f| (0..ny).any(|j| f.velocity[r.grid.index([nx - 1, j, 0])][0] < 0.0)).count();
    let change = r.frames.iter().map(|f| (f.mass / r.frames[0].mass - 1.0).abs()).fold(0.0, f64::max);
    println!(
        "Re 100 short channel for 20 s, convective outlet: {reversed} of {} frames with reversed outlet cells, largest mass change {change:.2e}",
        r.frames.len()
    );
    assert_eq!(r.stats.steps, 60_000);
    assert!(reversed > 10, "{reversed}");
    assert!(change < 0.05, "{change}");
}

/// The mean, and the largest of the fit by the first three harmonics of the shedding frequency
/// `f`, of `x(t)` over `t` (a whole number of periods): the drag's and lift's periodic part
/// without the pressure waves at other frequencies.
fn harmonic_fit(t: &[f64], x: &[f64], f: f64) -> (f64, f64) {
    let n = x.len().max(1) as f64;
    let mean = x.iter().sum::<f64>() / n;
    let w = 2.0 * std::f64::consts::PI * f;
    let coefficients: Vec<(f64, f64)> = (1..=3)
        .map(|k| {
            let k = f64::from(k);
            let a = 2.0 / n * t.iter().zip(x).map(|(t, x)| (x - mean) * (k * w * t).cos()).sum::<f64>();
            let b = 2.0 / n * t.iter().zip(x).map(|(t, x)| (x - mean) * (k * w * t).sin()).sum::<f64>();
            (a, b)
        })
        .collect();
    let fit = |t: f64| {
        mean + coefficients.iter().zip(1..).map(|((a, b), k)| a * (f64::from(k) * w * t).cos() + b * (f64::from(k) * w * t).sin()).sum::<f64>()
    };
    (mean, t.iter().map(|t| fit(*t)).fold(f64::MIN, f64::max))
}

/// Schäfer–Turek's benchmark 2D-2 (M. Schäfer, S. Turek, NNFM 52:547–566, 1996: the cylinder of
/// the channel at Re 100, parabolic inlet with peak 1.5, mean 1): vortex shedding with c_D max
/// 3.22–3.24, c_L max 0.99–1.01, St 0.295–0.305. A convective outlet 7.5 diameters behind the
/// cylinder (l = 1) gives the forces and frequency of the benchmark's channel (l = 2.2) within 3 %;
/// 20 cells across the staircase cylinder, quasi-2D, started at the mean inlet speed, the last
/// 40 % of 8 s (9 periods) measured.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn schafer_turek_2d2_does_not_depend_on_the_outlet_distance() {
    let (um, dx, d, end) = (1.5, 0.005, 0.1, 8.0);
    let ubar = 2.0 / 3.0 * um;
    let mut results = Vec::new();
    for l in [1.0, 2.2] {
        let r = channel(l, dx, um, Outflow::Convective, fluid(end, 0.05 * dx / um, 1), 1024, true).unwrap();
        let side = r.forces.iter().find(|f| f.name == "side").unwrap();
        // Per unit depth, over the dynamic pressure on the diameter.
        let k = 2.0 / (ubar * ubar * d * dx);
        let from = r.history_times.iter().position(|t| *t >= 0.6 * end).unwrap();
        let (t, f) = (&r.history_times[from..], &side.force[from..]);
        let cd: Vec<f64> = f.iter().map(|f| f[0] * k).collect();
        let cl: Vec<f64> = f.iter().map(|f| f[1] * k).collect();
        // The lift's upward crossings of its mean: the shedding period, and a whole number of
        // periods to measure over.
        let mean = cl.iter().sum::<f64>() / cl.len() as f64;
        let ups: Vec<usize> = (1..cl.len()).filter(|i| cl[i - 1] < mean && cl[*i] >= mean).collect();
        let (first, last) = (ups[0], ups[ups.len() - 1]);
        let periods = (ups.len() - 1) as f64;
        let frequency = periods / (t[last] - t[first]);
        let (cd_mean, cd_max) = harmonic_fit(&t[first..last], &cd[first..last], frequency);
        let (_, cl_max) = harmonic_fit(&t[first..last], &cl[first..last], frequency);
        let st = frequency * d / ubar;
        println!("Schäfer–Turek 2D-2, outlet at {l}: c_D max {cd_max:.4} (mean {cd_mean:.4}), c_L max {cl_max:.4}, St {st:.4} ({periods} periods)");
        results.push([cd_max, cl_max, st]);
    }
    // Measured (M4 Max): l = 1: c_D max 3.356, c_L max 1.072, St 0.2989; l = 2.2: 3.446, 1.102,
    // 0.2984. The staircase cylinder of 20 cells puts c_D max 4–7 % and c_L max 7–10 % above the
    // benchmark (2D-1 at this resolution: c_D 3.4 % above).
    for (near, far) in results[0].iter().zip(&results[1]) {
        assert!((near / far - 1.0).abs() < 0.04, "{results:?}");
    }
    for [cd, cl, st] in results {
        assert!((cd / 3.23 - 1.0).abs() < 0.08 && (cl - 1.0).abs() < 0.12 && (st / 0.30 - 1.0).abs() < 0.02, "{cd} {cl} {st}");
    }
}

/// Schäfer–Turek 2D-2 on a coarser grid: 10 cells across the cylinder (Δx = 0.01, Δt = 0.05 Δx /
/// 1.5, τ = 0.51), 1024 facets, the outlet at l = 1, started at the mean inlet speed, 8 s. The
/// uniform outlet diverges when the shed vortices reach it (measured: step 6699; at l = 2.2 step
/// 21055; started from rest instead, steps 3889 and 5471), the convective one runs to the end and
/// sheds at the benchmark's frequency.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn a_coarse_cylinder_at_re_100_runs_with_a_convective_outlet_where_the_uniform_one_diverges() {
    let (um, dx, end) = (1.5, 0.01, 8.0);
    let run = |outflow| channel(1.0, dx, um, outflow, fluid(end, 0.05 * dx / um, 1), 1024, true);
    let uniform = run(Outflow::Uniform);
    println!("coarse 2D-2, uniform outlet: {:?}", uniform.as_ref().map(|r| r.stats.steps).map_err(|e| e.to_string()));
    assert!(matches!(uniform, Err(fieldcraft_lbm::Error::Unstable(_))));
    let r = run(Outflow::Convective).unwrap();
    let side = r.forces.iter().find(|f| f.name == "side").unwrap();
    let from = r.history_times.iter().position(|t| *t >= 0.6 * end).unwrap();
    let cl: Vec<f64> = side.force[from..].iter().map(|f| f[1]).collect();
    let mean = cl.iter().sum::<f64>() / cl.len() as f64;
    let ups: Vec<usize> = (1..cl.len()).filter(|i| cl[i - 1] < mean && cl[*i] >= mean).collect();
    let t = &r.history_times[from..];
    let st = (ups.len() - 1) as f64 / (t[ups[ups.len() - 1]] - t[ups[0]]) * 0.1 / (2.0 / 3.0 * um);
    println!("coarse 2D-2, convective outlet: {} steps, St {st:.4}", r.stats.steps);
    assert_eq!(r.stats.steps, 24_000);
    assert!((st / 0.30 - 1.0).abs() < 0.02, "{st}");
}
