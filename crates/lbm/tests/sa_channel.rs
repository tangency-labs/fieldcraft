//! The Spalart–Allmaras model in a wall-resolved turbulent channel at Re_τ 178.12 (the "Re_τ 180" case), no wall
//! law, against (a) theory: the steady fully developed 1-D SA solution with the same constants, solved here; and
//! (b) a benchmark: Moser, Kim & Mansour's DNS (Phys. Fluids 11(4):943–945, 1999, `chan180.means`), printed and
//! not asserted: the gap is the model's own error (no published SA-against-DNS figure at Re_τ 180 gives a bar).
//!
//! **The 1-D solution.** In wall units, with χ = ν̃/ν, the momentum equation integrates to
//! `dU⁺/dy⁺ = (1 − y⁺/Re_τ) / (1 + χ f_v1)` and SA's equation is
//! `0 = P − D + (1/σ)[((1 + χ) χ')' + c_b2 χ'²]`, χ = 0 at the wall and χ' = 0 at the centre line, with the
//! closure of Spalart & Allmaras (1994) and the S̃ limiter of Allmaras, Johnson & Spalart (2012), written here
//! again from the papers (only the constants come from `fieldcraft_lbm::sa`). Second-order finite volumes,
//! pseudo-transient Newton (switched evolution relaxation, Mulder & van Leer, J. Comput. Phys. 59:232–246, 1985),
//! Richardson extrapolation from 1728 and 3456 cells.
//!
//! **The lattice run.** A box 1 × N × 1 cells, periodic in x and z, halfway walls, the body force u_τ²/δ, the
//! fluid's lattice viscosity 0.002 on every grid (τ⁺ = 0.506, just above the solver's floor 0.505), so Δt ∝ Δx².
//! That scaling is what makes the scheme second order here: the point-implicit sources of `SaField` leave a
//! steady splitting error of relative size Δt c_b1 S̃, which is ∝ Δx² only with a fixed lattice viscosity. The
//! flow has no gradient of ν̃ along the velocity, so the first-order upwind advection does not enter; the
//! diffusion is second order (`sa_mms.rs`) and so is the halfway wall. The accepted band for the observed order
//! is `sa_mms.rs`'s for a second-order scheme, [1.8, 2.2].
//!
//! Each test prints its numbers (`cargo test -p fieldcraft-lbm --test sa_channel -- --nocapture`).

mod common;

use common::*;
use fieldcraft_lbm::Flow;
use fieldcraft_lbm::sa::{C_B1, C_B2, C_T3, C_T4, C_V1, C_V2, C_V3, C_W2, C_W3, KAPPA, R_MAX, SIGMA, c_w1};
use fieldcraft_model::{Collision, FaceCondition, Procedure, Turbulence};

/// The friction Reynolds number u_τ δ / ν of Moser, Kim & Mansour's "Re_τ 180" data (the header of
/// `chan180.means`).
const RE_TAU: f64 = 178.12;
/// Their centre-line velocity over u_τ (the last row of `chan180.means`, y = 1).
const MKM_CENTRE: f64 = 18.301;
/// Their bulk velocity over u_τ (the trapezoid integral of `Umean` in `chan180.means`, 15.679).
const MKM_BULK: f64 = 15.68;

/// The fluid's lattice viscosity on every grid: τ⁺ = ½ + 3ν = 0.506, just above the solver's floor 0.505
/// (`setup::MIN_TAU`), which gives the lowest lattice speed that floor allows.
const LATTICE_NU: f64 = 0.002;
/// The longest run, in δ/u_τ.
const LONGEST: u64 = 150;
/// The decay ratio per δ/u_τ counts as settled when it changes by less than this share of itself: the
/// geometric tail it gives is then the distance to the steady state.
const RATIO_SETTLED: f64 = 0.01;
/// The accepted band for the observed order of a second-order scheme (as `sa_mms.rs`).
const SECOND_ORDER: std::ops::RangeInclusive<f64> = 1.8..=2.2;
/// The reference grids: every lattice grid's cell centres lie on faces of both.
const REFERENCE_CELLS: usize = 1728;
/// The most Newton iterations of the reference, and the update below which it has converged.
const NEWTON_CAP: usize = 200;
const NEWTON_STEP: f64 = 1e-10;

// ---------------------------------------------------------------------------------------------
// The 1-D reference.
// ---------------------------------------------------------------------------------------------

fn f_v1(chi: f64) -> f64 {
    let c3 = chi * chi * chi;
    c3 / (c3 + C_V1 * C_V1 * C_V1)
}

/// ν_t / ν = χ f_v1 (χ ≥ 0 in the reference).
fn nu_t(chi: f64) -> f64 {
    chi * f_v1(chi)
}

/// P − D in wall units at wall distance `y` with vorticity `omega`, for χ ≥ 0.
fn net_source(chi: f64, omega: f64, y: f64) -> f64 {
    let f_t2 = C_T3 * libm::exp(-C_T4 * chi * chi);
    let f_v2 = 1.0 - chi / (1.0 + chi * f_v1(chi));
    let k2y2 = KAPPA * KAPPA * y * y;
    let s_bar = chi * f_v2 / k2y2;
    let s_tilde = if s_bar >= -C_V2 * omega {
        omega + s_bar
    } else {
        omega + omega * (C_V2 * C_V2 * omega + C_V3 * s_bar) / ((C_V3 - 2.0 * C_V2) * omega - s_bar)
    };
    let r = if s_tilde * k2y2 > 0.0 { (chi / (s_tilde * k2y2)).min(R_MAX) } else { R_MAX };
    let g = r + C_W2 * (r.powi(6) - r);
    let c6 = C_W3.powi(6);
    let f_w = g * libm::pow((1.0 + c6) / (g.powi(6) + c6), 1.0 / 6.0);
    C_B1 * (1.0 - f_t2) * s_tilde * chi - (c_w1() * f_w - C_B1 / (KAPPA * KAPPA) * f_t2) * (chi / y) * (chi / y)
}

/// The residual of SA's steady equation in each of the cells of width `h` (see the module docs).
fn residual(chi: &[f64], h: f64) -> Vec<f64> {
    let m = chi.len();
    (0..m)
        .map(|j| {
            let v = chi[j];
            let y = (j as f64 + 0.5) * h;
            let omega = (1.0 - y / RE_TAU) / (1.0 + nu_t(v));
            let below = if j == 0 { 2.0 * v / h } else { 0.5 * (2.0 + chi[j - 1] + v) * (v - chi[j - 1]) / h };
            let above = if j + 1 == m { 0.0 } else { 0.5 * (2.0 + v + chi[j + 1]) * (chi[j + 1] - v) / h };
            let slope = if j == 0 {
                (chi[1] + 3.0 * v) / (3.0 * h)
            } else if j + 1 == m {
                (v - chi[j - 1]) / (2.0 * h)
            } else {
                (chi[j + 1] - chi[j - 1]) / (2.0 * h)
            };
            net_source(v, omega, y) + ((above - below) / h + C_B2 * slope * slope) / SIGMA
        })
        .collect()
}

/// The solution of the tridiagonal system with sub-diagonal `a`, diagonal `b` and super-diagonal `c`
/// (Thomas algorithm).
fn tridiagonal(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = b.len();
    let (mut cp, mut dp) = (vec![0.0; n], vec![0.0; n]);
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    for i in 1..n {
        let pivot = b[i] - a[i] * cp[i - 1];
        cp[i] = c[i] / pivot;
        dp[i] = (d[i] - a[i] * dp[i - 1]) / pivot;
    }
    let mut x = vec![0.0; n];
    x[n - 1] = dp[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }
    x
}

fn largest(v: &[f64]) -> f64 {
    v.iter().fold(0.0f64, |m, x| m.max(x.abs()))
}

/// The 1-D solution on `m` cells over the half channel.
struct Reference {
    /// χ at the cell centres.
    chi: Vec<f64>,
    /// U⁺ at the m + 1 faces, from the wall (0) to the centre line.
    u: Vec<f64>,
    bulk: f64,
    iterations: usize,
}

fn reference(m: usize) -> Reference {
    let h = RE_TAU / m as f64;
    let mut chi: Vec<f64> = (0..m)
        .map(|j| {
            let y = (j as f64 + 0.5) * h;
            KAPPA * y * (1.0 - y / (2.0 * RE_TAU))
        })
        .collect();
    let mut res = residual(&chi, h);
    let mut pseudo_step = 1.0;
    let bump = f64::EPSILON.sqrt();
    let mut iterations = 0;
    loop {
        iterations += 1;
        assert!(iterations <= NEWTON_CAP, "the 1-D solution on {m} cells did not converge in {NEWTON_CAP} iterations: residual {:e}", largest(&res));
        // The Jacobian's three diagonals by forward differences, three colours.
        let (mut sub, mut diagonal, mut sup) = (vec![0.0; m], vec![0.0; m], vec![0.0; m]);
        for colour in 0..3 {
            let mut bumped = chi.clone();
            let mut by = vec![0.0; m];
            for j in (colour..m).step_by(3) {
                by[j] = bump * chi[j].abs().max(1.0);
                bumped[j] += by[j];
            }
            let r = residual(&bumped, h);
            for j in (colour..m).step_by(3) {
                if j > 0 {
                    sup[j - 1] = (r[j - 1] - res[j - 1]) / by[j];
                }
                diagonal[j] = (r[j] - res[j]) / by[j];
                if j + 1 < m {
                    sub[j + 1] = (r[j + 1] - res[j + 1]) / by[j];
                }
            }
        }
        // (I/Δτ − J) δ = residual.
        let a: Vec<f64> = sub.iter().map(|x| -x).collect();
        let b: Vec<f64> = diagonal.iter().map(|x| 1.0 / pseudo_step - x).collect();
        let c: Vec<f64> = sup.iter().map(|x| -x).collect();
        let delta = tridiagonal(&a, &b, &c, &res);
        for (x, d) in chi.iter_mut().zip(&delta) {
            *x += d;
        }
        assert!(chi.iter().all(|x| x.is_finite() && *x >= 0.0), "the 1-D solution left χ ≥ 0 at iteration {iterations}");
        let new = residual(&chi, h);
        if largest(&new) > 0.0 {
            pseudo_step *= largest(&res) / largest(&new);
        }
        res = new;
        if largest(&delta) <= NEWTON_STEP {
            break;
        }
    }
    let mut u = Vec::with_capacity(m + 1);
    u.push(0.0);
    for (j, v) in chi.iter().enumerate() {
        let y = (j as f64 + 0.5) * h;
        let last = u[j];
        u.push(last + h * (1.0 - y / RE_TAU) / (1.0 + nu_t(*v)));
    }
    let bulk = (u.iter().sum::<f64>() - 0.5 * (u[0] + u[m])) * h / RE_TAU;
    Reference { chi, u, bulk, iterations }
}

/// The reference extrapolated (Richardson, second order) from `REFERENCE_CELLS` and twice as many.
struct Exact {
    coarse: Reference,
    fine: Reference,
}

impl Exact {
    fn new() -> Self {
        Exact { coarse: reference(REFERENCE_CELLS), fine: reference(2 * REFERENCE_CELLS) }
    }

    /// U⁺ at face `i` of the coarse grid.
    fn u(&self, i: usize) -> f64 {
        let (c, f) = (self.coarse.u[i], self.fine.u[2 * i]);
        f + (f - c) / 3.0
    }

    fn bulk(&self) -> f64 {
        self.fine.bulk + (self.fine.bulk - self.coarse.bulk) / 3.0
    }

    fn centre(&self) -> f64 {
        self.u(REFERENCE_CELLS)
    }

    /// The largest Richardson correction: the reference's own error estimate.
    fn correction(&self) -> f64 {
        (0..=REFERENCE_CELLS).map(|i| (self.fine.u[2 * i] - self.coarse.u[i]).abs() / 3.0).fold(0.0, f64::max)
    }

    /// U⁺ at the centre of cell `j` of a lattice `n` cells across the whole channel: y⁺ = (2j' + 1) Re_τ / n
    /// from the nearer wall, a face of the coarse grid.
    fn at_cell(&self, j: usize, n: usize) -> f64 {
        assert_eq!(REFERENCE_CELLS % n, 0, "a lattice of {n} cells shares no faces with the reference");
        let from_wall = j.min(n - 1 - j);
        self.u((2 * from_wall + 1) * (REFERENCE_CELLS / n))
    }
}

// ---------------------------------------------------------------------------------------------
// The lattice run.
// ---------------------------------------------------------------------------------------------

/// The channel on `n` cells across: δ = 1, u_τ = 1, ν = 1/Re_τ (so a velocity is u⁺ and a time is in δ/u_τ),
/// starting from Moser, Kim & Mansour's bulk velocity and ν̃ = 3ν. Returns the model and its time step.
fn channel_model(n: usize) -> (fieldcraft_model::Model, f64) {
    let dx = 2.0 / n as f64;
    let nu = 1.0 / RE_TAU;
    let dt = LATTICE_NU * dx * dx / nu;
    let mut d = domain([0.0; 3], [dx, 2.0, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.acceleration = [1.0, 0.0, 0.0];
    d.initial_velocity = [MKM_BULK, 0.0, 0.0];
    let mut p = collide(LONGEST as f64, dt, None, Collision::Trt);
    if let Procedure::Fluid { turbulence, .. } = &mut p {
        *turbulence = Some(Turbulence::SpalartAllmaras { prandtl: 0.9, wall_model: None, inflow_ratio: 3.0 });
    }
    (model(1.0, nu, d, p), dt)
}

/// A channel run to its steady state.
struct Settled {
    /// u / u_τ at each cell across the channel.
    u: Vec<f64>,
    /// ν̃ / ν at each cell.
    chi: Vec<f64>,
    /// The flow time it took, in δ/u_τ, and the lattice steps.
    time: u64,
    steps: u64,
    /// The distance left to the steady state, in u_τ (the geometric tail of the last changes).
    remaining: f64,
    /// The largest ν_t / ν.
    eddy_ratio: f64,
}

/// Run the channel on `n` cells until the distance left to its steady state is at most `tolerance` (in u_τ).
/// Per δ/u_τ the largest change of u⁺ is Δ_k; once its ratio ρ = Δ_k/Δ_{k−1} is below 1 and settled, what is
/// left is the geometric tail Δ_k ρ/(1 − ρ). Fails, with the last change, after `LONGEST` δ/u_τ.
fn settle(n: usize, tolerance: f64) -> Settled {
    let (m, dt) = channel_model(n);
    let (_, mut flow) = prepare(&m, &[]);
    let grid = *flow.grid();
    assert_eq!(grid.cells, [1, n, 1], "the box is one cell along x and z");
    let u_tau = dt * n as f64 / 2.0; // u_τ Δt/Δx, lattice units
    let per_time = (1.0 / dt).round() as u64; // steps in one δ/u_τ
    let column: Vec<usize> = (0..n).map(|j| grid.index([0, j, 0])).collect();
    let profile = |flow: &Flow| -> Vec<f64> { column.iter().map(|&c| flow.moments(c).map_or(f64::NAN, |(_, u)| u[0] / u_tau)).collect() };
    let mut last = profile(&flow);
    let (mut change_before, mut ratio_before) = (f64::INFINITY, f64::INFINITY);
    let mut done = None;
    for time in 1..=LONGEST {
        for _ in 0..per_time {
            flow.step().unwrap();
        }
        let now = profile(&flow);
        let change = now.iter().zip(&last).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        let ratio = change / change_before;
        last = now;
        if ratio < 1.0 && (ratio - ratio_before).abs() <= RATIO_SETTLED * ratio {
            let remaining = change * ratio / (1.0 - ratio);
            if remaining <= tolerance {
                done = Some((time, remaining));
                break;
            }
        }
        (change_before, ratio_before) = (change, ratio);
    }
    let Some((time, remaining)) = done else {
        panic!("N {n}: not steady after {LONGEST} δ/u_τ ({} steps): the last change per δ/u_τ is {change_before:.3e} u_τ", LONGEST * per_time);
    };
    let field = flow.spalart_allmaras().unwrap().nu_tilde();
    let chi: Vec<f64> = column.iter().map(|&c| field[c] / LATTICE_NU).collect();
    Settled { u: last, chi, time, steps: time * per_time, remaining, eddy_ratio: flow.eddy_viscosity_ratio() }
}

/// What a settled run measures against the reference.
struct Measured {
    l2: f64,
    linf: f64,
    /// y⁺ of the cell with the largest error.
    worst: f64,
    bulk: f64,
    /// The centre-line velocity: u at the two cells next to it and the two beyond, extrapolated to the centre
    /// line by the even parabola through them, (9 u₁ − u₂)/8.
    centre: f64,
}

fn measure(n: usize, run: &Settled, exact: &Exact) -> Measured {
    let errors: Vec<f64> = run.u.iter().enumerate().map(|(j, u)| u - exact.at_cell(j, n)).collect();
    let l2 = (errors.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();
    let (worst, linf) = errors.iter().enumerate().fold((0, 0.0f64), |(k, m), (j, e)| if e.abs() > m { (j, e.abs()) } else { (k, m) });
    let bulk = run.u.iter().sum::<f64>() / n as f64;
    let centre = (9.0 * run.u[n / 2 - 1] - run.u[n / 2 - 2]) / 8.0;
    Measured { l2, linf, worst: (2 * worst.min(n - 1 - worst) + 1) as f64 * RE_TAU / n as f64, bulk, centre }
}

fn report(n: usize, run: &Settled, got: &Measured, exact: &Exact) {
    let symmetry = (0..n / 2).map(|j| (run.u[j] - run.u[n - 1 - j]).abs()).fold(0.0, f64::max);
    println!(
        "N {n} (first cell y⁺ {:.2}): steady after {} δ/u_τ ({} steps), {:.2e} u_τ left; largest ν_t/ν {:.3}; χ at the centre {:.4} (1-D {:.4})",
        RE_TAU / n as f64,
        run.time,
        run.steps,
        run.remaining,
        run.eddy_ratio,
        run.chi[n / 2],
        exact.fine.chi[2 * REFERENCE_CELLS - 1]
    );
    println!(
        "  against the 1-D solution: L2 {:.4e}, L∞ {:.4e} u_τ at y⁺ {:.1}; U_b⁺ {:.5} (1-D {:.5}), U_c⁺ {:.5} (1-D {:.5}); bottom/top asymmetry {symmetry:.1e}",
        got.l2,
        got.linf,
        got.worst,
        got.bulk,
        exact.bulk(),
        got.centre,
        exact.centre()
    );
    println!(
        "  against Moser, Kim & Mansour: U_b⁺ {:+.4} ({:+.2} %), U_c⁺ {:+.4} ({:+.2} %), C_f {:.4e} against {:.4e} ({:+.2} %)",
        got.bulk - MKM_BULK,
        100.0 * (got.bulk / MKM_BULK - 1.0),
        got.centre - MKM_CENTRE,
        100.0 * (got.centre / MKM_CENTRE - 1.0),
        2.0 / (got.bulk * got.bulk),
        2.0 / (MKM_BULK * MKM_BULK),
        100.0 * ((MKM_BULK / got.bulk).powi(2) - 1.0)
    );
}

/// What holds on every grid: the turbulent branch (laminar flow has U_c⁺ = Re_τ/2 and ν_t = 0) and ν̃ ≥ 0
/// everywhere (the reference has no negative branch). The benchmark (b) is printed by `report`, not asserted:
/// the gap is the model's error, reported.
fn check(n: usize, run: &Settled, got: &Measured) {
    assert!(run.eddy_ratio > 0.0 && got.centre < 0.5 * RE_TAU, "N {n}: the flow is laminar (U_c⁺ {}, ν_t/ν {})", got.centre, run.eddy_ratio);
    assert!(run.chi.iter().all(|x| *x >= 0.0), "N {n}: ν̃ went below 0");
}

// ---------------------------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------------------------

/// The reference itself: the Newton iteration converges, the solution is the turbulent one (laminar flow has
/// U_c⁺ = Re_τ/2), and U⁺ at shared faces converges at second order over 864, 1728 and 3456 cells.
#[test]
fn the_reference_is_turbulent_and_second_order() {
    let grids = [REFERENCE_CELLS / 2, REFERENCE_CELLS, 2 * REFERENCE_CELLS];
    let solutions: Vec<Reference> = grids.iter().map(|&m| reference(m)).collect();
    for (m, s) in grids.iter().zip(&solutions) {
        let centre = s.u[*m];
        println!(
            "1-D SA, {m} cells: {} iterations, U_b⁺ {:.8}, U_c⁺ {:.8}, C_f {:.6e}, largest χ {:.5}",
            s.iterations,
            s.bulk,
            centre,
            2.0 / (s.bulk * s.bulk),
            s.chi.iter().copied().fold(0.0, f64::max)
        );
        assert!(centre < 0.5 * RE_TAU && s.chi.iter().all(|x| *x > 0.0), "{m} cells: not the turbulent solution");
    }
    let differences: Vec<f64> =
        solutions.windows(2).zip(&grids).map(|(w, &m)| (0..=m).map(|i| (w[1].u[2 * i] - w[0].u[i]).abs()).fold(0.0, f64::max)).collect();
    let order = libm::log2(differences[0] / differences[1]);
    println!("1-D SA: largest differences at shared faces {:.4e}, {:.4e}: order {order:.4}", differences[0], differences[1]);
    assert!(SECOND_ORDER.contains(&order), "{order}");
}

/// The channel on the coarsest grid (first cell at y⁺ 1.86), run until less than 0.01 u_τ is left: it settles
/// on the turbulent branch with ν̃ ≥ 0. Prints its distance from the 1-D solution and from the DNS.
#[test]
fn the_channel_settles_on_the_turbulent_branch() {
    let n = 96;
    let exact = Exact::new();
    let run = settle(n, 0.01);
    let got = measure(n, &run, &exact);
    report(n, &run, &got, &exact);
    check(n, &run, &got);
}

/// The refinement study: N = 96, 144, 216 (first cell at y⁺ 1.86, 1.24, 0.82), each run until less than
/// 1e-4 u_τ is left. The error of u⁺ against the 1-D solution falls at second order on the finest pair, in L2
/// and L∞, and falls on every refinement. The iteration and the reference are both converged well inside the
/// error they are compared with: what is left of the transient, and the reference's Richardson correction,
/// are each below 1 % of the lattice error.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_channel_converges_to_the_1d_solution_at_second_order() {
    let grids = [96usize, 144, 216];
    let exact = Exact::new();
    println!("1-D SA (Richardson): U_b⁺ {:.6}, U_c⁺ {:.6}; its largest correction {:.2e} u_τ", exact.bulk(), exact.centre(), exact.correction());
    let mut results = Vec::new();
    for &n in &grids {
        let run = settle(n, 1e-4);
        let got = measure(n, &run, &exact);
        report(n, &run, &got, &exact);
        check(n, &run, &got);
        assert!(run.remaining <= 0.01 * got.linf, "N {n}: {:.2e} u_τ left of the transient against an error of {:.2e}", run.remaining, got.linf);
        assert!(
            exact.correction() <= 0.01 * got.linf,
            "N {n}: the reference's correction {:.2e} against an error of {:.2e}",
            exact.correction(),
            got.linf
        );
        results.push(got);
    }
    let orders: Vec<(f64, f64)> = results
        .windows(2)
        .zip(grids.windows(2))
        .map(|(e, n)| {
            let refinement = libm::log(n[1] as f64 / n[0] as f64);
            (libm::log(e[0].l2 / e[1].l2) / refinement, libm::log(e[0].linf / e[1].linf) / refinement)
        })
        .collect();
    println!("observed orders (L2, L∞): {orders:.3?}");
    for (e, n) in results.windows(2).zip(grids.windows(2)) {
        assert!(e[1].l2 < e[0].l2 && e[1].linf < e[0].linf, "the error does not fall from N {} to N {}", n[0], n[1]);
    }
    let (l2, linf) = orders[orders.len() - 1];
    assert!(SECOND_ORDER.contains(&l2) && SECOND_ORDER.contains(&linf), "{orders:?}");
}
