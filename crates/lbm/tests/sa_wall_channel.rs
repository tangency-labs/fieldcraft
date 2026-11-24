//! The Spalart–Allmaras model with its own wall law (`wall_model: "spalart_allmaras_bounce"`) in a turbulent
//! channel at Re_τ 182.088 and Re_τ 1000.512, against (a) theory: the steady fully developed 1-D SA solution with
//! the law at its first node, solved here, and the exact identities of the lattice's steady column; and (b) Lee &
//! Moser (J. Fluid Mech. 774:395–415, 2015): their bulk and centre-line velocities, reported, not asserted.
//!
//! **The 1-D solutions.** In wall units, with χ = ν̃/ν, the momentum equation integrates to
//! `dU⁺/dy⁺ = (1 − y⁺/Re_τ) / (1 + χ f_v1)`, and SA's equation is `0 = P − D + (1/σ)[((1 + χ) χ')' + c_b2 χ'²]`
//! (Spalart & Allmaras, AIAA-92-0439, eqs. 7–9), with the S̃ limiter of Allmaras, Johnson & Spalart (2012); written here again from
//! the papers, only the constants come from `fieldcraft_lbm::sa`. Second-order finite volumes, pseudo-transient
//! Newton (switched evolution relaxation, Mulder & van Leer, J. Comput. Phys. 59:232–246, 1985), Richardson
//! extrapolation from m and 2m cells. The wall-modelled solution starts at the law's first node y₁ = Re_τ / N with
//! the law's stress share s and its ν̃ = κ s y₁ (Han, Ooka & Kikumoto 2021, eq. 14).
//!
//! **The lattice run.** A column 1 × (2N + 1) × 1 cells, δ = N cells, periodic in x and z, on-node walls at rows 0
//! and 2N, the body force u_τ²/δ, the regularized collision (TRT keeps the laminar floor 0.505 with this law and
//! cannot run these grids, #1063; HRR has a different steady identity and is not tested here), lattice bulk speed
//! 0.1. The walls take the whole body force of the column, the boundary nodes' half cells included (#1101), so the
//! law's share of the stress is 1. Each run is marched to steady state; then the identities of `check` are asserted, and `report` prints
//! the table of `docs/verification.md`.
//!
//! **The steady column's link identity (regularized collision).** Per row j (wall node at 0) let
//! `a_j = f(1,1,0) − f(−1,1,0)`, `b_j = f(1,−1,0) − f(−1,−1,0)`, `Π = a − b = Σ c_x c_y f` and `Q = a + b`. The
//! regularized update is `f* = f^eq + (1 − ω) f⁽¹⁾[P] + ½ F_q`. In the steady column (u = (u, 0, 0), F = (F, 0, 0),
//! ρ uniform, P_yy = 0) it gives `Π* = (1 − ω) Π` (the force's u_x F_y + u_y F_x is 0) and `Q* = ρu/3 + F/6` (the
//! equilibrium's ρu/3, half of the source's F/3, and a third-order part u_x P_yy + 2 u_y P_xy = 0: the odd moments
//! are rebuilt, not relaxed). Streaming gives `a_j = a*_{j−1}` and `b_j = b*_{j+1}`. The momentum flux through the
//! face j + ½ is `a*_j − b*_{j+1}`, and it falls by F per row; `Π_j + Π*_j` is the sum of the cell's two face fluxes,
//! so `(1 − ω_j/2)(−Π_j) = T_j = F (N − j)`: each cell's stress is the total stress, and the model's ν_t carries all of
//! it. Writing the face flux with `Q*` then gives, on every link j = 1 … 2N − 2,
//! `u_{j+1} − u_j = ½ (g_j + g_{j+1})`, `g_j = a (N − j) / (ν + ν_t,j)`: the trapezoid rule on the local shear, with no
//! further term. At row 1 eq. 14 gives `a_1 = b*_1 − τ_w` and the same form; the symmetric row N closes it with
//! `τ_w = F N` (the half cells' share included, #1101). The hybrid collision rebuilds from `σ P + (1 − σ) P_FD` and gets a five-row stencil that differs
//! by O(1 − σ), so it is not tested here.
//!
//! **Not in the trunk gate.** A settled lattice run is heavy (the five cases take about two minutes to reach 1e-6
//! u_τ): both lattice tests are `#[ignore]`. The trunk keeps the 1-D solutions; the
//! per-step law test is `the_sa_law_gives_the_walls_its_stress_and_the_first_node_its_nu_tilde` in `wall_bounce.rs`.
//! Two requested cases are absent: y⁺ 100 at Re_τ 180 needs fewer than 2 cells per δ, and Re_τ 2000 is not run:
//! its 1-D reference was not solved (from the plain start the iteration does not converge within its cap; the
//! scaled start of `resolved_start` was not tried there).
//!
//! Each test prints its numbers (`cargo nextest run -p fieldcraft-lbm --test sa_wall_channel --no-capture`).

mod common;

use common::*;
use fieldcraft_lbm::lattice::{C, Q};
use fieldcraft_lbm::sa::{C_B1, C_B2, C_T3, C_T4, C_V1, C_V2, C_V3, C_W2, C_W3, KAPPA, R_MAX, SIGMA, c_w1};
use fieldcraft_lbm::sa_law::u_plus;
use fieldcraft_lbm::{Flow, Setup};
use fieldcraft_model::{Collision, FaceCondition, Procedure, Turbulence, WallModel};

/// A Lee & Moser case: Re_τ, U_b⁺ = 1/u_τ and U_c⁺ where the repo cites it (headers and last rows of
/// `LM_Channel_0180/1000_mean_prof.dat`, the same constants as `crates/lbm-gpu/tests/verification.rs`), and the
/// reference grid's cell count (N must divide it).
struct Case {
    re_tau: f64,
    lm_bulk: f64,
    lm_centre: Option<f64>,
    reference_cells: usize,
}

const LM_180: Case = Case { re_tau: 182.088, lm_bulk: 1.0 / 6.37309e-2, lm_centre: None, reference_cells: 1728 };
const LM_1000: Case = Case { re_tau: 1000.512, lm_bulk: 1.0 / 5.00256e-2, lm_centre: Some(22.593), reference_cells: 3200 };
/// The lattice bulk speed the time step is set for (as `wall_bounce.rs`'s channels; the column's steady state
/// has no Mach-number error).
const LATTICE_SPEED: f64 = 0.1;
/// The longest run, in δ/u_τ: a cap (the change falls by about 0.91 per δ/u_τ, so 1e-6 u_τ takes about 120).
const LONGEST: u64 = 300;
/// The rows next to the wall whose modelled share of the stress is printed: the first cells, where wall-modelled
/// LES loses the resolved stress (Kawai & Larsson, Phys. Fluids 24:015105, 2012).
const SHARE_ROWS: usize = 3;
/// The decay ratio per δ/u_τ counts as settled when it changes by less than this share of itself (as `sa_channel.rs`).
const RATIO_SETTLED: f64 = 0.01;
/// The accepted band for the observed order of a second-order scheme (as `sa_mms.rs`).
const SECOND_ORDER: std::ops::RangeInclusive<f64> = 1.8..=2.2;
/// The most Newton iterations of the 1-D solution, and the update below which it has converged (as `sa_channel.rs`).
const NEWTON_CAP: usize = 200;
const NEWTON_STEP: f64 = 1e-10;
/// Simpson panels of the law's integral over [0, y₁] (reported only; the test prints the change from half as many).
const LAW_PANELS: usize = 4096;

// ---------------------------------------------------------------------------------------------
// The 1-D solutions.
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

/// SA's steady residual in each of the m cells of width `h` over [y0, re_tau] (wall units), with χ = `chi0` on the
/// lower face: the wall-resolved solution has y0 = chi0 = 0. Second-order fluxes; the lower face's flux and slope
/// are one-sided, as the wall-modelled solution's first cell needs.
fn residual(chi: &[f64], re_tau: f64, y0: f64, chi0: f64, h: f64) -> Vec<f64> {
    let m = chi.len();
    (0..m)
        .map(|j| {
            let v = chi[j];
            let y = y0 + (j as f64 + 0.5) * h;
            let omega = (1.0 - y / re_tau) / (1.0 + nu_t(v));
            let below = if j == 0 { (1.0 + chi0) * 2.0 * (v - chi0) / h } else { 0.5 * (2.0 + chi[j - 1] + v) * (v - chi[j - 1]) / h };
            let above = if j + 1 == m { 0.0 } else { 0.5 * (2.0 + v + chi[j + 1]) * (chi[j + 1] - v) / h };
            let slope = if j == 0 {
                (chi[1] + 3.0 * v - 4.0 * chi0) / (3.0 * h)
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

/// The 1-D solution on `m` cells over [y0, re_tau].
struct Reference {
    /// χ at the cell centres.
    chi: Vec<f64>,
    /// U⁺ at the m + 1 faces, from y0 (u0) to the centre line.
    u: Vec<f64>,
    iterations: usize,
}

/// The solution on `m` cells, Newton-started from `start` (χ at the cell centres) or, without one, from the profile
/// κ y (1 − y / 2 Re_τ). Prints the iteration count.
fn reference(re_tau: f64, m: usize, y0: f64, chi0: f64, u0: f64, start: Option<&[f64]>) -> Reference {
    let h = (re_tau - y0) / m as f64;
    let mut chi: Vec<f64> = match start {
        Some(s) => s.to_vec(),
        None => (0..m)
            .map(|j| {
                let y = (j as f64 + 0.5) * h;
                KAPPA * y * (1.0 - y / (2.0 * re_tau))
            })
            .collect(),
    };
    let mut res = residual(&chi, re_tau, y0, chi0, h);
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
            let r = residual(&bumped, re_tau, y0, chi0, h);
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
        let new = residual(&chi, re_tau, y0, chi0, h);
        if largest(&new) > 0.0 {
            pseudo_step *= largest(&res) / largest(&new);
        }
        res = new;
        if largest(&delta) <= NEWTON_STEP {
            break;
        }
    }
    let mut u = Vec::with_capacity(m + 1);
    u.push(u0);
    for (j, v) in chi.iter().enumerate() {
        let y = y0 + (j as f64 + 0.5) * h;
        let last = u[j];
        u.push(last + h * (1.0 - y / re_tau) / (1.0 + nu_t(*v)));
    }
    println!("1-D SA, {m} cells, y0 {y0:.3}, Re_τ {re_tau}: {iterations} iterations");
    Reference { chi, u, iterations }
}

/// The Re_τ of `LM_180`: the scale of the wall-resolved start.
const BASE_RE_TAU: f64 = 182.088;

/// The wall-resolved start at `re_tau` on `m` cells: the Re_τ 182.088 solution on the same cells, scaled by
/// Re_τ / 182.088 and capped at the inner law κ y⁺. The inner scaling keeps χ = κ y⁺ near the wall; in the outer
/// region the log-layer eddy viscosity ν_t = κ u_τ d scales with u_τ δ, so χ = ν_t/ν scales with Re_τ at fixed y/δ.
/// The start is the smaller of the two. (From the plain start the iteration needs 246 steps at Re_τ 1000.512 on
/// 1600 cells, above `NEWTON_CAP`; from this one 30.)
fn resolved_start(re_tau: f64, m: usize) -> Vec<f64> {
    let base = reference(BASE_RE_TAU, m, 0.0, 0.0, 0.0, None);
    let h = re_tau / m as f64;
    (0..m).map(|j| (KAPPA * (j as f64 + 0.5) * h).min(base.chi[j] * re_tau / BASE_RE_TAU)).collect()
}

/// The start of a wall-modelled grid of `m` cells from `y0` (spacing h = (Re_τ − y0)/m): the wall-resolved solution of
/// the same case on the same spacing, interpolated linearly to the grid's centres. The two differ by O(y0/δ) only.
fn modelled_start(re_tau: f64, m: usize, y0: f64) -> Vec<f64> {
    let h = (re_tau - y0) / m as f64;
    let donor_cells = (re_tau / h).round() as usize;
    let donor = reference(re_tau, donor_cells, 0.0, 0.0, 0.0, Some(&resolved_start(re_tau, donor_cells)));
    let donor_h = re_tau / donor_cells as f64;
    (0..m)
        .map(|j| {
            let p = (y0 + (j as f64 + 0.5) * h) / donor_h - 0.5;
            let i = (p.floor().max(0.0) as usize).min(donor_cells - 2);
            let t = (p - i as f64).clamp(0.0, 1.0);
            donor.chi[i] + t * (donor.chi[i + 1] - donor.chi[i])
        })
        .collect()
}

/// The Newton start of a solution on `m` cells from `y0`: the wall-resolved one, or the wall-modelled one.
fn start(re_tau: f64, m: usize, y0: f64) -> Vec<f64> {
    if y0 > 0.0 { modelled_start(re_tau, m, y0) } else { resolved_start(re_tau, m) }
}

/// The law's integral ∫₀^{y₁} s u⁺(s y) dy (composite Simpson, `LAW_PANELS` panels): the velocity the law gives over
/// the first cell when its stress share is `s` (u_τ = s in wall units).
fn law_integral(s: f64, y1: f64) -> f64 {
    let h = y1 / LAW_PANELS as f64;
    let f = |y: f64| s * u_plus(s * y);
    let odd: f64 = (1..LAW_PANELS).step_by(2).map(|k| f(k as f64 * h)).sum();
    let even: f64 = (2..LAW_PANELS).step_by(2).map(|k| f(k as f64 * h)).sum();
    h / 3.0 * (f(0.0) + 4.0 * odd + 2.0 * even + f(y1))
}

/// A 1-D solution extrapolated (Richardson, second order) from m and 2m cells, with the first cell's law integral
/// (0 for the wall-resolved solution).
struct Exact {
    /// U⁺ at the m + 1 coarse faces, from y0 to the centre line.
    u: Vec<f64>,
    /// The face spacing and the first face's height (wall units).
    h: f64,
    y0: f64,
    /// The law's integral over [0, y0] (0 when y0 = 0).
    law: f64,
    /// The largest Richardson correction: the reference's own error estimate.
    correction: f64,
}

fn exact(re_tau: f64, m: usize, y0: f64, chi0: f64, u0: f64, law: f64) -> Exact {
    let coarse = reference(re_tau, m, y0, chi0, u0, Some(&start(re_tau, m, y0)));
    let fine = reference(re_tau, 2 * m, y0, chi0, u0, Some(&start(re_tau, 2 * m, y0)));
    let u: Vec<f64> = (0..=m)
        .map(|i| {
            let (c, f) = (coarse.u[i], fine.u[2 * i]);
            f + (f - c) / 3.0
        })
        .collect();
    let correction = (0..=m).map(|i| (fine.u[2 * i] - coarse.u[i]).abs() / 3.0).fold(0.0, f64::max);
    Exact { u, h: (re_tau - y0) / m as f64, y0, law, correction }
}

impl Exact {
    /// U⁺ at height `y` (wall units), which must be a face of this solution (or its first node).
    fn at(&self, y: f64) -> f64 {
        let index = (y - self.y0) / self.h;
        assert!((index - index.round()).abs() < 1e-6, "y⁺ {y} is not a face of this solution");
        self.u[index.round() as usize]
    }

    /// U⁺ at the centre line.
    fn centre(&self) -> f64 {
        self.u[self.u.len() - 1]
    }

    /// The bulk velocity over the half channel, U_b⁺ = (law + trapezoid of the faces) / Re_τ.
    fn bulk(&self, re_tau: f64) -> f64 {
        let sum: f64 = self.u.iter().sum();
        (self.law + self.h * (sum - 0.5 * (self.u[0] + self.centre()))) / re_tau
    }
}

/// The wall-resolved 1-D solution of `case`.
fn wall_resolved(case: &Case) -> Exact {
    exact(case.re_tau, case.reference_cells, 0.0, 0.0, 0.0, 0.0)
}

/// The wall-modelled 1-D solution of `case` on `n` cells per δ with the law's stress share `s`: the first node is
/// y₁ = Re_τ / n; the grid is the wall-resolved one from y₁, so its faces are the lattice nodes.
fn wall_modelled(case: &Case, n: usize, s: f64) -> Exact {
    let (re, y1) = (case.re_tau, case.re_tau / n as f64);
    exact(re, case.reference_cells - case.reference_cells / n, y1, KAPPA * s * y1, s * u_plus(s * y1), law_integral(s, y1))
}

/// Solves on m/2, m and 2m cells, checks each is the turbulent solution, and returns the order of convergence of
/// the largest differences at the faces the grids share.
fn one_d_order(label: &str, re_tau: f64, m: usize, y0: f64, chi0: f64, u0: f64) -> f64 {
    let grids = [m / 2, m, 2 * m];
    let solutions: Vec<Reference> = grids.iter().map(|&k| reference(re_tau, k, y0, chi0, u0, Some(&start(re_tau, k, y0)))).collect();
    for (k, s) in grids.iter().zip(&solutions) {
        let centre = s.u[*k];
        println!("1-D SA, {label}, {k} cells: {} iterations, U_c⁺ {centre:.8}", s.iterations);
        assert!(centre < 0.5 * re_tau && s.chi.iter().all(|x| *x > 0.0), "{label}, {k} cells: not the turbulent solution");
    }
    let differences: Vec<f64> =
        solutions.windows(2).zip(&grids).map(|(w, &k)| (0..=k).map(|i| (w[1].u[2 * i] - w[0].u[i]).abs()).fold(0.0, f64::max)).collect();
    let order = libm::log2(differences[0] / differences[1]);
    println!("1-D SA, {label}: largest differences at shared faces {:.4e}, {:.4e}: order {order:.4}", differences[0], differences[1]);
    assert!(SECOND_ORDER.contains(&order), "{label}: order {order}");
    order
}

// ---------------------------------------------------------------------------------------------
// The lattice run.
// ---------------------------------------------------------------------------------------------

/// The walls' share of the body force, as a fraction of the law's stress (`s` in `wall_modelled`). The column's
/// body force is a·(M + ½Σρ) with M = 2N − 1 fluid cells and ½Σρ = 1 from the two nodes' half cells (#1101), so
/// a·2N, and the two walls carry half each: τ_w = a N = u_τ² (δ = 1, a = u_τ²/δ). The law's stress is the wall's
/// stress, so s = 1 at every N, and u_τw² = a N.
const WALL_SHARE: f64 = 1.0;

/// The lattice time step for `n` cells per δ.
fn time_step(case: &Case, n: usize) -> f64 {
    LATTICE_SPEED * (1.0 / n as f64) / case.lm_bulk
}

/// The column of `n` cells per δ: the model, δ = 1 and u_τ = 1 physically, ν = 1/Re_τ, the body force u_τ²/δ.
/// Returns the setup and the flow, started from the law's profile and ν̃ = κ u_τ y (1 − y/2) (wall units).
fn channel(case: &Case, n: usize) -> (Setup, Flow) {
    let dx = 1.0 / n as f64;
    let mut d = domain([0.0; 3], [dx, 2.0 + dx, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.acceleration = [1.0, 0.0, 0.0];
    let mut p = collide(LONGEST as f64, time_step(case, n), None, Collision::Regularized);
    if let Procedure::Fluid { turbulence, .. } = &mut p {
        *turbulence = Some(Turbulence::SpalartAllmaras { prandtl: 0.9, wall_model: Some(WallModel::SpalartAllmarasBounce), inflow_ratio: 3.0 });
    }
    let m = model(1.0, 1.0 / case.re_tau, d, p);
    let (setup, mut flow) = prepare(&m, &[]);
    let grid = *flow.grid();
    assert_eq!(grid.cells, [1, 2 * n + 1, 1], "the column is one cell along x and z");
    for w in &setup.warnings {
        println!("setup warning: {w}");
    }
    let ut = 1.0 / setup.velocity_scale();
    // The row's distance in cells from the nearer wall node: the walls are rows 0 and 2N.
    let from_wall = |j: usize| j.min(2 * n - j) as f64;
    let row_of = |cell: usize| ((grid.centre(cell)[1] - 0.5 * dx) / dx).round() as usize;
    flow.set_equilibrium(|cell| (1.0, [u_plus(case.re_tau * from_wall(row_of(cell)) / n as f64) * ut, 0.0, 0.0]));
    let mut nu_tilde = vec![0.0; grid.len()];
    for j in 0..=2 * n {
        let y = from_wall(j);
        nu_tilde[grid.index([0, j, 0])] = KAPPA * ut * y * (1.0 - y / (2.0 * n as f64));
    }
    flow.set_sa_nu_tilde(&nu_tilde).unwrap();
    (setup, flow)
}

/// A channel run to its steady state.
struct Settled {
    /// u⁺ at each row 0 … 2N (the walls' rows read 0: no velocity is stored there).
    u: Vec<f64>,
    /// ν̃ at each row (lattice units; 0 at the walls).
    nu_tilde: Vec<f64>,
    /// u_τ of the run (lattice), its ν (lattice), and the body acceleration (lattice).
    ut: f64,
    nu: f64,
    a: f64,
    /// The force the walls take in a step (lattice, x), the fluid mass and the two boundary nodes' half-cell mass
    /// (½ρ each, the interior density).
    wall: f64,
    mass: f64,
    half: f64,
    /// The friction velocity the walls carry, by row (lattice; 0 at the fluid rows).
    friction: Vec<f64>,
    /// The modelled share of the stress, rows 1 to 3 below the centre line: (1 − ω_j/2)(−Π_xy)/(ρ a (N − j)) from the
    /// stored populations, ω_j with the cell's eddy viscosity.
    share: Vec<f64>,
    /// The largest ν_t / ν.
    eddy_ratio: f64,
    /// The last change of u⁺ per δ/u_τ (u_τ), the distance left to steady state, the time in δ/u_τ and the steps.
    change: f64,
    remaining: f64,
    time: u64,
    steps: u64,
}

/// Run the column of `n` cells per δ until the distance left to its steady state is at most `tolerance` (in u_τ).
/// Per δ/u_τ the largest change of u⁺ is Δ_k; once its ratio ρ = Δ_k/Δ_{k−1} is below 1 and settled, what is left is
/// the geometric tail Δ_k ρ/(1 − ρ). Panics, with the last change and ratio, after `LONGEST` δ/u_τ.
fn settle(case: &Case, n: usize, tolerance: f64) -> Settled {
    let (setup, mut flow) = channel(case, n);
    let grid = *flow.grid();
    let ut = 1.0 / setup.velocity_scale();
    let nu = (setup.tau_plus - 0.5) / 3.0;
    let per_time = (1.0 / time_step(case, n)).round() as u64;
    let column: Vec<usize> = (0..=2 * n).map(|j| grid.index([0, j, 0])).collect();
    let profile = |flow: &Flow| -> Vec<f64> { column.iter().map(|&c| flow.moments(c).map_or(0.0, |(_, u)| u[0] / ut)).collect() };
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
                done = Some((time, change, remaining));
                break;
            }
        }
        (change_before, ratio_before) = (change, ratio);
    }
    let Some((time, change, remaining)) = done else {
        panic!(
            "N {n}: not steady after {LONGEST} δ/u_τ ({} steps): the last change per δ/u_τ is {change_before:.3e} u_τ, its ratio {ratio_before:.4}",
            LONGEST * per_time
        );
    };
    let field = flow.spalart_allmaras().unwrap().nu_tilde();
    let nu_tilde: Vec<f64> = column.iter().map(|&c| field[c]).collect();
    let mut friction = vec![0.0; 2 * n + 1];
    for w in flow.wall_nodes().unwrap().nodes() {
        if let Some(j) = column.iter().position(|&c| c == w.cell as usize) {
            friction[j] = w.friction;
        }
    }
    let a = flow.acceleration()[0];
    // Rows 1 to 3 below the centre line (whose stress is 0).
    let share: Vec<f64> = (1..=SHARE_ROWS.min(n - 1))
        .map(|j| {
            let cell = column[j];
            // Each cell relaxes with its own rate, its eddy viscosity included: τ⁺ = ½ + 3 (ν + ν_t).
            let omega = 1.0 / (0.5 + 3.0 * (nu + nu * nu_t(nu_tilde[j] / nu)));
            let pi_xy: f64 = (0..Q).map(|q| f64::from(C[q][0]) * f64::from(C[q][1]) * flow.population(q, cell)).sum();
            let rho = flow.moments(cell).map_or(f64::NAN, |(r, _)| r);
            (1.0 - 0.5 * omega) * (-pi_xy) / (rho * a * (n - j) as f64)
        })
        .collect();
    Settled {
        u: last,
        nu_tilde,
        ut,
        nu,
        a,
        wall: flow.forces().iter().map(|f| f[0]).sum(),
        mass: flow.mass(),
        half: flow.wall_nodes().unwrap().nodes().iter().map(|w| 0.5 * flow.moments(w.interior as usize).map_or(0.0, |m| m.0)).sum(),
        friction,
        share,
        eddy_ratio: flow.eddy_viscosity_ratio(),
        change,
        remaining,
        time,
        steps: time * per_time,
    }
}

/// The model-consistent bulk of a settled run over the half channel: the law's integral over [0, y₁] plus the
/// trapezoid over the nodes 1 … N (spacing y₁), divided by Re_τ.
fn lattice_bulk(run: &Settled, n: usize, y1: f64, re_tau: f64, s: f64) -> f64 {
    let nodes: f64 = (1..=n).map(|j| run.u[j]).sum::<f64>() - 0.5 * (run.u[1] + run.u[n]);
    (law_integral(s, y1) + y1 * nodes) / re_tau
}

/// The largest |U⁺_j − U⁺(j y₁)| over the nodes 1 … N of a run.
fn nodes_linf(run: &Settled, n: usize, y1: f64, reference: &Exact) -> f64 {
    (1..=n).map(|j| (run.u[j] - reference.at(j as f64 * y1)).abs()).fold(0.0, f64::max)
}

/// The identities of a settled run (the bars are the derivations in the comments; Δ is the run's last change).
fn check(case: &Case, n: usize, run: &Settled) {
    let delta = run.change;
    let eps = f64::EPSILON;
    let (ut, nu, a) = (run.ut, run.nu, run.a);
    let ne = |j: usize| nu + fieldcraft_lbm::sa::nu_t(run.nu_tilde[j], nu);

    // 1. The turbulent branch, ν̃ ≥ 0 on the fluid rows, and the symmetry of the profile.
    assert!(run.eddy_ratio > 0.0 && run.u[n] < 0.5 * case.re_tau, "N {n}: laminar (U_c⁺ {}, ν_t/ν {})", run.u[n], run.eddy_ratio);
    assert!((1..2 * n).all(|j| run.nu_tilde[j] >= 0.0), "N {n}: ν̃ below 0 on a fluid row");
    let asymmetry = (0..=n).map(|j| (run.u[j] - run.u[2 * n - j]).abs()).fold(0.0, f64::max);
    assert!(asymmetry <= delta, "N {n}: the profile is asymmetric by {asymmetry:.3e} against Δ {delta:.3e}");

    // 2. Momentum: the walls carry the body force on the fluid and on the boundary nodes' half cells, a (M + ½Σρ)
    // (#1101: the node's share, per step, is the wall's own; the fluid's budget is a M + Σ½ρa − F). The change of the
    // momentum over the last δ/u_τ moves at most M Δ u_τ, against a (M + ½Σρ) = (M + ½Σρ) u_τ²/δ per time N/u_τ:
    // a relative error of Δ. (Without the half cells the walls would carry a M, 1/(2N−1) short.)
    let balance = run.wall / (a * (run.mass + run.half)) - 1.0;
    println!("N {n}: momentum balance {balance:+.3e} (bar Δ {delta:.3e})");
    assert!(balance.abs() <= delta, "N {n}: momentum balance {balance:e}");

    // 3. The walls' friction velocity: τ_w = a N per wall (item 2), so u_τ = √(a N) = ut (δ = 1, a = u_τ²/δ): the
    // law's share of the stress is WALL_SHARE = 1 (the half cells' share is the walls' own, #1101).
    let share = WALL_SHARE.sqrt();
    for w in [0, 2 * n] {
        let kept = run.friction[w] / ut / share - 1.0;
        assert!(kept.abs() <= delta, "N {n}: row {w} friction off by {kept:e}");
    }

    // 4. The first node: U₁ = u_τ u⁺(u_τ y₁ / ν) with y₁ = 1 cell, and ν̃₁ = κ u_τ · 1 exactly.
    for (row, wall) in [(1, 0), (2 * n - 1, 2 * n)] {
        let f = run.friction[wall];
        let law = run.u[row] * ut / (f * u_plus(f / nu)) - 1.0;
        assert!(law.abs() <= delta, "N {n}: row {row} off the law by {law:e}");
        assert_eq!(run.nu_tilde[row], KAPPA * f * 1.0, "N {n}: row {row} ν̃");
    }

    // 5. The trapezoid identity on every link j = 1 … 2N − 2, in lattice velocity. The regularized update gives the
    // cell's stress as the total stress a ρ (N − j) and each face's flux as a ρ (N − j − ½), so
    // u_{j+1} − u_j = ½ (g_j + g_{j+1}) with g_j = a (N − j) / (ν + ν_t,j) (Spalart & Allmaras's ν_t as the collision
    // carries it; the derivation is in the module docs). The unsteady flux through a face is at most
    // Σ |du/dt| ≤ Δ u_τ², over the smallest ν + ν_t, row 1's.
    let g = |j: usize| a * (n as f64 - j as f64) / ne(j);
    let uu = |j: usize| run.u[j] * ut;
    let bar = delta * ut * ut / ne(1) + 64.0 * eps * uu(n);
    let mut worst = 0.0f64;
    for j in 1..2 * n - 1 {
        let residual = (uu(j + 1) - uu(j)) - 0.5 * (g(j) + g(j + 1));
        worst = worst.max(residual.abs());
    }
    println!("N {n}: link identity, largest residual {:.3e} u_τ against the bar {:.3e} u_τ", worst / ut, bar / ut);
    assert!(worst <= bar, "N {n}: link identity residual {worst:e} against the bar {bar:e}");

    // The LES check, for RANS: the modelled share of the stress in rows 1 to 3 (the identity makes it 1), and the
    // finite-difference reading (not the stress).
    for j in 1..=run.share.len() {
        let fd = if j == 1 { ne(1) * (uu(2) - uu(1)) / (a * (n - 1) as f64) } else { ne(j) * (uu(j + 1) - uu(j - 1)) / 2.0 / (a * (n - j) as f64) };
        println!("N {n}: row {j}: modelled share {:.6}; finite-difference reading {fd:.4} (not the stress)", run.share[j - 1]);
    }
}

/// What a settled run measures, against the wall-resolved solution, the law alone and the law with the stress as run.
fn report(case: &Case, n: usize, run: &Settled, wr: &Exact, law: &Exact, as_run: &Exact) {
    let re = case.re_tau;
    let y1 = re / n as f64;
    let s = WALL_SHARE.sqrt();
    let bulk = lattice_bulk(run, n, y1, re, s);
    let centre = run.u[n];
    let pct = |a: f64, b: f64| 100.0 * (a / b - 1.0);
    // C_f = 2/U_b²: C_f(a)/C_f(b) − 1.
    let cf_pct = |a: f64, b: f64| 100.0 * ((b / a).powi(2) - 1.0);
    let (wr_bulk, law_bulk, as_bulk) = (wr.bulk(re), law.bulk(re), as_run.bulk(re));
    println!(
        "Re_τ {re}, N {n} (y₁⁺ {y1:.3}): settled after {} δ/u_τ ({} steps), Δ {:.3e} u_τ, {:.2e} u_τ left",
        run.time, run.steps, run.change, run.remaining
    );
    println!("  lattice: U_b⁺ {bulk:.5}, U_c⁺ {centre:.5}, C_f {:.5e}", 2.0 / (bulk * bulk));
    println!("  1-D wall-resolved (Richardson, correction {:.1e}): U_b⁺ {wr_bulk:.5}, U_c⁺ {:.5}", wr.correction, wr.centre());
    println!(
        "  against wall-resolved: U_b {:+.3} %, U_c {:+.3} %, C_f {:+.2} %, L∞ over nodes {:.3e} u_τ",
        pct(bulk, wr_bulk),
        pct(centre, wr.centre()),
        cf_pct(bulk, wr_bulk),
        nodes_linf(run, n, y1, wr)
    );
    println!(
        "  law alone (s = 1): U_b {:+.3} %, C_f {:+.2} %, L∞ over nodes {:.3e} u_τ (against wall-resolved)",
        pct(law_bulk, wr_bulk),
        cf_pct(law_bulk, wr_bulk),
        (1..=n).map(|j| (law.at(j as f64 * y1) - wr.at(j as f64 * y1)).abs()).fold(0.0, f64::max)
    );
    println!(
        "  law with the wall stress as run (s = {s:.6}): U_b {:+.3} %, C_f {:+.2} %, L∞ over nodes {:.3e} u_τ (against wall-resolved)",
        pct(as_bulk, wr_bulk),
        cf_pct(as_bulk, wr_bulk),
        (1..=n).map(|j| (as_run.at(j as f64 * y1) - wr.at(j as f64 * y1)).abs()).fold(0.0, f64::max)
    );
    println!(
        "  lattice − law with the stress as run, L∞ over nodes {:.3e} u_τ (reported: the trapezoid rule's offset on the log profile)",
        nodes_linf(run, n, y1, as_run)
    );
    match case.lm_centre {
        Some(c) => println!(
            "  against Lee & Moser: U_b {:+.3} %, C_f {:+.2} %, U_c {:+.3} %",
            pct(bulk, case.lm_bulk),
            cf_pct(bulk, case.lm_bulk),
            pct(centre, c)
        ),
        None => {
            println!("  against Lee & Moser: U_b {:+.3} %, C_f {:+.2} % (no centre line given)", pct(bulk, case.lm_bulk), cf_pct(bulk, case.lm_bulk))
        }
    }
}

/// The three 1-D solutions of a case at `n` cells per δ: the wall-resolved one, the law alone (s = 1) and the law with
/// the stress share the lattice carries.
fn solutions(case: &Case, n: usize) -> (Exact, Exact, Exact) {
    assert_eq!(case.reference_cells % n, 0, "N {n} shares no faces with the reference grid of {}", case.reference_cells);
    (wall_resolved(case), wall_modelled(case, n, 1.0), wall_modelled(case, n, WALL_SHARE.sqrt()))
}

// ---------------------------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------------------------

/// The 1-D solutions are the turbulent ones and converge at second order at shared faces: the wall-resolved one on
/// m/2, m, 2m cells, and the wall-modelled one (N = 6, the law's share s = 1) on the same spacings.
#[test]
fn the_1d_solutions_are_turbulent_and_second_order() {
    let case = &LM_180;
    let n = 6;
    let m = case.reference_cells;
    one_d_order("wall-resolved Re_τ 182.088", case.re_tau, m, 0.0, 0.0, 0.0);
    let (re, y1) = (case.re_tau, case.re_tau / n as f64);
    one_d_order("wall law N 6 at Re_τ 182.088", re, m - m / n, y1, KAPPA * y1, u_plus(y1));
}

/// The wall law on the lattice at N = 6 (Re_τ 182.088) and N = 10 (Re_τ 1000.512), run until less than 1e-6 u_τ
/// is left: the identities of `check` hold, and the profile is reported against the 1-D solutions and Lee & Moser.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_wall_law_channel_settles_on_the_laws_stress() {
    for (case, n) in [(&LM_180, 6usize), (&LM_1000, 10)] {
        let run = settle(case, n, 1e-6);
        check(case, n, &run);
        let (wr, law, as_run) = solutions(case, n);
        report(case, n, &run, &wr, &law, &as_run);
    }
}

/// The wall law on the lattice against the wall-resolved solution: Re_τ 182.088 at N = 3, and Re_τ 1000.512 at N =
/// 16 and 32, each settled to 1e-6 u_τ and checked and reported as above; the 1-D solutions of Re_τ 1000.512 with the
/// law at N = 10 are checked for the second order first.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_wall_law_channel_against_the_wall_resolved_solution() {
    let case = &LM_1000;
    let (re, n) = (case.re_tau, 10usize);
    one_d_order("wall-resolved Re_τ 1000.512", re, case.reference_cells, 0.0, 0.0, 0.0);
    let y1 = re / n as f64;
    one_d_order("wall law N 10 at Re_τ 1000.512", re, case.reference_cells - case.reference_cells / n, y1, KAPPA * y1, u_plus(y1));
    for (case, n) in [(&LM_180, 3usize), (&LM_1000, 16), (&LM_1000, 32)] {
        let run = settle(case, n, 1e-6);
        check(case, n, &run);
        let (wr, law, as_run) = solutions(case, n);
        report(case, n, &run, &wr, &law, &as_run);
    }
}
