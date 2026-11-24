//! The Spalart–Allmaras RANS channel against Lee & Moser's DNS at Re_τ 182, 543 and 1000, and the
//! laminar channel that checks its measuring pipeline.
//!
//! Fully developed plane channel flow is one-dimensional in the mean, so a RANS run is a column:
//! one cell periodic in x and z, walls at y = 0 and 2δ, driven by a body force a = u_τ²/δ.
//! Units are u_τ = δ = 1 and ν = 1/Re_τ.
//!
//! What a run reads, each a time mean over the averaging window:
//! - **carried** wall stress: the momentum the flow hands to the walls, from the closed-box
//!   budget of `support/wall_balance.rs`, the body impulse less the change of the
//!   populations' momentum, per wall area;
//! - **applied**: the momentum exchange the walls report (`Flow::forces`), per area;
//! - **requested**: what the wall model asked for, ρ u_τ² per wall face: with SA, the wall
//!   function's u_τ at each first fluid node (`Flow::sa_wall_function`, Spalding's law on Han's
//!   on-node bounce, Musker's under the Musker wall); with a subgrid model under Musker's wall,
//!   Musker's request through `wall_balance::musker_request`;
//! - U_b, the bulk velocity, and C_f = τ_w / (½ ρ U_b²) from the carried stress. On on-node walls
//!   the interval between each wall and its first fluid node is integrated from the wall law
//!   (common::cases::WallInterval), with the trapezoid from u = 0 reported beside it.
//!
//! Settling ([`GATED`], common::cases): a flow-through is 2πδ/U_b, the length of the usual channel
//! box at the bulk speed. Two checks must both hold on means over the averaging window:
//! - the carried C_f moves by less than `drift` from the window before;
//! - the applied wall drag balances the body force to within `drift`, so a slow, steady approach to
//!   the fully developed state can't pass for settled.
//!
//! The averaging window is a flow-through, and at least two periods of the column's acoustic mode
//! between its walls (2 N_y / c_s steps).
//!
//! The reference is Lee & Moser's DNS at Re_τ = 1000.512: U_b⁺ = 1/0.0500256 = 19.990 and
//! C_f = 2 (u_τ/U_b)² = 5.005 × 10⁻³. The values come from the bulk velocity and friction velocity
//! in the header of their published mean profile, LM_Channel_1000_mean_prof.dat
//! (<https://turbulence.oden.utexas.edu/channel2015/>; M. Lee and R. D. Moser, *Direct numerical
//! simulation of turbulent channel flow up to Re_τ ≈ 5200*, J. Fluid Mech. 774:395–415, 2015,
//! doi:10.1017/jfm.2015.268).

mod common;
#[path = "support/sa_1d.rs"]
#[allow(dead_code)]
mod sa_1d;
#[path = "support/wall_balance.rs"]
mod wall_balance;

use common::cases::*;
use common::*;
use fieldcraft_lbm::Flow;
use fieldcraft_lbm::statistics::Blocks;
use fieldcraft_model::{Collision, FaceCondition, Procedure, Turbulence, WallModel};

/// Lee & Moser (2015), Re_τ = 1000.512: U_b⁺ and C_f on the bulk velocity (see the module docs).
const LEE_MOSER_1000_BULK: f64 = 1.0 / 0.0500256;
const LEE_MOSER_1000_CF: f64 = 2.0 * 0.0500256 * 0.0500256;

/// Lee & Moser (2015), Re_τ = 182.088: u_τ = 0.0637309 at U_b = 1 (ν = 3.5 × 10⁻⁴, δ = 1), from
/// the header of LM_Channel_0180_mean_prof.dat (same source as the Re_τ 1000 values).
const LEE_MOSER_180_RE_TAU: f64 = 182.088;
const LEE_MOSER_180_BULK: f64 = 1.0 / 0.0637309;
const LEE_MOSER_180_CF: f64 = 2.0 * 0.0637309 * 0.0637309;
/// The centreline U⁺ of Lee & Moser's Re_τ 182 profile, as the gate's window rounds it (17.75 to
/// 18.65 is 18.2 ± 2.5 %; the tabulated mean at y = δ reads about 18.27).
const LEE_MOSER_180_CENTRE: f64 = 18.2;

/// What the lattice's momentum transfer alone does to C_f on `cells` per half-height, without a flow solve (the first-rows run prints it
/// before its runs): ν̃ at the solution's profile on the lattice's nodes, row 1 at the solution's velocity, each later interval stepping by the
/// total stress t(y_mid) = 1 − y_mid over the harmonic mean of the two nodes' ν + ν_t (the lattice's transfer rule), the bulk velocity the
/// trapezoid sum over the nodes with the wall interval from the solution's own profile (the harness uses the law's, within 0.1 u⁺ of it), C_f = 2 / U_b².
/// Returns the 1→2 increment's excess over the solution's, and C_f's from the solution's with the rule as it is and with the 1→2 increment
/// taken from the solution (what `local_balance_node_one` sets, the harmonic sum floor of the layer), both as fractions.
fn harmonic_transfer_floor(y: &[f64], nt: &[f64], u: &[f64], cf_1d: f64, cells: usize) -> (f64, f64, f64) {
    let at = |yd: f64, f: &[f64]| {
        let k = y.partition_point(|v| *v < yd).clamp(1, y.len() - 1);
        f[k - 1] + (f[k] - f[k - 1]) * (yd - y[k - 1]) / (y[k] - y[k - 1])
    };
    let h = 1.0 / cells as f64;
    let visc = |j: usize| sa_1d::NU + fieldcraft_lbm::sa::nu_t(at(j as f64 * h, nt), sa_1d::NU);
    let step = |j: usize| (1.0 - (j as f64 + 0.5) * h) * h / (2.0 / (1.0 / visc(j) + 1.0 / visc(j + 1)));
    let wall_interval: f64 = (0..200).map(|i| at((i as f64 + 0.5) * h / 200.0, u) * h / 200.0).sum();
    let cf = |first_step: f64| {
        let mut v = at(h, u);
        let mut sum = 0.0;
        for j in 1..cells {
            let next = v + if j == 1 { first_step } else { step(j) };
            sum += 0.5 * (v + next) * h;
            v = next;
        }
        let bulk = wall_interval + sum;
        2.0 / (bulk * bulk) / cf_1d - 1.0
    };
    let newton_step = at(2.0 * h, u) - at(h, u);
    (step(1) / newton_step - 1.0, cf(step(1)), cf(newton_step))
}

/// A channel run's parameters.
#[derive(Debug, Clone)]
struct Channel {
    re_tau: f64,
    /// Cells per half-height δ.
    cells: usize,
    /// The expected U_b⁺, which sets the time step with `lattice_bulk`.
    bulk_guess: f64,
    /// U_b in lattice units (about).
    lattice_bulk: f64,
    collision: Collision,
    /// The hybrid collision's fixed σ (`None`: its default, or another collision).
    sigma: Option<f64>,
    turbulence: Option<Turbulence>,
    /// The walls on lattice nodes (the on-node bounce: one boundary node layer more) rather than
    /// halfway between cells.
    on_node: bool,
    max_steps: u64,
    sample: u64,
    settling: Settling,
}

/// A channel run's time means over its last averaging window.
#[derive(Debug, Clone)]
struct ChannelResult {
    /// U_b / u_τ, with u_τ from the carried stress. On on-node walls the interval between each
    /// wall and its first node is integrated from the wall law
    /// (common::cases::WallInterval); `bulk_trapezoid` is the same with the trapezoid from u = 0 at the wall.
    bulk: f64,
    bulk_trapezoid: f64,
    /// The centreline velocity over u_τ (the mean of the two middle rows when there is no row on
    /// the centreline).
    centre: f64,
    /// C_f from the carried, applied and requested wall stress, on the measured U_b (the law's);
    /// `carried_trapezoid` the carried C_f on the trapezoid's U_b.
    carried: f64,
    carried_trapezoid: f64,
    applied: f64,
    requested: Option<f64>,
    /// The applied wall drag over the body force (1 when fully developed).
    balance: f64,
    /// Every step's momentum budget closed (`wall_balance::closed_step`).
    budget_ok: bool,
    /// The smallest ν̃ of the run's last step (SA only).
    min_nu_tilde: Option<f64>,
    /// The carried C_f's block statistics over the whole run after the minimum flow-throughs,
    /// one block per averaging window.
    blocks: Blocks,
    settled: bool,
    steps: u64,
    /// The mean u and ν̃ (lattice units) at the first rows from the wall node (row 0 is the node),
    /// the two walls averaged, over the last averaging window.
    profile: Vec<[f64; 5]>,
    /// The friction velocity (lattice units) and the viscosity.
    u_tau: f64,
    nu: f64,
}

impl Channel {
    fn spacing(&self) -> f64 {
        1.0 / self.cells as f64
    }

    fn model(&self) -> fieldcraft_model::Model {
        let dx = self.spacing();
        let height = if self.on_node { 2.0 + dx } else { 2.0 };
        let mut d = domain([0.0; 3], [dx, height, dx], dx);
        for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
            *f = FaceCondition::Periodic {};
        }
        d.acceleration = [1.0, 0.0, 0.0];
        d.initial_velocity = [self.bulk_guess, 0.0, 0.0];
        let dt = self.lattice_bulk * dx / self.bulk_guess;
        let mut p = collide(self.max_steps as f64 * dt, dt, None, self.collision);
        if let Procedure::Fluid { turbulence, sigma, .. } = &mut p {
            *turbulence = self.turbulence;
            *sigma = self.sigma.map(fieldcraft_model::HrrSigma::Value);
        }
        model(1.0, 1.0 / self.re_tau, d, p)
    }

    /// Run to settling or the step cap.
    fn run(&self) -> ChannelResult {
        self.run_with(false)
    }

    /// [`Channel::run`] with `freeze` ν̃ held at the one-dimensional solution's profile.
    fn run_with(&self, freeze: bool) -> ChannelResult {
        let (setup, mut flow) = prepare(&self.model(), &[]);
        let grid = *flow.grid();
        if freeze {
            // ν̃/ν as a function of y/δ, y the cells from the nearer wall plane.
            let nu = (setup.tau_plus - 0.5) / 3.0;
            let g = sa_1d::grid(1200, 0.05 * sa_1d::NU, 1.0);
            let (nt, y, _, _) = sa_1d::solve_newton(sa_1d::NU, &g);
            let ratio = |yd: f64| {
                let k = y.partition_point(|v| *v < yd).clamp(1, y.len() - 1);
                (nt[k - 1] + (nt[k] - nt[k - 1]) * (yd - y[k - 1]) / (y[k] - y[k - 1])) / sa_1d::NU
            };
            let rows = grid.cells[1];
            let mut v = vec![0.0; grid.len()];
            for (c, out) in v.iter_mut().enumerate() {
                let r = grid.coords(c)[1];
                *out = nu * ratio(r.min(rows - 1 - r) as f64 / self.cells as f64);
            }
            flow.freeze_sa_nu_tilde(v).unwrap();
        }
        let rows = grid.cells[1];
        let a = flow.acceleration();
        let nu = (setup.tau_plus - 0.5) / 3.0;
        // Wall faces per wall: one cell in x and z.
        let area = 2.0 * (grid.cells[0] * grid.cells[2]) as f64;
        let musker = matches!(self.turbulence, Some(Turbulence::Smagorinsky { wall_model: Some(WallModel::Musker), .. }))
            || matches!(self.turbulence, Some(Turbulence::SpalartAllmaras { wall_model: Some(WallModel::Musker), .. }));
        // The bulk velocity: the mean over the rows between the walls (the on-node walls' rows,
        // on the walls themselves, at half weight), `[with the wall intervals from the wall law,
        // with the trapezoid]`: the two differ on on-node walls only.
        let sa_law = matches!(self.turbulence, Some(Turbulence::SpalartAllmaras { wall_model: Some(WallModel::SpalartAllmarasBounce), .. }));
        let interval = |u1: f64| if sa_law { WallInterval::new_sa_law(u1, 1.0, nu) } else { WallInterval::new(u1, 1.0, nu) };
        let bulk = |flow: &Flow| -> [f64; 2] {
            let row = |j: usize| {
                (0..grid.cells[0] * grid.cells[2])
                    .map(|k| flow.moments(grid.index([k % grid.cells[0], j, k / grid.cells[0]])).map_or(0.0, |(_, u)| u[0]))
                    .sum::<f64>()
                    / (grid.cells[0] * grid.cells[2]) as f64
            };
            let sum: f64 = (0..rows).map(|j| if self.on_node && (j == 0 || j + 1 == rows) { 0.5 * row(j) } else { row(j) }).sum();
            let height = if self.on_node { (rows - 1) as f64 } else { rows as f64 };
            // Each wall's interval, one cell (lattice units) from its node to its first fluid node.
            let law: f64 = if self.on_node {
                [row(1), row(rows - 2)].iter().filter_map(|&u1| interval(u1)).map(|w| w.law[0] - w.trapezoid[0]).sum()
            } else {
                0.0
            };
            [(sum + law) / height, sum / height]
        };
        let centre = |flow: &Flow| -> f64 {
            let row = |j: usize| {
                (0..grid.cells[0] * grid.cells[2])
                    .map(|k| flow.moments(grid.index([k % grid.cells[0], j, k / grid.cells[0]])).map_or(0.0, |(_, u)| u[0]))
                    .sum::<f64>()
                    / (grid.cells[0] * grid.cells[2]) as f64
            };
            if rows % 2 == 1 { row(rows / 2) } else { 0.5 * (row(rows / 2 - 1) + row(rows / 2)) }
        };
        let acoustic = 2.0 * rows as f64 * 3f64.sqrt();
        let bulk_lattice = self.lattice_bulk;
        let flow_through = 2.0 * std::f64::consts::PI * self.cells as f64 / bulk_lattice;
        let window = ((flow_through.max(2.0 * acoustic) / self.sample as f64).round() as usize).max(1);
        let min_steps = (self.settling.min_flow_throughs * flow_through).ceil() as u64;
        let mut blocks = Blocks::new(window, min_steps / self.sample).unwrap();
        // Per sample interval: [carried, applied, requested, body, U_b (law), U_b (trapezoid)]
        // summed, then averaged.
        let mut history: Vec<[f64; 7]> = Vec::new();
        let mut rows_history: Vec<Vec<[f64; 5]>> = Vec::new();
        let mut sums = [0.0; 7];
        let mut any_requested = false;
        let mut budget_ok = true;
        let mut p = wall_balance::population_momentum(&flow);
        let mut steps = 0;
        let mean =
            |rows: &[[f64; 7]]| -> [f64; 7] { [0, 1, 2, 3, 4, 5, 6].map(|c| rows.iter().map(|r| r[c]).sum::<f64>() / rows.len().max(1) as f64) };
        let cf = |stress: f64, ub: f64, rho: f64| 2.0 * stress / (rho * ub * ub);
        let rho = flow.mass() / (0..grid.len()).filter(|&c| flow.is_fluid(c)).count() as f64;
        let mut settled = false;
        while steps < self.max_steps {
            flow.step().unwrap();
            steps += 1;
            let next = wall_balance::population_momentum(&flow);
            let m = flow.mass();
            let wall: [f64; 3] = [0, 1, 2].map(|k| flow.forces().iter().map(|f| f[k]).sum());
            let density = |cell: usize| flow.moments(cell).map_or(f64::NAN, |(r, _)| r);
            // With SA under a wall model, the law's ρ u_τ² at each first fluid node (along the flow,
            // +x); otherwise Musker's request where the wall model is Musker.
            let requested = match (flow.sa_wall_function(), flow.walls()) {
                ((Some(_), samples), _) => [samples.iter().map(|&(cell, _, u_tau)| density(cell) * u_tau * u_tau).sum(), 0.0, 0.0],
                (_, Some(walls)) if musker => walls
                    .wall_cells()
                    .filter(|w| w.4.is_some())
                    .filter_map(|(cell, q, distance, _, _)| {
                        let target = walls.modelled(cell)?.0;
                        let normal = fieldcraft_lbm::lattice::C[q].map(|v| f64::from(v) * 0.5 / distance);
                        wall_balance::musker_request(target, normal, distance, nu).ok().map(|r| r.map(|r| r * density(cell)))
                    })
                    .fold([0.0; 3], |s, r| [0, 1, 2].map(|k| s[k] + r[k])),
                _ => [0.0; 3],
            };
            any_requested |= requested[0] != 0.0;
            let report =
                wall_balance::closed_step(p, next, a.map(|a| a * m), wall, requested.map(|r| r / area), area, m, 1e-12, f64::EPSILON).unwrap();
            budget_ok &= report.budget_ok;
            p = next;
            for (s, v) in sums.iter_mut().zip([report.carried[0], report.applied[0], requested[0] / area, a[0] * m, 0.0, 0.0, 0.0]) {
                *s += v;
            }
            if steps % self.sample != 0 {
                continue;
            }
            let n = self.sample as f64;
            let [law, trapezoid] = bulk(&flow);
            let row = [sums[0] / n, sums[1] / n, sums[2] / n, sums[3] / n, law, trapezoid, centre(&flow)];
            sums = [0.0; 7];
            assert!(row.iter().all(|v| v.is_finite()), "the channel is not finite after {steps} steps: {row:?}");
            // Past MAX_BLOCKS the stream refuses more samples; the statistics keep what they hold.
            blocks.push(cf(row[0], row[4], rho)).ok();
            history.push(row);
            let mean_row = |j: usize, f: &dyn Fn(usize) -> f64| -> f64 {
                (0..grid.cells[0] * grid.cells[2]).map(|k| f(grid.index([k % grid.cells[0], j, k / grid.cells[0]]))).sum::<f64>()
                    / (grid.cells[0] * grid.cells[2]) as f64
            };
            let nt = |c: usize| flow.spalart_allmaras().and_then(|f| f.nu_tilde().get(c).copied()).unwrap_or(0.0);
            let u = |c: usize| flow.moments(c).map_or(0.0, |(_, u)| u[0]);
            let upop = |c: usize| {
                let (m, r) = (0..fieldcraft_lbm::lattice::Q)
                    .fold((0.0, 0.0), |(m, r), q| (m + fieldcraft_lbm::lattice::C[q][0] as f64 * flow.population(q, c), r + flow.population(q, c)));
                if r > 0.0 { m / r } else { 0.0 }
            };
            let mom = |c: usize, cy: i32| {
                (0..fieldcraft_lbm::lattice::Q)
                    .filter(|&q| fieldcraft_lbm::lattice::C[q][1] == cy)
                    .map(|q| fieldcraft_lbm::lattice::C[q][0] as f64 * flow.population(q, c))
                    .sum::<f64>()
            };
            let pxy = |c: usize| {
                (0..fieldcraft_lbm::lattice::Q)
                    .map(|q| (fieldcraft_lbm::lattice::C[q][0] * fieldcraft_lbm::lattice::C[q][1]) as f64 * flow.population(q, c))
                    .sum::<f64>()
            };
            rows_history.push(
                (0..PROFILE_ROWS.min(rows / 2))
                    .map(|j| {
                        [
                            0.5 * (mean_row(j, &u) + mean_row(rows - 1 - j, &u)),
                            0.5 * (mean_row(j, &nt) + mean_row(rows - 1 - j, &nt)),
                            0.5 * (mean_row(j, &pxy) - mean_row(rows - 1 - j, &pxy)),
                            // The stress through the plane between rows j and j + 1 (j + ½ cells from the wall plane): the x-momentum the
                            // populations that arrived from the other side carry across it, bottom and top wall averaged.
                            0.5 * (mean_row(j, &|c| mom(c, -1)) - mean_row(j + 1, &|c| mom(c, 1)) + mean_row(rows - 1 - j, &|c| mom(c, 1))
                                - mean_row(rows - 2 - j, &|c| mom(c, -1))),
                            // The velocity the populations of the row carry (the wall node's row 0 has no fluid velocity otherwise).
                            0.5 * (mean_row(j, &upop) + mean_row(rows - 1 - j, &upop)),
                        ]
                    })
                    .collect(),
            );
            let carried: Vec<Vec<f64>> = history.iter().map(|r| vec![cf(r[0], r[4], rho)]).collect();
            let last = mean(&history[history.len().saturating_sub(window)..]);
            let balanced = (last[1] * area / last[3] - 1.0).abs() <= self.settling.drift;
            if balanced && self.settling.settled(&carried, window, steps as f64 / flow_through) {
                settled = true;
                break;
            }
        }
        let last = mean(&history[history.len().saturating_sub(window)..]);
        let u_tau = (last[0] / rho).sqrt();
        ChannelResult {
            bulk: last[4] / u_tau,
            bulk_trapezoid: last[5] / u_tau,
            centre: last[6] / u_tau,
            carried: cf(last[0], last[4], rho),
            carried_trapezoid: cf(last[0], last[5], rho),
            applied: cf(last[1], last[4], rho),
            requested: any_requested.then(|| cf(last[2], last[4], rho)),
            balance: last[1] * area / last[3],
            budget_ok,
            min_nu_tilde: flow
                .spalart_allmaras()
                .map(|f| (0..grid.len()).filter(|&c| flow.is_fluid(c)).filter_map(|c| f.nu_tilde().get(c).copied()).fold(f64::INFINITY, f64::min)),
            blocks,
            settled,
            steps,
            profile: {
                let tail = &rows_history[rows_history.len().saturating_sub(window)..];
                (0..PROFILE_ROWS.min(rows / 2))
                    .map(|j| [0, 1, 2, 3, 4].map(|c| tail.iter().filter_map(|r| r.get(j)).map(|r| r[c]).sum::<f64>() / tail.len().max(1) as f64))
                    .collect()
            },
            u_tau,
            nu,
        }
    }
}

/// The rows from each wall the first-rows run keeps (`Channel::run`'s profile).
const PROFILE_ROWS: usize = 8;

fn report(name: &str, r: &ChannelResult, bulk_reference: f64, cf_reference: f64) {
    let a = r.blocks.acceptance();
    println!(
        "{name}: {} steps, settled {}; U_b⁺ {:.4} (ref {bulk_reference:.4}, {:+.2} %; trapezoid {:.4}, {:+.2} %); C_f carried {:.5e} (ref {cf_reference:.5e}, {:+.2} %; trapezoid {:+.2} %), applied {:.5e}, requested {}; drag / body force {:.6}; budget closed {}; min ν̃ {:?}; {} blocks, half-width {:.2e}, accepted {}",
        r.steps,
        r.settled,
        r.bulk,
        100.0 * (r.bulk / bulk_reference - 1.0),
        r.bulk_trapezoid,
        100.0 * (r.bulk_trapezoid / bulk_reference - 1.0),
        r.carried,
        100.0 * (r.carried / cf_reference - 1.0),
        100.0 * (r.carried_trapezoid / cf_reference - 1.0),
        r.applied,
        r.requested.map_or("—".into(), |c| format!("{c:.5e}")),
        r.balance,
        r.budget_ok,
        r.min_nu_tilde,
        r.blocks.means().len(),
        r.blocks.half_width(),
        a.accepted,
    );
}

/// The measuring pipeline on the laminar channel: plane Poiseuille flow at Re_τ = 10, where
/// U_b⁺ = Re_τ/3 and C_f = 2/U_b⁺² = 18/Re_τ² exactly. On 16 cells per δ with halfway walls (TRT)
/// the discretisation error is small next to the 0.5 % asked of U_b⁺ and the carried C_f; the
/// budget closes every step and applied matches carried. The pipeline's settling is the SA
/// channel's own, [`GATED`].
#[test]
fn a_laminar_channel_reads_poiseuille_friction() {
    let re_tau = 10.0;
    let ch = Channel {
        re_tau,
        cells: 16,
        bulk_guess: re_tau / 3.0,
        lattice_bulk: 0.05,
        collision: Collision::Trt,
        turbulence: None,
        on_node: false,
        max_steps: 200_000,
        sample: 50,
        settling: GATED,
        sigma: None,
    };
    let r = ch.run();
    let (bulk, cf) = (re_tau / 3.0, 18.0 / (re_tau * re_tau));
    report("laminar channel, Re_τ 10, 16 cells per δ", &r, bulk, cf);
    assert!(r.settled, "the laminar channel did not settle in {} steps", r.steps);
    assert!(r.budget_ok);
    assert!((r.bulk / bulk - 1.0).abs() < 5e-3, "U_b⁺ {}", r.bulk);
    assert!((r.carried / cf - 1.0).abs() < 5e-3, "C_f {}", r.carried);
    assert!((r.applied / r.carried - 1.0).abs() < 1e-3, "applied {} vs carried {}", r.applied, r.carried);
}

/// The wall interval of an on-node wall (common::cases::WallInterval): its closed forms equal a
/// fine midpoint quadrature of Spalding's profile (inverted by bisection at each point) to 1e-9
/// relative over y₁⁺ 1–500. In the viscous sublayer the law is linear, so it equals the trapezoid
/// there. Unusable inputs give none.
#[test]
fn the_wall_interval_integrates_the_law_the_wall_applies() {
    use fieldcraft_lbm::wall::{BOUNCE_B, BOUNCE_KAPPA};
    let (k, e) = (BOUNCE_KAPPA, (-BOUNCE_KAPPA * BOUNCE_B).exp());
    let y_of = |u: f64| u + e * ((k * u).exp() - 1.0 - k * u - (k * u).powi(2) / 2.0 - (k * u).powi(3) / 6.0);
    let u_of = |y: f64| {
        let (mut lo, mut hi) = (0.0, y.max(1.0));
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if y_of(mid) < y { lo = mid } else { hi = mid }
        }
        0.5 * (lo + hi)
    };
    let nu = 1e-3;
    for y_plus in [1.0, 5.0, 11.0, 31.0, 50.0, 83.0, 200.0, 500.0] {
        let u_tau = 0.05;
        let y1 = y_plus * nu / u_tau;
        let u1 = u_tau * u_of(y_plus);
        let w = WallInterval::new(u1, y1, nu).unwrap();
        assert!((w.u_tau / u_tau - 1.0).abs() < 1e-12, "u_τ {} vs {u_tau}", w.u_tau);
        let n = 20_000;
        let h = y_plus / n as f64;
        let (mut i1, mut i2) = (0.0, 0.0);
        for j in 0..n {
            let u = u_of((j as f64 + 0.5) * h);
            i1 += u * h;
            i2 += u * u * h;
        }
        let want = [nu * i1, nu * u_tau * i2];
        for c in 0..2 {
            assert!((w.law[c] / want[c] - 1.0).abs() < 1e-7, "y⁺ {y_plus}: {:?} vs {want:?}", w.law);
        }
        println!("y⁺₁ {y_plus}: ∫u dy law / trapezoid {:.5}", w.law[0] / w.trapezoid[0]);
        if y_plus <= 1.0 {
            assert!((w.law[0] / w.trapezoid[0] - 1.0).abs() < 1e-3, "the sublayer is linear");
        }
    }
    for (u1, y1, nu) in [(0.0, 1.0, 1e-3), (f64::NAN, 1.0, 1e-3), (1.0, -1.0, 1e-3), (1.0, 1.0, 0.0), (1.0, f64::INFINITY, 1e-3)] {
        assert!(WallInterval::new(u1, y1, nu).is_none(), "{u1} {y1} {nu}");
    }
}

/// The reference numbers are Lee & Moser's header values (u_τ = 0.0500256 at U_b = 1).
#[test]
fn the_lee_moser_reference_is_their_published_bulk_and_friction_velocity() {
    assert!((LEE_MOSER_1000_BULK - 19.9898).abs() < 1e-4, "{LEE_MOSER_1000_BULK}");
    assert!((LEE_MOSER_1000_CF - 5.00512e-3).abs() < 1e-8, "{LEE_MOSER_1000_CF}");
    assert!((LEE_MOSER_180_BULK - 15.6910).abs() < 1e-4, "{LEE_MOSER_180_BULK}");
    assert!((LEE_MOSER_180_CF - 8.12326e-3).abs() < 1e-8, "{LEE_MOSER_180_CF}");
    assert!((LEE_MOSER_180_RE_TAU * 3.5e-4 - 0.0637309).abs() < 1e-6, "Re_τ ν = u_τ δ");
}

/// The model with wall model `wall`, inflow ν̃/ν 3 (the initial field).
fn sa(wall: WallModel) -> Turbulence {
    Turbulence::SpalartAllmaras { prandtl: 0.9, wall_model: Some(wall), inflow_ratio: 3.0 }
}

/// The frozen-ν̃ test mode holds the field after every step, the wall function's first-node values included, and refuses what it cannot hold.
#[test]
fn a_frozen_nu_tilde_is_held_and_bad_fields_are_refused() {
    let ch = Channel {
        re_tau: 1000.0,
        cells: 12,
        bulk_guess: LEE_MOSER_1000_BULK,
        lattice_bulk: 0.1,
        collision: Collision::Hrr,
        turbulence: Some(sa(WallModel::SpalartAllmarasBounce)),
        on_node: true,
        max_steps: 100,
        sample: 100,
        settling: GATED,
        sigma: Some(0.99),
    };
    let (setup, mut flow) = prepare(&ch.model(), &[]);
    let nu = (setup.tau_plus - 0.5) / 3.0;
    let n = flow.grid().len();
    assert!(flow.freeze_sa_nu_tilde(vec![0.0; n - 1]).is_err());
    assert!(flow.freeze_sa_nu_tilde(vec![f64::NAN; n]).is_err());
    let held: Vec<f64> = (0..n).map(|c| (3.0 + (c % 5) as f64) * nu).collect();
    flow.freeze_sa_nu_tilde(held.clone()).unwrap();
    for _ in 0..30 {
        flow.step().unwrap();
    }
    let got = flow.spalart_allmaras().map(|f| f.nu_tilde().to_vec()).unwrap_or_default();
    // The wall nodes are not fluid cells of the transport: ν̃ is 0 there.
    let nodes: Vec<usize> = flow.wall_nodes().map(|w| w.nodes().iter().map(|n| n.cell as usize).collect()).unwrap_or_default();
    assert!(!nodes.is_empty());
    for (c, (g, h)) in got.iter().zip(&held).enumerate() {
        assert!(if nodes.contains(&c) { *g == 0.0 } else { g == h }, "cell {c}: {g} vs {h}");
    }
}

/// With SA under a wall model the harness reads the wall function's request, ρ u_τ² at each first
/// fluid node: positive, finite and of the carried stress's order after a few hundred steps on
/// Han's on-node wall and on Musker's halfway one (a smoke test of the reading, not the gate).
/// Under a second in release.
#[test]
fn the_sa_wall_function_s_request_is_read() {
    for (wall, on_node, cells) in
        [(WallModel::SpaldingBounce, true, 20), (WallModel::SpalartAllmarasBounce, true, 20), (WallModel::Musker, false, 10)]
    {
        let ch = Channel {
            re_tau: 1000.0,
            cells,
            bulk_guess: LEE_MOSER_1000_BULK,
            lattice_bulk: 0.1,
            collision: Collision::Hrr,
            turbulence: Some(sa(wall)),
            on_node,
            max_steps: 400,
            sample: 100,
            settling: GATED,
            sigma: Some(0.99),
        };
        let r = ch.run();
        report(&format!("SA {wall:?}, {cells} cells per δ, 400 steps"), &r, LEE_MOSER_1000_BULK, LEE_MOSER_1000_CF);
        let requested = r.requested.unwrap_or(f64::NAN);
        assert!(requested > 0.0 && requested.is_finite(), "{wall:?}: requested {requested}");
        assert!(requested / r.carried > 0.1 && requested / r.carried < 10.0, "{wall:?}: requested {requested}, carried {}", r.carried);
        // The on-node wall's per-step budget leaves out its boundary nodes and is not exact in
        // the start-up; the halfway wall's closes every step.
        assert!(on_node || r.budget_ok, "{wall:?}");
    }
}

/// The Re_τ 1000 channel: SA under `wall` on Han's on-node wall, `cells` per δ, on `collision` (σ 0.99 with the hybrid one).
fn channel_at_re_tau_1000(wall: WallModel, cells: usize, collision: Collision, max_steps: u64) -> Channel {
    Channel {
        re_tau: 1000.0,
        cells,
        bulk_guess: LEE_MOSER_1000_BULK,
        lattice_bulk: 0.1,
        collision,
        sigma: (collision == Collision::Hrr).then_some(0.99),
        turbulence: Some(sa(wall)),
        on_node: true,
        max_steps,
        sample: 100,
        settling: GATED,
    }
}

/// The Spalart–Allmaras wall law (`spalart_allmaras_bounce`) on a periodic SA channel at
/// Re_τ 182, 12 and 24 cells per half-height (the first fluid node at y⁺ 15 and 7.6), on the regularized
/// collision and on the hybrid collision at σ 0.99 (MRT is below its stability bound at these viscosities). Spalding's law is biased low in the buffer layer (the finding
/// is in docs/verification.md); the gates are the centreline U⁺ within 2.5 % of Lee & Moser (the
/// window 17.75 to 18.65), carried C_f within 3.5 %, the request within 1 % of the carried
/// stress, settled, ν̃ ≥ 0. The C_f bar is the model's own error plus the grid's: the model's one-dimensional
/// solution is about 2 % low in C_f at this Re_τ and the runs sit within 1.5 % of it, so 3.5 % leaves no more than
/// the discretisation's share (the bar at Re_τ 1000, 5 %, is the 2.06 % the model is high there plus the same share).
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn sa_law_channel_at_re_tau_180_meets_lee_and_moser_at_12_and_24_cells() {
    let mut score = Score::default();
    let mut failures = Vec::new();
    for (collision, sigma) in [(Collision::Regularized, None), (Collision::Hrr, Some(0.99))] {
        for cells in [12, 24] {
            let r = Channel {
                re_tau: LEE_MOSER_180_RE_TAU,
                cells,
                bulk_guess: LEE_MOSER_180_BULK,
                lattice_bulk: 0.1,
                collision,
                sigma,
                turbulence: Some(sa(WallModel::SpalartAllmarasBounce)),
                on_node: true,
                max_steps: 2_000_000,
                sample: 100,
                settling: GATED,
            }
            .run();
            let name = format!("SA law, {collision:?}, {cells} cells per δ (first node y⁺ {:.1})", LEE_MOSER_180_RE_TAU / cells as f64);
            report(&name, &r, LEE_MOSER_180_BULK, LEE_MOSER_180_CF);
            println!("{name}: centre U⁺ {:.3} ({:+.2} % of {LEE_MOSER_180_CENTRE})", r.centre, 100.0 * (r.centre / LEE_MOSER_180_CENTRE - 1.0));
            score.push(Gate::new("SA law channel 180", &format!("{name}: centre U⁺"), r.centre, LEE_MOSER_180_CENTRE, 0.025));
            score.push(Gate::new("SA law channel 180", &format!("{name}: C_f carried"), r.carried, LEE_MOSER_180_CF, 0.035));
            score.push(Gate::new(
                "SA law channel 180",
                &format!("{name}: requested / carried"),
                r.requested.unwrap_or(f64::NAN) / r.carried,
                1.0,
                0.01,
            ));
            if !r.settled {
                failures.push(format!("{name}: not settled in {} steps", r.steps));
            }
            if !r.min_nu_tilde.is_some_and(|v| v >= 0.0) {
                failures.push(format!("{name}: ν̃ below 0 ({:?})", r.min_nu_tilde));
            }
        }
    }
    println!("{}", score.table());
    assert!(failures.is_empty(), "{failures:?}");
    assert!(score.passes());
}

/// The channel at Re_τ 1000 with the Spalart–Allmaras wall law, on the hybrid collision (σ 0.99), 20 and 32 cells per δ (the
/// first node at y⁺ 50 and 31): carried C_f within 5 % of Lee & Moser's, the request within 1 % of
/// the carried stress, settled, ν̃ ≥ 0.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn sa_law_channel_at_re_tau_1000_meets_lee_and_moser() {
    let mut score = Score::default();
    let mut failures = Vec::new();
    for cells in [20, 32] {
        let r = channel_at_re_tau_1000(WallModel::SpalartAllmarasBounce, cells, Collision::Hrr, 2_000_000).run();
        let name = format!("SA law, Hrr, {cells} cells per δ (first node y⁺ {:.0})", 1000.0 / cells as f64);
        report(&name, &r, LEE_MOSER_1000_BULK, LEE_MOSER_1000_CF);
        score.push(Gate::new("channel, SA law", &format!("{name}: C_f carried"), r.carried, LEE_MOSER_1000_CF, 0.05));
        score.push(Gate::new("channel, SA law", &format!("{name}: requested / carried"), r.requested.unwrap_or(f64::NAN) / r.carried, 1.0, 0.01));
        if !r.settled {
            failures.push(format!("{name}: not settled in {} steps", r.steps));
        }
        if !r.min_nu_tilde.is_some_and(|v| v >= 0.0) {
            failures.push(format!("{name}: ν̃ below 0 ({:?})", r.min_nu_tilde));
        }
    }
    println!("{}", score.table());
    assert!(failures.is_empty(), "{failures:?}");
    assert!(score.passes());
}

/// Lee & Moser's Re_τ 543.496 channel (header: ν = 1e-4, u_τ = 5.43496e-2, Re_b = 2 δ U_b / ν =
/// 20 000, so U_b = 1 in units of u_τ / 0.0543496 and C_f = 2 (u_τ / U_b)²); the centreline U⁺ of
/// their mean profile is 21.00.
const LEE_MOSER_550_RE_TAU: f64 = 543.496;
const LEE_MOSER_550_BULK: f64 = 1.0 / 0.0543496;
const LEE_MOSER_550_CF: f64 = 2.0 * 0.0543496 * 0.0543496;
const LEE_MOSER_550_CENTRE: f64 = 21.0;

/// The Spalart–Allmaras wall law at a third Reynolds number, Re_τ 543, with constants fixed
/// before this run: the regularized and the hybrid collision at σ 0.99, 16 and 24 cells per δ (the first node at
/// y⁺ 34 and 23). Gates: centreline U⁺ within 2.5 % of Lee & Moser, carried C_f within 5 %, the
/// request within 1 % of the carried stress.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn sa_law_channel_at_re_tau_550_meets_lee_and_moser() {
    let mut score = Score::default();
    let mut failures = Vec::new();
    for (collision, sigma) in [(Collision::Regularized, None), (Collision::Hrr, Some(0.99))] {
        for cells in [16, 24] {
            let r = Channel {
                re_tau: LEE_MOSER_550_RE_TAU,
                cells,
                bulk_guess: LEE_MOSER_550_BULK,
                lattice_bulk: 0.1,
                collision,
                sigma,
                turbulence: Some(sa(WallModel::SpalartAllmarasBounce)),
                on_node: true,
                max_steps: 2_000_000,
                sample: 100,
                settling: GATED,
            }
            .run();
            let name = format!("SA law, {collision:?}, {cells} cells per δ (first node y⁺ {:.0})", LEE_MOSER_550_RE_TAU / cells as f64);
            report(&name, &r, LEE_MOSER_550_BULK, LEE_MOSER_550_CF);
            println!("{name}: centre U⁺ {:.3} ({:+.2} % of {LEE_MOSER_550_CENTRE})", r.centre, 100.0 * (r.centre / LEE_MOSER_550_CENTRE - 1.0));
            score.push(Gate::new("SA law channel 550", &format!("{name}: centre U⁺"), r.centre, LEE_MOSER_550_CENTRE, 0.025));
            score.push(Gate::new("SA law channel 550", &format!("{name}: C_f carried"), r.carried, LEE_MOSER_550_CF, 0.05));
            score.push(Gate::new(
                "SA law channel 550",
                &format!("{name}: requested / carried"),
                r.requested.unwrap_or(f64::NAN) / r.carried,
                1.0,
                0.01,
            ));
            if !r.settled {
                failures.push(format!("{name}: not settled in {} steps", r.steps));
            }
            if !r.min_nu_tilde.is_some_and(|v| v >= 0.0) {
                failures.push(format!("{name}: ν̃ below 0 ({:?})", r.min_nu_tilde));
            }
        }
    }
    println!("{}", score.table());
    assert!(failures.is_empty(), "{failures:?}");
    assert!(score.passes());
}

/// The rows of the Re_τ 1000 channel beside the one-dimensional SA solution at the same y⁺: for each of the first rows from the wall
/// node, the lattice u⁺ and ν̃/(κ u_τ y) against the solution's, the carried C_f against the solution's, the stress rows, and the
/// velocity step of each of the first intervals against the harmonic prediction and the solution's. It prints, and asserts nothing: it
/// is how the C_f rows of the wall-layer corrections are produced, and runs only with `FIELDCRAFT_FIRST_ROWS=1` (hours at four grids):
///
/// ```text
/// FIELDCRAFT_FIRST_ROWS=1 [FIELDCRAFT_CHANNEL_CELLS=24,32,48] [FIELDCRAFT_CHANNEL_COLLISION=regularized|hrr]
///   [FIELDCRAFT_SA_FREEZE=newton]
///   cargo test --release -p fieldcraft-lbm --test sa_wall_law -- --ignored --nocapture first_rows
/// ```
///
/// `FIELDCRAFT_SA_FREEZE=newton` holds ν̃ at the solution's profile (the momentum equations then see the model's own eddy viscosity and
/// the layer's errors are told from the model's); `FIELDCRAFT_CHANNEL_COLLISION` defaults to `hrr` (σ 0.99).
#[test]
#[ignore = "local: FIELDCRAFT_FIRST_ROWS=1, hours"]
fn the_first_rows_of_the_sa_law_channel_beside_the_one_dimensional_solution() {
    const SWITCH: &str = "FIELDCRAFT_FIRST_ROWS";
    if std::env::var(SWITCH).is_err() {
        println!("the first-rows run prints only with {SWITCH}=1 (hours long): skipped");
        return;
    }
    let cells: Vec<usize> = std::env::var("FIELDCRAFT_CHANNEL_CELLS")
        .ok()
        .map(|v| v.split(',').filter_map(|c| c.trim().parse().ok()).collect())
        .unwrap_or_else(|| vec![16, 24, 32, 48]);
    let wall = WallModel::SpalartAllmarasBounce;
    let collision =
        if std::env::var("FIELDCRAFT_CHANNEL_COLLISION").is_ok_and(|v| v == "regularized") { Collision::Regularized } else { Collision::Hrr };
    let freeze = std::env::var("FIELDCRAFT_SA_FREEZE").is_ok_and(|v| v == "newton");
    let grid = sa_1d::grid(1200, 0.05 * sa_1d::NU, 1.0);
    let (nt, y, _, _) = sa_1d::solve_newton(sa_1d::NU, &grid);
    let u = sa_1d::velocity(&y, &nt, sa_1d::NU);
    let bulk: f64 = (1..y.len()).map(|i| 0.5 * (u[i] + u[i - 1]) * (y[i] - y[i - 1])).sum();
    let cf_1d = 2.0 / (bulk * bulk);
    let at = |yp: f64, f: &[f64]| {
        let target = yp * sa_1d::NU;
        let k = y.partition_point(|v| *v < target).clamp(1, y.len() - 1);
        f[k - 1] + (f[k] - f[k - 1]) * (target - y[k - 1]) / (y[k] - y[k - 1])
    };
    for &cells in &cells {
        let (excess, cf_rule, cf_floor) = harmonic_transfer_floor(&y, &nt, &u, cf_1d, cells);
        println!(
            "{cells} cells per δ, the momentum transfer alone: 1→2 increment {:+.1} % from the solution's; C_f {:+.2} % from the solution's with the harmonic rule, {:+.2} % with the 1→2 increment from the solution",
            100.0 * excess,
            100.0 * cf_rule,
            100.0 * cf_floor
        );
    }
    for cells in cells {
        let r = channel_at_re_tau_1000(wall, cells, collision, 2_000_000).run_with(freeze);
        println!(
            "{wall:?}, {cells} cells per δ: settled {}, C_f {:.5e}, {:+.2} % from the 1-D solution's {cf_1d:.5e}",
            r.settled,
            r.carried,
            100.0 * (r.carried / cf_1d - 1.0)
        );
        println!("  row   y+     u+ lattice   u+ 1-D     diff      nu~/(k y ut) lattice   1-D");
        for (j, [uj, ntj, _, _, _]) in r.profile.iter().enumerate().skip(1) {
            let yp = j as f64 * r.u_tau / r.nu;
            let (up, ratio) = (uj / r.u_tau, ntj / (fieldcraft_lbm::sa::KAPPA * j as f64 * r.u_tau));
            println!(
                "  {j:3} {yp:7.1} {up:11.3} {:10.3} {:+9.3} {ratio:16.3} {:10.3}",
                at(yp, &u),
                up - at(yp, &u),
                at(yp, &nt) / (fieldcraft_lbm::sa::KAPPA * yp * sa_1d::NU)
            );
        }
        // The stress rows: u_τ² is 1; expected 1 − y/δ with y the row's distance from the wall plane, δ the cells per half-height.
        println!("  row   expected   constitutive (nu+nu_t) dU/dy   population stress   population / expected");
        let cells_f = cells as f64;
        for j in 1..r.profile.len().saturating_sub(1) {
            let (u_tau2, nu) = (r.u_tau * r.u_tau, r.nu);
            let ntj = r.profile[j][1];
            let nu_t = fieldcraft_lbm::sa::nu_t(ntj, nu);
            let tau_bar = 3.0 * (nu + nu_t) + 0.5;
            let constitutive = (nu + nu_t) * 0.5 * (r.profile[j + 1][0] - r.profile[j - 1][0]) / u_tau2;
            let population = -(1.0 - 0.5 / tau_bar) * r.profile[j][2] / u_tau2;
            let expected = 1.0 - j as f64 / cells_f;
            println!("  {j:3} {expected:10.4} {constitutive:28.4} {population:19.4} {:22.4}", population / expected);
        }
        println!("  plane (cells from the wall plane)   stress / u_tau^2 through it (population flux)");
        for j in 0..r.profile.len().saturating_sub(1) {
            println!("  {:4.1} {:34.4}", j as f64 + 0.5, r.profile[j][3] / (r.u_tau * r.u_tau));
        }
        // The ledger in the reference solution's u_τ units: rows 0-4, the lattice u (row 0 from its populations) beside Newton's,
        // ν + ν_t per row, and per interval the lattice step, the harmonic prediction with τ_w (1 − y_mid/δ) and Newton's step.
        let u_ref = (cf_1d / 2.0).sqrt() * r.bulk * r.u_tau;
        let one = |yp: f64| at(yp, &u);
        println!("  ledger (u_ref {u_ref:.5e}, Newton's u_tau of the lattice bulk velocity)");
        println!("  row   u/u_ref   Newton u   (nu+nu_t)/nu");
        for j in 0..5.min(r.profile.len()) {
            let yp = j as f64 * u_ref / r.nu;
            let visc = 1.0 + fieldcraft_lbm::sa::nu_t(r.profile[j][1], r.nu) / r.nu;
            println!("  {j:3} {:9.4} {:10.4} {visc:14.3}", r.profile[j][4] / u_ref, one(yp));
        }
        println!("  interval   lattice step   harmonic prediction   Newton step");
        for j in 0..4.min(r.profile.len().saturating_sub(1)) {
            let (a, c) =
                (1.0 + fieldcraft_lbm::sa::nu_t(r.profile[j][1], r.nu) / r.nu, 1.0 + fieldcraft_lbm::sa::nu_t(r.profile[j + 1][1], r.nu) / r.nu);
            let harmonic = 2.0 / (1.0 / a + 1.0 / c);
            let dy_plus = u_ref / r.nu;
            let tau_mid = 1.0 - (j as f64 + 0.5) / cells as f64;
            let ys = [j as f64 * dy_plus, (j as f64 + 1.0) * dy_plus];
            println!(
                "  {j}-{}   {:12.4} {:21.4} {:13.4}",
                j + 1,
                (r.profile[j + 1][4] - r.profile[j][4]) / u_ref,
                tau_mid * dy_plus / harmonic,
                one(ys[1]) - one(ys[0])
            );
        }
        // The implied effective viscosity across the interval between rows j and j+1: the carried stress there
        // (the mean of the two rows' population stresses) over the measured du/dy, beside ν + ν_t at each row and their mean.
        println!("  interval   implied nu_eff/nu   (nu+nu_t)/nu at row j   at row j+1   mean");
        for j in 1..r.profile.len().saturating_sub(2) {
            let nu = r.nu;
            let stress = |k: usize| {
                let nu_t = fieldcraft_lbm::sa::nu_t(r.profile[k][1], nu);
                -(1.0 - 0.5 / (3.0 * (nu + nu_t) + 0.5)) * r.profile[k][2]
            };
            let implied = 0.5 * (stress(j) + stress(j + 1)) / (r.profile[j + 1][0] - r.profile[j][0]);
            let (a, c) = (1.0 + fieldcraft_lbm::sa::nu_t(r.profile[j][1], nu) / nu, 1.0 + fieldcraft_lbm::sa::nu_t(r.profile[j + 1][1], nu) / nu);
            println!("  {j}-{}   {:17.4} {a:23.4} {c:12.4} {:8.4}", j + 1, implied / nu, 0.5 * (a + c));
        }
    }
}
