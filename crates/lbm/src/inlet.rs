//! TRT velocity-inlet reconstruction for flows without a free surface (#820).
//!
//! **Base closure.** Zou & He (1997, doi:10.1063/1.869307) with the D3Q19 transverse-momentum
//! correction of Hecht & Harting (2010, doi:10.1088/1742-5468/2010/01/P01018) sets the unknown
//! populations from the prescribed velocity. Under TRT near τ⁺ = ½ its density/kinetic coupling
//! admits an inlet-local mode, alternating from cell to cell along the face, that grows.
//!
//! **Correction**, each part from a paper (`docs/solvers/literature.md`, C3.3a–C3.3e):
//! - Extrapolation (Guo, Zheng & Shi 2002, doi:10.1088/1009-1963/11/4/310): their
//!   equilibrium/non-equilibrium split and two-point extrapolation `2 x₁ − x₂` from the inner
//!   cells, applied to the non-equilibrium stress trace and to the density only. The prescribed
//!   velocity and the five deviatoric stress components stay Zou–He's.
//! - Residual filter (Shapiro 1970, doi:10.1029/RG008i002p00359; 1975,
//!   doi:10.1090/S0025-5718-1975-0389356-X): the n = 2 filter `S = [-1,4,10,4,-1]/16`, response
//!   `1 − sin⁴(k/2)`, along each tangential axis. Only the high-wavenumber part of the density
//!   residual is applied: `rho_b = rho_ZH + (I − S_t1 S_t2)(2 rho_1 − rho_2 − rho_ZH)`. Filtering
//!   only the unstable high wavenumbers follows Ricot, Marié, Sagaut & Bailly (2009,
//!   doi:10.1016/j.jcp.2009.03.030), who filter the bulk; here only the inlet residual.
//!
//! **What is new.** The cited inflow conditions do not combine these; the combination is FieldCraft's
//! own. The variants measured while building it fail: the density filter alone keeps the mode; the
//! trace and density both extrapolated outright fail on a broadband seed (step 2118 at τ⁺ = 0.505)
//! because long acoustic waves grow; the filter with the normal stress extrapolated instead of the
//! trace leaves oblique modes growing (about 2e-3 per step). The high-pass keeps
//! the long waves at Zou–He's density (its response falls as k⁴) and removes the two-cell wave.
//! A sponge or PML layer (Xu & Sagaut 2013; Najafi-Yazdi & Mongeau 2012) absorbs waves crossing
//! it but leaves the closure that feeds this mode; neither was implemented.
//!
//! **Consistency.** Mass: each tangential filter pass is symmetric with unit row sums (a wall side
//! reflects the stencil, a missing neighbour puts its weight on the centre), so the density
//! correction sums to zero over the face. This preserves Zou–He's net mass flux for uniform
//! prescribed normal velocity; for a varying profile the flux correction is `Σ δρᵢ u_n,ᵢ`, which need not vanish
//! (`the_density_correction_adds_no_mass_over_a_face`;
//! `a_duct_with_the_corrected_inlet_balances_its_mass_at_steady_state`). Momentum: the correction
//! is `equilibrium(δρ, u)`, carrying momentum `δρ u`, plus a trace term with zero mass and momentum,
//! so the prescribed velocity is exact (`correction_preserves_velocity_and_deviatoric_stress`).
//! Accuracy: second order on Poiseuille.
//!
//! **Against the published closures**, implemented from their papers and run on the same cases
//! (TRT, Λ = 3/16, U = 0.05): on the #820 16³ slip box at
//! τ⁺ = 0.505, 0.507, 0.51 and 0.52 this inlet holds 30 000 steps with checkerboard and broadband
//! seeds. Zou–He fails the broadband seed at every τ⁺ (from step 400), Guo et al.'s extrapolation fails at every
//! τ⁺ (from step 300), Latt et al.'s regularised inlet (2008, doi:10.1103/PhysRevE.77.056703) below
//! 0.52 (from step 400), and so do LODI characteristic inlets (Izquierdo
//! & Fueyo 2008, doi:10.1103/PhysRevE.78.046707; Heubes, Bartel & Ehrhardt 2014,
//! doi:10.1016/j.cam.2013.09.019) with bounce-back or regularised populations at every τ⁺ tried.
//!
//! CPU comparison (also in `docs/verification.md`), without a turbulence model. Box entries are
//! steps survived of 30 000 at τ⁺ 0.505 / 0.507 / 0.51 / 0.52, stopped on a solver error or
//! density departure above 0.1. Box survival is checked in 50-step chunks; a solver error reports
//! the last fully completed chunk. The walled 48 × 16 × 16 duct has τ⁺ = 0.52 and runs 60 000
//! steps in 1000-step chunks; drift and flux errors are measured over its last 10 000. Poiseuille has τ⁺ = 0.6,
//! H = 8/16/32; errors are relative to peak velocity, and order uses H = 16/32.
//!
//! | Closure | Box checkerboard | Box broadband | Duct steps; mass drift/step; flux vs ρ_b U A; in vs out | Poiseuille mid error at H=32, order; next to inlet |
//! |---|---|---|---|---|
//! | FieldCraft | all 30 000 | all 30 000 | 60 000; 1e-20; 8e-16; 7.9e-5 | 7.7e-5, 1.99; 5.6e-6 |
//! | Zou–He / Hecht–Harting | 2 800 / 12 200 / 30 000 / 30 000 | 400 / 450 / 600 / 5 700 | fails at 11 000 | 7.7e-5, 2.00; 6.8e-6 |
//! | Guo extrapolation | 300 / 350 / 500 / 23 950 | 300 / 350 / 500 / 27 050 | fails before 1 000 | 3.1e-4, 1.98; 2.7e-4 |
//! | Latt regularised | 400 / 550 / 950 / 30 000 | 450 / 550 / 1 000 / 30 000 | 60 000; 4e-20; 4e-16; 8.0e-5 | 7.7e-5, 2.00; 5.1e-5 |
//! | LODI bounce-back | 1 200 / 1 750 / 6 650 / 20 400 | 1 150 / 1 700 / 6 850 / 22 000 | 60 000; 0; −3.7e-4; 7.7e-5 | 7.6e-5, 1.97; 9.2e-5 |
//! | LODI regularised | 300 / 400 / 550 / 11 250 | 300 / 350 / 500 / 11 500 | fails before 1 000 | 9.3e-5, 2.08; 1.5e-4 |
//!
//! These are finite-grid checks, not a general stability proof. The comparison harness swapped
//! the closure in place and is not committed; the committed tests cover the adopted inlet.
//! No sponge or PML was built or measured.
//!
//! In the tested cases below τ⁺ = 0.505, reachable only with a turbulence model, none holds a 1e-6 disturbance at
//! Mach 0.087, this one included; with Λ = 3/16 the bulk is itself linearly unstable there for
//! τ⁺ ≤ 0.502 (`tests/trt_stability.rs`), and the setup warns below 0.505.

use crate::flow::{OpenCell, Target};
use crate::lattice::{C, Q, W, equilibrium, opposite};
use crate::outlet::Populations;
use crate::{Error, Grid, Result};

const MISSING: u32 = u32::MAX;
const WEIGHTS: [f64; 5] = [-1.0 / 16.0, 4.0 / 16.0, 10.0 / 16.0, 4.0 / 16.0, -1.0 / 16.0];

/// Indices of the two inner fluid cells and the two tangential five-point stencils.
/// Neighbours index the open-cell array. Missing/other-face neighbours use the centre.
#[derive(Debug, Clone, Copy)]
pub struct Stencil {
    pub inner: [u32; 2],
    pub neighbours: [[u32; 5]; 2],
}

impl Stencil {
    fn active(&self) -> bool {
        self.inner[0] != MISSING
    }
}

/// Geometry-derived stencils and scratch, overwritten before every boundary update.
#[derive(Debug, Clone, Default)]
pub(crate) struct Inlets {
    stencils: Vec<Stencil>,
    correction: Vec<[f64; 2]>,
    filtered: Vec<f64>,
}

fn filled<T: Clone>(len: usize, value: T) -> Result<Vec<T>> {
    let mut out = Vec::new();
    out.try_reserve_exact(len).map_err(|_| Error::Limit("not enough memory for inlet reconstruction".into()))?;
    out.resize(len, value);
    Ok(out)
}

impl Inlets {
    pub(crate) fn new(grid: &Grid, open: &[OpenCell], periodic: [bool; 3], fluid: &[bool]) -> Result<Self> {
        let mut stencils = filled(open.len(), Stencil { inner: [MISSING; 2], neighbours: [[0; 5]; 2] })?;
        let mut lookup = Vec::new();
        lookup.try_reserve_exact(open.len()).map_err(|_| Error::Limit("not enough memory for inlet lookup".into()))?;
        for (k, (oc, stencil)) in open.iter().zip(&mut stencils).enumerate() {
            let row = u32::try_from(k).map_err(|_| Error::Limit("too many inlet cells".into()))?;
            stencil.neighbours = [[row; 5]; 2];
            let Target::Velocity(velocity) = oc.target else { continue };
            if oc.axis >= 3 || velocity[oc.axis] * oc.inward < 0.0 {
                continue;
            }
            let mut xyz = grid.coords(oc.cell as usize);
            let mut cells = [MISSING; 2];
            for cell in &mut cells {
                let next = if oc.inward > 0.0 { xyz[oc.axis].checked_add(1) } else { xyz[oc.axis].checked_sub(1) };
                let Some(next) = next.filter(|v| *v < grid.cells[oc.axis]) else { break };
                xyz[oc.axis] = next;
                let inside = grid.index(xyz);
                if !fluid.get(inside).copied().unwrap_or(false) {
                    break;
                }
                *cell = u32::try_from(inside).map_err(|_| Error::Limit("inlet cell cannot be addressed".into()))?;
            }
            if cells[1] != MISSING {
                stencil.inner = cells;
                lookup.push((oc.cell, row));
            }
        }
        lookup.sort_unstable_by_key(|entry| entry.0);
        for (oc, stencil) in open.iter().zip(&mut stencils) {
            if !stencil.active() {
                continue;
            }
            let xyz = grid.coords(oc.cell as usize);
            for (pass, neighbours) in stencil.neighbours.iter_mut().enumerate() {
                let axis = (oc.axis + pass + 1) % 3;
                let len = grid.cells[axis] as i64;
                for (j, neighbour) in neighbours.iter_mut().enumerate() {
                    let coordinate = xyz[axis] as i64 + j as i64 - 2;
                    let coordinate = if periodic[axis] {
                        coordinate.rem_euclid(len)
                    } else {
                        let reflected = coordinate.rem_euclid(2 * len);
                        if reflected >= len { 2 * len - 1 - reflected } else { reflected }
                    };
                    let mut adjacent = xyz;
                    adjacent[axis] = coordinate as usize;
                    let cell = u32::try_from(grid.index(adjacent)).map_err(|_| Error::Limit("inlet cell cannot be addressed".into()))?;
                    if let Ok(found) = lookup.binary_search_by_key(&cell, |entry| entry.0)
                        && let Some(&(_, row)) = lookup.get(found)
                        && open.get(row as usize).is_some_and(|other| other.axis == oc.axis && other.inward == oc.inward)
                    {
                        *neighbour = row;
                    }
                }
            }
        }
        Ok(Self { stencils, correction: filled(open.len(), [0.0; 2])?, filtered: filled(open.len(), 0.0)? })
    }

    pub(crate) fn stencils(&self) -> &[Stencil] {
        &self.stencils
    }

    pub(crate) fn heap_bytes(&self) -> u64 {
        (self.stencils.capacity() * size_of::<Stencil>()
            + self.correction.capacity() * size_of::<[f64; 2]>()
            + self.filtered.capacity() * size_of::<f64>()) as u64
    }

    pub(crate) fn stage(&mut self, open: &[OpenCell], pop: &impl Populations) {
        for ((oc, stencil), correction) in open.iter().zip(&self.stencils).zip(&mut self.correction) {
            *correction = [0.0; 2];
            let Target::Velocity(u) = oc.target else { continue };
            if !stencil.active() {
                continue;
            }
            let cell = oc.cell as usize;
            let half_force = pop.half_force(cell);
            let bare: [f64; 3] = std::array::from_fn(|a| u[a] - half_force[a]);
            let mut known = 0.0;
            for (q, c) in C.iter().enumerate() {
                let cn = f64::from(c[oc.axis]) * oc.inward;
                if cn <= 0.0 {
                    known += pop.get(q, cell) * if cn < 0.0 { 2.0 } else { 1.0 };
                }
            }
            let rho = known / (1.0 - oc.inward * bare[oc.axis]);
            let mut trace = -rho * (1.0 + bare.iter().map(|v| v * v).sum::<f64>());
            for (q, c) in C.iter().enumerate() {
                let v = if f64::from(c[oc.axis]) * oc.inward > 0.0 {
                    let cu: f64 = c.iter().zip(bare).map(|(c, u)| f64::from(*c) * u).sum();
                    pop.get(opposite(q), cell) + 6.0 * W[q] * rho * cu
                } else {
                    pop.get(q, cell)
                };
                trace += v * c.iter().map(|v| f64::from(v * v)).sum::<f64>();
            }
            let Some((r1, t1)) = bulk_trace(pop, stencil.inner[0] as usize) else { continue };
            let Some((r2, t2)) = bulk_trace(pop, stencil.inner[1] as usize) else { continue };
            *correction = [2.0 * r1 - r2 - rho, 2.0 * t1 - t2 - trace];
        }
        for (stencil, filtered) in self.stencils.iter().zip(&mut self.filtered) {
            *filtered = stencil.neighbours[0]
                .iter()
                .zip(WEIGHTS)
                .map(|(row, weight)| weight * self.correction.get(*row as usize).map_or(0.0, |c| c[0]))
                .sum();
        }
        for (stencil, correction) in self.stencils.iter().zip(&mut self.correction) {
            if stencil.active() {
                correction[0] -= stencil.neighbours[1]
                    .iter()
                    .zip(WEIGHTS)
                    .map(|(row, weight)| weight * self.filtered.get(*row as usize).copied().unwrap_or(0.0))
                    .sum::<f64>();
            }
        }
    }

    pub(crate) fn apply(&self, k: usize, populations: &mut [f64; Q], bare: [f64; 3]) {
        if let Some(&[density, trace]) = self.correction.get(k) {
            for (q, value) in populations.iter_mut().enumerate() {
                *value += equilibrium(q, density, bare) + 1.5 * W[q] * (C[q].iter().map(|v| f64::from(v * v)).sum::<f64>() - 1.0) * trace;
            }
        }
    }
}

fn bulk_trace(pop: &impl Populations, cell: usize) -> Option<(f64, f64)> {
    let (mut rho, mut momentum, mut second) = (0.0, [0.0; 3], 0.0);
    for (q, c) in C.iter().enumerate() {
        let f = pop.get(q, cell);
        rho += f;
        for (m, c) in momentum.iter_mut().zip(c) {
            *m += f * f64::from(*c);
        }
        second += f * c.iter().map(|v| f64::from(v * v)).sum::<f64>();
    }
    (rho > 0.0 && rho.is_finite()).then(|| (rho, second - rho - momentum.iter().map(|v| v * v).sum::<f64>() / rho))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_inner_cells_and_velocity_outflow_keep_the_original_closure() {
        for nx in [2, 3] {
            let grid = Grid { origin: [0.0; 3], spacing: 1.0, cells: [nx, 1, 1] };
            let mut fluid = vec![true; grid.len()];
            for speed in [0.05, -0.05] {
                let open = [OpenCell { cell: 0, axis: 0, inward: 1.0, target: Target::Velocity([speed, 0.0, 0.0]) }];
                let inlet = Inlets::new(&grid, &open, [false; 3], &fluid).unwrap();
                assert_eq!(inlet.stencils[0].active(), nx == 3 && speed > 0.0);
                if nx == 3 {
                    fluid[2] = false;
                    let blocked = Inlets::new(&grid, &open, [false; 3], &fluid).unwrap();
                    assert!(!blocked.stencils[0].active());
                    fluid[2] = true;
                }
            }
        }
    }

    /// Cells at equilibrium with a density and velocity that vary from cell to cell.
    struct Noise(Vec<(f64, [f64; 3])>);

    impl Populations for Noise {
        fn get(&self, q: usize, cell: usize) -> f64 {
            self.0.get(cell).map_or(0.0, |(rho, u)| equilibrium(q, *rho, *u))
        }
        fn is_fluid(&self, _: usize) -> bool {
            true
        }
        fn half_force(&self, _: usize) -> [f64; 3] {
            [0.0; 3]
        }
    }

    /// Each tangential filter pass is symmetric with unit row sums, also where a wall side reflects
    /// the stencil and where a missing neighbour puts its weight on the centre, so the density
    /// correction sums to zero over the face. For uniform prescribed normal velocity, the mass flux Σ ρ_b u_n is
    /// therefore Zou–He's; this does not guarantee unchanged flux for a varying velocity profile.
    #[test]
    fn the_density_correction_adds_no_mass_over_a_face() {
        let grid = Grid { origin: [0.0; 3], spacing: 1.0, cells: [3, 7, 5] };
        let mut fluid = vec![true; grid.len()];
        // Two inlet cells without inner cells: holes in the filtered face.
        fluid[grid.index([1, 3, 2])] = false;
        fluid[grid.index([2, 5, 0])] = false;
        let pop = Noise(
            (0..grid.len())
                .map(|cell| {
                    let h = ((cell * 2_654_435_761) % 1000) as f64 / 1000.0 - 0.5;
                    (1.0 + 0.01 * h, [0.05 + 0.004 * h, 0.003 * h, -0.002 * h])
                })
                .collect(),
        );
        let open: Vec<OpenCell> = (0..7)
            .flat_map(|y| (0..5).map(move |z| [0, y, z]))
            .map(|xyz| OpenCell { cell: grid.index(xyz) as u32, axis: 0, inward: 1.0, target: Target::Velocity([0.05, 0.0, 0.0]) })
            .collect();
        for periodic in [[false; 3], [false, true, false], [false, true, true]] {
            let mut inlets = Inlets::new(&grid, &open, periodic, &fluid).unwrap();
            inlets.stage(&open, &pop);
            let largest = inlets.correction.iter().map(|c| c[0].abs()).fold(0.0, f64::max);
            let sum: f64 = inlets.correction.iter().map(|c| c[0]).sum();
            assert!(largest > 1e-5, "{periodic:?}: the residual must not vanish ({largest})");
            assert!(sum.abs() < 1e-15, "{periodic:?}: {sum}");
        }
    }

    /// The correction may change density/trace, but cannot change prescribed velocity or
    /// deviatoric non-equilibrium stress (including shear carried by the inlet).
    #[test]
    fn correction_preserves_velocity_and_deviatoric_stress() {
        let bare = [0.05, -0.01, 0.02];
        let mut f = std::array::from_fn(|q| equilibrium(q, 1.02, bare));
        // A nonzero shear and normal stress difference, with zero mass/momentum.
        for (q, value) in f.iter_mut().enumerate() {
            let c = C[q].map(f64::from);
            *value += W[q] * (0.0003 * c[0] * c[1] + 0.0002 * (c[1] * c[1] - c[2] * c[2]));
        }
        let moments = |f: &[f64; Q]| {
            let rho: f64 = f.iter().sum();
            let j: [f64; 3] = std::array::from_fn(|a| f.iter().zip(C).map(|(v, c)| v * f64::from(c[a])).sum());
            let stress: [[f64; 3]; 3] = std::array::from_fn(|a| {
                std::array::from_fn(|b| {
                    f.iter().zip(C).map(|(v, c)| v * f64::from(c[a] * c[b])).sum::<f64>() - j[a] * j[b] / rho - if a == b { rho / 3.0 } else { 0.0 }
                })
            });
            (rho, j.map(|v| v / rho), stress)
        };
        let (r0, u0, p0) = moments(&f);
        let inlets = Inlets { correction: vec![[0.001, -0.0007]], ..Inlets::default() };
        inlets.apply(0, &mut f, bare);
        let (r1, u1, p1) = moments(&f);
        assert!((r1 - r0 - 0.001).abs() < 1e-14);
        for a in 0..3 {
            assert!((u1[a] - u0[a]).abs() < 1e-14);
            for b in 0..3 {
                let expected = if a == b { -0.0007 / 3.0 } else { 0.0 };
                assert!((p1[a][b] - p0[a][b] - expected).abs() < 1e-14);
            }
        }
    }
}
