//! TRT velocity-inlet reconstruction for flows without a free surface (#820).
//!
//! The Zou–He density/kinetic coupling admits an alternating boundary mode near tau=1/2.
//! We extrapolate the non-equilibrium stress trace from two interior cells, preserving
//! deviatoric stress. Extrapolating density outright instead destabilizes long acoustic waves.
//! Only the high-frequency part of that density correction is therefore applied:
//! `rho_b = rho_ZH + (I - S_t1 S_t2)(2 rho_1 - rho_2 - rho_ZH)`.
//! Each `S` is Shapiro's fourth-order filter `[-1,4,10,4,-1]/16`.
//!
//! The combination is FieldCraft's own, not a published inflow condition: it borrows the
//! equilibrium/non-equilibrium split of Guo, Zheng & Shi's extrapolation (2002,
//! doi:10.1088/1009-1963/11/4/310) and Shapiro's selective filter (1970,
//! doi:10.1029/RG008i002p00359). Published closures were tried first on the #820 box and kept
//! failing: full non-equilibrium extrapolation (long acoustic waves grow), the regularised inlet of
//! Latt et al. (2008, doi:10.1103/PhysRevE.77.056703), an equilibrium inlet, and a characteristic
//! (LODI) density with each of those reconstructions, after Izquierdo & Fueyo (2008,
//! doi:10.1103/PhysRevE.78.046707) and Heubes, Bartel & Ehrhardt (2014,
//! doi:10.1016/j.cam.2013.09.019). Below τ⁺ = 0.505 (reachable only with a turbulence model)
//! none of them, this one included, holds a 1e-6 disturbance at Mach 0.087, and with Λ = 3/16
//! the bulk stream is itself linearly unstable there for τ⁺ ≤ 0.502 (`tests/trt_stability.rs`);
//! the setup warns below 0.505.
//! The filter acts on the reconstruction residual, not on bulk flow or imposed velocity.
//! Its zero-frequency response vanishes. The stress correction has zero mass/momentum
//! and an isotropic second moment; it leaves the five deviatoric stress components alone.

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
