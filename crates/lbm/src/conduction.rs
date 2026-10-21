//! Conjugate heat transfer: heat conducting through the obstacles that conduct, coupled to the
//! temperature lattice of the fluid at every face they share.
//!
//! # Scheme
//! - **In the solid:** a cell-centred finite-volume conduction on the lattice's own cells
//!   (Patankar 1980, *Numerical Heat Transfer and Fluid Flow*, ch. 4): the heat through the face
//!   between two cells is `g (T_b − T_a)` with `g` the harmonic mean of their conductivities
//!   along the face's axis (Patankar §4.2-4: the series resistance of the two half cells, which
//!   keeps the flux continuous across a change of material), a held box face adds `2 k (T_w − T)`
//!   (half a cell to the wall), a flux face its flux. Explicit in time, sub-cycled so that every
//!   cell's sum of conductances over its heat capacity stays below one (the explicit limit). The
//!   heat a held or flux box face takes from the solid is counted as that face's heat.
//!   Conductivities may differ per axis (a board's in-plane and through-thickness values):
//!   the faces are axis-aligned, so each face uses its own axis's value.
//! - **At the fluid:** each D3Q7 link from a fluid cell into a conducting cell is an
//!   anti-bounce-back held at the interface temperature
//!   `T_i = (k_s T_s + k_f T_f) / (k_s + k_f)`, the temperature halfway between the two cell
//!   centres at which the conductive fluxes of the two half cells are equal (the same series
//!   resistance as between two solids). The heat the link takes from the fluid, `g* − g'`, is
//!   given to the solid cell, so the heat leaving the fluid is the heat the solid receives to
//!   round-off. That heat is affine in the solid's temperature, `a − b T_s` with
//!   `b = 2 w_q k_s / (k_s + k_f)`, so the solid takes it point-implicitly: each sub-step solves
//!   its cell's balance with the link term at the new temperature, and the link's population is
//!   then set at the sub-steps' mean temperature, which gives the solid exactly the heat it
//!   took. The coupling is then stable for any heat capacity (an explicit link needs
//!   `Σ b < C`, which a corner cell of aluminium in water already breaks). Boundary and
//!   interface conditions for the thermal lattice equation, the
//!   conjugate interface among them: L. Li, R. Mei, J. F. Klausner, *Boundary conditions for
//!   thermal lattice Boltzmann equation method*, J. Comput. Phys. 237:366, 2013.
//! - **Heat sources:** a solid's power is spread evenly over its cells.
//! - **Steady state:** a step run to steady state scales each solid's heat capacity down to
//!   [`smallest_capacity`] when that is smaller: the steady temperatures do not depend on it,
//!   while a solid's own time constant (seconds to minutes for a heat sink, against
//!   milliseconds for the air around it) would otherwise set the run's length. The transient of
//!   such a run is not physical, and the setup says so.
//!
//! # Units
//! Heat is counted in the fluid's lattice heat: one unit raises one fluid cell by one degree
//! (`(ρ c)_f Δx³`). A solid's conductivity is `k Δt / ((ρ c)_f Δx²)`, its capacity
//! `ρ c / (ρ c)_f`, so a solid with the fluid's properties behaves exactly as the fluid at rest
//! (the fluid's lattice conductivity is its lattice diffusivity `(τ_T − ½) / 4`).
//!
//! # Determinism
//! Every cell gathers from its own faces in a fixed order and the interface links are visited in
//! cell order, so the result does not depend on the number of threads.

use fieldcraft_model::{BoxFace, ThermalCondition};

use crate::lattice::C;
use crate::setup::SolidSetup;
use crate::voxel::Grid;
use crate::{Error, Result};

/// No conducting cell.
pub const NONE: u32 = u32::MAX;

/// Most solid sub-steps per lattice step.
pub const MAX_SUBSTEPS: u32 = 64;

/// The share of the explicit limit one sub-step uses.
const SUBSTEP_MARGIN: f64 = 0.9;

/// Sub-steps a solid whose capacity is scaled for the steady state aims at.
pub const STEADY_SUBSTEPS: f64 = 4.0;

/// The conductance of two half cells of conductivities `a` and `b` in series (unit spacing).
pub fn series(a: f64, b: f64) -> f64 {
    if a + b > 0.0 { 2.0 * a * b / (a + b) } else { 0.0 }
}

/// The smallest heat capacity (fluid lattice units) a cell of a solid with lattice
/// conductivities `k` may take so that its conduction needs at most [`STEADY_SUBSTEPS`]
/// sub-steps: a cell has at most six faces, each at most `2 k` (a held wall). The fluid links
/// are implicit and set no limit.
pub fn smallest_capacity(k: [f64; 3]) -> f64 {
    4.0 * (k[0] + k[1] + k[2]) / (SUBSTEP_MARGIN * STEADY_SUBSTEPS)
}

/// What one face of a conducting cell touches.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Face {
    /// Nothing that conducts: a non-conducting obstacle, an adiabatic or open box face, or the
    /// fluid (whose heat comes through the interface links).
    None,
    /// Another conducting cell (its index among the conducting cells) through conductance `g`.
    Solid { other: u32, g: f64 },
    /// Box face `wall` (its index) held at temperature `t`, half a cell away: conductance
    /// `g = 2 k`.
    Held { t: f64, g: f64, wall: u8 },
    /// Box face `wall` with a heat flux into the cell (per step).
    Flux { q: f64, wall: u8 },
}

impl Face {
    fn conductance(&self) -> f64 {
        match *self {
            Face::Solid { g, .. } | Face::Held { g, .. } => g,
            Face::None | Face::Flux { .. } => 0.0,
        }
    }
}

/// A conducting cell: its heat capacity, its source per step and its six faces (+x, −x, +y, −y,
/// +z, −z: toward the neighbour at `c_q` for the D3Q19 axis directions 1..7).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cell {
    pub capacity: f64,
    pub source: f64,
    pub faces: [Face; 6],
}

/// A D3Q7 link from fluid cell `cell` (direction `q`, whose population comes from the solid)
/// into conducting cell `solid` (its index among the conducting cells), in force group `group`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interface {
    pub cell: u32,
    pub q: u8,
    pub solid: u32,
    pub group: u32,
    /// The link's index among the flow's links (`Flow::links`).
    pub link: u32,
    /// The series weight `k_s / (k_s + k_f)` of the interface temperature: the link's
    /// conductance is `2 k_f weight` (`series(k_s, k_f)`).
    pub weight: f64,
}

/// The conducting cells and their state.
#[derive(Debug, Clone)]
pub struct SolidHeat {
    /// The conducting cells (grid indices), in cell order.
    pub cells: Vec<u32>,
    /// Per solid (the setup's order), its cells (grid indices), in cell order.
    pub by_solid: Vec<Vec<u32>>,
    pub cell: Vec<Cell>,
    /// The fluid links into conducting cells, in fluid cell order.
    pub interface: Vec<Interface>,
    /// Sub-steps of the conduction per lattice step.
    pub substeps: u32,
    /// Per conducting cell, the heat its fluid links take this step, `a − b T` summed over them.
    link_a: Vec<f64>,
    link_b: Vec<f64>,
    /// Per conducting cell, the mean of its temperatures over the last step's sub-steps (what
    /// its fluid links take their heat at).
    pub effective: Vec<f64>,
    /// The heat each box face took from the solids in the last step.
    pub wall_heat: [f64; 6],
    /// Each conducting cell's temperature, and the sub-steps' second buffer.
    pub temperature: Vec<f64>,
    next: Vec<f64>,
}

/// What [`SolidHeat::build`] needs of the grid and the domain.
pub struct Layout<'a> {
    pub grid: &'a Grid,
    pub periodic: [bool; 3],
    /// Per grid cell its index among the conducting cells, or [`NONE`].
    pub conducting: &'a [u32],
    /// Each box face's thermal condition when it is a wall (`None`: open or periodic).
    pub walls: [Option<ThermalCondition>; 6],
    /// The fluid's lattice conductivity.
    pub kf: f64,
}

impl SolidHeat {
    /// The conducting cells `cells` (grid indices, in order) of solids `kind` (indices into
    /// `solids`), their faces from `layout`, and the fluid links into them (`weight` filled in
    /// here), all at temperature `initial`. Fails when a heat capacity is not a positive number
    /// or the conduction would need more than [`MAX_SUBSTEPS`].
    pub fn build(
        layout: &Layout,
        solids: &[SolidSetup],
        cells: Vec<u32>,
        kind: &[usize],
        mut interface: Vec<Interface>,
        initial: f64,
    ) -> Result<SolidHeat> {
        let k_of = |kind: usize, axis: usize| solids.get(kind).and_then(|s| s.conductivity.get(axis).copied()).unwrap_or(0.0);
        let mut by_solid = vec![Vec::new(); solids.len()];
        for (&cell, &k) in cells.iter().zip(kind) {
            if let Some(list) = by_solid.get_mut(k) {
                list.push(cell);
            }
        }
        let g = layout.grid;
        let dims = g.cells.map(|d| d as i64);
        let mut cell_data = Vec::with_capacity(cells.len());
        for (own, (&cell, &k)) in cells.iter().zip(kind).enumerate() {
            let ijk = g.coords(cell as usize);
            let p = ijk.map(|v| v as i64);
            let mut faces = [Face::None; 6];
            for (slot, cq) in faces.iter_mut().zip(C.iter().skip(1)) {
                let axis = cq.iter().position(|c| *c != 0).unwrap_or(0);
                let ks = k_of(k, axis);
                let mut nb = [0, 1, 2].map(|a| p[a] + i64::from(cq[a]));
                if !layout.periodic[axis] && !(0..dims[axis]).contains(&nb[axis]) {
                    let face = BoxFace::of(axis, nb[axis] >= dims[axis]);
                    let wall = face.index() as u8;
                    *slot = match layout.walls.get(face.index()).copied().flatten() {
                        Some(ThermalCondition::Temperature { value }) => Face::Held { t: value, g: 2.0 * ks, wall },
                        Some(ThermalCondition::Flux { value }) => Face::Flux { q: value, wall },
                        Some(ThermalCondition::Adiabatic {}) | None => Face::None,
                    };
                    continue;
                }
                nb[axis] = nb[axis].rem_euclid(dims[axis]);
                let at = g.index(nb.map(|v| v as usize));
                // (a periodic image of the cell itself exchanges nothing)
                if let Some(&other) = layout.conducting.get(at).filter(|o| **o != NONE && **o as usize != own) {
                    let kb = kind.get(other as usize).map_or(0.0, |kb| k_of(*kb, axis));
                    *slot = Face::Solid { other, g: series(ks, kb) };
                }
            }
            let s = solids.get(k);
            let count = by_solid.get(k).map_or(1, |l| l.len().max(1)) as f64;
            cell_data.push(Cell { capacity: s.map_or(1.0, |s| s.capacity), source: s.map_or(0.0, |s| s.power / count), faces });
        }
        for i in &mut interface {
            let axis = (i.q as usize).checked_sub(1).map_or(0, |d| d / 2);
            let ks = kind.get(i.solid as usize).map_or(0.0, |k| k_of(*k, axis));
            i.weight = ks / (ks + layout.kf);
        }

        let m = cells.len();
        let mut worst: f64 = 0.0;
        for c in &cell_data {
            let cap = c.capacity;
            if !(cap > 0.0 && cap.is_finite()) {
                return Err(Error::Invalid(format!("a conducting cell has heat capacity {cap}")));
            }
            worst = worst.max(c.faces.iter().map(Face::conductance).sum::<f64>() / cap);
        }
        // Explicit conduction is stable while Σ g Δt' / C ≤ 1; keep a margin.
        let substeps = (worst / SUBSTEP_MARGIN).ceil().max(1.0);
        if !(substeps <= MAX_SUBSTEPS as f64) {
            return Err(Error::Invalid(format!(
                "a conducting solid diffuses {worst:.1} cells per step on this lattice, needing {substeps} sub-steps (at most {MAX_SUBSTEPS}): shorten the time step"
            )));
        }
        Ok(SolidHeat {
            link_a: vec![0.0; m],
            link_b: vec![0.0; m],
            effective: vec![initial; m],
            wall_heat: [0.0; 6],
            temperature: vec![initial; m],
            next: vec![initial; m],
            cells,
            by_solid,
            cell: cell_data,
            interface,
            substeps: substeps as u32,
        })
    }

    /// Add a fluid link into conducting cell `solid` whose heat is `a − b T` at the cell's
    /// temperature `T`.
    pub fn link(&mut self, solid: u32, a: f64, b: f64) {
        if let (Some(x), Some(y)) = (self.link_a.get_mut(solid as usize), self.link_b.get_mut(solid as usize)) {
            *x += a;
            *y += b;
        }
    }

    /// Advance the conduction by one lattice step: the links implicit, the faces explicit, the
    /// sources and the link heat spread evenly over the sub-steps. Fills [`SolidHeat::effective`]
    /// and [`SolidHeat::wall_heat`]. Returns false when a temperature is not finite.
    pub fn step(&mut self) -> bool {
        let share = 1.0 / f64::from(self.substeps);
        self.effective.fill(0.0);
        self.wall_heat = [0.0; 6];
        let mut finite = true;
        for _ in 0..self.substeps {
            let old = &self.temperature;
            finite = true;
            let cells = self.next.iter_mut().zip(old).zip(&self.cell).zip(self.link_a.iter().zip(&self.link_b)).zip(&mut self.effective);
            for ((((slot, t), c), (a, b)), mean) in cells {
                let mut heat = a + c.source;
                for face in &c.faces {
                    heat += match *face {
                        Face::None => 0.0,
                        Face::Solid { other, g } => g * (old.get(other as usize).copied().unwrap_or(*t) - t),
                        Face::Held { t: tw, g, wall } => {
                            let into = g * (tw - t);
                            if let Some(w) = self.wall_heat.get_mut(wall as usize) {
                                *w -= into * share;
                            }
                            into
                        }
                        Face::Flux { q, wall } => {
                            if let Some(w) = self.wall_heat.get_mut(wall as usize) {
                                *w -= q * share;
                            }
                            q
                        }
                    };
                }
                // C (T' − T) = share (heat − b T')
                *slot = (t + heat * share / c.capacity) / (1.0 + b * share / c.capacity);
                *mean += *slot * share;
                finite &= slot.is_finite();
            }
            std::mem::swap(&mut self.temperature, &mut self.next);
        }
        self.link_a.fill(0.0);
        self.link_b.fill(0.0);
        finite
    }

    /// The total heat in the conducting cells above zero temperature (fluid lattice heat units).
    pub fn heat_content(&self) -> f64 {
        self.temperature.iter().zip(&self.cell).map(|(t, c)| t * c.capacity).sum()
    }
}
