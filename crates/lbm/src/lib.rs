//! FieldCraft's lattice Boltzmann flow solver (CPU).
//!
//! It solves a fluid step of the model ([`fieldcraft_model::Procedure::Fluid`]) in its fluid
//! domain: weakly compressible, isothermal, Newtonian flow on a uniform grid of cubic cells that
//! fills the domain's box, less the cells inside the obstacle bodies.
//!
//! - [`voxel`]: which cells the obstacles take, decided with exact predicates from the
//!   face-tagged boundary surface, and which face each wall link meets.
//! - [`setup`]: the step in lattice units, the time step, and the stability limits.
//! - [`flow`]: the D3Q19 two-relaxation-time scheme with Guo forcing, halfway bounce-back
//!   walls (moving walls after Ladd), slip faces, Zou–He inlets and outlets, and the
//!   momentum-exchange force on every wall; fans (a pressure–flow curve) and porous zones
//!   (Darcy–Forchheimer resistance) as body forces.
//! - [`outlet`]: convective pressure outlets, which carry a wake or vortices out of the box.
//! - [`solve`]: the time loop, output frames and force histories in model units.
//! - [`steady`]: the steady-state monitor and stop rule of a run.
//!
//! # Sources
//! - T. Krüger, H. Kusumaatmaja, A. Kuzmin, O. Shardt, G. Silva, E. M. Viggen, *The Lattice
//!   Boltzmann Method: Principles and Practice*, Springer 2017 (velocity sets, equilibria, units,
//!   stream–collide, bounce-back, forcing, boundary conditions).
//! - Y. H. Qian, D. d'Humières, P. Lallemand, *Lattice BGK models for Navier–Stokes equation*,
//!   EPL 17(6):479, 1992 (D3Q19).
//! - I. Ginzburg, F. Verhaeghe, D. d'Humières, *Two-relaxation-time lattice Boltzmann scheme:
//!   about parametrization, velocity, pressure and mixed boundary conditions*, Commun. Comput.
//!   Phys. 3(2):427–478, 2008 (TRT and the parameter Λ).
//! - Z. Guo, C. Zheng, B. Shi, *Discrete lattice effects on the forcing term in the lattice
//!   Boltzmann method*, PRE 65:046308, 2002 (body force).
//! - A. J. C. Ladd, *Numerical simulations of particulate suspensions via a discretized Boltzmann
//!   equation. Part 1*, JFM 271:285–309, 1994 (bounce-back, moving walls, momentum exchange).
//! - R. Mei, D. Yu, W. Shyy, L.-S. Luo, *Force evaluation in the lattice Boltzmann method
//!   involving curved geometry*, PRE 65:041203, 2002 (momentum exchange).
//! - Q. Zou, X. He, *On pressure and velocity boundary conditions for the lattice Boltzmann BGK
//!   model*, Phys. Fluids 9(6):1591–1598, 1997; M. Hecht, J. Harting, J. Stat. Mech. 2010,
//!   P01018 (the D3Q19 transverse correction).
//! - Z. Guo, T. S. Zhao, *Lattice Boltzmann model for incompressible flows through porous
//!   media*, PRE 66:036304, 2002 (the implicit velocity with a Darcy–Forchheimer force).
//! - I. Orlanski, *A simple boundary condition for unbounded hyperbolic flows*, J. Comput. Phys.
//!   21(3):251–269, 1976 (the convective outflow condition).
//! - P. Lallemand, L.-S. Luo, *Theory of the lattice Boltzmann method: dispersion, dissipation,
//!   isotropy, Galilean invariance, and stability*, PRE 61:6546, 2000 (stability limits).
//! - H. Edelsbrunner, E. P. Mücke, *Simulation of Simplicity*, ACM TOG 9(1):66–104, 1990;
//!   J. R. Shewchuk 1997 (exact predicates, through `fieldcraft_mesh::predicates`).
//!
//! No code of any lattice Boltzmann package was read; everything is written from these papers and
//! the book.
//!
//! # Determinism
//! Double precision and fixed loop orders. Collision and streaming run on threads by fixed blocks
//! of cells ([`fieldcraft_parallel`]), each cell's update independent of the others; the wall
//! forces and every other sum run in a fixed order on one thread. The same model gives the same
//! result bit for bit on any number of threads and on every platform, wasm included.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]
// Population arrays are indexed by cell and direction counts fixed when they are allocated.
#![allow(clippy::needless_range_loop)]

pub mod backend;
pub mod conduction;
pub mod flow;
pub mod lattice;
pub mod memory;
pub mod outlet;
pub mod regularized;
pub mod setup;
pub mod solve;
pub mod steady;
pub mod turbulence;
pub mod voxel;

pub use backend::{Backend, Choice, Fields, device_available, register_backend};
pub use flow::{Flow, Group, Link, OpenCell, Rule, Target, ThermalState, ZoneSample};
pub use memory::{DeviceBuffer, Memory};
pub use setup::{ContactSetup, SolidSetup, ThermalSetup};
pub use setup::{Obstacle, Setup, ZoneKind, ZoneSetup, fluid_step};
pub use solve::{
    ComponentCheck, FluidResult, FluidStats, ForceHistory, Frame, HeatHistory, Peak, SolidResult, ZoneHistory, run, solve_fluid, solve_fluid_on,
};
pub use steady::{MonitorSample, SteadyReport, SteadySetup, StopReason};
pub use voxel::Grid;

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum Error {
    #[error("invalid input: {0}")]
    Invalid(String),
    #[error("not found: {0}")]
    NotFound(String),
    #[error("limit reached: {0}")]
    Limit(String),
    /// The run diverged.
    #[error("unstable: {0}")]
    Unstable(String),
    #[error(transparent)]
    Model(#[from] fieldcraft_model::Error),
    #[error(transparent)]
    Mesh(#[from] fieldcraft_mesh::Error),
}

pub type Result<T> = std::result::Result<T, Error>;

/// Caps on a run's size and work; a run beyond any of them is refused before it starts. The
/// defaults hold on any machine that solves; `max_cells`, `max_cell_updates` and
/// `max_frame_bytes` may be raised above them up to what `max_memory` holds (the `*_ceiling`
/// methods), and a run with a raised memory cap must fit in `max_memory` ([`Memory::check`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Limits {
    pub max_cells: u64,
    pub max_steps: u64,
    /// Fluid cells times time steps.
    pub max_cell_updates: u64,
    /// The memory the output frames may take.
    pub max_frame_bytes: u64,
    /// The machine's bound on a run's memory in bytes; unbounded by default (an application
    /// sets it from the machine it runs on).
    pub max_memory: u64,
}

impl Limits {
    /// 16 million cells: about 5 GB of populations in double precision.
    pub const MAX_CELLS: u64 = 16_000_000;
    pub const MAX_STEPS: u64 = 10_000_000;
    pub const MAX_CELL_UPDATES: u64 = 1_000_000_000_000;
    pub const MAX_FRAME_BYTES: u64 = 2_000_000_000;
    /// The most cells a cell index addresses.
    pub const MAX_ADDRESSABLE_CELLS: u64 = u32::MAX as u64 - 1;

    /// Whether a memory cap (`max_cells`, `max_frame_bytes`) is above its default: the run then
    /// has to fit `max_memory`. (`max_cell_updates` is work, not memory.)
    pub fn raised(&self) -> bool {
        self.max_cells > Limits::MAX_CELLS || self.max_frame_bytes > Limits::MAX_FRAME_BYTES
    }

    /// Whether `max_memory` bounds anything (an application set it).
    pub fn bounded(&self) -> bool {
        self.max_memory != u64::MAX
    }

    /// The most `max_cells` may be: the cells `max_memory` holds at the fewest bytes a cell
    /// takes, never below the default nor above what a cell index addresses.
    pub fn max_cells_ceiling(&self) -> u64 {
        (self.max_memory / memory::MIN_BYTES_PER_CELL).clamp(Limits::MAX_CELLS, Limits::MAX_ADDRESSABLE_CELLS)
    }

    /// The most `max_cell_updates` may be: the most cells `max_memory` holds for the most steps
    /// (more updates than that cannot happen within the other caps).
    pub fn max_cell_updates_ceiling(&self) -> u64 {
        Limits::MAX_STEPS.saturating_mul(self.max_cells_ceiling()).max(Limits::MAX_CELL_UPDATES)
    }

    /// The most `max_frame_bytes` may be: `max_memory`, never below the default.
    pub fn max_frame_ceiling(&self) -> u64 {
        self.max_memory.max(Limits::MAX_FRAME_BYTES)
    }

    /// Refuse a grid of `cells` above `max_cells`, saying what to change.
    pub fn check_cells(&self, cells: [u64; 3]) -> Result<()> {
        let total = cells[0].saturating_mul(cells[1]).saturating_mul(cells[2]);
        if total > self.max_cells {
            return Err(Error::Limit(format!(
                "the domain has {total} cells ({} × {} × {}), above the cap of {}: coarsen the spacing{}",
                cells[0],
                cells[1],
                cells[2],
                self.max_cells,
                self.raise("max_cells", total, self.max_cells_ceiling())
            )));
        }
        Ok(())
    }

    /// `, or raise `key` to `needed` (up to `ceiling` on this machine)` when that would do,
    /// empty otherwise.
    pub(crate) fn raise(&self, key: &str, needed: u64, ceiling: u64) -> String {
        match needed <= ceiling {
            true if self.bounded() => format!(", or raise `{key}` to {needed} (this machine allows up to {ceiling})"),
            true => format!(", or raise `{key}` to {needed}"),
            false => String::new(),
        }
    }
}

impl Default for Limits {
    fn default() -> Self {
        Limits {
            max_cells: Limits::MAX_CELLS,
            max_steps: Limits::MAX_STEPS,
            max_cell_updates: Limits::MAX_CELL_UPDATES,
            max_frame_bytes: Limits::MAX_FRAME_BYTES,
            max_memory: u64::MAX,
        }
    }
}
