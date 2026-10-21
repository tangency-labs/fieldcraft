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
pub mod outlet;
pub mod regularized;
pub mod setup;
pub mod solve;
pub mod steady;
pub mod voxel;

pub use backend::{Backend, Choice, Fields, register_backend};
pub use flow::{Flow, Group, Link, OpenCell, Rule, Target, ThermalState, ZoneSample};
pub use setup::{Obstacle, Setup, ZoneKind, ZoneSetup, fluid_step};
pub use setup::{SolidSetup, ThermalSetup};
pub use solve::{FluidResult, FluidStats, ForceHistory, Frame, HeatHistory, SolidResult, ZoneHistory, run, solve_fluid, solve_fluid_on};
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

/// Caps on a run's size and work; a run beyond any of them is refused before it starts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Limits {
    pub max_cells: u64,
    pub max_steps: u64,
    /// Fluid cells times time steps.
    pub max_cell_updates: u64,
    /// The memory the output frames may take.
    pub max_frame_bytes: u64,
}

impl Limits {
    /// 16 million cells: about 5 GB of populations in double precision.
    pub const MAX_CELLS: u64 = 16_000_000;
    pub const MAX_STEPS: u64 = 10_000_000;
    pub const MAX_CELL_UPDATES: u64 = 1_000_000_000_000;
    pub const MAX_FRAME_BYTES: u64 = 2_000_000_000;

    /// A device backend's caps. The device keeps the populations in single precision and checks
    /// its own buffer limits; the host still builds the CPU flow first (about 440 bytes per cell
    /// at its peak while the device's is made), so 200 million cells take about 88 GB, which a
    /// 128 GB unified-memory machine holds. Steady runs at that size take some 10⁵ steps, and a
    /// frame of every cell about 20 bytes per cell.
    pub const GPU_MAX_CELLS: u64 = 200_000_000;
    pub const GPU_MAX_CELL_UPDATES: u64 = 100_000_000_000_000;
    pub const GPU_MAX_FRAME_BYTES: u64 = 8_000_000_000;

    /// The default caps for a solve on `choice`.
    pub fn for_backend(choice: Choice) -> Limits {
        match choice {
            Choice::Cpu => Limits::default(),
            Choice::Device => Limits {
                max_cells: Limits::GPU_MAX_CELLS,
                max_cell_updates: Limits::GPU_MAX_CELL_UPDATES,
                max_frame_bytes: Limits::GPU_MAX_FRAME_BYTES,
                ..Limits::default()
            },
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
        }
    }
}
