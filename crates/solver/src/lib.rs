//! FieldCraft's implicit solver: linear static and frequency analysis of a model's steps.
//!
//! The pipeline, for one step of a [`fieldcraft_model::Model`]:
//! 1. **System** ([`system`]): the model is resolved; every section's elements become solid
//!    elements with their material's elasticity and the section's formulation
//!    ([`fieldcraft_fem`] kernels); the step's boundary conditions give the prescribed degrees of
//!    freedom, its loads the load vector (forces and tractions as the model distributes them,
//!    consistent pressure loads, gravity and body forces).
//! 2. **Assembly** ([`sparse`]) on the node graph, elements in label order.
//! 3. **Ordering** ([`ordering`]): nested dissection of the nodes; the free degrees of freedom
//!    are numbered in that order and the prescribed ones eliminated (u_p moved to the right-hand
//!    side), with no penalty numbers.
//! 4. **Solution** ([`cholesky`]): sparse Cholesky; reactions K u − f at the prescribed degrees
//!    of freedom; stresses at the integration points, extrapolated and averaged at the nodes, and
//!    von Mises. Frequency analysis ([`eigen`]) runs shift-invert block Lanczos on the same
//!    kind of factorization. Indefinite matrices (tangents past a limit point, strained mixed
//!    elements, shifted operators) are factored by the symmetric indefinite LDLᵀ ([`ldlt`]) on
//!    the same symbolic factorization, which also counts their negative eigenvalues.
//!
//! Element matrices and stresses, the factorization's subtrees and dense tiles, the Lanczos
//! block solves and matrix products run on threads ([`fieldcraft_parallel`]), but every number
//! is computed by one task in an order fixed by the model alone and sums are taken in a fixed
//! order, so the same model gives the same numbers bit for bit on any number of threads and on
//! every platform; sizes and iterations are capped by [`Limits`].
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]
// Matrix kernels read as the textbook formulas with index loops.
#![allow(clippy::needless_range_loop)]

pub mod analysis;
pub mod buckle;
pub mod cholesky;
pub mod complex;
pub mod constraints;
pub mod contact;
pub mod discrete;
pub mod eigen;
pub mod foundation;
pub mod harmonic;
pub mod heat;
pub mod ldlt;
pub mod nonlinear;
pub mod ordering;
pub mod random;
pub mod sensitivity;
pub mod sparse;
pub mod spectrum;
pub mod system;
pub mod tie;
pub mod transient;

pub use analysis::{
    ElementStress, FrequencyOptions, FrequencyResult, Location, MassKind, Mode, PointResult, StaticResult, Stats, StepResult, solve_frequency,
    solve_static,
};
pub use buckle::{Buckling, solve_buckle};
pub use harmonic::{HarmonicResult, solve_harmonic};
pub use heat::{HeatResult, solve_heat};
pub use nonlinear::solve_nonlinear;
pub use random::{
    MilesResult, RandomFieldOutput, RandomInput, RandomResponse, distributed_force_modal_psd, integrate_psd, miles_base_acceleration_rms,
    modal_base_response, modal_coordinate_psd, modal_stress_field_psd, modal_stress_rms, sdof_base_acceleration_psd, sdof_force_bandpass_variance,
    von_mises_rms_from_components,
};

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum Error {
    #[error(transparent)]
    Model(#[from] fieldcraft_model::Error),
    #[error("invalid input: {0}")]
    Invalid(String),
    #[error("not supported yet: {0}")]
    NotSupported(String),
    /// An element kernel failed (inverted geometry, a formulation it lacks, …).
    #[error("element {label}: {source}")]
    Element { label: u32, source: fieldcraft_fem::Error },
    /// The stiffness (or shifted matrix) is singular.
    #[error("{0}")]
    Singular(String),
    #[error("limit reached: {0}")]
    Limit(String),
    #[error("did not converge: {0}")]
    NotConverged(String),
    /// A broken internal invariant (a bug), reported instead of a panic.
    #[error("internal error: {0}")]
    Internal(String),
}

pub type Result<T> = std::result::Result<T, Error>;

impl Error {
    /// Factor-value bytes alone; indices and numeric workspace need additional memory.
    pub(crate) fn factor_limit(entries: usize, limit: usize) -> Error {
        let bytes = entries as u128 * std::mem::size_of::<f64>() as u128;
        Error::Limit(format!(
            "the factorization needs {entries} entries ({bytes} bytes of factor values alone; indices and workspace are additional); the limit is {limit} entries (refine less or raise max_factor_entries)"
        ))
    }

    /// A kernel error, naming the element.
    pub fn element(label: u32, source: fieldcraft_fem::Error) -> Error {
        Error::Element { label, source }
    }
}

/// Caps on the size of a solution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Limits {
    /// Degrees of freedom of the mesh.
    pub max_dofs: usize,
    /// Factor storage in 8-byte real slots (two per native complex entry);
    /// also bounds the Lanczos basis and stored laminate ply-field memory (each separately).
    pub max_factor_entries: usize,
    /// Lanczos vectors (iterations × block size).
    pub max_lanczos_vectors: usize,
    /// Newton iterations of a nonlinear step, over all its increments.
    pub max_newton_iterations: usize,
    /// Optional whole-factor byte budget. This is checked after symbolic analysis and before
    /// numeric factor storage; the caller must keep separate budgets for matrices/results.
    pub max_factor_memory_bytes: Option<usize>,
}

impl Limits {
    pub const MAX_DOFS: usize = 6_000_000;
    pub const MAX_FACTOR_ENTRIES: usize = 400_000_000;
    pub const MAX_LANCZOS_VECTORS: usize = 3_000;
    pub const MAX_NEWTON_ITERATIONS: usize = 100_000;
}

impl Default for Limits {
    fn default() -> Limits {
        Limits {
            max_dofs: Limits::MAX_DOFS,
            max_factor_entries: Limits::MAX_FACTOR_ENTRIES,
            max_lanczos_vectors: Limits::MAX_LANCZOS_VECTORS,
            max_newton_iterations: Limits::MAX_NEWTON_ITERATIONS,
            max_factor_memory_bytes: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factor_limit_estimates_bytes_without_allocating_the_factor() {
        for (entries, bytes) in [(467_000_000, "3736000000"), (2_000_000_000, "16000000000")] {
            let e = Error::factor_limit(entries, 400_000_000).to_string();
            assert!(e.contains(&format!("{entries} entries ({bytes} bytes of factor values alone")), "{e}");
            assert!(e.contains("indices and workspace are additional"), "{e}");
            assert!(e.contains("limit is 400000000 entries"), "{e}");
            assert!(e.contains("raise max_factor_entries"), "{e}");
        }
    }
}
