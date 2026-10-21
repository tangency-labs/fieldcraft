//! Test tooling (a dev-dependency only, `cargo xtask layers`).
//!
//! [`corpus`] is the meshing test corpus: parts built with SolidCraft from recipes committed in
//! `tests/corpus/<name>/recipe.json`, and analytic shapes generated in code with closed-form
//! volumes. Meshing tests load their geometry from here; [`runner`] (the `corpus` binary,
//! `cargo xtask corpus`) builds and meshes every part, checks each part against its frozen
//! results ([`fixture`]), and writes docs/meshing/corpus.md.
//!
//! [`skill`] (the `skill` binary, `cargo xtask skill-results`) meshes every corpus part and
//! sample design by following the meshing skill and writes docs/meshing/skill-results.md.
//!
//! [`oracle`] (the `oracle` binary, `cargo xtask oracle`) runs the verification and validation
//! cases of `tests/oracle` against their references and fixtures and writes docs/oracle.md.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]

pub mod corpus;
pub mod fixture;
pub mod oracle;
pub mod runner;
pub mod skill;
