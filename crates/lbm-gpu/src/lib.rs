//! The lattice Boltzmann solver on the GPU: [`GpuFlow`] steps a flow that
//! [`fieldcraft_lbm::Flow`] prepared (grid, links, open cells, initial state) with the compute
//! shaders of `lbm.wgsl`, which are the CPU scheme's collision, streaming rules, open boundaries
//! and momentum exchange, and with heat its temperature lattice, written for the device; `fieldcraft_lbm::run` drives it through the
//! [`Backend`] trait in chunks of steps, so populations stay on the device and only frames,
//! force/zone histories and the stability check come back. Darcy–Forchheimer and loss-coefficient
//! zones use the CPU's implicit velocity correction; fans update their curve feedback before
//! every collision entirely on the device, independently of submission and readback chunks.
//!
//! Single precision: the CPU scheme is the reference, and this backend is deterministic on its
//! device and agrees with it within the tolerances its tests state (`tests/verification.rs`).
//! The device is the one `fieldcraft_gpu`'s adapter policy picks (`FIELDCRAFT_GPU`).
//!
//! Sources: the same as `fieldcraft_lbm` (Krüger et al. 2017; Ginzburg et al. 2008; Guo et al.
//! 2002; Guo & Zhao, PRE 66:036304, 2002; Ladd 1994; Mei et al. 2002; Zou & He 1997;
//! Hecht & Harting 2010; Orlanski 1976). No code of any
//! lattice Boltzmann package was read.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]

#[cfg(not(target_arch = "wasm32"))]
mod native;
#[cfg(not(target_arch = "wasm32"))]
mod zones;
#[cfg(all(test, not(target_arch = "wasm32")))]
mod zones_tests;
#[cfg(not(target_arch = "wasm32"))]
pub use native::{GpuFlow, KernelTime, Options, Streaming, register};

/// In the browser there is no device backend yet (WebGPU later): registering is a no-op and
/// `solve.fluid` with `backend: "gpu"` reports that none is available.
#[cfg(target_arch = "wasm32")]
pub fn register() {}
