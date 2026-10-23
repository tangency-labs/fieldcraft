//! What steps a prepared flow: the CPU scheme in [`crate::flow::Flow`], or a device backend an
//! application registers (a GPU crate above this layer implements the trait and hands a factory
//! to [`register_backend`]; `solve.fluid` picks it when asked). Every backend takes its grid,
//! links and open cells from the [`Flow`] that prepared them, so voxelisation and every boundary
//! rule exist once, here.
//!
//! The CPU result is the reference: a device backend is deterministic on its device and agrees
//! with the CPU within a stated tolerance (its tests say which).

use std::sync::OnceLock;

use crate::flow::{Flow, Group};
use crate::setup::Setup;
use crate::{Error, Result};

/// The density and velocity of every cell (lattice units; zero in solid cells), and the
/// temperature (model units) when the flow carries heat (empty otherwise).
#[derive(Debug, Clone, PartialEq)]
pub struct Fields {
    pub density: Vec<f64>,
    pub velocity: Vec<[f64; 3]>,
    pub temperature: Vec<f64>,
}

/// A stepper of one prepared flow.
pub trait Backend {
    /// What steps the flow, for reports: `cpu`, `gpu: <adapter>`.
    fn name(&self) -> String;
    /// Advance `steps` time steps. Fails with [`Error::Unstable`] when a density stops being a
    /// positive number or a lattice speed exceeds [`crate::setup::MAX_LATTICE_SPEED`] (a device
    /// backend may notice at the end of the chunk rather than the step).
    fn advance(&mut self, steps: u64) -> Result<()>;
    fn steps_done(&self) -> u64;
    /// The density and velocity (with the half force) of every cell.
    fn fields(&mut self) -> Result<Fields>;
    /// The total mass of the fluid cells (the sum of every population).
    fn mass(&mut self) -> Result<f64>;
    /// The momentum each force group received in the last step (lattice units: force per step).
    fn forces(&mut self) -> Result<Vec<[f64; 3]>>;
    /// The heat each force group received in the last step (lattice units: temperature times
    /// cell faces per step); empty without heat.
    fn heat(&mut self) -> Result<Vec<f64>>;
    /// Each fan's and porous zone's flow and pressure now (lattice units). A backend that does
    /// not step zones refuses a flow with any when it is made, so the default has none.
    fn zones(&mut self) -> Result<Vec<crate::flow::ZoneSample>> {
        Ok(Vec::new())
    }
    /// The largest eddy viscosity over the fluid's own viscosity in the last step. A backend
    /// without a turbulence model refuses a flow with one when it is made, so the default is 0.
    fn eddy_viscosity_ratio(&self) -> f64 {
        0.0
    }
    fn groups(&self) -> &[Group];
    fn fluid_cells(&self) -> usize;
    fn is_fluid(&self, cell: usize) -> bool;
    /// The cells of each obstacle that conducts heat (in the setup's order of solids), in cell
    /// order; a backend that does not conduct refuses such a flow when it is made and has none.
    fn conducting(&self) -> &[Vec<u32>];
    /// Obstacles that cover no cell centre.
    fn empty_obstacles(&self) -> &[String];
}

/// Makes a device backend from a prepared flow (its initial state included) and its setup.
pub type BackendFactory = dyn Fn(&Flow, &Setup) -> Result<Box<dyn Backend>> + Send + Sync;

/// The most bytes one buffer may hold on the device the backend would use (`None` when there is
/// no device or it cannot tell): what [`crate::Memory`] checks the population buffer against
/// before a flow is prepared.
pub type BufferLimit = dyn Fn() -> Option<u64> + Send + Sync;

struct Registered {
    make: Box<BackendFactory>,
    buffer_limit: Box<BufferLimit>,
}

static FACTORY: OnceLock<Registered> = OnceLock::new();

/// Register the process's device backend and its device's buffer limit (an application does
/// this once at start; a second registration is ignored).
pub fn register_backend(factory: Box<BackendFactory>, buffer_limit: Box<BufferLimit>) {
    let _ = FACTORY.set(Registered { make: factory, buffer_limit });
}

/// The registered device backend's factory, if an application registered one.
pub fn backend_factory() -> Option<&'static BackendFactory> {
    FACTORY.get().map(|b| b.make.as_ref())
}

/// The registered device's limit on one buffer, when there is one and it tells.
pub fn device_buffer_limit() -> Option<u64> {
    FACTORY.get().and_then(|b| (b.buffer_limit)())
}

/// Whether this machine has a device backend to solve on: one registered, with a device the
/// `FIELDCRAFT_GPU` policy picks; else why not.
pub fn device_available() -> std::result::Result<(), String> {
    if backend_factory().is_none() {
        return Err("this build has no GPU backend".into());
    }
    match device_buffer_limit() {
        Some(_) => Ok(()),
        None => Err("no GPU the FIELDCRAFT_GPU policy accepts was found".into()),
    }
}

/// Which backend a solve asks for.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Choice {
    /// The CPU scheme: the reference, bit for bit reproducible.
    #[default]
    Cpu,
    /// The registered device backend; an error when none is registered.
    Device,
}

impl Choice {
    /// `cpu` or `gpu`, exactly.
    pub fn parse(s: &str) -> Option<Choice> {
        match s {
            "cpu" => Some(Choice::Cpu),
            "gpu" => Some(Choice::Device),
            _ => None,
        }
    }
}

/// The backend `choice` names for `flow`.
pub fn backend(choice: Choice, flow: Flow, setup: &Setup) -> Result<Box<dyn Backend>> {
    match choice {
        Choice::Cpu => Ok(Box::new(flow)),
        Choice::Device => match backend_factory() {
            Some(make) => make(&flow, setup),
            None => Err(Error::Invalid("no GPU backend is available in this build: solve on the CPU (`backend: \"cpu\"`)".into())),
        },
    }
}
