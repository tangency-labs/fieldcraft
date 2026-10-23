//! What a run needs in memory, worked out from its setup (or from the model alone, for a form
//! before the run) before anything is allocated, and the checks against the caps and the
//! machine's bound, whose refusals say what to change.
//!
//! The figures follow the arrays the solver allocates, per cell: the momentum populations
//! (19 doubles, and the same again for the streaming targets on the CPU), the cell flags, the
//! temperature lattice with heat, the output frames (single precision: velocity and pressure,
//! and temperature with heat) and one frame's fields while it is written. Links at the box faces
//! are counted; links at obstacles depend on their surface and are not, so a run with obstacles
//! takes a little more (the tests state the margin).

use fieldcraft_model::{FaceCondition, FluidDomain, Model, Procedure};

use crate::backend::Choice;
use crate::flow::{Link, OpenCell};
use crate::lattice::Q;
use crate::setup::{Setup, fluid_step};
use crate::{Error, Limits, Result};

/// Bytes in a megabyte and a gigabyte, for messages.
const MEGABYTE: f64 = 1e6; // unit-factors: bytes of memory, not a physical quantity
const GIGABYTE: f64 = 1e9; // unit-factors: bytes of memory, not a physical quantity

/// Directions of the temperature lattice (D3Q7).
const QT: u64 = 7;
const F64: u64 = 8;
const F32: u64 = 4;

/// Per cell while a flow is prepared: its momentum populations, then its fluid flag, owner,
/// obstacle and open-cell marks.
const POPULATIONS: u64 = Q as u64 * F64;
const FLAGS: u64 = 1 + 4 + 1 + 1;
/// Per cell with heat: the temperature populations, the temperature and the velocity the scalar
/// equilibrium uses; with buoyancy also the cell's body force.
const HEAT: u64 = QT * F64 + F64 + 3 * F64;
const BUOYANCY: u64 = 3 * F64;
/// Per cell with fans or porous zones, and with conducting obstacles: an index.
const INDEX: u64 = 4;
/// Links per cell of a box face: the directions that cross it (5 of D3Q19).
const FACE_LINKS: u64 = 5;
/// Per cell of an output frame: velocity and pressure, and temperature with heat.
const FRAME: u64 = 4 * F32;
const FRAME_HEAT: u64 = F32;
/// Per cell of a frame being written: density and velocity in double precision (and the
/// temperature), and on the GPU the single-precision moments (and temperature populations) read
/// back first.
const READBACK: u64 = 4 * F64;
const READBACK_HEAT: u64 = F64;
const DEVICE_READBACK: u64 = 4 * F32;
const DEVICE_READBACK_HEAT: u64 = QT * F32;
/// Per cell of the result's solid flags, and of the GPU backend's host copy of them.
const RESULT: u64 = 1;
/// Per cell the GPU backend stages on the host while it uploads the flow: the flags twice and
/// one direction of populations.
const DEVICE_STAGING: u64 = 1 + 4 + F32;
/// Per cell of the GPU backend's population buffer (one buffer, single precision).
pub const DEVICE_POPULATIONS: u64 = Q as u64 * F32;

/// Words in a GPU fan/porous header, appended to its flags buffer (lbm-gpu/zones.rs).
pub const DEVICE_ZONE_WORDS: usize = 16;

/// The fewest bytes a cell takes on the host (a GPU run of a flow without heat while it is
/// prepared): what bounds the cells a machine's memory can hold.
pub const MIN_BYTES_PER_CELL: u64 = POPULATIONS + FLAGS + DEVICE_STAGING;

/// What a run needs in memory, in bytes (host memory unless named otherwise).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Memory {
    pub backend: Choice,
    pub cells: u64,
    /// The frames the run keeps, the initial one included.
    pub frames: u64,
    /// While the flow is prepared and handed to its backend.
    pub setup: u64,
    /// The flow's state while it steps.
    pub state: u64,
    /// The output frames.
    pub frame_bytes: u64,
    /// One frame's fields while it is written.
    pub readback: u64,
    /// On the GPU: the population buffer, and the most one buffer may hold on the device (when
    /// the backend reports it).
    pub device: Option<DeviceBuffer>,
}

/// The GPU backend's population buffer against the device's limit on one buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DeviceBuffer {
    pub populations: u64,
    pub limit: Option<u64>,
}

/// What sizes a run's memory.
#[derive(Debug, Clone, Copy)]
struct Shape {
    cells: [u64; 3],
    frames: u64,
    heat: bool,
    buoyancy: bool,
    zones: bool,
    zone_data: u64,
    zone_count: u64,
    conducting: bool,
    open_faces: [bool; 6],
}

impl Shape {
    /// The shape of a run over `domain` keeping `frames` frames after the initial one.
    fn new(domain: &FluidDomain, frames: u64) -> Result<Shape> {
        let thermal = domain.thermal.as_ref();
        Ok(Shape {
            cells: domain.cells()?,
            frames: frames.max(1) + 1,
            heat: thermal.is_some(),
            // The setup's buoyancy is −β g: nonzero with both.
            buoyancy: thermal.is_some_and(|t| t.gravity != [0.0; 3] && t.expansion.is_some_and(|b| b != 0.0)),
            zones: !domain.fans.is_empty() || !domain.porous.is_empty(),
            zone_count: (domain.fans.len() + domain.porous.len()) as u64,
            zone_data: ((domain.fans.len() + domain.porous.len()) as u64 * DEVICE_ZONE_WORDS as u64
                + domain.fans.iter().map(|f| 2 * f.curve.len() as u64).sum::<u64>())
                * F32,
            conducting: thermal.is_some_and(|t| !t.solids.is_empty()),
            open_faces: fieldcraft_model::BoxFace::ALL.map(|f| is_open(domain.faces.get(f))),
        })
    }
}

impl Memory {
    /// The memory of the run `setup` prepares on `backend`.
    pub fn of(setup: &Setup, backend: Choice) -> Result<Memory> {
        Ok(Memory::from_shape(&Shape::new(&setup.domain, u64::from(setup.frames).min(setup.steps))?, backend))
    }

    /// The memory of fluid step `step` of `model` on `backend`, from the model alone (no
    /// obstacle is triangulated): what a form shows before the run. The same as [`Memory::of`]
    /// but for a step with fewer time steps than frames.
    pub fn estimate(model: &Model, step: Option<&str>, backend: Choice) -> Result<Memory> {
        let (s, domain) = fluid_step(model, step)?;
        let Procedure::Fluid { frames, .. } = &s.procedure else {
            return Err(Error::Invalid(format!("step `{}` is not a fluid step", s.name)));
        };
        Ok(Memory::from_shape(&Shape::new(domain, u64::from(*frames))?, backend))
    }

    fn from_shape(s: &Shape, backend: Choice) -> Memory {
        let [nx, ny, nz] = s.cells;
        let n = nx.saturating_mul(ny).saturating_mul(nz);
        let per_cell = |bytes: u64| n.saturating_mul(bytes);
        // The cells of each box face, and the links and open cells they hold.
        let face_cells = [ny * nz, ny * nz, nx * nz, nx * nz, nx * ny, nx * ny].map(|c| c.min(n));
        let boundary: u64 = face_cells
            .iter()
            .zip(s.open_faces)
            .map(|(c, open)| {
                c.saturating_mul(
                    FACE_LINKS * size_of::<Link>() as u64
                        + if open { (size_of::<OpenCell>() + size_of::<[f64; crate::lattice::Q]>()) as u64 } else { 0 },
                )
            })
            .fold(0, u64::saturating_add);
        let mut flow = POPULATIONS + FLAGS;
        let mut targets = POPULATIONS;
        if s.heat {
            flow += HEAT;
            targets += QT * F64;
        }
        if s.buoyancy {
            flow += BUOYANCY;
        }
        if s.zones {
            // Both ownership per cell and the zones' lists of owned fluid cells.
            flow += 2 * INDEX;
        }
        if s.conducting {
            flow += INDEX;
        }
        let prepared = per_cell(flow).saturating_add(boundary);
        let (setup, state, readback) = match backend {
            // The CPU flow allocates its streaming targets on the first step and keeps
            // everything to the end.
            Choice::Cpu => {
                let state = prepared.saturating_add(per_cell(targets)).saturating_add(per_cell(RESULT));
                (state, state, per_cell(READBACK + if s.heat { READBACK_HEAT } else { 0 }))
            }
            // The device backend copies the flow up, and the flow is dropped before the steps.
            Choice::Device => (
                prepared.saturating_add(per_cell(DEVICE_STAGING)).saturating_add(s.zone_data.saturating_mul(3)),
                per_cell(2 * RESULT).saturating_add(s.zone_data),
                per_cell(READBACK + DEVICE_READBACK + if s.heat { READBACK_HEAT + DEVICE_READBACK_HEAT } else { 0 })
                    .saturating_add(s.zone_count.saturating_mul(4 * F32)),
            ),
        };
        let device = (backend == Choice::Device)
            .then(|| DeviceBuffer { populations: per_cell(DEVICE_POPULATIONS), limit: crate::backend::device_buffer_limit() });
        Memory {
            backend,
            cells: n,
            frames: s.frames,
            setup,
            state,
            frame_bytes: per_cell(FRAME + if s.heat { FRAME_HEAT } else { 0 }).saturating_mul(s.frames),
            readback,
            device,
        }
    }

    /// The most the run holds at once.
    pub fn total(&self) -> u64 {
        self.setup.max(self.state.saturating_add(self.frame_bytes).saturating_add(self.readback))
    }

    /// What the flow holds at its largest (while it is prepared or while it steps).
    pub fn flow_bytes(&self) -> u64 {
        self.state.max(self.setup)
    }

    /// The most frames (the initial one included) `bytes` hold.
    fn frames_in(&self, bytes: u64) -> u64 {
        bytes / (self.frame_bytes / self.frames.max(1)).max(1)
    }

    /// What the frames may take beside the flow within the machine's bound.
    fn frame_room(&self, limits: &Limits) -> u64 {
        limits.max_memory.saturating_sub(self.state.saturating_add(self.readback))
    }

    /// The most frames (the initial one included) that fit `limits`: the frame cap and, with a
    /// cap raised above its default, the machine's bound.
    pub fn frames_that_fit(&self, limits: &Limits) -> u64 {
        let mut fit = self.frames_in(limits.max_frame_bytes);
        if limits.raised() {
            fit = fit.min(self.frames_in(self.frame_room(limits)));
        }
        fit
    }

    /// Check the run against `limits`, before anything is allocated: the GPU's buffer, the frame
    /// cap and, when a cap is raised above its default, the machine's bound. Each refusal says
    /// what to change.
    pub fn check(&self, limits: &Limits) -> Result<()> {
        if let Some(DeviceBuffer { populations, limit: Some(limit) }) = self.device
            && populations > limit
        {
            return Err(Error::Limit(format!(
                "{} cells need {} of populations in one GPU buffer; this GPU allows {} per buffer, which holds {} cells: coarsen the spacing",
                self.cells,
                megabytes(populations),
                megabytes(limit),
                limit / DEVICE_POPULATIONS
            )));
        }
        self.check_frames(limits)?;
        if limits.raised() && self.total() > limits.max_memory {
            return Err(Error::Limit(format!(
                "the run needs {} ({} of flow state, {} for {} frames), above this machine's bound of {} for runs with raised caps: {}",
                megabytes(self.total()),
                megabytes(self.flow_bytes()),
                megabytes(self.frame_bytes),
                self.frames,
                megabytes(limits.max_memory),
                fewer_frames(self.frames_that_fit(limits))
            )));
        }
        Ok(())
    }

    /// Check the frames against the frame cap.
    pub fn check_frames(&self, limits: &Limits) -> Result<()> {
        if self.frame_bytes > limits.max_frame_bytes {
            // What the frames of this run may take on this machine, beside its flow.
            let room = self.frame_room(limits).clamp(Limits::MAX_FRAME_BYTES, limits.max_frame_ceiling());
            let raise = limits.raise("max_frame_bytes", self.frame_bytes, room);
            return Err(Error::Limit(format!(
                "{} frames of {} cells need {}, above the cap of {}: {}{raise}",
                self.frames,
                self.cells,
                megabytes(self.frame_bytes),
                megabytes(limits.max_frame_bytes),
                fewer_frames(self.frames_that_fit(limits))
            )));
        }
        Ok(())
    }
}

/// What to say when only `fit` frames (the initial one included) fit.
fn fewer_frames(fit: u64) -> String {
    if fit >= 2 {
        format!("{fit} frames fit: set the step's `frames` to at most {}", fit - 1)
    } else {
        "not even two frames fit: coarsen the spacing".into()
    }
}

/// Bytes for a message: megabytes below 10 GB, gigabytes above.
pub fn megabytes(bytes: u64) -> String {
    let b = bytes as f64;
    if b < 10.0 * GIGABYTE { format!("{:.0} MB", b / MEGABYTE) } else { format!("{:.1} GB", b / GIGABYTE) }
}

/// Whether a box face is an inlet or an outlet.
pub(crate) fn is_open(face: &FaceCondition) -> bool {
    matches!(face, FaceCondition::Velocity { .. } | FaceCondition::Pressure { .. })
}
