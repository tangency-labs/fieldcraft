//! The device backend (native targets).

use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

use fieldcraft_gpu::select::{self, Policy};
use fieldcraft_lbm::backend::{Backend, Fields};
use fieldcraft_lbm::lattice::{C, Q, W};
use fieldcraft_lbm::outlet::UNKNOWN;
use fieldcraft_lbm::setup::MAX_LATTICE_SPEED;
use fieldcraft_lbm::{Error, Flow, Group, Result, Rule, Target};
use fieldcraft_model::{BoxFace, ThermalCondition};
use wgpu::util::DeviceExt;

const SHADER: &str = include_str!("lbm.wgsl");
const WG: u32 = 256;
/// The widest row of workgroups. A dispatch wider than one row is exactly this wide ([`groups`]),
/// so the shader flattens every dispatch's workgroups with this one row width.
const MAX_GROUPS_X: u32 = 32768;
const KIND_PULL: u32 = 0;
const KIND_BOUNCE: u32 = 1;
const KIND_REFLECT: u32 = 2;
const KIND_OPEN: u32 = 3;
/// Marks a link into a cell whose scalar populations the open boundary sets (see lbm.wgsl).
const RESET: u32 = 1 << 24;
/// Marks a bounce link that interpolates where the obstacle's surface crosses it (see lbm.wgsl).
const INTERP: u32 = 1 << 25;
/// Floats per open cell in the outlet buffer (see lbm.wgsl).
const OUTLET: usize = fieldcraft_lbm::memory::DEVICE_OUTLET_WORDS;
/// The flag word of a boundary node of the on-node wall-function bounce: it streams, but it is not
/// fluid (lbm.wgsl, `is_fluid`).
const NODE: u32 = 2;
/// The floats of one wall model wall cell in `walls`: four rows (lbm.wgsl).
const WALL_ROW: usize = 16;
/// The scalar lattice's directions (the first seven of D3Q19) and weights.
const QT: usize = 7;
const WT: [f64; QT] = [0.25, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125];

/// Check before staging or creating the pressure buffer: in a two-cell-thick box it can be
/// larger than the population buffer, whose device-limit check is therefore insufficient.
fn outlet_words(cells: usize, cap: u64) -> Result<usize> {
    let words = cells.checked_mul(OUTLET).ok_or_else(|| Error::Limit("too many GPU outlet cells".into()))?;
    let bytes = (words as u64).saturating_mul(4);
    if bytes > cap {
        return Err(Error::Limit(format!("GPU pressure outlets need {bytes} bytes; this device allows {cap} per buffer")));
    }
    Ok(words)
}

fn gpu_error(e: impl std::fmt::Display) -> Error {
    update_diagnostic_phase(None);
    Error::Backend(format!("GPU: {e}"))
}

fn diagnostic_limits(limits: &wgpu::Limits) -> fieldcraft_model::progress::BackendLimits {
    fieldcraft_model::progress::BackendLimits {
        max_storage_buffers_per_shader_stage: limits.max_storage_buffers_per_shader_stage,
        max_buffer_size: limits.max_buffer_size,
        max_storage_buffer_binding_size: limits.max_storage_buffer_binding_size,
    }
}

/// Setup/failure boundaries only. A reporter belongs to one command, so concurrent solves
/// cannot overwrite one another's adapter details. Process counter scopes stay explicit.
fn update_diagnostic_phase(phase: Option<&str>) {
    if let Some(progress) = fieldcraft_model::progress::current()
        && let Some(mut context) = progress.backend_diagnostic()
    {
        if let Some(phase) = phase {
            context.phase = phase.to_string();
        }
        context.devices_opened = Some(DEVICES_OPENED.load(Ordering::Relaxed));
        context.instances_opened = Some(select::instances_opened());
        progress.set_backend_diagnostic(Some(context));
    }
}

/// Steps per submission: the stability check and a cancel take effect between submissions.
const MAX_CHUNK: u64 = 64;

/// How long one wait for the device may take before the solve gives up on it.
const DEVICE_TIMEOUT: std::time::Duration = std::time::Duration::from_secs(120);

/// How long a readback waits for its mapping callback before it polls the device again.
const CALLBACK_WAIT: std::time::Duration = std::time::Duration::from_millis(10);

/// Seconds in a nanosecond: GPU timestamp periods are in nanoseconds (wgpu).
const NANOSECOND: f64 = 1e-9; // unit-factors: wall-clock profiling time, not a model quantity

/// A bounded chunk's dispatches, shared by the normal and timestamped paths.
type Dispatch = (usize, &'static str, (u32, u32));

/// A vector of `len` copies of `value`, or [`Error::Limit`] when memory runs out (the host
/// arrays are sized by the grid: 16 M cells is over a gigabyte of populations).
fn filled<T: Clone>(len: usize, value: T, what: &str) -> Result<Vec<T>> {
    let mut v = Vec::new();
    v.try_reserve_exact(len).map_err(|_| Error::Limit(format!("not enough memory for the {what} ({len} entries)")))?;
    v.resize(len, value);
    Ok(v)
}

/// An empty vector with room for `len`, or [`Error::Limit`].
fn reserved<T>(len: usize, what: &str) -> Result<Vec<T>> {
    let mut v = Vec::new();
    v.try_reserve_exact(len).map_err(|_| Error::Limit(format!("not enough memory for the {what} ({len} entries)")))?;
    Ok(v)
}

/// Workgroup counts (x, y) covering `count` items: one row of up to [`MAX_GROUPS_X`], or rows of
/// exactly that many.
fn groups(count: usize) -> Result<(u32, u32)> {
    let total = u32::try_from(count.div_ceil(WG as usize).max(1)).map_err(|_| Error::Limit("too many cells for one dispatch".into()))?;
    let gx = total.min(MAX_GROUPS_X);
    let gy = total.div_ceil(gx);
    if gy > 65535 {
        return Err(Error::Limit("too many cells for one dispatch".into()));
    }
    Ok((gx, gy))
}

fn bytes_u32(v: &[u32]) -> Vec<u8> {
    v.iter().flat_map(|x| x.to_le_bytes()).collect()
}

fn bytes_f32(v: &[f32]) -> Vec<u8> {
    v.iter().flat_map(|x| x.to_le_bytes()).collect()
}

/// Register this backend with `fieldcraft_lbm` (an application does this once at start), so
/// `solve.fluid` with `backend: "gpu"` steps on the device `FIELDCRAFT_GPU` picks.
pub fn register() {
    fieldcraft_lbm::register_backend(
        Box::new(|flow, _setup| GpuFlow::from_env(flow).map(|g| Box::new(g) as Box<dyn Backend>)),
        Box::new(buffer_limit),
    );
}

/// The most bytes one storage buffer may hold under `limits`.
fn buffer_cap(limits: &wgpu::Limits) -> u64 {
    limits.max_storage_buffer_binding_size.min(limits.max_buffer_size)
}

/// The most bytes one storage buffer may hold on the adapter `FIELDCRAFT_GPU` picks (the
/// device is requested with the adapter's own limits, see [`GpuFlow::with_options`]): asked once
/// per process, since enumerating adapters is slow.
fn buffer_limit() -> Option<u64> {
    static LIMIT: std::sync::OnceLock<Option<u64>> = std::sync::OnceLock::new();
    *LIMIT.get_or_init(|| {
        let policy = Policy::from_env().ok()?;
        let (adapter, _) = select::adapter(policy).ok()??;
        Some(buffer_cap(&adapter.limits()))
    })
}

/// How the momentum populations stream (see `lbm.wgsl`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Streaming {
    /// In place: measured faster on every device tried (RTX 4090, GB10, M4 Max; #563).
    #[default]
    Auto,
    /// One population buffer, collided and streamed in place (the AA pattern).
    InPlace,
    /// Two population buffers: collide in place, then pull into the other (A–B).
    TwoBuffers,
}

/// How a [`GpuFlow`] runs; the defaults are what a solve uses.
#[derive(Clone, Copy, Debug, Default)]
pub struct Options {
    pub streaming: Streaming,
    /// Time every kernel with timestamp queries ([`GpuFlow::kernel_times`]). Each dispatch then
    /// gets its own pass and every chunk waits for its timestamps: for measurements, not solves.
    pub profile: bool,
}

/// GPU time per kernel, summed over the dispatches since the flow was created (with
/// [`Options::profile`]).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct KernelTime {
    pub dispatches: u64,
    pub seconds: f64,
}

/// A prepared flow on the device.
pub struct GpuFlow {
    device: wgpu::Device,
    queue: wgpu::Queue,
    /// With [`Options::profile`]: the GPU time of each kernel so far.
    kernel_times: Option<std::cell::RefCell<std::collections::BTreeMap<&'static str, KernelTime>>>,
    /// How long compiling the kernels took for this flow (zero when the process had compiled
    /// them on this device before).
    pipeline_seconds: f64,
    aa: bool,
    adapter: String,
    /// The kernels, shared with every flow on this device that uses the same variant.
    kernels: Arc<Kernels>,
    /// Each kernel's two bind groups: AB swaps populations; AA keeps one population buffer and
    /// changes its orientation. Both modes swap the separate scalar lattice buffers.
    bind: HashMap<&'static str, [wgpu::BindGroup; 2]>,
    /// Which bind group's `f` holds the state.
    current: usize,
    n: usize,
    cells: (u32, u32),
    link_count: usize,
    links_groups: (u32, u32),
    /// The wall model's sample dispatch, when there are modelled wall cells, and its buffer.
    wall_groups: Option<(u32, u32)>,
    /// Whether the wall rows are the halfway wall-function bounce's (`bounce_sample`).
    wall_bounce: bool,
    /// Whether they are the on-node bounce's boundary nodes (`bounce_sample`, `node_bounce`).
    on_node: bool,
    /// Whether the collision takes the finite-difference velocity gradient (the hybrid collision,
    /// Vreman's model): `gradient_velocity` runs before it.
    gradients: bool,
    walls: wgpu::Buffer,
    wall_count: usize,
    /// Van Driest's damping, when it has wall cells.
    damping: Option<Damping>,
    /// The uniforms of the two phases (the body acceleration can change between steps).
    params: [wgpu::Buffer; 2],
    open_count: usize,
    open_groups: (u32, u32),
    /// Whether a pressure face needs the preparatory outlet pass.
    prepare_outlets: bool,
    group_count: usize,
    force_groups: (u32, u32),
    stats: wgpu::Buffer,
    moments: wgpu::Buffer,
    forces: wgpu::Buffer,
    partials: wgpu::Buffer,
    partial_count: usize,
    /// Where the fans and porous zones are in the flags buffer ([`crate::zones`]).
    zones: crate::zones::Layout,
    /// The workgroups of `fan_partials` (one per 256 entries of the fans' cell lists).
    fan_groups: (u32, u32),
    /// The scalar populations remain A–B with heat; and T₀ they are shifted by.
    thermal: Option<[wgpu::Buffer; 2]>,
    initial_temperature: f64,
    fluid: Vec<bool>,
    fluid_cells: usize,
    groups: Vec<Group>,
    empty_obstacles: Vec<String>,
    /// What the half-way wall made of the obstacles' walls (the flow's report), for the solve's warnings.
    wall_layer: Option<fieldcraft_lbm::wall_layer::WallLayerReport>,
    /// With interpolated walls, how the obstacles' links bounce back (the flow's).
    wall_links: Option<fieldcraft_lbm::WallLinks>,
    steps_done: u64,
    // Last uploaded uniform acceleration, used only by the momentum diagnostic.
    diagnostic_acceleration: [f64; 3],
    /// With Smagorinsky's model: the fluid's own viscosity (lattice units), and the largest eddy
    /// viscosity of the last chunk over it.
    nu0: f64,
    eddy_ratio: f64,
    /// The momentum (xyz) and heat (w) each group received in the last step, once read back
    /// (zero before the first step).
    exchange_cache: Option<Vec<[f64; 4]>>,
    /// Each group's share of the momentum from the populations' weights, which the kernel leaves
    /// out (`link_forces` sums the deviations): added to what it reads back.
    weight_forces: Vec<[f64; 3]>,
    /// The conducting solids, when there are any.
    solid: Option<Solids>,
    /// The Spalart–Allmaras model, when the flow has it, and its per-cell vorticity and wall
    /// distance buffer.
    sa: Option<SaDevice>,
    sa_cells: Option<wgpu::Buffer>,
    /// With the model on Han's wall, its wall function's rows (`sa_wall_function`) and their
    /// dispatch.
    sa_wall_rows: u32,
    sa_wall_groups: Option<(u32, u32)>,
}

/// The Spalart–Allmaras model on the device (see `lbm.wgsl`): its ν̃ buffer and the chunk's
/// sub-step count.
struct SaDevice {
    nu_tilde: wgpu::Buffer,
    /// The sub-steps each lattice step of the next chunk takes; none until the first chunk
    /// probes the initial field.
    substeps: Option<u32>,
    /// Whether the count is fixed ([`GpuFlow::fix_sa_substeps`]) rather than planned per chunk.
    fixed: bool,
}

/// A boundary node's side row (lbm.wgsl `node_side`): its Han links, its links to the fluid (one bit
/// per direction) and its flags (flat, sampled, the stress axes).
fn node_roles(w: &fieldcraft_lbm::bounce::WallNode) -> (u32, u32, u32) {
    use fieldcraft_lbm::bounce::{Role, carries};
    let (mut han, mut to_fluid) = (0u32, 0u32);
    for (q, role) in w.roles.iter().enumerate().skip(1) {
        let bit = 1u32 << q;
        if *role == Role::Han {
            han |= bit;
        }
        if matches!(role, Role::Han | Role::Edge) {
            to_fluid |= bit;
        }
    }
    let (axis, inward) = (usize::from(w.axis), i32::from(w.inward));
    let mut flags = u32::from(w.flat) | (u32::from(w.sampled()) << 1);
    for b in 0..3 {
        if carries(&w.roles, axis, inward, b) {
            flags |= 1 << (2 + b);
        }
    }
    (han, to_fluid, flags)
}

/// The Spalart–Allmaras sub-steps per lattice step, at most (`fieldcraft_lbm::sa::MAX_SUBSTEPS`).
const SA_MAX_SUBSTEPS: u32 = fieldcraft_lbm::sa::MAX_SUBSTEPS;
/// The margin of the next chunk's count over the largest coefficient sum of the last.
const SA_MARGIN: f64 = 1.25;
/// Past this many times the monotone bound (dt · coefficient sum) in a chunk, the explicit
/// transport is no longer stable and the run stops. Between 1 and it the chunk has run past the
/// bound (ν̃ may dip below 0, which SA-neg handles; the CPU redoes the step instead) and the next
/// chunk takes more sub-steps.
const SA_STOP: f64 = 2.0;
/// `params.sa`'s byte offset in the uniforms, and the sub-step dt's within it.
const SA_PARAMS: u64 = 192;
const SA_DT: u64 = SA_PARAMS + 8;

/// The Spalart–Allmaras model's box-face kinds for `lbm.wgsl` (2 bits per face: 1 wall,
/// 2 inflow, 3 open), checked against every fluid cell's faces on the CPU: the device's rule is a
/// fluid neighbour, a boundary node of Han's on-node wall, a solid cell (a halfway wall) or,
/// beyond a face that is not periodic, that face's one kind. A face it cannot reproduce (an inlet
/// patch) is refused.
fn sa_faces(flow: &Flow, field: &fieldcraft_lbm::sa::SaField) -> Result<u32> {
    use fieldcraft_lbm::sa::Side;
    let refuse = |what: String| {
        Error::Invalid(format!("not supported yet: the Spalart–Allmaras model on the GPU with {what}; solve on the CPU (`backend: \"cpu\"`)"))
    };
    let grid = flow.grid();
    let periodic = flow.periodic();
    let edge = |p: [usize; 3], k: usize| {
        let axis = k / 2;
        let at = if k % 2 == 1 { p[axis] + 1 == grid.cells[axis] } else { p[axis] == 0 };
        at && !periodic[axis]
    };
    let mut nodes = vec![false; grid.len()];
    for w in flow.wall_nodes().map_or(&[][..], |n| n.nodes()) {
        if let Some(x) = nodes.get_mut(w.cell as usize) {
            *x = true;
        }
    }
    let mut kinds = [None::<u32>; 6];
    for (cell, sides) in field.sides().iter().enumerate() {
        if !flow.is_fluid(cell) {
            continue;
        }
        let p = grid.coords(cell);
        for (k, side) in sides.iter().enumerate() {
            if !edge(p, k) {
                continue;
            }
            let kind = match side {
                Side::Wall => 1,
                Side::Inflow => 2,
                Side::Open => 3,
                other => return Err(refuse(format!("a box face whose cells meet {other:?}"))),
            };
            match kinds.get_mut(k) {
                Some(slot @ None) => *slot = Some(kind),
                Some(Some(prev)) if *prev != kind => return Err(refuse("two kinds of face on one box face (an inlet patch)".into())),
                _ => {}
            }
        }
    }
    for (cell, sides) in field.sides().iter().enumerate() {
        if !flow.is_fluid(cell) {
            continue;
        }
        let p = grid.coords(cell);
        for (k, side) in sides.iter().enumerate() {
            let (axis, up) = (k / 2, k % 2 == 1);
            let device = if edge(p, k) {
                match kinds.get(k).copied().flatten() {
                    Some(1) => Side::Wall,
                    Some(2) => Side::Inflow,
                    _ => Side::Open,
                }
            } else {
                let n = grid.cells[axis];
                let mut q = p;
                q[axis] = if up { (p[axis] + 1) % n } else { (p[axis] + n - 1) % n };
                let nb = grid.index(q);
                if flow.is_fluid(nb) {
                    Side::Cell(u32::try_from(nb).map_err(|_| Error::Limit("too many cells for the GPU solver".into()))?)
                } else if nodes.get(nb).copied().unwrap_or(false) {
                    Side::Node
                } else {
                    Side::Wall
                }
            };
            if device != *side {
                return Err(refuse(format!(
                    "a face the device's neighbour rule does not reproduce (cell {cell}: {side:?}, the device's {device:?})"
                )));
            }
        }
    }
    Ok(kinds.iter().enumerate().fold(0u32, |bits, (k, kind)| bits | kind.unwrap_or(3) << (2 * k)))
}

/// Van Driest's damping on the device (see `lbm.wgsl`, `damping`).
struct Damping {
    /// The friction velocities first, then the wall cells' rows and the cells' distances.
    buffer: wgpu::Buffer,
    groups: (u32, u32),
    /// The wall cells' distances to their walls, and the fluid's viscosity (lattice units), for
    /// their y⁺.
    distances: Vec<f64>,
    nu: f64,
}

/// Van Driest's damping data in the layout of `lbm.wgsl`'s `damping`, and its wall cells'
/// distances: none without damping or without wall cells.
fn damping_data(flow: &Flow, n: usize) -> Result<Option<(Vec<f32>, Vec<f64>)>> {
    let Some(walls) = flow.walls().filter(|w| w.damped()) else { return Ok(None) };
    let distances: Vec<f64> = walls.wall_cells().map(|w| w.2).collect();
    let d = distances.len();
    if d == 0 {
        return Ok(None);
    }
    let too_many = || Error::Limit("too many cells for the GPU's wall distances".into());
    let len = n.checked_mul(2).and_then(|v| v.checked_add(9 * d)).filter(|v| u32::try_from(*v).is_ok()).ok_or_else(too_many)?;
    let mut data = reserved(len, "wall distances")?;
    data.resize(d, 0.0);
    for (cell, link, distance, velocity, _) in walls.wall_cells() {
        let c = fieldcraft_lbm::lattice::c(link);
        let norm = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
        data.extend(c.map(|v| (v / norm) as f32));
        data.push(distance as f32);
        data.extend(velocity.map(|v| v as f32));
        data.push(f32::from_bits(u32::try_from(cell).map_err(|_| too_many())?));
    }
    for cell in 0..n {
        match (walls.distance(cell), walls.owner(cell).and_then(|o| u32::try_from(o).ok())) {
            (Some(y), Some(owner)) => data.extend([y as f32, f32::from_bits(owner)]),
            _ => data.extend([0.0, f32::from_bits(u32::MAX)]),
        }
    }
    Ok(Some((data, distances)))
}

/// The conducting solids on the device (see `lbm.wgsl`, "Conducting solids").
struct Solids {
    /// Their state (the temperatures first, as θ = T − T₀).
    state: wgpu::Buffer,
    /// Their cells (grid indices) and each solid's cells.
    cells: Vec<u32>,
    by_solid: Vec<Vec<u32>>,
    substeps: u32,
    cell_groups: (u32, u32),
    link_groups: (u32, u32),
}

/// A conducting cell's face kinds, as `lbm.wgsl` names them.
const FACE_NONE: u32 = 0;
const FACE_SOLID: u32 = 1;
const FACE_HELD: u32 = 2;
const FACE_FLUX: u32 = 3;

/// The conducting solids on the device: their data in the layouts `lbm.wgsl` reads ("Conducting
/// solids"), its params, and what the host keeps of them.
struct SolidUpload {
    u: Vec<u32>,
    f: Vec<f32>,
    state: Vec<f32>,
    params: [u32; 4],
    cells: Vec<u32>,
    by_solid: Vec<Vec<u32>>,
    substeps: u32,
}

/// The conducting solids' data, temperatures as θ relative to `t0`.
fn solid_data(solid: &fieldcraft_lbm::conduction::SolidHeat, t0: f64) -> Result<SolidUpload> {
    use fieldcraft_lbm::conduction::Face;
    let u = |v: usize| u32::try_from(v).map_err(|_| Error::Limit("too many conducting cells for the GPU".into()));
    let m = solid.cells.len();
    let k = solid.interface.len();
    // The interface links ordered by conducting cell (stable: the CPU's order within a cell),
    // and each cell's first one.
    let mut links: Vec<&fieldcraft_lbm::conduction::Interface> = solid.interface.iter().collect();
    links.sort_by_key(|l| l.solid);
    let mut count = vec![0usize; m + 1];
    for l in &links {
        if let Some(c) = count.get_mut(l.solid as usize + 1) {
            *c += 1;
        }
    }
    let mut su = reserved(2 * m + 1 + 4 * k + 12 * m, "conducting cells")?;
    su.extend_from_slice(&solid.cells);
    let mut start = 0usize;
    for c in &count[..m + 1] {
        start += c;
        su.push(u(start)?);
    }
    for l in &links {
        su.extend_from_slice(&[l.cell, u32::from(l.q), l.solid, l.link]);
    }
    let mut sf = reserved(14 * m + k, "conducting cells")?;
    sf.extend(solid.cell.iter().map(|c| c.capacity as f32));
    sf.extend(solid.cell.iter().map(|c| c.source as f32));
    for c in &solid.cell {
        for face in &c.faces {
            // A wall's `other` is its box face, where the heat it takes is counted.
            let (kind, other, g, value) = match *face {
                Face::None => (FACE_NONE, 0, 0.0, 0.0),
                Face::Solid { other, g } => (FACE_SOLID, other, g, 0.0),
                Face::Held { t, g, wall } => (FACE_HELD, u32::from(wall), g, t - t0),
                Face::Flux { q, wall } => (FACE_FLUX, u32::from(wall), 0.0, q),
            };
            su.extend_from_slice(&[kind, other]);
            sf.extend_from_slice(&[g as f32, value as f32]);
        }
    }
    sf.extend(links.iter().map(|l| l.weight as f32));
    let mut state = filled(11 * m + 3 * k, 0.0f32, "conducting state")?;
    for (slot, temp) in state.iter_mut().zip(&solid.temperature) {
        *slot = (temp - t0) as f32;
    }
    let substeps = solid.substeps.max(1);
    Ok(SolidUpload {
        u: su,
        f: sf,
        state,
        params: [u(m)?, u(k)?, 0, (1.0 / substeps as f32).to_bits()],
        cells: solid.cells.clone(),
        by_solid: solid.by_solid.clone(),
        substeps,
    })
}

/// The entry points of `lbm.wgsl` a flow dispatches.
const KERNELS: [&str; 39] = [
    "clear_eddy",
    "gradient_velocity",
    "fan_partials",
    "fan_update",
    "zone_samples",
    "collide",
    "wall_sample",
    "damping_refresh",
    "bounce_sample",
    "node_bounce",
    "stream",
    "stream_links",
    "stream_links_exchange",
    "scatter_links",
    "outlet_convect",
    "open_boundaries",
    "link_forces",
    "wall_forces",
    "group_forces",
    "cell_moments",
    "wall_moments",
    "mass",
    "mean_velocity",
    "collide_thermal",
    "stream_thermal",
    "stream_links_thermal",
    "open_thermal",
    "thermal_check",
    "solid_ab",
    "solid_sub0",
    "solid_sub1",
    "solid_copy",
    "solid_scatter",
    "solid_contrib",
    "solid_wall_heat",
    "sa_prepare",
    "sa_rate",
    "sa_update",
    "sa_wall_function",
];

/// The bindings each kernel uses, read from the shader (naga's analysis of which globals an entry
/// point reaches): each kernel gets a bind group of only those. One layout of all 26 buffers
/// needs 25 storage buffers per compute stage, more than GL's 16 (#640); the busiest kernel
/// (`collide`, with the gradient velocities of the hybrid collision and Vreman's model) needs 8,
/// WebGPU's default.
struct Reflection {
    kernels: Vec<(&'static str, Vec<(u32, wgpu::BufferBindingType)>)>,
    /// The most storage buffers one kernel binds.
    storage_buffers: u32,
}

fn reflection() -> Result<&'static Reflection> {
    static REFLECTION: OnceLock<std::result::Result<Reflection, String>> = OnceLock::new();
    REFLECTION.get_or_init(reflect).as_ref().map_err(|e| gpu_error(format!("the fluid shader does not parse: {e}")))
}

fn reflect() -> std::result::Result<Reflection, String> {
    use wgpu::naga;
    let module = naga::front::wgsl::parse_str(SHADER).map_err(|e| e.emit_to_string(SHADER))?;
    let info = naga::valid::Validator::new(naga::valid::ValidationFlags::all(), naga::valid::Capabilities::all())
        .validate(&module)
        .map_err(|e| e.emit_to_string(SHADER))?;
    let mut kernels = Vec::with_capacity(KERNELS.len());
    for name in KERNELS {
        let index = module.entry_points.iter().position(|e| e.name == name).ok_or_else(|| format!("no entry point {name}"))?;
        let uses = info.get_entry_point(index);
        let mut bindings = Vec::new();
        for (handle, var) in module.global_variables.iter() {
            let (Some(binding), false) = (&var.binding, uses[handle].is_empty()) else { continue };
            let ty = match var.space {
                naga::AddressSpace::Uniform => wgpu::BufferBindingType::Uniform,
                naga::AddressSpace::Storage { access } => {
                    wgpu::BufferBindingType::Storage { read_only: !access.contains(naga::StorageAccess::STORE) }
                }
                _ => return Err(format!("{name}: binding {} is not a buffer", binding.binding)),
            };
            if binding.group != 0 {
                return Err(format!("{name}: binding {} is not in group 0", binding.binding));
            }
            bindings.push((binding.binding, ty));
        }
        bindings.sort_by_key(|(b, _)| *b);
        kernels.push((name, bindings));
    }
    let storage = |b: &[(u32, wgpu::BufferBindingType)]| b.iter().filter(|(_, t)| matches!(t, wgpu::BufferBindingType::Storage { .. })).count();
    let storage_buffers = kernels.iter().map(|(_, b)| storage(b)).max().unwrap_or(0) as u32;
    Ok(Reflection { kernels, storage_buffers })
}

/// Refuse an adapter that allows fewer storage buffers per compute stage than a kernel binds,
/// before anything is built, with what to do instead.
fn check_storage_buffers(backend: wgpu::Backend, have: u32, need: u32) -> Result<()> {
    if have >= need {
        return Ok(());
    }
    let backend = match backend {
        wgpu::Backend::Gl => "OpenGL",
        wgpu::Backend::Vulkan => "Vulkan",
        wgpu::Backend::Metal => "Metal",
        wgpu::Backend::Dx12 => "DX12",
        wgpu::Backend::BrowserWebGpu => "WebGPU",
        wgpu::Backend::Noop => "no-op",
    };
    Err(Error::Backend(format!(
        "this GPU's {backend} backend allows {have} storage buffers per compute stage and the fluid solver needs {need}: solve with `backend: \"cpu\"`, or use a Vulkan, Metal or DX12 adapter (`WGPU_BACKEND`)"
    )))
}

/// The override constants a flow's kernels are compiled with.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Variant {
    aa: bool,
    regularized: bool,
    mrt: bool,
    smagorinsky: bool,
    wall_model: bool,
    van_driest: bool,
    wall_bounce: bool,
    /// Han et al.'s on-node wall-function bounce (`node_bounce`).
    on_node: bool,
    /// The hybrid collision, Vreman's model, and the hybrid collision's dynamic σ.
    hrr: bool,
    vreman: bool,
    dynamic: bool,
    /// Fans or porous zones.
    zones: bool,
    /// Interpolated wall links.
    interpolated: bool,
    /// The Spalart–Allmaras model.
    sa: bool,
    /// Its wall function is the model's own constant-stress law (`sa_law`), not Spalding's.
    sa_law: bool,
}

/// A compiled kernel: its pipeline, its bind group's layout and the bindings in it.
struct Kernel {
    pipeline: wgpu::ComputePipeline,
    layout: wgpu::BindGroupLayout,
    bindings: Vec<u32>,
}

type Kernels = HashMap<&'static str, Kernel>;

/// The kernels of `variant` on `device`, each with the layout of the bindings it uses.
fn compile(device: &wgpu::Device, variant: Variant) -> Result<Kernels> {
    let module =
        device.create_shader_module(wgpu::ShaderModuleDescriptor { label: Some("fieldcraft lbm"), source: wgpu::ShaderSource::Wgsl(SHADER.into()) });
    let mut kernels = HashMap::with_capacity(KERNELS.len());
    for (name, bindings) in &reflection()?.kernels {
        let entries: Vec<wgpu::BindGroupLayoutEntry> = bindings
            .iter()
            .map(|&(binding, ty)| wgpu::BindGroupLayoutEntry {
                binding,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer { ty, has_dynamic_offset: false, min_binding_size: None },
                count: None,
            })
            .collect();
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor { label: Some(name), entries: &entries });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some(name),
            bind_group_layouts: &[Some(&layout)],
            immediate_size: 0,
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some(name),
            layout: Some(&pipeline_layout),
            module: &module,
            entry_point: Some(name),
            compilation_options: wgpu::PipelineCompilationOptions {
                constants: &[
                    ("AA_MODE", if variant.aa { 1.0 } else { 0.0 }),
                    ("REGULARIZED", if variant.regularized { 1.0 } else { 0.0 }),
                    ("MRT", if variant.mrt { 1.0 } else { 0.0 }),
                    ("SMAGORINSKY", if variant.smagorinsky { 1.0 } else { 0.0 }),
                    ("WALL_MODEL", if variant.wall_model { 1.0 } else { 0.0 }),
                    ("VAN_DRIEST", if variant.van_driest { 1.0 } else { 0.0 }),
                    ("WALL_BOUNCE", if variant.wall_bounce { 1.0 } else { 0.0 }),
                    ("ON_NODE", if variant.on_node { 1.0 } else { 0.0 }),
                    ("HRR", if variant.hrr { 1.0 } else { 0.0 }),
                    ("VREMAN", if variant.vreman { 1.0 } else { 0.0 }),
                    ("DYNAMIC", if variant.dynamic { 1.0 } else { 0.0 }),
                    ("ZONES", if variant.zones { 1.0 } else { 0.0 }),
                    ("INTERPOLATED", if variant.interpolated { 1.0 } else { 0.0 }),
                    ("SA", if variant.sa { 1.0 } else { 0.0 }),
                    ("SA_LAW", if variant.sa_law { 1.0 } else { 0.0 }),
                    ("MAX_SPEED", MAX_LATTICE_SPEED),
                ],
                ..Default::default()
            },
            cache: None,
        });
        kernels.insert(*name, Kernel { pipeline, layout, bindings: bindings.iter().map(|(b, _)| *b).collect() });
    }
    Ok(kernels)
}

/// A device opened for fluid solves, kept for the process with the kernels compiled on it, so a
/// solve after the first compiles nothing (#588: DX12 keeps its compiled shaders per device, and
/// recompiling them with FXC took 9 s per solve). Pipelines cannot outlive their device, so the
/// device is what is kept; wgpu's pipeline cache, which could keep them across runs, is unsafe
/// to create (and a stub on DX12).
struct Shared {
    /// The adapter and the features the device was opened with.
    key: String,
    device: wgpu::Device,
    queue: wgpu::Queue,
    /// Set by the device-lost callback: the next solve opens a new device.
    lost: Arc<AtomicBool>,
    kernels: Mutex<Vec<(Variant, Arc<Kernels>)>>,
}

/// The devices kept so far (one per adapter and feature set).
static DEVICES: Mutex<Vec<Arc<Shared>>> = Mutex::new(Vec::new());
/// Successful LBM compute-device creations, including devices subsequently lost/retired.
static DEVICES_OPENED: AtomicU64 = AtomicU64::new(0);

/// The kept device for `adapter` with `features`, or a new one, opened with the adapter's buffer
/// sizes and `storage_buffers` per compute stage.
fn shared_device(adapter: &wgpu::Adapter, features: wgpu::Features, storage_buffers: u32) -> Result<(Arc<Shared>, bool)> {
    let info = adapter.get_info();
    let key = format!("{:?} {} {:04x}:{:04x} {:?} {} {features:?}", info.backend, info.name, info.vendor, info.device, info.device_type, info.driver);
    let mut devices = DEVICES.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    devices.retain(|d| !d.lost.load(Ordering::Acquire));
    if let Some(d) = devices.iter().find(|d| d.key == key) {
        return Ok((d.clone(), false));
    }
    let have = adapter.limits();
    let mut limits = wgpu::Limits::default().using_resolution(have.clone());
    limits.max_buffer_size = have.max_buffer_size;
    limits.max_storage_buffer_binding_size = have.max_storage_buffer_binding_size;
    // Exactly what the busiest kernel binds, so every device (the CI's lavapipe too) proves the
    // kernels fit an adapter that allows no more.
    limits.max_storage_buffers_per_shader_stage = storage_buffers;
    update_diagnostic_phase(Some("device_requested"));
    let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
        label: Some("fieldcraft lbm"),
        required_features: features,
        required_limits: limits,
        ..Default::default()
    }))
    .map_err(gpu_error)?;
    DEVICES_OPENED.fetch_add(1, Ordering::Relaxed);
    let lost = Arc::new(AtomicBool::new(false));
    let flag = lost.clone();
    device.set_device_lost_callback(move |reason, message| {
        log::warn!("the fluid solver's GPU device was lost ({reason:?}): {message}");
        flag.store(true, Ordering::Release);
    });
    let shared = Arc::new(Shared { key, device, queue, lost, kernels: Mutex::new(Vec::new()) });
    devices.push(shared.clone());
    Ok((shared, true))
}

impl Shared {
    /// The kernels of `variant`, compiled on the first solve that needs them; and how long that
    /// took (zero when they were compiled before).
    fn kernels(&self, variant: Variant) -> Result<(Arc<Kernels>, f64)> {
        let mut built = self.kernels.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some((_, k)) = built.iter().find(|(v, _)| *v == variant) {
            return Ok((k.clone(), 0.0));
        }
        fieldcraft_model::progress::stage("Compiling the GPU kernels", None, "");
        let start = std::time::Instant::now();
        let scope = self.device.push_error_scope(wgpu::ErrorFilter::Validation);
        let kernels = compile(&self.device, variant)?;
        // Kept only when they compiled.
        if let Some(e) = pollster::block_on(scope.pop()) {
            return Err(gpu_error(e));
        }
        let kernels = Arc::new(kernels);
        built.push((variant, kernels.clone()));
        Ok((kernels, start.elapsed().as_secs_f64()))
    }
}

impl GpuFlow {
    /// `flow` on the device `FIELDCRAFT_GPU` picks.
    pub fn from_env(flow: &Flow) -> Result<GpuFlow> {
        let policy = Policy::from_env().map_err(gpu_error)?;
        GpuFlow::new(flow, policy)
    }

    /// `flow` (its grid, rules and initial state: everything the step needs) on the device
    /// `policy` picks.
    pub fn new(flow: &Flow, policy: Policy) -> Result<GpuFlow> {
        GpuFlow::with_options(flow, policy, Options::default())
    }

    /// [`GpuFlow::new`] with `options`.
    pub fn with_options(flow: &Flow, policy: Policy, options: Options) -> Result<GpuFlow> {
        if flow.has_baffles() {
            return Err(Error::Invalid("baffles are not supported on the GPU yet: solve on the CPU".into()));
        }
        if let Some(r) = flow.rheology() {
            return Err(Error::Invalid(format!(
                "the GPU backend does not step a shear-dependent viscosity yet (a {} fluid): solve on the CPU (`backend: \"cpu\"`)",
                r.fluid.name()
            )));
        }
        if flow.phase_change() {
            return Err(Error::Invalid(
                "melting and freezing (a phase change in the fluid's material) are not on the GPU backend yet: solve on the CPU (`backend: \"cpu\"`)"
                    .into(),
            ));
        }
        if flow.free_surface() {
            return Err(Error::Invalid("a free surface is not on the GPU backend yet: solve on the CPU (`backend: \"cpu\"`)".into()));
        }
        if flow.radiates() {
            return Err(Error::Invalid(
                "radiation from conducting solids is not on the GPU backend yet: solve on the CPU (`backend: \"cpu\"`)".into(),
            ));
        }
        if flow.wall_nodes().is_some_and(|n| !matches!(n.law(), fieldcraft_lbm::sa::WallLaw::Spalding | fieldcraft_lbm::sa::WallLaw::SpalartAllmaras))
        {
            return Err(Error::Invalid("not supported yet on the GPU: this wall law; solve on the CPU (`backend: \"cpu\"`)".into()));
        }
        // The Spalart–Allmaras model's parts not on the device yet, and its faces.
        let sa_faces = match flow.spalart_allmaras() {
            None => None,
            Some(field) => {
                let refuse = |what: &str| {
                    Err(Error::Invalid(format!(
                        "not supported yet: the Spalart–Allmaras model on the GPU with {what}; solve on the CPU (`backend: \"cpu\"`)"
                    )))
                };
                if flow.wall_model() {
                    return refuse("Malaspinas and Sagaut's wall model with Musker's profile (`musker`)");
                }
                if flow.wall_bounce() {
                    return refuse("the halfway wall-function bounce (`spalding_bounce_halfway`)");
                }
                if flow.thermal().is_some() {
                    return refuse("heat");
                }
                if flow.hybrid_sigma() == Some(fieldcraft_model::HrrSigma::Dynamic) {
                    return refuse("the hybrid collision's dynamic σ");
                }
                if flow.sa_wall_distance().is_none_or(|d| d.len() != flow.grid().len()) {
                    return refuse("no wall distance");
                }
                Some(sa_faces(flow, field)?)
            }
        };
        fieldcraft_model::progress::backend_diagnostic(fieldcraft_model::progress::BackendDiagnostic {
            phase: "selecting_adapter".into(),
            devices_opened: Some(DEVICES_OPENED.load(Ordering::Relaxed)),
            instances_opened: Some(select::instances_opened()),
            counter_scope: "devices: successful LBM compute creations, including retired devices; viewport, meshing and external/driver devices are not counted. instances: gpu::select creations; external viewport and driver-internal instances are not counted.".into(),
            ..Default::default()
        });
        // Validate zone values before opening a device or doing device work.
        let n = flow.grid().len();
        let mut fluid_u32 = filled(n, 0u32, "cell flags")?;
        for (cell, flag) in fluid_u32.iter_mut().enumerate() {
            *flag = u32::from(flow.is_fluid(cell));
        }
        let zones = crate::zones::append(flow, &mut fluid_u32)?;
        if sa_faces.is_some() && zones.zones > 0 {
            return Err(Error::Invalid(
                "not supported yet: the Spalart–Allmaras model on the GPU with fans or porous zones; solve on the CPU (`backend: \"cpu\"`)".into(),
            ));
        }
        let Some((adapter, tier)) = select::adapter(policy).map_err(gpu_error)? else {
            return Err(Error::Backend(format!("{}=cpu: no GPU device; solve with `backend: \"cpu\"`", select::ENV)));
        };
        let info = adapter.get_info();
        if let Some(progress) = fieldcraft_model::progress::current()
            && let Some(mut context) = progress.backend_diagnostic()
        {
            context.phase = "adapter_selected".into();
            context.backend = Some(format!("{:?}", info.backend));
            context.adapter = Some(info.name.clone());
            context.driver = Some(info.driver.clone());
            context.driver_info = Some(info.driver_info.clone());
            context.adapter_limits = Some(diagnostic_limits(&adapter.limits()));
            context.instances_opened = Some(select::instances_opened());
            progress.set_backend_diagnostic(Some(context));
        }
        let profile = options.profile;
        let features = if profile { wgpu::Features::TIMESTAMP_QUERY } else { wgpu::Features::empty() };
        if !adapter.features().contains(features) {
            return Err(gpu_error("profiling needs timestamp queries, which this adapter does not offer"));
        }
        let storage_buffers = reflection()?.storage_buffers;
        check_storage_buffers(info.backend, adapter.limits().max_storage_buffers_per_shader_stage, storage_buffers)?;
        let (shared, created) = shared_device(&adapter, features, storage_buffers)?;
        let (device, queue) = (&shared.device, &shared.queue);
        if let Some(progress) = fieldcraft_model::progress::current()
            && let Some(mut context) = progress.backend_diagnostic()
        {
            context.phase = if created { "device_created" } else { "device_reused" }.into();
            context.device_limits = Some(diagnostic_limits(&device.limits()));
            progress.set_backend_diagnostic(Some(context));
        }
        update_diagnostic_phase(None);
        let adapter_name = format!("{}{}", info.name, if tier == select::Tier::Software { " (software)" } else { "" });

        let aa = match options.streaming {
            Streaming::Auto | Streaming::InPlace => true,
            Streaming::TwoBuffers => false,
        };
        let elem = 4u64;
        let pop_bytes = elem * (Q * n) as u64;
        let cap = buffer_cap(&device.limits());
        if pop_bytes > cap {
            use fieldcraft_lbm::memory::megabytes;
            return Err(Error::Limit(format!(
                "{n} cells need {} of populations per buffer; this GPU allows {} per buffer: coarsen the spacing",
                megabytes(pop_bytes),
                megabytes(cap)
            )));
        }

        let flag_bytes = (fluid_u32.len() as u64).checked_mul(4).ok_or_else(|| Error::Limit("too much GPU zone data".into()))?;
        if flag_bytes > cap {
            return Err(Error::Limit(format!("GPU cell flags and fan/porous data need {flag_bytes} bytes; this device allows {cap} per buffer")));
        }

        // The state and the flags. The populations go up one direction at a time (below), so
        // the host never holds a whole single-precision copy of them beside the flow's.
        let thermal = flow.thermal();
        let t0 = thermal.map_or(0.0, |t| t.initial());
        let mut fluid = filled(n, false, "cell flags")?;
        for cell in 0..n {
            if let Some(b) = fluid.get_mut(cell) {
                *b = flow.is_fluid(cell);
            }
        }
        // The wall model's wall cells with a second cell: four rows each in `walls`, and in each
        // one's flag word its row (k + 1) past the fluid bit.
        let mut wall_rows: Vec<f32> = Vec::new();
        let nu0 = (1.0 / flow.omegas().0 - 0.5) / 3.0;
        // With the wall-function bounce instead, the wall cells on lattice planes (their distance
        // half a cell), with the cell itself where the second cell goes; their flags let link
        // streaming shift the populations reflected into them (`bounce_shift`).
        let bounce = flow.wall_bounce();
        if let Some(walls) = flow.walls().filter(|_| flow.wall_model() || bounce) {
            let links = flow.links();
            for (cell, link, distance, velocity, second) in walls.wall_cells() {
                let second = if bounce { (distance == 0.5).then_some(cell) } else { second };
                let Some(second) = second else { continue };
                let row = u32::try_from(wall_rows.len() / WALL_ROW + 1)
                    .ok()
                    .filter(|r| *r < 1 << 31)
                    .ok_or_else(|| Error::Limit("too many wall cells for the GPU's flags".into()))?;
                let second = u32::try_from(second).map_err(|_| Error::Limit("too many cells for the GPU solver".into()))?;
                if let Some(u) = fluid_u32.get_mut(cell) {
                    *u = 1 | row << 1;
                }
                let c = fieldcraft_lbm::lattice::c(link);
                let len = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
                wall_rows.extend(c.map(|v| (v / len) as f32));
                wall_rows.push(distance as f32);
                wall_rows.extend(velocity.map(|v| v as f32));
                wall_rows.push(f32::from_bits(second));
                wall_rows.extend(velocity.map(|v| v as f32));
                wall_rows.push((1.0 / (3.0 * nu0 + 0.5)) as f32);
                // Its nearest wall link (the links are in cell order), which takes the momentum
                // its reconstruction removes, as on the CPU.
                let first = links.partition_point(|l| (l.cell as usize) < cell);
                let nearest = links
                    .get(first..)
                    .unwrap_or_default()
                    .iter()
                    .take_while(|l| l.cell as usize == cell)
                    .position(|l| usize::from(l.q) == link)
                    .and_then(|i| u32::try_from(first + i).ok())
                    .ok_or_else(|| Error::Invalid(format!("wall cell {cell} has no link in its wall's direction")))?;
                wall_rows.extend([0.0, 0.0, 0.0, f32::from_bits(nearest)]);
            }
        }
        // Han et al.'s on-node bounce instead: its boundary nodes in `walls` (lbm.wgsl,
        // `node_bounce`), their flags NODE, and for each where what it took from the fluid goes:
        // on a plane box face, the link of its face's group from outside the box along the normal
        // (an open placeholder with no exchange of its own); anywhere else, a slot of its own
        // after the links. After the nodes' four rows, a side row pair per node (`node_side`):
        // its roles from the CPU (`bounce::roles`) as bit masks, its cell and flags.
        let mut node_links: Vec<(u32, u32)> = Vec::new();
        let mut node_sides: Vec<f32> = Vec::new();
        let mut node_slots = 0u32;
        if let Some(nodes) = flow.wall_nodes() {
            let links = flow.links();
            for w in nodes.nodes() {
                let (cell, axis) = (w.cell as usize, usize::from(w.axis).min(2));
                let mut normal = [0i32; 3];
                normal[axis] = i32::from(w.inward.signum());
                let q = C.iter().position(|c| *c == normal).ok_or_else(|| Error::Invalid(format!("boundary node {cell} has no normal")))?;
                if let Some(u) = fluid_u32.get_mut(cell) {
                    *u = NODE;
                }
                let first = links.partition_point(|l| (l.cell as usize) < cell);
                let open = links
                    .get(first..)
                    .unwrap_or_default()
                    .iter()
                    .take_while(|l| l.cell as usize == cell)
                    .position(|l| usize::from(l.q) == q && matches!(l.rule, Rule::Open))
                    .and_then(|i| u32::try_from(first + i).ok());
                let link = match open.filter(|_| w.flat && (w.group as usize) < BoxFace::ALL.len()) {
                    Some(l) => l,
                    None => {
                        let slot = u32::try_from(links.len())
                            .ok()
                            .and_then(|l| l.checked_add(node_slots))
                            .ok_or_else(|| Error::Limit("too many links and boundary nodes for the GPU solver".into()))?;
                        node_slots += 1;
                        slot
                    }
                };
                let (han, to_fluid, flags) = node_roles(w);
                node_sides.extend([f32::from_bits(han), f32::from_bits(to_fluid), f32::from_bits(w.cell), f32::from_bits(flags)]);
                node_sides.extend([0.0; 4]);
                let read = if nodes.second_cell() { w.sample_cell } else { w.interior };
                wall_rows.extend(normal.map(|v| v as f32));
                wall_rows.push(nodes.law_distance(w) as f32);
                wall_rows.extend([0.0, 0.0, 0.0, f32::from_bits(w.interior)]);
                wall_rows.extend([0.0, 0.0, 0.0, f32::from_bits(read)]);
                wall_rows.extend([0.0, 0.0, 0.0, f32::from_bits(link)]);
                node_links.push((w.group, link));
            }
        }
        let on_node = !node_links.is_empty();
        let wall_count = wall_rows.len() / WALL_ROW;
        wall_rows.extend_from_slice(&node_sides);
        drop(node_sides);

        // The links, by cell, and the ones with a force by group.
        let links = flow.links();
        let mut link_u = reserved(4 * links.len(), "links")?;
        // With interpolated walls, also each link's fraction and each group's wall velocity
        // (lbm.wgsl, `link_f`).
        let interpolations = flow.interpolations();
        let interpolated = !interpolations.is_empty();
        let mut link_q = reserved(if interpolated { links.len() } else { 0 }, "link fractions")?;
        let mut link_f = reserved(if interpolated { 2 * links.len() + 4 * flow.groups().len() } else { links.len() }, "links")?;
        let mut by_group: Vec<(u32, u32)> = reserved(links.len(), "links by group")?;
        let group_count = flow.groups().len();
        let mut weight_forces = filled(group_count, [0.0f64; 3], "force groups")?;
        let c = |q: u32| C.get(q as usize).copied().unwrap_or([0; 3]).map(f64::from);
        let w = |q: u32| W.get(q as usize).copied().unwrap_or(0.0);
        for (l, link) in links.iter().enumerate() {
            // The link's fraction q is ½ but for an interpolated link.
            let (kind, from, group, q2, du, frac) = match link.rule {
                Rule::Pull { from } => (KIND_PULL, from, 0, 0, 0.0, 0.5),
                Rule::Bounce { group, du } => (KIND_BOUNCE, 0, group, 0, du, 0.5),
                // A bounce flagged INTERP, its next cell in the second word.
                Rule::Interpolated { group, slot } => {
                    let i = interpolations.get(slot as usize).ok_or_else(|| Error::Invalid(format!("interpolated wall link {slot} is missing")))?;
                    (KIND_BOUNCE, i.next, group, 0, i.du, i.q)
                }
                Rule::Reflect { from, q, group } => (KIND_REFLECT, from, group, u32::from(q), 0.0, 0.5),
                Rule::Open => (KIND_OPEN, 0, 0, 0, 0.0, 0.5),
            };
            if interpolated {
                link_q.push(frac as f32);
            }
            // What the weights add to the link's exchange (see `link_forces`).
            let weight = match link.rule {
                // The interpolated exchange's weights cancel in out − back, as the halfway one's.
                Rule::Bounce { .. } | Rule::Interpolated { .. } => c(u32::from(link.q)).map(|x| -2.0 * w(u32::from(link.q)) * x),
                Rule::Reflect { .. } => std::array::from_fn(|a| (c(q2)[a] - c(u32::from(link.q))[a]) * w(q2)),
                _ => [0.0; 3],
            };
            if link.rule.wall_group().is_some() {
                by_group.push((group, l as u32));
                if let Some(g) = weight_forces.get_mut(group as usize) {
                    for (s, x) in g.iter_mut().zip(weight) {
                        *s += x;
                    }
                }
            }
            let reset = if thermal.is_some_and(|t| t.resets(link.cell)) { RESET } else { 0 };
            let interp = if matches!(link.rule, Rule::Interpolated { .. }) { INTERP } else { 0 };
            link_u.extend_from_slice(&[u32::from(link.q) | (kind << 8) | (q2 << 16) | reset | interp, from, group, link.cell]);
            link_f.push(du as f32);
        }
        if interpolated {
            link_f.extend_from_slice(&link_q);
            let velocities = flow.wall_velocities();
            for g in 0..group_count {
                let u = velocities.get(g).copied().unwrap_or([0.0; 3]);
                link_f.extend_from_slice(&[u[0] as f32, u[1] as f32, u[2] as f32, 0.0]);
            }
        }
        // The boundary nodes' links: the weights' share of what a node takes is −⅓ n (the
        // populations arriving from the fluid, −n/6, less those sent back, n/6).
        // On any other node, −2 Σ w_q c_q over its links to the fluid.
        for (&(group, l), node) in node_links.iter().zip(flow.wall_nodes().map_or(&[][..], |n| n.nodes())) {
            by_group.push((group, l));
            if let Some(g) = weight_forces.get_mut(group as usize) {
                if node.flat {
                    let axis = usize::from(node.axis).min(2);
                    g[axis] -= f64::from(node.inward.signum()) / 3.0;
                } else {
                    for (q, role) in node.roles.iter().enumerate() {
                        if matches!(role, fieldcraft_lbm::bounce::Role::Han | fieldcraft_lbm::bounce::Role::Edge) {
                            for (a, x) in g.iter_mut().enumerate() {
                                *x -= 2.0 * w(q as u32) * c(q as u32)[a];
                            }
                        }
                    }
                }
            }
        }
        by_group.sort();
        let mut group_start = filled(group_count + 1, 0u32, "force groups")?;
        for (group, _) in &by_group {
            if let Some(next) = group_start.get_mut(*group as usize + 1) {
                *next += 1;
            }
        }
        for g in 0..group_count {
            let prev = group_start.get(g).copied().unwrap_or(0);
            if let Some(next) = group_start.get_mut(g + 1) {
                *next += prev;
            }
        }
        let by_group: Vec<u32> = by_group.iter().map(|(_, l)| *l).collect();
        // Each group's thermal condition (held temperatures shifted by T₀).
        let mut group_t = reserved(4 * group_count, "thermal conditions")?;
        for g in 0..group_count {
            let (kind, value) = match thermal.and_then(|t| t.groups().get(g)).copied().unwrap_or_default() {
                ThermalCondition::Adiabatic {} => (0.0, 0.0),
                ThermalCondition::Temperature { value } => (1.0, (value - t0) as f32),
                ThermalCondition::Flux { value } => (2.0, value as f32),
            };
            group_t.extend_from_slice(&[kind, value, 0.0, 0.0]);
        }

        // Pressure cells each have a slot in the outlet buffer: convective history and
        // transient gathered populations for uniform-pressure extrapolation (see lbm.wgsl).
        let open = flow.open_cells();
        let mut open_u = reserved(4 * open.len(), "open cells")?;
        let mut open_f = reserved(4 * open.len(), "open cells")?;
        let outlet_words = outlet_words(open.iter().filter(|oc| !matches!(oc.target, Target::Velocity(_))).count(), cap)?;
        let mut outlet = reserved(outlet_words, "outlet cells")?;
        let states = flow.outlets().states();
        for (k, oc) in open.iter().enumerate() {
            let (kind, target) = match oc.target {
                Target::Velocity(u) => (0u32, [u[0] as f32, u[1] as f32, u[2] as f32]),
                Target::Density(rho) => (1u32, [rho as f32, 0.0, 0.0]),
                Target::Convective(rho) => (2u32, [rho as f32, 0.0, 0.0]),
            };
            let slot = u32::try_from(outlet.len() / OUTLET).map_err(|_| Error::Limit("too many open cells for the GPU solver".into()))?;
            open_u.extend_from_slice(&[oc.cell, oc.axis as u32, kind, slot]);
            open_f.extend_from_slice(&[target[0], target[1], target[2], oc.inward as f32]);
            if let Target::Convective(_) = oc.target {
                let state = states.get(k).copied().unwrap_or_default();
                let dirs = fieldcraft_lbm::outlet::unknowns(oc.axis, oc.inward);
                outlet.extend(state.unknowns.iter().zip(dirs).map(|(v, q)| (v - W.get(q).copied().unwrap_or(0.0)) as f32));
                outlet.push(state.speed as f32);
                outlet.extend([0.0; OUTLET - UNKNOWN - 1]);
            } else if matches!(oc.target, Target::Density(_)) {
                // Transient pressure populations, gathered before any open cell is changed.
                outlet.extend([0.0; OUTLET]);
            }
        }
        let prepare_outlets = !outlet.is_empty();
        let mut open_t = reserved(4 * open.len(), "open cell temperatures")?;
        for k in 0..open.len() {
            match thermal.and_then(|t| t.open().get(k)).copied().flatten() {
                Some(value) => open_t.extend_from_slice(&[1.0, (value - t0) as f32, 0.0, 0.0]),
                None => open_t.extend_from_slice(&[0.0; 4]),
            }
        }

        let (omega_plus, omega_minus) = flow.omegas();
        let regularized = flow.collision().is_regularized();
        let acc = flow.acceleration();
        let forced = acc != [0.0; 3];
        let cells = groups(n)?;
        let links_groups = groups(links.len())?;
        let open_groups = groups(open.len())?;
        let force_groups = groups((group_count * WG as usize).max(1))?;
        let partial_count = (cells.0 as usize) * (cells.1 as usize);
        let mut params = Vec::with_capacity(64);
        let u = |v: usize| u32::try_from(v).map_err(|_| Error::Limit("too many cells for the GPU solver".into()));
        for v in [
            u(n)?,
            u(flow.grid().cells[0])?,
            u(flow.grid().cells[1])?,
            u(flow.grid().cells[2])?,
            MAX_GROUPS_X,
            u(links.len())?,
            u(open.len())?,
            u(group_count)?,
        ] {
            params.extend_from_slice(&v.to_le_bytes());
        }
        for v in [
            omega_plus as f32,
            omega_minus as f32,
            (1.0 - 0.5 * omega_plus) as f32,
            (1.0 - 0.5 * omega_minus) as f32,
            acc[0] as f32,
            acc[1] as f32,
            acc[2] as f32,
            if forced { 1.0 } else { 0.0 },
        ] {
            params.extend_from_slice(&v.to_le_bytes());
        }
        let buoyancy = thermal.map_or([0.0; 3], |t| t.buoyancy());
        for v in [
            buoyancy[0] as f32,
            buoyancy[1] as f32,
            buoyancy[2] as f32,
            if buoyancy != [0.0; 3] { 1.0 } else { 0.0 },
            thermal.map_or(0.0, |t| t.omegas().1) as f32,
            t0 as f32,
            if thermal.is_some() { 1.0 } else { 0.0 },
            thermal.map_or(0.0, |t| t.omegas().0) as f32,
        ] {
            params.extend_from_slice(&v.to_le_bytes());
        }

        // Buffers, with the device's out-of-memory errors caught.
        let scopes = [wgpu::ErrorFilter::OutOfMemory, wgpu::ErrorFilter::Validation].map(|f| device.push_error_scope(f));
        let storage = |label: &str, bytes: &[u8]| {
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(label),
                contents: if bytes.is_empty() { &[0u8; 16][..] } else { bytes },
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
            })
        };
        let empty = |label: &str, size: u64| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: size.max(16),
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            })
        };
        let solid = thermal.and_then(|t| t.solid()).map(|sh| solid_data(sh, t0)).transpose()?;
        let solid_params = solid.as_ref().map_or([0; 4], |sd| sd.params);
        // Smagorinsky's model (lbm.wgsl's `smag`): C², 1 / Pr_t, Λ, ν.
        let c = flow.smagorinsky();
        let prandtl = flow.turbulent_prandtl();
        let smag = [
            (c * c) as f32,
            if flow.turbulent() && thermal.is_some() && prandtl > 0.0 { (1.0 / prandtl) as f32 } else { 0.0 },
            ((1.0 / omega_plus - 0.5) * (1.0 / omega_minus - 0.5)) as f32,
            ((1.0 / omega_plus - 0.5) / 3.0) as f32,
        ];
        let damping = damping_data(flow, n)?;
        let damping_count = damping.as_ref().map_or(0, |(_, distances)| distances.len());
        let a_plus = flow.walls().and_then(|w| w.a_plus()).unwrap_or(1.0);
        let wall_params = [
            u32::try_from(wall_count).map_err(|_| Error::Limit("too many wall cells for the GPU solver".into()))?,
            0,
            u32::try_from(damping_count).map_err(|_| Error::Limit("too many wall cells for the GPU solver".into()))?,
            ((1.0 / a_plus) as f32).to_bits(),
        ];
        // The hybrid collision's fixed σ and Vreman's c (`les`), and the periodic axes.
        // A fixed σ of 1 is the recursive regularised collision exactly: it compiles to those
        // kernels (no blend, and no velocity pass without Vreman's model), so the GPU's `hrr` at
        // σ 1 is its `regularized` bit for bit (the hybrid kernels skip the blend there too, but
        // the driver rounds the rest differently: 1.5e-8 of 0.03 on the GB10).
        use fieldcraft_model::HrrSigma;
        let hybrid = flow.hybrid_sigma().filter(|s| *s != HrrSigma::Value(1.0));
        let fixed_sigma = match hybrid {
            Some(HrrSigma::Value(s)) => s,
            _ => 1.0,
        };
        // Dynamic σ's cap on the von Kármán length: the box's shortest extent in cells.
        let l_max = flow.grid().cells.iter().copied().min().unwrap_or(1) as f32;
        let les = [fixed_sigma as f32, flow.vreman() as f32, l_max, 0.0];
        let periodic = flow.periodic().map(u32::from);
        let mut sa_wall_rows = 0u32;
        let mut sa_wall_total = 0u32;
        if sa_faces.is_some() && on_node {
            let nodes = flow.wall_nodes().map_or(&[][..], |n| n.nodes());
            let second = flow.wall_nodes().is_some_and(|n| n.second_cell());
            let mut passes: Vec<std::collections::BTreeMap<u32, (f32, Vec<u32>)>> = vec![std::collections::BTreeMap::new(); 1 + usize::from(second)];
            for (k, w) in nodes.iter().enumerate() {
                if !w.sampled() {
                    continue;
                }
                if let (Some(first), Some(n)) = (passes.first_mut(), flow.wall_nodes()) {
                    first.entry(w.interior).or_insert_with(|| (n.distance() as f32, Vec::new())).1.push(k as u32);
                }
                if second
                    && w.sample_cell != w.interior
                    && let (Some(extra), Some(n)) = (passes.get_mut(1), flow.wall_nodes())
                {
                    extra.entry(w.sample_cell).or_insert_with(|| (n.law_distance(w) as f32, Vec::new())).1.push(k as u32);
                }
            }
            for (pass, by_cell) in passes.iter().enumerate() {
                for (cell, (y, nodes)) in by_cell {
                    if nodes.len() > 6 {
                        return Err(Error::Invalid(format!("cell {cell} is the first fluid node of {} boundary nodes (at most 6)", nodes.len())));
                    }
                    let mut row = [0u32; 8];
                    row[0] = *cell;
                    row[1] = nodes.len() as u32;
                    for (slot, k) in row.iter_mut().skip(2).zip(nodes) {
                        *slot = *k;
                    }
                    wall_rows.extend(row.map(f32::from_bits));
                    wall_rows.extend([0.0, 0.0, *y, 0.0]);
                }
                let n = u32::try_from(by_cell.len()).map_err(|_| Error::Limit("too many first fluid nodes for the GPU solver".into()))?;
                if pass == 0 {
                    sa_wall_rows = n;
                }
                sa_wall_total = sa_wall_total.saturating_add(n);
            }
        }
        let sa_law = sa_faces.is_some() && on_node && flow.wall_nodes().is_some_and(|n| n.law() == fieldcraft_lbm::sa::WallLaw::SpalartAllmaras);
        if sa_law {
            let segments = fieldcraft_lbm::sa_law::segments();
            let count = u32::try_from(segments.len()).map_err(|_| Error::Limit("the wall law's table is too long for the GPU solver".into()))?;
            wall_rows.extend([f32::from_bits(count), 0.0, 0.0, 0.0]);
            for row in &segments {
                wall_rows.extend(row.map(|v| v as f32));
            }
        }
        let sa_inflow = flow.spalart_allmaras().map_or(0.0, |f| f.inflow());
        let params = [0u32, 1].map(|phase| {
            let mut bytes = params.clone();
            bytes.extend_from_slice(&bytes_u32(&[phase, zones.zones, zones.fans, zones.fan_table]));
            bytes.extend_from_slice(&bytes_u32(&solid_params));
            bytes.extend_from_slice(&bytes_f32(&smag));
            bytes.extend_from_slice(&bytes_u32(&wall_params));
            bytes.extend_from_slice(&bytes_f32(&les));
            bytes.extend_from_slice(&bytes_u32(&[periodic[0], periodic[1], periodic[2], 0]));
            bytes.extend_from_slice(&bytes_u32(&[sa_faces.unwrap_or(0), (sa_inflow as f32).to_bits(), 1.0f32.to_bits(), 1]));
            bytes.extend_from_slice(&bytes_u32(&[sa_wall_total, 0, 0, 0]));
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("lbm params"),
                contents: &bytes,
                usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            })
        });
        let fluid_buf = storage("lbm fluid", &bytes_u32(&fluid_u32));
        let walls_buf = storage("lbm walls", &bytes_f32(&wall_rows));
        let damping_buf = storage("lbm wall damping", &bytes_f32(damping.as_ref().map_or(&[][..], |(data, _)| &data[..])));
        // Only the wall cells' distances stay on the host.
        let damping = damping.map(|(_, distances)| distances);
        // AA reuses binding g for boundary scratch only, never another full population field.
        let pops = [empty("lbm f", pop_bytes), empty("lbm g / AA boundary scratch", if aa { 4 * links.len() as u64 } else { pop_bytes })];
        // Populations as deviations (see lbm.wgsl), one direction at a time: `value(q, cell)` less
        // `offset(q)`, zero in solid cells. Each direction's write is submitted and waited for
        // before the next: a write's staging copy lives until its submission completes, so the
        // directions would otherwise pile up into a whole copy of the populations.
        let mut direction = filled(n, 0.0f32, "populations")?;
        // (The on-node bounce's boundary nodes hold populations too.)
        let holds = |cell: usize| fluid_u32.get(cell).is_some_and(|u| *u != 0);
        let mut upload = |buffer: &wgpu::Buffer, count: usize, offset: &dyn Fn(usize) -> f64, value: &dyn Fn(usize, usize) -> f64| -> Result<()> {
            for q in 0..count {
                let w = offset(q);
                for (cell, slot) in direction.iter_mut().enumerate() {
                    *slot = if holds(cell) { (value(q, cell) - w) as f32 } else { 0.0 };
                }
                queue.write_buffer(buffer, 4 * (q * n) as u64, &bytes_f32(&direction));
                queue.submit([]);
                device
                    .poll(wgpu::PollType::Wait { submission_index: None, timeout: Some(DEVICE_TIMEOUT) })
                    .map_err(|e| gpu_error(format!("the device did not answer within {} s: {e}", DEVICE_TIMEOUT.as_secs())))?;
            }
            Ok(())
        };
        upload(&pops[0], Q, &|q| fieldcraft_lbm::lattice::W.get(q).copied().unwrap_or(0.0), &|q, cell| flow.population(q, cell))?;
        let link_u = storage("lbm links", &bytes_u32(&link_u));
        let link_f = storage("lbm link du", &bytes_f32(&link_f));
        let open_u = storage("lbm open cells", &bytes_u32(&open_u));
        let open_f = storage("lbm open targets", &bytes_f32(&open_f));
        let outlet = storage("lbm outlet cells", &bytes_f32(&outlet));
        let moments = empty("lbm moments", 16 * n as u64);
        // With the hybrid collision or Vreman's model, every cell's velocity before the collision
        // (three f32 per cell) for the finite-difference gradient.
        let gradients = hybrid.is_some() || flow.vreman() > 0.0 || sa_faces.is_some();
        let velocity = empty("lbm gradient velocities", if gradients { 12 * n as u64 } else { 16 });
        let stats = storage("lbm stats", &bytes_u32(&[0; 8]));
        let (nu_tilde, sa_rates, sa_cell) = match (flow.spalart_allmaras(), flow.sa_wall_distance()) {
            (Some(field), Some(distance)) if sa_faces.is_some() => {
                let mut v = filled(n, 0.0f32, "Spalart–Allmaras ν̃")?;
                for (out, x) in v.iter_mut().zip(field.nu_tilde()) {
                    *out = *x as f32;
                }
                let nu_tilde = storage("lbm sa nu_tilde", &bytes_f32(&v));
                drop(v);
                let mut cell = filled(2 * n, 0.0f32, "Spalart–Allmaras wall distances")?;
                for (out, d) in cell.as_chunks_mut::<2>().0.iter_mut().zip(distance) {
                    out[1] = *d as f32;
                }
                let sa_cell = storage("lbm sa cells", &bytes_f32(&cell));
                (nu_tilde, empty("lbm sa rates", 8 * n as u64), sa_cell)
            }
            _ => (empty("lbm sa nu_tilde", 16), empty("lbm sa rates", 16), empty("lbm sa cells", 16)),
        };
        let contrib = empty("lbm link forces", 16 * (links.len() as u64 + u64::from(node_slots)));
        let by_group = storage("lbm links by group", &bytes_u32(&by_group));
        let group_start = storage("lbm group start", &bytes_u32(&group_start));
        let forces = empty("lbm forces", 16 * group_count as u64);
        // Also each zone's sample (4 floats) and each fan workgroup's flow (`fan_partials`).
        let partial_slots = partial_count.max(4 * zones.zones as usize).max(zones.fan_groups as usize);
        let partials = empty("lbm partials", 4 * partial_slots as u64);
        let therm_bytes = elem * (QT * n) as u64;
        let therm = if thermal.is_some() {
            let tb = [empty("lbm t", therm_bytes), empty("lbm t2", therm_bytes)];
            // The scalar populations as deviations from w_q T₀ (see lbm.wgsl).
            upload(&tb[0], QT, &|q| WT.get(q).copied().unwrap_or(0.0) * t0, &|q, cell| flow.thermal_population(q, cell))?;
            tb
        } else {
            [empty("lbm t", 16), empty("lbm t2", 16)]
        };
        let group_t = storage("lbm thermal conditions", &bytes_f32(&group_t));
        let open_t = storage("lbm open temperatures", &bytes_f32(&open_t));
        let (solid_u, solid_f, solid_s) = solid.as_ref().map_or((&[][..], &[][..], &[][..]), |sd| (&sd.u[..], &sd.f[..], &sd.state[..]));
        let solid_u_buf = storage("lbm solid data", &bytes_u32(solid_u));
        let solid_f_buf = storage("lbm solid values", &bytes_f32(solid_f));
        let solid_s_buf = storage("lbm solid state", &bytes_f32(solid_s));

        let (kernels, pipeline_seconds) = shared.kernels(Variant {
            aa,
            regularized,
            mrt: flow.collision() == fieldcraft_model::Collision::Mrt,
            smagorinsky: flow.smagorinsky() > 0.0,
            wall_model: wall_count > 0 && flow.wall_model(),
            van_driest: damping_count > 0,
            wall_bounce: wall_count > 0 && flow.wall_bounce(),
            on_node,
            hrr: hybrid.is_some(),
            vreman: flow.vreman() > 0.0,
            dynamic: hybrid == Some(HrrSigma::Dynamic),
            zones: zones.zones > 0,
            interpolated,
            sa: sa_faces.is_some(),
            sa_law,
        })?;
        // Every binding's buffer in bind group k (the two differ in the phase and the A–B swap).
        let resource = |binding: u32, k: usize| -> Option<&wgpu::Buffer> {
            Some(match binding {
                0 => &params[k],
                1 => &fluid_buf,
                2 => &pops[if aa { 0 } else { k }],
                3 => &pops[if aa { 1 } else { 1 - k }],
                4 => &outlet,
                5 => &link_u,
                6 => &link_f,
                7 => &open_u,
                8 => &open_f,
                9 => &moments,
                10 => &stats,
                11 => &contrib,
                12 => &by_group,
                13 => &group_start,
                14 => &forces,
                15 => &partials,
                16 => &therm[k],
                17 => &therm[1 - k],
                18 => &group_t,
                19 => &open_t,
                20 => &solid_u_buf,
                21 => &solid_f_buf,
                22 => &solid_s_buf,
                23 => &walls_buf,
                24 => &damping_buf,
                25 => &velocity,
                // The flags buffer again, writable, for `fan_update` alone (lbm.wgsl).
                26 => &fluid_buf,
                27 => &nu_tilde,
                28 => &sa_rates,
                29 => &sa_cell,
                _ => return None,
            })
        };
        let mut bind = HashMap::with_capacity(kernels.len());
        for (name, kernel) in kernels.iter() {
            let group = |k: usize| -> Result<wgpu::BindGroup> {
                let entries = kernel
                    .bindings
                    .iter()
                    .map(|&b| {
                        resource(b, k)
                            .map(|buffer| wgpu::BindGroupEntry { binding: b, resource: buffer.as_entire_binding() })
                            .ok_or_else(|| gpu_error(format!("{name}: no buffer for binding {b}")))
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(device.create_bind_group(&wgpu::BindGroupDescriptor { label: Some(name), layout: &kernel.layout, entries: &entries }))
            };
            bind.insert(*name, [group(0)?, group(1)?]);
        }
        for scope in scopes.into_iter().rev() {
            if let Some(e) = pollster::block_on(scope.pop()) {
                return Err(gpu_error(e));
            }
        }
        Ok(GpuFlow {
            device: device.clone(),
            queue: queue.clone(),
            kernel_times: profile.then(Default::default),
            pipeline_seconds,
            aa,
            adapter: adapter_name,
            kernels,
            bind,
            current: 0,
            n,
            cells,
            link_count: links.len(),
            links_groups,
            wall_groups: if wall_count > 0 { Some(groups(wall_count)?) } else { None },
            wall_bounce: bounce,
            on_node,
            gradients,
            walls: walls_buf,
            wall_count,
            damping: match damping {
                Some(distances) => Some(Damping { buffer: damping_buf, groups: groups(distances.len())?, distances, nu: nu0 }),
                None => None,
            },
            params,
            open_count: open.len(),
            open_groups,
            prepare_outlets,
            group_count,
            force_groups,
            stats,
            moments,
            forces,
            partials,
            partial_count,
            zones,
            fan_groups: groups(zones.fan_groups as usize * WG as usize)?,
            thermal: thermal.map(|_| therm),
            initial_temperature: t0,
            fluid_cells: flow.fluid_cells(),
            fluid,
            groups: flow.groups().to_vec(),
            empty_obstacles: flow.empty_obstacles().to_vec(),
            wall_layer: fieldcraft_lbm::backend::Backend::wall_layer(flow),
            wall_links: flow.wall_links(),
            steps_done: 0,
            diagnostic_acceleration: flow.acceleration().map(|a| f64::from(a as f32)),
            nu0: (1.0 / omega_plus - 0.5) / 3.0,
            eddy_ratio: 0.0,
            exchange_cache: Some(vec![[0.0; 4]; group_count]),
            weight_forces,
            solid: match solid {
                Some(sd) => Some(Solids {
                    state: solid_s_buf,
                    cell_groups: groups(sd.cells.len())?,
                    link_groups: groups(sd.params[1] as usize)?,
                    cells: sd.cells,
                    by_solid: sd.by_solid,
                    substeps: sd.substeps,
                }),
                None => None,
            },
            sa: sa_faces.map(|_| SaDevice { nu_tilde, substeps: None, fixed: false }),
            sa_cells: sa_faces.map(|_| sa_cell),
            sa_wall_rows,
            sa_wall_groups: if sa_wall_total > 0 { Some(groups(sa_wall_total as usize)?) } else { None },
        })
    }

    /// With the Spalart–Allmaras model: take exactly `n` sub-steps of ν̃ per lattice step from now
    /// on (1 to `fieldcraft_lbm::sa::MAX_SUBSTEPS`), as the CPU's [`Flow::set_sa_substeps`] does,
    /// for parity with it; a chunk whose coefficients break the monotone bound then stops the run.
    /// `None` plans the count per chunk again.
    #[doc(hidden)]
    pub fn fix_sa_substeps(&mut self, n: Option<u32>) -> Result<()> {
        let Some(sa) = &mut self.sa else { return Err(Error::Invalid("the flow has no Spalart–Allmaras model".into())) };
        if let Some(n) = n.filter(|n| !(1..=SA_MAX_SUBSTEPS).contains(n)) {
            return Err(Error::Invalid(format!("Spalart–Allmaras: {n} sub-steps; give 1 to {SA_MAX_SUBSTEPS}")));
        }
        sa.fixed = n.is_some();
        sa.substeps = n;
        Ok(())
    }

    /// With the Spalart–Allmaras model: whether its sources act (on by default); off, ν̃ is only
    /// transported, as the CPU's [`Flow::set_sa_model`].
    #[doc(hidden)]
    pub fn set_sa_model(&mut self, on: bool) -> Result<()> {
        if self.sa.is_none() {
            return Err(Error::Invalid("the flow has no Spalart–Allmaras model".into()));
        }
        for p in &self.params {
            self.queue.write_buffer(p, SA_PARAMS + 12, &u32::from(on).to_le_bytes());
        }
        Ok(())
    }

    /// With the Spalart–Allmaras model, the wall distance each cell holds on the device (lattice
    /// units, uploaded once as f32).
    pub fn sa_wall_distance(&self) -> Result<Option<Vec<f64>>> {
        if self.sa.is_none() {
            return Ok(None);
        }
        let Some(cells) = self.sa_cells.as_ref() else { return Ok(None) };
        let v = self.read_f32(cells, 2 * self.n)?;
        Ok(Some(v.as_chunks::<2>().0.iter().map(|c| f64::from(c[1])).collect()))
    }

    /// With the Spalart–Allmaras model on Han's wall, the wall function's samples of the last step
    /// as the CPU's [`Flow::sa_wall_function`] gives them: each first fluid node's cell, its
    /// distance to the wall and the friction velocity it read (lattice units), in cell order.
    /// Empty without them.
    pub fn sa_wall_function(&self) -> Result<Vec<(usize, f64, f64)>> {
        if self.sa_wall_rows == 0 {
            return Ok(Vec::new());
        }
        let base = (WALL_ROW + 8) * self.wall_count;
        let rows = self.read_f32(&self.walls, base + 12 * self.sa_wall_rows as usize)?;
        Ok(rows
            .get(base..)
            .unwrap_or_default()
            .as_chunks::<12>()
            .0
            .iter()
            .map(|r| (r[0].to_bits() as usize, f64::from(r[10]), f64::from(r[8])))
            .collect())
    }

    /// With the Spalart–Allmaras model, the sub-steps per lattice step the next chunk takes (none
    /// before the first chunk plans them).
    pub fn sa_substeps(&self) -> Option<u32> {
        self.sa.as_ref().and_then(|s| s.substeps)
    }

    /// Write the sub-step dt = 1/n into both phases' uniforms.
    fn set_sa_dt(&self, n: u32) {
        let dt = (1.0 / f64::from(n.max(1))) as f32;
        for p in &self.params {
            self.queue.write_buffer(p, SA_DT, &dt.to_bits().to_le_bytes());
        }
    }

    /// Before the first chunk: the largest ν̃ coefficient sum of the initial field (one rate pass
    /// on the current velocities), so the first chunk's count is planned like the others'.
    fn probe_sa(&mut self) -> Result<f64> {
        self.queue.write_buffer(&self.stats, 0, &bytes_u32(&[0; 8]));
        self.submit("lbm sa probe", &[(self.current, "gradient_velocity", self.cells), (self.current, "sa_rate", self.cells)])?;
        let stats = self.read(&self.stats, 32)?;
        Ok(f64::from(f32::from_bits(stats.get(16..20).and_then(|b| b.try_into().ok()).map(u32::from_le_bytes).unwrap_or(0))))
    }

    /// The sub-step count for a largest coefficient sum `most`, as SaField::advance bounds it,
    /// with [`SA_MARGIN`]; the CPU's error past [`SA_MAX_SUBSTEPS`].
    fn plan_sa(most: f64) -> Result<u32> {
        if !most.is_finite() {
            return Err(Error::Unstable("Spalart–Allmaras: ν̃'s transport coefficients are not finite".into()));
        }
        if most > f64::from(SA_MAX_SUBSTEPS) {
            return Err(Error::Unstable(format!(
                "Spalart–Allmaras: ν̃'s transport needs more than {SA_MAX_SUBSTEPS} sub-steps per lattice step: refine the grid or lower the inflow ν̃"
            )));
        }
        Ok(((SA_MARGIN * most).ceil() as u32).clamp(1, SA_MAX_SUBSTEPS))
    }

    /// The adapter's name.
    pub fn adapter_name(&self) -> &str {
        &self.adapter
    }

    /// Whether the populations stream in place (see [`Streaming`]).
    pub fn streams_in_place(&self) -> bool {
        self.aa
    }

    /// How long compiling the kernels took for this flow, in seconds: zero when this process had
    /// compiled them on this device before.
    pub fn pipeline_seconds(&self) -> f64 {
        self.pipeline_seconds
    }

    /// With [`Options::profile`], the GPU time of each kernel so far, by kernel name.
    pub fn kernel_times(&self) -> Option<std::collections::BTreeMap<&'static str, KernelTime>> {
        self.kernel_times.as_ref().map(|t| t.borrow().clone())
    }

    fn dispatch(&self, pass: &mut wgpu::ComputePass<'_>, bind: usize, name: &str, (gx, gy): (u32, u32)) -> Result<()> {
        let kernel = self.kernels.get(name).ok_or_else(|| gpu_error(format!("no kernel {name}")))?;
        let group = self.bind.get(name).and_then(|b| b.get(bind)).ok_or_else(|| gpu_error(format!("no bind group for {name}")))?;
        pass.set_bind_group(0, group, &[]);
        pass.set_pipeline(&kernel.pipeline);
        pass.dispatch_workgroups(gx, gy, 1);
        Ok(())
    }

    /// Timestamped passes are diagnostic: their extra pass boundaries and readback are not
    /// used for throughput measurements. Normally every chunk remains one compute pass.
    fn submit(&self, label: &str, commands: &[Dispatch]) -> Result<()> {
        let mut enc = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some(label) });
        let Some(times) = &self.kernel_times else {
            {
                let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor { label: Some(label), timestamp_writes: None });
                for &(bind, name, groups) in commands {
                    self.dispatch(&mut pass, bind, name, groups)?;
                }
            }
            self.queue.submit([enc.finish()]);
            return Ok(());
        };
        let count = u32::try_from(commands.len() * 2).map_err(|_| gpu_error("too many timestamp queries"))?;
        if count == 0 || count > wgpu::QUERY_SET_MAX_QUERIES {
            return Err(gpu_error("timestamp query count is outside the device limit"));
        }
        let query = self.device.create_query_set(&wgpu::QuerySetDescriptor { label: Some(label), ty: wgpu::QueryType::Timestamp, count });
        let result = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("lbm timestamps"),
            size: u64::from(count) * 8,
            usage: wgpu::BufferUsages::QUERY_RESOLVE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        for (i, &(bind, name, groups)) in commands.iter().enumerate() {
            let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some(name),
                timestamp_writes: Some(wgpu::ComputePassTimestampWrites {
                    query_set: &query,
                    beginning_of_pass_write_index: Some(2 * i as u32),
                    end_of_pass_write_index: Some(2 * i as u32 + 1),
                }),
            });
            self.dispatch(&mut pass, bind, name, groups)?;
        }
        enc.resolve_query_set(&query, 0..count, &result, 0);
        self.queue.submit([enc.finish()]);
        let bytes = self.read(&result, u64::from(count) * 8)?;
        let period = f64::from(self.queue.get_timestamp_period());
        let mut totals = times.borrow_mut();
        for ((_, name, _), pair) in commands.iter().zip(bytes.as_chunks::<16>().0) {
            let [begin, end] = pair.as_chunks::<8>().0 else { return Err(gpu_error("invalid timestamp pair")) };
            let ticks = u64::from_le_bytes(*end).checked_sub(u64::from_le_bytes(*begin)).ok_or_else(|| gpu_error("invalid GPU timestamps"))?;
            let total = totals.entry(name).or_default();
            total.dispatches += 1;
            total.seconds += ticks as f64 * period * NANOSECOND;
        }
        Ok(())
    }

    /// Copy `size` bytes of `src` back to the host.
    fn read(&self, src: &wgpu::Buffer, size: u64) -> Result<Vec<u8>> {
        let size = size.max(16);
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("lbm readback"),
            size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let mut enc = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("lbm readback") });
        enc.copy_buffer_to_buffer(src, 0, &staging, 0, size);
        let copied = self.queue.submit([enc.finish()]);
        let slice = staging.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        // A bounded wait: a hung or lost device is an error, never a hang. The device is shared by
        // the solves of the process ([`DEVICES`]): another solve's poll can collect this mapping's
        // callback and run it just after this poll returns (#800), so the callback is waited for
        // too, polling again until it has run. Each poll waits for this copy's submission only, not
        // for whatever another solve submitted since.
        let deadline = std::time::Instant::now() + DEVICE_TIMEOUT;
        loop {
            let left = deadline.saturating_duration_since(std::time::Instant::now());
            if let Err(e) = self.device.poll(wgpu::PollType::Wait { submission_index: Some(copied.clone()), timeout: Some(left) }) {
                return Err(gpu_error(format!("the device did not answer within {} s: {e}", DEVICE_TIMEOUT.as_secs())));
            }
            match rx.recv_timeout(CALLBACK_WAIT.min(left)) {
                Ok(mapped) => break mapped.map_err(gpu_error)?,
                Err(std::sync::mpsc::RecvTimeoutError::Timeout) if std::time::Instant::now() < deadline => {}
                Err(_) => return Err(gpu_error("the readback did not complete")),
            }
        }
        let bytes = slice.get_mapped_range().map_err(gpu_error)?.to_vec();
        staging.unmap();
        Ok(bytes)
    }

    fn read_f32(&self, src: &wgpu::Buffer, count: usize) -> Result<Vec<f32>> {
        let bytes = self.read(src, 4 * count as u64)?;
        Ok(bytes.as_chunks::<4>().0.iter().take(count).map(|b| f32::from_le_bytes(*b)).collect())
    }

    /// The momentum (xyz) and heat (w) each group received in the last step.
    fn exchange(&mut self) -> Result<&[[f64; 4]]> {
        if self.exchange_cache.is_none() {
            let v = self.read_f32(&self.forces, 4 * self.group_count)?;
            self.exchange_cache = Some(
                v.as_chunks::<4>()
                    .0
                    .iter()
                    .zip(&self.weight_forces)
                    .map(|(g, w)| [f64::from(g[0]) + w[0], f64::from(g[1]) + w[1], f64::from(g[2]) + w[2], f64::from(g[3])])
                    .collect(),
            );
        }
        Ok(self.exchange_cache.as_deref().unwrap_or_default())
    }

    /// The wall model's last sample: each modelled wall cell's grid cell, its modelled velocity
    /// and relaxation rate (lattice units), in the order of [`fieldcraft_lbm::wall::Walls::wall_cells`]
    /// that have a second cell. Empty without the wall model.
    pub fn wall_targets(&self) -> Result<Vec<([f64; 3], f64)>> {
        if self.wall_count == 0 {
            return Ok(Vec::new());
        }
        let rows = self.read_f32(&self.walls, WALL_ROW * self.wall_count)?;
        Ok(rows.as_chunks::<WALL_ROW>().0.iter().map(|r| ([f64::from(r[8]), f64::from(r[9]), f64::from(r[10])], f64::from(r[11]))).collect())
    }

    /// Stored population momentum before any wall-target substitution. Diagnostic only;
    /// currently restricted to closed, unheated boxes without zones or on-node boundaries.
    /// Uses the same current AA/AB orientation as the normal moments readback.
    pub fn population_momentum(&self) -> Result<[f64; 3]> {
        if self.thermal.is_some() || self.zones.zones != 0 || self.on_node || self.open_count != 0 {
            return Err(Error::Invalid("population momentum diagnostic supports only closed, unheated boxes without zones or on-node walls".into()));
        }
        self.compute_moments()?;
        let values = self.read_f32(&self.moments, 4 * self.n)?;
        let mut momentum = [0.0; 3];
        for (cell, row) in values.as_chunks::<4>().0.iter().enumerate() {
            if self.fluid.get(cell).copied().unwrap_or(false) {
                for (a, value) in momentum.iter_mut().enumerate() {
                    *value += f64::from(row[3]) * (f64::from(row[a]) - 0.5 * self.diagnostic_acceleration[a]);
                }
            }
        }
        if !momentum.iter().all(|v| v.is_finite()) {
            return Err(Error::Invalid("non-finite population momentum".into()));
        }
        Ok(momentum)
    }
    /// The moments of the current state into the moments buffer.
    fn compute_moments(&self) -> Result<()> {
        self.submit("lbm moments", &[(self.current, "cell_moments", self.cells)])
    }
}

impl GpuFlow {
    /// The most dispatches one time step takes (the flow's, the scalar lattice's and the
    /// solids' with their sub-steps).
    fn dispatches_per_step(&self) -> usize {
        let solid = self.solid.as_ref().map_or(0, |sd| 3 + sd.substeps as usize);
        let sa = self.sa.as_ref().map_or(0, |s| 1 + 2 * s.substeps.unwrap_or(SA_MAX_SUBSTEPS) as usize);
        13 + solid + sa + if self.zones.fans > 0 { 2 } else { 0 }
    }

    /// One submission of `steps` (≤ MAX_CHUNK) time steps, then the stability check.
    fn advance_chunk(&mut self, steps: u64) -> Result<()> {
        let sa_substeps = match self.sa.as_ref().map(|s| s.substeps) {
            None => 0,
            Some(Some(n)) => n,
            Some(None) => {
                let n = GpuFlow::plan_sa(self.probe_sa()?)?;
                if let Some(sa) = &mut self.sa {
                    sa.substeps = Some(n);
                }
                n
            }
        };
        if sa_substeps > 0 {
            self.set_sa_dt(sa_substeps);
        }
        self.queue.write_buffer(&self.stats, 0, &bytes_u32(&[0; 8]));
        let mut commands = Vec::with_capacity(self.dispatches_per_step() * steps as usize + 7);
        for step in 0..steps {
            // Each fan's flow from the populations and its previous push, then its rise moved
            // toward its curve, before anything else in the step reads the push (as on the CPU).
            if self.zones.fans > 0 {
                if self.zones.fan_groups > 0 {
                    commands.push((self.current, "fan_partials", self.fan_groups));
                }
                commands.push((self.current, "fan_update", (self.zones.fans, 1)));
            }
            if step + 1 == steps {
                commands.push((self.current, "clear_eddy", (1, 1)));
            }
            if let Some(d) = &self.damping {
                commands.push((self.current, "damping_refresh", d.groups));
            }
            if let Some(groups) = self.wall_groups {
                commands.push((self.current, if self.wall_bounce || self.on_node { "bounce_sample" } else { "wall_sample" }, groups));
            }
            if self.gradients {
                commands.push((self.current, "gradient_velocity", self.cells));
            }
            if sa_substeps > 0 {
                commands.push((self.current, "sa_prepare", self.cells));
                for _ in 0..sa_substeps {
                    commands.push((self.current, "sa_rate", self.cells));
                    commands.push((self.current, "sa_update", self.cells));
                }
                if let Some(groups) = self.sa_wall_groups {
                    commands.push((self.current, "sa_wall_function", groups));
                }
            }
            commands.push((self.current, "collide", self.cells));
            if let Some(groups) = self.wall_groups.filter(|_| self.on_node) {
                commands.push((self.current, "node_bounce", groups));
            }
            if !self.aa {
                commands.push((self.current, "stream", self.cells));
            }
            if self.link_count > 0 {
                // The last step's momentum exchange is the one read back.
                let exchange = self.aa && step + 1 == steps;
                commands.push((self.current, if exchange { "stream_links_exchange" } else { "stream_links" }, self.links_groups));
                if self.aa {
                    commands.push((self.current, "scatter_links", self.links_groups));
                }
            }
            if self.prepare_outlets {
                commands.push((self.current, "outlet_convect", self.open_groups));
            }
            if self.open_count > 0 {
                commands.push((self.current, "open_boundaries", self.open_groups));
            }
            if self.thermal.is_some() {
                commands.push((self.current, "collide_thermal", self.cells));
                commands.push((self.current, "stream_thermal", self.cells));
                if self.link_count > 0 {
                    commands.push((self.current, "stream_links_thermal", self.links_groups));
                }
                if let Some(sd) = &self.solid {
                    commands.push((self.current, "solid_ab", sd.cell_groups));
                    for k in 0..sd.substeps {
                        commands.push((self.current, if k % 2 == 0 { "solid_sub0" } else { "solid_sub1" }, sd.cell_groups));
                    }
                    if sd.substeps % 2 == 1 {
                        commands.push((self.current, "solid_copy", sd.cell_groups));
                    }
                    commands.push((self.current, "solid_scatter", sd.link_groups));
                }
                // After the interface, as on the CPU: an inlet holds its temperature on its
                // links into a solid too, and an outlet copies its neighbour's final populations.
                if self.open_count > 0 {
                    commands.push((self.current, "open_thermal", self.open_groups));
                }
            }
            self.current ^= 1;
        }
        // The momentum (and heat) exchange of the last step, from its post-collision
        // populations; with solids also without fluid links, for their heat into box faces.
        if (self.link_count > 0 || self.solid.is_some()) && self.group_count > 0 {
            commands.push((self.current ^ 1, "link_forces", self.links_groups));
            if let Some(groups) = self.wall_groups {
                commands.push((self.current ^ 1, "wall_forces", groups));
            }
            if let Some(sd) = &self.solid {
                commands.push((self.current ^ 1, "solid_contrib", sd.link_groups));
            }
            commands.push((self.current ^ 1, "group_forces", self.force_groups));
            if self.solid.is_some() {
                commands.push((self.current ^ 1, "solid_wall_heat", (6, 1)));
            }
        }
        if self.thermal.is_some() {
            commands.push((self.current, "thermal_check", self.cells));
        }
        self.submit("lbm steps", &commands)?;
        self.steps_done += steps;
        self.exchange_cache = None;
        let stats = self.read(&self.stats, 32)?;
        let word = |k: usize| stats.get(4 * k..4 * k + 4).and_then(|b| b.try_into().ok()).map(u32::from_le_bytes).unwrap_or(0);
        let max_u2 = f32::from_bits(word(0));
        if word(1) != 0 || !(f64::from(max_u2) <= MAX_LATTICE_SPEED * MAX_LATTICE_SPEED) {
            return Err(Error::Unstable(format!(
                "the flow went unstable by time step {}: the lattice speed reached {:.3} (limit {MAX_LATTICE_SPEED}) or a density stopped being a positive number; lower the time step, refine the spacing or check the boundary conditions",
                self.steps_done,
                max_u2.sqrt()
            )));
        }
        if word(2) != 0 {
            return Err(Error::Unstable(format!(
                "the heat went unstable by time step {}: a temperature stopped being a number; lower the time step, refine the spacing or check the heat fluxes",
                self.steps_done
            )));
        }
        // The last step's largest eddy viscosity, as on the CPU.
        if self.nu0 > 0.0 {
            self.eddy_ratio = f64::from(f32::from_bits(word(3))) / self.nu0;
        }
        if let Some(sa) = &mut self.sa {
            if word(5) != 0 {
                return Err(Error::Unstable(format!("Spalart–Allmaras: ν̃ stopped being a finite number by time step {}", self.steps_done)));
            }
            let most = f64::from(f32::from_bits(word(4)));
            let worst = most / f64::from(sa_substeps.max(1));
            if sa.fixed {
                if !(worst <= 1.0) {
                    return Err(Error::Unstable(format!(
                        "Spalart–Allmaras: {sa_substeps} fixed sub-steps break ν̃'s monotone bound, which needs {} (by time step {})",
                        most.ceil(),
                        self.steps_done
                    )));
                }
            } else {
                if !(worst <= SA_STOP) {
                    return Err(Error::Unstable(format!(
                        "Spalart–Allmaras: ν̃'s transport went past its stability bound by time step {} ({worst:.2} times the monotone step with {sa_substeps} sub-steps): refine the grid or lower the inflow ν̃",
                        self.steps_done
                    )));
                }
                sa.substeps = Some(GpuFlow::plan_sa(most)?);
            }
        }
        Ok(())
    }
}

impl Backend for GpuFlow {
    fn eddy_viscosity_ratio(&self) -> f64 {
        self.eddy_ratio
    }

    fn nu_tilde(&mut self) -> Result<Option<Vec<f64>>> {
        let Some(sa) = &self.sa else { return Ok(None) };
        Ok(Some(self.read_f32(&sa.nu_tilde, self.n)?.into_iter().map(f64::from).collect()))
    }

    fn wall_y_plus(&self) -> Option<f64> {
        let d = self.damping.as_ref()?;
        // The friction velocities lead the buffer; a device that does not answer reports none.
        let ut = self.read_f32(&d.buffer, d.distances.len()).ok()?;
        let y = d.distances.iter().zip(ut).map(|(y, ut)| y * f64::from(ut) / d.nu).fold(0.0f64, f64::max);
        Some(if y.is_finite() { y } else { 0.0 })
    }

    fn set_acceleration(&mut self, a: [f64; 3]) -> Result<()> {
        if !a.iter().all(|v| v.is_finite()) {
            return Err(Error::Invalid(format!("the body acceleration {a:?} is not a number")));
        }
        // `acc` follows the eight counts and the four rates: 48 bytes in.
        let acc = [a[0] as f32, a[1] as f32, a[2] as f32, if a != [0.0; 3] { 1.0 } else { 0.0 }];
        self.diagnostic_acceleration = a.map(|v| f64::from(v as f32));
        for p in &self.params {
            self.queue.write_buffer(p, 48, &bytes_f32(&acc));
        }
        Ok(())
    }

    fn name(&self) -> String {
        format!("gpu: {}, {} streaming", self.adapter, if self.aa { "in-place" } else { "two-buffer" })
    }

    fn advance(&mut self, steps: u64) -> Result<()> {
        // Submissions of at most MAX_CHUNK steps: the stability check runs after each, so a
        // run that blows up is caught within that many steps whatever the caller asks for.
        // With timestamps every dispatch takes two queries of one query set.
        let most = if self.kernel_times.is_some() {
            ((u64::from(wgpu::QUERY_SET_MAX_QUERIES) / 2).saturating_sub(6) / self.dispatches_per_step() as u64).clamp(1, MAX_CHUNK)
        } else {
            MAX_CHUNK
        };
        let mut left = steps;
        while left > 0 {
            let chunk = left.min(most);
            self.advance_chunk(chunk)?;
            left -= chunk;
        }
        Ok(())
    }

    fn steps_done(&self) -> u64 {
        self.steps_done
    }

    fn fields(&mut self) -> Result<Fields> {
        self.compute_moments()?;
        // The wall model's cells at the velocity their reconstruction imposes (as the CPU's).
        if self.wall_count > 0 && self.steps_done > 0 {
            self.submit("lbm wall moments", &[(self.current, "wall_moments", self.cells)])?;
        }
        let m = self.read_f32(&self.moments, 4 * self.n)?;
        let density = m.as_chunks::<4>().0.iter().map(|v| f64::from(v[3])).collect();
        let velocity = m.as_chunks::<4>().0.iter().map(|v| [f64::from(v[0]), f64::from(v[1]), f64::from(v[2])]).collect();
        // The temperature from the scalar populations, summed on the host: T₀ + Σ t (T₀ in
        // solid cells, whose populations stay zero).
        let temperature = match &self.thermal {
            None => Vec::new(),
            Some(therm) => {
                let t = self.read_f32(therm.get(self.current).ok_or_else(|| gpu_error("no thermal buffer"))?, QT * self.n)?;
                let mut temperature: Vec<f64> = (0..self.n)
                    .map(|cell| self.initial_temperature + (0..QT).map(|q| t.get(q * self.n + cell).copied().map_or(0.0, f64::from)).sum::<f64>())
                    .collect();
                // The conducting cells' own temperatures.
                if let Some(sd) = &self.solid {
                    let theta = self.read_f32(&sd.state, sd.cells.len())?;
                    for (cell, th) in sd.cells.iter().zip(theta) {
                        if let Some(slot) = temperature.get_mut(*cell as usize) {
                            *slot = self.initial_temperature + f64::from(th);
                        }
                    }
                }
                temperature
            }
        };
        Ok(Fields { density, velocity, temperature, fill: Vec::new() })
    }

    fn heat(&mut self) -> Result<Vec<f64>> {
        if self.thermal.is_none() {
            return Ok(Vec::new());
        }
        Ok(self.exchange()?.iter().map(|g| g[3]).collect())
    }

    fn mass(&mut self) -> Result<f64> {
        self.compute_moments()?;
        self.submit("lbm mass", &[(self.current, "mass", self.cells)])?;
        let partials = self.read_f32(&self.partials, self.partial_count)?;
        // Summed in a fixed order on the host, so the result is deterministic on this device.
        Ok(partials.iter().map(|v| f64::from(*v)).sum())
    }

    fn forces(&mut self) -> Result<Vec<[f64; 3]>> {
        Ok(self.exchange()?.iter().map(|g| [g[0], g[1], g[2]]).collect())
    }

    fn mean_velocity(&mut self, axis: usize) -> Result<f64> {
        let axis = u32::try_from(axis).ok().filter(|a| *a < 3).ok_or_else(|| Error::Invalid(format!("axis {axis} is not x, y or z")))?;
        // `wall.y` follows the eight counts, the five vectors before `smag` and `smag` itself and
        // the wall cells' count: 148 bytes in.
        for p in &self.params {
            self.queue.write_buffer(p, 148, &bytes_u32(&[axis]));
        }
        self.compute_moments()?;
        self.submit("lbm mean velocity", &[(self.current, "mean_velocity", self.cells)])?;
        let partials = self.read_f32(&self.partials, self.partial_count)?;
        // Summed in a fixed order on the host, so the result is deterministic on this device.
        Ok(partials.iter().map(|v| f64::from(*v)).sum::<f64>() / self.fluid_cells.max(1) as f64)
    }

    fn zones(&mut self) -> Result<Vec<fieldcraft_lbm::ZoneSample>> {
        let count = self.zones.zones as usize;
        if count == 0 {
            return Ok(Vec::new());
        }
        self.compute_moments()?;
        self.submit("lbm zone samples", &[(self.current, "zone_samples", (self.zones.zones, 1))])?;
        let samples = self.read_f32(&self.partials, 4 * count)?;
        Ok(samples
            .as_chunks::<4>()
            .0
            .iter()
            .map(|v| fieldcraft_lbm::ZoneSample { flow: f64::from(v[0]), pressure: f64::from(v[1]), fan_rise: f64::from(v[2]) })
            .collect())
    }

    fn groups(&self) -> &[Group] {
        &self.groups
    }

    fn fluid_cells(&self) -> usize {
        self.fluid_cells
    }

    fn is_fluid(&self, cell: usize) -> bool {
        self.fluid.get(cell).copied().unwrap_or(false)
    }

    fn empty_obstacles(&self) -> &[String] {
        &self.empty_obstacles
    }

    fn wall_links(&self) -> Option<fieldcraft_lbm::WallLinks> {
        self.wall_links
    }

    fn wall_layer(&self) -> Option<fieldcraft_lbm::wall_layer::WallLayerReport> {
        self.wall_layer
    }

    fn conducting(&self) -> &[Vec<u32>] {
        self.solid.as_ref().map_or(&[], |s| &s.by_solid)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pressure_staging_refuses_oversized_buffers_before_allocation() {
        assert_eq!(outlet_words(2, 160).unwrap(), 40);
        assert!(outlet_words(2, 159).unwrap_err().to_string().contains("pressure outlets"));
        assert!(outlet_words(usize::MAX, u64::MAX).is_err());
    }

    /// #640: each kernel binds only the buffers it uses, so the solver fits GL's 16 storage
    /// buffers per compute stage, and WebGPU's default of 8.
    #[test]
    fn every_kernel_binds_at_most_eight_storage_buffers() {
        let r = reflection().unwrap();
        assert_eq!(r.kernels.len(), KERNELS.len());
        assert!(r.kernels.iter().all(|(_, b)| !b.is_empty()));
        println!("the busiest kernel binds {} storage buffers", r.storage_buffers);
        assert!((1..=8).contains(&r.storage_buffers), "{}", r.storage_buffers);
    }

    /// A boundary node's side row holds the CPU's roles exactly: bit q of the Han mask is a Han
    /// link, of the fluid mask a Han or Edge link (what `on_links` counts in the momentum the node
    /// takes), and the flags say flat, sampled and the axes it carries, for every normal (and
    /// none, a concave corner) and a spread of neighbourhoods.
    #[test]
    fn a_node_s_side_row_packs_the_cpu_s_roles() {
        use fieldcraft_lbm::bounce::{Neighbour, Role, WallNode, carries, roles};
        let mut seed = 0x2545_f491_4f6c_dd1du64;
        let mut checked = 0;
        for normal in [None, Some((0, 1)), Some((0, -1)), Some((1, 1)), Some((1, -1)), Some((2, 1)), Some((2, -1))] {
            for pattern in 0..64 {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                let kinds: [Neighbour; Q] = std::array::from_fn(|q| match (seed >> (2 * q)) & 3 {
                    _ if pattern == 0 => Neighbour::Fluid,
                    0 | 1 => Neighbour::Fluid,
                    2 => Neighbour::Node,
                    _ => Neighbour::Wall,
                });
                let r = roles(normal, |q| kinds[q]);
                let (axis, inward) = normal.unwrap_or((1, 1));
                let w = WallNode::new(7, 9, axis as u8, inward as i8, 0, r);
                let (han, to_fluid, flags) = node_roles(&w);
                for (q, role) in r.iter().enumerate().skip(1) {
                    assert_eq!(han >> q & 1 == 1, *role == Role::Han, "{normal:?} {pattern} q {q}");
                    assert_eq!(to_fluid >> q & 1 == 1, matches!(role, Role::Han | Role::Edge), "{normal:?} {pattern} q {q}");
                }
                assert_eq!(han & 1, 0);
                assert_eq!(flags & 1 == 1, w.flat);
                assert_eq!(flags >> 1 & 1 == 1, w.sampled());
                for b in 0..3 {
                    assert_eq!(flags >> (2 + b) & 1 == 1, carries(&r, axis, inward, b), "{normal:?} {pattern} axis {b}");
                }
                checked += 1;
            }
        }
        println!("{checked} neighbourhoods packed as the CPU's roles");
    }

    /// The wall function's constants in the shader are the CPU's (κ, B, c_v1), and its f32
    /// algorithm (`sa_exp_tail`, `sa_spalding_nu_t`, `sa_nu_tilde_for`, ported here line by line)
    /// gives the CPU's first-node ν̃ within 2e-5 relative over y⁺ from 0.1 to 3000 (measured:
    /// 9.7e-7), the viscous sublayer included, where the series keeps e^x − 1 − x − x²/2 from
    /// cancelling.
    #[test]
    fn the_wall_function_s_f32_roots_follow_the_cpu() {
        use fieldcraft_lbm::sa::{C_V1, WallLaw, law_eddy_viscosity, nu_tilde_for};
        use fieldcraft_lbm::wall::{BOUNCE_B, BOUNCE_KAPPA};
        let constant = |name: &str| -> f64 {
            let at = SHADER.find(&format!("const {name}: f32 = ")).unwrap() + name.len() + 14;
            SHADER[at..].split(';').next().unwrap().trim().parse().unwrap()
        };
        let (law_e, kappa, cv1) = (constant("SA_LAW_E"), constant("SA_KAPPA"), constant("SA_CV1"));
        assert!(((law_e - (-BOUNCE_KAPPA * BOUNCE_B).exp()) / law_e).abs() < 1e-12, "{law_e}");
        assert_eq!(kappa, BOUNCE_KAPPA);
        assert_eq!(cv1, C_V1);
        let (e, k, cv1) = (law_e as f32, kappa as f32, cv1 as f32);
        let tail = |x: f32, m: u32| -> f32 {
            if x < 0.5 {
                let x2 = x * x;
                let x4 = x2 * x2;
                let mut s = x4 / 24.0 + x4 * x / 120.0 + x4 * x2 / 720.0 + x4 * x2 * x / 5040.0 + x4 * x4 / 40320.0 + x4 * x4 * x / 362880.0;
                if m == 3 {
                    s += x2 * x / 6.0;
                }
                return s;
            }
            let mut s = x.exp() - 1.0 - x - 0.5 * x * x;
            if m == 4 {
                s -= x * x * x / 6.0;
            }
            s
        };
        let bisect = |mut lo: f32, mut hi: f32, below: &dyn Fn(f32) -> bool| -> f32 {
            for _ in 0..40 {
                let mid = 0.5 * (lo + hi);
                if hi - lo <= 1e-7 * hi || !(mid > lo && mid < hi) {
                    break;
                }
                if below(mid) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            hi
        };
        let fv1 = |chi: f32| {
            let c3 = chi * chi * chi;
            c3 / (c3 + 357.911)
        };
        let nu = 1e-4f32;
        let mut worst = 0.0f64;
        let mut y = 0.1f64;
        while y < 3000.0 {
            let yp = y as f32;
            let u = bisect(0.0, yp, &|u| u + e * tail(k * u, 4) < yp);
            let nu_t_plus = (k * e * tail(k * u, 3)).max(0.0);
            let t = nu_t_plus;
            let device = if t > 0.0 { f64::from(bisect(t, t + cv1, &|chi| chi * fv1(chi) < t) * nu) } else { 0.0 };
            let host = nu_tilde_for(f64::from(nu) * law_eddy_viscosity(WallLaw::Spalding, y).unwrap(), f64::from(nu)).unwrap();
            let rel = (device - host).abs() / host.max(1e-30);
            assert!(rel < 2e-5, "y⁺ {y}: device {device:e} host {host:e}");
            worst = worst.max(rel);
            y *= 1.05;
        }
        println!("the device's first-node ν̃ within {worst:.1e} of the CPU's over y⁺ 0.1–3000");
    }

    /// An adapter that allows fewer is refused before anything is built, saying what to do.
    #[test]
    fn an_adapter_with_too_few_storage_buffers_is_refused_clearly() {
        assert!(check_storage_buffers(wgpu::Backend::Gl, 16, 7).is_ok());
        let e = check_storage_buffers(wgpu::Backend::Gl, 4, 7).unwrap_err();
        assert!(matches!(&e, Error::Backend(_)));
        let e = e.to_string();
        println!("{e}");
        for part in ["OpenGL backend allows 4 storage buffers", "needs 7", "backend: \"cpu\"", "WGPU_BACKEND"] {
            assert!(e.contains(part), "{e}");
        }
    }

    /// The shader flattens every dispatch with one row width: a dispatch of several rows must be
    /// exactly [`MAX_GROUPS_X`] wide, whatever the cell count (a link dispatch wider than the
    /// cells' row once overlapped itself and skipped links), and cover its items.
    #[test]
    fn every_dispatch_of_several_rows_is_max_groups_x_wide() {
        let row = (MAX_GROUPS_X * WG) as usize;
        for count in [0, 1, 255, 256, 257, row - 1, row, row + 1, 3 * row - 7, 5_120_000, 15_360_000, 190_000_000] {
            let (gx, gy) = groups(count).unwrap();
            assert!(gy == 1 || gx == MAX_GROUPS_X, "{count}: {gx} × {gy}");
            assert!((gx * gy * WG) as usize >= count, "{count}: {gx} × {gy}");
        }
    }
}
