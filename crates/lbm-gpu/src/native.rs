//! The device backend (native targets).

use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

use fieldcraft_gpu::select::{self, Policy};
use fieldcraft_lbm::backend::{Backend, Fields};
use fieldcraft_lbm::lattice::{C, Q, W};
use fieldcraft_lbm::outlet::UNKNOWN;
use fieldcraft_lbm::setup::MAX_LATTICE_SPEED;
use fieldcraft_lbm::{Error, Flow, Group, Result, Rule, Target};
use fieldcraft_model::ThermalCondition;
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
/// Floats per open cell in the outlet buffer (see lbm.wgsl).
const OUTLET: usize = 16;
/// The scalar lattice's directions (the first seven of D3Q19) and weights.
const QT: usize = 7;
const WT: [f64; QT] = [0.25, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125];

fn gpu_error(e: impl std::fmt::Display) -> Error {
    Error::Invalid(format!("GPU: {e}"))
}

/// Steps per submission: the stability check and a cancel take effect between submissions.
const MAX_CHUNK: u64 = 64;

/// How long one wait for the device may take before the solve gives up on it.
const DEVICE_TIMEOUT: std::time::Duration = std::time::Duration::from_secs(120);

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
    open_count: usize,
    open_groups: (u32, u32),
    /// Whether a face is a convective outlet.
    convective: bool,
    group_count: usize,
    force_groups: (u32, u32),
    stats: wgpu::Buffer,
    moments: wgpu::Buffer,
    forces: wgpu::Buffer,
    partials: wgpu::Buffer,
    partial_count: usize,
    /// The scalar populations remain A–B with heat; and T₀ they are shifted by.
    thermal: Option<[wgpu::Buffer; 2]>,
    initial_temperature: f64,
    fluid: Vec<bool>,
    fluid_cells: usize,
    groups: Vec<Group>,
    empty_obstacles: Vec<String>,
    steps_done: u64,
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
const KERNELS: [&str; 23] = [
    "collide",
    "stream",
    "stream_links",
    "stream_links_exchange",
    "scatter_links",
    "outlet_convect",
    "open_boundaries",
    "link_forces",
    "group_forces",
    "cell_moments",
    "mass",
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
];

/// The bindings each kernel uses, read from the shader (naga's analysis of which globals an entry
/// point reaches): each kernel gets a bind group of only those. One layout of all 23 buffers
/// needs 22 storage buffers per compute stage, more than GL's 16 (#640); the busiest kernel
/// needs 7, within WebGPU's default of 8.
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
    Err(Error::Invalid(format!(
        "this GPU's {backend} backend allows {have} storage buffers per compute stage and the fluid solver needs {need}: solve with `backend: \"cpu\"`, or use a Vulkan, Metal or DX12 adapter (`WGPU_BACKEND`)"
    )))
}

/// The override constants a flow's kernels are compiled with.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Variant {
    aa: bool,
    regularized: bool,
    smagorinsky: bool,
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
                    ("SMAGORINSKY", if variant.smagorinsky { 1.0 } else { 0.0 }),
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

/// The kept device for `adapter` with `features`, or a new one, opened with the adapter's buffer
/// sizes and `storage_buffers` per compute stage.
fn shared_device(adapter: &wgpu::Adapter, features: wgpu::Features, storage_buffers: u32) -> Result<Arc<Shared>> {
    let info = adapter.get_info();
    let key = format!("{:?} {} {:04x}:{:04x} {:?} {} {features:?}", info.backend, info.name, info.vendor, info.device, info.device_type, info.driver);
    let mut devices = DEVICES.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    devices.retain(|d| !d.lost.load(Ordering::Acquire));
    if let Some(d) = devices.iter().find(|d| d.key == key) {
        return Ok(d.clone());
    }
    let have = adapter.limits();
    let mut limits = wgpu::Limits::default().using_resolution(have.clone());
    limits.max_buffer_size = have.max_buffer_size;
    limits.max_storage_buffer_binding_size = have.max_storage_buffer_binding_size;
    // Exactly what the busiest kernel binds, so every device (the CI's lavapipe too) proves the
    // kernels fit an adapter that allows no more.
    limits.max_storage_buffers_per_shader_stage = storage_buffers;
    let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
        label: Some("fieldcraft lbm"),
        required_features: features,
        required_limits: limits,
        ..Default::default()
    }))
    .map_err(gpu_error)?;
    let lost = Arc::new(AtomicBool::new(false));
    let flag = lost.clone();
    device.set_device_lost_callback(move |reason, message| {
        log::warn!("the fluid solver's GPU device was lost ({reason:?}): {message}");
        flag.store(true, Ordering::Release);
    });
    let shared = Arc::new(Shared { key, device, queue, lost, kernels: Mutex::new(Vec::new()) });
    devices.push(shared.clone());
    Ok(shared)
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
        if let Some(zone) = flow.zones().first() {
            return Err(Error::Invalid(format!(
                "the GPU backend does not step fans or porous zones yet (`{}`): solve on the CPU (`backend: \"cpu\"`)",
                zone.name
            )));
        }
        let Some((adapter, tier)) = select::adapter(policy).map_err(gpu_error)? else {
            return Err(Error::Invalid(format!("{}=cpu: no GPU device; solve with `backend: \"cpu\"`", select::ENV)));
        };
        let info = adapter.get_info();
        let profile = options.profile;
        let features = if profile { wgpu::Features::TIMESTAMP_QUERY } else { wgpu::Features::empty() };
        if !adapter.features().contains(features) {
            return Err(gpu_error("profiling needs timestamp queries, which this adapter does not offer"));
        }
        let storage_buffers = reflection()?.storage_buffers;
        check_storage_buffers(info.backend, adapter.limits().max_storage_buffers_per_shader_stage, storage_buffers)?;
        let shared = shared_device(&adapter, features, storage_buffers)?;
        let (device, queue) = (&shared.device, &shared.queue);
        let adapter_name = format!("{}{}", info.name, if tier == select::Tier::Software { " (software)" } else { "" });

        let aa = match options.streaming {
            Streaming::Auto | Streaming::InPlace => true,
            Streaming::TwoBuffers => false,
        };
        let n = flow.grid().len();
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

        // The state and the flags. The populations go up one direction at a time (below), so
        // the host never holds a whole single-precision copy of them beside the flow's.
        let thermal = flow.thermal();
        let t0 = thermal.map_or(0.0, |t| t.initial());
        let mut fluid = filled(n, false, "cell flags")?;
        let mut fluid_u32 = filled(n, 0u32, "cell flags")?;
        for cell in 0..n {
            if let (Some(b), Some(u)) = (fluid.get_mut(cell), fluid_u32.get_mut(cell)) {
                *b = flow.is_fluid(cell);
                *u = u32::from(*b);
            }
        }

        // The links, by cell, and the ones with a force by group.
        let links = flow.links();
        let mut link_u = reserved(4 * links.len(), "links")?;
        let mut link_f = reserved(links.len(), "links")?;
        let mut by_group: Vec<(u32, u32)> = reserved(links.len(), "links by group")?;
        let group_count = flow.groups().len();
        let mut weight_forces = filled(group_count, [0.0f64; 3], "force groups")?;
        let c = |q: u32| C.get(q as usize).copied().unwrap_or([0; 3]).map(f64::from);
        let w = |q: u32| W.get(q as usize).copied().unwrap_or(0.0);
        for (l, link) in links.iter().enumerate() {
            let (kind, from, group, q2, du) = match link.rule {
                Rule::Pull { from } => (KIND_PULL, from, 0, 0, 0.0),
                Rule::Bounce { group, du } => (KIND_BOUNCE, 0, group, 0, du),
                Rule::Reflect { from, q, group } => (KIND_REFLECT, from, group, u32::from(q), 0.0),
                Rule::Open => (KIND_OPEN, 0, 0, 0, 0.0),
            };
            // What the weights add to the link's exchange (see `link_forces`).
            let weight = match link.rule {
                Rule::Bounce { .. } => c(u32::from(link.q)).map(|x| -2.0 * w(u32::from(link.q)) * x),
                Rule::Reflect { .. } => std::array::from_fn(|a| (c(q2)[a] - c(u32::from(link.q))[a]) * w(q2)),
                _ => [0.0; 3],
            };
            if matches!(link.rule, Rule::Bounce { .. } | Rule::Reflect { .. }) {
                by_group.push((group, l as u32));
                if let Some(g) = weight_forces.get_mut(group as usize) {
                    for (s, x) in g.iter_mut().zip(weight) {
                        *s += x;
                    }
                }
            }
            let reset = if thermal.is_some_and(|t| t.resets(link.cell)) { RESET } else { 0 };
            link_u.extend_from_slice(&[u32::from(link.q) | (kind << 8) | (q2 << 16) | reset, from, group, link.cell]);
            link_f.push(du as f32);
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

        // The open cells, and what each convective outlet cell keeps between steps in its slot
        // (`open_u.w`) of the outlet buffer: its unknown populations as deviations from their
        // weights, and its convection speed (see lbm.wgsl).
        let open = flow.open_cells();
        let mut open_u = reserved(4 * open.len(), "open cells")?;
        let mut open_f = reserved(4 * open.len(), "open cells")?;
        let mut outlet = reserved(OUTLET * open.iter().filter(|oc| matches!(oc.target, Target::Convective(_))).count(), "outlet cells")?;
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
            }
        }
        let convective = !outlet.is_empty();
        let mut open_t = reserved(4 * open.len(), "open cell temperatures")?;
        for k in 0..open.len() {
            match thermal.and_then(|t| t.open().get(k)).copied().flatten() {
                Some(value) => open_t.extend_from_slice(&[1.0, (value - t0) as f32, 0.0, 0.0]),
                None => open_t.extend_from_slice(&[0.0; 4]),
            }
        }

        let (omega_plus, omega_minus) = flow.omegas();
        let regularized = flow.collision() == fieldcraft_model::Collision::Regularized;
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
        // Smagorinsky's model (lbm.wgsl's `smag`): C², 1 / Pr_t, Λ.
        let c = flow.smagorinsky();
        let prandtl = flow.turbulent_prandtl();
        let smag = [
            (c * c) as f32,
            if c > 0.0 && thermal.is_some() && prandtl > 0.0 { (1.0 / prandtl) as f32 } else { 0.0 },
            ((1.0 / omega_plus - 0.5) * (1.0 / omega_minus - 0.5)) as f32,
            0.0,
        ];
        let params = [0u32, 1].map(|phase| {
            let mut bytes = params.clone();
            bytes.extend_from_slice(&bytes_u32(&[phase, 0, 0, 0]));
            bytes.extend_from_slice(&bytes_u32(&solid_params));
            bytes.extend_from_slice(&bytes_f32(&smag));
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("lbm params"),
                contents: &bytes,
                usage: wgpu::BufferUsages::UNIFORM,
            })
        });
        let fluid_buf = storage("lbm fluid", &bytes_u32(&fluid_u32));
        // AA reuses binding g for boundary scratch only, never another full population field.
        let pops = [empty("lbm f", pop_bytes), empty("lbm g / AA boundary scratch", if aa { 4 * links.len() as u64 } else { pop_bytes })];
        // Populations as deviations (see lbm.wgsl), one direction at a time: `value(q, cell)` less
        // `offset(q)`, zero in solid cells. Each direction's write is submitted and waited for
        // before the next: a write's staging copy lives until its submission completes, so the
        // directions would otherwise pile up into a whole copy of the populations.
        let mut direction = filled(n, 0.0f32, "populations")?;
        let mut upload = |buffer: &wgpu::Buffer, count: usize, offset: &dyn Fn(usize) -> f64, value: &dyn Fn(usize, usize) -> f64| -> Result<()> {
            for q in 0..count {
                let w = offset(q);
                for (cell, slot) in direction.iter_mut().enumerate() {
                    *slot = if flow.is_fluid(cell) { (value(q, cell) - w) as f32 } else { 0.0 };
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
        let stats = storage("lbm stats", &bytes_u32(&[0, 0, 0, 0]));
        let contrib = empty("lbm link forces", 16 * links.len() as u64);
        let by_group = storage("lbm links by group", &bytes_u32(&by_group));
        let group_start = storage("lbm group start", &bytes_u32(&group_start));
        let forces = empty("lbm forces", 16 * group_count as u64);
        let partials = empty("lbm partials", 4 * partial_count as u64);
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

        let (kernels, pipeline_seconds) = shared.kernels(Variant { aa, regularized, smagorinsky: flow.smagorinsky() > 0.0 })?;
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
            open_count: open.len(),
            open_groups,
            convective,
            group_count,
            force_groups,
            stats,
            moments,
            forces,
            partials,
            partial_count,
            thermal: thermal.map(|_| therm),
            initial_temperature: t0,
            fluid_cells: flow.fluid_cells(),
            fluid,
            groups: flow.groups().to_vec(),
            empty_obstacles: flow.empty_obstacles().to_vec(),
            steps_done: 0,
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
        })
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
        let kernel = self.kernels.get(name).ok_or_else(|| Error::Invalid(format!("GPU: no kernel {name}")))?;
        let group = self.bind.get(name).and_then(|b| b.get(bind)).ok_or_else(|| Error::Invalid(format!("GPU: no bind group for {name}")))?;
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
        self.queue.submit([enc.finish()]);
        let slice = staging.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        // A bounded wait: a hung or lost device is an error, never a hang.
        match self.device.poll(wgpu::PollType::Wait { submission_index: None, timeout: Some(DEVICE_TIMEOUT) }) {
            Ok(_) => {}
            Err(e) => return Err(gpu_error(format!("the device did not answer within {} s: {e}", DEVICE_TIMEOUT.as_secs()))),
        }
        rx.try_recv().map_err(|_| gpu_error("the readback did not complete"))?.map_err(gpu_error)?;
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
        9 + solid
    }

    /// One submission of `steps` (≤ MAX_CHUNK) time steps, then the stability check.
    fn advance_chunk(&mut self, steps: u64) -> Result<()> {
        self.queue.write_buffer(&self.stats, 0, &bytes_u32(&[0, 0, 0, 0]));
        let mut commands = Vec::with_capacity(self.dispatches_per_step() * steps as usize + 6);
        for step in 0..steps {
            commands.push((self.current, "collide", self.cells));
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
            if self.convective {
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
        let stats = self.read(&self.stats, 16)?;
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
        // The chunk's largest eddy viscosity (the CPU reports its last step's).
        if self.nu0 > 0.0 {
            self.eddy_ratio = f64::from(f32::from_bits(word(3))) / self.nu0;
        }
        Ok(())
    }
}

impl Backend for GpuFlow {
    fn eddy_viscosity_ratio(&self) -> f64 {
        self.eddy_ratio
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
        let m = self.read_f32(&self.moments, 4 * self.n)?;
        let density = m.as_chunks::<4>().0.iter().map(|v| f64::from(v[3])).collect();
        let velocity = m.as_chunks::<4>().0.iter().map(|v| [f64::from(v[0]), f64::from(v[1]), f64::from(v[2])]).collect();
        // The temperature from the scalar populations, summed on the host: T₀ + Σ t (T₀ in
        // solid cells, whose populations stay zero).
        let temperature = match &self.thermal {
            None => Vec::new(),
            Some(therm) => {
                let t = self.read_f32(therm.get(self.current).ok_or_else(|| Error::Invalid("GPU: no thermal buffer".into()))?, QT * self.n)?;
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
        Ok(Fields { density, velocity, temperature })
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

    fn conducting(&self) -> &[Vec<u32>] {
        self.solid.as_ref().map_or(&[], |s| &s.by_solid)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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

    /// An adapter that allows fewer is refused before anything is built, saying what to do.
    #[test]
    fn an_adapter_with_too_few_storage_buffers_is_refused_clearly() {
        assert!(check_storage_buffers(wgpu::Backend::Gl, 16, 7).is_ok());
        let e = check_storage_buffers(wgpu::Backend::Gl, 4, 7).unwrap_err().to_string();
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
