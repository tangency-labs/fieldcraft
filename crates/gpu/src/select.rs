//! Which adapter draws a headless picture. `FIELDCRAFT_GPU` chooses the policy:
//!
//! - `auto` (the default) walks a hardware GPU, then a software adapter (lavapipe, llvmpipe,
//!   WARP, SwiftShader), then the CPU rasterizer;
//! - `hardware` and `software` are strict: an error when there is no such adapter;
//! - `cpu` forces the CPU rasterizer.
//!
//! Backends are Metal on macOS, DX12 and Vulkan on Windows, Vulkan and GL elsewhere;
//! `WGPU_BACKEND` narrows them (`WGPU_BACKEND=gl` forces GL). Where both Vulkan and GL are on
//! offer, Vulkan is asked first and GL is brought up only when Vulkan has no adapter of the tier
//! the policy wants: NVIDIA's EGL torn down beside a working Vulkan device in the same process
//! deadlocked its driver on a GB10 (#614). Each backend set has one `wgpu::Instance` for the
//! whole process: NVIDIA's driver stops showing the GPU to new instances once 25 of them have
//! opened a device, even after they are gone (#636), so a long session would lose the GPU at
//! its 26th solve. A discrete GPU outranks an integrated one, then a virtual one, then any
//! other; a CPU device is the software tier.
//!
//! wgpu's backend options follow its environment variables: on DX12, `WGPU_DX12_COMPILER`
//! (`fxc` or `dxc`; unset, DXC when `dxcompiler.dll` is on the DLL search path, else FXC, see
//! #561). `staticdxc` is refused: it needs wgpu's `static-dxc` feature, which FieldCraft does not
//! build, and wgpu would panic on it.

use std::sync::Mutex;

use crate::{GpuError, Result};

/// The adapter policy (`FIELDCRAFT_GPU`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Policy {
    #[default]
    Auto,
    Hardware,
    Software,
    Cpu,
}

/// The environment variable that chooses the [`Policy`].
pub const ENV: &str = "FIELDCRAFT_GPU";

impl Policy {
    pub fn parse(s: &str) -> Option<Policy> {
        Some(match s.trim().to_ascii_lowercase().as_str() {
            "" | "auto" => Policy::Auto,
            "hardware" => Policy::Hardware,
            "software" => Policy::Software,
            "cpu" => Policy::Cpu,
            _ => return None,
        })
    }

    /// The policy `FIELDCRAFT_GPU` names; [`Policy::Auto`] when it is unset.
    pub fn from_env() -> Result<Policy> {
        match std::env::var(ENV) {
            Ok(v) => Policy::parse(&v).ok_or_else(|| GpuError::Policy(v.chars().take(40).collect())),
            Err(_) => Ok(Policy::Auto),
        }
    }
}

pub use crate::Tier;

/// Set to `1` where a GPU adapter must exist (CI with lavapipe): tests that draw then fail
/// instead of skipping when none is found, so the GPU path is proven to draw.
pub const REQUIRE_ENV: &str = "FIELDCRAFT_REQUIRE_GPU";

/// Whether [`REQUIRE_ENV`] is set to `1`.
pub fn gpu_required() -> bool {
    std::env::var(REQUIRE_ENV).is_ok_and(|v| v.trim() == "1")
}

/// An adapter's tier and its rank within it (lower is better).
pub fn rank(t: wgpu::DeviceType) -> (Tier, u8) {
    match t {
        wgpu::DeviceType::DiscreteGpu => (Tier::Hardware, 0),
        wgpu::DeviceType::IntegratedGpu => (Tier::Hardware, 1),
        wgpu::DeviceType::VirtualGpu => (Tier::Hardware, 2),
        wgpu::DeviceType::Other => (Tier::Hardware, 3),
        wgpu::DeviceType::Cpu => (Tier::Software, 0),
    }
}

/// The index of the adapter `policy` picks among adapters of `types`: `Ok(None)` for the CPU
/// rasterizer, an error when a strict policy finds no adapter of its tier.
pub fn pick(types: &[wgpu::DeviceType], policy: Policy) -> Result<Option<usize>> {
    let best = |tier: Tier| types.iter().enumerate().filter(|(_, t)| rank(**t).0 == tier).min_by_key(|(_, t)| rank(**t).1).map(|(i, _)| i);
    match policy {
        Policy::Cpu => Ok(None),
        Policy::Auto => Ok(best(Tier::Hardware).or_else(|| best(Tier::Software))),
        Policy::Hardware => best(Tier::Hardware).map(Some).ok_or_else(|| GpuError::NoAdapter(format!("{ENV}=hardware: no hardware GPU adapter"))),
        Policy::Software => best(Tier::Software)
            .map(Some)
            .ok_or_else(|| GpuError::NoAdapter(format!("{ENV}=software: no software GPU adapter (lavapipe, llvmpipe, WARP)"))),
    }
}

/// The backends tried on this platform, narrowed by `WGPU_BACKEND`.
pub fn backends() -> wgpu::Backends {
    let platform = if cfg!(target_os = "macos") || cfg!(target_os = "ios") {
        wgpu::Backends::METAL
    } else if cfg!(windows) {
        wgpu::Backends::DX12 | wgpu::Backends::VULKAN
    } else {
        wgpu::Backends::VULKAN | wgpu::Backends::GL
    };
    match wgpu::Backends::from_env() {
        Some(env) => platform & env,
        None => platform,
    }
}

/// The adapter `policy` picks on this machine and its tier; `Ok(None)` for the CPU rasterizer.
pub fn adapter(policy: Policy) -> Result<Option<(wgpu::Adapter, Tier)>> {
    if policy == Policy::Cpu {
        return Ok(None);
    }
    let all = backends();
    let first = if all.contains(wgpu::Backends::VULKAN | wgpu::Backends::GL) { wgpu::Backends::VULKAN } else { all };
    let mut adapters = enumerate(first)?;
    if first != all && !settled(&device_types(&adapters), policy) {
        adapters.extend(enumerate(all - first)?);
    }
    let types = device_types(&adapters);
    let Some(i) = pick(&types, policy)? else { return Ok(None) };
    let tier = types.get(i).map_or(Tier::Hardware, |t| rank(*t).0);
    Ok((i < adapters.len()).then(|| (adapters.swap_remove(i), tier)))
}

/// The adapters of `backends`.
fn enumerate(backends: wgpu::Backends) -> Result<Vec<wgpu::Adapter>> {
    let instance = instance(backends)?;
    Ok(pollster::block_on(instance.enumerate_adapters(backends)))
}

/// The instances made so far, by backend set and backend options (`Dx12Compiler` has no
/// equality, so the options are compared by their debug form).
static INSTANCES: Mutex<Vec<(wgpu::Backends, String, wgpu::Instance)>> = Mutex::new(Vec::new());

/// Instances actually created by this module in the process. Clones/reuses do not count as
/// new instances; external viewport and driver-internal instances are outside this scope.
pub fn instances_opened() -> u64 {
    INSTANCES.lock().unwrap_or_else(|e| e.into_inner()).len() as u64
}

/// The process's instance of `backends` with the backend options the environment asks for: the
/// one made at its first use, so the driver sees one instance per backend set however many
/// devices are opened (#636).
fn instance(backends: wgpu::Backends) -> Result<wgpu::Instance> {
    let options = backend_options(wgpu::BackendOptions::from_env_or_default())?;
    let key = format!("{options:?}");
    let mut made = INSTANCES.lock().unwrap_or_else(|e| e.into_inner());
    if let Some((_, _, instance)) = made.iter().find(|(b, k, _)| *b == backends && *k == key) {
        return Ok(instance.clone());
    }
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
        backends,
        backend_options: options,
        ..wgpu::InstanceDescriptor::new_without_display_handle()
    });
    made.push((backends, key, instance.clone()));
    Ok(instance)
}

fn device_types(adapters: &[wgpu::Adapter]) -> Vec<wgpu::DeviceType> {
    adapters.iter().map(|a| a.get_info().device_type).collect()
}

/// Whether adapters of `types` already hold one of the tier `policy` wants, so the other
/// backends need not be brought up.
fn settled(types: &[wgpu::DeviceType], policy: Policy) -> bool {
    let want = if policy == Policy::Software { Tier::Software } else { Tier::Hardware };
    types.iter().any(|t| rank(*t).0 == want)
}

/// `options`, unless they ask for a DX12 compiler this build does not have.
fn backend_options(options: wgpu::BackendOptions) -> Result<wgpu::BackendOptions> {
    if matches!(options.dx12.shader_compiler, wgpu::Dx12Compiler::StaticDxc) {
        return Err(GpuError::Device("WGPU_DX12_COMPILER=staticdxc is not built in: use dxc with dxcompiler.dll, or fxc".into()));
    }
    Ok(options)
}

#[cfg(test)]
mod tests {
    use wgpu::DeviceType::*;

    use super::*;

    #[test]
    fn a_dx12_compiler_that_is_not_built_in_is_refused() {
        let with = |shader_compiler| {
            let mut o = wgpu::BackendOptions::default();
            o.dx12.shader_compiler = shader_compiler;
            backend_options(o)
        };
        assert!(with(wgpu::Dx12Compiler::Fxc).is_ok());
        assert!(with(wgpu::Dx12Compiler::Auto).is_ok());
        assert!(with(wgpu::Dx12Compiler::StaticDxc).is_err());
    }

    #[test]
    fn policies_walk_the_tiers() {
        assert_eq!(pick(&[Cpu, IntegratedGpu, DiscreteGpu, VirtualGpu], Policy::Auto), Ok(Some(2)));
        assert_eq!(pick(&[Cpu, Other, IntegratedGpu], Policy::Auto), Ok(Some(2)));
        assert_eq!(pick(&[Cpu, DiscreteGpu], Policy::Software), Ok(Some(0)));
        // Auto falls to a software adapter, then to the CPU rasterizer.
        assert_eq!(pick(&[Cpu], Policy::Auto), Ok(Some(0)));
        assert_eq!(pick(&[], Policy::Auto), Ok(None));
        assert_eq!(pick(&[DiscreteGpu], Policy::Cpu), Ok(None));
        // The strict policies refuse.
        assert!(pick(&[Cpu], Policy::Hardware).is_err());
        assert!(pick(&[DiscreteGpu], Policy::Software).is_err());
        assert!(pick(&[], Policy::Hardware).is_err());
    }

    /// #614: GL is brought up only when Vulkan has no adapter of the tier the policy wants.
    #[test]
    fn vulkan_settles_when_it_has_the_wanted_tier() {
        assert!(settled(&[Cpu, IntegratedGpu], Policy::Auto));
        assert!(settled(&[DiscreteGpu], Policy::Hardware));
        assert!(settled(&[Cpu], Policy::Software));
        // Only a software adapter: auto still asks GL for a hardware one, as before.
        assert!(!settled(&[Cpu], Policy::Auto));
        assert!(!settled(&[DiscreteGpu], Policy::Software));
        assert!(!settled(&[], Policy::Hardware));
    }

    /// #636: one instance per backend set for the process, however often adapters are asked for.
    #[test]
    fn a_backend_set_keeps_one_instance() {
        let count = |b: wgpu::Backends| INSTANCES.lock().unwrap().iter().filter(|(k, _, _)| *k == b).count();
        let b = wgpu::Backends::NOOP;
        for _ in 0..3 {
            let _ = instance(b).unwrap();
            let _ = enumerate(b).unwrap();
        }
        assert_eq!(count(b), 1);
    }

    #[test]
    fn policy_names() {
        assert_eq!(Policy::parse("Software"), Some(Policy::Software));
        assert_eq!(Policy::parse(""), Some(Policy::Auto));
        assert_eq!(Policy::parse("cpu"), Some(Policy::Cpu));
        assert_eq!(Policy::parse("gpu"), None);
    }
}
