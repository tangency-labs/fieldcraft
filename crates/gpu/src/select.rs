//! Which adapter draws a headless picture. `FIELDCRAFT_GPU` chooses the policy:
//!
//! - `auto` (the default) walks a hardware Vulkan GPU, then a software Vulkan adapter
//!   (lavapipe or SwiftShader), then the CPU rasterizer;
//! - `hardware` and `software` are strict: an error when there is no such Vulkan adapter;
//! - `cpu` forces the CPU rasterizer.
//!
//! Vulkan is the only native GPU backend. `WGPU_BACKEND` may narrow this choice but cannot
//! enable another backend. A discrete GPU outranks an integrated one, then a virtual one,
//! then any other; a CPU device is the software tier.
//!
//! One Vulkan instance serves the whole process: NVIDIA's driver stops showing the GPU to
//! new instances once 25 have opened a device (#636). The instance and the solvers' shared
//! devices keep repeated solves from recreating the driver resources.

use std::sync::OnceLock;

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
        Policy::Hardware => {
            best(Tier::Hardware).map(Some).ok_or_else(|| GpuError::NoAdapter(format!("{ENV}=hardware: no hardware Vulkan GPU adapter")))
        }
        Policy::Software => best(Tier::Software)
            .map(Some)
            .ok_or_else(|| GpuError::NoAdapter(format!("{ENV}=software: no software Vulkan GPU adapter (install lavapipe or SwiftShader)"))),
    }
}

/// Vulkan, unless `WGPU_BACKEND` explicitly excludes it. No override enables another backend.
pub fn backends() -> Result<wgpu::Backends> {
    supported_backends(wgpu::Backends::from_env())
}

fn supported_backends(requested: Option<wgpu::Backends>) -> Result<wgpu::Backends> {
    let enabled = requested.unwrap_or(wgpu::Backends::VULKAN) & wgpu::Backends::VULKAN;
    if enabled.is_empty() {
        return Err(GpuError::NoAdapter(
            "native GPU support requires Vulkan: unset WGPU_BACKEND or set WGPU_BACKEND=vulkan, and install a Vulkan driver".into(),
        ));
    }
    Ok(enabled)
}

/// The adapter `policy` picks on this machine and its tier; `Ok(None)` for the CPU rasterizer.
pub fn adapter(policy: Policy) -> Result<Option<(wgpu::Adapter, Tier)>> {
    if policy == Policy::Cpu {
        return Ok(None);
    }
    let mut adapters = enumerate(backends()?);
    let types = device_types(&adapters);
    let Some(i) = pick(&types, policy)? else { return Ok(None) };
    let tier = types.get(i).map_or(Tier::Hardware, |t| rank(*t).0);
    Ok((i < adapters.len()).then(|| (adapters.swap_remove(i), tier)))
}

/// The Vulkan adapters. Never initialize an unsupported backend as a fallback.
fn enumerate(backends: wgpu::Backends) -> Vec<wgpu::Adapter> {
    pollster::block_on(instance().enumerate_adapters(backends))
}

static INSTANCE: OnceLock<wgpu::Instance> = OnceLock::new();

fn instance() -> &'static wgpu::Instance {
    INSTANCE.get_or_init(|| {
        wgpu::Instance::new(wgpu::InstanceDescriptor { backends: wgpu::Backends::VULKAN, ..wgpu::InstanceDescriptor::new_without_display_handle() })
    })
}

fn device_types(adapters: &[wgpu::Adapter]) -> Vec<wgpu::DeviceType> {
    adapters.iter().map(|a| a.get_info().device_type).collect()
}

#[cfg(test)]
mod tests {
    use wgpu::DeviceType::*;

    use super::*;

    #[test]
    fn only_vulkan_is_compiled_as_a_native_backend() {
        assert_eq!(wgpu::Instance::enabled_backend_features(), wgpu::Backends::VULKAN);
    }

    #[test]
    fn backend_overrides_can_never_enable_other_native_backends() {
        use wgpu::Backends as B;
        for requested in [None, Some(B::VULKAN), Some(B::all()), Some(B::DX12 | B::VULKAN)] {
            assert_eq!(supported_backends(requested).unwrap(), B::VULKAN);
        }
        for requested in [B::empty(), B::DX12, B::METAL, B::GL, B::BROWSER_WEBGPU, B::NOOP] {
            let error = supported_backends(Some(requested)).unwrap_err().to_string();
            assert!(error.contains("Vulkan"));
            assert!(error.contains("WGPU_BACKEND=vulkan"));
        }
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

    /// #636: repeated requests reuse the same instance without a new driver context.
    #[test]
    fn vulkan_keeps_one_instance() {
        assert!(std::ptr::eq(instance(), instance()));
    }

    #[test]
    fn policy_names() {
        assert_eq!(Policy::parse("Software"), Some(Policy::Software));
        assert_eq!(Policy::parse(""), Some(Policy::Auto));
        assert_eq!(Policy::parse("cpu"), Some(Policy::Cpu));
        assert_eq!(Policy::parse("gpu"), None);
    }
}
