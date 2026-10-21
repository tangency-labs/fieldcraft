//! The machine a session runs on: its physical memory, which bounds how far a fluid solve's
//! memory caps may be raised above their defaults ([`Limits::memory`](crate::Limits::memory)).

/// The share of physical memory a run with raised caps may take: half, leaving the rest to the
/// system, the application, its results and other programs.
const MEMORY_SHARE: u64 = 2;

/// The machine's physical memory in bytes (on Linux the root cgroup's limit when lower, so a
/// container's), or `None` where it cannot be read (the browser). Read once per process.
#[cfg(not(target_arch = "wasm32"))]
pub fn physical_memory() -> Option<u64> {
    use sysinfo::{MemoryRefreshKind, RefreshKind, System};
    static MEMORY: std::sync::OnceLock<Option<u64>> = std::sync::OnceLock::new();
    *MEMORY.get_or_init(|| {
        if !sysinfo::IS_SUPPORTED_SYSTEM {
            return None;
        }
        let system = System::new_with_specifics(RefreshKind::nothing().with_memory(MemoryRefreshKind::nothing().with_ram()));
        let total = system.total_memory();
        let total = system.cgroup_limits().map_or(total, |c| if c.total_memory > 0 { total.min(c.total_memory) } else { total });
        (total > 0).then_some(total)
    })
}

/// In the browser the memory is unknown: the caps keep their defaults.
#[cfg(target_arch = "wasm32")]
pub fn physical_memory() -> Option<u64> {
    None
}

/// The bound on a run with raised caps: [`MEMORY_SHARE`] of the physical memory, or 0 where it
/// is unknown (the caps then cannot be raised).
pub fn memory_bound() -> u64 {
    physical_memory().map_or(0, |m| m / MEMORY_SHARE)
}
