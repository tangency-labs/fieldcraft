//! The GPU device a test holds while it runs: one test at a time on the device (`turn`), a flow on
//! the adapter `FIELDCRAFT_GPU` picks, shared by the verification suite and the channel study.

use std::ops::{Deref, DerefMut};
use std::sync::{Mutex, MutexGuard};

use fieldcraft_gpu::select::{Policy, gpu_required};
use fieldcraft_lbm::{Error, Flow};
use fieldcraft_lbm_gpu::{GpuFlow, Options, Streaming};

/// Whose turn it is on the GPU.
static TURN: Mutex<()> = Mutex::new(());

/// The next turn on the GPU (a test that failed while holding it leaves it to the next).
pub fn turn() -> MutexGuard<'static, ()> {
    TURN.lock().unwrap_or_else(|e| e.into_inner())
}

/// A flow on the GPU, holding the turn until it is dropped (the flow first, then the turn).
pub struct Gpu {
    flow: GpuFlow,
    _turn: MutexGuard<'static, ()>,
}

impl Deref for Gpu {
    type Target = GpuFlow;
    fn deref(&self) -> &GpuFlow {
        &self.flow
    }
}

impl DerefMut for Gpu {
    fn deref_mut(&mut self) -> &mut GpuFlow {
        &mut self.flow
    }
}

/// The flow on the GPU `FIELDCRAFT_GPU` picks, or `None` to skip (a failure when a GPU is
/// required).
pub fn gpu(flow: &Flow, streaming: Streaming) -> Option<Gpu> {
    let policy = Policy::from_env().expect("FIELDCRAFT_GPU names a policy");
    let turn = turn();
    match GpuFlow::with_options(flow, policy, Options { streaming, profile: false }) {
        Ok(g) => {
            println!("GPU: {} ({streaming:?})", g.adapter_name());
            Some(Gpu { flow: g, _turn: turn })
        }
        Err(e) if gpu_required() => panic!("a GPU is required: {e}"),
        // Only a missing device skips: a shader that does not compile is a failure anywhere.
        Err(e @ Error::NoDevice) => {
            println!("skipping: {e}");
            None
        }
        Err(e) => panic!("{e}"),
    }
}
