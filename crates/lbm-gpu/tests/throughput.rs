//! Throughput of the GPU backend (ignored; run with
//! `cargo test --release -p fieldcraft-lbm-gpu --test throughput -- --ignored --nocapture`),
//! on the cases of `fieldcraft-lbm`'s throughput test, streaming in place and with two buffers
//! (what `Streaming::Auto` picks on this device is marked). `FIELDCRAFT_LBM_GPU_PROFILE=1`
//! adds each kernel's GPU time (timestamp queries; the throughput is then not representative).

#[path = "../../lbm/tests/common/mod.rs"]
mod common;

use common::*;
use fieldcraft_gpu::select::Policy;
use fieldcraft_lbm::backend::Backend;
use fieldcraft_lbm_gpu::{GpuFlow, Options, Streaming};
use fieldcraft_model::{BoxFace, FaceCondition};

#[test]
#[ignore]
fn cell_updates_per_second() {
    for (name, l, sphere_radius) in [("periodic box", 64usize, 0.0), ("periodic box", 128, 0.0), ("box with a sphere", 128, 32.0)] {
        let side = l as f64;
        let mut d = domain([0.0; 3], [side; 3], 1.0);
        for f in BoxFace::ALL {
            *d.faces.get_mut(f) = FaceCondition::Periodic {};
        }
        d.acceleration = [1e-6, 0.0, 0.0];
        let mut obstacles = Vec::new();
        if sphere_radius > 0.0 {
            d.obstacles = vec!["Ball".into()];
            obstacles.push(obstacle("Ball", sphere([side / 2.0; 3], sphere_radius, 4), &["skin"]));
        }
        let steps = 200;
        let m = model(
            1.0,
            1.0 / 6.0,
            d,
            fieldcraft_model::Procedure::Fluid {
                domain: "Flow".into(),
                end_time: None,
                steps: Some(steps),
                time_step: Some(1.0),
                frames: 1,
                magic: None,
                steady: None,
            },
        );
        let (_, flow) = prepare(&m, &obstacles);
        let policy = Policy::from_env().expect("FIELDCRAFT_GPU names a policy");
        let profile = std::env::var("FIELDCRAFT_LBM_GPU_PROFILE").is_ok_and(|v| v == "1");
        let Ok(auto) = GpuFlow::new(&flow, policy) else {
            println!("skipping {name}, {l}³: no GPU or not enough memory");
            continue;
        };
        let picked = auto.streams_in_place();
        drop(auto);
        for (streaming, in_place) in [(Streaming::InPlace, true), (Streaming::TwoBuffers, false)] {
            let mut gpu = GpuFlow::with_options(&flow, policy, Options { streaming, profile }).unwrap();
            // One warm-up chunk, then the timed run: the chunk's steps are one submission.
            gpu.advance(1).unwrap();
            let start = std::time::Instant::now();
            gpu.advance(steps).unwrap();
            let _ = gpu.mass().unwrap();
            let s = start.elapsed().as_secs_f64();
            println!(
                "{name}, {l}³ on {}, {}{}: {:.1} M cell updates/s (pipelines {:.2} s)",
                gpu.adapter_name(),
                if in_place { "in place" } else { "two buffers" },
                if in_place == picked { " (auto)" } else { "" },
                (gpu.fluid_cells() as f64 * steps as f64) / s / 1e6,
                gpu.pipeline_seconds()
            );
            for (kernel, t) in gpu.kernel_times().unwrap_or_default() {
                println!("    {kernel}: {} dispatches, {:.3} ms each", t.dispatches, 1e3 * t.seconds / t.dispatches.max(1) as f64);
            }
        }
    }
}
