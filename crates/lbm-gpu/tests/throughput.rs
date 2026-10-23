//! Throughput of the GPU backend (ignored; run with
//! `cargo test --release -p fieldcraft-lbm-gpu --test throughput -- --ignored --nocapture`),
//! on the cases of `fieldcraft-lbm`'s throughput test, with both collisions, streaming in place
//! and with two buffers (what `Streaming::Auto` picks on this device is marked). `FIELDCRAFT_LBM_GPU_PROFILE=1`
//! adds each kernel's GPU time (timestamp queries; the throughput is then not representative).

#[path = "../../lbm/tests/common/mod.rs"]
mod common;

use common::*;
use fieldcraft_gpu::select::Policy;
use fieldcraft_lbm::backend::Backend;
use fieldcraft_lbm_gpu::{GpuFlow, Options, Streaming};
use fieldcraft_model::{BoxFace, FaceCondition};

/// Real-device cost of zone momentum forcing and of the optional diagnostic readback. Fan
/// feedback itself stays on the device; only this explicit zones() call reads its sample.
#[test]
#[ignore = "a real-device throughput measurement"]
fn zone_throughput_and_diagnostic_cost() {
    use fieldcraft_model::{Axis, Fan, FanDirection, PorousZone, Resistance};
    for kind in ["none", "porous", "fan + porous"] {
        let mut d = domain([0.0; 3], [64.0; 3], 1.0);
        for face in BoxFace::ALL {
            *d.faces.get_mut(face) = FaceCondition::Periodic {};
        }
        d.acceleration = [1e-6, 0.0, 0.0];
        if kind != "none" {
            d.porous = vec![PorousZone {
                name: "Vent".into(),
                min: [30.0, 0.0, 0.0],
                max: [34.0, 64.0, 64.0],
                axis: Axis::X,
                resistance: Resistance::DarcyForchheimer { permeability: 5.0, inertial: 2.0 },
            }];
        }
        if kind == "fan + porous" {
            d.fans = vec![Fan {
                name: "Fan".into(),
                min: [5.0, 0.0, 0.0],
                max: [7.0, 64.0, 64.0],
                direction: FanDirection::XPlus,
                curve: vec![[0.0, 0.002], [100.0, 0.0]],
                relaxation: Some(0.01),
            }];
        }
        let m = model(1.0, 0.1, d, step(232.0, 1.0, None));
        let (_, flow) = prepare(&m, &[]);
        let policy = Policy::from_env().expect("FIELDCRAFT_GPU names a policy");
        let Ok(mut device) = GpuFlow::with_options(&flow, policy, Options { streaming: Streaming::InPlace, profile: false }) else {
            println!("skipping {kind}: no GPU or not enough memory");
            return;
        };
        device.advance(32).unwrap();
        let start = std::time::Instant::now();
        device.advance(200).unwrap();
        let seconds = start.elapsed().as_secs_f64();
        let start = std::time::Instant::now();
        for _ in 0..10 {
            let _ = device.zones().unwrap();
        }
        let sample_ms = start.elapsed().as_secs_f64() * 100.0;
        println!(
            "{kind}: {}: {} cells, {:.2} M updates/s, {:.4} s/200 steps; zone diagnostic {:.3} ms/query; {} bytes/query",
            device.adapter_name(),
            device.fluid_cells(),
            device.fluid_cells() as f64 * 200.0 / seconds / 1e6,
            seconds,
            sample_ms,
            flow.zones().len() * 16
        );
    }
}

#[test]
#[ignore]
fn cell_updates_per_second() {
    let cases = [("periodic box", 64usize, 0.0), ("periodic box", 128, 0.0), ("box with a sphere", 128, 32.0)];
    for ((name, l, sphere_radius), collision) in cases.into_iter().flat_map(|c| COLLISIONS.map(move |k| (c, k))) {
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
                collision: Some(collision),
                steady: None,
            },
        );
        let (_, flow) = prepare(&m, &obstacles);
        let policy = Policy::from_env().expect("FIELDCRAFT_GPU names a policy");
        let profile = std::env::var("FIELDCRAFT_LBM_GPU_PROFILE").is_ok_and(|v| v == "1");
        let Ok(auto) = GpuFlow::new(&flow, policy) else {
            println!("skipping {name}, {l}³ ({collision:?}): no GPU or not enough memory");
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
                "{name}, {l}³, {collision:?} on {}, {}{}: {:.1} M cell updates/s (pipelines {:.2} s)",
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

/// The process's peak resident memory so far (Linux `VmHWM`), bytes; 0 elsewhere.
fn peak_resident() -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| s.lines().find(|l| l.starts_with("VmHWM:")).and_then(|l| l.split_whitespace().nth(1)).and_then(|v| v.parse::<u64>().ok()))
        .map_or(0, |kb| kb * 1024)
}

/// Conjugate heat transfer, CPU against GPU: air past a powered conducting base with six
/// conducting fin plates in a channel (`side` cells across, twice as long), heat on. Steady
/// state does not matter here: a fixed number of steps of the full CHT step on each backend,
/// and the host's peak memory around building the device flow. Set `FIELDCRAFT_CHT_SIDE`
/// (default 64) for a bigger grid.
#[test]
#[ignore]
fn conjugate_cell_updates_per_second() {
    let side: usize = std::env::var("FIELDCRAFT_CHT_SIDE").ok().and_then(|v| v.parse().ok()).unwrap_or(64);
    let (l, h) = (2.0 * side as f64, side as f64);
    let mut d = domain([0.0; 3], [l, h, h], 1.0);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: fieldcraft_model::Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: fieldcraft_model::Outflow::Convective };
    let mut solids = vec![conducting("Base", 1.0)];
    let mut obstacles = vec![obstacle("Base", cuboid([0.3 * l, 0.25 * h, 0.0], [0.6 * l, 0.75 * h, 0.1 * h]), &CUBOID_FACES)];
    for k in 0..6 {
        let y = 0.28 * h + k as f64 * 0.08 * h;
        let name = format!("Fin{k}");
        obstacles.push(obstacle(&name, cuboid([0.3 * l, y, 0.1 * h], [0.6 * l, y + 0.03 * h, 0.6 * h]), &CUBOID_FACES));
        solids.push(conducting(&name, 0.0));
    }
    d.obstacles = obstacles.iter().map(|o| o.body.clone()).collect();
    d.thermal = Some(fieldcraft_model::FluidThermal { initial_temperature: kelvin(20.0), solids, ..Default::default() });
    let steps = 100;
    let procedure = fieldcraft_model::Procedure::Fluid {
        domain: "Flow".into(),
        end_time: None,
        steps: Some(steps),
        time_step: Some(1.0),
        frames: 1,
        magic: None,
        collision: Some(fieldcraft_model::Collision::Regularized),
        steady: None,
    };
    let m = cht_model([1.0, 0.02, 0.03], (2.0, [2.0; 3], 1.0), d, procedure);
    let (_, mut cpu) = prepare(&m, &obstacles);
    let cells = cpu.grid().len();
    let conducting: usize = Backend::conducting(&cpu).iter().map(Vec::len).sum();
    let before = peak_resident();
    let policy = Policy::from_env().expect("FIELDCRAFT_GPU names a policy");
    let Ok(mut gpu) = GpuFlow::new(&cpu, policy) else {
        println!("skipping the conjugate case: no GPU or not enough memory");
        return;
    };
    let after = peak_resident();
    gpu.advance(1).unwrap();
    let start = std::time::Instant::now();
    gpu.advance(steps).unwrap();
    let _ = gpu.heat().unwrap();
    let g = start.elapsed().as_secs_f64();
    let cpu_steps = (steps / 10).max(5);
    let start = std::time::Instant::now();
    cpu.advance(cpu_steps).unwrap();
    let c = start.elapsed().as_secs_f64();
    let rate = |steps: u64, s: f64| cells as f64 * steps as f64 / s / 1e6;
    println!(
        "conjugate, {cells} cells ({conducting} conducting) on {}: GPU {:.1} M cell updates/s, CPU {:.1} ({} threads): {:.1}x; host peak {:.0} B/cell before the device flow, +{:.0} B/cell building it",
        gpu.adapter_name(),
        rate(steps, g),
        rate(cpu_steps, c),
        std::thread::available_parallelism().map_or(1, |n| n.get()),
        rate(steps, g) / rate(cpu_steps, c),
        before as f64 / cells as f64,
        after.saturating_sub(before) as f64 / cells as f64,
    );
}

/// The fixed cost of a GPU solve: three solves in a row of one step of a small heated channel
/// with a conducting block (every kernel), each building its device flow from the start as
/// `solve.fluid` does, timed from the build to the wall heat read back (#588).
#[test]
#[ignore]
fn fixed_cost_per_solve() {
    let (l, h) = (32.0, 16.0);
    let mut d = domain([0.0; 3], [l, h, h], 1.0);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: fieldcraft_model::Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: fieldcraft_model::Outflow::Convective };
    let obstacles = vec![obstacle("Block", cuboid([8.0, 4.0, 0.0], [16.0, 12.0, 4.0]), &CUBOID_FACES)];
    d.obstacles = vec!["Block".into()];
    d.thermal =
        Some(fieldcraft_model::FluidThermal { initial_temperature: kelvin(20.0), solids: vec![conducting("Block", 1.0)], ..Default::default() });
    let procedure = fieldcraft_model::Procedure::Fluid {
        domain: "Flow".into(),
        end_time: None,
        steps: Some(1),
        time_step: Some(1.0),
        frames: 1,
        magic: None,
        collision: None,
        steady: None,
    };
    let m = cht_model([1.0, 0.02, 0.03], (2.0, [2.0; 3], 1.0), d, procedure);
    let (_, flow) = prepare(&m, &obstacles);
    let policy = Policy::from_env().expect("FIELDCRAFT_GPU names a policy");
    for solve in 1..=3 {
        let start = std::time::Instant::now();
        let Ok(mut gpu) = GpuFlow::new(&flow, policy) else {
            println!("skipping: no GPU");
            return;
        };
        gpu.advance(1).unwrap();
        let _ = gpu.heat().unwrap();
        println!(
            "solve {solve} on {}: {:.3} s from the build to the first step's heat (pipelines {:.3} s)",
            gpu.adapter_name(),
            start.elapsed().as_secs_f64(),
            gpu.pipeline_seconds()
        );
    }
}
