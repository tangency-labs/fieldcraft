//! Throughput of the solver (ignored; run with
//! `cargo test --release -p fieldcraft-lbm --test throughput -- --ignored --nocapture`).

mod common;

use common::*;
use fieldcraft_model::{BoxFace, FaceCondition};

#[test]
#[ignore]
fn cell_updates_per_second() {
    let cases = [("periodic box", 64usize, 0.0), ("periodic box", 60, 0.0), ("box with a sphere", 64, 16.0)];
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
        let steps = 100;
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
                turbulence: None,
                steady: None,
            },
        );
        let (setup, mut flow) = prepare(&m, &obstacles);
        let start = std::time::Instant::now();
        run_all(&setup, &mut flow);
        let s = start.elapsed().as_secs_f64();
        println!("{name}, {l}³, {collision:?}: {:.1} M cell updates/s", (flow.fluid_cells() as f64 * steps as f64) / s / 1e6);
    }
}

/// How the solver scales with threads on the grid of the `flow-past` sample (96 × 48 × 40 cells:
/// a velocity inlet, a pressure outlet, slip sides and top, a floor, a bracket-sized block on
/// it), with each thread count's result checked bit for bit against one thread's. Set
/// `FIELDCRAFT_SCALING_STEPS` (default 200) for longer runs.
#[test]
#[ignore]
fn scaling_over_threads() {
    let steps: u64 = std::env::var("FIELDCRAFT_SCALING_STEPS").ok().and_then(|v| v.parse().ok()).unwrap_or(200);
    let mut d = domain([0.0; 3], [96.0, 48.0, 40.0], 1.0);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: fieldcraft_model::Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: fieldcraft_model::Outflow::Uniform };
    for f in [BoxFace::YMin, BoxFace::YMax, BoxFace::ZMax] {
        *d.faces.get_mut(f) = FaceCondition::Slip {};
    }
    d.initial_velocity = [0.05, 0.0, 0.0];
    d.obstacles = vec!["Bracket".into()];
    let obstacles = vec![obstacle("Bracket", cuboid([24.0, 16.0, 0.0], [28.0, 32.0, 20.0]), &CUBOID_FACES)];
    let m = model(
        1.0,
        0.02,
        d,
        fieldcraft_model::Procedure::Fluid {
            domain: "Flow".into(),
            end_time: None,
            steps: Some(steps),
            time_step: Some(1.0),
            frames: 1,
            magic: None,
            collision: None,
            steady: None,
            turbulence: None,
        },
    );
    let available = fieldcraft_parallel::available();
    let mut counts = vec![1, 2, 4, 8];
    counts.retain(|t| *t < available);
    counts.push(available);
    let mut reference: Option<(f64, Vec<u64>)> = None;
    for threads in counts {
        let (setup, mut flow) = prepare(&m, &obstacles);
        let start = std::time::Instant::now();
        fieldcraft_parallel::with_threads(threads, || run_all(&setup, &mut flow));
        let s = start.elapsed().as_secs_f64();
        let bits: Vec<u64> = (0..flow.grid().len()).filter_map(|c| flow.moments(c)).flat_map(|(rho, u)| [rho.to_bits(), u[0].to_bits()]).collect();
        let checksum = bits.iter().fold(0u64, |h, b| h.rotate_left(5) ^ b);
        let one = reference.get_or_insert((s, bits.clone()));
        assert!(one.1 == bits, "{threads} threads differ from one");
        println!(
            "{threads:>3} threads: {s:.2} s, {:.1} M cell updates/s, {:.2}x (state {checksum:016x})",
            (flow.fluid_cells() as f64 * steps as f64) / s / 1e6,
            one.0 / s
        );
    }
}
