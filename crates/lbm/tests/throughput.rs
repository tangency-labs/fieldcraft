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
