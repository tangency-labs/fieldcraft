//! Han, Ooka and Kikumoto's on-node wall on the GPU against the CPU where the wall is not a plane
//! box face: a block on a channel floor, so the device runs the general node path (a node's links
//! to the fluid and its link pairs, edge and corner nodes, interior cells of several nodes) with
//! Smagorinsky's model. Velocity within 1e-4 of the largest CPU speed, density 1e-5, every group's
//! force averaged over the last 200 steps 1e-3 of the largest group force; the bounds of the other
//! LES cases here. Prints what it measured.

use super::*;
use fieldcraft_model::{FaceCondition, WallModel};

/// A block on a forced channel's floor: 32 × 16 × 16 cells at Δ 1/8, periodic in x and z, Han's
/// walls at the y faces, the block 8 × 5 × 4 cells, started at the bulk speed 15.7 (lattice 0.1),
/// τ − ½ = `tau_half`, Smagorinsky's model with Han's wall.
fn floor_block(collision: Collision, tau_half: f64, steps: f64) -> (Setup, Flow) {
    let (dx, ub) = (1.0 / 8.0, 15.7);
    let mut d = domain([0.0; 3], [4.0, 2.0, 2.0], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    let obstacles = vec![obstacle("Block", cuboid([1.5, -1.0, 0.75], [2.5, 5.0 * dx, 1.25]), &CUBOID_FACES)];
    d.obstacles = vec!["Block".into()];
    d.acceleration = [1.0, 0.0, 0.0];
    d.initial_velocity = [ub, 0.0, 0.0];
    let dt = 0.1 * dx / ub;
    let nu = tau_half / 3.0 * dx * dx / dt;
    let mut p = collide(steps * dt, dt, None, collision);
    if let Procedure::Fluid { turbulence, .. } = &mut p {
        *turbulence = Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: Some(WallModel::SpaldingBounce) });
    }
    prepare(&model(1.0, nu, d, p), &obstacles)
}

/// The block on the floor, GPU against CPU, both layouts, the two-relaxation-time and the
/// multiple-relaxation-time collisions at τ − ½ 0.1, 400 steps.
#[test]
fn han_s_wall_follows_the_cpu_around_a_block() {
    let mut failed = Vec::new();
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        for collision in [Collision::Trt, Collision::Mrt] {
            let (setup, mut cpu) = floor_block(collision, 0.1, 400.0);
            let nodes = cpu.wall_nodes().map(|w| w.nodes().len()).unwrap_or(0);
            assert!(nodes > 1000, "the block's walls are boundary nodes: {nodes}");
            let Some(mut gpu) = gpu(&cpu, streaming) else { return };
            let tail = 200u64;
            cpu.advance(setup.steps - tail).unwrap();
            gpu.advance(setup.steps - tail).unwrap();
            let groups = cpu.groups().len();
            let (mut mc, mut mg) = (vec![[0.0f64; 3]; groups], vec![[0.0f64; 3]; groups]);
            for _ in 0..tail {
                cpu.advance(1).unwrap();
                gpu.advance(1).unwrap();
                let (fc, fg) = (Backend::forces(&mut cpu).unwrap(), gpu.forces().unwrap());
                assert_eq!(fg.len(), groups);
                for g in 0..groups {
                    for k in 0..3 {
                        mc[g][k] += fc[g][k] / tail as f64;
                        mg[g][k] += fg[g][k] / tail as f64;
                    }
                }
            }
            let norm = |v: &[f64; 3]| v.iter().map(|x| x * x).sum::<f64>().sqrt();
            let scale = mc.iter().map(norm).fold(0.0, f64::max).max(1e-300);
            let dforce = mc.iter().zip(&mg).map(|(a, b)| norm(&[a[0] - b[0], a[1] - b[1], a[2] - b[2]])).fold(0.0, f64::max) / scale;
            let (du, drho) = field_difference(&mut cpu, &mut gpu);
            let what = format!("Han's wall, block on a floor, {collision:?}, {streaming:?}, τ − ½ 0.1, 400 steps");
            println!("{what}: {nodes} nodes; velocity {du:.2e}, density {drho:.2e}, group forces, mean of 200 steps {dforce:.2e}");
            if !(du < 1e-4 && drho < 1e-5 && dforce < 1e-3) {
                failed.push(what);
            }
        }
    }
    assert!(failed.is_empty(), "{failed:?}");
}
