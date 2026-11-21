//! Han et al.'s on-node wall (`spalding_bounce`) at an obstacle's walls runs on the CPU only: the
//! device's `node_bounce` is eq. 14 on a plane box face, so its setup refuses, naming the wall.
#[path = "../../lbm/tests/common/mod.rs"]
mod common;

use common::*;
use fieldcraft_gpu::select::{Policy, gpu_required};
use fieldcraft_lbm_gpu::{GpuFlow, Options};
use fieldcraft_model::{Collision, FaceCondition, Procedure, Turbulence, WallModel};

#[test]
fn the_on_node_wall_at_an_obstacle_is_refused_on_the_device() {
    let dx = 1.0 / 16.0;
    let mut d = domain([0.0; 3], [24.0 * dx, 16.0 * dx, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    let obstacles = vec![obstacle("Plate", cuboid([8.0 * dx, 6.0 * dx, -1.0], [16.0 * dx, 7.0 * dx, 1.0]), &CUBOID_FACES)];
    d.obstacles = obstacles.iter().map(|o| o.body.clone()).collect();
    d.acceleration = [1.0, 0.0, 0.0];
    let dt = 0.05 * dx;
    let mut p = collide(10.0 * dt, dt, None, Collision::Trt);
    if let Procedure::Fluid { turbulence, .. } = &mut p {
        *turbulence = Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: Some(WallModel::SpaldingBounce) });
    }
    let (_setup, flow) = prepare(&model(1.0, 0.02, d, p), &obstacles);
    assert!(flow.wall_nodes().is_some_and(|n| n.nodes().iter().any(|w| !w.flat)), "the plate makes nodes that are not flat");
    match GpuFlow::with_options(&flow, Policy::from_env().unwrap(), Options::default()) {
        Ok(_) => panic!("the device accepted an on-node wall at an obstacle"),
        Err(e) if !gpu_required() && e.to_string().contains("no GPU device") => println!("skipping: {e}"),
        Err(e) => {
            let text = e.to_string();
            assert!(text.contains("not supported yet") && text.contains("obstacle wall, an edge or a corner"), "{text}");
        }
    }
}
