//! CPU/device agreement on the independently constructed curved conjugate benchmark.
#[path = "../../lbm/tests/common/mod.rs"]
mod common;

use fieldcraft_gpu::select::{Policy, gpu_required};
use fieldcraft_lbm::{Limits, Setup};
use fieldcraft_lbm_gpu::{GpuFlow, Options};

#[test]
fn the_passive_conducting_cylinder_matches_the_cpu_during_startup() {
    let (model, obstacles) = common::curved_cht::model_case(40, 1, 512, 5.0, 0.025);
    let setup = Setup::new(&model, None, &obstacles, &Limits::default()).unwrap();
    let mut cpu = fieldcraft_lbm::Flow::new(&setup, &obstacles).unwrap();
    let mut gpu = match GpuFlow::with_options(&cpu, Policy::from_env().unwrap(), Options::default()) {
        Ok(gpu) => gpu,
        Err(e) if !gpu_required() && e.to_string().contains("no GPU device") => {
            println!("skipping: {e}");
            return;
        }
        Err(e) => panic!("{e}"),
    };
    println!("curved conjugate startup device: {}", gpu.adapter_name());
    let a = fieldcraft_lbm::run(&setup, &mut cpu, &Limits::default()).unwrap();
    let b = fieldcraft_lbm::run(&setup, &mut gpu, &Limits::default()).unwrap();
    let (a_field, b_field) = (&a.frames.last().unwrap().temperature, &b.frames.last().unwrap().temperature);
    assert_eq!(a_field.len(), b_field.len());
    assert!(a_field.iter().chain(b_field).all(|t| t.is_finite()), "both temperature fields are finite before comparing them");
    let difference = a_field.iter().zip(b_field).map(|(a, b)| (a - b).abs()).fold(0.0_f32, f32::max);
    assert!(difference < 2e-5, "same 1 K field scale and backend tolerance as the existing thermal verification: {difference}");
    let hot = |r: &fieldcraft_lbm::FluidResult| {
        *r.heat.iter().find(|h| h.body.as_deref() == Some("RightWall") && h.name == "x_min").unwrap().power.last().unwrap()
    };
    assert!((hot(&a) - hot(&b)).abs() < 1e-4 * hot(&a).abs(), "wall heat uses the existing thermal backend tolerance");
}
