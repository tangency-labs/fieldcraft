//! The absorbing layers (`crate::absorbing`) in whole runs: a uniform stream through them keeps
//! its state, a plane pulse returns the stated reflection to the grid's and the viscosity's bounds,
//! the stream lowers the return as the round trip says, and the layers are counted, warned about
//! and refused where the solver does not run them.
//!
//! Each test prints its numbers (`cargo test -p fieldcraft-lbm --test absorbing -- --nocapture`).

mod common;

use common::*;
use fieldcraft_lbm::reflection::layer_reflection;
use fieldcraft_lbm::{Choice, Flow, Memory};
use fieldcraft_model::{AbsorbingLayer, BoxFace, Collision, FaceCondition, FluidDomain, Outflow, Profile};

const COLLISIONS: [Collision; 4] = [Collision::Trt, Collision::Regularized, Collision::Hrr, Collision::Mrt];

/// The 64 × 4 × 4 box of the characteristic outlet's tests: a uniform inlet at x_min, a uniform pressure
/// outlet at x_max, periodic sides, lattice speed 0.05 and ν = 0.05/3.
fn box_domain() -> FluidDomain {
    let mut d = domain([0.0; 3], [64.0, 4.0, 4.0], 1.0);
    for f in [&mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
    d.initial_velocity = [0.05, 0.0, 0.0];
    d
}

fn box_model(d: FluidDomain, collision: Collision) -> fieldcraft_model::Model {
    model(1.0, 0.05 / 3.0, d, collide(4000.0, 1.0, None, collision))
}

fn layer(face: BoxFace, cells: f64, reflection: f64, reference: Option<[f64; 3]>) -> AbsorbingLayer {
    AbsorbingLayer { face, width: cells, reflection, reference_velocity: reference }
}

/// The box from the uniform stream with a density perturbation `amplitude · cos(π x / 64)`.
fn perturbed(d: FluidDomain, collision: Collision, amplitude: f64) -> (fieldcraft_lbm::Setup, Flow) {
    let (setup, mut flow) = prepare(&box_model(d, collision), &[]);
    let grid = *flow.grid();
    flow.set_equilibrium(|cell| (1.0 + amplitude * libm::cos(std::f64::consts::PI * grid.centre(cell)[0] / 64.0), [0.05, 0.0, 0.0]));
    (setup, flow)
}

/// Layers on both ends of the box keep a uniform stream as the box without them keeps it: the
/// layers pull toward the state they find, each collision.
#[test]
fn layers_keep_a_uniform_stream() {
    for collision in COLLISIONS {
        let mut plain = perturbed(box_domain(), collision, 0.0).1;
        let mut d = box_domain();
        d.layers = vec![
            layer(BoxFace::XMin, 8.0, 0.01, None),
            layer(BoxFace::XMax, 8.0, 0.01, Some([0.05, 0.0, 0.0])),
        ];
        let mut layered = perturbed(d, collision, 0.0).1;
        for _ in 0..500 {
            plain.step().unwrap();
            layered.step().unwrap();
        }
        for cell in 0..plain.grid().len() {
            let ((r0, u0), (r1, u1)) = (plain.moments(cell).unwrap(), layered.moments(cell).unwrap());
            assert!((r0 - r1).abs() < 1e-13 && (0..3).all(|a| (u0[a] - u1[a]).abs() < 1e-13), "{collision:?} cell {cell}");
        }
    }
}

/// The layers' pulse returns the stated reflection, up to the viscous bound the setup reports: at
/// two reflections, and the slope of the log of the measured return against the log of the stated one.
#[test]
fn a_plane_pulse_returns_the_stated_reflection() {
    for collision in COLLISIONS {
        let a = layer_reflection(collision, 40, 0.1, 0.0).unwrap();
        let b = layer_reflection(collision, 40, 0.01, 0.0).unwrap();
        println!(
            "{collision:?}: R 0.1 measured {:.6e} (plain {:.4}, viscous bound {:.2e}); R 0.01 measured {:.6e} (plain {:.4}, viscous bound {:.2e})",
            a.measured, a.plain, a.report.viscous_reflection, b.measured, b.plain, b.report.viscous_reflection
        );
        assert!((a.measured - 0.1).abs() <= a.report.viscous_reflection, "{collision:?}: {}", a.measured);
        assert!((b.measured - 0.01).abs() <= b.report.viscous_reflection, "{collision:?}: {}", b.measured);
        let slope = libm::log(b.measured) / libm::log(a.measured);
        let bound = a.report.viscous_reflection / (0.1 * std::f64::consts::LN_10) + b.report.viscous_reflection / (0.01 * libm::log(100.0));
        assert!((slope - 2.0).abs() <= 2.0 * bound, "{collision:?}: slope {slope}");
    }
    // The grid's own part of the reflection, for the manual: printed, not asserted.
    for cells in [8, 16] {
        let r = layer_reflection(Collision::Trt, cells, 0.01, 0.0).unwrap();
        println!("TRT, {cells} cells, R 0.01: measured {:.6e}", r.measured);
    }
}

/// A stream of 0.1 along the layer's normal returns the round trip's value, R^(1 / (1 − M²)), and not the value at rest.
#[test]
fn a_stream_lowers_the_return_as_the_round_trip_says() {
    let r = layer_reflection(Collision::Regularized, 40, 0.1, 0.1).unwrap();
    let predicted = libm::exp(libm::log(0.1) / (1.0 - 0.01 / (1.0 / 3.0)));
    println!("stream 0.1: measured {:.6e}, round trip {predicted:.6e}, at rest 0.1", r.measured);
    assert!(((r.report.stream_reflection - predicted) / predicted).abs() < 1e-12);
    assert!((r.measured - predicted).abs() < 0.5 * (0.1 - predicted), "{}", r.measured);
}

/// The strongest layer in the box runs: one cell at R = 1e-12 (rate 7.976 per step, above the speed of
/// sound) stays stable for 2000 steps, and the setup says it cannot resolve that reflection.
#[test]
fn the_strongest_layer_runs_and_is_warned() {
    for collision in [Collision::Trt, Collision::Regularized] {
        let mut d = box_domain();
        d.layers = vec![layer(BoxFace::XMax, 1.0, 1e-12, Some([0.05, 0.0, 0.0]))];
        let (setup, mut flow) = perturbed(d, collision, 1e-3);
        for _ in 0..2000 {
            flow.step().unwrap();
        }
        for cell in 0..flow.grid().len() {
            let (rho, _) = flow.moments(cell).unwrap();
            assert!((rho - 1.0).abs() < 2e-3, "{collision:?} cell {cell}: density {rho}");
        }
        assert!(setup.warnings.iter().any(|w| w.contains("do not resolve")), "{:?}", setup.warnings);
        let peak = setup.layers[0].report.peak_rate;
        let expected = libm::log(1e12) * (1.0_f64 / 3.0).sqrt() / 2.0;
        assert!(((peak - expected) / expected).abs() < 1e-12, "{peak} against {expected}");
    }
}

/// The memory estimate counts the layers' cells as the flow holds them: 48 bytes each (an index, a half
/// factor, a density and a velocity).
#[test]
fn the_memory_estimate_counts_the_layers_cells() {
    for [ny, nz] in [[4_u64, 4], [3, 3]] {
        let mut d = box_domain();
        d.max = [64.0, ny as f64, nz as f64];
        let plain = box_model(d.clone(), Collision::Trt);
        d.layers = vec![layer(BoxFace::XMax, 8.0, 0.01, Some([0.05, 0.0, 0.0]))];
        let m = box_model(d, Collision::Trt);
        let (without, with) = (Memory::estimate(&plain, None, Choice::Cpu).unwrap(), Memory::estimate(&m, None, Choice::Cpu).unwrap());
        let extra = 8 * ny * nz * 48;
        assert_eq!((with.setup - without.setup, with.state - without.state), (extra, extra));
        let (setup, flow) = prepare(&m, &[]);
        assert_eq!(Memory::of(&setup, Choice::Cpu).unwrap(), with);
        assert_eq!(flow.heap_bytes() - prepare(&plain, &[]).1.heap_bytes(), extra);
        assert_eq!(flow.absorbing_cells() as u64, 8 * ny * nz);
    }
}

/// Layers are refused on the GPU, and with a free surface, by name.
#[test]
fn layers_are_refused_on_the_gpu_and_with_a_free_surface() {
    let mut d = box_domain();
    d.layers = vec![layer(BoxFace::XMax, 8.0, 0.01, Some([0.05, 0.0, 0.0]))];
    let (setup, flow) = prepare(&box_model(d, Collision::Trt), &[]);
    let refused = fieldcraft_lbm::backend::backend(Choice::Device, flow, &setup).err().unwrap().to_string();
    assert!(refused.contains("absorbing layers"), "{refused}");

    let mut d = box_domain();
    d.free_surface = Some(fieldcraft_model::FreeSurface { liquid: vec![] });
    d.layers = vec![layer(BoxFace::XMin, 4.0, 0.01, None)];
    let e = fieldcraft_lbm::Setup::new(&box_model(d, Collision::Trt), None, &[], &fieldcraft_lbm::Limits::default()).unwrap_err().to_string();
    assert!(e.contains("not supported yet") && e.contains("absorbing layers"), "{e}");
}

/// An obstacle inside a layer is counted: its four cells of the layer's inner column are not fluid, so
/// they are not absorbing cells, and the flow counts them as blocked.
#[test]
fn an_obstacle_inside_a_layer_is_counted() {
    let mut d = box_domain();
    d.obstacles = vec!["Post".into()];
    d.layers = vec![layer(BoxFace::XMax, 4.0, 0.01, Some([0.05, 0.0, 0.0]))];
    let post = [obstacle("Post", cuboid([62.0, 0.0, 0.0], [63.0, 2.0, 2.0]), &CUBOID_FACES)];
    let (_, flow) = prepare(&box_model(d, Collision::Trt), &post);
    assert_eq!(flow.absorbing_blocked(), 4);
    assert_eq!(flow.absorbing_cells(), 4 * 16 - 4);
}
