//! The characteristic (partially non-reflecting) pressure outlet (`Flow::characteristic_outlet`):
//! it keeps a uniform stream and lets the box's acoustic mode out where the uniform outlet,
//! holding its density, reflects it.
//!
//! Each test prints its numbers (`cargo test -p fieldcraft-lbm --test characteristic -- --nocapture`).

mod common;

use common::*;
use fieldcraft_lbm::Flow;
use fieldcraft_model::{BoxFace, Collision, FaceCondition, Outflow, Profile};

/// A 64 × 4 × 4 cell box (uniform inlet at x_min, periodic sides, TRT
/// τ 0.55, lattice speed 0.05) with a uniform pressure outlet at x_max and a density perturbation
/// `amplitude · cos(π x / L)`.
fn box_flow(amplitude: f64) -> Flow {
    let mut d = domain([0.0; 3], [64.0, 4.0, 4.0], 1.0);
    for f in [&mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.initial_velocity = [0.05, 0.0, 0.0];
    let (_, mut flow) = prepare(&model(1.0, 0.05 / 3.0, d, collide(4000.0, 1.0, None, Collision::Trt)), &[]);
    let grid = *flow.grid();
    flow.set_equilibrium(|cell| (1.0 + amplitude * (std::f64::consts::PI * grid.centre(cell)[0] / 64.0).cos(), [0.05, 0.0, 0.0]));
    flow
}

fn stats(flow: &Flow) -> (f64, f64) {
    let (mut sum, mut sq, mut n) = (0.0, 0.0, 0.0);
    for cell in 0..flow.grid().len() {
        if let Some((rho, _)) = flow.moments(cell) {
            sum += rho - 1.0;
            sq += (rho - 1.0) * (rho - 1.0);
            n += 1.0;
        }
    }
    ((sq / n).sqrt(), sum / n)
}

#[test]
fn the_characteristic_outlet_keeps_a_uniform_stream() {
    let mut plain = box_flow(0.0);
    let mut open = box_flow(0.0);
    assert_eq!(open.characteristic_outlet(BoxFace::XMax, 0.25).unwrap(), 16);
    for _ in 0..500 {
        plain.step().unwrap();
        open.step().unwrap();
    }
    for cell in 0..plain.grid().len() {
        let ((r0, u0), (r1, u1)) = (plain.moments(cell).unwrap(), open.moments(cell).unwrap());
        assert!((r0 - r1).abs() < 1e-13 && (0..3).all(|a| (u0[a] - u1[a]).abs() < 1e-13), "{cell}");
    }
    assert_eq!(open.characteristic_held(), 0);
    let mut f = box_flow(0.0);
    assert!(f.characteristic_outlet(BoxFace::XMin, 0.25).is_err());
    assert!(f.characteristic_outlet(BoxFace::XMax, 0.0).is_err());
    // One face: a second call is refused and keeps the first.
    assert_eq!(f.characteristic_outlet(BoxFace::XMax, 0.25).unwrap(), 16);
    assert!(f.characteristic_outlet(BoxFace::XMax, 0.5).unwrap_err().to_string().contains("more than one face"));
    assert_eq!(f.characteristic_outlet_cells(), 16);
    f.step().unwrap();
    assert!(f.characteristic_outlet(BoxFace::XMax, 0.25).unwrap_err().to_string().contains("before the first step"));
}

#[test]
fn the_characteristic_outlet_lets_the_box_s_acoustic_mode_out() {
    let run = |characteristic: bool| {
        let mut flow = box_flow(1e-3);
        if characteristic {
            flow.characteristic_outlet(BoxFace::XMax, 0.25).unwrap();
        }
        let mut history = Vec::new();
        for k in 0..2000 {
            flow.step().unwrap();
            if k % 250 == 249 {
                history.push(stats(&flow));
            }
        }
        history
    };
    let (plain, open) = (run(false), run(true));
    for (k, (a, b)) in plain.iter().zip(&open).enumerate() {
        println!(
            "step {}: density RMS {:.3e} (mean {:+.2e}) uniform outlet, {:.3e} (mean {:+.2e}) characteristic",
            250 * (k + 1),
            a.0,
            a.1,
            b.0,
            b.1
        );
    }
    // The characteristic outlet lets the wave out: what is left after 2000 steps is below a
    // hundredth of what the reflecting uniform outlet keeps (the oscillating part and the mean
    // alike), and the uniform outlet's mean does oscillate.
    let (a, b) = (plain.last().unwrap(), open.last().unwrap());
    assert!(b.0 < 0.01 * a.0, "{} vs {}", b.0, a.0);
    assert!(plain.windows(2).any(|w| w[0].1.signum() != w[1].1.signum()), "the uniform outlet's mean oscillates");
}
