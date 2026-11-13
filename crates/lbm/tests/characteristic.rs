//! The characteristic (partially non-reflecting) pressure outlet (`Flow::characteristic_outlet`):
//! it keeps a uniform stream and lets the box's acoustic mode out where the uniform outlet,
//! holding its density, reflects it.
//!
//! Each test prints its numbers (`cargo test -p fieldcraft-lbm --test characteristic -- --nocapture`).

mod common;

use common::*;
use fieldcraft_lbm::Flow;
use fieldcraft_model::{Axis, BoxFace, Collision, FaceCondition, FluidDomain, Model, Outflow, PorousZone, Profile, Resistance};

/// A 64 × 4 × 4 cell box (uniform inlet at x_min, periodic sides, TRT
/// τ 0.55, lattice speed 0.05) with a uniform pressure outlet at x_max and a density perturbation
/// `amplitude · cos(π x / L)`.
fn box_flow(amplitude: f64) -> Flow {
    flow_of(box_domain(), amplitude)
}

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

fn box_model(d: FluidDomain) -> Model {
    model(1.0, 0.05 / 3.0, d, collide(4000.0, 1.0, None, Collision::Trt))
}

fn flow_of(d: FluidDomain, amplitude: f64) -> Flow {
    let (_, mut flow) = prepare(&box_model(d), &[]);
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

/// `outflow: characteristic` on a model's pressure face sets the outlet up as the call does,
/// σ 0.25 when it gives none; a second such face is refused by name.
#[test]
fn a_model_s_characteristic_outflow_is_the_call_s_outlet() {
    let mut called = box_flow(1e-3);
    called.characteristic_outlet(BoxFace::XMax, 0.25).unwrap();
    let mut d = box_domain();
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Characteristic, sigma: None };
    let mut set = flow_of(d.clone(), 1e-3);
    assert_eq!(set.characteristic_outlet_cells(), 16);
    for _ in 0..200 {
        called.step().unwrap();
        set.step().unwrap();
    }
    for cell in 0..called.grid().len() {
        assert_eq!(called.moments(cell), set.moments(cell), "{cell}");
    }
    d.faces.x_min = d.faces.x_max;
    let e = fieldcraft_lbm::Setup::new(&box_model(d), None, &[], &fieldcraft_lbm::Limits::default()).unwrap_err().to_string();
    assert!(e.contains("more than one box face"), "{e}");
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

/// A stream through a porous zone whose last layer is the outlet's inner cell, Darcy's
/// resistance balanced by the body acceleration. Settled, the characteristic outlet holds the
/// face's density and the steady balance `L₅ = 0`: `u_b = u₁ + c_s (ρ₁ − ρ₀) / ρ₀` with u₁ the inner
/// cell's physical velocity, the zone's action included, as the closure reads it. (The uniform
/// outlet settles elsewhere: the open layers, outside the zone, leave ρ₁ off ρ₀ by 5e-5.)
#[test]
fn the_characteristic_outlet_reads_the_inner_velocity_of_a_porous_zone() {
    let mut d = domain([0.0; 3], [32.0, 2.0, 2.0], 1.0);
    for f in [&mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
    let permeability = 10.0;
    d.porous = vec![PorousZone {
        name: "Bed".into(),
        min: [1.0, 0.0, 0.0],
        max: [31.0, 2.0, 2.0],
        axis: Axis::X,
        resistance: Resistance::DarcyForchheimer { permeability, inertial: 0.0 },
    }];
    // ν u / K, the resistance per mass at the stream's speed.
    d.acceleration = [0.05 / 3.0 * 0.05 / permeability, 0.0, 0.0];
    let mut flow = flow_of(d, 0.0);
    flow.characteristic_outlet(BoxFace::XMax, 0.25).unwrap();
    for _ in 0..4000 {
        flow.step().unwrap();
    }
    let ((rho_b, u_b), (rho1, u1)) = (flow.moments(31).unwrap(), flow.moments(30).unwrap());
    let balance = u1[0] + (1.0_f64 / 3.0).sqrt() * (rho1 - 1.0) - u_b[0];
    println!("outlet ρ_b − ρ₀ {:.3e}, inner ρ₁ − ρ₀ {:.3e}, u_b − u₁ − c_s (ρ₁ − ρ₀) {balance:.3e}", rho_b - 1.0, rho1 - 1.0);
    assert!((rho_b - 1.0).abs() < 1e-9 && balance.abs() < 1e-9, "{rho_b} {balance}");
    assert_eq!(flow.characteristic_held(), 0);
}

/// The memory estimate counts the characteristic outlet's cells as the flow holds them.
#[test]
fn the_memory_estimate_counts_the_characteristic_outlet() {
    use fieldcraft_lbm::{Choice, Memory};
    // Include a non-power-of-two face, so spare vector capacity cannot escape the estimate.
    for [ny, nz] in [[4_u64, 4], [3, 3]] {
        let mut d = box_domain();
        d.max = [64.0, ny as f64, nz as f64];
        let plain = box_model(d.clone());
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Characteristic, sigma: None };
        let m = box_model(d);
        let (without, with) = (Memory::estimate(&plain, None, Choice::Cpu).unwrap(), Memory::estimate(&m, None, Choice::Cpu).unwrap());
        // Each outlet cell keeps its open-cell index, density and velocity (24 bytes).
        let extra = ny * nz * 24;
        assert_eq!((with.setup - without.setup, with.state - without.state), (extra, extra));
        let (setup, flow) = prepare(&m, &[]);
        assert_eq!(Memory::of(&setup, Choice::Cpu).unwrap(), with);
        assert_eq!(flow.heap_bytes() - prepare(&plain, &[]).1.heap_bytes(), extra);
    }
}
