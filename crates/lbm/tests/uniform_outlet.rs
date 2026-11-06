//! #819: an equilibrium stream is a fixed point even with an open pressure boundary.
//! Theory: constant density and velocity solve the unforced Navier–Stokes equations.
mod common;

use common::*;
use fieldcraft_lbm::{Flow, Limits, Setup};
use fieldcraft_model::{Collision, FluidDomain, Procedure, Turbulence};
use serde_json::json;

fn uniform(collision: Collision, tau: f64, size: [usize; 3], axis: usize, sign: f64, outlet: &str, steps: u64) -> Result<f64, String> {
    let mut dimensions = size;
    dimensions.swap(0, axis);
    let mut velocity = [0.0; 3];
    velocity[axis] = sign * 0.05;
    let periodic = json!({"type": "periodic"});
    let mut faces = json!({"x_min": periodic, "x_max": periodic, "y_min": periodic,
        "y_max": periodic, "z_min": periodic, "z_max": periodic});
    let ends = [["x_min", "x_max"], ["y_min", "y_max"], ["z_min", "z_max"]][axis];
    if outlet != "periodic" {
        let inlet = usize::from(sign < 0.0);
        faces[ends[inlet]] = json!({"type": "velocity", "velocity": velocity, "profile": "uniform"});
        faces[ends[1 - inlet]] = json!({"type": "pressure", "pressure": 0.0, "outflow": outlet});
    }
    let domain: FluidDomain = serde_json::from_value(json!({"material": "Fluid", "min": [0.0, 0.0, 0.0],
        "max": dimensions, "spacing": 1.0, "faces": faces, "initial_velocity": velocity}))
    .unwrap();
    let mut procedure = collide(steps as f64, 1.0, None, collision);
    if tau < 0.505
        && let Procedure::Fluid { turbulence, .. } = &mut procedure
    {
        // Exercise the supported turbulent relaxation floor, with no production-limit override.
        *turbulence = Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: None });
    }
    let setup = Setup::new(&model(1.0, (tau - 0.5) / 3.0, domain, procedure), None, &[], &Limits::default()).unwrap();
    assert!((setup.tau_plus - tau).abs() < 1e-14);
    let mut flow = Flow::new(&setup, &[]).unwrap();
    flow.set_equilibrium(|_| (1.0, velocity));
    let mut worst = 0.0f64;
    for step in 1..=steps {
        flow.step().map_err(|e| format!("step {step}: {e}"))?;
        for cell in 0..flow.grid().len() {
            let (rho, u) = flow.moments(cell).unwrap();
            worst = worst.max((rho - 1.0).abs());
            for a in 0..3 {
                worst = worst.max((u[a] - velocity[a]).abs());
            }
        }
    }
    Ok(worst)
}

#[test]
fn uniform_outlet_does_not_amplify_roundoff() {
    let worst = uniform(Collision::Trt, 0.518, [4, 1, 1], 0, 1.0, "uniform", 2000).unwrap();
    println!("TRT tau=.518 4x1x1 2000 steps: {worst:.4e}");
    assert!(worst <= 1e-12, "{worst}");
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn uniform_outlet_acceptance_sweep() {
    let mut failures = Vec::new();
    for collision in [Collision::Trt, Collision::Mrt, Collision::Regularized, Collision::Hrr] {
        for tau in [0.5001, 0.505, 0.518, 0.55] {
            for size in [[4, 1, 1], [16, 8, 2]] {
                for axis in 0..3 {
                    for sign in [-1.0, 1.0] {
                        let result = uniform(collision, tau, size, axis, sign, "uniform", 30_000);
                        let label = format!("{collision:?} tau={tau} size={size:?} axis={axis} sign={sign}");
                        println!("{label}: {result:?}");
                        if !matches!(result, Ok(worst) if worst <= 1e-12) {
                            failures.push((label, result));
                        }
                    }
                }
            }
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn convective_and_periodic_uniform_controls() {
    for collision in [Collision::Trt, Collision::Mrt, Collision::Regularized, Collision::Hrr] {
        for outlet in ["convective", "periodic"] {
            let worst = uniform(collision, 0.518, [4, 1, 1], 0, 1.0, outlet, 30_000).unwrap();
            println!("{collision:?} {outlet}: {worst:.4e}");
            assert!(worst <= 1e-12, "{worst}");
        }
    }
}

/// A zone may occupy the inner cell even though it cannot occupy the open layer.
/// Extrapolate its physical velocity, including the implicit porous or fan force.
#[test]
fn pressure_extrapolates_the_inner_zone_velocity() {
    for collision in [Collision::Trt, Collision::Mrt, Collision::Regularized, Collision::Hrr] {
        for fan in [false, true] {
            let (_, mut flow) = common::pressure_zone::neighbor(collision, fan);
            let mut worst = 0.0f64;
            for _ in 0..100 {
                flow.step().unwrap();
                let inside = flow.moments(flow.grid().index([6, 0, 0])).unwrap().1[0];
                let (rho, u) = flow.moments(flow.grid().index([7, 0, 0])).unwrap();
                worst = worst.max((u[0] - inside).abs());
                assert!((rho - 1.0).abs() < 1e-12);
                assert!(u[1].abs().max(u[2].abs()) < 1e-12);
            }
            println!("{collision:?} fan={fan}: physical normal-velocity mismatch {worst:.3e}");
            assert!(worst < 1e-12, "{worst}");
        }
    }
}
