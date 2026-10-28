//! Vreman's subgrid eddy viscosity (`turbulence: {"type": "vreman"}`, Vreman 2004) on the CPU:
//! laminar shear (plane Couette and Poiseuille flow, whose gradients have rank one) unchanged by
//! the model bit for bit with every collision, a decaying Taylor–Green vortex where it acts, the
//! wall models with it, its determinism across threads and the refusals of its constant; and the
//! hybrid collision's dynamic σ (Jacob et al. 2018, eqs. 48–50), which the model's eddy
//! viscosity sets.
//!
//! Each test prints its numbers (`cargo test -p fieldcraft-lbm --test vreman -- --nocapture`).

mod common;

use common::*;
use fieldcraft_lbm::{Flow, Limits, Setup};
use fieldcraft_model::{BoxFace, Collision, FaceCondition, HrrSigma, Model, Procedure, Step, Turbulence, WallModel};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};
const ALL: [Collision; 4] = [Collision::Trt, Collision::Regularized, Collision::Mrt, Collision::Hrr];

/// Vreman's model with Smagorinsky-equivalent constant `c`.
fn vreman(c: f64) -> Turbulence {
    Turbulence::Vreman { constant: c, prandtl: 0.9, wall_model: None }
}

/// `p` with turbulence model `t`.
fn with(mut p: Procedure, t: Option<Turbulence>) -> Procedure {
    if let Procedure::Fluid { turbulence, .. } = &mut p {
        *turbulence = t;
    }
    p
}

/// `p` with the hybrid collision and dynamic σ.
fn dynamic(mut p: Procedure) -> Procedure {
    if let Procedure::Fluid { collision, sigma, .. } = &mut p {
        *collision = Some(Collision::Hrr);
        *sigma = Some(HrrSigma::Dynamic);
    }
    p
}

/// Every fluid cell's density and velocity, as bits.
fn state(flow: &Flow) -> Vec<Option<(u64, [u64; 3])>> {
    (0..flow.grid().len()).map(|c| flow.moments(c).map(|(r, u)| (r.to_bits(), u.map(f64::to_bits)))).collect()
}

/// A plane channel `n` cells across (y), one cell along the periodic x and z: Couette (the top
/// wall moving) or Poiseuille (driven by a body force).
fn channel(n: usize, couette: bool) -> fieldcraft_model::FluidDomain {
    let (h, nu, umax) = (1.0, 0.05, 0.1);
    let dx = h / n as f64;
    let mut d = domain([0.0; 3], [dx, h, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    if couette {
        d.faces.y_max = FaceCondition::MovingWall { velocity: [umax, 0.0, 0.0] };
    } else {
        d.acceleration = [8.0 * nu * umax / (h * h), 0.0, 0.0];
    }
    d
}

/// In plane Couette and Poiseuille flow the only nonzero derivative is ∂_y u_x (the x and z
/// differences are of equal periodic images, exactly 0), so the gradient has one nonzero row and
/// Vreman's B_β is exactly 0: every cell keeps the fluid's own rate, and the flow is the one
/// without a model, bit for bit, with every collision. The reported largest ν_t/ν is exactly 0.
#[test]
fn couette_and_poiseuille_flow_are_unchanged_by_the_model_bit_for_bit() {
    let (n, nu) = (16, 0.05);
    let dx = 1.0 / n as f64;
    for couette in [true, false] {
        for collision in ALL {
            let run = |t: Option<Turbulence>| {
                let p = with(collide(2.0, dx * dx, None, collision), t);
                let (setup, mut flow) = prepare(&model(1.0, nu, channel(n, couette), p), &[]);
                run_all(&setup, &mut flow);
                flow
            };
            let (plain, modelled) = (run(None), run(Some(vreman(0.18))));
            assert_eq!(state(&plain), state(&modelled), "{collision:?}, Couette {couette}");
            assert_eq!(plain.forces(), modelled.forces());
            assert_eq!(modelled.eddy_viscosity_ratio(), 0.0);
            let u: Vec<f64> = (0..n).map(|j| modelled.moments(j).unwrap().1[0]).collect();
            println!("Vreman, {} ({collision:?}): equal bit for bit, ν_t = 0; u_x(y) {u:.4?}", if couette { "Couette" } else { "Poiseuille" });
        }
    }
}

/// A Taylor–Green vortex in a periodic box 2π (n³ cells, Re 50 so that the run without a model
/// is taken, lattice speed 0.05): the model
/// acts there (the gradient has rank two or three), its largest ν_t/ν is positive, and the
/// kinetic energy after a while is below the run's without a model (it only adds dissipation).
fn taylor_green(collision: Collision, t: Option<Turbulence>, n: usize, steps: usize) -> (f64, f64) {
    let (k, r, _) = taylor_green_with(|p| p, collision, t, n, steps);
    (k, r)
}

/// [`taylor_green`] with the procedure changed by `change`; also the dynamic σ's smallest and mean
/// value in the last step.
fn taylor_green_with(
    change: impl Fn(Procedure) -> Procedure,
    collision: Collision,
    t: Option<Turbulence>,
    n: usize,
    steps: usize,
) -> (f64, f64, Option<(f64, f64)>) {
    use std::f64::consts::PI;
    let side = 2.0 * PI;
    let dx = side / n as f64;
    let dt = 0.05 * dx;
    let mut d = domain([0.0; 3], [side; 3], dx);
    for f in BoxFace::ALL {
        *d.faces.get_mut(f) = PERIODIC;
    }
    let p = change(with(collide(steps as f64 * dt, dt, None, collision), t));
    let (setup, mut flow) = prepare(&model(1.0, 1.0 / 50.0, d, p), &[]);
    let grid = *flow.grid();
    let vs = setup.velocity_scale();
    flow.set_equilibrium(|cell| {
        let [x, y, z] = grid.centre(cell);
        let pressure = ((2.0 * x).cos() + (2.0 * y).cos()) * ((2.0 * z).cos() + 2.0) / 16.0;
        (1.0 + 3.0 * pressure / (vs * vs), [x.sin() * y.cos() * z.cos() / vs, -x.cos() * y.sin() * z.cos() / vs, 0.0])
    });
    run_all(&setup, &mut flow);
    let energy = (0..grid.len()).map(|c| flow.moments(c).map_or(0.0, |(_, u)| u.iter().map(|v| v * v).sum::<f64>())).sum::<f64>();
    (energy * vs * vs / (2.0 * grid.len() as f64), flow.eddy_viscosity_ratio(), flow.dynamic_sigma())
}

#[test]
fn the_taylor_green_vortex_takes_a_positive_eddy_viscosity_and_decays_faster() {
    for collision in [Collision::Trt, Collision::Hrr] {
        let (k0, r0) = taylor_green(collision, None, 16, 160);
        let (k, r) = taylor_green(collision, Some(vreman(0.18)), 16, 160);
        println!("Taylor–Green 16³ ({collision:?}), t 3.1: k {k0:.6} without a model, {k:.6} with Vreman's (largest ν_t/ν {r:.3})");
        assert_eq!(r0, 0.0);
        assert!(r > 0.0, "{r}");
        assert!(k < k0 && k > 0.0, "{k} vs {k0}");
    }
}

/// The same result bit for bit on one and on several threads, with the hybrid and the
/// two-relaxation-time collisions.
#[test]
fn results_are_bit_for_bit_the_same_on_one_and_many_threads() {
    for collision in [Collision::Trt, Collision::Hrr] {
        let (m, ball) = ball_channel(collision, None, Some(vreman(0.18)));
        let solve = |threads: usize| {
            fieldcraft_parallel::with_threads(threads, || {
                fieldcraft_lbm::solve_fluid(&m, None, std::slice::from_ref(&ball), &Limits::default()).unwrap()
            })
        };
        let one = solve(1);
        assert!(one.grid.len() > 3 * 4096, "several parallel blocks: {}", one.grid.len());
        assert!(one.stats.eddy_viscosity_ratio.is_some_and(|r| r > 0.0), "{:?}", one.stats.eddy_viscosity_ratio);
        assert_eq!(one, solve(7));
        assert_eq!(one, solve(3));
    }
}

/// Every collision runs with the model in a channel with a ball, an inlet, an outlet and a slip
/// face (every kind of neighbour), and the model acts.
#[test]
fn every_collision_runs_with_the_model() {
    for collision in ALL {
        let (m, ball) = ball_channel(collision, None, Some(vreman(0.18)));
        let (_, mut flow) = prepare(&m, std::slice::from_ref(&ball));
        assert!(flow.vreman() > 0.0 && flow.smagorinsky() == 0.0);
        for _ in 0..40 {
            flow.step().unwrap();
        }
        let ratio = flow.eddy_viscosity_ratio();
        println!("Vreman with {collision:?}, 40 steps: largest ν_t/ν {ratio:.4}");
        assert!(ratio > 0.0 && ratio.is_finite(), "{ratio}");
    }
}

/// The wall models take Vreman's model as they take Smagorinsky's: a coarse Couette flow under
/// Musker's wall model (the wall cells relax with the wall model's rate, so the largest ν_t/ν is
/// theirs) and a coarse forced channel under the wall-function bounce, a few hundred steps each.
#[test]
fn the_wall_models_run_with_the_model() {
    let n = 4;
    let dx = 1.0 / n as f64;
    for (couette, wall) in [(true, WallModel::Musker), (false, WallModel::SpaldingBounce)] {
        for collision in [Collision::Mrt, Collision::Regularized] {
            let t = Turbulence::Vreman { constant: 0.18, prandtl: 0.9, wall_model: Some(wall) };
            // The forced channel takes a shorter time step for its Mach number.
            let (dt, nu) = if couette { (dx, 1e-5) } else { (0.1 * dx, 1e-4) };
            let p = with(collide(300.0 * dt, dt, None, collision), Some(t));
            let (setup, mut flow) = prepare(&model(1.0, nu, channel(n, couette), p), &[]);
            assert_eq!((flow.wall_model(), flow.wall_bounce()), (wall == WallModel::Musker, wall == WallModel::SpaldingBounce));
            run_all(&setup, &mut flow);
            let u: Vec<f64> = (0..n).map(|j| flow.moments(j).unwrap().1[0] * setup.velocity_scale()).collect();
            println!("Vreman with {wall:?} ({collision:?}): u_x {u:.4?}, largest ν_t/ν {:.3}", flow.eddy_viscosity_ratio());
            assert!(u.iter().all(|v| v.is_finite()));
            if wall == WallModel::Musker {
                assert!(flow.eddy_viscosity_ratio() > 0.0);
            }
        }
    }
}

/// A constant out of (0, 0.5], or not a number, and a turbulent Prandtl number out of range are
/// refused, naming Vreman's model; Van Driest's damping is not a field of it. Nothing panics.
#[test]
fn bad_constants_are_refused() {
    let d = || {
        let mut d = domain([0.0; 3], [1.0, 1.0, 0.25], 0.25);
        d.faces.y_max = FaceCondition::MovingWall { velocity: [0.01, 0.0, 0.0] };
        d
    };
    for bad in [0.0, -0.1, 0.51, f64::NAN, f64::INFINITY] {
        let m = model(1.0, 0.1, d(), with(collide(0.1, 0.01, None, Collision::Trt), Some(vreman(bad))));
        let e = Setup::new(&m, None, &[], &Limits::default()).unwrap_err().to_string();
        println!("C_s = {bad}: {e}");
        assert!(e.contains("Vreman") && e.contains("at most 0.5"), "{e}");
    }
    let t = Turbulence::Vreman { constant: 0.18, prandtl: 0.0, wall_model: None };
    let e = Setup::new(&model(1.0, 0.1, d(), with(collide(0.1, 0.01, None, Collision::Trt), Some(t))), None, &[], &Limits::default())
        .unwrap_err()
        .to_string();
    assert!(e.contains("Prandtl"), "{e}");
    let setup =
        Setup::new(&model(1.0, 0.1, d(), with(collide(0.1, 0.01, None, Collision::Trt), Some(vreman(0.18)))), None, &[], &Limits::default()).unwrap();
    assert!((setup.vreman - 2.5 * 0.18 * 0.18).abs() < 1e-15 && setup.smagorinsky == 0.0, "c = 2.5 C_s²");
    let step = serde_json::json!({"type": "fluid", "domain": "Flow", "end_time": 0.1,
        "turbulence": {"type": "vreman", "constant": 0.18, "van_driest": 26}});
    assert!(serde_json::from_value::<Procedure>(step).is_err());
    let mut m = Model::default();
    assert!(m.add_step(Step::new("Run", with(collide(0.1, 0.01, None, Collision::Trt), Some(vreman(f64::NAN)))), None).is_err());
}

/// Dynamic σ where the model gives no eddy viscosity (Couette and Poiseuille flow with Vreman's
/// model) is σ = 1 in every cell: the recursive regularised collision, bit for bit. With the
/// subgrid viscosity out of τ, the flow is the regularised one without a model.
#[test]
fn dynamic_sigma_without_eddy_viscosity_is_the_regularised_collision_bit_for_bit() {
    let (n, nu) = (16, 0.05);
    let dx = 1.0 / n as f64;
    for couette in [true, false] {
        let run = |p: Procedure| {
            let (setup, mut flow) = prepare(&model(1.0, nu, channel(n, couette), p), &[]);
            run_all(&setup, &mut flow);
            flow
        };
        let regularized = run(collide(2.0, dx * dx, None, Collision::Regularized));
        let flow = run(dynamic(with(collide(2.0, dx * dx, None, Collision::Hrr), Some(vreman(0.18)))));
        assert_eq!(state(&regularized), state(&flow), "Couette {couette}");
        assert_eq!(regularized.forces(), flow.forces());
        assert_eq!(flow.dynamic_sigma(), Some((1.0, 1.0)));
        assert_eq!(regularized.dynamic_sigma(), None);
        println!("dynamic σ, {}: σ = 1 everywhere, the regularised collision bit for bit", if couette { "Couette" } else { "Poiseuille" });
    }
}

/// In the Taylor–Green vortex the model's eddy viscosity sets σ below 1, within [0, 1]; the
/// subgrid viscosity does not enter τ, and the run decays (Vreman's and Smagorinsky's models).
#[test]
fn dynamic_sigma_acts_in_the_taylor_green_vortex() {
    let (k0, _) = taylor_green(Collision::Regularized, None, 16, 160);
    let smagorinsky = Turbulence::Smagorinsky { constant: 0.17, prandtl: 0.9, van_driest: None, wall_model: None };
    for t in [vreman(0.18), smagorinsky] {
        let (k, ratio, sigma) = taylor_green_with(dynamic, Collision::Hrr, Some(t), 16, 160);
        let (min, mean) = sigma.unwrap();
        println!(
            "dynamic σ, Taylor–Green 16³, {}: k {k:.6} (regularised without a model {k0:.6}), σ min {min:.4}, mean {mean:.4}, ν_t/ν {ratio:.3}",
            t.name()
        );
        assert!((0.0..1.0).contains(&min) && min <= mean && mean < 1.0, "{min} {mean}");
        assert!(ratio > 0.0);
        assert!(k > 0.0 && k < k0 * 1.001, "{k} vs {k0}");
    }
}

/// From its equilibrium start the Taylor–Green vortex's flux is zero but for round-off, so
/// Smagorinsky's ν_t is 0 and dynamic σ is 1 in every cell of the first step, exactly. Before the
/// round-off floor ([`fieldcraft_lbm::turbulence::eddy_viscosity`]) the rounding gave ν_t > 0 where
/// ∇²u vanishes, and σ went towards 0 there.
#[test]
fn an_equilibrium_start_gives_dynamic_sigma_one_exactly() {
    let smagorinsky = Turbulence::Smagorinsky { constant: 0.17, prandtl: 0.9, van_driest: None, wall_model: None };
    let (_, ratio, sigma) = taylor_green_with(dynamic, Collision::Hrr, Some(smagorinsky), 16, 1);
    println!("dynamic σ, Taylor–Green 16³, first step from equilibrium: σ min and mean {sigma:?}, ν_t/ν {ratio:e}");
    assert_eq!(sigma, Some((1.0, 1.0)));
    assert_eq!(ratio, 0.0);
}

/// Dynamic σ gives the same result bit for bit on one and on several threads (its σ statistics
/// too: the blocks' sums are added in block order).
#[test]
fn dynamic_sigma_is_bit_for_bit_the_same_on_one_and_many_threads() {
    let smagorinsky = Turbulence::Smagorinsky { constant: 0.17, prandtl: 0.9, van_driest: None, wall_model: None };
    for t in [vreman(0.18), smagorinsky] {
        let (m, ball) = ball_channel(Collision::Hrr, Some(HrrSigma::Dynamic), Some(t));
        let run = |threads: usize| {
            fieldcraft_parallel::with_threads(threads, || {
                let (_, mut flow) = prepare(&m, std::slice::from_ref(&ball));
                for _ in 0..30 {
                    flow.step().unwrap();
                }
                (state(&flow), flow.dynamic_sigma().map(|(a, b)| (a.to_bits(), b.to_bits())))
            })
        };
        let one = run(1);
        assert!(one.1.is_some());
        assert_eq!(one, run(7));
        assert_eq!(one, run(3));
    }
}
