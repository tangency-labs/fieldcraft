//! The Spalart–Allmaras model coupled to the flow on the CPU (`turbulence: {"type":
//! "spalart_allmaras"}`): ν̃ transported on the cells that collide, its ν_t = ν̃ f_v1 in the
//! collision's rate. With no ν̃ the flow is the laminar one bit for bit; a forced channel takes a
//! positive ν_t that stays finite, with halfway and on-node walls; the result is the same bit for
//! bit on any number of threads; dynamic σ and the GPU are refused, naming the model.
//!
//! Each test prints its numbers (`cargo test -p fieldcraft-lbm --test spalart_allmaras -- --nocapture`).

mod common;

use common::*;
use fieldcraft_lbm::{Choice, Flow, Limits, Memory, Setup};
use fieldcraft_model::{Collision, FaceCondition, HrrSigma, Procedure, Turbulence, WallModel};

/// The model with inflow ν̃/ν `ratio` and wall model `wall_model`.
fn sa(ratio: f64, wall_model: Option<WallModel>) -> Turbulence {
    Turbulence::SpalartAllmaras { prandtl: 0.9, wall_model, inflow_ratio: ratio }
}

/// `p` with turbulence model `t`.
fn with(mut p: Procedure, t: Option<Turbulence>) -> Procedure {
    if let Procedure::Fluid { turbulence, .. } = &mut p {
        *turbulence = t;
    }
    p
}

/// Every fluid cell's density and velocity, as bits.
fn state(flow: &Flow) -> Vec<Option<(u64, [u64; 3])>> {
    (0..flow.grid().len()).map(|c| flow.moments(c).map(|(r, u)| (r.to_bits(), u.map(f64::to_bits)))).collect()
}

/// A plane Couette channel `n` cells across (y), one cell along the periodic x and z.
fn couette(n: usize) -> fieldcraft_model::FluidDomain {
    let dx = 1.0 / n as f64;
    let mut d = domain([0.0; 3], [dx, 1.0, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.faces.y_max = FaceCondition::MovingWall { velocity: [0.1, 0.0, 0.0] };
    d
}

/// A forced channel two high at Re_τ about 180, six cells per half-height (as the checkpoint
/// tests'), periodic in x and z; `on_node` adds the layer of boundary nodes the on-node bounce
/// takes.
fn forced_channel(collision: Collision, turbulence: Turbulence, on_node: bool) -> fieldcraft_model::Model {
    let dx = 1.0 / 6.0;
    let height = if on_node { 2.0 + dx } else { 2.0 };
    let mut d = domain([0.0; 3], [4.0 * dx, height, 4.0 * dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.acceleration = [1.0, 0.0, 0.0];
    d.initial_velocity = [15.7, 0.0, 0.0];
    let dt = 0.1 * dx / 15.7;
    model(1.0, 1.0 / 180.0, d, with(collide(60.0 * dt, dt, None, collision), Some(turbulence)))
}

/// With ν̃ = 0 at the start and no inflow, every source of ν̃ is 0 (production and destruction are
/// proportional to ν̃, diffusion and advection of a zero field vanish), so ν̃ stays exactly 0, ν_t
/// is 0 and the cells keep the fluid's own rates: the laminar flow bit for bit, with the
/// two-relaxation-time and hybrid collisions.
#[test]
fn with_no_nu_tilde_the_flow_is_laminar_bit_for_bit() {
    let n = 16;
    let dx = 1.0 / n as f64;
    for collision in [Collision::Trt, Collision::Hrr] {
        let run = |t: Option<Turbulence>| {
            let (setup, mut flow) = prepare(&model(1.0, 0.05, couette(n), with(collide(2.0, dx * dx, None, collision), t)), &[]);
            run_all(&setup, &mut flow);
            flow
        };
        let (plain, modelled) = (run(None), run(Some(sa(0.0, None))));
        assert!(modelled.turbulent() && !plain.turbulent());
        let field = modelled.spalart_allmaras().unwrap();
        assert!(field.nu_tilde().iter().all(|v| v.to_bits() == 0), "ν̃ stays exactly 0");
        assert_eq!(state(&plain), state(&modelled), "{collision:?}");
        assert_eq!(modelled.eddy_viscosity_ratio(), 0.0);
        println!("SA with ν̃ = 0 ({collision:?}): equal to the laminar flow bit for bit");
    }
}

/// A forced channel from ν̃/ν = 3 everywhere: ν̃ stays finite, the largest ν_t/ν is positive, and
/// the wall distances are the walls' (half a cell to the halfway wall, a full cell to the on-node
/// wall's boundary node), with each collision and wall treatment.
#[test]
fn a_forced_channel_takes_a_positive_finite_eddy_viscosity() {
    let cases = [(Collision::Trt, None, false), (Collision::Hrr, None, false), (Collision::Mrt, Some(WallModel::SpaldingBounce), true)];
    for (collision, wall_model, on_node) in cases {
        let (setup, mut flow) = prepare(&forced_channel(collision, sa(3.0, wall_model), on_node), &[]);
        run_all(&setup, &mut flow);
        let nu = flow.spalart_allmaras().unwrap().nu_tilde().to_vec();
        let d = flow.sa_wall_distance().unwrap();
        let grid = *flow.grid();
        let column: Vec<usize> = (0..grid.cells[1]).map(|j| grid.index([0, j, 0])).collect();
        let profile: Vec<f64> = column.iter().map(|&c| nu[c]).collect();
        let distance: Vec<f64> = column.iter().map(|&c| d[c]).collect();
        let ratio = flow.eddy_viscosity_ratio();
        println!("SA, {collision:?}, on-node {on_node}: largest ν_t/ν {ratio:.4}; ν̃ (lattice) {profile:?}; d {distance:.2?}");
        assert!(nu.iter().all(|v| v.is_finite()));
        assert!(ratio > 0.0 && ratio.is_finite(), "{ratio}");
        let first = if on_node { 1 } else { 0 };
        let expected = if on_node { 1.0 } else { 0.5 };
        assert_eq!(distance[first], expected, "the first fluid cell's wall distance");
        assert_eq!(distance[grid.cells[1] - 1 - first], expected);
    }
}

/// The same result bit for bit on one and on several threads, in a channel with a ball, an inlet
/// (ν̃ fixed), an outlet and a slip face, over several parallel blocks.
#[test]
fn results_are_bit_for_bit_the_same_on_one_and_many_threads() {
    let (m, ball) = ball_channel(Collision::Trt, None, Some(sa(3.0, None)));
    let solve = |threads: usize| {
        fieldcraft_parallel::with_threads(threads, || fieldcraft_lbm::solve_fluid(&m, None, std::slice::from_ref(&ball), &Limits::default()).unwrap())
    };
    let one = solve(1);
    assert!(one.grid.len() > 3 * 4096, "several parallel blocks: {}", one.grid.len());
    let ratio = one.stats.eddy_viscosity_ratio;
    println!("SA, ball channel: largest ν_t/ν {ratio:?}");
    assert!(ratio.is_some_and(|r| r > 0.0 && r.is_finite()), "{ratio:?}");
    assert_eq!(one, solve(7));
    assert_eq!(one, solve(3));
}

/// Dynamic σ, an inflow ratio out of range and the GPU are refused with the model's name and
/// what to do instead; nothing panics. The memory estimate counts the model's field.
#[test]
fn what_is_not_supported_is_refused_and_the_memory_counts_the_field() {
    let refused = |t: Turbulence, collision: Collision, sigma: Option<HrrSigma>| {
        let mut p = with(collide(0.1, 0.01, None, collision), Some(t));
        if let Procedure::Fluid { sigma: s, .. } = &mut p {
            *s = sigma;
        }
        Setup::new(&model(1.0, 0.1, couette(4), p), None, &[], &Limits::default()).unwrap_err().to_string()
    };
    let e = refused(sa(3.0, None), Collision::Hrr, Some(HrrSigma::Dynamic));
    assert!(e.contains("not supported yet") && e.contains("dynamic σ") && e.contains("Spalart–Allmaras"), "{e}");
    // Interpolated walls put the wall where the surface is, not at the half-way link the distance reads.
    let mut p = with(collide(0.1, 0.01, None, Collision::Trt), Some(sa(3.0, None)));
    if let Procedure::Fluid { walls, .. } = &mut p {
        *walls = fieldcraft_model::WallTreatment::Interpolated;
    }
    let e = Setup::new(&model(1.0, 0.1, couette(4), p), None, &[], &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("not supported yet") && e.contains("spalart_allmaras"), "{e}");
    for bad in [-1.0, 101.0, f64::NAN, f64::INFINITY] {
        let e = refused(sa(bad, None), Collision::Trt, None);
        assert!(e.contains("inflow ratio"), "{bad}: {e}");
    }
    // τ⁺ = 0.5019: the eddy viscosity of ν̃ > 0 lowers the floor to 0.5001; with no ν̃ anywhere there is none, and the grid is refused.
    let low = |ratio: f64| {
        Setup::new(
            &model(1.0, 0.004, couette(4), with(collide(0.1, 0.01, None, Collision::Trt), Some(sa(ratio, None)))),
            None,
            &[],
            &Limits::default(),
        )
    };
    assert!(low(3.0).is_ok());
    let e = low(0.0).unwrap_err().to_string();
    assert!(e.contains("below 0.505") && e.contains("turbulence"), "{e}");
    let m = model(1.0, 0.1, couette(4), with(collide(0.1, 0.01, None, Collision::Trt), Some(sa(3.0, None))));
    let (setup, flow) = prepare(&m, &[]);
    match fieldcraft_lbm::backend::backend(Choice::Device, flow, &setup) {
        Ok(_) => panic!("the GPU must refuse the model"),
        Err(e) => {
            let e = e.to_string();
            assert!(e.contains("not supported yet") && e.contains("Spalart–Allmaras") && e.contains("GPU"), "{e}");
        }
    }
    let plain = model(1.0, 0.1, couette(4), collide(0.1, 0.01, None, Collision::Trt));
    let (with_sa, without) = (Memory::estimate(&m, None, Choice::Cpu).unwrap(), Memory::estimate(&plain, None, Choice::Cpu).unwrap());
    let cells = with_sa.cells;
    println!("memory of {cells} cells: {} B with SA, {} B without", with_sa.total(), without.total());
    assert!(with_sa.total() >= without.total() + cells * fieldcraft_lbm::sa::BYTES_PER_CELL, "{with_sa:?} vs {without:?}");
    let of = Memory::of(&setup, Choice::Cpu).unwrap();
    assert_eq!(of.total(), with_sa.total(), "the setup's and the model's estimate agree");
}
