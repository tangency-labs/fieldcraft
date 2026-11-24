//! The Spalart–Allmaras model coupled to the flow on the CPU (`turbulence: {"type":
//! "spalart_allmaras"}`): ν̃ transported on the cells that collide, its ν_t = ν̃ f_v1 in the
//! collision's rate. With no ν̃ the flow is the laminar one bit for bit; a forced channel takes a
//! positive ν_t that stays finite; the result is the same bit for
//! bit on any number of threads; dynamic σ is refused, naming the model (the GPU's own
//! refusals are lbm-gpu's tests).
//!
//! Each test prints its numbers (`cargo test -p fieldcraft-lbm --test spalart_allmaras -- --nocapture`).

mod common;

use common::*;
use fieldcraft_lbm::sa::{WallLaw, law_eddy_viscosity};
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
    channel(n, true)
}

/// A plane channel `n` cells across (y), one cell along the periodic x and z: Couette (the top
/// wall moving) or Poiseuille (driven by a body force, peak speed 0.1 at ν 0.05).
fn channel(n: usize, couette: bool) -> fieldcraft_model::FluidDomain {
    let dx = 1.0 / n as f64;
    let mut d = domain([0.0; 3], [dx, 1.0, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    if couette {
        d.faces.y_max = FaceCondition::MovingWall { velocity: [0.1, 0.0, 0.0] };
    } else {
        d.acceleration = [8.0 * 0.05 * 0.1, 0.0, 0.0];
    }
    d
}

/// With ν̃ = 0 at the start and no inflow, every source of ν̃ is 0 (production and destruction are
/// proportional to ν̃, diffusion and advection of a zero field vanish), so ν̃ stays exactly 0, ν_t
/// is 0 and the cells keep the fluid's own rates: the laminar flow bit for bit, in plane Couette
/// and body-force Poiseuille flow, with the two-relaxation-time and hybrid collisions.
#[test]
fn with_no_nu_tilde_the_flow_is_laminar_bit_for_bit() {
    let n = 16;
    let dx = 1.0 / n as f64;
    for couette in [true, false] {
        for collision in [Collision::Trt, Collision::Hrr] {
            let run = |t: Option<Turbulence>| {
                let (setup, mut flow) = prepare(&model(1.0, 0.05, channel(n, couette), with(collide(2.0, dx * dx, None, collision), t)), &[]);
                run_all(&setup, &mut flow);
                flow
            };
            let (plain, modelled) = (run(None), run(Some(sa(0.0, None))));
            assert!(modelled.turbulent() && !plain.turbulent());
            let field = modelled.spalart_allmaras().unwrap();
            assert!(field.nu_tilde().iter().all(|v| v.to_bits() == 0), "ν̃ stays exactly 0");
            assert_eq!(state(&plain), state(&modelled), "{collision:?}, Couette {couette}");
            assert_eq!(modelled.eddy_viscosity_ratio(), 0.0);
            let flow = if couette { "Couette" } else { "Poiseuille" };
            println!("SA with ν̃ = 0, {flow} ({collision:?}): equal to the laminar flow bit for bit");
        }
    }
}

/// A forced, wall-resolved channel from ν̃/ν = 3 everywhere: ν̃ stays finite, the largest ν_t/ν
/// is positive, and the first fluid cell's wall distance is half a cell (the halfway wall), with
/// each collision.
#[test]
fn a_forced_channel_takes_a_positive_finite_eddy_viscosity() {
    for collision in [Collision::Trt, Collision::Hrr, Collision::Mrt] {
        let (setup, mut flow) = prepare(&resolved_channel(collision, sa(3.0, None)), &[]);
        run_all(&setup, &mut flow);
        let nu = flow.spalart_allmaras().unwrap().nu_tilde().to_vec();
        let d = flow.sa_wall_distance().unwrap();
        let grid = *flow.grid();
        let column: Vec<usize> = (0..grid.cells[1]).map(|j| grid.index([0, j, 0])).collect();
        let profile: Vec<f64> = column.iter().map(|&c| nu[c]).collect();
        let distance: Vec<f64> = column.iter().map(|&c| d[c]).collect();
        let ratio = flow.eddy_viscosity_ratio();
        println!("SA, {collision:?}: largest ν_t/ν {ratio:.4}; ν̃ (lattice) {:?}; d {:.2?}", &profile[..8], &distance[..8]);
        assert!(nu.iter().all(|v| v.is_finite()));
        assert!(ratio > 0.0 && ratio.is_finite(), "{ratio}");
        assert_eq!(distance[0], 0.5, "the first fluid cell's wall distance");
        assert_eq!(distance[grid.cells[1] - 1], 0.5);
    }
}

/// A channel at Re_τ 1000 (half-height 1, u_τ 1 from the body force, ν 10⁻³), 24 cells per
/// half-height so the first fluid node sits near y⁺ 40, started at its bulk speed.
fn re_tau_1000(collision: Collision, wall: WallModel, on_node: bool) -> fieldcraft_model::Model {
    let dx = 1.0 / 24.0;
    let height = if on_node { 2.0 + dx } else { 2.0 };
    let mut d = domain([0.0; 3], [4.0 * dx, height, 4.0 * dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.acceleration = [1.0, 0.0, 0.0];
    d.initial_velocity = [20.0, 0.0, 0.0];
    let dt = 0.1 * dx / 20.0;
    model(1.0, 1e-3, d, with(collide(100.0 * dt, dt, None, collision), Some(sa(3.0, Some(wall)))))
}

/// The wall function (Nichols & Nelson 2004): at Re_τ 1000 on Han, Ooka & Kikumoto's on-node
/// wall (Spalding's law) and under Musker's wall model (the control), every first fluid node's
/// ν_t is the law's at its y⁺ with the u_τ the wall applied this step (bar 1e-12; the
/// inverse of ν̃ f_v1 gives it to round-off), u_τ is positive and y⁺ past the buffer
/// layer. Without a wall model there are no samples.
#[test]
fn the_wall_function_sets_the_first_node_eddy_viscosity() {
    for (collision, wall, on_node) in [(Collision::Hrr, WallModel::SpaldingBounce, true), (Collision::Hrr, WallModel::Musker, false)] {
        let (setup, mut flow) = prepare(&re_tau_1000(collision, wall, on_node), &[]);
        run_all(&setup, &mut flow);
        let nu = (setup.tau_plus - 0.5) / 3.0;
        let (law, samples) = flow.sa_wall_function();
        let (law, samples) = (law.unwrap(), samples.to_vec());
        assert_eq!(law, if on_node { WallLaw::Spalding } else { WallLaw::Musker });
        let field = flow.spalart_allmaras().unwrap();
        assert!(!samples.is_empty());
        let grid = *flow.grid();
        let (mut worst, mut y_plus) = (0.0f64, (f64::INFINITY, 0.0f64));
        for &(cell, y, ut) in &samples {
            assert!(ut > 0.0 && ut.is_finite(), "u_τ {ut} at {cell}");
            let want = nu * law_eddy_viscosity(law, y * ut / nu).unwrap();
            if law == WallLaw::Musker {
                // The wall model's own eddy viscosity: wall_law at the speed whose u_τ this is.
                let speed = ut * fieldcraft_lbm::wall::musker(3.0 * y * ut / nu).0;
                let own = fieldcraft_lbm::wall::wall_law(speed, y, 3.0 * y, nu).1;
                assert!((want / own - 1.0).abs() < 1e-9, "wall_law's {own} vs {want}");
            }
            worst = worst.max((field.nu_t(cell, nu) / want - 1.0).abs());
            y_plus = (y_plus.0.min(y * ut / nu), y_plus.1.max(y * ut / nu));
            let j = grid.coords(cell)[1];
            assert!(j == usize::from(on_node) || j == grid.cells[1] - 1 - usize::from(on_node), "a first fluid node: row {j}");
        }
        println!(
            "SA wall function, {wall:?}: {} first nodes, y⁺ {:.1}–{:.1}, worst |ν_t/ν_t,law − 1| {worst:.2e}",
            samples.len(),
            y_plus.0,
            y_plus.1
        );
        assert!(worst < 1e-12, "{worst}");
        assert!(y_plus.0 > 20.0 && y_plus.1 < 100.0, "{y_plus:?}");
    }
    let (setup, mut flow) = prepare(&resolved_channel(Collision::Trt, sa(3.0, None)), &[]);
    flow.step().unwrap();
    assert!(flow.sa_wall_function().0.is_none() && flow.sa_wall_function().1.is_empty() && setup.steps > 0);
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

/// Dynamic σ and an inflow ratio out of range are refused with the model's name and what to do
/// instead; nothing panics. The memory estimate counts the model's field.
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
    // τ⁺ = 0.5019: SA keeps the laminar floor 0.505 whatever ν̃ is (its ν_t is 0 at a wall), so the grid is refused.
    let low = |ratio: f64| {
        Setup::new(
            &model(1.0, 0.004, couette(4), with(collide(0.1, 0.01, None, Collision::Trt), Some(sa(ratio, None)))),
            None,
            &[],
            &Limits::default(),
        )
    };
    for ratio in [0.0, 3.0] {
        let e = low(ratio).unwrap_err().to_string();
        assert!(e.contains("below 0.505") && !e.contains("add a turbulence model"), "{ratio}: {e}");
    }
    // The model's own wall law takes only Spalart–Allmaras turbulence; with a subgrid model it is refused.
    let e = refused(
        Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: Some(WallModel::SpalartAllmarasBounce) },
        Collision::Trt,
        None,
    );
    assert!(e.contains("spalart_allmaras") && e.contains("spalding_bounce"), "{e}");
    // `musker` beside an obstacle: the first-node wall function is built for the planar box faces.
    let mut d = slab_column(8, [300.0, 300.0], 0.0);
    d.thermal = None;
    let m = model(1.0, 0.1, d, with(collide(0.1, 0.01, None, Collision::Hrr), Some(sa(3.0, Some(WallModel::Musker)))));
    let obstacles = [slab_obstacle([0.25, 0.5])];
    let setup = Setup::new(&m, None, &obstacles, &Limits::default()).unwrap();
    let e = Flow::new(&setup, &obstacles).err().map(|e| e.to_string()).unwrap_or_default();
    assert!(e.contains("not supported yet") && e.contains("musker") && e.contains("obstacles"), "{e}");
    let m = model(1.0, 0.1, couette(4), with(collide(0.1, 0.01, None, Collision::Trt), Some(sa(3.0, None))));
    let (setup, _) = prepare(&m, &[]);
    let plain = model(1.0, 0.1, couette(4), collide(0.1, 0.01, None, Collision::Trt));
    let (with_sa, without) = (Memory::estimate(&m, None, Choice::Cpu).unwrap(), Memory::estimate(&plain, None, Choice::Cpu).unwrap());
    let cells = with_sa.cells;
    println!("memory of {cells} cells: {} B with SA, {} B without", with_sa.total(), without.total());
    assert!(with_sa.total() >= without.total() + cells * fieldcraft_lbm::sa::BYTES_PER_CELL, "{with_sa:?} vs {without:?}");
    // The wall distance's transform, on this grid one cell thick several times the cell count, on top of the prepared flow.
    let transform = fieldcraft_lbm::distance::transform_bytes(setup.grid.cells.map(|c| c as u64));
    assert!(transform > 3 * cells * 8, "{transform} bytes for {cells} cells");
    assert!(with_sa.setup >= transform + cells * fieldcraft_lbm::sa::BYTES_PER_CELL, "{with_sa:?}: {transform} bytes of transform");
    let of = Memory::of(&setup, Choice::Cpu).unwrap();
    assert_eq!(of.total(), with_sa.total(), "the setup's and the model's estimate agree");
}
