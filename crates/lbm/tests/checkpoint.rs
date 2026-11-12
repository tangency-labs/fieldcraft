//! Restart fidelity uses independent original analytic cases, not changed numerical targets.
mod common;

use common::*;
use fieldcraft_lbm::{Error, Limits, Obstacle, checkpoint::Store, solve::solve_fluid_checkpointed, solve_fluid};
use fieldcraft_model::{
    Averaging, BoxFace, Collision, FaceCondition, Fan, FanDirection, FlowMonitor, FluidBoundary, FluidThermal, HrrSigma, Material, Model, Outflow,
    PhaseChange, Procedure, Steady, Step, ThermalCondition, Turbulence, WallModel,
};

#[derive(Default)]
struct MemoryStore {
    bytes: Vec<u8>,
    resume: bool,
    stop_after: usize,
    writes: usize,
    interval: u64,
}
impl Store for MemoryStore {
    fn every(&self) -> u64 {
        if self.interval == 0 { 11 } else { self.interval }
    }
    fn load(&mut self, f: &mut dyn FnMut(&mut dyn std::io::Read) -> fieldcraft_lbm::Result<u64>) -> fieldcraft_lbm::Result<bool> {
        if !self.resume {
            return Ok(false);
        }
        f(&mut std::io::Cursor::new(&self.bytes))?;
        self.resume = false;
        Ok(true)
    }
    fn save(&mut self, _step: u64, f: &mut dyn FnMut(&mut dyn std::io::Write) -> fieldcraft_lbm::Result<()>) -> fieldcraft_lbm::Result<()> {
        let mut next = Vec::new();
        f(&mut next)?;
        self.bytes = next;
        self.writes += 1;
        if self.writes == self.stop_after {
            return Err(Error::Invalid("test interruption after an atomic save".into()));
        }
        Ok(())
    }
}

fn parity(m: &fieldcraft_model::Model, obstacles: &[Obstacle]) {
    parity_at(m, obstacles, 11);
}

fn parity_at(m: &fieldcraft_model::Model, obstacles: &[Obstacle], interval: u64) {
    let normal = fieldcraft_parallel::with_threads(1, || solve_fluid(m, None, obstacles, &Limits::default())).unwrap();
    let mut uninterrupted = MemoryStore { interval, ..Default::default() };
    assert_eq!(solve_fluid_checkpointed(m, None, obstacles, &Limits::default(), &mut uninterrupted).unwrap(), normal);
    let mut interrupted = MemoryStore { stop_after: 2, interval, ..Default::default() };
    assert!(
        solve_fluid_checkpointed(m, None, obstacles, &Limits::default(), &mut interrupted).unwrap_err().to_string().contains("test interruption")
    );
    interrupted.resume = true;
    interrupted.stop_after = 0;
    // Cadence is storage policy, not solver physics; it can change at the restart.
    interrupted.interval = interval + 2;
    let resumed =
        fieldcraft_parallel::with_threads(2, || solve_fluid_checkpointed(m, None, obstacles, &Limits::default(), &mut interrupted)).unwrap();
    assert_eq!(resumed, normal, "all saved frames, histories, statistics and stop decisions");
    assert_eq!(interrupted.bytes, uninterrupted.bytes, "every retained float bit in final state and histories");
    // A completed checkpoint returns its result without advancing or writing again.
    interrupted.resume = true;
    let writes = interrupted.writes;
    assert_eq!(solve_fluid_checkpointed(m, None, obstacles, &Limits::default(), &mut interrupted).unwrap(), normal);
    assert_eq!(interrupted.writes, writes);
}

#[test]
fn a_convective_outlet_and_relaxed_fan_resume_with_every_bit_and_history() {
    for collision in COLLISIONS {
        let mut d = domain([0.0; 3], [12.0, 2.0, 2.0], 1.0);
        d.faces.x_min = FaceCondition::Velocity { velocity: [0.02, 0.0, 0.0], profile: fieldcraft_model::Profile::Uniform };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective, sigma: None };
        for f in [&mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
            *f = FaceCondition::Periodic {};
        }
        d.fans.push(Fan {
            name: "Fan".into(),
            min: [4.0, 0.0, 0.0],
            max: [6.0, 2.0, 2.0],
            direction: FanDirection::XPlus,
            curve: vec![[0.0, 0.001], [1.0, 0.0]],
            relaxation: Some(0.05),
        });
        parity(&model(1.0, 0.1, d, collide(67.0, 1.0, None, collision)), &[]);
    }
}

#[test]
fn thermal_buoyancy_and_conducting_solid_state_resume_inside_a_steady_window() {
    let dx = 0.125;
    let mut d = domain([0.0; 3], [dx, 1.0, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.obstacles.push("Slab".into());
    let mut thermal = FluidThermal::default();
    thermal.faces.y_min = ThermalCondition::Temperature { value: kelvin(0.0) };
    thermal.faces.y_max = ThermalCondition::Temperature { value: kelvin(1.0) };
    thermal.solids.push(conducting("Slab", 0.01));
    thermal.gravity = [0.0, -0.001, 0.0];
    thermal.expansion = Some(0.01);
    d.thermal = Some(thermal);
    let mut p = step(0.67, 0.01, None);
    if let Procedure::Fluid { steady, .. } = &mut p {
        *steady = Some(Steady { window: Some(0.16), tolerance: 1e-6, min_time: None, divergence: None, monitors: vec![], averaging: None });
    }
    let m = cht_model([1.0, 0.1, 0.05], (1.0, [0.2; 3], 1.0), d, p);
    let o = obstacle("Slab", cuboid([-1.0, 0.25, -1.0], [2.0, 0.5, 2.0]), &CUBOID_FACES);
    parity(&m, &[o]);
}

#[test]
fn heated_smagorinsky_flow_retains_averages_before_and_after_the_start() {
    for collision in COLLISIONS {
        // Non-sample-aligned starts exercise the window shape change and None changes before
        // a full averaging window. The save at step22 is before or after the first average.
        for start in [0.0, 0.17, 0.4] {
            let dx = 0.125;
            let mut d = domain([0.0; 3], [dx, 1.0, dx], dx);
            for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
                *f = FaceCondition::Periodic {};
            }
            d.faces.y_max = FaceCondition::MovingWall { velocity: [0.1, 0.0, 0.0] };
            let mut thermal = FluidThermal::default();
            thermal.faces.y_min = ThermalCondition::Temperature { value: kelvin(0.0) };
            thermal.faces.y_max = ThermalCondition::Temperature { value: kelvin(1.0) };
            thermal.gravity = [0.0, -0.001, 0.0];
            thermal.expansion = Some(0.01);
            d.thermal = Some(thermal);
            let mut p = collide(0.67, 0.01, None, collision);
            if let Procedure::Fluid { turbulence, steady, .. } = &mut p {
                *turbulence = Some(Turbulence::Smagorinsky { constant: 0.2, prandtl: 0.9, van_driest: None, wall_model: None });
                *steady = Some(Steady {
                    window: Some(0.16),
                    tolerance: 1e-6,
                    min_time: None,
                    divergence: None,
                    averaging: Some(Averaging { start }),
                    monitors: vec![
                        FlowMonitor::Force { faces: vec![FluidBoundary::Box(BoxFace::YMax)] },
                        FlowMonitor::Heat { faces: vec![FluidBoundary::Box(BoxFace::YMax)] },
                        FlowMonitor::MeanTemperature {},
                    ],
                });
            }
            let m = hot_model(1.0, 0.1, 0.05, d, p);
            let normal = solve_fluid(&m, None, &[], &Limits::default()).unwrap();
            assert!(normal.stats.eddy_viscosity_ratio.unwrap() > 0.0);
            let report = normal.steady.as_ref().unwrap();
            assert_eq!(report.averages.len(), 5);
            assert!(report.averages.iter().all(|a| a.samples > 0));
            parity(&m, &[]);
        }
    }
}

#[test]
fn averaged_steady_stop_and_parallel_flow_resume_across_thread_counts() {
    let mut d = domain([0.0; 3], [8.0, 16.0, 8.0], 1.0);
    // 1024 cells cross two fixed collision blocks; a quiescent flow stops on two full
    // post-start windows, with an unscheduled final frame, on any thread count.
    for f in BoxFace::ALL {
        *d.faces.get_mut(f) = FaceCondition::Periodic {};
    }
    let mut p = step(137.0, 1.0, None);
    if let Procedure::Fluid { turbulence, steady, .. } = &mut p {
        *turbulence = Some(Turbulence::Smagorinsky { constant: 0.2, prandtl: 0.9, van_driest: None, wall_model: None });
        *steady = Some(Steady {
            window: Some(16.0),
            tolerance: 1e-6,
            min_time: None,
            divergence: None,
            monitors: vec![],
            averaging: Some(Averaging { start: 17.0 }),
        });
    }
    let m = model(1.0, 0.1, d, p);
    let normal = solve_fluid(&m, None, &[], &Limits::default()).unwrap();
    assert_eq!(normal.steady.as_ref().unwrap().reason, fieldcraft_lbm::StopReason::Steady);
    assert_eq!(normal.stats.steps, 64);
    assert_eq!(normal.steady.as_ref().unwrap().averages[0].samples, 24);
    parity(&m, &[]);
}

#[test]
fn truncated_checkpoint_and_wrong_shape_are_refused() {
    let m = model(1.0, 0.1, domain([0.0; 3], [2.0; 3], 1.0), step(67.0, 1.0, None));
    let mut store = MemoryStore { stop_after: 2, ..Default::default() };
    assert!(solve_fluid_checkpointed(&m, None, &[], &Limits::default(), &mut store).is_err());
    let saved = store.bytes.clone();
    for len in [0, 1, 16, saved.len() - 1] {
        store.bytes = saved[..len].to_vec();
        store.resume = true;
        assert!(solve_fluid_checkpointed(&m, None, &[], &Limits::default(), &mut store).is_err(), "truncation {len}");
    }
    store.bytes = saved;
    store.resume = true;
    // The absolute step is the second u64, before any retained arrays can be allocated.
    store.bytes[8..16].copy_from_slice(&u64::MAX.to_le_bytes());
    assert!(solve_fluid_checkpointed(&m, None, &[], &Limits::default(), &mut store).is_err());
}

#[test]
fn two_consecutive_steady_windows_and_the_unscheduled_final_frame_survive_resume() {
    let mut p = step(67.0, 1.0, None);
    if let Procedure::Fluid { steady, .. } = &mut p {
        *steady = Some(Steady { window: Some(16.0), tolerance: 1e-6, min_time: None, divergence: None, monitors: vec![], averaging: None });
    }
    let m = model(1.0, 0.1, domain([0.0; 3], [2.0; 3], 1.0), p);
    let normal = solve_fluid(&m, None, &[], &Limits::default()).unwrap();
    assert_eq!(normal.steady.as_ref().unwrap().reason, fieldcraft_lbm::StopReason::Steady);
    assert_eq!(normal.stats.steps, 32);
    assert_eq!(normal.frames.last().unwrap().step, 32);
    parity(&m, &[]);
}

#[test]
fn quiet_windows_before_a_late_minimum_are_retained_even_beyond_two() {
    let mut p = step(137.0, 1.0, None);
    if let Procedure::Fluid { steady, .. } = &mut p {
        *steady = Some(Steady { window: Some(16.0), tolerance: 1e-6, min_time: Some(80.0), divergence: None, monitors: vec![], averaging: None });
    }
    let m = model(1.0, 0.1, domain([0.0; 3], [2.0; 3], 1.0), p);
    assert_eq!(solve_fluid(&m, None, &[], &Limits::default()).unwrap().stats.steps, 80);
    parity(&m, &[]);
}

#[test]
fn checkpoints_split_existing_advance_chunks_without_changing_bits() {
    // History samples every two steps and progress every64. 1009 is coprime to both;
    // the restarted cadence1011 differs too. Only eight cells keep this fast.
    let mut d = domain([0.0; 3], [2.0; 3], 1.0);
    d.faces.y_max = FaceCondition::MovingWall { velocity: [0.01, 0.0, 0.0] };
    let m = model(1.0, 0.1, d, step(10003.0, 1.0, None));
    parity_at(&m, &[], 1009);
}

#[test]
fn an_advance_failure_keeps_the_uninterrupted_reporting_boundary() {
    let dx = 0.1;
    let mut d = domain([0.0; 3], [dx; 3], dx);
    for f in fieldcraft_model::BoxFace::ALL {
        *d.faces.get_mut(f) = FaceCondition::Periodic {};
    }
    d.acceleration = [1.0, 0.0, 0.0];
    let dt = (0.009 * dx).sqrt();
    let p = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: None,
        steps: Some(10003),
        time_step: Some(dt),
        frames: 20,
        magic: None,
        collision: None,
        sigma: None,
        turbulence: None,
        steady: Some(Steady { window: Some(160.0 * dt), tolerance: 1e-6, min_time: None, divergence: None, monitors: vec![], averaging: None }),
    };
    let m = model(1.0, 0.01, d, p);
    let normal = solve_fluid(&m, None, &[], &Limits::default()).unwrap();
    assert_eq!(normal.steady.as_ref().unwrap().reason, fieldcraft_lbm::StopReason::Diverged);
    let mut store = MemoryStore { interval: 55, ..Default::default() };
    assert_eq!(solve_fluid_checkpointed(&m, None, &[], &Limits::default(), &mut store).unwrap(), normal);
    assert!(!store.bytes.is_empty(), "sound state55 retained before failure56");
    store.resume = true;
    assert_eq!(solve_fluid_checkpointed(&m, None, &[], &Limits::default(), &mut store).unwrap(), normal);
}

#[test]
fn rounded_steady_windows_can_retain_more_than_the_history_sampling_target() {
    //169999/10000 rounds up to17; eight sub-samples round that window down to16.
    //This valid existing monitor therefore keeps10625 windows, not at most10000.
    let p = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: None,
        steps: Some(169999),
        time_step: Some(1.0),
        frames: 2,
        magic: None,
        collision: None,
        sigma: None,
        turbulence: None,
        steady: Some(Steady { window: Some(17.0), tolerance: 1e-6, min_time: Some(169999.0), divergence: None, monitors: vec![], averaging: None }),
    };
    let m = model(1.0, 0.1, domain([0.0; 3], [1.0; 3], 1.0), p);
    let normal = solve_fluid(&m, None, &[], &Limits::default()).unwrap();
    assert!(normal.steady.as_ref().unwrap().samples.len() > fieldcraft_lbm::solve::MAX_HISTORY as usize + 1);
    parity_at(&m, &[], 50000);
}

/// A channel periodic in x and z, its walls the y faces, driven along x from its bulk speed
/// (Re_τ about 180, six cells per half-height, as the wall-model tests use; the first cells off
/// the walls past y⁺ 5), with `turbulence`, `collision` and
/// `sigma`. `on_node` adds the cell layer the on-node bounce takes for its boundary nodes.
fn wall_channel(collision: Collision, sigma: Option<HrrSigma>, turbulence: Turbulence, on_node: bool) -> Model {
    let dx = 1.0 / 6.0;
    let height = if on_node { 2.0 + dx } else { 2.0 };
    let mut d = domain([0.0; 3], [4.0 * dx, height, 4.0 * dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.acceleration = [1.0, 0.0, 0.0];
    d.initial_velocity = [15.7, 0.0, 0.0];
    let dt = 0.1 * dx / 15.7;
    let mut p = collide(60.0 * dt, dt, None, collision);
    if let Procedure::Fluid { turbulence: t, sigma: s, .. } = &mut p {
        *t = Some(turbulence);
        *s = sigma;
    }
    model(1.0, 1.0 / 180.0, d, p)
}

/// #711's wall treatments and #696's collisions keep their per-step state across a restart: the
/// wall model (Musker) with Van Driest's damping under MRT, the
/// on-node bounce under the hybrid collision with a dynamic σ, and Vreman's model with the wall
/// model under the regularised collision. A resumed finished run reports the same wall y⁺.
#[test]
fn wall_models_and_every_collision_resume_bit_for_bit() {
    let smagorinsky = |van_driest, wall_model| Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest, wall_model };
    let cases = [
        (Collision::Mrt, None, smagorinsky(Some(26.0), Some(WallModel::Musker)), false),
        (Collision::Hrr, Some(HrrSigma::Dynamic), smagorinsky(Some(26.0), Some(WallModel::SpaldingBounce)), true),
        (Collision::Regularized, None, Turbulence::Vreman { constant: 0.18, prandtl: 0.9, wall_model: Some(WallModel::Musker) }, false),
    ];
    for (collision, sigma, turbulence, on_node) in cases {
        let m = wall_channel(collision, sigma, turbulence, on_node);
        let normal = solve_fluid(&m, None, &[], &Limits::default()).unwrap();
        assert!(normal.stats.eddy_viscosity_ratio.unwrap() > 0.0, "{collision:?}: the subgrid model is active");
        let damped = matches!(turbulence, Turbulence::Smagorinsky { van_driest: Some(_), .. });
        assert_eq!(normal.warnings.iter().any(|w| w.contains("y⁺")), damped, "{collision:?}: the wall y⁺ report {:?}", normal.warnings);
        parity(&m, &[]);
    }
}

/// A shear-thinning channel clamps its relaxation time in the middle, where it is at rest; the
/// clamped share over the run in the warning counts the steps before the restart too.
#[test]
fn shear_dependent_viscosity_clamps_count_the_whole_run_after_resume() {
    for collision in [Collision::Trt, Collision::Mrt] {
        let dx = 1.0 / 8.0;
        let mut d = domain([0.0; 3], [dx, 1.0, dx], dx);
        for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
            *f = FaceCondition::Periodic {};
        }
        // A power law of index ½ whose plane Poiseuille flow peaks at 0.1 with ν 1/6 at the walls.
        let (n, h, peak, wall_nu) = (0.5f64, 0.5f64, 0.1f64, 1.0 / 6.0);
        let root = peak * (n + 1.0) / n / h.powf((n + 1.0) / n);
        let k = wall_nu / (root * h.powf(1.0 / n)).powf(n - 1.0);
        d.acceleration = [k * root.powf(n), 0.0, 0.0];
        let mut m = Model::default();
        let fluid = fieldcraft_model::Fluid::PowerLaw { consistency: k, index: n, reference_shear_rate: None };
        m.add_material("Fluid", Material { density: Some(1.0), fluid: Some(fluid), ..Material::default() }).unwrap();
        m.add_fluid_domain("Flow", d).unwrap();
        m.add_step(Step::new("Run", collide(0.5, dx * dx, None, collision)), None).unwrap();
        let normal = solve_fluid(&m, None, &[], &Limits::default()).unwrap();
        assert!(normal.warnings.iter().any(|w| w.contains("capped")), "{:?}", normal.warnings);
        parity(&m, &[]);
    }
}

/// A channel between a wall held at −1 and one at +1 whose fluid melts at −0.5, driven along x:
/// the layer that freezes on the cold wall during the run holds still, and the liquid fractions
/// that the next collision's resistance reads are restored, not rebuilt from the initial state.
#[test]
fn a_freezing_channel_resumes_with_its_frozen_layer() {
    let cells = 16usize;
    let dx = 1.0 / cells as f64;
    let mut d = domain([0.0; 3], [dx, 1.0, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.acceleration = [1e-3, 0.0, 0.0];
    let mut t = FluidThermal { initial_temperature: kelvin(0.0), ..Default::default() };
    t.faces.y_min = ThermalCondition::Temperature { value: kelvin(-1.0) };
    t.faces.y_max = ThermalCondition::Temperature { value: kelvin(1.0) };
    d.thermal = Some(t);
    let alpha = 0.05;
    let dt = dx * dx / (8.0 * alpha);
    let mut m = Model::default();
    let material = Material {
        phase_change: Some(PhaseChange { latent_heat: 1.0, solidus: kelvin(-0.5), liquidus: kelvin(-0.5) }),
        ..hot_fluid(1.0, 0.1, alpha)
    };
    m.add_material("Fluid", material).unwrap();
    m.add_fluid_domain("Flow", d).unwrap();
    m.add_step(Step::new("Run", step(80.0 * dt, dt, None)), None).unwrap();
    let normal = solve_fluid(&m, None, &[], &Limits::default()).unwrap();
    assert!(normal.stats.frozen.is_some_and(|f| f > 0.0), "part of the fluid freezes: {:?}", normal.stats.frozen);
    parity(&m, &[]);
}

#[test]
fn a_collapsing_column_resumes_with_its_surface() {
    let dx = 0.025;
    let mut d = domain([0.0; 3], [0.4, dx, 0.3], dx);
    d.faces.y_min = FaceCondition::Periodic {};
    d.faces.y_max = FaceCondition::Periodic {};
    d.acceleration = [0.0, 0.0, -9.81];
    d.free_surface = Some(fieldcraft_model::FreeSurface { liquid: vec![fieldcraft_model::LiquidRegion { min: [0.0; 3], max: [0.15, dx, 0.2] }] });
    let dt = 0.04 * dx / (2.0 * 9.81 * 0.2f64).sqrt();
    let m = model(1.0, 0.02 * dx * dx / dt, d, step(60.0 * dt, dt, None));
    let normal = solve_fluid(&m, None, &[], &Limits::default()).unwrap();
    let fs = normal.stats.free_surface.unwrap();
    assert!(fs.filled > 0 && fs.emptied > 0, "cells convert within the run: {fs:?}");
    parity(&m, &[]);
}
