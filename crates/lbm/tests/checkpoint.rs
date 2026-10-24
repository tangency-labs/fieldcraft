//! Restart fidelity uses independent original analytic cases, not changed numerical targets.
mod common;

use common::*;
use fieldcraft_lbm::{Error, Limits, Obstacle, checkpoint::Store, solve::solve_fluid_checkpointed, solve_fluid};
use fieldcraft_model::{FaceCondition, Fan, FanDirection, FluidThermal, Outflow, Procedure, Steady, ThermalCondition};

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
    let normal = solve_fluid(m, None, obstacles, &Limits::default()).unwrap();
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
    let resumed = solve_fluid_checkpointed(m, None, obstacles, &Limits::default(), &mut interrupted).unwrap();
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
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective };
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
        *steady = Some(Steady { window: Some(0.16), tolerance: 1e-6, min_time: None, divergence: None, monitors: vec![] });
    }
    let m = cht_model([1.0, 0.1, 0.05], (1.0, [0.2; 3], 1.0), d, p);
    let o = obstacle("Slab", cuboid([-1.0, 0.25, -1.0], [2.0, 0.5, 2.0]), &CUBOID_FACES);
    parity(&m, &[o]);
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
        *steady = Some(Steady { window: Some(16.0), tolerance: 1e-6, min_time: None, divergence: None, monitors: vec![] });
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
        *steady = Some(Steady { window: Some(16.0), tolerance: 1e-6, min_time: Some(80.0), divergence: None, monitors: vec![] });
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
        steady: Some(Steady { window: Some(160.0 * dt), tolerance: 1e-6, min_time: None, divergence: None, monitors: vec![] }),
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
        steady: Some(Steady { window: Some(17.0), tolerance: 1e-6, min_time: Some(169999.0), divergence: None, monitors: vec![] }),
    };
    let m = model(1.0, 0.1, domain([0.0; 3], [1.0; 3], 1.0), p);
    let normal = solve_fluid(&m, None, &[], &Limits::default()).unwrap();
    assert!(normal.steady.as_ref().unwrap().samples.len() > fieldcraft_lbm::solve::MAX_HISTORY as usize + 1);
    parity_at(&m, &[], 50000);
}
