//! The same flow on one thread and on several, bit for bit (`fieldcraft-parallel`'s rules): the
//! plain two-relaxation-time collision (no model that changes a cell's state or rates) with
//! halfway, linear and quadratic interpolated walls, and heat with buoyancy (the temperatures and the buoyancy
//! forces worked out by blocks). Each grid spans several parallel blocks, and the wall links
//! more than one task of the wall forces.

mod common;

use common::*;
use fieldcraft_lbm::Flow;
use fieldcraft_model::{FaceCondition, FluidThermal, Outflow, Profile, ThermalCondition, WallTreatment};

/// Every bit of the state a step leaves: each cell's density and velocity (and temperature with
/// heat), each group's force and heat.
fn state(flow: &Flow) -> Vec<u64> {
    let grid = flow.grid();
    let cells = (0..grid.len()).flat_map(|c| {
        let moments = flow.moments(c).map_or([0.0; 4], |(rho, u)| [rho, u[0], u[1], u[2]]);
        let temperature = flow.temperature(c).unwrap_or(0.0);
        [moments[0], moments[1], moments[2], moments[3], temperature]
    });
    let forces = flow.forces().iter().flatten().copied();
    cells.chain(forces).chain(flow.heat().iter().copied()).map(f64::to_bits).collect()
}

/// `steps` steps of `m` around `obstacles` on each thread count, compared with one thread's.
fn same_on_every_thread_count(what: &str, m: &fieldcraft_model::Model, obstacles: &[fieldcraft_lbm::Obstacle], steps: usize) {
    let run = |threads: usize| {
        fieldcraft_parallel::with_threads(threads, || {
            let (_, mut flow) = prepare(m, obstacles);
            for _ in 0..steps {
                flow.step().unwrap();
            }
            flow
        })
    };
    let one = run(1);
    // Blocks of at least 512 cells (`flow::block_cells`), links by tasks of 4096.
    assert!(one.grid().len() >= 8 * 512, "{what}: several parallel blocks ({} cells)", one.grid().len());
    assert!(one.links().len() > 2 * 4096, "{what}: several tasks of wall links ({})", one.links().len());
    let reference = state(&one);
    assert!(reference.iter().any(|b| f64::from_bits(*b) != 0.0 && f64::from_bits(*b) != 1.0), "{what}: the flow moved");
    for threads in [2, 3, 7] {
        assert!(state(&run(threads)) == reference, "{what}: {threads} threads differ from one");
    }
}

/// A channel with a velocity inlet, a pressure outlet, a slip top, a body force and a block on
/// the floor: the plain collision with halfway, linear and quadratic interpolated walls.
#[test]
fn the_plain_collision_is_the_same_on_one_and_many_threads() {
    let dx = 0.05;
    let mut d = domain([0.0; 3], [1.6, 0.8, 0.7], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.2, 0.01, 0.0], profile: Profile::Parabolic };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
    d.faces.z_max = FaceCondition::Slip {};
    d.acceleration = [0.0, 0.0, -0.05];
    d.obstacles = vec!["Block".into()];
    let block = obstacle("Block", cuboid([0.43, 0.27, 0.0], [0.71, 0.52, 0.37]), &CUBOID_FACES);
    for walls in [WallTreatment::Halfway, WallTreatment::Interpolated, WallTreatment::InterpolatedQuadratic] {
        let m = model(1.0, 2e-3, d.clone(), with_walls(step(1.0, 0.01, None), walls));
        same_on_every_thread_count(&format!("{walls:?} walls"), &m, std::slice::from_ref(&block), 60);
    }
}

/// A cavity heated from one side with buoyancy (Boussinesq, through Guo's force term), so each
/// step works out the temperatures and the buoyancy forces cell by cell.
#[test]
fn heat_with_buoyancy_is_the_same_on_one_and_many_threads() {
    let (n, alpha, nu) = (32usize, 1.0, 0.71);
    let dx = 1.0 / n as f64;
    let mut d = domain([0.0; 3], [1.0, 1.0, 16.0 * dx], dx);
    let mut t =
        FluidThermal { initial_temperature: kelvin(0.5), gravity: [0.0, -1e4 * nu * alpha, 0.0], expansion: Some(1.0), ..FluidThermal::default() };
    t.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
    t.faces.x_max = ThermalCondition::Temperature { value: kelvin(0.0) };
    d.thermal = Some(t);
    let dt = dx * dx / (4.0 * alpha);
    let m = hot_model(1.0, nu, alpha, d, step(1.0, dt, None));
    same_on_every_thread_count("buoyant cavity", &m, &[], 60);
}
