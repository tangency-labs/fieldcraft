//! The free surface: a liquid column at rest under gravity (hydrostatics), a plug filling a
//! frictionless channel from an inlet (the front at U t, the mass balance), a collapsing column
//! (mass conserved to round-off, the front moving, the same on any thread count), and what is
//! refused.

mod common;

use common::*;
use fieldcraft_lbm::{Flow, Kind, Limits, Setup, solve_fluid};
use fieldcraft_model::{FaceCondition, FluidDomain, FreeSurface, LiquidRegion, Model, Profile};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};
const SLIP: FaceCondition = FaceCondition::Slip {};

/// A box `nx × 1 × nz` cells, periodic along y, with a free surface whose liquid fills
/// `[lo, hi]` (cell coordinates, inclusive).
fn tank(nx: usize, nz: usize, dx: f64, liquid: &[([f64; 3], [f64; 3])]) -> FluidDomain {
    let mut d = domain([0.0; 3], [nx as f64 * dx, dx, nz as f64 * dx], dx);
    d.faces.y_min = PERIODIC;
    d.faces.y_max = PERIODIC;
    d.free_surface = Some(FreeSurface { liquid: liquid.iter().map(|(min, max)| LiquidRegion { min: *min, max: *max }).collect() });
    d
}

/// The liquid's height in each column of cells along z (the sum of the fill levels times Δx).
fn heights(flow: &Flow, nx: usize, nz: usize, dx: f64) -> Vec<f64> {
    let s = flow.surface().unwrap();
    (0..nx).map(|i| (0..nz).map(|k| s.fill(flow.grid().index([i, 0, k]))).sum::<f64>() * dx).collect()
}

/// A column of liquid at rest under gravity in a closed box stays at rest: the surface keeps its
/// height, the pressure grows linearly with depth to ρ g h at the bottom, and the mass is the
/// column's to round-off.
#[test]
fn a_liquid_column_under_gravity_is_hydrostatic() {
    let (nx, nz, dx) = (4usize, 40usize, 0.025);
    let h = 0.5; // half the box
    let g = 9.81;
    let mut d = tank(nx, nz, dx, &[([0.0; 3], [nx as f64 * dx, dx, h])]);
    d.acceleration = [0.0, 0.0, -g];
    // A lattice gravity of 1e-4 (the density grows by 0.6 % over the column) and τ = 0.6.
    let g_lattice = 1.0e-4;
    let dt = (g_lattice * dx / g).sqrt();
    let nu = (0.6 - 0.5) / 3.0 * dx * dx / dt;
    let m = model(1000.0, nu, d, step(4000.0 * dt, dt, None));
    let (setup, mut flow) = prepare(&m, &[]);
    assert!(flow.free_surface());
    let m0 = flow.mass();
    assert!((m0 - (nx * nz / 2) as f64).abs() < 1e-9, "the column's cells at density 1: {m0}");
    for _ in 0..setup.steps {
        flow.step().unwrap();
    }
    let report = flow.surface_report().unwrap();
    println!("hydrostatic: {report:?}");
    assert!(report.mass_error.abs() < 1e-12, "mass to round-off: {:e}", report.mass_error);
    // The weakly compressible liquid is denser by `g (h − z) / c_s²` at depth: its mean density
    // is `1 + 1.5 g_lattice h_cells`, so the column, its mass kept, stands that much lower.
    let compressed = h / (1.0 + 1.5 * g_lattice * (nz / 2) as f64);
    for height in heights(&flow, nx, nz, dx) {
        assert!((height - compressed).abs() < 0.02 * dx, "the surface stays at {compressed:.5} (h = {h}): {height}");
    }
    // The largest speed over the liquid, against the free-fall speed over the column.
    let vs = setup.velocity_scale();
    let mut peak = 0.0f64;
    for c in 0..flow.grid().len() {
        if flow.surface().unwrap().kind(c) == Some(Kind::Liquid)
            && let Some((_, u)) = flow.moments(c)
        {
            peak = peak.max((u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt() * vs);
        }
    }
    let free_fall = (2.0 * g * h).sqrt();
    println!("peak speed {peak:.3e} of free fall {free_fall:.3}");
    // measured 2.7e-4 of it on this grid (spurious currents at the interface cells)
    assert!(peak < 1e-3 * free_fall, "the liquid holds still: {peak}");
    // The pressure at the cell centres from the surface down: ρ g (h − z).
    let mut worst = 0.0f64;
    for k in 0..nz / 2 {
        let cell = flow.grid().index([1, 0, k]);
        let (rho, _) = flow.moments(cell).unwrap();
        let p = setup.gauge_pressure(rho);
        let z = (k as f64 + 0.5) * dx;
        let want = 1000.0 * g * (h - z);
        worst = worst.max((p - want).abs() / (1000.0 * g * h));
    }
    println!("hydrostatic pressure error {worst:.2e} of ρ g h");
    assert!(worst < 0.03, "ρ g (h − z) within 3 % of ρ g h: {worst}");
}

/// A channel with frictionless walls, periodic across, filled from a velocity inlet at U: the
/// liquid is a plug, its front at `U t`, the volume behind it `U t A`, and the mass the inlet
/// injected equals what the box holds.
#[test]
fn a_plug_fills_a_frictionless_channel_at_the_inlet_speed() {
    let dx = 1.0 / 64.0;
    let u = 0.05;
    let mut d = domain([0.0; 3], [1.0, dx, dx], dx);
    d.faces.y_min = PERIODIC;
    d.faces.y_max = PERIODIC;
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d.faces.x_min = FaceCondition::Velocity { velocity: [u, 0.0, 0.0], profile: Profile::Uniform };
    d.free_surface = Some(FreeSurface { liquid: vec![] });
    let dt = 0.05 * dx / u; // lattice speed 0.05
    let end = 0.6 / u; // the front at x = 0.6
    let m = model(1.0, 1e-3, d, step(end, dt, None));
    let r = solve_fluid(&m, None, &[], &Limits::default()).unwrap();
    let fs = r.stats.free_surface.unwrap();
    println!("plug: {fs:?}");
    assert!(fs.mass_error.abs() < 1e-11, "the mass balance closes: {:e}", fs.mass_error);
    assert!(fs.initial_mass > 0.0 && fs.initial_mass < 3.0 * dx.powi(3), "the inlet layer alone starts liquid: {}", fs.initial_mass);
    // The volume behind the front: the fill levels summed along the channel.
    let last = r.frames.last().unwrap();
    let filled: f64 = last.fill.iter().map(|e| f64::from(*e)).sum::<f64>() * dx;
    let want = u * end + dx; // the inlet cell holds liquid from the start
    println!("front at {filled:.4} (fill sum), U t + Δx = {want:.4}");
    assert!((filled - want).abs() < 0.5 * dx, "the front at U t within half a cell: {filled} vs {want}");
    // Half of the fill lies where the front is.
    let half = last.fill.iter().position(|e| *e < 0.5).unwrap() as f64 * dx;
    assert!((half - want).abs() < 1.5 * dx, "the half-full cell at the front: {half} vs {want}");
    for w in &r.warnings {
        assert!(!w.contains("free surface lost"), "{w}");
    }
    // Between slip faces (two cells across z) the plug is the same and the balance closes: the
    // populations reflected along the walls are paired with their mirror cells.
    let mut d = domain([0.0; 3], [1.0, dx, 2.0 * dx], dx);
    d.faces.y_min = PERIODIC;
    d.faces.y_max = PERIODIC;
    d.faces.z_min = SLIP;
    d.faces.z_max = SLIP;
    d.faces.x_min = FaceCondition::Velocity { velocity: [u, 0.0, 0.0], profile: Profile::Uniform };
    d.free_surface = Some(FreeSurface { liquid: vec![] });
    let m = model(1.0, 1e-3, d, step(end, dt, None));
    let r = solve_fluid(&m, None, &[], &Limits::default()).unwrap();
    let fs = r.stats.free_surface.unwrap();
    println!("plug between slip faces: {fs:?}");
    assert!(fs.mass_error.abs() < 1e-11, "the mass balance closes between slip faces: {:e}", fs.mass_error);
    let filled: f64 = r.frames.last().unwrap().fill.iter().map(|e| f64::from(*e)).sum::<f64>() * dx / 2.0;
    assert!((filled - want).abs() < 0.5 * dx, "the front at U t: {filled} vs {want}");
}

/// A column collapsing under gravity in a closed box: the mass stays what it was to round-off, the
/// front runs along the floor, the column's top comes down, and the result is the same bit for bit
/// on one and four threads.
#[test]
fn a_collapsing_column_keeps_its_mass_on_any_thread_count() {
    let (nx, nz, dx) = (48usize, 36usize, 0.0125);
    let a = 0.15; // the column: a wide, 2a high, under a box 3a high
    let g = 9.81;
    let mut d = tank(nx, nz, dx, &[([0.0; 3], [a, dx, 2.0 * a])]);
    d.acceleration = [0.0, 0.0, -g];
    let free_fall = (2.0 * g * 2.0 * a).sqrt();
    let dt = 0.04 * dx / free_fall;
    let nu = 0.02 * dx * dx / dt; // τ = 0.56
    let end = 0.6 * (a / g).sqrt() * 2.0;
    let m = model(1.0, nu, d, step(end, dt, None));
    let run = |threads| fieldcraft_parallel::with_threads(threads, || solve_fluid(&m, None, &[], &Limits::default()).unwrap());
    let r = run(1);
    let fs = r.stats.free_surface.unwrap();
    println!("dam break: {fs:?}; steps {}", r.stats.steps);
    assert!(fs.mass_error.abs() < 1e-11, "mass to round-off: {:e}", fs.mass_error);
    assert!(fs.filled > 0 && fs.emptied > 0, "cells converted both ways: {fs:?}");
    let front = |frame: &fieldcraft_lbm::Frame| {
        (0..nx).rev().find(|i| (0..nz).any(|k| frame.fill[r.grid.index([*i, 0, k])] > 0.5)).map_or(0.0, |i| (i as f64 + 1.0) * dx)
    };
    let (first, last) = (r.frames.first().unwrap(), r.frames.last().unwrap());
    println!("front {:.4} → {:.4}, top {:.4} → {:.4}", front(first), front(last), top(first, &r.grid, dx), top(last, &r.grid, dx));
    assert!((front(first) - a).abs() < 0.5 * dx);
    assert!(front(last) > a + 4.0 * dx, "the front ran along the floor");
    assert!(top(last, &r.grid, dx) < 2.0 * a - dx, "the column came down");
    let history = r.surface.as_ref().unwrap();
    assert_eq!(history.max.len(), r.history_times.len());
    assert!((history.max.last().unwrap()[0] - front(last)).abs() < 1e-12, "the extent history ends at the front");
    // the volume follows the mass to the weakly compressible density's deviation from 1
    assert!((history.volume[0] - first.mass).abs() < 1e-3 * first.mass, "the liquid's volume is its mass at density 1");
    assert!(r.frames.iter().all(|f| (f.mass - first.mass).abs() < 1e-9 * first.mass), "every frame's mass is the first's");
    assert_eq!(run(4), r, "bit for bit on four threads");
}

/// The liquid's highest cell at the back wall (x = 0).
fn top(frame: &fieldcraft_lbm::Frame, grid: &fieldcraft_lbm::Grid, dx: f64) -> f64 {
    let nz = grid.cells[2];
    (0..nz).rev().find(|k| frame.fill[grid.index([0, 0, *k])] > 0.5).map_or(0.0, |k| (k as f64 + 1.0) * dx)
}

/// What a free surface refuses by name, and what the model refuses.
#[test]
fn unsupported_combinations_are_refused_by_name() {
    let refused = |m: &Model, what: &str| {
        let e = Setup::new(m, None, &[], &Limits::default()).unwrap_err().to_string();
        assert!(e.contains("free surface") && e.contains(what), "{what}: {e}");
    };
    let refuse = |edit: &dyn Fn(&mut FluidDomain), what: &str| {
        let mut d = tank(8, 8, 0.1, &[([0.0; 3], [0.8, 0.1, 0.4])]);
        edit(&mut d);
        refused(&model(1.0, 0.01, d, step(1.0, 0.1, None)), what);
    };
    refuse(&|d| d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Default::default(), sigma: None }, "pressure face");
    let mut hot = tank(8, 8, 0.1, &[([0.0; 3], [0.8, 0.1, 0.4])]);
    hot.thermal = Some(Default::default());
    refused(&hot_model(1.0, 0.01, 0.01, hot, step(1.0, 0.1, None)), "heat");
    refuse(
        &|d| {
            d.porous.push(fieldcraft_model::PorousZone {
                name: "Vent".into(),
                min: [0.2, 0.0, 0.0],
                max: [0.3, 0.1, 0.8],
                axis: fieldcraft_model::Axis::X,
                resistance: fieldcraft_model::Resistance::LossCoefficient { coefficient: 1.0 },
            })
        },
        "porous",
    );
    // No liquid anywhere: refused when the domain is defined.
    let mut m = Model::default();
    m.add_material(
        "Fluid",
        fieldcraft_model::Material { density: Some(1.0), fluid: Some(fieldcraft_model::Fluid::Newtonian { viscosity: 0.01 }), ..Default::default() },
    )
    .unwrap();
    let mut d = tank(8, 8, 0.1, &[]);
    d.free_surface = Some(FreeSurface { liquid: vec![] });
    let e = m.add_fluid_domain("Dry", d).unwrap_err().to_string();
    assert!(e.contains("no liquid"), "{e}");
    // A region too thin for a cell centre.
    let d = tank(8, 8, 0.1, &[([0.0, 0.0, 0.0], [0.8, 0.1, 0.01])]);
    let m = model(1.0, 0.01, d, step(1.0, 0.1, None));
    let e = Setup::new(&m, None, &[], &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("no cell centre"), "{e}");
    // Slip faces run.
    let mut d = tank(8, 8, 0.1, &[([0.0; 3], [0.8, 0.1, 0.4])]);
    d.faces.z_max = SLIP;
    let m = model(1.0, 0.01, d, step(1.0, 0.1, None));
    assert!(Setup::new(&m, None, &[], &Limits::default()).is_ok());
}
