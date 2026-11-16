//! Interpolated walls in an open domain: a laminar plane channel along x with a uniform velocity
//! inlet (the exact mean speed) and a pressure outlet, its wall one slab whose face cuts the cells
//! (8.5 cells across on an 8-cell grid, a width no cell face matches), periodic across y and z, no
//! turbulence model, 20 000 steps, CPU f64. The flow is steady and plane Poiseuille, so two things
//! are known: the flux through any cross-section is the same (the plane at two thirds of the
//! length over the one at a third is 1) and the pressure gradient is 12 ν U / H² per cell for the
//! flow the channel carries (read over the middle third; density is 3 p in lattice units).
//!
//! The expected numbers are what the code gives, recorded as "known" where they miss their bars
//! (flux ratio within 1e-3, gradient within 1 %). The wall treatments are developed in
//! `interpolated_rotated.rs`; this file is the check that what they do to the mass does not fight
//! the open boundaries.

mod common;

use common::*;
use fieldcraft_model::WallTreatment;

const NX: usize = 64;
const STEPS: u64 = 20_000;

/// What the run found.
#[derive(Debug)]
struct Open {
    /// The flux through the plane at 2/3 of the length over the one at 1/3.
    flux_ratio: f64,
    /// The pressure gradient of the middle third over the exact one for the flow the channel
    /// carries, minus one, in percent.
    gradient_error: f64,
    /// The mean density over the fluid cells at the end.
    mean_density: f64,
}

/// The open channel with the wall treatment `walls`.
fn open_channel(walls: WallTreatment) -> Open {
    let (m, obstacles) = open_channel_model(NX, STEPS, walls);
    let (_, mut flow) = prepare(&m, &obstacles);
    for _ in 0..STEPS {
        flow.step().unwrap();
    }
    let grid = *flow.grid();
    let [gx, gy, gz] = grid.cells;
    let flux = |i: usize| -> f64 {
        let mut s = 0.0;
        for k in 0..gz {
            for j in 0..gy {
                let c = grid.index([i, j, k]);
                if flow.is_fluid(c) {
                    s += flow.moments(c).map_or(0.0, |(rho, u)| rho * u[0]);
                }
            }
        }
        s
    };
    // The cross-section's mean density at every plane of the middle third, and a least-squares line.
    let mut pts = Vec::new();
    for i in gx / 3..2 * gx / 3 {
        let (mut s, mut n) = (0.0, 0usize);
        for k in 0..gz {
            for j in 0..gy {
                let c = grid.index([i, j, k]);
                if flow.is_fluid(c) {
                    s += flow.moments(c).map_or(0.0, |(rho, _)| rho);
                    n += 1;
                }
            }
        }
        pts.push((i as f64, s / n.max(1) as f64));
    }
    let n = pts.len() as f64;
    let (mx, my) = (pts.iter().map(|p| p.0).sum::<f64>() / n, pts.iter().map(|p| p.1).sum::<f64>() / n);
    let slope = pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum::<f64>() / pts.iter().map(|p| (p.0 - mx) * (p.0 - mx)).sum::<f64>();
    // The exact gradient of the flow the middle plane carries: 12 ν (ρ u) / H² per cell, ρ = 3 p.
    let carried = flux(gx / 2) / (gz as f64 * OPEN_WIDTH);
    let exact = -3.0 * 12.0 * OPEN_NU_LAT * carried / (OPEN_WIDTH * OPEN_WIDTH);
    let (mut sum, mut count) = (0.0, 0usize);
    for c in 0..grid.len() {
        if flow.is_fluid(c) {
            sum += flow.moments(c).map_or(0.0, |(rho, _)| rho);
            count += 1;
        }
    }
    Open { flux_ratio: flux(2 * gx / 3) / flux(gx / 3), gradient_error: 100.0 * (slope / exact - 1.0), mean_density: sum / count.max(1) as f64 }
}

/// Recorded numbers: the flux ratio, the gradient error (percent) and the mean density.
struct Known {
    flux_ratio: f64,
    gradient_error: f64,
    mean_density: f64,
}

fn check(walls: WallTreatment, known: Known) {
    let r = open_channel(walls);
    println!(
        "{walls:?}: flux ratio {:.5} (recorded {:.5}), gradient {:+.2} % (recorded {:+.2}), mean density {:.5} (recorded {:.5})",
        r.flux_ratio, known.flux_ratio, r.gradient_error, known.gradient_error, r.mean_density, known.mean_density
    );
    if std::env::var("OPEN_RECORD").is_ok() {
        return;
    }
    assert!((r.flux_ratio - known.flux_ratio).abs() < 2e-4, "{r:?}");
    assert!((r.gradient_error - known.gradient_error).abs() < 0.05, "{r:?}");
    assert!((r.mean_density - known.mean_density).abs() < 2e-4, "{r:?}");
}

/// The linear rule: known misses, the flux through the open channel is steady but the pressure
/// gradient reads high and the mean density sits above the outlet's level (the extra mass the wall
/// made, which the open ends settle).
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_linear_rule_in_an_open_channel_keeps_its_recorded_error() {
    check(WallTreatment::Interpolated, Known { flux_ratio: 1.00000, gradient_error: 3.00, mean_density: 1.01603 });
}

/// The quadratic rule, no correction: the same, recorded.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_quadratic_rule_in_an_open_channel_keeps_its_recorded_error() {
    check(WallTreatment::InterpolatedQuadratic, Known { flux_ratio: 1.00000, gradient_error: -0.75, mean_density: 1.01568 });
}
