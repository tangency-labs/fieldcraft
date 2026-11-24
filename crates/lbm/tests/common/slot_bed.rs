//! The slot bed: a plate with a slot, driven through (tests/trace_bulk.rs, the trace relaxation's gate).
//!
//! A plate two cells thick across a box that is periodic on every face, with a slot `slot` cells wide (three slots to the period),
//! the body force that drives a plane channel of the slot's width at the lattice speed `u_lat`.

use super::*;
use fieldcraft_lbm::{Limits, Obstacle, Setup};
use fieldcraft_model::{Collision, FaceCondition, HrrSigma, Procedure, Turbulence};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};
/// The box's length along x in cells, the plate's cells along x and its first cell.
pub const LX: usize = 32;
pub const THICK: usize = 2;
pub const PLATE: usize = 10;

pub struct Bed {
    pub setup: Setup,
    pub obstacles: Vec<Obstacle>,
    /// The acceleration of a plane channel of the slot's width at `u_lat`, in lattice units.
    pub g_lat: f64,
}

pub fn build(slot: usize, tau_minus_half: f64, u_lat: f64, sigma: f64, trace_bulk_one: bool) -> Result<Bed, String> {
    let nu_lat = tau_minus_half / 3.0;
    let dx = 1.0 / slot as f64;
    let dt = nu_lat * dx * dx;
    let period = 3 * slot;
    let (lx, ly, lz) = (LX as f64 * dx, period as f64 * dx, 4.0 * dx);
    let (x0, x1) = (PLATE as f64 * dx, (PLATE + THICK) as f64 * dx);
    let (ya, yb) = (slot as f64 * dx, 2.0 * slot as f64 * dx);
    let obstacles = vec![
        obstacle("PlateLow", cuboid([x0, -dx, -dx], [x1, ya, lz + dx]), &CUBOID_FACES),
        obstacle("PlateHigh", cuboid([x0, yb, -dx], [x1, ly + dx, lz + dx]), &CUBOID_FACES),
    ];
    let mut d = domain([0.0; 3], [lx, ly, lz], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.obstacles = obstacles.iter().map(|o| o.body.clone()).collect();
    let g_lat = 12.0 * nu_lat * u_lat / (slot * slot) as f64;
    d.acceleration = [g_lat * dx / (dt * dt), 0.0, 0.0];
    let mut p = collide(1e6 * dt, dt, None, Collision::Hrr);
    if let Procedure::Fluid { sigma: s, turbulence, trace_bulk_one: on, .. } = &mut p {
        *s = Some(HrrSigma::Value(sigma));
        *on = trace_bulk_one;
        if tau_minus_half < 5e-3 - 1e-12 {
            *turbulence = Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: None });
        }
    }
    let m = model(1.0, 1.0, d, p);
    let setup = Setup::new(&m, None, &obstacles, &Limits::default()).map_err(|e| e.to_string())?;
    Ok(Bed { setup, obstacles, g_lat })
}
