//! Shared CPU/GPU #820 reproductions and controls.
use super::*;
use fieldcraft_lbm::Backend;
use fieldcraft_model::{FaceCondition, Outflow, Profile};

pub fn box_flow(collision: Collision, tau: f64, sides: &str, seed: f64, size: [usize; 3], axis: usize, sign: f64) -> Flow {
    let mut d = domain([0.0; 3], size.map(|n| n as f64), 1.0);
    let side = if sides == "slip" { FaceCondition::Slip {} } else { FaceCondition::Periodic {} };
    for face in fieldcraft_model::BoxFace::ALL {
        *d.faces.get_mut(face) = side;
    }
    let velocity = std::array::from_fn(|a| if a == axis { 0.05 * sign } else { 0.0 });
    let pair = &fieldcraft_model::BoxFace::ALL[2 * axis..2 * axis + 2];
    if sides == "closed" {
        for face in pair {
            *d.faces.get_mut(*face) = FaceCondition::Periodic {};
        }
    } else {
        *d.faces.get_mut(pair[usize::from(sign < 0.0)]) = FaceCondition::Velocity { velocity, profile: Profile::Uniform };
        *d.faces.get_mut(pair[usize::from(sign > 0.0)]) = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective };
    }
    d.initial_velocity = velocity;
    let (_, mut flow) = prepare(&model(1.0, (tau - 0.5) / 3.0, d, collide(30_000.0, 1.0, None, collision)), &[]);
    let grid = *flow.grid();
    flow.set_equilibrium(|cell| {
        let mut u = velocity;
        u[axis] += if grid.coords(cell)[(axis + 1) % 3] % 2 == 0 { seed } else { -seed };
        (1.0, u)
    });
    flow
}

/// Deterministic broad-band seed in all three velocity components.
pub fn broadband(flow: &mut Flow, axis: usize, sign: f64, seed: f64) {
    flow.set_equilibrium(|cell| {
        (
            1.0,
            std::array::from_fn(|a| {
                let mut h = ((cell * 3 + a + 1) as u64).wrapping_mul(0x9e3779b97f4a7c15);
                h = (h ^ (h >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
                h = (h ^ (h >> 27)).wrapping_mul(0x94d049bb133111eb);
                let noise = ((h ^ (h >> 31)) >> 11) as f64 / ((1u64 << 53) as f64) * 2.0 - 1.0;
                (if a == axis { 0.05 * sign } else { 0.0 }) + seed * noise
            }),
        )
    });
}

pub fn departures(flow: &mut impl Backend, axis: usize, sign: f64) -> [f64; 2] {
    let fields = flow.fields().unwrap();
    let mut peak = [0.0f64; 2];
    for (rho, u) in fields.density.iter().zip(fields.velocity) {
        peak[0] = peak[0].max((rho - 1.0).abs());
        for (a, v) in u.into_iter().enumerate() {
            peak[1] = peak[1].max((v - if a == axis { 0.05 * sign } else { 0.0 }).abs());
        }
    }
    peak
}
