//! A velocity inlet under the MRT collision stays stable. The Zou-He closure leaves the
//! non-hydrodynamic part of the inlet cell's non-equilibrium unconstrained; the regularised
//! collisions project it every step, MRT relaxes it raw, so without the projection of the inlet
//! cell (`zou_he`, `project`) a plate-with-slot channel at low viscosity went unstable in
//! about 80 steps, in the first columns behind the inlet. This channel (uniform inlet at 0.05 in
//! lattice units with a cross-flow of 0.02, pressure outlet, no-slip walls on y) runs 3000 steps at
//! tau - 1/2 = 1e-2, 6e-3 and 5e-3 and must stay finite with a bounded speed (on develop's inlet
//! the first goes unstable at step 855). The other collisions do not change: the populations of the same channel
//! after 300 steps hash to what the code gave before the projection existed.

mod common;

use common::*;
use fieldcraft_model::{Collision, FaceCondition, Outflow, Profile};

fn flow(tau_minus_half: f64, collision: Collision, steps: u64) -> fieldcraft_lbm::Flow {
    let nu = tau_minus_half / 3.0;
    let u = 0.05;
    let (nx, ny, nz) = (40usize, 12usize, 4usize);
    let mut d = domain([0.0; 3], [nx as f64, ny as f64, nz as f64], 1.0);
    for f in [&mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.faces.x_min = FaceCondition::Velocity { velocity: [u, 0.02, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.initial_velocity = [u, 0.0, 0.0];
    let p = collide(steps as f64, 1.0, None, collision);
    let (_, mut flow) = prepare(&model(1.0, nu, d, p), &[]);
    for _ in 0..steps {
        flow.step().unwrap();
    }
    flow
}

fn peak_speed(flow: &fieldcraft_lbm::Flow) -> f64 {
    let mut peak = 0.0f64;
    for c in 0..flow.grid().len() {
        if flow.is_fluid(c)
            && let Some((_, v)) = flow.moments(c)
        {
            peak = peak.max((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt());
        }
    }
    peak
}

#[test]
fn a_velocity_inlet_under_mrt_stays_stable_at_low_viscosity() {
    for tau in [1e-2, 6e-3, 5e-3] {
        let peak = peak_speed(&flow(tau, Collision::Mrt, 3000));
        assert!(peak.is_finite() && peak < 0.2, "tau - 1/2 {tau}: peak speed {peak}");
    }
}

/// FNV-1a over the bits of every population of every cell.
fn hash(flow: &fieldcraft_lbm::Flow) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for c in 0..flow.grid().len() {
        for q in 0..19 {
            for b in flow.population(q, c).to_bits().to_le_bytes() {
                h = (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3);
            }
        }
    }
    h
}

#[test]
fn the_other_collisions_are_what_they_were() {
    let got: Vec<(Collision, u64)> =
        [Collision::Trt, Collision::Regularized, Collision::Hrr].into_iter().map(|c| (c, hash(&flow(1e-2, c, 300)))).collect();
    println!("{got:#018x?}");
    for (c, h) in got {
        let want = match c {
            Collision::Trt => TRT,
            Collision::Regularized => REGULARIZED,
            _ => HRR,
        };
        assert_eq!(h, want, "{c:?}: {h:#018x}");
    }
}

const TRT: u64 = 0x4c81fcc62090ed4d;
const REGULARIZED: u64 = 0x54a1868475433d25;
const HRR: u64 = 0x61336226143a4f8d;
