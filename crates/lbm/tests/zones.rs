//! Fans and porous zones against analytic values: a porous plate across a plug flow drops the
//! pressure by the Darcy–Forchheimer law (or its loss coefficient), and a fan pushing a periodic
//! loop through such a plate settles at the intersection of its curve with the plate's
//! resistance curve, conserving mass.
//!
//! The flows are one-dimensional (periodic across the channel), so the superficial velocity is
//! the same everywhere and the resistance law is the whole pressure drop. Model and lattice units
//! coincide (ρ = 1, Δx = Δt = 1).

mod common;

use common::{COLLISIONS, collide, domain, model, prepare, run_all, step};
use fieldcraft_lbm::{Flow, Limits, Setup, solve_fluid};
use fieldcraft_model::{Axis, Collision, FaceCondition, Fan, FanDirection, PorousZone, Profile, Resistance};

const NU: f64 = 0.1;

fn plate(x0: f64, x1: f64, resistance: Resistance) -> PorousZone {
    PorousZone { name: "Vent".into(), min: [x0, 0.0, 0.0], max: [x1, 2.0, 2.0], axis: Axis::X, resistance }
}

/// A 40 × 2 × 2 channel periodic across, a uniform inlet at `u`, a pressure outlet, and a porous
/// plate over x 18..22 (four cell layers): the measured pressure drop and flow through it.
fn plate_drop(u: f64, resistance: Resistance, collision: Collision) -> (f64, f64) {
    let mut d = domain([0.0; 3], [40.0, 2.0, 2.0], 1.0);
    d.faces.x_min = FaceCondition::Velocity { velocity: [u, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0 };
    for f in [&mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.porous = vec![plate(18.0, 22.0, resistance)];
    let m = model(1.0, NU, d, collide(6000.0, 1.0, None, collision));
    let (setup, mut flow) = prepare(&m, &[]);
    assert_eq!(setup.zones[0].layers(), 4);
    run_all(&setup, &mut flow);
    let s = flow.zone_samples()[0];
    (s.pressure, s.flow)
}

#[test]
fn a_porous_plate_drops_the_pressure_by_the_darcy_forchheimer_law() {
    // Both collisions: the resistance enters through Guo's force term, which each re-derives.
    for collision in COLLISIONS {
        let (k, c2, t) = (5.0, 2.0, 4.0);
        for u in [0.01, 0.03] {
            let (dp, q) = plate_drop(u, Resistance::DarcyForchheimer { permeability: k, inertial: c2 }, collision);
            assert!((q - 4.0 * u).abs() < 0.01 * 4.0 * u, "u {u}: flow {q}, expected {}", 4.0 * u);
            // The law at the superficial velocity and the density in the plate: the weakly
            // compressible scheme raises the density by p / c_s² (1 + 1.5 Δp halfway through the
            // plate, the outlet at 0), which moves both by about that fraction.
            let (up, rho) = (q / 4.0, 1.0 + 1.5 * dp);
            let expected = rho * (NU / k * up * t + 0.5 * c2 * up * up * t);
            assert!((dp - expected).abs() < 0.005 * expected, "u {u}: Δp {dp}, Darcy–Forchheimer {expected}");
        }
        // A loss coefficient: Δp = k ρ u²/2 across the plate, whatever its thickness.
        let (coefficient, u) = (3.0, 0.03);
        let (dp, q) = plate_drop(u, Resistance::LossCoefficient { coefficient }, collision);
        let expected = (1.0 + 1.5 * dp) * 0.5 * coefficient * (q / 4.0).powi(2);
        assert!((dp - expected).abs() < 0.005 * expected, "Δp {dp}, k ρ u²/2 {expected}");
    }
}

/// A periodic 40 × 2 × 2 loop: a fan over x 5..7 blowing +x with a linear curve, and the
/// Darcy–Forchheimer plate of the test above.
fn fan_loop(collision: Collision) -> fieldcraft_model::Model {
    let mut d = domain([0.0; 3], [40.0, 2.0, 2.0], 1.0);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.porous = vec![plate(18.0, 22.0, Resistance::DarcyForchheimer { permeability: 5.0, inertial: 2.0 })];
    d.fans = vec![Fan {
        name: "Fan".into(),
        min: [5.0, 0.0, 0.0],
        max: [7.0, 2.0, 2.0],
        direction: FanDirection::XPlus,
        curve: vec![[0.0, 0.004], [0.1, 0.0025], [0.2, 0.0]],
        relaxation: Some(0.05),
    }];
    model(1.0, NU, d, collide(5000.0, 1.0, None, collision))
}

#[test]
fn a_fan_settles_at_the_intersection_of_its_curve_and_the_system_resistance() {
    for collision in COLLISIONS {
        let m = fan_loop(collision);
        let (setup, mut flow) = prepare(&m, &[]);
        let mass0 = flow.mass();
        run_all(&setup, &mut flow);
        // Operating point: the curve's first segment Δp = 0.004 − 0.015 Q meets the plate's
        // Δp = (ν/K) u t + (C₂/2) u² t with u = Q / 4 (area 2 × 2), t = 4.
        let (a, b, c) = (2.0 / 16.0 * 4.0 / 2.0, NU / 5.0 * 4.0 / 4.0 + 0.015, -0.004);
        let q = (-b + (b * b - 4.0 * a * c).sqrt()) / (2.0 * a);
        let rise = 0.004 - 0.015 * q;
        assert!(q < 0.1, "the operating point {q} is on the first segment");
        let samples = flow.zone_samples();
        let (vent, fan) = (samples[0], samples[1]);
        assert!((fan.flow - q).abs() < 0.005 * q, "fan flow {}, operating point {q}", fan.flow);
        assert!((fan.fan_rise - rise).abs() < 0.01 * rise, "fan rise {}, operating point {rise}", fan.fan_rise);
        assert!((vent.flow - q).abs() < 0.005 * q, "vent flow {}", vent.flow);
        assert!((vent.pressure - rise).abs() < 0.005 * rise, "the plate drops {} of the fan's {rise}", vent.pressure);
        assert!((fan.pressure - rise).abs() < 0.005 * rise, "measured rise across the fan {}", fan.pressure);
        assert!(((flow.mass() - mass0) / mass0).abs() < 1e-12, "mass {} → {}", mass0, flow.mass());

        // The solve reports the same through its zone histories, in model units.
        let r = solve_fluid(&m, None, &[], &Limits::default()).unwrap();
        let z = r.zones.iter().find(|z| z.fan).unwrap();
        assert!((z.flow.last().unwrap() - q).abs() < 0.01 * q);
        assert!((z.fan_rise.last().unwrap() - rise).abs() < 0.01 * rise);
        assert_eq!(r.zones.len(), 2);
        assert_eq!(z.flow.len(), r.history_times.len());
    }
}

#[test]
fn zones_off_the_grid_or_on_an_open_layer_are_refused() {
    let base = || {
        let mut d = domain([0.0; 3], [10.0, 2.0, 2.0], 1.0);
        d.faces.x_min = FaceCondition::Velocity { velocity: [0.01, 0.0, 0.0], profile: Profile::Uniform };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0 };
        d
    };
    let resistance = Resistance::LossCoefficient { coefficient: 1.0 };
    for zone in [
        plate(20.0, 30.0, resistance), // outside the box
        plate(3.1, 3.4, resistance),   // thinner than a cell, no centre in it
        plate(0.0, 2.0, resistance),   // on the inlet layer
        plate(9.0, 10.0, resistance),  // on the outlet layer
    ] {
        let mut d = base();
        d.porous = vec![zone.clone()];
        let m = model(1.0, NU, d, step(10.0, 1.0, None));
        let made = Setup::new(&m, None, &[], &Limits::default()).and_then(|s| Flow::new(&s, &[]));
        let e = made.err().unwrap_or_else(|| panic!("{zone:?} was accepted")).to_string();
        assert!(e.contains("Vent"), "{e}");
    }
}
