//! Fan/porous parity, timing and analytical operating points. CPU cases: lbm/tests/zones.rs.
use super::*;
use fieldcraft_model::{Axis, Fan, FanDirection, PorousZone, Resistance};

fn vent(min: [f64; 3], max: [f64; 3], resistance: Resistance) -> PorousZone {
    PorousZone { name: "Vent".into(), min, max, axis: Axis::X, resistance }
}

fn periodic(d: &mut FluidDomain) {
    for face in BoxFace::ALL {
        *d.faces.get_mut(face) = PERIODIC;
    }
}

fn fan_loop(collision: Collision, steps: u64) -> (Setup, Flow) {
    let mut d = domain([0.0; 3], [40.0, 2.0, 2.0], 1.0);
    periodic(&mut d);
    d.porous = vec![vent([18.0, 0.0, 0.0], [22.0, 2.0, 2.0], Resistance::DarcyForchheimer { permeability: 5.0, inertial: 2.0 })];
    d.fans = vec![Fan {
        name: "Fan".into(),
        min: [5.0, 0.0, 0.0],
        max: [7.0, 2.0, 2.0],
        direction: FanDirection::XPlus,
        curve: vec![[0.0, 0.004], [0.1, 0.0025], [0.2, 0.0]],
        relaxation: Some(0.05),
    }];
    let m = model(1.0, 0.1, d, collide(steps as f64, 1.0, None, collision));
    prepare(&m, &[])
}

fn samples_match(cpu: &mut Flow, gpu: &mut GpuFlow) {
    let a = cpu.zone_samples();
    let b = gpu.zones().unwrap();
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(&b) {
        for (name, x, y) in [("flow", a.flow, b.flow), ("pressure", a.pressure, b.pressure), ("fan rise", a.fan_rise, b.fan_rise)] {
            let err = (x - y).abs() / x.abs().max(1e-8);
            println!("{name}: CPU {x:.9e}, GPU {y:.9e}, relative {err:.2e}");
            assert!(err < 2e-4, "{name}: {err}");
        }
    }
}

#[test]
fn fan_feedback_runs_every_step_and_is_independent_of_chunks_and_readbacks() {
    for (collision, streaming) in cases() {
        let (_, mut cpu) = fan_loop(collision, 129);
        let Some(mut device) = gpu(&cpu, streaming) else { continue };
        // First three updates make an off-by-one relaxation error readily visible.
        for k in [1, 1, 1, 14, 46, 1, 1, 64] {
            cpu.advance(k).unwrap();
            device.advance(k).unwrap();
            samples_match(&mut cpu, &mut device);
        }
        let (du, drho) = field_difference(&mut cpu, &mut device);
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
        let expected = device.fields().unwrap();
        let samples = device.zones().unwrap();
        drop(device); // One GPU flow at a time, including within this test.
        let (_, fresh) = fan_loop(collision, 129);
        let Some(mut single) = gpu(&fresh, streaming) else { continue };
        single.advance(129).unwrap();
        assert_eq!(single.fields().unwrap(), expected, "{collision:?} {streaming:?}: chunking changed fields");
        assert_eq!(single.zones().unwrap(), samples, "chunking changed zone diagnostics");
    }
}

#[test]
fn strong_porous_resistance_has_the_implicit_analytic_velocity() {
    for (collision, streaming) in cases() {
        let mut d = domain([0.0; 3], [2.0; 3], 1.0);
        periodic(&mut d);
        d.acceleration = [0.001, 0.0002, 0.0];
        d.initial_velocity = [0.01, -0.005, 0.0];
        d.porous = vec![vent([0.0; 3], [2.0; 3], Resistance::DarcyForchheimer { permeability: 0.005, inertial: 200.0 })];
        let m = model(1.0, 0.1, d, collide(96.0, 1.0, None, collision));
        let (setup, flow) = prepare(&m, &[]);
        let Some((mut cpu, mut device)) = both(&setup, flow, streaming) else { continue };
        let magnitude = (0.001f64.powi(2) + 0.0002f64.powi(2)).sqrt();
        let speed = 2.0 * magnitude / (20.0 + (400.0 + 400.0 * magnitude).sqrt());
        let velocity = device.fields().unwrap().velocity[0];
        let expected = [speed * 0.001 / magnitude, speed * 0.0002 / magnitude, 0.0];
        for (got, expected) in velocity.into_iter().zip(expected) {
            assert!((got - expected).abs() < 2e-5 * speed, "{got} vs {expected}");
        }
        let (du, drho) = field_difference(&mut cpu, &mut device);
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
    }
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn porous_pressure_drop_and_fan_operating_point_match_the_analytic_laws() {
    for (collision, streaming) in cases() {
        for resistance in [Resistance::DarcyForchheimer { permeability: 5.0, inertial: 2.0 }, Resistance::LossCoefficient { coefficient: 3.0 }] {
            let mut d = domain([0.0; 3], [40.0, 2.0, 2.0], 1.0);
            periodic(&mut d);
            d.faces.x_min = FaceCondition::Velocity { velocity: [0.03, 0.0, 0.0], profile: Profile::Uniform };
            d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
            d.porous = vec![vent([18.0, 0.0, 0.0], [22.0, 2.0, 2.0], resistance)];
            let m = model(1.0, 0.1, d, collide(6000.0, 1.0, None, collision));
            let (setup, flow) = prepare(&m, &[]);
            let Some((mut cpu, mut device)) = both(&setup, flow, streaming) else { continue };
            let sample = device.zones().unwrap()[0];
            let (u, rho) = (sample.flow / 4.0, 1.0 + 1.5 * sample.pressure);
            let exact = match resistance {
                Resistance::DarcyForchheimer { .. } => rho * (0.1 / 5.0 * u * 4.0 + u * u * 4.0),
                Resistance::LossCoefficient { .. } => 0.5 * 3.0 * rho * u * u,
            };
            println!("{collision:?} {streaming:?}: porous measured {}, analytic {exact}", sample.pressure);
            assert!((sample.pressure - exact).abs() < 0.005 * exact);
            samples_match(&mut cpu, &mut device);
        }
        let (setup, flow) = fan_loop(collision, 5000);
        let Some((mut cpu, mut device)) = both(&setup, flow, streaming) else { continue };
        let q = (-0.035f64 + (0.035f64.powi(2) + 4.0 * 0.25 * 0.004).sqrt()) / 0.5;
        let rise = 0.004 - 0.015 * q;
        let samples = device.zones().unwrap();
        println!("{collision:?} {streaming:?}: operating point Q={q}, dp={rise}; measured {samples:?}");
        assert!((samples[1].flow - q).abs() < 0.005 * q);
        assert!((samples[1].fan_rise - rise).abs() < 0.01 * rise);
        assert!((samples[0].pressure - rise).abs() < 0.005 * rise);
        samples_match(&mut cpu, &mut device);
        let mass = device.mass().unwrap();
        assert!((mass - 160.0).abs() < 1e-6 * 160.0, "{mass}");
    }
}

#[test]
fn a_powered_conductor_with_fan_porous_drag_and_buoyancy_matches_cpu() {
    for (collision, streaming) in cases() {
        let mut d = domain([0.0; 3], [16.0, 6.0, 2.0], 1.0);
        periodic(&mut d);
        d.faces.y_min = FaceCondition::Wall {};
        d.faces.y_max = FaceCondition::Wall {};
        d.obstacles = vec!["Plate".into()];
        d.porous = vec![vent([8.0, 1.0, 0.0], [10.0, 6.0, 2.0], Resistance::LossCoefficient { coefficient: 3.0 })];
        d.fans = vec![Fan {
            name: "Fan".into(),
            min: [2.0, 1.0, 0.0],
            max: [4.0, 6.0, 2.0],
            direction: FanDirection::XPlus,
            curve: vec![[0.0, 0.002], [1.0, 0.0]],
            relaxation: Some(0.05),
        }];
        let mut t = FluidThermal { initial_temperature: kelvin(1.0), gravity: [0.0, -0.001, 0.0], expansion: Some(0.001), ..Default::default() };
        t.faces.y_max = ThermalCondition::Temperature { value: kelvin(0.0) };
        t.solids = vec![conducting("Plate", 0.01)];
        d.thermal = Some(t);
        let m = cht_model([1.0, 0.1, 0.05], (2.0, [0.2; 3], 1.0), d, collide(129.0, 1.0, None, collision));
        let plate = [obstacle("Plate", cuboid([-1.0, -1.0, -1.0], [17.0, 1.0, 3.0]), &CUBOID_FACES)];
        let (setup, flow) = prepare(&m, &plate);
        let Some((mut cpu, mut device)) = both(&setup, flow, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut device);
        let (dt, dh) = heat_difference(&mut cpu, &mut device);
        println!("{collision:?} {streaming:?}: fan/porous/CHT/buoyancy velocity {du:e}, density {drho:e}, temperature {dt:e}, heat {dh:e}");
        assert!(du < 3e-5 && drho < 2e-6, "{du} {drho}");
        assert!(dt < 2e-5 && dh < 5e-5, "{dt} {dh}");
        samples_match(&mut cpu, &mut device);
    }
}

#[test]
fn fan_axes_signs_curve_endpoints_and_overlapping_ownership_match_cpu() {
    for (axis, direction) in [
        (0, FanDirection::XPlus),
        (0, FanDirection::XMinus),
        (1, FanDirection::YPlus),
        (1, FanDirection::YMinus),
        (2, FanDirection::ZPlus),
        (2, FanDirection::ZMinus),
    ] {
        let mut d = domain([0.0; 3], [8.0; 3], 1.0);
        periodic(&mut d);
        d.initial_velocity[axis] = 0.02;
        let (mut lo, mut hi) = ([0.0; 3], [8.0; 3]);
        lo[axis] = 2.0;
        hi[axis] = 6.0;
        let mut p = vent(lo, hi, Resistance::LossCoefficient { coefficient: 1.0 });
        p.axis = [Axis::X, Axis::Y, Axis::Z][axis];
        d.porous = vec![p];
        lo[axis] = 4.0;
        hi[axis] = 8.0;
        d.fans =
            vec![Fan { name: "Fan".into(), min: lo, max: hi, direction, curve: vec![[0.0, 0.003], [1.0, 0.002], [2.0, 0.0]], relaxation: Some(0.1) }];
        let collision = if axis == 1 { Collision::Regularized } else { Collision::Trt };
        let streaming = if axis == 2 { Streaming::TwoBuffers } else { Streaming::InPlace };
        let m = model(1.0, 0.1, d, collide(17.0, 1.0, None, collision));
        let (setup, flow) = prepare(&m, &[]);
        let Some((mut cpu, mut device)) = both(&setup, flow, streaming) else { continue };
        samples_match(&mut cpu, &mut device);
        let (du, drho) = field_difference(&mut cpu, &mut device);
        assert!(du < 2e-5 && drho < 1e-6, "{direction:?} {du} {drho}");
    }
}

#[test]
fn a_warm_cpu_flow_keeps_its_fan_action_when_uploaded() {
    let (_, mut flow) = fan_loop(Collision::Regularized, 17);
    flow.advance(17).unwrap();
    let Some(mut device) = gpu(&flow, Streaming::InPlace) else { return };
    samples_match(&mut flow, &mut device);
    flow.advance(3).unwrap();
    device.advance(3).unwrap();
    samples_match(&mut flow, &mut device);
}

#[test]
fn a_darcy_term_whose_square_overflows_f32_still_matches_cpu() {
    for (collision, streaming) in cases() {
        let mut d = domain([0.0; 3], [2.0; 3], 1.0);
        periodic(&mut d);
        d.initial_velocity = [0.01, 0.0, 0.0];
        d.porous = vec![vent([0.0; 3], [2.0; 3], Resistance::DarcyForchheimer { permeability: 1e-25, inertial: 0.0 })];
        let m = model(1.0, 0.1, d, collide(3.0, 1.0, None, collision));
        let (setup, flow) = prepare(&m, &[]);
        let Some((mut cpu, mut device)) = both(&setup, flow, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut device);
        assert!(du < 2e-5 && drho < 1e-6, "{collision:?} {streaming:?} {du} {drho}");
    }
}
