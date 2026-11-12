//! Open-face integration against independently summed, nonuniform fields.
mod common;

use common::*;
use fieldcraft_lbm::{FluidResult, HeatHistory, Limits, PortProperties, ports::inlet_temperature, run};
use fieldcraft_model::{FaceCondition, FluidThermal, Outflow, Profile, ThermalCondition};

fn saved() -> FluidResult {
    let mut d = domain([0.0; 3], [2.0; 3], 1.0);
    d.thermal = Some(Default::default());
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.1, 0.0, 0.0], profile: Default::default() };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Default::default(), sigma: None };
    let m = hot_model(2.0, 0.1, 0.1, d, step(0.1, 0.1, None));
    let (setup, mut flow) = prepare(&m, &[]);
    let mut r = run(&setup, &mut flow, &Limits::default()).unwrap();
    r.port_properties = Some(PortProperties { density: 2.0, specific_heat: Some(5.0) });
    r.temperature_reference = 273.15;
    let f = r.frames.last_mut().unwrap();
    f.velocity.fill([0.0; 3]);
    f.pressure.fill(10.0);
    f.temperature.fill(20.0);
    r
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-10, "{a} vs {b}");
}

#[test]
fn backflow_is_separate_and_temperatures_are_flow_weighted() {
    let mut r = saved();
    let f = r.frames.last_mut().unwrap();
    // x_max cells 1, 3, 5, 7: two outgoing, one incoming, one blocked.
    f.velocity[1][0] = 1.0;
    f.velocity[3][0] = 3.0;
    f.velocity[5][0] = -2.0;
    f.velocity[7][0] = 1000.0;
    f.temperature[1] = 10.0;
    f.temperature[3] = 30.0;
    f.temperature[5] = 40.0;
    f.pressure[1] = 4.0;
    f.pressure[3] = 8.0;
    r.solid[7] = true;
    let ports = r.ports(r.frames.len() - 1).unwrap();
    let p = &ports[1];
    close(p.area, 3.0);
    let mean = p.mean.unwrap();
    close(mean.pressure, 22.0 / 3.0);
    close(mean.temperature.unwrap(), 273.15 + 80.0 / 3.0);
    close(mean.velocity[0], 2.0 / 3.0);
    close(p.outflow.volume_flow, 4.0);
    close(p.inflow.volume_flow, 2.0);
    close(p.outflow.mass_flow.unwrap(), 8.0);
    close(p.inflow.mass_flow.unwrap(), 4.0);
    close(p.outflow.temperature.unwrap(), 298.15); // (1*10 + 3*30)/4, not (10+30)/2
    close(p.inflow.temperature.unwrap(), 313.15);
    close(p.outflow.pressure.unwrap(), 7.0);
    // Heat relative to the mixed inflow temperature, the backflow's 40 (x_min has no flow):
    // ρ c Σ Q (T − 40) = 10 (1 (10 − 40) + 3 (30 − 40)) out, 0 in.
    close(inlet_temperature(&ports).unwrap(), 313.15);
    close(p.outflow.sensible_heat.unwrap(), -600.0);
    close(p.inflow.sensible_heat.unwrap(), 0.0);
    // Reading an earlier frame must not reuse the final field.
    assert_ne!(r.ports(0).unwrap(), ports);
}

#[test]
fn every_axis_uses_the_outward_normal_and_cell_area() {
    let mut r = saved();
    r.grid.spacing = 0.5;
    r.open_faces = [true; 6];
    r.frames.last_mut().unwrap().velocity.fill([1.0, 2.0, 3.0]);
    for p in r.ports(r.frames.len() - 1).unwrap() {
        close(p.area, 1.0);
        let (stream, rest) = if p.face.is_max() { (p.outflow, p.inflow) } else { (p.inflow, p.outflow) };
        close(stream.volume_flow, (p.face.axis() + 1) as f64);
        close(rest.volume_flow, 0.0);
        assert_eq!(rest.temperature, None);
        assert_eq!(rest.pressure, None);
    }
}

#[test]
fn absent_properties_are_not_invented_and_isothermal_results_do_not_invent_heat() {
    let mut r = saved();
    r.frames.last_mut().unwrap().velocity.fill([1.0, 0.0, 0.0]);
    r.port_properties = None;
    let p = &r.ports(r.frames.len() - 1).unwrap()[1];
    close(p.outflow.volume_flow, 4.0);
    assert!(p.outflow.temperature.is_some());
    assert_eq!(p.outflow.mass_flow, None);
    assert_eq!(p.outflow.sensible_heat, None);
    r.frames.last_mut().unwrap().temperature.clear();
    r.port_properties = Some(PortProperties { density: 2.0, specific_heat: None });
    let p = &r.ports(r.frames.len() - 1).unwrap()[1];
    assert_eq!(p.outflow.temperature, None);
    assert_eq!(p.mean.unwrap().temperature, None);
    assert_eq!(p.outflow.sensible_heat, None);
    close(p.outflow.mass_flow.unwrap(), 8.0);
}

#[test]
fn malformed_saved_fields_are_refused() {
    let r = saved();
    assert!(r.ports(usize::MAX).is_err());
    for kind in 0..9 {
        let mut bad = r.clone();
        match kind {
            0 => bad.grid.cells = [usize::MAX; 3],
            1 => bad.grid.cells[0] = 0,
            2 => bad.grid.spacing = f64::INFINITY,
            3 => {
                bad.frames.last_mut().unwrap().temperature.pop();
            }
            4 => bad.frames.last_mut().unwrap().velocity[1][0] = f32::NAN,
            5 => bad.port_properties.as_mut().unwrap().density = -1.0,
            6 => bad.port_properties.as_mut().unwrap().specific_heat = Some(f64::NAN),
            7 => bad.grid.spacing = f64::MAX,
            _ => bad.frames.last_mut().unwrap().velocity[1][1] = f32::NAN,
        }
        assert!(bad.ports(bad.frames.len() - 1).is_err(), "kind {kind}");
    }
}

#[test]
fn stagnant_faces_have_area_means_but_blocked_faces_do_not() {
    let mut r = saved();
    let p = &r.ports(r.frames.len() - 1).unwrap()[0];
    let mean = p.mean.unwrap();
    close(mean.pressure, 10.0);
    close(mean.temperature.unwrap(), 293.15);
    assert_eq!(mean.velocity, [0.0; 3]);
    assert_eq!(p.inflow.temperature, None);
    r.solid.fill(true);
    for p in r.ports(r.frames.len() - 1).unwrap() {
        assert_eq!(p.area, 0.0);
        assert_eq!(p.mean, None);
    }
}

#[test]
fn wall_power_interpolates_between_recorded_samples_without_extrapolating() {
    let mut r = saved();
    let at = |r: &FluidResult, power: &[f64], time| r.wall_power(&HeatHistory { name: "wall".into(), body: None, power: power.to_vec() }, time);
    r.history_times = vec![1.0, 3.0];
    assert_eq!(at(&r, &[-4.0, 8.0], 2.0), Some(2.0));
    assert_eq!(at(&r, &[-4.0, 8.0], 1.0), Some(-4.0));
    assert_eq!(at(&r, &[-4.0, 8.0], 0.0), None);
    assert_eq!(at(&r, &[-4.0, 8.0], 4.0), None);
    assert_eq!(at(&r, &[], 1.0), None);
    r.history_times = vec![1.0, 1.0];
    assert_eq!(at(&r, &[0.0, 1.0], 1.0), None);
    r.history_times = vec![1.0];
    assert_eq!(at(&r, &[f64::NAN], 1.0), None);
}

#[test]
fn solve_records_the_fluid_properties_without_display_unit_conversions() {
    let mut d = domain([0.0; 3], [2.0; 3], 1.0);
    d.thermal = Some(Default::default());
    d.initial_velocity = [0.0; 3];
    let m = hot_model(7.0, 0.1, 0.1, d, step(0.1, 0.1, None));
    let (setup, mut flow) = prepare(&m, &[]);
    let r = run(&setup, &mut flow, &Limits::default()).unwrap();
    let props = r.port_properties.unwrap();
    close(props.density, 7.0);
    close(props.specific_heat.unwrap(), 1.0);
    assert!(r.ports(0).unwrap().is_empty());
}

/// The heat carried through the ports does not depend on the lattice's temperature zero, also
/// while the net volume flow is not zero: the same field written above another zero gives the
/// same heat (#704 review).
#[test]
fn port_heat_does_not_depend_on_the_temperature_zero() {
    let mut r = saved();
    let f = r.frames.last_mut().unwrap();
    // In through x_min at 20 and out through x_max at 25, twice as fast: a transient.
    for (c, v) in f.velocity.iter_mut().enumerate() {
        v[0] = if c % 2 == 0 { 1.0 } else { 2.0 };
    }
    for (c, t) in f.temperature.iter_mut().enumerate() {
        *t = if c % 2 == 0 { 20.0 } else { 25.0 };
    }
    let heat = |r: &FluidResult| {
        let ports = r.ports(r.frames.len() - 1).unwrap();
        ports.iter().map(|p| p.outflow.sensible_heat.unwrap() - p.inflow.sensible_heat.unwrap()).sum::<f64>()
    };
    let before = heat(&r);
    // ρ c Q_out (T_out − T_in) = 10 · 8 · 5.
    close(before, 400.0);
    let mut shifted = r.clone();
    shifted.temperature_reference -= 100.0;
    for t in &mut shifted.frames.last_mut().unwrap().temperature {
        *t += 100.0;
    }
    close(heat(&shifted), before);
}

/// A solved flow, end to end: a uniform inflow at T_in = 0 (periodic across, so the profile
/// stays uniform away from the block) past a conducting block that generates P, run to a steady
/// state (Pe 50 on the box's height). The inlet face carries the prescribed U A, the outlet
/// carries what enters, and the heat the ports carry out is the block's power: the plume is
/// downstream, so next to nothing conducts back into the held inlet. The lattice is weakly
/// compressible: the volume flow grows by the density the block's pressure drop takes away
/// (Δρ/ρ = Δp / (ρ c_s²)), 10 % at lattice speed 0.05 here, so the run is at 0.01 (0.4 %).
#[test]
fn a_solved_heated_flow_s_ports_carry_its_inflow_and_its_heat() {
    let (n, h, l, u, peclet, power) = (8usize, 1.0, 4.0, 0.1, 50.0, 0.01);
    let dx = h / n as f64;
    let alpha = u * h / peclet;
    let mut d = domain([0.0; 3], [l, h, dx], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [u, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
    for f in [&mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.obstacles.push("Block".into());
    let mut t = FluidThermal { initial_temperature: kelvin(0.0), ..FluidThermal::default() };
    t.faces.x_min = ThermalCondition::Temperature { value: kelvin(0.0) };
    t.solids.push(conducting("Block", power));
    d.thermal = Some(t);
    let block = obstacle("Block", cuboid([1.5, 0.375, -1.0], [2.0, 0.625, 2.0]), &CUBOID_FACES);
    let (rho, cp) = (1.0, 1.0);
    // Six flow-through times at lattice speed 0.01.
    let m = cht_model([rho, 0.05, alpha], (1.0, [0.5; 3], 1.0), d, step(6.0 * l / u, 0.01 * dx / u, None));
    let setup = fieldcraft_lbm::Setup::new(&m, None, std::slice::from_ref(&block), &Limits::default()).unwrap();
    let mut flow = fieldcraft_lbm::Flow::new(&setup, std::slice::from_ref(&block)).unwrap();
    let r = run(&setup, &mut flow, &Limits::default()).unwrap();
    let ports = r.ports(r.frames.len() - 1).unwrap();
    let (inlet, outlet) = (&ports[0], &ports[1]);
    assert_eq!((inlet.face.name(), outlet.face.name()), ("x_min", "x_max"));
    let expected = u * h * dx;
    let q_in = inlet.inflow.volume_flow;
    let carried: f64 = ports.iter().map(|p| p.outflow.sensible_heat.unwrap() - p.inflow.sensible_heat.unwrap()).sum();
    println!(
        "heated block: inflow {q_in:.6e} (U A {expected:.6e}), outflow {:.6e}; heat carried {carried:.6e} of the block's {power:e}; outlet {:.4} K over the inlet",
        outlet.outflow.volume_flow,
        outlet.outflow.temperature.unwrap() - inlet.inflow.temperature.unwrap()
    );
    assert!((q_in / expected - 1.0).abs() < 1e-6, "{q_in} vs {expected}");
    assert!(inlet.outflow.volume_flow == 0.0 && outlet.inflow.volume_flow == 0.0, "no backflow");
    assert!((outlet.outflow.volume_flow / q_in - 1.0).abs() < 0.01, "steady: what enters leaves");
    assert!((inlet_temperature(&ports).unwrap() - 273.15).abs() < 1e-4, "the inlet's held temperature");
    // The mixed outlet temperature is the power over ρ c Q.
    let rise = outlet.outflow.temperature.unwrap() - 273.15;
    assert!((rise * rho * cp * q_in / power - 1.0).abs() < 0.02, "outlet rise {rise}");
    assert!((carried / power - 1.0).abs() < 0.02, "{carried} carried vs {power} generated");
}
