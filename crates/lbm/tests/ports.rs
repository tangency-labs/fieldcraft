//! Open-face integration against independently summed, nonuniform fields.
mod common;

use common::*;
use fieldcraft_lbm::{FluidResult, Limits, PortProperties, run};
use fieldcraft_model::FaceCondition;

fn saved() -> FluidResult {
    let mut d = domain([0.0; 3], [2.0; 3], 1.0);
    d.thermal = Some(Default::default());
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.1, 0.0, 0.0], profile: Default::default() };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Default::default() };
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
    close(p.outflow.volume_flow, 4.0);
    close(p.inflow.volume_flow, 2.0);
    close(p.outflow.mass_flow.unwrap(), 8.0);
    close(p.inflow.mass_flow.unwrap(), 4.0);
    close(p.outflow.temperature.unwrap(), 298.15); // (1*10 + 3*30)/4, not (10+30)/2
    close(p.inflow.temperature.unwrap(), 313.15);
    close(p.outflow.pressure.unwrap(), 7.0);
    close(p.outflow.sensible_heat.unwrap(), 1000.0);
    close(p.inflow.sensible_heat.unwrap(), 800.0);
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
fn old_results_do_not_invent_properties_and_isothermal_results_do_not_invent_heat() {
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
    assert_eq!(p.outflow.sensible_heat, None);
    close(p.outflow.mass_flow.unwrap(), 8.0);
}

#[test]
fn malformed_saved_fields_are_refused() {
    let r = saved();
    assert!(r.ports(usize::MAX).is_err());
    for kind in 0..8 {
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
            _ => bad.grid.spacing = f64::MAX,
        }
        assert!(bad.ports(bad.frames.len() - 1).is_err(), "kind {kind}");
    }
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
