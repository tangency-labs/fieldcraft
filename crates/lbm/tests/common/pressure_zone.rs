//! A pressure face whose inner fluid cell belongs to a force zone.
use super::*;
use fieldcraft_model::{Axis, FaceCondition, Fan, FanDirection, Outflow, PorousZone, Profile, Resistance};

pub fn neighbor(collision: Collision, fan: bool) -> (Setup, Flow) {
    let mut d = domain([0.0; 3], [8.0, 2.0, 2.0], 1.0);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.03, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    for f in [&mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.initial_velocity = [0.03, 0.0, 0.0];
    if fan {
        d.fans.push(Fan {
            name: "Outlet neighbor".into(),
            min: [6.0, 0.0, 0.0],
            max: [7.0, 2.0, 2.0],
            direction: FanDirection::XPlus,
            curve: vec![[0.0, 0.001], [1.0, 0.001]],
            relaxation: Some(1.0),
        });
    } else {
        d.porous.push(PorousZone {
            name: "Outlet neighbor".into(),
            min: [6.0, 0.0, 0.0],
            max: [7.0, 2.0, 2.0],
            axis: Axis::X,
            resistance: Resistance::DarcyForchheimer { permeability: 5.0, inertial: 2.0 },
        });
    }
    prepare(&model(1.0, 0.1, d, collide(100.0, 1.0, None, collision)), &[])
}
