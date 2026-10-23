//! GPU zone estimates include immutable curve/header upload, sampled reductions, and CPU cell
//! ownership/lists while preparing. They agree before and after preparation.
mod common;
use common::*;
use fieldcraft_lbm::{Choice, Memory};
use fieldcraft_model::{BoxFace, FaceCondition, Fan, FanDirection};

#[test]
fn zones_and_curve_points_increase_the_memory_estimate() {
    let mut d = domain([0.0; 3], [8.0, 2.0, 2.0], 1.0);
    for f in BoxFace::ALL {
        *d.faces.get_mut(f) = FaceCondition::Periodic {};
    }
    let no_zones = model(1.0, 0.1, d.clone(), step(2.0, 1.0, None));
    let base = Memory::estimate(&no_zones, None, Choice::Device).unwrap();
    d.fans = vec![Fan {
        name: "Fan".into(),
        min: [2.0, 0.0, 0.0],
        max: [4.0, 2.0, 2.0],
        direction: FanDirection::XPlus,
        curve: vec![[0.0, 0.002], [0.1, 0.0]],
        relaxation: None,
    }];
    let m = model(1.0, 0.1, d.clone(), step(2.0, 1.0, None));
    let estimate = Memory::estimate(&m, None, Choice::Device).unwrap();
    let (setup, _) = prepare(&m, &[]);
    assert_eq!(estimate, Memory::of(&setup, Choice::Device).unwrap());
    // 32 cells * 2 u32 ownership entries; 16 header words + 4 curve words, staged 3 times.
    assert_eq!(estimate.setup - base.setup, 32 * 8 + 80 * 3);
    assert_eq!(estimate.state - base.state, 80);
    assert_eq!(estimate.readback - base.readback, 16);
    d.fans[0].curve.insert(1, [0.05, 0.001]);
    let m = model(1.0, 0.1, d, step(2.0, 1.0, None));
    let extra = Memory::estimate(&m, None, Choice::Device).unwrap();
    assert_eq!(extra.setup - estimate.setup, 3 * 8);
    assert_eq!(extra.state - estimate.state, 8);
}
