//! Upload checks that need no GPU device.
#[path = "../../lbm/tests/common/mod.rs"]
#[allow(clippy::unwrap_used)]
mod common;
use super::zones;
use common::*;
use fieldcraft_model::{BoxFace, FaceCondition, Fan, FanDirection};

#[test]
fn curves_that_collapse_in_single_precision_are_refused_before_device_work() {
    let mut d = domain([0.0; 3], [8.0, 2.0, 2.0], 1.0);
    for face in BoxFace::ALL {
        *d.faces.get_mut(face) = FaceCondition::Periodic {};
    }
    d.fans = vec![Fan {
        name: "Indistinguishable".into(),
        min: [2.0, 0.0, 0.0],
        max: [4.0, 2.0, 2.0],
        direction: FanDirection::XPlus,
        curve: vec![[0.0, 0.002], [0.1, 0.001], [0.1 + 1e-10, 0.0]],
        relaxation: None,
    }];
    let m = model(1.0, 0.1, d, step(1.0, 1.0, None));
    let (_, flow) = prepare(&m, &[]);
    let mut flags = vec![1; flow.grid().len()];
    let e = zones::append(&flow, &mut flags).unwrap_err().to_string();
    assert!(e.contains("Indistinguishable") && e.contains("single precision"), "{e}");
}

#[test]
fn unzoned_flags_are_unchanged_and_curve_offsets_are_absolute() {
    let mut d = domain([0.0; 3], [8.0, 2.0, 2.0], 1.0);
    for face in BoxFace::ALL {
        *d.faces.get_mut(face) = FaceCondition::Periodic {};
    }
    let m = model(1.0, 0.1, d.clone(), step(1.0, 1.0, None));
    let (_, flow) = prepare(&m, &[]);
    let mut flags = vec![1; 32];
    zones::append(&flow, &mut flags).unwrap();
    assert_eq!(flags, vec![1; 32]);
    d.fans = vec![Fan {
        name: "Fan".into(),
        min: [2.0, 0.0, 0.0],
        max: [4.0, 2.0, 2.0],
        direction: FanDirection::XPlus,
        curve: vec![[0.0, 0.002], [0.1, 0.0]],
        relaxation: None,
    }];
    let m = model(1.0, 0.1, d, step(1.0, 1.0, None));
    let (_, flow) = prepare(&m, &[]);
    zones::append(&flow, &mut flags).unwrap();
    assert_eq!(flags.len(), 32 + 16 + 4);
    let offset = flags[32 + 13] as usize;
    assert_eq!(offset, 48);
    assert_eq!(f32::from_bits(flags[offset]), 0.0);
    assert_eq!(f32::from_bits(flags[offset + 1]), 0.002f32);
    assert_eq!(flags[..32].iter().filter(|flag| **flag == 2).count(), 8);
}
