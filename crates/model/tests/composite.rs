use fieldcraft_model::{Layup, MaterialOrientation, Ply, PlyStrength};

#[test]
fn material_frames_are_right_handed_and_shell_directions_are_projected() {
    let rectangular = MaterialOrientation::Rectangular { axis1: [1.0, 0.0, 1.0], axis2: [0.0, 2.0, 0.0] };
    let shell = rectangular.shell_axes_at([0.0; 3], [0.0, 0.0, 1.0], 1).unwrap();
    assert_eq!(shell, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    for scale in [1e-300, 1.0, 1e300] {
        let r = MaterialOrientation::Rectangular { axis1: [scale, 0.0, 0.0], axis2: [scale, scale, 0.0] }.axes_at([0.0; 3]).unwrap();
        assert_eq!(r, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    }
    let cylinder = MaterialOrientation::Cylindrical { origin: [1.0, 2.0, 3.0], axis: [0.0, 0.0, 1.0] };
    assert_eq!(cylinder.axes_at([1.0, 4.0, 5.0]).unwrap(), [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
    assert!(cylinder.axes_at([1.0, 2.0, 5.0]).is_err());
    assert!(cylinder.shell_axes_at([1.0, 4.0, 5.0], [0.0, 1.0, 0.0], 1).is_err());
    assert!(cylinder.shell_axes_at([1.0, 4.0, 5.0], [0.0, 1.0, 0.0], 3).is_ok());
    assert!(cylinder.shell_axes_at([1.0, 4.0, 5.0], [0.0, 1.0, 0.0], 0).is_err());
    for bad in [[0.0; 3], [f64::NAN, 0.0, 0.0], [f64::INFINITY, 0.0, 0.0]] {
        assert!(MaterialOrientation::Rectangular { axis1: bad, axis2: [0.0, 1.0, 0.0] }.validate().is_err());
    }
    assert!(MaterialOrientation::Rectangular { axis1: [1.0, 0.0, 0.0], axis2: [2.0, 0.0, 0.0] }.validate().is_err());
}

#[test]
fn nearly_parallel_frames_remain_orthogonal_and_cylinder_origins_can_move_along_the_axis() {
    let frame = MaterialOrientation::Rectangular { axis1: [1.0; 3], axis2: [1.0 + 1e-11, 1.0, 1.0] }.axes_at([0.0; 3]).unwrap();
    for i in 0..3 {
        for j in 0..3 {
            let dot: f64 = frame[i].iter().zip(frame[j]).map(|(a, b)| a * b).sum();
            assert!((dot - if i == j { 1.0 } else { 0.0 }).abs() < 1e-14);
        }
    }
    let a = MaterialOrientation::Cylindrical { origin: [0.0; 3], axis: [0.0, 0.0, 1.0] };
    let b = MaterialOrientation::Cylindrical { origin: [0.0, 0.0, 1e8], axis: [0.0, 0.0, 1.0] };
    assert_eq!(a.axes_at([1e-6, 0.0, 1e8]).unwrap(), b.axes_at([1e-6, 0.0, 1e8]).unwrap());
}

#[test]
fn layup_and_strength_json_is_explicit_and_bounded() {
    let layup: Layup = serde_json::from_str(r#"{"plies":[{"material":"Carbon","thickness":0.001}]}"#).unwrap();
    layup.validate().unwrap();
    assert_eq!(layup.plies[0].points, 3);
    assert_eq!(layup.thickness(), 0.001);
    assert_eq!(serde_json::from_str::<Layup>(&serde_json::to_string(&layup).unwrap()).unwrap(), layup);
    let p = layup.plies[0].clone();
    for plies in [
        vec![],
        vec![p.clone(); 257],
        vec![Ply { points: 0, ..p.clone() }],
        vec![Ply { thickness: f64::NAN, ..p.clone() }],
        vec![Ply { angle: f64::INFINITY, ..p.clone() }],
    ] {
        assert!(Layup { plies, offset: 0.0 }.validate().is_err());
    }
    assert!(Layup { offset: f64::NAN, ..layup }.validate().is_err());
    let s: PlyStrength = serde_json::from_str(r#"{"tension":[1,2],"compression":[3,4],"shear12":5,"shear23":6}"#).unwrap();
    assert_eq!(s.tsai_wu_interaction, -0.5);
    s.validate().unwrap();
    assert!(PlyStrength { tsai_wu_interaction: 1.0, ..s }.validate().is_err());
    assert!(PlyStrength { shear23: 0.0, ..s }.validate().is_err());
    assert!(serde_json::from_str::<PlyStrength>(r#"{"tension":[1,2],"compression":[3,4],"shear12":5}"#).is_err());
    assert!(serde_json::from_str::<Ply>(r#"{"material":"Carbon","thickness":1,"typo":2}"#).is_err());
}
