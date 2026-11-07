use fieldcraft_fem::material::{elasticity, elasticity_in_axes, stress};
use fieldcraft_model::{Elastic, MaterialOrientation};

#[test]
fn rotated_solid_elasticity_preserves_energy_and_maps_material_axes() {
    let material = Elastic::EngineeringConstants { young: [100.0, 10.0, 20.0], poisson: [0.0; 3], shear: [5.0, 4.0, 3.0] };
    let axes = MaterialOrientation::Rectangular { axis1: [0.0, 1.0, 0.0], axis2: [-1.0, 0.0, 0.0] }.axes_at([0.0; 3]).unwrap();
    let d = elasticity_in_axes(&material, &axes).unwrap();
    for (value, expected) in [(d[0][0], 10.0), (d[1][1], 100.0), (d[4][4], 3.0), (d[5][5], 4.0)] {
        assert!((value - expected).abs() < 3e-14 * expected);
    }
    let axes = MaterialOrientation::Rectangular { axis1: [1.0, 1.0, 0.0], axis2: [-1.0, 1.0, 0.0] }.axes_at([0.0; 3]).unwrap();
    let d = elasticity_in_axes(&material, &axes).unwrap();
    let sigma = stress(&d, &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    for (s, reference) in sigma.iter().zip([32.5, 22.5, 0.0, 22.5, 0.0, 0.0]) {
        assert!((s - reference).abs() < 3e-14 * 100.0);
    }
    // Work in material axes: ε1=ε2=.5, γ12=-1, σ1=50, σ2=5, τ12=-5.
    assert!((sigma[0] - (0.5 * 50.0 + 0.5 * 5.0 + 5.0)).abs() < 3e-14 * 100.0);
    for (i, row) in d.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            assert_eq!(*value, d[j][i]);
        }
    }
    let iso = Elastic::Isotropic { young: 100.0, poisson: 0.3 };
    let a = elasticity(&iso).unwrap();
    let b = elasticity_in_axes(&iso, &axes).unwrap();
    for (x, y) in a.iter().flatten().zip(b.iter().flatten()) {
        assert!((x - y).abs() < 1e-13);
    }
    assert!(elasticity_in_axes(&material, &[[1.0; 3]; 3]).is_err());
    assert!(elasticity_in_axes(&material, &[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]).is_err());
}
