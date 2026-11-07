//! Independent affine-field energy and tensor checks for spatially varying solid axes.

use fieldcraft_fem::dense::Matrix;
use fieldcraft_fem::element::{Formulation, Solid};
use fieldcraft_fem::material::elasticity_in_axes;
use fieldcraft_model::{Elastic, ElementKind, MaterialOrientation};

const MATERIAL: Elastic = Elastic::EngineeringConstants { young: [120.0, 10.0, 20.0], poisson: [0.0; 3], shear: [5.0, 4.0, 3.0] };
const CYLINDER: MaterialOrientation = MaterialOrientation::Cylindrical { origin: [0.0; 3], axis: [0.0, 0.0, 1.0] };
const RECTANGULAR: MaterialOrientation = MaterialOrientation::Rectangular { axis1: [2.0, 1.0, 0.0], axis2: [-1.0, 2.0, 0.0] };
const STRAIN: [f64; 6] = [0.02, -0.01, 0.03, 0.04, -0.02, 0.01];

fn brick(kind: ElementKind) -> Vec<[f64; 3]> {
    // An affine 2 × 2 × 1 brick: det J = 1/2, volume = 4. Its material directions
    // vary substantially between Gauss points, while the whole brick stays off the axis.
    let corners: Vec<_> = fieldcraft_fem::shape::HEX_CORNERS.iter().map(|p| [2.0 + p[0], 0.2 + p[1], 0.5 * p[2]]).collect();
    let mut xyz = corners.clone();
    for edge in kind.edges() {
        if let [a, b, _] = **edge {
            xyz.push(std::array::from_fn(|k| (corners[a][k] + corners[b][k]) / 2.0));
        }
    }
    xyz
}

fn affine(x: [f64; 3]) -> [f64; 3] {
    let [xx, yy, zz, xy, xz, yz] = STRAIN;
    [xx * x[0] + xy * x[1] / 2.0 + xz * x[2] / 2.0, xy * x[0] / 2.0 + yy * x[1] + yz * x[2] / 2.0, xz * x[0] / 2.0 + yz * x[1] / 2.0 + zz * x[2]]
}

fn energy(k: &Matrix, u: &[[f64; 3]]) -> f64 {
    let flat: Vec<_> = u.iter().flatten().copied().collect();
    flat.iter().zip(k.mul_vec(&flat)).map(|(a, b)| a * b).sum()
}

/// Explicit planar tensor formula, independent of axes_at and elasticity_in_axes.
fn cylindrical_stress(x: [f64; 3]) -> [f64; 6] {
    let radius = x[0].hypot(x[1]);
    let (c, s) = (x[0] / radius, x[1] / radius);
    let [xx, yy, zz, xy, xz, yz] = STRAIN;
    let sr = 120.0 * (c * c * xx + s * s * yy + c * s * xy);
    let st = 10.0 * (s * s * xx + c * c * yy - c * s * xy);
    let sz = 20.0 * zz;
    let trt = 5.0 * (2.0 * c * s * (yy - xx) + (c * c - s * s) * xy);
    let trz = 4.0 * (c * xz + s * yz);
    let ttz = 3.0 * (-s * xz + c * yz);
    [
        c * c * sr + s * s * st - 2.0 * c * s * trt,
        s * s * sr + c * c * st + 2.0 * c * s * trt,
        sz,
        c * s * (sr - st) + (c * c - s * s) * trt,
        c * trz - s * ttz,
        s * trz + c * ttz,
    ]
}

#[test]
fn rectangular_axes_reproduce_constant_elasticity_and_initial_stress() {
    let d = elasticity_in_axes(&MATERIAL, &RECTANGULAR.axes_at([0.0; 3]).unwrap()).unwrap();
    for (kind, formulation) in [
        (ElementKind::Hex8, Formulation::Full),
        (ElementKind::Hex8, Formulation::IncompatibleModes),
        (ElementKind::Hex20, Formulation::Full),
        (ElementKind::Hex20, Formulation::Reduced),
    ] {
        let xyz = brick(kind);
        let solid = Solid::new(kind, &xyz, formulation).unwrap();
        let u: Vec<_> = xyz.iter().copied().map(affine).collect();
        let initial = vec![[0.4, -0.2, 0.1, 0.3, -0.5, 0.6]; xyz.len()];
        assert_eq!(solid.stiffness_in_material_axes(&MATERIAL, &RECTANGULAR).unwrap().data, solid.stiffness(&d).unwrap().data);
        assert_eq!(
            solid.stresses_in_material_axes(&MATERIAL, &RECTANGULAR, &u, Some(&initial)).unwrap(),
            solid.stresses_from(&d, &u, None, Some(&initial)).unwrap()
        );
    }
}

#[test]
fn cylindrical_axes_integrate_the_pointwise_energy_and_stress() {
    for (kind, formulation) in
        [(ElementKind::Hex8, Formulation::Full), (ElementKind::Hex20, Formulation::Full), (ElementKind::Hex20, Formulation::Reduced)]
    {
        let xyz = brick(kind);
        let solid = Solid::new(kind, &xyz, formulation).unwrap();
        let k = solid.stiffness_in_material_axes(&MATERIAL, &CYLINDER).unwrap();
        let u: Vec<_> = xyz.iter().copied().map(affine).collect();
        let expected: f64 = solid
            .rule()
            .points()
            .iter()
            .map(|p| {
                let sigma = cylindrical_stress(solid.position(p.xi));
                0.5 * p.weight * STRAIN.iter().zip(sigma).map(|(a, b)| a * b).sum::<f64>()
            })
            .sum();
        assert!((energy(&k, &u) - expected).abs() < 2e-13 * expected, "{kind:?} {formulation:?}");
        // A centroid-only frame would pass an affine strain test, but gives the wrong
        // energy for this direction field. This case detects that shortcut.
        let centroid: f64 = 4.0 * STRAIN.iter().zip(cylindrical_stress([2.0, 0.2, 0.0])).map(|(a, b)| a * b).sum::<f64>();
        assert!((centroid - expected).abs() > 0.01 * expected);
        let initial: Vec<_> = xyz.iter().map(|x| [x[0] + 2.0 * x[1], 0.0, 0.0, 0.0, x[2], 0.0]).collect();
        for p in solid.stresses_in_material_axes(&MATERIAL, &CYLINDER, &u, Some(&initial)).unwrap() {
            for (value, exact) in p.strain.iter().zip(STRAIN) {
                assert!((value - exact).abs() < 3e-15);
            }
            let mut exact = cylindrical_stress(p.xyz);
            exact[0] += p.xyz[0] + 2.0 * p.xyz[1];
            exact[4] += p.xyz[2];
            for (value, exact) in p.stress.iter().zip(exact) {
                assert!((value - exact).abs() < 2e-13);
            }
        }
    }
}

#[test]
fn cylindrical_stiffness_has_symmetric_blocks_and_six_rigid_modes() {
    let xyz = brick(ElementKind::Hex20);
    let solid = Solid::new(ElementKind::Hex20, &xyz, Formulation::Full).unwrap();
    let k = solid.stiffness_in_material_axes(&MATERIAL, &CYLINDER).unwrap();
    let scale = k.data.iter().map(|v| v.abs()).fold(0.0, f64::max);
    for i in 0..k.n {
        for j in 0..k.n {
            assert_eq!(k.get(i, j), k.get(j, i));
        }
    }
    for axis in 0..3 {
        let unit: [f64; 3] = std::array::from_fn(|k| if k == axis { 1.0 } else { 0.0 });
        let translation: Vec<_> = xyz.iter().flat_map(|_| unit).collect();
        let rotation: Vec<_> =
            xyz.iter().flat_map(|x| [unit[1] * x[2] - unit[2] * x[1], unit[2] * x[0] - unit[0] * x[2], unit[0] * x[1] - unit[1] * x[0]]).collect();
        for mode in [translation, rotation] {
            assert!(k.mul_vec(&mode).iter().all(|v| v.abs() < 2e-13 * scale));
        }
    }
}

const ROTATION: [[f64; 3]; 3] = [[0.6, -0.48, 0.64], [0.8, 0.36, -0.48], [0.0, 0.8, 0.6]];

fn rotate(x: [f64; 3]) -> [f64; 3] {
    ROTATION.map(|r| r.iter().zip(x).map(|(a, b)| a * b).sum())
}

#[test]
fn cylindrical_material_objectivity_rotates_coordinates_axis_and_tensor() {
    let xyz = brick(ElementKind::Hex8);
    let shift = [0.3, -0.8, 0.4];
    let transformed: Vec<_> = xyz.iter().map(|x| std::array::from_fn(|k| rotate(*x)[k] + shift[k])).collect();
    let rotated_axis = MaterialOrientation::Cylindrical { origin: shift, axis: rotate([0.0, 0.0, 1.0]) };
    let solid = Solid::new(ElementKind::Hex8, &xyz, Formulation::Full).unwrap();
    let turned = Solid::new(ElementKind::Hex8, &transformed, Formulation::Full).unwrap();
    let k = solid.stiffness_in_material_axes(&MATERIAL, &CYLINDER).unwrap();
    let rotated_k = turned.stiffness_in_material_axes(&MATERIAL, &rotated_axis).unwrap();
    let scale = k.data.iter().map(|v| v.abs()).fold(0.0, f64::max);
    for i in 0..xyz.len() {
        for j in 0..xyz.len() {
            for (a, row) in ROTATION.iter().enumerate() {
                for (b, col) in ROTATION.iter().enumerate() {
                    let expected: f64 = row
                        .iter()
                        .enumerate()
                        .map(|(c, r)| col.iter().enumerate().map(|(d, s)| r * s * k.get(3 * i + c, 3 * j + d)).sum::<f64>())
                        .sum();
                    assert!((rotated_k.get(3 * i + a, 3 * j + b) - expected).abs() < 2e-13 * scale);
                }
            }
        }
    }
    let u: Vec<_> = xyz.iter().copied().map(affine).collect();
    let turned_u: Vec<_> = u.iter().copied().map(rotate).collect();
    let points = solid.stresses_in_material_axes(&MATERIAL, &CYLINDER, &u, None).unwrap();
    let turned_points = turned.stresses_in_material_axes(&MATERIAL, &rotated_axis, &turned_u, None).unwrap();
    for (p, q) in points.iter().zip(turned_points) {
        let [xx, yy, zz, xy, xz, yz] = p.stress;
        let sigma = [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]];
        for (value, (a, b)) in q.stress.iter().zip([(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]) {
            let expected: f64 = (0..3).map(|c| (0..3).map(|d| ROTATION[a][c] * sigma[c][d] * ROTATION[b][d]).sum::<f64>()).sum();
            assert!((value - expected).abs() < 2e-13);
        }
    }
}

#[test]
fn bad_axes_fields_and_unsupported_formulations_are_errors() {
    let xyz = brick(ElementKind::Hex8);
    let solid = Solid::new(ElementKind::Hex8, &xyz, Formulation::Full).unwrap();
    let u = vec![[0.0; 3]; xyz.len()];
    let axis_through_point = MaterialOrientation::Cylindrical { origin: solid.position(solid.rule().points()[0].xi), axis: [0.0, 0.0, 1.0] };
    for bad in [
        axis_through_point,
        MaterialOrientation::Cylindrical { origin: [0.0; 3], axis: [0.0; 3] },
        MaterialOrientation::Rectangular { axis1: [1.0, 0.0, 0.0], axis2: [f64::NAN; 3] },
    ] {
        assert!(solid.stiffness_in_material_axes(&MATERIAL, &bad).is_err());
        assert!(solid.stresses_in_material_axes(&MATERIAL, &bad, &u, None).is_err());
    }
    for axes in [RECTANGULAR, CYLINDER] {
        assert!(solid.stresses_in_material_axes(&MATERIAL, &axes, &u[..7], None).is_err());
        assert!(solid.stresses_in_material_axes(&MATERIAL, &axes, &[[f64::NAN; 3]; 8], None).is_err());
        assert!(solid.stresses_in_material_axes(&MATERIAL, &axes, &u, Some(&[[0.0; 6]; 7])).is_err());
        assert!(solid.stresses_in_material_axes(&MATERIAL, &axes, &u, Some(&[[f64::INFINITY; 6]; 8])).is_err());
        let overflow: Vec<_> = xyz.iter().map(|x| x.map(|v| v * 1e307)).collect();
        assert!(solid.stresses_in_material_axes(&MATERIAL, &axes, &overflow, None).is_err());
    }
    let large_xyz: Vec<_> = xyz.iter().map(|x| x.map(|v| v * 1e6)).collect();
    let large = Solid::new(ElementKind::Hex8, &large_xyz, Formulation::Full).unwrap();
    let stiff = Elastic::Isotropic { young: 1e307, poisson: 0.0 };
    for axes in [RECTANGULAR, CYLINDER] {
        assert!(large.stiffness_in_material_axes(&stiff, &axes).is_err(), "overflow cannot return a nonfinite matrix");
    }
    for formulation in [Formulation::IncompatibleModes, Formulation::Hybrid { pressures: 1 }] {
        let solid = Solid::new(ElementKind::Hex8, &xyz, formulation).unwrap();
        let err = solid.stiffness_in_material_axes(&MATERIAL, &CYLINDER).unwrap_err().to_string();
        assert!(err.contains("spatially varying material axes"), "{err}");
        assert!(solid.stresses_in_material_axes(&MATERIAL, &CYLINDER, &u, None).is_err());
    }
}
