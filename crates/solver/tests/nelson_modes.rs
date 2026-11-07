mod common;
use common::*;
use fieldcraft_model::{
    Elastic, ElementKind, Procedure,
    design::{MaterialConstant, Parameter},
};
use fieldcraft_solver::{Limits, sensitivity::frequency_sensitivity_with_vectors, solve_frequency};

#[test]
fn nelson_direction_change_matches_sign_aligned_mass_normalized_modes() {
    let mut g = grid([3, 1, 1], ElementKind::Hex8, |p| [p[0] * 2.0, p[1] * 0.3, p[2] * 0.2]);
    g.define(
        Elastic::EngineeringConstants { young: [1000.0, 800.0, 600.0], poisson: [0.2, 0.15, 0.1], shear: [300.0, 250.0, 200.0] },
        Some(2.0),
        full(),
        Procedure::Frequency { modes: 3, min_hz: None, max_hz: None },
    );
    g.fixed("fixed", |p| near(p[0], 0.0));
    let parameter = Parameter::Material { name: "Mat".into(), constant: MaterialConstant::G12 };
    for h in [0.3, 0.03] {
        let base = frequency_sensitivity_with_vectors(&g.model, None, &parameter, h, &Default::default(), true, &Limits::default()).unwrap();
        let (mut plus, mut minus) = (g.model.clone(), g.model.clone());
        parameter.set(&mut plus, 300.0 + h).unwrap();
        parameter.set(&mut minus, 300.0 - h).unwrap();
        let a = solve_frequency(&plus, None, &Default::default(), &Limits::default()).unwrap();
        let b = solve_frequency(&minus, None, &Default::default(), &Limits::default()).unwrap();
        for ((mode, ap), am) in base.modes.iter().zip(&a.modes).zip(&b.modes) {
            let v = mode.vector.as_ref().unwrap();
            let sign = |shape: &[[f64; 3]]| {
                let dot: f64 = shape.iter().flatten().zip(v.shape_si.iter().flatten()).map(|(a, b)| a * b).sum();
                if dot < 0.0 { -1.0 } else { 1.0 }
            };
            let sa = sign(&ap.shape);
            let sb = sign(&am.shape);
            let fd: Vec<_> = ap.shape.iter().flatten().zip(am.shape.iter().flatten()).map(|(a, b)| (sa * a - sb * b) / (2.0 * h)).collect();
            let scale = fd.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
            assert!(scale > 1e-8);
            let error = v.derivative_si.iter().flatten().zip(fd).map(|(a, b)| (a - b).abs()).fold(0.0_f64, f64::max);
            // This weak direction derivative is about 2.5e-7 per Pa: retain a
            // 2e-11 absolute FD round-off floor as well as a 1e-4 relative bound.
            assert!(error < 1e-4 * scale + 2e-11, "mode {}, error {error}, scale {scale}", mode.frequency);
        }
    }
}
