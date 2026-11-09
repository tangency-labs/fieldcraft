//! Simply supported elastic column: omega^2 = EI/rho/A * (pi/L)^4
//! plus P/rho/A * (pi/L)^2, tension positive (Bathe 2014 §6.8.2;
//! Timoshenko & Gere 1961 §2.1). A half-model keeps the first symmetric pair.
mod common;
use common::*;
use fieldcraft_model::*;
use fieldcraft_solver::{FrequencyOptions, Limits, solve_frequency};

fn column(alpha: f64) -> Model {
    let mut g = grid([2, 2, 10], ElementKind::Hex20, |x| [x[0], x[1], 10.0 * x[2]]);
    g.define(iso(1000.0, 0.0), Some(1.0), Formulation::default(), Procedure::Static {});
    g.face_set("end", |x| near(x[2], 0.0));
    // L=20, I=1/12, A=1. Positive pressure at the bottom is compression.
    let critical = std::f64::consts::PI.powi(2) * 1000.0 / (12.0 * 20.0f64.powi(2));
    g.load("axial", Load::Pressure { set: "end".into(), magnitude: -alpha * critical });
    g.node_set("pin", |x| near(x[2], 0.0));
    g.bc("pin", BoundaryCondition::Displacement { target: Target::Set("pin".into()), values: [Some(0.0), Some(0.0), None, None, None, None] });
    g.symmetry("mid", Symmetry::Z, |x| near(x[2], 1.0));
    let mut vibration = Step::new("Modes", Procedure::Frequency { modes: 2, min_hz: None, max_hz: None });
    vibration.boundary_conditions = g.model.step("Step").unwrap().boundary_conditions.clone();
    vibration.base_state = Some("Step".into());
    g.model.add_step(vibration, None).unwrap();
    g.model
}

#[test]
fn axial_force_stiffens_and_softens_the_same_column() {
    let options = FrequencyOptions::default();
    let limits = Limits::default();
    let zero = solve_frequency(&column(0.0), Some("Modes"), &options, &limits).unwrap();
    for alpha in [-0.9, -0.5, 0.5] {
        let result = solve_frequency(&column(alpha), Some("Modes"), &options, &limits).unwrap();
        let ratio = result.modes[0].eigenvalue / zero.modes[0].eigenvalue;
        // 2% absolute in squared-frequency ratio allows 3D shear/rotary effects
        // while resolving the sign and the approach to the Euler load.
        assert!((ratio - (1.0 + alpha)).abs() < 0.02, "alpha={alpha}, ratio={ratio}");
        assert!((result.modes[1].eigenvalue / result.modes[0].eigenvalue - 1.0).abs() < 1e-5);
    }
}

/// Two bolt halves, both outer ends held. q = F(2L)/(EA), A=L=1.
fn cut_bolt() -> Model {
    let mut model = Model::default();
    let corner = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]];
    for half in 0..2u32 {
        for (i, x) in corner.iter().enumerate() {
            model.set_node(half * 8 + i as u32 + 1, [x[0], x[1], x[2] + f64::from(half)]).unwrap();
        }
        model.set_element(half + 1, ElementKind::Hex8, (half * 8 + 1..=half * 8 + 8).collect()).unwrap();
    }
    for (name, nodes) in [("first", vec![5, 6, 7, 8]), ("second", vec![9, 10, 11, 12]), ("ends", vec![1, 2, 3, 4, 13, 14, 15, 16])] {
        model.add_set(name, Set { nodes, ..Set::default() }).unwrap();
    }
    model.add_set("all", Set { elements: vec![1, 2], ..Set::default() }).unwrap();
    model.add_material("Steel", Material { elastic: Some(iso(1000., 0.)), density: Some(1.), ..Material::default() }).unwrap();
    model
        .add_section(
            "Bolt",
            Section {
                set: "all".into(),
                material: "Steel".into(),
                kind: SectionKind::Solid {},
                orientation: None,
                formulation: Formulation::default(),
            },
        )
        .unwrap();
    model.add_reference_point("Control", ReferencePoint { xyz: [0.5, 0.5, 1.], label: 100 }).unwrap();
    model
        .add_pretension(
            "Cut",
            PretensionSection { first: "first".into(), second: "second".into(), face: None, normal: [0., 0., 1.], control: "Control".into() },
        )
        .unwrap();
    model.add_boundary_condition("Ends", BoundaryCondition::Fixed { target: Target::Set("ends".into()) }).unwrap();
    let mut preload = Step::new("Preload", Procedure::Static {});
    preload.boundary_conditions.insert("Ends".into());
    preload.pretension.insert("Cut".into(), PretensionControl::Force { value: 1. });
    model.add_step(preload, None).unwrap();
    model
}

#[test]
fn cut_force_balance_and_locked_adjustment_are_work_conjugate() {
    let mut model = cut_bolt();
    let load = fieldcraft_solver::solve_static(&model, Some("Preload"), &Limits::default()).unwrap();
    let q = load.reference_points.iter().find(|p| p.name == "Control").unwrap().displacement[0];
    assert!((q / 0.002 - 1.).abs() < 1e-10, "q={q}");
    assert!((load.pretension[0].section_force - 1.).abs() < 1e-9);
    assert!((load.pretension[0].applied_force - 1.).abs() < 1e-12);
    for e in &load.elements {
        for p in &e.points {
            assert!((p.stress[2] - 1.).abs() < 1e-9);
        }
    }
    assert!(load.applied_total.iter().all(|v| v.abs() < 1e-12), "scalar cut force is an internal load pair");
    let mut locked = model.step("Preload").unwrap().clone();
    locked.name = "Locked".into();
    locked.pretension.insert("Cut".into(), PretensionControl::Adjustment { value: q });
    model.add_step(locked, None).unwrap();
    let held = fieldcraft_solver::solve_static(&model, Some("Locked"), &Limits::default()).unwrap();
    let force = held.reference_points.iter().find(|p| p.name == "Control").unwrap().reaction[0];
    assert!((force - 1.).abs() < 1e-9, "section force={force}");
    assert!((held.pretension[0].section_force - 1.).abs() < 1e-9);
    assert_eq!(held.pretension[0].applied_force, 0.);
    for (a, b) in load.displacement.iter().flatten().zip(held.displacement.iter().flatten()) {
        assert!((a - b).abs() < 1e-12);
    }
    let mut modes = Step::new("Modes", Procedure::Frequency { modes: 6, min_hz: None, max_hz: None });
    modes.boundary_conditions.insert("Ends".into());
    modes.base_state = Some("Locked".into());
    model.add_step(modes, None).unwrap();
    let vibration = solve_frequency(&model, Some("Modes"), &FrequencyOptions::default(), &Limits::default()).unwrap();
    assert_eq!(vibration.modes.len(), 6);
    assert!(vibration.modes.iter().all(|m| m.frequency > 0. && m.frequency.is_finite()));
}
