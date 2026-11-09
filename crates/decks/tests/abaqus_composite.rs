//! Original small meshes; no vendor example input or commercial solver output.
use fieldcraft_decks::abaqus;
use fieldcraft_model::units::Preset;
use fieldcraft_model::*;

fn plate() -> Model {
    let mut m = Model::default();
    for (i, p) in [[0.0, 0.0, 0.0], [2.0, 1.0, 0.0], [3.0, 3.0, 0.0], [0.0, 2.0, 0.0]].into_iter().enumerate() {
        m.set_node(i as u32 + 1, p).unwrap();
    }
    m.set_element(1, ElementKind::Quad4, vec![1, 2, 3, 4]).unwrap();
    m.add_set("all", Set { elements: vec![1], ..Set::default() }).unwrap();
    m.add_material(
        "Carbon",
        Material {
            elastic: Some(Elastic::EngineeringConstants { young: [120e9, 9e9, 7e9], poisson: [0.25, 0.2, 0.3], shear: [5e9, 4e9, 3e9] }),
            density: Some(1600.0),
            ..Material::default()
        },
    )
    .unwrap();
    m.add_section(
        "laminate",
        Section {
            set: "all".into(),
            material: String::new(),
            orientation: None,
            formulation: Formulation::default(),
            kind: SectionKind::CompositeShell {
                layup: Layup {
                    offset: 0.25,
                    plies: vec![
                        Ply { material: "Carbon".into(), thickness: 0.001, points: 3, angle: std::f64::consts::FRAC_PI_6 },
                        Ply { material: "Carbon".into(), thickness: 0.002, points: 5, angle: -std::f64::consts::FRAC_PI_4 },
                    ],
                },
            },
        },
    )
    .unwrap();
    m.add_boundary_condition("hold", BoundaryCondition::Fixed { target: Target::Set("all".into()) }).unwrap();
    let mut step = Step::new("Load", Procedure::Static {});
    step.boundary_conditions.insert("hold".into());
    m.add_step(step, None).unwrap();
    m
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-12 * b.abs().max(1.0), "{a} vs {b}");
}

#[test]
fn composite_roundtrip_preserves_stack_offset_angles_and_shear_order_in_each_unit_system() {
    let m = plate();
    for units in [Preset::MKgS, Preset::MmTonneS, Preset::MmKgMs, Preset::InLbfS] {
        let text = abaqus::write(&m, units).unwrap();
        assert!(text.contains("COMPOSITE, OFFSET=0.25") && text.contains("Carbon, 30.") && text.contains("Carbon, -45."), "{text}");
        let read = abaqus::read_in(&text, units).unwrap();
        assert!(
            read.unknown
                .iter()
                .all(|u| ["STEP", "STATIC", "BOUNDARY", "OUTPUT", "NODEOUTPUT", "ELEMENTOUTPUT", "ENDSTEP"].contains(&u.keyword.as_str())),
            "{:?}",
            read.unknown
        );
        assert_eq!(read.model.set("all").unwrap().elements, vec![1]);
        let section = read.model.sections().values().next().unwrap();
        assert!(section.material.is_empty());
        let SectionKind::CompositeShell { layup } = &section.kind else { panic!("composite") };
        close(layup.offset, 0.25);
        for (ply, (t, angle, points)) in layup.plies.iter().zip([(0.001, std::f64::consts::FRAC_PI_6, 3), (0.002, -std::f64::consts::FRAC_PI_4, 5)]) {
            close(ply.thickness, t);
            close(ply.angle, angle);
            assert_eq!(ply.points, points);
            assert_eq!(ply.material, "Carbon");
        }
        let Elastic::EngineeringConstants { shear, .. } = read.model.material("Carbon").unwrap().elastic.unwrap() else { panic!("orthotropic") };
        for (a, b) in shear.into_iter().zip([5e9, 4e9, 3e9]) {
            close(a, b);
        }
    }
}

#[test]
fn implicit_quad_frames_follow_opposite_edge_midpoints_after_node_permutations() {
    for (nodes, along, normal) in [
        (vec![1, 2, 3, 4], [5.0, 2.0, 0.0], [0.0, 0.0, 1.0]),
        (vec![2, 3, 4, 1], [1.0, 4.0, 0.0], [0.0, 0.0, 1.0]),
        (vec![1, 4, 3, 2], [1.0, 4.0, 0.0], [0.0, 0.0, -1.0]),
    ] {
        let mut m = plate();
        m.set_element(1, ElementKind::Quad4, nodes).unwrap();
        let text = abaqus::write(&m, Preset::MKgS).unwrap();
        let read = abaqus::read(&text).unwrap();
        let o = read.model.sections().values().next().unwrap().orientation.unwrap();
        let axes = o.frame.shell_axes_at([1.0, 1.0, 0.0], normal, o.shell_reference_axis).unwrap();
        let length = along.iter().map(|v| v * v).sum::<f64>().sqrt();
        for i in 0..3 {
            close(axes[0][i], along[i] / length);
            close(axes[2][i], normal[i]);
        }
        assert!(text.contains("3, 0.0"));
    }
    let mut m = plate();
    m.set_element(1, ElementKind::Tri3, vec![2, 3, 1]).unwrap();
    let text = abaqus::write(&m, Preset::MKgS).unwrap();
    let read = abaqus::read(&text).unwrap();
    let axes = read.model.sections().values().next().unwrap().orientation.unwrap().frame.axes_at([0.0; 3]).unwrap();
    close(axes[0][0], 1.0 / 5.0f64.sqrt());
    close(axes[0][1], 2.0 / 5.0f64.sqrt());
}

#[test]
fn rectangular_and_cylindrical_reference_axes_use_the_abaqus_cyclic_projection_rule() {
    for (frame, reference, rotation) in [
        (MaterialOrientation::Rectangular { axis1: [0.0, 0.0, 1.0], axis2: [1.0, 1.0, 0.0] }, 2, 1),
        (MaterialOrientation::Rectangular { axis1: [0.0, 0.0, 1.0], axis2: [1.0, 1.0, 0.0] }, 3, 2),
        (MaterialOrientation::Cylindrical { origin: [-1.0, -2.0, 0.0], axis: [0.0, 0.0, 1.0] }, 2, 1),
    ] {
        let mut m = plate();
        let mut s = m.section("laminate").unwrap().clone();
        s.orientation = Some(SectionOrientation { frame, shell_reference_axis: reference });
        m.replace_section("laminate", s).unwrap();
        let text = abaqus::write(&m, Preset::MmTonneS).unwrap();
        assert!(text.contains(&format!("\n{rotation}, 0.0\n")), "{text}");
        let read = abaqus::read_in(&text, Preset::MmTonneS).unwrap();
        let o = read.model.sections().values().next().unwrap().orientation.unwrap();
        assert_eq!(o.shell_reference_axis, reference);
        for p in [[0.5, 0.5, 0.0], [2.0, 1.5, 0.0]] {
            let want = frame.shell_axes_at([1.25, 1.5, 0.0], [0.0, 0.0, 1.0], reference).unwrap();
            let got = o.frame.shell_axes_at(p, [0.0, 0.0, 1.0], o.shell_reference_axis).unwrap();
            for (a, b) in got.into_iter().flatten().zip(want.into_iter().flatten()) {
                close(a, b);
            }
        }
    }
}

#[test]
fn unsupported_layer_and_orientation_semantics_are_refused_at_their_lines() {
    let text = abaqus::write(&plate(), Preset::MKgS).unwrap();
    for (from, to, why) in [
        ("COMPOSITE, OFFSET", "COMPOSITE=YES, OFFSET", "flag"),
        ("COMPOSITE, OFFSET", "COMPOSITE, SYMMETRIC, OFFSET", "SYMMETRIC"),
        ("COMPOSITE, OFFSET", "COMPOSITE, SECTION INTEGRATION=GAUSS, OFFSET", "Gauss"),
        ("DEFINITION=COORDINATES", "DEFINITION=NODES", "COORDINATES"),
        ("ORIENTATION=FC_ORI_1", "ORIENTATION=Missing", "Missing"),
        ("0.001, 3, Carbon", "0.001, 2, Carbon", "odd Simpson"),
    ] {
        let bad = text.replace(from, to);
        assert_ne!(bad, text, "replacement {from}");
        let error = abaqus::read(&bad).unwrap_err().to_string();
        assert!(error.contains("line ") && error.contains(why), "{error}");
    }
    let rotated = text.replace("3, 0.0\n*SHELL", "3, 10.0\n*SHELL");
    assert_ne!(rotated, text);
    let read = abaqus::read::read_analysis_limited(&rotated, Preset::MKgS, abaqus::read::Limits::default()).unwrap();
    assert!(read.analysis_complete);
    let normalized = abaqus::write(&read.model, Preset::MKgS).unwrap();
    let back = abaqus::read(&normalized).unwrap();
    let a = read.model.sections().values().next().unwrap().orientation.unwrap();
    let b = back.model.sections().values().next().unwrap().orientation.unwrap();
    for (x, y) in a.frame.axes_at([0.; 3]).unwrap().into_iter().flatten().zip(b.frame.axes_at([0.; 3]).unwrap().into_iter().flatten()) {
        close(x, y);
    }
    assert!(abaqus::read("*FAIL STRESS\n1,2,3,4,5,0\n").unwrap_err().to_string().contains("S23"));
    let mut one_point = plate();
    let mut s = one_point.section("laminate").unwrap().clone();
    if let SectionKind::CompositeShell { layup } = &mut s.kind {
        layup.plies[0].points = 1;
    }
    one_point.replace_section("laminate", s).unwrap();
    assert!(abaqus::write(&one_point, Preset::MKgS).unwrap_err().to_string().contains("ply 1"));
}

#[test]
fn a_homogeneous_orthotropic_shell_retains_its_material_frame() {
    let mut m = plate();
    let mut s = m.section("laminate").unwrap().clone();
    s.material = "Carbon".into();
    s.kind = SectionKind::Shell { thickness: 0.003, offset: -0.5, points: 5 };
    s.orientation = Some(SectionOrientation {
        frame: MaterialOrientation::Rectangular { axis1: [1.0, 2.0, 0.0], axis2: [-2.0, 1.0, 0.0] },
        shell_reference_axis: 1,
    });
    m.replace_section("laminate", s.clone()).unwrap();
    let text = abaqus::write(&m, Preset::MmTonneS).unwrap();
    let read = abaqus::read_in(&text, Preset::MmTonneS).unwrap();
    let loaded = read.model.sections().values().next().unwrap();
    assert_eq!(loaded.kind, s.kind);
    assert_eq!(loaded.orientation, s.orientation);
    assert_eq!(loaded.material, s.material);
}
