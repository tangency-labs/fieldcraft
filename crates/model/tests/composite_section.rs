//! Laminate definitions preserve every ply reference and dimensional input at the model boundary.
use fieldcraft_model::rename::Kind;
use fieldcraft_model::units::{DisplaySystem, Preset, Quantity};
use fieldcraft_model::*;

fn strength() -> PlyStrength {
    PlyStrength { tension: [1.5e9, 50e6], compression: [1.2e9, 150e6], shear12: 80e6, shear23: 60e6, tsai_wu_interaction: -0.25 }
}

fn plate() -> Model {
    let mut m = Model::default();
    for (i, xyz) in [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]].into_iter().enumerate() {
        m.set_node(i as u32 + 1, xyz).unwrap();
    }
    m.set_element(1, ElementKind::Quad4, vec![1, 2, 3, 4]).unwrap();
    m.add_set("all", Set { elements: vec![1], ..Set::default() }).unwrap();
    for name in ["Carbon", "Glass"] {
        m.add_material(
            name,
            Material {
                elastic: Some(Elastic::Isotropic { young: 70e9, poisson: 0.25 }),
                density: Some(1600.0),
                ply_strength: Some(strength()),
                ..Material::default()
            },
        )
        .unwrap();
    }
    m.add_section("laminate", section()).unwrap();
    m.add_boundary_condition("hold", BoundaryCondition::Fixed { target: Target::Set("all".into()) }).unwrap();
    let mut step = Step::new("Modes", Procedure::Frequency { modes: 1, min_hz: None, max_hz: None });
    step.boundary_conditions.insert("hold".into());
    m.add_step(step, None).unwrap();
    m
}

fn section() -> Section {
    Section {
        set: "all".into(),
        material: String::new(),
        formulation: Formulation::default(),
        orientation: Some(SectionOrientation {
            frame: MaterialOrientation::Cylindrical { origin: [0.1, 0.2, 0.3], axis: [0.0, 0.0, 2.0] },
            shell_reference_axis: 3,
        }),
        kind: SectionKind::CompositeShell {
            layup: Layup {
                plies: vec![
                    Ply { material: "Carbon".into(), thickness: 0.001, angle: std::f64::consts::FRAC_PI_4, points: 3 },
                    Ply { material: "Glass".into(), thickness: 0.002, angle: 0.0, points: 1 },
                    Ply { material: "Carbon".into(), thickness: 0.001, angle: -std::f64::consts::FRAC_PI_4, points: 3 },
                ],
                offset: 0.25,
            },
        },
    }
}

#[test]
fn composite_json_omits_the_section_material_and_keeps_the_stack() {
    let section = section();
    let json = serde_json::to_value(&section).unwrap();
    assert!(json.get("material").is_none());
    assert_eq!(json["kind"]["type"], "composite_shell");
    assert_eq!(serde_json::from_value::<Section>(json).unwrap(), section);
    let m = plate();
    m.resolve().unwrap();
    let loaded: Model = serde_json::from_str(&serde_json::to_string(&m).unwrap()).unwrap();
    assert_eq!(loaded, m);
    loaded.resolve().unwrap();
    let with_material = Section { material: "Carbon".into(), ..section.clone() };
    assert!(with_material.validate().unwrap_err().to_string().contains("omit the section material"));
    assert!(Section { kind: SectionKind::Solid {}, orientation: None, ..section }.validate().is_err());
    let default: SectionOrientation = serde_json::from_str(r#"{"frame":{"type":"rectangular","axis1":[1,0,0],"axis2":[0,1,0]}}"#).unwrap();
    assert_eq!(default.shell_reference_axis, 1);
    assert!(serde_json::from_str::<SectionOrientation>(r#"{"frame":{"type":"cylindrical","origin":[0,0,0],"axis":[0,0,1]},"extra":1}"#).is_err());
}

#[test]
fn material_references_rename_and_deletion_include_every_ply() {
    let mut m = plate();
    assert_eq!(m.material_users("Carbon"), ["section `laminate`"]);
    assert_eq!(m.material_users("Glass"), ["section `laminate`"]);
    assert!(m.remove_material("Carbon").unwrap_err().to_string().contains("laminate"));
    m.rename(Kind::Material, "Carbon", "Fibre").unwrap();
    assert_eq!(m.section("laminate").unwrap().material_names().collect::<Vec<_>>(), ["Fibre", "Glass", "Fibre"]);
    assert!(m.material_users("Carbon").is_empty());
    let mut missing = m.section("laminate").unwrap().clone();
    if let SectionKind::CompositeShell { layup } = &mut missing.kind {
        layup.plies[1].material = "Missing".into();
    }
    let error = m.replace_section("laminate", missing).unwrap_err().to_string();
    assert!(error.contains("laminate") && error.contains("Missing"), "{error}");
    m.remove_section("laminate").unwrap();
    m.remove_material("Fibre").unwrap();
    m.remove_material("Glass").unwrap();
}

#[test]
fn every_ply_material_is_checked_for_elasticity_and_mass() {
    let mut m = plate();
    let mut glass = m.material("Glass").unwrap().clone();
    glass.density = None;
    m.replace_material("Glass", glass).unwrap();
    let diagnostic = m.check().iter().map(ToString::to_string).collect::<Vec<_>>().join("\n");
    assert!(diagnostic.contains("material `Glass`: has no density"), "{diagnostic}");
    m.replace_material("Glass", Material { density: Some(1600.0), ..Material::default() }).unwrap();
    let diagnostic = m.check().iter().map(ToString::to_string).collect::<Vec<_>>().join("\n");
    assert!(diagnostic.contains("Glass") && diagnostic.contains("elastic"), "{diagnostic}");
    assert!(Material { density: Some(1600.0), ply_strength: Some(strength()), ..Material::default() }.validate().is_err());
}

#[test]
fn orientation_validates_its_kind_and_reference_axis() {
    let mut s = section();
    for axis in [0, 4, u8::MAX] {
        s.orientation.as_mut().unwrap().shell_reference_axis = axis;
        assert!(s.validate().is_err());
    }
    s.orientation.as_mut().unwrap().shell_reference_axis = 1;
    s.material = "Carbon".into();
    for kind in [SectionKind::Solid {}, SectionKind::Shell { thickness: 0.1, offset: 0.0, points: 3 }] {
        s.kind = kind;
        s.validate().unwrap();
    }
    s.kind = SectionKind::PlaneStress { thickness: 0.1 };
    assert!(s.validate().is_err());
    s.kind = SectionKind::Solid {};
    s.orientation.as_mut().unwrap().shell_reference_axis = 3;
    assert!(s.validate().is_err());
}

#[test]
fn conversion_scales_strengths_origins_thickness_and_angles_but_not_directions() {
    let m = plate();
    let system = DisplaySystem::new(Preset::MmTonneS).with_override(Quantity::Angle, "deg").unwrap();
    let expressed = m.converted(&|q, v| system.express(q, v).ok()).unwrap();
    let s = expressed.section("laminate").unwrap();
    let orientation = s.orientation.unwrap();
    assert_eq!(orientation.frame, MaterialOrientation::Cylindrical { origin: [100.0, 200.0, 300.0], axis: [0.0, 0.0, 2.0] });
    assert_eq!(orientation.shell_reference_axis, 3);
    let SectionKind::CompositeShell { layup } = &s.kind else { panic!("composite section") };
    assert_eq!(layup.thickness(), 4.0);
    assert!((layup.plies[0].angle - 45.0).abs() < 1e-12);
    assert_eq!(layup.offset, 0.25);
    assert_eq!(layup.plies[1].points, 1);
    let strength = expressed.material("Carbon").unwrap().ply_strength.unwrap();
    assert_eq!(strength.tension, [1500.0, 50.0]);
    assert_eq!(strength.compression, [1200.0, 150.0]);
    assert_eq!(strength.shear12, 80.0);
    assert_eq!(strength.shear23, 60.0);
    assert_eq!(strength.tsai_wu_interaction, -0.25);
    let back = expressed.from_converted(&expressed, &m, &|q, v| system.to_si(q, v).ok()).unwrap();
    assert_eq!(back, m);
}
