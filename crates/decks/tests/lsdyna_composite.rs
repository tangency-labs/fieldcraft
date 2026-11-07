//! Re-read our own decks against the R14 card definitions. No external solver is run.
mod common;
use fieldcraft_decks::lsdyna;
use fieldcraft_model::units::{DisplaySystem, Preset, Quantity as Q};
use fieldcraft_model::{ElementKind, MaterialOrientation, Model, OutputVariable, SectionKind};

fn model() -> Model {
    let mut m = common::laminate_model();
    let mut s = m.section("laminate").unwrap().clone();
    if let SectionKind::CompositeShell { layup } = &mut s.kind {
        for p in &mut layup.plies {
            p.points = 3;
        }
    }
    m.replace_section("laminate", s).unwrap();
    let mut step = m.steps()[0].clone();
    step.output.insert(OutputVariable::Strain);
    m.replace_step(step).unwrap();
    m
}

fn blocks(text: &str, keyword: &str) -> Vec<Vec<Vec<String>>> {
    let mut out = Vec::<Vec<Vec<String>>>::new();
    let mut active = false;
    for line in text.lines().filter(|l| !l.starts_with('$')) {
        if line.starts_with('*') {
            active = line == keyword;
            if active {
                out.push(Vec::new());
            }
        } else if active {
            out.last_mut().unwrap().push(line.as_bytes().chunks(20).map(|c| String::from_utf8(c.to_vec()).unwrap().trim().to_owned()).collect());
        }
    }
    out
}
fn real(fields: &[String], i: usize) -> f64 {
    fields[i].parse().unwrap()
}

#[test]
fn laminate_cards_preserve_axes_stack_units_and_exact_linear_thickness_moments() {
    let m = model();
    let original = serde_json::to_value(&m).unwrap();
    for units in [Preset::MKgS, Preset::MmTonneS, Preset::MmKgMs, Preset::InLbfS] {
        let deck = lsdyna::write(&m, units).unwrap().text;
        assert_eq!(deck, lsdyna::write(&m, units).unwrap().text);
        let d = DisplaySystem::new(units);
        let parts = blocks(&deck, "*PART_COMPOSITE_LONG");
        assert_eq!(parts.len(), 1);
        let p = &parts[0];
        assert_eq!(p[1], ["OPTCARD", "103"]);
        assert_eq!(real(&p[2], 2), 5.0 / 6.0);
        assert_eq!(real(&p[2], 3), 0.5);
        assert_eq!(p.len(), 5);
        assert!((real(&p[3], 2) - 30.0).abs() < 1e-12);
        assert!((real(&p[4], 2) + 45.0).abs() < 1e-12);
        let thickness: Vec<_> = p[3..].iter().map(|p| d.to_si(Q::Length, real(p, 1)).unwrap()).collect();
        assert!((thickness[0] - 0.001).abs() < 1e-15 && (thickness[1] - 0.002).abs() < 1e-15);
        // IRPL=103 evaluates each actual ply at bottom/middle/top, with Simpson weights.
        // Its zeroth/first/second moments must equal the exact linear laminate moments.
        let h: f64 = thickness.iter().sum();
        let mut bottom = -(0.5 + real(&p[2], 3) / 2.0) * h;
        for t in thickness {
            let top = bottom + t;
            for k in 0_i32..=2 {
                let moment = t * (bottom.powi(k) + 4.0 * ((bottom + top) / 2.0).powi(k) + top.powi(k)) / 6.0;
                let exact = (top.powi(k + 1) - bottom.powi(k + 1)) / f64::from(k + 1);
                assert!((moment - exact).abs() < 1e-12 * exact.abs().max(t.powi(k + 1)));
            }
            bottom = top;
        }
        let mats = blocks(&deck, "*MAT_ORTHOTROPIC_THERMAL_TITLE");
        assert_eq!(mats.len(), 1);
        let mat = &mats[0];
        assert_eq!(mat[1][0], p[3][0]);
        assert_eq!(p[3][0], p[4][0]);
        assert!((real(&mat[1], 5) - 0.25 * 9.0 / 120.0).abs() < 1e-15);
        for (i, g) in [5e9, 3e9, 4e9].into_iter().enumerate() {
            assert!((d.to_si(Q::Pressure, real(&mat[2], i)).unwrap() / g - 1.0).abs() < 1e-12);
        }
        assert_eq!(mat[2][3..].iter().map(|v| v.parse::<f64>().unwrap()).collect::<Vec<_>>(), vec![0., 0., 0., 2., 1.]);
        assert_eq!(mat[3][3..].iter().map(|v| v.parse::<f64>().unwrap()).collect::<Vec<_>>(), vec![0., 1., 0.]);
        let output = blocks(&deck, "*DATABASE_EXTENT_BINARY");
        assert_eq!(real(&output[0][0], 2), 6.0);
        assert_eq!(real(&output[0][1], 0), 1.0);
    }
    assert_eq!(serde_json::to_value(&m).unwrap(), original, "export must not mutate the model");
}

#[test]
fn default_and_cylindrical_shell_frames_are_explicit_and_solid_frames_are_not_flattened() {
    let mut m = model();
    let mut s = m.section("laminate").unwrap().clone();
    s.orientation = None;
    m.replace_section("laminate", s.clone()).unwrap();
    let text = lsdyna::write(&m, Preset::MKgS).unwrap().text;
    let mats = blocks(&text, "*MAT_ORTHOTROPIC_THERMAL_TITLE");
    assert!((real(&mats[0][3], 4) / real(&mats[0][3], 3) - 0.4).abs() < 1e-12);
    s.orientation = Some(fieldcraft_model::SectionOrientation {
        frame: MaterialOrientation::Cylindrical { origin: [0.; 3], axis: [0., 0., 1.] },
        shell_reference_axis: 2,
    });
    m.replace_section("laminate", s.clone()).unwrap();
    let text = lsdyna::write(&m, Preset::MKgS).unwrap().text;
    let mats = blocks(&text, "*MAT_ORTHOTROPIC_THERMAL_TITLE");
    assert!((real(&mats[0][3], 4) / real(&mats[0][3], 3) + 1.25 / 1.5).abs() < 1e-12);
    m.set_node(4, [0., 0., 1.]).unwrap();
    m.set_element(1, ElementKind::Tet4, vec![1, 2, 3, 4]).unwrap();
    s.kind = SectionKind::Solid {};
    s.orientation.as_mut().unwrap().shell_reference_axis = 1;
    s.material = "ply".into();
    m.replace_section("laminate", s.clone()).unwrap();
    assert!(lsdyna::write(&m, Preset::MKgS).unwrap_err().to_string().contains("integration point"));
    s.orientation.as_mut().unwrap().frame = MaterialOrientation::Rectangular { axis1: [0., 1., 0.], axis2: [-1., 0., 0.] };
    m.replace_section("laminate", s).unwrap();
    let text = lsdyna::write(&m, Preset::MKgS).unwrap().text;
    assert!(text.contains("*SECTION_SOLID"));
    let mats = blocks(&text, "*MAT_ORTHOTROPIC_THERMAL_TITLE");
    assert_eq!(real(&mats[0][3], 4), 1.0);
    assert_eq!(real(&mats[0][4], 3), -1.0);
}

#[test]
fn midpoint_integration_and_density_averaging_are_named_refusals() {
    let mut m = model();
    let mut s = m.section("laminate").unwrap().clone();
    if let SectionKind::CompositeShell { layup } = &mut s.kind {
        layup.plies[0].points = 1;
    }
    m.replace_section("laminate", s).unwrap();
    assert!(lsdyna::write(&m, Preset::MKgS).unwrap_err().to_string().contains("IRPL=103"));
    let mut m = model();
    let mut material = m.material("ply").unwrap().clone();
    material.density = Some(1200.);
    m.add_material("lighter", material).unwrap();
    let mut s = m.section("laminate").unwrap().clone();
    if let SectionKind::CompositeShell { layup } = &mut s.kind {
        layup.plies[1].material = "lighter".into();
    }
    m.replace_section("laminate", s).unwrap();
    assert!(lsdyna::write(&m, Preset::MKgS).unwrap_err().to_string().contains("unequal ply densities"));
}
