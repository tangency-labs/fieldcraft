use super::*;
use fieldcraft_model::{ElementKind as K, Integration};

#[test]
fn constant_property_collection_refuses_excess_at_its_actual_data_line() {
    let input = format!("*MATERIAL,NAME=M\n*ELASTIC\n{}\n{}\nBAD\n", "1,".repeat(16), "1,".repeat(16));
    let message = read(&input).unwrap_err().to_string();
    assert!(message.contains("line 5") && message.contains("32 numeric entries"), "{message}");
}

#[test]
fn malformed_history_never_leaks_data_back_into_model_scope() {
    for (text, line) in [("*STEP\n*STEP\n*END STEP\n*NODE\n99,1,2,3\n*END STEP\n", 2), ("*END STEP\n", 1), ("*STEP\n*STATIC\n1,1\n", 1)] {
        let message = read(text).unwrap_err().to_string();
        assert!(message.contains(&format!("line {line}")), "{message}");
    }
    let input = "*STEP\n*NODE\n99,1,2,3\n*END STEP\n*NODE\n1,0,0,0\n";
    let got = read(input).unwrap();
    assert!(got.model.node(99).is_none());
    assert!(got.model.node(1).is_some());
}

#[test]
fn pure_reader_empty_input_is_explicitly_empty_for_host_replacement_policy() {
    for text in ["", "** comment\n", "*HEADING\nAn original model title\n", "*UNKNOWN\nignored\n"] {
        let got = read(text).unwrap();
        assert!(got.model.nodes().is_empty() && got.model.elements().is_empty() && got.model.materials().is_empty());
    }
}

#[test]
fn scoped_cancellation_stops_the_lexer_and_periodic_member_expansion() {
    let progress = progress::Progress::new();
    let mut raw = RawSet::default();
    let mut budget = 0;
    // Some work has completed before cancellation; the next periodic range-expansion check
    // refuses to append another member. No timing or extra cancellation API is involved.
    for label in 1..256 {
        append_member(&mut raw, Member::Label(label, 1), &mut budget, Limits::default().members, 1).unwrap();
    }
    progress.cancel();
    progress::scoped(Some(progress), || {
        assert!(matches!(read(TET), Err(crate::Error::Model(fieldcraft_model::Error::Cancelled))));
        assert!(matches!(
            append_member(&mut raw, Member::Label(256, 1), &mut budget, Limits::default().members, 1),
            Err(crate::Error::Model(fieldcraft_model::Error::Cancelled))
        ));
    });
    assert_eq!(raw.members.len(), 255, "the cancelled range append did not allocate another member");
    assert!(read(TET).is_ok(), "scoped cancellation does not leak into another caller");
}

#[test]
fn global_system_cannot_make_imported_nodes_silently_use_local_coordinates() {
    let input = "*SYSTEM\n10,20,30\n*NODE\n1,0,0,0\n";
    let message = read(input).unwrap_err().to_string();
    assert!(message.contains("line 1") && message.contains("SYSTEM"), "{message}");
    assert_eq!(read("*NODE,SYSTEM=R\n1,,,\n").unwrap().model.node(1).unwrap().xyz, [0., 0., 0.]);
}

#[test]
fn bare_valued_parameters_and_valued_flags_never_fall_back_to_defaults() {
    for text in [
        "*MATERIAL,NAME=M\n*ELASTIC,TYPE\n2,0.3",
        "*NODE,NSET\n1,0,0,0",
        "*ELEMENT,TYPE=C3D4,ELSET\n1,1,2,3,4",
        "*NSET,NSET=A,GENERATE=YES\n1,4",
        "*NSET,NSET=A,GENERATE=\n1,4",
    ] {
        let message = read(text).unwrap_err().to_string();
        assert!(message.contains("parameter") || message.contains("flag"), "{message}");
    }
}

#[test]
fn unread_material_behaviors_preserve_the_current_material_and_physical_lines() {
    let text = "*MATERIAL,NAME=M\n*CONDUCTIVITY\n40\n*ELASTIC\n2,0.3\n*PLASTIC\n1,0\n*DENSITY\n1\n*MATERIAL,NAME=N\n*SPECIFIC HEAT\n5\n*DENSITY\n2\n";
    let got = read(text).unwrap();
    assert_eq!(got.model.material("M").unwrap().density, Some(1.));
    assert_eq!(got.model.material("M").unwrap().elastic, Some(Elastic::Isotropic { young: 2., poisson: 0.3 }));
    assert!(got.model.material("M").unwrap().plastic.is_none());
    assert_eq!(got.model.material("N").unwrap().density, Some(2.));
    assert!(got.model.material("N").unwrap().elastic.is_none());
    assert_eq!(
        got.unknown,
        vec![
            Unknown { keyword: "CONDUCTIVITY".into(), line: 2 },
            Unknown { keyword: "PLASTIC".into(), line: 6 },
            Unknown { keyword: "SPECIFICHEAT".into(), line: 11 }
        ]
    );
    assert_eq!(got.unknown_count, 3);
}

#[test]
fn nonmaterial_and_unknown_options_end_material_context() {
    for block in ["*NODE\n1,0,0,0", "*NSET,NSET=A\n1", "*CUSTOM\nignored", "*STEP\n*END STEP"] {
        let text = format!("*MATERIAL,NAME=M\n*ELASTIC\n2,0.3\n{block}\n*DENSITY\n1\n");
        let message = read(&text).unwrap_err().to_string();
        assert!(message.contains("outside a supported *MATERIAL block"), "{message}");
    }
}

#[test]
fn own_hyperelastic_writer_decks_import_supported_data_and_report_unread_physics() {
    use fieldcraft_model::{Hyperelastic, NonlinearStatic, Procedure, Step};
    for hyperelastic in [Hyperelastic::NeoHookean {}, Hyperelastic::MooneyRivlin { c01: 0.1 }, Hyperelastic::Yeoh { c20: 0.1, c30: 0.2 }] {
        let mut source = read_in(TET, Preset::MmTonneS).unwrap().model;
        let material = source.material("Steel").unwrap().clone();
        source.replace_material("Steel", Material { hyperelastic: Some(hyperelastic), ..material }).unwrap();
        source.add_step(Step::new("Load", Procedure::NonlinearStatic(NonlinearStatic::default())), None).unwrap();
        let written = super::super::write(&source, Preset::MmTonneS).unwrap();
        assert!(written.contains("*HYPERELASTIC") && !written.contains("*ELASTIC"));
        let got = read_in(&written, Preset::MmTonneS).unwrap();
        assert_eq!(got.model.nodes(), source.nodes());
        assert_eq!(got.model.elements(), source.elements());
        let imported = got.model.material("Steel").unwrap();
        assert_eq!(imported.density, source.material("Steel").unwrap().density);
        assert!(imported.elastic.is_none() && imported.hyperelastic.is_none());
        let line = written.lines().position(|line| line.starts_with("*HYPERELASTIC")).unwrap() + 1;
        assert!(got.unknown.contains(&Unknown { keyword: "HYPERELASTIC".into(), line }));
        assert!(got.model.steps().is_empty());
    }
}

#[test]
fn same_name_node_element_and_surface_sets_combine_without_losing_members() {
    let got = read(&format!("{TET}*NSET,NSET=Combined\n1,4\n*ELSET,ELSET=Combined\n10\n*SURFACE,NAME=Combined,TYPE=ELEMENT\n10,S1\n")).unwrap();
    let set = got.model.set("Combined").unwrap();
    assert_eq!(set.nodes, vec![1, 4]);
    assert_eq!(set.elements, vec![10]);
    assert_eq!(set.faces, vec![ElementFace { element: 10, face: 1 }]);
}

#[test]
fn generate_zero_increment_is_named_at_the_data_line() {
    for kind in ["NSET", "ELSET"] {
        let message = read(&format!("*{kind},{kind}=A,GENERATE\n1,4,0\n")).unwrap_err().to_string();
        assert!(message.contains("line 2") && message.contains("increment") && message.contains("zero"), "{message}");
    }
}

#[test]
fn all_writer_beam_profiles_and_engineering_elasticity_retain_their_dimensions() {
    let si = DisplaySystem::new(Preset::MKgS);
    let length = |s| Quantity::Length.parse(s, &si).unwrap();
    let cases = [
        ("RECT", "40,60", BeamProfile::Rectangle { width: length("40 mm"), height: length("60 mm") }),
        ("CIRC", "20", BeamProfile::Circle { radius: length("20 mm") }),
        ("PIPE", "20,2", BeamProfile::Pipe { radius: length("20 mm"), thickness: length("2 mm") }),
        ("I", "50,100,50,50,8,8,5", BeamProfile::I { height: length("100 mm"), width: length("50 mm"), web: length("5 mm"), flange: length("8 mm") }),
        ("BOX", "40,60,2,2,2,2", BeamProfile::Box { width: length("40 mm"), height: length("60 mm"), thickness: length("2 mm") }),
    ];
    for (shape, data, profile) in cases {
        let source = format!(
            "*NODE\n1,0,0,0\n2,1000,0,0\n*ELEMENT,TYPE=B31,ELSET=Beam\n1,1,2\n*MATERIAL,NAME=Steel\n*ELASTIC,TYPE=ENGINEERING CONSTANTS\n140000,10000,10000,0.3,0.3,0.4,5000,5000\n3800\n*BEAM SECTION,ELSET=Beam,MATERIAL=Steel,SECTION={shape}\n{data}\n0,1,0\n"
        );
        let mut got = read_in(&source, Preset::MmTonneS).unwrap();
        assert_eq!(got.model.sections().values().next().unwrap().kind, SectionKind::Beam { profile, orientation: [0., 1., 0.] });
        assert_eq!(got.model.node(2).unwrap().xyz, [length("1000 mm"), 0., 0.]);
        let material = got.model.material("Steel").unwrap();
        if let Some(Elastic::EngineeringConstants { young, poisson, shear }) = material.elastic {
            assert_eq!(
                young,
                [
                    Quantity::Pressure.parse("140 GPa", &si).unwrap(),
                    Quantity::Pressure.parse("10 GPa", &si).unwrap(),
                    Quantity::Pressure.parse("10 GPa", &si).unwrap()
                ]
            );
            assert_eq!(poisson, [0.3, 0.3, 0.4]);
            assert_eq!(shear[2], Quantity::Pressure.parse("3.8 GPa", &si).unwrap());
        } else {
            panic!("engineering elasticity missing");
        }
        // Orthotropic beam export is deliberately refused by the existing writer; use its
        // supported isotropic material to check geometric roundtrip independently.
        got.model
            .replace_material(
                "Steel",
                Material {
                    elastic: Some(Elastic::Isotropic { young: Quantity::Pressure.parse("210 GPa", &si).unwrap(), poisson: 0.3 }),
                    ..Material::default()
                },
            )
            .unwrap();
        got.model.add_step(fieldcraft_model::Step::new("Load", fieldcraft_model::Procedure::Static {}), None).unwrap();
        let written = super::super::write(&got.model, Preset::MmTonneS).unwrap();
        let again = read_in(&written, Preset::MmTonneS).unwrap();
        assert_eq!(got.model.sections(), again.model.sections());
    }
}

const TET: &str = "*NODE\n1,0,0,0\n2,1,0,0\n3,0,1,0\n4,0,0,1\n*ELEMENT,TYPE=C3D4,ELSET=All\n10,1,2,3,4\n*NSET,NSET=Root\n1,3,4\n*SURFACE,NAME=Bottom,TYPE=ELEMENT\nAll,S1\n*MATERIAL,NAME=Steel\n*ELASTIC\n210000,0.3\n*DENSITY\n7.85D-9\n*SOLID SECTION,ELSET=All,MATERIAL=Steel\n,\n";

#[test]
fn original_tetra_preserves_labels_surfaces_and_explicit_units() {
    let mut read = read_in(TET, Preset::MmTonneS).unwrap();
    assert_eq!(read.model.element(10).unwrap().nodes, vec![1, 2, 3, 4]);
    let si = DisplaySystem::new(Preset::MKgS);
    assert_eq!(read.model.node(2).unwrap().xyz, [Quantity::Length.parse("1 mm", &si).unwrap(), 0., 0.]);
    assert_eq!(read.model.set("Bottom").unwrap().faces, vec![ElementFace { element: 10, face: 1 }]);
    assert_eq!(read.model.material("Steel").unwrap().density, Some(Quantity::Density.parse("7.85e-9 t/mm^3", &si).unwrap()));
    assert_eq!(read.unknown_count, 0);
    read.model.add_step(fieldcraft_model::Step::new("Load", fieldcraft_model::Procedure::Static {}), None).unwrap();
    let written = super::super::write(&read.model, Preset::MmTonneS).unwrap();
    let mut again = read_in(&written, Preset::MmTonneS).unwrap();
    again.model.add_step(fieldcraft_model::Step::new("Load", fieldcraft_model::Procedure::Static {}), None).unwrap();
    assert_eq!(read.model.nodes(), again.model.nodes());
    assert_eq!(read.model.elements(), again.model.elements());
    assert_eq!(written, super::super::write(&again.model, Preset::MmTonneS).unwrap());
}

#[test]
fn supported_writer_table_has_canonical_inverse_formulations() {
    for (name, kind) in [("C3D4", K::Tet4), ("C3D10", K::Tet10), ("C3D6", K::Wedge6), ("C3D5", K::Pyramid5), ("C3D8", K::Hex8), ("C3D20", K::Hex20)] {
        for suffix in ["", "H"] {
            let t = Type::parse(&format!("{name}{suffix}"), 2).unwrap();
            assert_eq!(t.kind, kind);
            assert_eq!(t.formulation.hybrid, suffix == "H");
            assert_eq!(t.formulation.integration, Integration::Full);
        }
    }
    for name in ["C3D8R", "C3D8RH", "C3D20R", "C3D20RH", "S3R", "S4R", "S8R"] {
        assert_eq!(Type::parse(name, 1).unwrap().formulation.integration, Integration::Reduced);
    }
    for name in ["C3D8I", "C3D8IH"] {
        assert!(Type::parse(name, 1).unwrap().formulation.incompatible_modes);
    }
    for prefix in ["CPS", "CPE", "CAX"] {
        for suffix in ["3", "6", "4", "4R", "4I", "8", "8R"] {
            assert!(Type::parse(&format!("{prefix}{suffix}"), 1).is_ok(), "{prefix}{suffix}");
        }
    }
    for (name, kind) in
        [("S3", K::Tri3), ("S4", K::Quad4), ("B31", K::Beam2), ("T3D2", K::Rod2), ("T3D3", K::Rod3), ("T2D2", K::Rod2), ("T2D3", K::Rod3)]
    {
        assert_eq!(Type::parse(name, 1).unwrap().kind, kind);
    }
    for name in ["C3D4R", "C3D10M", "S8", "S6", "B32", "CPS4H", "T3D2H"] {
        assert!(Type::parse(name, 1).is_err(), "{name}");
    }
}

#[test]
fn connectivity_and_keywords_continue_without_losing_physical_lines() {
    let mut deck = String::from("*NODE\n");
    for i in 1..=20 {
        deck.push_str(&format!("{i},{i},0,0\n"));
    }
    deck.push_str("*ELEMENT,\n TYPE=C3D20, ELSET=Long\n99,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,\n** original comment\n\n16,17,18,19,20\n");
    let got = read(&deck).unwrap();
    assert_eq!(got.model.element(99).unwrap().nodes, (1..=20).collect::<Vec<_>>());
    let bad = deck.replace("16,17,18,19,20", "16,17");
    assert!(read(&bad).unwrap_err().to_string().contains("line 27"));
}

#[test]
fn case_insensitive_references_quoted_identifiers_and_finite_defaults() {
    let text = TET.replace("ELSET=All", "ELSET=\"All\"").replace("All,S1", "aLL,S1").replace("ELSET=All,MATERIAL=Steel", "ELSET=all,MATERIAL=sTEel");
    let got = read(&text).unwrap();
    assert!(got.model.set("All").is_some());
    assert!(read(&TET.replace("NAME=Steel", "NAME=\"Steel Space\"")).is_err());
    assert!(read("*NSET,NSET=Root\n1\n*NSET,NSET=root\n1").unwrap_err().to_string().contains("differ only in case"));
    assert_eq!(read("*NODE\n1,2,,\n").unwrap().model.node(1).unwrap().xyz, [2., 0., 0.]);
}

#[test]
fn nested_sets_generate_and_node_from_element_set_are_bounded() {
    let deck = format!(
        "{TET}*ELSET,ELSET=Pair\nAll\n*NSET,NSET=Generated,GENERATE\n1,4\n*NSET,NSET=Nested\nGenerated,Root\n*NSET,NSET=FromElements,ELSET=Pair\n"
    );
    let got = read(&deck).unwrap();
    assert_eq!(got.model.set("Nested").unwrap().nodes, vec![1, 2, 3, 4]);
    assert_eq!(got.model.set("FromElements").unwrap().nodes, vec![1, 2, 3, 4]);
    for bad in [
        "*NSET,NSET=A\nA",
        "*NSET,NSET=A\nB\n*NSET,NSET=B\nA",
        "*NSET,NSET=A,GENERATE\n1,99999999,1",
        "*NSET,NSET=A,GENERATE\n1,4,0",
        "*NSET,NSET=A,GENERATE\n4,1,1",
    ] {
        assert!(read(bad).is_err(), "{bad}");
    }
    assert!(parse(&deck, Preset::MKgS, Limits { members: 8, ..Limits::default() }).is_err());
}

#[test]
fn unknown_history_options_all_keep_their_original_lines() {
    let deck = "*NODE\n1,0,0,0\n*UNKNOWN,FOO=BAR\ndata\n*STEP\n*STATIC\n1,1\n*NODE\n99,1,2,3\n*CLOAD\n1,1,4\n*END STEP\n";
    let got = read(deck).unwrap();
    assert!(got.model.node(99).is_none());
    assert_eq!(
        got.unknown.iter().map(|u| (u.keyword.as_str(), u.line)).collect::<Vec<_>>(),
        vec![("UNKNOWN", 3), ("STEP", 5), ("STATIC", 6), ("NODE", 8), ("CLOAD", 10), ("ENDSTEP", 12)]
    );
    let many = "*CUSTOM\n".repeat(120);
    let got = read(&many).unwrap();
    assert_eq!(got.unknown_count, 120);
    assert_eq!(got.unknown.len(), MAX_UNKNOWN_SHOWN);
}

#[test]
fn named_refusals_cover_unsupported_model_data_and_missing_references() {
    for bad in [
        "*INCLUDE,INPUT=x.inp",
        "*PART,NAME=P",
        "*ASSEMBLY,NAME=A",
        "*INSTANCE,NAME=I",
        "*ELEMENT,TYPE=C3D10M\n1,1",
        "*ELASTIC\n2,0.3",
        "*DENSITY\n1",
        "*NODE,SYSTEM=C\n1,1,0,0",
        "*NODE,INPUT=x.inp",
        "*NSET,NSET=A\n900",
        "*ELSET,ELSET=A\n900",
    ] {
        let e = read(bad).unwrap_err().to_string();
        assert!(e.contains("line"), "{e}");
    }
    assert!(read(&TET.replace("All,S1", "All,S6")).is_err());
    assert!(read(&TET.replace("MATERIAL=Steel", "MATERIAL=Absent")).is_err());
    assert!(read(&TET.replace("210000,0.3", "210000,0.3,300")).is_err());
    let mixed = TET
        .replace("*ELEMENT,TYPE=C3D4", "*NODE\n5,1,1,0\n6,1,0,1\n7,0,1,1\n8,1,1,1\n*ELEMENT,TYPE=C3D4")
        .replace("*NSET,NSET=Root", "*ELEMENT,TYPE=C3D8R,ELSET=All\n11,1,2,3,4,5,6,7,8\n*NSET,NSET=Root");
    assert!(read(&mixed).unwrap_err().to_string().contains("mixes element types/formulations"));
}

#[test]
fn hostile_numbers_continuations_and_every_prefix_never_panic() {
    for bad in [
        "*NODE\n0,0,0,0",
        "*NODE\n-1,0,0,0",
        "*NODE\n1000000000,0,0,0",
        "*NODE\n1,NaN,0,0",
        "*NODE\n1,1e400,0,0",
        "*NODE\n1,0,0,0\n1,0,0,0",
        "*ELEMENT,TYPE=C3D4\n1,1,2,",
        "*NODE,\n*NSET,NSET=A",
    ] {
        assert!(read(bad).is_err(), "{bad}");
    }
    assert!(read(&"x".repeat(1024 * 1024)).is_err());
    assert!(read(&format!("*NODE,\n{}NSET=A", format!("{},\n", ",".repeat(250)).repeat(20))).is_err());
    for at in 0..=TET.len() {
        assert!(std::panic::catch_unwind(|| read(&TET[..at])).is_ok(), "prefix {at}");
    }
}

#[test]
fn caller_limits_refuse_nodes_elements_and_members_at_their_line() {
    let limited = |nodes, elements, members| read_limited(TET, Preset::MmTonneS, Limits { nodes, elements, members });
    assert!(limited(4, 1, 16).is_ok());
    let message = limited(3, 1, 16).unwrap_err().to_string();
    assert!(message.contains("line 5") && message.contains("node limit of 3"), "{message}");
    let message = limited(4, 0, 16).unwrap_err().to_string();
    assert!(message.contains("line 7") && message.contains("element limit of 0"), "{message}");
    let message = read_limited("*NSET,NSET=A,GENERATE\n1,1000\n", Preset::MKgS, Limits { members: 10, ..Limits::default() }).unwrap_err().to_string();
    assert!(message.contains("line 2") && message.contains("member limit"), "{message}");
}
