//! Original close-frequency rods: public equal-damping CQC with signed components.
use fieldcraft_engine::Session;
use serde_json::json;

const INPUT: &str = "*NODE\n1,0,0,0\n2,-0.7071067811865476,-0.7071067811865476,0\n3,0.7071067811865476,-0.7071067811865476,0\n*ELEMENT,TYPE=T3D2,ELSET=Bar1\n1,2,1\n*ELEMENT,TYPE=T3D2,ELSET=Bar2\n2,3,1\n*MATERIAL,NAME=E1\n*ELASTIC\n192.3769,.25\n*DENSITY\n1.5\n*MATERIAL,NAME=E2\n*ELASTIC\n194.0449,.25\n*DENSITY\n1.5\n*SOLID SECTION,ELSET=Bar1,MATERIAL=E1\n1\n*SOLID SECTION,ELSET=Bar2,MATERIAL=E2\n1\n*BOUNDARY\n2,1,3\n3,1,3\n1,3,3\n*SPECTRUM,NAME=S,TYPE=ACCELERATION\n1,1,.05\n1,5,.05\n*STEP,NAME=Modes,PERTURBATION\n*FREQUENCY,NORMALIZATION=MASS\n2\n*END STEP\n*STEP,NAME=Spectrum,PERTURBATION\n*RESPONSE SPECTRUM,SUM=CQC,COMP=SRSS\nS,1,0,0,1\n*MODAL DAMPING,DEFINITION=MODE NUMBERS\n1,2,.05\n*OUTPUT,FIELD\n*NODE OUTPUT\nU,A\n*ELEMENT OUTPUT\nS\n*END STEP\n";

#[test]
fn cqc_import_solves_and_survives_native_and_deck_roundtrips() {
    let dir = std::env::temp_dir().join(format!("fieldcraft-cqc-import-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut peaks = Vec::new();
    for combination in ["CQC", "SRSS"] {
        let text = INPUT.replace("SUM=CQC", &format!("SUM={combination}"));
        let input = dir.join("input.inp");
        std::fs::write(&input, &text).unwrap();
        let mut s = Session::default();
        let imported = s.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap();
        assert_eq!(imported["analysis_complete"], true, "{imported}");
        assert!(imported["warnings"].as_array().unwrap().iter().any(|w| w.as_str().unwrap().contains("explicit mass=lumped")));
        s.execute("solve.spectrum", &json!({"step":"Spectrum","mass":"lumped","threads":1})).unwrap();
        let probe = s.execute("solve.probe", &json!({"step":"Spectrum","nodes":[1]})).unwrap();
        let r = (192.3769_f64 / 194.0449).sqrt();
        let z = 0.05_f64;
        let rho = if combination == "CQC" {
            8.0 * z * z * (1.0 + r) * r.powf(1.5) / ((1.0 - r * r).powi(2) + 4.0 * z * z * r * (1.0 + r).powi(2))
        } else {
            0.0
        };
        let a = 0.75 / 192.3769;
        let b = 0.75 / 194.0449;
        let actual = probe["nodes"][0]["samples"][0]["displacement"][1].as_f64().unwrap();
        assert!((actual - (a * a + b * b - 2.0 * rho * a * b).sqrt()).abs() < 1e-12);
        peaks.push(actual);
        let saved = dir.join("saved.fieldcraft");
        s.execute("file.save_as", &json!({"path":saved,"results":true})).unwrap();
        s.execute("file.open", &json!({"path":saved})).unwrap();
        assert_eq!(s.execute("solve.probe", &json!({"step":"Spectrum","nodes":[1]})).unwrap(), probe);
        let replay = dir.join(format!("replay-{combination}.inp"));
        s.execute("deck.source.write", &json!({"path":replay})).unwrap();
        assert_eq!(std::fs::read_to_string(replay).unwrap(), text);
        let denied = s
            .execute(
                "deck.write",
                &json!({"format":"abaqus","path":dir.join("consistent.inp"),"steps":["Spectrum"],"mass":"consistent","units":"m-kg-s-N-Pa"}),
            )
            .unwrap_err();
        assert!(denied.to_string().contains("stored mass Lumped differs from export Consistent"), "{denied}");
        let normalized = dir.join("normalized.inp");
        s.execute("deck.write", &json!({"format":"abaqus","path":normalized,"steps":["Spectrum"],"mass":"lumped","units":"m-kg-s-N-Pa"})).unwrap();
        let imported = s.execute("deck.read", &json!({"format":"abaqus","path":normalized,"units":"m-kg-s-N-Pa"})).unwrap();
        assert_eq!(imported["analysis_complete"], true);
        s.execute("solve.spectrum", &json!({"step":"Spectrum","mass":"lumped","threads":1})).unwrap();
        assert_eq!(s.execute("solve.probe", &json!({"step":"Spectrum","nodes":[1]})).unwrap(), probe);
    }
    assert!(peaks[0] < 0.1 * peaks[1], "CQC must not silently use SRSS");
    std::fs::remove_dir_all(dir).unwrap();
}
