//! Original tetrahedron on a normal elastic bed. No external solver data.
use fieldcraft_engine::Session;
use serde_json::{Value, json};

const INPUT: &str = "*NODE\n1,0,0,0\n2,1,0,0\n3,0,1,0\n4,.3333333333333333,.3333333333333333,1\n*ELEMENT,TYPE=C3D4,ELSET=Body\n1,1,2,3,4\n*NSET,NSET=All\n1,2,3,4\n*MATERIAL,NAME=M\n*ELASTIC\n1000,0\n*SOLID SECTION,ELSET=Body,MATERIAL=M\n*FOUNDATION\nBody,F1,1000\n*BOUNDARY\nAll,1,2\n*STEP,NAME=Load\n*STATIC\n1,1\n*CLOAD\n4,3,-1\n*END STEP\n*STEP,NAME=Double\n*STATIC\n1,1\n*CLOAD,OP=NEW\n4,3,-2\n*END STEP\n";

fn solve(s: &mut Session) -> Vec<Value> {
    ["Load", "Double"]
        .iter()
        .enumerate()
        .map(|(i, step)| {
            let r = s.execute("solve.static", &json!({"step":step,"threads":1})).unwrap();
            // Area=1/2, height=1, nu=0: bed translation F/(k A), body shortening 3Fh/(EA).
            assert!((r["max_displacement"]["value"].as_f64().unwrap() - 0.008 * (i + 1) as f64).abs() < 1e-12);
            s.execute("solve.probe", &json!({"step":step,"nodes":[1,2,3,4]})).unwrap()
        })
        .collect()
}

#[test]
fn foundation_is_solved_persisted_and_written_as_model_data() {
    let dir = std::env::temp_dir().join(format!("fieldcraft-foundation-import-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("input.inp");
    std::fs::write(&input, INPUT).unwrap();
    let mut s = Session::default();
    let read = s.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(read["analysis_complete"], true, "{read}");
    let expected = solve(&mut s);
    let saved = dir.join("saved.fieldcraft");
    s.execute("file.save_as", &json!({"path":saved,"results":true})).unwrap();
    s.execute("file.open", &json!({"path":saved})).unwrap();
    for (i, step) in ["Load", "Double"].iter().enumerate() {
        assert_eq!(s.execute("solve.probe", &json!({"step":step,"nodes":[1,2,3,4]})).unwrap(), expected[i]);
    }
    let replay = dir.join("replay.inp");
    s.execute("deck.source.write", &json!({"path":replay})).unwrap();
    assert_eq!(std::fs::read_to_string(replay).unwrap(), INPUT);
    let normal = dir.join("normalized.inp");
    s.execute("deck.write", &json!({"format":"abaqus","path":normal,"units":"m-kg-s-N-Pa"})).unwrap();
    s.execute("deck.read", &json!({"format":"abaqus","path":normal,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(solve(&mut s), expected);
    assert_eq!(s.model.foundations().values().next().unwrap().normal, 1000.);
    // Input identifiers are case-insensitive even when the generated native names are mixed case.
    let extra = INPUT.replace("*FOUNDATION", "*NSET,NSET=PLACEHOLDER\n1\n*FOUNDATION");
    let row = extra.lines().position(|l| l == "Body,F1,1000").unwrap() + 1;
    let collision = extra.replace("PLACEHOLDER", &format!("INPFOUNDATIONFACES{row}_0"));
    std::fs::write(&input, collision).unwrap();
    s.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap();
    assert!(s.model.foundations().values().next().unwrap().set.ends_with("_1"));
    assert_eq!(solve(&mut s), expected);
    for (old, new, reason) in [
        ("*FOUNDATION\n", "*FOUNDATION,AMPLITUDE=Ramp\n", "AMPLITUDE"),
        ("Body,F1,1000", "Body,F,1000", "FOUNDATION type F"),
        ("*STEP,NAME=Load", "*STEP,NAME=Load,NLGEOM=YES", "FOUNDATION requires"),
    ] {
        std::fs::write(&input, INPUT.replace(old, new)).unwrap();
        let r = s.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap();
        assert_eq!(r["analysis_complete"], false, "{r}");
        assert!(r.to_string().contains(reason), "{r}");
        assert!(s.execute("solve.static", &json!({"step":"Load"})).is_err());
    }
    std::fs::remove_dir_all(dir).unwrap();
}
