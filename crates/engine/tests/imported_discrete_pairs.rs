//! Independent two-mass chain, imported through the production commands.
use fieldcraft_engine::Session;
use serde_json::{Value, json};

const INPUT: &str = "*NODE\n1,0,0,0\n2,1,0,0\n3,2,0,0\n*ELEMENT,TYPE=MASS,ELSET=M2\n10,2\n*ELEMENT,TYPE=MASS,ELSET=M3\n11,3\n*MASS,ELSET=M2\n2\n*MASS,ELSET=M3\n3\n*ELEMENT,TYPE=SPRING2,ELSET=K1\n20,1,2\n*ELEMENT,TYPE=SPRING2,ELSET=K2\n21,2,3\n*SPRING,ELSET=K1\n1,1\n100\n*SPRING,ELSET=K2\n1,1\n200\n*BOUNDARY\n1,1,3\n2,2,3\n3,2,3\n*STEP,NAME=Static\n*STATIC\n1,1\n*CLOAD\n3,1,6\n*END STEP\n*STEP,NAME=Modes\n*FREQUENCY\n2\n*END STEP\n";

fn solve(s: &mut Session) -> Value {
    s.execute("solve.static", &json!({"step":"Static","threads":1})).unwrap();
    let probe = s.execute("solve.probe", &json!({"step":"Static","nodes":[1,2,3]})).unwrap();
    for (i, u) in [0., 0.06, 0.09].iter().enumerate() {
        assert!((probe["nodes"][i]["displacement"][0].as_f64().unwrap() - u).abs() < 1e-12);
    }
    let modes = s.execute("solve.frequency", &json!({"step":"Modes","mass":"consistent","threads":1})).unwrap();
    // det(K-lambda M)=6 lambda^2-1300 lambda+20000, roots 50/3 and 200.
    for (i, lambda) in [50. / 3., 200.].iter().enumerate() {
        let f = modes["modes"][i]["frequency"].as_f64().unwrap();
        assert!((f - f64::sqrt(*lambda) / std::f64::consts::TAU).abs() < 1e-10);
    }
    probe
}

#[test]
fn spring2_chain_solves_reopens_and_roundtrips() {
    let dir = std::env::temp_dir().join(format!("fieldcraft-import-pairs-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("input.inp");
    std::fs::write(&input, INPUT).unwrap();
    let mut s = Session::default();
    let r = s.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(r["analysis_complete"], true, "{r}");
    s.execute("model.units", &json!({"system":"m-kg-s-N-Pa"})).unwrap();
    let actual = solve(&mut s);
    let saved = dir.join("saved.fieldcraft");
    s.execute("file.save_as", &json!({"path":saved,"results":true})).unwrap();
    s.execute("file.open", &json!({"path":saved})).unwrap();
    assert_eq!(s.execute("solve.probe", &json!({"step":"Static","nodes":[1,2,3]})).unwrap(), actual);
    let replay = dir.join("replay.inp");
    s.execute("deck.source.write", &json!({"path":replay})).unwrap();
    assert_eq!(std::fs::read_to_string(replay).unwrap(), INPUT);
    let normal = dir.join("normalized.inp");
    s.execute("deck.write", &json!({"format":"abaqus","path":normal,"units":"m-kg-s-N-Pa"})).unwrap();
    s.execute("deck.read", &json!({"format":"abaqus","path":normal,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(solve(&mut s), actual);
    std::fs::remove_dir_all(dir).unwrap();
}
