//! Original diagonal three-axis oscillator; global principal masses stay distinct.
use fieldcraft_engine::Session;
use serde_json::{Value, json};
const INPUT: &str = r#"** Original global directional point mass, SI
*NODE
7,0,0,0
*NSET,NSET=Point
7
*ELEMENT,TYPE=MASS,ELSET=M
10,7
*MASS,ELSET=M,TYPE=ANISOTROPIC
2,3,5
*ELEMENT,TYPE=SPRING1,ELSET=Kx
11,7
*ELEMENT,TYPE=SPRING1,ELSET=Ky
12,7
*ELEMENT,TYPE=SPRING1,ELSET=Kz
13,7
*SPRING,ELSET=Kx
1
200
*SPRING,ELSET=Ky
2
600
*SPRING,ELSET=Kz
3
1500
*STEP,NAME=Modes
*FREQUENCY
3
*END STEP
*STEP,NAME=ResponseA
*STEADY STATE DYNAMICS,DIRECT
1,0
*CLOAD
Point,1,2
Point,2,3
Point,3,5
*END STEP
*STEP,NAME=ResponseB
*STEADY STATE DYNAMICS,DIRECT
1,0
*CLOAD
Point,1,4
*END STEP
"#;
fn solve(s: &mut Session) -> Vec<Value> {
    let modes = s.execute("solve.frequency", &json!({"step":"Modes","mass":"consistent"})).unwrap();
    for (i, lambda) in [100_f64, 200., 300.].into_iter().enumerate() {
        assert!((modes["modes"][i]["frequency"].as_f64().unwrap() - lambda.sqrt() / (2. * std::f64::consts::PI)).abs() < 1e-12);
    }
    let mut probes = vec![s.execute("solve.probe", &json!({"step":"Modes","nodes":[7]})).unwrap()];
    for (name, force) in [("ResponseA", [2., 3., 5.]), ("ResponseB", [4., 0., 0.])] {
        s.execute("solve.harmonic", &json!({"step":name,"mass":"consistent"})).unwrap();
        let probe = s.execute("solve.probe", &json!({"step":name,"nodes":[7]})).unwrap();
        let sample = &probe["nodes"][0]["samples"][0];
        for (axis, (k, m)) in [200., 600., 1500.].into_iter().zip([2., 3., 5.]).enumerate() {
            let expected = force[axis] / (k - m * (2. * std::f64::consts::PI).powi(2));
            assert!((sample["displacement_real"][axis].as_f64().unwrap() - expected).abs() < 1e-12);
            assert!(sample["displacement_imaginary"][axis].as_f64().unwrap().abs() < 1e-12);
        }
        probes.push(probe);
    }
    probes
}
#[test]
fn directional_mass_solves_reopens_replays_and_roundtrips() {
    let dir = std::env::temp_dir().join(format!("fieldcraft-directional-mass-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("input.inp");
    std::fs::write(&input, INPUT).unwrap();
    let mut s = Session::default();
    let read = s.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(read["analysis_complete"], true, "{read}");
    let probes = solve(&mut s);
    let file = dir.join("saved.fieldcraft");
    s.execute("file.save_as", &json!({"path":file,"results":true})).unwrap();
    s.execute("file.new", &json!({})).unwrap();
    s.execute("file.open", &json!({"path":file})).unwrap();
    assert_eq!(s.model.masses().values().next().unwrap().mass, Some(fieldcraft_model::Mass::Directional([2., 3., 5.])));
    for (step, expected) in ["Modes", "ResponseA", "ResponseB"].into_iter().zip(&probes) {
        assert_eq!(&s.execute("solve.probe", &json!({"step":step,"nodes":[7]})).unwrap(), expected);
    }
    let replay = dir.join("replay.inp");
    s.execute("deck.source.write", &json!({"path":replay})).unwrap();
    assert_eq!(std::fs::read_to_string(replay).unwrap(), INPUT);
    let output = dir.join("normalized.inp");
    s.execute("deck.write", &json!({"format":"abaqus","path":output,"units":"m-kg-s-N-Pa"})).unwrap();
    assert!(std::fs::read_to_string(&output).unwrap().contains("TYPE=ANISOTROPIC"));
    let read = s.execute("deck.read", &json!({"format":"abaqus","path":output,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(read["analysis_complete"], true, "{read}");
    assert_eq!(solve(&mut s), probes);
    std::fs::remove_dir_all(dir).unwrap();
}
