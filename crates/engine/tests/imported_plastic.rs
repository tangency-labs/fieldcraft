//! Original single tetrahedron in uniaxial stress beyond first yield.
use fieldcraft_engine::Session;
use serde_json::{Value, json};
const INPUT: &str = r#"** Original constant J2 hardening tetrahedron; SI, small strain
*NODE
1,0,0,0
2,1,0,0
3,0,1,0
4,0,0,1
*ELEMENT,TYPE=C3D4,ELSET=Body
1,1,2,3,4
*NSET,NSET=All
1,2,3,4
*MATERIAL,NAME=M
*ELASTIC
1000,.25
*PLASTIC,HARDENING=ISOTROPIC
1,0
2,.01
*SOLID SECTION,ELSET=Body,MATERIAL=M
*STEP,NAME=Pull,NLGEOM=NO,INC=100
*STATIC
.1,1,.00001,.1
*BOUNDARY
1,1,3
2,2,3
3,1,1
3,3,3
4,1,2
*CLOAD
2,1,.25
*NODE PRINT,NSET=All
U,RF
*EL PRINT,ELSET=Body
S,PEEQ
*END STEP
"#;
fn solve(s: &mut Session) -> Value {
    let result = s.execute("solve.nonlinear", &json!({"step":"Pull"})).unwrap();
    assert!((result["max_plastic_strain"]["value"].as_f64().unwrap() - 0.005).abs() < 1e-10);
    let probe = s.execute("solve.probe", &json!({"step":"Pull","nodes":[1,2,3,4]})).unwrap();
    let expected = [[0., 0., 0.], [0.0065, 0., 0.], [0., -0.002875, 0.], [0., 0., -0.002875]];
    for (i, node) in probe["nodes"].as_array().unwrap().iter().enumerate() {
        for axis in 0..3 {
            assert!((node["displacement"][axis].as_f64().unwrap() - expected[i][axis]).abs() < 1e-10);
            let reaction = if i == 0 && axis == 0 { -0.25 } else { 0. };
            assert!((node["reaction"][axis].as_f64().unwrap() - reaction).abs() < 1e-10);
        }
        for axis in 0..6 {
            let stress = if axis == 0 { 1.5 } else { 0. };
            assert!((node["stress"][axis].as_f64().unwrap() - stress).abs() < 1e-10);
        }
    }
    probe
}
#[test]
fn small_strain_plastic_table_solves_and_roundtrips_without_losing_history_scope() {
    let dir = std::env::temp_dir().join(format!("fieldcraft-import-plastic-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("input.inp");
    std::fs::write(&input, INPUT).unwrap();
    let mut s = Session::default();
    let read = s.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(read["analysis_complete"], true, "{read}");
    let probe = solve(&mut s);
    let file = dir.join("saved.fieldcraft");
    s.execute("file.save_as", &json!({"path":file,"results":true})).unwrap();
    s.execute("file.new", &json!({})).unwrap();
    s.execute("file.open", &json!({"path":file})).unwrap();
    assert_eq!(s.execute("solve.probe", &json!({"step":"Pull","nodes":[1,2,3,4]})).unwrap(), probe);
    let replay = dir.join("replay.inp");
    s.execute("deck.source.write", &json!({"path":replay})).unwrap();
    assert_eq!(std::fs::read_to_string(&replay).unwrap(), INPUT);
    let normalized = dir.join("normalized.inp");
    s.execute("deck.write", &json!({"format":"abaqus","path":normalized,"units":"m-kg-s-N-Pa"})).unwrap();
    let read = s.execute("deck.read", &json!({"format":"abaqus","path":normalized,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(read["analysis_complete"], true, "{read}");
    assert_eq!(solve(&mut s), probe);
    for text in [
        INPUT.replace("NLGEOM=NO", "NLGEOM=YES"),
        INPUT.replace("*STATIC\n.1,1,.00001,.1", "*FREQUENCY\n1"),
        format!("{INPUT}*STEP,NAME=Next\n*STATIC\n*END STEP\n"),
    ] {
        std::fs::write(&input, &text).unwrap();
        let read = s.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap();
        assert_eq!(read["analysis_complete"], false, "{read}");
        assert!(s.execute("solve.nonlinear", &json!({"step":"Pull"})).is_err());
        s.execute("deck.source.write", &json!({"path":replay})).unwrap();
        assert_eq!(std::fs::read_to_string(&replay).unwrap(), text);
    }
    std::fs::remove_dir_all(dir).unwrap();
}
