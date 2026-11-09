//! Original orthotropic flat-shell membrane patch with a rotated material frame.
use fieldcraft_engine::Session;
use serde_json::{Value, json};
const INPUT: &str = r#"** Original orthotropic membrane patch; SI; material rotation 90 degrees
*NODE
1,0,0,0
2,1,0,0
3,1,1,0
4,0,1,0
*ELEMENT,TYPE=S4,ELSET=Body
1,1,2,3,4
*NSET,NSET=All
1,2,3,4
*NSET,NSET=Left
1,4
*NSET,NSET=Right
2,3
*MATERIAL,NAME=M
*ELASTIC,TYPE=ENGINEERING CONSTANTS
2000,1000,1500,.2,.1,.15,400,500
600
*ORIENTATION,NAME=Axes
1,0,0,0,1,0
3,90
*SHELL SECTION,ELSET=Body,MATERIAL=M,ORIENTATION=Axes
.1,5
*STEP,NAME=Pull
*STATIC
*BOUNDARY
All,3,6
Left,1,1
1,2,2
*CLOAD
Right,1,.5
*NODE PRINT,NSET=All
U,RF
*EL PRINT,ELSET=Body
S
*END STEP
"#;
fn solve(s: &mut Session) -> Value {
    s.execute("solve.static", &json!({"step":"Pull"})).unwrap();
    let p = s.execute("solve.probe", &json!({"step":"Pull","nodes":[1,2,3,4]})).unwrap();
    let u = [[0., 0., 0.], [0.01, 0., 0.], [0.01, -0.001, 0.], [0., -0.001, 0.]];
    for (i, node) in p["nodes"].as_array().unwrap().iter().enumerate() {
        for j in 0..3 {
            assert!((node["displacement"][j].as_f64().unwrap() - u[i][j]).abs() < 1e-10);
            let rf = if (i == 0 || i == 3) && j == 0 { -0.5 } else { 0. };
            assert!((node["reaction"][j].as_f64().unwrap() - rf).abs() < 1e-10);
        }
        for j in 0..6 {
            assert!((node["stress"][j].as_f64().unwrap() - if j == 0 { 10. } else { 0. }).abs() < 1e-10);
        }
    }
    p
}
#[test]
fn rotated_shell_frame_survives_solution_save_and_canonical_round_trip() {
    let dir = std::env::temp_dir().join(format!("fieldcraft-shell-rotation-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("input.inp");
    std::fs::write(&input, INPUT).unwrap();
    let mut s = Session::default();
    let r = s.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(r["analysis_complete"], true, "{r}");
    let probe = solve(&mut s);
    let file = dir.join("saved.fieldcraft");
    s.execute("file.save_as", &json!({"path":file,"results":true})).unwrap();
    s.execute("file.new", &json!({})).unwrap();
    s.execute("file.open", &json!({"path":file})).unwrap();
    assert_eq!(s.execute("solve.probe", &json!({"step":"Pull","nodes":[1,2,3,4]})).unwrap(), probe);
    let replay = dir.join("replay.inp");
    s.execute("deck.source.write", &json!({"path":replay})).unwrap();
    assert_eq!(std::fs::read_to_string(replay).unwrap(), INPUT);
    let path = dir.join("normalized.inp");
    s.execute("deck.write", &json!({"format":"abaqus","path":path,"units":"m-kg-s-N-Pa"})).unwrap();
    let r = s.execute("deck.read", &json!({"format":"abaqus","path":path,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(r["analysis_complete"], true, "{r}");
    solve(&mut s);
    std::fs::remove_dir_all(dir).unwrap();
}
