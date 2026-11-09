use fieldcraft_engine::Session;
use serde_json::json;
#[test]
fn frequency_mass_kind12_roundtrip_and_distinct_mass() {
    let mut s = Session::default();
    for (name, mass) in [("C", "consistent"), ("L", "lumped")] {
        s.execute("file.sample", &json!({"name":"cantilever"})).unwrap();
        s.execute("model.step.create", &json!({"name":"Pulse","procedure":{"type":"frequency","modes":1},"boundary_conditions":[],"loads":[]}))
            .unwrap();
        let r = s.execute("solve.frequency", &json!({"step":"Pulse","mass":mass})).unwrap();
        let frequency = r["modes"][0]["frequency"].as_f64().unwrap();
        assert!(frequency.is_finite());
        let path = std::env::temp_dir().join(format!("kind12-{name}-{}.fieldcraft", std::process::id()));
        s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
        let mut o = Session::default();
        o.execute("file.open", &json!({"path":path})).unwrap();
        let p = o.execute("solve.probe", &json!({"step":"Pulse","nodes":[1]})).unwrap();
        assert_eq!(p["mass_kind"], mass);
        assert!(p["section_mass_policy"].as_str().unwrap().contains("declared section"));
        std::fs::remove_file(path).unwrap();
    }
}
