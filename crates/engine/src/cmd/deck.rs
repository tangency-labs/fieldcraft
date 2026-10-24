//! Solver decks: write model data, or read a flat Abaqus deck into a new session.

use std::path::Path;

use fieldcraft_io::roots::Access;
use serde_json::{Value, json};

use super::{opt, req};
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const WRITE: CommandSpec = CommandSpec {
    id: "deck.write",
    label: "Write Solver Deck",
    params: r#"{"format": "abaqus|nastran|lsdyna", "path": "part.inp", "steps": ["Load"], "units": "mm-t-s-N-MPa|m-kg-s-N-Pa|mm-kg-ms-kN-GPa|mm-g-ms-N-MPa|in-lbf-s-psi"} — writes the model as a solver deck in one consistent unit system, stated in its header (default: the display system's, ignoring its overrides) (formats: fieldcraft_decks::FORMATS); `steps` (optional) writes only those steps, as a Nastran run needs (static or frequency steps, not both); LS-DYNA writes several steps as *CASE analyses (run it with CASE); fails naming every error of the model check and everything the format cannot express"#,
    form: Some(super::forms::DECK),
    kind: CommandKind::Session,
    requires: &[Requirement::Analysable, Requirement::ElementMaterial, Requirement::ElementSection],
    run: write,
};

pub const READ: CommandSpec = CommandSpec {
    id: "deck.read",
    label: "Read Abaqus Deck",
    params: r#"{"format":"abaqus","path":"part.inp","units"?:"mm-t-s-N-MPa|m-kg-s-N-Pa|mm-kg-ms-kN-GPa|mm-g-ms-N-MPa|in-lbf-s-psi"} — replaces the model with a flat Abaqus deck's mesh, sets, element surfaces, elastic/density materials and sections; input numbers use this consistent preset (default: the current display preset, ignoring overrides), then store SI; the reply states units, unread keywords with original lines (first 100 and total count), and warnings. No steps, loads, conditions or unknown text are retained; INCLUDE/parts/assemblies and unsupported physical variants are refused. Empty supported data is refused. Success clears old results/design/path/undo and is unsaved; errors leave the whole session intact."#,
    form: Some(super::forms::DECK_READ),
    kind: CommandKind::Session,
    requires: &[],
    run: read,
};

fn read(s: &mut Session, p: &Value) -> Result<Value> {
    let format: &str = req(p, "format")?;
    if format != "abaqus" {
        return Err(EngineError::invalid("format", "only abaqus is supported by deck.read"));
    }
    let units = match opt::<String>(p, "units")? {
        Some(label) => fieldcraft_model::units::Preset::from_label(&label)
            .ok_or_else(|| EngineError::invalid("units", "choose a known consistent unit preset"))?,
        None => s.units().preset(),
    };
    let path = super::file_param(s, p, "path", Access::Read)?;
    let bytes = super::file::read_file(&path, fieldcraft_decks::abaqus::read::MAX_INPUT_BYTES as u64)?;
    let text = std::str::from_utf8(&bytes).map_err(|_| EngineError::invalid("path", "the input deck is not UTF-8 text"))?;
    let imported = fieldcraft_decks::abaqus::read_in(text, units)?;
    let model = &imported.model;
    if model.nodes().is_empty()
        && model.elements().is_empty()
        && model.sets().values().all(fieldcraft_model::Set::is_empty)
        && model.materials().is_empty()
        && model.sections().is_empty()
    {
        return Err(EngineError::invalid(
            "path",
            format!("the deck contains no supported model data ({} unread keywords); the current model is unchanged", imported.unknown_count),
        ));
    }
    let mut kinds = std::collections::BTreeMap::<&str, usize>::new();
    for element in model.elements().values() {
        *kinds.entry(element.kind.name()).or_default() += 1;
    }
    let unknown: Vec<Value> = imported.unknown.iter().map(|u| json!({"keyword":u.keyword,"line":u.line})).collect();
    let reply = json!({"format":"abaqus","path":path,"units":units.label(),"nodes":model.nodes().len(),"elements":kinds,
        "sets":model.sets().len(),"materials":model.materials().len(),"sections":model.sections().len(),"unknown":unknown,
        "unknown_count":imported.unknown_count,"warnings":imported.warnings,"history_imported":false});
    // The new, unsaved document has no journal, but still needs its own dirty revision.
    let opened = Session { model: imported.model, revision: crate::next_revision(), ..Session::default() };
    s.replace_with(opened);
    Ok(reply)
}

fn write(s: &mut Session, p: &Value) -> Result<Value> {
    let name: &str = req(p, "format")?;
    let format = fieldcraft_decks::format(name).ok_or_else(|| {
        let shown: String = name.chars().take(32).collect();
        let known: Vec<&str> = fieldcraft_decks::FORMATS.iter().map(|f| f.name).collect();
        EngineError::Params(format!("unknown deck format `{shown}`; known: {}", known.join(", ")))
    })?;
    let path = super::file_param(s, p, "path", Access::Write)?;
    let steps: Option<Vec<String>> = super::opt(p, "steps")?;
    let units = match p.get("units").and_then(Value::as_str) {
        None => s.units().preset(),
        Some(label) => fieldcraft_model::units::Preset::from_label(label).ok_or_else(|| {
            let known: Vec<&str> = fieldcraft_model::units::Preset::ALL.iter().map(|u| u.label()).collect();
            EngineError::invalid("units", format!("must be one of {}", known.join(", ")))
        })?,
    };
    let deck = format.write_steps(&s.model, steps.as_deref(), units)?;
    fieldcraft_io::fs::write_output(Path::new(&path), deck.text.as_bytes())?;
    Ok(json!({"format": format.name, "path": path, "bytes": deck.text.len(), "units": units.label(), "warnings": deck.warnings}))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_the_defined_model_and_refuses_what_it_cannot() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-deck-write-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("bar.inp");
        let path_str = path.to_string_lossy().into_owned();

        let mut s = Session::default();
        crate::cmd::define::tests::defined_bar(&mut s);
        let journal = s.journal.len();
        let out = s.execute("deck.write", &json!({"format": "abaqus", "path": path_str})).unwrap();
        let written = std::fs::read_to_string(&path).unwrap();
        assert_eq!(written, fieldcraft_decks::abaqus::write(&s.model, fieldcraft_model::units::Preset::MmTonneS).unwrap());
        assert_eq!(out["bytes"], written.len());
        assert_eq!(out["warnings"].as_array().map(Vec::len), Some(1), "`press` is in no step: {out}");
        assert!(written.contains("*ELEMENT, TYPE=C3D8I"), "{written}");
        assert_eq!(s.journal.len(), journal, "writing a deck does not change the model");
        assert!(std::fs::read_dir(&dir).unwrap().count() == 1, "no partial file is left behind");

        // Unknown formats, missing parameters and models that are not ready fail without writing.
        let other = dir.join("other.inp").to_string_lossy().into_owned();
        assert!(s.execute("deck.write", &json!({"format": "nope", "path": other})).unwrap_err().to_string().contains("known: abaqus"));
        assert!(s.execute("deck.write", &json!({"path": other})).is_err());
        assert!(s.execute("deck.write", &json!({"format": "abaqus"})).is_err());
        let mut empty = Session::default();
        assert!(matches!(empty.execute("deck.write", &json!({"format": "abaqus", "path": other})), Err(EngineError::Unavailable(_))));
        assert!(!dir.join("other.inp").exists());
        assert!(matches!(s.execute("deck.write", &json!({"format": "abaqus", "path": dir.to_string_lossy()})), Err(EngineError::File(_))));

        // Nastran: one run is static or frequency steps, so the bar's two steps need two decks.
        let bdf = dir.join("bar.bdf").to_string_lossy().into_owned();
        let e = s.execute("deck.write", &json!({"format": "nastran", "path": bdf})).unwrap_err().to_string();
        assert!(e.contains("not supported by Nastran") && e.contains("separate decks"), "{e}");
        for (step, sol) in [("Bend", "SOL 101"), ("Modes", "SOL 103")] {
            s.execute("deck.write", &json!({"format": "nastran", "path": bdf, "steps": [step]})).unwrap();
            let deck = std::fs::read_to_string(&bdf).unwrap();
            assert!(deck.contains(sol) && deck.contains("RBE2") && deck.ends_with("ENDDATA\n"), "{deck}");
        }
        assert_eq!(s.journal.len(), journal);
        for steps in [json!([]), json!(["Nope"]), json!("Bend"), json!([1, 2]), json!([["Bend"]])] {
            assert!(s.execute("deck.write", &json!({"format": "nastran", "path": bdf, "steps": steps})).is_err(), "{steps}");
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn writes_lsdyna_keyword_decks() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-deck-lsdyna-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("bar.k");
        let mut s = Session::default();
        crate::cmd::define::tests::defined_bar(&mut s);
        let out = s.execute("deck.write", &json!({"format": "lsdyna", "path": path.to_string_lossy()})).unwrap();
        let written = std::fs::read_to_string(&path).unwrap();
        assert_eq!(written, fieldcraft_decks::lsdyna::write(&s.model, fieldcraft_model::units::Preset::MmTonneS).unwrap().text);
        assert_eq!(out["format"], "lsdyna");
        assert!(written.starts_with("*KEYWORD long=y\n") && written.contains("*CASE_BEGIN_2"), "two steps, two cases");
        // One chosen step: a single analysis, no cases.
        s.execute("deck.write", &json!({"format": "lsdyna", "path": path.to_string_lossy(), "steps": ["Modes"]})).unwrap();
        let one = std::fs::read_to_string(&path).unwrap();
        assert!(one.contains("*CONTROL_IMPLICIT_EIGENVALUE") && !one.contains("*CASE_BEGIN") && !one.contains("*CONTROL_IMPLICIT_SOLUTION"), "{one}");
        // What LS-DYNA cannot take is named and nothing is written.
        s.execute("model.load.create", &json!({"name": "heat", "type": "body_force", "set": "all", "vector": [0, 0, 1]})).unwrap();
        s.execute("model.step.activate", &json!({"step": "Bend", "loads": ["heat"]})).unwrap();
        let other = dir.join("heat.k");
        let e = s.execute("deck.write", &json!({"format": "lsdyna", "path": other.to_string_lossy()})).unwrap_err().to_string();
        assert!(e.contains("not supported by LS-DYNA") && e.contains("load `heat`"), "{e}");
        assert!(!other.exists());
        let _ = std::fs::remove_dir_all(&dir);
    }
}

#[cfg(test)]
#[path = "deck/read_tests.rs"]
mod read_tests;
