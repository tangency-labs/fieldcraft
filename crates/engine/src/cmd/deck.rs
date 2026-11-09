//! Solver decks: interpret Abaqus model data and retain original files for explicit replay.

mod source;
pub use source::WRITE_SOURCE;

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
    params: r#"{"format":"abaqus","path":"part.inp","units"?:"mm-t-s-N-MPa|m-kg-s-N-Pa|mm-kg-ms-kN-GPa|mm-g-ms-N-MPa|in-lbf-s-psi"} — replaces the model with an Abaqus deck's mesh, sets, element surfaces, elastic/density materials and sections; expands parts and instances and bounded relative INCLUDE files; input numbers use this consistent preset (default: the current display preset, ignoring overrides), then store SI; the reply states units, unread keywords with original lines (first 100 and total count), and warnings. Original files are stored with the project for deck.source.write. Supported linear static/frequency/buckling and direct single-frequency-list harmonic steps, structural boundaries, concentrated forces/moments, pressure and gravity are mapped; tabular/smooth-step static amplitudes are evaluated at the step end. Unsupported physical options retain source and block affected analysis mappings; inspect analysis_complete, imported_steps and warnings before solving. Constant isotropic MASS and grounded SPRING1 map to native point masses/springs with generated node targets; original discrete ELSET membership remains in the source. Material/global constant damping is mapped; complex load quadratures, frequency ranges and modal harmonic extraction remain unsupported. Unsupported physical variants are refused. Empty supported data is refused. Success clears old results/design/path/undo and is unsaved; errors leave the whole session intact."#,
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
    let limits = fieldcraft_decks::abaqus::Limits {
        nodes: s.limits.max_import_nodes,
        elements: s.limits.max_import_elements,
        members: s.limits.max_import_members,
    };
    let bundle = source::load(s, &path, text)?;
    let imported = fieldcraft_decks::abaqus::read::read_analysis_limited(&bundle.expanded, units, limits).map_err(|e| bundle.source_error(e))?;
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
    let unknown: Vec<Value> = imported
        .unknown
        .iter()
        .map(|u| {
            let at = bundle.location(u.line);
            let original = at.and_then(|at| {
                bundle.files.get(&at.file).and_then(|doc| doc.blocks().iter().position(|b| b.line == at.line).and_then(|index| doc.block_text(index)))
            });
            let preview = original.map(|text| text.chars().take(4096).collect::<String>());
            let truncated = original.zip(preview.as_deref()).is_some_and(|(full, shown)| full.len() != shown.len());
            json!({"keyword":u.keyword,"line":at.map_or(u.line, |at|at.line),"file":at.map(|at| &at.file),
                "original":preview,"original_truncated":truncated})
        })
        .collect();
    let instances: Vec<Value> = imported
        .instances
        .iter()
        .map(|i| {
            json!({"name":i.name,"part":i.part,
        "nodes":i.nodes,"elements":i.elements,"names":i.names})
        })
        .collect();
    let requirements: Vec<Value> = imported
        .mass_requirements
        .iter()
        .map(|m| {
            json!({
                "element_type":m.element_type,"elset":m.elset,"choice":m.choice,"rule":m.rule,
                "conditions":m.conditions,"native_gap":m.native_gap
            })
        })
        .collect();
    let choices: std::collections::BTreeSet<_> = imported.mass_requirements.iter().filter_map(|m| m.choice).collect();
    let mass_report = json!({"applied":false,"conditions_verified":false,
        "evidence_scope":"unconstrained MATRIX GENERATE,MASS; no cross-procedure inference",
        "mixed_choices_require_per_type_or_section":choices.len()>1,"requirements":requirements,
        "warning":"No imported mass default is applied. Native solver defaults are not a verified source mass mapping; consistent/lumped alone does not identify quadrature or lumping formula."});
    let reply = json!({"format":"abaqus","path":path,"units":units.label(),"nodes":model.nodes().len(),"elements":kinds,
        "sets":model.sets().len(),"materials":model.materials().len(),"sections":model.sections().len(),"unknown":unknown,
        "masses":model.masses().len(),"springs":model.springs().len(),"unknown_count":imported.unknown_count,"warnings":imported.warnings,"history_imported":!model.steps().is_empty(),"analysis_complete":imported.analysis_complete,"steps":model.steps().len(),"imported_steps":model.steps().iter().map(|s|s.name.as_str()).collect::<Vec<_>>(),
        "native_lumped_mass_policies":model.sections().iter().filter(|(_,s)|!s.formulation.lumped_mass.is_default()).map(|(n,s)|json!({"section":n,"set":s.set,"strategy":s.formulation.lumped_mass,"activation":"explicit mass=lumped solve"})).collect::<Vec<_>>(),"mass_semantics":mass_report,"source_retained":true,"source_storage":"project","source_files":bundle.files.keys().collect::<Vec<_>>(),"instances":instances});
    // The new, unsaved document has no journal, but still needs its own dirty revision.
    let archive = fieldcraft_io::DeckSource {
        root: bundle.root.clone(),
        files: bundle.files.iter().map(|(name, doc)| (name.clone(), doc.original().to_string())).collect(),
    };
    let opened =
        Session { model: imported.model, deck_source: Some(std::sync::Arc::new(archive)), revision: crate::next_revision(), ..Session::default() };
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
