//! Explicit, test-only transfer of hash-identified native results between frozen IO builds.
//! Not a fallback reader: ordinary file.open must reject overlapping old result kinds.
use fieldcraft_io::{Extras, ModelFile, codec, results};
use fieldcraft_solver::{StepResult, transient::Response};
use serde_json::{Value, json};
use std::path::Path;

fn write(file: &ModelFile, include_results: bool) -> Vec<u8> {
    let extras = Extras {
        deck_source: file.deck_source.as_ref(),
        results: if include_results {
            file.results.iter().map(|r| (r.result.step(), results::StoredRef { definition: r.definition, result: r.result.as_ref() })).collect()
        } else {
            Vec::new()
        },
        embed: file.embedded.iter().map(|(k, v)| (k.clone(), v.as_slice())).collect(),
        run: file.run.clone(),
    };
    fieldcraft_io::write_with(&file.model, &file.geometry, &file.journal, &file.journal_header, &file.view, &extras).unwrap()
}
#[test]
#[ignore = "explicit source/SHA-approved local evidence transfer; no solving"]
fn verified_source_transfer() {
    let config = std::fs::read(std::env::var("FIELDCRAFT_RESULT_TRANSFER").unwrap()).unwrap();
    let config: Value = serde_json::from_slice(&config).unwrap();
    let root = Path::new(config["root"].as_str().unwrap());
    let out = root.join(config["output"].as_str().unwrap());
    let exporting = config["mode"] == "export";
    let mut report = Vec::new();
    for case in config["cases"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let original = root.join(case["path"].as_str().unwrap());
        let original_bytes = std::fs::read(&original).unwrap();
        assert_eq!(fieldcraft_io::sha256_hex(&original_bytes), case["sha256"].as_str().unwrap());
        let folder = out.join(name);
        if exporting {
            assert_eq!(results::MAGIC, b"MCRSLT05", "export must use the verified pre-repair IO build");
            std::fs::create_dir_all(&folder).unwrap();
            let file = fieldcraft_io::read(&original_bytes).unwrap();
            assert_eq!(file.results.len(), 1);
            let stored = &file.results[0];
            let results::Solved::Implicit(StepResult::Dynamics(response)) = &stored.result else { panic!("expected identified dynamics sample") };
            let payload = codec::encode(response).unwrap();
            let model = write(&file, false);
            std::fs::write(folder.join("model-only.fieldcraft"), &model).unwrap();
            std::fs::write(folder.join("response.payload"), &payload).unwrap();
            std::fs::write(folder.join("response.json"), serde_json::to_vec_pretty(response).unwrap()).unwrap();
            let packet = json!({"layout":"Response_cf1c04cc", "source":case,"definition":stored.definition,
                "payload_sha256":fieldcraft_io::sha256_hex(&payload),"model_sha256":fieldcraft_io::sha256_hex(&model),
                "reader_source":config["reader_source"],"reader_io_sha256":config["reader_io_sha256"]});
            std::fs::write(folder.join("packet.json"), serde_json::to_vec_pretty(&packet).unwrap()).unwrap();
            report.push(json!({"name":name,"exported":true,"frames":response.frames.len()}));
        } else {
            assert_eq!(results::MAGIC, b"MCRSLT06");
            let error = fieldcraft_io::read(&original_bytes).unwrap_err().to_string();
            assert!(error.contains("conflicting producer layouts"), "{error}");
            let packet: Value = serde_json::from_slice(&std::fs::read(folder.join("packet.json")).unwrap()).unwrap();
            assert_eq!(packet["source"], *case);
            assert_eq!(packet["layout"], "Response_cf1c04cc");
            assert_eq!(packet["reader_source"], config["reader_source"]);
            assert_eq!(packet["reader_io_sha256"], config["reader_io_sha256"]);
            let model = std::fs::read(folder.join("model-only.fieldcraft")).unwrap();
            let payload = std::fs::read(folder.join("response.payload")).unwrap();
            assert_eq!(fieldcraft_io::sha256_hex(&model), packet["model_sha256"]);
            assert_eq!(fieldcraft_io::sha256_hex(&payload), packet["payload_sha256"]);
            let response: Response = codec::decode(&payload).unwrap();
            let definition: [u8; 32] = serde_json::from_value(packet["definition"].clone()).unwrap();
            let mut file = fieldcraft_io::read(&model).unwrap();
            assert_eq!(results::definition_hash(&file.model, &response.step, &file.geometry), Some(definition));
            file.results.push(results::Stored { definition, result: results::Solved::Implicit(StepResult::Dynamics(response)) });
            let bytes = write(&file, true);
            let restored = fieldcraft_io::read(&bytes).unwrap();
            assert_eq!(restored.results, file.results, "all typed result values retained exactly");
            assert_eq!(restored.model, file.model);
            std::fs::write(folder.join("restored.fieldcraft"), &bytes).unwrap();
            report.push(json!({"name":name,"original_refusal":error,"typed_values_exact":true,"native_sha256":fieldcraft_io::sha256_hex(&bytes)}));
        }
        assert_eq!(std::fs::read(&original).unwrap(), original_bytes, "original unchanged");
    }
    std::fs::write(out.join(if exporting { "export-report.json" } else { "restore-report.json" }), serde_json::to_vec_pretty(&report).unwrap())
        .unwrap();
}

#[test]
#[ignore = "write version 06 fixtures from immutable version 05 kinds 0-5 without solving"]
fn rewrite_unambiguous_fixtures() {
    assert_eq!(results::MAGIC, b"MCRSLT06");
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../engine/tests/results");
    let out = root.join("MCRSLT06");
    std::fs::create_dir_all(&out).unwrap();
    for name in ["cantilever", "couette", "impact", "blocks"] {
        let file = fieldcraft_io::read(&std::fs::read(root.join("MCRSLT05").join(format!("{name}.fieldcraft"))).unwrap()).unwrap();
        let bytes = write(&file, true);
        assert_eq!(fieldcraft_io::read(&bytes).unwrap().results, file.results);
        std::fs::write(out.join(format!("{name}.fieldcraft")), bytes).unwrap();
    }
}
