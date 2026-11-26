//! The marker `deck.openfoam` writes beside a case (`fieldcraft.json`), frozen. The marker is the
//! contract between the writer and `results.openfoam`: a case written by an earlier FieldCraft
//! must still read back. `tests/openfoam/marker-1129.json` is the text the writer produces for the
//! sample air case in the #1129 format (no gas key), frozen so that a change to the writer's
//! marker fails here, and so that the frozen text still passes `Marker::check`.

use serde_json::json;

use super::MARKER;
use super::Marker;
use super::tests_gas_write::air;

/// The #1129 marker of an incompressible case, as the writer produces it for the sample air case.
const FROZEN: &str = include_str!("../../../tests/openfoam/marker-1129.json");

#[test]
fn the_writer_still_produces_the_frozen_1129_marker() {
    let mut a = air(false, false);
    let case = a.dir.join("marker-1129");
    a.s.execute("deck.openfoam", &json!({"path": case.to_str().unwrap()})).unwrap();
    let written = std::fs::read_to_string(case.join(MARKER)).unwrap();
    assert_eq!(written, FROZEN, "the writer's marker moved from the frozen #1129 format");
}

#[test]
fn the_frozen_1129_marker_still_reads() {
    let read: Marker = serde_json::from_str(FROZEN).unwrap();
    read.check().unwrap();
    assert_eq!((read.step.as_str(), read.solver.as_str(), read.time_precision), ("Flow", "pimpleFoam", 9));
    assert!(read.gas.is_none(), "a #1129 marker records no gas");
}
