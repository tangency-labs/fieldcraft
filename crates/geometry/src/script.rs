//! Designs from SolidCraft command scripts (feature `scripts`): the recipes of the test corpus.
//!
//! A script is SolidCraft's own script form, a JSON array of `{"command", "params"}` or an object
//! with a `commands` array, the same file `solidcraft-cli run` takes. [`script_design`] runs it in
//! a SolidCraft session and saves the result as `.solidcraft` file bytes, which then open through
//! [`crate::open_design`] like any design a person saved in SolidCraft.

use solidcraft_kernel::{KernelError, guard};

use crate::{Error, Result};

/// Largest script read (SolidCraft's design limit is far larger than any recipe needs).
pub const MAX_SCRIPT_BYTES: usize = 4 << 20;

/// Run a SolidCraft command script and return the design it builds as `.solidcraft` file bytes.
/// A command that fails stops the script and is the error.
pub fn script_design(script: &str) -> Result<Vec<u8>> {
    if script.len() > MAX_SCRIPT_BYTES {
        return Err(Error::Invalid(format!("script is larger than {MAX_SCRIPT_BYTES} bytes")));
    }
    let value: serde_json::Value = serde_json::from_str(script).map_err(|e| Error::Invalid(format!("script is not JSON: {e}")))?;
    // A failing command is the recipe's error; only a panic is the kernel's.
    let doc = guard("run script", || {
        let mut session = solidcraft_engine::Session::default();
        Ok(session.run_script(&value).map(|_| session.doc.clone()))
    })?
    .map_err(|e| Error::Invalid(format!("script: {e}")))?;
    guard("save design", || solidcraft_io::write_design_container(&doc).map_err(|e| KernelError::Internal(format!("save design: {e}"))))
        .map_err(Error::from)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::open_design;

    #[test]
    fn a_script_builds_a_design_that_opens() {
        let script = r#"{"commands": [
            {"command": "solid.box", "params": {"length": 10, "width": 20, "height": 30, "body_name": "B"}}
        ]}"#;
        let design = open_design(&script_design(script).unwrap()).unwrap();
        assert!(design.errors.is_empty(), "{:?}", design.errors);
        let m = design.body("B").unwrap().measure().unwrap();
        assert!((m.volume - 6000.0).abs() < 1e-6, "{}", m.volume);
    }

    /// Bodies of components come out where the assembly places them: through every occurrence
    /// (instances are bodies of their own) and through rotations as well as translations.
    #[test]
    fn component_bodies_are_placed_by_their_occurrences() {
        let script = r#"{"commands": [
            {"command": "component.create", "params": {"name": "Block"}},
            {"command": "solid.box", "params": {"length": 10, "width": 4, "height": 2, "body_name": "B"}},
            {"command": "component.activate", "params": {"component": "root"}},
            {"command": "occurrence.copy", "params": {"component": "Block", "translate": [100, 0, 0]}},
            {"command": "occurrence.ground", "params": {"occurrence": "Block:1", "grounded": false}},
            {"command": "occurrence.move", "params": {"occurrence": "Block:1", "axis": [0, 0, 1], "angle": "90 deg", "translate": [0, 50, 0]}}
        ]}"#;
        let design = open_design(&script_design(script).unwrap()).unwrap();
        assert!(design.errors.is_empty(), "{:?}", design.errors);
        let mut boxes: Vec<([f64; 3], [f64; 3], f64)> = design
            .bodies
            .iter()
            .map(|b| {
                let m = b.measure().unwrap();
                (m.bbox_min, m.bbox_max, m.volume)
            })
            .collect();
        boxes.sort_by(|a, b| a.0[0].total_cmp(&b.0[0]));
        assert_eq!(boxes.len(), 2, "{:?}", design.bodies.iter().map(|b| &b.name).collect::<Vec<_>>());
        let close = |a: [f64; 3], b: [f64; 3]| a.iter().zip(&b).all(|(x, y)| (x - y).abs() < 1e-6);
        // The first occurrence turned 90 degrees about Z, then moved 50 along Y.
        assert!(close(boxes[0].0, [-4.0, 50.0, 0.0]) && close(boxes[0].1, [0.0, 60.0, 2.0]), "{boxes:?}");
        // The instance, moved 100 along X.
        assert!(close(boxes[1].0, [100.0, 0.0, 0.0]) && close(boxes[1].1, [110.0, 4.0, 2.0]), "{boxes:?}");
        assert!(boxes.iter().all(|b| (b.2 - 80.0).abs() < 1e-6));
        // Placed bodies keep their persistent face names.
        assert!(design.bodies.iter().all(|b| b.faces().len() == 6), "{:?}", design.bodies.iter().map(|b| b.faces().len()).collect::<Vec<_>>());
    }

    #[test]
    fn hostile_scripts_are_an_error() {
        for s in ["", "{", "[1]", "{}", r#"[{"command": "no.such"}]"#, r#"[{"command": "solid.box", "params": {"length": -1e400}}]"#] {
            assert!(script_design(s).is_err(), "{s}");
        }
        assert!(script_design(&" ".repeat(MAX_SCRIPT_BYTES + 1)).is_err());
    }
}
