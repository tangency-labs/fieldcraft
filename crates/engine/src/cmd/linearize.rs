//! `results.linearize`: membrane, bending and peak stresses along a stress classification line
//! through a wall, by the stress-integration method of ASME VIII-2 Annex 5-A
//! ([`crate::linearise`]). The line is sampled with the solid shape functions of
//! [`crate::path`]; the stresses are each element's own nodal stresses (extrapolated from its
//! integration points), not the average over neighbouring elements, which would smear the
//! distribution at a material or thickness change.

use serde_json::{Value, json};

use super::path::{current_result, path_index};
use super::results::{Want, step_for};
use super::xyz;
use crate::linearise::{COMPONENTS, MAX_POINTS, MIN_POINTS, frame, linearise};
use crate::requirement::Requirement;
use crate::results::{Kind, shown_all};
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};
use fieldcraft_model::units::Quantity;

pub const LINEARIZE: CommandSpec = CommandSpec {
    id: "results.linearize",
    label: "Linearize",
    params: r#"{"step": "Pressure", "from": [x, y, z], "to": [x, y, z], "points": 41} — linearises the stresses of a solved linear static step along the line from `from` to `to` (a stress classification line through a wall, in undeformed reference coordinates, both ends on its surfaces), following ASME VIII-2 Annex 5-A: each containing solid element's own nodal stresses (extrapolated from its integration points, not averaged with its neighbours) are interpolated with its shape functions at `points` equally spaced stations (odd, 5 to 1001, default 41) and integrated with Simpson's rule. Replies, in the line's local frame (`frame`: x along the line, z the global axis least aligned with it, Z preferred, y completing the right-handed frame; tensors as `components`, xx, yy, zz, xy, xz, yz), the `membrane` tensor, the `bending` tensor at `from` and at `to` (the two normal stresses in the section and the shear between them carry bending; the normal stress along the line and the shears on its planes do not), the `peak` and `total` tensors at both ends, the `equivalent` (von Mises) stresses of the membrane, of membrane plus bending and of the total at both ends, the `length`, and the `samples` (`x` from `from`, `stress`). Refused by name: a stale result; a frequency, nonlinear static, explicit, heat or fluid step; axisymmetric and other 2D models (not supported yet); coincident ends or a length under 1e-9 of the model's size; an end or station outside the solid elements or across a gap (the first such point is named); an even number of points. Allowable stresses are not judged: they depend on the design code edition and the material. A query: not journalled"#,
    form: Some(super::forms::RESULTS_LINEARIZE),
    kind: CommandKind::Session,
    requires: &[Requirement::Result],
    run: linearize,
};

/// The shortest line linearised, as a share of the model's size (the diagonal of its solid
/// elements' bounds): anything shorter is taken for coincident ends.
const MIN_LENGTH: f64 = 1e-9; // unit-factors: a relative tolerance of the path length, not a conversion

fn linearize(s: &mut Session, p: &Value) -> Result<Value> {
    fieldcraft_model::progress::check()?;
    let step = step_for(s, p, Want::Any)?;
    let result = current_result(s, &step)?;
    let Kind::Static(r) = result.kind() else {
        return Err(EngineError::invalid(
            "step",
            format!("step `{step}` is not a linear static step: stresses are linearised for linear static steps only"),
        ));
    };
    if r.nonlinear.is_some() {
        return Err(EngineError::invalid("step", "not supported yet: linearisation of nonlinear static steps"));
    }
    if let Some(extent) = s.model.sections().values().find_map(|section| section.kind.extent()) {
        let what = match extent {
            fieldcraft_model::Extent::Revolved => "axisymmetric",
            fieldcraft_model::Extent::Thickness(_) => "2D plane stress and plane strain",
        };
        return Err(EngineError::invalid("step", format!("not supported yet: linearisation of {what} models")));
    }
    if r.elements.is_empty() {
        return Err(EngineError::invalid("step", format!("step `{step}` kept no element stresses to linearise")));
    }
    if !r.elements.is_sorted_by_key(|e| e.label) {
        return Err(EngineError::invalid("step", format!("step `{step}` keeps its element stresses out of label order; solve it again")));
    }
    let from = xyz(p, "from")?;
    let to = xyz(p, "to")?;
    let points = super::mesh::optional_count(p, "points", 41, MIN_POINTS as u64..=MAX_POINTS as u64)?;
    if points.is_multiple_of(2) {
        return Err(EngineError::invalid("points", "must be odd: Simpson's rule integrates pairs of intervals"));
    }
    let index = path_index(s, &step, &result)?;
    let size = index.bounds().map_or(0.0, |b| fieldcraft_fem::frame::norm(fieldcraft_fem::frame::sub(b.max, b.min)));
    let (axes, length) = frame(from, to).filter(|(_, length)| *length >= MIN_LENGTH * size).ok_or_else(|| {
        EngineError::invalid(
            "to",
            "`from` and `to` coincide or are not finite (a line shorter than 1e-9 of the model's size has no thickness to linearise over)",
        )
    })?;
    let (distance, stations) = crate::path::stations(&[from, to], points)?;
    let sampled = index.tensors(&s.model, &stations, |label| {
        r.elements.binary_search_by_key(&label, |e| e.label).ok().and_then(|k| r.elements.get(k)).map(|e| e.nodal.as_slice())
    })?;
    let mut samples = Vec::with_capacity(points);
    let mut elements = Vec::with_capacity(points);
    for (i, (point, sample)) in stations.iter().zip(&sampled).enumerate() {
        let Some((label, stress)) = sample else {
            let shown = shown_all(Quantity::Length, *point, s.units());
            let (key, which) = match i {
                0 => ("from", "`from`".to_string()),
                _ if i + 1 == points => ("to", "`to`".to_string()),
                _ => ("to", format!("station {i}")),
            };
            return Err(EngineError::invalid(
                key,
                format!(
                    "{which} at [{}, {}, {}] {} is outside the solid elements of the mesh (or in a gap between them); keep the line inside one wall",
                    shown[0],
                    shown[1],
                    shown[2],
                    s.units().symbol(Quantity::Length)
                ),
            ));
        };
        samples.push(*stress);
        elements.push(*label);
    }
    let out = linearise(&samples, length, axes)?;
    Ok(json!({
        "step": step,
        "from": from,
        "to": to,
        "length": length,
        "points": points,
        "frame": {"x": axes[0], "y": axes[1], "z": axes[2]},
        "components": COMPONENTS,
        "membrane": out.membrane,
        "bending": {"from": out.bending, "to": out.bending.map(|v| -v)},
        "peak": {"from": out.peak[0], "to": out.peak[1]},
        "total": {"from": out.total[0], "to": out.total[1]},
        "equivalent": {
            "membrane": out.membrane_equivalent,
            "membrane_plus_bending": {"from": out.membrane_plus_bending_equivalent[0], "to": out.membrane_plus_bending_equivalent[1]},
            "total": {"from": out.total_equivalent[0], "to": out.total_equivalent[1]},
        },
        "samples": distance.iter().zip(&samples).map(|(x, stress)| json!({"x": x, "stress": stress})).collect::<Vec<_>>(),
        "elements": elements,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_fem::rotation::{Mat3, apply, rotate, transpose};

    const RECIPE: &str = include_str!("../../../../tests/oracle/lame-cylinder/recipe.json");

    fn about_z(degrees: f64) -> Mat3 {
        let a = degrees.to_radians();
        [[a.cos(), -a.sin(), 0.0], [a.sin(), a.cos(), 0.0], [0.0, 0.0, 1.0]]
    }

    /// Lamé's thick-walled cylinder: a = 10, b = 20, p = 100, ν = 0.3, plane strain, from the
    /// oracle case, rotated by `angle` about Z (the symmetry planes follow the rotation).
    fn lame(angle_degrees: f64) -> Session {
        lame_divisions(angle_degrees, [8, 16])
    }

    fn lame_divisions(angle_degrees: f64, divisions: [u64; 2]) -> Session {
        let q = about_z(angle_degrees);
        let mut recipe: Value = serde_json::from_str(RECIPE).unwrap();
        for entry in recipe.as_array_mut().unwrap() {
            let params = &mut entry["params"];
            if let Some(corners) = params.get_mut("corners") {
                for corner in corners.as_array_mut().unwrap() {
                    let p: [f64; 3] = serde_json::from_value(corner.clone()).unwrap();
                    *corner = json!(apply(&q, &p));
                }
            }
            if let Some(d) = params.get_mut("divisions") {
                *d = json!([divisions[0], divisions[1], 1]);
            }
            if let Some(plane) = params.get_mut("plane") {
                let normal = match plane.as_str().unwrap() {
                    "x" => [1.0, 0.0, 0.0],
                    "y" => [0.0, 1.0, 0.0],
                    _ => [0.0, 0.0, 1.0],
                };
                *plane = json!({"normal": apply(&q, &normal)});
            }
        }
        let mut s = Session::default();
        s.run_script(&recipe).unwrap();
        s
    }

    fn line(angle_degrees: f64) -> Value {
        let q = about_z(angle_degrees);
        json!({"step": "Pressure", "from": apply(&q, &[10.0, 0.0, 0.5]), "to": apply(&q, &[20.0, 0.0, 0.5])})
    }

    /// The cylinder rotated by `angle`, linearised along its rotated radial line, checked
    /// against the closed forms.
    fn linearised(angle_degrees: f64) -> (Session, Value) {
        let mut s = lame(angle_degrees);
        let reply = s.execute("results.linearize", &line(angle_degrees)).unwrap();
        check_lame(&reply);
        (s, reply)
    }

    /// The closed forms of the issue: A₀ = p a² / (b² − a²), B₀ = A₀ b², σ_θ = A₀ + B₀/r²,
    /// σ_r = A₀ − B₀/r², σ_z = 2νA₀, linearised over the wall t = b − a.
    struct Reference {
        hoop_membrane: f64,
        hoop_bending_a: f64,
        hoop_peak: [f64; 2],
        radial_membrane: f64,
        axial_membrane: f64,
        membrane_von_mises: f64,
    }

    fn reference() -> Reference {
        let (a, b, p, nu) = (10.0_f64, 20.0_f64, 100.0_f64, 0.3_f64);
        let t = b - a;
        let a0 = p * a * a / (b * b - a * a);
        let b0 = a0 * b * b;
        let c = (a + b) / 2.0;
        let hoop = |r: f64| a0 + b0 / (r * r);
        let hoop_membrane = p * a / (b - a);
        let hoop_bending_a = 6.0 * b0 / (t * t) * (c * (1.0 / a - 1.0 / b) - (b / a).ln());
        let radial_membrane = -p * a / (a + b);
        let axial_membrane = 2.0 * nu * a0;
        Reference {
            hoop_membrane,
            hoop_bending_a,
            hoop_peak: [hoop(a) - hoop_membrane - hoop_bending_a, hoop(b) - hoop_membrane + hoop_bending_a],
            radial_membrane,
            axial_membrane,
            membrane_von_mises: fieldcraft_fem::von_mises(&[radial_membrane, hoop_membrane, axial_membrane, 0.0, 0.0, 0.0]),
        }
    }

    fn within(reply: &Value, pointer: &str, expected: f64, relative: f64) {
        let got = reply.pointer(pointer).and_then(Value::as_f64).unwrap_or_else(|| panic!("{pointer} in {reply}"));
        assert!((got - expected).abs() <= relative * expected.abs(), "{pointer}: {got} vs {expected} ({relative:e})");
    }

    fn check_lame(reply: &Value) {
        let r = reference();
        assert_eq!(reply["components"], json!(["xx", "yy", "zz", "xy", "xz", "yz"]));
        within(reply, "/length", 10.0, 1e-12);
        within(reply, "/membrane/1", r.hoop_membrane, 0.01);
        within(reply, "/bending/from/1", r.hoop_bending_a, 0.01);
        within(reply, "/bending/to/1", -r.hoop_bending_a, 0.01);
        // The peak is the small difference of the surface stress and the membrane and bending
        // parts, each 0.2 to 0.4% off on this mesh (0.3, 0.23 and 0.18 in stress units): added
        // up, 3.4% of the 21.2 peak at the bore. `finer_mesh_brings_the_peak_within_one_percent`
        // shows it converging.
        within(reply, "/peak/from/1", r.hoop_peak[0], 0.05);
        within(reply, "/peak/to/1", r.hoop_peak[1], 0.05);
        within(reply, "/membrane/0", r.radial_membrane, 0.01);
        within(reply, "/membrane/2", r.axial_membrane, 0.01);
        within(reply, "/equivalent/membrane", r.membrane_von_mises, 0.01);
        // The components without a bending part report none, and the shears vanish.
        assert_eq!(reply["bending"]["from"][0], 0.0);
        assert_eq!(reply["bending"]["from"][3], 0.0);
        for c in 3..6 {
            assert!(reply["membrane"][c].as_f64().unwrap().abs() < 0.5, "{reply}");
        }
    }

    #[test]
    fn lame_cylinder_matches_the_closed_forms_converges_and_is_not_an_edit() {
        let mut s = lame(0.0);
        let journal = s.journal.len();
        let revision = s.revision();
        let reply = s.execute("results.linearize", &line(0.0)).unwrap();
        check_lame(&reply);
        assert_eq!(reply["points"], 41);
        assert_eq!(reply["samples"].as_array().unwrap().len(), 41);
        assert_eq!(reply["samples"][40]["x"], 10.0);
        assert_eq!(reply["frame"], json!({"x": [1.0, 0.0, 0.0], "y": [0.0, 1.0, 0.0], "z": [0.0, 0.0, 1.0]}));
        assert_eq!(reply["units"]["pressure"], "MPa");
        assert_eq!(reply["units"]["length"], "mm");
        assert_eq!(s.journal.len(), journal);
        assert_eq!(s.revision(), revision);
        let mut fine = line(0.0);
        fine["points"] = json!(81);
        let finer = s.execute("results.linearize", &fine).unwrap();
        check_lame(&finer);
        for pointer in ["/membrane/0", "/membrane/1", "/membrane/2", "/bending/from/1"] {
            within(&finer, pointer, reply.pointer(pointer).and_then(Value::as_f64).unwrap(), 1e-3);
        }
    }

    /// Quadratic elements converge the stresses at second order: twice the divisions in both
    /// directions bring the surface stress and the integrals within about a quarter of their
    /// error, and the peak at the bore from 3.4% to under 1%.
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn finer_mesh_brings_the_peak_within_one_percent() {
        let mut s = lame_divisions(0.0, [16, 32]);
        let reply = s.execute("results.linearize", &line(0.0)).unwrap();
        eprintln!("16 × 32 hex20: membrane {} bending {} peak {} total {}", reply["membrane"], reply["bending"], reply["peak"], reply["total"]);
        let r = reference();
        within(&reply, "/membrane/1", r.hoop_membrane, 0.001);
        within(&reply, "/bending/from/1", r.hoop_bending_a, 0.0025);
        within(&reply, "/peak/from/1", r.hoop_peak[0], 0.01);
        within(&reply, "/peak/to/1", r.hoop_peak[1], 0.01);
    }

    #[test]
    fn a_cylinder_rotated_about_z_gives_the_same_local_frame_values() {
        let (_, reply) = linearised(30.0);
        let frame: [[f64; 3]; 3] = ["x", "y", "z"].map(|k| serde_json::from_value(reply["frame"][k].clone()).unwrap());
        let q = about_z(30.0);
        let expected = transpose(&q);
        for i in 0..3 {
            for k in 0..3 {
                assert!((frame[i][k] - expected[i][k]).abs() < 1e-12, "{frame:?}");
            }
        }
        // The same frame, transposed, carries the local membrane back to the global tensor the
        // unrotated cylinder reports.
        let (_, unrotated) = linearised(0.0);
        let local: [f64; 6] = serde_json::from_value(reply["membrane"].clone()).unwrap();
        let global = rotate(&transpose(&frame), &local);
        let reference: [f64; 6] = serde_json::from_value(unrotated["membrane"].clone()).unwrap();
        let rotated = rotate(&q, &reference);
        for k in 0..6 {
            assert!((global[k] - rotated[k]).abs() < 1e-2, "{global:?} vs {rotated:?}");
        }
    }

    #[test]
    fn refusals_are_named_and_change_nothing() {
        let (mut s, reply) = linearised(0.0);
        let journal = s.journal.len();
        let revision = s.revision();
        let refuse = |s: &mut Session, change: Value, expect: &str| {
            let mut p = line(0.0);
            for (k, v) in change.as_object().unwrap() {
                p[k] = v.clone();
            }
            let err = s.execute("results.linearize", &p).unwrap_err().to_string();
            assert!(err.contains(expect), "{p}: {err}");
            assert_eq!(s.journal.len(), journal);
            assert_eq!(s.revision(), revision);
        };
        refuse(&mut s, json!({"from": [10, 0, 0.5], "to": [10, 0, 0.5]}), "coincide");
        refuse(&mut s, json!({"from": [10, 0, 0.5], "to": [10.000000000001, 0, 0.5]}), "coincide");
        refuse(&mut s, json!({"points": 40}), "odd");
        refuse(&mut s, json!({"points": 3}), "5..=1001");
        refuse(&mut s, json!({"points": 1003}), "5..=1001");
        refuse(&mut s, json!({"points": 1.5}), "5..=1001");
        refuse(&mut s, json!({"points": -1}), "5..=1001");
        refuse(&mut s, json!({"from": [9, 0, 0.5]}), "`from` at [9, 0, 0.5] mm is outside");
        refuse(&mut s, json!({"to": [21, 0, 0.5]}), "station 37 at [20.175, 0, 0.5] mm is outside");
        refuse(&mut s, json!({"from": [19, 0, 0.5], "to": [20.02, 0, 0.5]}), "`to` at [20.02, 0, 0.5] mm is outside");
        refuse(&mut s, json!({"from": [15, 0, 0.5], "to": [15, 0, 3]}), "station 9 at [15, 0, 1.0625] mm is outside");
        refuse(&mut s, json!({"step": "Nope"}), "no result");
        // A gap: remove an element along the line, keep the result's node values.
        let mut gapped = lame(0.0);
        gapped.execute("mesh.element.remove", &json!({"label": reply["elements"][20]})).unwrap();
        let result = gapped.latest_result("Pressure").unwrap().0.as_ref().clone();
        gapped.store_result(result);
        let err = gapped.execute("results.linearize", &line(0.0)).unwrap_err().to_string();
        assert!(err.contains("station") && err.contains("gap"), "{err}");
        // A stale result.
        s.execute("model.material.edit", &json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 100000.0, "poisson": 0.3}})).unwrap();
        let err = s.execute("results.linearize", &line(0.0)).unwrap_err().to_string();
        assert!(err.contains("stale"), "{err}");
    }

    #[test]
    fn frequency_nonlinear_and_axisymmetric_steps_are_refused_by_name() {
        let mut s = Session::default();
        crate::cmd::solve::tests::cantilever(&mut s);
        s.execute("solve.frequency", &json!({"step": "Modes"})).unwrap();
        let err = s.execute("results.linearize", &json!({"step": "Modes", "from": [1, 0.5, 0.5], "to": [9, 0.5, 0.5]})).unwrap_err().to_string();
        assert!(err.contains("not a linear static step"), "{err}");
        s.execute("model.step.create", &json!({"name": "NL", "procedure": {"type": "nonlinear_static", "initial_increment": 0.5}, "boundary_conditions": ["clamp"], "loads": ["push"]}))
            .unwrap();
        s.execute("solve.nonlinear", &json!({"step": "NL"})).unwrap();
        let err = s.execute("results.linearize", &json!({"step": "NL", "from": [1, 0.5, 0.5], "to": [9, 0.5, 0.5]})).unwrap_err().to_string();
        assert!(err.contains("not supported yet: linearisation of nonlinear static steps"), "{err}");

        for (recipe, what) in [
            (include_str!("../../../../tests/oracle/lame-cylinder-axisymmetric/recipe.json"), "axisymmetric"),
            (include_str!("../../../../tests/oracle/lame-cylinder-plane-strain/recipe.json"), "2D plane stress and plane strain"),
        ] {
            let mut planar = Session::default();
            planar.run_script(&serde_json::from_str(recipe).unwrap()).unwrap();
            let step = planar.result_steps().next().unwrap().0.to_string();
            let err = planar.execute("results.linearize", &json!({"step": step, "from": [10, 0, 0], "to": [20, 0, 0]})).unwrap_err().to_string();
            assert!(err.contains(&format!("not supported yet: linearisation of {what} models")), "{what}: {err}");
        }
    }
}
