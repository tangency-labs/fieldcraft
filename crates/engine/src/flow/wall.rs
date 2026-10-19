//! Face-averaged heat received by obstacle walls, in SI: the faces' areas come from the
//! design's tessellation, in millimetres, converted once.
//! The history is interpolated only between recorded samples: frame zero predates it.

use std::collections::BTreeMap;

use fieldcraft_lbm::FluidResult;
use fieldcraft_render::{DisplayMesh, Vec3};
use serde::Serialize;

use crate::results::Pick;
use crate::{EngineError, Result, Session};

#[derive(Clone, Debug, Serialize)]
pub struct FaceHeat {
    pub body: String,
    pub face: String,
    pub area: f64,
    /// Positive into the wall. None before the first recorded history sample.
    pub power: Option<f64>,
    pub flux: Option<f64>,
    #[serde(skip)]
    pub index: u32,
}

fn interpolate(times: &[f64], values: &[f64], time: f64) -> Option<f64> {
    if !time.is_finite()
        || times.len() != values.len()
        || times.is_empty()
        || times.iter().chain(values).any(|v| !v.is_finite())
        || times.windows(2).any(|p| p[0] >= p[1])
    {
        return None;
    }
    let k = times.partition_point(|t| *t < time);
    if times.get(k) == Some(&time) {
        return values.get(k).copied();
    }
    let lo = k.checked_sub(1)?;
    let (a, b) = (*times.get(lo)?, *times.get(k)?);
    let weight = (time - a) / (b - a);
    let value = (1.0 - weight) * values.get(lo)? + weight * values.get(k)?;
    value.is_finite().then_some(value)
}

/// All recorded obstacle faces, using their complete tessellated geometric area (not the
/// visible triangle area or a lattice-link count). Missing geometry is an error, never zero.
pub fn faces(s: &Session, r: &FluidResult, pick: Pick) -> Result<Vec<FaceHeat>> {
    if s.latest_result(&r.step).is_some_and(|(_, current)| !current) {
        return Err(EngineError::Params("wall heat requires current geometry: solve the edited model again before dividing by face area".into()));
    }
    let time = super::frame(r, pick)?.time;
    let mut areas = BTreeMap::new();
    let mut output = Vec::new();
    for h in r.heat.iter().filter(|h| h.body.is_some()) {
        let Some(name) = h.body.as_deref() else { continue };
        let body = s
            .design
            .iter()
            .flat_map(|d| &d.bodies)
            .find(|b| b.name == name)
            .ok_or_else(|| EngineError::Params(format!("wall heat: obstacle `{name}` is no longer in the design")))?;
        let index = body
            .face_index(&h.name)
            .ok_or_else(|| EngineError::Params(format!("wall heat: face `{}` of `{name}` is no longer in the design", h.name)))?;
        if let std::collections::btree_map::Entry::Vacant(entry) = areas.entry(name.to_owned()) {
            let mesh = body.tessellate((body.size()? * 0.001).max(1e-4))?;
            let mut by_face = BTreeMap::<u32, f64>::new();
            for (tri, face) in mesh.triangles.iter().zip(&mesh.tri_face) {
                let point = |k: u32| mesh.positions.get(k as usize).copied().map(Vec3::from);
                let (Some(a), Some(b), Some(c)) = (point(tri[0]), point(tri[1]), point(tri[2])) else { continue };
                *by_face.entry(*face).or_default() += 0.5 * (b - a).cross(c - a).len();
            }
            entry.insert(by_face);
        }
        let area = areas
            .get(name)
            .and_then(|a| a.get(&index))
            .copied()
            .and_then(|a| crate::units::geometry().to_si(fieldcraft_model::units::Quantity::Area, a).ok())
            .filter(|a| a.is_finite() && *a > 0.0)
            .ok_or_else(|| EngineError::Params(format!("wall heat: face `{}` of `{name}` has no finite positive area", h.name)))?;
        let power = interpolate(&r.history_times, &h.power, time);
        let flux = power.map(|q| q / area).filter(|q| q.is_finite());
        output.push(FaceHeat { body: name.into(), face: h.name.clone(), index, area, power, flux });
    }
    Ok(output)
}

/// Faces as [`faces`] lists them (JSON) with their areas, powers and fluxes in `display`.
pub fn show(faces: &mut serde_json::Value, display: &fieldcraft_model::DisplaySystem) {
    use crate::results::{Quantity, shown};
    for f in faces.as_array_mut().into_iter().flatten() {
        for (key, q) in [("area", Quantity::Area), ("power", Quantity::Power), ("flux", Quantity::HeatFlux)] {
            if let Some(v) = f.get(key).and_then(serde_json::Value::as_f64) {
                f[key] = serde_json::json!(shown(q, v, display));
            }
        }
    }
}

/// The units [`show`] gives faces in.
pub fn units(display: &fieldcraft_model::DisplaySystem) -> serde_json::Value {
    use fieldcraft_model::units::Quantity;
    serde_json::json!({"area": display.symbol(Quantity::Area), "power": display.symbol(Quantity::Power), "flux": display.symbol(Quantity::HeatFlux)})
}

/// One optional value per displayed surface triangle. Imported/unrecorded faces remain plain.
pub fn triangle_values(faces: &[FaceHeat], base: &DisplayMesh, bodies: &[&str]) -> Vec<Option<f64>> {
    let by_face: BTreeMap<_, _> = faces.iter().map(|f| ((f.body.as_str(), f.index), f.flux)).collect();
    base.surface_tris
        .iter()
        .map(|t| {
            let body = bodies.get(t.body as usize)?;
            by_face.get(&(*body, t.face)).copied().flatten()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn wall_heat_uses_face_identity_area_sign_and_frame_in_both_views() {
        let mut session = crate::flow::tests::plate_in_a_tank();
        let crate::StepResult::Fluid(mut r) = session.latest_result("Flow").unwrap().0.as_ref().clone() else { panic!() };
        let names: Vec<String> = session.design.as_ref().unwrap().bodies[0].faces().iter().map(|f| f.name.clone()).collect();
        let end = r.frames.last().unwrap().time;
        r.history_times = vec![end * 0.5, end];
        r.heat = names
            .iter()
            .enumerate()
            .map(|(k, name)| fieldcraft_lbm::HeatHistory {
                name: name.clone(),
                body: Some("Body1".into()),
                power: vec![0.0, if k == 0 { -120.0 } else { 0.0 }],
            })
            .collect();
        for f in &mut r.frames {
            f.temperature = vec![20.0; r.grid.len()];
        }
        let face_values = faces(&session, &r, Pick::default()).unwrap();
        assert!(face_values.iter().all(|f| f.area > 0.0));
        // Areas in m² (the design's 40 × 10 mm is 4e-4 m²).
        assert!(face_values.iter().any(|f| (f.area - 4e-4).abs() < 1e-14), "40 by 10 mm outer face");
        assert!(face_values.iter().any(|f| (f.area - 3e-4).abs() < 1e-14), "30 by 10 mm outer face");
        assert!(face_values[0].flux.unwrap() < 0.0);
        assert!((face_values.iter().map(|f| f.flux.unwrap() * f.area).sum::<f64>() + 120.0).abs() < 1e-10);
        assert!(face_values[1..].iter().all(|f| f.flux == Some(0.0)), "adiabatic faces are zero");
        session.store_result(crate::StepResult::Fluid(r));
        let report = session.execute("results.surface", &json!({"field": "wall_heat_flux"})).unwrap();
        assert_eq!(report["wall_heat"]["positive"], "into wall");
        // 120 W, shown in the default display system's mW.
        assert_eq!(report["wall_heat"]["faces"][0]["power"], -120_000.0);
        assert_eq!(report["wall_heat"]["units"]["power"], "mW");
        assert!(report["slice"].is_null(), "wall heat is a surface-only field");
        session.execute("results.streamlines", &json!({"show": false})).unwrap();
        let (scene, _, shown) = crate::view::scene(&session, &crate::view::Palette::LIGHT);
        assert!(shown.contour.as_ref().is_some_and(|c| c.0 == "wall_heat_flux"));
        assert!(scene.tris.windows(2).any(|t| t[0].color != t[1].color));
        let report = session.execute("results.frame", &json!({"frame": 0})).unwrap();
        assert!(report["wall_heat"]["faces"].as_array().unwrap().iter().all(|f| f["flux"].is_null()));
        assert!(session.execute("results.slice", &json!({"field": "wall_heat_flux"})).is_err());
        assert!(session.execute("results.streamlines", &json!({"field": "wall_heat_flux"})).is_err());
    }

    #[test]
    fn histories_interpolate_without_fabricating_initial_heat() {
        assert_eq!(interpolate(&[1.0, 3.0], &[-4.0, 8.0], 2.0), Some(2.0));
        assert_eq!(interpolate(&[1.0, 3.0], &[-4.0, 8.0], 1.0), Some(-4.0));
        assert_eq!(interpolate(&[1.0, 3.0], &[-4.0, 8.0], 0.0), None);
        assert_eq!(interpolate(&[1.0, 3.0], &[-4.0, 8.0], 4.0), None);
        assert_eq!(interpolate(&[1.0, 1.0], &[0.0, 1.0], 1.0), None);
        assert_eq!(interpolate(&[1.0], &[f64::NAN], 1.0), None);
        assert_eq!(interpolate(&[1.0], &[], 1.0), None);
    }
}
