//! Flat STEP export of the current world-space B-rep bodies. The pinned SolidCraft writer
//! returns the exact ADVANCED_FACE entity for each source record face, including cavity faces.

use std::collections::BTreeMap;

use solidcraft_kernel::{ExportBody, KernelError, StepHeader, guard, step_export_bodies_with_face_ids, step_string};

use crate::{Result, Shape};

/// One placed body, its STEP solid name, and boundary names keyed by current FieldCraft face
/// index. Any face without an override keeps its persistent name, when it has one.
pub struct StepBody<'a> {
    pub shape: &'a Shape,
    pub name: &'a str,
    pub face_names: BTreeMap<u32, &'a str>,
}

/// Write the bodies as one AP242 product with their occurrence placements already applied.
/// A virtual merge labels each of its original faces. A virtual split has no B-rep edge to
/// delimit its two regions in STEP and is refused.
pub fn export_named_step(bodies: &[StepBody<'_>]) -> Result<String> {
    Ok(guard("export named STEP", || {
        let mut export = Vec::with_capacity(bodies.len());
        for body in bodies {
            if !body.shape.cuts.is_empty() {
                return Err(KernelError::Invalid(format!("body `{}` has a virtual face split that cannot be represented in STEP", body.shape.name)));
            }
            check_label(body.name)?;
            export.push(ExportBody { name: body.name.into(), body: &body.shape.body, color: None });
        }
        let header = StepHeader { file_name: "FieldCraft.step".into(), description: "FieldCraft analysis geometry".into(), ..StepHeader::default() };
        let written = step_export_bodies_with_face_ids(&export, &header)?;
        if written.face_ids.len() != bodies.len() {
            return Err(KernelError::Failed("STEP writer returned an incomplete body map".into()));
        }
        let mut labels = BTreeMap::new();
        for (body, ids) in bodies.iter().zip(&written.face_ids) {
            collect_labels(body, ids, &mut labels)?;
        }
        label_entities(&written.text, labels)
    })?)
}

fn check_label(name: &str) -> solidcraft_kernel::Result<()> {
    if name.trim().is_empty() || name.chars().any(char::is_control) {
        return Err(KernelError::Invalid("STEP region names must contain text and cannot include control characters".into()));
    }
    Ok(())
}

fn collect_labels<'a>(body: &'a StepBody<'_>, ids: &[u64], labels: &mut BTreeMap<u64, &'a str>) -> solidcraft_kernel::Result<()> {
    let shape = body.shape;
    let current_count = if shape.faces.is_empty() { ids.len() } else { shape.faces.len() };
    if let Some((index, _)) = body.face_names.last_key_value().filter(|(index, _)| **index as usize >= current_count) {
        return Err(KernelError::Invalid(format!("body `{}` has no face at index {index}", shape.name)));
    }
    let count = if shape.face_map.is_empty() { current_count } else { shape.face_map.len() };
    if count != ids.len() {
        return Err(KernelError::Failed(format!("STEP writer returned the wrong number of faces for body `{}`", shape.name)));
    }
    for (raw, &id) in ids.iter().enumerate() {
        let index = if let Some(&mapped) = shape.face_map.get(raw) {
            mapped
        } else {
            u32::try_from(raw).map_err(|_| KernelError::Invalid("too many faces for STEP export".into()))?
        };
        let persistent = shape.faces.get(index as usize).map(|face| face.name.as_str()).filter(|s| !s.is_empty());
        if let Some(name) = body.face_names.get(&index).copied().or(persistent) {
            check_label(name)?;
            if id == 0 || labels.insert(id, name).is_some() {
                return Err(KernelError::Failed(format!("body `{}`, face {index} has no unique written STEP entity", shape.name)));
            }
        }
    }
    Ok(())
}

/// Only process the kernel writer's own output: one entity per line, a blank first argument
/// on each ADVANCED_FACE. Check every requested id and copy every other byte verbatim.
fn label_entities(text: &str, mut labels: BTreeMap<u64, &str>) -> solidcraft_kernel::Result<String> {
    let mut output_size = text.len();
    let too_large = || KernelError::Invalid("named STEP file exceeds the supported file size".into());
    if output_size > solidcraft_io::MAX_STEP_BYTES {
        return Err(too_large());
    }
    let mut named = String::with_capacity(output_size);
    let mut literals = BTreeMap::new();
    for line in text.split_inclusive('\n') {
        let entity = line.strip_prefix('#').and_then(|l| l.split_once('='));
        if let Some((number, tail)) = entity
            && let Ok(id) = number.parse::<u64>()
            && let Some(name) = labels.remove(&id)
        {
            let rest =
                tail.strip_prefix("ADVANCED_FACE('',").ok_or_else(|| KernelError::Failed(format!("STEP entity #{id} is not a writable face")))?;
            let literal = literals.entry(name).or_insert_with(|| step_string(name));
            output_size = output_size.saturating_add(literal.len().saturating_sub(2));
            if output_size > solidcraft_io::MAX_STEP_BYTES {
                return Err(too_large());
            }
            named.push('#');
            named.push_str(&id.to_string());
            named.push_str("=ADVANCED_FACE(");
            named.push_str(literal);
            named.push(',');
            named.push_str(rest);
        } else {
            named.push_str(line);
        }
    }
    if let Some((id, _)) = labels.first_key_value() {
        return Err(KernelError::Failed(format!("STEP writer omitted named face #{id}")));
    }
    Ok(named)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn placed_bodies_keep_their_positions_and_independent_face_labels() {
        let bytes = crate::sample::step_assembly(10.0).unwrap();
        let (design, _) = crate::open_geometry("pair.step", &bytes).unwrap();
        let names = ["air", "steel"];
        let labels = ["outer", "interface"];
        let bodies: Vec<_> = design
            .bodies
            .iter()
            .zip(names)
            .zip(labels)
            .map(|((shape, name), label)| StepBody { shape, name, face_names: [(0, label), (2, label)].into() })
            .collect();
        let step = export_named_step(&bodies).unwrap();
        assert!(!step.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE("));
        for label in labels {
            assert_eq!(step.matches(&format!("=ADVANCED_FACE('{label}',")).count(), 2);
        }
        let read = guard("read named STEP for test", || solidcraft_kernel::step_import(&step)).unwrap();
        assert_eq!(read.bodies.len(), 2);
        for ((before, after), name) in design.bodies.iter().zip(read.bodies).zip(names) {
            assert_eq!(after.name, name);
            let measure = guard("measure imported solid", || solidcraft_kernel::measure(&after.body)).unwrap();
            assert!((measure.centroid.x - before.measure().unwrap().centroid[0]).abs() < 1e-8);
        }
    }

    #[test]
    fn merged_faces_keep_their_boundary_and_step_strings_are_escaped() {
        let (design, _) = crate::open_geometry("pair.step", &crate::sample::step_assembly(10.0).unwrap()).unwrap();
        let source = &design.bodies[0];
        let members: Vec<_> = source.faces().iter().map(|face| face.name.clone()).collect();
        let shape = source.merge_faces("enclosure", &members).unwrap();
        let name = r"air's\ä";
        let one = |index, label| export_named_step(&[StepBody { shape: &shape, name: "part", face_names: [(index, label)].into() }]);
        let step = one(0, name).unwrap();
        assert_eq!(step.matches(r"=ADVANCED_FACE('air''s\\\X2\00E4\X0\',").count(), 6);
        assert!(one(99, name).is_err());
        assert!(one(0, "line\nbreak").is_err());
    }
}
