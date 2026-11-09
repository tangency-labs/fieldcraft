//! Constant isotropic MASS and grounded SPRING1 import.
//! Syntax: Abaqus public Keywords Reference, *MASS and *SPRING:
//! https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-mass.htm
//! https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-spring.htm
//! Discrete labels become native definitions, not fictitious mesh elements. Original labels
//! and ELSET memberships remain in the retained source; generated target sets are explicit.
use super::{Block, canonical, check_work, error, key, label, quantity, row_fields, values};
use crate::Result;
use fieldcraft_model::units::Quantity;
use fieldcraft_model::{Direction, DisplaySystem, Mass, Model, PointMass, Set, Spring, Target};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum Kind {
    Mass,
    Spring,
}
impl Kind {
    pub(super) fn parse(s: &str) -> Option<Self> {
        match canonical(s).as_str() {
            "MASS" => Some(Self::Mass),
            "SPRING1" => Some(Self::Spring),
            _ => None,
        }
    }
}
pub(super) struct Record {
    pub kind: Kind,
    pub node: u32,
    pub line: usize,
}

pub(super) fn apply(
    model: &mut Model,
    records: &BTreeMap<u32, Record>,
    properties: &[&Block<'_>],
    elsets: &BTreeMap<String, Vec<u32>>,
    units: &DisplaySystem,
) -> Result<()> {
    let mut assigned = BTreeSet::new();
    let mut names: BTreeSet<String> = model.sets().keys().map(|s| s.to_ascii_uppercase()).collect();
    let mut serial = 0usize;
    for b in properties {
        let line = b.line;
        let mass = b.keyword == "MASS";
        b.only(if mass { &["ELSET", "TYPE"] } else { &["ELSET"] })?;
        if mass && b.params.contains_key("TYPE") && canonical(b.required("TYPE")?) != "ISOTROPIC" {
            return Err(error(line, "only constant isotropic MASS is supported"));
        }
        let entries = values(b)?;
        let (kind, dof, value) = if mass {
            let [value] = entries.as_slice() else {
                return Err(error(line, "constant isotropic MASS needs exactly one value"));
            };
            (Kind::Mass, 0, quantity(units, Quantity::Mass, *value, line)?)
        } else {
            if b.rows.len() != 2 {
                return Err(error(line, "constant SPRING1 needs a degree-of-freedom row and a stiffness row"));
            }
            let first = b.rows.first().ok_or_else(|| error(line, "missing spring degree of freedom"))?;
            let fields = row_fields(first.text, first.line)?;
            if fields.len() != 1 || entries.len() != 2 {
                return Err(error(
                    line,
                    "SPRING1 requires one degree of freedom and one constant stiffness; frequency/temperature dependence is not mapped",
                ));
            }
            let dof = label(fields.first().map(String::as_str).unwrap_or_default(), first.line)?;
            if !(1..=6).contains(&dof) {
                return Err(error(first.line, "spring degree of freedom must be 1 through 6"));
            }
            let q = if dof <= 3 { Quantity::Stiffness } else { Quantity::RotationalStiffness };
            (Kind::Spring, dof as u8, quantity(units, q, *entries.get(1).ok_or_else(|| error(line, "missing stiffness"))?, line)?)
        };
        let set = key(b.required("ELSET")?, line)?;
        let labels = elsets.get(&set).filter(|v| !v.is_empty()).ok_or_else(|| error(line, format!("missing/empty discrete element set `{set}`")))?;
        for (index, id) in labels.iter().enumerate() {
            check_work(index)?;
            let record =
                records.get(id).filter(|r| r.kind == kind).ok_or_else(|| error(line, format!("discrete property does not match element {id}")))?;
            if !assigned.insert(*id) {
                return Err(error(line, format!("overlapping discrete properties for element {id}")));
            }
            if model.node(record.node).is_none() {
                return Err(error(record.line, format!("discrete element {id} names missing node {}", record.node)));
            }
            // Each iteration consumes either an existing name or a new definition. The bound
            // includes all input sets, so hostile names cannot cause an unbounded search.
            let mut target = None;
            for _ in 0..=names.len() {
                serial = serial.checked_add(1).ok_or_else(|| error(line, "discrete name limit reached"))?;
                check_work(serial)?;
                let name = format!("InpDiscreteTarget{serial}");
                if names.insert(name.to_ascii_uppercase()) {
                    target = Some(name);
                    break;
                }
            }
            let target = target.ok_or_else(|| error(line, "no free discrete target name"))?;
            model.add_set(&target, Set { nodes: vec![record.node], ..Set::default() }).map_err(|e| error(line, e.to_string()))?;
            let target = Target::Set(target);
            let name = format!("InpDiscreteElement{id}");
            if mass {
                model.add_mass(&name, PointMass { target, mass: Some(Mass::Isotropic(value)), inertia: None })
            } else {
                model.add_spring(&name, Spring::Spring { first: target, second: None, direction: Direction::Dof(dof), stiffness: value })
            }
            .map_err(|e| error(line, e.to_string()))?;
        }
    }
    if let Some((id, record)) = records.iter().find(|(id, _)| !assigned.contains(*id)) {
        return Err(error(record.line, format!("discrete element {id} has no property")));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::abaqus::read::{Limits, read_analysis_limited, read_in};
    use fieldcraft_model::units::Preset;

    // Original one-mass oscillator: k=200 N/m, m=2 kg, f=sqrt(k/m)/(2*pi).
    const SOURCE: &str = "*NODE\n7,0,0,0\n*ELEMENT,TYPE=MASS,ELSET=Mass\n10,7\n*ELEMENT,TYPE=SPRING1,ELSET=Spring\n11,7\n*NSET,NSET=InpDiscreteTarget1\n7\n*NSET,NSET=FromMass,ELSET=Mass\n*MASS,ELSET=Mass\n2\n*SPRING,ELSET=Spring\n1\n200\n*BOUNDARY\n7,2,3\n*STEP,NAME=Modes\n*FREQUENCY\n1\n*END STEP\n";

    #[test]
    fn oscillator_maps_properties_references_and_normalized_roundtrip() {
        let got = read_analysis_limited(SOURCE, Preset::MKgS, Limits::default()).unwrap();
        assert!(got.analysis_complete, "{:?}", got.warnings);
        assert!(got.model.elements().is_empty());
        assert_eq!(got.model.set("FromMass").unwrap().nodes, vec![7]);
        let resolved = got.model.resolve().unwrap();
        assert_eq!(resolved.masses.len(), 1);
        assert_eq!(resolved.springs.len(), 1);
        assert_eq!(got.model.mass("InpDiscreteElement10").unwrap().mass, Some(Mass::Isotropic(2.0)));
        assert_eq!(got.model.set("InpDiscreteTarget1").unwrap().nodes, vec![7]);
        let normalized = crate::abaqus::write(&got.model, Preset::MKgS).unwrap();
        let back = read_analysis_limited(&normalized, Preset::MKgS, Limits::default()).unwrap();
        assert!(back.analysis_complete, "{:?}", back.warnings);
        assert_eq!(back.model.resolve().unwrap().masses.len(), 1);
        assert_eq!(back.model.resolve().unwrap().springs.len(), 1);
        assert_eq!(back.model.masses().values().next().unwrap().mass, Some(Mass::Isotropic(2.0)));
    }

    #[test]
    fn discrete_properties_convert_units_and_reject_missing_or_ambiguous_semantics() {
        let m = read_in(SOURCE, Preset::MmTonneS).unwrap().model;
        assert_eq!(m.masses().values().next().unwrap().mass, Some(Mass::Isotropic(2000.0)));
        let Spring::Spring { stiffness, .. } = m.springs().values().next().unwrap() else { panic!("expected spring") };
        assert_eq!(*stiffness, 200000.0);
        for text in [
            SOURCE.replace("*MASS,ELSET=Mass\n2\n", ""),
            SOURCE.replace("*MASS,ELSET=Mass", "*MASS,ELSET=Mass,ALPHA=0.1"),
            SOURCE.replace("*SPRING,ELSET=Spring", "*SPRING,ELSET=Spring,NONLINEAR"),
            SOURCE.replace("1\n200\n", "1\n200,5\n"),
            SOURCE.replace("10,7", "10,99"),
            SOURCE.replace("11,7", "10,7"),
            SOURCE.replace("*BOUNDARY", "*MASS,ELSET=Mass\n3\n*BOUNDARY"),
        ] {
            assert!(read_in(&text, Preset::MKgS).is_err(), "{text}");
        }
        assert!(super::super::read_limited(SOURCE, Preset::MKgS, Limits { elements: 1, ..Limits::default() }).is_err());
    }

    #[test]
    fn gravity_does_not_silently_omit_discrete_mass() {
        let source = SOURCE.replace("*FREQUENCY\n1\n", "*STATIC\n1,1\n*DLOAD\n,GRAV,9.81,1,0,0\n");
        let got = read_analysis_limited(&source, Preset::MKgS, Limits::default()).unwrap();
        assert!(!got.analysis_complete);
        assert!(got.model.steps().is_empty());
        assert!(got.warnings.iter().any(|w| w.contains("discrete load membership")), "{:?}", got.warnings);
    }

    #[test]
    fn translated_instance_keeps_discrete_references_and_rotated_spring_is_refused() {
        let part = SOURCE.split("*BOUNDARY").next().unwrap();
        let source = format!(
            "*PART,NAME=Oscillator\n{part}*END PART\n*ASSEMBLY,NAME=A\n*INSTANCE,NAME=One,PART=Oscillator\n1,2,3\n*END INSTANCE\n*END ASSEMBLY\n"
        );
        let got = read_in(&source, Preset::MKgS).unwrap();
        assert_eq!(got.instances.len(), 1);
        assert_eq!(got.model.masses().len(), 1);
        assert_eq!(got.model.springs().len(), 1);
        let rotated = source.replace("1,2,3\n*END INSTANCE", "1,2,3\n0,0,0,0,0,1,90\n*END INSTANCE");
        assert!(read_in(&rotated, Preset::MKgS).unwrap_err().to_string().contains("rotated instance"));
    }
}
