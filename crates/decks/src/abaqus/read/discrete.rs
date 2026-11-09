//! Constant global MASS/ROTARYI and constant global SPRING1/2 and DASHPOT1/2 import.
//! Syntax: Abaqus public Keywords Reference, *MASS, *ROTARY INERTIA and *SPRING:
//! https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-mass.htm
//! https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-rotaryinertia.htm
//! https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-spring.htm
//! https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-dashpot.htm
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
    Inertia,
    Spring,
    Spring2,
    Dashpot,
    Dashpot2,
}
impl Kind {
    pub(super) fn node_count(self) -> usize {
        if matches!(self, Self::Spring2 | Self::Dashpot2) { 2 } else { 1 }
    }
    pub(super) fn parse(s: &str) -> Option<Self> {
        match canonical(s).as_str() {
            "MASS" => Some(Self::Mass),
            "ROTARYI" => Some(Self::Inertia),
            "SPRING1" => Some(Self::Spring),
            "SPRING2" => Some(Self::Spring2),
            "DASHPOT1" => Some(Self::Dashpot),
            "DASHPOT2" => Some(Self::Dashpot2),
            _ => None,
        }
    }
}
pub(super) struct Record {
    pub kind: Kind,
    pub node: u32,
    pub second: Option<u32>,
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
        let rotary = b.keyword == "ROTARYINERTIA";
        b.only(if mass { &["ELSET", "TYPE"] } else { &["ELSET"] })?;
        let entries = values(b)?;
        let mut point_mass = None;
        let mut inertia = None;
        let (kind, dof, value) = if mass {
            if b.rows.len() != 1 {
                return Err(error(line, "constant MASS requires exactly one data row"));
            }
            let convert = |v| quantity(units, Quantity::Mass, v, line);
            let typ = if b.params.contains_key("TYPE") { b.required("TYPE")? } else { "ISOTROPIC" };
            point_mass = Some(match (canonical(typ).as_str(), entries.as_slice()) {
                ("ISOTROPIC", [m]) => Mass::Isotropic(convert(*m)?),
                ("ANISOTROPIC", [x, y, z]) => Mass::Directional([convert(*x)?, convert(*y)?, convert(*z)?]),
                _ => {
                    return Err(error(
                        line,
                        "constant MASS needs TYPE=ISOTROPIC with one value or TYPE=ANISOTROPIC with three global principal values",
                    ));
                }
            });
            (Kind::Mass, 0, 0.0)
        } else if rotary {
            if b.rows.len() != 1 || !matches!(entries.len(), 3 | 6) {
                return Err(error(line, "ROTARY INERTIA needs one row: I11,I22,I33 and optional I12,I13,I23"));
            }
            let mut tensor = [0.0; 6];
            for (value, input) in tensor.iter_mut().zip(&entries) {
                *value = quantity(units, Quantity::RotaryInertia, *input, line)?;
            }
            inertia = Some(tensor);
            (Kind::Inertia, 0, 0.0)
        } else {
            if b.rows.len() != 2 {
                return Err(error(line, format!("{} requires exactly a DOF row and a constant coefficient row", b.keyword)));
            }
            let first = b.rows.first().ok_or_else(|| error(line, "missing discrete degree of freedom"))?;
            let fields = row_fields(first.text, first.line)?;
            if !matches!(fields.len(), 1 | 2) || entries.len() != fields.len() + 1 {
                return Err(error(
                    line,
                    format!(
                        "{} requires one or two DOFs and one constant coefficient; frequency/temperature/field dependence is not mapped",
                        b.keyword
                    ),
                ));
            }
            let dof = label(fields.first().map(String::as_str).unwrap_or_default(), first.line)?;
            if !(1..=6).contains(&dof) {
                return Err(error(first.line, "spring/dashpot degree of freedom must be 1 through 6"));
            }
            let pair = fields.len() == 2;
            if pair && fields.get(1).map(|v| label(v, first.line)).transpose()? != Some(dof) {
                return Err(error(
                    first.line,
                    format!("{}2 different endpoint DOFs are not expressible by the native shared global direction", b.keyword),
                ));
            }
            let dashpot = b.keyword == "DASHPOT";
            let q = match (dashpot, dof <= 3) {
                (false, true) => Quantity::Stiffness,
                (false, false) => Quantity::RotationalStiffness,
                (true, true) => Quantity::Damping,
                (true, false) => Quantity::RotationalDamping,
            };
            let kind = match (dashpot, pair) {
                (false, false) => Kind::Spring,
                (false, true) => Kind::Spring2,
                (true, false) => Kind::Dashpot,
                (true, true) => Kind::Dashpot2,
            };
            (kind, dof as u8, quantity(units, q, *entries.last().ok_or_else(|| error(line, "missing coefficient"))?, line)?)
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
            let mut targets = Vec::new();
            for node in std::iter::once(record.node).chain(record.second) {
                if model.node(node).is_none() {
                    return Err(error(record.line, format!("discrete element {id} names missing node {node}")));
                }
                // Bounded collision-free target names, including both endpoints.
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
                model.add_set(&target, Set { nodes: vec![node], ..Set::default() }).map_err(|e| error(line, e.to_string()))?;
                targets.push(Target::Set(target));
            }
            let target = targets.first().cloned().ok_or_else(|| error(line, "missing discrete target"))?;
            let second = targets.get(1).cloned();
            let name = format!("InpDiscreteElement{id}");
            if mass || rotary {
                model.add_mass(&name, PointMass { target, mass: point_mass, inertia })
            } else if b.keyword == "DASHPOT" {
                model.add_spring(&name, Spring::Dashpot { first: target, second, direction: Direction::Dof(dof), coefficient: value })
            } else {
                model.add_spring(&name, Spring::Spring { first: target, second, direction: Direction::Dof(dof), stiffness: value })
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

    #[test]
    fn rotary_inertia_keeps_signed_global_tensor_and_units() {
        let input = "*NODE\n7,0,0,0\n*ELEMENT,TYPE=ROTARYI,ELSET=J\n10,7\n*ROTARY INERTIA,ELSET=J\n100,200,250,10,-20,30\n";
        let r = read_in(input, Preset::MmTonneS).unwrap();
        assert_eq!(r.model.masses().values().next().unwrap().inertia, Some([0.1, 0.2, 0.25, 0.01, -0.02, 0.03]));
        for (a, b) in [
            ("ELSET=J\n100", "ELSET=J,ORIENTATION=Local\n100"),
            ("ELSET=J\n100", "ELSET=J,ALPHA=1\n100"),
            ("ELSET=J\n100", "ELSET=J,COMPOSITE=.05\n100"),
            ("100,200,250,10,-20,30", "1,1,1,2,0,0"),
        ] {
            assert!(read_in(&input.replace(a, b), Preset::MmTonneS).is_err());
        }
    }

    #[test]
    fn two_point_properties_keep_endpoints_units_and_refuse_distinct_dofs() {
        let source = "*NODE\n1,0,0,0\n2,1,0,0\n*ELEMENT,TYPE=SPRING2,ELSET=K\n10,1,2\n*ELEMENT,TYPE=DASHPOT2,ELSET=C\n11,2,1\n*SPRING,ELSET=K\n1,1\n200\n*DASHPOT,ELSET=C\n1,1\n6\n*NSET,NSET=Ends,ELSET=K\n";
        let mut model = read_in(source, Preset::MmTonneS).unwrap().model;
        model.add_step(fieldcraft_model::Step::new("Check", fieldcraft_model::Procedure::Static {}), None).unwrap();
        assert_eq!(model.set("Ends").unwrap().nodes, vec![1, 2]);
        let resolved = model.resolve().unwrap();
        assert!(resolved.springs.iter().any(|s| !s.dashpot && s.first == 1 && s.second == Some(2) && s.value == 200000.));
        assert!(resolved.springs.iter().any(|s| s.dashpot && s.first == 2 && s.second == Some(1) && s.value == 6000.));
        for text in [
            source.replace("1,1\n200", "1,2\n200"),
            source.replace("1,1\n6", "1,2\n6"),
            source.replace("11,2,1", "11,2,99"),
            source.replace("11,2,1", "11,2,2"),
            source.replace("*DASHPOT,ELSET=C", "*DASHPOT,ELSET=C,NONLINEAR"),
            source.replace("1,1\n6", "1,1\n6,10"),
        ] {
            assert!(read_in(&text, Preset::MKgS).is_err(), "{text}");
        }
    }

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
            SOURCE.replace("*MASS,ELSET=Mass", "*MASS,ELSET=Mass,TYPE"),
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
    fn global_directional_mass_preserves_axes_and_refuses_unmapped_options() {
        let input = SOURCE.replace("*MASS,ELSET=Mass\n2", "*MASS,ELSET=Mass,TYPE=ANISOTROPIC\n2,3,5");
        let got = read_analysis_limited(&input, Preset::MKgS, Limits::default()).unwrap();
        assert!(got.analysis_complete);
        assert_eq!(got.model.masses().values().next().unwrap().mass, Some(Mass::Directional([2., 3., 5.])));
        for (from, to) in [("TYPE=ANISOTROPIC", "TYPE=ANISOTROPIC,ORIENTATION=Axes"), ("2,3,5", "2,3,5,20"), ("2,3,5", "2,-3,5")] {
            assert!(read_in(&input.replace(from, to), Preset::MKgS).is_err());
        }
        let changed = input.replace("*FREQUENCY\n1", "*FREQUENCY\n1\n*MASS,ELSET=Mass,TYPE=ANISOTROPIC\n3,4,6");
        let denied = read_analysis_limited(&changed, Preset::MKgS, Limits::default()).unwrap();
        assert!(!denied.analysis_complete);
        assert_eq!(denied.model.masses().values().next().unwrap().mass, Some(Mass::Directional([2., 3., 5.])));
        assert!(denied.unknown.iter().any(|u| u.keyword == "MASS"));
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
