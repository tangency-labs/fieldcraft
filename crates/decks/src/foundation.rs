//! Abaqus reference-surface foundations and sparse matrix fallback.
//! Public format: Abaqus 2025 Keywords *FOUNDATION, *MATRIX INPUT and
//! *MATRIX ASSEMBLE; Analysis Guide "Using Matrices". No vendor model is used.
use crate::{Error, Result};
use fieldcraft_model::{ElementKind, Model, Procedure, Resolved};
use std::collections::BTreeMap;

pub(crate) fn refuse(format: &'static str, model: &Model) -> Result<()> {
    if let Some((name, _)) = model.foundations().iter().next() {
        return Err(Error::Unsupported {
            format,
            what: format!("surface foundation `{name}`: equivalent matrix export is not implemented by this writer"),
        });
    }
    Ok(())
}
pub(crate) fn abaqus(model: &Model, r: &Resolved) -> Result<String> {
    let mut out = String::new();
    let mut native: BTreeMap<(u32, String), f64> = BTreeMap::new();
    for (index, (name, f)) in model.foundations().iter().enumerate() {
        let error = |reason: &str| Error::Unsupported { format: "Abaqus", what: format!("foundation `{name}`: {reason}") };
        if r.is_planar() {
            return Err(error("planar/axisymmetric foundations"));
        }
        // Matrix data stays fixed in a nonlinear step (Analysis Guide, "Using
        // Matrices"). This preserves our reference normal and area even for a
        // normal-only foundation; no follower support stiffness is introduced.
        let finite = model.steps().iter().any(|s| matches!(&s.procedure, Procedure::NonlinearStatic(n) if n.nonlinear_geometry));
        let matrix = finite || f.tangential != 0. || f.normal_damping != 0. || f.tangential_damping != 0. || f.loss_factor != 0.;
        for s in model.steps() {
            match &s.procedure {
                Procedure::Static {} | Procedure::NonlinearStatic(_) | Procedure::Frequency { .. } => (),
                Procedure::Harmonic(h) if !matrix || matches!(h.method, fieldcraft_model::harmonic::Method::Direct) => (),
                _ => {
                    return Err(error(
                        "only static, nonlinear static, frequency and direct harmonic matrix usage are supported; modal matrix usage requires SIM architecture",
                    ));
                }
            }
        }
        let set = r.sets.get(&f.set).ok_or_else(|| error("missing face set"))?;
        if !matrix {
            for face in &set.faces {
                let e = model.element(face.element).ok_or_else(|| error("missing element"))?;
                let kind = if matches!(e.kind, ElementKind::Tri3 | ElementKind::Tri6 | ElementKind::Quad4 | ElementKind::Quad8) {
                    "F".into()
                } else {
                    format!("F{}", face.face)
                };
                let value = native.entry((face.element, kind)).or_default();
                *value += f.normal;
                if !value.is_finite() {
                    return Err(error("overlapping foundation stiffness overflows"));
                }
            }
            continue;
        }
        for (operator, normal, tangential) in [
            ("STIFFNESS", f.normal, f.tangential),
            ("VISCOUS DAMPING", f.normal_damping, f.tangential_damping),
            ("STRUCTURAL DAMPING", f.normal * f.loss_factor, f.tangential * f.loss_factor),
        ] {
            if normal == 0. && tangential == 0. {
                continue;
            }
            let mut entries: BTreeMap<(u32, usize, u32, usize), f64> = BTreeMap::new();
            for face in &set.faces {
                let nodes = model.face_nodes(*face).ok_or_else(|| error("invalid face"))?;
                let xyz = nodes.iter().map(|n| model.label_xyz(*n).ok_or_else(|| error("missing node"))).collect::<Result<Vec<_>>>()?;
                let mut invalid = false;
                fieldcraft_model::nodal::integrate_face(&xyz, |shape, normal_area| {
                    let area = normal_area[0].hypot(normal_area[1]).hypot(normal_area[2]);
                    if area <= 0. || !area.is_finite() {
                        invalid = true;
                        return;
                    }
                    let n = normal_area.map(|x| x / area);
                    for ((&a, na), (&b, nb)) in nodes.iter().zip(shape).flat_map(|a| nodes.iter().zip(shape).map(move |b| (a, b))) {
                        for i in 0..3 {
                            for j in 0..3 {
                                if (a, i) < (b, j) {
                                    continue;
                                }
                                let v = na * nb * area * ((normal - tangential) * n[i] * n[j] + if i == j { tangential } else { 0. });
                                *entries.entry((a, i + 1, b, j + 1)).or_default() += v;
                            }
                        }
                    }
                })?;
                if invalid {
                    return Err(error("degenerate face"));
                }
            }
            let tag = format!("FC_FOUNDATION_{index}_{}", operator.replace(' ', "_"));
            out.push_str(&format!("*MATRIX INPUT, NAME={tag}, TYPE=SYMMETRIC\n"));
            for ((a, i, b, j), v) in entries {
                if !v.is_finite() {
                    return Err(error("matrix overflow"));
                }
                if v != 0. {
                    out.push_str(&format!("{a}, {i}, {b}, {j}, {v:.17e}\n"));
                }
            }
            out.push_str(&format!("*MATRIX ASSEMBLE, {operator}={tag}\n"));
        }
    }
    if !native.is_empty() {
        out.push_str("*FOUNDATION\n");
    }
    for ((element, kind), value) in native {
        out.push_str(&format!("{element}, {kind}, {value:.17e}\n"));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_model::*;
    fn model() -> Model {
        let mut m = Model::default();
        for (i, p) in [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]].into_iter().enumerate() {
            m.set_node(i as u32 + 1, p).unwrap();
        }
        m.set_element(1, ElementKind::Tri3, vec![1, 2, 3]).unwrap();
        m.add_set("all", Set { elements: vec![1], ..Set::default() }).unwrap();
        m.add_set("bed", Set { faces: vec![ElementFace { element: 1, face: 1 }], ..Set::default() }).unwrap();
        m.add_material("M", Material { elastic: Some(Elastic::Isotropic { young: 1000., poisson: 0.3 }), ..Default::default() }).unwrap();
        m.add_section(
            "S",
            Section {
                set: "all".into(),
                material: "M".into(),
                kind: SectionKind::Shell { thickness: 0.1, offset: 0., points: 3 },
                orientation: None,
                formulation: Default::default(),
            },
        )
        .unwrap();
        m.add_step(Step::new("S", Procedure::Static {}), None).unwrap();
        m.add_foundation(
            "Bed",
            Foundation { set: "bed".into(), normal: 100., tangential: 0., normal_damping: 0., tangential_damping: 0., loss_factor: 0. },
        )
        .unwrap();
        m
    }
    #[test]
    fn pure_normal_and_coupled_matrix_exports_preserve_coefficients() {
        let mut m = model();
        let r = m.resolve().unwrap();
        let text = abaqus(&m, &r).unwrap();
        assert!(text.contains("*FOUNDATION\n1, F, 1.00000000000000000e2"));
        let duplicate = m.foundation("Bed").unwrap().clone();
        m.add_foundation("Second", duplicate).unwrap();
        let doubled = abaqus(&m, &m.resolve().unwrap()).unwrap();
        assert!(doubled.contains("2.00000000000000000e2"));
        m.remove_foundation("Second").unwrap();
        let mut f = m.foundation("Bed").unwrap().clone();
        f.tangential = 50.;
        f.normal_damping = 10.;
        f.loss_factor = 0.2;
        m.replace_foundation("Bed", f).unwrap();
        let text = abaqus(&m, &m.resolve().unwrap()).unwrap();
        assert!(!text.contains("*FOUNDATION\n"));
        assert!(text.contains("VISCOUS DAMPING="));
        assert!(text.contains("STRUCTURAL DAMPING="));
        let mut sums = Vec::new();
        let mut sum = 0.;
        for line in text.lines() {
            if line.starts_with("*MATRIX INPUT") {
                sum = 0.;
            } else if line.starts_with("*MATRIX ASSEMBLE") {
                sums.push(sum);
            } else {
                let v: Vec<_> = line.split(',').map(str::trim).collect();
                if v.len() == 5 && v[1] == "3" && v[3] == "3" {
                    sum += v[4].parse::<f64>().unwrap() * if v[0] == v[2] { 1. } else { 2. };
                }
            }
        }
        for (a, b) in sums.iter().zip([50., 5., 10.]) {
            assert!((a - b).abs() < 1e-12);
        }
        assert!(refuse("Nastran", &m).unwrap_err().to_string().contains("Bed"));
    }

    #[test]
    fn nonlinear_foundation_export_keeps_the_reference_matrix() {
        for tangential in [0., 50.] {
            let mut m = model();
            let mut f = m.foundation("Bed").unwrap().clone();
            f.tangential = tangential;
            m.replace_foundation("Bed", f).unwrap();
            let linear = abaqus(&m, &m.resolve().unwrap()).unwrap();
            let mut step = m.step("S").unwrap().clone();
            step.procedure = Procedure::NonlinearStatic(NonlinearStatic { nonlinear_geometry: true, ..NonlinearStatic::default() });
            m.replace_step(step).unwrap();
            let finite = abaqus(&m, &m.resolve().unwrap()).unwrap();
            assert!(finite.contains("*MATRIX ASSEMBLE, STIFFNESS="));
            assert!(!finite.contains("*FOUNDATION\n"));
            if tangential != 0. {
                assert_eq!(linear, finite);
            }
            // A constant unit normal displacement of a triangle of area 1/2
            // produces total force k_n/2, including symmetric off-diagonals.
            let sum: f64 = finite
                .lines()
                .filter_map(|line| {
                    let v: Vec<_> = line.split(',').map(str::trim).collect();
                    (v.len() == 5 && v[1] == "3" && v[3] == "3").then(|| v[4].parse::<f64>().unwrap() * if v[0] == v[2] { 1. } else { 2. })
                })
                .sum();
            assert!((sum - 50.).abs() < 1e-12);
        }
    }
}
