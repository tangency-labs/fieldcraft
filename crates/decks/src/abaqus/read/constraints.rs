//! Constraint mapping from public *TIE, *COUPLING and *KINEMATIC keyword syntax.
//! Numerical enforcement is FieldCraft's native mortar / rigid coupling, not vendor coefficients.
use super::{Block, Limits, MAX_RECORDS, Read, canonical, error, history, lex, number, row_fields};
use crate::{Error, Result};
use fieldcraft_model::units::{Preset, Quantity};
use fieldcraft_model::{Coupling, CouplingKind, DisplaySystem, Dofs, Model, ReferencePoint, Set, Tie};
use std::collections::BTreeSet;

fn set_name(model: &Model, value: &str, line: usize) -> Result<String> {
    let name = lex::unquote(value, line)?;
    model.sets().keys().find(|n| n.eq_ignore_ascii_case(&name)).cloned().ok_or_else(|| error(line, format!("missing constraint set/surface {name}")))
}
fn nodes(model: &Model, value: &str, line: usize) -> Result<Vec<u32>> {
    if let Ok(id) = value.parse::<u32>() {
        return model.node(id).map(|_| vec![id]).ok_or_else(|| error(line, format!("missing constraint node {id}")));
    }
    let name = set_name(model, value, line)?;
    let set = model.set(&name).ok_or_else(|| error(line, "missing constraint set"))?;
    if set.nodes.is_empty() {
        return Err(history::unsupported(line, format!("constraint target {name} requires a node-only set")));
    }
    Ok(set.nodes.clone())
}
fn tie(b: &Block<'_>, model: &mut Model, preset: Preset) -> Result<()> {
    history::allow(b, &["NAME", "ADJUST", "POSITIONTOLERANCE", "TYPE"])?;
    if canonical(history::value(b, "ADJUST", "YES")?) != "NO" {
        return Err(history::unsupported(b.line, "*TIE ADJUST (requires explicit NO)"));
    }
    if canonical(history::value(b, "TYPE", "SURFACETOSURFACE")?) != "SURFACETOSURFACE" {
        return Err(history::unsupported(b.line, "*TIE TYPE (requires SURFACE TO SURFACE)"));
    }
    let Some(tolerance) = b.value("POSITIONTOLERANCE") else {
        return Err(history::unsupported(b.line, "*TIE default POSITION TOLERANCE; specify a positive length"));
    };
    if b.rows.len() != 1 {
        return Err(history::unsupported(b.line, "*TIE requires one surface pair per named constraint"));
    }
    let row = b.rows.first().ok_or_else(|| error(b.line, "missing tie pair"))?;
    let pair = row_fields(row.text, row.line)?;
    if pair.len() != 2 {
        return Err(error(row.line, "TIE needs secondary, main surface"));
    }
    let mut names = Vec::new();
    for value in pair {
        let name = set_name(model, &value, row.line)?;
        let set = model.set(&name).ok_or_else(|| error(row.line, "missing tie surface"))?;
        if set.faces.is_empty() || set.faces.iter().any(|f| model.element(f.element).is_none_or(|e| e.kind.dimension() != 3)) {
            return Err(history::unsupported(row.line, format!("*TIE surface {name}: only solid element faces")));
        }
        names.push(name);
    }
    let [secondary, main]: [String; 2] = names.try_into().map_err(|_| error(row.line, "invalid tie pair"))?;
    let tolerance = super::quantity(&DisplaySystem::new(preset), Quantity::Length, number(tolerance, b.line)?, b.line)?;
    model
        .add_tie(&history::name(b.required("NAME")?, b.line)?, Tie { main, secondary, position_tolerance: Some(tolerance) })
        .map_err(|e| error(b.line, e.to_string()))
}

fn coupling(b: &Block<'_>, sub: Option<&Block<'_>>, blocks: &[Block<'_>], model: &mut Model, limits: Limits) -> Result<Vec<usize>> {
    history::allow(b, &["CONSTRAINTNAME", "REFNODE", "SURFACE"])?;
    if !b.rows.is_empty() {
        return Err(error(b.line, "COUPLING has no data rows"));
    }
    let sub = sub.ok_or_else(|| error(b.line, "COUPLING needs KINEMATIC"))?;
    let kind = match sub.keyword.as_str() {
        "KINEMATIC" => {
            history::allow(sub, &[])?;
            CouplingKind::Kinematic
        }
        "DISTRIBUTING" => {
            history::allow(sub, &["COUPLING", "WEIGHTINGMETHOD", "ROTATIONALCOUPLING"])?;
            if canonical(history::value(sub, "COUPLING", "CONTINUUM")?) != "CONTINUUM" {
                return Err(history::unsupported(sub.line, "*DISTRIBUTING COUPLING=STRUCTURAL"));
            }
            if canonical(history::value(sub, "WEIGHTINGMETHOD", "UNIFORM")?) != "UNIFORM" {
                return Err(history::unsupported(sub.line, "*DISTRIBUTING non-UNIFORM WEIGHTING METHOD"));
            }
            if sub.params.contains_key("ROTATIONALCOUPLING") && canonical(history::value(sub, "ROTATIONALCOUPLING", "")?) != "CONTINUUM" {
                return Err(history::unsupported(sub.line, "*DISTRIBUTING ROTATIONAL COUPLING=STRUCTURAL"));
            }
            CouplingKind::Distributing
        }
        _ => return Err(history::unsupported(sub.line, format!("*COUPLING subcard *{}", sub.keyword))),
    };
    let mut bits = if sub.rows.is_empty() { 63 } else { 0 };
    for row in &sub.rows {
        let f = row_fields(row.text, row.line)?;
        if f.len() > 2 {
            return Err(error(row.line, "coupling DOF row has more than two fields"));
        }
        if f.first().is_none_or(|v| v.is_empty()) {
            bits = 63;
            continue;
        }
        let start =
            f.first().ok_or_else(|| error(row.line, "missing kinematic DOF"))?.parse::<u8>().map_err(|_| error(row.line, "invalid kinematic DOF"))?;
        let end = f.get(1).filter(|v| !v.is_empty()).map_or(Ok(start), |v| v.parse::<u8>().map_err(|_| error(row.line, "invalid kinematic DOF")))?;
        if f.len() > 2 || start < 1 || end > 6 || end < start {
            return Err(error(row.line, "kinematic DOF range must be within 1..6"));
        }
        for d in start..=end {
            bits |= 1 << (d - 1);
        }
    }
    if (kind == CouplingKind::Distributing && bits != 63 && bits != 7) || bits & 7 == 0 {
        return Err(history::unsupported(
            sub.line,
            format!(
                "*{} DOFs: distributing requires explicit 1..3 or 1..6; partial rotational release is not mapped; kinematic requires at least one translation",
                sub.keyword
            ),
        ));
    }
    let labels = nodes(model, b.required("REFNODE")?, b.line)?;
    let [label]: [u32; 1] = labels.try_into().map_err(|_| history::unsupported(b.line, "*COUPLING REF NODE requires a singleton"))?;
    if model.elements().values().any(|e| e.nodes.contains(&label)) {
        return Err(history::unsupported(b.line, "*COUPLING REF NODE belongs to mesh connectivity"));
    }
    let contains_reference =
        |target: &fieldcraft_model::Target| target.set().and_then(|name| model.set(name)).is_some_and(|set| set.nodes.contains(&label));
    if model.masses().values().any(|m| contains_reference(&m.target)) || model.springs().values().any(|s| s.targets().any(contains_reference)) {
        return Err(history::unsupported(b.line, "*COUPLING REF NODE with MASS/SPRING target requires discrete reference-point mapping"));
    }
    let surf = history::name(b.required("SURFACE")?, b.line)?;
    let surfaces: Vec<_> = blocks
        .iter()
        .filter(|s| s.keyword == "SURFACE" && s.value("NAME").is_some_and(|v| history::name(v, s.line).is_ok_and(|v| v == surf)))
        .collect();
    let [surface]: [&Block<'_>; 1] = surfaces.try_into().map_err(|_| error(b.line, format!("missing or duplicate coupling surface {surf}")))?;
    history::allow(surface, &["NAME", "TYPE"])?;
    let element_surface = match canonical(history::value(surface, "TYPE", "ELEMENT")?).as_str() {
        "NODE" => false,
        "ELEMENT" if kind == CouplingKind::Kinematic => true,
        "ELEMENT" => return Err(history::unsupported(surface.line, "*DISTRIBUTING element-surface tributary weights are not mapped")),
        _ => return Err(history::unsupported(surface.line, "*COUPLING surface TYPE requires NODE or kinematic ELEMENT")),
    };
    let mut members = BTreeSet::new();
    if element_surface {
        let name = set_name(model, &surf, surface.line)?;
        let resolved = model.resolve_set(&name).map_err(|e| error(surface.line, e.to_string()))?;
        if resolved.nodes.len() > limits.members {
            return Err(error(surface.line, "coupling member limit reached"));
        }
        members.extend(resolved.nodes);
    } else {
        for row in &surface.rows {
            let f = row_fields(row.text, row.line)?;
            if f.len() > 2 || f.is_empty() {
                return Err(error(row.line, "node surface needs node/set and optional weight"));
            }
            if let Some(weight) = f.get(1)
                && number(weight, row.line)? != 1.0
            {
                return Err(history::unsupported(row.line, "*SURFACE non-unit node weight"));
            }
            for n in nodes(model, f.first().map(String::as_str).unwrap_or_default(), row.line)? {
                if members.len() >= limits.members && !members.contains(&n) {
                    return Err(error(row.line, "coupling member limit reached"));
                }
                members.insert(n);
            }
        }
    }
    if members.is_empty() || members.contains(&label) {
        return Err(error(b.line, "coupling surface is empty or includes its reference node"));
    }
    // Native continuum distribution averages translations only. Restrict this subset to
    // solid clouds, where the vendor's default rotational coupling has no nodal rotations.
    let solid_nodes: BTreeSet<_> = model.elements().values().filter(|e| e.kind.dimension() == 3).flat_map(|e| e.nodes.iter().copied()).collect();
    let rotating_nodes: BTreeSet<_> = model.elements().values().filter(|e| e.kind.dimension() != 3).flat_map(|e| e.nodes.iter().copied()).collect();
    if members.iter().any(|n| !solid_nodes.contains(n) || rotating_nodes.contains(n)) {
        return Err(history::unsupported(b.line, "*COUPLING cloud requires nodes belonging only to solid elements"));
    }
    if kind == CouplingKind::Distributing && bits == 7 {
        // The new translation-only subset is the exact equal-weight mean at the
        // cloud centroid. Do not infer a released-rotation rule for an offset RP.
        let mut centroid = [0.0; 3];
        let mut scale = 0.0_f64;
        for node in &members {
            let xyz = model.node(*node).ok_or_else(|| error(b.line, "missing distributing cloud node"))?.xyz;
            for axis in 0..3 {
                centroid[axis] += xyz[axis] / members.len() as f64;
                scale = scale.max(xyz[axis].abs());
            }
        }
        let xyz = model.node(label).ok_or_else(|| error(b.line, "missing distributing reference node"))?.xyz;
        if xyz.iter().zip(centroid).any(|(x, c)| (*x - c).abs() > 64.0 * f64::EPSILON * scale) {
            return Err(history::unsupported(
                b.line,
                "*DISTRIBUTING DOFs 1..3 require REF NODE at the equal-weight cloud centroid; offset reference points are not mapped",
            ));
        }
    }
    let rp = format!("InpReference{label}");
    if model.reference_point(&rp).is_none() {
        let xyz = model.node(label).ok_or_else(|| error(b.line, "missing reference node"))?.xyz;
        model.add_reference_point(&rp, ReferencePoint { xyz, label }).map_err(|e| error(b.line, e.to_string()))?;
    }
    let coupling_set = if element_surface {
        // Keep the actual spelling and face identity of the resolved model set.
        set_name(model, &surf, surface.line)?
    } else {
        let set = Set { nodes: members.into_iter().collect(), ..Set::default() };
        if let Some(existing) = model.set(&surf) {
            if existing != &set {
                return Err(error(surface.line, "coupling surface name collides with another set"));
            }
        } else {
            model.add_set(&surf, set).map_err(|e| error(surface.line, e.to_string()))?;
        }
        surf
    };
    model
        .add_coupling(
            &history::name(b.required("CONSTRAINTNAME")?, b.line)?,
            Coupling { reference_point: rp, set: coupling_set, kind, dofs: Dofs::from_bits(bits) },
        )
        .map_err(|e| error(b.line, e.to_string()))?;
    Ok(vec![b.line, sub.line, surface.line])
}

pub(super) fn import(text: &str, read: &mut Read, preset: Preset, limits: Limits) -> Result<()> {
    let blocks = lex::blocks(text, MAX_RECORDS)?;
    if !blocks.iter().any(|b| matches!(b.keyword.as_str(), "TIE" | "COUPLING")) {
        return Ok(());
    }
    let mut model = read.model.clone();
    let mut consumed = BTreeSet::new();
    let mut in_step = false;
    for (i, b) in blocks.iter().enumerate() {
        if b.keyword == "STEP" {
            in_step = true;
        }
        if in_step {
            continue;
        }
        let result = match b.keyword.as_str() {
            "TIE" => tie(b, &mut model, preset).map(|()| vec![b.line]),
            "COUPLING" => coupling(b, blocks.get(i + 1), &blocks, &mut model, limits),
            _ => continue,
        };
        match result {
            Ok(lines) => consumed.extend(lines),
            Err(Error::Unsupported { what, .. }) => {
                read.warnings.push(what);
                return Ok(());
            }
            Err(e) => return Err(e),
        }
    }
    read.warnings.retain(|w| !consumed.iter().any(|line| w.starts_with(&format!("line {line}: *SURFACE TYPE="))));
    read.model = model;
    let before = read.unknown.len();
    read.unknown.retain(|u| !consumed.contains(&u.line));
    read.unknown_count -= before - read.unknown.len();
    Ok(())
}

// Keep original reference labels available while history expands node sets; native targets are
// materialized separately. Then remove their temporary mesh-node aliases and set memberships.
pub(super) fn finish(model: &mut Model) -> Result<()> {
    let labels: BTreeSet<_> = model.reference_points().values().map(|p| p.label).collect();
    if labels.is_empty() {
        return Ok(());
    }
    let mut sets = model.sets().clone();
    for (name, set) in &mut sets {
        set.nodes.retain(|n| !labels.contains(n));
        if set.is_empty() {
            model.remove_set(name).map_err(|e| error(1, e.to_string()))?;
        } else {
            model.replace_set(name, set.clone()).map_err(|e| error(1, e.to_string()))?;
        }
    }
    for label in labels {
        model.remove_node(label).map_err(|e| error(1, e.to_string()))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn translation_only_distributing_requires_the_explicit_supported_subset() {
        let input = r#"** Original equal-weight translation-only distributing coupling, SI
*NODE
1,0,0,0
2,1,0,0
3,0,1,0
4,0,0,1
99,.3333333333333333,0,.3333333333333333
*ELEMENT,TYPE=C3D4,ELSET=Solid
1,1,2,3,4
*NSET,NSET=Base
1,2,3
*NSET,NSET=Cloud
1,2,4
*NSET,NSET=RP
99
*SURFACE,TYPE=NODE,NAME=Pad
Cloud,1
*MATERIAL,NAME=Elastic
*ELASTIC
120000,.25
*SOLID SECTION,ELSET=Solid,MATERIAL=Elastic
*COUPLING,CONSTRAINT NAME=Mean,REF NODE=RP,SURFACE=Pad
*DISTRIBUTING,COUPLING=CONTINUUM,WEIGHTING METHOD=UNIFORM
1,3
*STEP,NAME=Pull
*STATIC
*BOUNDARY
Base,1,3
RP,4,6
*CLOAD
RP,1,3
*END STEP
"#;
        let read = |text: &str| super::super::read_analysis_limited(text, Preset::MKgS, Limits::default()).unwrap();
        let got = read(input);
        assert!(got.analysis_complete, "{:?}", got.warnings);
        assert_eq!(got.model.couplings().values().next().unwrap().dofs.list(), vec![1, 2, 3]);
        for (old, new, reason) in [
            ("1,3\n*STEP", "1,5\n*STEP", "partial rotational"),
            ("99,.3333333333333333", "99,.4", "centroid"),
            ("WEIGHTING METHOD=UNIFORM", "WEIGHTING METHOD=LINEAR", "WEIGHTING"),
            ("*STATIC\n", "*FREQUENCY\n1\n", "linear STATIC"),
        ] {
            let denied = read(&input.replace(old, new));
            assert!(!denied.analysis_complete);
            assert!(denied.model.steps().is_empty());
            assert!(denied.warnings.iter().any(|w| w.contains(reason)), "{:?}", denied.warnings);
        }
    }
}
