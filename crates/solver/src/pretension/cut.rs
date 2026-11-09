//! Private mesh preparation for a continuous planar pretension section.
//!
//! The duplicated trace uses the same constant-jump elimination as an explicit
//! node-pair cut (Cook et al., 2002, §9.2). No stiffness or load is invented at
//! the face: element virtual work supplies the force conjugate to q.
use crate::{Error, Result};
use fieldcraft_model::{Model, Resolved, Set};
use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};

pub(crate) struct Prepared<'a> {
    pub model: Cow<'a, Model>,
    pub copies: BTreeMap<u32, u32>,
}

/// Called before building Mesh/Setup. The caller's model and its labels never
/// change. The negative-normal side is duplicated, so q = n·(u_minus-u_plus)
/// is positive tightening; source-node output retains the positive-side trace.
pub(crate) fn prepare<'a>(model: &'a Model, resolved: &Resolved, max_dofs: usize) -> Result<Prepared<'a>> {
    if model.pretensions().values().all(|p| p.face.is_none()) {
        return Ok(Prepared { model: Cow::Borrowed(model), copies: BTreeMap::new() });
    }
    let mut out = model.clone();
    // Replacing connectivity invalidates geometric classification. Freeze the
    // resolved regions first, keeping element/face membership (not old traces).
    for (name, set) in &resolved.sets {
        if set.is_empty() {
            continue;
        }
        let mut covered = BTreeSet::new();
        for label in &set.elements {
            if let Some(e) = model.element(*label) {
                covered.extend(e.nodes.iter().copied());
            }
        }
        for face in &set.faces {
            covered.extend(model.face_nodes(*face).unwrap_or_default());
        }
        for edge in &set.edges {
            covered.extend(model.edge_nodes(*edge).unwrap_or_default());
        }
        out.replace_set(
            name,
            Set {
                nodes: set.nodes.iter().filter(|n| !covered.contains(n)).copied().collect(),
                elements: set.elements.clone(),
                faces: set.faces.clone(),
                edges: set.edges.clone(),
                ..Set::default()
            },
        )?;
    }
    let mut next = model.nodes().keys().copied().chain(resolved.reference_point_nodes.values().copied()).max().unwrap_or(0);
    let mut copies = BTreeMap::new();
    let mut cut_nodes = BTreeSet::new();
    for (index, (name, section)) in model.pretensions().iter().enumerate() {
        let Some(face_name) = &section.face else { continue };
        let invalid = |why: &str| Error::Invalid(format!("pretension `{name}` internal face `{face_name}`: {why}"));
        let faces = &resolved.sets.get(face_name).ok_or_else(|| invalid("missing set"))?.faces;
        if faces.is_empty() {
            return Err(invalid("requires an element-face set"));
        }
        let mut targets = BTreeSet::new();
        let mut nodes = BTreeSet::new();
        for face in faces {
            let e = model.element(face.element).ok_or_else(|| invalid("missing element"))?;
            if e.kind.dimension() != 3 {
                return Err(invalid("requires solid elements"));
            }
            let mut labels = model.face_nodes(*face).ok_or_else(|| invalid("invalid element face"))?;
            labels.sort_unstable();
            if !targets.insert(labels.clone()) {
                return Err(invalid("both sides of the same face were selected"));
            }
            nodes.extend(labels);
        }
        if nodes.len() < 3 || nodes.len() > 10_000 {
            return Err(invalid("requires 3..10000 cut nodes"));
        }
        if !cut_nodes.is_disjoint(&nodes) {
            return Err(invalid("intersecting pretension cuts are unsupported"));
        }
        cut_nodes.extend(nodes.iter().copied());
        let xyz = |label| model.node(label).map(|n| n.xyz).ok_or_else(|| invalid("missing node"));
        let origin = xyz(*nodes.first().ok_or_else(|| invalid("empty cut"))?)?;
        let mut extent = 0.0_f64;
        for node in &nodes {
            let x = xyz(*node)?;
            extent = extent.max((0..3).map(|c| (x[c] - origin[c]).powi(2)).sum::<f64>().sqrt());
        }
        if !extent.is_finite() || extent <= 0.0 {
            return Err(invalid("zero cut extent"));
        }
        const PLANAR_TOLERANCE: f64 = 1e-9; // unit-factors: dimensionless geometric tolerance
        let tolerance = PLANAR_TOLERANCE * extent;
        let distance = |x: [f64; 3]| (0..3).map(|c| (x[c] - origin[c]) * section.normal[c]).sum::<f64>();
        for node in &nodes {
            if distance(xyz(*node)?).abs() > tolerance {
                return Err(invalid("cut is not planar normal to its axis"));
            }
        }
        let mut incidence: BTreeMap<Vec<u32>, Vec<bool>> = targets.into_iter().map(|f| (f, Vec::new())).collect();
        let mut negative = Vec::new();
        let mut minus_nodes = BTreeSet::new();
        let mut plus_nodes = BTreeSet::new();
        // Visit every element in the cut-node stars, not just the selected face
        // owners. Edge/corner neighbours on that side must be rewired too.
        for (&label, e) in model.elements() {
            if !e.nodes.iter().any(|n| nodes.contains(n)) {
                continue;
            }
            if e.kind.dimension() != 3 {
                return Err(invalid("cut node is also used by a non-solid element"));
            }
            let mut below = false;
            let mut above = false;
            for node in &e.nodes {
                let d = distance(xyz(*node)?);
                below |= d < -tolerance;
                above |= d > tolerance;
            }
            if below == above {
                return Err(invalid("an incident element crosses or lies in the cut plane"));
            }
            if above && faces.iter().any(|face| face.element == label) {
                return Err(invalid("normal must point away from the selected face's underlying element"));
            }
            if below {
                negative.push(label);
                minus_nodes.extend(e.nodes.iter().copied());
            } else {
                plus_nodes.extend(e.nodes.iter().copied());
            }
            for (i, _) in e.kind.faces().iter().enumerate() {
                let face = fieldcraft_model::ElementFace { element: label, face: (i + 1) as u8 };
                let mut key = model.face_nodes(face).ok_or_else(|| invalid("invalid adjacent face"))?;
                key.sort_unstable();
                if let Some(sides) = incidence.get_mut(&key) {
                    sides.push(below);
                }
            }
        }
        for sides in incidence.values() {
            if sides.len() != 2 || sides.first() == sides.get(1) {
                return Err(invalid("each face must have exactly two connected solid neighbours on opposite sides"));
            }
        }
        if minus_nodes.intersection(&plus_nodes).any(|n| !nodes.contains(n)) {
            return Err(invalid("incomplete cut: the two sides share unselected nodes"));
        }
        if model.nodes().len().saturating_add(copies.len()).saturating_add(nodes.len()) > max_dofs / 3 {
            return Err(Error::Invalid(format!("pretension `{name}` internal copies exceed the degree-of-freedom limit")));
        }
        let mut pairs = BTreeMap::new();
        for node in &nodes {
            next = next.checked_add(1).filter(|n| *n <= fieldcraft_model::MAX_LABEL).ok_or_else(|| invalid("no labels left for internal copies"))?;
            out.set_node(next, xyz(*node)?)?;
            pairs.insert(*node, next);
            copies.insert(next, *node);
        }
        for label in negative {
            let e = out.element(label).ok_or_else(|| invalid("missing internal element"))?;
            let kind = e.kind;
            let connectivity = e.nodes.iter().map(|n| pairs.get(n).copied().unwrap_or(*n)).collect();
            out.set_element(label, kind, connectivity)?;
        }
        let mut private = section.clone();
        let names = [format!("InternalPretension{index}First"), format!("InternalPretension{index}Second")];
        if names.iter().any(|n| out.set(n).is_some()) {
            return Err(invalid("reserved internal cut-set name is already used"));
        }
        out.add_set(&names[0], Set { nodes: pairs.values().copied().collect(), ..Set::default() })?;
        out.add_set(&names[1], Set { nodes: pairs.keys().copied().collect(), ..Set::default() })?;
        private.first = names[0].clone();
        private.second = names[1].clone();
        private.face = None;
        out.replace_pretension(name, private)?;
    }
    Ok(Prepared { model: Cow::Owned(out), copies })
}

fn retain<T>(values: &mut Vec<T>, keep: &[bool]) {
    let mut i = 0;
    values.retain(|_| {
        let yes = keep.get(i).copied().unwrap_or(false);
        i += 1;
        yes
    });
}

impl Prepared<'_> {
    pub fn project_frequency(&self, result: &mut crate::FrequencyResult) {
        let keep: Vec<bool> = result.nodes.iter().map(|n| !self.copies.contains_key(n)).collect();
        for mode in &mut result.modes {
            retain(&mut mode.shape, &keep);
            retain(&mut mode.rotation, &keep);
        }
        retain(&mut result.nodes, &keep);
    }

    pub fn project_harmonic(&self, result: &mut crate::HarmonicResult) {
        let keep: Vec<bool> = result.nodes.iter().map(|n| !self.copies.contains_key(n)).collect();
        for frame in &mut result.frames {
            retain(&mut frame.real, &keep);
            retain(&mut frame.imaginary, &keep);
            retain(&mut frame.rotation_real, &keep);
            retain(&mut frame.rotation_imaginary, &keep);
            if let Some(stress) = &mut frame.stress_real {
                retain(stress, &keep);
            }
            if let Some(stress) = &mut frame.stress_imaginary {
                retain(stress, &keep);
            }
        }
        retain(&mut result.nodes, &keep);
    }

    /// Strip only private trace nodes at the public result boundary. Element
    /// stresses and the scalar jump remain available; no averaged displacement
    /// is invented for a discontinuity.
    pub fn project(&self, result: &mut crate::StaticResult) {
        if self.copies.is_empty() {
            return;
        }
        let keep: Vec<bool> = result.nodes.iter().map(|n| !self.copies.contains_key(n)).collect();
        let owners: Vec<usize> = result
            .nodes
            .iter()
            .enumerate()
            .map(|(i, n)| self.copies.get(n).and_then(|original| result.nodes.binary_search(original).ok()).unwrap_or(i))
            .collect();
        // Reactions from either trace act on the same source node. Unlike
        // displacement/stress traces, they must be summed, never discarded.
        for (i, owner) in owners.iter().copied().enumerate() {
            if i == owner {
                continue;
            }
            for values in [&mut result.reaction, &mut result.reaction_moment] {
                if let Some(force) = values.get(i).copied()
                    && let Some(total) = values.get_mut(owner)
                {
                    for c in 0..3 {
                        total[c] += force[c];
                    }
                }
            }
        }
        retain(&mut result.nodes, &keep);
        retain(&mut result.displacement, &keep);
        retain(&mut result.rotation, &keep);
        retain(&mut result.reaction, &keep);
        retain(&mut result.reaction_moment, &keep);
        retain(&mut result.stress, &keep);
        retain(&mut result.von_mises, &keep);
        if let Some(t) = &mut result.temperature {
            retain(t, &keep);
        }
        if let Some(history) = &mut result.nonlinear {
            let mut next = 0;
            let indices: Vec<Option<usize>> = keep
                .iter()
                .map(|keep| {
                    if *keep {
                        let i = next;
                        next += 1;
                        Some(i)
                    } else {
                        None
                    }
                })
                .collect();
            for frame in &mut history.frames {
                retain(&mut frame.displacement, &keep);
                retain(&mut frame.stress, &keep);
                let mut reactions = BTreeMap::<usize, [f64; 3]>::new();
                for (i, f) in &frame.reaction {
                    if let Some(j) = owners.get(*i).and_then(|owner| indices.get(*owner)).copied().flatten() {
                        let total = reactions.entry(j).or_default();
                        for c in 0..3 {
                            total[c] += f[c];
                        }
                    }
                }
                frame.reaction = reactions.into_iter().collect();
            }
        }
    }
}
