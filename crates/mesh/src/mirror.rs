//! Mirroring a meshed region across a plane.
//!
//! Each element of the region gets a copy reflected across the plane, `x' = x − 2 ((x − p)·n) n`.
//! A reflection turns a right-handed element left-handed (its Jacobian determinant changes sign),
//! so the copy lists its nodes in the reflected order of [`reflected_order`]: corners swapped so
//! the element is right-handed again, mid-side nodes following their edges. Nodes the region has
//! on the plane are their own mirror images; the copies there are merged with them
//! ([`crate::join::merge_coincident`]), so the region and its image share the plane's nodes and
//! form one conforming mesh.

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_model::{ElementKind, MAX_LABEL, Model};

use crate::{Error, Result};

/// The reflected node order of an element kind: the copy of an element with nodes `n` lists
/// `n[order[i]]` (mirrored) as its node `i`. `None` for kinds that are not mirrored yet.
pub fn reflected_order(kind: ElementKind) -> Option<&'static [usize]> {
    Some(match kind {
        // Corners 1 and 2 swap; the mid-edge nodes follow: 0-1 ↔ 0-2 (4 ↔ 6), 1-3 ↔ 2-3 (8 ↔ 9).
        ElementKind::Tet4 => &[0, 2, 1, 3],
        ElementKind::Tet10 => &[0, 2, 1, 3, 6, 5, 4, 7, 9, 8],
        // Each triangle's last two corners swap.
        ElementKind::Wedge6 => &[0, 2, 1, 3, 5, 4],
        // The base's corners in reverse cyclic order; the apex stays.
        ElementKind::Pyramid5 => &[0, 3, 2, 1, 4],
        // Each quadrilateral's corners in reverse cyclic order from its first.
        ElementKind::Hex8 => &[0, 3, 2, 1, 4, 7, 6, 5],
        ElementKind::Hex20 => &[0, 3, 2, 1, 4, 7, 6, 5, 11, 10, 9, 8, 15, 14, 13, 12, 16, 19, 18, 17],
        _ => return None,
    })
}

/// For each face (1-based) of an element, the face of its reflected copy over the same nodes.
pub fn reflected_faces(kind: ElementKind) -> Option<Vec<u8>> {
    let order = reflected_order(kind)?;
    let corners = kind.corner_count();
    // Original node position → the copy's position of its image.
    let mut inverse = vec![0usize; order.len()];
    for (new, &old) in order.iter().enumerate() {
        *inverse.get_mut(old)? = new;
    }
    let key = |face: &[usize]| {
        let mut k: Vec<usize> = face.iter().copied().filter(|&i| i < corners).collect();
        k.sort_unstable();
        k
    };
    kind.faces()
        .iter()
        .map(|face| {
            let image: Vec<usize> = face.iter().filter_map(|&i| inverse.get(i).copied()).collect();
            let want = key(&image);
            let k = kind.faces().iter().position(|f| key(f) == want)?;
            u8::try_from(k + 1).ok()
        })
        .collect()
}

/// The same for edges.
pub fn reflected_edges(kind: ElementKind) -> Option<Vec<u8>> {
    let order = reflected_order(kind)?;
    let mut inverse = vec![0usize; order.len()];
    for (new, &old) in order.iter().enumerate() {
        *inverse.get_mut(old)? = new;
    }
    let key = |e: &[usize]| {
        let (a, b) = (*e.first().unwrap_or(&0), *e.get(1).unwrap_or(&0));
        (a.min(b), a.max(b))
    };
    kind.edges()
        .iter()
        .map(|edge| {
            let image: Vec<usize> = edge.iter().take(2).filter_map(|&i| inverse.get(i).copied()).collect();
            let k = kind.edges().iter().position(|e| key(e) == key(&image))?;
            u8::try_from(k + 1).ok()
        })
        .collect()
}

/// What [`mirror`] added.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Mirrored {
    /// Each mirrored element and the label of its copy.
    pub elements: BTreeMap<u32, u32>,
    /// Each node of the region and the label of its image (itself for a node on the plane).
    pub nodes: BTreeMap<u32, u32>,
    /// Nodes on the plane, shared by the region and its image.
    pub on_plane: usize,
}

/// Reflect the elements `region` across the plane through `point` with normal `normal` (any
/// length but zero), adding the copies after the model's largest labels, and merge the copies of
/// the nodes within `tol` of the plane with the originals. Geometry tags are not copied (the
/// caller names the image). Refused, with the model unchanged, for an element kind that is not
/// mirrored yet ([`reflected_order`]) or a node too far for the coordinate range.
pub fn mirror(model: &mut Model, region: &BTreeSet<u32>, point: [f64; 3], normal: [f64; 3], tol: f64) -> Result<Mirrored> {
    let len = normal.iter().map(|x| x * x).sum::<f64>().sqrt();
    if !(len > 1e-12 && len.is_finite() && point.iter().all(|x| x.is_finite())) {
        return Err(Error::Invalid("the mirror plane needs a finite point and a non-zero normal".into()));
    }
    if !(tol.is_finite() && tol > 0.0) {
        return Err(Error::invalid(crate::message!["the merge distance must be a positive number, not ", crate::message::length(tol)]));
    }
    let n = normal.map(|x| x / len);
    let mut elements = Vec::with_capacity(region.len());
    for &e in region {
        let el = model.element(e).ok_or_else(|| Error::Invalid(format!("no element {e}")))?;
        let order = reflected_order(el.kind)
            .ok_or_else(|| Error::Invalid(format!("not supported yet: mirroring {} elements (element {e})", el.kind.name())))?;
        elements.push((e, el.kind, order, el.nodes.clone()));
    }
    let nodes: BTreeSet<u32> = elements.iter().flat_map(|(.., ns)| ns.iter().copied()).collect();
    let first_node = model.nodes().keys().next_back().map_or(Some(1), |l| l.checked_add(1));
    let first_element = model.elements().keys().next_back().map_or(Some(1), |l| l.checked_add(1));
    let (Some(first_node), Some(first_element)) = (first_node, first_element) else {
        return Err(Error::Invalid(format!("the mirror needs labels beyond {MAX_LABEL}")));
    };
    let fits = |first: u32, count: usize| u64::from(first) + count as u64 <= u64::from(MAX_LABEL) + 1;
    if !fits(first_node, nodes.len()) || !fits(first_element, elements.len()) {
        return Err(Error::Invalid(format!("the mirror needs labels beyond {MAX_LABEL}")));
    }
    let mut work = model.clone();
    let mut image = BTreeMap::new();
    let mut plane_nodes = BTreeSet::new();
    for (k, &node) in nodes.iter().enumerate() {
        let x = work.node(node).map(|n| n.xyz).ok_or_else(|| Error::Invalid(format!("no node {node}")))?;
        let d = (0..3).map(|i| (x[i] - point[i]) * n[i]).sum::<f64>();
        let label = first_node + k as u32;
        work.set_node(label, [x[0] - 2.0 * d * n[0], x[1] - 2.0 * d * n[1], x[2] - 2.0 * d * n[2]]).map_err(|e| Error::Invalid(e.to_string()))?;
        if d.abs() <= tol {
            plane_nodes.insert(node);
            plane_nodes.insert(label);
        }
        image.insert(node, label);
    }
    let mut copies = BTreeMap::new();
    for (k, (e, kind, order, ns)) in elements.iter().enumerate() {
        let label = first_element + k as u32;
        let copy: Vec<u32> = order.iter().map(|&i| ns.get(i).and_then(|n| image.get(n)).copied().unwrap_or(0)).collect();
        work.set_element(label, *kind, copy).map_err(|e| Error::Invalid(e.to_string()))?;
        copies.insert(*e, label);
    }
    // The images of the nodes on the plane are the nodes themselves.
    let merged = crate::join::merge_coincident(&mut work, Some(&plane_nodes), tol)?;
    for to in image.values_mut() {
        if let Some(&survivor) = merged.get(to) {
            *to = survivor;
        }
    }
    let on_plane = merged.len();
    *model = work;
    Ok(Mirrored { elements: copies, nodes: image, on_plane })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_model::ElementFace;

    /// The determinant of an element's Jacobian at its first corner, from its three corner edges
    /// there (the sign of its orientation).
    fn orientation(m: &Model, e: u32) -> f64 {
        let el = m.element(e).unwrap();
        let p = |i: usize| m.node(el.nodes[i]).unwrap().xyz;
        let (a, b, c) = match el.kind {
            ElementKind::Tet4 | ElementKind::Tet10 => (1, 2, 3),
            ElementKind::Wedge6 => (1, 2, 3),
            _ => (1, 3, 4),
        };
        let o = p(0);
        let d = |i: usize| [p(i)[0] - o[0], p(i)[1] - o[1], p(i)[2] - o[2]];
        let (u, v, w) = (d(a), d(b), d(c));
        u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0])
    }

    /// One element of `kind` with its corners at the unit reference shape, offset by 1 in x,
    /// mid-side nodes at the edge middles.
    fn one(kind: ElementKind) -> Model {
        let corners: &[[f64; 3]] = match kind {
            ElementKind::Tet4 | ElementKind::Tet10 => &[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            ElementKind::Wedge6 => &[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [0., 1., 1.]],
            _ => &[[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]],
        };
        let mut pts: Vec<[f64; 3]> = corners.to_vec();
        if kind.is_quadratic() {
            let mut mids = vec![[0.0; 3]; kind.node_count() - kind.corner_count()];
            for e in kind.edges() {
                let (a, b, m) = (corners[e[0]], corners[e[1]], e[2]);
                mids[m - kind.corner_count()] = [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, (a[2] + b[2]) / 2.0];
            }
            pts.extend(mids);
        }
        let mut m = Model::default();
        for (i, p) in pts.iter().enumerate() {
            m.set_node(i as u32 + 1, [p[0] + 1.0, p[1], p[2]]).unwrap();
        }
        m.set_element(1, kind, (1..=pts.len() as u32).collect()).unwrap();
        m
    }

    #[test]
    fn every_kind_mirrors_right_handed_with_faces_and_mid_side_nodes_in_place() {
        for kind in [ElementKind::Tet4, ElementKind::Tet10, ElementKind::Wedge6, ElementKind::Hex8, ElementKind::Hex20] {
            let mut m = one(kind);
            assert!(orientation(&m, 1) > 0.0, "{kind:?}");
            // Across x = 1: the element's face at x = 1 lies on the plane.
            let r = mirror(&mut m, &[1].into(), [1.0, 0.0, 0.0], [-2.0, 0.0, 0.0], 1e-9).unwrap();
            let copy = r.elements[&1];
            assert!(orientation(&m, copy) > 0.0, "{kind:?} copy is left-handed");
            // Every mid-side node of the copy is the middle of its edge.
            let el = m.element(copy).unwrap().clone();
            for e in kind.edges().iter().filter(|e| e.len() == 3) {
                let p = |i: usize| m.node(el.nodes[e[i]]).unwrap().xyz;
                for k in 0..3 {
                    assert!((p(2)[k] - (p(0)[k] + p(1)[k]) / 2.0).abs() < 1e-12, "{kind:?} edge {e:?}");
                }
            }
            // Each face maps to the face of the copy over the images of its nodes.
            let faces = reflected_faces(kind).unwrap();
            for (f, &g) in faces.iter().enumerate() {
                let orig = m.face_nodes(ElementFace { element: 1, face: f as u8 + 1 }).unwrap();
                let mut want: Vec<u32> = orig.iter().map(|n| r.nodes[n]).collect();
                let mut got = m.face_nodes(ElementFace { element: copy, face: g }).unwrap();
                want.sort_unstable();
                got.sort_unstable();
                assert_eq!(got, want, "{kind:?} face {}", f + 1);
            }
            assert_eq!(reflected_edges(kind).unwrap().len(), kind.edges().len());
            // The nodes on the plane are shared.
            // The face at x = 1 (with its mid-side nodes).
            let shared = match kind {
                ElementKind::Tet4 => 3,
                ElementKind::Tet10 => 6,
                ElementKind::Hex20 => 8,
                _ => 4,
            };
            assert_eq!(r.on_plane, shared, "{kind:?}");
            assert_eq!(m.nodes().len(), 2 * kind.node_count() - shared, "{kind:?}");
        }
    }

    #[test]
    fn unsupported_kinds_and_hostile_planes_change_nothing() {
        let mut m = one(ElementKind::Hex8);
        let before = m.clone();
        for (p, n) in [([0.0; 3], [0.0; 3]), ([f64::NAN, 0.0, 0.0], [1.0, 0.0, 0.0]), ([0.0; 3], [f64::INFINITY, 0.0, 0.0])] {
            assert!(mirror(&mut m, &[1].into(), p, n, 1e-9).is_err());
        }
        assert!(mirror(&mut m, &[1].into(), [0.0; 3], [1.0, 0.0, 0.0], 0.0).is_err());
        assert!(mirror(&mut m, &[2].into(), [0.0; 3], [1.0, 0.0, 0.0], 1e-9).is_err());
        assert_eq!(m, before);
        let mut q = Model::default();
        for (l, p) in [(1, [0., 0., 0.]), (2, [1., 0., 0.]), (3, [0., 1., 0.])] {
            q.set_node(l, p).unwrap();
        }
        q.set_element(1, ElementKind::Tri3, vec![1, 2, 3]).unwrap();
        let e = mirror(&mut q, &[1].into(), [0.0; 3], [1.0, 0.0, 0.0], 1e-9).unwrap_err();
        assert!(e.to_string().contains("not supported yet"), "{e}");
    }
}
