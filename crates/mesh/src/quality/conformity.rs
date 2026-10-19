//! Whether a solid mesh is conforming: the checks that find cracks.
//!
//! A mesh of solids without cracks has
//! - free faces (element faces used by one element) that form a closed 2-manifold: every corner
//!   edge of a free face lies on exactly two free faces ([`Conformity::open_edges`]),
//! - no two distinct nodes at the same place ([`Conformity::coincident`]): a duplicate node
//!   splits the elements around it, and the faces between them become internal free faces,
//! - one mid-side node per element edge for quadratic kinds ([`Conformity::midside_mismatches`]):
//!   two elements on the same corner edge with different mid-side nodes are joined at the
//!   corners and open along the edge.
//!
//! The report is deterministic: every list is sorted.

use fieldcraft_model::ElementKind;

use crate::{Error, Result};

/// What [`conformity`] found.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Conformity {
    /// Free faces, as (element, face number from 1 as [`ElementKind::faces`] numbers them).
    pub free_faces: Vec<(u32, u8)>,
    /// Corner edges (sorted node pair) of free faces that are not on exactly two free faces.
    pub open_edges: Vec<[u32; 2]>,
    /// Pairs of distinct nodes (lower first) closer than the tolerance.
    pub coincident: Vec<[u32; 2]>,
    /// Corner edges (sorted node pair) that elements give different mid-side nodes.
    pub midside_mismatches: Vec<[u32; 2]>,
}

impl Conformity {
    /// No open edges, no coincident nodes and no mismatched mid-side nodes.
    pub fn is_conforming(&self) -> bool {
        self.open_edges.is_empty() && self.coincident.is_empty() && self.midside_mismatches.is_empty()
    }
}

/// Check the solid elements (`(label, kind, nodes)`; shells and line elements are skipped) over
/// the nodes `(label, xyz)`: nodes closer than `tolerance` count as coincident. Elements that
/// name a node missing from `nodes` are an error, as are non-finite coordinates or tolerance.
pub fn conformity<'a>(
    nodes: &[(u32, [f64; 3])],
    elements: impl IntoIterator<Item = (u32, ElementKind, &'a [u32])>,
    tolerance: f64,
) -> Result<Conformity> {
    if !(tolerance >= 0.0 && tolerance.is_finite()) {
        return Err(Error::invalid(crate::message!["the tolerance must be a non-negative number, got ", crate::message::length(tolerance)]));
    }
    let mut labels: Vec<u32> = nodes.iter().map(|n| n.0).collect();
    labels.sort_unstable();
    // (sorted corners, corners in order, element, face), triangles padded with `u32::MAX`.
    let mut faces: Vec<([u32; 4], [u32; 4], u32, u8)> = Vec::new();
    // (corner pair, mid-side node) of every edge of a quadratic element.
    let mut mids: Vec<([u32; 2], u32)> = Vec::new();
    for (label, kind, conn) in elements {
        if kind.dimension() != 3 {
            continue;
        }
        if conn.len() != kind.node_count() {
            return Err(Error::Invalid(format!("element {label} has {} nodes; a {} has {}", conn.len(), kind.name(), kind.node_count())));
        }
        if let Some(n) = conn.iter().find(|n| labels.binary_search(n).is_err()) {
            return Err(Error::Invalid(format!("element {label} names node {n}, which does not exist")));
        }
        let node = |i: usize| conn.get(i).copied().unwrap_or(u32::MAX);
        for (k, f) in kind.faces().iter().enumerate() {
            let mut ordered = [u32::MAX; 4];
            for (slot, &i) in ordered.iter_mut().zip(f.iter().take(face_corners(f.len()))) {
                *slot = node(i);
            }
            let mut key = ordered;
            key.sort_unstable();
            faces.push((key, ordered, label, u8::try_from(k + 1).unwrap_or(u8::MAX)));
        }
        if kind.is_quadratic() {
            for e in kind.edges() {
                if let (Some(&a), Some(&b), Some(&m)) = (e.first(), e.get(1), e.get(2)) {
                    let (a, b) = (node(a), node(b));
                    mids.push(([a.min(b), a.max(b)], node(m)));
                }
            }
        }
    }
    faces.sort_unstable();
    let mut out = Conformity::default();
    let mut rim: Vec<[u32; 2]> = Vec::new();
    for run in faces.chunk_by(|a, b| a.0 == b.0) {
        let [(_, ordered, element, face)] = run else { continue };
        out.free_faces.push((*element, *face));
        let ordered: Vec<u32> = ordered.iter().copied().filter(|&n| n != u32::MAX).collect();
        for (j, &a) in ordered.iter().enumerate() {
            let b = ordered.get((j + 1) % ordered.len()).copied().unwrap_or(a);
            rim.push([a.min(b), a.max(b)]);
        }
    }
    out.free_faces.sort_unstable();
    rim.sort_unstable();
    out.open_edges = rim.chunk_by(|a, b| a == b).filter(|r| r.len() != 2).filter_map(|r| r.first().copied()).collect();
    mids.sort_unstable();
    mids.dedup();
    out.midside_mismatches = mids.chunk_by(|a, b| a.0 == b.0).filter(|r| r.len() > 1).filter_map(|r| r.first().map(|m| m.0)).collect();
    out.coincident = coincident(nodes, tolerance)?;
    Ok(out)
}

/// Corners of a face with `n` node positions: a quadratic face lists its mid-side nodes after
/// its corners, as many of each.
fn face_corners(n: usize) -> usize {
    match n {
        6 | 8 => n / 2,
        n => n,
    }
}

/// Pairs of distinct nodes closer than `tol`: nodes are bucketed in a grid of cells `tol` wide
/// (sorted, so the work stays O(n log n) on any input) and compared with the 27 cells around.
fn coincident(nodes: &[(u32, [f64; 3])], tol: f64) -> Result<Vec<[u32; 2]>> {
    if let Some((l, _)) = nodes.iter().find(|(_, p)| !p.iter().all(|v| v.is_finite())) {
        return Err(Error::Invalid(format!("node {l} has a non-finite coordinate")));
    }
    // A zero tolerance still needs a cell size; exact copies share a cell at any size.
    let h = if tol > 0.0 { tol } else { 1.0 };
    let cell = |p: &[f64; 3]| p.map(|v| (v / h).floor() as i64);
    let mut grid: Vec<([i64; 3], u32, [f64; 3])> = nodes.iter().map(|(l, p)| (cell(p), *l, *p)).collect();
    grid.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut out = Vec::new();
    for (c, la, pa) in &grid {
        for d in 0..27i64 {
            let near = [c[0].saturating_add(d % 3 - 1), c[1].saturating_add((d / 3) % 3 - 1), c[2].saturating_add(d / 9 - 1)];
            let start = grid.partition_point(|g| g.0 < near);
            for (_, lb, pb) in grid.get(start..).unwrap_or(&[]).iter().take_while(|g| g.0 == near) {
                let d2: f64 = (0..3).map(|k| (pa[k] - pb[k]).powi(2)).sum();
                if la < lb && d2 <= tol * tol {
                    out.push([*la, *lb]);
                }
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two unit hex8 cells side by side along x (nodes 1–12), elements 1 and 2.
    fn pair() -> (Vec<(u32, [f64; 3])>, Vec<Vec<u32>>) {
        let mut nodes = Vec::new();
        for k in 0..2 {
            for j in 0..2 {
                for i in 0..3 {
                    nodes.push((1 + i + 3 * j + 6 * k, [f64::from(i), f64::from(j), f64::from(k)]));
                }
            }
        }
        let hex = |i: u32| vec![1 + i, 2 + i, 5 + i, 4 + i, 7 + i, 8 + i, 11 + i, 10 + i];
        (nodes, vec![hex(0), hex(1)])
    }

    fn check(nodes: &[(u32, [f64; 3])], elements: &[Vec<u32>], kind: ElementKind) -> Conformity {
        conformity(nodes, elements.iter().enumerate().map(|(e, n)| (e as u32 + 1, kind, n.as_slice())), 1e-9).unwrap()
    }

    #[test]
    fn a_conforming_pair_has_a_closed_skin() {
        let (nodes, elements) = pair();
        let c = check(&nodes, &elements, ElementKind::Hex8);
        assert!(c.is_conforming(), "{c:?}");
        assert_eq!(c.free_faces.len(), 10);
        assert!(!c.free_faces.contains(&(1, 4)) && !c.free_faces.contains(&(2, 6)), "the shared face is inside");
    }

    #[test]
    fn a_duplicated_node_is_a_crack() {
        let (mut nodes, mut elements) = pair();
        // Element 2 takes its own copy of node 2 (at x = 1).
        nodes.push((99, [1.0, 0.0, 0.0]));
        elements[1][0] = 99;
        let c = check(&nodes, &elements, ElementKind::Hex8);
        assert_eq!(c.coincident, vec![[2, 99]]);
        assert!(!c.open_edges.is_empty(), "the faces beside the copy no longer close up");
        assert_eq!(c.free_faces.len(), 12, "the faces between the elements are free");
    }

    #[test]
    fn different_mid_side_nodes_on_one_edge_are_a_crack() {
        let (mut nodes, elements) = pair();
        let mut next = 100;
        let mut quadratic = Vec::new();
        for (e, n) in elements.iter().enumerate() {
            let mut q = n.clone();
            for edge in ElementKind::Hex20.edges() {
                next += 1;
                let (a, b) = (n[edge[0]], n[edge[1]]);
                let (pa, pb) = (nodes.iter().find(|x| x.0 == a).unwrap().1, nodes.iter().find(|x| x.0 == b).unwrap().1);
                nodes.push((next, [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2]) + e as f64 * 1e-3]));
                q.push(next);
            }
            quadratic.push(q);
        }
        let c = check(&nodes, &quadratic, ElementKind::Hex20);
        // The four edges of the shared face got a mid-side node from each element.
        assert_eq!(c.midside_mismatches.len(), 4);
        assert!(c.open_edges.is_empty(), "the corners still close up");
    }

    #[test]
    fn hostile_input_is_an_error() {
        let (nodes, elements) = pair();
        assert!(conformity(&nodes, std::iter::once((1, ElementKind::Hex8, &elements[0][..4])), 0.0).is_err());
        assert!(conformity(&nodes, std::iter::once((1, ElementKind::Hex8, &[1, 2, 3, 4, 5, 6, 7, 1000][..])), 0.0).is_err());
        assert!(conformity(&nodes, std::iter::empty(), f64::NAN).is_err());
        assert!(conformity(&[(1, [f64::INFINITY, 0.0, 0.0])], std::iter::empty(), 0.0).is_err());
        let far = [(1, [1e300, 0.0, 0.0]), (2, [1e300, 0.0, 0.0]), (3, [-1e300, 0.0, 0.0])];
        assert_eq!(conformity(&far, std::iter::empty(), 1e-300).unwrap().coincident, vec![[1, 2]]);
    }
}
