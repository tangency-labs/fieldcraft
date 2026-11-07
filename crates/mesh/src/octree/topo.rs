//! Hexahedral mesh topology for the boundary stages: boundary faces, the manifold clean-up of
//! the core, padding a boundary with a layer of hexahedra, and pillowing a shrink set.
//!
//! *Padding* puts one hexahedron on every boundary quad, between the quad and a copy of it, so
//! no hexahedron keeps more than one face on the boundary. *Pillowing* a shrink set `S` inserts
//! a layer of hexahedra along the faces between `S` and the rest of the mesh; the layer reaches
//! the boundary along the loop where `S`'s boundary faces meet the others' (S. A. Mitchell,
//! T. J. Tautges, *Pillowing doublets: refining a mesh to ensure that faces share at most one
//! edge*, IMR 1995; Gao, Shen, Panozzo 2019, §3.4 use the same insertion around features).

use std::collections::{BTreeMap, HashMap};

use crate::{Error, Result};

/// The six faces of a hex8 as local node positions, right-hand normals outward (the model's
/// `ElementKind::Hex8.faces()` order).
pub(crate) const FACES: [[usize; 4]; 6] = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]];

/// A face key: the vertices sorted.
pub(crate) fn key(q: [u32; 4]) -> [u32; 4] {
    let mut k = q;
    k.sort_unstable();
    k
}

/// The root of `i` in the union-find forest `c`, compressing the path from `i`.
pub(crate) fn find(c: &mut [usize], i: usize) -> usize {
    let mut r = i;
    while let Some(&p) = c.get(r).filter(|&&p| p != r) {
        r = p;
    }
    if let Some(x) = c.get_mut(i) {
        *x = r;
    }
    r
}

/// A face of an element: (hexahedron, local face 0..6).
pub(crate) type FaceRef = (u32, u8);

/// The vertices of a local face, outward.
pub(crate) fn face(h: &[u32; 8], f: u8) -> [u32; 4] {
    FACES.get(f as usize).map_or([0; 4], |l| l.map(|i| h[i]))
}

/// Every face's elements, keyed by sorted vertices (deterministic map).
pub(crate) fn face_map(hexes: &[[u32; 8]]) -> BTreeMap<[u32; 4], Vec<FaceRef>> {
    let mut m: BTreeMap<[u32; 4], Vec<FaceRef>> = BTreeMap::new();
    for (i, h) in hexes.iter().enumerate() {
        for f in 0..6u8 {
            m.entry(key(face(h, f))).or_default().push((i as u32, f));
        }
    }
    m
}

/// The boundary faces (on one element only), outward, in face-map order.
pub(crate) fn boundary(hexes: &[[u32; 8]]) -> Vec<(FaceRef, [u32; 4])> {
    face_map(hexes)
        .into_values()
        .filter(|v| v.len() == 1)
        .filter_map(|v| v.first().copied())
        .map(|(h, f)| ((h, f), hexes.get(h as usize).map_or([0; 4], |x| face(x, f))))
        .collect()
}

/// Edges of a quad list with the quads using each (ascending vertex pair → quad indices).
pub(crate) fn quad_edges(quads: &[[u32; 4]]) -> BTreeMap<(u32, u32), Vec<u32>> {
    let mut m: BTreeMap<(u32, u32), Vec<u32>> = BTreeMap::new();
    for (i, q) in quads.iter().enumerate() {
        for k in 0..4 {
            let (a, b) = (q[k], q[(k + 1) % 4]);
            m.entry((a.min(b), a.max(b))).or_default().push(i as u32);
        }
    }
    m
}

/// The hexahedra to drop so that the boundary of `keep` is a 2-manifold: elements at edges with
/// four boundary faces (elements touching along an edge only) and at vertices whose boundary
/// faces form more than one fan. Returns the indices to drop (empty when manifold).
pub(crate) fn nonmanifold(hexes: &[[u32; 8]], keep: &[bool]) -> Vec<u32> {
    let kept: Vec<u32> = (0..hexes.len() as u32).filter(|&i| keep.get(i as usize).copied().unwrap_or(false)).collect();
    let sub: Vec<[u32; 8]> = kept.iter().filter_map(|&i| hexes.get(i as usize).copied()).collect();
    let b = boundary(&sub);
    let quads: Vec<[u32; 4]> = b.iter().map(|(_, q)| *q).collect();
    let mut bad_vertices: Vec<u32> = Vec::new();
    for ((a, c), qs) in quad_edges(&quads) {
        if qs.len() != 2 {
            bad_vertices.extend([a, c]);
        }
    }
    // Fans: around each vertex, boundary quads connected through shared edges at the vertex.
    let mut at: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for (i, q) in quads.iter().enumerate() {
        for &v in q {
            at.entry(v).or_default().push(i as u32);
        }
    }
    for (v, qs) in &at {
        let mut comp: Vec<usize> = (0..qs.len()).collect();
        for i in 0..qs.len() {
            for j in i + 1..qs.len() {
                let (qa, qb) = (quads[qs[i] as usize], quads[qs[j] as usize]);
                // Share an edge through v: another common vertex adjacent to v in both.
                let nb = |q: [u32; 4]| {
                    let p = q.iter().position(|&x| x == *v).unwrap_or(0);
                    [q[(p + 1) % 4], q[(p + 3) % 4]]
                };
                if nb(qa).iter().any(|x| nb(qb).contains(x)) {
                    let (ri, rj) = (find(&mut comp, i), find(&mut comp, j));
                    comp[ri] = rj;
                }
            }
        }
        let roots = (0..qs.len()).filter(|&i| find(&mut comp, i) == i).count();
        if roots > 1 {
            bad_vertices.push(*v);
        }
    }
    bad_vertices.sort_unstable();
    bad_vertices.dedup();
    if bad_vertices.is_empty() {
        return Vec::new();
    }
    kept.into_iter().filter(|&i| hexes.get(i as usize).is_some_and(|h| h.iter().any(|v| bad_vertices.binary_search(v).is_ok()))).collect()
}

/// Pad: a hexahedron on each listed boundary face (outward), between it and a copy; the copies
/// of the faces' vertices are new points (appended, initially at the original positions).
/// Returns the copy of each original vertex.
pub(crate) fn pad(points: &mut Vec<[f64; 3]>, hexes: &mut Vec<[u32; 8]>, faces: &[[u32; 4]]) -> Result<BTreeMap<u32, u32>> {
    let mut copy: BTreeMap<u32, u32> = BTreeMap::new();
    for q in faces {
        for &v in q {
            if let std::collections::btree_map::Entry::Vacant(slot) = copy.entry(v) {
                let n = u32::try_from(points.len()).map_err(|_| Error::Limit("too many points".into()))?;
                points.push(points.get(v as usize).copied().unwrap_or([0.0; 3]));
                slot.insert(n);
            }
        }
    }
    for q in faces {
        let c = q.map(|v| copy.get(&v).copied().unwrap_or(v));
        hexes.push([q[0], q[1], q[2], q[3], c[0], c[1], c[2], c[3]]);
    }
    Ok(copy)
}

/// One vertex copy made by [`pillow`]: the original, the copy, and the shrink-set elements
/// around the original that now use the copy.
#[derive(Debug, Clone)]
pub(crate) struct PillowCopy {
    pub original: u32,
    pub copy: u32,
    pub elements: Vec<u32>,
}

/// Pillow the shrink set `inside` (one flag per element): a new layer of hexahedra on every face
/// between an element of the set and one outside it. Vertices of those faces get one copy per
/// group of set elements around them (connected through faces inside the set), which the set's
/// elements now use; the new layer joins the originals to the copies. Faces on the boundary are
/// not pillowed: where the layer meets the boundary, the set's boundary faces move to the copies,
/// so the boundary gains a ring of quads along the loop. New points are appended at the
/// originals' positions. Returns the copies made.
pub(crate) fn pillow(points: &mut Vec<[f64; 3]>, hexes: &mut Vec<[u32; 8]>, inside: &[bool]) -> Result<Vec<PillowCopy>> {
    let fm = face_map(hexes);
    let in_set = |h: u32| inside.get(h as usize).copied().unwrap_or(false);
    // Faces to pillow, oriented from the outside element into the set.
    let mut sheet: Vec<(u32, [u32; 4])> = Vec::new();
    for refs in fm.values() {
        if let [(a, fa), (b, fb)] = refs.as_slice() {
            match (in_set(*a), in_set(*b)) {
                (true, false) => sheet.push((*a, hexes.get(*a as usize).map_or([0; 4], |h| face(h, *fa)))),
                (false, true) => sheet.push((*b, hexes.get(*b as usize).map_or([0; 4], |h| face(h, *fb)))),
                _ => {}
            }
        }
    }
    if sheet.is_empty() {
        return Ok(Vec::new());
    }
    let mut on_sheet: Vec<u32> = sheet.iter().flat_map(|(_, q)| *q).collect();
    on_sheet.sort_unstable();
    on_sheet.dedup();
    // Set elements around each sheet vertex.
    let mut around: HashMap<u32, Vec<u32>> = HashMap::new();
    for (i, h) in hexes.iter().enumerate() {
        if !in_set(i as u32) {
            continue;
        }
        for &v in h {
            if on_sheet.binary_search(&v).is_ok() {
                let e = around.entry(v).or_default();
                if e.last() != Some(&(i as u32)) {
                    e.push(i as u32);
                }
            }
        }
    }
    let mut copies: Vec<PillowCopy> = Vec::new();
    // (vertex, element) → copy.
    let mut copy_of: HashMap<(u32, u32), u32> = HashMap::new();
    for &v in &on_sheet {
        let mut els = around.get(&v).cloned().unwrap_or_default();
        els.sort_unstable();
        els.dedup();
        // Group by faces inside the set that contain v.
        let mut group: Vec<usize> = (0..els.len()).collect();
        for i in 0..els.len() {
            for j in i + 1..els.len() {
                let (Some(hi), Some(hj)) = (hexes.get(els[i] as usize), hexes.get(els[j] as usize)) else { continue };
                let shares = (0..6u8).any(|f| {
                    let q = face(hi, f);
                    q.contains(&v) && (0..6u8).any(|g| key(face(hj, g)) == key(q))
                });
                if shares {
                    let (a, b) = (find(&mut group, i), find(&mut group, j));
                    group[a] = b;
                }
            }
        }
        // Per group: its copy and that copy's index in `copies`.
        let mut made: BTreeMap<usize, (u32, usize)> = BTreeMap::new();
        for (i, &e) in els.iter().enumerate() {
            let r = find(&mut group, i);
            let (c, at) = match made.get(&r) {
                Some(&m) => m,
                None => {
                    let c = u32::try_from(points.len()).map_err(|_| Error::Limit("too many points".into()))?;
                    points.push(points.get(v as usize).copied().unwrap_or([0.0; 3]));
                    copies.push(PillowCopy { original: v, copy: c, elements: Vec::new() });
                    made.insert(r, (c, copies.len() - 1));
                    (c, copies.len() - 1)
                }
            };
            copy_of.insert((v, e), c);
            if let Some(pc) = copies.get_mut(at) {
                pc.elements.push(e);
            }
        }
    }
    // The new layer: original face (seen from outside, normal into the set) below, copies above.
    let mut layer: Vec<[u32; 8]> = Vec::with_capacity(sheet.len());
    for (h, q) in &sheet {
        let r = [q[0], q[3], q[2], q[1]];
        let c = r.map(|v| copy_of.get(&(v, *h)).copied().unwrap_or(v));
        layer.push([r[0], r[1], r[2], r[3], c[0], c[1], c[2], c[3]]);
    }
    for (i, h) in hexes.iter_mut().enumerate() {
        if in_set(i as u32) {
            for v in h.iter_mut() {
                if let Some(&c) = copy_of.get(&(*v, i as u32)) {
                    *v = c;
                }
            }
        }
    }
    hexes.extend(layer);
    Ok(copies)
}

/// Hexahedra to drop so that every vertex of the boundary of `keep` can be pushed outward: the
/// normals of its boundary quads must lie strictly inside one half-space (their mean has a
/// positive product with each, at least `margin`). Where this fails, padding from the vertex
/// to the part's surface inverts elements, so for each such vertex the kept element at it with
/// the most boundary faces is dropped (the buffer-zone rule of H. Tong, E. Halilaj, Y. J. Zhang,
/// *HybridOctree_Hex*, arXiv 2401.05984, 2024, §2.3, in a mean-normal form).
pub(crate) fn hidden_vertices(hexes: &[[u32; 8]], keep: &[bool], xyz: &dyn Fn(u32) -> [f64; 3], margin: f64) -> Vec<u32> {
    let kept: Vec<u32> = (0..hexes.len() as u32).filter(|&i| keep.get(i as usize).copied().unwrap_or(false)).collect();
    let sub: Vec<[u32; 8]> = kept.iter().filter_map(|&i| hexes.get(i as usize).copied()).collect();
    let b = boundary(&sub);
    let mut normals: BTreeMap<u32, Vec<[f64; 3]>> = BTreeMap::new();
    let mut faces_of: BTreeMap<u32, usize> = BTreeMap::new();
    for ((h, _), q) in &b {
        let p = q.map(xyz);
        let d1 = [0, 1, 2].map(|k| p[2][k] - p[0][k]);
        let d2 = [0, 1, 2].map(|k| p[3][k] - p[1][k]);
        let n = [d1[1] * d2[2] - d1[2] * d2[1], d1[2] * d2[0] - d1[0] * d2[2], d1[0] * d2[1] - d1[1] * d2[0]];
        let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if l > 0.0 {
            for &v in q {
                normals.entry(v).or_default().push(n.map(|c| c / l));
            }
        }
        *faces_of.entry(*h).or_default() += 1;
    }
    let mut drop: Vec<u32> = Vec::new();
    for (v, ns) in &normals {
        let mut m = [0.0; 3];
        for n in ns {
            for k in 0..3 {
                m[k] += n[k];
            }
        }
        let l = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt();
        let ok = l > 0.0 && ns.iter().all(|n| (n[0] * m[0] + n[1] * m[1] + n[2] * m[2]) / l > margin);
        if ok {
            continue;
        }
        // The kept element at v with the most boundary faces (lowest index on ties).
        let best = sub
            .iter()
            .enumerate()
            .filter(|(_, h)| h.contains(v))
            .max_by(|a, b| faces_of.get(&(a.0 as u32)).unwrap_or(&0).cmp(faces_of.get(&(b.0 as u32)).unwrap_or(&0)).then(b.0.cmp(&a.0)))
            .and_then(|(i, _)| kept.get(i).copied());
        if let Some(e) = best {
            drop.push(e);
        }
    }
    drop.sort_unstable();
    drop.dedup();
    drop
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    #[test]
    fn faces_follow_the_model_hex8_face_order() {
        let model = fieldcraft_model::ElementKind::Hex8.faces();
        assert_eq!(model.len(), FACES.len());
        for (m, f) in model.iter().zip(FACES) {
            assert_eq!(*m, &f[..]);
        }
    }

    /// An `n × n × n` block of unit cubes.
    pub(crate) fn block(n: u32) -> (Vec<[f64; 3]>, Vec<[u32; 8]>) {
        let id = |x: u32, y: u32, z: u32| x + (n + 1) * (y + (n + 1) * z);
        let mut p = Vec::new();
        for z in 0..=n {
            for y in 0..=n {
                for x in 0..=n {
                    p.push([x as f64, y as f64, z as f64]);
                }
            }
        }
        let mut h = Vec::new();
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    h.push([
                        id(x, y, z),
                        id(x + 1, y, z),
                        id(x + 1, y + 1, z),
                        id(x, y + 1, z),
                        id(x, y, z + 1),
                        id(x + 1, y, z + 1),
                        id(x + 1, y + 1, z + 1),
                        id(x, y + 1, z + 1),
                    ]);
                }
            }
        }
        (p, h)
    }

    fn closed(hexes: &[[u32; 8]]) {
        for (k, v) in face_map(hexes) {
            assert!(v.len() <= 2, "face {k:?} on {} elements", v.len());
        }
        let b = boundary(hexes);
        for (_, qs) in quad_edges(&b.iter().map(|x| x.1).collect::<Vec<_>>()) {
            assert_eq!(qs.len(), 2);
        }
    }

    #[test]
    fn padding_and_pillowing_keep_a_closed_conforming_mesh() {
        let (mut p, mut h) = block(3);
        let b: Vec<[u32; 4]> = boundary(&h).into_iter().map(|x| x.1).collect();
        assert_eq!(b.len(), 54);
        let copy = pad(&mut p, &mut h, &b).unwrap();
        assert_eq!(copy.len(), 56);
        assert_eq!(h.len(), 27 + 54);
        closed(&h);
        // Pillow the padding hexes on the face z = 0 (outer faces down).
        let pad_start = 27;
        let inside: Vec<bool> = (0..h.len()).map(|i| i >= pad_start && (0..4).all(|k| p[h[i][k] as usize][2] == 0.0)).collect();
        assert_eq!(inside.iter().filter(|x| **x).count(), 9);
        let before = h.len();
        let copies = pillow(&mut p, &mut h, &inside).unwrap();
        closed(&h);
        // Sheet: the 9 inner faces plus the 12 side faces to the neighbouring padding.
        assert_eq!(h.len(), before + 9 + 12);
        // The boundary gains the ring: 54 + 12 quads.
        assert_eq!(boundary(&h).len(), 54 + 12);
        assert!(!copies.is_empty());
    }

    #[test]
    fn nonmanifold_edges_and_vertices_are_found() {
        let (_, h) = block(2);
        // Two cubes touching along an edge.
        let keep = |ids: &[usize]| (0..8).map(|i| ids.contains(&i)).collect::<Vec<bool>>();
        assert!(!nonmanifold(&h, &keep(&[0, 3])).is_empty(), "edge contact");
        assert!(!nonmanifold(&h, &keep(&[0, 7])).is_empty(), "vertex contact");
        assert!(nonmanifold(&h, &keep(&[0, 1, 3])).is_empty());
        assert!(nonmanifold(&h, &keep(&[0, 1, 2, 3, 4, 5, 6, 7])).is_empty());
    }
}
