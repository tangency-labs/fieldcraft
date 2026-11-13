//! Conforming hex–tet transitions: a layer of pyramids on a hexahedral mesh's interface faces,
//! and the boundary surface of the tetrahedral region that meets them.
//!
//! A hexahedron's face is a quadrilateral and a tetrahedron's a triangle, so a tet mesh cannot
//! share a hex mesh's faces. A pyramid can: its base is the quadrilateral and its four sides are
//! triangles. [`pyramid_surface`] puts one pyramid on every interface quadrilateral (its apex
//! inside the tet region, [`APEX_HEIGHT`] of the quadrilateral's mean edge above its centre) and builds
//! the surface the tetrahedra fill: the tet region's own surface without the interface faces,
//! closed by the pyramids' sides. Where the interface meets the tet region's other faces (the
//! rim), those faces' triangles are re-stitched to the hex mesh's rim nodes: each rim vertex of
//! the tessellation moves to the nearest hex rim node and each rim edge becomes a fan over the hex
//! rim nodes it spans, so the rim of the surface is exactly the hex mesh's. The tet mesher then
//! keeps the pyramids' sides unchanged ([`crate::tet::QualityOptions::frozen_faces`]) and
//! [`apex_points`] checks that it did.
//!
//! The result is one conforming mesh: hexahedra, pyramids on their faces, tetrahedra on the
//! pyramids' sides and on the rest of the region. The surface is checked by
//! [`BoundarySurface::new`] like any other (closed, manifold, no self-intersection): a pyramid
//! that reaches through a thin region, or a rim that cannot be stitched, is an error, and the
//! caller falls back to another join.
//!
//! # Sources
//! - S. J. Owen, S. A. Canann, S. Saigal, *Pyramid elements for maintaining tetrahedra to
//!   hexahedra conformability*, AMD Vol. 220, Trends in Unstructured Mesh Generation, ASME 1997,
//!   pp. 123–129: a pyramid on each quadrilateral of the interface, its apex placed off the face
//!   inside the tetrahedral region.
//! - L. Bedrosian, *Shape functions and integration formulas for three-dimensional finite
//!   element analysis*, IJNME 35(1):95–108, 1992: the pyramid5 element that makes the
//!   transition conforming (`fieldcraft_fem::basis`).

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use crate::surface::BoundarySurface;
use crate::tet::QualityMesh;
use crate::tet::complex::face_toward;
use crate::{Error, Result};

/// The apex's height above the quadrilateral's centre, over the quadrilateral's mean edge. Below
/// one half, so the sides of two pyramids on faces that meet at a right angle (a box-shaped tet
/// cell's corner between hex cells) lean apart instead of crossing.
pub const APEX_HEIGHT: f64 = 0.4;

/// The tet region's surface for a pyramid transition ([`pyramid_surface`]).
#[derive(Debug, Clone)]
pub struct Transition {
    pub positions: Vec<[f64; 3]>,
    pub triangles: Vec<[u32; 3]>,
    pub tri_face: Vec<u32>,
    /// One apex per interface quadrilateral, in their order.
    pub apexes: Vec<[f64; 3]>,
    /// Each pyramid's four side triangles (positions), outward from the tet region, quad by quad.
    pub sides: Vec<[[f64; 3]; 3]>,
}

fn key(p: [f64; 3]) -> [u64; 3] {
    p.map(f64::to_bits)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub(a, b);
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

/// The tet region's surface with a pyramid on every quadrilateral of `quads` (corners in order,
/// their right-hand normal pointing into the tet region, i.e. out of the hexahedra), in place of
/// the triangles of `surface` on the B-rep faces `cut`. The pyramids' sides are tagged `tag`.
pub fn pyramid_surface(surface: &BoundarySurface, cut: &[u32], quads: &[[[f64; 3]; 4]], tag: u32) -> Result<Transition> {
    if quads.is_empty() {
        return Err(Error::Invalid("no interface quadrilaterals".into()));
    }
    let cut: BTreeSet<u32> = cut.iter().copied().collect();
    if cut.contains(&tag) || surface.tri_face().contains(&tag) {
        return Err(Error::Invalid(format!("the pyramid tag {tag} is a face of the region")));
    }
    let mut positions: Vec<[f64; 3]> = surface.points().to_vec();
    let (mut kept, mut removed): (Vec<([u32; 3], u32)>, Vec<[u32; 3]>) = (Vec::new(), Vec::new());
    for (t, &f) in surface.triangles().iter().zip(surface.tri_face()) {
        if cut.contains(&f) {
            removed.push(*t);
        } else {
            kept.push((*t, f));
        }
    }
    if removed.is_empty() {
        return Err(Error::Invalid("the interface faces have no triangles".into()));
    }
    // The hex nodes of the interface, after the surface's points.
    let mut hex_index: BTreeMap<[u64; 3], u32> = BTreeMap::new();
    let mut quad_nodes: Vec<[u32; 4]> = Vec::with_capacity(quads.len());
    for q in quads {
        if q.iter().flatten().any(|c| !c.is_finite()) {
            return Err(Error::Invalid("an interface quadrilateral has a non-finite corner".into()));
        }
        let mut ids = [0u32; 4];
        for (slot, p) in ids.iter_mut().zip(q) {
            let next = u32::try_from(positions.len()).map_err(|_| Error::Limit("too many points".into()))?;
            *slot = *hex_index.entry(key(*p)).or_insert_with(|| {
                positions.push(*p);
                next
            });
        }
        quad_nodes.push(ids);
    }
    // The rim of the quadrilaterals: edges of one quadrilateral, as a graph.
    let mut edge_use: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for q in &quad_nodes {
        for i in 0..4 {
            let (a, b) = (q[i], q[(i + 1) % 4]);
            *edge_use.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    if edge_use.values().any(|&n| n > 2) {
        return Err(Error::Invalid("an edge of the interface is shared by more than two quadrilaterals".into()));
    }
    let mut rim: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for (&(a, b), _) in edge_use.iter().filter(|(_, n)| **n == 1) {
        rim.entry(a).or_default().push(b);
        rim.entry(b).or_default().push(a);
    }
    if rim.values().any(|n| n.len() != 2) {
        return Err(Error::Unsupported("the interface's rim is not a set of simple loops (it touches itself at a node)".into()));
    }
    // The surface's rim: the kept triangles' edges whose reverse is a removed triangle's edge.
    let removed_edges: BTreeSet<(u32, u32)> = removed.iter().flat_map(|t| [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])]).collect();
    let removed_vertices: BTreeSet<u32> = removed.iter().flatten().copied().collect();
    let rim_nodes: Vec<u32> = rim.keys().copied().collect();
    let mut nearest: BTreeMap<u32, u32> = BTreeMap::new();
    for &(t, _) in &kept {
        for v in t {
            if !removed_vertices.contains(&v) || nearest.contains_key(&v) {
                continue;
            }
            let p = positions.get(v as usize).copied().ok_or_else(|| Error::Failed("a surface vertex is out of range".into()))?;
            let mut best: Option<(f64, u32)> = None;
            for &h in &rim_nodes {
                let d = positions.get(h as usize).map_or(f64::INFINITY, |q| dist(p, *q));
                if best.is_none_or(|b| d < b.0) {
                    best = Some((d, h));
                }
            }
            if let Some((_, h)) = best {
                nearest.insert(v, h);
            }
        }
    }
    let at = |v: u32| nearest.get(&v).copied().unwrap_or(v);
    let mut triangles: Vec<[u32; 3]> = Vec::with_capacity(kept.len() + 4 * quads.len());
    let mut tri_face: Vec<u32> = Vec::with_capacity(triangles.capacity());
    for &(t, f) in &kept {
        let rim_edges: Vec<usize> = (0..3).filter(|&i| removed_edges.contains(&(t[(i + 1) % 3], t[i]))).collect();
        match rim_edges.as_slice() {
            [] => {
                let m = t.map(at);
                if m[0] == m[1] || m[1] == m[2] || m[2] == m[0] {
                    return Err(Error::Unsupported(
                        "a triangle beside the interface collapses onto its rim (the tessellation is coarser than the hex mesh there)".into(),
                    ));
                }
                triangles.push(m);
                tri_face.push(f);
            }
            [i] => {
                let (u, v, w) = (t[*i], t[(i + 1) % 3], t[(i + 2) % 3]);
                let path = rim_path(&rim, at(u), at(v))?;
                let apex = at(w);
                for pair in path.windows(2) {
                    if let [a, b] = *pair {
                        if apex == a || apex == b {
                            return Err(Error::Unsupported("a triangle beside the interface lies along the rim on two sides".into()));
                        }
                        triangles.push([a, b, apex]);
                        tri_face.push(f);
                    }
                }
            }
            [i, j] => {
                // The triangle spans a corner of the rim: its corner is the vertex the two rim
                // edges share, its edges a → corner → b both lie on the rim and b → a is inside.
                let corner = if j - i == 1 { *j } else { *i };
                let (a, c, b) = (t[(corner + 2) % 3], t[corner], t[(corner + 1) % 3]);
                let mut toward_a = rim_path(&rim, at(a), at(c))?;
                toward_a.reverse();
                let toward_b = rim_path(&rim, at(c), at(b))?;
                corner_strip(&toward_a, &toward_b, &positions, &mut triangles)?;
                tri_face.resize(triangles.len(), f);
            }
            _ => return Err(Error::Unsupported("a triangle beside the interface lies along the rim on all three edges".into())),
        }
    }
    // The pyramids: an apex above each quadrilateral, the sides toward the tet region.
    let mut apexes = Vec::with_capacity(quads.len());
    let mut sides = Vec::with_capacity(4 * quads.len());
    for (q, ids) in quads.iter().zip(&quad_nodes) {
        // Newell normal (into the tet region) and the mean edge.
        let mut n = [0.0; 3];
        let mut edge = 0.0;
        for i in 0..4 {
            let (p, r) = (q[i], q[(i + 1) % 4]);
            n = [n[0] + (p[1] - r[1]) * (p[2] + r[2]), n[1] + (p[2] - r[2]) * (p[0] + r[0]), n[2] + (p[0] - r[0]) * (p[1] + r[1])];
            edge += dist(p, r) / 4.0;
        }
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if !(len > 0.0 && edge > 0.0 && len.is_finite()) {
            return Err(Error::Invalid("an interface quadrilateral is degenerate".into()));
        }
        let c = [0, 1, 2].map(|k| q.iter().map(|p| p[k]).sum::<f64>() / 4.0);
        let apex = [0, 1, 2].map(|k| c[k] + n[k] / len * APEX_HEIGHT * edge);
        let a = u32::try_from(positions.len()).map_err(|_| Error::Limit("too many points".into()))?;
        positions.push(apex);
        apexes.push(apex);
        for i in 0..4 {
            // The pyramid's side on edge i → i+1 faces out of the pyramid; the region's boundary
            // there faces into it.
            let (p, r) = (ids[i], ids[(i + 1) % 4]);
            triangles.push([r, p, a]);
            tri_face.push(tag);
            sides.push([q[(i + 1) % 4], q[i], apex]);
        }
    }
    // Drop the points no triangle uses (the interface's own and the moved rim vertices).
    let mut used = vec![u32::MAX; positions.len()];
    let mut compact = Vec::new();
    for t in &mut triangles {
        for v in t.iter_mut() {
            let slot = used.get_mut(*v as usize).ok_or_else(|| Error::Failed("a vertex is out of range".into()))?;
            if *slot == u32::MAX {
                *slot = compact.len() as u32;
                compact.push(positions.get(*v as usize).copied().unwrap_or([0.0; 3]));
            }
            *v = *slot;
        }
    }
    Ok(Transition { positions: compact, triangles, tri_face, apexes, sides })
}

/// Triangulate the region between two rim paths that leave the same corner node (`left[0] ==
/// right[0]`) and end at the two ends of an inside edge: the strip of triangles from the corner
/// out to that edge, each step along the path whose next diagonal is shorter. The triangles turn
/// the way of the triangle (left end, corner, right end) they replace, and the last one's edge
/// right end → left end is the inside edge its neighbour shares.
fn corner_strip(left: &[u32], right: &[u32], positions: &[[f64; 3]], triangles: &mut Vec<[u32; 3]>) -> Result<()> {
    let at = |v: u32| positions.get(v as usize).copied().ok_or_else(|| Error::Failed("a rim node is out of range".into()));
    if left.len() < 2 || right.len() < 2 || left.first() != right.first() {
        return Err(Error::Unsupported(
            "a triangle at a corner of the interface collapses onto its rim (the tessellation is coarser than the hex mesh there)".into(),
        ));
    }
    let mut nodes: BTreeSet<u32> = left.iter().copied().collect();
    if right.iter().skip(1).any(|v| !nodes.insert(*v)) {
        return Err(Error::Unsupported("a triangle at a corner of the interface spans the rim on both sides of the corner".into()));
    }
    let (corner, (mut i, mut j)) = (left[0], (1, 1));
    triangles.push([left[1], corner, right[1]]);
    while i + 1 < left.len() || j + 1 < right.len() {
        let along_left = match (left.get(i + 1), right.get(j + 1)) {
            (Some(&l), Some(&r)) => dist(at(l)?, at(right[j])?) <= dist(at(left[i])?, at(r)?),
            (l, _) => l.is_some(),
        };
        if along_left {
            triangles.push([left[i], right[j], left[i + 1]]);
            i += 1;
        } else {
            triangles.push([left[i], right[j], right[j + 1]]);
            j += 1;
        }
    }
    Ok(())
}

/// The rim nodes from `from` to `to` along the shorter way round their loop.
fn rim_path(rim: &BTreeMap<u32, Vec<u32>>, from: u32, to: u32) -> Result<Vec<u32>> {
    if from == to {
        return Ok(vec![from]);
    }
    // Breadth-first: the loop's two ways meet at the far side, so this is the shorter one.
    let mut previous: BTreeMap<u32, u32> = BTreeMap::new();
    let mut queue = VecDeque::from([from]);
    previous.insert(from, from);
    while let Some(v) = queue.pop_front() {
        if v == to {
            break;
        }
        for &w in rim.get(&v).map(Vec::as_slice).unwrap_or_default() {
            if let std::collections::btree_map::Entry::Vacant(e) = previous.entry(w) {
                e.insert(v);
                queue.push_back(w);
            }
        }
    }
    let mut path = vec![to];
    let mut v = to;
    while v != from {
        v = *previous.get(&v).ok_or_else(|| Error::Unsupported("the tessellation's rim and the hex mesh's rim are on different loops".into()))?;
        path.push(v);
        if path.len() > rim.len() + 1 {
            return Err(Error::Failed("the rim walk did not end".into()));
        }
    }
    path.reverse();
    Ok(path)
}

/// The tet mesh's point index of each pyramid's apex, after checking that the mesh kept every
/// pyramid side as one of its boundary faces (so the pyramids and tetrahedra share them).
pub fn apex_points(mesh: &QualityMesh, t: &Transition) -> Result<Vec<u32>> {
    let index: BTreeMap<[u64; 3], u32> = mesh.points.iter().enumerate().map(|(i, p)| (key(*p), i as u32)).collect();
    let find = |p: [f64; 3]| index.get(&key(p)).copied();
    let mut faces: BTreeSet<[u32; 3]> = BTreeSet::new();
    for b in &mesh.boundary {
        let Some(e) = mesh.elements.get(b.tet as usize) else { continue };
        let Some(vs) = e.get(..4).and_then(|s| <[u32; 4]>::try_from(s).ok()) else { continue };
        let mut f = face_toward(vs, usize::from(b.opposite));
        f.sort_unstable();
        faces.insert(f);
    }
    let mut missing = 0usize;
    for s in &t.sides {
        match (find(s[0]), find(s[1]), find(s[2])) {
            (Some(a), Some(b), Some(c)) => {
                let mut f = [a, b, c];
                f.sort_unstable();
                missing += usize::from(!faces.contains(&f));
            }
            _ => missing += 1,
        }
    }
    if missing > 0 {
        return Err(Error::Failed(format!("the tet mesh changed {missing} of the {} pyramid sides it had to keep", t.sides.len())));
    }
    t.apexes.iter().map(|p| find(*p).ok_or_else(|| Error::Failed("a pyramid's apex is not a vertex of the tet mesh".into()))).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cdt::{CdtOptions, tetrahedralize};
    use crate::surface::SurfaceOptions;
    use crate::tet::{QualityOptions, quality_mesh};

    /// A box [0, 2] × [0, 2] × [0, height] as twelve triangles (two per side, outward), sides tagged
    /// 0..6 (−x, +x, −y, +y, −z, +z).
    fn slab(height: f64) -> BoundarySurface {
        let p: Vec<[f64; 3]> = (0..8).map(|i| [2.0 * (i & 1) as f64, 2.0 * ((i >> 1) & 1) as f64, height * ((i >> 2) & 1) as f64]).collect();
        let quads = [([0, 4, 6, 2], 0), ([1, 3, 7, 5], 1), ([0, 1, 5, 4], 2), ([2, 6, 7, 3], 3), ([0, 2, 3, 1], 4), ([4, 5, 7, 6], 5)];
        let (mut t, mut f) = (Vec::new(), Vec::new());
        for (q, tag) in quads {
            t.push([q[0], q[1], q[2]]);
            t.push([q[0], q[2], q[3]]);
            f.extend([tag, tag]);
        }
        BoundarySurface::new(&p, &t, &f, &SurfaceOptions::default()).unwrap()
    }

    /// The bottom (z = 0) as a 4 × 4 grid of hex faces, their normals up into the slab.
    fn bottom_quads() -> Vec<[[f64; 3]; 4]> {
        let mut q = Vec::new();
        for j in 0..4 {
            for i in 0..4 {
                let (x, y, h) = (0.5 * i as f64, 0.5 * j as f64, 0.5);
                q.push([[x, y, 0.0], [x + h, y, 0.0], [x + h, y + h, 0.0], [x, y + h, 0.0]]);
            }
        }
        q
    }

    #[test]
    fn pyramids_on_a_face_are_kept_by_the_tet_mesh() {
        let s = slab(1.0);
        let quads = bottom_quads();
        let t = pyramid_surface(&s, &[4], &quads, 6).unwrap();
        assert_eq!(t.apexes.len(), 16);
        let surface = BoundarySurface::new(&t.positions, &t.triangles, &t.tri_face, &SurfaceOptions::default()).unwrap();
        // The slab less the pyramids (each 0.5² × 0.3 / 3).
        let pyramids = 16.0 * 0.25 * (APEX_HEIGHT * 0.5) / 3.0;
        assert!((surface.volume() - (4.0 - pyramids)).abs() < 1e-12, "{}", surface.volume());
        let cdt = tetrahedralize(&surface, &CdtOptions::default()).unwrap();
        let options = QualityOptions { size: Some(0.3), frozen_faces: vec![6], ..QualityOptions::default() };
        let mesh = quality_mesh(&surface, &cdt, Vec::new(), &options).unwrap();
        let apexes = apex_points(&mesh, &t).unwrap();
        assert_eq!(apexes.len(), 16);
        // Without freezing, refinement splits the pyramid sides and the check says so.
        let free = quality_mesh(&surface, &cdt, Vec::new(), &QualityOptions { size: Some(0.1), ..QualityOptions::default() }).unwrap();
        assert!(apex_points(&free, &t).is_err());
    }

    #[test]
    fn a_triangle_across_a_corner_of_the_rim_is_stitched_along_both_edges() {
        // The cube [0, 2]³ with hex faces on its bottom and its −x side. The −y and +y sides are
        // split by the diagonal that leaves a triangle along both of those faces at their corner.
        let p: Vec<[f64; 3]> = (0..8).map(|i| [2.0 * (i & 1) as f64, 2.0 * ((i >> 1) & 1) as f64, 2.0 * ((i >> 2) & 1) as f64]).collect();
        let t = [[0, 4, 6], [0, 6, 2], [1, 3, 7], [1, 7, 5], [0, 1, 4], [1, 5, 4], [2, 6, 3], [6, 7, 3], [0, 2, 3], [0, 3, 1], [4, 5, 7], [4, 7, 6]];
        let s = BoundarySurface::new(&p, &t, &[0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5], &SurfaceOptions::default()).unwrap();
        // The −x side's faces, their normals along +x into the cube.
        let side = (0..16)
            .map(|k| ((k % 4) as f64 * 0.5, (k / 4) as f64 * 0.5))
            .map(|(y, z)| [[0.0, y, z], [0.0, y + 0.5, z], [0.0, y + 0.5, z + 0.5], [0.0, y, z + 0.5]]);
        let quads: Vec<[[f64; 3]; 4]> = bottom_quads().into_iter().chain(side).collect();
        let t = pyramid_surface(&s, &[4, 0], &quads, 6).unwrap();
        let surface = BoundarySurface::new(&t.positions, &t.triangles, &t.tri_face, &SurfaceOptions::default()).unwrap();
        let pyramids = 32.0 * 0.25 * (APEX_HEIGHT * 0.5) / 3.0;
        assert!((surface.volume() - (8.0 - pyramids)).abs() < 1e-12, "{}", surface.volume());
        let cdt = tetrahedralize(&surface, &CdtOptions::default()).unwrap();
        let options = QualityOptions { size: Some(0.4), frozen_faces: vec![6], ..QualityOptions::default() };
        let mesh = quality_mesh(&surface, &cdt, Vec::new(), &options).unwrap();
        assert_eq!(apex_points(&mesh, &t).unwrap().len(), 32);
    }

    #[test]
    fn a_pyramid_through_a_thin_region_is_refused() {
        // 1 × 1 quads: apexes 0.4 up, inside a 1 mm slab, through a 0.3 mm one.
        let quads: Vec<[[f64; 3]; 4]> = (0..4)
            .map(|k| ((k % 2) as f64, (k / 2) as f64))
            .map(|(x, y)| [[x, y, 0.0], [x + 1.0, y, 0.0], [x + 1.0, y + 1.0, 0.0], [x, y + 1.0, 0.0]])
            .collect();
        for (height, ok) in [(1.0, true), (0.3, false)] {
            let t = pyramid_surface(&slab(height), &[4], &quads, 6).unwrap();
            assert_eq!(BoundarySurface::new(&t.positions, &t.triangles, &t.tri_face, &SurfaceOptions::default()).is_ok(), ok, "{height}");
        }
    }
}
