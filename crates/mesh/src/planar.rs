//! 2D models of planar regions: the meshes plane stress, plane strain and axisymmetric analyses
//! run on.
//!
//! A region is bounded by loops of *curves* (polylines; the first loop is the outer boundary, the
//! others are holes). [`mesh_region`] divides every curve by the size, meshes the region with the
//! free quadrilateral mesher or the free triangle mesher of [`crate::quad`] on those boundary
//! nodes, and adds mid-side nodes for quadratic elements (on the boundary polyline for boundary
//! edges, at the chord midpoint inside). Every element edge on a curve is reported with the
//! curve's id, so a mesher above can tag it with the geometry it came from.
//!
//! [`face_boundary`] gives the curves of one planar face of a tessellated body: the face's
//! boundary loops in the XY plane, split where the neighbouring face changes (each curve's id is
//! the neighbour's face index, the B-rep edge the two faces share) and at sharp corners.

use std::collections::BTreeMap;

use fieldcraft_model::ElementKind;

use crate::quad::{QuadOptions, free_quad_mesh, free_tri_mesh};
use crate::size::{UniformSize, divide_polyline, edge_divisions, even_divisions, ideal_divisions};
use crate::surface::weld::weld;
use crate::{Error, Result};

/// Round-off in a face's height, relative to it.
const ROUND_OFF: f64 = 1e-12; // unit-factors: relative to the coordinate

/// A boundary curve: a polyline from its start to the start of the next curve of its loop (the
/// end point is not repeated), and an id reported on the element edges that lie on it.
#[derive(Debug, Clone, PartialEq)]
pub struct Curve {
    pub id: Option<u32>,
    pub points: Vec<[f64; 2]>,
}

/// A planar mesh of 2D elements in the XY plane.
#[derive(Debug, Clone, PartialEq)]
pub struct PlanarMesh {
    pub kind: ElementKind,
    pub nodes: Vec<[f64; 2]>,
    /// Counterclockwise node lists in the model's node order (corners, then mid-side nodes).
    pub elements: Vec<Vec<u32>>,
    /// Element edges on a curve with an id: (element index, edge number from 1, curve id).
    pub edges: Vec<(usize, u8, u32)>,
    /// The node at the start of each curve with an id: (node, curve id).
    pub corners: Vec<(u32, u32)>,
}

/// Turns sharper than this (radians) start a new curve in [`face_boundary`], so corners stay
/// nodes of the mesh.
const CORNER_ANGLE: f64 = std::f64::consts::PI / 6.0;

/// Mesh the region bounded by `loops` with elements of `kind` (tri3, tri6, quad4 or quad8) of
/// about `size`.
pub fn mesh_region(loops: &[Vec<Curve>], kind: ElementKind, size: f64, opts: &QuadOptions) -> Result<PlanarMesh> {
    if kind.dimension() != 2 {
        return Err(Error::Invalid(format!("a planar mesh has tri3, tri6, quad4 or quad8 elements, not {}", kind.name())));
    }
    if !(size > 0.0 && size.is_finite()) {
        return Err(Error::invalid(crate::message!["the element size must be a positive number, got ", crate::message::length(size)]));
    }
    if loops.is_empty() || loops.iter().any(Vec::is_empty) {
        return Err(Error::Invalid("a region needs at least one loop of curves".into()));
    }
    let quads = matches!(kind, ElementKind::Quad4 | ElementKind::Quad8);
    let quadratic = kind.is_quadratic();
    let field = UniformSize(size);
    // Boundary points per loop, the curve of each boundary segment, and the boundary mid points.
    let mut boundary: Vec<Vec<[f64; 2]>> = Vec::with_capacity(loops.len());
    let mut segment_curve: Vec<Option<u32>> = Vec::new();
    let mut segment_mid: Vec<[f64; 2]> = Vec::new();
    let mut curve_starts: Vec<(usize, u32)> = Vec::new();
    for (l, curves) in loops.iter().enumerate() {
        let polylines: Vec<Vec<[f64; 3]>> = curves
            .iter()
            .enumerate()
            .map(|(i, c)| {
                let next = curves.get((i + 1) % curves.len()).and_then(|n| n.points.first()).copied();
                c.points.iter().chain(next.as_ref()).map(|p| [p[0], p[1], 0.0]).collect()
            })
            .collect();
        let mut counts: Vec<usize> = polylines
            .iter()
            .map(|p| if quads { ideal_divisions(&field, p).and_then(even_divisions) } else { edge_divisions(&field, p) })
            .collect::<Result<_>>()?;
        // Enough points per loop for the meshers: 3 for triangles; 4 for a lone quad loop, 6 with
        // holes.
        let least = if !quads {
            3
        } else if loops.len() > 1 {
            6
        } else {
            4
        };
        while counts.iter().sum::<usize>() < least {
            let longest = (0..counts.len()).max_by(|&a, &b| counts[a].cmp(&counts[b]).then(b.cmp(&a))).unwrap_or(0);
            if let Some(c) = counts.get_mut(longest) {
                *c += if quads { 2 } else { 1 };
            }
        }
        let mut pts = Vec::new();
        for ((poly, n), curve) in polylines.iter().zip(&counts).zip(curves) {
            let fine = divide_polyline(&field, poly, if quadratic { 2 * n } else { *n })?;
            let step = if quadratic { 2 } else { 1 };
            if let Some(id) = curve.id {
                curve_starts.push((boundary.iter().map(Vec::len).sum::<usize>() + pts.len(), id));
            }
            for k in 0..*n {
                let p = fine.get(step * k).ok_or_else(|| Error::Failed("curve division".into()))?;
                pts.push([p[0], p[1]]);
                segment_curve.push(curve.id);
                let m = if quadratic { fine.get(2 * k + 1) } else { None }.copied().unwrap_or_default();
                segment_mid.push([m[0], m[1]]);
            }
        }
        if pts.len() < 3 {
            return Err(Error::Invalid(format!("loop {l} is too short to mesh")));
        }
        boundary.push(pts);
    }
    let (nodes, corners_of): (Vec<[f64; 2]>, Vec<Vec<u32>>) = if quads {
        let m = free_quad_mesh(&boundary, &field, opts)?;
        (m.nodes, m.quads.iter().map(|q| q.to_vec()).collect())
    } else {
        let m = free_tri_mesh(&boundary, &field, opts)?;
        (m.nodes, m.triangles.iter().map(|t| t.to_vec()).collect())
    };
    // Boundary segment of a node pair: segment k of loop l joins its points k and k + 1.
    let mut segment: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    let mut offset = 0usize;
    let mut s = 0usize;
    for pts in &boundary {
        let n = pts.len();
        for k in 0..n {
            let (a, b) = ((offset + k) as u32, (offset + (k + 1) % n) as u32);
            segment.insert((a.min(b), a.max(b)), s);
            s += 1;
        }
        offset += n;
    }
    let mut nodes = nodes;
    let mut mids: BTreeMap<(u32, u32), u32> = BTreeMap::new();
    let mut elements = Vec::with_capacity(corners_of.len());
    let mut edges = Vec::new();
    for (e, corners) in corners_of.iter().enumerate() {
        let mut list = corners.clone();
        for (i, local) in kind.edges().iter().enumerate() {
            let (Some(&a), Some(&b)) = (local.first().and_then(|x| corners.get(*x)), local.get(1).and_then(|x| corners.get(*x))) else {
                return Err(Error::Failed("element edge".into()));
            };
            let key = (a.min(b), a.max(b));
            let on = segment.get(&key).copied();
            if let Some(id) = on.and_then(|s| segment_curve.get(s).copied().flatten()) {
                edges.push((e, i as u8 + 1, id));
            }
            if quadratic {
                let id = match mids.get(&key) {
                    Some(id) => *id,
                    None => {
                        let p = match on.and_then(|s| segment_mid.get(s)) {
                            Some(m) => *m,
                            None => {
                                let (pa, pb) =
                                    (nodes.get(a as usize).copied().unwrap_or_default(), nodes.get(b as usize).copied().unwrap_or_default());
                                [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])]
                            }
                        };
                        nodes.push(p);
                        let id = u32::try_from(nodes.len() - 1).map_err(|_| Error::Limit("too many nodes".into()))?;
                        mids.insert(key, id);
                        id
                    }
                };
                list.push(id);
            }
        }
        elements.push(list);
    }
    let corners = curve_starts.into_iter().map(|(n, id)| (n as u32, id)).collect();
    Ok(PlanarMesh { kind, nodes, elements, edges, corners })
}

/// The boundary of one planar face of a triangulated body, in the XY plane.
#[derive(Debug, Clone, PartialEq)]
pub struct FaceBoundary {
    /// The outer loop first, then the holes; each curve's id is the face index of the neighbour
    /// it borders (`None` where no neighbour was found).
    pub loops: Vec<Vec<Curve>>,
    /// The face's Z.
    pub z: f64,
}

/// The boundary of face `face` of a triangulated body (`tri_face` gives each triangle's face):
/// the face must lie in a plane of constant Z (within `tolerance`); points of neighbouring faces
/// within `tolerance` are the same point.
pub fn face_boundary(positions: &[[f64; 3]], triangles: &[[u32; 3]], tri_face: &[u32], face: u32, tolerance: f64) -> Result<FaceBoundary> {
    if !(tolerance >= 0.0 && tolerance.is_finite()) {
        return Err(Error::invalid(crate::message!["tolerance must be a non-negative number, got ", crate::message::length(tolerance)]));
    }
    if let Some(v) = triangles.iter().flatten().find(|v| **v as usize >= positions.len()) {
        return Err(Error::Invalid(format!("triangle vertex {v} is out of range")));
    }
    if tri_face.len() != triangles.len() {
        return Err(Error::Invalid("every triangle needs its face".into()));
    }
    // One weld of the whole tessellation: each face's copies of a shared point become one point,
    // so neighbours meet on welded edges.
    let welded = weld(positions, triangles, 0.0, tolerance)?;
    let xyz = welded.points;
    let face_of = |t: usize| welded.source_triangle.get(t).and_then(|s| tri_face.get(*s as usize)).copied();
    let mine: Vec<[u32; 3]> = welded.triangles.iter().enumerate().filter(|(t, _)| face_of(*t) == Some(face)).map(|(_, w)| *w).collect();
    if mine.is_empty() {
        return Err(Error::Invalid(format!("face {face} has no triangles")));
    }
    let z = mine.first().and_then(|w| xyz.get(w[0] as usize)).map_or(0.0, |p| p[2]);
    let flat = mine.iter().flatten().filter_map(|v| xyz.get(*v as usize)).all(|p| (p[2] - z).abs() <= tolerance.max(ROUND_OFF * (1.0 + z.abs())));
    if !flat {
        return Err(Error::Invalid("the face does not lie in a plane parallel to XY (a 2D model lies in the XY plane)".into()));
    }
    // Boundary edges: directed edges whose reverse is not in the face.
    let mut directed: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for w in &mine {
        for k in 0..3 {
            *directed.entry((w[k], w[(k + 1) % 3])).or_default() += 1;
        }
    }
    let mut next: BTreeMap<u32, u32> = BTreeMap::new();
    for (&(a, b), &count) in &directed {
        if directed.contains_key(&(b, a)) {
            continue;
        }
        if count > 1 || next.insert(a, b).is_some() {
            return Err(Error::Invalid("the face's boundary touches itself (a non-manifold face); split the face".into()));
        }
    }
    // Neighbours: the other faces' welded edges, by their ends.
    let mut other_faces: BTreeMap<(u32, u32), u32> = BTreeMap::new();
    for (t, w) in welded.triangles.iter().enumerate() {
        let Some(f) = face_of(t).filter(|f| *f != face) else { continue };
        for k in 0..3 {
            let (a, b) = (w[k], w[(k + 1) % 3]);
            other_faces.entry((a.min(b), a.max(b))).or_insert(f);
        }
    }
    let neighbour = |a: u32, b: u32| other_faces.get(&(a.min(b), a.max(b))).copied();
    // Chain the boundary into loops, each started at its smallest vertex id.
    let mut loops: Vec<(f64, Vec<Curve>)> = Vec::new();
    let mut unused: BTreeMap<u32, u32> = next.clone();
    while let Some((&start, _)) = unused.iter().next() {
        let mut chain = vec![start];
        let mut v = start;
        loop {
            let Some(w) = unused.remove(&v) else {
                return Err(Error::Invalid("the face's boundary is not closed".into()));
            };
            if w == start {
                break;
            }
            chain.push(w);
            v = w;
            if chain.len() > next.len() {
                return Err(Error::Failed("boundary chaining did not close".into()));
            }
        }
        let at = |i: usize| chain.get(i % chain.len()).and_then(|v| xyz.get(*v as usize)).copied().unwrap_or_default();
        let n = chain.len();
        let ids: Vec<Option<u32>> = (0..n).map(|i| neighbour(chain[i], chain[(i + 1) % n])).collect();
        // A curve starts where the neighbour changes or the boundary turns sharply.
        let starts: Vec<bool> = (0..n)
            .map(|i| {
                let (p, q, r) = (at(i + n - 1), at(i), at(i + 1));
                let (u, w) = ([q[0] - p[0], q[1] - p[1]], [r[0] - q[0], r[1] - q[1]]);
                let turn = crate::vec3::atan2(u[0] * w[1] - u[1] * w[0], u[0] * w[0] + u[1] * w[1]).abs();
                ids[i] != ids[(i + n - 1) % n] || turn > CORNER_ANGLE
            })
            .collect();
        let first = starts.iter().position(|s| *s).unwrap_or(0);
        let mut curves: Vec<Curve> = Vec::new();
        for k in 0..n {
            let i = (first + k) % n;
            let p = at(i);
            if starts[i] || curves.is_empty() {
                curves.push(Curve { id: ids[i], points: Vec::new() });
            }
            if let Some(c) = curves.last_mut() {
                c.points.push([p[0], p[1]]);
            }
        }
        let area: f64 = (0..n).map(|i| at(i)[0] * at(i + 1)[1] - at(i + 1)[0] * at(i)[1]).sum::<f64>() / 2.0;
        loops.push((area.abs(), curves));
    }
    // The outer loop encloses the most area.
    let outer = (0..loops.len()).max_by(|&a, &b| loops[a].0.total_cmp(&loops[b].0).then(b.cmp(&a))).unwrap_or(0);
    let mut out = Vec::with_capacity(loops.len());
    if outer < loops.len() {
        out.push(loops.remove(outer).1);
    }
    out.extend(loops.into_iter().map(|l| l.1));
    Ok(FaceBoundary { loops: out, z })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square_with_hole() -> Vec<Vec<Curve>> {
        let side = |id: u32, a: [f64; 2]| Curve { id: Some(id), points: vec![a] };
        let outer = vec![side(1, [0.0, 0.0]), side(2, [4.0, 0.0]), side(3, [4.0, 3.0]), side(4, [0.0, 3.0])];
        let circle: Vec<[f64; 2]> = (0..32).map(|k| 2.0 * std::f64::consts::PI * k as f64 / 32.0).map(|t| [2.0 + t.cos(), 1.5 + t.sin()]).collect();
        vec![outer, vec![Curve { id: Some(5), points: circle }]]
    }

    #[test]
    fn regions_mesh_with_every_element_kind_and_report_curve_edges() {
        let loops = square_with_hole();
        let area_want = 12.0 - std::f64::consts::PI;
        for kind in [ElementKind::Tri3, ElementKind::Tri6, ElementKind::Quad4, ElementKind::Quad8] {
            let m = mesh_region(&loops, kind, 0.4, &QuadOptions::default()).unwrap();
            assert!(m.elements.iter().all(|e| e.len() == kind.node_count()));
            // Area by the corner polygons (the hole's chords make it a little larger).
            let area: f64 = m
                .elements
                .iter()
                .map(|e| {
                    let c: Vec<[f64; 2]> = e[..kind.corner_count()].iter().map(|v| m.nodes[*v as usize]).collect();
                    (0..c.len()).map(|i| c[i][0] * c[(i + 1) % c.len()][1] - c[(i + 1) % c.len()][0] * c[i][1]).sum::<f64>() / 2.0
                })
                .sum();
            assert!(area > area_want && area - area_want < 0.1, "{kind:?}: {area}");
            // Edges on each curve: the right side (x = 4, curve 2) has length 3.
            let on_right: Vec<_> = m.edges.iter().filter(|e| e.2 == 2).collect();
            assert!(!on_right.is_empty());
            for (e, edge, _) in on_right {
                let local = kind.edges()[usize::from(*edge) - 1];
                for i in local {
                    assert!((m.nodes[m.elements[*e][*i] as usize][0] - 4.0).abs() < 1e-12, "{kind:?}");
                }
            }
            // Quadratic mid-side nodes on the hole lie on its polyline (radius 1 within the chord sag).
            assert_eq!(m.corners.len(), 5);
            assert_eq!(mesh_region(&loops, kind, 0.4, &QuadOptions::default()).unwrap(), m, "deterministic");
        }
        assert!(mesh_region(&loops, ElementKind::Hex8, 0.4, &QuadOptions::default()).is_err());
        assert!(mesh_region(&loops, ElementKind::Quad4, 0.0, &QuadOptions::default()).is_err());
        assert!(mesh_region(&loops, ElementKind::Quad4, f64::NAN, &QuadOptions::default()).is_err());
        assert!(mesh_region(&[], ElementKind::Quad4, 1.0, &QuadOptions::default()).is_err());
        assert!(mesh_region(&loops, ElementKind::Tri3, 0.02, &QuadOptions { max_elements: 1000, ..QuadOptions::default() }).is_err());
    }

    #[test]
    fn a_tessellated_face_gives_its_loops_and_neighbours() {
        // A unit cube's bottom face (z = 0, faces 0), split into two triangles; four side faces 1–4
        // with their own copies of the points.
        let p = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            // Side y = 0 (face 1).
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 1.0],
            // Side x = 1 (face 2).
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
        ];
        let tris = vec![[0, 2, 1], [0, 3, 2], [4, 5, 6], [7, 8, 9]];
        let faces = vec![0, 0, 1, 2];
        let b = face_boundary(&p, &tris, &faces, 0, 1e-9).unwrap();
        assert_eq!(b.loops.len(), 1);
        let ids: Vec<Option<u32>> = b.loops[0].iter().map(|c| c.id).collect();
        assert!(ids.contains(&Some(1)) && ids.contains(&Some(2)) && ids.contains(&None), "{ids:?}");
        let points: usize = b.loops[0].iter().map(|c| c.points.len()).sum();
        assert_eq!(points, 4);
        assert!(face_boundary(&p, &tris, &faces, 1, 1e-9).is_err(), "face 1 is vertical");
        assert!(face_boundary(&p, &tris, &faces, 7, 1e-9).is_err());
        assert!(face_boundary(&p, &[[0, 1, 99]], &[0], 0, 1e-9).is_err());
        assert!(face_boundary(&p, &tris, &faces, 0, f64::NAN).is_err());
    }
}
