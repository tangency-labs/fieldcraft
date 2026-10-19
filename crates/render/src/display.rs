//! What a model looks like on screen.
//!
//! Finite elements: the faces nobody shares (the skin of the solid elements, and every shell),
//! the element edges on those faces, line elements and nodes. Faces shared by two solid elements
//! are inside the mesh and are left out, so a model of a million elements draws its surface, not
//! its interior. Higher-order elements (tet10, hex20, tri6, quad8) draw through all their nodes:
//! a face is cut into triangles whose vertices are its corners and mid-side nodes, and each
//! edge is drawn as its two halves through its mid-side node, so the faces, the edges, the
//! feature lines and a highlight all share the same vertices and meet exactly at any zoom.
//! Every element primitive keeps its node labels, so results (displacements, contours) can be
//! drawn on it ([`DisplayMesh::displaced`]).
//!
//! Geometry: triangulated boundary surfaces (tessellated bodies, imported surfaces) added with
//! [`DisplayMesh::add_surface`], drawn with the edges between their faces.
//!
//! The output is deterministic: elements come out in label order, edges sorted.

use fieldcraft_model::{ElementKind, Model};

use crate::pick::{Hit, first_hit};
use crate::{Aabb3, Vec3};

/// One drawn triangle of an element face.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DisplayTri {
    pub p: [Vec3; 3],
    /// The node labels at the corners.
    pub nodes: [u32; 3],
    /// The face's unit normal, pointing out of its element (a shell's own orientation).
    pub normal: Vec3,
    pub kind: ElementKind,
    /// The element it is a face of, and the face's number (from 1, as
    /// [`ElementKind::faces`] numbers them; a shell's is 1, its top side).
    pub element: u32,
    pub face: u8,
}

impl DisplayTri {
    /// The element face seen looking along `dir`: a shell's top side (1) from in front and its
    /// bottom side (2) from behind; a solid's face as it is.
    pub fn seen_face(&self, dir: Vec3) -> u8 {
        if self.kind.dimension() == 2 && self.normal.dot(dir) > 0.0 { 2 } else { self.face }
    }
}

/// A drawn segment of an element edge: the whole edge of a linear element, a half (corner to
/// mid-side node) of a quadratic one, so edges bend with the faces they bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DisplayEdge {
    pub p: [Vec3; 2],
    pub nodes: [u32; 2],
    /// The corner nodes of the element edge it is (part of), lower first.
    pub corners: [u32; 2],
    /// The element edge lies on a geometry edge: a feature line of the mesh.
    pub feature: bool,
}

/// A line element.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DisplayLine {
    pub a: Vec3,
    pub b: Vec3,
    pub nodes: [u32; 2],
    pub kind: ElementKind,
    pub element: u32,
}

/// A node.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DisplayNode {
    pub label: u32,
    pub p: Vec3,
}

/// A triangle of a geometry surface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceTri {
    pub p: [Vec3; 3],
    pub normal: Vec3,
    /// Which surface (body) it belongs to: the order [`DisplayMesh::add_surface`] got them in.
    pub body: u32,
    /// The face of that surface it lies on (its `tri_face` tag; `u32::MAX` when it has none).
    pub face: u32,
}

/// A segment of a geometry surface's B-rep edge (a named design edge), for picking and marking.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceCurve {
    pub a: Vec3,
    pub b: Vec3,
    /// The surface (body) it belongs to, as [`SurfaceTri::body`].
    pub body: u32,
    /// The edge's index in that body's tessellation.
    pub edge: u32,
}

/// Display data for a model and its geometry.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DisplayMesh {
    /// Visible element faces as triangles (quads are split in two).
    pub tris: Vec<DisplayTri>,
    /// Element edges on the visible faces, each once.
    pub edges: Vec<DisplayEdge>,
    /// Line elements (rods, beams).
    pub lines: Vec<DisplayLine>,
    /// Every node, in label order.
    pub nodes: Vec<DisplayNode>,
    /// Geometry surface triangles.
    pub surface_tris: Vec<SurfaceTri>,
    /// Geometry edges: where two faces meet, and open boundaries.
    pub surface_edges: Vec<[Vec3; 2]>,
    /// The B-rep edges of the geometry surfaces that have them, segment by segment.
    pub curves: Vec<SurfaceCurve>,
    /// Geometry surfaces (bodies) added.
    pub bodies: u32,
    pub bounds: Aabb3,
}

/// The faces drawn of an element, as positions in its node list ([`ElementKind::faces`]: the
/// corners counter-clockwise seen from outside, then the mid-side nodes of quadratic kinds): a
/// solid's faces in their numbering, a shell's top side (face 1).
fn drawn_faces(kind: ElementKind) -> &'static [&'static [usize]] {
    let all = kind.faces();
    if kind.dimension() == 2 { all.get(..1).unwrap_or(&[]) } else { all }
}

/// The corners of a face listed with `n` nodes: a quadratic face (tri6, quad8) has as many
/// mid-side nodes as corners.
fn face_corner_count(n: usize) -> usize {
    if n >= 6 { n / 2 } else { n }
}

/// Triangles (positions in the face's boundary ring) that fill a face whose boundary ring has
/// `n` points. A quadratic face's ring alternates corner and mid-side node, and every triangle's
/// vertices are nodes: one triangle cut off at each corner, and the polygon of the mid-side
/// nodes fanned, so the drawn face passes through every node and its sides are the drawn edges.
fn face_triangles(n: usize, quadratic: bool) -> Vec<[usize; 3]> {
    if !quadratic {
        return (1..n.saturating_sub(1)).map(|k| [0, k, k + 1]).collect();
    }
    let mut out: Vec<[usize; 3]> = (0..n).step_by(2).map(|c| [c, c + 1, (c + n - 1) % n]).collect();
    let mids: Vec<usize> = (1..n).step_by(2).collect();
    for k in 1..mids.len().saturating_sub(1) {
        if let (Some(&a), Some(&b), Some(&c)) = (mids.first(), mids.get(k), mids.get(k + 1)) {
            out.push([a, b, c]);
        }
    }
    out
}

/// A face's identity regardless of winding and starting corner: its sorted corner labels, a
/// triangle padded with `u32::MAX`.
fn face_key(e: &[u32], face: &[usize]) -> Option<[u32; 4]> {
    let mut k = [u32::MAX; 4];
    for (slot, &i) in k.iter_mut().zip(face.iter().take(face_corner_count(face.len()))) {
        *slot = *e.get(i)?;
    }
    k.sort_unstable();
    Some(k)
}

/// Newell's normal of a polygon (robust for non-planar quads); `None` when degenerate.
pub(crate) fn polygon_normal(p: &[Vec3]) -> Option<Vec3> {
    let mut n = Vec3::ZERO;
    for (i, a) in p.iter().enumerate() {
        let b = p.get((i + 1) % p.len()).copied().unwrap_or(*a);
        n = n + Vec3::new((a.y - b.y) * (a.z + b.z), (a.z - b.z) * (a.x + b.x), (a.x - b.x) * (a.y + b.y));
    }
    n.normalized()
}

pub(crate) fn centre(p: &[Vec3]) -> Vec3 {
    p.iter().fold(Vec3::ZERO, |s, q| s + *q) * (1.0 / p.len().max(1) as f64)
}

impl DisplayMesh {
    pub fn from_model(model: &Model) -> DisplayMesh {
        DisplayMesh::from_model_hiding(model, &std::collections::BTreeSet::new())
    }

    /// [`DisplayMesh::from_model`] without the `hidden` elements (a hidden body's): their faces
    /// are not drawn, a face they share with a shown element is drawn as its boundary, and nodes
    /// only they use are left out (so nothing picks or frames them).
    pub fn from_model_hiding(model: &Model, hidden: &std::collections::BTreeSet<u32>) -> DisplayMesh {
        let mut out = DisplayMesh::default();
        let pos = |label: u32| model.node(label).map(|n| Vec3::from(n.xyz));
        let shown = |label: &u32| !hidden.contains(label);
        // With elements hidden, the nodes only they use go too.
        let hidden_nodes: std::collections::BTreeSet<u32> = if hidden.is_empty() {
            Default::default()
        } else {
            let used: std::collections::BTreeSet<u32> =
                model.elements().iter().filter(|(l, _)| shown(l)).flat_map(|(_, e)| e.nodes.iter().copied()).collect();
            hidden.iter().filter_map(|l| model.element(*l)).flat_map(|e| e.nodes.iter().copied()).filter(|n| !used.contains(n)).collect()
        };
        for (&label, n) in model.nodes().iter().filter(|(l, _)| !hidden_nodes.contains(l)) {
            let p = Vec3::from(n.xyz);
            out.nodes.push(DisplayNode { label, p });
            out.bounds.add(p);
        }
        // Faces two solid elements share are inside the mesh: sort every solid face's key and
        // keep the ones that occur more than once.
        let mut keys: Vec<[u32; 4]> = model
            .elements()
            .iter()
            .filter(|(l, _)| shown(l))
            .map(|(_, e)| e)
            .filter(|e| e.kind.dimension() == 3)
            .flat_map(|e| drawn_faces(e.kind).iter().filter_map(|f| face_key(&e.nodes, f)))
            .collect();
        keys.sort_unstable();
        let mut inner: Vec<[u32; 4]> = keys.windows(2).filter(|w| w.first() == w.get(1)).filter_map(|w| w.first().copied()).collect();
        inner.dedup();
        drop(keys);
        // Element edges on geometry edges (by their corner nodes): the mesh's feature lines.
        let classification = model.classification();
        let mut features: Vec<(u32, u32)> = Vec::new();
        for (&element, tags) in classification.elements().iter().filter(|(l, _)| shown(l)) {
            let Some(e) = model.element(element) else { continue };
            for (&edge, &tag) in &tags.edges {
                if classification.entity(tag).is_some_and(|g| g.entity.dimension() == 1)
                    && let Some(local) = e.kind.edges().get(usize::from(edge).wrapping_sub(1))
                    && let (Some(&a), Some(&b)) = (local.first().and_then(|&i| e.nodes.get(i)), local.get(1).and_then(|&i| e.nodes.get(i)))
                {
                    features.push((a.min(b), a.max(b)));
                }
            }
        }
        features.sort_unstable();
        features.dedup();
        // (segment's nodes, the element edge's corners), each sorted.
        let mut edges: Vec<((u32, u32), (u32, u32))> = Vec::new();
        for (&label, e) in model.elements().iter().filter(|(l, _)| shown(l)) {
            if e.kind.dimension() == 1 {
                if let (Some(&na), Some(&nb)) = (e.nodes.first(), e.nodes.get(1))
                    && let (Some(a), Some(b)) = (pos(na), pos(nb))
                {
                    out.lines.push(DisplayLine { a, b, nodes: [na, nb], kind: e.kind, element: label });
                }
                continue;
            }
            let Some(at) = e.nodes.iter().map(|&n| pos(n)).collect::<Option<Vec<Vec3>>>() else { continue };
            if at.len() != e.kind.node_count() {
                continue;
            }
            let solid = e.kind.dimension() == 3;
            let element_centre = centre(at.get(..e.kind.corner_count()).unwrap_or(&at));
            for (number, f) in drawn_faces(e.kind).iter().enumerate() {
                let face = u8::try_from(number + 1).unwrap_or(u8::MAX);
                if solid && face_key(&e.nodes, f).is_some_and(|k| inner.binary_search(&k).is_ok()) {
                    continue;
                }
                // The face's boundary: corner, mid-side node, corner, … (corners only when linear).
                let corners = face_corner_count(f.len());
                let ring: Vec<usize> = (0..corners).flat_map(|j| std::iter::once(j).chain((f.len() > corners).then_some(corners + j))).collect();
                let Some(ring) =
                    ring.iter().map(|&k| f.get(k).and_then(|&i| Some((*at.get(i)?, *e.nodes.get(i)?)))).collect::<Option<Vec<(Vec3, u32)>>>()
                else {
                    continue;
                };
                let p: Vec<Vec3> = ring.iter().map(|r| r.0).collect();
                let Some(mut normal) = polygon_normal(&p) else { continue };
                if solid && normal.dot(centre(&p) - element_centre) < 0.0 {
                    normal = -normal;
                }
                for [a, b, c] in face_triangles(ring.len(), corners < ring.len()) {
                    if let (Some(&(pa, la)), Some(&(pb, lb)), Some(&(pc, lc))) = (ring.get(a), ring.get(b), ring.get(c)) {
                        out.tris.push(DisplayTri { p: [pa, pb, pc], nodes: [la, lb, lc], normal, kind: e.kind, element: label, face });
                    }
                }
                let step = ring.len() / corners.max(1);
                for (k, &(_, a)) in ring.iter().enumerate() {
                    let b = ring.get((k + 1) % ring.len()).map_or(a, |r| r.1);
                    let (c0, c1) = (ring.get(k - k % step).map_or(a, |r| r.1), ring.get((k - k % step + step) % ring.len()).map_or(b, |r| r.1));
                    edges.push(((a.min(b), a.max(b)), (c0.min(c1), c0.max(c1))));
                }
            }
        }
        edges.sort_unstable();
        edges.dedup();
        out.edges = edges
            .into_iter()
            .filter_map(|((a, b), corners)| {
                let feature = features.binary_search(&corners).is_ok();
                Some(DisplayEdge { p: [pos(a)?, pos(b)?], nodes: [a, b], corners: [corners.0, corners.1], feature })
            })
            .collect();
        out
    }

    /// Add a triangulated geometry surface: `positions`, counter-clockwise `triangles` seen from
    /// outside, and the face each triangle lies on. Edges are drawn where triangles of different
    /// faces meet and along open boundaries; degenerate triangles and out-of-range indices are
    /// skipped.
    pub fn add_surface(&mut self, positions: &[[f64; 3]], triangles: &[[u32; 3]], tri_face: &[u32]) {
        let pos = |i: u32| positions.get(i as usize).map(|p| Vec3::from(*p)).filter(|p| p.is_finite());
        // Tessellations may repeat a point for each triangle or face: copies with identical
        // coordinates are one vertex (the lowest index of the copies).
        let mut order: Vec<u32> = (0..u32::try_from(positions.len()).unwrap_or(u32::MAX)).collect();
        let coords = |i: &u32| positions.get(*i as usize).copied().unwrap_or([f64::NAN; 3]);
        order.sort_by(|a, b| {
            let (p, q) = (coords(a), coords(b));
            p[0].total_cmp(&q[0]).then(p[1].total_cmp(&q[1])).then(p[2].total_cmp(&q[2])).then(a.cmp(b))
        });
        let mut vertex: Vec<u32> = (0..order.len()).map(|i| i as u32).collect();
        for run in order.chunk_by(|a, b| coords(a).map(f64::to_bits) == coords(b).map(f64::to_bits)) {
            let first = run.first().copied().unwrap_or(0);
            for &i in run {
                if let Some(v) = vertex.get_mut(i as usize) {
                    *v = first;
                }
            }
        }
        let vertex = |i: u32| vertex.get(i as usize).copied().unwrap_or(i);
        // (low vertex, high vertex, face) per triangle side; a side alone, or sides on
        // different faces, make an edge to draw.
        let mut sides: Vec<(u32, u32, u32)> = Vec::with_capacity(triangles.len().saturating_mul(3));
        for (t, &[i0, i1, i2]) in triangles.iter().enumerate() {
            let (Some(a), Some(b), Some(c)) = (pos(i0), pos(i1), pos(i2)) else { continue };
            let Some(normal) = (b - a).cross(c - a).normalized() else { continue };
            for p in [a, b, c] {
                self.bounds.add(p);
            }
            let face = tri_face.get(t).copied().unwrap_or(u32::MAX);
            self.surface_tris.push(SurfaceTri { p: [a, b, c], normal, body: self.bodies, face });
            let [v0, v1, v2] = [vertex(i0), vertex(i1), vertex(i2)];
            for (u, v) in [(v0, v1), (v1, v2), (v2, v0)] {
                sides.push((u.min(v), u.max(v), face));
            }
        }
        sides.sort_unstable();
        for run in sides.chunk_by(|a, b| a.0 == b.0 && a.1 == b.1) {
            let Some(&(u, v, face)) = run.first() else { continue };
            if (run.len() == 1 || run.iter().any(|s| s.2 != face))
                && let (Some(a), Some(b)) = (pos(u), pos(v))
            {
                self.surface_edges.push([a, b]);
            }
        }
        self.bodies = self.bodies.saturating_add(1);
    }

    /// Add the B-rep edges (polylines, by edge index) of the surface [`DisplayMesh::add_surface`]
    /// adds next, so they are picked as that body's edges; non-finite points are skipped.
    pub fn add_curves(&mut self, edges: &[Vec<[f64; 3]>]) {
        let body = self.bodies;
        for (e, points) in edges.iter().enumerate() {
            let Ok(edge) = u32::try_from(e) else { break };
            for w in points.windows(2) {
                let (a, b) = (Vec3::from(w[0]), Vec3::from(w[1]));
                if a.is_finite() && b.is_finite() {
                    self.curves.push(SurfaceCurve { a, b, body, edge });
                }
            }
        }
    }

    /// The colour of a geometry triangle: `single` when there is one body, else body `k`'s
    /// colour from `palette` (cycling), so the bodies of a design, such as the cells of a
    /// partitioned part, tell apart.
    pub fn body_color<C: Copy>(&self, tri: &SurfaceTri, single: C, palette: &[C]) -> C {
        if self.bodies <= 1 || palette.is_empty() {
            return single;
        }
        palette.get(tri.body as usize % palette.len()).copied().unwrap_or(single)
    }

    /// Where node `label` is drawn, if it is (nodes are kept in label order).
    pub fn node_position(&self, label: u32) -> Option<Vec3> {
        self.nodes.binary_search_by_key(&label, |n| n.label).ok().and_then(|i| self.nodes.get(i)).map(|n| n.p)
    }

    /// The element face triangle ([`DisplayMesh::tris`]) a ray (`origin`, `dir`) hits first
    /// ([`pick::first_hit`]).
    pub fn hit_element_face(&self, origin: Vec3, dir: Vec3) -> Option<Hit> {
        first_hit(origin, dir, self.tris.iter().map(|t| t.p))
    }

    /// The geometry surface triangle ([`DisplayMesh::surface_tris`]) a ray hits first.
    pub fn hit_surface(&self, origin: Vec3, dir: Vec3) -> Option<Hit> {
        first_hit(origin, dir, self.surface_tris.iter().map(|t| t.p))
    }

    /// The node nearest to where a ray (`origin`, `dir`) first hits an element face ahead of
    /// its origin: the hit face's closest corner, and the hit point. Line elements and geometry
    /// are not hit.
    pub fn pick(&self, origin: Vec3, dir: Vec3) -> Option<(u32, Vec3)> {
        let h = self.hit_element_face(origin, dir)?;
        let t = self.tris.get(h.index)?;
        let hit = origin + dir * h.t;
        let (i, _) = t.p.iter().enumerate().min_by(|a, b| a.1.dist(hit).total_cmp(&b.1.dist(hit)))?;
        Some((*t.nodes.get(i)?, hit))
    }

    /// The elements moved by `offset(label)` (no offset for a node it has none for): faces,
    /// edges, line elements and nodes move, face normals follow the deformation. Geometry
    /// surfaces are left out (they do not deform).
    pub fn displaced(&self, offset: impl Fn(u32) -> Option<Vec3>) -> DisplayMesh {
        let at = |p: Vec3, label: u32| p + offset(label).filter(|d| d.is_finite()).unwrap_or(Vec3::ZERO);
        let tris = self
            .tris
            .iter()
            .map(|t| {
                let p = [at(t.p[0], t.nodes[0]), at(t.p[1], t.nodes[1]), at(t.p[2], t.nodes[2])];
                // Unmoved, it is as it was (a normal worked out again could differ in its last bit).
                if p == t.p {
                    return *t;
                }
                // The triangle's own normal, on the side its face's normal was.
                let before = (t.p[1] - t.p[0]).cross(t.p[2] - t.p[0]);
                let side = if before.dot(t.normal) < 0.0 { -1.0 } else { 1.0 };
                let normal = ((p[1] - p[0]).cross(p[2] - p[0]) * side).normalized().unwrap_or(t.normal);
                DisplayTri { p, normal, ..*t }
            })
            .collect();
        let mut out = DisplayMesh {
            tris,
            edges: self.edges.iter().map(|e| DisplayEdge { p: [at(e.p[0], e.nodes[0]), at(e.p[1], e.nodes[1])], ..*e }).collect(),
            lines: self.lines.iter().map(|l| DisplayLine { a: at(l.a, l.nodes[0]), b: at(l.b, l.nodes[1]), ..*l }).collect(),
            nodes: self.nodes.iter().map(|n| DisplayNode { label: n.label, p: at(n.p, n.label) }).collect(),
            ..Default::default()
        };
        for n in &out.nodes {
            out.bounds.add(n.p);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `nx × ny × nz` unit hex8 cells.
    fn block(nx: u32, ny: u32, nz: u32) -> Model {
        let mut m = Model::default();
        let id = |i: u32, j: u32, k: u32| 1 + i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
        for k in 0..=nz {
            for j in 0..=ny {
                for i in 0..=nx {
                    m.set_node(id(i, j, k), [i as f64, j as f64, k as f64]).unwrap();
                }
            }
        }
        let mut label = 1;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let n = vec![
                        id(i, j, k),
                        id(i + 1, j, k),
                        id(i + 1, j + 1, k),
                        id(i, j + 1, k),
                        id(i, j, k + 1),
                        id(i + 1, j, k + 1),
                        id(i + 1, j + 1, k + 1),
                        id(i, j + 1, k + 1),
                    ];
                    m.set_element(label, ElementKind::Hex8, n).unwrap();
                    label += 1;
                }
            }
        }
        m
    }

    #[test]
    fn a_block_shows_its_skin_only() {
        let d = DisplayMesh::from_model(&block(3, 2, 2));
        // Skin quads: 2·(3·2 + 3·2 + 2·2) = 32, two triangles each.
        assert_eq!(d.tris.len(), 64);
        assert_eq!(d.nodes.len(), 4 * 3 * 3);
        assert_eq!(d.bounds.max, Vec3::new(3.0, 2.0, 2.0));
        // Every normal points out of the block.
        let c = d.bounds.center();
        for t in &d.tris {
            let fc = (t.p[0] + t.p[1] + t.p[2]) * (1.0 / 3.0);
            assert!(t.normal.dot(fc - c) > 0.0, "{t:?}");
        }
        // Euler on the skin: V - E + F = 2 (two nodes are inside).
        let skin_nodes = 4 * 3 * 3 - 2;
        assert_eq!(skin_nodes as i64 - d.edges.len() as i64 + 32, 2);
    }

    /// Hidden elements (a hidden body's) are not drawn; the face they shared with a shown one
    /// is its boundary now, and nodes only they used are gone (#365).
    #[test]
    fn hidden_elements_are_left_out() {
        let m = block(2, 1, 1);
        let d = DisplayMesh::from_model_hiding(&m, &[2].into_iter().collect());
        assert_eq!(d.tris.len(), 12, "one cube's six faces");
        assert_eq!(d.nodes.len(), 8);
        assert_eq!(d.bounds.max, Vec3::new(1.0, 1.0, 1.0));
    }

    #[test]
    fn shells_lines_and_tets_draw() {
        let mut m = Model::default();
        for (i, p) in [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.], [1., 1., 0.]].into_iter().enumerate() {
            m.set_node(i as u32 + 1, p).unwrap();
        }
        m.set_element(1, ElementKind::Tet4, vec![1, 2, 3, 4]).unwrap();
        m.set_element(2, ElementKind::Tri3, vec![2, 5, 3]).unwrap();
        m.set_element(3, ElementKind::Beam2, vec![4, 5]).unwrap();
        let d = DisplayMesh::from_model(&m);
        assert_eq!(d.tris.len(), 5);
        assert_eq!(d.lines.len(), 1);
        assert_eq!(d.edges.len(), 6 + 2);
        let tet_normals: Vec<Vec3> = d.tris.iter().filter(|t| t.kind == ElementKind::Tet4).map(|t| t.normal).collect();
        assert!(tet_normals.contains(&Vec3::new(0.0, 0.0, -1.0)));
        assert_eq!(d.lines[0].nodes, [4, 5]);
        assert!(d.tris.iter().all(|t| t.nodes.iter().zip(t.p).all(|(&l, p)| Vec3::from(m.node(l).unwrap().xyz) == p)));
    }

    #[test]
    fn surfaces_draw_their_face_boundaries() {
        // A unit square in two faces (one triangle each), plus a degenerate and an out-of-range
        // triangle, which are skipped.
        let p = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]];
        let mut d = DisplayMesh::default();
        // Its B-rep edges: the diagonal (two segments) and a side; a non-finite point is skipped.
        d.add_curves(&[vec![p[0], [0.5, 0.5, 0.0], p[2]], vec![p[0], p[1]], vec![[f64::NAN; 3], p[1]]]);
        d.add_surface(&p, &[[0, 1, 2], [0, 2, 3], [0, 0, 1], [0, 1, 9]], &[0, 1, 1, 1]);
        assert_eq!(d.curves.iter().map(|c| (c.body, c.edge)).collect::<Vec<_>>(), vec![(0, 0), (0, 0), (0, 1)]);
        assert_eq!(d.surface_tris.len(), 2);
        // Four open sides and the diagonal between the two faces.
        assert_eq!(d.surface_edges.len(), 5);
        let mut one_face = DisplayMesh::default();
        one_face.add_surface(&p, &[[0, 1, 2], [0, 2, 3]], &[0, 0]);
        assert_eq!(one_face.surface_edges.len(), 4);
        assert_eq!(one_face.bounds.max, Vec3::new(1.0, 1.0, 0.0));
        // A triangle soup (each triangle with its own copies of the corners) is welded first.
        let soup = [p[0], p[1], p[2], p[0], p[2], p[3]];
        let mut welded = DisplayMesh::default();
        welded.add_surface(&soup, &[[0, 1, 2], [3, 4, 5]], &[0, 0]);
        assert_eq!(welded.surface_edges.len(), 4);
        let mut hostile = DisplayMesh::default();
        let mut two = DisplayMesh::default();
        two.add_surface(&p, &[[0, 1, 2]], &[0]);
        two.add_surface(&p, &[[0, 2, 3]], &[0]);
        assert_eq!(two.bodies, 2);
        let colors: Vec<u8> = two.surface_tris.iter().map(|t| two.body_color(t, 0u8, &[7, 8])).collect();
        assert_eq!(colors, vec![7, 8]);
        assert_eq!(one_face.body_color(&one_face.surface_tris[0], 0u8, &[7, 8]), 0);
        hostile.add_surface(&[[f64::NAN, 0., 0.], [1., 0., 0.], [0., 1., 0.]], &[[0, 1, 2], [u32::MAX, 0, 1]], &[]);
        assert!(hostile.surface_tris.is_empty() && hostile.surface_edges.is_empty());
    }

    #[test]
    fn displacing_moves_elements_by_node_and_keeps_normals_outward() {
        let d = DisplayMesh::from_model(&block(2, 1, 1));
        // Stretch along x by 50 %: node labels carry the offset to every primitive.
        let m = block(2, 1, 1);
        let moved = d.displaced(|label| m.node(label).map(|n| Vec3::new(n.xyz[0] * 0.5, 0.0, 0.0)));
        assert_eq!(moved.tris.len(), d.tris.len());
        assert_eq!(moved.bounds.max, Vec3::new(3.0, 1.0, 1.0));
        assert!(moved.edges.iter().all(|e| e.p[0].x == 1.5 * m.node(e.nodes[0]).unwrap().xyz[0]));
        let c = moved.bounds.center();
        for t in &moved.tris {
            let fc = (t.p[0] + t.p[1] + t.p[2]) * (1.0 / 3.0);
            assert!(t.normal.dot(fc - c) > 0.0, "{t:?}");
        }
        // No offset, or a non-finite one, leaves a node where it is.
        let same = d.displaced(|label| (label == 1).then_some(Vec3::new(f64::NAN, 0.0, 0.0)));
        assert_eq!(same.tris, d.tris);
    }

    #[test]
    fn picking_finds_the_nearest_corner_of_the_face_hit_first() {
        let m = block(2, 1, 1);
        let d = DisplayMesh::from_model(&m);
        // From above, aimed near the corner (2, 1, 1): the top face is hit, not the bottom.
        let (label, hit) = d.pick(Vec3::new(1.9, 0.9, 10.0), Vec3::new(0.0, 0.0, -1.0)).unwrap();
        assert_eq!(m.node(label).unwrap().xyz, [2.0, 1.0, 1.0]);
        assert!((hit.z - 1.0).abs() < 1e-12);
        // From below the bottom face wins.
        let (label, _) = d.pick(Vec3::new(0.1, 0.1, -10.0), Vec3::new(0.0, 0.0, 1.0)).unwrap();
        assert_eq!(m.node(label).unwrap().xyz, [0.0, 0.0, 0.0]);
        assert_eq!(d.pick(Vec3::new(5.0, 5.0, 10.0), Vec3::new(0.0, 0.0, -1.0)), None);
        assert_eq!(d.pick(Vec3::new(f64::NAN, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0)), None);
    }

    /// A solid's drawn faces are numbered as the model numbers them, so a picked triangle names
    /// the element face a set or load refers to.
    #[test]
    fn drawn_faces_carry_the_models_face_numbers() {
        let m = block(2, 1, 1);
        let d = DisplayMesh::from_model(&m);
        for t in &d.tris {
            let face = m.face_nodes(fieldcraft_model::ElementFace { element: t.element, face: t.face }).unwrap();
            assert!(t.nodes.iter().all(|n| face.contains(n)), "{t:?} {face:?}");
        }
        // The ray from +x hits element 2's face 4 (corners 1, 2, 6, 5: x = 2).
        let h = d.hit_element_face(Vec3::new(10.0, 0.3, 0.6), Vec3::new(-1.0, 0.0, 0.0)).unwrap();
        let t = d.tris[h.index];
        assert_eq!((t.element, t.face, h.t), (2, 4, 8.0));
    }

    /// `block(2, 1, 1)` as hex20: a mid-side node on every edge, the ones on the top face
    /// raised by 0.1 (a curved top, as a mesher puts them on a curved surface), and element 1's
    /// edge 5 (corners 5 → 6, along the top) on a geometry edge.
    fn curved_hex20_block() -> Model {
        let linear = block(2, 1, 1);
        let mut m = Model::default();
        for (&l, n) in linear.nodes() {
            m.set_node(l, n.xyz).unwrap();
        }
        let mut mids = std::collections::BTreeMap::new();
        let mut next = 100;
        for (&label, e) in linear.elements() {
            let mut nodes = e.nodes.clone();
            for edge in ElementKind::Hex20.edges() {
                let (a, b) = (e.nodes[edge[0]], e.nodes[edge[1]]);
                let id = *mids.entry((a.min(b), a.max(b))).or_insert_with(|| {
                    let (p, q) = (linear.node(a).unwrap().xyz, linear.node(b).unwrap().xyz);
                    let mut mid = [0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1]), 0.5 * (p[2] + q[2])];
                    if p[2] == 1.0 && q[2] == 1.0 {
                        mid[2] += 0.1;
                    }
                    next += 1;
                    m.set_node(next, mid).unwrap();
                    next
                });
                nodes.push(id);
            }
            m.set_element(label, ElementKind::Hex20, nodes).unwrap();
        }
        let edge = fieldcraft_model::GeometryRef::new("B", fieldcraft_model::Entity::Edge("top".into())).unwrap();
        let tag = m.geometry_tag(&edge).unwrap();
        m.tag_element_edge(1, 5, tag).unwrap();
        m
    }

    /// A curved hex20 face is drawn through its mid-side nodes, and its faces, edges and feature
    /// lines share their vertices exactly (#260: lines met the faces nowhere when zoomed in).
    #[test]
    fn quadratic_faces_and_edges_share_their_vertices() {
        let m = curved_hex20_block();
        let d = DisplayMesh::from_model(&m);
        let at = |l: u32| Vec3::from(m.node(l).unwrap().xyz);
        // Every vertex drawn is a node, exactly where the node is.
        for t in &d.tris {
            assert!(t.nodes.iter().zip(t.p).all(|(&l, p)| at(l) == p), "{t:?}");
        }
        for e in &d.edges {
            assert!(e.nodes.iter().zip(e.p).all(|(&l, p)| at(l) == p), "{e:?}");
        }
        // 10 skin faces, 6 triangles each; the 20 element edges in two halves each.
        assert_eq!(d.tris.len(), 60);
        assert_eq!(d.edges.len(), 40);
        // The drawn edges are exactly the sides of the faces' outlines.
        let mut sides: Vec<((u32, u8), [u32; 2])> = d
            .tris
            .iter()
            .flat_map(|t| {
                [0, 1, 2].map(|k| {
                    ((t.element, t.face), {
                        let mut s = [t.nodes[k], t.nodes[(k + 1) % 3]];
                        s.sort_unstable();
                        s
                    })
                })
            })
            .collect();
        sides.sort_unstable();
        let mut outline: Vec<[u32; 2]> = sides.chunk_by(|a, b| a == b).filter(|r| r.len() == 1).map(|r| r[0].1).collect();
        outline.sort_unstable();
        outline.dedup();
        let mut drawn: Vec<[u32; 2]> = d.edges.iter().map(|e| e.nodes).collect();
        drawn.sort_unstable();
        assert_eq!(outline, drawn);
        // A top edge bows up through its mid-side node.
        assert!(d.edges.iter().any(|e| e.p.iter().any(|p| (p.z - 1.1).abs() < 1e-12)));
        // The edge on the geometry edge is a feature line, both halves; the rest are not.
        let feature: Vec<[u32; 2]> = d.edges.iter().filter(|e| e.feature).map(|e| e.corners).collect();
        let e1 = &m.element(1).unwrap().nodes;
        let top = [e1[4].min(e1[5]), e1[4].max(e1[5])];
        assert_eq!(feature, vec![top, top]);
    }

    #[test]
    fn empty_and_deterministic() {
        assert_eq!(DisplayMesh::from_model(&Model::default()), DisplayMesh::default());
        let m = block(2, 2, 1);
        assert_eq!(DisplayMesh::from_model(&m), DisplayMesh::from_model(&m.clone()));
    }
}
