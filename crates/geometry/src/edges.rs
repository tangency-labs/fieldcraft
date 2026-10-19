//! Persistent names of a body's B-rep edges.
//!
//! An edge is named after the faces it separates: the persistent names of its faces, sorted and
//! joined by `|` (`F2:side:l1|F2:top`), the scheme SolidCraft names edges by. Where several edges
//! separate the same faces (the top and bottom rim of a hole's half-cylinder against two caps
//! cannot repeat, but the two edges between a fillet and a face it splits can), they are numbered
//! `~0`, `~1`, … in order of position along the axis they spread most along, so the names do not
//! depend on the order the kernel lists edges in. Because the names come from the face names, an
//! edge keeps its name through an edit of the design that keeps its faces' names, and the cut
//! edges of a partition are named after the `partition:<name>:a|b` faces.
//!
//! A vertex is named the same way after the faces that meet there (the faces of the edges that
//! end at it): `F1:top|F2:side:l1|F2:side:l2`, numbered `~0`, `~1`, … where several vertices
//! touch the same faces.
//!
//! A shape built in code has no face names; its edges are named after the face indices.

use crate::{Result, Shape, SurfaceMesh};

/// The grid points are matched on, relative to their spread.
const GRID: f64 = 1e-6; // unit-factors: relative to the spread

/// Tessellation tolerance (a fraction of the body's size) of [`Shape::edges`]: the edge polylines
/// are for display, picking and matching edges by position.
pub const EDGE_TOLERANCE: f64 = 1e-3; // unit-factors: relative to the body's size

/// A B-rep edge of a body.
#[derive(Debug, Clone, PartialEq)]
pub struct Edge {
    /// The persistent name (see the module docs).
    pub name: String,
    /// The faces it bounds, by face index (one or two; both the same for a seam).
    pub faces: Vec<u32>,
    /// A polyline along it, start to end.
    pub points: Vec<[f64; 3]>,
}

impl Edge {
    /// The length of its polyline.
    pub fn length(&self) -> f64 {
        self.points.windows(2).map(|w| distance(w[0], w[1])).sum()
    }

    /// The point `fraction` of its length along it from its start (0: the start, 1: the end);
    /// `None` for a fraction outside 0..=1 or an edge without points.
    pub fn point_at(&self, fraction: f64) -> Option<[f64; 3]> {
        if !(0.0..=1.0).contains(&fraction) {
            return None;
        }
        let mut left = fraction * self.length();
        for w in self.points.windows(2) {
            let l = distance(w[0], w[1]);
            if left <= l && l > 0.0 {
                let t = left / l;
                return Some([0, 1, 2].map(|k| w[0][k] + t * (w[1][k] - w[0][k])));
            }
            left -= l;
        }
        self.points.last().copied()
    }
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// A B-rep vertex of a body: where edges end.
#[derive(Debug, Clone, PartialEq)]
pub struct Vertex {
    /// The persistent name (see the module docs).
    pub name: String,
    pub point: [f64; 3],
}

/// Most edge ends looked at for vertices (each against the vertices found so far).
const MAX_EDGE_ENDS: usize = 20_000;

/// The length-weighted centre of a polyline (its first point when it has no length).
fn centre(points: &[[f64; 3]]) -> [f64; 3] {
    let mut sum = [0.0; 3];
    let mut weight = 0.0;
    for w in points.windows(2) {
        let l = distance(w[0], w[1]);
        for (k, s) in sum.iter_mut().enumerate() {
            *s += l * 0.5 * (w[0][k] + w[1][k]);
        }
        weight += l;
    }
    if weight > 0.0 { sum.map(|s| s / weight) } else { points.first().copied().unwrap_or_default() }
}

/// Number the repeated names `~0`, `~1`, … by position: along the axis their centres spread most
/// along, then the other two. Positions are compared on a grid of a millionth of that spread, so
/// a tessellation within a different tolerance numbers them the same.
fn number_repeats(names: &mut [String], centres: &[[f64; 3]]) {
    let mut groups: std::collections::BTreeMap<String, Vec<usize>> = std::collections::BTreeMap::new();
    for (i, n) in names.iter().enumerate() {
        groups.entry(n.clone()).or_default().push(i);
    }
    for (name, mut same) in groups {
        if same.len() < 2 {
            continue;
        }
        let at = |i: usize| centres.get(i).copied().unwrap_or([0.0; 3]);
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for c in same.iter().map(|&i| at(i)) {
            for k in 0..3 {
                lo[k] = lo[k].min(c[k]);
                hi[k] = hi[k].max(c[k]);
            }
        }
        let spread = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
        let axes = if spread[0] >= spread[1] && spread[0] >= spread[2] {
            [0, 1, 2]
        } else if spread[1] >= spread[2] {
            [1, 2, 0]
        } else {
            [2, 0, 1]
        };
        let grid = spread.iter().fold(0.0f64, |m, s| m.max(*s)).max(fieldcraft_model::units::min_geometry_length()) * GRID;
        let key = |i: usize| axes.map(|k| (at(i)[k] / grid).round() as i64);
        same.sort_by_key(|&i| (key(i), i));
        for (k, i) in same.into_iter().enumerate() {
            if let Some(slot) = names.get_mut(i) {
                *slot = format!("{name}~{k}");
            }
        }
    }
}

impl Shape {
    /// The name of face `f`: its persistent name, or its index for a shape without names.
    pub(crate) fn face_label(&self, f: u32) -> String {
        match self.faces.get(f as usize) {
            Some(face) if !face.name.is_empty() => face.name.clone(),
            _ => f.to_string(),
        }
    }

    /// The persistent names of the edges of `mesh` (a tessellation of this shape), by edge index.
    pub fn edge_names(&self, mesh: &SurfaceMesh) -> Vec<String> {
        let mut names: Vec<String> = mesh
            .edge_faces
            .iter()
            .map(|faces| {
                let mut n: Vec<String> = faces.iter().map(|&f| self.face_label(f)).collect();
                n.sort();
                n.join("|")
            })
            .collect();
        names.resize(mesh.edges.len(), String::new());
        let centres: Vec<[f64; 3]> = mesh.edges.iter().map(|e| centre(e)).collect();
        number_repeats(&mut names, &centres);
        names
    }

    /// The vertices of `mesh` (a tessellation of this shape) with their persistent names: the
    /// ends of its edges, those within the tessellation's weld distance being one vertex. A
    /// closed edge with no other edge on it (a full circle) has no vertex. Empty for a shape with
    /// more than [`MAX_EDGE_ENDS`] edge ends.
    pub fn vertex_names(&self, mesh: &SurfaceMesh) -> Vec<Vertex> {
        let ends = mesh.edges.iter().enumerate().flat_map(|(i, e)| [e.first(), e.last()].into_iter().flatten().map(move |p| (i, *p)));
        if mesh.edges.len().saturating_mul(2) > MAX_EDGE_ENDS {
            return Vec::new();
        }
        let close = (2.0 * mesh.weld_distance).max(fieldcraft_model::units::min_geometry_length());
        // Per vertex: its point and the edges that end at it.
        let mut found: Vec<([f64; 3], Vec<usize>)> = Vec::new();
        for (edge, p) in ends {
            match found.iter_mut().find(|(q, _)| distance(*q, p) <= close) {
                Some((_, edges)) => edges.push(edge),
                None => found.push((p, vec![edge])),
            }
        }
        // A closed edge alone (its two ends one vertex, no other edge there) has no corner.
        found.retain(|(_, edges)| edges.windows(2).any(|w| w[0] != w[1]) || edges.len() > 2);
        let mut names: Vec<String> = found
            .iter()
            .map(|(_, edges)| {
                let mut faces: Vec<u32> = edges.iter().filter_map(|e| mesh.edge_faces.get(*e)).flatten().copied().collect();
                faces.sort_unstable();
                faces.dedup();
                let mut n: Vec<String> = faces.into_iter().map(|f| self.face_label(f)).collect();
                n.sort();
                n.join("|")
            })
            .collect();
        let points: Vec<[f64; 3]> = found.iter().map(|(p, _)| *p).collect();
        number_repeats(&mut names, &points);
        names.into_iter().zip(points).map(|(name, point)| Vertex { name, point }).collect()
    }

    /// The body's vertices with their persistent names.
    pub fn vertices(&self) -> Result<Vec<Vertex>> {
        let mesh = self.tessellate((self.size()? * EDGE_TOLERANCE).max(1e-4))?;
        Ok(self.vertex_names(&mesh))
    }

    /// The body's edges with their persistent names, from a tessellation fine enough to place them.
    pub fn edges(&self) -> Result<Vec<Edge>> {
        let mesh = self.tessellate((self.size()? * EDGE_TOLERANCE).max(1e-4))?;
        Ok(self.edges_of(&mesh))
    }

    /// The edges of `mesh` (a tessellation of this shape) with their persistent names.
    pub(crate) fn edges_of(&self, mesh: &SurfaceMesh) -> Vec<Edge> {
        let names = self.edge_names(mesh);
        names
            .into_iter()
            .zip(&mesh.edges)
            .enumerate()
            .map(|(i, (name, points))| Edge { name, faces: mesh.edge_faces.get(i).cloned().unwrap_or_default(), points: points.clone() })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::sample::{self, Plate};
    use crate::{Cut, open_design};

    fn names(edges: &[super::Edge]) -> Vec<String> {
        let mut n: Vec<String> = edges.iter().map(|e| e.name.clone()).collect();
        n.sort();
        n
    }

    #[test]
    fn every_edge_of_the_plate_has_its_own_name_from_its_faces() {
        let d = open_design(&sample::plate_design(&Plate::default()).unwrap()).unwrap();
        let body = &d.bodies[0];
        let edges = body.edges().unwrap();
        // A box with a through hole split in two half-cylinders: 12 box edges, 2 × 2 rim arcs and
        // 2 seams.
        assert_eq!(edges.len(), 18, "{:?}", names(&edges));
        let mut unique = names(&edges);
        unique.dedup();
        assert_eq!(unique.len(), edges.len(), "{unique:?}");
        let faces: Vec<&str> = body.faces().iter().map(|f| f.name.as_str()).collect();
        for e in &edges {
            let base = e.name.split('~').next().unwrap();
            assert!(base.split('|').all(|f| faces.contains(&f)), "{} is not named after the faces {faces:?}", e.name);
            assert!(e.length() > 0.0, "{}", e.name);
        }
    }

    #[test]
    fn edge_names_survive_a_parameter_change_and_a_finer_tessellation() {
        let a = open_design(&sample::plate_design(&Plate::default()).unwrap()).unwrap();
        let b = open_design(&sample::plate_design(&Plate { width: 70.0, thickness: 15.0, ..Plate::default() }).unwrap()).unwrap();
        let (ea, eb) = (a.bodies[0].edges().unwrap(), b.bodies[0].edges().unwrap());
        assert_eq!(names(&ea), names(&eb));
        // The through-thickness edges are as long as the plate is thick, before and after.
        for e in &ea {
            let other = eb.iter().find(|x| x.name == e.name).unwrap();
            if (e.length() - 10.0).abs() < 1e-6 {
                assert!((other.length() - 15.0).abs() < 1e-6, "{}: {}", e.name, other.length());
            }
        }
        let fine = a.bodies[0].tessellate(1e-3).unwrap();
        let mut n = a.bodies[0].edge_names(&fine);
        n.sort();
        assert_eq!(n, names(&ea));
    }

    #[test]
    fn every_corner_of_the_plate_has_its_own_name_from_its_faces() {
        let d = open_design(&sample::plate_design(&Plate::default()).unwrap()).unwrap();
        let body = &d.bodies[0];
        let vertices = body.vertices().unwrap();
        // The box's 8 corners and the ends of the hole's two seams (2 × 2).
        assert_eq!(vertices.len(), 12, "{vertices:?}");
        let mut unique: Vec<&str> = vertices.iter().map(|v| v.name.as_str()).collect();
        unique.sort_unstable();
        unique.dedup();
        assert_eq!(unique.len(), vertices.len(), "{unique:?}");
        // A box corner is where three faces meet.
        let corner = vertices.iter().find(|v| v.point == [0.0, 0.0, 0.0]).unwrap();
        assert_eq!(corner.name.split('|').count(), 3, "{}", corner.name);
        // The names survive an edit that keeps the faces.
        let b = open_design(&sample::plate_design(&Plate { width: 70.0, ..Plate::default() }).unwrap()).unwrap();
        let mut other: Vec<String> = b.bodies[0].vertices().unwrap().into_iter().map(|v| v.name).collect();
        other.sort_unstable();
        assert_eq!(other, unique);
    }

    #[test]
    fn points_along_an_edge() {
        let e = super::Edge { name: "e".into(), faces: vec![], points: vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 2.0, 0.0]] };
        assert_eq!(e.point_at(0.0), Some([0.0, 0.0, 0.0]));
        assert_eq!(e.point_at(0.5), Some([2.0, 0.0, 0.0]));
        assert_eq!(e.point_at(0.75), Some([2.0, 1.0, 0.0]));
        assert_eq!(e.point_at(1.0), Some([2.0, 2.0, 0.0]));
        for bad in [-0.1, 1.5, f64::NAN] {
            assert_eq!(e.point_at(bad), None);
        }
        assert_eq!(super::Edge { points: vec![], ..e }.point_at(0.5), None);
    }

    #[test]
    fn partition_cut_edges_are_named_after_the_cut_faces() {
        let d = open_design(&sample::plate_design(&Plate::default()).unwrap()).unwrap();
        let cut = Cut { name: "Mid".into(), point: [5.0, 15.0, 0.0], normal: [1.0, 0.0, 0.0], bodies: Vec::new() };
        let p = d.partitioned(&[cut]).unwrap();
        for (cell, side) in p.bodies.iter().zip(["a", "b"]) {
            let edges = cell.edges().unwrap();
            let cut_edges: Vec<&str> = edges.iter().map(|e| e.name.as_str()).filter(|n| n.contains(&format!("partition:Mid:{side}"))).collect();
            // The cut face is a rectangle: four edges.
            assert_eq!(cut_edges.len(), 4, "{} {cut_edges:?}", cell.name);
        }
    }
}
