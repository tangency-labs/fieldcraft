//! Virtual topology: faces merged for meshing only.
//!
//! SolidCraft's kernel cannot heal a B-rep (no face merge or split), and a mesher should not have
//! to: a fillet the kernel left in pieces, or a face an imprint split, only needs to *mesh* as one
//! face. A [`Merge`] makes adjacent faces of a body one virtual face in what the meshers take from
//! the body ([`Shape::tessellate`], [`Shape::faces`], [`Shape::face_surfaces`], edge and vertex
//! names): the triangles of its B-rep faces carry one tag, the edges between them are gone (an
//! edge they leave split in pieces along a neighbour is joined into one), and its surface is the
//! one its B-rep faces share (none when they lie on different surfaces: the projection mesher and
//! node relaxation then keep to each B-rep piece's own surface, [`Shape::tessellate_pieces`],
//! the other meshers to the tessellation) — except that an absorb keeps the surface of the face
//! it was absorbed into, however the absorbed one lies (#546): the cleaning must not cost the
//! neighbour its classification. The B-rep itself is unchanged.
//!
//! Merges live in the FieldCraft model, by persistent face names, and are applied again after the
//! partitions whenever the design is opened, so an edit of the design followed by a re-open
//! replays them. On a body the partitions cut (a mesh plan's symmetry and corner cuts), a merge
//! applies to each cell: the pieces of its faces that touch there are one face of the cell. A
//! merge whose faces are gone or no longer touch is listed in [`Design::errors`]
//! (feature `merge <name>`), as a feature that failed to rebuild; the others still apply.
//!
//! The merged face is named `merge:<name>`. A face name in a merge selects the face of that name,
//! or the pieces a partition cut it into (`<name>#0`, `<name>#1`, …), as sets do.

use std::collections::BTreeSet;

use crate::edges::EDGE_TOLERANCE;
use crate::{Design, Error, Face, FaceSurface, FeatureError, Result, Shape, SurfaceMesh};

/// Surfaces this close, relative to the body's size, are one.
const SAME_SURFACE: f64 = 1e-6; // unit-factors: relative to the body's size

/// Most merges a design takes.
pub const MAX_MERGES: usize = 1024;
/// Most faces one merge names.
pub const MAX_MERGE_FACES: usize = 1024;
/// Most edges of one pair of faces joined end to end (more are left as they are).
const MAX_JOINED: usize = 256;

/// Faces of one body made one face for meshing.
#[derive(Debug, Clone, PartialEq)]
pub struct Merge {
    /// Names the merged face (`merge:<name>`).
    pub name: String,
    pub body: String,
    /// The faces, by persistent name (a name also selects the pieces a partition cut it into).
    pub faces: Vec<String>,
    /// Absorb: the merged face keeps the name of this one of `faces` instead of `merge:<name>`
    /// (a tiny face absorbed into a neighbour leaves the neighbour, by its own name).
    pub into: Option<String>,
}

/// The name of the face merge `name` makes.
pub fn merged_face_name(name: &str) -> String {
    format!("merge:{name}")
}

/// Whether a feature error comes from the model's partitions or merges (applied again from the
/// design as opened) rather than from the design's own timeline.
pub(crate) fn derived(feature: &str) -> bool {
    feature.starts_with("partition ") || feature.starts_with("merge ") || feature.starts_with("split ")
}

impl Design {
    /// The design with `merges` applied in order, on top of its partitions. A merge on a body
    /// the partitions cut applies to each of its cells: the faces it names that a cell has
    /// (their pieces included) and that touch there are made one face of that cell ([`cell_merges`]).
    /// A merge that fails on this design (a body or face that is not there, faces that do not
    /// touch, in any cell) is listed in [`Design::errors`] (feature `merge <name>`) and changes
    /// nothing; the rest still apply. The bodies as opened are merged too
    /// ([`Design::source_body`], the feature graph), whether partitions cut them or not.
    pub fn merged(mut self, merges: &[Merge]) -> Result<Design> {
        if merges.len() > MAX_MERGES {
            return Err(Error::Invalid(format!("a design takes at most {MAX_MERGES} face merges")));
        }
        let mut any = false;
        for m in merges {
            let feature = format!("merge {}", m.name);
            let fail = |message: String| FeatureError { feature: feature.clone(), message, bodies: vec![m.body.clone()] };
            if let Some(i) = m.into.as_ref().filter(|i| !m.faces.contains(i)) {
                self.errors.push(fail(format!("the face it keeps, `{i}`, is not one it merges")));
                continue;
            }
            if let Some(body) = self.bodies.iter_mut().find(|b| b.name == m.body) {
                let name = m.into.clone().unwrap_or_else(|| merged_face_name(&m.name));
                match body.merge_faces(&name, &m.faces) {
                    Ok(merged) => {
                        *body = merged;
                        any = true;
                    }
                    Err(e) => self.errors.push(fail(e.to_string())),
                }
                continue;
            }
            let prefix = format!("{}/", m.body);
            let cells: Vec<usize> = (0..self.bodies.len()).filter(|&i| self.bodies.get(i).is_some_and(|b| b.name.starts_with(&prefix))).collect();
            if cells.is_empty() {
                self.errors.push(fail(format!("no body `{}`", m.body)));
                continue;
            }
            // Each cell merged on its own; the merge fails as a whole when no cell has two of its
            // faces touching, or one of them refuses it.
            let mut next = Vec::with_capacity(cells.len());
            let mut error = None;
            for &i in &cells {
                let Some(cell) = self.bodies.get(i) else { continue };
                match cell_merges(cell, m) {
                    Ok(Some(merged)) => next.push((i, merged)),
                    Ok(None) => {}
                    Err(e) => {
                        error = Some(format!("cell `{}`: {e}", cell.name));
                        break;
                    }
                }
            }
            match (error, next.is_empty()) {
                (Some(e), _) => self.errors.push(fail(e)),
                (None, true) => self.errors.push(fail(format!("no cell of `{}` has two of its faces touching", m.body))),
                (None, false) => {
                    for (i, merged) in next {
                        if let Some(slot) = self.bodies.get_mut(i) {
                            *slot = merged;
                        }
                    }
                    any = true;
                }
            }
        }
        if any {
            let opened = self
                .source
                .iter()
                .map(|s| match self.bodies.iter().find(|b| b.name == s.name) {
                    Some(b) => b.clone(),
                    // Cut into cells: the body as opened with its own merges (each one that
                    // applies; the cells reported the rest).
                    None => merges.iter().filter(|m| m.body == s.name).fold(s.clone(), |shape, m| {
                        let name = m.into.clone().unwrap_or_else(|| merged_face_name(&m.name));
                        shape.merge_faces(&name, &m.faces).unwrap_or(shape)
                    }),
                })
                .collect();
            self.opened = std::sync::Arc::new(opened);
        }
        Ok(self)
    }
}

/// Merge `m` on a cell partitions cut from its body: the faces it names that the cell has (or
/// pieces of them), in the groups that touch there, each group of two or more made one face —
/// named after the face it keeps (`into`; that face's piece in the cell), else `merge:<name>`
/// (`merge:<name>#0`, `#1`, … when the cell has several groups, as a partition names pieces).
/// `None` when the cell has no two of them touching.
fn cell_merges(cell: &Shape, m: &Merge) -> Result<Option<Shape>> {
    let mut members = BTreeSet::new();
    for f in &m.faces {
        members.extend(cell.faces_named(f));
    }
    let groups: Vec<BTreeSet<u32>> = cell.components(&members)?.into_iter().filter(|g| g.len() > 1).collect();
    if groups.is_empty() {
        return Ok(None);
    }
    let kept: BTreeSet<u32> = m.into.as_deref().map(|i| cell.faces_named(i).into_iter().collect()).unwrap_or_default();
    // Names first (indices change as faces merge), then the merges by those names.
    let unnamed = groups.iter().filter(|g| g.is_disjoint(&kept)).count();
    let mut k = 0;
    let mut plan = Vec::with_capacity(groups.len());
    for g in &groups {
        let names: Vec<String> = g.iter().map(|&f| cell.face_label(f)).collect();
        let name = match g.intersection(&kept).next() {
            Some(&f) => cell.face_label(f),
            None if unnamed == 1 => merged_face_name(&m.name),
            None => {
                k += 1;
                format!("{}#{}", merged_face_name(&m.name), k - 1)
            }
        };
        plan.push((name, names));
    }
    let mut shape = cell.clone();
    for (name, faces) in plan {
        shape = shape.merge_faces(&name, &faces)?;
    }
    Ok(Some(shape))
}

impl Shape {
    /// The faces with this name: the face itself, else the pieces a partition cut it into
    /// (`<name>#k`).
    pub fn faces_named(&self, name: &str) -> Vec<u32> {
        let exact: Vec<u32> =
            (0..self.faces.len()).filter(|&i| self.faces.get(i).is_some_and(|f| f.name == name)).filter_map(|i| u32::try_from(i).ok()).collect();
        if !exact.is_empty() {
            return exact;
        }
        (0..self.faces.len())
            .filter(|&i| self.faces.get(i).is_some_and(|f| crate::partition::base_name(&f.name) == name && f.name != name))
            .filter_map(|i| u32::try_from(i).ok())
            .collect()
    }

    /// The faces made of several B-rep faces (by face index), each with how many.
    pub fn merged_faces(&self) -> Vec<(u32, usize)> {
        let mut count = vec![0usize; self.faces.len()];
        for &f in &self.face_map {
            if let Some(c) = count.get_mut(f as usize) {
                *c += 1;
            }
        }
        count.into_iter().enumerate().filter(|(_, c)| *c > 1).filter_map(|(i, c)| Some((u32::try_from(i).ok()?, c))).collect()
    }

    /// This shape with the faces named `faces` (two or more, touching along edges) made one face
    /// named `name` (a new name, or the name of one of them), in place of the first of them.
    pub fn merge_faces(&self, name: &str, faces: &[String]) -> Result<Shape> {
        if faces.len() > MAX_MERGE_FACES {
            return Err(Error::Invalid(format!("a merge takes at most {MAX_MERGE_FACES} faces")));
        }
        if !self.cuts.is_empty() {
            return Err(Error::Invalid("merges come before splits: a merge cannot take a split face".into()));
        }
        let mut members = BTreeSet::new();
        for f in faces {
            let found = self.faces_named(f);
            if found.is_empty() {
                return Err(Error::Invalid(format!("body `{}` has no face `{f}`", self.name)));
            }
            members.extend(found);
        }
        if self.faces.iter().zip(0u32..).any(|(f, i)| f.name == name && !members.contains(&i)) {
            return Err(Error::Invalid(format!("body `{}` already has a face `{name}`", self.name)));
        }
        if members.len() < 2 {
            return Err(Error::Invalid("a merge needs two or more faces".into()));
        }
        self.connected(&members)?;
        // An absorb keeps the face it is named after: its classification survives the merge
        // (#546) even where the other members lie on different surfaces. A merge that keeps no
        // face folds its members' surfaces as before.
        let kept_surface = match members.iter().copied().find(|&m| self.faces.get(m as usize).is_some_and(|f| f.name == name)) {
            Some(m) => self.face_surfaces()?.get(m as usize).copied().flatten(),
            None => None,
        };
        // The faces left, the merged one where its first face was; the old index → the new one.
        let mut kept: Vec<Face> = Vec::with_capacity(self.faces.len());
        let mut renumber = Vec::with_capacity(self.faces.len());
        let mut merged: Option<u32> = None;
        for (i, f) in self.faces.iter().enumerate() {
            let index = u32::try_from(kept.len()).map_err(|_| Error::Invalid("too many faces".into()))?;
            if !u32::try_from(i).is_ok_and(|i| members.contains(&i)) {
                renumber.push(index);
                kept.push(f.clone());
                continue;
            }
            match merged {
                Some(m) => {
                    renumber.push(m);
                    if let Some(face) = kept.get_mut(m as usize) {
                        face.exact &= f.exact;
                        face.members.extend(members_of(f));
                    }
                }
                None => {
                    merged = Some(index);
                    renumber.push(index);
                    kept.push(Face { name: name.to_string(), exact: f.exact, members: members_of(f), kept_surface });
                }
            }
        }
        let brep = if self.face_map.is_empty() { (0..u32::try_from(self.faces.len()).unwrap_or(0)).collect() } else { self.face_map.clone() };
        let face_map = brep.iter().map(|&v| renumber.get(v as usize).copied().unwrap_or(v)).collect();
        Ok(Shape { name: self.name.clone(), body: self.body.clone(), faces: kept, face_map, cuts: Vec::new() })
    }

    /// The neighbour face `face` shares the longest edges with (by name): where a tiny or sliver
    /// face is best absorbed. An error for a face the body lacks or one with no neighbour.
    pub fn absorbing_neighbour(&self, face: &str) -> Result<String> {
        let own = self.faces_named(face);
        let [f] = own.as_slice() else {
            return Err(Error::Invalid(format!("body `{}` has no face `{face}` (or several pieces of it: name one)", self.name)));
        };
        let mesh = self.tessellate((self.size()? * EDGE_TOLERANCE).max(1e-4))?;
        let mut shared: std::collections::BTreeMap<u32, f64> = std::collections::BTreeMap::new();
        for (points, faces) in mesh.edges.iter().zip(&mesh.edge_faces) {
            if !faces.contains(f) {
                continue;
            }
            let length: f64 = points.windows(2).map(|w| (0..3).map(|k| (w[1][k] - w[0][k]).powi(2)).sum::<f64>().sqrt()).sum();
            for g in faces.iter().filter(|g| *g != f) {
                *shared.entry(*g).or_default() += length;
            }
        }
        let best = shared.into_iter().max_by(|a, b| a.1.total_cmp(&b.1).then(b.0.cmp(&a.0))).map(|(g, _)| g);
        best.map(|g| self.face_label(g)).ok_or_else(|| Error::Invalid(format!("face `{face}` of body `{}` has no neighbour to absorb it", self.name)))
    }

    /// An error unless `members` touch one another along edges, all in one piece.
    fn connected(&self, members: &BTreeSet<u32>) -> Result<()> {
        let groups = self.components(members)?;
        match groups.get(1).and_then(|g| g.first()) {
            None => Ok(()),
            Some(&f) => Err(Error::Invalid(format!(
                "face `{}` of body `{}` does not touch the others along an edge: a merge makes adjacent faces one",
                self.face_label(f),
                self.name
            ))),
        }
    }

    /// `members` in groups that touch along edges (the group of the first member first).
    fn components(&self, members: &BTreeSet<u32>) -> Result<Vec<BTreeSet<u32>>> {
        if members.is_empty() {
            return Ok(Vec::new());
        }
        let mesh = self.tessellate((self.size()? * EDGE_TOLERANCE).max(1e-4))?;
        let mut groups: Vec<BTreeSet<u32>> = Vec::new();
        for &start in members {
            if groups.iter().any(|g| g.contains(&start)) {
                continue;
            }
            let mut reached = BTreeSet::new();
            let mut todo = vec![start];
            while let Some(f) = todo.pop() {
                if !reached.insert(f) {
                    continue;
                }
                for faces in mesh.edge_faces.iter().filter(|fs| fs.contains(&f)) {
                    todo.extend(faces.iter().filter(|g| members.contains(g) && !reached.contains(g)));
                }
            }
            groups.push(reached);
        }
        Ok(groups)
    }

    /// The surfaces of the B-rep faces face `face` is made of (one for a face not merged); `None`
    /// for a free-form one.
    pub fn member_surfaces(&self, face: u32) -> Result<Vec<Option<FaceSurface>>> {
        let brep = self.brep_face_surfaces()?;
        // The b side of a split is made of the face it was split from.
        let face = self.cuts.iter().rev().fold(face, |f, c| if c.faces().1 == f { c.faces().0 } else { f });
        if self.face_map.is_empty() {
            return Ok(brep.get(face as usize).copied().into_iter().collect());
        }
        Ok(self.face_map.iter().zip(brep).filter(|(v, _)| **v == face).map(|(_, s)| s).collect())
    }

    /// `brep` (surfaces by B-rep face) by face: a merged face has its B-rep faces' surface when
    /// they share one, else none.
    pub(crate) fn virtual_surfaces(&self, brep: Vec<Option<FaceSurface>>) -> Result<Vec<Option<FaceSurface>>> {
        let mut out = self.merged_surfaces(brep)?;
        for c in &self.cuts {
            let (a, b) = c.faces();
            let s = out.get(a as usize).copied().flatten();
            if out.len() <= b as usize {
                out.resize(b as usize + 1, None);
            }
            if let Some(slot) = out.get_mut(b as usize) {
                *slot = s;
            }
        }
        Ok(out)
    }

    fn merged_surfaces(&self, brep: Vec<Option<FaceSurface>>) -> Result<Vec<Option<FaceSurface>>> {
        if self.face_map.is_empty() {
            return Ok(brep);
        }
        let tol = (self.size()? * SAME_SURFACE).max(fieldcraft_model::units::min_geometry_length());
        let mut out: Vec<Option<Option<FaceSurface>>> = vec![None; self.faces.len()];
        for (b, &v) in self.face_map.iter().enumerate() {
            let surface = brep.get(b).copied().flatten();
            if let Some(slot) = out.get_mut(v as usize) {
                *slot = match *slot {
                    None => Some(surface),
                    Some(s) if crate::partition::same_surface(s, surface, tol) && s.is_some() => Some(s),
                    Some(_) => Some(None),
                };
            }
        }
        // An absorbed face keeps the classification of the face it was merged into (#546), however
        // its other members lie; a merge that kept no face keeps the fold's verdict.
        for (v, face) in self.faces.iter().enumerate() {
            if let Some(s) = face.kept_surface
                && let Some(slot) = out.get_mut(v)
            {
                *slot = Some(Some(s));
            }
        }
        Ok(out.into_iter().map(Option::flatten).collect())
    }

    /// `mesh` (a tessellation of the B-rep) with the merged faces: triangles tagged by face,
    /// the edges inside a merged face dropped and the pieces of an edge it leaves joined.
    pub(crate) fn virtual_mesh(&self, mesh: SurfaceMesh) -> SurfaceMesh {
        let mut mesh = self.merged_mesh(mesh);
        crate::face_split::apply(&mut mesh, &self.cuts);
        mesh
    }

    fn merged_mesh(&self, mut mesh: SurfaceMesh) -> SurfaceMesh {
        if self.face_map.is_empty() {
            return mesh;
        }
        let map = |f: u32| self.face_map.get(f as usize).copied().unwrap_or(f);
        for f in &mut mesh.tri_face {
            *f = map(*f);
        }
        let mut edges = Vec::with_capacity(mesh.edges.len());
        let mut edge_faces = Vec::with_capacity(mesh.edges.len());
        for (points, faces) in mesh.edges.into_iter().zip(mesh.edge_faces) {
            let two = faces.windows(2).any(|w| w[0] != w[1]);
            let mut mapped: Vec<u32> = faces.iter().map(|&f| map(f)).collect();
            mapped.dedup();
            // Between two B-rep faces of one merged face: inside it.
            if two && mapped.windows(2).all(|w| w[0] == w[1]) {
                continue;
            }
            edges.push(points);
            edge_faces.push(mapped);
        }
        join_pieces(&mut edges, &mut edge_faces, (2.0 * mesh.weld_distance).max(fieldcraft_model::units::min_geometry_length()));
        mesh.edges = edges;
        mesh.edge_faces = edge_faces;
        mesh
    }
}

/// The B-rep faces `f` is made of, by name: itself, or a merged face's.
fn members_of(f: &Face) -> Vec<String> {
    if f.members.is_empty() { vec![f.name.clone()] } else { f.members.clone() }
}

fn close(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
    (0..3).map(|k| (a[k] - b[k]).powi(2)).sum::<f64>().sqrt() <= tol
}

/// Join edges between the same faces end to end where they meet and no other edge ends: the
/// pieces of one edge a merge left (a fillet's pieces each had their own edge along the face
/// beside them).
fn join_pieces(edges: &mut Vec<Vec<[f64; 3]>>, faces: &mut Vec<Vec<u32>>, tol: f64) {
    let key = |f: &Vec<u32>| {
        let mut k = f.clone();
        k.sort_unstable();
        k.dedup();
        k
    };
    // Every edge end, to count the edges meeting at a point.
    let ends =
        |edges: &[Vec<[f64; 3]>], p: [f64; 3]| edges.iter().flat_map(|e| [e.first(), e.last()]).flatten().filter(|q| close(**q, p, tol)).count();
    let mut joined = 0;
    'again: while joined < MAX_JOINED {
        for i in 0..edges.len() {
            for j in (i + 1)..edges.len() {
                let (Some(a), Some(b)) = (edges.get(i), edges.get(j)) else { continue };
                if faces.get(i).map(key) != faces.get(j).map(key) || a.len() < 2 || b.len() < 2 {
                    continue;
                }
                let (a0, a1, b0, b1) = (a[0], a[a.len() - 1], b[0], b[b.len() - 1]);
                // Closed already, or would close into a loop at both ends: leave it.
                if close(a0, a1, tol) || close(b0, b1, tol) {
                    continue;
                }
                let at = if close(a1, b0, tol) {
                    Some((a1, false, false))
                } else if close(a1, b1, tol) {
                    Some((a1, false, true))
                } else if close(a0, b1, tol) {
                    Some((a0, true, true))
                } else if close(a0, b0, tol) {
                    Some((a0, true, false))
                } else {
                    None
                };
                let Some((p, flip_a, flip_b)) = at else { continue };
                if ends(edges, p) != 2 {
                    continue;
                }
                let mut b = edges.remove(j);
                faces.remove(j);
                let Some(a) = edges.get_mut(i) else { return };
                if flip_a {
                    a.reverse();
                }
                if flip_b {
                    b.reverse();
                }
                a.extend(b.into_iter().skip(1));
                joined += 1;
                continue 'again;
            }
        }
        break;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sample::{self, Plate};
    use crate::{FaceSurface, open_design};

    fn plate() -> Design {
        open_design(&sample::plate_design(&Plate::default()).unwrap()).unwrap()
    }

    /// The plate's hole is two half-cylinders.
    fn halves(d: &Design) -> Vec<String> {
        let body = &d.bodies[0];
        let s = body.face_surfaces().unwrap();
        body.faces().iter().zip(s).filter(|(_, s)| matches!(s, Some(FaceSurface::Cylinder { .. }))).map(|(f, _)| f.name.clone()).collect()
    }

    /// Merging the hole's halves makes one cylindrical face: one tag on its triangles, its seams
    /// gone from the edges, the rim arcs joined into one circle each, and the edge and face names
    /// follow; the B-rep keeps its faces.
    #[test]
    fn merged_halves_are_one_face_without_seams() {
        let d = plate();
        let hole = halves(&d);
        assert_eq!(hole.len(), 2, "{hole:?}");
        let before = d.bodies[0].tessellate(0.01).unwrap();
        let body = d.bodies[0].name.clone();
        let d = d.merged(&[Merge { name: "Hole".into(), body, faces: hole.clone(), into: None }]).unwrap();
        assert!(d.errors.is_empty(), "{:?}", d.errors);
        let body = &d.bodies[0];
        assert_eq!(body.faces().len(), 7);
        let f = body.face_index("merge:Hole").unwrap();
        assert!(body.face_index(&hole[0]).is_none());
        assert!(matches!(body.face_surfaces().unwrap()[f as usize], Some(FaceSurface::Cylinder { .. })));
        assert_eq!(body.merged_faces(), vec![(f, 2)]);
        let mesh = body.tessellate(0.01).unwrap();
        assert_eq!(mesh.triangles, before.triangles);
        let tags: BTreeSet<u32> = mesh.tri_face.iter().copied().collect();
        assert_eq!(tags.len(), 7);
        // 12 box edges and the two rims, each now one circle.
        let edges = body.edges().unwrap();
        assert_eq!(edges.len(), 14, "{:?}", edges.iter().map(|e| &e.name).collect::<Vec<_>>());
        assert!(edges.iter().all(|e| e.faces.iter().all(|&g| (g as usize) < 7)));
        let rims: Vec<_> = edges.iter().filter(|e| e.faces.contains(&f)).collect();
        assert_eq!(rims.len(), 2);
        for r in rims {
            assert!(r.name.contains("merge:Hole"), "{}", r.name);
            assert!(close(r.points[0], *r.points.last().unwrap(), 1e-6), "a rim is a closed circle");
        }
        assert_eq!(body.measure().unwrap().faces, 7);
    }

    /// A merge of faces that do not touch, of one face, of a missing face or body is an error
    /// listed with the design, which keeps its faces.
    #[test]
    fn bad_merges_are_listed_and_change_nothing() {
        let d = plate();
        let body = d.bodies[0].name.clone();
        let names: Vec<String> = d.bodies[0].faces().iter().map(|f| f.name.clone()).collect();
        let s = d.bodies[0].face_surfaces().unwrap();
        let planes: Vec<&String> = names
            .iter()
            .zip(&s)
            .filter(|(_, s)| matches!(s, Some(FaceSurface::Plane { normal, .. }) if normal[2].abs() > 0.5))
            .map(|(n, _)| n)
            .collect();
        assert_eq!(planes.len(), 2, "top and bottom");
        let m = |name: &str, body: &str, faces: Vec<String>| Merge { name: name.into(), body: body.into(), faces, into: None };
        let merged = d
            .clone()
            .merged(&[
                m("Apart", &body, planes.iter().map(|s| s.to_string()).collect()),
                m("One", &body, vec![names[0].clone()]),
                m("Missing", &body, vec![names[0].clone(), "nope".into()]),
                m("Nobody", "Nobody", vec![names[0].clone(), names[1].clone()]),
            ])
            .unwrap();
        let failed: Vec<&str> = merged.errors.iter().map(|e| e.feature.as_str()).collect();
        assert_eq!(failed, ["merge Apart", "merge One", "merge Missing", "merge Nobody"]);
        assert!(merged.errors[0].message.contains("does not touch"), "{}", merged.errors[0].message);
        assert_eq!(merged.bodies[0].faces().len(), 8);
        // A second merge into a name already there is refused too.
        let twice = d.merged(&[m("H", &body, halves(&plate())), m("H", &body, vec![names[0].clone(), names[1].clone()])]).unwrap();
        assert_eq!(twice.errors.len(), 1);
    }

    /// Merges stack: a merged face merges again with its neighbour, and the surfaces differ then.
    #[test]
    fn a_merged_face_merges_again() {
        let d = plate();
        let body = d.bodies[0].name.clone();
        let hole = halves(&d);
        let d = d.merged(&[Merge { name: "Hole".into(), body: body.clone(), faces: hole, into: None }]).unwrap();
        let b = &d.bodies[0];
        let f = b.face_index("merge:Hole").unwrap();
        let mesh = b.tessellate(0.01).unwrap();
        let next = mesh.edge_faces.iter().find(|fs| fs.contains(&f)).and_then(|fs| fs.iter().find(|&&g| g != f)).copied().unwrap();
        let next = b.faces()[next as usize].name.clone();
        let d = d.merged(&[Merge { name: "More".into(), body, faces: vec!["merge:Hole".into(), next], into: None }]).unwrap();
        assert!(d.errors.is_empty(), "{:?}", d.errors);
        let b = &d.bodies[0];
        let f = b.face_index("merge:More").unwrap();
        assert_eq!(b.faces().len(), 6);
        assert_eq!(b.merged_faces(), vec![(f, 3)]);
        assert_eq!(b.face_surfaces().unwrap()[f as usize], None);
    }

    /// A merge on a body the partitions cut applies to each cell: the plate cut across its hole
    /// leaves a piece of each half in each cell, and each cell's two pieces are one face there,
    /// named after the face kept; the body as opened keeps the merge too.
    #[test]
    fn a_merge_applies_to_the_cells_of_a_partitioned_body() {
        let d = plate();
        let body = d.bodies[0].name.clone();
        let hole = halves(&d);
        let cut = crate::partition::Cut { name: "Mid".into(), point: [0.0, 0.0, 5.0], normal: [0.0, 0.0, 1.0], bodies: Vec::new() };
        let absorb = Merge { name: "A".into(), body: body.clone(), faces: hole.clone(), into: Some(hole[1].clone()) };
        let d = d.partitioned(&[cut]).unwrap().merged(&[absorb]).unwrap();
        assert!(d.errors.is_empty(), "{:?}", d.errors);
        assert_eq!(d.bodies.len(), 2);
        for cell in &d.bodies {
            let [(f, 2)] = cell.merged_faces()[..] else { panic!("{}: {:?}", cell.name, cell.merged_faces()) };
            assert_eq!(cell.faces()[f as usize].name, hole[1], "{}", cell.name);
            assert!(cell.faces_named(&hole[0]).is_empty(), "{}", cell.name);
        }
        let opened = d.source_body(&body).unwrap();
        assert_eq!(opened.merged_faces().len(), 1);
        // A merge whose faces touch in no cell is listed.
        let names: Vec<String> = plate().bodies[0].faces().iter().map(|f| f.name.clone()).collect();
        let cut = crate::partition::Cut { name: "Mid".into(), point: [0.0, 0.0, 5.0], normal: [0.0, 0.0, 1.0], bodies: Vec::new() };
        let s = plate().bodies[0].face_surfaces().unwrap();
        let caps: Vec<String> = names
            .iter()
            .zip(&s)
            .filter(|(_, s)| matches!(s, Some(FaceSurface::Plane { normal, .. }) if normal[2].abs() > 0.5))
            .map(|(n, _)| n.clone())
            .collect();
        let apart = Merge { name: "Apart".into(), body, faces: caps, into: None };
        let d = plate().partitioned(&[cut]).unwrap().merged(&[apart]).unwrap();
        assert_eq!(d.errors.iter().map(|e| e.feature.as_str()).collect::<Vec<_>>(), ["merge Apart"]);
    }

    /// The neighbour a half of the hole is best absorbed into is the other half (two 10 mm seams,
    /// against a 15.7 mm rim arc with each cap); absorbing it keeps the other half's name.
    #[test]
    fn an_absorbed_face_leaves_its_neighbour_by_name() {
        let d = plate();
        let body = d.bodies[0].name.clone();
        let hole = halves(&d);
        let best = d.bodies[0].absorbing_neighbour(&hole[0]).unwrap();
        assert_eq!(best, hole[1]);
        let d = d.merged(&[Merge { name: "A".into(), body, faces: vec![hole[0].clone(), hole[1].clone()], into: Some(hole[1].clone()) }]).unwrap();
        assert!(d.errors.is_empty(), "{:?}", d.errors);
        let b = &d.bodies[0];
        assert_eq!(b.faces().len(), 7);
        assert!(b.face_index(&hole[1]).is_some() && b.face_index(&hole[0]).is_none());
        assert!(b.absorbing_neighbour("nope").is_err());
    }

    /// An absorbed face keeps the classification of the face it was absorbed into (#546): the
    /// hole's halves merged into the top plane leave one planar face, not a free-form one, so the
    /// meshers keep the route they had before the absorb. A plain merge of the same faces stays
    /// free-form: its faces lie on different surfaces and no face was kept.
    #[test]
    fn an_absorbed_face_keeps_the_surface_of_the_face_absorbed_into() {
        let d = plate();
        let body = d.bodies[0].name.clone();
        let hole = halves(&d);
        let top = d.bodies[0]
            .face_surfaces()
            .unwrap()
            .iter()
            .zip(&d.bodies[0].faces().iter().map(|f| f.name.clone()).collect::<Vec<_>>())
            .find(|(s, _)| matches!(s, Some(FaceSurface::Plane { normal, .. }) if normal[2] > 0.5))
            .map(|(_, n)| n.clone())
            .unwrap();
        // First absorb the hole's halves into one another (one cylindrical face), that face into
        // the top plane: the chain ends planar.
        let d = d
            .merged(&[Merge { name: "Hole".into(), body: body.clone(), faces: hole, into: None }])
            .unwrap()
            .merged(&[Merge { name: "A".into(), body: body.clone(), faces: vec!["merge:Hole".into(), top.clone()], into: Some(top.clone()) }])
            .unwrap();
        assert!(d.errors.is_empty(), "{:?}", d.errors);
        let b = &d.bodies[0];
        let f = b.face_index(&top).unwrap();
        assert!(matches!(b.face_surfaces().unwrap()[f as usize], Some(FaceSurface::Plane { .. })), "{:?}", b.face_surfaces().unwrap()[f as usize]);
        // The same two merges without a face kept stay free-form: a plain merge of a cylinder and
        // a plane lies on no one surface.
        let plain = plate()
            .merged(&[Merge { name: "Hole".into(), body: body.clone(), faces: halves(&plate()), into: None }])
            .unwrap()
            .merged(&[Merge { name: "A".into(), body, faces: vec!["merge:Hole".into(), top], into: None }])
            .unwrap();
        assert!(plain.errors.is_empty(), "{:?}", plain.errors);
        let f = plain.bodies[0].face_index("merge:A").unwrap();
        assert_eq!(plain.bodies[0].face_surfaces().unwrap()[f as usize], None);
    }
}
