//! Global validity of a finished tetrahedral mesh, checked on the test side, independently of the
//! mesher's own gate (`fieldcraft_mesh::tet::linear_gate`): the corpus checks every tet mesh it
//! makes against these invariants (docs/meshing/corpus.md).
//!
//! - every tetrahedron has a positive volume (the exact sign of `det[b − a, c − a, d − a]`);
//! - every triangle is a face of one tetrahedron (on the boundary) or two (inside), and the two
//!   lie on opposite sides of it;
//! - the boundary faces form a closed, consistently oriented 2-manifold: every boundary edge is
//!   in exactly two boundary faces that run it in opposite directions, and the boundary faces
//!   round each boundary vertex form a single fan;
//! - no tetrahedron is joined to the rest only by vertices or edges: the tetrahedra round each
//!   vertex are connected through faces that hold the vertex, and the mesh has one
//!   face-connected piece per solid;
//! - the boundary faces the mesher lists are the faces held by one tetrahedron.
//!
//! These are topological and sign checks, so a fold (a face held by three tetrahedra, two layers
//! on one side of a face, a pinched skin) shows even when the volumes add up.

use std::collections::{BTreeMap, HashMap};

use fieldcraft_mesh::predicates::{Point3, Sign, tet_volume_sign};
use fieldcraft_mesh::tet::TaggedFace;

/// The faces of a tet4 by its local vertices, each ordered so its right-hand normal points out of
/// a positive tetrahedron (the model's tet4 face order).
const FACES: [[usize; 3]; 4] = [[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]];
/// Violations listed with their place; the rest are counted.
pub const MAX_LISTED: usize = 20;

/// One broken invariant and where.
#[derive(Debug, Clone, PartialEq)]
pub struct Violation {
    /// The invariant ([`Report::counts`] key).
    pub kind: &'static str,
    pub what: String,
    /// A point of the place: the element's, face's or edge's centroid, or the vertex.
    pub at: [f64; 3],
}

impl std::fmt::Display for Violation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let [x, y, z] = self.at;
        write!(f, "{} at ({x:.6}, {y:.6}, {z:.6})", self.what)
    }
}

/// The result of [`check`].
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Report {
    /// Violations by invariant.
    pub counts: BTreeMap<&'static str, usize>,
    /// The first [`MAX_LISTED`] violations.
    pub listed: Vec<Violation>,
}

impl Report {
    pub fn is_valid(&self) -> bool {
        self.counts.is_empty()
    }

    pub fn total(&self) -> usize {
        self.counts.values().sum()
    }

    /// Add `other`'s violations, naming their place by `scope` (a body, a stage).
    pub fn absorb(&mut self, scope: &str, other: Report) {
        for (k, n) in other.counts {
            *self.counts.entry(k).or_default() += n;
        }
        for v in other.listed {
            if self.listed.len() < MAX_LISTED {
                self.listed.push(Violation { what: format!("{scope}: {}", v.what), ..v });
            }
        }
    }

    fn add(&mut self, kind: &'static str, what: String, at: [f64; 3]) {
        *self.counts.entry(kind).or_default() += 1;
        if self.listed.len() < MAX_LISTED {
            self.listed.push(Violation { kind, what, at });
        }
    }

    /// Counts by invariant, then the listed places.
    pub fn summary(&self) -> String {
        let counts: Vec<String> = self.counts.iter().map(|(k, n)| format!("{n} {k}")).collect();
        let places: Vec<String> = self.listed.iter().map(Violation::to_string).collect();
        format!("{}: {}", counts.join(", "), places.join("; "))
    }
}

/// Invariant names ([`Violation::kind`]).
pub mod kind {
    pub const ELEMENT: &str = "malformed element";
    pub const VOLUME: &str = "tetrahedron without positive volume";
    pub const OVERFULL_FACE: &str = "face held by more than two tetrahedra";
    pub const SAME_SIDE: &str = "two tetrahedra on one side of a face";
    pub const BOUNDARY_EDGE: &str = "boundary edge not in exactly two boundary faces";
    pub const ORIENTATION: &str = "boundary oriented inconsistently";
    pub const PINCHED_VERTEX: &str = "boundary pinched at a vertex";
    pub const VERTEX_STAR: &str = "tetrahedra joined by a vertex or an edge only";
    pub const PIECES: &str = "face-connected pieces other than the solids";
    pub const BOUNDARY_LIST: &str = "listed boundary differs from the faces held once";
}

fn centroid(ps: &[[f64; 3]]) -> [f64; 3] {
    let n = ps.len().max(1) as f64;
    let mut c = [0.0; 3];
    for p in ps {
        for (c, x) in c.iter_mut().zip(p) {
            *c += x / n;
        }
    }
    c
}

/// A union-find over `n` items.
struct Sets(Vec<usize>);

impl Sets {
    fn new(n: usize) -> Sets {
        Sets((0..n).collect())
    }

    fn root(&mut self, mut i: usize) -> usize {
        while let Some(&p) = self.0.get(i) {
            if p == i {
                break;
            }
            let grand = self.0.get(p).copied().unwrap_or(p);
            if let Some(slot) = self.0.get_mut(i) {
                *slot = grand;
            }
            i = grand;
        }
        i
    }

    fn join(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.root(a), self.root(b));
        if ra != rb
            && let Some(slot) = self.0.get_mut(ra.max(rb))
        {
            *slot = ra.min(rb);
        }
    }
}

/// Where a triangle (by its sorted vertices) is a face: (tetrahedron, its local face).
type Holders = HashMap<[u32; 3], Vec<(usize, usize)>>;

/// Check a tet mesh: `elements` give their four corners first (tet4 or tet10 in the model's node
/// order), `boundary` is the mesher's list of boundary faces, and `solids` is the number of
/// separate solids the mesh fills (the boundary surface's outer shells).
pub fn check(points: &[[f64; 3]], elements: &[Vec<u32>], boundary: &[TaggedFace], solids: usize) -> Report {
    let mut r = Report::default();
    let at = |v: u32| points.get(v as usize).copied().unwrap_or([f64::NAN; 3]);
    // The corners of the well-formed elements, by element index.
    let mut tets: Vec<Option<[u32; 4]>> = Vec::with_capacity(elements.len());
    for (i, e) in elements.iter().enumerate() {
        let corners = match e.get(..4) {
            Some(&[a, b, c, d]) => [a, b, c, d],
            _ => {
                r.add(kind::ELEMENT, format!("element {i} has {} nodes", e.len()), [f64::NAN; 3]);
                tets.push(None);
                continue;
            }
        };
        let mut sorted = corners;
        sorted.sort_unstable();
        if corners.iter().any(|&v| v as usize >= points.len()) || sorted.windows(2).any(|w| w[0] == w[1]) {
            r.add(kind::ELEMENT, format!("element {i} has corners {corners:?}"), centroid(&corners.map(at)));
            tets.push(None);
            continue;
        }
        let p = corners.map(|v| Point3::new(at(v)[0], at(v)[1], at(v)[2]));
        let sign = match p {
            [Ok(a), Ok(b), Ok(c), Ok(d)] => Some(tet_volume_sign(&a, &b, &c, &d)),
            _ => None,
        };
        if sign != Some(Sign::Positive) {
            let what = match sign {
                Some(Sign::Negative) => "inverted",
                Some(_) => "flat",
                None => "outside the exact-predicate domain",
            };
            r.add(kind::VOLUME, format!("tetrahedron {i} is {what}"), centroid(&corners.map(at)));
        }
        tets.push(Some(corners));
    }

    // Every face and the tetrahedra that hold it.
    let mut holders: Holders = HashMap::with_capacity(tets.len() * 2);
    for (t, corners) in tets.iter().enumerate() {
        let Some(c) = corners else { continue };
        for (k, f) in FACES.iter().enumerate() {
            let mut key = f.map(|j| c[j]);
            key.sort_unstable();
            holders.entry(key).or_default().push((t, k));
        }
    }
    let oriented = |t: usize, k: usize| -> Option<[u32; 3]> {
        let c = tets.get(t).copied().flatten()?;
        FACES.get(k).map(|f| f.map(|j| c[j]))
    };
    // Even permutations of a face keep its orientation.
    let same_turn = |a: [u32; 3], b: [u32; 3]| a == b || a == [b[1], b[2], b[0]] || a == [b[2], b[0], b[1]];

    let mut faces: Vec<(&[u32; 3], &Vec<(usize, usize)>)> = holders.iter().collect();
    faces.sort_unstable_by_key(|(key, _)| **key);
    let mut skin: Vec<(usize, usize, [u32; 3])> = Vec::new();
    // (vertex, tet) pairs joined through a face holding the vertex; index t * 4 + local corner.
    let mut star = Sets::new(tets.len() * 4);
    let mut pieces = Sets::new(tets.len());
    let local = |t: usize, v: u32| tets.get(t).copied().flatten().and_then(|c| c.iter().position(|&x| x == v)).map(|j| t * 4 + j);
    for (key, held) in &faces {
        let place = centroid(&key.map(at));
        match held.as_slice() {
            [(t, k)] => {
                if let Some(f) = oriented(*t, *k) {
                    skin.push((*t, *k, f));
                }
            }
            [(t1, k1), (t2, k2)] => {
                if let (Some(a), Some(b)) = (oriented(*t1, *k1), oriented(*t2, *k2))
                    && same_turn(a, b)
                {
                    r.add(kind::SAME_SIDE, format!("tetrahedra {t1} and {t2} lie on the same side of face {key:?}"), place);
                }
                pieces.join(*t1, *t2);
                for &v in key.iter() {
                    if let (Some(a), Some(b)) = (local(*t1, v), local(*t2, v)) {
                        star.join(a, b);
                    }
                }
            }
            more => {
                let list: Vec<usize> = more.iter().map(|(t, _)| *t).collect();
                r.add(kind::OVERFULL_FACE, format!("face {key:?} is held by {} tetrahedra {list:?}", more.len()), place);
            }
        }
    }

    // The boundary as a surface: directed edges and the fan round each vertex.
    let mut directed: HashMap<(u32, u32), usize> = HashMap::with_capacity(skin.len() * 3);
    for (_, _, f) in &skin {
        for j in 0..3 {
            *directed.entry((f[j], f[(j + 1) % 3])).or_default() += 1;
        }
    }
    let mut edges: Vec<(u32, u32)> = directed.keys().map(|&(a, b)| (a.min(b), a.max(b))).collect();
    edges.sort_unstable();
    edges.dedup();
    for (a, b) in edges {
        let (ab, ba) = (directed.get(&(a, b)).copied().unwrap_or(0), directed.get(&(b, a)).copied().unwrap_or(0));
        let place = centroid(&[at(a), at(b)]);
        if ab + ba != 2 {
            r.add(kind::BOUNDARY_EDGE, format!("boundary edge {a}–{b} is in {} boundary faces", ab + ba), place);
        } else if ab != 1 {
            r.add(kind::ORIENTATION, format!("both boundary faces on edge {a}–{b} run it the same way"), place);
        }
    }
    // Round each boundary vertex, its boundary faces chain into fans through their shared edges;
    // more than one fan is a pinch. A vertex's fan link: face (v, a, b) joins a and b.
    let mut link: BTreeMap<u32, Vec<(u32, u32)>> = BTreeMap::new();
    for (_, _, f) in &skin {
        for j in 0..3 {
            link.entry(f[j]).or_default().push((f[(j + 1) % 3], f[(j + 2) % 3]));
        }
    }
    for (v, around) in &link {
        let mut ids: Vec<u32> = around.iter().flat_map(|&(a, b)| [a, b]).collect();
        ids.sort_unstable();
        ids.dedup();
        let index = |x: u32| ids.binary_search(&x).unwrap_or(0);
        let mut fans = Sets::new(ids.len());
        for &(a, b) in around {
            fans.join(index(a), index(b));
        }
        let count = (0..ids.len()).filter(|&i| fans.root(i) == i).count();
        if count > 1 {
            r.add(kind::PINCHED_VERTEX, format!("the boundary round vertex {v} makes {count} separate fans"), at(*v));
        }
    }

    // The tetrahedra round each vertex, connected through faces holding it.
    let mut round: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for (t, corners) in tets.iter().enumerate() {
        if let Some(c) = corners {
            for (j, &v) in c.iter().enumerate() {
                round.entry(v).or_default().push(t * 4 + j);
            }
        }
    }
    for (v, slots) in &round {
        let mut roots: Vec<usize> = slots.iter().map(|&s| star.root(s)).collect();
        roots.sort_unstable();
        roots.dedup();
        if roots.len() > 1 {
            r.add(kind::VERTEX_STAR, format!("the tetrahedra round vertex {v} make {} groups joined by no face", roots.len()), at(*v));
        }
    }
    let live: Vec<usize> = (0..tets.len()).filter(|&t| tets.get(t).is_some_and(Option::is_some)).collect();
    let mut piece_roots: BTreeMap<usize, usize> = BTreeMap::new();
    for &t in &live {
        *piece_roots.entry(pieces.root(t)).or_default() += 1;
    }
    if !live.is_empty() && piece_roots.len() != solids {
        // The smallest piece is the likely stray.
        let smallest = piece_roots.iter().min_by_key(|&(_, n)| *n).map(|(&root, &n)| (root, n));
        let place = smallest.and_then(|(root, _)| tets.get(root).copied().flatten()).map_or([f64::NAN; 3], |c| centroid(&c.map(at)));
        let n = smallest.map_or(0, |(_, n)| n);
        r.add(
            kind::PIECES,
            format!("the tetrahedra make {} face-connected pieces for {solids} solid(s); the smallest has {n}", piece_roots.len()),
            place,
        );
    }

    // The mesher's boundary list against the faces held once.
    let mut listed: Vec<(usize, usize)> = boundary.iter().map(|b| (b.tet as usize, usize::from(b.opposite))).collect();
    listed.sort_unstable();
    let mut held_once: Vec<(usize, usize)> = skin.iter().map(|&(t, k, _)| (t, opposite_corner(k))).collect();
    held_once.sort_unstable();
    if listed != held_once {
        let extra = listed.iter().find(|x| held_once.binary_search(x).is_err());
        let missing = held_once.iter().find(|x| listed.binary_search(x).is_err());
        let first = extra.or(missing).copied();
        let place = first.and_then(|(t, _)| tets.get(t).copied().flatten()).map_or([f64::NAN; 3], |c| centroid(&c.map(at)));
        r.add(
            kind::BOUNDARY_LIST,
            format!(
                "the mesher lists {} boundary faces, {} faces are held once; first difference (tetrahedron, opposite corner) {first:?}",
                listed.len(),
                held_once.len()
            ),
            place,
        );
    }
    r
}

/// The local corner a tet4 face (by [`FACES`] index) leaves out.
fn opposite_corner(k: usize) -> usize {
    FACES.get(k).and_then(|f| (0..4).find(|j| !f.contains(j))).unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The unit cube's corners, then its centre.
    fn cube() -> Vec<[f64; 3]> {
        let mut p: Vec<[f64; 3]> = (0..8).map(|i| [f64::from(i & 1), f64::from((i >> 1) & 1), f64::from((i >> 2) & 1)]).collect();
        p.push([0.5, 0.5, 0.5]);
        p
    }

    /// The cube as 12 tetrahedra from its centre to its faces' halves, all positive.
    fn cube_tets() -> Vec<Vec<u32>> {
        // Each face as two triangles, outward.
        let tris: [[u32; 3]; 12] =
            [[0, 2, 1], [1, 2, 3], [4, 5, 6], [5, 7, 6], [0, 1, 4], [1, 5, 4], [2, 6, 3], [3, 6, 7], [0, 4, 2], [2, 4, 6], [1, 3, 5], [3, 7, 5]];
        tris.iter().map(|t| vec![t[0], t[2], t[1], 8]).collect()
    }

    fn boundary_of(elements: &[Vec<u32>]) -> Vec<TaggedFace> {
        // Faces held once, listed as the mesher does.
        let mut count: HashMap<[u32; 3], usize> = HashMap::new();
        for e in elements {
            for f in FACES {
                let mut k = f.map(|j| e[j]);
                k.sort_unstable();
                *count.entry(k).or_default() += 1;
            }
        }
        let mut out = Vec::new();
        for (t, e) in elements.iter().enumerate() {
            for (k, f) in FACES.iter().enumerate() {
                let mut key = f.map(|j| e[j]);
                key.sort_unstable();
                if count[&key] == 1 {
                    out.push(TaggedFace { tet: t as u32, opposite: opposite_corner(k) as u8, face: 0 });
                }
            }
        }
        out
    }

    #[test]
    fn a_valid_mesh_passes() {
        let (p, e) = (cube(), cube_tets());
        let r = check(&p, &e, &boundary_of(&e), 1);
        assert!(r.is_valid(), "{}", r.summary());
        // The faces match the model's tet4 face order and opposite corners.
        assert_eq!((0..4).map(opposite_corner).collect::<Vec<_>>(), vec![3, 2, 0, 1]);
    }

    #[test]
    fn an_inverted_tetrahedron_is_found() {
        let (p, mut e) = (cube(), cube_tets());
        e[0].swap(0, 1);
        let r = check(&p, &e, &boundary_of(&e), 1);
        assert_eq!(r.counts.get(kind::VOLUME), Some(&1), "{}", r.summary());
    }

    /// A tetrahedron held by its vertices alone (#659): a separate piece and a split vertex star.
    #[test]
    fn a_tetrahedron_held_by_a_vertex_is_found() {
        let mut p = cube();
        let mut e = cube_tets();
        let base = p.len() as u32;
        p.extend([[2.0, 0.0, 0.0], [1.0, 1.0, -1.0], [2.0, 1.0, 0.0]]);
        e.push(vec![1, base + 2, base, base + 1]);
        let r = check(&p, &e, &boundary_of(&e), 1);
        assert!(r.counts.contains_key(kind::PIECES) && r.counts.contains_key(kind::VERTEX_STAR), "{}", r.summary());
        assert!(r.counts.contains_key(kind::PINCHED_VERTEX), "{}", r.summary());
    }

    /// A fold (#666): a third tetrahedron laid on an inner face two already hold.
    #[test]
    fn a_fold_is_found() {
        let mut p = cube();
        let mut e = cube_tets();
        p.push([0.5, 0.2, 0.2]);
        e.push(vec![0, 1, 8, 9]);
        let r = check(&p, &e, &boundary_of(&e), 1);
        assert!(r.counts.contains_key(kind::OVERFULL_FACE), "{}", r.summary());
        assert!(!r.listed[0].to_string().is_empty());
    }

    /// Two tetrahedra stacked on the same side of a face (a fold whose volumes add up).
    #[test]
    fn two_layers_on_one_side_are_found() {
        let p = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 2.0]];
        let e = vec![vec![0, 1, 2, 3], vec![0, 1, 2, 4]];
        let r = check(&p, &e, &boundary_of(&e), 1);
        assert!(r.counts.contains_key(kind::SAME_SIDE), "{}", r.summary());
    }

    #[test]
    fn a_wrong_boundary_list_is_found() {
        let (p, e) = (cube(), cube_tets());
        let mut b = boundary_of(&e);
        b.pop();
        let r = check(&p, &e, &b, 1);
        assert_eq!(r.counts.keys().copied().collect::<Vec<_>>(), vec![kind::BOUNDARY_LIST], "{}", r.summary());
    }

    #[test]
    fn malformed_elements_are_reported_not_panicked_on() {
        let p = cube();
        let e = vec![vec![0, 1, 2], vec![0, 1, 2, 99], vec![0, 0, 1, 2]];
        let r = check(&p, &e, &[], 1);
        assert_eq!(r.counts.get(kind::ELEMENT), Some(&3), "{}", r.summary());
    }
}
