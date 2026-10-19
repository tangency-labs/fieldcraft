//! Constrained Delaunay tetrahedralization of a closed boundary surface: every input triangle is
//! a union of tetrahedron faces, with no tolerance anywhere in the decisions.
//!
//! # Method
//! L. Diazzi, D. Panozzo, A. Vaxman, M. Attene, *Constrained Delaunay Tetrahedrization: A Robust
//! and Practical Approach*, ACM TOG 42(6):181, 2023 (arXiv 2309.09805; paper only, the authors'
//! code is GPL/LGPL and is not read), on the theory of H. Si, K. Gärtner, *Meshing Piecewise
//! Linear Complexes by Constrained Delaunay Tetrahedralizations*, IMR 2005, and J. R. Shewchuk,
//! *Constrained Delaunay Tetrahedralizations and Provably Good Boundary Recovery*, IMR 2002 (also
//! papers only; TetGen is AGPL).
//!
//! 1. The Delaunay tetrahedralization ([`crate::delaunay`]) of the surface's vertices and the
//!    eight corners of a box around them, so every Steiner point lies strictly inside the hull.
//! 2. **Segment recovery** ([`segments`], §3.3): every edge of the surface is a segment; a missing
//!    one is split at an implicit Steiner point (an [`Lnc`](crate::predicates::Lnc) point on the
//!    original segment, §4.3) placed by Si's rule and inserted by Bowyer–Watson with the same
//!    perturbed in-sphere test as the Delaunay step, until all segments are edges.
//! 3. **Facet recovery** ([`facets`], §3.4, §4.4): every input triangle is a facet; the
//!    tetrahedra crossing a missing facet are deleted and the two half-cavities refilled by the
//!    modified gift-wrapping, which builds the constrained Delaunay tetrahedralization of each
//!    half with exact predicates on implicit points.
//! 4. **Inside/outside** (§4.5): a flood fill from the box corners that flips at every facet
//!    triangle; the orientation of every facet triangle against its inside tetrahedron is then
//!    checked exactly.
//!
//! Every decision is an exact predicate ([`crate::predicates`]): `orient3d`, `orient2d` and
//! coordinate comparison on explicit and LNC points, and `insphere` with the Devillers–Teillaud
//! symbolic perturbation. The result is checked before it is returned: every segment is an edge,
//! no edge crosses any facet, every tetrahedron is positive (exactly), and the labelling is
//! consistent with the surface's orientation.
//!
//! # Output
//! [`TetMesh`]: the inside tetrahedra in tet4 order, the surface's vertices first (same ids) and
//! then the Steiner points, each with its segment and parameter, rounded to doubles; and for each
//! boundary tetrahedron face the input triangle it lies on (so the B-rep face tag carries over).
//! Rounding a Steiner point to the nearest doubles moves it off its segment by under an ulp and
//! can, on a nearly flat tetrahedron, flip it (§5.1); nearly cospherical input also gives
//! tetrahedra that are positive yet flat to rounding. [`repair`] removes them by flips that keep
//! the boundary (§5.1); [`TetMesh::gate`] runs the validity gate on the rounded mesh and any
//! failure left is reported, never hidden.
//!
//! # Determinism and bounds
//! Same surface → same mesh: the Delaunay step is deterministic, segments and facets are
//! processed in a fixed order, and hash maps are only looked up. Steiner points, tetrahedra,
//! passes and cavities are capped ([`CdtOptions`]); hitting a cap is an [`Error::Limit`].

mod complex;
mod facets;
mod meet;
mod repair;
mod segments;

#[cfg(test)]
pub(crate) mod tests;

use std::collections::HashMap;

use fieldcraft_model::ElementKind;

use crate::delaunay::Triangulation;
use crate::predicates::{GenericPoint, Point3, Sign};
use crate::quality::{GateOptions, GateReport, validity_gate_elements};
use crate::surface::{BoundarySurface, Box3};
use crate::{Error, Result};

use complex::{Complex, NONE, corrupt, face_toward, sorted3};
use segments::Segment;

/// Caps on the work of [`tetrahedralize`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CdtOptions {
    /// Most Steiner points added on segments. A boundary whose triangles are long next to
    /// features much closer than their length (needle triangles) needs Steiner points at that
    /// spacing all along them; the default is well above what the corpus needs (30k at most) and
    /// stops such input in seconds.
    pub max_steiner: usize,
    /// Most tetrahedra (inside and outside, while building).
    pub max_tets: usize,
}

impl Default for CdtOptions {
    fn default() -> Self {
        CdtOptions { max_steiner: 50_000, max_tets: 20_000_000 }
    }
}

/// A Steiner point: `points[a] + t (points[b] − points[a])` of the surface edge `[a, b]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SteinerPoint {
    pub segment: [u32; 2],
    pub t: f64,
}

/// A tetrahedron face on the boundary: tetrahedron `tet`'s face opposite its vertex
/// `opposite` (0..4) lies on surface triangle `triangle`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BoundaryFace {
    pub tet: u32,
    pub opposite: u8,
    pub triangle: u32,
}

impl BoundaryFace {
    /// The face's 1-based number in the model's tet4 face order ([`ElementKind::faces`]): the
    /// face without the opposite vertex.
    pub fn tet4_face(&self) -> Option<u8> {
        let k = usize::from(self.opposite);
        let n = ElementKind::Tet4.faces().iter().position(|f| !f.contains(&k))?;
        u8::try_from(n + 1).ok()
    }
}

/// What the recovery did.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CdtStats {
    /// Tetrahedra of the initial Delaunay tetrahedralization (with the box).
    pub delaunay_tets: usize,
    /// Missing subsegments split.
    pub missing_subsegments: usize,
    /// Cavities refilled by gift-wrapping.
    pub cavities: usize,
    /// Tetrahedra made by gift-wrapping.
    pub gift_wrapped_tets: usize,
    /// Tetrahedra in the box at the end.
    pub box_tets: usize,
    /// Tetrahedra that were flat, inverted or inadmissible on the rounded coordinates and were
    /// removed by flips ([`repair`]).
    pub repaired: usize,
}

/// The constrained Delaunay tetrahedralization of a closed surface (see the module docs).
#[derive(Debug, Clone, PartialEq)]
pub struct TetMesh {
    /// The surface's vertices (same ids as [`BoundarySurface::points`]), then the Steiner
    /// points rounded to doubles.
    pub points: Vec<[f64; 3]>,
    /// The Steiner points, for `points[surface vertices..]`.
    pub steiner: Vec<SteinerPoint>,
    /// The inside tetrahedra, tet4 order (positive volume in exact arithmetic).
    pub tets: Vec<[u32; 4]>,
    /// Every tetrahedron face on the surface, with its input triangle.
    pub boundary: Vec<BoundaryFace>,
    pub stats: CdtStats,
}

impl TetMesh {
    /// The validity gate on the rounded coordinates: every tetrahedron positive, and the volume
    /// equal to `expected_volume` (the surface's, [`BoundarySurface::volume`]) within 1e-9.
    pub fn gate(&self, expected_volume: f64) -> Result<GateReport> {
        let coords = |t: &[u32; 4]| -> Result<[[f64; 3]; 4]> {
            let p = |v: u32| self.points.get(v as usize).copied().ok_or_else(|| corrupt("output vertex id"));
            Ok([p(t[0])?, p(t[1])?, p(t[2])?, p(t[3])?])
        };
        let elements = self.tets.iter().enumerate().map(|(i, t)| Ok((i as u32 + 1, ElementKind::Tet4, coords(t)?))).collect::<Result<Vec<_>>>()?;
        validity_gate_elements(elements, &GateOptions { expected_volume: Some(expected_volume), volume_tolerance: 1e-9 })
    }
}

/// The state of one recovery.
pub(crate) struct Recovery {
    cx: Complex,
    /// The surface's triangles (the facets), vertex ids as in the complex.
    triangles: Vec<[u32; 3]>,
    segments: Vec<Segment>,
    /// Sorted edge → segment.
    edge_segment: HashMap<(u32, u32), usize>,
    /// Recovered facet triangles (sorted vertex ids) → facet.
    constraint: HashMap<[u32; 3], u32>,
    /// The surface's vertices are the complex's first `surface` vertices, then 8 box corners,
    /// then the Steiner points.
    surface: usize,
    /// Steiner points in creation order: (segment, t).
    steiner: Vec<(usize, f64)>,
    options: CdtOptions,
    stats: CdtStats,
}

/// The constrained Delaunay tetrahedralization of `surface` (see the module docs).
pub fn tetrahedralize(surface: &BoundarySurface, options: &CdtOptions) -> Result<TetMesh> {
    let mut rec = Recovery::new(surface, options)?;
    rec.classify_acute()?;
    rec.recover_segments()?;
    for f in 0..rec.triangles.len() {
        rec.recover_facet(f)?;
    }
    rec.verify()?;
    let inside = rec.label()?;
    rec.output(&inside)
}

impl Recovery {
    fn new(surface: &BoundarySurface, options: &CdtOptions) -> Result<Self> {
        let n = surface.points().len();
        let triangles = surface.triangles().to_vec();
        if triangles.is_empty() {
            return Err(Error::Invalid("the surface has no triangles".into()));
        }
        let mut points: Vec<[f64; 3]> = surface.points().iter().map(|p| p.map(|x| x + 0.0)).collect();
        let mut bbox = Box3::EMPTY;
        for p in &points {
            bbox.add_point(*p);
        }
        let size = (0..3).map(|k| bbox.hi[k] - bbox.lo[k]).fold(0.0, f64::max).max(1.0);
        for corner in 0..8 {
            points.push(std::array::from_fn(|k| {
                let x = if corner >> k & 1 == 0 { bbox.lo[k] - size } else { bbox.hi[k] + size };
                crate::predicates::snap_to_domain(x) + 0.0
            }));
        }
        let tri = Triangulation::new(&points)?;
        if tri.dimension() != 3 || (0..points.len() as u32).any(|v| tri.duplicate_of(v).is_some()) {
            return Err(Error::Invalid("the surface's vertices must be distinct".into()));
        }
        let generic = points.iter().map(|&p| Point3::try_from(p).map(GenericPoint::Explicit)).collect::<Result<Vec<_>>>()?;
        let cx = Complex::from_delaunay(&tri, generic, options.max_tets)?;
        let mut edges: Vec<(u32, u32)> =
            triangles.iter().flat_map(|t| (0..3).map(move |k| (t[k].min(t[(k + 1) % 3]), t[k].max(t[(k + 1) % 3])))).collect();
        edges.sort_unstable();
        edges.dedup();
        if let Some(&(a, b)) = edges.iter().find(|(a, b)| a == b || *b as usize >= n) {
            return Err(Error::Invalid(format!("surface edge {a}–{b} is degenerate or has no vertex")));
        }
        let segments: Vec<Segment> = edges.iter().map(|&(a, b)| Segment::new(a, b)).collect();
        let edge_segment = edges.iter().enumerate().map(|(i, &e)| (e, i)).collect();
        let stats = CdtStats { delaunay_tets: cx.live_count(), ..CdtStats::default() };
        Ok(Recovery { cx, triangles, segments, edge_segment, constraint: HashMap::new(), surface: n, steiner: Vec::new(), options: *options, stats })
    }

    /// Every segment is an edge, no edge crosses a facet, every facet has triangles, every
    /// tetrahedron is positive.
    fn verify(&mut self) -> Result<()> {
        let mut scratch = Vec::new();
        for s in 0..self.segments.len() {
            let chain: Vec<u32> = self.segment(s)?.splits.iter().map(|&(_, v)| v).collect();
            for w in chain.windows(2) {
                if let &[u, v] = w
                    && !self.cx.has_edge(u, v, &mut scratch)?
                {
                    let at = self.cx.approx(u)?;
                    return Err(Error::failed(crate::message![
                        "a boundary edge near ",
                        crate::message::point(at),
                        " was lost during facet recovery"
                    ]));
                }
            }
        }
        let mut has_faces = vec![false; self.triangles.len()];
        for &f in self.constraint.values() {
            if let Some(h) = has_faces.get_mut(f as usize) {
                *h = true;
            }
        }
        for f in 0..self.triangles.len() {
            if !has_faces.get(f).copied().unwrap_or(false) || !self.facet_present(f)? {
                return Err(Error::Failed(format!("input triangle {f} is not a union of tetrahedron faces")));
            }
        }
        let live: Vec<[u32; 4]> = self.cx.live_tets().map(|(_, vs)| vs).collect();
        for vs in live {
            if self.cx.volume(vs)? != Sign::Positive {
                let at = self.cx.approx(vs[0])?;
                return Err(Error::failed(crate::message![
                    "constrained tetrahedralization is inconsistent: a tetrahedron near ",
                    crate::message::point(at),
                    " is not positive"
                ]));
            }
        }
        Ok(())
    }

    /// Inside (`true`) or outside, per tetrahedron slot: a flood fill from the box corners that
    /// flips across every facet triangle (Diazzi et al. §4.5), then the exact check that each
    /// facet triangle has its inside tetrahedron on the side its normal points away from.
    fn label(&mut self) -> Result<Vec<bool>> {
        let mut label: Vec<Option<bool>> = vec![None; self.cx.slots()];
        let corners = self.surface..self.surface + 8;
        let corner = self.cx.live_tets().find(|(_, vs)| vs.iter().any(|&v| corners.contains(&(v as usize))));
        let (start, _) = corner.ok_or_else(|| corrupt("no tetrahedron touches the bounding box"))?;
        if let Some(l) = label.get_mut(start as usize) {
            *l = Some(false);
        }
        let mut queue = vec![start];
        while let Some(t) = queue.pop() {
            let inside = label.get(t as usize).copied().flatten().unwrap_or(false);
            let vs = self.cx.tet(t)?;
            for (k, nb) in self.cx.adj(t)?.into_iter().enumerate() {
                let want = inside != self.constraint.contains_key(&sorted3(face_toward(vs, k)));
                match label.get_mut(nb as usize) {
                    None => {} // the hull
                    Some(l @ None) => {
                        *l = Some(want);
                        queue.push(nb);
                    }
                    Some(Some(got)) if *got != want => {
                        let at = self.cx.approx(vs[k])?;
                        return Err(Error::failed(crate::message![
                            "inside/outside labelling is inconsistent near ",
                            crate::message::point(at),
                            ": the surface does not bound a volume"
                        ]));
                    }
                    Some(_) => {}
                }
            }
        }
        let live: Vec<(u32, [u32; 4])> = self.cx.live_tets().collect();
        for (t, vs) in live {
            match label.get(t as usize).copied().flatten() {
                None => return Err(corrupt("a tetrahedron the flood fill did not reach")),
                Some(false) => {}
                Some(true) => {
                    for k in 0..4 {
                        let Some(&f) = self.constraint.get(&sorted3(face_toward(vs, k))) else { continue };
                        let [a, b, c] = *self.triangles.get(f as usize).ok_or_else(|| corrupt("facet id"))?;
                        let x = vs.get(k).copied().unwrap_or(NONE);
                        if self.cx.volume([a, b, c, x])? != Sign::Negative {
                            let at = self.cx.approx(a)?;
                            return Err(Error::failed(crate::message![
                                format!("input triangle {f} near "),
                                crate::message::point(at),
                                " faces into the solid: the surface is inverted"
                            ]));
                        }
                    }
                }
            }
        }
        Ok(label.into_iter().map(|l| l == Some(true)).collect())
    }

    fn output(&self, inside: &[bool]) -> Result<TetMesh> {
        let n = self.surface;
        // Surface vertices keep their ids; box corners disappear; Steiner points follow.
        let total = self.cx.points().len();
        let mut id = vec![NONE; total];
        let mut points = Vec::with_capacity(n + self.steiner.len());
        for v in 0..n {
            if let Some(s) = id.get_mut(v) {
                *s = v as u32;
            }
            points.push(self.cx.approx(v as u32)?);
        }
        let mut steiner = Vec::with_capacity(self.steiner.len());
        for (k, &(s, t)) in self.steiner.iter().enumerate() {
            let v = n + 8 + k;
            if let Some(slot) = id.get_mut(v) {
                *slot = points.len() as u32;
            }
            points.push(self.cx.approx(v as u32)?);
            steiner.push(SteinerPoint { segment: self.segment(s)?.ends, t });
        }
        let mut tets = Vec::new();
        let mut boundary = Vec::new();
        for (t, vs) in self.cx.live_tets() {
            if !inside.get(t as usize).copied().unwrap_or(false) {
                continue;
            }
            let mapped = vs.map(|v| id.get(v as usize).copied().unwrap_or(NONE));
            if mapped.contains(&NONE) {
                return Err(corrupt("an inside tetrahedron uses a bounding-box corner"));
            }
            for k in 0..4 {
                if let Some(&f) = self.constraint.get(&sorted3(face_toward(vs, k))) {
                    boundary.push(BoundaryFace { tet: tets.len() as u32, opposite: k as u8, triangle: f });
                }
            }
            tets.push(mapped);
        }
        let mut stats = self.stats;
        stats.box_tets = self.cx.live_count();
        stats.repaired = repair::repair(&points, &mut tets, &mut boundary, self.options.max_tets)?;
        Ok(TetMesh { points, steiner, tets, boundary, stats })
    }
}
