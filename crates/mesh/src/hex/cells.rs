//! Conforming structured hex meshes of a part partitioned into cells (multi-block meshing).
//!
//! An analyst meshes a part no single sweep fits by splitting it into cells that each are a
//! logical cube (six four-sided faces, possibly composite), meshing every cell as a mapped block
//! and joining the blocks so that they conform: the multi-block grids of J. F. Thompson, B. K.
//! Soni, N. P. Weatherill (eds.), *Handbook of Grid Generation*, CRC Press 1999, ch. 3 and 13.
//! [`hex_mesh_cells`] does it in five steps:
//!
//! 1. **Topology.** Each cell's B-rep is read with smoothly meeting faces joined
//!    ([`super::brep::Brep::composite`]: the halves of a hole, a plane running into a fillet) and
//!    recognized as a logical cube ([`super::mapped::recognize`]). Corners of different cells at
//!    the same point (within the merge tolerance, the one geometric decision) are one vertex;
//!    cell edges between the same two vertices along the same curve are one edge; cell faces
//!    bounded by the same four edges are one face, shared by at most two cells (an interface).
//!    Everything after this works on these identities, never on coordinates.
//! 2. **Interval assignment** over the edges of all cells together ([`super::interval`]): the
//!    four edges along each direction of every cell are equal, so a count is the same on both
//!    sides of every interface. Edges that run away from a curved face (radial lines off a hole,
//!    the diagonal off a fillet) are then graded: the spacing at that end matches the spacing
//!    along the curved face there and grows by at most [`GROWTH`] per element towards the target
//!    size; their counts are reassigned with the grading's goals.
//! 3. **Edges** are divided once, on the geometry, and their nodes belong to the edge.
//! 4. **Faces** get one transfinite grid each (Gordon and Hall), from their four edges and
//!    projected onto the face; an interface's grid is made once and both cells use its nodes.
//! 5. **Cells** fill their interior by trivariate transfinite interpolation (Gordon and Thiel)
//!    of their six face grids, blended by the edges' arc-length distribution.
//!
//! Nodes are merged by identity (a vertex, an edge's k-th node, a face grid's node), so the
//! cells conform exactly: no node is duplicated on an interface and none hangs. The result goes
//! through the same finishing as a single body ([`super::finish`]): the certified validity gate,
//! the tiling and fidelity checks, a closed 2-manifold boundary, and no boundary face lying on
//! another cell (a free internal face, where cells touch without sharing a face).

use std::collections::BTreeMap;

use super::brep::Brep;
use super::geom::{TriangleTree, add, dist, scale};
use super::mapped::{self, Cube, face_ring, other};
use super::{Cells, Ctx, HexMesh, HexMethod, HexOptions, Raw, Surface, finish, project_to_face, project_to_faces};
use crate::message;
use crate::quad::{lattice_index, transfinite_grid, transfinite_volume};
use crate::seeds::EdgeSeed;
use crate::size::{SizeField, divide_polyline, ideal_divisions};
use crate::surface::BoundarySurface;
use crate::vec3::{P, cross, dot, norm, sub};

use crate::{Error, Result};

/// One cell: its name (for messages), its boundary surface (face-tagged) and the analytic
/// surfaces of its faces by tag.
#[derive(Clone, Copy)]
pub struct Cell<'a> {
    pub name: &'a str,
    pub surface: &'a BoundarySurface,
    pub faces: &'a [Option<Surface>],
}

/// A conforming mesh of several cells.
#[derive(Debug, Clone, PartialEq)]
pub struct CellsMesh {
    /// The mesh; boundary faces are tagged with the B-rep faces of their element's cell.
    pub mesh: HexMesh,
    /// The cell of each element (index into the cells given).
    pub element_cell: Vec<u32>,
    /// Each cell's block: divisions along its three logical directions.
    pub blocks: Vec<[usize; 3]>,
    /// Faces shared by two cells.
    pub interfaces: usize,
}

/// The most cells meshed together.
pub const MAX_CELLS: usize = 4096;
/// Spacing growth from one element to the next along a graded edge (1.2: each at most 20 %
/// longer than the one before).
pub const GROWTH: f64 = 1.2;
/// Segments a graded edge's polyline is cut into for the integration of its grading (the size
/// helpers take the size constant along a sample).
const GRADED_SEGMENTS: usize = 64;

/// An edge shared by cells: its end vertices (`a < b`), the cell that divides it and how.
struct Edge {
    a: usize,
    b: usize,
    /// The owning cell, its edge and whether that edge runs from `a` to `b`.
    owner: (usize, usize, bool),
    length: f64,
}

/// A face shared by one or two cells.
struct Face {
    /// The owner's ring of vertices (the order of [`face_ring`]).
    ring: [usize; 4],
    /// The cells, with the direction and side of the face in each; the first owns the grid.
    uses: Vec<(usize, usize, usize)>,
}

/// Mesh `cells` (bodies that together make the part, meeting face to face) into one conforming
/// all-hex mesh, every cell a mapped block. `opts.method` must be [`HexMethod::Auto`] or
/// [`HexMethod::Mapped`]; a cell that is not a logical cube, or cells that meet without sharing
/// whole faces, are [`Error::Unsupported`] with the reason.
pub fn hex_mesh_cells(cells: &[Cell], field: &dyn SizeField, opts: &HexOptions) -> Result<CellsMesh> {
    if !matches!(opts.method, HexMethod::Auto | HexMethod::Mapped) {
        return Err(Error::Invalid(format!("conforming cells are meshed as mapped blocks; method `{}` does not apply", opts.method.name())));
    }
    if !(opts.chordal_tolerance >= 0.0 && opts.chordal_tolerance.is_finite()) {
        return Err(Error::invalid(crate::message![
            "the chordal tolerance must be a non-negative number, got ",
            crate::message::length(opts.chordal_tolerance)
        ]));
    }
    if cells.is_empty() || cells.len() > MAX_CELLS {
        return Err(Error::Invalid(format!("{} cells; 1 to {MAX_CELLS} are meshed together", cells.len())));
    }
    if let Some(c) = cells.iter().find(|c| c.surface.shells().len() != 1) {
        return Err(Error::Unsupported(format!("cell `{}` has {} shells; a cell is one solid without voids", c.name, c.surface.shells().len())));
    }
    let surfaces: Vec<&BoundarySurface> = cells.iter().map(|c| c.surface).collect();
    super::check_estimate(&surfaces, field, opts)?;
    let breps: Vec<Brep> = cells.iter().map(|c| Brep::composite(c.surface, c.faces)).collect::<Result<_>>()?;
    let trees: Vec<TriangleTree> = cells.iter().map(|c| TriangleTree::new(c.surface)).collect::<Result<_>>()?;
    let ctxs: Vec<Ctx> = breps.iter().zip(&trees).map(|(brep, tree)| Ctx::new(brep, tree, field, opts)).collect::<Result<_>>()?;
    let names: Vec<&str> = cells.iter().map(|c| c.name).collect();
    let cubes: Vec<Cube> = breps
        .iter()
        .zip(&names)
        .map(|(b, name)| mapped::recognize(b).map_err(|e| Error::Unsupported(format!("cell `{name}`: {}", super::plain(&e)))))
        .collect::<Result<_>>()?;
    let topo = Topology::new(&breps, &cubes, &names, opts.chordal_tolerance)?;
    let divisions = topo.intervals(&ctxs, &cubes, field)?;
    let mut builder =
        Builder { topo: &topo, ctxs: &ctxs, cubes: &cubes, divisions: &divisions, field, opts, nodes: Vec::new(), edge_nodes: Vec::new() };
    let (raw, element_cell, blocks) = builder.build()?;
    let shared = topo.interfaces();
    let mesh = finish(&Cells { ctxs: &ctxs, of: &element_cell, names: &names }, raw, HexMethod::Mapped)?;
    Ok(CellsMesh { mesh, element_cell, blocks, interfaces: shared })
}

/// Mesh `cells` (bodies that together make a part, meeting face to face) into one conforming
/// all-hex mesh by sweeping them together, along one direction or about one axis (`opts.method`:
/// [`HexMethod::Extrude`], [`HexMethod::Revolve`], or [`HexMethod::Auto`] / [`HexMethod::Swept`]
/// for the first that fits).
///
/// This is many-to-many sweeping ([`super::layered`]) over the cells at once: the walls of every
/// cell (its cut faces included) are imprinted into one cross-section and the levels of every cell
/// are levels of all, so one interval assignment gives every curve and every level interval one
/// count, and each region of each interval is filled by the cell it lies in. Where two cells meet
/// they share the curves and layers there, so their elements meet node for node with no seeds;
/// a cut that ends on another cut (a T-junction) is a curve ending on a curve like any other.
/// Each element belongs to one cell and boundary faces are tagged with that cell's faces. Cells
/// that do not all sweep in one frame are [`Error::Unsupported`] with the first reason.
pub fn hex_mesh_swept_cells(cells: &[Cell], field: &dyn SizeField, opts: &HexOptions) -> Result<CellsMesh> {
    let methods: &[HexMethod] = match opts.method {
        HexMethod::Auto | HexMethod::Swept => &[HexMethod::Extrude, HexMethod::Revolve],
        HexMethod::Extrude => &[HexMethod::Extrude],
        HexMethod::Revolve => &[HexMethod::Revolve],
        m => return Err(Error::Invalid(format!("cells are swept together by extrusion or revolution; method `{}` does not apply", m.name()))),
    };
    if !(opts.chordal_tolerance >= 0.0 && opts.chordal_tolerance.is_finite()) {
        return Err(Error::invalid(crate::message![
            "the chordal tolerance must be a non-negative number, got ",
            crate::message::length(opts.chordal_tolerance)
        ]));
    }
    if cells.is_empty() || cells.len() > MAX_CELLS {
        return Err(Error::Invalid(format!("{} cells; 1 to {MAX_CELLS} are meshed together", cells.len())));
    }
    if let Some(c) = cells.iter().find(|c| c.surface.shells().len() != 1) {
        return Err(Error::Unsupported(format!("cell `{}` has {} shells; a cell is one solid without voids", c.name, c.surface.shells().len())));
    }
    let surfaces: Vec<&BoundarySurface> = cells.iter().map(|c| c.surface).collect();
    super::check_estimate(&surfaces, field, opts)?;
    let breps: Vec<Brep> = cells.iter().map(|c| Brep::new(c.surface, c.faces)).collect::<Result<_>>()?;
    let trees: Vec<TriangleTree> = cells.iter().map(|c| TriangleTree::new(c.surface)).collect::<Result<_>>()?;
    // One classification tolerance for all, from their joint extent.
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for p in cells.iter().flat_map(|c| c.surface.points()) {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let tol = super::RELATIVE_TOLERANCE * dist(lo, hi);
    if !(tol > 0.0 && tol.is_finite()) {
        return Err(Error::Invalid("the cells have no extent".into()));
    }
    let ctxs: Vec<Ctx> = breps.iter().zip(&trees).map(|(brep, tree)| Ctx::with_tolerance(brep, tree, field, opts, tol)).collect::<Result<_>>()?;
    let names: Vec<&str> = cells.iter().map(|c| c.name).collect();
    let mut reasons: Vec<String> = Vec::new();
    for &method in methods {
        let setups = match super::candidate_setups(&ctxs, method) {
            Ok(s) => s,
            Err(why) => {
                reasons.push(format!("{}: {why}", method.name()));
                continue;
            }
        };
        let mut last = None;
        for setup in setups {
            fieldcraft_model::progress::check()?;
            let prepared = match super::layered::prepare(&ctxs, setup) {
                Ok(p) => p,
                Err(e) if e.kind() == Some(crate::Kind::Limit) => return Err(e),
                Err(e) => {
                    last = Some(e);
                    continue;
                }
            };
            // As a single body's sweep: on a face mesher's failure, again with finer curves.
            for refine in 0..=super::MAX_REFINE {
                fieldcraft_model::progress::check()?;
                let made = super::layered::mesh(&ctxs, &prepared, refine)
                    .and_then(|(raw, of)| Ok((finish(&Cells { ctxs: &ctxs, of: &of, names: &names }, raw, method)?, of)));
                match made {
                    Ok((mesh, element_cell)) => {
                        let interfaces = shared_faces(&mesh, &element_cell);
                        return Ok(CellsMesh { mesh, element_cell, blocks: Vec::new(), interfaces });
                    }
                    Err(e) if e.kind() == Some(crate::Kind::Limit) => return Err(e),
                    Err(e) if e.kind() == Some(crate::Kind::Failed) => last = Some(e),
                    Err(e) => {
                        last = Some(e);
                        break;
                    }
                }
            }
        }
        if let Some(e) = last {
            reasons.push(format!("{}: {}", method.name(), super::plain(&e)));
        }
    }
    Err(Error::Unsupported(format!("the cells do not sweep together ({})", reasons.join("; "))))
}

/// The pairs of cells that share element faces in `mesh`, counted once each.
fn shared_faces(mesh: &HexMesh, element_cell: &[u32]) -> usize {
    let mut owner: std::collections::HashMap<[u32; 4], u32> = std::collections::HashMap::new();
    let mut pairs = std::collections::BTreeSet::new();
    for e in 0..mesh.element_count() {
        let (Some(h), Some(&c)) = (mesh.element(e), element_cell.get(e)) else { continue };
        for k in 0..6 {
            let mut key = super::hex_face(k).map(|i| h.get(i).copied().unwrap_or(0));
            key.sort_unstable();
            match owner.insert(key, c) {
                Some(o) if o != c => {
                    pairs.insert((o.min(c), o.max(c)));
                }
                _ => {}
            }
        }
    }
    pairs.len()
}

/// Vertices, edges and faces shared by the cells.
struct Topology {
    /// Shared vertices: where they are.
    vertices: Vec<P>,
    edges: Vec<Edge>,
    faces: Vec<Face>,
    /// Per cell: the vertex at each corner, `corner[i][j][k]`.
    corners: Vec<[[[usize; 2]; 2]; 2]>,
    /// Per cell, per direction: the four edges (in the order of [`Cube::groups`]).
    cell_edges: Vec<[[usize; 4]; 3]>,
    /// Per cell: the face on each side of each direction.
    cell_faces: Vec<[[usize; 2]; 3]>,
}

/// Points at most this far apart (relative to the diagonal of all cells) are one vertex.
const RELATIVE_MERGE: f64 = 1e-6;

impl Topology {
    fn new(breps: &[Brep], cubes: &[Cube], names: &[&str], chordal: f64) -> Result<Topology> {
        let name = |c: usize| names.get(c).copied().unwrap_or("?");
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for b in breps {
            for p in b.surface.points() {
                for k in 0..3 {
                    lo[k] = lo[k].min(p[k]);
                    hi[k] = hi[k].max(p[k]);
                }
            }
        }
        let merge = RELATIVE_MERGE * dist(lo, hi);
        if !(merge > 0.0 && merge.is_finite()) {
            return Err(Error::Invalid("the cells have no extent".into()));
        }
        // Vertices: corners within `merge` of each other, found through a grid of that spacing.
        let mut vertices: Vec<P> = Vec::new();
        let mut grid: BTreeMap<[i64; 3], Vec<usize>> = BTreeMap::new();
        let cellkey = |p: P| p.map(|x| (x / merge).floor() as i64);
        let mut corners = Vec::with_capacity(cubes.len());
        for (c, (brep, cube)) in breps.iter().zip(cubes).enumerate() {
            let mut ids = [[[0usize; 2]; 2]; 2];
            let mut seen: Vec<usize> = Vec::with_capacity(8);
            for (i, plane) in ids.iter_mut().enumerate() {
                for (j, row) in plane.iter_mut().enumerate() {
                    for (k, slot) in row.iter_mut().enumerate() {
                        let p = brep.point(cube.at(i, j, k));
                        let key = cellkey(p);
                        let mut found: Option<usize> = None;
                        for dx in -1..=1 {
                            for dy in -1..=1 {
                                for dz in -1..=1 {
                                    let near = [key[0] + dx, key[1] + dy, key[2] + dz];
                                    for &v in grid.get(&near).map(Vec::as_slice).unwrap_or(&[]) {
                                        if vertices.get(v).is_some_and(|&x| dist(x, p) <= merge) {
                                            if found.is_some_and(|f| f != v) {
                                                return Err(Error::Measured(
                                                    crate::Kind::Unsupported,
                                                    crate::message![
                                                        format!("cell `{}`: a corner at ", name(c)),
                                                        crate::message::point(p),
                                                        " is within the merge distance of two vertices"
                                                    ],
                                                ));
                                            }
                                            found = Some(v);
                                        }
                                    }
                                }
                            }
                        }
                        let v = match found {
                            Some(v) => v,
                            None => {
                                vertices.push(p);
                                grid.entry(key).or_default().push(vertices.len() - 1);
                                vertices.len() - 1
                            }
                        };
                        if seen.contains(&v) {
                            return Err(Error::Measured(
                                crate::Kind::Unsupported,
                                crate::message![format!("cell `{}`: two corners at the same point ", name(c)), crate::message::point(p)],
                            ));
                        }
                        seen.push(v);
                        *slot = v;
                    }
                }
            }
            corners.push(ids);
        }
        // Edges: the same two vertices and the same curve (midpoints within the merge distance
        // and the tessellation's chordal deviation).
        let edge_tol = merge + 4.0 * chordal;
        let mut edges: Vec<Edge> = Vec::new();
        let mut by_ends: BTreeMap<(usize, usize), Vec<(usize, P)>> = BTreeMap::new();
        let mut cell_edges = Vec::with_capacity(cubes.len());
        for (c, (brep, cube)) in breps.iter().zip(cubes).enumerate() {
            let ids = corners.get(c).copied().ok_or_else(|| Error::Failed("cell corners".into()))?;
            let mut out = [[0usize; 4]; 3];
            for (d, row) in out.iter_mut().enumerate() {
                for (m, slot) in row.iter_mut().enumerate() {
                    let (u, v) = edge_ends(&ids, d, m);
                    let (e, fwd) = cube.groups.get(d).and_then(|g| g.get(m)).copied().ok_or_else(|| Error::Failed("cube groups".into()))?;
                    let pts = brep.edge_points(e, fwd);
                    let (mid, length) = midpoint(&pts);
                    let (a, b) = (u.min(v), u.max(v));
                    let list = by_ends.entry((a, b)).or_default();
                    let id = match list.iter().find(|(_, q)| dist(*q, mid) <= edge_tol) {
                        Some(&(id, _)) => id,
                        None => {
                            edges.push(Edge { a, b, owner: (c, e, fwd == (u == a)), length });
                            list.push((edges.len() - 1, mid));
                            edges.len() - 1
                        }
                    };
                    *slot = id;
                }
            }
            cell_edges.push(out);
        }
        // Faces: the same four edges, on the same surface, with the cells on opposite sides.
        let mut faces: Vec<Face> = Vec::new();
        let mut by_edges: BTreeMap<[usize; 4], usize> = BTreeMap::new();
        let mut cell_faces = Vec::with_capacity(cubes.len());
        for (c, (brep, cube)) in breps.iter().zip(cubes).enumerate() {
            let ids = corners.get(c).copied().ok_or_else(|| Error::Failed("cell corners".into()))?;
            let mut out = [[0usize; 2]; 3];
            for (d, row) in out.iter_mut().enumerate() {
                for (side, slot) in row.iter_mut().enumerate() {
                    let ring = face_ring(d, side).map(|[i, j, k]| ids[i][j][k]);
                    let mut key = [0usize; 4];
                    for (m, k) in key.iter_mut().enumerate() {
                        let (u, v) = (ring[m], ring[(m + 1) % 4]);
                        let groups = cell_edges.get(c).ok_or_else(|| Error::Failed("cell edges".into()))?;
                        *k = edge_between(&edges, groups, u, v).ok_or_else(|| Error::Failed("a face side is no cell edge".into()))?;
                    }
                    key.sort_unstable();
                    let face = cube.faces.get(d).and_then(|r| r.get(side)).copied().ok_or_else(|| Error::Failed("cube faces".into()))?;
                    let id = match by_edges.get(&key) {
                        Some(&id) => {
                            let f = faces.get_mut(id).ok_or_else(|| Error::Failed("face index".into()))?;
                            if f.uses.len() >= 2 {
                                return Err(Error::Unsupported(format!("cell `{}`: a face is shared by more than two cells", name(c))));
                            }
                            if let Some(&(c0, d0, s0)) = f.uses.first() {
                                if c0 == c {
                                    return Err(Error::Unsupported(format!("cell `{}`: two faces of the cell have the same four edges", name(c))));
                                }
                                let f0 = cubes.get(c0).and_then(|q| q.faces.get(d0)).and_then(|r| r.get(s0)).copied().unwrap_or(0);
                                let b0 = breps.get(c0).ok_or_else(|| Error::Failed("cell index".into()))?;
                                let (centre0, area0) = face_centre(b0, f0);
                                let (centre1, area1) = face_centre(brep, face);
                                if dist(centre0, centre1) > edge_tol + 1e-3 * norm(area0).sqrt() || dot(area0, area1) >= 0.0 {
                                    return Err(Error::Unsupported(format!(
                                        "cells `{}` and `{}` share the four edges of a face but not the face (different surfaces, or on the same side)",
                                        name(c0),
                                        name(c)
                                    )));
                                }
                            }
                            f.uses.push((c, d, side));
                            id
                        }
                        None => {
                            faces.push(Face { ring, uses: vec![(c, d, side)] });
                            by_edges.insert(key, faces.len() - 1);
                            faces.len() - 1
                        }
                    };
                    *slot = id;
                }
            }
            cell_faces.push(out);
        }
        Ok(Topology { vertices, edges, faces, corners, cell_edges, cell_faces })
    }

    /// Faces shared by two cells.
    fn interfaces(&self) -> usize {
        self.faces.iter().filter(|f| f.uses.len() == 2).count()
    }

    /// The edge between vertices `u` and `v` among cell `c`'s edges, and whether it runs from
    /// `u` to `v`.
    fn edge(&self, c: usize, u: usize, v: usize) -> Option<(usize, bool)> {
        let e = edge_between(&self.edges, self.cell_edges.get(c)?, u, v)?;
        Some((e, self.edges.get(e)?.a == u))
    }

    /// How every edge is divided, all cells together; graded edges get the counts their
    /// grading asks for.
    fn intervals(&self, ctxs: &[Ctx], cubes: &[Cube], field: &dyn SizeField) -> Result<Vec<Division>> {
        let polylines: Vec<Vec<P>> = self.edges.iter().map(|e| self.polyline(ctxs, e)).collect::<Result<_>>()?;
        let goals: Vec<f64> = polylines.iter().map(|p| ideal_divisions(field, p)).collect::<Result<_>>()?;
        // The seed on each edge (a shared edge is one edge, so its seed holds in every cell it
        // bounds), oriented from `a` to `b`.
        let seeds: Vec<Option<(&EdgeSeed, bool)>> = self
            .edges
            .iter()
            .map(|e| {
                let (c, be, fwd) = e.owner;
                ctxs.get(c).and_then(|x| x.seed_of(be)).map(|(seed, rev)| (seed, rev == fwd))
            })
            .collect();
        // A biased seed's spacing carries to the unseeded edges parallel to it, through every
        // cell they bound, oriented along each cell's direction (#257); the counts come from the
        // edges' own seeds alone.
        let carried = crate::seeds::carry(&seeds, |s: &&EdgeSeed| s.bias.is_some(), &self.parallel_ties());
        let first = self.solve(&goals, &seeds)?;
        // Ends of edges that leave a curved face: the spacing along that face there.
        let ends = self.graded_ends(ctxs, cubes, &first);
        let mut graded = goals;
        let mut divisions = Vec::with_capacity(polylines.len());
        let spacing: Vec<Option<(&EdgeSeed, bool)>> = seeds.iter().zip(carried).map(|(own, carried)| own.or(carried)).collect();
        for (((poly, ends), goal), seed) in polylines.into_iter().zip(ends).zip(graded.iter_mut()).zip(&spacing) {
            // A seeded edge keeps the seed's spacing, not a grading.
            let ends = if seed.is_some() { [None, None] } else { ends };
            let poly = if ends.iter().any(Option::is_some) {
                let dense = densify(&poly);
                *goal = goal.max(ideal_divisions(&Grading { field, ends }, &dense)?);
                dense
            } else {
                poly
            };
            divisions.push(Division { count: 1, ends, polyline: poly, seed: seed.map(|(s, rev)| (s.clone(), rev)) });
        }
        for (d, n) in divisions.iter_mut().zip(self.solve(&graded, &seeds)?) {
            d.count = n;
        }
        Ok(divisions)
    }

    /// The four edges along each direction of every cell, tied pairwise, each traversed from the
    /// cell's low corner to its high one.
    fn parallel_ties(&self) -> Vec<[crate::seeds::Use; 2]> {
        let mut ties = Vec::new();
        for (groups, ids) in self.cell_edges.iter().zip(&self.corners) {
            for (d, g) in groups.iter().enumerate() {
                let uses: Vec<crate::seeds::Use> = g
                    .iter()
                    .enumerate()
                    .map(|(m, &e)| {
                        let (lo, _) = edge_ends(ids, d, m);
                        (e, self.edges.get(e).is_some_and(|x| x.a == lo))
                    })
                    .collect();
                if let Some((&first, rest)) = uses.split_first() {
                    ties.extend(rest.iter().map(|&u| [first, u]));
                }
            }
        }
        ties
    }

    fn solve(&self, goals: &[f64], seeds: &[Option<(&EdgeSeed, bool)>]) -> Result<Vec<usize>> {
        let mut ia = super::interval::Intervals::default();
        for (&g, seed) in goals.iter().zip(seeds) {
            let v = ia.var(g, 1);
            if let Some((seed, _)) = seed {
                ia.seed(v, seed);
            }
        }
        for groups in &self.cell_edges {
            for g in groups {
                for &e in g.iter().skip(1) {
                    ia.equal.push((vec![g[0]], vec![e]));
                }
            }
        }
        ia.solve()
    }

    /// For every edge, at each end (`a`, `b`): the point and the spacing to start from when the
    /// edge leaves a curved face there (ends on it without lying on it) — the smallest spacing,
    /// by `counts`, of that face's edges at the end — or `None`.
    fn graded_ends(&self, ctxs: &[Ctx], cubes: &[Cube], counts: &[usize]) -> Vec<GradedEnds> {
        // The edges of every curved face, by the vertices they meet.
        let mut curved: Vec<[usize; 4]> = Vec::new();
        let mut curved_at: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for f in &self.faces {
            let Some(&(c, d, side)) = f.uses.first() else { continue };
            let planar = cubes
                .get(c)
                .and_then(|q| q.faces.get(d))
                .and_then(|r| r.get(side))
                .and_then(|&bf| ctxs.get(c).and_then(|x| x.brep.faces.get(bf)))
                .is_none_or(|bf| bf.plane.is_some());
            if planar {
                continue;
            }
            let mut es = [usize::MAX; 4];
            for (m, slot) in es.iter_mut().enumerate() {
                if let Some((e, _)) = self.edge(c, f.ring[m], f.ring[(m + 1) % 4]) {
                    *slot = e;
                }
            }
            for &v in &f.ring {
                curved_at.entry(v).or_default().push(curved.len());
            }
            curved.push(es);
        }
        let spacing = |e: usize| -> f64 { self.edges.get(e).map_or(f64::INFINITY, |x| x.length / counts.get(e).copied().unwrap_or(1).max(1) as f64) };
        self.edges
            .iter()
            .enumerate()
            .map(|(e, edge)| {
                [edge.a, edge.b].map(|v| {
                    let mut best: Option<f64> = None;
                    for face in curved_at.get(&v).map(Vec::as_slice).unwrap_or(&[]).iter().filter_map(|&k| curved.get(k)) {
                        if face.contains(&e) {
                            continue;
                        }
                        for &o in face {
                            if self.edges.get(o).is_some_and(|x| x.a == v || x.b == v) {
                                let s = spacing(o);
                                best = Some(best.map_or(s, |b| b.min(s)));
                            }
                        }
                    }
                    best.and_then(|s| self.vertices.get(v).map(|&x| (x, s)))
                })
            })
            .collect()
    }

    /// The owner's polyline of an edge, from `a` to `b`.
    fn polyline(&self, ctxs: &[Ctx], e: &Edge) -> Result<Vec<P>> {
        let (c, be, fwd) = e.owner;
        let ctx = ctxs.get(c).ok_or_else(|| Error::Failed("edge owner".into()))?;
        let mut pts = ctx.brep.edge_points(be, fwd);
        // The ends are the shared vertices.
        if let (Some(first), Some(a)) = (pts.first_mut(), self.vertices.get(e.a)) {
            *first = *a;
        }
        if let (Some(last), Some(b)) = (pts.last_mut(), self.vertices.get(e.b)) {
            *last = *b;
        }
        if pts.len() < 2 {
            return Err(Error::Failed("an edge without points".into()));
        }
        Ok(pts)
    }
}

/// The edge among a cell's edges (`groups`) between vertices `u` and `v`.
fn edge_between(edges: &[Edge], groups: &[[usize; 4]; 3], u: usize, v: usize) -> Option<usize> {
    let key = (u.min(v), u.max(v));
    groups.iter().flatten().copied().find(|&e| edges.get(e).is_some_and(|x| (x.a, x.b) == key))
}

/// The corners (low, high) of edge `m` of direction `d` of a cell whose corners are `ids`, in
/// the order of [`Cube::groups`].
fn edge_ends(ids: &[[[usize; 2]; 2]; 2], d: usize, m: usize) -> (usize, usize) {
    let (s, t) = ((m / 2) & 1, m & 1);
    match d {
        0 => (ids[0][s][t], ids[1][s][t]),
        1 => (ids[s][0][t], ids[s][1][t]),
        _ => (ids[s][t][0], ids[s][t][1]),
    }
}

/// The point halfway along a polyline (by arc length) and its length.
fn midpoint(pts: &[P]) -> (P, f64) {
    let length: f64 = pts.windows(2).map(|w| dist(w[0], w[1])).sum();
    let mut walked = 0.0;
    for w in pts.windows(2) {
        let l = dist(w[0], w[1]);
        if walked + l >= 0.5 * length && l > 0.0 {
            let t = (0.5 * length - walked) / l;
            return (add(w[0], scale(sub(w[1], w[0]), t)), length);
        }
        walked += l;
    }
    (pts.first().copied().unwrap_or([0.0; 3]), length)
}

/// A B-rep face's area-weighted centre and its area vector (outward).
fn face_centre(brep: &Brep, f: usize) -> (P, P) {
    let (mut c, mut n, mut total) = ([0.0; 3], [0.0; 3], 0.0);
    for &t in brep.faces.get(f).map(|x| x.triangles.as_slice()).unwrap_or(&[]) {
        let Some(tri) = brep.surface.triangles().get(t as usize) else { continue };
        let [a, b, d] = tri.map(|v| brep.point(v));
        let area = cross(sub(b, a), sub(d, a));
        let w = 0.5 * norm(area);
        c = add(c, scale(add(add(a, b), d), w / 3.0));
        n = add(n, scale(area, 0.5));
        total += w;
    }
    (if total > 0.0 { scale(c, 1.0 / total) } else { c }, n)
}

/// A polyline cut into at least [`GRADED_SEGMENTS`] pieces (for the integration of a grading).
fn densify(pts: &[P]) -> Vec<P> {
    let length: f64 = pts.windows(2).map(|w| dist(w[0], w[1])).sum();
    let step = length / GRADED_SEGMENTS as f64;
    let mut out = Vec::with_capacity(pts.len() + GRADED_SEGMENTS);
    for w in pts.windows(2) {
        let l = dist(w[0], w[1]);
        let k = if step > 0.0 { ((l / step).ceil() as usize).clamp(1, GRADED_SEGMENTS) } else { 1 };
        for s in 0..k {
            out.push(add(w[0], scale(sub(w[1], w[0]), s as f64 / k as f64)));
        }
    }
    if let Some(&last) = pts.last() {
        out.push(last);
    }
    out
}

/// Each end of an edge (`a`, `b`): where grading starts from and the spacing there.
type GradedEnds = [Option<(P, f64)>; 2];

/// How an edge is divided: its count, its graded ends and the polyline (from `a` to `b`; cut
/// finer when graded) it is divided along.
struct Division {
    count: usize,
    ends: GradedEnds,
    polyline: Vec<P>,
    /// The edge's seed, and whether it runs from `b` to `a`.
    seed: Option<(EdgeSeed, bool)>,
}

/// The size along a graded edge: the field's, capped near each graded end by the spacing there
/// growing by [`GROWTH`] per element (linearly with the distance from that end).
struct Grading<'a> {
    field: &'a dyn SizeField,
    ends: GradedEnds,
}

impl SizeField for Grading<'_> {
    fn size(&self, p: [f64; 3]) -> f64 {
        let mut h = self.field.size(p);
        for (q, s) in self.ends.iter().flatten() {
            h = h.min(s + (GROWTH - 1.0) * dist(p, *q));
        }
        h
    }
}

/// The nodes of a cell's lattice, or of a face's grid, by identity.
struct Builder<'t, 'c, 'a> {
    topo: &'t Topology,
    ctxs: &'c [Ctx<'a>],
    cubes: &'c [Cube],
    divisions: &'c [Division],
    field: &'c dyn SizeField,
    opts: &'c HexOptions,
    nodes: Vec<P>,
    /// Per edge: its node ids from `a` to `b`.
    edge_nodes: Vec<Vec<u32>>,
}

impl Builder<'_, '_, '_> {
    fn count(&self, e: usize) -> usize {
        self.divisions.get(e).map_or(1, |d| d.count)
    }

    fn push(&mut self, p: P) -> Result<u32> {
        let id = u32::try_from(self.nodes.len()).map_err(|_| Error::Limit("too many nodes".into()))?;
        self.nodes.push(p);
        Ok(id)
    }

    fn build(&mut self) -> Result<(Raw, Vec<u32>, Vec<[usize; 3]>)> {
        let (topo, field, max_elements) = (self.topo, self.field, self.opts.max_elements);
        // Block sizes and the element budget before any node is made.
        let mut blocks = Vec::with_capacity(topo.cell_edges.len());
        let mut total = 0usize;
        for groups in &topo.cell_edges {
            let n = groups.map(|g| self.count(g[0]));
            let nodes = (n[0] + 1).checked_mul(n[1] + 1).and_then(|x| x.checked_mul(n[2] + 1)).filter(|&x| x <= mapped::MAX_NODES);
            if nodes.is_none() {
                return Err(Error::Limit(format!("a {} × {} × {} block exceeds {} nodes", n[0], n[1], n[2], mapped::MAX_NODES)));
            }
            total = total.saturating_add(n[0] * n[1] * n[2]);
            blocks.push(n);
        }
        if total > max_elements {
            return Err(Error::Limit(format!("the mesh would have {total} hexahedra; the limit is {max_elements}")));
        }
        // Vertices.
        for &v in &topo.vertices {
            self.push(v)?;
        }
        // Edges: divided once, on the geometry.
        for (edge, division) in topo.edges.iter().zip(self.divisions) {
            let n = division.count;
            let pts = if let Some((seed, rev)) = &division.seed {
                crate::seeds::place(&division.polyline, &seed.fractions(n, *rev))?
            } else if division.ends.iter().any(Option::is_some) {
                divide_polyline(&Grading { field, ends: division.ends }, &division.polyline, n)?
            } else {
                divide_polyline(field, &division.polyline, n)?
            };
            let (c, be, _) = edge.owner;
            let ctx = self.ctxs.get(c).ok_or_else(|| Error::Failed("edge owner".into()))?;
            let faces = ctx.brep.edges.get(be).map(|x| x.faces.to_vec()).unwrap_or_default();
            let mut ids = Vec::with_capacity(n + 1);
            ids.push(edge.a as u32);
            for q in pts.iter().take(n).skip(1) {
                let p = project_to_faces(ctx, &faces, *q);
                ids.push(self.push(p)?);
            }
            ids.push(edge.b as u32);
            self.edge_nodes.push(ids);
        }
        // Faces: one grid each, from the owner's side.
        let mut face_grids: Vec<FaceGrid> = Vec::with_capacity(topo.faces.len());
        for face in &topo.faces {
            face_grids.push(self.face_grid(face)?);
        }
        // Cells.
        let mut hexes: Vec<[u32; 8]> = Vec::with_capacity(total);
        let mut element_cell: Vec<u32> = Vec::with_capacity(total);
        for (c, &n) in blocks.iter().enumerate() {
            let ids = self.cell_lattice(c, n, &face_grids)?;
            let [ni, nj, nk] = n;
            let id = |i: usize, j: usize, k: usize| ids.get(lattice_index(n, [i, j, k])).copied().unwrap_or(0);
            for k in 0..nk {
                for j in 0..nj {
                    for i in 0..ni {
                        hexes.push([
                            id(i, j, k),
                            id(i + 1, j, k),
                            id(i + 1, j + 1, k),
                            id(i, j + 1, k),
                            id(i, j, k + 1),
                            id(i + 1, j, k + 1),
                            id(i + 1, j + 1, k + 1),
                            id(i, j + 1, k + 1),
                        ]);
                        element_cell.push(c as u32);
                    }
                }
            }
        }
        let shared = topo.interfaces();
        let description = message![format!(
            "{} cell{} as conforming mapped blocks, {} shared face{}",
            blocks.len(),
            if blocks.len() == 1 { "" } else { "s" },
            shared,
            if shared == 1 { "" } else { "s" }
        )];
        let nodes = std::mem::take(&mut self.nodes);
        Ok((Raw { nodes, hexes, description }, element_cell, blocks))
    }

    /// The nodes of edge `e` from vertex `u`.
    fn edge_from(&self, e: usize, u: usize) -> Result<Vec<u32>> {
        let edge = self.topo.edges.get(e).ok_or_else(|| Error::Failed("edge index".into()))?;
        let mut ids = self.edge_nodes.get(e).cloned().ok_or_else(|| Error::Failed("edge nodes".into()))?;
        if edge.b == u {
            ids.reverse();
        } else if edge.a != u {
            return Err(Error::Failed("an edge does not start at the vertex asked".into()));
        }
        Ok(ids)
    }

    /// A face's transfinite grid in its owner's frame.
    fn face_grid(&mut self, face: &Face) -> Result<FaceGrid> {
        let &(c, d, side) = face.uses.first().ok_or_else(|| Error::Failed("a face without a cell".into()))?;
        let r = face.ring;
        let side_nodes = |b: &Self, u: usize, v: usize| -> Result<Vec<u32>> {
            let (e, _) = b.topo.edge(c, u, v).ok_or_else(|| Error::Failed("a face side is no edge".into()))?;
            b.edge_from(e, u)
        };
        let bottom = side_nodes(self, r[0], r[1])?;
        let right = side_nodes(self, r[1], r[2])?;
        let top = side_nodes(self, r[3], r[2])?;
        let left = side_nodes(self, r[0], r[3])?;
        let (na, nb) = (bottom.len() - 1, left.len() - 1);
        if top.len() != na + 1 || right.len() != nb + 1 {
            return Err(Error::Failed("opposite sides of a face have different counts".into()));
        }
        let at = |b: &Self, ids: &[u32]| -> Vec<P> { ids.iter().map(|&i| b.nodes.get(i as usize).copied().unwrap_or([0.0; 3])).collect() };
        let grid = transfinite_grid(&at(self, &bottom), &at(self, &right), &at(self, &top), &at(self, &left))?;
        let ctx = self.ctxs.get(c).ok_or_else(|| Error::Failed("face owner".into()))?;
        let brep_face = self.cubes.get(c).and_then(|q| q.faces.get(d)).and_then(|row| row.get(side)).copied();
        let mut ids = vec![0u32; (na + 1) * (nb + 1)];
        for b in 0..=nb {
            for a in 0..=na {
                let id = if b == 0 {
                    bottom.get(a).copied()
                } else if b == nb {
                    top.get(a).copied()
                } else if a == 0 {
                    left.get(b).copied()
                } else if a == na {
                    right.get(b).copied()
                } else {
                    let q = grid.get(b * (na + 1) + a).copied().ok_or_else(|| Error::Failed("face grid".into()))?;
                    let q = brep_face.and_then(|f| project_to_face(ctx, f, q)).unwrap_or(q);
                    Some(self.push(q)?)
                };
                if let Some(slot) = ids.get_mut(b * (na + 1) + a) {
                    *slot = id.ok_or_else(|| Error::Failed("face side node".into()))?;
                }
            }
        }
        Ok(FaceGrid { ring: r, na, nb, ids })
    }

    /// Cell `c`'s lattice of node ids: boundary from the shared face grids, interior new.
    fn cell_lattice(&mut self, c: usize, n: [usize; 3], face_grids: &[FaceGrid]) -> Result<Vec<u32>> {
        let topo = self.topo;
        let corners = topo.corners.get(c).copied().ok_or_else(|| Error::Failed("cell corners".into()))?;
        let faces = topo.cell_faces.get(c).copied().ok_or_else(|| Error::Failed("cell faces".into()))?;
        let total = (n[0] + 1) * (n[1] + 1) * (n[2] + 1);
        let mut ids: Vec<Option<u32>> = vec![None; total];
        for d in 0..3 {
            let (da, db) = other(d);
            for side in 0..2 {
                let gf = faces.get(d).and_then(|r| r.get(side)).copied().ok_or_else(|| Error::Failed("cell face".into()))?;
                let grid = face_grids.get(gf).ok_or_else(|| Error::Failed("face grid".into()))?;
                let ring = face_ring(d, side).map(|[i, j, k]| corners[i][j][k]);
                let map = grid.map_from(ring, n[da], n[db])?;
                for b in 0..=n[db] {
                    for a in 0..=n[da] {
                        let mut l = [0usize; 3];
                        l[d] = side * n[d];
                        l[da] = a;
                        l[db] = b;
                        if let Some(slot) = ids.get_mut(lattice_index(n, l)) {
                            *slot = Some(map(a, b)?);
                        }
                    }
                }
            }
        }
        // Interior: trivariate interpolation, blended by the edges' arc-length distribution.
        let mut x: Vec<Option<P>> = ids.iter().map(|i| i.and_then(|i| self.nodes.get(i as usize).copied())).collect();
        let groups = topo.cell_edges.get(c).copied().ok_or_else(|| Error::Failed("cell edges".into()))?;
        let mut blend: [Vec<f64>; 3] = [Vec::new(), Vec::new(), Vec::new()];
        for (d, g) in groups.iter().enumerate() {
            let mut sum = vec![0.0; n[d] + 1];
            for (m, &e) in g.iter().enumerate() {
                let pts: Vec<P> = self
                    .edge_from(e, edge_ends(&corners, d, m).0)?
                    .iter()
                    .map(|&i| self.nodes.get(i as usize).copied().unwrap_or([0.0; 3]))
                    .collect();
                let params = crate::quad::tfi::params(&pts, "edge")?;
                if params.len() != sum.len() {
                    return Err(Error::Failed("an edge's count differs from its block's".into()));
                }
                for (acc, p) in sum.iter_mut().zip(params) {
                    *acc += 0.25 * p;
                }
            }
            if let Some(first) = sum.first_mut() {
                *first = 0.0;
            }
            if let Some(last) = sum.last_mut() {
                *last = 1.0;
            }
            blend[d] = sum;
        }
        transfinite_volume(&mut x, n, [&blend[0], &blend[1], &blend[2]])?;
        let mut out = Vec::with_capacity(total);
        for (slot, p) in ids.iter().zip(&x) {
            out.push(match (slot, p) {
                (Some(id), _) => *id,
                (None, Some(p)) => self.push(*p)?,
                (None, None) => return Err(Error::Failed("a block node is missing".into())),
            });
        }
        Ok(out)
    }
}

/// A face's grid of node ids in its owner's frame: `(na + 1) × (nb + 1)`, node `(a, b)` at
/// `b·(na + 1) + a`, corners at the owner's ring.
struct FaceGrid {
    ring: [usize; 4],
    na: usize,
    nb: usize,
    ids: Vec<u32>,
}

impl FaceGrid {
    /// The lookup of node `(a, b)` of a cell's view of this face, whose ring (the order of
    /// [`face_ring`]) is `ring` and counts `na × nb`.
    fn map_from(&self, ring: [usize; 4], na: usize, nb: usize) -> Result<impl Fn(usize, usize) -> Result<u32> + '_> {
        // The owner's lattice position of each of the cell's corners.
        let corner = |v: usize| -> Option<[i64; 2]> {
            let k = self.ring.iter().position(|&x| x == v)?;
            let (a, b) = [(0, 0), (self.na, 0), (self.na, self.nb), (0, self.nb)].get(k).copied()?;
            Some([a as i64, b as i64])
        };
        let (Some(o), Some(pa), Some(pb)) = (corner(ring[0]), corner(ring[1]), corner(ring[3])) else {
            return Err(Error::Failed("a cell's face does not have the shared face's corners".into()));
        };
        let step = |p: [i64; 2], n: usize| -> Option<[i64; 2]> {
            let d = [p[0] - o[0], p[1] - o[1]];
            let n = n as i64;
            // One axis of the owner's lattice, `n` long.
            (n > 0 && d.iter().filter(|&&x| x != 0).count() == 1 && d.iter().all(|&x| x == 0 || x.abs() == n)).then(|| d.map(|x| x / n))
        };
        let (Some(du), Some(dv)) = (step(pa, na), step(pb, nb)) else {
            return Err(Error::Failed("a cell's view of a shared face does not match its grid".into()));
        };
        let (na_owner, ids) = (self.na, &self.ids);
        Ok(move |a: usize, b: usize| -> Result<u32> {
            let (a, b) = (a as i64, b as i64);
            let u = o[0] + a * du[0] + b * dv[0];
            let v = o[1] + a * du[1] + b * dv[1];
            usize::try_from(v)
                .ok()
                .zip(usize::try_from(u).ok())
                .and_then(|(v, u)| ids.get(v * (na_owner + 1) + u))
                .copied()
                .ok_or_else(|| Error::Failed("shared face grid index".into()))
        })
    }
}
