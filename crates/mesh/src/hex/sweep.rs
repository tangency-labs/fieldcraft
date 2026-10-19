//! One-to-one sweeping from a planar source face to a planar target face through rings of
//! linking faces.
//!
//! **Recognition** (Blacker, *The Cooper tool*, Proc. 5th International Meshing Roundtable,
//! 1996). Starting from the source's loops, each loop edge leads to the linking face across it.
//! A linking face has one loop made of a *bottom* chain (the edges it shares with the previous
//! ring), one *side* edge, a *top* chain and another side edge; its top chains form the next ring.
//! Rings are followed until every top edge of a loop borders one planar face, the target. Every
//! face but the source and the target must be used exactly once.
//!
//! **Intervals** ([`super::interval`]): on each linking face the bottom and top chains have equal
//! sums and the two sides are equal; the source is mapped when it is a four-cornered single loop,
//! otherwise its loops are even for the free all-quad mesher. Every ring has the same number of
//! layers per stage.
//!
//! **Placement.** Edge nodes lie on both faces of their edge; each linking face gets the
//! transfinite grid of its four sides (Gordon and Hall 1973) projected onto its surface. The
//! interior nodes of each layer follow from the source mesh by the least-squares affine map that
//! carries the source boundary ring onto the layer's ring (Roca and Sarrate, *An automatic and
//! general least-squares projection procedure for sweep meshing*, Engineering with Computers
//! 26(4):391–406, 2010), plus the map's residuals at the ring interpolated by inverse-distance
//! weights (the boundary-residual correction of Knupp, *Next-generation sweep tool: a method for
//! generating all-hex meshes on two-and-one-half dimensional geometries*, Proc. 7th International
//! Meshing Roundtable, 1998). The last layer is projected onto the target face.

use std::collections::{BTreeMap, BTreeSet};

use super::brep::Brep;
use super::geom::{add, scale};
use super::interval::Intervals;
use super::layered::{Frame, PlaneField};
use super::region;
use super::{Ctx, Raw, project_to_face, project_to_faces};
use crate::quad::transfinite_grid;
use crate::size::ideal_divisions;
use crate::vec3::{P, dot, sub};
use crate::{Error, Result};

/// The most source faces tried.
const MAX_SOURCES: usize = 8;
/// The most interior-node × ring-node × layer products the residual correction may spend.
const MAX_CORRECTION_WORK: u64 = 2_000_000_000;

/// A linking face of one stage: its bottom chain, sides and top chain (oriented edges as the face
/// traverses them).
#[derive(Debug, Clone)]
struct Link {
    face: usize,
    bottom: Vec<(usize, bool)>,
    right: (usize, bool),
    top: Vec<(usize, bool)>,
    left: (usize, bool),
}

/// The sweep's topology: per source loop, its stages, each a ring of linking faces in loop order.
struct Topology {
    source: usize,
    target: usize,
    rings: Vec<Vec<Vec<Link>>>,
}

/// Sweep the body, with the source's curves `2^refine` times finer than the size asks.
pub(crate) fn mesh(ctx: &Ctx, refine: u32) -> Result<Raw> {
    let brep = ctx.brep;
    let mut planar: Vec<usize> = (0..brep.faces.len()).filter(|&f| brep.faces.get(f).is_some_and(|x| x.plane.is_some())).collect();
    if planar.is_empty() {
        return Err(Error::Unsupported("no planar face to sweep from".into()));
    }
    let area = |f: usize| brep.faces.get(f).map_or(0.0, |x| x.area);
    planar.sort_by(|&a, &b| area(b).total_cmp(&area(a)).then(a.cmp(&b)));
    let mut reasons = Vec::new();
    for &s in planar.iter().take(MAX_SOURCES) {
        match topology(brep, s) {
            Ok(topo) => return build(ctx, &topo, refine),
            Err(e) => reasons.push(format!("from face {}: {e}", brep.face_tag(s))),
        }
    }
    Err(Error::Unsupported(reasons.into_iter().take(2).collect::<Vec<_>>().join("; ")))
}

/// Follow the rings of linking faces from source face `s`.
fn topology(brep: &Brep, s: usize) -> std::result::Result<Topology, String> {
    let source = brep.faces.get(s).ok_or("no such face")?;
    let mut used: BTreeSet<usize> = BTreeSet::from([s]);
    let mut rings = Vec::new();
    let mut target: Option<usize> = None;
    for front0 in &source.loops {
        // Front edges as the face below them traverses them, with that face.
        let mut front: Vec<(usize, bool, usize)> = front0.iter().map(|&(e, f)| (e, f, s)).collect();
        let mut stages: Vec<Vec<Link>> = Vec::new();
        loop {
            if stages.len() > brep.faces.len() {
                return Err("the rings do not end".into());
            }
            let across: Vec<usize> = front.iter().map(|&(e, _, below)| brep.edges.get(e).map_or(usize::MAX, |x| x.other(below))).collect();
            let first = across.first().copied().ok_or("an empty loop")?;
            if across.iter().all(|&f| f == first) && brep.faces.get(first).is_some_and(|x| x.plane.is_some()) && is_loop_of(brep, first, &front) {
                match target {
                    None => target = Some(first),
                    Some(t) if t == first => {}
                    Some(_) => return Err("the loops end on different faces (not one-to-one)".into()),
                }
                break;
            }
            // Rotate the front to start where the face across changes, then take runs.
            let n = front.len();
            let start = (0..n).find(|&k| across.get((k + n - 1) % n) != across.get(k)).unwrap_or(0);
            let mut ring: Vec<Link> = Vec::new();
            let mut k = 0;
            while k < n {
                let f = across.get((start + k) % n).copied().ok_or("front index")?;
                if used.contains(&f) || ring.iter().any(|l| l.face == f) {
                    return Err(format!("face {} is reached twice (not a one-to-one sweep)", brep.face_tag(f)));
                }
                let mut run = Vec::new();
                while k < n && across.get((start + k) % n) == Some(&f) {
                    let &(e, fwd, _) = front.get((start + k) % n).ok_or("front index")?;
                    run.push((e, fwd));
                    k += 1;
                }
                ring.push(link(brep, f, &run).map_err(|e| match cap(brep, s, f) {
                    Some(why) => why,
                    None => e,
                })?);
            }
            used.extend(ring.iter().map(|l| l.face));
            front = ring.iter().flat_map(|l| l.top.iter().map(move |&(e, fwd)| (e, fwd, l.face))).collect();
            stages.push(ring);
        }
        if stages.is_empty() {
            return Err("the source touches the target".into());
        }
        rings.push(stages);
    }
    let target = target.ok_or("no target face")?;
    used.insert(target);
    if used.len() != brep.faces.len() {
        return Err(format!("{} faces are not on the sweep path", brep.faces.len() - used.len()));
    }
    if brep.faces.get(target).map(|f| f.loops.len()) != Some(source.loops.len()) {
        return Err("the source and target have different numbers of loops".into());
    }
    Ok(Topology { source: s, target, rings })
}

/// Why face `f`, reached from a front swept from source `s`, cannot be a linking face when it is
/// a plane parallel to the source: a cap that only part of the ring reaches. The body then ends
/// at several levels (a step, a boss, a pocket floor), which no one-to-one sweep takes; saying so
/// names the reason where the linking-face split would only report the symptom (#444).
fn cap(brep: &Brep, s: usize, f: usize) -> Option<String> {
    let (_, ns) = brep.faces.get(s)?.plane?;
    let (_, nf) = brep.faces.get(f)?.plane?;
    (dot(ns, nf).abs() > 1.0 - 1e-9).then(|| {
        format!(
            "face {} is parallel to the source and only part of the ring reaches it: the body ends at more than one level (a step, a boss or a pocket), and a one-to-one sweep has one target face; method `extrude` meshes such bodies level by level",
            brep.face_tag(f)
        )
    })
}

/// Whether the front's edges are exactly one loop of face `t`.
fn is_loop_of(brep: &Brep, t: usize, front: &[(usize, bool, usize)]) -> bool {
    let edges: BTreeSet<usize> = front.iter().map(|x| x.0).collect();
    brep.faces.get(t).is_some_and(|f| f.loops.iter().any(|l| l.iter().map(|x| x.0).collect::<BTreeSet<_>>() == edges))
}

/// Split linking face `f`'s loop into bottom (the front `run`), sides and top.
fn link(brep: &Brep, f: usize, run: &[(usize, bool)]) -> std::result::Result<Link, String> {
    let face = brep.faces.get(f).ok_or("face index")?;
    let tag = face.id;
    let [lp] = face.loops.as_slice() else {
        return Err(format!("linking face {tag} has {} boundary loops (periodic or holed linking faces are not swept)", face.loops.len()));
    };
    let n = lp.len();
    let bottom_edges: BTreeSet<usize> = run.iter().map(|x| x.0).collect();
    let in_bottom = |k: usize| lp.get(k % n).is_some_and(|x| bottom_edges.contains(&x.0));
    let start = (0..n)
        .find(|&k| in_bottom(k) && !in_bottom(k + n - 1))
        .ok_or_else(|| format!("linking face {tag}: the shared edges do not form one chain"))?;
    let len = (0..n).take_while(|&k| in_bottom(start + k)).count();
    if len != run.len() || (0..n).filter(|&k| in_bottom(k)).count() != len {
        return Err(format!("linking face {tag}: the shared edges do not form one chain"));
    }
    if n < len + 3 {
        return Err(format!("linking face {tag} has no top side"));
    }
    let at = |k: usize| lp.get(k % n).copied().ok_or_else(|| "loop index".to_string());
    let bottom: Vec<(usize, bool)> = (0..len).map(|k| at(start + k)).collect::<std::result::Result<_, _>>()?;
    let right = at(start + len)?;
    let left = at(start + n - 1)?;
    let top: Vec<(usize, bool)> = (len + 1..n - 1).map(|k| at(start + k)).collect::<std::result::Result<_, _>>()?;
    Ok(Link { face: f, bottom, right, top, left })
}

/// Node positions while building.
struct Nodes {
    xyz: Vec<P>,
}

impl Nodes {
    fn add(&mut self, p: P) -> u32 {
        self.xyz.push(p);
        self.xyz.len() as u32 - 1
    }
    fn get(&self, i: u32) -> P {
        self.xyz.get(i as usize).copied().unwrap_or([0.0; 3])
    }
}

/// Node ids along oriented edge uses, end to end (each use's last node is the next one's first).
fn chain(edge_nodes: &[Vec<u32>], uses: &[(usize, bool)]) -> Vec<u32> {
    let mut out: Vec<u32> = Vec::new();
    for &(e, fwd) in uses {
        let mut ids = edge_nodes.get(e).cloned().unwrap_or_default();
        if !fwd {
            ids.reverse();
        }
        if !out.is_empty() && !ids.is_empty() {
            ids.remove(0);
        }
        out.extend(ids);
    }
    out
}

fn build(ctx: &Ctx, topo: &Topology, refine: u32) -> Result<Raw> {
    let brep = ctx.brep;
    let source = brep.faces.get(topo.source).ok_or_else(|| Error::Failed("source face".into()))?;
    let (s_point, s_normal) = source.plane.ok_or_else(|| Error::Failed("the source is not planar".into()))?;
    // Interval assignment over the B-rep edges.
    let mut ia = Intervals::default();
    let sides: BTreeSet<usize> = topo.rings.iter().flatten().flatten().flat_map(|l| [l.right.0, l.left.0]).collect();
    let fine = f64::from(1u32 << refine.min(8));
    for e in 0..brep.edges.len() {
        let goal = ideal_divisions(ctx.field, &brep.edge_points(e, true))?;
        ctx.edge_var(&mut ia, e, if sides.contains(&e) { goal } else { goal * fine });
    }
    let vars = |uses: &[(usize, bool)]| uses.iter().map(|x| x.0).collect::<Vec<_>>();
    for stages in &topo.rings {
        for ring in stages {
            for l in ring {
                ia.equal.push((vars(&l.bottom), vars(&l.top)));
                ia.equal.push((vec![l.right.0], vec![l.left.0]));
            }
        }
    }
    // The same layers on every ring: stage by stage when the rings have as many stages.
    let side_of = |stages: &Vec<Vec<Link>>, k: usize| stages.get(k).and_then(|r| r.first()).map(|l| l.right.0);
    if let Some(first) = topo.rings.first() {
        for other in topo.rings.iter().skip(1) {
            if other.len() == first.len() {
                for k in 0..first.len() {
                    if let (Some(a), Some(b)) = (side_of(first, k), side_of(other, k)) {
                        ia.equal.push((vec![a], vec![b]));
                    }
                }
            } else {
                let all = |st: &Vec<Vec<Link>>| (0..st.len()).filter_map(|k| side_of(st, k)).collect::<Vec<_>>();
                ia.equal.push((all(first), all(other)));
            }
        }
    }
    // The source's frame: u × v along the sweep (into the solid).
    let frame = Frame::translate(s_point, scale(s_normal, -1.0)).ok_or_else(|| Error::Failed("source normal".into()))?;
    let q = |p: P| frame.coords(p).0;
    let polylines: Vec<Vec<Vec<[f64; 2]>>> =
        source.loops.iter().map(|lp| lp.iter().map(|&(e, f)| brep.edge_points(e, f).into_iter().map(q).collect()).collect()).collect();
    let scheme = region::classify(&polylines);
    let loop_vars: Vec<Vec<usize>> = source.loops.iter().map(|lp| vars(lp)).collect();
    region::constrain(&mut ia, &loop_vars, scheme);
    let counts = ia.solve()?;
    // A biased seed along the sweep carries its layers to every side edge of its stage: a
    // linking face's right side runs up as its left side runs down, and the first side of each
    // stage runs up on every ring.
    let mut ties: Vec<[(usize, bool); 2]> = Vec::new();
    for stages in &topo.rings {
        for ring in stages {
            for l in ring {
                ties.push([l.right, (l.left.0, !l.left.1)]);
            }
        }
    }
    if let Some(first) = topo.rings.first() {
        for other in topo.rings.iter().skip(1).filter(|o| o.len() == first.len()) {
            for (a, b) in first.iter().zip(other) {
                if let (Some(a), Some(b)) = (a.first(), b.first()) {
                    ties.push([a.right, b.right]);
                }
            }
        }
    }
    let carried = ctx.carried(&ties);
    // Edge nodes: B-rep vertices once, interiors on both faces of the edge.
    let mut nodes = Nodes { xyz: Vec::new() };
    let mut vertex_node: BTreeMap<u32, u32> = BTreeMap::new();
    let mut edge_nodes: Vec<Vec<u32>> = Vec::with_capacity(brep.edges.len());
    for (e, edge) in brep.edges.iter().enumerate() {
        let n = counts.get(e).copied().unwrap_or(1);
        let pts = ctx.divide_edge(e, true, n, &carried)?;
        let last = pts.len().saturating_sub(1);
        let mut ids = Vec::with_capacity(pts.len());
        for (k, p) in pts.iter().enumerate() {
            let end = if k == 0 {
                Some(edge.first())
            } else if k == last {
                Some(edge.last())
            } else {
                None
            };
            ids.push(match end {
                Some(v) => *vertex_node.entry(v).or_insert_with(|| nodes.add(brep.point(v))),
                None => nodes.add(project_to_faces(ctx, &edge.faces, *p)),
            });
        }
        edge_nodes.push(ids);
    }
    // Linking-face grids and the ring of every layer: rings_by_layer[loop][layer].
    let mut rings_by_layer: Vec<Vec<Vec<u32>>> = Vec::with_capacity(topo.rings.len());
    for stages in &topo.rings {
        let mut layers: Vec<Vec<u32>> = Vec::new();
        for (k, ring) in stages.iter().enumerate() {
            let n_side = ring.first().map_or(0, |l| chain(&edge_nodes, &[l.right]).len().saturating_sub(1));
            let mut rows: Vec<Vec<u32>> = vec![Vec::new(); n_side + 1];
            for l in ring {
                let grid = face_grid(ctx, &mut nodes, &edge_nodes, l)?;
                if grid.len() != n_side + 1 {
                    return Err(Error::Failed(format!("linking face {} has sides of different lengths", brep.face_tag(l.face))));
                }
                // Rows run along the face's bottom traversal, which is the front's direction.
                for (row, slot) in grid.iter().zip(rows.iter_mut()) {
                    let mut r = row.clone();
                    r.pop();
                    slot.extend(r);
                }
            }
            // A later stage's rings may start elsewhere around the loop: turn them to continue the
            // previous stage's columns.
            if let (Some(prev), Some(first)) = (layers.last(), rows.first()) {
                let n = first.len();
                let turn = (0..n)
                    .find(|&o| (0..n).all(|i| first.get((i + o) % n) == prev.get(i)))
                    .ok_or_else(|| Error::Failed(format!("stage {k} of the sweep does not continue the previous stage's ring")))?;
                for r in rows.iter_mut() {
                    r.rotate_left(turn);
                }
            }
            layers.extend(rows.into_iter().enumerate().filter(|(j, _)| k == 0 || *j > 0).map(|(_, r)| r));
        }
        rings_by_layer.push(layers);
    }
    let n_layers = rings_by_layer.first().map_or(0, Vec::len).saturating_sub(1);
    if n_layers == 0 || rings_by_layer.iter().any(|l| l.len() != n_layers + 1) {
        return Err(Error::Failed("the rings have different numbers of layers".into()));
    }
    // Source mesh on the layer-0 rings, which run as the source's loops do.
    let ring0: Vec<Vec<u32>> = rings_by_layer.iter().map(|l| l.first().cloned().unwrap_or_default()).collect();
    let loops2: Vec<Vec<[f64; 2]>> = ring0.iter().map(|r| r.iter().map(|&i| q(nodes.get(i))).collect()).collect();
    let plane_field = PlaneField { field: ctx.field, frame, h: 0.0 };
    // Where each loop's edge uses start in its ring: at the edges' first vertices.
    let starts: Vec<Vec<usize>> = source
        .loops
        .iter()
        .zip(&ring0)
        .map(|(lp, ring)| {
            lp.iter()
                .filter_map(|&(e, fwd)| {
                    brep.edge_vertices(e, fwd).first().and_then(|v| vertex_node.get(v)).and_then(|id| ring.iter().position(|x| x == id))
                })
                .collect()
        })
        .collect();
    // The rings may start mid-loop: give the mapped scheme its corners as ring positions.
    // Corners (curve uses) as ring positions, in order; `None` unless every one is found.
    let positions = |c: &[usize]| {
        let s0 = starts.first().cloned().unwrap_or_default();
        let mut pos: Vec<usize> = c.iter().filter_map(|&k| s0.get(k).copied()).collect();
        pos.sort_unstable();
        (pos.len() == c.len()).then_some(pos)
    };
    let (scheme, starts) = match scheme {
        region::Scheme::Mapped(c) => match positions(&c) {
            Some(pos) => (region::Scheme::Mapped([0, 1, 2, 3]), vec![pos]),
            None => (region::Scheme::Free, starts),
        },
        region::Scheme::Triangle(c) => match positions(&c) {
            Some(pos) => (region::Scheme::Triangle([0, 1, 2]), vec![pos]),
            None => (region::Scheme::Free, starts),
        },
        other => (other, starts),
    };
    let (mesh, used_scheme) = region::mesh(&loops2, &starts, scheme, &plane_field, &ctx.opts.quad)?;
    let boundary_ids: Vec<u32> = ring0.concat();
    let n_boundary = boundary_ids.len();
    if mesh.boundary_nodes != n_boundary {
        return Err(Error::Failed("the source mesh changed its boundary".into()));
    }
    let interior2: Vec<[f64; 2]> = mesh.nodes.iter().skip(n_boundary).copied().collect();
    let count = mesh.quads.len().saturating_mul(n_layers);
    if count > ctx.opts.max_elements {
        return Err(Error::Limit(format!("the mesh would have {count} hexahedra; the limit is {}", ctx.opts.max_elements)));
    }
    let work = (interior2.len() as u64).saturating_mul(n_boundary as u64).saturating_mul(n_layers as u64);
    if work > MAX_CORRECTION_WORK {
        return Err(Error::Limit(format!("placing the sweep's interior nodes would take {work} steps; the limit is {MAX_CORRECTION_WORK}")));
    }
    let s_ring: Vec<[f64; 2]> = boundary_ids.iter().map(|&i| q(nodes.get(i))).collect();
    // Interior nodes, layer by layer. The residual weights depend only on the source plane.
    let weights: Vec<Vec<(usize, f64)>> = interior2.iter().map(|p| shepard_weights(*p, &s_ring)).collect();
    let mut layer_nodes: Vec<Vec<u32>> = Vec::with_capacity(n_layers + 1);
    for j in 0..=n_layers {
        let ring: Vec<u32> = rings_by_layer.iter().flat_map(|l| l.get(j).cloned().unwrap_or_default()).collect();
        if ring.len() != n_boundary {
            return Err(Error::Failed(format!("layer {j}'s ring has {} nodes, the source's {n_boundary}", ring.len())));
        }
        let mut ids = ring.clone();
        if j == 0 {
            ids.extend(interior2.iter().map(|p| nodes.add(frame.point(*p, 0.0))));
        } else {
            let x: Vec<P> = ring.iter().map(|&i| nodes.get(i)).collect();
            let map = Affine::fit(&s_ring, &x)?;
            let residual: Vec<P> = s_ring.iter().zip(&x).map(|(s, xi)| sub(*xi, map.at(*s))).collect();
            for (p, w) in interior2.iter().zip(&weights) {
                let correction = w.iter().fold([0.0; 3], |acc, &(k, wk)| add(acc, scale(residual.get(k).copied().unwrap_or([0.0; 3]), wk)));
                let mut y = add(map.at(*p), correction);
                if j == n_layers {
                    y = project_to_face(ctx, topo.target, y).unwrap_or(y);
                }
                ids.push(nodes.add(y));
            }
        }
        layer_nodes.push(ids);
    }
    let mut hexes = Vec::with_capacity(count);
    for quad in &mesh.quads {
        for w in layer_nodes.windows(2) {
            let (lo, hi) = (&w[0], &w[1]);
            let b = quad.map(|v| lo.get(v as usize).copied().unwrap_or(0));
            let t = quad.map(|v| hi.get(v as usize).copied().unwrap_or(0));
            hexes.push([b[0], b[1], b[2], b[3], t[0], t[1], t[2], t[3]]);
        }
    }
    let stages: usize = topo.rings.first().map_or(0, Vec::len);
    let description = crate::message![format!(
        "swept from face {} to face {} through {stages} stage{} of linking faces: {n_layers} layers of {} {} source quadrilaterals",
        brep.face_tag(topo.source),
        brep.face_tag(topo.target),
        if stages == 1 { "" } else { "s" },
        mesh.quads.len(),
        match used_scheme {
            region::Scheme::Mapped(_) => "mapped",
            region::Scheme::Triangle(_) => "triangle",
            region::Scheme::Annulus => "annulus",
            region::Scheme::Disk => "disk",
            region::Scheme::Free => "free",
        }
    )];
    Ok(Raw { nodes: nodes.xyz, hexes, description })
}

/// The transfinite grid of a linking face projected onto it: `grid[j][i]`, `i` along the bottom
/// as the front runs, `j` from bottom to top.
fn face_grid(ctx: &Ctx, nodes: &mut Nodes, edge_nodes: &[Vec<u32>], l: &Link) -> Result<Vec<Vec<u32>>> {
    // The face traverses bottom → right → top → left; the front runs against its bottom.
    let mut bottom = chain(edge_nodes, &l.bottom);
    bottom.reverse();
    let mut right = chain(edge_nodes, &[l.left]);
    right.reverse();
    let top = chain(edge_nodes, &l.top);
    let left = chain(edge_nodes, &[l.right]);
    let pts = |ids: &[u32]| ids.iter().map(|&i| nodes.get(i)).collect::<Vec<P>>();
    let grid = transfinite_grid(&pts(&bottom), &pts(&right), &pts(&top), &pts(&left))
        .map_err(|e| Error::Failed(format!("linking face {}: {e}", ctx.brep.face_tag(l.face))))?;
    let (nu, nv) = (bottom.len().saturating_sub(1), right.len().saturating_sub(1));
    let mut out = Vec::with_capacity(nv + 1);
    for j in 0..=nv {
        let mut row = Vec::with_capacity(nu + 1);
        for i in 0..=nu {
            let id = if j == 0 {
                bottom.get(i).copied()
            } else if j == nv {
                top.get(i).copied()
            } else if i == 0 {
                left.get(j).copied()
            } else if i == nu {
                right.get(j).copied()
            } else {
                let p = grid.get(j * (nu + 1) + i).copied().unwrap_or([0.0; 3]);
                Some(nodes.add(project_to_face(ctx, l.face, p).unwrap_or(p)))
            };
            row.push(id.ok_or_else(|| Error::Failed("grid side".into()))?);
        }
        out.push(row);
    }
    Ok(out)
}

/// `x ≈ c + M s`: the least-squares affine map from the source plane into space.
struct Affine {
    c: P,
    m: [P; 2],
}

impl Affine {
    fn fit(s: &[[f64; 2]], x: &[P]) -> Result<Affine> {
        if s.len() != x.len() || s.len() < 3 {
            return Err(Error::Failed("a sweep layer needs at least three ring nodes".into()));
        }
        // Normal equations of [1, s_u, s_v], centred for conditioning.
        let n = s.len() as f64;
        let (mu, mv) = (s.iter().map(|p| p[0]).sum::<f64>() / n, s.iter().map(|p| p[1]).sum::<f64>() / n);
        let (mut suu, mut suv, mut svv) = (0.0, 0.0, 0.0);
        let (mut xu, mut xv, mut xm) = ([0.0; 3], [0.0; 3], [0.0; 3]);
        for (p, xi) in s.iter().zip(x) {
            let (u, v) = (p[0] - mu, p[1] - mv);
            suu += u * u;
            suv += u * v;
            svv += v * v;
            xu = add(xu, scale(*xi, u));
            xv = add(xv, scale(*xi, v));
            xm = add(xm, scale(*xi, 1.0 / n));
        }
        let det = suu * svv - suv * suv;
        if !(det > 1e-12 * suu * svv) {
            return Err(Error::Failed("the source ring is collinear; the sweep's layers cannot be mapped".into()));
        }
        let mu_vec = scale(sub(scale(xu, svv), scale(xv, suv)), 1.0 / det);
        let mv_vec = scale(sub(scale(xv, suu), scale(xu, suv)), 1.0 / det);
        let c = sub(xm, add(scale(mu_vec, mu), scale(mv_vec, mv)));
        Ok(Affine { c, m: [mu_vec, mv_vec] })
    }

    fn at(&self, s: [f64; 2]) -> P {
        add(self.c, add(scale(self.m[0], s[0]), scale(self.m[1], s[1])))
    }
}

/// Inverse-distance (Shepard, power 4) weights of the ring nodes at `p`, normalized: the
/// residual correction at `p` is their weighted sum of the ring's residuals.
pub(super) fn shepard_weights(p: [f64; 2], ring: &[[f64; 2]]) -> Vec<(usize, f64)> {
    let mut w: Vec<(usize, f64)> = Vec::with_capacity(ring.len());
    for (k, s) in ring.iter().enumerate() {
        let d2 = (p[0] - s[0]).powi(2) + (p[1] - s[1]).powi(2);
        if d2 == 0.0 {
            return vec![(k, 1.0)];
        }
        w.push((k, 1.0 / (d2 * d2)));
    }
    let sum: f64 = w.iter().map(|x| x.1).sum();
    if sum > 0.0 {
        for x in &mut w {
            x.1 /= sum;
        }
    }
    w
}
