//! The free all-quadrilateral mesher: half-resolution constrained Delaunay triangulation,
//! recombination into quadrilaterals, one uniform split that lands on the prescribed boundary
//! nodes, then constrained smoothing. See the module docs of [`super`] for the method and its
//! sources.

use std::collections::BTreeMap;

use super::cdt::{CLEARANCE, Cdt, Triangulation};
use super::domain::{self, Region};
use super::guard::Guard;
use super::{QuadMesh, QuadOptions, finish};
use crate::predicates::{Sign, orient2d};
use crate::quality::{planar_corner, turns_left};
use crate::size::SizeField;
use crate::{Error, Result};

/// Mesh the region bounded by `loops` (the first is the outer boundary, the rest are holes;
/// either orientation) with quadrilaterals only. See [`super::free_quad_mesh`].
pub(crate) fn mesh(loops: &[Vec<[f64; 2]>], size: &dyn SizeField, opts: &QuadOptions) -> Result<QuadMesh> {
    let region = domain::region(loops)?;
    for (k, l) in region.loops.iter().enumerate() {
        if l.len() % 2 == 1 {
            return Err(Error::Invalid(format!(
                "loop {k} has an odd number ({}) of boundary points; an all-quadrilateral mesh needs an even number on every loop \
                 (choose the edge divisions with fieldcraft_mesh::size::loop_divisions or even_divisions)",
                l.len()
            )));
        }
    }
    if region.loops.len() == 1 && region.loops.first().is_some_and(|l| l.len() == 4) {
        let l = region.loops.first().cloned().unwrap_or_default();
        let quad = [l.first(), l.get(1), l.get(2), l.get(3)];
        let quad = match quad {
            [Some(&a), Some(&b), Some(&c), Some(&d)] => [a, b, c, d],
            _ => return Err(Error::Failed("quad loop".into())),
        };
        return finish(region.pts.clone(), vec![quad], region.pts.len(), opts, Default::default());
    }
    if let Some((k, l)) = region.loops.iter().enumerate().find(|(_, l)| l.len() < 6) {
        return Err(Error::Invalid(format!(
            "loop {k} has {} boundary points; the free quad mesher needs at least 6 on every loop of a region with more than one loop \
             (split its sides)",
            l.len()
        )));
    }

    // Both parity classes, preferred first; the second is tried when the first cannot be made
    // valid.
    let mut last = None;
    for coarse in coarse_regions(&region)? {
        match attempt(&region, coarse, size, opts) {
            Ok(mesh) => return Ok(mesh),
            Err(e) if e.kind() == Some(crate::Kind::Failed) => last = Some(e),
            Err(e) => return Err(e),
        }
    }
    Err(last.unwrap_or_else(|| Error::Failed("no half-resolution boundary".into())))
}

/// One run of the pipeline on one half-resolution region.
fn attempt(region: &Region, coarse: Coarse, size: &dyn SizeField, opts: &QuadOptions) -> Result<QuadMesh> {
    let Coarse { region: coarse, to_input, skipped, before } = coarse;
    let chords = chords(&coarse, &skipped)?;
    // For each coarse boundary edge (coarse ids, smaller first), the input point it skips.
    let boundary_mid: BTreeMap<(u32, u32), u32> = chords.iter().map(|&(a, b, m)| ((a.min(b), a.max(b)), m)).collect();
    let guard =
        Guard::new(&chords.iter().map(|&(a, b, m)| Ok((*coarse.point(a)?, *coarse.point(b)?, *region.point(m)?))).collect::<Result<Vec<_>>>()?)?;
    let need = fan_needs(region, &to_input, &before, &skipped)?;
    let max_tris = opts.max_elements / 2 + 4 * region.pts.len() + 16;
    let mut cdt = Cdt::build(&coarse, max_tris)?;
    cdt.refine(size, 2.0, &|p, h| guard.allows(p, CLEARANCE * h))?;
    cdt.ensure_fans(&need)?;
    let tri = cdt.into_triangulation();
    let (quads, leftovers) = recombine(&tri, &need)?;
    let count = 4 * quads.len() + 3 * leftovers.len();
    if count > opts.max_elements {
        return Err(Error::Limit(format!("the mesh would have {count} quadrilaterals; the limit is {}", opts.max_elements)));
    }
    let (nodes, elements) = split(region, &tri, &to_input, &boundary_mid, &quads, &leftovers)?;
    let stats = super::QuadStats { recombined: quads.len(), leftover_triangles: leftovers.len(), ..Default::default() };
    finish(nodes, elements, region.pts.len(), opts, stats)
}

/// `1 − cos θ` of the turn at each loop position (0 on a straight run, 2 on a full reversal).
fn turns(region: &Region, l: &[u32]) -> Result<Vec<f64>> {
    let n = l.len();
    let at = |i: usize| -> Result<[f64; 2]> {
        let v = l.get(i % n).copied().ok_or_else(|| Error::Failed("loop index".into()))?;
        region.pts.get(v as usize).copied().ok_or_else(|| Error::Failed("point index".into()))
    };
    (0..n)
        .map(|i| {
            let (a, b, c) = (at(i + n - 1)?, at(i)?, at(i + 1)?);
            let (u, w) = ([b[0] - a[0], b[1] - a[1]], [c[0] - b[0], c[1] - b[1]]);
            let den = (u[0] * u[0] + u[1] * u[1]).sqrt() * (w[0] * w[0] + w[1] * w[1]).sqrt();
            let cos = (u[0] * w[0] + u[1] * w[1]) / den;
            Ok(if den > 0.0 && cos.is_finite() { 1.0 - cos } else { 0.0 })
        })
        .collect()
}

/// The half-resolution regions to try: every other point of each loop, first keeping the parity
/// class with the sharper turns (so loop corners survive when they sit at an even spacing), then
/// the other class of every loop. Choices whose chords cross or flip are dropped.
fn coarse_regions(region: &Region) -> Result<Vec<Coarse>> {
    let mut preferred = Vec::with_capacity(region.loops.len());
    for l in &region.loops {
        let t = turns(region, l)?;
        let even: f64 = t.iter().step_by(2).sum();
        let odd: f64 = t.iter().skip(1).step_by(2).sum();
        preferred.push(usize::from(odd > even));
    }
    let mut out = Vec::new();
    let mut last = None;
    for flip in [0, 1] {
        let offsets: Vec<usize> = preferred.iter().map(|&o| o ^ flip).collect();
        let coarse = pick(region, &offsets)?;
        match check_coarse(&coarse.region) {
            Ok(()) => out.push(coarse),
            Err(e) => last = Some(e),
        }
    }
    if out.is_empty() {
        return Err(Error::Failed(format!(
            "the boundary is too coarse for the free quad mesher: the chords joining every other boundary point cross or flip ({}); \
             refine the boundary divisions",
            last.map(|e| e.to_string()).unwrap_or_default()
        )));
    }
    Ok(out)
}

/// The half-resolution region, with each coarse point's input id and the input id of the
/// boundary point after it (the one the coarse edge skips).
struct Coarse {
    region: Region,
    to_input: Vec<u32>,
    skipped: Vec<u32>,
    /// The input point before each coarse point.
    before: Vec<u32>,
}

fn pick(region: &Region, offsets: &[usize]) -> Result<Coarse> {
    let (mut pts, mut p2, mut loops, mut to_input, mut skipped, mut before) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for (l, &off) in region.loops.iter().zip(offsets) {
        let mut ids = Vec::with_capacity(l.len() / 2);
        for pos in (off..l.len()).step_by(2) {
            let (v, after, prior) = match (l.get(pos), l.get((pos + 1) % l.len()), l.get((pos + l.len() - 1) % l.len())) {
                (Some(&v), Some(&a), Some(&b)) => (v, a, b),
                _ => return Err(Error::Failed("loop index".into())),
            };
            ids.push(pts.len() as u32);
            pts.push(region.pts.get(v as usize).copied().ok_or_else(|| Error::Failed("point index".into()))?);
            p2.push(*region.point(v)?);
            to_input.push(v);
            skipped.push(after);
            before.push(prior);
        }
        loops.push(ids);
    }
    Ok(Coarse { region: Region { pts, p2, loops }, to_input, skipped, before })
}

/// The coarse loops must be simple, keep their orientation and nesting.
fn check_coarse(coarse: &Region) -> Result<()> {
    domain::check_simple(coarse)?;
    domain::check_orientation(coarse)?;
    domain::check_nesting(coarse)
}

/// The coarse boundary edges `(a, b, m)`: coarse ids in loop order and the input point `m` the
/// edge skips.
fn chords(coarse: &Region, skipped: &[u32]) -> Result<Vec<(u32, u32, u32)>> {
    let mut out = Vec::with_capacity(coarse.pts.len());
    for (a, b) in domain::segments(&coarse.loops) {
        out.push((a, b, skipped.get(a as usize).copied().ok_or_else(|| Error::Failed("coarse id".into()))?));
    }
    Ok(out)
}

/// How many coarse elements each coarse boundary vertex needs around it. A coarse vertex `v`
/// gets one final quadrilateral per incident coarse element, and the corner of the first and last
/// of them at `v` opens towards the prescribed neighbours of `v`. A single element would span the
/// whole boundary angle at `v`, so a straight or nearly straight (interior angle above 150°) or
/// reflex vertex needs two, and one above 300° needs three. Decided with the exact `orient2d` and
/// the cosine of the turn.
fn fan_needs(region: &Region, to_input: &[u32], before: &[u32], after: &[u32]) -> Result<Vec<usize>> {
    let mut need = Vec::with_capacity(to_input.len());
    for ((&v, &a), &b) in to_input.iter().zip(before).zip(after) {
        let (pa, pv, pb) = (region.point(a)?, region.point(v)?, region.point(b)?);
        let o = orient2d(pa, pv, pb);
        let (a, v, b) = (pa.coords(), pv.coords(), pb.coords());
        let (u, w) = ([v[0] - a[0], v[1] - a[1]], [b[0] - v[0], b[1] - v[1]]);
        let cos = (u[0] * w[0] + u[1] * w[1]) / ((u[0] * u[0] + u[1] * u[1]).sqrt() * (w[0] * w[0] + w[1] * w[1]).sqrt());
        // Interior angle 180° − turn: above 150° when the turn is below 30° (cos > √3/2).
        need.push(match o {
            Sign::Positive if cos > 0.866_025_403_784_438_6 => 2,
            Sign::Positive => 1,
            Sign::Zero => 2,
            Sign::Negative if cos < -0.5 => 3,
            Sign::Negative => 2,
        });
    }
    Ok(need)
}

/// The smallest corner scaled Jacobian of a planar quadrilateral (corner `i`:
/// `cross(p[i+1] − p[i], p[i−1] − p[i]) / (|…| |…|)`), in `[-1, 1]`.
/// A degenerate corner (zero-length edge) scores −1, the worst, so optimizers move away from it.
pub(crate) fn min_scaled_jacobian(p: [[f64; 2]; 4]) -> f64 {
    (0..4).fold(f64::INFINITY, |m, i| m.min(planar_corner(&p, i).1.unwrap_or(-1.0)))
}

/// The smallest corner scaled Jacobian of a recombined quadrilateral (sin 15°): a flatter corner
/// would survive the split as a nearly straight angle, which rounding of the mid nodes can turn
/// degenerate.
const MIN_RECOMBINED_QUALITY: f64 = 0.258_819_045_102_520_74;

/// Greedy recombination of adjacent triangle pairs into strictly convex quadrilaterals, best
/// first by the smallest corner scaled Jacobian (Lee and Lo 1994; Knupp 2001). Returns the
/// quadrilaterals (coarse ids, counterclockwise) and the triangles left over.
///
/// A merge removes one element around the two shared vertices; it is skipped when a boundary
/// vertex would drop below the elements it needs ([`fan_needs`]).
fn recombine(tri: &Triangulation, need: &[usize]) -> Result<(Vec<[u32; 4]>, Vec<[u32; 3]>)> {
    let mut count = vec![0usize; tri.pts.len()];
    for t in &tri.tris {
        for &v in t {
            if let Some(c) = count.get_mut(v as usize) {
                *c += 1;
            }
        }
    }
    let mut edges: Vec<(u32, u32, u32, usize)> = Vec::with_capacity(3 * tri.tris.len());
    for (t, v) in tri.tris.iter().enumerate() {
        for i in 0..3 {
            let (a, b) = (v[(i + 1) % 3], v[(i + 2) % 3]);
            edges.push((a.min(b), a.max(b), t as u32, i));
        }
    }
    edges.sort_unstable();
    let pt = |v: u32| tri.p2.get(v as usize).ok_or_else(|| Error::Failed("point index".into()));
    let xy = |v: u32| tri.pts.get(v as usize).copied().ok_or_else(|| Error::Failed("point index".into()));
    let mut candidates: Vec<(f64, u32, u32, [u32; 4], (u32, u32))> = Vec::new();
    for w in edges.windows(2) {
        let [(a0, b0, t1, i), (a1, b1, t2, j)] = w else { continue };
        if (a0, b0) != (a1, b1) {
            continue;
        }
        let (v1, v2) = match (tri.tris.get(*t1 as usize), tri.tris.get(*t2 as usize)) {
            (Some(x), Some(y)) => (x, y),
            _ => continue,
        };
        let (p, q, r, s) = (v1[*i], v1[(i + 1) % 3], v1[(i + 2) % 3], v2[*j]);
        let quad = [p, q, s, r];
        if !turns_left([pt(p)?, pt(q)?, pt(s)?, pt(r)?]) {
            continue;
        }
        let quality = min_scaled_jacobian([xy(p)?, xy(q)?, xy(s)?, xy(r)?]);
        if quality < MIN_RECOMBINED_QUALITY {
            continue;
        }
        candidates.push((quality, *t1, *t2, quad, (q, r)));
    }
    candidates.sort_by(|x, y| y.0.total_cmp(&x.0).then((x.1, x.2).cmp(&(y.1, y.2))));
    let mut used = vec![false; tri.tris.len()];
    let mut quads = Vec::new();
    for (_, t1, t2, quad, (q, r)) in candidates {
        let (Some(&u1), Some(&u2)) = (used.get(t1 as usize), used.get(t2 as usize)) else { continue };
        if u1 || u2 {
            continue;
        }
        let keeps = |v: u32| count.get(v as usize).copied().unwrap_or(0) > need.get(v as usize).copied().unwrap_or(0);
        if !keeps(q) || !keeps(r) {
            continue;
        }
        for v in [q, r] {
            if let Some(c) = count.get_mut(v as usize) {
                *c -= 1;
            }
        }
        for t in [t1, t2] {
            if let Some(u) = used.get_mut(t as usize) {
                *u = true;
            }
        }
        quads.push(quad);
    }
    let leftovers = tri.tris.iter().zip(&used).filter(|(_, u)| !**u).map(|(t, _)| *t).collect();
    Ok((quads, leftovers))
}

/// The uniform split: every coarse edge gets a mid node (on the boundary, the prescribed point
/// it skipped), every element a centre node; a quadrilateral becomes four quadrilaterals, a
/// triangle three. Node order: the input boundary points (input order), the interior coarse
/// vertices, the interior edge mid nodes, the centres.
fn split(
    region: &Region,
    tri: &Triangulation,
    to_input: &[u32],
    boundary_mid: &BTreeMap<(u32, u32), u32>,
    quads: &[[u32; 4]],
    leftovers: &[[u32; 3]],
) -> Result<(Vec<[f64; 2]>, Vec<[u32; 4]>)> {
    let mut nodes = region.pts.clone();
    let mut final_id = vec![u32::MAX; tri.pts.len()];
    for (c, &input) in to_input.iter().enumerate() {
        if let Some(f) = final_id.get_mut(c) {
            *f = input;
        }
    }
    let xy = |v: u32| tri.pts.get(v as usize).copied().ok_or_else(|| Error::Failed("point index".into()));
    let elements: Vec<&[u32]> = quads.iter().map(|q| q.as_slice()).chain(leftovers.iter().map(|t| t.as_slice())).collect();
    for e in &elements {
        for &v in e.iter() {
            let f = final_id.get_mut(v as usize).ok_or_else(|| Error::Failed("vertex id".into()))?;
            if *f == u32::MAX {
                *f = nodes.len() as u32;
                nodes.push(xy(v)?);
            }
        }
    }
    let fid = |v: u32| final_id.get(v as usize).copied().filter(|&f| f != u32::MAX).ok_or_else(|| Error::Failed("vertex id".into()));
    let mut mids: BTreeMap<(u32, u32), u32> = BTreeMap::new();
    for e in &elements {
        for (i, &a) in e.iter().enumerate() {
            let b = e.get((i + 1) % e.len()).copied().ok_or_else(|| Error::Failed("element".into()))?;
            let key = (a.min(b), a.max(b));
            if let std::collections::btree_map::Entry::Vacant(slot) = mids.entry(key) {
                let id = match boundary_mid.get(&key) {
                    Some(&m) => m,
                    None => {
                        let (pa, pb) = (xy(a)?, xy(b)?);
                        nodes.push([0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])]);
                        nodes.len() as u32 - 1
                    }
                };
                slot.insert(id);
            }
        }
    }
    let mid = |a: u32, b: u32| mids.get(&(a.min(b), a.max(b))).copied().ok_or_else(|| Error::Failed("mid node".into()));
    let mut out = Vec::with_capacity(4 * quads.len() + 3 * leftovers.len());
    for e in &elements {
        let k = e.len();
        let mut c = [0.0, 0.0];
        for &v in e.iter() {
            let p = xy(v)?;
            c = [c[0] + p[0], c[1] + p[1]];
        }
        nodes.push([c[0] / k as f64, c[1] / k as f64]);
        let centre = nodes.len() as u32 - 1;
        for (i, &v) in e.iter().enumerate() {
            let next = e.get((i + 1) % k).copied().ok_or_else(|| Error::Failed("element".into()))?;
            let prev = e.get((i + k - 1) % k).copied().ok_or_else(|| Error::Failed("element".into()))?;
            out.push([fid(v)?, mid(v, next)?, centre, mid(prev, v)?]);
        }
    }
    Ok((nodes, out))
}
