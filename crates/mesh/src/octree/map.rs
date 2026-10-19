//! The topological map from the boundary's feature graph (corners, curves, patches) onto the
//! core's boundary quads, after X. Gao, H. Shen, D. Panozzo, *Feature preserving octree-based
//! hexahedral meshing*, CGF 38(5), 2019, §3.3 (paper only; the authors' code has no licence and
//! is not read). Gao et al. trace the curves first and read the patches off them; here the
//! patches come first and the curves are their borders, which keeps every region on the side of
//! the part it faces:
//!
//! 1. **Label.** Each quad takes the patch its outward normal meets first (within two quad
//!    sizes), or else the patch nearest its centre.
//! 2. **Clean.** A quad whose label differs from all its edge neighbours takes their majority
//!    label; then every patch keeps only its largest connected group of quads, the other
//!    groups joining the label they border most. Both repeat a bounded number of times.
//! 3. **Check.** The map must be a bijection: every patch has a region, and two regions border
//!    each other exactly where the part has a curve between their patches.
//! 4. **Features.** Region borders are the curves (each border edge takes the curve between its
//!    two patches nearest to it); a vertex where three or more regions meet is a corner and takes
//!    the nearest corner of the part around those patches. Corners that are not region junctions
//!    (where one curve turns sharply) go to the nearest vertex on their curve's border.
//!
//! Any failure names where the grid is too coarse for the features; the caller refines there and
//! tries again (§3.3, "Bijective map").

use std::collections::{BTreeMap, BTreeSet};

use super::geom::{Geometry, closest_on_segment, dist2};
use super::smooth::Constraint;
use super::topo::quad_edges;
use crate::vec3::{P, cross, norm, sub};

/// Rounds of label cleaning.
const CLEAN_ROUNDS: usize = 16;

/// Where and why the map failed.
#[derive(Debug, Clone)]
pub(crate) struct MapFailure {
    pub at: Vec<P>,
    pub why: crate::message::Message,
}

/// The map: a constraint per quad vertex and a patch per quad.
#[derive(Debug, Clone)]
pub(crate) struct FeatureMap {
    pub role: BTreeMap<u32, Constraint>,
    pub quad_patch: Vec<u32>,
}

/// Map the feature graph onto `quads` (the core's boundary, outward) over `xyz`.
pub(crate) fn map_features(g: &Geometry, quads: &[[u32; 4]], xyz: &dyn Fn(u32) -> P) -> Result<FeatureMap, MapFailure> {
    let n = quads.len();
    let centre = |q: usize| {
        let mut c = [0.0; 3];
        for &v in &quads[q] {
            let p = xyz(v);
            for k in 0..3 {
                c[k] += 0.25 * p[k];
            }
        }
        c
    };
    let edges = quad_edges(quads);
    let mut nbr: Vec<Vec<u32>> = vec![Vec::new(); n];
    for qs in edges.values() {
        if let [a, b] = qs.as_slice() {
            nbr[*a as usize].push(*b);
            nbr[*b as usize].push(*a);
        }
    }

    // 1. Label: the patch the quad faces (a ray along its outward normal), else the nearest.
    let mut label: Vec<u32> = (0..n)
        .map(|q| {
            let c = centre(q);
            let p = quads[q].map(xyz);
            let nrm = cross(sub(p[2], p[0]), sub(p[3], p[1]));
            let l = norm(nrm);
            let reach = 2.0 * (norm(sub(p[2], p[0])) + norm(sub(p[3], p[1])));
            let hit = if l > 0.0 { g.ray(c, nrm.map(|x| x / l), reach, u32::MAX).map(|(t, _)| t) } else { None };
            hit.or_else(|| g.nearest(c, |_| true).map(|(t, _, _)| t)).and_then(|t| g.tri_patch.get(t as usize).copied()).unwrap_or(0)
        })
        .collect();

    // A patch no quad faces (a narrow face the core's boundary passes by) takes the quad
    // nearest to it, which the cleaning keeps.
    let mut seeded: BTreeSet<usize> = BTreeSet::new();
    let have: BTreeSet<u32> = label.iter().copied().collect();
    for p in 0..g.patches.len() as u32 {
        if have.contains(&p) {
            continue;
        }
        let best = (0..n)
            .filter(|q| !seeded.contains(q))
            .filter_map(|q| g.nearest(centre(q), |t| g.tri_patch.get(t as usize) == Some(&p)).map(|x| (q, x.2)))
            .min_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
        if let Some((q, _)) = best {
            label[q] = p;
            seeded.insert(q);
        }
    }
    // 2. Clean.
    for _ in 0..CLEAN_ROUNDS {
        let mut changed = false;
        for q in 0..n {
            let mut count: BTreeMap<u32, usize> = BTreeMap::new();
            for &m in &nbr[q] {
                *count.entry(label[m as usize]).or_default() += 1;
            }
            let same = count.get(&label[q]).copied().unwrap_or(0);
            if same == 0
                && !nbr[q].is_empty()
                && !seeded.contains(&q)
                && let Some((&best, _)) = count.iter().filter(|(l, _)| **l != label[q]).max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))
            {
                label[q] = best;
                changed = true;
            }
        }
        let mut comp = vec![usize::MAX; n];
        let mut comps: Vec<(u32, Vec<usize>)> = Vec::new();
        for s in 0..n {
            if comp[s] != usize::MAX {
                continue;
            }
            let id = comps.len();
            let mut members = vec![s];
            comp[s] = id;
            let mut i = 0;
            while i < members.len() {
                let q = members[i];
                i += 1;
                for &m in &nbr[q] {
                    let m = m as usize;
                    if comp[m] == usize::MAX && label[m] == label[s] {
                        comp[m] = id;
                        members.push(m);
                    }
                }
            }
            comps.push((label[s], members));
        }
        let mut largest: BTreeMap<u32, usize> = BTreeMap::new();
        for (i, (l, m)) in comps.iter().enumerate() {
            let e = largest.entry(*l).or_insert(i);
            if comps[*e].1.len() < m.len() {
                *e = i;
            }
        }
        for (i, (l, members)) in comps.iter().enumerate() {
            if largest.get(l) == Some(&i) {
                continue;
            }
            let mut count: BTreeMap<u32, usize> = BTreeMap::new();
            for &q in members {
                for &m in &nbr[q] {
                    if comp[m as usize] != i {
                        *count.entry(label[m as usize]).or_default() += 1;
                    }
                }
            }
            if let Some((&best, _)) = count.iter().max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0))) {
                for &q in members {
                    label[q] = best;
                }
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }

    // 3. Check: every patch has a region; borders match curves.
    let present: BTreeSet<u32> = label.iter().copied().collect();
    for (p, patch) in g.patches.iter().enumerate() {
        if !present.contains(&(p as u32)) {
            let at: Vec<P> = patch.triangles.iter().take(64).map(|&t| g.triangle_centre(t)).collect();
            return Err(MapFailure { at, why: crate::message![format!("patch {p} (B-rep faces {:?}) has no region", patch.faces)] });
        }
    }
    let mut curves_between: BTreeMap<(u32, u32), Vec<u32>> = BTreeMap::new();
    for (i, c) in g.curves.iter().enumerate() {
        match c.patches.as_slice() {
            [a, b] => curves_between.entry((*a, *b)).or_default().push(i as u32),
            [a] => curves_between.entry((*a, *a)).or_default().push(i as u32),
            _ => {}
        }
    }
    let mut role: BTreeMap<u32, Constraint> = BTreeMap::new();
    let mut border_pairs: BTreeSet<(u32, u32)> = BTreeSet::new();
    let mut vertex_labels: BTreeMap<u32, BTreeSet<u32>> = BTreeMap::new();
    for (q, quad) in quads.iter().enumerate() {
        for &v in quad {
            vertex_labels.entry(v).or_default().insert(label[q]);
        }
    }
    for ((a, b), qs) in &edges {
        let (q0, q1) = match qs.as_slice() {
            [x, y] => (*x as usize, *y as usize),
            _ => continue,
        };
        let (la, lb) = (label[q0].min(label[q1]), label[q0].max(label[q1]));
        if la == lb {
            continue;
        }
        let Some(cs) = curves_between.get(&(la, lb)) else {
            return Err(MapFailure {
                at: vec![xyz(*a), xyz(*b)],
                why: crate::message![format!("regions of patches {la} and {lb} meet where the part has no curve")],
            });
        };
        border_pairs.insert((la, lb));
        let mid = [0, 1, 2].map(|k| 0.5 * (xyz(*a)[k] + xyz(*b)[k]));
        let c = nearest_curve(g, cs, mid);
        for v in [*a, *b] {
            role.entry(v).or_insert(Constraint::Curve(c));
        }
    }
    for (pair, cs) in &curves_between {
        if pair.0 != pair.1 && !border_pairs.contains(pair) {
            let at: Vec<P> = cs.iter().flat_map(|&c| g.curves[c as usize].vertices.iter().map(|&v| g.xyz(v)).collect::<Vec<_>>()).collect();
            return Err(MapFailure {
                at,
                why: crate::message![format!("the regions of patches {} and {} do not meet along their curve", pair.0, pair.1)],
            });
        }
    }

    // 4. Corners.
    let mut taken: BTreeSet<u32> = BTreeSet::new();
    let junctions: Vec<u32> = vertex_labels.iter().filter(|(_, l)| l.len() >= 3).map(|(v, _)| *v).collect();
    for &v in &junctions {
        let p = xyz(v);
        let labels = &vertex_labels[&v];
        let best = g
            .corners
            .iter()
            .copied()
            .filter(|c| !taken.contains(c) && corner_patches(g, *c).is_superset(labels))
            .min_by(|&a, &b| dist2(g.xyz(a), p).total_cmp(&dist2(g.xyz(b), p)).then(a.cmp(&b)));
        match best {
            Some(c) => {
                taken.insert(c);
                role.insert(v, Constraint::Corner(c));
            }
            None => return Err(MapFailure { at: vec![p], why: crate::message![format!("regions {labels:?} meet where the part has no corner")] }),
        }
    }
    for &c in &g.corners {
        if taken.contains(&c) {
            continue;
        }
        let p = g.xyz(c);
        let cp = corner_patches(g, c);
        let best = role
            .iter()
            .filter(|(_, r)| matches!(r, Constraint::Curve(_)))
            .map(|(v, _)| *v)
            .filter(|v| vertex_labels.get(v).is_some_and(|l| l.is_subset(&cp)))
            .min_by(|&a, &b| dist2(xyz(a), p).total_cmp(&dist2(xyz(b), p)).then(a.cmp(&b)));
        match best {
            Some(v) => {
                role.insert(v, Constraint::Corner(c));
            }
            None => {
                return Err(MapFailure { at: vec![p], why: crate::message!["the feature corner at ", crate::message::point(p), " has no vertex"] });
            }
        }
    }
    for (q, quad) in quads.iter().enumerate() {
        for &v in quad {
            role.entry(v).or_insert(Constraint::Patch(label[q]));
        }
    }
    Ok(FeatureMap { role, quad_patch: label })
}

/// The patches around a corner of the part (on either side of the curves that end there).
fn corner_patches(g: &Geometry, c: u32) -> BTreeSet<u32> {
    g.curves.iter().filter(|cv| cv.vertices.first() == Some(&c) || cv.vertices.last() == Some(&c)).flat_map(|cv| cv.patches.iter().copied()).collect()
}

/// The curve of `cs` nearest to `p`.
fn nearest_curve(g: &Geometry, cs: &[u32], p: P) -> u32 {
    cs.iter().copied().min_by(|&a, &b| curve_distance(g, a, p).total_cmp(&curve_distance(g, b, p)).then(a.cmp(&b))).unwrap_or(0)
}

fn curve_distance(g: &Geometry, c: u32, p: P) -> f64 {
    g.curves.get(c as usize).map_or(f64::INFINITY, |cv| {
        cv.vertices.windows(2).map(|w| dist2(p, closest_on_segment(p, g.xyz(w[0]), g.xyz(w[1])).0)).fold(f64::INFINITY, f64::min)
    })
}
