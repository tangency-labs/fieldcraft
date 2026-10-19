//! The model's geometry names against a design opened again (#605).
//!
//! Definitions name the design's faces, edges and vertices by persistent name
//! ([`fieldcraft_model::Model::name_uses`]). A design rebuilt by another SolidCraft can name the
//! same face differently (SolidCraft 0.16 names the pieces of a split face from their history:
//! `F9:face:4` became `F9:section`), and a name the design no longer has selects nothing: a set
//! comes out empty and a seed is skipped. Opening never lets that pass in silence:
//! - [`rename_stale`] finds a stale face again where the model has a mesh on it: the mesh's
//!   boundary facets tagged with the old name are located on the new design's tessellation, and
//!   the name is rewritten only when every sample lands, within tolerance, on one face whose area
//!   agrees with the facets', and no other tag's mesh, stale or not, lands there (else the sets
//!   on it would grow in silence by a sliver another tag held). The facets are grouped by tag in
//!   one pass; a facet longer than the body is refused, and at most [`MAX_STALE`] names are
//!   looked for (the rest are reported).
//! - The mesh of every name that still resolves is checked the same way (every facet's centre on
//!   a body with a stale face, a sample of them elsewhere): one lying mostly on another face is
//!   reported (`mismatched`), the sign of SolidCraft giving an old name to another face. Edges and vertices follow
//!   their faces: an edge `A|B` becomes the edge between the faces `A` and `B` became, when the
//!   design has exactly that edge. Anything else is refused and reported, never guessed.
//! - [`stale_uses`] lists every use still naming what the design lacks, so the opener names each
//!   one and asks for it to be picked again.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use fieldcraft_geometry::{Design, Shape};
use fieldcraft_model::{ElementFace, Entity, GeometryRef, Model};
use serde_json::{Value, json};

/// A use of a name the design lacks.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Stale {
    pub(crate) owner: String,
    pub(crate) reference: GeometryRef,
}

impl Stale {
    pub(crate) fn message(&self) -> String {
        format!("{} names the {}, which the design no longer has: pick it again", self.owner, self.reference)
    }

    pub(crate) fn json(&self) -> Value {
        json!({"owner": self.owner, "reference": self.reference})
    }
}

/// A name moved to the one the design now gives the same entity.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Renamed {
    pub(crate) from: GeometryRef,
    pub(crate) to: GeometryRef,
    /// Definitions that followed it (the mesh's tags follow too).
    pub(crate) uses: usize,
    /// How it was found: where the mesh lies, or between the faces it separates.
    pub(crate) how: &'static str,
}

impl Renamed {
    pub(crate) fn message(&self) -> String {
        let name = |r: &GeometryRef| r.entity.name().unwrap_or_default().to_string();
        let followed = match self.uses {
            0 => "the mesh's tags follow it".to_string(),
            1 => "1 definition follows it".to_string(),
            n => format!("{n} definitions follow it"),
        };
        format!(
            "{} `{}` of body `{}` is now `{}` in the design (found {}); {followed}: save the model to keep the new name",
            self.from.entity.kind(),
            name(&self.from),
            self.from.body,
            name(&self.to),
            self.how
        )
    }

    pub(crate) fn json(&self) -> Value {
        json!({"from": self.from, "to": self.to, "uses": self.uses, "found": self.how})
    }
}

/// A stale name not rewritten, and why.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Refused {
    pub(crate) reference: GeometryRef,
    pub(crate) why: String,
}

/// The bodies of `design` a reference to `body` reaches: the body, or the cells partitions cut
/// it into.
fn bodies_for<'a>(design: &'a Design, body: &str) -> Vec<&'a Shape> {
    design.bodies.iter().filter(|b| b.name == body || (design.source_of(&b.name) == Some(body) && b.name != body)).collect()
}

/// Whether `name` (with its body and kind taken from `r`) selects `have`.
fn selects(r: &GeometryRef, have: &str) -> bool {
    GeometryRef::new(r.body.clone(), r.entity.clone().with_name(have)).is_ok_and(|g| r.selects(&g))
}

/// The names of a design's bodies, found once per body.
struct Names<'a> {
    design: &'a Design,
    /// Per body: its faces' names (and a merged face's members'), its edges' and its vertices'.
    bodies: BTreeMap<String, [Vec<String>; 3]>,
}

impl<'a> Names<'a> {
    fn new(design: &'a Design) -> Self {
        Names { design, bodies: BTreeMap::new() }
    }

    /// Whether the design has what `r` names.
    fn resolves(&mut self, r: &GeometryRef) -> bool {
        let bodies = bodies_for(self.design, &r.body);
        if r.entity == Entity::Body {
            return !bodies.is_empty() || self.design.source_of(&r.body).is_some();
        }
        let k = match r.entity {
            Entity::Face(_) => 0,
            Entity::Edge(_) => 1,
            _ => 2,
        };
        bodies.iter().any(|b| {
            let names = self.bodies.entry(b.name.clone()).or_insert_with(|| {
                let faces = b.faces().iter().flat_map(|f| std::iter::once(f.name.clone()).chain(f.members.iter().cloned())).collect();
                let edges = b.edges().map(|es| es.into_iter().map(|e| e.name).collect()).unwrap_or_default();
                let vertices = b.vertices().map(|vs| vs.into_iter().map(|v| v.name).collect()).unwrap_or_default();
                [faces, edges, vertices]
            });
            names.get(k).is_some_and(|n| n.iter().any(|have| selects(r, have)))
        })
    }
}

/// Every definition's use of a name `design` lacks on a body it has, in
/// [`Model::name_uses`]'s order (a body the design lacks altogether may be an imported surface's,
/// and is the open's `stale.sets` and `stale.mesh_bodies`).
pub(crate) fn stale_uses(model: &Model, design: &Design) -> Vec<Stale> {
    let mut known: BTreeMap<GeometryRef, bool> = BTreeMap::new();
    let mut names = Names::new(design);
    let mut out = Vec::new();
    for u in model.name_uses().into_iter().filter(|u| !bodies_for(design, &u.reference.body).is_empty()) {
        let ok = *known.entry(u.reference.clone()).or_insert_with(|| names.resolves(&u.reference));
        if !ok {
            out.push(Stale { owner: u.owner, reference: u.reference });
        }
    }
    out
}

/// Samples of a facet: its centre and its corners pulled halfway to it.
fn samples(corners: &[[f64; 3]]) -> Vec<[f64; 3]> {
    let n = corners.len() as f64;
    let c = corners.iter().fold([0.0; 3], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n, a[2] + p[2] / n]);
    std::iter::once(c).chain(corners.iter().map(|p| [0.5 * (p[0] + c[0]), 0.5 * (p[1] + c[1]), 0.5 * (p[2] + c[2])])).collect()
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn triangle_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let n = cross(sub(b, a), sub(c, a));
    0.5 * dot(n, n).sqrt()
}

/// Distance from `p` to triangle `abc` (Ericson, *Real-Time Collision Detection*, 5.1.5).
fn point_triangle(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let (ab, ac, ap) = (sub(b, a), sub(c, a), sub(p, a));
    let at = |u: f64, v: f64| [a[0] + u * ab[0] + v * ac[0], a[1] + u * ab[1] + v * ac[1], a[2] + u * ab[2] + v * ac[2]];
    let dist = |q: [f64; 3]| dot(sub(p, q), sub(p, q)).sqrt();
    let (d1, d2) = (dot(ab, ap), dot(ac, ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return dist(a);
    }
    let bp = sub(p, b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0.0 && d4 <= d3 {
        return dist(b);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return dist(at(d1 / (d1 - d3), 0.0));
    }
    let cp = sub(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 {
        return dist(c);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return dist(at(0.0, d2 / (d2 - d6)));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        let q = [b[0] + w * (c[0] - b[0]), b[1] + w * (c[1] - b[1]), b[2] + w * (c[2] - b[2])];
        return dist(q);
    }
    let denom = 1.0 / (va + vb + vc);
    dist(at(vb * denom, vc * denom))
}

/// Most cells along the longest side of a body's search grid.
const GRID: f64 = 64.0;

/// Most stale names looked for again in one open; the rest are reported stale.
pub(crate) const MAX_STALE: usize = 1_000;

/// Facets of a face tag sampled to check the tag still names the face its mesh lies on, on a
/// body with no stale face (all of them on a body with one: a rename must see every tag).
const CHECKED_FACETS: usize = 64;

fn finite(p: [f64; 3]) -> bool {
    p.iter().all(|x| x.is_finite())
}

/// A body's tessellation, bucketed for "which face is near this point".
struct Located {
    triangles: Vec<([[f64; 3]; 3], u32)>,
    cells: HashMap<[i64; 3], Vec<usize>>,
    cell: f64,
    /// The body's size: no facet of its mesh is longer.
    size: f64,
    /// Area per face index.
    areas: Vec<f64>,
    /// The tessellation's own chordal tolerance.
    tolerance: f64,
    names: Vec<String>,
}

impl Located {
    fn new(shape: &Shape) -> Option<Located> {
        let size = shape.size().ok().filter(|s| s.is_finite() && *s > 0.0)?;
        let tolerance = (size * fieldcraft_geometry::EDGE_TOLERANCE).max(1e-4);
        let mesh = shape.tessellate(tolerance).ok()?;
        let cell = size / GRID;
        let key = |p: [f64; 3]| p.map(|x| (x / cell).floor() as i64);
        let mut triangles = Vec::with_capacity(mesh.triangles.len());
        let mut cells: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
        let mut areas = vec![0.0; shape.faces().len()];
        for (t, &f) in mesh.triangles.iter().zip(&mesh.tri_face) {
            let p = |k: usize| t.get(k).and_then(|&i| mesh.positions.get(i as usize)).copied().filter(|p| finite(*p));
            let (Some(a), Some(b), Some(c)) = (p(0), p(1), p(2)) else { continue };
            if let Some(x) = areas.get_mut(f as usize) {
                *x += triangle_area(a, b, c);
            }
            let (lo, hi) = (
                key([a[0].min(b[0]).min(c[0]), a[1].min(b[1]).min(c[1]), a[2].min(b[2]).min(c[2])]),
                key([a[0].max(b[0]).max(c[0]), a[1].max(b[1]).max(c[1]), a[2].max(b[2]).max(c[2])]),
            );
            let i = triangles.len();
            triangles.push(([a, b, c], f));
            for x in lo[0]..=hi[0] {
                for y in lo[1]..=hi[1] {
                    for z in lo[2]..=hi[2] {
                        cells.entry([x, y, z]).or_default().push(i);
                    }
                }
            }
        }
        Some(Located { triangles, cells, cell, size, areas, tolerance, names: shape.faces().iter().map(|f| f.name.clone()).collect() })
    }

    fn name(&self, face: u32) -> String {
        self.names.get(face as usize).cloned().unwrap_or_default()
    }

    /// The face nearest `p` within `reach`, and its distance; `None` when no face is that close.
    fn nearest(&self, p: [f64; 3], reach: f64) -> Option<(u32, f64)> {
        if !finite(p) || !(reach.is_finite() && reach >= 0.0) {
            return None;
        }
        let r = (reach / self.cell).ceil().clamp(0.0, GRID) as i64;
        let k = p.map(|x| (x / self.cell).floor() as i64);
        let span = |c: i64| c.saturating_sub(r)..=c.saturating_add(r);
        let mut best: Option<(u32, f64)> = None;
        let mut seen = BTreeSet::new();
        for x in span(k[0]) {
            for y in span(k[1]) {
                for z in span(k[2]) {
                    for &i in self.cells.get(&[x, y, z]).into_iter().flatten() {
                        if !seen.insert(i) {
                            continue;
                        }
                        let Some(([a, b, c], f)) = self.triangles.get(i) else { continue };
                        let d = point_triangle(p, *a, *b, *c);
                        if d <= reach && best.is_none_or(|(_, bd)| d < bd) {
                            best = Some((*f, d));
                        }
                    }
                }
            }
        }
        best
    }

    /// How far a sample of a facet may lie from its face: a facet of a curved face stands off it
    /// by its sagitta, the tessellation by its own. `None` for a facet longer than the body (no
    /// mesh of it has one; a hostile file would make every sample search the whole grid).
    fn reach(&self, facet: &[[f64; 3]]) -> Option<f64> {
        let n = facet.len();
        let longest = (0..n).filter_map(|i| Some(sub(*facet.get(i)?, *facet.get((i + 1) % n)?))).map(|d| dot(d, d).sqrt()).fold(0.0, f64::max);
        (longest.is_finite() && longest <= self.size).then_some(0.1 * longest + 2.0 * self.tolerance)
    }
}

fn facet_area(f: &[[f64; 3]]) -> f64 {
    match f {
        [a, b, c] => triangle_area(*a, *b, *c),
        [a, b, c, d] => triangle_area(*a, *b, *c) + triangle_area(*a, *c, *d),
        _ => 0.0,
    }
}

/// The mesh's boundary facets per face tag (their corners), in one pass over the elements.
fn facets_by_tag(model: &Model) -> BTreeMap<u32, Vec<Vec<[f64; 3]>>> {
    let faces: BTreeSet<u32> = model
        .classification()
        .entities()
        .iter()
        .enumerate()
        .filter(|(_, e)| matches!(e.entity, Entity::Face(_)))
        .filter_map(|(i, _)| u32::try_from(i).ok())
        .collect();
    let mut out: BTreeMap<u32, Vec<Vec<[f64; 3]>>> = BTreeMap::new();
    for (&element, tags) in model.classification().elements() {
        for (&face, &t) in &tags.faces {
            if !faces.contains(&t) {
                continue;
            }
            let nodes = model.face_nodes(ElementFace { element, face }).unwrap_or_default();
            // The corners: the first three of a triangle (tri3, tri6), four of a quadrilateral.
            let corners = if matches!(nodes.len(), 3 | 6) { 3 } else { 4 };
            let pts: Vec<[f64; 3]> = nodes.iter().take(corners).filter_map(|&n| model.node(n).map(|n| n.xyz)).filter(|p| finite(*p)).collect();
            if pts.len() == corners {
                out.entry(t).or_default().push(pts);
            }
        }
    }
    out
}

/// A face tag of the mesh that lies on another face of the design than the one its name
/// names: the design may give the name to another face now.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Mismatch {
    pub(crate) reference: GeometryRef,
    /// The faces of the design its mesh lies on.
    pub(crate) lies_on: Vec<String>,
}

impl Mismatch {
    pub(crate) fn message(&self) -> String {
        format!(
            "the mesh on the {} lies on {} of the design: the design may give that name to another face now; check what names it, and remesh",
            self.reference,
            self.lies_on.iter().map(|f| format!("`{f}`")).collect::<Vec<_>>().join(", ")
        )
    }
}

/// Where one body's face tags lie: per tag, the faces its sampled facets' centres land on.
struct Footprints {
    /// Tag → (design face → facets landing there).
    by_tag: BTreeMap<u32, BTreeMap<u32, usize>>,
}

impl Footprints {
    /// `every`: sample every facet (a body with a stale face), else at most [`CHECKED_FACETS`]
    /// per tag, evenly spread.
    fn of(tags: &[u32], facets: &BTreeMap<u32, Vec<Vec<[f64; 3]>>>, at: &Located, every: bool) -> Footprints {
        let mut by_tag = BTreeMap::new();
        for &t in tags {
            let list = facets.get(&t).map(Vec::as_slice).unwrap_or_default();
            let step = if every { 1 } else { list.len().div_ceil(CHECKED_FACETS).max(1) };
            let mut landed: BTreeMap<u32, usize> = BTreeMap::new();
            for f in list.iter().step_by(step) {
                let Some(reach) = at.reach(f) else { continue };
                if let Some(c) = samples(f).first()
                    && let Some((g, _)) = at.nearest(*c, reach)
                {
                    *landed.entry(g).or_default() += 1;
                }
            }
            by_tag.insert(t, landed);
        }
        Footprints { by_tag }
    }

    /// The other tags whose mesh lies on design face `g`.
    fn others_on(&self, g: u32, tag: u32) -> Vec<u32> {
        self.by_tag.iter().filter(|(t, landed)| **t != tag && landed.contains_key(&g)).map(|(t, _)| *t).collect()
    }
}

/// Rewrite the model's stale face, edge and vertex names to the design's names for the same
/// entities, where they are found unambiguously (see the module docs); the mesh's tags follow.
/// Returns the names rewritten, the stale names that were not (with why), and the face tags
/// whose mesh lies on another face than the one their name names.
///
/// `framed` is the same model in the design's millimetres when `model` is the stored one (in
/// SI): the mesh is located from it, the names move in `model`. `None`: `model` is in
/// millimetres itself.
pub(crate) fn rename_stale(model: &mut Model, framed: Option<&Model>, design: &Design) -> (Vec<Renamed>, Vec<Refused>, Vec<Mismatch>) {
    // Every stale name, from the definitions and from the mesh's tags.
    let mut stale: BTreeSet<GeometryRef> = BTreeSet::new();
    let mut names = Names::new(design);
    let candidates = model.name_uses().into_iter().map(|u| u.reference).chain(model.classification().entities().iter().cloned());
    for r in candidates {
        if r.entity != Entity::Body && design.body(&r.body).is_some() && !stale.contains(&r) && !names.resolves(&r) {
            stale.insert(r);
        }
    }
    let mut refused = Vec::new();
    // Past the cap, a name is only reported.
    let over: Vec<GeometryRef> = stale.iter().skip(MAX_STALE).cloned().collect();
    for r in over {
        stale.remove(&r);
        refused.push(Refused { reference: r, why: format!("more than {MAX_STALE} names are stale; pick it again") });
    }
    let entities: Vec<GeometryRef> = model.classification().entities().to_vec();
    let tag_of = |r: &GeometryRef| entities.iter().position(|e| e == r).and_then(|i| u32::try_from(i).ok());
    let facets = facets_by_tag(framed.unwrap_or(model));
    // Per body with face tags: where each face tag's mesh lies.
    let mut found: Vec<(GeometryRef, GeometryRef)> = Vec::new();
    let mut mismatched = Vec::new();
    let bodies: BTreeSet<&str> = entities.iter().filter(|e| matches!(e.entity, Entity::Face(_))).map(|e| e.body.as_str()).collect();
    for body in bodies {
        let Some(shape) = design.body(body) else { continue };
        let tags: Vec<u32> = entities
            .iter()
            .enumerate()
            .filter(|(_, e)| e.body == body && matches!(e.entity, Entity::Face(_)))
            .filter_map(|(i, _)| u32::try_from(i).ok())
            .collect();
        let stale_here: Vec<&GeometryRef> = stale.iter().filter(|r| r.body == body && matches!(r.entity, Entity::Face(_))).collect();
        let Some(at) = Located::new(shape) else {
            refused.extend(stale_here.into_iter().map(|r| Refused { reference: r.clone(), why: "the body could not be tessellated".into() }));
            continue;
        };
        let prints = Footprints::of(&tags, &facets, &at, !stale_here.is_empty());
        // Tags that still resolve: their mesh should lie on the face of that name.
        for &t in &tags {
            let Some(r) = entities.get(t as usize) else { continue };
            let Some(landed) = prints.by_tag.get(&t) else { continue };
            if stale.contains(r) || landed.is_empty() {
                continue;
            }
            let total: usize = landed.values().sum();
            let own: usize = landed.iter().filter(|(g, _)| selects(r, &at.name(**g))).map(|(_, n)| *n).sum();
            if 2 * own < total {
                mismatched.push(Mismatch { reference: r.clone(), lies_on: landed.keys().map(|g| at.name(*g)).collect() });
            }
        }
        for r in stale_here {
            let Some(tag) = tag_of(r) else {
                refused.push(Refused { reference: r.clone(), why: "no mesh lies on it to find it by".into() });
                continue;
            };
            let result = locate_face(facets.get(&tag).map(Vec::as_slice).unwrap_or_default(), &at).and_then(|g| {
                let others = prints.others_on(g, tag);
                match others.first().and_then(|t| entities.get(*t as usize)) {
                    Some(other) => Err(format!("the mesh on the {other} lies on face {} too", at.name(g))),
                    None => GeometryRef::face(body, at.name(g)).map_err(|e| e.to_string()),
                }
            });
            match result {
                Ok(to) => found.push((r.clone(), to)),
                Err(why) => refused.push(Refused { reference: r.clone(), why }),
            }
        }
    }
    // A stale face no mesh is tagged with on a body the loop above did not see.
    for r in stale.iter().filter(|r| matches!(r.entity, Entity::Face(_))) {
        if !found.iter().any(|(f, _)| f == r) && !refused.iter().any(|x| &x.reference == r) {
            refused.push(Refused { reference: r.clone(), why: "no mesh lies on it to find it by".into() });
        }
    }
    // Edges and vertices through their faces.
    let faces: BTreeMap<(String, String), String> =
        found.iter().filter_map(|(f, t)| Some(((f.body.clone(), f.entity.name()?.to_string()), t.entity.name()?.to_string()))).collect();
    for r in stale.iter().filter(|r| matches!(r.entity, Entity::Edge(_) | Entity::Vertex(_))) {
        match through_faces(design, r, &faces) {
            Ok(to) => found.push((r.clone(), to)),
            Err(why) => refused.push(Refused { reference: r.clone(), why }),
        }
    }
    let mut renamed = Vec::new();
    for (from, to) in found {
        let how = if matches!(from.entity, Entity::Face(_)) { "where the mesh lies on it" } else { "between the faces it separates" };
        match model.rename_entity(&from, &to) {
            Ok(uses) => renamed.push(Renamed { from, to, uses, how }),
            Err(e) => refused.push(Refused { reference: from, why: e.to_string() }),
        }
    }
    (renamed, refused, mismatched)
}

/// Relative difference allowed between the area of a face's facets and the face's.
const AREA_TOLERANCE: f64 = 0.05;

/// The design face (index) a stale face tag's `facets` lie on, or why there is no one such face:
/// every sample of every facet on that face, and the facets' area the face's.
fn locate_face(facets: &[Vec<[f64; 3]>], at: &Located) -> Result<u32, String> {
    if facets.is_empty() {
        return Err("no mesh lies on it to find it by".into());
    }
    let mut face: Option<u32> = None;
    let mut area = 0.0;
    for f in facets {
        let reach = at.reach(f).ok_or("its mesh has a facet larger than the body")?;
        for p in samples(f) {
            // Its point is in the geometry's millimetres, which the user may not read: not named.
            let (g, _) = at.nearest(p, reach).ok_or("the mesh on it lies off every face of the design")?;
            match face {
                None => face = Some(g),
                Some(h) if h != g => {
                    return Err(format!("the mesh on it lies on two faces of the design, {} and {}", at.name(h), at.name(g)));
                }
                Some(_) => {}
            }
        }
        area += facet_area(f);
    }
    let g = face.ok_or("no mesh lies on it to find it by")?;
    let want = at.areas.get(g as usize).copied().unwrap_or(0.0);
    if !(want > 0.0) || (area - want).abs() > AREA_TOLERANCE * want {
        // As a share of the face: its areas are in the geometry's mm².
        let share = if want > 0.0 { 100.0 * area / want } else { f64::INFINITY };
        return Err(format!("the mesh on it covers {share:.1} % of the area of face {}: a part of a face, or more than one", at.name(g)));
    }
    Ok(g)
}

/// The edge or vertex the design names after the faces `r`'s faces became (`faces`: stale face
/// → new name, per body), when the design has exactly one such.
fn through_faces(design: &Design, r: &GeometryRef, faces: &BTreeMap<(String, String), String>) -> Result<GeometryRef, String> {
    let shape = design.body(&r.body).ok_or("the design has no such body")?;
    let name = r.entity.name().unwrap_or_default();
    let (base, number) = match name.rsplit_once('~') {
        Some((b, k)) if !k.is_empty() && k.bytes().all(|c| c.is_ascii_digit()) => (b, Some(k)),
        _ => (name, None),
    };
    if number.is_some() {
        return Err("one of several between the same faces (`~k`): pick it again".into());
    }
    let design_faces: BTreeSet<&str> = shape.faces().iter().map(|f| f.name.as_str()).collect();
    let known = |n: &str| -> Option<String> {
        faces.get(&(r.body.clone(), n.to_string())).cloned().or_else(|| design_faces.contains(n).then(|| n.to_string()))
    };
    // A face name may hold `|` itself: every way of cutting the name into known faces.
    let parts: Vec<&str> = base.split('|').collect();
    let mut results: BTreeSet<String> = BTreeSet::new();
    let mut stack: Vec<(usize, Vec<String>)> = vec![(0, Vec::new())];
    let mut budget = 10_000usize;
    while let Some((at, mut got)) = stack.pop() {
        budget = budget.saturating_sub(1);
        if budget == 0 {
            return Err("its name could not be read as faces".into());
        }
        if at == parts.len() {
            if got.len() >= 2 {
                got.sort();
                results.insert(got.join("|"));
            }
            continue;
        }
        for end in at + 1..=parts.len() {
            if let Some(n) = parts.get(at..end).map(|p| p.join("|")).and_then(|n| known(&n)) {
                let mut next = got.clone();
                next.push(n);
                stack.push((end, next));
            }
        }
    }
    let new = match results.len() {
        0 => return Err("a face it lies on was not found again".into()),
        1 => results.into_iter().next().unwrap_or_default(),
        _ => return Err("its name reads as more than one set of faces".into()),
    };
    let names: Vec<String> = match r.entity {
        Entity::Edge(_) => shape.edges().map_err(|e| e.to_string())?.into_iter().map(|e| e.name).collect(),
        _ => shape.vertices().map_err(|e| e.to_string())?.into_iter().map(|v| v.name).collect(),
    };
    if !names.contains(&new) {
        return Err(format!("the design has no single {} `{new}` between those faces", r.entity.kind()));
    }
    GeometryRef::new(r.body.clone(), r.entity.clone().with_name(&new)).map_err(|e| e.to_string())
}

/// What opening a design did to the model's names.
#[derive(Debug, Clone, Default)]
pub(crate) struct Report {
    pub(crate) renamed: Vec<Renamed>,
    pub(crate) refused: Vec<Refused>,
    /// Uses still naming what the design lacks, after the renames.
    pub(crate) stale: Vec<Stale>,
    /// Face tags whose mesh lies on another face than their name's.
    pub(crate) mismatched: Vec<Mismatch>,
}

impl Report {
    /// One line per rename and per stale use (with why it was not found again).
    pub(crate) fn warnings(&self) -> Vec<String> {
        let why = |r: &GeometryRef| self.refused.iter().find(|x| &x.reference == r).map(|x| x.why.as_str());
        self.renamed
            .iter()
            .map(Renamed::message)
            .chain(self.stale.iter().map(|s| match why(&s.reference) {
                Some(w) => format!("{} (not found again: {w})", s.message()),
                None => s.message(),
            }))
            .chain(self.mismatched.iter().map(Mismatch::message))
            .collect()
    }

    pub(crate) fn json(&self) -> Value {
        json!({
            "renamed": self.renamed.iter().map(Renamed::json).collect::<Vec<_>>(),
            "references": self.stale.iter().map(|s| {
                let mut v = s.json();
                if let Some(r) = self.refused.iter().find(|x| x.reference == s.reference) {
                    v["why"] = json!(r.why);
                }
                v
            }).collect::<Vec<_>>(),
            "mismatched": self.mismatched.iter().map(|m| json!({"reference": m.reference, "lies_on": m.lies_on})).collect::<Vec<_>>(),
        })
    }

    pub(crate) fn changed_model(&self) -> bool {
        !self.renamed.is_empty()
    }
}

/// Bring the model's names to `design` (the design as opened, `raw`, cut by the model's
/// partitions): stale names found again are rewritten ([`rename_stale`]), the partitions, merges
/// and splits applied again when a name they hold moved, and the uses still stale listed.
///
/// `model` is the stored model (SI) and `framed` the same in millimetres, as the geometry is: the
/// mesh is located and the design cut by `framed`, the names move in `model`, and `framed` is
/// made from it again after a rename. `None`: `model` is in millimetres itself (a command run in
/// the geometry's frame).
pub(crate) fn reconcile(model: &mut Model, framed: Option<&mut Model>, raw: &Design, design: Design) -> crate::Result<(Design, Report)> {
    let (renamed, refused, mismatched) = rename_stale(model, framed.as_deref(), &design);
    let design = if renamed.is_empty() {
        design
    } else if let Some(f) = framed {
        *f = super::file::frame_geometry(model).map_err(crate::EngineError::Params)?;
        super::datum::partitioned(raw, f)?
    } else {
        super::datum::partitioned(raw, model)?
    };
    let stale = stale_uses(model, &design);
    Ok((design, Report { renamed, refused, stale, mismatched }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Session;

    /// Models FieldCraft 0.7.0 wrote on the bracket SolidCraft 0.11.0 saved (`tests/stale-names`):
    /// set `cut` on `F9:face:4` of `Bracket.1` with a pressure load on it and a seed on the edge
    /// `F2:side:l4|F9:face:4`, meshed or not. SolidCraft 0.16 names that face `F9:section`.
    fn open(file: &str) -> (Session, Value) {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/stale-names").join(file);
        let mut s = Session::default();
        let r = s.execute("file.open", &json!({"path": path.to_string_lossy()})).unwrap();
        (s, r)
    }

    fn face(n: &str) -> GeometryRef {
        GeometryRef::face("Bracket.1", n).unwrap()
    }

    fn renamed(r: &Value) -> Vec<(String, String)> {
        let name = |v: &Value| ["face", "edge", "vertex"].iter().find_map(|k| v[*k].as_str()).unwrap_or_default().to_string();
        r["renamed"].as_array().unwrap().iter().map(|x| (name(&x["from"]), name(&x["to"]))).collect()
    }

    #[test]
    fn a_0_7_0_model_follows_its_renamed_faces_where_the_mesh_lies() {
        let (mut s, r) = open("meshed.fieldcraft");
        let moved = renamed(&r);
        assert!(moved.contains(&("F9:face:4".into(), "F9:section".into())), "{r}");
        assert!(moved.contains(&("F2:side:l4|F9:face:4".into(), "F2:side:l4|F9:section".into())), "{r}");
        assert_eq!(r["stale"], json!({"mesh_bodies": [], "sets": [], "seeds": [], "references": []}), "{r}");
        assert_eq!(r["mismatched"], json!([]), "{r}");
        let warnings = r["warnings"].to_string();
        assert!(warnings.contains("`F9:face:4` of body `Bracket.1` is now `F9:section`"), "{warnings}");
        assert_eq!(s.model.sets().get("cut").map(|c| c.geometry.clone()), Some(vec![face("F9:section")]));
        assert!(s.is_dirty(), "the rename is a change to save");
        // The remesh keeps the seed and the set, and says nothing was skipped.
        let m = s.execute("mesh.tet.generate", &json!({"body": "Bracket.1", "size": 3})).unwrap();
        assert_eq!(m["bodies"][0]["seeds"][0]["mesh_edges"], 7, "{m}");
        assert!(m.get("warnings").is_none(), "{m}");
        let check = s.execute("model.check", &json!({})).unwrap();
        assert!(!check["diagnostics"].to_string().contains("set `cut`"), "{check}");
    }

    #[test]
    fn a_0_7_0_model_without_a_mesh_asks_for_its_stale_names_to_be_picked_again() {
        let (mut s, r) = open("unmeshed.fieldcraft");
        assert_eq!(r["renamed"], json!([]), "{r}");
        let owners: Vec<&str> = r["stale"]["references"].as_array().unwrap().iter().filter_map(|x| x["owner"].as_str()).collect();
        assert_eq!(owners, ["set `cut`", "edge seed"], "{r}");
        let warnings = r["warnings"].to_string();
        assert!(warnings.contains("face F9:face:4 of body Bracket.1, which the design no longer has: pick it again"), "{warnings}");
        assert!(!s.is_dirty());
        let m = s.execute("mesh.tet.generate", &json!({"body": "Bracket.1", "size": 3})).unwrap();
        assert!(m["warnings"].to_string().contains("edge seed was not applied"), "{m}");
        let check = s.execute("model.check", &json!({})).unwrap().to_string();
        assert!(check.contains("the mesh has no face F9:face:4 of body Bracket.1"), "{check}");
    }

    #[test]
    fn a_stale_face_the_mesh_does_not_cover_alone_is_refused() {
        let (mut s, _) = open("meshed.fieldcraft");
        let design = s.design.clone().unwrap();
        // The mesh on F9:section back under its old name, with the facets of F7:section on it too.
        s.model.rename_entity(&face("F9:section"), &face("F9:face:4")).unwrap();
        let (old, other) = (s.model.geometry_tag(&face("F9:face:4")).unwrap(), s.model.geometry_tag(&face("F7:section")).unwrap());
        let on_other: Vec<(u32, u8)> = s
            .model
            .classification()
            .elements()
            .iter()
            .flat_map(|(&e, t)| t.faces.iter().filter(|(_, g)| **g == other).map(move |(&f, _)| (e, f)))
            .collect();
        for (e, f) in on_other {
            s.model.tag_element_face(e, f, old).unwrap();
        }
        let (renamed, refused, _) = {
            let framed = super::super::file::frame_geometry(&s.model).unwrap();
            rename_stale(&mut s.model, Some(&framed), &design)
        };
        assert!(renamed.iter().all(|r| r.from != face("F9:face:4")), "{renamed:?}");
        let why = refused.iter().find(|r| r.reference == face("F9:face:4")).map(|r| r.why.clone()).unwrap_or_default();
        assert!(why.contains("two faces of the design"), "{why}");
        // Half of a face is not the face: its area gives it away.
        let (mut s, _) = open("meshed.fieldcraft");
        s.model.rename_entity(&face("F9:section"), &face("F9:face:4")).unwrap();
        let (old, half) = (s.model.geometry_tag(&face("F9:face:4")).unwrap(), s.model.geometry_tag(&face("Gone")).unwrap());
        let on_old: Vec<(u32, u8)> = s
            .model
            .classification()
            .elements()
            .iter()
            .flat_map(|(&e, t)| t.faces.iter().filter(|(_, g)| **g == old).map(move |(&f, _)| (e, f)))
            .collect();
        for &(e, f) in on_old.iter().step_by(2) {
            s.model.tag_element_face(e, f, half).unwrap();
        }
        let (renamed, refused, _) = {
            let framed = super::super::file::frame_geometry(&s.model).unwrap();
            rename_stale(&mut s.model, Some(&framed), &design)
        };
        assert!(renamed.iter().all(|r| r.from != face("F9:face:4") && r.from != face("Gone")), "{renamed:?}");
        for r in [face("F9:face:4"), face("Gone")] {
            let why = refused.iter().find(|x| x.reference == r).map(|x| x.why.clone()).unwrap_or_default();
            assert!(why.contains("a part of a face"), "{r}: {why}");
        }
    }

    /// The facets tagged with a face, by tag.
    fn facets_of(s: &Session, r: &GeometryRef) -> Vec<(u32, u8)> {
        let tag = s.model.classification().entities().iter().position(|e| e == r).unwrap() as u32;
        s.model.classification().elements().iter().flat_map(|(&e, t)| t.faces.iter().filter(|(_, g)| **g == tag).map(move |(&f, _)| (e, f))).collect()
    }

    /// A stale face whose mesh is most of a face, the rest of it under another tag: renaming it
    /// would grow the sets on it in silence (the 5 % area tolerance lets it through), so another
    /// tag's mesh on the face refuses it.
    #[test]
    fn a_face_another_tag_also_lies_on_is_not_taken() {
        let (mut s, _) = open("meshed.fieldcraft");
        let design = s.design.clone().unwrap();
        s.model.rename_entity(&face("F9:section"), &face("F9:face:4")).unwrap();
        let on = facets_of(&s, &face("F9:face:4"));
        let other = s.model.geometry_tag(&face("F7:section")).unwrap();
        // A sliver: 3 % of the facets.
        for &(e, f) in on.iter().take(on.len() * 3 / 100 + 1) {
            s.model.tag_element_face(e, f, other).unwrap();
        }
        let (renamed, refused, _) = {
            let framed = super::super::file::frame_geometry(&s.model).unwrap();
            rename_stale(&mut s.model, Some(&framed), &design)
        };
        assert!(renamed.iter().all(|r| r.from != face("F9:face:4")), "{renamed:?}");
        let why = refused.iter().find(|r| r.reference == face("F9:face:4")).map(|r| r.why.clone()).unwrap_or_default();
        assert!(why.contains("F7:section of body Bracket.1 lies on face F9:section too"), "{why}");
    }

    /// Names that still exist but now name another face (SolidCraft reusing a name) are reported.
    #[test]
    fn a_name_now_on_another_face_is_reported() {
        let (mut s, _) = open("meshed.fieldcraft");
        let design = s.design.clone().unwrap();
        // Swap the tags of F9:section and F7:section.
        s.model.rename_entity(&face("F9:section"), &face("Tmp")).unwrap();
        s.model.rename_entity(&face("F7:section"), &face("F9:section")).unwrap();
        s.model.rename_entity(&face("Tmp"), &face("F7:section")).unwrap();
        let (renamed, _, mismatched) = {
            let framed = super::super::file::frame_geometry(&s.model).unwrap();
            rename_stale(&mut s.model, Some(&framed), &design)
        };
        assert!(renamed.is_empty(), "{renamed:?}");
        let names: Vec<(GeometryRef, Vec<String>)> = mismatched.into_iter().map(|m| (m.reference, m.lies_on)).collect();
        assert!(names.contains(&(face("F9:section"), vec!["F7:section".to_string()])), "{names:?}");
        assert!(names.contains(&(face("F7:section"), vec!["F9:section".to_string()])), "{names:?}");
    }

    /// A hostile file's mesh: a facet longer than the body is refused before it is searched for,
    /// and points far out or not finite find nothing (and overflow nothing).
    #[test]
    fn hostile_facets_cost_nothing_and_find_nothing() {
        let (s, _) = open("meshed.fieldcraft");
        let design = s.design.clone().unwrap();
        let at = Located::new(design.body("Bracket.1").unwrap()).unwrap();
        let huge = vec![vec![[0.0, 0.0, 0.0], [1e12, 0.0, 0.0], [0.0, 1e12, 0.0]]];
        assert!(locate_face(&huge, &at).unwrap_err().contains("larger than the body"));
        for p in [[f64::NAN, 0.0, 0.0], [f64::INFINITY, 0.0, 0.0], [9e18, -9e18, 9e18], [1e300, 0.0, 0.0]] {
            assert_eq!(at.nearest(p, 1.0), None, "{p:?}");
        }
        assert_eq!(at.nearest([0.0; 3], f64::NAN), None);
    }
}
