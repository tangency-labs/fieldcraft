//! Selection by picking in the viewport.
//!
//! A click picks what is under the pointer through the selection filter: a geometry face, edge,
//! vertex or body, or a mesh node, element face or element. A box drag picks every visible node,
//! element face or element inside it. What is picked is kept in the model's persistent terms
//! ([`Item`]): geometry by name ([`GeometryRef`], body and face name), so a selection survives a
//! remesh, and the mesh by label. A [`Selection`]'s members are a set definition
//! ([`Selection::params`]), which is how command forms and `model.set.create` take it.
//!
//! Picking runs on the CPU against the display data the viewport draws (the GPU draws the same
//! triangles, so both renderers show what a pick finds): a ray through the pointer's pixel finds
//! the first face it hits with the watertight test of [`fieldcraft_render::pick`]; nodes, vertices
//! and edges are the nearest on screen within a few pixels whose own ray reaches them before any
//! face. Every choice is the nearest, ties going to the lower index, so a pick is deterministic.
//!
//! Geometry is named through the mesh's geometry tags ([`fieldcraft_model::Classification`]) once
//! the model is meshed, and through the drawn bodies' face names before. A design's edges are
//! picked on its drawn B-rep edges by their persistent names (`<face>|<face>`, before and after
//! meshing); other edges and vertices have names only on a mesh that carries them (block and
//! planar meshes).

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_model::{ElementEdge, ElementFace, Entity, GeometryRef, Model, Set, check_name};
use fieldcraft_render::pick::{DepthMap, Projector, segment_distance};
use fieldcraft_render::{DisplayMesh, Rgb, Scene, SceneLine, ScenePoint, SceneTri, Vec3};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

/// Most items a selection holds; a box or a list beyond it is cut there (and says so).
pub const MAX_SELECTION: usize = 100_000;
/// How far (pixels) from a node or vertex, and from an edge, a click still picks it.
const POINT_RADIUS: f64 = 8.0;
const EDGE_RADIUS: f64 = 6.0;
/// The nearest candidates on screen whose visibility is checked with a ray of their own.
const VISIBILITY_CHECKS: usize = 16;
/// Width of a highlighted edge and size of a highlighted node (points).
const EDGE_WIDTH: f64 = 4.0;
const NODE_SIZE: f64 = 9.0;
/// Width of the outline around highlighted faces (points).
const OUTLINE_WIDTH: f64 = 2.5;
/// How far a point may lie behind the first face along its ray and still be on it, as a fraction
/// of the model's diagonal, and at least (metres).
const SLACK: f64 = 1e-4;
const MIN_SLACK: f64 = 1e-9; // unit-factors: a floor in metres (SI)

/// What a pick selects.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Filter {
    /// Geometry faces.
    #[default]
    Faces,
    Edges,
    Vertices,
    Bodies,
    /// Mesh nodes.
    Nodes,
    ElementFaces,
    Elements,
}

impl Filter {
    pub const ALL: [Filter; 7] =
        [Filter::Faces, Filter::Edges, Filter::Vertices, Filter::Bodies, Filter::Nodes, Filter::ElementFaces, Filter::Elements];

    pub fn label(self) -> &'static str {
        match self {
            Filter::Faces => "Faces",
            Filter::Edges => "Edges",
            Filter::Vertices => "Vertices",
            Filter::Bodies => "Bodies",
            Filter::Nodes => "Nodes",
            Filter::ElementFaces => "Element faces",
            Filter::Elements => "Elements",
        }
    }

    /// The filter that picks `item`.
    pub fn of(item: &Item) -> Filter {
        match item {
            Item::Geometry(g) => match g.entity {
                Entity::Body => Filter::Bodies,
                Entity::Face(_) => Filter::Faces,
                Entity::Edge(_) => Filter::Edges,
                Entity::Vertex(_) => Filter::Vertices,
            },
            Item::Node(_) => Filter::Nodes,
            Item::Element(_) => Filter::Elements,
            Item::Face(_) => Filter::ElementFaces,
        }
    }

    /// Geometry (named) rather than mesh (labelled) entities.
    pub fn geometry(self) -> bool {
        matches!(self, Filter::Faces | Filter::Edges | Filter::Vertices | Filter::Bodies)
    }

    /// Whether a box drag selects with this filter.
    pub fn boxes(self) -> bool {
        matches!(self, Filter::Nodes | Filter::ElementFaces | Filter::Elements)
    }
}

/// One selected entity, in the model's persistent terms.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Item {
    Geometry(GeometryRef),
    Node(u32),
    Element(u32),
    Face(ElementFace),
}

impl Item {
    /// The item in words: `face F2:top of body Body1`, `node 5`, `element 3 face 2`.
    pub fn describe(&self) -> String {
        match self {
            Item::Geometry(g) => g.to_string(),
            Item::Node(n) => format!("node {n}"),
            Item::Element(e) => format!("element {e}"),
            Item::Face(f) => format!("element {} face {}", f.element, f.face),
        }
    }

    /// The item as JSON: a geometry reference (`{"body", "face"}`), `{"node"}`, `{"element"}` or
    /// `{"element", "face"}`.
    pub fn to_json(&self) -> Value {
        match self {
            Item::Geometry(g) => json!(g),
            Item::Node(n) => json!({"node": n}),
            Item::Element(e) => json!({"element": e}),
            Item::Face(f) => json!(f),
        }
    }
}

/// How picked items change the selection.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Mode {
    /// The picked items become the selection (a plain click).
    #[default]
    Replace,
    /// Added (Shift).
    Add,
    Remove,
    /// Added, or removed when already selected (Cmd/Ctrl).
    Toggle,
}

/// The selected items, sorted and without repeats, at most [`MAX_SELECTION`].
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Selection {
    items: BTreeSet<Item>,
    /// Items were left out to stay within [`MAX_SELECTION`].
    capped: bool,
}

impl Selection {
    pub fn items(&self) -> &BTreeSet<Item> {
        &self.items
    }

    /// The one filter that picks every selected item (`None` when empty or mixed).
    pub fn filter(&self) -> Option<Filter> {
        let mut kinds = self.items.iter().map(Filter::of);
        let first = kinds.next()?;
        kinds.all(|k| k == first).then_some(first)
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    pub fn capped(&self) -> bool {
        self.capped
    }

    pub fn clear(&mut self) {
        self.items.clear();
        self.capped = false;
    }

    /// Change the selection by `items` as `mode` says; anything past [`MAX_SELECTION`] is left
    /// out and marks the selection capped. Returns whether the items changed.
    pub fn apply(&mut self, items: impl IntoIterator<Item = Item>, mode: Mode) -> bool {
        let mut changed = false;
        if mode == Mode::Replace && !self.is_empty() {
            // A replacement counts as a change, even by the same items.
            changed = true;
            self.clear();
        }
        for item in items {
            changed |= match mode {
                Mode::Remove => self.items.remove(&item),
                Mode::Toggle if self.items.remove(&item) => true,
                _ if self.items.contains(&item) => false,
                _ if self.items.len() >= MAX_SELECTION => {
                    self.capped = true;
                    false
                }
                _ => self.items.insert(item),
            };
        }
        changed
    }

    /// Keep only the items `keep` accepts.
    pub fn retain(&mut self, keep: impl FnMut(&Item) -> bool) {
        self.items.retain(keep);
    }

    /// The selection as a set definition's members (`model.set.create`'s parameters without the
    /// name): `geometry`, `nodes`, `elements` and `faces`, each only when it has members.
    pub fn params(&self) -> Value {
        let mut set = Set::default();
        for item in &self.items {
            match item {
                Item::Geometry(g) => set.geometry.push(g.clone()),
                Item::Node(n) => set.nodes.push(*n),
                Item::Element(e) => set.elements.push(*e),
                Item::Face(f) => set.faces.push(*f),
            }
        }
        // A set's JSON form leaves out its empty lists.
        serde_json::to_value(set).unwrap_or_else(|_| json!({}))
    }

    /// What is selected in a few words: `2 faces, 15 nodes`.
    pub fn summary(&self) -> String {
        let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
        for item in &self.items {
            let kind = match item {
                Item::Geometry(g) => g.entity.kind(),
                Item::Node(_) => "node",
                Item::Element(_) => "element",
                Item::Face(_) => "element face",
            };
            *counts.entry(kind).or_default() += 1;
        }
        if counts.is_empty() {
            return "nothing selected".into();
        }
        let plural = |kind: &str, n: usize| match (kind, n) {
            (k, 1) => k.to_string(),
            ("body", _) => "bodies".into(),
            ("vertex", _) => "vertices".into(),
            (k, _) => format!("{k}s"),
        };
        counts.iter().map(|(k, n)| format!("{n} {}", plural(k, *n))).collect::<Vec<_>>().join(", ")
    }
}

/// The members a selection is made of: a set definition's (`geometry`, `nodes`, `elements`,
/// `faces`; element edges cannot be picked yet).
const MEMBERS: [&str; 4] = ["geometry", "nodes", "elements", "faces"];

/// Read set members (`geometry`, `nodes`, `elements`, `faces`, as [`Selection::params`] writes
/// them) as items, checked as `model.set.create` checks a set. Each list holds at most
/// [`MAX_SELECTION`] entries.
pub fn items_from_params(p: &Value) -> Result<Vec<Item>, String> {
    let mut members = serde_json::Map::new();
    for key in MEMBERS {
        match p.get(key) {
            None | Some(Value::Null) => {}
            Some(Value::Array(a)) if a.len() <= MAX_SELECTION => {
                members.insert(key.into(), Value::Array(a.clone()));
            }
            Some(_) => return Err(format!("`{key}` must be a list of at most {MAX_SELECTION} entries")),
        }
    }
    if members.is_empty() {
        return Ok(Vec::new());
    }
    let set: Set = serde_json::from_value(Value::Object(members)).map_err(|e| e.to_string())?;
    // Empty lists are no members; anything else is checked as a set's members are.
    if set.is_empty() {
        return Ok(Vec::new());
    }
    let set = set.normalized().map_err(|e| e.to_string())?;
    let mut items: Vec<Item> = set.geometry.into_iter().map(Item::Geometry).collect();
    items.extend(set.nodes.into_iter().map(Item::Node));
    items.extend(set.elements.into_iter().map(Item::Element));
    items.extend(set.faces.into_iter().map(Item::Face));
    Ok(items)
}

/// A drawn geometry surface's names: its body and, by face tag, its faces; by edge index, its
/// B-rep edges ([`DisplayMesh::curves`]); its B-rep vertices with their points.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BodyNames {
    pub name: String,
    pub faces: Vec<String>,
    pub edges: Vec<String>,
    pub vertices: Vec<(String, Vec3)>,
}

impl BodyNames {
    /// The reference to face `face` of this body, when it has a name.
    fn face(&self, face: u32) -> Option<GeometryRef> {
        GeometryRef::face(self.name.clone(), self.faces.get(face as usize)?.clone()).ok()
    }

    /// The reference to B-rep vertex `name` of this body.
    fn vertex(&self, name: &str) -> Option<GeometryRef> {
        GeometryRef::new(self.name.clone(), Entity::Vertex(name.to_string())).ok()
    }

    /// The reference to B-rep edge `edge` of this body, when it has a name.
    pub fn edge(&self, edge: u32) -> Option<GeometryRef> {
        GeometryRef::new(self.name.clone(), Entity::Edge(self.edges.get(edge as usize)?.clone())).ok()
    }
}

/// Where the mesh's drawn primitives lie on the geometry: built once per model revision from the
/// geometry tags.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Index {
    /// Per element face triangle: the geometry face it lies on (a boundary face's tag, or the
    /// face a shell fills).
    face_tag: Vec<Option<u32>>,
    /// Per element face triangle: the entity its element fills (a body for a solid).
    region_tag: Vec<Option<u32>>,
    /// Per line element: the entity it fills (an edge).
    line_tag: Vec<Option<u32>>,
    /// Element edges on geometry edges: their corner nodes ([`edge_key`]) → tag.
    edge_tag: BTreeMap<(u32, u32), u32>,
    /// Nodes on geometry vertices, and the vertex's tag.
    vertices: Vec<(u32, u32)>,
}

/// An edge's key, whichever way round its corner nodes come.
fn edge_key(a: u32, b: u32) -> (u32, u32) {
    (a.min(b), a.max(b))
}

impl Index {
    pub fn new(model: &Model, display: &DisplayMesh) -> Index {
        let c = model.classification();
        let dimension = |tag: u32| c.entity(tag).map(|e| e.entity.dimension());
        let tags = |element: u32| c.elements().get(&element);
        let region = |element: u32| tags(element).and_then(|t| t.region);
        let face_tag = display
            .tris
            .iter()
            .map(|t| tags(t.element).and_then(|g| g.faces.get(&t.face).copied().or(g.region.filter(|r| dimension(*r) == Some(2)))))
            .collect();
        let region_tag = display.tris.iter().map(|t| region(t.element)).collect();
        let line_tag = display.lines.iter().map(|l| region(l.element).filter(|r| dimension(*r) == Some(1))).collect();
        let mut edge_tag = BTreeMap::new();
        for (&element, t) in c.elements() {
            for (&edge, &tag) in &t.edges {
                if let Some(nodes) = model.edge_nodes(ElementEdge { element, edge })
                    && let (Some(&a), Some(&b)) = (nodes.first(), nodes.get(1))
                {
                    edge_tag.entry(edge_key(a, b)).or_insert(tag);
                }
            }
        }
        let vertices = c.nodes().iter().filter(|(_, t)| dimension(**t) == Some(0)).map(|(n, t)| (*n, *t)).collect();
        Index { face_tag, region_tag, line_tag, edge_tag, vertices }
    }
}

/// What picking looks at: the drawn display data and the names behind it.
pub struct Picker<'a> {
    /// What the viewport draws: the model's display data, or a results plot's.
    pub display: &'a DisplayMesh,
    pub names: &'a [BodyNames],
    pub model: &'a Model,
    pub index: &'a Index,
    /// Element faces are drawn (and picked).
    pub elements: bool,
    /// Geometry surfaces are drawn (and picked).
    pub geometry: bool,
    /// The geometry's B-rep edges are drawn (and picked), over the mesh too.
    pub curves: bool,
}

/// The candidate nearest on screen within `radius` whose point passes `visible`: candidates are
/// (index, screen distance, the point to check); the first of the nearest few that passes wins.
fn nearest_passing(mut near: Vec<(usize, f64, Vec3)>, radius: f64, visible: impl Fn(Vec3) -> bool) -> Option<usize> {
    near.retain(|c| c.1 <= radius);
    near.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    near.into_iter().take(VISIBILITY_CHECKS).find(|c| visible(c.2)).map(|c| c.0)
}

/// The first drawn face along a ray.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Surface {
    /// An element face triangle ([`DisplayMesh::tris`]) and the distance.
    Element(usize, f64),
    /// A geometry surface triangle ([`DisplayMesh::surface_tris`]).
    Geometry(usize, f64),
}

impl Surface {
    fn t(self) -> f64 {
        match self {
            Surface::Element(_, t) | Surface::Geometry(_, t) => t,
        }
    }
}

impl Picker<'_> {
    fn entity(&self, tag: u32) -> Option<&GeometryRef> {
        self.model.classification().entity(tag)
    }

    /// The first drawn face a ray hits; an element face wins a tie with a geometry face.
    fn surface(&self, o: Vec3, d: Vec3) -> Option<Surface> {
        let e = self.elements.then(|| self.display.hit_element_face(o, d)).flatten().map(|h| Surface::Element(h.index, h.t));
        let g = self.geometry.then(|| self.display.hit_surface(o, d)).flatten().map(|h| Surface::Geometry(h.index, h.t));
        match (e, g) {
            (Some(e), Some(g)) => Some(if g.t() < e.t() { g } else { e }),
            (e, g) => e.or(g),
        }
    }

    /// The drawn faces' triangles (for a depth map).
    fn drawn(&self) -> impl Iterator<Item = [Vec3; 3]> + '_ {
        let e = self.display.tris.iter().filter(|_| self.elements).map(|t| t.p);
        e.chain(self.display.surface_tris.iter().filter(|_| self.geometry).map(|t| t.p))
    }

    /// How far a point may lie behind the first face along its ray and still be on it.
    fn slack(&self) -> f64 {
        let d = self.display.bounds.diagonal();
        if d.is_finite() { (d * SLACK).max(MIN_SLACK) } else { MIN_SLACK }
    }

    /// Whether `p` is in front: the ray through its pixel reaches it before any drawn face.
    fn visible(&self, proj: &Projector, p: Vec3) -> bool {
        self.visible_within(proj, p, self.slack())
    }

    /// Whether `p` is in front, counting it on the first face along its ray within `slack`.
    fn visible_within(&self, proj: &Projector, p: Vec3, slack: f64) -> bool {
        let Some([x, y, _]) = proj.screen(p) else { return false };
        let (o, d) = proj.ray(x, y);
        let along = (p - o).dot(d);
        self.surface(o, d).is_none_or(|s| along <= s.t() + slack)
    }

    /// The candidate nearest on screen to `at` within `radius` that is visible ([`nearest_passing`]).
    fn nearest_visible(&self, proj: &Projector, near: Vec<(usize, f64, Vec3)>, radius: f64) -> Option<usize> {
        nearest_passing(near, radius, |p| self.visible(proj, p))
    }

    /// Whether a point of a B-rep edge is in front. The edge's polyline and the drawn faces are
    /// two approximations of the same curve, which can lie apart by the tessellation's chord
    /// error: a point still counts as on the first face within a hundredth of the model's size.
    fn curve_visible(&self, proj: &Projector, p: Vec3) -> bool {
        self.visible_within(proj, p, self.slack() * 100.0)
    }

    /// The design edge (a [`DisplayMesh::curves`] segment's) nearest `at` on screen, visible.
    fn pick_curve(&self, proj: &Projector, at: [f64; 2]) -> Option<GeometryRef> {
        if !self.curves {
            return None;
        }
        let near = self.display.curves.iter().enumerate().filter_map(|(i, c)| {
            let (sa, sb) = (proj.screen(c.a)?, proj.screen(c.b)?);
            let (distance, k) = segment_distance(at, [sa[0], sa[1]], [sb[0], sb[1]]);
            (distance <= EDGE_RADIUS).then(|| (i, distance, c.a + (c.b - c.a) * k))
        });
        let i = nearest_passing(near.collect(), EDGE_RADIUS, |p| self.curve_visible(proj, p))?;
        let c = self.display.curves.get(i)?;
        self.names.get(c.body as usize)?.edge(c.edge)
    }

    /// The geometry face of an element face triangle or a geometry triangle.
    fn geometry_face(&self, s: Surface) -> Option<GeometryRef> {
        match s {
            Surface::Element(i, _) => self.entity((*self.index.face_tag.get(i)?)?).cloned(),
            Surface::Geometry(i, _) => {
                let t = self.display.surface_tris.get(i)?;
                self.names.get(t.body as usize)?.face(t.face)
            }
        }
    }

    /// The body an element face triangle or a geometry triangle belongs to.
    fn body(&self, s: Surface) -> Option<GeometryRef> {
        let name = match s {
            Surface::Element(i, _) => {
                let tag = self.index.region_tag.get(i).copied().flatten().or(self.index.face_tag.get(i).copied().flatten())?;
                self.entity(tag)?.body.clone()
            }
            Surface::Geometry(i, _) => self.names.get(self.display.surface_tris.get(i)?.body as usize)?.name.clone(),
        };
        GeometryRef::body(name).ok()
    }

    /// What is under pixel `at` with `filter`.
    pub fn pick(&self, proj: &Projector, at: [f64; 2], filter: Filter) -> Option<Item> {
        if !(at[0].is_finite() && at[1].is_finite()) {
            return None;
        }
        let (o, d) = proj.ray(at[0], at[1]);
        let screen = |p: Vec3| proj.screen(p).map(|s| ((s[0] - at[0]).powi(2) + (s[1] - at[1]).powi(2)).sqrt());
        match filter {
            Filter::Faces => self.geometry_face(self.surface(o, d)?).map(Item::Geometry),
            Filter::Bodies => self.body(self.surface(o, d)?).map(Item::Geometry),
            Filter::ElementFaces => match self.surface(o, d)? {
                Surface::Element(i, _) => {
                    let t = self.display.tris.get(i)?;
                    Some(Item::Face(ElementFace { element: t.element, face: t.seen_face(d) }))
                }
                Surface::Geometry(..) => None,
            },
            Filter::Elements => {
                let lines = self.display.lines.iter().enumerate().filter_map(|(i, l)| {
                    let (a, b) = (proj.screen(l.a)?, proj.screen(l.b)?);
                    let (distance, k) = segment_distance(at, [a[0], a[1]], [b[0], b[1]]);
                    Some((i, distance, l.a + (l.b - l.a) * k))
                });
                if let Some(l) = self.nearest_visible(proj, lines.collect(), EDGE_RADIUS).and_then(|i| self.display.lines.get(i)) {
                    return Some(Item::Element(l.element));
                }
                match self.surface(o, d)? {
                    Surface::Element(i, _) => self.display.tris.get(i).map(|t| Item::Element(t.element)),
                    Surface::Geometry(..) => None,
                }
            }
            Filter::Nodes => {
                let near = self.display.nodes.iter().enumerate().filter_map(|(i, n)| Some((i, screen(n.p)?, n.p))).collect();
                self.nearest_visible(proj, near, POINT_RADIUS).and_then(|i| self.display.nodes.get(i)).map(|n| Item::Node(n.label))
            }
            Filter::Vertices => {
                // A design's B-rep vertices first, by their persistent names.
                if self.curves {
                    let design: Vec<(&BodyNames, &(String, Vec3))> = self.names.iter().flat_map(|b| b.vertices.iter().map(move |v| (b, v))).collect();
                    let near = design.iter().enumerate().filter_map(|(i, (_, (_, p)))| Some((i, screen(*p)?, *p))).collect();
                    if let Some((b, (name, _))) = nearest_passing(near, POINT_RADIUS, |p| self.curve_visible(proj, p)).and_then(|i| design.get(i)) {
                        return b.vertex(name).map(Item::Geometry);
                    }
                }
                let tagged = &self.index.vertices;
                let near = tagged.iter().enumerate().filter_map(|(i, (n, _))| {
                    let p = self.display.node_position(*n)?;
                    Some((i, screen(p)?, p))
                });
                let (_, tag) = tagged.get(self.nearest_visible(proj, near.collect(), POINT_RADIUS)?)?;
                self.entity(*tag).cloned().map(Item::Geometry)
            }
            Filter::Edges => {
                // A design's B-rep edges first: they are what an edge is named by before meshing,
                // and what seeds and sets refer to after.
                if let Some(g) = self.pick_curve(proj, at) {
                    return Some(Item::Geometry(g));
                }
                // Element edges and line elements on geometry edges.
                let mut tagged: Vec<(Vec3, Vec3, u32)> = self
                    .display
                    .edges
                    .iter()
                    .filter_map(|e| Some((e.p[0], e.p[1], *self.index.edge_tag.get(&edge_key(e.corners[0], e.corners[1]))?)))
                    .collect();
                tagged.extend(self.display.lines.iter().zip(&self.index.line_tag).filter_map(|(l, t)| Some((l.a, l.b, (*t)?))));
                let near = tagged.iter().enumerate().filter_map(|(i, (a, b, _))| {
                    let (sa, sb) = (proj.screen(*a)?, proj.screen(*b)?);
                    // The edge's point nearest the pointer on screen is the one checked.
                    let (distance, k) = segment_distance(at, [sa[0], sa[1]], [sb[0], sb[1]]);
                    Some((i, distance, *a + (*b - *a) * k))
                });
                let (_, _, tag) = tagged.get(self.nearest_visible(proj, near.collect(), EDGE_RADIUS)?)?;
                self.entity(*tag).cloned().map(Item::Geometry)
            }
        }
    }

    /// The visible nodes, element faces or elements (with `filter`) inside the screen box
    /// `rect` ([x0, y0, x1, y1], any two opposite corners), sorted; other filters select
    /// nothing. ([`Selection::apply`] keeps at most [`MAX_SELECTION`].)
    pub fn pick_box(&self, proj: &Projector, rect: [f64; 4], filter: Filter) -> Vec<Item> {
        if !filter.boxes() || !rect.iter().all(|v| v.is_finite()) {
            return Vec::new();
        }
        let (x0, x1, y0, y1) = (rect[0].min(rect[2]), rect[0].max(rect[2]), rect[1].min(rect[3]), rect[1].max(rect[3]));
        let inside = |p: Vec3| proj.screen(p).is_some_and(|s| (x0..=x1).contains(&s[0]) && (y0..=y1).contains(&s[1]));
        let Some(depth) = DepthMap::new(proj, self.drawn()) else { return Vec::new() };
        let mut out = BTreeSet::new();
        let mut add = |item: Item| {
            out.insert(item);
        };
        match filter {
            Filter::Nodes => {
                for n in self.display.nodes.iter().filter(|n| inside(n.p) && depth.visible(n.p)) {
                    add(Item::Node(n.label));
                }
            }
            Filter::ElementFaces | Filter::Elements => {
                let faces = filter == Filter::ElementFaces;
                for t in self.display.tris.iter().filter(|_| self.elements) {
                    let centre = (t.p[0] + t.p[1] + t.p[2]) * (1.0 / 3.0);
                    if t.p.iter().all(|p| inside(*p)) && depth.visible(centre) {
                        add(if faces {
                            let dir = proj.screen(centre).map_or(Vec3::ZERO, |s| proj.ray(s[0], s[1]).1);
                            Item::Face(ElementFace { element: t.element, face: t.seen_face(dir) })
                        } else {
                            Item::Element(t.element)
                        });
                    }
                }
                if !faces {
                    for l in &self.display.lines {
                        if inside(l.a) && inside(l.b) && depth.visible((l.a + l.b) * 0.5) {
                            add(Item::Element(l.element));
                        }
                    }
                }
            }
            _ => {}
        }
        out.into_iter().collect()
    }

    /// The selection and the hovered item drawn over the scene: their faces in `selected` (or
    /// `hover`) colour, their edges as wide lines and their nodes as large points.
    pub fn overlay(&self, selection: &Selection, hover: Option<&Item>, selected: Rgb, hover_color: Rgb) -> Scene {
        let mut s = Scene::default();
        self.mark(&Marks::of(self, selection.items().iter()), selected, &mut s);
        if let Some(h) = hover {
            self.mark(&Marks::of(self, std::iter::once(h)), hover_color, &mut s);
        }
        s
    }

    /// `items` (a set's members) drawn over the scene in `color`.
    pub fn items_overlay(&self, items: &[Item], color: Rgb) -> Scene {
        let mut s = Scene::default();
        self.mark(&Marks::of(self, items.iter()), color, &mut s);
        s
    }

    /// The hovered item drawn over the scene in `color`. With `fill` off (while result contours
    /// are shown) its faces are left out and only their outline, edges and nodes are drawn, so
    /// the hover never hides the colours under the pointer.
    pub fn hover(&self, item: &Item, color: Rgb, fill: bool) -> Scene {
        let mut s = Scene::default();
        self.mark(&Marks::of(self, std::iter::once(item)), color, &mut s);
        if !fill {
            s.tris.clear();
        }
        s
    }

    fn mark(&self, m: &Marks, color: Rgb, s: &mut Scene) {
        if m.is_empty() {
            return;
        }
        let first = s.tris.len();
        let flag = |flags: &[bool], tag: Option<u32>| tag.is_some_and(|t| flags.get(t as usize).copied().unwrap_or(false));
        let on = |tag: Option<u32>| flag(&m.tags, tag);
        let in_body = |tag: Option<u32>| flag(&m.in_bodies, tag);
        if self.elements {
            for (i, t) in self.display.tris.iter().enumerate() {
                let (face, region) = (self.index.face_tag.get(i).copied().flatten(), self.index.region_tag.get(i).copied().flatten());
                // A shell's triangle draws both of its sides.
                let shell_back = t.kind.dimension() == 2 && m.faces.contains(&ElementFace { element: t.element, face: 2 });
                if m.elements.contains(&t.element)
                    || m.faces.contains(&ElementFace { element: t.element, face: t.face })
                    || shell_back
                    || on(face)
                    || in_body(region)
                    || in_body(face)
                {
                    s.tris.push(SceneTri { p: t.p, normal: t.normal, color });
                }
            }
            for (l, tag) in self.display.lines.iter().zip(&self.index.line_tag) {
                if m.elements.contains(&l.element) || on(*tag) {
                    s.lines.push(SceneLine { a: l.a, b: l.b, color, width: EDGE_WIDTH });
                }
            }
        }
        if self.geometry && !m.refs.is_empty() {
            // Which faces of each drawn body are marked, worked out once per face.
            let marked: Vec<Vec<bool>> = self
                .names
                .iter()
                .map(|body| {
                    (0..body.faces.len())
                        .map(|f| {
                            let face = body.face(u32::try_from(f).unwrap_or(u32::MAX));
                            m.refs.iter().any(|r| match &r.entity {
                                Entity::Body => r.body == body.name,
                                _ => face.as_ref().is_some_and(|g| r.selects(g)),
                            })
                        })
                        .collect()
                })
                .collect();
            let whole = |b: usize| self.names.get(b).is_some_and(|n| m.bodies.contains(&n.name));
            for t in &self.display.surface_tris {
                let b = t.body as usize;
                if whole(b) || marked.get(b).and_then(|f| f.get(t.face as usize)).copied().unwrap_or(false) {
                    s.tris.push(SceneTri { p: t.p, normal: t.normal, color });
                }
            }
        }
        outline(s, first, color);
        // A design edge is marked along its own polyline, which the scene draws over the mesh too
        // (most meshes do not tag their element edges with the design's edge names), and along
        // the mesh's edges a mesh does tag with it (below).
        if m.edges && self.curves {
            // Which edges of each drawn body are marked, worked out once per edge.
            let marked: Vec<Vec<bool>> = self
                .names
                .iter()
                .map(|body| {
                    (0..body.edges.len())
                        .map(|e| body.edge(u32::try_from(e).unwrap_or(u32::MAX)).is_some_and(|g| m.refs.iter().any(|r| r.selects(&g))))
                        .collect()
                })
                .collect();
            for c in &self.display.curves {
                if marked.get(c.body as usize).and_then(|b| b.get(c.edge as usize)).copied().unwrap_or(false) {
                    s.lines.push(SceneLine { a: c.a, b: c.b, color, width: EDGE_WIDTH });
                }
            }
        }
        if m.edges {
            for e in &self.display.edges {
                if on(self.index.edge_tag.get(&edge_key(e.corners[0], e.corners[1])).copied()) {
                    s.lines.push(SceneLine { a: e.p[0], b: e.p[1], color, width: EDGE_WIDTH });
                }
            }
        }
        if m.vertices && self.curves {
            for body in self.names {
                for (name, p) in &body.vertices {
                    if body.vertex(name).is_some_and(|g| m.refs.iter().any(|r| r.selects(&g))) {
                        s.points.push(ScenePoint { p: *p, color, size: NODE_SIZE });
                    }
                }
            }
        }
        let vertices = self.index.vertices.iter().filter(|_| m.vertices).filter(|(_, t)| on(Some(*t))).map(|(n, _)| *n);
        for n in m.nodes.iter().copied().chain(vertices) {
            if let Some(p) = self.display.node_position(n) {
                s.points.push(ScenePoint { p, color, size: NODE_SIZE });
            }
        }
    }
}

/// The outline of the triangles `s.tris[first..]`: their sides no other of them shares, drawn as
/// lines so a highlight reads at a glance.
fn outline(s: &mut Scene, first: usize, color: Rgb) {
    let key = |p: Vec3| [p.x, p.y, p.z].map(f64::to_bits);
    let mut sides: Vec<([[u64; 3]; 2], Vec3, Vec3)> = Vec::new();
    for t in s.tris.get(first..).unwrap_or(&[]) {
        for (a, b) in [(t.p[0], t.p[1]), (t.p[1], t.p[2]), (t.p[2], t.p[0])] {
            let (ka, kb) = (key(a), key(b));
            sides.push((if ka < kb { [ka, kb] } else { [kb, ka] }, a, b));
        }
    }
    sides.sort_by_key(|x| x.0);
    for run in sides.chunk_by(|x, y| x.0 == y.0) {
        if let [(_, a, b)] = run {
            s.lines.push(SceneLine { a: *a, b: *b, color, width: OUTLINE_WIDTH });
        }
    }
}

/// What a group of items marks, for drawing.
#[derive(Default)]
struct Marks {
    /// Geometry tags selected (by index).
    tags: Vec<bool>,
    /// Geometry tags of an entity of a body selected whole.
    in_bodies: Vec<bool>,
    /// Bodies selected whole.
    bodies: BTreeSet<String>,
    /// The geometry references (for the drawn geometry surfaces).
    refs: Vec<GeometryRef>,
    /// Whether an edge or a vertex is among the references.
    edges: bool,
    vertices: bool,
    nodes: BTreeSet<u32>,
    elements: BTreeSet<u32>,
    faces: BTreeSet<ElementFace>,
}

impl Marks {
    fn of<'a>(p: &Picker, items: impl Iterator<Item = &'a Item>) -> Marks {
        let c = p.model.classification();
        let n = c.entities().len();
        let mut m = Marks { tags: vec![false; n], in_bodies: vec![false; n], ..Marks::default() };
        for item in items {
            match item {
                Item::Geometry(g) => {
                    match g.entity {
                        Entity::Body => {
                            m.bodies.insert(g.body.clone());
                        }
                        Entity::Edge(_) => m.edges = true,
                        Entity::Vertex(_) => m.vertices = true,
                        Entity::Face(_) => {}
                    }
                    for t in c.tags_selected_by(g) {
                        if let Some(slot) = m.tags.get_mut(t as usize) {
                            *slot = true;
                        }
                    }
                    m.refs.push(g.clone());
                }
                Item::Node(n) => {
                    m.nodes.insert(*n);
                }
                Item::Element(e) => {
                    m.elements.insert(*e);
                }
                Item::Face(f) => {
                    m.faces.insert(*f);
                }
            }
        }
        if !m.bodies.is_empty() {
            for (slot, e) in m.in_bodies.iter_mut().zip(c.entities()) {
                *slot = m.bodies.contains(&e.body);
            }
        }
        m
    }

    fn is_empty(&self) -> bool {
        self.refs.is_empty() && self.nodes.is_empty() && self.elements.is_empty() && self.faces.is_empty()
    }
}

/// Where a command form puts the selection.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Target {
    /// The members themselves (a set's definition).
    Members,
    /// A set made from the selection, named in `target: {"set"}` (boundary conditions, forces).
    TargetSet,
    /// A set made from the selection, named in `set` (sections, couplings, surface and body
    /// loads).
    Set,
    /// The selected geometry edges, as `edges: [{body, edge}]` (edge seeds).
    Edges,
    /// A point at what is selected (`geometry.centre`: a vertex, a node, a hole's centre, a
    /// centroid), as the reference point's `xyz`.
    Position,
    /// A datum point on what is selected: at a vertex, along an edge, at a face's centre, or at
    /// the point `geometry.centre` finds (a hole's centre, nodes, several items).
    DatumPoint,
    /// The selected bodies (or the bodies of the selected faces, edges and vertices) as the
    /// `bodies` list (the cells of a partitioned mesh).
    Bodies,
    /// The selected faces as an extrusion's `source` while it has none, then as its
    /// `to_surface`: pick the source face, then the target surface.
    ExtrudeFaces,
    /// The selected geometry faces as `faces: [{body, face}]`, in place of a sphere or labels
    /// (the faces Relax Nodes relaxes, the faces a local size refines).
    Faces,
}

/// Where command `id`'s form takes the selection, if it does. `params` (the form's current
/// parameters) tells a force (target) from a pressure (set).
pub fn target(id: &str, params: &Value) -> Option<Target> {
    match id {
        "model.set.create" | "model.set.edit" => Some(Target::Members),
        "mesh.seed.edge" | "mesh.seed.clear" => Some(Target::Edges),
        "model.reference_point.create" | "model.reference_point.edit" => Some(Target::Position),
        "geometry.datum.point" | "geometry.datum.point.edit" => Some(Target::DatumPoint),
        "mesh.cells" => Some(Target::Bodies),
        "mesh.hex.extrude" => Some(Target::ExtrudeFaces),
        "mesh.node.relax" | "mesh.size.set" => Some(Target::Faces),
        "model.boundary_condition.create" | "model.boundary_condition.edit" => Some(Target::TargetSet),
        "model.section.create" | "model.section.edit" | "model.coupling.create" | "model.coupling.edit" | "interaction.rigid_wall" => {
            Some(Target::Set)
        }
        "model.load.create" | "model.load.edit" => {
            let kind = params.get("type").and_then(Value::as_str);
            Some(if matches!(kind, None | Some("force" | "moment")) { Target::TargetSet } else { Target::Set })
        }
        _ => None,
    }
}

/// The filter command `id`'s form picks with, when its target is one kind of entity: edges for
/// edge seeds, and the vertex, edge or face a datum point is defined on.
pub fn form_filter(id: &str, params: &Value) -> Option<Filter> {
    match target(id, params) {
        Some(Target::Edges) => return Some(Filter::Edges),
        Some(Target::ExtrudeFaces) => return Some(Filter::Faces),
        _ => {}
    }
    match id {
        "geometry.datum.point" | "geometry.datum.point.edit" => match params.pointer("/by/type").and_then(Value::as_str) {
            Some("vertex") => Some(Filter::Vertices),
            Some("edge") => Some(Filter::Edges),
            Some("face_centre") => Some(Filter::Faces),
            _ => None,
        },
        _ => None,
    }
}

/// `params` with the selection's members put in (replacing any members it had).
pub fn with_members(params: &Value, selection: &Selection) -> Value {
    let mut out = params.as_object().cloned().unwrap_or_default();
    for key in ["geometry", "nodes", "elements", "faces", "edges"] {
        out.remove(key);
    }
    if let Value::Object(members) = selection.params() {
        out.extend(members);
    }
    Value::Object(out)
}

/// `params` with the selection's geometry edges as `edges` (replacing any it had); `None` when no
/// edge is selected.
pub fn with_edges(params: &Value, selection: &Selection) -> Option<Value> {
    let edges: Vec<Value> = selection
        .items()
        .iter()
        .filter_map(|i| match i {
            Item::Geometry(g) if matches!(g.entity, Entity::Edge(_)) => Some(json!(g)),
            _ => None,
        })
        .collect();
    if edges.is_empty() {
        return None;
    }
    let mut out = params.as_object().cloned().unwrap_or_default();
    out.insert("edges".into(), Value::Array(edges));
    Some(Value::Object(out))
}

/// `params` with the selection's geometry faces as `faces`, and no sphere, nodes or elements
/// (what Relax Nodes relaxes); `None` when no face is selected.
pub fn with_faces(params: &Value, selection: &Selection) -> Option<Value> {
    let faces: Vec<Value> = selection
        .items()
        .iter()
        .filter_map(|i| match i {
            Item::Geometry(g) if matches!(g.entity, Entity::Face(_)) => Some(json!(g)),
            _ => None,
        })
        .collect();
    if faces.is_empty() {
        return None;
    }
    let mut out = params.as_object().cloned().unwrap_or_default();
    for key in ["center", "radius", "around_element", "nodes", "elements"] {
        out.remove(key);
    }
    out.insert("faces".into(), Value::Array(faces));
    Some(Value::Object(out))
}

/// The key of an extrusion's form the selected faces go into: `source` while it names no face
/// (and no `source_set`), then `to_surface`.
pub fn extrude_key(params: &Value) -> &'static str {
    let named = |k: &str| {
        params
            .get(k)
            .is_some_and(|v| v.as_array().is_some_and(|a| a.iter().any(|f| f.get("face").and_then(Value::as_str).is_some_and(|n| !n.is_empty()))))
    };
    if named("source") || params.get("source_set").and_then(Value::as_str).is_some_and(|s| !s.is_empty()) { "to_surface" } else { "source" }
}

/// `params` with the selection's faces as an extrusion's source or target surface
/// ([`extrude_key`]); `None` when no face is selected.
pub fn with_extrude_faces(params: &Value, selection: &Selection) -> Option<Value> {
    let faces: Vec<Value> = selection
        .items()
        .iter()
        .filter_map(|i| match i {
            Item::Geometry(g) if matches!(g.entity, Entity::Face(_)) => Some(json!(g)),
            _ => None,
        })
        .collect();
    if faces.is_empty() {
        return None;
    }
    let mut out = params.as_object().cloned().unwrap_or_default();
    out.insert(extrude_key(params).into(), Value::Array(faces));
    Some(Value::Object(out))
}

/// `params` with the set `set` as the target.
pub fn with_set(params: &Value, target: Target, set: &str) -> Value {
    let mut out = params.as_object().cloned().unwrap_or_default();
    match target {
        Target::TargetSet => {
            out.insert("target".into(), json!({"set": set}));
        }
        Target::Set => {
            out.insert("set".into(), json!(set));
        }
        Target::Members | Target::Edges | Target::Position | Target::DatumPoint | Target::Bodies | Target::ExtrudeFaces | Target::Faces => {}
    }
    Value::Object(out)
}

/// The set of `model` whose members are exactly those of `set` (`model.set.create`'s
/// parameters), if any.
pub fn same_set(model: &Model, set: &Value) -> Option<String> {
    let mut members = set.as_object().cloned()?;
    members.remove("name");
    let wanted = serde_json::from_value::<Set>(Value::Object(members)).ok()?.normalized().ok()?;
    model.sets().iter().find(|(_, s)| **s == wanted).map(|(n, _)| n.clone())
}

/// The bodies of the selection's geometry (a body, or the body of a face, edge or vertex), each
/// once; `None` when no geometry is selected.
pub fn bodies(selection: &Selection) -> Option<Vec<String>> {
    let mut out: Vec<String> = Vec::new();
    for item in selection.items() {
        if let Item::Geometry(g) = item
            && !out.contains(&g.body)
        {
            out.push(g.body.clone());
        }
    }
    (!out.is_empty()).then_some(out)
}

/// A datum point's `by` for the selection, given the point `geometry.centre` found for it
/// (`centre`: its reply): one vertex, one edge or one face by name, unless the centre is a
/// hole's (`from: circle`); anything else at the centre's coordinates.
pub fn datum_point_by(selection: &Selection, centre: &Value) -> Value {
    let circle = centre.get("from").and_then(Value::as_str) == Some("circle");
    let point = centre.get("point").cloned().unwrap_or(Value::Null);
    let items = selection.items();
    match items.iter().next() {
        Some(Item::Geometry(g)) if items.len() == 1 && !circle => match &g.entity {
            Entity::Vertex(v) => json!({"type": "vertex", "body": g.body, "vertex": v}),
            Entity::Edge(e) => json!({"type": "edge", "body": g.body, "edge": e, "fraction": 0.5}),
            Entity::Face(f) => json!({"type": "face_centre", "body": g.body, "face": f}),
            Entity::Body => json!({"type": "coordinates", "point": point}),
        },
        _ => json!({"type": "coordinates", "point": point}),
    }
}

/// A set name not taken in `model`: `<base>`, then `<base>-2`, `<base>-3`, …; `Set` when
/// `base` is not a valid name.
pub fn fresh_set_name(model: &Model, base: &str) -> String {
    let base = if check_name("set", base).is_ok() { base } else { "Set" };
    std::iter::once(base.to_string())
        .chain((2..100_000).map(|k| format!("{base}-{k}")))
        .find(|n| check_name("set", n).is_ok() && model.set(n).is_none())
        .unwrap_or_else(|| "Set".into())
}

#[cfg(test)]
mod tests {
    use fieldcraft_model::ElementKind;
    use fieldcraft_render::{Camera, StandardView};

    use super::*;

    /// The filter a selection holds: one kind names it, a mix or nothing names none (#341).
    #[test]
    fn a_selection_names_the_filter_of_what_it_holds() {
        let edge = Item::Geometry(GeometryRef::new("Bar", Entity::Edge("top".into())).unwrap());
        let mut sel = Selection::default();
        assert_eq!(sel.filter(), None);
        sel.apply([edge], Mode::Replace);
        assert_eq!(sel.filter(), Some(Filter::Edges));
        sel.apply([Item::Node(3)], Mode::Add);
        assert_eq!(sel.filter(), None);
        sel.apply([Item::Node(4)], Mode::Replace);
        assert_eq!(sel.filter(), Some(Filter::Nodes));
    }

    /// `nx × ny × nz` unit hex8 cells, the faces at x = nx tagged as face `end` of body `Bar`, the
    /// solids as the body, the edge along x = nx, z = nz as edge `top` and the corner node at
    /// (nx, 0, 0) as vertex `v`.
    fn bar(nx: u32, ny: u32, nz: u32) -> Model {
        let mut m = Model::default();
        let id = |i: u32, j: u32, k: u32| 1 + i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
        for k in 0..=nz {
            for j in 0..=ny {
                for i in 0..=nx {
                    m.set_node(id(i, j, k), [i as f64, j as f64, k as f64]).unwrap();
                }
            }
        }
        let body = m.geometry_tag(&GeometryRef::body("Bar").unwrap()).unwrap();
        let end = m.geometry_tag(&GeometryRef::face("Bar", "end").unwrap()).unwrap();
        let top = m.geometry_tag(&GeometryRef::new("Bar", Entity::Edge("top".into())).unwrap()).unwrap();
        let v = m.geometry_tag(&GeometryRef::new("Bar", Entity::Vertex("v".into())).unwrap()).unwrap();
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
                    m.tag_element(label, body).unwrap();
                    if i + 1 == nx {
                        // Hex8 face 4: corners 1, 2, 6, 5 (x = i + 1).
                        m.tag_element_face(label, 4, end).unwrap();
                        if k + 1 == nz {
                            // Hex8 edge 6 joins corners 5 and 6 (x = i + 1, z = k + 1).
                            m.tag_element_edge(label, 6, top).unwrap();
                        }
                    }
                    label += 1;
                }
            }
        }
        m.tag_node(id(nx, 0, 0), v).unwrap();
        m
    }

    struct Scene3 {
        model: Model,
        display: DisplayMesh,
        index: Index,
    }

    impl Scene3 {
        fn new(model: Model) -> Scene3 {
            let display = DisplayMesh::from_model(&model);
            let index = Index::new(&model, &display);
            Scene3 { model, display, index }
        }
        fn picker(&self) -> Picker<'_> {
            Picker { display: &self.display, names: &[], model: &self.model, index: &self.index, elements: true, geometry: false, curves: false }
        }
    }

    /// A front-on view of the bar's end (looking along −x from +x), 400 × 300 pixels.
    fn end_view(d: &DisplayMesh) -> Projector {
        let mut cam = Camera::default();
        cam.set_view(StandardView::Right);
        cam.fit(&d.bounds);
        Projector::new(&cam, 400.0, 300.0, d.bounds.diagonal())
    }

    fn pixel(p: &Projector, at: Vec3) -> [f64; 2] {
        let s = p.screen(at).unwrap();
        [s[0], s[1]]
    }

    #[test]
    fn a_click_on_the_end_picks_its_face_body_element_face_element_and_node() {
        let s = Scene3::new(bar(4, 2, 2));
        let p = s.picker();
        let proj = end_view(&s.display);
        // The centre of the end cell at y ∈ [0, 1], z ∈ [1, 2]: element 1 + 3 + 0·4 + 1·8 = 12.
        let at = pixel(&proj, Vec3::new(4.0, 0.5, 1.5));
        assert_eq!(p.pick(&proj, at, Filter::Faces), Some(Item::Geometry(GeometryRef::face("Bar", "end").unwrap())));
        assert_eq!(p.pick(&proj, at, Filter::Bodies), Some(Item::Geometry(GeometryRef::body("Bar").unwrap())));
        assert_eq!(p.pick(&proj, at, Filter::ElementFaces), Some(Item::Face(ElementFace { element: 12, face: 4 })));
        assert_eq!(p.pick(&proj, at, Filter::Elements), Some(Item::Element(12)));
        // Near the corner (4, 1, 2), its node; far from every node, nothing.
        let corner = pixel(&proj, Vec3::new(4.0, 1.0, 2.0));
        let n = p.pick(&proj, [corner[0] + 3.0, corner[1] - 2.0], Filter::Nodes);
        assert_eq!(n, Some(Item::Node(1 + 4 + 5 + 2 * 15)));
        assert_eq!(p.pick(&proj, at, Filter::Nodes), None);
        // The tagged vertex and edge.
        let v = pixel(&proj, Vec3::new(4.0, 0.0, 0.0));
        assert_eq!(p.pick(&proj, v, Filter::Vertices), Some(Item::Geometry(GeometryRef::new("Bar", Entity::Vertex("v".into())).unwrap())));
        let e = pixel(&proj, Vec3::new(4.0, 1.5, 2.0));
        assert_eq!(
            p.pick(&proj, [e[0], e[1] + 2.0], Filter::Edges),
            Some(Item::Geometry(GeometryRef::new("Bar", Entity::Edge("top".into())).unwrap()))
        );
        // Off the model, and hostile pixels.
        assert_eq!(p.pick(&proj, [2.0, 2.0], Filter::Faces), None);
        assert_eq!(p.pick(&proj, [f64::NAN, 1.0], Filter::Nodes), None);
        assert_eq!(p.pick(&proj, [1e300, -1e300], Filter::Edges), None);
    }

    /// Nodes behind the end face are not picked through it, even right under the pointer.
    #[test]
    fn hidden_nodes_are_not_picked() {
        let s = Scene3::new(bar(4, 2, 2));
        let p = s.picker();
        let mut cam = Camera::default();
        cam.set_view(StandardView::Right);
        cam.fit(&s.display.bounds);
        // Orthographic along −x: the node (3, 1, 1) projects onto (4, 1, 1), which is in front.
        let proj = Projector::new(&cam, 400.0, 300.0, s.display.bounds.diagonal());
        let at = pixel(&proj, Vec3::new(4.0, 1.0, 1.0));
        let front = 1 + 4 + 5 + 15;
        assert_eq!(p.pick(&proj, at, Filter::Nodes), Some(Item::Node(front)));
        let items = p.pick_box(&proj, [0.0, 0.0, 400.0, 300.0], Filter::Nodes);
        assert_eq!(items.len(), 3 * 3, "the end's nodes only: {items:?}");
        assert!(items.contains(&Item::Node(front)));
        // Element faces and elements in the box: the end's eight cells.
        let faces = p.pick_box(&proj, [0.0, 0.0, 400.0, 300.0], Filter::ElementFaces);
        assert_eq!(faces.len(), 4);
        assert!(faces.iter().all(|f| matches!(f, Item::Face(ElementFace { face: 4, .. }))));
        let elements = p.pick_box(&proj, [400.0, 300.0, 0.0, 0.0], Filter::Elements);
        assert_eq!(elements, vec![Item::Element(4), Item::Element(8), Item::Element(12), Item::Element(16)]);
        // Filters without boxes, empty and hostile boxes.
        assert_eq!(p.pick_box(&proj, [0.0, 0.0, 400.0, 300.0], Filter::Faces), vec![]);
        assert_eq!(p.pick_box(&proj, [5.0, 5.0, 6.0, 6.0], Filter::Nodes), vec![]);
        assert_eq!(p.pick_box(&proj, [f64::NAN, 0.0, 1.0, 1.0], Filter::Nodes), vec![]);
    }

    #[test]
    fn an_empty_scene_picks_nothing() {
        let s = Scene3::new(Model::default());
        let p = s.picker();
        let proj = Projector::new(&Camera::default(), 400.0, 300.0, 1.0);
        for f in Filter::ALL {
            assert_eq!(p.pick(&proj, [200.0, 150.0], f), None);
            assert_eq!(p.pick_box(&proj, [0.0, 0.0, 400.0, 300.0], f), vec![]);
        }
        assert!(p.overlay(&Selection::default(), None, Rgb(1, 2, 3), Rgb(4, 5, 6)).tris.is_empty());
    }

    /// Geometry surfaces (no mesh) are named by their bodies' face names.
    #[test]
    fn drawn_geometry_names_its_faces() {
        let mut display = DisplayMesh::default();
        let square = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]];
        // Its one B-rep edge: the diagonal between the two faces.
        display.add_curves(&[vec![square[0], square[2]]]);
        display.add_surface(&square, &[[0, 1, 2], [0, 2, 3]], &[0, 1]);
        let names =
            [BodyNames { name: "Plate".into(), faces: vec!["F1:a".into(), "F1:b".into()], edges: vec!["F1:a|F1:b".into()], vertices: Vec::new() }];
        let model = Model::default();
        let index = Index::new(&model, &display);
        let p = Picker { display: &display, names: &names, model: &model, index: &index, elements: false, geometry: true, curves: true };
        let mut cam = Camera::default();
        cam.set_view(StandardView::Top);
        cam.fit(&display.bounds);
        let proj = Projector::new(&cam, 400.0, 300.0, 2.0);
        let at = pixel(&proj, Vec3::new(0.8, 0.2, 0.0));
        let face = GeometryRef::face("Plate", "F1:a").unwrap();
        assert_eq!(p.pick(&proj, at, Filter::Faces), Some(Item::Geometry(face.clone())));
        assert_eq!(p.pick(&proj, at, Filter::Bodies), Some(Item::Geometry(GeometryRef::body("Plate").unwrap())));
        assert_eq!(p.pick(&proj, at, Filter::ElementFaces), None);
        // The face is highlighted: its one triangle in the selection colour.
        let mut sel = Selection::default();
        sel.apply([Item::Geometry(face)], Mode::Replace);
        let o = p.overlay(&sel, None, Rgb(9, 9, 9), Rgb(0, 0, 0));
        assert_eq!(o.tris.len(), 1);
        assert_eq!(o.tris[0].color, Rgb(9, 9, 9));
        // The design edge is picked by its persistent name, and drawn wide when selected.
        let edge = GeometryRef::new("Plate", Entity::Edge("F1:a|F1:b".into())).unwrap();
        let on = pixel(&proj, Vec3::new(0.5, 0.5, 0.0));
        assert_eq!(p.pick(&proj, on, Filter::Edges), Some(Item::Geometry(edge.clone())));
        assert_eq!(p.pick(&proj, at, Filter::Edges), None);
        sel.apply([Item::Geometry(edge)], Mode::Replace);
        let o = p.overlay(&sel, None, Rgb(9, 9, 9), Rgb(0, 0, 0));
        assert_eq!((o.tris.len(), o.lines.len()), (0, 1));
        // Without the curves drawn, the edge is not picked.
        let hidden = Picker { curves: false, ..p };
        assert_eq!(hidden.pick(&proj, on, Filter::Edges), None);
    }

    #[test]
    fn the_overlay_marks_what_is_selected_and_hovered() {
        let s = Scene3::new(bar(4, 2, 2));
        let p = s.picker();
        let (sel_c, hov_c) = (Rgb(200, 0, 0), Rgb(0, 200, 0));
        let mut sel = Selection::default();
        sel.apply([Item::Geometry(GeometryRef::face("Bar", "end").unwrap())], Mode::Replace);
        let o = p.overlay(&sel, None, sel_c, hov_c);
        // Four quads at the end, two triangles each.
        assert_eq!(o.tris.len(), 8);
        assert!(o.tris.iter().all(|t| t.p.iter().all(|q| q.x == 4.0)));
        // Hovering an element adds its skin triangles in the hover colour, drawn after.
        let o = p.overlay(&sel, Some(&Item::Element(1)), sel_c, hov_c);
        assert_eq!(o.tris.iter().filter(|t| t.color == hov_c).count(), 6);
        // Over result contours the hover is its outline only: no face fills.
        let filled = p.hover(&Item::Element(1), hov_c, true);
        assert_eq!(filled.tris.len(), 6);
        let outlined = p.hover(&Item::Element(1), hov_c, false);
        assert!(outlined.tris.is_empty());
        assert_eq!(outlined.lines, filled.lines);
        assert!(!outlined.lines.is_empty());
        // Nodes, the edge and the vertex.
        sel.apply(
            [
                Item::Node(1),
                Item::Geometry(GeometryRef::new("Bar", Entity::Edge("top".into())).unwrap()),
                Item::Geometry(GeometryRef::new("Bar", Entity::Vertex("v".into())).unwrap()),
            ],
            Mode::Replace,
        );
        let o = p.overlay(&sel, None, sel_c, hov_c);
        assert_eq!((o.tris.len(), o.lines.len(), o.points.len()), (0, 2, 2));
        // A whole body: every skin triangle.
        sel.apply([Item::Geometry(GeometryRef::body("Bar").unwrap())], Mode::Replace);
        assert_eq!(p.overlay(&sel, None, sel_c, hov_c).tris.len(), s.display.tris.len());
    }

    /// An extrusion's form takes the source face first, then the target surface (#400).
    #[test]
    fn an_extrusion_takes_its_source_then_its_target_surface() {
        assert_eq!(target("mesh.hex.extrude", &json!({})), Some(Target::ExtrudeFaces));
        assert_eq!(form_filter("mesh.hex.extrude", &json!({})), Some(Filter::Faces));
        let mut s = Selection::default();
        s.apply([Item::Geometry(GeometryRef::face("Post", "F1:-z").unwrap()), Item::Node(4)], Mode::Replace);
        let form = with_extrude_faces(&json!({"layers": 4, "source": [{"body": "", "face": ""}]}), &s).unwrap();
        assert_eq!(form, json!({"layers": 4, "source": [{"body": "Post", "face": "F1:-z"}]}));
        s.apply([Item::Geometry(GeometryRef::face("Wall", "F7").unwrap())], Mode::Replace);
        let form = with_extrude_faces(&form, &s).unwrap();
        assert_eq!(form["to_surface"], json!([{"body": "Wall", "face": "F7"}]));
        assert_eq!(form["source"], json!([{"body": "Post", "face": "F1:-z"}]), "the source stays");
        // Picking again replaces the target surface; a source set counts as a source.
        assert_eq!(extrude_key(&form), "to_surface");
        assert_eq!(extrude_key(&json!({"source_set": "Section"})), "to_surface");
        s.apply([Item::Node(4)], Mode::Replace);
        assert!(with_extrude_faces(&form, &s).is_none(), "no face selected");
    }

    #[test]
    fn selections_add_remove_toggle_and_cap() {
        let mut s = Selection::default();
        s.apply([Item::Node(3), Item::Node(1), Item::Node(3)], Mode::Replace);
        assert_eq!(s.len(), 2);
        s.apply([Item::Element(7)], Mode::Add);
        s.apply([Item::Node(1), Item::Node(9)], Mode::Toggle);
        assert_eq!(s.items().iter().cloned().collect::<Vec<_>>(), vec![Item::Node(3), Item::Node(9), Item::Element(7)]);
        s.apply([Item::Node(3)], Mode::Remove);
        assert_eq!(s.summary(), "1 element, 1 node");
        s.apply((1..=MAX_SELECTION as u32 + 10).map(Item::Node), Mode::Replace);
        assert!(s.len() == MAX_SELECTION && s.capped());
        s.clear();
        assert!(s.is_empty() && !s.capped() && s.summary() == "nothing selected");
    }

    /// The selection becomes a set's members, and reads back from them.
    #[test]
    fn the_selection_maps_to_set_members_and_command_parameters() {
        let mut s = Selection::default();
        let face = GeometryRef::face("Body1", "F2:top").unwrap();
        s.apply([Item::Geometry(face.clone()), Item::Node(5), Item::Element(2), Item::Face(ElementFace { element: 2, face: 3 })], Mode::Replace);
        let p = s.params();
        assert_eq!(
            p,
            json!({"geometry": [{"body": "Body1", "face": "F2:top"}], "nodes": [5], "elements": [2], "faces": [{"element": 2, "face": 3}]})
        );
        assert_eq!(items_from_params(&p).unwrap().into_iter().collect::<BTreeSet<_>>(), *s.items());
        assert_eq!(s.summary(), "1 element, 1 element face, 1 face, 1 node");
        // A set's form: the members replace the old ones, the name stays.
        let form = with_members(&json!({"name": "end", "nodes": [1, 2], "edges": [{"element": 1, "edge": 1}]}), &s);
        assert_eq!(form["name"], "end");
        assert_eq!(form["nodes"], json!([5]));
        assert!(form.get("edges").is_none());
        // Definitions that take a set.
        assert_eq!(target("model.set.create", &json!({})), Some(Target::Members));
        assert_eq!(target("model.boundary_condition.create", &json!({})), Some(Target::TargetSet));
        assert_eq!(target("model.load.create", &json!({"type": "force"})), Some(Target::TargetSet));
        assert_eq!(target("model.load.create", &json!({"type": "pressure"})), Some(Target::Set));
        assert_eq!(target("model.section.edit", &json!({})), Some(Target::Set));
        assert_eq!(target("model.coupling.create", &json!({})), Some(Target::Set));
        assert_eq!(target("mesh.node.add", &json!({})), None);
        assert_eq!(form_filter("mesh.seed.edge", &json!({})), Some(Filter::Edges));
        assert_eq!(form_filter("geometry.datum.point", &json!({"by": {"type": "vertex"}})), Some(Filter::Vertices));
        assert_eq!(form_filter("model.load.create", &json!({})), None);
        assert_eq!(with_set(&json!({"name": "fix"}), Target::TargetSet, "fix-set"), json!({"name": "fix", "target": {"set": "fix-set"}}));
        assert_eq!(with_set(&json!({"name": "web"}), Target::Set, "web-set"), json!({"name": "web", "set": "web-set"}));
        let mut m = Model::default();
        assert_eq!(fresh_set_name(&m, "fix-set"), "fix-set");
        m.add_set("fix-set", fieldcraft_model::Set { nodes: vec![1], ..Default::default() }).unwrap();
        assert_eq!(fresh_set_name(&m, "fix-set"), "fix-set-2");
        assert_eq!(fresh_set_name(&m, " \u{7}"), "Set");
        assert_eq!(fresh_set_name(&m, "my clamp-set"), "Set");
    }

    #[test]
    fn hostile_members_are_refused() {
        for bad in [
            json!({"nodes": [0]}),
            json!({"nodes": "x"}),
            json!({"nodes": [-1]}),
            json!({"elements": [1.5]}),
            json!({"faces": [{"element": 1, "face": 0}]}),
            json!({"faces": [{"element": 0, "face": 1}]}),
            json!({"geometry": [{"body": ""}]}),
            json!({"geometry": [{"body": "B", "face": "f", "edge": "e"}]}),
            json!({"nodes": vec![1; MAX_SELECTION + 1]}),
        ] {
            assert!(items_from_params(&bad).is_err(), "{bad}");
        }
        assert_eq!(items_from_params(&json!({})).unwrap(), vec![]);
        assert_eq!(items_from_params(&json!(null)).unwrap(), vec![]);
    }
}
