//! Degrees of freedom: what a node carries and how the global unknowns are numbered.
//!
//! Every unknown of a discrete system is a ([`DofKind`], node) pair. An element declares the
//! kinds each of its nodes uses ([`ElementDofs`]): three translations for solids, six
//! (translations and rotations) for shells and beams, two for plane elements, three (Ux, Uy, Rz)
//! for the beams of a 2D model, a temperature for
//! heat transfer, a pressure for mixed u–p elements, and the internal unknowns it condenses out
//! itself (incompatible modes, a discontinuous pressure), which never reach the global system.
//! A node carries the union of what its elements, constraints and boundary conditions need.
//!
//! [`DofLayout`] numbers them: mesh nodes in label order, each node's kinds in [`DofKind`] order,
//! then the mesh's edges and faces that carry unknowns (edge elements,
//! docs/rfc/0001-electromagnetics.md §3.1), each by its nodes' sorted labels, then extra points
//! (reference points) in the order they are added, so the same model always gives the same
//! numbering (the convention of T. J. R. Hughes, *The Finite Element Method*,
//! Dover 2000, §2.8–2.10, with the ID array replaced by per-node offsets). Assembly reads a
//! node's unknowns as one contiguous block, so a node pair is a dense block of the global matrix.
//!
//! Boundary conditions and loads address a node and a [`DofKind`]; the model's 1-based degree
//! of freedom numbers (1–3 translations, 4–6 rotations, 11 temperature, 8 pressure, the common
//! convention of the deck formats) convert with [`DofKind::from_number`].

use std::collections::BTreeMap;
use std::ops::Range;

/// One kind of nodal unknown. The order is the numbering order within a node.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DofKind {
    /// Translation along X.
    Ux,
    /// Translation along Y.
    Uy,
    /// Translation along Z.
    Uz,
    /// Rotation about X.
    Rx,
    /// Rotation about Y.
    Ry,
    /// Rotation about Z.
    Rz,
    /// Temperature.
    T,
    /// Pressure (hydrostatic, of a mixed u–p element).
    P,
    /// The first field unknown of an edge or a face (an edge's Whitney function, a face's first
    /// function: `fieldcraft_fem::nedelec`).
    E1,
    /// The second (an edge's gradient function, a face's second function).
    E2,
}

impl DofKind {
    pub const ALL: [DofKind; 10] =
        [DofKind::Ux, DofKind::Uy, DofKind::Uz, DofKind::Rx, DofKind::Ry, DofKind::Rz, DofKind::T, DofKind::P, DofKind::E1, DofKind::E2];

    /// Position in [`DofKind::ALL`].
    pub const fn index(self) -> usize {
        self as usize
    }

    /// The translation along axis `axis` (0–2).
    pub fn translation(axis: usize) -> Option<DofKind> {
        [DofKind::Ux, DofKind::Uy, DofKind::Uz].get(axis).copied()
    }

    /// The rotation about axis `axis` (0–2).
    pub fn rotation(axis: usize) -> Option<DofKind> {
        [DofKind::Rx, DofKind::Ry, DofKind::Rz].get(axis).copied()
    }

    /// The structural direction 0–5 (translations, then rotations), `None` for T and P.
    pub fn direction(self) -> Option<usize> {
        (self.index() < 6).then_some(self.index())
    }

    /// The kind of a structural direction 0–5.
    pub fn from_direction(direction: usize) -> Option<DofKind> {
        DofKind::ALL.get(direction).copied().filter(|k| k.index() < 6)
    }

    /// The model's and decks' 1-based number: 1–3 translations, 4–6 rotations, 8 pressure,
    /// 11 temperature; 0 for the unknowns of edges and faces, which have none.
    pub fn number(self) -> u8 {
        match self {
            DofKind::T => 11,
            DofKind::P => 8,
            DofKind::E1 | DofKind::E2 => 0,
            k => k.index() as u8 + 1,
        }
    }

    /// The kind of a model degree of freedom number ([`DofKind::number`]).
    pub fn from_number(number: u8) -> Option<DofKind> {
        DofKind::ALL.into_iter().find(|k| number != 0 && k.number() == number)
    }

    /// What it is, in the user's terms.
    pub fn describe(self) -> &'static str {
        match self {
            DofKind::Ux => "X",
            DofKind::Uy => "Y",
            DofKind::Uz => "Z",
            DofKind::Rx => "rotation about X",
            DofKind::Ry => "rotation about Y",
            DofKind::Rz => "rotation about Z",
            DofKind::T => "temperature",
            DofKind::P => "pressure",
            DofKind::E1 => "its first field unknown",
            DofKind::E2 => "its second field unknown",
        }
    }

    const fn bit(self) -> u16 {
        1 << self as u16
    }
}

/// A set of [`DofKind`]s, iterated in numbering order.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct DofSet(u16);

impl DofSet {
    pub const EMPTY: DofSet = DofSet(0);
    /// Ux, Uy, Uz: solids, trusses, membranes in 3D.
    pub const TRANSLATIONS: DofSet = DofSet(0b0000_0111);
    /// Ux, Uy: plane elements.
    pub const PLANE: DofSet = DofSet(0b0000_0011);
    /// Ux, Uy, Rz: beams in a 2D model (the nodes of a plane frame).
    pub const PLANE_FRAME: DofSet = DofSet(0b0010_0011);
    /// Rx, Ry, Rz.
    pub const ROTATIONS: DofSet = DofSet(0b0011_1000);
    /// Translations and rotations: shells, beams, reference points.
    pub const STRUCTURAL: DofSet = DofSet(0b0011_1111);
    pub const TEMPERATURE: DofSet = DofSet(DofKind::T.bit());
    pub const PRESSURE: DofSet = DofSet(DofKind::P.bit());
    /// E1: an edge of a first-order edge element.
    pub const FIELD: DofSet = DofSet(DofKind::E1.bit());
    /// E1, E2: the edges and faces of a second-order one.
    pub const FIELDS: DofSet = DofSet(DofKind::E1.bit() | DofKind::E2.bit());

    pub const fn of(kinds: &[DofKind]) -> DofSet {
        let mut bits = 0;
        let mut i = 0;
        while i < kinds.len() {
            bits |= kinds[i].bit();
            i += 1;
        }
        DofSet(bits)
    }

    pub const fn contains(self, kind: DofKind) -> bool {
        self.0 & kind.bit() != 0
    }

    /// Whether every kind of `other` is in this set.
    pub const fn includes(self, other: DofSet) -> bool {
        self.0 & other.0 == other.0
    }

    pub const fn union(self, other: DofSet) -> DofSet {
        DofSet(self.0 | other.0)
    }

    pub const fn len(self) -> usize {
        self.0.count_ones() as usize
    }

    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }

    /// The position of `kind` among this set's kinds (its offset in a node's block).
    pub const fn position(self, kind: DofKind) -> Option<usize> {
        if self.contains(kind) { Some((self.0 & (kind.bit() - 1)).count_ones() as usize) } else { None }
    }

    /// The kinds, in numbering order.
    pub fn iter(self) -> impl Iterator<Item = DofKind> {
        DofKind::ALL.into_iter().filter(move |k| self.contains(*k))
    }
}

impl FromIterator<DofKind> for DofSet {
    fn from_iter<I: IntoIterator<Item = DofKind>>(iter: I) -> DofSet {
        DofSet(iter.into_iter().fold(0, |m, k| m | k.bit()))
    }
}

/// What an element type contributes to the unknowns: the kinds every one of its nodes carries
/// (its element matrices are ordered node by node, these kinds in order within a node) and the
/// internal unknowns it condenses out before assembly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ElementDofs {
    pub node: DofSet,
    pub internal: Internal,
}

/// Element-internal unknowns, condensed at the element level (static condensation; Hughes
/// 2000, §4.7 for incompatible modes, §4.4 for discontinuous pressure).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Internal {
    /// Incompatible (bubble) displacement modes.
    pub incompatible_modes: usize,
    /// Discontinuous pressure unknowns of a mixed element.
    pub pressure: usize,
}

impl Internal {
    pub const NONE: Internal = Internal { incompatible_modes: 0, pressure: 0 };

    pub const fn count(self) -> usize {
        self.incompatible_modes + self.pressure
    }
}

impl ElementDofs {
    /// Three translations per node, nothing internal.
    pub const SOLID: ElementDofs = ElementDofs { node: DofSet::TRANSLATIONS, internal: Internal::NONE };

    pub const fn new(node: DofSet) -> ElementDofs {
        ElementDofs { node, internal: Internal::NONE }
    }

    /// Unknowns of an element with `nodes` nodes, before condensation.
    pub const fn total(self, nodes: usize) -> usize {
        nodes * self.node.len() + self.internal.count()
    }
}

/// What a slot of a layout numbers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Entity {
    /// A mesh node, by label.
    Node(u32),
    /// A mesh edge, by its nodes' labels, ascending.
    Edge([u32; 2]),
    /// A mesh face, by its nodes' labels, ascending.
    Face([u32; 3]),
    /// An extra point (a reference point), by label.
    Point(u32),
}

/// One numbered entity: what it is, its kinds and its first global unknown.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Slot {
    entity: Entity,
    kinds: DofSet,
    first: usize,
}

/// The global numbering of the unknowns.
///
/// Slots (positions) `0..nodes()` are the mesh nodes in label order; then the edges and the
/// faces that carry unknowns, each in the order of their sorted labels; then the extra points in
/// the order they were added. A slot's unknowns are contiguous, in [`DofKind`] order.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct DofLayout {
    slots: Vec<Slot>,
    nodes: usize,
    edges: usize,
    faces: usize,
    total: usize,
}

/// Collects what every node, edge and face needs before numbering: the union of the kinds asked
/// for.
#[derive(Debug, Clone, Default)]
pub struct LayoutBuilder {
    nodes: BTreeMap<u32, DofSet>,
    edges: BTreeMap<[u32; 2], DofSet>,
    faces: BTreeMap<[u32; 3], DofSet>,
    points: Vec<(u32, DofSet)>,
}

impl LayoutBuilder {
    /// Node `label` needs `kinds` (added to what it needs already).
    pub fn node(&mut self, label: u32, kinds: DofSet) -> &mut Self {
        let entry = self.nodes.entry(label).or_default();
        *entry = entry.union(kinds);
        self
    }

    /// Every node of an element of this kind needs the element's kinds.
    pub fn element(&mut self, labels: &[u32], dofs: ElementDofs) -> &mut Self {
        for &l in labels {
            self.node(l, dofs.node);
        }
        self
    }

    /// The edge between nodes `a` and `b` (either order) needs `kinds`; no kinds declare
    /// nothing (an edge or a face is a slot only when it carries unknowns), and neither does a
    /// node twice, which is no edge (it has no slot, and [`DofLayout::edge`] finds none).
    pub fn edge(&mut self, [a, b]: [u32; 2], kinds: DofSet) -> &mut Self {
        if kinds.is_empty() || a == b {
            return self;
        }
        let entry = self.edges.entry([a.min(b), a.max(b)]).or_default();
        *entry = entry.union(kinds);
        self
    }

    /// The face of nodes `nodes` (any order) needs `kinds`; no kinds declare nothing, and
    /// neither does a node repeated, which is no face.
    pub fn face(&mut self, mut nodes: [u32; 3], kinds: DofSet) -> &mut Self {
        nodes.sort_unstable();
        let [a, b, c] = nodes;
        if kinds.is_empty() || a == b || b == c {
            return self;
        }
        let entry = self.faces.entry(nodes).or_default();
        *entry = entry.union(kinds);
        self
    }

    /// An extra point (a reference point) numbered after every mesh entity, in the order added;
    /// each call adds one.
    pub fn point(&mut self, label: u32, kinds: DofSet) -> &mut Self {
        self.points.push((label, kinds));
        self
    }

    pub fn build(&self) -> DofLayout {
        let entities = (self.nodes.iter().map(|(l, k)| (Entity::Node(*l), *k)))
            .chain(self.edges.iter().map(|(e, k)| (Entity::Edge(*e), *k)))
            .chain(self.faces.iter().map(|(f, k)| (Entity::Face(*f), *k)))
            .chain(self.points.iter().map(|(l, k)| (Entity::Point(*l), *k)));
        let mut total = 0;
        let slots = entities
            .map(|(entity, kinds)| {
                let slot = Slot { entity, kinds, first: total };
                total += kinds.len();
                slot
            })
            .collect();
        DofLayout { slots, nodes: self.nodes.len(), edges: self.edges.len(), faces: self.faces.len(), total }
    }
}

impl DofLayout {
    pub fn builder() -> LayoutBuilder {
        LayoutBuilder::default()
    }

    /// Every label in `labels` (sorted or not) with the same kinds.
    pub fn uniform(labels: impl IntoIterator<Item = u32>, kinds: DofSet) -> DofLayout {
        let mut b = LayoutBuilder::default();
        for l in labels {
            b.node(l, kinds);
        }
        b.build()
    }

    /// Total unknowns.
    pub fn dofs(&self) -> usize {
        self.total
    }

    /// Mesh nodes (slots `0..nodes()`).
    pub fn nodes(&self) -> usize {
        self.nodes
    }

    /// Mesh edges with unknowns (slots after the nodes).
    pub fn edges(&self) -> usize {
        self.edges
    }

    /// Mesh faces with unknowns (slots after the edges).
    pub fn faces(&self) -> usize {
        self.faces
    }

    /// Mesh entities and extra points.
    pub fn slots(&self) -> usize {
        self.slots.len()
    }

    /// The first slot after the mesh's entities: the first extra point's.
    fn mesh_slots(&self) -> usize {
        self.nodes + self.edges + self.faces
    }

    /// Unknowns of the mesh's entities (nodes, edges, faces); the extra points' follow them.
    pub fn node_dofs(&self) -> usize {
        self.slots.get(self.mesh_slots()).map_or(self.total, |s| s.first)
    }

    /// The slot of extra point `k`.
    pub fn point_slot(&self, k: usize) -> usize {
        self.mesh_slots() + k
    }

    /// The slot of a mesh node label.
    pub fn node(&self, label: u32) -> Option<usize> {
        self.slots.get(..self.nodes)?.binary_search_by(|s| s.entity.cmp(&Entity::Node(label))).ok()
    }

    /// The slot of the edge between nodes `a` and `b` (either order), if it carries unknowns.
    pub fn edge(&self, [a, b]: [u32; 2]) -> Option<usize> {
        let edges = self.slots.get(self.nodes..self.nodes + self.edges)?;
        edges.binary_search_by(|s| s.entity.cmp(&Entity::Edge([a.min(b), a.max(b)]))).ok().map(|k| self.nodes + k)
    }

    /// The slot of the face of nodes `nodes` (any order), if it carries unknowns.
    pub fn face(&self, mut nodes: [u32; 3]) -> Option<usize> {
        nodes.sort_unstable();
        let start = self.nodes + self.edges;
        let faces = self.slots.get(start..start + self.faces)?;
        faces.binary_search_by(|s| s.entity.cmp(&Entity::Face(nodes))).ok().map(|k| start + k)
    }

    /// The slot of a label: a mesh node, else an extra point.
    pub fn slot(&self, label: u32) -> Option<usize> {
        let points = self.mesh_slots();
        self.node(label).or_else(|| self.slots.get(points..)?.iter().position(|s| s.entity == Entity::Point(label)).map(|k| points + k))
    }

    /// What a slot numbers.
    pub fn entity(&self, slot: usize) -> Option<Entity> {
        self.slots.get(slot).map(|s| s.entity)
    }

    /// The label of a node's or an extra point's slot (`None` for an edge or a face).
    pub fn label(&self, slot: usize) -> Option<u32> {
        match self.entity(slot)? {
            Entity::Node(l) | Entity::Point(l) => Some(l),
            Entity::Edge(_) | Entity::Face(_) => None,
        }
    }

    /// The kinds a slot carries.
    pub fn kinds(&self, slot: usize) -> DofSet {
        self.slots.get(slot).map_or(DofSet::EMPTY, |s| s.kinds)
    }

    /// A slot's unknowns (empty for a slot that does not exist).
    pub fn range(&self, slot: usize) -> Range<usize> {
        self.slots.get(slot).map_or(0..0, |s| s.first..s.first + s.kinds.len())
    }

    /// The global unknown of `kind` at `slot`, if the slot carries it.
    pub fn dof(&self, slot: usize, kind: DofKind) -> Option<usize> {
        let s = self.slots.get(slot)?;
        s.kinds.position(kind).map(|p| s.first + p)
    }

    /// The slot and kind of a global unknown.
    pub fn locate(&self, dof: usize) -> Option<(usize, DofKind)> {
        if dof >= self.total {
            return None;
        }
        let slot = self.slots.partition_point(|s| s.first <= dof).checked_sub(1)?;
        let s = self.slots.get(slot)?;
        let kind = s.kinds.iter().nth(dof - s.first)?;
        Some((slot, kind))
    }

    /// The global unknowns of an element: its slots' `kinds`, node by node, `None` when a slot
    /// lacks one of them.
    pub fn element_dofs(&self, slots: &[usize], kinds: DofSet) -> Option<Vec<usize>> {
        let mut out = Vec::with_capacity(slots.len() * kinds.len());
        for &slot in slots {
            for k in kinds.iter() {
                out.push(self.dof(slot, k)?);
            }
        }
        Some(out)
    }

    /// `Some(s)` when every slot carries the same `s` kinds, so unknown `s · slot + offset` is
    /// the `offset`-th kind of `slot` (node-major vectors of `s`-vectors).
    pub fn stride(&self) -> Option<usize> {
        let first = self.slots.first()?.kinds;
        self.slots.iter().all(|s| s.kinds == first).then_some(first.len())
    }

    /// What a global unknown is, in the user's terms.
    pub fn describe(&self, dof: usize) -> String {
        match self.locate(dof) {
            Some((slot, kind)) => {
                let what = match self.entity(slot) {
                    Some(Entity::Node(l)) => format!("node {l}"),
                    Some(Entity::Point(l)) => format!("point {l}"),
                    Some(Entity::Edge([a, b])) => format!("edge {a}–{b}"),
                    Some(Entity::Face([a, b, c])) => format!("face {a}–{b}–{c}"),
                    None => String::new(),
                };
                format!("{what} in {}", kind.describe())
            }
            None => format!("unknown degree of freedom {dof}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kinds_numbers_and_sets() {
        for k in DofKind::ALL.into_iter().filter(|k| k.number() != 0) {
            assert_eq!(DofKind::from_number(k.number()), Some(k));
        }
        assert_eq!(DofKind::from_number(1), Some(DofKind::Ux));
        assert_eq!(DofKind::from_number(6), Some(DofKind::Rz));
        assert_eq!(DofKind::from_number(11), Some(DofKind::T));
        assert_eq!(DofKind::from_number(0), None);
        assert_eq!(DofKind::from_number(7), None);
        assert_eq!(DofKind::from_direction(4), Some(DofKind::Ry));
        assert_eq!(DofKind::from_direction(6), None);
        let s = DofSet::of(&[DofKind::T, DofKind::Ux, DofKind::Uy, DofKind::Uz]);
        assert_eq!(s.len(), 4);
        assert_eq!(s.iter().collect::<Vec<_>>(), vec![DofKind::Ux, DofKind::Uy, DofKind::Uz, DofKind::T]);
        assert_eq!(s.position(DofKind::T), Some(3));
        assert_eq!(s.position(DofKind::Rx), None);
        assert!(DofSet::STRUCTURAL.includes(DofSet::TRANSLATIONS));
        assert!(!DofSet::TRANSLATIONS.includes(DofSet::STRUCTURAL));
        assert_eq!(DofSet::TRANSLATIONS.union(DofSet::ROTATIONS), DofSet::STRUCTURAL);
        assert_eq!(DofSet::PLANE_FRAME.iter().collect::<Vec<_>>(), vec![DofKind::Ux, DofKind::Uy, DofKind::Rz]);
        // The edge and face unknowns have no model number, and number 0 is none.
        assert_eq!((DofKind::E1.number(), DofKind::E2.number()), (0, 0));
        assert_eq!(DofSet::FIELDS.iter().collect::<Vec<_>>(), vec![DofKind::E1, DofKind::E2]);
        assert_eq!(DofSet::STRUCTURAL.union(DofSet::FIELDS).position(DofKind::E2), Some(7));
        assert_eq!(ElementDofs::SOLID.total(8), 24);
        let hybrid = ElementDofs { node: DofSet::TRANSLATIONS, internal: Internal { incompatible_modes: 0, pressure: 1 } };
        assert_eq!(hybrid.total(8), 25);
    }

    #[test]
    fn layout_numbers_nodes_by_label_then_kind_then_points() {
        let mut b = DofLayout::builder();
        // A shell on 7 and 3, a solid on 3, 5 and 7, a heat element on 5; a reference point.
        b.element(&[7, 3], ElementDofs::new(DofSet::STRUCTURAL));
        b.element(&[3, 5, 7], ElementDofs::SOLID);
        b.node(5, DofSet::TEMPERATURE);
        b.point(100, DofSet::STRUCTURAL);
        let l = b.build();
        assert_eq!((l.nodes(), l.slots(), l.dofs(), l.node_dofs()), (3, 4, 22, 16));
        assert_eq!(l.node(3), Some(0));
        assert_eq!(l.node(5), Some(1));
        assert_eq!(l.node(100), None);
        assert_eq!(l.slot(100), Some(3));
        assert_eq!(l.range(1), 6..10);
        assert_eq!(l.dof(1, DofKind::T), Some(9));
        assert_eq!(l.dof(1, DofKind::Rx), None);
        assert_eq!(l.dof(2, DofKind::Rz), Some(15));
        assert_eq!(l.dof(3, DofKind::Ux), Some(16));
        assert_eq!(l.locate(9), Some((1, DofKind::T)));
        assert_eq!(l.locate(16), Some((3, DofKind::Ux)));
        assert_eq!(l.locate(22), None);
        for d in 0..l.dofs() {
            let (slot, kind) = l.locate(d).unwrap();
            assert_eq!(l.dof(slot, kind), Some(d));
        }
        assert_eq!(l.element_dofs(&[0, 1], DofSet::TRANSLATIONS), Some(vec![0, 1, 2, 6, 7, 8]));
        assert_eq!(l.element_dofs(&[0, 1], DofSet::STRUCTURAL), None);
        assert_eq!(l.describe(9), "node 5 in temperature");
        assert_eq!(l.describe(19), "point 100 in rotation about X");
        // The same needs in another order give the same layout.
        let mut c = DofLayout::builder();
        c.point(100, DofSet::STRUCTURAL).node(5, DofSet::TEMPERATURE).element(&[5, 7, 3], ElementDofs::SOLID);
        c.element(&[3, 7], ElementDofs::new(DofSet::STRUCTURAL));
        assert_eq!(c.build(), l);
        // A uniform translational layout is the classic 3 × node numbering.
        let u = DofLayout::uniform([4, 2, 9], DofSet::TRANSLATIONS);
        assert_eq!(u.dof(2, DofKind::Uy), Some(7));
        assert_eq!(u.dofs(), 9);
        assert_eq!(u.range(7), 0..0);
        assert_eq!(u.locate(usize::MAX), None);
        assert_eq!(u.stride(), Some(3));
        assert_eq!(l.stride(), None);
        assert_eq!(DofLayout::default().stride(), None);
    }

    #[test]
    fn edges_and_faces_are_numbered_after_the_nodes_by_their_sorted_labels() {
        // Two tetrahedra of nodes 1, 4, 6, 9 and 4, 6, 9, 12 with second-order edge unknowns, a
        // node with a temperature and a reference point.
        let mut b = DofLayout::builder();
        b.point(500, DofSet::STRUCTURAL).node(6, DofSet::TEMPERATURE);
        for tet in [[9, 1, 6, 4], [12, 6, 4, 9]] {
            for (i, a) in tet.iter().enumerate() {
                for c in tet.iter().skip(i + 1) {
                    b.edge([*a, *c], DofSet::FIELDS);
                }
            }
            for skip in 0..4 {
                let f: Vec<u32> = tet.iter().enumerate().filter(|(k, _)| *k != skip).map(|(_, l)| *l).collect();
                b.face([f[0], f[1], f[2]], DofSet::FIELDS);
            }
        }
        let l = b.build();
        // One node, 9 edges (the shared face's 3 once), 7 faces, the point.
        assert_eq!((l.nodes(), l.edges(), l.faces(), l.slots()), (1, 9, 7, 18));
        assert_eq!(l.edge([4, 1]), Some(1));
        assert_eq!(l.edge([1, 4]), Some(1));
        assert_eq!(l.entity(1), Some(Entity::Edge([1, 4])));
        assert_eq!(l.edge([1, 12]), None);
        assert_eq!(l.face([9, 4, 6]), l.face([4, 6, 9]));
        assert_eq!(l.entity(l.face([9, 6, 4]).unwrap()), Some(Entity::Face([4, 6, 9])));
        assert_eq!(l.face([1, 4, 12]), None);
        // Each slot's unknowns follow the previous one's; the node's temperature comes first.
        assert_eq!(l.range(1), 1..3);
        assert_eq!(l.dof(1, DofKind::E2), Some(2));
        assert_eq!((l.node_dofs(), l.dofs()), (1 + 2 * 16, 1 + 2 * 16 + 6));
        assert_eq!(l.slot(500), Some(17));
        assert_eq!(l.point_slot(0), 17);
        assert_eq!((l.label(1), l.label(0), l.label(17)), (None, Some(6), Some(500)));
        assert_eq!(l.describe(2), "edge 1–4 in its second field unknown");
        let last_face = l.face([6, 9, 12]).unwrap();
        assert_eq!(l.describe(l.range(last_face).start), "face 6–9–12 in its first field unknown");
        for d in 0..l.dofs() {
            let (slot, kind) = l.locate(d).unwrap();
            assert_eq!(l.dof(slot, kind), Some(d));
        }
        // The same needs in another order give the same layout.
        let mut c = DofLayout::builder();
        for f in [[12, 9, 6], [1, 6, 9], [4, 12, 6], [9, 6, 4], [4, 1, 9], [4, 9, 12], [1, 4, 6]] {
            c.face(f, DofSet::FIELD).face(f, DofSet::FIELDS);
        }
        for e in [[12, 9], [6, 1], [9, 4], [12, 6], [4, 6], [9, 1], [4, 12], [1, 4], [6, 9]] {
            c.edge(e, DofSet::FIELDS);
        }
        c.node(6, DofSet::TEMPERATURE).point(500, DofSet::STRUCTURAL);
        // An edge or a face with no unknowns is no slot, and a repeated node is no edge or face.
        c.edge([1, 12], DofSet::EMPTY).face([1, 6, 12], DofSet::EMPTY);
        c.edge([4, 4], DofSet::FIELDS).face([6, 1, 6], DofSet::FIELDS).face([9, 9, 9], DofSet::FIELDS);
        assert_eq!(c.build(), l);
        assert_eq!((l.edge([4, 4]), l.face([1, 6, 6])), (None, None));
    }
}
