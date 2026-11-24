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
//! then extra points (reference points) in the order they are added, so the same model always
//! gives the same numbering (the convention of T. J. R. Hughes, *The Finite Element Method*,
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
}

impl DofKind {
    pub const ALL: [DofKind; 8] = [DofKind::Ux, DofKind::Uy, DofKind::Uz, DofKind::Rx, DofKind::Ry, DofKind::Rz, DofKind::T, DofKind::P];

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
    /// 11 temperature.
    pub fn number(self) -> u8 {
        match self {
            DofKind::T => 11,
            DofKind::P => 8,
            k => k.index() as u8 + 1,
        }
    }

    /// The kind of a model degree of freedom number ([`DofKind::number`]).
    pub fn from_number(number: u8) -> Option<DofKind> {
        DofKind::ALL.into_iter().find(|k| k.number() == number)
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
        }
    }

    const fn bit(self) -> u8 {
        1 << self as u8
    }
}

/// A set of [`DofKind`]s, iterated in numbering order.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct DofSet(u8);

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

/// One numbered node: its label, its kinds and its first global unknown.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Slot {
    label: u32,
    kinds: DofSet,
    first: usize,
}

/// The global numbering of the unknowns.
///
/// Slots (positions) `0..nodes()` are the mesh nodes in label order; `nodes()..slots()` are the
/// extra points in the order they were added. A slot's unknowns are contiguous, in
/// [`DofKind`] order.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct DofLayout {
    slots: Vec<Slot>,
    nodes: usize,
    total: usize,
}

/// Collects what every node needs before numbering: the union of the kinds asked for.
#[derive(Debug, Clone, Default)]
pub struct LayoutBuilder {
    nodes: BTreeMap<u32, DofSet>,
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

    /// An extra point (a reference point) numbered after every mesh node, in the order added;
    /// each call adds one.
    pub fn point(&mut self, label: u32, kinds: DofSet) -> &mut Self {
        self.points.push((label, kinds));
        self
    }

    pub fn build(&self) -> DofLayout {
        let mut slots = Vec::with_capacity(self.nodes.len() + self.points.len());
        let mut total = 0;
        for (&label, &kinds) in self.nodes.iter().chain(self.points.iter().map(|p| (&p.0, &p.1))) {
            slots.push(Slot { label, kinds, first: total });
            total += kinds.len();
        }
        DofLayout { slots, nodes: self.nodes.len(), total }
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

    /// Mesh nodes and extra points.
    pub fn slots(&self) -> usize {
        self.slots.len()
    }

    /// Unknowns of the mesh nodes; the extra points' follow them.
    pub fn node_dofs(&self) -> usize {
        self.slots.get(self.nodes).map_or(self.total, |s| s.first)
    }

    /// The slot of extra point `k`.
    pub fn point_slot(&self, k: usize) -> usize {
        self.nodes + k
    }

    /// The slot of a mesh node label.
    pub fn node(&self, label: u32) -> Option<usize> {
        self.slots.get(..self.nodes)?.binary_search_by_key(&label, |s| s.label).ok()
    }

    /// The slot of a label: a mesh node, else an extra point.
    pub fn slot(&self, label: u32) -> Option<usize> {
        self.node(label).or_else(|| self.slots.get(self.nodes..)?.iter().position(|s| s.label == label).map(|k| self.nodes + k))
    }

    pub fn label(&self, slot: usize) -> Option<u32> {
        self.slots.get(slot).map(|s| s.label)
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
                let what = if slot < self.nodes { "node" } else { "point" };
                format!("{what} {} in {}", self.label(slot).unwrap_or(0), kind.describe())
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
        for k in DofKind::ALL {
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
}
