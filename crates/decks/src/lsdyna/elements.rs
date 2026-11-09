//! Element connectivity in LS-DYNA's node order, and the element formulation (ELFORM) for each
//! element kind and section formulation.
//!
//! Node order (from the element figures of the public keyword manual):
//! - hexahedra: nodes 1–4 on one face, 5–8 over them, as in the model; the 20-node hexahedron's
//!   mid-side nodes are on edges 1-2, 2-3, 3-4, 4-1, then 5-6 … 8-5, then 1-5 … 4-8, also as in
//!   the model;
//! - tetrahedra are degenerate hexahedra `n1 n2 n3 n4 n4 n4 n4 n4`;
//! - wedges are degenerate hexahedra `n1 n2 n3 n4 n5 n5 n6 n6` whose nodes 1–4 are a
//!   quadrilateral face, 5 is over edge 1-2 and 6 over edge 3-4. The model's wedge (triangles
//!   a1-a2-a3 and a4-a5-a6) becomes `a2 a1 a4 a5 a3 a3 a6 a6`: the side a2-a1-a4-a5, then the
//!   ridge a3–a6;
//! - the 10-node tetrahedron has its mid-side nodes 5, 6, 7 on edges 1-2, 2-3, 3-1, and 8, 9, 10
//!   on edges 2-4, 3-4, 1-4 (the model puts 8, 9, 10 on 1-4, 2-4, 3-4);
//! - shells: counter-clockwise corners (the top surface), triangles as `n1 n2 n3 n3` and six-node
//!   triangles as `n1 n2 n3 n3 n4 n5 n6 n6`; eight-node quadrilaterals as in the model.
//!
//! Every mapping is checked by the sign of the element Jacobian in the tests.

use fieldcraft_model::{ElementKind, Formulation, Integration, SectionKind};

/// How an element is written: its keyword and its node list in LS-DYNA order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Keyword {
    /// `*ELEMENT_SOLID`: up to 10 nodes on the second card.
    Solid,
    /// `*ELEMENT_SOLID_H20`: 20 nodes on two cards.
    SolidH20,
    /// `*ELEMENT_SHELL`: 4 or 8 nodes.
    Shell,
    /// `*ELEMENT_BEAM`: 2 nodes (trusses).
    Beam,
    /// `*ELEMENT_BEAM_ORIENTATION`: 2 nodes, then the section's orientation vector on a second
    /// card (beams).
    BeamOrientation,
}

impl Keyword {
    pub(crate) fn of(kind: ElementKind) -> Keyword {
        match kind {
            ElementKind::Hex20 => Keyword::SolidH20,
            ElementKind::Tet4 | ElementKind::Tet10 | ElementKind::Wedge6 | ElementKind::Pyramid5 | ElementKind::Hex8 => Keyword::Solid,
            ElementKind::Tri3 | ElementKind::Tri6 | ElementKind::Quad4 | ElementKind::Quad8 => Keyword::Shell,
            ElementKind::Rod2 | ElementKind::Rod3 => Keyword::Beam,
            ElementKind::Beam2 => Keyword::BeamOrientation,
        }
    }

    pub(crate) fn name(self) -> &'static str {
        match self {
            Keyword::Solid => "*ELEMENT_SOLID",
            Keyword::SolidH20 => "*ELEMENT_SOLID_H20",
            Keyword::Shell => "*ELEMENT_SHELL",
            Keyword::Beam => "*ELEMENT_BEAM",
            Keyword::BeamOrientation => "*ELEMENT_BEAM_ORIENTATION",
        }
    }
}

/// Positions in the model's node list, in LS-DYNA's order (repeats for degenerate forms).
pub(crate) fn order(kind: ElementKind) -> &'static [usize] {
    match kind {
        ElementKind::Rod2 | ElementKind::Beam2 => &[0, 1],
        ElementKind::Rod3 => &[0, 1, 2],
        ElementKind::Tri3 => &[0, 1, 2, 2],
        ElementKind::Tri6 => &[0, 1, 2, 2, 3, 4, 5, 5],
        ElementKind::Quad4 => &[0, 1, 2, 3],
        ElementKind::Quad8 => &[0, 1, 2, 3, 4, 5, 6, 7],
        ElementKind::Tet4 => &[0, 1, 2, 3, 3, 3, 3, 3],
        ElementKind::Tet10 => &[0, 1, 2, 3, 4, 5, 6, 8, 9, 7],
        ElementKind::Wedge6 => &[1, 0, 3, 4, 2, 2, 5, 5],
        // The degenerate hexahedron: the base, then the apex four times.
        ElementKind::Pyramid5 => &[0, 1, 2, 3, 4, 4, 4, 4],
        ElementKind::Hex8 => &[0, 1, 2, 3, 4, 5, 6, 7],
        ElementKind::Hex20 => &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
    }
}

/// The element's node labels in LS-DYNA's order.
pub(crate) fn connectivity(kind: ElementKind, nodes: &[u32]) -> Option<Vec<u32>> {
    order(kind).iter().map(|&i| nodes.get(i).copied()).collect()
}

/// Hourglass control a reduced-integration formulation needs (`*HOURGLASS`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Hourglass {
    /// Solids: the assumed-strain stiffness form (type 6) with coefficient 1, which gives
    /// accurate coarse-mesh bending for elastic materials.
    Solid,
    /// Quadrilateral shells: the stiffness form (type 4) with the default coefficient.
    Shell,
}

/// A part's element formulation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Elform {
    pub(crate) elform: i32,
    pub(crate) hourglass: Option<Hourglass>,
}

/// The formulation of a section's elements of one kind.
///
/// Solids: hex8 2 (eight-point), reduced 1 (one-point) with hourglass control, incompatible
/// modes −18 (eight-point enhanced strain with incompatible modes); wedges 2 (as degenerate
/// hexahedra), reduced 15 (two-point pentahedron); pyramids 2 (as degenerate hexahedra); tet4 10 (one point, exact for a linear
/// tetrahedron), hybrid 13 (nodal pressure); tet10 16; hex20 23. Shells: quad4 16 (fully
/// integrated), reduced 2 (one-point) with hourglass control; tri3 4 (C0 triangle); quad8 23;
/// tri6 24. Rods: 3 (truss). Beams: 2 (the Belytschko–Schwer resultant beam, whose section is
/// its area, second moments, torsion constant and shear area; no integration rule applies). 2D (in the XY plane): plane stress 12, plane strain 13, both one-point
/// forms on quad4 (written for reduced integration, with hourglass control) and tri3 (exact at one
/// point); quadratic 2D elements and axisymmetric sections have no form written here.
/// What it cannot express is the reason, naming the section.
pub(crate) fn elform(kind: ElementKind, section: &SectionKind, f: Formulation, section_name: &str) -> Result<Elform, String> {
    let plain = |elform| Ok(Elform { elform, hourglass: None });
    let reduced = f.integration == Integration::Reduced;
    let no = |what: &str| Err(format!("section `{section_name}`: {what} {}", kind.name()));
    match section {
        SectionKind::Solid {} => {
            if f.hybrid {
                return match kind {
                    ElementKind::Tet4 => plain(13),
                    _ => no(
                        "the hybrid (mixed displacement-pressure) formulation has a counterpart only for 4-node tetrahedra (nodal pressure, ELFORM 13), not",
                    ),
                };
            }
            match kind {
                ElementKind::Hex8 if f.incompatible_modes => plain(-18),
                ElementKind::Hex8 if reduced => Ok(Elform { elform: 1, hourglass: Some(Hourglass::Solid) }),
                ElementKind::Hex8 => plain(2),
                ElementKind::Wedge6 if reduced => plain(15),
                ElementKind::Wedge6 => plain(2),
                ElementKind::Pyramid5 if reduced => no("there is no reduced-integration form of"),
                ElementKind::Pyramid5 => plain(2),
                ElementKind::Tet4 => plain(10),
                _ if reduced && kind.is_quadratic() => no("there is no reduced-integration form of"),
                ElementKind::Tet10 => plain(16),
                ElementKind::Hex20 => plain(23),
                _ => no("a solid section cannot hold"),
            }
        }
        SectionKind::Shell { .. } => match kind {
            ElementKind::Quad4 if reduced => Ok(Elform { elform: 2, hourglass: Some(Hourglass::Shell) }),
            ElementKind::Quad4 => plain(16),
            ElementKind::Tri3 => plain(4),
            _ if reduced && kind.is_quadratic() => {
                no("reduced integration of quadratic shells is a whole-model setting in LS-DYNA, not written per section, for")
            }
            ElementKind::Quad8 => plain(23),
            ElementKind::Tri6 => plain(24),
            _ => no("a shell section cannot hold"),
        },
        SectionKind::Rod { .. } => match kind {
            ElementKind::Rod2 => plain(3),
            ElementKind::Rod3 => no("LS-DYNA's trusses have two nodes; there is no counterpart for"),
            _ => no("a rod section cannot hold"),
        },
        SectionKind::Beam { .. } => match kind {
            ElementKind::Beam2 => plain(2),
            _ => no("a beam section cannot hold"),
        },
        SectionKind::Axisymmetric {} => {
            no("axisymmetric sections are not written to LS-DYNA decks yet (its axisymmetric forms count loads per radian; not verified):")
        }
        SectionKind::PlaneStress { .. } | SectionKind::PlaneStrain { .. } => {
            let form = if matches!(section, SectionKind::PlaneStress { .. }) { 12 } else { 13 };
            match kind {
                ElementKind::Tri3 => plain(form),
                ElementKind::Quad4 if reduced => Ok(Elform { elform: form, hourglass: Some(Hourglass::Shell) }),
                ElementKind::Quad4 => no("LS-DYNA's 2D forms integrate at one point (set the section's integration to reduced) for"),
                _ => no("LS-DYNA's 2D forms are linear; there is no form for"),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_kind_maps_every_node_once_except_degenerate_repeats() {
        for kind in ElementKind::ALL {
            let o = order(kind);
            let mut seen: Vec<usize> = o.to_vec();
            seen.sort_unstable();
            seen.dedup();
            assert_eq!(seen, (0..kind.node_count()).collect::<Vec<_>>(), "{}", kind.name());
            let nodes: Vec<u32> = (1..=kind.node_count() as u32).collect();
            assert_eq!(connectivity(kind, &nodes).unwrap().len(), o.len());
            assert!(connectivity(kind, &nodes[..1]).is_none());
        }
    }

    #[test]
    fn formulations_map_or_say_why_not() {
        let solid = SectionKind::Solid {};
        let full = Formulation::default();
        let reduced = Formulation { integration: Integration::Reduced, ..full };
        let im = Formulation { incompatible_modes: true, ..full };
        let hybrid = Formulation { hybrid: true, ..full };
        assert_eq!(elform(ElementKind::Hex8, &solid, full, "s").unwrap().elform, 2);
        assert_eq!(elform(ElementKind::Hex8, &solid, reduced, "s").unwrap(), Elform { elform: 1, hourglass: Some(Hourglass::Solid) });
        assert_eq!(elform(ElementKind::Hex8, &solid, im, "s").unwrap().elform, -18);
        assert_eq!(elform(ElementKind::Tet4, &solid, hybrid, "s").unwrap().elform, 13);
        assert_eq!(elform(ElementKind::Tet10, &solid, full, "s").unwrap().elform, 16);
        assert_eq!(elform(ElementKind::Hex20, &solid, full, "s").unwrap().elform, 23);
        for (k, f) in [(ElementKind::Hex8, hybrid), (ElementKind::Tet10, reduced), (ElementKind::Hex20, reduced)] {
            let e = elform(k, &solid, f, "web").unwrap_err();
            assert!(e.contains("section `web`") && e.contains(k.name()), "{e}");
        }
        let shell = SectionKind::Shell { thickness: 1.0, offset: 0.0, points: 5 };
        assert_eq!(elform(ElementKind::Quad4, &shell, full, "s").unwrap().elform, 16);
        assert_eq!(elform(ElementKind::Quad4, &shell, reduced, "s").unwrap().hourglass, Some(Hourglass::Shell));
        assert_eq!(elform(ElementKind::Tri6, &shell, full, "s").unwrap().elform, 24);
        assert!(elform(ElementKind::Quad8, &shell, reduced, "s").is_err());
        assert!(elform(ElementKind::Hex8, &shell, full, "s").is_err());
        assert_eq!(elform(ElementKind::Rod2, &SectionKind::Rod { area: 1.0 }, full, "s").unwrap().elform, 3);
        let beam = SectionKind::Beam { profile: fieldcraft_model::BeamProfile::Circle { radius: 1.0 }, orientation: [0.0, 0.0, 1.0] };
        assert_eq!(elform(ElementKind::Beam2, &beam, full, "s").unwrap(), Elform { elform: 2, hourglass: None });
        assert_eq!(elform(ElementKind::Beam2, &beam, reduced, "s").unwrap().elform, 2, "a resultant beam has no integration rule");
        assert!(elform(ElementKind::Rod2, &beam, full, "s").unwrap_err().contains("a beam section cannot hold rod2"));
    }
}
