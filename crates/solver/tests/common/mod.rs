//! Structured meshes generated in code for the analytic verification suite: a block of
//! hexahedra (hex8, hex20), of wedges (two per cell) or of Kuhn-split tetrahedra (tet4, tet10) on [0, 1]³, mapped to the part
//! by a function, and helpers to define the model on it.
#![allow(dead_code)]

use std::collections::BTreeMap;

use fieldcraft_model::{
    BoundaryCondition, Coupling, CouplingKind, Dofs, Elastic, ElementFace, ElementKind, Formulation, Integration, Load, Material, Model, Procedure,
    ReferencePoint, Section, SectionKind, Set, Step, Symmetry, Target,
};

pub struct Grid {
    pub model: Model,
    /// Parametric coordinates of every node.
    pub param: BTreeMap<u32, [f64; 3]>,
    pub kind: ElementKind,
}

/// A structured mesh with `n` cells per direction. Nodes live on the doubled lattice
/// (2n + 1 points per direction) so quadratic elements get their mid-side nodes from the map
/// too (curved edges follow it).
pub fn grid(n: [usize; 3], kind: ElementKind, map: impl Fn([f64; 3]) -> [f64; 3]) -> Grid {
    let size = n.map(|k| 2 * k + 1);
    let label = |i: [usize; 3]| (1 + i[0] + size[0] * (i[1] + size[1] * i[2])) as u32;
    let param_of = |i: [usize; 3]| [0, 1, 2].map(|k| i[k] as f64 / (2 * n[k]) as f64);
    let mut model = Model::default();
    let mut param = BTreeMap::new();
    let mut node = |model: &mut Model, i: [usize; 3]| -> u32 {
        let l = label(i);
        if let std::collections::btree_map::Entry::Vacant(slot) = param.entry(l) {
            let p = param_of(i);
            model.set_node(l, map(p)).unwrap();
            slot.insert(p);
        }
        l
    };
    let mut next = 1u32;
    for cz in 0..n[2] {
        for cy in 0..n[1] {
            for cx in 0..n[0] {
                let base = [2 * cx, 2 * cy, 2 * cz];
                let at = |d: [usize; 3]| [base[0] + d[0], base[1] + d[1], base[2] + d[2]];
                match kind {
                    ElementKind::Hex8 | ElementKind::Hex20 => {
                        let corners = [[0, 0, 0], [2, 0, 0], [2, 2, 0], [0, 2, 0], [0, 0, 2], [2, 0, 2], [2, 2, 2], [0, 2, 2]];
                        let mut nodes: Vec<u32> = corners.iter().map(|c| node(&mut model, at(*c))).collect();
                        if kind == ElementKind::Hex20 {
                            for e in kind.edges() {
                                let (a, b) = (corners[e[0]], corners[e[1]]);
                                nodes.push(node(&mut model, at([0, 1, 2].map(|k| (a[k] + b[k]) / 2))));
                            }
                        }
                        model.set_element(next, kind, nodes).unwrap();
                        next += 1;
                    }
                    ElementKind::Wedge6 => {
                        // Two prisms along z, split by the cell's (x, y) diagonal, both counter-clockwise from above.
                        for tri in [[[0, 0], [2, 0], [2, 2]], [[0, 0], [2, 2], [0, 2]]] {
                            let nodes: Vec<u32> =
                                [0, 2].iter().flat_map(|&z| tri.map(|t| [t[0], t[1], z])).map(|c| node(&mut model, at(c))).collect();
                            model.set_element(next, kind, nodes).unwrap();
                            next += 1;
                        }
                    }
                    _ => {
                        // Six tetrahedra around the cell's main diagonal (Kuhn), conforming across cells.
                        for perm in [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]] {
                            let mut v = [[0usize; 3]; 4];
                            for s in 0..3 {
                                v[s + 1] = v[s];
                                v[s + 1][perm[s]] = 2;
                            }
                            // Orientation from the parametric volume.
                            let p = v.map(|c| c.map(|x| x as f64));
                            let e = |i: usize| [0, 1, 2].map(|k| p[i][k] - p[0][k]);
                            let (a, b, c) = (e(1), e(2), e(3));
                            let vol = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
                            if vol < 0.0 {
                                v.swap(1, 2);
                            }
                            let mut nodes: Vec<u32> = v.iter().map(|c| node(&mut model, at(*c))).collect();
                            if kind == ElementKind::Tet10 {
                                for e in kind.edges() {
                                    let (a, b) = (v[e[0]], v[e[1]]);
                                    nodes.push(node(&mut model, at([0, 1, 2].map(|k| (a[k] + b[k]) / 2))));
                                }
                            }
                            model.set_element(next, kind, nodes).unwrap();
                            next += 1;
                        }
                    }
                }
            }
        }
    }
    Grid { model, param, kind }
}

/// A hex–tet transition mesh on [0, 1]³ with `n` cells per direction: the cells with x below
/// one half are hex8, every other cell is split from its centre into six pyramids over its faces,
/// and each pyramid whose face is not shared with a hex cell into two tetrahedra (split along
/// the diagonal through the face's lowest node label, so neighbours agree). The mesh is
/// conforming: hexahedra meet pyramids on whole quadrilaterals, everything else on triangles.
pub fn transition_grid(n: [usize; 3], map: impl Fn([f64; 3]) -> [f64; 3]) -> Grid {
    let size = n.map(|k| 2 * k + 1);
    let label = |i: [usize; 3]| (1 + i[0] + size[0] * (i[1] + size[1] * i[2])) as u32;
    let param_of = |i: [usize; 3]| [0, 1, 2].map(|k| i[k] as f64 / (2 * n[k]) as f64);
    let mut model = Model::default();
    let mut param = BTreeMap::new();
    let mut node = |model: &mut Model, i: [usize; 3]| -> u32 {
        let l = label(i);
        if let std::collections::btree_map::Entry::Vacant(slot) = param.entry(l) {
            let p = param_of(i);
            model.set_node(l, map(p)).unwrap();
            slot.insert(p);
        }
        l
    };
    let is_hex = |c: [isize; 3]| (0..3).all(|k| c[k] >= 0 && (c[k] as usize) < n[k]) && 2 * c[0] < n[0] as isize;
    const CORNERS: [[usize; 3]; 8] = [[0, 0, 0], [2, 0, 0], [2, 2, 0], [0, 2, 0], [0, 0, 2], [2, 0, 2], [2, 2, 2], [0, 2, 2]];
    // Hex8's faces (outward) and the cell across each.
    let across: [[isize; 3]; 6] = [[0, 0, -1], [0, 0, 1], [0, -1, 0], [1, 0, 0], [0, 1, 0], [-1, 0, 0]];
    let mut next = 1u32;
    for cz in 0..n[2] {
        for cy in 0..n[1] {
            for cx in 0..n[0] {
                let base = [2 * cx, 2 * cy, 2 * cz];
                let at = |d: [usize; 3]| [base[0] + d[0], base[1] + d[1], base[2] + d[2]];
                let cell = [cx as isize, cy as isize, cz as isize];
                let corners: Vec<u32> = CORNERS.iter().map(|c| node(&mut model, at(*c))).collect();
                if is_hex(cell) {
                    model.set_element(next, ElementKind::Hex8, corners).unwrap();
                    next += 1;
                    continue;
                }
                let apex = node(&mut model, at([1, 1, 1]));
                for (face, d) in ElementKind::Hex8.faces().iter().zip(across) {
                    // Reversed, the outward face's normal points into the cell, at the apex.
                    let q: Vec<u32> = face.iter().rev().map(|&i| corners[i]).collect();
                    if is_hex([cell[0] + d[0], cell[1] + d[1], cell[2] + d[2]]) {
                        model.set_element(next, ElementKind::Pyramid5, vec![q[0], q[1], q[2], q[3], apex]).unwrap();
                        next += 1;
                        continue;
                    }
                    let lowest = (0..4).min_by_key(|&i| q[i]).unwrap();
                    let (a, b, c, e) = if lowest % 2 == 0 { (q[0], q[1], q[2], q[3]) } else { (q[1], q[2], q[3], q[0]) };
                    for tet in [[a, b, c, apex], [a, c, e, apex]] {
                        model.set_element(next, ElementKind::Tet4, tet.to_vec()).unwrap();
                        next += 1;
                    }
                }
            }
        }
    }
    Grid { model, param, kind: ElementKind::Pyramid5 }
}

pub fn near(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

impl Grid {
    pub fn nodes_where(&self, f: impl Fn([f64; 3]) -> bool) -> Vec<u32> {
        self.param.iter().filter(|(_, p)| f(**p)).map(|(l, _)| *l).collect()
    }

    /// Element faces whose nodes all satisfy `f` (on parametric coordinates).
    pub fn faces_where(&self, f: impl Fn([f64; 3]) -> bool) -> Vec<ElementFace> {
        let mut out = Vec::new();
        for (label, e) in self.model.elements() {
            for (i, face) in e.kind.faces().iter().enumerate() {
                if face.iter().all(|k| f(self.param[&e.nodes[*k]])) {
                    out.push(ElementFace { element: *label, face: i as u8 + 1 });
                }
            }
        }
        out
    }

    /// Material, section over every element and an empty step of the given procedure.
    pub fn define(&mut self, elastic: Elastic, density: Option<f64>, formulation: Formulation, procedure: Procedure) {
        let all: Vec<u32> = self.model.elements().keys().copied().collect();
        self.model.add_set("all", Set { elements: all, ..Set::default() }).unwrap();
        self.model.add_material("Mat", Material { elastic: Some(elastic), density, ..Material::default() }).unwrap();
        self.model
            .add_section("solid", Section { orientation: None, set: "all".into(), material: "Mat".into(), kind: SectionKind::Solid {}, formulation })
            .unwrap();
        self.model.add_step(Step::new("Step", procedure), None).unwrap();
    }

    pub fn node_set(&mut self, name: &str, f: impl Fn([f64; 3]) -> bool) {
        let nodes = self.nodes_where(f);
        self.model.add_set(name, Set { nodes, ..Set::default() }).unwrap();
    }

    pub fn face_set(&mut self, name: &str, f: impl Fn([f64; 3]) -> bool) {
        let faces = self.faces_where(f);
        self.model.add_set(name, Set { faces, ..Set::default() }).unwrap();
    }

    pub fn bc(&mut self, name: &str, bc: BoundaryCondition) {
        self.model.add_boundary_condition(name, bc).unwrap();
        self.model.set_boundary_condition_active("Step", name, true).unwrap();
    }

    pub fn load(&mut self, name: &str, load: Load) {
        self.model.add_load(name, load).unwrap();
        self.model.set_load_active("Step", name, true).unwrap();
    }

    pub fn symmetry(&mut self, name: &str, plane: Symmetry, f: impl Fn([f64; 3]) -> bool) {
        self.node_set(name, f);
        self.bc(name, BoundaryCondition::Symmetry { target: Target::Set(name.into()), plane });
    }

    /// A reference point `name` at `xyz`, coupled to the nodes where `f` holds.
    pub fn couple(&mut self, name: &str, xyz: [f64; 3], kind: CouplingKind, dofs: &[u8], f: impl Fn([f64; 3]) -> bool) {
        let set = format!("{name}_nodes");
        self.node_set(&set, f);
        self.model.add_reference_point(name, ReferencePoint { xyz, label: 0 }).unwrap();
        let dofs = Dofs::try_from(dofs.to_vec()).unwrap();
        self.model.add_coupling(&format!("{name}_coupling"), Coupling { reference_point: name.into(), set, kind, dofs }).unwrap();
    }

    pub fn fixed(&mut self, name: &str, f: impl Fn([f64; 3]) -> bool) {
        self.node_set(name, f);
        self.bc(name, BoundaryCondition::Fixed { target: Target::Set(name.into()) });
    }
}

pub fn iso(young: f64, poisson: f64) -> Elastic {
    Elastic::Isotropic { young, poisson }
}

pub fn full() -> Formulation {
    Formulation::default()
}

pub fn reduced() -> Formulation {
    Formulation { integration: Integration::Reduced, ..Formulation::default() }
}

pub fn incompatible() -> Formulation {
    Formulation { incompatible_modes: true, ..Formulation::default() }
}

/// The element kinds and formulations under test.
pub fn variants() -> Vec<(ElementKind, Formulation, &'static str)> {
    vec![
        (ElementKind::Tet4, full(), "tet4"),
        (ElementKind::Tet10, full(), "tet10"),
        (ElementKind::Wedge6, full(), "wedge6"),
        (ElementKind::Hex8, full(), "hex8"),
        (ElementKind::Hex8, incompatible(), "hex8 incompatible"),
        (ElementKind::Hex20, full(), "hex20"),
        (ElementKind::Hex20, reduced(), "hex20 reduced"),
    ]
}
