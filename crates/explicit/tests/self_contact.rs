//! Self contact (single-surface contact) in explicit dynamics: general contact on a body's own
//! exterior faces, every face against every other one of the same surface. A body that only bends
//! or squeezes must never feel a force from its own faces.
//!
//! - A layer of elements squeezed through its thickness moves exactly as it does without contact
//!   (issue #1034: a layer one hex8 thick squeezed by 60 % took 193 N, one two elements thick
//!   squeezed by 80 %, 1 216 N, because the depth limit used the undeformed element thickness and a
//!   hexahedron's opposite faces share no node).

use std::collections::BTreeMap;

use fieldcraft_explicit::{ExplicitResult, Limits, solve_explicit};
use fieldcraft_model::*;

type V3 = [f64; 3];

/// A model built from hexahedra whose nodes are merged by position within a body tag: blocks of
/// one tag that touch share their nodes (one body), blocks of different tags do not.
struct Mesh {
    model: Model,
    nodes: BTreeMap<(u8, [i64; 3]), u32>,
    next_node: u32,
    next_element: u32,
}

const QUANTUM: f64 = 1e-9;

impl Mesh {
    fn new() -> Mesh {
        Mesh { model: Model::default(), nodes: BTreeMap::new(), next_node: 1, next_element: 1 }
    }

    fn node(&mut self, body: u8, x: V3) -> u32 {
        let key = (body, x.map(|c| (c / QUANTUM).round() as i64));
        if let Some(n) = self.nodes.get(&key) {
            return *n;
        }
        let l = self.next_node;
        self.model.set_node(l, x).unwrap();
        self.nodes.insert(key, l);
        self.next_node += 1;
        l
    }

    /// One hex8 of corners (body tag, position) in the element's order.
    fn hex(&mut self, corners: [(u8, V3); 8]) -> u32 {
        let nodes: Vec<u32> = corners.iter().map(|(b, x)| self.node(*b, *x)).collect();
        let e = self.next_element;
        self.model.set_element(e, ElementKind::Hex8, nodes).unwrap();
        self.next_element += 1;
        e
    }

    /// A box from `origin` of `size` with `n` hex8 along each axis, positions mapped by `map`.
    fn block(&mut self, body: u8, origin: V3, size: V3, n: [usize; 3], map: impl Fn(V3) -> V3) -> Vec<u32> {
        let at = |i: usize, j: usize, k: usize| map([0, 1, 2].map(|d| origin[d] + size[d] * [i, j, k][d] as f64 / n[d] as f64));
        let mut out = Vec::new();
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    out.push(self.hex([
                        (body, at(i, j, k)),
                        (body, at(i + 1, j, k)),
                        (body, at(i + 1, j + 1, k)),
                        (body, at(i, j + 1, k)),
                        (body, at(i, j, k + 1)),
                        (body, at(i + 1, j, k + 1)),
                        (body, at(i + 1, j + 1, k + 1)),
                        (body, at(i, j + 1, k + 1)),
                    ]));
                }
            }
        }
        out
    }

    fn xyz(&self, node: u32) -> V3 {
        self.model.node(node).unwrap().xyz
    }

    fn nodes_where(&self, elements: &[u32], f: impl Fn(V3) -> bool) -> Vec<u32> {
        let mut out: Vec<u32> = elements.iter().flat_map(|e| self.model.element(*e).unwrap().nodes.clone()).filter(|n| f(self.xyz(*n))).collect();
        out.sort_unstable();
        out.dedup();
        out
    }

    fn material(&mut self, name: &str, young: f64, poisson: f64, density: f64, hyper: bool) {
        let m = Material {
            elastic: Some(Elastic::Isotropic { young, poisson }),
            density: Some(density),
            hyperelastic: hyper.then_some(Hyperelastic::NeoHookean {}),
            ..Material::default()
        };
        self.model.add_material(name, m).unwrap();
    }

    fn section(&mut self, name: &str, elements: &[u32], material: &str) {
        self.model.add_set(name, Set { elements: elements.to_vec(), ..Set::default() }).unwrap();
        let formulation = Formulation { integration: Integration::Reduced, ..Formulation::default() };
        self.model.add_section(name, Section { set: name.into(), material: material.into(), kind: SectionKind::Solid {}, formulation }).unwrap();
    }

    fn node_set(&mut self, name: &str, nodes: Vec<u32>) {
        assert!(!nodes.is_empty(), "{name}");
        self.model.add_set(name, Set { nodes, ..Set::default() }).unwrap();
    }

    fn general(&mut self, name: &str, surface: Option<&str>, friction: Friction) {
        let c = Contact { friction, ..Contact::new(ContactSurfaces::General { surface: surface.map(String::from) }) };
        self.model.add_contact(name, c).unwrap();
    }

    fn step(&mut self, step: ExplicitDynamic, bcs: &[(&str, BoundaryCondition)]) {
        let mut s = Step::new("Run", Procedure::ExplicitDynamic(step));
        for (name, bc) in bcs {
            self.model.add_boundary_condition(name, bc.clone()).unwrap();
            s.boundary_conditions.insert(name.to_string());
        }
        self.model.add_step(s, None).unwrap();
    }

    fn solve(&self) -> ExplicitResult {
        let r = solve_explicit(&self.model, Some("Run"), &Limits::default()).unwrap();
        assert!(r.completed, "{:?}", r.warnings);
        r
    }
}

fn held(set: &str, dofs: &[usize], value: f64) -> BoundaryCondition {
    let mut values = [None; 6];
    for k in dofs {
        values[*k] = Some(value);
    }
    BoundaryCondition::Displacement { target: Target::Set(set.into()), values }
}

/// Run `m` as it is and without its contacts: the same motion, bit for bit, and no contact force.
fn same_without_contact(mut m: Mesh, what: &str) {
    let with = m.solve();
    let names: Vec<String> = m.model.contacts().keys().cloned().collect();
    for n in names {
        m.model.remove_contact(&n).unwrap();
    }
    let without = m.solve();
    let peak = with.contacts.iter().map(|s| s.peak_force).fold(0.0, f64::max);
    eprintln!("{what}: peak contact force {peak:.3e}");
    assert!(with.contacts.iter().all(|s| s.first_contact.is_none()), "{what}: {:?}", with.contacts);
    assert_eq!(with.frames, without.frames, "{what}");
}

/// A layer of `layers` elements through its thickness squeezed by the fraction `squeeze` in a
/// closed die (no lateral motion, uniaxial strain): the far sides of its elements come closer than
/// half an undeformed element, but they are still the far sides, not a surface the near side has
/// passed. It moves as it does without contact and takes no contact force.
fn squeezed_layer(layers: usize, squeeze: f64) {
    let mut m = Mesh::new();
    let t = 0.1;
    let layer = m.block(0, [0.0; 3], [0.4, 0.4, t], [4, 4, layers], |x| x);
    m.material("Foam", 1.0e3, 0.0, 1.0, true);
    m.section("layer", &layer, "Foam");
    m.node_set("bottom", m.nodes_where(&layer, |x| x[2].abs() < 1e-9));
    m.node_set("top", m.nodes_where(&layer, |x| (x[2] - t).abs() < 1e-9));
    m.node_set("all_nodes", m.nodes_where(&layer, |_| true));
    m.general("self", None, Friction::NONE);
    let c = (1.0e3f64).sqrt();
    let end = 200.0 * t / c;
    let step = ExplicitDynamic::new(end);
    m.step(step, &[("die", held("all_nodes", &[0, 1], 0.0)), ("bottom", held("bottom", &[2], 0.0)), ("press", held("top", &[2], -squeeze * t))]);
    same_without_contact(m, &format!("{layers} layer(s) squeezed {squeeze}"));
}

/// One hex8 thick, squeezed by 60 %: before the depth limit followed the current element, its far
/// face took a contact force of 193 N (issue #1034).
#[test]
fn a_squeezed_one_element_layer_feels_no_self_contact() {
    squeezed_layer(1, 0.6);
}

/// Two hex8 thick, squeezed by 80 %: before the fix, 1 216 N (issue #1034).
#[test]
fn a_squeezed_two_element_layer_feels_no_self_contact() {
    squeezed_layer(2, 0.8);
}
