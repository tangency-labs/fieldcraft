//! Contact in nonlinear static steps against closed forms: the contact patch test on
//! non-matching meshes, Newton's convergence with contact, and what implicit contact refuses;
//! ties of non-matching meshes for displacement and for heat.
//! Run with `--nocapture` to see the iterations.

use std::collections::BTreeMap;

use fieldcraft_model::*;
use fieldcraft_solver::contact::Status;
use fieldcraft_solver::{Error, Limits, StaticResult, solve_nonlinear};

/// A model built block by block.
#[derive(Default)]
struct Blocks {
    model: Model,
    xyz: BTreeMap<u32, [f64; 3]>,
    /// Each corner node's point in its block's unit cube.
    params: BTreeMap<u32, [f64; 3]>,
    next_node: u32,
    next_element: u32,
}

impl Blocks {
    /// A structured block of `n` cells over the unit cube mapped by `map`, hex8 or tet4 (six
    /// Kuhn tetrahedra per cell), or hex20 or tet10 with straight edges (mid-side nodes at the
    /// middles of the mapped corners); its elements' labels.
    fn block(&mut self, kind: ElementKind, n: [usize; 3], map: impl Fn([f64; 3]) -> [f64; 3]) -> Vec<u32> {
        self.quadratic_block(kind, n, false, map)
    }

    /// As [`Blocks::block`], with `curved` the mid-side nodes of hex20 and tet10 at the image of
    /// their edges' middles in the unit cube, so the faces follow `map` (a mapped circle stays
    /// on the circle at every node).
    fn quadratic_block(&mut self, kind: ElementKind, n: [usize; 3], curved: bool, map: impl Fn([f64; 3]) -> [f64; 3]) -> Vec<u32> {
        let linear = match kind {
            ElementKind::Hex20 => ElementKind::Hex8,
            ElementKind::Tet10 => ElementKind::Tet4,
            k => k,
        };
        let out = self.linear_block(linear, n, &map);
        if linear != kind {
            let mut mids: BTreeMap<(u32, u32), u32> = BTreeMap::new();
            for l in &out {
                let mut nodes = self.model.element(*l).unwrap().nodes.clone();
                for e in kind.edges() {
                    let (a, b) = (nodes[e[0]], nodes[e[1]]);
                    let m = *mids.entry((a.min(b), a.max(b))).or_insert_with(|| {
                        let x = if curved {
                            map([0, 1, 2].map(|k| 0.5 * (self.params[&a][k] + self.params[&b][k])))
                        } else {
                            [0, 1, 2].map(|k| 0.5 * (self.xyz[&a][k] + self.xyz[&b][k]))
                        };
                        self.add_node(x)
                    });
                    nodes.push(m);
                }
                self.model.set_element(*l, kind, nodes).unwrap();
            }
        }
        out
    }

    /// A new node at `x`; its label.
    fn add_node(&mut self, x: [f64; 3]) -> u32 {
        self.next_node += 1;
        self.model.set_node(self.next_node, x).unwrap();
        self.xyz.insert(self.next_node, x);
        self.next_node
    }

    /// The first node (lowest label) at `x`.
    fn node_at(&self, x: [f64; 3]) -> Option<u32> {
        self.nodes(|y| (0..3).all(|k| near(x[k], y[k]))).first().copied()
    }

    fn linear_block(&mut self, kind: ElementKind, n: [usize; 3], map: impl Fn([f64; 3]) -> [f64; 3]) -> Vec<u32> {
        let first = self.next_node + 1;
        let label = |i: [usize; 3]| first + (i[0] + (n[0] + 1) * (i[1] + (n[1] + 1) * i[2])) as u32;
        for k in 0..=n[2] {
            for j in 0..=n[1] {
                for i in 0..=n[0] {
                    let p = [i as f64 / n[0] as f64, j as f64 / n[1] as f64, k as f64 / n[2] as f64];
                    let x = map(p);
                    let l = label([i, j, k]);
                    self.params.insert(l, p);
                    self.model.set_node(l, x).unwrap();
                    self.xyz.insert(l, x);
                    self.next_node = self.next_node.max(l);
                }
            }
        }
        let mut out = Vec::new();
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    let at = |d: [usize; 3]| label([i + d[0], j + d[1], k + d[2]]);
                    let cells: Vec<Vec<u32>> = match kind {
                        ElementKind::Hex8 => vec![
                            [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]].iter().map(|d| at(*d)).collect(),
                        ],
                        _ => [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]
                            .iter()
                            .map(|perm| {
                                let mut v = [[0usize; 3]; 4];
                                for s in 0..3 {
                                    v[s + 1] = v[s];
                                    v[s + 1][perm[s]] = 1;
                                }
                                let p = v.map(|c| c.map(|x| x as f64));
                                let e = |i: usize| [0, 1, 2].map(|k| p[i][k] - p[0][k]);
                                let (a, b, c) = (e(1), e(2), e(3));
                                let vol =
                                    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
                                if vol < 0.0 {
                                    v.swap(1, 2);
                                }
                                v.iter().map(|d| at(*d)).collect()
                            })
                            .collect(),
                    };
                    for nodes in cells {
                        self.next_element += 1;
                        self.model.set_element(self.next_element, kind, nodes).unwrap();
                        out.push(self.next_element);
                    }
                }
            }
        }
        out
    }

    /// The faces of `elements` whose nodes all satisfy `f`.
    fn faces(&self, elements: &[u32], f: impl Fn([f64; 3]) -> bool) -> Vec<ElementFace> {
        let mut out = Vec::new();
        for l in elements {
            let e = self.model.element(*l).unwrap();
            for (i, face) in e.kind.faces().iter().enumerate() {
                if face.iter().all(|k| f(self.xyz[&e.nodes[*k]])) {
                    out.push(ElementFace { element: *l, face: i as u8 + 1 });
                }
            }
        }
        out
    }

    fn nodes(&self, f: impl Fn([f64; 3]) -> bool) -> Vec<u32> {
        self.xyz.iter().filter(|(_, x)| f(**x)).map(|(l, _)| *l).collect()
    }

    fn node_set(&mut self, name: &str, f: impl Fn([f64; 3]) -> bool) {
        let nodes = self.nodes(f);
        self.model.add_set(name, Set { nodes, ..Set::default() }).unwrap();
    }

    fn face_set(&mut self, name: &str, elements: &[u32], f: impl Fn([f64; 3]) -> bool) {
        let faces = self.faces(elements, f);
        assert!(!faces.is_empty(), "no faces for `{name}`");
        self.model.add_set(name, Set { faces, ..Set::default() }).unwrap();
    }

    fn bc(&mut self, name: &str, values: [Option<f64>; 3], f: impl Fn([f64; 3]) -> bool) {
        self.node_set(name, f);
        let [a, b, c] = values;
        self.model
            .add_boundary_condition(name, BoundaryCondition::Displacement { target: Target::Set(name.into()), values: [a, b, c, None, None, None] })
            .unwrap();
        self.model.set_boundary_condition_active("Step", name, true).unwrap();
    }

    /// One elastic material and section over every element, and the nonlinear step.
    fn define(&mut self, young: f64, poisson: f64, settings: NonlinearStatic) {
        let all: Vec<u32> = self.model.elements().keys().copied().collect();
        self.model.add_set("all", Set { elements: all, ..Set::default() }).unwrap();
        self.model.add_material("Mat", Material { elastic: Some(Elastic::Isotropic { young, poisson }), ..Material::default() }).unwrap();
        self.model
            .add_section(
                "solid",
                Section { set: "all".into(), material: "Mat".into(), kind: SectionKind::Solid {}, formulation: Formulation::default() },
            )
            .unwrap();
        self.model.add_step(Step::new("Step", Procedure::NonlinearStatic(settings)), None).unwrap();
    }

    fn pair(&mut self, name: &str, main: &str, secondary: &str) {
        let c = Contact::new(ContactSurfaces::Pair { main: main.into(), secondary: secondary.into() });
        self.model.add_contact(name, c).unwrap();
    }
}

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

/// A flat contact whose secondary surface also contains the cube's uncontacted sides, of linear
/// and quadratic faces, in small and finite sliding. Those sides must not tilt the normal of the
/// nodes the bottom face carries: the normal force is the assembled force, every bottom node
/// (corners included) carries pressure and the others none, and each bottom node's pressure is the force
/// over the unit area.
#[test]
fn uncontacted_corner_faces_do_not_tilt_reported_pressure() {
    for (kind, nonlinear_geometry) in [(ElementKind::Hex8, false), (ElementKind::Hex8, true), (ElementKind::Hex20, false), (ElementKind::Hex20, true)]
    {
        let mut b = Blocks::default();
        let low = b.block(ElementKind::Hex8, [1, 1, 1], |p| p);
        let upper = b.block(kind, [1, 1, 1], |p| [p[0], p[1], p[2] + 1.0]);
        b.define(1000.0, 0.0, NonlinearStatic { nonlinear_geometry, initial_increment: 1.0, ..NonlinearStatic::default() });
        b.face_set("main", &low, |p| near(p[2], 1.0));
        b.face_set("secondary", &upper, |_| true);
        b.pair("Interface", "main", "secondary");
        let base = b.model.element(low[0]).unwrap().nodes.clone();
        let rest: Vec<u32> = b.nodes(|_| true).into_iter().filter(|l| !base.contains(l)).collect();
        let bottom: Vec<u32> = rest.iter().copied().filter(|l| near(b.xyz[l][2], 1.0)).collect();
        let top: Vec<u32> = rest.iter().copied().filter(|l| near(b.xyz[l][2], 2.0)).collect();
        for (name, nodes, values) in [
            ("base", base, [Some(0.0), Some(0.0), Some(0.0), None, None, None]),
            ("upper", rest, [Some(0.0), Some(0.0), None, None, None, None]),
            ("top", top, [None, None, Some(-0.01), None, None, None]),
        ] {
            b.model.add_set(name, Set { nodes, ..Set::default() }).unwrap();
            b.model.add_boundary_condition(name, BoundaryCondition::Displacement { target: Target::Set(name.into()), values }).unwrap();
            b.model.set_boundary_condition_active("Step", name, true).unwrap();
        }
        let r = solve_nonlinear(&b.model, Some("Step"), &Limits::default()).unwrap();
        let contact = &r.contacts[0];
        let normal = contact.force[2];
        let case = format!("{} (nonlinear geometry {nonlinear_geometry})", kind.name());
        // About E·ε = 10 (less under finite strain, and by the augmentation's remaining gap).
        assert!(normal > 9.5 && normal <= 10.0, "{case}: physical compression force: {normal}");
        // The uncontacted sides would tilt the corners' normals by about 55° (a third of their
        // force lost); a quadratic face pressed to slightly different gaps bends by far less.
        assert!((contact.normal_force - normal).abs() < 1e-6 * normal, "{case}: normal force {} against {normal}", contact.normal_force);
        for (label, pressure) in contact.nodes.iter().zip(&contact.pressure) {
            if !bottom.contains(label) {
                assert_eq!(*pressure, 0.0, "{case}: node {label}");
            } else {
                // Uniform pressure on a unit-area interface equals its actual assembled force;
                // this checks reporting independently of the augmentation's stopping gap. Under
                // finite strain a hex20's corner and mid-side nodes differ by a few parts in 10⁷,
                // as much with the bottom face alone (compared below).
                let tolerance = if nonlinear_geometry { 1e-6 } else { 1e-9 };
                assert!((pressure - normal).abs() < tolerance * normal, "{case}: bottom node {label}: pressure {pressure}, physical force {normal}");
            }
        }
        // Selecting only the bottom face must give the same solution, including the gap
        // that controls augmentation and the pressure-field coefficients.
        let faces = b.faces(&upper, |p| near(p[2], 1.0));
        b.model.replace_set("secondary", Set { faces, ..Set::default() }).unwrap();
        let only_bottom = solve_nonlinear(&b.model, Some("Step"), &Limits::default()).unwrap();
        let other = &only_bottom.contacts[0];
        assert!((other.force[2] - normal).abs() < 1e-9 * normal, "{case}: bottom-only force {} against {normal}", other.force[2]);
        for (i, label) in contact.nodes.iter().enumerate().filter(|(_, l)| bottom.contains(l)) {
            let j = other.nodes.iter().position(|l| l == label).unwrap();
            assert!((contact.pressure[i] - other.pressure[j]).abs() < 1e-9 * normal, "{case}: pressure at {label}");
            assert!((contact.multiplier[i] - other.multiplier[j]).abs() < 1e-9 * normal, "{case}: multiplier at {label}");
            assert!((contact.gap[i].unwrap() - other.gap[j].unwrap()).abs() < 1e-12, "{case}: gap at {label}");
        }
    }
}

fn report(r: &StaticResult) {
    let h = r.nonlinear.as_ref().unwrap();
    for inc in &h.increments {
        let forces: Vec<String> = inc.iterations.iter().map(|i| format!("{:.1e}", i.force)).collect();
        eprintln!("increment {} λ {:.3} aug {} [{}]", inc.number, inc.load_factor, inc.augmentations, forces.join(" "));
    }
}

/// Two blocks pressed together across non-matching meshes (`lower` × `upper` cells, each of
/// `lower_kind` / `upper_kind`), between frictionless symmetry planes: a uniaxial stress.
fn patch(lower_kind: ElementKind, upper_kind: ElementKind, geometry: bool, distorted: bool) -> (Blocks, StaticResult) {
    let mut b = Blocks::default();
    // A graded lower block, its cells boxes of different sizes, or distorted (its faces stay the
    // planes of the cube either way).
    let low = if distorted {
        b.block(lower_kind, [3, 3, 2], |p| [p[0] + 0.3 * p[0] * (1.0 - p[0]) * (p[1] - 0.5), p[1] + 0.3 * p[1] * (1.0 - p[1]) * (p[0] - 0.4), p[2]])
    } else {
        b.block(lower_kind, [3, 3, 2], |p| [p[0].powf(1.4), p[1].powf(0.8), p[2]])
    };
    let upper = [2, 2, 2];
    let up = b.block(upper_kind, upper, |p| [p[0], p[1], 1.0 + p[2]]);
    b.define(
        1000.0,
        0.3,
        NonlinearStatic {
            nonlinear_geometry: geometry,
            initial_increment: 0.5,
            // Tight for small sliding (linear: exact); finite sliding converges to the default ones.
            tolerances: if geometry { Tolerances::default() } else { Tolerances { force: 1e-10, displacement: 1e-8, energy: 1e-16 } },
            ..NonlinearStatic::default()
        },
    );
    b.face_set("lower_top", &low, |x| near(x[2], 1.0));
    b.face_set("upper_bottom", &up, |x| near(x[2], 1.0));
    b.pair("interface", "lower_top", "upper_bottom");
    b.bc("x0", [Some(0.0), None, None], |x| near(x[0], 0.0));
    b.bc("y0", [None, Some(0.0), None], |x| near(x[1], 0.0));
    b.bc("z0", [None, None, Some(0.0)], |x| near(x[2], 0.0));
    b.bc("push", [None, None, Some(-0.02)], |x| near(x[2], 2.0));
    let r = solve_nonlinear(&b.model, Some("Step"), &Limits::default()).unwrap();
    report(&r);
    (b, r)
}

/// The contact patch test (Taylor & Papadopoulos 1991; Crisfield 2000): a uniform pressure
/// crosses a non-matching interface exactly, so every element of both blocks has the same
/// uniaxial stress, with nothing across.
#[test]
fn the_mortar_interface_passes_the_contact_patch_test() {
    for (lower, upper, name) in [
        (ElementKind::Hex8, ElementKind::Hex8, "hex8 on hex8"),
        (ElementKind::Hex8, ElementKind::Tet4, "tet4 on hex8"),
        (ElementKind::Tet4, ElementKind::Hex8, "hex8 on tet4"),
    ] {
        for geometry in [false, true] {
            // Small sliding is linear and exact; finite sliding is exact up to the gaps the
            // augmentation leaves (their geometric terms, GAP_TOLERANCE) and the convergence.
            let tol = if geometry { 1e-4 } else { 1e-9 };
            let (b, r) = patch(lower, upper, geometry, false);
            let szz: Vec<f64> = r.elements.iter().flat_map(|e| e.nodal.iter().map(|s| s[2])).collect();
            let mean = szz.iter().sum::<f64>() / szz.len() as f64;
            assert!(mean < 0.0, "{name}: {mean}");
            for (e, s) in r.elements.iter().flat_map(|e| e.nodal.iter().map(move |s| (e.label, s))) {
                assert!((s[2] - mean).abs() < tol * mean.abs(), "{name} ({geometry}) element {e}: σzz {} vs {mean}", s[2]);
                assert!(s[0].abs() < tol * mean.abs() && s[1].abs() < tol * mean.abs(), "{name} ({geometry}) element {e}: {s:?}");
            }
            let c = &r.contacts[0];
            // The pressure is uniform, equal to the stress, at every secondary node.
            assert!(c.pressure.iter().all(|p| (p + mean).abs() < tol * mean.abs()), "{name}: {:?} vs {mean}", c.pressure);
            // The contact force on the upper block balances the reaction at its top.
            let top: f64 = b.nodes(|x| near(x[2], 2.0)).iter().map(|l| r.reaction[r.nodes.binary_search(l).unwrap()][2]).sum();
            assert!((c.force[2] + top).abs() < tol * top.abs(), "{name}: {:?} vs {top}", c.force);
            assert!(c.force[0].abs() < tol * top.abs() && c.force[1].abs() < tol * top.abs(), "{name}: {:?}", c.force);
            assert!((c.normal_force + top).abs() < tol * top.abs(), "{name}: {} vs {top}", c.normal_force);
            // Augmented to a gap tolerance far below the penetration a penalty alone leaves.
            assert!(c.penetration < 1e-5, "{name}: {}", c.penetration);
        }
    }
}

/// The contact patch test with quadratic faces: hex20 and tet10 blocks pressed together, against
/// each other and against linear ones, the upper block secondary, with small (`finite`: finite)
/// sliding. The transformed multiplier basis of the quadratic secondary faces carries the
/// uniform pressure at corners and mid-side nodes alike, and the facets of the segmentation are
/// the faces exactly (straight edges).
fn quadratic_patch(finite: bool) {
    // Small sliding is linear and exact; finite sliding is exact up to the gaps the augmentation
    // leaves and the convergence, as with linear faces.
    let tol = if finite { 2e-4 } else { 1e-9 };
    for (lower, upper, name) in [
        (ElementKind::Hex20, ElementKind::Hex20, "hex20 on hex20"),
        (ElementKind::Tet10, ElementKind::Tet10, "tet10 on tet10"),
        (ElementKind::Tet10, ElementKind::Hex20, "hex20 on tet10"),
        (ElementKind::Hex8, ElementKind::Hex20, "hex20 on hex8"),
        (ElementKind::Hex20, ElementKind::Tet4, "tet4 on hex20"),
    ] {
        let (b, r) = patch(lower, upper, finite, false);
        let szz: Vec<f64> = r.elements.iter().flat_map(|e| e.nodal.iter().map(|s| s[2])).collect();
        let mean = szz.iter().sum::<f64>() / szz.len() as f64;
        assert!(mean < 0.0, "{name}: {mean}");
        // σzz uniform and σxx, σyy zero in every element, over |σzz|.
        let spread =
            r.elements.iter().flat_map(|e| &e.nodal).map(|s| (s[2] - mean).abs().max(s[0].abs()).max(s[1].abs())).fold(0.0, f64::max) / mean.abs();
        let c = &r.contacts[0];
        // The pressure equal to the stress at every secondary node.
        let pressure = c.pressure.iter().map(|p| (p + mean).abs()).fold(0.0, f64::max) / mean.abs();
        let top: f64 = b.nodes(|x| near(x[2], 2.0)).iter().map(|l| r.reaction[r.nodes.binary_search(l).unwrap()][2]).sum();
        let force = (c.force[2] + top).abs().max((c.normal_force + top).abs()) / top.abs();
        eprintln!(
            "{name} (finite {finite}): σzz {mean:.6}, spread {spread:.1e}, pressure {pressure:.1e}, force {force:.1e}, penetration {:.1e}",
            c.penetration
        );
        assert!(spread < tol && pressure < tol && force < tol, "{name} (finite {finite}): {spread:e} {pressure:e} {force:e}");
        assert!(c.force[0].abs() < tol * top.abs() && c.force[1].abs() < tol * top.abs(), "{name}: {:?}", c.force);
        assert!(c.penetration < 1e-5, "{name}: {}", c.penetration);
    }
}

#[test]
fn quadratic_faces_pass_the_contact_patch_test() {
    quadratic_patch(false);
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn quadratic_faces_pass_the_contact_patch_test_with_finite_sliding() {
    quadratic_patch(true);
}

/// What quadratic contact refuses by name: ties of quadratic faces (dual mortar takes linear
/// faces), and a secondary face collapsed by a mid-side node moved onto its corner.
#[test]
fn what_quadratic_contact_refuses() {
    let mut q = Blocks::default();
    let low = q.block(ElementKind::Hex20, [1, 1, 1], |p| p);
    let up = q.block(ElementKind::Hex20, [1, 1, 1], |p| [p[0], p[1], 1.0 + p[2]]);
    q.define(1000.0, 0.3, NonlinearStatic::default());
    q.face_set("lower_top", &low, |x| near(x[2], 1.0));
    q.face_set("upper_bottom", &up, |x| near(x[2], 1.0));
    q.bc("z0", [Some(0.0), Some(0.0), Some(0.0)], |x| near(x[2], 0.0));
    q.bc("push", [Some(0.0), Some(0.0), Some(-0.01)], |x| near(x[2], 2.0));
    let mut tied = q.model.clone();
    tied.add_tie("glue", Tie::new("lower_top", "upper_bottom")).unwrap();
    let e = solve_nonlinear(&tied, Some("Step"), &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("hex20") && e.contains("linear faces"), "{e}");
    // A mid-side node of the secondary face moved onto a corner.
    q.pair("interface", "lower_top", "upper_bottom");
    let face = q.faces(&up, |x| near(x[2], 1.0))[0];
    let nodes = q.model.face_nodes(face).unwrap();
    let corner = q.xyz[&nodes[0]];
    q.model.set_node(nodes[4], corner).unwrap();
    // (The element under it is refused first, by name.)
    let e = solve_nonlinear(&q.model, Some("Step"), &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("hex20") && e.contains("degenerate"), "{e}");
}

/// A secondary surface whose tri6 and quad8 faces share a mid-side node is refused naming the
/// node: the two faces transform the multiplier basis by different α there.
#[test]
fn a_mid_side_node_of_a_tri6_and_a_quad8_face_is_refused() {
    let mut q = Blocks::default();
    // The upper block first: its nodes are the first found at the interface.
    let up = q.block(ElementKind::Hex20, [1, 1, 1], |p| [p[0], p[1], 1.0 + p[2]]);
    let low = q.block(ElementKind::Hex8, [2, 1, 1], |p| [2.0 * p[0], p[1], p[2]]);
    // A tet10 beside the hex20, sharing its edge x = 1, z = 1 (corners and mid-side node).
    // Positive volume: (1, 0, 1), (2, 0, 1), (1, 1, 1), (1.5, 0, 2).
    let corners = [[1.0, 0.0, 1.0], [2.0, 0.0, 1.0], [1.0, 1.0, 1.0], [1.5, 0.0, 2.0]].map(|x| q.node_at(x).unwrap_or_else(|| q.add_node(x)));
    let mut nodes = corners.to_vec();
    for e in ElementKind::Tet10.edges() {
        let mid = [0, 1, 2].map(|k| 0.5 * (q.xyz[&corners[e[0]]][k] + q.xyz[&corners[e[1]]][k]));
        nodes.push(q.node_at(mid).unwrap_or_else(|| q.add_node(mid)));
    }
    let shared = q.node_at([1.0, 0.5, 1.0]).unwrap();
    assert!(nodes.contains(&shared));
    q.next_element += 1;
    let tet = q.next_element;
    q.model.set_element(tet, ElementKind::Tet10, nodes).unwrap();
    q.define(1000.0, 0.3, NonlinearStatic::default());
    q.face_set("lower_top", &low, |x| near(x[2], 1.0));
    let secondary: Vec<u32> = up.iter().copied().chain([tet]).collect();
    q.face_set("upper_bottom", &secondary, |x| near(x[2], 1.0));
    q.pair("interface", "lower_top", "upper_bottom");
    q.bc("z0", [Some(0.0), Some(0.0), Some(0.0)], |x| near(x[2], 0.0));
    q.bc("push", [Some(0.0), Some(0.0), Some(-0.01)], |x| near(x[2], 2.0));
    let e = solve_nonlinear(&q.model, Some("Step"), &Limits::default()).unwrap_err();
    let named = |m: &str| m.contains(&format!("mid-side node {shared} ")) && m.contains("tri6") && m.contains("quad8");
    assert!(matches!(&e, Error::NotSupported(m) if named(m)), "{e}");
}

/// Hertz contact of a cylinder on a flat (plane strain): a quarter model of an indenter whose
/// bottom is a circle of radius `r` on a block, symmetric about x = 0, one element thick in z
/// with z held (plane strain), of `kind` elements with `1 / coarsen` the cells of the hex8 model
/// in x and y; hex20 and tet10 have their mid-side nodes on the circle (curved faces). Returns the
/// model and the result.
fn hertz_cylinder(kind: ElementKind, coarsen: usize, finite: bool, push: f64) -> (Blocks, StaticResult) {
    let r = 10.0;
    let mut b = Blocks::default();
    let grade = |s: f64| s.powf(1.6);
    let foundation = b.quadratic_block(kind, [24 / coarsen, 14 / coarsen, 1], true, |p| [5.0 * grade(p[0]), -5.0 * grade(1.0 - p[1]), 0.2 * p[2]]);
    let indenter = b.quadratic_block(kind, [20 / coarsen, 10 / coarsen, 1], true, |p| {
        let x = 2.5 * grade(p[0]);
        [x, r - (r * r - x * x).sqrt() + 2.5 * grade(p[1]), 0.2 * p[2]]
    });
    b.define(1000.0, 0.3, NonlinearStatic { nonlinear_geometry: finite, initial_increment: 0.25, ..NonlinearStatic::default() });
    b.face_set("flat", &foundation, |x| near(x[1], 0.0));
    let bottom = |x: [f64; 3]| (x[1] - (r - (r * r - x[0] * x[0]).sqrt())).abs() < 1e-9;
    b.face_set("round", &indenter, bottom);
    b.pair("hertz", "flat", "round");
    b.bc("symmetry", [Some(0.0), None, None], |x| near(x[0], 0.0));
    b.bc("plane", [None, None, Some(0.0)], |_| true);
    b.bc("base", [None, Some(0.0), None], |x| near(x[1], -5.0));
    let top = 2.5 + r - (r * r - 2.5 * 2.5f64).sqrt();
    b.bc("push", [None, Some(-push), None], move |x| {
        x[1] > 2.5 && (x[1] - (r - (r * r - x[0] * x[0]).sqrt()) - 2.5).abs() < 1e-9 && x[1] <= top + 1e-9
    });
    let r = solve_nonlinear(&b.model, Some("Step"), &Limits::default()).unwrap();
    report(&r);
    (b, r)
}

/// Hertz's closed form for the model of [`hertz_cylinder`]: the load per length P (the reaction
/// at the base over the slab's thickness, both halves), the half-width b and the peak pressure
/// p₀, with the largest nodal pressure found.
fn hertz_line(b: &Blocks, r: &StaticResult) -> (f64, f64, f64, f64) {
    let base: f64 = b.nodes(|x| near(x[1], -5.0)).iter().map(|l| r.reaction[r.nodes.binary_search(l).unwrap()][1]).sum();
    let p = 2.0 * base / 0.2;
    let e_star = 1000.0 / (2.0 * (1.0 - 0.09));
    let half_width = (4.0 * p * 10.0 / (std::f64::consts::PI * e_star)).sqrt();
    let peak = 2.0 * p / (std::f64::consts::PI * half_width);
    let found = r.contacts[0].pressure.iter().fold(0.0f64, |m, v| m.max(*v));
    (p, half_width, peak, found)
}

/// Hertz's line contact (H. Hertz, *Über die Berührung fester elastischer Körper*, J. reine
/// angew. Math. 92:156–171, 1882; K. L. Johnson, *Contact Mechanics*, Cambridge 1985, §4.3):
/// for a load P per length, b = √(4 P R / (π E*)) and p₀ = 2 P / (π b), E* = E / (2 (1 − ν²))
/// for two bodies of one material.
#[test]
fn hertz_cylinder_on_a_flat() {
    for finite in [false, true] {
        let (b, r) = hertz_cylinder(ElementKind::Hex8, 1, finite, 0.03);
        let (p, half_width, peak, found) = hertz_line(&b, &r);
        let (before, next) = edge(&b, &r.contacts[0], 2);
        eprintln!("finite {finite}: P {p:.4} b {half_width:.4} ({before:.4}..{next:.4}) p0 {peak:.4} found {found:.4}");
        assert!((found - peak).abs() < 0.03 * peak, "p0 {found} vs {peak}");
        assert!(before <= half_width && half_width <= next, "b {half_width} not in ({before}, {next})");
        if finite {
            quadratic(&r);
        }
    }
}

/// Hertz's line contact on quadratic faces (the closed form of [`hertz_cylinder_on_a_flat`]):
/// `kind` (hex20 or tet10) with half the cells in x and y (the hex8 model's nodes), their mid-side
/// nodes on the circle, so the secondary faces are curved and the facets of the segmentation are
/// not the faces. The nodal pressures alternate between corners and mid-side nodes by up to about
/// a tenth around the ellipse (the transformed basis spreads each node's force over its own
/// area), so the peak is the largest nodal value and the profile is checked within 15 % inside.
fn hertz_on_quadratic_faces(kind: ElementKind, peak_tol: f64) {
    let (b, r) = hertz_cylinder(kind, 2, false, 0.03);
    let c = &r.contacts[0];
    let (p, half_width, peak, found) = hertz_line(&b, &r);
    let (before, next) = edge(&b, c, 2);
    eprintln!("{kind:?}: P {p:.4} b {half_width:.4} ({before:.4}..{next:.4}) p0 {peak:.4} found {found:.4}");
    assert!((found - peak).abs() < peak_tol * peak, "{kind:?}: p0 {found} vs {peak}");
    assert!(before <= half_width && half_width <= next, "{kind:?}: b {half_width} not in ({before}, {next})");
    // The contact carries the whole load (P is per length of both halves).
    let base = p * 0.2 / 2.0;
    assert!((c.normal_force.abs() - base.abs()).abs() < 1e-3 * base.abs(), "{kind:?}: {} vs {base}", c.normal_force);
    for (l, v) in c.nodes.iter().zip(&c.pressure) {
        let x = b.xyz[l][0];
        if near(b.xyz[l][2], 0.0) && x < 0.6 * half_width {
            let hertz = peak * (1.0 - (x / half_width).powi(2)).sqrt();
            assert!((v - hertz).abs() < 0.15 * peak, "{kind:?}: p({x}) {v} vs {hertz}");
        }
    }
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn hertz_cylinder_on_hex20_faces() {
    hertz_on_quadratic_faces(ElementKind::Hex20, 0.02);
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn hertz_cylinder_on_tet10_faces() {
    hertz_on_quadratic_faces(ElementKind::Tet10, 0.05);
}

/// The contact's edge along the x axis (the nodes on the plane where coordinate `across` is 0):
/// the node before the last closed one and the first open one after it, so that the edge lies
/// within an element of the last closed node (whose status the augmentation's tolerance decides).
fn edge(b: &Blocks, c: &fieldcraft_solver::contact::ContactResult, across: usize) -> (f64, f64) {
    let mut along: Vec<(f64, Status)> =
        c.nodes.iter().zip(&c.status).filter(|(l, _)| near(b.xyz[l][across], 0.0)).map(|(l, s)| (b.xyz[l][0], *s)).collect();
    along.sort_by(|x, y| x.0.total_cmp(&y.0));
    let last = along.iter().rposition(|(_, s)| *s == Status::Closed).unwrap_or(0);
    let before = last.checked_sub(1).and_then(|i| along.get(i)).map_or(0.0, |p| p.0);
    let next = along.get(last + 1).map_or(f64::INFINITY, |p| p.0);
    (before, next)
}

/// Newton's rate with contact: within each solution of an increment (its first, and each
/// augmentation's), once the residual is small every correction squares it, up to a constant
/// (the consistent tangent of the mortar integrals, segmentation and normal included).
fn quadratic(r: &StaticResult) {
    let mut checked = 0;
    for inc in &r.nonlinear.as_ref().unwrap().increments {
        for w in inc.iterations.windows(2) {
            let (a, b) = (w[0], w[1]);
            // Below 1e-10 of the force scale the sums' round-off takes over.
            if a.augmentation != b.augmentation || a.force > 1e-3 || b.force < 1e-10 {
                continue;
            }
            assert!(b.force <= 100.0 * a.force * a.force, "increment {}: {:e} after {:e}", inc.number, b.force, a.force);
            checked += 1;
        }
    }
    assert!(checked >= 3, "too few iterations near the solution to judge ({checked})");
}

/// Two unit blocks, the upper one starting `overlap` into the lower one, pressed together.
fn stacked(overlap: f64, penetration: InitialPenetration) -> Blocks {
    let mut b = Blocks::default();
    let low = b.block(ElementKind::Hex8, [2, 2, 1], |p| p);
    let up = b.block(ElementKind::Hex8, [3, 3, 1], |p| [p[0], p[1], 1.0 - overlap + p[2]]);
    b.define(1000.0, 0.3, NonlinearStatic { nonlinear_geometry: true, initial_increment: 1.0, ..NonlinearStatic::default() });
    b.face_set("lower_top", &low, |x| near(x[2], 1.0));
    b.face_set("upper_bottom", &up, |x| near(x[2], 1.0 - overlap));
    let mut c = Contact::new(ContactSurfaces::Pair { main: "lower_top".into(), secondary: "upper_bottom".into() });
    c.initial_penetration = penetration;
    b.model.add_contact("interface", c).unwrap();
    b.bc("x0", [Some(0.0), None, None], |x| near(x[0], 0.0));
    b.bc("y0", [None, Some(0.0), None], |x| near(x[1], 0.0));
    b.bc("z0", [None, None, Some(0.0)], |x| near(x[2], 0.0));
    b.bc("push", [None, None, Some(-0.01)], move |x| near(x[2], 2.0 - overlap));
    b
}

#[test]
fn initial_overclosure_is_adjusted_or_refused() {
    let b = stacked(1e-3, InitialPenetration::Adjust);
    let r = solve_nonlinear(&b.model, Some("Step"), &Limits::default()).unwrap();
    assert!(r.warnings.iter().any(|w| w.contains("start inside the main surface")), "{:?}", r.warnings);
    // Adjusted, the blocks behave as if they started touching: the same stress as without the
    // overlap.
    let touching = solve_nonlinear(&stacked(0.0, InitialPenetration::Adjust).model, Some("Step"), &Limits::default()).unwrap();
    let (a, t) = (r.elements[0].nodal[0][2], touching.elements[0].nodal[0][2]);
    assert!((a - t).abs() < 1e-4 * t.abs(), "{a} vs {t}");
    let e = solve_nonlinear(&stacked(1e-3, InitialPenetration::Refuse).model, Some("Step"), &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("starts") && e.contains("inside"), "{e}");
}

#[test]
fn contact_results_do_not_depend_on_the_thread_count() {
    let b = stacked(0.0, InitialPenetration::Adjust);
    let run = |threads| fieldcraft_parallel::with_threads(threads, || solve_nonlinear(&b.model, Some("Step"), &Limits::default()).unwrap());
    let (one, many) = (run(1), run(4));
    assert_eq!(one.displacement, many.displacement);
    assert_eq!(one.contacts, many.contacts);
}

#[test]
fn what_implicit_contact_refuses() {
    let solve = |m: &Model| solve_nonlinear(m, Some("Step"), &Limits::default()).unwrap_err().to_string();
    // Linear analyses refuse every contact by name.
    let mut b = stacked(0.0, InitialPenetration::Adjust);
    b.model.add_step(Step::new("Linear", Procedure::Static {}), None).unwrap();
    let e = fieldcraft_solver::solve_static(&b.model, Some("Linear"), &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("interface") && e.contains("linear static"), "{e}");
    // Arc-length control.
    let mut b = stacked(0.0, InitialPenetration::Adjust);
    let arc = NonlinearStatic { arc_length: Some(ArcLength::default()), ..NonlinearStatic::default() };
    b.model.replace_step(Step { procedure: Procedure::NonlinearStatic(arc), ..b.model.step("Step").unwrap().clone() }).unwrap();
    assert!(solve(&b.model).contains("arc-length"));
}

/// Hertz's point contact (Johnson 1985, §4.2): a sphere of radius R on a flat of the same
/// material under a load P has the contact radius a = (3 P R / (4 E*))^(1/3) and the peak
/// pressure p₀ = 3 P / (2 π a²), E* = E / (2 (1 − ν²)). A quarter model, displacement controlled;
/// the load is the reaction at the base.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn hertz_sphere_on_a_flat() {
    let (r, e, nu) = (10.0, 1000.0, 0.3);
    let mut b = Blocks::default();
    let grade = |s: f64| s.powf(1.5);
    let foundation = b.block(ElementKind::Hex8, [16, 16, 10], |p| [5.0 * grade(p[0]), 5.0 * grade(p[1]), -5.0 * grade(1.0 - p[2])]);
    let cap = move |x: f64, y: f64| r - (r * r - x * x - y * y).sqrt();
    let indenter = b.block(ElementKind::Hex8, [14, 14, 8], |p| {
        let (x, y) = (2.5 * grade(p[0]), 2.5 * grade(p[1]));
        let bottom = cap(x, y);
        [x, y, bottom + (3.0 - bottom) * grade(p[2])]
    });
    b.define(e, nu, NonlinearStatic { nonlinear_geometry: false, initial_increment: 0.25, ..NonlinearStatic::default() });
    b.face_set("flat", &foundation, |x| near(x[2], 0.0));
    b.face_set("round", &indenter, move |x| (x[2] - cap(x[0], x[1])).abs() < 1e-9 && x[2] < 1.0);
    b.pair("hertz", "flat", "round");
    b.bc("x0", [Some(0.0), None, None], |x| near(x[0], 0.0));
    b.bc("y0", [None, Some(0.0), None], |x| near(x[1], 0.0));
    b.bc("base", [None, None, Some(0.0)], |x| near(x[2], -5.0));
    b.bc("push", [None, None, Some(-0.04)], |x| near(x[2], 3.0));
    let res = solve_nonlinear(&b.model, Some("Step"), &Limits::default()).unwrap();
    report(&res);
    let c = &res.contacts[0];
    let base: f64 = b.nodes(|x| near(x[2], -5.0)).iter().map(|l| res.reaction[res.nodes.binary_search(l).unwrap()][2]).sum();
    let p = 4.0 * base;
    let e_star = e / (2.0 * (1.0 - nu * nu));
    let a = (3.0 * p * r / (4.0 * e_star)).cbrt();
    let peak = 3.0 * p / (2.0 * std::f64::consts::PI * a * a);
    let found = c.pressure.iter().fold(0.0f64, |m, v| m.max(*v));
    let (before, next) = edge(&b, c, 1);
    eprintln!("P {p:.4} a {a:.4} ({before:.4}..{next:.4}) p0 {peak:.4} found {found:.4}");
    assert!((found - peak).abs() < 0.05 * peak, "p0 {found} vs {peak}");
    assert!(before <= a && a <= next, "a {a} not in ({before}, {next})");
}

/// Two blocks, the upper one pressed down by `down` and moved sideways by `side` at its top
/// (both in proportion over the step), with Coulomb friction `friction`; the reactions at the
/// top (normal, tangential) and the result.
fn sheared(down: f64, side: f64, friction: Friction, finite: bool) -> (f64, f64, StaticResult) {
    sheared_blocks(ElementKind::Hex8, down, side, friction, finite, false)
}

/// [`sheared`] with blocks of `kind`; `whole`: every face of the upper block is on the secondary
/// surface (its sides and top too), not only its bottom.
fn sheared_blocks(kind: ElementKind, down: f64, side: f64, friction: Friction, finite: bool, whole: bool) -> (f64, f64, StaticResult) {
    let mut b = Blocks::default();
    let low = b.block(kind, [3, 3, 1], |p| p);
    let up = b.block(kind, [2, 2, 1], |p| [p[0], p[1], 1.0 + p[2]]);
    b.define(1000.0, 0.3, NonlinearStatic { nonlinear_geometry: finite, initial_increment: 0.25, ..NonlinearStatic::default() });
    b.face_set("lower_top", &low, |x| near(x[2], 1.0));
    // The upper block's bottom, or its whole boundary (each side's faces).
    let mut faces = b.faces(&up, |x| near(x[2], 1.0));
    if whole {
        for (k, at) in [(0, 0.0), (0, 1.0), (1, 0.0), (1, 1.0), (2, 2.0)] {
            faces.extend(b.faces(&up, |x| near(x[k], at)));
        }
    }
    b.model.add_set("upper_bottom", Set { faces, ..Set::default() }).unwrap();
    let mut c = Contact::new(ContactSurfaces::Pair { main: "lower_top".into(), secondary: "upper_bottom".into() });
    c.friction = friction;
    b.model.add_contact("interface", c).unwrap();
    b.bc("base", [Some(0.0), Some(0.0), Some(0.0)], |x| near(x[2], 0.0));
    b.bc("top", [Some(side), Some(0.0), Some(-down)], |x| near(x[2], 2.0));
    let r = solve_nonlinear(&b.model, Some("Step"), &Limits::default()).unwrap();
    report(&r);
    let sum = |k: usize| -> f64 { b.nodes(|x| near(x[2], 2.0)).iter().map(|l| r.reaction[r.nodes.binary_search(l).unwrap()][k]).sum() };
    (-sum(2), sum(0), r)
}

/// Coulomb's law across a frictional interface: a tangential demand below μ N sticks (the blocks
/// shear as one, nothing slips), above it the interface slips and carries exactly μ N; with a
/// kinetic coefficient the slipping interface carries μ_k N.
#[test]
fn coulomb_friction_sticks_below_the_cone_and_slips_on_it() {
    for finite in [false, true] {
        let mu = 0.3;
        // Stuck as one body, the shear the top displacement asks for: a frictionless trial gives
        // the normal force, a stuck one with a very large μ the tangential force.
        let (n, t, _) = sheared(0.01, 0.001, Friction::coulomb(5.0), finite);
        let ratio = t / n;
        eprintln!("finite {finite}: stuck N {n:.4} T {t:.4} ratio {ratio:.4}");
        // Well below the cone (the shear concentrates at the edges of the interface, which slip
        // first): sticks, the same forces as glued.
        let below = 0.4 * mu / ratio;
        let (n1, t1, r1) = sheared(0.01, 0.001 * below, Friction::coulomb(mu), finite);
        assert!((t1 / n1 - 0.4 * mu).abs() < 1e-3 * mu, "finite {finite}: {} vs {}", t1 / n1, 0.4 * mu);
        // The middle of the interface sticks (its corners, where the shear concentrates, may slip).
        // (The upper block's 3 × 3 bottom nodes come first, in x then y order: the fifth.)
        assert_eq!(r1.contacts[0].status[4], Status::Stick, "finite {finite}: {:?}", r1.contacts[0].status);
        // Above it: slips, carrying μ N.
        let (n2, t2, r2) = sheared(0.01, 0.001 * 2.0 * mu / ratio, Friction::coulomb(mu), finite);
        // With finite sliding the upper block turns a little as it slides, and the contact acts
        // along its face's normal, a few thousandths of a radian from the vertical.
        let tol = if finite { 0.03 } else { 1e-3 };
        assert!((t2 - mu * n2).abs() < tol * mu * n2, "finite {finite}: T {t2} vs μN {}", mu * n2);
        assert!(r2.contacts[0].status.iter().all(|s| *s == Status::Slip), "finite {finite}: {:?}", r2.contacts[0].status);
        let shear: f64 = r2.contacts[0].shear.iter().map(|s| s[0]).zip(&r2.contacts[0].pressure).map(|(s, p)| s / p).sum::<f64>()
            / r2.contacts[0].nodes.len() as f64;
        // Node by node the nodal forces' ratio scatters a little about μ (the smoothing of the
        // nodal pressure and shear differs at the interface's edges); the mean is μ.
        assert!((shear + mu).abs() < 0.05 * mu, "finite {finite}: shear over pressure {shear}");
        // Kinetic friction once slipping.
        let decaying = Friction { static_coefficient: mu, kinetic: Some(0.2), decay: 1.0 };
        let (n3, t3, _) = sheared(0.01, 0.001 * 2.0 * mu / ratio, decaying, finite);
        assert!((t3 - 0.2 * n3).abs() < tol * 0.2 * n3, "finite {finite}: T {t3} vs μk N {}", 0.2 * n3);
    }
}

/// Coulomb friction on quadratic faces: hex20 blocks slipping carry μ N, every secondary node
/// (corners and mid-side nodes) slipping.
#[test]
fn quadratic_faces_slip_carrying_mu_n() {
    let mu = 0.3;
    let (n, t, r) = sheared_blocks(ElementKind::Hex20, 0.01, 0.03, Friction::coulomb(mu), false, false);
    eprintln!("hex20: N {n:.4} T {t:.4}");
    assert!((t - mu * n).abs() < 1e-3 * mu * n, "T {t} vs μN {}", mu * n);
    assert!(r.contacts[0].status.iter().all(|s| *s == Status::Slip), "{:?}", r.contacts[0].status);
}

/// The node normal at the corners of a frictional contact: with the upper block's sides on the
/// secondary surface too, its slipping bottom carries μ N, and every node's reported shear is μ
/// times its pressure, as the friction law's cone asks (with the sides' area vectors in the
/// normal the corners reported up to 25 μ p).
#[test]
fn uncontacted_sides_do_not_tilt_the_friction_plane() {
    let mu = 0.3;
    for finite in [false, true] {
        let (n, t, r) = sheared_blocks(ElementKind::Hex8, 0.01, 0.03, Friction::coulomb(mu), finite, true);
        let tol = if finite { 0.03 } else { 1e-3 };
        assert!((t - mu * n).abs() < tol * mu * n, "finite {finite}: T {t} vs μN {}", mu * n);
        let c = &r.contacts[0];
        for (((label, p), s), q) in c.nodes.iter().zip(&c.pressure).zip(&c.status).zip(&c.shear) {
            assert_eq!(*p > 0.0, *s == Status::Slip, "finite {finite}: node {label}: pressure {p}, {s:?}");
            if *p > 0.0 {
                let ratio = q.iter().map(|v| v * v).sum::<f64>().sqrt() / (mu * p);
                assert!((ratio - 1.0).abs() < 5e-3, "finite {finite}: node {label}: shear over μ p {ratio}");
            }
        }
    }
}

/// Two blocks of non-matching meshes tied at z = 1 (`lower_kind` 3 × 3 × 2, graded or distorted,
/// below; `upper_kind` 2 × 2 × 2 above), every outer boundary node moved by the linear field
/// u = G x: the tie patch test.
fn tied(lower_kind: ElementKind, upper_kind: ElementKind, distorted: bool, procedure: Procedure) -> (Blocks, [[f64; 3]; 3]) {
    let mut b = Blocks::default();
    let low = if distorted {
        b.block(lower_kind, [3, 3, 2], |p| [p[0] + 0.3 * p[0] * (1.0 - p[0]) * (p[1] - 0.5), p[1] + 0.3 * p[1] * (1.0 - p[1]) * (p[0] - 0.4), p[2]])
    } else {
        b.block(lower_kind, [3, 3, 2], |p| [p[0].powf(1.4), p[1].powf(0.8), p[2]])
    };
    let up = b.block(upper_kind, [2, 2, 2], |p| [p[0], p[1], 1.0 + p[2]]);
    let all: Vec<u32> = b.model.elements().keys().copied().collect();
    b.model.add_set("all", Set { elements: all, ..Set::default() }).unwrap();
    b.model.add_material("Mat", Material { elastic: Some(Elastic::Isotropic { young: 1000.0, poisson: 0.3 }), ..Material::default() }).unwrap();
    b.model
        .add_section("solid", Section { set: "all".into(), material: "Mat".into(), kind: SectionKind::Solid {}, formulation: Formulation::default() })
        .unwrap();
    b.model.add_step(Step::new("Step", procedure), None).unwrap();
    b.face_set("lower_top", &low, |x| near(x[2], 1.0));
    b.face_set("upper_bottom", &up, |x| near(x[2], 1.0));
    b.model.add_tie("glue", Tie::new("lower_top", "upper_bottom")).unwrap();
    let g = [[1e-3, 2e-4, -3e-4], [4e-4, -5e-4, 1e-4], [-2e-4, 3e-4, 2e-3]];
    let outer = |x: [f64; 3]| near(x[0], 0.0) || near(x[0], 1.0) || near(x[1], 0.0) || near(x[1], 1.0) || near(x[2], 0.0) || near(x[2], 2.0);
    for l in b.nodes(outer) {
        let x = b.xyz[&l];
        let u: [f64; 3] = std::array::from_fn(|i| (0..3).map(|j| g[i][j] * x[j]).sum());
        let name = format!("n{l}");
        b.model.add_set(&name, Set { nodes: vec![l], ..Set::default() }).unwrap();
        b.model
            .add_boundary_condition(
                &name,
                BoundaryCondition::Displacement { target: Target::Set(name.clone()), values: [Some(u[0]), Some(u[1]), Some(u[2]), None, None, None] },
            )
            .unwrap();
        b.model.set_boundary_condition_active("Step", &name, true).unwrap();
    }
    (b, g)
}

/// The tie patch test (dual mortar: Wohlmuth 2000): a linear displacement field across a tie
/// of non-matching meshes is reproduced exactly, so every element has the same stress, in linear
/// static, nonlinear static and through the interface nodes the tie moves.
#[test]
fn a_tie_passes_the_patch_test() {
    for (lower, upper, distorted, name) in [
        (ElementKind::Hex8, ElementKind::Hex8, false, "hex8 on graded hex8"),
        (ElementKind::Hex8, ElementKind::Tet4, false, "tet4 on hex8"),
        (ElementKind::Tet4, ElementKind::Hex8, false, "hex8 on tet4"),
        (ElementKind::Hex8, ElementKind::Hex8, true, "hex8 on distorted hex8"),
    ] {
        let (b, g) = tied(lower, upper, distorted, Procedure::Static {});
        let r = fieldcraft_solver::solve_static(&b.model, Some("Step"), &Limits::default()).unwrap();
        // The interface nodes follow the field too.
        let tol = if distorted { 1e-6 } else { 1e-11 };
        for l in b.nodes(|x| near(x[2], 1.0)) {
            let x = b.xyz[&l];
            let u = r.displacement[r.nodes.binary_search(&l).unwrap()];
            for i in 0..3 {
                let want: f64 = (0..3).map(|j| g[i][j] * x[j]).sum();
                assert!((u[i] - want).abs() < tol * 1e-2, "{name}: node {l} u {u:?}");
            }
        }
        let first = r.elements[0].nodal[0];
        for e in &r.elements {
            for s in &e.nodal {
                for k in 0..6 {
                    assert!((s[k] - first[k]).abs() < tol * 10.0, "{name}: element {} {s:?} vs {first:?}", e.label);
                }
            }
        }
        assert!(r.warnings.iter().all(|w| !w.contains("stay free")), "{name}: {:?}", r.warnings);
    }
    // The same in a nonlinear step with small strain.
    let (b, _) = tied(
        ElementKind::Hex8,
        ElementKind::Hex8,
        false,
        Procedure::NonlinearStatic(NonlinearStatic { nonlinear_geometry: false, ..NonlinearStatic::default() }),
    );
    let r = solve_nonlinear(&b.model, Some("Step"), &Limits::default()).unwrap();
    let first = r.elements[0].nodal[0];
    assert!(r.elements.iter().flat_map(|e| &e.nodal).all(|s| (0..6).all(|k| (s[k] - first[k]).abs() < 1e-8)));
}

#[test]
fn what_ties_refuse() {
    let (mut b, _) = tied(ElementKind::Hex8, ElementKind::Hex8, false, Procedure::Static {});
    // (Touching surfaces are tied within any tolerance.)
    b.model.replace_tie("glue", Tie { position_tolerance: Some(1e-9), ..Tie::new("lower_top", "upper_bottom") }).unwrap();
    let r = fieldcraft_solver::solve_static(&b.model, Some("Step"), &Limits::default()).unwrap();
    assert!(r.warnings.iter().all(|w| !w.contains("stay free")), "{:?}", r.warnings);
}

/// Two blocks tied at z = 1, conducting heat: `lower` (`lower_n` cells, graded unless `n`
/// divides evenly, conductivity `k_lower`) below and `upper` (2 × 2 × 2, `k_upper`) above, the
/// lower block secondary, in step "Step" of `procedure`. ρc = 1 and α = 1e-5 (reference 0) in
/// both.
fn tied_conductor(lower: ElementKind, lower_n: [usize; 3], graded: bool, upper: ElementKind, k: [f64; 2], procedure: Procedure) -> Blocks {
    let mut b = Blocks::default();
    let low = b.block(lower, lower_n, |p| if graded { [p[0].powf(1.4), p[1].powf(0.8), p[2]] } else { p });
    let up = b.block(upper, [2, 2, 2], |p| [p[0], p[1], 1.0 + p[2]]);
    for (name, elements, k) in [("low", &low, k[0]), ("up", &up, k[1])] {
        b.model.add_set(name, Set { elements: elements.clone(), ..Set::default() }).unwrap();
        let thermal = Thermal {
            conductivity: Some(Directional::Isotropic(k)),
            specific_heat: Some(1.0),
            expansion: Some(Directional::Isotropic(1e-5)),
            // Temperatures here are above a reference of 0 (the default is 0 °C, 273.15 K).
            reference_temperature: 0.0,
        };
        let material = Material {
            elastic: Some(Elastic::Isotropic { young: 1000.0, poisson: 0.3 }),
            density: Some(1.0),
            thermal: Some(thermal),
            ..Material::default()
        };
        b.model.add_material(name, material).unwrap();
        b.model
            .add_section(name, Section { set: name.into(), material: name.into(), kind: SectionKind::Solid {}, formulation: Formulation::default() })
            .unwrap();
    }
    b.model.add_step(Step::new("Step", procedure), None).unwrap();
    b.face_set("lower_top", &low, |x| near(x[2], 1.0));
    b.face_set("upper_bottom", &up, |x| near(x[2], 1.0));
    b.model.add_tie("glue", Tie::new("upper_bottom", "lower_top")).unwrap();
    b
}

/// Hold the nodes satisfying `f` at temperature `value` in step "Step".
fn hold_temperature(b: &mut Blocks, name: &str, value: f64, f: impl Fn([f64; 3]) -> bool) {
    b.node_set(name, f);
    b.model.add_boundary_condition(name, BoundaryCondition::Temperature { target: Target::Set(name.into()), value }).unwrap();
    b.model.set_boundary_condition_active("Step", name, true).unwrap();
}

fn temperature(r: &fieldcraft_solver::HeatResult, label: u32) -> f64 {
    r.temperature()[r.nodes.binary_search(&label).unwrap()]
}

/// Steady conduction through two tied blocks of non-matching meshes, 20 at z = 0 and 100 at
/// z = 2: with one conductivity the profile is the linear T = 20 + 40 z everywhere, the interface
/// nodes included, the flux −k ∂T/∂z in every element and the heat flow k A ΔT / L in at the hot
/// end and out at the cold one (A = 1); with k₁ below and k₂ above, the series-resistance
/// solution, linear in each block with the kink at the interface temperature
/// T_i = (k₁ T₀ + k₂ T₁)/(k₁ + k₂) and the heat flow ΔT/(1/k₁ + 1/k₂).
#[test]
fn heat_crosses_a_tie_exactly() {
    for (lower, upper, n, graded, name) in [
        (ElementKind::Hex8, ElementKind::Hex8, [2, 2, 2], false, "matching hex8"),
        (ElementKind::Hex8, ElementKind::Hex8, [3, 3, 2], true, "hex8 on graded hex8"),
        (ElementKind::Tet4, ElementKind::Hex8, [3, 3, 2], false, "hex8 on tet4"),
        (ElementKind::Hex8, ElementKind::Tet4, [3, 3, 2], true, "tet4 on graded hex8"),
        (ElementKind::Tet4, ElementKind::Tet4, [3, 3, 2], false, "tet4 on tet4"),
    ] {
        for k in [[45.0, 45.0], [10.0, 40.0]] {
            let mut b = tied_conductor(lower, n, graded, upper, k, Procedure::HeatTransfer { transient: None });
            hold_temperature(&mut b, "cold", 20.0, |x| near(x[2], 0.0));
            hold_temperature(&mut b, "hot", 100.0, |x| near(x[2], 2.0));
            let r = fieldcraft_solver::solve_heat(&b.model, Some("Step"), &Limits::default()).unwrap();
            assert!(r.warnings.iter().all(|w| !w.contains("stay free")), "{name}: {:?}", r.warnings);
            let interface = (k[0] * 20.0 + k[1] * 100.0) / (k[0] + k[1]);
            let q = 80.0 / (1.0 / k[0] + 1.0 / k[1]);
            let want = |z: f64| if z <= 1.0 { 20.0 + (interface - 20.0) * z } else { interface + (100.0 - interface) * (z - 1.0) };
            for (l, x) in &b.xyz {
                let t = temperature(&r, *l);
                assert!((t - want(x[2])).abs() < 1e-9 * 100.0, "{name} {k:?}: node {l} at {x:?}: {t} vs {}", want(x[2]));
            }
            for e in &r.elements {
                for p in &e.points {
                    assert!(
                        (p.flux[2] + q).abs() < 1e-9 * q && p.flux[0].abs() < 1e-9 * q && p.flux[1].abs() < 1e-9 * q,
                        "{name} {k:?}: element {} {:?}",
                        e.label,
                        p.flux
                    );
                }
            }
            // The heat flow: in at the hot end, out at the cold one, nothing lost at the tie.
            let flow = |f: &dyn Fn([f64; 3]) -> bool| -> f64 { b.nodes(f).iter().map(|l| r.reaction[r.nodes.binary_search(l).unwrap()]).sum() };
            assert!((flow(&|x| near(x[2], 2.0)) - q).abs() < 1e-9 * q, "{name} {k:?}: {}", flow(&|x| near(x[2], 2.0)));
            assert!((flow(&|x| near(x[2], 0.0)) + q).abs() < 1e-9 * q, "{name} {k:?}: {}", flow(&|x| near(x[2], 0.0)));
            assert!(r.reaction_total.abs() < 1e-9 * q && r.applied_total.abs() < 1e-9 * q, "{name} {k:?}");
        }
    }
}

/// Transient conduction into a bar of two tied blocks (backward Euler, 0 initially, the top held
/// at 100 from the start): the single-body bar's history to round-off, with matching meshes and
/// with a graded 3 × 3 lower block (the field varies along the bar alone, which both hexahedral
/// meshes represent alike, so the tie is the only difference).
#[test]
fn transient_heat_through_a_tie_is_the_single_body_s() {
    let tr = Procedure::HeatTransfer { transient: Some(Transient { theta: 1.0, initial_temperature: 0.0, ..Transient::new(1.0, 0.05) }) };
    // The single body: one 2 × 2 × 4 block of the same material.
    let mut single = Blocks::default();
    let all = single.block(ElementKind::Hex8, [2, 2, 4], |p| [p[0], p[1], 2.0 * p[2]]);
    single.model.add_set("all", Set { elements: all, ..Set::default() }).unwrap();
    let thermal = Thermal { conductivity: Some(Directional::Isotropic(1.0)), specific_heat: Some(1.0), ..Thermal::default() };
    single.model.add_material("Mat", Material { density: Some(1.0), thermal: Some(thermal), ..Material::default() }).unwrap();
    single
        .model
        .add_section("solid", Section { set: "all".into(), material: "Mat".into(), kind: SectionKind::Solid {}, formulation: Formulation::default() })
        .unwrap();
    single.model.add_step(Step::new("Step", tr.clone()), None).unwrap();
    hold_temperature(&mut single, "hot", 100.0, |x| near(x[2], 2.0));
    let reference = fieldcraft_solver::solve_heat(&single.model, Some("Step"), &Limits::default()).unwrap();
    let at = |r: &fieldcraft_solver::HeatResult, xyz: &BTreeMap<u32, [f64; 3]>, z: f64| -> f64 {
        let nodes: Vec<u32> = xyz.iter().filter(|(_, x)| near(x[2], z)).map(|(l, _)| *l).collect();
        nodes.iter().map(|l| temperature(r, *l)).sum::<f64>() / nodes.len() as f64
    };
    for (n, graded) in [([2, 2, 2], false), ([3, 3, 2], true)] {
        let mut b = tied_conductor(ElementKind::Hex8, n, graded, ElementKind::Hex8, [1.0, 1.0], tr.clone());
        hold_temperature(&mut b, "hot", 100.0, |x| near(x[2], 2.0));
        let r = fieldcraft_solver::solve_heat(&b.model, Some("Step"), &Limits::default()).unwrap();
        assert_eq!(r.times.len(), reference.times.len());
        for z in [0.0, 0.5, 1.0, 1.5] {
            let (a, want) = (at(&r, &b.xyz, z), at(&reference, &single.xyz, z));
            assert!(want > 1.0, "z {z}: the heat reaches it ({want})");
            assert!((a - want).abs() < 1e-9 * want, "{n:?} z {z}: {a} vs {want}");
        }
    }
}

/// A coupled temperature–displacement step with a tie, held at a uniform 100 above the reference
/// 0: the temperature is uniform through the tie and both blocks expand freely, u = α ΔT x, with
/// no stress.
#[test]
fn a_tied_body_expands_freely_under_a_uniform_temperature() {
    let mut b = tied_conductor(ElementKind::Hex8, [3, 3, 2], true, ElementKind::Tet4, [10.0, 40.0], Procedure::CoupledTemperatureDisplacement {});
    hold_temperature(&mut b, "cold", 100.0, |x| near(x[2], 0.0));
    hold_temperature(&mut b, "hot", 100.0, |x| near(x[2], 2.0));
    // Statically determinate supports: the origin fixed, (1, 0, 0) held in y and z, (0, 1, 0) in z.
    b.bc("o", [Some(0.0), Some(0.0), Some(0.0)], |x| near(x[0], 0.0) && near(x[1], 0.0) && near(x[2], 0.0));
    b.bc("a", [None, Some(0.0), Some(0.0)], |x| near(x[0], 1.0) && near(x[1], 0.0) && near(x[2], 0.0));
    b.bc("c", [None, None, Some(0.0)], |x| near(x[0], 0.0) && near(x[1], 1.0) && near(x[2], 0.0));
    let r = fieldcraft_solver::solve_static(&b.model, Some("Step"), &Limits::default()).unwrap();
    let temps = r.temperature.clone().unwrap();
    assert!(temps.iter().all(|t| (t - 100.0).abs() < 1e-9), "{temps:?}");
    let strain = 1e-5 * 100.0;
    for (i, l) in r.nodes.iter().enumerate() {
        let x = b.xyz[l];
        for (u, x) in r.displacement[i].iter().zip(x) {
            assert!((u - strain * x).abs() < 1e-9 * strain, "node {l}: {:?}", r.displacement[i]);
        }
    }
    assert!(r.elements.iter().flat_map(|e| &e.nodal).flatten().all(|s| s.abs() < 1e-6), "stress-free");
}

/// A C-shaped ring of rubber (an annulus of radii 1 and 1.5 with a 10° opening at +x, one hex8
/// layer thick in z with z held: plane strain) whose two ends are pressed together through the
/// outer nodes two elements from each end, against general (single-surface) contact on every
/// exterior face of the body: Mooney–Rivlin, nearly incompressible, hybrid elements.
fn c_ring(finite: bool, closing: f64) -> (Blocks, StaticResult) {
    let mut b = Blocks::default();
    let (n, open) = (36usize, 5f64.to_radians());
    let angle = move |i: f64| open + i * (2.0 * std::f64::consts::PI - 2.0 * open);
    b.block(ElementKind::Hex8, [n, 2, 1], |p| {
        let (t, r) = (angle(p[0]), 1.5 - 0.5 * p[1]);
        [r * t.cos(), r * t.sin(), 0.2 * p[2]]
    });
    let all: Vec<u32> = b.model.elements().keys().copied().collect();
    b.model.add_set("all", Set { elements: all, ..Set::default() }).unwrap();
    let rubber = Material {
        elastic: Some(Elastic::Isotropic { young: 3.0, poisson: 0.499 }),
        hyperelastic: Some(Hyperelastic::MooneyRivlin { c01: 0.1 }),
        ..Material::default()
    };
    b.model.add_material("Rubber", rubber).unwrap();
    b.model
        .add_section(
            "solid",
            Section {
                set: "all".into(),
                material: "Rubber".into(),
                kind: SectionKind::Solid {},
                formulation: Formulation { hybrid: true, ..Formulation::default() },
            },
        )
        .unwrap();
    let settings = NonlinearStatic { nonlinear_geometry: finite, initial_increment: 0.1, max_increments: 200, ..NonlinearStatic::default() };
    b.model.add_step(Step::new("Step", Procedure::NonlinearStatic(settings)), None).unwrap();
    b.model.add_contact("self", Contact::new(ContactSurfaces::General { surface: None })).unwrap();
    // The outer nodes two elements from each end close towards the axis of symmetry y = 0; the
    // back of the ring is held in its plane; every node in z.
    let step = 2.0 * std::f64::consts::PI - 2.0 * open;
    let at = |x: [f64; 3], i: f64| {
        ((x[0] * x[0] + x[1] * x[1]).sqrt() - 1.5).abs() < 1e-9
            && (x[1].atan2(x[0]).rem_euclid(2.0 * std::f64::consts::PI) - angle(i)).abs() < 1e-6 * step
    };
    let near_end = 2.0 / n as f64;
    b.bc("upper", [None, Some(-closing), None], move |x| at(x, near_end));
    b.bc("lower", [None, Some(closing), None], move |x| at(x, 1.0 - near_end));
    b.bc("back", [Some(0.0), Some(0.0), None], move |x| at(x, 0.5));
    b.bc("plane", [None, None, Some(0.0)], |_| true);
    let r = solve_nonlinear(&b.model, Some("Step"), &Limits::default()).unwrap();
    report(&r);
    (b, r)
}

/// Self-contact: the ring's two ends, faces of one body and one general contact, meet and press
/// on each other without passing through, and the result keeps the model's mirror symmetry
/// about y = 0.
#[test]
fn a_rubber_ring_closes_on_itself_without_penetration() {
    for finite in [true] {
        let (b, r) = c_ring(finite, 0.15);
        let c = &r.contacts[0];
        assert!(c.closed() > 0, "finite {finite}: the ends touch");
        // No penetration beyond the augmented Lagrangian's tolerance.
        assert!(c.penetration < 1e-4, "finite {finite}: {}", c.penetration);
        let deformed = |l: u32| -> [f64; 3] {
            let u = r.displacement[r.nodes.binary_search(&l).unwrap()];
            let x = b.xyz[&l];
            [x[0] + u[0], x[1] + u[1], x[2] + u[2]]
        };
        // The ends' nodes stay on their own side of the plane of symmetry, up to the overlap the
        // mortar's weak (weighted) enforcement allows at the edge of a contact zone: a hundredth
        // of the elements' size (0.25).
        for (l, x) in &b.xyz {
            let y = x[1];
            if x[0] > 0.9 && y.abs() < 0.2 {
                let d = deformed(*l);
                assert!(d[1] * y.signum() >= -2.5e-3, "finite {finite}: node {l} crossed: {y} → {}", d[1]);
            }
        }
        // Mirror symmetry: the node at (x, −y) moves as the mirror of the node at (x, y).
        let mut checked = 0;
        for (l, x) in &b.xyz {
            let Some((m, _)) = b.xyz.iter().find(|(_, y)| (y[0] - x[0]).abs() < 1e-9 && (y[1] + x[1]).abs() < 1e-9 && (y[2] - x[2]).abs() < 1e-9)
            else {
                continue;
            };
            let (u, v) = (r.displacement[r.nodes.binary_search(l).unwrap()], r.displacement[r.nodes.binary_search(m).unwrap()]);
            assert!((u[0] - v[0]).abs() < 1e-6 && (u[1] + v[1]).abs() < 1e-6, "finite {finite}: nodes {l}, {m}: {u:?} {v:?}");
            checked += 1;
        }
        assert!(checked > 50, "{checked}");
        eprintln!("finite {finite}: closed {} penetration {:.2e} force {:?}", c.closed(), c.penetration, c.force);
    }
}

/// A rubber bellows (a quarter of a corrugated tube: three convolutions, mean radius 2, fold
/// depth 0.5, wall 0.15, height 3) squeezed axially by 0.8 with general contact on all its faces:
/// its folds close and press on one another, as in a boot seal. (Squeezed to half its height
/// the step stalls at 57%: deeply folded rubber needs better increment control; see M7.)
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn a_rubber_bellows_folds_onto_itself() {
    let mut b = Blocks::default();
    let (folds, height, depth, wall) = (3.0, 3.0, 0.5, 0.15);
    let quarter = std::f64::consts::FRAC_PI_2;
    b.block(ElementKind::Hex8, [4, 48, 1], move |p| {
        let theta = quarter * p[0];
        let s = p[1] * folds;
        let tri = 1.0 - 2.0 * (s - s.floor() - 0.5).abs(); // 0 → 1 → 0 per convolution
        let r = 2.0 + depth * tri + wall * p[2];
        [r * theta.cos(), r * theta.sin(), height * p[1]]
    });
    let all: Vec<u32> = b.model.elements().keys().copied().collect();
    b.model.add_set("all", Set { elements: all, ..Set::default() }).unwrap();
    let rubber = Material {
        elastic: Some(Elastic::Isotropic { young: 3.0, poisson: 0.499 }),
        hyperelastic: Some(Hyperelastic::MooneyRivlin { c01: 0.1 }),
        ..Material::default()
    };
    b.model.add_material("Rubber", rubber).unwrap();
    b.model
        .add_section(
            "solid",
            Section {
                set: "all".into(),
                material: "Rubber".into(),
                kind: SectionKind::Solid {},
                formulation: Formulation { hybrid: true, ..Formulation::default() },
            },
        )
        .unwrap();
    let settings = NonlinearStatic { initial_increment: 0.05, max_increment: 0.1, max_increments: 400, ..NonlinearStatic::default() };
    b.model.add_step(Step::new("Step", Procedure::NonlinearStatic(settings)), None).unwrap();
    b.model.add_contact("folds", Contact::new(ContactSurfaces::General { surface: None })).unwrap();
    b.bc("x0", [Some(0.0), None, None], |x| near(x[0], 0.0));
    b.bc("y0", [None, Some(0.0), None], |x| near(x[1], 0.0));
    b.bc("base", [None, None, Some(0.0)], |x| near(x[2], 0.0));
    b.bc("squeeze", [None, None, Some(-0.8)], move |x| near(x[2], height));
    let r = solve_nonlinear(&b.model, Some("Step"), &Limits::default()).unwrap();
    report(&r);
    let c = &r.contacts[0];
    eprintln!("closed {} penetration {:.2e} augmentations {}", c.closed(), c.penetration, c.augmentations);
    assert!(c.closed() > 0, "the folds touch");
    assert!(c.penetration < 1e-3 * wall, "{}", c.penetration);
}

/// The C-ring of [`c_ring`] as a plane strain model (quad4, hybrid, 0.2 thick: the hex8 layer's
/// thickness) with general contact on every exterior edge.
fn planar_c_ring(closing: f64, threads: usize) -> (Blocks, StaticResult) {
    let mut b = Blocks::default();
    let (n, open) = (36usize, 5f64.to_radians());
    let angle = move |i: f64| open + i * (2.0 * std::f64::consts::PI - 2.0 * open);
    let label = |i: usize, j: usize| 1 + (i + (n + 1) * j) as u32;
    for j in 0..=2 {
        for i in 0..=n {
            let (t, r) = (angle(i as f64 / n as f64), 1.5 - 0.25 * j as f64);
            let x = [r * t.cos(), r * t.sin(), 0.0];
            b.model.set_node(label(i, j), x).unwrap();
            b.xyz.insert(label(i, j), x);
        }
    }
    let mut next = 0;
    for j in 0..2 {
        for i in 0..n {
            next += 1;
            // Counter-clockwise in the plane: along the outer arc, then inwards and back.
            b.model.set_element(next, ElementKind::Quad4, vec![label(i, j), label(i + 1, j), label(i + 1, j + 1), label(i, j + 1)]).unwrap();
        }
    }
    let all: Vec<u32> = b.model.elements().keys().copied().collect();
    b.model.add_set("all", Set { elements: all, ..Set::default() }).unwrap();
    let rubber = Material {
        elastic: Some(Elastic::Isotropic { young: 3.0, poisson: 0.499 }),
        hyperelastic: Some(Hyperelastic::MooneyRivlin { c01: 0.1 }),
        ..Material::default()
    };
    b.model.add_material("Rubber", rubber).unwrap();
    let section = Section {
        set: "all".into(),
        material: "Rubber".into(),
        kind: SectionKind::PlaneStrain { thickness: 0.2 },
        formulation: Formulation { hybrid: true, ..Formulation::default() },
    };
    b.model.add_section("ring", section).unwrap();
    let settings = NonlinearStatic { initial_increment: 0.1, max_increments: 200, ..NonlinearStatic::default() };
    b.model.add_step(Step::new("Step", Procedure::NonlinearStatic(settings)), None).unwrap();
    b.model.add_contact("self", Contact::new(ContactSurfaces::General { surface: None })).unwrap();
    let step = 2.0 * std::f64::consts::PI - 2.0 * open;
    let at = |x: [f64; 3], i: f64| {
        ((x[0] * x[0] + x[1] * x[1]).sqrt() - 1.5).abs() < 1e-9
            && (x[1].atan2(x[0]).rem_euclid(2.0 * std::f64::consts::PI) - angle(i)).abs() < 1e-6 * step
    };
    let near_end = 2.0 / n as f64;
    b.bc("upper", [None, Some(-closing), None], move |x| at(x, near_end));
    b.bc("lower", [None, Some(closing), None], move |x| at(x, 1.0 - near_end));
    b.bc("back", [Some(0.0), Some(0.0), None], move |x| at(x, 0.5));
    let r = fieldcraft_parallel::with_threads(threads, || solve_nonlinear(&b.model, Some("Step"), &Limits::default()).unwrap());
    (b, r)
}

/// 2D self-contact: the plane strain ring closes on itself like the one-layer solid ring with
/// its z held (the same problem), the edges' mortar pressing as the faces' does; bit for bit on
/// any thread count.
#[test]
fn a_plane_strain_ring_closes_on_itself_as_the_solid_ring_does() {
    let (b2, r2) = planar_c_ring(0.15, 1);
    let (b3, r3) = c_ring(true, 0.15);
    let (c2, c3) = (&r2.contacts[0], &r3.contacts[0]);
    assert!(c2.closed() > 0 && c2.penetration < 1e-4, "{} {}", c2.closed(), c2.penetration);
    // Each edge node closes where its face nodes do (two per edge node, at z = 0 and 0.2).
    assert_eq!(2 * c2.closed(), c3.closed());
    // The same motion: every node of the plane against the solid's node at z = 0.
    let peak = c2.pressure.iter().copied().fold(0.0, f64::max);
    assert!(peak > 0.0);
    // Both independently solved models report the same contact, including where an end face
    // meets an uncontacted outer or Z face. Compare fields, not a frozen value from either
    // dimensional discretisation; their motions already agree to the tolerance below.
    assert!((c2.normal_force - c3.normal_force).abs() < 1e-5 * c2.normal_force.abs());
    let pressure =
        |contact: &fieldcraft_solver::contact::ContactResult, node| contact.nodes.binary_search(&node).ok().map_or(0.0, |i| contact.pressure[i]);
    let mut checked = 0;
    for (l, x) in &b2.xyz {
        let Some((m, _)) = b3.xyz.iter().find(|(_, y)| (y[0] - x[0]).abs() < 1e-9 && (y[1] - x[1]).abs() < 1e-9 && y[2] == 0.0) else { continue };
        let (u, v) = (r2.displacement[r2.nodes.binary_search(l).unwrap()], r3.displacement[r3.nodes.binary_search(m).unwrap()]);
        assert!((u[0] - v[0]).abs() < 1e-5 && (u[1] - v[1]).abs() < 1e-5, "node {l}: {u:?} {v:?}");
        for (solid, _) in b3.xyz.iter().filter(|(_, y)| (y[0] - x[0]).abs() < 1e-9 && (y[1] - x[1]).abs() < 1e-9) {
            assert!((pressure(c2, *l) - pressure(c3, *solid)).abs() < 1e-5 * peak, "pressure at plane node {l}, solid node {solid}");
        }
        checked += 1;
    }
    assert_eq!(checked, b2.xyz.len());
    let (_, many) = planar_c_ring(0.15, 4);
    assert_eq!(r2.displacement, many.displacement);
    assert_eq!(r2.contacts, many.contacts);
}

#[test]
fn what_2d_contact_refuses() {
    // Quadratic edges.
    let (mut b, _) = planar_c_ring(0.0, 1);
    let mut m = Model::default();
    for (i, x) in [[0., 0.], [1., 0.], [1., 1.], [0., 1.], [0.5, 0.], [1., 0.5], [0.5, 1.], [0., 0.5]].iter().enumerate() {
        m.set_node(i as u32 + 1, [x[0], x[1], 0.0]).unwrap();
    }
    m.set_element(1, ElementKind::Quad8, (1..=8).collect()).unwrap();
    m.add_set("all", Set { elements: vec![1], ..Set::default() }).unwrap();
    m.add_material("Mat", Material { elastic: Some(Elastic::Isotropic { young: 1.0, poisson: 0.3 }), ..Material::default() }).unwrap();
    m.add_section(
        "s",
        Section { set: "all".into(), material: "Mat".into(), kind: SectionKind::PlaneStrain { thickness: 1.0 }, formulation: Formulation::default() },
    )
    .unwrap();
    m.add_step(Step::new("Step", Procedure::NonlinearStatic(NonlinearStatic::default())), None).unwrap();
    m.add_contact("self", Contact::new(ContactSurfaces::General { surface: None })).unwrap();
    let e = solve_nonlinear(&m, Some("Step"), &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("edges of quad8") && e.contains("tri3 and quad4"), "{e}");
    // A plane strain ring with arc-length control.
    let arc = NonlinearStatic { arc_length: Some(ArcLength::default()), ..NonlinearStatic::default() };
    b.model.replace_step(Step { procedure: Procedure::NonlinearStatic(arc), ..b.model.step("Step").unwrap().clone() }).unwrap();
    let e = solve_nonlinear(&b.model, Some("Step"), &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("arc-length"), "{e}");
}
