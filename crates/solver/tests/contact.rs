//! Contact in nonlinear static steps against closed forms: the contact patch test on
//! non-matching meshes, Newton's convergence with contact, and what implicit contact refuses;
//! ties of non-matching meshes for displacement and for heat.
//! Run with `--nocapture` to see the iterations.

use std::collections::BTreeMap;

use fieldcraft_model::*;
use fieldcraft_solver::contact::Status;
use fieldcraft_solver::{Limits, StaticResult, solve_nonlinear};

/// A model built block by block.
#[derive(Default)]
struct Blocks {
    model: Model,
    xyz: BTreeMap<u32, [f64; 3]>,
    next_node: u32,
    next_element: u32,
}

impl Blocks {
    /// A structured block of `n` cells over the unit cube mapped by `map`, hex8 or tet4 (six
    /// Kuhn tetrahedra per cell); its elements' labels.
    fn block(&mut self, kind: ElementKind, n: [usize; 3], map: impl Fn([f64; 3]) -> [f64; 3]) -> Vec<u32> {
        let first = self.next_node + 1;
        let label = |i: [usize; 3]| first + (i[0] + (n[0] + 1) * (i[1] + (n[1] + 1) * i[2])) as u32;
        for k in 0..=n[2] {
            for j in 0..=n[1] {
                for i in 0..=n[0] {
                    let x = map([i as f64 / n[0] as f64, j as f64 / n[1] as f64, k as f64 / n[2] as f64]);
                    let l = label([i, j, k]);
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

/// A flat contact whose secondary surface also contains the cube's uncontacted sides. Those
/// sides must not tilt the normal used to report the pressure carried by the bottom face.
#[test]
fn uncontacted_corner_faces_do_not_tilt_reported_pressure() {
    let mut b = Blocks::default();
    let low = b.block(ElementKind::Hex8, [1, 1, 1], |p| p);
    let upper = b.block(ElementKind::Hex8, [1, 1, 1], |p| [p[0], p[1], p[2] + 1.0]);
    b.define(1000.0, 0.0, NonlinearStatic { nonlinear_geometry: false, initial_increment: 1.0, ..NonlinearStatic::default() });
    b.face_set("main", &low, |p| near(p[2], 1.0));
    b.face_set("secondary", &upper, |_| true);
    b.pair("Interface", "main", "secondary");
    for (name, nodes, values) in [
        ("base", (1..=8).collect(), [Some(0.0), Some(0.0), Some(0.0), None, None, None]),
        ("upper", (9..=16).collect(), [Some(0.0), Some(0.0), None, None, None, None]),
        ("top", (13..=16).collect(), [None, None, Some(-0.01), None, None, None]),
    ] {
        b.model.add_set(name, Set { nodes, ..Set::default() }).unwrap();
        b.model.add_boundary_condition(name, BoundaryCondition::Displacement { target: Target::Set(name.into()), values }).unwrap();
        b.model.set_boundary_condition_active("Step", name, true).unwrap();
    }
    let r = solve_nonlinear(&b.model, Some("Step"), &Limits::default()).unwrap();
    let contact = &r.contacts[0];
    let normal = contact.force[2];
    assert!(normal > 9.9 && normal <= 10.0, "physical compression force: {normal}");
    for (label, pressure) in contact.nodes.iter().zip(&contact.pressure) {
        if (9..=12).contains(label) {
            // Uniform pressure on a unit-area interface equals its actual assembled force;
            // this checks reporting independently of the unchanged augmentation stopping gap.
            assert!((pressure - normal).abs() < 1e-9, "bottom node {label}: pressure {pressure}, physical force {normal}");
        } else {
            assert_eq!(*pressure, 0.0);
        }
    }
    assert!((contact.normal_force - normal).abs() < 1e-10, "reported normal part must follow the contacted face");
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

/// Hertz contact of a cylinder on a flat (plane strain): a quarter model of an indenter whose
/// bottom is a circle of radius `r` on a block, symmetric about x = 0, one element thick in z
/// with z held (plane strain). Returns the model and the result.
fn hertz_cylinder(finite: bool, push: f64) -> (Blocks, StaticResult) {
    let r = 10.0;
    let mut b = Blocks::default();
    let grade = |s: f64| s.powf(1.6);
    let foundation = b.block(ElementKind::Hex8, [24, 14, 1], |p| [5.0 * grade(p[0]), -5.0 * grade(1.0 - p[1]), 0.2 * p[2]]);
    let indenter = b.block(ElementKind::Hex8, [20, 10, 1], |p| {
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

/// Hertz's line contact (H. Hertz, *Über die Berührung fester elastischer Körper*, J. reine
/// angew. Math. 92:156–171, 1882; K. L. Johnson, *Contact Mechanics*, Cambridge 1985, §4.3):
/// for a load P per length, b = √(4 P R / (π E*)) and p₀ = 2 P / (π b), E* = E / (2 (1 − ν²))
/// for two bodies of one material.
#[test]
fn hertz_cylinder_on_a_flat() {
    for finite in [false, true] {
        let (b, r) = hertz_cylinder(finite, 0.03);
        let c = &r.contacts[0];
        let thickness = 0.2;
        // The load per length: the reaction at the base over the slab's thickness, both halves.
        let base: f64 = b.nodes(|x| near(x[1], -5.0)).iter().map(|l| r.reaction[r.nodes.binary_search(l).unwrap()][1]).sum();
        let p = 2.0 * base / thickness;
        let e_star = 1000.0 / (2.0 * (1.0 - 0.09));
        let half_width = (4.0 * p * 10.0 / (std::f64::consts::PI * e_star)).sqrt();
        let peak = 2.0 * p / (std::f64::consts::PI * half_width);
        let found = c.pressure.iter().fold(0.0f64, |m, v| m.max(*v));
        let (before, next) = edge(&b, c, 2);
        eprintln!("finite {finite}: P {p:.4} b {half_width:.4} ({before:.4}..{next:.4}) p0 {peak:.4} found {found:.4}");
        assert!((found - peak).abs() < 0.03 * peak, "p0 {found} vs {peak}");
        assert!(before <= half_width && half_width <= next, "b {half_width} not in ({before}, {next})");
        if finite {
            quadratic(&r);
        }
    }
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
    let mut b = Blocks::default();
    let low = b.block(ElementKind::Hex8, [3, 3, 1], |p| p);
    let up = b.block(ElementKind::Hex8, [2, 2, 1], |p| [p[0], p[1], 1.0 + p[2]]);
    b.define(1000.0, 0.3, NonlinearStatic { nonlinear_geometry: finite, initial_increment: 0.25, ..NonlinearStatic::default() });
    b.face_set("lower_top", &low, |x| near(x[2], 1.0));
    b.face_set("upper_bottom", &up, |x| near(x[2], 1.0));
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
