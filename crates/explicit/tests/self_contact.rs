//! Self contact (single-surface contact) in explicit dynamics: general contact on a body's own
//! exterior faces, every face against every other one of the same surface. The contact never asks
//! which body a face belongs to, so a body that closes on itself must behave as two bodies that
//! collide, keep its momentum and its energy balance, and never feel a force from faces that only
//! bend or squeeze past each other.
//!
//! - Two bars joined by a soft hinge (one body) collide end to end like the same bars as two bodies
//!   with a contact pair, and both like one-dimensional wave theory (K. F. Graff, *Wave Motion in
//!   Elastic Solids*, Dover 1991, §2.1–2.3): contact for 2L/c at the force ρ c A v₀/2, velocities
//!   exchanged.
//! - A hairpin whose arms close on each other and slide with friction keeps its linear momentum to
//!   round-off and its angular momentum: the contact forces come in equal and opposite pairs
//!   along the main face's normal (T. Belytschko, W. K. Liu, B. Moran, K. I. Elkhodary,
//!   *Nonlinear Finite Elements for Continua and Structures*, 2nd ed., Wiley 2014, §10.3–10.4:
//!   the penalty contact forces of a pair of points act along one line), and the energy balance
//!   closes with the friction dissipated.
//! - A ring crushed flat between a rigid wall and a stiff platen: its inner surfaces meet and carry
//!   the load, with the penetration the penalty gives, p/ε (the penalty method of
//!   single-surface contact: D. J. Benson, J. O. Hallquist, *A single surface contact algorithm
//!   for the post-buckling analysis of shell structures*, CMAME 78:141–163, 1990).
//! - No false contact: a cantilever in large-amplitude bending, a thick tube with a small bore
//!   ovalising, and a layer of elements squeezed through its thickness move exactly as they do
//!   without contact.
//!
//! Tolerances (each is asserted in its test): the one-body and two-body bars agree within 1 % in
//! peak force and contact time and 0.005 v₀ in the rebound velocities; against wave theory the
//! contact lasts 2L/c within 2 %, the peak force is ρ c A v₀/2 within 5 % and the velocities are
//! exchanged within 0.03 v₀ (a penalty contact on 200 cells per bar is looser than the two-body
//! comparison); linear momentum drifts by less than 1e-10 of itself; the energy balance closes
//! within 0.005 (0.01 for the ring); the hairpin's angular momentum changes by less than 1e-6 of
//! itself frictionless and, with friction, by no more than the couple's bound δ_max μ ∫ F_n dt;
//! the crushed ring's walls cross by more than 0.2 t without contact and by less than 0.1 t with
//! it, and its inner surface carries the punch's force within 3 %; a body that only bends or
//! squeezes moves bit for bit as without contact.

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
    /// Density per element (for the lumped masses).
    density: BTreeMap<u32, f64>,
}

const QUANTUM: f64 = 1e-9;

impl Mesh {
    fn new() -> Mesh {
        Mesh { model: Model::default(), nodes: BTreeMap::new(), next_node: 1, next_element: 1, density: BTreeMap::new() }
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

    fn faces_where(&self, elements: &[u32], f: impl Fn(V3) -> bool) -> Vec<ElementFace> {
        let exterior = self.model.exterior_faces();
        let mut out = Vec::new();
        for &e in elements {
            let el = self.model.element(e).unwrap();
            for (i, local) in el.kind.faces().iter().enumerate() {
                let c = [0, 1, 2].map(|d| local.iter().map(|n| self.xyz(el.nodes[*n])[d]).sum::<f64>() / local.len() as f64);
                let face = ElementFace { element: e, face: i as u8 + 1 };
                if f(c) && exterior.binary_search(&face).is_ok() {
                    out.push(face);
                }
            }
        }
        out
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
        let rho = self.model.material(material).and_then(|m| m.density).unwrap();
        for e in elements {
            self.density.insert(*e, rho);
        }
        self.model.add_set(name, Set { elements: elements.to_vec(), ..Set::default() }).unwrap();
        let formulation = Formulation { integration: Integration::Reduced, ..Formulation::default() };
        self.model.add_section(name, Section { set: name.into(), material: material.into(), kind: SectionKind::Solid {}, formulation }).unwrap();
    }

    fn face_set(&mut self, name: &str, faces: Vec<ElementFace>) {
        assert!(!faces.is_empty(), "{name}");
        self.model.add_set(name, Set { faces, ..Set::default() }).unwrap();
    }

    fn node_set(&mut self, name: &str, nodes: Vec<u32>) {
        assert!(!nodes.is_empty(), "{name}");
        self.model.add_set(name, Set { nodes, ..Set::default() }).unwrap();
    }

    fn general(&mut self, name: &str, surface: Option<&str>, friction: Friction) {
        let c = Contact { friction, ..Contact::new(ContactSurfaces::General { surface: surface.map(String::from) }) };
        self.model.add_contact(name, c).unwrap();
    }

    fn pair(&mut self, name: &str, main: &str, secondary: &str, friction: Friction) {
        let c = Contact { friction, ..Contact::new(ContactSurfaces::Pair { main: main.into(), secondary: secondary.into() }) };
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

    /// The lumped (row-sum) mass of every node, as the solver lumps it.
    fn masses(&self) -> BTreeMap<u32, f64> {
        let mut out = BTreeMap::new();
        for (label, e) in self.model.elements() {
            let x: [V3; 8] = std::array::from_fn(|i| self.xyz(e.nodes[i]));
            let k = fieldcraft_fem::uniform::Hex8::hex8(&x).unwrap();
            for (n, v) in e.nodes.iter().zip(k.nodal_volume) {
                *out.entry(*n).or_insert(0.0) += self.density[label] * v;
            }
        }
        out
    }
}

fn held(set: &str, dofs: &[usize], value: f64) -> BoundaryCondition {
    let mut values = [None; 6];
    for k in dofs {
        values[*k] = Some(value);
    }
    BoundaryCondition::Displacement { target: Target::Set(set.into()), values }
}

fn index(r: &ExplicitResult, label: u32) -> usize {
    r.nodes.binary_search(&label).unwrap()
}

fn mean_velocity(r: &ExplicitResult, nodes: &[u32]) -> V3 {
    let f = r.frames.last().unwrap();
    let mut v = [0.0; 3];
    for n in nodes {
        let x = f.velocity[index(r, *n)];
        for k in 0..3 {
            v[k] += x[k] / nodes.len() as f64;
        }
    }
    v
}

fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn norm(a: V3) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

// ---------------------------------------------------------------------------------------------
// Two bars, one body or two.

const BAR_E: f64 = 1.0e4;
const BAR_RHO: f64 = 1.0;
const BAR_L: f64 = 10.0;
const BAR_W: f64 = 0.1;
const BAR_CELLS: usize = 200;

/// Two bars end to end, the first moving at v₀ into the second at rest, each with one hinge cell
/// on top of its end. `joined`: the two hinge cells share their top edge at the joint (one body,
/// general contact on all its exterior faces); otherwise they do not (two bodies, a contact pair
/// between the bars' ends). (mesh, bar A with its hinge cell, bar B with its hinge cell.)
fn hinged_bars(joined: bool, v0: f64, end: f64) -> (Mesh, Vec<u32>, Vec<u32>) {
    let mut m = Mesh::new();
    let dx = BAR_L / BAR_CELLS as f64;
    let (l, w) = (BAR_L, BAR_W);
    let mut a = m.block(0, [0.0; 3], [l, w, w], [BAR_CELLS, 1, 1], |x| x);
    let mut b = m.block(1, [l, 0.0, 0.0], [l, w, w], [BAR_CELLS, 1, 1], |x| x);
    // The hinge: a cell on each bar's end; their top edges at x = L are one line of nodes
    // (body tag 2 on both) when joined.
    let (pa, pb) = (2u8, if joined { 2u8 } else { 3u8 });
    let p = m.hex([
        (0, [l - dx, 0.0, w]),
        (0, [l, 0.0, w]),
        (0, [l, w, w]),
        (0, [l - dx, w, w]),
        (pa, [l - dx, 0.0, w + dx]),
        (pa, [l, 0.0, w + dx]),
        (pa, [l, w, w + dx]),
        (pa, [l - dx, w, w + dx]),
    ]);
    let q = m.hex([
        (1, [l, 0.0, w]),
        (1, [l + dx, 0.0, w]),
        (1, [l + dx, w, w]),
        (1, [l, w, w]),
        (pb, [l, 0.0, w + dx]),
        (3, [l + dx, 0.0, w + dx]),
        (3, [l + dx, w, w + dx]),
        (pb, [l, w, w + dx]),
    ]);
    m.material("Bar", BAR_E, 0.0, BAR_RHO, false);
    // Soft: its stiffness takes a ten-thousandth of the bars'; the same density, so its nodes
    // weigh what the bars' do per area and the penalty is the pair's.
    m.material("Hinge", BAR_E * 1e-4, 0.0, BAR_RHO, false);
    m.section("A", &a, "Bar");
    m.section("B", &b, "Bar");
    m.section("hinge", &[p, q], "Hinge");
    if joined {
        m.general("self", None, Friction::NONE);
    } else {
        let end_a = m.faces_where(&a, |c| (c[0] - l).abs() < 1e-9);
        let start_b = m.faces_where(&b, |c| (c[0] - l).abs() < 1e-9);
        m.face_set("end_a", end_a);
        m.face_set("start_b", start_b);
        m.pair("impact", "start_b", "end_a", Friction::NONE);
    }
    a.push(p);
    b.push(q);
    m.model.add_set("moving", Set { elements: a.clone(), ..Set::default() }).unwrap();
    let mut step = ExplicitDynamic::new(end);
    step.initial_velocities = vec![InitialVelocity { set: "moving".into(), velocity: [v0, 0.0, 0.0] }];
    step.bulk_viscosity = BulkViscosity { linear: 0.0, quadratic: 1.5 };
    step.history_interval = Some(end / 4000.0);
    m.step(step, &[]);
    (m, a, b)
}

/// The bars joined by a hinge (one body, general contact) collide as the same bars as two bodies
/// with a contact pair: the same peak force, contact time and rebound, and both as wave theory
/// says; momentum is kept to round-off and the energy balance closes.
#[test]
fn hinged_bars_collide_like_two_bodies() {
    let c = (BAR_E / BAR_RHO).sqrt();
    let v0 = 0.1;
    let duration = 2.0 * BAR_L / c;
    let plateau = BAR_RHO * c * BAR_W * BAR_W * v0 / 2.0;
    let mut out = Vec::new();
    for joined in [false, true] {
        let (m, a, b) = hinged_bars(joined, v0, 1.5 * duration);
        let r = m.solve();
        let s = &r.contacts[0];
        let lasted = s.last_contact.unwrap() - s.first_contact.unwrap();
        let va = mean_velocity(&r, &m.nodes_where(&a, |_| true))[0];
        let vb = mean_velocity(&r, &m.nodes_where(&b, |_| true))[0];
        let p0 = r.history.momentum[0][0];
        let drift = r.history.momentum.iter().map(|p| norm([p[0] - p0, p[1], p[2]])).fold(0.0, f64::max) / p0;
        let balance = r.energy.relative_error(r.initial_kinetic);
        eprintln!(
            "{}: contact {lasted:.5} vs 2L/c {duration:.5}; peak {:.5e} vs ρcAv₀/2 {plateau:.5e}; v_A {va:.5} v_B {vb:.5}; KE/KE₀ {:.5}; momentum drift {drift:.2e}; balance {balance:.2e}; penetration {:.3e}; {} increments",
            if joined { "one body, self contact" } else { "two bodies, pair" },
            s.peak_force,
            r.energy.kinetic / r.initial_kinetic,
            s.max_penetration,
            r.increments
        );
        assert!(((lasted - duration) / duration).abs() < 0.02, "{lasted} vs {duration}");
        assert!(((s.peak_force - plateau) / plateau).abs() < 0.05, "peak {} vs ρcAv₀/2 {plateau}", s.peak_force);
        assert!(va.abs() < 0.03 * v0 && (vb - v0).abs() < 0.03 * v0, "{va} {vb}");
        assert!(drift < 1e-10, "momentum drift {drift}");
        assert!(balance < 0.005, "{:?}", r.energy);
        out.push((s.peak_force, lasted, va, vb));
    }
    let (two, one) = (out[0], out[1]);
    eprintln!(
        "one body / two bodies: peak force {:+.3}%, contact time {:+.3}%, v_A {:+.2e} v₀, v_B {:+.2e} v₀",
        100.0 * (one.0 / two.0 - 1.0),
        100.0 * (one.1 / two.1 - 1.0),
        (one.2 - two.2) / v0,
        (one.3 - two.3) / v0
    );
    assert!((one.0 / two.0 - 1.0).abs() < 0.01, "peak {} vs {}", one.0, two.0);
    assert!((one.1 / two.1 - 1.0).abs() < 0.01, "time {} vs {}", one.1, two.1);
    assert!((one.2 - two.2).abs() < 0.005 * v0 && (one.3 - two.3).abs() < 0.005 * v0);
}

// ---------------------------------------------------------------------------------------------
// A hairpin closing on itself.

/// A U-shaped body: a spine x ∈ [−0.1, 0], z ∈ [−0.15, 0.15] and two arms x ∈ [0, 1] of thickness
/// 0.1 with a gap of 0.1 between them; the upper arm thrown down and along x, the lower arm up.
fn hairpin(friction: Friction, end: f64) -> (Mesh, Vec<u32>) {
    let mut m = Mesh::new();
    let mut all = m.block(0, [-0.1, 0.0, -0.15], [0.1, 0.1, 0.3], [2, 2, 6], |x| x);
    let upper = m.block(0, [0.0, 0.0, 0.05], [1.0, 0.1, 0.1], [20, 2, 2], |x| x);
    let lower = m.block(0, [0.0, 0.0, -0.15], [1.0, 0.1, 0.1], [20, 2, 2], |x| x);
    m.model.add_set("upper", Set { elements: upper.clone(), ..Set::default() }).unwrap();
    m.model.add_set("lower", Set { elements: lower.clone(), ..Set::default() }).unwrap();
    all.extend(&upper);
    all.extend(&lower);
    m.material("Mat", 1.0e4, 0.3, 1.0, false);
    m.section("all", &all, "Mat");
    m.general("self", None, friction);
    let mut step = ExplicitDynamic::new(end);
    step.initial_velocities =
        vec![InitialVelocity { set: "upper".into(), velocity: [0.5, 0.0, -2.0] }, InitialVelocity { set: "lower".into(), velocity: [0.0, 0.0, 2.0] }];
    step.output_interval = Some(end / 40.0);
    step.history_interval = Some(end / 2000.0);
    m.step(step, &[]);
    (m, all)
}

/// The total angular momentum about the origin, Σ m (X + u) × v, in each frame.
fn angular_momentum(m: &Mesh, r: &ExplicitResult) -> Vec<V3> {
    let masses = m.masses();
    r.frames
        .iter()
        .map(|f| {
            let mut l = [0.0; 3];
            for (i, n) in r.nodes.iter().enumerate() {
                let x0 = m.xyz(*n);
                let x = [0, 1, 2].map(|k| x0[k] + f.displacement[i][k]);
                let p = f.velocity[i].map(|v| masses[n] * v);
                let c = cross(x, p);
                for k in 0..3 {
                    l[k] += c[k];
                }
            }
            l
        })
        .collect()
}

/// How far the upper arm's lower surface has gone below the lower arm's upper surface, at its
/// deepest over the saved frames (negative while the arms are apart). The surfaces are the nodes
/// that start at z = ±0.05 (x ≥ 0), the lower one interpolated linearly along x at the same y.
fn arm_overlap(m: &Mesh, r: &ExplicitResult, all: &[u32]) -> f64 {
    let upper = m.nodes_where(all, |x| x[0] >= -1e-9 && (x[2] - 0.05).abs() < 1e-9);
    let lower = m.nodes_where(all, |x| x[0] >= -1e-9 && (x[2] + 0.05).abs() < 1e-9);
    let mut deepest = f64::NEG_INFINITY;
    for f in &r.frames {
        let at = |n: u32| {
            let (x0, i) = (m.xyz(n), index(r, n));
            [0, 1, 2].map(|k| x0[k] + f.displacement[i][k])
        };
        for &u in &upper {
            let pu = at(u);
            let y0 = m.xyz(u)[1];
            let mut line: Vec<V3> = lower.iter().filter(|n| (m.xyz(**n)[1] - y0).abs() < 1e-9).map(|n| at(*n)).collect();
            line.sort_by(|a, b| a[0].total_cmp(&b[0]));
            for w in line.windows(2) {
                if w[0][0] <= pu[0] && pu[0] <= w[1][0] && w[1][0] > w[0][0] {
                    let z = w[0][2] + (w[1][2] - w[0][2]) * (pu[0] - w[0][0]) / (w[1][0] - w[0][0]);
                    deepest = deepest.max(z - pu[2]);
                }
            }
        }
    }
    deepest
}

/// The arms of a hairpin hit each other: linear momentum kept to round-off; angular momentum kept
/// to round-off without friction (each normal force and its reaction act along one line, the
/// main face's normal through both points) and within the bound friction allows with it (the
/// tangential force at a node and its reaction at the point of the face under it are a couple of
/// arm at most the penetration δ, so |ΔL| ≤ δ_max μ ∫ F_n dt); the energy balance closed with
/// the friction dissipated; no pass-through (the upper arm's lower surface never goes below the
/// lower arm's upper surface by more than a tenth of an arm's thickness, where the same arms without
/// contact pass through each other by more than their own thickness); the same numbers on one thread and four.
#[test]
fn a_hairpin_closing_on_itself_keeps_momentum_and_energy() {
    let end = 0.12;
    let mu = 0.3;
    for friction in [Friction::NONE, Friction::coulomb(mu)] {
        let (m, all) = hairpin(friction, end);
        let r = m.solve();
        let crossing = arm_overlap(&m, &r, &all);
        eprintln!("hairpin μ = {}: the arms overlap by at most {crossing:.3e}", friction.static_coefficient);
        assert!(crossing < 0.01, "the arms cross: {crossing}");
        if friction.is_frictionless() {
            // The same arms without the contact do pass through each other.
            let (mut free, all) = hairpin(friction, end);
            free.model.remove_contact("self").unwrap();
            let through = arm_overlap(&free, &free.solve(), &all);
            eprintln!("hairpin without contact: the arms overlap by {through:.3e}");
            assert!(through > 0.1, "{through}");
        }
        let s = &r.contacts[0];
        let p0 = r.history.momentum[0];
        let drift = r.history.momentum.iter().map(|p| norm([p[0] - p0[0], p[1] - p0[1], p[2] - p0[2]])).fold(0.0, f64::max) / norm(p0);
        let l = angular_momentum(&m, &r);
        let change = l.iter().map(|x| norm([x[0] - l[0][0], x[1] - l[0][1], x[2] - l[0][2]])).fold(0.0, f64::max);
        let (t, f) = (&r.history.time, &r.history.contact_force);
        let impulse: f64 = t.windows(2).zip(f.windows(2)).map(|(t, f)| 0.5 * (f[0][0] + f[1][0]) * (t[1] - t[0])).sum();
        let bound = s.max_penetration * friction.static_coefficient * impulse;
        let balance = r.energy.relative_error(r.initial_kinetic);
        eprintln!(
            "hairpin μ = {}: contact {:.4}..{:.4}, peak {:.4e}; penetration {:.3e}; normal impulse {impulse:.4e}; friction {:.4e} of KE₀ {:.4e}; |Δp|/|p₀| {drift:.2e}; |ΔL| {change:.2e} = {:.2e} |L₀| (bound δ μ ∫F_n dt {bound:.2e}); balance {balance:.2e}; KE end / KE₀ {:.4}; {} increments",
            friction.static_coefficient,
            s.first_contact.unwrap_or(f64::NAN),
            s.last_contact.unwrap_or(f64::NAN),
            s.peak_force,
            s.max_penetration,
            r.energy.friction,
            r.initial_kinetic,
            change / norm(l[0]),
            r.energy.kinetic / r.initial_kinetic,
            r.increments
        );
        assert!(s.first_contact.is_some() && s.last_contact < Some(0.8 * end), "{s:?}");
        assert!(s.max_penetration < 0.1 * 0.05, "penetration {}", s.max_penetration);
        assert!(drift < 1e-10, "momentum drift {drift}");
        if friction.is_frictionless() {
            assert!(change < 1e-6 * norm(l[0]), "{change}");
        } else {
            assert!(r.energy.friction > 0.0);
            assert!(change <= bound, "{change} vs {bound}");
        }
        assert!(balance < 0.005, "{:?}", r.energy);
        if friction.is_frictionless() {
            continue;
        }
        // With friction (the most state): the same numbers on one thread and four.
        let one = fieldcraft_parallel::with_threads(1, || solve_explicit(&m.model, None, &Limits::default()).unwrap());
        let four = fieldcraft_parallel::with_threads(4, || solve_explicit(&m.model, None, &Limits::default()).unwrap());
        assert_eq!(one.frames, four.frames);
        assert_eq!(one.history, four.history);
        assert_eq!(one.contacts, four.contacts);
        assert_eq!(one.frames, r.frames);
    }
}

// ---------------------------------------------------------------------------------------------
// A ring crushed flat.

const RING_E: f64 = 1.0e4;
const RING_R: f64 = 1.0;
const RING_T: f64 = 0.1;
const RING_W: f64 = 0.1;

/// Half a ring (x ≥ 0, symmetric about x = 0) of mean radius 1 and wall 0.1, two neo-Hookean
/// elements through the wall, in plane strain, held at the bottom and pushed down `crush` over
/// `end` at the top, each over a punch x < 0.1 on its outer surface. `contact`: "none", "inner"
/// (general contact on the ring's inner surface) or "all" (general contact on every exterior face).
fn crushed_ring(contact: &str, cells: usize, crush: f64, end: f64) -> (Mesh, Vec<u32>) {
    use std::f64::consts::FRAC_PI_2;
    let mut m = Mesh::new();
    let (ri, ro) = (RING_R - RING_T / 2.0, RING_R + RING_T / 2.0);
    let ring = m.block(0, [-FRAC_PI_2, ri, 0.0], [2.0 * FRAC_PI_2, RING_T, RING_W], [cells, 2, 1], |p| [p[1] * p[0].cos(), p[2], p[1] * p[0].sin()]);
    m.material("Rubber", RING_E, 0.3, 1.0, true);
    m.section("ring", &ring, "Rubber");
    let radius = |x: V3| (x[0] * x[0] + x[2] * x[2]).sqrt();
    m.node_set("x0", m.nodes_where(&ring, |x| x[0].abs() < 1e-9));
    m.node_set("everything", m.nodes_where(&ring, |_| true));
    m.node_set("top", m.nodes_where(&ring, |x| (radius(x) - ro).abs() < 1e-9 && x[0] < 0.1 && x[2] > 0.0));
    m.node_set("bottom", m.nodes_where(&ring, |x| (radius(x) - ro).abs() < 1e-9 && x[0] < 0.1 && x[2] < 0.0));
    let inner = m.faces_where(&ring, |c| radius(c) < RING_R && c[1].abs() > 1e-9 && (c[1] - RING_W).abs() > 1e-9);
    m.face_set("inner", inner);
    match contact {
        "all" => m.general("self", None, Friction::NONE),
        "inner" => m.general("self", Some("inner"), Friction::NONE),
        _ => {}
    }
    let mut step = ExplicitDynamic::new(end);
    step.output_interval = Some(end / 20.0);
    step.history_interval = Some(end / 2000.0);
    m.step(
        step,
        &[
            ("sym", held("x0", &[0], 0.0)),
            ("plane", held("everything", &[1], 0.0)),
            ("support", held("bottom", &[2], 0.0)),
            ("press", held("top", &[2], -crush)),
        ],
    );
    (m, ring)
}

/// The inner surface's nodes in the last frame, upper (z > 0 at rest) or lower half, in the row
/// at y = `y`, as (x, z) sorted by x.
fn inner_profile(m: &Mesh, r: &ExplicitResult, ring: &[u32], upper: bool, y: f64) -> Vec<[f64; 2]> {
    let ri = RING_R - RING_T / 2.0;
    let f = r.frames.last().unwrap();
    let mut out: Vec<[f64; 2]> = m
        .nodes_where(ring, |x| ((x[0] * x[0] + x[2] * x[2]).sqrt() - ri).abs() < 1e-9 && (x[1] - y).abs() < 1e-9 && (x[2] > 0.0) == upper)
        .iter()
        .map(|n| {
            let (x0, u) = (m.xyz(*n), f.displacement[index(r, *n)]);
            [x0[0] + u[0], x0[2] + u[2]]
        })
        .collect();
    out.sort_by(|a, b| a[0].total_cmp(&b[0]));
    out
}

/// z of a polyline at x (None outside it).
fn height(line: &[[f64; 2]], x: f64) -> Option<f64> {
    line.windows(2)
        .find(|w| w[0][0] <= x && x <= w[1][0] && w[1][0] > w[0][0])
        .map(|w| w[0][1] + (x - w[0][0]) / (w[1][0] - w[0][0]) * (w[1][1] - w[0][1]))
}

/// How far each node of the upper inner surface (row y) is below the lower one (positive:
/// crossed), with its x.
fn crossings(m: &Mesh, r: &ExplicitResult, ring: &[u32], y: f64) -> Vec<[f64; 2]> {
    let (up, down) = (inner_profile(m, r, ring, true, y), inner_profile(m, r, ring, false, y));
    up.iter().filter_map(|p| height(&down, p[0]).map(|z| [p[0], z - p[1]])).collect()
}

/// ∫ max(d, 0) dx of a piecewise-linear d(x), cut exactly where it changes sign.
fn positive_integral(d: &[[f64; 2]]) -> f64 {
    d.windows(2)
        .map(|w| {
            let ([x0, a], [x1, b]) = (w[0], w[1]);
            match (a > 0.0, b > 0.0) {
                (true, true) => 0.5 * (a + b) * (x1 - x0),
                (true, false) => 0.5 * a * (x1 - x0) * a / (a - b),
                (false, true) => 0.5 * b * (x1 - x0) * b / (b - a),
                _ => 0.0,
            }
        })
        .sum()
}

/// The punch's mean force over the last `samples` history samples: the rate of the constraints'
/// work over the punch's speed (the supports and symmetry planes do no work).
fn punch_force(r: &ExplicitResult, crush: f64, end: f64, samples: usize) -> f64 {
    let (t, e) = (&r.history.time, &r.history.energies);
    let (k, n) = (t.len() - 1 - samples, t.len() - 1);
    (e[n].constraint - e[k].constraint) / (t[n] - t[k]) / (crush / end)
}

/// The mean of a contact's normal force over the last `samples` history samples.
fn mean_contact_force(r: &ExplicitResult, samples: usize) -> f64 {
    let f = &r.history.contact_force;
    f.iter().rev().take(samples).map(|x| x.first().copied().unwrap_or(0.0)).sum::<f64>() / samples as f64
}

/// min over the surface's nodes of m_i / A_i, the area A_i = κ_i max_d Σ_f a_f |n_f,d| with
/// κ_i = max_f Σ_d |n_f,d| (the contact module's penalty, recomputed from the mesh).
fn mass_per_area(m: &Mesh, faces: &[ElementFace]) -> f64 {
    let masses = m.masses();
    let mut directional: BTreeMap<u32, (V3, f64)> = BTreeMap::new();
    for f in faces {
        let nodes = m.model.face_nodes(*f).unwrap();
        let x: Vec<V3> = nodes.iter().map(|n| m.xyz(*n)).collect();
        let sub = |a: V3, b: V3| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let a = cross(sub(x[2], x[0]), sub(x[3], x[1])).map(|c| 0.5 * c);
        let kappa = a.iter().map(|c| c.abs()).sum::<f64>() / norm(a);
        for n in &nodes {
            let e = directional.entry(*n).or_insert(([0.0; 3], 1.0));
            for (sum, c) in e.0.iter_mut().zip(a) {
                *sum += c.abs() / 4.0;
            }
            e.1 = e.1.max(kappa);
        }
    }
    directional.iter().map(|(n, (d, k))| masses[n] / (k * d.iter().copied().fold(0.0, f64::max))).fold(f64::INFINITY, f64::min)
}

/// A ring crushed flat between two punches: without self contact its inner surfaces pass
/// through each other; with general contact (on the inner surface, or on every exterior face)
/// they meet and carry the punch's load, no deeper than the penalty's p/ε. (Three quasi-static
/// runs of 31 000 increments: the trunk gate has the hinged bars and the hairpin.)
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn a_ring_crushed_flat_carries_load_on_its_inner_surface() {
    let (cells, crush, end) = (48, 1.92, 8.0);
    let mut bending = f64::NAN;
    for kind in ["none", "inner", "all"] {
        let (m, ring) = crushed_ring(kind, cells, crush, end);
        let r = m.solve();
        let rows = [crossings(&m, &r, &ring, 0.0), crossings(&m, &r, &ring, RING_W)];
        let deepest = rows.iter().flatten().map(|c| c[1]).fold(f64::NEG_INFINITY, f64::max);
        let punch = punch_force(&r, crush, end, 10);
        let carried = mean_contact_force(&r, 10);
        let s = r.contacts.first().cloned().unwrap_or_default();
        let last = r.history.contact_force.last().and_then(|f| f.first().copied()).unwrap_or(0.0);
        let balance = r.energy.relative_error(r.initial_kinetic);
        eprintln!(
            "ring {kind}: deepest crossing {deepest:.4e}; punch force {punch:.4e}, contact force {carried:.4e} (means over the last 10 samples; last {last:.4e}, peak {:.4e}, first at {:?}, penetration {:.3e}); balance {balance:.2e}; hourglass {:.3e} of internal {:.3e}; {} increments",
            s.peak_force, s.first_contact, s.max_penetration, r.energy.hourglass, r.energy.internal, r.increments,
        );
        assert!(balance < 0.01, "{kind}: {:?}", r.energy);
        if kind == "none" {
            // The walls pass through each other; only the ring's bending resists.
            assert!(deepest > 0.2 * RING_T, "{deepest}");
            bending = punch;
            continue;
        }
        // The walls stop at each other, no deeper than a tenth of the wall, and the contact
        // carries the punch's load, many times what bending alone took.
        assert!(deepest < 0.1 * RING_T && s.max_penetration < 0.1 * RING_T, "{kind}: {deepest} {}", s.max_penetration);
        assert!((carried - punch).abs() < 0.03 * punch, "{kind}: contact {carried} vs punch {punch}");
        assert!(punch > 10.0 * bending, "{kind}: {punch} vs {bending}");
        if kind == "inner" {
            // The penalty law at the last increment: the force is ε W ∫ d dx over the overlap d
            // (the mean of the two rows of nodes), ε = γ min(m/A) / Δt² with the stable
            // increment of that evaluation, which lies among those of the last samples.
            let inner: Vec<ElementFace> = m.model.set("inner").unwrap().faces.clone();
            let n = r.history.increment.len();
            let recent = &r.history.increment[n - 11..n - 1];
            let (lo, hi) = recent.iter().fold((f64::INFINITY, 0.0f64), |(a, b), d| (a.min(*d), b.max(*d)));
            let ma = mass_per_area(&m, &inner);
            let integral = 0.5 * (positive_integral(&rows[0]) + positive_integral(&rows[1]));
            let force = |dt: f64| 0.3 * ma / (dt * dt) * RING_W * integral;
            eprintln!(
                "ring inner: min m/A {ma:.4e}; Δt {lo:.4e}..{hi:.4e}: ε W ∫ overlap {:.4e}..{:.4e} vs the contact force {last:.4e}; mean overlap {:.4e} vs mean pressure / ε {:.4e}..{:.4e}",
                force(hi),
                force(lo),
                integral / 0.09,
                last / (RING_W * 0.09) / (0.3 * ma / (lo * lo)),
                last / (RING_W * 0.09) / (0.3 * ma / (hi * hi)),
            );
            assert!(0.95 * force(hi) <= last && last <= 1.05 * force(lo), "{last} vs {}..{}", force(hi), force(lo));
        }
    }
}

// ---------------------------------------------------------------------------------------------
// No false contact.

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

/// A cantilever flung into large-amplitude bending: its compressed side never touches itself.
#[test]
fn a_bending_cantilever_feels_no_self_contact() {
    let mut m = Mesh::new();
    let beam = m.block(0, [0.0; 3], [1.0, 0.1, 0.1], [20, 2, 2], |x| x);
    m.material("Mat", 1.0e4, 0.3, 1.0, true);
    m.section("beam", &beam, "Mat");
    m.node_set("root", m.nodes_where(&beam, |x| x[0].abs() < 1e-9));
    m.general("self", None, Friction::coulomb(0.2));
    let mut step = ExplicitDynamic::new(0.3);
    step.initial_velocities = vec![InitialVelocity { set: "beam".into(), velocity: [0.0, 0.0, 5.0] }];
    m.step(step, &[("root", BoundaryCondition::Fixed { target: Target::Set("root".into()) })]);
    let tip = m.nodes_where(&beam, |x| (x[0] - 1.0).abs() < 1e-9);
    let r = m.solve();
    let f = r.frames.iter().map(|f| tip.iter().map(|n| f.displacement[index(&r, *n)][2]).fold(0.0, f64::max)).fold(0.0, f64::max);
    eprintln!("cantilever: largest tip deflection {f:.3}");
    assert!(f > 0.4, "{f}");
    same_without_contact(m, "cantilever");
}

/// A thick tube with a bore of a third of its outer diameter, squeezed into an oval by its own
/// motion: the faces across the bore face each other and never touch.
#[test]
fn an_ovalising_tube_feels_no_self_contact() {
    use std::f64::consts::PI;
    let mut m = Mesh::new();
    let tube = m.block(0, [0.0, 0.05, 0.0], [2.0 * PI, 0.1, 0.1], [24, 2, 1], |p| [p[1] * p[0].cos(), p[2], p[1] * p[0].sin()]);
    m.material("Mat", 1.0e4, 0.3, 1.0, true);
    m.section("tube", &tube, "Mat");
    let top: Vec<u32> = tube.iter().copied().filter(|e| m.model.element(*e).unwrap().nodes.iter().all(|n| m.xyz(*n)[2] > 0.0)).collect();
    let bottom: Vec<u32> = tube.iter().copied().filter(|e| m.model.element(*e).unwrap().nodes.iter().all(|n| m.xyz(*n)[2] < 0.0)).collect();
    m.model.add_set("top", Set { elements: top, ..Set::default() }).unwrap();
    m.model.add_set("bottom", Set { elements: bottom, ..Set::default() }).unwrap();
    m.general("self", None, Friction::NONE);
    let mut step = ExplicitDynamic::new(0.05);
    step.initial_velocities =
        vec![InitialVelocity { set: "top".into(), velocity: [0.0, 0.0, -0.5] }, InitialVelocity { set: "bottom".into(), velocity: [0.0, 0.0, 0.5] }];
    m.step(step, &[]);
    same_without_contact(m, "tube");
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

/// A stack of `layers` unit hex8 whose nodes are all prescribed to a layer squeezed unevenly: the
/// node at level k over corner c ends at height k h_c, so the stack is `h_c` thick over corner c
/// (the corners of the unit square in the element's order). Every element keeps a positive volume
/// and nothing crosses, but where `h_c` is small the far exterior face is only `layers h_c` away
/// from the near one, however thick the element is on average. The stack moves as without
/// contact and takes no contact force.
fn tapered_stack(layers: usize, h: [f64; 4]) {
    let mut m = Mesh::new();
    let stack = m.block(0, [0.0; 3], [1.0, 1.0, layers as f64], [1, 1, layers], |x| x);
    m.material("Foam", 1.0e3, 0.0, 1.0, false);
    m.section("stack", &stack, "Foam");
    m.node_set("all_nodes", m.nodes_where(&stack, |_| true));
    m.general("self", None, Friction::NONE);
    let mut bcs: Vec<(String, BoundaryCondition)> = vec![("die".into(), held("all_nodes", &[0, 1], 0.0))];
    for n in m.nodes_where(&stack, |_| true) {
        let x = m.xyz(n);
        let corner = match [x[0] > 0.5, x[1] > 0.5] {
            [false, false] => 0,
            [true, false] => 1,
            [true, true] => 2,
            [false, true] => 3,
        };
        let level = x[2];
        m.node_set(&format!("n{n}"), vec![n]);
        bcs.push((format!("z{n}"), held(&format!("n{n}"), &[2], level * h[corner] - level)));
    }
    let bcs: Vec<(&str, BoundaryCondition)> = bcs.iter().map(|(n, b)| (n.as_str(), b.clone())).collect();
    m.step(ExplicitDynamic::new(0.2), &bcs);
    same_without_contact(m, &format!("{layers} layers tapered {h:?}"));
}

/// Two hex8 whose middle surface is 0.001 high over three corners and 1.0 over the fourth, and
/// whose top surface is twice as high (automated review of #1120): the lower element is 0.251
/// thick on average, but over the thin corners its far side is 0.001 away and the next element's
/// exterior face 0.002 away, so a limit of half the average thickness reached it.
#[test]
fn a_tapered_two_element_stack_feels_no_self_contact() {
    tapered_stack(2, [0.001, 1.0, 0.001, 0.001]);
}

/// The same with three elements through the thickness: the far exterior face is three thin
/// elements away over the thin corners.
#[test]
fn a_tapered_three_element_stack_feels_no_self_contact() {
    tapered_stack(3, [0.001, 1.0, 0.001, 0.001]);
}

/// What self contact costs: a slab of n × n × 2 hex8 tumbling freely, run for 200 increments with
/// general contact on all its exterior faces and without, for growing n. The contact's time per
/// increment over the number of exterior faces stays about constant (the search and the
/// integration are linear in the faces: a bucket grid, a fixed number of candidates per face).
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn self_contact_cost_scales_with_the_exterior_faces() {
    for n in [8usize, 16, 32, 64] {
        let run = |contact: bool| {
            let mut m = Mesh::new();
            let slab = m.block(0, [0.0; 3], [1.0, 1.0, 2.0 / n as f64], [n, n, 2], |x| x);
            m.material("Mat", 1.0e4, 0.3, 1.0, false);
            m.section("slab", &slab, "Mat");
            if contact {
                m.general("self", None, Friction::coulomb(0.2));
            }
            let c = (1.0e4f64).sqrt();
            let dt = 0.9 * (1.0 / n as f64) / c;
            let mut step = ExplicitDynamic::new(200.0 * dt * 0.6);
            step.max_increments = 200;
            step.initial_velocities = vec![InitialVelocity { set: "slab".into(), velocity: [1.0, -0.5, 0.2] }];
            step.output_interval = Some(1.0);
            m.step(step, &[]);
            let faces = m.model.exterior_faces().len();
            let start = std::time::Instant::now();
            let r = solve_explicit(&m.model, Some("Run"), &Limits::default()).unwrap();
            (start.elapsed().as_secs_f64(), r.increments, faces)
        };
        let (with, increments, faces) = run(true);
        let (without, _, _) = run(false);
        let per = (with - without) / increments as f64 / faces as f64;
        eprintln!(
            "slab {n} × {n} × 2: {faces} exterior faces, {increments} increments; {with:.3} s with self contact, {without:.3} s without; contact {:.3} µs per face per increment",
            per * 1e6
        );
    }
}
