//! Deformable contact and friction in explicit dynamics against closed-form solutions:
//! - two identical elastic bars colliding exchange their velocities after a contact of 2L/c, with
//!   the interface force ρ c A v₀/2 (one-dimensional wave theory: K. F. Graff, *Wave Motion in
//!   Elastic Solids*, Dover 1991, §2.1–2.3; W. Johnson, *Impact Strength of Materials*, Arnold
//!   1972, §1.6), and keep their energy (frictionless elastic impact);
//! - a block sliding on a plane with Coulomb friction decelerates at μg and stops after
//!   v₀²/(2μg), on a rigid wall and on a deformable plate (any dynamics text, e.g. J. L. Meriam,
//!   L. G. Kraige, *Engineering Mechanics: Dynamics*, §3/5), and with the exponential decay law
//!   after the distance the decaying coefficient gives;
//! - two elastic spheres colliding: Hertz's contact duration and peak force (H. Hertz, *Über die
//!   Berührung fester elastischer Körper*, J. reine angew. Math. 92:156–171, 1882; K. L. Johnson,
//!   *Contact Mechanics*, Cambridge 1985, §11.4);
//! - a contact patch test: two blocks with non-matching meshes pressed together carry a uniform
//!   stress across the interface (R. L. Taylor, P. Papadopoulos, *On a patch test for contact
//!   problems in two dimensions*, in Nonlinear Computational Mechanics, Springer 1991, 690–702;
//!   T. Y. Crisfield-style uniform-pressure form for three dimensions);
//! - and what a contact must never do: hold back a body in free flight, depend on the thread
//!   count, or crash on hostile input.

use fieldcraft_explicit::{ExplicitResult, Limits, solve_explicit};
use fieldcraft_model::*;

/// A model built block by block, each block of hex8 with its own nodes.
struct Mesh {
    model: Model,
    next_node: u32,
    next_element: u32,
}

impl Mesh {
    fn new() -> Mesh {
        Mesh { model: Model::default(), next_node: 1, next_element: 1 }
    }

    /// A box from `origin` of `size` with `n` hex8 along each axis, positions mapped by `map`;
    /// its elements.
    fn block_mapped(&mut self, origin: [f64; 3], size: [f64; 3], n: [usize; 3], map: impl Fn([f64; 3]) -> [f64; 3]) -> Vec<u32> {
        let base = self.next_node;
        let at = |i: usize, j: usize, k: usize| base + (i + (n[0] + 1) * (j + (n[1] + 1) * k)) as u32;
        for k in 0..=n[2] {
            for j in 0..=n[1] {
                for i in 0..=n[0] {
                    let x = [0, 1, 2].map(|d| origin[d] + size[d] * [i, j, k][d] as f64 / n[d] as f64);
                    self.model.set_node(at(i, j, k), map(x)).unwrap();
                }
            }
        }
        self.next_node = at(n[0], n[1], n[2]) + 1;
        let mut out = Vec::new();
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    let nodes = vec![
                        at(i, j, k),
                        at(i + 1, j, k),
                        at(i + 1, j + 1, k),
                        at(i, j + 1, k),
                        at(i, j, k + 1),
                        at(i + 1, j, k + 1),
                        at(i + 1, j + 1, k + 1),
                        at(i, j + 1, k + 1),
                    ];
                    self.model.set_element(self.next_element, ElementKind::Hex8, nodes).unwrap();
                    out.push(self.next_element);
                    self.next_element += 1;
                }
            }
        }
        out
    }

    fn block(&mut self, origin: [f64; 3], size: [f64; 3], n: [usize; 3]) -> Vec<u32> {
        self.block_mapped(origin, size, n, |x| x)
    }

    fn xyz(&self, node: u32) -> [f64; 3] {
        self.model.node(node).unwrap().xyz
    }

    /// The faces of `elements` whose centre satisfies `f`.
    fn faces_where(&self, elements: &[u32], f: impl Fn([f64; 3]) -> bool) -> Vec<ElementFace> {
        let mut out = Vec::new();
        for &e in elements {
            let el = self.model.element(e).unwrap();
            for (i, local) in el.kind.faces().iter().enumerate() {
                let c = [0, 1, 2].map(|d| local.iter().map(|n| self.xyz(el.nodes[*n])[d]).sum::<f64>() / local.len() as f64);
                if f(c) {
                    out.push(ElementFace { element: e, face: i as u8 + 1 });
                }
            }
        }
        out
    }

    fn nodes_where(&self, elements: &[u32], f: impl Fn([f64; 3]) -> bool) -> Vec<u32> {
        let mut out: Vec<u32> = elements.iter().flat_map(|e| self.model.element(*e).unwrap().nodes.clone()).filter(|n| f(self.xyz(*n))).collect();
        out.sort_unstable();
        out.dedup();
        out
    }

    fn material(&mut self, name: &str, young: f64, poisson: f64, density: f64) {
        let m = Material { elastic: Some(Elastic::Isotropic { young, poisson }), density: Some(density), ..Material::default() };
        self.model.add_material(name, m).unwrap();
    }

    /// A solid section of one-point hex8 of `material` on `elements`, under the set `name`.
    fn section(&mut self, name: &str, elements: &[u32], material: &str) {
        self.model.add_set(name, Set { elements: elements.to_vec(), ..Set::default() }).unwrap();
        let formulation = Formulation { integration: Integration::Reduced, ..Formulation::default() };
        self.model
            .add_section(name, Section { orientation: None, set: name.into(), material: material.into(), kind: SectionKind::Solid {}, formulation })
            .unwrap();
    }

    fn face_set(&mut self, name: &str, faces: Vec<ElementFace>) {
        assert!(!faces.is_empty(), "{name}");
        self.model.add_set(name, Set { faces, ..Set::default() }).unwrap();
    }

    fn node_set(&mut self, name: &str, nodes: Vec<u32>) {
        assert!(!nodes.is_empty(), "{name}");
        self.model.add_set(name, Set { nodes, ..Set::default() }).unwrap();
    }

    fn pair(&mut self, name: &str, main: &str, secondary: &str, friction: Friction) {
        let c = Contact { friction, ..Contact::new(ContactSurfaces::Pair { main: main.into(), secondary: secondary.into() }) };
        self.model.add_contact(name, c).unwrap();
    }

    fn step(&mut self, step: ExplicitDynamic, bcs: &[(&str, BoundaryCondition)], loads: &[(&str, Load)]) {
        let mut s = Step::new("Run", Procedure::ExplicitDynamic(step));
        for (name, bc) in bcs {
            self.model.add_boundary_condition(name, bc.clone()).unwrap();
            s.boundary_conditions.insert(name.to_string());
        }
        for (name, load) in loads {
            self.model.add_load(name, load.clone()).unwrap();
            s.loads.insert(name.to_string());
        }
        self.model.add_step(s, None).unwrap();
    }

    fn solve(&self) -> ExplicitResult {
        let r = solve_explicit(&self.model, Some("Run"), &Limits::default()).unwrap();
        assert!(r.completed, "{:?}", r.warnings);
        r
    }
}

fn node_index(r: &ExplicitResult, label: u32) -> usize {
    r.nodes.binary_search(&label).unwrap()
}

/// The mean velocity of some nodes in the last frame.
fn mean_velocity(r: &ExplicitResult, nodes: &[u32]) -> [f64; 3] {
    let f = r.frames.last().unwrap();
    let mut v = [0.0; 3];
    for n in nodes {
        let x = f.velocity[node_index(r, *n)];
        for k in 0..3 {
            v[k] += x[k] / nodes.len() as f64;
        }
    }
    v
}

const BAR_E: f64 = 1.0e4;
const BAR_RHO: f64 = 1.0;
const BAR_L: f64 = 10.0;
const BAR_W: f64 = 0.1;

/// Two bars of `cells` hex8, the first moving at v₀ along x into the second at rest (touching
/// end to end, separate nodes), with `contact` between their ends: (mesh, bar A, bar B).
fn bars(cells: usize, v0: f64, end: f64, contact: impl FnOnce(&mut Mesh, Vec<ElementFace>, Vec<ElementFace>)) -> (Mesh, Vec<u32>, Vec<u32>) {
    let mut m = Mesh::new();
    let a = m.block([0.0; 3], [BAR_L, BAR_W, BAR_W], [cells, 1, 1]);
    let b = m.block([BAR_L, 0.0, 0.0], [BAR_L, BAR_W, BAR_W], [cells, 1, 1]);
    m.material("Mat", BAR_E, 0.0, BAR_RHO);
    m.section("A", &a, "Mat");
    m.section("B", &b, "Mat");
    let end_a = m.faces_where(&a, |c| (c[0] - BAR_L).abs() < 1e-9);
    let start_b = m.faces_where(&b, |c| (c[0] - BAR_L).abs() < 1e-9);
    contact(&mut m, end_a, start_b);
    let mut step = ExplicitDynamic::new(end);
    step.initial_velocities = vec![InitialVelocity { set: "A".into(), velocity: [v0, 0.0, 0.0] }];
    // The elastic solution: no linear bulk viscosity (it dissipates at the wave fronts).
    step.bulk_viscosity = BulkViscosity { linear: 0.0, quadratic: 1.5 };
    step.history_interval = Some(end / 4000.0);
    m.step(step, &[], &[]);
    (m, a, b)
}

/// Velocity exchange, contact duration 2L/c, interface force ρ c A v₀ / 2 and the energy of a
/// frictionless elastic impact, for a pair (both ways) and for general contact.
#[test]
fn colliding_bars_exchange_their_velocities() {
    let c = (BAR_E / BAR_RHO).sqrt();
    let v0 = 0.1;
    let duration = 2.0 * BAR_L / c;
    for kind in ["pair", "general"] {
        let (m, a, b) = bars(200, v0, 1.5 * duration, |m, end_a, start_b| match kind {
            "pair" => {
                m.face_set("end_a", end_a);
                m.face_set("start_b", start_b);
                m.pair("impact", "start_b", "end_a", Friction::NONE);
            }
            _ => m.model.add_contact("everything", Contact::new(ContactSurfaces::General { surface: None })).unwrap(),
        });
        let r = m.solve();
        let s = &r.contacts[0];
        let (first, last) = (s.first_contact.unwrap(), s.last_contact.unwrap());
        let plateau = BAR_RHO * c * BAR_W * BAR_W * v0 / 2.0;
        let va = mean_velocity(&r, &m.nodes_where(&a, |_| true))[0];
        let vb = mean_velocity(&r, &m.nodes_where(&b, |_| true))[0];
        let ke = r.energy.kinetic / r.initial_kinetic;
        // The force across the interface while the waves cross the bars.
        let mid: Vec<f64> = r
            .history
            .time
            .iter()
            .zip(&r.history.contact_force)
            .filter(|(t, _)| **t > 0.2 * duration && **t < 0.8 * duration)
            .map(|(_, f)| f[0])
            .collect();
        let mean_force = mid.iter().sum::<f64>() / mid.len() as f64;
        eprintln!(
            "{kind}: contact {first:.5}..{last:.5} ({:.5}) vs 2L/c = {duration:.5}; force {mean_force:.5e} (peak {:.5e}) vs ρcAv₀/2 = {plateau:.5e}; v_A {va:.5} v_B {vb:.5} vs 0, {v0}; KE/KE₀ {ke:.5}; contact energy {:.3e} of KE₀ {:.3e}; balance {:.2e}; max penetration {:.3e}; {} increments",
            last - first,
            s.peak_force,
            r.energy.contact,
            r.initial_kinetic,
            r.energy.relative_error(r.initial_kinetic),
            s.max_penetration,
            r.increments
        );
        assert!(first < 0.01 * duration, "{kind}: contact from the start, not {first}");
        assert!(((last - first) - duration).abs() < 0.015 * duration, "{kind}: contact for {} vs {duration}", last - first);
        assert!((mean_force - plateau).abs() < 0.02 * plateau, "{kind}: {mean_force} vs {plateau}");
        assert!(va.abs() < 0.02 * v0 && (vb - v0).abs() < 0.02 * v0, "{kind}: {va} {vb}");
        assert!(ke > 0.97 && ke <= 1.0 + 1e-6, "{kind}: KE ratio {ke}");
        // After the separation the contact springs are empty again: what they took they gave back.
        assert!(r.energy.contact.abs() < 0.005 * r.initial_kinetic, "{kind}: {:?}", r.energy);
        assert!(r.energy.relative_error(r.initial_kinetic) < 0.005, "{kind}: {:?}", r.energy);
        assert!(s.max_penetration < 0.01 * BAR_L / 200.0 * 10.0, "{kind}: penetration {}", s.max_penetration);
        assert_eq!(r.energy.friction, 0.0);
    }
}

/// A body in free flight with general contact on all its faces feels no contact force (no node
/// grabs a face of its own body), and its motion is that of the same run without contact.
#[test]
fn general_contact_leaves_a_free_body_alone() {
    let run = |contact: bool| {
        let mut m = Mesh::new();
        let e = m.block([0.0; 3], [1.0, 0.6, 0.4], [5, 3, 2]);
        m.material("Mat", 1e4, 0.3, 1.0);
        m.section("all", &e, "Mat");
        if contact {
            m.model.add_contact("all", Contact::new(ContactSurfaces::General { surface: None })).unwrap();
        }
        let mut step = ExplicitDynamic::new(0.05);
        step.initial_velocities = vec![InitialVelocity { set: "all".into(), velocity: [1.0, -0.5, 0.2] }];
        m.step(step, &[], &[]);
        m.solve()
    };
    let (with, without) = (run(true), run(false));
    assert!(with.history.contact_force.iter().all(|f| f[0] == 0.0));
    assert_eq!(with.contacts[0].first_contact, None);
    assert_eq!(with.frames.last().unwrap().displacement, without.frames.last().unwrap().displacement);
    assert_eq!(with.energy.contact, 0.0);
}

/// What kind of plane a block slides on.
#[derive(Clone, Copy, PartialEq, Debug)]
enum Plane {
    RigidWall,
    /// A plate fixed at its base, in contact with the block.
    Plate,
}

/// A block (1 × 1 × 0.5, 4 × 4 × 2 hex8) launched along x at v₀ on a horizontal plane under
/// gravity g with `friction`: the distance its centre of mass travels and the time it stops.
fn slide(plane: Plane, friction: Friction, v0: f64, g: f64, end: f64) -> (f64, f64, ExplicitResult) {
    let mut m = Mesh::new();
    let block = m.block([0.0, 0.0, 0.0], [1.0, 1.0, 0.5], [4, 4, 2]);
    m.material("Mat", 1.0e5, 0.3, 1.0);
    m.section("block", &block, "Mat");
    let bottom = m.nodes_where(&block, |x| x[2].abs() < 1e-9);
    let mut step = ExplicitDynamic::new(end);
    step.initial_velocities = vec![InitialVelocity { set: "block".into(), velocity: [v0, 0.0, 0.0] }];
    step.history_interval = Some(end / 2000.0);
    let mut bcs = Vec::new();
    match plane {
        Plane::RigidWall => {
            m.node_set("bottom", bottom);
            step.rigid_walls = vec![RigidWall { set: "bottom".into(), point: [0.0; 3], normal: [0.0, 0.0, 1.0], friction }];
        }
        Plane::Plate => {
            let plate = m.block([-0.5, -0.5, -0.5], [4.5, 2.0, 0.5], [9, 4, 1]);
            m.section("plate", &plate, "Mat");
            let top = m.faces_where(&plate, |c| c[2].abs() < 1e-9);
            let under = m.faces_where(&block, |c| c[2].abs() < 1e-9);
            m.face_set("top", top);
            m.face_set("under", under);
            m.pair("slide", "top", "under", friction);
            m.node_set("base", m.nodes_where(&plate, |x| (x[2] + 0.5).abs() < 1e-9));
            bcs.push(("hold", BoundaryCondition::Fixed { target: Target::Set("base".into()) }));
        }
    }
    let gravity = Load::Gravity { acceleration: [0.0, 0.0, -g], set: Some("block".into()) };
    m.step(step, &bcs, &[("weight", gravity)]);
    let r = m.solve();
    let nodes = m.nodes_where(&block, |_| true);
    let f = r.frames.last().unwrap();
    let distance = nodes.iter().map(|n| f.displacement[node_index(&r, *n)][0]).sum::<f64>() / nodes.len() as f64;
    // The block's momentum per mass (the plate's is nothing along x when it is held).
    let mass = match plane {
        Plane::RigidWall => r.total_mass,
        Plane::Plate => 0.5,
    };
    let stop = r.history.time.iter().zip(&r.history.momentum).find(|(_, p)| p[0] / mass < 0.005 * v0).map_or(f64::NAN, |(t, _)| *t);
    (distance, stop, r)
}

/// Coulomb friction μ: deceleration μg, stopping after v₀/(μg) at v₀²/(2μg); on a rigid wall
/// (kinematic friction impulses) and on a deformable plate (penalty contact, return map).
#[test]
fn a_block_slides_to_a_stop_under_coulomb_friction() {
    let (v0, g, mu) = (3.0, 10.0, 0.3);
    let (t_stop, d_stop) = (v0 / (mu * g), v0 * v0 / (2.0 * mu * g));
    for plane in [Plane::RigidWall, Plane::Plate] {
        let (d, t, r) = slide(plane, Friction::coulomb(mu), v0, g, 1.3 * t_stop);
        eprintln!(
            "{plane:?}: stops at {t:.4} vs v₀/(μg) = {t_stop:.4}, after {d:.4} vs v₀²/(2μg) = {d_stop:.4}; friction dissipated {:.4} vs ½mv₀² = {:.4}; balance {:.2e}; {} increments",
            r.energy.friction,
            0.5 * 0.5 * v0 * v0,
            r.energy.relative_error(r.initial_kinetic),
            r.increments
        );
        assert!((t - t_stop).abs() < 0.02 * t_stop, "{plane:?}: {t} vs {t_stop}");
        assert!((d - d_stop).abs() < 0.02 * d_stop, "{plane:?}: {d} vs {d_stop}");
        // The kinetic energy went into friction (and a little into the bouncing on the plane).
        let ke0 = 0.5 * 0.5 * v0 * v0;
        assert!((r.energy.friction - ke0).abs() < 0.03 * ke0, "{plane:?}: {:?}", r.energy);
        assert!(r.energy.relative_error(r.initial_kinetic) < 0.01, "{plane:?}: {:?}", r.energy);
    }
}

/// The exponential decay law μ(v) = μ_k + (μ_s − μ_k) e^{−d v}: the block decelerates at
/// g μ(v), so it stops after t = ∫ dv / (g μ(v)) = [v₀ + ln((μ_k + Δμ e^{−d v₀}) / μ_s) / d] / (g μ_k)
/// and travels ∫ v dv / (g μ(v)) (by quadrature).
#[test]
fn a_block_slides_with_a_decaying_friction_coefficient() {
    let (v0, g) = (3.0, 10.0);
    let f = Friction { static_coefficient: 0.5, kinetic: Some(0.2), decay: 1.5 };
    let (a, b, d): (f64, f64, f64) = (0.2, 0.3, 1.5);
    let t_stop = (v0 + ((a + b * (-d * v0).exp()) / (a + b)).ln() / d) / (g * a);
    let n = 20_000;
    let d_stop: f64 = (0..n)
        .map(|i| {
            let v = (i as f64 + 0.5) * v0 / n as f64;
            v / (g * f.coefficient(v)) * v0 / n as f64
        })
        .sum();
    for plane in [Plane::RigidWall, Plane::Plate] {
        let (dist, t, _) = slide(plane, f, v0, g, 1.3 * t_stop);
        eprintln!("{plane:?}, decaying μ: stops at {t:.4} vs {t_stop:.4}, after {dist:.4} vs {d_stop:.4}");
        assert!((t - t_stop).abs() < 0.03 * t_stop, "{plane:?}: {t} vs {t_stop}");
        assert!((dist - d_stop).abs() < 0.03 * d_stop, "{plane:?}: {dist} vs {d_stop}");
    }
}

/// Two blocks with non-matching meshes (3 × 3 below, 2 × 2 above) pressed together slowly
/// between symmetry planes: the stress is uniform σ_zz on both sides of the interface and the
/// force across it is σ_zz A (a contact patch test in uniaxial stress).
#[test]
fn a_contact_patch_test_carries_a_uniform_stress_across_non_matching_meshes() {
    let mut m = Mesh::new();
    let lower = m.block([0.0; 3], [1.0, 1.0, 1.0], [3, 3, 2]);
    let upper = m.block([0.0, 0.0, 1.0], [1.0, 1.0, 1.0], [2, 2, 2]);
    let (young, rho) = (1000.0, 1.0);
    m.material("Mat", young, 0.3, rho);
    m.section("lower", &lower, "Mat");
    m.section("upper", &upper, "Mat");
    let both: Vec<u32> = lower.iter().chain(&upper).copied().collect();
    m.face_set("lower_top", m.faces_where(&lower, |c| (c[2] - 1.0).abs() < 1e-9));
    m.face_set("upper_bottom", m.faces_where(&upper, |c| (c[2] - 1.0).abs() < 1e-9));
    m.pair("interface", "lower_top", "upper_bottom", Friction::NONE);
    m.node_set("x0", m.nodes_where(&both, |x| x[0].abs() < 1e-9));
    m.node_set("y0", m.nodes_where(&both, |x| x[1].abs() < 1e-9));
    m.node_set("z0", m.nodes_where(&both, |x| x[2].abs() < 1e-9));
    m.node_set("top", m.nodes_where(&both, |x| (x[2] - 2.0).abs() < 1e-9));
    let squeeze = 0.01;
    let held = |set: &str, k: usize| {
        let mut values = [None; 6];
        values[k] = Some(0.0);
        BoundaryCondition::Displacement { target: Target::Set(set.into()), values }
    };
    let mut push_values = [None; 6];
    push_values[2] = Some(-squeeze);
    let push = BoundaryCondition::Displacement { target: Target::Set("top".into()), values: push_values };
    // A hundred wave transits: quasi-static.
    let c = (young * 0.7 / (1.3 * 0.4) / rho).sqrt();
    let end = 100.0 * 2.0 / c;
    let mut step = ExplicitDynamic::new(end);
    step.output_interval = Some(end / 4.0);
    m.step(step, &[("x0", held("x0", 0)), ("y0", held("y0", 1)), ("z0", held("z0", 2)), ("push", push)], &[]);
    let r = m.solve();
    let f = r.frames.last().unwrap();
    let szz: Vec<f64> = f.stress.iter().map(|s| s[2]).collect();
    let mean = szz.iter().sum::<f64>() / szz.len() as f64;
    let spread = szz.iter().map(|s| (s - mean).abs()).fold(0.0, f64::max);
    let lateral = f.stress.iter().flat_map(|s| [s[0], s[1], s[3], s[4], s[5]]).map(f64::abs).fold(0.0, f64::max);
    let force = *r.history.contact_force.last().unwrap().first().unwrap();
    eprintln!(
        "patch: σ_zz mean {mean:.5} (E ε = {:.5}), spread {spread:.3e}, largest other component {lateral:.3e}; force across {force:.5} vs |σ_zz| A = {:.5}",
        -young * squeeze / 2.0,
        mean.abs()
    );
    assert!(spread < 0.02 * mean.abs(), "σ_zz {szz:?}");
    assert!(lateral < 0.02 * mean.abs(), "{lateral}");
    assert!((force - mean.abs()).abs() < 0.02 * mean.abs(), "{force} vs {mean}");
    // The interface's penalty takes some of the squeeze (uniformly): each block is strained by
    // what is left.
    let closure = r.contacts[0].max_penetration;
    let expected = young * (squeeze - closure) / 2.0;
    eprintln!("patch: interface closure {closure:.4e}, σ_zz vs E (squeeze − closure) / 2 = {expected:.5}");
    assert!((mean.abs() - expected).abs() < 0.01 * expected, "{mean} vs {expected}");
}

/// The same contact run gives the same numbers, bit for bit, on one thread and on four.
#[test]
fn contact_results_do_not_depend_on_the_thread_count() {
    let c = (BAR_E / BAR_RHO).sqrt();
    let (m, _, _) = bars(60, 0.1, 2.0 * BAR_L / c, |m, _, _| {
        m.model.add_contact("all", Contact { friction: Friction::coulomb(0.2), ..Contact::new(ContactSurfaces::General { surface: None }) }).unwrap()
    });
    let one = fieldcraft_parallel::with_threads(1, || solve_explicit(&m.model, None, &Limits::default()).unwrap());
    let four = fieldcraft_parallel::with_threads(4, || solve_explicit(&m.model, None, &Limits::default()).unwrap());
    assert_eq!(one.frames, four.frames);
    assert_eq!(one.history, four.history);
    assert_eq!(one.energy, four.energy);
    assert_eq!(one.contacts, four.contacts);
}

/// Nodes that start inside the other surface are moved onto it and reported, or refused with
/// the node, the face and the depth; a cap on contact faces is honoured; contact on the faces
/// of tetrahedra works too.
#[test]
fn initial_penetration_caps_and_tetrahedra() {
    let build = |mode: InitialPenetration, overlap: f64| {
        let mut m = Mesh::new();
        let a = m.block([0.0; 3], [1.0, 1.0, 1.0], [2, 2, 2]);
        let b = m.block([0.0, 0.0, 1.0 - overlap], [1.0, 1.0, 1.0], [3, 3, 2]);
        m.material("Mat", 1000.0, 0.3, 1.0);
        m.section("a", &a, "Mat");
        m.section("b", &b, "Mat");
        m.face_set("a_top", m.faces_where(&a, |c| (c[2] - 1.0).abs() < 1e-9));
        m.face_set("b_bottom", m.faces_where(&b, |c| (c[2] - 1.0 + overlap).abs() < 1e-9));
        let c = Contact { initial_penetration: mode, ..Contact::new(ContactSurfaces::Pair { main: "a_top".into(), secondary: "b_bottom".into() }) };
        m.model.add_contact("stack", c).unwrap();
        m.step(ExplicitDynamic::new(0.01), &[], &[]);
        m
    };
    let r = build(InitialPenetration::Adjust, 0.01).solve();
    let s = &r.contacts[0];
    // Both passes: b's 16 bottom nodes into a, and a's 9 top nodes into b (they are inside b
    // after b's nodes moved down? no: b's nodes are moved onto a's top first, then a's top nodes
    // are on b's bottom).
    assert!(!s.adjusted.is_empty(), "{s:?}");
    assert!(s.adjusted.iter().all(|a| (a.depth - 0.01).abs() < 1e-9), "{:?}", s.adjusted);
    assert!(r.warnings.iter().any(|w| w.contains("contact `stack`") && w.contains("moved onto it")), "{:?}", r.warnings);
    // After the adjustment nothing is pushed apart: no contact force worth the name.
    assert!(s.peak_force < 1e-6, "{s:?}");
    let e = solve_explicit(&build(InitialPenetration::Refuse, 0.01).model, None, &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("contact `stack`") && e.contains("starts 1.0000e-2 inside face 2 of element"), "{e}");
    // Touching exactly is not a penetration.
    assert!(build(InitialPenetration::Refuse, 0.0).solve().contacts[0].adjusted.is_empty());
    let capped = Limits { max_contact_faces: 10, ..Limits::default() };
    let e = solve_explicit(&build(InitialPenetration::Adjust, 0.0).model, None, &capped).unwrap_err().to_string();
    assert!(e.contains("limit") && e.contains("contact `stack`"), "{e}");
}

/// Interval ends on [0, 1] in `n` intervals growing from the first by `ratio`.
fn graded(n: usize, ratio: f64) -> Vec<f64> {
    let widths: Vec<f64> = (0..n).map(|k| ratio.powi(k as i32)).collect();
    let total: f64 = widths.iter().sum();
    let mut out = vec![0.0];
    let mut at = 0.0;
    for w in widths {
        at += w / total;
        out.push(at);
    }
    if let Some(last) = out.last_mut() {
        *last = 1.0;
    }
    out
}

/// A quarter (x ≥ 0, y ≥ 0) of a ball of radius `radius` about `center`, meshed as a cubed sphere
/// of hex8: an inner cube of half-size 0.4 R with `n` uniform elements per edge, and six blocks
/// from its faces to the sphere with `n` × `n` elements across (equiangular, clustered towards
/// each face's centre by `lateral`) and `layers` elements through (thinner towards the surface
/// by `radial`). Its elements and the faces on the sphere within 30° of the pole below (−z) or
/// above (+z) the centre (`towards`).
fn quarter_ball(
    m: &mut Mesh,
    center: [f64; 3],
    radius: f64,
    n: usize,
    layers: usize,
    lateral: f64,
    radial: f64,
    towards: f64,
) -> (Vec<u32>, Vec<ElementFace>) {
    use std::collections::BTreeMap;
    assert!(n.is_multiple_of(2));
    let inner = 0.4 * radius;
    // Clustered lateral parameter on [−1, 1]: symmetric, finest at 0.
    let half = graded(n / 2, lateral);
    let mut clustered: Vec<f64> = half.iter().rev().map(|s| -s).collect();
    clustered.extend(half.iter().skip(1));
    let uniform: Vec<f64> = (0..=n).map(|i| -1.0 + 2.0 * i as f64 / n as f64).collect();
    // Radial blend from the cube (0) to the sphere (1), finest at the sphere.
    let through: Vec<f64> = graded(layers, radial).iter().rev().map(|s| 1.0 - s).collect();
    let mut labels: BTreeMap<[i64; 3], u32> = BTreeMap::new();
    let quantum = 1e-9 * radius;
    let mut node = |m: &mut Mesh, x: [f64; 3]| -> u32 {
        let key = x.map(|c| (c / quantum).round() as i64);
        *labels.entry(key).or_insert_with(|| {
            let l = m.next_node;
            m.model.set_node(l, [center[0] + x[0], center[1] + x[1], center[2] + x[2]]).unwrap();
            m.next_node += 1;
            l
        })
    };
    let kept = |x: &[[f64; 3]]| {
        let c = [0, 1].map(|d| x.iter().map(|p| p[d]).sum::<f64>() / x.len() as f64);
        c[0] > 0.0 && c[1] > 0.0
    };
    let mut hexes: Vec<[[f64; 3]; 8]> = Vec::new();
    // The inner cube.
    let cube = |i: usize, j: usize, k: usize| [uniform[i], uniform[j], uniform[k]].map(|s| inner * s);
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                hexes.push([
                    cube(i, j, k),
                    cube(i + 1, j, k),
                    cube(i + 1, j + 1, k),
                    cube(i, j + 1, k),
                    cube(i, j, k + 1),
                    cube(i + 1, j, k + 1),
                    cube(i + 1, j + 1, k + 1),
                    cube(i, j + 1, k + 1),
                ]);
            }
        }
    }
    // The six blocks: face axis a, sign σ; (u, v) along the other two axes.
    for axis in 0..3 {
        for sign in [-1.0, 1.0] {
            let place = |u: f64, v: f64, w: f64| -> [f64; 3] {
                let mut p = [0.0; 3];
                p[axis] = w;
                p[(axis + 1) % 3] = u;
                p[(axis + 2) % 3] = v;
                p
            };
            let point = |i: usize, j: usize, k: usize| -> [f64; 3] {
                let t = through[k];
                let p_in = place(uniform[i], uniform[j], sign).map(|c| inner * c);
                let angle = |s: f64| (s * std::f64::consts::FRAC_PI_4).tan();
                let d = place(angle(clustered[i]), angle(clustered[j]), sign);
                let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                let p_out = d.map(|c| radius * c / len);
                [0, 1, 2].map(|c| (1.0 - t) * p_in[c] + t * p_out[c])
            };
            for k in 0..layers {
                for j in 0..n {
                    for i in 0..n {
                        hexes.push([
                            point(i, j, k),
                            point(i + 1, j, k),
                            point(i + 1, j + 1, k),
                            point(i, j + 1, k),
                            point(i, j, k + 1),
                            point(i + 1, j, k + 1),
                            point(i + 1, j + 1, k + 1),
                            point(i, j + 1, k + 1),
                        ]);
                    }
                }
            }
        }
    }
    let mut elements = Vec::new();
    for mut h in hexes.into_iter().filter(|h| kept(h)) {
        // Positive volume: the bottom face counter-clockwise seen from the top.
        let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let (a, b, c) = (sub(h[1], h[0]), sub(h[3], h[0]), sub(h[4], h[0]));
        let det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
        if det < 0.0 {
            h = [h[4], h[5], h[6], h[7], h[0], h[1], h[2], h[3]];
        }
        let nodes: Vec<u32> = h.iter().map(|x| node(m, *x)).collect();
        m.model.set_element(m.next_element, ElementKind::Hex8, nodes).unwrap();
        elements.push(m.next_element);
        m.next_element += 1;
    }
    let pole = [0.0, 0.0, towards];
    let exterior = m.model.exterior_faces();
    let faces: Vec<ElementFace> = m
        .faces_where(&elements, |c| {
            let r = [c[0] - center[0], c[1] - center[1], c[2] - center[2]];
            let len = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
            len > 0.9 * radius && (r[0] * pole[0] + r[1] * pole[1] + r[2] * pole[2]) / len > (std::f64::consts::PI / 6.0).cos()
        })
        .into_iter()
        .filter(|f| exterior.binary_search(f).is_ok())
        .collect();
    (elements, faces)
}

/// Hertz's impact of two identical elastic spheres (radius R, mass m, E, ν) at relative speed v:
/// with R* = R/2, m* = m/2, E* = E / (2(1 − ν²)), the approach reaches
/// δ = (15 m* v² / (16 E* √R*))^{2/5}, the force F = (4/3) E* √R* δ^{3/2}, and the contact lasts
/// 2.9432 δ / v (K. L. Johnson, *Contact Mechanics*, Cambridge 1985, §11.4).
fn hertz(radius: f64, young: f64, poisson: f64, density: f64, v: f64) -> (f64, f64, f64) {
    let mass = 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3) * density;
    let (r_eff, m_eff, e_eff) = (radius / 2.0, mass / 2.0, young / (2.0 * (1.0 - poisson * poisson)));
    let delta = (15.0 * m_eff * v * v / (16.0 * e_eff * r_eff.sqrt())).powf(0.4);
    (delta, 4.0 / 3.0 * e_eff * r_eff.sqrt() * delta.powf(1.5), 2.9432 * delta / v)
}

/// Two spheres in a head-on collision (a quarter of each, between the symmetry planes x = 0 and
/// y = 0), with a contact pair between their facing caps: the peak force and the contact
/// duration against Hertz. `n` and `layers` set the mesh; `a_over_r` the contact radius the speed
/// is chosen for.
fn sphere_impact(n: usize, layers: usize, lateral: f64, radial: f64, a_over_r: f64) -> (f64, f64, f64, f64, ExplicitResult) {
    let (radius, young, poisson, density) = (1.0, 1.0, 0.3, 1.0);
    // The speed that gives a contact radius a = √(R* δ).
    let delta = a_over_r * a_over_r * radius * radius / (radius / 2.0);
    let mass = 4.0 / 3.0 * std::f64::consts::PI * density;
    let e_eff = young / (2.0 * (1.0 - poisson * poisson));
    let v = (delta.powf(2.5) * 16.0 * e_eff * (radius / 2.0f64).sqrt() / (15.0 * mass / 2.0)).sqrt();
    let (delta, force, duration) = hertz(radius, young, poisson, density, v);
    let mut m = Mesh::new();
    let (upper, upper_cap) = quarter_ball(&mut m, [0.0, 0.0, radius], radius, n, layers, lateral, radial, -1.0);
    let (lower, lower_cap) = quarter_ball(&mut m, [0.0, 0.0, -radius], radius, n, layers, lateral, radial, 1.0);
    m.material("Mat", young, poisson, density);
    m.section("upper", &upper, "Mat");
    m.section("lower", &lower, "Mat");
    m.face_set("upper_cap", upper_cap);
    m.face_set("lower_cap", lower_cap);
    m.pair("impact", "lower_cap", "upper_cap", Friction::NONE);
    let both: Vec<u32> = upper.iter().chain(&lower).copied().collect();
    m.node_set("x0", m.nodes_where(&both, |x| x[0].abs() < 1e-9));
    m.node_set("y0", m.nodes_where(&both, |x| x[1].abs() < 1e-9));
    let held = |set: &str, k: usize| {
        let mut values = [None; 6];
        values[k] = Some(0.0);
        BoundaryCondition::Displacement { target: Target::Set(set.into()), values }
    };
    let mut step = ExplicitDynamic::new(1.15 * duration);
    step.initial_velocities = vec![
        InitialVelocity { set: "upper".into(), velocity: [0.0, 0.0, -v / 2.0] },
        InitialVelocity { set: "lower".into(), velocity: [0.0, 0.0, v / 2.0] },
    ];
    step.bulk_viscosity = BulkViscosity { linear: 0.0, quadratic: 1.5 };
    step.history_interval = Some(duration / 2000.0);
    step.output_interval = Some(1.15 * duration);
    m.step(step, &[("x0", held("x0", 0)), ("y0", held("y0", 1))], &[]);
    let r = m.solve();
    let s = &r.contacts[0];
    // A quarter of each sphere carries a quarter of the force.
    let peak = 4.0 * s.peak_force;
    let lasted = s.last_contact.unwrap_or(0.0) - s.first_contact.unwrap_or(0.0);
    eprintln!(
        "Hertz ({} elements, a/R = {a_over_r}): v = {v:.4e}, δ = {delta:.4e}; peak force {peak:.5e} vs {force:.5e} ({:+.2}%), duration {lasted:.4} vs {duration:.4} ({:+.2}%); max penetration {:.2e}; KE after/before {:.4}; balance {:.2e}; {} increments",
        both.len(),
        100.0 * (peak / force - 1.0),
        100.0 * (lasted / duration - 1.0),
        s.max_penetration,
        r.energy.kinetic / r.initial_kinetic,
        r.energy.relative_error(r.initial_kinetic),
        r.increments
    );
    (peak, force, lasted, duration, r)
}

/// A coarse sphere impact in the main gate: Hertz's force and duration to a few percent.
#[test]
fn spheres_collide_as_hertz_says_coarse() {
    let (peak, force, lasted, duration, _) = sphere_impact(12, 6, 1.25, 1.25, 0.15);
    assert!((peak - force).abs() < 0.1 * force, "{peak} vs {force}");
    assert!((lasted - duration).abs() < 0.1 * duration, "{lasted} vs {duration}");
}

/// The fine sphere impact: contact radius a tenth of the radius, about six elements across it.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn spheres_collide_as_hertz_says() {
    let (peak, force, lasted, duration, r) = sphere_impact(24, 12, 1.23, 1.185, 0.1);
    assert!((peak - force).abs() < 0.03 * force, "{peak} vs {force}");
    assert!((lasted - duration).abs() < 0.03 * duration, "{lasted} vs {duration}");
    assert!(r.energy.kinetic > 0.95 * r.initial_kinetic);
}
