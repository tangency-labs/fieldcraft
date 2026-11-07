//! Trusses, point masses, rotary inertia and springs against closed-form solutions. Run with
//! `--nocapture` to see the errors.
#![allow(clippy::needless_range_loop)]
//!
//! - Statically determinate trusses (a plane Pratt truss and a space tripod): the member forces
//!   are those of the method of joints, whatever the members' stiffness.
//! - The three-bar truss (statically indeterminate): N = P cos²θ / (1 + 2 cos³θ) in the inclined
//!   bars, P / (1 + 2 cos³θ) in the vertical one.
//! - rod3: constant-stress patch with off-centre middle nodes; a bar hanging under its own
//!   weight (quadratic displacement, exact at the nodes); curved members are checked in
//!   `fieldcraft_fem::truss` (a ring's radial expansion gives the hoop strain δ/R exactly).
//! - Spring–mass ω = √(k/m), springs along three skew directions, the two-mass chain, a bar with
//!   a tip mass ω = √(EA/(LM)), the torsional pendulum ω = √(k_t/J).

use std::f64::consts::PI;

use fieldcraft_fem::truss::axial;
use fieldcraft_model::*;
use fieldcraft_solver::{FrequencyOptions, FrequencyResult, Limits, MassKind, StaticResult, solve_frequency, solve_static};

const ZERO_Z: [Option<f64>; 6] = [None, None, Some(0.0), None, None, None];

/// A model built point by point.
struct Build {
    m: Model,
    step: &'static str,
}

impl Build {
    fn new(procedure: Procedure) -> Build {
        let mut m = Model::default();
        m.add_step(Step::new("Step", procedure), None).unwrap();
        Build { m, step: "Step" }
    }

    fn statics() -> Build {
        Build::new(Procedure::Static {})
    }

    fn modes(n: u32) -> Build {
        Build::new(Procedure::Frequency { modes: n, min_hz: None, max_hz: None })
    }

    fn node(&mut self, label: u32, xyz: [f64; 3]) {
        self.m.set_node(label, xyz).unwrap();
    }

    fn nodes(&mut self, name: &str, nodes: &[u32]) -> Target {
        self.m.add_set(name, Set { nodes: nodes.to_vec(), ..Set::default() }).unwrap();
        Target::Set(name.into())
    }

    /// Members with one area and material.
    fn members(&mut self, name: &str, kind: ElementKind, members: &[(u32, Vec<u32>)], area: f64, young: f64, density: Option<f64>) {
        for (label, nodes) in members {
            self.m.set_element(*label, kind, nodes.clone()).unwrap();
        }
        self.m.add_set(name, Set { elements: members.iter().map(|m| m.0).collect(), ..Set::default() }).unwrap();
        self.m.add_material(name, Material { elastic: Some(Elastic::Isotropic { young, poisson: 0.3 }), density, ..Material::default() }).unwrap();
        self.m
            .add_section(
                name,
                Section {
                    orientation: None,
                    set: name.into(),
                    material: name.into(),
                    kind: SectionKind::Rod { area },
                    formulation: Formulation::default(),
                },
            )
            .unwrap();
    }

    fn bc(&mut self, name: &str, target: Target, values: [Option<f64>; 6]) {
        self.m.add_boundary_condition(name, BoundaryCondition::Displacement { target, values }).unwrap();
        self.m.set_boundary_condition_active(self.step, name, true).unwrap();
    }

    fn fix(&mut self, name: &str, target: Target) {
        self.m.add_boundary_condition(name, BoundaryCondition::Fixed { target }).unwrap();
        self.m.set_boundary_condition_active(self.step, name, true).unwrap();
    }

    fn load(&mut self, name: &str, load: Load) {
        self.m.add_load(name, load).unwrap();
        self.m.set_load_active(self.step, name, true).unwrap();
    }

    fn force(&mut self, name: &str, node: u32, vector: [f64; 3]) {
        let target = self.nodes(name, &[node]);
        self.load(name, Load::Force { target, vector, distribution: Distribution::Total });
    }

    fn spring(&mut self, name: &str, first: Target, second: Option<Target>, direction: Direction, stiffness: f64) {
        self.m.add_spring(name, Spring::Spring { first, second, direction, stiffness }).unwrap();
    }

    fn mass(&mut self, name: &str, target: Target, mass: Option<Mass>, inertia: Option<[f64; 6]>) {
        self.m.add_mass(name, PointMass { target, mass, inertia }).unwrap();
    }

    fn solve(&self) -> StaticResult {
        solve_static(&self.m, None, &Limits::default()).unwrap()
    }

    fn modes_of(&self, mass: MassKind) -> FrequencyResult {
        solve_frequency(&self.m, None, &FrequencyOptions { mass, shift: None }, &Limits::default()).unwrap()
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn unit(a: [f64; 3]) -> [f64; 3] {
    let l = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    a.map(|c| c / l)
}

/// The axial force of each member at its first integration point: N = A tᵀ σ t.
fn member_forces(m: &Model, r: &StaticResult, area_of: impl Fn(u32) -> f64) -> Vec<(u32, f64)> {
    r.elements
        .iter()
        .map(|e| {
            let el = m.element(e.label).unwrap();
            let t = unit(sub(m.node(el.nodes[1]).unwrap().xyz, m.node(el.nodes[0]).unwrap().xyz));
            (e.label, area_of(e.label) * axial(&e.points[0].stress, t, false))
        })
        .collect()
}

/// Gaussian elimination with partial pivoting (the method of joints' equations).
fn solve_dense(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for c in 0..n {
        let p = (c..n).max_by(|i, j| a[*i][c].abs().total_cmp(&a[*j][c].abs())).unwrap();
        a.swap(c, p);
        b.swap(c, p);
        assert!(a[c][c].abs() > 1e-12, "the truss is not statically determinate");
        for r in c + 1..n {
            let f = a[r][c] / a[c][c];
            for k in c..n {
                a[r][k] -= f * a[c][k];
            }
            b[r] -= f * b[c];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        x[r] = (b[r] - (r + 1..n).map(|k| a[r][k] * x[k]).sum::<f64>()) / a[r][r];
    }
    x
}

/// Member forces by the method of joints: at every joint and direction, Σ N t + F + R = 0, with
/// the reactions `supports` (node, direction) as unknowns too.
fn method_of_joints(
    xyz: &[(u32, [f64; 3])],
    members: &[(u32, [u32; 2])],
    supports: &[(u32, usize)],
    loads: &[(u32, [f64; 3])],
    dims: usize,
) -> Vec<(u32, f64)> {
    let row = |node: u32, d: usize| xyz.iter().position(|(l, _)| *l == node).unwrap() * dims + d;
    let pos = |node: u32| xyz.iter().find(|(l, _)| *l == node).unwrap().1;
    let n = xyz.len() * dims;
    assert_eq!(members.len() + supports.len(), n, "determinate");
    let mut a = vec![vec![0.0; n]; n];
    let mut b = vec![0.0; n];
    for (k, (_, [i, j])) in members.iter().enumerate() {
        let t = unit(sub(pos(*j), pos(*i)));
        for d in 0..dims {
            // Tension pulls joint i towards j and j towards i.
            a[row(*i, d)][k] += t[d];
            a[row(*j, d)][k] -= t[d];
        }
    }
    for (k, (node, d)) in supports.iter().enumerate() {
        a[row(*node, *d)][members.len() + k] = 1.0;
    }
    for (node, f) in loads {
        for d in 0..dims {
            b[row(*node, d)] -= f[d];
        }
    }
    let x = solve_dense(a, b);
    members.iter().zip(x).map(|((l, _), n)| (*l, n)).collect()
}

fn worst_relative(fe: &[(u32, f64)], exact: &[(u32, f64)]) -> f64 {
    let scale = exact.iter().fold(0.0f64, |m, (_, n)| m.max(n.abs()));
    exact.iter().map(|(l, n)| (fe.iter().find(|(m, _)| m == l).unwrap().1 - n).abs() / scale).fold(0.0, f64::max)
}

/// A plane Pratt truss of four panels, pinned at one end and on a roller at the other, loaded at
/// the bottom joints; chords and web members have different areas. The out-of-plane direction is
/// held. Member forces equal the method of joints'.
#[test]
fn pratt_truss_member_forces_by_the_method_of_joints() {
    let (a, h, p) = (3.0, 4.0, 10.0);
    let mut b = Build::statics();
    let mut xyz = Vec::new();
    for i in 0..5u32 {
        xyz.push((1 + i, [a * f64::from(i), 0.0, 0.0]));
    }
    for i in 1..4u32 {
        xyz.push((10 + i, [a * f64::from(i), h, 0.0]));
    }
    for (l, x) in &xyz {
        b.node(*l, *x);
    }
    let chords: Vec<(u32, [u32; 2])> = vec![(1, [1, 2]), (2, [2, 3]), (3, [3, 4]), (4, [4, 5]), (5, [11, 12]), (6, [12, 13])];
    let web: Vec<(u32, [u32; 2])> = vec![(7, [2, 11]), (8, [3, 12]), (9, [4, 13]), (10, [1, 11]), (11, [13, 5]), (12, [11, 3]), (13, [13, 3])];
    let as_members = |list: &[(u32, [u32; 2])]| list.iter().map(|(l, n)| (*l, n.to_vec())).collect::<Vec<_>>();
    b.members("chords", ElementKind::Rod2, &as_members(&chords), 2.0, 2.0e5, None);
    b.members("web", ElementKind::Rod2, &as_members(&web), 0.7, 7.0e4, None);
    let all: Vec<u32> = xyz.iter().map(|x| x.0).collect();
    let plane = b.nodes("all_nodes", &all);
    b.bc("plane", plane, ZERO_Z);
    let pin = b.nodes("pin", &[1]);
    b.bc("pin", pin, [Some(0.0), Some(0.0), None, None, None, None]);
    let roller = b.nodes("roller", &[5]);
    b.bc("roller", roller, [None, Some(0.0), None, None, None, None]);
    let loads = [(2, [0.0, -p, 0.0]), (3, [0.0, -1.5 * p, 0.0]), (4, [0.2 * p, -p, 0.0])];
    for (node, f) in loads {
        b.force(&format!("P{node}"), node, f);
    }
    let r = b.solve();
    let fe = member_forces(&b.m, &r, |l| if l <= 6 { 2.0 } else { 0.7 });
    let members: Vec<(u32, [u32; 2])> = chords.iter().chain(&web).copied().collect();
    let exact = method_of_joints(&xyz, &members, &[(1, 0), (1, 1), (5, 1)], &loads, 2);
    let err = worst_relative(&fe, &exact);
    eprintln!("Pratt truss: worst member force error {err:.2e} of the largest");
    assert!(err < 1e-10, "{fe:?} vs {exact:?}");
    // Reactions balance the loads.
    for d in 0..2 {
        assert!((r.applied_total[d] + r.reaction_total[d]).abs() < 1e-9 * p);
    }
}

/// A space tripod: three members from an apex to three fixed points, a skew load at the apex.
#[test]
fn space_tripod_member_forces() {
    let (rad, h) = (2.0, 3.0);
    let mut b = Build::statics();
    let apex = (4, [0.1, -0.2, h]);
    let feet: Vec<(u32, [f64; 3])> = (0..3)
        .map(|i| (1 + i as u32, [rad * (PI / 2.0 + 2.0 * PI * i as f64 / 3.0).cos(), rad * (PI / 2.0 + 2.0 * PI * i as f64 / 3.0).sin(), 0.0]))
        .collect();
    for (l, x) in feet.iter().chain([&apex]) {
        b.node(*l, *x);
    }
    let members: Vec<(u32, Vec<u32>)> = (1..=3).map(|i| (i, vec![i, 4])).collect();
    b.members("legs", ElementKind::Rod2, &members, 1.5, 1.0e5, None);
    let feet_set = b.nodes("feet", &[1, 2, 3]);
    b.fix("feet", feet_set);
    let load = [3.0, -2.0, -10.0];
    b.force("P", 4, load);
    let r = b.solve();
    let fe = member_forces(&b.m, &r, |_| 1.5);
    // Equilibrium of the apex: Σ N_i t_i + P = 0 with t_i from the apex to foot i (tension
    // pulls the apex towards the foot).
    let t: Vec<[f64; 3]> = feet.iter().map(|(_, x)| unit(sub(*x, apex.1))).collect();
    let a = (0..3).map(|d| (0..3).map(|i| t[i][d]).collect()).collect();
    let exact: Vec<(u32, f64)> = (1..=3).zip(solve_dense(a, load.map(|c| -c).to_vec())).collect();
    let err = worst_relative(&fe, &exact);
    eprintln!("space tripod: worst member force error {err:.2e}");
    assert!(err < 1e-10, "{fe:?} vs {exact:?}");
}

/// The statically indeterminate three-bar truss: a vertical bar of length L and two at ±θ
/// meeting at a loaded joint.
#[test]
fn three_bar_truss() {
    let (l, theta, p, ea) = (2.0, 35f64.to_radians(), 50.0, 3.0e4);
    for kind in [ElementKind::Rod2, ElementKind::Rod3] {
        let mut b = Build::statics();
        b.node(1, [-l * theta.tan(), l, 0.0]);
        b.node(2, [0.0, l, 0.0]);
        b.node(3, [l * theta.tan(), l, 0.0]);
        b.node(4, [0.0, 0.0, 0.0]);
        let mut members = Vec::new();
        for (i, top) in [1u32, 2, 3].into_iter().enumerate() {
            let mut nodes = vec![top, 4];
            if kind == ElementKind::Rod3 {
                // Middle nodes a third of the way along: still straight, not centred.
                let (x0, x1) = (b.m.node(top).unwrap().xyz, [0.0; 3]);
                let mid = 10 + i as u32;
                b.node(mid, [0, 1, 2].map(|k| x0[k] + (x1[k] - x0[k]) / 3.0));
                nodes.push(mid);
                // A straight member's middle node has no stiffness across it: hold it there.
                let across = unit([x0[1], -x0[0], 0.0]);
                let set = b.nodes(&format!("mid{mid}"), &[mid]);
                b.spring(&format!("mid{mid}"), set, None, Direction::Vector(across), 1.0);
            }
            members.push((1 + i as u32, nodes));
        }
        b.members("bars", kind, &members, 1.0, ea, None);
        let top = b.nodes("top", &[1, 2, 3]);
        b.fix("top", top);
        let all: Vec<u32> = b.m.nodes().keys().copied().collect();
        let plane = b.nodes("plane_nodes", &all);
        b.bc("plane", plane, ZERO_Z);
        b.force("P", 4, [0.0, -p, 0.0]);
        let r = b.solve();
        let c3 = theta.cos().powi(3);
        let vertical = p / (1.0 + 2.0 * c3);
        let inclined = p * theta.cos().powi(2) / (1.0 + 2.0 * c3);
        let fe = member_forces(&b.m, &r, |_| 1.0);
        let exact = [(1, inclined), (2, vertical), (3, inclined)];
        let err = worst_relative(&fe, &exact);
        let delta = r.displacement[r.nodes.binary_search(&4).unwrap()][1];
        let want = -vertical * l / ea;
        eprintln!("three-bar truss ({kind:?}): force error {err:.2e}, deflection {delta:.6e} vs {want:.6e}");
        assert!(err < 1e-10);
        assert!((delta - want).abs() < 1e-10 * want.abs());
    }
}

/// A rod3 bar with off-centre middle nodes under end tension carries σ = F/A exactly, at every
/// integration point and extrapolated to every node.
#[test]
fn rod3_constant_stress_patch() {
    let (area, young, f) = (2.0, 1000.0, 6.0);
    let mut b = Build::statics();
    let dir = unit([1.0, 2.0, 2.0]);
    let at = |s: f64| dir.map(|c| c * s);
    // Three elements of unequal length, middle nodes at 30%, 50% and 65% of each.
    let ends = [0.0, 1.0, 2.5, 3.2];
    let fractions = [0.3, 0.5, 0.65];
    for (i, s) in ends.iter().enumerate() {
        b.node(1 + i as u32, at(*s));
    }
    let mut members = Vec::new();
    for i in 0..3 {
        let s = ends[i] + fractions[i] * (ends[i + 1] - ends[i]);
        b.node(10 + i as u32, at(s));
        members.push((1 + i as u32, vec![1 + i as u32, 2 + i as u32, 10 + i as u32]));
    }
    b.members("bar", ElementKind::Rod3, &members, area, young, None);
    let start = b.nodes("start", &[1]);
    b.fix("start", start);
    // The other nodes slide along the bar only: hold them with springs across it.
    let others: Vec<u32> = b.m.nodes().keys().copied().filter(|l| *l != 1).collect();
    let rest = b.nodes("rest", &others);
    let across = unit([2.0, -1.0, 0.0]);
    let third = [dir[1] * across[2] - dir[2] * across[1], dir[2] * across[0] - dir[0] * across[2], dir[0] * across[1] - dir[1] * across[0]];
    b.spring("across1", rest.clone(), None, Direction::Vector(across), 1.0e6);
    b.spring("across2", rest, None, Direction::Vector(third), 1.0e6);
    b.force("pull", 4, dir.map(|c| c * f));
    let r = b.solve();
    let mut worst: f64 = 0.0;
    for e in &r.elements {
        for p in &e.points {
            worst = worst.max((axial(&p.stress, dir, false) - f / area).abs());
        }
        for s in &e.nodal {
            worst = worst.max((axial(s, dir, false) - f / area).abs());
        }
    }
    let tip = r.displacement[r.nodes.binary_search(&4).unwrap()];
    let stretch = tip.iter().zip(dir).map(|(u, d)| u * d).sum::<f64>();
    eprintln!("rod3 patch: stress error {worst:.2e}; tip stretch {stretch:.10} vs FL/EA {:.10}", f * 3.2 / (young * area));
    assert!(worst < 1e-10 * f / area);
    assert!((stretch - f * 3.2 / (young * area)).abs() < 1e-12);
}

/// A bar hanging under its own weight: u(z) = ρ g (L z − z²/2) / E from the support. Two-node
/// bars are exact at the nodes in one dimension; three-node bars reproduce the quadratic exactly.
#[test]
fn hanging_bar_under_gravity() {
    let (l, rho, g, young, area) = (6.0, 2.0, 9.81, 1.0e4, 0.5);
    for (kind, cells) in [(ElementKind::Rod2, 6usize), (ElementKind::Rod3, 3)] {
        let mut b = Build::statics();
        let h = l / cells as f64;
        for i in 0..=cells {
            b.node(1 + i as u32, [0.0, 0.0, -(i as f64) * h]);
        }
        let mut members = Vec::new();
        for i in 0..cells {
            let mut nodes = vec![1 + i as u32, 2 + i as u32];
            if kind == ElementKind::Rod3 {
                b.node(100 + i as u32, [0.0, 0.0, -(i as f64 + 0.5) * h]);
                nodes.push(100 + i as u32);
            }
            members.push((1 + i as u32, nodes));
        }
        b.members("bar", kind, &members, area, young, Some(rho));
        let top = b.nodes("top", &[1]);
        b.fix("top", top);
        let all: Vec<u32> = b.m.nodes().keys().copied().collect();
        let line = b.nodes("line", &all);
        b.bc("line", line, [Some(0.0), Some(0.0), None, None, None, None]);
        b.load("weight", Load::Gravity { acceleration: [0.0, 0.0, -g], set: None });
        let r = b.solve();
        let mut worst: f64 = 0.0;
        for (label, u) in r.nodes.iter().zip(&r.displacement) {
            let z = -b.m.node(*label).unwrap().xyz[2];
            let exact = -rho * g * (l * z - z * z / 2.0) / young;
            worst = worst.max((u[2] - exact).abs());
        }
        let tip = rho * g * l * l / (2.0 * young);
        eprintln!("hanging bar ({kind:?}): nodal displacement error {:.2e} of the tip's", worst / tip);
        assert!(worst < 1e-12 * tip.max(1.0), "{worst}");
        assert!((r.reaction_total[2] - rho * g * area * l).abs() < 1e-10);
    }
}

/// One node, a spring to the ground in X and a mass: ω = √(k/m) (consistent and lumped mass
/// alike); a static force gives u = F/k.
#[test]
fn spring_mass_natural_frequency_and_static_stretch() {
    let (k, m) = (400.0, 2.5);
    for mass in [MassKind::Consistent, MassKind::Lumped] {
        let mut b = Build::modes(1);
        b.node(1, [1.0, 2.0, 3.0]);
        let n = b.nodes("n", &[1]);
        b.spring("k", n.clone(), None, Direction::Dof(1), k);
        b.mass("m", n.clone(), Some(Mass::Isotropic(m)), None);
        b.bc("yz", n, [None, Some(0.0), Some(0.0), None, None, None]);
        let r = b.modes_of(mass);
        let f = (k / m).sqrt() / (2.0 * PI);
        let err = (r.modes[0].frequency - f).abs() / f;
        eprintln!("spring–mass ({mass:?}): f {:.12} vs √(k/m)/2π {f:.12}, error {err:.1e}", r.modes[0].frequency);
        assert!(err < 1e-12);
        assert!((r.total_mass - m).abs() < 1e-12);
        assert!((r.modes[0].effective_mass[0] - m).abs() < 1e-9);
    }
    let mut b = Build::statics();
    b.node(1, [0.0; 3]);
    let n = b.nodes("n", &[1]);
    b.spring("k", n.clone(), None, Direction::Dof(1), k);
    b.bc("yz", n, [None, Some(0.0), Some(0.0), None, None, None]);
    b.force("F", 1, [8.0, 0.0, 0.0]);
    let r = b.solve();
    assert!((r.displacement[0][0] - 8.0 / k).abs() < 1e-15);
    // The spring carries the load to the ground: no support reaction, and no imbalance.
    assert_eq!(r.reaction_total[0], 0.0);
    assert!(r.warnings.iter().all(|w| !w.contains("differ")), "{:?}", r.warnings);
}

/// Springs to the ground along three orthogonal skew vectors: frequencies √(kᵢ/m) with the mode
/// shapes along the vectors.
#[test]
fn springs_in_skew_directions() {
    let m = 3.0;
    let a = unit([1.0, 1.0, 0.5]);
    let b2 = unit([1.0, -1.0, 0.0]);
    let c = [a[1] * b2[2] - a[2] * b2[1], a[2] * b2[0] - a[0] * b2[2], a[0] * b2[1] - a[1] * b2[0]];
    let ks = [30.0, 120.0, 270.0];
    let mut b = Build::modes(3);
    b.node(7, [0.0; 3]);
    let n = b.nodes("n", &[7]);
    for (i, (v, k)) in [a, b2, c].iter().zip(ks).enumerate() {
        b.spring(&format!("k{i}"), n.clone(), None, Direction::Vector(v.map(|x| 2.0 * x)), k);
    }
    b.mass("m", n, Some(Mass::Isotropic(m)), None);
    let r = b.modes_of(MassKind::Consistent);
    for ((mode, k), v) in r.modes.iter().zip(ks).zip([a, b2, c]) {
        let f = (k / m).sqrt() / (2.0 * PI);
        let s = mode.shape[0];
        let along = (s[0] * v[0] + s[1] * v[1] + s[2] * v[2]).abs() / (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
        eprintln!("skew springs: f {:.10} vs {f:.10}; shape along the spring {along:.12}", mode.frequency);
        assert!((mode.frequency - f).abs() < 1e-10 * f);
        assert!((along - 1.0).abs() < 1e-9);
    }
}

/// Ground — k₁ — m₁ — k₂ — m₂, the two-point spring acting along the line of its nodes:
/// ω² are the roots of m₁m₂ω⁴ − (k₁m₂ + k₂m₁ + k₂m₂)ω² + k₁k₂ = 0.
#[test]
fn two_mass_chain() {
    let (k1, k2, m1, m2) = (100.0, 60.0, 2.0, 0.5);
    let mut b = Build::modes(2);
    let dir = unit([3.0, 0.0, 4.0]);
    b.node(1, [0.0; 3]);
    b.node(2, dir.map(|c| 2.0 * c));
    let n1 = b.nodes("n1", &[1]);
    let n2 = b.nodes("n2", &[2]);
    b.spring("k1", n1.clone(), None, Direction::Vector(dir), k1);
    b.spring("k2", n1.clone(), Some(n2.clone()), Direction::Line, k2);
    b.mass("m1", n1, Some(Mass::Isotropic(m1)), None);
    b.mass("m2", n2, Some(Mass::Directional([m2, 7.0, m2])), None);
    // Only the motion along the line: hold the others with stiff springs across it.
    let both = b.nodes("both", &[1, 2]);
    b.spring("hold_y", both.clone(), None, Direction::Dof(2), 1.0e9);
    b.spring("hold_n", both, None, Direction::Vector([dir[2], 0.0, -dir[0]]), 1.0e9);
    let r = b.modes_of(MassKind::Lumped);
    let (pa, pb, pc) = (m1 * m2, -(k1 * m2 + k2 * m1 + k2 * m2), k1 * k2);
    let disc = (pb * pb - 4.0 * pa * pc).sqrt();
    let exact = [(-pb - disc) / (2.0 * pa), (-pb + disc) / (2.0 * pa)];
    for (mode, w2) in r.modes.iter().zip(exact) {
        let err = (mode.eigenvalue - w2).abs() / w2;
        eprintln!("two-mass chain: ω² {:.10} vs {w2:.10}, error {err:.1e}", mode.eigenvalue);
        assert!(err < 1e-8, "the holding springs perturb at 1e-9");
    }
}

/// A massless bar with a tip mass M: ω = √(EA/(L M)); the bar's own mass ρAL adds the
/// Rayleigh-quotient correction ω² ≈ EA/(L(M + ρAL/3)) for a light bar.
#[test]
fn bar_with_a_tip_mass() {
    let (l, young, area, tip) = (2.0, 2.0e5, 1e-2, 50.0);
    for kind in [ElementKind::Rod2, ElementKind::Rod3] {
        let mut b = Build::modes(1);
        b.node(1, [0.0; 3]);
        b.node(2, [l, 0.0, 0.0]);
        let mut nodes = vec![1, 2];
        if kind == ElementKind::Rod3 {
            b.node(3, [l / 2.0, 0.0, 0.0]);
            nodes.push(3);
        }
        b.members("bar", kind, &[(1, nodes)], area, young, Some(1e-9));
        let root = b.nodes("root", &[1]);
        b.fix("root", root);
        let all: Vec<u32> = b.m.nodes().keys().copied().collect();
        let line = b.nodes("line", &all);
        b.bc("line", line, [None, Some(0.0), Some(0.0), None, None, None]);
        let end = b.nodes("end", &[2]);
        b.mass("tip", end, Some(Mass::Isotropic(tip)), None);
        let r = b.modes_of(MassKind::Consistent);
        let f = (young * area / (l * tip)).sqrt() / (2.0 * PI);
        let err = (r.modes[0].frequency - f).abs() / f;
        eprintln!("bar with tip mass ({kind:?}): f {:.10} vs {f:.10}, error {err:.1e}", r.modes[0].frequency);
        assert!(err < 1e-9);
    }
}

/// A torsional pendulum: a reference point with polar inertia J on a torsional spring k_t to
/// the ground; ω = √(k_t/J). The other directions are held.
#[test]
fn torsional_pendulum() {
    let (kt, j) = (5.0e3, 0.8);
    let mut b = Build::modes(1);
    b.m.add_reference_point("Disc", ReferencePoint { xyz: [1.0, 1.0, 0.0], label: 0 }).unwrap();
    let disc = Target::ReferencePoint("Disc".into());
    b.spring("shaft", disc.clone(), None, Direction::Dof(6), kt);
    // Products of inertia about the other axes do not couple into the held rotations.
    b.mass("disc", disc.clone(), Some(Mass::Isotropic(4.0)), Some([0.5, 0.5, j, 0.0, 0.1, 0.1]));
    b.bc("held", disc, [Some(0.0), Some(0.0), Some(0.0), Some(0.0), Some(0.0), None]);
    let r = b.modes_of(MassKind::Lumped);
    let f = (kt / j).sqrt() / (2.0 * PI);
    let err = (r.modes[0].frequency - f).abs() / f;
    eprintln!("torsional pendulum: f {:.12} vs √(k_t/J)/2π {f:.12}, error {err:.1e}", r.modes[0].frequency);
    assert!(err < 1e-12);
    assert_eq!(r.reference_points.len(), 1);
    assert!(r.modes[0].point_shape[0][5].abs() > 0.0);
}

/// Gravity acts on point masses (whole-model gravity); dashpots are ignored with a warning in
/// a static step; the runs are bit-for-bit repeatable.
#[test]
fn point_mass_gravity_dashpots_and_determinism() {
    let mut b = Build::statics();
    b.node(1, [0.0; 3]);
    b.node(2, [0.0, 0.0, -1.0]);
    b.members("rod", ElementKind::Rod2, &[(1, vec![1, 2])], 1.0, 1.0e3, Some(0.0 + 1e-6));
    let top = b.nodes("top", &[1]);
    b.fix("top", top);
    let bottom = b.nodes("bottom", &[2]);
    b.bc("line", bottom.clone(), [Some(0.0), Some(0.0), None, None, None, None]);
    b.mass("weight", bottom.clone(), Some(Mass::Isotropic(3.0)), None);
    b.m.add_spring("damper", Spring::Dashpot { first: bottom, second: None, direction: Direction::Dof(3), coefficient: 2.0 }).unwrap();
    b.load("g", Load::Gravity { acceleration: [0.0, 0.0, -10.0], set: None });
    let r = b.solve();
    let total = 3.0 * 10.0 + 1e-6 * 10.0;
    assert!((r.reaction_total[2] - total).abs() < 1e-10, "{:?}", r.reaction_total);
    assert!((r.displacement[r.nodes.binary_search(&2).unwrap()][2] + total * 1.0 / 1.0e3).abs() < 1e-6);
    assert!(r.warnings.iter().any(|w| w.contains("dashpots `damper`")), "{:?}", r.warnings);
    assert_eq!(r, b.solve());
}

/// What the solver cannot do is an error that says so.
#[test]
fn errors_name_the_problem() {
    // A point mass with nothing holding it.
    let mut b = Build::modes(1);
    b.node(1, [0.0; 3]);
    let n = b.nodes("n", &[1]);
    b.mass("m", n, Some(Mass::Isotropic(1.0)), None);
    // A free particle has no stiffness at all to shift against.
    assert!(solve_frequency(&b.m, None, &FrequencyOptions::default(), &Limits::default()).is_err());
    // A truss of an orthotropic material.
    let mut b = Build::statics();
    b.node(1, [0.0; 3]);
    b.node(2, [1.0, 0.0, 0.0]);
    b.m.set_element(1, ElementKind::Rod2, vec![1, 2]).unwrap();
    b.m.add_set("rod", Set { elements: vec![1], ..Set::default() }).unwrap();
    b.m.add_material(
        "Ortho",
        Material {
            elastic: Some(Elastic::EngineeringConstants { young: [1.0; 3], poisson: [0.0; 3], shear: [1.0; 3] }),
            density: None,
            ..Material::default()
        },
    )
    .unwrap();
    b.m.add_section(
        "rod",
        Section {
            orientation: None,
            set: "rod".into(),
            material: "Ortho".into(),
            kind: SectionKind::Rod { area: 1.0 },
            formulation: Formulation::default(),
        },
    )
    .unwrap();
    let e = solve_static(&b.m, None, &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("orthotropic"), "{e}");
}

/// A 2D (plane stress) plate stiffened by trusses along its loaded edges: the truss takes the
/// model's two in-plane unknowns per node, and the two act as parallel springs,
/// u = F L / (E t w + E_b A) (exact: uniform tension, ν = 0).
#[test]
fn plane_truss_stiffens_a_plane_stress_plate() {
    let (l, w, t, young, area, young_bar, force) = (4.0, 1.0, 0.2, 1000.0, 0.05, 2.0e4, 30.0);
    for (kind, mid) in [(ElementKind::Rod2, false), (ElementKind::Rod3, true)] {
        let mut b = Build::statics();
        for (label, xyz) in [(1, [0.0, 0.0, 0.0]), (2, [l, 0.0, 0.0]), (3, [l, w, 0.0]), (4, [0.0, w, 0.0])] {
            b.node(label, xyz);
        }
        b.m.set_element(1, ElementKind::Quad4, vec![1, 2, 3, 4]).unwrap();
        b.m.add_set("plate", Set { elements: vec![1], ..Set::default() }).unwrap();
        b.m.add_material("Plate", Material { elastic: Some(Elastic::Isotropic { young, poisson: 0.0 }), density: None, ..Material::default() })
            .unwrap();
        b.m.add_section(
            "plate",
            Section {
                orientation: None,
                set: "plate".into(),
                material: "Plate".into(),
                kind: SectionKind::PlaneStress { thickness: t },
                formulation: Formulation::default(),
            },
        )
        .unwrap();
        // Two bars of half the area, along the bottom and top edges, keep the tension uniform.
        let mut members = vec![(2, vec![1, 2]), (3, vec![4, 3])];
        if mid {
            b.node(5, [l / 2.0, 0.0, 0.0]);
            b.node(6, [l / 2.0, w, 0.0]);
            members[0].1.push(5);
            members[1].1.push(6);
        }
        b.members("bar", kind, &members, area / 2.0, young_bar, None);
        let left = b.nodes("left", &[1, 4]);
        b.bc("left", left, [Some(0.0), None, None, None, None, None]);
        let corner = b.nodes("corner", &[1]);
        b.bc("corner", corner, [None, Some(0.0), None, None, None, None]);
        if mid {
            let middle = b.nodes("middle", &[5, 6]);
            b.bc("middle", middle, [None, Some(0.0), None, None, None, None]);
        }
        let right = b.nodes("right", &[2, 3]);
        b.load("pull", Load::Force { target: right, vector: [force, 0.0, 0.0], distribution: Distribution::PerNode });
        let r = b.solve();
        let u = r.displacement[r.nodes.binary_search(&2).unwrap()][0];
        let exact = 2.0 * force * l / (young * t * w + young_bar * area);
        eprintln!("plane truss ({kind:?}): u {u:.12} vs {exact:.12}; {} unknowns", r.stats.dofs);
        assert!((u - exact).abs() < 1e-12 * exact);
        assert_eq!(r.stats.dofs, 2 * b.m.nodes().len(), "two unknowns per node in a 2D model");
    }
}
