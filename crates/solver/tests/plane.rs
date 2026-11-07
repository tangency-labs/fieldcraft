//! The analytic verification suite of the 2D elements (plane stress, plane strain,
//! axisymmetric): patch tests, a traction patch, a cantilever, Kirsch's plate, Lamé's cylinder
//! and rotating discs and cylinders against closed-form solutions. Run with `--nocapture` to see
//! the errors.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use fieldcraft_model::{
    BoundaryCondition, Elastic, ElementEdge, ElementKind, Formulation, Integration, Load, Material, Model, Procedure, Section, SectionKind, Set,
    Step, Symmetry, Target,
};
use fieldcraft_solver::{FrequencyOptions, Limits, StaticResult, solve_frequency, solve_static};

/// A structured 2D mesh of `n` cells on [0, 1]², mapped by `map`; nodes on the doubled lattice so
/// quadratic elements take their mid-side nodes from the map.
struct Grid {
    model: Model,
    param: BTreeMap<u32, [f64; 2]>,
}

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

fn grid(n: [usize; 2], kind: ElementKind, map: impl Fn([f64; 2]) -> [f64; 2]) -> Grid {
    let size = [2 * n[0] + 1, 2 * n[1] + 1];
    let mut model = Model::default();
    let mut param = BTreeMap::new();
    let mut node = |model: &mut Model, i: [usize; 2]| -> u32 {
        let l = (1 + i[0] + size[0] * i[1]) as u32;
        param.entry(l).or_insert_with(|| {
            let p = [i[0] as f64 / (2 * n[0]) as f64, i[1] as f64 / (2 * n[1]) as f64];
            let x = map(p);
            model.set_node(l, [x[0], x[1], 0.0]).unwrap();
            p
        });
        l
    };
    let mut next = 1u32;
    for cy in 0..n[1] {
        for cx in 0..n[0] {
            let at = |d: [usize; 2]| [2 * cx + d[0], 2 * cy + d[1]];
            let polygons: Vec<Vec<[usize; 2]>> = match kind {
                ElementKind::Quad4 | ElementKind::Quad8 => vec![vec![[0, 0], [2, 0], [2, 2], [0, 2]]],
                _ => vec![vec![[0, 0], [2, 0], [2, 2]], vec![[0, 0], [2, 2], [0, 2]]],
            };
            for corners in polygons {
                let mut nodes: Vec<u32> = corners.iter().map(|c| node(&mut model, at(*c))).collect();
                if kind.is_quadratic() {
                    for e in kind.edges() {
                        let (a, b) = (corners[e[0]], corners[e[1]]);
                        nodes.push(node(&mut model, at([(a[0] + b[0]) / 2, (a[1] + b[1]) / 2])));
                    }
                }
                model.set_element(next, kind, nodes).unwrap();
                next += 1;
            }
        }
    }
    Grid { model, param }
}

impl Grid {
    fn nodes_where(&self, f: impl Fn([f64; 2]) -> bool) -> Vec<u32> {
        self.param.iter().filter(|(_, p)| f(**p)).map(|(l, _)| *l).collect()
    }

    fn edges_where(&self, f: impl Fn([f64; 2]) -> bool) -> Vec<ElementEdge> {
        let mut out = Vec::new();
        for (label, e) in self.model.elements() {
            for (i, edge) in e.kind.edges().iter().enumerate() {
                if edge.iter().all(|k| f(self.param[&e.nodes[*k]])) {
                    out.push(ElementEdge { element: *label, edge: i as u8 + 1 });
                }
            }
        }
        out
    }

    fn define(&mut self, kind: SectionKind, young: f64, nu: f64, density: Option<f64>, formulation: Formulation, procedure: Procedure) {
        let all: Vec<u32> = self.model.elements().keys().copied().collect();
        self.model.add_set("all", Set { elements: all, ..Set::default() }).unwrap();
        self.model
            .add_material("Mat", Material { elastic: Some(Elastic::Isotropic { young, poisson: nu }), density, ..Material::default() })
            .unwrap();
        self.model.add_section("plane", Section { orientation: None, set: "all".into(), material: "Mat".into(), kind, formulation }).unwrap();
        self.model.add_step(Step::new("Step", procedure), None).unwrap();
    }

    fn bc(&mut self, name: &str, bc: BoundaryCondition) {
        self.model.add_boundary_condition(name, bc).unwrap();
        self.model.set_boundary_condition_active("Step", name, true).unwrap();
    }

    fn load(&mut self, name: &str, load: Load) {
        self.model.add_load(name, load).unwrap();
        self.model.set_load_active("Step", name, true).unwrap();
    }

    fn symmetry(&mut self, name: &str, plane: Symmetry, f: impl Fn([f64; 2]) -> bool) {
        let nodes = self.nodes_where(f);
        self.model.add_set(name, Set { nodes, ..Set::default() }).unwrap();
        self.bc(name, BoundaryCondition::Symmetry { target: Target::Set(name.into()), plane });
    }

    fn fixed(&mut self, name: &str, f: impl Fn([f64; 2]) -> bool) {
        let nodes = self.nodes_where(f);
        self.model.add_set(name, Set { nodes, ..Set::default() }).unwrap();
        self.bc(name, BoundaryCondition::Fixed { target: Target::Set(name.into()) });
    }

    fn edge_set(&mut self, name: &str, f: impl Fn([f64; 2]) -> bool) {
        let edges = self.edges_where(f);
        self.model.add_set(name, Set { edges, ..Set::default() }).unwrap();
    }

    fn solve(&self) -> StaticResult {
        solve_static(&self.model, Some("Step"), &Limits::default()).unwrap()
    }
}

fn at(r: &StaticResult, label: u32) -> usize {
    r.nodes.binary_search(&label).unwrap()
}

fn full() -> Formulation {
    Formulation::default()
}

fn reduced() -> Formulation {
    Formulation { integration: Integration::Reduced, ..Formulation::default() }
}

fn incompatible() -> Formulation {
    Formulation { incompatible_modes: true, ..Formulation::default() }
}

fn variants() -> Vec<(ElementKind, Formulation, &'static str)> {
    vec![
        (ElementKind::Tri3, full(), "tri3"),
        (ElementKind::Tri6, full(), "tri6"),
        (ElementKind::Quad4, full(), "quad4"),
        (ElementKind::Quad4, incompatible(), "quad4 incompatible"),
        (ElementKind::Quad8, full(), "quad8"),
        (ElementKind::Quad8, reduced(), "quad8 reduced"),
    ]
}

const STRESS: SectionKind = SectionKind::PlaneStress { thickness: 0.5 };
const STRAIN: SectionKind = SectionKind::PlaneStrain { thickness: 1.0 };
const AXI: SectionKind = SectionKind::Axisymmetric {};

/// Displacement patch test on a distorted 2 × 2 patch over [1, 2] × [0, 1]: a linear
/// displacement prescribed on the boundary is reproduced inside, with exactly constant strain
/// and stress (plane: a general linear field; axisymmetric: u_r = a r, u_z = b z, the
/// constant-strain state a body of revolution can be in without body forces; a constant shear
/// γ_rz would need an axial body force τ/r).
#[test]
fn patch_test_every_element_and_theory() {
    let (young, nu) = (1000.0, 0.25);
    let map = |p: [f64; 2]| {
        let b = (PI * p[0]).sin() * (PI * p[1]).sin();
        [1.0 + p[0] + 0.15 * b, p[1] - 0.1 * b]
    };
    let lambda = young * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = young / (2.0 * (1.0 + nu));
    for theory in [STRESS, STRAIN, AXI] {
        // u = field(x); strain (xx, yy, zz, xy).
        let (field, strain): (fn([f64; 3]) -> [f64; 2], [f64; 4]) = if theory == AXI {
            (|x| [1e-3 * x[0], -5e-4 * x[1]], [1e-3, -5e-4, 1e-3, 0.0])
        } else {
            (|x| [1e-3 * x[0] + 2e-4 * x[1], 5e-4 * x[0] - 2e-3 * x[1]], [1e-3, -2e-3, 0.0, 7e-4])
        };
        let mut e = strain;
        if let SectionKind::PlaneStress { .. } = theory {
            e[2] = -nu / (1.0 - nu) * (e[0] + e[1]);
        }
        let tr = e[0] + e[1] + e[2];
        let stress = [lambda * tr + 2.0 * mu * e[0], lambda * tr + 2.0 * mu * e[1], lambda * tr + 2.0 * mu * e[2], mu * e[3]];
        for (kind, formulation, name) in variants() {
            let mut g = grid([2, 2], kind, map);
            if theory == AXI {
                // Straight sides: on a curved quadratic element the axisymmetric integrand
                // (r ∂N/∂r + N) det J is rational, so no Gauss rule integrates it exactly and the
                // patch test holds only to the rule's accuracy.
                straighten(&mut g.model);
            }
            g.define(theory.clone(), young, nu, None, formulation, Procedure::Static {});
            let boundary = g.nodes_where(|p| near(p[0], 0.0) || near(p[0], 1.0) || near(p[1], 0.0) || near(p[1], 1.0));
            for (i, label) in boundary.into_iter().enumerate() {
                let u = field(g.model.node(label).unwrap().xyz);
                let set = format!("b{i}");
                g.model.add_set(&set, Set { nodes: vec![label], ..Set::default() }).unwrap();
                g.bc(
                    &set,
                    BoundaryCondition::Displacement { target: Target::Set(set.clone()), values: [Some(u[0]), Some(u[1]), None, None, None, None] },
                );
            }
            let r = g.solve();
            let (mut eu, mut es): (f64, f64) = (0.0, 0.0);
            for (label, u) in r.nodes.iter().zip(&r.displacement) {
                let want = field(g.model.node(*label).unwrap().xyz);
                eu = eu.max((u[0] - want[0]).abs().max((u[1] - want[1]).abs()));
                assert_eq!(u[2], 0.0);
            }
            for el in &r.elements {
                for p in &el.points {
                    for (got, want) in p.stress.iter().zip(&stress) {
                        es = es.max((got - want).abs() / young);
                    }
                }
            }
            eprintln!("patch {} {name}: displacement error {eu:.1e}, stress error {es:.1e}", theory.name());
            assert!(eu < 1e-15 && es < 1e-14, "{} {name}: {eu:e} {es:e}", theory.name());
        }
    }
}

/// Put every mid-side node at the middle of its edge's ends.
fn straighten(model: &mut Model) {
    let elements: Vec<(ElementKind, Vec<u32>)> = model.elements().values().map(|e| (e.kind, e.nodes.clone())).collect();
    for (kind, nodes) in elements {
        for e in kind.edges() {
            if let [a, b, m] = **e {
                let (pa, pb) = (model.node(nodes[a]).unwrap().xyz, model.node(nodes[b]).unwrap().xyz);
                model.set_node(nodes[m], [0, 1, 2].map(|k| (pa[k] + pb[k]) / 2.0)).unwrap();
            }
        }
    }
}

/// A uniform traction on the free edge of a plate (plane stress), a slab (plane strain) and a
/// tube's end (axisymmetric), held by symmetry: the exact stress and displacement everywhere.
#[test]
fn uniform_traction_is_exact() {
    let (young, nu, s0) = (1000.0, 0.3, 2.0);
    for theory in [STRESS, STRAIN, AXI] {
        for (kind, formulation, name) in variants() {
            // Axisymmetric: a tube r ∈ [1, 2], pulled along its axis (Y); plane: pulled along X.
            let mut g = grid([3, 2], kind, |p| [1.0 + p[0], 2.0 * p[1]]);
            g.define(theory.clone(), young, nu, None, formulation, Procedure::Static {});
            if theory == AXI {
                g.symmetry("bottom", Symmetry::Y, |p| near(p[1], 0.0));
                g.edge_set("top", |p| near(p[1], 1.0));
                g.load("pull", Load::Traction { set: "top".into(), vector: [0.0, s0, 0.0] });
            } else {
                g.symmetry("left", Symmetry::X, |p| near(p[0], 0.0));
                g.symmetry("bottom", Symmetry::Y, |p| near(p[1], 0.0));
                g.edge_set("right", |p| near(p[0], 1.0));
                g.load("pull", Load::Pressure { set: "right".into(), magnitude: -s0 });
            }
            let r = g.solve();
            // Strains: along the load s/E' and across −ν' s/E'.
            let (along, across) = match theory {
                SectionKind::PlaneStrain { .. } => (s0 * (1.0 - nu * nu) / young, -nu * (1.0 + nu) * s0 / young),
                _ => (s0 / young, -nu * s0 / young),
            };
            let mut eu: f64 = 0.0;
            for (label, u) in r.nodes.iter().zip(&r.displacement) {
                let x = g.model.node(*label).unwrap().xyz;
                let want = if theory == AXI { [across * x[0], along * x[1]] } else { [along * (x[0] - 1.0), across * x[1]] };
                eu = eu.max((u[0] - want[0]).abs().max((u[1] - want[1]).abs()));
            }
            let total: f64 = r.reaction_total.iter().map(|c| c.abs()).sum();
            let area = if theory == AXI {
                PI * (4.0 - 1.0)
            } else {
                2.0 * if let SectionKind::PlaneStress { thickness } = theory { thickness } else { 1.0 }
            };
            eprintln!("traction {} {name}: displacement error {eu:.1e}, reaction {total:.6} vs {:.6}", theory.name(), s0 * area);
            assert!(eu < 1e-13, "{} {name}: {eu:e}", theory.name());
            assert!((total - s0 * area).abs() < 1e-11 * s0 * area, "{} {name}: {total} vs {}", theory.name(), s0 * area);
        }
    }
}

/// Cantilever L × h of thickness t in plane stress, clamped at x = 0, a shear traction on the free
/// end; Timoshenko's deflection δ = P L³/(3 E I) + P L/(κ G A), κ = 5/6, with ν = 0.
fn cantilever_tip(kind: ElementKind, formulation: Formulation, n: [usize; 2]) -> f64 {
    let (l, h, t, young, p) = (10.0, 1.0, 1.0, 1000.0, 1.0);
    let mut g = grid(n, kind, |q| [l * q[0], h * q[1]]);
    g.define(SectionKind::PlaneStress { thickness: t }, young, 0.0, None, formulation, Procedure::Static {});
    g.fixed("root", |q| near(q[0], 0.0));
    g.edge_set("end", |q| near(q[0], 1.0));
    g.load("shear", Load::Traction { set: "end".into(), vector: [0.0, -p / (h * t), 0.0] });
    let r = g.solve();
    let end = g.nodes_where(|q| near(q[0], 1.0));
    -end.iter().map(|l| r.displacement[at(&r, *l)][1]).sum::<f64>() / end.len() as f64
}

#[test]
fn cantilever_in_plane_stress() {
    let (l, h, young, p) = (10.0f64, 1.0f64, 1000.0, 1.0);
    let exact = p * l.powi(3) / (3.0 * young * h.powi(3) / 12.0) + p * l / (5.0 / 6.0 * young / 2.0 * h);
    let err = |d: f64| (d - exact).abs() / exact;
    let study = |kind, f, meshes: &[[usize; 2]], name: &str| -> Vec<f64> {
        meshes
            .iter()
            .map(|n| {
                let d = cantilever_tip(kind, f, *n);
                eprintln!("cantilever {name} {n:?}: {d:.5} vs {exact:.5} (error {:.2}%)", 100.0 * err(d));
                err(d)
            })
            .collect()
    };
    let q8 = study(ElementKind::Quad8, full(), &[[10, 1], [20, 2]], "quad8");
    assert!(q8[1] < 0.005, "{q8:?}");
    let q8r = study(ElementKind::Quad8, reduced(), &[[10, 2]], "quad8 reduced");
    assert!(q8r[0] < 0.005, "{q8r:?}");
    let im = study(ElementKind::Quad4, incompatible(), &[[10, 1], [20, 2]], "quad4 incompatible");
    assert!(im[0] < 0.03 && im[1] < 0.01, "{im:?}");
    let t6 = study(ElementKind::Tri6, full(), &[[10, 2]], "tri6");
    assert!(t6[0] < 0.01, "{t6:?}");
    let meshes = [[5, 1], [10, 2], [20, 4], [40, 8]];
    for (kind, name) in [(ElementKind::Quad4, "quad4"), (ElementKind::Tri3, "tri3")] {
        let e = study(kind, full(), &meshes, name);
        let rates: Vec<f64> = e.windows(2).map(|w| (w[0] / w[1]).log2()).collect();
        eprintln!("  {name} convergence rates {rates:.2?}");
        assert!(e.windows(2).all(|w| w[1] < w[0]), "{name} converges monotonically");
        assert!(*rates.last().unwrap() > 1.5, "{name} approaches second order: {rates:?}");
    }
}

/// Kirsch's plate with a hole of radius a under remote tension, plane stress, quarter model of
/// half-width 20a: the hoop stress at the top of the hole tends to 3σ₀.
#[test]
fn kirsch_plate_with_hole() {
    let (a, w, s0) = (1.0f64, 20.0f64, 1.0);
    let ratio: f64 = 40.0;
    let map = |q: [f64; 2]| {
        let th = PI / 2.0 * q[1];
        let outer = if th <= PI / 4.0 { [w, w * th.tan()] } else { [w / th.tan(), w] };
        let s = (ratio.powf(q[0]) - 1.0) / (ratio - 1.0);
        let inner = [a * th.cos(), a * th.sin()];
        [inner[0] + s * (outer[0] - inner[0]), inner[1] + s * (outer[1] - inner[1])]
    };
    let mut results = Vec::new();
    for (kind, formulation, n, name) in [
        (ElementKind::Quad8, full(), [6, 8], "quad8 coarse"),
        (ElementKind::Quad8, full(), [12, 16], "quad8 fine"),
        (ElementKind::Tri6, full(), [12, 16], "tri6"),
        (ElementKind::Quad4, incompatible(), [24, 32], "quad4 incompatible"),
    ] {
        let mut g = grid(n, kind, map);
        g.define(SectionKind::PlaneStress { thickness: 0.1 }, 1000.0, 0.3, None, formulation, Procedure::Static {});
        g.symmetry("x0", Symmetry::X, |q| near(q[1], 1.0));
        g.symmetry("y0", Symmetry::Y, |q| near(q[1], 0.0));
        g.edge_set("right", |q| near(q[0], 1.0) && q[1] <= 0.5 + 1e-9);
        g.load("pull", Load::Traction { set: "right".into(), vector: [s0, 0.0, 0.0] });
        let r = g.solve();
        let top = g.nodes_where(|q| near(q[0], 0.0) && near(q[1], 1.0));
        let kt = top.iter().map(|l| r.stress[at(&r, *l)][0]).sum::<f64>() / top.len() as f64 / s0;
        eprintln!("kirsch {name} {n:?}: Kt = {kt:.4} (error vs 3: {:.2}%)", 100.0 * (kt - 3.0).abs() / 3.0);
        assert!(r.stress.iter().all(|s| s[2] == 0.0), "plane stress has no σzz");
        results.push(kt);
    }
    assert!((results[1] - 3.0).abs() < (results[0] - 3.0).abs(), "{results:?}");
    assert!((results[1] - 3.0).abs() / 3.0 < 0.02, "{results:?}");
    assert!((results[2] - 3.0).abs() / 3.0 < 0.03, "{results:?}");
    assert!((results[3] - 3.0).abs() / 3.0 < 0.04, "{results:?}");
}

/// Lamé's thick cylinder a ≤ r ≤ b under internal pressure p in plane strain:
/// σ_r = A − B/r², σ_θ = A + B/r², σ_z = 2νA, u_r = (1 + ν)/E ((1 − 2ν) A r + B/r),
/// A = p a²/(b² − a²), B = p a² b²/(b² − a²). Solved as a quarter section in plane strain and as
/// an axisymmetric slice held axially on both faces.
#[test]
fn lame_thick_cylinder() {
    let (a, b, p, young, nu) = (1.0f64, 2.0f64, 10.0, 1000.0, 0.3);
    let big_a = p * a * a / (b * b - a * a);
    let big_b = p * a * a * b * b / (b * b - a * a);
    let ur = |r: f64| (1.0 + nu) / young * ((1.0 - 2.0 * nu) * big_a * r + big_b / r);
    let sigma_t = |r: f64| big_a + big_b / (r * r);
    for (kind, formulation, n, name) in [
        (ElementKind::Quad8, full(), [6, 8], "quad8"),
        (ElementKind::Quad4, incompatible(), [12, 16], "quad4 incompatible"),
        (ElementKind::Tri6, full(), [6, 8], "tri6"),
    ] {
        // Plane strain quarter: (r, θ) mapped.
        let mut g = grid(n, kind, |q| {
            let (r, th) = (a + (b - a) * q[0], PI / 2.0 * q[1]);
            [r * th.cos(), r * th.sin()]
        });
        g.define(STRAIN, young, nu, None, formulation, Procedure::Static {});
        g.symmetry("y0", Symmetry::Y, |q| near(q[1], 0.0));
        g.symmetry("x0", Symmetry::X, |q| near(q[1], 1.0));
        g.edge_set("bore", |q| near(q[0], 0.0));
        g.load("p", Load::Pressure { set: "bore".into(), magnitude: p });
        let r = g.solve();
        let mut eu: f64 = 0.0;
        for (label, u) in r.nodes.iter().zip(&r.displacement) {
            let x = g.model.node(*label).unwrap().xyz;
            let rad = x[0].hypot(x[1]);
            eu = eu.max(((u[0] * x[0] + u[1] * x[1]) / rad - ur(rad)).abs() / ur(a));
        }
        let bore = g.nodes_where(|q| near(q[0], 0.0) && near(q[1], 0.5));
        let st = bore
            .iter()
            .map(|l| {
                let s = r.stress[at(&r, *l)];
                let x = g.model.node(*l).unwrap().xyz;
                let (c, sn) = (x[0] / x[0].hypot(x[1]), x[1] / x[0].hypot(x[1]));
                s[0] * sn * sn - 2.0 * s[3] * sn * c + s[1] * c * c
            })
            .sum::<f64>()
            / bore.len() as f64;
        let es = (st - sigma_t(a)).abs() / sigma_t(a);
        eprintln!("lamé plane strain {name} {n:?}: u_r error {:.3}%, bore hoop {st:.4} vs {:.4} ({:.2}%)", 100.0 * eu, sigma_t(a), 100.0 * es);
        assert!(eu < 0.01 && es < 0.03, "{name}");

        // Axisymmetric slice r ∈ [a, b], 0 ≤ y ≤ 0.5, held axially on both faces.
        let mut g = grid([n[1], 2], kind, |q| [a + (b - a) * q[0], 0.5 * q[1]]);
        g.define(AXI, young, nu, None, formulation, Procedure::Static {});
        g.symmetry("ends", Symmetry::Y, |q| near(q[1], 0.0) || near(q[1], 1.0));
        g.edge_set("bore", |q| near(q[0], 0.0));
        g.load("p", Load::Pressure { set: "bore".into(), magnitude: p });
        let r = g.solve();
        let mut eu: f64 = 0.0;
        for (label, u) in r.nodes.iter().zip(&r.displacement) {
            let x = g.model.node(*label).unwrap().xyz;
            eu = eu.max((u[0] - ur(x[0])).abs() / ur(a));
        }
        let bore = g.nodes_where(|q| near(q[0], 0.0));
        let st = bore.iter().map(|l| r.stress[at(&r, *l)][2]).sum::<f64>() / bore.len() as f64;
        let sz = r.stress.iter().zip(&r.nodes).map(|(s, _)| s[1]).sum::<f64>() / r.nodes.len() as f64;
        let es = (st - sigma_t(a)).abs() / sigma_t(a);
        eprintln!(
            "lamé axisymmetric {name}: u_r error {:.3}%, bore hoop {st:.4} vs {:.4} ({:.2}%), mean σ_z {sz:.4} vs {:.4}",
            100.0 * eu,
            sigma_t(a),
            100.0 * es,
            2.0 * nu * big_a
        );
        assert!(eu < 0.01 && es < 0.03 && (sz - 2.0 * nu * big_a).abs() < 0.02 * sigma_t(a), "{name}");
    }
}

/// Rotating discs and cylinders (S. Timoshenko, J. N. Goodier, *Theory of Elasticity*, 3rd ed.,
/// McGraw-Hill 1970, §32), solid, radius b, speed ω, density ρ:
/// - a thin disc (plane stress, rotating about Z): σ_r = (3 + ν)/8 ρω² (b² − r²),
///   σ_θ = ρω²/8 ((3 + ν) b² − (1 + 3ν) r²);
/// - a long cylinder (axisymmetric, axial strain held at zero):
///   σ_r = (3 − 2ν)/(8(1 − ν)) ρω² (b² − r²), σ_θ = ρω²/(8(1 − ν)) ((3 − 2ν) b² − (1 + 2ν) r²).
#[test]
fn rotating_solid_disc_and_cylinder() {
    let (b, rho, w, young, nu) = (1.0f64, 1e-3, 30.0, 1000.0, 0.3);
    let k = rho * w * w;
    let disc = |r: f64| [(3.0 + nu) / 8.0 * k * (b * b - r * r), k / 8.0 * ((3.0 + nu) * b * b - (1.0 + 3.0 * nu) * r * r)];
    let cylinder = |r: f64| {
        let c = k / (8.0 * (1.0 - nu));
        [c * (3.0 - 2.0 * nu) * (b * b - r * r), c * ((3.0 - 2.0 * nu) * b * b - (1.0 + 2.0 * nu) * r * r)]
    };
    let peak = disc(0.0)[0];
    // Linear elements resolve the quadratic stress field to about their size (h/b of the peak
    // in nodal averages); they get a finer mesh.
    for (kind, formulation, n, tol, name) in [
        (ElementKind::Quad8, full(), 8, 0.01, "quad8"),
        (ElementKind::Tri6, full(), 8, 0.01, "tri6"),
        (ElementKind::Quad4, incompatible(), 32, 0.04, "quad4 incompatible"),
    ] {
        // Plane stress quarter disc, polar-mapped from a pin hole of radius b/50 (a polar map
        // folds at the centre); the hole's disturbance decays as (r₀/r)², so the solid disc's
        // field is compared at r ≥ b/2.
        let mut g = grid([n, 8], kind, |q| {
            let (r, th) = (b * (0.02 + 0.98 * q[0]), PI / 2.0 * q[1]);
            [r * th.cos(), r * th.sin()]
        });
        g.define(SectionKind::PlaneStress { thickness: 0.05 }, young, nu, Some(rho), formulation, Procedure::Static {});
        g.symmetry("y0", Symmetry::Y, |q| near(q[1], 0.0));
        g.symmetry("x0", Symmetry::X, |q| near(q[1], 1.0));
        g.load("spin", Load::Centrifugal { angular_velocity: w, point: [0.0; 3], axis: [0.0, 0.0, 1.0], set: None });
        let r = g.solve();
        let mut e: f64 = 0.0;
        for (label, s) in r.nodes.iter().zip(&r.stress) {
            let x = g.model.node(*label).unwrap().xyz;
            let rad = x[0].hypot(x[1]);
            if rad < 0.5 * b {
                continue;
            }
            let (c, sn) = (x[0] / rad, x[1] / rad);
            let sr = s[0] * c * c + 2.0 * s[3] * sn * c + s[1] * sn * sn;
            let st = s[0] * sn * sn - 2.0 * s[3] * sn * c + s[1] * c * c;
            let want = disc(rad);
            e = e.max((sr - want[0]).abs().max((st - want[1]).abs()) / peak);
        }
        eprintln!("rotating disc (plane stress) {name}: stress error {:.3}% of the peak (r ≥ b/2)", 100.0 * e);
        assert!(e < tol, "{name}: {e}");

        // Axisymmetric solid cylinder 0 ≤ r ≤ b, a slice held axially on both faces: nodes on the
        // axis are held radially by the solver.
        let mut g = grid([n, 2], kind, |q| [b * q[0], 0.2 * q[1]]);
        g.define(AXI, young, nu, Some(rho), formulation, Procedure::Static {});
        g.symmetry("ends", Symmetry::Y, |q| near(q[1], 0.0) || near(q[1], 1.0));
        g.load("spin", Load::Centrifugal { angular_velocity: w, point: [0.0; 3], axis: [0.0, 1.0, 0.0], set: None });
        let r = g.solve();
        let cyl_peak = cylinder(0.0)[0];
        let mut e: f64 = 0.0;
        for (label, s) in r.nodes.iter().zip(&r.stress) {
            let x = g.model.node(*label).unwrap().xyz[0];
            let want = cylinder(x);
            e = e.max((s[0] - want[0]).abs().max((s[2] - want[1]).abs()) / cyl_peak);
        }
        let on_axis = g.nodes_where(|q| near(q[0], 0.0));
        assert!(on_axis.iter().all(|l| r.displacement[at(&r, *l)][0] == 0.0), "the axis does not move radially");
        eprintln!("rotating cylinder (axisymmetric) {name}: stress error {:.3}% of the peak", 100.0 * e);
        assert!(e < 2.0 * tol, "{name}: {e}");
    }
}

/// The total mass of a 2D model is its area times its thickness, or the volume of revolution
/// (Pappus) when axisymmetric; a free plate has three rigid-body modes in its plane.
#[test]
fn masses_and_rigid_modes() {
    let rho = 2.0;
    for (theory, mass) in [(STRESS, rho * 2.0 * 0.5), (AXI, rho * 2.0 * PI * 1.5 * 2.0)] {
        let mut g = grid([4, 2], ElementKind::Quad8, |p| [1.0 + p[0], 2.0 * p[1]]);
        g.define(theory.clone(), 1000.0, 0.3, Some(rho), full(), Procedure::Frequency { modes: 4, min_hz: None, max_hz: None });
        let f = solve_frequency(&g.model, None, &FrequencyOptions::default(), &Limits::default()).unwrap();
        eprintln!("{}: mass {} vs {mass}, frequencies {:?}", theory.name(), f.total_mass, f.modes.iter().map(|m| m.frequency).collect::<Vec<_>>());
        assert!((f.total_mass - mass).abs() < 1e-12 * mass);
        let rigid = f.modes.iter().filter(|m| m.frequency < 1e-3 * f.modes.last().unwrap().frequency).count();
        // In-plane: two translations and a rotation; axisymmetric: the axial translation only.
        assert_eq!(rigid, if theory == AXI { 1 } else { 3 }, "{}", theory.name());
    }
}

#[test]
fn two_dimensional_results_are_deterministic_and_errors_name_the_problem() {
    let make = || {
        let mut g = grid([4, 2], ElementKind::Quad4, |p| [1.0 + p[0] + 0.1 * p[1] * p[0], p[1]]);
        g.define(AXI, 1000.0, 0.3, Some(1.0), incompatible(), Procedure::Static {});
        g.symmetry("bottom", Symmetry::Y, |p| near(p[1], 0.0));
        g.edge_set("bore", |p| near(p[0], 0.0));
        g.load("p", Load::Pressure { set: "bore".into(), magnitude: 3.0 });
        g.load("g", Load::Gravity { acceleration: [0.0, -9.81, 0.0], set: None });
        g
    };
    let g = make();
    assert_eq!(g.solve(), g.solve());
    // A coupling is not supported in 2D yet.
    let mut g = make();
    g.model.add_reference_point("RP", fieldcraft_model::ReferencePoint { xyz: [0.0, 2.0, 0.0], label: 0 }).unwrap();
    g.model
        .add_coupling(
            "c",
            fieldcraft_model::Coupling {
                reference_point: "RP".into(),
                set: "bottom".into(),
                kind: fieldcraft_model::CouplingKind::Kinematic,
                dofs: fieldcraft_model::Dofs::all(),
            },
        )
        .unwrap();
    let e = solve_static(&g.model, None, &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("couplings in 2D models"), "{e}");
    // A load with a Z component, and a rotation about a wrong axis, are refused by the check.
    let mut g = make();
    g.load("lift", Load::Gravity { acceleration: [0.0, 0.0, 1.0], set: None });
    let e = solve_static(&g.model, None, &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("Z component"), "{e}");
    // No support: singular, naming a node and direction.
    let mut g = grid([2, 1], ElementKind::Quad4, |p| [p[0], p[1]]);
    g.define(STRESS, 1000.0, 0.3, None, full(), Procedure::Static {});
    g.edge_set("right", |p| near(p[0], 1.0));
    g.load("p", Load::Pressure { set: "right".into(), magnitude: 1.0 });
    let e = solve_static(&g.model, None, &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("singular") && e.contains("node"), "{e}");
}

/// Give the material an expansion coefficient and heat every node by `dt` above the reference.
fn heat(g: &mut Grid, alpha: f64, dt: f64) {
    let mut m = g.model.material("Mat").unwrap().clone();
    m.thermal = Some(fieldcraft_model::Thermal {
        expansion: Some(fieldcraft_model::Directional::Isotropic(alpha)),
        reference_temperature: 20.0,
        ..Default::default()
    });
    g.model.replace_material("Mat", m).unwrap();
    let nodes: Vec<u32> = g.param.keys().copied().collect();
    g.model.add_set("hot", Set { nodes, ..Set::default() }).unwrap();
    g.load("heat", Load::Temperature { target: Target::Set("hot".into()), value: 20.0 + dt });
}

#[test]
fn thermal_expansion_in_every_theory() {
    // Uniform ΔT, held only by symmetry: plane stress and axisymmetric expand freely
    // (u = αΔT x, no stress; εzz or the hoop strain αΔT), plane strain holds εzz = 0, so
    // σzz = −EαΔT and the in-plane strain is (1 + ν)αΔT.
    let (e, nu, alpha, dt) = (2e5, 0.3, 1.2e-5, 80.0);
    for (kind, f, name) in variants() {
        for section in [SectionKind::PlaneStress { thickness: 0.5 }, SectionKind::PlaneStrain { thickness: 0.5 }, SectionKind::Axisymmetric {}] {
            let axisymmetric = matches!(section, SectionKind::Axisymmetric {});
            let mut g = grid([2, 2], kind, |p| [0.5 + 2.0 * p[0] + 0.2 * p[1], p[1]]);
            g.define(section.clone(), e, nu, None, f, Procedure::Static {});
            heat(&mut g, alpha, dt);
            g.symmetry("bottom", Symmetry::Y, |p| near(p[1], 0.0) && near(p[0], 0.0));
            if !axisymmetric {
                g.symmetry("left", Symmetry::X, |p| near(p[0], 0.0) && near(p[1], 0.0));
                let corner = g.nodes_where(|p| near(p[0], 1.0) && near(p[1], 0.0));
                g.model.add_set("roller", Set { nodes: corner, ..Set::default() }).unwrap();
                g.bc("roller", BoundaryCondition::Symmetry { target: Target::Set("roller".into()), plane: Symmetry::Y });
            }
            let r = g.solve();
            let plane_strain = matches!(section, SectionKind::PlaneStrain { .. });
            let strain = if plane_strain { (1.0 + nu) * alpha * dt } else { alpha * dt };
            for (l, i) in r.nodes.iter().zip(0..) {
                let x = g.model.node(*l).unwrap().xyz;
                let x0 = g.model.node(g.nodes_where(|p| near(p[0], 0.0) && near(p[1], 0.0))[0]).unwrap().xyz;
                // Axisymmetric: the radial displacement is αΔT r (no rigid radial motion).
                let ux = if axisymmetric { strain * x[0] } else { strain * (x[0] - x0[0]) };
                assert!((r.displacement[i][0] - ux).abs() < 1e-9, "{name} {section:?} node {l}: {:?}", r.displacement[i]);
                let uy = if axisymmetric { alpha * dt * x[1] } else { strain * x[1] };
                assert!((r.displacement[i][1] - uy).abs() < 1e-9, "{name} {section:?} node {l}: {:?}", r.displacement[i]);
            }
            for s in &r.stress {
                let zz = if plane_strain { -e * alpha * dt } else { 0.0 };
                assert!((s[0].abs() + s[1].abs() + s[3].abs()) < 1e-6 && (s[2] - zz).abs() < 1e-6, "{name} {section:?}: {s:?}");
            }
        }
    }
}

#[test]
fn a_plate_held_on_every_edge_carries_the_biaxial_thermal_stress() {
    // Plane stress, normal displacements held on all four edges: σxx = σyy = −EαΔT / (1 − ν).
    let (e, nu, alpha, dt) = (2e5, 0.3, 1.2e-5, 50.0);
    let want = -e * alpha * dt / (1.0 - nu);
    for (kind, f, name) in variants() {
        let mut g = grid([2, 2], kind, |p| [p[0], p[1]]);
        g.define(SectionKind::PlaneStress { thickness: 1.0 }, e, nu, None, f, Procedure::Static {});
        heat(&mut g, alpha, dt);
        g.symmetry("x0", Symmetry::X, |p| near(p[0], 0.0));
        g.symmetry("x1", Symmetry::X, |p| near(p[0], 1.0));
        g.symmetry("y0", Symmetry::Y, |p| near(p[1], 0.0));
        g.symmetry("y1", Symmetry::Y, |p| near(p[1], 1.0));
        let r = g.solve();
        for s in &r.stress {
            assert!((s[0] - want).abs() < 1e-8 * want.abs() && (s[1] - want).abs() < 1e-8 * want.abs() && s[3].abs() < 1e-6, "{name}: {s:?}");
        }
    }
}
