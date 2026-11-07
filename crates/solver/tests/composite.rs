//! Independent laminate solver checks. Constitutive and plate conventions are documented in
//! docs/solvers/laminates.md. Navier's separated sine/cosine amplitudes give the small 3 × 3
//! FSDT eigenproblem; the strip is the exact cylindrical pure-bending solution, not a fixture
//! from another solver. All dimensions and material values here are SI.

use std::f64::consts::PI;

mod common;

use fieldcraft_fem::laminate::{Laminate, resolve_plies};
use fieldcraft_model::composite::MAX_PLIES;
use fieldcraft_model::*;
use fieldcraft_solver::{FrequencyOptions, Limits, StaticResult, solve_buckle, solve_frequency, solve_static};

struct Plate {
    model: Model,
    nx: usize,
    ny: usize,
}

impl Plate {
    fn new(nx: usize, ny: usize, length: [f64; 2], thickness: f64, angles: &[f64], procedure: Procedure) -> Self {
        let mut model = Model::default();
        model
            .add_material(
                "ply",
                Material {
                    elastic: Some(Elastic::EngineeringConstants { young: [100e9, 10e9, 10e9], poisson: [0.0; 3], shear: [5e9, 3.5e9, 3.5e9] }),
                    density: Some(1600.0),
                    ply_strength: Some(PlyStrength {
                        tension: [500e6, 5e6],
                        compression: [350e6, 75e6],
                        shear12: 35e6,
                        shear23: 35e6,
                        tsai_wu_interaction: -0.5,
                    }),
                    ..Material::default()
                },
            )
            .unwrap();
        let mut p = Self { model, nx, ny };
        for j in 0..=ny {
            for i in 0..=nx {
                p.model.set_node(p.node(i, j), [length[0] * i as f64 / nx as f64, length[1] * j as f64 / ny as f64, 0.0]).unwrap();
            }
        }
        let mut elements = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                let label = (1 + i + nx * j) as u32;
                p.model.set_element(label, ElementKind::Quad4, vec![p.node(i, j), p.node(i + 1, j), p.node(i + 1, j + 1), p.node(i, j + 1)]).unwrap();
                elements.push(label);
            }
        }
        p.model.add_set("plate", Set { elements, ..Set::default() }).unwrap();
        let layup = Layup {
            offset: 0.0,
            plies: angles.iter().map(|a| Ply { material: "ply".into(), thickness: thickness / angles.len() as f64, angle: *a, points: 3 }).collect(),
        };
        p.model
            .add_section(
                "plate",
                Section {
                    set: "plate".into(),
                    material: String::new(),
                    kind: SectionKind::CompositeShell { layup },
                    orientation: Some(SectionOrientation {
                        frame: MaterialOrientation::Rectangular { axis1: [1.0, 0.0, 0.0], axis2: [0.0, 1.0, 0.0] },
                        shell_reference_axis: 1,
                    }),
                    formulation: Formulation::default(),
                },
            )
            .unwrap();
        p.model.add_step(Step::new("step", procedure), None).unwrap();
        p
    }

    fn node(&self, i: usize, j: usize) -> u32 {
        (1 + i + (self.nx + 1) * j) as u32
    }

    fn hold(&mut self, name: &str, values: [Option<f64>; 6], choose: impl Fn(usize, usize) -> bool) {
        let nodes =
            (0..=self.ny).flat_map(|j| (0..=self.nx).map(move |i| (i, j))).filter(|(i, j)| choose(*i, *j)).map(|(i, j)| self.node(i, j)).collect();
        self.model.add_set(name, Set { nodes, ..Set::default() }).unwrap();
        self.model.add_boundary_condition(name, BoundaryCondition::Displacement { target: Target::Set(name.into()), values }).unwrap();
        self.model.set_boundary_condition_active("step", name, true).unwrap();
    }

    fn load(&mut self, name: &str, load: Load) {
        self.model.add_load(name, load).unwrap();
        self.model.set_load_active("step", name, true).unwrap();
    }

    fn properties(&self) -> Laminate {
        let SectionKind::CompositeShell { layup } = &self.model.section("plate").unwrap().kind else { panic!("test layup") };
        Laminate::new(&resolve_plies(layup, &self.model).unwrap(), layup.offset).unwrap()
    }

    fn target_node(&mut self, label: u32) -> Target {
        let name = format!("node{label}");
        if self.model.set(&name).is_none() {
            self.model.add_set(&name, Set { nodes: vec![label], ..Set::default() }).unwrap();
        }
        Target::Set(name)
    }

    fn supported(&mut self, membrane_fixed: bool) {
        let (nx, ny) = (self.nx, self.ny);
        // Hard simply supported Navier edges: w=0; rotation about the edge normal is zero.
        self.hold("edges", [None, None, Some(0.0), None, None, None], |i, j| i == 0 || i == nx || j == 0 || j == ny);
        self.hold("x_edges", [None, None, None, Some(0.0), None, None], |i, _| i == 0 || i == nx);
        self.hold("y_edges", [None, None, None, None, Some(0.0), None], |_, j| j == 0 || j == ny);
        if membrane_fixed {
            self.hold("membrane", [Some(0.0), Some(0.0), None, None, None, Some(0.0)], |_, _| true);
        } else {
            self.hold("drill", [None, None, None, None, None, Some(0.0)], |_, _| true);
            self.hold("left", [Some(0.0), None, None, None, None, None], |i, _| i == 0);
            self.hold("origin", [None, Some(0.0), None, None, None, None], |i, j| i == 0 && j == 0);
        }
    }
}

fn at(result: &StaticResult, label: u32) -> usize {
    result.nodes.binary_search(&label).unwrap()
}

#[test]
fn laminate_strip_pure_bending_and_layer_fields() {
    let mut p = Plate::new(8, 2, [1.0, 0.1], 0.02, &[0.0, PI / 2.0, PI / 2.0, 0.0], Procedure::Static {});
    p.hold("clamp", [Some(0.0); 6], |i, _| i == 0);
    let end: Vec<u32> = (0..=p.ny).map(|j| p.node(p.nx, j)).collect();
    for (j, label) in end.iter().enumerate() {
        let moment = 0.1 / p.ny as f64 * if j == 0 || j == p.ny { 0.5 } else { 1.0 };
        let target = p.target_node(*label);
        p.load(&format!("moment{j}"), Load::Moment { target, vector: [0.0, moment, 0.0], distribution: Distribution::PerNode });
    }
    let r = solve_static(&p.model, None, &Limits::default()).unwrap();
    let kappa = 1.0 / p.properties().abd()[3][3];
    for (label, u) in r.nodes.iter().zip(&r.displacement) {
        let x = p.model.node(*label).unwrap().xyz[0];
        assert!((u[2] + 0.5 * kappa * x * x).abs() < 2e-9 * kappa, "{label}: {u:?}, kappa {kappa}");
    }
    assert_eq!(r.shells.len(), 16);
    for e in &r.shells {
        assert_eq!(e.stress.plies.len(), 4);
        for layers in &e.stress.plies {
            assert_eq!(layers.len(), 4);
            for (k, layer) in layers.iter().enumerate() {
                let s = layer.middle.stress;
                let axial = if k == 0 || k == 3 { s[0] } else { s[1] };
                let modulus = if k == 0 || k == 3 { 100e9 } else { 10e9 };
                assert!((axial - modulus * kappa * layer.middle.z).abs() < 2e-8 * 100e9 * kappa * 0.02);
                assert!(layer.middle.failure.is_some());
                assert_eq!(layer.points.len(), 3);
            }
        }
    }
    eprintln!("composite strip tip={:.12e}, exact={:.12e}", r.displacement[at(&r, end[1])][2], -kappa / 2.0);
}

#[test]
fn ply_result_storage_is_bounded_before_recovery() {
    let mut plate = Plate::new(1, 1, [1.0, 0.1], 0.01, &[0.0], Procedure::Static {});
    plate.hold("clamp", [Some(0.0); 6], |i, _| i == 0);
    let target = plate.target_node(plate.node(1, 0));
    plate.load("tip", Load::Force { target, vector: [0.0, 0.0, -1.0], distribution: Distribution::PerNode });
    let limits = Limits { max_factor_entries: 1000, ..Limits::default() };
    assert!(solve_static(&plate.model, None, &limits).is_ok(), "one-ply fields and the same small sparse factor fit");
    let mut section = plate.model.section("plate").unwrap().clone();
    if let SectionKind::CompositeShell { layup } = &mut section.kind {
        layup.plies = vec![Ply { material: "ply".into(), thickness: 0.01 / MAX_PLIES as f64, angle: 0.0, points: 15 }; MAX_PLIES];
    }
    plate.model.replace_section("plate", section).unwrap();
    let error = solve_static(&plate.model, None, &limits).unwrap_err();
    assert!(matches!(error, fieldcraft_solver::Error::Limit(_)));
    assert!(error.to_string().contains("laminate ply results"));
}

/// First flexural eigenvalue of the Navier FSDT amplitude system with rotary inertia.
/// Krr=(D11 a²+D66 b²+Sx, (D12+D66)ab; ..., D66 a²+D22 b²+Sy),
/// Krw=(Sx a,Sy b), Kww=Sx a²+Sy b²; M=diag(I2,I2,I0).
fn navier_lambda(p: &Laminate, m: usize, n: usize) -> f64 {
    let abd = p.abd();
    let s = p.shear();
    let mass = p.mass_moments().unwrap();
    let a = m as f64 * PI;
    let b = n as f64 * PI;
    let r1 = abd[3][3] * a * a + abd[5][5] * b * b + s[0][0];
    let r2 = abd[5][5] * a * a + abd[4][4] * b * b + s[1][1];
    let r12 = (abd[3][4] + abd[5][5]) * a * b;
    let g1 = s[0][0] * a;
    let g2 = s[1][1] * b;
    let mut low = 0.0;
    let mut high = (abd[3][3] * a.powi(4) + 2.0 * (abd[3][4] + 2.0 * abd[5][5]) * a * a * b * b + abd[4][4] * b.powi(4)) / mass[0];
    for _ in 0..100 {
        let v = 0.5 * (low + high);
        let x = r1 - v * mass[2];
        let y = r2 - v * mass[2];
        let f = s[0][0] * a * a + s[1][1] * b * b - v * mass[0] - (g1 * g1 * y + g2 * g2 * x - 2.0 * g1 * g2 * r12) / (x * y - r12 * r12);
        if f > 0.0 {
            low = v;
        } else {
            high = v;
        }
    }
    0.5 * (low + high)
}

#[test]
fn laminate_plate_navier_frequencies_thin_and_thick() {
    for h in [0.01, 0.1] {
        let mut errors = Vec::new();
        for cells in [4, 8, 16] {
            let mut p =
                Plate::new(cells, cells, [1.0; 2], h, &[0.0, PI / 2.0, PI / 2.0, 0.0], Procedure::Frequency { modes: 3, min_hz: None, max_hz: None });
            p.supported(true);
            let properties = p.properties();
            let mut reference: Vec<f64> =
                (1..=3).flat_map(|m| (1..=3).map(move |n| (m, n))).map(|(m, n)| navier_lambda(&properties, m, n).sqrt() / (2.0 * PI)).collect();
            reference.sort_by(f64::total_cmp);
            let r = solve_frequency(&p.model, None, &FrequencyOptions::default(), &Limits::default()).unwrap();
            assert_eq!(r.modes.len(), 3);
            let err = (r.modes[0].frequency / reference[0] - 1.0).abs();
            errors.push(err);
            eprintln!(
                "composite frequency h={h}, n={cells}, got={:?}, Navier={:?}",
                r.modes.iter().map(|m| m.frequency).collect::<Vec<_>>(),
                &reference[..3]
            );
            if cells == 16 {
                for threads in [1, 2, 7] {
                    let repeated = fieldcraft_parallel::with_threads(threads, || {
                        solve_frequency(&p.model, None, &FrequencyOptions::default(), &Limits::default()).unwrap()
                    });
                    assert_eq!(repeated, r, "modal threads={threads}");
                }
                for (mode, want) in r.modes.iter().zip(&reference) {
                    assert!((mode.frequency / want - 1.0).abs() < 0.05);
                }
                assert!(err < 0.02);
            }
        }
        assert!(errors[1] < 0.4 * errors[0]);
        assert!(errors[2] < 0.4 * errors[1]);
    }
}

#[test]
fn laminate_flat_plate_buckling_and_thread_determinism() {
    let mut p = Plate::new(8, 8, [1.0; 2], 0.01, &[0.0, PI / 2.0, PI / 2.0, 0.0], Procedure::Buckle { modes: 3, shift: None });
    p.supported(false);
    for j in 0..=p.ny {
        let f = -(if j == 0 || j == p.ny { 0.5 } else { 1.0 }) / p.ny as f64;
        let target = p.target_node(p.node(p.nx, j));
        p.load(&format!("compression{j}"), Load::Force { target, vector: [f, 0.0, 0.0], distribution: Distribution::PerNode });
    }
    let mut step = p.model.step("step").unwrap().clone();
    step.name = "preload".into();
    step.procedure = Procedure::Static {};
    p.model.add_step(step, None).unwrap();
    let d = p.properties().abd();
    let want = PI * PI * (d[3][3] + 2.0 * (d[3][4] + 2.0 * d[5][5]) + d[4][4]);
    let run = |threads| {
        fieldcraft_parallel::with_threads(threads, || {
            let s = solve_static(&p.model, Some("preload"), &Limits::default()).unwrap();
            let b = solve_buckle(&p.model, Some("step"), &Limits::default()).unwrap();
            (s, b)
        })
    };
    let first = run(1);
    let got = first.1.modes[0].eigenvalue;
    eprintln!("composite buckling first={got:.12e}, thin CLT={want:.12e}");
    assert!((got / want - 1.0).abs() < 0.04);
    for threads in [2, 7] {
        assert_eq!(run(threads), first, "threads={threads}");
    }
    // An unsymmetric plate bends under membrane loading; the flat-plate operator must not
    // silently ignore that pre-buckling deformation.
    let mut section = p.model.section("plate").unwrap().clone();
    if let SectionKind::CompositeShell { layup } = &mut section.kind {
        layup.plies = vec![layup.plies[0].clone(), layup.plies[1].clone()];
    }
    p.model.replace_section("plate", section).unwrap();
    let e = solve_buckle(&p.model, Some("step"), &Limits::default()).unwrap_err();
    assert!(e.to_string().contains("non-flat shell preload"), "{e}");
}

#[test]
fn laminate_missing_density_is_not_filled_for_modes_or_gravity() {
    let mut p = Plate::new(2, 2, [1.0; 2], 0.1, &[0.0, PI / 2.0], Procedure::Static {});
    p.hold("clamp", [Some(0.0); 6], |i, _| i == 0);
    let mut m = p.model.material("ply").unwrap().clone();
    m.density = None;
    p.model.replace_material("ply", m).unwrap();
    p.load("gravity", Load::Gravity { acceleration: [0.0, 0.0, -1.0], set: None });
    let e = solve_static(&p.model, None, &Limits::default()).unwrap_err();
    assert!(e.to_string().contains("no density"), "{e}");
    let mut step = p.model.step("step").unwrap().clone();
    step.procedure = Procedure::Frequency { modes: 1, min_hz: None, max_hz: None };
    p.model.replace_step(step).unwrap();
    let e = solve_frequency(&p.model, None, &FrequencyOptions::default(), &Limits::default()).unwrap_err();
    assert!(e.to_string().contains("no density"), "{e}");
}

#[test]
fn ply_density_first_moment_reaches_gravity_reactions() {
    let mut p = Plate::new(3, 1, [1.0, 0.1], 0.02, &[0.0, PI / 2.0], Procedure::Static {});
    p.hold("clamp", [Some(0.0); 6], |i, _| i == 0);
    let mut material = p.model.material("ply").unwrap().clone();
    material.density = Some(3200.0);
    p.model.add_material("heavy", material).unwrap();
    let mut section = p.model.section("plate").unwrap().clone();
    if let SectionKind::CompositeShell { layup } = &mut section.kind {
        layup.plies[1].material = "heavy".into();
        layup.offset = 0.2;
    }
    p.model.replace_section("plate", section).unwrap();
    let mass = p.properties().mass_moments().unwrap();
    p.load("gravity", Load::Gravity { acceleration: [1.0, 0.0, 0.0], set: None });
    let r = solve_static(&p.model, None, &Limits::default()).unwrap();
    let fx: f64 = r.reaction.iter().map(|r| r[0]).sum();
    let my: f64 = r.reaction_moment.iter().map(|r| r[1]).sum();
    assert!((fx + 0.1 * mass[0]).abs() < 1e-10 * mass[0]);
    assert!((my + 0.1 * mass[1]).abs() < 1e-9 * mass[1].abs(), "moment {my}, I1 {}", mass[1]);
}

#[test]
fn solid_material_direction_reaches_assembly_and_recovery() {
    let mut g = common::grid([2, 1, 1], ElementKind::Hex8, |x| x);
    g.define(
        Elastic::EngineeringConstants { young: [100e9, 10e9, 20e9], poisson: [0.0; 3], shear: [5e9, 4e9, 3e9] },
        Some(1600.0),
        Formulation::default(),
        Procedure::Static {},
    );
    let mut section = g.model.section("solid").unwrap().clone();
    section.orientation = Some(SectionOrientation {
        frame: MaterialOrientation::Rectangular { axis1: [0.0, 1.0, 0.0], axis2: [-1.0, 0.0, 0.0] },
        shell_reference_axis: 1,
    });
    g.model.replace_section("solid", section).unwrap();
    g.fixed("left", |x| x[0] == 0.0);
    g.node_set("right", |x| x[0] == 1.0);
    g.load("pull", Load::Force { target: Target::Set("right".into()), vector: [1.0, 0.0, 0.0], distribution: Distribution::Total });
    let r = solve_static(&g.model, None, &Limits::default()).unwrap();
    for (label, u) in r.nodes.iter().zip(&r.displacement) {
        let x = g.model.node(*label).unwrap().xyz[0];
        assert!((u[0] - x / 10e9).abs() < 1e-22);
    }
    for e in &r.elements {
        for p in &e.points {
            assert!((p.stress[0] - 1.0).abs() < 1e-12);
        }
    }
    let mut step = g.model.step("Step").unwrap().clone();
    step.procedure = Procedure::NonlinearStatic(NonlinearStatic::default());
    g.model.replace_step(step).unwrap();
    let e = fieldcraft_solver::solve_nonlinear(&g.model, None, &Limits::default()).unwrap_err();
    assert!(e.to_string().contains("material orientation"), "{e}");
}

#[test]
fn harmonic_oriented_solid_keeps_material_damping_and_stress_axes() {
    use fieldcraft_model::harmonic::{Damping, Frequencies, Harmonic, Method};
    use fieldcraft_solver::harmonic::solve_harmonic;
    let make = |material_damping: bool| {
        let mut g = common::grid([2, 1, 1], ElementKind::Hex8, |x| x);
        g.define(
            Elastic::EngineeringConstants { young: [100e9, 10e9, 20e9], poisson: [0.0; 3], shear: [5e9, 4e9, 3e9] },
            Some(1600.0),
            Formulation::default(),
            Procedure::Static {},
        );
        let mut section = g.model.section("solid").unwrap().clone();
        section.orientation = Some(SectionOrientation {
            frame: MaterialOrientation::Rectangular { axis1: [0.0, 1.0, 0.0], axis2: [-1.0, 0.0, 0.0] },
            shell_reference_axis: 1,
        });
        let damping = Damping { alpha: 0.2, beta: 0.0001, loss_factor: 0.03 };
        if material_damping {
            let mut material = g.model.materials()[&section.material].clone();
            material.damping = Some(damping);
            g.model.replace_material(&section.material, material).unwrap();
        }
        g.model.replace_section("solid", section).unwrap();
        g.fixed("left", |x| x[0] == 0.0);
        g.node_set("right", |x| x[0] == 1.0);
        g.load("pull", Load::Force { target: Target::Set("right".into()), vector: [1.0, 0.0, 0.0], distribution: Distribution::Total });
        let mut step = g.model.step("Step").unwrap().clone();
        step.procedure = Procedure::Harmonic(Harmonic {
            frequencies: Frequencies::List(vec![0.0, 20.0]),
            method: Method::Direct,
            damping: if material_damping { Damping::default() } else { damping },
            load_factors: Default::default(),
        });
        step.output.insert(OutputVariable::Stress);
        g.model.replace_step(step).unwrap();
        solve_harmonic(&g.model, None, &Limits::default()).unwrap()
    };
    let a = make(true);
    let b = make(false);
    for (a, b) in a.frames.iter().zip(&b.frames) {
        for (a, b) in a.real.iter().flatten().zip(b.real.iter().flatten()) {
            assert!((a - b).abs() < 1e-20, "{a} != {b}");
        }
        for (a, b) in a.imaginary.iter().flatten().zip(b.imaginary.iter().flatten()) {
            assert!((a - b).abs() < 1e-20, "{a} != {b}");
        }
    }
    // At zero frequency with a uniform loss factor, (1+i eta) K u = f.
    // Unit traction gives sigma_xx = 1/(1+i eta), independently of the rotated modulus.
    let zero = &a.frames[0];
    for s in zero.stress_real.as_ref().unwrap() {
        assert!((s[0] - 1.0 / (1.0 + 0.03_f64.powi(2))).abs() < 1e-10);
    }
}

#[test]
fn harmonic_laminate_static_limit_and_named_ply_damping_refusal() {
    use fieldcraft_model::harmonic::{Damping, Frequencies, Harmonic, Method};
    use fieldcraft_solver::harmonic::solve_harmonic;
    let mut p = Plate::new(2, 1, [1.0, 0.2], 0.02, &[0.0, 90.0], Procedure::Static {});
    p.hold("root", [Some(0.0); 6], |i, _| i == 0);
    let target = p.target_node(p.node(2, 1));
    p.load("tip", Load::Force { target, vector: [0.0, 0.0, -1.0], distribution: Distribution::Total });
    let static_result = solve_static(&p.model, None, &Limits::default()).unwrap();
    let mut step = p.model.step("step").unwrap().clone();
    step.procedure = Procedure::Harmonic(Harmonic {
        frequencies: Frequencies::List(vec![0.0]),
        method: Method::Direct,
        damping: Damping::default(),
        load_factors: Default::default(),
    });
    step.output.clear();
    step.output.insert(OutputVariable::Displacement);
    p.model.replace_step(step).unwrap();
    let h = solve_harmonic(&p.model, None, &Limits::default()).unwrap();
    for (u, v) in h.frames[0].real.iter().flatten().zip(static_result.displacement.iter().flatten()) {
        assert!((u - v).abs() < 1e-12);
    }
    let mut material = p.model.materials()["ply"].clone();
    material.damping = Some(Damping { alpha: 0.2, beta: 0.0, loss_factor: 0.0 });
    p.model.replace_material("ply", material).unwrap();
    let error = solve_harmonic(&p.model, None, &Limits::default()).unwrap_err().to_string();
    assert!(error.contains("per-ply material damping"), "{error}");
}
