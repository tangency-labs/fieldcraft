//! Hyperelastic materials in explicit dynamics: the stress of the shared finite-strain material
//! point (neo-Hookean, Mooney–Rivlin, Yeoh) from the deformation gradient, on one-point hex8 and
//! F-bar tet4, against the closed form of uniaxial strain reached quasi-statically, and the
//! refusals.
//!
//! In uniaxial strain F = diag(λ, 1, 1): the hyperelastic Cauchy stress is the material point's
//! Kirchhoff stress over J = λ (fieldcraft_fem::constitutive, whose laws are checked there
//! against their strain energies), so the explicit elements must reproduce it once the motion has
//! settled; the wave dynamics of a slow pull leave a fraction of a percent.

use fieldcraft_explicit::{ExplicitResult, Limits, solve_explicit};
use fieldcraft_fem::constitutive::MaterialModel;
use fieldcraft_model::*;

const N: usize = 2;

/// A unit cube of N³ hex8 (or six tet4 per cell), held on its faces x = 0, y = 0, y = 1, z = 0,
/// z = 1 in their normal directions and pulled on x = 1 to a stretch λ over `end`.
fn block(tet: bool, hyper: Hyperelastic, stretch: f64, end: f64, nonlinear: bool) -> Model {
    let mut m = Model::default();
    let h = 1.0 / N as f64;
    let label = |i: usize, j: usize, k: usize| (1 + i + (N + 1) * (j + (N + 1) * k)) as u32;
    for k in 0..=N {
        for j in 0..=N {
            for i in 0..=N {
                m.set_node(label(i, j, k), [i as f64 * h, j as f64 * h, k as f64 * h]).unwrap();
            }
        }
    }
    let corners = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
    let mut next = 1;
    for k in 0..N {
        for j in 0..N {
            for i in 0..N {
                let at = |c: [usize; 3]| label(i + c[0], j + c[1], k + c[2]);
                if !tet {
                    m.set_element(next, ElementKind::Hex8, corners.iter().map(|c| at(*c)).collect()).unwrap();
                    next += 1;
                    continue;
                }
                for perm in [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]] {
                    let mut v = [[0usize; 3]; 4];
                    for s in 0..3 {
                        v[s + 1] = v[s];
                        v[s + 1][perm[s]] = 1;
                    }
                    let e = |a: usize| [0, 1, 2].map(|d| v[a][d] as f64 - v[0][d] as f64);
                    let (a, b, c) = (e(1), e(2), e(3));
                    let vol = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
                    if vol < 0.0 {
                        v.swap(1, 2);
                    }
                    m.set_element(next, ElementKind::Tet4, v.iter().map(|c| at(*c)).collect()).unwrap();
                    next += 1;
                }
            }
        }
    }
    let all: Vec<u32> = m.elements().keys().copied().collect();
    m.add_set("all", Set { elements: all, ..Set::default() }).unwrap();
    let nodes = |f: &dyn Fn([f64; 3]) -> bool| -> Vec<u32> { m.nodes().iter().filter(|(_, n)| f(n.xyz)).map(|(l, _)| *l).collect() };
    let sets = [
        ("x0", nodes(&|x| x[0] == 0.0)),
        ("x1", nodes(&|x| x[0] == 1.0)),
        ("y", nodes(&|x| x[1] == 0.0 || x[1] == 1.0)),
        ("z", nodes(&|x| x[2] == 0.0 || x[2] == 1.0)),
    ];
    for (name, list) in sets {
        m.add_set(name, Set { nodes: list, ..Set::default() }).unwrap();
    }
    let elastic = Elastic::Isotropic { young: 1.0, poisson: 0.45 };
    m.add_material("Rubber", Material { elastic: Some(elastic), density: Some(1.0), hyperelastic: Some(hyper), ..Material::default() }).unwrap();
    let integration = if tet { Integration::Full } else { Integration::Reduced };
    let formulation = Formulation { integration, ..Formulation::default() };
    m.add_section("all", Section { orientation: None, set: "all".into(), material: "Rubber".into(), kind: SectionKind::Solid {}, formulation })
        .unwrap();
    let held = |set: &str, k: usize, v: f64| {
        let mut values = [None; 6];
        values[k] = Some(v);
        BoundaryCondition::Displacement { target: Target::Set(set.into()), values }
    };
    let mut step = ExplicitDynamic::new(end);
    step.nonlinear_geometry = nonlinear;
    step.output_interval = Some(end / 2.0);
    let mut s = Step::new("Pull", Procedure::ExplicitDynamic(step));
    for (name, bc) in [("x0", held("x0", 0, 0.0)), ("x1", held("x1", 0, stretch - 1.0)), ("y", held("y", 1, 0.0)), ("z", held("z", 2, 0.0))] {
        m.add_boundary_condition(name, bc).unwrap();
        s.boundary_conditions.insert(name.into());
    }
    m.add_step(s, None).unwrap();
    m
}

/// σ of uniaxial strain λ from the material point.
fn closed_form(hyper: Hyperelastic, stretch: f64) -> [f64; 6] {
    let model = MaterialModel::from_model(&Elastic::Isotropic { young: 1.0, poisson: 0.45 }, Some(&hyper)).unwrap();
    let f = [[stretch, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    model.response(&f, None, &[], &mut []).unwrap().stress.map(|t| t / stretch)
}

fn mean_stress(r: &ExplicitResult) -> [f64; 6] {
    let f = r.frames.last().unwrap();
    let mut s = [0.0; 6];
    for e in &f.stress {
        for k in 0..6 {
            s[k] += e[k] / f.stress.len() as f64;
        }
    }
    s
}

#[test]
fn hyperelastic_uniaxial_strain_matches_the_material_point() {
    let stretch = 1.5;
    // Dilatational wave speed √((κ + 4μ/3)/ρ) ≈ 2.0: a pull over 400 transits of the cube.
    let end = 200.0;
    let materials = [Hyperelastic::NeoHookean {}, Hyperelastic::MooneyRivlin { c01: 0.08 }, Hyperelastic::Yeoh { c20: -0.01, c30: 0.002 }];
    for hyper in materials {
        let reference = closed_form(hyper, stretch);
        for tet in [false, true] {
            let r = solve_explicit(&block(tet, hyper, stretch, end, true), None, &Limits::default()).unwrap();
            assert!(r.completed);
            let s = mean_stress(&r);
            let spread = r.frames.last().unwrap().stress.iter().map(|e| (e[0] - s[0]).abs()).fold(0.0, f64::max);
            eprintln!(
                "{} {}: σ_xx {:.5} vs {:.5}, σ_yy {:.5} vs {:.5}; spread {spread:.2e}; balance {:.2e}; {} increments",
                hyper.name(),
                if tet { "tet4" } else { "hex8" },
                s[0],
                reference[0],
                s[1],
                reference[1],
                r.energy.relative_error(r.initial_kinetic),
                r.increments
            );
            for k in 0..3 {
                assert!((s[k] - reference[k]).abs() < 0.01 * reference[0].abs(), "{} tet {tet}: {s:?} vs {reference:?}", hyper.name());
            }
            assert!(spread < 0.01 * reference[0].abs(), "uniform state: {spread}");
            assert!(r.energy.relative_error(r.initial_kinetic) < 0.01, "{:?}", r.energy);
        }
    }
}

#[test]
fn hyperelastic_materials_need_nonlinear_geometry() {
    let e = solve_explicit(&block(false, Hyperelastic::NeoHookean {}, 1.1, 1.0, false), None, &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("nonlinear geometry") && e.contains("Rubber"), "{e}");
}
