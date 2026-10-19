//! Thermal strain in explicit dynamics: a temperature load reached slowly over the step expands
//! a free bar by α ΔT L without stress; a bar held at both ends carries σ = −E α ΔT, and with
//! perfectly plastic J2 the stress stops at the yield stress.

use fieldcraft_explicit::{Limits, solve_explicit};
use fieldcraft_model::*;

const YOUNG: f64 = 1.0e4;
const ALPHA: f64 = 1.0e-4;
const DT: f64 = 20.0;
const LENGTH: f64 = 4.0;

/// A 4 × 1 × 1 bar of hex8 (one-point) or Kuhn-split tet4 cells, ν = 0, heated by DT over a
/// step much longer than the wave transit, held at x = 0 (and at x = L when `clamped`) along x
/// and on the planes y = 0, z = 0.
fn bar(kind: ElementKind, clamped: bool, plastic: Option<Plastic>) -> Model {
    let mut m = Model::default();
    let n = 4usize;
    let label = |i: usize, j: usize, k: usize| (1 + i + (n + 1) * (j + 2 * k)) as u32;
    for k in 0..2 {
        for j in 0..2 {
            for i in 0..=n {
                m.set_node(label(i, j, k), [i as f64, j as f64, k as f64]).unwrap();
            }
        }
    }
    let corners = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
    let mut next = 1;
    for i in 0..n {
        let at = |c: [usize; 3]| label(i + c[0], c[1], c[2]);
        if kind == ElementKind::Hex8 {
            m.set_element(next, kind, corners.iter().map(|c| at(*c)).collect()).unwrap();
            next += 1;
            continue;
        }
        for perm in [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]] {
            let mut v = [[0usize; 3]; 4];
            for s in 0..3 {
                v[s + 1] = v[s];
                v[s + 1][perm[s]] = 1;
            }
            let x = v.map(|c| c.map(|k| k as f64));
            let e = |a: usize| [0, 1, 2].map(|k| x[a][k] - x[0][k]);
            let (a, b, c) = (e(1), e(2), e(3));
            let vol = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
            if vol < 0.0 {
                v.swap(1, 2);
            }
            m.set_element(next, kind, v.iter().map(|c| at(*c)).collect()).unwrap();
            next += 1;
        }
    }
    let all: Vec<u32> = m.elements().keys().copied().collect();
    let nodes: Vec<u32> = m.nodes().keys().copied().collect();
    let pick = |f: &dyn Fn([f64; 3]) -> bool| -> Vec<u32> { nodes.iter().copied().filter(|l| f(m.node(*l).unwrap().xyz)).collect() };
    let (x0, x1, y0, z0) = (pick(&|p| p[0] == 0.0), pick(&|p| p[0] == LENGTH), pick(&|p| p[1] == 0.0), pick(&|p| p[2] == 0.0));
    m.add_set("all", Set { elements: all, nodes: nodes.clone(), ..Set::default() }).unwrap();
    for (name, set) in [("x0", x0), ("x1", x1), ("y0", y0), ("z0", z0)] {
        m.add_set(name, Set { nodes: set, ..Set::default() }).unwrap();
    }
    let thermal = Thermal { expansion: Some(Directional::Isotropic(ALPHA)), ..Thermal::default() };
    m.add_material(
        "Mat",
        Material {
            elastic: Some(Elastic::Isotropic { young: YOUNG, poisson: 0.0 }),
            density: Some(1.0),
            thermal: Some(thermal),
            plastic,
            ..Material::default()
        },
    )
    .unwrap();
    let integration = if kind == ElementKind::Hex8 { Integration::Reduced } else { Integration::Full };
    let formulation = Formulation { integration, ..Formulation::default() };
    m.add_section("bar", Section { set: "all".into(), material: "Mat".into(), kind: SectionKind::Solid {}, formulation }).unwrap();
    let mut step = Step::new("Heat", Procedure::ExplicitDynamic(ExplicitDynamic::new(20.0)));
    let mut bcs = vec![("x0", Symmetry::X), ("y0", Symmetry::Y), ("z0", Symmetry::Z)];
    if clamped {
        bcs.push(("x1", Symmetry::X));
    }
    for (set, plane) in bcs {
        let name = format!("hold_{set}");
        m.add_boundary_condition(&name, BoundaryCondition::Symmetry { target: Target::Set(set.into()), plane }).unwrap();
        step.boundary_conditions.insert(name);
    }
    // DT above the reference temperature (0 °C by default; temperatures are kelvin).
    m.add_load("warm", Load::Temperature { target: Target::Set("all".into()), value: fieldcraft_model::units::celsius_to_kelvin(DT) }).unwrap();
    step.loads.insert("warm".into());
    m.add_step(step, None).unwrap();
    m
}

#[test]
fn a_free_bar_expands_without_stress() {
    for kind in [ElementKind::Hex8, ElementKind::Tet4] {
        let r = solve_explicit(&bar(kind, false, None), None, &Limits::default()).unwrap();
        let last = r.frames.last().unwrap();
        let tip = r.nodes.iter().position(|l| *l == 5).unwrap();
        let want = ALPHA * DT * LENGTH;
        assert!((last.displacement[tip][0] - want).abs() < 0.02 * want, "{kind:?}: {:?} vs {want}", last.displacement[tip]);
        let worst = last.stress.iter().flatten().fold(0.0f64, |m, s| m.max(s.abs()));
        assert!(worst < 0.02 * YOUNG * ALPHA * DT, "{kind:?}: {worst}");
    }
}

#[test]
fn a_clamped_bar_carries_minus_e_alpha_dt_up_to_yield() {
    let sigma = -YOUNG * ALPHA * DT;
    for kind in [ElementKind::Hex8, ElementKind::Tet4] {
        let r = solve_explicit(&bar(kind, true, None), None, &Limits::default()).unwrap();
        for s in &r.frames.last().unwrap().stress {
            assert!((s[0] - sigma).abs() < 0.02 * sigma.abs(), "{kind:?}: {s:?}");
        }
        // Perfectly plastic at half that stress: the bar yields and the stress stays at −σ_y.
        let plastic = Plastic::Bilinear { yield_stress: 0.5 * sigma.abs(), hardening_modulus: 0.0, hardening: Hardening::Isotropic };
        let r = solve_explicit(&bar(kind, true, Some(plastic)), None, &Limits::default()).unwrap();
        let last = r.frames.last().unwrap();
        for (s, ep) in last.stress.iter().zip(&last.plastic_strain) {
            let mises = fieldcraft_fem::von_mises(s);
            assert!((mises - 0.5 * sigma.abs()).abs() < 0.02 * sigma.abs(), "{kind:?}: {s:?}");
            assert!(*ep > 0.0, "{kind:?}: plastic strain {ep}");
        }
        // A bilinear material follows no temperature: none is reported (not 0 K).
        assert!(!last.temperature.is_empty() && last.temperature.iter().all(|t| t.is_nan()), "{kind:?}: {:?}", last.temperature);
    }
}

#[test]
fn orthotropic_expansion_and_heat_loads_are_refused() {
    let mut m = bar(ElementKind::Hex8, false, None);
    let thermal = Thermal { expansion: Some(Directional::Orthotropic([1e-5, 2e-5, 1e-5])), ..Thermal::default() };
    m.replace_material(
        "Mat",
        Material {
            elastic: Some(Elastic::Isotropic { young: YOUNG, poisson: 0.0 }),
            density: Some(1.0),
            thermal: Some(thermal),
            ..Material::default()
        },
    )
    .unwrap();
    assert!(solve_explicit(&m, None, &Limits::default()).unwrap_err().to_string().contains("orthotropic thermal expansion"));
}
