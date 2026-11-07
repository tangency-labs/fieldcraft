//! Analytic checks of static steps that start from initial fields (#160): a residual stress
//! the supports hold stays where it is and moves nothing; a residual stress linear through the
//! thickness of a free beam is released entirely into a uniform curvature, as is an initial
//! temperature linear through the thickness when the beam is cooled to a uniform temperature
//! (a compatible eigenstrain is stress-free: S. P. Timoshenko, J. N. Goodier, *Theory of
//! Elasticity*, 3rd ed., McGraw-Hill 1970, §150).

mod common;

use std::collections::BTreeMap;

use common::*;
use fieldcraft_model::*;
use fieldcraft_solver::{Limits, StaticResult, solve_frequency, solve_static};

const E: f64 = 2.0e11;
const NU: f64 = 0.3;
const L: f64 = 0.4;
const B: f64 = 0.04;
const H: f64 = 0.02;

fn beam(kind: ElementKind, formulation: Formulation) -> Grid {
    let mut g = grid([8, 2, 2], kind, |p| [L * p[0], B * p[1], H * p[2]]);
    g.define(iso(E, NU), None, formulation, Procedure::Static {});
    g
}

/// Hold the beam statically determinately (3-2-1 at its x = 0 end), so it deforms freely.
fn hold_321(g: &mut Grid) {
    let point = |name: &str, at: [f64; 3], values: [Option<f64>; 6], g: &mut Grid| {
        g.node_set(name, |p| near(p[0], at[0]) && near(p[1], at[1]) && near(p[2], at[2]));
        g.bc(name, BoundaryCondition::Displacement { target: Target::Set(name.into()), values });
    };
    let held = Some(0.0);
    point("a", [0.0, 0.0, 0.0], [held, held, held, None, None, None], g);
    point("b", [0.0, 1.0, 0.0], [held, None, held, None, None, None], g);
    point("c", [0.0, 0.0, 1.0], [held, None, None, None, None, None], g);
}

fn stress_field(g: &mut Grid, name: &str, f: impl Fn([f64; 3]) -> [f64; 6]) {
    let values: BTreeMap<u32, [f64; 6]> = g.model.nodes().iter().map(|(l, n)| (*l, f(n.xyz))).collect();
    g.model.set_initial_field(name, InitialField::Stress { values, source: String::new() }).unwrap();
}

fn displacement(r: &StaticResult, label: u32) -> [f64; 3] {
    r.displacement[r.nodes.binary_search(&label).unwrap()]
}

/// The curvature w″ of the beam's mid-line (y = B/2, z = H/2) from its deflection at both ends
/// and the middle: w(0) − 2 w(L/2) + w(L) = w″ L²/4 for a uniform curvature.
fn curvature(g: &Grid, r: &StaticResult) -> f64 {
    let w = |x: f64| displacement(r, g.nodes_where(|p| near(p[0], x) && near(p[1], 0.5) && near(p[2], 0.5))[0])[2];
    (w(0.0) - 2.0 * w(0.5) + w(1.0)) / (L * L / 4.0)
}

fn max_point_stress(r: &StaticResult) -> f64 {
    r.elements.iter().flat_map(|e| &e.points).flat_map(|p| p.stress).fold(0.0, |m, s| m.max(s.abs()))
}

#[test]
fn a_self_equilibrated_residual_stress_held_at_the_ends_moves_nothing() {
    // σ_xx = s (1 − 12 ζ²), ζ = z/H − ½: no net force or moment over the section, no traction
    // on the free faces, ∇·σ = 0. Holding u_x at both ends (the stress's own support), the
    // exact solution is u = 0 with σ = σ₀, and the ends react no net force. (The hexahedra's
    // product rules integrate the load exactly here; the tet10 four-point rule only nearly.)
    let s = 50e6;
    let sigma = |x: [f64; 3]| {
        let z = x[2] / H - 0.5;
        [s * (1.0 - 12.0 * z * z), 0.0, 0.0, 0.0, 0.0, 0.0]
    };
    for (kind, f, name) in [
        (ElementKind::Hex8, full(), "hex8"),
        (ElementKind::Hex8, incompatible(), "hex8 incompatible"),
        (ElementKind::Hex20, full(), "hex20"),
        (ElementKind::Hex20, reduced(), "hex20 reduced"),
    ] {
        let mut g = beam(kind, f);
        g.node_set("ends", |p| near(p[0], 0.0) || near(p[0], 1.0));
        g.bc("ends", BoundaryCondition::Displacement { target: Target::Set("ends".into()), values: [Some(0.0), None, None, None, None, None] });
        g.node_set("a", |p| near(p[0], 0.0) && near(p[1], 0.0) && near(p[2], 0.0));
        g.bc("a", BoundaryCondition::Displacement { target: Target::Set("a".into()), values: [None, Some(0.0), Some(0.0), None, None, None] });
        g.node_set("b", |p| near(p[0], 0.0) && near(p[1], 1.0) && near(p[2], 0.0));
        g.bc("b", BoundaryCondition::Displacement { target: Target::Set("b".into()), values: [None, None, Some(0.0), None, None, None] });
        stress_field(&mut g, "residual", sigma);
        let r = solve_static(&g.model, None, &Limits::default()).unwrap();
        let moved = r.displacement.iter().flatten().fold(0.0f64, |m, u| m.max(u.abs()));
        assert!(moved < 1e-9 * s / E * L, "{name}: moved {moved:e}");
        // Quadratic elements interpolate σ₀ exactly: the stress at every point is σ₀ there.
        if kind == ElementKind::Hex20 {
            for p in r.elements.iter().flat_map(|e| &e.points) {
                let want = sigma(p.xyz);
                assert!(p.stress.iter().zip(want).all(|(a, b)| (a - b).abs() < 1e-6 * s), "{name} at {:?}: {:?} vs {want:?}", p.xyz, p.stress);
            }
        }
        assert!(r.reaction_total.iter().all(|f| f.abs() < 1e-9 * s * B * H), "{name}: {:?}", r.reaction_total);
        assert!(r.warnings.iter().all(|w| !w.contains("initial")), "{name}: {:?}", r.warnings);
    }
}

#[test]
fn a_free_beam_releases_a_linear_residual_stress_into_a_uniform_curvature() {
    // σ_xx = s (2z/H − 1) is the stress of the eigenstrain ε* = −C σ₀, linear in z and so
    // compatible: released freely it leaves no stress and bends the beam to w″ = 2 s / (E H).
    let s = 30e6;
    for (kind, f, name) in [(ElementKind::Hex20, full(), "hex20"), (ElementKind::Hex8, incompatible(), "hex8 incompatible")] {
        let mut g = beam(kind, f);
        hold_321(&mut g);
        stress_field(&mut g, "residual", |x| [s * (2.0 * x[2] / H - 1.0), 0.0, 0.0, 0.0, 0.0, 0.0]);
        let r = solve_static(&g.model, None, &Limits::default()).unwrap();
        let want = 2.0 * s / (E * H);
        let got = curvature(&g, &r);
        assert!((got - want).abs() < 1e-8 * want, "{name}: {got} vs {want}");
        assert!(max_point_stress(&r) < 1e-6 * s, "{name}: {}", max_point_stress(&r));
    }
}

#[test]
fn thermal_strain_is_measured_from_the_initial_temperature() {
    // Ejected with T₀ = 120 + g z (a hotter top), cooled to a uniform 20: ε_th = α (20 − T₀) is
    // linear in z, so the free beam bends stress-free to w″ = α g (top shrinks more). The
    // material's reference temperature (0) no longer matters.
    let (alpha, gradient) = (1.2e-5, 2000.0);
    let mut g = beam(ElementKind::Hex20, full());
    g.model
        .replace_material(
            "Mat",
            Material {
                elastic: Some(iso(E, NU)),
                thermal: Some(Thermal { expansion: Some(Directional::Isotropic(alpha)), ..Thermal::default() }),
                ..Material::default()
            },
        )
        .unwrap();
    hold_321(&mut g);
    let all: Vec<u32> = g.model.nodes().keys().copied().collect();
    g.model.add_set("nodes", Set { nodes: all, ..Set::default() }).unwrap();
    g.load("room", Load::Temperature { target: Target::Set("nodes".into()), value: 20.0 });
    let values = g.model.nodes().iter().map(|(l, n)| (*l, 120.0 + gradient * n.xyz[2])).collect();
    g.model.set_initial_field("ejected", InitialField::Temperature { values, source: "mould".into() }).unwrap();
    let r = solve_static(&g.model, None, &Limits::default()).unwrap();
    let want = alpha * gradient;
    let got = curvature(&g, &r);
    assert!((got - want).abs() < 1e-8 * want, "{got} vs {want}");
    assert!(max_point_stress(&r) < 1e-6 * E * alpha * gradient * H, "{}", max_point_stress(&r));
    // The mean shrinkage is α (20 − T₀) at mid-height, along the length.
    let ux = |x: f64| displacement(&r, g.nodes_where(|p| near(p[0], x) && near(p[1], 0.5) && near(p[2], 0.5))[0])[0];
    let strain = (ux(1.0) - ux(0.0)) / L;
    let mean = alpha * (20.0 - (120.0 + gradient * H / 2.0));
    assert!((strain - mean).abs() < 1e-9 * mean.abs(), "{strain} vs {mean}");
}

#[test]
fn steps_and_elements_that_take_no_initial_fields_say_so() {
    let mut g = beam(ElementKind::Hex8, full());
    hold_321(&mut g);
    stress_field(&mut g, "residual", |_| [1e6, 0.0, 0.0, 0.0, 0.0, 0.0]);
    // A field naming a node the mesh lacks: left out, with a warning.
    g.model.set_initial_field("stale", InitialField::Temperature { values: [(99_999, 300.0)].into(), source: String::new() }).unwrap();
    let r = solve_static(&g.model, None, &Limits::default()).unwrap();
    assert!(r.warnings.iter().any(|w| w.contains("1 initial field values are at nodes the mesh does not have")), "{:?}", r.warnings);
    g.model.add_material("Steel", Material { elastic: Some(iso(E, NU)), density: Some(7800.0), ..Material::default() }).unwrap();
    g.model
        .replace_section(
            "solid",
            Section { orientation: None, set: "all".into(), material: "Steel".into(), kind: SectionKind::Solid {}, formulation: full() },
        )
        .unwrap();
    g.model.add_step(Step::new("Modes", Procedure::Frequency { modes: 3, min_hz: None, max_hz: None }), None).unwrap();
    let modes = solve_frequency(&g.model, Some("Modes"), &Default::default(), &Limits::default()).unwrap();
    assert!(modes.warnings.iter().any(|w| w.contains("step `Modes` starts without them")), "{:?}", modes.warnings);
    // Two fields disagreeing at a node are refused.
    let first = *g.model.nodes().keys().next().unwrap();
    g.model
        .set_initial_field("other", InitialField::Stress { values: [(first, [2e6, 0.0, 0.0, 0.0, 0.0, 0.0])].into(), source: String::new() })
        .unwrap();
    assert!(solve_static(&g.model, Some("Step"), &Limits::default()).unwrap_err().to_string().contains("is given both"));
}
