//! Hybrid (mixed u–p) solids for nearly incompressible materials, against closed forms: a block
//! in uniaxial strain (exact for every element) and Lamé's thick cylinder in plane strain at
//! ν = 0.4999, where the displacement elements lock and the hybrid ones do not. Run with
//! `--nocapture` to see the errors.

mod common;

use common::*;
use fieldcraft_model::*;
use fieldcraft_solver::{Limits, solve_static};

const NU: f64 = 0.4999;

fn hybrid() -> Formulation {
    Formulation { hybrid: true, ..Formulation::default() }
}

fn hybrids() -> Vec<(ElementKind, &'static str)> {
    vec![(ElementKind::Hex8, "hex8 hybrid"), (ElementKind::Hex20, "hex20 hybrid"), (ElementKind::Tet10, "tet10 hybrid")]
}

#[test]
fn block_in_uniaxial_strain() {
    // Held normal to the bottom and the four sides, pressed on top: ε_zz = −p/M with the
    // constrained modulus M = E(1 − ν)/((1 + ν)(1 − 2ν)); σ_xx = σ_yy = −p ν/(1 − ν).
    let (e, p) = (10.0, 1.0);
    let m = e * (1.0 - NU) / ((1.0 + NU) * (1.0 - 2.0 * NU));
    for (kind, name) in hybrids() {
        let mut g = grid([2, 2, 2], kind, |x| x);
        g.define(iso(e, NU), None, hybrid(), Procedure::Static {});
        g.symmetry("x0", Symmetry::X, |x| near(x[0], 0.0));
        g.symmetry("x1", Symmetry::X, |x| near(x[0], 1.0));
        g.symmetry("y0", Symmetry::Y, |x| near(x[1], 0.0));
        g.symmetry("y1", Symmetry::Y, |x| near(x[1], 1.0));
        g.symmetry("z0", Symmetry::Z, |x| near(x[2], 0.0));
        g.face_set("top", |x| near(x[2], 1.0));
        g.load("press", Load::Pressure { set: "top".into(), magnitude: p });
        let r = solve_static(&g.model, None, &Limits::default()).unwrap();
        for (i, l) in r.nodes.iter().enumerate() {
            let z = g.model.node(*l).unwrap().xyz[2];
            assert!((r.displacement[i][2] + p * z / m).abs() < 1e-9 * p / m, "{name} node {l}: {:?}", r.displacement[i]);
        }
        for s in &r.stress {
            assert!((s[2] + p).abs() < 1e-6 && (s[0] + p * NU / (1.0 - NU)).abs() < 1e-6 && (s[1] - s[0]).abs() < 1e-6, "{name}: {s:?}");
        }
    }
}

/// The inner radial displacement of a quarter of a thick cylinder a = 1, b = 2 in plane strain
/// under internal pressure 1, and Lamé's value.
fn cylinder(kind: ElementKind, formulation: Formulation, n: [usize; 3]) -> (f64, f64) {
    let (e, p, a, b) = (1000.0, 1.0, 1.0, 2.0);
    let mut g = grid(n, kind, |x| {
        let (r, t) = (a + (b - a) * x[0], std::f64::consts::FRAC_PI_2 * x[1]);
        [r * t.cos(), r * t.sin(), 0.2 * x[2]]
    });
    g.define(iso(e, NU), None, formulation, Procedure::Static {});
    g.symmetry("y0", Symmetry::Y, |x| near(x[1], 0.0));
    g.symmetry("x0", Symmetry::X, |x| near(x[1], 1.0));
    g.symmetry("z0", Symmetry::Z, |x| near(x[2], 0.0));
    g.symmetry("z1", Symmetry::Z, |x| near(x[2], 1.0));
    g.face_set("bore", |x| near(x[0], 0.0));
    g.load("pressure", Load::Pressure { set: "bore".into(), magnitude: p });
    let r = solve_static(&g.model, None, &Limits::default()).unwrap();
    let inner = g.nodes_where(|x| near(x[0], 0.0) && near(x[1], 0.0) && near(x[2], 0.0))[0];
    let u = r.displacement[r.nodes.binary_search(&inner).unwrap()][0];
    let exact = (1.0 + NU) * p * a * a / (e * (b * b - a * a)) * ((1.0 - 2.0 * NU) * a + b * b / a);
    (u, exact)
}

#[test]
fn nearly_incompressible_thick_cylinder_does_not_lock() {
    for (kind, name, n, tol) in [
        (ElementKind::Hex8, "hex8 hybrid", [8, 12, 1], 0.01),
        (ElementKind::Hex20, "hex20 hybrid", [4, 6, 1], 1e-4),
        (ElementKind::Tet10, "tet10 hybrid", [4, 6, 1], 5e-3),
    ] {
        let (u, exact) = cylinder(kind, hybrid(), n);
        eprintln!("{name}: u(a) = {u:.6e}, Lamé {exact:.6e}, error {:.2e}", (u - exact) / exact);
        assert!(((u - exact) / exact).abs() < tol, "{name}: {u} vs {exact}");
    }
    // The displacement hex8 locks: far too stiff on the same mesh.
    let (u, exact) = cylinder(ElementKind::Hex8, Formulation::default(), [8, 12, 1]);
    eprintln!("hex8 (displacement): u(a) = {u:.6e}, Lamé {exact:.6e}, error {:.2e}", (u - exact) / exact);
    assert!(u < 0.7 * exact, "full hex8 should lock at ν = {NU}: {u} vs {exact}");
}

/// Why a one-cell hybrid model of `kind` and `elastic` cannot be solved.
fn refusal(kind: ElementKind, elastic: Elastic) -> String {
    let mut g = grid([1, 1, 1], kind, |x| x);
    g.define(elastic, None, hybrid(), Procedure::Static {});
    g.fixed("base", |x| near(x[2], 0.0));
    g.face_set("top", |x| near(x[2], 1.0));
    g.load("press", Load::Pressure { set: "top".into(), magnitude: 1.0 });
    solve_static(&g.model, None, &Limits::default()).unwrap_err().to_string()
}

#[test]
fn hybrid_refusals() {
    for (kind, needle) in [(ElementKind::Tet4, "hybrid tet4"), (ElementKind::Wedge6, "hybrid wedge6")] {
        let e = refusal(kind, iso(10.0, NU));
        assert!(e.contains(needle), "{e}");
    }
    let ortho = Elastic::EngineeringConstants { young: [10.0, 20.0, 30.0], poisson: [0.3, 0.2, 0.1], shear: [5.0, 6.0, 7.0] };
    let e = refusal(ElementKind::Hex8, ortho);
    assert!(e.contains("orthotropic material `Mat`"), "{e}");
}

#[test]
fn heated_confined_block_carries_the_bulk_thermal_stress() {
    // Held normal to all six faces and heated by ΔT: σ = −3 K α ΔT in every direction, the
    // pressure the hybrid elements carry as their own unknown.
    let (e, alpha, dt) = (10.0, 1e-4, 30.0);
    let bulk = e / (3.0 * (1.0 - 2.0 * NU));
    for (kind, name) in hybrids() {
        let mut g = grid([2, 2, 2], kind, |x| x);
        g.define(iso(e, NU), None, hybrid(), Procedure::Static {});
        let thermal = Thermal { expansion: Some(Directional::Isotropic(alpha)), reference_temperature: 0.0, ..Thermal::default() };
        g.model.replace_material("Mat", Material { elastic: Some(iso(e, NU)), thermal: Some(thermal), ..Material::default() }).unwrap();
        for (set, plane, k, v) in [
            ("x0", Symmetry::X, 0, 0.0),
            ("x1", Symmetry::X, 0, 1.0),
            ("y0", Symmetry::Y, 1, 0.0),
            ("y1", Symmetry::Y, 1, 1.0),
            ("z0", Symmetry::Z, 2, 0.0),
            ("z1", Symmetry::Z, 2, 1.0),
        ] {
            g.symmetry(set, plane, |x| near(x[k], v));
        }
        g.node_set("everything", |_| true);
        g.load("warm", Load::Temperature { target: Target::Set("everything".into()), value: dt });
        let r = solve_static(&g.model, None, &Limits::default()).unwrap();
        let want = -3.0 * bulk * alpha * dt;
        for s in &r.stress {
            for k in 0..3 {
                assert!((s[k] - want).abs() < 1e-6 * want.abs(), "{name}: {s:?} vs {want}");
            }
        }
    }
}
