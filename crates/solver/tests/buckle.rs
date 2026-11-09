//! Linear buckling against closed forms: Euler columns of hex20 solids (pinned–pinned by a
//! symmetry plane at midspan, and clamped–free), with Engesser's shear correction
//! P = P_e / (1 + P_e/(κ G A)) (S. P. Timoshenko, J. M. Gere, *Theory of Elastic Stability*,
//! 2nd ed., McGraw-Hill 1961, §2.1, §2.17), the shifted (LDLᵀ) form against the plain one, and
//! what buckling steps refuse.

mod common;

use common::*;
use fieldcraft_model::*;
use fieldcraft_solver::{Limits, solve_buckle};

const E: f64 = 1000.0;
const B: f64 = 1.0;
const L: f64 = 20.0;

fn buckle(modes: u32, shift: Option<f64>) -> Procedure {
    Procedure::Buckle { modes, shift, follower_pressure: false }
}

/// Engesser's critical load of a square column of side B and effective length `le`, ν = 0.
fn engesser(le: f64) -> f64 {
    let (a, i, g) = (B * B, B.powi(4) / 12.0, E / 2.0);
    let pe = std::f64::consts::PI.powi(2) * E * i / (le * le);
    pe / (1.0 + pe / (5.0 / 6.0 * g * a))
}

/// A column of `length` along z, 2 × 2 × `cells` hex20, pressed by a unit pressure on its top.
fn column(length: f64, cells: usize, shift: Option<f64>) -> Grid {
    let mut g = grid([2, 2, cells], ElementKind::Hex20, |x| [B * x[0], B * x[1], length * x[2]]);
    g.define(iso(E, 0.0), None, Formulation::default(), buckle(3, shift));
    g.face_set("top", |x| near(x[2], 1.0));
    g.load("press", Load::Pressure { set: "top".into(), magnitude: 1.0 });
    g
}

#[test]
fn euler_columns() {
    // Pinned–pinned: half the column, the pinned end (z = 0 here is the loaded top) held
    // laterally over its face (which leaves its rotation free), midspan a symmetry plane.
    let mut g = grid([2, 2, 10], ElementKind::Hex20, |x| [B * x[0], B * x[1], L / 2.0 * x[2]]);
    g.define(iso(E, 0.0), None, Formulation::default(), buckle(2, None));
    g.face_set("end", |x| near(x[2], 0.0));
    g.load("press", Load::Pressure { set: "end".into(), magnitude: 1.0 });
    g.node_set("pin", |x| near(x[2], 0.0));
    g.bc("pin", BoundaryCondition::Displacement { target: Target::Set("pin".into()), values: [Some(0.0), Some(0.0), None, None, None, None] });
    g.symmetry("mid", Symmetry::Z, |x| near(x[2], 1.0));
    let r = solve_buckle(&g.model, Some("Step"), &Limits::default()).unwrap();
    let want = engesser(L) / (B * B);
    let got = r.modes[0].eigenvalue;
    eprintln!("pinned: λ = {got:.6} (Engesser {want:.6}, {:+.3} %), second {:.6}", 100.0 * (got / want - 1.0), r.modes[1].eigenvalue);
    assert!((got / want - 1.0).abs() < 5e-3);
    assert!((r.modes[1].eigenvalue / got - 1.0).abs() < 1e-6, "the square section's two planes");
    assert!(r.buckling.as_ref().unwrap().negative.is_empty() || r.buckling.as_ref().unwrap().negative[0] < -got);

    // Clamped–free: the base held, the top free: le = 2L.
    let mut g = column(L / 2.0, 10, None);
    g.fixed("base", |x| near(x[2], 0.0));
    let r = solve_buckle(&g.model, Some("Step"), &Limits::default()).unwrap();
    let want = engesser(L) / (B * B);
    let got = r.modes[0].eigenvalue;
    eprintln!("clamped–free: λ = {got:.6} (Engesser {want:.6}, {:+.3} %)", 100.0 * (got / want - 1.0));
    assert!((got / want - 1.0).abs() < 5e-3);
    // Scaled to a largest component of 1, which is lateral at the free top.
    let big = r.modes[0].shape.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
    assert!((big - 1.0).abs() < 1e-12);

    // Shifted between the first and the higher modes (LDLᵀ of K + σ K_σ): its inertia counts
    // the two load factors below the shift, and the modes above it come out.
    let mut s = column(L / 2.0, 10, Some(2.0 * got));
    s.fixed("base", |x| near(x[2], 0.0));
    let shifted = solve_buckle(&s.model, Some("Step"), &Limits::default()).unwrap();
    let info = shifted.buckling.as_ref().unwrap();
    assert_eq!(info.below_shift, Some(2));
    let third = r.modes[2].eigenvalue;
    assert!((shifted.modes[0].eigenvalue / third - 1.0).abs() < 1e-7, "{} vs {third}", shifted.modes[0].eigenvalue);
}

#[test]
fn buckling_steps_refuse_what_they_cannot_do() {
    // No load: nothing to buckle.
    let mut g = column(L, 4, None);
    g.model.set_load_active("Step", "press", false).unwrap();
    g.fixed("base", |x| near(x[2], 0.0));
    assert!(solve_buckle(&g.model, Some("Step"), &Limits::default()).is_err());
    // Unsupported: a located singular stiffness.
    let g = column(L, 4, None);
    let e = solve_buckle(&g.model, Some("Step"), &Limits::default()).unwrap_err();
    assert!(e.to_string().contains("singular"), "{e}");
    // A shift that is not finite, and a step of another kind.
    assert!(Procedure::Buckle { modes: 0, shift: None, follower_pressure: false }.validate().is_err());
    assert!(Procedure::Buckle { modes: 1, shift: Some(f64::NAN), follower_pressure: false }.validate().is_err());
}

/// A clamped–free column given its first buckling mode as an imperfection δ₀ at the top, under a
/// dead end force of half its buckling load: the top's lateral deflection grows to
/// δ₀ / (1 − P/P_cr) = 2 δ₀ (Timoshenko & Gere 1961, §1.12, for the mode-shaped imperfection).
#[test]
fn an_imperfection_from_the_buckling_mode() {
    let mut g = column(L / 2.0, 10, None);
    g.fixed("base", |x| near(x[2], 0.0));
    let critical = solve_buckle(&g.model, Some("Step"), &Limits::default()).unwrap().modes[0].eigenvalue * B * B;
    let delta0 = 0.01;
    g.node_set("tip", |x| near(x[2], 1.0));
    g.model
        .add_load("push", Load::Force { target: Target::Set("tip".into()), vector: [0.0, 0.0, -0.5 * critical], distribution: Distribution::Total })
        .unwrap();
    let settings = NonlinearStatic { imperfection: Some(Imperfection { step: "Step".into(), mode: 1, scale: delta0 }), ..NonlinearStatic::default() };
    let mut nl = Step::new("Push", Procedure::NonlinearStatic(settings));
    nl.boundary_conditions.insert("base".into());
    nl.loads.insert("push".into());
    g.model.add_step(nl, None).unwrap();
    let r = fieldcraft_solver::solve_nonlinear(&g.model, Some("Push"), &Limits::default()).unwrap();
    assert!(r.warnings.iter().any(|w| w.contains("imperfect geometry")), "{:?}", r.warnings);
    let lateral = r.displacement.iter().map(|u| u[0].hypot(u[1])).fold(0.0, f64::max);
    // The deflection added to δ₀: δ₀ α/(1 − α) = δ₀ at α = ½.
    eprintln!("imperfection {delta0}: added top deflection {lateral:.6e}");
    assert!((lateral / delta0 - 1.0).abs() < 0.05, "{lateral}");
    // A mode the buckling step does not have, or a step that does not buckle.
    let mut bad = g.model.clone();
    let mut s = bad.step("Push").unwrap().clone();
    if let Procedure::NonlinearStatic(n) = &mut s.procedure {
        n.imperfection = Some(Imperfection { step: "Step".into(), mode: 9, scale: delta0 });
    }
    bad.replace_step(s).unwrap();
    assert!(fieldcraft_solver::solve_nonlinear(&bad, Some("Push"), &Limits::default()).is_err());
}
