//! Reference points, couplings and moments against closed-form solutions. Run with
//! `--nocapture` to see the errors.

mod common;

use std::f64::consts::PI;

use common::*;
use fieldcraft_model::{
    BoundaryCondition, Coupling, CouplingKind, Distribution, Dofs, ElementKind, Load, Procedure, ReferencePoint, Set, Symmetry, Target,
};
use fieldcraft_solver::{FrequencyOptions, Limits, StaticResult, solve_frequency, solve_static};

const ALL: [u8; 6] = [1, 2, 3, 4, 5, 6];

fn solve(g: &Grid) -> StaticResult {
    solve_static(&g.model, Some("Step"), &Limits::default()).unwrap()
}

fn at(r: &StaticResult, label: u32) -> usize {
    r.nodes.binary_search(&label).unwrap()
}

fn rp(name: &str) -> Target {
    Target::ReferencePoint(name.into())
}

fn force(target: Target, vector: [f64; 3]) -> Load {
    Load::Force { target, vector, distribution: Distribution::default() }
}

fn moment(target: Target, vector: [f64; 3]) -> Load {
    Load::Moment { target, vector, distribution: Distribution::default() }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/// A block held only through a reference point coupled to one face, the point given a small
/// rigid motion: the whole block follows it rigidly, u = u_r + θ × (x − x_r), with no stress.
/// Kinematic couplings eliminate the face nodes; a distributing coupling with a prescribed point
/// eliminates face nodes chosen by the constraint itself (the point is dependent).
#[test]
fn prescribed_reference_point_moves_the_block_rigidly() {
    let motion = [2e-3, -1e-3, 3e-3, 1e-3, -2e-3, 1.5e-3];
    let xr = [0.0, 0.6, 0.4];
    for (kind, formulation, name) in [
        (ElementKind::Hex8, incompatible(), "hex8 incompatible"),
        (ElementKind::Hex20, full(), "hex20"),
        (ElementKind::Tet10, full(), "tet10"),
        (ElementKind::Wedge6, full(), "wedge6"),
    ] {
        for coupling in [CouplingKind::Kinematic, CouplingKind::Distributing] {
            let mut g = grid([3, 2, 2], kind, |q| [3.0 * q[0], 1.2 * q[1], 0.8 * q[2]]);
            g.define(iso(1000.0, 0.3), None, formulation, Procedure::Static {});
            g.couple("RP", xr, coupling, &ALL, |q| near(q[0], 0.0));
            let values = motion.map(Some);
            g.bc("move", BoundaryCondition::Displacement { target: rp("RP"), values });
            let r = solve(&g);
            let (t, w) = ([motion[0], motion[1], motion[2]], [motion[3], motion[4], motion[5]]);
            let mut worst: f64 = 0.0;
            for (label, u) in r.nodes.iter().zip(&r.displacement) {
                let x = g.model.node(*label).unwrap().xyz;
                let rot = cross(w, [x[0] - xr[0], x[1] - xr[1], x[2] - xr[2]]);
                worst = worst.max((0..3).map(|k| (u[k] - t[k] - rot[k]).abs()).fold(0.0, f64::max));
            }
            let stress = r.stress.iter().flatten().fold(0.0f64, |m, s| m.max(s.abs()));
            let p = &r.reference_points[0];
            eprintln!("rigid {name} {coupling:?}: displacement error {worst:.1e}, largest stress {stress:.1e}");
            assert!(worst < 1e-13, "{name} {coupling:?}: {worst}");
            assert!(stress < 1e-9, "{name} {coupling:?}: {stress}");
            assert_eq!(p.displacement, motion);
            assert!(p.reaction.iter().all(|x| x.abs() < 1e-9), "{:?}", p.reaction);
        }
    }
}

/// A cantilever clamped through a fixed reference point kinematically coupled to its root face
/// is the cantilever clamped directly: same displacements and natural frequencies, and the
/// point's reaction carries the whole load (force and moment about the point).
#[test]
fn clamping_through_a_coupled_reference_point_is_clamping() {
    let build = |through_point: bool, procedure: Procedure| {
        let mut g = grid([8, 1, 1], ElementKind::Hex20, |q| [8.0 * q[0], q[1], q[2]]);
        g.define(iso(1000.0, 0.3), Some(1.0), full(), procedure);
        if through_point {
            g.couple("Root", [0.0, 0.5, 0.5], CouplingKind::Kinematic, &[1, 2, 3], |q| near(q[0], 0.0));
            g.bc("clamp", BoundaryCondition::Fixed { target: rp("Root") });
        } else {
            g.fixed("clamp", |q| near(q[0], 0.0));
        }
        g
    };
    let (mut direct, mut coupled) = (build(false, Procedure::Static {}), build(true, Procedure::Static {}));
    for g in [&mut direct, &mut coupled] {
        g.face_set("end", |q| near(q[0], 1.0));
        g.load("shear", Load::Traction { set: "end".into(), vector: [0.0, 0.0, -1.0] });
    }
    let (a, b) = (solve(&direct), solve(&coupled));
    let scale = a.displacement.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
    let diff = a.displacement.iter().flatten().zip(b.displacement.iter().flatten()).fold(0.0f64, |m, (x, y)| m.max((x - y).abs()));
    let p = &b.reference_points[0];
    eprintln!("clamp through a point: displacement difference {:.1e} (relative), reaction {:?}", diff / scale, p.reaction);
    assert!(diff < 1e-10 * scale, "{diff}");
    // The tip load −1 in z at x = 8 (centroid line y = z = 0.5): reaction +1 in z and the moment
    // about the point M_y = −8 (r × F with r = (8, 0, 0), F = (0, 0, 1) → (0, −8, 0)).
    assert!((p.reaction[2] - 1.0).abs() < 1e-9 && p.reaction[0].abs() < 1e-9 && p.reaction[1].abs() < 1e-9, "{:?}", p.reaction);
    assert!((p.reaction[4] + 8.0).abs() < 1e-8 && p.reaction[3].abs() < 1e-8 && p.reaction[5].abs() < 1e-8, "{:?}", p.reaction);
    assert!((b.reaction_total[2] - 1.0).abs() < 1e-9, "{:?}", b.reaction_total);

    let modes = Procedure::Frequency { modes: 3, min_hz: None, max_hz: None };
    let (direct, coupled) = (build(false, modes.clone()), build(true, modes));
    let f = |g: &Grid| solve_frequency(&g.model, Some("Step"), &FrequencyOptions::default(), &Limits::default()).unwrap();
    let (fa, fb) = (f(&direct), f(&coupled));
    for (x, y) in fa.modes.iter().zip(&fb.modes) {
        assert!((x.frequency - y.frequency).abs() < 1e-9 * x.frequency, "{} vs {}", x.frequency, y.frequency);
    }
    assert_eq!(fb.reference_points.len(), 1);
    assert!(fb.modes.iter().all(|m| m.point_shape[0] == [0.0; 6]), "the fixed point does not move");
}

/// A bar in tension loaded through a distributing coupling at the centroid of its end face
/// against the same total force as a uniform pressure: away from the end (Saint-Venant) both give
/// σxx = F/A and u_x = F x / (E A); the coupled load reaches the face as equal nodal forces.
#[test]
fn distributing_force_equals_the_equivalent_pressure_away_from_the_end() {
    let (l, b, h, young, nu, total) = (20.0, 1.0, 1.0, 1000.0, 0.3, 4.0);
    let area = b * h;
    let mut errors = Vec::new();
    for (kind, formulation, name) in [(ElementKind::Hex20, full(), "hex20"), (ElementKind::Hex8, incompatible(), "hex8 incompatible")] {
        let mut results = Vec::new();
        for coupled in [false, true] {
            let mut g = grid([20, 2, 2], kind, |q| [l * q[0], b * q[1], h * q[2]]);
            g.define(iso(young, nu), None, formulation, Procedure::Static {});
            g.symmetry("x0", Symmetry::X, |q| near(q[0], 0.0));
            g.symmetry("y0", Symmetry::Y, |q| near(q[1], 0.0));
            g.symmetry("z0", Symmetry::Z, |q| near(q[2], 0.0));
            if coupled {
                // The quarter model's end face: its centroid is at (l, b/2, h/2).
                g.couple("End", [l, b / 2.0, h / 2.0], CouplingKind::Distributing, &ALL, |q| near(q[0], 1.0));
                g.load("pull", force(rp("End"), [total, 0.0, 0.0]));
            } else {
                g.face_set("end", |q| near(q[0], 1.0));
                g.load("pull", Load::Pressure { set: "end".into(), magnitude: -total / area });
            }
            let r = solve(&g);
            // At least eight widths of the quarter (four of the whole bar) from the loaded end.
            let middle: Vec<u32> = g.nodes_where(|q| q[0] > 0.25 && q[0] < 0.6);
            let es = middle.iter().map(|n| (r.stress[at(&r, *n)][0] - total / area).abs() / (total / area)).fold(0.0, f64::max);
            let eu = middle
                .iter()
                .map(|n| {
                    let x = g.model.node(*n).unwrap().xyz[0];
                    let want = total * x / (young * area);
                    (r.displacement[at(&r, *n)][0] - want).abs() / want
                })
                .fold(0.0, f64::max);
            eprintln!("tension {name} {}: mid-length σxx error {:.2e}, u_x error {:.2e}", if coupled { "coupled" } else { "pressure" }, es, eu);
            assert!((r.applied_total[0] - total).abs() < 1e-12 && (r.reaction_total[0] + total).abs() < 1e-9, "{name}");
            errors.push((es, eu));
            results.push(r);
        }
        // The point moves with the face's mean: between the pressure-loaded face's extremes.
        let end = x_range(&results[0]);
        let ux = results[1].reference_points[0].displacement[0];
        assert!(ux > 0.9 * end.0 && ux < 1.1 * end.1, "{ux} {end:?}");
    }
    assert!(errors.iter().all(|(es, eu)| *es < 1e-4 && *eu < 1e-4), "{errors:?}");
}

/// The smallest and largest x displacement.
fn x_range(r: &StaticResult) -> (f64, f64) {
    r.displacement.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), u| (lo.min(u[0]), hi.max(u[0])))
}

/// Pure bending by a moment on a reference point with a distributing coupling on the end face:
/// away from the end σxx = M z / I and the curvature κ = M / (E I) (u_z = κ x² / 2 on the axis).
#[test]
fn bending_moment_through_a_distributing_coupling() {
    let (l, b, h, young, nu, m) = (10.0, 1.0, 2.0, 1000.0, 0.3, 3.0);
    let i = b * h * h * h / 12.0;
    let mut g = grid([10, 2, 4], ElementKind::Hex20, |q| [l * q[0], b * (q[1] - 0.5), h * (q[2] - 0.5)]);
    g.define(iso(young, nu), None, full(), Procedure::Static {});
    g.symmetry("x0", Symmetry::X, |q| near(q[0], 0.0));
    // The pure-bending field has u_y = −ν κ y z and u_z = κ (x² + ν (z² − y²)) / 2: zero u_y on
    // y = 0 at x = 0, and zero u_z at the origin.
    g.node_set("axis", |q| near(q[0], 0.0) && near(q[1], 0.5) && near(q[2], 0.5));
    g.bc("axis", BoundaryCondition::Displacement { target: Target::Set("axis".into()), values: [None, Some(0.0), Some(0.0), None, None, None] });
    g.node_set("edge", |q| near(q[0], 0.0) && near(q[1], 0.5) && near(q[2], 0.0));
    g.bc("edge", BoundaryCondition::Displacement { target: Target::Set("edge".into()), values: [None, Some(0.0), None, None, None, None] });
    g.couple("End", [l, 0.0, 0.0], CouplingKind::Distributing, &ALL, |q| near(q[0], 1.0));
    g.load("bend", moment(rp("End"), [0.0, m, 0.0]));
    let r = solve(&g);
    let mut es: f64 = 0.0;
    for n in g.nodes_where(|q| q[0] > 0.35 && q[0] < 0.65) {
        let z = g.model.node(n).unwrap().xyz[2];
        es = es.max((r.stress[at(&r, n)][0] - m * z / i).abs() / (m * h / 2.0 / i));
    }
    let kappa = m / (young * i);
    let mid = g.nodes_where(|q| near(q[0], 0.5) && near(q[1], 0.5) && near(q[2], 0.5))[0];
    let x = g.model.node(mid).unwrap().xyz[0];
    // Positive M_y stretches +z: the axis bends toward −z.
    let want = -kappa * x * x / 2.0;
    let eu = (r.displacement[at(&r, mid)][2] - want).abs() / want.abs();
    eprintln!("bending through a distributing coupling: σxx error {es:.2e} (of the peak), mid-span deflection error {eu:.2e}");
    assert!(es < 5e-3 && eu < 5e-3, "{es} {eu}");
    for k in 0..3 {
        assert!(r.applied_total[k].abs() < 1e-12 && r.reaction_total[k].abs() < 1e-9, "a moment adds no net force");
    }
}

/// A bar in tension, its end face coupled distributing through the face's own element faces: the
/// consistent nodal loads of the uniform traction are exact, so every node carries σxx = F/A and
/// u_x = F x / (E A) to round-off, and the point moves with the face's mean (Timoshenko and
/// Goodier, *Theory of Elasticity*, 3rd ed., 1970, ch. 1).
#[test]
fn tension_through_a_distributing_coupling_on_faces_is_uniform() {
    let (young, nu, total, area, l) = (1000.0, 0.3, 4.0, 1.0, 4.0);
    for (kind, formulation, name) in [(ElementKind::Hex20, full(), "hex20"), (ElementKind::Hex8, full(), "hex8")] {
        let mut g = grid([4, 2, 2], kind, |q| [l * q[0], q[1], q[2]]);
        g.define(iso(young, nu), None, formulation, Procedure::Static {});
        g.symmetry("x0", Symmetry::X, |q| near(q[0], 0.0));
        g.symmetry("y0", Symmetry::Y, |q| near(q[1], 0.0));
        g.symmetry("z0", Symmetry::Z, |q| near(q[2], 0.0));
        g.couple_faces("End", [l, 0.5, 0.5], CouplingKind::Distributing, &ALL, |q| near(q[0], 1.0));
        g.load("pull", force(rp("End"), [total, 0.0, 0.0]));
        let r = solve(&g);
        let (mut es, mut eu) = (0.0f64, 0.0f64);
        for (i, label) in r.nodes.iter().enumerate() {
            let x = g.model.node(*label).unwrap().xyz[0];
            es = es.max((r.stress[i][0] - total / area).abs() / (total / area));
            eu = eu.max((r.displacement[i][0] - total * x / (young * area)).abs());
            assert!((r.stress[i][0] - total / area).abs() < 1e-9 * total / area, "{name}: σxx at node {label}: {}", r.stress[i][0]);
            assert!((r.displacement[i][0] - total * x / (young * area)).abs() < 1e-9 * total * l / (young * area), "{name}: u_x at node {label}");
        }
        let ux = r.reference_points[0].displacement[0];
        let want = total * l / (young * area);
        eprintln!("tension through faces {name}: σxx error {es:.2e}, u_x error {eu:.2e}, point u_x error {:.2e}", (ux - want).abs() / want);
        assert!((ux - want).abs() < 1e-9 * want, "{name}: point u_x {ux}");
    }
}

/// Pure bending through a distributing coupling on the end face's element faces: the linear
/// traction the fit carries is exact, so σxx = M z / I and the end rotation M L / (E I) hold at
/// every node and at the point.
#[test]
fn bending_through_a_distributing_coupling_on_faces_is_exact() {
    let (l, b, h, young, nu, m) = (10.0, 1.0, 2.0, 1000.0, 0.3, 3.0);
    let i = b * h * h * h / 12.0;
    let mut g = grid([10, 2, 4], ElementKind::Hex20, |q| [l * q[0], b * (q[1] - 0.5), h * (q[2] - 0.5)]);
    g.define(iso(young, nu), None, full(), Procedure::Static {});
    g.symmetry("x0", Symmetry::X, |q| near(q[0], 0.0));
    g.node_set("axis", |q| near(q[0], 0.0) && near(q[1], 0.5) && near(q[2], 0.5));
    g.bc("axis", BoundaryCondition::Displacement { target: Target::Set("axis".into()), values: [None, Some(0.0), Some(0.0), None, None, None] });
    g.node_set("edge", |q| near(q[0], 0.0) && near(q[1], 0.5) && near(q[2], 0.0));
    g.bc("edge", BoundaryCondition::Displacement { target: Target::Set("edge".into()), values: [None, Some(0.0), None, None, None, None] });
    g.couple_faces("End", [l, 0.0, 0.0], CouplingKind::Distributing, &ALL, |q| near(q[0], 1.0));
    g.load("bend", moment(rp("End"), [0.0, m, 0.0]));
    let r = solve(&g);
    let mut es: f64 = 0.0;
    for (k, label) in r.nodes.iter().enumerate() {
        let z = g.model.node(*label).unwrap().xyz[2];
        let err = (r.stress[k][0] - m * z / i).abs();
        es = es.max(err);
        assert!(err < 1e-8 * (m * h / 2.0 / i), "node {label}: σxx error {err:.2e}");
    }
    let want = m * l / (young * i);
    let got = r.reference_points[0].displacement[4];
    eprintln!("bending through faces: σxx error {es:.2e}, point rotation error {:.2e}", (got - want).abs() / want);
    assert!((got - want).abs() < 1e-8 * want, "{got} vs {want}");
}

/// A set with element faces and extra nodes off them: the surface does not reach those nodes, so
/// the coupling refuses rather than leave them out of the fit.
#[test]
fn a_distributing_coupling_refuses_nodes_off_its_faces() {
    let mut g = grid([4, 2, 2], ElementKind::Hex8, |q| [4.0 * q[0], q[1], q[2]]);
    g.define(iso(1000.0, 0.3), None, full(), Procedure::Static {});
    g.symmetry("x0", Symmetry::X, |q| near(q[0], 0.0));
    g.symmetry("y0", Symmetry::Y, |q| near(q[1], 0.0));
    g.symmetry("z0", Symmetry::Z, |q| near(q[2], 0.0));
    let faces = g.faces_where(|q| near(q[0], 1.0));
    let nodes = g.nodes_where(|q| near(q[0], 0.0));
    g.model.add_set("mixed", Set { faces, nodes, ..Set::default() }).unwrap();
    g.model.add_reference_point("End", ReferencePoint { xyz: [4.0, 0.5, 0.5], label: 0 }).unwrap();
    let dofs = Dofs::try_from(ALL.to_vec()).unwrap();
    g.model
        .add_coupling("End_coupling", Coupling { reference_point: "End".into(), set: "mixed".into(), kind: CouplingKind::Distributing, dofs })
        .unwrap();
    g.load("pull", force(rp("End"), [4.0, 0.0, 0.0]));
    let e = solve_static(&g.model, Some("Step"), &Limits::default()).unwrap_err().to_string();
    eprintln!("refusal: {e}");
    assert!(e.contains("on none of the set's element faces"), "{e}");
}

/// A tube of hex20 elements in (r, θ, z), the seam at θ = 0 shared.
fn tube(nr: usize, nt: usize, nz: usize, ri: f64, ro: f64, length: f64) -> Grid {
    let mut g = grid([nr, nt, nz], ElementKind::Hex20, |q| {
        let (r, t) = (ri + (ro - ri) * q[0], 2.0 * PI * q[1]);
        [r * t.cos(), r * t.sin(), length * q[2]]
    });
    // Merge the seam: nodes at θ = 2π take the labels of those at θ = 0.
    let seam: Vec<(u32, u32)> = g
        .param
        .iter()
        .filter(|(_, p)| near(p[1], 1.0))
        .map(|(l, p)| (*l, *g.param.iter().find(|(_, o)| near(o[1], 0.0) && near(o[0], p[0]) && near(o[2], p[2])).unwrap().0))
        .collect();
    let elements: Vec<(u32, fieldcraft_model::Element)> = g.model.elements().iter().map(|(l, e)| (*l, e.clone())).collect();
    for (label, e) in elements {
        let nodes: Vec<u32> = e.nodes.iter().map(|n| seam.iter().find(|(a, _)| a == n).map_or(*n, |(_, b)| *b)).collect();
        g.model.set_element(label, e.kind, nodes).unwrap();
    }
    for (old, _) in &seam {
        g.model.remove_node(*old).unwrap();
        g.param.remove(old);
    }
    g
}

/// Torsion of a tube by a moment on a reference point kinematically coupled to its end face, the
/// other end clamped: twist φ = T L / (G J), shear stress τ = T r / J, J = π (r_o⁴ − r_i⁴) / 2
/// (circular sections do not warp, so the clamp and the rigid end are exact).
#[test]
fn torsion_of_a_tube_through_a_reference_point() {
    let (ri, ro, length, young, nu, torque) = (0.6f64, 1.0f64, 4.0, 1000.0, 0.3, 5.0);
    let shear = young / (2.0 * (1.0 + nu));
    let j = PI / 2.0 * (ro.powi(4) - ri.powi(4));
    let mut g = tube(2, 16, 4, ri, ro, length);
    g.define(iso(young, nu), None, full(), Procedure::Static {});
    g.fixed("root", |q| near(q[2], 0.0));
    g.couple("Tip", [0.0, 0.0, length], CouplingKind::Kinematic, &ALL, |q| near(q[2], 1.0));
    g.load("twist", moment(rp("Tip"), [0.0, 0.0, torque]));
    let r = solve(&g);
    let phi = r.reference_points[0].displacement[5];
    let want = torque * length / (shear * j);
    let e_phi = (phi - want).abs() / want;
    // Shear stress on the outer surface at mid-length, τ_θz = −sin θ σxz + cos θ σyz.
    let outer = g.nodes_where(|q| near(q[0], 1.0) && near(q[2], 0.5));
    let tau = outer
        .iter()
        .map(|n| {
            let x = g.model.node(*n).unwrap().xyz;
            let t = x[1].atan2(x[0]);
            let s = r.stress[at(&r, *n)];
            -t.sin() * s[4] + t.cos() * s[5]
        })
        .sum::<f64>()
        / outer.len() as f64;
    let e_tau = (tau - torque * ro / j).abs() / (torque * ro / j);
    eprintln!("torsion: twist {phi:.6e} vs {want:.6e} ({:.3}%), τ {tau:.5} vs {:.5} ({:.3}%)", 100.0 * e_phi, torque * ro / j, 100.0 * e_tau);
    assert!(e_phi < 1e-3 && e_tau < 5e-3, "{e_phi} {e_tau}");
    // The rotation is about the axis only; the reaction moment at the root balances the torque
    // (reaction forces sum to zero).
    let p = &r.reference_points[0];
    assert!(p.displacement[..5].iter().all(|v| v.abs() < 1e-9 * phi), "{:?}", p.displacement);
    assert!(r.reaction_total.iter().all(|v| v.abs() < 1e-9));
}

#[test]
fn coupled_results_are_deterministic() {
    let build = || {
        let mut g = grid([4, 2, 2], ElementKind::Tet10, |q| [4.0 * q[0], q[1], q[2]]);
        g.define(iso(1000.0, 0.3), None, full(), Procedure::Static {});
        g.couple("Root", [0.0, 0.5, 0.5], CouplingKind::Kinematic, &ALL, |q| near(q[0], 0.0));
        g.bc("clamp", BoundaryCondition::Fixed { target: rp("Root") });
        g.couple("Tip", [4.5, 0.5, 0.5], CouplingKind::Distributing, &ALL, |q| near(q[0], 1.0));
        g.load("f", force(rp("Tip"), [0.0, 1.0, -2.0]));
        g.load("m", moment(rp("Tip"), [0.5, 0.0, 0.0]));
        g
    };
    let (a, b) = (solve(&build()), solve(&build()));
    assert_eq!(a, b);
    // Balance: forces and the reactions at the clamped point.
    for k in 0..3 {
        assert!((a.applied_total[k] + a.reaction_total[k]).abs() < 1e-9, "{:?} {:?}", a.applied_total, a.reaction_total);
    }
}

#[test]
fn coupling_errors_name_the_problem() {
    let base = || {
        let mut g = grid([2, 1, 1], ElementKind::Hex8, |q| [2.0 * q[0], q[1], q[2]]);
        g.define(iso(1000.0, 0.3), None, full(), Procedure::Static {});
        g.fixed("root", |q| near(q[0], 0.0));
        g
    };
    let err = |g: &Grid| solve_static(&g.model, Some("Step"), &Limits::default()).unwrap_err().to_string();

    // A moment on a point whose distributing coupling ties translations only.
    let mut g = base();
    g.couple("Tip", [2.0, 0.5, 0.5], CouplingKind::Distributing, &[1, 2, 3], |q| near(q[0], 1.0));
    g.load("m", moment(rp("Tip"), [1.0, 0.0, 0.0]));
    assert!(err(&g).contains("rotation about X is not coupled"), "{}", err(&g));

    let loose = |g: &mut Grid| g.model.add_reference_point("Loose", fieldcraft_model::ReferencePoint { xyz: [3.0, 0.0, 0.0], label: 0 }).unwrap();

    // A force on a reference point nothing couples.
    let mut g = base();
    loose(&mut g);
    g.couple("Tip", [2.0, 0.5, 0.5], CouplingKind::Kinematic, &ALL, |q| near(q[0], 1.0));
    g.load("f", force(rp("Loose"), [1.0, 0.0, 0.0]));
    assert!(err(&g).contains("no coupling ties"), "{}", err(&g));

    // Rotations through nodes on one line.
    let mut g = base();
    g.couple("Edge", [2.0, 0.0, 0.5], CouplingKind::Distributing, &ALL, |q| near(q[0], 1.0) && near(q[1], 0.0));
    g.load("f", force(rp("Edge"), [1.0, 0.0, 0.0]));
    assert!(err(&g).contains("one line"), "{}", err(&g));

    // Only rotations coupled kinematically.
    let mut g = base();
    g.couple("Tip", [2.0, 0.5, 0.5], CouplingKind::Kinematic, &[4, 5, 6], |q| near(q[0], 1.0));
    assert!(err(&g).contains("only rotations"), "{}", err(&g));

    // A rigid face whose nodes boundary conditions pull apart.
    let mut g = grid([2, 1, 1], ElementKind::Hex8, |q| [2.0 * q[0], q[1], q[2]]);
    g.define(iso(1000.0, 0.3), None, full(), Procedure::Static {});
    g.couple("Root", [0.0, 0.5, 0.5], CouplingKind::Kinematic, &ALL, |q| near(q[0], 0.0));
    for (name, y, v) in [("a", 0.0, 0.1), ("b", 1.0, 0.0)] {
        g.node_set(name, |q| near(q[0], 0.0) && near(q[1], y) && near(q[2], 0.0));
        g.bc(name, BoundaryCondition::Displacement { target: Target::Set(name.into()), values: [None, Some(v), None, None, None, None] });
    }
    assert!(err(&g).contains("contradicts"), "{}", err(&g));

    // A free reference point rotation: one coupled node gives the point's rotations no stiffness.
    let mut g = base();
    g.couple("One", [2.5, 0.5, 0.5], CouplingKind::Kinematic, &[1, 2, 3], |q| near(q[0], 1.0) && near(q[1], 1.0) && near(q[2], 1.0));
    g.load("f", force(rp("One"), [1.0, 0.0, 0.0]));
    assert!(err(&g).contains("reference point `One` in rotation"), "{}", err(&g));

    // A boundary condition on an uncoupled point holds nothing and says so.
    let mut g = base();
    loose(&mut g);
    g.bc("hold", BoundaryCondition::Fixed { target: rp("Loose") });
    g.face_set("end", |q| near(q[0], 1.0));
    g.load("p", Load::Pressure { set: "end".into(), magnitude: 1.0 });
    let r = solve(&g);
    assert!(r.warnings.iter().any(|w| w.contains("holds nothing")), "{:?}", r.warnings);

    // Hostile geometry: a reference point far away still solves or fails cleanly.
    for kind in [CouplingKind::Kinematic, CouplingKind::Distributing] {
        let mut g = base();
        g.couple("Far", [1e9, -1e9, 1e9], kind, &ALL, |q| near(q[0], 1.0));
        g.load("m", moment(rp("Far"), [1.0, 1.0, 1.0]));
        let _ = solve_static(&g.model, Some("Step"), &Limits::default());
    }
}

/// A distributing coupling on a large face whose reference point is held ties every coupled
/// node to every other (a dense block of (3N)² entries): refused up front, before the reduced
/// matrix is built, with the size it would need. Any other reduction that would outgrow the
/// budget stops while counting, before allocating.
#[test]
fn oversized_reductions_are_refused_before_allocating() {
    let mut g = grid([110, 110, 1], ElementKind::Hex8, |q| [110.0 * q[0], 110.0 * q[1], q[2]]);
    g.define(iso(1000.0, 0.3), None, full(), Procedure::Static {});
    g.couple("Base", [55.0, 55.0, 0.0], CouplingKind::Distributing, &ALL, |q| near(q[2], 0.0));
    g.bc("hold", BoundaryCondition::Fixed { target: rp("Base") });
    g.face_set("top", |q| near(q[2], 1.0));
    g.load("p", Load::Pressure { set: "top".into(), magnitude: 1.0 });
    let start = std::time::Instant::now();
    let raised = Limits { max_factor_entries: 2_000_000_000, ..Limits::default() };
    let e = solve_static(&g.model, Some("Step"), &raised).unwrap_err().to_string();
    eprintln!("refused in {:?}: {e}", start.elapsed());
    assert!(e.contains("couples 12321 nodes") && e.contains("GB") && e.contains("kinematic coupling"), "{e}");

    let mut g = grid([2, 1, 1], ElementKind::Hex8, |q| [2.0 * q[0], q[1], q[2]]);
    g.define(iso(1000.0, 0.3), None, full(), Procedure::Static {});
    g.couple("Root", [0.0, 0.5, 0.5], CouplingKind::Kinematic, &ALL, |q| near(q[0], 0.0));
    g.bc("clamp", BoundaryCondition::Fixed { target: rp("Root") });
    g.face_set("end", |q| near(q[0], 1.0));
    g.load("p", Load::Pressure { set: "end".into(), magnitude: 1.0 });
    let tight = Limits { max_factor_entries: 4, ..Limits::default() };
    let e = solve_static(&g.model, Some("Step"), &tight).unwrap_err().to_string();
    assert!(e.contains("reduced matrices need more than"), "{e}");
    assert!(solve_static(&g.model, Some("Step"), &Limits::default()).is_ok());
}

/// Contact belongs to explicit dynamics so far: the implicit solves refuse it by name.
#[test]
fn implicit_solves_refuse_contact_by_name() {
    let mut g = grid([1, 1, 1], ElementKind::Hex8, |x| x);
    g.define(iso(1000.0, 0.3), Some(1.0), full(), Procedure::Static {});
    g.face_set("top", |x| near(x[2], 1.0));
    g.model.add_contact("touch", fieldcraft_model::Contact::new(fieldcraft_model::ContactSurfaces::General { surface: Some("top".into()) })).unwrap();
    let e = solve_static(&g.model, Some("Step"), &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("contact (`touch`)") && e.contains("not supported"), "{e}");
}
