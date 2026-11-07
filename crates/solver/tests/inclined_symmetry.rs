//! Symmetry about an inclined plane (#242): n · u = 0 on the plane, and zero rotation about the
//! two axes in it at shell nodes, eliminated as linear constraints.
//!
//! The checks need no reference program: a finite element model turned by a rotation R solves
//! to R times the displacements of the unturned one (the elements are frame-indifferent), so an
//! inclined plane is checked against the axis-aligned plane it is the image of; and half of a
//! part whose mesh is mirror-symmetric about the plane, with the condition and half the load,
//! solves exactly as the whole does at the nodes they share.

mod common;

use std::collections::BTreeMap;

use common::{full, grid, iso, near, reduced};
use fieldcraft_model::{
    BoundaryCondition, Elastic, ElementFace, ElementKind, Formulation, Load, Material, Model, Procedure, Section, SectionKind, Set, Step, Symmetry,
    Target,
};
use fieldcraft_solver::{FrequencyOptions, Limits, StaticResult, solve_frequency, solve_static};

/// A rotation by `angle` about the unit axis `a` (Rodrigues).
fn rotation(a: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
    let l = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    let a = a.map(|c| c / l);
    let (s, c) = angle.sin_cos();
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let k = [[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]][i][j];
            r[i][j] = if i == j { c } else { 0.0 } + (1.0 - c) * a[i] * a[j] + s * k;
        }
    }
    r
}

fn apply(r: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [0, 1, 2].map(|i| (0..3).map(|j| r[i][j] * v[j]).sum())
}

fn norm(v: [f64; 3]) -> f64 {
    v.iter().map(|c| c * c).sum::<f64>().sqrt()
}

fn key(x: [f64; 3]) -> [i64; 3] {
    x.map(|c| (c * 1e6).round() as i64)
}

/// The displacements of a result by node position (rounded).
fn by_position(model: &Model, r: &StaticResult) -> BTreeMap<[i64; 3], ([f64; 3], f64)> {
    r.nodes.iter().enumerate().map(|(i, l)| (key(model.nodes()[l].xyz), (r.displacement[i], r.von_mises[i]))).collect()
}

/// A cantilever L × 2b × h along x (y from −b to b when `whole`, 0 to b otherwise) clamped at
/// x = 0, under a downward traction on its end; the half has the symmetry condition on y = 0.
/// Everything is turned by `turn` (nodes, the traction and the plane's normal).
fn cantilever(kind: ElementKind, formulation: Formulation, whole: bool, turn: &[[f64; 3]; 3]) -> (Model, StaticResult) {
    let (l, b, h) = (10.0, 1.0, 1.0);
    let n = if whole { [8, 4, 2] } else { [8, 2, 2] };
    let mut g = grid(n, kind, |p| apply(turn, [l * p[0], if whole { b * (2.0 * p[1] - 1.0) } else { b * p[1] }, h * p[2]]));
    g.define(iso(210e3, 0.3), None, formulation, Procedure::Static {});
    g.fixed("root", |p| near(p[0], 0.0));
    g.face_set("tip", |p| near(p[0], 1.0));
    g.load("push", Load::Traction { set: "tip".into(), vector: apply(turn, [0.0, 0.0, -10.0]) });
    if !whole {
        let normal = apply(turn, [0.0, 1.0, 0.0]);
        g.symmetry("mid", Symmetry::from_normal(normal).unwrap(), |p| near(p[1], 0.0));
    }
    let r = solve_static(&g.model, None, &Limits::default()).unwrap();
    (g.model, r)
}

/// The turned half cantilever with its inclined plane solves to the turned displacements of the
/// unturned half (and, for the mirror-symmetric hexahedral meshes, of the whole), at every node;
/// the von Mises stress is the same; the plane's reactions balance the load.
#[test]
fn a_turned_half_solid_matches_the_axis_aligned_half_and_the_whole() {
    let turn = rotation([1.0, 2.0, 3.0], 0.7);
    let identity = rotation([0.0, 0.0, 1.0], 0.0);
    for (kind, formulation, mirrored) in [
        (ElementKind::Hex8, full(), true),
        (ElementKind::Hex20, full(), true),
        (ElementKind::Hex20, reduced(), true),
        (ElementKind::Tet10, full(), false),
        (ElementKind::Wedge6, full(), false),
    ] {
        let (m_axis, axis) = cantilever(kind, formulation, false, &identity);
        let (m_turned, turned) = cantilever(kind, formulation, false, &turn);
        let size = axis.displacement.iter().fold(0.0f64, |m, u| m.max(norm(*u)));
        let stress = axis.von_mises.iter().fold(0.0f64, |m, s| m.max(*s));
        let reference = by_position(&m_axis, &axis);
        let mut worst = (0.0f64, 0.0f64);
        let unturn = |x: [f64; 3]| [0, 1, 2].map(|i| (0..3).map(|j| turn[j][i] * x[j]).sum::<f64>());
        for (i, l) in turned.nodes.iter().enumerate() {
            let (u, s) = reference[&key(unturn(m_turned.nodes()[l].xyz))];
            let want = apply(&turn, u);
            let got = turned.displacement[i];
            worst.0 = worst.0.max(norm([0, 1, 2].map(|k| got[k] - want[k])) / size);
            worst.1 = worst.1.max((turned.von_mises[i] - s).abs() / stress);
        }
        eprintln!(
            "{kind:?} {formulation:?}: tip deflection {size:.6e}, turned vs axis-aligned: displacement {:.1e}, von Mises {:.1e}",
            worst.0, worst.1
        );
        assert!(worst.0 < 1e-9 && worst.1 < 1e-8, "{kind:?}: {worst:?}");
        for r in [&axis, &turned] {
            let imbalance = (0..3).map(|c| (r.applied_total[c] + r.reaction_total[c]).abs()).fold(0.0, f64::max);
            assert!(imbalance < 1e-8 * norm(r.applied_total), "{kind:?}: {:?} vs {:?}", r.applied_total, r.reaction_total);
            assert!(!r.warnings.iter().any(|w| w.contains("reactions")), "{:?}", r.warnings);
        }
        if mirrored {
            let (m_whole, whole) = cantilever(kind, formulation, true, &identity);
            let all = by_position(&m_whole, &whole);
            let mut worst = 0.0f64;
            for (x, (u, _)) in &reference {
                let w = all[x].0;
                worst = worst.max(norm([0, 1, 2].map(|k| u[k] - w[k])) / size);
            }
            eprintln!("{kind:?}: half with the plane vs whole: displacement {worst:.1e}");
            assert!(worst < 1e-9, "{kind:?}: {worst}");
        }
    }
}

/// A shell model: one material and step; `n` × `n` cells of the quarter of a simply supported
/// square plate (side 1), two triangles per cell split along x = y, so the mesh is mirror-
/// symmetric about the diagonal plane x = y; `eighth` keeps the cells below the diagonal and
/// puts the inclined symmetry plane (normal (1, −1, 0)) on it.
fn plate(n: usize, eighth: bool, procedure: Procedure) -> (Model, BTreeMap<u32, [f64; 3]>) {
    let (a, t) = (1.0, 0.01);
    let mut m = Model::default();
    m.add_material("Mat", Material { elastic: Some(Elastic::Isotropic { young: 1e7, poisson: 0.3 }), density: Some(1.0), ..Material::default() })
        .unwrap();
    m.add_step(Step::new("Step", procedure), None).unwrap();
    let label = |i: usize, j: usize| (1 + i + (n + 1) * j) as u32;
    let mut xyz = BTreeMap::new();
    for j in 0..=n {
        for i in 0..=n {
            if eighth && j > i {
                continue;
            }
            let p = [i as f64 / n as f64 * a / 2.0, j as f64 / n as f64 * a / 2.0, 0.0];
            m.set_node(label(i, j), p).unwrap();
            xyz.insert(label(i, j), p);
        }
    }
    let mut elements = Vec::new();
    for j in 0..n {
        for i in 0..n {
            let (p, q, r, s) = (label(i, j), label(i + 1, j), label(i + 1, j + 1), label(i, j + 1));
            for (tri, below) in [([p, q, r], true), ([p, r, s], false)] {
                if eighth && (j > i || (i == j && !below)) {
                    continue;
                }
                let e = elements.len() as u32 + 1;
                m.set_element(e, ElementKind::Tri3, tri.to_vec()).unwrap();
                elements.push(e);
            }
        }
    }
    let faces = elements.iter().map(|e| ElementFace { element: *e, face: 1 }).collect();
    m.add_set("plate", Set { elements, ..Set::default() }).unwrap();
    m.add_section(
        "plate",
        Section {
            orientation: None,
            set: "plate".into(),
            material: "Mat".into(),
            kind: SectionKind::Shell { thickness: t, offset: 0.0, points: 5 },
            formulation: Formulation::default(),
        },
    )
    .unwrap();
    let on = |x: f64, c: f64| (x - c).abs() < 1e-9;
    let free = None;
    let hold = |m: &mut Model, name: &str, bc: BoundaryCondition, f: &dyn Fn([f64; 3]) -> bool| {
        let nodes = xyz.iter().filter(|(_, x)| f(**x)).map(|(l, _)| *l).collect();
        m.add_set(name, Set { nodes, ..Set::default() }).unwrap();
        m.add_boundary_condition(name, bc).unwrap();
        m.set_boundary_condition_active("Step", name, true).unwrap();
    };
    let held = |values| move |name: &str| BoundaryCondition::Displacement { target: Target::Set(name.into()), values };
    // Hard simple support on x = 0 and y = 0; symmetry on x = a/2 and y = a/2.
    let x0 = held([free, free, Some(0.0), Some(0.0), free, free]);
    let y0 = held([free, free, Some(0.0), free, Some(0.0), free]);
    let sym = |name: &str, plane| BoundaryCondition::Symmetry { target: Target::Set(name.into()), plane };
    if !eighth {
        hold(&mut m, "x0", x0("x0"), &|p| on(p[0], 0.0));
        hold(&mut m, "ys", sym("ys", Symmetry::Y), &|p| on(p[1], a / 2.0));
    }
    hold(&mut m, "y0", y0("y0"), &|p| on(p[1], 0.0));
    hold(&mut m, "xs", sym("xs", Symmetry::X), &|p| on(p[0], a / 2.0));
    if eighth {
        hold(&mut m, "diagonal", sym("diagonal", Symmetry::from_normal([1.0, -1.0, 0.0]).unwrap()), &|p| on(p[0], p[1]));
    }
    if matches!(m.steps()[0].procedure, Procedure::Static {}) {
        m.add_set("top", Set { faces, ..Set::default() }).unwrap();
        m.add_load("q", Load::Pressure { set: "top".into(), magnitude: 1.0 }).unwrap();
        m.set_load_active("Step", "q", true).unwrap();
    }
    (m, xyz)
}

/// An eighth of the simply supported plate, cut on the diagonal with the inclined plane, solves
/// as the quarter (deflections and rotations at every shared node, to round-off; the first
/// frequency, whose mode is symmetric about the diagonal), and both match Navier's series.
#[test]
fn an_eighth_plate_on_the_diagonal_plane_matches_the_quarter() {
    let n = 16;
    let (mq, quarter) = plate(n, false, Procedure::Static {});
    let (me, eighth) = plate(n, true, Procedure::Static {});
    let rq = solve_static(&mq, None, &Limits::default()).unwrap();
    let re = solve_static(&me, None, &Limits::default()).unwrap();
    let at = |r: &StaticResult, l: u32| r.nodes.binary_search(&l).unwrap();
    let size = rq.displacement.iter().fold(0.0f64, |m, u| m.max(norm(*u)));
    let turn = rq.rotation.iter().fold(0.0f64, |m, u| m.max(norm(*u)));
    let (mut du, mut dr) = (0.0f64, 0.0f64);
    for l in eighth.keys() {
        assert!(quarter.contains_key(l));
        let (a, b) = (at(&rq, *l), at(&re, *l));
        du = du.max(norm([0, 1, 2].map(|k| rq.displacement[a][k] - re.displacement[b][k])) / size);
        dr = dr.max(norm([0, 1, 2].map(|k| rq.rotation[a][k] - re.rotation[b][k])) / turn);
    }
    // The diagonal turns only about its normal: θx + θy = 0 and θz = 0 there.
    for (l, x) in &eighth {
        if (x[0] - x[1]).abs() < 1e-9 {
            let r = re.rotation[at(&re, *l)];
            assert!((r[0] + r[1]).abs() < 1e-12 * turn && r[2].abs() < 1e-12 * turn, "node {l}: {r:?}");
        }
    }
    let centre = *eighth.iter().find(|(_, x)| (x[0] - 0.5).abs() < 1e-9 && (x[1] - 0.5).abs() < 1e-9).unwrap().0;
    let w = -re.displacement[at(&re, centre)][2];
    let d = 1e7 * 0.01f64.powi(3) / (12.0 * (1.0 - 0.09));
    let navier = navier(1.0, 1.0, d);
    eprintln!("eighth vs quarter plate: displacement {du:.1e}, rotation {dr:.1e}; centre deflection {w:.6e} vs Navier {navier:.6e}");
    assert!(du < 1e-9 && dr < 1e-9, "{du} {dr}");
    assert!((w - navier).abs() < 0.01 * navier, "{w} vs {navier}");
    let imbalance = (0..3).map(|c| (re.applied_total[c] + re.reaction_total[c]).abs()).fold(0.0, f64::max);
    assert!(imbalance < 1e-9 * norm(re.applied_total), "{:?} vs {:?}", re.applied_total, re.reaction_total);
    assert!(!re.warnings.iter().any(|w| w.contains("reactions")), "{:?}", re.warnings);
    // An eighth carries an eighth of the plate's load (the quarter a quarter).
    assert!((2.0 * re.applied_total[2] - rq.applied_total[2]).abs() < 1e-12 * rq.applied_total[2].abs());

    let frequency = |eighth| {
        let (m, _) = plate(n, eighth, Procedure::Frequency { modes: 1, min_hz: None, max_hz: None });
        solve_frequency(&m, None, &FrequencyOptions::default(), &Limits::default()).unwrap().modes[0].frequency
    };
    let (fq, fe) = (frequency(false), frequency(true));
    eprintln!("first frequency: quarter {fq:.9}, eighth {fe:.9}");
    assert!((fq - fe).abs() < 1e-8 * fq, "{fq} vs {fe}");
}

/// Navier's series for the centre deflection of a simply supported square plate under uniform
/// pressure (S. Timoshenko, S. Woinowsky-Krieger, *Theory of Plates and Shells*, 2nd ed.,
/// McGraw-Hill 1959, §30).
fn navier(q: f64, a: f64, d: f64) -> f64 {
    let mut s = 0.0;
    for m in (1..400).step_by(2) {
        for n in (1..400).step_by(2) {
            let (mf, nf) = (f64::from(m), f64::from(n));
            let sign = if ((m + n) / 2 - 1) % 2 == 0 { 1.0 } else { -1.0 };
            s += sign / (mf * nf * (mf * mf + nf * nf).powi(2));
        }
    }
    16.0 * q * a.powi(4) / (std::f64::consts::PI.powi(6) * d) * s
}

/// A node that is also held, or listed in two planes, is no conflict; a prescribed motion off
/// the plane is.
#[test]
fn conflicts() {
    // The clamp's nodes on y = 0 are also on the plane: held there, the plane's rows are redundant.
    let (_, r) = cantilever(ElementKind::Hex8, full(), false, &rotation([0.0, 0.0, 1.0], 0.3));
    assert!(!r.warnings.iter().any(|w| w.contains("reactions")), "{:?}", r.warnings);
    let mut g = grid([4, 1, 1], ElementKind::Hex8, |p| [4.0 * p[0], p[1], p[2]]);
    g.define(iso(1e3, 0.3), None, full(), Procedure::Static {});
    g.fixed("root", |p| near(p[0], 0.0));
    g.symmetry("a", Symmetry::Normal([1.0, 1.0, 0.0]), |p| near(p[1], 0.0) && near(p[0], 1.0));
    g.symmetry("b", Symmetry::Normal([2.0, 2.0, 0.0]), |p| near(p[1], 0.0) && near(p[0], 1.0));
    g.face_set("tip", |p| near(p[0], 1.0));
    g.load("push", Load::Traction { set: "tip".into(), vector: [0.0, 0.0, -1.0] });
    let r = solve_static(&g.model, None, &Limits::default()).unwrap();
    let u = r.displacement[r.nodes.binary_search(&g.nodes_where(|p| near(p[1], 0.0) && near(p[0], 1.0))[0]).unwrap()];
    assert!((u[0] + u[1]).abs() < 1e-12 * norm(u), "{u:?}");
    // A prescribed motion the plane cannot follow is a contradiction, reported by name.
    let nodes = g.nodes_where(|p| near(p[1], 0.0) && near(p[0], 1.0));
    g.model.add_set("pulled", Set { nodes, ..Set::default() }).unwrap();
    g.bc("pull", BoundaryCondition::Displacement { target: Target::Set("pulled".into()), values: [Some(0.1), Some(0.1), None, None, None, None] });
    let e = solve_static(&g.model, None, &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("boundary condition `") && e.contains("contradicts"), "{e}");
}
