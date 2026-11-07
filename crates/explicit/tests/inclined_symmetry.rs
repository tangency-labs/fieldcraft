//! Symmetry about inclined planes in explicit dynamics (#242): the velocity along the plane's
//! normal is removed at its nodes every increment. A bar turned by a rotation R moves as R times
//! the unturned bar on the axis-aligned planes; and a half (or quarter) of a bar whose mesh is
//! mirror-symmetric moves as the whole at the nodes they share. One plane, and two planes meeting
//! at an edge.

use fieldcraft_explicit::{ExplicitResult, Limits, solve_explicit};
use fieldcraft_model::*;

/// A rotation by `angle` about the axis `a` (Rodrigues).
fn rotation(a: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
    let l = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    let a = a.map(|c| c / l);
    let (s, c) = angle.sin_cos();
    let k = [[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]];
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = if i == j { c } else { 0.0 } + (1.0 - c) * a[i] * a[j] + s * k[i][j];
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

/// A 6 × 1 × 1 bar of one-point hex8 cells along x, clamped at x = 0, under the body force
/// `force` for 0.2 time units. `whole[a]` spans −1..1 across axis a (y or z) with two cells,
/// otherwise 0..1 with one cell and a symmetry plane at 0. Everything is turned by `turn`.
fn bar(whole: [bool; 2], force: [f64; 3], turn: &[[f64; 3]; 3]) -> (Model, ExplicitResult) {
    let mut m = Model::default();
    let cells = [6usize, if whole[0] { 2 } else { 1 }, if whole[1] { 2 } else { 1 }];
    let start = [0.0, if whole[0] { -1.0 } else { 0.0 }, if whole[1] { -1.0 } else { 0.0 }];
    let label = |i: usize, j: usize, k: usize| (1 + i + (cells[0] + 1) * (j + (cells[1] + 1) * k)) as u32;
    for k in 0..=cells[2] {
        for j in 0..=cells[1] {
            for i in 0..=cells[0] {
                m.set_node(label(i, j, k), apply(turn, [i as f64, start[1] + j as f64, start[2] + k as f64])).unwrap();
            }
        }
    }
    let corners = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
    let mut next = 1;
    for k in 0..cells[2] {
        for j in 0..cells[1] {
            for i in 0..cells[0] {
                m.set_element(next, ElementKind::Hex8, corners.iter().map(|c| label(i + c[0], j + c[1], k + c[2])).collect()).unwrap();
                next += 1;
            }
        }
    }
    let all: Vec<u32> = m.elements().keys().copied().collect();
    m.add_set("all", Set { elements: all, ..Set::default() }).unwrap();
    let nodes = |f: &dyn Fn(usize, usize, usize) -> bool| -> Vec<u32> {
        let mut out = Vec::new();
        for k in 0..=cells[2] {
            for j in 0..=cells[1] {
                for i in 0..=cells[0] {
                    if f(i, j, k) {
                        out.push(label(i, j, k));
                    }
                }
            }
        }
        out.sort_unstable();
        out
    };
    m.add_set("root", Set { nodes: nodes(&|i, _, _| i == 0), ..Set::default() }).unwrap();
    m.add_set("tip", Set { nodes: nodes(&|i, _, _| i == cells[0]), ..Set::default() }).unwrap();
    m.add_material("Mat", Material { elastic: Some(Elastic::Isotropic { young: 1.0e4, poisson: 0.3 }), density: Some(1.0), ..Material::default() })
        .unwrap();
    let formulation = Formulation { integration: Integration::Reduced, ..Formulation::default() };
    m.add_section("bar", Section { orientation: None, set: "all".into(), material: "Mat".into(), kind: SectionKind::Solid {}, formulation }).unwrap();
    let settings = ExplicitDynamic::new(0.2);
    let mut step = Step::new("Pull", Procedure::ExplicitDynamic(settings));
    m.add_boundary_condition("clamp", BoundaryCondition::Fixed { target: Target::Set("root".into()) }).unwrap();
    step.boundary_conditions.insert("clamp".into());
    for (a, name) in [(0usize, "y0"), (1, "z0")] {
        if whole[a] {
            continue;
        }
        m.add_set(name, Set { nodes: nodes(&|_, j, k| [j, k][a] == 0), ..Set::default() }).unwrap();
        let mut normal = [0.0; 3];
        normal[a + 1] = 1.0;
        let plane = Symmetry::from_normal(apply(turn, normal)).unwrap();
        m.add_boundary_condition(name, BoundaryCondition::Symmetry { target: Target::Set(name.into()), plane }).unwrap();
        step.boundary_conditions.insert(name.into());
    }
    m.add_load("pull", Load::BodyForce { set: "all".into(), vector: apply(turn, force) }).unwrap();
    step.loads.insert("pull".into());
    m.add_step(step, None).unwrap();
    let r = solve_explicit(&m, None, &Limits::default()).unwrap();
    (m, r)
}

fn key(x: [f64; 3]) -> [i64; 3] {
    x.map(|c| (c * 1e6).round() as i64)
}

/// The largest difference between `b`'s last-frame displacements and `turn` times `a`'s at the
/// same (unturned) places, over the nodes of `b`, relative to `a`'s largest.
fn difference(a: &(Model, ExplicitResult), b: &(Model, ExplicitResult), turn: &[[f64; 3]; 3]) -> f64 {
    let ua = &a.1.frames.last().unwrap().displacement;
    let at: std::collections::BTreeMap<[i64; 3], [f64; 3]> = a.1.nodes.iter().zip(ua).map(|(l, u)| (key(a.0.node(*l).unwrap().xyz), *u)).collect();
    let size = ua.iter().fold(0.0f64, |m, u| m.max(norm(*u)));
    let unturn = |x: [f64; 3]| [0, 1, 2].map(|i| (0..3).map(|j| turn[j][i] * x[j]).sum::<f64>());
    let ub = &b.1.frames.last().unwrap().displacement;
    b.1.nodes
        .iter()
        .zip(ub)
        .map(|(l, u)| {
            let want = apply(turn, at[&key(unturn(b.0.node(*l).unwrap().xyz))]);
            norm([0, 1, 2].map(|k| u[k] - want[k])) / size
        })
        .fold(0.0, f64::max)
}

#[test]
fn turned_planes_move_as_the_axis_aligned_ones_and_the_whole() {
    let identity = rotation([0.0, 0.0, 1.0], 0.0);
    let turn = rotation([1.0, 2.0, 3.0], 0.7);
    // Bending about one plane; an axial pull and a bend about both, two planes meeting at an edge.
    for (whole, force) in [([false, true], [0.0, 0.0, -1.0]), ([false, false], [4.0, 0.0, 0.0])] {
        let axis = bar(whole, force, &identity);
        let turned = bar(whole, force, &turn);
        let full = bar([true, true], force, &identity);
        let full_turned = bar([true, true], force, &turn);
        // The part moves as the whole, turned or not, to round-off: the planes hold exactly.
        let (mirrored, mirrored_turned) = (difference(&full, &axis, &identity), difference(&full_turned, &turned, &identity));
        // Turning the model changes the stable increment in its last digits, which the central
        // difference scheme carries over its increments: the turned part may differ from the
        // unturned by as much as the turned whole does (its increments may end a step apart), no more.
        let (rotated, own) = (difference(&axis, &turned, &turn), difference(&full, &full_turned, &turn));
        let tip = axis.1.frames.last().unwrap().displacement.iter().fold(0.0f64, |m, u| m.max(norm(*u)));
        eprintln!(
            "planes {whole:?}: tip {tip:.4e}; part vs whole {mirrored:.1e} (turned {mirrored_turned:.1e}); turned vs axis-aligned {rotated:.1e} (the whole: {own:.1e})"
        );
        assert!(tip > 1e-4, "{tip}");
        assert!(mirrored < 1e-12 && mirrored_turned < 1e-12, "{mirrored} {mirrored_turned}");
        assert!(rotated < 1.01 * own + 1e-12 && own < 1e-3, "{rotated} {own}");
        assert!(turned.1.completed);
    }
}

/// A prescribed velocity that leaves no free component along the normal must move along the
/// plane; one that does not is refused.
#[test]
fn a_prescribed_motion_off_the_plane_is_refused() {
    let turn = rotation([0.0, 0.0, 1.0], 0.5);
    let (mut m, _) = bar([false, true], [0.0, 0.0, -1.0], &turn);
    let nodes = m.set("y0").unwrap().nodes.iter().copied().filter(|l| m.set("tip").unwrap().nodes.contains(l)).collect();
    m.add_set("pushed", Set { nodes, ..Set::default() }).unwrap();
    let push = BoundaryCondition::Displacement { target: Target::Set("pushed".into()), values: [Some(0.01), Some(0.01), None, None, None, None] };
    m.add_boundary_condition("push", push).unwrap();
    m.set_boundary_condition_active("Pull", "push", true).unwrap();
    let e = solve_explicit(&m, None, &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("moves off a symmetry plane"), "{e}");
}
