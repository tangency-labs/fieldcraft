//! Interpolated walls on a curved wall that moves: annular Couette flow between two coaxial
//! cylinders, the inner one sliding along its axis (theory: the closed-form logarithmic profile
//! u(r) = U ln(r/R₂) / ln(R₁/R₂), and the drag on the rod 2π μ U L / ln(R₂/R₁)). Both walls are
//! curved, the inner one convex and moving, so the case checks Ladd's moving-wall term and the
//! momentum exchange at the interpolated wall as well as the wall's position. Halfway bounce-back
//! puts both walls on staircases (about first order); interpolated walls put them on the
//! cylinders (about second order).
//!
//! Each test prints its errors
//! (`cargo test -p fieldcraft-lbm --test interpolated_couette -- --nocapture`).

mod common;

use common::*;
use fieldcraft_lbm::WallLinks;
use fieldcraft_mesh::surface::BoundarySurface;
use fieldcraft_mesh::tet::geometry::AnalyticSurface;
use fieldcraft_model::{FaceCondition, GeometryRef, SurfaceBoundary, SurfaceCondition, WallTreatment};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};
const HALFWAY: WallTreatment = WallTreatment::Halfway;
const INTERPOLATED: WallTreatment = WallTreatment::Interpolated;

/// The outer cylinder's radius (the inner one's is half of it), the kinematic viscosity and the
/// rod's speed.
const R2: f64 = 1.0;
const R1: f64 = 0.5;
const NU: f64 = 1.0 / 6.0;
const SPEED: f64 = 0.05;

/// A closed rod along x from `x0` to `x1` about the x axis, radius `r`, `m` segments. Faces: 0 the
/// side, 1 the end at `x0`, 2 the end at `x1`.
fn rod(x0: f64, x1: f64, r: f64, m: usize) -> BoundarySurface {
    let mut positions = Vec::new();
    for x in [x0, x1] {
        for k in 0..m {
            let phi = 2.0 * std::f64::consts::PI * k as f64 / m as f64;
            positions.push([x, r * phi.cos(), r * phi.sin()]);
        }
    }
    positions.push([x0, 0.0, 0.0]);
    positions.push([x1, 0.0, 0.0]);
    let id = |e: usize, s: usize| (e * m + s % m) as u32;
    let (start, end) = ((2 * m) as u32, (2 * m + 1) as u32);
    let mut tris = Vec::new();
    let mut tags = Vec::new();
    for s in 0..m {
        tris.push([id(0, s), id(1, s), id(1, s + 1)]);
        tris.push([id(0, s), id(1, s + 1), id(0, s + 1)]);
        tags.extend([0, 0]);
        tris.push([start, id(0, s), id(0, s + 1)]);
        tags.push(1);
        tris.push([end, id(1, s), id(1, s + 1)]);
        tags.push(2);
    }
    surface(positions, tris, tags, |p, tag| match tag {
        0 => [0.0, p[1], p[2]],
        1 => [-1.0, 0.0, 0.0],
        _ => [1.0, 0.0, 0.0],
    })
}

/// The closed-form axial velocity at radius `r` in the gap.
fn exact(r: f64) -> f64 {
    SPEED * (r / R2).ln() / (R1 / R2).ln()
}

/// The closed-form force of the fluid on a length `l` of the rod, along x (against its motion).
fn exact_drag(l: f64) -> f64 {
    -2.0 * std::f64::consts::PI * NU * SPEED * l / (R2 / R1).ln()
}

/// What one run measured.
struct Couette {
    /// The relative L2 error of the velocity over the fluid cells.
    l2: f64,
    /// The rod's drag over the closed form's, less 1.
    drag: f64,
    /// The rod's drag plus the outer wall's over the rod's (steady: the two balance).
    balance: f64,
    links: Option<WallLinks>,
}

/// Annular Couette flow with `k` cells per outer radius, one cell deep and periodic along x, run
/// from rest for five viscous times of the gap (τ = 1 on every grid).
fn couette(k: usize, w: WallTreatment) -> Couette {
    let dx = R2 / k as f64;
    let half = R2 + 2.0 * dx;
    let mut d = domain([0.0, -half, -half], [dx, half, half], dx);
    d.faces.x_min = PERIODIC;
    d.faces.x_max = PERIODIC;
    d.obstacles = vec!["Tube".into(), "Rod".into()];
    d.surfaces = vec![SurfaceBoundary {
        face: GeometryRef::face("Rod", "side").unwrap(),
        condition: SurfaceCondition::MovingWall { velocity: [SPEED, 0.0, 0.0] },
    }];
    let axis = |r: f64| Some(AnalyticSurface::Cylinder { origin: [0.0; 3], axis: [1.0, 0.0, 0.0], radius: r });
    let mut tube = obstacle("Tube", tube(-1.0, dx + 1.0, [0.0, 0.0], R2, 3.0 * half, 512), &["inner", "outer", "start", "end"]);
    tube.exact = vec![axis(R2), None, None, None];
    let mut inner = obstacle("Rod", rod(-1.0, dx + 1.0, R1, 512), &["side", "start", "end"]);
    inner.exact = vec![axis(R1), None, None];
    let gap = R2 - R1;
    let m = model(1.0, NU, d, with_walls(step(5.0 * gap * gap / NU, dx * dx, None), w));
    let (setup, mut flow) = prepare(&m, &[tube, inner]);
    run_all(&setup, &mut flow);
    let grid = *flow.grid();
    let (mut e2, mut u2) = (0.0, 0.0);
    for c in 0..grid.len() {
        let Some((_, u)) = flow.moments(c) else { continue };
        let [_, y, z] = grid.centre(c);
        let reference = exact((y * y + z * z).sqrt());
        let u = u[0] * setup.velocity_scale();
        e2 += (u - reference) * (u - reference);
        u2 += reference * reference;
    }
    let force = |name: &str| -> f64 {
        flow.groups().iter().zip(flow.forces()).filter(|(g, _)| g.name == name).map(|(_, f)| f[0]).sum::<f64>() * setup.force_scale()
    };
    let (on_rod, on_tube) = (force("side"), force("inner"));
    Couette { l2: (e2 / u2).sqrt(), drag: on_rod / exact_drag(dx) - 1.0, balance: (on_rod + on_tube) / on_rod, links: flow.wall_links() }
}

/// The orders of the velocity's L2 error and of the drag error, each wall treatment, over grids:
/// [halfway L2, interpolated L2, halfway drag, interpolated drag].
fn study(cells: &[usize]) -> [f64; 4] {
    let mut runs: [Vec<f64>; 4] = Default::default();
    for &k in cells {
        let h = couette(k, HALFWAY);
        let i = couette(k, INTERPOLATED);
        let links = i.links.unwrap_or_default();
        println!(
            "couette, {k} cells per radius: L2 halfway {:.3e}, interpolated {:.3e}; drag halfway {:+.3e}, interpolated {:+.3e}; balance {:.1e}, {:.1e}; links {links:?}",
            h.l2, i.l2, h.drag, i.drag, h.balance, i.balance
        );
        assert!(h.links.is_none());
        // Every link into either cylinder finds its surface and has a fluid cell behind it.
        assert!(links.interpolated > 0 && links.halfway() == 0, "{links:?}");
        // At steady state the outer wall takes back what the rod gives the fluid.
        assert!(h.balance.abs() < 1e-6 && i.balance.abs() < 1e-6, "{} {}", h.balance, i.balance);
        // The interpolated walls are closer on every grid (measured: 7 to 31 times), and so is
        // the rod's drag: within 0.5 % (measured 0.11 % at most), where the staircase's moves
        // between 0.5 and 14 % as the grid cuts the rod differently.
        assert!(i.l2 < 0.5 * h.l2, "{k}: {} {}", i.l2, h.l2);
        assert!(i.drag.abs() < 5e-3 && i.drag.abs() < 0.25 * h.drag.abs(), "{k}: interpolated {}, halfway {}", i.drag, h.drag);
        for (run, v) in runs.iter_mut().zip([h.l2, i.l2, h.drag, i.drag]) {
            run.push(v);
        }
    }
    let orders = runs.each_ref().map(|r| fitted_order(cells, r));
    println!(
        "couette, {cells:?} cells per radius: L2 halfway {} (order {:.2}), interpolated {} (order {:.2}); drag halfway {} (order {:.2}), interpolated {} (order {:.2})",
        sci(&runs[0]),
        orders[0],
        sci(&runs[1]),
        orders[1],
        sci(&runs[2]),
        orders[2],
        sci(&runs[3]),
        orders[3]
    );
    orders
}

/// The trunk gate's representative: the sliding rod at 8 cells per outer radius, halfway against
/// interpolated (measured: L2 5.3e-2 against 3.9e-3, drag +1.5 % against −0.11 %).
#[test]
fn a_sliding_rod_on_eight_cells_per_radius_is_closer_with_interpolated_walls() {
    let (h, i) = (couette(8, HALFWAY), couette(8, INTERPOLATED));
    println!(
        "couette, 8 cells per radius: L2 halfway {:.3e}, interpolated {:.3e}; drag halfway {:+.3e}, interpolated {:+.3e}",
        h.l2, i.l2, h.drag, i.drag
    );
    let links = i.links.unwrap_or_default();
    assert!(links.interpolated > 0 && links.halfway() == 0, "{links:?}");
    assert!(h.balance.abs() < 1e-6 && i.balance.abs() < 1e-6, "{} {}", h.balance, i.balance);
    assert!(i.l2 < 0.1 * h.l2, "{} {}", i.l2, h.l2);
    assert!(i.drag.abs() < 5e-3 && i.drag.abs() < 0.25 * h.drag.abs(), "interpolated {}, halfway {}", i.drag, h.drag);
}

/// The sliding rod on 6 to 16 cells per outer radius (3 to 8 across the rod's radius).
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn annular_couette_with_a_sliding_rod_converges_at_second_order_with_interpolated_walls() {
    let [halfway, interpolated, ..] = study(&[6, 8, 12, 16]);
    assert!(interpolated > 1.6, "interpolated order {interpolated}");
    assert!(halfway < 1.3, "halfway order {halfway}");
}

/// The same study on seven grids up to 48 cells per outer radius.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn annular_couette_orders_on_seven_grids() {
    let [halfway, interpolated, ..] = study(&[6, 8, 12, 16, 24, 32, 48]);
    assert!(interpolated > 1.8, "interpolated order {interpolated}");
    assert!(halfway < 1.3, "halfway order {halfway}");
}
