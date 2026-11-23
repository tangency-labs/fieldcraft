//! The inlet corner: with an on-node wall law, a cell in a velocity inlet's layer that touches two walls (its in-plane
//! neighbours: a solid cell, or a box face that is a wall) is a wall cell (plain bounce-back, no inlet condition, no
//! wall law). Without the rule the density at such a cell grew until the run died (a slotted channel on 16,000 cells,
//! wall law on, around step 480). The inlet's velocity is raised so its volume flow stays. Other wall treatments are unchanged.
mod common;

use common::*;
use fieldcraft_lbm::backend::Backend;
use fieldcraft_lbm::{Limits, Setup};
use fieldcraft_model::{Collision, FaceCondition, HrrSigma, Outflow, Procedure, Profile, Turbulence, WallModel};

const DX: f64 = 0.05;

fn procedure(wall_model: Option<WallModel>, steps: u64, time_step: f64) -> Procedure {
    let sa = Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model };
    Procedure::Fluid {
        domain: "Flow".into(),
        end_time: None,
        steps: Some(steps),
        time_step: Some(time_step),
        frames: 2,
        magic: None,
        sigma: Some(HrrSigma::Value(0.95)),
        collision: Some(Collision::Hrr),
        turbulence: Some(sa),
        steady: None,
        walls: fieldcraft_model::WallTreatment::Halfway,
    }
}

/// A channel `nx` x `ny` x `nz` cells with a parabolic inlet at x min, a pressure outlet at x max, and the y and z faces
/// as given; Smagorinsky turbulence with `wall_model`, with the given obstacles.
fn channel(
    dims: (usize, usize, usize),
    y_periodic: bool,
    z_max_slip: bool,
    wall_model: Option<WallModel>,
    names: Vec<String>,
    obstacles: Vec<fieldcraft_lbm::Obstacle>,
) -> (Setup, Vec<fieldcraft_lbm::Obstacle>) {
    let m = channel_model(dims, y_periodic, z_max_slip, wall_model, names);
    let setup = Setup::new(&m, None, &obstacles, &Limits::default()).unwrap();
    (setup, obstacles)
}

fn channel_model(
    (nx, ny, nz): (usize, usize, usize),
    y_periodic: bool,
    z_max_slip: bool,
    wall_model: Option<WallModel>,
    names: Vec<String>,
) -> fieldcraft_model::Model {
    let mut d = domain([0.0; 3], [nx as f64 * DX, ny as f64 * DX, nz as f64 * DX], DX);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.1, 0.0, 0.0], profile: Profile::Parabolic };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
    if y_periodic {
        d.faces.y_min = FaceCondition::Periodic {};
        d.faces.y_max = FaceCondition::Periodic {};
    }
    if z_max_slip {
        d.faces.z_max = FaceCondition::Slip {};
    }
    d.obstacles = names;
    model(1.0, 1.0e-3, d, procedure(wall_model, 10, 0.01))
}

/// All four sides walls: the four edge cells of the inlet face are wall cells, reported by the flow; the cells next to
/// them (one wall each) stay fluid and keep the inlet's rule.
#[test]
fn an_inlet_cell_touching_two_walls_is_a_wall_cell() {
    let (setup, obstacles) = channel((16, 8, 8), false, false, Some(WallModel::SpaldingBounce), vec![], vec![]);
    let flow = fieldcraft_lbm::Flow::new(&setup, &obstacles).unwrap();
    let report = flow.inlet_corner_report();
    assert_eq!((report.cells, report.inlet_cells), (4, 64), "{report:?}");
    // Only the cells on no wall's row deliver: the inner 6 x 6 of the 8 x 8 layer. The rescale is the whole parabolic
    // flux over that of the inner cells, per axis squared.
    let w = |i: usize| {
        let s = (i as f64 + 0.5) / 8.0;
        4.0 * s * (1.0 - s)
    };
    let (all, inner): (f64, f64) = ((0..8).map(w).sum(), (1..7).map(w).sum());
    let expected = (all / inner).powi(2);
    assert!((report.flux_rescale / expected - 1.0).abs() < 1e-12, "{report:?}, expected {expected}");
    let g = flow.grid();
    for (j, k) in [(0, 0), (7, 0), (0, 7), (7, 7)] {
        assert!(!Backend::is_fluid(&flow, g.index([0, j, k])), "corner ({j},{k}) is a wall cell");
    }
    for (j, k) in [(1, 0), (0, 1), (3, 0), (7, 3), (3, 3)] {
        assert!(Backend::is_fluid(&flow, g.index([0, j, k])), "({j},{k}) keeps the inlet rule");
    }
}

/// A solved run reports the count and the factor the inlet was raised by, and a run without the rule reports nothing.
#[test]
fn a_run_prints_the_reclassified_cells_and_the_raise() {
    let m = channel_model((16, 8, 8), false, false, Some(WallModel::SpaldingBounce), vec![]);
    let r = fieldcraft_lbm::solve_fluid(&m, None, &[], &Limits::default()).unwrap();
    let line = r.warnings.iter().find(|w| w.starts_with("inlet corner rule:")).expect("the run reports the corner rule");
    assert!(line.contains("4 of 64 inlet-layer cells") && line.contains("raised by up to 1.") && line.contains("by flux"), "{line}");
    let m = channel_model((16, 8, 8), false, false, None, vec![]);
    let r = fieldcraft_lbm::solve_fluid(&m, None, &[], &Limits::default()).unwrap();
    assert!(r.warnings.iter().all(|w| !w.starts_with("inlet corner rule:")), "{:?}", r.warnings);
}

/// Walls on one side only (periodic span, z max slip): the inlet's cells on the z min wall's row touch one wall, none
/// touches two, so none is reclassified; the row takes the wall's velocity (at rest) and delivers nothing, so the inlet's
/// velocity is raised and the requested flow is delivered to round-off. (Before, the inlet delivered less.)
#[test]
fn an_inlet_meeting_a_wall_without_a_corner_cell_delivers_its_flow() {
    let (nx, ny, nz) = (16, 8, 8);
    let (setup, obstacles) = channel((nx, ny, nz), true, true, Some(WallModel::SpaldingBounce), vec![], vec![]);
    let flow = fieldcraft_lbm::Flow::new(&setup, &obstacles).unwrap();
    let report = flow.inlet_corner_report();
    assert_eq!(report.cells, 0, "{report:?}");
    assert!(report.flux_rescale > 1.0, "the wall row's flux is made up: {report:?}");
    let g = flow.grid();
    assert!(Backend::is_fluid(&flow, g.index([0, 0, 0])));
    assert!(Backend::is_fluid(&flow, g.index([0, 3, 0])));
    let u0 = setup.to_lattice([0.1, 0.0, 0.0])[0];
    let profile = |k: usize| {
        let s = (k as f64 + 0.5) / nz as f64;
        4.0 * s * (1.0 - s)
    };
    let requested: f64 = (0..ny).flat_map(|j| (0..nz).map(move |k| (j, k))).map(|(_, k)| u0 * profile(k)).sum();
    let delivered: f64 = flow
        .open_cells()
        .iter()
        .filter(|o| o.axis == 0 && o.inward > 0.0)
        .map(|o| match o.target {
            fieldcraft_lbm::Target::Velocity(v) => v[0],
            _ => 0.0,
        })
        .sum();
    assert!((delivered / requested - 1.0).abs() < 1e-12, "delivered {delivered}, requested {requested}");
}

/// Without a wall law nothing changes: the corner cells are inlet cells as before.
#[test]
fn without_a_wall_law_the_corner_cells_are_unchanged() {
    let (setup, obstacles) = channel((16, 8, 8), false, false, None, vec![], vec![]);
    let flow = fieldcraft_lbm::Flow::new(&setup, &obstacles).unwrap();
    assert_eq!(flow.inlet_corner_report(), Default::default());
    assert!(Backend::is_fluid(&flow, flow.grid().index([0, 0, 0])));
}

/// The section's shape: a triangle of solid cells in the inlet's corner (j + k < 6), built as one cuboid per row, leaves a
/// diagonal staircase rim. Every cell on the rim touches two walls (a solid cell below and one beside it), so the 7 cells
/// on j + k = 6 are wall cells, with the 3 free box corners of the 12 x 12 section; the other cells keep the inlet.
/// A 12 x 12 inlet whose lower rim is a staircase of six steps (21 solid cells in the layer).
fn staircase() -> (Setup, Vec<fieldcraft_lbm::Obstacle>, usize) {
    let section = 12;
    let mut names = Vec::new();
    let mut obstacles = Vec::new();
    for j in 0..6usize {
        let name = format!("Step{j}");
        let (y0, y1) = (j as f64 * DX, (j + 1) as f64 * DX);
        obstacles.push(obstacle(&name, cuboid([0.0, y0, 0.0], [0.5, y1, (6 - j) as f64 * DX]), &CUBOID_FACES));
        names.push(name);
    }
    let (setup, obstacles) = channel((16, section, section), false, false, Some(WallModel::SpaldingBounce), names, obstacles);
    (setup, obstacles, section)
}

#[test]
fn a_staircased_inlet_rim_is_wall_cells_with_the_flow_rescaled() {
    let (setup, obstacles, section) = staircase();
    let flow = fieldcraft_lbm::Flow::new(&setup, &obstacles).unwrap();
    let report = flow.inlet_corner_report();
    let fluid_in_layer = section * section - 21;
    assert_eq!((report.cells, report.inlet_cells), (10, fluid_in_layer), "{report:?}");
    assert!(report.flux_rescale > 1.0, "the inlet's volume flow is kept: {report:?}");
    let g = flow.grid();
    for j in 0..=6usize {
        assert!(!Backend::is_fluid(&flow, g.index([0, j, 6 - j])), "rim cell ({j},{}) is a wall cell", 6 - j);
    }
    for (j, k) in [(section - 1, 0), (0, section - 1), (section - 1, section - 1)] {
        assert!(!Backend::is_fluid(&flow, g.index([0, j, k])), "box corner ({j},{k})");
    }
    for (j, k) in [(1, 6), (2, 5), (3, 4), (4, 4), (6, 6), (6, 1), (8, 8)] {
        assert!(Backend::is_fluid(&flow, g.index([0, j, k])), "({j},{k}) keeps the inlet rule");
    }
}

/// The slotted channel the corner failed in (2 x 1 x 1 m on 40 x 20 x 20 cells, a plate with a 4-cell slot, hybrid collision sigma 0.95, tau - 1/2 = 2.4e-4, inlet mean 0.035 lattice, walls on y min/max and
/// z min, slip or a wall on z max): without the rule the density at the inlet's wall corners grows to 1 +- 1.1 and the
/// run goes unstable at step 481 with the Spalding law (474 on four walls); with it the run holds 600 steps.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_slotted_channel_with_wall_laws_holds_with_the_corner_rule() {
    for (law, four_walls) in [(WallModel::SpaldingBounce, false), (WallModel::SpaldingBounce, true)] {
        let dx = 0.05;
        let dt = 0.035 * dx / 0.1;
        let nu = 2.4e-4 / 3.0 * dx * dx / dt;
        let mut d = domain([0.0; 3], [2.0, 1.0, 1.0], dx);
        d.faces.x_min = FaceCondition::Velocity { velocity: [0.1, 0.02, 0.0], profile: Profile::Parabolic };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
        if !four_walls {
            d.faces.z_max = FaceCondition::Slip {};
        }
        d.acceleration = [0.0, 0.0, -0.01];
        d.obstacles = vec!["PlateLow".into(), "PlateHigh".into()];
        let obstacles = vec![
            obstacle("PlateLow", cuboid([0.6, 0.0, 0.0], [0.8, 1.0, 0.4]), &CUBOID_FACES),
            obstacle("PlateHigh", cuboid([0.6, 0.0, 0.6], [0.8, 1.0, 1.0]), &CUBOID_FACES),
        ];
        let m = model(1.0, nu, d, procedure(Some(law), 600, dt));
        let setup = Setup::new(&m, None, &obstacles, &Limits::default()).unwrap();
        let mut flow = fieldcraft_lbm::Flow::new(&setup, &obstacles).unwrap();
        Backend::advance(&mut flow, 600).unwrap_or_else(|e| panic!("{law:?}, four walls {four_walls}: {e}"));
        let f = Backend::fields(&mut flow).unwrap();
        let ok = (0..f.density.len()).filter(|&c| Backend::is_fluid(&flow, c)).all(|c| f.density[c].is_finite() && (f.density[c] - 1.0).abs() < 0.3);
        // What this guards: the density mode the rule removes (1 +- 1.1 before the run stopped). A healthy run at this
        // lattice speed varies the density by about Ma^2, a few thousandths, so 0.3 only has to separate the two.
        assert!(ok, "{law:?}, four walls {four_walls}: a density left 1 +- 0.3");
    }
}

/// The flow the inlet layer delivers (the x component of the velocity every inlet cell is driven with, lattice units).
fn delivered(flow: &fieldcraft_lbm::Flow) -> f64 {
    flow.open_cells()
        .iter()
        .filter(|o| o.axis == 0 && o.inward > 0.0)
        .map(|o| match o.target {
            fieldcraft_lbm::Target::Velocity(v) => v[0],
            _ => 0.0,
        })
        .sum()
}

/// What bounds the inlet's y axis: walls, a periodic span, or a wall at y min and a symmetry plane at y max (the half of a
/// channel, whose profile peaks at the plane).
#[derive(Clone, Copy, PartialEq)]
enum Span {
    Walls,
    Periodic,
    SymmetryAtMax,
}

/// The flow the setup asked of an inlet on a `ny` x `nz` section with a parabolic profile across the axes that are not
/// periodic, over the cells that are not `solid` in the layer (lattice units, the profile counted at cell centres).
fn requested(setup: &Setup, (ny, nz): (usize, usize), y: Span, solid: impl Fn(usize, usize) -> bool) -> f64 {
    let u0 = setup.to_lattice([0.1, 0.0, 0.0])[0];
    let profile = |i: usize, n: usize| {
        let s = (i as f64 + 0.5) / n as f64;
        4.0 * s * (1.0 - s)
    };
    let profile_y = |j: usize| match y {
        Span::Walls => profile(j, ny),
        Span::Periodic => 1.0,
        Span::SymmetryAtMax => profile(j, 2 * ny) * 1.0,
    };
    (0..ny).flat_map(|j| (0..nz).map(move |k| (j, k))).filter(|&(j, k)| !solid(j, k)).map(|(j, k)| u0 * profile_y(j) * profile(k, nz)).sum()
}

/// The velocity the inlet delivers is the one it was asked for, to round-off: with the edge cells walls and the row
/// beside them on the walls' plane (at rest), the cells that keep the inlet carry the whole requested flow, parabolic
/// profile included. Four walls (corner cells), a duct with a periodic span (wall rows and no corner cell), and a
/// staircased rim (corner cells made by corner cells).
#[test]
fn the_inlet_delivers_the_requested_flow() {
    let n = 8;
    for (what, y_periodic) in [("four walls", false), ("duct with a periodic span", true)] {
        let (setup, obstacles) = channel((16, n, n), y_periodic, false, Some(WallModel::SpaldingBounce), vec![], vec![]);
        let flow = fieldcraft_lbm::Flow::new(&setup, &obstacles).unwrap();
        let (d, r) = (delivered(&flow), requested(&setup, (n, n), if y_periodic { Span::Periodic } else { Span::Walls }, |_, _| false));
        assert!((d / r - 1.0).abs() < 1e-12, "{what}: delivered {d}, requested {r}");
    }
    // A symmetry plane at y max: the profile is the half of a channel (s / 2), so the plane's end is its peak.
    let mut d = domain([0.0; 3], [16.0 * DX, n as f64 * DX, n as f64 * DX], DX);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.1, 0.0, 0.0], profile: Profile::Parabolic };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
    d.faces.y_max = FaceCondition::Symmetry {};
    let m = model(1.0, 1.0e-3, d, procedure(Some(WallModel::SpaldingBounce), 10, 0.01));
    let setup = Setup::new(&m, None, &[], &Limits::default()).unwrap();
    let flow = fieldcraft_lbm::Flow::new(&setup, &[]).unwrap();
    let (d, r) = (delivered(&flow), requested(&setup, (n, n), Span::SymmetryAtMax, |_, _| false));
    assert!((d / r - 1.0).abs() < 1e-12, "symmetry plane: delivered {d}, requested {r}");
    let (setup, obstacles, section) = staircase();
    let flow = fieldcraft_lbm::Flow::new(&setup, &obstacles).unwrap();
    let (d, r) = (delivered(&flow), requested(&setup, (section, section), Span::Walls, |j, k| j < 6 && k < 6 - j));
    assert!((d / r - 1.0).abs() < 1e-12, "staircase: delivered {d}, requested {r}");
}

/// A uniform inlet that meets one wall (periodic span, z max slip) delivers the requested flow, to round-off: the wall's
/// row takes the wall's velocity and delivers nothing, an eighth of an 8-cell section, so without the raise the inlet
/// delivered 0.875 of it (a parabolic one, 0.956).
#[test]
fn a_uniform_inlet_meeting_a_wall_delivers_its_flow() {
    let n = 8;
    let mut d = domain([0.0; 3], [16.0 * DX, n as f64 * DX, n as f64 * DX], DX);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.1, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
    d.faces.y_min = FaceCondition::Periodic {};
    d.faces.y_max = FaceCondition::Periodic {};
    d.faces.z_max = FaceCondition::Slip {};
    let m = model(1.0, 1.0e-3, d, procedure(Some(WallModel::SpaldingBounce), 10, 0.01));
    let setup = Setup::new(&m, None, &[], &Limits::default()).unwrap();
    let flow = fieldcraft_lbm::Flow::new(&setup, &[]).unwrap();
    let requested = setup.to_lattice([0.1, 0.0, 0.0])[0] * (n * n) as f64;
    let d = delivered(&flow);
    assert!((d / requested - 1.0).abs() < 1e-12, "delivered {d}, requested {requested}");
    let report = flow.inlet_corner_report();
    assert!((report.flux_rescale - 8.0 / 7.0).abs() < 1e-12, "{report:?}");
}

/// An inlet with no cell that can deliver flow (a bar across the middle of a 3-cell section leaves only cells that touch
/// the box wall and the bar) is refused: its flow cannot be delivered at any velocity.
#[test]
fn an_inlet_with_no_cell_left_to_deliver_its_flow_is_refused() {
    let bar = obstacle("Bar", cuboid([0.0, DX, -DX], [0.5, 2.0 * DX, 4.0 * DX]), &CUBOID_FACES);
    let m = channel_model((16, 3, 3), false, false, Some(WallModel::SpaldingBounce), vec!["Bar".to_string()]);
    let obstacles = vec![bar];
    let setup = Setup::new(&m, None, &obstacles, &Limits::default()).unwrap();
    let e = match fieldcraft_lbm::Flow::new(&setup, &obstacles) {
        Ok(_) => panic!("no inlet cell can deliver flow"),
        Err(e) => e.to_string(),
    };
    assert!(e.contains("none left to deliver its flow"), "{e}");
}

/// A wall cell the rule makes is a wall of its neighbours' force groups: no link bounces into the inlet's own group.
#[test]
fn the_corner_cells_bounce_into_the_walls_they_touch() {
    let (setup, obstacles) = channel((16, 8, 8), false, false, Some(WallModel::SpaldingBounce), vec![], vec![]);
    let flow = fieldcraft_lbm::Flow::new(&setup, &obstacles).unwrap();
    let bounces = flow.links().iter().filter(|l| matches!(l.rule, fieldcraft_lbm::Rule::Bounce { .. })).count();
    assert!(bounces > 0);
    for l in flow.links() {
        if let fieldcraft_lbm::Rule::Bounce { group, .. } = l.rule {
            assert!(group != 0, "cell {} direction {} bounces into the inlet's group", l.cell, l.q);
        }
    }
}

/// The Mach check covers the raised velocity. On four walls the 8 x 8 parabolic inlet is raised by 1.2 (see
/// `an_inlet_cell_touching_two_walls_is_a_wall_cell`), so a lattice speed of 0.16 (Mach 0.277, inside the limit of
/// 0.3) is refused, naming the raise, and 0.14 (raised to 0.168, Mach 0.291) is accepted.
#[test]
fn a_raise_that_takes_the_inlet_past_the_mach_limit_is_refused() {
    // A velocity of v m/s with a time step of 0.01 s on 0.05 m cells is v / 5 in lattice units.
    let build = |v: f64| {
        let mut d = domain([0.0; 3], [16.0 * DX, 8.0 * DX, 8.0 * DX], DX);
        d.faces.x_min = FaceCondition::Velocity { velocity: [v, 0.0, 0.0], profile: Profile::Parabolic };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
        let m = model(1.0, 1.0e-3, d, procedure(Some(WallModel::SpaldingBounce), 10, 0.01));
        let setup = Setup::new(&m, None, &[], &Limits::default()).unwrap();
        fieldcraft_lbm::Flow::new(&setup, &[])
    };
    let e = match build(0.8) {
        Ok(_) => panic!("the raised inlet is past the Mach limit"),
        Err(e) => e.to_string(),
    };
    assert!(e.contains("must be raised by") && e.contains("Mach"), "{e}");
    build(0.7).unwrap();
}
