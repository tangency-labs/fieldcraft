//! The inlet corner: with an on-node wall law, a cell in a velocity inlet's layer that touches two walls (its in-plane
//! neighbours: a solid cell, or a box face that is a wall) is a wall cell (plain bounce-back, no inlet condition, no
//! wall law). Without the rule the density at such a cell grew until the run died (a slotted channel on 16,000 cells,
//! wall law on, around step 480). The inlet's velocity is raised so its volume flow stays, and an inlet that would
//! lose more than a tenth of its cells is refused. Other wall treatments are unchanged.
mod common;

use common::*;
use fieldcraft_lbm::backend::Backend;
use fieldcraft_lbm::{Limits, Setup};
use fieldcraft_model::{Collision, FaceCondition, HrrSigma, Outflow, Procedure, Profile, Turbulence, WallModel};

const DX: f64 = 0.05;

fn procedure(wall_model: Option<WallModel>, steps: u64, time_step: f64) -> Procedure {
    let sa = Turbulence::SpalartAllmaras { prandtl: 0.9, wall_model, inflow_ratio: 3.0 };
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
/// as given; spalart_allmaras turbulence with `wall_model`, with the given obstacles.
fn channel(
    (nx, ny, nz): (usize, usize, usize),
    y_periodic: bool,
    z_max_slip: bool,
    wall_model: Option<WallModel>,
    names: Vec<String>,
    obstacles: Vec<fieldcraft_lbm::Obstacle>,
) -> (Setup, Vec<fieldcraft_lbm::Obstacle>) {
    let mut d = domain([0.0; 3], [nx as f64 * DX, ny as f64 * DX, nz as f64 * DX], DX);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.1, 0.0, 0.0], profile: Profile::Parabolic };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    if y_periodic {
        d.faces.y_min = FaceCondition::Periodic {};
        d.faces.y_max = FaceCondition::Periodic {};
    }
    if z_max_slip {
        d.faces.z_max = FaceCondition::Slip {};
    }
    d.obstacles = names;
    let m = model(1.0, 1.0e-3, d, procedure(wall_model, 10, 0.01));
    let setup = Setup::new(&m, None, &obstacles, &Limits::default()).unwrap();
    (setup, obstacles)
}

/// All four sides walls: the four edge cells of the inlet face are wall cells, reported by the flow; the cells next to
/// them (one wall each) stay fluid and keep the inlet's rule.
#[test]
fn an_inlet_cell_touching_two_walls_is_a_wall_cell() {
    let (setup, obstacles) = channel((16, 8, 8), false, false, Some(WallModel::SpaldingBounce), vec![], vec![]);
    let flow = fieldcraft_lbm::Flow::new(&setup, &obstacles).unwrap();
    let report = flow.inlet_corner_report();
    assert_eq!((report.cells, report.inlet_cells), (4, 64), "{report:?}");
    assert!((report.area_lost - 4.0 / 64.0).abs() < 1e-12);
    assert!(report.flux_rescale > 1.0 && report.flux_rescale < 1.3, "{report:?}");
    let g = flow.grid();
    for (j, k) in [(0, 0), (7, 0), (0, 7), (7, 7)] {
        assert!(!Backend::is_fluid(&flow, g.index([0, j, k])), "corner ({j},{k}) is a wall cell");
    }
    for (j, k) in [(1, 0), (0, 1), (3, 0), (7, 3), (3, 3)] {
        assert!(Backend::is_fluid(&flow, g.index([0, j, k])), "({j},{k}) keeps the inlet rule");
    }
}

/// Walls on one side only (periodic span, z max slip): the inlet cells touch one wall, none is reclassified, and the
/// inlet is not rescaled.
#[test]
fn an_inlet_cell_touching_one_wall_keeps_the_inlet() {
    let (setup, obstacles) = channel((16, 8, 8), true, true, Some(WallModel::SpaldingBounce), vec![], vec![]);
    let flow = fieldcraft_lbm::Flow::new(&setup, &obstacles).unwrap();
    assert_eq!(flow.inlet_corner_report(), Default::default());
    let g = flow.grid();
    assert!(Backend::is_fluid(&flow, g.index([0, 0, 0])));
    assert!(Backend::is_fluid(&flow, g.index([0, 3, 0])));
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
#[test]
fn a_staircased_inlet_rim_is_wall_cells_with_the_flow_rescaled() {
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

/// An inlet too coarse for its shape (a 5 x 5 section loses its 4 corners, 16 %) is refused with a clear error.
#[test]
fn an_inlet_that_loses_more_than_a_tenth_is_refused() {
    let (setup, obstacles) = channel((16, 5, 5), false, false, Some(WallModel::SpaldingBounce), vec![], vec![]);
    let e = match fieldcraft_lbm::Flow::new(&setup, &obstacles) {
        Ok(_) => panic!("the inlet loses 16 % of its cells"),
        Err(e) => e.to_string(),
    };
    assert!(e.contains("too coarse") && e.contains("16.0 %"), "{e}");
    let (setup, obstacles) = channel((16, 5, 5), false, false, None, vec![], vec![]);
    assert!(fieldcraft_lbm::Flow::new(&setup, &obstacles).is_ok());
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
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
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
        assert!(ok, "{law:?}, four walls {four_walls}: a density left 1 +- 0.3");
    }
}

/// The velocity the inlet delivers is the one it was asked for: with the edge cells walls and the row beside them on the
/// walls' plane (at rest), the cells that keep the inlet carry the whole requested flow, parabolic profile included.
#[test]
fn the_inlet_delivers_the_requested_flow() {
    let (nx, n) = (16, 8);
    let (setup, obstacles) = channel((nx, n, n), false, false, Some(WallModel::SpaldingBounce), vec![], vec![]);
    let flow = fieldcraft_lbm::Flow::new(&setup, &obstacles).unwrap();
    let u0 = setup.to_lattice([0.1, 0.0, 0.0])[0];
    let profile = |i: usize| {
        let s = (i as f64 + 0.5) / n as f64;
        4.0 * s * (1.0 - s)
    };
    let requested: f64 = (0..n).flat_map(|j| (0..n).map(move |k| (j, k))).map(|(j, k)| u0 * profile(j) * profile(k)).sum();
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
