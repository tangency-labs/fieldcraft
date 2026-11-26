//! The refined setup: one `Setup` per level, the ghost layers, the viscosity-preserving τ scaling, the refusals and
//! the caps over every level (Qi, Klimach & Roller 2019; Schornbaum & Rüde 2016).

mod common;

use common::*;
use fieldcraft_lbm::setup::AUTO_TAU_MAX;
use fieldcraft_lbm::{Error, GHOST_LAYERS, Limits, Result, Setup};
use fieldcraft_model::{BoxFace, Collision, FaceCondition, FluidDomain, Outflow, Procedure, Profile, Refinement};

/// The 16 × 8 × 4 channel of these tests: walls on every face but z, which is periodic.
fn channel() -> FluidDomain {
    let mut d = domain([0.0; 3], [16.0, 8.0, 4.0], 1.0);
    d.faces.z_min = FaceCondition::Periodic {};
    d.faces.z_max = FaceCondition::Periodic {};
    d
}

fn boxed(levels: Vec<Refinement>) -> FluidDomain {
    let mut d = channel();
    d.refinement = levels;
    d
}

fn level(l: u32, min: [f64; 3], max: [f64; 3]) -> Refinement {
    Refinement { level: l, min, max }
}

/// Two nested boxes: level 1 is cells 4..12, 2..6, 0..4 of the domain, level 2 sits inside it.
fn two_levels() -> FluidDomain {
    boxed(vec![level(1, [4.0, 2.0, 0.0], [12.0, 6.0, 4.0]), level(2, [6.0, 3.0, 0.0], [10.0, 5.0, 4.0])])
}

/// ν = 0.02 and Δt = 1 on level 0 (τ⁺ = 0.56), with the magic parameter 3/16.
fn run(end: f64, dt: f64) -> Procedure {
    step(end, dt, Some(3.0 / 16.0))
}

fn levels_of(d: FluidDomain, procedure: Procedure) -> Result<Vec<Setup>> {
    Setup::levels(&model(1.0, 0.02, d, procedure), None, &[], &Limits::default())
}

fn refusal(d: FluidDomain, procedure: Procedure) -> String {
    match levels_of(d, procedure) {
        Ok(_) => String::from("no refusal"),
        Err(e) => e.to_string(),
    }
}

#[test]
fn the_fluid_solve_refuses_a_refined_domain_as_not_supported_yet() {
    let m = model(1.0, 0.02, two_levels(), run(1.0, 1.0));
    let e = Setup::new(&m, None, &[], &Limits::default()).unwrap_err();
    assert!(matches!(e, Error::Invalid(_)), "{e}");
    assert!(e.to_string().contains("not supported yet: a refined fluid domain"), "{e}");
    let e = fieldcraft_lbm::solve_fluid(&m, None, &[], &Limits::default()).unwrap_err();
    assert!(e.to_string().contains("not supported yet: a refined fluid domain"), "{e}");
}

#[test]
fn a_domain_without_refinement_gives_one_setup_without_a_level() {
    let levels = levels_of(channel(), run(1.0, 1.0)).unwrap();
    assert_eq!(levels.len(), 1);
    assert!(levels[0].level.is_none());
}

#[test]
fn levels_keep_the_viscosity_and_scale_the_lattice() {
    let mut d = two_levels();
    d.acceleration = [0.0, -0.001, 0.0];
    d.initial_velocity = [0.05, 0.0, 0.0];
    let levels = levels_of(d, run(1.0, 1.0)).unwrap();
    assert_eq!(levels.len(), 3);
    let base = &levels[0];
    for (l, s) in levels.iter().enumerate() {
        let scale = f64::from(1u32 << l);
        assert_eq!(s.grid.spacing, 1.0 / scale, "level {l}");
        assert_eq!(s.time_step, 1.0 / scale, "level {l}");
        assert_eq!(s.tau_plus, 3.0 * (0.02 * scale) + 0.5, "level {l}");
        assert_eq!(s.tau_minus, 0.5 + (3.0 / 16.0) / (s.tau_plus - 0.5), "level {l}");
        for a in 0..3 {
            assert_eq!(s.acceleration[a], base.acceleration[a] / scale, "level {l} axis {a}");
            assert_eq!(s.initial_velocity[a], base.initial_velocity[a], "level {l} axis {a}");
        }
        assert_eq!(s.viscosity, 0.02, "level {l}");
        assert_eq!(s.steps, base.steps << l, "level {l}");
    }
}

#[test]
fn a_level_grid_has_two_ghost_layers_inside_each_interface() {
    let levels = levels_of(two_levels(), run(1.0, 1.0)).unwrap();
    let one = levels[1].level.unwrap();
    assert_eq!(GHOST_LAYERS, 2);
    assert_eq!(levels[1].grid.cells, [20, 12, 8]);
    assert_eq!(levels[1].grid.origin, [3.0, 1.0, 0.0]);
    assert_eq!(one.interface, [true, true, true, true, false, false]);
    assert_eq!(levels[1].domain.faces.z_min, FaceCondition::Periodic {});
    assert_eq!(levels[0].level.unwrap().covered, Some(([4, 2, 0], [11, 5, 3])));
    assert_eq!(one.covered, Some(([6, 4, 0], [13, 7, 7])));
    assert_eq!(levels[2].level.unwrap().covered, None);
}

#[test]
fn a_level_that_touches_a_wall_has_no_ghost_layer_there() {
    let levels = levels_of(boxed(vec![level(1, [4.0, 0.0, 0.0], [12.0, 4.0, 4.0])]), run(1.0, 1.0)).unwrap();
    let one = levels[1].level.unwrap();
    assert!(!one.interface[BoxFace::YMin.index()]);
    assert_eq!(levels[1].grid.origin[1], 0.0);
    assert_eq!(levels[1].grid.cells[1], 8 + 2);
}

#[test]
fn an_open_face_touched_by_a_box_is_refused() {
    let mut d = boxed(vec![level(1, [0.0, 0.0, 0.0], [4.0, 4.0, 4.0])]);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.1, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
    let msg = refusal(d, run(1.0, 1.0));
    assert!(msg.contains("not supported yet") && msg.contains("open face `x_min`"), "{msg}");
}

#[test]
fn a_periodic_face_touched_on_one_side_is_refused() {
    let msg = refusal(boxed(vec![level(1, [4.0, 2.0, 0.0], [12.0, 6.0, 2.0])]), run(1.0, 1.0));
    assert!(msg.contains("touches the periodic face `z_min` without spanning the box along z"), "{msg}");
}

#[test]
fn the_mrt_collision_is_refused_on_a_refined_domain() {
    let msg = refusal(two_levels(), collide(1.0, 1.0, Some(3.0 / 16.0), Collision::Mrt));
    assert!(msg.contains("not supported yet: the MRT collision on a refined domain"), "{msg}");
}

#[test]
fn an_obstacle_crossing_an_interface_is_refused_and_one_inside_is_accepted() {
    let mut d = boxed(vec![level(1, [4.0, 2.0, 0.0], [12.0, 6.0, 4.0])]);
    d.obstacles = vec!["Ball".into()];
    let m = model(1.0, 0.02, d, run(1.0, 1.0));
    let across = obstacle("Ball", sphere([4.0, 4.0, 2.0], 1.0, 2), &["skin"]);
    let msg = Setup::levels(&m, None, &[across], &Limits::default()).unwrap_err().to_string();
    assert!(msg.contains("obstacle `Ball` crosses the interface of refinement level 1"), "{msg}");
    let inside = obstacle("Ball", sphere([8.0, 4.0, 2.0], 1.0, 2), &["skin"]);
    assert_eq!(Setup::levels(&m, None, &[inside], &Limits::default()).unwrap().len(), 2);
    let far = obstacle("Ball", sphere([1.5, 4.0, 2.0], 0.5, 2), &["skin"]);
    assert_eq!(Setup::levels(&m, None, &[far], &Limits::default()).unwrap().len(), 2);
}

#[test]
fn the_caps_count_every_level() {
    let level_zero_cells: u64 = 16 * 8 * 4;
    let m = model(1.0, 0.02, two_levels(), run(1.0, 1.0));
    let tight = Limits { max_cells: level_zero_cells + 1, ..Limits::default() };
    let msg = Setup::levels(&m, None, &[], &tight).unwrap_err().to_string();
    assert!(msg.contains("refinement levels have") && msg.contains("cells"), "{msg}");
    let updates = Limits { max_cell_updates: 1, ..Limits::default() };
    let msg = Setup::levels(&m, None, &[], &updates).unwrap_err().to_string();
    assert!(msg.contains("cell updates"), "{msg}");
}

#[test]
fn the_automatic_time_step_keeps_tau_at_most_one_on_the_finest_level() {
    let mut p = run(1.0, 1.0);
    if let Procedure::Fluid { end_time, steps, time_step, .. } = &mut p {
        *end_time = None;
        *steps = Some(8);
        *time_step = None;
    }
    let levels = levels_of(two_levels(), p).unwrap();
    let finest = levels.last().unwrap();
    assert!(finest.tau_plus - AUTO_TAU_MAX <= 8.0 * f64::EPSILON, "{}", finest.tau_plus);
    assert!(finest.tau_plus > 0.99, "{}", finest.tau_plus);
}
