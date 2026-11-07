//! A symmetry plane: half of a symmetric model, solved with the plane as a box face, is the
//! whole model's half, with the inputs of the whole model: a powered solid or a fan the plane
//! cuts is entered with its whole power or curve.

mod common;

use common::*;
use fieldcraft_lbm::{Backend, Flow};
use fieldcraft_model::{FaceCondition, Fan, FanDirection, FluidThermal, Outflow, Profile, ThermalCondition};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};

/// A channel 3 × `height` (8 cells per unit), a parabolic inlet at 1 °C, an outlet, the floor
/// held at 1 °C, and a conducting block generating `power` across y ∈ [0.375, `top`]: the whole
/// channel (`height` 1, the roof a wall held at 1 °C too, `top` 0.625) or its lower half
/// (`height` 0.5, the roof the symmetry plane through the block, which reaches past it to
/// `top` 0.625 or is cut at it, `top` 0.5: either way half of one block). Run part way to steady
/// state.
fn heated_channel(height: f64, roof: FaceCondition, power: f64, top: f64) -> (fieldcraft_lbm::Setup, Flow) {
    let (n, nu, alpha) = (8usize, 0.02, 0.02);
    let dx = 1.0 / n as f64;
    let mut d = domain([0.0; 3], [3.0, height, dx], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Parabolic };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.faces.y_max = roof;
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d.obstacles = vec!["Block".into()];
    let mut t = FluidThermal { initial_temperature: kelvin(1.0), ..FluidThermal::default() };
    t.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
    t.faces.y_min = ThermalCondition::Temperature { value: kelvin(1.0) };
    if roof == (FaceCondition::Wall {}) {
        t.faces.y_max = ThermalCondition::Temperature { value: kelvin(1.0) };
    }
    t.solids = vec![conducting("Block", power)];
    d.thermal = Some(t);
    let dt = 0.25 * dx * dx / nu;
    let m = cht_model([1.0, nu, alpha], (2.0, [0.5; 3], 1.0), d, step(40.0, dt, None));
    let block = [obstacle("Block", cuboid([0.75, 0.375, -1.0], [1.25, top, 1.0]), &CUBOID_FACES)];
    let (setup, mut flow) = prepare(&m, &block);
    run_all(&setup, &mut flow);
    (setup, flow)
}

/// The heated channel solved whole and as its lower half with a symmetry plane through the
/// block, both given the block's whole power, the half's block reaching past the plane or cut at
/// it (the usual half model): every cell of the half (fluid and solid) carries the whole
/// channel's velocity and temperature, the block's drag on the half is half the whole's, and its
/// heat into the floor is the whole's into each of floor and roof (each the same to round-off,
/// 1e-12 of its scale asserted (4e-15 measured): the specular reflection and the adiabatic
/// bounce-back are the mirror image's populations, the parabolic inlet peaks at the plane, and
/// the block's power is spread over its cells and their mirror images, so half of it is in the
/// box). A slip roof alone, with its parabola ending at the plane, misses: that is the inlet, not
/// the plane.
#[test]
fn a_half_channel_with_a_symmetry_plane_is_half_the_whole() {
    let power = 0.01;
    let (_, whole) = heated_channel(1.0, FaceCondition::Wall {}, power, 0.625);
    let (half_setup, half) = heated_channel(0.5, FaceCondition::Symmetry {}, power, 0.625);
    assert_eq!(half.power_share(), [0.5], "the plane cuts the block in half");
    assert_eq!(whole.power_share(), [1.0]);
    let grid = *half.grid();
    // The largest differences from the whole channel's cells in velocity and temperature, and
    // the whole's peak speed.
    let difference = |half: &Flow| {
        let (mut du, mut dt, mut u_scale) = (0.0f64, 0.0f64, 0.0f64);
        for cell in 0..grid.cells.iter().product::<usize>() {
            let same = whole.grid().index(grid.coords(cell));
            let (a, b) = (half.moments(cell).map(|m| m.1).unwrap_or_default(), whole.moments(same).map(|m| m.1).unwrap_or_default());
            du = du.max((0..3).map(|k| (a[k] - b[k]).abs()).fold(0.0, f64::max));
            u_scale = u_scale.max(b[0].abs());
            dt = dt.max((half.temperature(cell).unwrap() - whole.temperature(same).unwrap()).abs());
        }
        (du, dt, u_scale)
    };
    let (du, dt, u_scale) = difference(&half);
    let group = |f: &Flow, name: &str, body: Option<&str>| f.groups().iter().position(|g| g.name == name && g.body.as_deref() == body).unwrap();
    let block = |f: &Flow, axis: usize| CUBOID_FACES.iter().map(|name| f.forces()[group(f, name, Some("Block"))][axis]).sum::<f64>();
    let heat = |f: &Flow, name: &str| f.heat()[group(f, name, None)];
    let (whole_drag, half_drag) = (block(&whole, 0), block(&half, 0));
    let (floor, roof, half_floor) = (heat(&whole, "y_min"), heat(&whole, "y_max"), heat(&half, "y_min"));
    // The block's temperature rise over the inlet, the scale of the temperatures.
    let rise = (0..grid.cells.iter().product::<usize>()).map(|c| half.temperature(c).unwrap() - 1.0).fold(0.0, f64::max);
    println!(
        "heated channel, whole vs half with a symmetry plane: velocity {:.2e} of the peak, temperature {:.2e} of the rise ({rise:.4}); \
         block drag whole {whole_drag:.6e}, 2 × half {:.6e} (lift whole {:.2e}); floor heat whole {floor:.6e}, roof {roof:.6e}, half {half_floor:.6e}",
        du / u_scale,
        dt / rise,
        2.0 * half_drag,
        block(&whole, 1)
    );
    assert!(du < 1e-12 * u_scale && dt < 1e-12 * rise, "{du} {dt}");
    assert!((2.0 * half_drag - whole_drag).abs() < 1e-12 * whole_drag.abs(), "{half_drag} {whole_drag}");
    assert!(block(&whole, 1).abs() < 1e-12 * whole_drag.abs(), "the whole block has no lift");
    assert!((half_floor - floor).abs() < 1e-12 * floor.abs() && (roof - floor).abs() < 1e-12 * floor.abs(), "{floor} {roof} {half_floor}");
    assert!(!half_setup.warnings.iter().any(|w| w.starts_with("symmetry plane")), "{:?}", half_setup.warnings);

    // The block cut at the plane, its geometry ending there: the same half of one block.
    let (cut_setup, cut) = heated_channel(0.5, FaceCondition::Symmetry {}, power, 0.5);
    assert_eq!(cut.power_share(), [0.5], "a block ending at the plane is half of one block");
    let (du, dt, _) = difference(&cut);
    println!("the block cut at the plane: velocity {:.2e} of the peak, temperature {:.2e} of the rise", du / u_scale, dt / rise);
    assert!(du < 1e-12 * u_scale && dt < 1e-12 * rise, "{du} {dt}");
    assert!((2.0 * block(&cut, 0) - whole_drag).abs() < 1e-12 * whole_drag.abs() && (heat(&cut, "y_min") - floor).abs() < 1e-12 * floor.abs());
    assert!(!cut_setup.warnings.iter().any(|w| w.starts_with("symmetry plane")), "{:?}", cut_setup.warnings);

    let (_, slip) = heated_channel(0.5, FaceCondition::Slip {}, power, 0.625);
    let (du, _, _) = difference(&slip);
    println!("a slip roof with the inlet's parabola ending at it: velocity {:.2e} of the peak", du / u_scale);
    assert!(du > 0.1 * u_scale, "{du}");
}

/// A loop 40 × `height` × 1 (lattice units, periodic along x and z), the floor a wall, driven by
/// a fan across x ∈ [5, 7] and y ∈ [0, 8] with its whole curve: the whole loop (`height` 8, the roof a
/// wall) or its lower half (`height` 4, the roof the symmetry plane, which cuts the fan).
fn fan_loop(height: f64, roof: FaceCondition) -> Flow {
    let mut d = domain([0.0; 3], [40.0, height, 1.0], 1.0);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.faces.y_max = roof;
    d.fans = vec![Fan {
        name: "Fan".into(),
        min: [5.0, 0.0, 0.0],
        max: [7.0, 8.0, 1.0],
        direction: FanDirection::XPlus,
        curve: vec![[0.0, 0.004], [0.1, 0.0]],
        relaxation: Some(0.05),
    }];
    let m = model(1.0, 0.1, d, step(400.0, 1.0, None));
    let (setup, mut flow) = prepare(&m, &[]);
    run_all(&setup, &mut flow);
    flow
}

/// A fan the symmetry plane cuts, entered with the whole fan's curve: the half carries half the
/// whole loop's flow through it at the whole's pressure rise, every cell the whole's velocity
/// (its curve is read at the whole fan's flow, twice that through its half).
#[test]
fn a_fan_cut_by_a_symmetry_plane_runs_on_the_whole_fans_curve() {
    let (whole, half) = (fan_loop(8.0, FaceCondition::Wall {}), fan_loop(4.0, FaceCondition::Symmetry {}));
    let (w, h) = (whole.zone_samples()[0], half.zone_samples()[0]);
    let grid = *half.grid();
    let du = (0..grid.cells.iter().product::<usize>())
        .map(|c| {
            let (a, b) = (
                half.moments(c).map(|m| m.1).unwrap_or_default(),
                whole.moments(whole.grid().index(grid.coords(c))).map(|m| m.1).unwrap_or_default(),
            );
            (0..3).map(|k| (a[k] - b[k]).abs()).fold(0.0, f64::max)
        })
        .fold(0.0, f64::max);
    println!(
        "fan loop, whole vs half: flow {:.6e} and 2 × {:.6e}, rise {:.6e} and {:.6e}, velocity {du:.2e}",
        w.flow, h.flow, w.fan_rise, h.fan_rise
    );
    assert!(w.flow > 0.01, "the fan drives the loop: {}", w.flow);
    assert!((2.0 * h.flow - w.flow).abs() < 1e-12 * w.flow, "{} {}", h.flow, w.flow);
    assert!((h.fan_rise - w.fan_rise).abs() < 1e-12 * w.fan_rise, "{} {}", h.fan_rise, w.fan_rise);
    assert!(du < 1e-12 * w.flow, "{du}");
}

/// A half model whose inputs cross its symmetry plane is told so when it is set up.
#[test]
fn an_asymmetric_half_model_is_warned_about_at_setup() {
    let mut d = domain([0.0; 3], [2.0, 0.5, 0.125], 0.125);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.01, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.faces.y_max = FaceCondition::Symmetry {};
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    let m = model(1.0, 0.02, d, step(1.0, 0.01, None));
    let (setup, _) = prepare(&m, &[]);
    let found: Vec<_> = setup.warnings.iter().filter(|w| w.starts_with("symmetry plane y_max")).collect();
    assert_eq!(found.len(), 1, "{:?}", setup.warnings);
    assert!(found[0].contains("inlet velocity on x_min"), "{}", found[0]);
}

/// A square duct 3 × 1 × 1 (8 cells per unit), a uniform inlet at 1 °C, an outlet, its four sides
/// held at 1 °C, and a conducting block generating `power` across y, z ∈ [0.375, 0.625]: the whole
/// duct or its quarter y, z ∈ [0, 0.5] with symmetry planes y_max and z_max through the block
/// (a quarter of one block). Run part way to steady state.
fn heated_duct(quarter: bool, power: f64) -> Flow {
    let (n, nu, alpha) = (8usize, 0.02, 0.02);
    let dx = 1.0 / n as f64;
    let side = if quarter { 0.5 } else { 1.0 };
    let mut d = domain([0.0; 3], [3.0, side, side], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.obstacles = vec!["Block".into()];
    let mut t = FluidThermal { initial_temperature: kelvin(1.0), ..FluidThermal::default() };
    t.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
    t.faces.y_min = ThermalCondition::Temperature { value: kelvin(1.0) };
    t.faces.z_min = ThermalCondition::Temperature { value: kelvin(1.0) };
    if quarter {
        d.faces.y_max = FaceCondition::Symmetry {};
        d.faces.z_max = FaceCondition::Symmetry {};
    } else {
        t.faces.y_max = ThermalCondition::Temperature { value: kelvin(1.0) };
        t.faces.z_max = ThermalCondition::Temperature { value: kelvin(1.0) };
    }
    t.solids = vec![conducting("Block", power)];
    d.thermal = Some(t);
    let dt = 0.25 * dx * dx / nu;
    let m = cht_model([1.0, nu, alpha], (2.0, [0.5; 3], 1.0), d, step(20.0, dt, None));
    let block = [obstacle("Block", cuboid([0.75, 0.375, 0.375], [1.25, 0.625, 0.625]), &CUBOID_FACES)];
    let (setup, mut flow) = prepare(&m, &block);
    run_all(&setup, &mut flow);
    flow
}

/// The largest differences in velocity and temperature (when both carry heat) between every
/// cell of `part` and the same cell of `whole`.
fn differences(part: &Flow, whole: &Flow) -> (f64, f64) {
    let grid = *part.grid();
    let (mut du, mut dt) = (0.0f64, 0.0f64);
    for cell in 0..grid.cells.iter().product::<usize>() {
        let same = whole.grid().index(grid.coords(cell));
        let (a, b) = (part.moments(cell).map(|m| m.1).unwrap_or_default(), whole.moments(same).map(|m| m.1).unwrap_or_default());
        du = du.max((0..3).map(|k| (a[k] - b[k]).abs()).fold(0.0, f64::max));
        if let (Some(ta), Some(tb)) = (part.temperature(cell), whole.temperature(same)) {
            dt = dt.max((ta - tb).abs());
        }
    }
    (du, dt)
}

/// A quarter model, both planes through the heated block, given the block's whole power: a
/// quarter of it is in the box, and every cell carries the whole duct's velocity and temperature.
#[test]
fn a_quarter_duct_is_a_quarter_of_the_whole() {
    let power = 0.01;
    let (whole, quarter) = (heated_duct(false, power), heated_duct(true, power));
    assert_eq!(quarter.power_share(), [0.25], "both planes cut the block");
    let (du, dt) = differences(&quarter, &whole);
    let grid = *quarter.grid();
    let rise = (0..grid.cells.iter().product::<usize>()).filter_map(|c| quarter.temperature(c)).map(|t| t - 1.0).fold(0.0, f64::max);
    println!("heated duct, whole vs quarter: velocity {du:.2e}, temperature {:.2e} of the rise ({rise:.4})", dt / rise);
    assert!(rise > 0.0 && du < 1e-12 * 0.05 && dt < 1e-12 * rise, "{du} {dt} {rise}");
}

/// The fan loop in 3D: 40 × `side` × `side` (lattice units, periodic along x), its four sides
/// walls, driven by a fan across x ∈ [5, 7] and y, z ∈ [0, 8] with its whole curve: the whole
/// loop (`side` 8) or its quarter (`side` 4, symmetry planes y_max and z_max through the fan).
fn fan_duct(quarter: bool) -> Flow {
    let side = if quarter { 4.0 } else { 8.0 };
    let mut d = domain([0.0; 3], [40.0, side, side], 1.0);
    d.faces.x_min = PERIODIC;
    d.faces.x_max = PERIODIC;
    if quarter {
        d.faces.y_max = FaceCondition::Symmetry {};
        d.faces.z_max = FaceCondition::Symmetry {};
    }
    d.fans = vec![Fan {
        name: "Fan".into(),
        min: [5.0, 0.0, 0.0],
        max: [7.0, 8.0, 8.0],
        direction: FanDirection::XPlus,
        curve: vec![[0.0, 0.004], [0.8, 0.0]],
        relaxation: Some(0.05),
    }];
    let m = model(1.0, 0.1, d, step(300.0, 1.0, None));
    let (setup, mut flow) = prepare(&m, &[]);
    run_all(&setup, &mut flow);
    flow
}

/// A fan both planes cut, given the whole fan's curve: the quarter carries a quarter of the
/// whole duct's flow at the whole's pressure rise.
#[test]
fn a_fan_in_a_quarter_duct_runs_on_the_whole_fans_curve() {
    let (whole, quarter) = (fan_duct(false), fan_duct(true));
    let (w, q) = (whole.zone_samples()[0], quarter.zone_samples()[0]);
    let (du, _) = differences(&quarter, &whole);
    println!(
        "fan duct, whole vs quarter: flow {:.6e} and 4 × {:.6e}, rise {:.6e} and {:.6e}, velocity {du:.2e}",
        w.flow, q.flow, w.fan_rise, q.fan_rise
    );
    assert!(w.flow > 0.01, "the fan drives the loop: {}", w.flow);
    assert!((4.0 * q.flow - w.flow).abs() < 1e-12 * w.flow, "{} {}", q.flow, w.flow);
    assert!((q.fan_rise - w.fan_rise).abs() < 1e-12 * w.fan_rise, "{} {}", q.fan_rise, w.fan_rise);
    assert!(du < 1e-12 * w.flow, "{du}");
}

/// An obstacle the plane cuts off its middle is warned about at setup; one centred on the plane
/// or ending at it is not.
#[test]
fn an_obstacle_cut_off_its_middle_is_warned_about_at_setup() {
    let mut d = domain([0.0; 3], [2.0, 0.5, 0.125], 0.125);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.faces.y_max = FaceCondition::Symmetry {};
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d.obstacles = vec!["Centred".into(), "Ending".into(), "Lopsided".into()];
    let m = model(1.0, 0.02, d, step(1.0, 0.01, None));
    let obstacles = [
        obstacle("Centred", cuboid([0.5, 0.25, -1.0], [0.75, 0.75, 1.0]), &CUBOID_FACES),
        obstacle("Ending", cuboid([1.0, 0.25, -1.0], [1.25, 0.5, 1.0]), &CUBOID_FACES),
        obstacle("Lopsided", cuboid([1.5, 0.25, -1.0], [1.75, 0.875, 1.0]), &CUBOID_FACES),
    ];
    let (setup, _) = prepare(&m, &obstacles);
    let found: Vec<_> = setup.warnings.iter().filter(|w| w.starts_with("symmetry plane")).collect();
    assert_eq!(found.len(), 1, "{:?}", setup.warnings);
    assert!(found[0].starts_with("symmetry plane y_max: it cuts obstacle `Lopsided` off its middle"), "{}", found[0]);
}
