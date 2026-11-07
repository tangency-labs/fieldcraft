//! A symmetry plane: half of a symmetric model, solved with the plane as a box face, is the
//! whole model's half.

mod common;

use common::*;
use fieldcraft_lbm::Flow;
use fieldcraft_model::{FaceCondition, FluidThermal, Outflow, Profile, ThermalCondition};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};

/// A channel 3 × `height` (8 cells per unit), a parabolic inlet at 1 °C, an outlet, the floor
/// held at 1 °C, and a powered conducting block across y ∈ [0.375, 0.625]: the whole channel
/// (`height` 1, the roof a wall held at 1 °C too) or its lower half (`height` 0.5, the roof the
/// symmetry plane through the block, which cuts it and takes half its power). Run part way to
/// steady state.
fn heated_channel(height: f64, roof: FaceCondition, power: f64) -> (fieldcraft_lbm::Setup, Flow) {
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
    let block = [obstacle("Block", cuboid([0.75, 0.375, -1.0], [1.25, 0.625, 1.0]), &CUBOID_FACES)];
    let (setup, mut flow) = prepare(&m, &block);
    run_all(&setup, &mut flow);
    (setup, flow)
}

/// The heated channel solved whole and as its lower half with a symmetry plane through the
/// block: every cell of the half (fluid and solid) carries the whole channel's velocity and
/// temperature, the block's drag on the half is half the whole's, and its heat into the floor is
/// the whole's into each of floor and roof (each the same to round-off, 1e-12 of its scale asserted (4e-15 measured): the specular
/// reflection and the adiabatic bounce-back are the mirror image's populations, and the
/// parabolic inlet peaks at the plane). A slip roof alone, with its parabola ending at the plane,
/// misses: that is the inlet, not the plane.
#[test]
fn a_half_channel_with_a_symmetry_plane_is_half_the_whole() {
    let power = 0.01;
    let (_, whole) = heated_channel(1.0, FaceCondition::Wall {}, power);
    let (half_setup, half) = heated_channel(0.5, FaceCondition::Symmetry {}, power / 2.0);
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

    let (_, slip) = heated_channel(0.5, FaceCondition::Slip {}, power / 2.0);
    let (du, _, _) = difference(&slip);
    println!("a slip roof with the inlet's parabola ending at it: velocity {:.2e} of the peak", du / u_scale);
    assert!(du > 0.1 * u_scale, "{du}");
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
