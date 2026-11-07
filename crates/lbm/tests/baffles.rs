//! Our own analytic parallel-plate and isolated-fluid models; no external decks.
mod common;

use common::*;
use fieldcraft_lbm::{Flow, backend::Backend};
use fieldcraft_model::{Axis, BoxFace, FaceCondition, FluidThermal, ThermalCondition, ThinWall};

fn wall(axis: Axis, conductance: f64) -> ThinWall {
    ThinWall { name: "Sheet".into(), axis, position: 0.5, conductance }
}

fn plates(axis: Axis, conductance: f64) -> fieldcraft_model::Model {
    plates_with(axis, conductance, |_| {})
}

fn plates_with(axis: Axis, conductance: f64, edit: impl Fn(&mut fieldcraft_model::FluidDomain)) -> fieldcraft_model::Model {
    let a = axis.index();
    let mut max = [0.125; 3];
    max[a] = 1.0;
    let mut d = domain([0.0; 3], max, 0.125);
    for f in BoxFace::ALL {
        if f.axis() != a {
            *d.faces.get_mut(f) = FaceCondition::Periodic {};
        }
    }
    let mut t = FluidThermal { baffles: vec![wall(axis, conductance)], ..Default::default() };
    *t.faces.get_mut(BoxFace::of(a, false)) = ThermalCondition::Temperature { value: kelvin(1.0) };
    *t.faces.get_mut(BoxFace::of(a, true)) = ThermalCondition::Temperature { value: kelvin(0.0) };
    d.thermal = Some(t);
    edit(&mut d);
    hot_model(1.0, 0.1, 0.05, d, step(100.0, 0.125 * 0.125 / (8.0 * 0.05), None))
}

/// The steady temperature at `y` between plates at 1 (y = 0) and 0 (y = 1), k_f = 0.05, with
/// a sheet at y = ½ passing the flux `q`.
fn series(y: f64, q: f64) -> f64 {
    if y < 0.5 { 1.0 - q * y / 0.05 } else { q * (1.0 - y) / 0.05 }
}

/// Fourier resistances in series: q = ΔT / (L/k_f + 1/G), independent of
/// spacing and sheet thickness. Side temperatures differ by q/G.
#[test]
fn thin_wall_matches_series_conduction_on_every_axis() {
    for axis in [Axis::X, Axis::Y, Axis::Z] {
        for conductance in [0.0, 0.05, 1e12] {
            let m = plates(axis, conductance);
            let (setup, mut flow) = prepare(&m, &[]);
            run_all(&setup, &mut flow);
            let q = if conductance == 0.0 { 0.0 } else { 1.0 / (20.0 + 1.0 / conductance) };
            for cell in 0..flow.grid().len() {
                let x = flow.grid().centre(cell)[axis.index()];
                assert!((flow.temperature(cell).unwrap() - series(x, q)).abs() < 1e-9);
            }
            let report = flow.baffles().remove(0);
            let expected = [kelvin(1.0 - q * 10.0), kelvin(q * 10.0)];
            for (got, exact) in report.temperature.into_iter().zip(expected) {
                assert!((got - exact).abs() < 1e-9, "{got} vs {exact}");
            }
            assert!((report.heat_flux[0] - q).abs() < 1e-10, "{report:?}");
            assert_eq!(report.heat_flux[0], -report.heat_flux[1]);
            assert_eq!(flow.fluid_cells(), 8, "sheet consumes no fluid cells");
            println!("{axis:?}, G={conductance}: flux={}, exact={q}, temperatures={:?}", report.heat_flux[0], report.temperature);
        }
    }
}

#[test]
fn baffle_conserves_heat_and_each_passages_mass() {
    for conductance in [0.0, 0.05, 1e12] {
        let mut d = domain([0.0; 3], [0.5, 1.0, 0.25], 0.125);
        d.faces.x_min = FaceCondition::Periodic {};
        d.faces.x_max = FaceCondition::Periodic {};
        d.faces.z_min = FaceCondition::Slip {};
        d.faces.z_max = FaceCondition::Slip {};
        d.thermal = Some(FluidThermal { baffles: vec![wall(Axis::Y, conductance)], ..Default::default() });
        let m = hot_model(1.0, 0.1, 0.05, d, step(1.0, 0.01, None));
        let (_, mut flow) = prepare(&m, &[]);
        let grid = *flow.grid();
        flow.set_temperature(|c| if grid.centre(c)[1] < 0.5 { 1.0 } else { 0.0 });
        flow.set_equilibrium(|c| if grid.centre(c)[1] < 0.5 { (1.01, [0.03, 0.01, 0.0]) } else { (0.99, [-0.03, -0.01, 0.0]) });
        let sums = |f: &Flow| {
            let mut mass = [0.0; 2];
            let mut heat = [0.0; 2];
            for c in 0..grid.len() {
                let side = usize::from(grid.centre(c)[1] >= 0.5);
                mass[side] += f.moments(c).unwrap().0;
                heat[side] += f.temperature(c).unwrap();
            }
            (mass, heat)
        };
        let (mass, energy) = sums(&flow);
        for _ in 0..100 {
            flow.step().unwrap();
            let (now_mass, now_energy) = sums(&flow);
            for i in 0..2 {
                assert!((mass[i] - now_mass[i]).abs() < 1e-11);
            }
            assert!((energy.iter().sum::<f64>() - now_energy.iter().sum::<f64>()).abs() < 1e-11);
            if conductance == 0.0 {
                assert_eq!(flow.baffles()[0].heat_flux, [0.0, -0.0]);
                for i in 0..2 {
                    assert!((energy[i] - now_energy[i]).abs() < 1e-11);
                }
            }
        }
    }
}

/// Two counter-flowing streams either side of the sheet: each passage is plane Couette flow
/// between its moving box wall and the no-slip sheet, u = U (1 − 2y) on both sides, and with the flow along the sheet the temperature is still the series conduction of
/// `thin_wall_matches_series_conduction_on_every_axis` (no viscous heating in the model).
#[test]
fn counter_flowing_streams_either_side_of_a_conducting_sheet() {
    let speed = 0.1;
    let conductance = 0.05;
    let m = plates_with(Axis::Y, conductance, |d| {
        d.faces.y_min = FaceCondition::MovingWall { velocity: [speed, 0.0, 0.0] };
        d.faces.y_max = FaceCondition::MovingWall { velocity: [-speed, 0.0, 0.0] };
    });
    let (setup, mut flow) = prepare(&m, &[]);
    run_all(&setup, &mut flow);
    let lattice = setup.to_lattice([speed, 0.0, 0.0])[0];
    let q = 1.0 / (20.0 + 1.0 / conductance);
    for cell in 0..flow.grid().len() {
        let y = flow.grid().centre(cell)[1];
        let u = 1.0 - 2.0 * y;
        let got = flow.moments(cell).unwrap().1;
        assert!((got[0] - u * lattice).abs() < 1e-9 * lattice && got[1].abs() < 1e-12 && got[2].abs() < 1e-12, "{y}: {got:?}");
        assert!((flow.temperature(cell).unwrap() - series(y, q)).abs() < 1e-9, "{y}");
    }
    let report = flow.baffles().remove(0);
    assert!((report.heat_flux[0] - q).abs() < 1e-10, "{report:?}");
    // the sheet's no-slip force: equal and opposite shear on its two sides
    let sides: Vec<usize> = flow.groups().iter().enumerate().filter(|(_, g)| g.name.starts_with("baffle/Sheet/")).map(|(i, _)| i).collect();
    let [low, high] = [sides[0], sides[1]].map(|g| flow.forces()[g][0]);
    assert!(low.abs() > 0.0 && (low + high).abs() < 1e-9 * low.abs(), "{low} {high}");
}

/// The same wall resolved: a conducting slab one cell thick (1/32), of conductivity
/// G × thickness, between the same held plates. The thin wall leaves that cell to the fluid,
/// so the two differ by the slab's fluid-layer resistance, thickness / k_f over the total
/// (1/32 / 0.05 / 40 = 1.6 %): the thin wall's flux and side temperatures match the slab's
/// within 2 %.
#[test]
fn thin_wall_matches_a_resolved_conducting_slab() {
    let (n, alpha, conductance) = (32, 0.05, 0.05);
    let dx = 1.0 / n as f64;
    let d = slab_column(n, [1.0, 0.0], 0.0);
    let m = cht_model([1.0, 0.1, alpha], (1.0, [conductance * dx; 3], 1.0), d, step(8.0 / alpha, 0.125 * dx * dx / alpha, None));
    let obstacles = [slab_obstacle([0.5, 0.5 + dx])];
    let (setup, mut flow) = prepare(&m, &obstacles);
    run_all(&setup, &mut flow);
    let hs = setup.thermal.as_ref().unwrap().heat_scale;
    // the heat the hot plate receives is negative: it gives the flux to the fluid
    let slab_flux = -flow.heat()[BoxFace::YMin.index()] * hs / (dx * dx);
    let solid = flow.thermal().unwrap().solid().unwrap();
    assert_eq!(solid.cells.len(), 1);
    // the slab's face temperatures, from the fluid next to it and the flux through the half
    // cells (k_f = α with ρ = c = 1)
    let below = flow.temperature(n / 2 - 1).unwrap() - slab_flux * 0.5 * dx / alpha;
    let above = flow.temperature(n / 2 + 1).unwrap() + slab_flux * 0.5 * dx / alpha;

    let (setup, mut thin) = prepare(&plates(Axis::Y, conductance), &[]);
    run_all(&setup, &mut thin);
    let report = thin.baffles().remove(0);
    let sides = report.temperature.map(|t| t - kelvin(0.0));
    println!("slab: flux {slab_flux:.6}, faces {below:.6} {above:.6}; thin wall: flux {:.6}, sides {sides:?}", report.heat_flux[0]);
    assert!((report.heat_flux[0] / slab_flux - 1.0).abs() < 0.02, "{} vs {slab_flux}", report.heat_flux[0]);
    for (got, slab) in sides.into_iter().zip([below, above]) {
        assert!((got - slab).abs() < 0.02, "{got} vs {slab}");
    }
}
