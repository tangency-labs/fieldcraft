//! `inlet_ramp` in an empty tube against the ideal tube of `common::ramp_tube` (isothermal flow by characteristics):
//! the pressure rms left after a ramp of one and of two periods. This is the sharp check of the ramp's residual at
//! the solver's Mach number; `inlet_ramp_cylinder.rs` (heavy) is the bluff-body case of the manual's table.

mod common;

use common::ramp_tube::{REFERENCE_ERROR, RampTube, linear_rms, rest_period, richardson, rms};
use common::*;
use fieldcraft_model::{Collision, FaceCondition, Outflow, Procedure, Profile, Turbulence, WallTreatment};

/// The tube: 192 cells long (12 D at the cylinder case's 16 cells per D), 4 cells high, one cell of span.
const CELLS: usize = 192;
const CELLS_PER_D: f64 = 16.0;
const RE: f64 = 3900.0;
const LATTICE_SPEED: f64 = 0.05;
/// The probe: this cell from the inlet (an eighth of the tube, as the cylinder case's far-field probe).
const PROBE_CELL: usize = 24;
/// Whole rest periods in the rms window, which starts at the ramp's end.
const WINDOW_PERIODS: f64 = 4.0;

/// The pressure coefficient per unit of density: c_s² / (½ U²).
fn cp_per_density(speed: f64) -> f64 {
    1.0 / 3.0 / (0.5 * speed * speed)
}

/// The acoustic length is between two bounds. The velocity and pressure faces set their cells' moments at the cell
/// centres (wet nodes), CELLS − 1 apart; a face that acted half a cell outside its centre, where the box's face is,
/// would add half a cell at each end, CELLS. The reference is the middle, and half the distance between the two
/// bounds' results is a term of the tolerance. The probe moves with the inlet's half cell.
fn tube(placement: f64, speed: f64, ramp_steps: u64, intervals: usize) -> RampTube {
    RampTube { length: (CELLS - 1) as f64 + placement, speed, ramp_steps: ramp_steps as f64, probe: PROBE_CELL as f64 + 0.5 * placement, intervals }
}

/// The ramp and the rms window, in steps, for a ramp of `m` rest periods of the middle length.
fn ramp_and_window(m: f64) -> (u64, u64) {
    let period = rest_period(CELLS as f64 - 0.5);
    ((m * period).round() as u64, (WINDOW_PERIODS * period).round() as u64)
}

/// The ideal tube's density at the probe over the window after the ramp.
fn ideal_density(placement: f64, speed: f64, m: f64, intervals: usize) -> Vec<f64> {
    let (ramp, window) = ramp_and_window(m);
    tube(placement, speed, ramp, intervals).density(ramp + window).split_off(ramp as usize)
}

/// The ideal tube's C_p rms over the window after the ramp.
fn ideal_rms(placement: f64, speed: f64, m: f64, intervals: usize) -> f64 {
    rms(&ideal_density(placement, speed, m, intervals)) * cp_per_density(speed)
}

/// The acoustic length a density signal shows, in cells: its period T, the mean spacing of the upward crossings of its
/// mean, is 4 L / (c_s (1 − Ma²)) with the inlet at its speed (down at c_s + U, back at c_s − U), so
/// L = T c_s (1 − Ma²) / 4.
fn acoustic_length(rho: &[f64], speed: f64) -> f64 {
    let mean = rho.iter().sum::<f64>() / rho.len() as f64;
    let mut ups = Vec::new();
    for i in 1..rho.len() {
        let (a, b) = (rho[i - 1] - mean, rho[i] - mean);
        if a < 0.0 && b >= 0.0 {
            ups.push(i as f64 - 1.0 + a / (a - b));
        }
    }
    assert!(ups.len() >= 2, "the signal has fewer than two upward crossings");
    let period = (ups[ups.len() - 1] - ups[0]) / (ups.len() - 1) as f64;
    period / 3f64.sqrt() * (1.0 - 3.0 * speed * speed) / 4.0
}

/// The solver's density at the probe over the same window in the empty tube, with the cylinder case's collision and
/// outlet.
fn lbm_density(m: f64) -> Vec<f64> {
    let (ramp, window) = ramp_and_window(m);
    let dx = 1.0 / CELLS_PER_D;
    let dt = LATTICE_SPEED * dx;
    let mut d = domain([0.0; 3], [CELLS as f64 * dx, 4.0 * dx, dx], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [1.0, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective, sigma: None };
    d.faces.y_min = FaceCondition::Slip {};
    d.faces.y_max = FaceCondition::Slip {};
    d.faces.z_min = FaceCondition::Periodic {};
    d.faces.z_max = FaceCondition::Periodic {};
    let procedure = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some((ramp + window) as f64 * dt),
        steps: None,
        time_step: Some(dt),
        frames: 1,
        magic: None,
        sigma: None,
        collision: Some(Collision::Hrr),
        turbulence: Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: None }),
        steady: None,
        walls: WallTreatment::Halfway,
        inlet_ramp: Some(ramp as f64 * dt),
        trace_bulk_one: false,
    };
    let m = model(1.0, 1.0 / RE, d, procedure);
    let (setup, mut flow) = prepare(&m, &[]);
    assert_eq!(flow.inlet_ramp_steps(), ramp, "the ramp is not the requested number of steps");
    assert_eq!(setup.steps, ramp + window, "the run is not the ramp and the window");
    let grid = *flow.grid();
    let probe = grid.cell_at([(PROBE_CELL as f64 + 0.5) * dx, 2.5 * dx, dx / 2.0]).expect("the probe is in the tube");
    let mut rho = Vec::with_capacity(setup.steps as usize);
    for _ in 0..setup.steps {
        flow.step().unwrap();
        rho.push(flow.moments(probe).expect("the probe is a fluid cell").0);
    }
    rho.split_off(ramp as usize)
}

/// The reference converges at its observed order, and the Richardson bound on its 192-interval result is within
/// `REFERENCE_ERROR`: the fast tests take that constant as the reference's error. Heavy: the 768-interval runs take
/// about 20 s.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_reference_converges() {
    for m in [1.0, 2.0] {
        let c = |intervals| ideal_rms(0.5, LATTICE_SPEED, m, intervals);
        let (order, bound) = richardson(c(CELLS), c(2 * CELLS), c(4 * CELLS));
        println!(
            "m = {m}: ideal tube rms {:.5} on {CELLS}, {:.5} on {}, {:.5} on {}; observed order {order:.2}; Richardson bound {bound:.2e}",
            c(CELLS),
            c(2 * CELLS),
            c(4 * CELLS),
            2 * CELLS,
            4 * CELLS
        );
        assert!(bound <= REFERENCE_ERROR * c(CELLS), "m = {m}: the reference's error bound {bound} is above {REFERENCE_ERROR}");
    }
}

/// At a vanishing Mach number the ideal tube is linear acoustics: the sum of the odd harmonics the ramp leaves.
#[test]
fn the_reference_is_linear_acoustics_at_a_vanishing_mach_number() {
    let speed = 1e-7;
    for m in [1.0, 2.0] {
        let (ramp, window) = ramp_and_window(m);
        let t = tube(0.5, speed, ramp, CELLS);
        let ideal = ideal_rms(0.5, speed, m, CELLS);
        let linear = linear_rms(t.length, speed, ramp as f64, t.probe, 50);
        // The window is four periods rounded to a whole step: half a step of a signal whose square is at most three
        // times its mean square moves the mean square by 1 / window, the rms by half of that.
        let tolerance = REFERENCE_ERROR + 0.5 / window as f64;
        println!("m = {m}: ideal tube rms {ideal:.6e}, linear acoustics {linear:.6e}, relative tolerance {tolerance:.2e}");
        assert!((ideal - linear).abs() <= tolerance * linear, "m = {m}: {ideal} against {linear}");
    }
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn an_empty_tube_leaves_the_ideal_tubes_residual() {
    for m in [1.0, 2.0] {
        let (ramp, window) = ramp_and_window(m);
        let rho = lbm_density(m);
        let measured = rms(&rho) * cp_per_density(LATTICE_SPEED);
        let ideal = ideal_rms(0.5, LATTICE_SPEED, m, CELLS);
        let (short, long) = (ideal_rms(0.0, LATTICE_SPEED, m, CELLS), ideal_rms(1.0, LATTICE_SPEED, m, CELLS));
        // Where the faces act: half the distance between the two bounds of the acoustic length.
        let placement = 0.5 * (long - short).abs() / ideal;
        // The reference's own error: `the_reference_converges` asserts its Richardson bound is within this constant.
        let reference = REFERENCE_ERROR;
        // The lattice's sound speed is second-order accurate in k Δx: a relative error of at most (k Δx)² (Lallemand &
        // Luo, Phys. Rev. E 61, 6546 (2000), the dispersion relation's expansion for the BGK-type single-relaxation
        // operator; the Hrr collision of this run has a different second-order coefficient, so this is an estimate of
        // the term, not a bound). It moves m by the same share and the residual by 8 m² / (4 m² − 1) of it (the slope of
        // 1 / (4 m² − 1); the cosine is flat at a whole m).
        let k = std::f64::consts::PI / (2.0 * (CELLS as f64 - 0.5));
        let dispersion = 8.0 * m * m / (4.0 * m * m - 1.0) * k * k;
        // Viscous decay of a sound wave over the run: exp(−k² ν t) with the lattice's bulk viscosity 2ν/3
        // (Landau & Lifshitz, Fluid Mechanics, §79: ½ k² (4ν/3 + ν_bulk)). This is the attenuation the ideal tube does
        // not have; the outlet's reflection loss is not separately bounded here.
        let nu = LATTICE_SPEED * CELLS_PER_D / RE;
        let viscous = k * k * nu * (ramp + window) as f64;
        let tolerance = placement + reference + dispersion + viscous;
        println!(
            "m = {m}: empty tube C_p rms {measured:.4}, ideal tube {ideal:.4} (bounds {short:.4} and {long:.4}); relative difference {:.2e}, tolerance {tolerance:.2e} (placement {placement:.2e}, reference {reference:.2e}, dispersion {dispersion:.2e}, viscous {viscous:.2e})",
            (measured - ideal) / ideal
        );
        assert!((measured - ideal).abs() <= tolerance * ideal, "m = {m}: {measured} against the ideal tube's {ideal}");
        // The acoustic length the run shows is between the two bounds. The estimate's own error is what it makes on
        // the ideal tube, whose length is known.
        let length = acoustic_length(&rho, LATTICE_SPEED);
        let estimate_error = (acoustic_length(&ideal_density(0.5, LATTICE_SPEED, m, CELLS), LATTICE_SPEED) - (CELLS as f64 - 0.5)).abs();
        println!("m = {m}: acoustic length {length:.2} cells of {CELLS} (the estimate's error on the ideal tube {estimate_error:.3})");
        assert!(
            length >= (CELLS - 1) as f64 - estimate_error && length <= CELLS as f64 + estimate_error,
            "m = {m}: the acoustic length {length} is outside the {} to {CELLS} cells between the faces' cell centres and the box's faces",
            CELLS - 1
        );
    }
}
