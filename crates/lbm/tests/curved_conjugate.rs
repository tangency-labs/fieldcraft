//! Measurement-only curved conjugate benchmark; the true-surface algorithm remains parked.
mod common;

use common::*;
use fieldcraft_lbm::{Flow, Limits};

#[test]
fn aligned_slots_and_the_cylinder_match_the_analytic_cell_mask() {
    for n in [20, 40, 80, 160] {
        let (model, obstacles) = curved_cht::model_case(n, 1, 512, 1.0, 0.025);
        let (setup, flow) = prepare(&model, &obstacles);
        let grid = &setup.grid;
        for c in 0..grid.len() {
            let [x, y, _] = grid.centre(c);
            let wall = (x < 0.0 && y > 0.1) || (x > 1.0 && y < 0.9);
            let cylinder = (x - 0.5).powi(2) + (y - 0.5).powi(2) < 0.1_f64.powi(2);
            assert_eq!(flow.is_fluid(c), !wall && !cylinder, "{n}: {x},{y}");
        }
        assert_eq!(flow.open_cells().len(), 2 * n / 10, "the prescribed 0.1 L slots stay aligned");
        let solid = flow.thermal().unwrap().solid().unwrap();
        assert!(solid.cell.iter().all(|c| c.source == 0.0), "the reference cylinder generates no heat");
    }
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn rahman_curved_conjugate_coarse_measurement() {
    let (model, obstacles) = curved_cht::model_case(40, 1, 512, 100.0, 0.025);
    let setup = fieldcraft_lbm::Setup::new(&model, None, &obstacles, &Limits::default()).unwrap();
    let mut flow = Flow::new(&setup, &obstacles).unwrap();
    let result = fieldcraft_lbm::run(&setup, &mut flow, &Limits::default()).unwrap();
    assert_eq!(result.steady.as_ref().unwrap().reason.name(), "steady", "the bounded run must settle, not merely reach its cap");
    let observed = curved_cht::observation(&setup, &result, 1, 512);
    assert!(observed.heat_spread < 1e-3, "the wall heat is stable over the final tenth of the recorded history");
}
