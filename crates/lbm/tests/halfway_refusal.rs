//! #794: legacy files remain readable, but unstable halfway walls never start a solve.
mod common;
use common::*;
use fieldcraft_lbm::{Choice, Limits, Setup, checkpoint::Store, solve::solve_fluid_checkpointed, solve_fluid_on};
use fieldcraft_model::{Collision, Model, Procedure, Turbulence, WallModel};

struct UnusedStore;
impl Store for UnusedStore {
    fn every(&self) -> u64 {
        panic!("refused setup reached checkpoint cadence")
    }
    fn checked_memory(&mut self, _: u64) {
        panic!("refused setup reached checkpoint allocation")
    }
    fn load(&mut self, _: &mut dyn FnMut(&mut dyn std::io::Read) -> fieldcraft_lbm::Result<u64>) -> fieldcraft_lbm::Result<bool> {
        panic!("refused setup reached checkpoint load")
    }
    fn save(&mut self, _: u64, _: &mut dyn FnMut(&mut dyn std::io::Write) -> fieldcraft_lbm::Result<()>) -> fieldcraft_lbm::Result<()> {
        panic!("refused setup reached checkpoint save")
    }
}

#[test]
fn legacy_halfway_files_are_readable_but_every_solve_route_refuses_before_allocation() {
    for collision in [Collision::Trt, Collision::Mrt, Collision::Regularized, Collision::Hrr] {
        for t in [
            Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: Some(WallModel::SpaldingBounceHalfway) },
            Turbulence::Vreman { constant: 0.18, prandtl: 0.9, wall_model: Some(WallModel::SpaldingBounceHalfway) },
        ] {
            let mut p = collide(1.0, 1.0, None, collision);
            if let Procedure::Fluid { turbulence, .. } = &mut p {
                *turbulence = Some(t);
            }
            let m = model(1.0, 0.01, domain([0.0; 3], [4.0; 3], 1.0), p);
            let text = serde_json::to_string(&m).unwrap();
            assert!(text.contains("spalding_bounce_halfway"));
            let loaded: Model = serde_json::from_str(&text).unwrap();
            assert_eq!(loaded, m, "the old wall selection is retained, never mapped to on-node");
            let limits = Limits { max_cells: 0, max_memory: 1, ..Limits::default() };
            let check = |error: fieldcraft_lbm::Error| {
                let e = error.to_string();
                assert!(e.contains("spalding_bounce_halfway") && e.contains("disabled") && e.contains("remove `wall_model`"), "{e}");
            };
            check(Setup::new(&loaded, None, &[], &limits).unwrap_err());
            for backend in [Choice::Cpu, Choice::Device] {
                check(solve_fluid_on(&loaded, None, &[], &limits, backend).unwrap_err());
            }
            check(solve_fluid_checkpointed(&loaded, None, &[], &limits, &mut UnusedStore).unwrap_err());
        }
    }
}
