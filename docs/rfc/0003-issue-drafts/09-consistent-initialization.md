# Test consistent population initialization on a decaying Taylor–Green vortex

Conditional CPU prototype under [#1153](https://github.com/tangency-labs/fieldcraft/issues/1153) if
an established-flow start is needed after the quiescent ramp comparison. Prescribing velocity and
placing every population at equilibrium leaves the initial non-equilibrium stress unspecified by the
velocity gradient. Mei et al. give a way to obtain consistent pressure and populations for a
prescribed divergence-free field.
[Mei et al.](https://doi.org/10.1016/j.compfluid.2005.08.008)

## Suggested change

- Prototype the paper's iteration inside an LBM crate test module with access to the populations.
  Keep the prescribed periodic velocity field fixed during initialization. Give the iteration a
  residual criterion and a step cap, and report both. The existing HRR Taylor–Green case already
  supplies an analytical pressure field, velocity mode and amplitude check.
  [Taylor–Green fixture](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm/tests/hrr.rs#L110-L145) ·
  [Equilibrium initializer](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm/src/flow.rs#L2279-L2316)
- Use three starts with identical velocity: uniform density at equilibrium, analytical pressure at
  equilibrium, and the paper's initialization. Before comparing physical time histories, account for
  the arbitrary mean pressure in the periodic box. Record the initial stress and pressure transient,
  mass drift and velocity-amplitude decay with the same physical stepper.

## Acceptance

- [ ] On the periodic control, publish initialization residuals, number of iterations and errors
  against analytical pressure and velocity decay for each start. A case reaching the iteration cap
  reports an inconclusive initialization.
- [ ] Decide from the control whether to try the initializer in the established-flow cylinder case.
  The prototype is complete with the analytic comparison and decision. A product initializer or a
  GPU port would be separate work.
