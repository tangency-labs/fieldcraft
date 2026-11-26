# Measure the pressure pulse from a vortex crossing one CPU refinement interface

Test for the two-level coupling in [#1107](https://github.com/tangency-labs/fieldcraft/issues/1107).
The planned Couette case measures steady stress continuity. A passing steady case can still radiate
sound when a vortex crosses a 2:1 interface. This test depends on the CPU multilevel stepper
becoming available.
[Refinement design and implementation status](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/docs/rfc/0003-lbm-les-roadmap.md#L131-L140) ·
[Astoul et al.](https://arxiv.org/abs/2004.11863)

## Suggested change

- Add a bounded CPU study next to the LBM verification tests. Initialize the convected vortex from
  Astoul et al. on a uniform fine grid and through the real two-level stepper using identical
  physical coordinates, vortex strength, background velocity, viscosity and pressure field. The
  existing periodic vortex test shows how to initialize pressure and project a velocity field in
  lattice units.
  [Periodic HRR example](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm/tests/hrr.rs#L110-L145) ·
  [Field readback](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm/src/backend.rs#L19-L50)
- Sample pressure at fixed physical probes and integrated vorticity over the same physical vortex
  core on both grids at synchronized coarse-step times. Compare fine-to-coarse and coarse-to-fine
  crossings separately. Select the window before a reflected or periodic image reaches the probe.
  Report excess pressure relative to the uniform trace and the change in transmitted circulation.

## Acceptance

- [ ] A uniform-stream control crossing the interface retains its density and velocity within the
  coupling scheme's documented error budget. Use a uniform vortex control to quantify the sampling
  and discretization floor.
- [ ] Report both crossing directions, interface position, grid spacings, time window and metrics.
  State any proposed regression tolerances from analysis and control accuracy before measuring the
  interface. A numerical failure should identify which observable and crossing failed.
