# Close the force budget over the matched-wall GPU averaging window

Follow-up to [#1108](https://github.com/tangency-labs/fieldcraft/issues/1108). The short CPU test
already checks the matched Spalding wall's requested shear and half-cell body-force accounting. All
six long GPU runs that started reached the time cap. A friction coefficient inferred from bulk speed
alone cannot show whether the force balance had settled.
[CPU balance](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm/tests/wall_matching.rs#L132-L177) ·
[Long-run results](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/docs/verification.md#L595-L611)

## Suggested change

- Add an opt-in on-node diagnostic to the matched-channel study. Accumulate requested tangential
  shear, applied force on both walls, body impulse and the change in fluid momentum over the same
  averaging window. For the uniform drive in lattice units the total applied body force is `a × (M +
  ½ Σρ_sample)`. Include the half-cell contribution exactly once. Use the CPU control's
  stored-momentum definition, including populations streamed onto wall nodes.
  [CPU momentum inventory](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm/tests/wall_matching.rs#L78-L111) ·
  [GPU node bounce](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/src/lbm.wgsl#L815-L852)
- Preserve AA parity and sum the actual per-step impulses over each interval. A single endpoint
  force multiplied by the readback stride cannot verify a transient balance. Reduce scalar totals on
  the device between readbacks. The current GPU population-momentum diagnostic explicitly excludes
  on-node walls.
  [GPU diagnostic](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/src/native.rs#L1903-L1927) ·
  [Sampling loop](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/tests/matched_channel.rs#L177-L216)

## Acceptance

- [ ] A short driven CPU/GPU case exercises even and odd steps. Check requested versus applied shear
  and `ΔP = body impulse − wall impulse` against documented precision-aware tolerances.
- [ ] Print the averaging interval, mean wall and driving forces, momentum change, residual and stop
  reason on one long GPU grid. Keep the existing time-cap verdict and physical acceptance criteria.
