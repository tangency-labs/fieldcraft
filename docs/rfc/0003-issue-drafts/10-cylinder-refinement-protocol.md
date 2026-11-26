# Specify the D/40 evidence and D/80 budget for cylinder acceptance

Protocol slice of [#1107](https://github.com/tangency-labs/fieldcraft/issues/1107). Agree on what
the D/40 run can establish before reserving the GB10. A uniform run can measure wake statistics.
Checking coarse-level instability requires the planned hierarchy or a multilevel control at its
operating conditions. The default uniform D/40 box also exceeds the GB10's single population-buffer
limit.
[RFC costs and staging](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/docs/rfc/0003-lbm-les-roadmap.md#L115-L129) ·
[Buffer split, #699](https://github.com/tangency-labs/fieldcraft/issues/699)

## Suggested change

- Prepare a reproducible run sheet on the tracking issue. Derive dimensions, steps and sample
  cadence from the same parameters as the driver. For 24 D × 12 D × 3.125 D at D/40, list 57,600,000
  cells and 4,377,600,000 population bytes, versus the GB10's 4,294,967,292-byte limit. Use the CPU
  for a uniform run, or wait for #699 for that GPU configuration. For the refined case, preflight
  each level's buffers and record spacing, τ, characteristic speed, interface location and
  observation window.
  [Driver defaults and scaling](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/tests/bluff_body_les.rs#L64-L83) ·
  [GPU allocation check](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/src/native.rs#L1038-L1057)
- Specify experimental bands, their sources and the estimator uncertainty for drag, shedding and
  rear-surface recirculation length before running. Use block statistics, probe and lift spectra,
  and the same physical averaging windows for grid comparisons. Mark a short record as inconclusive.
  Keep long modulation sampling on D/16 initially. At the RFC estimate, four periods at refined D/80
  would cost roughly 21 GB10-hours for sampling alone, versus five hours budgeted for 150 D/U total.
  [Current measurements](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/tests/bluff_body_les.rs#L337-L379) ·
  [Parnaudeau et al.](https://doi.org/10.1063/1.2957018) ·
  [Suss & Mary](https://doi.org/10.2514/6.2026-3485)

## Acceptance

- [ ] A maintainer can reproduce the memory and time estimates from the run sheet and identify the
  backend and unresolved dependencies for each run.
- [ ] The sheet names the reference and tolerance for each observable, stationarity and stop rules,
  and which multilevel evidence is needed for D/80. This issue delivers the protocol for approval.
  Execution and hardware allocation remain with #1107.
