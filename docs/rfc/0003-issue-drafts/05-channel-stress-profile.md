# Measure the shear-stress partition across the matched-wall sampling row

Diagnostic slice for [#1108](https://github.com/tangency-labs/fieldcraft/issues/1108). At Reτ = 1000
on 16 cells per half-height, the velocity profile errors on rows 3, 4 and 5 are −5.65%, +0.05% and
+5.56%. Row 4 is the sample and the first row outside the mixing-length layer. The sign change
suggests a location to investigate.
[Measured rows](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/docs/verification.md#L624-L637)

## Suggested change

- Reuse the GPU test's row sums for mean velocity and resolved covariance. Its current
  subgrid-stress estimate uses a velocity gradient and a bulk Smagorinsky formula. For the matched
  wall, derive the molecular and modeled shear from the force-corrected non-equilibrium flux and
  **the relaxation rate actually used at that cell**, with a stated sign, time level and units. Use
  the layer rate on rows below the sample and the bulk rate above it.
  [Existing budget](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/tests/verification.rs#L2717-L2768) ·
  [Device layer rate](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/src/lbm.wgsl#L1136-L1173) ·
  [Layer derivation](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm/src/bounce.rs#L375-L394)
- Average the three shear components over the same blocks as the velocity profile. Fold the two
  channel halves with the wall-normal sign corrected. In steady fully developed flow compare their
  sum with the line from measured wall shear to zero at the centre. Use the long-run force balance
  to identify a drifting window before interpreting the profile.
  [Liu, Xu & Huang](https://doi.org/10.1016/j.jcp.2025.114029)

## Acceptance

- [ ] Check zero resolved Reynolds shear on a laminar control and the discrete shear balance on a
  fully developed driven channel. Document whether the diagnostic samples a cell centre or a face.
- [ ] For one matched GPU grid, publish row, physical height, effective rate, all three shear terms,
  their sum, target and window status. A discrepancy can motivate a separate closure change.
