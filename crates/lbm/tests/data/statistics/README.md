# Statistical rejection controls

Original FieldCraft `cb976117` audit, NVIDIA GB10, in-place GPU streaming.
Each line is the signed drag coefficient sampled every 0.1 D/U. These are
failure controls, not successful physical-validation oracles.

- `cylinder-early.csv`: Re3900 N12 HRR cylinder, 100–200 D/U inclusive (1001 samples).
- `cylinder-late.csv`: same trajectory, 400–600 D/U inclusive (2001 samples).
- `cylinder-potential-late.csv`: divergence-free potential-flow start, same case,
  400–600 D/U inclusive (2001 samples).

The analytic geometry, model parameters, seeds, scaling, recipes and reported
statistics are in [issue #796](https://github.com/tangency-labs/fieldcraft/issues/796).
Only the Cd column of the original force histories is retained here. The tests
first reproduce the published whole-window mean, then apply acceptance to
complete blocks; the remaining final partial block is explicitly excluded.
