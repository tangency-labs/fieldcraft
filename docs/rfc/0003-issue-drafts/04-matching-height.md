# Specify the matched Spalding wall’s sampling distance in grid rows

Implementation slice for [#1108](https://github.com/tangency-labs/fieldcraft/issues/1108). Sampling
at four rows changes the physical matching height when the grid is refined. A controlled study at
fixed height relative to the channel half-height needs one row count shared by setup, collision and
reporting.
[Matching-height rationale](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm/src/bounce.rs#L238-L251) ·
[Kawai & Larsson](https://doi.org/10.1063/1.3678331)

## Suggested change

- Carry an optional, validated integer row count through the fluid step and `model.step.create` for
  the matched wall. Keep the existing serialized wall choice and resolve an omitted count to four.
  Use checked arithmetic when validating that both wall layers and their sample nodes fit the
  domain, before allocating their maps.
  [Step representation](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/model/src/steps.rs#L383-L412) ·
  [Step command](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/engine/src/cmd/define.rs#L454-L463) ·
  [CPU node setup](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm/src/flow.rs#L1143-L1207)
- The GPU packer already derives its layer count from `nodes.rows()`. Pass the resolved count
  through that path. Change the channel driver's printed height and lower bound for its profile
  check to use the same count.
  [GPU layer layout](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/src/native.rs#L1117-L1148) ·
  [Profile comparison](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/tests/matched_channel.rs#L274-L292)

## Acceptance

- [ ] An existing model with the field omitted and a model specifying four produce the same node map
  and results. Invalid counts, incompatible wall selections and grids too narrow for the request
  give actionable errors.
- [ ] On a supported alternate count, verify the sampled node, every mixing-length row and CPU/GPU
  parity. Show two grid/row pairs with the same physical matching fraction, e.g. 4/40 and 5/50, in
  setup/reporting tests. A long DNS comparison is a subsequent study.
