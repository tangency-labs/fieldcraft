# Compare quiescent and stream-function starts in the ramped cylinder test

Follow-up to [#1153](https://github.com/tangency-labs/fieldcraft/issues/1153), scoped to the short
GPU reproducer. At fixed σ = 0.99, the Vreman case fails around 3.6 D/U while the 2.5-period inlet
ramp is still rising. The initial field approaches free-stream speed away from the body. Changing
just that field can distinguish a startup mismatch from instability shared by both starts.

## Suggested change

- Add a benchmark option such as `init=stream|rest` to `Case`. The existing `start` value sets the
  start of averaging. Preserve the stream-function default for existing recipes.
  [Case and parser](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/tests/bluff_body_les.rs#L44-L126)
- At the existing `set_equilibrium` call, select either the current field or unit density and zero
  velocity. Keep the ramp runs unseeded. With `dump` enabled, record the last valid field's extrema,
  locations and readback step. If the device stops between readbacks, report that step interval.
  [Initialization](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/tests/bluff_body_les.rs#L245-L261) ·
  [Advance and readback loop](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm-gpu/tests/bluff_body_les.rs#L287-L326)

For the paired run, use the existing ignored test and vary only the new option. The following recipe
is the proposed invocation for its quiescent member. Use `init=stream` for the control.

```sh
FIELDCRAFT_BLUFF_CASE='r:shape=cylinder,sgs=vreman,cs=0.18,sigma=0.99,ramp=2.5,start=5,end=40,perturb=0,init=rest,dump=1' \
  cargo test --release -p fieldcraft-lbm-gpu --test bluff_body_les bluff_body_les_rerun -- --ignored --nocapture
```

The existing one-cell-span CPU cylinder provides a quiescent-start example. Its halfway walls and
Smagorinsky closure differ from the failing 3D case.
[CPU ramp control](https://github.com/tangency-labs/fieldcraft/blob/8802ec4a845556afc5af282004a6b0a75fb04b8a/crates/lbm/tests/inlet_ramp_cylinder.rs#L1-L17)

## Acceptance

- [ ] Check the legacy recipe still selects the stream-function start and a small quiescent case
  starts at rest throughout the fluid.
- [ ] Run the #1153 recipe twice on the same device with only `init` changed. Report the inputs,
  ramp duration, maximum lattice speed, density range and first sampled excursion or failure
  interval for each. A run surviving the 40 D/U reproducer establishes that observation window only.
