# GPU fan and porous cooling demonstration

The 39-second recording and its still belong in attachments on
[PR #705](https://github.com/tangency-labs/fieldcraft/pull/705), outside git.
The recipes, reports and recording provenance below describe the demonstration without requiring a video download.

The 39-second clip shows the actual FieldCraft window: domain properties, choosing
the GPU backend, solve progress, saved transient frames, the component limit report,
and velocity/temperature plots. Captures use `ui.screenshot`; captions are added
below the unchanged window with the repository's Inter font. H.264, 1440 × 1048,
10 fps; playback of the saved frames is accelerated. The numerical evidence below
comes from command replies, not image measurements.

## Reproduce

Build the app and CLI, then start `fieldcraft --control 7096`. Execute the command
array [gpu-fan-porous.json](../tools/demo-video/gpu-fan-porous.json) through
`engine.script`, or run it with `fieldcraft-cli run`. To reproduce the reference,
replace only `solve.fluid`'s backend with `cpu` (four threads were used here).
The original `heated-flow` sample generates the bracket from the editable
`UPRIGHT_BRACKET` recipe in `crates/engine/src/sample.rs`; no imported CAD or
external result data is used. A GUI recording can select the domain with
`ui.select {"kind":"fluid_domain","name":"WarmAir"}`, open the solve form with
`ui.start {"command":"solve.fluid"}`, and fill `backend: "gpu"` with `ui.dialog`.
Select frames 0–10 with `results.frame` and capture the presented window.

The 130 × 80 × 70 mm air box has 5 mm cells: 5824 total, 5504 fluid, 320 conducting.
Both open ends hold zero gauge pressure. Air starts at 20 °C; the inlet holds
20 °C. The bracket uses aluminium (2700 kg/m³, 200 W/(m K), 900 J/(kg K)), generates
80 W, and has an 80 °C limit. Its real heat capacity is used. The upstream fan
curve runs from (0, 0.0001 Pa) to (0.00005 m³/s, 0), relaxed by 0.01 each step;
the downstream vent has loss coefficient 3. The curve is an illustrative low-speed
input, not measured fan data. This coarse transient demonstrates the implementation;
it is not a grid-converged electronics design or a two-fluid exchanger.

Regularized collision, automatic in-place streaming, 0.01 s steps, 30 s physical
time, 3000 steps, 11 saved frames. The temperature legend is fixed at 20–50 °C.
The final 80 W mostly enters solid thermal storage: final heat to the fluid is
only about 0.3955 W. The clip does not claim a converged operating temperature.

## Recorded comparison

2026-10-23, Windows, NVIDIA GeForce RTX 4090 / Vulkan, Ryzen 9 7950X.
Full command replies, including units and warnings:
[GPU](verification/gpu-zones/cooling-gpu.json),
[CPU](verification/gpu-zones/cooling-cpu.json).

| quantity at 30 s | CPU | GPU |
|---|---:|---:|
| bracket peak, °C | 44.63342285 | 44.63268280 |
| fan volume flow, mm³/s | 21631.783071 | 21631.708741 |
| applied fan rise, Pa | 0.000057001158 | 0.000057001312 |
| heat to fluid, W | 0.395496608 | 0.395483728 |

The peak differs by 0.0007401 °C, flow by 0.000344%, and heat by 0.00326%.
Applied rise comes from the relaxed fan curve; the separate measured pressure
rise samples the adjacent layers and includes the transient flow. The lattice
Mach number is 0.04472. The GPU reply has only its normal single-precision notice.

## Numerical verification

`cargo test -p fieldcraft-lbm-gpu -p fieldcraft-lbm` checks CPU parity, warm-flow
upload, feedback on the first steps, chunk/readback independence, all six fan
directions, curve endpoints, overlapping ownership, strong resistance, heat and
buoyancy, and the existing unzoned cases. GPU test skips are disabled for the
reported run with `FIELDCRAFT_REQUIRE_GPU=1` and `FIELDCRAFT_GPU=hardware`.

Run the analytical pressure/operating-point test with:

```text
cargo test -p fieldcraft-lbm-gpu --test verification porous_pressure_drop_and_fan_operating_point_match_the_analytic_laws -- --ignored --nocapture
```

It checks Darcy–Forchheimer and loss-coefficient drops to 0.5%, and the fan-system
intersection to 0.5% flow / 1% rise, against closed forms, plus CPU parity to 0.02%,
for both collisions and both layouts. The periodic fan reference is
Q = 0.074568322948, Δp = 0.002881475156 in lattice units. Representative regularized
A–B output is Q = 0.07451981306, fan rise = 0.002882201457, porous drop =
0.002885186113. The pressure sampling and finite grid explain the small difference.

The [implementation and hardware throughput](lbm-gpu-performance.md#fans-and-porous-resistance)
document buffer costs, binding limits, precision refusals, and feedback overhead.
