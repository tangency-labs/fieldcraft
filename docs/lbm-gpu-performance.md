# GPU lattice Boltzmann: streaming layout and performance

The momentum populations stream in one of two layouts, chosen per flow with
`fieldcraft_lbm_gpu::Options::streaming`:

- **In place** (the AA pattern: Bailey, Myre, Walsh, Lilja and Saar, "Accelerating lattice
  Boltzmann fluid flow simulations using graphics processors", ICPP 2009; Wittmann, Zeiser, Hager
  and Wellein, Comput. Math. Appl. 65 (2013) 924–935, §6). One population buffer. Even steps read
  population `q` of cell `x` at `(q, x)` and write the collided one at `(opposite(q), x)`; odd steps
  read it at `(opposite(q), x − c_q)` and write it at `(q, x + c_q)`. In both, a cell writes exactly
  the addresses it read, so cells never race, and streaming happens inside the collision: 152
  population bytes per cell update.
- **Two buffers** (A–B). Collide in place, then pull into the other buffer: 304 bytes per update.

The wall, reflection and open-face rules are the same per-link rules in both. In place, a gather
pass reads every link's post-collision source into a scratch value before a scatter pass writes,
so no link overwrites another's source; the momentum exchange is taken in the gather of the step
whose forces are read back (a submission's last). The temperature lattice is two-buffer in both.

`Streaming::Auto`, what a solve uses, streams in place: it was faster on every device measured
(the ignored throughput test, million cell updates per second, means of three runs; the 4090 row
is from the first version of this layout, before its addressing was made cheaper):

| device | case | in place | two buffers |
|---|---|---|---|
| RTX 4090, Vulkan | periodic 128³ | 3897 | 2578 |
| RTX 4090, Vulkan | sphere 128³ | 3572 | 2441 |
| GB10 (unified memory), Vulkan | periodic 128³ | 1224 | 656 |
| GB10 (unified memory), Vulkan | sphere 128³ | 1146 | 630 |
| M4 Max, Metal | periodic 128³ | 1364 | 937 |
| M4 Max, Metal | sphere 128³ | 1256 | 872 |

On the 4090 and the GB10 the two-buffer step is bound by memory bandwidth (~0.8 TB/s and ~0.2 of
the GB10's 0.27), so halving the bytes pays directly; on the M4 Max it gains once the in-place
addressing is cheap (a cell's coordinates once, neighbours wrapped by compare and add). With
128³ cells in place also needs one 152 MiB population buffer less. 64³ fits in cache, where the
two layouts are closer and vary more between runs.

## Measuring

```
cargo test --release -p fieldcraft-lbm-gpu --test throughput -- --ignored --nocapture
```

runs each case in both layouts on the device `FIELDCRAFT_GPU` picks, marks the one `Auto` chooses
and prints how long creating the pipelines took (on Metal this reads near zero: the driver
compiles lazily or from its cache, so the cost lands in the first dispatches instead). With `FIELDCRAFT_LBM_GPU_PROFILE=1` it also prints
each kernel's GPU time (`Options::profile`, timestamp queries; the throughput is then not
representative, as every dispatch gets its own pass).

## Fans and porous resistance

`solve.fluid {"backend":"gpu"}` accepts the existing fan curves, Darcy–Forchheimer
and loss-coefficient boxes. It preserves the CPU's ownership of overlapping zones,
direction conventions, curve endpoint clamping, and under-relaxation per **lattice
step**. Before each collision, a workgroup per fan reduces its owned cells' flow
using the current populations and previous fan action, then updates the rise. There
is no host readback in that feedback loop. This runs in AA and A–B, with TRT and
regularized collision, including buoyancy and conducting solids.

The porous velocity solves the implicit half-force quadratic as on the CPU, after
Guo and Zhao, *Lattice Boltzmann model for incompressible flows through porous media*,
PRE 66:036304 (2002), [DOI](https://doi.org/10.1103/PhysRevE.66.036304). The root is
scaled before squaring, so a representable strong Darcy term does not overflow its
single-precision square. Tests include a Darcy coefficient of 10²⁴.

The flags buffer holds ownership followed by 64 bytes per zone and 8 per curve
point, including mutable fan rise. No storage binding is added; reflection tests
still bound every kernel to at most seven storage buffers. Uniform phase data holds
the zone count. The existing partial-sum buffer is reused for requested zone
diagnostics (16 bytes read per zone); pressure uses density deviations directly
to retain small differences. Histories never change subsequent fan feedback.
The memory planner accounts for the metadata, staging and diagnostic readback.
Unrepresentable coefficients and curves whose flow points collapse in f32 are
refused by zone name before selecting a device; buffer caps are checked before
upload. Unzoned collision expressions and memory estimates stay the same.

Measured on 2026-10-23, RTX 4090 / Vulkan, Ryzen 9 7950X, 64³ periodic cells,
TRT / AA, 32 warm-up then 200 timed steps, the normal optimized development profile:

| zones | million cell updates/s | seconds / 200 steps | requested zone report, ms | bytes read / report |
|---|---:|---:|---:|---:|
| none | 4601.76 | 0.0114 | 0.000 | 0 |
| porous | 2678.59 | 0.0196 | 0.800 | 16 |
| fan + porous | 661.44 | 0.0793 | 0.702 | 32 |

Reports average ten queries and include the moments pass, reduction and readback.
Timings are one run, not a portable performance guarantee; fan reductions have a
visible cost even without host feedback. Small thermal models can favor the CPU.
Reproduce with `cargo test -p fieldcraft-lbm-gpu --test throughput
zone_throughput_and_diagnostic_cost -- --ignored --nocapture` and
`FIELDCRAFT_GPU=hardware`. The [cooling demonstration](gpu-fan-porous-video.md)
includes real app captures and CPU/GPU reports.
