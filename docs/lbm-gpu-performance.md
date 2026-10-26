# GPU lattice Boltzmann: streaming layout and performance

Native GPU execution requires Vulkan. The Metal measurements below are historical evidence
for the streaming-layout decision, not a currently supported backend.

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
and prints how long creating the pipelines took. With `FIELDCRAFT_LBM_GPU_PROFILE=1` it also prints
each kernel's GPU time (`Options::profile`, timestamp queries; the throughput is then not
representative, as every dispatch gets its own pass).
