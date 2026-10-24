# CPU fluid checkpoints

`solve.fluid` can periodically save complete CPU restart state to an explicit native file.
This protects a long run against a stopped process; it does not change the fluid scheme.
A project's saved display results are different: they hold output frames, without the
populations and internal histories needed to continue a solve.

Save the original model before starting, then choose the CPU, a fixed thread count and a
restart path. The Solve Fluid form exposes the same options as the command:

```json
{"command":"solve.fluid","params":{"backend":"cpu","threads":3,"checkpoint":"cooling.restart","checkpoint_every":1000}}
```

This writes every 1000 lattice steps and at a sound final stop. The default interval is 64.
The reply reports `checkpoint_every`, the actual `last_checkpoint.path` and `.step`, and
`resumed_from` (null for a fresh run). A returned failure after a successful save names that
last checkpoint while preserving the original error. Cancellation retains the latest
completed periodic save; it does not try to save a partially advanced or unstable step.
`checkpoint_solver_bytes` gives the conservative solver-state bound checked for this run,
including restart/history/monitor workspace beyond the ordinary frame estimate.

After restarting the app or CLI, reopen the original model and resume explicitly:

```json
{"command":"solve.fluid","params":{"backend":"cpu","threads":3,"checkpoint":"cooling.restart","checkpoint_every":2000,"resume":true}}
```

The model's relevant definitions, step settings, geometry content, actual triangulated
obstacles, CPU caps and thread count must match. A changed identity is refused before loading
the state; choose a new path for a changed model. Storage cadence may change on resume because
it is not a physical setting. Relative paths follow the normal command working directory;
automation write roots apply, and resuming or inspecting an existing destination also needs
read authority. No file is discovered automatically. An unrelated output, model or CAD file
is not replaced. Checkpoints are kept after completion; remove them deliberately when no
longer needed. A completed checkpoint returns its saved result without advancing again.
Use a separate checkpoint path for each simultaneous job; cross-process writer locking is
outside this first slice.

The state retains both momentum and thermal population buffers, cell forces and velocities,
convective-outlet history, relaxed fan state, conducting-solid temperatures and interface
working state, the absolute step and output cursors, all retained frames/force/heat/zone
histories and the steady monitor's partial windows and previous samples. Static topology is
rebuilt and checked. A resumed steady calculation keeps its original stopping rule and any
configured thermal-capacity scaling; its iterations are not relabelled as a physical
transient.

Numbers are stored as little-endian IEEE bits. The versioned `FCLBM001` container includes the
input/settings identity and a SHA-256 of the streamed state; lengths, finite values and
expected shapes are validated, and truncation, a bad checksum or trailing bytes are refused.
The restart schema must change when state layout or compatible CPU solver semantics change.
This is not a cross-version migration format.

Native IO uses one 64 KiB buffer and atomically replaces the destination after a successful
flush. Restore fills a single rebuilt flow rather than allocating a second population
snapshot or loading a whole restart into memory. The final checkpoint run checks the rebuilt
flow, saved frames, retained histories, monitor windows and IO buffer against the machine's
memory bound; the preliminary frame estimate explains that it excludes the extra history and
monitor accounting. Individual state arrays are checked against the rebuilt model before
allocation. A failed replacement leaves the prior checkpoint intact.

This is bounded solver-state accounting, not an exact whole-process peak: application,
geometry preparation and allocator/OS overhead remain separate. Each checkpoint rewrites the
complete retained state and can be large; choose its cadence for the available storage and
acceptable lost work. No throughput improvement is claimed.

This first slice supports CPU fluid only. GPU restart, the explicit/nonlinear/transient-heat
solvers, browser persistence, wall-clock scheduling and automatic resume offers remain open
under #653. The existing fluid and conjugate-heat approximations are unchanged.
