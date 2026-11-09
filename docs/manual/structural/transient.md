# Linear transient and response-spectrum analysis

The original rod regression exercises the finite-element command path, analytical displacement
and stress, modal equivalence, persistence, display units, damping interpolation and direction
combination. Independent benchmark validation and large-model modal truncation studies remain pending.

Create a `transient` step with `model.step.create`, then run `solve.transient` or `solve.run`.
The original two-node rod recipe is `crates/engine/tests/data/transient-rod.json`.
It defines a complete finite-element model, rather than supplying matrices to a kernel test.

```json
{"name":"Pulse","procedure":{"type":"transient","time_step":"0.001 s",
 "increments":1000,"output_every":100,"integration":{"type":"average_acceleration"},
 "method":{"type":"direct"},"load_amplitudes":{"force":{"type":"tabular",
 "points":[["0 s",1],["1 s",1]]}}},"boundary_conditions":["clamp","guide"],
 "loads":["force"],"output":["displacement","stress"]}
```

`time_step` and table abscissae use display time units at command boundaries; storage is SI.
Load amplitudes are dimensionless signed multipliers of existing spatial loads. A missing
curve means a constant multiplier of one. Tables are strictly ascending and linearly
interpolated, holding endpoint values outside their domain. `smooth_step` replaces linear
interpolation by the quintic `10s³-15s⁴+6s⁵`. `periodic` uses `origin`, `offset`,
`angular_frequency` (angular velocity units, internally rad/s), and `harmonics` containing
`[cosine coefficient, sine coefficient]` for harmonics starting at one.
All times are local to the step; states are not inherited from preceding steps.

The solver assembles consistent M, K (including surface foundations), and C (dashpots,
material and foundation viscous damping). Step Rayleigh damping adds `alpha M + beta K`.
Structural loss factors are refused, not converted into a time-domain damping constant.
`initial_displacement` and `initial_velocity` map node labels to XYZ translations.
They default to zero and must satisfy active constraints. Initial acceleration is obtained
from `M a0 = f(0+) - C v0 - K u0`; it is not assumed zero.

The fixed-step methods are Newmark average acceleration (`beta=1/4, gamma=1/2`),
linear acceleration (`beta=1/6, gamma=1/2`), and HHT with `alpha` in `[-1/3,0]`,
`beta=(1-alpha)²/4`, `gamma=1/2-alpha`. HHT weights the stiffness, damping and forcing
residuals at both endpoints; mass acceleration remains at the new endpoint.
The effective operator is factored once. Linear acceleration is conditionally stable:
`dt * highest angular frequency <= sqrt(12)` must be checked by the caller.

`method: {type:"modal",modes:20,damping_ratios:[0.02]}` uses the existing mass-normalized
modal basis. Damping ratios are additional to projected viscous damping; off-diagonal
projected operators are retained. Initial states are projected with the mass inner product.
Truncation requires a mode-count convergence check; no claim of full-equation convergence is
made from the small projected system. Residual-vector augmentation is not yet exposed here.

`base_acceleration: {direction:[1,0,0],magnitude:"1 m/s^2",amplitude:{...}}` adds
`-M r a_ground(t)` using the full mass matrix before constraint reduction. Direction must
be a unit vector. Displacement, velocity and acceleration outputs are relative to the base;
absolute acceleration is not included. Nonzero prescribed motion and reference-point
output currently return named errors. Contact and thermal loads are refused. Massless
uncondensed DOFs cannot be initialized by the SPD mass solve and are refused.

`solve.probe {step:"Pulse",nodes:[2]}` returns stored time samples, or select `frame`.
Times, displacement, rotation, velocity, relative acceleration and stress use display units.
Stress is a global nodal average ordered `xx,yy,zz,xy,xz,yz`; shell stress uses the top surface.
Beam section stress is not supported. Frames also contain kinetic and elastic strain energy;
these alone are not a full energy-balance ledger. Output storage has a bounded value budget;
use `output_every` to reduce it. The final state and initial state are always stored.
`.fieldcraft` save/open retains SI frames through a dedicated result kind, checking lengths,
finite values and time ordering on read. Dynamics contours, VTU, and XY UI paths currently
refuse the result; nodal probes are the supported reader.

## Response spectra

```json
{"name":"Spectrum","procedure":{"type":"response_spectrum","modes":1,
 "damping_ratio":0.05,"combination":"cqc","direction":[1,0,0],
 "pseudo_acceleration":[["0.1 Hz","2 m/s^2"],["1 Hz","2 m/s^2"]]},
 "boundary_conditions":["clamp","guide"],"loads":[],"output":["displacement","stress"]}
```

Run `solve.spectrum` and read `solve.probe`. The table is pseudo-acceleration at the declared
modal damping, interpolated linearly in frequency without extrapolation. `Gamma=Phiᵀ M r`
is calculated with full consistent mass; `q_peak=Gamma Sa/omega²`. Component signs remain
intact until SRSS/CQC combination. CQC uses the Der Kiureghian correlation for broadband
stationary excitation; damping must be positive. Rigid modes have no spectral displacement
and are refused. Spatial loads and assembled viscous damping are refused for this procedure;
the spectrum's damping ratio must describe the input spectrum.

A peak frame has no time and no velocity or energy. Component peaks are not simultaneous:
do not compute a signed stress history or fatigue from them. Acceleration is a
pseudo-acceleration peak, not the maximum true absolute acceleration. Residual corrections remain unimplemented.

All three deck writers explicitly refuse these new procedures pending complete export
semantics. The project handoff decks are separate original inputs, not verified exports.

## References and pending verification

- N. M. Newmark (1959), [A Method of Computation for Structural Dynamics](https://ascelibrary.org/doi/10.1061/JMCEA3.0000098).
- H. M. Hilber, T. J. R. Hughes and R. L. Taylor (1977), [Improved numerical dissipation](https://doi.org/10.1002/eqe.4290050306).
- A. Der Kiureghian (1981), *A response spectrum method for random vibration analysis of MDF systems*, Earthquake Engineering & Structural Dynamics 9, 419–435.

Pending tests cover the original consistent-mass rod's step response, stress recovery,
modal/direct agreement, nonzero initial conditions, spectrum participation, corrupt result
rejection, save/open and display-unit changes. Closed-form rod values derive from its reduced
FE equations, not continuum longitudinal-wave solutions. Independent published benchmarks,
HHT dissipation, BASE performance and third-party same-mesh comparison remain outstanding.

## Damping interpolation and direction combinations (unverified candidate)

Use either `pseudo_acceleration` (one frequency table at `damping_ratio`), `tables`
(frequency tables) or `period_tables` (period tables), never multiple forms together.
Each entry of `tables` or `period_tables` is `{damping_ratio, points:[[x,Sa],...]}`.
Tables must ascend by damping ratio, and their points by positive frequency or period.
Frequencies use display frequency units; periods use display time units; Sa uses
acceleration units. Interpolation is linear in the chosen axis first, then linear in
damping ratio. This is an explicit convention, not logarithmic interpolation or an
implicit damping-correction factor. Axis and damping extrapolation are refused.

`modal_damping_ratios` can contain one common ratio or one per requested mode;
empty uses `damping_ratio`. The same ratio selects the spectrum ordinate and enters
CQC correlation, preventing a mismatch between modal damping and the supplied spectrum.

`directions` takes 1–3 `{direction:[x,y,z],scale}` entries with orthogonal unit directions
and nonnegative dimensionless scales. Empty uses the existing single `direction`.
Each direction uses the same declared spectrum with its own scale. Modal signed responses
are combined separately within each direction. `direction_combination:"srss"` then takes
the root sum of squares of direction peaks; `"percent1004040"` takes the envelope of all
100/40/40 permutations, `0.4*sum(peaks)+0.6*max(peaks)`. Missing axes contribute zero.
These are independent-direction design combination rules, not a simulation of correlated
multi-axis time histories. The command still returns component peak estimates only.

The end-to-end rod test uses two oblique orthogonal axes, damping-table interpolation,
period inputs in ms, and millimetre output. It is written but has not run on Dell yet.
