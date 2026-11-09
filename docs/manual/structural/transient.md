# Linear transient and response-spectrum analysis

The original rod regression exercises the finite-element command path, analytical displacement
and stress, modal equivalence, persistence, display units, damping interpolation and direction
combination. Additional regressions cover Duhamel pulse loading, prescribed base displacement, a cantilever
modal impact series, HHT high-frequency dissipation and a published CQC correlation table.
Large-model modal truncation studies and external same-mesh comparison remain separate validation tasks.

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
absolute acceleration is not included.

For uniform prescribed base translation, use
`base_displacement: {direction:[1,0,0],magnitude:"1 m",amplitude:{...},initial_velocity:"0 m/s",initial_acceleration:"0 m/s^2"}`.
It is mutually exclusive with `base_acceleration`. The curve specifies base positions at each
integration time. Base acceleration is recovered using the chosen Newmark beta/gamma
kinematics, then applied through the same full `-M r` vector. Supply both initial base
derivatives explicitly; the displacement curve alone does not define them. Derivative jumps
or coarse displacement sampling may create acceleration oscillations: inspect the time step
and prefer a supplied acceleration history when available. This is uniform rigid base motion,
not independently prescribed support motions. Initial nodal states and all outputs remain
relative. Arbitrary nonzero displacement boundary conditions and reference-point output
currently return named errors. Contact and thermal loads are refused. Massless
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

## References and verification

- N. M. Newmark (1959), [A Method of Computation for Structural Dynamics](https://ascelibrary.org/doi/10.1061/JMCEA3.0000098).
- H. M. Hilber, T. J. R. Hughes and R. L. Taylor (1977), [Improved numerical dissipation](https://doi.org/10.1002/eqe.4290050306).
- A. Der Kiureghian (1981), *A response spectrum method for random vibration analysis of MDF systems*, Earthquake Engineering & Structural Dynamics 9, 419–435.

The command-level regressions in `crates/engine/tests/transient.rs` cover the original
consistent-mass rod, file round trips, units, base displacement/acceleration equivalence,
triangular-pulse Duhamel and Fourier-load responses, and a 40-element cantilever retaining
four bending modes. The cantilever compares against an independently evaluated four-mode
Euler–Bernoulli series (1% peak-relative allowance for FE mesh, shear/rotary inertia and
time integration); this does not establish convergence in retained mode count.

Core regressions check both Newmark discrete periods and the negative-alpha HHT
high-frequency amplification polynomial. The CQC reference is E. L. Wilson's published
five-mode correlation table, including two close pairs, accounting for the printed
frequency intervals (±0.005 rad/s) and correlation intervals (±0.0005); the
command test uses an original rotated two-spring model to verify signed modal cancellation.
No commercial model files or outputs are used in these tests. Closed-form rod values derive
from its reduced consistent-mass FE equations, not continuum longitudinal-wave solutions:
with E=4, rho=3, A=L=1, Mff=1 and Kff=4, a constant 1 N gives
`u=(1-cos(2t))/4`, `v=sin(2t)/2`, `a=cos(2t)` in SI.
The full per-direction mass is `[[1,0.5],[0.5,1]] kg`, total physical mass 3 kg;
constraint elimination keeps Mff=1 rather than redistributing the fixed node's mass.

- E. L. Wilson, [Chapter 15, Table 15.1](https://edwilson.neocities.org/bookshelf/2005/2005%20Chapter%2015%20MODIFIED%202014.pdf), 2014 revision, p.15-11.
- Sandia, [Sierra/SolidMechanics 5.20 Theory Manual](https://www.sandia.gov/app/uploads/sites/315/2025/07/SM_Theory_5_20.pdf), equations 12.5–12.6 (Newmark/HHT convention).
- TU Graz, [Cantilever beam modes](https://lampz.tugraz.at/~hadley/memm/mechanics/cbeam.php): eigenfunctions and frequency equation, used with Duhamel integration for the original pulse test.

## Damping interpolation and direction combinations

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
period inputs in ms, and millimetre output. It exercises the full command, assembly, combination and probe path.

## Scope limits

This is a linear, fixed-step solver. The reusable reduced-matrix integration interface is a
building block for future tangent updates, not a nonlinear equilibrium iteration. Adaptive
stepping, residual vectors, independent support motions and rotational initial states are
not implemented. These limits are separate from the available amplitude, direct/modal,
SRSS/CQC and damping/direction interpolation features. A successful solve does not establish
modal truncation accuracy or an external benchmark agreement. Performance reports must name
the actual integration method, retained modes, time step, increment count and saved frames;
command timing does not include later file serialization.
