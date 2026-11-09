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
They default to zero and must satisfy active constraints. Node existence is checked against
the resolved solve mesh; native files load their definition before the separate mesh part. Initial acceleration is obtained
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

Abaqus native export uses `deck.write` for supported dynamics; Nastran and LS-DYNA
still refuse these procedures by name. Earlier handoff decks were separate original inputs;
they are not evidence of native export. See the export scope below.

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

## Finite-step HHT damping and period

For an undamped, unforced linear oscillator, set `Omega = omega * dt`. Using the negative
alpha convention of Sierra equations 12.5–12.6, the discrete equilibrium is

`M*a[n+1] + (1+alpha)*(C*v[n+1]+K*u[n+1]-f[n+1]) - alpha*(C*v[n]+K*u[n]-f[n]) = 0`.

Substitution of a discrete mode `u[n] ~ lambda^n` into this equation and the Newmark
kinematics gives the characteristic polynomial (derived here, not a quoted formula):

`lambda*(lambda-1)^2 + Omega^2*((1+alpha)*lambda-alpha) *
 (beta*lambda^2 + (1/2-2*beta+gamma)*lambda + (1/2+beta-gamma)) = 0`.

For the oscillatory root `lambda=rho*exp(i*theta)`, the numerical period divided by the
physical period is `Omega/theta`, amplitude retained per numerical cycle is
`rho^(2*pi/theta)`, and the equivalent algorithmic damping ratio is
`-ln(rho)/sqrt(ln(rho)^2+theta^2)`. These quantities describe the oscillatory pair after
the parasitic root has decayed; they are not physical material damping. For alpha zero,
`rho=1` and `theta=2*atan(Omega/2)` recover average-acceleration Newmark. At high frequency
the three roots approach `alpha/(1+alpha)` and the repeated
`-(1+alpha)/(1-alpha)`, with spectral radius `(1+alpha)/(1-alpha)` in the permitted range.
Finite-step phase error must be checked separately: high-frequency stability alone does
not establish time-history accuracy. These expressions apply to a fixed-step, linear,
unforced oscillator with zero physical damping; the full damped or forced response is
not a single decaying sinusoid.

## Original command-chain examples

[Modal history](examples/modal-close-modes.json) and
[CQC with 100/40/40 directions](examples/spectrum-close-modes.json) contain complete
original model-building commands in SI: a 1 kg point mass and two orthogonal oblique
springs, with angular frequencies 13.87 and 13.93 rad/s. Run the recipe through the CLI
or a session, then `solve.transient {step:"ModalPulse"}` or
`solve.spectrum {step:"Spectrum"}` and `solve.probe {nodes:[1],step:...}`.
Use `file.save_as {path:...,results:true}` and `file.open` to retain and read the actual fields.

The modal recipe uses two modes, damping ratios 0.02 and 0.05, a constant 1 N X force,
zero initial displacement/velocity, dt=0.0005 s and 2000 increments, storing 201 samples.
Its maximum component errors against the damped step-load closed form are
1.785e-7 m, 2.373e-6 m/s and 3.474e-5 m/s²; allowances are 1e-6, 1e-5 and 1e-4 in SI.
For modal unit-force responses q1/q2, physical X/Y are their half sum/difference.
The response-spectrum recipe uses constant Sa=1 m/s², damping 0.05, X scale 1 and Y scale
0.6. SRSS and CQC, each with single-axis, directional SRSS and 100/40/40 variants, have
been checked through solve, save, open and probe against the two-mode closed form
with a 1e-10 m absolute allowance. This is a linear discrete-model check, not evidence of
continuum mesh convergence or an external software comparison.

A native command-chain check uses m=1 kg, k=4 N/m, u0=0.001 m, v0=0, no force or
physical damping, and 200 increments for each alpha in {0,-0.1,-1/3} and Omega in
{0.2,0.5,1}. From saved and reopened displacement frames 20–200, the oscillatory recurrence
is fitted and compared with the independent polynomial. The largest discrepancy among
rho, theta, period ratio, equivalent damping and per-cycle amplitude is 2.905e-13
(absolute allowance 1e-9); all three discrete equilibrium/kinematics residuals are below
3e-18 in their SI units (allowance 1e-12).

| alpha | Omega | measured period ratio | measured algorithmic damping | amplitude per cycle |
|---:|---:|---:|---:|---:|
| 0 | 0.2 | 1.00332449071 | 4.622742896e-14 | 1 |
| 0 | 0.5 | 1.02049703762 | 1.155638894e-14 | 1 |
| 0 | 1.0 | 1.07840521615 | -4.789081203e-16 | 1 |
| -0.1 | 0.2 | 1.00416336495 | 7.947059974e-05 | 0.999500796138 |
| -0.1 | 0.5 | 1.02540631365 | 0.001128788865 | 0.992932697545 |
| -0.1 | 1.0 | 1.09450000982 | 0.006754745798 | 0.95844577793 |
| -0.333333333 | 0.2 | 1.00496977752 | 0.0001450047039 | 0.999089323483 |
| -0.333333333 | 0.5 | 1.03014766172 | 0.002037813892 | 0.987277633442 |
| -0.333333333 | 1.0 | 1.11052104369 | 0.01193140211 | 0.927768955893 |

Damping of order 1e-14 at alpha zero is floating-point roundoff. In particular, damping
high frequencies does not remove period elongation: at alpha=-1/3 and Omega=1 the
measured period ratio is 1.110521, so the step is unsuitable if small phase error is required.

## Native Abaqus dynamics export

`deck.write {format:"abaqus",path:"pulse.inp",steps:["Pulse"],units:"m-kg-s-N-Pa"}`
writes one independent dynamics step per deck through the same units, mass, spring, section
and resolved-load paths as other native exports. Direct transient writes `*DYNAMIC, DIRECT`
with explicit alpha/beta/gamma (including zero algorithmic damping for average acceleration).
Modal transient writes its own mass-normalized frequency extraction and `*MODAL DYNAMIC`,
with explicit mode damping. A constant pseudo-acceleration spectrum writes frequency
extraction, `*SPECTRUM` and `*RESPONSE SPECTRUM`, preserving SRSS/CQC and SRSS/R40 directions.

Time/load units are converted with the model; amplitude tables are not copied in display
units. A common load amplitude can be tabular, smooth-step or a zero-origin Fourier series.
Absent amplitudes become an explicit constant curve so the target's default ramp cannot
change the initial load. All three exports require zero initial displacement and velocity;
the target's initial output frame still needs independent interpretation. Transient files
request U/V/A with an explicit node set at the model output stride; spectrum files request
U/A peaks. Requested stress is written as element S; this is not the native solver's global
nodal average and must be compared at a common recovery location.

Named refusals cover mixed independent steps, nonzero initial states/support motion, step
Rayleigh damping, different amplitudes per load, modal HHT or linear-acceleration integration,
and unrepresented projected damping. Spectrum period/damping tables, unequal modal damping
and nonconstant ordinates are refused: Abaqus interpolates log-log and extends endpoints,
whereas FieldCraft interpolates linearly and rejects extrapolation. Even for constant spectra,
check extracted modes lie inside the declared band. Modal time integration belongs to the
target solver, not FieldCraft's Newmark recurrence; use the stated time-discretization tolerance.

This writes input only. It neither invokes a commercial solver nor certifies its preprocessing
or numerical response. Public syntax sources:
[DYNAMIC](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-dynamic.htm),
[AMPLITUDE](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-amplitude.htm),
[MODAL DYNAMIC](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-modaldynamic.htm),
[response-spectrum interpolation](https://abaqus.uclouvain.be/English/SIMACAEANLRefMap/simaanl-c-responsespectrum.htm).


### Energy ledger

New time histories expose `kinetic_energy`, `strain_energy`, cumulative
`damping_dissipation`, `external_work` and `energy_balance` in `solve.probe.frames`.
Values are in the active energy unit (SI storage: J), and survive save/open.
Old result files retain their original fields; their missing work history is not
reconstructed from sparse output. Spectra have no energy ledger.

At every integration step, including steps between output frames, the solver uses
`D += dt * v_mid^T C v_mid` and `W += (f_old+f_new)^T (u_new-u_old)/2`.
The balance is `T+V-(T0+V0)+D-W`, with `T=v^T M v/2`, `V=u^T K u/2`.
These are discrete midpoint-power/trapezoidal-work quantities, not exact continuous
integrals. The average-acceleration Newmark identities imply roundoff-level balance
for constant linear symmetric operators. HHT and linear-acceleration residuals also
include algorithmic and quadrature effects; a nonzero residual is not by itself a
physical instability. Nonzero initial energy is subtracted, never counted as load
work. For modal histories the integrated projected damping includes the specified
modal ratios. For base excitation this is the relative-coordinate effective-force
work, not absolute ground input energy. The load interpolation and output cadence
are unchanged.
