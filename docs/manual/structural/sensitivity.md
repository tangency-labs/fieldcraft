# Linear model sensitivities

## Synopsis

`solve.sensitivity.static` differentiates a scalar displacement or rotation,
external-load work, and optionally continuum/rod integration-point or shell
node stress components. `solve.sensitivity.frequency` differentiates simple natural
frequencies and eigenvalues. Both commands use the same registry in scripts,
CLI and MCP, leave the model unchanged and do not replace stored results.

## Parameters

```json
{"command":"solve.sensitivity.static","params":{
  "step":"Load","parameter":{"type":"material","name":"Steel","constant":"young"},
  "node":1,"dof":1,"stress_element":1,"relative_increment":0.00001}}
```

For frequencies use `solve.sensitivity.frequency`, the frequency step name and
the same `parameter` and `relative_increment`; omit node/dof/stress_element.

Parameters identify existing properties: `material` with `young`, `poisson`,
`density`, `e1`, `e2`, `e3`, `nu12`, `nu13`, `nu23`, `g12`, `g13`, `g23`;
`shell_thickness` with the section name; or `spring_stiffness` with the spring
name. Grounded springs represent compliant supports. Changing the parameter
uses the model's ordinary material/section/connector validation. Orthotropic
reciprocal Poisson ratios therefore follow the existing elasticity model.

The central increment is the positive `relative_increment` (default 1e-5,
maximum 0.1) times the absolute parameter value. At zero Poisson ratio the
scale is one. Both perturbed values must remain valid; no one-sided fallback
silently changes the method. `dof` is 1–6; omitted means x translation.
`threads` and work caps are those of `solve.static`.

## Theory and algorithm

The implementation is semi-analytical: central differences of assembled
stiffness and loads followed by differentiation of equilibrium. It does not
subtract independently solved static outputs. For fixed constraints,
K u'=f'−K'u; the reduced right-hand side includes prescribed-displacement
elimination. One original stiffness factor serves the primal, tangent and
transpose-adjoint solves. Density-dependent body loads are differentiated.

Stress recovery evaluates the constitutive operator and the tangent motion
together; it includes the explicit constitutive derivative. External-load
work is W=fᵀu and W'=f'ᵀu+fᵀu'. For moving supports W is not twice the total
strain energy. Thermal and initial fields are currently refused.

For a simple mode, λ'=φᵀ(K'−λM')φ/(φᵀMφ) and f'=f λ'/(2λ). The model is
solved with one additional mode to check the upper gap. Relative gaps at or
below 1e-5 are refused; rigid/near-zero modes and frequency windows are also
refused. Modal mass is recomputed, rather than assumed to equal one.

## Output and units

Values and derivatives use the display system, with an explicit unit for
each field. A derivative's denominator is the display parameter unit: m/Pa
becomes mm/MPa when those are selected. Static `direct` and `adjoint` refer
to the same observation; `stress` and `stress_derivative` follow the requested
element's sampling order (reported in `stress_sampling`), with tensor stress components
`xx, yy, zz, xy, xz, yz`. `primal_relative_residual` and
`tangent_relative_residual` refer to their respective reduced equations.

## Verification

The `sensitivity-spring` oracle differentiates the exact spring/mass equations.
Solver tests cover a traction bar (displacement, continuum stress and work),
Poisson ratio at zero, density-dependent gravity, prescribed support motion
and material scaling of actual FE natural modes. Repeat with several
increments: truncation and cancellation have competing step-size dependence.

## Limitations

This is a linear-model entry point. It does not provide clustered-subspace
derivatives, beam section stress derivatives or
nonlinear/thermal sensitivity. Harmonic sensitivities support direct response or a fixed modal basis; derivatives of the truncation basis are not included. Modal analysis uses consistent mass at the command boundary.
Assembly differencing can lose small derivatives when unrelated contributions
dominate; the residual checks alone do not bound that error. No BASE-scale
performance or third-party agreement is claimed by this page.

## Implementation and references

`model::design` owns parameter identity and validated updates;
`solver::sensitivity::model` owns derivatives; engine commands and the units
boundary handle display values. Existing FE kernels and sparse solvers are reused.

- Adelman and Haftka, *Sensitivity Analysis for Discrete Structural Systems—A
  Survey*, NASA TM-86333 (1984), static-response section and eigenvalue section:
  https://ntrs.nasa.gov/citations/19850006877.
- Fox and Kapoor, *Rates of Change of Eigenvalues and Eigenvectors*, AIAA J.
  6:2426–2429 (1968), doi:10.2514/3.5008.
- Nelson, *Simplified Calculation of Eigenvector Derivatives*, AIAA J.
  14(9):1201–1205 (1976), doi:10.2514/3.7211. See also NASA TM-89127
  (1987), https://ntrs.nasa.gov/api/citations/19870012842/downloads/19870012842.pdf
  for normalization and the limitations near repeated eigenvalues.

## Bounded one-parameter frequency inversion

`solve.inverse.frequency` fits one named parameter to the ordered simple modes
of a frequency step. This is the initial inversion interface, not joint
multi-parameter identification. By default each candidate recomputes actual FE modes and semi-analytical
derivatives. The explicit uniform-Young reuse option below has stricter
model conditions and reports a separate FE evaluation count.

```json
{"command":"solve.inverse.frequency","params":{
  "step":"Modes","parameter":{"type":"material","name":"Steel","constant":"young"},
  "lower":"100 GPa","upper":"300 GPa","observed":["10 Hz","30 Hz"],
  "standard_deviation":["0.1 Hz","0.1 Hz"],"max_iterations":16}}
```

The current parameter is the starting estimate and must lie inside the bounds.
Bare bounds use the parameter display unit; frequencies and known independent
Gaussian standard deviations use the frequency display unit. Typed units are
accepted. Observation count must equal the requested mode count (at most 256).
`max_iterations` is 1–64; each trial can be an expensive full model evaluation.
`max_lanczos_vectors` independently caps each modal solve (default 3000).
`relative_increment` has the sensitivity command's meaning.

The scalar Levenberg–Marquardt equation is solved in scaled parameter
coordinates, projected to bounds and accepted only when the objective decreases.
`regularization` (default zero) adds λ[(p−p_initial)/scale]². Missing observation
standard deviations select relative residuals. Output includes `estimate`,
`predicted`, `residuals`, normalized `objective`, iteration/evaluation counts
and `at_bound`; the input model is not changed. A local stationary point can
retain a nonzero residual, so inspect the residuals as well as termination.

A local 95% normal interval is returned only for known observation standard
deviations, an interior estimate, zero regularization and nonzero local
information. It is conditional on the model and measurement error assumptions,
not a model-discrepancy or global-identifiability certificate. Otherwise it is
null with a warning. No residual-estimated covariance or mode-crossing tracking is implemented. Complex harmonic inversion is described below.

The oracle recovers k=9 N/m from its exact frequency starting at k=4 N/m.
Command tests also perturb that observation by 0.2 of its known standard
deviation and check that the interval contains the generating value; this is
one deterministic coverage example, not a statistical coverage study.

Reference: K. Levenberg (1944), *A Method for the Solution of Certain Non-Linear
Problems in Least Squares*, Q. Appl. Math. 2:164–168, doi:10.1090/qam/10666.

## Complex harmonic response

`solve.sensitivity.harmonic` takes the same named parameter and relative
increment, a direct or modal harmonic step, and an observation `node` / `dof` (1–6).
Each frame returns complex displacement, direct and transpose-adjoint
parameter derivatives, magnitude and locally unwrapped phase derivatives.
Phase is in radians; all derivatives divide by the displayed parameter unit.
The assembly differences include mass, material Rayleigh/loss damping, loads,
and elimination of real prescribed support motion. Step damping is held fixed for physical parameters; explicit damping parameters differentiate it analytically. Try several increments near resonance. A zero response has undefined
phase; unresolved cancellation or a singular operator returns an error.
Modal response uses a fixed modal/Ritz basis, with applicability described below.
The direct primal, tangent and adjoint currently factor separately.

`solve.inverse.harmonic` fits one named parameter within `lower` / `upper`.
`observed` is a flattened sequence of real, imaginary displacement pairs in
excitation-frequency order. Typed displacement or rotation units are accepted;
`standard_deviation`, when provided, has one known independent Gaussian error
per component. Without known errors the fit uses the magnitude of each observed
phasor as a relative scale and reports no confidence interval. Cartesian
residuals avoid phase wrapping, but do not remove local minima near resonance.
The initial value comes from the model; neither command changes it. The bounds,
regularization and local interval conditions follow the frequency fit above.
This is a scalar fit, not joint material identification.

For simple modes, `solve.sensitivity.frequency` also accepts `mode_shapes:true`.
It uses Nelson's constrained particular solution and a mass-normalization
correction. The returned `vector` contains normalized mode translations,
rotations, reference-point components and their derivatives. Keys ending in
`_si` deliberately stay in the original SI equation basis; vector derivatives
are per **SI** parameter unit, independent of the display-unit scalar outputs.
The normalization is φᵀMφ=1, so density and thickness changes can change the
normalized vector even when its geometric direction is unchanged. The local
sign follows the primal mode. Differentiated equation and normalization
residuals accompany each vector; clustered modes remain unsupported.

For a uniform isotropic Young-modulus fit only,
`solve.inverse.frequency` accepts `uniform_young_scaling:true`. It verifies
that every section uses the named isotropic material and refuses springs,
initial fields, other parameter kinds and nonpositive bounds. With fixed
Poisson ratio, geometry and mass, K(E)=E/E₀ K(E₀) exactly, so frequencies scale
as √(E/E₀). One actual FE sensitivity solve supplies the baseline; bounded LM
then evaluates this exact relation. `model_evaluations:1` distinguishes this
from the default repeated-FE evaluator, and `baseline` includes the actual
frequency derivatives. This option does not apply to mixed-material models or
spring-support identification. It reports no new FE solve for each iteration.

Static `stress_element` also accepts rods and shells. Continuum and rod
samples are global tensors at integration points. Shell samples are global
tensors at element nodes, concatenated bottom, middle, then top. The
`stress_sampling` string identifies this ordering. Differentiation includes
the explicit material/thickness change and the displacement derivative.
Beam section stress recovery is still unsupported.

## Joint bounded identification

`solve.inverse.joint` fits two or three unique existing parameters together,
using mixed static displacement, simple frequency, or direct complex-response
observations. The existing one-parameter commands remain available. Every
number in this command is explicitly SI, independent of the display system:

```json
{"parameters":[
  {"parameter":{"type":"material","name":"Steel","constant":"young"},
   "lower_si":180e9,"upper_si":260e9,"scale_si":210e9},
  {"parameter":{"type":"material","name":"Steel","constant":"density"},
   "lower_si":7000,"upper_si":9000,"scale_si":7850}],
 "observations":[
  {"response":{"kind":"static","step":"Load","node":1,"dof":3},
   "value_si":-1e-6,"standard_deviation_si":1e-9},
  {"response":{"kind":"frequency","step":"Modes","mode":0},
   "value_si":100,"standard_deviation_si":0.1}],
 "regularization":0,"max_iterations":32,"relative_increment":0.00001}
```

The model supplies the initial parameters. Bounds and positive `scale_si`
normalize parameter coordinates; known positive observation standard deviations
whiten the residuals. Frequency `mode` and harmonic `frame` indices start at
zero; `harmonic_real` and `harmonic_imaginary` also take step, node and dof.
The values above illustrate the schema and are not a physical benchmark.

Before adding regularization or LM damping, a one-sided Jacobi SVD checks the
scaled data Jacobian. Rank deficiency or condition number above 1e8 is refused,
including when a prior penalty was requested. Repeated frequencies of a uniform
body constrain E/rho, not E and rho separately. More observations help only if
they add independent parameter information; a mass-independent static response
can supply the missing direction. Local full rank does not prove global uniqueness.

Output includes `estimate_si`, `predicted_si`, `residuals_si`, `jacobian_si`,
parameter scales, scaled singular values, rank and condition number. For an
interior unregularized solution, `covariance_si` is the inverse data information
in original parameter units. The 95% joint confidence region is
(delta p)^T covariance^-1 (delta p) <= `joint_chi_square` (5.99146 for two,
7.81473 for three). `joint_projection_intervals_si` are projections of that
ellipsoid, not separate marginal intervals; correlations matter when checking
whether a parameter vector is inside it. Known independent Gaussian errors and
local linearization are assumptions; model discrepancy is excluded. Regularized
or bound-active solutions return no covariance.

Optional `uniform_material_scaling:true` accepts exactly [Young, density] for
one isotropic material, no springs/point masses, homogeneous prescribed motion,
mass-independent loads and static/frequency observations. Actual FE primal
baselines are followed by exact u~1/E and f~sqrt(E/rho) evaluation, not repeated
FE solves. Mixed materials or elastic supports require the general evaluator.
The general path reuses existing derivatives for every Jacobian column and can
be expensive. This command does not introduce grouped support stiffnesses,
shape variables, nonlinear sensitivity or a differentiated harmonic modal basis.

For a grounded spring targeting a node set, the named stiffness is applied to
each target node. It is a per-node stiffness, not the sum over the set. A common
value on N nodes gives total parallel stiffness N*k; arbitrary weighted groups
are not an additional parameter type. The weakest-direction refusal is in
scaled parameter coordinates; multiply components by the supplied scales to
interpret a physical parameter perturbation.

For inverse commands, `max_iterations` limits optimization only; the independent
`max_lanczos_vectors` option limits each modal basis (default 3000, also bounded
by factor-memory capacity). Sensitivity and joint modal baselines that reach a
truncated basis cap are refused rather than treating capped Ritz values as
converged observations. Increase the modal budget, not the optimization count.

### Measured real modes and reference-normalized components

`solve.inverse.joint` optionally accepts `modal_observations`, a list of groups:

```json
{"id":"test","step":"Modes",
 "sensors":[{"node":10,"dof":3},{"node":20,"dof":3}],
 "measured_modes":[[1,0.7]],"min_mac":0.8,"min_assignment_margin":0.01,
 "normalization":{"kind":"reference","reference_sensor":0}}
```

Use `{"kind":"matched_frequency","group":"test","measured_mode":0}`
for a frequency, or `{"kind":"mode_shape","group":"test","measured_mode":0,"sensor":1}`
for a dimensionless normalized component. `value_si` and `standard_deviation_si`
are the measured frequency in Hz or the **reference-normalized ratio** and its
known standard deviation. Indices start at zero. The reference component itself
is constant and is refused as an observation. The measured vectors define mode
identity; their amplitudes need not be mass calibrated. The response value is
still explicitly supplied in the observation list.

Each trial recomputes real sensor-space MAC and a maximum-total-MAC one-to-one
assignment within each group. Groups define independent measurement sets; combine
modes in one group when they must have distinct computed identities. Returned
`modal_pairings` refer to the final estimate, with computed indices, MAC, sign
alignment and the total-MAC gap to the best distinct feasible assignment. Weak
or ambiguous assignments are errors, not a frequency-order fallback. Thresholds
are caller policy, not confidence probabilities. Sensors must be unique and have
consistent translation/rotation dimensions. At most 64 groups, 64 candidate modes
per frequency step and 4096 sensors per group are accepted. Complex measured
modes and projected oblique sensor directions are not supported by this interface.

For r=φ_i/φ_ref the Jacobian uses
r'=(φ_i'−r φ_ref')/φ_ref, with Nelson's mass-normalized derivative; sign and
arbitrary modal scale cancel. A reference component below 1e-8 of the largest
sensor magnitude is refused. The existing simple-eigenvalue gap checks still
apply: MAC tracking across an exchange does not make individual derivatives at
a repeated eigenvalue unique. These groups cannot use uniform-material reuse.

The covariance assumes independent Gaussian errors **in the normalized observation
coordinates**. Dividing noisy measured components by the same noisy reference
usually produces correlated, non-Gaussian errors; supplying their marginal sigmas
does not validate that covariance. Reference selection, sensor mapping and measurement
noise modelling remain explicit input responsibilities.


### Step damping and fixed-basis modal harmonic derivatives

A design parameter may be `{"type":"damping","step":"H","component":{"kind":"alpha"}}`.
Kinds are `alpha` (s^-1), `beta` (s), `loss_factor` (dimensionless), or
`modal_ratio` with a zero-based `mode` (dimensionless). Alpha is inverse time,
not a cyclic or angular natural frequency. Values are finite and nonnegative.
Updating one modal ratio expands any common ratio to per-mode values while retaining
all other modes. Step and material damping contributions remain additive.

With exp(iωt), step derivatives are Z_α=iωM, Z_β=iωK and Z_η=iK.
A modal-ratio derivative is i 2ωω_j Mφ_jφ_jᵀM. These derivatives are analytic,
including at a zero coefficient; they do not centrally perturb damping below zero.
The existing bounded joint and scalar harmonic fits accept these identities.
An alpha/beta pair may be unidentifiable on a single-mode structure; the data-rank
check still applies, and priors cannot supply the missing information.

For a modal step, a fresh basis with the requested residual-vector cycles is built
inside the sensitivity call. Physical operator derivatives are projected once,
then all frequency points solve the small Galerkin derivative system. No basis
cache or out-of-call preprocessing is required. The reported derivative holds the
modal/Ritz basis fixed. Full-basis and Ritz-completed two-coordinate tests compare
to direct derivatives with absolute/relative tolerance 1e-8; a truncated large model
requires a separate truncation assessment. Residual fields are projected-equation
residuals, not full-space error bounds. Prescribed nonzero support motion and
physical-parameter derivatives with nonzero modal ratios are refused in this method.

`solve.sensitivity.harmonic` computes in SI internally. The engine converts its reply
to display response units per display parameter unit, with an explicit `units` map;
phase remains radians. `solve.inverse.joint` observations retain their explicit SI contract. Modal and
direct complex outputs remain observable through the same command.

### Bounded harmonic batches

`solve.sensitivity.harmonic` also accepts `requests`, an array of 1–16 single-sweep requests (each with `step`, `parameter`, `node`, optional `dof` and `relative_increment`). `threads` and work limits remain on the outer request. The reply contains `results` in request order, per-result display units, and `modal_basis_builds`.

A batch retains modal/Ritz vectors only during that immutable-model call. Exact matching of supports, active loads, complex load factors, mode count and Ritz cycles is required for reuse; excitation frequencies and step damping may differ because they do not enter the real basis extraction. Every derivative column still assembles and projects its own physical derivative operators. All preprocessing is inside the command, so end-to-end timing includes the first extraction. Nothing survives into the next command. Different load factors require a distinct basis even if a complete basis would happen to give the same answer.
