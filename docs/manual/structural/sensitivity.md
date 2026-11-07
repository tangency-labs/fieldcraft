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
derivatives, joint-parameter inversion, beam section stress derivatives or
nonlinear/thermal sensitivity. Direct harmonic response is required for
harmonic sensitivities; a truncated modal basis is not differentiated. Modal analysis uses consistent mass at the command boundary.
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
null with a warning. No residual-estimated covariance, joint-parameter fit,
or mode-crossing tracking is implemented. Complex harmonic inversion is described below.

The oracle recovers k=9 N/m from its exact frequency starting at k=4 N/m.
Command tests also perturb that observation by 0.2 of its known standard
deviation and check that the interval contains the generating value; this is
one deterministic coverage example, not a statistical coverage study.

Reference: K. Levenberg (1944), *A Method for the Solution of Certain Non-Linear
Problems in Least Squares*, Q. Appl. Math. 2:164–168, doi:10.1090/qam/10666.

## Complex harmonic response

`solve.sensitivity.harmonic` takes the same named parameter and relative
increment, a direct harmonic step, and an observation `node` / `dof` (1–6).
Each frame returns complex displacement, direct and transpose-adjoint
parameter derivatives, magnitude and locally unwrapped phase derivatives.
Phase is in radians; all derivatives divide by the displayed parameter unit.
The assembly differences include mass, material Rayleigh/loss damping, loads,
and elimination of real prescribed support motion. Global damping is held
fixed. Try several increments near resonance. A zero response has undefined
phase; unresolved cancellation or a singular operator returns an error.
Modal response and derivatives of a truncated modal basis are not supported.
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
