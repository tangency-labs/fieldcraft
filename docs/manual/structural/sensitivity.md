# Linear model sensitivities

## Synopsis

`solve.sensitivity.static` differentiates a scalar displacement or rotation,
external-load work, and optionally a continuum element's integration-point
stress components. `solve.sensitivity.frequency` differentiates simple natural
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
element's integration-point order, with tensor stress components
`xx, yy, zz, xy, xz, yz`. `primal_relative_residual` and
`tangent_relative_residual` refer to their respective reduced equations.

## Verification

The `sensitivity-spring` oracle differentiates the exact spring/mass equations.
Solver tests cover a traction bar (displacement, continuum stress and work),
Poisson ratio at zero, density-dependent gravity, prescribed support motion
and material scaling of actual FE natural modes. Repeat with several
increments: truncation and cancellation have competing step-size dependence.

## Limitations

This is a linear-model entry point. It does not provide Nelson eigenvector
or clustered-subspace derivatives, harmonic model parameter derivatives,
joint-parameter inversion, shell/beam stress derivatives or nonlinear/thermal
sensitivity. Modal analysis uses consistent mass at the command boundary.
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

## Bounded one-parameter frequency inversion

`solve.inverse.frequency` fits one named parameter to the ordered simple modes
of a frequency step. This is the initial inversion interface, not joint
multi-parameter identification. Each candidate recomputes actual FE modes and
semi-analytical derivatives; there is no uniform-material scaling shortcut.

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
mode-crossing tracking or harmonic inverse is implemented yet.

The oracle recovers k=9 N/m from its exact frequency starting at k=4 N/m.
Command tests also perturb that observation by 0.2 of its known standard
deviation and check that the interval contains the generating value; this is
one deterministic coverage example, not a statistical coverage study.

Reference: K. Levenberg (1944), *A Method for the Solution of Certain Non-Linear
Problems in Least Squares*, Q. Appl. Math. 2:164–168, doi:10.1090/qam/10666.
