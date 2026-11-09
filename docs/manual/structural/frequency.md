# Natural frequencies

Affine swept `hex8` continuum shells use the opt-in `solid_shell` formulation.
See [geometry, physical ply layers, supported analyses and named refusals](../../solvers/solid-shell.md).

Six-node `tri6` shells use [MITC6-a](../../solvers/mitc6.md), with shared curved
directors and layered constitutive/mass integration. Its documented analysis restrictions and coarse-mesh accuracy limits apply.

Curved `quad8` MITC8 supports consistent and positive-diagonal lumped mass for
frequency and harmonic analysis, including ply densities and rotary inertia.
Consistent mass retains offset/density first-moment coupling; lumped mass omits it.
Both include an explicitly artificial drilling inertia. Transient and spectrum steps
remain refused; curved-shell prestressed frequency has the limited
[reference-geometry modal path](../../solvers/curved-shell-modes.md).
Curved prestressed harmonic response remains refused.
See [MITC8 scope](../../solvers/mitc8.md).

`solve.frequency` computes the lowest natural frequencies and mode shapes of a structure,
K φ = ω² M φ, with every boundary condition of the step held at zero, together with each mode's
participation factors and effective masses. The eigenproblem is solved by the shifted block
Lanczos method on the spectral transformation (K − σM)⁻¹M, with full reorthogonalization, on the
sparse Cholesky factor of the shifted stiffness. A structure that is free (or partly free) gets a
small negative shift automatically, so its rigid-body modes come out at zero frequency.

| | |
|---|---|
| **Commands** | `solve.frequency` (and `solve.run`); `solve.probe` (`mode`), `results.modes`, `results.animate`, `results.export.vtu` |
| **Step procedure** | `{"type": "frequency", "modes", "min_hz", "max_hz"}` |
| **Solver** | `crates/solver` (CPU, f64) |
| **Method** | spectral transformation Lanczos (Ericsson & Ruhe 1980), shifted block Lanczos (Grimes, Lewis & Simon 1994), full reorthogonalization (Parlett 1998, ch. 13), dense tridiagonal QR (Golub & Van Loan 2013, §8.3) |
| **Verified by** | `cantilever-frequencies`, `free-beam-frequencies`, `nafems-fv32`; `cantilever_natural_frequencies`, `free_free_beam_has_six_rigid_body_modes`, `simply_supported_plate_deflection_and_frequency`, `cantilever_beam_frequencies`, `crates/solver/tests/one_d.rs`, `eigen::tests` |

## Synopsis

```json
[
  {"command": "model.step.create", "params": {"name": "Modes",
    "procedure": {"type": "frequency", "modes": 10, "max_hz": 500}, "boundary_conditions": ["clamp"]}},
  {"command": "solve.frequency", "params": {"step": "Modes", "mass": "consistent"}},
  {"command": "solve.probe", "params": {"step": "Modes", "mode": 1, "set": "tip"}}
]
```

## Parameters

### The step

| Field | Value | Default | Meaning |
|---|---|---|---|
| `modes` | 1 to 100 000 | — | the number of modes to extract (lowest first, within the range) |
| `min_hz` | ≥ 0 | none | lower bound of the range, cycles per time unit |
| `max_hz` | ≥ 0, above `min_hz` | none | upper bound of the range |
| `boundary_conditions` | names | none | held at zero in the eigenproblem; prescribed values are ignored |
| `loads` | names | none | no effect |

### `solve.frequency`

| Field | Value | Default | Meaning |
|---|---|---|---|
| `step` | step name | the only frequency step | |
| `mass` | `consistent`, `lumped` | `consistent` | the mass matrix |
| `shift` | finite, eigenvalue units (rad per time unit)² | 0 for a supported structure, −10⁻⁶ Σ K_ii / Σ M_ii for a free one | the spectral shift σ; it must lie below every eigenvalue of the free equations |
| `max_dofs` | 1 to 6 000 000 | 6 000 000 | may only be lowered |
| `max_factor_entries` | 1 to 2 000 000 000 | 400 000 000 | factor entries; also bounds the Lanczos basis ([Linear solvers](linear-solvers.md#memory-limits)) |
| `max_iterations` | 1 to 3000 | 3000 | cap on basis vectors, not block iterations (the final block may be partial); may only be lowered |
| `threads` | 1 to 256 | the session's | changes the time, never the result |

Every element with mass needs a density in its material ("element N has no density" otherwise).

## Theory

The undamped free vibration of the discretised structure is the generalized symmetric
eigenproblem

```text
K φ = λ M φ,   λ = ω²,   f = √λ / 2π
```

with K the stiffness of the [element library](elements.md) and M its mass: **consistent**
(M = ∫ ρ Nᵀ N dV on each element's mass rule; shells and beams with their rotary inertia) or
**lumped** (the diagonal scaling of Hinton, Rock & Zienkiewicz 1976: each element's consistent
diagonal scaled to its mass). Point masses and rotary inertia add the same terms to both. Springs
add stiffness; dashpots are ignored with a warning.

**Spectral transformation** (Ericsson & Ruhe 1980). With a shift σ below the wanted eigenvalues,
the operator

```text
(K − σM)⁻¹ M,   θ = 1 / (λ − σ)
```

is symmetric in the M inner product, and the lowest λ become its largest, well separated θ.
**Block Lanczos** (Grimes, Lewis & Simon 1994) builds an M-orthonormal basis Q block by block,
W = (K − σM)⁻¹ M Q_j, A_j = Q_jᵀ M W, and projects the operator onto the block tridiagonal T with
diagonal blocks A_j and off-diagonal blocks B_j. A block of six vectors finds eigenvalues repeated
up to six times, such as the six rigid-body modes of a free body. Every new block is
orthogonalised against all earlier ones and itself, twice, by classical Gram–Schmidt in the M
inner product (full reorthogonalization; Parlett 1998, ch. 13), so no spurious copies appear. The
eigenpairs (θ, s) of T come from Householder reduction to tridiagonal form and the implicit
symmetric QR algorithm with Wilkinson shifts (Golub & Van Loan 2013, algorithms 8.3.1 and 8.3.3).
A Ritz pair has converged when its residual bound ‖B_{j+1} s_last‖ ≤ 10⁻⁹ θ, s_last the last block
of the Ritz vector.

**Participation.** For a mass-normalised mode (φᵀ M φ = 1) the participation factor along
direction d is Γ_d = φᵀ M r_d, r_d the unit translation along d on the free equations (mesh and
reference point translations), and the effective mass is Γ_d². The rigid translations of a free
body carry the whole mass as effective mass.

![The spectral transformation maps the lowest eigenvalues to the largest, well separated θ](../../media/manual/frequency-1.svg)

## Algorithm

1. **Select and check.** The step named, or the only frequency step; a non-finite `shift` is
   refused; a model with contact is refused.
2. **Mesh and reduction.** As for [linear static](linear-static.md#algorithm), with every
   prescribed value held at zero. Every unknown held: "every degree of freedom is held; there is
   nothing to vibrate". The reduction's memory is budgeted for two matrices.
3. **Assemble K and M** on the same node pattern; the total mass is the sum of the
   translation–translation mass entries divided by the number of translation directions (3, or
   2 in a 2D model); rotary inertia and translation–rotation coupling are left out of it.
4. **Factor.** With `shift` given, K − σM is factored by Cholesky. Without one, K is factored;
   if it is not positive definite, σ = −10⁻⁶ Σ diag K / Σ diag M is taken, K − σM factored, and
   the warning "the structure is not fully supported; rigid-body modes appear at (near) zero
   frequency" added.
5. **Lanczos.** The basis may hold at most min(`max_iterations`, `max_factor_entries` / (2 n))
   vectors (n the free equations), and at least min(6, n). The start block is pseudo-random from a
   fixed seed (xorshift64*), so results are reproducible. Each block's six solves run in parallel.
   A vector that vanishes under orthogonalization is replaced by a fresh random direction (up to
   four attempts).
6. **Convergence check** on a schedule, since the dense eigenproblem costs m³: first at
   max(modes + 6, 12) vectors, then every max(m/4, 6) more. Ritz values are taken in decreasing θ
   (increasing λ) until one has not converged; values below `min_hz` are skipped; the iteration
   stops when `modes` values have converged or one passes `max_hz`.
7. **Stop.** When the basis reaches its cap before that: "found k of N modes after m Lanczos
   vectors (the cap); raise the iteration limit or ask for fewer modes". When the basis spans every
   free equation, the converged modes are returned; fewer than `modes` without a `max_hz` gives
   the warning "found k of the N modes asked for".
8. **Modes.** Ritz vectors expanded to every node; frequency √max(λ, 0)/2π; participation and
   effective mass.

## Output

| Key | Meaning |
|---|---|
| `step`, `dofs`, `equations` | the step, its unknowns, the free ones |
| `lanczos_vectors` | the basis size used |
| `total_mass` | the model's mass |
| `modes` | per mode: `mode` (from 1), `frequency` (cycles per display time unit; 0 for a rigid-body mode), `eigenvalue` λ = ω², `effective_mass` [x, y, z] |
| `warnings` | free structure, modes not found, ignored dashpots |

The stored result has each mode's mass-normalised shape per node (translations, and rotations at
shell, beam and reference point nodes) and its participation factors. `solve.probe` with `mode`
gives the shape at nodes; `results.modes` shows a mode and lists every frequency;
`results.animate` swings it; `results.export.vtu` writes `mode_shape` with the mode, frequency and
eigenvalue.

## Worked example

A slender steel cantilever of square section, clamped at one end: its first eight modes against
Euler–Bernoulli beam theory. The square section gives every bending frequency twice (bending in
y and in z), and the block Lanczos iteration returns each pair together.

| | |
|---|---|
| **Geometry** | 200 mm along x, 10 × 10 mm section; one block (`mesh.block`), face `root` at x = 0 |
| **Mesh** | 80 hex20 (20 × 2 × 2), 621 nodes, 1 863 unknowns (1 800 free) |
| **Material** | steel: E = 210 000 MPa, ν = 0.3, ρ = 7.85 × 10⁻⁹ t/mm³ |
| **Boundary conditions** | `root` fixed |
| **Step** | frequency, `modes` 8; `solve.frequency` with the default consistent mass |
| **Units** | mm, t, s, N, MPa |

```json
[
  {"command": "model.units", "params": {"system": "mm-t-s-N-MPa"}},
  {"command": "mesh.block", "params": {"body": "Beam", "kind": "hex20",
    "corners": [[0,0,0],[200,0,0],[200,10,0],[0,10,0],[0,0,10],[200,0,10],[200,10,10],[0,10,10]],
    "divisions": [20, 2, 2], "faces": {"6": "root", "4": "tip"}}},
  {"command": "model.set.create", "params": {"name": "all", "geometry": [{"body": "Beam"}]}},
  {"command": "model.set.create", "params": {"name": "root", "geometry": [{"body": "Beam", "face": "root"}]}},
  {"command": "model.material.create", "params": {"name": "Steel",
    "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}, "density": 7.85e-9}},
  {"command": "model.section.create", "params": {"name": "Beam", "set": "all", "material": "Steel", "kind": {"type": "solid"}}},
  {"command": "model.boundary_condition.create", "params": {"name": "Clamp", "type": "fixed", "target": {"set": "root"}}},
  {"command": "model.step.create", "params": {"name": "Modes",
    "procedure": {"type": "frequency", "modes": 8}, "boundary_conditions": ["Clamp"]}},
  {"command": "solve.frequency", "params": {"step": "Modes"}},
  {"command": "results.modes", "params": {"step": "Modes", "mode": 1}},
  {"command": "results.contour", "params": {"step": "Modes", "field": "displacement"}},
  {"command": "results.deformed", "params": {"step": "Modes", "show": true, "scale": 0.4}}
]
```

Results from the `solve.frequency` reply (54 Lanczos vectors, 0.1 s, Apple M4 Max, CPU). The
references: bending f = (βₙL)² / (2π L²) · √(E I / ρ A) with βₙL = 1.8751, 4.6941, 7.8548;
torsion f = (1/4L) · √(G J / ρ I_p) with J = 0.1406 a⁴ for the square (Timoshenko & Goodier 1970,
§109) and I_p = a⁴/6; axial f = (1/4L) · √(E/ρ).

| Mode | Shape | FieldCraft (Hz) | Reference (Hz) | Difference |
|---|---|---|---|---|
| 1, 2 | first bending, y and z | 209.515, 209.515 | 208.879 | +0.30 % |
| 3, 4 | second bending | 1 298.345, 1 298.345 | 1 309.023 | −0.82 % |
| 5, 6 | third bending | 3 572.927, 3 572.927 | 3 665.303 | −2.52 % |
| 7 | first torsion | 3 717.734 | 3 682.701 | +0.95 % |
| 8 | first axial | 6 480.322 | 6 465.243 | +0.23 % |

| Quantity | FieldCraft | Reference | Difference |
|---|---|---|---|
| `total_mass` | 1.570 × 10⁻⁴ t | ρ A L = 1.570 × 10⁻⁴ t | < 10⁻¹² |
| Effective mass, modes 1, 3, 5 (in their direction) / total mass | 0.611, 0.184, 0.065 | 0.613, 0.188, 0.065 (Euler–Bernoulli mode shapes) | −0.3 %, −2.4 %, −0.1 % |
| Effective mass, mode 8 (x) / total mass | 0.809 | 8/π² = 0.811 (axial bar) | −0.2 % |

With ν = 0 the same mesh gives f₁ = 208.549 Hz (−0.16 %, the value of the
`cantilever-frequencies` case): at ν = 0.3 the fixed root face also holds the section's Poisson
contraction, which stiffens the first mode by 0.46 %.

![First bending mode of the clamped steel cantilever, seen from above, coloured by displacement magnitude](../../media/manual/examples/structural-frequency.svg)

![Second bending mode (1 298 Hz) with its node near x = 0.78 L](../../media/manual/examples/structural-frequency-2.svg)

![The eight natural frequencies of the worked example against beam theory: Euler–Bernoulli bending and axial, Saint-Venant torsion](../../media/manual/frequency-2.svg)

**What to look at.** The first two frequencies agree to nine digits, and their effective masses
lie along y and z: one bending pair, as K and M have the section's symmetry. The second bending
mode has its node at 0.78 L, where the beam crosses its undeformed axis. The bending frequencies
fall below Euler–Bernoulli as the mode number rises (−0.8 %, −2.5 %): at L/a = 20 shear
deformation and rotary inertia, which the solid elements carry and the Euler–Bernoulli beam
omits, lower the higher modes. The mode shapes are mass-normalised (φᵀ M φ = 1), so their
amplitude has no physical size; the pictures use deformation scales of 0.4 and 0.15.

## Verification

| Case | Reference | Error | Tolerance | Test |
|---|---|---|---|---|
| Cantilever, first and second bending frequencies, hex20 | theory (Euler–Bernoulli, Timoshenko, Young & Weaver 1974, §5.3) | −0.16 % (both first modes), −1.08 % (second) | 1 %, 2 % | `cantilever-frequencies` |
| Free–free beam: six rigid-body modes and the first bending mode; total mass | theory | rigid modes ≈ 0; −0.80 %; mass < 10⁻⁴ % | 1.5 %, 10⁻¹² | `free-beam-frequencies` |
| NAFEMS FV32 tapered membrane, six modes (hex20 layer) | benchmark | ≤ 0.0036 % | 1 % | `nafems-fv32` |
| Cantilever, hex20 consistent and lumped, tet10 | theory | first pair < 1 %, second pair < 3 %; mass exact | as stated | `cantilever_natural_frequencies` |
| Free–free beam | theory | rigid modes < 10⁻³ f₁; bending < 1.5 %; effective mass of the rigid translations = total mass to 10⁻⁶ | as stated | `free_free_beam_has_six_rigid_body_modes` |
| Simply supported plate (MITC4, MITC3), first frequency | theory (Navier) | < 1 % at 16 × 16 cells | 1 % | `simply_supported_plate_deflection_and_frequency` |
| Slender pipe cantilever, beam elements, three modes | theory (Euler–Bernoulli) | < 5·10⁻⁴ consistent, < 5·10⁻³ lumped; mass to 10⁻¹² | as stated | `cantilever_beam_frequencies` |
| Spring–mass, skew springs, two-mass chain, bar with a tip mass, torsional pendulum | theory | 10⁻¹² to 10⁻⁸ | as stated | `crates/solver/tests/one_d.rs` |
| 2D total mass (area × thickness, Pappus) and three in-plane rigid modes | theory | 10⁻¹² | exact | `masses_and_rigid_modes` |
| Repeated and zero eigenvalues; dense symmetric eigenproblem | theory | — | — | `eigen::tests::lanczos_finds_repeated_and_zero_eigenvalues`, `dense_symmetric_eigen` |
| Two runs of the same model | determinism | identical | exact | `results_are_deterministic` |

## Limitations

- Undamped, real modes only; no damping or complex eigenvalues. Equilibrium-base
  perturbations are limited to the explicitly supported element/tangent combinations.
- `shift` is factored by Cholesky: a shift above the lowest eigenvalue makes K − σM indefinite and
  the solve stops with the singular-matrix error naming a degree of freedom.
- Contact is refused; dashpots are ignored.
- Rotary inertia and translation–rotation mass coupling are not part of `total_mass`.

## Implementation

| File | Contents |
|---|---|
| `crates/engine/src/cmd/solve.rs` | `FREQUENCY` spec, `run_frequency` |
| `crates/solver/src/analysis.rs` | `solve_frequency`, `FrequencyOptions`, `MassKind`, `Mode`, `FrequencyResult`, `BLOCK`, `lanczos_basis`, `RIGID_SHIFT`, `TOLERANCE` |
| `crates/solver/src/eigen.rs` | `symmetric_eigen`, `lanczos`, `Operator`, `Wanted`, the xorshift64* generator |
| `crates/fem/src/element.rs`, `shell.rs`, `beam.rs`, `truss.rs` | consistent and lumped element masses |
| `crates/solver/src/discrete.rs` | point masses and rotary inertia |

## References

- T. Ericsson, A. Ruhe. *The spectral transformation Lanczos method for the numerical solution of large sparse generalized symmetric eigenvalue problems.* Math. Comp. 35(152):1251–1268, 1980.
- R. G. Grimes, J. G. Lewis, H. D. Simon. *A shifted block Lanczos algorithm for solving sparse symmetric generalized eigenproblems.* SIAM J. Matrix Anal. Appl. 15(1):228–272, 1994.
- B. N. Parlett. *The Symmetric Eigenvalue Problem.* SIAM, 1998.
- G. H. Golub, C. F. Van Loan. *Matrix Computations*, 4th ed. Johns Hopkins, 2013.
- S. Vigna. *An experimental exploration of Marsaglia's xorshift generators, scrambled.* ACM TOMS 42(4), 2016.
- E. Hinton, T. Rock, O. C. Zienkiewicz. *A note on mass lumping and related processes in the finite element method.* Earthquake Eng. Struct. Dyn. 4(3):245–249, 1976.
- S. Timoshenko, D. H. Young, W. Weaver. *Vibration Problems in Engineering*, 4th ed. Wiley, 1974.
- S. P. Timoshenko, J. N. Goodier. *Theory of Elasticity*, 3rd ed. McGraw-Hill, 1970.
- NAFEMS. *The Standard NAFEMS Benchmarks*, TNSB Rev. 3, October 1990.

## Axisymmetric ring shells

The two-ring axisymmetric formulation and its analysis/export restrictions are
documented in [Axisymmetric meridional shell](../../solvers/axisymmetric-shell.md).
Use the explicit `axisymmetric_shell` section flag; `beam2` alone remains a beam.
### Explicit section policies for a lumped solve

The command contract remains `mass=consistent|lumped`, default `consistent`. A
solid or planar-continuum section can set
`formulation.lumped_mass` to `hrz` (default), `row_sum`, or `consistent`, through
`model.section.create` / `model.section.edit`. This policy is active **only when
the solve explicitly selects `mass=lumped`**. A consistent solve ignores it.
The frequency reply reports the requested mass and section policies. The policy
is stored with the model; the global solve option remains a per-call choice.
Reopening a project does not make a later omitted `mass` argument lumped.

For `row_sum`, the diagonal is `L_ii = sum_j C_ij`, where `C` is the native
consistent element mass at its existing mass quadrature. No stiffness, density,
load or integration point for stiffness is changed. Supported positive-basis
linear continuum types are Tet4, Wedge6, Hex8, and planar Tri3/Quad4 (including
axisymmetric radius weighting). Higher-order sign-changing bases, pyramids,
shells and beams are refused for this override. Every nodal value must be finite
and strictly positive. `hrz` keeps the existing consistent-diagonal scaling and
its existing applicability. `consistent` retains the full consistent matrix of
a continuum section even within a globally lumped solve; this allows a linear
row-sum section and a quadratic consistent section in one model. Split sections
when their element types need different policies; one section has one policy.

An analytic tapered Hex8 (`x=(1+z)xi`, `y=eta`, unit parameter cube, density 1)
has mass 3/2. Row sums are 1/6 at each bottom node and 5/24 at each top node;
HRZ gives 5/32 and 7/32. A mixed Hex8/Tet4 assembly test also checks the tetrahedral
consistent diagonal `V/10`, off-diagonal `V/20`, and unchanged stiffness entries.
These are formula and implementation checks, not a claim of another solver's
mass defaults or stiffness equivalence.

Harmonic and native modal random commands accept the shared `mass` choice and
apply these section policies to inertia and material mass-proportional damping.
Default remains consistent. Projection of stored harmonic frames inherits that
result's mass choice and rejects an override; the project retains section policies.

## Convergence diagnostics

`solve.frequency` returns `convergence.iterations`, the actual number of block
recurrences performed. `lanczos_vectors` / `convergence.basis_vectors` count basis
vectors; these are not iteration counts. The legacy parameter `max_iterations`
remains a **basis-vector capacity**, not a count of block recurrences. A final
partial block counts as one iteration. `basis_limit` is the effective capacity
(user setting, memory budget and equation count, with at least one start block); `basis_limit_reached` reports
whether that capacity was reached, and `full_space_reached` distinguishes a complete
basis from a truncated user/memory cap. A truncated cap alone never certifies
convergence. An unconverged cap returns an error with the actual block iterations,
vector count and `basis_limit_reached=true`, and stores no replacement result.

Each returned mode has `equation_residual`:

- `residual_norm` = ‖Kφ − λMφ‖₂;
- `stiffness_action_norm` = ‖Kφ‖₂;
- `relative_residual` = their quotient, or `null` with `undefined_reason` when the
  denominator is zero (or the quotient overflows).

These use the actual **unshifted, constrained** stiffness and the selected
consistent/lumped mass, with the returned mode in the reduced equations. With a
base state, K is its actual terminal tangent, including its stress stiffness and
frozen constraints. Held-equation reactions are not counted as residual equations.
The Euclidean norms are evaluated in the solver's SI equation coordinates; they
are independent of mode normalization, but not of arbitrary equation rescaling.
This is not the projected Lanczos residual and does not replace it in the algorithm.
A small transformed Ritz bound does not by itself certify a small physical ratio.

For an exact rigid mode Kφ=0, the requested ratio is undefined, including 0/0: it
is never reported as zero. Near a rigid mode, roundoff can make ‖Kφ‖ nonzero;
`near_zero_stiffness_action` warns when it is at most √ε ‖K‖∞ ‖φ‖₂. A finite ratio
in that region is ill-conditioned, not a reliable rigid-mode acceptance criterion.
The threshold is a diagnostic warning only; it does not suppress a finite ratio.

Diagnostics are returned by this solve command, and by the Rust API
`solve_frequency_report`, separately from the numeric `FrequencyResult`. Save the
command's JSON reply when retaining this evidence. The unchanged native results
payload does **not** persist this new report; opening an old solved model cannot
supply historical iteration counts or residuals. Computing the report retains K
through Lanczos and performs two sparse products per mode inside the complete
solve command; its time and memory belong in performance comparisons.

For an opt-in orthogonality check, set `FIELDCRAFT_MODAL_ORTHOGONALITY=1`.
The solve evaluates the returned modes against its actual reduced mass matrix and
writes `FIELDCRAFT_MASS_ORTHOGONALITY` to stderr: mode count and
max|ΦᵀMΦ−I|. With solver timing enabled, `frequency_mass_orthogonality` reports
this additional work, which is included in the complete solve command time.
This diagnostic neither changes the modes nor adds fields to the stored result.

`FIELDCRAFT_PROJECTION_WORK` records successful projection stages when solver
timing is enabled: implementation, equations, basis columns, candidate columns,
passes, vector pairings, and dot-product elements. The subtraction traverses the
same number of vector elements. These are logical work counts, not cache or
DRAM counters; timings include workspace preparation and worker scheduling.

The block projection uses two reorthogonalization passes, fixed eight-column
panels and fixed 512-row partial dot products combined in row-block order.
Independent row blocks have one writer each. The worker count does not choose
the partition or the reduction order. This can round differently from scalar
Gram–Schmidt; thread-count reproducibility and agreement with an older
algorithm within its physical-result tolerance are separate checks. Repeated
eigenvalues require a subspace comparison, since individual vectors within a
repeated eigenspace are not unique.

`FIELDCRAFT_PROJECTION=scalar` selects the scalar reference for diagnostics in
the same executable. Both paths retain the convergence tests and random
replacement of vanished directions. Block projection additionally stores
packed basis and mass-basis panels; account for this workspace when comparing
peak process memory. This selection does not change the archive format.

## Curved-shell equilibrium bases

Curved MITC8 shells support small-preload linear static bases and
optional conservative pressure load stiffness through
`perturbation_follower_pressure` on the frequency step. See the
[formula, cylinder boundary conditions and refusals](../../solvers/curved-shell-modes.md).

Curved MITC6-a prestressed frequency remains refused: the new cylinder study
failed coarse-preload accuracy and fine-mesh actual-residual checks. Its existing
ordinary modes and verified buckling support are unchanged.
