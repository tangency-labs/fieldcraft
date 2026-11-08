# Steady-state harmonic response

`solve.harmonic` solves a linear structure's response to sinusoidal loads at specified
frequencies. Create a `harmonic` step, activate its loads and supports, then solve it by
name. The command is available to scripts, the CLI and MCP; `solve.run` also dispatches it.

## Parameters

```json
{
  "name": "Sweep",
  "procedure": {
    "type": "harmonic",
    "frequencies": {"start": "0 Hz", "end": "200 Hz", "count": 41},
    "method": {"type": "direct"},
    "damping": {"alpha": 0.2, "beta": 0.00001, "loss_factor": 0.01},
    "load_factors": {"Force": {"amplitude": 1, "phase": "30 deg"}}
  },
  "boundary_conditions": ["Clamp"],
  "loads": ["Force"],
  "output": ["displacement", "stress"]
}
```

Pass this object to `model.step.create`. Then call `solve.harmonic` with
`{"step":"Sweep"}` and `solve.probe` with `{"step":"Sweep","nodes":[42]}`.
Names must refer to definitions in the model. Each load keeps its existing spatial
distribution; its factor multiplies that distribution. A missing factor means unit amplitude
and zero phase. Factors may only name loads active in the step.

| Parameter | Meaning |
|---|---|
| `frequencies` | An ascending list, or `{start,end,count,logarithmic}`; 1–10,000 finite nonnegative points. Logarithmic ranges require positive endpoints. Frequencies are in the display frequency unit. |
| `method` | `{"type":"direct"}` or `{"type":"modal","modes":40,"damping_ratios":[0.02]}`. |
| `damping.alpha` | Mass-proportional viscous coefficient, inverse time. |
| `damping.beta` | Stiffness-proportional viscous coefficient, time. |
| `damping.loss_factor` | Dimensionless structural loss factor, not a damping ratio. |
| `damping_ratios` | Modal method only: empty, one ratio for every mode, or one per requested mode. |
| `method.residual_vectors` | Optional 0–16 load-dependent static/Ritz cycles per real and imaginary load component; default 0. |
| `load_factors` | Nonnegative dimensionless amplitude and phase in the display angle unit. Stored phases are radians. |
| `threads` | A solve parameter, using the same deterministic CPU thread pool as the other solvers. |

Material definitions accept the same `damping` object. Material and step coefficients are
**additive**. All coefficients and ratios must be finite and nonnegative. Consistent mass
is used, including point masses and rotary inertia. Material density is required where an
element carries inertia. Dashpots add their assembled viscous matrix.

Oriented solids use the same material-axis stiffness and stress recovery as static
analysis, including material damping. Laminated shells support analysis-step damping;
nonzero damping on individual ply materials is explicitly refused because ply-weighted
damping assembly is not yet implemented.

## Equations and damping

Physical quantities follow `u(t) = Re(U exp(i omega t))`, with `omega = 2 pi f`:

```text
[K - omega² M + i (omega C + G)] U = F
C = alpha_step M + beta_step K + sum(alpha_e M_e + beta_e K_e) + C_dashpot
G = eta_step K + sum(eta_e K_e)
F_load = amplitude (cos(phase) + i sin(phase)) F_spatial
```

`G` remains at zero frequency. Thus a constant structural loss gives complex compliance
there. To test the static limit, set all loss factors to zero; viscous damping then vanishes
as frequency tends to zero. No frequency offset, stiffness shift or damping floor is added
to a direct solve. An undamped resonance or an unsupported singular system returns an error.
Near a lightly damped pole, refine the frequency list and compare the full complex response:
small frequency or mesh errors can move the peak substantially.

The matrix is complex **symmetric**, not Hermitian. The direct solver uses native
complex supernodal LDLᵀ, reusing the existing sparse ordering and symbolic structure.
Local rook pivoting compares complex moduli and chooses one- or two-column blocks;
no positive-definiteness assumption is made near resonance. A real-only operator uses
the existing real factor. The interleaved real symmetric formulation
`[[A,-B],[-B,-A]] [Ur,Ui] = [Fr,-Fi]` remains a matrix-level verification path.

A native complex factor stores two real slots per original entry, versus approximately
four for the doubled real system. `max_factor_entries` counts 8-byte real slots, including
both components of every complex value. It defaults to 400 million and can be explicitly
raised to 4 billion for harmonic solves (32 GB of factor values alone). Indices, matrices
and workspace are additional. The symbolic cap is checked before numeric allocation.
Bounded iterative refinement checks the original complex equations. A failed local pivot
is reported; no pivot is delayed into a later supernode.

The modal method obtains real mass-normalized modes from the existing block Lanczos solver.
It projects material damping in full: `Phiᵀ C Phi` and `Phiᵀ G Phi`. Step Rayleigh damping
and loss are diagonal in this basis. The additional modal viscous coefficient is
`2 zeta_j omega_j`. Proportional damping uses scalar complex divisions; nonproportional
damping retains the off-diagonal terms and solves the small complex system. Increasing the
mode count is a necessary convergence check. A small projected residual alone does not
certify the untruncated physical equations.

For loads poorly represented by the retained modes, `residual_vectors` augments the basis
with static load responses and successive `K⁻¹ M` responses, mass-orthogonalized against
the retained basis. This follows Wilson, Yuan and Dickens (1982), DOI
[10.1002/eqe.4290100606](https://doi.org/10.1002/eqe.4290100606), and the force-truncation
construction in [SAND2025-03995, §2.22](https://www.sandia.gov/app/uploads/sites/315/2025/07/SD_Theory_5_24.pdf).
The recurrence and mass orthogonalization are also described in
[Abdallah (1990), §2.3, pp. 10–12](https://ntrs.nasa.gov/api/citations/19900013497/downloads/19900013497.pdf).
It reuses the extraction's static factor. Positive-definite unshifted stiffness and
homogeneous supports are required. Zero or dependent load vectors add no vectors.
The augmented operators retain their full stiffness, mass and damping couplings;
modal damping ratios still apply only to the requested physical modes. The reply reports
the actual `residual_vectors` count separately from `natural_frequencies`.
Increase the basis and check the full residual and response convergence; augmentation
does not guarantee accuracy near an omitted resonance. Deck export refuses this option
by name because an equivalent residual-vector construction is not implemented there.

Free-body mode extraction uses a spectral shift, which does not supply a physical
static stiffness. A zero-frequency modal solve is refused if extraction required
that shift; use the direct method to check the original system. For an undamped
diagonal mode, a forcing frequency indistinguishable from its pole at the relative
pivot tolerance is also refused. Finite neighboring frequencies and damped poles
retain their specified values.

## Output

The solve summary gives natural frequencies for modal runs, matrix sizes, warnings and,
per forcing frequency, `relative_residual`, `projected_residual` when applicable, input power
and dissipated power. Residuals use the infinity norm of the full reduced equations divided
by the force norm. With modal ratios, the full residual includes the reconstructed modal
damping force. Direct residuals above `1e-8` are refused; modal truncation residuals are
reported without concealing omitted response.

`residual_norm_si` also reports the algebraic Euclidean norm of the same full residual,
in the internal SI equation scaling. It is a numerical diagnostic and is not converted
with display units: a model with both translations and rotations mixes force and moment
equations. For a translational-only solid model, it is the force-imbalance norm in newtons.

`solve.probe` returns each node's frequency `samples`: displacement, rotation when present,
and requested stress components, each with `_real`, `_imaginary`, `_amplitude` and `_phase`.
An optional zero-based `frame` selects one frequency. Component amplitudes are nonnegative;
phase is `atan2(imaginary, real)`, with zero for exactly zero amplitude. These component
amplitudes do not define a von Mises stress amplitude. Units convert at the command boundary,
including phase, frequency and power. Probe CSV exports one row per node and frequency.
Saved `.fieldcraft` results retain the SI phasors.

For homogeneous supports the mean powers are
`P_in = omega/2 Im(F dot conj(U))` and
`P_loss = omega/2 U* (omega C + G) U`, including modal damping when requested.
For nonzero real prescribed displacements the reported powers refer to the reduced effective
force and response; they exclude a separate accounting of support work.

## Scope and verification

This is an undeformed linear elastic response: no preload stress stiffening, contact, thermal
harmonic loads or nonlinear material evolution. Stress is elastic nodal stress for solids
and rods, and the existing shell top-surface nodal recovery. Beam section stress, harmonic
contours/animation and VTK export are currently refused by name; use displacement-only output
for beams and probe/CSV for curves. Prescribed support displacements are real phasors; complex
base motion is not a separate input type.

The Abaqus writer emits direct or SIM modal steady-state steps, material/global damping,
modal ratios and separate real/imaginary loads, from public keyword documentation. A modal
step with nonzero prescribed motion, or a harmonic step sharing a nonlinear base state, is
refused. Nastran and LS-DYNA harmonic procedures are refused with the step name; these writers
do not silently export a static substitute.

Verification recipes are `harmonic-sdof`, `harmonic-reciprocity`, `harmonic-cantilever` and
`harmonic-nafems-5h`. Their source values and tolerances are declared in `expected.json`;
execution and fixture status appear in [the oracle report](../../oracle.md). Unit and model
tests also cover nonproportional projection, prescribed displacement, hostile input, the
static limit and deterministic thread counts. Published reference values do not imply that
a commercial solver was run.

## References and implementation

- J. K. Vandiver, MIT 2.003SC (2011), [lectures 24–25](https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-25-modal-analysis-response-to-ics-and-to-harmonic-forces/): harmonic forcing and modal response.
- J. R. Bunch (1982), [Stable algorithms for solving symmetric and skew-symmetric systems](https://doi.org/10.1017/S0004972700005621), §4: complex symmetric diagonal pivoting. C. Ashcraft, R. G. Grimes and J. G. Lewis (1998), [Accurate symmetric indefinite linear equation solvers](https://doi.org/10.1137/S0895479896296921): bounded rook search.
- D. Day and M. A. Heroux (2001), [Solving complex-valued linear systems via equivalent real formulations](https://doi.org/10.1137/S1064827500372262), SIAM J. Sci. Comput. 23(2), 480–498. The symmetric block identity above is derived directly from the real and imaginary equations.
- R. R. Craig and A. J. Kurdila, *Fundamentals of Structural Dynamics*, 2nd ed., 2006: real-mode projection and nonproportional damping.
- [Abaqus public damping reference](https://docs.software.vt.edu/abaqusv2026/English/SIMACAEMATRefMap/simamat-c-dampingopt.htm) and [steady-state keyword](https://docs.software.vt.edu/abaqusv2026/English/SIMACAEKEYRefMap/simakey-r-steadystatedynamics.htm): independent coefficient and export definitions.

The model contract is in `crates/model/src/harmonic.rs`; the matrix kernels in
`crates/solver/src/{complex,harmonic}.rs`; model assembly and recovery in
`crates/solver/src/harmonic/model.rs`; commands in `crates/engine/src/cmd/solve.rs`.

## Projecting an input PSD onto stored frames

After a unit-amplitude spatial load has been solved, `solve.random` accepts the stored
harmonic step and a one-sided PSD table in SI units, per Hz:

```json
{"command":"solve.random","params":{"step":"Sweep","input_psd":[[0,1],[2000,1]]}}
```

The table must cover every stored frequency. Linear interpolation supplies the input PSD;
complex transfer magnitudes squared give response PSD, and trapezoidal frequency integration
gives RMS. Returned nodal displacement spectra sum the three translational component spectra.
The optional stress spectrum sums six component spectra: it is neither an individual normal
stress spectrum nor a von Mises statistic. Recover individual components from the complex
harmonic stress frames when those are required. The returned arrays follow the stored result's
node order. This path does not assemble a general random-vibration analysis step or automatically
interpret a base-motion excitation.

Without `base_direction`, the acceleration calculation is `omega^4` times the displacement PSD. It represents absolute
acceleration for a force-loaded structure with stationary supports. If harmonic frames instead
represent relative displacement under a unit base-acceleration inertial load `-M r`, that output
is relative acceleration. Pass `"base_direction":[rx,ry,rz]` (a unit vector) when the stored frames are the response to
the matching unit inertial load `-M r`. The command forms absolute acceleration coherently
from `r - omega^2 H_relative`, before squaring and multiplying by the base PSD. Adding the
base and relative acceleration PSDs would lose their cross term. A white acceleration spectrum
extending to zero frequency does not define a finite stationary absolute base-displacement RMS;
label relative displacement explicitly.

For uniform translational support motion, an existing whole-model gravity load with acceleration
opposite the unit base direction supplies the consistent inertial spatial load. This relative-motion
construction requires supports moving together and damping formulated for relative motion; it
does not cover independent support motions. Modal truncation, residual correction, frequency
sampling and stress recovery still require checks for the particular model.

`nodes:[label,...]` selects solved nodes in the requested order (no duplicates). `probes`
contains per-node `samples` with frequency, displacement/acceleration component PSDs and,
when recovered by harmonic, six stress component PSDs (xx,yy,zz,xy,xz,yz). Per-component RMS
and vector RMS are both reported. Returned random quantities are explicitly SI, irrespective
of display units; the input PSD must match the unit spatial harmonic load's SI amplitude.
Base motion returns relative displacement and absolute acceleration, labeled in the result.
The user must supply the matching inertial spatial load: the command does not infer it from
a PSD or change the model. At most four million selected node-frame samples are allowed.

Executable command oracles in `tests/oracle/random-*` cover exact finite-band SDOF integrals,
Miles' approximation separately, low/narrow-band limits, coherent cantilever modal response,
and a source-bound 5R known discrepancy. They do not certify arbitrary load normalization,
independent multiple inputs, or a general `RandomVibration` step assembler.
