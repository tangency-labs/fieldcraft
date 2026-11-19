# RFC 0001: High-frequency electromagnetics

- **Status:** proposed, open for comment
- **Discussion:** [#971](https://github.com/tangency-labs/fieldcraft/issues/971)
- **Scope:** time-harmonic Maxwell (frequency domain) on tetrahedral meshes: cavity eigenmodes,
  driven problems with waveguide and lumped ports, S-parameters, radiation and scattering
- **Out of scope here:** time-domain EM (FDTD, DGTD), low-frequency magnetostatics and eddy
  currents (ROADMAP M19 "low-frequency electromagnetics", a later RFC), fast integral-equation
  methods, and coupled EM–thermal–structural analysis (once this lands, it becomes a coupling
  through the existing heat solver)

This RFC proposes how electromagnetics fits FieldCraft's existing layers, what it reuses, what it
adds, and how it is verified. It is a plan to comment on; the decisions marked **Open** are the
owner's.

## 1. What we solve

The vector wave equation for the electric field **E**, at angular frequency ω:

  ∇ × (μᵣ⁻¹ ∇ × **E**) − k₀² εᵣ **E** = −j k₀ Z₀ **J**,  k₀ = ω √(μ₀ε₀)

with complex, possibly anisotropic εᵣ and μᵣ (losses, PML). The weak form in H(curl) gives

- **eigenmodes** (lossless cavity): K e = k² M e, K and M real symmetric, M positive definite,
  K positive semi-definite;
- **driven problems**: (K − k₀² M + j k₀ B) e = f, complex **symmetric** (A = Aᵀ, not Hermitian)
  when the materials are reciprocal.

Theory: J.-M. Jin, *The Finite Element Method in Electromagnetics*, 3rd ed., Wiley 2014; P. Monk,
*Finite Element Methods for Maxwell's Equations*, Oxford 2003.

## 2. What exists today, and what it means for EM

| Piece | Today | For EM |
|---|---|---|
| Unknowns (`fem::dof`) | Nodal only: a `DofKind` per node, numbered node by node in label order | Needs unknowns on **edges** (and on faces at second order) |
| Matrix (`solver::sparse::Csr`) | `f64` values | Needs complex values |
| Direct solvers (`solver::cholesky`, `solver::ldlt`) | Supernodal Cholesky and LDLᵀ with bounded Bunch–Kaufman pivoting, multithreaded, fixed operation order | Reused: LDLᵀ for indefinite real systems as is; its numeric kernels made generic for complex-symmetric |
| Ordering (`solver::ordering`) | Geometric nested dissection on a coordinate per unknown | Reused unchanged: an edge unknown takes its edge midpoint |
| Eigen (`solver::eigen`) | Shifted block Lanczos, K φ = λ M φ, M inner product, full reorthogonalization | Reused unchanged for cavity modes and for port modes |
| Constraints (`solver::constraints`) | Prescribed values and linear constraint rows, eliminated from the system | PEC (tangential **E** = 0) is a prescribed value on boundary edges; periodic boundaries are constraint rows |
| Geometry map (`fem::basis`) | Isoparametric map, Jacobians on every tet family | Reused for the covariant Piola map of edge elements |
| Results, jobs, cancellation, determinism (`fieldcraft-parallel`) | In place | Reused |
| UI results | Contours, mode animation, path plot | Complex fields fit contours and animation; S-parameters, Smith chart and patterns need new panels |

No iterative solver exists. Direct factorization covers small and medium models (a few million
unknowns within the memory caps the solvers already enforce).

## 3. Proposal

### 3.1 Edge degrees of freedom

Extend `DofLayout` so a degree of freedom belongs to a **mesh entity**: a node (today), an edge or
a face. Edges and faces are numbered after the nodes, in a fixed order derived from sorted vertex
labels, so the same mesh always gives the same numbering (the determinism rule).

- **Orientation:** an edge runs from its lower to its higher global vertex label; an element's
  local edge takes sign ±1 against it. Face DOFs of second-order elements use the same rule on the
  sorted vertex triple. Sign conventions are the main source of silent errors in edge elements, so
  they get their own unit tests (a constant curl field reproduced exactly, the tangential
  continuity across shared faces).
- **Elements:** Nédélec first kind on tetrahedra (J.-C. Nédélec, *Mixed finite elements in ℝ³*,
  Numer. Math. 35:315–341, 1980): first order (6 edge unknowns, Whitney forms) and second order
  (hierarchical, 20 unknowns; J. P. Webb, *Hierarchal vector basis functions of arbitrary order
  for triangular and tetrahedral finite elements*, IEEE Trans. Antennas Propag. 47(8):1244–1253,
  1999). Curved (quadratic-geometry) tets through the covariant Piola map.

**Open:** a separate `EdgeLayout` beside `DofLayout`, or one layout keyed by entity. The proposal
is one layout, because constraints, results and the control channel then see a single numbering.

### 3.2 Complex arithmetic

Make `Csr`, assembly and the LDLᵀ numeric phase generic over a scalar trait implemented for `f64`
and a complex type. The symbolic phase (elimination tree, supernodes, column counts) depends only
on the pattern and is shared unchanged.

- Complex-symmetric factorization uses the transpose, not the conjugate transpose:
  A = L D Lᵀ with the same bounded Bunch–Kaufman pivot search on |aᵢⱼ|.
- The `f64` path must stay bit-for-bit identical: the oracle fixtures prove it.
- **Open:** a hand-written complex type in `fieldcraft-fem` or the `num-complex` crate
  (MIT/Apache-2.0). The proposal is `num-complex`, recorded in docs/adr/0001-libraries.md.

### 3.3 Eigenmodes

A cavity step reuses `solver::eigen` as it is. The curl–curl operator has a large null space
(every gradient field has k² = 0). With a shift σ between 0 and the first cutoff, K − σM is
indefinite, so the step factors it with LDLᵀ. The null-space modes map to θ = −1/σ, away from the
wanted largest θ. The LDLᵀ inertia (a Sturm count) confirms how many modes lie below a frequency.

**Open:** whether to also remove the null space explicitly (a discrete-gradient projection,
Monk 2003 §13) for robustness on large models.

### 3.4 Ports, boundaries and sources

- **PEC and PMC:** prescribed tangential **E** (PEC) and the natural condition (PMC).
- **Waveguide ports:** a 2D eigenproblem on the port face (edge elements on triangles, the same
  Lanczos) gives the modal fields and propagation constants. A modal boundary condition excites
  and absorbs them (Jin 2014, ch. 8).
- **Lumped ports:** a gap source with a reference impedance.
- **Open boundaries:** a first-order absorbing boundary condition, then a PML as an anisotropic
  material (Z. S. Sacks, D. M. Kingsland, R. Lee, J.-F. Lee, *A perfectly matched anisotropic
  absorber for use as an absorbing boundary condition*, IEEE Trans. Antennas Propag.
  43(12):1460–1463, 1995).
- **Periodic (Floquet) boundaries:** constraint rows with a phase, for unit cells.

### 3.5 Post-processing (commands first, panels second)

Every view reads its numbers from a command, so scripts, MCP and exports get exactly what the UI
draws.

- `results.sparameters` — S-matrix per frequency, reference impedances, Touchstone export
- `results.far_field` — near-to-far-field transformation on a closed surface (equivalence
  principle; C. A. Balanis, *Antenna Theory*, 4th ed., Wiley 2016, ch. 12): gain, directivity,
  RCS, polarization
- `results.port_modes` — cutoff, propagation constant and impedance per mode
- Complex nodal or element fields: magnitude, real and imaginary parts, phase at a chosen angle

UI, on SolidCraft's UI kit like every panel:
- complex contours and an ωt animation through the existing contour and animation path;
- an S-parameter chart (magnitude in dB, phase);
- a Smith chart;
- radiation patterns: polar cuts and a 3D pattern in the viewport.

### 3.6 Commands and model

The usual pattern, each a `CommandSpec` with tests:

- `model.material.create` gains εᵣ, μᵣ and σ (with loss tangent)
- `model.boundary_condition.create` gains `pec`, `pmc`, `absorbing`, `periodic`
- `model.port.create` (waveguide and lumped)
- `model.step.create` gains `em_eigen` and `em_driven` (a frequency list or sweep)
- `solve.em`

Units go through `crates/model/src/units` (SI stored, display units at the boundaries).

## 4. Verification

Each case is our own model, built from its published description, compared with theory, measured
data or a published benchmark, and frozen as an oracle fixture (docs/verification.md).

| Case | Reference | Checks |
|---|---|---|
| Rectangular cavity TE₁₀₁ and higher modes | D. M. Pozar, *Microwave Engineering*, 4th ed., Wiley 2012, §6.4 | Resonant frequencies; no spurious modes; convergence order of first- vs second-order elements |
| Partially dielectric-filled cavity | Jin 2014 | Frequencies from the transcendental equation |
| WR-90 waveguide section | Pozar 2012 §3.3 | β of TE₁₀; \|S₂₁\| = 1, S₁₁ ≈ 0 for a matched line |
| Dielectric slab in a waveguide | Pozar 2012, ch. 2–3 | S₁₁ and S₂₁ against the closed-form transmission-line result |
| Coaxial line | Pozar 2012 §3.5 | Characteristic impedance |
| PEC sphere scattering | Mie series (J. A. Stratton, *Electromagnetic Theory*, McGraw-Hill 1941, §9.25) | Bistatic RCS; PML reflection |
| Half-wave dipole | Balanis 2016 §4.6 | Pattern, directivity 1.64, input impedance ≈ 73 + j42.5 Ω |

Cheap cases run in the trunk gate; sweeps and convergence studies are `heavy`.

## 5. Plan

Each phase lands as reviewed pull requests with its tests and its verification cases.

| Phase | Work | Proposed owner |
|---|---|---|
| 0 | This RFC agreed; open questions decided | Maintainer and EM contributor |
| 1 | Entity-keyed `DofLayout` (edges, faces); scalar-generic `Csr`, assembly and LDLᵀ; `f64` path bit-identical | Core |
| 2 | Nédélec elements (order 1, then 2); PEC; cavity eigen step; cavity cases | EM contributor |
| 3 | Complex driven step; waveguide and lumped ports; S-parameters; waveguide cases | EM contributor |
| 4 | Absorbing boundary and PML; far field; sphere and dipole cases | EM contributor |
| 5 | UI: complex contours and animation, S-parameter chart, Smith chart, radiation patterns | Core and EM contributor |
| Later | Periodic boundaries, iterative solver with an auxiliary-space preconditioner for large models (R. Hiptmair, J. Xu, *Nodal auxiliary space preconditioning in H(curl) and H(div) spaces*, SIAM J. Numer. Anal. 45(6):2483–2509, 2007), EM–thermal coupling | — |

## 6. Rules that apply (from AGENTS.md)

- Methods come from the cited papers and textbooks; no code from other EM codes is read, ported or
  adapted, and GPL/LGPL/AGPL code is never read.
- Results are verified against theory, measurement and published benchmarks, never against
  another solver's output.
- No panics; every iteration and size is bounded; input-derived numbers are checked.
- Deterministic results on CPU, bit for bit.
- Everything is a command; the UI, scripts, CLI and MCP reach the same commands.

## 7. Open questions

1. One entity-keyed `DofLayout`, or a separate edge layout (§3.1)?
2. `num-complex` or our own complex type (§3.2)?
3. Explicit null-space projection in the eigen step (§3.3)?
4. First-order only for the first release, or first and second order together?
5. Which first application drives the priorities: filters and waveguide components (ports and
   S-parameters first) or antennas (PML and far field first)?
