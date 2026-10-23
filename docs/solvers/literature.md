# Solver literature survey for FieldCraft: explicit dynamics, equations of state, lattice Boltzmann, FSI

Compiled 2026-10-13. Scope: three new solver families on top of the planned implicit linear solver (M6):
(A) an explicit Lagrangian finite element solver of the LS-DYNA / Abaqus-Explicit class, (B) equations of state for
that solver's shock and detonation work, (C) a lattice Boltzmann CFD solver of the Palabos class, (D) fluid–structure
interaction coupling them, and (E) the alternatives (finite volume, SPH) and GPU compute in pure Rust. The plan is at
the end (Part F).

## How this was verified

- **DOIs.** Every DOI below resolved through the DOI handle API (`doi.org/api/handles/<doi>` → responseCode 1),
  and its first author, year, title and venue were checked against Crossref metadata. Volume and page numbers are
  from memory and were not each checked; title, first author, year and venue were. Four DOIs I had from memory were
  wrong and were corrected from Crossref search before inclusion:
  - Flanagan–Belytschko 1984 eigenvalue paper is `10.1115/1.3167594` (`…598` is an unrelated crack paper).
  - Johnson–Stryk 1987 eroding interface is `10.1016/0734-743X(87)90057-1`.
  - He–Luo 1997 incompressible LBM is `10.1023/B:JOSS.0000015179.12689.e4` (`10.1007/BF02732106` is a plant-biology paper).
  - Belytschko–Lin 1987 erosion is `10.1016/0734-743X(87)90033-9`.
- **Other URLs** returned HTTP 200 on 2026-10-13 unless marked. Items without a DOI or working URL are marked
  **(unverified)**.
- **Licences** come from GitHub/GitLab licence metadata (`api.github.com/repos/<r>`, GitLab `?license=true`), the
  repository's LICENSE/COPYING file header, or a README licence sentence when the API reported nothing. **No source
  code was opened.** Where the API said `NOASSERTION`, the first lines of the licence file decided it.
- No section numbers are given: I did not read the PDFs for this survey, so I give none rather than guess.

### Licence legend for associated code

| Status | Meaning |
|---|---|
| **READ OK** | MIT, BSD, Apache-2.0, ISC, MPL-2.0, public domain. May be read for reference; still re-implement from the paper and record anything adapted in ATTRIBUTION.md. |
| **FORBIDDEN** | GPL, LGPL or AGPL, or a non-commercial / research-only licence. Do not open. |
| **NO LICENCE** | No licence found in metadata: all rights reserved. Do not open. |

### Associated code bases (licence from metadata only)

| Code | What it is | Licence (source) | Status |
|---|---|---|---|
| Palabos `gitlab.com/unigespc/palabos` | LBM framework | AGPL-3.0 (GitLab API) | FORBIDDEN |
| OpenLB `gitlab.com/openlb/release` | LBM framework | GPL-2.0+ (GitLab API) | FORBIDDEN |
| waLBerla `i10git.cs.fau.de/walberla/walberla` | LBM framework | GPL-3.0 (GitLab API) | FORBIDDEN |
| lbmpy `i10git.cs.fau.de/pycodegen/lbmpy` | LBM code generator (cumulant, MRT) | GPL-3.0+ (GitLab API) | FORBIDDEN |
| VirtualFluids `git.rz.tu-bs.de/irmb/VirtualFluids` | cumulant LBM (Geier/Krafczyk group) | GPL-3.0 (GitLab API) | FORBIDDEN |
| TCLB `CFD-GO/TCLB` | GPU LBM | GPL-3.0 | FORBIDDEN |
| HemeLB `hemelb-codes/hemelb` | sparse-geometry LBM | LGPL-3.0 | FORBIDDEN |
| Sailfish `sailfish-team/sailfish` | GPU LBM | LGPL-3 (README; API: none) | FORBIDDEN |
| FluidX3D `ProjectPhysX/FluidX3D` | GPU LBM, esoteric pull | custom: "public research, education or personal use" (LICENSE.md) | FORBIDDEN (non-commercial) |
| ESPResSo `espressomd/espresso` | MD + LBM | GPL-3.0 | FORBIDDEN |
| OpenFOAM `OpenFOAM/OpenFOAM-dev` | FV CFD | GPL (COPYING) | FORBIDDEN |
| SU2 `su2code/SU2` | FV CFD | LGPL-2.1 (LICENSE.md) | FORBIDDEN |
| Fluidity | FE/FV CFD | LGPL-2.1 | FORBIDDEN |
| OpenRadioss | explicit FE | AGPL-3.0 (ADR 0001; GitHub mirrors report AGPL-3.0; the `OpenRadioss/OpenRadioss` repo URL returned 404 today) | FORBIDDEN |
| preCICE `precice/precice` | partitioned coupling, IQN, RBF mapping | LGPL-3.0 | FORBIDDEN |
| 4C (ex-BACI) `4C-multiphysics/4C` | FSI, Aitken (Küttler–Wall group) | LGPL-3.0 | FORBIDDEN |
| MOOSE `idaholab/moose` | multiphysics FE | LGPL-2.1 | FORBIDDEN |
| FEniCS `FEniCS/dolfinx` | FE library | LGPL-3.0 | FORBIDDEN |
| deal.II `dealii/dealii` | FE library | now Apache-2.0 WITH LLVM-exception **and** LGPL-2.1 (LICENSE.md) | FORBIDDEN by ADR 0001 (it lists deal.II as LGPL). The relicensing is an owner decision; until the ADR changes, don't open. |
| DualSPHysics `DualSPHysics/DualSPHysics` | GPU SPH | LGPL-2.1 | FORBIDDEN |
| IB2d `nickabattista/IB2d` | 2D immersed boundary | GPL-3.0 | FORBIDDEN |
| CalculiX, Code_Aster, Elmer | implicit FE | GPL family (ADR 0001) | FORBIDDEN |
| STLBM `gitlab.com/unigehpfs/stlbm` | Latt et al. 2021 portable C++ LBM (AA-pattern, BGK/TRT/regularized) | MIT (GitLab API) | **READ OK** — the one readable code by Palabos's authors |
| XLB `Autodesk/XLB` | differentiable LBM (JAX/Warp) | Apache-2.0 (LICENSE) | READ OK |
| lettuce `lettucecfd/lettuce` | PyTorch LBM | MIT | READ OK |
| pylbm `pylbm/pylbm` | LBM schemes generator | BSD-3 ("New BSD", LICENSE) | READ OK |
| Ludwig `ludwig-cf/ludwig` | LBM complex fluids | BSD-3-style (LICENSE: redistribution conditions) | READ OK |
| Laghos `CEED/Laghos` | high-order Lagrangian hydro, artificial viscosity | BSD-2-Clause | READ OK |
| MFEM `mfem/mfem` | FE library | BSD-3-Clause | READ OK |
| LULESH `LLNL/LULESH` | Lagrangian hydro proxy app (one-point hex + hourglass) | none in metadata, no licence file, README has no licence sentence | NO LICENCE (until a licence file is confirmed; do not open) |
| Kratos `KratosMultiphysics/Kratos` | multiphysics incl. FSI, IQN, mapping | BSD-style (kratos/license.txt: "Redistribution … permitted"; per-application licences may differ) | READ OK for the core; check each application folder's own licence file first |
| IBAMR `IBAMR/IBAMR` | immersed boundary (Peskin family) | BSD-3-style (COPYRIGHT file) | READ OK |
| MUI `MxUI/MUI` | multiscale coupling interface, RBF | Apache-2.0 (GitHub API; README says dual GPL/Apache, take Apache) | READ OK |
| FEBio `febiosoftware/FEBio` | nonlinear implicit FE, FSI | MIT | READ OK |
| WARP3D `rhdodds/warp3d` | nonlinear implicit FE | NCSA/MIT-style (LICENSE) | READ OK |
| SPlisHSPlasH `InteractiveComputerGraphics/SPlisHSPlasH` | SPH fluids | MIT | READ OK |
| SPHinXsys `Xiangyu-Hu/SPHinXsys` | SPH multiphysics incl. FSI | Apache-2.0 | READ OK |
| PySPH `pypr/pysph` | SPH | BSD-3-style (LICENSE.txt) | READ OK |
| Cabana `ECP-copa/Cabana` | particle toolkit | BSD-3-Clause | READ OK |
| Nek5000 | spectral-element CFD | BSD-3-style (LICENSE) | READ OK |
| fenris `InteractiveComputerGraphics/fenris` | Rust FE library | Apache-2.0 | READ OK |
| wgpu `gfx-rs/wgpu` | Rust WebGPU | MIT OR Apache-2.0 | READ OK, dependency candidate |
| rust-gpu `Rust-GPU/rust-gpu` | Rust → SPIR-V | MIT OR Apache-2.0 (API: Apache-2.0) | READ OK |
| faer `sarah-quinones/faer-rs` | Rust dense/sparse LA | MIT | READ OK |

---

## Part A: Explicit finite element dynamics

### A1. Time integration and stability

**A1.1 Belytschko, Liu, Moran, Elkhodary. *Nonlinear Finite Elements for Continua and Structures*, 2nd ed., Wiley 2014.** **Verdict: build on (the spine of the explicit solver).**
- **Links.** ISBN 978-1-118-63270-3; publisher page https://www.wiley.com/en-us/Nonlinear+Finite+Elements+for+Continua+and+Structures%2C+2nd+Edition-p-9781118632703 (200). No DOI checked.
- **Relevant chapters** (by topic; chapter numbers from memory of the 2nd edition, check before citing in module docs):
  updated Lagrangian formulation and the internal-force loop; explicit central difference with the
  "flowchart" for one step, energy balance check, critical time step and element eigenvalue bounds; constitutive
  models (hypoelastic–plastic, objective rates: Jaumann, Truesdell, Green–Naghdi); element technology (one-point
  quadrature, hourglass control, B-bar, assumed strain, shells including Belytschko–Lin–Tsay); contact-impact
  (penalty, Lagrange multipliers, explicit contact algorithms); stability (linearised stability, energy method).
- **Method.** Central difference in the half-step form: `v^{n+1/2} = v^{n-1/2} + Δt^n M⁻¹(f_ext − f_int)^n`,
  `x^{n+1} = x^n + Δt^{n+1/2} v^{n+1/2}`, with a diagonal (lumped) mass so each step is a gather of element
  internal forces and a nodal division. Stable for `Δt ≤ 2/ω_max` (with damping ratio ξ:
  `Δt ≤ (2/ω_max)(√(1+ξ²) − ξ)`). ω_max is bounded by the largest element eigenvalue (element eigenvalue inequality),
  giving `Δt = α · min_e L_e / c_e` with a safety factor α ≈ 0.9 and `c = √((λ+2μ)/ρ)` for solids.
- **Solves.** Everything in Part A; the book is the integration framework all the papers below plug into.
- **Limits.** Conditionally stable; cost ∝ (number of elements) × (end time / Δt), and Δt is set by the smallest
  or stiffest element anywhere in the model. One distorted sliver element sets the step for the whole model.
- **Complexity.** The core loop is small (≈1–2k lines); everything hard is in elements, materials and contact.
- **Code.** None.

**A1.2 Courant, Friedrichs, Lewy 1928.** **Verdict: reference (the CFL condition).**
- **Citation.** R. Courant, K. Friedrichs, H. Lewy. *Über die partiellen Differenzengleichungen der mathematischen Physik.* Math. Ann. 100:32–74, 1928. DOI [10.1007/BF01448839](https://doi.org/10.1007/BF01448839).
- **Use.** The origin of `Δt ≤ L/c`; cite for the explicit FE and LBM time-step rules.

**A1.3 Flanagan & Belytschko 1984, eigenvalues and stable time steps.** **Verdict: build on.**
- **Citation.** D. P. Flanagan, T. Belytschko. *Eigenvalues and Stable Time Steps for the Uniform Strain Hexahedron and Quadrilateral.* J. Appl. Mech. 51(1):35–40, 1984. DOI [10.1115/1.3167594](https://doi.org/10.1115/1.3167594).
- **Method.** Closed-form upper bound on the maximum eigenvalue of the one-point hexahedron/quadrilateral in terms of
  its gradient vectors `b_iI` and volume; it yields a characteristic length that is tighter and cheaper than
  "shortest edge" and is exact for a cube.
- **Solves.** The per-element critical time step for the element in A2.1. Tets and shells use their own lengths
  (tet: volume / largest face area times a factor; shell: area / longest diagonal or side).
- **Limits.** A bound per element, so it is conservative; nodal estimates (below) can be larger.
- **Complexity.** Low.
- **Note on nodal estimates.** Nodal time-step estimates (from nodal mass and an assembled nodal stiffness bound) give
  larger stable steps on graded meshes. I found no single verified primary paper to cite for them; treat as
  "reference later" and keep the element estimate as the guaranteed-stable default. **(unverified)**

**A1.4 Lumped mass.** **Verdict: build on.** Row-sum lumping for linear elements (hex8, tet4, quad4, tri3), which is
positive; for quadratic elements row-sum gives zero or negative corner masses, so use HRZ (diagonal scaling)
lumping or avoid quadratic elements in explicit (the usual choice). Source: A1.1 book; Hughes, *The Finite Element
Method* (Dover 2000, no DOI) for HRZ. **(book, unverified page)**

**A1.5 Olovsson, Simonsson, Unosson 2005, selective mass scaling.** **Verdict: build on (after plain mass scaling).**
- **Citation.** L. Olovsson, K. Simonsson, M. Unosson. *Selective mass scaling for explicit finite element analyses.* IJNME 63(10):1436–1445, 2005. DOI [10.1002/nme.1293](https://doi.org/10.1002/nme.1293).
- **Method.** Conventional mass scaling adds nodal mass to small elements until each reaches a target Δt (adds
  inertia to rigid-body modes too, so it changes the dynamics). Selective mass scaling adds a non-diagonal mass term
  that leaves translational (rigid-body) inertia unchanged and lowers only the high frequencies, at the cost of a cheap
  iterative solve per step.
- **Solves.** Quasi-static and thin-structure explicit runs whose step is limited by a few small elements.
- **Limits.** Conventional scaling must report the added mass (% of total) — silent mass scaling is a classic source
  of wrong answers. Selective scaling needs a solve each step (CG on a well-conditioned system).
- **Complexity.** Conventional: trivial. Selective: moderate.
- **Code.** None.

### A2. Solid element technology

**A2.1 Flanagan & Belytschko 1981, uniform strain hexahedron with orthogonal hourglass control.** **Verdict: build on (the default solid element).**
- **Citation.** D. P. Flanagan, T. Belytschko. *A uniform strain hexahedron and quadrilateral with orthogonal hourglass control.* IJNME 17(5):679–706, 1981. DOI [10.1002/nme.1620170504](https://doi.org/10.1002/nme.1620170504).
- **Method.** One-point integration with the mean (volume-averaged) gradient operator, so the element passes the patch
  test and costs one stress update per step. The four hourglass modes per direction are made orthogonal to the
  linear velocity field (the γ vectors), and a viscous or stiffness-type hourglass resisting force with a
  user coefficient (≈0.05–0.15) suppresses them without touching constant-strain modes.
- **Solves.** Fast, lock-free (under one-point integration) hex8 for impact and crash.
- **Limits.** Hourglass coefficient is a tuning knob; stiffness hourglass is too stiff in bending on coarse meshes,
  viscous hourglass lets hourglass energy grow in long runs. Monitor hourglass energy against internal energy (rule of
  thumb < 10%) and report it.
- **Complexity.** Moderate (≈500 lines for the element with both hourglass types).
- **Code.** LULESH implements this element but has NO LICENCE in its metadata: don't open. Laghos (BSD-2) is a
  high-order hydro and does not use this element.

**A2.2 Belytschko, Ong, Liu, Kennedy 1984, hourglass control in linear and nonlinear problems.** **Verdict: reference.**
- **Citation.** CMAME 43(3):251–276, 1984. DOI [10.1016/0045-7825(84)90067-7](https://doi.org/10.1016/0045-7825(84)90067-7).
- **Method.** Unifies hourglass control as stabilisation of the one-point element, gives the stabilisation
  parameters and their consistency for nonlinear problems.
- **Use.** Background for choosing coefficients; read with A2.1.

**A2.3 Belytschko & Bindeman 1993, assumed strain stabilisation of the hex8.** **Verdict: build on (second solid formulation).**
- **Citation.** T. Belytschko, L. P. Bindeman. *Assumed strain stabilization of the eight node hexahedral element.* CMAME 105(2):225–260, 1993. DOI [10.1016/0045-7825(93)90124-G](https://doi.org/10.1016/0045-7825(93)90124-G).
- **Method.** Derive the hourglass (stabilisation) forces from an assumed-strain field instead of an arbitrary
  coefficient, so the stabilisation stiffness is physical: correct bending with one element through the thickness,
  no shear locking, and an incompressible variant that avoids volumetric locking. Co-rotational frame.
- **Solves.** Coarse-mesh bending accuracy that A2.1 lacks, without a tuning coefficient.
- **Limits.** Physical stabilisation is elastic: it uses the elastic moduli (or a tangent estimate) in plasticity,
  so it is approximate in large plastic flow.
- **Complexity.** Moderate–high.

**A2.4 Puso 2000, enhanced assumed strain physically stabilised hexahedron.** **Verdict: reference (later upgrade over A2.3).**
- **Citation.** M. A. Puso. *A highly efficient enhanced assumed strain physically stabilized hexahedral element.* IJNME 49(8):1029–1064, 2000. DOI [10.1002/1097-0207(20001120)49:8<1029::AID-NME990>3.0.CO;2-3](https://doi.org/10.1002/1097-0207(20001120)49:8%3C1029::AID-NME990%3E3.0.CO;2-3).
- **Method.** Takes Taylor-expanded EAS (Simo–Rifai) stabilisation to one-point integration; better on distorted
  (trapezoidal) meshes than Belytschko–Bindeman, still ≈ one-point cost.
- **Limits.** More algebra; same elastic-stabilisation caveat.
- **Complexity.** High.

**A2.5 Simo & Rifai 1990, enhanced assumed strain.** **Verdict: reference (implicit solver first).**
- **Citation.** IJNME 29(8):1595–1638, 1990. DOI [10.1002/nme.1620290802](https://doi.org/10.1002/nme.1620290802).
- **Use.** EAS hex8 for the implicit solver (M6); background for A2.4. Static condensation per element is awkward in
  explicit.

**A2.6 Selective / reduced integration and B-bar.** **Verdict: build on (fully integrated hex8 option).**
- **Citations.**
  - D. S. Malkus, T. J. R. Hughes. *Mixed finite element methods — reduced and selective integration techniques: a unification of concepts.* CMAME 15(1):63–81, 1978. DOI [10.1016/0045-7825(78)90005-1](https://doi.org/10.1016/0045-7825(78)90005-1).
  - T. J. R. Hughes. *Generalization of selective integration procedures to anisotropic and nonlinear media.* IJNME 15(9):1413–1418, 1980. DOI [10.1002/nme.1620150914](https://doi.org/10.1002/nme.1620150914).
- **Method.** Integrate the deviatoric part at 2×2×2 points and the volumetric part at one point (B-bar: replace the
  dilatational part of B with its element mean). No hourglass modes, no volumetric locking.
- **Limits.** Eight stress updates per element (≈4–8× the cost of A2.1) and shear locking in bending.
- **Complexity.** Low given the implicit solver's hex8.
- **Shared with M6.** The implicit solver needs exactly this element; write it once in a shared crate (Part F).

**A2.7 Bonet & Burton 1998, average nodal pressure tetrahedron.** **Verdict: build on (explicit on tets before hexes exist).**
- **Citations.**
  - J. Bonet, A. J. Burton. *A simple average nodal pressure tetrahedral element for incompressible and nearly incompressible dynamic explicit applications.* Commun. Numer. Methods Eng. 14(5):437–449, 1998. DOI [10.1002/(SICI)1099-0887(199805)14:5<437::AID-CNM162>3.0.CO;2-W](https://doi.org/10.1002/(SICI)1099-0887(199805)14:5%3C437::AID-CNM162%3E3.0.CO;2-W).
  - J. Bonet, H. Marriott, O. Hassan. *An averaged nodal deformation gradient linear tetrahedral element for large strain explicit dynamic applications.* CNME 17:551–561, 2001. DOI [10.1002/cnm.429](https://doi.org/10.1002/cnm.429).
- **Method.** Linear tet with the volumetric strain (pressure) averaged over the elements sharing a node, which
  removes the volumetric locking of tet4 in nearly incompressible (plastic, rubber) flow. The 2001 paper averages the
  whole deformation gradient.
- **Solves.** Explicit on the tet meshes M2 delivers before the hex mesher (M3) exists.
- **Limits.** Known pressure-checkerboard/spurious low-energy modes in some settings (nodal-integration family);
  needs a small stabilisation. Plain tet4 must be offered only with a warning for plasticity.
- **Complexity.** Low–moderate (a nodal gather pass per step).

### A3. Shell elements

**A3.1 Belytschko, Lin, Tsay 1984, explicit shells (BLT).** **Verdict: build on (default shell).**
- **Citation.** T. Belytschko, J. I. Lin, C.-S. Tsay. *Explicit algorithms for the nonlinear dynamics of shells.* CMAME 42(2):225–251, 1984. DOI [10.1016/0045-7825(84)90026-4](https://doi.org/10.1016/0045-7825(84)90026-4).
- **Method.** Four-node quadrilateral with one in-plane integration point, Mindlin plate kinematics in a co-rotational
  local frame, several through-thickness points for plasticity, and hourglass control. Cost is a small multiple of a
  solid element's per thickness point. The workhorse of crash simulation.
- **Limits.** Fails the patch test for warped elements (the co-rotational flat-element assumption); hourglass control
  again; poor for strongly warped geometry.
- **Complexity.** Moderate.

**A3.2 Belytschko, Wong, Chiang 1992, advances in one-point quadrature shells.** **Verdict: build on (warping fix for BLT).**
- **Citation.** CMAME 96(1):93–107, 1992. DOI [10.1016/0045-7825(92)90100-X](https://doi.org/10.1016/0045-7825(92)90100-X).
- **Method.** Corrections to BLT for warped geometry (and physical hourglass stabilisation), at little extra cost.

**A3.3 Hughes & Liu 1981, degenerated continuum shells.** **Verdict: reference (second shell, later).**
- **Citation.** T. J. R. Hughes, W. K. Liu. *Nonlinear finite element analysis of shells: Part I. Three-dimensional shells.* CMAME 26(3):331–362, 1981. DOI [10.1016/0045-7825(81)90121-3](https://doi.org/10.1016/0045-7825(81)90121-3).
- **Method.** Shell derived from a 3D solid by constraining fibres (degenerated continuum); handles warping
  naturally and thickness change, about 2–4× the cost of BLT with one in-plane point.
- **Use.** Accuracy option once BLT is verified.

### A4. Large deformation constitutive integration

**A4.1 Simo & Hughes, *Computational Inelasticity*, Springer 1998.** **Verdict: build on.**
- **Links.** DOI [10.1007/b98904](https://doi.org/10.1007/b98904).
- **Method.** Return mapping (closest point projection) for J2 plasticity; radial return for von Mises with
  isotropic/kinematic hardening is closed form for linear hardening and a scalar Newton iteration for nonlinear
  hardening; incrementally objective integration of hypoelastic–plastic models (rotate stress to the configuration at
  the midstep, update, rotate back); consistent tangents (needed by implicit only).
- **Solves.** The stress update at every integration point.
- **Limits.** Radial return in explicit is unconditionally stable but first-order accurate in the strain increment;
  fine since explicit steps are tiny.
- **Complexity.** Low for J2 radial return; moderate for the objective integration.

**A4.2 Krieg & Krieg 1977, accuracy of radial return.** **Verdict: reference.** J. Press. Vessel Technol. 99(4):510–515, 1977. DOI [10.1115/1.3454568](https://doi.org/10.1115/1.3454568). Shows radial return is accurate for the elastic–perfectly plastic model; the standard justification.

**A4.3 Hughes & Winget 1980, incrementally objective integration.** **Verdict: build on.**
- **Citation.** IJNME 15(12):1862–1867, 1980. DOI [10.1002/nme.1620151210](https://doi.org/10.1002/nme.1620151210).
- **Method.** Rotate the stress with the incremental rotation from the midpoint velocity gradient (a Cayley transform
  of the spin increment), so a rigid rotation leaves the stress exact. This is the Jaumann-rate update most explicit
  codes use.
- **Complexity.** Low.

**A4.4 Dienes 1979, rotation and stress rate.** **Verdict: build on (as a verification case).**
- **Citation.** J. K. Dienes. *On the analysis of rotation and stress rate in deforming bodies.* Acta Mech. 32:217–232, 1979. DOI [10.1007/BF01379008](https://doi.org/10.1007/BF01379008).
- **Method.** Simple shear of a hypoelastic material: the Jaumann rate gives an oscillating shear stress (unphysical),
  the Green–Naghdi rate a monotonic one. Both have closed-form solutions.
- **Use.** Analytic verification of the stress-rate implementation (exactly the kind of test the project wants).

**A4.5 Flanagan & Taylor 1987, Green–Naghdi stress integration.** **Verdict: build on (second rate).**
- **Citation.** D. P. Flanagan, L. M. Taylor. *An accurate numerical algorithm for stress integration with finite rotations.* CMAME 62(3):305–320, 1987. DOI [10.1016/0045-7825(87)90065-X](https://doi.org/10.1016/0045-7825(87)90065-X).
- **Method.** Integrates the rotation of the polar decomposition incrementally (no full polar decomposition each step)
  and updates stress in the unrotated frame: Green–Naghdi rate at low cost. Used by the PRONTO family.
- **Limits.** Drift in the rotation over very long runs (re-orthonormalise).

### A5. Contact

**A5.1 Hallquist, Goudreau, Benson 1985, sliding interfaces with contact-impact.** **Verdict: build on.**
- **Citation.** CMAME 51(1–3):107–137, 1985. DOI [10.1016/0045-7825(85)90030-1](https://doi.org/10.1016/0045-7825(85)90030-1).
- **Method.** Node-to-segment penalty contact for explicit: each slave node finds its nearest master segment, a
  penetration produces a penalty force with stiffness scaled from the segment's bulk modulus, area and volume (so it
  does not reduce the stable step), symmetric (two-pass) treatment, Coulomb friction with a trial-force return,
  incremental search tracking the previous master segment. This is the original DYNA3D contact.
- **Solves.** General impact between deformable bodies.
- **Limits.** Penalty contact allows small penetration and can be noisy; node-to-segment fails at segment edges and
  corners (snagging, nodes slipping through when the search loses track); initial penetrations must be removed or
  reported; very soft-against-stiff pairs need soft-constraint penalty (scale on nodal mass/Δt²).
- **Complexity.** Moderate for the force; high for robust search.

**A5.2 Benson & Hallquist 1990, single-surface contact.** **Verdict: build on (self-contact for buckling).**
- **Citation.** CMAME 78(2):141–163, 1990. DOI [10.1016/0045-7825(90)90098-7](https://doi.org/10.1016/0045-7825(90)90098-7).
- **Method.** Self-contact of one surface (folding shells) with bucket-sort global search; no master/slave pairs needed.
- **Solves.** Crush and post-buckling.

**A5.3 Belytschko & Neal 1991, pinball contact.** **Verdict: build on (simplest robust option first).**
- **Citation.** T. Belytschko, M. O. Neal. *Contact-impact by the pinball algorithm with penalty and Lagrangian methods.* IJNME 31(3):547–572, 1991. DOI [10.1002/nme.1620310309](https://doi.org/10.1002/nme.1620310309).
- **Method.** Embed spheres ("pinballs") in surface elements; contact is sphere–sphere overlap, which is trivial and
  vectorises; penalty or Lagrange multiplier enforcement; hierarchical pinballs refine near contact.
- **Solves.** Robust contact detection without the edge/corner cases of node-to-segment; ideal for erosion
  (A7) where the surface keeps changing.
- **Limits.** Contact surface is lumpy; accuracy of the contact pressure is low.
- **Complexity.** Low.

**A5.4 Contact search: Oldenburg & Nilsson 1994 position code; Yang & Laursen 2006 BVH.** **Verdict: build on (BVH).**
- **Citations.**
  - M. Oldenburg, L. Nilsson. *The position code algorithm for contact searching.* IJNME 37(3):359–386, 1994. DOI [10.1002/nme.1620370302](https://doi.org/10.1002/nme.1620370302).
  - B. Yang, T. A. Laursen. *A contact searching algorithm including bounding volume trees applied to finite sliding mortar formulations.* Comput. Mech. 41:189–205, 2008 (online 2006). DOI [10.1007/s00466-006-0116-z](https://doi.org/10.1007/s00466-006-0116-z).
- **Method.** Global search finds candidate pairs: bucket/position-code sorting into a spatial grid (O(n log n), good for
  uniform meshes), or a bounding-volume hierarchy refitted each step and rebuilt occasionally (robust to non-uniform
  meshes and large motion). Local search then computes the closest point on the candidate segment.
- **Choice.** A BVH with AABBs inflated by `v·Δt·k` and refitting, deterministic in build order. The mesher already
  needs spatial indexing; share it if its API fits.
- **Complexity.** Moderate.

**A5.5 Books: Wriggers; Laursen; Hughes et al. 1976.** **Verdict: reference.**
- P. Wriggers. *Computational Contact Mechanics*, 2nd ed., Springer 2006. DOI [10.1007/978-3-540-32609-0](https://doi.org/10.1007/978-3-540-32609-0).
- T. A. Laursen. *Computational Contact and Impact Mechanics*, Springer 2003. DOI [10.1007/978-3-662-04864-1](https://doi.org/10.1007/978-3-662-04864-1).
- T. J. R. Hughes et al. *A finite element method for a class of contact-impact problems.* CMAME 8(3):249–276, 1976. DOI [10.1016/0045-7825(76)90018-9](https://doi.org/10.1016/0045-7825(76)90018-9). (The impact–release conditions; analytic bar-impact test.)
- Use: mortar/segment-to-segment contact for the implicit solver later; theory of penalty regularisation.

### A6. Johnson–Cook plasticity and damage

**A6.1 Johnson & Cook 1983 (strength) and 1985 (fracture).** **Verdict: build on.**
- **Citations.**
  - G. R. Johnson, W. H. Cook. *A constitutive model and data for metals subjected to large strains, high strain rates and high temperatures.* Proc. 7th Int. Symp. on Ballistics, The Hague, 541–547, 1983. **No DOI; (unverified)** — conference proceedings, not online through Crossref.
  - G. R. Johnson, W. H. Cook. *Fracture characteristics of three metals subjected to various strains, strain rates, temperatures and pressures.* Eng. Fract. Mech. 21(1):31–48, 1985. DOI [10.1016/0013-7944(85)90052-9](https://doi.org/10.1016/0013-7944(85)90052-9).
- **Method.** Flow stress `σ_y = (A + B ε_p^n)(1 + C ln ε̇*)(1 − T*^m)`, a multiplicative split of hardening, rate and
  thermal softening (homologous temperature T*; adiabatic heating from plastic work with a Taylor–Quinney factor).
  Fracture strain `ε_f = (D1 + D2 e^{D3 σ*})(1 + D4 ln ε̇*)(1 + D5 T*)` with stress triaxiality σ*; damage
  `D = Σ Δε_p / ε_f` accumulates and the point fails at D = 1.
- **Solves.** Metals under impact, machining, ballistics.
- **Limits.** Rate term `ln ε̇*` is singular at ε̇* → 0 (clamp at 1); T* outside [0, 1] must be clamped; radial return
  with rate dependence needs a scalar Newton iteration (bounded). Fracture is mesh-dependent without regularisation.
  Parameter sets are taken from published papers (the 1983/1985 tables), never from a commercial material library.
- **Complexity.** Low on top of A4.1.

### A7. Element erosion

**A7.1 Johnson & Stryk 1987; Belytschko & Lin 1987.** **Verdict: build on (simple erosion), reference (eroding contact).**
- **Citations.**
  - G. R. Johnson, R. A. Stryk. *Eroding interface and improved tetrahedral element algorithms for high-velocity impact computations in three dimensions.* Int. J. Impact Eng. 5:411–421, 1987. DOI [10.1016/0734-743X(87)90057-1](https://doi.org/10.1016/0734-743X(87)90057-1).
  - T. Belytschko, J. I. Lin. *A three-dimensional impact-penetration algorithm with erosion.* Int. J. Impact Eng. 5:111–127, 1987. DOI [10.1016/0734-743X(87)90033-9](https://doi.org/10.1016/0734-743X(87)90033-9).
- **Method.** Delete an element when a criterion is met (damage D = 1, effective plastic strain, minimum Δt); keep its
  nodal mass and momentum (optionally as free particles); rebuild the exposed contact surface so the projectile keeps
  contacting interior faces.
- **Limits.** Not mass/energy conservative unless the eroded energy is accounted for (report it); results depend on
  the erosion threshold. Contact surface update is the hard part.
- **Complexity.** Low for deletion, moderate for the surface update.

### A8. Artificial bulk viscosity and shock capturing

**A8.1 von Neumann & Richtmyer 1950.** **Verdict: build on.**
- **Citation.** J. von Neumann, R. D. Richtmyer. *A Method for the Numerical Calculation of Hydrodynamic Shocks.* J. Appl. Phys. 21(3):232–237, 1950. DOI [10.1063/1.1699639](https://doi.org/10.1063/1.1699639).
- **Method.** Add a pressure `q` (only in compression) quadratic in the velocity jump across the element, spreading
  shocks over a few elements so the Rankine–Hugoniot jump is reached without oscillation.

**A8.2 Wilkins 1980, artificial viscosity in multidimensions.** **Verdict: build on.**
- **Citation.** M. L. Wilkins. *Use of artificial viscosity in multidimensional fluid dynamic calculations.* J. Comput. Phys. 36(3):281–303, 1980. DOI [10.1016/0021-9991(80)90161-8](https://doi.org/10.1016/0021-9991(80)90161-8).
- **Method.** The combined form used by explicit FE codes: `q = ρ L (C_q L (ε̇_kk)² − C_l c ε̇_kk)` for `ε̇_kk < 0`,
  with a quadratic (C_q ≈ 1.5) and a linear (C_l ≈ 0.06) term; the stable step becomes
  `Δt = L / (Q + √(Q² + c²))`, `Q = C_l c + C_q L |ε̇_kk|`. Default on in every explicit solid element: without it,
  impact runs ring.
- **Limits.** Smears shocks over 3–4 elements; spurious heating in shockless compression (wall heating, see Noh).

**A8.3 Caramana, Shashkov, Whalen 1998.** **Verdict: reference (later, for hydro).**
- **Citation.** JCP 144(1):70–97, 1998. DOI [10.1006/jcph.1998.5989](https://doi.org/10.1006/jcph.1998.5989).
- **Method.** Edge-based tensor viscosity with limiters for multidimensional Lagrangian hydro; avoids spurious
  viscosity in uniform compression and preserves symmetry.
- **Code.** Laghos (BSD-2) has related high-order viscosities: READ OK.

### A9. Hydrocode references

**A9.1 Wilkins, *Computer Simulation of Dynamic Phenomena*, Springer 1999.** **Verdict: build on (hydro side of the explicit solver).**
- **Links.** DOI [10.1007/978-3-662-03885-7](https://doi.org/10.1007/978-3-662-03885-7).
- **Method.** The HEMP lineage: staggered Lagrangian grid, the stress/energy update order, artificial viscosity,
  elastic–plastic flow with a pressure from an EOS and a deviator from a hypoelastic law with radial return, energy
  equation coupling, programmed-burn detonation, sliding interfaces; with many worked impact problems.
- **Complexity.** Moderate, and the order of operations is exactly what the explicit step needs.

**A9.2 Benson 1992, Lagrangian and Eulerian hydrocodes.** **Verdict: build on (survey and design reference).**
- **Citation.** D. J. Benson. *Computational methods in Lagrangian and Eulerian hydrocodes.* CMAME 99(2–3):235–394, 1992. DOI [10.1016/0045-7825(92)90042-I](https://doi.org/10.1016/0045-7825(92)90042-I).
- **Method.** 160-page survey: Lagrangian step, artificial viscosity, EOS coupling and the energy update, Eulerian and
  ALE advection (operator split: Lagrangian step then remap), mixed-material cells, interface reconstruction.
- **Use.** The map for whether and how to add ALE/Eulerian multi-material (blast, fluids with EOS). Read before
  designing anything Eulerian.

**A9.3 Benson 1989, simple ALE.** **Verdict: reference (later).** CMAME 72(3):305–350, 1989. DOI [10.1016/0045-7825(89)90003-0](https://doi.org/10.1016/0045-7825(89)90003-0). An operator-split ALE added to a Lagrangian FE code with a donor-cell/van Leer advection pass.

### A10. Verification problems (analytic or published experimental, no commercial results)

| Case | Reference | Expected value |
|---|---|---|
| Elastic wave in a bar (step load on a free or fixed end) | d'Alembert; A1.1 book | Wave front at `x = c t`, `c = √(E/ρ)` (rod) or `√((λ+2μ)/ρ)` (constrained); stress `σ = ρ c v`; reflection doubling at a fixed end. Measures dispersion/overshoot of the element and q. |
| Elastic bar impact on a rigid wall | Hughes et al. 1976 | Contact duration `2L/c`, release velocity equal and opposite (energy conserved). |
| Patch test, rigid rotation | A2.1, A4.3 | Constant strain reproduced; zero stress after a 90° rigid rotation. |
| Simple shear, hypoelastic | Dienes 1979 | Closed-form shear stress for Jaumann (oscillating) and Green–Naghdi rates. |
| Uniaxial J2 / Johnson–Cook at constant rate | A4.1, A6.1 | Closed-form flow curve from the model's own equations. |
| Cantilever free vibration | Euler–Bernoulli / M6 modal solver | Period from the analytic first frequency; cross-check with the implicit modal solver. |
| Energy balance | A1.1 | `KE + IE + hourglass + contact + eroded − W_ext` constant within 1%. |
| Taylor bar impact | G. I. Taylor 1948, Proc. R. Soc. A 194:289–299, DOI [10.1098/rspa.1948.0081](https://doi.org/10.1098/rspa.1948.0081); M. L. Wilkins, M. W. Guinan 1973, J. Appl. Phys. 44(3):1200–1206, DOI [10.1063/1.1662328](https://doi.org/10.1063/1.1662328) | Taylor's 1D estimate is analytic (final length ratio vs ρv²/σ_y). Wilkins–Guinan experimental final lengths are published experiments: allowed by the "no commercial data" rule, but **owner decision** whether experiment-derived expected values may enter the reference suite. |
| Sod, Noh, Sedov (hydro) | Sod 1978 JCP 27(1):1–31, DOI [10.1016/0021-9991(78)90023-2](https://doi.org/10.1016/0021-9991(78)90023-2); Noh 1987 JCP 72(1):78–120, DOI [10.1016/0021-9991(87)90074-X](https://doi.org/10.1016/0021-9991(87)90074-X) | Exact Riemann solution; Noh's exact post-shock density 64 (spherical, γ=5/3); Sedov self-similar solution. Tests q and EOS coupling. |

---

## Part B: Equations of state

All EOS models below give pressure `p(ρ, e)` (e: specific internal energy, or energy per reference volume E);
the deviatoric stress comes from the strength model (A4/A6). Every EOS in the solver needs: finite-input checks,
the sound speed `c² = (∂p/∂ρ)_s = (∂p/∂ρ)_e + (p/ρ²)(∂p/∂e)_ρ` (for Δt and q), a tension cut-off, and a guard on
compression singularities.

**B1. Mie–Grüneisen with a linear Us–up Hugoniot.** **Verdict: build on.**
- **Citations.**
  - R. G. McQueen, S. P. Marsh, J. W. Taylor, J. N. Fritz, W. J. Carter. *The equation of state of solids from shock wave studies.* In *High-Velocity Impact Phenomena* (R. Kinslow, ed.), Academic Press 1970, 293–417. DOI [10.1016/B978-0-12-408950-1.50012-4](https://doi.org/10.1016/B978-0-12-408950-1.50012-4).
  - M. A. Meyers. *Dynamic Behavior of Materials.* Wiley 1994. DOI [10.1002/9780470172278](https://doi.org/10.1002/9780470172278). (Textbook derivation of Rankine–Hugoniot and Mie–Grüneisen.)
- **Method.** The Hugoniot `U_s = c_0 + s u_p` and Rankine–Hugoniot give, with `η = 1 − ρ_0/ρ`,
  `p_H = ρ_0 c_0² η / (1 − s η)²` and `e_H = p_H η / (2ρ_0)`; off the Hugoniot,
  `p = p_H + Γ ρ (e − e_H)` with `Γ ρ = Γ_0 ρ_0` (or `Γ = Γ_0 + a η`). Quadratic/cubic Hugoniot fits extend it.
- **Solves.** Metals and water/polymers under shock (impact, penetration, underwater explosion water).
- **Limits.** Singular at `η = 1/s` (guard: error or cap, never NaN); tension branch needs a separate form (linear
  `p = ρ_0 c_0² η`) and a spall cut-off; invalid near phase transitions and in expansion to vapour.
- **Complexity.** Low.

**B2. JWL for detonation products.** **Verdict: build on.**
- **Citations.**
  - E. L. Lee, H. C. Hornig, J. W. Kury. *Adiabatic Expansion of High Explosive Detonation Products.* LLNL report UCRL-50422, 1968. DOI [10.2172/4783904](https://doi.org/10.2172/4783904) (OSTI).
  - B. M. Dobratz. *LLNL explosives handbook: properties of chemical explosives and explosive simulants.* UCRL-52997, 1981. DOI [10.2172/6530310](https://doi.org/10.2172/6530310) (published JWL parameter sets).
- **Method.** `p = A(1 − ω/(R_1 V)) e^{−R_1 V} + B(1 − ω/(R_2 V)) e^{−R_2 V} + ω E / V`, V = relative volume, E =
  energy per initial volume; fitted to cylinder-test expansion. Used with programmed burn: each element's lighting
  time = distance from the detonation point / detonation speed D, a burn fraction ramps the pressure in
  (Wilkins A9.1 describes this).
- **Solves.** High-explosive products; with B1 water, underwater explosion; with an Eulerian/ALE or SPH carrier, blast.
- **Limits.** Programmed burn ignores shock initiation and shadowing (needs a ray/eikonal lighting time around
  obstacles); in a Lagrangian mesh the products' expansion tangles the mesh quickly, so practical blast needs
  ALE/Eulerian (A9.2/A9.3) or SPH (E2).
- **Complexity.** Low for the EOS; programmed burn moderate.

**B3. Ideal gas / gamma law, polynomial, stiffened gas.** **Verdict: build on.**
- Ideal gas `p = (γ − 1) ρ e`. Linear polynomial `p = C_0 + C_1 μ + C_2 μ² + C_3 μ³ + (C_4 + C_5 μ + C_6 μ²) E`,
  `μ = ρ/ρ_0 − 1` (a general form that contains ideal gas and linear elastic as special cases; a hydrocode
  convention, described in A9.2 — write our own description). Stiffened gas `p = (γ − 1) ρ e − γ p_∞`:
  O. Le Métayer, J. Massoni, R. Saurel, Int. J. Therm. Sci. 43(3):265–276, 2004, DOI [10.1016/j.ijthermalsci.2003.09.002](https://doi.org/10.1016/j.ijthermalsci.2003.09.002) (French title in Crossref).
- **Limits.** Polynomial may give negative sound speed in tension for arbitrary coefficients: validate c² > 0 over the
  admissible range at model check time.
- **Complexity.** Trivial.

**B4. Tabulated EOS.** **Verdict: build on (simple bilinear table), reference (SESAME-class).**
- **Method.** `p(ρ, e)` (or `p(ρ, T)`, `e(ρ, T)`) on a rectilinear grid with bilinear or monotone bicubic
  interpolation; for a 1D compaction table `p(ε_v)` with unloading slopes (porous materials, soils).
- **Limits.** Thermodynamic consistency (positive c², convexity) is not guaranteed by the user's table: check at
  model-check time and report the bad cells. The LANL SESAME tables carry their own distribution terms: do not
  bundle any.
- **Reference.** R. Menikoff, B. J. Plohr. *The Riemann problem for fluid flow of real materials.* Rev. Mod. Phys. 61(1):75–130, 1989. DOI [10.1103/RevModPhys.61.75](https://doi.org/10.1103/RevModPhys.61.75) — convexity (fundamental derivative) conditions an EOS must meet for well-posed shocks.

**B5. Tillotson.** **Verdict: reference (later; for hypervelocity impact).**
- **Citations.**
  - J. H. Tillotson. *Metallic equations of state for hypervelocity impact.* General Atomic report GA-3216, 1962. **(unverified: no DOI or OSTI record found by API search)**.
  - A. L. Brundage. *Implementation of Tillotson Equation of State for Hypervelocity Impact of Metals, Geologic Materials, and Liquids.* Procedia Eng. 58:461–470, 2013. DOI [10.1016/j.proeng.2013.05.053](https://doi.org/10.1016/j.proeng.2013.05.053).
- **Method.** Separate forms for compressed, cold-expanded and hot-expanded (vapour) states with an interpolation
  region between incipient and complete vaporisation energies.
- **Limits.** The region switch is discontinuous in published forms; Brundage documents the fix. Needed only for
  hypervelocity (> ~5 km/s) work.
- **Complexity.** Low–moderate.

**B6. Energy equation coupling and pressure–energy iteration.** **Verdict: build on (A9.1, A9.2).**
- **Method.** Per element, per step: `e^{n+1} = e^n − ((p^{n+1} + p^n)/2 + q) ΔV/m + (deviatoric work)/m`. Since
  `p^{n+1}` depends on `e^{n+1}`, the update is implicit. For EOS linear in energy (`p = f(ρ) + g(ρ) e`: ideal gas,
  polynomial, Mie–Grüneisen, JWL) it is solved in closed form; for tabulated or Tillotson, a bounded scalar
  fixed-point or Newton iteration (≤ 10 iterations, report non-convergence).
- **Limits.** Energy conservation of the whole scheme depends on using the same ΔV for work and volume change, and on
  adding q-work to e (shock heating).

**B7. Shock capturing.** Lagrangian explicit captures shocks with artificial viscosity (A8). Godunov/Riemann-solver
shock capturing belongs to Eulerian finite volume (E1); not needed for the Lagrangian solver.

---

## Part C: Lattice Boltzmann method

### C1. Foundations

**C1.1 Chen & Doolen 1998.** **Verdict: reference (orientation).**
- **Citation.** S. Chen, G. D. Doolen. *Lattice Boltzmann method for fluid flows.* Annu. Rev. Fluid Mech. 30:329–364, 1998. DOI [10.1146/annurev.fluid.30.1.329](https://doi.org/10.1146/annurev.fluid.30.1.329).
- **Method.** Review: from lattice-gas automata to BGK LBM, Chapman–Enskog recovery of Navier–Stokes, boundary
  conditions, multiphase and thermal extensions.

**C1.2 Krüger, Kusumaatmaja, Kuzmin, Shardt, Silva, Viggen. *The Lattice Boltzmann Method: Principles and Practice*, Springer 2017.** **Verdict: build on (the spine of the LBM solver).**
- **Links.** DOI [10.1007/978-3-319-44649-3](https://doi.org/10.1007/978-3-319-44649-3).
- **Content used.** Velocity sets (D2Q9, D3Q15/19/27) and equilibria; Chapman–Enskog; unit conversion and the choice
  of lattice parameters; stream-collide algorithm; BGK, TRT and MRT; bounce-back, interpolated bounce-back, Zou–He
  and other boundary conditions; forces (Guo); stability and accuracy analysis; sound waves. The book has
  companion example code on its site: **licence not checked — do not open until checked**.
- **Method in one line.** `f_i(x + c_i Δt, t + Δt) = f_i(x, t) + Ω_i(f)`; density and momentum are moments of f;
  viscosity `ν = c_s² (τ − ½) Δt` with `c_s² = 1/3` in lattice units.
- **Limits** (general, LBM): weakly compressible — the error is O(Ma²), so lattice velocity `u_lat ≲ 0.1`
  (Ma ≲ 0.17); `τ → ½` (high Re on a coarse grid) destabilises BGK; uniform Cartesian grid, so curved walls are
  staircased unless interpolated BCs are used; acoustic time step `Δt ∝ Δx` (acoustic scaling) or `∝ Δx²`
  (diffusive scaling) makes incompressible low-speed flows costly.
- **Complexity.** BGK D3Q19 with bounce-back on a uniform grid: ≈1k lines. Everything production-grade (refinement,
  turbulence, BCs, GPU) is the real work.

**C1.3 Bhatnagar, Gross, Krook 1954; Qian, d'Humières, Lallemand 1992.** **Verdict: build on (BGK as the first collision).**
- P. L. Bhatnagar, E. P. Gross, M. Krook. Phys. Rev. 94(3):511–525, 1954. DOI [10.1103/PhysRev.94.511](https://doi.org/10.1103/PhysRev.94.511).
- Y. H. Qian, D. d'Humières, P. Lallemand. *Lattice BGK Models for Navier-Stokes Equation.* EPL 17(6):479, 1992. DOI [10.1209/0295-5075/17/6/001](https://doi.org/10.1209/0295-5075/17/6/001).
- **Method.** Single relaxation toward a second-order Hermite equilibrium; the DnQm lattices.
- **Limits.** Viscosity-dependent wall location with bounce-back; unstable at low viscosity.

**C1.4 He & Luo 1997, incompressible LBM.** **Verdict: reference.** J. Stat. Phys. 88:927–944, 1997. DOI [10.1023/B:JOSS.0000015179.12689.e4](https://doi.org/10.1023/B:JOSS.0000015179.12689.e4). An equilibrium that removes O(Ma²) compressibility error in steady flows; useful for steady validation.

**C1.5 Sterling & Chen 1996; Lallemand & Luo 2000 — stability.** **Verdict: reference.**
- J. D. Sterling, S. Chen. *Stability Analysis of Lattice Boltzmann Methods.* JCP 123(1):196–206, 1996. DOI [10.1006/jcph.1996.0016](https://doi.org/10.1006/jcph.1996.0016).
- P. Lallemand, L.-S. Luo. *Theory of the lattice Boltzmann method: Dispersion, dissipation, isotropy, Galilean invariance, and stability.* PRE 61:6546–6562, 2000. DOI [10.1103/PhysRevE.61.6546](https://doi.org/10.1103/PhysRevE.61.6546).
- **Use.** Linear stability (von Neumann) analysis of BGK/MRT; the source for the stable (u, τ) region and the case for MRT.

**C1.6 Guo, Zheng, Shi 2002, forcing.** **Verdict: build on.** PRE 65:046308, 2002. DOI [10.1103/PhysRevE.65.046308](https://doi.org/10.1103/PhysRevE.65.046308). The second-order body-force scheme every IB-LBM coupling (D2) needs.

### C2. Collision operators

**C2.1 d'Humières et al. 2002, MRT in 3D.** **Verdict: build on (second collision).**
- **Citation.** D. d'Humières, I. Ginzburg, M. Krafczyk, P. Lallemand, L.-S. Luo. *Multiple-relaxation-time lattice Boltzmann models in three dimensions.* Phil. Trans. R. Soc. A 360(1792):437–451, 2002. DOI [10.1098/rsta.2001.0955](https://doi.org/10.1098/rsta.2001.0955).
- **Method.** Collide in moment space with a separate relaxation rate per moment group (D3Q15/D3Q19 matrices given);
  free rates damp non-hydrodynamic modes, raising stability and fixing the wall location.
- **Limits.** Still loses stability at very high Re; parameters need care.
- **Complexity.** Low–moderate (a constant 19×19 transform).

**C2.2 Ginzburg, Verhaeghe, d'Humières 2008, TRT.** **Verdict: build on (default over BGK).**
- **Citation.** Commun. Comput. Phys. 3(2):427–478, 2008. DOI [10.4208/cicp.2008.v3.p427](https://doi.org/10.4208/cicp.2008.v3.p427).
- **Method.** Two relaxation rates (symmetric/antisymmetric parts); the "magic parameter" Λ = 3/16 puts the
  bounce-back wall exactly halfway and makes the solution viscosity-independent in steady flow. BGK cost, much of
  MRT's benefit.
- **Complexity.** Low.

**C2.3 Geier et al. 2015, cumulant LBM.** **Verdict: build on (high-Re collision, D3Q27).**
- **Citations.**
  - M. Geier, M. Schönherr, A. Pasquali, M. Krafczyk. *The cumulant lattice Boltzmann equation in three dimensions: Theory and validation.* Comput. Math. Appl. 70(4):507–547, 2015. DOI [10.1016/j.camwa.2015.05.001](https://doi.org/10.1016/j.camwa.2015.05.001).
  - M. Geier, A. Pasquali, M. Schönherr. *Parametrization of the cumulant lattice Boltzmann method for fourth order accurate diffusion part I.* JCP 348:862–888, 2017. DOI [10.1016/j.jcp.2017.05.040](https://doi.org/10.1016/j.jcp.2017.05.040).
- **Method.** Relax statistically independent cumulants (computed via a central-moment chain transform on D3Q27);
  Galilean invariant and very stable; the 2017 parametrisation chooses the free rates so the scheme is fourth-order
  accurate in diffusion, and acts as an implicit LES at high Re without an explicit sub-grid model.
- **Limits.** D3Q27 (27 populations: ≈40% more memory and bandwidth than D3Q19); the transform is long (hundreds of
  flops) but bandwidth, not flops, limits LBM; a limiter is needed for under-resolved flows.
- **Complexity.** Moderate–high (careful transcription of the transforms; test each against the moment definitions).
- **Code.** VirtualFluids and lbmpy implement it: both GPL-3.0, FORBIDDEN. Implement from the paper.

**C2.4 Regularised and recursive-regularised collision.** **Verdict: built (recursive regularised on D3Q19, `collision: "regularized"`, #557; `docs/lbm-collision.md`).**
- **Citations.**
  - J. Latt, B. Chopard. *Lattice Boltzmann method with regularized pre-collision distribution functions.* Math. Comput. Simul. 72(2–6):165–168, 2006. DOI [10.1016/j.matcom.2006.05.017](https://doi.org/10.1016/j.matcom.2006.05.017).
  - O. Malaspinas. *Increasing stability and accuracy of the lattice Boltzmann scheme: recursivity and regularization.* arXiv:1505.06900, 2015. https://arxiv.org/abs/1505.06900 (200).
    (A preprint, not on Crossref: the implementation cites Coreixas et al. 2017 for the recursive form.)
  - C. Coreixas, G. Wissocq, G. Puigt, J.-F. Boussuge, P. Sagaut. *Recursive regularization step for high-order lattice Boltzmann methods.* PRE 96:033306, 2017. DOI [10.1103/PhysRevE.96.033306](https://doi.org/10.1103/PhysRevE.96.033306).
  - C. Coreixas, B. Chopard, J. Latt. *Comprehensive comparison of collision models in the lattice Boltzmann framework: Theory and investigations.* PRE 100:033305, 2019. DOI [10.1103/PhysRevE.100.033305](https://doi.org/10.1103/PhysRevE.100.033305).
- **Method.** Before collision, project the non-equilibrium part onto its second-order Hermite component (regularised),
  or reconstruct higher-order Hermite coefficients recursively from it (recursive regularised): filters the ghost modes
  that destabilise BGK. Coreixas 2019 compares BGK, MRT (raw/central/Hermite moments), cumulant, regularised and
  recursive-regularised on the same linear stability and benchmark cases — **the paper to pick the operator with**.
- **Limits.** Regularised adds numerical dissipation; RR is near-cumulant in stability on D3Q19/27.
- **Code.** STLBM (MIT) has BGK/TRT/regularised kernels: READ OK.

**C2.5 Velocity sets.** D3Q19 for laminar and moderate Re (cheapest isotropic-enough set); D3Q27 for cumulant and
wall-bounded turbulence. Kang & Hassan 2013 show D3Q19 artefacts in wall-bounded turbulence: S. K. Kang, Y. A. Hassan,
JCP 232(1):100–117, 2013, DOI [10.1016/j.jcp.2012.07.023](https://doi.org/10.1016/j.jcp.2012.07.023). **Verdict: reference.**

### C3. Boundary conditions

**C3.1 Bounce-back; Ladd 1994 moving walls.** **Verdict: build on.**
- A. J. C. Ladd. *Numerical simulations of particulate suspensions via a discretized Boltzmann equation. Part 1.* JFM 271:285–309, 1994. DOI [10.1017/S0022112094001771](https://doi.org/10.1017/S0022112094001771).
- **Method.** Halfway bounce-back reflects populations crossing a wall link: second-order for a wall halfway between
  nodes, staircase otherwise. Ladd adds the moving-wall momentum term (`−2 w_i ρ c_i·u_w / c_s²`) and the
  momentum-exchange force on the body, which is how LBM computes drag/lift and FSI loads with bounce-back.
- **Complexity.** Low.

**C3.2 Bouzidi, Firdaouss, Lallemand 2001, interpolated bounce-back.** **Verdict: build on.**
- **Citation.** Phys. Fluids 13(11):3452–3459, 2001. DOI [10.1063/1.1399290](https://doi.org/10.1063/1.1399290).
- **Method.** For each wall link, the fraction q ∈ (0, 1] of the link inside the fluid; linear or quadratic
  interpolation of the reflected population, second-order for curved walls.
- **Needs from meshing.** q comes from ray–triangle intersection of the lattice link with the face-tagged boundary
  surface — robust predicates (fieldcraft-mesh) give the inside/outside decision; q itself is a constructed value.
- **Limits.** Not mass conservative; small q (< ~0.01) can destabilise: treat as q = 0 or use the Mei variant.
- **Also.** R. Mei, L.-S. Luo, W. Shyy, JCP 155(2):307–330, 1999, DOI [10.1006/jcph.1999.6334](https://doi.org/10.1006/jcph.1999.6334) — the alternative curved-boundary scheme. **Reference.**

**C3.3 Zou & He 1997, pressure and velocity BCs.** **Verdict: build on (inlet/outlet).**
- **Citation.** Phys. Fluids 9(6):1591–1598, 1997. DOI [10.1063/1.869307](https://doi.org/10.1063/1.869307).
- **Method.** Unknown incoming populations from the prescribed velocity (or density) plus bounce-back of the
  non-equilibrium part normal to the boundary.
- **Limits.** Defined for flat, axis-aligned boundaries; corners need special treatment; unstable at high Re outlets —
  use regularised or convective/sponge outlets there.

**C3.4 Orlanski 1976, convective outflow condition.** **Verdict: build on (convective outlets, #599).**
- **Citation.** I. Orlanski, *A simple boundary condition for unbounded hyperbolic flows*, J. Comput. Phys.
  21(3):251–269, 1976 (DOI not checked).
- **Method.** At an open boundary every quantity obeys `∂φ/∂t + U ∂φ/∂n = 0`, discretised upwind and implicitly
  in time, `φ_N^{t+1} = (φ_N^t + λ φ_{N−1}^{t+1}) / (1 + λ)`, `λ = U Δt/Δx`; Orlanski estimates `U` from the
  interior, others prescribe it. For lattice Boltzmann it is applied to the unknown populations: Q. Lou, Z. Guo,
  B. Shi, PRE 87:063301, 2013; Z. Yang, Comput. Math. Appl. 65(2):160–171, 2013 (DOIs not checked).
- **Used as.** `Outflow::Convective` on a pressure face: the unknown populations convected with the cell's own
  outward speed, the cell then held at the face's density (`w_q δρ` on every population, the momentum unchanged),
  Zou–He where the flow comes back in (`fieldcraft_lbm::outlet`).
- **Limits.** Holding each cell at the face's pressure still reflects part of a strong vortex's pressure field;
  a characteristic (non-reflecting) outlet would pass it.

### C4. Grid refinement

**C4.1 Filippova & Hänel 1998; Dupuis & Chopard 2003.** **Verdict: build on (Dupuis–Chopard cell-vertex refinement).**
- O. Filippova, D. Hänel. *Grid Refinement for Lattice-BGK Models.* JCP 147(1):219–228, 1998. DOI [10.1006/jcph.1998.6089](https://doi.org/10.1006/jcph.1998.6089).
- A. Dupuis, B. Chopard. *Theory and applications of an alternative lattice Boltzmann grid refinement algorithm.* PRE 67:066707, 2003. DOI [10.1103/PhysRevE.67.066707](https://doi.org/10.1103/PhysRevE.67.066707).
- **Method.** Nested grids at ratio 2 with time-step ratio 2 (acoustic scaling); τ rescaled per level to keep ν;
  populations at the interface are rescaled non-equilibrium parts with spatial and temporal interpolation from the
  coarse level.
- **Limits.** Interface artefacts (spurious acoustic reflections); mass conservation not exact in the vertex scheme.

**C4.2 Rohde et al. 2006, mass-conservative cell-centred refinement.** **Verdict: build on (with the octree).**
- **Citation.** M. Rohde, D. Kandhai, J. J. Derksen, H. E. A. van den Akker. *A generic, mass conservative local grid refinement technique for lattice-Boltzmann schemes.* IJNMF 51(4):439–468, 2006. DOI [10.1002/fld.1140](https://doi.org/10.1002/fld.1140).
- **Method.** Volumetric (cell-centred) formulation: a coarse cell is split into 8 fine cells, populations are
  distributed homogeneously (explode) and gathered (coalesce) across levels; exactly mass conservative and naturally
  fits an octree with 2:1 balance.
- **Fit.** The M3 octree (2:1-balanced, adapted to the surface) is the same data structure. Reuse it.

**C4.3 Lagrava et al. 2012.** **Verdict: reference.** JCP 231(14):4808–4822, 2012. DOI [10.1016/j.jcp.2012.03.015](https://doi.org/10.1016/j.jcp.2012.03.015). Multi-domain refinement with filtering of the coarse-to-fine transfer (the Palabos approach; the paper is fine to read, the code is AGPL).

### C5. Turbulence

**C5.1 Smagorinsky LES in LBM: Smagorinsky 1963; Hou et al. 1995; Krafczyk, Tölke, Luo 2003.** **Verdict: build on.**
- J. Smagorinsky. Mon. Weather Rev. 91(3):99–164, 1963. DOI [10.1175/1520-0493(1963)091<0099:GCEWTP>2.3.CO;2](https://doi.org/10.1175/1520-0493(1963)091%3C0099:GCEWTP%3E2.3.CO;2).
- S. Hou, J. Sterling, S. Chen, G. D. Doolen. *A lattice Boltzmann subgrid model for high Reynolds number flows.* Fields Inst. Commun. 6:151–166, 1996 (Crossref year 1995). DOI [10.1090/fic/006/12](https://doi.org/10.1090/fic/006/12).
- M. Krafczyk, J. Tölke, L.-S. Luo. *Large-eddy simulations with a multiple-relaxation-time LBE model.* Int. J. Mod. Phys. B 17:33–39, 2003. DOI [10.1142/S0217979203017059](https://doi.org/10.1142/S0217979203017059).
- **Method.** The strain rate is local in LBM (from the non-equilibrium second moment), so `ν_t = (C_s Δ)² |S|` is
  computed without finite differences and enters as a local τ. Cheap and local.
- **Limits.** Smagorinsky is over-dissipative near walls (needs damping or a wall model) and in transition.
- **As implemented** (`crates/lbm/src/turbulence.rs`, #663/#664): Hou et al.'s closed-form local τ from the
  non-equilibrium flux `P = Π^neq + ½(uF + Fu)` of Guo's forcing, in both the regularised and the two-relaxation-time
  collision (τ⁻ follows τ⁺, keeping Λ); with heat, the turbulent diffusivity `α_t = ν_t / Pr_t` adds to each cell's
  τ_T. No wall damping yet. The GPU backend (`crates/lbm-gpu/src/lbm.wgsl`, #665) runs the same closed form in f32
  behind an override constant, so with the model off its kernels and results are bit-for-bit the earlier ones; the
  cell's ν_t reaches the scalar lattice through the moments buffer's spare component (no new buffer). Measured on an
  RTX 4070 Laptop GPU, 128³: within ±2.3 % of the model-off throughput with or without heat (the step is
  memory-bound).

**C5.2 Malaspinas & Sagaut 2012, 2014 — consistent SGS and wall models.** **Verdict: reference (wall model later).**
- *Consistent subgrid scale modelling for lattice Boltzmann methods.* JFM 700:514–542, 2012. DOI [10.1017/jfm.2012.155](https://doi.org/10.1017/jfm.2012.155).
- *Wall model for large-eddy simulation based on the lattice Boltzmann method.* JCP 275:25–40, 2014. DOI [10.1016/j.jcp.2014.06.020](https://doi.org/10.1016/j.jcp.2014.06.020).
- **Method.** Derives the SGS closure consistently from the filtered LB equation; the wall model imposes a log-law
  (Musker profile) wall shear stress through a regularised BC at the first off-wall node.
- **Use.** Wall-modelled LES is what makes industrial-Re external flow tractable; it comes after laminar validation.

### C6. Memory layouts and GPU

**C6.1 Bailey et al. 2009, AA pattern.** **Verdict: build on.**
- **Citation.** P. Bailey, J. Myre, S. D. C. Walsh, D. J. Lilja, M. O. Saar. *Accelerating Lattice Boltzmann Fluid Flow Simulations Using Graphics Processors.* ICPP 2009, 550–557. DOI [10.1109/ICPP.2009.38](https://doi.org/10.1109/ICPP.2009.38).
- **Method.** One population array updated in place by alternating two kernels (even step: read/write own cell
  reversed; odd step: read/write neighbours), halving memory against the two-array (A–B) scheme with no extra traffic.
- **Limits.** Even/odd steps have different index maps; boundaries must be written for both.

**C6.2 Geier & Schönherr 2017, esoteric twist.** **Verdict: reference (alternative to AA).**
- **Citation.** *Esoteric Twist: An Efficient in-Place Streaming Algorithmus for the Lattice Boltzmann Method on Massively Parallel Hardware.* Computation 5(2):19, 2017. DOI [10.3390/computation5020019](https://doi.org/10.3390/computation5020019).
- **Method.** In-place streaming with a single kernel by storing populations at "twisted" locations (each cell owns
  the populations pointing in negative directions of itself and its neighbours), with pointer swaps instead of
  alternating kernels; works on sparse/block grids.
- **Also.** Lehmann 2022, *Esoteric Pull and Esoteric Push*, Computation 10(6):92, DOI [10.3390/computation10060092](https://doi.org/10.3390/computation10060092): simpler variants. The paper is open access and readable; its author's code (FluidX3D) is non-commercial: FORBIDDEN.
- Wittmann et al. 2013, *Comparison of different propagation steps for lattice Boltzmann methods*, Comput. Math. Appl. 65(6):924–935, DOI [10.1016/j.camwa.2012.05.002](https://doi.org/10.1016/j.camwa.2012.05.002): measured comparison of A–B, AA, esoteric twist and others on CPUs. **Reference.**

**C6.3 Tölke & Krafczyk 2008; Latt et al. 2021.** **Verdict: reference.**
- J. Tölke, M. Krafczyk. *TeraFLOP computing on a desktop PC with GPUs for 3D CFD.* Int. J. Comput. Fluid Dyn. 22(7):443–456, 2008. DOI [10.1080/10618560802238275](https://doi.org/10.1080/10618560802238275). (DOI resolves; publisher page 403 to a script.)
- J. Latt, C. Coreixas, J. Beny. *Cross-platform programming model for many-core lattice Boltzmann simulations.* PLoS ONE 16(4):e0250306, 2021. DOI [10.1371/journal.pone.0250306](https://doi.org/10.1371/journal.pone.0250306). Code STLBM: MIT, READ OK — the best readable reference for data layouts (AoS vs SoA, AA) with measured performance.
- J. Latt et al. *Palabos: Parallel Lattice Boltzmann Solver.* Comput. Math. Appl. 81:334–350, 2021. DOI [10.1016/j.camwa.2020.03.022](https://doi.org/10.1016/j.camwa.2020.03.022). Paper: read for architecture (block-structured multi-block, data processors). Code: AGPL, FORBIDDEN.

**Cost model.** LBM is memory-bound: a D3Q19 update moves 2 × 19 populations. In f32 that is 152 B/cell/step,
so a 500 GB/s GPU tops out near 3 G cell-updates/s and a 50 GB/s CPU near 0.3 G (half in f64). Memory: AA-pattern
D3Q19 f32 ≈ 80 B/cell (with flags): 10⁸ cells ≈ 8 GB. These numbers decide the GPU question (E3).

### C7. Free-surface LBM

**C7.1 Körner et al. 2005; Thürey 2007; Thürey & Rüde 2009.** **Verdict: reference (a later milestone).**
- C. Körner, M. Thies, T. Hofmann, N. Thürey, U. Rüde. *Lattice Boltzmann Model for Free Surface Flow for Modeling Foaming.* J. Stat. Phys. 121:179–196, 2005. DOI [10.1007/s10955-005-8879-8](https://doi.org/10.1007/s10955-005-8879-8).
- N. Thürey. *Physically based animation of free surface flows with the lattice Boltzmann method.* PhD thesis, FAU Erlangen-Nürnberg, 2007. PDF https://www10.cs.fau.de/publications/dissertations/Diss_2007-Thuerey.pdf (200, application/pdf; title page not machine-checked).
- N. Thürey, U. Rüde. *Stable free surface flows with the lattice Boltzmann method on adaptively coarsened grids.* Comput. Vis. Sci. (Crossref issue year 2008; volume/pages not checked). DOI [10.1007/s00791-008-0090-4](https://doi.org/10.1007/s00791-008-0090-4).
- **Method.** Volume-of-fluid-like cell types (fluid, interface, gas) with a fill level; missing populations at the
  interface reconstructed from the gas pressure; mass exchange through links; cell conversion with flag
  propagation.
- **Solves.** Sloshing, slamming, filling — the FSI cases engineers ask for with water.
- **Limits.** Single-phase (gas pressure only); mass loss without care; many special cases in cell conversion.
- **Complexity.** High.

### C8. Validation cases

| Case | Reference | Expected |
|---|---|---|
| Poiseuille channel/pipe | analytic | parabolic profile; second-order convergence with halfway BB (TRT magic Λ gives exact steady result) |
| Couette, Womersley | analytic | linear profile; Womersley oscillatory solution |
| Taylor–Green vortex (decay) | analytic | `exp(−2νk²t)` decay; measures viscosity and order |
| Lid-driven cavity Re 100–1000 | U. Ghia, K. N. Ghia, C. T. Shin, JCP 48(3):387–411, 1982, DOI [10.1016/0021-9991(82)90058-4](https://doi.org/10.1016/0021-9991(82)90058-4) | published numerical reference values (not analytic: owner decision, same as Taylor bar) |
| Flow past a cylinder 2D-1/2D-2/3D | M. Schäfer, S. Turek, *Benchmark computations of laminar flow around a cylinder*, NNFM 52:547–566, 1996, DOI [10.1007/978-3-322-89849-4_39](https://doi.org/10.1007/978-3-322-89849-4_39); LBM results: J. Geller et al., Comput. Fluids 35:888–897, 2006, DOI [10.1016/j.compfluid.2005.08.009](https://doi.org/10.1016/j.compfluid.2005.08.009) | drag/lift coefficient intervals and Strouhal number given by the benchmark (published numerical reference) |
| Sphere Stokes drag | analytic | `F = 6πμRU` (with confinement correction) |
| Differentially heated cavity Ra 10³–10⁶ | G. de Vahl Davis, *Natural convection of air in a square cavity: a bench mark numerical solution*, Int. J. Numer. Meth. Fluids 3:249–264, 1983 | hot-wall Nusselt numbers (published numerical reference) |
| Cavity with a centred conducting body (conjugate) | J. M. House, C. Beckermann, T. F. Smith, *Effect of a centered conducting body on natural convection heat transfer in an enclosure*, Numer. Heat Transfer A 18(2):213–225, 1990, DOI [10.1080/10407789008944791](https://doi.org/10.1080/10407789008944791) | hot-wall Nusselt numbers at Ra 10⁵ for conductivity ratios 0.2, 1, 5, and the conduction limit (Table 3; published numerical reference) |
| Square array of cylinders, conduction (conjugate, curved interfaces) | W. T. Perrins, D. R. McKenzie, R. C. McPhedran, *Transport properties of regular arrays of cylinders*, Proc. R. Soc. Lond. A 369:207–225, 1979; the coefficients from the lattice sums of Y. A. Godin, J. Math. Phys. 53:063703, 2012 | effective conductivity series at a given area fraction (analytic, truncated series); Keller's reciprocity k(κ) k(1/κ) = 1 |

### C9. Heat: the scalar lattice and conjugate heat transfer, as implemented

- **Temperature lattice.** D3Q7 passive scalar on the flow's links, equilibrium linear in u
  (T. Krüger et al., *The Lattice Boltzmann Method: Principles and Practice*, Springer 2017,
  ch. 8), two relaxation times with Λ_T fixed (I. Ginzburg, Adv. Water Resour. 28:1171, 2005:
  the antisymmetric rate sets the diffusivity, and with Λ_T fixed the steady solution and its
  anti-bounce-back walls do not depend on it; #629); its velocity-dependent numerical diffusion −(τ_T − ½) u² (B. Chopard, J.-L. Falcone,
  J. Latt, Eur. Phys. J. Special Topics 171:245–249, 2009) bounds the lattice speed of a heat
  problem. Held temperatures by anti-bounce-back, fluxes added to the bounced population.
- **Conduction in solids.** Cell-centred finite volumes on the lattice's cells with the
  harmonic mean of the two cells' conductivities on each face (S. V. Patankar, *Numerical Heat
  Transfer and Fluid Flow*, Hemisphere 1980, ch. 4), explicit and sub-cycled under its stability
  limit (`fieldcraft_lbm::conduction`).
- **Fluid–solid interface.** Each scalar link into a conducting cell is an anti-bounce-back at
  the temperature that equalises the conductive fluxes of the two half cells, and the heat it
  takes is given to the solid cell: continuity of temperature and flux with heat conserved to
  round-off. The link's heat is affine in the solid's temperature, so the solid takes it
  point-implicitly (an explicit link needs `Σ 2 w_q k_s / (k_s + k_f) < ρ c / (ρ c)_f`, which a
  corner cell of aluminium in water breaks). Interface conditions for the thermal lattice equation, conjugate ones among them:
  L. Li, R. Mei, J. F. Klausner, *Boundary conditions for thermal lattice Boltzmann equation
  method*, J. Comput. Phys. 237:366–395, 2013.
- **Thermal contacts.** A paste, pad or gap filler between two conducting solids is not
  resolved: its resistance per area adds in series with the two half cells on every face
  between the bodies' cells, `g = 1 / (1/(2 k_a) + 1/(2 k_b) + R (ρ c)_f Δx / Δt)` in lattice
  units (the contact resistance of a composite wall's series circuit: F. P. Incropera,
  D. P. DeWitt, T. L. Bergman, A. S. Lavine, *Fundamentals of Heat and Mass Transfer*, 7th ed.,
  Wiley 2011, §3.1.4). Both cells see the same conductance, so heat is conserved, and a
  resistance only lowers it, so the sub-cycling limit holds. Between a solid and the fluid
  there is none yet.
- **Limits.** Staircase interfaces (cell faces, not the true surface): curved interfaces lose
  accuracy, which cut cells with the flux on the true surface recover; a steady run scales the
  solids' heat capacity, so its transient is not physical.

---

## Part D: Fluid–structure interaction

### D1. Immersed boundary

**D1.1 Peskin 2002.** **Verdict: build on.**
- **Citation.** C. S. Peskin. *The immersed boundary method.* Acta Numerica 11:479–517, 2002. DOI [10.1017/S0962492902000077](https://doi.org/10.1017/S0962492902000077).
- **Method.** The structure is a Lagrangian set of points on a fixed Eulerian fluid grid; regularised delta functions
  (2-, 3-, 4-point kernels with stated moment conditions) spread structure forces to the grid and interpolate fluid
  velocity back. No body-fitted fluid mesh, so arbitrarily large motion is free.
- **Limits.** The interface is smeared over ~2 cells (first-order near the wall); force is computed from the
  structure's constitutive law (classic) which brings a stiff time-step restriction for stiff structures.
- **Code.** IBAMR (BSD-3-style): READ OK. IB2d (GPL-3.0): FORBIDDEN.

**D1.2 Mittal & Iaccarino 2005; Uhlmann 2005.** **Verdict: reference / build on.**
- R. Mittal, G. Iaccarino. *Immersed boundary methods.* Annu. Rev. Fluid Mech. 37:239–261, 2005. DOI [10.1146/annurev.fluid.37.061903.175743](https://doi.org/10.1146/annurev.fluid.37.061903.175743). (Review: continuous vs discrete forcing, sharp-interface methods.)
- M. Uhlmann. *An immersed boundary method with direct forcing for the simulation of particulate flows.* JCP 209(2):448–476, 2005. DOI [10.1016/j.jcp.2005.03.017](https://doi.org/10.1016/j.jcp.2005.03.017). Direct forcing: the force is whatever makes the interpolated fluid velocity equal the structure velocity, so no stiff spring. **Build on.**

**D1.3 Feng & Michaelides 2004; Wu & Shu 2009 — IB-LBM.** **Verdict: build on (Wu–Shu).**
- Z.-G. Feng, E. E. Michaelides. *The immersed boundary-lattice Boltzmann method for solving fluid–particles interaction problems.* JCP 195(2):602–628, 2004. DOI [10.1016/j.jcp.2003.10.013](https://doi.org/10.1016/j.jcp.2003.10.013).
- J. Wu, C. Shu. *Implicit velocity correction-based immersed boundary-lattice Boltzmann method and its applications.* JCP 228(6):1963–1979, 2009. DOI [10.1016/j.jcp.2008.11.019](https://doi.org/10.1016/j.jcp.2008.11.019).
- **Method.** Feng–Michaelides: penalty (feedback) forcing in LBM via a body-force term. Wu–Shu: treat the velocity
  correction at all Lagrangian points implicitly — one small linear system per step (matrix of delta-kernel
  overlaps) — so the no-slip condition is satisfied exactly at the markers with no penetration.
- **Limits.** Wu–Shu solves a dense-ish system in the number of markers each step (factor once for rigid bodies,
  re-solve for deforming ones; use CG). Smeared interface, so wall shear is first-order.

**D1.4 Momentum-exchange bounce-back for moving bodies (Ladd C3.1, Bouzidi C3.2)** is the sharp-interface
alternative to IB in LBM: cells are uncovered/covered as the body moves (refill rule needed). **Verdict: reference** —
IB is simpler with a deforming FE surface; bounce-back is better for rigid or slowly deforming walls.

### D2. LBM coupled to finite elements

**D2.1 Kollmannsberger et al. 2009.** **Verdict: reference.** S. Kollmannsberger, S. Geller, A. Düster, J. Tölke, C. Sorger, M. Krafczyk, E. Rank. *Fixed-grid fluid–structure interaction in two dimensions based on a partitioned Lattice Boltzmann and p-FEM approach.* IJNME 79(7):817–845, 2009. DOI [10.1002/nme.2581](https://doi.org/10.1002/nme.2581). Partitioned implicit coupling (structure solved implicitly with high-order p-FEM) with interpolated bounce-back on a fixed LBM grid; includes Turek–Hron-type validation.

**D2.2 De Rosis et al. 2013/2014.** **Verdict: reference.** A. De Rosis, G. Falcucci, S. Ubertini, F. Ubertini. *A coupled lattice Boltzmann-finite element approach for two-dimensional fluid–structure interaction.* Comput. Fluids 86:558–568, 2013. DOI [10.1016/j.compfluid.2013.08.004](https://doi.org/10.1016/j.compfluid.2013.08.004). Explicit-in-time staggered LBM–FE coupling for flexible plates.

**D2.3 Li & Favier 2017.** **Verdict: build on (closest to our design).** Z. Li, J. Favier. *A non-staggered coupling of finite element and lattice Boltzmann methods via an immersed boundary scheme for fluid-structure interaction.* Comput. Fluids 143:90–102, 2017. DOI [10.1016/j.compfluid.2016.11.008](https://doi.org/10.1016/j.compfluid.2016.11.008). IB-LBM fluid coupled to a nonlinear FE structure with synchronous time stepping; validated on flexible-flag benchmarks.

**D2.4 Owen, Leonardi, Feng 2011.** **Verdict: reference.** IJNME 87(1–5):66–95, 2011 (online 2010). DOI [10.1002/nme.2985](https://doi.org/10.1002/nme.2985). LBM with immersed moving boundaries (partial-volume bounce-back) coupled to explicit discrete/finite elements — evidence that explicit structure + LBM with subcycling works.

### D3. Partitioned coupling and stability

**D3.1 Causin, Gerbeau, Nobile 2005, added mass.** **Verdict: build on (as the design constraint).**
- **Citation.** P. Causin, J.-F. Gerbeau, F. Nobile. *Added-mass effect in the design of partitioned algorithms for fluid–structure problems.* CMAME 194(42–44):4506–4527, 2005. DOI [10.1016/j.cma.2004.12.005](https://doi.org/10.1016/j.cma.2004.12.005).
- **Result.** For incompressible flow, explicit (loosely coupled) partitioned schemes are unconditionally unstable when
  the structure-to-fluid mass ratio falls below a geometry-dependent threshold — and refining Δt makes it *worse*.
  Strongly coupled (iterated) schemes converge only with relaxation whose required strength grows as the mass ratio
  drops.
- **Consequence.** Light structures in water (membranes, valves, Turek–Hron FSI3 with ρ_s/ρ_f = 1) need strong
  coupling with IQN (D3.3). Heavy structures (steel plates in water, anything in air) work loosely coupled.
- Also: C. Förster, W. A. Wall, E. Ramm, CMAME 196(7):1278–1293, 2007, DOI [10.1016/j.cma.2006.09.002](https://doi.org/10.1016/j.cma.2006.09.002) — the same instability for sequential staggered coupling, and why predictors don't cure it. **Reference.**
- Serino, Banks, Henshaw, Schwendeman 2019, *A stable added-mass partitioned (AMP) algorithm for elastic solids and incompressible flow*, JCP 399:108923, DOI [10.1016/j.jcp.2019.108923](https://doi.org/10.1016/j.jcp.2019.108923): Robin-type interface conditions derived from the solid's characteristics give a stable loosely coupled scheme even for light solids. **Reference (research; promising for explicit solids, whose characteristic impedance is cheap).**

**D3.2 Küttler & Wall 2008, Aitken relaxation.** **Verdict: build on (first strong coupling).**
- **Citation.** U. Küttler, W. A. Wall. *Fixed-point fluid–structure interaction solvers with dynamic relaxation.* Comput. Mech. 43:61–72, 2008. DOI [10.1007/s00466-008-0255-5](https://doi.org/10.1007/s00466-008-0255-5).
- **Method.** Fixed-point (Dirichlet–Neumann) iteration on the interface displacement with Aitken's Δ² dynamic
  relaxation factor computed from the last two residuals; also steepest-descent variants.
- **Limits.** Convergence degrades with strong added mass; still the cheapest robust starting point.
- **Complexity.** Low (≈100 lines).
- **Code.** 4C (LGPL-3.0): FORBIDDEN.

**D3.3 Degroote, Bathe, Vierendeels 2009, IQN-ILS.** **Verdict: build on (production strong coupling).**
- **Citation.** J. Degroote, K.-J. Bathe, J. Vierendeels. *Performance of a new partitioned procedure versus a monolithic procedure in fluid–structure interaction.* Comput. Struct. 87(11–12):793–801, 2009. DOI [10.1016/j.compstruc.2008.11.013](https://doi.org/10.1016/j.compstruc.2008.11.013).
- **Method.** Interface quasi-Newton with an approximate inverse Jacobian from least squares on the differences of
  previous residuals and outputs (black-box: needs only the two solvers' input/output); within ~5–10 iterations even
  for strong added mass; competitive with monolithic.
- **Limits.** The least-squares system becomes ill-conditioned as columns accumulate: filter (Haelterman et al. 2016,
  Comput. Struct. 171:9–17, DOI [10.1016/j.compstruc.2016.04.001](https://doi.org/10.1016/j.compstruc.2016.04.001), **build on**) and reuse columns from earlier steps with care.
- **Complexity.** Low–moderate (QR updates on a tall skinny matrix).
- **Code.** preCICE (LGPL-3.0) FORBIDDEN; Kratos core (BSD-style) READ OK — check the CoSimulation application's own licence file first.

**D3.4 Farhat et al. 2010, explicit–explicit staggered for compressible FSI.** **Verdict: reference (the explicit–explicit path).**
- **Citation.** C. Farhat, A. Rallu, K. Wang, T. Belytschko. *Robust and provably second-order explicit–explicit and implicit–explicit staggered time-integrators for highly non-linear compressible fluid–structure interaction problems.* IJNME 84(1):73–107, 2010. DOI [10.1002/nme.2883](https://doi.org/10.1002/nme.2883).
- **Method.** Staggered schemes with a carefully placed fluid/structure time offset that are second-order and
  stable for compressible flow + explicit structure (the added-mass problem is far milder for compressible fluids).
- **Use.** LBM is a weakly compressible explicit scheme; this is the closest analysis for explicit LBM + explicit FE.

### D4. ALE

**D4.1 Donea et al. 2004; Hughes, Liu, Zimmermann 1981.** **Verdict: reference (not first).**
- J. Donea, A. Huerta, J.-Ph. Ponthot, A. Rodríguez-Ferran. *Arbitrary Lagrangian–Eulerian methods.* Encyclopedia of Computational Mechanics, Wiley 2004. DOI [10.1002/0470091355.ecm009](https://doi.org/10.1002/0470091355.ecm009).
- T. J. R. Hughes, W. K. Liu, T. K. Zimmermann. *Lagrangian-Eulerian finite element formulation for incompressible viscous flows.* CMAME 29(3):329–349, 1981. DOI [10.1016/0045-7825(81)90049-9](https://doi.org/10.1016/0045-7825(81)90049-9).
- **Method.** The mesh moves with an arbitrary velocity: Lagrangian at the FSI interface, Eulerian far away, with mesh
  smoothing and the convective term relative to mesh motion; geometric conservation law.
- **Why not first.** LBM is fixed-grid; ALE needs a body-fitted, moving fluid mesh and a FV/FE fluid solver. ALE in the
  *explicit solver* (A9.3) is a different use (multi-material hydro) and comes after Lagrangian.

### D5. Interface data transfer

**D5.1 de Boer, van Zuijlen, Bijl 2007; Farhat, Lesoinne, Le Tallec 1998.** **Verdict: build on.**
- A. de Boer, A. H. van Zuijlen, H. Bijl. *Review of coupling methods for non-matching meshes.* CMAME 196(8):1515–1525, 2007. DOI [10.1016/j.cma.2006.03.017](https://doi.org/10.1016/j.cma.2006.03.017).
- C. Farhat, M. Lesoinne, P. Le Tallec. *Load and motion transfer algorithms for fluid/structure interaction problems with non-matching discrete interfaces.* CMAME 157(1–2):95–114, 1998. DOI [10.1016/S0045-7825(97)00216-8](https://doi.org/10.1016/S0045-7825(97)00216-8).
- **Method.** Nearest neighbour, nearest projection (project each fluid point onto the structural surface and use the
  element's shape functions), and RBF; conservative transfer: if displacements map with H (`u_f = H u_s`) then forces
  map with Hᵀ (`f_s = Hᵀ f_f`), which conserves energy over the interface.
- **Use.** With IB-LBM the Lagrangian markers *are* the structural surface nodes (or quadrature points of its faces),
  so mapping reduces to the IB kernel plus Hᵀ. Projection mapping is needed only for body-fitted coupling.

**D5.2 Beckert & Wendland 2001, RBF.** **Verdict: reference (needed only for non-IB couplings).** Aerosp. Sci. Technol. 5(2):125–134, 2001. DOI [10.1016/S1270-9638(00)01087-7](https://doi.org/10.1016/S1270-9638(00)01087-7). Compactly supported Wendland RBF interpolation for mesh-free interface mapping.

### D6. Verification

| Case | Reference | Expected |
|---|---|---|
| Turek–Hron FSI1/FSI2/FSI3 | S. Turek, J. Hron, *Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow*, in *Fluid-Structure Interaction*, LNCSE 53:371–385, Springer 2006, DOI [10.1007/3-540-34596-5_15](https://doi.org/10.1007/3-540-34596-5_15) | Tip displacement amplitude/frequency, drag/lift (published numerical reference values). FSI1 steady; FSI2 large periodic motion (ρ_s/ρ_f = 10); FSI3 strong added mass (ρ_s/ρ_f = 1) — the IQN test. |
| Rigid cylinder in a fluid, added-mass coefficient | potential flow (analytic) | C_a = 1 for a circular cylinder; tests the coupling's inertia, not just loads |
| 1D piston / elastic slab on a fluid column | analytic | linear acoustic–structure solution; checks energy conservation of the transfer |
| Sphere settling (terminal velocity, Stokes regime) | analytic | `U = 2(ρ_s − ρ_f) g R² / (9μ)` with wall correction |

---

## Part E: Alternatives and GPU compute

**E1. Finite-volume CFD (SIMPLE/PISO).** **Verdict: reference (not first).**
- S. V. Patankar, D. B. Spalding. *A calculation procedure for heat, mass and momentum transfer in three-dimensional parabolic flows.* Int. J. Heat Mass Transf. 15(10):1787–1806, 1972. DOI [10.1016/0017-9310(72)90054-3](https://doi.org/10.1016/0017-9310(72)90054-3).
- R. I. Issa. *Solution of the implicitly discretised fluid flow equations by operator-splitting.* JCP 62(1):40–65, 1986. DOI [10.1016/0021-9991(86)90099-9](https://doi.org/10.1016/0021-9991(86)90099-9).
- C. M. Rhie, W. L. Chow. AIAA J. 21(11):1525–1532, 1983. DOI [10.2514/3.8284](https://doi.org/10.2514/3.8284). (Collocated-grid pressure–velocity coupling.)
- J. H. Ferziger, M. Perić, R. L. Street. *Computational Methods for Fluid Dynamics*, 4th ed., Springer 2020. DOI [10.1007/978-3-319-99693-6](https://doi.org/10.1007/978-3-319-99693-6).
- **Why LBM first anyway.** FV needs a body-fitted (polyhedral or hex-dominant) fluid mesh, a pressure Poisson solve
  (algebraic multigrid in pure Rust: real work) and turbulence-model wall functions tuned to y+; its robustness lives
  in the mesh quality, which is exactly what our meshers will be judged on. LBM needs only a voxel/octree grid from the
  surface (no fluid volume mesh to fail), is explicit and local (trivially parallel, GPU-friendly, wasm-friendly), and
  its FSI with IB needs no mesh motion. LBM is the right first fluid solver for incompressible, low-Mach, transient flow
  around complex geometry. FV is the right choice for steady RANS at high Re, compressible high-Mach flow and heat
  transfer with conjugate solids, none of which LBM does well. **If wrong:** if owners mostly want steady internal
  flows with pressure drop at high Re (pipes, manifolds), LBM spends hours of transient time-stepping to reach a steady
  state that a SIMPLE solver reaches in minutes; that is the trigger to add FV.

**E2. SPH for fluids and impact.** **Verdict: reference (a later, separate path for fragmentation and blast).**
- J. J. Monaghan. *Smoothed particle hydrodynamics.* Rep. Prog. Phys. 68(8):1703–1759, 2005. DOI [10.1088/0034-4885/68/8/R01](https://doi.org/10.1088/0034-4885/68/8/R01).
- L. D. Libersky, A. G. Petschek, T. C. Carney, J. R. Hipp, F. A. Allahdadi. *High Strain Lagrangian Hydrodynamics.* JCP 109(1):67–75, 1993. DOI [10.1006/jcph.1993.1199](https://doi.org/10.1006/jcph.1993.1199). (SPH with material strength.)
- J. W. Swegle, D. L. Hicks, S. W. Attaway. *Smoothed Particle Hydrodynamics Stability Analysis.* JCP 116(1):123–134, 1995. DOI [10.1006/jcph.1995.1010](https://doi.org/10.1006/jcph.1995.1010). (The tensile instability.)
- G. R. Johnson, R. A. Stryk, S. R. Beissel. *SPH for high velocity impact computations.* CMAME 139:347–373, 1996. DOI [10.1016/S0045-7825(96)01089-4](https://doi.org/10.1016/S0045-7825(96)01089-4).
- S. W. Attaway, M. W. Heinstein, J. W. Swegle. *Coupling of smooth particle hydrodynamics with the finite element method.* Nucl. Eng. Des. 150:199–205, 1994. DOI [10.1016/0029-5493(94)90136-8](https://doi.org/10.1016/0029-5493(94)90136-8).
- C. Antoci, M. Gallati, S. Sibilla. *Numerical simulation of fluid–structure interaction by SPH.* Comput. Struct. 85(11–14):879–890, 2007. DOI [10.1016/j.compstruc.2007.01.002](https://doi.org/10.1016/j.compstruc.2007.01.002).
- **Assessment.** SPH shares the explicit solver's whole stack (time step, EOS, Johnson–Cook, artificial viscosity,
  contact via the same penalty), so it slots *into* `fieldcraft-explicit` as another discretisation — the natural
  home for hypervelocity impact, fragmentation and blast with JWL, which LBM cannot do. As a general fluid solver it is
  noisy (pressure), needs many particles for accuracy and has the tensile instability; LBM beats it for viscous
  low-Mach flow. Code: SPHinXsys (Apache-2.0), SPlisHSPlasH (MIT), PySPH (BSD-3-style), Cabana (BSD-3): READ OK;
  DualSPHysics (LGPL-2.1): FORBIDDEN.

**E3. GPU compute in pure Rust.** **Verdict: CPU first for both solvers; wgpu compute for LBM second; explicit FE on GPU only on evidence.**
- **Stack.** `wgpu` (MIT OR Apache-2.0) compute pipelines with WGSL shaders. The API we call is safe Rust, so
  `unsafe_code = "forbid"` in our crates holds (wgpu's own internals use unsafe, as every dependency on the OS does).
  `rust-gpu` (Rust → SPIR-V) would let kernels be written in Rust, but it needs a pinned nightly toolchain to build
  shaders: reject for now, revisit when it builds on stable.
- **Precision.** WGSL (https://www.w3.org/TR/WGSL/, 200) has no `f64` type; wgpu exposes `Features::SHADER_F64`
  (https://docs.rs/wgpu/latest/wgpu/struct.Features.html, 200) only on native backends that support it, and not on
  Metal or the browser. **So GPU kernels are f32.** For LBM that is standard (store `f_i − w_i` to recover
  precision). For explicit FE, f32 is workable for stresses and velocities (single precision is common in explicit
  codes) but coordinates/displacements over long runs need f64 accumulation; that rules out a naive GPU port.
- **Layering.** `xtask/src/layers.rs` lists `wgpu` in `UI_CRATES`: nothing below L5 may depend on it. A GPU backend
  therefore needs either an owner decision to allow `wgpu` in a dedicated compute crate below L5 (headless compute is
  not UI), or a `fieldcraft-gpu` crate at L5 that implements a backend trait the L1 solver crates define. I recommend
  the second (no rule change): the solver crate defines a `LatticeBackend` trait and a CPU implementation; the GPU
  implementation lives at L5. **Owner decision either way.**
- **Determinism.** The project requires bit-identical output on every platform. GPU f32 results differ from CPU and
  between vendors (FMA contraction, transcendental precision). LBM kernels are per-cell with no reductions, so a GPU
  run is deterministic on one device, but not equal to the CPU. Explicit FE assembly on GPU needs scatter-add with
  float atomics, which is order-dependent (non-deterministic even on one device) unless assembly is a gather over a
  node-to-element list. **The rule needs an explicit exception: CPU path is the bit-identical reference; GPU results
  are deterministic per device and verified against CPU to a tolerance.** Owner decision.
- **CPU parallelism.** `rayon` (MIT OR Apache-2.0) over cell blocks (LBM) and element blocks with a deterministic
  gather assembly (explicit FE) is deterministic by construction. On `wasm32-unknown-unknown` there are no threads by
  default: rayon must be behind a feature, with a serial fallback (same results).
- **Browser.** WebGPU's default limits (128 MiB per storage binding, 256 MiB per buffer, unless the adapter grants
  more) cap a browser LBM domain at roughly 1.5 M D3Q19 cells per buffer unless populations are split across buffers.
  Fine for demos, not for production. **(limits from memory of the WebGPU spec; check the adapter's reported limits at
  run time, never assume)**
- **Feasibility verdict.** LBM on wgpu: high feasibility, 10–30× over a multi-core CPU on a discrete GPU since it is
  bandwidth-bound and kernels are simple. Explicit FE on wgpu: feasible for hex8 one-point + J2 (a single kernel) but
  contact search, erosion and shells make it a large porting effort; LS-DYNA-class codes remain CPU/MPI-first for good
  reasons. Do it last, and only if profiling shows the element loop dominates.

---

## Part F: Plan (opinionated)

### F1. Strategic call

Build **explicit FE first, LBM second, FSI third.** Explicit FE reuses the M6 element and material work, runs on the
tet meshes M2 produces (Bonet–Burton) before hexes exist, writes the same LS-DYNA keyword concepts M5 targets, and
its verification is analytic. LBM is independent of everything except the octree (M3) and is the best place to spend
GPU effort. FSI needs both to be verified first; a coupling bug on top of an unverified solver is undiagnosable.

One correction to the brief: **LBM + EOS is not a pairing.** LBM is a weakly compressible, low-Mach method; EOS
matters for shocks, detonations and blast, which LBM cannot represent. EOS-driven fluids (water under underwater
explosion, detonation products, air blast) belong to the explicit solver (Lagrangian first, then ALE/Eulerian per
Benson 1992 or SPH). The LBM–explicit FSI covers low-Mach problems: flexible structures in flow, valves, sloshing,
slamming with free-surface LBM. If the owner's motivating use cases are blast/impact on structures, the order changes:
explicit FE → EOS → SPH or ALE multi-material, and LBM drops a step.

### F2. Crate layout

Current layers: L0 `model` → L1 `geometry`, `mesh`, `solver`, `decks`, `render` → L2 `doc` → L3 `io` → L4 `engine` → L5.
Implicit and explicit share element kinematics (shape functions, quadrature, B-matrices, hex8/tet4/quad4), materials
(elastic, J2, Johnson–Cook) and EOS. Two L1 crates cannot depend on each other, and FSI depends on both solvers. So:

| Layer | Crate | Contents |
|---|---|---|
| L0 | `fieldcraft-model` | + materials, EOS, sections/formulations, explicit and fluid steps, contact, fluid domains, coupling interfaces (F3) |
| L1 | `fieldcraft-fem` (new) | element kernels (shape functions, quadrature, B/B-bar, one-point hex with FB hourglass and Belytschko–Bindeman stabilisation, Bonet–Burton tet, BLT shell), constitutive updates (hypoelastic, objective rates, J2 radial return, Johnson–Cook, damage), EOS evaluation. Pure functions over `model` types; no solver loop. |
| L2 | `fieldcraft-lbm` (new) | lattices, collision operators, BCs, voxelisation from the face-tagged surface and the octree from `mesh`\*, stepping on CPU, momentum-exchange forces, IB kernels; a backend trait for GPU |
| L2 | `fieldcraft-solver` (moved up from L1) | implicit linear static/modal (M6), later nonlinear implicit; uses `fem` |
| L2 | `fieldcraft-explicit` (new) | central-difference driver, time step, mass scaling, q, contact (search + penalty), erosion, energy accounting, output; uses `fem` (and later SPH) |
| L3 | `fieldcraft-fsi` (new) | IB coupling, interface mapping, Aitken/IQN-ILS, subcycling, coupled step driver; uses `explicit`, `solver`, `lbm` |
| L4… | `doc`, `io`, `engine`, then L5+ shifts up by two | results (time histories, field frames) in `doc`/`io` |
| L5+ | `fieldcraft-gpu` (new, with the UI layer or an owner-approved compute layer) | wgpu implementations of the solver backend traits |

\* `fieldcraft-lbm` needs the octree and predicates from `fieldcraft-mesh` (L1), so it sits at L2. **`lbm` at L2** next to `solver` and `explicit` (it may then use `mesh`'s octree and
predicates), and `fem` at L1. That is one coherent rule: L1 = geometry and math kernels, L2 = solvers, L3 = coupling.

Cost of this being wrong: renumbering the layer table is a mechanical edit to `xtask/src/layers.rs` and
`docs/architecture.md`; the risk is churn, not breakage. The alternative (one `fieldcraft-solver` crate with
`implicit`/`explicit`/`lbm`/`fsi` modules) needs no layer change, but a 50k-line crate compiles slowly and lets the
LBM depend on FE internals by accident. If the owner prefers no renumbering, take the single-crate option and enforce
module boundaries with `pub(crate)` discipline. This layout change is an owner decision.

### F3. Solver-neutral model extensions in `fieldcraft-model`

- **Materials** (named, unit-consistent, no product numbering): density; isotropic elastic; J2 plasticity with
  isotropic/kinematic hardening (bilinear, power law, table); Johnson–Cook strength (A, B, n, C, m, ε̇₀, T_room,
  T_melt, specific heat, Taylor–Quinney); damage (JC D1–D5); failure/erosion criteria; Newtonian fluid (ρ, μ).
- **EOS** as a separate attachable record: `ideal_gas{γ}`, `polynomial{c0..c6}`, `mie_gruneisen{c0, s1[,s2,s3], γ0, a}`,
  `jwl{A, B, R1, R2, ω, e0, D, p_cj}`, `tillotson{…}`, `tabulated{…}`; a material with an EOS takes its pressure from
  it and only the deviatoric part from the strength model.
- **Sections / formulations:** solid formulation enum (`one_point_hourglass{kind, coeff}`, `assumed_strain`,
  `selective_reduced`, `full`, `nodal_pressure_tet`), shell sections (thickness, integration points, BLT/Hughes–Liu).
- **Steps:** `explicit_dynamic{end_time, dt_scale (default 0.9), mass_scaling: none | target_dt, output_interval}`,
  `fluid{end_time or steps, lattice, collision, les}`, `coupled{explicit/implicit structure, fluid, scheme: loose |
  aitken | iqn_ils, tolerances, max_iters}`. Every iterative or time-marching step carries caps (max steps, max
  iterations, max cells) per the never-crash rules.
- **Loads/conditions:** initial velocity, prescribed motion, rigid walls, gravity, detonation points with lighting
  times.
- **Contact:** surface pairs and self-contact on face sets, penalty scale, friction, search interval.
- **Fluid domain:** bounding box, base spacing, refinement regions (by distance to named faces), boundary conditions on
  box faces (velocity inlet, pressure outlet, slip/no-slip) and on named geometry faces (wall, moving wall); fluid
  material.
- **Coupling interface:** named structural face set ↔ fluid domain.
- **Output requests:** history (energies, reaction forces, contact forces, probe points) and field frames.

All of these must map to the LS-DYNA keyword and Abaqus input concepts M5 writes, described in our own words.

### F4. Milestones

| # | Milestone | Needs from meshing | Verification (analytic first) |
|---|---|---|---|
| X0 | `fieldcraft-fem` extraction with M6: hex8/tet4 kinematics, elastic and J2 updates, shared with the implicit solver | M2 tets, generated boxes | patch tests, J2 uniaxial curve, rigid rotation |
| X1 | Explicit core: central difference, lumped mass, element Δt (FB 1984), one-point hex8 + FB hourglass, Bonet–Burton tet4, artificial viscosity, energy accounting, prescribed BCs, initial velocity, history output | structured hex8 boxes generated in code (no M3 needed), M2 tets | bar wave (front position and amplitude), bar–rigid wall contact time 2L/c, energy balance, cantilever period vs M6 modal |
| X2 | Large deformation materials: Hughes–Winget Jaumann, Flanagan–Taylor Green–Naghdi, J2 radial return, Johnson–Cook strength + damage, erosion | same | Dienes simple shear (both rates, closed form), Taylor bar (Taylor 1948 analytic estimate; experimental comparison only if the owner allows) |
| X3 | Contact: BVH global search, node-to-segment penalty two-pass + friction, pinball for eroding surfaces, rigid walls | surface extraction of volume meshes (have it) | two bars colliding (momentum exchange, analytic), block sliding with friction on a slope (analytic), Hertz sphere-on-plane at low speed (analytic contact time and force) |
| X4 | EOS: polynomial, ideal gas, Mie–Grüneisen, JWL + programmed burn, energy coupling | hex8 strongly preferred (tets smear shocks) — M3 octree hex | 1D shock in Mie–Grüneisen material vs Rankine–Hugoniot with `U_s = c0 + s u_p`, Sod, Noh, Sedov (Lagrangian 1D/axisymmetric) |
| X5 | Shells: BLT + Belytschko–Wong–Chiang warping, self-contact (Benson–Hallquist) | quad4 surface meshes (a quad-dominant surface mesher is not yet on the roadmap — gap) | plate free vibration (analytic), cylinder crush (energy only) |
| L1 | LBM core on CPU: D3Q19, TRT, halfway bounce-back, Zou–He, uniform grid from voxelising the face-tagged surface with exact predicates, momentum-exchange forces | voxelisation (inside/outside with `mesh` predicates); no volume mesh | Poiseuille (2nd order), Couette, Taylor–Green decay, sphere Stokes drag |
| L2 | Accuracy and range: Bouzidi interpolated walls (link fractions from ray–triangle tests), regularised/RR and cumulant D3Q27, Smagorinsky LES, Schäfer–Turek, cavity | link-surface intersections | Schäfer–Turek 2D-1/2D-2/3D (published numerical bounds, owner sign-off), Womersley |
| L3 | Octree refinement: Rohde cell-centred on the M3 2:1 octree; then wgpu backend at L5 (AA pattern) | M3 octree (shared data structure) | refinement across a channel (Poiseuille exact across levels, mass conservation to round-off), GPU-vs-CPU tolerance check |
| F1 | FSI loose: IB-LBM (Wu–Shu) with markers on the structure's surface nodes, explicit FE structure with subcycling, conservative Hᵀ force transfer | structural surface of a hex/tet mesh | rigid cylinder added-mass C_a = 1, settling sphere, Turek–Hron FSI1 |
| F2 | FSI strong: Aitken, then IQN-ILS with filtering; implicit (M6 linear dynamic) structure as well as explicit | same | Turek–Hron FSI2 and FSI3 |
| F3 | Free-surface LBM (Körner/Thürey) with FSI | same | dam break front position (published experiment — owner sign-off), hydrostatic column (analytic) |

**Time-step mismatch, concretely.** Steel structure (c ≈ 5,900 m/s) with 1 mm elements: Δt_s ≈ 1.7×10⁻⁷ s. Water LBM
at 1 mm spacing with lattice velocity 0.1 for a 1 m/s flow: Δt_f = 0.1·Δx/u = 10⁻⁴ s. The structure needs ~600
substeps per fluid step: subcycle the structure (cheap, it is a thin interface region) or apply mass scaling with
reported added mass. For soft structures (Turek–Hron), the implicit structure solver with one step per fluid step is
the better pairing — this is why FSI should support both structure solvers.

### F5. Biggest risks, and what happens if wrong

1. **Hourglass and contact robustness, not the time integrator, decide whether people trust the explicit solver.**
   If hourglass control is weak, coarse-mesh results look plausible and are wrong; if contact misses penetrations,
   runs explode at hour 3. Mitigation: hourglass and contact energies reported every run with thresholds that warn;
   contact verification cases from day one; a crash-report replay path for solver blow-ups.
2. **Silent mass scaling and time step collapse.** One bad element sets Δt for the whole model. Mitigation: report the
   ten controlling elements and the added mass; mesh-quality gate before a run. If wrong, users see runs that never
   finish rather than wrong answers — annoying, not dangerous.
3. **Determinism vs GPU vs parallelism.** The determinism rule and GPU f32 are in direct conflict (E3). If no
   exception is agreed, the GPU work cannot ship; if the exception is too loose, regression tests flake. Decide before
   L3.
4. **Added-mass instability in FSI.** Loose coupling will diverge for light structures in water no matter how small
   the step (Causin 2005). If F1 ships only loose coupling, the first realistic marine/biomedical case fails; F2 (IQN)
   must follow immediately, and the loose scheme must detect divergence (interface energy growth) and stop with an
   error.
5. **LBM at industrial Reynolds numbers.** Laminar validation is easy; turbulent external flow needs D3Q27 cumulant
   or RR, LES, wall models and refinement before results are useful. If underestimated, the LBM looks like a toy for a
   year. Budget L2–L3 as the main LBM work, not L1.
6. **Scope creep toward LS-DYNA's material library.** LS-DYNA has hundreds of materials; each is a paper. Ship
   elastic, J2, Johnson–Cook, a hyperelastic rubber, a crushable foam, and EOS; add on demand. Never mirror any
   product's material numbering in the model — only in the deck writer's mapping.
7. **Expected-value policy.** Taylor bar, cavity, Schäfer–Turek, Turek–Hron, dam break all compare with published
   experiments or published numerical references, not closed forms. The project rule says expected values are
   analytic. **Owner decision:** allow published, non-commercial reference values for benchmarks (I recommend yes,
   tagged as "benchmark", kept separate from the analytic suite).
8. **Quad shell meshing is missing from the roadmap.** X5 needs quad-dominant surface meshes; without them shells
   run on tri3 (BLT has a triangular variant, less accurate). Flag for the roadmap.

## Part G: Implicit solver (M4), as implemented

The papers and books `fieldcraft-fem` and `fieldcraft-solver` are written from (cited again in
the module docs):

- Elements: O. C. Zienkiewicz, R. L. Taylor, J. Z. Zhu, *The Finite Element Method: Its Basis
  and Fundamentals*, 7th ed., 2013; T. J. R. Hughes, *The Finite Element Method*, Dover 2000.
- Incompatible modes: Wilson, Taylor, Doherty, Ghaboussi 1973; R. L. Taylor, P. J. Beresford,
  E. L. Wilson, IJNME 10(6):1211–1219, 1976 (the det J₀/det J correction that passes the patch
  test).
- Tetrahedron quadrature: A. H. Stroud, *Approximate Calculation of Multiple Integrals*, 1971
  (conical product rules).
- Mass lumping: E. Hinton, T. Rock, O. C. Zienkiewicz, Earthquake Eng. Struct. Dyn.
  4(3):245–249, 1976.
- Stress extrapolation: E. Hinton, J. S. Campbell, IJNME 8(3):461–480, 1974.
- Principal stresses (results display): O. K. Smith, *Eigenvalues of a symmetric 3 × 3 matrix*,
  Commun. ACM 4(4):168, 1961 (closed form).
- Ordering: A. George, SIAM J. Numer. Anal. 10(2):345–363, 1973; H. D. Simon, Computing Systems
  in Engineering 2:135–148, 1991.
- Sparse Cholesky: T. A. Davis, *Direct Methods for Sparse Linear Systems*, SIAM 2006; J. W. H.
  Liu, SIAM J. Matrix Anal. Appl. 11(1):134–172, 1990; supernodes: J. W. H. Liu, E. G. Ng, B. W.
  Peyton, SIAM J. Matrix Anal. Appl. 14(1):242–252, 1993; left-looking supernodal factorization:
  E. Ng, B. W. Peyton, SIAM J. Sci. Comput. 14(5):1034–1056, 1993; parallel block factorization:
  E. Rothberg, A. Gupta, SIAM J. Sci. Comput. 15(6):1413–1439, 1994.
- Eigenvalues: T. Ericsson, A. Ruhe, Math. Comp. 35(152):1251–1268, 1980; R. G. Grimes, J. G.
  Lewis, H. D. Simon, SIAM J. Matrix Anal. Appl. 15(1):228–272, 1994; B. N. Parlett, *The
  Symmetric Eigenvalue Problem*, SIAM 1998; G. H. Golub, C. F. Van Loan, *Matrix Computations*,
  4th ed., 2013.
- 2D elements (plane stress, plane strain, axisymmetric): Zienkiewicz, Taylor & Zhu 2013, ch. 6
  and §6.6; Hughes 2000, §2.7 and §3.5. Quad4 incompatible modes (QM6): Wilson et al. 1973 and
  Taylor, Beresford & Wilson 1976 as for hex8; in axisymmetric elements the modes' strain
  operator is shifted by its volume mean, the patch-test correction of E. L. Wilson and A.
  Ibrahimbegovic, *Use of incompatible displacement modes for the calculation of element
  stiffnesses or stresses*, Finite Elements in Analysis and Design 7(3):229–241, 1990. Rotating
  discs and cylinders for verification: S. Timoshenko, J. N. Goodier, *Theory of Elasticity*,
  3rd ed., McGraw-Hill 1970, §32.
- Wedge6: the linear triangle times the linear interval; the symmetric three-point triangle rule
  (Zienkiewicz, Taylor & Zhu 2013, table 6.3) times two-point Gauss–Legendre.
- Couplings and constraint elimination: rigid links and master–slave elimination by
  transformation, R. D. Cook, D. S. Malkus, M. E. Plesha, R. J. Witt, *Concepts and Applications
  of Finite Element Analysis*, 4th ed., Wiley 2002, §9.2, and C. A. Felippa, *Introduction to
  Finite Element Methods*, ch. 8 (University of Colorado lecture notes); automatic choice of the
  dependent degree of freedom, J. F. Abel, M. S. Shephard, IJNME 14(3):464–467, 1979; the
  distributing (interpolation) constraint as a least-squares rigid fit with loads spread by
  virtual work, G. H. K. Heirman, W. Desmet, Multibody Syst. Dyn. 24(2):219–234, 2010; the
  reduced matrix product T_fᵀ K T_f row by row, F. G. Gustavson, ACM TOMS 4(3):250–269, 1978.
- Trusses and discrete elements: the isoparametric bar along its arc, R. D. Cook, D. S. Malkus,
  M. E. Plesha, R. J. Witt, *Concepts and Applications of Finite Element Analysis*, 4th ed.,
  Wiley 2002, §3.6 and §6.2, and Hughes 2000, §1.10; springs k [n nᵀ, −n nᵀ; −n nᵀ, n nᵀ],
  Cook et al. 2002, §2.2 and §2.10; HRZ lumping as above. Verified against the method of joints
  (statically determinate trusses), the three-bar truss, the hanging bar, √(k/m), the two-mass
  chain and the torsional pendulum √(k_t/J).
- Shells: E. N. Dvorkin, K.-J. Bathe, *A continuum mechanics based four-node shell element for
  general nonlinear analysis*, Engineering Computations 1(1):77–88, 1984 (MITC4 shear tying);
  K.-J. Bathe, E. N. Dvorkin, IJNME 21(2):367–383, 1985 (the flat-plate form); P.-S. Lee,
  K.-J. Bathe, *Development of MITC isotropic triangular shell finite elements*, Computers &
  Structures 82(11–13):945–962, 2004 (MITC3); T. J. R. Hughes, F. Brezzi, *On drilling degrees of
  freedom*, CMAME 72(1):105–121, 1989 (drilling penalty); facet shells, O. C. Zienkiewicz,
  R. L. Taylor, D. D. Fox, *The Finite Element Method for Solid and Structural Mechanics*, 7th
  ed., 2014, ch. 9–10.
- Beams: J. S. Przemieniecki, *Theory of Matrix Structural Analysis*, McGraw-Hill 1968, §5.6;
  Z. Friedman, J. B. Kosmatka, Computers & Structures 47(3):473–481, 1993 (interdependent
  Timoshenko interpolation); G. R. Cowper, J. Appl. Mech. 33(2):335–340, 1966 (shear
  coefficients); S. P. Timoshenko, J. N. Goodier, *Theory of Elasticity*, 3rd ed., 1970, §109
  (torsion of a rectangle); W. D. Pilkey, *Analysis and Design of Elastic Beams*, Wiley 2002
  (thin-walled shear areas and torsion).
- Shell benchmarks (published values, tagged, not the analytic gate): Scordelis–Lo roof 0.3024
  (MacNeal & Harder 1985); pinched cylinder 1.8248 · 10⁻⁵ (T. Belytschko, H. Stolarski, W. K.
  Liu, N. Carpenter, J. S.-J. Ong, CMAME 51(1–3):221–258, 1985). Plate: S. Timoshenko,
  S. Woinowsky-Krieger, *Theory of Plates and Shells*, 2nd ed., 1959, §30 (Navier).
- Verification: R. H. MacNeal, R. L. Harder, *A proposed standard set of problems to test finite
  element accuracy*, Finite Elem. Anal. Des. 1(1):3–20, 1985 (patch tests); Lamé, Kirsch,
  Timoshenko and Euler–Bernoulli closed forms; Saint-Venant's principle (tension through a
  coupling), pure bending σ = M z / I, torsion of a circular tube φ = T L / (G J), τ = T r / J.

## Part H: Explicit solver core (X1), as implemented

The papers and books `fieldcraft-explicit` and `fieldcraft_fem::uniform` are written from (cited
again in the module docs):

- Time integration, energy balance, element eigenvalue inequality: T. Belytschko, W. K. Liu,
  B. Moran, K. Elkhodary, *Nonlinear Finite Elements for Continua and Structures*, 2nd ed., Wiley
  2014 (ch. 6: central differences in the half-step form with variable increments, the energy
  balance check, the critical time step).
- One-point hexahedron and hourglass control: D. P. Flanagan, T. Belytschko, IJNME
  17(5):679–706, 1981 (integrated gradients b, mean strain, the γ hourglass shape vectors
  orthogonal to linear fields, viscous and stiffness resisting forces). Our scaling of the
  coefficient with ρ, c and the element size is our own: κ is the fraction of critical damping
  (viscous) or the squared frequency ratio (stiffness) of a cube's hourglass modes at 2c/h.
- Element frequency bound: in the spirit of D. P. Flanagan, T. Belytschko, J. Appl. Mech.
  51(1):35–40, 1984; we bound uᵀKu ≤ (k/V) λ_max(Σ b bᵀ) ‖u‖² with k = max(3λ + 2μ, 2μ)
  (isotropic) or 2 λ_max(D), which is exact for a cube; the 3 × 3 eigenvalue in closed form after
  O. K. Smith, Commun. ACM 4(4):168, 1961.
- Average nodal pressure tetrahedron: J. Bonet, A. J. Burton, Commun. Numer. Methods Eng.
  14(5):437–449, 1998.
- Bulk viscosity and its time-step correction: J. von Neumann, R. D. Richtmyer, J. Appl. Phys.
  21(3):232–237, 1950; M. L. Wilkins, J. Comput. Phys. 36(3):281–303, 1980.
- Springs and dashpots: lagged-velocity damping and the damped stable increment
  Δt ≤ (2/ω)(√(1 + ξ²) − ξ), Belytschko et al. 2014, §6.2–6.3; Gershgorin's bound on
  M^{-1/2} K M^{-1/2} per degree of freedom. Verified against the damped oscillator's decay e^{−ξωt}
  and period 2π/ω_d.
- Rigid walls as kinematic velocity constraints: T. J. R. Hughes, R. L. Taylor, J. L. Sackman,
  A. Curnier, W. Kanoknukulchai, CMAME 8(3):249–276, 1976.
- Verification: d'Alembert's solution of the 1D wave equation (front at c t, σ = ρ c v), the
  elastic bar impact on a rigid wall (contact 2L/c, rebound at the impact speed; Hughes et al.
  1976), Euler–Bernoulli and the implicit modal solver for the cantilever period. The one-point
  element keeps 1 − 1/n² of the bending energy with n elements through the depth, which predicts
  its period with viscous hourglass control to 0.03%.

## Part I: Lattice Boltzmann core (M12, plan L1), as implemented

`fieldcraft-lbm` is written from these papers and the book (cited again in its module docs). No
lattice Boltzmann package's code was opened.

- Lattice, equilibrium, units, stream–collide, bounce-back, forcing and boundary conditions:
  Krüger et al. 2017 (C1.2); D3Q19 from Qian, d'Humières & Lallemand 1992 (C1.3).
- Collision: two relaxation times, Ginzburg, Verhaeghe & d'Humières 2008 (C2.2). The parameter
  Λ = (τ⁺ − ½)(τ⁻ − ½) is the step's `magic` (default 3/16: walls exactly halfway).
- Body force: Guo, Zheng & Shi 2002 (C1.6), split into the symmetric and antisymmetric parts
  relaxed with τ⁺ and τ⁻.
- Walls: halfway bounce-back, Ladd's moving-wall term and momentum exchange (Ladd 1994, C3.1);
  momentum exchange also after R. Mei, D. Yu, W. Shyy, L.-S. Luo, *Force evaluation in the
  lattice Boltzmann method involving curved geometry*, PRE 65:041203, 2002 (DOI not checked).
- Inlets and outlets: Zou & He 1997 (C3.3), with the D3Q19 transverse-momentum correction of
  M. Hecht, J. Harting, *Implementation of on-site velocity boundary conditions for D3Q19
  lattice Boltzmann simulations*, J. Stat. Mech. 2010, P01018 (DOI not checked); convective
  outlets after Orlanski 1976 (C3.4).
- Stability limits (Mach ≤ 0.3, τ⁺ ≥ 0.505, runtime speed ≤ 0.5): Sterling & Chen 1996 and
  Lallemand & Luo 2000 (C1.5) for the shrinking stable region as τ → ½.
- Voxelisation: Shewchuk's exact `orient2d`/`orient3d` (`fieldcraft_mesh::predicates`) with a
  symbolic perturbation after H. Edelsbrunner, E. P. Mücke, *Simulation of Simplicity*, ACM TOG
  9(1):66–104, 1990; link faces by nearest triangle, C. Ericson, *Real-Time Collision Detection*,
  Morgan Kaufmann 2005, §5.1.5.
- Verification: analytic Poiseuille, Couette and Taylor–Green solutions (C8); the Stokes drag of
  a simple cubic array of spheres, H. Hasimoto, *On the periodic fundamental solutions of the
  Stokes equations and their application to viscous flow past a cubic array of spheres*, JFM
  5(2):317–328, 1959, with the higher terms of A. S. Sangani, A. Acrivos, Int. J. Multiphase Flow
  8(4):343–360, 1982 (DOIs not checked). Benchmark (published numerical reference values, tagged,
  not part of the analytic gate): Schäfer & Turek 1996, 2D-1 (C8).

## Part J: Large deformation and inelastic materials (X2), as implemented

The papers and books `fieldcraft_fem::{rotation, plasticity}` and `fieldcraft-explicit`'s X2
update are written from (cited again in the module docs):

- Updated Lagrangian one-point elements, strain increments on the midpoint configuration:
  T. Belytschko, W. K. Liu, B. Moran, K. Elkhodary, *Nonlinear Finite Elements for Continua and
  Structures*, 2nd ed., Wiley 2014 (ch. 4–6); the closed-form integrated gradients b = ∂V/∂x of
  the hexahedron, D. P. Flanagan, T. Belytschko, IJNME 17(5):679–706, 1981 (we tabulate the 24
  coefficients of the volume form and check them against their exact integration).
- Jaumann rate, incrementally objective: T. J. R. Hughes, J. Winget, IJNME 15(12):1862–1867,
  1980 (Cayley rotation of the midpoint spin increment; a rigid rotation of any size is exact).
- Green–Naghdi rate: D. P. Flanagan, L. M. Taylor, CMAME 62(3):305–320, 1987 (incremental polar
  rotation from V̇ = L V − V Ω; we solve [(tr V) I − V] x = axial(D V − V D) and integrate V with
  the midpoint rule).
- Radial return for J2 with isotropic and linear kinematic (Prager) hardening: J. C. Simo, T. J.
  R. Hughes, *Computational Inelasticity*, Springer 1998, §3.3; R. D. Krieg, D. B. Krieg, J.
  Press. Vessel Technol. 99(4):510–515, 1977.
- Johnson–Cook flow stress (rate clamped at ε̇* ≥ 1, T* in [0, 1], adiabatic heating with a
  Taylor–Quinney fraction) and fracture strain D1–D5 with D = Σ Δε_p/ε_f: G. R. Johnson, W. H.
  Cook 1983 (7th Int. Symp. Ballistics) and Eng. Fract. Mech. 21(1):31–48, 1985. The fracture
  strain is floored at 10⁻⁶ (our choice).
- Erosion: G. R. Johnson, R. A. Stryk, Int. J. Impact Eng. 5:411–421, 1987; T. Belytschko, J. I.
  Lin, ibid. 5:111–127, 1987 (delete, keep nodal mass and momentum, report the energy).
- Verification: J. K. Dienes, Acta Mech. 32:217–232, 1979 (simple shear, closed forms for both
  rates); G. I. Taylor, Proc. R. Soc. A 194:289–299, 1948 (the yield stress from the final length
  and undeformed length of a flat-ended projectile). Published experimental Taylor-bar lengths
  (Wilkins & Guinan 1973) are not used: that needs the owner's sign-off for a tagged benchmark
  suite.

## Part K: Heat transfer and thermal stress, as implemented

- **Conduction.** Galerkin form of ρc Ṫ = ∇·(k∇T) + Q with surface flux and convection:
  O. C. Zienkiewicz, R. L. Taylor, J. Z. Zhu, *The Finite Element Method: Its Basis and
  Fundamentals*, 7th ed., Butterworth-Heinemann 2013, §3.2; T. J. R. Hughes, *The Finite Element
  Method*, Dover 2000, ch. 1. Kernels in `fieldcraft_fem::heat`, solver in `fieldcraft_solver::heat`.
- **Time integration.** The generalized trapezoidal (θ) family, unconditionally stable for
  θ ≥ ½ (Hughes 2000, §8.2): backward Euler (θ = 1, the default) and Crank–Nicolson (θ = ½).
- **Thermal strain.** Initial-strain load ∫ Bᵀ D ε_th dV and stress D(ε − ε_th)
  (Zienkiewicz, Taylor & Zhu 2013, §6.2.3); incompatible modes condensed with the same load.
- **Coupled temperature–displacement.** One-way coupling (temperature strains the solid) solved
  in turn, which is the exact monolithic solution of linear quasi-static thermoelasticity
  without thermoelastic dissipation (B. A. Boley, J. H. Weiner, *Theory of Thermal Stresses*,
  Wiley 1960, ch. 1).
- **Verification** (`crates/solver/tests/thermal.rs`, analytic): steady slab conduction, a slab
  with a heat source and convection or flux (quadratic profile), radial conduction through a
  thick cylinder (T ∝ ln r), the semi-infinite solid T = erfc(x / 2√(κt)) (H. S. Carslaw,
  J. C. Jaeger, *Conduction of Heat in Solids*, 2nd ed., Oxford 1959, §2.4), free and
  constrained thermal expansion (σ = −EαΔT), the bimetal strip (S. Timoshenko, *Analysis of
  bi-metal thermostats*, J. Opt. Soc. Am. 11(3):233–255, 1925), and the coupled step against
  the heat step followed by the static step.

## Part L: Hybrid (mixed u–p) solids, as implemented

- Perturbed-Lagrangian u–p formulation with element-wise condensed discontinuous pressure:
  T. J. R. Hughes, *The Finite Element Method*, Dover 2000, §4.3–4.5 (B-bar and mean
  dilatation as the constant-pressure case); J. C. Nagtegaal, D. M. Parks, J. R. Rice, CMAME
  4(2):153–177, 1974 (locking of displacement elements near incompressibility).
- Pressure spaces and the inf-sup condition: F. Brezzi, M. Fortin, *Mixed and Hybrid Finite
  Element Methods*, Springer 1991, ch. VI; K.-J. Bathe, *Finite Element Procedures*, 2nd ed.,
  2014, §4.5 (Q2–P1disc and P2–P0 stable, Q1–P0 checkerboard-prone, P1–P0 locks).
- Verification (`crates/solver/tests/hybrid.rs`, analytic): a block in uniaxial strain at
  ν = 0.4999 (exact), Lamé's thick cylinder in plane strain at ν = 0.4999 (hex20 1e-5,
  tet10 1e-3, hex8 3e-3 relative error in the bore displacement, where the displacement hex8
  is 80 % too stiff). Cook's membrane is a published benchmark, not part of the analytic gate.

## Part M: Nonlinear static analysis (finite strain), as implemented

- Total Lagrangian formulation written with spatial quantities (Kirchhoff stress, spatial
  tangent, current gradients): J. Bonet, R. D. Wood, *Nonlinear Continuum Mechanics for Finite
  Element Analysis*, 2nd ed., Cambridge 2008, ch. 5–9; G. A. Holzapfel, *Nonlinear Solid
  Mechanics*, Wiley 2000, §6.6 (the isochoric–volumetric split and its tangents); K.-J. Bathe,
  *Finite Element Procedures*, 2nd ed., 2014, §6.3; T. Belytschko, W. K. Liu, B. Moran,
  K. Elkhodary, *Nonlinear Finite Elements for Continua and Structures*, 2nd ed., Wiley 2014,
  ch. 4–6.
- Decoupled neo-Hookean energy: P. J. Flory, Trans. Faraday Soc. 57:829–838, 1961.
- Mixed u–p at finite strain: J. C. Simo, R. L. Taylor, K. S. Pister, CMAME 51:177–208, 1985;
  mean dilatation: J. C. Nagtegaal, D. M. Parks, J. R. Rice, CMAME 4:153–177, 1974.
- Enhanced assumed strain at finite strain (why incompatible modes are refused for now): J. C.
  Simo, F. Armero, IJNME 33:1413–1449, 1992; their hourglassing in compression: P. Wriggers,
  S. Reese, CMAME 135:201–209, 1996.
- Follower pressure and its load stiffness: K. Schweizerhof, E. Ramm, Comput. Struct.
  18(6):1099–1114, 1984.
- Newton–Raphson, line search and convergence criteria: M. A. Crisfield, *Non-linear Finite
  Element Analysis of Solids and Structures*, vol. 1, Wiley 1991, ch. 9; K.-J. Bathe,
  A. P. Cimento, CMAME 22:59–85, 1980.
- Symmetric eigen-decomposition for the logarithmic strain: G. H. Golub, C. F. Van Loan,
  *Matrix Computations*, 4th ed., Johns Hopkins 2013, §8.5 (cyclic Jacobi).
- Verification references: the elastica, K. E. Bisshopp, D. C. Drucker, Q. Appl. Math.
  3(3):272–275, 1945, with K. Mattiasson, IJNME 17(1):145–153, 1981; Cook's membrane (pending its
  transcribed finite-strain values), R. D. Cook, J. Struct. Div. ASCE 100(ST9):1851–1863, 1974,
  and J. C. Simo, F. Armero, IJNME 33:1413–1449, 1992; tube inflation, R. S. Rivlin, Phil. Trans. R. Soc. A
  242:173–195, 1949, and R. W. Ogden, *Non-Linear Elastic Deformations*, Dover 1997.

## Part N: Implicit plasticity (finite and small strain), as implemented

- Multiplicative J2 plasticity with a Hencky elastic energy and the exponential return map in
  the principal logarithmic strains, which keeps the classical radial return of the
  infinitesimal theory: J. C. Simo, *Algorithms for static and dynamic multiplicative plasticity
  that preserve the classical return mapping schemes of the infinitesimal theory*, CMAME
  99(1):61–112, 1992, DOI [10.1016/0045-7825(92)90123-2](https://doi.org/10.1016/0045-7825(92)90123-2);
  J. C. Simo, T. J. R. Hughes, *Computational Inelasticity*, Springer 1998, ch. 9.
- The radial return (shared with explicit dynamics) and its consistent tangent
  K I⊗I + 2Gθ ℙ − 2Gθ̄ n⊗n: Simo & Hughes 1998, §3.3, Box 3.2 (A4.1); J. C. Simo, R. L. Taylor,
  *Consistent tangent operators for rate-independent elastoplasticity*, CMAME 48(1):101–118,
  1985, DOI [10.1016/0045-7825(85)90070-2](https://doi.org/10.1016/0045-7825(85)90070-2).
- The derivative of the matrix logarithm by divided differences (Daleckii–Krein), smooth through
  equal stretches: R. Bhatia, *Matrix Analysis*, Springer 1997, §V.3.
- Kinematic hardening at finite strain (refused for now; the route when it comes): additive
  plasticity in the Lagrangian logarithmic strain space, C. Miehe, N. Apel, M. Lambrecht, CMAME
  191(47–48):5383–5425, 2002, DOI [10.1016/S0045-7825(02)00438-3](https://doi.org/10.1016/S0045-7825(02)00438-3).
- Prager's linear kinematic hardening: W. Prager, *The theory of plasticity: a survey of recent
  achievements*, Proc. Inst. Mech. Eng. 169:41–57, 1955.
- Verification references: the thick cylinder's plastic zone and collapse pressure, R. Hill,
  *The Mathematical Theory of Plasticity*, Oxford 1950, ch. V; the cantilever's plastic hinge,
  B. G. Neal, *The Plastic Methods of Structural Analysis*, 3rd ed., Chapman & Hall 1977, and the
  continuum's velocity fields, W. Prager, P. G. Hodge, *Theory of Perfectly Plastic Solids*,
  Wiley 1951.

## Part O: Polynomial hyperelasticity and thermal strain in nonlinear steps, as implemented

- Decoupled invariant-based hyperelasticity, the fictitious stress and elasticity and their
  projection: G. A. Holzapfel, *Nonlinear Solid Mechanics*, Wiley 2000, §6.4–6.6 (eqs. 6.168,
  6.193).
- Mooney–Rivlin: M. Mooney, *A theory of large elastic deformation*, J. Appl. Phys.
  11(9):582–592, 1940; R. S. Rivlin, *Large elastic deformations of isotropic materials IV*,
  Phil. Trans. R. Soc. A 241:379–397, 1948. Yeoh: O. H. Yeoh, *Some forms of the strain energy
  function for rubber*, Rubber Chem. Technol. 66(5):754–771, 1993.
- The homogeneous modes of rubber testing: L. R. G. Treloar, *The Physics of Rubber
  Elasticity*, 3rd ed., Oxford 1975, ch. 10; the inflation and extension of a tube: R. S. Rivlin,
  Phil. Trans. R. Soc. A 242:173–195, 1949, and R. W. Ogden, *Non-Linear Elastic Deformations*,
  Dover 1997, §5.3.
- Thermal strain by the multiplicative split F = F_m F_θ: V. A. Lubarda, *Constitutive
  theories based on the multiplicative decomposition of deformation gradient: thermoelasticity,
  elastoplasticity, and biomechanics*, Appl. Mech. Rev. 57(2):95–108, 2004.

## Part P: Symmetric indefinite factorization, linear buckling and arc-length, as implemented

- Symmetric indefinite pivoting: J. R. Bunch, L. Kaufman, *Some stable methods for calculating
  inertia and solving symmetric linear systems*, Math. Comp. 31(137):163–179, 1977; the bounded
  (rook) search, C. Ashcraft, R. G. Grimes, J. G. Lewis, *Accurate symmetric indefinite linear
  equation solvers*, SIAM J. Matrix Anal. Appl. 20(2):513–561, 1998; the panel's partial
  factorization of LAPACK's `dsytrf`, E. Anderson et al., *LAPACK Users' Guide*, 3rd ed., SIAM
  1999.
- Pivoting within supernodes on a static structure: O. Schenk, K. Gärtner, *On fast
  factorization pivoting methods for sparse symmetric indefinite systems*, Electron. Trans.
  Numer. Anal. 23:158–179, 2006. Delayed pivots (the multifrontal alternative, not used: it
  changes the structure as it goes): I. S. Duff, J. K. Reid, *The multifrontal solution of
  indefinite sparse symmetric linear equations*, ACM TOMS 9(3):302–325, 1983; J. W. H. Liu, *The
  multifrontal method for sparse matrix solution: theory and practice*, SIAM Review
  34(1):82–109, 1992.
- Arc-length: E. Riks, *An incremental approach to the solution of snapping and buckling
  problems*, Int. J. Solids Struct. 15(7):529–551, 1979; M. A. Crisfield, *A fast
  incremental/iterative solution procedure that handles "snap-through"*, Computers & Structures
  13:55–62, 1981; E. Ramm, *Strategies for tracing the nonlinear response near limit points*, in
  *Nonlinear Finite Element Analysis in Structural Mechanics*, Springer 1981, 63–89; the
  predictor's sign, M. Feng, D. Perić, D. R. J. Owen, Computers & Structures 58(3):479–485, 1996.
  Snap-through of a shallow sinusoidal arch: Y. C. Fung, A. Kaplan, NACA TN 2840, 1952;
  Timoshenko & Gere 1961, §7.6.
- Iterative refinement: N. J. Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed.,
  SIAM 2002, §12.2.
- Linear buckling and the initial-stress stiffness: K.-J. Bathe, *Finite Element Procedures*,
  2nd ed., 2014, §6.8.2; the buckling spectral transformations in the K inner product: R. G.
  Grimes, J. G. Lewis, H. D. Simon, SIAM J. Matrix Anal. Appl. 15(1):228–272, 1994, §3; Euler's
  columns with Engesser's shear correction: S. P. Timoshenko, J. M. Gere, *Theory of Elastic
  Stability*, 2nd ed., McGraw-Hill 1961, §2.1, §2.17.

## Part Q: Contact and friction (explicit), as implemented

The papers and books `fieldcraft_fem::contact` and `fieldcraft_explicit::contact` are written
from (cited again in the module docs); docs/fem.md describes the formulation.

- Closest point on a triangle by its Voronoi regions: C. Ericson, *Real-Time Collision
  Detection*, Morgan Kaufmann 2005, §5.1.5; on a bilinear quadrilateral by Newton's method on the
  stationarity conditions, projected on the parameter square: the local search of node-to-segment
  contact, J. O. Hallquist, G. L. Goudreau, D. J. Benson, CMAME 51:107–137, 1985 (A5.1), and P.
  Wriggers, *Computational Contact Mechanics*, 2nd ed., Springer 2006, §4.1.
- Global search by bucket sorting into a uniform grid: D. J. Benson, J. O. Hallquist, CMAME
  78:141–163, 1990 (A5.2). We chose a grid over the BVH of A5.4: one pass per increment, memory
  bounded by caps whatever the mesh, candidates in face order (determinism); the cell grows when a
  few huge faces would blow the entry cap.
- Segment-to-segment integration (projection on the secondary face's plane, clipping,
  sub-triangle quadrature): M. A. Puso, T. A. Laursen, *A mortar segment-to-segment contact
  method for large deformation solid mechanics*, CMAME 193:601–629, 2004; with the penalty applied
  at the integration points, K. A. Fischer, P. Wriggers, Comput. Mech. 36:226–244, 2005; polygon
  clipping, I. E. Sutherland, G. W. Hodgman, Commun. ACM 17(1):32–42, 1974. Chosen over
  node-to-segment because the latter fails the contact patch test on non-matching meshes (R. L.
  Taylor, P. Papadopoulos, in *Nonlinear Computational Mechanics*, Springer 1991, 690–702; M. A.
  Crisfield, IJNME 48:435–449, 2000), which our first node-to-segment version showed too (12%
  spread of σ_zz in the patch test, 0.1% now).
- Penalty from the masses and the stable increment (the "soft constraint" idea: J. O. Hallquist,
  *LS-DYNA Theory Manual*, LSTC 2006, the contact chapter's soft-constraint option; the bound is ours: Gershgorin's row sums of the
  integrated contact stiffness with the node's directional area, added to the element
  frequencies by Weyl's inequality).
- Friction by a regularised (penalty) return map: P. Wriggers, *Computational Contact
  Mechanics*, §5.2 and §10.3; J. C. Simo, T. A. Laursen, Comput. Struct. 42(1):97–116, 1992;
  rate-dependent coefficients reviewed by J. T. Oden, J. A. C. Martins, CMAME 52:527–634, 1985.
  The exponential decay μ = μ_k + (μ_s − μ_k) e^{−d γ̇} is the form the LS-DYNA and Abaqus
  manuals document, so the decks carry it unchanged.
- Kinematic rigid walls with Coulomb impulses: T. J. R. Hughes et al., CMAME 8(3):249–276, 1976
  (A5.5), the tangential impulse capped at μ times the normal one.
- Verification: the impact of equal bars, K. F. Graff, *Wave Motion in Elastic Solids*, Dover
  1991, §2.1–2.3, and W. Johnson, *Impact Strength of Materials*, Arnold 1972, §1.6; Hertz's
  impact of spheres, H. Hertz, J. reine angew. Math. 92:156–171, 1882, and K. L. Johnson,
  *Contact Mechanics*, Cambridge 1985, §11.4; the sliding block, any dynamics text (J. L. Meriam,
  L. G. Kraige, *Engineering Mechanics: Dynamics*, §3/5).

## Part R: Contact in nonlinear static steps (M7), as implemented

The papers `fieldcraft_fem::mortar`, `fieldcraft_fem::dual` and `fieldcraft_solver::contact` are
written from (cited again in the module docs); docs/fem.md describes the formulation.

- Mortar contact, segment to segment, with nodal weighted gaps and the pressure interpolated by
  the non-mortar side's shape functions: M. A. Puso, T. A. Laursen, *A mortar segment-to-segment
  contact method for large deformation solid mechanics*, CMAME 193:601–629, 2004; T. W. McDevitt,
  T. A. Laursen, *A mortar-finite element formulation for frictional contact problems*, IJNME
  48:1525–1547, 2000; the segmentation (projection on the non-mortar plane, clipping, integration
  cells) and the linearisation it needs, A. Popp, M. W. Gee, W. A. Wall, *A finite deformation
  mortar contact formulation using a primal–dual active set strategy*, IJNME 79:1354–1391, 2009.
  Chosen over two-pass node-to-surface because mortar passes the contact patch test on
  non-matching meshes (R. L. Taylor, P. Papadopoulos 1991; M. A. Crisfield 2000, as in Part Q).
- The consistent tangent by automatic differentiation of the mortar integrals rather than by
  hand: J. Lengiewicz, J. Korelc, S. Stupkiewicz, *Automation of finite element formulations for
  large deformation contact problems*, IJNME 85:1252–1279, 2011; forward-mode (truncated Taylor)
  arithmetic, A. Griewank, A. Walther, *Evaluating Derivatives*, 2nd ed., SIAM 2008.
- Augmented Lagrangian with Uzawa updates after each converged increment: J. C. Simo, T. A.
  Laursen, *An augmented Lagrangian treatment of contact problems involving friction*, Comput.
  Struct. 42(1):97–116, 1992; P. Wriggers, *Computational Contact Mechanics*, 2nd ed., Springer
  2006, ch. 6 and 10; T. A. Laursen, *Computational Contact and Impact Mechanics*, Springer 2003.
- Friction: the regularised Coulomb return map at the nodes, J. C. Simo, T. A. Laursen 1992 and
  P. Wriggers 2006, §5.2 (as in Part Q), on the mortar slip in its rate form (the weighted
  relative motion projected on the tangent plane), T. W. McDevitt, T. A. Laursen 2000 and M. A.
  Puso, T. A. Laursen, *A mortar segment-to-segment frictional contact method for large
  deformations*, CMAME 193:4891–4913, 2004; C. A. Coulomb's law, *Théorie des machines simples*,
  1785. The tangent is symmetrised like the follower pressure's (the non-symmetric coupling to
  the pressure is left out; the residual is exact).
- Ties by dual mortar constraints: B. I. Wohlmuth, *A mortar finite element method using dual
  spaces for the Lagrange multiplier*, SIAM J. Numer. Anal. 38(3):989–1012, 2000; A. Popp, B. I.
  Wohlmuth, M. W. Gee, W. A. Wall, *Dual quadratic mortar finite element methods for 3D finite
  deformation contact*, SIAM J. Sci. Comput. 34(4):B421–B446, 2012; mesh tying in solid
  mechanics, M. A. Puso, *A 3D mortar method for solid mechanics*, IJNME 59:315–336, 2004.
- The active set inside Newton (the max() of the regularised pressure, a semismooth Newton
  method): P. Alart, A. Curnier, *A mixed formulation for frictional contact problems prone to
  Newton like solution methods*, CMAME 92:353–375, 1991; Popp et al. 2009.
- Integration of the overlap cells: J. Radon's seven-point degree-5 rule on triangles (A. H.
  Stroud, *Approximate Calculation of Multiple Integrals*, Prentice-Hall 1971, T2:5-1).
- Verification: Hertz's contact of cylinders and spheres, H. Hertz, J. reine angew. Math.
  92:156–171, 1882; K. L. Johnson, *Contact Mechanics*, Cambridge 1985, §4.2–4.3 (eqs. 4.45–4.49
  for the line contact's pressure and subsurface stresses).
