# Eulerian solids on the lattice: the reference map technique

Research notes for issue #594: solids described on a fixed grid (the lattice the LB fluid
already uses) instead of a separate FEM mesh, and the lattice Boltzmann coupling of such solids
(LBRMT). This page is our own summary, written from the papers listed at the end; the papers
themselves are kept locally and never committed. The prototype that follows these notes is
`crates/eulerian` (solid only, 2D plane strain), with `crates/eulerian-gpu` for its deep
interior on a GPU; fluid coupling is a later step. **The prototype's code is not on `develop`:**
it is research code with no command, kept on the research branch of PR #597 until fluid–structure
coupling needs it, so the `crates/eulerian…` paths, tests and examples named below refer to that
branch. This page keeps the results.

## Summary

- **The published scheme is biased at loaded surfaces** (§8). Its blurred surface makes the
  Lamé cylinder 10.9 % too stiff with ten nodes across the wall, and the cause is the linear
  extrapolation of the reference map into the band, not the blur width: with the exact map
  prescribed there the error is +0.30 %. Its static answer also depends on the relaxation
  damping (−7.5 % to −20.7 % for γ × 0.5 to × 4), because band nodes integrate a velocity
  they never act on.
- **Mapped cut cells remove the bias** (§10) and are the default surface: momentum only in
  the body's cut control volumes, the traction imposed on the curved surface, the
  deformation gradient near it from a surface-fitted map. Lamé u_r(a) within 0.023 % at
  h = 1 (the hoop stress at the bore within 1.4 %), Kirsch RMS 0.53 % of 3S at h = 0.5, second
  order on stresses, independent of damping and viscosity.
  A ghost-value fill (§9) gets second order too but is less accurate and slower per step; it
  and the blurred surface stay only as options for the comparison studies.
- **Statics no longer need relaxation** (§11): Newton–Krylov with a geometric multigrid
  preconditioner solves the mapped equilibrium in O(n) with GMRES counts flat in h, about 430
  times faster than dynamic relaxation at h = 0.125 (1.55 s against 661 s for the Lamé cylinder).
  It is still about 200 times slower than FieldCraft's FEM for a solid alone.
- **Local refinement near the surface** (§12) works (nested levels, refluxing to round-off
  momentum, second-order subcycling) but pays modestly: at equal surface accuracy 1.9–6.7
  times fewer node updates and 0.8–2.4 times less CPU, and only for one or two levels,
  because the interior error stays at the coarse spacing.
- **Threads and a GPU** (§13): the step is deterministic on any thread count (bit for bit),
  but scales only 2.1–3.3 times on 8–10 cores because the surface work and the full-grid scans
  around it dominate. On the GPU the deep interior runs at 1.7–2.0 × 10⁹ node updates/s,
  70–80 times ten CPU threads, accurate in single precision once the stress is formed from
  small quantities; with the surface still on the CPU the hybrid step gains at most about
  12 %, and nothing on the smallest grid (1 % slower in the table's run).
- **What blocks FSI** (§15): the surface work on the CPU, and the time step of a stiff solid,
  hundreds of times below the lattice's; then the line-search stalls of the static solve after
  topology changes, and everything in 3D.

## 1. The idea

A Lagrangian solver follows material points. The reference map technique (RMT) does the
opposite: every grid point x stores the reference (undeformed) position ξ(x, t) of the material
that sits at x now. ξ is the inverse of the motion χ(X, t), and it is carried along with the
material, so it obeys the pure advection equation

    ∂ξ/∂t + (v · ∇) ξ = 0.

Its spatial gradient is the inverse deformation gradient, so

    F = (∇ξ)⁻¹,   J = det F,   ρ = ρ₀ / J,

and any finite-strain constitutive law σ = σ(F) (hyperelastic, or with internal variables
stored and advected beside ξ) gives the Cauchy stress at every grid point. The velocity follows
from the Eulerian momentum balance

    ρ (∂v/∂t + (v · ∇) v) = ∇ · σ + ρ b.

Everything lives on one fixed grid: no mesh, no remeshing, no Lagrangian markers. Displacement
boundary conditions are natural (the boundary's ξ is the boundary's reference position), and
a pre-strained state is just a non-identity initial ξ. Unlike Eulerian methods that advect F
itself, F stays compatible (it is a gradient), so no curl constraint has to be enforced.

The price: ξ is only defined where material is, the boundary of the solid is not on the grid,
and the stress needs a derivative of ξ, so every stencil that touches the boundary needs values
of ξ outside the solid.

## 2. Tracking the solid region

- Kamrin, Rycroft & Nave (2012) use a level set φ whose zero contour is the deformed boundary,
  with a sharp-interface treatment: values just outside the body are extrapolated linearly
  along grid lines from inside values and the boundary condition, and points that cross into
  the body during a step are given values by interpolation along the normal. Cut cells smaller
  than about a tenth of a cell are handled the same way.
- Valkov, Rycroft & Kamrin (2015) observe that the boundary is already encoded in ξ: if
  φ₀(X) is a signed distance to the body in the reference configuration, φ(x) = φ₀(ξ(x)) marks
  the deformed body without advecting a separate level set. The incompressible RMT (Rycroft,
  Wu, Yu & Kamrin 2020) and the LBRMT use this and reinitialise φ₀(ξ) to a signed distance in a
  narrow band each step (straddling points by a bicubic interpolant and Newton iterations, the
  rest by fast marching).
- From Valkov et al. (2015) on, the interface is blurred: a smoothed Heaviside H_ε(φ) over a
  transition zone of half-width ε (one to two cells) mixes solid and fluid stress and density,
  σ = H σ_f + (1 − H) σ_s, ρ = H ρ_f + (1 − H) ρ_s. Interface jump conditions are then satisfied
  implicitly in the zone rather than by explicit ghost-value fixes; only ξ is extrapolated
  (two scalar fields per step instead of the twenty of the sharp-interface method, by the
  authors' count).

## 3. Extrapolating ξ: the hard part

The solid stress is needed wherever H < 1, which includes the outer half of the blur zone
(0 < φ < ε), and its finite differences reach one more cell; ξ must therefore exist in a band
of a few cells outside the solid. ξ there is not advected (that region behaves like fluid and
would shear without bound, making the solid stress meaningless); it is rebuilt from the solid
each step.

- Sharp-interface RMT (2012): linear extrapolation along grid lines and PDE-based extension
  (Aslam 2004) of every needed field — many fields, and the papers report that accumulated
  extrapolation error could destabilise the interface (pressure spikes when points enter the
  solid).
- Valkov et al. (2015) introduce a constrained extrapolation that enforces physical and
  smoothness requirements on the extrapolated map.
- Rycroft et al. (2020) and the LBRMT replace this by a local least-squares fit: for each
  outside point, in order of increasing distance (2020: increasing φ; LBRMT: layer by layer
  outwards), fit a linear map ξ ≈ A x + c to all known values within a square window
  (half-width 3 in 2020, 2 in the LBRMT), weighted by distance, widening the window when the
  fit is degenerate (collinear points). The averaging blurs ξ far from the interface, which the
  authors argue is harmless because only the first cells matter, and it damps high-frequency
  noise.

Why it is hard: the extrapolated ξ decides (a) the stress in the blur zone, hence the traction
the surface carries, (b) where the interface is next step (φ₀(ξ) at the band points), so
material entering a cell gets an extrapolated, not an advected, reference position, and (c)
the deformation gradient computed across the interface. A linear fit is exact for homogeneous
deformation and first-order wrong where the strain varies, and it is the main source of the
first-order error near the boundary that the blurred methods show. Points where the fit is
degenerate (thin features, corners, two bodies touching) need special handling.

## 4. Numerical schemes in the papers

| | Kamrin, Rycroft & Nave 2012 | Valkov et al. 2015 | Rycroft et al. 2020 (IncRMT) | Sun & Rycroft 2025 (LBRMT) |
|---|---|---|---|---|
| fluid | compressible, explicit | compressible Navier–Stokes, explicit | incompressible, projection (Chorin), MAC + FE projections, multigrid | LB, D2Q9 BGK, Guo forcing |
| solid law | general hyperelastic, plane strain | compressible neo-Hookean | incompressible neo-Hookean (deviatoric) | incompressible neo-Hookean (deviatoric) |
| grid | v, ξ at cell corners; F, σ at cell centres | v, ξ corners, σ centres | cell-centred fields, stress on cell edges | LB nodes; F, σ on half-edges |
| advection | WENO (Liu, Osher & Chan 1994) | second-order upwinded ENO | second-order explicit Godunov scheme with Colella's limited slopes | second-order one-sided upwind differences |
| time | Euler step inside TVD RK3 (Shu & Osher 1988) | explicit | second order (projection) | LB step |
| interface | sharp, level set, extrapolation | blurred, φ = φ₀(ξ), one extrapolation | blurred, least-squares extrapolation, reinitialised φ | as IncRMT, blur zone + extrapolation zone of width ≥ 1.5Δx + 2ε |
| stabilisation | solid viscosity | artificial viscosity | artificial viscosity μ_e = κ ρ_s c_s h in the solid (κ = 0.4) | the LB fluid's own viscosity acts inside the solid |

Common points:
- **Stress stencil.** The deformation gradient is evaluated on cell edges (or centres) from
  compact differences of ξ: on a vertical edge ∂ξ/∂x from the two nodes it joins, ∂ξ/∂y as the
  average of the four transverse differences; the stress divergence at a node is then the
  difference of the edge stresses. This is compact (a 3 × 3 node stencil) and avoids the
  odd–even decoupling of a wide centred stencil.
- **CFL.** Explicit, limited by the solid wave speed: Δt ≲ α h / c with c = √(G/ρ) (shear) in
  the incompressible papers, the P-wave speed for compressible solids, plus the viscous limits
  Δt ≲ ρ h² / (2 μ (1 + q) d) of the fluid and the artificial viscosity.
- **Statics** are reached by relaxation: add damping (viscosity) and run the dynamics to rest
  (2012, §6.2: the washer, the flower). A good initial guess for ξ (a pre-strained map that
  already satisfies the boundary displacements) shortens the relaxation.
- **Verification** in the 2012 paper: a finite-strain travelling wave in a periodic domain
  (second-order convergence of ξ, v and stress), the Levinson–Burgess washer under finite
  rotation (closed form), uniform motion of a pre-deformed disk (exact to round-off because
  every field is linear). The incompressible paper's FSI tests converge at about order 1.8 in
  the L1 and L2 norms, with the largest errors at the interface (the blur zone alters the
  equations over a width O(h)) and first order in L∞ at the corners of a square solid.

## 5. How the LBRMT couples to the lattice

Sun & Rycroft (J. Comput. Phys. 2025; preprint 2024; Sun's 2020 master's thesis is an earlier
version of the method) put both phases on the D2Q9 lattice with one velocity field:

1. Advect ξ inside the solid (φ < 0) with the LB velocity; extrapolate ξ layer by layer into
   the extrapolation zone; relabel nodes with φ = φ₀(ξ).
2. Compute the solid's deviatoric stress on half-edges from ∇ξ, then the force density
   F = ∇ · [(1 − H_ε(φ)) τ_s] + (1 − H_ε) f_ext,s + H_ε f_ext,f.
3. Hand F to the LB collision as an external force with Guo's forcing (and the half-force
   velocity correction). The fluid stress is never computed explicitly: it is in the LB
   populations, and it also acts inside the solid, where it serves as the stabilising
   viscosity the IncRMT added by hand.
4. A "smooth flux correction" keeps the solid–fluid density difference: the equilibrium is
   shifted by the target density field ρ_t = H ρ_f + (1 − H) ρ_s, moving the excess flux back
   onto the rest population so mass and momentum stay conserved.

No interface mesh, no force spreading between Lagrangian and Eulerian points, no coupling
iteration; multi-body contact is a repulsive stress between overlapping blur zones (a separate
ξ field per body). Validation: a deformable disk in a lid-driven cavity against published
results, rotating and settling solids, suspensions mixing.

For FieldCraft this maps onto `crates/lbm` (D3Q19 TRT with Guo forcing already): the solid
supplies a force density per lattice cell, the lattice supplies the velocity; the prototype
keeps velocity and ξ at grid nodes for that reason.

## 6. Open problems (ours and the papers')

- **Accuracy at the boundary.** The blurred interface is formally first order at the
  interface (O(h) zone); the free surface sits
  somewhere in a zone ±ε wide, and ξ in that zone is extrapolated. Stresses at the surface (the
  numbers an analyst wants, e.g. the hoop stress at a bore) are the least accurate.
  §8 measures the bias; §10 removes it.
- **Extrapolation robustness**: thin parts a cell or two thick, sharp corners, near-contact,
  points that enter the solid with an extrapolated reference position.
- **Mixing / large distortion.** ξ keeps the link to the reference configuration, so extreme
  shearing (stirring) degrades ∇ξ; the 2012 paper suggests periodic reinitialisation of ξ with
  a stored F.
- **Level-set reinitialisation** cost and robustness (the IncRMT falls back to a first-order
  method when Newton fails).
- **Stability** relies on artificial viscosity or the LB viscosity; the solid's damping is
  then tied to grid or fluid parameters.
- **Stiff solids.** The explicit step is set by the solid's wave speed; for engineering metals
  in a slow flow, c_solid ≫ flow speed and the solid dominates the step count. The LBRMT's
  examples use soft solids (G comparable with fluid stresses).
- **Incompressible vs compressible.** The LB fluid is weakly compressible; the LBRMT's solid is
  incompressible neo-Hookean with the pressure from the lattice. A compressible solid on the
  lattice needs its volumetric stress in F.
- **3D** is direct in principle (3 × 3 ∇ξ, D3Q19) but the band and extrapolation cost grow
  with the surface.
- **Cost vs FEM** for a solid alone: an RMT grid resolves the whole bounding box and needs a
  dynamic relaxation to reach a static answer; it is not a substitute for the FEM solver for
  statics, only a route to FSI with large solid motion.
  §11 replaces the relaxation by a static solve; the FEM still wins.

## 7. The prototype and its verification

`crates/eulerian` (`fieldcraft-eulerian`, layer 2, no command yet): 2D plane strain on a
regular grid, f64, deterministic on any number of threads (§13), no panics, capped sizes and
step counts.
This section describes the first version, which follows the published scheme with its blurred
surface (`Boundary::Blurred`).
The default surface is now the mapped one of §10; the blurred one stays for the comparison
studies.

- Nodes carry ζ = ξ − x (periodic when the motion is) and v — the nodes a lattice will share.
- F on cell edges from compact differences; ∇·σ at nodes from the four edges.
- Blurred surface (ε = 1.5 h, the smoothed Heaviside of the 2020 paper); outside the body a
  prescribed *ambient stress* (0 for a free surface, −p I in a pressurised bore, any field), so
  surface tractions are σ_ambient·n with no surface tracking.
- φ = φ₀(ξ) (Valkov et al.), band of four 8-neighbour rings rebuilt each step by a weighted
  linear least-squares fit using only earlier rings (independent of node order, mirror
  symmetric). Band nodes with φ < 0 become solid (material entering a cell).
- Second-order one-sided upwind advection (first order next to unknown nodes); central
  differences in time (v at half steps), Δt = 0.4 h / c_P; artificial viscosity
  μ_e = κ ρ c_P h (κ = 0.1 by default); implicit mass-proportional damping for statics.
- Compressible neo-Hookean (Bonet & Wood) and small-strain linear elasticity.
- Later additions: the decoupled neo-Hookean of FieldCraft's FEM (μ/2 (Ī₁ − 3) +
  κ/2 (J − 1)²) for finite-strain comparisons, and a body force per unit volume.

Verification (`crates/eulerian/tests/verification.rs`):

| check | reference | resolution | result | tolerance |
|---|---|---|---|---|
| rigid translation, body crossing > 5 cells | ξ = x − vt, σ = 0 | h = 0.05 | ξ, v exact to 1e-11, ‖σ‖ < 1e-8 | 1e-11 / 1e-8 |
| rigid rotation by π/4 (prescribed velocity) | ξ = R(−θ)x, σ = 0 | h = 0.05 / 0.025 | max ‖ξ error‖ 4.7e-4 → 2.4e-4; max ‖σ‖/μ 3.4e-3 → 1.7e-3 (first order: the explicit update of ξ) | < 2e-3, ratio < 0.6 |
| uniaxial stretch s = 1.25, neo-Hookean | σxx = μ(s²−1)/s + λ ln s/s, σyy = λ ln s/s | h = 0.05 | exact to 1e-9 relative, block stays at rest | 1e-9 |
| simple shear γ = 0.3, neo-Hookean | σxx = μγ², σxy = μγ, σyy = 0 | h = 0.05 | exact to 1e-9, at rest | 1e-9 |
| free block under traction T, T/E = 1e-4 | εxx = T(1−ν²)/E, εyy = −Tν(1+ν)/E | h = 0.1 | linear: 9.0992e-5 / −3.9002e-5 (exact 9.1e-5 / −3.9e-5); neo-Hookean within 5.2e-5 relative | 1e-4 / 5e-4 |
| P-wave, periodic standing wave | √((λ+2μ)/ρ) = 11.6024 | 16 / 32 nodes per wavelength | 11.5395 / 11.5867 (−0.54 % / −0.14 %, second order) | 1 % / 0.25 %, ratio < 0.35 |
| S-wave | √(μ/ρ) = 6.2017 | 16 / 32 | 6.1638 / 6.1922 (−0.61 % / −0.15 %) | same |
| thick cylinder, plane strain (a = 10, b = 20, p = 100, E = 210000, ν = 0.3) | Lamé; FieldCraft FEM (quad8, oracle `lame-cylinder-plane-strain`) u_r(a) = 0.0090793 | h = 1 / 0.5 (10 / 20 nodes across the wall) | u_r(a) 0.008088 / 0.008789 (−10.9 % / −3.2 %); u_r(15) −11.2 % / −3.3 %; σθ(15) 81.9 / 89.3 vs 92.6 (−11.6 % / −3.5 %) | 15 % / 4 %, ratio < 0.4 |

u_r(a) is evaluated from a fit u_r = A r + B/r over the nodes with φ < −2h; the FEM gets the
same number to 3e-6 on a 8 × 16 quad8 quarter.
The checks except the Lamé ones now run with the blurred and the ghost-filled surfaces (§9).
The rigid motions, patch tests and wave speeds also run with the mapped surface
(`tests/mapped.rs`, §10).

Where it falls short:
- **The free surface.** Homogeneous states are exact (linear fields are reproduced by every
  stencil and by the extrapolation). With a curved surface and a strain gradient across it,
  the blur zone makes the ring too stiff: −11 % with ten nodes across the wall, −3.2 % with
  twenty (observed order ≈ 1.8, as the incompressible paper reports for its FSI tests). A wider
  blur is stiffer (ε = 1.0 h: −8.5 %; 2.5 h: −16.8 % at h = 1). Stresses at the surface itself
  are not usable; the numbers above are a node layer or more inside.
  §8 finds the cause.
- **Extrapolation** is linear: first-order wrong where the strain varies, and its fit needs
  three non-collinear nodes, so features thinner than about two cells are not represented.
- **Statics by relaxation.** The cylinder takes 1200 steps (2601 nodes) at h = 1 and about
  4000 steps (10201 nodes) at h = 0.5 to come to rest; the FEM solves the same problem with one
  sparse factorisation of under a thousand unknowns, more accurately. For a solid alone the RMT is
  not competitive with the FEM; its case is FSI with large solid motion on the fluid's grid.
- **Rotation** is transported with a first-order time error (explicit update of ξ with the
  current ∇ξ); second order in time would need a predictor–corrector or an RK stage.
- Not done: 3D, contact between bodies, prescribed moving boundaries, the incompressible law
  the LBRMT uses, the LB coupling (step 3).

## 8. Why the blurred surface is biased

Two suspects: the blur itself (the boundary condition spread over a band 2ε wide) and the
band's ζ, a linear extrapolation that continues the interior strain past the surface. The
`blur_diagnosis` example separates them on the Lamé cylinder (neo-Hookean, p = 100; errors of
u_r(a) from the fit over φ < −2h, relative), with two switches the solid keeps for this:
`Solid::set_exterior_map` prescribes the exact Lamé ζ in the band instead of extrapolating
it, and `Boundary::BlurredBody` integrates the momentum only at the body's nodes (φ < 0) with
v extrapolated into the band ("sm" in the study notes).

| band ζ | momentum at | ε / h | κ | γ × | h = 1 | h = 0.5 | h = 0.25 |
|---|---|---|---|---|---|---|---|
| extrapolated | φ < ε (published) | 1.5 | 0.1 | 1 | −10.9 % | −3.17 % | −0.85 % |
| exact | φ < ε | 1.5 | 0.1 | 1 | +0.30 % | +0.11 % | – (unstable) |
| extrapolated | φ < 0 | 1.5 | 0.1 | 1 | −39.1 % | −24.2 % | −13.8 % |
| exact | φ < 0 | 1.5 | 0.1 | 1 | +0.31 % | +0.088 % | +0.017 % |
| exact | φ < ε | 1.0 / 2.0 | 0.1 | 1 | +0.15 % / +0.53 % | +0.075 % / +0.18 % | |
| extrapolated | φ < ε | 1.0 / 2.0 | 0.1 | 1 | −8.5 % / −13.9 % | −2.3 % / −4.5 % | |
| extrapolated | φ < ε | 1.5 | 0 | 1 | −39.2 % | −24.2 % | |
| extrapolated | φ < ε | 1.5 | 0.1 | 0.5 / 2 / 4 | −7.5 % / −15.4 % / −20.7 % | | |
| exact | φ < 0 | 1.5 | 0.1 | 0.5 / 2 / 4 | +0.31 % / +0.31 % / +0.31 % | | |

What this shows:
- **The extrapolation is the bias.** With the exact ζ in the band the published scheme is
  within 0.3 % at h = 1; the blur width alone costs a few tenths of a percent (+0.15 / +0.30 /
  +0.53 % for ε = 1.0 / 1.5 / 2.0 h).
- **The published −10.9 % is a cancellation.** The band nodes with 0 < φ < ε integrate
  momentum, but their ζ is overwritten by the extrapolation every step, so they never come to
  rest, and the artificial viscosity turns their velocity into a surface stress that offsets
  part of the extrapolation error. Without the
  band's momentum (φ < 0), or without the viscosity, the true extrapolation error shows:
  −39 % at h = 1, order 0.4–1.0.
- **The static answer depends on the damping** used to reach it: −7.5 % to −20.7 % for
  γ × 0.5 to × 4 at h = 1, and on κ. With the momentum in the body and the exact band, both
  drop out.
- Measured and dropped (not committed): quadratic least-squares extrapolation (−3.7 % at h = 1,
  +0.16 % at h = 0.25, but noisy and divergent without viscosity), harmonic extensions (worse
  and unstable), other blur kernels (none removes the bias; the best with the exact band gives
  +0.10 % at h = 1), averaging over grid rotations (the angular spread also comes from the
  extrapolation).

A flat face isolates the same effect without curvature or staircase (a slab |x| < L under a
body force, σxx falling linearly to zero at the faces): the blurred surface puts the face at
an effective position off by a fixed fraction of h that depends mostly on where the face cuts
the grid, so the stress error is O(1) · (dσ/dn) · h, first order.

| | ε / h | face θ h past a node: 0 | 0.25 | 0.5 | 0.75 |
|---|---|---|---|---|---|
| blurred, effective face shift δ / h | 1.0 | +0.06 | +0.25 | +0.82 | +2.51 |
| | 1.5 | +0.13 | +0.30 | +0.66 | +1.47 |
| | 2.5 | +0.24 | +0.39 | +0.63 | +1.01 |
| blurred, σxx error at the last node / (ρ g h) | 1.5 | 0.82 | 0.91 | 1.21 | 1.98 |

δ / h is the same at h = 1, 0.5 and 0.25 (0.658 h at θ = 0.5, ε = 1.5 h): first order.

So a sharp surface needs two things: the traction imposed where the surface is, by the
solid's own stress (no blending), and no stress computed from extrapolated ζ.

## 9. A sharp boundary by ghost values (kept for comparison)

`Boundary::Ghost` (`crates/eulerian/src/sharp.rs`) was the first sharp surface. It is kept as
an option for the comparison studies; §10's mapped cut cells replace it.

Momentum is integrated at the body's nodes only; every edge with a body node at either end
uses the solid stress alone (no blending); the first 8-neighbour ring of nodes outside (the
ghosts) gets ζ from a local constrained least-squares fit: a polynomial ζ(d) per component
fitted to the body nodes in a 7 × 7 (cubic) or 5 × 5 (quadratic) window, subject to the two
equations [σ(G₀) + C:(∇ζ − G₀)]·n = σ_ambient·n at the surface point nearest the ghost
(the stress linearised about the previous fit, C by differences of the material law; two
passes, so the traction is exact at equilibrium for any law). The surface point is the
projection d ← d − φ∇φ/|∇φ|² onto φ₀(x + ζ(x)) = 0; the normal is ∇φ there. Outer band
layers stay linear (they only feed φ and the upwind advection).

Rules for the hard cases:
- small cut cells: nothing divides by the cut distance, so a surface arbitrarily close to a
  node gives a regular system;
- a surface that comes to rest on a node: body nodes leave at φ ≥ 0, band nodes enter only at
  φ < −h/10 (without this hysteresis the finite-strain cylinder rang forever, 8 bore nodes
  flipping in and out); a ghost may therefore lie up to h/10 inside, and the projection runs
  from either side;
- fits the data barely determine (the inside of a convex corner): the KKT system is
  equilibrated and a pivot below 1e-6 refuses the fit; then the quadratic, then a linear fit
  with the same traction, then the plain extrapolation (counted; zero in every case below);
- other things measured and dropped: other fit weights (no change: the error is truncation),
  extra traction rows at the window's edge crossings (doubled the error: they over-constrain
  the fit's gradient), the quadratic fit (exact for linear stress, but 3–5 times the cubic's
  error at a curved surface).

Stability: with no artificial viscosity (κ = 0) a growing surface mode appears after a few
thousand steps in the cylinder, at any time step (CFL 0.4 and 0.2 alike): the closure is not
symmetric. κ = 0.02–0.1 damps it, with or without the relaxation damping, and the static
answer does not depend on κ or γ.

Verification (linear elastic at strain 1e-6 unless noted; `examples/boundary_study.rs`). The
slab of §8 is exact to the relaxation tolerance wherever its faces cut the grid (effective
face shift −7e-7 to −6e-6 h). Lamé cylinder (a = 10, b = 20; relative errors; σ at the bore
extrapolated from the three nearest nodes on y = 0):

| h | blurred u_r(a) | order | ghost quadratic u_r(a) | order | ghost cubic u_r(a) | order | cubic σθ(15) | cubic σrr(a)/p (−1) | cubic σθ(a) |
|---|---|---|---|---|---|---|---|---|---|
| 2 | −2.61e-1 | – | −1.07e-1 | – | −1.42e-1 | – | −1.33e-1 | −0.838 | −1.80e-1 |
| 1 | −1.12e-1 | 1.22 | −3.69e-2 | 1.54 | −2.57e-2 | 2.47 | −2.71e-2 | −0.986 | −4.61e-2 |
| 0.5 | −3.45e-2 | 1.70 | −1.11e-2 | 1.73 | −4.01e-3 | 2.68 | −5.32e-3 | −1.0015 | −1.01e-2 |
| 0.25 | −1.01e-2 | 1.78 | −3.12e-3 | 1.83 | −6.29e-4 | 2.67 | −1.25e-3 | −1.0011 | −2.15e-3 |
| 0.125 | −2.75e-3 | 1.87 | −7.64e-4 | 2.03 | −5.89e-5 | 3.42 | +6.6e-5 | −1.0004 | −4.14e-4 |

The maximum nodal stress error of the cubic fit converges at order 2.2 (1.5e-2 p at h = 0.5).
The blurred surface's σrr at the bore is −0.03 / −0.46 / −0.72 / −0.86 / −0.93 p for h = 2 …
0.125; the ghost fill reaches −p within 0.2 % from h = 0.5.

Kirsch hole (a = 5 in a disc R = 20 with the Kirsch traction on r = R, so the finite domain
is exact; σθ at the hole edge extrapolated along the axis):

| h | blurred σθ(a, 90°)/S (3) | σθ(a, 0°)/S (−1) | ghost quadratic | | ghost cubic | |
|---|---|---|---|---|---|---|
| 1 | 2.184 | −0.254 | 3.169 | −1.338 | 3.166 | −1.321 |
| 0.5 | 2.693 | −0.668 | 3.056 | −1.111 | 3.010 | −1.054 |
| 0.25 | 2.938 | −0.897 | 3.020 | −1.036 | 2.997 | −1.007 |
| 0.125 | 3.002 | −0.972 | 3.006 | −1.011 | 2.998 | −1.000 |

The blurred value at 90° crosses 3 S on its way from below (+0.05 % at h = 0.125 is a
coincidence: its error around the hole oscillates between −0.02 S and +0.07 S at h = 0.25,
against ±0.02 S for the cubic fit). Maximum nodal stress error: blurred 0.135 S, cubic 5.4e-3 S
at h = 0.125.

Rotating annulus (a = 10, b = 20, free surfaces, centrifugal body force, closed form
u = C₁ r + C₂/r − K r³; non-harmonic):

| h | blurred u_r(a) | σθ(a) | ghost cubic u_r(a) | σθ(a) | max σ error / σθ(a) |
|---|---|---|---|---|---|
| 1 | −7.8e-2 | −2.3e-2 | −1.96e-2 | −2.8e-2 | 2.5e-2 |
| 0.5 | −2.4e-2 | +1.8e-2 | −3.0e-3 | −6.2e-3 | 5.6e-3 |
| 0.25 | −6.8e-3 | +1.9e-2 | −4.8e-4 | −1.3e-3 | 1.2e-3 |

(orders of the ghost fill 2.7 / 2.2 / 2.2; the blurred σθ(a) does not converge). The static
answer no longer depends on the damping: u_r(a) at h = 1 equal to 2e-5 for γ × 0.5 and × 4
(blurred: −7.9 % against −20.9 %).

**Finite strain: neo-Hookean cylinder.** The decoupled neo-Hookean law of FieldCraft's FEM
(μ, κ from E = 210000, ν = 0.3), bore pressure p = 8000 (p/E = 0.038; the linear solution
scaled to it is u_r(a)/a = 7.3 %, the nonlinear one 8.2 %), applied in ten increments.
Reference: FieldCraft's FEM, nonlinear static with follower pressure, plane-strain quad8
quarter, converged by refinement (16 × 32 → 32 × 64: u_r(a) 0.815032 → 0.815032, σθ at r = 15
in the deformed body 8777.3 → 8777.5).

| | u_r(a) | σθ(r = 15) | σθ at the bore | σrr at the bore |
|---|---|---|---|---|
| linear Lamé scaled to p | 0.726 | 7407 | 13333 | −8000 |
| FEM 32 × 64 | 0.815032 | 8777.5 | 14700 | −7998 |
| ghost, h = 2 | 0.6440 | 7114 | 11548 | −7188 |
| ghost, h = 1 | 0.7924 | 8471 | 13972 | −8151 |
| ghost, h = 0.5 | 0.8108 | 8714 | 14558 | −8040 |
| ghost, h = 0.25 | 0.8144 | 8765 | 14677 | −8009 |
| ghost, h = 0.125 | 0.81499 (−4.7e-5) | 8777.9 | 14697 | −8004 |
| Richardson (0.5, 0.25, 0.125) | 0.81510 (order 2.7) | 8782.5 (order 1.95) | 14701 (order 2.6) | |

The blurred surface folds (det ∇ξ < 0) during the load ramp at every resolution at
p = 8000 and at p = 4000, and at some resolutions at p = 2000–3000.

**Cost** (Apple M4 Max, release build, one thread per run, cylinder, linear):

| h | nodes | blurred steps | ms/step | wall s | ghost steps | ms/step | wall s | ghost u_r(a) error |
|---|---|---|---|---|---|---|---|---|
| 2 | 1089 | 780 | 0.15 | 0.1 | 1200 | 1.1 | 1.3 | 1.4e-1 |
| 1 | 2809 | 1860 | 0.35 | 0.6 | 2140 | 2.3 | 4.9 | 2.6e-2 |
| 0.5 | 8649 | 4240 | 0.93 | 3.9 | 4300 | 4.8 | 20.5 | 4.0e-3 |
| 0.25 | 29929 | 7820 | 2.9 | 22.4 | 9120 | 10.6 | 96.5 | 6.3e-4 |
| 0.125 | 110889 | 16480 | 9.8 | 161 | 17220 | 34.4 | 591 | 5.9e-5 |

A profile at h = 0.25: ghost fits 65 %, grid neighbour lookups 12 %, the linear band fill 11 %,
edge stresses and momentum 6 %. Steps to rest grow as 1/h, so a static answer by relaxation
costs O(h⁻³) in 2D: 10³–10⁵ times the FEM, which reaches u_r(a) within 6e-3 on 1 × 2 quad8
elements, 6e-4 on 2 × 4 and 4e-5 on 4 × 8 (121 nodes; σθ(15) within 8e-4), each in about
1 ms, and the finite-strain case in 0.04–0.9 s (8 × 16 to 32 × 64).

## 10. Mapped cut cells (the default surface)

`Boundary::Mapped` (`crates/eulerian/src/{cut,mapped}.rs`) is the default surface. It removes
both causes of §8: no blending, no stress from extrapolated ξ, momentum only inside the body.

- **Control volumes.** A node inside the body (φ < 0) carries momentum for the solid part of
  its dual cell. The surface in each cut grid cell is a quadratic curve through its two edge
  crossings and a third point on φ₀(ξ) = 0 above the chord. Areas, centroids, the solid
  pieces of the faces and the surface pieces of each quarter cell come from Green's theorem
  over that curve (three Gauss points per side, exact for it), so every control volume closes
  exactly and a constant stress puts no net force on any cut cell. The lumped mass is ρ times
  the true solid area. Solid in the dual cell of an outside node is shared between the inside
  neighbours across the surface along x and along y, with smooth weights in |∂φ/∂x| against
  |∂φ/∂y|; its net force and mass go with it. A control volume below a quarter of h² would
  join its deepest inside neighbour (never needed in the cases below).
- **Fluxes.** Faces between full cells keep the compact interior stencil. Faces next to a cut
  cell integrate σ(F)·n with two Gauss points on each solid piece; on the curved side the flux
  is the prescribed traction σ_ambient·n with the curve's line Jacobian.
- **The deformation gradient near the surface** comes from a surface-fitted map. Where a grid
  line along axis d crosses the surface (up to about 72° from the normal) there is a surface
  point with its own value ξ_S. Between neighbouring lines, T(u, η) = (x_d′ = (l + u) h,
  x_d = Y(u) − σ η h) maps the unit square onto a cell whose side η = 0 is the surface, Y(u)
  the crossings interpolated with the same Lagrange shape functions in u as the field;
  ζ(u, η) = Σ L_l(u) p_l(η), p_l quadratic through ξ_S and the two inside nodes on line l;
  ∇ζ = (∂ζ/∂u, ∂ζ/∂η) T′⁻¹ and F = (I + ∇ζ)⁻¹. Near the diagonal both families of lines give
  a map and the two are blended smoothly.
- **Surface values** satisfy σ(F)·n = t at every surface point: Newton on the conditions of
  each chain of points along the surface, banded finite-difference Jacobian, LU with partial
  pivoting, warm-started from the previous step (two or three iterations).
- **The band** is still filled by the linear least-squares extrapolation, but only the ξ
  advection (upwind neighbours) and the crossing positions use it.

Three things the first version got wrong and the final one fixes: a discrete choice that can
flip between steps makes the force discontinuous and the explicit scheme chatters instead of
coming to rest (the family of lines chosen at 45°, the inside node that takes an outside
node's solid); Gauss–Seidel along the surface diverges at the ends of a family of lines on a
tight curve (hence the banded Newton); the cut state must belong to the configuration it is
used with (it is rebuilt at the end of each step).

Verification (`crates/eulerian/tests/mapped.rs`, and the `rmt_cases` example for the
studies; same cases and metrics as before, u_r(a) from the fit over φ < −2h):

| check | result |
|---|---|
| uniaxial stretch 1.25, simple shear 0.3 (neo-Hookean disc held by the matching ambient stress) | exact to 1e-9 at every node and surface point, at rest |
| rigid translation of a disc across the grid | exact to round-off; the surface map amplifies it where the surface passes very close to a node (1/s weights): v to 5e-10, σ to 4e-9 (blurred: 1e-11) |
| rigid rotation of a free disc by 90° (dynamic) | ‖σ‖/μ 3.0e-3 (the centrifugal response ρω²R²/μ = 6.5e-3), ξ within 2e-3 |
| P- and S-wave speeds (no surface) | the plain stencil: same as §7 |

Lamé cylinder, |u_r(a) error| (%):

| h | blurred | blur + exact ξ in band (§8) | mapped, p = 100 neo-Hookean | mapped, small strain (load × 0.01) |
|---|---|---|---|---|
| 2 | 30.3 | – | 0.85 | 0.99 |
| 1 | 10.9 | 0.31 | 0.154 | 0.023 |
| 0.5 | 3.17 | 0.088 | 0.181 | 0.0004 |
| 0.25 | 0.85 | 0.017 | 0.184 | 0.0018 |
| 0.125 | 0.11 | – | 0.185 | 0.0016 |

The mapped error at p = 100 levels at +0.18 %: that is the finite-strain answer (the pressure
acts on the expanded bore, u is measured at the current position), not a discretisation
error; with the load cut a hundredfold it falls to 0.0018 %, the same floor scaled. Below h =
0.5 the small-strain error sits on that floor, so the order is read off the stresses: the
worst σθ over the bore surface points (small strain) is 6.7, 1.4, 0.32, 0.082, 0.026 % for h =
2 … 0.125 (orders 2.3, 2.2, 2.0, 1.7). Rotating disc (ρω² = 0.5, both surfaces free): u_r(a)
small strain −6.9, −0.015, +0.003, +0.002, +0.001 %; worst bore σθ 4.2, 0.89, 0.20, 0.052,
0.017 % (blurred u_r(a): −47, −7.6, −2.2, −0.54, −0.04 %). Kirsch hole (a = 5 in a disc of
radius 20 with the exact traction outside): RMS error of σθ around the hole edge 3.4, 3.7,
0.53, 0.12, 0.027 % of 3S (orders 2.8, 2.2, 2.1 from h = 1); σθ(a, 90°) within 0.04 % of 3S
from h = 0.5; the blurred surface's worst node next to the hole is 30, 18, 6.2, 3.6, 2.1 % off.

- **Damping** no longer moves the answer: u_r(a) at h = 0.5 is 0.00022 / 0.00036 / 0.00041 %
  for γ × 0.5 / 1 / 4 (blurred: −1.9 / −3.2 / −7.7 %), nor does the artificial viscosity
  (κ = 0.02: +0.150 % against +0.154 % at h = 1; blurred −22.6 % against −10.9 %).
- **The band** decides nothing measurable: 2, 4 or 6 extrapolation layers give identical
  results; prescribing the exact ξ there moves u_r(a) by 0.004 % (h = 1) and 0.015 % (0.5).
- **Time step.** With κ = 0.02 (so the viscous limit does not cap Δt) the mapped scheme relaxes
  up to Δt = 0.7 h/c_P at h = 1 and 0.8 at 0.5, the blurred one up to 0.8; the default 0.4 has
  margin. The largest eigenvalue of the mapped operator is 0.22–0.26 of the explicit limit at
  0.4. The operator is not symmetric near the surface (asymmetry 0.17 of its largest entry at
  h = 1); its complex eigenvalues grow by up to 0.5 % per step without damping, so undamped
  dynamics need the artificial viscosity (κ = 0.1 suffices; the blurred scheme without it
  diverges at t = 0.5).
- **Cost** per step on one core (Lamé; blurred / mapped, ms): h = 2: 0.10 / 0.66; 1: 0.31 /
  1.5; 0.5: 1.0 / 3.5; 0.25: 3.6 / 9.0; 0.125: 13.8 / 26.6. The cut cells (≈ 0.9 ms at h =
  0.5) and the surface solve (≈ 0.6 ms) grow with the surface, about × 2.2 and × 1.5 per
  halving of h; the rest is the interior's O(n) work, so the ratio falls from 6.6 to 1.9.

**The three surfaces side by side** (small strain; the blurred and ghost columns from §9's
study at strain 1e-6, the mapped one from this section's at load × 0.01; per-step costs from
the two studies' own runs on one core):

| | h | blurred | ghost fill (cubic) | mapped |
|---|---|---|---|---|
| Lamé u_r(a) error | 1 | −11.2 % | −2.57 % | 0.023 % |
| | 0.5 | −3.45 % | −0.40 % | 0.0004 % |
| | 0.25 | −1.01 % | −0.063 % | 0.0018 % (floor) |
| Kirsch σθ(a, 90°)/S (exact 3) | 0.5 | 2.693 | 3.010 | within 0.04 % of 3 |
| | 0.25 | 2.938 | 2.997 | within 0.04 % of 3 |
| rotating disc u_r(a) error | 1 | −7.8 % | −1.96 % | −0.015 % |
| | 0.5 | −2.4 % | −0.30 % | +0.003 % |
| static answer, γ × 0.5 against × 4 | 1 / 0.5 | −7.9 / −20.9 % (h = 1) | equal to 2e-5 (h = 1) | 0.00022 / 0.00041 % (h = 0.5) |
| order of the surface stress | | none (σθ(a) does not converge) | 2.2 | 2.0–2.3 |
| ms per step, h = 0.125 | | 9.8–13.8 | 34.4 | 26.6 |

## 11. Statics without relaxation: Newton–Krylov with multigrid

Dynamic relaxation needs a number of steps that grows as 1/h (the explicit step is bounded by
h/c_P while the slowest mode's period is fixed), so its cost is O(n h⁻¹) = O(h⁻³) on n grid
nodes: 661 s of CPU for the Lamé cylinder at h = 0.125. `Solid::solve_static`
(`crates/eulerian/src/{statics,multigrid}.rs`) solves the equilibrium of the mapped boundary
directly. In statics the reference map is the inverse motion, ξ(x) = X(x); nothing is advected.

**The problem.** Unknowns: ζ = ξ − x at the body nodes that hold a control volume (a merged
node moves with the node it joined). Residual R(ζ): the net force on every control volume of
§10 (face fluxes of σ(F), the traction on the curved surface, the body force) at a load factor
λ, with the surface values ζ_S solved from their traction conditions and the band refilled from
the body inside every evaluation. R = 0 is exactly the rest state of the dynamic scheme.

**The body depends on the unknown.** {φ₀(ξ) < 0} moves while Newton converges, so the unknown
set and the cut structure change. Each linearisation and its line search *hold the topology*:
the node kinds, which surface points carry values and how they are matched, the merged control
volumes. With it held, crossings, areas, normals and surface values move continuously with ζ.
A trial step that moves the surface across a node is judged on the rebuilt topology (a node
that entered the body starts from its extrapolated ζ); the topology is also re-derived at
every converged state, and convergence counts only when the two agree. The load is ramped by a
load factor with a secant predictor; an increment whose Newton fails (line search below 1/64)
is halved. The structure that depends only on the node kinds (cut-cell list, crossing edges,
full-cell areas, face list, the band fill as a linear map) is built once per topology; a held
evaluation recomputes only the near-surface geometry.

**Rigid motions.** The four bodies are free in the plane. The equilibrium the relaxation reaches
is R ∈ span(MΦ) (Φ: two translations and the rotation about the centre of mass; M: lumped
masses), because the discrete net load drives a rigid drift. The solver measures QR with
Q = I − MΦ(ΦᵀMΦ)⁻¹Φᵀ (inertia relief) and keeps updates free of rigid motion with
P = I − Φ(ΦᵀMΦ)⁻¹ΦᵀM, so GMRES solves the consistent singular system QJP δ = −QR.

**Jacobian: matrix-free (JFNK).** J v ≈ (R(ζ + εv) − R(ζ))/ε with ε‖v‖∞ = 10⁻⁷ h. The residual
contains an implicit solve (the surface values), root finds (crossings, arc midpoints) and
least-squares fits (the band); an analytic Jacobian would have to differentiate all of them,
the surface values through the implicit function theorem, and follow every change of the
discretisation. The difference quotient is the exact Jacobian's action at one residual
evaluation per Krylov iteration, so Newton keeps its quadratic convergence; only the
preconditioner needs a matrix, and an approximate one does. Krylov: restarted GMRES(60), right
preconditioned (the operator is not symmetric near the surface), Eisenstat–Walker forcing
(choice 2, η ≤ 10⁻²). Converged when ‖QR‖∞ ≤ 10⁻⁹ times the larger of the largest external
force on a control volume and the largest force through a face.

**Preconditioner.** A geometric V-cycle on Galerkin coarse operators A_{l+1} = Pᵀ A_l P over the
grid hierarchy (every second node), point-block Gauss–Seidel (2 sweeps forward before, 2
backward after), a dense LU on the coarsest level (≤ 400 nodes) bordered with the rigid
motions. *Coarsening rule at the cut cells:* a coarse node is a coarse unknown when the fine
node on it is a fine unknown, plus, for a fine unknown left without one among its bilinear
parents, its parent of largest weight; each fine unknown interpolates from the parents it has,
weights renormalised. Constants are reproduced everywhere, linear fields exactly in the
interior and to O(h) at the surface, and P has full column rank. (Taking every parent of a fine
unknown keeps linear fields exact but makes P rank deficient where several coarse nodes outside
the body parent only the same cut-cell node: the first version's Galerkin operator had a fourth
null vector and its coarse LU failed.) Geometric multigrid works at the cut cells; neither AMG
nor ILU was needed. The fine matrix A₀ is a *ghost-cell* approximation of the Jacobian,
assembled in O(n) without a residual evaluation: every face with solid on it carries the
tangent of its compact-stencil flux (∂f/∂(∇ζ) by differences of the constitutive law) times its
solid fraction, ζ at band nodes expanded through the band fill, and each dual cell's rows go to
the control volumes holding its solid with their shares. Away from the surface that is the
exact Jacobian; at the surface it ignores the cut geometry, the surface map and the surface
values. The alternative implemented beside it probes the rows near the surface (2 s²
residual evaluations perturbing one colour of unknowns at a time, s = 5 or 7, interior rows
from the tangent). It is closer at the surface but costs 50–100 evaluations per build, and the
ghost cells give the same iteration counts (table below), so they are the default.

**Reuse.** The hierarchy is kept while the topology is unchanged, across Newton steps and load
increments, and rebuilt when it is stale (GMRES needs more than twice the iterations it took
when fresh). The factorised banded Jacobian of every surface chain is kept the same way, and
the surface solve inside each residual evaluation becomes a chord iteration (refactorised when
it contracts slowly). A kept hierarchy *adapted* to a changed topology (unknowns matched by
node, new ones on their probed diagonal blocks) cost 45–225 GMRES iterations instead of about
20 even for 12–19 new nodes, growing with 1/h, so a topology change always rebuilds; with the
ghost-cell matrix a rebuild is cheap.

Verification (static against relaxed; relaxed with the §10 settings and tolerance, then
compared node by node; "non-rigid" removes the best-fit rigid motion of the difference, which
the free relaxed body accumulates):

| case | h | max \|ξ_static − ξ_relaxed\| | non-rigid | rigid drift | max \|u\| | body nodes differing |
|---|---|---|---|---|---|---|
| Lamé | 2 / 1 / 0.5 / 0.25 | 2.7e-7 / 7.7e-7 / 1.8e-6 / 5.8e-7 | 2.6e-7 / 7.7e-7 / 1.8e-6 / 5.3e-7 | ≤ 2e-8 | 9.1e-3 | 0 |
| rotating disc | 2 / 1 / 0.5 / 0.25 | 4.3e-7 / 7.6e-7 / 1.8e-6 / 1.0e-6 | 2.6e-7 / 7.6e-7 / 1.8e-6 / 1.0e-6 | ≤ 2e-7 | 7.7e-3 | 0 |
| Kirsch | 2 / 1 / 0.5 / 0.25 | 3.4e-7 / 7.8e-6 / 8.4e-6 / 9.8e-8 | 1.4e-7 / 1.6e-7 / 7.8e-8 / 8.7e-8 | up to 6e-6 | 1.0e-2 | 0 |
| neo-Hookean, p = 8000 | 2 | 2.0e-2 | 9.9e-4 | 1.1e-2, turn 4.4e-4 | 0.79 | 0 |

The differences fall with the relaxation's tolerance (Lamé, h = 1: 7.7e-7, 1.1e-7, 3.9e-9 for
the tolerance × 1, 0.1, 0.01), so the static solution is the relaxation's limit; the tests
(`crates/eulerian/tests/statics.rs`) check Lamé at h = 1 against a relaxation ten times
tighter than §10's (max |Δu| below 10⁻⁴ of max |u|) and a finite homogeneous stretch to 10⁻⁹.

Neo-Hookean cylinder (decoupled law, p = 8000, FEM u_r(a) = 0.815032), static solve:

| h | 1 | 0.5 | 0.25 | 0.125 | 0.0625 |
|---|---|---|---|---|---|
| u_r(a) | 0.811736 | 0.814304 | 0.814869 | 0.814995 | 0.815024 |
| error | −4.0e-3 | −8.9e-4 | −2.0e-4 | −4.5e-5 | −9.8e-6 |

(observed order 2.1–2.2). The static/relaxed comparison of this case exists only at h = 2,
and it is weak there: the relaxed body turned rigidly by 4.4e-4 rad and moved 1.1e-2, and a
body turned on the grid has a discrete equilibrium that differs by about the discretisation
error (the non-rigid difference, 1.3e-3 of u, is of that size at h = 2). The ramped relaxation
(ten increments, at most 300 000 steps each) took 303 240 steps at h = 2 and had not finished
at h = 1 after more than 20 min of CPU, so the finer grids are unverified against it; the
static answers converge to the FEM's at second order instead.

**Iterations** (one core, CPU s; N Newton steps; K GMRES iterations, mean and largest per step):

| h | nodes | Lamé N / K | disc N / K | Kirsch N / K | NH N / K (increments) |
|---|---|---|---|---|---|
| 2 | 676 | 3 / 16, 25 | 3 / 17, 25 | 4 / 11, 15 | 10 / 14, 26 (3) |
| 1 | 2601 | 4 / 13, 23 | 4 / 13, 21 | 3 / 15, 21 | 11 / 9, 20 (3) |
| 0.5 | 10201 | 3 / 15, 22 | 3 / 15, 20 | 3 / 14, 21 | 12 / 10, 20 (3) |
| 0.25 | 40401 | 3 / 15, 23 | 3 / 12, 20 | 3 / 14, 22 | 17 / 10, 24 (4) |
| 0.125 | 160801 | 3 / 13, 18 | 3 / 14, 20 | 3 / 13, 19 | 20 / 10, 23 (4 + 1 failed) |
| 0.0625 | 641601 | 4 / 12, 22 | 3 / 15, 22 | 3 / 11, 18 | 27 / 9, 19 (5 + 2 failed) |

GMRES iterations to 10⁻⁸ on the first Newton system (Lamé; the other cases with the ghost-cell
matrix within 20–28 at every h):

| h | 2 | 1 | 0.5 | 0.25 | 0.125 | 0.0625 |
|---|---|---|---|---|---|---|
| unknown nodes | 236 | 928 | 3756 | 15044 | 60288 | 241264 |
| multigrid, ghost-cell matrix (default) | 24 | 20 | 26 | 23 | 24 | 25 |
| multigrid, probed matrix, s = 5 | 15 | 21 | 19 | 26 | 27 | 31 |
| multigrid, probed matrix, s = 7 | 14 | 18 | 17 | 17 | 18 | 21 |
| probed matrix as the operator (the V-cycle alone) | 11 | 20 | 17 | 22 | 22 | 25 |
| no preconditioner | 181 | 360 | 756 | 3107 | – | – |

GMRES counts per Newton step are h-independent, and so are the Newton counts of the
linear-load cases; the neo-Hookean cylinder's total grows (10 → 27) because finer grids see
more topology changes and failed increments during the ramp (about 5 Newton steps per accepted
increment throughout). With the probed matrix the counts creep up with s = 5 (its lumping of
the surface coupling) and much less with s = 7; the V-cycle alone on its own matrix creeps up
slowly too (11 → 25), mostly from the coarse levels at the surface (a W-cycle gave 14, 14, 19
at h = 1, 0.25, 0.0625; stopping the coarsening at about 240 instead of 60 nodes took h =
0.0625 from 28 to 25).

**Cost** (CPU s, one core, Apple M4 Max with other jobs running; the relaxation is the §10
driver to its default tolerance):

| h | Lamé static | relaxation | ratio | disc static | relaxation | Kirsch static | relaxation | NH static |
|---|---|---|---|---|---|---|---|---|
| 2 | 0.04 | 0.85 | 21 | 0.05 | 0.92 | 0.04 | 0.78 | 0.13 |
| 1 | 0.09 | 3.6 | 40 | 0.10 | 4.4 | 0.09 | 3.9 | 0.24 |
| 0.5 | 0.18 | 16.3 | 91 | 0.20 | 17.7 | 0.20 | 21.9 | 0.61 |
| 0.25 | 0.56 | 102 | 182 | 0.46 | 105 | 0.54 | – | 2.5 |
| 0.125 | 1.55 | 661 | 426 | 1.60 | – | 1.62 | – | 9.4 |
| 0.0625 | 6.2 | – | | 5.6 | – | 5.0 | – | 37 |

The static cost is O(n) = O(h⁻²): × 2.8 and × 4.0 per halving for Lamé, × 3.8 and × 4.0 for
the neo-Hookean cylinder from h = 0.25; the relaxation is × 6.3 to × 6.5 per halving, O(h^−2.7)
here (steps ∝ h⁻¹ times O(n) per step, the surface work making up the rest). The FEM still wins
by orders of magnitude: the quad8 quarter reaches u_r(a) within 4e-5 on 4 × 8 elements in about
1 ms; the static RMT matches that at h = 0.5 (§10's small-strain error) in 0.18 s, about 200 ×
the FEM, and the bore stress to 0.3 % only from h = 0.25. The neo-Hookean cylinder takes the
FEM 0.04–0.9 s (8 × 16 to 32 × 64); the static RMT 9.4 s at h = 0.125 for −4.5e-5.

**Reuse, measured.** Keeping the surface chains' factorisations cuts the chain factorisations
15–20 × and the traction evaluations 4–5 ×, about halving the cut-and-surface stage of an
evaluation; with the ghost-cell matrix, keeping the hierarchy saves one build of two or three.
Together they change the total by up to 20 % (Lamé h = 0.125: 1.40 s rebuilding everything,
1.28 keeping the hierarchy, 1.46 keeping the LUs, 1.50 keeping both; the neo-Hookean cylinder:
10.6, 10.2, 11.9, 9.1 s); the spread between runs of a different Newton path is as large. With
the probed matrix (50 evaluations per build) keeping the hierarchy matters: Lamé 0.60 / 1.60 /
5.38 s rebuilding against 0.46 / 1.32 / 4.01 s keeping it at h = 0.5 / 0.25 / 0.125 (−25 %),
the neo-Hookean cylinder 6.1 / 10.4 / 63.8 s against 3.7 / 7.9 s at h = 0.5 / 0.25 (−25 to
−40 %). The ghost-cell matrix made that reuse mostly unnecessary.

**Where the time goes** (share of the solve at h = 0.125 / 0.0625): face forces (the interior's
compact fluxes and the surface map at the cut faces) 26–37 %; V-cycles 15–24 %; the ghost-cell
assembly and multigrid setup 15–24 %; the cut cells and surface solve 6–16 %, falling with h;
topology rebuilds 4–8 % (16 % for the neo-Hookean cylinder, whose bore moves 6.5 cells at
h = 0.125 and changes the topology 12–16 times); the rest is GMRES orthogonalisation and vector
work. Every residual evaluation costs O(n) (interior fluxes) plus O(h⁻¹) with a large constant
(the surface map, the cut cells, the surface solve), and a solve takes 40–60 of them for the
linear cases and 120–310 for the neo-Hookean cylinder.

Where it falls short:
- **Not faster than the FEM.** The cost order is now that of the FEM's iterative solvers, but
  the constant is not: a residual evaluation does the cut geometry, the surface solve and the
  surface map, and a solve needs 40–300 of them. For a solid alone the FEM remains the tool; the
  static solve matters for the RMT inside a fluid grid (an initial equilibrium, quasi-static
  loading, an implicit solid step) and as a fast verification path.
- **Finite strain.** The surface sweeps over nodes during the load ramp, so the topology changes
  every few Newton steps (16 times at h = 0.0625) and each change rebuilds the preconditioner.
  Right after such a change the line search sometimes stalls (step below 1/64 with ‖QR‖ flat),
  most likely at surface points a node only just behind the surface (depth s → 0, map weights
  ∝ 1/s); the increment is then halved (1–2 failures per solve at h ≤ 0.125). A surface stencil
  that skips nodes too close to the surface would remove it, but changes the discretisation.
- **The neo-Hookean cylinder is verified against the FEM, not against the relaxation**, below
  h = 2 (the ramped relaxation did not finish at h = 1); the three linear-load cases are
  verified against it to h = 0.25.
- **Multigrid at the surface** is good, not perfect: the V-cycle alone creeps from 11 to 25
  iterations over 1000 × more unknowns; GMRES absorbs it.
- Not done: 3D, a held (Dirichlet) part beside the free bodies (supported, untested here beyond
  the code path).

## 12. Local refinement near the surface

The mapped surface is second order, but on a uniform grid its accuracy costs the same
resolution everywhere. `crates/eulerian/src/refine.rs` refines a band around φ = 0 instead: a
quadtree stored level by level, each level a sparse uniform lattice (`Grid::sparse`: only the
active nodes are stored and walked, so the whole `Solid` and the mapped boundary run on a level
unchanged), the surface on the finest level only.

- **Regions.** The finest region is the cells of level L − 1 within `band` (default 6) finest
  cells of φ = 0, each split in four. Every coarser region contains the next one dilated by two
  of its own cells (proper nesting, which implies the 2:1 balance across edges and corners).
  A level assembles only its *leaf* cells (those not refined further); the base level keeps
  only the cells touching the body.
- **Roles.** On level l a node is interior (all its cells are leaves here), a *master* (on the
  interface with level l − 1, at a coarse position; updated here), *hanging* (on that
  interface between two masters), a *ghost* (two layers outside the region), or *held*
  (on the interface with level l + 1, or under it in a buffer of two cells; copied from there).
- **Hanging nodes.** Value: cubic along the interface through the two masters and the next
  nodes on the line (quadratic at a corner). Mass and force: half to each master. With the
  half-faces of the cells this keeps a constant stress free of net force (the patch test) for
  any interpolation exact on linear fields; linear values (the transpose of the sharing) give
  an O(h²) value error and a first-order stress error in the next cells (1.2 % at h₀ = 1 next
  to the interface against 0.29 % with cubic values).
- **Ghosts.** Biquadratic in the coarser level's advanced nodes and this level's masters (a
  one-sided 3 × 3 window where the centred one reaches under the finer level), linear in time
  between the coarse states before and after its step.
- **Fluxes and the register.** Each half-face of a leaf carries the flux between two of its
  corners, so every flux is computed once and given to both sides. A dual face that the
  interface halves is integrated over the leaf's half with the midpoint rule (the stress a
  quarter cell from the edge, ¾ of the edge's plus ¼ of the parallel edge's). The force of the
  coarse leaves on an interface node, and their solid area, complete the finer level's master
  there (the flux register): summed over every degree of freedom of every level, the impulses
  of the internal fluxes cancel to round-off (`coarse_fine_fluxes_conserve_momentum`: 1e-16 of
  Σ|impulse|).
- **Time.** Berger–Oliger subcycling: a step of level l, then two of level l + 1 with Δt/2,
  then the finer values copied onto the held nodes (velocities averaged over the two
  substeps, which puts them at the coarse half step). With leapfrog, the coarse kick at tⁿ
  stands for [tⁿ − Δt/2, tⁿ + Δt/2] and the fine kicks for [tⁿ − Δt/4, tⁿ + 3Δt/4]; holding
  Fⁿ in the register through both substeps lags the coupling by a quarter coarse step, a
  first-order error that a P wave shows as a reflection ∝ h/λ. The register is extrapolated
  to each fine kick, Fⁿ + s(Fⁿ − Fⁿ⁻¹); the extra impulse, (Δt/4)(Fⁿ − Fⁿ⁻¹) per step,
  telescopes (momentum in transit between the levels' clocks) and is reported by the audit.
- **Statics.** `Hierarchy::residual` evaluates the composite equilibrium residual (all levels
  at one time, ghosts and held nodes made consistent, the registers added): the operator a
  static solver needs, independent of the time loop.
- **Regridding.** `Hierarchy::regrid` (or every N coarse steps) re-flags from the current
  surface and moves ζ by interpolation from the finest old level holding each point (exact for
  linear maps) and v by the overlap of old and new control volumes, rasterised on quarter cells
  of the finer of the two hierarchies (every piece of old mass finds new nodes; Σ m v is kept
  for full cells; the cut cells next to the surface are rebuilt with the new band and move in
  the eighth digit).

Verification (`crates/eulerian/tests/refine.rs`; studies `examples/refine_cases.rs` and
`examples/refine_wave.rs`):

| check | result |
|---|---|
| patch tests: uniaxial stretch 1.25, simple shear 0.3, neo-Hookean disc, 1 and 2 levels | composite residual ≤ 1.1e-13 at every degree of freedom; at rest after 60 steps; σ exact to 1e-7 at every node of every level and at the surface |
| no refinement | the hierarchy on one level equals the mapped solid on the full grid to 1e-15 |
| control volumes | the composite solid area equals the uniform finest grid's to 1e-12 |
| momentum | Σ impulses over all degrees of freedom of a free, stressed body, 2 levels subcycled: 1e-16 of Σ\|impulse\| per coarse step (after the register's in-transit part) |
| rigid translation of a disc through 30 regrids (2 levels) | ξ to 6e-15, v to 2e-12, σ to 1.4e-11 |
| regrid | every piece of old mass finds new nodes; Σ m v changes by 1.2e-8 (the rebuilt cut cells) |

**A P wave across a refined strip** (periodic box, Gaussian pulse of half-width w coarse cells,
linear elastic, κ = 0; reflected amplitude in % of the incident one, entry / exit):

| w | 1 level | 2 levels | 3 levels | 1 level, register held at Fⁿ |
|---|---|---|---|---|
| 2 | 1.5 / 17 | 1.3 / 14 | 1.2 / 14 | 2.7 / 18 |
| 4 | 0.43 / 0.43 | 0.49 / 0.50 | 0.44 / 0.45 | 1.25 / 1.46 |
| 8 | 0.11 / 0.11 | 0.13 / 0.13 | 0.13 / 0.13 | 0.55 / 0.58 |
| 16 | 0.027 / 0.027 | 0.033 / 0.033 | 0.035 / 0.035 | 0.27 / 0.27 |

Second order in h/w (first order with the held register); every level at the finest Δt gives
0.57, 0.15, 0.037 % (1 level), so the subcycled coupling reflects no more than the spatial one.
At w = 2 the pulse is not resolved on the coarse grid and what the fine grid carries cannot
leave it (17 % at the exit).

**Stability.** Without artificial viscosity the coupled operator has slowly growing modes:
the pulse's amplitude reaches 80 × (1 level, after 5–7 trips round the box) to 300 × (2
levels) and saturates there; every level at the finest Δt grows too (10–40 × in 40 trips), so
it is the spatial coupling, not the subcycling. The ghost-cell coupling is conservative but not
variational (a fine stencil reads interpolated coarse values whose reaction no coarse node
receives), as the mapped surface already is (§10: its eigenvalues grow by 0.5 % per step
without damping). κ = 0.01 (1 level) and κ = 0.02 (2 and 3 levels) keep 40 trips decaying;
the default κ = 0.1 has a wide margin. The explicit step with κ = 0.1 relaxes the refined
Lamé cylinder at Δt = 0.6 h/c_P (the viscous cap is 0.625); with κ = 0.02 it relaxes at 0.5
and diverges at 0.6 (late, t = 40) and 0.7, where the uniform mapped scheme still relaxes.

**Statics.** Linear elastic, load × 0.01 (strain 1e-5), relaxed until the deformation speed
stays below 1e-6 of its scale for a full period 2π b/c_P (a single mode passes through zero
velocity everywhere at once, and stopping at the first quiet check left the free disc's
displacements 10 × less accurate). Errors in % of σθ(a) exact (Lamé, disc) or 3S (Kirsch).
Surface: the worst σθ at the bore's surface points (Kirsch: RMS of σθ around the hole);
interior: the worst nodal stress component at nodes with a + 3 < r < b − 3. CPU on one core
of a shared, heavily loaded machine (indicative only).

| grid | finest h | Lamé surface | interior | Kirsch surface | interior | disc surface | interior | updates (Lamé) | CPU s (Lamé) |
|---|---|---|---|---|---|---|---|---|---|
| uniform | 1 | 1.403 | 0.476 | 3.651 | 0.698 | 0.894 | 0.316 | 2.78e+06 | 6 |
| uniform | 0.5 | 0.316 | 0.132 | 0.535 | 0.194 | 0.201 | 0.086 | 2.31e+07 | 30 |
| uniform | 0.25 | 0.082 | 0.035 | 0.119 | 0.051 | 0.052 | 0.023 | 1.84e+08 | 158 |
| uniform | 0.125 | 0.027 | 0.010 | 0.027 | 0.013 | 0.017 | 0.006 | 1.47e+09 | 794 |
| uniform | 0.0625 | 0.009 | 0.003 | 0.007 | 0.003 | 0.006 | 0.002 | 1.18e+10 | 4299 |
| h₀ = 1, 1 level | 0.5 | 0.228 | 0.560 | 0.527 | 0.455 | 0.153 | 0.402 | 1.70e+07 | 30 |
| h₀ = 1, 2 levels | 0.25 | 0.074 | 0.779 | 0.115 | 0.459 | 0.030 | 0.620 | 7.36e+07 | 142 |
| h₀ = 1, 3 levels | 0.125 | 0.109 | 0.731 | 0.056 | 0.459 | 0.050 | 0.315 | 2.98e+08 | 544 |
| h₀ = 0.5, 1 level | 0.25 | 0.058 | 0.151 | 0.120 | 0.164 | 0.039 | 0.101 | 8.05e+07 | 141 |
| h₀ = 0.5, 2 levels | 0.125 | 0.040 | 0.149 | 0.038 | 0.168 | 0.021 | 0.099 | 3.05e+08 | 544 |
| h₀ = 0.5, 3 levels | 0.0625 | 0.076 | 0.154 | 0.043 | 0.165 | 0.040 | 0.102 | 1.21e+09 | 2066 |

- **The surface** reaches the uniform finest-grid accuracy while it is above the coarse floor,
  and then stops: at h₀ = 0.5 two levels give 0.040 / 0.038 / 0.021 % (Lamé / Kirsch / disc),
  uniform h = 0.125 0.027 / 0.027 / 0.017 %; a third level gives 0.076 / 0.043 / 0.040 %.
  At h₀ = 1 the Lamé and disc surfaces stall at the second level (0.074 → 0.109 %, 0.030 →
  0.050 %); Kirsch still gains at the third (0.115 → 0.056 %). The surface stress is set by the
  whole wall's equilibrium, and the coarse interior's error reaches it: around the Kirsch hole
  the refined runs carry a smooth cos 2θ offset of 0.05 % (h₀ = 1) and 0.025 % (h₀ = 0.5) on top
  of the point-to-point scatter of the finest grid.
- **The interior** keeps the error of h₀: 0.15 / 0.17 / 0.10 % at h₀ = 0.5 for every number of
  levels (uniform h₀ = 0.5: 0.13 / 0.19 / 0.086 %), against 0.0095 / 0.013 / 0.006 % on the
  uniform h = 0.125 grid and 0.003 / 0.003 / 0.002 % at 0.0625. At h₀ = 1 the refined bands fill
  most of the 10 mm wall and the interior maximum (0.7 %) sits on level-1 nodes next to the
  level-0 interface, about 1.5 × the coarse grid's own error. The payoff therefore exists only
  where the interior error at h₀ is below the target: for these cases h₀ = 0.5 holds 0.1–0.2 %
  in the interior, so refinement pays for a surface target between about 0.04 % and that level,
  and not below.

**Cost** in node updates (a degree of freedom advanced by one step of its level) to rest, and
in CPU time:

| grid | finest h | updates, uniform same h / refined (Lamé, Kirsch, disc) | CPU, same | uniform updates at the refined surface error / refined | CPU, same |
|---|---|---|---|---|---|
| h₀ = 1, 1 level | 0.5 | 1.4, 1.8, 1.4 | 1.0, 1.0, 1.0 | 2.2, 1.9, 2.1 | 1.5, 1.1, 1.4 |
| h₀ = 1, 2 levels | 0.25 | 2.5, 3.7, 2.5 | 1.1, 1.3, 1.1 | 3.0, 3.9, 6.7 | 1.3, 1.4, 2.4 |
| h₀ = 1, 3 levels | 0.125 | 4.9, 7.5, 5.0 | 1.5, 1.8, 1.5 | 0.4, 2.7, 0.7 | 0.2, 0.8, 0.3 |
| h₀ = 0.5, 1 level | 0.25 | 2.3, 3.0, 2.3 | 1.1, 1.3, 1.1 | 4.4, 3.0, 3.8 | 1.9, 1.3, 1.7 |
| h₀ = 0.5, 2 levels | 0.125 | 4.8, 7.0, 4.8 | 1.5, 1.8, 1.5 | 2.3, 4.4, 3.2 | 0.8, 1.2, 1.1 |
| h₀ = 0.5, 3 levels | 0.0625 | 9.7, 14.6, 9.8 | 2.1, 2.7, 2.1 | 0.2, 1.0, 0.2 | 0.1, 0.3, 0.1 |

At the same finest h the refined surface error is 0.6–1.0 × the uniform one with one level
(and two at h₀ = 1), 1.3–1.5 × with two levels on h₀ = 0.5, and 2–9 × once it has stalled;
the interior error is 2–77 × the uniform one (it is h₀'s).

An estimate made before it was built (two levels near both surfaces on h₀ = 0.5, a band of ±4
fine cells: about 6e4 node updates per coarse step against 4.4e5 for a uniform h = 0.125 grid,
roughly 7 ×) holds for the count: Lamé at h₀ = 0.5 with two levels makes 5.0e4 updates per
coarse step, the uniform h = 0.125 grid 2.4e5 (its 60300 degrees of freedom, not the 110889
nodes of the box), 4.8 × at the same finest h (Kirsch 7.0 ×). At equal surface accuracy the
saving in updates is 1.9–6.7 × before the stall and below 1 after it. In CPU time it is 0.8–2.4
× before the stall: the hierarchy walks 3–4 active nodes per degree of freedom on the refined
levels (level 2: 31584 active for 10972 degrees of freedom; level 1: 14320 for 1884: the outer
half of the band, ghosts, hanging and held nodes, the covered buffers), and every per-node loop
of the mapped solid (band fill, cut build, fluxes) visits them; per active node-step it costs
the same as the uniform grid (0.56–0.76 µs). A narrower band (4 finest cells instead of 6)
cuts 23 % of the updates but doubles the Kirsch surface error (RMS 0.069 % against 0.038 %).

Where it falls short:
- The interior error stays at h₀, so local refinement cannot reach a uniform accuracy target
  cheaper unless h₀ already meets it, and the surface stalls after one or two levels, at a
  floor set by the coarse interior (about 0.04 % here at h₀ = 0.5).
- CPU savings are small (≤ 2.4 ×) because inactive-but-walked nodes dominate the refined
  levels; skipping frozen and empty nodes in the solid's per-node loops, and not activating the
  outside of the body beyond the band, would bring the CPU ratio towards the update ratio.
- The coupling needs artificial viscosity (κ ≥ 0.02) for long undamped dynamics, and the
  stable Δt with little viscosity is somewhat lower than on a uniform grid.
- Static answers still come from dynamic relaxation with Δt set by the finest level's
  subcycles; for statics the composite residual (`Hierarchy::residual`) is the operator the
  Newton–Krylov solver of §11 should use, and its multigrid could reuse the levels.
- Not done: regridding a body that moves more than the band between regrids (the transfer
  needs the old band to hold the new surface), refinement on criteria other than the distance
  to the surface, 3D.

## 13. Threads and the GPU

**Threads.** Every loop of the explicit step over nodes, cells, faces and surface points runs on
`fieldcraft-parallel`'s pool (`map_indexed`: results written in place, at least a fixed grain
of indices per task). Each value is computed by the operations of the sequential loop and every
sum a node gathers (face fluxes in face order, dual-cell areas in the order of their cells, the
shares of control volumes and forces) is formed in the sequential scatter's order, so the step
gives the sequential numbers bit for bit on any thread count (`tests/threads.rs`: 1, 2 and 7
threads, mapped, blurred and ghost-filled; checked also against checksums of the sequential
code for the static solve and the refined hierarchy). The surface chains are solved one after
another; within a chain the residuals and the Jacobian's columns run in parallel.

One explicit step of the Lamé cylinder by phase (`examples/rmt_profile.rs`; Apple M4 Max, ms,
one thread):

| h | nodes | interior faces | mapped faces | surface forces | momentum | advection | band | cut cells | surface points | Newton | step |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.125 | 160801 | 6.1 | 3.7 | 3.0 | 3.2 | 0.9 | 5.2 | 7.9 | 0.5 | 2.6 | 33.1 |
| 0.0625 | 641601 | 23.8 | 8.3 | 11.6 | 11.8 | 2.8 | 16.2 | 32.4 | 1.2 | 4.7 | 112.7 |

Speed-up of the whole step over one thread (10 performance and 4 efficiency cores):

| threads | 2 | 4 | 8 | 10 | 14 |
|---|---|---|---|---|---|
| h = 0.125 | 1.53 | 1.89 | 2.12 | 2.03 | 1.95 |
| h = 0.0625 | 1.72 | 2.64 | 3.34 | 3.31 | 3.23 |

The pure stencil phases scale 6–7 × on 8 threads (interior faces, advection), the surface
ones barely: the Newton chains not at all, the mapped faces and surface points 1–2 ×, and
the band and cut-cell scans, which walk the whole grid to find the few nodes near the surface,
1.6–2.3 × and 2.5–3.9 ×. These measurements were made on the mapped step before it was
merged with the static solve's held-topology structures; the merged step gives the same
numbers bit for bit and costs about 18 % more per step (38–39 against 32–33 ms at h = 0.125
on one thread, measured back to back): the static solve's per-topology structure, its face
list with the compact stencils' nodes, is rebuilt every explicit step.

**The deep interior on the GPU** (`crates/eulerian-gpu`, layer 7, research only). Far from the
surface the mapped step is a pure stencil on full cells (the compact edge stress, the momentum
of a full control volume, the upwind advection of ζ). `Solid::set_remote` lets a device take
the solid nodes more than six rings inside (`CPU_RINGS`: the band's extrapolation windows, the
surface map and the cut cells' corners need the CPU's own values that deep); `GpuSolid` steps
them in a wgpu compute shader while the CPU steps everything else at the same time, and each
step exchanges the values within two rings of the seam. With no surface (a periodic wave)
every node is on the device and steps go in chunks of 64 without exchange.

Single precision needs care: the stress formed as written loses the strain to cancellation
(F ≈ I). The shader forms it from small quantities: H = F − I = −F G with G = ∇ζ, B − I =
H + Hᵀ + H Hᵀ, ln J = −log1p(tr G + det G) by its series. Lamé cylinder at h = 0.25 after
8000 steps, largest stress difference from the CPU (double precision) over p, and the largest
stress error against the closed form over the nodes 2h inside:

| load (strain) | GPU, stable stress | GPU, stress as written | error vs closed form: CPU / stable / as written |
|---|---|---|---|
| p = 100 (1e-3) | 7.7e-6 | 3.1e-4 | 5.25e-3 / 5.25e-3 / 5.30e-3 |
| p = 1 (1e-5) | 1.2e-5 | 2.9e-2 | 1.49e-3 / 1.49e-3 / 2.89e-2 |

The device's error is two to three orders below the discretisation error (`tests/agreement.rs`:
a periodic shear wave, the pressurised ring, a disc translating across the grid while the
split follows it, against the CPU; skipped without an adapter unless `FIELDCRAFT_REQUIRE_GPU=1`).

Throughput (`examples/rmt_gpu.rs`; Apple M4 Max; ms per step):

| shear wave, nodes | 4096 | 16384 | 65536 | 262144 | 1048576 | 3686400 |
|---|---|---|---|---|---|---|
| CPU, 1 thread | 0.84 | 3.3 | 13.3 | 55.9 | 240 | 846 |
| CPU, 10 threads | 0.78 | 1.8 | 3.1 | 11.0 | 39.4 | 150 |
| GPU, every node | 0.034 | 0.041 | 0.063 | 0.15 | 0.53 | 1.92 |
| GPU, 10⁹ node updates/s | 0.12 | 0.40 | 1.0 | 1.7 | 2.0 | 1.9 |

The interior stencil runs 70–80 × faster on the GPU than on ten CPU threads from a quarter of a
million nodes up, at 1.7–2.0 × 10⁹ node updates per second (the lattice Boltzmann GPU backend
on the same machine: 0.87–1.31 × 10⁹ D3Q19 cell updates per second).

With a surface it does not pay yet. Lamé cylinder, ten CPU threads (ms per step):

| h | nodes | nodes on the GPU | CPU alone | GPU interior + CPU surface | CPU part of it | wait for the GPU |
|---|---|---|---|---|---|---|
| 0.25 | 40401 | 7364 | 8.4 | 8.5 | 8.3 | 0.05 |
| 0.125 | 160801 | 44928 | 17.9 | 16.6 | 16.4 | 0.07 |
| 0.0625 | 641601 | 210544 | 36.2 | 32.9 | 32.5 | 0.12 |
| 0.03125 | 2563201 | 903620 | 93.5 | 83.7 | 82.9 | 0.24 |

The device finishes long before the CPU (the wait is under 0.3 ms), so the hybrid step is the
CPU's surface work: 1 % slower than the CPU alone at h = 0.25 in the table's run (1 % faster in
a second), 2–7 %, 8–9 % and 10–13 % faster at h = 0.125, 0.0625 and 0.03125 over the two
runs. Marking the deep interior remote
leaves the CPU step at 16–17 / 33–35 ms on ten threads (h = 0.125 / 0.0625), the same as with
every node on the CPU (16.3 / 34.1): what remains are the band, the cut cells and the faces
next to the surface, and the full-grid scans that find them (widening the grid around the same
cylinder from 160801 to 1442401 nodes takes the CPU step from 16 to 41 ms with or without the
remote interior).

What it would take, in order: active-node lists and persistent buffers so the CPU's work
scales with the surface rather than the grid; the cut cells and mapped faces on the device;
the surface chains last (their Newton is sequential along a chain).

## 14. Coupling to the lattice, and the momentum risk

In the blurred LBRMT the fluid stress never appears: the solid's force density
∇·[(1 − H) τ_s] goes to the LB collision, and momentum is conserved because the force density
is the divergence of one stress field (the edge sums telescope). With a sharp solid surface the
solid needs the fluid's traction σ_f·n on the surface: the pressure from the density and the
viscous stress from the non-equilibrium moments (σ'_f = −(1 − 1/(2τ)) Π^neq, local at every
lattice node), at the surface's Gauss points; that is the ambient stress field the mapped
surface integrates now (`Solid::set_ambient_field`). The fluid must then see the solid as a
moving wall. The risk is momentum: the solid's momentum change is the balance of its own
control volumes, and the fluid's comes from its own boundary scheme; if the two compute the
exchanged force separately they differ by a truncation error, a free body drifts and the total
momentum is not conserved to round-off as in the blend. The cut cells hold what a conservative
coupling needs: the exact solid fraction of every cell (the input of partially saturated cells,
Noble & Torczynski 1998) and the surface pieces of every cut cell (the links a momentum
exchange crosses, Ladd 1994; Aidun et al. 1998). Computing the exchanged force once and
applying it to both sides makes the coupling conservative by construction. Not measured here
(no fluid yet).

## 15. What still falls short

- **For a solid alone the FEM wins.** The static solve (§11) is about 200 times slower than
  FieldCraft's FEM at equal accuracy (0.18 s against about 1 ms for the Lamé cylinder). The
  RMT's case is FSI with large solid motion on the fluid's grid, not structural analysis.
- **FSI blockers.**
  - *The surface is on the CPU.* The deep interior runs at GPU speed (§13), but the surface
    work (band, cut cells, surface map, Newton chains) and the full-grid scans that find it set
    the step; a lattice on the GPU with a solid surface on the CPU would wait on the CPU every
    step.
  - *Stiff solids.* The explicit solid step is bounded by 0.4 h/c_P, while the LBM steps at the
    lattice's own (low-Mach) speed scale: for an engineering metal in a slow flow the solid's
    stable step is hundreds of times below the lattice's (about 700 times in the estimate made
    for the GPU study). A soft solid (the LBRMT's examples) avoids it; a stiff one needs an
    implicit or subcycled solid step (the static solve's Newton–Krylov machinery is the
    starting point) or a rigid-body treatment.
- **Line-search stalls after topology changes** in the finite-strain static solve (§11): the
  increment is halved 1–2 times per solve at h ≤ 0.125, most likely at surface points a node
  only just behind the surface.
- **Refinement pays modestly** (§12): the interior error stays at the coarse spacing.
- **The ghost-filled and mapped surfaces are not symmetric** and need artificial viscosity
  (κ ≥ 0.02–0.1) for long undamped dynamics.
- Not done: 3D; contact between bodies; prescribed moving boundaries; the incompressible law the
  LBRMT uses; thin features below about two cells; the LB coupling itself.

## 16. Related work and what appears new

A literature search (Crossref and the open literature to October 2026, about 110 papers) found
the following. It is a search, not a proof: some relevant texts could not be read (the Valkov
MIT thesis, Lin's Harvard thesis, both Lara & Kamrin CMAME papers on RMT contact, the full text
of Richter 2013), and they may contradict the claims below.

**Closest published methods.**
- *Sharp traction on Cartesian grids, small strain*: cut-cell finite volumes with the traction
  integrated over the true boundary segment (Theillard et al. 2013) and ghost points with
  σ·n = t imposed at the closest boundary point plus geometric multigrid (Coco et al. 2014),
  both second order in displacement. They are the closest analogues of the mapped cut cells and
  the static solve, at small strain. The explicit simplified interface method (Lombard &
  Piraux 2004/2006) builds the interface condition into a least-squares extension and keeps
  second order only with extension order ≥ 3, as our ghost fill found for its cubic fit.
- *Finite-strain Eulerian solids with sharp interfaces* (the Godunov line of Barton, Drikakis,
  Pullin; Gorsse, Iollo & Milcent's backward characteristics, which are the reference map)
  impose interface conditions through Riemann ghost states but report no interface order.
- *Fully Eulerian FSI*: the finite element line (Dunne 2006; Richter & Wick 2010; Richter 2013;
  Frei & Richter 2014) defines the solid by the Initial Point Set, which is the reference map
  under another name, and solves it implicitly by Newton with multigrid; it reaches optimal
  order only by fitting the interface (locally modified parametric elements, the closest
  precedent for our surface-fitted map). The finite-difference line (Sugiyama et al. 2011;
  the RMT papers of Kamrin, Rycroft, Jain, Lin, Wang, Sun) is explicit with a blended band and
  first order at the interface: Sugiyama et al. measure it against an exact solution, the LBRMT
  reports about 1.4 in self-convergence, Wang et al. 2022 show second order only towards the
  model with a frozen blur width.
- *Unfitted finite elements*: merged control volumes ≈ aggregated unfitted FEM (Badia et al.
  2018); the ghost-cell preconditioner matrix ≈ in spirit the ghost penalty (Burman 2010);
  cut-cell-aware coarsening ≈ immersed multigrid with Schwarz smoothers (de Prenter et al.
  2020); the shifted boundary method (Main & Scovazzi 2018; Atallah et al. 2021) is the
  strongest alternative to cut cells, optimal for linear elasticity, with no finite-strain
  version found.
- *Static equilibrium in the inverse motion* on body-fitted meshes: Govindjee & Mihalic 1996.
- *GPU*: no GPU implementation of the RMT, the LBRMT, the full-Eulerian method or a Godunov
  Eulerian solid was found; the only GPU Eulerian solids are unvalidated graphics papers
  (Levin et al. 2011; Fan et al. 2013).

**What appears new (hedged).**
1. The stiffness bias of the linear least-squares extrapolation of ξ, measured against an
   analytic static solution (§8). The ingredients are published (first-order ξ in the band,
   Valkov et al. 2015; a grid-independent det F error of about 14 % in the band, Jain et al.
   2019; second order only towards the frozen-blur model, Wang et al. 2022), but we found no
   paper linking them to a stiffness error, and every published extrapolation test uses an
   affine ξ, for which a linear fit is exact. The Valkov thesis is the most likely place for an
   earlier static test.
2. The dependence of a blurred RMT's static answer on the damping, through band nodes that
   integrate momentum while their ξ is overwritten (§8). Not found.
3. An implicit static RMT: Newton–Krylov with geometric multigrid and h-independent iteration
   counts (§11). New for the RMT as far as we can tell, but not the first implicit Eulerian
   solid (Richter & Wick 2010; Richter 2013) nor the first static solve in the inverse motion
   (Govindjee & Mihalic 1996).
4. A measured second-order surface stress for a finite-strain reference-map solid with sharp
   traction on a fixed grid (§10, §11). The individual techniques have precedents (above); their
   combination in an RMT with the traction on the true surface, and the measured order at
   finite strain, appear unreported.
5. Verification against exact solutions instead of self-convergence: the blurred RMT is −11 %
   at ten nodes across a wall while its self-convergence looks healthy. Wang et al. 2022's
   frozen-blur appendix comes closest to saying so.
6. An RMT solid on a GPU (§13), with the caveat that only its deep interior runs there.

## 17. Next steps

In order of payoff for the goal (FSI with large solid motion on the GPU lattice):

1. **The LB coupling** (step 3 of #594): monolithic and conservative through the cut-cell data
   (§14): partially saturated cells driven by the cut cells' solid fractions, or a momentum
   exchange on the cut links with covered and uncovered nodes accounted for (Wen et al. 2014's
   Galilean-invariant form for moving walls), one velocity field as in the LBRMT. Risk: how a
   sharp solid traction and a one-cell-diffuse PSM force agree at second order is unproven;
   the fallback is LBRMT-style forcing on the true solid fraction (first order, still
   conservative).
2. **The stiff-solid time step**: an implicit or subcycled solid step built on the static
   solve's Newton–Krylov and multigrid, or a rigid-body limit for very stiff parts (Wang et
   al. 2022's projection; Millet & Liberge 2026).
3. **The surface on the GPU**: active-node lists so the CPU's work scales with the surface;
   the cut cells and mapped faces on the device; then the chains. A head-to-head with the
   shifted boundary method (whole-cell stencils, no cut integration, no chains) is the
   alternative if the chains do not port.
4. **A condition-built, higher-order extension at the surface** (in the spirit of Lombard &
   Piraux and Coco et al.): the traction built into a quadratic or cubic fit replacing the 1/s
   weights near nodes just behind the surface, to remove the line-search stalls of §11 and give
   a band accurate for dynamics. Schwarz-patch smoothing (de Prenter et al. 2020) for the
   V-cycle's creep.
5. **3D** once the coupling works in 2D.

## Papers

All read from the open-access versions listed; journal details checked against the published
DOI records (Crossref) where marked.

1. K. Kamrin, J.-C. Nave, *An Eulerian approach to the simulation of deformable solids:
   Application to finite-strain elasticity*, arXiv:0901.3799 (2009).
   <https://arxiv.org/abs/0901.3799> — the first RMT paper (journal version not checked).
2. K. Kamrin, C. H. Rycroft, J.-C. Nave, *Reference map technique for finite-strain elasticity
   and fluid–solid interaction*, J. Mech. Phys. Solids 60(11):1952–1969 (2012).
   doi:[10.1016/j.jmps.2012.06.003](https://doi.org/10.1016/j.jmps.2012.06.003) (verified).
   Read from the author copy <https://www.math.mcgill.ca/jcnave/papers/elasticity.pdf>.
3. B. Valkov, C. H. Rycroft, K. Kamrin, *Eulerian method for multiphase interactions of soft
   solid bodies in fluids*, J. Appl. Mech. 82(4):041011 (2015).
   doi:[10.1115/1.4029765](https://doi.org/10.1115/1.4029765) (verified). Preprint
   arXiv:1409.6183, titled *Eulerian method for fluid–structure interaction and submerged
   solid–solid contact problems*.
4. C. H. Rycroft, C.-H. Wu, Y. Yu, K. Kamrin, *Reference map technique for incompressible
   fluid–structure interaction*, J. Fluid Mech. 898:A9 (2020).
   doi:[10.1017/jfm.2020.353](https://doi.org/10.1017/jfm.2020.353) (verified). Preprint
   arXiv:1810.03015.
5. Y. Sun, C. H. Rycroft, *A fully-integrated lattice Boltzmann method for fluid–structure
   interaction*, J. Comput. Phys. 526:113774 (2025).
   doi:[10.1016/j.jcp.2025.113774](https://doi.org/10.1016/j.jcp.2025.113774) (verified).
   Read from the preprint arXiv:2402.12696 (2024), <https://arxiv.org/abs/2402.12696>; section
   and equation numbers cited in the code are the preprint's.
6. Y. Sun, *A lattice Boltzmann implementation of the reference map technique*, master's
   thesis, Harvard University (2020).
   <https://people.math.wisc.edu/~chr/publish/2020_Sun.pdf>.

For the sharp boundaries (§9, §10), journal records checked against Crossref (metadata only;
the methods were written from the equations, these are pointers to the established ideas, not
sources of code or text):

7. R. Mittal, H. Dong, M. Bozkurttas, F. M. Najjar, A. Vargas, A. von Loebbecke, *A versatile
   sharp interface immersed boundary method for incompressible flows with complex
   boundaries*, J. Comput. Phys. 227(10):4825–4852 (2008). doi:10.1016/j.jcp.2008.01.028
   (verified): ghost cells with the boundary condition at the boundary intercept.
8. P. T. Barton, D. Drikakis, *An Eulerian method for multi-component problems in non-linear
   elasticity with sliding interfaces*, J. Comput. Phys. 229(15):5518–5540 (2010).
   doi:10.1016/j.jcp.2010.04.012 (record verified; its ghost-fluid treatment not checked).
9. T. D. Aslam, *A partial differential equation approach to multidimensional extrapolation*,
   J. Comput. Phys. 193(1):349–355 (2004). doi:10.1016/j.jcp.2003.08.001 (verified).
10. H. Johansen, P. Colella, *A Cartesian grid embedded boundary method for Poisson's equation
    on irregular domains*, J. Comput. Phys. 147(1):60–85 (1998). doi:10.1006/jcph.1998.5965
    (verified): conservative cut cells.
11. C. S. Peskin, *The immersed boundary method*, Acta Numerica 11:479–517 (2002).
    doi:10.1017/S0962492902000077 (verified).

Related work (§16; DOIs checked on Crossref by the literature search; read in full where open
copies existed, else from abstracts):

12. M. Theillard, L. Djodom, J. Vié, F. Gibou, *A second-order sharp numerical method for
    solving the linear elasticity equations on irregular domains and adaptive grids –
    Application to shape optimization*, J. Comput. Phys. 233:430–448 (2013).
    doi:10.1016/j.jcp.2012.09.002.
13. A. Coco, G. Currenti, C. Del Negro, G. Russo, *A second order finite-difference ghost-point
    method for elasticity problems on unbounded domains with applications to volcanology*,
    Commun. Comput. Phys. 16(4):983–1009 (2014). doi:10.4208/cicp.210713.010414a.
14. B. Lombard, J. Piraux, *Numerical treatment of two-dimensional interfaces for acoustic and
    elastic waves*, J. Comput. Phys. 195(1):90–116 (2004). doi:10.1016/j.jcp.2003.09.024.
15. S. Jain, K. Kamrin, A. Mani, *A conservative and non-dissipative Eulerian formulation for
    the simulation of soft solids in fluids*, J. Comput. Phys. 399:108922 (2019).
    doi:10.1016/j.jcp.2019.108922.
16. Y. Lin, N. Derr, C. Rycroft, *Eulerian simulation of complex suspensions and
    biolocomotion in three dimensions*, Proc. Natl. Acad. Sci. 119:e2105338118 (2022).
    doi:10.1073/pnas.2105338118.
17. X. Wang, K. Kamrin, C. Rycroft, *An incompressible Eulerian method for fluid–structure
    interaction with mixed soft and rigid solids*, Phys. Fluids 34:033604 (2022).
    doi:10.1063/5.0082233.
18. K. Sugiyama, S. Ii, S. Takeuchi, S. Takagi, Y. Matsumoto, *A full Eulerian finite
    difference approach for solving fluid–structure coupling problems*, J. Comput. Phys.
    230(3):596–627 (2011). doi:10.1016/j.jcp.2010.09.032.
19. T. Richter, T. Wick, *Finite elements for fluid–structure interaction in ALE and fully
    Eulerian coordinates*, Comput. Methods Appl. Mech. Eng. 199:2633–2642 (2010).
    doi:10.1016/j.cma.2010.04.016.
20. T. Richter, *A fully Eulerian formulation for fluid–structure-interaction problems*, J.
    Comput. Phys. 233:227–240 (2013). doi:10.1016/j.jcp.2012.08.047.
21. S. Frei, T. Richter, *A locally modified parametric finite element method for interface
    problems*, SIAM J. Numer. Anal. 52(5):2315–2334 (2014). doi:10.1137/130919489.
22. S. Govindjee, P. Mihalic, *Computational methods for inverse finite elastostatics*,
    Comput. Methods Appl. Mech. Eng. 136:47–57 (1996). doi:10.1016/0045-7825(96)01045-6.
23. S. Badia, F. Verdugo, A. Martín, *The aggregated unfitted finite element method for
    elliptic problems*, Comput. Methods Appl. Mech. Eng. 336:533–553 (2018).
    doi:10.1016/j.cma.2018.03.022.
24. E. Burman, *Ghost penalty*, C. R. Math. 348:1217–1220 (2010). doi:10.1016/j.crma.2010.10.006.
25. F. de Prenter, C. Verhoosel, E. van Brummelen, J. Evans, C. Messe, J. Benzaken, K. Maute,
    *Multigrid solvers for immersed finite element methods and immersed isogeometric analysis*, Comput. Mech. 65:807–838 (2020). doi:10.1007/s00466-019-01796-y.
26. A. Main, G. Scovazzi, *The shifted boundary method for embedded domain computations. Part
    I: Poisson and Stokes problems*, J. Comput. Phys. 372:972–995 (2018).
    doi:10.1016/j.jcp.2017.10.026.
27. N. Atallah, C. Canuto, G. Scovazzi, *The shifted boundary method for solid mechanics*,
    Int. J. Numer. Methods Eng. 122(20):5935–5970 (2021). doi:10.1002/nme.6779.
28. D. Noble, J. Torczynski, *A lattice-Boltzmann method for partially saturated
    computational cells*, Int. J. Mod. Phys. C 9(8):1189–1201 (1998).
    doi:10.1142/s0129183198001084.
29. A. Ladd, *Numerical simulations of particulate suspensions via a discretized
    Boltzmann equation. Part 1. Theoretical foundation*, J. Fluid Mech. 271:285–309 (1994).
    doi:10.1017/s0022112094001771.
30. C. Aidun, Y. Lu, E. Ding, *Direct analysis of particulate suspensions with inertia
    using the discrete Boltzmann equation*, J. Fluid Mech. 373:287–311 (1998).
    doi:10.1017/s0022112098002493.
31. B. Wen, C. Zhang, Y. Tu, C. Wang, H. Fang, *Galilean invariant fluid–solid interfacial
    dynamics in lattice Boltzmann simulations*, J. Comput. Phys. 266:161–170 (2014).
    doi:10.1016/j.jcp.2014.02.018.
32. D. Levin, J. Litven, G. Jones, S. Sueda, D. Pai, *Eulerian solid simulation with
    contact*, ACM Trans. Graph. 30(4) (2011). doi:10.1145/2010324.1964931.
33. Y. Fan, J. Litven, D. Levin, D. Pai, *Eulerian-on-Lagrangian simulation*, ACM
    Trans. Graph. 32(3) (2013). doi:10.1145/2487228.2487230.
34. T. Millet, E. Liberge, *A monolithic lattice Boltzmann–fictitious domain method for particle
    sedimentation*, Results Eng. 30:110942 (2026). doi:10.1016/j.rineng.2026.110942.

Background cited by these papers and used in the prototype: J. Bonet, R. D. Wood, *Nonlinear
Continuum Mechanics for Finite Element Analysis*, 2nd ed., Cambridge 2008 (the compressible
neo-Hookean law σ = (μ/J)(B − I) + (λ ln J / J) I); T. D. Aslam, *A partial differential
equation approach to multidimensional extrapolation*, J. Comput. Phys. 193(1):349–355 (2004).
