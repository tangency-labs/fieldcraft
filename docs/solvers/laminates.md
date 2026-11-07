# Laminate kernels

The constitutive and MITC kernels are implemented in `fem::laminate`, `fem::failure` and
`fem::shell`. Integration with model sections, commands and solver decks is in progress;
the tests below are kernel tests, not a completed command oracle or cross-program validation.

## Coordinates and thickness

The node order defines the shell's right-hand normal. Plies are listed **bottom to top**, from
the negative to the positive normal. A rectangular material frame uses its first direction
and a reference direction in the 1–2 plane. A cylindrical frame is radial, circumferential,
axial about its specified axis. A shell projects an explicitly selected frame axis into its
plane at the element centre, then rotates by the ply angle. Use the hoop or axial reference
on a cylindrical surface; a radial reference parallel to its normal is an error. The frame is
constant within each facet. Angles stored by the model and FEM kernels are radians.

With total thickness h and fractional offset o, the laminate centre is e = −oh, and its
bottom is z = e − h/2 relative to the nodal surface. The kinematics remain those of the existing
MITC4/MITC3: ε(z) = ε₀ + zκ, with engineering shear γ12 = 2ε12. In-plane vectors are ordered
11,22,12; full vectors are 11,22,33,12,13,23. Tension is positive.

Following [Nettles, NASA RP-1351 (1994), §§II–IV](https://ntrs.nasa.gov/citations/19950009349),
each lamina is reduced to plane stress. With T the engineering-strain transformation from
shell to material axes, Qbar = TᵀQT. A, B and D are its zeroth, first and second thickness
moments. Homogeneous ply moments are integrated analytically, including each ply's t³/12
term, so using one section sample does not remove bending stiffness. One to fifteen
Gauss–Legendre section samples supplement each ply's bottom, middle and top fields.

The two transverse shear directions have a full rotated 2×2 stiffness matrix, with the
existing 5/6 correction. Ply densities supply I0 = ∫ρdz, I1 = ∫ρz dz and I2 = ∫ρz²dz.
Consistent mass and gravity retain the I1 translation–rotation coupling. Diagonal lumping
omits I1, as it already does for a homogeneous offset shell. Missing density is an error for
mass or gravity, not zero mass.

## Fields and failure

Ply fields remain separate on either side of an interface; they are not averaged across
materials. They include plane-stress ε33 from the normal compliance. The transverse stress
entries are the effective FSDT values (5/6 Gγ). They do not constitute an equilibrated
interlaminar stress recovery, do not vanish at free surfaces, and must not be used to claim
agreement with a three-dimensional interlaminar peak. Homogeneous-shell stress output retains
its older parabolic shear reconstruction; equivalent isotropic laminates agree in stiffness
and mass, but intentionally use the documented laminate field convention.

Strengths Xt, Xc, Yt, Yc, S12 and S23 are positive magnitudes. The four criteria apply to
material-axis σ11, σ22, τ12, with failure at an index of 1:

- Maximum stress is the largest stress/strength magnitude, choosing tension or compression
  strength by the corresponding normal stress sign.
- Tsai–Hill is the signed-product polynomial `(σ1/X)^2 − σ1σ2/X^2 + (σ2/Y)^2 + (τ12/S12)^2`,
  where X and Y follow the respective stress signs. It is not Azzi–Tsai–Hill.
- Tsai–Wu uses tensile/compressive linear terms and the quadratic interaction
  `F12 = f12 sqrt(F11 F22)`. The input f12 is strictly between −1 and 1. Its default −0.5 is
  a stated modelling assumption. The polynomial is retained even when it is negative.
- Hashin uses the 1980 plane-stress equations with α=1 and squared fibre compression;
  the inactive sign branches are absent. Matrix compression uses the separately supplied
  S23. No unreported substitution of S12 is made.

These are first-ply failure indicators, without degradation or progressive damage. Quadratic
indices are neither square roots nor inverse reserve factors. References:
[Tsai & Wu (1971)](https://doi.org/10.1177/002199837100500106),
[Hashin (1980)](https://doi.org/10.1115/1.3153664), and the public
[failure-measures](https://abaqus.uclouvain.be/English/SIMACAEMATRefMap/simamat-c-failuremeasures.htm)
and [Hashin](https://abaqus.uclouvain.be/English/SIMACAEMATRefMap/simamat-c-damageinitfibercomposite.htm)
reference manuals. Only documentation and published numbers are used.

## Verification and limits

`cargo test -p fieldcraft-model -p fieldcraft-fem` checks closed-form ABD for [0/90]s,
[±45]s and [0/90], off-axis work conjugacy, asymmetric density and offset moments,
one-sided fields, failure thresholds, hostile inputs and coordinate objectivity. A balanced
symmetric angle-ply stack has zero A16/A26 and B; it need not have zero D16/D26.

`fem/tests/laminated_shell.rs` independently assembles a square plate under a consistently
integrated sinusoidal load. The symmetric cross-ply case uses hard simply supported edges,
and compares to the 3×3 Navier FSDT amplitude system with the same shear correction. Its thin
limit also checks CLT. Refining 4×4 to 8×8 must reduce the displacement error by more than a
factor of 1/0.35; the fine-grid error must be below 3%.

The same test file constructs NAFEMS R0093/1 from its published physical description:
[Hopkins (2005), public reproduction, Tables 1–3 and 8–12](https://doc.comsol.com/6.3/doc/com.comsol.help.models.compmat.failure_prediction_in_a_laminated_composite_shell/failure_prediction_in_a_laminated_composite_shell.html).
It solves the prescribed single quad and compares its centre (the mean of the extrapolated
nodal values) at each ply midplane, not every nodal field. The predeclared tolerance is 1%
per stress component plus 0.0005 MPa for published rounding. The independent failure test
uses the published material-axis stresses and a 0.0002 absolute FI tolerance.

`Shell::plate_geometric_stiffness` supplies only the small-slope, flat-plate initial-membrane-
stress term. It uses the preload at stiffness integration points, with compression negative.
It is not a curved-shell or finite-rotation tangent and includes no follower-load stiffness.
The solver must enforce a planar pre-buckling state before exposing this kernel. See
[MIT 2.080J, Lecture 11, §§11.1–11.2](https://ocw.mit.edu/courses/2-080j-structural-mechanics-fall-2013/345c462fe40421ea0f795bd35d2e4af1_MIT2_080JF13_Lecture11.pdf)
for the flat-plate stability equations (the lecture uses the opposite compression sign).

The new calculations have fixed iteration and accumulation order. Thread-count and platform
claims must be checked separately; these local tests alone are not a cross-platform bitwise
determinism certificate.
