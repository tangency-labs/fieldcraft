# Eight-node curved MITC8

`quad8` with an elastic `shell` or `composite_shell` section and **full** integration uses
`fem::mitc8::Mitc8` in static, frequency and harmonic solves. This is an eight-node serendipity shell, not MITC9:
four corners counterclockwise from (-1,-1), followed by the four edge midpoints. There
is no centre node, internal centre displacement, or condensation. Each node shares the
existing global `[ux, uy, uz, rx, ry, rz]` structural layout.

## Formulation and geometry

The implementation follows Bathe and Dvorkin (1986), *A formulation of general shell
elements—the use of mixed interpolation of tensorial components*, IJNME 22:697–722,
[doi:10.1002/nme.1620220312](https://doi.org/10.1002/nme.1620220312), Fig. 7 and equations
(9)–(22). The original paper numbers its corners from (+1,+1); the implementation uses
the model's corner order throughout, including its tying tensors.

| Formula | Code | Meaning |
|---|---|---|
| BD86 (9); Lee–Bathe 2004 (1),(2) | `Basis::Quad8`, `Mitc8::map`, `displacement_gradient` | quadratic surface, nodal directors, director rotations |
| Lee–Bathe (3)–(8); BD86 (10),(12) | `covariant`, `in_layer` | covariant small strain and its physical tensor |
| BD86 (11)–(18) | `assumed_layer` | eight scaled in-layer tying tensors; edge reconstruction in orthogonalized covariant bases |
| BD86 (19)–(22) | `shear` | edge shear tensors and centre dyad times averaged RA/RB or SA/SB covariant shear |
| BD86 §4.3 | `stiffness` | 3 x 3 area integration; three thickness Gauss points |
| BD86 §4.4 | `stresses` | direct displacement-strain nodal stress samples |

Lee and Bathe (2004), *Development of MITC isotropic triangular shell finite elements*,
Computers & Structures 82:945–962,
[doi:10.1016/j.compstruc.2004.02.004](https://doi.org/10.1016/j.compstruc.2004.02.004),
provides the general geometry/strain description; its triangular MITC6 interpolation
is not substituted for the eight-node formulation.

At node i the director is the normalized cross product of the reference surface's
quadratic tangents. With h the constant section thickness and e=-offset*h,
`x(r,s,t) = sum N_i(r,s) [x_i + (e+h*t/2) n_i]`. Director vectors are interpolated
without renormalization inside the element. Their variation is `theta_i cross n_i`.
This preserves infinitesimal rigid motion even on a curved surface and with offset.
The solver aligns incident nodal normals before averaging and normalizing them. Each
shared node then has one director axis, signed along each element's own normal. All
nodes retain global rotation vectors: `delta n = theta cross n`, so a shared quadratic
edge has the same thickened position and displacement on both sides at a common
physical thickness position. This follows the nodal-normal discussion in Lee–Bathe,
*Computers & Structures* 83 (2005), 69–90, section 3,
[doi:10.1016/j.compstruc.2004.07.005](https://doi.org/10.1016/j.compstruc.2004.07.005).
It does not replace the eight-node BD86 strain interpolation with MITC9.

The implementation rejects incident normal axes separated by 60 degrees or more
and mixed-order shell junctions. Common edges must have the same midside node;
nonconforming and nonmanifold edge connections are rejected. It is a smooth-patch construction, not a folded
junction model. The normal-axis screen is local: it cannot distinguish every
fold-back or overlapping patch from a reversed connectivity, and is not a global
self-intersection or orientability certificate. Supply a valid smooth surface mesh. Reversing connectivity reverses the signed normal. To describe the
same asymmetric physical laminate, also reverse ply order, negate offset and ply
angles, and keep the physical section reference direction fixed. Rotating local
corner numbering preserves an explicitly oriented material; an unspecified material
orientation intentionally follows the element r direction. Different thicknesses or
offsets on adjacent sections do not imply matching outer faces.

The Jacobian uses the actual three covariant bases through the thickness. Zero,
nonfinite or negative scaled determinants are rejected. Both thickness faces, all
nodes, the centre, quadrature and tying points are checked. Nodal normals turning
90 degrees or more from the centre normal are rejected; refine such a patch. These
sampled checks are not a proof of positivity everywhere in an arbitrarily distorted map.

Layer strain tensors carry membrane and bending together through the thickness.
At each integration point they are projected onto an orthonormal frame whose normal
is the interpolated director. Plane stress eliminates the normal stress; engineering
shears include the factor two. Transverse shear uses the fixed 5/6 correction of the
existing shell model. Thin-shell locking is checked independently of patch accuracy.
The section `points` value remains a deck sampling parameter; this elastic kernel
uses three thickness Gauss points per ply regardless of that value. These integrate
flat elastic layer moments exactly; curved-shell integration remains numerical.

## Sixth rotation and scope

The paper has five physical unknowns per node. The extra drilling rotation in this
six-DOF interface is tied to the reference surface's tangential infinitesimal spin
with the Hughes–Brezzi penalty `G*h/100`, as in the existing flat shell. It is an
explicit additional stabilization, not a claim of algebraic identity with the
original five-DOF element or with an Abaqus element. Dynamics adds the explicit drilling inertia described below.

Homogeneous and layered elastic linear statics are enabled.
Curved frequency and direct/modal harmonic analyses use the mass operator below.
Linear buckling is limited to coplanar plates with a flat preload. Curved-shell
buckling, non-flat preloads (including laminate extension–bending coupling), transient
and response-spectrum procedures are refused by name. Nonlinear steps do not support these shells.
The existing MITC3/MITC4 analysis paths are unchanged.

Reduced-integration quad8 shells are refused in the internal solver: Abaqus S8R is
not this full-integration MITC8. The deck writer also refuses full-integration quad8
shells rather than silently naming them S8R. A separately authored S8R benchmark deck
can compare the same physical problem and mesh, but cannot establish formula identity.
SC8R and SAX1 are not supplied by this implementation.

## Results and verification

Global bottom/middle/top stresses are sampled from the direct displacement field at
each of the eight nodes; they are neither Gauss extrapolations nor equilibrated
interlaminar stresses. Transverse shear is the constitutive shear stress, without the
flat kernel's parabolic 3/2 peak reconstruction. The existing single `axes` field is
the centre frame. Section resultants are thickness integrals projected onto that
fixed frame, **not** local tangent components at each node of a curved shell. No
result fields or binary format numbers are added.

`fem/tests/mitc8.rs` checks curved/offset rigid motions, constant axial extension of a
parabolic cylinder, a distorted straight-sided membrane patch, analytical pure-bending
and shear energies down to h/length=1e-4, consistent body force and offset moment,
folded geometry refusal, symmetry, and positive definiteness after removing six rigid
unknowns. Polynomial energies are computed independently of the element operators.
The command-level oracle cases `mitc8-pure-bending`, `mitc8-curved-extension` and
`scordelis-lo-roof-mitc8` provide actual static solutions with fixed tolerances.
The curved axial patch is not a bending-dominated validation; the roof is checked
separately. These checks do not validate every shell geometry or the unsupported
procedures. The axisymmetric LE7 thin-shell discrepancy remains unresolved.

## Layered elastic sections

The law follows Nettles, NASA RP-1351 (1994), sections II–IV, using the existing
`Lamina` energy-conjugate engineering-strain transformations. `Mitc8::thickness_rule`
splits at every ply interface; `material` evaluates the rotated plane-stress and
transverse-shear law in each ply. Thus the actual layer positions produce A, B and D
coupling, rather than replacing the stack by an averaged modulus. Constant layer
thicknesses and the fixed 5/6 shear correction are used.

At every surface point, `ply_angle` projects the supplied rectangular or cylindrical
material reference axis onto the interpolated-director plane, then applies the ply
angle. This physical fibre direction is held along the director through thickness;
it does not follow the varying covariant tangent at each depth. A singular reference
projection is rejected. With no section orientation, the reference is the local r
direction. Material axes and ply angles are independent of the shared global DOFs.

Existing result fields retain separate one-sided ply interface stresses and sampling
points. Ply stress and strain components are in material axes; each nodal `angle`
is measured in that node's reference-surface frame. For a curved patch, the single
reported centre `axes` does not reconstruct every nodal material basis. Use the model
geometry, shared directors and section orientation to reconstruct those bases. No
new result type or field is added. The fields are an equivalent-single-layer recovery,
not equilibrated interlaminar stresses or delamination predictions. Static gravity
integrates the actual density of each ply, including its first moment, as does
consistent dynamic mass.

Independent tests check isotropic single-ply degeneration on flat and curved
surfaces, nonsymmetric layer membrane/bending energy, one-sided interface stress,
explicit fibre directions under local renumbering, and reversal of normals/layer
order/offset. Command cases `mitc8-single-ply-bending`, `mitc8-navier-thin` and
`mitc8-navier-thick-quarter16` and `mitc8-coupled-strip` check the composite solve path. Navier references apply to the
symmetric cross-ply plate with B=D16=D26=0; they do not validate a general unsymmetric
curved laminate. Thick-plate full 4/8/16 meshes missed the fixed 3% band
(14.2/8.72/3.43% errors); the equivalent full 32 mesh passed at 0.999%.
Those coarse failures remain explicit in the thick-case README.


## Inertia and mass choices

`Mitc8::mass` inserts the Lee–Bathe (2004), equations (1)–(8), displacement
kinematics into kinetic energy; it does not interpolate inertia with the MITC strain
operator:

```
v = sum N_i (v_i + z omega_i cross n_i)
T_physical = 1/2 integral rho |v|² dV
T_drilling = 1/2 integral rho z² (sum N_i omega_i dot n_i)² dV
```

The second term is artificial drilling inertia, following the existing flat-shell
six-DOF convention. It is not physical inertia of a director spinning about itself,
not a massless drilling constraint, and not an assertion that S8R uses this model.
It makes the six-DOF mass positive definite. Rigid-body rotation energy tests include
this explicitly reported addition. Modes dominated by drilling are numerical modes;
no physical validation of their high-frequency spectrum is claimed.

Three Gauss points in each surface direction and per ply integrate the actual
thickness Jacobian and each layer density. Consistent mass retains inter-node blocks,
rotary inertia and the translation–rotation first moment (including offset). Missing
or nonpositive ply density is refused, not replaced by an average. Unit translation
inertia equals the consistent gravity vector.

`Mitc8::lumped_mass` uses positive diagonal scaling after Hinton, Rock & Zienkiewicz,
Earthquake Engineering & Structural Dynamics 4 (1976), 245–249,
[doi:10.1002/eqe.4290040305](https://doi.org/10.1002/eqe.4290040305).
The diagonal weights integral rho N_i² and integral rho z² N_i² are separately
scaled to integral rho and integral rho z². Each supplies three global components;
the rotational diagonal includes the same artificial drilling convention. Lumping
preserves total translation mass and these scalar rotary moments in the reference
thickness coordinate z. On curved elements this is not the rotary energy of the
interpolated, non-unit interior director. It omits I1,
inter-node coupling and exact rigid-rotation energy. It is not row-sum lumping:
quad8 corner row sums can be negative. Do not equate consistent and lumped spectra.

Both choices enter the existing assembly, Lanczos and direct/modal harmonic paths.
The Navier tests hold drilling and in-plane DOFs to isolate the symmetric laminate's
flexural branch. They compare the first three frequencies and complex phasors below
resonance against the independent separated 3x3 FSDT equations, for both mass choices.
The curved free-body command test checks six rigid modes, and shared-node tests check
assembled mass under cyclic/reversed numbering. General curved spectra still need
problem-specific convergence checks. Per-ply damping and unrepresented analysis
combinations retain the existing named refusals.

## Linear plate buckling

`Mitc8::plate_geometric_stiffness` uses the existing small-slope flat-plate theory
(Bathe, Finite Element Procedures, 2014, section6.8.2):
`K_sigma = integral grad(w)^T N grad(w) dA`.
At each point, N is obtained by integrating the actual MITC preload strains and each
ply law through thickness; compression is negative. `buckle::initial_stress_stiffness`
assembles this operator from the solved static preload, then the existing Lanczos
buckling path solves `(K + lambda K_sigma) phi = 0`.

Every shell reference surface/director must be coplanar. The preload must remain flat,
including director tilt; otherwise the solver refuses the combination. This supports
flat symmetric laminates, or other layers whose actual constrained preload remains
flat; it does not silently drop pre-bending from an unsymmetric laminate. No curved
shell geometric tangent, finite rotations, nonlinear material tangent, follower
pressure or load stiffness is provided. The same scope checks protect prestressed
linear perturbations. The Bryan square-plate command case checks the entire preload,
material/geometric assembly and eigenproblem against `Ncr=4*pi²*D/a²`.

Frequency and buckling retain the existing frequency-result representation; harmonic
uses the existing mass-aware harmonic result. Native save/reopen tests compare modal,
buckling and complex displacement probes exactly, with no result format extension.
