# Eight-node curved MITC8 statics

`quad8` with an elastic `shell` or `composite_shell` section and **full** integration uses
`fem::mitc8::Mitc8` in `solve.static`. This is an eight-node serendipity shell, not MITC9:
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
original five-DOF element or with an Abaqus element. No artificial mass is introduced.

Homogeneous and layered elastic linear statics are enabled.
Modal/frequency/harmonic/transient analyses, initial-stress stiffness and buckling are
not implemented for MITC8 and are refused. Nonlinear steps do not support these shells.
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
integrates the actual density of each ply, including its first moment; a dynamic
mass matrix is still not supplied.

Independent tests check isotropic single-ply degeneration on flat and curved
surfaces, nonsymmetric layer membrane/bending energy, one-sided interface stress,
explicit fibre directions under local renumbering, and reversal of normals/layer
order/offset. Command cases `mitc8-single-ply-bending`, `mitc8-navier-thin` and
`mitc8-navier-thick-quarter16` and `mitc8-coupled-strip` check the composite solve path. Navier references apply to the
symmetric cross-ply plate with B=D16=D26=0; they do not validate a general unsymmetric
curved laminate. Thick-plate full 4/8/16 meshes missed the fixed 3% band
(14.2/8.72/3.43% errors); the equivalent full 32 mesh passed at 0.999%.
Those coarse failures remain explicit in the thick-case README.
