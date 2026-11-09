# Six-node curved MITC6-a

`tri6` with a full-integration homogeneous or composite shell section selects MITC6-a.
It has exactly six geometric/displacement nodes: corners `(0,0),(1,0),(0,1)`, followed
by edge midpoints `0–1,1–2,2–0`. Each carries global `[ux,uy,uz,rx,ry,rz]`; there are
36 element DOFs. The centroid supplies two shear samples, never an extra node or
condensed unknown. This is not MITC6-b, MITC7, MITC8 or MITC9.

The source is P.-S. Lee and K.-J. Bathe, *Development of MITC isotropic triangular
shell finite elements*, Computers & Structures 82 (2004), 945–962,
[doi:10.1016/j.compstruc.2004.02.004](https://doi.org/10.1016/j.compstruc.2004.02.004).
The word isotropic describes invariance under local node numbering; it does not
restrict the elastic material to an isotropic law.

| Public formula | Implementation in `crates/fem/src/mitc8.rs` |
|---|---|
| (1)–(8), quadratic geometry and displacement, covariant strains | `QuadraticShell<6,36>::map`, `displacement_gradient`, `covariant` |
| (29)–(34), Fig.6, linear in-plane covariant strain components | `triangle_layer`; directional qq strain reconstructs rs |
| (35)–(38), Fig.7a, edge-linear / interior-quadratic shear | `triangle_shear`; six edge points plus two centroid components |
| Plane stress constitutive contraction | `strain_operator`, `material`, `stiffness` |

Membrane and bending coexist through `z = -offset*h + h*t/2`. Shape functions
interpolate the unit nodal directors without renormalizing their interpolation;
rotational displacement is `z*(theta cross director)`. This preserves infinitesimal
rigid motion. The actual signed thickness Jacobian is used in stiffness, mass and
body loads; inverted or singular sampled maps are errors. Checks are not a global
positivity certificate for arbitrary curved geometry. Thickness is constant per
section; rotations and strains are small.

Three thickness Gauss points per ply resolve the material law, directions and
membrane/bending coupling. Surface integration averages all six barycentric
permutations of the positive 3x3 conical product (degree four), preserving numbering
invariance on curved maps. This quadrature choice is separate from the paper's
strain interpolation. Shared directors and global rotations use the same smooth
patch interface as [MITC8](mitc8.md); nonconforming midsides, nonmanifold edges and
sharp folded junctions are refused. Tri6/quad8 smooth interfaces can share directors;
linear/quadratic shell junctions remain refused.

The additional drilling penalty and drilling inertia are the documented MITC8
six-DOF extensions, not part of the paper's physical five-DOF kinematics. Consistent
mass retains actual ply density, translation/rotation coupling and rotary inertia.
Positive diagonal scaling preserves total translational mass and scalar thickness
moment; it removes coupling and is not row-sum lumping or exact rigid-rotation
energy. Ply stress output is direct nodal strain recovery in material coordinates,
not integration-point extrapolation or equilibrated interlaminar stress.

`solve.static`, `solve.frequency`, and direct/modal `solve.harmonic` use the production
assembly. `solve.buckle` uses the solved reference membrane stresses, but only for
coplanar plates with flat preloads. Curved-shell buckling, prebent configurations,
finite-rotation/nonlinear shells, transient and response-spectrum MITC6-a analyses
are not validated and remain refused. Native result save/open uses the existing
static, frequency/buckling and harmonic result layouts.

Verification includes independent rigid motion, membrane, bending and shear energy
integrals; thickness ratios down to 0.0002 on a unit triangle; positive constrained
stiffness and mass; shared curved edges and cyclic/reversed numbering; one-ply
reduction and laminate coupling. Production cases use pure bending, a layered strip,
Navier plate vibration/response, Bryan buckling and the Scordelis–Lo roof. The
32-triangle coarse harmonic mesh misses the fixed 3% band (up to 17.4%) and the
coarse Bryan mesh misses 2% (-3.54%). These failures remain applicability limits;
they do not qualify as validated coarse meshes. The 128-triangle harmonic mesh
still misses 3% (up to 4.81%); 288 triangles reduce all four solution/mass variants
to at most 2.24%. Bryan buckling with 128 triangles has -1.08% error; the
128-triangle curved roof has -0.053% displacement error. These are separate
accepted meshes, not replacements for the failed coarse results. The source paper itself reports
residual locking for some thin curved-shell problems: no claim of uniformly
locking-free performance is made.

Abaqus export refuses an automatic MITC6-a mapping. `STRI65` is a six-node thin shell
with five active DOFs at most nodes and a different transverse-shear treatment;
`S6` is not a supported mapping. Neither name establishes formula equivalence.
See the public [shell selection manual](https://docs.software.vt.edu/abaqusv2026/English/SIMACAEELMRefMap/simaelm-c-shellelem.htm).
A separately authored thin-shell comparison input must state these semantic limits;
it is not a lossless native export. No general thick-shell STRI65 comparison is claimed.
