# Affine eight-node continuum solid-shell

Use a `solid` section with `formulation: {"solid_shell": true}` on a swept `hex8`
mesh. The first four nodes form the lower face, the next four the upper face;
1–5, 2–6, 3–7 and 4–8 are the positive thickness pairs. Each node has only global
`ux, uy, uz`: 24 external coordinates. Thickness, eccentricity and density volume
come from the actual geometry, not a midsurface thickness property. There is no
independent rotational or drilling inertia.

This is the small-strain, orthogonal affine specialization of Harnau,
Schweizerhof and Hauptmann, *On solid-shell elements with linear and quadratic
shape functions for small and large deformations*, ECCOMAS 2000,
[original paper](https://congress2.cimne.com/eccomas/proceedings/eccomas2000/pdf/826.pdf).
`crates/fem/src/solid_shell.rs` maps equations (1),(3) to paired trilinear
geometry/displacement, (5) to transverse shear ANS, (16) to thickness ANS, and
(6),(15),(10) to six membrane and four thickness EAS parameters. In local
engineering-strain order `[xx, yy, zz, xy, xz, yz]`, the membrane modes are
`xx:r, yy:s, xy:r, xy:s, xx:rs, yy:rs`; the thickness modes are
`zz:t, zz:rt, zz:st, zz:rst`. Constant Jacobian and basis scale factors are
absorbed into the local parameters on the admitted affine geometry.

Eight Gauss points integrate the full 3D elastic law. The ten EAS parameters
are statically condensed; their recovered values also enter the reported strains
and stresses. Stiffness integrates the relaxed strain energy rather than subtracting
two large stiffness blocks. There is no plane-stress constitutive reduction,
5/6 shear factor, reduced quadrature or hourglass penalty. The development did
not use the different Ko–Bathe 2018 3D-MITC8 strain interpolation.

## Current limits

Geometry must be an orthogonal affine brick, possibly rigidly rotated. Curved,
warped, oblique, variable-thickness and inverted sweeps are refused. Maximum
edge aspect ratio is 100. More extreme ratios exposed loss of bending energy
accuracy when storing the translational stiffness in f64; they are refused,
not presented as passing thin-shell verification. Refine the in-plane mesh.
This limit is per element, not the whole structure's length/thickness ratio.

Linear static, consistent/lumped mass, modal and direct/modal harmonic assembly
use the existing production solvers. Mass is the ordinary 3D translational
integral; rotational kinetic energy follows from differential upper/lower motion.
Lumped mass uses positive diagonal scaling. Static condensation introduces no
additional physical inertia. Buckling/prestressed perturbation, nonlinear,
transient/response-spectrum, thermal and initial-stress loads, and spatially
varying cylindrical material axes are explicitly refused. Existing curved
MITC8/MITC6 buckling refusals are unchanged.

A laminate is currently represented by **separate physical swept layers** with
shared interface nodes and a solid section/material orientation for each layer.
Rectangular orthotropic axes are supported. This retains thickness-normal
stress and uses the actual density of each layer. It is not single-element
through-thickness composite-section integration: `composite_shell` on a hex8 is
unsupported. Existing surface-shell layups are unchanged.

No faithful SC8R export is claimed: SC8R uses a different reduced-integration
continuum-shell formulation. Native Abaqus, Nastran and LS-DYNA writers refuse
the new formulation by name rather than silently writing an ordinary solid.
An independently prepared same-mesh comparison deck must state its formulation
and thickness/material semantics explicitly.

## Verification

The independent element tests cover rigid zero strain, all constant strains,
constant-curvature plate energy including Poisson thickness relaxation, uniform
shear energy, positive constrained stiffness, and hostile/unsupported geometry.
Command recipes `solid-shell-pure-bending`, `solid-shell-thick-bending` and
`solid-shell-coupled-layers` compare with closed-form bending and extension–bending
coupling. `solid-shell-ko-cantilever` uses the regular mesh dimensions and
reference in Ko & Bathe (2018), doi:10.1016/j.compstruc.2018.02.015, Figure 6(a),
Table 1; it does not assert that the two element formulations are identical.
Native command tests additionally check total mass, the slender cantilever's
Euler frequency, the low-frequency static limit and saved/reopened results.
