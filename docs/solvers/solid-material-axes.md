# Material directions in linear solids

`Solid::stiffness_in_material_axes` and `Solid::stresses_in_material_axes` use the
section's material directions to rotate the engineering-constant elasticity into
global coordinates. A rectangular frame is constant. A cylindrical frame is radial,
circumferential, axial, evaluated at **each stiffness integration point** using its
physical position. An undefined radial direction on the axis is an error.

The stiffness is the usual displacement integral `K = ∫ Bᵀ D(x) B dV`; direction
variation changes `D(x)`, not the element geometry or its quadrature rule. This is
the formulation of Hughes, *The Finite Element Method* (2000), chapter 3, with the
tensor transformation of Jones, *Mechanics of Composite Materials*, second edition
(1999), §2.4. The public [Abaqus orientation reference](https://docs.software.vt.edu/abaqusv2026/English/SIMACAEKEYRefMap/simakey-r-orientation.htm)
also defines material calculations at integration points; no product was run and
no input deck, code or reference output was copied.

Both kernels return their matrices and fields in **global axes**. Stress and strain
components remain `xx, yy, zz, xy, xz, yz`, with engineering shear strains. Optional
initial stresses are global nodal stresses, interpolated with the shape functions:
`σ = D(x) ε + σ₀(x)`. These solid fields are distinct from the laminate ply fields,
which are reported in each ply's material axes.

Rectangular directions use the existing constant-elasticity formulations, including
incompatible modes. A hybrid formulation still requires isotropic elasticity.
Cylindrical directions currently support the full and reduced displacement
formulations; incompatible-mode and hybrid formulations give named errors. These
methods do not implement thermal strains, nonlinear material evolution or rotating
material axes at finite strain.

`crates/fem/tests/solid_material_axes.rs` checks constant-frame equivalence, six rigid
motions, exact matrix symmetry and invariance under a general rotation and
translation. Its affine brick independently evaluates the cylindrical constitutive
law at the Gauss points, compares both stress and strain energy, and ensures a
centroid-only frame would fail. Invalid axes, nonfinite fields and unsupported
formulations are tested as errors. The affine check establishes quadrature
consistency, not exact integration of an arbitrary cylindrical direction field.
