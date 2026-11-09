# Thirteen-variable enhanced-gradient Hex8

The native `enhanced_strain13` section option is independent of the existing
nine-mode formulation. It does not authorize a C3D8I/IH import/export mapping.

Let J map parent coordinates to physical coordinates, J0 its value at the centre,
and r = det(J0)/det(J). For each parent direction k and physical component c,
use the symmetric gradient of r e_c tensor (J0^-T e_k xi_k). These nine fields
are supplemented by r times xy I, xz I, yz I and xyz I in parent coordinates.
Engineering shear components form a 6 by 13 matrix A. Integrate B^T D B,
L = B^T D A and H = A^T D A with full 2 by 2 by 2 quadrature. The nodal stiffness
is Kuu - L H^-1 L^T; internal amplitudes solve H a = -L^T u. Stress recovery uses
D(Bu + A a). Loads act on nodal displacements, not internal strain variables.
Cholesky solves replace explicit inverses. The determinant ratio makes the
integrated enhanced field zero on valid distorted cells. Patch consistency does
not guarantee arbitrary-distortion bending accuracy.

Sources: [Abaqus public Theory Guide, linear incompatible formulation](https://docs.software.vt.edu/abaqusv2025/English/SIMACAETHERefMap/simathe-c-incompatible.htm)
provides the explicit three-dimensional basis. [Simo and Rifai (1990),
IJNME 29, 1595–1638](https://doi.org/10.1002/nme.1620290802), sections 2.3 and 3.1,
provides the variational elimination and patch-test requirements; its broader
nonlinear theory is not implemented by this option. Historical nine-mode work
by [Taylor, Beresford and Wilson (1976)](https://doi.org/10.1002/nme.1620100602)
does not establish equality of nine- and thirteen-dimensional enhancement spaces.

Scope: linear isotropic, small-strain static, full Hex8 with positive Jacobians
at the centre and integration points. Other procedures, material axes,
thermal/initial stress, hybrid pressure and other enhancement combinations are
refused. Fully incompressible material is outside this displacement formulation.
Native save/open preserves the option and existing static results. No result
kind or binary version is added. Verification covers rigid modes, constant
strain on affine/non-affine cells, zero-mean enhancement, absence of extra zero
energy modes, rotation covariance, continuum pure bending and integrated energy
versus condensed nodal energy. Independent distorted-case equivalence remains
an explicit validation gate, separate from these theory checks.

The nonhybrid thirteen-variable basis is explicit. The hybrid case is not inferred
from it: the public theory discussion omits the extra dilatational modes for
hybrids, while the element-library table lists thirteen incompatible variables
for C3D8IH and separately describes four linear-pressure variables. This does not
supply a verified complete mixed block for this implementation. C3D8IH therefore
remains outside scope even if the nonhybrid candidate later passes comparisons.
Stress here is the constitutive value at each Gauss point, followed by the
existing trilinear nodal extrapolation; it is not a claim to implement every
assumed-stress projection discussed in the paper.
