# Axisymmetric meridional shell

`formulation.axisymmetric_shell: true` on a shell or composite-shell section
selects a two-ring, linear conical shell on `beam2` connectivity. It does not
select a beam constitutive law. Coordinates are `(r,z,0)` in model XYZ; active
DOFs per node are `(Ux,Uy,Rz) = (u_r,u_z,beta)`. Each node represents an entire
ring. No circumferential displacement, torsion or non-axisymmetric harmonic is
represented. Ring loads and reactions are full-circumference integrals in N and
N m. A radial ring resultant is not a net Cartesian vector force; the static
balance diagnostic checks the axial component only.

## Formula and code correspondence

The source is Simões Moita et al., *Applied Sciences* 13 (2023), 5325,
[doi:10.3390/app13095325](https://doi.org/10.3390/app13095325), equations 5–6
(kinematics) and 8–11 (laminate resultants). See `fem::axishell::Axishell`:
`operator` implements the strains, `stiffness` integrates their energy, `mass`
integrates the velocity field, and `stresses` recovers section results.

For meridional arclength s, tangent t=(c,s_z), normal n=(-s_z,c), linear nodal
shape functions N and beta positive about model Z:

- eps_s = c u_r' + s_z u_z'; eps_h = u_r/r;
- k_s = -beta'; k_h = -c beta/r;
- gamma_sn = -s_z u_r' + c u_z' - beta.

Rows in the existing laminate operator are `(eps_s,eps_h,0,k_s,k_h,0,gamma,0)`;
columns are node-major `(u_r,u_z,beta)`. The radius is interpolated, ds/dxi=L/2,
and integration measure is 2 pi r ds. Membrane and bending use two Gauss points;
shear uses its midpoint value and one-point integration. This **selective** rule
differs from the fully reduced rule in the cited paper. It is also not a claim
that this discretization is Abaqus SAX1. No MITC or solid B-bar operator is used.
The thin-shell measure omits through-thickness curvature Jacobian corrections;
nonnegative radii and the hoop-curvature thickness screen below are required.

Section material axes are `(t,negative model Z,n)`, a right-handed triad. Positive
thickness follows n: inward for a meridian running upward on a cylinder, upward
in model Y for an annulus running radially outward. Bottom-to-top plies follow
this normal; reversing connectivity also reverses it, so an unsymmetric physical
layup must reverse its layer stack and offset. Shared beta is a global rotation,
not an independently chosen element-local director angle. Only 0/90-degree plies
(modulo pi) are accepted. Off-axis extension/twist coupling requires a missing
circumferential DOF and is rejected. Existing laminate ABD, transverse shear and
bottom/middle/top or ply-point stress routines are reused. Reported global stress
components use the model meridian frame `(radial,axial,hoop)` at the represented
meridian, not one spatially uniform Cartesian tensor around the entire ring.

For density moments I_j=integral rho zeta^j dzeta, the sectional mass block is
`[[I0,0,-c I1],[0,I0,-s_z I1],[-c I1,-s_z I1,I2]]`.
`consistent` retains this translation/rotation coupling and two-point spatial
integration. `lumped` uses separate translational and rotational HRZ diagonal
scaling, preserving integrated I0 and I2 but dropping the I1 coupling; it is not
a row sum of a coupled mass matrix. Existing Rayleigh damping uses the selected
assembled M and K. The rigid-motion mass and total mass tests include ring measure.

## Analysis scope and refusals

Production commands `solve.static`, `solve.frequency`, and direct/modal
`solve.harmonic` assemble this stiffness and mass through the common constrained
solver. Nodal ring forces, moments, density-integrated gravity and body force
are supported. Pressure and thermal loads are explicitly rejected; callers must
supply independently integrated mechanical ring forces. Point masses/springs,
mixed 3D sections, negative-radius axis crossing, explicit material frames, transient/response
spectrum, explicit dynamics and sensitivity analyses are not supported.

`solve.buckle` supports **flat annular plates with a flat membrane preload only**.
The axisymmetric von Karman term is delta eps_s = (w')^2/2, giving
KG = integral 2 pi r N_s (dw/dr)^2 dr. `plate_geometric_stiffness` implements it.
The existing generalized eigensolver supplies load factors. Cylindrical, conical
and curved-shell buckling, prebent states and non-axisymmetric modes are named
refusals. A membrane-preloaded simply supported annulus (r=1..2, nu=.3) compares
against the Bessel solution w=A+B log r+C J0(kr)+D Y0(kr), w=0 and
w''+nu w'/r=0 at both edges, k^2=P/D. The first critical load is
93.80004147396585 for E=1e8, h=.01; 32 elements give 93.91324277507324
(+0.12068%). A huge negative factor warning from the singular geometric-stiffness
pencil is retained, not reinterpreted as a physical reverse-load mode.

Abaqus SAX1 export is explicitly refused: active DOF correspondence (1,2,6)
alone does not establish equivalence of shear interpolation, thin-shell measure,
normal convention and section semantics. Nastran and LS-DYNA export likewise
refuse this section rather than silently emit a beam. Existing native 07 result
kinds suffice; no new identifier or archive version is introduced.

## Evidence

`tests/oracle/axishell-cylinder` checks the exact radial membrane solution.
The engine test `axishell` additionally checks its 50.329212104487 Hz breathing
mode, 62.831853071796 kg total mass, and 30 Hz direct/modal complex response for
both mass choices against the scalar oscillator, then saves and reopens results.
`axishell-annular-laminate` checks an unsymmetric two-layer spherical-bending
patch and membrane/bending coupling. `axishell-annular-buckling` checks the
independent Bessel root. FEM tests check axial rigid motion, conical membrane
strain, reversed connectivity, curvature, thin limits and mass. These are
restricted-scope checks, not evidence for arbitrary curved laminates or poles.

The public [Abaqus axisymmetric-shell theory](https://docs.software.vt.edu/abaqusv2025/English/SIMACAETHERefMap/simathe-c-axishells.htm)
describes SAX1 with lumped mass and one-point distributed-load integration,
and a finite-strain theory. Those choices do not match the present consistent
mass option or small-strain selective rule. The public
[shell overview](https://docs.software.vt.edu/abaqusv2026/English/SIMACAEELMRefMap/simaelm-c-shelloverview.htm)
defines its axisymmetric positive normal by a counterclockwise quarter turn
from node 1 toward node 2, and excludes explicit orientations for SAX1/SAX2.
These conventions are checked individually; matching a normal convention does
not establish an element-formula mapping.

## Regular poles and LE7

A meridian may end at r=0; negative radius and an edge entirely on the axis are
refused. NASA CR-2568 (1975), section 4.1.11, supplies the regularity conditions
and endpoint limiting procedure, not this element's interpolation. The common
production constraint path imposes u_r=beta=0 at every pole, leaves axial motion
free, and rejects conflicting prescribed values. For smooth fields,
`lim u_r/r = u_r'/r'`, `lim beta/r = beta'/r'`; hence
`lim k_h = -beta'`. The constrained pole shape contributes zero to the hoop
operator. The non-pole linear shape uses `N'/r'` at the endpoint. Interior Gauss
points retain their positive radii and the measure 2*pi*r; no epsilon radius,
new cap node, mass, or fictitious stiffness is introduced. Stress recovery also
rejects nonregular pole coordinates. Faceted tangent/director errors near a
smooth pole converge with meridian refinement; a sharp conical apex is not a
smooth-pole accuracy claim.

The thin-shell screen now uses `h*abs(dz/ds)/r_mid <= .05` (hoop curvature),
not h divided by the smallest nodal ring radius. A sphere has finite curvature
at its pole although its ring radius vanishes. Meridional curvature is resolved
by conical facets; this screen alone is not a universal thin-shell certificate.

`axishell_pole` builds a distinct, SI-unit NAFEMS LE7 **shell** model: R=1 m,
h=.025 m, E=210 GPa, nu=.3; full-circumference loads integrate 1 MPa normal
pressure on the mid-surface. F fixes axial motion and symmetry rotation; A has
the automatic pole constraints. D is a node at (1,1.4034,0), with its physical
outer surface at (1.0125,1.4034,0). Connectivity runs from F toward A, making
positive thickness inward: **bottom**, not the default nodal top stress, is the
outer surface. Both adjacent element endpoint sigma_YY values are checked
against 25.86 MPa at the unchanged 1% limit; 50/100/200 elements per metre give
outer endpoint pairs (25.61413,25.93402), (25.71127,25.88083),
(25.76013,25.84542) MPa. Native static results reopen exactly. This does not
replace the old axisymmetric-solid LE7 recipe, its millimetre-scale geometry,
or its -2.233% shell-reference failure.
