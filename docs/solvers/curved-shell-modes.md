# Curved quadratic-shell modal perturbations

MITC8 can use an earlier **linear static** `base_state` for
`solve.frequency`. The reference-geometry tangent is
`K_T = K_elastic + K_sigma + K_pressure`, with reference mass. The stress
operator is the same Green-strain/director Hessian used in curved-shell
buckling, evaluated from the actual preload displacements at the original
layer quadrature. Nodal averaged output stress is not a restart state.

Set the frequency step's `perturbation_follower_pressure: true` to add
`K_pressure = -d f_pressure / d u` for pressure face loads in the base step.
The default is false and retains dead loads. This model option is distinct
from the buckling procedure's existing `follower_pressure` option. A nodal
force is always a dead load; assigning the same initial resultant does not
make its incremental loading law a follower pressure.

The pressure operator reuses Schweizerhof and Ramm (1984), *Displacement
dependent pressure loads in nonlinear finite element analyses*, Computers &
Structures 18:1099–1114, doi:10.1016/0045-7949(84)90154-8, through
`fem::face::follower_pressure` and `solver::buckle::pressure_stiffness`.
`solver::prestress` checks the **reduced** pressure tangent for symmetry
before removing roundoff asymmetry. Nonconservative pressure boundaries
are refused, not solved using a symmetrized substitute.

## Small-preload cylinder check

The reproducible command test `curved_shell_prestress` uses a quarter cylinder
with R=1 m, h=.02 m, L=.2 m, E=1e7 Pa, nu=0, rho=1000 kg/m³. All axial
translations are held, both ends remain radially free, and the meridian
planes have mirror supports. It is the axially uniform plane-strain ring
limit, **not** a simply supported finite cylinder or a capped pressure vessel.
There is no end-cap axial force. For nu=0 the membrane resultants are
N_theta=pR and N_z=0; static expansion is u_r=pR²/(Eh). Positive internal p
uses negative native pressure on the outward reference surface.

NASA MCR-68-87, vol. II (1969), p. II-1,
[public report](https://ntrs.nasa.gov/api/citations/19690031776/downloads/19690031776.pdf),
gives the classical inextensional ring/cylinder frequency. The following
pressure extension is derived independently from its displacement ansatz
and the public virtual-work/load-potential formulation, not transcribed from
an inaccessible pressure-frequency paper. With radial w=A cos(n theta),
tangential v=−A sin(n theta)/n and bending rigidity D=Eh³/[12(1−nu²)], the
integrated mass coefficient per axial length is rho*h*R*pi*(1+1/n²).
The three stiffness coefficients are:

- elastic bending: pi*D*(n²−1)²/R³;
- hoop initial stress: pi*p*(n−1/n)²;
- internal follower-pressure potential: pi*p*(1−1/n²).

The last follows by differentiating the enclosed area:
`delta² area = integral(w²+v²+2*w*v') dtheta` and potential `−p*area`.
Therefore

`omega_n² = [D*n²*(n²−1)²/R⁴ + p*n²*(n²−1)/R] / [rho*h*(n²+1)]`.

For dead radial loads omit the last energy: replace the numerator's pressure
factor `n²*(n²−1)` by `(n²−1)²`. At n=2 these slopes are 12/5 and 9/5,
respectively. Continuing to negative p gives follower/dead critical external
pressures 3D/R³ and 4D/R³, agreeing with the separately checked buckling
limits. This is not a fluid-coupled mode: fluid added mass, compressibility
and acoustic dynamics are absent.

Four- and eight-element-per-quadrant MITC6-a meshes passed the frequency
gate but failed the independent 1% pointwise static expansion check (about
4.82% and 3.81% at the first checked node under 5 Pa). The eight-element
field has opposite end-corner deviations of about +3.81% and -3.85%, while
the mid-length nodes are within 0.18%. These coarse-preload failures are
retained; the sixteen-element first-node error was still about 2.26%. They
are not accepted by loosening a gate or fitting pressure. No altered director or strain interpolation is used to fit this benchmark.

The command test compares zero-pressure degeneration, p=5 and 10 Pa and
the dead-load alternative, checks the actual reduced-pencil residuals,
total mass and native save/reopen. The analytic reference omits rotary
inertia and transverse shear; the FE mass retains both translational and
rotary inertia. At h/R=.02 the a-priori frequency gate is 3% for the coarse
four MITC8 quadratic elements per quadrant, with a separate 1e−6 algebraic residual
gate. No material, density, pressure or tolerance is fitted to a result.

MITC6-a curved prestressed modes and its follower-pressure modal tangent are
**refused by name**. At 64 circumferential elements per quadrant the two-mode
run failed the second-mode actual residual (1.04e-5 at zero pressure); asking
only for the first mode also failed (2.43e-5). This is not evidence that those
modes are zero-energy mechanisms: the actual residual is inadequate and the
fine-mesh preload check was not reached. The coarse preload discrepancies and
fine residual failures remain open. Ordinary MITC6-a modes and its separately
verified static/buckling paths retain their previous scope. No shared eigen
solver tolerance or frequency result identity was changed to accept this case.

## Scope and persistence

Only MITC8 surface shells with unchanged direct supports are admitted
for this curved modal path. Mixed solids/beams/ring shells, pretension cuts,
MPCs/inclined supports, initial stress/temperature fields, thermal preload,
offset follower faces and nonconservative pressure
are refused. Curved prestressed harmonic response remains unsupported.
The pressure-frequency example verifies a homogeneous isotropic shell with
consistent mass; it does not validate pressure-dependent laminate spectra or
preloaded lumped-mass accuracy.
Finite-preload deformation, finite rotations, material-tangent updates,
nonlinear shell equilibrium and postbuckling are outside the formulation.
The linear preload and reference mass are a small-strain approximation;
a frequency increase is not evidence for accuracy at large pressure.

The model saves the pressure-tangent option and equilibrium dependency.
Numeric modes retain the existing result representation. This change adds
no result kind and does not edit frequency mass-provenance or codec fields;
retain the solve command's explicit mass selection and section policy in
external evidence when using a writer whose frequency payload lacks that
identity. Faithful S8R/STRI65/SAX1/SC8R mappings are not established.
