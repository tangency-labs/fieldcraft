# Geometric nonlinear transient

`solve.transient.nonlinear {"step":"Motion","mass":"consistent",
"max_iterations":30,"tolerance":1e-8}` integrates an existing transient step with
full Newton iterations. It is an experimental **rod2** path: finite rotations and
St Venant–Kirchhoff axial elasticity, homogeneous fixed supports, dead nodal
forces (including their amplitude curves), consistent or lumped reference mass,
optional point masses and mass-proportional physical viscous damping.

The reply contains SI nodal displacement, velocity, acceleration and T, V, W,
physical damping dissipation, energy balance, Newton correction counts and
residuals. Successful solves store a distinct **MCRSLT07 kind 11** result,
which `file.save_as` and `file.open` preserve without solving. Use
`solve.probe {"step":"Motion","nodes":[2]}` to retrieve the complete nonlinear
SI history (all nodes and elements, including reference identities). The node
selector is checked against the model but does not reduce this full-history reply.
This is not a linear kind-6 response. A failed solve restores the previous session
result; no rejected trial state is stored as an accepted frame. The error reports
the failed increment and last committed time. All accepted increments, including
those omitted by `output_every`, retain correction counts and relative residuals.

The model's existing transient step supplies time controls; the explicit nonlinear
command selects the rod formulation and stores its actual settings with the result.
Saving a result-free model or exporting via `deck.write` does not select nonlinear
analysis. Reproduce with `solve.transient.nonlinear`; use `deck.write.nonlinear`
for the separately supported user-element export. No archive-version upgrade is
required; older unified-07 readers reject unknown kind 11; pre-07 readers retain their
unsupported-magic result-drop warning. Kinds 0–8 are unchanged; 9/10
remain unused, and old 05/06 kinds above 5 remain refused.
See [the complete layout](../../results-kind-layout.md#nonlinear-rod-kind-11).

Solids, shells, beams, rod3, contacts, couplings, foundations, springs/dashpots,
plasticity/hyperelasticity, material damping, follower loads, base motion,
nonzero prescribed displacements, requested reaction/thermal field output, modal reduction, adaptive steps and
stiffness-proportional damping are refused. Point masses must provide a positive
mass operator on the free coordinates. This rod formulation is not a tension-only
cable law: compression is elastic, and collapse is refused. The nonlinear static solver also supports pure elastic rod2 equilibria through
this shared force/tangent kernel; these equilibria can supply the frozen tangent
to frequency and direct/modal perturbation steps (see [base states](transient.md)).

## Equations and verified sources

Newmark uses u*=u_n+h v_n+h²(1/2−β)a_n and v*=v_n+h(1−γ)a_n,
then a=(u−u*)/(βh²), v=v*+γh a. Average acceleration uses β=1/4,
γ=1/2. HHT uses α∈[−1/3,0], β=(1−α)²/4 and γ=1/2−α.
For constant Cartesian M and C and displacement-independent applied forces,

```
r = M a + (1+α)[f_int(u)+C v−f_(n+1)]
          − α[f_int(u_n)+C v_n−f_n]
K_eff = M/(βh²) + (1+α)[K_t(u)+γ C/(βh)]
K_eff Δu = −r
```

Every iteration recomputes internal force and its consistent tangent, including
geometric stiffness. The implementation solves the equivalent βh²-scaled
system, using symmetric pivoting when necessary. A time step commits u/v/a only
after the relative free-force residual passes; a failed iteration leaves the
previous state intact and the command leaves the session unchanged. Initial
acceleration solves M a0=f(0)−f_int(u0)−C v0 without HHT weighting. Initial
translations and velocities must satisfy the same constraints as the solve.

For a reference rod length L, current chord d, and axial rigidity EA,
E_G=(d·d−L²)/(2L²), V=EA L E_G²/2, f=EA E_G d/L, and
K=EA/L(d⊗d/L²+E_G I), assembled with the usual opposite endpoint signs.
Thus finite rigid rotation has zero strain energy, and the tangent is the exact
force derivative. M is assembled on the reference geometry with the explicitly
selected mass convention; the same matrix defines initial acceleration and
T=vᵀMv/2. Here C=a_M M is fixed and nonnegative.

Work W uses endpoint-force trapezoidal work in displacement; physical
D uses trapezoidal integration of vᵀCv in time. The reported balance is
T+V−(T0+V0)+D−W. It includes time integration and quadrature error, and must not
be labelled physical damping or a pure measurement of HHT algorithmic loss.
HHT high-frequency dissipation is examined separately with C=0 and step
refinement. Linear unconditional stability claims do not guarantee nonlinear
convergence or conservation for an arbitrarily large step.

Sources actually checked:

- Bathe & Cimento (1980), §2 and §3.1, equations (2)–(4), (14)–(17),
  [author-hosted paper](https://web.mit.edu/kjb/www/Publications_Prior_to_1998/Some_Practical_Procedures_for_the_Solution_of_Nonlinear_Finite_Element_Equations.pdf):
  nonlinear dynamic balance, trapezoidal kinematics, current tangent and the
  effective mass contribution 4M/h². Full Newton is used here, not their BFGS
  or modified-Newton alternatives.
- Bathe's [MIT lecture 13](https://ocw.mit.edu/courses/res-2-002-finite-element-procedures-for-solids-and-structures-spring-2010/resources/lecture-13/)
  and [lecture 14 study guide](https://ocw.mit.edu/courses/res-2-002-finite-element-procedures-for-solids-and-structures-spring-2010/e55f2c75fd9a8f8632f85c662b184167_MITRES2_002S10_lec14_22.pdf),
  topic 14 pp. 14-30–14-32: large-displacement pendulum and full Newton at each
  step; iteration convergence alone does not prevent large-step instability.
  The original elastic pendulum here is not a reproduction of that guide's
  rigid-limit parameter set or its plotted numerical values.
- Hughes, *The Finite Element Method*, ch. 9, is the requested integration
  reference; the book's complete chapter has not been available for a new
  page-by-page check. HHT formulas were checked against the publicly available
  [Sierra Solid Mechanics 5.20 theory manual](https://www.sandia.gov/app/uploads/sites/315/2025/07/SM_Theory_5_20.pdf),
  equations (12.5)–(12.6), using its negative-alpha convention. The original
  Hilber–Hughes–Taylor (1977), [doi:10.1002/eqe.4290050306](https://doi.org/10.1002/eqe.4290050306),
  is the bibliographic origin; its complete paper has not been newly audited. The nonlinear effective tangent above is
  explicitly differentiated from that residual, not attributed to an unchecked
  book equation number. Bathe textbook edition-specific section numbering also
  remains distinct from the publicly checked lecture/paper equations.

Verification is the actual model/command path for an elastic pendulum and a
sagged two-segment tensioned string against independent Cartesian RK4 histories,
three time steps, both masses, with and without physical damping. Separate
force/energy finite differences check the finite-rotation tangent. A forced
Newton failure checks state restoration and preservation of stored linear
results. These checks do not validate unsupported solid/contact/material laws.

For the original validation models, L=1 m, A=1 m², E=1000 Pa and
rho=3 kg/m³. The pendulum uses a fixed 9.81 N downward **nodal force**,
not an automatically mass-scaled body gravity; free mass is 1 or 1.5 kg.
Its initial rotation is 60 degrees and its duration is 1 s. The string has two
segments, free mass 2 or 3 kg, a constant 8 N transverse nodal force, and
initial midpoint displacement 0.25 m; its static equilibrium is 0.2 m.
The physical coefficients are identical between mass choices. Validation
uses dt=0.01, 0.005, 0.0025 s; HHT alpha=-0.1, and a_M=0 or 0.2 /s.
Initial velocities are zero. Every step is returned, with no interpolation.
The independent ODE reference uses RK4 with 100 substeps per output step.
Tests require successive displacement errors to fall below 0.4 of the previous
error, finest-grid component errors below 1% of their reference peaks, energy
balance below 0.1% of the characteristic mechanical/work scale, and Newton
relative residual below 1e-8. HHT acceleration is the algorithm's acceleration
at the nominal time, not a replacement obtained by imposing instantaneous
physical force balance after the weighted HHT solve.


## Native nonlinear black-box input

`deck.write.nonlinear {"path":"motion.inp","step":"Motion",
"mass":"consistent","representation":"uel"}` writes **new** `.inp` and `.c`
files in SI. Both are required. It does not run an external solver. Existing
paths are refused. The generated C is an original user-element implementation,
not a copy of a vendor example. The executable application remains Rust.

A native finite-strain T3D2 with `*ELASTIC` is not this Green-strain/second-Piola
StVK law and is explicitly refused. The generated general user element uses
the same reference nodes, connectivity, E, A and density, and the selected full
consistent or lumped mass. It supplies current internal force and tangent,
initial-acceleration mass/residual, HHT residual/Jacobian and elastic/kinetic
energy through the public UEL interface. This checks a target solver's assembly,
constraints and time integration with our specified constitutive law; it is not
an independent validation against the target's built-in material implementation.

The supported export is deliberately limited to constant dead nodal forces,
homogeneous supports, zero initial velocity, no physical damping and rod2.
Load-amplitude curves, adaptive/modal/base excitation, other elements, point
masses, nonlinear materials and interactions are refused by name. Other options
are not silently linearized. The solver's broader diagnostic support is distinct
from this export subset.

`FC_INITIALIZE` is a fully prescribed geometrically nonlinear static step that
establishes the original initial displacement and the constant applied forces.
`FC_MOTION` retains the physical supports and releases preparation constraints
with step amplitude, using `NLGEOM=YES`, `DYNAMIC,DIRECT,NOHAF`, and explicit
alpha/beta/gamma. It does not disable initial acceleration. Compare dynamic step
time and the unchanged reference displacement definition, not the preparation
pseudo-time. Initial stresses/strain energy are retained. Subtract dynamic-start
external work when comparing a cumulative multi-step target energy record.

Run the generated source through the target's user-subroutine build with
`user=motion.c output_precision=full`. A site-configured C compiler and the
vendor-supplied interoperability header/link environment are required. A separately
validated Linux object can also be supplied via `user=...o`; it is not portable
across platforms or automatically proof of target ABI compatibility. The
standalone C test macro only substitutes the calling-name convention for tests;
it does not include or reproduce a commercial header.

Nodal U/V/A/RF are requested on `FC_NODES`, at the specified output frequency.
The six element state variables are resisting forces at end 1 XYZ, then end 2
XYZ, in N; `EL FILE,SDV` requests them in the results file. `ENERGY FILE` requests
total energies with no data line. User-element stress output in ODB is not
claimed. Independently recovered axial measures are Green strain, second Piola
stress S, nominal stress P=lambda S and axial force A0 P. They are not Cauchy
stress; this one-dimensional law has no independent lateral/volume constitutive
rule. The solver's `recover_rod` API evaluates its real force kernel from saved
nodal displacement without advancing time or inventing a prior stored field.

The full-output execution option selects double-precision Standard field data;
history data and timestamps have separate precision limits. Preserve raw values,
zero frames and failures of a strict time-matching gate. Do not change physical
initial acceleration to imitate a target's zero-output-frame convention.

Public contracts actually consulted (2025 edition; not a claim of a 2024 runtime
check): [truss theory](https://docs.software.vt.edu/abaqusv2025/English/SIMACAETHERefMap/simathe-c-trusses.htm),
[UEL](https://docs.software.vt.edu/abaqusv2025/English/SIMACAESUBRefMap/simasub-c-uel.htm),
[DYNAMIC](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-dynamic.htm),
[user-subroutine C interoperability](https://docs.software.vt.edu/abaqusv2025/English/SIMACAESUBRefMap/simasub-c-subroutineover.htm),
[execution precision](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEEXCRefMap/simaexc-c-analysisproc.htm),
[ENERGY FILE](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-energyfile.htm).
The user-element export itself does not change the unified
07 envelope or kinds 0–8. Native nonlinear histories use the separately allocated kind 11 described above.
