# Geometric nonlinear transient diagnostics

`solve.transient.nonlinear {"step":"Motion","mass":"consistent",
"max_iterations":30,"tolerance":1e-8}` integrates an existing transient step with
full Newton iterations. It is an experimental **rod2** path: finite rotations and
St Venant–Kirchhoff axial elasticity, homogeneous fixed supports, dead nodal
forces (including their amplitude curves), consistent or lumped reference mass,
optional point masses and mass-proportional physical viscous damping.

The reply contains SI nodal displacement, velocity, acceleration and T, V, W,
physical damping dissipation, energy balance, Newton correction counts and
residuals. It is ephemeral diagnostics: save the command reply separately.
It neither creates nor replaces a stored linear result. `file.save_as` does not
persist these nonlinear frames. No native result kind or payload is added.
Native nonlinear persistence requires a separately approved version contract.
The existing transient step supplies time controls; the nonlinear selection is
the explicit command above. Reproduce with the full command script. Saving
only the model or exporting that step with `deck.write` does not turn it into
a nonlinear analysis; a native nonlinear deck-export path is not delivered here.
The coordinated fatigue merge supplies unified MCRSLT07 for existing linear
results. This diagnostic command adds no persisted nonlinear kind or payload;
existing 06 archives are not 07 deliverables.

Solids, shells, beams, rod3, contacts, couplings, foundations, springs/dashpots,
plasticity/hyperelasticity, material damping, follower loads, base motion,
nonzero prescribed displacements, requested stress/other field output, modal reduction, adaptive steps and
stiffness-proportional damping are refused. Point masses must provide a positive
mass operator on the free coordinates. This rod formulation is not a tension-only
cable law: compression is elastic, and collapse is refused. It does not extend
the supported elements of the separate nonlinear static solver.

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
  page-by-page check. HHT formulas use the original Hilber–Hughes–Taylor (1977),
  [doi:10.1002/eqe.4290050306](https://doi.org/10.1002/eqe.4290050306), already
  used by the linear integrator. The nonlinear effective tangent above is
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
