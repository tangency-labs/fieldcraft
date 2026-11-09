# Prestress and bolt pretension

A frequency or harmonic step may name an earlier static or nonlinear static step
with `base_state`. The solver recomputes that equilibrium from current inputs;
changing its definition, loads, boundary conditions or sets invalidates the
perturbation result. Perturbation boundary conditions remain explicit.

```json
{"command":"model.step.create","params":{"name":"Loaded modes","base_state":"Preload","procedure":{"type":"frequency","modes":10},"boundary_conditions":["Supports"]}}
```

A linear static base supplies integration-point stress stiffness Kσ, added to K
on the reference mesh (Bathe 2014 §6.8.2). Loads are dead loads. Solids and the
existing flat-shell operator are supported; beams and trusses are refused.
A finite-deformation nonlinear base supplies its converged consistent tangent
while material history is still available; Kσ is not added twice. Mass uses the
reference mesh. Small-strain nonlinear bases, geometric imperfections and plastic
unloading tangents are currently refused.

At a nonlinear terminal state, mortar pairing and the positive-pressure active
set are frozen. Closed constraints become bonded in all translations, including
nodes that previously slipped; open constraints are absent. The original contact
penalty tangent is excluded. This is an explicit bonded-when-closed policy, not a
Coulomb slip tangent. Existing ties keep the same reference mapping used by the
base solver. The capture is used within the solve; a reopened model recomputes
its equilibrium rather than interpreting output stress averages as restart history.

## Pretension cuts

`model.pretension.create` defines two coincident, conforming planar cut-node sets,
a fixed unit normal and a dedicated existing reference point. Cut nodes must be
duplicated across the cut. For every paired node,
`u_first − u_second = normal × q`; the reference point's Ux is the scalar q,
not a physical global X translation. Its other directions are inactive.
A positive force conjugate to q tightens the bolt halves.

For a continuous solid mesh, the optional `face` names a planar internal
element-face set. Its normal must point away from the selected underlying
elements. The source keeps shared nodes; static and perturbation solvers prepare a
private copy before assembly, duplicate the negative-normal trace (including
quadratic mid-side nodes), and rewire every incident element on that side.
The existing node-pair constraints and exact elimination then enforce tangential
continuity and `u_minus - u_plus = normal * q`. `first` and `second` remain named
sets for compatibility but do not impose end-face displacement constraints in
this mode. In particular, q can be nonzero with both outer ends fixed.

```json
{"command":"model.pretension.create","params":{"name":"Bolt","first":"Bottom","second":"Top","face":"InternalCut","normal":[0,0,1],"control":"BoltControl"}}
```

Frequency and harmonic analyses assemble mass, capture the equilibrium tangent
and freeze constraints on that same private mesh. The public mode and phasor
fields are projected back to source nodes only after solving. Boundary faces, intersecting cuts,
incomplete sections and nonplanar cuts are refused. Public nodal output retains
the positive-side trace at source cut nodes; q reports the jump separately, and
element stresses are retained on both sides. Private duplicate labels are never
written into the source model or its native deck. The native face export uses
the shared connectivity and `*PRE-TENSION SECTION`, without transverse equations.

```json
[
 {"command":"model.pretension.create","params":{"name":"Bolt","first":"cut1","second":"cut2","normal":[0,0,1],"control":"BoltControl"}},
 {"command":"model.pretension.force","params":{"name":"Bolt","step":"Preload","value":1000}},
 {"command":"solve.nonlinear","params":{"step":"Preload"}},
 {"command":"model.pretension.fix","params":{"name":"Bolt","step":"Service"}},
 {"command":"solve.nonlinear","params":{"step":"Service"}},
 {"command":"solve.pretension","params":{"name":"Bolt","step":"Service"}}
]
```

The named steps, cut sets, reference point and other model definitions must
already exist. `model.pretension.adjust` prescribes q directly (`value`, length).
A section without an explicit control has q=0; a perturbation holds the increment
of q at zero. `fix` resolves the actual converged q of the preceding static or
nonlinear static step from current model inputs during each solve, skipping
intervening perturbation steps. Chains of current-value controls are evaluated
in analysis order, with at most 1024 predecessor equilibria. The source keeps
the `fixed` control, and native export writes `*BOUNDARY, FIXED`; no old result
or prescribed numerical q is substituted in the source. Changes to the predecessor
invalidate the locked equilibrium and its perturbations, including saved results.
This elastic re-equilibration does not transfer plastic or frictional history;
those current-value chains are refused. Frictionless terminal contact is supported.

Alternatively, `lock` snapshots the adjustment of a current solved earlier step
into the target. After changing the preload, solve it and lock again. The journal
records the numerical adjustment with its unit, so replay needs no transient
results. Locking does not transfer plastic history; plastic materials are refused
by `lock`. Elastic later steps are solved anew with the locked adjustment and
their own explicit loads and boundary conditions.

`solve.pretension` reports adjustment, applied scalar force and section force
recovered from equilibrium by virtual work through the constraint transformation.
The scalar load is an internal force pair and is excluded from global force sums.
Definitions, controls and section results are serialized with the model/results.

Current limitations: conforming coincident planar cuts only, fixed normal and
constant jump; no automatic selection of a bolt cut, nonconforming pretension cut, or evolving
cut normal. Keep the cut control point separate from ordinary loads, springs,
couplings and boundary conditions. The Abaqus writer emits native `*PRE-TENSION SECTION` for conforming cuts whose element faces can be resolved, writes the control-node force in the preload general step, and carries a locked adjustment with a nonzero control-node boundary value. Frequency and harmonic steps with `base_state` are emitted as perturbation steps; a general static step feeding such a perturbation is written with `NLGEOM=YES` so Abaqus includes preload/initial-stress stiffness. Nonconforming or face-less cuts, unsupported native procedures, and pretension on a perturbation step are refused by name. Nastran and LS-DYNA writers continue to refuse pretension and `base_state`.
Sensitivity paths that assume unloaded stiffness refuse prestressed steps.

Verification sources: Bathe, *Finite Element Procedures*, 2nd ed., 2014 §6.8.2;
Cook et al., *Concepts and Applications of Finite Element Analysis*, 2002 §9.2;
Puso & Laursen, CMAME 193:601–629, 2004 (mortar integration). Regression scenarios
cover axial-force beam frequencies, cut-force balance and locking, open/closed
contact, persistence/journal replay and a bolt/clamp command sequence. The axial column, cut-force/lock, open/closed contact and command-chain regressions
have passed on Linux. Full-suite/oracle acceptance and BASE-scale validation are
separate checks. Contact cut closure includes the permitted terminal penetration;
report this gap separately when comparing with an ideal hard-contact formula.

`internal_pretension_chain` checks a uniform axial mode with one free coordinate:
two unit cubes with E=1e6 Pa, nu=0, rho=1 kg/m³ give consistent reduced mass 2/3 kg.
The linear prestress stiffness is 2(E+F); the finite StVK tangent is E(3 lambda²−1),
where lambda=1+q/2. Direct and modal phasors are compared with
1/(K−omega²M+i omega alpha M). This checks the stress contribution without adding
K_sigma twice. Zero preload recovers the unloaded operators.

Frequency results currently retain the Lanczos basis size, not individual modal
residuals or iteration histories. Missing diagnostics must be reported as missing;
absence of warnings is not a residual measurement. Harmonic frames retain full
and projected residuals and power balance; nonlinear history retains Newton
residual measures, iterations, cutbacks and converged increments. Modal derivative
paths requiring an unprojected physical mass operator refuse private cuts.
