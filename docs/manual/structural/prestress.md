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

```json
[
 {"command":"model.pretension.create","params":{"name":"Bolt","first":"cut1","second":"cut2","normal":[0,0,1],"control":"BoltControl"}},
 {"command":"model.pretension.force","params":{"name":"Bolt","step":"Preload","value":1000}},
 {"command":"solve.nonlinear","params":{"step":"Preload"}},
 {"command":"model.pretension.lock","params":{"name":"Bolt","step":"Service","from":"Preload"}},
 {"command":"solve.nonlinear","params":{"step":"Service"}},
 {"command":"solve.pretension","params":{"name":"Bolt","step":"Service"}}
]
```

The named steps, cut sets, reference point and other model definitions must
already exist. `model.pretension.adjust` prescribes q directly (`value`, length).
A section without an explicit control has q=0; a perturbation holds the increment
of q at zero. `lock` snapshots the adjustment of a current solved earlier step
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
constant jump; no automatic bolt cutting, nonconforming pretension cut, or evolving
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
