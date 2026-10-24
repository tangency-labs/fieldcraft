---
name: fieldcraft-modelling
description: "Choose and check a FieldCraft structural, heat or fluid analysis: units, material evidence, element idealisation, supports, contacts, steps and quantitative result checks. Use when setting up or judging an analysis, before treating a contour as an engineering conclusion."
---

# Model for the question

Start with the quantity the user needs, its operating conditions and acceptable uncertainty:
deflection, frequency, stress away from a fixture, component temperature, pressure loss or heat
rate. Record the idealisations that could change that answer. Preserve the user's model; make
comparisons in separate saved variants.

Discover the running version's commands and schemas with MCP `list_commands` or
`fieldcraft-cli commands --json`. Use `inspect_model` / `model.info`, `model.activity` and
`model.check` to inspect the existing setup. A command being available does not establish that
every combination of physics, element and material is supported. Read a refusal; do not bypass
it by changing the intended physics. Repository links below refer to the matching checkout.

## Establish the physical model

- **Units:** set or read `model.units` before entering bare numbers. The default is
  `mm-t-s-N-MPa` with Celsius, not SI input. Prefer dimensioned values for material properties;
  read reply `units` and legend/CSV labels. Stored values are SI; changing display units does
  not rescale the physical model. Distinguish absolute temperature, temperature difference,
  dynamic viscosity and kinematic viscosity. See [units](../../docs/units.md).
- **Materials:** record the source, grade/condition, temperature range and uncertainty of each
  property. Example sample values are assumptions, not a certified material library. Elastic
  constants suffice for a static elastic model; density matters for inertia, gravity and modal
  analysis. Transient heat needs density and heat capacity as well as conductivity. Check
  anisotropic directions, thermal expansion/reference temperature, fluid density and dynamic
  viscosity against the operating conditions. Do not infer strength or fatigue allowables
  from Young's modulus. Use public data or user-authorised inputs; no customer portals or
  copied proprietary libraries. [NIST material-data guidance](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication960-11.pdf)
  explains why provenance and conditions matter.

Choose the lowest-dimensional idealisation that retains the quantity of interest:

| Behaviour needed | Candidate and check |
|---|---|
| Axial member | `rod2` with rod area; it cannot represent bending. |
| Slender member bending/torsion | `beam2`; verify profile, orientation and load eccentricity. |
| Thin-wall membrane/bending | Shell section; verify thickness, offset, normals and section stresses. |
| 2D continuum | Plane stress for a thin in-plane-loaded sheet; plane strain for justified out-of-plane restraint; axisymmetry only for axisymmetric geometry **and** loading (FieldCraft: X is radius, Y is axis). |
| Local 3D stress, thick parts, interfaces | Solid elements; resolve the region of interest. Linear tet4 is not a reliable default for bending stress; compare a suitable hex or quadratic formulation. |
| Solid conduction | Heat-transfer step with supported elements/materials; mechanical restraints do not supply thermal boundary conditions. |
| Fluid motion/advected heat | A fluid domain and fluid step on the LBM grid; a structural volume mesh alone does not define the flow grid. |

Use the [meshing skill](../fieldcraft-meshing/SKILL.md) for mesh decisions. Check element/section
compatibility and formulation in the actual schema; reduced integration, incompatible modes
and hybrid behaviour are choices to justify, not interchangeable cures for a bad mesh.

## Make loads and constraints describe the hardware

Draw the load path and free-body balance first. Restrain the actual rigid-body freedoms of
each disconnected structural part without blocking intended Poisson contraction, sliding or
thermal expansion. A singular stiffness matrix asks for a missing connection/support or an
invalid idealisation; adding arbitrary clamps can produce a solvable but different problem.

Apply force over the real attachment footprint when local solid stresses matter. A force
load's default `distribution: "total"` preserves its resultant across a node set; `per_node`
changes the total when the mesh changes. Pressure acts per area (positive inward), so check
face normals, force and moment. Keep support/contact edges and ideal sharp corners separate
from the stress quantity used for assessment: refinement can increase an ideal singular peak
without converging it. Do not hide it by smoothing or by declaring an arbitrary exclusion
distance. Resolve the real fixture/radius, or justify and report an appropriate remote quantity.
[Wood et al., NAFEMS 2015](https://pure.strath.ac.uk/ws/portalfiles/portal/41158000/Wood_etal_NAFEMS_World_Congress_2015_Theoretical_elastic_stress_singularities_much_maligned_and_misunderstood_Jun_2015.pdf)
discuss this dependence on the idealisation.

Use symmetry only when geometry, material directions, loads, supports and the sought response
have that symmetry. `model.boundary_condition.symmetry` constrains normal displacement and
the relevant rotations; it is not a fixed face. Scale total forces to the retained fraction,
but retain pressure/traction intensity. Symmetric modal constraints omit antisymmetric modes;
symmetry can also suppress a physical buckling or flow instability.

Shared conforming nodes represent continuity. `interaction.tie` represents a permanent bond
between nonmatching surfaces; it cannot open or slip. Use supported nonlinear/explicit
`interaction.contact.*` for unilateral closure/friction, checking gap, penetration, normal
direction and contact status. A mechanical tie/contact is not a measured thermal contact
resistance. Fluid conducting-solid interfaces have their own thermal definitions.

## Choose the step, then check its activity

Read `model.activity` and resolve named sets after remeshing. A defined load or boundary
condition acts only where the step activates it (`model.step.activate`).

- Small elastic equilibrium: static; use nonlinear static for relevant large deformation,
  material nonlinearity or changing contact, and inspect converged increments/residuals.
- Natural modes: frequency with justified mass and restraints. Eigenmode display amplitudes
  are arbitrary, not predicted operating displacement. Linear buckling factors are not a
  nonlinear collapse prediction.
- Dynamics/impact: explicit with physical duration, output interval and energy account.
  For a quasi-static explicit approximation, demonstrate negligible inertia for the chosen
  quantity; report mass scaling and added energy rather than assuming a slow-looking movie
  establishes it.
- Conduction: steady heat or transient heat with initial conditions and thermal storage.
  Coupled temperature-displacement is distinct from fluid-solid conjugate heat transfer.
- Fluid: justify domain extent, spacing, walls/inlet/outlet and laminar/turbulent assumptions.
  Inspect the actual backend, Mach/relaxation warnings and stop reason. Reaching a step cap
  does not establish steady state. Accelerated steady iterations with scaled solid heat
  capacity are not physical transient time.

A held-hot-wall, single-air-domain example is not a coupled two-fluid exchanger. Conducting
solids in one fluid domain may exchange heat through the supported conjugate model; separate
hot/cold fluid domains coupled through a wall require an actually available, verified coupling
workflow. Check the current implementation and its verification before promising exchanger
effectiveness. Likewise, do not infer a GPU LES/fan/porous capability from the CPU option or
an open pull request. Use [current verification](../../docs/verification.md) and
[oracle coverage](../../docs/oracle.md), not roadmap aspirations.

## Decide whether the answer is usable

Before solving, resolve `model.check` errors and explain each relevant warning. After solving:

1. Compare an independent estimate with a consistently measured quantity; use
   [small verification exercises](references/checks.md) when learning the workflow.
2. Check structural force **and moment** balance including reactions; check heat sources,
   boundary flows and storage; check fluid mass/energy balance including inlet/outlet fluxes.
   A transient need not have equal inlet and outlet heat while energy is accumulating.
3. Inspect `solve.probe`, reactions/histories and units as well as the contour. State the
   deformation scale, time/frame, averaging and location. A fluid wall-heat value is positive
   **into the wall**; a thermal input boundary's positive flux is into the model/fluid.
4. Refine the mesh/grid and, for transients, time resolution independently while retaining
   physical loads, geometry and duration. Compare the same quantity and location; agree an
   error target before accepting the change. One refinement is a sensitivity check, not proof
   of asymptotic convergence; a singular peak cannot be rescued by a percentage criterion.

Report the quantity with units, source assumptions, balance/refinement evidence, warnings and
remaining uncertainty. Code/oracle agreement verifies an implementation; matching a physical
experiment under matched conditions validates that modelling use. The
[NASA V&V tutorial](https://www.grc.nasa.gov/www/wind/valid/tutorial/tutorial.html) distinguishes
these checks. If evidence is insufficient, say which decision the current result can support
and what measurement or model refinement would resolve the uncertainty.
