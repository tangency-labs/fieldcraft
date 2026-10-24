# Small checks before a design claim

These are original, synthetic verification exercises, not characterised engineering materials
or validation against experiment. Open a fresh scratch model for each recipe; running one in
an existing model adds definitions and may collide with names. JSON files use the normal
FieldCraft script format (`fieldcraft-cli run <recipe.json>` or individual MCP `execute`
commands). Read the command replies, not only the picture.

## Axial force, units and supports

Run [axial-rod.json](axial-rod.json): a 0.3 m rod, area 0.0001 m², assumed Young's modulus
70 GPa, and a total 1000 N axial force. The rod idealisation has no transverse stiffness;
transverse degrees of freedom are restrained, while its root alone is held axially. This
would be an inappropriate model for a bending question.

For small linear extension, `u = F L / (E A) = 4.285714285714286e-5 m`,
`sigma = F/A = 10 MPa`. Expect the root reaction to be −1000 N along X, balancing the load.
Check `solve.probe` at both nodes. Change only the display to `mm-t-s-N-MPa`: extension must
read 0.04285714285714286 mm and the same physical stress 10 MPa. A change of units must not
change the physical model. This is an exact 1D check, not evidence of a converged local 3D
stress field. The [University of Illinois axial-loading reference](https://mechref.engr.illinois.edu/sol/axial.html)
gives the underlying relation and its assumptions.

## Conduction, heat-flow sign and spatial sensitivity

Run [steady-conductor.json](steady-conductor.json): an insulated-sided 0.12 m long prism with
0.02 × 0.03 m cross-section, assumed constant conductivity 15 W/(m K), and prescribed end
temperatures 293.15 K and 353.15 K. No heat capacity or mechanical clamp is required for this
steady pure-conduction problem.

Fourier conduction gives a linear temperature: midpoint 323.15 K,
`q_z = -k (T_hot - T_cold)/L = -7500 W/m²`, and heat-rate magnitude 4.5 W. Check the probe's
actual node coordinates and distance. Probe the `Cold` and `Hot` sets separately and sum their
nodal heat reactions: they must be opposite, each of magnitude 4.5 W. A net reaction near
zero is a balance check, not the through-flow magnitude.

Repeat in a fresh model with block divisions changed from `[2,3,6]` to `[4,6,12]`; compare
temperature, flux and signed end rates. This linear patch is exactly representable on both
meshes, so equality checks the setup/implementation, not an observed convergence order for
curved geometry or thermal boundary layers. For a design, use a meaningful quantity on a
mesh family; the [NASA spatial-convergence discussion](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html)
explains why three levels and an asymptotic-range check are more informative than one change.

## Heated flow: inspect the approximation before solving

`file.sample {"name":"heated-flow"}` opens the built-in bracket example. Inspect the
`WarmAir` fluid domain and `HeatedFlow` step with the matching command schemas: it has one
air domain, an imposed inlet and held-hot obstacle faces, 5 mm lattice spacing and 180 steps.
The solid is not cooling, there is no second fluid circuit, and this sample has no conjugate
solid conduction or buoyancy. It demonstrates transient heat carried by air around an
imposed-temperature obstacle; it does not establish exchanger effectiveness or steady state.

Set a deliberate display system, run `model.check`, then `solve.fluid` with an explicit
supported backend. Inspect the returned time step, physical end time, warnings and stop
status. Use `results.slice` with `field: "temperature"` and `results.surface` with
`field: "wall_heat_flux"`. Negative heat **into** the held-hot bracket is expected when it
supplies heat to cooler air. A temperature range within the prescribed bounds is a smoke
check; it does not establish an accurate heat rate, mass/energy balance or grid convergence.
Do not use a coarse demonstration to certify a cooling design.

For each use, retain the actual version/backend, settings, units and measured quantities.
Choose design-specific error targets and uncertainty sources before refinement. None of
these exercises establishes contact, plasticity, turbulence, fatigue or multi-fluid coupling
accuracy merely because a solver offers those options.
