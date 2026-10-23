# Two-fluid heat exchanger: implementation contract for 0.7.5

This is the proposed minimum for [#700](https://github.com/tangency-labs/Fieldcraft/issues/700),
under [release acceptance #697](https://github.com/tangency-labs/Fieldcraft/issues/697).
It specifies work still to implement. The existing `cooling-aio` is an air-side study with a held
radiator temperature. Two domain definitions and two independent solves do not exchange heat.

## First coupled case

Use two straight, counterflow parallel-plate passages separated by one conducting plate.
All dimensions and stored quantities below are SI. These are deliberately prescribed benchmark
inputs, not a material-library claim or an experiment. There is no fan, pump, buoyancy, radiation,
phase change, turbulence closure or contact resistance in this first case.

| Item | Definition |
|---|---|
| Length and repeated width | x = 0–0.05 m; z = 0–0.002 m, periodic in z |
| Air passage | y = 0–0.002 m; flow along +x; inlet 293.15 K |
| Separator | y = 0.002–0.003 m; one body named `Separator` |
| Coolant passage | y = 0.003–0.005 m; flow along −x; inlet 323.15 K |
| Air properties | ρ = 1.2 kg/m³, μ = 1.8e−5 Pa·s, cp = 1000 J/(kg·K), k = 0.024 W/(m·K) |
| Coolant properties | ρ = 1000 kg/m³, μ = 0.001 Pa·s, cp = 4200 J/(kg·K), k = 0.6 W/(m·K) |
| Separator properties | ρ = 8000 kg/m³, cp = 500 J/(kg·K), k = 15 W/(m·K) |
| Mean inlet speeds | air 0.05 m/s; coolant 0.0005 m/s; fully developed parabolic profiles |
| Other boundaries | no-slip, insulated outer plates; zero gauge pressure at each outlet; insulated separator end faces |
| Refinement | cubic spacing 0.00025, 0.000125 and 0.0000625 m; all interfaces grid aligned |

Using hydraulic diameter 2H = 0.004 m gives Re = 13.33 and 2 respectively: the reference is
laminar. Prescribe the *mean* speeds above; a parabolic inlet's peak command value is 1.5 times
its mean. Keep each circuit's mass, momentum, density scale and material independent. The shared
plate is not a coolant region and cannot carry either fluid's momentum populations.

For an independent hydraulic check use two interior stations 0.03 m apart in the developed
region: Δp = 12 μ U_mean L/H² gives 0.081 Pa for air and 0.045 Pa for coolant. Exclude entrance
and exit losses; report station positions and extrapolation/interpolation. The saved port-layer
static pressures are useful diagnostics but do not replace interior pressure taps.

## Ownership and conservative transfer

Start with aligned, equal-spacing CPU grids and a common admissible fluid time step. A global
cell key is the integer world-grid coordinate. Each fluid cell has exactly one circuit owner;
each separator cell has exactly one thermal owner. Reject overlapping fluid occupancy, gaps
in the intended interface, incompatible grid origins/spacing, duplicate wall definitions or
different materials/powers assigned to the same wall. A separator's capacity and source are
allocated once. Each circuit has its own stable lattice scaling: reject a common time step that
fails either fluid's existing stability checks, with both domains named.

Extend the existing conservative CHT link operation to exchange *energy in joules* across circuit
boundaries. For every fluid–wall link, accumulate the transfer once and apply equal and opposite
increments to fluid and wall. Convert the fluid increment using that circuit's ρ cp ΔV, and the
wall increment using its own ρ cp ΔV. Do not copy a lattice temperature or a lattice heat flux
between differently scaled circuits. Use the existing half-cell harmonic conduction resistance,
with both sides' conductivities and the actual link area; verify the discrete transfer against
the existing single-fluid CHT implementation before joining the domains.

The first orchestrator advances both fluid states and the *single* wall state together. Bound
the explicit wall increment by its capacity divided by the sum of incident conductances and
sources' existing constraints; reduce the common time step when necessary. No unbounded outer
iteration, manual flux handoff, duplicate solid solve or hidden relaxation is permitted. Initially
use relaxation factor 1 and physical heat capacities. Any later under-relaxation must preserve
the equal/opposite transferred energy, report its factor, and pass a tighter-tolerance study.

The first public support is **CPU + CPU, steady state only**. Reject transient requests, GPU
and mixed backends by name until verified. Iteration playback is labelled convergence toward
steady state. Accelerated solid capacity is disabled initially; if added later it is saved and
labelled, and is never presented as a physical warm-up.

## Stop rules, budgets and reference hierarchy

Require both circuits' field-change monitors and their mixed outlet temperatures to settle over
two successive windows. Require relative outlet-temperature change below 1e−4 of the inlet
temperature difference and each circuit's relative reference-density flow mismatch below 1e−3.
Require the combined energy residual below 1e−3 of the larger circuit duty; use a documented
absolute floor for the zero-duty case. Set finite time/iteration caps and retain results with an
explicit `max_time` verdict when those checks do not pass. A tolerance is an acceptance proposal,
not a measured accuracy claim.

Compute inlet/outlet mixed temperatures from the directional integral of ρ cp u·n T, not an
arithmetic average of cells. Report backflow separately. For constant cp and balanced flow:

```
Q_hot  = m_dot_hot  cp_hot  (T_hot_in  - T_hot_out)
Q_cold = m_dot_cold cp_cold (T_cold_out - T_cold_in)
```

When flow is imbalanced or reversed, use the separately integrated directional sensible-heat
rates and disclose their temperature reference. Expose Q_hot and Q_cold independently, wall
storage rate, applied power, external wall losses, conductive flux through open ports, and the
signed residual. Sum each internal link once per owner so it cancels from the combined budget.
The new `ports` reports provide **advective transport only**, using reference-density mass flow
from saved fields. They do not yet measure port conduction, wall storage or a full energy budget.

Check the coupling arithmetic first against an independently prescribed-conductance network.
For a counterflow network with C_hot = 2 W/K, C_cold = 1 W/K, UA = 1 W/K and inlet temperatures
343.15 K / 293.15 K, let r = C_min/C_max and N = UA/C_min:

```
epsilon = (1 - exp(-N*(1-r))) / (1 - r*exp(-N*(1-r)))
Q = epsilon * C_min * (T_hot_in - T_cold_in)
T_hot_out = T_hot_in - Q/C_hot
T_cold_out = T_cold_in + Q/C_cold
```

For these inputs epsilon = 0.5647334016064162, Q = 28.23667008032081 W,
T_hot_out = 329.0316649598396 K and T_cold_out = 321.3866700803208 K.

The independently specified network must match the closed form within 0.1% and conserve energy
to round-off per exchange. Also cover r = 1 (epsilon = N/(1+N)), UA = 0 and equal inlet
temperatures. This is a coupling-unit reference, **not** a CFD fixture: do not choose UA by fitting
the CFD solution and then call the comparison validation. The planar CFD case additionally needs
spatial and time-step refinement, tighter stop tolerances, the pressure-gradient check, existing
single-domain regressions and the closed insulated two-fluid/wall energy test. Keep #623's
experimental forced-convection validation and #631's curved-interface uncertainty separate.

## Commands, results and release evidence

Add an explicit coupled-system model definition, validation and solve command through the normal
engine catalogue. Existing `solve.fluid` remains a single-domain solve. Before allocating, show
the sum of both fluid states, the wall owned once, coupling buffers, staging and output frames
in the memory estimate; all allocations and work have checked caps.

A coupled result groups both named circuits and one wall under one solve revision, retains
solve-time material properties, common output times, stop settings/verdict and budget histories.
Save/reopen that group without solving. Probes/CSV identify their circuit and wall, their frame
and units. Selecting temperature, velocity, streamlines or a pressure section must name which
circuit it displays; stale results keep their original properties.

The release model must include an original geometry recipe, exact source commit, model and
saved results, resolution/time step, convergence settings, hardware and wall time. Its in-app
video shows setup, both results, budgets and save/reopen from that same run. Each implementation
PR includes a shorter in-app clip of its actual feature; a clip of port diagnostics or the existing
AIO does not claim to demonstrate a coupled exchanger. Track remaining acceptance in #697/#700.

Sources: [MIT Unified Engineering, §18.5](https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node131.html)
(two-fluid energy accounting and wall resistance), and the
[EnergyPlus engineering reference, Fluid-to-fluid heat exchangers](https://energyplus.readthedocs.io/en/stable/guides/engineering-reference/16.8-heat-exchangers.html)
(effectiveness–NTU equations). Equations and benchmark inputs here are our own reference setup;
no external solver input or commercial geometry is used.
