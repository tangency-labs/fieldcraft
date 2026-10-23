# Inlet and outlet transport

[In-app demonstration recipe, settings and measured report](port-transport-video.md).
The recording and still are intended as attachments on [PR #704](https://github.com/tangency-labs/fieldcraft/pull/704), outside git.

After a fluid solve, expand **Results → Inlets and outlets**. Each open box face shows separate
inflow and outflow rates, mixed temperature and flow-weighted static gauge pressure. It follows
the displayed frame, including when a saved solved model is reopened. Backflow is retained on
its own side of the table. No flow has no mixed temperature or pressure (shown as a dash).

`solve.fluid` returns the final frame under `ports`; `results.frame`, `results.contour` and the
other results-view replies return the displayed frame there. Its `units` map describes the
display units, and all conversions use the model's display system, including absolute
temperature offsets. `faces` is in box-face order and contains `face`, fluid `area`, `inflow`
and `outflow`; each direction has `volume_flow`, `mass_flow`, `temperature`, `pressure` and
`sensible_heat`. Net rates are positive outward.

Each fluid cell on a boundary layer contributes |u·n| Δx². Solid/obstacle cells do not contribute
area or transport. Mass flow uses the fluid's **reference density**, consistent with the
incompressible interpretation; it is not an exact weakly compressible lattice mass-flux budget.
Temperature and pressure are weighted by the directional flow. Material properties are retained
from solve time, so an edit that makes a result stale cannot change its reported mass or heat.
Older saved results retain volume transport and mixed temperatures but show missing mass/heat
as null/dashes until solved again.

Sensible heat is ρ cp ∑ |u·n| Δx² (T − T_reference), with the reference explicitly displayed.
`net_sensible_heat` subtracts incoming from outgoing transport. This is an advective diagnostic:
it excludes port conduction, storage, wall heat and sources, and makes no steady-state or
energy-balance pass/fail claim. A boundary-layer pressure can include an outlet artifact; measure
pressure loss at appropriate interior stations. No turbulent time averaging is implied by a
single saved frame.

For the two-fluid exchanger contract and remaining implementation, see
[heat-exchanger.md](solvers/heat-exchanger.md) and [#700](https://github.com/tangency-labs/Fieldcraft/issues/700).
