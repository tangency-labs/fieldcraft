# Free surface

A fluid domain with a free surface holds liquid in part of its box and gas elsewhere, and the
liquid moves through the gas: a column collapsing (a dam break), a cavity filling from a gate, a
wave running up a wall. The gas is at rest at the gauge pressure 0 and carries no momentum, so the
surface is a free-slip, constant-pressure boundary that moves with the liquid. Velocity inlets
inject liquid, so a filling can start from an empty box. The solve tracks the liquid's mass exactly
and reports the balance, the fill level of every cell, and the liquid's extent over time.

| | |
|---|---|
| **Commands** | `model.fluid_domain.create` / `.edit` (`free_surface`), `solve.fluid` (`free_surface` in the reply), `solve.probe` (`line` with `field: "fill"`, `front`, `depth`) |
| **Step procedure** | `"fluid"` on a domain with `free_surface` settings |
| **Solver** | `crates/lbm/src/surface.rs` (CPU, f64); not on the GPU backend yet |
| **Method** | Volume of fluid on the lattice with interface cells: mass exchange, reconstruction from the gas and cell conversion (Körner et al. 2005; Thürey 2007) |
| **Verified by** | oracle `free-surface-hydrostatic-column`, `free-surface-plug-filling` (theory), `dam-break-martin-moyce` (experiment); `crates/lbm/tests/free_surface.rs`, `checkpoint.rs` |

## Synopsis

```json
[
  {"command": "model.units", "params": {"system": "m-kg-s-N-Pa"}},
  {"command": "model.material.create", "params": {"name": "Liquid", "density": 1000,
    "fluid": {"type": "newtonian", "viscosity": 0.2}}},
  {"command": "model.fluid_domain.create", "params": {"name": "Tank", "material": "Liquid",
    "min": [0, 0, 0], "max": [0.286, 0.00114, 0.171], "spacing": 0.00114,
    "acceleration": [0, 0, -9.81],
    "faces": {"y_min": {"type": "periodic"}, "y_max": {"type": "periodic"}},
    "free_surface": {"liquid": [{"min": [0, 0, 0], "max": [0.057, 0.00114, 0.114]}]}}},
  {"command": "model.step.create", "params": {"name": "Collapse", "procedure":
    {"type": "fluid", "domain": "Tank", "end_time": 0.19, "time_step": 3.1e-5, "frames": 6}}},
  {"command": "solve.fluid", "params": {"step": "Collapse"}},
  {"command": "solve.probe", "params": {"step": "Collapse",
    "front": {"axis": "x", "at": [[0.082, 0, 0], [0.171, 0, 0]]}}}
]
```

A column a = 57 mm wide and 2a high collapses in a tank 5a long under gravity, one cell deep and
periodic across (a plane problem). The probe gives the surge front's position over time and when
it passed 1.44a and 3a.

## Parameters

### The domain (`model.fluid_domain.create`, `.edit`)

| Field | Value | Default | Meaning |
|---|---|---|---|
| `free_surface` | `{…}` | none | A free surface. Left out, the box is full of one fluid. |
| `free_surface.liquid` | `[{"min": [x, y, z], "max": [x, y, z]}, …]` | `[]` | The regions that hold liquid at the start: the cells whose centres they contain. Everything else is gas. At most 1000; each must contain a cell centre. Empty when every drop comes in through an inlet. |

With a free surface the domain takes walls, moving walls, slip faces, periodic faces and velocity
inlets on its box faces, obstacles with walls and moving walls, an acceleration (gravity) and an
initial velocity, and a Newtonian fluid. Everything else is refused by name as not supported yet:
pressure faces, heat (`thermal`), a turbulence model, porous zones and fans, a shear-dependent
viscosity, and the GPU backend. Surface tension is not modelled.

A velocity inlet injects liquid: its cells are liquid for good and the Zou–He condition sets the
velocity. Nothing leaves the box (there are no outlets yet), so a filling ends when the box is full
or the step ends.

### The reply of `solve.fluid`

`free_surface`: the liquid's `initial_mass`, its `mass` at the end, what the inlets `injected`, the
`mass_error` of that balance, `(mass − initial_mass − injected)` over the larger of the initial and
injected mass (round-off, about 10⁻¹³, when all is well: a warning names a larger one), the `fill`
share of the fluid cells' volume that is liquid, the counts of `liquid_cells`, `interface_cells`
and `gas_cells` at the end, and how many cells `filled` and `emptied` over the run.

### Probing (`solve.probe`)

- `{"step", "line": {…, "field": "fill"}}`: the fill level along a line, 1 in the liquid, 0 in
  the gas, between in the interface cells.
- `{"step", "front": {"axis": "x", "end": "max", "at": [[x, y, z], …]}}`: the liquid's extent
  along the axis at every history time (`times`, `extent`: the `max` or `min` face of the cells
  at least half full), its `volume`, and per point in `at` the time the extent first reached the
  point's coordinate along the axis, coming or going (`arrivals`; between two samples by linear
  interpolation; null when it never did). The surge front of a dam break is `axis: "x"`,
  `end: "max"`; the column's falling top is `axis: "z"`, `end: "max"` with points at the heights
  of interest.
- `{"step", "depth": {"axis": "z", "at": [[x, y, z], …]}, "frame": k}`: per point, the liquid in
  the column of cells through it along the axis in that frame, the fill levels summed times the
  spacing: a wave gauge's depth, or the height of a collapsing column at its wall.
- `points` and the other line fields read the liquid's and the interface cells' velocity and
  pressure; a gas cell has velocity 0 and pressure 0.

The frames' fill levels are saved with the results. The viewport does not draw the surface yet:
`results.slice` with `pressure` or `velocity` shows the liquid, since the gas has neither.

## Theory

Every fluid cell is **gas**, **interface** or **liquid**. Liquid cells run the ordinary lattice
Boltzmann collision and streaming, and their mass is their density. Gas cells hold no populations
and do nothing. Interface cells, the one-cell layer between the two, collide and stream too, and
carry a mass m apart from their density ρ; their fill level ε = m / ρ says how much of the cell is
liquid (Körner et al. 2005; Thürey 2007, ch. 4).

**Mass exchange.** After streaming, an interface cell's mass changes by
Δm = Σ_i (f_ī(x + c_i) − f_i(x)) over its liquid neighbours and by the same weighted with
(ε(x) + ε(x + c_i)) / 2 over its interface neighbours; nothing passes to a gas cell or through a
wall. Each exchange is antisymmetric between its two cells, so the total mass, the liquid cells'
densities plus the interface cells' masses, is conserved to round-off. An inlet's mass enters
through its face: the incoming populations the boundary sets less the outgoing ones that left.

**Reconstruction.** The populations an interface cell would receive from a gas cell at x + c_i do
not exist. They are set to f_ī = f_i^eq(ρ_A, u) + f_ī^eq(ρ_A, u) − f_i(x): the equilibria at the
gas density ρ_A (the gauge pressure 0) and the cell's velocity, less the population leaving toward
the gas. That makes the surface a free-slip boundary at the gas pressure. The populations arriving
from the gas side of the surface normal (n · c_i > 0, n pointing from the liquid into the gas, from
the central difference of the fill levels) are reconstructed too, so the balance at the surface
does not depend on which cells happen to be gas.

**Conversion.** An interface cell whose mass exceeds (1 + κ) ρ, or that has no gas neighbour
left, becomes liquid; one whose mass falls below −κ ρ, or that has no liquid neighbour left,
becomes gas (κ = 10⁻³ keeps a cell from converting back in the next step). A cell that fills turns
its gas neighbours into interface cells (mass 0, populations at the equilibrium of their populated
neighbours' mean density and velocity); a cell that empties turns its liquid neighbours into
interface cells (mass: their density) and drops its populations. A cell next to one that fills does
not empty in the same step. The mass a conversion leaves over, m − ρ of a filled cell and m of an
emptied one, goes to the interface neighbours weighted by n · c_i toward the gas for a filled cell
and away from it for an emptied one, evenly when no weight is positive, and into a reservoir shared
out over every interface cell when there is none. Cells next to an inlet never empty.

**Time step.** With a free surface and an acceleration the solver's characteristic velocity is the
free-fall speed over the box's largest extent, √(2 a L), so the automatic Δt keeps a falling
liquid at lattice speed 0.05. Give `time_step` to control τ: the viscosity the lattice can hold
is ν ≥ (τ_min − ½) Δx² / (3 Δt) with τ_min = 0.505, which for water on millimetre cells is a
hundred times its viscosity. Dam breaks and fillings are driven by gravity and inertia; their
fronts are insensitive to the Reynolds number once it is in the hundreds, but the grid refinement
in the dam-break case's README shows what the viscosity and the resolution cost.

**Determinism.** The mass exchange and the reconstruction run by fixed blocks of cells in
parallel, each cell writing only its own mass and populations; the conversions run on one thread
in cell order. The result is the same bit for bit on any number of threads, and a checkpointed run
resumes to the same bits.

## Verification

| Case | Kind | What it checks |
|---|---|---|
| `free-surface-hydrostatic-column` | theory | A column at rest under gravity: ρ g (h − z) at the bottom, the surface staying, the liquid holding still, the mass to round-off. |
| `free-surface-plug-filling` | theory | A frictionless channel filled from an inlet: the front at U t + Δx, the liquid's volume, the mass balance. |
| `dam-break-martin-moyce` | experiment | Martin and Moyce's collapsing column (1952): the surge front's travel times between stations and the column's fall, against their Table 2 and Table 6. |

## Sources

- C. Körner, M. Thies, T. Hofmann, N. Thürey, U. Rüde, *Lattice Boltzmann model for free surface
  flow for modeling foaming*, J. Stat. Phys. 121(1–2):179–196, 2005.
- N. Thürey, *Physically based animation of free surface flows with the lattice Boltzmann method*,
  PhD thesis, Universität Erlangen-Nürnberg, 2007, ch. 4.
- N. Thürey, U. Rüde, *Stable free surface flows with the lattice Boltzmann method on adaptively
  coarsened grids*, Comput. Visual. Sci. 12(5):247–263, 2009.
- J. C. Martin, W. J. Moyce, *Part IV. An experimental study of the collapse of liquid columns on a
  rigid horizontal plane*, Phil. Trans. R. Soc. A 244:312–324, 1952.
