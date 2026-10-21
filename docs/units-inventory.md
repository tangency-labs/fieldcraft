# Units inventory (#579)

*The state before SI storage, kept as the map the conversion followed; docs/units.md describes
how units work now.*

Every place a dimensional number enters or leaves FieldCraft today, the unit it is assumed to be
in, and where the SI conversion goes. Today nothing converts: every stored `f64` means "whatever
consistent system the user chose", labelled by `model.units` (`UnitSystem`,
`crates/model/src/units/mod.rs`). The SI table that later steps convert through is
`crates/model/src/units/` (`Quantity`, `Unit`, `DisplaySystem`). **Risky** marks an implicit
assumption, a hard-coded factor or a missing dimension.

## Model storage (`crates/model/src`): all bare `f64`, no unit tag
Conversion goes at the command/file boundary; storage becomes SI.

| Where | Fields | Quantity |
|---|---|---|
| `lib.rs` | `units: UnitSystem` (serde default `unspecified`); `MAX_COORD = 1e9` "model units" | label; length |
| `mesh.rs`, `datums.rs`, `constraints.rs`, `relaxations.rs`, `mesh_plans.rs` | node/datum/reference point xyz, plan centre/radius/size/tolerance | length |
| `sizes.rs`, `seeds.rs`, `virtual_faces.rs`, `tie.rs` | size regions, seed sizes, tolerances, position_tolerance | length (bias, growth: dimensionless) |
| `materials.rs` | young, shear, hyperelastic c01/c20/c30, yield_stress, hardening_modulus, hardening table stress, JC a/b | pressure |
| | poisson, strains, JC n/c/m, taylor_quinney, damage d1–d5 | dimensionless |
| | JC reference_strain_rate | frequency (1/s) |
| | JC room/melt temperature, thermal reference_temperature | temperature (absolute) |
| | conductivity, specific_heat, expansion | conductivity, specific heat, expansion (1/ΔT) |
| | density; fluid viscosity (doc: "1e-9 MPa·s for water") | density; dynamic viscosity |
| `sections.rs` | profile dims, shell thickness/offset, 2D thickness, area | length, area |
| `loads.rs` | force, moment, traction, pressure, gravity, body force, centrifugal ω/point, surface flux, film, sink, heat source power, temperature | force, moment, pressure, acceleration, force/volume, angular velocity, heat flux, HTC, temperature, heat generation |
| `bcs.rs` | displacement `[Option<f64>; 6]` (3 lengths + 3 rotations rad), temperature | length, angle, temperature |
| `steps.rs` | frequency min/max_hz, buckle shift, end_time/increments/intervals, initial velocity, rigid wall point, hourglass/bulk viscosity coefs, heat theta/initial_temperature, imperfection scale, displacement control | frequency, time, velocity, length, temperature (coefs dimensionless) |
| `contact.rs` | μ static/kinetic, penalty_scale; **decay (time per length)** | dimensionless; s/m |
| `discrete.rs` | point mass, inertia, spring stiffness, dashpot coefficient | mass, mass·length², stiffness, damping |
| `fluids.rs` | wall/inlet velocity, pressure, thermal temperature/flux, initial temperature, gravity, β, domain min/max/spacing/tolerance, a conducting solid's power, a thermal contact's resistance or layer thickness and conductivity | velocity, pressure, temperature, heat flux, acceleration, expansion, length, power, thermal contact resistance, thermal conductivity |

**Risky:** `units` missing in a file means `unspecified`, so a migration cannot know the system of
old files without one; the step cache key includes `units`.

## Commands and forms (`crates/engine/src/cmd`, `schema.rs`, `cmd/forms.rs`)
Definition commands (`model.{material,section,load,boundary_condition,…}.{create,edit}`) deserialize
JSON straight into the model structs (`cmd/define.rs`), so every number is a bare `f64`. Each form
field carries a `Dim` (`schema.rs:299-418`); `Dim::unit(system)` only builds a label. Conversion
goes where the form value is read: `Field` → `Quantity` and `ValueInput::to_si`.

**Risky:**
- `Dim::Millimetre` is a fixed "mm" for meshing and geometry commands (tet/hex size, tolerances,
  datums, hex project, defeature, mesh plan, geometry check), while the same quantities in mesh
  editing commands (node add/move/relax, seeds, size regions, block, mirror, planar, stitch) are
  `Dim::Length` and follow the system. One physical quantity, two assumptions.
  *Resolved:* `Dim::Millimetre` is gone; every length is `Dim::Length` and follows the display
  system, the geometry frame converting at the command boundary (docs/units.md). Journals from
  before SI storage keep their millimetres through a frozen table.
- Missing dimensions: load heat_source power (`forms.rs:561`, really W/m³), frequency shift
  (`forms.rs:2024`, really 1/s²), contact decay (`forms.rs:1015`, s/m), JC reference strain
  rate (`forms.rs:137`, 1/s).
- `Dim::Temperature` has no unit; `Dim::Expansion` is fixed "1/K" while temperatures are unlabelled.
- JSON Schema (`schema.rs:1260-1365`, list_commands, `fieldcraft-cli commands --json`, MCP) emits
  only `x-dimension` names and none for vectors and number_or_auto (forces, positions, gravity,
  velocities lose their dimension).
- Help/defaults assume mm-t-s: `forms.rs:520` (gravity −9810), `define.rs:98` (210000, 7.85e-9),
  `define.rs:291`; samples `engine/src/sample.rs:362` `GRAVITY = 9810.0` and aluminium
  70000 / 2.7e-9.
- `model.units` (`cmd/model.rs:39-54`) only sets the label; it warns about mm geometry only while a
  design is open.

## UI display (`crates/ui-egui`)
Forms (`form.rs`) label through `Field::unit_label/help_text`; descriptions (`describe.rs`),
replies (`reply.rs`), results panel and probes (`results_panel.rs`, `lib.rs` `probe_rows`) through
`Quantity::unit(system)`. Display conversion goes here (`Quantity::format` with the display
system). **Risky:** literal "mm", "mm²", "mm³" in `reply.rs:398, 708, 720, 723, 770` and
`feature_graph.rs:52-257`.

## Results (`crates/engine/src/results.rs`, `flow.rs`, `cmd/xy.rs`, `cmd/results.rs`)
`Field::quantity()` maps fields to the engine's `Quantity` (Length, Stress, Force, Moment, Mass,
Time, Velocity, Temperature, Energy, HeatFlow, HeatFlux); `Quantity::unit` builds the label, with
per-system strings for energy/power/flux. Reactions and probe replies are numbers only; contour
replies carry `"unit"`. **Risky:** temperature is always `"°C or K"` (also in in-lbf-s);
`ContactGap` (a length) maps to `None`. Conversion: results stay SI; legends and replies format.

## Files
| Format | Today | Conversion |
|---|---|---|
| `.fieldcraft` (`crates/io`, `FORMAT = "fieldcraft/1"`) | `model.json` numbers as entered plus `"units"`; `mesh.bin` raw coordinates | scale to SI once on open by the declared system. **Risky:** FieldCraft has no pending-format mechanism (only the `FORMAT` const and its refusal); SI storage needs a format version step. Files with `unspecified` need a rule (assume mm-t-s-N-MPa, since geometry is mm). |
| Abaqus / Nastran / LS-DYNA (`crates/decks`, write only) | unit label in a comment only | write in a chosen consistent system (`DisplaySystem` factors). **Risky:** no Nastran `PARAM,WTMASS`, no LS-DYNA `*CONTROL_UNITS`; there are no deck readers yet. |
| SolidCraft geometry (`crates/geometry`) | "lengths are millimetres" (`lib.rs:16, 96`) | ×1e-3 once at import |
| 3MF (`geometry/src/tagged_mesh.rs:151-161`) | **the one real conversion outside the table**: micron/mm/cm/inch/foot/metre → mm; writes `unit="millimeter"` | move factors to the table, convert to m |
| STEP | no import exists | — |
| STL export (`io/src/stl.rs`) | raw coordinates | write in the display length (STL has no unit field) |
| VTU/PVD (`io/src/vtk.rs`) | no unit metadata | add units in field data |
| CSV (`cmd/results.rs:1023-1075`) | bare headers; reply has `"units"` | units in headers (`fx [N]`) |

## CLI and MCP
`crates/mcp/src/server.rs:26-27` instructions: "unit-consistent and never converts … use
mm-t-s-N-MPa"; `tools.rs:105` example `model.units`; screenshot distance "in model units". The
CLI adds no units itself (commands schema as above). Replies switch to values with units.

## Solvers
- LBM (`crates/lbm/src/setup.rs`): lattice units Δx = Δt = ρ₀ = 1; ν, velocities, accelerations,
  pressure, force, heat scale through h, dt, ρ, c. Self-consistent in any system, no fixed
  factors; unchanged by SI storage. Temperatures pass unscaled (only differences matter).
- Implicit, explicit, fem: no unit literals; thermal strain uses α(T − T_ref), Johnson–Cook uses
  room/melt in one scale: correct in kelvin. No 273.15 or Stefan–Boltzmann anywhere (radiation
  would need absolute kelvin, which SI storage gives).

## Literal factors
`cargo xtask unit-factors` lists float literals equal to common factors outside the table: 548
today (mesh 187, engine 175, geometry 34, render 33, fem 25, …). Almost all are tolerances,
limits or bytes→MB; the real conversions are the 3MF table and the mm-t-s sample constants.

## Implicit "geometry in mm"
Mesh coordinates come from mm geometry whatever `model.units` says: an m-kg-s or in-lbf-s model
reads mm numbers as metres or inches; the only guard is the `model.units` warning. SI storage
fixes this at the import boundary. `MAX_COORD` is "1e9 model units" in code, "±10⁹ mm" in
`docs/formats/tagged-mesh.md:45`.
