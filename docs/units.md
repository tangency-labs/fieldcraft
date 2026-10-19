# Units

FieldCraft stores every quantity in SI and converts only where a number crosses into or out of
the program, through one table: `crates/model/src/units/` (#579).

## The rule

**Store SI, convert only at boundaries through `crates/model/src/units`.** The model, the result
types and every solver work in metres, kilograms, seconds, newtons, pascals, kelvin and radians.
A conversion factor written anywhere else is a bug (`cargo xtask unit-factors` lists the literals
that look like one; most of what it lists today are tolerances, see the end of this page).

## The display system

Numbers are shown and typed in a *display system*: a preset with optional per-quantity overrides
(`fieldcraft_model::DisplaySystem`).

| Preset | Length | Mass | Time | Force | Stress | Temperature |
|---|---|---|---|---|---|---|
| `mm-t-s-N-MPa` (default) | mm | t | s | N | MPa | °C |
| `m-kg-s-N-Pa` (SI) | m | kg | s | N | Pa | K |
| `mm-kg-ms-kN-GPa` | mm | kg | ms | kN | GPa | °C |
| `mm-g-ms-N-MPa` | mm | g | ms | N | MPa | °C |
| `in-lbf-s-psi` | in | lbf·s²/in | s | lbf | psi | °F |

Every other quantity's unit follows from the preset (`Quantity` lists them; energy in mJ for
mm-t-s, J for SI, …). Angles are radians in every preset; an override such as
`{"preset": "mm-t-s-N-MPa", "overrides": {"angle": "deg", "energy": "J"}}` changes one quantity.

Which system applies, in order:

1. the model's own (`model.units {system, overrides}`, saved in the model as `display_units`),
2. the user's setting (Help › Settings › Units, `prefs.json`; `ui.prompt settings set {units}`),
3. `mm-t-s-N-MPa` with °C.

A model migrated from format 2 keeps the system its file declared (`unspecified` displays as
mm-t-s-N-MPa), so it shows the numbers it showed before. The defaults are pending the owner's
confirmation.

## Boundaries: where conversion happens

- **Command parameters** (`crates/engine/src/units.rs`): a bare number is in the display unit of
  its field's quantity (the form schema knows it, `schema::Dim`); `"210 GPa"` or
  `{"value": 210, "unit": "GPa"}` carries its own unit. The parameters are converted field by
  field into the system the command works in ([`Frame`](../crates/engine/src/units.rs)): SI for
  the model and the solvers, millimetres for the geometry commands and the meshers. A value of the
  wrong dimension or an unknown unit is refused, naming the field. Fields of a fixed unit
  (`Dim::Millimetre` for datum coordinates and mesh sizes typed against the geometry,
  `Dim::Degree`) take that unit always.
- **Replies**: each command's dimensional reply fields are listed in `units::REPLIES` and
  converted into the display system, with a `units` map naming the unit of each. A fluid
  step's steady-state monitor history follows its monitors, so `solve.fluid` converts it per
  monitor (velocity, temperature change, mass flow, force, heat, mean temperature) and names
  each one's unit in `stop.units`.
- **Forms, panels, legends, probes, XY plots**: the UI reads the model and results in SI and
  shows them through the display system (`describe.rs`, `results.rs`, `results_panel.rs`); an
  edit form is filled with the definition converted to the display system, so OK without a
  change gives back the same model.
- **Viewport**: the scene is drawn in metres; a click that places something (a relax sphere, a
  node drag) writes the display units into the form. `ui.camera`'s distance and target and
  `ui.worldToScreen`'s points are in the display units.
- **Exports**: CSV headers name their unit (`fx [N]`); VTU files carry a
  `<!-- units: … -->` comment and their values are in the display units; an STL file (which
  has no unit of its own) is in the display length, and `mesh.export.stl` replies with it.
- **Solver decks**: written in one consistent system (`deck.write {units}`, default the display
  preset), stated in the header; Nastran gets `PARAM,WTMASS,1.0`, LS-DYNA `*CONTROL_UNITS`.
  Written in a model's old system, a deck is byte-identical to before apart from the unit lines.
- **Files**: `.fieldcraft` format 3 stores SI (see below).
- **Geometry**: SolidCraft designs and imported surfaces are in millimetres; their
  tessellations become metres where they meet the model (`units::geometry_point`), and the
  meshers run in millimetres inside the geometry frame. 3MF files are read in their declared unit
  through the table.

## Where conversion never happens

The model (`fieldcraft-model`), the solvers (`fieldcraft-solver`, `-fem`, `-explicit`, `-lbm`,
`-lbm-gpu`) and the result types never convert. Solver tolerances are relative (to the model's
size, the load, the stiffness) so a model in metres converges as it did in millimetres; the
camera's clipping and minimum sizes are relative or in metres (`render::MIN_SIZE`).

Relative is not enough where two dimensions meet in one norm. A nonlinear step's convergence
test measures each kind of freedom against its own scale (docs/fem.md, *Convergence norms*):
residual forces against the largest force, residual moments against the largest moment (at
least the force scale times the model's size L), translations of a correction against the
increment's translation, rotations in radians against the increment's rotation (at least the
translation over L). One norm over newtons and newton-metres, or metres and radians, would
change with the unit system and the part's size, since a moment scales as a force times a
length and a rotation not at all. The arc length of an arc-length step weighs rotations by L² for
the same reason.

## Exactness

A typed value is scaled on its shortest decimal form (`units/decimal.rs`): `3.97` mm is stored as
the double nearest 0.00397 m and shown again as `3.97`, for any value of up to 15 significant
digits. Longer values (computed coordinates, results) are scaled by one correctly rounded
operation by an exact power of ten. Offsets (°C, °F) are added in decimal: 20 °C is 293.15 K
exactly, and an absolute °F, (F + 459.67) × 5/9, is worked out in decimal where the result is a
short decimal (68 °F is 293.15 K). Factors that are not powers of ten (inch, pound) are applied
in binary.

A command on the geometry runs on the model expressed in millimetres: the millimetre image of the
stored metres, a function of the stored model alone. So a session, its undo, a saved and reopened
model, a journal replay and a crash replay all give a mesher the same input, and the same result.
The meshers alone (tet, hex, cells, `mesh.plan.apply`) give what FieldCraft 0.7 gave, and a
stored model is unchanged by the move to SI. But a mesher-computed coordinate does not always
survive mm → m → mm: 6–18% of them come back one ulp off. A relaxation of such a mesh (guarded or
strong) amplifies that, as a 1-ulp nudge does on 0.7 itself: on the samples, relaxed nodes land up
to about 2 mm from where 0.7 put them (the crank: 1.87 mm guarded, 1.09 mm strong). Results after
a relaxation are not those of 0.7, bit for bit or to the micrometre.

## Temperatures

Absolute temperatures are kelvin; a temperature *difference* (a film coefficient's per-degree, an
expansion coefficient) has no offset. Left out, an initial or reference temperature is 0 °C
(273.15 K), and a plane-strain thickness 1 mm: the defaults of the default display system. The
lattice-Boltzmann thermal solver runs its temperature lattice relative to the zero of the
temperature unit the model is shown in (its own display units, else °C): 0 °C for the
millimetre systems, 0 K for SI, 0 °F for in-lbf-s-psi. Its walls are not invariant under a shift
of the temperature scale, and before SI storage it ran on the numbers as typed, so every model
gives the results it gave.

## The file format and old files

Format 3 (`fieldcraft/3`, docs/formats/fieldcraft.md) stores SI: `model.json` in SI, `mesh.bin`
as `MCMESH02` with coordinates in metres. `migrate::step` 2 → 3 (`crates/io/src/migrate/v2_to_v3.rs`)
reads the declared `units` of the file and scales every dimensional field with its own frozen
field → quantity map and its own copy of the MCMESH01 codec, so a later change to the model types
cannot change how an old file migrates; a number it does not know is refused. Datums, mesh
plans, seeds and size controls migrate from millimetres in a file built on a design or 3MF meshes
(the notes below). The journal's numbers are not
rewritten: the step gives it a header naming the system they were typed in, and replay reads
them through it, so replaying the journal gives the migrated model (but for the last bits of
coordinates a mesher computed on another mesher's mesh, above). Where such a journal's
commands left out a field whose format 2 default meant a number in that system (a reference
temperature of 0, a plane-strain thickness of 1, a Johnson–Cook reference strain rate of 1, an
initial temperature of 0), replay writes that default in, from a table frozen with the step
(`migrate::format_2_defaults`). In a model built on a design or 3MF meshes the geometry and mesher
commands' numbers are read as millimetres, as they were then; one that acted on a block typed in
another unit (beside a mesh of the design) cannot be replayed exactly and is refused, naming the
block.

Notes on migrating older files:

- An m-kg-s-N-Pa file's temperatures are read as kelvin, the system's own unit (the release
  that wrote it labelled them "°C or K"; its solvers depended only on differences, and its
  fluid solver's heat lattice keeps counting from 0 K, as it did, through the domain's
  `temperature_zero`). An in-lbf-s-psi file's are °F, the millimetre systems' °C.
- A model that declared metres or inches but holds a mesh made from a SolidCraft design or a
  3MF file: format 2 never scaled such a mesh, so its numbers are millimetres and the part
  solved was 1000× (25.4×) the design. Each body is migrated by its origin, traced through the
  journal: a mesh of the geometry from millimetres (so it overlays the design again and agrees
  with the datums, mesh plans, seeds and size controls), a block (`mesh.block`, tagged
  elements, and what was mirrored or extruded from them) from the declared unit, and a
  relaxation field with the mesh it was made on (its corners, else its bodies, else where it
  lies). A body the journal does not explain is read as the geometry's, and the note says so.
  Opening such a file says so (in the open reply, the status line and the document's
  migrated-from note): results will differ from the old ones. The first save keeps the old file
  beside the new one.
- A node shared by a mesh of the geometry and a block typed in metres or inches (a block joined
  onto a face of the design's mesh) has no consistent scale: such a file is refused, naming the
  bodies, and is not changed; it still opens in FieldCraft 0.7.x.
- Datums, mesh plans, seeds and size controls were typed against the geometry: millimetres in a
  file built on a design or 3MF meshes, the declared unit in one that was not.

## Literal factors

`cargo xtask unit-factors` reports float literals equal to common conversion factors outside the
table. What remains are tolerances and limits (1e-12, 1e-3 relative tolerances, 1e12 caps), the
frozen factors of the 2 → 3 migration step (which must not follow the table), and the
mm-t-s sample scripts' typed values (`sample::GRAVITY`, 9810 mm/s²): typed input, converted like
any other. None converts a stored value.
