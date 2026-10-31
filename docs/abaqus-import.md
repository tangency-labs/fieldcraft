# Read an Abaqus input deck

`deck.read` imports flat model data for inspection and a new FieldCraft analysis:

```json
{"command":"deck.read","params":{"format":"abaqus","path":"part.inp","units":"mm-t-s-N-MPa"}}
```

Choose the input's consistent unit system. A deck has no intrinsic units; FieldCraft does not
infer them from its numbers, title or comments. Without `units`, the command uses the current
display **preset**, ignoring individual display overrides. The reply states that preset.
Coordinates, elastic moduli, density and section dimensions are converted through the central
quantity table and stored in SI. The pure library `abaqus::read` defaults explicitly to SI;
`abaqus::read_in` takes a chosen preset, and `abaqus::read_limited` also takes the caller's
node, element and set-member caps.

The file is read and its new model constructed before the current session changes. A read,
syntax or validation error keeps the existing model, results, path, undo and unsaved work.
An empty, heading-only or unknown-only file cannot replace the current project. Success
replaces its mesh/definitions, clears old geometry/results/history/path, and leaves an unsaved
model. Save it as a new `.fieldcraft` project.

In the desktop app, open **Read Abaqus Deck** from the command search or File commands,
choose the path and consistent units, then confirm. The reader runs as a cancellable background
job. Replacing unsaved work asks Save / Don't Save / Cancel; an untouched imported model also
needs saving even though its old command history has been cleared. The reply panel lists the
unread keywords and their lines. A programmatic command never opens a dialog.

## Reproduce an original import and solve

The reusable [beam recipe](../tools/abaqus-import/beam.json) starts a new document, constructs
an original 80 x 10 x 10 mm steel beam (81 nodes, 32 hex8 elements with incompatible modes),
clamps its root and applies a total 100 N downward tip force. It solves that source model,
writes its own Abaqus deck, imports the model data with explicit mm-t-s-N-MPa units, then
recreates the intended constraint/load/static step and solves again. Its final project includes
the newly solved results. The generated deck's analysis history is reported as unread.

Run the built CLI from an empty work directory of your choice, using the recipe's absolute path:

```text
fieldcraft-cli run /absolute/path/to/Fieldcraft/tools/abaqus-import/beam.json
```

Outputs `abaqus-beam.inp` and `abaqus-imported.fieldcraft` are relative to that work directory;
change their `path` fields to choose other locations. Neither generated output is source data.
The same command list can be replayed through the app's normal command interface. For an
interactive demonstration, pause before `deck.read`, submit its path/units in Read Abaqus Deck,
review the unread-history report, then continue the remaining commands. The recipe replaces
the current document; save existing work before running it in a session you are using.

## What is imported

- `NODE`, with optional node set; labels and coordinates.
- `ELEMENT`, with type, connectivity and optional element set. Long connectivity may continue.
- `NSET` and `ELSET`: label lists, other set names and bounded `GENERATE` ranges. A node set may
  also name an element set. References resolve without changing their labels.
- `SURFACE, TYPE=ELEMENT`: solid faces, shell sides and planar element edges, using an element
  label or element set.
- `MATERIAL`, constant `ELASTIC` (isotropic or engineering constants) and constant `DENSITY`.
- `SOLID SECTION`, `SHELL SECTION` and the writer's `BEAM SECTION` profiles: rectangle, circle,
  pipe, centred symmetric I, and equal-thickness box. Thickness, rod area, shell offset/points
  and beam orientation are retained where representable.

The reader inverts the writer's supported solid, shell, planar, beam and truss element types:
tet4/tet10/wedge6/pyramid5/hex8/hex20; tri3/quad4/quad8 shells; plane stress, plane strain and
axisymmetric triangles/quadrilaterals; B31 and two-/three-node trusses. Encoded reduced,
incompatible-mode and hybrid formulations are assigned to sections. A section mixing full and
reduced hex elements is refused by name instead of being guessed or split.

Names are matched case-insensitively and retain their original spelling. Definitions differing
only in case are refused. An `NSET`, `ELSET` and element `SURFACE` with the same exact
name combine into one FieldCraft set holding nodes, elements and faces, as in the writer. Quoted identifiers must still fit FieldCraft's name rule:
a letter followed by letters, digits, `_` or `-`, up to 64 ASCII characters. Unsupported quoted
spaces/characters and longer names are refused rather than silently renamed. Mesh labels must
fit `1..=99_999_999`.

Comments and headings are ignored. Abaqus sections do not carry a semantic section name;
the importer generates `Section1`, `Section2`, and so on. It does not recover names from comments.
Canonical roundtrip tests compare the supported keyword/data definitions, mesh and physical
section data, not comment text or discarded analysis history. Some writer encodings cannot
recover extra metadata: a reduced tet4, for example, is written as C3D4 without a separate
integration marker.

## Unread data and refusals

The reply lists unread keywords with their original one-based physical lines. `unknown` holds
the first 100 records; `unknown_count` counts all of them. Read the warnings and this list before
creating an analysis. No unknown text is stored in the model or re-emitted on export.

Steps, boundary conditions, loads, contact, coupling and other history/model extensions are not
imported. Every keyword inside an unread step stays reported, even if that keyword is supported
at model level. Nested, unmatched or truncated step scopes are refused. Node-based surfaces are
reported as unread; they are not converted into element surfaces.

`INCLUDE`, parts, assemblies and instances are refused. Global `SYSTEM` records are refused
because later Cartesian nodes could otherwise be misread in local coordinates. `NODE, SYSTEM=R`
is accepted; cylindrical/spherical node coordinates are refused. Alternate input files,
coordinate-system transforms, unsupported physical parameters and temperature-/field-dependent property tables
are not silently interpreted as their default variants. A bare valued parameter such as `TYPE`
is an error; `GENERATE` is a flag and must not have a value.

Recognized material behaviors keep the current `MATERIAL` context even when they are unread.
This includes `HYPERELASTIC`, `PLASTIC`, `RATE DEPENDENT`, `SPECIFIC HEAT`, `INELASTIC HEAT
FRACTION`, `CONDUCTIVITY`, `EXPANSION`, `CREEP`, `VISCOELASTIC`, `DAMPING`, `DAMAGE INITIATION`,
`DAMAGE EVOLUTION`, `DEPVAR` and `USER MATERIAL`. Constant elasticity and density before or
after them are imported into that material; every unread option is still reported with its
line number. For example, a writer-generated hyperelastic material retains its density, but
its hyperelastic law is unread and its stiffness is absent. Review and supply the intended
material behavior before solving. No plastic, hyperelastic or thermal law is inferred.

Another model option or an unrecognized keyword ends material context. A later `ELASTIC` or
`DENSITY` then needs a new `MATERIAL` definition, preventing association across unknown scopes.
The rule follows the "Material data definition" chapter of the public *Abaqus 2025 Materials
Guide*.

Element family/formulation is stored on sections, not on mesh elements. Sectionless elements
therefore produce an explicit warning: the mesh can be inspected, but those physical semantics
are not retained until compatible sections are assigned. T2D-only trusses without a planar
continuum cannot retain the model's required planar context and are refused. Mixing T3D trusses
with a planar continuum is refused rather than silently imposing different constraints.

Supported records must remain finite and consistent: missing referenced labels/sets/materials,
overlapping section assignments, invalid connectivity/faces and unrepresentable section shapes
are refused with line numbers. This reader does not promise that an arbitrary imported mesh is
ready to solve. Add the intended constraints, loads and analysis steps through normal commands,
then run Check Model.

## Bounds and automation

Files are checked against the session's read roots and capped at 64 MiB. Nodes, elements and
set and surface memberships are capped by the session's import limits, the same ones the
orphan-mesh import uses, so a host or test with a smaller budget is refused at once. Physical
lines are capped at 256 bytes; logical keyword continuation, total record count, property-value collection,
set expansion/work and nesting are bounded separately. A short malicious range cannot allocate
an enormous set. Input must be UTF-8 with ASCII keyword/data records. The parser has no filesystem
access and builds for wasm; browser files still pass through the existing host file store.

The format implementation is independently written from Dassault's public keyword documentation
and FieldCraft's original writer/model mappings. Tests use original handwritten decks and
FieldCraft-generated decks, never a commercial CAE program's output or a vendor example model.

## Format references

The syntax and supported options follow Dassault Systemes' public documentation, cited by title
and version: the "Input syntax rules" chapter of the *Abaqus 2017 Analysis User's Guide*, and
the `*NODE`, `*ELEMENT`, `*NSET`, `*SURFACE`, `*ELASTIC` and `*SOLID SECTION` entries of the
*Abaqus 2017 Keywords Reference Guide* and its `*BEAM SECTION` entry in the 2025 edition.
