# Read an Abaqus input deck

`deck.read` imports flat model data and the linear analysis history emitted by the writer:

```json
{"command":"deck.read","params":{"format":"abaqus","path":"part.inp","units":"mm-t-s-N-MPa"}}
```

Choose the input's consistent unit system. A deck has no intrinsic units; FieldCraft does not
infer them from its numbers, title or comments. Without `units`, the command uses the current
display **preset**, ignoring individual display overrides. The reply states that preset.
Coordinates, elastic moduli, density and section dimensions are converted through the central
quantity table and stored in SI. The pure library `abaqus::read` defaults explicitly to SI;
`abaqus::read_in` takes a chosen preset, and `abaqus::read_limited` also takes the caller's
node, element and set-member caps. `abaqus::read_analysis_limited` additionally imports the
writer's static, frequency and buckling steps with boundary conditions and concentrated loads;
other history keywords stay in the unread report.

The file is read and its new model constructed before the current session changes. A read,
syntax or validation error keeps the existing model, results, path, undo and unsaved work.
History mapping problems are handled consistently: an unsupported or malformed history block
stays unread with a line-numbered warning, and that step and later steps are not imported. The
flat model data still succeeds unchanged when the base reader can read it; only lexical errors
that the base reader rejects fail the whole read.
An empty, heading-only or unknown-only file cannot replace the current project. Success
replaces its mesh/definitions, clears old geometry/results/history/path, and leaves an unsaved
model. Save it as a new `.fieldcraft` project. The command reply reports the counts of nodes,
elements by kind, sets, materials, sections, steps, boundary conditions and loads alongside the
unit preset, unread keywords and warnings.

In the desktop app, open **Read Abaqus Deck** from the command search or File commands,
choose the path and consistent units, then confirm. The reader runs as a cancellable background
job. Replacing unsaved work asks Save / Don't Save / Cancel; an untouched imported model also
needs saving even though its old command history has been cleared. The reply panel lists the
unread keywords and their lines. A programmatic command never opens a dialog.

## Reproduce an original import and solve

The reusable [beam recipe](../tools/abaqus-import/beam.json) starts a new document, constructs
an original 80 x 10 x 10 mm steel beam (81 nodes, 32 hex8 elements with incompatible modes),
clamps its root and applies a total 100 N downward tip force. It solves that source model,
writes its own Abaqus deck, imports the model and its static history with explicit
mm-t-s-N-MPa units, then solves the imported step. Its final project includes the newly solved
results.

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
- The writer's linear history: `STEP` with `STATIC`, `FREQUENCY` or `BUCKLE`, plus `BOUNDARY`
  and `CLOAD`. Boundary targets may be node labels or node sets; supported fixed, symmetry,
  displacement and concentrated force/moment forms are converted to the corresponding model
  objects.

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
Canonical roundtrip tests compare the supported keyword/data definitions, mesh, physical section
data, and supported linear history, not comment text or unsupported history. Some writer encodings cannot
recover extra metadata: a reduced tet4, for example, is written as C3D4 without a separate
integration marker.

Automatically generated sets, boundary conditions and loads are named `InpAnalysisN`. Input step
names remain the step names, so `*STEP, NAME=Bend` creates a step named `Bend`; step names are
case-insensitive and a name differing only by case from an earlier imported step leaves that step
and all following history unread. Unnamed steps receive `ImportedStepN`, with the number chosen to
avoid every explicit step name in the input (including names in later steps). The generated object
names do not replace the step name. A named node-set target remains its own boundary condition or
load object: rows for that target merge by type and component, while different targets remain
separate and their common-node loads add in the model. Rows addressed by node labels still group
nodes with equal written values or force/moment vectors. A node-label group reuses an input-named
node set when its resolved membership exactly matches, otherwise it reuses another matching set or
creates one; the canonical target spelling is retained in the load identity, so a node-label load
and a set-targeted load remain separate even when they use the same set and vector. The same
node-label definition is reused across later steps. A history row naming a set uses the node
members listed by that set, even when an element set with the same name contributes extra resolved
nodes; a name with only element members is not a node target.

For history state, `OP=NEW` clears a family once and all blocks in that family must use the same
mode; multiple `OP=NEW` blocks therefore accumulate. Concentrated-load `OP=MOD` replaces only
the same target spelling and DOF, while different target spellings add. A node and DOF may be
held through several targets when every prescribed value is zero; introducing a non-zero value
through a second target is unread. Reusing one target modifies it. A type-format boundary
(`ENCASTRE`, `PINNED` or a symmetry label) and a direct-format boundary on the same target and
degree of freedom are different definitions and require `OP=NEW` to change between them; a
type-format change to another type also requires `OP=NEW`. A direct row on a different target is
still governed only by node-overlap and value rules, because that constraint is the same and is
also emitted by the writer. Static state carries forward, frequency and buckling perturbations do not. Frequency
boundary magnitudes remain on the imported frequency step, but are treated as zero only while
checking target overlap; the preceding and following static state keeps its original values. A
buckling step after any general step remains unread, even when that step carried no load. A
`BUCKLE` block without `EIGENSOLVER` is read when its data line holds only the number of modes: the format's default solver for
buckling is not Lanczos, but with that single field the request is the same. An explicitly named
solver other than Lanczos stays unread.

## Unread data and refusals

The reply lists unread keywords with their original one-based physical lines. `unknown` holds
the first 100 records; `unknown_count` counts all of them. Read the warnings and this list before
creating an analysis. No unknown text is stored in the model or re-emitted on export.

Unsupported history, contact, coupling and other model extensions are not imported. Every
unsupported keyword inside a step stays reported, even if that keyword is supported at model
level. A semantic history problem (unsupported mapping, missing target, duplicate/conflicting
entry, or non-finite value) leaves that step and later steps unread with a line-numbered warning.
An unsupported model-level boundary leaves all history unread. Nested, unmatched or truncated
step scopes remain lexical errors. Node-based surfaces are reported as unread; they are not
converted into element surfaces.

`INCLUDE`, parts, assemblies and instances are refused. Global `SYSTEM` records are refused
because later Cartesian nodes could otherwise be misread in local coordinates. `NODE, SYSTEM=R`
is accepted; cylindrical/spherical node coordinates are refused. Alternate input files,
coordinate-system transforms, unsupported physical parameters and temperature-/field-dependent
property tables are not silently interpreted as their default variants. A bare valued parameter
such as `TYPE` is an error; `GENERATE` is a flag and must not have a value.

A model-level `*TRANSFORM` stays in the unread report. Concentrated loads and boundary conditions
at its nodes are given in the nodes' local system, which the reader does not convert. A step with
a `*BOUNDARY` or `*CLOAD` row on one of those nodes is therefore left unread, with the steps after
it; a model-level `*BOUNDARY` on one of them leaves all history unread, and so does a `*TRANSFORM`
whose `NSET` is missing or cannot be resolved. Each case gives a line-numbered warning.

A model-level `*BOUNDARY` is accepted only before the first `*STEP`. If it appears after a
successful step, that block and all following steps remain unread and no boundary is added to the
model.

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
ready to solve. Add the history that was left unread through normal commands, then run Check
Model.

### Information not recoverable from the deck

The writer expands a model `Total` force into a value at each node when it writes the deck. The
reader therefore groups nodes by their written per-node values and does not infer that they came
from one total load. Small floating-point spelling differences can create separate objects: the
beam example's one written load can read back as six loads because values such as
`-6.249999999999999` and `-6.249999999999998` are distinct. They produce the same solve, but the
reader preserves the values instead of guessing a total-load intent.

## Bounds and automation

Files are checked against the session's read roots and capped at 64 MiB. Nodes, elements and
set and surface memberships are capped by the session's import limits, the same ones the
orphan-mesh import uses, so a host or test with a smaller budget is refused at once. Physical
lines are capped at 256 bytes; logical keyword continuation, total record count, property-value collection,
set expansion/work and nesting are bounded separately. A short malicious range cannot allocate
an enormous set. History entries use the same global set-member and expansion-work cap; if one
step would exceed the model's definition limit, that step and later steps remain unread. Input must
be UTF-8 with ASCII keyword/data records. The parser has no filesystem access and builds for wasm;
browser files still pass through the existing host file store.

The format implementation is independently written from Dassault's public keyword documentation
and FieldCraft's original writer/model mappings. Tests use original handwritten decks and
FieldCraft-generated decks, never a commercial CAE program's output or a vendor example model.

## Format references

The syntax and supported options follow Dassault Systemes' public documentation, cited by title
and version: the "Input syntax rules" chapter of the *Abaqus 2017 Analysis User's Guide*, and
the `*NODE`, `*ELEMENT`, `*NSET`, `*SURFACE`, `*ELASTIC` and `*SOLID SECTION` entries of the
*Abaqus 2017 Keywords Reference Guide* and its `*BEAM SECTION` entry in the 2025 edition.
