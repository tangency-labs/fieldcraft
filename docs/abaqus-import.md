# Read an Abaqus input deck

`deck.read` imports model data, including part instances, for inspection and a new FieldCraft analysis:

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
model. Save it as a new `.fieldcraft` project. Original INP files are retained separately in
the project as `deck-source.json`, protected by the package size and SHA-256 checks.
Reopening validates paths and bounds but never reads the original input files from disk.

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
solves the imported static step directly, without recreating constraints or loads. Its final project includes
the newly solved results. The source bundle and solved model are both saved.

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
- Constant isotropic `MASS` and grounded `SPRING1` elements with `MASS`/`SPRING` properties.
  These become native point masses and springs named `InpDiscreteElement<label>` with explicit,
  collision-free node target sets. Translational and rotational spring stiffness use their own
  unit quantities. Discrete ELSETs resolve properties and `NSET, ELSET=...`; discrete labels are
  excluded from native mesh element sets, with a warning. Their original labels and memberships
  remain in the retained source. Translated instances remap labels and property set references;
  springs in rotated instances are explicitly refused. Anisotropic/damped masses, nonlinear or
  frequency-dependent springs, SPRING2/A, dashpots and rotary inertia are not mapped yet.
  Gravity history with imported point masses is refused until discrete load membership is mapped.
- `NSET` and `ELSET`: label lists, other set names and bounded `GENERATE` ranges. A node set may
  also name an element set. Flat model labels are unchanged; instance labels are mapped to
  deterministic global labels and reported with their original identity.
- `SURFACE, TYPE=ELEMENT`: solid faces, shell sides and planar element edges, using an element
  label or element set.
- `MATERIAL`, constant `ELASTIC` (isotropic or engineering constants) and constant `DENSITY`.
- Constant rectangular/cylindrical `ORIENTATION` records and `SHELL SECTION, COMPOSITE`;
  see [composite deck mappings and refusals](composite-decks.md).
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

Comments and headings do not affect the interpreted model but remain in the source bundle. Abaqus sections do not carry a semantic section name;
the importer generates `Section1`, `Section2`, and so on. It does not recover names from comments.
Canonical roundtrip tests compare the supported keyword/data definitions, mesh and physical
section data, not comment text or uninterpreted analysis history. Source-replay tests separately require exact bytes. Some writer encodings cannot
recover extra metadata: a reduced tet4, for example, is written as C3D4 without a separate
integration marker.

## Parts, instances and original source

`PART`/`END PART`, one `ASSEMBLY`/`END ASSEMBLY`, and its `INSTANCE`/`END INSTANCE`
scopes are expanded before model construction. Mesh data may live in the part or the instance.
Local node/element labels are allocated globally without collisions. Local sets, surfaces and
orientations receive deterministic `InpLocal…` names; the `instances` reply carries the original
instance/part names and label/name maps. `INSTANCE=` sets, `GENERATE`, `Instance.Set`,
`Instance.Label` and assembly-qualified references are resolved. Material names stay global.
Uninstanced parts remain in the source bundle and are listed in a warning. Rotated orthotropic
solids must supply an explicit part orientation; implicit material axes are not inferred.
Translation precedes axis rotation. Coordinates, constant orientation points and beam directions
follow that transform. Missing parts/references, duplicate labels/instances and malformed scope
endings are errors. More general coordinate systems and imported-result instances are refused.

`INCLUDE, INPUT=…` resolves relative to the main INP directory, used as the submission directory,
including when the directive occurs in a nested file. Data-only fragments are supported. Absolute
paths, `.`/`..` components, backslash/drive paths, encrypted inputs and symlink escapes are refused.
The directory boundary is enforced even in an unconfined interactive session, in addition to
any automation read root. Cycles, depth over 32, more than 1024 files, source/expanded text over
64 MiB and excessive expansion work are refused. Errors identify the source file and line.

To recover the imported source without normalizing its spelling or interpreting its history:

```json
{"command":"deck.source.write","params":{"path":"copy/original.inp"}}
```

The destination parent must exist. The command writes the retained root file and its dependencies
at their original relative paths, preserving CRLF/LF, comments, unknown options and missing final
newlines. It never rereads the inputs. It checks all destinations before writing, writes dependencies
first and the root last; individual writes are atomic, the whole bundle is not a transaction.
A root filename colliding with a dependency is refused. The reply explicitly states that model
edits were not exported. Use `deck.write` for an edited model. Starting/opening another project
replaces the retained source with that project’s saved source, or clears it if absent.
Save the project or export the source before discarding an unsaved import.

These syntax mappings follow the public [Assembly Definition](https://abaqus.uclouvain.be/English/SIMACAEMODRefMap/simamod-c-partassy.htm),
[INSTANCE](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-instance.htm),
[INCLUDE](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-include.htm) and
[Input Syntax Rules](https://abaqus.uclouvain.be/English/SIMACAEMODRefMap/simamod-c-inputsyntax.htm).
Fixtures are original handwritten inputs, not manual examples.

## Unread data and refusals

The reply lists unread keywords with their original one-based physical lines. `unknown` holds
the first 100 records, with file/line and up to 4096 characters of original text each;
`original_truncated` indicates a preview limit, while the source bundle retains the full text.
`unknown_count` counts all unread occurrences. Read the warnings and this list before
creating an analysis. Unknown text is retained in the original source bundle, separate from the
model; `deck.write` exports only the interpreted model. `deck.source.write` replays originals.

The command maps linear `STATIC`, `FREQUENCY` and a first-base-state `BUCKLE` to native
steps. Structural `BOUNDARY` accepts numeric DOFs, ENCASTRE, PINNED, X/Y/ZSYMM and X/Y/ZASYMM;
`CLOAD` maps per-node forces/moments, `DLOAD` maps element pressures and gravity, and
`DSLOAD` maps element-surface pressure. `OP=NEW` clears the named family once per step;
`OP=MOD` modifies the same original target/load type and retains other definitions. Different
node/set identities remain additive, even where their expanded nodes overlap.

Relative TABULAR/SMOOTH STEP amplitudes are evaluated at the linear static equilibrium end.
STEP TIME values freeze at the preceding step end when carried forward; TOTAL TIME curves
continue along accumulated general-step time. This does not validate any transient trajectory.
The full curve remains in the original source. Native frequency steps use homogeneous support
DOFs and no applied loads; native modal normalization/output defaults apply.

Contact, unsupported coupling variants, unsupported material laws, nonlinear/dynamic steps and unsupported load
variants are not silently omitted from runnable analyses. Unhandled model physics blocks analysis
mapping; an unsupported step stops that and later mappings. `history_imported` says whether any
step was imported; `analysis_complete` says whether the entire supported history was mapped.
It is not a solver-success flag. Output requests remain in the unread report; native procedure
default outputs are used. Every keyword inside an unread step stays reported. Nested, unmatched or truncated step scopes are refused. Node-based surfaces used by supported couplings become native node sets; other node-based
surfaces remain unread. They are not converted into element surfaces.

`PARAMETER` and global `SYSTEM` records are refused
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
then run Check Model. The pure `abaqus::read_limited` retains its model-data-only contract;
`abaqus::read::read_analysis_limited` enables the analysis mapping used by `deck.read`.

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


Direct harmonic import also maps constant material/global Rayleigh and structural damping,
explicit single-frequency rows (`upper=0`) and real concentrated/distributed loads. Frequency
ranges with bias, modal/subspace harmonic methods, imaginary loads, prescribed harmonic motion
and prestressed harmonic base states remain reported as unsupported. The `imported_steps` reply
lists the steps actually created; `analysis_complete=false` still means the entire source is
not mapped. A later matrix-export step is a separate native matrix-generation procedure
(see below); `solve.harmonic` executes only the selected harmonic step, while
`solve.matrices` assembles and optionally exports the operators without repeating
the harmonic sweep.

Native solver limits still apply after a complete physical mapping: B31 harmonic stress
recovery is currently unsupported and is reported by `solve.harmonic`; the importer does not
remove requested stress to force a solve. The direct-harmonic solve regression uses an original
solid beam, and a separate B31 regression requires the explicit unsupported-stress error.


### Modal random response and operator export

`deck.read` maps a preceding unwindowed `FREQUENCY` extraction followed by
`RANDOM RESPONSE` into a native modal random step. Current scope is one real,
constant `PSD-DEFINITION,TYPE=FORCE` multiplier, one fully correlated load case
with unit correlation weight, homogeneous nodal constraints, and mode-number
`MODAL DAMPING`. Distributed pressure is case 1; concentrated forces/moments
require an explicit other case. Spatial load magnitudes remain physical; the
scalar spectrum multiplies the coherent response power once. Interval sampling
uses the extracted natural frequencies, the requested points **per interval**,
bias and linear/logarithmic spacing. `solve.random {"step":"Random","nodes":[1]}`
executes this model definition. It returns SI component PSD and trapezoidal RMS
in its reply. Save that reply separately: native project persistence retains the
random definition and input source, not random response frames. The default nodal stress components
use native averaged harmonic recovery; the aggregate stress RMS is a component
norm, not von Mises RMS. Optional `elements` selects unaveraged solid recovery as
described below. A small projected
residual does not certify a truncated modal basis. Nonconstant spectra, multiple
or complex correlations, base excitation and preloaded modal histories remain
explicitly unsupported by this import path.

`MATRIX GENERATE` maps stiffness, consistent mass, viscous and structural damping
to a native matrix-generation step. `MATRIX OUTPUT,FORMAT=MATRIX INPUT` is retained;
`solve.matrices {"step":"Matrices","path":"operators.json"}` exports a manifest
and sibling five-column operator files. Each row identifies native model node/DOF
pairs (flattened global labels after instance mapping), with a symmetric lower triangle by node/DOF labels; the manifest separately records independent-equation order. Operators are
always SI (rotations in radians), with `K - omega² M + i(omega C + G)` convention;
material Rayleigh damping contributes to C and loss factors to G. Homogeneous
supports and native MPC transformations are eliminated. The manifest gives the
independent coordinate ordering, units and warnings. No response calculation or
factorization is needed. With no path the command assembles and reports metadata
only; `solve.run` does not silently create files. Loads, initial fields, preload,
nonlinear material behavior and nondefault matrix options are refused. Files are
atomic individually, not a transaction. A native project saves the procedure;
operator files are separate artifacts. Normalized deck writers currently refuse
these procedures; `deck.source.write` replays their original text exactly.

The SDOF random test compares every sampled displacement PSD with the analytic
transfer function, then reopens and executes the saved definition. The axial-bar
matrix test checks K, M, C and G numerically and after reopening. Source replay
checks byte equality only; it does not establish solver equivalence for unknown
keywords or output requests. Syntax sources: the public Abaqus keyword references
for RANDOM RESPONSE, PSD-DEFINITION, CORRELATION, MODAL DAMPING, MATRIX GENERATE and
MATRIX OUTPUT, and its random-response analysis guide. No commercial outputs are
used in implementation or these checks.


### Element stress recovery and truncation checks

`solve.random {"step":"Random","nodes":[1],"elements":[10,11]}` also returns
`element_stress` for each selected solid, both at its own extrapolated nodes and
at its integration points. The same option works with stored harmonic frames and
an explicit `input_psd` table. Recovery uses the solved complex displacements,
the element constitutive law and the existing quadrature/extrapolation kernels;
it never reconstructs an element field from averaged nodal stresses. Each location
includes the element label, local point index, global physical coordinates and,
where applicable, node label or natural integration coordinates. Match quadrature
points by coordinates, not by assuming another program's point numbering.

Component ordering is `[xx, yy, zz, xy, xz, yz]` in global Cartesian axes. The
reply contains SI component PSD and trapezoidal component RMS; optional
`element_phasors:true` includes real and imaginary stresses in Pa for the
deterministic coherent spatial-load transfer, before PSD multiplication; these
are not instantaneous random stresses. Coordinates
remain metres even when the display system differs. There is no claimed von
Mises RMS. Only solid elements are supported by this optional recovery; shell,
beam, rod and 2D selections fail explicitly. Selection is bounded to 1024 unique
elements and one million location-frame samples. Harmonic project results retain
the displacements needed to repeat recovery after reopening; random replies must
still be saved separately.

Original `OUTPUT` selection/frequency requests remain retained source, not native
output-selection definitions. In particular, the caller must explicitly select
the requested element labels. `POSITION=NODES` denotes each element's unaveraged
extrapolated values, whereas integration-point requests use the other location
rows. A successful source replay or `analysis_complete` does not establish that
all original output requests were executed. See the public
[*ELEMENT OUTPUT reference](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-elementoutput.htm).

The full-equation modal residual is computed on all independent equations after
constraint elimination, including the lifted modal damping force. Its relative
value is the maximum complex residual magnitude divided by the maximum complex
load magnitude. `residual_norm_si` separately reports the residual L2 norm.
The projected residual checks only the equations in the retained modal subspace,
relative to their projected loads. A small projected residual demonstrates that
the reduced equations were solved, not that omitted modes make negligible
contributions to displacement or stress. Neither number is a response-error
percentage. Compare responses at representative frequencies and the static limit;
keep the spatial load, PSD and damping assumptions explicit. Truncating through a
repeated-frequency cluster can also make a mode-number damping boundary ambiguous.
Do not silently extend a damping range or replace an undamped pole with an arbitrary
damping value to obtain a result.

Validation separates exact source replay, solved-project reopening, constitutive
recovery and response accuracy. The solid harmonic regression checks a uniform
axial stress, component PSD/RMS and physical/natural recovery positions, then
reopens the solved project and repeats recovery. A manufactured affine field with
two distinct moduli at shared nodes verifies all six components and the absence
of cross-element averaging. These tests do not certify modal truncation for a
particular engineering model.

Named symmetry and antisymmetry conditions expand to global structural DOFs; nodal
transforms remain unsupported and block analysis mapping. The command regression
checks all eight named conditions against explicit DOF lists through project
save/open, exact source replay, and normalized deck write/read. Reference:
[Boundary conditions](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEPRCRefMap/simaprc-c-boundary.htm).

Command-level refusal regressions also reopen source bundles containing unsupported `TIE`
defaults, coupling orientations, `PLASTIC` or nodal `TRANSFORM`. They require the original
text to survive and no runnable step to be created. This verifies conservative
retention, not support for those physical definitions. Named antisymmetry is a
linear structural mapping; it does not enable finite-rotation or transformed-node
analysis.


### Bonded surfaces and reference-point couplings

`TIE` maps one secondary/main pair of solid element surfaces to a native `Tie`.
It requires `ADJUST=NO` and an explicit positive `POSITION TOLERANCE`, converted
with the selected length units. Omitted `TYPE` or `TYPE=SURFACE TO SURFACE` is
accepted. FieldCraft enforces its native dual-mortar constraint; this is not a
claim that interpolation coefficients, selection tolerances or nonmatching-mesh
responses equal another solver. Default adjustment/tolerance, `TYPE=NODE TO SURFACE`,
`TIED NSET`, `NO ROTATION`, shell ties and multiple pairs per name are not mapped.

`COUPLING` maps an isolated numeric reference node or singleton node set, a
`SURFACE, TYPE=NODE` with unit node weights, and `KINEMATIC` or uniform
`DISTRIBUTING, COUPLING=CONTINUUM` to native reference-point constraints.
The supported cloud consists only of solid-element nodes, in global axes, with
all six reference DOFs selected (omitted data, blank first DOF, or ranges covering
1–6). Continuum distribution uses equal nodal weights and translation-derived
cloud rotation. The reference point retains its original label and coordinates;
its temporary mesh-node alias is removed. Original reference node sets are
expanded before this removal. Boundary, force and moment groups split reference
points from ordinary mesh nodes, preserving per-node values and step activity.
Source-only reference sets remain recoverable in exact replay; normalized export
writes reference-point targets as numeric labels.

Part/instance expansion rewrites local constraint names, tie surface pairs,
coupling surfaces and reference node/set names before mapping. Invalid references
fail the read transaction. Unsupported parameters/subcards produce named warnings,
retain the source and block runnable history: e.g. `ORIENTATION`, `INFLUENCE RADIUS`,
`KINEMATIC, ALPHA`, partial DOFs, nonuniform weights, `DISTRIBUTING, COUPLING=STRUCTURAL`,
explicit structural rotational coupling, `FREE NODES`, `SINGLE FACET`, `UNIFORM`,
and reference nodes already targeted by `MASS`/`SPRING` (discrete-target retargeting
is not yet mapped).
Some of these may be expressible with further native mapping; they are not claimed
implemented or validated here. Unknown physics never becomes a successful solve
merely because its source is retained.

Original two-cube command tests compare imported constraints to independently built
native commands, including a translated instance and reference-point force/moment
transfer. The axial patch has area 1 m², length 2 m, E=1000 Pa, nu=0 and load 1 N:
its theoretical end displacement is 0.002 m. Tests compare displacement, nodal stress
and total support reaction, then save/open and normalized write/read and solve again.
Exact source replay checks original bytes separately. These small native comparisons
validate mapping and persistence, not commercial parity, arbitrary mesh mortar accuracy,
nonlinear motion, shell coupling or modal truncation.

Syntax references: [TIE](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-tie.htm),
[COUPLING](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-coupling.htm),
[KINEMATIC](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-kinematic.htm),
[DISTRIBUTING](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-distributing.htm).

### Mass semantics report: requirements are not applied settings

`deck.read` returns `mass_semantics` separately from `analysis_complete`. It lists
source element types and ELEMENT-card ELSET hints (not resolved section assignments), a conditional `consistent`/`lumped` choice when
identified, the actual formula/quadrature requirement, applicability conditions,
and native implementation gaps. `applied=false` and `conditions_verified=false`
are explicit: importing a runnable history does not establish a matching source
mass matrix. Source replay after project reopening preserves the original types;
normalized export represents the native model, not a verified vendor default.

The current evidence scope is uniform-property, unconstrained `MATRIX GENERATE,MASS`
small cases. It does not establish frequency, harmonic or random-response defaults.
Do not change an engineering response calculation on that evidence alone. A global
`mass=lumped` switch is insufficient: native solid/plane lumping uses HRZ diagonal
scaling, whereas a row-sum requirement means summing every consistent matrix row.
These coincide for affine uniform simplices and straight uniform two-node rods,
but not generally for distorted hexes/wedges or radius-weighted axisymmetric elements.
A mixed linear/quadratic mesh may need lumped and consistent mass by section/type;
`mixed_choices_require_per_type_or_section` identifies that conflict. It does not
apply or persist such a policy. Quadratic axisymmetric triangles additionally need
an identified integration rule: six-point degree-four quadrature is distinct from
native `TriConical(4)`. Even `consistent` does not encode this difference.

Pyramid distribution (`U-PYR`: C3D5/H) and shell projection (`U-SHELL`: S3/R, S4/R)
remain named unknowns. S8R translational agreement does not establish its rotational
or drilling inertia. Beam equivalence is limited to uniform straight prismatic,
unoffset conditions. No unidentified option is assigned a guessed default. Commercial
matrix originals are not repository fixtures; formula tests and original command
round trips are separate from local, authorized observational evidence.

### Native section mass policy round trips

An explicit native section setting is separate from inferred source defaults.
The writer retains nondefault policies using a FieldCraft-only comment:

```text
** FIELDCRAFT LUMPED MASS, ELSET=Body, STRATEGY=ROW_SUM
```

`STRATEGY` accepts `HRZ`, `ROW_SUM`, or `CONSISTENT`; its target must have a section.
The reader restores `formulation.lumped_mass` and reports it separately as
`native_lumped_mass_policies`, with activation `explicit mass=lumped solve`.
Duplicate/unknown strategies and missing targets are errors. Metadata on a deck
with PART/INSTANCE expansion is refused until scoped policy references are
supported; the native writer emits a flat deck. Nastran and LS-DYNA exports refuse
nondefault policies instead of discarding them. The comment is ignored by Abaqus:
it is not an Abaqus keyword, a portable mass control, or a way to force a vendor
analysis to use FieldCraft's rule. Exact source replay remains byte-preserving.

Model save/open preserves the policy; a native normalized write/read restores it
and reproduces the explicitly lumped frequency solution. Default/explicit consistent
solutions stay identical. The imported mass-evidence report is still conditional
and unapplied: it does not silently select this native policy for arbitrary input.
S8R rotational inertia, CAX6's particular six-point rule and the named unknown
pyramid/shell distributions remain separate gaps.

**Stiffness boundary:** the current `C3D8` mapping uses FieldCraft's standard
full-integration Hex8. It is not established equivalent to a volume-averaged
B-bar formulation. Consequently writing this native element as `C3D8` does not
certify identical stiffness semantics in another program. The distributing
coupling's kinematic/virtual-work identity and the element stiffness discretization
must be checked separately. Selecting a mass policy does not change or resolve
this stiffness boundary. A B-bar implementation requires its own public-formula
validation and mapping work; changing a label alone is insufficient.
