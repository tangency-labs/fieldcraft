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
an original 80 x 10 x 10 mm steel beam (81 nodes, 32 Hex8 elements with explicit B-bar),
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
- Constant isotropic or global diagonal anisotropic `MASS` and grounded `SPRING1` elements with `MASS`/`SPRING` properties.
  These become native point masses and springs named `InpDiscreteElement<label>` with explicit,
  collision-free node target sets. Translational and rotational spring stiffness use their own
  unit quantities. Discrete ELSETs resolve properties and `NSET, ELSET=...`; discrete labels are
  excluded from native mesh element sets, with a warning. Their original labels and memberships
  remain in the retained source. Translated instances remap labels and property set references;
  springs and anisotropic masses in rotated instances are explicitly refused. Oriented/damped masses, nonlinear or
  frequency-dependent springs, SPRING2/A, dashpots and rotary inertia are not mapped yet.
  Gravity may select ordinary mesh elements despite unrelated point masses; targets containing
  discrete elements and whole-model gravity with point masses remain refused.
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
planar incompatible-mode and hybrid formulations are assigned to sections; C3D8 selects B-bar, while C3D8I selects the thirteen-variable formulation in the limited static scope below; C3D8IH remains refused. A section mixing full and
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
material behavior before solving. No hyperelastic or thermal law is inferred. Constant isotropic plastic hardening has the explicit limited mapping below.

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
ranges with bias, modal/subspace harmonic methods, nonempty imaginary loads and prescribed
harmonic motion remain reported as unsupported. Static-history base-state mapping and its
stiffness policy are described below. The `imported_steps` reply
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
defaults, coupling orientations, unsupported plastic histories or nodal `TRANSFORM`. They require the original
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
kinematic global DOF ranges containing at least one translation, or all six
reference DOFs for distributing (omitted data, blank first DOF, or ranges covering
1–6). Explicit distributing DOFs 1–3 additionally map all translations with all
rotations released, limited to small-displacement linear static steps and a reference
point at the equal-weight cloud centroid (roundoff tolerance only; no coordinate
adjustment). Partial rotational release and offset references for this subset remain
refused. Kinematic element surfaces resolve all face nodes (including midside nodes);
element-based distributing weights remain refused. Continuum distribution uses equal nodal weights and translation-derived
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
`KINEMATIC, ALPHA`, rotation-only solid kinematic DOFs, partial rotational distributing DOFs, nonuniform weights, `DISTRIBUTING, COUPLING=STRUCTURAL`,
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

A translation-only original Tet4 fixture fixes the three base nodes and reference
rotations and applies a 3 N reference force in x. With E=120000 Pa, nu=0.25 and
unit coordinate edges, the free node displacement is 0.000125 m and shear stress
is 6 Pa: the centroid force distributes equally to the three cloud nodes. The
same model with all six DOFs coupled gives x displacement 1/104000 m and z
displacement -2/104000 m, distinguishing released rotations from a full coupling.
Both native commands and normalized input round trips check this behavior; the
comparison is a discrete analytical reference, not external solver validation.

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

**C3D8 stiffness mapping:** `C3D8` now selects `formulation.b_bar=true`:
2x2x2 volume-weighted mean dilatation, with the deviatoric/shear strains retained.
Native Hex8 defaults to the original standard full-integration formulation;
its Abaqus export is refused by name rather than silently relabelled `C3D8`.
The flag is saved in the model and recovered on normalized deck read. This is a
physical formulation choice, independent of `lumped_mass`.

The current B-bar scope is Hex8, linear isotropic elasticity, no material
orientation, small-strain static, frequency, harmonic, correlated modal random
or matrix generation. Reduced,
hybrid, incompatible, non-Hex8, nonlinear/material-plastic and other analysis
combinations are refused. Thermal strain and initial-stress recovery are not
implemented. Positive Jacobians are required at every integration point;
tests cover rectangular and moderately warped elements, not arbitrary mesh
quality or a version-wide commercial equivalence claim.

`C3D8I` maps to the distinct thirteen-variable enhanced-strain option for
nonhybrid, small-strain, constant-isotropic static analysis. The original native
nine-mode option is still not exported as C3D8I. C3D8IH remains a named refusal.
Source replay preserves original text separately from semantic export.

Formula: Hughes (1980), DOI 10.1002/nme.1620150914, mean-dilatation specialization;
Nagtegaal, Parks and Rice (1974), CMAME 4, 153–178. The public 2025 solid theory
(https://docs.software.vt.edu/abaqusv2025/English/SIMACAETHERefMap/simathe-c-solidisoquadhex.htm)
provides the C3D8 mapping rationale. This does not certify the unavailable 2024
manual version or finite-strain extensions. No commercial response values are
used to fit the implementation.


### Constant orthotropic stiffness and additional solid-frame rotation

`*ELASTIC, TYPE=ORTHOTROPIC` accepts the nine constant stiffness entries in the
[public keyword order](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-elastic.htm).
The symmetric normal 3×3 block is inverted to reciprocal engineering constants;
shear moduli are retained independently. Positive definiteness, finite numbers,
units and the native material validation are enforced. Temperature/dependency
columns and a plane-stress `LAMINA` completion are not inferred.

For rectangular solid frames, the optional [orientation rotation](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-orientation.htm)
rotates the orthonormal material axes about local axis 1, 2 or 3. The resulting
frame is persisted and exported explicitly; source replay retains the original
axis/angle lines. Rectangular shell frames use the same rotation followed by
projection, as described below. Cylindrical twist and nodal/distributed frames
remain named refusals. This does not relax the isotropic, unoriented B-bar material restriction.

An original rotated tetrahedron tests `Kxx = E2/6`, `ux = 0.003 m` for a 1 N load,
against independently constructed native commands, then save/open and normalized
write/read. Partial kinematic and element-face coupling tests likewise compare
native constraints, displacement, nodal recovered stress and reactions through
these phases. Exact source replay is a separate byte equality assertion; it is
not evidence that every retained keyword is solvable. Tie adjustment/default
tolerances and local coupling axes are still not mapped. Plasticity, thermal
material cards and other previously unimplemented physics remain explicitly
outside this constant-elastic increment.

### Internal-mode discrimination remains a validation gate

Native Hex8 incompatible modes use nine condensed displacement variables and the
centre-Jacobian correction cited in `fem::element` ([Taylor, Beresford and Wilson,
1976](https://doi.org/10.1002/nme.1620100602)). The public
[C3D8I library entry](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEELMRefMap/simaelm-r-3delem.htm)
specifies thirteen internal variables. Condensation removes internal DOFs but
does not prove equality of the enhanced strain spaces. The public
[distortion discussion](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEGSARefMap/simagsa-c-ctmincompat.htm)
also distinguishes parallel and trapezoidal distortions. Matching one rectangular
nu=0 bending case cannot establish equivalence for either distortion or nonzero
Poisson ratio. The nine-mode option remains refused for C3D8I export. Native
nine-mode candidate results, if supplied for independent black-box comparison,
are labelled as such and never as successful C3D8I imports.


### Modal random mass and fixed sampling

Native `solve.random` accepts `mass=consistent|lumped`, default `consistent`.
With `lumped`, an explicit continuum section policy selects row-sum, HRZ or
consistent inertia; the same assembled operator enters mode extraction, harmonic
inertia, Rayleigh material damping and modal damping. Frequency results without
a mass identity are not reused as the modal source. The reply reports the global
selection and section policies. Stored harmonic-frame PSD projection continues
to inherit its source mass and rejects an override. No result-kind allocation or
binary layout changes are involved.

Optional `sampling_frequencies` gives an increasing, bounded SI-Hz replay grid
inside the input modal bands. It changes sampling only; modes, damping and loads
remain defined by the imported model. Its reply is labelled `explicit_replay_grid`;
without it, `source_modal_bands` constructs the usual biased grid from the newly
computed natural frequencies. This distinction matters when comparing material or
inertia choices: fixing only band endpoints does not fix resonance-biased points.

Linear isotropic B-bar now permits harmonic and correlated modal random steps.
Its stiffness and assumed-strain stress recovery are shared with static analysis;
thermal/initial stress, nonlinear materials, orientation and unsupported analysis
combinations remain refused. A constant-simplex one-DOF test independently checks
consistent versus row-sum modal frequencies and every sampled force PSD, including
an incompatible cached source, save/reopen and exact source replay. Normalized
random-step export remains a named refusal; source replay is not semantic export. A mixed
B-bar Hex8/consistent Tet4 matrix test checks inertia integrals and `C=alpha M+beta K`.


### Thirteen-variable Hex8 mapping

Native sections can select `formulation: {enhanced_strain13: true}`. This is
separate from the original nine-variable `incompatible_modes` option and B-bar.
It implements nine enhanced-gradient and four dilatational variables with the
centre-Jacobian determinant correction, full eight-point integration and local
static condensation. The strain used for integration-point stress recovery
includes the recovered internal amplitudes. The formulation is stored in the
native project; it does not change the results binary layout.

The current scope is small-strain linear static analysis, positive-Jacobian Hex8
and constant isotropic elasticity. Reduced integration, hybrid pressure, B-bar,
nine-mode combinations, material axes, thermal/initial stress and other analysis
procedures are refused. Passing patch and energy tests is not proof of accuracy
for arbitrarily distorted bending meshes. `C3D8I` import and normalized export
select this formulation in the stated scope. `C3D8IH` remains refused; neither
frequency nor nonlinear use is inferred from static verification. Exact
source replay, where a source document is retained, is not semantic conversion.
See [the formula and verification scope](enhanced-strain13.md).

### First nonlinear rod equilibrium and frequency about its base state

`*STEP, NLGEOM=YES` (or the flag `NLGEOM`) with automatic `*STATIC` maps to native
nonlinear static for constant isotropic elastic T3D2/T2D2 rods or C3D4 solids.
The current nonlinear solver refuses rods: their steps/base-state links survive
round-trip, but neither nonlinear equilibrium nor dependent frequency is claimed
solved. C3D4 supplies the executable nonlinear static round-trip check. The input increment
and minimum/maximum are divided by the step period; `INC` bounds increments.
Zero/omitted entries follow the public [STATIC defaults](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-static.htm).
The source period remains in the source bundle; native static time is normalized.
A subsequent FREQUENCY step references the preceding equilibrium through `base_state`;
ordinary equilibrium-to-frequency export keeps the sequence and NLGEOM flag.
The executable elastic continuum uses total-Lagrangian kinematics. This mapping
preserves the declared geometric nonlinearity, not a claim of arbitrary finite-strain
commercial material equivalence or identical convergence controls.

Only the first nonlinear equilibrium is imported: sequential state continuation,
other nonlinear solid/shell imports, RIKS, DIRECT, STABILIZE, named nonlinear amplitudes,
instantaneous STEP loading and unsupported controls are named refusals. Their source
text remains available for exact replay; they are not replaced by linear steps.
This ordinary rod path does not implement pretension-section import/export.

### Multi-step endpoint loads and boundary identity

For supported linear static steps, the importer evaluates relative amplitudes at
step end. STEP TIME uses that step's duration and then holds the attained value
if the condition is inherited; TOTAL TIME continues along accumulated general-step
time. A new OP=MOD reference replaces the previous magnitude/amplitude, while
OP=NEW removes the whole corresponding load or boundary family before restating it.
Removing a displacement restraint is different from setting its value to zero.
Normalized export represents these effective equilibrium states; it does not
reconstruct the original amplitude curves or time sampling. Source replay retains
those original cards and all INCLUDE bytes exactly.

Boundary history also keeps its original node/set identity and named/direct format.
Changing a set-based boundary to a node-based one, or changing a named boundary to
direct format without OP=NEW, is refused rather than silently overwriting the same
expanded DOF. Conflicting repeated definitions within one step are refused too.
The [public boundary propagation reference](https://abaqus.uclouvain.be/English/SIMACAEPRCRefMap/simaprc-c-boundary.htm)
explains the identity and removal requirements. Unsupported steps remain in the
source bundle and the import report names the first refusal; preceding supported
steps are not a certificate that the remaining history can be solved.

The command regression uses two instances of an original axial rod and two levels
of INCLUDE, five general steps of durations 2, 3, 1, 2, 1, separate step/total-time
curves, OP=MOD and OP=NEW, and a support transferred from a root to a tip. It checks
each step's expanded constraints and nodal force vectors, analytic displacement
and reaction, saved results, exact source replay and normalized readback. INCLUDE
paths in this workflow are relative to the root input directory, including nested
files; a missing dependency is a named error, never silently skipped. This does
not extend nonlinear paths, transient loading, FIXED or follower-load support.

### Independent direct time history

`DYNAMIC, DIRECT` now maps to the native transient procedure for explicit point
`MASS` and linear grounded `SPRING1` models. It requires explicit `ALPHA=0`,
`BETA=.25`, `GAMMA=.5`, `APPLICATION=TRANSIENT FIDELITY` (or omitted application),
a positive increment and period containing a whole number of increments, and
sufficient `STEP, INC`. The omitted HHT alpha default is **not** average
acceleration and is refused rather than guessed. Initial displacement/velocity
and physical damping are zero. Only homogeneous fixed boundaries and translational
`CLOAD`s sharing one explicit relative TABULAR/SMOOTH STEP amplitude are mapped.
The reference force remains unscaled; the native solver evaluates the complete
curve at each time, including endpoint holding. This independent first step starts
at total time zero, so STEP TIME and TOTAL TIME coincide. Output is retained every
increment; original output-request text remains in source replay.

Run `solve.transient` with `mass=consistent` (the shared default) or `lumped`:
explicit point masses are unchanged by this choice. This does **not** establish
continuous-element mass equivalence. Structural elements, continuation/base state,
nonzero or amplitude-driven boundaries, automatic increments, DIRECT=NO STOP,
other HHT parameters, multiple amplitude names and distributed loads are refused.
Modal time history and response-spectrum imports have the separate restricted
contract below; other combinations remain source-preserved and non-runnable.
No step is silently converted into static equilibrium or another procedure.

Verification uses an original 2 kg / 200 N/m oscillator driven by a triangular
1 N pulse on [0,.2] s and followed to .4 s. Native-command construction, saved
results reopened and normalized INP reread must reproduce all sampled U/V/A.
The independent continuous convolution supplies the physical reference; finite
integration error is distinguished from import error. Byte-exact source replay
checks text preservation separately. Commercial same-input comparison is pending.
Public syntax: [DYNAMIC](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-dynamic.htm)
and [AMPLITUDE](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-amplitude.htm).

### Modal time history and response spectrum

The importer accepts `MODAL DYNAMIC, CONTINUE=NO` (also the default) after one
immediately preceding, unwindowed `FREQUENCY` step with no equilibrium base state.
All extracted modes are used, with unchanged fixed boundary conditions. The mode
count and those conditions compile into native modal re-extraction, **not cached
frequency-result reuse**. The original frequency step and its name remain in the
model/source. Editing that step later does not change the compiled dynamic mode
count. This is currently restricted to explicit MASS/SPRING1 models, whose masses
are identical under `mass=consistent|lumped`; distributed-element mass mappings
are refused. Modal selection, changed boundaries, prestress, continuation and
initial conditions are refused. No nonlinear or alternate procedure is inferred.

`MODAL DAMPING` maps non-overlapping mode-number ranges of fractions of critical
viscous damping; omitted modes have zero damping. Frequency-range, structural,
composite and Rayleigh damping definitions are refused. Independent modal histories
retain the complete common relative TABULAR/SMOOTH STEP amplitude, translational
CLOAD magnitudes and directions, zero initial state and the entire fixed time
interval. STEP TIME and TOTAL TIME coincide in this restricted first perturbation
history. **Native integration is Newmark average acceleration, not Abaqus's modal
integration algorithm.** Continuous-response error must be checked separately;
matching the imported and native-command models does not prove equal time discretization.

`RESPONSE SPECTRUM` accepts explicit `SUM=SRSS|CQC` and `COMP=SRSS|R40` (default COMP
is SRSS). Omitted SUM defaults to ABS in the source format and is refused, not
silently changed. One shared positive constant `SPECTRUM, TYPE=ACCELERATION`
table supplies pseudo-acceleration in the declared length/time-squared units,
frequency in cycles/time, and one critical damping ratio matching every selected
mode. Constant ordinates avoid replacing logarithmic interpolation with native
linear interpolation. The solver refuses frequencies outside the supplied band.
One to three global orthonormal direction vectors and nonnegative scale factors
map to independent directional excitation. Direction combination follows modal
combination. Spatial loads are refused; output is a peak envelope, not a time
history or simultaneous signed state, and acceleration is pseudo-acceleration.
Nonconstant tables, TYPE=G/VELOCITY/DISPLACEMENT, damping tables, distinct spectra
by direction, ABS/DSC and rigid-response additions are explicitly refused.

Response-spectrum T3D2 rods additionally require explicit `mass=lumped` in both
solve and export; import does not select or change the solver mass policy.
Consistent rod mass is not an equivalent T3D2 dynamic export and remains refused.
The modal time-history importer retains its discrete MASS/SPRING1 restriction.

CQC uses the native signed modal quadratic combination with the supplied common
positive critical damping ratio; it is never substituted by SRSS. Zero damping
is refused for CQC. Unequal modal damping remains refused by this single-damping
spectrum-table importer, even though the native CQC kernel supports it. Close-mode
verification compares an original two-mass rod model with the public equal-damping
CQC correlation formula, an explicit SRSS control, save/reopen and normalized replay.

Use `deck.write steps=["<dynamic step>"]` to write the independent analysis;
the native writer supplies the equivalent modal extraction step. All transient
increments are stored even when original output requests selected fewer frames;
the exact original requests remain in source replay. MCRSLT07 kind 6 retains
nodal U/V/A and spectrum peaks as f64. Reopening must preserve these fields exactly;
this does not assert the precision of unrelated result kinds or commercial files.

Validation uses an original two-direction MASS/SPRING1 model, unequal modal
damping and signed loads for the transient, and oblique orthogonal directions
with unequal scales for SRSS/R40 spectra. Checks compare imported/native-command
models, independent damped step-response and peak-combination formulas, native
save/reopen, source replay, and normalized write/read. Commercial comparison is
pending and uses the supplied time-step/peak tolerances.
Sources: [MODAL DYNAMIC](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-modaldynamic.htm),
[MODAL DAMPING](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-modaldamping.htm),
[SPECTRUM](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-spectrum.htm),
[RESPONSE SPECTRUM](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-responsespectrum.htm).

### Pretension and cross-format formulation boundaries

The Abaqus writer checks which procedures can reference equilibrium base states
separately from pretension control history. Frequency and harmonic perturbations
may reference an earlier locked equilibrium; missing lock controls, force-only
bases and unbased frequency steps after pretension are refused. Chain validation
runs before element-type mapping, so an unrelated element refusal cannot conceal
an invalid control history. Standard full Hex8 still has no equivalent C3D8
export: it is refused rather than silently replaced with B-bar. Selecting B-bar
for a linear export fixture does not establish nonlinear pretension equivalence.
Contact closure, force control, locked adjustment persistence and frozen contact
or tie mappings remain native solver checks on their original formulation.

The current Nastran and LS-DYNA writers do not encode the volume-averaged B-bar
section contract. Both reject the named section through the formulation-policy
check. Their ordinary solid cards must not be used as silent substitutes. This
states the implemented writer boundary, not that those products lack alternative
solid formulations.

### Gravity on selected mesh elements in a model with point masses

`DLOAD, GRAV` on an explicit ordinary mesh element label or element set does not
load unrelated point masses. Set aliases are resolved before selecting the loaded elements.
Zero acceleration components remain in the history state for replacement and inheritance,
but do not create invalid zero-vector native loads. Relative amplitudes and `OP=MOD/NEW`
use the same history rules as other supported static loads.

A target containing any imported discrete element (including a mixed mesh/discrete set or
an alias), and a blank whole-model target when point masses exist, are explicitly refused.
The source membership is checked before discrete labels are removed from mesh sets; filtering
cannot silently turn a mixed target into a partial load. The original text is still retained.
This does not implement gravity on selected point masses or expand the supported analysis types.

The regression model is an axial rod with an independent mass and grounded spring: the
rod's uniform body load gives `u_tip = rho*g*L^2/(2*E)` and `RF_root = -rho*A*L*g`,
while the unloaded spring's displacement is zero. These checks distinguish scoped gravity
from whole-model acceleration and cover amplitude scaling, inheritance and replacement.
Syntax reference: [Abaqus public DLOAD reference](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-dload.htm).

### Fixed surface traction

`DSLOAD` with `TRVEC` maps to native surface traction for Tet4 triangular element faces in
small-displacement linear static steps. The vector is given in global coordinates, with
explicit `FOLLOWER=NO` and `CONSTANT RESULTANT=YES`. The reference magnitude multiplies
the normalized direction and the relative amplitude. `OP=NEW` clears both pressure and
traction in the DSLOAD family; an inherited STEP TIME amplitude holds its previous endpoint,
while TOTAL TIME continues to use accumulated step duration. Replacing a direction also
replaces its zero components, so a previous component cannot remain active accidentally.

Other elements, local `ORIENTATION`, follower or current-area controls, `TRSHR`, nonuniform
traction and non-static procedures are named refusals. Defaults that do not explicitly select
the fixed/reference-area form are also refused. This mapping does not add these capabilities
to the nonlinear or transient importer.

The native model keeps the face selection and integrates the constant traction with the
linear triangle shape functions. The existing canonical writer emits the equivalent consistent
nodal forces (area times traction divided by three at each face node), which reproduce the
same linear static FE equations exactly in this scope. It does not emit a general follower
traction card. Source replay preserves the original cards, amplitudes and time data unchanged.
Static results and canonical output represent step endpoints, not the intervening load history.

An original tetrahedron test uses four steps to distinguish STEP TIME endpoint inheritance
from TOTAL TIME evaluation and load replacement; U, RF and constant element stress are
checked against the same native model and `K u = f` with isotropic elasticity. Units are
converted through the shared quantity table. Syntax follows the
[public DSLOAD reference](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-dsload.htm).

### Constant global volume force

`DLOAD` labels `BX`, `BY` and `BZ` map to native `BodyForce` on explicitly selected Tet4
mesh elements in small-displacement linear static steps. Their magnitude is force per volume,
converted with `Quantity::ForcePerVolume`; density is not a multiplier. This is distinct from
`GRAV` acceleration and surface traction. Element-set aliases are resolved; blank targets and
sets containing discrete elements are refused instead of silently changing membership.

Each global component has its own history identity. `OP=MOD` changes that component while
retaining the others; `OP=NEW` removes all prior DLOAD body forces, gravity and element pressure.
Relative STEP TIME and TOTAL TIME amplitudes follow the existing static endpoint rules.
The canonical writer retains the BX/BY/BZ representation in the selected unit system.

Nonuniform labels such as `BXNU`, local `ORIENTATION`, other element types and non-static
procedures are named refusals. Original cards remain in the retained source. An original
four-step tetrahedron regression checks component replacement, inheritance, load-family clearing,
SI and millimetre units, and the exact consistent volume load `f_i = b*V/4` against static U/RF/S.
The load definition follows the [public three-dimensional element reference](https://abaqus.uclouvain.be/English/SIMACAEELMRefMap/simaelm-r-3delem.htm).


Element pressure `DLOAD, Pn` requires an explicit element or element set. Like scoped gravity
and body force, its target is checked against original discrete membership before mesh-set
filtering. A mixed solid/MASS or solid/SPRING1 target is refused, including set aliases; it
cannot be accepted as pressure on only the surviving solid elements. A blank target is refused
instead of being silently interpreted as the whole model. The affected step is not runnable,
and the original source remains available through project reopening and source replay.


### Static history followed by direct harmonic response

A frequency or direct harmonic step after an ordinary `STATIC` step retains that
latest equilibrium as `base_state`, with `linear_base: true`. This explicitly uses
the reference-geometry elastic stiffness without initial-stress stiffness, matching
an input whose general step has `NLGEOM=NO` (including the default). The native
`base_state` default remains stress stiffened; `linear_base` is opt-in and requires
an ordinary linear static base without pretension. It is stored with the model,
included in result invalidation, and exposed by `model.step.create/edit`.
Normalized export retains `NLGEOM=NO`; it refuses a nonlatest base or a base shared
by linear and stress-stiffened consumers because one source general step cannot
express both choices. Nonlinear equilibrium remains a distinct procedure: no YES
flag is changed to NO and unsupported nonlinear rod solving still refuses.

Frequency and harmonic perturbations retain the base support DOFs with zero
incremental values. Inherited static amplitudes are cleared before validating the
new procedure; explicitly supplied unsupported amplitudes still refuse. Each
harmonic perturbation starts with no loads: only loads written in that perturbation act.
The general history remains unchanged for subsequent general steps. Nonzero
prescribed harmonic motion and frequency-dependent amplitudes remain named
refusals. An empty `IMAGINARY, OP=NEW` load block maps an explicitly zero imaginary
component without clearing the real component; nonempty imaginary loads remain
unsupported. No time-domain amplitude is reinterpreted as a frequency curve.

The original axial oscillator uses K=4 N/m and a consistently assembled free mass
of 1 kg. Its static displacement is 0.25 m, frequency 1/pi Hz and response to a
1 N load at 0.1 Hz is 1/(4-(2*pi*0.1)^2) m. Import, solving, saved-result reopening,
source replay and normalized writing/reading retain the complete three-step order.
A second input adds a force only to the static step to distinguish load histories.
Mass is explicitly selected as `consistent` for these comparisons; this does not
infer or apply a source solver's default lumping policy. A separate prestressed
column check distinguishes elastic stiffness from initial-stress stiffening.

Public procedure semantics: [general and perturbation procedures](https://abaqus.uclouvain.be/English/SIMACAEANLRefMap/simaanl-c-linearnonlinear.htm),
[direct steady-state dynamics](https://abaqus.uclouvain.be/English/SIMACAEANLRefMap/simaanl-c-steadystdyndirect.htm),
and [loads in perturbation steps](https://abaqus.uclouvain.be/English/SIMACAEPRCRefMap/simaprc-c-loading.htm).


### Global directional point mass

`MASS, TYPE=ANISOTROPIC` with exactly one row of three constant principal values
maps to the native global diagonal `Mass::Directional([mx,my,mz])`. Values are masses,
not weights: each is converted through the mass unit. Without `ORIENTATION` the
principal directions are the global axes, per the public [MASS keyword](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-mass.htm).
Normalized export writes the three values and `TYPE=ANISOTROPIC`; equal components
may normalize to an isotropic mass with the same physical tensor. Exact source replay
retains the original keyword and numbers separately.

The mapped range is linear static, frequency and direct harmonic analysis. Local
orientation, a rotated instance, mass-proportional/composite damping parameters,
extra temperature/field columns and other analysis combinations are refused by name.
A mass redefinition inside a step is retained as unsupported history, never applied
retroactively to earlier steps. The mass is model data and remains fixed across
supported steps. No constitutive law, continuum mass quadrature or lumping default
is inferred from this point-mass mapping.

An original three-axis oscillator uses masses [2,3,5] kg and grounded spring stiffnesses
[200,600,1500] N/m. Independent formulas are f_i=sqrt(k_i/m_i)/(2*pi) and
U_i=F_i/(k_i-m_i*(2*pi*f)^2). Two direct response steps use forces [2,3,5] N and
[4,0,0] N at 1 Hz. Tests compare native commands, SI/mm-tonne-second inputs,
saved-result reopening, exact source replay and normalized export/read. There is
no continuum stress in this discrete model; harmonic reaction forces are not
exposed by the current probe. This is a mapping/analytical check, not commercial
validation or evidence for general rotated mass tensors.

### Constant isotropic plastic hardening in one small-strain equilibrium

`*PLASTIC` (default or explicit `HARDENING=ISOTROPIC`) accepts constant pairs of
true yield stress and equivalent plastic strain. Stress uses the chosen deck pressure unit;
strain is dimensionless. The native J2 table interpolates between entries and holds the last
yield stress beyond the final entry. The first plastic strain must be zero, strains strictly
increase, stresses are positive and nondecreasing. No guessed elastic strain is subtracted.
[Public PLASTIC keyword](https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-plastic.htm)

Analysis import is limited to exactly one automatic `STATIC`, `NLGEOM=NO`, isotropic C3D4
step with no material orientation. It selects the existing material-nonlinear static solver
with small-displacement kinematics, not the linear elastic solver. Increment sizes are
normalized by the step period; this rate-independent solve uses a proportional ramp from
zero initial plastic strain. Material state advances through converged increments within
that step. It does not transfer plastic history to another step. Preserved results are
endpoint results, not a restart file for integration-point history.

Temperature/field/rate dependence, softening, other hardening types, local orientations,
finite-strain plastic imports, perturbations, multiple steps and named load amplitudes are
refused. Unsupported source history remains replayable; malformed material cards reject
`deck.read` atomically and the original input file remains the evidence. No constitutive
kernel or result layout is added. A positive import is not validation of other plastic laws.

The original one-tetrahedron check has E=1000 Pa, nu=0.25, yield=1 Pa and hardening slope
100 Pa, with 1.5 Pa uniaxial stress. The independent small-strain solution is plastic strain
0.005, axial displacement 0.0065 m over 1 m, and lateral displacement -0.002875 m.
The SI and mm-tonne-second inputs use the same forces, axes and strain table. Native commands,
source replay, saved-result reopening and canonical reimport are separate checks.

### Rectangular shell frames with an additional rotation

A constant rectangular `ORIENTATION` now applies its scalar angle (degrees) about the
selected local axis before the existing shell projection. The cyclic successor axis is
projected into the element plane; the positive element normal completes the right-handed
frame. This order is significant on tilted surfaces. The rotated frame and reference-axis
index are native model data and survive save/open. Canonical output writes that frame with
zero additional angle; source replay instead retains the original angle and text.
[Public orientation algorithm](https://abaqus.uclouvain.be/English/SIMACAEMODRefMap/simamod-c-orientation.htm)

This reuses existing shell material axes; it adds no constitutive law or nonlinear shell
analysis. Cylindrical twist, distributed angles, node-based frames and unsupported coordinate
systems remain named refusals. A degenerate projected direction is rejected by the native
frame validation. The executable check is a flat orthotropic S4 membrane in linear statics:
E1=2000 Pa, E2=1000 Pa, nu12=0.2, thickness=0.1 m, unit edge, total axial force=1 N.
A 90-degree material rotation gives axial displacement 0.01 m and transverse contraction
-0.001 m; the unrotated native control gives axial displacement 0.005 m. SI and mm-tonne-second
inputs are compared with independently constructed native commands, source replay, reopened
endpoint results and a canonical write/read/solve. This is not general shell bending or
finite-rotation verification. Tensor comparisons must state whether components are global,
element-local or material-local; a local stress component is not a global one by name alone.

Plastic state still starts from zero in each independent native nonlinear equilibrium.
A single step's load path is not state transfer between general steps. Multiple-step plastic
import remains refused; the limited single-step hardening mapping above is unchanged.
