# Composite deck mappings

The model stores SI, ply angles in radians, layers from the negative to the positive shell
normal, and the reference surface at `offset * total_thickness` above the laminate centre.
Directions and layer names refer to the material model; the adapters do not tune them against
an external solver result. All tests use original meshes and inputs, without commercial CAE.

| Adapter | Elastic laminate and directions | First-ply strengths |
| --- | --- | --- |
| Abaqus write/read | Conventional shells, rectangular/cylindrical frames, orthotropic engineering constants | Explicitly refused; no damage law is inferred |
| Nastran write | PCOMP/MAT8, rectangular solid frames, explicit rectangular/cylindrical shell frames; one mid-ply sample | Refused |
| LS-DYNA write | MAT_021/PART_COMPOSITE, rectangular solid frames, explicit shell frames; three Simpson points per ply | Refused |

## Abaqus

`SHELL SECTION, COMPOSITE` keeps the layer order and offset sign. Each layer writes its
thickness, section-point count, material and angle in degrees. `ELASTIC, TYPE=ENGINEERING
CONSTANTS` retains E1/E2/E3, ν12/ν13/ν23 and G12/G13/G23. The adapter uses Simpson integration:
an odd count of at least three points per ply keeps the model's exact linear-elastic thickness
moments. The extra sample count is retained, but Simpson's equally spaced output locations
are different from FieldCraft's Gauss samples; only the explicitly reported bottom/middle/top
locations are common. One-point/even-count plies and Gauss section
integration are refused. [Abaqus 2026 SHELL SECTION](https://docs.software.vt.edu/abaqusv2026/English/SIMACAEKEYRefMap/simakey-r-shellsection.htm)

`ORIENTATION, DEFINITION=COORDINATES` carries rectangular frames and cylindrical radial/hoop/
axial frames. A rectangular origin is eliminated by taking coordinate differences. For solids, cylindrical endpoints retain the axis line. Shell cylindrical frames are
evaluated at the element node-average centre and written as a constant rectangular frame,
matching FieldCraft's current discrete shell definition; refinement reevaluates those frames.
The arbitrary scale of a direction vector is not retained.
Nonzero additional rotations, distributions, node-based definitions and unsupported coordinate
systems are refused on import. Solids retain the whole frame. [Abaqus 2026 ORIENTATION](https://docs.software.vt.edu/abaqusv2026/English/SIMACAEKEYRefMap/simakey-r-orientation.htm)

For shells, Abaqus projects the cyclic successor of the specified rotation axis. The adapter
therefore maps `shell_reference_axis` 1/2/3 to rotation axis 3/1/2, with zero additional rotation.
This preserves the projected material 1 direction and the positive element normal.
[Abaqus 2026 Orientations](https://docs.software.vt.edu/abaqusv2026/English/SIMACAEMODRefMap/simamod-c-orientation.htm)

An absent model orientation uses FieldCraft's element frame. Tri3 uses the first edge; Quad4
uses the vector between the midpoints of edges 4–1 and 2–3, with the normal from the crossed
diagonals. These are generally different from the global-axis default in Abaqus. The writer
splits an un-oriented composite or orthotropic shell section into explicit per-element frames
and section sets, while preserving every original load/support set. This is supported for
Tri3/Quad4; higher-order shells need an explicit frame. A default orthotropic solid receives
the explicit global rectangular frame. Import does not guess an omitted orientation for an
orthotropic/composite shell. Generated section names/grouping may differ; physical assignments
and directions are what round-trip.

Strength metadata is not silently discarded. Abaqus `FAIL STRESS` has no separate transverse
shear strength S23. `DAMAGE INITIATION, CRITERION=HASHIN` has six strengths and a selectable
fiber-shear coefficient, but no Tsai–Wu interaction parameter; it is a damage-initiation
material option. This adapter does not select or combine these laws from the model's aggregate
postprocessing strengths. Export of `ply_strength` is refused by material name; import of
`FAIL STRESS`/`FAIL STRAIN` is refused by line. Existing unread damage options remain in the
reader's explicit unread-keyword report.
[Abaqus FailStress data](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKERRefMap/simaker-c-failstresspyc.htm),
[Abaqus 2026 DAMAGE INITIATION](https://docs.software.vt.edu/abaqusv2026/English/SIMACAEKEYRefMap/simakey-r-damageinitiation.htm)

Additional unsupported inputs include symmetric half-stack declarations, named/distributed
ply angles, ply-name metadata, nonlinear ply laws and unrepresentable shell formulations.
The reader rejects unknown physical parameters rather than treating them as defaults. Its
existing unread-history policy is unchanged; it never executes or preserves arbitrary keywords.
These are format mapping checks, not evidence that another solver has run or that its shell
formulation, shear correction, drilling stiffness or output conventions match FieldCraft.

The regression tests are in `crates/decks/tests/abaqus_composite.rs`: four unit presets, material
shear ordering, layer/offset round-trip, rectangular and cylindrical axis selection, skew quads,
cyclic/reversed node order and refusals with source lines. Run:

```text
cargo test -p fieldcraft-decks
```

## Nastran

The Nastran writer emits PCOMP plies in bottom-to-top order with `Z0 = -(0.5 + offset) h`, MAT8 engineering constants and explicit CQUAD/CTRIA material coordinate IDs. CORD2R preserves rectangular frames; a shell cylindrical frame is frozen at its element centre, matching the facet convention used by the kernel. Default orthotropic shell frames are also explicit. Solid rectangular directions use the PSOLID material coordinate system. The public source is the [MSC Nastran 2022.1 Quick Reference Guide](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf), entries PCOMP, MAT8, CQUAD4, CORD2R and PSOLID.

This writer supports the single mid-ply sampling request; it refuses additional through-ply sample counts, independent S23/Hashin strengths, nonlinear plies and spatial cylindrical solid directions instead of dropping them. PCOMP has no field for FieldCraft's extra sample count. It integrates the linear laminate stiffness; the internal one-sample request does not remove each ply's bending stiffness. These are writer/readback checks, not a commercial-solver execution certificate.

## LS-DYNA

The R14 writer uses plain MAT_021 (`MAT_ORTHOTROPIC_THERMAL`) with zero thermal expansion for oriented linear materials, and PART_COMPOSITE_LONG with IRPL=103 for laminates. MAT_002 does not accept a separate ply angle at each integration point; MAT_021 does. Each oriented section receives distinct material definitions, preserving shared input materials without changing the model or its load sets. AOPT=2 carries the explicit reference vector, with bottom-to-top ply angles in degrees, NLOC=2×offset and SHRF=5/6. Solid rectangular frames use both material vectors; cylindrical shell frames are frozen at each facet centre.

This export currently requires three points per ply: LS-DYNA uses bottom/middle/top Simpson points, whereas FieldCraft additionally reports Gauss samples. Their positions must not be compared as identical. Simpson integrates each linear ply's zeroth, first and second stiffness moments exactly; default one-point integration would omit the ply's own t³/12 term. Unequal ply densities, other point counts, cylindrical solid frames, nonlinear oriented materials and the four simultaneous failure measures are named refusals. No averaging of I1/I2 is silently substituted. D3PLOT requests all ply stress points and material-axis output; its strain output remains the format's outermost-surface convention, not FieldCraft's per-ply strain collection.

References: [R14 Volume I](https://ansyshelp.ansys.com/public/Views/Secured/Doc_Assets/Release/v242/LS-DYNA_Manual_Volume_I_R14.pdf), PART_COMPOSITE (35-19–35-30), SECTION_SHELL (39-59), DATABASE_EXTENT_BINARY (16-61–16-63); [R14 Volume II](https://ansyshelp.ansys.com/public/Views/Secured/Doc_Assets/Release/v242/LS-DYNA_Manual_Volume_II_R14.pdf), MAT_021 (2-227 onwards). Tests independently re-read generated cards across four unit systems and verify thickness moments, material/shear order, frame selection and explicit refusals. No commercial solver was run.
