# Carbon/epoxy cantilever panel: measured frequencies

The original 16×16 MITC4 recipe uses the M2 specimen in Tables 1 and 5 of
[Sharma, Mahapatra & Panda (2017)](https://doi.org/10.1590/1679-78253830).
Its 150×150×3.75 mm dimensions and four equal [45/-45/-45/45] plies use the
reported engineering constants and density, with material 1 along x. The x=0
edge is clamped in all six DOFs, the other edges are free, offset is zero,
and the normal is +z. The model uses consistent mass and fixed 5/6 shear correction.

Only the **measured** Table 5 frequencies [62.5,152,378,485,561] Hz are references.
Table 2 describes another material and is not used. No vendor input or numerical
result is used as our model or fixture. The 15% validation band was fixed before
the first solve; it allows physical support/model differences, not just mesh error.
The paper reports imperfect repeatability of boundary conditions but no uncertainty
intervals, so this band is not a confidence interval or general accuracy guarantee.
`solver/tests/composite.rs` independently checks 4/8/16 meshes, with each 8→16
frequency change required below 3%. No effective span, clamp compliance or property
has been inferred by matching the measurements. Frequencies pair by ascending
order: the publication does not supply measured full-field modes for MAC.

Run `cargo xtask oracle --case composite-carbon-epoxy-experiment`.
Bless only after passing: `cargo xtask oracle --bless composite-carbon-epoxy-experiment`.
This does not validate acoustic radiation, layer stress or failure.
