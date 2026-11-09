# Unified result kinds and producer identity

`results.bin` now writes **MCRSLT06** (8 ASCII bytes), followed by little-endian
`u32 step_count`. Each name-sorted record is `u64 name_bytes`, UTF-8 name,
32 definition-hash bytes, `u8 kind`, `u64 payload_bytes`, then the existing bounded
serde binary codec payload. The outer `.fieldcraft` ZIP manifest still records the
part size and SHA256. This change does not alter model/deck physics or solver data.

| Kind | Version 06 payload |
|---|---|
| 0 | StaticResult without fatigue |
| 1 | FrequencyResult |
| 2 | HeatResult |
| 3 | ExplicitResult |
| 4 | FluidResult |
| 5 | Pre-mass HarmonicResult field order, without fatigue |
| 6 | Reserved for lane 16 linear dynamics; this unmerged fatigue build refuses it by name |
| 7 | Tuple (StaticResult, FatigueField) |
| 8 | Tuple (pre-mass HarmonicResult, FatigueField) |
| 9 | HarmonicResult with mass and operator_diagonal_bounds |
| 10 | Tuple (HarmonicResult with mass and operator_diagonal_bounds, FatigueField) |

The pre-mass harmonic order is step, nodes, reference_points, frames,
natural_frequencies, residual_vectors, stats, warnings. It implies consistent
mass, with no operator-diagonal diagnostic. The writer uses kinds 5/8 only when
mass is consistent and the diagnostic is absent; otherwise it uses 9/10, never
losing lumped identity. The ordinary harmonic solver produces the diagnostic, so
new solves use 9/10. Fatigue is omitted from the base result's codec and encoded
once as the second tuple member for 7/8/10. No optional trailer is introduced:
separate registered kinds plus an explicit format version are the smaller repair.

## Old files: do not infer a producer from a number

MCRSLT05 kinds **0–5** remain readable with their unchanged payload layouts.
MCRSLT05 kinds **6 and above are refused before payload decoding**, even when a
payload could successfully decode as one of the known types. Historically lane
18 used 6=static+fatigue, 7=pre-mass harmonic+fatigue, 8=mass harmonic,
9=mass harmonic+fatigue, while lane 16 used overlapping dynamics numbers. Numeric
shape, decode success and file extension do not identify the producer.

Keep the original file and SHA256. Recover needed results with its independently
verified old source/binary identity, export via that program's supported commands,
and regenerate with the new program. A verified producer can export the model
without results and its existing probes for comparison; only small necessary
examples need solving again. Never edit a kind or magic byte as a migration.
Original BASE results need no unsolicited recomputation. Old versions of the
application do not understand version 06 results; keep originals when exchanging
files across versions.

## Cross-lane integration

Kind 6 belongs exclusively to linear dynamics. The user-authorized intermediate
merge must preserve the lane 16 decoder/writer at 6 under the same version 06
header; this branch does not import its solver. All later kind assignments require
Router allocation before implementation. Do not add a fallback decoder for old
version 05 kind 6 in a merged build: the original producer identity is still needed.

The inp-reader work separates random-frame evaluation and element-local stress
recovery and fixes modal damping symmetry. An integration must retain these paths
alongside fatigue's stored-result mass identity; its mass-requirement report does
not select a solver option. Lane 17's 5R calculation waits for the user-directed
intermediate merge and an explicit `solve.harmonic mass` choice. Do not treat
`lumped` as automatic agreement with every imported mass-distribution formula.
