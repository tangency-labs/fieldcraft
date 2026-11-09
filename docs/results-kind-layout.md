# Unified version 07 results

The only current writer uses **MCRSLT07**. The little-endian envelope is eight
ASCII magic bytes, `u32 step_count` (at most 10,000), then name-sorted records:
`u64 name_bytes`, UTF-8 name, 32 definition-hash bytes, `u8 kind`,
`u64 payload_bytes`, bounded serde binary codec payload. No layout sentinel,
optional trailer, alternate development-kind aliases or trial decoders are used.
The ZIP manifest verifies part sizes and SHA256.

| Kind | Version 07 payload |
|---|---|
| 0 | StaticResult without fatigue, including pretension and nested nonlinear/contact results |
| 1 | FrequencyResult |
| 2 | HeatResult |
| 3 | ExplicitResult |
| 4 | FluidResult |
| 5 | HarmonicResult including mass and operator_diagonal_bounds |
| 6 | Full Dynamics Response: frames, spectrum flag, natural frequencies, warnings, energy ledger and mass |
| 7 | Tuple (StaticResult, FatigueField) |
| 8 | Tuple (HarmonicResult including mass and operator_diagonal_bounds, FatigueField) |

Fatigue is skipped in the base result's codec and encoded exactly once as the
second tuple member. Kinds 5 and 8 always carry mass, including consistent results
without an operator diagnostic. Kind 6 stores the complete direct/modal history
or response-spectrum result; it is not a placeholder refusal. There are no current
kinds 9 or 10. Additional kind assignments require Router approval.

For decoders, type 5 field order is `mass, operator_diagonal_bounds, step, nodes,
reference_points, frames, natural_frequencies, residual_vectors, stats, warnings`.
Type 6 order is `mass, step, nodes, frames, spectrum, natural_frequencies, warnings,
energy`. The fatigue member of 7/8 is `damage, life_seconds, inputs_json,
dependencies`. The codec writes structure-field presence bytes, little-endian
fixed-width numbers, `u64` list/string lengths, and `u32` enum variants; skipped
fatigue fields in the base types emit no member. A complete decoder must consume
the entire length-delimited payload, including optional values and nested fields.

## Old version 05 and 06

Only kinds **0–5** are accepted using their explicit old schemas:

- 05/static uses the frozen pre-pretension static structure and restores an empty
  pretension vector. 06/static uses the pre-merge prestress structure including
  pretension (the frozen `fusion-p4a` producer).
- 1–4 retain their original payloads in both versions.
- 05/06 harmonic kind 5 has field order step, nodes, reference_points, frames,
  natural_frequencies, residual_vectors, stats, warnings. It restores consistent
  mass and no operator diagnostic; it is not read as version 07 kind 5.
- Every old kind **6 or above** is refused before decoding its payload, with the
  literal message **统一前的开发文件，请用原程序重新生成**.

The pre-merge fatigue producer also wrote static kind 0 under 06 without the
prestress field. That incompatible file is not guessed or retried: retain it and
recover using its verified original program. Accepting version 06 does not mean
all incompatible development producers had the same static schema. Historical
06/7 or 06/10 dynamics and fatigue collisions never receive aliases in this reader.
The withdrawn intermediate unified-06 sentinel experiment is not an accepted
schema. Other unsupported magic versions retain the existing result-drop note.

Keep each original file, source/program identity and SHA256. Where a result is
still needed, use its verified original program to export a result-free model and
reference replies, or use an independently retained original command recipe;
regenerate through the version 07 native solver/writer and compare all fields.
Never change old magic/kind bytes as migration. Current built-in test fixtures
are generated under MCRSLT07; old05/06 fixtures stay in their historical folders.
Do not recalculate large BASE responses merely to repair metadata.

## Mass and cross-lane integration

`solve.transient` (direct/modal), `solve.spectrum` and `solve.harmonic` share
`MassKind`: `consistent|lumped`, default `consistent` on each solve. The selection
feeds mass assembly, damping and modal projection and is stored in types 5/6/8.
The random-frame pure evaluator and unaveraged element recovery coexist with
session-bound fatigue persistence. Random/fatigue inherit their solved source's
mass; a separate random `mass` override is refused. Imported modal random steps
retain their supported default-consistent contract and refuse external overrides.
No PSD is turned into a signed tensor history or nonlinear principal-stress spectrum.

`solve.run mass` retains its scoped contract: only frequency, transient and
response-spectrum steps together; other mixed cases use individual commands.
The integration does not start lane 17's separately scheduled 5R or any BASE run.

## Native evidence tool

`cargo run --release -p fieldcraft-engine --example results_roundtrip -- input.fieldcraft new-directory`
opens through the native reader, rewrites through the native writer, compares the
entire decoded result objects, and exports complete fields plus archive hashes.
It does not solve, infer a legacy producer, repair invalid data or replace a
physical validation. The built-in fixture generator and IO synthetic-layout tests
separately cover old compatibility, all current kinds and explicit refusals.
