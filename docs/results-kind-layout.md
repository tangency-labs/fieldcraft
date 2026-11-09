# Result encoding and kind allocation

The current header is **MCRSLT06**. The enclosing `.fieldcraft` ZIP and its part
checksums are unchanged. An eight-byte header precedes a little-endian `u32`
record count. Each name-sorted record contains `u64` UTF-8 name length, name,
32 definition-hash bytes, `u8 kind`, `u64 payload length`, and the bounded serde
codec payload. Both the container checksums and the decoded result are checked.

Kinds 0–5 keep their existing layouts (static, frequency, heat, explicit, fluid,
and pre-mass harmonic). **Kind 6 is linear dynamics**, including direct/modal
histories and response spectra. Kinds 7–10 are allocated to the fatigue/harmonic
extension: static+fatigue, pre-mass harmonic+fatigue, mass-aware harmonic, and
mass-aware harmonic+fatigue, respectively. This unmerged dynamics branch refuses
7–10 by name; it does not implement or trial-decode those payloads. Any new kind
requires Router allocation before implementation.

## Version 06, kind 6: frozen field order

The payload is `transient::Response`, in this order:

| Field | Encoding / meaning |
|---|---|
| mass | `MassKind`, LE u32 enum index: 0 consistent, 1 lumped |
| step | UTF-8 string |
| nodes | sequence of u32 labels, strictly ascending |
| frames | sequence of Frame |
| spectrum | bool; false for histories, true for component peak estimates |
| natural_frequencies | sequence of f64 Hz |
| warnings | sequence of strings |
| energy | optional sequence of Energy, aligned with frames; absent for spectra or recovered old results without a ledger |

Every struct field has a presence byte before its value. Thus a normally written
payload begins with `01`, followed by the four-byte mass index. Strings and
sequences have LE u64 lengths; arrays have no length; f64 values are IEEE 754 LE.
An Option uses a one-byte discriminant. The energy field's `skip_serializing_if`
can instead make its struct presence byte zero. There is no optional trailer or
payload sniffing. Missing or extra fields are not an extension mechanism.

Frame order: optional time (seconds; absent for spectra), displacement (node×3),
rotation (node×3), velocity (node×3, empty for spectra), acceleration (node×3),
optional stress (node×6), optional kinetic energy, optional strain energy. Energy
order: damping_dissipation, external_work, energy_balance, all f64 joules. Values
are SI; stress order is xx, yy, zz, xy, xz, yz. The mass choice is stored with the
result, not inferred from the model or a default on reopening. A dynamics result
is not sufficient by itself to distinguish direct from modal integration: that
procedure remains in the model step, bound by the definition hash.

## Historical boundaries

| Header/kind | Historical lane-16 payload | Current reader |
|---|---|---|
| 05 / 0–5 | Existing non-dynamics payloads | Read unchanged |
| 05 / 6 | step, nodes, frames, spectrum, natural_frequencies, warnings; implicit consistent mass, no cumulative energy ledger | Refused |
| 05 / 7 | Above plus optional energy; implicit consistent mass | Refused |
| 05 / 10 | mass first, then step, nodes, frames, spectrum, natural_frequencies, warnings, optional energy | Refused |
| 06 / 6 | Full mass-aware Response described above | Read/write |
| 06 / 7–10 | Allocated foreign extension layouts | Reserved; unavailable in this branch |

The old 05 numbers overlapped across independent branches. Consequently **all
05 kinds above 5 are refused before payload decoding**, including a byte stream
that would happen to decode successfully. Version 06 kind 6 is intentionally
not the original 05 kind 6 payload. Its mass-aware payload matches lane 16's
former 05 kind 10 field order, but the new header and allocated kind remove the
ambiguous interpretation. Never patch a magic/kind byte to claim migration.

Keep original archives, probes, source/binary identities and hashes. Recover
identified old files with a verified old reader, exporting the model without
results and typed result data, then validate the model definition hash and
re-encode with the current native writer. The ignored `result_transfer` test
helper supports this explicitly approved evidence workflow; it is not an
application fallback or a producer detector. Unknown producers remain refused.
No solve is necessary for an identified typed transfer. Readers older than 06
may open the model with an unsupported-result warning and omit the results.

Updating this branch is **not a cross-branch merge**. The independent fatigue
branch can reserve/refuse kind 6 until the user-directed integration brings in
its exact writer/reader. Any later field-order or variant-order change requires
an explicit version contract and coordination, not another unallocated kind.

## Interim fusion-p4a merge (not the unified 07 delivery)

The merged static payload includes lane 15's pretension field. Version 05 kind 0
uses its explicit pre-pretension decoder; current 06 static fixtures come from
that merge. Lane 16's earlier 06 static archives do not acquire that field by
guessing. Preserve them with their producer identity. Linear dynamics kind 6
is unchanged in this interim merge. The coordinator has superseded the earlier
allocation with a unified 07 contract to be delivered by lane 18; no archive
written here is represented as 07, and no new kind is allocated locally.
