Designs saved by SolidCraft, opened by `tests/void_shells.rs` and `tests/design_formats.rs` as a
person's file would be. Each `<name>.solidcraft` here is `solidcraft-cli run <name>.json --out
<name>.solidcraft` (design format 21); the volumes the tests expect are analytic.

`solidcraft-<version>/` holds designs saved by that SolidCraft release, as the desktop app saves
them (a container): `cd solidcraft-<version> && solidcraft-cli run <name>.json`, where
`solidcraft-<version>/<name>.json` is the recipe `<name>.json` of this folder with a `file.save`
step. From 0.17 on a saved design also carries SolidCraft's model cache (`cache/` in the
container), which FieldCraft never reads: opening rebuilds the timeline.

`solidcraft-<version>/resaved/` holds the format 21 designs of this folder opened and saved
again by that release, as a person saves an older design again:
`cd solidcraft-<version>/resaved && solidcraft-cli run <name>.json`, where `<name>.json` opens
`../../<name>.solidcraft` and saves it. A design saved again keeps each face name on its face. A
recipe built afresh may name the pieces of a split face in another order: from 0.30 SolidCraft
numbers them in its kernel's order for features made in it, while a migrated feature keeps the
old order (`kernel_piece_order` in the saved document).
`tests/design_formats.rs` opens every design here and checks that each saved copy of a recipe gives the bodies and face names of
`<name>.solidcraft` (each name on the same face for a design saved again; the same names for a
recipe built afresh), so a newer SolidCraft's designs keep opening as the older ones did
(docs/architecture.md, "SolidCraft pin"). The scheduled `SolidCraft compatibility` workflow saves
`solidcraft-latest/` and `solidcraft-latest/resaved/` with SolidCraft's latest release the same
way before running these tests.
