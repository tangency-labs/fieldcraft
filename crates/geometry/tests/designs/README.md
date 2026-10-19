Designs saved by SolidCraft, opened by `tests/void_shells.rs` and `tests/design_formats.rs` as a
person's file would be. Each `<name>.solidcraft` here is `solidcraft-cli run <name>.json --out
<name>.solidcraft` (design format 21); the volumes the tests expect are analytic.

`solidcraft-<version>/` holds designs saved by that SolidCraft release, as the desktop app saves
them (a container): `cd solidcraft-<version> && solidcraft-cli run <name>.json`, where
`solidcraft-<version>/<name>.json` is the recipe `<name>.json` of this folder with a `file.save`
step. `tests/design_formats.rs` opens every
design here and checks that each saved copy of a recipe gives the bodies and face names of
`<name>.solidcraft`, so a newer SolidCraft's designs keep opening as the older ones did
(docs/architecture.md, "SolidCraft pin"). The scheduled `SolidCraft compatibility` workflow saves
`solidcraft-latest/` with SolidCraft's latest release the same way before running these tests.
