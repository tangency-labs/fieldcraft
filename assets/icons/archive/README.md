# Retired toolbar artwork

These original recipes and VectorCraft/SVG pairs are retained unchanged for history.
`model.orientation` and `model.section.composite` were planned catalogue entries;
material directions and layups are now exposed through `model.section.create`,
with laminate properties through `model.section.properties`.

Only SVGs directly in `assets/icons/` are embedded in the active toolbar, and
`cargo xtask icons` reads recipes directly in `assets/icons/src/`. Retired artwork
lives here so it is neither embedded as a nonexistent command nor regenerated into
the active directory. `cargo xtask assets` still checks its attribution.
