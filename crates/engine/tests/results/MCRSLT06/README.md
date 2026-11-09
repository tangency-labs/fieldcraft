# Version 06 result fixtures

Generated on Linux by the repository's existing original small-model generator:

```
cargo test -p fieldcraft-engine --lib cmd::file::tests::write_the_results_fixtures -- --ignored --exact
```

These are newly solved/serialized examples, not edited version 05 files:
cantilever (static/frequency), couette (fluid), impact (explicit), and blocks
(nonlinear static/contact). Keep the version 05 originals unchanged. The engine
reading test checks both directories. IO's `result_kinds_tests` additionally
covers kinds 0–5 and 7–10, reserved 6 and rejection of ambiguous version 05 kinds.
See `docs/results-kind-layout.md` for the exact envelope and registered kinds.
