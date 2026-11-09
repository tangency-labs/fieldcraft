# LE7 auxiliary axisymmetric solid comparison

This preserves the original quad8 axisymmetric **solid** recipe and frozen result
byte for byte. It checks only the matching ESRD converged solid reference,
25.27 MPa with the original 1% tolerance. The frozen result is
25.282545696283336 MPa (+0.04965%). It is not verification of a shell formula.

The original recipe omitted `model.units`: its numeric radius 1 and thickness
0.025 are read in the default mm-tonne-s-N-MPa system, so the geometry is 1e-3
of the published metre geometry. E and pressure are 210 GPa and 1 MPa.
Uniform scaling preserves this pressure-only linear stress comparison, but
not displacements, reactions, mass or frequencies. The probe is at the legacy
numeric outer coordinate (1.0125, 1.4034, 0) in millimetres.

The unmodified [original reference contract](history/expected.json) and
[original description and refinement table](history/README.md) are retained
as historical evidence. The old solid-versus-shell comparison was
25.282545696283336 versus 25.86 MPa, -2.2330%, outside 1%. It remains a
historical applicability failure, not a current oracle known failure:
the current [LE7 shell case](../nafems-le7/) uses the validated ring-shell
model. The archived `known_failure` text describes the old state, not current
capabilities. Only this directory's top-level expected.json is executed.

Run `cargo xtask oracle --case nafems-le7-solid`. Do not bless or reinterpret
the solid fixture as a shell result.
