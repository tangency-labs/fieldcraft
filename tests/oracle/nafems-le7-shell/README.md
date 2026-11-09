# NAFEMS LE7: axisymmetric shell with a regular pole

Independent command model, explicit m-kg-s-N-Pa units. R=1 m, h=.025 m,
E=210 GPa, nu=.3, density7800kg/m³ (irrelevant to this static result), cylinder
length1.5m capped by a hemisphere. Two-ring conical shell elements on beam2
connectivity use 100 intervals per metre (rounded upward per region).
D=(1,1.4034,0) is a midsurface node; the target outer face is at R=1.0125.
Connectivity runs upward from F to A, so the section normal points inward and
bottom is the outer surface. The measured value is element D's first endpoint
sigma_YY, without averaging it with the adjacent element. Both endpoint values
and a 50/100/200 refinement check are documented in the axisymmetric shell manual.

F fixes axial translation and symmetry rotation. Pole radial translation and
rotation are automatically constrained by the production solver. Pressure1MPa
acts on the thin-shell midsurface as independently integrated consistent ring
nodal forces: 2*pi*p*integral N*r*(dz,-dr). No pressure-face follower behavior is
needed in this linear static problem. Its axial resultant is pi*p*R².

The 25.86MPa target and 1% tolerance come from the same public shell reference as
the old LE7 case (ESRD2018, p.23). This distinct shell input does not alter the
old solid recipe, its default-millimetre geometry or its shell-reference failure.
See docs/solvers/axisymmetric-shell.md for formulas, pole limits and restrictions.
