# Joint E and k recovery on a beam

The original composite reference-surface mesh and constraints are retained.
A unit normal point force acts at the centre; two transverse displacements
are frozen from a separate truth solve. Truth E=1.2 GPa, k=13000.0 N/m^3;
initial E=1 GPa, k=10000.0 N/m^3. Tangential stiffness is zero.
All observations and parameter bounds/scales are explicit in the recipe.
Replay identifies parameters without repeating the truth solve.

The acceptance is relative parameter error <=1e-6 with data-only rank 2.
This tests the actual command and composite-shell joint sensitivity path.
It is a synthetic consistency test, not independent physical validation.

The beam initially failed its iteration cap with sigma=1e-11 m. That failed
run remains in the lane evidence. The same observations were then whitened
with sigma=1e-6 m, giving the successful recorded setup used here. Sigma is
a weighting scale for this noiseless test, not a claimed measurement noise.
