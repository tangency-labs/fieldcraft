# Spring stiffness sensitivity

One node, a unit point mass and a grounded spring k=4 N/m, guided in y,z.
A unit x load gives u=1/4 m and du/dk=d(Fu)/dk=-1/16. The eigenvalue is
k/m and df/dk=f/(2k)=1/(8π). The model and numerical references are original.
No spatial discretization error exists. The 1e-8 relative tolerance covers
central assembly differencing at relative increment 1e-5, not FE error.

This verifies the command route and quotient units in SI. Solid stress,
material constants and nonzero prescribed motion have solver integration tests;
The bounded scalar frequency fit starts at k=4 N/m and recovers k=9 N/m
from its synthetic exact frequency, with the required 1e-6 relative tolerance.
None of this certifies clustered-mode derivatives, joint-parameter or harmonic
inversion, or a large-model benchmark.
